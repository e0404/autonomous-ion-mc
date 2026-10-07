"""Simulation driver, results and capabilities.

``Simulation(config)`` validates at construction (nothing is simulated if validation fails)
and ``run()`` returns a :class:`Result`. All reported quantities are per primary history:
energy in MeV, dose in Gy (``1 MeV/g = 1.602176634e-10 Gy``), computed with the exact
voxel-mass overlap of :mod:`ionmc.scoring`.

Fail-closed behaviour (decision 0037): a nonzero transport-limit counter makes the result
invalid; ``run()`` then raises :class:`~ionmc.errors.TransportLimitError` (carrying the invalid
result) unless ``RunOptions.allow_invalid_result`` is set, in which case the returned result has
``valid = False`` and every grid result is flagged ``valid = False``.
"""

from __future__ import annotations

import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.config import (
    BACKENDS,
    CHANNEL_BACKENDS,
    DEFAULT_CHUNK_HISTORIES,
    DELTA_ELECTRON_MODELS,
    MAX_CPU_WORKERS,
    MCS_MODELS,
    MIN_CHUNK_HISTORIES,
    PRECISIONS,
    STRAGGLING_MODELS,
    EffectiveConfig,
    SimulationConfig,
    cuda_available,
    validate,
)
from ionmc.environment import describe_environment
from ionmc.errors import TransportLimitError
from ionmc.lookup import AXES as LOOKUP_AXES
from ionmc.lookup import AXIS_SPACINGS as LOOKUP_AXIS_SPACINGS
from ionmc.physics.projectiles import PROTON
from ionmc.scoring import (
    GENERATION_CHOICES,
    MEV_PER_G_TO_GY,
    TALLY_QUANTITIES,
    QuantityResult,
    ScoringGrid,
    reduce_batches,
    reduce_ratio,
    voxel_mass_g,
)
from ionmc.species import producible
from ionmc.transport.channels import (
    CLASS_LOCAL,
    CLASS_STEP,
    EXCLUDED_CHANNEL_NAME,
    FE_F_MAX_EXPONENT,
    FE_F_MIN_EXPONENT,
    LOCAL_PIECE_COUNT_NAME,
    MAX_CHANNELS,
    MAX_PARTICLES,
    MAX_SPECTRUM_BINS,
    PIECE_COUNT_NAME,
)
from ionmc.transport.run import run_transport
from ionmc.transport.tally import NUCLEAR_TALLY_NAMES, ChannelRaw, RawTransport


@dataclass(frozen=True)
class TransportCounters:
    """Transport-limit counters; any nonzero value invalidates the result.

    ``step_truncation``: histories stopped by ``max_steps``; ``stall``: histories stopped after
    more than three consecutive zero-length steps; ``straggling_rejection``: straggling samples
    that exceeded 64 attempts; ``genealogy_overflow`` and ``queue_overflow``: secondary
    bookkeeping limits (always 0 until secondaries exist); ``source_energy_out_of_range``:
    sampled source energies outside ``[E_cut, table maximum]``; ``energy_inversion``: steps in
    which the inverse-range energy exceeded the initial energy (clamped to it);
    ``accumulator_overflow``: a fixed-point voxel accumulator reached its capacity;
    ``scoring_pieces_overflow``: a leg of the track-length scoring crossed more voxel pieces than
    the validated bound (the remainder is deposited in the last voxel and the result is invalid).
    """

    step_truncation: int = 0
    stall: int = 0
    straggling_rejection: int = 0
    genealogy_overflow: int = 0
    queue_overflow: int = 0
    source_energy_out_of_range: int = 0
    energy_inversion: int = 0
    accumulator_overflow: int = 0
    scoring_pieces_overflow: int = 0

    @property
    def any_nonzero(self) -> bool:
        """True if any counter is nonzero."""
        return any(self.as_dict().values())

    def as_dict(self) -> dict[str, int]:
        """The counters as a dictionary."""
        return {
            "step_truncation": self.step_truncation,
            "stall": self.stall,
            "straggling_rejection": self.straggling_rejection,
            "genealogy_overflow": self.genealogy_overflow,
            "queue_overflow": self.queue_overflow,
            "source_energy_out_of_range": self.source_energy_out_of_range,
            "energy_inversion": self.energy_inversion,
            "accumulator_overflow": self.accumulator_overflow,
            "scoring_pieces_overflow": self.scoring_pieces_overflow,
        }


@dataclass(frozen=True)
class NuclearTransportCounters(TransportCounters):
    """The counters of a ``nuclear=True`` run (decision 0041 section 5): the transport-limit
    counters plus the conditional block. ``majorant_violation``: candidates with
    ``Sigma(E1) > S^(E0)``; ``nuclear_rejection_limit``: events that exhausted 64 attempts;
    ``nuclear_conservation``: events whose ledger ``T1 = sum T_lab + T_r + binding + imbalance``
    missed by more than 1e-9 T1 (should never fire). Any nonzero value invalidates the result."""

    majorant_violation: int = 0
    nuclear_rejection_limit: int = 0
    nuclear_conservation: int = 0

    def as_dict(self) -> dict[str, int]:
        """The counters as a dictionary (with the nuclear block)."""
        return {
            **super().as_dict(),
            "majorant_violation": self.majorant_violation,
            "nuclear_rejection_limit": self.nuclear_rejection_limit,
            "nuclear_conservation": self.nuclear_conservation,
        }


@dataclass(frozen=True)
class EnergyBalance:
    """Energy bookkeeping summed over all histories [MeV] (not per primary).

    Every tally is accumulated independently. ``cutoff_mev`` (energy deposited locally below
    ``E_cut``) is part of the grid and outside deposits; the closure identity is
    ``initial = step_deposit + cutoff + escaped + truncated + unaccounted`` and, for every
    grid, ``in_grid + quantization + outside = step_deposit + cutoff`` (``quantization`` is the
    summed rounding residual of the fixed-point grids, at most q/2 per deposited piece).
    """

    initial_mev: float
    step_deposit_mev: float
    cutoff_mev: float
    escaped_mev: float
    truncated_mev: float
    unaccounted_mev: float
    in_grid_mev: tuple[float, ...]
    outside_mev: tuple[float, ...]
    quantization_mev: tuple[float, ...] = ()

    @property
    def nuclear_destinations_mev(self) -> float:
        """Energy in the conditional nuclear destinations (0 here; see the nuclear subclass)."""
        return 0.0

    @property
    def local_deposit_mev(self) -> float:
        """Local deposits counted in the grid identity besides step deposit and cutoff (0 here)."""
        return 0.0

    @property
    def closure_residual_mev(self) -> float:
        """``initial`` minus the sum of all independently tallied destinations."""
        return self.initial_mev - (
            self.step_deposit_mev
            + self.cutoff_mev
            + self.escaped_mev
            + self.truncated_mev
            + self.unaccounted_mev
            + self.nuclear_destinations_mev
        )

    @property
    def relative_residual(self) -> float:
        """``|closure residual| / initial`` (0 if nothing was simulated)."""
        if self.initial_mev == 0.0:
            return 0.0
        return abs(self.closure_residual_mev) / self.initial_mev

    def grid_relative_residual(self, grid: int) -> float:
        """``|in_grid + quantization + outside - (step_deposit + cutoff)| / initial``."""
        if self.initial_mev == 0.0:
            return 0.0
        total = self.step_deposit_mev + self.cutoff_mev + self.local_deposit_mev
        quant = self.quantization_mev[grid] if self.quantization_mev else 0.0
        return abs(self.in_grid_mev[grid] + quant + self.outside_mev[grid] - total) / (
            self.initial_mev
        )


@dataclass(frozen=True)
class NuclearEnergyBalance(EnergyBalance):
    """Energy balance of a ``nuclear=True`` run (decision 0041 section 5, amended 2026-10-07):
    ``initial = step_deposit + cutoff + nuclear_local + escaped + nuclear_escaped_neutron +
    nuclear_escaped_gamma + nuclear_binding + nuclear_imbalance + truncated + unaccounted`` and,
    for every grid, ``in_grid + quantization + outside = step_deposit + cutoff + nuclear_local``.
    ``nuclear_mev`` holds the six ``NUCLEAR_TALLY_NAMES`` (``nuclear_alpha_local`` is the alpha
    part of ``nuclear_local``, not a separate destination; ``nuclear_binding`` and
    ``nuclear_imbalance`` are signed)."""

    nuclear_mev: Mapping[str, float] = field(default_factory=dict)

    @property
    def nuclear_destinations_mev(self) -> float:
        n = self.nuclear_mev
        return (
            n["nuclear_local"]
            + n["nuclear_escaped_neutron"]
            + n["nuclear_escaped_gamma"]
            + n["nuclear_binding"]
            + n["nuclear_imbalance"]
        )

    @property
    def local_deposit_mev(self) -> float:
        return float(self.nuclear_mev["nuclear_local"])


@dataclass(frozen=True, eq=False)
class GridResult:
    """Results of one scoring grid; arrays have shape ``(nx, ny, nz)`` (batch array
    ``(B, nx, ny, nz)``). Everything is per primary.

    ``energy_mev``/``energy_std_mev``: mean deposited energy per primary and its standard
    error [MeV]; ``dose_gy``/``dose_std_gy``: dose per primary [Gy] and standard error;
    ``relative_uncertainty``: ``std / mean`` (NaN where the mean is 0: an undefined value is
    never reported as zero); ``defined_mask``: voxels with positive mean energy and mass;
    ``mass_g``: exact voxel mass [g] (0 outside the geometry; dose is NaN there);
    ``n_nonzero_batches``: batches with a nonzero sum; ``sparse_mask``: defined voxels with
    fewer than ``B/2`` nonzero batches (their uncertainty is unreliable);
    ``low_batch_count``: ``B < 10`` (variance estimate itself is poor); ``batch_energy_mev``:
    per-batch estimates ``E_b / (N/B)``; ``valid``: False if the run was invalid (the arrays
    are then not physically meaningful).
    """

    name: str
    energy_mev: NDArray[np.float64]
    energy_std_mev: NDArray[np.float64]
    dose_gy: NDArray[np.float64]
    dose_std_gy: NDArray[np.float64]
    relative_uncertainty: NDArray[np.float64]
    defined_mask: NDArray[np.bool_]
    mass_g: NDArray[np.float64]
    n_nonzero_batches: NDArray[np.int32]
    sparse_mask: NDArray[np.bool_]
    low_batch_count: bool
    batch_energy_mev: NDArray[np.float64]
    valid: bool
    quantities: Mapping[str, QuantityResult] = field(default_factory=lambda: MappingProxyType({}))


@dataclass(frozen=True, eq=False)
class Result:
    """Result of a transport run (see :class:`GridResult`, :class:`EnergyBalance`).

    ``valid`` is False if any transport-limit counter is nonzero. ``timings`` in seconds:
    ``setup``, ``compile`` (0 for the Python backend), ``transport``, ``reduce``, ``total``.
    ``diagnostics`` holds the optional diagnostics requested in the configuration;
    ``transport_report`` holds backend facts (workers, per-worker device, compile seconds,
    chunking, register count) that are not part of the physics result.
    """

    valid: bool
    requested_config: SimulationConfig
    effective_config: EffectiveConfig
    backend: str
    device: str
    precision: str
    seed: int
    rng: dict[str, Any]
    n_histories: int
    n_batches: int
    grids: tuple[GridResult, ...]
    energy_balance: EnergyBalance
    counters: TransportCounters
    timings: dict[str, float]
    environment: dict[str, Any]
    diagnostics: dict[str, Any] = field(default_factory=dict)
    transport_report: dict[str, Any] = field(default_factory=dict)
    channel_raw: ChannelRaw | None = None
    lookups: tuple[dict[str, Any], ...] = ()

    def channel_batches(self, index: int) -> NDArray[np.int64]:
        """Integer per-batch accumulators ``[B, size]`` of channel ``index`` of the compiled plan
        (``effective_config.channels``), in units of its quantum ``2^-k``."""
        plan = self.effective_config.channels
        if plan is None or self.channel_raw is None:
            raise KeyError("this result has no scoring channels")
        c = plan.channels[index]
        out: NDArray[np.int64] = self.channel_raw.acc[:, c.offset : c.offset + c.size].copy()
        return out

    def grid(self, name: str) -> GridResult:
        """The grid result named ``name``."""
        for g in self.grids:
            if g.name == name:
                return g
        raise KeyError(f"no scoring grid named {name!r}")


def _tally_capabilities(nuclear: bool = False) -> dict[str, Any]:
    """The scoring-channel part of the capability report (decision 0040). Every entry is derived
    from the objects that enforce it (``CHANNEL_BACKENDS``, ``TALLY_QUANTITIES``, ``AXES``,
    ``producible``), so the report cannot drift from the validation. With ``nuclear=True`` the
    producible species, the generations and their note follow the nuclear engine (secondary
    protons and deuterons, generation >= 1, are transported and accepted)."""
    pairs = sorted(producible(PROTON, nuclear=nuclear))
    producible_generations = sorted({g for _, g in pairs})
    # nuclear=True runs on the Python reference only (warp kernels are V3-005B): the tally
    # backends are limited accordingly and the warp ones listed as deferred
    tally_backends = ("python",) if nuclear else CHANNEL_BACKENDS
    deferred = {"deferred_backends": {b: "V3-005B" for b in CHANNEL_BACKENDS if b != "python"}}
    return {
        "quantities": list(TALLY_QUANTITIES),
        "backends": {
            b: {
                "quantities": list(TALLY_QUANTITIES),
                "available": b != "warp-cuda" or cuda_available(),
            }
            for b in tally_backends
        },
        **(deferred if nuclear else {}),
        "producible": [{"species": n, "generation": g} for n, g in pairs],
        "generations": {
            "accepted": [
                g for g in GENERATION_CHOICES if g == "all" or g in producible_generations
            ],
            "rejected": [
                g for g in GENERATION_CHOICES if g != "all" and g not in producible_generations
            ],
            "note": (
                "secondary protons and deuterons (generation >= 1) are transported and accepted "
                "(nuclear=True, decision 0041)"
                if nuclear
                else "secondary particles are not transported yet (V3-005A); no secondary species"
            ),
        },
        "let_medium": ["water"],
        "dose_reference": ["medium"],
        "lookup": {
            "axes": list(LOOKUP_AXES),
            "axis_spacings": list(LOOKUP_AXIS_SPACINGS),
            "uniform_axis_required": True,
            "non_uniform": "rejected; resample_uniform() is the explicit, recorded alternative",
            "values": "per species, finite and non-negative; clinical tables and RBE models "
            "do not ship",
            "max_value_range": [2.0**FE_F_MIN_EXPONENT, 2.0**FE_F_MAX_EXPONENT],
            "max_value_range_note": "the largest value of a requested table must lie in "
            f"[2^{FE_F_MIN_EXPONENT}, 2^{FE_F_MAX_EXPONENT}] (fixed-point scale and bound "
            "representable); other magnitudes fail closed",
        },
        "fluence_spectrum": {
            "edges": "uniform in energy or in ln energy (MeV per nucleon), under/overflow bins",
            "max_bins": MAX_SPECTRUM_BINS,
        },
        "automatic_channels": [PIECE_COUNT_NAME, LOCAL_PIECE_COUNT_NAME, EXCLUDED_CHANNEL_NAME],
        "max_channels": MAX_CHANNELS,
        "fail_closed": [
            (
                "unknown or unproducible species, generation 'secondary' (accepted with "
                "nuclear=True: secondary protons and deuterons are transported)"
                if nuclear
                else "unknown or unproducible species, generation 'secondary'"
            ),
            "let_medium other than water, dose_reference other than medium",
            "unknown grid or lookup, duplicate or reserved request names, unused lookups",
            "lookup species gaps, axis coverage gaps, sha256 mismatch, non-uniform tables or "
            "spectrum edges, negative or non-finite lookup values",
            "a quantum above the precision floor, accumulator memory above the budget",
            "a backend without channels"
            + (" (warp backends with nuclear=True until V3-005B)" if nuclear else ""),
        ],
    }


E_CUT_DEUTERON_DEFAULT_MEV = 4.0  # PhysicsOptions.e_cut_deuteron_mev (checked by a test)


def _nuclear_capabilities() -> dict[str, Any]:
    """The ``nuclear`` section of the capability report (only with ``nuclear=True``)."""
    pairs = sorted(producible(PROTON, nuclear=True))
    return {
        "backends": ["python"],
        "warp_backends": "not before V3-005B (rejected before transport)",
        "source": "proton only; E0 + 6 sigma_E <= 250 MeV",
        "stopping": "analytic Bethe (deuteron table needed); 'nist-star' is rejected",
        "table": "derived nuclear-proton table by id (loaded, re-hashed, source pins checked)",
        "elements": "every element of every material needs an evaluated or surrogate entry",
        "producible": [{"species": n, "generation": g} for n, g in pairs],
        "secondaries": ["proton", "deuteron"],
        "local_deposit": "nuclear_local (alpha, residual recoil; generation = parent + 1)",
        "deuteron_cutoff_mev_default": E_CUT_DEUTERON_DEFAULT_MEV,
        "max_particles_per_history": MAX_PARTICLES,
        "capacity": "B_L and the E bound follow decision 0041 section 5 (amended 2026-10-08)",
    }


def capabilities(nuclear: bool = False) -> dict[str, Any]:
    """What the engine supports (reported, not negotiated: unsupported requests raise).
    ``nuclear=True`` adds the ``nuclear`` section; the default report is unchanged."""
    report = _capabilities_base()
    if nuclear:
        report["nuclear"] = _nuclear_capabilities()
        report["physics"]["nuclear"] = True
        report["tallies"] = _tally_capabilities(nuclear=True)
    return report


def _capabilities_base() -> dict[str, Any]:
    return {
        "species": ["proton"],
        "backends": {
            "python": "available (float64 reference; cpu_workers > 1 splits over processes)",
            "warp-cpu": "available (float32 production, float64 validation; cpu_workers > 1)",
            "warp-cuda": (
                "available (cuda:0, float32 and float64; chunked launches)"
                if cuda_available()
                else "not available (no usable CUDA device; no fallback)"
            ),
        },
        "trace": "float64 only (float32 with trace_histories > 0 is rejected)",
        "chunk_histories": {"default": DEFAULT_CHUNK_HISTORIES, "minimum": MIN_CHUNK_HISTORIES},
        "max_cpu_workers": MAX_CPU_WORKERS,
        "backend_names": list(BACKENDS),
        "precisions": list(PRECISIONS),
        "physics": {
            "nuclear": False,
            "stopping": "tabulated or analytic electronic stopping (CSDA, inverse range)",
            "straggling_models": list(STRAGGLING_MODELS),
            "mcs_models": list(MCS_MODELS),
            "delta_electrons": list(DELTA_ELECTRON_MODELS),
        },
        "geometry": ["VoxelGeometry", "BoxPhantom"],
        "sources": ["PencilBeamSource"],
        "scoring": [
            "ScoringGrid (energy, dose)",
            "TallyRequest (edep, dose, fluence, LET, lookup averages, fluence spectra; "
            "species and generation channels; see 'tallies')",
        ],
        "tallies": _tally_capabilities(),
        "max_scoring_grids": 4,
    }


class Simulation:
    """A validated simulation. Construction raises before any transport if the request is
    unsupported (``UnsupportedCombinationError`` / ``BackendUnavailableError``)."""

    def __init__(self, config: SimulationConfig) -> None:
        t0 = time.perf_counter()
        self.config = config
        self.effective: EffectiveConfig = validate(config)
        self._setup_s = time.perf_counter() - t0

    def run(self) -> Result:
        """Run the transport and return the result (raises on invalid results unless allowed)."""
        t_start = time.perf_counter()
        eff = self.effective
        raw = run_transport(eff)
        t_transport = time.perf_counter()
        result = _assemble(eff, raw)
        t_end = time.perf_counter()
        timings = {
            "setup": self._setup_s,
            "compile": _compile_seconds(raw),
            "transport": t_transport - t_start,
            "reduce": t_end - t_transport,
            "total": self._setup_s + (t_end - t_start),
        }
        result = _with_timings(result, timings)
        if not result.valid and not self.config.run.allow_invalid_result:
            nonzero = {k: v for k, v in result.counters.as_dict().items() if v}
            if result.channel_raw is not None and result.channel_raw.lookup_out_of_domain:
                nonzero["lookup_out_of_domain"] = result.channel_raw.lookup_out_of_domain
            if result.channel_raw is not None and result.channel_raw.path_bound_exceeded:
                nonzero["path_bound_exceeded"] = result.channel_raw.path_bound_exceeded
            raise TransportLimitError(
                f"transport limits violated {nonzero}; the result is invalid "
                "(set RunOptions.allow_invalid_result to receive it anyway)",
                result,
            )
        return result


def _compile_seconds(raw: RawTransport) -> float:
    """Kernel compile/load seconds: the slowest worker's load plus, for a pool, the parent's
    compile before the workers were spawned (0 for the Python backend)."""
    parts = raw.meta.get("partials", [])
    worker = max((float(p.get("compile_s", 0.0)) for p in parts), default=0.0)
    parent = max((float(p.get("parent_compile_s", 0.0)) for p in parts), default=0.0)
    return worker + parent


def _device_description(eff: EffectiveConfig, raw: RawTransport) -> str:
    parts = raw.meta.get("partials", [])
    if eff.backend == "python":
        base = "python (float64)"
    else:
        base = str(parts[0].get("device", eff.backend)) if parts else eff.backend
    workers = eff.requested.run.cpu_workers
    return f"{base} x {workers} processes" if workers > 1 else base


def _with_timings(result: Result, timings: dict[str, float]) -> Result:
    return Result(
        valid=result.valid,
        requested_config=result.requested_config,
        effective_config=result.effective_config,
        backend=result.backend,
        device=result.device,
        precision=result.precision,
        seed=result.seed,
        rng=result.rng,
        n_histories=result.n_histories,
        n_batches=result.n_batches,
        grids=result.grids,
        energy_balance=result.energy_balance,
        counters=result.counters,
        timings=timings,
        environment=result.environment,
        diagnostics=result.diagnostics,
        transport_report=result.transport_report,
        channel_raw=result.channel_raw,
        lookups=result.lookups,
    )


def _channel_sums(eff: EffectiveConfig, ch: ChannelRaw, index: int) -> NDArray[np.float64]:
    """Per-batch sums ``[B, size]`` of channel ``index`` in the channel unit (float64)."""
    assert eff.channels is not None
    c = eff.channels.channels[index]
    block = ch.acc[:, c.offset : c.offset + c.size].astype(np.float64)
    return block * c.quantum


def ratio_rounding_bound(
    xhat: NDArray[np.float64],
    a: NDArray[np.float64],
    yhat: NDArray[np.float64],
    b: NDArray[np.float64],
    rhat: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Deterministic fixed-point bound of a ratio ``R = X / Y`` by interval arithmetic: with the
    quantized sums ``X in [xhat - a, xhat + a]``, ``Y in [yhat - b, yhat + b]``, ``xhat, yhat >= 0``
    the bound is ``max(|(xhat + a)/(yhat - b) - rhat|, |(xhat - a)/(yhat + b) - rhat|)``. If
    ``yhat - b <= 0`` the true denominator may be zero and the bound is ``+inf`` (it never
    silently becomes 0). All arguments are per primary (the level where the counts are defined)."""
    lower = yhat - b
    with np.errstate(divide="ignore", invalid="ignore"):
        hi = np.abs((xhat + a) / lower - rhat)
        lo = np.abs((xhat - a) / (yhat + b) - rhat)
    return np.where(lower > 0.0, np.maximum(hi, lo), np.inf)


def _grid_quantities(gi: int, eff: EffectiveConfig, raw: RawTransport) -> dict[str, QuantityResult]:
    """Quantities of grid ``gi`` from the channels (decision 0040): linear ones through
    ``reduce_batches``, ratios through ``reduce_ratio``, plus the automatic channels."""
    plan, ch = eff.channels, raw.channels
    assert plan is not None and ch is not None
    cfg = eff.requested
    grid = cfg.scoring[gi]
    run = cfg.run
    hpb = run.n_histories // run.n_batches
    mass = voxel_mass_g(grid, eff.geometry).reshape(-1)
    vol = float(np.prod(grid.spacing_mm))
    nvox = grid.n_voxels
    n_step_ci = plan.count_channel(gi, CLASS_STEP)
    n_local_ci = plan.count_channel(gi, CLASS_LOCAL)
    # pieces per primary and voxel, step and local class (exact integer channels)
    n_step_pp = _channel_sums(eff, ch, n_step_ci).mean(axis=0) / hpb
    n_local_pp = _channel_sums(eff, ch, n_local_ci).mean(axis=0) / hpb

    def n_pp(ci: int) -> NDArray[np.float64]:
        """Pieces per primary and voxel that channel ``ci`` can receive (by its class mask)."""
        total = np.zeros(nvox)
        for ni in plan.piece_count_indices(ci):
            total = total + (n_step_pp if ni == n_step_ci else n_local_pp)
        return total

    out: dict[str, QuantityResult] = {}

    def shape_of(a: NDArray[Any], extra: int) -> NDArray[Any]:
        return a.reshape((*grid.shape, extra)) if extra else a.reshape(grid.shape)

    def linear(
        name: str, quantity: str, ci: int, units: str, definition: str, scale: NDArray[Any],
        lookup: dict[str, Any] | None = None, positive_scale_only: bool = False,
    ) -> QuantityResult:  # fmt: skip
        c = plan.channels[ci]
        sums = _channel_sums(eff, ch, ci)
        stats = reduce_batches(sums, hpb)
        bins = c.size // nvox
        sc = np.repeat(scale, bins) if bins > 1 else scale
        mean = stats.mean * sc
        std = np.sqrt(stats.variance_of_mean) * sc
        bound = np.repeat(n_pp(ci), bins) * (c.quantum / 2.0) * sc
        defined = mean > 0.0
        if positive_scale_only:
            defined = defined & np.repeat(mass > 0.0, bins)
        extra = bins if bins > 1 else 0
        return QuantityResult(
            name, quantity, "linear", units, definition,
            shape_of(mean, extra), shape_of(std, extra), shape_of(defined, extra),
            shape_of(stats.n_nonzero, extra),
            shape_of(bound, extra), (c.k,), lookup,
        )  # fmt: skip

    ones = np.ones(nvox)
    for qd in plan.quantities:
        req = qd.request
        if req.grid != grid.name:
            continue
        lk = None if req.lookup is None else plan.lookups[
            [t.name for t in plan.lookups].index(req.lookup)
        ].provenance()  # fmt: skip
        if qd.kind == "linear":
            if req.quantity == "dose":
                with np.errstate(divide="ignore"):
                    sc = np.where(
                        mass > 0.0, MEV_PER_G_TO_GY / np.where(mass > 0, mass, 1.0), np.nan
                    )
                out[req.name] = linear(
                    req.name, "dose", qd.numerator, "Gy", qd.definition, sc, lk, True
                )
            elif req.quantity in ("fluence", "fluence_spectrum"):
                out[req.name] = linear(
                    req.name, req.quantity, qd.numerator, "mm^-2", qd.definition + " per volume",
                    ones / vol, lk,
                )  # fmt: skip
            else:
                out[req.name] = linear(
                    req.name, req.quantity, qd.numerator, qd.units, qd.definition, ones, lk
                )  # fmt: skip
        else:
            assert qd.denominator is not None
            xs = _channel_sums(eff, ch, qd.numerator)
            ys = _channel_sums(eff, ch, qd.denominator)
            st = reduce_ratio(xs, ys)
            cn, cd = plan.channels[qd.numerator], plan.channels[qd.denominator]
            xbar = xs.mean(axis=0) / hpb
            ybar = ys.mean(axis=0) / hpb
            bound = ratio_rounding_bound(
                xbar, n_pp(qd.numerator) * (cn.quantum / 2.0),
                ybar, n_pp(qd.denominator) * (cd.quantum / 2.0), st.mean,
            )  # fmt: skip
            bound = np.where(st.defined_mask, bound, np.nan)  # inf stays inf where defined
            out[req.name] = QuantityResult(
                req.name, req.quantity, "ratio", qd.units, qd.definition,
                st.mean.reshape(grid.shape), np.sqrt(st.variance_of_mean).reshape(grid.shape),
                st.defined_mask.reshape(grid.shape), st.n_nonzero.reshape(grid.shape),
                bound.reshape(grid.shape), (cn.k, cd.k), lk,
            )  # fmt: skip
    for ci in plan.channel_index("E", gi):
        if plan.channels[ci].class_mask == CLASS_LOCAL:
            out[EXCLUDED_CHANNEL_NAME] = linear(
                EXCLUDED_CHANNEL_NAME, "edep", ci, "MeV",
                "energy of cutoff and zero-length deposits, excluded from LET and lookups", ones,
            )  # fmt: skip
    out[PIECE_COUNT_NAME] = linear(
        PIECE_COUNT_NAME, "count", n_step_ci, "1", "scoring pieces of class step per primary", ones
    )
    out[LOCAL_PIECE_COUNT_NAME] = linear(
        LOCAL_PIECE_COUNT_NAME, "count", n_local_ci, "1",
        "scoring pieces of class local (point deposits) per primary", ones,
    )  # fmt: skip
    return out


def _grid_result(
    grid: ScoringGrid,
    batch_sums: NDArray[np.float64],
    eff: EffectiveConfig,
    valid: bool,
    quantities: Mapping[str, QuantityResult] | None = None,
) -> GridResult:
    run = eff.requested.run
    stats = reduce_batches(batch_sums, run.n_histories // run.n_batches)
    shape = grid.shape
    mass = voxel_mass_g(grid, eff.geometry)
    flat_mass = mass.reshape(-1)
    mean = stats.mean
    std = np.sqrt(stats.variance_of_mean)
    defined = (mean > 0.0) & (flat_mass > 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(mean > 0.0, std / mean, np.nan)
        scale = np.where(flat_mass > 0.0, MEV_PER_G_TO_GY / flat_mass, np.nan)
    dose = np.where(flat_mass > 0.0, mean * scale, np.nan)
    dose_std = np.where(flat_mass > 0.0, std * scale, np.nan)
    sparse = defined & (stats.n_nonzero < 0.5 * run.n_batches)
    per_batch = (batch_sums / float(run.n_histories // run.n_batches)).reshape(
        (run.n_batches, *shape)
    )
    return GridResult(
        name=grid.name,
        energy_mev=mean.reshape(shape),
        energy_std_mev=std.reshape(shape),
        dose_gy=dose.reshape(shape),
        dose_std_gy=dose_std.reshape(shape),
        relative_uncertainty=rel.reshape(shape),
        defined_mask=defined.reshape(shape),
        mass_g=mass,
        n_nonzero_batches=stats.n_nonzero.reshape(shape),
        sparse_mask=sparse.reshape(shape),
        low_batch_count=run.n_batches < 10,
        batch_energy_mev=per_batch,
        valid=valid,
        quantities=MappingProxyType(dict(quantities or {})),
    )


def _assemble(eff: EffectiveConfig, raw: RawTransport) -> Result:
    cfg = eff.requested
    counters = (
        NuclearTransportCounters(**raw.counters)
        if eff.nuclear is not None
        else TransportCounters(**raw.counters)
    )
    ood = 0 if raw.channels is None else raw.channels.lookup_out_of_domain
    pbe = 0 if raw.channels is None else raw.channels.path_bound_exceeded
    valid = not counters.any_nonzero and ood == 0 and pbe == 0
    used = set() if eff.channels is None else {c.grid for c in eff.channels.channels}
    grids = tuple(
        _grid_result(
            g,
            raw.edep_mev[i],
            eff,
            valid,
            _grid_quantities(i, eff, raw) if i in used else None,
        )  # fmt: skip
        for i, g in enumerate(cfg.scoring)
    )
    t = raw.tallies
    balance_cls = EnergyBalance if eff.nuclear is None else NuclearEnergyBalance
    extra: dict[str, Any] = (
        {}
        if eff.nuclear is None
        else {"nuclear_mev": {k: float(t[k]) for k in NUCLEAR_TALLY_NAMES}}
    )
    balance = balance_cls(
        initial_mev=t["initial"],
        step_deposit_mev=t["step_deposit"],
        cutoff_mev=t["cutoff"],
        escaped_mev=t["escaped"],
        truncated_mev=t["truncated"],
        unaccounted_mev=t["unaccounted"],
        in_grid_mev=tuple(float(a.sum()) for a in raw.edep_mev),
        outside_mev=tuple(float(x) for x in raw.outside_mev),
        quantization_mev=tuple(float(x) for x in raw.quantization_mev),
        **extra,
    )
    return Result(
        valid=valid,
        requested_config=cfg,
        effective_config=eff,
        backend=eff.backend,
        device=_device_description(eff, raw),
        precision=eff.precision,
        seed=cfg.run.seed,
        rng=dict(eff.rng),
        n_histories=cfg.run.n_histories,
        n_batches=cfg.run.n_batches,
        grids=grids,
        energy_balance=balance,
        counters=counters,
        timings={},
        environment=describe_environment(),
        diagnostics=raw.diagnostics,
        transport_report=raw.meta,
        channel_raw=raw.channels,
        lookups=()
        if eff.channels is None
        else tuple(lk.provenance() for lk in eff.channels.lookups),
    )
