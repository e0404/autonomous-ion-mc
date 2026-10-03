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
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.config import (
    BACKENDS,
    DELTA_ELECTRON_MODELS,
    MCS_MODELS,
    PRECISIONS,
    STRAGGLING_MODELS,
    EffectiveConfig,
    SimulationConfig,
    validate,
)
from ionmc.environment import describe_environment
from ionmc.errors import TransportLimitError
from ionmc.scoring import MEV_PER_G_TO_GY, ScoringGrid, reduce_batches, voxel_mass_g
from ionmc.transport.reference import RawTransport, run_reference


@dataclass(frozen=True)
class TransportCounters:
    """Transport-limit counters; any nonzero value invalidates the result.

    ``step_truncation``: histories stopped by ``max_steps``; ``stall``: histories stopped after
    more than three consecutive zero-length steps; ``straggling_rejection``: straggling samples
    that exceeded 64 attempts; ``genealogy_overflow`` and ``queue_overflow``: secondary
    bookkeeping limits (always 0 until secondaries exist); ``source_energy_out_of_range``:
    sampled source energies outside ``[E_cut, table maximum]``; ``energy_inversion``: steps in
    which the inverse-range energy exceeded the initial energy (clamped to it).
    """

    step_truncation: int = 0
    stall: int = 0
    straggling_rejection: int = 0
    genealogy_overflow: int = 0
    queue_overflow: int = 0
    source_energy_out_of_range: int = 0
    energy_inversion: int = 0

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
        }


@dataclass(frozen=True)
class EnergyBalance:
    """Energy bookkeeping summed over all histories [MeV] (not per primary).

    Every tally is accumulated independently. ``cutoff_mev`` (energy deposited locally below
    ``E_cut``) is part of the grid and outside deposits; the closure identity is
    ``initial = step_deposit + cutoff + escaped + truncated + unaccounted`` and, for every
    grid, ``in_grid + outside = step_deposit + cutoff``.
    """

    initial_mev: float
    step_deposit_mev: float
    cutoff_mev: float
    escaped_mev: float
    truncated_mev: float
    unaccounted_mev: float
    in_grid_mev: tuple[float, ...]
    outside_mev: tuple[float, ...]

    @property
    def closure_residual_mev(self) -> float:
        """``initial`` minus the sum of all independently tallied destinations."""
        return self.initial_mev - (
            self.step_deposit_mev
            + self.cutoff_mev
            + self.escaped_mev
            + self.truncated_mev
            + self.unaccounted_mev
        )

    @property
    def relative_residual(self) -> float:
        """``|closure residual| / initial`` (0 if nothing was simulated)."""
        if self.initial_mev == 0.0:
            return 0.0
        return abs(self.closure_residual_mev) / self.initial_mev

    def grid_relative_residual(self, grid: int) -> float:
        """``|in_grid + outside - (step_deposit + cutoff)| / initial`` for grid ``grid``."""
        if self.initial_mev == 0.0:
            return 0.0
        total = self.step_deposit_mev + self.cutoff_mev
        return abs(self.in_grid_mev[grid] + self.outside_mev[grid] - total) / self.initial_mev


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


@dataclass(frozen=True, eq=False)
class Result:
    """Result of a transport run (see :class:`GridResult`, :class:`EnergyBalance`).

    ``valid`` is False if any transport-limit counter is nonzero. ``timings`` in seconds:
    ``setup``, ``compile`` (0 for the Python backend), ``transport``, ``reduce``, ``total``.
    ``diagnostics`` holds the optional diagnostics requested in the configuration.
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

    def grid(self, name: str) -> GridResult:
        """The grid result named ``name``."""
        for g in self.grids:
            if g.name == name:
                return g
        raise KeyError(f"no scoring grid named {name!r}")


def capabilities() -> dict[str, Any]:
    """What the engine supports (reported, not negotiated: unsupported requests raise)."""
    return {
        "species": ["proton"],
        "backends": {
            "python": "available (float64 reference)",
            "warp-cpu": "not available (V3-003B)",
            "warp-cuda": "not available (V3-003B)",
        },
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
        "scoring": ["ScoringGrid (energy, dose)"],
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
        raw = run_reference(eff)
        t_transport = time.perf_counter()
        result = _assemble(eff, raw)
        t_end = time.perf_counter()
        timings = {
            "setup": self._setup_s,
            "compile": 0.0,
            "transport": t_transport - t_start,
            "reduce": t_end - t_transport,
            "total": self._setup_s + (t_end - t_start),
        }
        result = _with_timings(result, timings)
        if not result.valid and not self.config.run.allow_invalid_result:
            nonzero = {k: v for k, v in result.counters.as_dict().items() if v}
            raise TransportLimitError(
                f"transport limits violated {nonzero}; the result is invalid "
                "(set RunOptions.allow_invalid_result to receive it anyway)",
                result,
            )
        return result


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
    )


def _grid_result(
    grid: ScoringGrid, batch_sums: NDArray[np.float64], eff: EffectiveConfig, valid: bool
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
    )


def _assemble(eff: EffectiveConfig, raw: RawTransport) -> Result:
    cfg = eff.requested
    counters = TransportCounters(**raw.counters)
    valid = not counters.any_nonzero
    grids = tuple(_grid_result(g, raw.edep_mev[i], eff, valid) for i, g in enumerate(cfg.scoring))
    t = raw.tallies
    balance = EnergyBalance(
        initial_mev=t["initial"],
        step_deposit_mev=t["step_deposit"],
        cutoff_mev=t["cutoff"],
        escaped_mev=t["escaped"],
        truncated_mev=t["truncated"],
        unaccounted_mev=t["unaccounted"],
        in_grid_mev=tuple(float(a.sum()) for a in raw.edep_mev),
        outside_mev=tuple(float(x) for x in raw.outside_mev),
    )
    return Result(
        valid=valid,
        requested_config=cfg,
        effective_config=eff,
        backend=eff.backend,
        device="python (float64)",
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
    )
