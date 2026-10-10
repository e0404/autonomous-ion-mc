"""Simulation configuration, fail-closed validation and the effective configuration.

``validate(config)`` implements the fail-closed capability contract of decision 0039: every
request outside the supported envelope raises ``UnsupportedCombinationError`` (or
``BackendUnavailableError`` for backends that are not available) before any transport, and
nothing is silently clamped, defaulted or substituted. The result of ``validate`` is the
:class:`EffectiveConfig`, which records what will actually run.

Units follow decision 0037 (MeV, mm, g/cm3, rad); names carry the unit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Any, Literal

from ionmc._validate import choice, fail, integer, real
from ionmc.errors import BackendUnavailableError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.lookup import LookupTable
from ionmc.materials import WATER, Material
from ionmc.physics.em import E_S_MEV
from ionmc.physics.projectiles import DEUTERON, PROTON, Projectile
from ionmc.physics.scattering import scattering_length_g_cm2
from ionmc.physics.stopping import StoppingSource, StoppingTable
from ionmc.scoring import MAX_SCORING_GRIDS, ScoringGrid, TallyRequest
from ionmc.sources import PencilBeamSource
from ionmc.species import producible
from ionmc.transport.channels import ChannelPlan, NuclearCapacity, compile_channels
from ionmc.transport.tables import TransportTables, thaw
from ionmc.transport.tally import (
    MAX_QUANTA,
    QUANTUM_MEV,
    TRACE_N_CONTINUOUS,
    TRACE_N_DISCRETE,
)

STRAGGLING_MODELS = ("bohr_gamma_v1", "bohr_gauss_clamped_gamma_v1")
MCS_MODELS = ("differential_moliere",)
DELTA_ELECTRON_MODELS = ("local",)
BACKENDS = ("python", "warp-cpu", "warp-cuda")
PRECISIONS = ("float32", "float64")
MAX_REJECTION_ATTEMPTS = 64
DEFAULT_CHUNK_HISTORIES = 2**18
MIN_CHUNK_HISTORIES = 2**10
MAX_CPU_WORKERS = 256
MAX_TRACE_BUFFER_BYTES = 2**30
MAX_SCORING_PIECES = 4096
MAX_ENERGY_SIGMA_FRACTION = 0.05
MAX_ENERGY_LOSS_FRACTION = 0.2
CHANNEL_BACKENDS: tuple[str, ...] = ("python", "warp-cpu", "warp-cuda")
"""Backends that implement the scoring channels of decision 0040; a configuration with tallies is
rejected on every other backend, so a requested tally is never silently ignored."""
NUCLEAR_MAX_ENERGY_MEV = 250.0

U01_MAPPING = {
    "float64": "((w >> 8) + 0.5) * 2**-24",
    "float32": "((w >> 9) + 0.5) * 2**-23",
}


@dataclass(frozen=True)
class PhysicsOptions:
    """Physics switches and numerical parameters.

    ``nuclear`` has no default (the caller must state it). ``True`` (decision 0041; python
    backend, proton source only) additionally needs ``nuclear_table_id`` (the derived nuclear
    table, loaded and verified by ``validate()``) and uses ``e_cut_deuteron_mev`` (the deuteron
    cutoff, an engineering default of 4 MeV = twice the 2 MeV total-energy Bethe table floor).
    ``elastic`` (default ``True``, V3-005C, decision 0041 slice C) adds the hadronic elastic channel
    (p-p and p + A, all backends since C4) to a ``nuclear=True`` run: ``Sigma_tot =
    Sigma_nonel + Sigma_el`` in the thinning, the p-p slower proton as a transported secondary and
    the p + A recoil deposited locally (``elastic_recoil_local``). ``elastic=False`` reproduces the
    pre-D2 non-elastic-only behaviour bit for bit; without ``nuclear=True`` the switch has no
    effect. ``elastic_table_id`` names the derived elastic table (default: the id pinned in
    ``src/ionmc/data/elastic_table_pin.json``); it is loaded, re-hashed and qualification-checked
    fail-closed. ``elastic_only`` (internal, row V11) switches the non-elastic channel off so that
    only EM and elastic scattering act; it needs ``nuclear=True`` and ``elastic=True``.
    ``stopping`` supplies the electronic stopping tables (for example an
    analytic :class:`~ionmc.physics.stopping.BetheStoppingSource`). Below ``e_cut_mev`` the
    remaining kinetic energy is deposited locally; ``max_step_mm`` bounds the step length;
    ``max_energy_loss_fraction`` bounds the mean energy loss per step as a fraction of the
    kinetic energy; ``range_alpha`` and ``range_rho_f_mm`` parametrise the Geant4 range step
    function; steps shorter than ``short_step_fraction`` of the residual range use the
    linear loss ``S(E_mid) t`` with ``E_mid = E - S(E) t / 2`` (midpoint rule).
    ``truncated_hinge_diagnostic`` (default off, a diagnostic for the T14 negative control, not a
    physics model) is "truncate-first": the planned step ends at the first
    plane the straight line reaches, the angle is sampled for that length, leg 2 is not cut again
    and the end point is snapped onto the plane (the displacements are recorded per history; the
    direction is never changed).
    """

    nuclear: bool
    stopping: StoppingSource
    straggling: bool = True
    straggling_model: str = "bohr_gamma_v1"
    multiple_scattering: bool = True
    mcs_model: str = "differential_moliere"
    delta_electrons: str = "local"
    e_cut_mev: float = 2.0
    max_step_mm: float = 1.0
    max_energy_loss_fraction: float = 0.02
    range_alpha: float = 0.2
    range_rho_f_mm: float = 0.1
    short_step_fraction: float = 1.0e-2
    truncated_hinge_diagnostic: bool = False
    e_cut_deuteron_mev: float = 4.0
    nuclear_table_id: str | None = None
    elastic: bool = True
    elastic_only: bool = False
    elastic_table_id: str | None = None

    def __post_init__(self) -> None:
        for name in (
            "nuclear",
            "straggling",
            "multiple_scattering",
            "truncated_hinge_diagnostic",
            "elastic",
            "elastic_only",
        ):
            if not isinstance(getattr(self, name), bool):
                raise fail(f"{name} must be a bool, got {getattr(self, name)!r}")
        for name in ("straggling_model", "mcs_model", "delta_electrons"):
            if not isinstance(getattr(self, name), str):
                raise fail(f"{name} must be a string")
        if not (hasattr(self.stopping, "table") and hasattr(self.stopping, "name")):
            raise fail("stopping must be a StoppingSource (with .name and .table())")
        real("e_cut_mev", self.e_cut_mev, positive=True)
        real("e_cut_deuteron_mev", self.e_cut_deuteron_mev, positive=True)
        if self.nuclear_table_id is not None and not isinstance(self.nuclear_table_id, str):
            raise fail("nuclear_table_id must be a string or None")
        if self.elastic_table_id is not None and not isinstance(self.elastic_table_id, str):
            raise fail("elastic_table_id must be a string or None")
        if self.elastic_only and not (self.nuclear and self.elastic):
            raise fail("elastic_only needs nuclear=True and elastic=True")
        real("max_step_mm", self.max_step_mm, positive=True)
        f = real("max_energy_loss_fraction", self.max_energy_loss_fraction, positive=True)
        if f > MAX_ENERGY_LOSS_FRACTION:
            raise fail(f"max_energy_loss_fraction must be <= {MAX_ENERGY_LOSS_FRACTION}, got {f}")
        a = real("range_alpha", self.range_alpha, positive=True)
        if a > 1.0:
            raise fail(f"range_alpha must be in (0, 1], got {a}")
        real("range_rho_f_mm", self.range_rho_f_mm, positive=True)
        s = real("short_step_fraction", self.short_step_fraction, positive=True)
        if s > 0.01:
            raise fail(f"short_step_fraction must be in (0, 0.01], got {s}")


@dataclass(frozen=True)
class RunOptions:
    """Execution options. ``seed`` is a 64-bit unsigned integer; ``n_histories`` primaries are
    split into ``n_batches`` interleaved batches (history ``h`` belongs to batch ``h mod B``).
    ``max_steps`` (per history) defaults to a computed bound that is recorded. A nonzero
    transport-limit counter raises ``TransportLimitError`` unless ``allow_invalid_result``.
    ``memory_budget_bytes`` bounds the energy-deposit accumulators (of all worker processes
    together). ``cpu_workers`` > 1 splits the histories over that many spawned worker processes
    on the ``python`` and ``warp-cpu`` backends (``warp-cuda`` rejects it); a worker that fails
    or exceeds ``worker_timeout_s`` terminates the run without a partial result.
    ``chunk_histories`` (a power of two, at least 2**10, recorded in the effective
    configuration) bounds the histories per Warp launch; it affects memory only: the deposit grids
    are int64 fixed-point quanta with integer atomic adds, so the result is bit-identical for
    any chunk size within a precision."""

    backend: Literal["python", "warp-cpu", "warp-cuda"]
    precision: Literal["float32", "float64"]
    seed: int
    n_histories: int
    n_batches: int = 20
    cpu_workers: int = 1
    max_steps: int | None = None
    allow_invalid_result: bool = False
    worker_timeout_s: float | None = None
    memory_budget_bytes: int = 2**31
    chunk_histories: int = DEFAULT_CHUNK_HISTORIES

    def __post_init__(self) -> None:
        choice("backend", self.backend, BACKENDS)
        choice("precision", self.precision, PRECISIONS)
        integer("seed", self.seed)
        integer("n_histories", self.n_histories)
        integer("n_batches", self.n_batches)
        integer("cpu_workers", self.cpu_workers, minimum=1)
        if self.max_steps is not None:
            integer("max_steps", self.max_steps, minimum=1)
        if not isinstance(self.allow_invalid_result, bool):
            raise fail("allow_invalid_result must be a bool")
        if self.worker_timeout_s is not None:
            real("worker_timeout_s", self.worker_timeout_s, positive=True)
        integer("memory_budget_bytes", self.memory_budget_bytes, minimum=1)
        chunk = integer("chunk_histories", self.chunk_histories, minimum=MIN_CHUNK_HISTORIES)
        if chunk & (chunk - 1) != 0:
            raise fail(f"chunk_histories must be a power of two, got {chunk}")


@dataclass(frozen=True)
class DiagnosticsOptions:
    """Optional diagnostics (all off by default): the end position and reason of every track,
    records of escaping particles, and a per-step trace of the first ``trace_histories``
    histories (used for trajectory-level parity)."""

    track_end_positions: bool = False
    escape_records: bool = False
    trace_histories: int = 0

    def __post_init__(self) -> None:
        for name in ("track_end_positions", "escape_records"):
            if not isinstance(getattr(self, name), bool):
                raise fail(f"{name} must be a bool")
        integer("trace_histories", self.trace_histories, minimum=0)


@dataclass(frozen=True)
class SimulationConfig:
    """Complete simulation request."""

    source: PencilBeamSource
    geometry: VoxelGeometry | BoxPhantom
    scoring: tuple[ScoringGrid, ...]
    physics: PhysicsOptions
    run: RunOptions
    diagnostics: DiagnosticsOptions = DiagnosticsOptions()
    tallies: tuple[TallyRequest, ...] = ()
    lookups: tuple[LookupTable, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.source, PencilBeamSource):
            raise fail("source must be a PencilBeamSource")
        if not isinstance(self.geometry, VoxelGeometry | BoxPhantom):
            raise fail("geometry must be a VoxelGeometry or BoxPhantom")
        if not isinstance(self.scoring, tuple) or not all(
            isinstance(g, ScoringGrid) for g in self.scoring
        ):
            raise fail("scoring must be a tuple of ScoringGrid")
        if not isinstance(self.physics, PhysicsOptions):
            raise fail("physics must be PhysicsOptions")
        if not isinstance(self.run, RunOptions):
            raise fail("run must be RunOptions")
        if not isinstance(self.diagnostics, DiagnosticsOptions):
            raise fail("diagnostics must be DiagnosticsOptions")
        if not isinstance(self.tallies, tuple) or not all(
            isinstance(t, TallyRequest) for t in self.tallies
        ):
            raise fail("tallies must be a tuple of TallyRequest")
        if not isinstance(self.lookups, tuple) or not all(
            isinstance(t, LookupTable) for t in self.lookups
        ):
            raise fail("lookups must be a tuple of LookupTable")


@dataclass(frozen=True, eq=False)
class EffectiveConfig:
    """What will actually run (recorded in every result).

    ``unit_direction`` is the normalised beam direction; ``max_steps`` the per-history step
    bound (``max_steps_origin`` is ``"user"`` or ``"computed"``); ``tables`` the transport
    tables; ``scattering_length_g_cm2`` the Gottschalk ``X_S`` of every material;
    ``production`` is True only for float32 Warp backends; ``scoring_pieces`` the bound of the
    voxel pieces per leg of the track-length scoring; ``rng`` describes the generator;
    ``channels`` the compiled scoring channels (None without tallies).
    """

    requested: SimulationConfig
    geometry: VoxelGeometry
    unit_direction: tuple[float, float, float]
    tables: TransportTables
    max_steps: int
    max_steps_origin: str
    backend: str
    precision: str
    production: bool
    scoring_pieces: int
    scattering_length_g_cm2: tuple[float, ...]
    rng: dict[str, Any]
    channels: ChannelPlan | None = None
    nuclear: NuclearSetup | None = None

    def summary(self) -> dict[str, Any]:
        """JSON-serialisable summary of the effective configuration."""
        p = self.requested.physics
        r = self.requested.run
        s = self.requested.source
        out: dict[str, Any] = {
            "backend": self.backend,
            "precision": self.precision,
            "production": self.production,
            "seed": r.seed,
            "n_histories": r.n_histories,
            "n_batches": r.n_batches,
            "cpu_workers": r.cpu_workers,
            "chunk_histories": r.chunk_histories,
            "max_steps": self.max_steps,
            "scoring_pieces": self.scoring_pieces,
            "max_steps_origin": self.max_steps_origin,
            "unit_direction": list(self.unit_direction),
            "source": {
                "projectile": s.projectile.name,
                "position_mm": list(s.position_mm),
                "kinetic_energy_mev": s.kinetic_energy_mev,
                "energy_sigma_mev": s.energy_sigma_mev,
                "lateral_sigma_mm": s.lateral_sigma_mm,
            },
            "physics": {
                "nuclear": p.nuclear,
                "straggling": p.straggling,
                "straggling_model": p.straggling_model,
                "multiple_scattering": p.multiple_scattering,
                "mcs_model": p.mcs_model,
                "delta_electrons": p.delta_electrons,
                "e_cut_mev": p.e_cut_mev,
                "max_step_mm": p.max_step_mm,
                "max_energy_loss_fraction": p.max_energy_loss_fraction,
                "range_alpha": p.range_alpha,
                "range_rho_f_mm": p.range_rho_f_mm,
                "short_step_fraction": p.short_step_fraction,
                "truncated_hinge_diagnostic": p.truncated_hinge_diagnostic,
                "scattering_E_s_mev": E_S_MEV,
                "max_rejection_attempts": MAX_REJECTION_ATTEMPTS,
            },
            "geometry": {
                "origin_mm": list(self.geometry.origin_mm),
                "spacing_mm": list(self.geometry.spacing_mm),
                "shape": list(self.geometry.shape),
                "materials": [m.name for m in self.geometry.materials],
                "z_exit_mm": self.geometry.z_exit_mm,
            },
            "scoring": [
                {
                    "name": g.name,
                    "origin_mm": list(g.origin_mm),
                    "spacing_mm": list(g.spacing_mm),
                    "shape": list(g.shape),
                }
                for g in self.requested.scoring
            ],
            "tables": {
                "sha256": self.tables.sha256,
                "stopping_source": p.stopping.name,
                "materials": thaw(self.tables.identity),
                "n_e": self.tables.n_e,
                "n_r": self.tables.n_r,
                **(
                    {"water_row": thaw(self.tables.water_identity)} if self.tables.has_water else {}
                ),
            },
            "scattering_length_g_cm2": list(self.scattering_length_g_cm2),
            "rng": dict(self.rng),
            "tallies": None if self.channels is None else self.channels.summary(),
        }
        if self.nuclear is not None:  # nuclear=False summaries are unchanged
            out["physics"]["e_cut_deuteron_mev"] = p.e_cut_deuteron_mev
            out["nuclear"] = self.nuclear.summary(self.geometry)
            if self.nuclear.elastic is not None:  # elastic=False summaries are unchanged
                out["elastic"] = self.nuclear.elastic.summary()
            if self.backend != "python":  # the uploaded arrays of the Warp kernels (V3-005B)
                from ionmc.transport.nuclear_device import host_sha256, pack_nuclear

                out["nuclear"]["device_sha256"] = host_sha256(
                    pack_nuclear(self.nuclear.table, self.geometry.materials), self.precision
                )
                if self.nuclear.elastic is not None:
                    from ionmc.transport import elastic_device as ed

                    el = self.nuclear.elastic
                    out["elastic"]["device_sha256"] = ed.host_sha256(
                        ed.pack_elastic(el.table, self.geometry.materials, rows=el.rows),
                        self.precision,
                    )
        return out


def _resolve_geometry(geometry: VoxelGeometry | BoxPhantom) -> VoxelGeometry:
    return geometry.to_geometry() if isinstance(geometry, BoxPhantom) else geometry


def _check_table_identity(
    table: StoppingTable, material: Material, source_name: str, projectile: Projectile = PROTON
) -> None:
    """Fail closed unless ``table`` describes exactly the requested projectile and material."""
    if not isinstance(table, StoppingTable):
        raise fail(f"source {source_name!r} returned {type(table).__name__}, not a StoppingTable")
    if table.projectile != projectile:
        raise fail(
            f"source {source_name!r} returned a table for projectile {table.projectile.name!r}, "
            f"not the requested {projectile.name}"
        )
    got = table.material
    if got != material:
        raise fail(
            f"source {source_name!r} returned a table for material {got.name!r} "
            f"(density {got.density_g_cm3}, Z/A {got.z_over_a}, I {got.mean_excitation_eV} eV, "
            f"composition {dict(got.mass_fractions)}) that differs from the requested "
            f"{material.name!r} (density {material.density_g_cm3}, Z/A {material.z_over_a}, "
            f"I {material.mean_excitation_eV} eV, composition {dict(material.mass_fractions)}); "
            "the full material definition must be identical"
        )


def _stopping_tables(
    source: StoppingSource, materials: tuple[Material, ...], projectile: Projectile = PROTON
) -> list[StoppingTable]:
    tables: list[StoppingTable] = []
    for material in materials:
        try:
            table = source.table(material, projectile)
        except (ValueError, KeyError) as exc:
            raise fail(
                f"no stopping table for material {material.name!r} from source "
                f"{source.name!r}: {exc}"
            ) from exc
        _check_table_identity(table, material, source.name, projectile)
        tables.append(table)
    return tables


def _computed_max_steps(
    geometry: VoxelGeometry, tables: TransportTables, config: SimulationConfig
) -> int:
    """Conservative bound on the number of steps of one history (recorded, not a physics cut).

    Path length bound ``L`` from the CSDA range at the highest source energy in the least
    dense voxel of each material; steps are bounded by ``L / max_step`` (maximum-step
    steps), ``L (1/dx + 1/dy + 1/dz)`` (voxel crossings), the number of energy-loss-limited
    steps and of range-limited steps, with a safety factor of 4.
    """
    src, ph = config.source, config.physics
    e_hi = min(src.kinetic_energy_mev + 6.0 * src.energy_sigma_mev, float(tables.e_max_mev.min()))
    dens = geometry.densities_g_cm3()
    path_mm = 0.0
    for m in range(len(geometry.materials)):
        mask = geometry.material_index == m
        if not mask.any():
            continue
        rho_min = float(dens[mask].min())
        path_mm = max(path_mm, tables.range_g_cm2(m, e_hi) * 10.0 / rho_min)
    sp = geometry.spacing_mm
    n_len = math.ceil(path_mm / ph.max_step_mm)
    n_cross = math.ceil(path_mm * (1.0 / sp[0] + 1.0 / sp[1] + 1.0 / sp[2])) + 3
    n_eloss = math.ceil(1.5 * math.log(max(e_hi / ph.e_cut_mev, 1.0)) / ph.max_energy_loss_fraction)
    n_range = math.ceil(math.log(max(path_mm / ph.range_rho_f_mm, 1.0)) / ph.range_alpha) + 10
    return int(4 * (n_len + n_cross + n_eloss + n_range) + 100)


def scoring_pieces_bound(max_step_mm: float, grid: ScoringGrid) -> int:
    """Bound ``3 ceil(max_step / spacing_min) + 4`` of the voxel pieces of one leg.

    A leg of length ``L <= max_step`` meets at most ``floor(L / d) + 1 <= ceil(max_step / d) + 1``
    planes of an axis with spacing ``d`` (a start exactly on a plane counts: it is a
    zero-length hop), so at most ``3 ceil(.) + 3`` planes in all and one more piece than
    planes. A leg that exceeds the bound at run time increments ``scoring_pieces_overflow`` and
    invalidates the result."""
    return 3 * math.ceil(max_step_mm / min(grid.spacing_mm)) + 4


def trace_buffer_bytes(trace_histories: int, max_steps: int) -> int:
    """Bytes of the trace buffers (int32 discrete and float64 continuous columns per step)."""
    return trace_histories * max_steps * (4 * TRACE_N_DISCRETE + 8 * TRACE_N_CONTINUOUS)


def _compile_tallies(
    config: SimulationConfig,
    geometry: VoxelGeometry,
    tables: TransportTables,
    e_hi_mev: float,
    max_steps: int,
    scoring_pieces: int,
    water: StoppingTable | None,
    nuclear: NuclearSetup | None = None,
) -> ChannelPlan | None:
    """Fail-closed validation and compilation of the scoring requests (None without tallies)."""
    if not config.tallies:
        if config.lookups:
            raise fail("lookup tables were given but no tally request uses them")
        return None
    assert water is not None
    src, ph, run = config.source, config.physics, config.run
    plan = compile_channels(
        config.tallies,
        config.lookups,
        config.scoring,
        geometry=geometry,
        tables=tables,
        water=water,
        projectile=src.projectile,
        e_cut_mev=ph.e_cut_mev,
        e_hi_mev=e_hi_mev,
        n_histories=run.n_histories,
        n_batches=run.n_batches,
        cpu_workers=run.cpu_workers,
        memory_budget_bytes=run.memory_budget_bytes,
        max_steps=max_steps,
        scoring_pieces=scoring_pieces,
        producible=producible(src.projectile, ph.nuclear),
        max_step_mm=ph.max_step_mm,
        nuclear=None if nuclear is None else nuclear.capacity,
    )
    unused = {lk.name for lk in config.lookups} - {t.lookup for t in config.tallies}
    if unused:
        raise fail(f"lookup tables {sorted(unused)} are not used by any tally request")
    return plan


def _water_table(
    config: SimulationConfig, projectile: Projectile | None = None
) -> StoppingTable | None:
    """The water stopping table of the LET definition from the run's own stopping source (None
    without tallies); fail closed if the source has none. ``projectile`` defaults to the
    source's (a nuclear run also asks for the deuteron)."""
    if not config.tallies:
        return None
    ph = config.physics
    proj = config.source.projectile if projectile is None else projectile
    try:
        water = ph.stopping.table(WATER, proj)
    except (ValueError, KeyError) as exc:
        raise fail(
            f"no water stopping table for the LET definition from source {ph.stopping.name!r}: "
            f"{exc}"
        ) from exc
    _check_table_identity(water, WATER, ph.stopping.name, proj)
    return water


@dataclass(frozen=True, eq=False)
class ElasticSetup:
    """The verified elastic inputs of a run with the hadronic elastic channel (V3-005C): the loaded
    ``table`` (an ``ElasticTable``), its per-material ``rows`` (``MaterialElastic``),
    ``elastic_only`` (non-elastic channel off), the per-material lower domain bounds used by the
    diagnostic counters (``e_min_pa``: largest ``e_min_shape`` over the material's p + A targets;
    ``e_min_pp``: ``E_min,pp`` if the material contains hydrogen, else 0) and whether the table is
    the one pinned in ``elastic_table_pin.json``."""

    table: Any
    rows: tuple[Any, ...]
    elastic_only: bool
    e_min_pa: tuple[float, ...]
    e_min_pp: tuple[float, ...]
    pinned: bool

    def summary(self) -> dict[str, Any]:
        info = self.table.info
        return {
            "table_id": self.table.table_id,
            "pinned": self.pinned,
            "schema": info["schema"],
            "builder_version": info["builder_version"],
            "npz_sha256": info["npz_sha256"],
            "elastic_only": self.elastic_only,
            "elastic_domain": {k: list(v) for k, v in self.table.elastic_domain().items()},
        }


def _pinned_elastic_table_id() -> str:
    import json
    from pathlib import Path

    path = Path(__file__).resolve().parent / "data" / "elastic_table_pin.json"
    return str(json.loads(path.read_text(encoding="utf-8"))["table_id"])


def _elastic_setup(config: SimulationConfig, geometry: VoxelGeometry) -> ElasticSetup:
    """Fail-closed rules of the elastic channel (all backends since C4): a loaded, re-hashed
    and qualified table, every element of every material covered (UnsupportedCombinationError)."""
    from ionmc.nuclear.elastic_tables import ElasticTable

    ph = config.physics
    pin = _pinned_elastic_table_id()
    tid = ph.elastic_table_id if ph.elastic_table_id is not None else pin
    table = ElasticTable.load(None, tid)
    rows = tuple(table.material_rows(m) for m in geometry.materials)
    dom = table.elastic_domain()
    e_pa, e_pp = [], []
    for m in geometry.materials:
        pa = [
            dom[table.info["elements"][sym]["target"]][0]
            for sym in m.mass_fractions
            if table.info["elements"][sym]["target"] != "H-1"
        ]
        e_pa.append(max(pa, default=0.0))
        e_pp.append(dom["H-1"][0] if "H" in m.mass_fractions else 0.0)
    return ElasticSetup(table, rows, ph.elastic_only, tuple(e_pa), tuple(e_pp), tid == pin)


@dataclass(frozen=True, eq=False)
class NuclearSetup:
    """The verified nuclear inputs of a ``nuclear=True`` run (decision 0041 sections 5, 6):
    the loaded ``table`` (a ``NuclearTable``), its per-material ``rows`` (``MaterialNuclear``, one
    per geometry material), the deuteron transport tables (also ``tables.deuteron``),
    ``deuteron_water`` (the deuteron water stopping table, tallies only), the deuteron cutoff,
    the table's per-particle energy bound, the producible set and, with tallies, the
    ``capacity`` inputs of the channel compiler."""

    table: Any
    rows: tuple[Any, ...]
    deuteron_tables: TransportTables
    deuteron_water: StoppingTable | None
    e_cut_deuteron_mev: float
    transport_energy_bound_mev: float
    history_energy_bound_mev: float
    producible: frozenset[tuple[str, str]]
    capacity: NuclearCapacity | None = None
    elastic: ElasticSetup | None = None

    def summary(self, geometry: VoxelGeometry) -> dict[str, Any]:
        """JSON-serialisable record of the nuclear configuration (the effective config)."""
        from ionmc.nuclear.tables import DEFAULT_F_E, MAJORANT_FACTOR

        info = self.table.info
        elements = sorted({e for m in geometry.materials for e in m.mass_fractions})
        return {
            "table_id": self.table.table_id,
            "schema": info["schema"],
            "builder_version": info["builder_version"],
            "npz_sha256": info["npz_sha256"],
            "source_sha256": {k: v["sha256"] for k, v in sorted(info["sources"].items())},
            "elements": {
                e: {
                    "target": info["elements"][e]["target"],
                    "surrogate": bool(info["elements"][e]["surrogate"]),
                    "sigma_scale": float(info["elements"][e]["sigma_scale"]),
                }
                for e in elements
                if e in info["elements"]
            },
            "q_plus_table_mev": float(info["q_plus_table_mev"]),
            "transport_energy_bound_mev": self.transport_energy_bound_mev,
            "history_energy_bound_mev": self.history_energy_bound_mev,
            "recoil_t_max_mev": float(info["recoil_t_max_mev"]),
            "transport_path_bound_terms": thaw(info["transport_path_bound_terms"]),
            "grid_size": int(info["grid"]["n_points"]),
            "multiplicity_model": info["multiplicity"]["model"],
            "majorant_factor": float(MAJORANT_FACTOR),
            "majorant_window_f_e": float(self.rows[0].f_e) if self.rows else DEFAULT_F_E,
            "e_cut_deuteron_mev": self.e_cut_deuteron_mev,
            "deuteron_tables_sha256": self.deuteron_tables.sha256,
            "producible": [{"species": n, "generation": g} for n, g in sorted(self.producible)],
        }


def _nuclear_setup_checks(
    config: SimulationConfig, geometry: VoxelGeometry
) -> tuple[Any, tuple[Any, ...]]:
    """Fail-closed rules of ``nuclear=True`` that need no tables (decision 0041 section 5): any
    ``nist-star`` (no deuteron table), a missing table id, a missing / stale /
    mis-pinned table (the loader's ``NuclearTableError``, an ``UnsupportedCombinationError``),
    ``E0 + 6 sigma_E > 250 MeV`` and an element of a material without a table entry. Returns
    the loaded table and the per-material rows."""
    from ionmc.nuclear.tables import NuclearTable  # lazy: heavy and only for nuclear runs

    src, ph = config.source, config.physics
    if ph.stopping.name == "nist-star":
        raise fail(
            "nuclear=True is not available with the 'nist-star' stopping source: secondary "
            "deuterons need a deuteron stopping table, which NIST PSTAR/ASTAR do not provide; "
            "use an analytic Bethe source"
        )
    if ph.nuclear_table_id is None:
        raise fail("nuclear=True needs PhysicsOptions.nuclear_table_id (the derived nuclear table)")
    table = NuclearTable.load(None, ph.nuclear_table_id)
    e_max = src.kinetic_energy_mev + 6.0 * src.energy_sigma_mev
    if e_max > NUCLEAR_MAX_ENERGY_MEV:
        raise fail(
            f"nuclear=True supports source energies up to {NUCLEAR_MAX_ENERGY_MEV} MeV "
            f"(E0 + 6 sigma_E = {e_max} MeV; the cross sections are extended to 250 MeV only; at "
            "runtime every sampled source energy above the limit is counted in "
            "source_energy_out_of_range and invalidates the result)"
        )
    rows = tuple(table.material_rows(m) for m in geometry.materials)
    return table, rows


def cuda_available() -> bool:
    """True if Warp sees a usable CUDA device (``cuda:0``)."""
    import warp as wp

    try:
        return bool(wp.is_cuda_available())
    except Exception:  # pragma: no cover - driver probing differs between hosts
        return False


def _nuclear_tables(
    config: SimulationConfig,
    geometry: VoxelGeometry,
    tables: TransportTables,
    nuc_table: Any,
    nuc_rows: tuple[Any, ...],
) -> tuple[TransportTables, NuclearSetup]:
    """Deuteron transport tables (``tables.deuteron``), the energy-bound and cutoff rules and the
    capacity inputs of a nuclear run; fail closed if the stopping source has no deuteron table
    or the table's energy bound exceeds a stopping-table maximum."""
    ph = config.physics
    d_stop = _stopping_tables(ph.stopping, geometry.materials, DEUTERON)
    d_water = _water_table(config, DEUTERON)
    d_tables = TransportTables.from_stopping_tables(
        d_stop, water=d_water, projectile=DEUTERON, energy_axis="total_kinetic_mev"
    )
    t_bound = float(nuc_table.info["transport_energy_bound_mev"])
    if t_bound > float(tables.e_max_mev.min()) or t_bound > float(d_tables.e_max_mev.min()):
        raise fail(
            f"the nuclear table's per-particle energy bound {t_bound:g} MeV exceeds the maximum "
            "energy of the proton or deuteron stopping tables"
        )
    d_floor = float(d_tables.e_min_mev.max())
    if d_floor > 0.5 * ph.e_cut_deuteron_mev:
        raise fail(
            f"the deuteron stopping tables start at {d_floor} MeV (total), above half of "
            f"e_cut_deuteron_mev ({ph.e_cut_deuteron_mev}); the deuteron cutoff must be at "
            "least twice the table floor"
        )
    h_bound = float(nuc_table.info["history_energy_bound_mev"])
    terms = {  # the transported species (proton, deuteron) enter the path bound
        k: (float(v["n_max"]), float(v["t_lab_max_mev"]))
        for k, v in nuc_table.info["transport_path_bound_terms"].items()
        if k in ("p", "d")
    }
    capacity = (
        None
        if d_water is None
        else NuclearCapacity(d_tables, d_water, ph.e_cut_deuteron_mev, t_bound, h_bound, terms)
    )
    setup = NuclearSetup(
        table=nuc_table,
        rows=nuc_rows,
        deuteron_tables=d_tables,
        deuteron_water=d_water,
        e_cut_deuteron_mev=ph.e_cut_deuteron_mev,
        transport_energy_bound_mev=t_bound,
        history_energy_bound_mev=h_bound,
        producible=producible(config.source.projectile, True),
        capacity=capacity,
    )
    return replace(tables, deuteron=d_tables), setup


def validate(config: SimulationConfig) -> EffectiveConfig:
    """Check every fail-closed rule and return the effective configuration.

    Raises ``UnsupportedCombinationError`` for an unsupported or invalid request and
    ``BackendUnavailableError`` for a backend that is not available (nothing falls back).
    """
    if not isinstance(config, SimulationConfig):
        raise fail("config must be a SimulationConfig")
    src, ph, run = config.source, config.physics, config.run
    geometry = _resolve_geometry(config.geometry)

    if src.projectile != PROTON:
        raise fail(
            f"projectile {src.projectile!r} is not supported: the transport engine implements "
            f"exactly the canonical proton {PROTON!r} (ions arrive with V3-008)"
        )
    nuc_table = None
    nuc_rows: tuple[Any, ...] = ()
    if ph.nuclear:
        nuc_table, nuc_rows = _nuclear_setup_checks(config, geometry)
    choice("straggling_model", ph.straggling_model, STRAGGLING_MODELS)
    choice("mcs_model", ph.mcs_model, MCS_MODELS)
    choice("delta_electrons", ph.delta_electrons, DELTA_ELECTRON_MODELS)

    if run.backend == "python" and run.precision == "float32":
        raise fail("backend 'python' runs in float64 only; use precision='float64'")
    if run.backend == "warp-cuda" and run.cpu_workers > 1:
        raise fail("cpu_workers > 1 is meaningless on backend 'warp-cuda'")
    if not 0 <= run.seed < 2**64:
        raise fail(f"seed must be in [0, 2**64), got {run.seed}")
    if not 1 <= run.n_histories <= 2**32:
        raise fail(f"n_histories must be in [1, 2**32], got {run.n_histories}")
    if run.n_batches < 2:
        raise fail(f"n_batches must be >= 2 for a variance estimate, got {run.n_batches}")
    if run.n_histories % run.n_batches != 0:
        raise fail(
            f"n_histories ({run.n_histories}) must be a multiple of n_batches ({run.n_batches})"
        )
    if run.cpu_workers > MAX_CPU_WORKERS:
        raise fail(f"cpu_workers must be <= {MAX_CPU_WORKERS}, got {run.cpu_workers}")
    if run.cpu_workers > run.n_histories:
        raise fail(
            f"cpu_workers ({run.cpu_workers}) must not exceed n_histories ({run.n_histories})"
        )

    if not 1 <= len(config.scoring) <= MAX_SCORING_GRIDS:
        raise fail(f"between 1 and {MAX_SCORING_GRIDS} scoring grids are required")
    names = [g.name for g in config.scoring]
    if len(set(names)) != len(names):
        raise fail(f"scoring grid names must be unique, got {names}")
    bytes_per = 8  # int64 fixed-point accumulators on every backend
    accumulator_bytes = (
        run.n_batches * sum(g.n_voxels for g in config.scoring) * bytes_per * run.cpu_workers
    )
    if accumulator_bytes > run.memory_budget_bytes:
        raise fail(
            f"per-batch accumulators need {accumulator_bytes} bytes (including one private copy "
            f"per worker process), above the memory budget of {run.memory_budget_bytes} bytes "
            "(reduce n_batches, the grids or cpu_workers)"
        )
    scoring_pieces = max(scoring_pieces_bound(ph.max_step_mm, g) for g in config.scoring)
    if scoring_pieces > MAX_SCORING_PIECES:
        raise fail(
            f"max_step_mm ({ph.max_step_mm}) against the finest scoring spacing would allow "
            f"{scoring_pieces} voxel pieces per leg in the track-length scoring, above the "
            f"limit of {MAX_SCORING_PIECES}; increase the scoring spacing or reduce max_step_mm"
        )

    stopping_tables = _stopping_tables(ph.stopping, geometry.materials)
    water_table = _water_table(config)
    tables = TransportTables.from_stopping_tables(stopping_tables, water=water_table)
    nuclear: NuclearSetup | None = None
    if nuc_table is not None:
        tables, nuclear = _nuclear_tables(config, geometry, tables, nuc_table, nuc_rows)
        if ph.elastic:
            nuclear = replace(nuclear, elastic=_elastic_setup(config, geometry))
    e_lo = float(tables.e_min_mev.max())
    e_hi_table = float(tables.e_max_mev.min())
    if e_lo > 0.5 * ph.e_cut_mev:
        raise fail(
            f"the stopping tables start at {e_lo} MeV, above half of e_cut_mev "
            f"({ph.e_cut_mev}); the cutoff energy must be at least twice the table floor"
        )
    e0 = src.kinetic_energy_mev
    if not 2.0 * ph.e_cut_mev <= e0 <= e_hi_table:
        raise fail(
            f"source energy {e0} MeV outside [2 e_cut, table maximum] = "
            f"[{2.0 * ph.e_cut_mev}, {e_hi_table}] MeV"
        )
    if src.energy_sigma_mev > MAX_ENERGY_SIGMA_FRACTION * e0:
        raise fail(
            f"energy_sigma_mev ({src.energy_sigma_mev}) exceeds "
            f"{MAX_ENERGY_SIGMA_FRACTION} of the source energy"
        )

    e_hi = min(e0 + 6.0 * src.energy_sigma_mev, e_hi_table)
    # the energy bound of a history: a nuclear history may deposit more than its initial energy
    # (decision 0041 section 5, amended 2026-10-08): max(E_hi, the table's per-history bound)
    e_cap = e_hi if nuclear is None else max(e_hi, nuclear.history_energy_bound_mev)
    per_batch = run.n_histories // run.n_batches
    if per_batch * e_cap / QUANTUM_MEV >= MAX_QUANTA:
        raise fail(
            f"{per_batch} histories per batch of up to {e_cap:g} MeV exceed the capacity "
            f"{MAX_QUANTA * QUANTUM_MEV:g} MeV of a fixed-point voxel accumulator (quantum "
            f"{QUANTUM_MEV:g} MeV); increase n_batches or reduce n_histories"
        )

    max_steps_origin = "user" if run.max_steps is not None else "computed"
    max_steps = (
        run.max_steps
        if run.max_steps is not None
        else _computed_max_steps(geometry, tables, config)
    )
    if max_steps * (1 + MAX_REJECTION_ATTEMPTS) >= 2**32:
        raise fail(
            f"max_steps ({max_steps}) times {1 + MAX_REJECTION_ATTEMPTS} draw blocks per step "
            "reaches the 2**32 block-counter bound"
        )

    diag = config.diagnostics
    if diag.trace_histories > 0:
        if run.precision != "float64":
            raise fail(
                "the per-step trace is only available in float64 (trace_histories > 0 with "
                f"precision {run.precision!r}); the trace is the trajectory-parity reference"
            )
        if diag.trace_histories > run.n_histories:
            raise fail(
                f"trace_histories ({diag.trace_histories}) exceeds n_histories ({run.n_histories})"
            )
        trace_bytes = trace_buffer_bytes(diag.trace_histories, max_steps)
        if trace_bytes > MAX_TRACE_BUFFER_BYTES:
            raise fail(
                f"the trace buffers ({diag.trace_histories} histories x {max_steps} steps) need "
                f"{trace_bytes} bytes, above the limit of {MAX_TRACE_BUFFER_BYTES} bytes"
            )

    channels = _compile_tallies(
        config, geometry, tables, e_hi, max_steps, scoring_pieces, water_table, nuclear
    )
    if channels is not None and run.backend not in CHANNEL_BACKENDS:
        raise fail(
            f"backend {run.backend!r} does not implement scoring channels yet "
            f"(supported: {list(CHANNEL_BACKENDS)}); tallies are never silently ignored"
        )

    if run.backend == "warp-cuda" and not cuda_available():
        raise BackendUnavailableError(
            "backend 'warp-cuda' requires a CUDA device (cuda:0) that Warp can use; none is "
            "available and nothing falls back to another backend"
        )

    return EffectiveConfig(
        requested=config,
        geometry=geometry,
        unit_direction=src.unit_direction(),
        tables=tables,
        max_steps=max_steps,
        max_steps_origin=max_steps_origin,
        backend=run.backend,
        precision=run.precision,
        production=run.backend != "python" and run.precision == "float32",
        scoring_pieces=scoring_pieces,
        scattering_length_g_cm2=tuple(scattering_length_g_cm2(m) for m in geometry.materials),
        channels=channels,
        nuclear=nuclear,
        rng={
            "generator": "philox4x32-10",
            "key": "seed low/high 32 bits",
            "counter": "(history, genealogy id, block, purpose)",
            "u01": U01_MAPPING[run.precision],
        },
    )
