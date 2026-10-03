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
from dataclasses import dataclass
from typing import Any, Literal

from ionmc._validate import choice, fail, integer, real
from ionmc.errors import BackendUnavailableError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.materials import Material
from ionmc.physics.em import E_S_MEV
from ionmc.physics.projectiles import PROTON
from ionmc.physics.scattering import scattering_length_g_cm2
from ionmc.physics.stopping import StoppingSource, StoppingTable
from ionmc.scoring import MAX_SCORING_GRIDS, ScoringGrid
from ionmc.sources import PencilBeamSource
from ionmc.transport.tables import TransportTables, thaw

STRAGGLING_MODELS = ("bohr_gauss_clamped_gamma_v1",)
MCS_MODELS = ("differential_moliere",)
DELTA_ELECTRON_MODELS = ("local",)
BACKENDS = ("python", "warp-cpu", "warp-cuda")
PRECISIONS = ("float32", "float64")
MAX_REJECTION_ATTEMPTS = 64
MAX_ENERGY_SIGMA_FRACTION = 0.05
MAX_ENERGY_LOSS_FRACTION = 0.2
U01_MAPPING = {
    "float64": "((w >> 8) + 0.5) * 2**-24",
    "float32": "((w >> 9) + 0.5) * 2**-23",
}


@dataclass(frozen=True)
class PhysicsOptions:
    """Physics switches and numerical parameters.

    ``nuclear`` has no default (the caller must state it); ``True`` is rejected until nuclear
    interactions exist. ``stopping`` supplies the electronic stopping tables (for example an
    analytic :class:`~ionmc.physics.stopping.BetheStoppingSource`). Below ``e_cut_mev`` the
    remaining kinetic energy is deposited locally; ``max_step_mm`` bounds the step length;
    ``max_energy_loss_fraction`` bounds the mean energy loss per step as a fraction of the
    kinetic energy; ``range_alpha`` and ``range_rho_f_mm`` parametrise the Geant4 range step
    function; steps shorter than ``short_step_fraction`` of the residual range use the
    linear loss ``S t``.
    """

    nuclear: bool
    stopping: StoppingSource
    straggling: bool = True
    straggling_model: str = "bohr_gauss_clamped_gamma_v1"
    multiple_scattering: bool = True
    mcs_model: str = "differential_moliere"
    delta_electrons: str = "local"
    e_cut_mev: float = 2.0
    max_step_mm: float = 1.0
    max_energy_loss_fraction: float = 0.02
    range_alpha: float = 0.2
    range_rho_f_mm: float = 0.1
    short_step_fraction: float = 1.0e-3

    def __post_init__(self) -> None:
        for name in ("nuclear", "straggling", "multiple_scattering"):
            if not isinstance(getattr(self, name), bool):
                raise fail(f"{name} must be a bool, got {getattr(self, name)!r}")
        for name in ("straggling_model", "mcs_model", "delta_electrons"):
            if not isinstance(getattr(self, name), str):
                raise fail(f"{name} must be a string")
        if not (hasattr(self.stopping, "table") and hasattr(self.stopping, "name")):
            raise fail("stopping must be a StoppingSource (with .name and .table())")
        real("e_cut_mev", self.e_cut_mev, positive=True)
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
    ``memory_budget_bytes`` bounds the per-batch accumulators. ``cpu_workers`` (> 1) and
    ``worker_timeout_s`` belong to the multi-process Warp CPU driver (V3-003B); the python
    backend runs in-process and rejects ``cpu_workers > 1``."""

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


@dataclass(frozen=True, eq=False)
class EffectiveConfig:
    """What will actually run (recorded in every result).

    ``unit_direction`` is the normalised beam direction; ``max_steps`` the per-history step
    bound (``max_steps_origin`` is ``"user"`` or ``"computed"``); ``tables`` the transport
    tables; ``scattering_length_g_cm2`` the Gottschalk ``X_S`` of every material;
    ``production`` is True only for float32 Warp backends; ``rng`` describes the generator.
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
    scattering_length_g_cm2: tuple[float, ...]
    rng: dict[str, Any]

    def summary(self) -> dict[str, Any]:
        """JSON-serialisable summary of the effective configuration."""
        p = self.requested.physics
        r = self.requested.run
        s = self.requested.source
        return {
            "backend": self.backend,
            "precision": self.precision,
            "production": self.production,
            "seed": r.seed,
            "n_histories": r.n_histories,
            "n_batches": r.n_batches,
            "cpu_workers": r.cpu_workers,
            "max_steps": self.max_steps,
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
                "scattering_E_s_mev": E_S_MEV,
                "max_rejection_attempts": MAX_REJECTION_ATTEMPTS,
            },
            "geometry": {
                "origin_mm": list(self.geometry.origin_mm),
                "spacing_mm": list(self.geometry.spacing_mm),
                "shape": list(self.geometry.shape),
                "materials": [m.name for m in self.geometry.materials],
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
            },
            "scattering_length_g_cm2": list(self.scattering_length_g_cm2),
            "rng": dict(self.rng),
        }


def _resolve_geometry(geometry: VoxelGeometry | BoxPhantom) -> VoxelGeometry:
    return geometry.to_geometry() if isinstance(geometry, BoxPhantom) else geometry


def _check_table_identity(table: StoppingTable, material: Material, source_name: str) -> None:
    """Fail closed unless ``table`` describes exactly the requested proton and material."""
    if not isinstance(table, StoppingTable):
        raise fail(f"source {source_name!r} returned {type(table).__name__}, not a StoppingTable")
    if table.projectile != PROTON:
        raise fail(
            f"source {source_name!r} returned a table for projectile {table.projectile.name!r}, "
            f"not the requested proton"
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
    source: StoppingSource, materials: tuple[Material, ...]
) -> list[StoppingTable]:
    tables: list[StoppingTable] = []
    for material in materials:
        try:
            table = source.table(material, PROTON)
        except (ValueError, KeyError) as exc:
            raise fail(
                f"no stopping table for material {material.name!r} from source "
                f"{source.name!r}: {exc}"
            ) from exc
        _check_table_identity(table, material, source.name)
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
    if ph.nuclear:
        raise fail(
            "nuclear=True is not supported: nuclear interactions are not implemented yet "
            "(V3-005); pass nuclear=False explicitly for an electromagnetic-only run"
        )
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
    if run.cpu_workers > 1 and run.backend == "python":
        raise BackendUnavailableError(
            "multi-process execution (cpu_workers > 1) is implemented in V3-003B"
        )

    if not 1 <= len(config.scoring) <= MAX_SCORING_GRIDS:
        raise fail(f"between 1 and {MAX_SCORING_GRIDS} scoring grids are required")
    names = [g.name for g in config.scoring]
    if len(set(names)) != len(names):
        raise fail(f"scoring grid names must be unique, got {names}")
    bytes_per = 4 if run.precision == "float32" else 8
    accumulator_bytes = run.n_batches * sum(g.n_voxels for g in config.scoring) * bytes_per
    if accumulator_bytes > run.memory_budget_bytes:
        raise fail(
            f"per-batch accumulators need {accumulator_bytes} bytes, above the memory budget "
            f"of {run.memory_budget_bytes} bytes (reduce n_batches or the grids)"
        )
    min_spacing = min(min(g.spacing_mm) for g in config.scoring)
    if ph.max_step_mm > min_spacing:
        raise fail(
            f"max_step_mm ({ph.max_step_mm}) exceeds the smallest scoring spacing "
            f"({min_spacing} mm): the midpoint deposit needs steps no longer than one scoring "
            "voxel; reduce max_step_mm (nothing is clamped silently)"
        )

    stopping_tables = _stopping_tables(ph.stopping, geometry.materials)
    tables = TransportTables.from_stopping_tables(stopping_tables)
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

    if run.backend in ("warp-cpu", "warp-cuda"):
        raise BackendUnavailableError(
            f"backend {run.backend!r} is implemented in V3-003B; only backend 'python' is "
            "available (no fallback is performed)"
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
        scattering_length_g_cm2=tuple(scattering_length_g_cm2(m) for m in geometry.materials),
        rng={
            "generator": "philox4x32-10",
            "key": "seed low/high 32 bits",
            "counter": "(history, genealogy id, block, purpose)",
            "u01": U01_MAPPING[run.precision],
        },
    )
