"""Transport engine: reference Python and Warp execution of CSDA depth dose.

The engine runs a monoenergetic proton pencil beam through a water slab or a 1-D
voxelized, materially-heterogeneous phantom (decisions 0014, 0015) and returns
the integral depth dose (energy deposited per depth bin). The
**reference Python** path is a scalar per-history loop calling the float64
shared-source step functions (:mod:`ionmc.physics.transport`); the **warp** path
launches the equivalent kernel on the CPU or CUDA device (decision ``0009``).
Both paths use the same step algorithm and the same tabulated stopping power,
so the reference path is a faithful oracle for the accelerated one.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from ionmc.backend import mathlib, reference
from ionmc.constants import AVOGADRO
from ionmc.data.stopping_tables import StoppingTable
from ionmc.materials import WATER
from ionmc.particles import PROTON, Particle
from ionmc.physics.secondaries import (
    F_HEAVY as SECONDARY_F_HEAVY,
)
from ionmc.physics.secondaries import (
    F_SECONDARY as SECONDARY_F_PROTON,
)
from ionmc.physics.secondaries import (
    generate_secondaries,
)
from ionmc.rng import RandomState
from ionmc.stopping_power import mass_stopping_power_ratio
from ionmc.transport.depth_dose import DepthDoseGrid, DepthLateralGrid
from ionmc.transport.geometry import VoxelSlab, WaterSlab
from ionmc.transport.source import PencilBeamSource
from ionmc.transport.state import ParticleState, Status

#: Default fractional energy loss per step and the absolute step cap [mm].
DEFAULT_MAX_FRACTION: float = 0.02
DEFAULT_MAX_STEP_MM: float = 1.0
#: Default energy cutoff [MeV]; below it the residual energy is deposited locally.
DEFAULT_ENERGY_CUT_MEV: float = 0.5
#: Hard per-history step-count cap (a truncated history is counted, not dropped).
DEFAULT_MAX_STEPS: int = 100_000
#: Below this energy [MeV] energy-loss straggling is switched off and the step is
#: deterministic, so the terminal-step Gaussian clamp cannot bias the range; the
#: residual range below it (<~0.1 mm in water) is negligible (decision 0010).
DEFAULT_STRAGGLING_FLOOR_MEV: float = 2.0
#: Fraction of a reacting primary's kinetic energy deposited locally at the
#: nonelastic reaction vertex (short-range recoils and heavy fragments); the
#: remainder is booked to the escaping/deferred channel for secondary transport
#: in Stage 2 (decision 0012).
DEFAULT_NUCLEAR_LOCAL_FRACTION: float = 0.30


def _merge_voxels(
    z_boundaries: np.ndarray, *per_voxel: np.ndarray
) -> tuple[np.ndarray, ...]:
    """Collapse consecutive voxels sharing **all** their per-voxel physics
    quantities (decisions 0014, 0015, 0017): the water-equivalent density, the
    oxygen-equivalent nuclear density, the physical density and the radiation
    length. Merging on the full set keeps every returned array aligned and safe.

    Returns ``(z_boundaries, *merged_arrays)`` in the same order as the inputs. A
    uniform single-material slab (or a homogeneous ``WaterSlab``) collapses to a
    single voxel, so the voxelized transport reproduces the homogeneous transport
    exactly; step limiting then only happens at genuine material/density
    interfaces.
    """
    arrays = [np.ascontiguousarray(a, dtype=np.float64) for a in per_voxel]
    z = np.ascontiguousarray(z_boundaries, dtype=np.float64)
    merged = [[float(a[0])] for a in arrays]
    merged_z = [float(z[0]), float(z[1])]
    for v in range(1, arrays[0].shape[0]):
        if all(a[v] == m[-1] for a, m in zip(arrays, merged, strict=True)):
            merged_z[-1] = float(z[v + 1])  # extend the current merged voxel
        else:
            for a, m in zip(arrays, merged, strict=True):
                m.append(float(a[v]))
            merged_z.append(float(z[v + 1]))
    return (
        np.asarray(merged_z, dtype=np.float64),
        *(np.asarray(m, dtype=np.float64) for m in merged),
    )


@dataclass(frozen=True)
class DepthDoseResult:
    """Outcome of a depth-dose transport run."""

    edep_mev: np.ndarray  # energy deposited per depth bin [MeV]
    grid: DepthDoseGrid
    n_histories: int
    energy_in_mev: float  # total kinetic energy started
    truncated: int  # histories that hit the step cap
    path: str
    device: str | None
    range_mean_mm: float = float("nan")  # mean stopping depth of stopped histories
    range_sigma_mm: float = float("nan")  # range straggling (std of stopping depth)
    n_stopped: int = 0
    #: Kinetic energy [MeV] carried off by nonelastic reaction products (secondary
    #: protons, neutrons, gammas) and booked for later secondary transport
    #: (decision 0012). Zero unless the run had ``nuclear=True``.
    escaped_mev: float = 0.0
    #: Number of primaries removed by a nonelastic nuclear reaction.
    n_reactions: int = 0
    #: Per-depth-bin dose from transported secondary protons and their sub-cut
    #: local deposits [MeV]; included in ``edep_mev``. Zero unless the run had
    #: ``secondaries=True`` (decision 0013).
    secondary_edep_mev: np.ndarray | None = None
    #: Number of secondary protons transported.
    n_secondaries: int = 0

    @property
    def energy_deposited_mev(self) -> float:
        return float(np.sum(self.edep_mev))

    @property
    def secondary_dose_fraction(self) -> np.ndarray:
        """Per-bin fraction of the deposited dose that is of secondary origin.

        Zero where no secondary array was recorded or the total dose is zero.
        """
        out = np.zeros_like(self.edep_mev)
        if self.secondary_edep_mev is None:
            return out
        nz = self.edep_mev > 0.0
        out[nz] = self.secondary_edep_mev[nz] / self.edep_mev[nz]
        return out

    @property
    def energy_balance(self) -> float:
        """Relative energy imbalance ``(deposited + escaped - in) / in``.

        With ``nuclear=True`` a reacting primary's energy is split between the
        locally deposited fraction and the escaping channel; both must be
        counted for the balance to close (decision 0012).
        """
        accounted = self.energy_deposited_mev + self.escaped_mev
        return (accounted - self.energy_in_mev) / self.energy_in_mev

    def r80_mm(self) -> float:
        """Distal 80 %-of-maximum depth of the integral depth dose [mm]."""
        return _distal_level_depth(self.edep_mev, self.grid, 0.80)


@dataclass(frozen=True)
class BatchedDepthDoseResult:
    """Mean depth dose and its statistical uncertainty over independent batches."""

    mean_edep_mev: np.ndarray
    standard_error_mev: np.ndarray
    grid: DepthDoseGrid
    n_batches: int
    histories_per_batch: int
    truncated: int
    path: str
    device: str | None

    @property
    def relative_standard_error(self) -> np.ndarray:
        """Per-bin relative standard error where the mean is positive."""
        out = np.zeros_like(self.mean_edep_mev)
        nz = self.mean_edep_mev > 0.0
        out[nz] = self.standard_error_mev[nz] / self.mean_edep_mev[nz]
        return out


@dataclass(frozen=True)
class ScatteringResult:
    """Outcome of a 3-D transport run with multiple Coulomb scattering."""

    edep_zx_mev: np.ndarray  # (n_depth, n_lateral) energy [MeV]
    grid: DepthLateralGrid
    n_histories: int
    energy_in_mev: float
    truncated: int
    path: str
    device: str | None
    range_mean_mm: float = float("nan")  # mean projected stopping depth
    n_stopped: int = 0

    @property
    def depth_dose_mev(self) -> np.ndarray:
        return self.grid.depth_dose(self.edep_zx_mev)

    @property
    def sigma_x_mm(self) -> np.ndarray:
        return self.grid.sigma_x_mm(self.edep_zx_mev)

    @property
    def energy_deposited_mev(self) -> float:
        return float(np.sum(self.edep_zx_mev))

    @property
    def energy_balance(self) -> float:
        return (self.energy_deposited_mev - self.energy_in_mev) / self.energy_in_mev

    def sigma_x_at_depth(self, depth_mm: float) -> float:
        """Lateral RMS ``sigma_x`` [mm] interpolated at ``depth_mm``."""
        centers = self.grid.depth_centers_mm
        sigma = self.sigma_x_mm
        good = np.isfinite(sigma)
        return float(np.interp(depth_mm, centers[good], sigma[good]))


def _scatter_direction(
    dx: float, dy: float, dz: float, theta_x: float, theta_y: float
) -> tuple[float, float, float]:
    """Tilt a unit direction by small projected angles in its transverse frame.

    Builds an orthonormal frame perpendicular to ``d`` (reference axis = the
    least-aligned world axis, avoiding degeneracy), applies
    ``d' = normalize(d + theta_x e1 + theta_y e2)`` (the two-plane sampler, space
    angle variance ``2 theta0^2``). Mirrors the Warp kernel exactly; plain float
    arithmetic keeps the two paths numerically identical (decision 0011).
    """
    ax, ay, az = abs(dx), abs(dy), abs(dz)
    if ax <= ay and ax <= az:
        rx, ry, rz = 1.0, 0.0, 0.0
    elif ay <= az:
        rx, ry, rz = 0.0, 1.0, 0.0
    else:
        rx, ry, rz = 0.0, 0.0, 1.0
    c1x, c1y, c1z = dy * rz - dz * ry, dz * rx - dx * rz, dx * ry - dy * rx
    inv1 = 1.0 / math.sqrt(c1x * c1x + c1y * c1y + c1z * c1z)
    e1x, e1y, e1z = c1x * inv1, c1y * inv1, c1z * inv1
    e2x, e2y, e2z = (
        dy * e1z - dz * e1y,
        dz * e1x - dx * e1z,
        dx * e1y - dy * e1x,
    )
    nx = dx + theta_x * e1x + theta_y * e2x
    ny = dy + theta_x * e1y + theta_y * e2y
    nz = dz + theta_x * e1z + theta_y * e2z
    inv = 1.0 / math.sqrt(nx * nx + ny * ny + nz * nz)
    return nx * inv, ny * inv, nz * inv


def _deposit_zx(
    edep: np.ndarray,
    z0: float,
    z1: float,
    x_mid: float,
    energy: float,
    dz: float,
    half_width: float,
    dx_bin: float,
    n_depth: int,
    n_lateral: int,
) -> None:
    """Deposit ``energy`` across the depth bins spanned by ``[z0, z1]`` (by depth
    overlap) at the lateral bin containing ``x_mid`` (decision 0011).

    Mirrors the Warp kernel deposition. Energy outside the grid is dropped.
    """
    xb = math.floor((x_mid + half_width) / dx_bin)
    if not (0 <= xb < n_lateral):
        return
    span = z1 - z0
    if span <= 0.0:
        b = math.floor(z0 / dz)
        if 0 <= b < n_depth:
            edep[b, xb] += energy
        return
    inv = 1.0 / span
    pos = z0
    b = math.floor(z0 / dz)
    while pos < z1 - 1.0e-12:
        bin_end = (b + 1) * dz
        seg_end = min(bin_end, z1)
        if 0 <= b < n_depth:
            edep[b, xb] += energy * (seg_end - pos) * inv
        pos = seg_end
        b += 1


class TransportEngine:
    """Runs CSDA depth-dose transport for one stopping-power table and geometry."""

    def __init__(
        self,
        table: StoppingTable,
        slab: WaterSlab | VoxelSlab,
        grid: DepthDoseGrid,
        max_fraction: float = DEFAULT_MAX_FRACTION,
        max_step_mm: float = DEFAULT_MAX_STEP_MM,
        energy_cut_mev: float = DEFAULT_ENERGY_CUT_MEV,
        max_steps: int = DEFAULT_MAX_STEPS,
        particle: Particle = PROTON,
        straggling: bool = True,
        straggling_floor_mev: float = DEFAULT_STRAGGLING_FLOOR_MEV,
        nuclear: bool = False,
        nuclear_local_fraction: float = DEFAULT_NUCLEAR_LOCAL_FRACTION,
        secondaries: bool = False,
        secondary_heavy_fraction: float = SECONDARY_F_HEAVY,
        secondary_proton_fraction: float = SECONDARY_F_PROTON,
        material_reference_energy_mev: float = 150.0,
    ) -> None:
        if not (0.0 < max_fraction < 1.0):
            raise ValueError("max_fraction must be in (0, 1)")
        if not (0.0 <= nuclear_local_fraction <= 1.0):
            raise ValueError("nuclear_local_fraction must be in [0, 1]")
        if secondaries and not nuclear:
            raise ValueError("secondaries=True requires nuclear=True (decision 0013)")
        if secondary_heavy_fraction + secondary_proton_fraction > 1.0:
            raise ValueError("heavy + proton secondary fractions must be <= 1")
        self.table = table
        self.slab = slab
        self.grid = grid
        self.max_fraction = max_fraction
        self.max_step_mm = max_step_mm
        self.energy_cut_mev = max(energy_cut_mev, float(table.energy_mev[0]))
        self.max_steps = max_steps
        self.particle = particle
        self.straggling = straggling
        self.straggling_floor_mev = straggling_floor_mev
        self.nuclear = nuclear
        self.nuclear_local_fraction = nuclear_local_fraction
        self.secondaries = secondaries
        self.secondary_heavy_fraction = secondary_heavy_fraction
        self.secondary_proton_fraction = secondary_proton_fraction
        self.material_reference_energy_mev = material_reference_energy_mev
        #: Both transport paths operate in a *water-equivalent* frame (water table
        #: x water-equivalent density; decisions 0015, 0017), so the Bohr
        #: straggling prefactor <Z/A> is water's on every path; the material's
        #: <Z/A> enters through the stopping-power ratio.
        self.depth_dose_za_ratio = WATER.electrons_per_gram_ratio
        #: Radiation length [g/cm^2] of the front material (MCS; decision 0011).
        self.radiation_length_g_per_cm2 = slab.material.radiation_length_g_per_cm2
        #: Voxel geometry along the beam axis (decisions 0014, 0015). Each voxel
        #: is transported as water at its **water-equivalent density** rho_we =
        #: SPR(material) x rho_phys, with a composition-scaled oxygen-equivalent
        #: nuclear density. Consecutive voxels sharing both transported quantities
        #: are merged, so a uniform single-material slab (and a homogeneous
        #: WaterSlab) collapse to a single voxel and reproduce the homogeneous
        #: transport exactly.
        raw_z, raw_rho = slab.voxel_profile()
        raw_mats = slab.materials_profile()
        n_phys = int(raw_rho.shape[0])
        we_density = np.empty(n_phys, dtype=np.float64)
        ox_density = np.empty(n_phys, dtype=np.float64)
        phys_density = np.ascontiguousarray(raw_rho, dtype=np.float64)
        radlen = np.empty(n_phys, dtype=np.float64)
        spr_cache: dict[int, float] = {}
        for v in range(n_phys):
            mat_v = raw_mats[v]
            key = id(mat_v)
            if key not in spr_cache:
                spr_cache[key] = mass_stopping_power_ratio(
                    mat_v, self.material_reference_energy_mev, particle
                )
            we_density[v] = spr_cache[key] * raw_rho[v]
            ox_density[v] = AVOGADRO * mat_v.oxygen_equivalent_per_gram * raw_rho[v]
            radlen[v] = mat_v.radiation_length_g_per_cm2
        (
            self.voxel_z_mm,
            self.voxel_density,
            self.voxel_oxygen_density,
            self.voxel_physical_density,
            self.voxel_radiation_length,
        ) = _merge_voxels(raw_z, we_density, ox_density, phys_density, radlen)
        self.n_voxels = int(self.voxel_density.shape[0])
        self.depth_mm = float(self.voxel_z_mm[-1])
        #: Physical front-voxel density for scalar consumers.
        self.density_g_per_cm3 = float(raw_rho[0])
        #: Lazily-built engine that transports secondary protons (nuclear off, no
        #: further secondaries), reusing this engine's geometry (decision 0013).
        self._sec_engine: TransportEngine | None = None

    def run(
        self,
        source: PencilBeamSource,
        n_histories: int = 1,
        seed: int = 12345,
        path: str = "python",
        device: str = "cpu",
    ) -> DepthDoseResult:
        state = source.sample(n_histories, seed)
        energy_in = float(np.sum(state.energy_mev * state.weight))
        (
            edep,
            truncated,
            final_z,
            final_status,
            escaped,
            n_reactions,
            react_z,
            react_e,
        ) = self._transport(state, path, device)
        # second pass: transport the secondary protons produced at the reaction
        # vertices and fold their dose into the grid (decision 0013).
        secondary_edep: np.ndarray | None = None
        n_secondaries = 0
        if self.secondaries and n_reactions > 0:
            edep, secondary_edep, escaped, n_secondaries, sec_truncated = (
                self._transport_secondaries(
                    react_z, react_e, state.weight, edep, escaped, seed, path, device
                )
            )
            truncated += sec_truncated
        stopped = final_z[final_status == int(Status.STOPPED)]
        n_stopped = int(stopped.shape[0])
        range_mean = float(np.mean(stopped)) if n_stopped else float("nan")
        range_sigma = float(np.std(stopped, ddof=1)) if n_stopped > 1 else float("nan")
        return DepthDoseResult(
            edep_mev=edep,
            grid=self.grid,
            n_histories=n_histories,
            energy_in_mev=energy_in,
            truncated=truncated,
            path=path,
            device=device if path == "warp" else None,
            range_mean_mm=range_mean,
            range_sigma_mm=range_sigma,
            n_stopped=n_stopped,
            escaped_mev=escaped,
            n_reactions=n_reactions,
            secondary_edep_mev=secondary_edep,
            n_secondaries=n_secondaries,
        )

    def _transport(
        self, state: ParticleState, path: str, device: str
    ) -> tuple[
        np.ndarray, int, np.ndarray, np.ndarray, float, int, np.ndarray, np.ndarray
    ]:
        """Dispatch one transport pass to the reference or Warp driver."""
        if path == "warp":
            return self._run_warp(state, device)
        if path == "python":
            return self._run_reference(state)
        raise ValueError(f"unknown transport path {path!r} (use 'python' or 'warp')")

    def _secondary_engine(self) -> TransportEngine:
        """Engine that transports secondary protons: same geometry, nuclear off,
        no further secondaries (decision 0013)."""
        if self._sec_engine is None:
            self._sec_engine = TransportEngine(
                self.table,
                self.slab,
                self.grid,
                max_fraction=self.max_fraction,
                max_step_mm=self.max_step_mm,
                energy_cut_mev=self.energy_cut_mev,
                max_steps=self.max_steps,
                particle=self.particle,
                straggling=self.straggling,
                straggling_floor_mev=self.straggling_floor_mev,
                nuclear=False,
                secondaries=False,
            )
        return self._sec_engine

    def _transport_secondaries(
        self,
        react_z: np.ndarray,
        react_e: np.ndarray,
        weight: np.ndarray,
        primary_edep: np.ndarray,
        escaped: float,
        seed: int,
        path: str,
        device: str,
    ) -> tuple[np.ndarray, np.ndarray, float, int, int]:
        """Generate and transport secondary protons; return the combined dose,
        the secondary-only dose, the updated escaping energy, the secondary
        count, and any secondary truncations (decision 0013)."""
        batch = generate_secondaries(
            react_z,
            react_e,
            weight,
            seed=seed,
            f_heavy=self.secondary_heavy_fraction,
            f_secondary=self.secondary_proton_fraction,
        )
        dz = self.grid.bin_width_mm
        n_bins = self.grid.n_bins
        secondary_edep = self.grid.empty()
        # short-range (sub-cut) secondaries deposit locally at their vertex
        for z, weighted_energy in batch.local_deposit:
            b = math.floor(z / dz)
            if 0 <= b < n_bins:
                secondary_edep[b] += weighted_energy
        sec_truncated = 0
        if batch.state is not None:
            sec_engine = self._secondary_engine()
            s_edep, s_trunc, _, _, _, _, _, _ = sec_engine._transport(
                batch.state, path, device
            )
            secondary_edep = secondary_edep + s_edep
            sec_truncated += s_trunc
        # Every MeV a secondary deposits was pulled from the primary's escaping
        # channel; whatever the secondaries do not deposit (sub-cut energy or
        # transported energy that leaves the grid/geometry) simply stays in the
        # escaping channel. So subtracting the deposited secondary dose closes
        # deposited + escaped = energy_in exactly for any geometry (decision 0013;
        # the remaining escaping energy is the neutron/gamma/binding fraction plus
        # any secondary that left the scored region).
        escaped -= float(np.sum(secondary_edep))
        return (
            primary_edep + secondary_edep,
            secondary_edep,
            escaped,
            batch.n_secondaries,
            sec_truncated,
        )

    def run_batched(
        self,
        source: PencilBeamSource,
        n_histories: int,
        n_batches: int = 10,
        seed: int = 12345,
        path: str = "warp",
        device: str = "cpu",
    ) -> BatchedDepthDoseResult:
        """Run ``n_batches`` independent batches and estimate the statistical
        uncertainty of the mean depth dose (decision 0010).

        Each batch uses a distinct seed so the batches are independent Monte
        Carlo samples. The per-bin standard error of the mean is the standard
        deviation of the batch profiles divided by ``sqrt(n_batches)``.
        """
        if n_batches < 2:
            raise ValueError("need at least two batches for an uncertainty estimate")
        per_batch = max(1, n_histories // n_batches)
        profiles = np.empty((n_batches, self.grid.n_bins), dtype=np.float64)
        truncated = 0
        for b in range(n_batches):
            res = self.run(
                source, per_batch, seed=seed + 1 + b, path=path, device=device
            )
            profiles[b] = res.edep_mev
            truncated += res.truncated
        mean = profiles.mean(axis=0)
        standard_error = profiles.std(axis=0, ddof=1) / math.sqrt(n_batches)
        return BatchedDepthDoseResult(
            mean_edep_mev=mean,
            standard_error_mev=standard_error,
            grid=self.grid,
            n_batches=n_batches,
            histories_per_batch=per_batch,
            truncated=truncated,
            path=path,
            device=device if path == "warp" else None,
        )

    def run_scattering(
        self,
        source: PencilBeamSource,
        grid: DepthLateralGrid,
        n_histories: int = 1,
        seed: int = 12345,
        path: str = "warp",
        device: str = "cpu",
    ) -> ScatteringResult:
        """Run 3-D transport with multiple Coulomb scattering into a 2-D
        (depth, lateral-x) grid (decision 0011).

        Scattering is always applied (above the 2 MeV floor); energy-loss
        straggling independently follows the engine's ``straggling`` flag, so a
        scattering-only study (``straggling=False``) is possible. The medium's
        radiation length must be known (> 0).

        The 3-D scattering path supports **1-D voxelized heterogeneous materials**
        (decisions 0016, 0017): the energy loss uses the per-voxel water-
        equivalent density and the multiple scattering uses the per-voxel physical
        density and material radiation length, looked up by depth. Every voxel
        material must have a known radiation length (> 0).
        """
        if np.any(self.voxel_radiation_length <= 0.0):
            raise ValueError(
                "a voxel material has no radiation length; multiple scattering is "
                "unavailable (populate radiation_length_g_per_cm2)"
            )
        state = source.sample(n_histories, seed)
        energy_in = float(np.sum(state.energy_mev * state.weight))
        if path == "warp":
            edep, truncated, final_z, final_status = self._run_scattering_warp(
                state, grid, device
            )
        elif path == "python":
            edep, truncated, final_z, final_status = self._run_scattering_reference(
                state, grid
            )
        else:
            raise ValueError(
                f"unknown transport path {path!r} (use 'python' or 'warp')"
            )
        stopped = final_z[final_status == int(Status.STOPPED)]
        n_stopped = int(stopped.shape[0])
        range_mean = float(np.mean(stopped)) if n_stopped else float("nan")
        return ScatteringResult(
            edep_zx_mev=edep,
            grid=grid,
            n_histories=n_histories,
            energy_in_mev=energy_in,
            truncated=truncated,
            path=path,
            device=device if path == "warp" else None,
            range_mean_mm=range_mean,
            n_stopped=n_stopped,
        )

    def _run_scattering_warp(
        self, state: ParticleState, grid: DepthLateralGrid, device: str
    ) -> tuple[np.ndarray, int, np.ndarray, np.ndarray]:
        if not mathlib.HAVE_WARP:
            raise ImportError(
                "warp-lang is not installed; the warp path is unavailable"
            )
        from ionmc.backend.warp_transport import ScatteringKernel

        kernel = ScatteringKernel(self.table, device)
        return kernel.run(
            state=state,
            grid=grid,
            voxel_z_mm=self.voxel_z_mm,
            voxel_density=self.voxel_density,
            voxel_physical_density=self.voxel_physical_density,
            voxel_radiation_length=self.voxel_radiation_length,
            za_ratio=self.depth_dose_za_ratio,
            max_fraction=self.max_fraction,
            max_step_mm=self.max_step_mm,
            geom_depth_mm=self.depth_mm,
            energy_cut_mev=self.energy_cut_mev,
            max_steps=self.max_steps,
            rest_energy_mev=self.particle.rest_energy_mev,
            charge=self.particle.charge,
            straggling=self.straggling,
            straggling_floor_mev=self.straggling_floor_mev,
        )

    def _run_scattering_reference(
        self, state: ParticleState, grid: DepthLateralGrid
    ) -> tuple[np.ndarray, int, np.ndarray, np.ndarray]:
        tp = reference.load_bound_module(
            "ionmc.physics.transport",
            "python",
            rebind_dependencies=["ionmc.physics.tabulated"],
        )
        t = self.table
        args: tuple[Any, ...] = (
            t.energy_mev,
            t.stopping_mev_cm2_per_g,
            t.slope,
            t.size,
            t.bisection_steps,
        )
        # per-voxel physics along the beam axis (decisions 0016, 0017): the
        # energy loss uses the water-equivalent density, the multiple scattering
        # the physical density and material radiation length.
        voxel_z = self.voxel_z_mm
        voxel_density = self.voxel_density  # water-equivalent (stopping)
        voxel_phys = self.voxel_physical_density  # physical (MCS)
        voxel_radlen = self.voxel_radiation_length  # material X0 (MCS)
        n_vox = self.n_voxels
        rest_energy = self.particle.rest_energy_mev
        charge = self.particle.charge
        za = self.depth_dose_za_ratio  # water-equivalent frame (decision 0015)
        straggling = self.straggling
        floor = self.straggling_floor_mev
        geom_depth = self.depth_mm
        dz = grid.depth_bin_mm
        half_w = grid.half_width_mm
        dxb = grid.lateral_bin_mm
        nz, nx = grid.n_depth, grid.n_lateral
        edep = grid.empty()
        truncated = 0
        final_pz = np.zeros(state.size, dtype=np.float64)
        final_status = np.zeros(state.size, dtype=np.int32)
        for h in range(state.size):
            e = float(state.energy_mev[h])
            px, py, pz = (float(v) for v in state.position_mm[h])
            dx, dy, dzr = 0.0, 0.0, 1.0
            w = float(state.weight[h])
            rng = RandomState.from_state(int(state.rng_state[h]))
            voxel = min(int(np.searchsorted(voxel_z, pz, side="right")) - 1, n_vox - 1)
            voxel = max(voxel, 0)
            step = 0
            status = Status.ALIVE
            while status == Status.ALIVE and step < self.max_steps:
                density = float(voxel_density[voxel])
                dzr_pos = max(dzr, 1.0e-6)
                s = tp.energy_loss_step_length(
                    e, self.max_fraction, self.max_step_mm, density, *args
                )
                # limit the step so its depth advance stays within the voxel
                s = min(
                    s,
                    (geom_depth - pz) / dzr_pos,
                    (float(voxel_z[voxel + 1]) - pz) / dzr_pos,
                )
                de = tp.midpoint_energy_loss(e, s, density, *args)
                if straggling and e > floor:
                    sigma = tp.bohr_straggling_sigma(
                        e, rest_energy, s, density, za, charge
                    )
                    de = tp.straggled_energy_loss(de, sigma, rng.randn(), e)
                # random hinge: straight a, scatter, straight (s - a)
                a = rng.randf() * s
                z_start = pz
                x_start = px
                px, py, pz = px + a * dx, py + a * dy, pz + a * dzr
                # scattering is independent of energy-loss straggling: it is
                # always applied above the floor in a scattering run (decision 0011).
                if e > floor:
                    # MCS uses the physical density and material radiation length
                    theta0 = tp.highland_theta0(
                        e,
                        rest_energy,
                        charge,
                        s,
                        float(voxel_phys[voxel]),
                        float(voxel_radlen[voxel]),
                    )
                    dx, dy, dzr = _scatter_direction(
                        dx, dy, dzr, rng.randn() * theta0, rng.randn() * theta0
                    )
                px, py, pz = px + (s - a) * dx, py + (s - a) * dy, pz + (s - a) * dzr
                _deposit_zx(
                    edep,
                    z_start,
                    pz,
                    0.5 * (x_start + px),
                    w * de,
                    dz,
                    half_w,
                    dxb,
                    nz,
                    nx,
                )
                e -= de
                step += 1
                # advance the voxel index across depth boundaries the step reached
                while voxel + 1 < n_vox and pz >= float(voxel_z[voxel + 1]) - 1.0e-9:
                    voxel += 1
                if e <= self.energy_cut_mev:
                    _deposit_zx(edep, pz, pz, px, w * e, dz, half_w, dxb, nz, nx)
                    e = 0.0
                    status = Status.STOPPED
                if pz >= geom_depth:
                    status = Status.ESCAPED
            final_pz[h] = pz
            if step >= self.max_steps and status == Status.ALIVE:
                truncated += 1
                final_status[h] = 3
            else:
                final_status[h] = int(status)
        return edep, truncated, final_pz, final_status

    # -- reference Python path ------------------------------------------------

    def _reaction_local_fraction(self) -> float:
        """Fraction of a reacting primary's energy deposited locally at the
        vertex. With secondary transport on this is the heavy-fragment fraction
        (the secondary protons carry the rest); otherwise the DEV-007 lumped
        local fraction (decision 0013)."""
        if self.secondaries:
            return self.secondary_heavy_fraction
        return self.nuclear_local_fraction

    def _run_reference(
        self, state: ParticleState
    ) -> tuple[
        np.ndarray, int, np.ndarray, np.ndarray, float, int, np.ndarray, np.ndarray
    ]:
        tp = reference.load_bound_module(
            "ionmc.physics.transport",
            "python",
            rebind_dependencies=["ionmc.physics.tabulated"],
        )
        nuc = reference.load_bound_module("ionmc.physics.nuclear", "python")
        t = self.table
        args: tuple[Any, ...] = (
            t.energy_mev,
            t.stopping_mev_cm2_per_g,
            t.slope,
            t.size,
            t.bisection_steps,
        )
        dz = self.grid.bin_width_mm
        geom_depth = self.depth_mm  # transport is bounded by the medium...
        n_bins = self.grid.n_bins  # ...scoring only within the grid extent
        # per-voxel density profile along the beam axis (decision 0014)
        voxel_z = self.voxel_z_mm
        voxel_density = self.voxel_density
        voxel_oxygen = self.voxel_oxygen_density
        n_vox = self.n_voxels
        rest_energy = self.particle.rest_energy_mev
        charge = self.particle.charge
        za = self.depth_dose_za_ratio
        straggling = self.straggling
        floor = self.straggling_floor_mev
        nuclear = self.nuclear
        local_fraction = self._reaction_local_fraction()
        edep = self.grid.empty()
        truncated = 0
        escaped = 0.0
        n_reactions = 0
        final_z = np.zeros(state.size, dtype=np.float64)
        final_status = np.zeros(state.size, dtype=np.int32)
        # per-reacting-history record for secondary generation (decision 0013):
        # vertex depth and residual energy, zero for histories that do not react.
        react_z = np.zeros(state.size, dtype=np.float64)
        react_e = np.zeros(state.size, dtype=np.float64)
        for h in range(state.size):
            e = float(state.energy_mev[h])
            z = float(state.position_mm[h, 2])
            w = float(state.weight[h])
            rng = RandomState.from_state(int(state.rng_state[h]))
            # current voxel index (non-decreasing; forward transport, decision 0014)
            voxel = min(int(np.searchsorted(voxel_z, z, side="right")) - 1, n_vox - 1)
            voxel = max(voxel, 0)
            step = 0
            status = Status.ALIVE
            while status == Status.ALIVE and step < self.max_steps:
                density = float(voxel_density[voxel])
                oxygen_density = float(voxel_oxygen[voxel])
                # physics step: fractional energy loss, not bin-limited, so the
                # straggling and clamp are unbiased (decision 0010). The step is
                # additionally limited to the current voxel boundary so its
                # density is unambiguous (decision 0014). Deposition is split
                # across the bins the step spans (below).
                dl = tp.energy_loss_step_length(
                    e, self.max_fraction, self.max_step_mm, density, *args
                )
                dl = min(dl, geom_depth - z, float(voxel_z[voxel + 1]) - z)
                de = tp.midpoint_energy_loss(e, dl, density, *args)
                if straggling and e > floor:
                    sigma = tp.bohr_straggling_sigma(
                        e, rest_energy, dl, density, za, charge
                    )
                    de = tp.straggled_energy_loss(de, sigma, rng.randn(), e)
                if nuclear:
                    # one uniform per alive step keeps the reference and Warp
                    # streams aligned; the probability is zero below threshold
                    # so the draw never triggers there (decision 0012).
                    p_nuc = nuc.nonelastic_step_probability(e, dl, oxygen_density)
                    if rng.randf() < p_nuc:
                        dep_bin = math.floor(z / dz)
                        if 0 <= dep_bin < n_bins:
                            edep[dep_bin] += w * local_fraction * e
                        escaped += w * (1.0 - local_fraction) * e
                        react_z[h] = z
                        react_e[h] = e
                        n_reactions += 1
                        status = Status.REACTED
                        break
                _deposit_along_step(edep, z, dl, w * de, dz, n_bins)
                z += dl
                e -= de
                step += 1
                # advance the voxel index across any boundaries the step reached
                # (z is non-decreasing, so the index only moves forward)
                while voxel + 1 < n_vox and z >= voxel_z[voxel + 1] - 1.0e-9:
                    voxel += 1
                if e <= self.energy_cut_mev:
                    dep_bin = math.floor(z / dz)
                    if 0 <= dep_bin < n_bins:
                        edep[dep_bin] += w * e
                    e = 0.0
                    status = Status.STOPPED
                if z >= geom_depth:
                    status = Status.ESCAPED
            final_z[h] = z
            if step >= self.max_steps and status == Status.ALIVE:
                truncated += 1
                final_status[h] = 3
            else:
                final_status[h] = int(status)
        return (
            edep,
            truncated,
            final_z,
            final_status,
            escaped,
            n_reactions,
            react_z,
            react_e,
        )

    # -- Warp path ------------------------------------------------------------

    def _run_warp(
        self, state: ParticleState, device: str
    ) -> tuple[
        np.ndarray, int, np.ndarray, np.ndarray, float, int, np.ndarray, np.ndarray
    ]:
        if not mathlib.HAVE_WARP:
            raise ImportError(
                "warp-lang is not installed; the warp path is unavailable"
            )
        from ionmc.backend.warp_transport import DepthDoseKernel

        kernel = DepthDoseKernel(self.table, device)
        return kernel.run(
            energy0=state.energy_mev,
            z0=state.position_mm[:, 2],
            weight=state.weight,
            rng_state=state.rng_state,
            voxel_z_mm=self.voxel_z_mm,
            voxel_density=self.voxel_density,
            voxel_oxygen_density=self.voxel_oxygen_density,
            max_fraction=self.max_fraction,
            max_step_mm=self.max_step_mm,
            bin_width_mm=self.grid.bin_width_mm,
            geom_depth_mm=self.depth_mm,
            energy_cut_mev=self.energy_cut_mev,
            max_steps=self.max_steps,
            n_bins=self.grid.n_bins,
            rest_energy_mev=self.particle.rest_energy_mev,
            charge=self.particle.charge,
            za_ratio=self.depth_dose_za_ratio,
            straggling=self.straggling,
            straggling_floor_mev=self.straggling_floor_mev,
            nuclear=self.nuclear,
            nuclear_local_fraction=self._reaction_local_fraction(),
        )


def _deposit_along_step(
    edep: np.ndarray, z0: float, dl: float, energy: float, dz: float, n_bins: int
) -> None:
    """Distribute ``energy`` uniformly over the depth bins the step ``[z0, z0+dl]``
    spans, proportional to the path length in each bin (decision 0010).

    Mirrors the Warp kernel's deposition loop. Energy that falls outside the
    scoring grid is dropped (the ``0 <= b < n_bins`` guard).
    """
    if dl <= 0.0:
        b = math.floor(z0 / dz)
        if 0 <= b < n_bins:
            edep[b] += energy
        return
    z1 = z0 + dl
    pos = z0
    b = math.floor(z0 / dz)
    inv_dl = 1.0 / dl
    while pos < z1 - 1.0e-12:
        bin_end = (b + 1) * dz
        seg_end = min(bin_end, z1)
        if 0 <= b < n_bins:
            edep[b] += energy * (seg_end - pos) * inv_dl
        pos = seg_end
        b += 1


def _distal_level_depth(edep: np.ndarray, grid: DepthDoseGrid, level: float) -> float:
    """Depth [mm] where the profile falls through ``level * max`` on the distal side.

    Linear interpolation between the two bin centres bracketing the crossing
    beyond the maximum. Returns ``nan`` if the profile is all zero.
    """
    if not np.any(edep > 0.0):
        return float("nan")
    peak = int(np.argmax(edep))
    threshold = level * edep[peak]
    centers = grid.centers_mm
    for i in range(peak, len(edep) - 1):
        if edep[i] >= threshold >= edep[i + 1]:
            e0, e1 = edep[i], edep[i + 1]
            c0, c1 = centers[i], centers[i + 1]
            if e0 == e1:
                return float(c0)
            frac = (e0 - threshold) / (e0 - e1)
            return float(c0 + frac * (c1 - c0))
    return float(centers[-1])
