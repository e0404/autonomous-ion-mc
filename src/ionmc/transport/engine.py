"""Transport engine: reference Python and Warp execution of CSDA depth dose.

The engine runs a monoenergetic proton pencil beam through a homogeneous water
slab and returns the integral depth dose (energy deposited per depth bin). The
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
from ionmc.particles import PROTON, Particle
from ionmc.rng import RandomState
from ionmc.transport.depth_dose import DepthDoseGrid, DepthLateralGrid
from ionmc.transport.geometry import WaterSlab
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

    @property
    def energy_deposited_mev(self) -> float:
        return float(np.sum(self.edep_mev))

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
        slab: WaterSlab,
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
    ) -> None:
        if not (0.0 < max_fraction < 1.0):
            raise ValueError("max_fraction must be in (0, 1)")
        if not (0.0 <= nuclear_local_fraction <= 1.0):
            raise ValueError("nuclear_local_fraction must be in [0, 1]")
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
        #: Electrons per gram <Z/A> of the medium (Bohr straggling; decision 0010).
        self.za_ratio = slab.material.electrons_per_gram_ratio
        #: Radiation length [g/cm^2] of the medium (MCS; decision 0011).
        self.radiation_length_g_per_cm2 = slab.material.radiation_length_g_per_cm2
        #: Oxygen number density [1/cm^3] of the medium: only oxygen contributes
        #: catastrophic nonelastic removal of primaries (decision 0012). Zero if
        #: the medium has no oxygen, which disables nuclear removal there.
        material = slab.material
        if "O" in material.mass_fractions:
            self.oxygen_density_per_cm3 = (
                AVOGADRO * material.atoms_per_gram("O") * material.density_g_per_cm3
            )
        else:
            self.oxygen_density_per_cm3 = 0.0

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
        if path == "warp":
            edep, truncated, final_z, final_status, escaped, n_reactions = (
                self._run_warp(state, device)
            )
        elif path == "python":
            edep, truncated, final_z, final_status, escaped, n_reactions = (
                self._run_reference(state)
            )
        else:
            raise ValueError(
                f"unknown transport path {path!r} (use 'python' or 'warp')"
            )
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
        """
        if self.radiation_length_g_per_cm2 <= 0.0:
            raise ValueError(
                "the medium has no radiation length; multiple scattering is unavailable"
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
            density=self.slab.density_g_per_cm3,
            radiation_length_g_per_cm2=self.radiation_length_g_per_cm2,
            max_fraction=self.max_fraction,
            max_step_mm=self.max_step_mm,
            geom_depth_mm=self.slab.depth_mm,
            energy_cut_mev=self.energy_cut_mev,
            max_steps=self.max_steps,
            rest_energy_mev=self.particle.rest_energy_mev,
            charge=self.particle.charge,
            za_ratio=self.za_ratio,
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
        density = self.slab.density_g_per_cm3
        radlen = self.radiation_length_g_per_cm2
        rest_energy = self.particle.rest_energy_mev
        charge = self.particle.charge
        za = self.za_ratio
        straggling = self.straggling
        floor = self.straggling_floor_mev
        geom_depth = self.slab.depth_mm
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
            step = 0
            status = Status.ALIVE
            while status == Status.ALIVE and step < self.max_steps:
                s = tp.energy_loss_step_length(
                    e, self.max_fraction, self.max_step_mm, density, *args
                )
                s = min(s, (geom_depth - pz) / max(dzr, 1.0e-6))
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
                    theta0 = tp.highland_theta0(
                        e, rest_energy, charge, s, density, radlen
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

    def _run_reference(
        self, state: ParticleState
    ) -> tuple[np.ndarray, int, np.ndarray, np.ndarray, float, int]:
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
        density = self.slab.density_g_per_cm3
        dz = self.grid.bin_width_mm
        geom_depth = self.slab.depth_mm  # transport is bounded by the medium...
        n_bins = self.grid.n_bins  # ...scoring only within the grid extent
        rest_energy = self.particle.rest_energy_mev
        charge = self.particle.charge
        za = self.za_ratio
        straggling = self.straggling
        floor = self.straggling_floor_mev
        nuclear = self.nuclear
        oxygen_density = self.oxygen_density_per_cm3
        local_fraction = self.nuclear_local_fraction
        edep = self.grid.empty()
        truncated = 0
        escaped = 0.0
        n_reactions = 0
        final_z = np.zeros(state.size, dtype=np.float64)
        final_status = np.zeros(state.size, dtype=np.int32)
        for h in range(state.size):
            e = float(state.energy_mev[h])
            z = float(state.position_mm[h, 2])
            w = float(state.weight[h])
            rng = RandomState.from_state(int(state.rng_state[h]))
            step = 0
            status = Status.ALIVE
            while status == Status.ALIVE and step < self.max_steps:
                # physics step: fractional energy loss, not bin-limited, so the
                # straggling and clamp are unbiased (decision 0010). Deposition
                # is split across the bins the step spans (below).
                dl = tp.energy_loss_step_length(
                    e, self.max_fraction, self.max_step_mm, density, *args
                )
                dl = min(dl, geom_depth - z)
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
                        n_reactions += 1
                        status = Status.REACTED
                        break
                _deposit_along_step(edep, z, dl, w * de, dz, n_bins)
                z += dl
                e -= de
                step += 1
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
        return edep, truncated, final_z, final_status, escaped, n_reactions

    # -- Warp path ------------------------------------------------------------

    def _run_warp(
        self, state: ParticleState, device: str
    ) -> tuple[np.ndarray, int, np.ndarray, np.ndarray, float, int]:
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
            density=self.slab.density_g_per_cm3,
            max_fraction=self.max_fraction,
            max_step_mm=self.max_step_mm,
            bin_width_mm=self.grid.bin_width_mm,
            geom_depth_mm=self.slab.depth_mm,
            energy_cut_mev=self.energy_cut_mev,
            max_steps=self.max_steps,
            n_bins=self.grid.n_bins,
            rest_energy_mev=self.particle.rest_energy_mev,
            charge=self.particle.charge,
            za_ratio=self.za_ratio,
            straggling=self.straggling,
            straggling_floor_mev=self.straggling_floor_mev,
            nuclear=self.nuclear,
            oxygen_density_per_cm3=self.oxygen_density_per_cm3,
            nuclear_local_fraction=self.nuclear_local_fraction,
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
