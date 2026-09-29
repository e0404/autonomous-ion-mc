"""Reference Python transport backend (float64, one particle at a time).

This backend is the correctness oracle of decision 0038: it executes the
shared step physics (``step_physics.py`` bound to Python's ``math``) in an
explicit, readable per-particle loop with NumPy random streams, independent
of the Warp kernels' scheduling, compaction and atomics. It is slow by
design (order 10⁴ steps/s) and intended for validation cases of 10²–10⁴
histories.

Step algorithm (class-II-like condensed history, electromagnetic only):

1. locate the voxel; leave if outside (energy accounted as escaped);
2. candidate step = min(distance to the voxel face, max_step, step that
   loses ``energy_step_fraction`` of the energy per the CSDA range table);
3. mean loss from range-table inversion, then a Gamma-distributed loss with
   Bohr variance (``straggling``);
4. multiple scattering: polar angle from Gottschalk's differential Molière
   scattering power, azimuth uniform, applied at a random fraction of the
   step (random hinge) so lateral displacement is unbiased;
5. deposit the energy along the step into the scoring grid (path-length
   weighted); particles below the cutoff deposit their remaining energy
   locally (counted separately).

Random streams: one ``numpy.random.Generator`` per batch seeded from
``(seed, batch)`` via SeedSequence, so batches are independent and
reproducible.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from types import ModuleType

import numpy as np

from ionmc.config import SimulationConfig
from ionmc.geometry import VoxelGeometry
from ionmc.scoring import ScoringGrid, Tally
from ionmc.transport.shared import physics
from ionmc.transport.tables import TableSet


@dataclass
class EnergyAccounting:
    """Per-primary energy bookkeeping (MeV), summed over all histories."""

    initial: float = 0.0
    deposited_continuous: float = 0.0
    deposited_cutoff: float = 0.0
    escaped: float = 0.0
    truncated_max_steps: float = 0.0
    histories: int = 0
    steps: int = 0

    def per_primary(self) -> dict[str, float]:
        n = max(self.histories, 1)
        return {
            "initial_mev": self.initial / n,
            "deposited_continuous_mev": self.deposited_continuous / n,
            "deposited_cutoff_mev": self.deposited_cutoff / n,
            "escaped_mev": self.escaped / n,
            "truncated_max_steps_mev": self.truncated_max_steps / n,
            "steps_per_history": self.steps / n,
        }


def deposit_along_segment(
    edep: np.ndarray,
    grid: ScoringGrid,
    x0: float,
    y0: float,
    z0: float,
    x1: float,
    y1: float,
    z1: float,
    energy: float,
) -> None:
    """Distribute ``energy`` over the scoring voxels crossed by the segment (path-length weighted)."""
    length = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)
    if energy <= 0.0:
        return
    ox, oy, oz = grid.origin_mm
    dx, dy, dz = grid.spacing_mm
    nx, ny, nz = grid.shape
    if length == 0.0:
        i, j, k = (
            math.floor((x0 - ox) / dx),
            math.floor((y0 - oy) / dy),
            math.floor((z0 - oz) / dz),
        )
        if 0 <= i < nx and 0 <= j < ny and 0 <= k < nz:
            edep[i, j, k] += energy
        return
    ux, uy, uz = (x1 - x0) / length, (y1 - y0) / length, (z1 - z0) / length
    s = 0.0
    x, y, z = x0, y0, z0
    guard = 0
    while s < length - 1e-9 and guard < 10000:
        guard += 1
        i, j, k = (
            math.floor((x - ox) / dx),
            math.floor((y - oy) / dy),
            math.floor((z - oz) / dz),
        )
        # distance to the next scoring face along the direction
        best = length - s
        if ux > 0:
            best = min(best, (ox + (i + 1) * dx - x) / ux)
        elif ux < 0:
            best = min(best, (ox + i * dx - x) / ux)
        if uy > 0:
            best = min(best, (oy + (j + 1) * dy - y) / uy)
        elif uy < 0:
            best = min(best, (oy + j * dy - y) / uy)
        if uz > 0:
            best = min(best, (oz + (k + 1) * dz - z) / uz)
        elif uz < 0:
            best = min(best, (oz + k * dz - z) / uz)
        best = max(best, 0.0)
        seg = best + 1e-9  # nudge across the face
        seg = min(seg, length - s)
        if 0 <= i < nx and 0 <= j < ny and 0 <= k < nz:
            edep[i, j, k] += energy * seg / length
        s += seg
        x, y, z = x0 + ux * s, y0 + uy * s, z0 + uz * s


def transport_history(
    py: ModuleType,
    tables: TableSet,
    geometry: VoxelGeometry,
    grid: ScoringGrid,
    edep: np.ndarray,
    acc: EnergyAccounting,
    rng: np.random.Generator,
    species_index: int,
    x: float,
    y: float,
    z: float,
    ux: float,
    uy: float,
    uz: float,
    t: float,
    cfg: SimulationConfig,
    max_steps: int = 100000,
) -> None:
    ph = cfg.physics
    a = float(tables.mass_number[species_index])
    mass = float(tables.mass_mev[species_index])
    z_charge = float(tables.charge[species_index])
    mass_per_nucleon = mass / a
    ox, oy, oz = geometry.origin_mm
    dx, dy, dz = geometry.spacing_mm
    nx, ny, nz = geometry.shape
    nt = tables.t_grid.size
    t_log_min, t_inv = tables.t_log_min, tables.t_inv_dlog
    p1v1 = py.pv_mev(t, mass, a)
    acc.initial += t * a
    steps = 0
    # vacuum flight to the geometry box (no interactions outside the phantom)
    d_entry = py.distance_to_box_entry(
        x, y, z, ux, uy, uz, ox, oy, oz, ox + nx * dx, oy + ny * dy, oz + nz * dz
    )
    if d_entry < 0.0:
        acc.escaped += t * a
        acc.histories += 1
        return
    x, y, z = x + ux * d_entry, y + uy * d_entry, z + uz * d_entry
    while True:
        if t <= ph.cutoff_mev_per_u:
            deposit_along_segment(edep, grid, x, y, z, x, y, z, t * a)
            acc.deposited_cutoff += t * a
            break
        if steps >= max_steps:
            acc.truncated_max_steps += t * a
            break
        i = py.voxel_axis_index(x, ox, dx, nx)
        j = py.voxel_axis_index(y, oy, dy, ny)
        k = py.voxel_axis_index(z, oz, dz, nz)
        if i < 0 or j < 0 or k < 0:
            acc.escaped += t * a
            break
        m = int(geometry.material_index[i, j, k])
        rho = float(geometry.density_g_cm3[i, j, k])
        # candidate step (mm)
        d_face = py.distance_to_voxel_boundary(
            x, y, z, ux, uy, uz, ox, oy, oz, dx, dy, dz
        )
        r_now = py.range_lookup(
            tables.log_range, species_index, m, t, t_log_min, t_inv, nt
        )
        r_frac = py.range_lookup(
            tables.log_range,
            species_index,
            m,
            t * (1.0 - ph.energy_step_fraction),
            t_log_min,
            t_inv,
            nt,
        )
        d_energy = (r_now - r_frac) / rho * 10.0  # g/cm² -> mm at density rho
        step = min(d_face + 1e-6, ph.max_step_mm, max(d_energy, 1e-3))
        rho_path = rho * step * 0.1  # g/cm²
        # mean energy loss via range inversion
        t_new = py.energy_after_path(
            tables.log_range, species_index, m, t, rho_path, t_log_min, t_inv, nt
        )
        de = max((t - t_new) * a, 0.0)
        if ph.straggling and t_new > 0.0:
            beta2 = py.beta_squared(t, mass_per_nucleon)
            z_eff = py.effective_charge(z_charge, math.sqrt(beta2))
            var = py.straggling_variance_mev2(
                z_eff, float(tables.z_over_a[m]), beta2, rho_path
            )
            de = py.sample_energy_loss(rng, de, var)
            de = min(de, t * a)
        # multiple scattering with random hinge
        hinge = rng.random() if ph.multiple_scattering else 1.0
        x0, y0, z0 = x, y, z
        xh, yh, zh = x + ux * step * hinge, y + uy * step * hinge, z + uz * step * hinge
        if ph.multiple_scattering:
            beta2 = py.beta_squared(t, mass_per_nucleon)
            pv = py.pv_mev(t, mass, a)
            theta0_sq = py.mcs_theta0_squared(
                z_charge, pv, p1v1, float(tables.rho_x_s[m]), rho_path
            )
            theta = math.sqrt(theta0_sq) * math.sqrt(
                -2.0 * math.log(max(rng.random(), 1e-300))
            )
            theta = min(theta, 1.0)
            phi = 2.0 * math.pi * rng.random()
            ct, st = math.cos(theta), math.sin(theta)
            nux = py.rotate_x(ux, uy, uz, ct, st, phi)
            nuy = py.rotate_y(ux, uy, uz, ct, st, phi)
            nuz = py.rotate_z(ux, uy, uz, ct, st, phi)
            norm = math.sqrt(nux * nux + nuy * nuy + nuz * nuz)
            ux, uy, uz = nux / norm, nuy / norm, nuz / norm
        x1, y1, z1 = (
            xh + ux * step * (1.0 - hinge),
            yh + uy * step * (1.0 - hinge),
            zh + uz * step * (1.0 - hinge),
        )
        # deposit along the two legs of the hinged step
        deposit_along_segment(edep, grid, x0, y0, z0, xh, yh, zh, de * hinge)
        deposit_along_segment(edep, grid, xh, yh, zh, x1, y1, z1, de * (1.0 - hinge))
        acc.deposited_continuous += de
        t = t - de / a
        x, y, z = x1, y1, z1
        steps += 1
    acc.steps += steps
    acc.histories += 1


def run_reference(config: SimulationConfig, tables: TableSet) -> dict:
    """Run all batches on the reference backend; returns tallies and accounting."""
    py = physics("python")
    cfg = config.with_defaults()
    grid = cfg.scoring
    assert grid is not None
    sp_index = tables.species_index(cfg.source.species)
    energy = Tally("energy", grid.shape, cfg.batches)
    acc = EnergyAccounting()
    seeds = np.random.SeedSequence(cfg.seed).spawn(cfg.batches)
    per_batch = [
        cfg.histories // cfg.batches + (1 if b < cfg.histories % cfg.batches else 0)
        for b in range(cfg.batches)
    ]
    t0 = time.perf_counter()
    for b in range(cfg.batches):
        rng = np.random.default_rng(seeds[b])
        primaries = cfg.source.sample(per_batch[b], rng)
        for n in range(per_batch[b]):
            px, py_, pz = primaries["position_mm"][n]
            ux, uy, uz = primaries["direction"][n]
            transport_history(
                py,
                tables,
                cfg.geometry,
                grid,
                energy.current,
                acc,
                rng,
                sp_index,
                float(px),
                float(py_),
                float(pz),
                float(ux),
                float(uy),
                float(uz),
                float(primaries["energy_mev_per_u"][n]),
                cfg,
            )
        energy.close_batch(per_batch[b])
    wall = time.perf_counter() - t0
    return {
        "energy": energy,
        "accounting": acc,
        "histories_per_batch": per_batch,
        "wall_seconds": wall,
        "rng": "numpy PCG64 via SeedSequence(seed).spawn(batches)",
    }
