"""Warp kernel for continuous-slowing-down proton depth-dose transport.

One history per thread with an in-kernel step loop and atomic scoring into a
1-D depth grid (decision ``0009``). The kernel contains no physics: the step
length and the midpoint energy loss come from the shared-source functions in
:mod:`ionmc.physics.transport`, which also run as the float64 reference driver.
Runs in float32 on the CPU and CUDA devices.
"""

# mypy: disable-error-code="valid-type"

from __future__ import annotations

from typing import Any

import numpy as np
import warp as wp

from ionmc.data.stopping_tables import StoppingTable
from ionmc.physics import nuclear as nuclear_phys
from ionmc.physics import transport


@wp.kernel
def csda_depth_dose_kernel(
    energy0: wp.array(dtype=float),
    z0: wp.array(dtype=float),
    weight: wp.array(dtype=float),
    rng_state0: wp.array(dtype=wp.uint32),
    voxel_z: wp.array(dtype=float),
    voxel_density: wp.array(dtype=float),
    voxel_oxygen: wp.array(dtype=float),
    n_vox: int,
    max_fraction: float,
    max_step_mm: float,
    bin_width_mm: float,
    geom_depth_mm: float,
    energy_cut_mev: float,
    max_steps: int,
    rest_energy_mev: float,
    charge: float,
    za_ratio: float,
    straggling: int,
    straggling_floor_mev: float,
    nuclear: int,
    nuclear_local_fraction: float,
    table_e: wp.array(dtype=float),
    table_s: wp.array(dtype=float),
    table_d: wp.array(dtype=float),
    n: int,
    n_steps: int,
    n_bins: int,
    edep: wp.array(dtype=wp.float64),
    truncated: wp.array(dtype=int),
    final_z: wp.array(dtype=float),
    final_status: wp.array(dtype=int),
    escaped: wp.array(dtype=wp.float64),
    reactions: wp.array(dtype=int),
    react_z: wp.array(dtype=float),
    react_e: wp.array(dtype=float),
):
    i = wp.tid()
    e = energy0[i]
    z = z0[i]
    w = weight[i]
    rng = rng_state0[i]
    # dynamic (mutable) loop variables: int(...) keeps Warp from treating them
    # as constants inside the while loop (noqa keeps ruff from stripping int()).
    step = int(0)  # noqa: UP018, RUF046
    alive = int(1)  # noqa: UP018, RUF046
    reacted = int(0)  # noqa: UP018, RUF046
    # current voxel index for the per-voxel density profile (decision 0014):
    # advance to the voxel containing z0 (non-decreasing thereafter)
    voxel = int(0)  # noqa: UP018, RUF046
    while voxel + 1 < n_vox and z >= voxel_z[voxel + 1]:
        voxel = voxel + 1
    while alive == 1 and step < max_steps:
        density = voxel_density[voxel]
        ox = voxel_oxygen[voxel]
        # physics step: fractional energy loss, not bin-limited (decision 0010),
        # so the straggling and clamp are unbiased; additionally limited to the
        # current voxel boundary so the density is unambiguous (decision 0014);
        # deposition is split across the bins the step spans.
        dl_e = transport.energy_loss_step_length(
            e, max_fraction, max_step_mm, density, table_e, table_s, table_d, n, n_steps
        )
        dl = wp.min(wp.min(dl_e, geom_depth_mm - z), voxel_z[voxel + 1] - z)
        de = transport.midpoint_energy_loss(
            e, dl, density, table_e, table_s, table_d, n, n_steps
        )
        if straggling == 1 and e > straggling_floor_mev:
            sigma = transport.bohr_straggling_sigma(
                e, rest_energy_mev, dl, density, za_ratio, charge
            )
            variate = wp.randn(rng)
            de = transport.straggled_energy_loss(de, sigma, variate, e)
        if nuclear == 1:
            # one uniform per alive step keeps the Warp and reference streams
            # aligned; the probability is zero below threshold (decision 0012).
            p_nuc = nuclear_phys.nonelastic_step_probability(e, dl, ox)
            if wp.randf(rng) < p_nuc:
                rbin = int(wp.floor(z / bin_width_mm))
                if rbin >= 0 and rbin < n_bins:
                    wp.atomic_add(
                        edep, rbin, wp.float64(w * nuclear_local_fraction * e)
                    )
                wp.atomic_add(
                    escaped, 0, wp.float64(w * (1.0 - nuclear_local_fraction) * e)
                )
                wp.atomic_add(reactions, 0, 1)
                react_z[i] = z
                react_e[i] = e
                reacted = 1
                alive = 0
                break
        deposit = w * de
        z1 = z + dl
        pos = z
        b = int(wp.floor(z / bin_width_mm))
        inv_dl = 1.0 / wp.max(dl, 1.0e-12)
        while pos < z1 - 1.0e-12:
            bin_end = float(b + 1) * bin_width_mm
            seg_end = wp.min(bin_end, z1)
            if b >= 0 and b < n_bins:
                wp.atomic_add(edep, b, wp.float64(deposit * (seg_end - pos) * inv_dl))
            pos = seg_end
            b = b + 1
        z = z + dl
        e = e - de
        step = step + 1
        # advance across any voxel boundaries the step reached (z non-decreasing)
        while voxel + 1 < n_vox and z >= voxel_z[voxel + 1] - 1.0e-9:
            voxel = voxel + 1
        if e <= energy_cut_mev:
            dep_bin = int(wp.floor(z / bin_width_mm))
            if dep_bin >= 0 and dep_bin < n_bins:
                wp.atomic_add(edep, dep_bin, wp.float64(w * e))
            e = 0.0
            alive = 0
        if z >= geom_depth_mm:
            alive = 0
    # per-history outcome: status 1 stopped, 2 escaped, 3 truncated, 4 reacted
    final_z[i] = z
    if reacted == 1:
        final_status[i] = 4
    elif step >= max_steps and alive == 1:
        wp.atomic_add(truncated, 0, 1)
        final_status[i] = 3
    elif z >= geom_depth_mm:
        final_status[i] = 2
    else:
        final_status[i] = 1


class DepthDoseKernel:
    """Launches the CSDA depth-dose kernel for one table on one device."""

    def __init__(self, table: StoppingTable, device: str = "cpu") -> None:
        wp.init()
        self.device = device
        self.table = table
        self.n = int(table.size)
        self.n_steps = int(table.bisection_steps)
        self.tables: dict[str, Any] = {
            name: wp.array(
                np.asarray(getattr(table, attr), dtype=np.float32),
                dtype=float,
                device=device,
            )
            for name, attr in (
                ("e", "energy_mev"),
                ("s", "stopping_mev_cm2_per_g"),
                ("d", "slope"),
            )
        }

    def run(
        self,
        energy0: np.ndarray,
        z0: np.ndarray,
        weight: np.ndarray,
        rng_state: np.ndarray,
        voxel_z_mm: np.ndarray,
        voxel_density: np.ndarray,
        voxel_oxygen_density: np.ndarray,
        max_fraction: float,
        max_step_mm: float,
        bin_width_mm: float,
        geom_depth_mm: float,
        energy_cut_mev: float,
        max_steps: int,
        n_bins: int,
        rest_energy_mev: float,
        charge: float,
        za_ratio: float,
        straggling: bool,
        straggling_floor_mev: float,
        nuclear: bool = False,
        nuclear_local_fraction: float = 0.0,
    ) -> tuple[
        np.ndarray, int, np.ndarray, np.ndarray, float, int, np.ndarray, np.ndarray
    ]:
        """Return (edep per bin [MeV], truncated, final z, final status, escaped
        energy [MeV], number of nonelastic reactions, per-history reaction vertex
        depth [mm], per-history reaction residual energy [MeV])."""
        d = self.device
        n_hist = int(energy0.shape[0])
        e0: Any = wp.array(
            np.ascontiguousarray(energy0, dtype=np.float32), dtype=float, device=d
        )
        zz: Any = wp.array(
            np.ascontiguousarray(z0, dtype=np.float32), dtype=float, device=d
        )
        ww: Any = wp.array(
            np.ascontiguousarray(weight, dtype=np.float32), dtype=float, device=d
        )
        rng: Any = wp.array(
            np.ascontiguousarray(rng_state, dtype=np.uint32),
            dtype=wp.uint32,
            device=d,
        )
        vz: Any = wp.array(
            np.ascontiguousarray(voxel_z_mm, dtype=np.float32), dtype=float, device=d
        )
        vrho: Any = wp.array(
            np.ascontiguousarray(voxel_density, dtype=np.float32), dtype=float, device=d
        )
        vox: Any = wp.array(
            np.ascontiguousarray(voxel_oxygen_density, dtype=np.float32),
            dtype=float,
            device=d,
        )
        n_vox = int(voxel_density.shape[0])
        edep = wp.zeros(n_bins, dtype=wp.float64, device=d)
        truncated = wp.zeros(1, dtype=int, device=d)
        final_z = wp.zeros(n_hist, dtype=float, device=d)
        final_status = wp.zeros(n_hist, dtype=int, device=d)
        escaped = wp.zeros(1, dtype=wp.float64, device=d)
        reactions = wp.zeros(1, dtype=int, device=d)
        react_z = wp.zeros(n_hist, dtype=float, device=d)
        react_e = wp.zeros(n_hist, dtype=float, device=d)
        t = self.tables
        wp.launch(
            csda_depth_dose_kernel,
            dim=n_hist,
            inputs=[
                e0,
                zz,
                ww,
                rng,
                vz,
                vrho,
                vox,
                n_vox,
                float(max_fraction),
                float(max_step_mm),
                float(bin_width_mm),
                float(geom_depth_mm),
                float(energy_cut_mev),
                int(max_steps),
                float(rest_energy_mev),
                float(charge),
                float(za_ratio),
                (1 if straggling else 0),
                float(straggling_floor_mev),
                (1 if nuclear else 0),
                float(nuclear_local_fraction),
                t["e"],
                t["s"],
                t["d"],
                self.n,
                self.n_steps,
                int(n_bins),
                edep,
                truncated,
                final_z,
                final_status,
                escaped,
                reactions,
                react_z,
                react_e,
            ],
            device=d,
        )
        wp.synchronize_device(d)
        return (
            edep.numpy().astype(np.float64),  # already float64 on device
            int(truncated.numpy()[0]),
            final_z.numpy().astype(np.float64),
            final_status.numpy().astype(np.int32),
            float(escaped.numpy()[0]),
            int(reactions.numpy()[0]),
            react_z.numpy().astype(np.float64),
            react_e.numpy().astype(np.float64),
        )


@wp.func
def _scatter_dir(dx: float, dy: float, dz: float, tx: float, ty: float) -> wp.vec3:
    """Tilt a unit direction by projected angles ``tx, ty`` in its transverse
    frame (mirrors ``engine._scatter_direction``; two-plane sampler)."""
    ax = wp.abs(dx)
    ay = wp.abs(dy)
    az = wp.abs(dz)
    rx = float(0.0)  # noqa: UP018
    ry = float(0.0)  # noqa: UP018
    rz = float(0.0)  # noqa: UP018
    if ax <= ay and ax <= az:
        rx = 1.0
    elif ay <= az:
        ry = 1.0
    else:
        rz = 1.0
    c1x = dy * rz - dz * ry
    c1y = dz * rx - dx * rz
    c1z = dx * ry - dy * rx
    inv1 = 1.0 / wp.sqrt(c1x * c1x + c1y * c1y + c1z * c1z)
    e1x = c1x * inv1
    e1y = c1y * inv1
    e1z = c1z * inv1
    e2x = dy * e1z - dz * e1y
    e2y = dz * e1x - dx * e1z
    e2z = dx * e1y - dy * e1x
    nx = dx + tx * e1x + ty * e2x
    ny = dy + tx * e1y + ty * e2y
    nz = dz + tx * e1z + ty * e2z
    inv = 1.0 / wp.sqrt(nx * nx + ny * ny + nz * nz)
    return wp.vec3(nx * inv, ny * inv, nz * inv)


@wp.kernel
def csda_scattering_kernel(
    energy0: wp.array(dtype=float),
    weight: wp.array(dtype=float),
    rng_state0: wp.array(dtype=wp.uint32),
    density: float,
    radlen: float,
    max_fraction: float,
    max_step_mm: float,
    geom_depth_mm: float,
    energy_cut_mev: float,
    max_steps: int,
    rest_energy_mev: float,
    charge: float,
    za_ratio: float,
    straggling: int,
    straggling_floor_mev: float,
    depth_bin_mm: float,
    half_width_mm: float,
    lateral_bin_mm: float,
    n_depth: int,
    n_lateral: int,
    table_e: wp.array(dtype=float),
    table_s: wp.array(dtype=float),
    table_d: wp.array(dtype=float),
    n: int,
    n_steps: int,
    edep: wp.array(dtype=wp.float64),
    truncated: wp.array(dtype=int),
    final_z: wp.array(dtype=float),
    final_status: wp.array(dtype=int),
):
    i = wp.tid()
    e = energy0[i]
    w = weight[i]
    rng = rng_state0[i]
    px = float(0.0)  # noqa: UP018
    py = float(0.0)  # noqa: UP018
    pz = float(0.0)  # noqa: UP018
    dx = float(0.0)  # noqa: UP018
    dy = float(0.0)  # noqa: UP018
    dz = float(1.0)  # noqa: UP018
    step = int(0)  # noqa: UP018, RUF046
    alive = int(1)  # noqa: UP018, RUF046
    while alive == 1 and step < max_steps:
        s = transport.energy_loss_step_length(
            e, max_fraction, max_step_mm, density, table_e, table_s, table_d, n, n_steps
        )
        s = wp.min(s, (geom_depth_mm - pz) / wp.max(dz, 1.0e-6))
        de = transport.midpoint_energy_loss(
            e, s, density, table_e, table_s, table_d, n, n_steps
        )
        if straggling == 1 and e > straggling_floor_mev:
            sigma = transport.bohr_straggling_sigma(
                e, rest_energy_mev, s, density, za_ratio, charge
            )
            de = transport.straggled_energy_loss(de, sigma, wp.randn(rng), e)
        a = wp.randf(rng) * s
        z_start = pz
        x_start = px
        px = px + a * dx
        py = py + a * dy
        pz = pz + a * dz
        if e > straggling_floor_mev:
            theta0 = transport.highland_theta0(
                e, rest_energy_mev, charge, s, density, radlen
            )
            nd = _scatter_dir(
                dx, dy, dz, wp.randn(rng) * theta0, wp.randn(rng) * theta0
            )
            dx = nd[0]
            dy = nd[1]
            dz = nd[2]
        px = px + (s - a) * dx
        py = py + (s - a) * dy
        pz = pz + (s - a) * dz
        # deposit w*de across depth bins [z_start, pz] at lateral bin of x mid
        x_mid = 0.5 * (x_start + px)
        xb = int(wp.floor((x_mid + half_width_mm) / lateral_bin_mm))
        if xb >= 0 and xb < n_lateral:
            span = pz - z_start
            if span <= 0.0:
                b0 = int(wp.floor(z_start / depth_bin_mm))
                if b0 >= 0 and b0 < n_depth:
                    wp.atomic_add(edep, b0 * n_lateral + xb, wp.float64(w * de))
            else:
                invs = 1.0 / span
                pos = z_start
                b = int(wp.floor(z_start / depth_bin_mm))
                while pos < pz - 1.0e-12:
                    bin_end = float(b + 1) * depth_bin_mm
                    seg_end = wp.min(bin_end, pz)
                    if b >= 0 and b < n_depth:
                        wp.atomic_add(
                            edep,
                            b * n_lateral + xb,
                            wp.float64(w * de * (seg_end - pos) * invs),
                        )
                    pos = seg_end
                    b = b + 1
        e = e - de
        step = step + 1
        if e <= energy_cut_mev:
            xb2 = int(wp.floor((px + half_width_mm) / lateral_bin_mm))
            b2 = int(wp.floor(pz / depth_bin_mm))
            if xb2 >= 0 and xb2 < n_lateral and b2 >= 0 and b2 < n_depth:
                wp.atomic_add(edep, b2 * n_lateral + xb2, wp.float64(w * e))
            e = 0.0
            alive = 0
        if pz >= geom_depth_mm:
            alive = 0
    final_z[i] = pz
    if step >= max_steps and alive == 1:
        wp.atomic_add(truncated, 0, 1)
        final_status[i] = 3
    elif pz >= geom_depth_mm:
        final_status[i] = 2
    else:
        final_status[i] = 1


class ScatteringKernel:
    """Launches the 3-D scattering depth-lateral kernel for one table/device."""

    def __init__(self, table: StoppingTable, device: str = "cpu") -> None:
        wp.init()
        self.device = device
        self.table = table
        self.n = int(table.size)
        self.n_steps = int(table.bisection_steps)
        self.tables: dict[str, Any] = {
            name: wp.array(
                np.asarray(getattr(table, attr), dtype=np.float32),
                dtype=float,
                device=device,
            )
            for name, attr in (
                ("e", "energy_mev"),
                ("s", "stopping_mev_cm2_per_g"),
                ("d", "slope"),
            )
        }

    def run(
        self,
        state: Any,
        grid: Any,
        density: float,
        radiation_length_g_per_cm2: float,
        max_fraction: float,
        max_step_mm: float,
        geom_depth_mm: float,
        energy_cut_mev: float,
        max_steps: int,
        rest_energy_mev: float,
        charge: float,
        za_ratio: float,
        straggling: bool,
        straggling_floor_mev: float,
    ) -> tuple[np.ndarray, int, np.ndarray, np.ndarray]:
        d = self.device
        n_hist = int(state.energy_mev.shape[0])
        nz, nx = int(grid.n_depth), int(grid.n_lateral)
        e0: Any = wp.array(
            np.ascontiguousarray(state.energy_mev, dtype=np.float32),
            dtype=float,
            device=d,
        )
        ww: Any = wp.array(
            np.ascontiguousarray(state.weight, dtype=np.float32), dtype=float, device=d
        )
        rng: Any = wp.array(
            np.ascontiguousarray(state.rng_state, dtype=np.uint32),
            dtype=wp.uint32,
            device=d,
        )
        edep = wp.zeros(nz * nx, dtype=wp.float64, device=d)
        truncated = wp.zeros(1, dtype=int, device=d)
        final_z = wp.zeros(n_hist, dtype=float, device=d)
        final_status = wp.zeros(n_hist, dtype=int, device=d)
        t = self.tables
        wp.launch(
            csda_scattering_kernel,
            dim=n_hist,
            inputs=[
                e0,
                ww,
                rng,
                float(density),
                float(radiation_length_g_per_cm2),
                float(max_fraction),
                float(max_step_mm),
                float(geom_depth_mm),
                float(energy_cut_mev),
                int(max_steps),
                float(rest_energy_mev),
                float(charge),
                float(za_ratio),
                (1 if straggling else 0),
                float(straggling_floor_mev),
                float(grid.depth_bin_mm),
                float(grid.half_width_mm),
                float(grid.lateral_bin_mm),
                nz,
                nx,
                t["e"],
                t["s"],
                t["d"],
                self.n,
                self.n_steps,
                edep,
                truncated,
                final_z,
                final_status,
            ],
            device=d,
        )
        wp.synchronize_device(d)
        return (
            edep.numpy().astype(np.float64).reshape(nz, nx),
            int(truncated.numpy()[0]),
            final_z.numpy().astype(np.float64),
            final_status.numpy().astype(np.int32),
        )
