"""NVIDIA Warp backends (CPU and CUDA) executing the shared step physics.

One thread transports one primary from its birth to its death (a bounded
``while`` loop inside the kernel; no per-step host involvement, V1-MUST-035).
The kernel is *generic* in the floating-point type: the same source is
instantiated for float32 (default) and float64 (numerical falsification).
Scoring uses atomic adds into a per-batch device grid that is folded into a
float64 host tally after each batch; energy accounting uses a small atomic
array. Random streams: ``wp.rand_init(seed, batch_offset + tid)`` gives every
history its own counter-based stream, independent of thread scheduling
(V1-MUST-028). Primaries are sampled on the host with the same
``PencilBeam.sample`` definitions and NumPy streams as the reference backend.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np
import warp as wp

from ionmc.config import SimulationConfig, UnsupportedConfigurationError
from ionmc.scoring import Tally
from ionmc.transport.reference import EnergyAccounting
from ionmc.transport.shared import physics
from ionmc.transport.tables import TableSet

wpv = physics("warp")

MAX_STEPS = 200000


@wp.kernel
def transport_kernel(
    pos: wp.array2d(dtype=Any),
    dirs: wp.array2d(dtype=Any),
    energy: wp.array(dtype=Any),
    seed: int,
    batch_offset: int,
    species: int,
    charge: Any,
    mass_number: Any,
    mass_mev: Any,
    log_range: wp.array3d(dtype=Any),
    t_log_min: Any,
    t_inv_dlog: Any,
    nt: int,
    z_over_a: wp.array(dtype=Any),
    rho_x_s: wp.array(dtype=Any),
    mat_index: wp.array3d(dtype=wp.int32),
    density: wp.array3d(dtype=Any),
    ox: Any,
    oy: Any,
    oz: Any,
    dx: Any,
    dy: Any,
    dz: Any,
    nx: int,
    ny: int,
    nz: int,
    edep: wp.array3d(dtype=Any),
    sox: Any,
    soy: Any,
    soz: Any,
    sdx: Any,
    sdy: Any,
    sdz: Any,
    snx: int,
    sny: int,
    snz: int,
    straggling: int,
    mcs: int,
    max_step_mm: Any,
    energy_step_fraction: Any,
    cutoff: Any,
    acc: wp.array(dtype=wp.float64),
    steps_out: wp.array(dtype=wp.int32),
):
    tid = wp.tid()
    state = wp.rand_init(seed, batch_offset + tid)
    t = energy[tid]
    x = pos[tid, 0]
    y = pos[tid, 1]
    z = pos[tid, 2]
    ux = dirs[tid, 0]
    uy = dirs[tid, 1]
    uz = dirs[tid, 2]
    zero = type(t)(0.0)
    one = type(t)(1.0)
    mass_per_nucleon = mass_mev / mass_number
    p1v1 = wpv.pv_mev(t, mass_mev, mass_number)
    # vacuum flight to the geometry box
    d_entry = wpv.distance_to_box_entry(
        x,
        y,
        z,
        ux,
        uy,
        uz,
        ox,
        oy,
        oz,
        ox + type(t)(nx) * dx,
        oy + type(t)(ny) * dy,
        oz + type(t)(nz) * dz,
    )
    if d_entry < zero:
        wp.atomic_add(acc, 2, wp.float64(t * mass_number))
        steps_out[tid] = 0
        return
    x = x + ux * d_entry
    y = y + uy * d_entry
    z = z + uz * d_entry
    steps = int(0)
    alive = int(1)
    while alive == 1:
        if t <= cutoff:
            wpv.deposit_segment(
                edep,
                sox,
                soy,
                soz,
                sdx,
                sdy,
                sdz,
                snx,
                sny,
                snz,
                x,
                y,
                z,
                x,
                y,
                z,
                t * mass_number,
            )
            wp.atomic_add(acc, 1, wp.float64(t * mass_number))
            alive = 0
            break
        if steps >= MAX_STEPS:
            wp.atomic_add(acc, 3, wp.float64(t * mass_number))
            alive = 0
            break
        i = wpv.voxel_axis_index(x, ox, dx, nx)
        j = wpv.voxel_axis_index(y, oy, dy, ny)
        k = wpv.voxel_axis_index(z, oz, dz, nz)
        if i < 0 or j < 0 or k < 0:
            wp.atomic_add(acc, 2, wp.float64(t * mass_number))
            alive = 0
            break
        m = int(mat_index[i, j, k])
        rho = density[i, j, k]
        d_face = wpv.distance_to_voxel_boundary(
            x, y, z, ux, uy, uz, ox, oy, oz, dx, dy, dz
        )
        r_now = wpv.range_lookup(log_range, species, m, t, t_log_min, t_inv_dlog, nt)
        r_frac = wpv.range_lookup(
            log_range,
            species,
            m,
            t * (one - energy_step_fraction),
            t_log_min,
            t_inv_dlog,
            nt,
        )
        d_energy = (r_now - r_frac) / rho * type(t)(10.0)
        step = MATH_min3(
            d_face + wpv.boundary_overshoot(t),
            max_step_mm,
            MATH_max(d_energy, type(t)(1.0e-3)),
        )
        rho_path = rho * step * type(t)(0.1)
        t_new = wpv.energy_after_path(
            log_range, species, m, t, rho_path, t_log_min, t_inv_dlog, nt
        )
        de = MATH_max((t - t_new) * mass_number, zero)
        if straggling == 1 and t_new > zero:
            beta2 = wpv.beta_squared(t, mass_per_nucleon)
            z_eff = wpv.effective_charge(charge, wp.sqrt(beta2))
            var = wpv.straggling_variance_mev2(z_eff, z_over_a[m], beta2, rho_path)
            de, state = wpv.sample_energy_loss(state, de, var)
            de = MATH_min(de, t * mass_number)
        hinge = type(t)(1.0)
        if mcs == 1:
            hinge = type(t)(wp.randf(state))
        x0 = x
        y0 = y
        z0 = z
        xh = x + ux * step * hinge
        yh = y + uy * step * hinge
        zh = z + uz * step * hinge
        if mcs == 1:
            pv = wpv.pv_mev(t, mass_mev, mass_number)
            theta0_sq = wpv.mcs_theta0_squared(charge, pv, p1v1, rho_x_s[m], rho_path)
            u = MATH_max(type(t)(wp.randf(state)), type(t)(1.0e-30))
            theta = wp.sqrt(theta0_sq) * wp.sqrt(-type(t)(2.0) * wp.log(u))
            theta = MATH_min(theta, one)
            phi = type(t)(6.283185307179586) * type(t)(wp.randf(state))
            ct = wp.cos(theta)
            st = wp.sin(theta)
            nux = wpv.rotate_x(ux, uy, uz, ct, st, phi)
            nuy = wpv.rotate_y(ux, uy, uz, ct, st, phi)
            nuz = wpv.rotate_z(ux, uy, uz, ct, st, phi)
            norm = wp.sqrt(nux * nux + nuy * nuy + nuz * nuz)
            ux = nux / norm
            uy = nuy / norm
            uz = nuz / norm
        x1 = xh + ux * step * (one - hinge)
        y1 = yh + uy * step * (one - hinge)
        z1 = zh + uz * step * (one - hinge)
        wpv.deposit_segment(
            edep,
            sox,
            soy,
            soz,
            sdx,
            sdy,
            sdz,
            snx,
            sny,
            snz,
            x0,
            y0,
            z0,
            xh,
            yh,
            zh,
            de * hinge,
        )
        wpv.deposit_segment(
            edep,
            sox,
            soy,
            soz,
            sdx,
            sdy,
            sdz,
            snx,
            sny,
            snz,
            xh,
            yh,
            zh,
            x1,
            y1,
            z1,
            de * (one - hinge),
        )
        wp.atomic_add(acc, 0, wp.float64(de))
        t = t - de / mass_number
        x = x1
        y = y1
        z = z1
        steps = steps + 1
    steps_out[tid] = steps


@wp.func
def MATH_min(a: Any, b: Any) -> Any:
    return wp.min(a, b)


@wp.func
def MATH_max(a: Any, b: Any) -> Any:
    return wp.max(a, b)


@wp.func
def MATH_min3(a: Any, b: Any, c: Any) -> Any:
    return wp.min(wp.min(a, b), c)


def _device_name(backend: str) -> str:
    if backend == "warp-cpu":
        return "cpu"
    if backend == "warp-cuda":
        if not wp.is_cuda_available():
            raise UnsupportedConfigurationError(
                "backend 'warp-cuda' requested but no CUDA device is available"
            )
        return "cuda:0"
    raise UnsupportedConfigurationError(f"unknown Warp backend {backend!r}")


def run_warp(config: SimulationConfig, tables: TableSet) -> dict:
    """Run all batches on a Warp device; returns tallies, accounting and timing."""
    cfg = config.with_defaults()
    grid = cfg.scoring
    assert grid is not None
    device = _device_name(cfg.backend)
    np_dtype = np.float32 if cfg.precision == "float32" else np.float64
    wp_dtype = wp.float32 if cfg.precision == "float32" else wp.float64
    ph = cfg.physics
    sp_index = tables.species_index(cfg.source.species)
    geo = cfg.geometry
    t_cold = time.perf_counter()
    with wp.ScopedDevice(device):
        log_range = wp.array(tables.log_range.astype(np_dtype), dtype=wp_dtype)
        z_over_a = wp.array(tables.z_over_a.astype(np_dtype), dtype=wp_dtype)
        rho_x_s = wp.array(tables.rho_x_s.astype(np_dtype), dtype=wp_dtype)
        mat_index = wp.array(
            np.ascontiguousarray(geo.material_index, dtype=np.int32), dtype=wp.int32
        )
        density = wp.array(
            np.ascontiguousarray(geo.density_g_cm3, dtype=np_dtype), dtype=wp_dtype
        )
        edep = wp.zeros(grid.shape, dtype=wp_dtype)
        acc = wp.zeros(4, dtype=wp.float64)
    upload_seconds = time.perf_counter() - t_cold
    energy = Tally("energy", grid.shape, cfg.batches)
    accounting = EnergyAccounting()
    seeds = np.random.SeedSequence(cfg.seed).spawn(cfg.batches)
    per_batch = [
        cfg.histories // cfg.batches + (1 if b < cfg.histories % cfg.batches else 0)
        for b in range(cfg.batches)
    ]
    kernel_seconds = 0.0
    host_seconds = 0.0
    t_wall = time.perf_counter()
    first_launch_seconds = None
    offset = 0
    for b in range(cfg.batches):
        n = per_batch[b]
        rng = np.random.default_rng(seeds[b])
        prim = cfg.source.sample(n, rng)
        t_host = time.perf_counter()
        with wp.ScopedDevice(device):
            pos = wp.array(prim["position_mm"].astype(np_dtype), dtype=wp_dtype)
            dirs = wp.array(prim["direction"].astype(np_dtype), dtype=wp_dtype)
            en = wp.array(prim["energy_mev_per_u"].astype(np_dtype), dtype=wp_dtype)
            steps_out = wp.zeros(n, dtype=wp.int32)
            edep.zero_()
            acc.zero_()
            wp.synchronize_device(device)
            t0 = time.perf_counter()
            wp.launch(
                transport_kernel,
                dim=n,
                inputs=[
                    pos,
                    dirs,
                    en,
                    int(cfg.seed),
                    int(offset),
                    int(sp_index),
                    wp_dtype(float(tables.charge[sp_index])),
                    wp_dtype(float(tables.mass_number[sp_index])),
                    wp_dtype(float(tables.mass_mev[sp_index])),
                    log_range,
                    wp_dtype(tables.t_log_min),
                    wp_dtype(tables.t_inv_dlog),
                    int(tables.t_grid.size),
                    z_over_a,
                    rho_x_s,
                    mat_index,
                    density,
                    wp_dtype(geo.origin_mm[0]),
                    wp_dtype(geo.origin_mm[1]),
                    wp_dtype(geo.origin_mm[2]),
                    wp_dtype(geo.spacing_mm[0]),
                    wp_dtype(geo.spacing_mm[1]),
                    wp_dtype(geo.spacing_mm[2]),
                    int(geo.shape[0]),
                    int(geo.shape[1]),
                    int(geo.shape[2]),
                    edep,
                    wp_dtype(grid.origin_mm[0]),
                    wp_dtype(grid.origin_mm[1]),
                    wp_dtype(grid.origin_mm[2]),
                    wp_dtype(grid.spacing_mm[0]),
                    wp_dtype(grid.spacing_mm[1]),
                    wp_dtype(grid.spacing_mm[2]),
                    int(grid.shape[0]),
                    int(grid.shape[1]),
                    int(grid.shape[2]),
                    int(ph.straggling),
                    int(ph.multiple_scattering),
                    wp_dtype(ph.max_step_mm),
                    wp_dtype(ph.energy_step_fraction),
                    wp_dtype(ph.cutoff_mev_per_u),
                    acc,
                    steps_out,
                ],
            )
            wp.synchronize_device(device)
            dt = time.perf_counter() - t0
            if first_launch_seconds is None:
                first_launch_seconds = dt
            kernel_seconds += dt
            energy.current[...] = edep.numpy().astype(np.float64)
            a = acc.numpy().astype(np.float64)
            accounting.deposited_continuous += float(a[0])
            accounting.deposited_cutoff += float(a[1])
            accounting.escaped += float(a[2])
            accounting.truncated_max_steps += float(a[3])
            accounting.steps += int(steps_out.numpy().sum())
        accounting.initial += float(
            np.sum(prim["energy_mev_per_u"]) * tables.mass_number[sp_index]
        )
        accounting.histories += n
        energy.close_batch(n)
        host_seconds += time.perf_counter() - t_host - dt
        offset += n
    wall = time.perf_counter() - t_wall
    return {
        "energy": energy,
        "accounting": accounting,
        "histories_per_batch": per_batch,
        "wall_seconds": wall,
        "kernel_seconds": kernel_seconds,
        "host_seconds": host_seconds,
        "upload_seconds": upload_seconds,
        "first_launch_seconds": first_launch_seconds,
        "device": device,
        "rng": "warp rand_init(seed, batch_offset + history) per history; primaries from numpy SeedSequence(seed).spawn(batches)",
    }
