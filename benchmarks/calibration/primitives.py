#!/usr/bin/env python3
"""Calibration primitives for freezing v2 performance targets (decision 0039).

This script measures what the workstation and NVIDIA Warp can deliver for the
*building blocks* of condensed-history ion transport, before any real physics
exists, so that the numeric targets in ``validation/release-plan.json`` are
derived from measurements rather than guessed or copied. It is NOT a physics
benchmark and its numbers are never release evidence by themselves.

Measured primitives (each repeated ``--repeats`` times, per device):

* ``init_seconds``: wall time of ``warp.init()`` in a fresh process.
* ``compile_seconds``: wall time of the first launch of a representative
  transport-like module after ``warp.clear_kernel_cache()`` (cold compile) and
  of a warm launch in a new process with the cache populated.
* ``step_kernel``: a synthetic condensed-history step loop per particle with
  the memory and arithmetic pattern of real transport — three RNG draws,
  a log-spaced table lookup with linear interpolation, an energy decrement, a
  direction rotation using sqrt/log/cos/sin, a voxel index computation and an
  atomic dose deposit into a 60 x 60 x 150 grid — in float32 and float64, with
  and without the atomic deposit, for 1e4..1e6 particles and a fixed number of
  steps per particle. Reported as steps per second and particle-steps per
  launch.
* ``memory``: attributable device memory of the particle state and scoring
  grid, and the host process peak resident set size.

Everything runs with only ``numpy`` and ``warp`` so the script executes
unchanged inside the controlled host runner.
"""

import argparse
import json
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

GRID = (60, 60, 150)
TABLE_SIZE = 256


def hardware_info() -> dict[str, Any]:
    cpu_model = None
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu_model = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    mem_total_kib = None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal"):
                mem_total_kib = int(line.split()[1])
                break
    except OSError:
        pass
    return {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "cpu_model": cpu_model,
        "logical_cpus": platform.os.cpu_count(),
        "host_memory_kib": mem_total_kib,
    }


def code_sha() -> dict[str, Any]:
    root = Path(__file__).resolve().parents[2]
    try:
        sha = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        dirty = bool(
            subprocess.check_output(
                ["git", "-C", str(root), "status", "--porcelain"],
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        )
        return {"git_sha": sha, "git_dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"git_sha": None, "git_dirty": None}


def build_kernel(wp, dtype):
    """Build the synthetic step kernel for one scalar type (float32/float64).

    Kernels are built lazily (the module imports without warp) and separately
    per precision through a closure so every literal is cast explicitly.
    """
    two_pi = dtype(6.283185307179586)

    @wp.func
    def table_lookup(
        table: wp.array(dtype=dtype), energy: dtype, e_min: dtype, inv_dlog: dtype
    ):
        # log-spaced lookup with linear interpolation, clamped to the table
        x = wp.log(energy / e_min) * inv_dlog
        x = wp.max(x, dtype(0.0))
        i = int(x)
        n = table.shape[0]
        if i >= n - 1:
            return table[n - 1]
        f = x - dtype(i)
        return table[i] * (dtype(1.0) - f) + table[i + 1] * f

    @wp.kernel
    def step_kernel(
        pos: wp.array(dtype=dtype, ndim=2),
        dirs: wp.array(dtype=dtype, ndim=2),
        energy: wp.array(dtype=dtype),
        table: wp.array(dtype=dtype),
        grid: wp.array3d(dtype=dtype),
        seed: int,
        n_steps: int,
        deposit: int,
        voxel_mm: dtype,
        e_min: dtype,
        inv_dlog: dtype,
    ):
        tid = wp.tid()
        state = wp.rand_init(seed, tid)
        e = energy[tid]
        x = pos[tid, 0]
        y = pos[tid, 1]
        z = pos[tid, 2]
        ux = dirs[tid, 0]
        uy = dirs[tid, 1]
        uz = dirs[tid, 2]
        for _ in range(n_steps):
            if e <= e_min:
                break
            u1 = dtype(wp.randf(state))
            u2 = dtype(wp.randf(state))
            u3 = dtype(wp.randf(state))
            stopping = table_lookup(table, e, e_min, inv_dlog)
            step = voxel_mm
            de = stopping * step * (dtype(1.0) + dtype(0.05) * (u1 - dtype(0.5)))
            if de > e:
                de = e
            # Gaussian-like small polar angle from two uniforms (Box-Muller)
            theta = dtype(0.02) * wp.sqrt(-dtype(2.0) * wp.log(u2 + dtype(1.0e-7)))
            phi = two_pi * u3
            sin_t = wp.sin(theta)
            cos_t = wp.cos(theta)
            # rotate the direction about a perpendicular axis
            if wp.abs(uz) < dtype(0.99):
                px = -uy
                py = ux
                pz = dtype(0.0)
            else:
                px = dtype(1.0)
                py = dtype(0.0)
                pz = dtype(0.0)
            pn = wp.sqrt(px * px + py * py + pz * pz)
            px = px / pn
            py = py / pn
            pz = pz / pn
            qx = uy * pz - uz * py
            qy = uz * px - ux * pz
            qz = ux * py - uy * px
            cp = wp.cos(phi)
            sp = wp.sin(phi)
            nx = cos_t * ux + sin_t * (cp * px + sp * qx)
            ny = cos_t * uy + sin_t * (cp * py + sp * qy)
            nz = cos_t * uz + sin_t * (cp * pz + sp * qz)
            nn = wp.sqrt(nx * nx + ny * ny + nz * nz)
            ux = nx / nn
            uy = ny / nn
            uz = nz / nn
            x = x + ux * step
            y = y + uy * step
            z = z + uz * step
            e = e - de
            ix = int((x + dtype(60.0)) / dtype(2.0))
            iy = int((y + dtype(60.0)) / dtype(2.0))
            iz = int(z / dtype(2.0))
            if ix >= 0 and ix < 60 and iy >= 0 and iy < 60 and iz >= 0 and iz < 150:
                if deposit == 1:
                    wp.atomic_add(grid, ix, iy, iz, de)
            else:
                break
        energy[tid] = e
        pos[tid, 0] = x
        pos[tid, 1] = y
        pos[tid, 2] = z
        dirs[tid, 0] = ux
        dirs[tid, 1] = uy
        dirs[tid, 2] = uz

    return step_kernel


def make_table(dtype_np) -> np.ndarray:
    # Proton-like linear stopping power in water, MeV/mm, shaped as S ~ E^-0.8 and
    # scaled so that S(150 MeV) ~ 0.55 MeV/mm; a 150 MeV particle then stops after
    # roughly 150 one-millimetre steps, exercising the early-exit paths.
    e = np.geomspace(0.5, 300.0, TABLE_SIZE)
    return (30.0 * e**-0.8).astype(dtype_np)


def run_step_measurement(
    wp, kernels, device, n_particles, n_steps, dtype, deposit, repeats
):
    dtype_np = np.float32 if dtype == "float32" else np.float64
    wp_dtype = wp.float32 if dtype == "float32" else wp.float64
    kernel = kernels[dtype]
    rng = np.random.default_rng(12345)
    pos_np = np.zeros((n_particles, 3), dtype=dtype_np)
    pos_np[:, 0] = rng.normal(0.0, 3.0, n_particles)
    pos_np[:, 1] = rng.normal(0.0, 3.0, n_particles)
    dir_np = np.zeros((n_particles, 3), dtype=dtype_np)
    dir_np[:, 2] = 1.0
    energy_np = np.full(n_particles, 150.0, dtype=dtype_np)
    table = wp.array(make_table(dtype_np), dtype=wp_dtype, device=device)
    grid = wp.zeros(GRID, dtype=wp_dtype, device=device)
    timings = []
    free_before = None
    if device.startswith("cuda"):
        free_before = wp.get_device(device).free_memory
    pos = wp.array(pos_np, dtype=wp_dtype, device=device)
    dirs = wp.array(dir_np, dtype=wp_dtype, device=device)
    energy = wp.array(energy_np, dtype=wp_dtype, device=device)
    device_bytes = None
    if free_before is not None:
        device_bytes = free_before - wp.get_device(device).free_memory
    e_min = wp_dtype(0.5)
    inv_dlog = wp_dtype((TABLE_SIZE - 1) / np.log(300.0 / 0.5))
    voxel = wp_dtype(1.0)
    for _ in range(repeats):
        pos.assign(pos_np)
        dirs.assign(dir_np)
        energy.assign(energy_np)
        grid.zero_()
        wp.synchronize_device(device)
        t0 = time.perf_counter()
        wp.launch(
            kernel,
            dim=n_particles,
            inputs=[
                pos,
                dirs,
                energy,
                table,
                grid,
                7,
                n_steps,
                deposit,
                voxel,
                e_min,
                inv_dlog,
            ],
            device=device,
        )
        wp.synchronize_device(device)
        timings.append(time.perf_counter() - t0)
    total_deposit = float(grid.numpy().sum()) if deposit else None
    energy_after = energy.numpy()
    steps_executed = int(n_particles * n_steps)  # upper bound; particles may stop early
    best = min(timings)
    return {
        "device": device,
        "dtype": dtype,
        "particles": n_particles,
        "steps_per_particle": n_steps,
        "deposit": bool(deposit),
        "seconds": timings,
        "best_seconds": best,
        "particle_steps_max": steps_executed,
        "particle_steps_per_second_upper": steps_executed / best,
        "mean_final_energy": float(energy_after.mean()),
        "grid_total_deposit": total_deposit,
        "attributable_device_bytes": device_bytes,
    }


def measure(args) -> dict[str, Any]:
    import warp as wp

    report: dict[str, Any] = {
        "schema_version": 1,
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "hardware": hardware_info(),
        "code": code_sha(),
        "warp_version": wp.__version__,
        "measurements": {},
    }
    if args.clear_cache:
        wp.clear_kernel_cache()
    t0 = time.perf_counter()
    wp.init()
    report["measurements"]["init_seconds"] = time.perf_counter() - t0
    devices = [d for d in args.devices if d == "cpu" or wp.is_cuda_available()]
    report["devices"] = {}
    for alias in devices:
        d = wp.get_device(alias)
        report["devices"][alias] = {
            "name": d.name,
            "arch": getattr(d, "arch", None),
            "total_memory_bytes": getattr(d, "total_memory", None),
        }
    kernel = {
        "float32": build_kernel(wp, wp.float32),
        "float64": build_kernel(wp, wp.float64),
    }
    report["measurements"]["compile_seconds"] = {}
    for alias in devices:
        # First launch triggers module load/compile for that device.
        t0 = time.perf_counter()
        run_step_measurement(wp, kernel, alias, 64, 4, "float32", 1, 1)
        run_step_measurement(wp, kernel, alias, 64, 4, "float64", 1, 1)
        report["measurements"]["compile_seconds"][alias] = time.perf_counter() - t0
    results = []
    for alias in devices:
        for n in args.particles:
            if alias == "cpu" and n > args.cpu_max_particles:
                continue
            for dtype in args.dtypes:
                for deposit in (1, 0):
                    results.append(
                        run_step_measurement(
                            wp,
                            kernel,
                            alias,
                            n,
                            args.steps,
                            dtype,
                            deposit,
                            args.repeats,
                        )
                    )
    report["measurements"]["step_kernel"] = results
    report["measurements"]["host_peak_rss_kib"] = resource.getrusage(
        resource.RUSAGE_SELF
    ).ru_maxrss
    return report


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--devices", nargs="+", default=["cpu", "cuda:0"])
    p.add_argument(
        "--particles", nargs="+", type=int, default=[10_000, 100_000, 1_000_000]
    )
    p.add_argument("--cpu-max-particles", type=int, default=100_000)
    p.add_argument("--steps", type=int, default=200)
    p.add_argument("--dtypes", nargs="+", default=["float32", "float64"])
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--clear-cache", action="store_true", help="measure a cold compile")
    p.add_argument("--output", type=Path)
    args = p.parse_args(argv)
    report = measure(args)
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
