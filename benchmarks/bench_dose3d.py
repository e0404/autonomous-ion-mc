"""Benchmark: deterministic 3-D voxel-grid dose (Stage 6, decision 0030).

The representative treatment-planning workload: a proton pencil beam transported
through a 3-D voxel geometry with ray/voxel DDA traversal, depositing into an
independent lab-frame 3-D dose grid by per-step atomic scatter (decisions 0020,
0021). Scattering and straggling are off, so the workload is deterministic and its
per-history 3-D dose is a fixed digest that is gated cross-backend to the established
float32 budget (decision 0021 / validation/v4_dose3d.py). This kernel is
expected to be memory-bound (atomic voxel accumulation), in contrast to the more
compute-bound 1-D depth dose of `bench_depth_dose.py` (a characterisation to be
confirmed by the deferred occupancy/roofline profiling, not asserted here).

Beyond the fixed-size timing, this driver runs a **history-count scaling sweep** on
the Warp backends, reporting throughput (histories/s) at each size — characterising
GPU utilisation and saturation, and quantifying the under-utilisation caveat noted
for the small-workload 1-D benchmark. Wall-clock time is recorded, never asserted;
only the physics is gated.

Exit 0 if the cross-backend physics gate passes, 3 otherwise, 4 if the dataset is
missing.
"""

from __future__ import annotations

import argparse
import sys
from typing import Any

import numpy as np

from ionmc.backend import mathlib
from ionmc.benchmarking import (
    array_digest,
    backend_label,
    dump_report,
    make_sync,
    measure,
    new_report,
    relative_agreement,
    timing_to_dict,
)
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid

# deterministic 3-D dose float32 budget (decision 0021 / validation/v4_dose3d.py):
# per-history dose, reference-vs-Warp and CUDA-vs-CPU tight on the integral and the
# worst voxel. Comparisons carry absolute magnitude (never unit-normalised).
REF_VS_WARP_TOTAL = 1.0e-5
REF_VS_WARP_MAXVOXEL = 5.0e-3
CUDA_VS_CPU_TOTAL = 1.0e-5
CUDA_VS_CPU_MAXVOXEL = 5.0e-3

# Warp scaling-sweep history counts (CPU runs serially, so it sweeps smaller sizes).
CUDA_SWEEP = (1_000, 10_000, 100_000, 1_000_000)
CPU_SWEEP = (1_000, 5_000, 20_000)


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(60, 60, 150), origin_mm=(-60.0, -60.0, 0.0), spacing_mm=(2.0, 2.0, 2.0)
    )


def _build(
    cache_dir: str | None,
) -> tuple[TransportEngine, PencilBeamSource, DepthLateralGrid]:
    table = load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, cache_dir))
    box = VoxelGrid3D.uniform(
        (120, 120, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-60.0, -60.0, 0.0)
    )
    eng = TransportEngine(
        table, box, DepthDoseGrid(300.0, 10), straggling=False, scattering=False
    )
    return eng, PencilBeamSource(150.0), DepthLateralGrid(300.0, 600, 30.0, 400)


def _run(
    eng: TransportEngine,
    src: PencilBeamSource,
    lat: DepthLateralGrid,
    n: int,
    path: str,
    device: str,
) -> np.ndarray:
    """One 3-D dose run, returning the per-history dose (fresh dose grid each call so
    accumulation never leaks between runs)."""
    res = eng.run_scattering(
        src, lat, n, seed=1, path=path, device=device, dose_grid=_dose()
    )
    return np.asarray(res.dose3d_mev, dtype=np.float64) / float(n)


def _digest(per_hist: np.ndarray) -> dict[str, Any]:
    return {
        "integral_per_history_mev": float(per_hist.sum()),
        "peak_voxel_mev": float(per_hist.max()),
        "n_nonzero_voxels": int(np.count_nonzero(per_hist)),
        "shape_digest": array_digest(per_hist.ravel()),
    }


def _sweep(
    eng: TransportEngine,
    src: PencilBeamSource,
    lat: DepthLateralGrid,
    path: str,
    device: str,
    sizes: tuple[int, ...],
    repeats: int,
    warmup: int,
) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    label = backend_label(path, device)
    sync = make_sync(path, device)
    for n in sizes:
        timing = measure(
            lambda m=n: eng.run_scattering(
                src, lat, m, seed=1, path=path, device=device, dose_grid=_dose()
            ),
            label="dose3d",
            backend=label,
            repeats=repeats,
            warmup=warmup,
            sync=sync,
            work_units=float(n),
        )
        row = timing_to_dict(timing)
        row["histories"] = n
        out.append(row)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-cuda", action="store_true")
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--corr-histories", type=int, default=1000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    args = parser.parse_args()

    report = new_report("dose3d")
    report["config"] = {
        "energy_mev": 150.0,
        "voxel_box": "VoxelGrid3D.uniform((120,120,300),1mm)",
        "dose_grid": "DoseGrid3D((60,60,150),2mm)",
        "straggling": False,
        "scattering": False,
        "corr_histories": args.corr_histories,
        "cuda_sweep": list(CUDA_SWEEP),
        "cpu_sweep": list(CPU_SWEEP),
        "repeats": args.repeats,
        "warmup": args.warmup,
        "reference_precision": "float64",
        "warp_precision": "float32",
    }
    try:
        eng, src, lat = _build(args.cache_dir)
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        dump_report(report, sys.stdout)
        return 4

    # -- correctness: per-history deterministic 3-D dose, cross-backend gate ----
    ref_perhist = _run(eng, src, lat, 1, "python", "cpu")
    report["digest"] = {"reference": _digest(ref_perhist)}
    report["cross_backend"] = {}

    gate_ok = True
    report["scaling"] = []
    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        if args.require_cuda:
            gate_ok = False
    else:
        wp = mathlib.warp_module()
        wp.config.log_level = wp.LOG_WARNING
        wp.init()
        devices = [dv.alias for dv in wp.get_devices()]
        report["warp"] = {
            "available": True,
            "version": str(wp.config.version),
            "devices": devices,
        }
        cpu_perhist: np.ndarray | None = None
        # process "cpu" before any CUDA device so the CUDA-vs-CPU gate always has its
        # baseline, regardless of the order wp.get_devices() returns
        for device in sorted(devices, key=lambda d: 0 if d == "cpu" else 1):
            perhist = _run(eng, src, lat, args.corr_histories, "warp", device)
            label = backend_label("warp", device)
            report["digest"][label] = _digest(perhist)
            agree = relative_agreement(ref_perhist.ravel(), perhist.ravel())
            report["cross_backend"][f"{label}_vs_reference"] = agree
            ok = (
                agree["total_rel_diff"] <= REF_VS_WARP_TOTAL
                and agree["max_bin_rel_diff"] <= REF_VS_WARP_MAXVOXEL
            )
            if device == "cpu":
                cpu_perhist = perhist
            elif cpu_perhist is not None:
                cagree = relative_agreement(cpu_perhist.ravel(), perhist.ravel())
                report["cross_backend"][f"{label}_vs_cpu"] = cagree
                ok = ok and (
                    cagree["total_rel_diff"] <= CUDA_VS_CPU_TOTAL
                    and cagree["max_bin_rel_diff"] <= CUDA_VS_CPU_MAXVOXEL
                )
            gate_ok = gate_ok and ok

        # -- scaling sweep -----------------------------------------------------
        for device in devices:
            sizes = CPU_SWEEP if device == "cpu" else CUDA_SWEEP
            report["scaling"].extend(
                _sweep(eng, src, lat, "warp", device, sizes, args.repeats, args.warmup)
            )

        if args.require_cuda and not any(d != "cpu" for d in devices):
            gate_ok = False

    # -- per-backend peak throughput (max over the sweep) ----------------------
    peak: dict[str, float] = {}
    for row in report["scaling"]:
        b = row["backend"]
        if row["throughput_per_s"] is not None:
            peak[b] = max(peak.get(b, 0.0), row["throughput_per_s"])
    report["peak_throughput_per_s"] = peak

    report["physics_gate_passed"] = bool(gate_ok)
    dump_report(report, sys.stdout)
    return 0 if gate_ok else 3


if __name__ == "__main__":
    raise SystemExit(main())
