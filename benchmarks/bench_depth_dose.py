"""Benchmark: deterministic proton CSDA depth dose (Stage 6, decision 0030).

The first reproducible benchmark of the ``ionmc`` performance suite. It times the
core longitudinal transport kernel (150 MeV proton pencil beam, 0.5 mm bins, energy
straggling off so the workload is deterministic and its scientific output is a
fixed digest) on the reference, Warp CPU and Warp CUDA backends, and records:

* full run provenance (hardware, versions, devices);
* per-backend timing statistics and throughput (histories/s), with device
  synchronisation so the GPU is timed honestly;
* a **scientific digest** of the depth-dose curve and the cross-backend relative
  agreement, so an optimisation can be shown to leave the physics unchanged
  (milestone V6) and a benchmark that silently broke the physics fails.

The reference backend is a scalar Python oracle (~10^2 histories/s), so it is timed
at a small history count; the Warp backends are timed at a larger count. Throughput
is a per-history rate, so the numbers remain comparable (the GPU is under-utilised
at these sizes — a scaling sweep is a later task). Wall-clock time is recorded, never
asserted against an absolute threshold; only the physics is gated.

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
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine
from ionmc.transport.geometry import WaterSlab

# float32 budget for the deterministic cross-backend depth dose, from the
# established V1 depth-dose gate (decision 0009 / validation/v1_depth_dose_csda.py):
# the edge-aware cumulative_rel_diff, tight to 1e-4 reference-vs-Warp and 1e-5
# CUDA-vs-CPU. Comparisons are on per-history-normalised curves so a uniform-scale
# divergence is caught (not hidden by unit normalisation).
REF_VS_WARP_CUM = 1.0e-4
CUDA_VS_CPU_CUM = 1.0e-5
CUDA_VS_CPU_MAXBIN = 5.0e-3


def _build(cache_dir: str | None) -> tuple[TransportEngine, PencilBeamSource]:
    proton = load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, cache_dir))
    grid = DepthDoseGrid(250.0, 500)  # 0-25 cm, 0.5 mm bins
    eng = TransportEngine(
        proton, WaterSlab(250.0), grid, straggling=False, nuclear=False
    )
    return eng, PencilBeamSource(150.0)


def _per_history(edep: np.ndarray, n_histories: int) -> np.ndarray:
    """Per-history depth dose. For the deterministic (straggling-off) CSDA workload
    every history is identical, so this is directly comparable across backends timed
    at different history counts, while retaining absolute magnitude."""
    return edep / float(n_histories)


def _digest(per_hist: np.ndarray, grid: DepthDoseGrid) -> dict[str, Any]:
    """A reproducible scientific fingerprint of the per-history depth-dose curve
    (independent of the timing history count and carrying absolute magnitude)."""
    total = float(per_hist.sum())
    shape = per_hist / total if total > 0.0 else per_hist
    kpk = int(per_hist.argmax())
    return {
        "peak_depth_mm": float(grid.centers_mm[kpk]),
        "integral_per_history_mev": total,
        "shape_digest": array_digest(shape),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-cuda", action="store_true")
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--histories", type=int, default=20000, help="Warp workload")
    parser.add_argument("--ref-histories", type=int, default=200)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmup", type=int, default=1)
    args = parser.parse_args()

    report = new_report("depth_dose_csda")
    report["config"] = {
        "energy_mev": 150.0,
        "grid": "DepthDoseGrid(250.0, 500)",
        "straggling": False,
        "histories_warp": args.histories,
        "histories_reference": args.ref_histories,
        "repeats": args.repeats,
        "warmup": args.warmup,
        # the cross-backend gate is a float32-vs-float64 check: Warp accumulates in
        # float32, the reference oracle in float64 (decision 0009 budget)
        "reference_precision": "float64",
        "warp_precision": "float32",
    }
    try:
        eng, src = _build(args.cache_dir)
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        dump_report(report, sys.stdout)
        return 4

    grid = eng.grid
    # -- reference: small history count (scalar oracle), correctness digest -----
    # per-history dose retains absolute magnitude, so the cross-backend gate catches
    # a uniform-scale divergence (a units / normalisation / accumulation-constant bug)
    # that a unit-normalised comparison would hide.
    ref_res = eng.run(src, args.ref_histories, seed=1, path="python")
    ref_perhist = _per_history(ref_res.edep_mev, args.ref_histories)
    report["digest"] = {"reference": _digest(ref_perhist, grid)}
    report["cross_backend"] = {}

    ref_timing = measure(
        lambda: eng.run(src, args.ref_histories, seed=1, path="python"),
        label="depth_dose_csda",
        backend="reference",
        repeats=args.repeats,
        warmup=args.warmup,
        work_units=float(args.ref_histories),
    )
    report["timings"].append(timing_to_dict(ref_timing))

    gate_ok = True
    # -- Warp backends ---------------------------------------------------------
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
        for device in devices:
            n = args.histories
            res = eng.run(src, n, seed=1, path="warp", device=device)
            perhist = _per_history(res.edep_mev, n)
            label = backend_label("warp", device)
            report["digest"][label] = _digest(perhist, grid)
            agree = relative_agreement(ref_perhist, perhist)
            report["cross_backend"][f"{label}_vs_reference"] = agree
            ok = agree["cumulative_rel_diff"] <= REF_VS_WARP_CUM
            if device == "cpu":
                cpu_perhist = perhist
            elif cpu_perhist is not None:
                cagree = relative_agreement(cpu_perhist, perhist)
                report["cross_backend"][f"{label}_vs_cpu"] = cagree
                ok = (
                    ok
                    and cagree["cumulative_rel_diff"] <= CUDA_VS_CPU_CUM
                    and (cagree["max_bin_rel_diff"] <= CUDA_VS_CPU_MAXBIN)
                )
            gate_ok = gate_ok and ok
            timing = measure(
                lambda d=device, m=n: eng.run(src, m, seed=1, path="warp", device=d),
                label="depth_dose_csda",
                backend=label,
                repeats=args.repeats,
                warmup=args.warmup,
                sync=make_sync("warp", device),
                work_units=float(n),
            )
            report["timings"].append(timing_to_dict(timing))

        if args.require_cuda and not any(d != "cpu" for d in devices):
            gate_ok = False

    # -- speedups relative to the reference per-history rate -------------------
    rates = {t["backend"]: t["throughput_per_s"] for t in report["timings"]}
    ref_rate = rates.get("reference")
    if ref_rate:
        report["speedup_vs_reference"] = {
            b: (r / ref_rate) for b, r in rates.items() if r is not None
        }

    report["physics_gate_passed"] = bool(gate_ok)
    dump_report(report, sys.stdout)
    return 0 if gate_ok else 3


if __name__ == "__main__":
    raise SystemExit(main())
