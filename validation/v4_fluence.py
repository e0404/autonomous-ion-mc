"""V4 (fluence): energy-resolved fluence-spectrum scoring + lookup accumulation.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0026 gates:

* ``lookup_reproduces_offline`` - the on-the-fly lookup accumulator equals the
                            offline post-processing of the scored spectrum to
                            round-off (the V4 milestone gate);
* ``stopping_lookup_reproduces_energy`` - with w = S_lin the accumulator recovers
                            the total step energy deposited (deterministic mode);
* ``fluence_sanity``      - counts >= 0 and the spectrum support is physical;
* ``warp_cpu_vs_reference`` - the deterministic per-bin spectrum + accumulator
                            match the reference (float32 budget);
* ``warp_cpu_vs_cuda``    - the deterministic per-bin spectrum + accumulator agree
                            CUDA-vs-CPU.

Exit 0 if every gate passes, 3 otherwise, 4 if the dataset is missing.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from typing import Any

import numpy as np

from ionmc import __version__
from ionmc.backend import mathlib
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    FluenceSpectrum,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid


def _lat() -> DepthLateralGrid:
    return DepthLateralGrid(200.0, 400, 20.0, 40)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (11, 11, 400), (1.0, 1.0, 0.5), 1.0, origin_mm=(-5.5, -5.5, 0.0)
    )


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(11, 11, 100), origin_mm=(-5.5, -5.5, 0.0), spacing_mm=(1.0, 1.0, 2.0)
    )


def _jsonable(o: Any) -> Any:
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-cuda", action="store_true")
    parser.add_argument("--cache-dir", default=None)
    args = parser.parse_args()

    report: dict[str, Any] = {
        "schema_version": 1,
        "ionmc_version": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "gates": {},
    }
    gates = report["gates"]
    try:
        table = load_stopping_table(
            cache.load_path(MCSQUARE_PSTAR_WATER, args.cache_dir)
        )
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        json.dump(report, sys.stdout, indent=2, default=_jsonable)
        sys.stdout.write("\n")
        return 4

    src = PencilBeamSource(150.0)
    lat = _lat()
    fl = FluenceSpectrum()

    # -- exact lookup identity (scattering on, reference) --------------------
    eng = TransportEngine(table, _box(), DepthDoseGrid(200.0, 10))
    r = eng.run_scattering(
        src, lat, 200, seed=7, path="python", dose_grid=_dose(), fluence=fl
    )
    a_gate = float(r.fluence_lookup_sum)
    a_off = fl.postprocess_lookup(r.fluence_counts_mm, r.fluence_lookup_table)
    ident = abs(a_gate - a_off) / abs(a_off)
    report["lookup"] = {"a_gate": a_gate, "a_offline": a_off, "rel_diff": ident}
    gates["lookup_reproduces_offline"] = ident <= 1e-12

    # -- stopping-power lookup reproduces energy (deterministic) -------------
    eng_det = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    rd = eng_det.run_scattering(
        src, lat, 200, seed=5, path="python", dose_grid=_dose(), fluence=fl
    )
    dose_total = float(rd.dose3d_mev.sum())
    ratio = float(rd.fluence_lookup_sum) / dose_total
    report["energy"] = {
        "a_gate": float(rd.fluence_lookup_sum),
        "dose_total": dose_total,
        "ratio": ratio,
    }
    gates["stopping_lookup_reproduces_energy"] = 0.99 <= ratio <= 1.0

    # -- fluence sanity ------------------------------------------------------
    counts = rd.fluence_counts_mm
    occ = np.nonzero(counts > 0.0)[0]
    sane = bool(
        np.all(counts >= 0.0)
        and fl.centers_mev[occ].max() <= 150.0
        and fl.centers_mev[occ].min() >= 0.0
        and counts.sum() > 0.0
    )
    report["fluence"] = {
        "total_track_length_mm": float(counts.sum()),
        "support_lo_mev": float(fl.centers_mev[occ].min()),
        "support_hi_mev": float(fl.centers_mev[occ].max()),
    }
    gates["fluence_sanity"] = sane

    # -- Warp cross-backend (deterministic per-bin) --------------------------
    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        gates["warp_cpu_vs_reference"] = False
        if args.require_cuda:
            gates["warp_cpu_vs_cuda"] = False
    else:
        wp = mathlib.warp_module()
        wp.config.log_level = wp.LOG_WARNING
        wp.init()
        devices = [dv.alias for dv in wp.get_devices()]
        report["warp"] = {
            "available": True,
            "version": wp.config.version,
            "devices": devices,
        }
        ref_d = eng_det.run_scattering(
            src, lat, 50, seed=4, path="python", dose_grid=_dose(), fluence=fl
        )
        cpu_d = eng_det.run_scattering(
            src,
            lat,
            50,
            seed=4,
            path="warp",
            device="cpu",
            dose_grid=_dose(),
            fluence=fl,
        )
        peak = float(ref_d.fluence_counts_mm.max())
        cpu_bin = float(
            np.max(np.abs(cpu_d.fluence_counts_mm - ref_d.fluence_counts_mm)) / peak
        )
        cpu_acc = abs(cpu_d.fluence_lookup_sum - ref_d.fluence_lookup_sum) / abs(
            ref_d.fluence_lookup_sum
        )
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "max_bin_rel_diff": cpu_bin,
            "accumulator_rel_diff": cpu_acc,
        }
        gates["warp_cpu_vs_reference"] = cpu_bin <= 5e-3 and cpu_acc <= 5e-3

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            ok = True
            for device in cuda_devices:
                gpu_d = eng_det.run_scattering(
                    src,
                    lat,
                    50,
                    seed=4,
                    path="warp",
                    device=device,
                    dose_grid=_dose(),
                    fluence=fl,
                )
                gbin = float(
                    np.max(np.abs(gpu_d.fluence_counts_mm - cpu_d.fluence_counts_mm))
                    / peak
                )
                gacc = abs(gpu_d.fluence_lookup_sum - cpu_d.fluence_lookup_sum) / abs(
                    cpu_d.fluence_lookup_sum
                )
                cb = report.setdefault("cross_backend", {})
                cb[f"{device}_vs_cpu"] = {
                    "max_bin_rel_diff": gbin,
                    "accumulator_rel_diff": gacc,
                }
                ok = ok and gbin <= 5e-3 and gacc <= 5e-3
            gates["warp_cpu_vs_cuda"] = ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
