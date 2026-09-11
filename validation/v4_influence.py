"""V4 (influence): beamlet-resolved scoring and sparse dose-influence matrices.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0022 gates:

* ``beamlet_sum_equals_broadfield`` - the summed influence-matrix dose equals the
                            batched broad-field dose to round-off (the first V4
                            gate: sum of beamlet doses == broad-field dose);
* ``sparse_vs_dense``     - a 1 % thresholded matrix keeps >= 99 % of the energy;
* ``energy_conservation`` - the summed influence dose equals the deposited energy;
* ``warp_cpu_vs_reference`` - the CPU broad-field dose matches the reference;
* ``warp_cpu_vs_cuda``    - the CUDA broad-field dose agrees with CPU.

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
from ionmc.influence import assemble_influence_matrix
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid

N_LARGE = 20000


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (120, 120, 250), (1.0, 1.0, 1.0), 1.0, origin_mm=(-60.0, -60.0, 0.0)
    )


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(40, 40, 125), origin_mm=(-60.0, -60.0, 0.0), spacing_mm=(3.0, 3.0, 2.0)
    )


def _lat() -> DepthLateralGrid:
    return DepthLateralGrid(250.0, 500, 30.0, 400)


def _beamlets() -> list[PencilBeamSource]:
    return [
        PencilBeamSource(150.0, position_mm=(x, 0.0, 0.0), beamlet=b)
        for b, x in enumerate([-20.0, -7.0, 7.0, 20.0])
    ]


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

    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    srcs = _beamlets()
    lat = _lat()

    # -- beamlet sum == broad field (reference) ------------------------------
    m = assemble_influence_matrix(
        eng, srcs, lat, _dose(), n_histories=400, seed=1000, path="python"
    )
    combined = eng.run_scattering_multi(
        srcs, lat, n_histories=400, seed=1000, path="python", dose_grid=_dose()
    )
    tot = m.total_dose()
    cb = combined.dose3d_mev
    denom = float(cb.sum())
    sum_diff = float(np.max(np.abs(tot - cb)) / denom)
    report["beamlet_sum"] = {
        "max_rel_diff": sum_diff,
        "influence_total_mev": float(tot.sum()),
        "broadfield_total_mev": denom,
        "n_beamlets": m.n_beamlets,
        "nnz": m.nnz,
    }
    gates["beamlet_sum_equals_broadfield"] = sum_diff <= 1e-9

    # -- sparse vs dense (thresholding) --------------------------------------
    mt = assemble_influence_matrix(
        eng,
        srcs,
        lat,
        _dose(),
        n_histories=400,
        seed=1000,
        path="python",
        threshold_frac=0.01,
    )
    kept = float(mt.total_dose().sum()) / float(tot.sum())
    report["sparse"] = {
        "full_nnz": m.nnz,
        "thresholded_nnz": mt.nnz,
        "energy_kept": kept,
    }
    gates["sparse_vs_dense"] = mt.nnz <= m.nnz and kept >= 0.99

    # -- energy conservation -------------------------------------------------
    econ = abs(float(tot.sum()) - float(400 * len(srcs) * 150.0)) / float(
        400 * len(srcs) * 150.0
    )
    report["energy_conservation_rel"] = econ
    gates["energy_conservation"] = econ <= 1e-9

    # -- Warp cross-backend --------------------------------------------------
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
        # Cross-backend broad-field parity is checked DETERMINISTICALLY (scattering
        # and straggling off): with scattering on, float32-vs-float64 DDA
        # face-flip decorrelation makes the per-voxel dose two independent MC
        # estimates (large per-voxel diff at finite N), so the deterministic run
        # isolates the float32 arithmetic -- the tight, discriminating spatial
        # metric (the stochastic DoseGrid3D cross-backend was validated in
        # DEV-016). The total energy stays tight either way.
        eng_det = TransportEngine(
            table, _box(), DepthDoseGrid(250.0, 10), straggling=False, scattering=False
        )
        ref_bf = eng_det.run_scattering_multi(
            srcs, lat, n_histories=1, seed=5, path="python", dose_grid=_dose()
        )
        cpu_bf = eng_det.run_scattering_multi(
            srcs,
            lat,
            n_histories=1,
            seed=5,
            path="warp",
            device="cpu",
            dose_grid=_dose(),
        )
        rd, cd = ref_bf.dose3d_mev, cpu_bf.dose3d_mev
        cpu_total = abs(float(cd.sum()) - float(rd.sum())) / float(rd.sum())
        cpu_max = float(np.max(np.abs(cd - rd)) / rd.max())
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "total_rel_diff": cpu_total,
            "max_voxel_rel_diff": cpu_max,
        }
        gates["warp_cpu_vs_reference"] = cpu_total <= 1e-5 and cpu_max <= 5e-3

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            cpu_l = eng.run_scattering_multi(
                srcs,
                lat,
                n_histories=N_LARGE,
                seed=9,
                path="warp",
                device="cpu",
                dose_grid=_dose(),
            )
            ok = True
            for device in cuda_devices:
                cuda_l = eng.run_scattering_multi(
                    srcs,
                    lat,
                    n_histories=N_LARGE,
                    seed=9,
                    path="warp",
                    device=device,
                    dose_grid=_dose(),
                )
                d0 = float(cpu_l.dose3d_mev.sum())
                tot_d = abs(float(cuda_l.dose3d_mev.sum()) - d0) / d0
                report.setdefault("cross_backend", {})[f"{device}_vs_cpu_total"] = tot_d
                ok = ok and tot_d <= 1e-4
            gates["warp_cpu_vs_cuda"] = ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
