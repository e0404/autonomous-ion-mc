"""V4 (uncertainty): batch-based statistical uncertainty for 3-D dose scoring.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0024 gates:

* ``se_scaling``          - the high-dose mean relative SEM scales as 1/sqrt(N)
                            (4x histories -> ~half the SEM), over a shared mask;
* ``mean_energy_conservation`` - the batch-mean total dose equals the per-batch
                            input energy for a contained beam;
* ``uncertainty_sanity``  - SEM > 0 in high-dose voxels, == 0 where no dose;
* ``warp_cpu_vs_reference`` - the batched-mean total dose agrees reference-vs-CPU;
* ``warp_cpu_vs_cuda``    - the batched-mean total dose agrees CUDA-vs-CPU.

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
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid


def _lat() -> DepthLateralGrid:
    return DepthLateralGrid(200.0, 400, 30.0, 60)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (41, 41, 400), (1.0, 1.0, 0.5), 1.0, origin_mm=(-20.5, -20.5, 0.0)
    )


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(41, 41, 100), origin_mm=(-20.5, -20.5, 0.0), spacing_mm=(1.0, 1.0, 2.0)
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

    # -- mean energy conservation (deterministic, reference) -----------------
    eng_det = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    rdet = eng_det.run_scattering_batched(
        src, lat, _dose(), 40, n_batches=4, seed=3, path="python"
    )
    expected = rdet.histories_per_batch * 150.0
    econ = abs(float(rdet.mean_dose3d_mev.sum()) - expected) / expected
    report["mean_energy_conservation_rel"] = econ
    gates["mean_energy_conservation"] = econ <= 1e-9

    # -- uncertainty sanity (stochastic, reference small) --------------------
    eng = TransportEngine(table, _box(), DepthDoseGrid(200.0, 10))  # scattering on
    rs = eng.run_scattering_batched(
        src, lat, _dose(), 400, n_batches=8, seed=5, path="python"
    )
    hi = rs.mean_dose3d_mev > 0.5 * rs.mean_dose3d_mev.max()
    sane = bool(
        np.all(rs.standard_error_mev[hi] > 0.0)
        and np.all(rs.standard_error_mev[rs.mean_dose3d_mev == 0.0] == 0.0)
    )
    report["uncertainty_sanity"] = {
        "mean_rel_uncertainty_hi": rs.mean_relative_uncertainty(0.5),
        "n_hi_voxels": int(hi.sum()),
    }
    gates["uncertainty_sanity"] = sane

    # -- Warp: 1/sqrt(N) scaling + cross-backend -----------------------------
    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        gates["se_scaling"] = False
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

        # 1/sqrt(N): 4x histories -> ~0.5x the high-dose mean relative SEM
        rN = eng.run_scattering_batched(
            src, lat, _dose(), 40000, n_batches=20, seed=100, path="warp", device="cpu"
        )
        r4 = eng.run_scattering_batched(
            src,
            lat,
            _dose(),
            160000,
            n_batches=20,
            seed=100,
            path="warp",
            device="cpu",
        )
        mask = r4.mean_dose3d_mev > 0.2 * r4.mean_dose3d_mev.max()
        uN = float(np.mean(rN.relative_standard_error[mask]))
        u4 = float(np.mean(r4.relative_standard_error[mask]))
        ratio = u4 / uN
        report["se_scaling"] = {
            "rel_sem_N": uN,
            "rel_sem_4N": u4,
            "ratio": ratio,
            "n_mask_voxels": int(mask.sum()),
        }
        gates["se_scaling"] = bool(0.35 <= ratio <= 0.71)

        # cross-backend: the batched-mean TOTAL dose agrees (per-voxel mean under
        # scattering decorrelates like dose, decision 0022; the total is tight)
        ref_b = eng.run_scattering_batched(
            src, lat, _dose(), 4000, n_batches=8, seed=9, path="python"
        )
        cpu_b = eng.run_scattering_batched(
            src, lat, _dose(), 4000, n_batches=8, seed=9, path="warp", device="cpu"
        )
        rt, ct = float(ref_b.mean_dose3d_mev.sum()), float(cpu_b.mean_dose3d_mev.sum())
        cpu_tot = abs(ct - rt) / rt
        report.setdefault("cross_backend", {})["cpu_vs_reference_total"] = cpu_tot
        gates["warp_cpu_vs_reference"] = cpu_tot <= 1e-4

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            ok = True
            d0 = float(cpu_b.mean_dose3d_mev.sum())
            for device in cuda_devices:
                cuda_b = eng.run_scattering_batched(
                    src,
                    lat,
                    _dose(),
                    4000,
                    n_batches=8,
                    seed=9,
                    path="warp",
                    device=device,
                )
                tot = abs(float(cuda_b.mean_dose3d_mev.sum()) - d0) / d0
                report.setdefault("cross_backend", {})[f"{device}_vs_cpu_total"] = tot
                ok = ok and tot <= 1e-4
            gates["warp_cpu_vs_cuda"] = ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
