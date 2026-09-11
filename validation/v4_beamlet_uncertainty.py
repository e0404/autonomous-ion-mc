"""V4 (beamlet uncertainty): beamlet-resolved statistical uncertainty.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0025 gates:

* ``sigma_alignment``     - data_sigma is aligned with data (same nnz), finite,
                            and >= 0;
* ``per_beamlet_se_scaling`` - a beamlet's high-dose mean relative SEM scales as
                            1/sqrt(N);
* ``beamlet_sum_within_statistics`` - the summed per-beamlet mean dose agrees with
                            a high-statistics broad-field dose within the combined
                            statistical uncertainty (the V4 "within statistics"
                            gate);
* ``warp_cpu_vs_reference`` - the DETERMINISTIC per-beamlet mean dose agrees
                            reference-vs-CPU per voxel (float32 budget);
* ``warp_cpu_vs_cuda``    - the deterministic per-beamlet mean dose agrees
                            CUDA-vs-CPU per voxel.

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
from ionmc.influence import assemble_influence_matrix_batched
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid


def _lat() -> DepthLateralGrid:
    return DepthLateralGrid(250.0, 500, 30.0, 400)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (60, 60, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-30.0, -30.0, 0.0)
    )


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(30, 30, 125), origin_mm=(-30.0, -30.0, 0.0), spacing_mm=(2.0, 2.0, 2.0)
    )


def _beamlets() -> list[PencilBeamSource]:
    return [
        PencilBeamSource(150.0, position_mm=(x, 0.0, 0.0), beamlet=b)
        for b, x in enumerate([-8.0, 0.0, 8.0])
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

    srcs = _beamlets()
    lat = _lat()
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))  # scattering on
    eng_det = TransportEngine(
        table, _box(), DepthDoseGrid(250.0, 10), straggling=False, scattering=False
    )

    if not mathlib.HAVE_WARP:
        # the estimator itself is backend-agnostic; run the reference-only gates
        m = assemble_influence_matrix_batched(
            eng, srcs, lat, _dose(), 400, n_batches=8, seed=1000, path="python"
        )
        gates["sigma_alignment"] = bool(
            m.data_sigma is not None
            and m.data_sigma.shape == m.data.shape
            and np.all(np.isfinite(m.data_sigma))
            and np.all(m.data_sigma >= 0.0)
        )
        report["warp"] = {"available": False}
        for g in (
            "per_beamlet_se_scaling",
            "beamlet_sum_within_statistics",
            "warp_cpu_vs_reference",
        ):
            gates[g] = False
        if args.require_cuda:
            gates["warp_cpu_vs_cuda"] = False
        report["all_gates_passed"] = all(gates.values())
        json.dump(report, sys.stdout, indent=2, default=_jsonable)
        sys.stdout.write("\n")
        return 0 if report["all_gates_passed"] else 3

    wp = mathlib.warp_module()
    wp.config.log_level = wp.LOG_WARNING
    wp.init()
    devices = [dv.alias for dv in wp.get_devices()]
    report["warp"] = {
        "available": True,
        "version": wp.config.version,
        "devices": devices,
    }

    # -- sigma alignment + per-beamlet SE scaling ----------------------------
    n_batches = 16
    mN = assemble_influence_matrix_batched(
        eng,
        srcs,
        lat,
        _dose(),
        8000,
        n_batches=n_batches,
        seed=100,
        path="warp",
        device="cpu",
    )
    m4 = assemble_influence_matrix_batched(
        eng,
        srcs,
        lat,
        _dose(),
        32000,
        n_batches=n_batches,
        seed=100,
        path="warp",
        device="cpu",
    )
    assert mN.data_sigma is not None
    gates["sigma_alignment"] = bool(
        mN.data_sigma.shape == mN.data.shape
        and np.all(np.isfinite(mN.data_sigma))
        and np.all(mN.data_sigma >= 0.0)
    )
    uN = mN.beamlet_relative_uncertainty(1, 0.2)
    u4 = m4.beamlet_relative_uncertainty(1, 0.2)
    report["per_beamlet"] = {"rel_sem_N": uN, "rel_sem_4N": u4, "ratio": u4 / uN}
    gates["per_beamlet_se_scaling"] = bool(0.35 <= u4 / uN <= 0.71)

    # -- beamlet sum within statistics ---------------------------------------
    # summed per-beamlet means (per-history normalised) vs a high-statistics
    # broad-field dose; check per-voxel agreement within the combined SEM.
    per_batch = 8000 // n_batches
    d_beam = mN.total_dose() / per_batch  # per-history-per-beamlet mean dose
    sig_beam = mN.total_sigma() / per_batch
    n_hi = 200000
    hi = eng.run_scattering_multi(
        srcs,
        lat,
        n_histories=n_hi,
        seed=777,
        path="warp",
        device="cpu",
        dose_grid=_dose(),
    )
    d_ref = hi.dose3d_mev / n_hi  # high-stat, negligible variance
    mask = d_ref > 0.2 * d_ref.max()
    z = np.abs(d_beam[mask] - d_ref[mask]) / np.maximum(sig_beam[mask], 1e-30)
    coverage_4sig = float(np.mean(z <= 4.0))
    report["within_statistics"] = {
        "coverage_within_4sigma": coverage_4sig,
        "max_z": float(z.max()),
        "n_hi_voxels": int(mask.sum()),
    }
    gates["beamlet_sum_within_statistics"] = coverage_4sig >= 0.95

    # -- cross-backend: deterministic per-beamlet mean per voxel -------------
    # compare dense per-beamlet doses (robust to float32/float64 nonzero-set
    # differences that would misalign the raw CSR data arrays)
    md_ref = assemble_influence_matrix_batched(
        eng_det, srcs, lat, _dose(), 8, n_batches=4, seed=5, path="python"
    )
    md_cpu = assemble_influence_matrix_batched(
        eng_det, srcs, lat, _dose(), 8, n_batches=4, seed=5, path="warp", device="cpu"
    )
    rr = np.stack([md_ref.beamlet_dose_flat(b) for b in range(md_ref.n_beamlets)])
    cc = np.stack([md_cpu.beamlet_dose_flat(b) for b in range(md_cpu.n_beamlets)])
    peak = float(rr.max())
    cpu_max = float(np.max(np.abs(cc - rr)) / peak)
    report.setdefault("cross_backend", {})["cpu_vs_reference_max_voxel"] = cpu_max
    gates["warp_cpu_vs_reference"] = cpu_max <= 5e-3

    cuda_devices = [dv for dv in devices if dv != "cpu"]
    if args.require_cuda and not cuda_devices:
        gates["warp_cpu_vs_cuda"] = False
    elif cuda_devices:
        ok = True
        for device in cuda_devices:
            md_gpu = assemble_influence_matrix_batched(
                eng_det,
                srcs,
                lat,
                _dose(),
                8,
                n_batches=4,
                seed=5,
                path="warp",
                device=device,
            )
            gg = np.stack(
                [md_gpu.beamlet_dose_flat(b) for b in range(md_gpu.n_beamlets)]
            )
            dmax = float(np.max(np.abs(gg - cc)) / peak)
            report["cross_backend"][f"{device}_vs_cpu_max_voxel"] = dmax
            ok = ok and dmax <= 5e-3
        gates["warp_cpu_vs_cuda"] = ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
