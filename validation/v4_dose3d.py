"""V4 (dose scoring): lab-frame 3-D dose scoring on the voxel-grid path.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0021 gates:

* ``energy_conservation`` - a contained beam deposits all its energy into the 3-D
                            dose grid (dose sum == deposited energy);
* ``depth_dose_consistency`` - the 3-D dose z-marginal reproduces the beam-frame
                            depth-dose R80 to within a dose voxel;
* ``grid_independence``   - the total 3-D dose is invariant to the dose grid's
                            resolution and alignment;
* ``warp_cpu_vs_reference`` - the CPU dose grid matches the reference grid driver;
* ``warp_cpu_vs_cuda``    - the CUDA dose grid agrees with CPU (same precision).

Exit 0 if every gate passes, 3 otherwise, 4 if the dataset is missing.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from typing import Any

import numpy as np

from ionmc import __version__, materials, particles
from ionmc.backend import mathlib
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid

REF_WARP_DD_TOL = 5e-4
N_LARGE = 40000


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(depth, int(depth * 2), 30.0, 400)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (120, 120, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-60.0, -60.0, 0.0)
    )


def _dose(nzv: int = 150, sp: float = 2.0) -> DoseGrid3D:
    n = int(120 / sp)
    return DoseGrid3D(
        shape=(n, n, nzv),
        origin_mm=(-60.0, -60.0, 0.0),
        spacing_mm=(sp, sp, 300.0 / nzv),
    )


def _r80(centers: np.ndarray, prof: np.ndarray) -> float:
    i = int(prof.argmax())
    lvl = 0.8 * float(prof.max())
    for k in range(i, len(prof) - 1):
        if prof[k] >= lvl >= prof[k + 1]:
            f = (prof[k] - lvl) / (prof[k] - prof[k + 1])
            return float(centers[k] + f * (centers[k + 1] - centers[k]))
    return float(centers[i])


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

    _ = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, "numpy")
    src = PencilBeamSource(150.0)
    lat = _lat(300.0)
    det = dict(straggling=False, scattering=False)

    # -- energy conservation (deterministic + scattering) --------------------
    eng_det = TransportEngine(table, _box(), DepthDoseGrid(300.0, 10), **det)
    rd = eng_det.run_scattering(src, lat, 1, seed=4, path="python", dose_grid=_dose())
    eng = TransportEngine(table, _box(), DepthDoseGrid(300.0, 10))
    rs = eng.run_scattering(src, lat, 200, seed=7, path="python", dose_grid=_dose())
    econ_det = abs(float(rd.dose3d_mev.sum()) - rd.energy_deposited_mev)
    econ_sc = abs(float(rs.dose3d_mev.sum()) - rs.energy_deposited_mev)
    report["energy_conservation"] = {"deterministic": econ_det, "scattering": econ_sc}
    gates["energy_conservation"] = econ_det <= 1e-9 and econ_sc <= 1e-9

    # -- depth-dose consistency vs beam-frame --------------------------------
    dose = _dose(nzv=300, sp=2.0)  # 1 mm z-voxels
    rc = eng_det.run_scattering(src, lat, 1, seed=4, path="python", dose_grid=dose)
    zc, prof = dose.axis_marginals(rc.dose3d_mev)
    r80_3d = _r80(zc, prof)
    r80_bf = _r80(rc.grid.depth_centers_mm, rc.depth_dose_mev)
    report["depth_dose"] = {"r80_3d_mm": r80_3d, "r80_beamframe_mm": r80_bf}
    gates["depth_dose_consistency"] = abs(r80_3d - r80_bf) <= 1.0

    # -- grid independence ---------------------------------------------------
    a = eng_det.run_scattering(
        src, lat, 1, seed=4, path="python", dose_grid=_dose(150, 2.0)
    )
    b = eng_det.run_scattering(
        src,
        lat,
        1,
        seed=4,
        path="python",
        dose_grid=DoseGrid3D(
            shape=(30, 30, 100),
            origin_mm=(-45.0, -45.0, -10.0),
            spacing_mm=(3.0, 3.0, 3.0),
        ),
    )
    ta, tb = float(a.dose3d_mev.sum()), float(b.dose3d_mev.sum())
    gi = abs(ta - tb) / ta
    report["grid_independence"] = {"total_rel_diff": gi}
    gates["grid_independence"] = gi <= 1e-12

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
        # deterministic dose grid: CPU vs reference (float32 point deposition)
        ref_d = eng_det.run_scattering(
            src, lat, 1, seed=4, path="python", dose_grid=_dose()
        )
        cpu_d = eng_det.run_scattering(
            src, lat, 1, seed=4, path="warp", device="cpu", dose_grid=_dose()
        )
        denom = float(ref_d.dose3d_mev.sum())
        cpu_total = abs(float(cpu_d.dose3d_mev.sum()) - denom) / denom
        cpu_max = float(
            np.max(np.abs(cpu_d.dose3d_mev - ref_d.dose3d_mev)) / ref_d.dose3d_mev.max()
        )
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "total_rel_diff": cpu_total,
            "max_voxel_rel_diff": cpu_max,
        }
        gates["warp_cpu_vs_reference"] = cpu_total <= 1e-5 and cpu_max <= 5e-3

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            cpu_l = eng.run_scattering(
                src, lat, N_LARGE, seed=9, path="warp", device="cpu", dose_grid=_dose()
            )
            ok = True
            for device in cuda_devices:
                cuda_l = eng.run_scattering(
                    src,
                    lat,
                    N_LARGE,
                    seed=9,
                    path="warp",
                    device=device,
                    dose_grid=_dose(),
                )
                d0 = float(cpu_l.dose3d_mev.sum())
                tot = abs(float(cuda_l.dose3d_mev.sum()) - d0) / d0
                report.setdefault("cross_backend", {})[f"{device}_vs_cpu_total"] = tot
                ok = ok and tot <= 1e-4
            gates["warp_cpu_vs_cuda"] = ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
