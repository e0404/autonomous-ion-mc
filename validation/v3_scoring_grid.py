"""V3 (scoring grid): grid-independence of integral dose under scoring-grid
resolution and alignment changes (decision 0019; closes milestone V3).

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached). Emits one JSON document with these gates:

* ``integral_dose_resolution`` - straight ray (scattering off): total deposited
                               energy equal on a coarse and a 4x-finer grid, and
                               the fine per-bin dose summed 4:1 equals the coarse;
* ``integral_dose_alignment`` - total deposited energy invariant under a depth-
                               origin and lateral-centre shift of the grid;
* ``lateral_origin``         - scattering off: a lateral window centred on the
                               beam captures the dose, one shifted off it captures
                               none, and the energy-weighted mean lateral position
                               is the beam axis -- a *discriminating* check of the
                               lateral origin (sigma_x is translation-invariant);
* ``lateral_shift_sigma_x``  - scattering on: sigma_x(z) invariant under a
                               lateral-centre shift (an invariance, not a
                               discriminating check);
* ``partial_coverage``       - a grid starting past the entrance captures strictly
                               less energy (the depth origin shifts the window);
* ``warp_cpu_vs_reference`` - CPU agrees with the reference on a shifted grid;
* ``warp_cuda_vs_oracle`` / ``warp_cpu_vs_cuda`` - CUDA reproduces the Fermi-Eyges
                               sigma_x on a shifted grid and CUDA vs CPU agree.

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
from ionmc.physics import fermi_eyges as fe
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.depth_dose import DepthLateralGrid

X0 = materials.WATER.radiation_length_g_per_cm2
SIGMA_TOL = 0.03
SIGMA_BACKEND_TOL_MM = 0.05
REF_WARP_DD_TOL = 5e-4
N_LARGE = 40000
DEPTH_MM = 300.0


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

    model = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, "numpy")
    src = PencilBeamSource(150.0)
    det = dict(straggling=False, scattering=False)
    eng_det = TransportEngine(
        table, WaterSlab(DEPTH_MM), DepthDoseGrid(DEPTH_MM, 10), **det
    )
    eng = TransportEngine(table, WaterSlab(DEPTH_MM), DepthDoseGrid(DEPTH_MM, 10))

    def _det(grid: DepthLateralGrid) -> Any:
        return eng_det.run_scattering(src, grid, 1, seed=4, path="python")

    # -- resolution independence ---------------------------------------------
    coarse = _det(DepthLateralGrid(DEPTH_MM, 150, 25.0, 200))
    fine = _det(DepthLateralGrid(DEPTH_MM, 600, 25.0, 200))
    rebinned = fine.depth_dose_mev.reshape(150, 4).sum(axis=1)
    denom = float(np.sum(coarse.depth_dose_mev))
    res_total = abs(coarse.energy_deposited_mev - fine.energy_deposited_mev) / (
        fine.energy_deposited_mev
    )
    res_bin = float(np.max(np.abs(rebinned - coarse.depth_dose_mev)) / denom)
    report["resolution"] = {"total_rel_diff": res_total, "per_bin_rel_diff": res_bin}
    gates["integral_dose_resolution"] = res_total <= 1e-12 and res_bin <= 1e-12

    # -- alignment independence ----------------------------------------------
    base = _det(DepthLateralGrid(DEPTH_MM, 300, 25.0, 200))
    shifted = _det(
        DepthLateralGrid(
            320.0, 320, 30.0, 240, depth_origin_mm=-20.0, lateral_center_mm=4.0
        )
    )
    align_diff = abs(shifted.energy_deposited_mev - base.energy_deposited_mev) / (
        base.energy_deposited_mev
    )
    report["alignment"] = {"total_rel_diff": align_diff}
    gates["integral_dose_alignment"] = align_diff <= 1e-12

    # -- partial coverage -----------------------------------------------------
    downstream = _det(DepthLateralGrid(200.0, 200, 25.0, 200, depth_origin_mm=100.0))
    report["partial_coverage"] = {
        "full_mev": base.energy_deposited_mev,
        "downstream_mev": downstream.energy_deposited_mev,
    }
    gates["partial_coverage"] = (
        0.0 < downstream.energy_deposited_mev < base.energy_deposited_mev
    )

    # -- lateral origin discriminated: window selection + first moment --------
    # scattering off: the beam stays on the x=0 axis. A window centred on the beam
    # captures the dose; one shifted off it captures none; and on a wide shifted
    # grid the energy-weighted mean lateral position is the beam axis (x=0), not
    # the grid centre. sigma_x alone cannot show this (it is translation-invariant).
    centred_lat = _det(DepthLateralGrid(DEPTH_MM, 300, 3.0, 60))
    off_lat = _det(DepthLateralGrid(DEPTH_MM, 300, 3.0, 60, lateral_center_mm=10.0))
    wide_shift = _det(DepthLateralGrid(DEPTH_MM, 300, 30.0, 600, lateral_center_mm=8.0))
    w = wide_shift.edep_zx_mev.sum(axis=0)
    mean_x = float((w * wide_shift.grid.lateral_centers_mm).sum() / w.sum())
    report["lateral_origin"] = {
        "centred_window_mev": centred_lat.energy_deposited_mev,
        "off_window_mev": off_lat.energy_deposited_mev,
        "mean_x_on_shifted_grid_mm": mean_x,
    }
    gates["lateral_origin"] = (
        centred_lat.energy_deposited_mev > 100.0
        and off_lat.energy_deposited_mev == 0.0
        and abs(mean_x) < 0.2
    )

    # -- lateral-shift sigma_x invariance (reference, scattering on) ----------
    lat_a = DepthLateralGrid(250.0, 500, 30.0, 600)
    lat_b = DepthLateralGrid(250.0, 500, 30.0, 600, lateral_center_mm=5.0)
    ref_a = eng.run_scattering(src, lat_a, 400, seed=7, path="python")
    ref_b = eng.run_scattering(src, lat_b, 400, seed=7, path="python")
    sig_shift = max(
        abs(ref_a.sigma_x_at_depth(z) - ref_b.sigma_x_at_depth(z))
        for z in (90.0, 130.0)
    )
    report["lateral_shift"] = {"max_sigma_x_diff_mm": sig_shift}
    gates["lateral_shift_sigma_x"] = sig_shift <= 1e-6

    # -- Warp cross-backend on a shifted grid --------------------------------
    r0 = float(model.csda_range(150.0)[0]) * 10.0
    shifted_grid = DepthLateralGrid(
        250.0, 500, 30.0, 600, depth_origin_mm=-10.0, lateral_center_mm=4.0
    )

    def _sigma_vs_oracle(res: Any, tag: str, detail: dict[str, Any]) -> bool:
        ok = True
        for frac in (0.5, 0.8):
            z = frac * r0
            mc = res.sigma_x_at_depth(z)
            oracle = float(
                fe.lateral_sigma_x_mm(model, 150.0, np.array([z]), X0, 1.0)[0]
            )
            detail.setdefault(tag, {})[f"{frac}R"] = {"mc": mc, "oracle": oracle}
            ok = ok and abs(mc / oracle - 1.0) <= SIGMA_TOL
        return ok and abs(res.energy_balance) <= 1e-5

    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        gates["warp_cpu_vs_reference"] = False
        if args.require_cuda:
            gates["warp_cuda_vs_oracle"] = False
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
        detail: dict[str, Any] = {}
        profiles: dict[str, Any] = {}
        for device in devices:
            profiles[device] = eng.run_scattering(
                src, shifted_grid, N_LARGE, seed=11, path="warp", device=device
            )
        # CPU vs reference on the shifted grid at a tractable N
        n_match = 4000
        ref_match = eng.run_scattering(
            src, shifted_grid, n_match, seed=13, path="python"
        )
        cpu_match = eng.run_scattering(
            src, shifted_grid, n_match, seed=13, path="warp", device="cpu"
        )
        dd_ref = ref_match.depth_dose_mev
        cum_cpu_ref = float(
            np.max(np.abs(np.cumsum(cpu_match.depth_dose_mev) - np.cumsum(dd_ref)))
            / np.sum(dd_ref)
        )
        sig_cpu_ref = all(
            abs(cpu_match.sigma_x_at_depth(f * r0) - ref_match.sigma_x_at_depth(f * r0))
            <= SIGMA_BACKEND_TOL_MM
            for f in (0.5, 0.8)
        )
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "depth_dose_cumulative": cum_cpu_ref,
            "sigma_x_ok": sig_cpu_ref,
        }
        gates["warp_cpu_vs_reference"] = cum_cpu_ref <= REF_WARP_DD_TOL and sig_cpu_ref

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cuda_vs_oracle"] = False
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            cpu_cuda_ok = True
            cuda_oracle_ok = True
            for device in cuda_devices:
                cuda_oracle_ok = cuda_oracle_ok and _sigma_vs_oracle(
                    profiles[device], device, detail
                )
                cum = float(
                    np.max(
                        np.abs(
                            np.cumsum(profiles[device].depth_dose_mev)
                            - np.cumsum(profiles["cpu"].depth_dose_mev)
                        )
                    )
                    / np.sum(profiles["cpu"].depth_dose_mev)
                )
                report.setdefault("cross_backend", {})[f"{device}_vs_cpu"] = cum
                cpu_cuda_ok = cpu_cuda_ok and cum <= REF_WARP_DD_TOL
            report["sigma_x_vs_oracle"] = detail
            gates["warp_cuda_vs_oracle"] = cuda_oracle_ok
            gates["warp_cpu_vs_cuda"] = cpu_cuda_ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
