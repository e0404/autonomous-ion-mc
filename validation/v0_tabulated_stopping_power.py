"""V0 closure: tabulated proton stopping power and range on all execution paths.

Runs on the controlled host runner (numpy + Warp only, no network; the
datasets must already be in the cache, see ``python -m ionmc.data acquire``)
and prints one JSON document with the gates fixed in decision 0008:

* tabulated (PCHIP interpolation of the MCsquare PSTAR water table) vs the
  analytic model of decision 0006 for 10-400 MeV: <= 1.0 percent;
* tabulated CSDA range (from the 0.5 MeV table floor, plus the PSTAR residual
  range below 0.5 MeV) vs NIST PSTAR at 100/150/200/250 MeV: <= 0.1 percent;
* analytic CSDA range (from 1 MeV, plus the PSTAR residual below 1 MeV) vs
  NIST PSTAR at the same energies: <= 0.5 percent;
* numpy path vs Python path: bitwise equal;
* Warp CPU and Warp CUDA (float32) vs the float64 reference: S rtol 1e-5,
  range rtol 2e-5;
* Warp CPU vs Warp CUDA (decision 0001 mixed criterion, class transcendental
  for both lookups): rtol 4e-6, atol 1e-6.

Exit code 0 if every gate passes, 3 otherwise; 4 if a dataset is missing.
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
from ionmc.data import MCSQUARE_G4_WATER, MCSQUARE_PSTAR_WATER, DatasetNotCached
from ionmc.stopping_power import AnalyticStoppingPower
from ionmc.tabulated_stopping_power import TabulatedStoppingPower

# NIST PSTAR, liquid water (material 276), retrieved 2026-09-10; g/cm^2.
PSTAR_CSDA_RANGE = {100.0: 7.718, 150.0: 15.77, 200.0: 25.96, 250.0: 37.94}
PSTAR_CSDA_AT_0_5_MEV = 8.869e-4
PSTAR_CSDA_AT_1_MEV = 2.458e-3

TOL_TABLE_VS_ANALYTIC = 0.010
TOL_TABLE_RANGE = 0.001
TOL_ANALYTIC_RANGE = 0.005
REF_VS_WARP = {"s": (1.0e-5, 0.0), "range": (2.0e-5, 0.0)}
CPU_VS_CUDA = {"s": (4.0e-6, 1.0e-6), "range": (4.0e-6, 1.0e-6)}

GRID = np.linspace(0.5, 400.0, 800)
GRID_ABOVE_10 = np.linspace(10.0, 400.0, 781)
RANGE_ENERGIES = np.array(sorted(PSTAR_CSDA_RANGE))


def normalized_diff(candidate: np.ndarray, ref: np.ndarray, rtol: float, atol: float) -> float:
    return float(np.max(np.abs(candidate - ref) / (atol + rtol * np.abs(ref))))


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
        ref = TabulatedStoppingPower.from_dataset(
            MCSQUARE_PSTAR_WATER, materials.WATER, particles.PROTON, args.cache_dir, path="python"
        )
        npy = TabulatedStoppingPower.from_dataset(
            MCSQUARE_PSTAR_WATER, materials.WATER, particles.PROTON, args.cache_dir, path="numpy"
        )
        g4 = TabulatedStoppingPower.from_dataset(
            MCSQUARE_G4_WATER, materials.WATER, particles.PROTON, args.cache_dir, path="numpy"
        )
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 4
    report["provenance"] = ref.provenance()
    analytic = AnalyticStoppingPower(materials.WATER, particles.PROTON, path="numpy")

    # 1. tabulated vs analytic above 10 MeV
    ratio = npy.mass_stopping_power(GRID_ABOVE_10) / analytic.mass_stopping_power(GRID_ABOVE_10) - 1.0
    report["tabulated_vs_analytic_10_400_mev"] = {
        "max_abs_relative_difference": float(np.max(np.abs(ratio))),
        "mean_relative_difference": float(np.mean(ratio)),
        "tolerance": TOL_TABLE_VS_ANALYTIC,
    }
    gates["tabulated_vs_analytic_above_10_mev"] = bool(np.max(np.abs(ratio)) <= TOL_TABLE_VS_ANALYTIC)

    # 2. ranges vs NIST PSTAR
    pstar = np.array([PSTAR_CSDA_RANGE[e] for e in RANGE_ENERGIES])
    r_tab = npy.csda_range(RANGE_ENERGIES) + PSTAR_CSDA_AT_0_5_MEV
    r_ana = analytic.csda_range(RANGE_ENERGIES) + PSTAR_CSDA_AT_1_MEV
    report["csda_range_vs_pstar"] = {
        "energies_mev": RANGE_ENERGIES.tolist(),
        "pstar_g_per_cm2": pstar.tolist(),
        "tabulated_plus_residual": r_tab.tolist(),
        "tabulated_relative_difference": (r_tab / pstar - 1.0).tolist(),
        "analytic_plus_residual": r_ana.tolist(),
        "analytic_relative_difference": (r_ana / pstar - 1.0).tolist(),
        "residual_below_0_5_mev": PSTAR_CSDA_AT_0_5_MEV,
        "residual_below_1_mev": PSTAR_CSDA_AT_1_MEV,
        "tolerances": {"tabulated": TOL_TABLE_RANGE, "analytic": TOL_ANALYTIC_RANGE},
    }
    gates["tabulated_range_vs_pstar"] = bool(np.all(np.abs(r_tab / pstar - 1.0) <= TOL_TABLE_RANGE))
    gates["analytic_range_vs_pstar"] = bool(np.all(np.abs(r_ana / pstar - 1.0) <= TOL_ANALYTIC_RANGE))

    # 3. numpy vs python
    s_py = ref.mass_stopping_power(GRID)
    s_np = npy.mass_stopping_power(GRID)
    r_py = ref.csda_range(RANGE_ENERGIES)
    gates["numpy_vs_python_exact"] = bool(
        np.array_equal(s_np, s_py) and np.array_equal(npy.csda_range(RANGE_ENERGIES), r_py)
    )

    # 4. second data source (informational)
    e4 = np.array([20.0, 50.0, 100.0, 200.0])
    report["g4_table_vs_pstar_table_percent"] = (
        100.0 * (g4.mass_stopping_power(e4) / npy.mass_stopping_power(e4) - 1.0)
    ).tolist()

    # 5. Warp paths
    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        gates["warp_cpu_vs_reference"] = False
        if args.require_cuda:
            gates["warp_cuda_vs_reference"] = False
            gates["warp_cpu_vs_cuda"] = False
    else:
        wp = mathlib.warp_module()
        wp.config.log_level = wp.LOG_WARNING
        wp.init()
        devices = [d.alias for d in wp.get_devices()]
        report["warp"] = {"available": True, "version": wp.config.version, "devices": devices}
        results: dict[str, dict[str, np.ndarray]] = {}
        for device in devices:
            model = TabulatedStoppingPower(ref.table, materials.WATER, particles.PROTON, path="warp", device=device)
            s_dev = model.mass_stopping_power(GRID)
            r_dev = model.csda_range(RANGE_ENERGIES)
            results[device] = {"s": s_dev, "range": r_dev}
            nd_s = normalized_diff(s_dev, s_py, *REF_VS_WARP["s"])
            nd_r = normalized_diff(r_dev, r_py, *REF_VS_WARP["range"])
            report["warp"][device] = {
                "s_vs_reference": {
                    "max_rel_diff": float(np.max(np.abs(s_dev / s_py - 1.0))),
                    "max_normalized_diff": nd_s,
                    "rtol": REF_VS_WARP["s"][0],
                },
                "range_vs_reference": {
                    "values": r_dev.tolist(),
                    "max_rel_diff": float(np.max(np.abs(r_dev / r_py - 1.0))),
                    "max_normalized_diff": nd_r,
                    "rtol": REF_VS_WARP["range"][0],
                },
            }
            key = "warp_cpu_vs_reference" if device == "cpu" else "warp_cuda_vs_reference"
            gates[key] = bool(nd_s <= 1.0 and nd_r <= 1.0)
        cuda_devices = [d for d in devices if d != "cpu"]
        if not cuda_devices:
            report["warp"]["cuda_note"] = "no CUDA device available"
            if args.require_cuda:
                gates["warp_cuda_vs_reference"] = False
                gates["warp_cpu_vs_cuda"] = False
        for device in cuda_devices:
            nd_s = normalized_diff(results[device]["s"], results["cpu"]["s"], *CPU_VS_CUDA["s"])
            nd_r = normalized_diff(results[device]["range"], results["cpu"]["range"], *CPU_VS_CUDA["range"])
            report["warp"][f"{device}_vs_cpu"] = {
                "class": "transcendental",
                "s_max_abs_diff": float(np.max(np.abs(results[device]["s"] - results["cpu"]["s"]))),
                "s_max_normalized_diff": nd_s,
                "range_max_abs_diff": float(np.max(np.abs(results[device]["range"] - results["cpu"]["range"]))),
                "range_max_normalized_diff": nd_r,
                "rtol": CPU_VS_CUDA["s"][0],
                "atol": CPU_VS_CUDA["s"][1],
            }
            gates["warp_cpu_vs_cuda"] = bool(nd_s <= 1.0 and nd_r <= 1.0)

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
