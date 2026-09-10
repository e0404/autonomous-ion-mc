"""V0 validation: analytical proton stopping power on all execution paths.

Runs on the controlled host runner (numpy + Warp only) and prints one JSON
document with:

* analytic (reference Python, float64) vs PSTAR-derived water table
  (decision 0006 criterion: <= 1.0 percent for E >= 10 MeV);
* CSDA ranges at 100/150/200/250 MeV (informational vs recalled PSTAR values,
  0.5 percent guard) and the R(200)-R(100) difference;
* numpy path vs Python path (expected bitwise equal);
* Warp CPU and Warp CUDA (float32) vs the float64 reference
  (decision 0005: rtol 1e-5 for S, 2e-5 for ranges);
* Warp CPU vs Warp CUDA under the decision 0001 mixed criterion with the
  per-kernel classes fixed in decision 0006 (S: rtol 4e-6/atol 1e-6;
  range: rtol 1e-5/atol 1e-6);
* the ICRU 90 (78 eV) offset as a separately reported model difference.

Exit code 0 if every gate passes, 3 otherwise. Producing the report never
fails by itself; missing CUDA is reported and fails the gate only when
``--require-cuda`` is given.
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
from ionmc.reference_data.pstar_water import PROVENANCE, PSTAR_WATER_STOPPING_POWER
from ionmc.stopping_power import AnalyticStoppingPower

TOL_PSTAR = 0.010
TOL_RANGE_RECALLED = 0.005
RECALLED_PSTAR_CSDA = {100.0: 7.718, 150.0: 15.77, 200.0: 25.96, 250.0: 37.94}
REF_VS_WARP = {"s": (1.0e-5, 0.0), "range": (2.0e-5, 0.0)}
CPU_VS_CUDA = {"s": (4.0e-6, 1.0e-6), "range": (1.0e-5, 1.0e-6)}

ENERGIES = np.array(sorted(PSTAR_WATER_STOPPING_POWER))
RANGE_ENERGIES = np.array(sorted(RECALLED_PSTAR_CSDA))
GRID = np.linspace(2.0, 400.0, 400)


def normalized_diff(
    candidate: np.ndarray, ref: np.ndarray, rtol: float, atol: float
) -> float:
    return float(np.max(np.abs(candidate - ref) / (atol + rtol * np.abs(ref))))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-cuda", action="store_true")
    args = parser.parse_args()

    report: dict[str, Any] = {
        "schema_version": 1,
        "ionmc_version": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "warp": None,
        "reference_table": PROVENANCE,
        "gates": {},
    }
    gates = report["gates"]

    ref = AnalyticStoppingPower(materials.WATER, particles.PROTON, path="python")
    npy = AnalyticStoppingPower(materials.WATER, particles.PROTON, path="numpy")
    report["provenance"] = ref.provenance()

    # 1. analytic vs PSTAR
    s_ref = ref.mass_stopping_power(ENERGIES)
    table = np.array([PSTAR_WATER_STOPPING_POWER[e] for e in ENERGIES])
    rel = s_ref / table - 1.0
    above = ENERGIES >= 10.0
    report["analytic_vs_pstar"] = {
        "energies_mev": ENERGIES.tolist(),
        "analytic": s_ref.tolist(),
        "pstar": table.tolist(),
        "relative_difference": rel.tolist(),
        "max_abs_relative_difference_above_10_mev": float(np.max(np.abs(rel[above]))),
        "tolerance": TOL_PSTAR,
    }
    gates["analytic_vs_pstar_above_10_mev"] = bool(
        np.all(np.abs(rel[above]) <= TOL_PSTAR)
    )

    # 2. ranges
    r_ref = ref.csda_range(RANGE_ENERGIES)
    recalled = np.array([RECALLED_PSTAR_CSDA[e] for e in RANGE_ENERGIES])
    report["csda_range"] = {
        "energies_mev": RANGE_ENERGIES.tolist(),
        "analytic_g_per_cm2": r_ref.tolist(),
        "recalled_pstar_g_per_cm2": recalled.tolist(),
        "relative_difference": (r_ref / recalled - 1.0).tolist(),
        "r200_minus_r100": float(r_ref[2] - r_ref[0]),
        "note": "recalled PSTAR values are informational (not fetched in-session)",
    }
    gates["csda_range_recalled_guard"] = bool(
        np.all(np.abs(r_ref / recalled - 1.0) <= TOL_RANGE_RECALLED)
    )

    # 3. numpy vs python
    s_np = npy.mass_stopping_power(GRID)
    s_py = ref.mass_stopping_power(GRID)
    report["numpy_vs_python"] = {"max_abs_diff": float(np.max(np.abs(s_np - s_py)))}
    gates["numpy_vs_python_exact"] = bool(np.array_equal(s_np, s_py))

    # 4. ICRU 90 offset (informational)
    s90 = AnalyticStoppingPower(materials.WATER_ICRU90, particles.PROTON, path="numpy")
    ratio90 = s90.mass_stopping_power(ENERGIES) / npy.mass_stopping_power(ENERGIES)
    report["icru90_offset_percent"] = (100.0 * (ratio90 - 1.0)).tolist()

    # 5. Warp paths
    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        gates["warp_cpu_vs_reference"] = False
        if args.require_cuda:
            gates["warp_cuda_vs_reference"] = False
            gates["warp_cpu_vs_cuda"] = False
    else:
        wp = mathlib.warp_module()
        # Keep Warp's banner and module-load notices off stdout: this script's
        # stdout is one JSON document (Warp >= 1.17 API; ``config.quiet`` is
        # deprecated there).
        wp.config.log_level = wp.LOG_WARNING
        wp.init()
        devices = [d.alias for d in wp.get_devices()]
        report["warp"] = {
            "available": True,
            "version": wp.config.version,
            "devices": devices,
        }
        r_py_grid = ref.csda_range(RANGE_ENERGIES)
        results: dict[str, dict[str, np.ndarray]] = {}
        for device in devices:
            model = AnalyticStoppingPower(
                materials.WATER, particles.PROTON, path="warp", device=device
            )
            s_dev = model.mass_stopping_power(GRID)
            r_dev = model.csda_range(RANGE_ENERGIES)
            results[device] = {"s": s_dev, "range": r_dev}
            nd_s = normalized_diff(s_dev, s_py, *REF_VS_WARP["s"])
            nd_r = normalized_diff(r_dev, r_py_grid, *REF_VS_WARP["range"])
            report["warp"][device] = {
                "s_vs_reference": {
                    "max_rel_diff": float(np.max(np.abs(s_dev / s_py - 1.0))),
                    "max_normalized_diff": nd_s,
                    "rtol": REF_VS_WARP["s"][0],
                },
                "range_vs_reference": {
                    "values": r_dev.tolist(),
                    "max_rel_diff": float(np.max(np.abs(r_dev / r_py_grid - 1.0))),
                    "max_normalized_diff": nd_r,
                    "rtol": REF_VS_WARP["range"][0],
                },
            }
            key = (
                "warp_cpu_vs_reference"
                if device == "cpu"
                else "warp_cuda_vs_reference"
            )
            gates[key] = bool(nd_s <= 1.0 and nd_r <= 1.0)
        cuda_devices = [d for d in devices if d != "cpu"]
        if not cuda_devices:
            report["warp"]["cuda_note"] = "no CUDA device available"
            if args.require_cuda:
                gates["warp_cuda_vs_reference"] = False
                gates["warp_cpu_vs_cuda"] = False
        for device in cuda_devices:
            s_dev, s_cpu = results[device]["s"], results["cpu"]["s"]
            nd_s = normalized_diff(s_dev, s_cpu, *CPU_VS_CUDA["s"])
            r_dev, r_cpu = results[device]["range"], results["cpu"]["range"]
            nd_r = normalized_diff(r_dev, r_cpu, *CPU_VS_CUDA["range"])
            report["warp"][f"{device}_vs_cpu"] = {
                "s": {
                    "class": "transcendental",
                    "max_abs_diff": float(np.max(np.abs(s_dev - s_cpu))),
                    "max_normalized_diff": nd_s,
                    "rtol": CPU_VS_CUDA["s"][0],
                    "atol": CPU_VS_CUDA["s"][1],
                },
                "range": {
                    "class": "iterative_accumulation",
                    "max_abs_diff": float(np.max(np.abs(r_dev - r_cpu))),
                    "max_normalized_diff": nd_r,
                    "rtol": CPU_VS_CUDA["range"][0],
                    "atol": CPU_VS_CUDA["range"][1],
                },
            }
            gates["warp_cpu_vs_cuda"] = bool(nd_s <= 1.0 and nd_r <= 1.0)

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
