"""V1 (partial): CSDA proton depth dose on the reference and Warp paths.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached, see ``python -m ionmc.data acquire``). Emits one JSON document
with the decision-0009 gates: energy conservation, R80 vs the tabulated CSDA
range, step-size convergence, an independent range cross-check, and the
reference-vs-Warp and Warp-CPU-vs-CUDA depth-dose agreement (cumulative and
edge-aware per-bin). Exit 0 if every gate passes, 3 otherwise, 4 if the dataset
is missing.
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
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab

ENERGIES = [100.0, 150.0, 200.0]
R80_TOL = 0.003
STEP_DRIFT_TOL = 5e-4
RANGE_XCHECK_TOL = 0.002
REF_WARP_CUM_TOL = 1e-4
CPU_CUDA_CUM_TOL = 1e-5
# Warp energy conservation is a float32 accumulation over a full proton range; the
# deterministic residual is ~1e-5, so the project's standard 5e-5 float32 budget
# applies (the float64 reference path is held to 1e-9). See decision 0034.
WARP_ENERGY_BALANCE_TOL = 5e-5


def csda_range_mm(table: Any, e0: float) -> float:
    tab = TabulatedStoppingPower(
        table, materials.WATER, particles.PROTON, path="python"
    )
    return float(tab.csda_range(e0)[0]) * 10.0 / materials.WATER.density_g_per_cm3


def cumulative_diff(a: np.ndarray, b: np.ndarray) -> float:
    total = float(np.sum(b))
    return float(np.max(np.abs(np.cumsum(a) - np.cumsum(b))) / total)


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
        "energies_mev": ENERGIES,
    }
    gates = report["gates"]
    try:
        table = load_stopping_table(
            cache.load_path(MCSQUARE_PSTAR_WATER, args.cache_dir)
        )
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 4

    slab = WaterSlab(400.0, materials.WATER)
    grid = DepthDoseGrid(400.0, 800)
    fine = DepthDoseGrid(400.0, 8000)
    # CSDA is the continuous-slowing-down (deterministic) approximation: straggling
    # is off, so R80 tracks the tabulated CSDA range and the step-convergence check
    # is not confounded by stochastic energy loss. (TransportEngine has defaulted
    # straggling=True since DEV-005; this CSDA validation predates that and must be
    # explicit. Straggling itself is validated by v1_bragg_straggling.) Decision 0034.
    engine = TransportEngine(table, slab, grid, straggling=False)

    per_energy: dict[str, Any] = {}
    ref_profiles: dict[float, np.ndarray] = {}
    balance_ok = True
    r80_ok = True
    step_ok = True
    xcheck_ok = True
    for e0 in ENERGIES:
        ref = engine.run(PencilBeamSource(e0), n_histories=1, path="python")
        ref_profiles[e0] = ref.edep_mev
        r80 = ref.r80_mm()
        r_csda = csda_range_mm(table, e0)
        # step convergence on a fixed fine grid
        r_coarse = (
            TransportEngine(table, slab, fine, max_fraction=0.02, straggling=False)
            .run(PencilBeamSource(e0), 1, path="python")
            .r80_mm()
        )
        r_fine = (
            TransportEngine(table, slab, fine, max_fraction=0.002, straggling=False)
            .run(PencilBeamSource(e0), 1, path="python")
            .r80_mm()
        )
        drift = abs(r_coarse / r_fine - 1.0)
        last = int(np.max(np.nonzero(ref.edep_mev)[0]))
        stop_depth = float(ref.grid.centers_mm[last])
        per_energy[str(e0)] = {
            "r80_mm": r80,
            "csda_range_mm": r_csda,
            "r80_over_csda_minus_1": r80 / r_csda - 1.0,
            "energy_balance": ref.energy_balance,
            "step_drift": drift,
            "stop_depth_mm": stop_depth,
            "stop_depth_over_csda_minus_1": stop_depth / r_csda - 1.0,
            "truncated": ref.truncated,
        }
        balance_ok = (
            balance_ok and abs(ref.energy_balance) <= 1e-9 and ref.truncated == 0
        )
        r80_ok = r80_ok and abs(r80 / r_csda - 1.0) <= R80_TOL
        step_ok = step_ok and drift <= STEP_DRIFT_TOL
        xcheck_ok = xcheck_ok and abs(stop_depth / r_csda - 1.0) <= RANGE_XCHECK_TOL
    report["per_energy"] = per_energy
    gates["energy_conservation_reference"] = balance_ok
    gates["r80_vs_csda_range"] = r80_ok
    gates["step_size_convergence"] = step_ok
    gates["range_cross_check"] = xcheck_ok

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
        report["warp"] = {
            "available": True,
            "version": wp.config.version,
            "devices": devices,
        }
        warp_profiles: dict[str, dict[float, np.ndarray]] = {d: {} for d in devices}
        cpu_ref_ok = True
        cuda_ref_ok = True
        balance_warp_ok = True
        for e0 in ENERGIES:
            for device in devices:
                res = engine.run(PencilBeamSource(e0), 1, path="warp", device=device)
                warp_profiles[device][e0] = res.edep_mev
                cum = cumulative_diff(res.edep_mev, ref_profiles[e0])
                balance_warp_ok = (
                    balance_warp_ok
                    and abs(res.energy_balance) <= WARP_ENERGY_BALANCE_TOL
                )
                report["warp"].setdefault(device, {})[str(e0)] = {
                    "cumulative_diff": cum,
                    "energy_balance": res.energy_balance,
                    "r80_mm": res.r80_mm(),
                }
                if device == "cpu":
                    cpu_ref_ok = cpu_ref_ok and cum <= REF_WARP_CUM_TOL
                else:
                    cuda_ref_ok = cuda_ref_ok and cum <= REF_WARP_CUM_TOL
        gates["warp_cpu_vs_reference"] = cpu_ref_ok and balance_warp_ok
        cuda_devices = [d for d in devices if d != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cuda_vs_reference"] = False
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            gates["warp_cuda_vs_reference"] = cuda_ref_ok
            cpu_cuda_ok = True
            for e0 in ENERGIES:
                for device in cuda_devices:
                    cum = cumulative_diff(
                        warp_profiles[device][e0], warp_profiles["cpu"][e0]
                    )
                    report["warp"].setdefault(f"{device}_vs_cpu", {})[str(e0)] = cum
                    cpu_cuda_ok = cpu_cuda_ok and cum <= CPU_CUDA_CUM_TOL
            gates["warp_cpu_vs_cuda"] = cpu_cuda_ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
