"""V1 (scattering): proton lateral spread from multiple Coulomb scattering.

Runs on the controlled host runner (numpy + Warp; the PSTAR water table must be
cached). Emits one JSON document with the decision-0011 gates: theta0/pv hand
values, MC sigma_x(z) vs the Fermi-Eyges oracle and vs published values,
energy conservation, the depth-dose Bragg peak and detour factor, and the
Warp-CPU-vs-CUDA lateral spread. Exit 0 if every gate passes, 3 otherwise, 4 if
the dataset is missing.
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
from ionmc.constants import PROTON_MASS_MEV
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import transport as tphys
from ionmc.physics.fermi_eyges import lateral_sigma_x_mm
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import (
    DepthDoseGrid,
    DepthLateralGrid,
    PencilBeamSource,
    TransportEngine,
    WaterSlab,
)

ENERGIES = [150.0, 200.0]
PUBLISHED_SIGMA_08R = {150.0: 2.4, 200.0: 3.9}
X0 = 36.08


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-cuda", action="store_true")
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--histories", type=int, default=40000)
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
        table = load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, args.cache_dir))
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 4

    if mathlib.HAVE_WARP:
        _wp = mathlib.warp_module()
        _wp.config.log_level = _wp.LOG_WARNING
        _wp.init()

    model = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, path="numpy")
    engine = TransportEngine(table, WaterSlab(400.0, materials.WATER), DepthDoseGrid(400.0, 10))

    # theta0 / pv hand values
    pv_ok = abs(tphys.momentum_times_velocity(100.0, PROTON_MASS_MEV) / 190.369 - 1.0) <= 1e-3
    th_ok = abs(tphys.highland_theta0(100.0, PROTON_MASS_MEV, 1.0, 1.0, 1.0, X0) * 1000.0 / 3.761 - 1.0) <= 1e-3
    gates["theta0_pv_hand_values"] = bool(pv_ok and th_ok)

    if not mathlib.HAVE_WARP:
        gates["sigma_vs_fermi_eyges"] = False
        gates["sigma_vs_published"] = False
        report["all_gates_passed"] = all(gates.values())
        json.dump(report, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 3

    per_energy: dict[str, Any] = {}
    fe_ok = pub_ok = balance_ok = peak_ok = detour_ok = True
    profiles: dict[float, Any] = {}
    for e0 in ENERGIES:
        r = float(model.csda_range(e0)[0]) * 10.0
        grid = DepthLateralGrid(r * 1.05, int(r * 1.05 / 0.5), 25.0, 500)
        res = engine.run_scattering(PencilBeamSource(e0), grid, args.histories, seed=42, path="warp")
        profiles[e0] = (grid, res)
        depths = np.array([0.5 * r, 0.8 * r])
        fe = lateral_sigma_x_mm(model, e0, depths, X0)
        mc = np.array([res.sigma_x_at_depth(d) for d in depths])
        dd = res.depth_dose_mev
        centers = grid.depth_centers_mm
        peak_bin = int(np.argmax(dd))
        entrance = float(dd[centers < 20.0].mean())
        per_energy[str(e0)] = {
            "range_mm": r,
            "sigma_mc_05R": float(mc[0]),
            "sigma_mc_08R": float(mc[1]),
            "sigma_fe_05R": float(fe[0]),
            "sigma_fe_08R": float(fe[1]),
            "mc_over_fe_05R": float(mc[0] / fe[0]),
            "mc_over_fe_08R": float(mc[1] / fe[1]),
            "sigma_published_08R": PUBLISHED_SIGMA_08R[e0],
            "mc_over_published_08R": float(mc[1] / PUBLISHED_SIGMA_08R[e0]),
            "peak_depth_mm": float(centers[peak_bin]),
            "peak_over_entrance": float(dd[peak_bin] / entrance),
            "energy_balance": res.energy_balance,
        }
        fe_ok = fe_ok and bool(np.all(np.abs(mc / fe - 1.0) <= 0.03))
        pub_ok = pub_ok and bool(abs(mc[1] / PUBLISHED_SIGMA_08R[e0] - 1.0) <= 0.08)
        balance_ok = balance_ok and bool(abs(res.energy_balance) <= 5e-5)
        peak_ok = peak_ok and bool(centers[peak_bin] > 0.9 * r and dd[peak_bin] > 3.0 * entrance)
        detour_ok = detour_ok and bool(centers[peak_bin] <= r * 1.001)
    report["per_energy"] = per_energy
    gates["sigma_vs_fermi_eyges"] = fe_ok
    gates["sigma_vs_published"] = pub_ok
    gates["energy_conservation"] = balance_ok
    gates["bragg_peak_marginal"] = peak_ok
    gates["detour_factor_small"] = detour_ok

    devices = [d.alias for d in mathlib.warp_module().get_devices()]
    report["warp_devices"] = devices
    cuda_devices = [d for d in devices if d != "cpu"]
    if cuda_devices:
        r = float(model.csda_range(150.0)[0]) * 10.0
        grid = DepthLateralGrid(r * 1.05, int(r * 1.05 / 0.5), 25.0, 500)
        cpu = engine.run_scattering(PencilBeamSource(150.0), grid, args.histories, seed=9, path="warp", device="cpu")
        cuda = engine.run_scattering(PencilBeamSource(150.0), grid, args.histories, seed=9, path="warp", device="cuda:0")
        d08 = 0.8 * r
        diff = abs(cpu.sigma_x_at_depth(d08) - cuda.sigma_x_at_depth(d08))
        report["cpu_vs_cuda_sigma_08R_mm"] = [cpu.sigma_x_at_depth(d08), cuda.sigma_x_at_depth(d08)]
        gates["warp_cpu_vs_cuda"] = bool(diff <= 0.02)
    elif args.require_cuda:
        gates["warp_cpu_vs_cuda"] = False

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
