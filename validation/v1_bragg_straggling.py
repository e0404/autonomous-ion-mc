"""V1 (straggling): proton Bragg peak and range straggling on all paths.

Runs on the controlled host runner (numpy + Warp; the PSTAR water table must be
cached). Emits one JSON document with the decision-0010 gates: unbiased mean
range, range straggling sigma_R vs the analytic Bohr integral and Bortfeld,
energy conservation, a Bragg peak, batch-based statistical consistency of two
seeds, and the reference-vs-Warp and Warp-CPU-vs-CUDA cumulative depth dose.
Exit 0 if every gate passes, 3 otherwise, 4 if the dataset is missing.
"""

from __future__ import annotations

import argparse
import json
import math
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
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab

ENERGIES = [100.0, 150.0, 200.0]


def analytic_sigma_r_mm(
    tab: TabulatedStoppingPower, e0: float, floor_mev: float = 2.0
) -> float:
    # integrate from the MC straggling floor so the reference matches the MC
    e = np.linspace(floor_mev, e0, 4000)
    tau = e / PROTON_MASS_MEV
    gamma = 1.0 + tau
    beta2 = tau * (tau + 2.0) / (gamma * gamma)
    f_rel = (1.0 - 0.5 * beta2) / (1.0 - beta2)
    za = materials.WATER.electrons_per_gram_ratio
    d_omega2_dx = tphys.BOHR_K_MEV2_CM2_PER_MOL * za * f_rel
    s_lin = tab.mass_stopping_power(e)
    return math.sqrt(np.trapezoid(d_omega2_dx / s_lin**3, e)) * 10.0


def cumulative_diff(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(np.cumsum(a) - np.cumsum(b))) / float(np.sum(b)))


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

    tab = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, path="numpy")
    tab_py = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, path="python")
    slab = WaterSlab(400.0, materials.WATER)
    grid = DepthDoseGrid(400.0, 2000)
    engine = TransportEngine(table, slab, grid, straggling=True)

    per_energy: dict[str, Any] = {}
    mean_ok = sigma_ana_ok = sigma_bort_ok = frac_ok = peak_ok = True
    for e0 in ENERGIES:
        res = engine.run(PencilBeamSource(e0), n_histories=args.histories, seed=7, path="warp")
        r_csda = float(tab_py.csda_range(e0)[0]) * 10.0
        sig_ana = analytic_sigma_r_mm(tab, e0)
        sig_bort = 0.012 * (r_csda / 10.0) ** 0.935 * 10.0
        frac = res.range_sigma_mm / res.range_mean_mm
        centers = grid.centers_mm
        peak_bin = int(np.argmax(res.edep_mev))
        entrance = float(res.edep_mev[centers < 20.0].mean())
        per_energy[str(e0)] = {
            "mean_range_mm": res.range_mean_mm,
            "csda_range_mm": r_csda,
            "mean_over_csda_minus_1": res.range_mean_mm / r_csda - 1.0,
            "sigma_r_mm": res.range_sigma_mm,
            "sigma_analytic_mm": sig_ana,
            "sigma_bortfeld_mm": sig_bort,
            "sigma_over_analytic": res.range_sigma_mm / sig_ana,
            "sigma_over_bortfeld": res.range_sigma_mm / sig_bort,
            "sigma_over_range_percent": 100.0 * frac,
            "peak_depth_mm": float(centers[peak_bin]),
            "peak_over_entrance": float(res.edep_mev[peak_bin] / entrance),
            "energy_balance": res.energy_balance,
        }
        mean_ok = mean_ok and bool(abs(res.range_mean_mm / r_csda - 1.0) <= 1e-3)
        sigma_ana_ok = sigma_ana_ok and bool(abs(res.range_sigma_mm / sig_ana - 1.0) <= 0.03)
        sigma_bort_ok = sigma_bort_ok and bool(abs(res.range_sigma_mm / sig_bort - 1.0) <= 0.10)
        frac_ok = frac_ok and bool(0.009 <= frac <= 0.012)
        peak_ok = peak_ok and bool(
            centers[peak_bin] > 0.9 * r_csda
            and res.edep_mev[peak_bin] > 3.0 * entrance
        )
    report["per_energy"] = per_energy
    gates["mean_range_unbiased"] = bool(mean_ok)
    gates["sigma_vs_analytic"] = bool(sigma_ana_ok)
    gates["sigma_vs_bortfeld"] = bool(sigma_bort_ok)
    gates["sigma_over_range_band"] = bool(frac_ok)
    gates["bragg_peak_present"] = bool(peak_ok)

    # statistical consistency of two independent seeds at 150 MeV
    a = engine.run_batched(PencilBeamSource(150.0), args.histories, 10, seed=100, path="warp")
    b = engine.run_batched(PencilBeamSource(150.0), args.histories, 10, seed=500, path="warp")
    mask = a.mean_edep_mev > 0.01 * a.mean_edep_mev.max()
    combined = np.sqrt(a.standard_error_mev**2 + b.standard_error_mev**2)
    t = (a.mean_edep_mev[mask] - b.mean_edep_mev[mask]) / combined[mask]
    rms_t = float(math.sqrt(np.mean(t**2)))
    report["statistical_consistency"] = {"rms_t": rms_t, "max_abs_t": float(np.max(np.abs(t)))}
    gates["statistical_consistency"] = bool(0.5 <= rms_t <= 1.6 and np.max(np.abs(t)) <= 5.0)

    # cross-backend cumulative (shared streams), reference vs Warp and CPU vs CUDA
    if not mathlib.HAVE_WARP:
        gates["warp_cpu_vs_reference"] = False
        if args.require_cuda:
            gates["warp_cpu_vs_cuda"] = False
    else:
        wp = mathlib.warp_module()
        devices = [d.alias for d in wp.get_devices()]
        report["warp_devices"] = devices
        ref = engine.run(PencilBeamSource(150.0), 300, seed=11, path="python")
        cpu = engine.run(PencilBeamSource(150.0), 300, seed=11, path="warp", device="cpu")
        cum_cpu_ref = cumulative_diff(cpu.edep_mev, ref.edep_mev)
        report["cumulative_cpu_vs_reference"] = cum_cpu_ref
        gates["warp_cpu_vs_reference"] = bool(cum_cpu_ref <= 1e-4)
        cuda_devices = [d for d in devices if d != "cpu"]
        if cuda_devices:
            big_cpu = engine.run(PencilBeamSource(150.0), args.histories, seed=11, path="warp", device="cpu")
            big_cuda = engine.run(PencilBeamSource(150.0), args.histories, seed=11, path="warp", device="cuda:0")
            cum = cumulative_diff(big_cuda.edep_mev, big_cpu.edep_mev)
            report["cumulative_cpu_vs_cuda"] = cum
            report["sigma_r_cpu_vs_cuda_mm"] = [big_cpu.range_sigma_mm, big_cuda.range_sigma_mm]
            gates["warp_cpu_vs_cuda"] = bool(cum <= 1e-5)
        elif args.require_cuda:
            gates["warp_cpu_vs_cuda"] = False

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
