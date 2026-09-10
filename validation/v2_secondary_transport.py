"""V2 (closing): secondary charged-particle transport from nonelastic reactions.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached, see ``python -m ionmc.data acquire``). Emits one JSON document
with the decision-0013 gates:

* ``secondary_dose_fraction`` - transported secondary protons contribute ~1-2 %
                                of the local dose at entrance and a few percent of
                                the total dose (Paganetti 2002);
* ``secondary_plateau_shape`` - the secondary-dose fraction rises from ~1-2 % at
                                entrance to ~5-10 % across the plateau, then drops
                                at the sharp Bragg peak (Paganetti 2002);
* ``energy_conservation``     - deposited + escaped = energy_in with the second
                                pass (reference exact, Warp float32 < 1e-5);
* ``secondaries_off_regression`` - secondaries=False reproduces the DEV-007
                                nuclear-only path (no second pass, zero secondaries);
* ``warp_cpu_vs_reference`` / ``warp_cuda_vs_reference`` / ``warp_cpu_vs_cuda``
                                - the backends produce the same secondary set and
                                agree cumulatively (decision 0001 tolerances).

Large-statistics quantities use the Warp CPU/CUDA path; the reference path is run
at small N for the exact budget and cross-backend parity. Exit 0 if every gate
passes, 3 otherwise, 4 if the dataset is missing.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from typing import Any

import numpy as np

from ionmc import __version__, materials
from ionmc.backend import mathlib
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab

ENERGIES = [150.0, 200.0]
N_LARGE = 40000  # Warp path
N_REFERENCE = 1000  # reference path (Python loop) for exact budget and parity
REF_WARP_CUM_TOL = 1e-4
CPU_CUDA_CUM_TOL = 1e-5


def _jsonable(o: Any) -> Any:
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")


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
        json.dump(report, sys.stdout, indent=2, default=_jsonable)
        sys.stdout.write("\n")
        return 4

    slab = WaterSlab(400.0, materials.WATER)
    grid = DepthDoseGrid(400.0, 800)
    eng_on = TransportEngine(table, slab, grid, nuclear=True, secondaries=True)
    eng_off = TransportEngine(table, slab, grid, nuclear=True, secondaries=False)

    # -- reference path: exact budget + off-regression ------------------------
    per_energy: dict[str, Any] = {}
    budget_ok = True
    regression_ok = True
    for e0 in ENERGIES:
        ref = eng_on.run(PencilBeamSource(e0), N_REFERENCE, seed=2024, path="python")
        budget_ok = budget_ok and abs(ref.energy_balance) <= 1e-9
        off = eng_off.run(PencilBeamSource(e0), N_REFERENCE, seed=2024, path="python")
        regression_ok = (
            regression_ok
            and off.n_secondaries == 0
            and off.secondary_edep_mev is None
            and abs(off.energy_balance) <= 1e-9
        )
        per_energy[str(e0)] = {
            "reference_reactions": ref.n_reactions,
            "reference_secondaries": ref.n_secondaries,
            "reference_secondary_dose_mev": float(np.sum(ref.secondary_edep_mev)),
            "reference_energy_balance": ref.energy_balance,
        }
    gates["energy_conservation_reference"] = budget_ok
    gates["secondaries_off_regression"] = regression_ok

    # -- Warp CPU/CUDA: dose fraction, lift, parity ---------------------------
    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        for g in (
            "secondary_dose_fraction",
            "secondary_plateau_shape",
            "warp_cpu_vs_reference",
        ):
            gates[g] = False
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
        fraction_ok = True
        shape_ok = True
        cpu_ref_ok = True
        cuda_ref_ok = True
        profiles: dict[str, dict[float, np.ndarray]] = {d: {} for d in devices}
        for e0 in ENERGIES:
            for device in devices:
                on = eng_on.run(
                    PencilBeamSource(e0), N_LARGE, seed=11, path="warp", device=device
                )
                profiles[device][e0] = on.edep_mev
                frac = on.secondary_dose_fraction
                entrance = float(np.mean(frac[5:25]))
                total_fraction = float(
                    np.sum(on.secondary_edep_mev) / on.energy_deposited_mev
                )
                entry: dict[str, Any] = {
                    "secondaries": on.n_secondaries,
                    "entrance_secondary_fraction": entrance,
                    "total_secondary_fraction": total_fraction,
                    "energy_balance": on.energy_balance,
                }
                report["warp"].setdefault(device, {})[str(e0)] = entry
                if device == "cpu":
                    fraction_ok = (
                        fraction_ok
                        and 0.005 <= entrance <= 0.04
                        and 0.02 <= total_fraction <= 0.12
                    )
                    # secondary-fraction depth shape (Paganetti 2002): the
                    # fraction rises from ~1-2 % at entrance to ~5-10 % across the
                    # plateau proximal to the peak, then collapses at the sharp
                    # Bragg peak where the primary dose dominates. Gate the mean
                    # plateau fraction (entrance..just before the peak) into the
                    # published band, require it above entrance (it rises), and
                    # require it to drop at the peak (secondaries don't peak).
                    peak_bin = int(np.argmax(on.edep_mev))
                    frac_at_peak = float(frac[peak_bin])
                    plateau = frac[25 : max(26, peak_bin - 10)]
                    plateau_mean = float(np.mean(plateau))
                    entry["secondary_fraction_at_peak"] = frac_at_peak
                    entry["plateau_mean_fraction"] = plateau_mean
                    shape_ok = (
                        shape_ok
                        and 0.03 <= plateau_mean <= 0.12
                        and plateau_mean > entrance
                        and frac_at_peak < entrance
                    )
                    # parity vs reference at matched N
                    ref_small = eng_on.run(
                        PencilBeamSource(e0), N_REFERENCE, seed=2024, path="python"
                    )
                    cpu_small = eng_on.run(
                        PencilBeamSource(e0),
                        N_REFERENCE,
                        seed=2024,
                        path="warp",
                        device="cpu",
                    )
                    parity_cum = cumulative_diff(cpu_small.edep_mev, ref_small.edep_mev)
                    entry["parity_secondaries_ref"] = ref_small.n_secondaries
                    entry["parity_secondaries_cpu"] = cpu_small.n_secondaries
                    entry["parity_cumulative_diff"] = parity_cum
                    # the secondary count derives from the float32-vs-float64
                    # reaction energy, so it is compared within a small tolerance;
                    # the cumulative depth dose is the real cross-backend gate
                    secondary_tol = max(1, round(0.005 * ref_small.n_secondaries))
                    cpu_ref_ok = (
                        cpu_ref_ok
                        and abs(on.energy_balance) <= 1e-5
                        and abs(cpu_small.n_secondaries - ref_small.n_secondaries)
                        <= secondary_tol
                        and parity_cum <= REF_WARP_CUM_TOL
                    )
                else:
                    cuda_ref_ok = cuda_ref_ok and abs(on.energy_balance) <= 1e-5
        gates["secondary_dose_fraction"] = fraction_ok
        gates["secondary_plateau_shape"] = shape_ok
        gates["warp_cpu_vs_reference"] = cpu_ref_ok

        cuda_devices = [d for d in devices if d != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cuda_vs_reference"] = False
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            gates["warp_cuda_vs_reference"] = cuda_ref_ok
            cpu_cuda_ok = True
            for e0 in ENERGIES:
                for device in cuda_devices:
                    cum = cumulative_diff(profiles[device][e0], profiles["cpu"][e0])
                    report["warp"].setdefault(f"{device}_vs_cpu", {})[str(e0)] = cum
                    cpu_cuda_ok = cpu_cuda_ok and cum <= CPU_CUDA_CUM_TOL
            gates["warp_cpu_vs_cuda"] = cpu_cuda_ok

    report["per_energy"] = per_energy
    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
