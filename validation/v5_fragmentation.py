"""V5 (fragmentation): the carbon-12 distal dose tail from nuclear fragmentation.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0029 gates:

* ``fragment_tail_present`` - a 290 MeV/u carbon beam deposits 8-20 % of the peak
                              dose just distal to the Bragg peak (canonical ~15 %),
                              and a primary-only run has ~0 dose there;
* ``tail_reach``            - fragment dose stays > 1 % of the peak out to >= 1.5x
                              the carbon range;
* ``primary_survival``      - the surviving primary fraction at the peak matches
                              S(R) = exp(-Sigma*R) ~ 0.47 at 290 MeV/u;
* ``energy_conservation``   - deposited + escaped == E0 exactly, escaped >= 0;
* ``warp_cpu_vs_reference`` - the Warp CPU total depth dose matches the reference;
* ``warp_cpu_vs_cuda``      - the CUDA total depth dose agrees with CPU.

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
from ionmc.data.stopping_tables import load_stopping_table, scale_ion_stopping_table
from ionmc.fragmentation import (
    FRAGMENT_SPECIES,
    carbon_fragmentation_depth_dose,
    macroscopic_carbon_reaction_per_cm,
)
from ionmc.particles import CARBON_12
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine
from ionmc.transport.geometry import WaterSlab

CARBON_RANGE_MM = 163.0  # 290 MeV/u carbon range in water


def _jsonable(o: Any) -> Any:
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")


def _run(proton, slab, grid, src, path, device="cpu"):
    return carbon_fragmentation_depth_dose(
        proton,
        slab,
        grid,
        src,
        n_histories=1,
        seed=1,
        path=path,
        device=device,
        straggling=False,
    )


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
        proton = load_stopping_table(
            cache.load_path(MCSQUARE_PSTAR_WATER, args.cache_dir)
        )
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        json.dump(report, sys.stdout, indent=2, default=_jsonable)
        sys.stdout.write("\n")
        return 4

    slab = WaterSlab(550.0)
    grid = DepthDoseGrid(550.0, 275)  # 2 mm bins
    src = PencilBeamSource(3480.0)
    frag = _run(proton, slab, grid, src, "python")

    # -- fragment tail present (and absent for primary-only) -----------------
    # gate on the resolution-robust integrated distal-dose fraction; tail_to_peak
    # (a single-bin point ratio at the pinned 2 mm resolution) is a diagnostic only.
    distal = {f"{d:.0f}mm": frag.distal_dose_fraction(d) for d in (5.0, 10.0, 20.0)}
    tail = {f"{d:.0f}mm": frag.tail_to_peak(d) for d in (5.0, 10.0, 20.0)}
    carbon = scale_ion_stopping_table(proton, CARBON_12)
    eng = TransportEngine(
        carbon, slab, grid, particle=CARBON_12, straggling=False, nuclear=False
    )
    prim = eng.run(src, n_histories=1, seed=1, path="python").edep_mev
    cprim = grid.centers_mm
    kpp = int(prim.argmax())
    primary_only_distal = float(prim[cprim > cprim[kpp] + 10.0].sum() / prim.sum())
    # resolution invariance of the gate metric across an 8x range of bin widths
    robust = []
    for n_bins in (138, 275, 550, 1100):  # 4, 2, 1, 0.5 mm
        r = _run(proton, slab, DepthDoseGrid(550.0, n_bins), src, "python")
        robust.append(r.distal_dose_fraction(10.0))
    resolution_spread = float(max(robust) - min(robust))
    report["tail"] = {
        "distal_dose_fraction": distal,
        "tail_to_peak_diagnostic_2mm": tail,
        "primary_only_distal_fraction_+10mm": primary_only_distal,
        "resolution_spread_distal_frac_+10mm": resolution_spread,
    }
    gates["fragment_tail_present"] = bool(
        0.10 <= distal["10mm"] <= 0.25
        and distal["5mm"] > distal["10mm"] > distal["20mm"] > 0.10
        and primary_only_distal < 1e-3
        and resolution_spread < 0.01
    )

    # -- tail reach ----------------------------------------------------------
    total = frag.total_edep_mev
    centers = grid.centers_mm
    kpk = int(total.argmax())
    sig = np.where(total > 0.01 * total[kpk])[0]
    reach_mm = float(centers[sig[-1]])
    report["reach"] = {
        "reach_mm": reach_mm,
        "reach_over_range": reach_mm / CARBON_RANGE_MM,
    }
    gates["tail_reach"] = bool(reach_mm >= 1.5 * CARBON_RANGE_MM)

    # -- primary survival ----------------------------------------------------
    sigma = macroscopic_carbon_reaction_per_cm(1.0)
    expected = float(np.exp(-sigma * CARBON_RANGE_MM / 10.0))
    report["survival"] = {
        "at_peak": frag.primary_survival_at_peak,
        "expected": expected,
    }
    gates["primary_survival"] = bool(
        abs(frag.primary_survival_at_peak - expected) / expected < 0.05
        and 0.4 <= frag.primary_survival_at_peak <= 0.55
    )

    # -- energy conservation (independently reconciled, not tautological) -----
    deposited = float(frag.total_edep_mev.sum())
    # independent reconstruction of the injected fragment KE from the reaction
    # weights, multiplicities and residual carbon energy (decision 0029)
    sigma = macroscopic_carbon_reaction_per_cm(1.0)
    s_edge = np.exp(-sigma * grid.edges_mm / 10.0)
    react = s_edge[:-1] - s_edge[1:]
    r_tot = float(np.interp(3480.0, carbon.energy_mev, carbon.csda_range_g_per_cm2))
    resid = r_tot - grid.centers_mm / 10.0
    e_carbon = np.where(
        resid > carbon.csda_range_g_per_cm2[0],
        np.interp(resid, carbon.csda_range_g_per_cm2, carbon.energy_mev),
        0.0,
    )
    n_a_over_12 = sum(m * sp.mass_number for sp, m in FRAGMENT_SPECIES) / 12.0
    injected_independent = float(np.sum(react * e_carbon * n_a_over_12))
    frag_deposited = float(frag.fragment_edep_mev.sum())
    report["energy"] = {
        "energy_in_mev": frag.energy_in_mev,
        "deposited_mev": deposited,
        "escaped_mev": frag.escaped_mev,
        "fragment_fraction": frag_deposited / frag.energy_in_mev,
        "fragment_energy_injected_mev": frag.fragment_energy_injected_mev,
        "fragment_energy_injected_independent_mev": injected_independent,
        "fragment_deposited_mev": frag_deposited,
    }
    gates["energy_conservation"] = bool(
        abs(deposited + frag.escaped_mev - frag.energy_in_mev) / frag.energy_in_mev
        < 1e-9
        and frag.escaped_mev >= 0.0
        # injected KE matches the independent reconstruction (ties the bookkeeping)
        and abs(frag.fragment_energy_injected_mev - injected_independent)
        / injected_independent
        < 1e-6
        # and is deposited in full (fragments stay inside the grid)
        and abs(frag_deposited - frag.fragment_energy_injected_mev)
        / frag.fragment_energy_injected_mev
        < 1e-4
    )

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
        cpu = _run(proton, slab, grid, src, "warp", device="cpu")
        rt, ct = frag.total_edep_mev, cpu.total_edep_mev
        cpu_tot = abs(float(ct.sum()) - float(rt.sum())) / float(rt.sum())
        cpu_max = float(np.max(np.abs(ct - rt)) / rt.max())
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "total_rel_diff": cpu_tot,
            "max_bin_rel_diff": cpu_max,
        }
        gates["warp_cpu_vs_reference"] = bool(cpu_tot <= 5e-5 and cpu_max <= 5e-3)

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            ok = True
            d0 = float(ct.sum())
            for device in cuda_devices:
                cu = _run(proton, slab, grid, src, "warp", device=device)
                tot = abs(float(cu.total_edep_mev.sum()) - d0) / d0
                mx = float(np.max(np.abs(cu.total_edep_mev - ct)) / ct.max())
                cb = report.setdefault("cross_backend", {})
                cb[f"{device}_vs_cpu"] = {"total_rel_diff": tot, "max_bin_rel_diff": mx}
                ok = ok and tot <= 1e-5 and mx <= 5e-3
            gates["warp_cpu_vs_cuda"] = ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
