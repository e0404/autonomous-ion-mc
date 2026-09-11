"""V5 (helium): helium-4 transport via equal-velocity z-squared scaling.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0027 gates:

* ``stopping_vs_bethe``   - the scaled-PSTAR helium table matches the independent
                            analytic Bethe model to ~1.5% (E/A = 10-250 MeV/u);
* ``range_identity``      - the helium CSDA range obeys the proton-range identity
                            and R_He(600 MeV) ~ 15.86 g/cm^2;
* ``helium_bragg_range``  - a 600 MeV helium beam stops at ~158 mm, conserving
                            energy (the Bragg-curve target);
* ``warp_cpu_vs_reference`` - the Warp CPU helium depth dose matches the reference;
* ``warp_cpu_vs_cuda``    - the CUDA helium depth dose agrees with CPU.

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
from ionmc.constants import PROTON_MASS_MEV
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table, scale_ion_stopping_table
from ionmc.materials import WATER
from ionmc.particles import ALPHA
from ionmc.stopping_power import AnalyticStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine
from ionmc.transport.geometry import WaterSlab


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
        proton = load_stopping_table(
            cache.load_path(MCSQUARE_PSTAR_WATER, args.cache_dir)
        )
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        json.dump(report, sys.stdout, indent=2, default=_jsonable)
        sys.stdout.write("\n")
        return 4

    helium = scale_ion_stopping_table(proton, ALPHA)
    mr = ALPHA.rest_energy_mev / PROTON_MASS_MEV

    # -- stopping vs independent Bethe model ---------------------------------
    bethe = AnalyticStoppingPower(WATER, ALPHA, path="numpy")
    stop: dict[str, Any] = {}
    ok_stop = True
    for e_per_u in (10.0, 50.0, 100.0, 150.0, 250.0):
        e_he = e_per_u * 4.0
        s_tab = float(np.interp(e_he, helium.energy_mev, helium.stopping_mev_cm2_per_g))
        s_bethe = float(bethe.mass_stopping_power(e_he)[0])
        rel = abs(s_tab - s_bethe) / s_bethe
        stop[f"{e_per_u:.0f}MeV_u"] = {
            "s_table": s_tab,
            "s_bethe": s_bethe,
            "rel": rel,
        }
        ok_stop = ok_stop and rel <= 0.015
    report["stopping"] = stop
    gates["stopping_vs_bethe"] = ok_stop

    # -- range identity + Bragg-range target ---------------------------------
    prefactor = ALPHA.rest_energy_mev / (4.0 * PROTON_MASS_MEV)
    ident_ok = True
    for e_he in (200.0, 400.0, 600.0):
        r_he = float(np.interp(e_he, helium.energy_mev, helium.csda_range_g_per_cm2))
        r_id = prefactor * float(
            np.interp(e_he / mr, proton.energy_mev, proton.csda_range_g_per_cm2)
        )
        ident_ok = ident_ok and abs(r_he - r_id) / r_id <= 1e-6
    r600 = float(np.interp(600.0, helium.energy_mev, helium.csda_range_g_per_cm2))
    report["range"] = {"r_he_600_g_per_cm2": r600, "target": 15.86}
    gates["range_identity"] = bool(ident_ok and abs(r600 - 15.86) <= 0.15)

    # -- helium Bragg range (deterministic transport) ------------------------
    grid = DepthDoseGrid(200.0, 400)
    eng = TransportEngine(
        helium,
        WaterSlab(200.0),
        grid,
        particle=ALPHA,
        straggling=False,
        nuclear=False,
    )
    r = eng.run(PencilBeamSource(600.0), n_histories=1, seed=1, path="python")
    centers = grid.centers_mm
    kpk = int(r.edep_mev.argmax())
    report["bragg"] = {
        "peak_depth_mm": float(centers[kpk]),
        "energy_balance": float(r.energy_balance),
    }
    gates["helium_bragg_range"] = bool(
        155.0 <= centers[kpk] <= 161.0 and abs(r.energy_balance) < 1e-9
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
        gd = DepthDoseGrid(200.0, 200)
        eng_w = TransportEngine(
            helium,
            WaterSlab(200.0),
            gd,
            particle=ALPHA,
            straggling=False,
            nuclear=False,
        )
        src = PencilBeamSource(600.0)
        ref = eng_w.run(src, n_histories=1, seed=2, path="python")
        cpu = eng_w.run(src, n_histories=1, seed=2, path="warp", device="cpu")
        rd, cd = ref.edep_mev, cpu.edep_mev
        cpu_tot = abs(float(cd.sum()) - float(rd.sum())) / float(rd.sum())
        cpu_max = float(np.max(np.abs(cd - rd)) / rd.max())
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "total_rel_diff": cpu_tot,
            "max_bin_rel_diff": cpu_max,
        }
        gates["warp_cpu_vs_reference"] = cpu_tot <= 1e-5 and cpu_max <= 5e-3

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            ok = True
            d0 = float(cd.sum())
            for device in cuda_devices:
                cu = eng_w.run(src, n_histories=1, seed=2, path="warp", device=device)
                tot = abs(float(cu.edep_mev.sum()) - d0) / d0
                mx = float(np.max(np.abs(cu.edep_mev - cd)) / cd.max())
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
