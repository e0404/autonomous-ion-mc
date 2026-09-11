"""V5 (carbon/oxygen): carbon-12 and oxygen-16 primary transport via z^2 scaling.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0028 gates:

* ``stopping_vs_bethe``   - the scaled carbon/oxygen table matches the independent
                            Bethe model within the ion-specific tolerance (carbon
                            <= 2 %, oxygen <= 3 %); the Z-dependent error is
                            recorded;
* ``range_identity``      - the carbon/oxygen CSDA range obeys the proton-range
                            identity;
* ``carbon_bragg_range``  - a 290 MeV/u carbon beam stops at ~163 mm, conserving
                            energy (the primary Bragg curve; fragment tail
                            deferred);
* ``warp_cpu_vs_reference`` - the Warp CPU carbon depth dose matches the reference;
* ``warp_cpu_vs_cuda``    - the CUDA carbon depth dose agrees with CPU.

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
from ionmc.particles import CARBON_12, OXYGEN_16
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

    # -- stopping vs independent Bethe (Z-dependent tolerance) ----------------
    stop: dict[str, Any] = {}
    ok_stop = True
    for particle, tol in ((CARBON_12, 0.02), (OXYGEN_16, 0.03)):
        tab = scale_ion_stopping_table(proton, particle)
        bethe = AnalyticStoppingPower(WATER, particle, path="numpy")
        worst = 0.0
        for e_per_u in (10.0, 50.0, 100.0, 150.0, 250.0, 400.0):
            e = e_per_u * particle.mass_number
            s_tab = float(np.interp(e, tab.energy_mev, tab.stopping_mev_cm2_per_g))
            s_bethe = float(bethe.mass_stopping_power(e)[0])
            worst = max(worst, abs(s_tab - s_bethe) / s_bethe)
        stop[particle.name] = {"worst_rel": worst, "tol": tol}
        ok_stop = ok_stop and worst <= tol
    report["stopping"] = stop
    gates["stopping_vs_bethe"] = ok_stop

    # -- range identity (carbon + oxygen) ------------------------------------
    ident_ok = True
    for particle in (CARBON_12, OXYGEN_16):
        tab = scale_ion_stopping_table(proton, particle)
        mr = particle.rest_energy_mev / PROTON_MASS_MEV
        pref = particle.rest_energy_mev / (particle.charge**2 * PROTON_MASS_MEV)
        for e_per_u in (100.0, 290.0):
            e = e_per_u * particle.mass_number
            r_ion = float(np.interp(e, tab.energy_mev, tab.csda_range_g_per_cm2))
            r_id = pref * float(
                np.interp(e / mr, proton.energy_mev, proton.csda_range_g_per_cm2)
            )
            ident_ok = ident_ok and abs(r_ion - r_id) / r_id <= 1e-6
    gates["range_identity"] = ident_ok

    # -- carbon Bragg range (deterministic) ----------------------------------
    carbon = scale_ion_stopping_table(proton, CARBON_12)
    grid = DepthDoseGrid(250.0, 500)
    eng = TransportEngine(
        carbon,
        WaterSlab(250.0),
        grid,
        particle=CARBON_12,
        straggling=False,
        nuclear=False,
    )
    r = eng.run(PencilBeamSource(3480.0), n_histories=1, seed=1, path="python")
    kpk = int(r.edep_mev.argmax())
    report["bragg"] = {
        "peak_depth_mm": float(grid.centers_mm[kpk]),
        "energy_balance": float(r.energy_balance),
    }
    gates["carbon_bragg_range"] = bool(
        159.0 <= grid.centers_mm[kpk] <= 167.0 and abs(r.energy_balance) < 1e-9
    )

    # -- Warp cross-backend (carbon) -----------------------------------------
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
        gd = DepthDoseGrid(250.0, 250)
        eng_w = TransportEngine(
            carbon,
            WaterSlab(250.0),
            gd,
            particle=CARBON_12,
            straggling=False,
            nuclear=False,
        )
        src = PencilBeamSource(3480.0)
        ref = eng_w.run(src, n_histories=1, seed=2, path="python")
        cpu = eng_w.run(src, n_histories=1, seed=2, path="warp", device="cpu")
        rd, cd = ref.edep_mev, cpu.edep_mev
        # heavier-ion float32 budget (~500 steps at 3480 MeV): total <= 5e-5
        cpu_tot = abs(float(cd.sum()) - float(rd.sum())) / float(rd.sum())
        cpu_max = float(np.max(np.abs(cd - rd)) / rd.max())
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "total_rel_diff": cpu_tot,
            "max_bin_rel_diff": cpu_max,
        }
        gates["warp_cpu_vs_reference"] = cpu_tot <= 5e-5 and cpu_max <= 5e-3

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
                ok = ok and tot <= 5e-5 and mx <= 5e-3
            gates["warp_cpu_vs_cuda"] = ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
