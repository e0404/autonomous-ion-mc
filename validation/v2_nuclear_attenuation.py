"""V2 (opening): proton nonelastic nuclear attenuation and local deposition.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached, see ``python -m ionmc.data acquire``). Emits one JSON document
with the decision-0012 gates:

* ``primary_survival``      - primary survival to the Bragg peak matches the
                              published 0.80 +- 0.03 (150 MeV) / 0.73 +- 0.04
                              (200 MeV) (Paganetti 2002; Gottschalk);
* ``energy_conservation``   - deposited + escaped = energy_in (reference exact,
                              Warp float32 < 1e-5);
* ``reaction_self_consistency`` - the MC reaction fraction agrees with the
                              analytic 1 - exp(-integral Sigma/S dE);
* ``nuclear_off_regression`` - nuclear=False reproduces the EM baseline exactly
                              (zero reactions, bit-identical depth dose);
* ``warp_cpu_vs_reference`` / ``warp_cuda_vs_reference`` / ``warp_cpu_vs_cuda``
                              - the backends remove the identical primary set and
                              agree cumulatively (decision 0001 tolerances).

Large-statistics survival uses the Warp CPU path (fast); the reference path is
run at small N for the exact energy budget and cross-backend parity. Exit 0 if
every gate passes, 3 otherwise, 4 if the dataset is missing.
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
from ionmc.constants import AVOGADRO
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import nuclear as nuclear_ref
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab

ENERGIES = [150.0, 200.0]
PUBLISHED_SURVIVAL = {150.0: (0.80, 0.03), 200.0: (0.73, 0.04)}
N_SURVIVAL = 40000  # Warp CPU, statistical error on the fraction ~0.002
N_REFERENCE = 2000  # reference path (Python loop) for exact budget and parity
REF_WARP_CUM_TOL = 1e-4
CPU_CUDA_CUM_TOL = 1e-5


def _jsonable(o: Any) -> Any:
    """Coerce numpy scalars (bool_/int/float) to native Python for json.dump."""
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


def analytic_reaction_fraction(table: Any, e0: float, oxygen_density: float) -> float:
    """1 - exp(-integral Sigma/S_lin dE) along the CSDA track (cutoff..E0)."""
    density = materials.WATER.density_g_per_cm3
    model = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, "numpy")
    e_grid = np.linspace(1.0, e0, 8000)
    s_lin = model.mass_stopping_power(e_grid) * density  # MeV/cm
    sigma = np.array(
        [nuclear_ref.macroscopic_nonelastic(float(e), oxygen_density) for e in e_grid]
    )
    integral = float(np.trapezoid(sigma / s_lin, e_grid))
    return 1.0 - math.exp(-integral)


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
    grid = DepthDoseGrid(400.0, 2000)
    engine = TransportEngine(table, slab, grid, nuclear=True)
    engine_off = TransportEngine(table, slab, grid, nuclear=False)
    oxygen_density = (
        AVOGADRO
        * materials.WATER.atoms_per_gram("O")
        * materials.WATER.density_g_per_cm3
    )
    report["oxygen_density_per_cm3"] = oxygen_density

    # -- reference path: exact energy budget + parity baseline ----------------
    per_energy: dict[str, Any] = {}
    ref_profiles: dict[float, np.ndarray] = {}
    ref_reactions: dict[float, int] = {}
    ref_escaped: dict[float, float] = {}
    budget_ok = True
    regression_ok = True
    for e0 in ENERGIES:
        ref = engine.run(
            PencilBeamSource(e0), N_REFERENCE, seed=2024, path="python"
        )
        ref_profiles[e0] = ref.edep_mev
        ref_reactions[e0] = ref.n_reactions
        ref_escaped[e0] = ref.escaped_mev
        budget_ok = budget_ok and abs(ref.energy_balance) <= 1e-9
        # nuclear=False draws no extra uniform: zero reactions, exact EM budget.
        off = engine_off.run(
            PencilBeamSource(e0), N_REFERENCE, seed=2024, path="python"
        )
        regression_ok = (
            regression_ok
            and off.n_reactions == 0
            and off.escaped_mev == 0.0
            and abs(off.energy_balance) <= 1e-9
        )
        per_energy[str(e0)] = {
            "reference_reactions": ref.n_reactions,
            "reference_survival": 1.0 - ref.n_reactions / N_REFERENCE,
            "reference_escaped_mev": ref.escaped_mev,
            "reference_energy_balance": ref.energy_balance,
        }
    gates["energy_conservation_reference"] = budget_ok
    gates["nuclear_off_regression"] = regression_ok

    # -- Warp CPU / CUDA: survival, self-consistency, parity ------------------
    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        gates["primary_survival"] = False
        gates["reaction_self_consistency"] = False
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
        survival_ok = True
        selfcons_ok = True
        cpu_ref_ok = True
        cuda_ref_ok = True
        warp_profiles: dict[str, dict[float, np.ndarray]] = {d: {} for d in devices}
        for e0 in ENERGIES:
            analytic = analytic_reaction_fraction(table, e0, oxygen_density)
            for device in devices:
                res = engine.run(
                    PencilBeamSource(e0),
                    N_SURVIVAL,
                    seed=99,
                    path="warp",
                    device=device,
                )
                warp_profiles[device][e0] = res.edep_mev
                survival = 1.0 - res.n_reactions / N_SURVIVAL
                mc_fraction = res.n_reactions / N_SURVIVAL
                stat = math.sqrt(mc_fraction * (1.0 - mc_fraction) / N_SURVIVAL)
                target, tol = PUBLISHED_SURVIVAL[e0]
                entry = {
                    "survival": survival,
                    "reactions": res.n_reactions,
                    "escaped_mev": res.escaped_mev,
                    "energy_balance": res.energy_balance,
                    "mc_reaction_fraction": mc_fraction,
                    "analytic_reaction_fraction": analytic,
                }
                report["warp"].setdefault(device, {})[str(e0)] = entry
                if device == "cpu":
                    survival_ok = survival_ok and abs(survival - target) <= tol
                    selfcons_ok = selfcons_ok and abs(
                        mc_fraction - analytic
                    ) <= 5.0 * stat + 0.01
                    # exact-parity: reference vs Warp CPU at matched N and seed;
                    # aligned counter RNG must remove the identical primary set.
                    cpu_small = engine.run(
                        PencilBeamSource(e0),
                        N_REFERENCE,
                        seed=2024,
                        path="warp",
                        device="cpu",
                    )
                    parity_cum = cumulative_diff(
                        cpu_small.edep_mev, ref_profiles[e0]
                    )
                    entry["parity_reactions_ref"] = ref_reactions[e0]
                    entry["parity_reactions_cpu"] = cpu_small.n_reactions
                    entry["parity_cumulative_diff"] = parity_cum
                    cpu_ref_ok = (
                        cpu_ref_ok
                        and abs(res.energy_balance) <= 1e-5
                        and cpu_small.n_reactions == ref_reactions[e0]
                        and math.isclose(
                            cpu_small.escaped_mev, ref_escaped[e0], rel_tol=1e-4
                        )
                        and parity_cum <= REF_WARP_CUM_TOL
                    )
                else:
                    cuda_ref_ok = cuda_ref_ok and abs(res.energy_balance) <= 1e-5
        gates["primary_survival"] = survival_ok
        gates["reaction_self_consistency"] = selfcons_ok
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
                    cum = cumulative_diff(
                        warp_profiles[device][e0], warp_profiles["cpu"][e0]
                    )
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
