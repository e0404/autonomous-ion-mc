"""V3 (opening): 1-D voxelized density-heterogeneous depth-dose transport.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached, see ``python -m ionmc.data acquire``). Emits one JSON document
with the decision-0014 gates:

* ``homogeneous_equivalence`` - a uniform voxel slab (density 1.0) reproduces the
                                homogeneous WaterSlab depth dose (reference exact,
                                Warp within the decision-0001 cumulative tolerance);
* ``wet_density_scaling``     - a uniform slab of density rho puts R80 at
                                R80_water / rho (water-equivalent thickness) at
                                150 and 200 MeV for rho in {0.5, 1.2};
* ``layered_interface``       - a dense layer shifts the Bragg peak proximally by
                                its extra water-equivalent thickness;
* ``energy_conservation``     - deposited + escaped = energy_in in a heterogeneous
                                phantom with nuclear + secondary transport;
* ``warp_cpu_vs_reference`` / ``warp_cuda_vs_reference`` / ``warp_cpu_vs_cuda``
                                - the backends agree across a density interface.

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
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.geometry import VoxelSlab

ENERGIES = [150.0, 200.0]
DENSITIES = [0.5, 1.2]
WET_DEPTH_MM = 700.0  # deep enough for the 1/rho-scaled range at 200 MeV
WET_TOL = 3e-3
SHIFT_TOL_MM = 0.5
REF_WARP_CUM_TOL = 1e-4
CPU_CUDA_CUM_TOL = 1e-5
#: dense interface layer for the interface / conservation / cross-backend gates.
LAYERS = [(50.0, 1.0), (20.0, 1.85), (330.0, 1.0)]
LAYER_SHIFT_MM = 20.0 * (1.85 - 1.0)  # extra water-equivalent thickness [mm]


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

    grid = DepthDoseGrid(400.0, 800)
    fine = DepthDoseGrid(WET_DEPTH_MM, 7000)  # 0.1 mm bins for the range shift

    # -- homogeneous equivalence (reference, bit-exact) -----------------------
    src = PencilBeamSource(150.0)
    homo = TransportEngine(table, WaterSlab(400.0), grid).run(
        src, 200, seed=3, path="python"
    )
    vox = TransportEngine(table, VoxelSlab.uniform(400.0, 1.0, 200), grid).run(
        src, 200, seed=3, path="python"
    )
    equiv_max_diff = float(np.max(np.abs(homo.edep_mev - vox.edep_mev)))
    gates["homogeneous_equivalence"] = equiv_max_diff == 0.0
    report["homogeneous_equivalence_max_diff"] = equiv_max_diff

    # -- water-equivalent-thickness scaling (reference, no straggling) --------
    wet: dict[str, Any] = {}
    wet_ok = True
    for e0 in ENERGIES:
        water = TransportEngine(
            table, WaterSlab(WET_DEPTH_MM), fine, straggling=False
        ).run(PencilBeamSource(e0), 1, path="python")
        r80_water = water.r80_mm()
        for rho in DENSITIES:
            dense = TransportEngine(
                table, VoxelSlab.uniform(WET_DEPTH_MM, rho, 1), fine, straggling=False
            ).run(PencilBeamSource(e0), 1, path="python")
            expected = r80_water / rho
            rel = abs(dense.r80_mm() / expected - 1.0)
            wet[f"{e0}/{rho}"] = {
                "r80_mm": dense.r80_mm(),
                "expected_mm": expected,
                "rel_err": rel,
            }
            wet_ok = wet_ok and rel <= WET_TOL
    report["wet_scaling"] = wet
    gates["wet_density_scaling"] = wet_ok

    # -- layered interface (reference, no straggling) -------------------------
    shift_ok = True
    interface: dict[str, Any] = {}
    layered_geom = VoxelSlab.from_layers(LAYERS)
    for e0 in ENERGIES:
        water = TransportEngine(table, WaterSlab(400.0), fine, straggling=False).run(
            PencilBeamSource(e0), 1, path="python"
        )
        lay = TransportEngine(table, layered_geom, fine, straggling=False).run(
            PencilBeamSource(e0), 1, path="python"
        )
        shift = water.r80_mm() - lay.r80_mm()
        interface[str(e0)] = {"shift_mm": shift, "expected_mm": LAYER_SHIFT_MM}
        shift_ok = shift_ok and abs(shift - LAYER_SHIFT_MM) <= SHIFT_TOL_MM
    report["interface"] = interface
    gates["layered_interface"] = shift_ok

    # -- energy conservation in a heterogeneous phantom (reference) -----------
    eng_het = TransportEngine(
        table, layered_geom, grid, nuclear=True, secondaries=True
    )
    ref_het = eng_het.run(PencilBeamSource(150.0), 1000, seed=7, path="python")
    gates["energy_conservation_reference"] = (
        abs(ref_het.energy_balance) <= 1e-9
        and ref_het.n_reactions > 0
        and ref_het.n_secondaries > 0
    )
    report["heterogeneous_reference"] = {
        "energy_balance": ref_het.energy_balance,
        "reactions": ref_het.n_reactions,
        "secondaries": ref_het.n_secondaries,
        "r80_mm": ref_het.r80_mm(),
    }

    # -- cross-backend across the density interface ---------------------------
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
        profiles: dict[str, np.ndarray] = {}
        cpu_ref_ok = True
        cuda_ref_ok = True
        for device in devices:
            res = eng_het.run(
                PencilBeamSource(150.0), 4000, seed=7, path="warp", device=device
            )
            profiles[device] = res.edep_mev
            # cumulative difference vs the reference at matched N and seed
            ref_match = eng_het.run(
                PencilBeamSource(150.0), 4000, seed=7, path="python"
            )
            cum = cumulative_diff(res.edep_mev, ref_match.edep_mev)
            entry = {
                "reactions": res.n_reactions,
                "reactions_ref": ref_match.n_reactions,
                "energy_balance": res.energy_balance,
                "cumulative_diff_vs_reference": cum,
                "r80_mm": res.r80_mm(),
            }
            report["warp"].setdefault(device, {})["150.0"] = entry
            if device == "cpu":
                cpu_ref_ok = (
                    cpu_ref_ok
                    and abs(res.energy_balance) <= 1e-5
                    and res.n_reactions == ref_match.n_reactions
                    and cum <= REF_WARP_CUM_TOL
                )
            else:
                cuda_ref_ok = cuda_ref_ok and cum <= REF_WARP_CUM_TOL
        gates["warp_cpu_vs_reference"] = cpu_ref_ok
        cuda_devices = [d for d in devices if d != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cuda_vs_reference"] = False
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            gates["warp_cuda_vs_reference"] = cuda_ref_ok
            cpu_cuda_ok = True
            for device in cuda_devices:
                cum = cumulative_diff(profiles[device], profiles["cpu"])
                report["warp"].setdefault(f"{device}_vs_cpu", {})["150.0"] = cum
                cpu_cuda_ok = cpu_cuda_ok and cum <= CPU_CUDA_CUM_TOL
            gates["warp_cpu_vs_cuda"] = cpu_cuda_ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
