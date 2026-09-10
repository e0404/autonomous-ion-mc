"""V3 (materials): per-voxel tissue materials via stopping-power ratios.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached). Emits one JSON document with the decision-0015 gates:

* ``water_backcompat``       - a water VoxelSlab reproduces the homogeneous
                               WaterSlab depth dose bit-for-bit;
* ``wer_published_bands``    - each tissue's water-equivalent ratio SPR*rho falls
                               in the published Schneider/ICRU band;
* ``wer_range_scaling``      - a tissue slab puts R80 at R80_water / WER at 150
                               and 200 MeV;
* ``nuclear_composition``    - the composition-scaled nonelastic content reduces
                               to n_O for water and is ~2x (bone) / ~3.4x
                               (adipose) larger than oxygen-only;
* ``interface_conservation`` - deposited + escaped = energy_in in a multi-material
                               phantom with nuclear + secondary transport;
* ``warp_cpu_vs_reference`` / ``warp_cuda_vs_reference`` / ``warp_cpu_vs_cuda``
                               - the backends agree across material interfaces.

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
from ionmc import materials as M
from ionmc.backend import mathlib
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.particles import PROTON
from ionmc.stopping_power import mass_stopping_power_ratio
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.geometry import VoxelSlab

ENERGIES = [150.0, 200.0]
WET_DEPTH_MM = 700.0
WET_TOL = 3e-3
REF_WARP_CUM_TOL = 1e-4
CPU_CUDA_CUM_TOL = 1e-5
#: material -> published WER band (Schneider 2000 / ICRU).
BANDS = {
    "cortical_bone": (M.CORTICAL_BONE, 1.60, 1.72),
    "adipose": (M.ADIPOSE, 0.95, 0.98),
    "soft_tissue": (M.SOFT_TISSUE, 1.01, 1.04),
    "skeletal_muscle": (M.SKELETAL_MUSCLE, 1.02, 1.05),
    "lung_tissue": (M.LUNG_TISSUE, 1.02, 1.05),
}
INTERFACE_LAYERS = [
    (40.0, M.WATER),
    (20.0, M.CORTICAL_BONE),
    (30.0, M.ADIPOSE),
    (310.0, M.WATER),
]


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
    fine = DepthDoseGrid(WET_DEPTH_MM, 7000)

    # -- water back-compatibility (reference, bit-exact) ----------------------
    src = PencilBeamSource(150.0)
    homo = TransportEngine(table, WaterSlab(400.0), grid).run(
        src, 200, seed=3, path="python"
    )
    vox = TransportEngine(
        table, VoxelSlab.from_material_layers([(400.0, M.WATER)]), grid
    ).run(src, 200, seed=3, path="python")
    backcompat_diff = float(np.max(np.abs(homo.edep_mev - vox.edep_mev)))
    gates["water_backcompat"] = backcompat_diff == 0.0
    report["water_backcompat_max_diff"] = backcompat_diff

    # -- published WER bands + nuclear composition scaling --------------------
    band_ok = True
    nuc_ok = True
    wer_report: dict[str, Any] = {}
    for name, (mat, lo, hi) in BANDS.items():
        spr = mass_stopping_power_ratio(mat, 150.0, PROTON)
        wer = spr * mat.density_g_per_cm3
        nuc_ratio = mat.oxygen_equivalent_per_gram / mat.atoms_per_gram("O")
        wer_report[name] = {"spr": spr, "wer": wer, "nuclear_over_oxygen": nuc_ratio}
        band_ok = band_ok and lo <= wer <= hi
    water_nuc = M.WATER.oxygen_equivalent_per_gram / M.WATER.atoms_per_gram("O")
    bone_nuc = wer_report["cortical_bone"]["nuclear_over_oxygen"]
    adipose_nuc = wer_report["adipose"]["nuclear_over_oxygen"]
    nuc_ok = (
        abs(water_nuc - 1.0) < 1e-9
        and 1.9 <= bone_nuc <= 2.3
        and 3.0 <= adipose_nuc <= 3.8
    )
    report["wer"] = wer_report
    gates["wer_published_bands"] = band_ok
    gates["nuclear_composition"] = nuc_ok

    # -- WER range scaling (reference, no straggling) -------------------------
    scaling_ok = True
    scaling: dict[str, Any] = {}
    for e0 in ENERGIES:
        water = TransportEngine(
            table, WaterSlab(WET_DEPTH_MM), fine, straggling=False
        ).run(PencilBeamSource(e0), 1, path="python")
        r80_water = water.r80_mm()
        for name in ("cortical_bone", "adipose", "skeletal_muscle"):
            mat = BANDS[name][0]
            slab = VoxelSlab.from_material_layers([(WET_DEPTH_MM, mat)])
            tissue = TransportEngine(table, slab, fine, straggling=False).run(
                PencilBeamSource(e0), 1, path="python"
            )
            wer = mass_stopping_power_ratio(mat, 150.0, PROTON) * mat.density_g_per_cm3
            expected = r80_water / wer
            rel = abs(tissue.r80_mm() / expected - 1.0)
            scaling[f"{e0}/{name}"] = {
                "r80_mm": tissue.r80_mm(),
                "expected_mm": expected,
                "rel_err": rel,
            }
            scaling_ok = scaling_ok and rel <= WET_TOL
    report["wer_range_scaling"] = scaling
    gates["wer_range_scaling"] = scaling_ok

    # -- interface conservation (reference) -----------------------------------
    eng_het = TransportEngine(
        table, VoxelSlab.from_material_layers(INTERFACE_LAYERS), grid,
        nuclear=True, secondaries=True,
    )
    ref_het = eng_het.run(PencilBeamSource(150.0), 1000, seed=7, path="python")
    gates["interface_conservation_reference"] = (
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

    # -- cross-backend across material interfaces -----------------------------
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
        ref_match = eng_het.run(PencilBeamSource(150.0), 4000, seed=7, path="python")
        for device in devices:
            res = eng_het.run(
                PencilBeamSource(150.0), 4000, seed=7, path="warp", device=device
            )
            profiles[device] = res.edep_mev
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
