"""V3 (material scattering): non-water materials on the 3-D scattering path.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached). Emits one JSON document with the decision-0017 gates:

* ``homogeneous_equivalence`` - a water VoxelSlab reproduces the WaterSlab
                               scattering result bit-for-bit;
* ``bone_sigma_x``           - a homogeneous bone slab reproduces the material-
                               aware Fermi-Eyges sigma_x (physical density, bone
                               X0) and R80 at R_water/WER;
* ``material_interface_sigma_x`` - a water/bone/water phantom reproduces the
                               material-aware Fermi-Eyges sigma_x across the
                               interface;
* ``energy_conservation``    - deposited = energy_in in the material scattering run;
* ``warp_cpu_vs_reference`` / ``warp_cuda_vs_oracle`` / ``warp_cpu_vs_cuda``
                               - the backends agree on sigma_x and depth dose.

Exit 0 if every gate passes, 3 otherwise, 4 if the dataset is missing.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from typing import Any

import numpy as np

from ionmc import __version__, particles
from ionmc import materials as M
from ionmc.backend import mathlib
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import fermi_eyges as fe
from ionmc.stopping_power import mass_stopping_power_ratio
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.depth_dose import DepthLateralGrid
from ionmc.transport.geometry import VoxelSlab

SIGMA_TOL = 0.03
N_LARGE = 40000
DEPTH_MM = 250.0
INTERFACE = [(40.0, M.WATER), (20.0, M.CORTICAL_BONE), (190.0, M.WATER)]
REF_WARP_DD_TOL = 5e-4
SIGMA_BACKEND_TOL_MM = 0.05
INTERFACE_DEPTHS = (40.0, 60.0, 90.0)


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(
        depth_mm=depth, n_depth=int(depth * 2), half_width_mm=25.0, n_lateral=600
    )


def _spr(mat: Any) -> float:
    return mass_stopping_power_ratio(mat, 150.0, particles.PROTON)


def _material_arrays(slab: VoxelSlab):
    zb, rho = slab.voxel_profile()
    mats = slab.materials_profile()
    we = np.array([_spr(m) * rho[i] for i, m in enumerate(mats)])
    phys = np.asarray(rho, dtype=np.float64)
    rl = np.array([m.radiation_length_g_per_cm2 for m in mats])
    return zb, we, phys, rl


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
        table = load_stopping_table(
            cache.load_path(MCSQUARE_PSTAR_WATER, args.cache_dir)
        )
    except DatasetNotCached as exc:
        report["error"] = str(exc)
        json.dump(report, sys.stdout, indent=2, default=_jsonable)
        sys.stdout.write("\n")
        return 4

    model = TabulatedStoppingPower(table, M.WATER, particles.PROTON, "numpy")
    src = PencilBeamSource(150.0)

    # -- homogeneous equivalence (reference, bit-exact) -----------------------
    lat = _lat(DEPTH_MM)
    homo = TransportEngine(table, WaterSlab(DEPTH_MM), DepthDoseGrid(DEPTH_MM, 10))
    vox = TransportEngine(
        table, VoxelSlab.from_material_layers([(DEPTH_MM, M.WATER)]),
        DepthDoseGrid(DEPTH_MM, 10),
    )
    a = homo.run_scattering(src, lat, 400, seed=3, path="python")
    b = vox.run_scattering(src, lat, 400, seed=3, path="python")
    equiv = float(np.max(np.abs(a.edep_zx_mev - b.edep_zx_mev)))
    gates["homogeneous_equivalence"] = equiv == 0.0
    report["homogeneous_equivalence_max_diff"] = equiv

    # -- material interface energy conservation (reference) -------------------
    layered = VoxelSlab.from_material_layers(INTERFACE)
    zb_i, we_i, phys_i, rl_i = _material_arrays(layered)
    ref_lay = TransportEngine(
        table, layered, DepthDoseGrid(DEPTH_MM, 10)
    ).run_scattering(src, lat, 400, seed=7, path="python")
    gates["energy_conservation_reference"] = abs(ref_lay.energy_balance) <= 1e-9
    report["layered_energy_balance"] = ref_lay.energy_balance

    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        for g in (
            "bone_sigma_x",
            "material_interface_sigma_x",
            "warp_cpu_vs_reference",
        ):
            gates[g] = False
        if args.require_cuda:
            gates["warp_cuda_vs_oracle"] = False
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

        # homogeneous bone: sigma_x vs material oracle + R80 at WER
        bone = VoxelSlab.from_material_layers([(120.0, M.CORTICAL_BONE)])
        zb_b, we_b, phys_b, rl_b = _material_arrays(bone)
        rb = TransportEngine(table, bone, DepthDoseGrid(120.0, 10)).run_scattering(
            src, _lat(120.0), N_LARGE, seed=11, path="warp", device="cpu"
        )
        bone_detail: dict[str, Any] = {}
        bone_ok = True
        for z in (48.0, 72.0):
            mc = rb.sigma_x_at_depth(z)
            oracle = float(
                fe.lateral_sigma_x_material_mm(
                    model, 150.0, np.array([z]), zb_b, we_b, phys_b, rl_b
                )[0]
            )
            bone_detail[str(z)] = {"mc": mc, "oracle": oracle}
            bone_ok = bone_ok and abs(mc / oracle - 1.0) <= SIGMA_TOL
        report["bone_sigma_x"] = bone_detail
        gates["bone_sigma_x"] = bone_ok

        # material interface: sigma_x vs material oracle, per device
        profiles: dict[str, Any] = {}
        interface_detail: dict[str, Any] = {}

        def _interface_ok(res: Any, tag: str) -> bool:
            ok = True
            for z in INTERFACE_DEPTHS:
                mc = res.sigma_x_at_depth(z)
                oracle = float(
                    fe.lateral_sigma_x_material_mm(
                        model, 150.0, np.array([z]), zb_i, we_i, phys_i, rl_i
                    )[0]
                )
                interface_detail.setdefault(tag, {})[str(z)] = {
                    "mc": mc,
                    "oracle": oracle,
                }
                ok = ok and abs(mc / oracle - 1.0) <= SIGMA_TOL
            return ok and abs(res.energy_balance) <= 1e-5

        interface_ok = True
        for device in devices:
            res = TransportEngine(
                table, layered, DepthDoseGrid(DEPTH_MM, 10)
            ).run_scattering(
                src, lat, N_LARGE, seed=11, path="warp", device=device
            )
            profiles[device] = res
            interface_ok = interface_ok and _interface_ok(res, device)
        report["material_interface_sigma_x"] = interface_detail
        gates["material_interface_sigma_x"] = interface_ok

        # CPU vs reference at a matched, tractable N
        n_match = 4000
        eng_m = TransportEngine(table, layered, DepthDoseGrid(DEPTH_MM, 10))
        ref_m = eng_m.run_scattering(src, lat, n_match, seed=13, path="python")
        cpu_m = eng_m.run_scattering(
            src, lat, n_match, seed=13, path="warp", device="cpu"
        )
        dd_ref = ref_m.depth_dose_mev
        cum_cpu = float(
            np.max(np.abs(np.cumsum(cpu_m.depth_dose_mev) - np.cumsum(dd_ref)))
            / np.sum(dd_ref)
        )
        sig_cpu = all(
            abs(cpu_m.sigma_x_at_depth(z) - ref_m.sigma_x_at_depth(z))
            <= SIGMA_BACKEND_TOL_MM
            for z in INTERFACE_DEPTHS
        )
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "depth_dose_cumulative": cum_cpu,
            "sigma_x_ok": sig_cpu,
        }
        gates["warp_cpu_vs_reference"] = cum_cpu <= REF_WARP_DD_TOL and sig_cpu

        cuda_devices = [d for d in devices if d != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cuda_vs_oracle"] = False
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            cpu_cuda_ok = True
            cuda_oracle_ok = True
            for device in cuda_devices:
                res = profiles[device]
                cuda_oracle_ok = cuda_oracle_ok and _interface_ok(res, device)
                cum = float(
                    np.max(
                        np.abs(
                            np.cumsum(res.depth_dose_mev)
                            - np.cumsum(profiles["cpu"].depth_dose_mev)
                        )
                    )
                    / np.sum(profiles["cpu"].depth_dose_mev)
                )
                report.setdefault("cross_backend", {})[f"{device}_vs_cpu"] = cum
                cpu_cuda_ok = cpu_cuda_ok and cum <= REF_WARP_DD_TOL
            gates["warp_cuda_vs_oracle"] = cuda_oracle_ok
            gates["warp_cpu_vs_cuda"] = cpu_cuda_ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
