"""V3 (scattering): density-heterogeneous 3-D multiple-scattering transport.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached). Emits one JSON document with the decision-0016 gates:

* ``homogeneous_equivalence`` - a uniform water VoxelSlab reproduces the
                               homogeneous WaterSlab scattering result;
* ``density_scaling_sigma_x`` - a uniform slab of density rho reproduces the
                               Fermi-Eyges lateral sigma_x at that density;
* ``layered_sigma_x``        - a water/dense/water phantom reproduces the
                               piecewise-density Fermi-Eyges sigma_x;
* ``energy_conservation``    - deposited = energy_in in the layered scattering run;
* ``warp_cpu_vs_reference`` - CPU agrees with the reference Python transport on
                               sigma_x (tight) and the depth dose (3-D float32
                               budget) across the interface;
* ``warp_cuda_vs_oracle`` / ``warp_cpu_vs_cuda`` - CUDA reproduces the Fermi-Eyges
                               sigma_x (its own physics gate) and CUDA vs CPU
                               agree on the depth dose.

Exit 0 if every gate passes, 3 otherwise, 4 if the dataset is missing.
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
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import fermi_eyges as fe
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.depth_dose import DepthLateralGrid
from ionmc.transport.geometry import VoxelSlab

X0 = materials.WATER.radiation_length_g_per_cm2
SIGMA_TOL = 0.03  # 3 % vs Fermi-Eyges
N_LARGE = 40000
DEPTH_MM = 250.0
LAYERS = [(60.0, 1.0), (30.0, 1.4), (160.0, 1.0)]
REF_WARP_DD_TOL = 5e-4  # 3-D depth-dose float32 budget
SIGMA_BACKEND_TOL_MM = 0.05


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(
        depth_mm=depth, n_depth=int(depth * 2), half_width_mm=25.0, n_lateral=600
    )


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

    model = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, "numpy")
    src = PencilBeamSource(150.0)
    grid = DepthDoseGrid(DEPTH_MM, 10)
    lat = _lat(DEPTH_MM)

    # -- homogeneous equivalence (reference, bit-exact) -----------------------
    homo = TransportEngine(table, WaterSlab(DEPTH_MM), grid).run_scattering(
        src, lat, 400, seed=3, path="python"
    )
    vox = TransportEngine(
        table, VoxelSlab.uniform(DEPTH_MM, 1.0, 100), grid
    ).run_scattering(src, lat, 400, seed=3, path="python")
    equiv_diff = float(np.max(np.abs(homo.edep_zx_mev - vox.edep_zx_mev)))
    gates["homogeneous_equivalence"] = equiv_diff == 0.0
    report["homogeneous_equivalence_max_diff"] = equiv_diff

    # -- layered energy conservation (reference) ------------------------------
    layered = VoxelSlab.from_layers(LAYERS)
    zb, dens = layered.voxel_profile()
    ref_lay = TransportEngine(table, layered, grid).run_scattering(
        src, lat, 400, seed=7, path="python"
    )
    gates["energy_conservation_reference"] = abs(ref_lay.energy_balance) <= 1e-9
    report["layered_energy_balance"] = ref_lay.energy_balance

    # -- Warp: sigma_x vs Fermi-Eyges + cross-backend -------------------------
    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        for g in (
            "density_scaling_sigma_x",
            "layered_sigma_x",
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
        # density scaling: uniform rho=1.2
        rho = 1.2
        uni = TransportEngine(
            table, VoxelSlab.uniform(140.0, rho, 1), DepthDoseGrid(140.0, 10)
        ).run_scattering(src, _lat(140.0), N_LARGE, seed=11, path="warp", device="cpu")
        r0 = float(model.csda_range(150.0)[0]) * 10.0
        density_ok = True
        density_detail: dict[str, Any] = {}
        for frac in (0.5, 0.8):
            z = frac * r0 / rho
            mc = uni.sigma_x_at_depth(z)
            oracle = float(
                fe.lateral_sigma_x_mm(model, 150.0, np.array([z]), X0, rho)[0]
            )
            density_detail[f"{frac}R"] = {"mc": mc, "oracle": oracle}
            density_ok = density_ok and abs(mc / oracle - 1.0) <= SIGMA_TOL
        report["density_scaling"] = density_detail
        gates["density_scaling_sigma_x"] = density_ok

        # each device's layered sigma_x is checked against the heterogeneous
        # Fermi-Eyges oracle (a genuine per-device physics gate), then the
        # backends are compared to each other and (CPU) to the reference.
        profiles: dict[str, Any] = {}
        layered_detail: dict[str, Any] = {}

        def _sigma_x_vs_oracle_ok(res: Any, tag: str) -> bool:
            ok = True
            for z in (60.0, 90.0, 120.0):
                mc = res.sigma_x_at_depth(z)
                oracle = float(
                    fe.lateral_sigma_x_heterogeneous_mm(
                        model, 150.0, np.array([z]), X0, zb, dens
                    )[0]
                )
                layered_detail.setdefault(tag, {})[str(z)] = {
                    "mc": mc,
                    "oracle": oracle,
                }
                ok = ok and abs(mc / oracle - 1.0) <= SIGMA_TOL
            return ok and abs(res.energy_balance) <= 1e-5

        layered_ok = True
        for device in devices:
            res = TransportEngine(table, layered, grid).run_scattering(
                src, lat, N_LARGE, seed=11, path="warp", device=device
            )
            profiles[device] = res
            layered_ok = layered_ok and _sigma_x_vs_oracle_ok(res, device)
        report["layered_sigma_x"] = layered_detail
        gates["layered_sigma_x"] = layered_ok

        # CPU vs the reference Python transport at a matched, tractable N (the
        # slow reference loop is small; the sigma_x-vs-oracle check above uses the
        # fast large-N Warp runs).
        n_match = 4000
        eng_match = TransportEngine(table, layered, grid)
        ref_match = eng_match.run_scattering(src, lat, n_match, seed=13, path="python")
        cpu_match = eng_match.run_scattering(
            src, lat, n_match, seed=13, path="warp", device="cpu"
        )
        dd_ref = ref_match.depth_dose_mev
        cum_cpu_ref = float(
            np.max(np.abs(np.cumsum(cpu_match.depth_dose_mev) - np.cumsum(dd_ref)))
            / np.sum(dd_ref)
        )
        sig_cpu_ref = all(
            abs(cpu_match.sigma_x_at_depth(z) - ref_match.sigma_x_at_depth(z))
            <= SIGMA_BACKEND_TOL_MM
            for z in (60.0, 90.0, 120.0)
        )
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "depth_dose_cumulative": cum_cpu_ref,
            "sigma_x_ok": sig_cpu_ref,
        }
        gates["warp_cpu_vs_reference"] = (
            cum_cpu_ref <= REF_WARP_DD_TOL and sig_cpu_ref
        )

        cuda_devices = [d for d in devices if d != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cuda_vs_oracle"] = False
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            # CUDA gets its own physics gate (sigma_x vs oracle, folded into
            # layered_sigma_x above); here CUDA is compared to CPU cumulatively.
            cpu_cuda_ok = True
            cuda_oracle_ok = True
            for device in cuda_devices:
                res = profiles[device]
                cuda_oracle_ok = cuda_oracle_ok and _sigma_x_vs_oracle_ok(res, device)
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
