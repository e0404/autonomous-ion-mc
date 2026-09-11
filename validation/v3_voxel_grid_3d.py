"""V3 (voxel grid): 3-D voxel geometry with ray/voxel DDA traversal.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0020 gates:

* ``grid_reduction_to_slab`` - a single-column grid reproduces the 1-D VoxelSlab
                               depth dose (deterministic straight ray);
* ``homogeneous_box_vs_slab`` - a homogeneous 3-D box reproduces WaterSlab to the
                               step-partition discretization level;
* ``oblique_wet_vs_siddon``  - an oblique beam through an off-axis dense insert
                               stops where an independent Siddon WET integral
                               reaches the water CSDA range (interior crossing);
* ``energy_conservation``    - a beam contained in a wide box deposits all energy;
* ``warp_cpu_vs_reference`` - CPU grid kernel agrees with the reference grid
                               driver (depth dose + sigma_x);
* ``warp_cuda_vs_oracle`` / ``warp_cpu_vs_cuda`` - CUDA reproduces the Fermi-Eyges
                               sigma_x on the grid and CUDA vs CPU agree.

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
from ionmc.transport import (
    DepthDoseGrid,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
    WaterSlab,
)
from ionmc.transport.depth_dose import DepthLateralGrid
from ionmc.transport.geometry import VoxelSlab

X0 = materials.WATER.radiation_length_g_per_cm2
SIGMA_TOL = 0.03
SIGMA_BACKEND_TOL_MM = 0.05
# Same-precision (CPU vs CUDA) depth-dose cumulative budget, as on the 1-D path.
REF_WARP_DD_TOL = 5e-4
# Reference(float64) vs Warp(float32) depth-dose cumulative on the DDA grid path:
# the DDA clips at every voxel face (~160 z-crossings plus lateral crossings per
# history, vs the 1-D path's few material interfaces), so float32-vs-float64
# near-corner face-decision flips decorrelate proportionally more histories and
# accumulate ~1e-3 into the depth dose. sigma_x (physics) stays the tight
# cross-backend metric; CPU-vs-CUDA (same precision) stays at REF_WARP_DD_TOL.
GRID_REF_WARP_DD_TOL = 3e-3
N_LARGE = 40000


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(depth, int(depth * 2), 30.0, 400)


def _r80(res: Any) -> float:
    dd = res.depth_dose_mev
    c = res.grid.depth_centers_mm
    i = int(dd.argmax())
    lvl = 0.8 * float(dd.max())
    for k in range(i, len(dd) - 1):
        if dd[k] >= lvl >= dd[k + 1]:
            f = (dd[k] - lvl) / (dd[k] - dd[k + 1])
            return float(c[k] + f * (c[k + 1] - c[k]))
    return float(c[i])


def _siddon_wet_mm(grid: VoxelGrid3D, p0: Any, d: Any, path_mm: float) -> float:
    o = np.asarray(grid.origin_mm)
    h = np.asarray(grid.spacing_mm)
    n = np.array([grid.nx, grid.ny, grid.nz])
    rho = grid.density_g_per_cm3
    p0 = np.asarray(p0, dtype=float)
    d = np.asarray(d, dtype=float)
    d = d / np.linalg.norm(d)
    ds = 1.0e-3
    wet = 0.0
    s = 0.5 * ds
    while s < path_mm:
        idx = np.floor((p0 + s * d - o) / h).astype(int)
        if np.all(idx >= 0) and np.all(idx < n):
            wet += float(rho[idx[0], idx[1], idx[2]]) * ds
        s += ds
    return wet


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
    r_water = float(model.csda_range(150.0)[0]) * 10.0
    src = PencilBeamSource(150.0)
    det = dict(straggling=False, scattering=False)

    # -- reduction to VoxelSlab ----------------------------------------------
    z = np.arange(0, 301, 10.0)
    rho = 1.0 + 0.3 * np.cos(np.arange(len(z) - 1))
    slab = VoxelSlab(z_boundaries_mm=z, density_g_per_cm3=rho)
    grid_col = VoxelGrid3D.from_voxel_slab(slab, transverse_extent_mm=60.0)
    lat = _lat(300.0)
    a = TransportEngine(table, slab, DepthDoseGrid(300.0, 10), **det).run_scattering(
        src, lat, 1, seed=4, path="python"
    )
    b = TransportEngine(
        table, grid_col, DepthDoseGrid(300.0, 10), **det
    ).run_scattering(src, lat, 1, seed=4, path="python")
    red_cum = float(
        np.max(np.abs(np.cumsum(b.depth_dose_mev) - np.cumsum(a.depth_dose_mev)))
        / np.sum(a.depth_dose_mev)
    )
    report["reduction"] = {"depth_dose_cumulative": red_cum}
    gates["grid_reduction_to_slab"] = red_cum <= 1e-9

    # -- homogeneous box vs WaterSlab ----------------------------------------
    box = VoxelGrid3D.uniform(
        (80, 80, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-40.0, -40.0, 0.0)
    )
    hb = TransportEngine(table, box, DepthDoseGrid(300.0, 10), **det).run_scattering(
        src, lat, 1, seed=4, path="python"
    )
    ws = TransportEngine(
        table, WaterSlab(300.0), DepthDoseGrid(300.0, 10), **det
    ).run_scattering(src, lat, 1, seed=4, path="python")
    hb_e = (
        abs(hb.energy_deposited_mev - ws.energy_deposited_mev) / ws.energy_deposited_mev
    )
    hb_r80 = abs(_r80(hb) - _r80(ws))
    report["homogeneous_box"] = {"energy_rel_diff": hb_e, "r80_diff_mm": hb_r80}
    gates["homogeneous_box_vs_slab"] = hb_e <= 1e-4 and hb_r80 <= 1e-2

    # -- oblique WET through a dense insert vs Siddon oracle ------------------
    dens = np.ones((160, 160, 300))
    dens[80:, :, 80:180] = 1.7
    insert = VoxelGrid3D(
        dens, origin_mm=(-80.0, -80.0, 0.0), spacing_mm=(1.0, 1.0, 1.0)
    )
    ang = np.radians(18.0)
    d_obl = (float(np.sin(ang)), 0.0, float(np.cos(ang)))
    obl = TransportEngine(
        table, insert, DepthDoseGrid(300.0, 10), **det
    ).run_scattering(
        PencilBeamSource(150.0, direction=d_obl), lat, 1, seed=3, path="python"
    )
    lo, hi = 0.0, 300.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _siddon_wet_mm(insert, (0.0, 0.0, 0.0), d_obl, mid) < r_water:
            lo = mid
        else:
            hi = mid
    siddon_range = 0.5 * (lo + hi)
    obl_diff = abs(obl.range_mean_mm - siddon_range) / siddon_range
    report["oblique_wet"] = {
        "engine_range_mm": obl.range_mean_mm,
        "siddon_range_mm": siddon_range,
        "rel_diff": obl_diff,
        "n_stopped": obl.n_stopped,
    }
    gates["oblique_wet_vs_siddon"] = obl.n_stopped == 1 and obl_diff <= 1e-2

    # -- energy conservation (contained, scattering on) ----------------------
    wide = VoxelGrid3D.uniform(
        (120, 120, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-60.0, -60.0, 0.0)
    )
    econ = TransportEngine(table, wide, DepthDoseGrid(300.0, 10)).run_scattering(
        src, lat, 400, seed=7, path="python"
    )
    report["energy_balance"] = econ.energy_balance
    gates["energy_conservation"] = abs(econ.energy_balance) <= 1e-9

    # -- Warp cross-backend ---------------------------------------------------
    scatter_box = VoxelGrid3D.uniform(
        (120, 120, 250), (1.0, 1.0, 1.0), 1.0, origin_mm=(-60.0, -60.0, 0.0)
    )
    lat2 = _lat(250.0)
    r0 = r_water

    def _sigma_vs_oracle(res: Any, tag: str, detail: dict[str, Any]) -> bool:
        ok = True
        for frac in (0.5, 0.8):
            zc = frac * r0
            mc = res.sigma_x_at_depth(zc)
            oracle = float(
                fe.lateral_sigma_x_mm(model, 150.0, np.array([zc]), X0, 1.0)[0]
            )
            detail.setdefault(tag, {})[f"{frac}R"] = {"mc": mc, "oracle": oracle}
            ok = ok and abs(mc / oracle - 1.0) <= SIGMA_TOL
        return ok and abs(res.energy_balance) <= 1e-5

    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        gates["warp_cpu_vs_reference"] = False
        if args.require_cuda:
            gates["warp_cuda_vs_oracle"] = False
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
        eng = TransportEngine(table, scatter_box, DepthDoseGrid(250.0, 10))
        profiles: dict[str, Any] = {}
        for device in devices:
            profiles[device] = eng.run_scattering(
                src, lat2, N_LARGE, seed=11, path="warp", device=device
            )
        n_match = 4000
        ref_match = eng.run_scattering(src, lat2, n_match, seed=13, path="python")
        cpu_match = eng.run_scattering(
            src, lat2, n_match, seed=13, path="warp", device="cpu"
        )
        dd_ref = ref_match.depth_dose_mev
        cum_cpu_ref = float(
            np.max(np.abs(np.cumsum(cpu_match.depth_dose_mev) - np.cumsum(dd_ref)))
            / np.sum(dd_ref)
        )
        # no-drift check: the *signed* cumulative at the last bin is the total
        # deposited-energy relative difference; for two energy-conserving contained
        # runs it must be ~0, proving the ~1e-3 max excursion is a zero-mean profile
        # wiggle (float32 face-flip decorrelation), not a systematic range/energy
        # bias -- the discriminator against a subtle DDA bug (decision 0020).
        signed_drift = float(
            (np.sum(cpu_match.depth_dose_mev) - np.sum(dd_ref)) / np.sum(dd_ref)
        )
        sig_cpu_ref = all(
            abs(cpu_match.sigma_x_at_depth(f * r0) - ref_match.sigma_x_at_depth(f * r0))
            <= SIGMA_BACKEND_TOL_MM
            for f in (0.5, 0.8)
        )
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "depth_dose_cumulative": cum_cpu_ref,
            "signed_total_drift": signed_drift,
            "sigma_x_ok": sig_cpu_ref,
        }
        gates["warp_cpu_vs_reference"] = (
            cum_cpu_ref <= GRID_REF_WARP_DD_TOL
            and abs(signed_drift) <= 1e-4
            and sig_cpu_ref
        )

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cuda_vs_oracle"] = False
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            detail: dict[str, Any] = {}
            cpu_cuda_ok = True
            cuda_oracle_ok = True
            for device in cuda_devices:
                cuda_oracle_ok = cuda_oracle_ok and _sigma_vs_oracle(
                    profiles[device], device, detail
                )
                cum = float(
                    np.max(
                        np.abs(
                            np.cumsum(profiles[device].depth_dose_mev)
                            - np.cumsum(profiles["cpu"].depth_dose_mev)
                        )
                    )
                    / np.sum(profiles["cpu"].depth_dose_mev)
                )
                report.setdefault("cross_backend", {})[f"{device}_vs_cpu"] = cum
                cpu_cuda_ok = cpu_cuda_ok and cum <= REF_WARP_DD_TOL
            report["sigma_x_vs_oracle"] = detail
            gates["warp_cuda_vs_oracle"] = cuda_oracle_ok
            gates["warp_cpu_vs_cuda"] = cpu_cuda_ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
