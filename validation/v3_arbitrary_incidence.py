"""V3 (arbitrary incidence): beam-frame transport and rotated-vs-axis-aligned
beam equivalence.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table
must be cached). Emits one JSON document with the decision-0018 gates:

* ``axis_aligned_energy_conservation`` - the default +z beam-frame run conserves
                               energy (sanity that the beam-frame path is intact);
* ``deterministic_rotation_equivalence`` - scattering off, a rigidly rotated
                               scene (beam + slab normal) reproduces the
                               axis-aligned beam-frame depth dose to round-off;
* ``statistical_rotation_covariance`` - scattering on, a rotated homogeneous slab
                               reproduces the Fermi-Eyges lateral sigma_x' and
                               the axis-aligned R80 within statistics;
* ``oblique_wet_traversal``  - a beam tilted by theta through a lab-fixed slab of
                               thickness D deposits the same energy as a normal
                               beam through D/cos(theta) (same water-equivalent
                               path), exercising the m_hat != z' plane traversal;
* ``warp_cpu_vs_reference`` - CPU agrees with the reference Python transport on a
                               rotated config (sigma_x' and the 3-D depth dose);
* ``warp_cuda_vs_oracle`` / ``warp_cpu_vs_cuda`` - CUDA reproduces the Fermi-Eyges
                               sigma_x' and CUDA vs CPU agree on the depth dose.

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

X0 = materials.WATER.radiation_length_g_per_cm2
SIGMA_TOL = 0.03  # 3 % vs Fermi-Eyges
SIGMA_BACKEND_TOL_MM = 0.05
REF_WARP_DD_TOL = 5e-4  # 3-D depth-dose float32 budget
N_LARGE = 40000
DEPTH_MM = 250.0
# a non-trivial rigid rotation of the whole scene (yaw about y, pitch about x)
ROT_Y_DEG = 28.0
ROT_X_DEG = -17.0


def _rot_x(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _rot_y(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _scene_direction() -> tuple[float, float, float]:
    d = _rot_y(ROT_Y_DEG) @ _rot_x(ROT_X_DEG) @ np.array([0.0, 0.0, 1.0])
    return tuple(float(v) for v in d)


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(
        depth_mm=depth, n_depth=int(depth * 2), half_width_mm=25.0, n_lateral=600
    )


def _r80_mm(res: Any) -> float:
    dd = res.depth_dose_mev
    centers = res.grid.depth_centers_mm
    peak = float(np.max(dd))
    i_peak = int(np.argmax(dd))
    level = 0.8 * peak
    for i in range(i_peak, len(dd) - 1):
        if dd[i] >= level >= dd[i + 1]:
            f = (dd[i] - level) / (dd[i] - dd[i + 1])
            return float(centers[i] + f * (centers[i + 1] - centers[i]))
    return float(centers[i_peak])


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
        "rotation_deg": {"yaw_y": ROT_Y_DEG, "pitch_x": ROT_X_DEG},
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
    grid = DepthDoseGrid(DEPTH_MM, 10)
    lat = _lat(DEPTH_MM)
    d = _scene_direction()
    report["scene_direction"] = list(d)

    # -- axis-aligned energy conservation (reference) -------------------------
    axis_ref = TransportEngine(table, WaterSlab(DEPTH_MM), grid).run_scattering(
        PencilBeamSource(150.0), lat, 400, seed=9, path="python"
    )
    gates["axis_aligned_energy_conservation"] = abs(axis_ref.energy_balance) <= 1e-9
    report["axis_aligned_energy_balance"] = axis_ref.energy_balance

    # -- deterministic rotation equivalence (scattering off, reference) -------
    kw_det = dict(straggling=False, scattering=False)
    axis_det = TransportEngine(
        table, WaterSlab(DEPTH_MM), grid, **kw_det
    ).run_scattering(PencilBeamSource(150.0), lat, 1, seed=5, path="python")
    rot_det = TransportEngine(
        table, WaterSlab(DEPTH_MM, normal=d), grid, **kw_det
    ).run_scattering(
        PencilBeamSource(150.0, direction=d), lat, 1, seed=5, path="python"
    )
    e0 = axis_det.energy_deposited_mev
    d_energy = abs(rot_det.energy_deposited_mev - e0) / e0
    dd_cum = float(
        np.max(
            np.abs(
                np.cumsum(rot_det.depth_dose_mev) - np.cumsum(axis_det.depth_dose_mev)
            )
        )
        / np.sum(axis_det.depth_dose_mev)
    )
    d_r80 = abs(_r80_mm(rot_det) - _r80_mm(axis_det))
    report["deterministic_rotation"] = {
        "energy_rel_diff": d_energy,
        "depth_dose_cumulative": dd_cum,
        "r80_diff_mm": d_r80,
    }
    gates["deterministic_rotation_equivalence"] = (
        d_energy <= 1e-12 and dd_cum <= 1e-9 and d_r80 <= 1e-4
    )

    # -- oblique WET traversal (analytic, reference, scattering off) ----------
    theta = 30.0
    cos_t = float(np.cos(np.radians(theta)))
    slab_d = 80.0  # thin: the 150 MeV proton escapes the back face
    d_obl = tuple(float(v) for v in (_rot_y(theta) @ np.array([0.0, 0.0, 1.0])))
    lat_obl = _lat(300.0)
    normal_eff = TransportEngine(
        table, WaterSlab(slab_d / cos_t), DepthDoseGrid(300.0, 10), **kw_det
    ).run_scattering(PencilBeamSource(150.0), lat_obl, 1, seed=3, path="python")
    oblique = TransportEngine(
        table, WaterSlab(slab_d), DepthDoseGrid(300.0, 10), **kw_det
    ).run_scattering(
        PencilBeamSource(150.0, direction=d_obl), lat_obl, 1, seed=3, path="python"
    )
    e_ref = normal_eff.energy_deposited_mev
    obl_diff = abs(oblique.energy_deposited_mev - e_ref) / e_ref
    report["oblique_wet"] = {
        "theta_deg": theta,
        "deposited_oblique_mev": oblique.energy_deposited_mev,
        "deposited_normal_effective_mev": e_ref,
        "rel_diff": obl_diff,
    }
    gates["oblique_wet_traversal"] = (
        0.0 < oblique.energy_deposited_mev < 150.0 and obl_diff <= 1e-3
    )

    # -- Warp: statistical covariance + cross-backend -------------------------
    r0 = float(model.csda_range(150.0)[0]) * 10.0

    def _sigma_vs_oracle(res: Any, tag: str, detail: dict[str, Any]) -> bool:
        ok = True
        for frac in (0.5, 0.8):
            z = frac * r0
            mc = res.sigma_x_at_depth(z)
            oracle = float(
                fe.lateral_sigma_x_mm(model, 150.0, np.array([z]), X0, 1.0)[0]
            )
            detail.setdefault(tag, {})[f"{frac}R"] = {"mc": mc, "oracle": oracle}
            ok = ok and abs(mc / oracle - 1.0) <= SIGMA_TOL
        return ok and abs(res.energy_balance) <= 1e-5

    if not mathlib.HAVE_WARP:
        report["warp"] = {"available": False}
        gates["statistical_rotation_covariance"] = False
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
        rot_eng = TransportEngine(table, WaterSlab(DEPTH_MM, normal=d), grid)
        axis_eng = TransportEngine(table, WaterSlab(DEPTH_MM), grid)
        rot_src = PencilBeamSource(150.0, direction=d)
        axis_src = PencilBeamSource(150.0)

        # statistical rotation covariance: rotated homogeneous slab sigma_x' vs
        # the Fermi-Eyges oracle, and axis-vs-rotated R80 within statistics.
        cov_detail: dict[str, Any] = {}
        profiles: dict[str, Any] = {}
        cov_ok = True
        for device in devices:
            res = rot_eng.run_scattering(
                rot_src, lat, N_LARGE, seed=11, path="warp", device=device
            )
            profiles[device] = res
            cov_ok = cov_ok and _sigma_vs_oracle(res, device, cov_detail)
        axis_large = axis_eng.run_scattering(
            axis_src, lat, N_LARGE, seed=11, path="warp", device="cpu"
        )
        r80_diff = abs(_r80_mm(profiles["cpu"]) - _r80_mm(axis_large))
        report["statistical_covariance"] = cov_detail
        report["statistical_covariance"]["r80_axis_vs_rot_diff_mm"] = r80_diff
        gates["statistical_rotation_covariance"] = cov_ok and r80_diff <= 1.0

        # CPU vs reference on the rotated config at a tractable N
        n_match = 4000
        ref_match = rot_eng.run_scattering(
            rot_src, lat, n_match, seed=13, path="python"
        )
        cpu_match = rot_eng.run_scattering(
            rot_src, lat, n_match, seed=13, path="warp", device="cpu"
        )
        dd_ref = ref_match.depth_dose_mev
        cum_cpu_ref = float(
            np.max(np.abs(np.cumsum(cpu_match.depth_dose_mev) - np.cumsum(dd_ref)))
            / np.sum(dd_ref)
        )
        sig_cpu_ref = all(
            abs(cpu_match.sigma_x_at_depth(f * r0) - ref_match.sigma_x_at_depth(f * r0))
            <= SIGMA_BACKEND_TOL_MM
            for f in (0.5, 0.8)
        )
        report.setdefault("cross_backend", {})["cpu_vs_reference"] = {
            "depth_dose_cumulative": cum_cpu_ref,
            "sigma_x_ok": sig_cpu_ref,
        }
        gates["warp_cpu_vs_reference"] = cum_cpu_ref <= REF_WARP_DD_TOL and sig_cpu_ref

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cuda_vs_oracle"] = False
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            cpu_cuda_ok = True
            cuda_oracle_ok = True
            for device in cuda_devices:
                res = profiles[device]
                cuda_oracle_ok = cuda_oracle_ok and _sigma_vs_oracle(
                    res, device, cov_detail
                )
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
