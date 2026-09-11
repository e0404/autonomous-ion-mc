"""V4 (LET): dose-averaged LET (LET_d) scoring on the voxel-grid path.

Runs on the controlled host runner (numpy + Warp only; the PSTAR water table must
be cached). Emits one JSON document with the decision-0023 gates:

* ``thin_voxel_analytic`` - a monoenergetic proton's thin-voxel LET_d equals the
                            tabulated electronic linear stopping power S(E)
                            [keV/um] at several energies (NIST PSTAR water);
* ``dose_weight_consistency`` - the LET_d denominator is the scored dose energy
                            (num > 0 exactly where dose > 0);
* ``distal_let_peak``     - LET_d rises with depth, its peak lies at/distal to the
                            Bragg dose peak, entrance ~0.5 keV/um at 150 MeV;
* ``voxel_size_convergence`` - the plateau LET_d is stable across dose-voxel
                            z-sizes (the property that motivates "Method C");
* ``warp_cpu_vs_reference`` - the deterministic CPU LET_d matches the reference
                            per voxel (float32 budget);
* ``warp_cpu_vs_cuda``    - the deterministic CUDA LET_d agrees with CPU per voxel
                            and the stochastic total agrees (per-voxel stochastic
                            diff recorded as a diagnostic, decision 0022).

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
from ionmc.backend import mathlib, reference
from ionmc.data import MCSQUARE_PSTAR_WATER, DatasetNotCached, cache
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
)
from ionmc.transport.depth_dose import DepthLateralGrid

N_LARGE = 20000
#: NIST PSTAR electronic linear stopping power of protons in water [keV/um]
#: (== S_mass [MeV cm^2/g] * rho, and 1 MeV/mm == 1 keV/um), the thin-voxel targets.
PSTAR_LET_KEV_UM = {250.0: 0.3911, 150.0: 0.5445, 100.0: 0.7289}


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(depth, int(depth * 2), 30.0, 400)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (21, 21, 400), (1.0, 1.0, 0.5), 1.0, origin_mm=(-10.5, -10.5, 0.0)
    )


def _dose(nzv: int = 100, sp: float = 2.0) -> DoseGrid3D:
    return DoseGrid3D(
        shape=(21, 21, nzv),
        origin_mm=(-10.5, -10.5, 0.0),
        spacing_mm=(1.0, 1.0, 200.0 / nzv),
    )


def _s_lin(table, energy: float) -> float:
    tp = reference.load_bound_module(
        "ionmc.physics.transport",
        "python",
        rebind_dependencies=["ionmc.physics.tabulated"],
    )
    return float(
        tp.linear_stopping_power(
            energy,
            1.0,
            table.energy_mev,
            table.stopping_mev_cm2_per_g,
            table.slope,
            table.size,
            table.bisection_steps,
        )
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

    # -- thin-voxel analytic limit vs NIST PSTAR -----------------------------
    eng_thin = TransportEngine(
        table,
        _box(),
        DepthDoseGrid(200.0, 10),
        straggling=False,
        scattering=False,
        max_step_mm=0.1,
    )
    thin: dict[str, Any] = {}
    ok_thin = True
    for e, pstar in PSTAR_LET_KEV_UM.items():
        dose = DoseGrid3D(
            shape=(1, 1, 1), origin_mm=(-1.5, -1.5, 0.0), spacing_mm=(3.0, 3.0, 1.0)
        )
        r = eng_thin.run_scattering(
            PencilBeamSource(e),
            _lat(200.0),
            1,
            seed=1,
            path="python",
            dose_grid=dose,
            score_let=True,
        )
        letd = float(dose.let_d_kev_um(r.let3d_num_mev_per_mm, r.dose3d_mev)[0, 0, 0])
        rel_pstar = abs(letd - pstar) / pstar
        rel_slin = abs(letd - _s_lin(table, e)) / _s_lin(table, e)
        thin[f"{e:.0f}MeV"] = {
            "let_d_kev_um": letd,
            "pstar_kev_um": pstar,
            "rel_vs_pstar": rel_pstar,
            "rel_vs_s_lin": rel_slin,
        }
        ok_thin = ok_thin and rel_pstar <= 0.02 and letd >= _s_lin(table, e)
    report["thin_voxel"] = thin
    gates["thin_voxel_analytic"] = ok_thin

    # -- dose-weight consistency & distal LET trend --------------------------
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    dose1 = _dose(nzv=200)
    rr = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(200.0),
        1,
        seed=5,
        path="python",
        dose_grid=dose1,
        score_let=True,
    )
    pos = rr.dose3d_mev > 0.0
    consistent = bool(
        np.all(rr.let3d_num_mev_per_mm[pos] > 0.0)
        and np.all(rr.let3d_num_mev_per_mm[~pos] == 0.0)
    )
    report["dose_weight_consistent"] = consistent
    gates["dose_weight_consistency"] = consistent

    zc, ddz = dose1.axis_marginals(rr.dose3d_mev)
    _, dnz = dose1.axis_marginals(rr.let3d_num_mev_per_mm)
    letz = np.divide(dnz, ddz, out=np.zeros_like(ddz), where=ddz > 0)
    sig = ddz > 0.01 * ddz.max()
    k_dose = int(np.argmax(ddz))
    k_let = int(np.argmax(np.where(sig, letz, 0.0)))
    entrance = float(letz[int(np.argmax(sig))])
    report["distal"] = {
        "dose_peak_z_mm": float(zc[k_dose]),
        "let_peak_z_mm": float(zc[k_let]),
        "entrance_let_kev_um": entrance,
        "peak_let_kev_um": float(letz[k_let]),
    }
    gates["distal_let_peak"] = bool(
        zc[k_let] >= zc[k_dose]
        and abs(entrance - 0.5445) / 0.5445 <= 0.05
        and letz[k_let] > 3.0
    )

    # -- voxel-size convergence (the property that motivates Method C) --------
    z0 = 40.0  # plateau depth, well before the ~157 mm Bragg peak
    conv_vals = []
    for sp in (0.5, 1.0, 2.0):
        dg = DoseGrid3D(
            shape=(21, 21, int(200.0 / sp)),
            origin_mm=(-10.5, -10.5, 0.0),
            spacing_mm=(1.0, 1.0, sp),
        )
        rc = eng.run_scattering(
            PencilBeamSource(150.0),
            _lat(200.0),
            1,
            seed=5,
            path="python",
            dose_grid=dg,
            score_let=True,
        )
        zcc, dd = dg.axis_marginals(rc.dose3d_mev)
        _, nn = dg.axis_marginals(rc.let3d_num_mev_per_mm)
        lz = np.divide(nn, dd, out=np.zeros_like(dd), where=dd > 0)
        conv_vals.append(float(lz[int(np.argmin(np.abs(zcc - z0)))]))
    conv_spread = (max(conv_vals) - min(conv_vals)) / min(conv_vals)
    report["voxel_size_convergence"] = {
        "z0_mm": z0,
        "let_d_kev_um_by_spacing": dict(
            zip(["0.5", "1.0", "2.0"], conv_vals, strict=True)
        ),
        "rel_spread": conv_spread,
    }
    gates["voxel_size_convergence"] = bool(
        all(v > 0.0 for v in conv_vals) and conv_spread < 0.03
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

        def _letd(res: Any) -> np.ndarray:
            return _dose().let_d_kev_um(res.let3d_num_mev_per_mm, res.dose3d_mev)

        # deterministic reference vs Warp CPU LET_d per voxel (float32 budget)
        ref_d = eng.run_scattering(
            PencilBeamSource(150.0),
            _lat(200.0),
            1,
            seed=5,
            path="python",
            dose_grid=_dose(),
            score_let=True,
        )
        cpu_d = eng.run_scattering(
            PencilBeamSource(150.0),
            _lat(200.0),
            1,
            seed=5,
            path="warp",
            device="cpu",
            dose_grid=_dose(),
            score_let=True,
        )
        lr, lc = _letd(ref_d), _letd(cpu_d)
        mask = lr > 0.0
        cpu_max = float(np.max(np.abs(lc[mask] - lr[mask]) / lr[mask]))
        report.setdefault("cross_backend", {})["cpu_vs_reference_max_voxel"] = cpu_max
        gates["warp_cpu_vs_reference"] = cpu_max <= 5e-3

        cuda_devices = [dv for dv in devices if dv != "cpu"]
        if args.require_cuda and not cuda_devices:
            gates["warp_cpu_vs_cuda"] = False
        elif cuda_devices:
            dev0 = cuda_devices[0]
            # deterministic per-voxel LET_d parity (the discriminating spatial gate)
            gpu_det = eng.run_scattering(
                PencilBeamSource(150.0),
                _lat(200.0),
                1,
                seed=5,
                path="warp",
                device=dev0,
                dose_grid=_dose(),
                score_let=True,
            )
            lg = _letd(gpu_det)
            det_max = float(np.max(np.abs(lg[mask] - lc[mask]) / lc[mask]))
            # stochastic: total LET numerator agrees; per-voxel LET_d is a diagnostic
            eng_s = TransportEngine(table, _box(), DepthDoseGrid(200.0, 10))
            cpu_s = eng_s.run_scattering(
                PencilBeamSource(150.0),
                _lat(200.0),
                N_LARGE,
                seed=9,
                path="warp",
                device="cpu",
                dose_grid=_dose(),
                score_let=True,
            )
            ok = det_max <= 5e-3
            cb = report.setdefault("cross_backend", {})
            cb[f"{dev0}_vs_cpu_deterministic_max_voxel"] = det_max
            ns = float(cpu_s.let3d_num_mev_per_mm.sum())
            ls_c = _letd(cpu_s)
            for device in cuda_devices:
                gpu_s = eng_s.run_scattering(
                    PencilBeamSource(150.0),
                    _lat(200.0),
                    N_LARGE,
                    seed=9,
                    path="warp",
                    device=device,
                    dose_grid=_dose(),
                    score_let=True,
                )
                num_tot = abs(float(gpu_s.let3d_num_mev_per_mm.sum()) - ns) / ns
                ls_g = _letd(gpu_s)
                m2 = ls_c > 0.0
                stoch_max = float(np.max(np.abs(ls_g[m2] - ls_c[m2]) / ls_c[m2]))
                cb[f"{device}_vs_cpu_num_total"] = num_tot
                cb[f"{device}_vs_cpu_max_voxel_stochastic"] = stoch_max
                ok = ok and num_tot <= 1e-4
            gates["warp_cpu_vs_cuda"] = ok

    report["all_gates_passed"] = all(gates.values())
    json.dump(report, sys.stdout, indent=2, default=_jsonable)
    sys.stdout.write("\n")
    return 0 if report["all_gates_passed"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
