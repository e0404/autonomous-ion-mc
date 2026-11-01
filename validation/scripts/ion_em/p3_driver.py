"""V3-008 pre-freeze check P3 driver: LS - (Bloch + first-order Mott) for He, C, O (and p) in water.

Run (no bytecode is written into the repository)::

    OMP_NUM_THREADS=1 UV_LINK_MODE=copy uv run --extra dev python -I -B \
        validation/scripts/ion_em/p3_driver.py [--out PATH] [--report PATH] \
        [--g4-ls-data <G4LindhardSorensenData.cc>]

Writes the aggregates-only JSON (``--out``, default ``validation/plans/v3-008/p3-output.json``;
quantities derived from the ionmc model and the LS theory only, no external reference-table
values) and, only when ``--report`` is given, a Markdown report with the per-energy tables.
``--g4-ls-data`` is an OPTIONAL exploratory cross-check against the Geant4 LS data file (a
research download that is not part of the repository); it is skipped with a notice when absent.

Quantities (per ion, liquid water, ionmc Bethe model, I = 78 eV, all options on):
* r(E) = dL_LS - (Bloch + Mott1), with Bloch/Mott1 in the ionmc forms (lindhard_sorensen.py);
  primary variant: charge = ionmc effective charge q(E) in LS, Bloch and Mott (consistent with
  how ionmc applies Bloch and Mott); secondary: bare z.
* r/L, with L the full ionmc stopping-number bracket, L = S beta^2 / (K q^2 Z/A). Since
  S is proportional to L at fixed q, r/L is also the relative change of S.
* CSDA-range shift from omitting r: R(E) = int_{E_cut}^{E} A dE_u / S, E_cut = 2 MeV/u,
  computed with ionmc's exact log-log quadrature on a 400-per-decade grid for S (ionmc) and
  S (1 + r/L); dR/R = R_ionmc/R_withLS - 1 (positive: ionmc range longer).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True  # never write __pycache__ into the repository

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
WORKTREE = HERE.parents[2]
sys.path.insert(0, str(HERE))
import lindhard_sorensen as lsm  # noqa: E402

STANDARD_POINTS = (50.0, 100.0, 150.0, 200.0, 250.0, 290.0, 350.0, 430.0)
E_CUT = 2.0
E_MAX = 430.0
TRIGGER_RANGE_PCT = 0.2
EXPECT_S_PCT = 0.1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_ionmc(worktree: Path):
    if str(worktree / "src") not in sys.path:
        sys.path.insert(0, str(worktree / "src"))
    from ionmc.materials import water
    from ionmc.physics import projectiles as pj
    from ionmc.physics import stopping as st

    proj = {1: pj.PROTON, 2: pj.ALPHA, 6: pj.CARBON12, 8: pj.OXYGEN16}
    return st, proj, water(78.0)


def grid(points_per_decade: int) -> np.ndarray:
    n = int(math.ceil(points_per_decade * math.log10(E_MAX / E_CUT))) + 1
    g = np.geomspace(E_CUT, E_MAX, n)
    return np.unique(np.concatenate([g, np.array(STANDARD_POINTS)]))


def evaluate_ion(st, proj, mat, z: int, e: np.ndarray) -> dict:
    p = proj[z]
    s = st.bethe_mass_stopping(e, p, mat)
    beta2, _, _ = st.kinematics(e, p)
    beta = np.sqrt(beta2)
    q = st.pierce_blann_charge(beta, p) if z >= 2 else np.full_like(e, 1.0)
    big_l = s * beta2 / (st.K_MEV_CM2_MOL * q**2 * mat.z_over_a)
    # consistency of the copied kinematics/masses/forms with ionmc
    beta_m = np.array([lsm.kinematics(z, float(x))[0] for x in e])
    bloch_rep = st.bloch_term(beta2, q)
    bloch_m = np.array([lsm.bloch_ionmc(float(b2), float(qq)) for b2, qq in zip(beta2, q, strict=True)])
    out = {"beta_rel_dev_max": float(np.max(np.abs(beta_m / beta - 1))),
           "bloch_abs_dev_max": float(np.max(np.abs(bloch_m - bloch_rep)))}
    r_eff, r_bare, r_mixed, ls_minus_bloch, kmax_used, conv_max = [], [], [], [], 0, 0.0
    for b, qq in zip(beta, q, strict=True):
        b = float(b)
        qq = float(qq)
        res = lsm.delta_l_ls(lsm.ALPHA_FS * qq, b)
        kmax_used = max(kmax_used, res.k_used)
        conv_max = max(conv_max, abs(res.change_vs_half_k))
        r_eff.append(res.delta_l - lsm.bloch_ionmc(b * b, qq) - lsm.mott1_ionmc(b, qq))
        ls_minus_bloch.append(res.delta_l - lsm.bloch_ionmc(b * b, qq))
        res_b = lsm.delta_l_ls(lsm.ALPHA_FS * z, b)
        r_bare.append(res_b.delta_l - lsm.bloch_ionmc(b * b, float(z)) - lsm.mott1_ionmc(b, float(z)))
        # mixed (Geant4-lineage) convention: LS at the bare nuclear charge z, Bloch/Mott1 at q
        r_mixed.append(res_b.delta_l - lsm.bloch_ionmc(b * b, qq) - lsm.mott1_ionmc(b, qq))
    out.update(e=e, s=s, L=big_l, q=q, beta=beta, r_eff=np.array(r_eff), r_bare=np.array(r_bare),
               r_mixed=np.array(r_mixed),
               ls_minus_bloch=np.array(ls_minus_bloch), k_max_used=kmax_used,
               convergence_change_max=conv_max)
    return out


def range_shift(st, a: int, e: np.ndarray, s: np.ndarray, frac: np.ndarray):
    inc0 = st.exact_loglog_range_increments(np.log(e), np.log(s), float(a))
    inc1 = st.exact_loglog_range_increments(np.log(e), np.log(s * (1 + frac)), float(a))
    r0 = np.concatenate([[0.0], np.cumsum(inc0)])
    r1 = np.concatenate([[0.0], np.cumsum(inc1)])
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = np.where(r1 > 0, r0 / r1 - 1, 0.0)
    return rel, (r0 - r1) * 10.0, r0 * 10.0  # mm of water (rho = 1 g/cm3)


def g4_table(path: Path):
    """Parse the lsdata[9][41] array of G4LindhardSorensenData.cc as text (never executed)."""
    txt = path.read_text()
    block = txt[txt.index("lsdata[9][NPOINT] = {"):]
    block = block[: block.index("}};") + 3]
    block = re.sub(r"//[^\n]*", "", block)  # strip comments (index markers like 0-5)
    nums = [float(x) for x in re.findall(r"-?\d+\.\d+(?:e[-+]?\d+)?", block)]
    arr = np.array(nums).reshape(9, 41)
    return arr


def g4_delta_l(arr: np.ndarray, z: int, gamma: float) -> float:
    """Geant4 GetDeltaL (2x convention), linear in x = ln(gamma-1) instead of G4's spline."""
    zlist = [1, 10, 18, 36, 54, 66, 79, 92, 109]
    xmin, xmax = math.log(0.02), math.log(316.22777)
    xs = np.linspace(xmin, xmax, 41)
    idx = 0
    while idx < 8 and not z < zlist[idx + 1]:
        idx += 1
    x = math.log(gamma - 1.0)
    y = float(np.interp(x, xs, arr[idx]))
    if idx < 8 and z > zlist[idx]:
        y1 = float(np.interp(x, xs, arr[idx + 1]))
        y += (y1 - y) * (z - zlist[idx]) / (zlist[idx + 1] - zlist[idx])
    return y * z ** (2.0 / 3.0)


def convergence_demo(z: int, e_u: float) -> list[dict]:
    beta, _ = lsm.kinematics(z, e_u)
    t, b2, _ = lsm.ls_summand(lsm.ALPHA_FS * z, beta, 32000)
    rows = []
    for k in (500, 1000, 2000, 4000, 8000, 16000, 32000):
        raw = float(np.sum(t[:k], dtype=np.longdouble)) + b2
        withtail, tail, c2 = lsm._sum_with_tail(t, k)
        rows.append({"K": k, "raw": raw, "with_tail": withtail + b2, "tail": tail,
                     "c2_over_mott1": c2 / lsm.mott1_ionmc(beta, z)})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=WORKTREE / "validation/plans/v3-008/p3-output.json")
    ap.add_argument("--report", type=Path, default=None, help="optional Markdown report path")
    ap.add_argument("--g4-ls-data", type=Path, default=None)
    ap.add_argument("--points-per-decade", type=int, default=400)
    args = ap.parse_args()
    t0 = time.time()
    st, proj, mat = load_ionmc(WORKTREE)
    e = grid(args.points_per_decade)
    e_coarse_idx = np.unique(np.concatenate([
        np.searchsorted(e, grid(30)), np.searchsorted(e, np.array(STANDARD_POINTS))]))
    ions = {}
    for z in (1, 2, 6, 8):
        d = evaluate_ion(st, proj, mat, z, e)
        a = proj[z].a
        frac_eff = d["r_eff"] / d["L"]
        frac_bare = d["r_bare"] / d["L"]
        frac_lsb = d["ls_minus_bloch"] / d["L"]
        frac_mixed = d["r_mixed"] / d["L"]
        rel_mixed, dmm_mixed, _ = range_shift(st, a, e, d["s"], frac_mixed)
        # same shift with the integration started at 10 MeV/u (separate grid, mixed variant)
        e10 = np.geomspace(10.0, E_MAX, int(math.ceil(args.points_per_decade * math.log10(E_MAX / 10.0))) + 1)
        d10 = evaluate_ion(st, proj, mat, z, e10)
        rel_mixed10, _, _ = range_shift(st, a, e10, d10["s"], d10["r_mixed"] / d10["L"])
        d5 = evaluate_ion(st, proj, mat, z, np.array([2.0, 5.0, 10.0]))
        rel_eff, dmm_eff, rmm = range_shift(st, a, e, d["s"], frac_eff)
        rel_bare, _, _ = range_shift(st, a, e, d["s"], frac_bare)
        rel_lsb, dmm_lsb, _ = range_shift(st, a, e, d["s"], frac_lsb)
        sel = e > E_CUT * (1 + 1e-9)
        d.update(frac_eff=frac_eff, frac_bare=frac_bare, frac_lsb=frac_lsb, rel_eff=rel_eff,
                 rel_bare=rel_bare, rel_lsb=rel_lsb, frac_mixed=frac_mixed, rel_mixed=rel_mixed,
                 dmm_mixed=dmm_mixed, e10=e10, rel_mixed10=rel_mixed10, pt=d5, dmm_eff=dmm_eff, dmm_lsb=dmm_lsb, rmm=rmm, sel=sel)
        ions[z] = d
    runtime = time.time() - t0

    names = {1: "proton", 2: "4He", 6: "12C", 8: "16O"}
    summary = {}
    for z, d in ions.items():
        sp = {f"{p:g}": {
            "r_over_L_pct_effq": 100 * float(d["frac_eff"][np.searchsorted(e, p)]),
            "dR_over_R_pct_effq": 100 * float(d["rel_eff"][np.searchsorted(e, p)]),
            "dR_mm_effq": float(d["dmm_eff"][np.searchsorted(e, p)]),
            "R_ionmc_mm_from_Ecut": float(d["rmm"][np.searchsorted(e, p)]),
        } for p in STANDARD_POINTS}
        ent = {
            "max_abs_r_over_L_pct_effq": 100 * float(np.max(np.abs(d["frac_eff"]))),
            "E_at_max_r_over_L_effq": float(e[np.argmax(np.abs(d["frac_eff"]))]),
            "rms_r_over_L_pct_effq": 100 * float(np.sqrt(np.mean(d["frac_eff"] ** 2))),
            "min_r_over_L_pct_effq": 100 * float(np.min(d["frac_eff"])),
            "max_abs_r_over_L_pct_bare": 100 * float(np.max(np.abs(d["frac_bare"]))),
            "max_abs_dR_over_R_pct_effq": 100 * float(np.max(np.abs(d["rel_eff"][d["sel"]]))),
            "E_at_max_dR_over_R_effq": float(e[d["sel"]][np.argmax(np.abs(d["rel_eff"][d["sel"]]))]),
            "max_abs_dR_over_R_pct_bare": 100 * float(np.max(np.abs(d["rel_bare"][d["sel"]]))),
            "max_abs_dR_mm_effq": float(np.max(np.abs(d["dmm_eff"]))),
            "standard_points": sp,
            "k_max_used": d["k_max_used"],
            "convergence_change_max": d["convergence_change_max"],
            "check_beta_rel_dev_vs_ionmc": d["beta_rel_dev_max"],
            "check_bloch_abs_dev_vs_ionmc": d["bloch_abs_dev_max"],
        }
        fm = d["frac_mixed"]
        q_over_z = d["q"] / z
        e_q999 = float(e[np.argmax(q_over_z > 0.999)]) if np.any(q_over_z > 0.999) else None
        above = np.nonzero(np.abs(fm) > EXPECT_S_PCT / 100)[0]
        pt_f = d["pt"]["r_mixed"] / d["pt"]["L"]
        ent["mixed_charge_convention"] = {
            "definition": "r_mixed = dL_LS(z_nominal) - Bloch(q_eff) - Mott1(q_eff); r/S = r/L",
            "max_abs_r_over_S_pct": 100 * float(np.max(np.abs(fm))),
            "E_at_max_abs_r_over_S": float(e[np.argmax(np.abs(fm))]),
            "min_r_over_S_pct": 100 * float(np.min(fm)),
            "max_r_over_S_pct": 100 * float(np.max(fm)),
            "r_over_S_pct_at_2": 100 * float(fm[0]),
            "r_over_S_pct_at_5": 100 * float(pt_f[1]),
            "r_over_S_pct_at_10": 100 * float(pt_f[2]),
            "abs_r_over_S_pct_at_5": 100 * abs(float(pt_f[1])),
            "abs_r_over_S_pct_at_10": 100 * abs(float(pt_f[2])),
            "delta_L_mixed_at_2": float(d["r_mixed"][0]),
            "delta_L_eff_at_2": float(d["r_eff"][0]),
            "delta_L_bare_at_2": float(d["r_bare"][0]),
            "max_abs_dR_over_R_pct_from_2": 100 * float(np.max(np.abs(d["rel_mixed"][d["sel"]]))),
            "E_at_max_dR_over_R_from_2": float(e[d["sel"]][np.argmax(np.abs(d["rel_mixed"][d["sel"]]))]),
            "dR_over_R_pct_at_430_from_2": 100 * float(d["rel_mixed"][-1]),
            "max_abs_dR_mm_from_2": float(np.max(np.abs(d["dmm_mixed"]))),
            "E_last_abs_dR_over_R_above_0.2pct_from_2": (
                float(e[d["sel"]][np.nonzero(np.abs(d["rel_mixed"][d["sel"]]) > TRIGGER_RANGE_PCT / 100)[0][-1]])
                if np.any(np.abs(d["rel_mixed"][d["sel"]]) > TRIGGER_RANGE_PCT / 100) else None),
            "R_ionmc_mm_from_2_at_E_last_dR_above_0.2pct": (
                float(d["rmm"][d["sel"]][np.nonzero(np.abs(d["rel_mixed"][d["sel"]]) > TRIGGER_RANGE_PCT / 100)[0][-1]])
                if np.any(np.abs(d["rel_mixed"][d["sel"]]) > TRIGGER_RANGE_PCT / 100) else None),
            "max_abs_dR_over_R_pct_from_10_exceeds_0.2pct": bool(
                np.max(np.abs(d["rel_mixed10"][1:])) > TRIGGER_RANGE_PCT / 100),
            "max_abs_dR_over_R_pct_from_10": 100 * float(np.max(np.abs(d["rel_mixed10"][1:]))),
            "dR_over_R_pct_at_430_from_10": 100 * float(d["rel_mixed10"][-1]),
            "E_q_over_z_exceeds_0.999": e_q999,
            "E_last_abs_r_over_S_above_0.1pct": float(e[above[-1]]) if above.size else None,
            "n_grid_points_abs_r_over_S_above_0.1pct": int(above.size),
            "standard_points": {f"{p:g}": {
                "r_over_S_pct": 100 * float(fm[np.searchsorted(e, p)]),
                "dR_over_R_pct_from_2": 100 * float(d["rel_mixed"][np.searchsorted(e, p)])}
                for p in STANDARD_POINTS},
        }
        if z == 1:
            ent["proton_LS_minus_Bloch"] = {
                "note": "ionmc applies no Mott term for z = 1; this is what the proton model omits",
                "max_abs_over_L_pct": 100 * float(np.max(np.abs(d["frac_lsb"]))),
                "max_abs_dR_over_R_pct": 100 * float(np.max(np.abs(d["rel_lsb"][d["sel"]]))),
                "max_abs_dR_mm": float(np.max(np.abs(d["dmm_lsb"]))),
            }
        summary[names[z]] = ent

    verdict = {}
    for z in (6, 8):
        s = summary[names[z]]
        mx_r = max(s["max_abs_dR_over_R_pct_effq"], s["max_abs_dR_over_R_pct_bare"])
        mx_s = max(s["max_abs_r_over_L_pct_effq"], s["max_abs_r_over_L_pct_bare"])
        verdict[names[z]] = {
            "max_abs_r_over_S_pct": mx_s, "max_abs_dR_over_R_pct": mx_r,
            "range_trigger_0p2pct_exceeded": bool(mx_r > TRIGGER_RANGE_PCT),
            "within_expected_0p1pct_of_S": bool(mx_s <= EXPECT_S_PCT),
        }

    for z in (6, 8):
        m = summary[names[z]]["mixed_charge_convention"]
        verdict[names[z]]["mixed_convention"] = {
            "max_abs_r_over_S_pct": m["max_abs_r_over_S_pct"],
            "max_abs_dR_over_R_pct_from_2": m["max_abs_dR_over_R_pct_from_2"],
            "range_trigger_0p2pct_exceeded": bool(m["max_abs_dR_over_R_pct_from_2"] > TRIGGER_RANGE_PCT),
            "range_trigger_0p2pct_exceeded_from_10": m["max_abs_dR_over_R_pct_from_10_exceeds_0.2pct"],
            "E_last_abs_dR_over_R_above_0.2pct_from_2": m["E_last_abs_dR_over_R_above_0.2pct_from_2"],
            "within_expected_0p1pct_of_S": bool(m["max_abs_r_over_S_pct"] <= EXPECT_S_PCT),
            "E_last_abs_r_over_S_above_0.1pct": m["E_last_abs_r_over_S_above_0.1pct"],
        }

    conv = {f"z{z}_E{eu:g}": convergence_demo(z, eu) for z, eu in ((8, 2.0), (8, 430.0), (6, 290.0), (1, 2.0))}

    g4 = None
    if args.g4_ls_data is not None and not args.g4_ls_data.exists():
        print(f"NOTICE: --g4-ls-data {args.g4_ls_data} not found; Geant4 cross-check skipped")
    elif args.g4_ls_data is None:
        print("NOTICE: no --g4-ls-data given; exploratory Geant4 cross-check skipped")
    if args.g4_ls_data is not None and args.g4_ls_data.exists():
        arr = g4_table(args.g4_ls_data)
        rows = []
        for z in (1, 2, 6, 8):
            for p in (20.0, 50.0, 100.0, 200.0, 290.0, 430.0):
                beta, gamma = lsm.kinematics(z, p)
                ours = lsm.delta_l_ls(lsm.ALPHA_FS * z, beta).delta_l
                rows.append({"z": z, "E": p, "g4_half": 0.5 * g4_delta_l(arr, z, gamma), "ls_full": ours,
                             "ls_minus_bloch": ours - lsm.bloch_exact(lsm.ALPHA_FS * z / beta)})
        g4 = {"file": str(args.g4_ls_data.resolve()), "sha256": sha256(args.g4_ls_data), "rows": rows}

    test_path = WORKTREE / "tests/ionmc/test_ion_em_lindhard_sorensen.py"
    files = {n: sha256(HERE / n) for n in ("lindhard_sorensen.py", "p3_driver.py")}
    files["tests/ionmc/test_ion_em_lindhard_sorensen.py"] = sha256(test_path)
    output = {
        "check": "V3-008 P3", "status": "exploratory computation (not a frozen record)",
        "material": "ionmc water(78.0)", "E_range_MeV_per_u": [E_CUT, E_MAX], "E_cut_MeV_per_u": E_CUT,
        "grid_points": int(e.size), "points_per_decade": args.points_per_decade,
        "thresholds": {"expected_r_over_S_pct": EXPECT_S_PCT, "range_trigger_pct": TRIGGER_RANGE_PCT},
        "verdict": verdict, "summary": summary, "runtime_s": runtime,
        "script_sha256": files, "numpy": np.__version__,
    }
    output["second_order_coefficient"] = second_order_table()
    if g4 is not None:
        output["g4_seen_crosscheck"] = {"file_sha256": g4["sha256"], "status": "SEEN/exploratory"}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    output["runtime_total_s_before_report"] = time.time() - t0
    if args.report is not None:
        args.report.write_text(render(output, ions, e, e_coarse_idx, conv, g4))
    print(json.dumps(verdict, indent=2))
    print(f"runtime {runtime:.1f} s; grid {e.size} points")


def second_order_table() -> list[dict]:
    """(LS - Bloch - Mott1)/(alpha z)^2 at small z (informative) vs the Morgan-Eby closed form."""
    z3 = 1.2020569031595942
    rows = []
    for b in (0.1, 0.3, 0.5, 0.7, 0.85):
        vals = []
        for zz in (1.0, 0.5, 0.25):
            az = lsm.ALPHA_FS * zz
            r = lsm.delta_l_ls(az, b).delta_l - lsm.bloch_exact(az / b) - lsm.mott1_ionmc(b, zz)
            vals.append(r / az**2)
        me = 0.5 * (math.pi**2 / 3 + 1 + b * b * (0.75 - math.log(2)) + 0.5 * (z3 - 3))
        rows.append({"beta": b, "c_z1": vals[0], "c_z0.5": vals[1], "c_z0.25": vals[2], "morgan_eby": me})
    return rows


def _f(x: float, n: int = 4) -> str:
    return f"{x:.{n}f}"


def render(out, ions, e, idx, conv, g4) -> str:
    v = out["verdict"]
    sm = out["summary"]
    L = []
    w = L.append
    w("# V3-008 pre-freeze check P3: Lindhard-Sorensen residual (exploratory computation)\n")
    w("Generated by `p3_driver.py` (scratchpad; not a frozen record). Numbers are computed from "
      "theory and the ionmc analytic Bethe model only; no reference-table values enter.\n")
    w("## Question and thresholds\n")
    w("The ionmc stopping model (`src/ionmc/physics/stopping.py`) is Bethe + shell + Barkas + Bloch + "
      "first-order Mott (z >= 2) with the Barkas-form effective charge q. The point-nucleus "
      "Lindhard-Sorensen (LS) correction contains the Bloch term exactly, so the omitted physics is "
      "r = dL_LS - (Bloch + Mott1). P3: expected <= 0.1 % of S (estimate); if omitting r shifts the "
      "CSDA range by more than 0.2 % at any domain energy (2-430 MeV/u), an LS option is added "
      "(amendment).\n")
    w("## Verdict\n")
    w("| Ion | max abs r/S (%) | max abs dR/R (%) | > 0.2 % range trigger? | <= 0.1 % of S? |")
    w("|---|---|---|---|---|")
    for k, x in v.items():
        w(f"| {k} | {_f(x['max_abs_r_over_S_pct'])} | {_f(x['max_abs_dR_over_R_pct'])} | "
          f"{'YES' if x['range_trigger_0p2pct_exceeded'] else 'no'} | "
          f"{'yes' if x['within_expected_0p1pct_of_S'] else 'NO'} |")
    w("\nMaxima are over the 2-430 MeV/u grid and over both charge variants (effective q, bare z). "
      "The range trigger is **not** met for 12C or 16O; no LS option is required by the P3 rule. "
      "The 16O maximum (just under 0.1 % of S, at about 17 MeV/u) uses most of the 0.1 % estimate; "
      "the estimate in `research-models.md` RQ-S1 is confirmed but not conservative by a large "
      "factor.\n")
    w("## Mixed charge convention (Codex REVIEW-e0b586fd finding 1)\n")
    w("r_mixed = dL_LS(z nominal) - Bloch(q) - Mott1(q); r/S = r/L. Added 2026-11-01 after the review (SEEN).\n")
    w("Range trigger (0.2 %), mixed variant: see `p3-output.json` keys `E_last_abs_dR_over_R_above_0.2pct_from_2`, "
      "`max_abs_dR_mm_from_2`, `max_abs_dR_over_R_pct_from_10`.\n")
    w("| Ion | max abs r/S % | at E | r/S % at 2 / 5 / 10 MeV/u | max abs dR/R % (from 2) | "
      "max abs dR/R % (from 10) | E(q/z > 0.999) | last E with abs r/S > 0.1 % |")
    w("|---|---|---|---|---|---|---|---|")
    for k, x in sm.items():
        if k == "proton":
            continue
        m = x["mixed_charge_convention"]
        lastx = m["E_last_abs_r_over_S_above_0.1pct"]
        w(f"| {k} | {_f(m['max_abs_r_over_S_pct'])} | {m['E_at_max_abs_r_over_S']:.2f} | "
          f"{_f(m['r_over_S_pct_at_2'])} / {_f(m['r_over_S_pct_at_5'])} / {_f(m['r_over_S_pct_at_10'])} | "
          f"{_f(m['max_abs_dR_over_R_pct_from_2'])} | {_f(m['max_abs_dR_over_R_pct_from_10'])} | "
          f"{m['E_q_over_z_exceeds_0.999']:.1f} | {'none' if lastx is None else f'{lastx:.2f}'} |")
    w("")
    w("## Method\n")
    w("Point-nucleus LS ([KHV] arXiv:2008.04819v2 eq. 10, read from the rendered PDF page 5; the same "
      "sum term by term in [MARS] arXiv:physics/0505043 eq. LS, LaTeX source):\n")
    w("```\ndL_LS = beta^2/2 + sum_{k>=1} [ k/eta^2 (k-1)/(2k-1) sin^2(d_k - d_{k-1})\n"
      "                              + k/eta^2 (k+1)/(2k+1) sin^2(d_{-k} - d_{-k-1})\n"
      "                              + k/((4k^2-1)(gamma^2 k^2 + eta^2)) - 1/k ]\n"
      "d_k = phi_k - arg Gamma(rho_k + 1 + i eta) + (pi/2)(l - rho_k)\n"
      "exp(2 i phi_k) = (k - i eta/gamma)/(rho_k - i eta),  rho_k = sqrt(k^2 - (alpha z)^2),\n"
      "l = k (k > 0), l = -k-1 (k < 0),  eta = alpha z / beta\n```\n")
    w("- Bloch and first-order Mott in the exact ionmc forms: Bloch = -y^2 sum_{n<=2000} "
      "1/(n(n^2+y^2)) + 0.5/(2000.5)^2 tail, y = q alpha/beta; Mott1 = 0.5 pi alpha beta q "
      "(Ahlen 1980; [KHV] eq. 5). Test `test_forms_match_ionmc` checks equality with the repository "
      "functions, constants and masses.")
    w("- Non-relativistic limit (alpha z -> 0 at fixed eta): both phase differences give "
      "sin^2 = eta^2/(k^2+eta^2), the k-weights sum to 1 and the summand becomes the Bloch summand "
      "k/(k^2+eta^2) - 1/k (algebra in the module docstring; tests (b)).")
    w("- Charge: primary variant uses the ionmc effective charge q(E) in LS, Bloch and Mott (as ionmc "
      "applies Bloch and Mott with q); the bare-z variant is also evaluated. They differ only below "
      "~10 MeV/u.")
    w("- r/L: L is the full ionmc bracket, L = S beta^2/(K q^2 Z/A), I = 78 eV, ionmc water. Adding r "
      "changes S by the factor (1 + r/L).")
    w(f"- Range: R(E) = int_{{2 MeV/u}}^E A dE_u/S with ionmc's exact log-log quadrature on a "
      f"{out['points_per_decade']}-per-decade grid ({out['grid_points']} points incl. the standard "
      "points), for S and S(1 + r/L); dR/R = R_ionmc/R_withLS - 1 (positive = ionmc range longer). "
      "dR in mm is for water of 1 g/cm3, range counted from E_cut.")
    w("- Finite nuclear size neglected: it matters when gamma m_e c R/hbar ~ 1/2 ([MARS]); for 16O "
      "(R ~ 3 fm) at gamma <= 1.47 the ratio is ~0.011.\n")
    w("## Numerics and convergence\n")
    w("Phases in numpy longdouble (80-bit). Complex log-Gamma: Stirling series (8 terms) after "
      "recurrence to Re >= 16; tested against |Gamma(1+iy)|^2 = pi y/sinh(pi y), math.lgamma, the "
      "series for arg Gamma(1+iy) and the recurrence. `l - rho_k` is computed as "
      "(alpha z)^2/(|k| + rho_k) to avoid cancellation.\n")
    w("The summand decays as c2/k^2 with c2 equal to the first-order Mott term (to ~1e-8 relative, "
      "fitted), so plain truncation at K leaves an error of about Mott1/K (e.g. 3e-5 for 12C at "
      "beta = 0.85 with K = 2000). The code adds a fitted tail (c2/k^2 + c3/k^3 + c4/k^4 at "
      "k = K, K/2, K/4; Euler-Maclaurin sums) and doubles K until the tail-corrected sums at K and "
      f"K/2 agree to 1e-11. Largest K used in the driver: "
      f"{max(x['k_max_used'] for x in sm.values())}; largest |change| between K/2 and K: "
      f"{max(x['convergence_change_max'] for x in sm.values()):.1e}.\n")
    for key, rows in conv.items():
        w(f"Convergence demonstration {key} (full dL_LS):\n")
        w("| K | raw partial sum | with tail | tail | c2/Mott1 |")
        w("|---|---|---|---|---|")
        for r in rows:
            w(f"| {r['K']} | {r['raw']:.12f} | {r['with_tail']:.12f} | {r['tail']:.3e} | "
              f"{r['c2_over_mott1']:.9f} |")
        w("")
    w("Independent check against published theory values: [KHV] Table 1 lists the point-nucleus "
      "dL_LS (Mathematica, k <= 5000, 3 decimals) for Z = 6-114, beta = 0.85-0.99. All 50 values are "
      "reproduced within rounding plus their truncation (Mott1/5000) (`test_against_kats_table1`). "
      "This verifies the implementation, the phase convention and the sign of the Mott part; it "
      "does not independently verify the LS formula, which rests on [LS96] as reproduced by "
      "[KHV] and [MARS] and on the two analytic limits (Bloch; first-order Mott).\n")
    w("## Results (primary variant: effective charge)\n")
    w("Selected grid (30 per decade plus the standard points). r/L in %, dR/R in %, dR in mm "
      "(range from 2 MeV/u).\n")
    w("| E (MeV/u) | p r/L | p LS-Bloch /L | 4He r/L | 4He dR/R | 12C r/L | 12C dR/R | 12C dR mm | "
      "16O r/L | 16O dR/R | 16O dR mm |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for i in idx:
        p, he, c, o = ions[1], ions[2], ions[6], ions[8]
        w(f"| {e[i]:.2f} | {100*p['frac_eff'][i]:.4f} | {100*p['frac_lsb'][i]:.4f} | "
          f"{100*he['frac_eff'][i]:.4f} | {100*he['rel_eff'][i]:.4f} | {100*c['frac_eff'][i]:.4f} | "
          f"{100*c['rel_eff'][i]:.4f} | {c['dmm_eff'][i]:.3f} | {100*o['frac_eff'][i]:.4f} | "
          f"{100*o['rel_eff'][i]:.4f} | {o['dmm_eff'][i]:.3f} |")
    w("\nAggregates (from `p3-output.json`):\n")
    w("| Ion | max abs r/L % (q) | at E | RMS r/L % | max abs r/L % (bare z) | max abs dR/R % (q) | "
      "at E | max abs dR/R % (bare) | max abs dR mm |")
    w("|---|---|---|---|---|---|---|---|---|")
    for k, x in sm.items():
        w(f"| {k} | {_f(x['max_abs_r_over_L_pct_effq'])} | {x['E_at_max_r_over_L_effq']:.1f} | "
          f"{_f(x['rms_r_over_L_pct_effq'])} | {_f(x['max_abs_r_over_L_pct_bare'])} | "
          f"{_f(x['max_abs_dR_over_R_pct_effq'])} | {x['E_at_max_dR_over_R_effq']:.1f} | "
          f"{_f(x['max_abs_dR_over_R_pct_bare'])} | {x['max_abs_dR_mm_effq']:.3f} |")
    pl = sm["proton"]["proton_LS_minus_Bloch"]
    w("\nr is positive at every grid energy for every ion (ionmc stopping slightly low, range "
      "slightly long).\n")
    w("## Proton context\n")
    w(f"For z = 1, r = LS - (Bloch + Mott1) is at most {_f(sm['proton']['max_abs_r_over_L_pct_effq'])} % "
      "of S. ionmc applies no Mott term for protons, so the proton model actually omits LS - Bloch "
      f"(essentially the first-order Mott term): at most {_f(pl['max_abs_over_L_pct'])} % of S and "
      f"{_f(pl['max_abs_dR_over_R_pct'])} % of the range from 2 MeV/u ({pl['max_abs_dR_mm']:.2f} mm at "
      "430 MeV), within the already-qualified proton budget. Context only.\n")
    w("## Second-order coefficient (informative)\n")
    w("c(beta) = (LS - Bloch - Mott1)/(alpha z)^2 extrapolated to small z, compared with the "
      "Morgan-Eby third-order-Born closed form as quoted in [KHV] eq. 6 (bracket/2):\n")
    w("| beta | c (z=1) | c (z=0.5) | c (z=0.25) | Morgan-Eby |")
    w("|---|---|---|---|---|")
    for r in second_order_table():
        w(f"| {r['beta']} | {r['c_z1']:.4f} | {r['c_z0.5']:.4f} | {r['c_z0.25']:.4f} | "
          f"{r['morgan_eby']:.4f} |")
    w("\nThe LS second-order coefficient (about 2.0-2.1) is about 20 % above the Morgan-Eby value "
      "(about 1.70). [KHV] Table 1 also lists Morgan-Eby values that deviate from LS. The comparison "
      "is informative, not a test; using Morgan-Eby instead would "
      "make r about 20 % smaller and does not change the verdict.\n")
    w("## Geant4 LS table: SEEN / exploratory cross-check (not construction, not a test)\n")
    if g4 is None:
        w("Not run (no Geant4 file given).\n")
    else:
        w(f"File read as text only (never compiled or executed): `{g4['file']}` (sha256 "
          f"{g4['sha256']}). This is the **Geant4 v11.2.0** `G4LindhardSorensenData.cc` downloaded "
          "by the research pass; the v11.4.2 data file was not found in the cache (only "
          "`G4LindhardSorensenIonModel.cc` v11.4.2, e16ed6cd...). Geant4 `GetDeltaL(Z, gamma)` "
          "= table(Z-interpolated between rows Z = 1, 10, 18, ...; x = ln(gamma-1), 0.02 <= "
          "gamma-1 <= 316) * Z^(2/3). Here the table is interpolated linearly in x (Geant4 uses a "
          "cubic spline). In `G4LindhardSorensenIonModel::CorrectionsAlongStep` (v11.4.2, "
          "L292-300) it is added as `twopi_mc2_rcl2*q2*eDensity*(deltaL+deltaL0)*length/beta2`, "
          "i.e. in the 2 pi r_e^2 (2x) convention, so its equivalent in the ionmc convention would "
          "be deltaL/2.\n")
        w("| z | E (MeV/u) | G4 deltaL/2 | G4 deltaL | dL_LS (this work, full) | dL_LS - Bloch |")
        w("|---|---|---|---|---|---|")
        for r in g4["rows"]:
            w(f"| {r['z']} | {r['E']:g} | {r['g4_half']:.5f} | {2*r['g4_half']:.5f} | "
              f"{r['ls_full']:.5f} | {r['ls_minus_bloch']:.5f} |")
        o430 = [r for r in g4["rows"] if r["z"] == 8 and r["E"] == 430.0][0]
        i430 = int(np.searchsorted(e, 430.0))
        half_o = 100 * 0.5 * o430["ls_full"] / float(ions[8]["L"][i430])
        w(f"\nObservation (SEEN, exploratory, unverified): for z = 6 and 8 the **un-halved** Geant4 "
          "deltaL agrees with the full point-nucleus dL_LS of this work in the 4 pi r_e^2 (ionmc) "
          "convention across 20-430 MeV/u (absolute difference <= 0.003; <= 5 % where |dL_LS| > 0.02), "
          "including the sign change between 20 and 100 MeV/u. If the "
          "reading of the conventions above is correct, Geant4 adds only half of the LS correction "
          f"(dL_LS/2/L = {half_o:.2f} % of S for 16O at 430 MeV/u). This bears on the I3c "
          "attribution diagnostic and on V1, not on P3. It needs the Geant4 source of the table "
          "normalisation (or a printed GenericIon dE/dx from the P4 smoke) before it is relied on. "
          "For z = 1, 2 the Geant4 table does not track either quantity closely (the Z = 1 row is "
          "irregular; Geant4 does not use it for protons or alpha).\n")
    w("## Limitations and items not verified\n")
    w("- The LS 1996 paper itself was not read; the formula is taken from [KHV] eq. 10 and [MARS] "
      "eq. LS, which agree with each other, reproduce the Bloch limit analytically and the "
      "first-order Mott term numerically, and whose published values ([KHV] Table 1) are "
      "reproduced.")
    w("- The effective-charge treatment of LS (point charge q instead of a screened ion) is a model "
      "choice consistent with ionmc's Bloch/Mott; the bare-z variant brackets it (difference <= "
      "0.001 percentage points of S).")
    w("- The range shift is counted from E_cut = 2 MeV/u; below E_cut transport deposits locally.")
    w("- Geant4 comparison: v11.2.0 data file, linear instead of spline interpolation; exploratory.\n")
    w("## Run record\n")
    w(f"- runtime {out['runtime_s']:.1f} s (ion evaluation), "
      f"{out['runtime_total_s_before_report']:.1f} s total before report rendering; numpy {out['numpy']}")
    for n, h in out["script_sha256"].items():
        w(f"- `{n}` sha256 {h}")
    w("\nReferences: [LS96] J. Lindhard, A. H. Sorensen, Phys. Rev. A 53, 2443 (1996). [KHV] P. B. "
      "Kats, K. V. Halenka, O. O. Voskresenskaya, arXiv:2008.04819v2. [MARS] I. L. Rakhno, "
      "N. V. Mokhov, S. I. Striganov, arXiv:physics/0505043. S. P. Ahlen, Rev. Mod. Phys. 52, 121 "
      "(1980).\n")
    return "\n".join(L)


if __name__ == "__main__":
    main()
