"""D-E8 (exploratory, report-only, never evidence): the go/no-go screen of Amendment 20 (c) as carried by Amendment 21 (e) for the model under
test M2' (V-M2' = ``m2prime`` arm); the same quantities for the report-only arms ``m2`` and ``rb`` (V-RB).

Reads the D-E3 harness documents (``de3-<E>-<variant>.json``, CRN seeds 20990100 / 20990101), the D-E5 variant records, the committed de4 ratios,
the committed V11 / V5 archive (``validation/results/transport/lv5c-dd976c5-single-process.json``, steps 16 and 21), the committed de7 record and
the committed de6 report. Nothing is run here and there is no re-run path: a missing document or a failed binding is a SystemExit; a run that
published an invalid document is recorded as executed with verdict STOP ("not established").

Quantities (Amendment 20 (c)1), computed exactly as quoted. ``F = (peak/plateau)_EM+el / (peak/plateau)_EM`` and ``gain = plateau_EM+el /
plateau_EM - 1`` with the plateau the 20-60 mm mean of the 0.5 mm IDD (``compare_idd_v5.PLATEAU_MM``); the EM-only run is shared.
  D1  c_E = (F_arm - F_base) / Delta_obs,  Delta_obs = F_TOPAS,arch - F_ionmc,arch (V11, step 16)
  D1b q_E = PP_ionmc,arch (F_arm / F_base) / PP_TOPAS,arch - 1 (V5, step 21; multiplicative composition, Amendment 19 (b)1)
  D2  d_g = g_ionmc,arch + (g_arm - g_base) - g_TOPAS,arch (additive)
  D3  p_E = (I_ionmc,arch + delta_E) / I_TOPAS,arch - 1, delta_E = (sum over 20 <= z <= 60 mm of IDD_arm - IDD_base) x bin_g_cm2 [MeV/primary] =
      the change of ``compare_idd_v5.curve_metrics`` ``plateau_integral_mev`` of the elastic-only runs (additive composition)
  D4  exceedance ratio (arm / table) from the D-E5 sampler x (table / Geant4) from de4 = P_arm / P_G4 at theta_2mm, water cases 0-4, p+A class
  D5  tail-moment and theta_lab median ratios (arm / table); D6 the D-E0 kernels K* and G of the arm table against the pinned table; D7 recoil share.
The paired (common-random-number) delete-one-batch jackknife SE is reported for D1-D3 and does NOT enter any rule.

CONTINUE iff (i), (i-b), (ii), (iii) hold at both energies for the M2' arm and (v) the committed de7 construction record is in order (triggers 2-4
and the builder-vs-scratch 2 % rule; referenced by sha256) and all M2' documents are valid. Condition (iv) is REPORT-ONLY (Amendment 21 review A21-3):
recorded at all five energies with the 30 / 50 MeV cases flagged "identical by construction" for tables that keep the pinned black disk below 80 MeV,
and the in-domain form (100 / 150 / 200 MeV, 3 of 3); it is never part of the verdict. Every band is recomputed from ``compare_idd_v5.t95`` (floored
Welch df) and the archived SE / df; the conditions use the exact values (the printed bands are their roundings).

Run (one execution, after the harness documents exist): ``python -m validation.scripts.elastic_attribution.de8_diagnostic --reviewed-head-review-id
REVIEW-xxxxxxxx``. The evaluator refuses a dirty tree; HEAD, the clean flag and the sha256 of the code files are recorded, and the review id is bound
to HEAD in ``de8-manifest.json``.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import numpy as np

from . import common as c

CRN_SEEDS = {150: 20990100, 200: 20990101}
ENERGIES = (150, 200)
ARMS = ("m2prime", "m2", "rb")
GATED_ARM = "m2prime"
D4_ENERGIES = (30, 50, 100, 150, 200)
WATER_CASES = {30: 0, 50: 1, 100: 2, 150: 3, 200: 4}
DOMAIN_EMIN_MEV = 80.0  # M2' / R-B keep the pinned black disk below this node energy (Amendment 21 (c)7)
IDENTICAL_BELOW = {"m2prime": DOMAIN_EMIN_MEV, "rb": DOMAIN_EMIN_MEV}
D4_SEED_BASES = {"m2prime": 20991300, "m2": 20991200}
TOL = {"F": 0.02, "PP": 0.02, "gain": 0.005, "I": 0.10}
N_HIST, N_BATCH = 100_000, 20
ARCHIVE_JSON = c.REPO / "validation" / "results" / "transport" / "lv5c-dd976c5-single-process.json"
RES = c.REPO / "validation" / "results" / "elastic_attribution"
H5_JSON = RES / "de0h5" / "h5.json"
DE7_JSON = RES / "de7" / "de7.json"
DE6_REPORT = RES / "de6" / "de6-report.md"
DE4_RATIOS = RES / "de4" / "de4-ratios.json"
DE5_COMMITTED = RES / "de5" / "de5.json"
REVIEW_ID = re.compile(r"^REVIEW-[0-9a-f]{8}$")
DISCLOSURE = (
    "Go/no-go screen, Amendment 20 (c) as carried by Amendment 21 (e) and its addendum (a); report-only, never evidence. It uses the V11 reference "
    "(archived TOPAS F and gain), the V5 reference (TOPAS ΔIDD and peak/plateau) and the Geant4 thin-target ratios to decide whether the single "
    "pre-committed model M2′ is transported at all. It is not a selection among models: M2′ has no free parameter, plan §4.2 forbids switching to M3 "
    "on such outcomes, and STOP leads to \"no model change\". It is a stopping rule on gating quantities. The lv5d V11 and V5 runs are therefore not "
    "blind to these estimates; they replicate them on independent seeds with the frozen comparators. The M2′ construction outcomes were known before "
    "registration (Amendment 21 (a)3) and the de7 record before this addendum; the screen's transport quantities were not. Condition (iv) was revised "
    "after de7 to 100 / 150 / 200 MeV and made report-only: its outcome there is predictable from the disclosed construction-level P values. The decisive "
    "blind conditions are (i), (i-b), (ii) and (iii). Executed once; recorded whichever way it falls."
)  # verbatim from the acceptance plan, Addendum (a) to Amendment 21, item 5 (a test compares it with the plan text)

GIT_ALLOW = (
    "disclosure", "verdict", "established", "failing_conditions/*", "iv_report_only_note", "conditions/**",
    "bands/*/*", "archived_se_df/*/*", "arms/**", "d4/**", "d5/**", "d6/**", "d7/**", "runs/**", "provenance/**",
    "leak_check/*",
)  # fmt: skip
GIT_FILES = {"de8.json": (GIT_ALLOW, ())}


# ---------------------------------------------------------------------------------------------------------------------------------
# archived values and bands (recomputed from the frozen rule)
# ---------------------------------------------------------------------------------------------------------------------------------
def _cmp() -> Any:
    from .de3_probes import _cmp as cmp_module

    return cmp_module()


def load_archive(path: Path = ARCHIVE_JSON) -> dict[int, dict[str, float]]:
    """Archived V11 (step 16) and V5 (step 21, TOPAS engine) values with SE and df, per energy; read from the committed single-process record."""
    d = json.loads(Path(path).read_text())["steps"]
    v11 = d["16-v11-compare"]["verdict"]["energies"]
    v5 = d["21-v5-compare"]["verdict"]["engines"]["topas"]["energies"]
    out: dict[int, dict[str, float]] = {}
    for e in ENERGIES:
        f, g = v11[str(e)]["F"], v11[str(e)]["plateau_gain"]
        pp, ii = v5[str(e)]["peak_over_plateau"], v5[str(e)]["delta_idd_plateau_integral"]
        out[e] = {
            "F_ion": f["ionmc"]["value"], "F_top": f["reference"]["value"], "F_se": f["se"], "F_df": f["df"],
            "g_ion": g["ionmc"]["value"], "g_top": g["reference"]["value"], "g_se": g["se"], "g_df": g["df"],
            "PP_ion": pp["ionmc"]["value"], "PP_top": pp["reference"]["value"], "PP_se": pp["se"], "PP_df": pp["df"],
            "I_ion": ii["ionmc"]["value"], "I_top": ii["reference"]["value"], "I_se": ii["se"], "I_df": ii["df"],
        }  # fmt: skip
    return out


def crosscheck_h5(arch: dict[int, dict[str, float]], h5_path: Path | None = None) -> None:
    """The archive's V11 values equal the committed H5 ``archived_v11`` block (to 1e-12); fail closed."""
    p = Path(h5_path or H5_JSON)
    blk = json.loads(p.read_text())["archived_v11"]["energies"]
    for e in ENERGIES:
        for a, b in (("F_ion", "F_ionmc"), ("F_top", "F_topas"), ("g_ion", "gain_ionmc"), ("g_top", "gain_topas"), ("F_se", "F_relative_diff_se")):
            if abs(arch[e][a] - blk[str(e)][b]) > 1e-12:
                raise SystemExit(f"archived {a} at {e} MeV differs from the committed H5 record ({b})")


def allowed(cmp: Any, tol: float, se: float, df: float) -> tuple[float, float]:
    """``(allowed |estimate|, t95)``: the point estimates whose 90 % interval ``d +- t95 SE`` would lie inside the tolerance at the archived SE."""
    t = cmp.t95(df)
    return tol - t * se, t


def compute_bands(cmp: Any, arch: dict[int, dict[str, float]], f_base: dict[int, float]) -> dict[int, dict[str, float]]:
    """Every clause band recomputed from ``t95`` (nothing printed is hard-coded). Closure bands: V11 F and V5 peak/plateau in units of Delta_obs."""
    out: dict[int, dict[str, float]] = {}
    for e in ENERGIES:
        a = arch[e]
        aF, tF = allowed(cmp, TOL["F"], a["F_se"], a["F_df"])
        aP, tP = allowed(cmp, TOL["PP"], a["PP_se"], a["PP_df"])
        aG, tG = allowed(cmp, TOL["gain"], a["g_se"], a["g_df"])
        aI, tI = allowed(cmp, TOL["I"], a["I_se"], a["I_df"])
        delta = a["F_top"] - a["F_ion"]
        r = a["PP_top"] / a["PP_ion"]
        out[e] = {
            "delta_obs": delta, "a_F": aF, "t95_F": tF, "a_PP": aP, "t95_PP": tP, "a_gain": aG, "t95_gain": tG, "a_I": aI, "t95_I": tI,
            "F_closure_lo": (a["F_top"] * (1 - aF) - a["F_ion"]) / delta, "F_closure_hi": (a["F_top"] * (1 + aF) - a["F_ion"]) / delta,
            "PP_closure_lo": ((1 - aP) * r - 1.0) * f_base[e] / delta, "PP_closure_hi": ((1 + aP) * r - 1.0) * f_base[e] / delta,
        }  # fmt: skip
    return out


# ---------------------------------------------------------------------------------------------------------------------------------
# run documents
# ---------------------------------------------------------------------------------------------------------------------------------
def sha256_bytes_of(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reused_hashes(report: Path = DE6_REPORT) -> dict[str, str]:
    """``{file name: sha256}`` of the input list of the committed de6 report (the base and EM-only documents are reused only against these)."""
    out = {}
    for m in re.finditer(r"^- `([^`]+)`\s+([0-9a-f]{64})\s*$", Path(report).read_text(), re.M):
        out[Path(m.group(1)).name] = m.group(2)
    return out


def load_reused(de3_dir: Path, report: Path = DE6_REPORT) -> dict[tuple[int, str], dict[str, Any]]:
    """The four base / EM-only documents, each bound by sha256 to ``de6-report.md``. A missing file or a mismatch is a SystemExit (no re-run here)."""
    want = reused_hashes(report)
    out = {}
    for e in ENERGIES:
        for v in ("all", "emonly"):
            p = Path(de3_dir) / f"de3-{e}-{v}.json"
            if not p.is_file():
                raise SystemExit(f"{p} is missing: the base / EM-only document is required (the evaluator never runs anything)")
            if want.get(p.name) != sha256_bytes_of(p):
                raise SystemExit(f"{p.name}: sha256 differs from the input list of {Path(report).name}")
            out[(e, v)] = json.loads(p.read_text())
    return out


def valid_document(doc: dict[str, Any], energy: int, mode: str, *, arm: bool) -> bool:
    """Frozen run conditions; for arm documents also the elastic-only configuration, the table binding and the recorded code provenance."""
    ok = (
        doc.get("n") == N_HIST and doc.get("n_batches") == N_BATCH and doc.get("seed") == CRN_SEEDS[energy] and doc.get("energy") == energy
        and doc.get("backend") == "warp-cpu" and doc.get("precision") == "float64" and doc.get("mode") == mode
        and bool(doc.get("valid")) and doc.get("relative_residual", float("inf")) <= 1e-12
    )  # fmt: skip
    if arm:
        ok = ok and doc.get("configuration") == "elastic_only" and bool(doc.get("elastic_table_id")) and bool(doc.get("elastic_table_npz_sha256"))
        ok = ok and "recoil_idd_batches" in doc
    return bool(ok)


REGISTERED_HARNESS_SEEDS = frozenset(CRN_SEEDS.values())
DE7_VARIANT = {"m2prime": "M2prime", "m2": "M2", "rb": "RB"}


def check_seeds(doc: dict[str, Any], what: str, *, d4_arm: str | None = None) -> None:
    """Refuse a seed outside the registered set: harness documents 20990100 / 20990101; a D-E5 record of an arm only its D4 block
    (20991300-20991309 for M2', 20991200-20991209 for M2): the seed base and every case seed."""
    if d4_arm is None:
        if doc.get("seed") not in REGISTERED_HARNESS_SEEDS:
            raise SystemExit(f"{what}: seed {doc.get('seed')} is not a registered harness seed {sorted(REGISTERED_HARNESS_SEEDS)}")
        return
    base = D4_SEED_BASES[d4_arm]
    seeds = [cs.get("seed") for cs in doc.get("cases", {}).values()]
    if doc.get("seed_base") != base or not all(isinstance(x, int) for x in seeds) or sorted(seeds) != list(range(base, base + 10)):
        raise SystemExit(f"{what}: D4 seeds are not the registered block {base}-{base + 9} of the {d4_arm} arm")


def check_expected_tables(expected: dict[str, dict[str, str]], de7_doc: dict[str, Any]) -> None:
    """The expected table identities are an explicit manifest input checked against the committed de7 record (M2' and M2 are recorded there; the
    R-B table is taken from the manifest as given)."""
    for arm, e in expected.items():
        if not e.get("table_id") or not e.get("npz_sha256"):
            raise SystemExit(f"expected table of the {arm} arm needs table_id and npz_sha256")
        if arm in DE7_VARIANT:
            rec = de7_doc["variants"][DE7_VARIANT[arm]]["identity"]
            if (e["table_id"], e["npz_sha256"]) != (rec.get("table_id"), rec.get("npz_sha256")):
                raise SystemExit(f"expected {arm} table differs from the committed de7 record")


def check_table(doc: dict[str, Any], expected: dict[str, str] | None, arm: str, what: str) -> None:
    """The run document's elastic table id and npz sha256 must equal the expected (de7-checked) table of the arm."""
    if expected is None:
        raise SystemExit(f"no expected table for the {arm} arm (manifest input required)")
    if (doc.get("elastic_table_id") or doc.get("table_id")) != expected["table_id"] or doc.get("elastic_table_npz_sha256") != expected["npz_sha256"]:
        raise SystemExit(f"{what}: elastic table id / npz sha256 differ from the expected {arm} table")


def check_code(doc: dict[str, Any], reviewed_sha: str, scripts: dict[str, str], what: str, ancestor: Any = None, head: str | None = None) -> None:
    """A run document must record a clean tree, the evaluator's HEAD (every document of the screen was produced at one commit), a HEAD that contains
    the reviewed commit, and the code hashes of that commit (SystemExit otherwise)."""
    cd = doc.get("code") or {}
    anc = ancestor or is_ancestor
    if head is not None and cd.get("git_head") != head:
        raise SystemExit(f"{what}: recorded HEAD {str(cd.get('git_head'))[:12]} differs from the evaluator HEAD {head[:12]} (all documents must come from one commit)")
    if cd.get("git_clean") is not True or not isinstance(cd.get("git_head"), str) or not anc(reviewed_sha, cd["git_head"]):
        raise SystemExit(f"{what}: recorded HEAD does not contain the reviewed commit {reviewed_sha[:12]}, or the tree was not clean (clean tree required)")
    for n, s in scripts.items():
        if (cd.get("scripts_sha256") or {}).get(n) != s:
            raise SystemExit(f"{what}: sha256 of {n} recorded in the document differs from the current file")


# ---------------------------------------------------------------------------------------------------------------------------------
# D1-D3, D5-D7
# ---------------------------------------------------------------------------------------------------------------------------------
def _jk(stat: Any, arrays: list[np.ndarray], den: np.ndarray) -> tuple[float, float]:
    from .de3_probes import _jackknife

    return _jackknife(stat, arrays, den)


def arm_quantities(cmp: Any, arch: dict[str, float], delta_obs: float, arm: np.ndarray, base: np.ndarray, em: np.ndarray) -> dict[str, Any]:
    """D1, D1b, D2, D3 of one arm at one energy with their paired-jackknife SE (reported, never used in a rule). ``arm``, ``base``, ``em`` are the
    ``[batches, nz]`` IDD arrays of the elastic-only arm run, the elastic-only base run and the shared EM-only run."""
    m = cmp.curve_metrics
    pop = lambda x: m(x)["peak_over_plateau"]  # noqa: E731
    f = lambda n, d: pop(n) / pop(d)  # noqa: E731
    gain = lambda n, d: m(n)["plateau"] / m(d)["plateau"] - 1.0  # noqa: E731
    integ = lambda n: m(n)["plateau_integral_mev"]  # noqa: E731
    f_base, _ = _jk(lambda n, d: f(n[0], d), [base], em)
    f_arm, _ = _jk(lambda n, d: f(n[0], d), [arm], em)
    c1, c1se = _jk(lambda n, d: (f(n[0], d) - f(n[1], d)) / delta_obs, [arm, base], em)
    q, qse = _jk(lambda n, d: arch["PP_ion"] * (f(n[0], d) / f(n[1], d)) / arch["PP_top"] - 1.0, [arm, base], em)
    dg, dgse = _jk(lambda n, d: arch["g_ion"] + (gain(n[0], d) - gain(n[1], d)) - arch["g_top"], [arm, base], em)
    dl, dlse = _jk(lambda n, d: integ(n[0]) - integ(n[1]), [arm, base], em)
    p = (arch["I_ion"] + dl) / arch["I_top"] - 1.0
    return {
        "F_arm": f_arm, "F_base": f_base, "D1_closure": c1, "D1_se": c1se, "D1b_q": q, "D1b_se": qse, "D2_dg": dg, "D2_se": dgse,
        "delta_plateau_integral": dl, "delta_se": dlse, "D3_p": p, "D3_se": dlse / arch["I_top"],
        "gain_arm": gain(arm.mean(axis=0), em.mean(axis=0)), "gain_base": gain(base.mean(axis=0), em.mean(axis=0)),
    }  # fmt: skip


def conditions(q: dict[str, Any], b: dict[str, float]) -> dict[str, dict[str, Any]]:
    """(i), (i-b), (ii), (iii) of one energy against the recomputed bands."""
    return {
        "i": {"value": q["D1_closure"], "lo": b["F_closure_lo"], "hi": b["F_closure_hi"], "holds": b["F_closure_lo"] <= q["D1_closure"] <= b["F_closure_hi"]},
        "i_b": {"value": q["D1b_q"], "allowed": b["a_PP"], "holds": abs(q["D1b_q"]) <= b["a_PP"]},
        "ii": {"value": q["D2_dg"], "allowed": b["a_gain"], "holds": abs(q["D2_dg"]) <= b["a_gain"]},
        "iii": {"value": q["D3_p"], "allowed": b["a_I"], "holds": abs(q["D3_p"]) <= b["a_I"]},
    }  # fmt: skip


def d4_report(arm: str, de5_arm: dict[str, Any], de5_table: dict[str, Any], de4: dict[str, Any]) -> dict[str, Any]:
    """D4 per energy (water cases 0-4, p+A class O-16, P(theta_lab > theta_2mm)): ratio arm / Geant4 = (arm / table) x (table / Geant4), and (iv) as
    REPORT-ONLY: |ratio_arm/G4 - 1| < |ratio_table/G4 - 1|, at all five energies and in the in-domain form (100 / 150 / 200 MeV), with the energies
    below the domain boundary flagged identical by construction. Also D5: tail-moment and theta_lab median ratios (arm / table)."""
    rows: dict[str, Any] = {}
    for e in D4_ENERGIES:
        k = str(WATER_CASES[e])
        a, t = de5_arm["cases"][k]["classes"]["O-16"], de5_table["cases"][k]["classes"]["O-16"]
        tg = de4["cases"][k]["O-16"]["signatures"]["P_theta_lab_gt_theta_2mm"]["ratio"]
        arm_over_table = a["P_theta_lab_gt_theta_2mm"] / t["P_theta_lab_gt_theta_2mm"]
        r = arm_over_table * tg
        ident = e < IDENTICAL_BELOW.get(arm, 0.0)
        rows[str(e)] = {
            "arm_over_table": arm_over_table, "table_over_geant4": tg, "arm_over_geant4": r, "moves_toward_1": bool(abs(r - 1.0) < abs(tg - 1.0)),
            "identical_by_construction": bool(ident),
            "D5_tail_moment_ratio": a["T_tail_moment"] / t["T_tail_moment"], "D5_theta_lab_median_ratio": a["theta_lab_median_deg"] / t["theta_lab_median_deg"],
        }  # fmt: skip
    toward = [v["moves_toward_1"] for v in rows.values()]
    dom = [rows[str(e)]["moves_toward_1"] for e in (100, 150, 200)]
    return {"rows": rows, "iv_toward_count_5": int(sum(toward)), "iv_5_of_5_form_ge_4": bool(sum(toward) >= 4),
            "iv_in_domain_3_of_3": bool(all(dom)), "report_only": True}  # fmt: skip


def d6_kernels(arm_rec: dict[str, Any], table_rec: dict[str, Any]) -> dict[str, Any]:
    """D6: ratios (arm / pinned table) of the numeric D-E0 kernel entries K* (O-16, table shape) and G at 150 / 200 MeV (report-only)."""
    out: dict[str, Any] = {}
    for e in ENERGIES:
        out[str(e)] = {}
        for grp in ("O-16", "G"):
            a = arm_rec["kernels"][str(e)]["table"].get(grp, {})
            t = table_rec["kernels"][str(e)]["table"].get(grp, {})
            out[str(e)][grp] = {k: a[k] / t[k] for k in sorted(a) if k in t and isinstance(a[k], (int, float)) and isinstance(t[k], (int, float)) and t[k]}
    return out


def d7_recoil(arm_doc: dict[str, Any], base_doc: dict[str, Any], em: np.ndarray, gain_excess: float) -> dict[str, Any]:
    """D7: the recoil channel's plateau contribution of the arm against the base and its share of the observed gain excess (H4 form of de6_swap)."""
    from .de6_swap import recoil_h4

    den = float(em.mean(axis=0)[40:121].mean())
    h = recoil_h4({"all": np.asarray(base_doc["recoil_idd_batches"]), "arm": np.asarray(arm_doc["recoil_idd_batches"])}, den, gain_excess)
    return h["arm"]  # type: ignore[no-any-return]


def de7_condition(de7_path: Path = DE7_JSON) -> dict[str, Any]:
    """(v): the committed de7 record of the M2' construction (triggers 2-4 and the builder-vs-scratch 2 % rule), referenced by sha256."""
    p = Path(de7_path)
    if not p.is_file():
        raise SystemExit(f"{p} is missing: condition (v) refers to the committed de7 record")
    v = json.loads(p.read_text())
    m = v["variants"]["M2prime"]
    flags = {"trigger2": bool(m["trigger2"]["passed"]), "trigger3": bool(m["trigger3"]["passes"]), "trigger4": bool(m["trigger4"]["passes"]),
             "reproduction_2pct": bool(str(v["reproduction"]["rule_2pct_verdict"]).startswith("all six cells within 2 %"))}  # fmt: skip
    return {**flags, "holds": all(flags.values()), "de7_sha256": sha256_bytes_of(p), "table_id": m["identity"]["table_id"]}


# ---------------------------------------------------------------------------------------------------------------------------------
# the evaluator
# ---------------------------------------------------------------------------------------------------------------------------------
SCRIPT_REL = "validation/scripts/elastic_attribution"
SHA40 = re.compile(r"^[0-9a-f]{40}$")


def _git(repo: Path | None, *args: str) -> Any:
    import subprocess

    return subprocess.run(["git", *args], cwd=Path(repo or c.REPO), capture_output=True, text=True)


def is_ancestor(sha: str, head: str, repo: Path | None = None) -> bool:
    """``git merge-base --is-ancestor sha head`` (true when equal); an unknown commit is a SystemExit."""
    r = _git(repo, "merge-base", "--is-ancestor", sha, head)
    if r.returncode not in (0, 1):
        raise SystemExit(f"git merge-base --is-ancestor {sha[:12]} {head[:12]} failed: {r.stderr.strip()}")
    return bool(r.returncode == 0)


def hashes_at(sha: str, repo: Path | None = None) -> dict[str, str]:
    """sha256 of the four code files as committed in ``sha`` (``git show <sha>:<path>``, not the working tree)."""
    out = {}
    for n in c.CODE_FILES:
        r = _git(repo, "show", f"{sha}:{SCRIPT_REL}/{n}")
        if r.returncode != 0:
            raise SystemExit(f"{n} is not in the reviewed commit {sha[:12]}: {r.stderr.strip()}")
        out[n] = hashlib.sha256(r.stdout.encode("utf-8")).hexdigest()
    return out


def working_hashes(repo: Path | None = None) -> dict[str, str]:
    root = Path(repo or c.REPO) / SCRIPT_REL
    for n in c.CODE_FILES:
        if not (root / n).is_file():
            raise SystemExit(f"{root / n} is missing")
    return {n: hashlib.sha256((root / n).read_bytes()).hexdigest() for n in c.CODE_FILES}


def check_review(review_id: str | None, reviewed_sha: str | None, state: dict[str, Any], manifest: Path, repo: Path | None = None) -> dict[str, Any]:
    """Freeze enforcement: a clean tree; a well-formed review id and reviewed commit SHA (40 hex); the reviewed commit is an ancestor of HEAD; the
    sha256 of the four code files at THAT commit (``git show``) equals the current files; no conflicting review id bound to that commit in the
    manifest. A Codex pass cannot be checked offline: the id and the SHA are recorded for the audit."""
    if not state["git_clean"]:
        raise SystemExit("the working tree is not clean (git status --porcelain --untracked-files=all): the evaluator refuses to run")
    if not review_id or not REVIEW_ID.match(review_id):
        raise SystemExit("--reviewed-head-review-id REVIEW-xxxxxxxx (the Codex review that passed on the freezing commit) is required")
    if not reviewed_sha or not SHA40.match(reviewed_sha):
        raise SystemExit("--reviewed-sha (40 hex: the commit that froze the four code files) is required")
    if not is_ancestor(reviewed_sha, state["git_head"], repo):
        raise SystemExit(f"the reviewed commit {reviewed_sha[:12]} is not an ancestor of HEAD {state['git_head'][:12]}")
    if hashes_at(reviewed_sha, repo) != working_hashes(repo):
        raise SystemExit("the sha256 of the code files differ between the reviewed commit and the working tree")
    prev = json.loads(manifest.read_text()) if Path(manifest).is_file() else {}
    bound = prev.get("reviews", {}).get(reviewed_sha)
    if bound is not None and bound != review_id:
        raise SystemExit(f"the manifest binds {bound} to the reviewed commit {reviewed_sha[:12]}, not {review_id}")
    return prev


def evaluate(de3_dir: Path, de5_dir: Path, de0_dir: Path | None, state: dict[str, Any], scripts: dict[str, str], *, reviewed_sha: str, ancestor: Any = None, archive: Path = ARCHIVE_JSON,
             de7_path: Path = DE7_JSON, de4_path: Path = DE4_RATIOS, de5_table_path: Path = DE5_COMMITTED, de6_report: Path = DE6_REPORT,
             h5_path: Path | None = None, expected_tables: dict[str, dict[str, str]] | None = None) -> dict[str, Any]:
    """The screen. ``de5_dir`` holds ``de5-<arm>/de5.json`` records, ``de0_dir`` the optional ``de0-<arm>/kernels.json``."""
    cmp = _cmp()
    expected = expected_tables or {}
    check_expected_tables(expected, json.loads(Path(de7_path).read_text()))
    arch = load_archive(archive)
    crosscheck_h5(arch, h5_path)
    reused = load_reused(de3_dir, de6_report)
    em = {e: np.asarray(reused[(e, "emonly")]["idd_batches"]) for e in ENERGIES}
    base = {e: np.asarray(reused[(e, "all")]["idd_batches"]) for e in ENERGIES}
    f_base = {e: cmp.curve_metrics(base[e].mean(axis=0))["peak_over_plateau"] / cmp.curve_metrics(em[e].mean(axis=0))["peak_over_plateau"] for e in ENERGIES}
    bands = compute_bands(cmp, arch, f_base)
    de4 = json.loads(Path(de4_path).read_text())
    de5_table = json.loads(Path(de5_table_path).read_text())
    for (e, v), d in reused.items():
        check_seeds(d, f"de3-{e}-{v}.json")
    runs: dict[str, Any] = {"reused": {f"de3-{e}-{v}.json": sha256_bytes_of(Path(de3_dir) / f"de3-{e}-{v}.json") for (e, v) in sorted(reused)}}
    arms: dict[str, Any] = {}
    valid_all: dict[str, bool] = {}
    for arm in ARMS:
        docs: dict[int, dict[str, Any]] = {}
        for e in ENERGIES:
            p = Path(de3_dir) / f"de3-{e}-{arm}.json"
            if p.is_file():
                docs[e] = json.loads(p.read_text())
                runs.setdefault("arms", {}).setdefault(arm, {})[str(e)] = {"sha256": sha256_bytes_of(p), "executed": True}
        if len(docs) != len(ENERGIES):
            if arm == GATED_ARM:
                raise SystemExit(f"the {arm} run documents are missing (150 and 200 MeV are required; the evaluator never runs anything)")
            arms[arm] = {"status": "not run"}
            continue
        for e, d in docs.items():
            check_seeds(d, f"de3-{e}-{arm}.json")
            check_table(d, expected.get(arm), arm, f"de3-{e}-{arm}.json")
            check_code(d, reviewed_sha, scripts, f"de3-{e}-{arm}.json", ancestor, state["git_head"])
        ok = {e: valid_document(docs[e], e, "emel", arm=True) for e in ENERGIES}
        for e in ENERGIES:
            runs["arms"][arm][str(e)]["valid"] = bool(ok[e])
        valid_all[arm] = all(ok.values())
        rec: dict[str, Any] = {"table_id": docs[150]["elastic_table_id"], "valid": valid_all[arm], "per_energy": {}}
        if valid_all[arm]:
            for e in ENERGIES:
                q = arm_quantities(cmp, arch[e], bands[e]["delta_obs"], np.asarray(docs[e]["idd_batches"]), base[e], em[e])
                q["D7_recoil"] = d7_recoil(docs[e], reused[(e, "all")], em[e], arch[e]["g_ion"] - arch[e]["g_top"])
                rec["per_energy"][str(e)] = {**q, "conditions": conditions(q, bands[e])}
        if arm in D4_SEED_BASES:
            p5 = Path(de5_dir) / f"de5-{arm}" / "de5.json"
            if p5.is_file():
                d5 = json.loads(p5.read_text())
                check_seeds(d5, str(p5), d4_arm=arm)
                check_table(d5, expected.get(arm), arm, str(p5))
                check_code(d5, reviewed_sha, scripts, str(p5), ancestor, state["git_head"])
                rec["D4"] = d4_report(arm, d5, de5_table, de4)
                runs["arms"][arm]["de5_sha256"] = sha256_bytes_of(p5)
            elif arm == GATED_ARM:
                raise SystemExit(f"{p5} is missing (D4 of the model under test)")
        if de0_dir is not None and (Path(de0_dir) / f"de0-{arm}" / "kernels.json").is_file():
            k0 = Path(de0_dir) / f"de0-{arm}" / "kernels.json"
            kd = json.loads(k0.read_text())
            check_table(kd, expected.get(arm), arm, str(k0))
            check_code(kd, reviewed_sha, scripts, str(k0), ancestor, state["git_head"])
            rec["D6"] = d6_kernels(kd, json.loads((RES / "de0" / "kernels.json").read_text()))
            runs["arms"][arm]["de0_sha256"] = sha256_bytes_of(k0)
        arms[arm] = rec
    v7 = de7_condition(de7_path)
    g = arms[GATED_ARM]
    failing: list[str] = []
    if not valid_all[GATED_ARM]:
        failing.append("run invalid: not established")
    else:
        for e in ENERGIES:
            for k, cond in g["per_energy"][str(e)]["conditions"].items():
                if not cond["holds"]:
                    failing.append(f"({k.replace('_', '-')}) at {e} MeV")
    if not v7["holds"]:
        failing.append("(v) construction record")
    established = bool(valid_all[GATED_ARM])
    return {
        "banner": c.BANNER, "disclosure": DISCLOSURE, "verdict": "CONTINUE" if not failing else "STOP", "established": established, "failing_conditions": failing,
        "iv_report_only_note": "condition (iv) is recorded, never part of the verdict; the 30 / 50 MeV cases are identical to the pinned table by construction for M2' and R-B",
        "conditions": {"v": {"construction_record": v7}},
        "bands": {str(e): bands[e] for e in ENERGIES}, "archived_se_df": {str(e): {k: arch[e][k] for k in arch[e] if k.endswith(("_se", "_df"))} for e in ENERGIES},
        "arms": arms, "runs": runs, "provenance": {**state, "reviewed_sha": reviewed_sha, "scripts_sha256": scripts, "archive_sha256": sha256_bytes_of(archive), "de7_sha256": v7["de7_sha256"]},
    }  # fmt: skip


def project_for_git(doc: dict[str, Any]) -> dict[str, Any]:
    """The Git projection: the evaluator output re-keyed for the allowlist (``arms/<arm>/per_energy/<E>/...`` and the D4 blocks)."""
    out = copy.deepcopy(doc)
    for arm in out["arms"].values():
        for key in ("D4",):
            if key in arm:
                out.setdefault("d4", {})[arm["table_id"]] = arm.pop(key)
        if "D6" in arm:
            out.setdefault("d6", {})[arm["table_id"]] = arm.pop("D6")
    return out


def render_report(doc: dict[str, Any]) -> str:
    f = lambda x, n=4: "-" if x is None else f"{x:.{n}f}"  # noqa: E731
    out = ["", doc["disclosure"], "", f"## Verdict: **{doc['verdict']}**" + ("" if doc.get("established") else " (not established: a run published an invalid document)"), ""]
    out += ["Failing conditions: " + (", ".join(doc.get("failing_conditions", [])) or "none"), "", doc["iv_report_only_note"] + ".", ""]
    out += ["## Bands (recomputed from the frozen rule)", "", "| E | F closure window | allowed |q| (i-b) | allowed |d_g| (ii) | allowed |p| (iii) |", "|---|---|---|---|---|"]
    for e, b in doc["bands"].items():
        out.append(f"| {e} | {f(b['F_closure_lo'], 3)} to {f(b['F_closure_hi'], 3)} | {f(b['a_PP'], 5)} | {f(b['a_gain'], 5)} | {f(b['a_I'], 4)} |")
    for name, a in doc["arms"].items():
        out += ["", f"## Arm {name}" + ("  (model under test)" if name == GATED_ARM else "  (report-only)"), ""]
        if a.get("status"):
            out.append(a["status"] + ".")
            continue
        out.append(f"Table id `{a['table_id']}`; documents valid: {a['valid']}.")
        out += ["", "| E | D1 closure (SE) | D1b q (SE) | D2 d_g (SE) | D3 p (SE) | (i) | (i-b) | (ii) | (iii) |", "|---|---|---|---|---|---|---|---|---|"]
        for e, q in a.get("per_energy", {}).items():
            cd = q["conditions"]
            out.append(f"| {e} | {f(q['D1_closure'], 3)} ({f(q['D1_se'], 3)}) | {f(q['D1b_q'], 4)} ({f(q['D1b_se'], 4)}) | {f(q['D2_dg'], 5)} ({f(q['D2_se'], 5)}) | "
                       f"{f(q['D3_p'], 4)} ({f(q['D3_se'], 4)}) | " + " | ".join("pass" if cd[k]["holds"] else "fail" for k in ("i", "i_b", "ii", "iii")) + " |")
    out += ["", "## D4 (report-only for (iv)): P(theta_lab > theta_2mm), arm / Geant4", ""]
    for tid, d in doc.get("d4", {}).items():
        out += [f"Table `{tid}`: moves toward 1 at {d['iv_toward_count_5']} of 5; in-domain 3 of 3: {d['iv_in_domain_3_of_3']}.", "", "| E | arm/table | table/G4 | arm/G4 | toward 1 | identical by construction |", "|---|---|---|---|---|---|"]
        for e, r in d["rows"].items():
            out.append(f"| {e} | {f(r['arm_over_table'])} | {f(r['table_over_geant4'])} | {f(r['arm_over_geant4'])} | {r['moves_toward_1']} | {r['identical_by_construction']} |")
        out.append("")
    cv = doc["conditions"]["v"]["construction_record"]
    out += [f"## (v) construction record: {'holds' if cv['holds'] else 'fails'} (de7 sha256 `{cv['de7_sha256']}`)", "",
            f"Provenance: HEAD `{doc['provenance']['git_head']}`, clean tree {doc['provenance']['git_clean']}.", "",
            f"Leak check (no absolute cross section, no measured value): passed {doc['leak_check']['passed']}."]  # fmt: skip
    return "\n".join(out) + "\n"


def absolutes_of(de5_dir: Path) -> list[float]:
    """Absolute table cross sections of the D-E5 variant records (mb, barn, fm^2): none may appear in the projection."""
    vals: list[float] = []
    for p in sorted(Path(de5_dir).glob("de5-*/de5.json")):
        for cs in json.loads(p.read_text())["cases"].values():
            vals += [float(x) for x in cs.get("sigma_el_mb", {}).values()]
    return [x * f for x in vals for f in (1.0, 1e-3, 0.1)]


def run(de3_dir: Path, de5_dir: Path, de0_dir: Path | None, out: Path, review_id: str | None, expected_tables: Path | None = None, reviewed_sha: str | None = None) -> dict[str, Any]:
    from . import de7_m2prime_record as d7

    state = c.repo_state()
    manifest_p = Path(out) / "de8-manifest.json"
    prev = check_review(review_id, reviewed_sha, state, manifest_p)
    assert reviewed_sha is not None
    published = [p for p in (Path(out) / "de8.json", Path(out) / "de8-report.md", Path(out) / "git" / "de8.json", Path(out) / "git" / "de8-report.md") if p.exists()]
    if prev.get("evaluation") or published:  # the screen is evaluated once; nothing already published is rewritten or re-derived
        raise SystemExit(f"a D-E8 evaluation is already published ({'manifest' if prev.get('evaluation') else published[0].name}): it is never overwritten")
    scripts = hashes_at(reviewed_sha)
    if expected_tables is None:
        raise SystemExit("--expected-tables FILE (the expected table id and npz sha256 per arm, checked against de7.json) is required")
    doc = evaluate(de3_dir, de5_dir, de0_dir, state, scripts, reviewed_sha=reviewed_sha, archive=ARCHIVE_JSON, de7_path=DE7_JSON, de4_path=DE4_RATIOS, de5_table_path=DE5_COMMITTED,
                   de6_report=DE6_REPORT, h5_path=H5_JSON, expected_tables=json.loads(Path(expected_tables).read_text()))
    digest = hashlib.sha256(json.dumps({k: doc[k] for k in ("verdict", "failing_conditions", "arms", "bands")}, sort_keys=True, default=float).encode()).hexdigest()
    proj = c.git_projection(json.loads(json.dumps(project_for_git(doc), default=float)), GIT_ALLOW)
    bad = d7.leak_violations(proj, absolutes_of(de5_dir))
    if bad:
        raise SystemExit(f"leak check failed: {bad[:6]}")
    doc["leak_check"] = {"passed": True}
    doc = project_for_git(doc)
    hdr = c.header("D-E8: go/no-go screen of the model under test M2' (Amendment 20 (c), Amendment 21 (e))",
                   "python -m validation.scripts.elastic_attribution.de8_diagnostic --reviewed-head-review-id " + str(review_id) + " --reviewed-sha " + str(reviewed_sha), [ARCHIVE_JSON, DE7_JSON, DE6_REPORT, DE4_RATIOS, DE5_COMMITTED])  # fmt: skip
    c.write_outputs(Path(out), "de8.json", doc, "de8-report.md", hdr, render_report, GIT_ALLOW)
    manifest = {"reviews": {**prev.get("reviews", {}), reviewed_sha: review_id}, "reviewed_sha": reviewed_sha, "head": state["git_head"], "clean": state["git_clean"],
                "scripts_sha256": scripts, "evaluation": {"digest": digest, "verdict": doc["verdict"], "established": doc["established"]},
                "runs": doc["runs"]}  # fmt: skip
    manifest_p.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    return doc


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--de3-dir", type=Path, default=c.OUT / "de3")
    ap.add_argument("--de5-dir", type=Path, default=c.OUT)
    ap.add_argument("--de0-dir", type=Path, default=c.OUT)
    ap.add_argument("--out", type=Path, default=c.OUT / "de8")
    ap.add_argument("--reviewed-head-review-id", required=True)
    ap.add_argument("--reviewed-sha", required=True, help="40-hex commit that froze the four code files (an ancestor of HEAD)")
    ap.add_argument("--expected-tables", type=Path, required=True, help='JSON {"m2prime": {"table_id": ..., "npz_sha256": ...}, "m2": {...}, "rb": {...}}')
    a = ap.parse_args(argv)
    run(a.de3_dir, a.de5_dir, a.de0_dir, a.out, a.reviewed_head_review_id, a.expected_tables, a.reviewed_sha)


if __name__ == "__main__":
    main()
