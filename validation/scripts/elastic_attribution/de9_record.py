"""D-E9 record (exploratory, report-only, never evidence): the projection of B0, B1, B2, B3 and A1 and the Lead A sensitivity table
(``validation/plans`` D-E9 plan, rev 2). Leads identified after the D-E8 STOP; this is NOT a re-evaluation: D-E8 is neither re-run nor re-derived and
no recorded verdict, band, seed or table of record changes. Every number is a gain, a ratio, a z value or a flag; no absolute cross section.

Inputs (``de3-<E>-<variant>[-r<k>].json`` harness documents; CRN seeds 20990100 / 20990101 unless a realisation seed is recorded):
* B1  the decomposition of the elastic plateau gain: ``all`` / ``pp-only`` / ``pA-only`` against ``emonly`` and the additivity;
* B3  ``pp-level-0.90-ge150``, ``pp-level-0.95-ge150``, ``pp-level-monotone`` against ``all``: the gain change of each variant and the slope
  dg / d ln s of the scaled family, as gains only. **No-inversion rule:** no function of this module computes or reports a scale or a level that
  would close any residual, or relates a variant's change to a residual;
* A1  K independent pinned-table realisations (EM+elastic, ``-r<k>``) with their EM-only partners: mean, SD, SE of the mean, the SE-calibration
  ratio rho = s / mean(per-run jackknife SE) against the pre-registered range, z of the archived baseline (step 16, ionmc) and of the D-E3 base,
  and the report-only d(mu) labelled "sensitivity; frozen d = <frozen value> stands" and never compared with a band. K is fixed (8 at 200 MeV,
  4 at 150 MeV) and is never extended;
* B2  the share of the elastic plateau gain carried by the generation-1 (p-p secondary) proton dose and by the recoil channel, from the existing
  species / generation tally channels (``--genealogy-split`` of ``de3_probes``);
* B0  the projected record of ``de9_pp_level``.

Plateau = 20-60 mm mean of the 0.5 mm IDD; gain = plateau_EM+el / plateau_EM - 1; paired delete-one-batch jackknife SE (common random numbers).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np

from . import common as c

ENERGIES = (150, 200)
A1_K = {200: 8, 150: 4}
A1_SEED_BASE = {200: 20992000, 150: 20992010}
RHO_RANGE = {8: (0.556, 1.418), 4: (0.342, 1.614)}  # 90 % range of sqrt(chi2_{K-1} / (K-1)), fixed in the plan before A1
B2_SEEDS = {150: 20992100, 200: 20992101}
B2_REALISATION = 100
SCALES = {"pp-level-0.90-ge150": 0.90, "pp-level-0.95-ge150": 0.95}
MONOTONE = "pp-level-monotone"
SENS_LABEL = "sensitivity; frozen d = {d:.5f} stands"
ARCHIVE_JSON = c.REPO / "validation" / "results" / "transport" / "lv5c-dd976c5-single-process.json"
RES = c.REPO / "validation" / "results" / "elastic_attribution"
DE9_CODE_FILES = ("de3_probes.py", "make_variant_tables.py", "de9_pp_level.py", "de9_record.py")
REVIEW_ID = re.compile(r"^REVIEW-[0-9a-f]{8}$")
SHA40 = re.compile(r"^[0-9a-f]{40}$")
SCRIPT_REL = "validation/scripts/elastic_attribution"
DISCLOSURE = "D-E9: exploratory, report-only; leads identified after the D-E8 STOP; not a re-evaluation. The D-E8 STOP stands; no band, seed or table of record changes."

GIT_ALLOW = (
    "banner", "disclosure", "b0/*", "b0/*/*", "b0/*/*/*", "b1/*/*", "b3/*/*/*", "b3/slope/*/*", "a1/*/*", "a1/*/*/*", "a1_sensitivity/**", "b2/*/*", "b2/*/*/*",
    "provenance/*", "provenance/*/*", "leak_check/*",
)  # fmt: skip
GIT_FILES = {"de9.json": (GIT_ALLOW, ())}


# ---------------------------------------------------------------------------------------------------------------------------------
def _cmp() -> Any:
    from .de3_probes import _cmp as m

    return m()


def _jk(stat: Any, arrays: list[np.ndarray], den: np.ndarray) -> tuple[float, float]:
    from .de3_probes import _jackknife

    return _jackknife(stat, arrays, den)


def plateau(cmp: Any, x: np.ndarray) -> float:
    """Mean of the 0.5 mm IDD over the 20-60 mm window (``compare_idd_v5.PLATEAU_MM``); no R80 is needed, so channels without a distal fall-off work."""
    d = cmp.depth_centres(x.size)
    return float(x[(d >= cmp.PLATEAU_MM[0]) & (d <= cmp.PLATEAU_MM[1])].mean())


def _idd(doc: dict[str, Any], key: str = "idd_batches") -> np.ndarray:
    return np.asarray(doc[key], dtype=float)


def gain(cmp: Any, emel: np.ndarray, em: np.ndarray) -> tuple[float, float]:
    """Plateau gain of the EM+elastic run over the EM-only run (same seed) and its paired jackknife SE."""
    return _jk(lambda n, d: plateau(cmp, n[0]) / plateau(cmp, d) - 1.0, [emel], em)


def gain_difference(cmp: Any, arm: np.ndarray, base: np.ndarray, em: np.ndarray) -> tuple[float, float]:
    """``g_arm - g_base`` (common random numbers: the batches are shared) and its paired jackknife SE."""
    g = lambda x, d: plateau(cmp, x) / plateau(cmp, d) - 1.0  # noqa: E731
    return _jk(lambda n, d: g(n[0], d) - g(n[1], d), [arm, base], em)


def b1(cmp: Any, docs: dict[tuple[int, str], dict[str, Any]]) -> dict[str, Any]:
    """The decomposition of the plateau gain (all / pp-only / pA-only) and its additivity, per energy."""
    out: dict[str, Any] = {}
    for e in ENERGIES:
        em = _idd(docs[(e, "emonly")])
        g = {v: gain(cmp, _idd(docs[(e, v)]), em) for v in ("all", "pp-only", "pA-only")}
        add = g["pp-only"][0] + g["pA-only"][0]
        out[str(e)] = {"gain_all": g["all"][0], "gain_all_se": g["all"][1], "gain_pp_only": g["pp-only"][0], "gain_pA_only": g["pA-only"][0],
                       "pp_share_of_sum": g["pp-only"][0] / add, "additivity_gap": g["all"][0] - add, "additivity_gap_over_all": (g["all"][0] - add) / g["all"][0]}  # fmt: skip
    return out


def b3(cmp: Any, docs: dict[tuple[int, str], dict[str, Any]]) -> dict[str, Any]:
    """Gain change of every B3 variant against ``all`` (CRN) and the slope dg/d ln s of the scaled family; gains only, no inversion."""
    out: dict[str, Any] = {}
    slope: dict[str, Any] = {}
    for e in ENERGIES:
        em, base = _idd(docs[(e, "emonly")]), _idd(docs[(e, "all")])
        dg: dict[str, tuple[float, float]] = {}
        for v in (*SCALES, MONOTONE):
            dg[v] = gain_difference(cmp, _idd(docs[(e, v)]), base, em)
            out.setdefault(v, {})[str(e)] = {"delta_gain": dg[v][0], "delta_gain_se": dg[v][1]}
        s90, s95 = SCALES["pp-level-0.90-ge150"], SCALES["pp-level-0.95-ge150"]
        slope[str(e)] = {"dg_dln_s_0.90_to_0.95": (dg["pp-level-0.95-ge150"][0] - dg["pp-level-0.90-ge150"][0]) / math.log(s95 / s90),
                         "dg_dln_s_0.95_to_1.00": (0.0 - dg["pp-level-0.95-ge150"][0]) / math.log(1.0 / s95)}  # fmt: skip
    out["slope"] = slope
    return out


def a1(cmp: Any, docs: dict[tuple[int, str, int], dict[str, Any]], arch: dict[int, dict[str, float]], de3_base: dict[int, tuple[float, float]]) -> dict[str, Any]:
    """A1 reading rule (fixed in the plan): mean, SD, SE of the mean, SE calibration rho, z of the archived and of the D-E3 baseline."""
    out: dict[str, Any] = {}
    for e in ENERGIES:
        k = A1_K[e]
        gs, ses = [], []
        for r in range(k):
            em = _idd(docs[(e, "emonly", r)])
            g, se = gain(cmp, _idd(docs[(e, "pinned", r)]), em)
            gs.append(g)
            ses.append(se)
        mu, s = float(np.mean(gs)), float(np.std(gs, ddof=1))
        se_mu = s / math.sqrt(k)
        rho = s / float(np.mean(ses))
        lo, hi = RHO_RANGE[k]
        za = (arch[e]["g_ion"] - mu) / math.hypot(arch[e]["g_ion_se"], se_mu)
        zd = (de3_base[e][0] - mu) / math.hypot(de3_base[e][1], se_mu)
        out[str(e)] = {"k": k, "mean_gain": mu, "sd": s, "se_mean": se_mu, "rho": rho, "rho_range_lo": lo, "rho_range_hi": hi,
                       "per_run_se_calibration": "per-run SE mis-calibrated" if not lo <= rho <= hi else "per-run SE consistent",
                       "z_archived": za, "z_de3_base": zd,
                       "archived_reading": "consistent with the realisation distribution" if abs(za) <= 2 else "outlier at 2 sigma",
                       "de3_reading": "consistent with the realisation distribution" if abs(zd) <= 2 else "outlier at 2 sigma"}  # fmt: skip
    return out


def a1_sensitivity(arch200: dict[str, float], delta: tuple[float, float], de3_base: tuple[float, float], a1_200: dict[str, Any] | None) -> dict[str, Any]:
    """The Lead A table at 200 MeV (report-only). Only the first row is the frozen prediction; every other row carries the label
    ``sensitivity; frozen d = <value> stands`` and is never compared with the band (``compared_with_band`` is false by construction)."""
    g_arch, se_arch, g_top = arch200["g_ion"], arch200["g_ion_se"], arch200["g_top"]
    d_frozen = g_arch + delta[0] - g_top
    lab = SENS_LABEL.format(d=d_frozen)

    def row(name: str, d: float, label: str, se: float | None = None) -> dict[str, Any]:
        return {"reading": name, "d": d, "se": se, "label": label, "compared_with_band": False}

    w_a, w_b = 1.0 / se_arch**2, 1.0 / de3_base[1] ** 2
    pooled, pooled_se = (w_a * g_arch + w_b * de3_base[0]) / (w_a + w_b), 1.0 / math.sqrt(w_a + w_b)
    rows = {
        "frozen": row("the prediction as frozen", d_frozen, "the frozen value; the STOP stands", math.hypot(se_arch, delta[1])),
        "archived_minus_2se": row("archived baseline - 2 SE", g_arch - 2 * se_arch + delta[0] - g_top, lab),
        "archived_plus_2se": row("archived baseline + 2 SE", g_arch + 2 * se_arch + delta[0] - g_top, lab),
        "de3_baseline": row("D-E3 baseline", de3_base[0] + delta[0] - g_top, lab),
        "pooled_baseline": row("inverse-variance pooled baseline", pooled + delta[0] - g_top, lab, pooled_se),
    }
    if a1_200 is not None:
        rows["a1_mean"] = row("A1 realisation mean (K fixed)", a1_200["mean_gain"] + delta[0] - g_top, lab, math.hypot(a1_200["se_mean"], delta[1]))
    return rows


def b2(cmp: Any, docs: dict[tuple[int, str], dict[str, Any]], em_docs: dict[int, dict[str, Any]]) -> dict[str, Any]:
    """Share of the elastic plateau gain carried by the secondary-proton (p-p) dose and by the recoil channel; None when the split was not recorded."""
    out: dict[str, Any] = {}
    for e in ENERGIES:
        d = docs.get((e, "pinned"))
        if d is None or "genealogy_split" not in d:
            out[str(e)] = {"status": "dropped: the genealogy split was not recorded"}
            continue
        em = _idd(em_docs[e]).mean(axis=0)
        pm = lambda x: plateau(cmp, np.asarray(x, dtype=float).mean(axis=0))  # noqa: E731
        excess = pm(_idd(d)) - plateau(cmp, em)
        sec = pm(d["genealogy_split"]["pp_secondary_proton"])
        rec = pm(d["recoil_idd_batches"])
        out[str(e)] = {"secondary_proton_share": sec / excess, "recoil_share": rec / excess, "remainder_share": 1.0 - (sec + rec) / excess, "status": "recorded"}
    return out


# ---------------------------------------------------------------------------------------------------------------------------------
def load_archived_ionmc_gain(path: Path = ARCHIVE_JSON) -> dict[int, dict[str, float]]:
    """Step 16 (V11 compare) archived ionmc gain and its own 20-batch SE, the TOPAS gain, per energy (committed single-process record)."""
    v = json.loads(Path(path).read_text())["steps"]["16-v11-compare"]["verdict"]["energies"]
    return {e: {"g_ion": v[str(e)]["plateau_gain"]["ionmc"]["value"], "g_ion_se": v[str(e)]["plateau_gain"]["ionmc"]["se"], "g_top": v[str(e)]["plateau_gain"]["reference"]["value"]} for e in ENERGIES}


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def leak_violations(proj: Any, absolutes: list[float]) -> list[str]:
    from .de7_m2prime_record import _numbers

    return [f"{p}: equals an absolute cross section" for p, x in _numbers(proj) if x not in (0.9, 0.95, 1.0) and any(a != 0.0 and abs(x - a) <= 1e-9 * abs(a) for a in absolutes)] + c.git_safety_violations(proj)


def hashes_at(sha: str, repo: Path | None = None) -> dict[str, str]:
    """sha256 of the four D-E9 code files as committed in ``sha`` (``git show``, not the working tree)."""
    from .de8_diagnostic import _git

    out = {}
    for n in DE9_CODE_FILES:
        r = _git(repo, "show", f"{sha}:{SCRIPT_REL}/{n}")
        if r.returncode != 0:
            raise SystemExit(f"{n} is not in the reviewed commit {sha[:12]}: {r.stderr.strip()}")
        out[n] = hashlib.sha256(r.stdout.encode("utf-8")).hexdigest()
    return out


def working_hashes(repo: Path | None = None) -> dict[str, str]:
    root = Path(repo or c.REPO) / SCRIPT_REL
    for n in DE9_CODE_FILES:
        if not (root / n).is_file():
            raise SystemExit(f"{root / n} is missing")
    return {n: _sha(root / n) for n in DE9_CODE_FILES}


def is_ancestor(sha: str, head: str, repo: Path | None = None) -> bool:
    from .de8_diagnostic import is_ancestor as anc

    return anc(sha, head, repo)


def check_freeze(review_id: str | None, reviewed_sha: str | None, state: dict[str, Any], manifest: Path, repo: Path | None = None) -> dict[str, Any]:
    """Freeze binding as for D-E8: a clean tree, a well-formed review id and reviewed SHA, the reviewed commit an ancestor of HEAD, the four code
    files identical between that commit and the working tree, and no conflicting review id bound to that commit in the manifest."""
    if not state["git_clean"]:
        raise SystemExit("the working tree is not clean (git status --porcelain --untracked-files=all): the D-E9 record refuses to run")
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


def check_document_code(doc: dict[str, Any], what: str, freeze: dict[str, Any]) -> None:
    """A run document must record the evaluator's HEAD, a clean tree, and (for the files it hashes) the code of the reviewed commit."""
    cd = doc.get("code") or doc.get("provenance") or {}
    if cd.get("git_head") != freeze["head"] or cd.get("git_clean") is not True:
        raise SystemExit(f"{what}: recorded HEAD / clean flag differ from the evaluator HEAD {freeze['head'][:12]} (all documents must come from one commit)")
    if not freeze["ancestor"](freeze["sha"], cd["git_head"]):
        raise SystemExit(f"{what}: recorded HEAD does not contain the reviewed commit {freeze['sha'][:12]}")
    for n, s in (cd.get("scripts_sha256") or {}).items():
        if n in freeze["scripts"] and freeze["scripts"][n] != s:
            raise SystemExit(f"{what}: sha256 of {n} recorded in the document differs from the reviewed commit")


def build_document(de3_dir: Path, de9_dir: Path, b0_path: Path | None, archive: Path = ARCHIVE_JSON, absolutes: list[float] | None = None,
                   freeze: dict[str, Any] | None = None) -> dict[str, Any]:
    """Assemble the record from the harness documents (missing documents are a SystemExit; nothing is run here)."""
    cmp = _cmp()
    arch = load_archived_ionmc_gain(archive)

    def load(p: Path, run_document: bool = False) -> dict[str, Any]:
        if not p.is_file():
            raise SystemExit(f"{p} is missing: the evaluator never runs anything")
        d = json.loads(p.read_text())
        if run_document and freeze is not None:
            check_document_code(d, p.name, freeze)
        return d  # type: ignore[no-any-return]

    shas: dict[str, str] = {}
    d3: dict[tuple[int, str], dict[str, Any]] = {}
    for e in ENERGIES:
        for v in ("all", "pp-only", "pA-only", "emonly", "m2prime", *SCALES, MONOTONE):
            p = Path(de3_dir) / f"de3-{e}-{v}.json"
            if v == "m2prime" and not p.is_file():
                continue
            d3[(e, v)] = load(p, run_document=v in ("m2prime",))
            shas[p.name] = _sha(p)
    a1d: dict[tuple[int, str, int], dict[str, Any]] = {}
    for e in ENERGIES:
        for r in range(A1_K[e]):
            for v in ("pinned", "emonly"):
                p = Path(de9_dir) / f"de3-{e}-{v}-r{r}.json"
                doc = load(p, run_document=True)
                if doc.get("seed") != A1_SEED_BASE[e] + r or doc.get("realisation") != r:
                    raise SystemExit(f"{p.name}: seed / realisation differ from the registered A1 block")
                a1d[(e, v, r)] = doc
                shas[p.name] = _sha(p)
    b2docs: dict[tuple[int, str], dict[str, Any]] = {}
    for e in ENERGIES:
        p = Path(de9_dir) / f"de3-{e}-pinned-r{B2_REALISATION}.json"
        if p.is_file():
            doc = load(p, run_document=True)
            if doc.get("seed") != B2_SEEDS[e]:
                raise SystemExit(f"{p.name}: seed differs from the registered B2 seed")
            b2docs[(e, "pinned")] = doc
            shas[p.name] = _sha(p)
    base_gain = {e: gain(cmp, _idd(d3[(e, "all")]), _idd(d3[(e, "emonly")])) for e in ENERGIES}
    a1r = a1(cmp, a1d, arch, base_gain)
    delta200 = (gain_difference(cmp, _idd(d3[(200, "m2prime")]), _idd(d3[(200, "all")]), _idd(d3[(200, "emonly")])) if (200, "m2prime") in d3 else None)
    doc: dict[str, Any] = {
        "banner": c.BANNER, "disclosure": DISCLOSURE, "b1": b1(cmp, d3), "b3": b3(cmp, d3), "a1": a1r,
        "b2": b2(cmp, b2docs, {e: d3[(e, "emonly")] for e in ENERGIES}),
    }  # fmt: skip
    if delta200 is not None:
        doc["a1_sensitivity"] = a1_sensitivity(arch[200], delta200, base_gain[200], a1r["200"])
    if b0_path is not None and Path(b0_path).is_file():
        b0 = json.loads(Path(b0_path).read_text())
        if freeze is not None:
            check_document_code(b0, Path(b0_path).name, freeze)
        doc["b0"] = {"match": b0["match"], "dip": b0["dip"], "reading": b0["reading"], "named_nodes": b0["named_nodes"]}
        shas[Path(b0_path).name] = _sha(b0_path)
    here = Path(__file__).resolve().parent
    st = c.repo_state()
    doc["provenance"] = {**st, "scripts_sha256": {n: _sha(here / n) for n in DE9_CODE_FILES if (here / n).is_file()}, "documents_sha256": shas,
                         "archive_sha256": _sha(archive), "attempts_sha256": _sha(c.ATTEMPTS) if c.ATTEMPTS.is_file() else None}  # fmt: skip
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), GIT_ALLOW)
    bad = leak_violations(proj, absolutes or [])
    if bad:
        raise SystemExit(f"leak check failed: {bad[:6]}")
    doc["leak_check"] = {"passed": True}
    return doc


def render_report(doc: dict[str, Any]) -> str:
    f = lambda x, n=5: "-" if x is None else f"{x:.{n}f}"  # noqa: E731
    out = ["", doc["disclosure"], ""]
    if "b0" in doc:
        out += [f"## B0: independent p-p level reconstruction: match {doc['b0']['match']['all_within']}, dip persists {doc['b0']['dip']['persists']}; reading: {doc['b0']['reading']}.", ""]
    out += ["## B1: decomposition of the elastic plateau gain", "", "| E | gain all | pp-only | pA-only | pp share of sum | additivity gap / all |", "|---|---|---|---|---|---|"]
    for e, r in doc["b1"].items():
        out.append(f"| {e} | {f(r['gain_all'])} | {f(r['gain_pp_only'])} | {f(r['gain_pA_only'])} | {f(r['pp_share_of_sum'], 3)} | {f(r['additivity_gap_over_all'], 3)} |")
    out += ["", "## B3: report-only variants (gain change against all, CRN; never candidates)", "", "| variant | E | delta gain | SE |", "|---|---|---|---|"]
    for v, rows in doc["b3"].items():
        if v == "slope":
            continue
        for e, r in rows.items():
            out.append(f"| {v} | {e} | {f(r['delta_gain'])} | {f(r['delta_gain_se'])} |")
    out += ["", "Slopes dg / d ln s of the scaled family: " + "; ".join(f"{e} MeV {f(s['dg_dln_s_0.90_to_0.95'])} (0.90-0.95), {f(s['dg_dln_s_0.95_to_1.00'])} (0.95-1.00)" for e, s in doc["b3"]["slope"].items()) + ".", ""]
    out += ["## A1: realisations of the pinned-table gain", "", "| E | K | mean | SD | SE mean | rho | SE calibration | z archived | z D-E3 |", "|---|---|---|---|---|---|---|---|---|"]
    for e, r in doc["a1"].items():
        out.append(f"| {e} | {r['k']} | {f(r['mean_gain'])} | {f(r['sd'])} | {f(r['se_mean'])} | {f(r['rho'], 3)} | {r['per_run_se_calibration']} | {f(r['z_archived'], 2)} ({r['archived_reading']}) | {f(r['z_de3_base'], 2)} ({r['de3_reading']}) |")
    if "a1_sensitivity" in doc:
        out += ["", "## Lead A sensitivity table at 200 MeV (report-only; never compared with the band)", "", "| reading | d | label |", "|---|---|---|"]
        for r in doc["a1_sensitivity"].values():
            out.append(f"| {r['reading']} | {f(r['d'])} | {r['label']} |")
    out += ["", "## B2: share of the elastic plateau gain", ""]
    for e, r in doc["b2"].items():
        out.append(f"- {e} MeV: " + (f"secondary-proton share {f(r['secondary_proton_share'], 3)}, recoil share {f(r['recoil_share'], 3)}, remainder {f(r['remainder_share'], 3)}" if r["status"] == "recorded" else r["status"]))
    out += ["", f"Leak check (no absolute cross section): passed {doc['leak_check']['passed']}."]
    return "\n".join(out) + "\n"


def run(de3_dir: Path, de9_dir: Path, out: Path, b0_path: Path | None = None, review_id: str | None = None, reviewed_sha: str | None = None,
        review_id_full: str | None = None) -> dict[str, Any]:
    out = Path(out)
    state = c.repo_state()
    manifest_p = out / "de9-manifest.json"
    prev = check_freeze(review_id, reviewed_sha, state, manifest_p)
    assert reviewed_sha is not None
    if any((out / n).exists() for n in ("de9.json", "git/de9.json", "git/de9-report.md")) or prev.get("record"):
        raise SystemExit("a D-E9 record is already published: it is never overwritten")
    freeze = {"sha": reviewed_sha, "head": state["git_head"], "scripts": hashes_at(reviewed_sha), "ancestor": is_ancestor}
    doc = build_document(de3_dir, de9_dir, b0_path, freeze=freeze)
    doc["provenance"]["reviewed_sha"] = reviewed_sha
    doc["provenance"]["review_id"] = review_id
    hdr = c.header("D-E9 record: B0 / B1 / B2 / B3 / A1 and the Lead A sensitivity table", "python -m validation.scripts.elastic_attribution.de9_record", [ARCHIVE_JSON])
    c.write_outputs(out, "de9.json", doc, "de9-report.md", hdr, render_report, GIT_ALLOW)
    manifest = {"reviews": {**prev.get("reviews", {}), reviewed_sha: review_id}, "review_id_full_text": review_id_full, "reviewed_sha": reviewed_sha, "head": state["git_head"],
                "scripts_sha256": freeze["scripts"], "record": {"documents_sha256": doc["provenance"]["documents_sha256"]}}  # fmt: skip
    manifest_p.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    return doc


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--de3-dir", type=Path, default=c.OUT / "de3")
    ap.add_argument("--de9-dir", type=Path, default=c.OUT / "de9")
    ap.add_argument("--b0", type=Path, default=c.OUT / "de9" / "git" / "de9-b0.json")
    ap.add_argument("--out", type=Path, default=c.OUT / "de9")
    ap.add_argument("--reviewed-head-review-id", required=True, help="REVIEW-xxxxxxxx (8 hex)")
    ap.add_argument("--reviewed-sha", required=True, help="40-hex commit that froze the four code files (an ancestor of HEAD)")
    ap.add_argument("--review-id-full", default=None, help="free text: the full Codex review id")
    a = ap.parse_args(argv)
    run(a.de3_dir, a.de9_dir, a.out, a.b0, a.reviewed_head_review_id, a.reviewed_sha, a.review_id_full)


if __name__ == "__main__":
    main()
