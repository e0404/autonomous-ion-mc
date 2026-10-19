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
    "complete", "incomplete_reasons/*", "invalid_documents/*", "provenance/*", "provenance/*/*", "leak_check/*",
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
    if "code" in doc:  # a harness document: the code hashes are REQUIRED and must equal those of the reviewed commit
        sh = cd.get("scripts_sha256")
        if not isinstance(sh, dict) or not sh:
            raise SystemExit(f"{what}: no scripts_sha256 recorded (the code binding cannot be verified)")
        for n, want in freeze["scripts"].items():
            if n in c.CODE_FILES and sh.get(n) != want:
                raise SystemExit(f"{what}: sha256 of {n} recorded in the document is missing or differs from the reviewed commit")
    elif "de9_pp_level.py" in freeze["scripts"] and cd.get("code_sha256") != freeze["scripts"]["de9_pp_level.py"]:  # the B0 record
        raise SystemExit(f"{what}: code_sha256 is missing or differs from the reviewed de9_pp_level.py")


REUSED = ("all", "pp-only", "pA-only", "emonly")
CRN = {150: 20990100, 200: 20990101}
N_HIST, N_BATCH = 100_000, 20
HASH_SOURCES = (RES / "de6" / "de6-report.md", RES / "de3" / "de3-summary.md")  # committed input lists that bind the reused D-E3 documents by sha256


def committed_hashes(sources: tuple[Path, ...] | None = None) -> dict[str, str]:
    """``{file name: sha256}`` of the committed input lists (de6-report.md, de3-summary.md); a name listed with two different hashes is a SystemExit."""
    out: dict[str, str] = {}
    for p in sources or HASH_SOURCES:
        for m in re.finditer(r"^- `([^`]+)`\s+([0-9a-f]{64})\s*$", Path(p).read_text(), re.M):
            n = Path(m.group(1)).name
            if out.setdefault(n, m.group(2)) != m.group(2):
                raise SystemExit(f"{n} is listed with different sha256 values in the committed input lists")
    return out


ATTEMPTS = c.ATTEMPTS
ATTEMPTS_DEFAULT = c.ATTEMPTS
DE8_RECORD = RES / "de8" / "de8.json"
PIN_JSON = c.REPO / "src" / "ionmc" / "data" / "elastic_table_pin.json"
UNIT = "MeV/(g/cm^2)/primary"


def expected_geometry(energy: int) -> dict[str, Any]:
    """The frozen V5 scoring geometry of the harness (r = 20 cm, 0.5 mm depth bins) at ``energy``, from the transport step module (read only)."""
    v5b = c.import_transport_module("steps_v5b")
    _geo, grid, nz, _r80 = v5b.v5_geometry(float(energy))  # the scoring grid: lateral half-width 200 mm (r = 20 cm), 0.5 mm depth bins
    return {"origin_mm": [float(x) for x in grid.origin_mm], "spacing_mm": [float(x) for x in grid.spacing_mm], "shape": [int(x) for x in grid.shape], "name": grid.name,
            "nz": int(nz), "bin_mm": float(v5b.V5_DZ_MM)}  # fmt: skip


def check_arrays(doc: dict[str, Any], what: str, energy: int) -> None:
    """Frozen configuration of the arrays: the scoring geometry, 0.5 mm bins and the bin count, the units, the tally names, the array shapes (batches x bins)
    and finiteness (SystemExit otherwise)."""
    g = expected_geometry(energy)
    sc = doc.get("scoring") or {}
    for k in ("origin_mm", "spacing_mm", "shape", "name"):
        if sc.get(k) != g[k]:
            raise SystemExit(f"{what}: scoring {k} {sc.get(k)!r} differs from the frozen V5 geometry {g[k]!r}")
    if doc.get("bin_mm") != g["bin_mm"] or doc.get("unit") != UNIT:
        raise SystemExit(f"{what}: bin size / unit differ from the frozen configuration ({g['bin_mm']} mm, {UNIT})")
    split = doc.get("genealogy_split")
    names = [] if doc.get("mode") == "emonly" else ["recoil", *(["pp_secondary_proton", "primary_proton"] if split else [])]
    if doc.get("tallies") != names:
        raise SystemExit(f"{what}: tally names {doc.get('tallies')!r} differ from the expected {names!r}")
    arrays = {"idd_batches": doc.get("idd_batches")}
    if doc.get("mode") != "emonly":
        arrays["recoil_idd_batches"] = doc.get("recoil_idd_batches")
    for k, v in (split or {}).items():
        arrays[f"genealogy_split/{k}"] = v
    for k, v in arrays.items():
        a = np.asarray(v, dtype=float) if v is not None else None
        if a is None or a.shape != (N_BATCH, g["nz"]) or not np.all(np.isfinite(a)):
            raise SystemExit(f"{what}: array {k} is missing, not shaped (batches x bins) = ({N_BATCH}, {g['nz']}), or not finite")


def attempt_log_sha(log: list[dict[str, Any]], name: str) -> str:
    """The sha256 recorded in the attempt log for the published document ``name`` (exactly one published attempt; SystemExit otherwise)."""
    ids = {r["id"] for r in log if r.get("outcome") == "started" and Path(str(r.get("output"))).name == name}
    pub = [r for r in log if r.get("attempt") in ids and str(r.get("outcome", "")).startswith("published")]
    if len(pub) != 1 or not pub[0].get("document_sha256"):
        raise SystemExit(f"{name}: the attempt log has no (unique) published entry with a document sha256")
    return str(pub[0]["document_sha256"])


BUILDER_IDENTITY = {  # exploratory model / Pauli variant / Pauli E_min of the official builder tables (Amendments 20, 21; the D-E7 record)
    "m2prime": ("m2prime-pauli-blocked", "clementel-villi", 80.0),
    "m2": ("m2-optical-glauber", None, None),
    "rb": ("m2prime-pauli-blocked", "bdc-pn", 80.0),
}
LA150_ZIP_SHA256 = "27bcafb89cf0444c53c6b9f3dd17618c62e2e6b3694f70d31c4f502399f103f7"  # the registered ENDF/B-VIII.0 proton sublibrary (``ionmc.data.registry``)
LA150_H1_MEMBER_SHA256 = "a1c3bc637ba242f263d88a3b3848b5a0a717728d6ec19395d24225e1b96a4fb5"  # ``p-001_H_001.endf`` inside it (frozen pin)


def verify_expected_tables(expected: dict[str, dict[str, str]], scratch: Path) -> None:
    """The expected-tables map is verified against the scratch cache SEMANTICALLY: for every variant the installed npz sha256, the sidecar's recorded hash and the
    table id recomputed from the sidecar must agree, and the sidecar itself must define the variant: ``exploratory.d_e3_variant`` equal to the name and its
    parameters equal to the frozen definition of ``make_variant_tables`` (and the pinned source id) for the D-E3 variants; the exploratory model, Pauli variant and
    Pauli E_min of the builder tables (``BUILDER_IDENTITY``); ``pinned`` against the committed pin and no exploratory block. ``variants.json`` is only a consistency
    hint, never the authority. ``emonly`` runs no elastic table and is exempt."""
    from ionmc.nuclear.build import table_id as tid_fn

    from .make_variant_tables import _parameters

    vj = json.loads((Path(scratch) / "variants.json").read_text())
    ids = {**{v: r["table_id"] for v, r in vj["variants"].items()}, **{v: r["table_id"] for v, r in vj.get("builders", {}).items()}}
    pin = json.loads(PIN_JSON.read_text())
    der = Path(scratch) / "derived"
    for v, e in expected.items():
        if v == "emonly":
            continue
        tid = e.get("table_id")
        if v != "pinned" and ids.get(v) != tid:
            raise SystemExit(f"expected table of {v!r}: the variant is not in the scratch cache variants.json with this table id")
        npz, sidecar = der / f"elastic-proton-{tid}.npz", der / f"elastic-proton-{tid}.json"
        if not npz.is_file() or not sidecar.is_file():
            raise SystemExit(f"expected table of {v!r}: the table files are not installed in the scratch cache")
        side = json.loads(sidecar.read_text())
        sha = _sha(npz)
        if e.get("npz_sha256") != sha or side.get("npz_sha256") != sha or tid_fn({**side, "npz_sha256": sha}) != tid:
            raise SystemExit(f"expected table of {v!r}: table id / npz sha256 are not consistent with the scratch cache (npz and sidecar)")
        ex = side.get("exploratory")
        if v == "pinned":
            if tid != c.TABLE_ID or tid != pin["table_id"] or sha != pin["npz_sha256"] or ex is not None:
                raise SystemExit("expected pinned table differs from the committed pin (or carries an exploratory block)")
        elif v in BUILDER_IDENTITY:
            model, pauli, emin = BUILDER_IDENTITY[v]
            if not isinstance(ex, dict) or (ex.get("model"), ex.get("pauli_variant"), ex.get("pauli_emin_mev")) != (model, pauli, emin):
                raise SystemExit(f"expected table of {v!r}: the sidecar does not define the {v} construction (model / Pauli variant / Pauli E_min)")
        else:
            if not isinstance(ex, dict) or ex.get("d_e3_variant") != v or ex.get("parameters") != _parameters(v) or ex.get("source_table_id") != c.TABLE_ID:
                raise SystemExit(f"expected table of {v!r}: the sidecar does not define this variant (d_e3_variant / frozen parameters / source table)")


B0_KEYS = ("match", "dip", "reading", "carry_up_check", "quadrature_convergence", "provenance", "named_nodes")


def validate_b0(b0: dict[str, Any], what: str, expected: dict[str, dict[str, str]]) -> None:
    """Schema and provenance of the (projected) B0 record: the required blocks, ``quadrature_convergence.accepted`` as a boolean, the code hash, the pinned table (id and
    npz sha256 against the expected / pinned map), the raw LA150 zip and H-1 member hashes against the frozen pins, and the single-look disclosure."""
    from .de9_pp_level import PRE_RUN_LOOKS

    missing = [k for k in B0_KEYS if k not in b0]
    if missing:
        raise SystemExit(f"{what}: the B0 record lacks {missing}")
    if not isinstance(b0["match"].get("all_within"), bool) or not isinstance(b0["dip"].get("persists"), bool) or not isinstance(b0["reading"], str) or "match" not in b0["carry_up_check"]:
        raise SystemExit(f"{what}: the match / dip / reading / carry-up blocks are malformed")
    if not isinstance(b0["quadrature_convergence"].get("accepted"), bool):
        raise SystemExit(f"{what}: quadrature_convergence.accepted is missing or not a boolean")
    pv = b0["provenance"]
    pinned = expected.get("pinned")
    if pinned is None:
        raise SystemExit(f"{what}: no expected pinned table to bind the B0 record to")
    for k in ("code_sha256", "table_id", "npz_sha256", "source_sha256", "member_sha256", "pre_run_looks"):
        if not pv.get(k):
            raise SystemExit(f"{what}: provenance lacks {k}")
    if (pv["table_id"], pv["npz_sha256"]) != (pinned["table_id"], pinned["npz_sha256"]):
        raise SystemExit(f"{what}: the pinned table of the B0 record differs from the expected pinned table")
    if pv["source_sha256"] != LA150_ZIP_SHA256 or pv["member_sha256"] != LA150_H1_MEMBER_SHA256:
        raise SystemExit(f"{what}: the raw LA150 sublibrary / H-1 member hashes differ from the frozen pins")
    if pv["pre_run_looks"] != PRE_RUN_LOOKS:
        raise SystemExit(f"{what}: the pre-run disclosure differs from the frozen text")


def validate_document(doc: dict[str, Any], what: str, *, energy: int, variant: str, seed: int, expected: dict[str, dict[str, str]], reused: dict[str, str] | None = None,
                      path: Path | None = None, freeze: dict[str, Any] | None = None) -> bool:
    """Fail-closed validation of one harness document: frozen configuration (energy, variant, 100000 histories, 20 batches, mode, warp-cpu, float64, seed),
    elastic table identity against the expected map (id and npz sha256), the committed sha256 for the reused D-E3 documents, the code provenance for new
    documents. Returns the document's own validity (``valid`` and the energy balance): an invalid published run is returned as False and is never treated as valid."""
    mode = "emonly" if variant == "emonly" else "emel"
    want = {"energy": energy, "variant": variant, "n": N_HIST, "n_batches": N_BATCH, "mode": mode, "backend": "warp-cpu", "precision": "float64", "seed": seed}
    for k, v in want.items():
        if doc.get(k) != v:
            raise SystemExit(f"{what}: {k} {doc.get(k)!r} differs from the frozen configuration {v!r}")
    is_reused = reused is not None and variant in REUSED and path is not None and path.name in reused
    check_arrays(doc, what, energy)
    if variant == "emonly":  # no elastic table: exempt from the table identity, but a new document must say so explicitly
        if not is_reused and doc.get("elastic") is not False:
            raise SystemExit(f"{what}: an EM-only document must carry the explicit marker elastic = false")
    else:
        exp = expected.get(variant)
        if exp is None or not exp.get("table_id"):
            raise SystemExit(f"{what}: no expected table for variant {variant!r} (manifest input required)")
        if not is_reused and doc.get("elastic") is not True:
            raise SystemExit(f"{what}: an elastic document must carry the explicit marker elastic = true")
        if doc.get("elastic_table_id") != exp["table_id"]:
            raise SystemExit(f"{what}: elastic table id differs from the expected {variant} table")
        if "elastic_table_npz_sha256" in doc or not is_reused:
            if doc.get("elastic_table_npz_sha256") != exp.get("npz_sha256") or not exp.get("npz_sha256"):
                raise SystemExit(f"{what}: elastic table npz sha256 differs from (or is missing against) the expected {variant} table")
    if reused is not None and variant in REUSED and path is not None:
        if reused.get(path.name) != _sha(path):
            raise SystemExit(f"{what}: sha256 differs from (or is not listed in) the committed input lists")
    elif freeze is not None:
        check_document_code(doc, what, freeze)
    return bool(doc.get("valid")) and doc.get("relative_residual", float("inf")) <= 1e-12


def build_document(de3_dir: Path, de9_dir: Path, b0_path: Path | None, archive: Path = ARCHIVE_JSON, absolutes: list[float] | None = None,
                   freeze: dict[str, Any] | None = None, expected_tables: dict[str, dict[str, str]] | None = None, hash_sources: tuple[Path, ...] | None = None,
                   allow_incomplete: bool = False, scratch_cache: Path | None = None, attempts_path: Path | None = None, de8_record: Path | None = None,
                   expected_sha256: str | None = None) -> dict[str, Any]:
    """Assemble the record from the harness documents. Every input is validated fail-closed (``validate_document``). A missing required input (the
    ``m2prime`` document of the Lead A table, the B0 result) is a SystemExit naming it, unless ``allow_incomplete``: then the section is omitted and the
    record is labelled INCOMPLETE (``complete`` false) and cannot be read as the D-E9 synthesis. A published-invalid document is recorded and no section
    is computed from invalid data. Nothing is run here."""
    cmp = _cmp()
    arch = load_archived_ionmc_gain(archive)
    expected = expected_tables or {}
    if scratch_cache is None:
        raise SystemExit("the scratch cache (--scratch-cache) is required: the expected-tables map is verified against it")
    verify_expected_tables(expected, scratch_cache)
    reused = committed_hashes(hash_sources)
    log = json.loads(Path(attempts_path or ATTEMPTS_DEFAULT).read_text()) if Path(attempts_path or ATTEMPTS_DEFAULT).is_file() else None
    de8_sha: dict[str, str] | None = None
    reasons: list[str] = []
    invalid: dict[str, str] = {}
    shas: dict[str, str] = {}

    def load(p: Path, *, energy: int, variant: str, seed: int, required: bool = True, bound: bool = False) -> dict[str, Any] | None:
        if not p.is_file():
            if required:
                raise SystemExit(f"{p} is missing: a required input of the D-E9 record (the evaluator never runs anything)")
            return None
        d = json.loads(p.read_text())
        if not bound:  # authenticate the generated document: D-E8 arm documents against the committed D-E8 record, new runs against the attempt log
            nonlocal de8_sha
            if variant == "m2prime":
                if de8_sha is None:
                    rec8 = json.loads(Path(de8_record or DE8_RECORD).read_text())["runs"]["arms"]["m2prime"]
                    de8_sha = {str(k): v["sha256"] for k, v in rec8.items() if isinstance(v, dict) and "sha256" in v}
                if de8_sha.get(str(energy)) != _sha(p):
                    raise SystemExit(f"{p.name}: sha256 differs from the committed D-E8 record")
            else:
                if log is None:
                    raise SystemExit("the attempt log is missing: generated documents cannot be authenticated")
                if attempt_log_sha(log, p.name) != _sha(p):
                    raise SystemExit(f"{p.name}: sha256 differs from the attempt-log entry of its published run")
        if not validate_document(d, p.name, energy=energy, variant=variant, seed=seed, expected=expected, reused=reused if bound else None, path=p, freeze=freeze):
            invalid[p.name] = "published-invalid"
        shas[p.name] = _sha(p)
        return d

    d3: dict[tuple[int, str], dict[str, Any]] = {}
    for e in ENERGIES:
        for v in (*REUSED, *SCALES, MONOTONE):
            dd = load(Path(de3_dir) / f"de3-{e}-{v}.json", energy=e, variant=v, seed=CRN[e], bound=v in REUSED)
            assert dd is not None
            d3[(e, v)] = dd
    m2p = load(Path(de3_dir) / "de3-200-m2prime.json", energy=200, variant="m2prime", seed=CRN[200], required=not allow_incomplete)
    if m2p is None:
        reasons.append("missing input de3-200-m2prime.json: the Lead A sensitivity table is not computed")
    a1d: dict[tuple[int, str, int], dict[str, Any]] = {}
    for e in ENERGIES:
        for r in range(A1_K[e]):
            for v in ("pinned", "emonly"):
                p = Path(de9_dir) / f"de3-{e}-{v}-r{r}.json"
                dd = load(p, energy=e, variant=v, seed=A1_SEED_BASE[e] + r)
                assert dd is not None
                if dd.get("realisation") != r:
                    raise SystemExit(f"{p.name}: realisation differs from the registered A1 block")
                a1d[(e, v, r)] = dd
    b2docs: dict[tuple[int, str], dict[str, Any]] = {}
    for e in ENERGIES:
        dd = load(Path(de9_dir) / f"de3-{e}-pinned-r{B2_REALISATION}.json", energy=e, variant="pinned", seed=B2_SEEDS[e], required=False)
        if dd is not None:
            b2docs[(e, "pinned")] = dd
    b0 = None
    if b0_path is not None and Path(b0_path).is_file():
        b0 = json.loads(Path(b0_path).read_text())
        if log is None:
            raise SystemExit("the attempt log is missing: the B0 result cannot be authenticated")
        if attempt_log_sha(log, Path(b0_path).name) != _sha(Path(b0_path)):
            raise SystemExit(f"{Path(b0_path).name}: sha256 differs from the attempt-log entry of its published run")
        validate_b0(b0, Path(b0_path).name, expected)
        if freeze is not None:
            check_document_code(b0, Path(b0_path).name, freeze)
        shas[Path(b0_path).name] = _sha(b0_path)
        if b0["quadrature_convergence"]["accepted"] is not True:
            reasons.append("B0 quadrature convergence not accepted (pre-declared acceptance: relative difference < 1e-6)")
    elif not allow_incomplete:
        raise SystemExit(f"{b0_path} is missing: the B0 result is a required input of the D-E9 record")
    else:
        reasons.append("missing input: the B0 result is not in the record")
    doc: dict[str, Any] = {"banner": c.BANNER, "disclosure": DISCLOSURE}
    if invalid:
        reasons.append("published-invalid run documents: no section is computed from invalid data")
    else:
        base_gain = {e: gain(cmp, _idd(d3[(e, "all")]), _idd(d3[(e, "emonly")])) for e in ENERGIES}
        a1r = a1(cmp, a1d, arch, base_gain)
        doc.update({"b1": b1(cmp, d3), "b3": b3(cmp, d3), "a1": a1r, "b2": b2(cmp, b2docs, {e: d3[(e, "emonly")] for e in ENERGIES})})
        if m2p is not None:
            delta200 = gain_difference(cmp, _idd(m2p), _idd(d3[(200, "all")]), _idd(d3[(200, "emonly")]))
            doc["a1_sensitivity"] = a1_sensitivity(arch[200], delta200, base_gain[200], a1r["200"])
        if b0 is not None:
            doc["b0"] = {"match": b0["match"], "dip": b0["dip"], "reading": b0["reading"], "named_nodes": b0["named_nodes"]}
            doc["b0"]["carry_up_match"] = b0["carry_up_check"]["match"]
            doc["b0"]["quadrature_accepted"] = b0["quadrature_convergence"]["accepted"]
    doc["complete"] = not reasons
    doc["incomplete_reasons"] = reasons
    if invalid:
        doc["invalid_documents"] = invalid
    here = Path(__file__).resolve().parent
    st = c.repo_state()
    doc["provenance"] = {**st, "scripts_sha256": {n: _sha(here / n) for n in DE9_CODE_FILES if (here / n).is_file()}, "documents_sha256": shas,
                         "archive_sha256": _sha(archive), "attempts_sha256": _sha(Path(attempts_path or ATTEMPTS_DEFAULT)) if Path(attempts_path or ATTEMPTS_DEFAULT).is_file() else None,
                         "expected_tables_sha256": expected_sha256}  # fmt: skip
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), GIT_ALLOW)
    bad = leak_violations(proj, absolutes or [])
    if bad:
        raise SystemExit(f"leak check failed: {bad[:6]}")
    doc["leak_check"] = {"passed": True}
    return doc


def render_report(doc: dict[str, Any]) -> str:
    f = lambda x, n=5: "-" if x is None else f"{x:.{n}f}"  # noqa: E731
    out = ["", doc["disclosure"], ""]
    if not doc.get("complete", False):
        out += ["**INCOMPLETE — not the D-E9 synthesis.** " + "; ".join(doc.get("incomplete_reasons", [])), ""]
    for n in doc.get("invalid_documents", {}):
        out.append(f"- {n}: published-invalid (recorded as executed; nothing is computed from it)")
    if "b1" not in doc:
        out += ["", f"Leak check (no absolute cross section): passed {doc['leak_check']['passed']}."]
        return "\n".join(out) + "\n"
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
        review_id_full: str | None = None, expected_tables: Path | None = None, allow_incomplete: bool = False, scratch_cache: Path | None = None) -> dict[str, Any]:
    out = Path(out)
    state = c.repo_state()
    manifest_p = out / "de9-manifest.json"
    prev = check_freeze(review_id, reviewed_sha, state, manifest_p)
    assert reviewed_sha is not None
    if any((out / n).exists() for n in ("de9.json", "git/de9.json", "git/de9-report.md")) or prev.get("record"):
        raise SystemExit("a D-E9 record is already published: it is never overwritten")
    freeze = {"sha": reviewed_sha, "head": state["git_head"], "scripts": hashes_at(reviewed_sha), "ancestor": is_ancestor}
    if expected_tables is None:
        raise SystemExit("--expected-tables FILE (variant -> table_id / npz_sha256 of the scratch cache) is required")
    doc = build_document(de3_dir, de9_dir, b0_path, freeze=freeze, expected_tables=json.loads(Path(expected_tables).read_text()), allow_incomplete=allow_incomplete,
                         scratch_cache=scratch_cache, attempts_path=ATTEMPTS_DEFAULT, de8_record=DE8_RECORD, expected_sha256=_sha(Path(expected_tables)))
    doc["provenance"]["reviewed_sha"] = reviewed_sha
    doc["provenance"]["review_id"] = review_id
    hdr = c.header("D-E9 record: B0 / B1 / B2 / B3 / A1 and the Lead A sensitivity table" if doc["complete"] else "INCOMPLETE — not the D-E9 synthesis", "python -m validation.scripts.elastic_attribution.de9_record", [ARCHIVE_JSON])
    c.write_outputs(out, "de9.json", doc, "de9-report.md", hdr, render_report, GIT_ALLOW)
    manifest = {"reviews": {**prev.get("reviews", {}), reviewed_sha: review_id}, "review_id_full_text": review_id_full, "reviewed_sha": reviewed_sha, "head": state["git_head"],
                "scripts_sha256": freeze["scripts"], "complete": doc["complete"], "record": {"documents_sha256": doc["provenance"]["documents_sha256"]}}  # fmt: skip
    manifest_p.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    return doc


def main(argv: list[str] | None = None) -> None:
    from .make_variant_tables import resolve_scratch

    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--de3-dir", type=Path, default=c.OUT / "de3")
    ap.add_argument("--de9-dir", type=Path, default=c.OUT / "de9")
    ap.add_argument("--b0", type=Path, default=c.OUT / "de9" / "git" / "de9-b0.json")
    ap.add_argument("--out", type=Path, default=c.OUT / "de9")
    ap.add_argument("--reviewed-head-review-id", required=True, help="REVIEW-xxxxxxxx (8 hex)")
    ap.add_argument("--reviewed-sha", required=True, help="40-hex commit that froze the four code files (an ancestor of HEAD)")
    ap.add_argument("--review-id-full", default=None, help="free text: the full Codex review id")
    ap.add_argument("--expected-tables", type=Path, required=True, help='JSON {"all": {"table_id": ..., "npz_sha256": ...}, "pinned": ..., "emonly": ..., "pp-only": ..., "pA-only": ..., "pp-level-...": ..., "m2prime": ...}')
    ap.add_argument("--scratch-cache", type=Path, default=None, help="scratch cache of make_variant_tables (default: IONMC_DE3_SCRATCH or the repository default)")
    ap.add_argument("--allow-incomplete", action="store_true", help="omit a missing m2prime document / B0 result and label the record INCOMPLETE (not the D-E9 synthesis)")
    a = ap.parse_args(argv)
    run(a.de3_dir, a.de9_dir, a.out, a.b0, a.reviewed_head_review_id, a.reviewed_sha, a.review_id_full, a.expected_tables, a.allow_incomplete, resolve_scratch(a.scratch_cache))


if __name__ == "__main__":
    main()
