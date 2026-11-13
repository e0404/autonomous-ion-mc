"""V3-008 row V5-SS (branch B): seal, run validity and the single grade.

Entry points::

    python validation/scripts/ion_em/ss_evaluate.py seal --implementation-sha SHA
        --host-records DIR --samples-root DIR [--out validation/results/transport/v8-ss/...]
    python validation/scripts/ion_em/ss_evaluate.py check-run --run-id REF-... --bundle NAME
        [--records-dir DIR] [--outputs DIR] [--out REPORT]
    python validation/scripts/ion_em/ss_evaluate.py grade [--outputs DIR] [--samples-root DIR]
        [--out validation/results/transport/v8-ss/ss-report.json]

Order of operations, file names and decisions D1-D6: "V5-SS execution note (2026-11-13)" of
``validation/plans/v3-008-acceptance.md``. Record verification reuses ``r2_io`` (committed reader,
Ed25519 attestation via ``infrastructure/experiment_v3/attest.py``, artifact listing, published
host-run output manifest, step-document parsing) and the R2-TOPAS-4 chronology check of
``r2_evaluate.check_run_validity``; ``ss_budget`` supplies the window factor ``WINDOW_K`` and the
gate (both frozen in P9).

**Statistic (both sides).** Pooled windowed rms about zero of the projected exit angles,
theta_w(W) = sqrt(sum_i [tx_i^2 1(|tx_i|<=W) + ty_i^2 1(|ty_i|<=W)] /
sum_i [1(|tx_i|<=W) + 1(|ty_i|<=W)]), with W_c = WINDOW_K * theta_anchor,c and theta_anchor,c the
untruncated pooled rms of the sealed ionmc sample of case c (decision D2). TOPAS:
theta = atan2(u, w), w = +sqrt(1 - u^2 - v^2) (the phase space stores u, v and a "w negative"
flag that must be false); ionmc: atan2(u_x, u_z) (``steps_v8.py ss-sample``).

**Pairs.** dev = (theta_model,A / theta_model,B) / (theta_ref,A / theta_ref,B) - 1 for He/C 150
(22/25), He/O 150 (22/26), He/C 250 (28/27), He/O 250 (28/31), C/O 430 (29/30). sigma = SD over
B = 2000 bootstrap replicates of dev (decision D5): TOPAS resampled with replacement within each
part (stratified by part), ionmc over histories, every (case, side) resampled independently with
its own child generator of ``np.random.SeedSequence(ANALYSIS_SEED)`` (an ANALYSIS seed, never a
transport seed). Verdict per pair ``ss_budget.gate(100 dev, 100 sigma)``; the row passes iff all
five pairs PASS (INCONCLUSIVE counts as FAIL).

**Binding SS-SEAL-1 (mirror of R2-TOPAS-4).** A V5-SS TOPAS run (smoke or evidence) is valid
only if (a) the SS seal commit (``ss_seal_sha`` of the plan's single fenced json block holding
exactly ``ss_seal_sha`` / ``ss_seal_review_id``) is an ancestor of the run's signed ``code_sha``
(``git merge-base --is-ancestor``) and (b) the SS seal review's ``finished_at`` (review record
:data:`REVIEW_PATH`, pinned at the run's ``code_sha``) precedes the run's signed
``started_at`` (= controlled ``reference_started``). Missing records fail closed (ERROR); a
violation is a BLINDING FAIL. The SS frozen manifest at the run's ``code_sha`` must equal HEAD's
(smoke: its own entry only).

Phase-space layout (header column names, Z of the exit surface in world coordinates) and log
string forms marked TO VERIFY below are checked by the scorer smoke; an unexpected column list is
a fail-closed error naming the observed list. Reports carry statistics, ids and digests only,
never per-primary data.
"""

from __future__ import annotations

import argparse
import hashlib
import io as _io
import json
import math
import re
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r2_evaluate as r2e  # noqa: E402
import r2_io as io  # noqa: E402
import ss_budget  # noqa: E402
from ss_wrapper_template import (  # noqa: E402,F401
    PLACEHOLDER,
    SS_WRAPPER_SHA256,
    SS_WRAPPER_TEMPLATE,
    render_wrapper,
)

RecordError = io.RecordError

# --- constants (D1-D6, plan "V5-SS execution note (2026-11-13)") -------------------------------
SS_CASES: tuple[int, ...] = (22, 25, 26, 27, 28, 29, 30, 31)
SS_R_INDEX = {c: i for i, c in enumerate(SS_CASES)}
SS_SEED_BASE = 20631004  # RERUN_BASE_2, first V5-SS run (section-C note 2026-11-13)
SS_HISTORIES = 400_000
SS_SCALE = "1"
TOPAS_BASE = 20619000
PARTS = (0, 1)
PART_HISTORIES = 30_000
# case -> (ion tag, MeV/u, PDG code of the primary, Geant4 particle of the EM printout)
CASE_INFO: dict[int, tuple[str, int, int, str]] = {
    22: ("he", 150, 1000020040, "alpha"),
    25: ("carbon", 150, 1000060120, "GenericIon"),
    26: ("oxygen", 150, 1000080160, "GenericIon"),
    27: ("carbon", 250, 1000060120, "GenericIon"),
    28: ("he", 250, 1000020040, "alpha"),
    29: ("carbon", 430, 1000060120, "GenericIon"),
    30: ("oxygen", 430, 1000080160, "GenericIon"),
    31: ("oxygen", 250, 1000080160, "GenericIon"),
}
# gated pairs (A, B) with the pre-registered expected dev (P10 section 7.5.1, percent)
PAIRS: tuple[tuple[str, int, int, float], ...] = (
    ("He/C 150", 22, 25, -0.49),
    ("He/O 150", 22, 26, -0.69),
    ("He/C 250", 28, 27, -0.38),
    ("He/O 250", 28, 31, -0.57),
    ("C/O 430", 29, 30, -0.17),
)
EXPECTED_ABS_FACTOR_PCT = (3.1, 4.4)  # reference-model factor of theta_ionmc/theta_SS (report only)
WINDOW_K = ss_budget.WINDOW_K
N_BOOT = 2000
ANALYSIS_SEED = 20261113  # analysis (bootstrap) seed, NOT a transport seed (decision D5)
SMOKE_BUNDLE = "v5ss-helium-water-250mevu-scorer-smoke"
SMOKE_SEEDS = (20600438, 20600439)  # 20600439: one re-smoke if a parameter name is wrong

SS_DIR = "validation/results/transport/v8-ss"
SEAL_PATH = f"{SS_DIR}/ss-sample-seal.json"
SEAL_SCHEMA = "ionmc-v8-ss-sample-seal-1"
REVIEW_PATH = f"{SS_DIR}/ss-seal-review.json"
HOST_RECORDS_DIR = f"{SS_DIR}/host-records"
TOPAS_RECORDS = f"{SS_DIR}/topas-runs"
SMOKE_RECORDS = f"{SS_DIR}/smoke-runs"
REPORT_PATH = f"{SS_DIR}/ss-report.json"
REPORT_SCHEMA = "ionmc-v8-ss-report-1"
CHECK_SCHEMA = "ionmc-v8-ss-check-run-1"
MANIFEST_PATH = "validation/reference_cases/topas/v5ss-frozen-manifest.json"
MANIFEST_SCHEMA = "ionmc-v5ss-frozen-bundles-1"
NPZ_SUBDIR = "v5ss-samples"
NOTE_KEYS = ("ss_seal_sha", "ss_seal_review_id")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")

# phase space (TO VERIFY on the scorer smoke: column names and order)
EXPECTED_COLUMNS: tuple[str, ...] = (
    "Position X [cm]",
    "Position Y [cm]",
    "Position Z [cm]",
    "Direction Cosine X",
    "Direction Cosine Y",
    "Energy [MeV]",
    "Weight",
    "Particle Type (in PDG Format)",
    "Flag to tell if Third Direction Cosine is Negative (1 means true)",
    "Flag to tell if this is the First Scored Particle from this History (1 means true)",
)
EXIT_Z_CM = 0.1  # far face of the slab in world coordinates (TO VERIFY: world vs local frame)
Z_TOL_CM = 1e-6
E_CUT_MIN_MEV = 1.16  # e- production energy threshold >= T_max of the fastest case
_COLUMN = re.compile(r"^\s*(\d+):\s*(.*?)\s*$")
_ORIGINAL = re.compile(r"^\s*Number of Original Histories:\s*(\d+)\s*$", re.M)
_SCORED = re.compile(r"^\s*Number of Scored Particles:\s*(\d+)\s*$", re.M)
_CS_HEAD = re.compile(r"^CoulombScat:\s+for\s+(\S+)\s", re.M)
_MSC_HEAD = re.compile(r"^msc:\s+for\s+(\S+)", re.M)  # per-particle msc block header
_ETHRESH = re.compile(r"Energy thresholds\s*:.*?\be-\s+([0-9.eE+-]+)\s*(eV|keV|MeV|GeV|TeV)")
_UNIT_MEV = {"eV": 1e-6, "keV": 1e-3, "MeV": 1.0, "GeV": 1e3, "TeV": 1e6}
SLAB_MATERIAL = "G4_WATER"  # the slab material of every V5-SS bundle (see gen_ss_bundles.py)
_COUPLE = re.compile(
    r"^\s*Index\s*:\s*(\d+)\s+used in the geometry\s*:\s*(\w+)\s*\n\s*Material\s*:\s*(\S+)"
    r"[^\n]*\n(?:(?!\s*Index\s*:)[^\n]*\n)*?\s*Energy thresholds\s*:([^\n]*)",
    re.M,
)
_THRESH_ONE = re.compile(r"(gamma|e-|e\+|proton)\s+([0-9.eE+-]+)\s*(eV|keV|MeV|GeV|TeV)")


# =============================================================================================
# names, seeds
# =============================================================================================


def topas_seed(case: int, part: int) -> int:
    return TOPAS_BASE + 100 * part + case


def ionmc_seed(case: int, k: int = 0) -> int:
    return SS_SEED_BASE + 100 * SS_R_INDEX[case] + k


def bundle_name(case: int, part: int) -> str:
    ion, e_u, _pdg, _p = CASE_INFO[case]
    return f"v5ss-{ion}-water-{e_u}mevu-case{case:02d}-part{part}"


EVIDENCE_BUNDLES = {bundle_name(c, p): (c, p) for c in SS_CASES for p in PARTS}


def parse_bundle(bundle: str) -> dict[str, Any]:
    """{"kind": "part"|"smoke", "case", "part"} of one of the 16 bundles or the scorer smoke."""
    if bundle in EVIDENCE_BUNDLES:
        case, part = EVIDENCE_BUNDLES[bundle]
        return {"kind": "part", "case": case, "part": part}
    if bundle == SMOKE_BUNDLE:
        return {"kind": "smoke", "case": 28, "part": None}
    raise RecordError(f"{bundle!r} is not a V5-SS bundle")


# =============================================================================================
# statistic and bootstrap
# =============================================================================================


def contributions(tx: np.ndarray, ty: np.ndarray, window: float) -> tuple[np.ndarray, np.ndarray]:
    """Per primary: s = tx^2 1(|tx|<=W) + ty^2 1(|ty|<=W) and c = 1(|tx|<=W) + 1(|ty|<=W)."""
    tx, ty = np.asarray(tx, np.float64), np.asarray(ty, np.float64)
    ix, iy = np.abs(tx) <= window, np.abs(ty) <= window
    return np.where(ix, tx * tx, 0.0) + np.where(iy, ty * ty, 0.0), ix.astype(np.int64) + iy


def windowed_rms(tx: np.ndarray, ty: np.ndarray, window: float | None) -> tuple[float, int]:
    """Pooled (x and y) rms about zero in |theta| <= window (None: untruncated); (rms, n used)."""
    w = math.inf if window is None else window
    s, c = contributions(tx, ty, w)
    n = int(c.sum())
    if n == 0:
        raise ValueError("no projection inside the window")
    return float(math.sqrt(s.sum() / n)), n


def boot_rms(
    s: np.ndarray, c: np.ndarray, strata: list[int] | None, n_boot: int, rng: np.random.Generator
) -> np.ndarray:
    """``n_boot`` bootstrap replicates of sqrt(sum s / sum c), resampling primaries with
    replacement within each stratum (``strata``: consecutive stratum sizes; None = one)."""
    sizes = strata or [len(s)]
    if sum(sizes) != len(s):
        raise ValueError("strata do not partition the sample")
    out = np.empty(n_boot)
    starts = np.cumsum([0, *sizes[:-1]])
    for b in range(n_boot):
        ss = cc = 0.0
        for st, n in zip(starts, sizes, strict=True):
            idx = st + rng.integers(0, n, size=n)
            ss += float(s[idx].sum())
            cc += float(c[idx].sum())
        out[b] = math.sqrt(ss / cc)
    return out


def pair_dev(model: dict[int, float], ref: dict[int, float], a: int, b: int) -> float:
    return (model[a] / model[b]) / (ref[a] / ref[b]) - 1.0


def grade_pairs(
    ref: dict[int, tuple[list[np.ndarray], list[np.ndarray]]],
    model: dict[int, tuple[np.ndarray, np.ndarray]],
    windows: dict[int, float],
    n_boot: int = N_BOOT,
    seed: int = ANALYSIS_SEED,
) -> dict[str, Any]:
    """Statistics of the five pairs. ``ref[c]`` = (tx per part, ty per part) of TOPAS,
    ``model[c]`` = (tx, ty) of ionmc, ``windows[c]`` = W_c. Returns per-case widths and per-pair
    dev, sigma, gate verdict and the row verdict."""
    cases = sorted({x for _, a, b, _ in PAIRS for x in (a, b)})
    kids = np.random.SeedSequence(seed).spawn(2 * len(cases))
    theta_ref, theta_model, rep_ref, rep_model, per_case = {}, {}, {}, {}, {}
    for i, c in enumerate(cases):
        w = windows[c]
        txs, tys = ref[c]
        tx, ty = np.concatenate(txs), np.concatenate(tys)
        s, n = contributions(tx, ty, w)
        theta_ref[c] = math.sqrt(s.sum() / n.sum())
        rep_ref[c] = boot_rms(
            s, n, [len(x) for x in txs], n_boot, np.random.default_rng(kids[2 * i])
        )
        ms, mn = contributions(*model[c], w)
        theta_model[c] = math.sqrt(ms.sum() / mn.sum())
        rep_model[c] = boot_rms(ms, mn, None, n_boot, np.random.default_rng(kids[2 * i + 1]))
        per_case[c] = {
            "window_rad": w,
            "theta_ref_rad": theta_ref[c],
            "theta_model_rad": theta_model[c],
            "n_ref_primaries": int(len(tx)),
            "n_model_histories": int(len(model[c][0])),
            "p_window_ref": float(n.sum() / (2 * len(tx))),
            "p_window_model": float(mn.sum() / (2 * len(model[c][0]))),
        }
    pairs = []
    for name, a, b, expected in PAIRS:
        dev = pair_dev(theta_model, theta_ref, a, b)
        reps = (rep_model[a] / rep_model[b]) / (rep_ref[a] / rep_ref[b]) - 1.0
        sigma = float(np.std(reps, ddof=1))
        verdict = ss_budget.gate(100.0 * dev, 100.0 * sigma)
        pairs.append(
            {
                "pair": name,
                "case_a": a,
                "case_b": b,
                "dev_pct": 100.0 * dev,
                "sigma_pct": 100.0 * sigma,
                "dev_plus_1645_sigma_pct": 100.0 * (abs(dev) + 1.645 * sigma),
                "verdict": verdict,
                "expected_dev_pct": expected,
            }
        )
    return {
        "cases": per_case,
        "pairs": pairs,
        "row_verdict": row_verdict(pairs),
        "n_boot": n_boot,
        "analysis_seed": seed,
    }


def row_verdict(pairs: list[dict[str, Any]]) -> str:
    """PASS iff every gated pair PASSes (INCONCLUSIVE and FAIL both fail the row)."""
    if len(pairs) != len(PAIRS):
        return "FAIL"
    return "PASS" if all(p["verdict"] == "PASS" for p in pairs) else "FAIL"


# =============================================================================================
# phase space (ASCII)
# =============================================================================================


def parse_header(text: str) -> dict[str, Any]:
    """Columns and counts of a TOPAS ASCII phase-space ``.header``. Fails closed on a column list
    other than :data:`EXPECTED_COLUMNS` (the error names the observed list) or missing counts."""
    cols = []
    inside = False
    for ln in text.splitlines():
        if ln.strip().startswith("Columns of data are as follows"):
            inside = True
            continue
        if inside:
            m = _COLUMN.match(ln)
            if m:
                cols.append((int(m.group(1)), m.group(2)))
            elif cols and not ln.strip():
                break
    names = [n for _, n in cols]
    if [i for i, _ in cols] != list(range(1, len(cols) + 1)) or tuple(names) != EXPECTED_COLUMNS:
        raise RecordError(
            f"phase-space header: observed column list {names} != expected {list(EXPECTED_COLUMNS)}"
        )
    orig, scored = _ORIGINAL.findall(text), _SCORED.findall(text)
    if len(orig) != 1 or len(scored) != 1:
        raise RecordError("phase-space header lacks the original-history or scored-particle count")
    return {"columns": names, "original_histories": int(orig[0]), "scored": int(scored[0])}


def parse_phasespace(header: str, rows: str, histories: int, pdg: int) -> dict[str, Any]:
    """Check one ASCII phase space (every primary scored exactly once at the far face) and return
    ``{"theta_x", "theta_y", "n"}`` (radians, float64). Problems raise RecordError."""
    h = parse_header(header)
    if h["original_histories"] != histories or h["scored"] != histories:
        raise RecordError(
            f"phase space: original histories {h['original_histories']} / scored {h['scored']} "
            f"!= {histories}"
        )
    try:
        data = np.loadtxt(_io.StringIO(rows), dtype=np.float64, ndmin=2)
    except ValueError as exc:
        raise RecordError(f"phase space rows are not numeric: {exc}") from exc
    if data.shape != (histories, len(EXPECTED_COLUMNS)):
        raise RecordError(f"phase space has shape {data.shape}, not ({histories}, 10)")
    if not np.all(np.isfinite(data)):
        raise RecordError("phase space holds non-finite values")
    z, u, v = data[:, 2], data[:, 3], data[:, 4]
    if not np.all(data[:, 7] == pdg):
        raise RecordError(f"phase space holds particles other than PDG {pdg}")
    if not np.all(data[:, 8] == 0):
        raise RecordError("phase space: a particle with negative third direction cosine")
    if not np.all(data[:, 9] == 1):
        raise RecordError("phase space: a row is not the first scored particle of its history")
    if not np.all(np.abs(z - EXIT_Z_CM) <= Z_TOL_CM):
        raise RecordError(f"phase space: Z differs from {EXIT_Z_CM} cm")
    w2 = 1.0 - u * u - v * v
    if not np.all(w2 > 0.0):
        raise RecordError("phase space: direction cosines with u^2 + v^2 >= 1")
    w = np.sqrt(w2)
    return {"theta_x": np.arctan2(u, w), "theta_y": np.arctan2(v, w), "n": histories}


# =============================================================================================
# plan note, review, manifest
# =============================================================================================


def parse_ss_note(text: str) -> dict[str, Any]:
    """The SS seal note: the one fenced ```json block containing ``"ss_seal_sha"`` with exactly
    the keys :data:`NOTE_KEYS`; returned normalized as ``{"seal_sha", "seal_review_id"}``."""
    blocks = [b for b in re.findall(r"```json[ \t]*\n(.*?)```", text, re.S) if '"ss_seal_sha"' in b]
    if len(blocks) != 1:
        raise RecordError(
            f"{io.PLAN_PATH}: expected exactly one SS seal-note block, found {len(blocks)}"
        )
    try:
        note = json.loads(blocks[0])
    except ValueError as exc:
        raise RecordError(f"{io.PLAN_PATH}: the SS seal-note block is not JSON: {exc}") from exc
    if not isinstance(note, dict) or sorted(note) != sorted(NOTE_KEYS):
        raise RecordError(
            f"{io.PLAN_PATH}: the SS seal note must hold exactly the keys {NOTE_KEYS}"
        )
    if not isinstance(note["ss_seal_sha"], str) or not SHA_RE.match(note["ss_seal_sha"]):
        raise RecordError("SS seal note: ss_seal_sha is not a 40-hex commit SHA")
    return {"seal_sha": note["ss_seal_sha"], "seal_review_id": note["ss_seal_review_id"]}


def load_note_and_review(reader: io.CommittedReader) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = parse_ss_note(reader.read(io.PLAN_PATH).decode("utf-8"))
    review = io.check_review(reader.json(REVIEW_PATH, "SS seal review record"), plan)
    return plan, review


def pinned_chronology(repo: Path, code_sha: str) -> dict[str, Any]:
    """SS seal note and review record as committed at the run's signed ``code_sha``."""
    plan_blob = io.git_show(repo, code_sha, io.PLAN_PATH)
    plan = parse_ss_note(plan_blob.decode("utf-8"))
    rev_blob = io.git_show(repo, code_sha, REVIEW_PATH)
    try:
        review = json.loads(rev_blob.decode("utf-8"))
    except ValueError as exc:
        raise RecordError(f"{REVIEW_PATH} at {code_sha} is not JSON: {exc}") from exc
    io.check_review(review, plan)
    return {
        "plan": plan,
        "review": review,
        "plan_sha256": hashlib.sha256(plan_blob).hexdigest(),
        "review_sha256": hashlib.sha256(rev_blob).hexdigest(),
    }


def ss_manifest(reader: io.CommittedReader, code_sha: str, bundle: str, smoke: bool) -> dict:
    """The SS frozen manifest at HEAD bound to the run's ``code_sha``: byte-identical (evidence)
    or, for the scorer smoke, the smoke's own entry identical."""
    at_run = io.git_show(reader.repo, code_sha, MANIFEST_PATH)
    at_head = reader.read(MANIFEST_PATH)
    docs = []
    for blob, where in ((at_head, "HEAD"), (at_run, f"code_sha {code_sha}")):
        try:
            doc = json.loads(blob.decode("utf-8"))
        except ValueError as exc:
            raise RecordError(f"{MANIFEST_PATH} at {where} is not JSON: {exc}") from exc
        if not isinstance(doc, dict) or doc.get("schema") != MANIFEST_SCHEMA:
            raise RecordError(f"{MANIFEST_PATH} at {where}: schema is not {MANIFEST_SCHEMA}")
        docs.append(doc)
    head, run = docs
    if not smoke and at_run != at_head:
        raise RecordError(
            f"{MANIFEST_PATH} differs between the run's code_sha {code_sha} and HEAD {reader.head}"
        )
    old, new = run["bundles"].get(bundle), head["bundles"].get(bundle)
    if not isinstance(old, dict) or old != new:
        raise RecordError(f"bundle {bundle}: SS manifest entry missing or changed after the run")
    return head


# =============================================================================================
# signed TOPAS runs
# =============================================================================================


def load_run(
    repo: Path,
    reader: io.CommittedReader,
    rdir: Path,
    outputs: Path,
    run_id: str,
    bundle: str,
    trust: dict[str, Any],
) -> io.TopasRun:
    """Load and bind one signed reference-runner record of a V5-SS bundle (``r2_io`` record
    rules: verified attestation, schema 3, run id from configuration id, success, listing,
    output_files_sha256, request binding to the bundle committed at ``code_sha``, SS manifest)."""
    info = parse_bundle(bundle)
    record_path = rdir / "record.json"
    result = reader.json(record_path, f"signed run record {run_id}")
    status, why = io.attest_module().classify_record(result, trust, io.REFERENCE_SIGNER)
    if status != "verified":
        raise RecordError(f"{run_id}: record.json attestation {status}: {why}")
    for key in io.RESULT_FIELDS:
        io.require(result, key, f"result record {run_id}")
    if result["schema_version"] != 3 or result.get("timed_out") is not False:
        raise RecordError(f"{run_id}: schema_version/timed_out not 3/false")
    if (
        not re.fullmatch(
            r"REF-" + re.escape(str(result["configuration_id"])[:20]) + r"-[0-9a-f]{8}", run_id
        )
        or result["run_id"] != run_id
    ):
        raise RecordError(f"{run_id}: run id not bound to the record/configuration")
    if result["engine"] != "topas" or not (
        result["succeeded"] is True and result["exit_code"] == 0 and result["missing_outputs"] == []
    ):
        raise RecordError(f"{run_id}: not a successful topas run")
    started = io._parse_time(result["started_at"], f"{run_id} started_at")
    if started > io._parse_time(result["finished_at"], f"{run_id} finished_at"):
        raise RecordError(f"{run_id}: started_at is later than finished_at")
    code_sha = str(result["code_sha"])
    sizes: dict[str, Any] = {}
    listing = io._listing(
        reader.json(rdir / "artifacts.json", "artifact listing"), run_id, code_sha, sizes
    )
    work = {}
    for p, d in listing.items():
        if p.startswith("work/"):
            b = sizes.get(p)
            if not isinstance(b, int) or isinstance(b, bool) or b < 0:
                raise RecordError(f"{run_id}: artifact listing lacks the byte size of {p}")
            work[p[len("work/") :]] = {"sha256": d, "bytes": b}
    if io.canonical_digest(work) != result["output_files_sha256"] or result["output_count"] != len(
        work
    ):
        raise RecordError(f"{run_id}: listed work/ outputs do not reproduce the signed digest")
    for name in ("request.json", "stdout.txt", "stderr.txt", "record.json"):
        if name not in listing:
            raise RecordError(f"{run_id}: artifact listing lacks {name}")
    if reader.reads[reader.rel(record_path)] != listing["record.json"]:
        raise RecordError(f"{run_id}: record.json differs from the listed record.json")
    for name, key in (("stdout.txt", "stdout_sha256"), ("stderr.txt", "stderr_sha256")):
        if result[key] != listing[name]:
            raise RecordError(f"{run_id}: signed {key} differs from the listing")
    run = io.TopasRun(
        run_id=run_id,
        case=info["case"],
        config="v5ss",
        part=info["part"],
        bundle=bundle,
        configuration_id=str(result["configuration_id"]),
        code_sha=code_sha,
        finished_at=str(result["finished_at"]),
        seconds=float(result["seconds"]),
        reference_started=str(result["started_at"]),
        histories=0,
        expected_outputs=[],
        listing=listing,
        outputs_dir=outputs / run_id,
    )
    for name in ("stdout.txt", "stderr.txt"):
        io.read_verified(run, name)
    request_blob = io.read_verified(run, "request.json")
    if hashlib.sha256(request_blob).hexdigest() != result["request_sha256"]:
        raise RecordError(f"{run_id}: request.json does not hash to the signed request_sha256")
    request = json.loads(request_blob.decode("utf-8"))
    if io.canonical_digest(request) != run.configuration_id or request.get("code_sha") != code_sha:
        raise RecordError(f"{run_id}: request.json not bound to the record")
    bdir = f"{io.BUNDLE_ROOT}/{bundle}"
    committed = {}
    for rel in io.git_ls(repo, code_sha, bdir):
        blob = io.git_show(repo, code_sha, rel)
        committed[rel[len(bdir) + 1 :]] = {
            "sha256": hashlib.sha256(blob).hexdigest(),
            "bytes": len(blob),
        }
    if not committed or request.get("input_files") != committed:
        raise RecordError(f"{run_id}: request input_files differ from the bundle at code_sha")
    if not io.is_ancestor(repo, code_sha, reader.head):
        raise RecordError(f"{run_id}: code_sha {code_sha} is not an ancestor of HEAD")
    manifest = ss_manifest(reader, code_sha, bundle, info["kind"] == "smoke")
    entry = manifest["bundles"][bundle]
    for name in ("case.json", "input.txt"):
        if committed.get(name, {}).get("sha256") != entry.get(name):
            raise RecordError(f"{run_id}: {bundle}/{name} differs from the SS frozen manifest")
    case_doc = json.loads(io.git_show(repo, code_sha, f"{bdir}/case.json").decode("utf-8"))
    if request.get("case") != case_doc:
        raise RecordError(f"{run_id}: request case differs from the committed case.json")
    seeds = case_doc.get("seeds")
    if info["kind"] == "part":
        if seeds != [topas_seed(info["case"], info["part"])]:
            raise RecordError(f"{run_id}: seed {seeds} breaks the formula")
        if case_doc.get("histories") != PART_HISTORIES:
            raise RecordError(f"{run_id}: histories {case_doc.get('histories')} != 30000")
    elif not (isinstance(seeds, list) and len(seeds) == 1 and seeds[0] in SMOKE_SEEDS):
        raise RecordError(f"{run_id}: smoke seed {seeds} is not a declared smoke seed")
    run.histories = int(case_doc["histories"])
    run.case_doc, run.request = case_doc, request
    run.native = io.git_show(repo, code_sha, f"{bdir}/input.txt").decode("utf-8", "replace")
    run.pinned = pinned_chronology(repo, code_sha)
    run.pinned["manifest_bundle_set_digest"] = manifest.get("bundle_set_digest")
    return run


def chronology(repo: Path, run: io.TopasRun, plan: dict[str, Any]) -> dict[str, Any]:
    """SS-SEAL-1 (a) ancestry and (b) review-before-start, pinned at the run's code_sha."""
    if run.pinned["plan"] != plan:
        raise RecordError(f"{run.run_id}: the SS seal note at code_sha differs from HEAD's")
    v = r2e.check_run_validity(
        repo, run.pinned["plan"], run.pinned["review"], run.chronology_record(), run.run_id
    )
    return {
        **v,
        "review_finished_at": run.pinned["review"]["finished_at"],
        "review_sha256_at_code_sha": run.pinned["review_sha256"],
    }


def coulomb_models(stdout: str, particle: str) -> list[str]:
    models: list[str] = []
    for m in _CS_HEAD.finditer(stdout):
        if m.group(1) == particle:
            models.extend(io._MODEL.findall(stdout[m.end() :].split("\n\n", 1)[0]))
    return models


def electron_threshold_mev(stdout: str) -> list[float]:
    """e- production energy thresholds of the couple table (TO VERIFY: string form)."""
    return [float(v) * _UNIT_MEV[u] for v, u in _ETHRESH.findall(stdout)]


def couples(stdout: str) -> list[dict[str, Any]]:
    """Couples of the Geant4 "Table of registered couples": index, material, used_in_geometry and
    the production energy thresholds in MeV (gamma, e-, e+, proton; absent if not printed)."""
    out: list[dict[str, Any]] = []
    for m in _COUPLE.finditer(stdout):
        thr = {n: float(v) * _UNIT_MEV[u] for n, v, u in _THRESH_ONE.findall(m.group(4))}
        out.append(
            {
                "index": int(m.group(1)),
                "material": m.group(3),
                "used_in_geometry": m.group(2) == "Yes",
                "gamma_threshold_mev": thr.get("gamma"),
                "e_minus_threshold_mev": thr.get("e-"),
                "e_plus_threshold_mev": thr.get("e+"),
                "proton_threshold_mev": thr.get("proton"),
            }
        )
    return out


def msc_info(stdout: str) -> dict[str, Any]:
    """Multiple-scattering facts of a log: particles that have an ``msc:  for <particle>`` block and
    whether the global "Registered Physics Processes" table lists ``msc`` (that table names
    processes only, not the particle, so a listed ``msc`` is informational: the physics list may
    attach it to e-/e+ while the beam particle uses single scattering only)."""
    names = io.registered_processes(stdout)
    return {
        "msc_blocks_for": sorted(set(_MSC_HEAD.findall(stdout))),
        "registered_msc_entry": None if names is None else "msc" in names,
    }


def log_checks(stdout: str, stderr: bytes, histories: int, particle: str) -> list[str]:
    """Run-validity problems of a TOPAS log (empty = valid). Single-scattering rule is per beam
    particle: fail if an ``msc:  for <beam particle>`` block exists, if ``CoulombScat`` is not in
    the registered list, or if the beam particle's CoulombScat models are not exactly
    ``IonCoulombScattering``. An msc block for another particle (e-, e+) or a global ``msc`` entry
    in the registered table is recorded by :func:`msc_info`, not a failure."""
    problems = []
    if stderr.strip():
        problems.append("stderr is not empty")
    bad = [ln for ln in stdout.splitlines() if io._KILLED.search(ln) and not io._BEGIN.search(ln)]
    if bad:
        problems.append(f"stdout reports killed/stuck tracks ({len(bad)} line(s))")
    found = [int(v) for v in io._HISTORIES.findall(stdout)]
    if not found or any(v != histories for v in found):
        problems.append(f"histories completed {found} != {histories}")
    models = coulomb_models(stdout, particle)
    if models != ["IonCoulombScattering"]:
        problems.append(f"CoulombScat models for {particle} are {models}, not IonCoulombScattering")
    if particle in msc_info(stdout)["msc_blocks_for"]:
        problems.append(f"an msc process block is printed for the beam particle {particle}")
    names = io.registered_processes(stdout)
    if names is None or "CoulombScat" not in names:
        problems.append("process list absent or lacks CoulombScat")
    slab = [c for c in couples(stdout) if c["material"] == SLAB_MATERIAL]
    if len(slab) != 1:
        problems.append(f"{len(slab)} couples with material {SLAB_MATERIAL} (exactly one required)")
    else:
        if not slab[0]["used_in_geometry"]:
            problems.append(f"the {SLAB_MATERIAL} couple is not used in the geometry")
        thr = slab[0]["e_minus_threshold_mev"]
        if thr is None or thr < E_CUT_MIN_MEV:
            problems.append(
                f"{SLAB_MATERIAL} e- energy threshold {thr} MeV is not >= {E_CUT_MIN_MEV} MeV"
            )
    return problems


def run_angles(run: io.TopasRun) -> dict[str, Any]:
    """Validated exit angles of one run's phase space (``work/ss_exit.*``)."""
    header = io.read_verified(run, "work/ss_exit.header").decode("utf-8", "replace")
    rows = io.read_verified(run, "work/ss_exit.phsp").decode("utf-8", "replace")
    return parse_phasespace(header, rows, run.histories, CASE_INFO[run.case][2])


def check_one(
    repo: Path,
    reader: io.CommittedReader,
    plan: dict[str, Any],
    run: io.TopasRun,
    read_angles: bool = True,
) -> dict[str, Any]:
    chron = chronology(repo, run, plan)
    stdout = io.read_verified(run, "stdout.txt").decode("utf-8", "replace")
    stderr = io.read_verified(run, "stderr.txt")
    problems = log_checks(stdout, stderr, run.histories, CASE_INFO[run.case][3])
    angles = None
    if read_angles:
        try:
            angles = run_angles(run)
        except RecordError as exc:
            problems.append(f"phase space: {exc}")
    verdict = (
        "BLINDING FAIL" if not chron["valid"] else "RUN-VALIDITY FAIL" if problems else "VALID"
    )
    return {
        "verdict": verdict,
        "chronology": chron,
        "problems": problems,
        "angles": angles,
        **msc_info(stdout),
        "couples": couples(stdout),
    }


def check_run(
    run_id: str,
    bundle: str,
    records_dir: Path | None = None,
    outputs: Path | None = None,
    *,
    repo: Path | None = None,
) -> dict[str, Any]:
    """Run validity of ONE signed V5-SS TOPAS run (evidence part or scorer smoke); nothing graded.
    VALID / BLINDING FAIL / RUN-VALIDITY FAIL; record gaps raise (ERROR)."""
    repo = Path(repo).resolve() if repo is not None else r2e.evaluator_repo()
    reader = io.CommittedReader(repo)
    info = parse_bundle(bundle)
    if records_dir is None:
        records_dir = Path(TOPAS_RECORDS if info["kind"] == "part" else SMOKE_RECORDS)
    records_dir = records_dir if records_dir.is_absolute() else repo / records_dir
    plan, review = load_note_and_review(reader)
    trust, trust_info = io.committed_trust(reader)
    outputs = outputs or repo / ".ionmc-cache" / "reference-runs"
    run = load_run(repo, reader, records_dir / run_id, outputs, run_id, bundle, trust)
    res = check_one(repo, reader, plan, run)
    return {
        "schema": CHECK_SCHEMA,
        "run_id": run_id,
        "bundle": bundle,
        "bundle_kind": info["kind"],
        "case": info["case"],
        "part": info["part"],
        "verdict": res["verdict"],
        "problems": res["problems"],
        "msc_blocks_for": res["msc_blocks_for"],
        "registered_msc_entry": res["registered_msc_entry"],
        "chronology": res["chronology"],
        "ss_seal_sha": plan["seal_sha"],
        "ss_seal_review_id": review["review_id"],
        "code_sha": run.code_sha,
        "reference_started": run.reference_started,
        "seconds": run.seconds,
        "histories": run.histories,
        "seeds": run.case_doc.get("seeds"),
        "files_sha256": dict(sorted(run.read.items())),
        "trust_file": trust_info,
        "evaluator": {
            "head_sha": reader.head,
            "committed_inputs": dict(sorted(reader.reads.items())),
        },
    }


# =============================================================================================
# ionmc samples: host-run records and the seal
# =============================================================================================


def _exact(got: Any, want: Any) -> bool:
    return io._exact(got, want)


def case_entry(case: int, doc: dict[str, Any], npz_path: Path, impl: str) -> dict[str, Any]:
    """Seal entry of one case from its step document and npz (fails closed on any mismatch)."""
    seed = ionmc_seed(case)
    expected = {
        "row": "v5ss",
        "case": case,
        "r_index": SS_R_INDEX[case],
        "k": 0,
        "seed": seed,
        "seed_base": SS_SEED_BASE,
        "histories": SS_HISTORIES,
        "frozen_histories": SS_HISTORIES,
        "reduced": False,
        "backend": "warp-cpu",
        "precision": "float64",
        "valid": True,
        "pass": True,
        "git_sha": impl,
        "suite": "v5ss",
        "n_escaped": SS_HISTORIES,
        "npz": f"ionmc-ss-case{case}-k0.npz",
    }
    for key, want in expected.items():
        if key not in doc or not _exact(doc[key], want):
            raise RecordError(f"case {case}: step document {key} {doc.get(key)!r} != {want!r}")
    counters = doc.get("counters")
    if (
        not isinstance(counters, dict)
        or not counters
        or not all(_exact(v, 0) for v in counters.values())
    ):
        raise RecordError(f"case {case}: counters are not all zero")
    eff = doc.get("effective_config")
    if not isinstance(eff, dict) or not all(
        _exact(eff.get(k), v)
        for k, v in {"cpu_workers": 1, "seed": seed, "n_histories": SS_HISTORIES}.items()
    ):
        raise RecordError(f"case {case}: effective_config is not the frozen run")
    if npz_path.is_symlink() or not npz_path.is_file():
        raise RecordError(f"case {case}: npz {npz_path} missing")
    blob = npz_path.read_bytes()
    digest = hashlib.sha256(blob).hexdigest()
    if digest != doc.get("npz_sha256"):
        raise RecordError(f"case {case}: npz sha256 differs from the step document")
    with np.load(_io.BytesIO(blob)) as z:
        tx, ty = np.asarray(z["theta_x"]), np.asarray(z["theta_y"])
        hist, zexit = np.asarray(z["escape_history"]), np.asarray(z["z_exit_mm"])
        stored = float(z["theta_rms_untruncated"])
    if tx.shape != (SS_HISTORIES,) or ty.shape != (SS_HISTORIES,) or tx.dtype != np.float64:
        raise RecordError(f"case {case}: npz angles are not {SS_HISTORIES} float64 values")
    if not np.array_equal(hist, np.arange(SS_HISTORIES)) or not np.all(np.abs(zexit - 1.0) <= 1e-9):
        raise RecordError(f"case {case}: npz escape records are not one far-face exit per history")
    anchor, _ = windowed_rms(tx, ty, None)
    if not math.isclose(anchor, stored, rel_tol=1e-12) or not math.isclose(
        anchor, float(doc.get("theta_anchor_rad", math.nan)), rel_tol=1e-12
    ):
        raise RecordError(f"case {case}: theta_anchor differs between npz, document and recompute")
    window = WINDOW_K * anchor
    theta_w, n_in = windowed_rms(tx, ty, window)
    return {
        "case": case,
        "seed": seed,
        "r_index": SS_R_INDEX[case],
        "k": 0,
        "histories": SS_HISTORIES,
        "npz": npz_path.name,
        "npz_sha256": digest,
        "npz_bytes": len(blob),
        "theta_anchor_rad": anchor,
        "window_k": WINDOW_K,
        "window_rad": window,
        "theta_w_rad": theta_w,
        "n_in_window": n_in,
        "p_window": n_in / (2.0 * SS_HISTORIES),
    }


def verify_sample_run(
    rec: Any, trust: dict[str, Any], impl: str, samples_root: Path, record_sha256: str
) -> dict[str, Any]:
    """Bind one signed host-run record of an ionmc sample: verified ``host-runner`` attestation,
    clean success at ``impl``, argv == ["python", "-c", render_wrapper(case, impl, NPZ_SUBDIR,
    "1")], published outputs reproduce the signed manifest and are exactly the npz, stdout and
    stderr; the step document binds the case (:func:`case_entry`). Returns the seal entry."""
    if not isinstance(rec, dict):
        raise RecordError("host-run record is not a JSON object")
    run_id = io.require(rec, "run_id", "host-run record")
    status, why = io.attest_module().classify_record(rec, trust, io.HOST_SIGNER)
    if status != "verified":
        raise RecordError(f"host run {run_id}: attestation {status}: {why}")
    io._bind(
        run_id,
        "signed",
        rec,
        {
            "sha": impl,
            "task_id": io.HOST_TASK_ID,
            "succeeded": True,
            "exit_code": 0,
            "timed_out": False,
            "dirty_before": False,
            "dirty_after": False,
            "committed_state_unchanged": True,
            "output_publication_error": None,
        },
    )
    argv = rec.get("argv")
    cases = [
        c
        for c in SS_CASES
        if argv == ["python", "-c", render_wrapper(c, impl, NPZ_SUBDIR, SS_SCALE)]
    ]
    if len(cases) != 1:
        raise RecordError(f"host run {run_id}: argv is not the reviewed SS wrapper for any case")
    case = cases[0]
    digest, files = io.published_manifest(samples_root, run_id)
    if digest != rec.get("output_manifest_sha256"):
        raise RecordError(f"host run {run_id}: published outputs do not reproduce the manifest")
    base = f"validation/generated/{NPZ_SUBDIR}/case{case}"
    names = (
        f"ionmc-ss-case{case}-k0.npz",
        f"case{case}-k0.stdout.txt",
        f"case{case}-k0.stderr.txt",
    )
    if set(files) != {f"{base}/{n}" for n in names}:
        raise RecordError(f"host run {run_id}: published outputs are not exactly {names}")
    top = Path(samples_root) / "validation" / "generated" / run_id / NPZ_SUBDIR / f"case{case}"
    stdout_blob = (top / names[1]).read_bytes()
    doc = io.step_document(stdout_blob.decode("utf-8"), run_id)
    entry = case_entry(case, doc, top / names[0], impl)
    if files[f"{base}/{names[0]}"] != entry["npz_sha256"]:
        raise RecordError(f"host run {run_id}: manifest npz digest differs")
    att = rec["attestation"]
    return {
        **entry,
        "run_id": run_id,
        "run_started_at": rec.get("started_at"),
        "run_finished_at": rec.get("finished_at"),
        "output_manifest_sha256": digest,
        "wrapper_stdout_sha256": rec.get("stdout_sha256"),
        "step_stdout_sha256": hashlib.sha256(stdout_blob).hexdigest(),
        "host_record_sha256": record_sha256,
        "wall_s": doc.get("wall_s"),
        "host_runner_attestation": att,
    }


def build_seal(entries: list[dict[str, Any]], impl: str) -> dict[str, Any]:
    """The seal record from one verified entry per case (all eight, each once)."""
    if not isinstance(impl, str) or not SHA_RE.match(impl):
        raise RecordError(f"implementation_sha {impl!r} is not a 40-hex SHA")
    by_case: dict[int, dict[str, Any]] = {}
    for e in entries:
        if e["case"] in by_case:
            raise RecordError(f"case {e['case']}: more than one sample run")
        if e["seed"] != ionmc_seed(e["case"]) or e["histories"] != SS_HISTORIES:
            raise RecordError(f"case {e['case']}: seed/histories break the declared design")
        if not math.isclose(e["window_rad"], WINDOW_K * e["theta_anchor_rad"], rel_tol=1e-15):
            raise RecordError(f"case {e['case']}: W != {WINDOW_K} theta_anchor")
        by_case[e["case"]] = e
    if sorted(by_case) != sorted(SS_CASES):
        raise RecordError(f"seal needs exactly the cases {SS_CASES}, got {sorted(by_case)}")
    return {
        "schema": SEAL_SCHEMA,
        "task_id": "V3-008",
        "row": "V5-SS",
        "implementation_sha": impl,
        "seed_base": SS_SEED_BASE,
        "seed_formula": "20631004 + 100 * r_index + k, k = 0",
        "histories": SS_HISTORIES,
        "window_k": WINDOW_K,
        "npz_subdir": NPZ_SUBDIR,
        "wrapper_template_sha256": SS_WRAPPER_SHA256,
        "backend": "warp-cpu",
        "precision": "float64",
        "single_process": True,
        "cases": [by_case[c] for c in SS_CASES],
        "host_runs": [
            {
                "run_id": by_case[c]["run_id"],
                "output_manifest_sha256": by_case[c]["output_manifest_sha256"],
            }
            for c in SS_CASES
        ],
    }


def seal(impl: str, host_records: Path, samples_root: Path, repo: Path | None = None) -> dict:
    repo = Path(repo).resolve() if repo is not None else r2e.evaluator_repo()
    reader = io.CommittedReader(repo)
    trust, _ = io.committed_trust(reader)
    if not io.is_ancestor(repo, impl, reader.head):
        raise RecordError(f"implementation_sha {impl} is not an ancestor of HEAD")
    entries = []
    for path in sorted(Path(host_records).glob("RUN-*.json")):
        blob = path.read_bytes()
        rec = json.loads(blob.decode("utf-8"))
        if path.stem != rec.get("run_id"):
            raise RecordError(f"{path.name} does not hold run {path.stem}")
        entries.append(
            verify_sample_run(rec, trust, impl, samples_root, hashlib.sha256(blob).hexdigest())
        )
    return build_seal(entries, impl)


def load_sealed_samples(
    repo: Path,
    reader: io.CommittedReader,
    plan: dict[str, Any],
    trust: dict[str, Any],
    samples_root: Path,
) -> tuple[dict, dict]:
    """The committed seal (byte-identical to the one at ss_seal_sha), every case re-verified
    against its committed signed host-run record; returns (seal, {case: (tx, ty)})."""
    blob = reader.read(SEAL_PATH)
    if blob != io.git_show(repo, plan["seal_sha"], SEAL_PATH):
        raise RecordError(f"{SEAL_PATH} differs from the one committed at the SS seal")
    sealed = json.loads(blob.decode("utf-8"))
    if sealed.get("schema") != SEAL_SCHEMA:
        raise RecordError("SS seal schema mismatch")
    impl = sealed.get("implementation_sha")
    if (
        not isinstance(impl, str)
        or not SHA_RE.match(impl)
        or not io.is_ancestor(repo, impl, plan["seal_sha"])
    ):
        raise RecordError("SS seal implementation_sha is not an ancestor of the seal commit")
    entries, samples = [], {}
    for c in sealed.get("cases", []):
        rel = f"{HOST_RECORDS_DIR}/{c.get('run_id')}.json"
        rec = json.loads(reader.read(rel).decode("utf-8"))
        e = verify_sample_run(rec, trust, impl, samples_root, reader.reads[rel])
        if e != c:
            raise RecordError(f"case {c.get('case')}: re-verified entry differs from the seal")
        entries.append(e)
        top = samples_root / "validation" / "generated" / e["run_id"] / NPZ_SUBDIR
        with np.load(top / f"case{e['case']}" / e["npz"]) as z:
            samples[e["case"]] = (np.asarray(z["theta_x"]), np.asarray(z["theta_y"]))
    if build_seal(entries, impl) != sealed:
        raise RecordError("the committed seal is not the seal of its verified entries")
    return sealed, samples


# =============================================================================================
# grade
# =============================================================================================


def load_index(reader: io.CommittedReader) -> dict[tuple[int, int], str]:
    idx = reader.json(f"{TOPAS_RECORDS}/index.json", "V5-SS run index")
    runs = io.require(idx, "runs", "V5-SS run index")
    out: dict[tuple[int, int], str] = {}
    for e in runs:
        key = (io.require(e, "case", "index entry"), io.require(e, "part", "index entry"))
        if key not in {(c, p) for c in SS_CASES for p in PARTS}:
            raise RecordError(f"index entry {key} is not a V5-SS (case, part)")
        if key in out:
            raise RecordError(f"(case, part) {key} listed twice")
        out[key] = str(io.require(e, "run_id", "index entry"))
    if len(out) != 16:
        raise RecordError(f"index lists {len(out)} of 16 parts")
    return out


def _kurtosis_corr(tx: np.ndarray, ty: np.ndarray, w: float) -> dict[str, float]:
    t = np.concatenate([tx[np.abs(tx) <= w], ty[np.abs(ty) <= w]])
    m2 = float(np.mean(t * t))
    both = (np.abs(tx) <= w) & (np.abs(ty) <= w)
    corr = float(np.corrcoef(tx[both] ** 2, ty[both] ** 2)[0, 1])
    return {"window_kurtosis": float(np.mean(t**4)) / (m2 * m2), "xy_corr_sq": corr}


def grade(
    outputs: Path | None = None,
    samples_root: Path | None = None,
    *,
    repo: Path | None = None,
    n_boot: int = N_BOOT,
) -> dict[str, Any]:
    repo = Path(repo).resolve() if repo is not None else r2e.evaluator_repo()
    reader = io.CommittedReader(repo)
    plan, review = load_note_and_review(reader)
    trust, trust_info = io.committed_trust(reader)
    outputs = outputs or repo / ".ionmc-cache" / "reference-runs"
    samples_root = samples_root or repo
    index = load_index(reader)
    runs, checks = {}, {}
    for (case, part), run_id in sorted(index.items()):
        run = load_run(
            repo,
            reader,
            repo / TOPAS_RECORDS / run_id,
            outputs,
            run_id,
            bundle_name(case, part),
            trust,
        )
        runs[(case, part)] = run
        checks[(case, part)] = check_one(repo, reader, plan, run)
    base = {
        "schema": REPORT_SCHEMA,
        "row": "V5-SS",
        "ss_seal_sha": plan["seal_sha"],
        "ss_seal_review_id": review["review_id"],
        "runs": [
            {
                "case": c,
                "part": p,
                "run_id": r.run_id,
                "code_sha": r.code_sha,
                "reference_started": r.reference_started,
                "seconds": r.seconds,
                "verdict": checks[(c, p)]["verdict"],
                "problems": checks[(c, p)]["problems"],
                "files_sha256": dict(sorted(r.read.items())),
            }
            for (c, p), r in sorted(runs.items())
        ],
        "trust_file": trust_info,
    }
    bad = [k for k, v in checks.items() if v["verdict"] != "VALID"]
    if bad:
        blind = any(checks[k]["verdict"] == "BLINDING FAIL" for k in bad)
        return {
            **base,
            "row_verdict": "FAIL",
            "category": "blinding control" if blind else "run validity",
            "invalid_parts": sorted(bad),
            "note": "nothing graded; never re-graded",
        }
    sealed, model = load_sealed_samples(repo, reader, plan, trust, samples_root)
    windows = {c["case"]: c["window_rad"] for c in sealed["cases"]}
    ref = {
        c: (
            [checks[(c, p)]["angles"]["theta_x"] for p in PARTS],
            [checks[(c, p)]["angles"]["theta_y"] for p in PARTS],
        )
        for c in SS_CASES
    }
    res = grade_pairs(ref, model, windows, n_boot=n_boot)
    report_only = {}
    for c in SS_CASES:
        w = windows[c]
        parts = [
            windowed_rms(
                checks[(c, p)]["angles"]["theta_x"], checks[(c, p)]["angles"]["theta_y"], w
            )[0]
            for p in PARTS
        ]
        tx, ty = np.concatenate(ref[c][0]), np.concatenate(ref[c][1])
        report_only[str(c)] = {
            "abs_ratio_model_over_ref_pct": 100.0
            * (res["cases"][c]["theta_model_rad"] / res["cases"][c]["theta_ref_rad"] - 1.0),
            "expected_reference_model_factor_pct": list(EXPECTED_ABS_FACTOR_PCT),
            "theta_ref_part_rad": parts,
            "ref": _kurtosis_corr(tx, ty, w),
            "model": _kurtosis_corr(*model[c], w),
            "ref_wall_s_per_part": [runs[(c, p)].seconds for p in PARTS],
        }
    return {
        **base,
        "row_verdict": res["row_verdict"],
        "pairs": res["pairs"],
        "cases": {str(c): v for c, v in res["cases"].items()},
        "bootstrap": {
            "n_boot": res["n_boot"],
            "analysis_seed": res["analysis_seed"],
            "scheme": "TOPAS within part, ionmc over histories, SD of dev",
        },
        "seal": {
            "path": SEAL_PATH,
            "sha256": reader.reads[SEAL_PATH],
            "implementation_sha": sealed["implementation_sha"],
            "cases": [
                {
                    k: e[k]
                    for k in (
                        "case",
                        "run_id",
                        "npz_sha256",
                        "seed",
                        "theta_anchor_rad",
                        "window_rad",
                    )
                }
                for e in sealed["cases"]
            ],
        },
        "report_only": report_only,
        "evaluator": {
            "head_sha": reader.head,
            "committed_inputs": dict(sorted(reader.reads.items())),
        },
    }


# =============================================================================================
# CLI
# =============================================================================================

_FAIL_CLOSED = (RecordError, r2e.ChronologyError, OSError, ValueError, KeyError, TypeError)


def _json_safe(v: Any) -> Any:
    if isinstance(v, dict):
        return {str(k): _json_safe(x) for k, x in v.items() if k != "angles"}
    if isinstance(v, list | tuple):
        return [_json_safe(x) for x in v]
    if isinstance(v, np.ndarray):
        raise TypeError("raw arrays are never written to a V5-SS report")
    if isinstance(v, np.generic):
        v = v.item()
    if isinstance(v, float) and not math.isfinite(v):
        return None
    return v


def _write(doc: dict[str, Any], out: Path | None) -> None:
    doc = _json_safe({**doc, "generated_at": datetime.now(UTC).isoformat()})
    text = json.dumps(doc, indent=1, sort_keys=True, allow_nan=False) + "\n"
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="V3-008 V5-SS seal / check-run / grade")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("seal")
    s.add_argument("--implementation-sha", required=True)
    s.add_argument("--host-records", type=Path, default=Path(HOST_RECORDS_DIR))
    s.add_argument("--samples-root", type=Path, required=True)
    s.add_argument("--out", type=Path, default=None)
    c = sub.add_parser("check-run")
    c.add_argument("--run-id", required=True)
    c.add_argument("--bundle", required=True)
    c.add_argument("--records-dir", type=Path, default=None)
    c.add_argument("--outputs", type=Path, default=None)
    c.add_argument("--out", type=Path, default=None)
    g = sub.add_parser("grade")
    g.add_argument("--outputs", type=Path, default=None)
    g.add_argument("--samples-root", type=Path, default=None)
    g.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    try:
        if a.cmd == "seal":
            doc = seal(a.implementation_sha, a.host_records, a.samples_root)
            text = json.dumps(_json_safe(doc), indent=1, sort_keys=True, allow_nan=False) + "\n"
            out = a.out or r2e.evaluator_repo() / SEAL_PATH
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(text, encoding="utf-8")
            return 0
        if a.cmd == "check-run":
            doc = check_run(a.run_id, a.bundle, a.records_dir, a.outputs)
            code = 0 if doc["verdict"] == "VALID" else 1
        else:
            doc = grade(a.outputs, a.samples_root)
            code = 0 if doc["row_verdict"] == "PASS" else 1
    except _FAIL_CLOSED as exc:
        doc = {
            "cmd": a.cmd,
            "status": "ERROR",
            "error": f"{type(exc).__name__}: {exc}",
            "note": "fail closed: nothing established; not a pass",
        }
        print(f"ss_evaluate {a.cmd}: ERROR (fail closed): {exc}", file=sys.stderr)
        code = 2
    _write(doc, a.out)
    return code


if __name__ == "__main__":
    sys.exit(main())
