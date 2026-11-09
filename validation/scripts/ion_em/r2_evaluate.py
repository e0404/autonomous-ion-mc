"""V3-008 revision 2: evaluator of the TOPAS-agreement rows (TD-EM, TD-FULL, TL-EM, TL-FULL).

Entry point::

    python validation/scripts/ion_em/r2_evaluate.py grade --case <n>
        --row <TD-EM|TD-FULL|TL-EM|TL-FULL> [--records <committed dir>] [--seal <committed path>]
        [--outputs <dir>] [--samples-root <dir>] [--out <report.json>]
    python validation/scripts/ion_em/r2_evaluate.py check-run --run-id <REF-...>
        --bundle <bundle dir name> [--records-dir <committed dir>] [--outputs <dir>]
        [--out <report.json>]

``check-run`` establishes the run validity of ONE signed run (evidence part or smoke; nothing
graded, see :func:`check_run`): exit 0 VALID, 1 BLINDING FAIL / RUN-VALIDITY FAIL, 2 ERROR.

The repository is always the git checkout that contains this file; plan seal note, review record,
run index, signed run records, listings, seal and trust file are read as committed at its HEAD
(a differing working copy fails closed; Codex REVIEW-a50374b1).

Design record ``validation/plans/v3-008/r2-topas-target.md``: section 4 estimators and section 5
criteria are in ``r2_statistics`` (numpy only), the fail-closed loaders of the controlled records,
the sealed ionmc chunks and the TOPAS outputs in ``r2_io`` (read its docstring for the evidence
directory layout and the record bindings). Exit status:
0 PASS, 1 FAIL / NOT EVALUABLE, 2 ERROR (fail closed; a report with ``status: ERROR`` is still
written). The report holds statistics, record ids and sha256 only, never raw arrays.

Chronology (section 6, "Run validity against the seal"), unchanged:

1. Ancestry: ``git merge-base --is-ancestor <seal> <code_sha>`` succeeds (the seal itself counts).
2. Review before start: the controlled review record proves that the seal commit itself passed
   its independent review (``head_sha`` equals the plan's ``seal_sha``, ``status`` is ``passed``
   and the review ``report`` carries ``verdict == "pass"``; a missing report or verdict fails
   closed) and its ``finished_at`` is strictly earlier than the run's ``reference_started``
   (records of the controlled services, never Git timestamps; the run's own ``finished_at`` is
   not used).

Every TOPAS run the grade uses (including the IDD-EM runs that fix R80_T,EM for the lateral depths
and the low-end case of TL (b)) is checked first; one invalid run makes the row a FAIL of the
blinding control and nothing is graded ("never re-graded").

Record shapes for the chronology checks (dicts; ISO-8601 timestamps with an explicit offset or Z):
``plan`` ``{"seal_sha", "seal_review_id"}``; ``review`` ``{"review_id", "head_sha", "status",
"finished_at", "report": {"verdict"}}``; ``run`` ``{"run_id", "code_sha", "reference_started"}``
(built by ``r2_io.TopasRun.chronology_record`` from the signed run record bound by run_id).

Statuses: PASS; FAIL with ``category`` "blinding control" (chronology, never re-graded), "run
validity" (a TOPAS part fails ``r2_io.run_checks`` -- non-empty stderr, killed/stuck tracks,
incomplete histories, msc not UrbanMsc for GenericIon/alpha -- or a group fails the replicate gate
``r2_io.replicate_gate``; ``rerun_allowed``: the part is re-run at a fresh re-run base per section
6) or a graded criterion FAIL; NOT EVALUABLE (orchestrator decision 2026-11-09): not a pass and not
a model FAIL but a recorded evidence gap (TL (a) without a gated point, TL (b) without common gated
depths); exit 1, so it blocks an integration criterion that needs PASS; ERROR: a missing or
malformed record, nothing graded.

The run-validity check includes the executed module list (``r2_io.module_check``: printout markers
derived from the bundle's ``sv:Ph/Default/Modules``; interim rule: a module without runtime
evidence is a run-validity FAIL, ``g4decay`` is evidenced by a registered ``Decay`` process in the
``Ph/ListProcesses`` printout). TOPAS outputs are bound to the signed record through its
``output_files_sha256``. The chronology uses the signed ``started_at`` of
the signed reference-runner record.json only (see ``r2_io``). Not done here: the blinding ORDER of
grading (orchestrator); the report-only items listed under ``not_computed`` in every report.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r2_io as io  # noqa: E402
import r2_statistics as st  # noqa: E402

SHA_RE = re.compile(r"^[0-9a-f]{40}$")
REPORT_SCHEMA = "ionmc-v3-008-r2-grade-1"
ROWS = ("TD-EM", "TD-FULL", "TL-EM", "TL-FULL")
DEPTH_FRACTIONS = (0.25, 0.5, 0.75, 0.9)  # section 4.2
DEXT_FRACTIONS = (0.75, 0.9)  # section 5.3 (b)
A_FRAG = {"he": 0.0015, "carbon": 0.0010, "oxygen": 0.0010}  # section 5.2
TL_BAND = (-0.05, 0.40)  # section 5.3 (a)
DEXT_TOL = 0.12  # section 5.3 (b)
NOT_COMPUTED = [
    "TD-EM report-only: ratio curve E_i/E_T over [0, R80_T + 5 sigma_R]",
    "TD-EM report-only: C/O 430 equal-velocity R80 ratio",
    "TD-FULL report-only: tail dose beyond R80_T + 3 sigma_R",
    "TD-FULL report-only: primary fluence vs depth (and TF-DIAG primary survival at R80)",
    "TL report-only: D_std; C/O 430 equal-velocity sigma ratio",
    "TF-DIAG: FULL/EM plateau and peak ratios",
]
DECISIONS = [
    "D1 pooled point values (estimator on pooled batches / summed TOPAS parts)",
    "D2 batch-means / part substitution uncertainties; sigma_T = 1.430 SD/sqrt(4), K = 4 only",
    "D3 relative tolerances graded as |x - r| + 1.645 sigma_comb <= frac * r_TOPAS",
    "D4 whole-bin windows by bin centre; fit window = bins with |centre| <= 2 sigma",
    "D5 range classes of section 4.1 from R_plan (section 3)",
    "D6 chunks combined per batch: history-weighted mean of the chunks' batch b (B = 40)",
    "D7 r and <r> of TL use the pooled TOPAS sigma per point; D_ext sigma from the same "
    "substitution (depth correlations kept)",
    "D8 TL (a) with no gated point and TL (b) without common gated depths: NOT EVALUABLE",
]


class ChronologyError(Exception):
    """A chronology record is missing or malformed: the check fails closed (never passes)."""


def _require(record: dict[str, Any] | None, key: str, what: str) -> Any:
    if not isinstance(record, dict) or record.get(key) in (None, ""):
        raise ChronologyError(f"missing record: {what} ({key})")
    return record[key]


def _sha(value: Any, what: str) -> str:
    if not isinstance(value, str) or not SHA_RE.match(value):
        raise ChronologyError(f"{what} is not a full 40-hex commit SHA: {value!r}")
    return value


def _time(value: Any, what: str) -> datetime:
    if not isinstance(value, str):
        raise ChronologyError(f"{what} is not a timestamp string: {value!r}")
    try:
        t = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ChronologyError(f"{what} is not an ISO-8601 timestamp: {value!r}") from exc
    if t.tzinfo is None:
        raise ChronologyError(f"{what} has no UTC offset: {value!r}")
    return t.astimezone(UTC)


def check_ancestry(repo: Path | str, seal_sha: Any, code_sha: Any) -> bool:
    """True iff ``seal_sha`` is an ancestor of (or equal to) ``code_sha`` in ``repo``; raises
    :class:`ChronologyError` for a missing/malformed SHA or a git failure other than "not
    ancestor"."""
    seal = _sha(seal_sha, "seal SHA")
    code = _sha(code_sha, "run code_sha")
    proc = subprocess.run(  # noqa: S603 - fixed argv, validated hex SHAs
        ["git", "-C", str(repo), "merge-base", "--is-ancestor", seal, code],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode == 0:
        return True
    if proc.returncode == 1:
        return False
    raise ChronologyError(f"git merge-base failed (rc={proc.returncode}): {proc.stderr.strip()}")


def check_review_before_start(
    review: dict[str, Any] | None, plan: dict[str, Any] | None, run: dict[str, Any] | None
) -> bool:
    """True iff the seal review recorded in ``plan`` passed for the seal commit and finished
    strictly before the run started. Any missing, malformed, mismatched or non-passing field raises
    :class:`ChronologyError`."""
    review_id = _require(plan, "seal_review_id", "plan seal note")
    seal = _sha(_require(plan, "seal_sha", "plan seal note"), "plan seal_sha")
    if _require(review, "review_id", "seal review record") != review_id:
        raise ChronologyError("the review record is not the review named by the plan's seal note")
    head = _sha(_require(review, "head_sha", "seal review head_sha"), "review head_sha")
    if head != seal:
        raise ChronologyError(f"the review covers commit {head}, not the plan's seal commit {seal}")
    status = _require(review, "status", "seal review status")
    if status != "passed":
        raise ChronologyError(f"the seal review did not pass (status={status!r})")
    report = _require(review, "report", "seal review report")
    if not isinstance(report, dict):
        raise ChronologyError(f"the review report is malformed: {report!r}")
    verdict = _require(report, "verdict", "seal review report verdict")
    if verdict != "pass":
        raise ChronologyError(f"the seal review verdict is not pass: {verdict!r}")
    finished = _time(_require(review, "finished_at", "seal review finished_at"), "finished_at")
    started = _time(
        _require(run, "reference_started", "run reference_started"), "reference_started"
    )
    return finished < started


def check_run_validity(
    repo: Path | str,
    plan: dict[str, Any] | None,
    review: dict[str, Any] | None,
    run: dict[str, Any] | None,
    expected_run_id: str,
) -> dict[str, Any]:
    """Both chronology checks of one TOPAS run. ``valid`` is True only if both hold; any missing
    record raises :class:`ChronologyError` (fail closed). ``run`` must be the run's own controlled
    record (see the module docstring): it must carry ``run_id``, ``code_sha`` and
    ``reference_started``, and its ``run_id`` must equal ``expected_run_id`` (mandatory)."""
    if not isinstance(expected_run_id, str) or not expected_run_id:
        raise ChronologyError(f"expected_run_id is required: {expected_run_id!r}")
    run_id = _require(run, "run_id", "run record run_id")
    code_sha = _require(run, "code_sha", "run record code_sha")
    _require(run, "reference_started", "run record reference_started")
    if run_id != expected_run_id:
        raise ChronologyError(
            f"the run record ({run_id!r}) is not bound to the run being evaluated "
            f"({expected_run_id!r})"
        )
    seal = _require(plan, "seal_sha", "plan seal note")
    ancestry = check_ancestry(repo, seal, code_sha)
    reviewed = check_review_before_start(review, plan, run)
    return {"ancestry": ancestry, "review_before_start": reviewed, "valid": ancestry and reviewed}


# =============================================================================================
# TD statistics (sections 4.1, 5.1, 5.2)
# =============================================================================================


def _safe(fn: Callable[[], Any]) -> tuple[Any, str | None]:
    try:
        return fn(), None
    except (st.StatisticsError, ValueError, ZeroDivisionError, FloatingPointError) as exc:
        return None, str(exc) or type(exc).__name__


def _quantities(e: np.ndarray, hz: float, r_plan: float) -> dict[str, Any]:
    """Per-quantity values; a failing quantity is None with its reason under ``err:<key>``."""
    out: dict[str, Any] = {}
    pk, err = _safe(lambda: st.peak(e, hz))
    if pk is None:
        for k in ("z_peak", "d_peak", "r80", "r20", "r90", "w", "pe"):
            out[k], out["err:" + k] = None, err
    else:
        out["z_peak"], out["d_peak"] = pk
        for key, lev in (("r80", 0.8), ("r20", 0.2), ("r90", 0.9)):
            out[key], out["err:" + key] = _safe(
                lambda lev=lev: st.distal_level(e, hz, pk[0], lev * pk[1])
            )
        out["w"] = (
            out["r20"] - out["r80"] if out["r20"] is not None and out["r80"] is not None else None
        )
    out["d_ent"], out["err:d_ent"] = _safe(lambda: st.entrance_dose(e, hz, r_plan))
    if out.get("d_peak") is not None and out["d_ent"]:
        out["pe"] = out["d_peak"] / out["d_ent"]
    else:
        out["pe"] = None
    return out


def _criterion(
    name: str,
    key: str,
    qi: dict[str, Any],
    qib: list[dict[str, Any]],
    qt: dict[str, Any],
    qtp: list[dict[str, Any]],
    tol_of: Callable[[float], float],
    gated: bool = True,
) -> dict[str, Any]:
    vals = [qi.get(key), qt.get(key)] + [q.get(key) for q in qib] + [q.get(key) for q in qtp]
    if any(v is None for v in vals):
        reasons = sorted(
            {q.get("err:" + key) or "undefined" for q in [qi, qt, *qib, *qtp] if q.get(key) is None}
        )
        res = st.tost(math.nan, math.nan, math.nan, math.nan, math.nan)
        res["note"] = "invalid: " + "; ".join(reasons)[:300]
    else:
        res = st.tost(
            qi[key],
            qt[key],
            st.batch_se([q[key] for q in qib]),
            st.topas_se([q[key] for q in qtp]),
            tol_of(qt[key]),
        )
    res["criterion"], res["observable"], res["gated"] = name, key, gated
    if not gated:
        res["verdict"] = "REPORT-ONLY"
    return res


def _r80_tol(r80_t: float) -> float:
    """TD-EM (a): max(min(0.35 % R80_T, 1.0 mm), 0.05 mm)."""
    return max(min(0.0035 * r80_t, 1.0), 0.05)


def td_statistics(
    row: str,
    case: int,
    ion_e: np.ndarray,
    topas_parts: list[tuple[np.ndarray, int]],
    hz: float,
    topas_em_parts: list[tuple[np.ndarray, int]] | None = None,
    topas_primary_parts: list[tuple[np.ndarray, int]] | None = None,
) -> dict[str, Any]:
    """TD-EM (section 5.1, every criterion (a)-(e) gating, (b) only where W_T >= 1.0 mm) or TD-FULL
    (section 5.2, gating on (a) R80 of the total IDD with tolerance TD-EM (a) + A_frag, and (c)
    z_peak 0.5 mm; PE, plateau, W and the entrance-dose ratio report-only).

    ``ion_e``: ionmc combined batches, energy per primary per depth bin [MeV], shape (B, nz).
    ``topas_parts``: K = 4 (run-total energy per depth bin [MeV], histories). Same edges.
    """
    info = io.case_info(case)
    r_plan = info["r_plan_mm"]

    def curve_t(parts: list[tuple[np.ndarray, int]]) -> tuple[np.ndarray, list[np.ndarray]]:
        pooled = sum(p[0] for p in parts) / float(sum(p[1] for p in parts)) / hz
        return pooled, [p[0] / float(p[1]) / hz for p in parts]

    ei_b = [np.asarray(b) / hz for b in ion_e]
    ei = np.mean(np.asarray(ion_e), axis=0) / hz
    et, et_p = curve_t(topas_parts)
    qi, qt = _quantities(ei, hz, r_plan), _quantities(et, hz, r_plan)
    qib = [_quantities(b, hz, r_plan) for b in ei_b]
    qtp = [_quantities(p, hz, r_plan) for p in et_p]
    crit: list[dict[str, Any]] = []
    report: dict[str, Any] = {}
    full = row == "TD-FULL"
    afrag = A_FRAG[info["ion"]]
    tol_a = (lambda r: _r80_tol(r) + afrag * r) if full else _r80_tol
    crit.append(_criterion("(a)", "r80", qi, qib, qt, qtp, tol_a))
    w_t = qt.get("w")
    w_gated = (not full) and (w_t is None or w_t >= 1.0)  # W_T undefined: fail closed (gated)
    crit.append(_criterion("(b)", "w", qi, qib, qt, qtp, lambda r: 0.10 * r, gated=w_gated))
    crit.append(_criterion("(c)", "z_peak", qi, qib, qt, qtp, lambda r: 0.5))
    crit.append(_criterion("(d)", "pe", qi, qib, qt, qtp, lambda r: 0.03 * r, gated=not full))
    r80_t = qt.get("r80")
    if r80_t is None:
        res = st.tost(math.nan, math.nan, math.nan, math.nan, math.nan)
        res.update(
            criterion="(e)", observable="plateau", gated=not full, note="invalid: R80_T undefined"
        )
        crit.append(res)
    else:
        groups = [
            _safe(lambda e=e: st.plateau_groups(e, hz, r_plan, r80_t))
            for e in [ei, et, *ei_b, *et_p]
        ]
        if any(g[0] is None for g in groups):
            res = st.tost(math.nan, math.nan, math.nan, math.nan, math.nan)
            res.update(
                criterion="(e)",
                observable="plateau",
                gated=not full,
                note="invalid: " + str(groups[0][1] or groups[1][1]),
            )
            crit.append(res)
        else:
            gi, gt = groups[0][0], groups[1][0]
            gib = np.array([g[0] for g in groups[2 : 2 + len(ei_b)]])
            gtp = np.array([g[0] for g in groups[2 + len(ei_b) :]])
            per = [
                st.tost(
                    float(gi[j]),
                    float(gt[j]),
                    st.batch_se(gib[:, j]),
                    st.topas_se(gtp[:, j]),
                    0.02 * float(gt[j]),
                )
                for j in range(gi.size)
            ]
            ok = all(p["verdict"] == "PASS" for p in per)
            worst = max(per, key=lambda p: (p["statistic"] or math.inf) - p["tol"])
            crit.append(
                {
                    "criterion": "(e)",
                    "observable": "plateau",
                    "gated": not full,
                    "groups": len(per),
                    "group_mm": st.range_class(r_plan)["group_mm"],
                    "start_mm": st.range_class(r_plan)["ent_lo"],
                    "end_mm": 0.8 * r80_t,
                    "failing_groups": [j for j, p in enumerate(per) if p["verdict"] != "PASS"],
                    "worst_group": worst,
                    "per_group": per,
                    "verdict": "PASS" if ok else "FAIL",
                }
            )
            if full:
                crit[-1]["verdict"] = "REPORT-ONLY"
    if full:
        crit[-1]["verdict"] = "REPORT-ONLY"
        report["entrance_dose_ratio"] = (
            (qi["d_ent"] / qt["d_ent"]) if qi.get("d_ent") and qt.get("d_ent") else None
        )
        if topas_primary_parts:
            qprim = _quantities(curve_t(topas_primary_parts)[0], hz, r_plan)
            report["r80_topas_iddprimary_mm"] = qprim.get("r80")
        if topas_em_parts:
            em_p, em_parts = curve_t(topas_em_parts)
            r80_em = _quantities(em_p, hz, r_plan).get("r80")
            em_r80s = [_quantities(p, hz, r_plan).get("r80") for p in em_parts]
            full_r80s = [q.get("r80") for q in qtp]
            if r80_em is not None and r80_t is not None and None not in em_r80s + full_r80s:
                report["delta_frag_mm"] = {
                    "value": r80_t - r80_em,
                    "sigma": math.hypot(st.topas_se(em_r80s), st.topas_se(full_r80s)),
                    "a_frag_mm": afrag * r80_t,
                }
    else:
        rep = _criterion("R90", "r90", qi, qib, qt, qtp, lambda r: math.inf, gated=False)
        report["r90"] = rep
    for c in crit:
        if c.get("gated") is False:
            c["verdict"] = "REPORT-ONLY"
    gating = [c for c in crit if c["verdict"] != "REPORT-ONLY"]
    status = "PASS" if gating and all(c["verdict"] == "PASS" for c in gating) else "FAIL"
    pooled = {
        k: {"ionmc": qi.get(k), "topas": qt.get(k)}
        for k in ("z_peak", "d_peak", "r80", "r20", "w", "d_ent", "pe")
    }
    return {"status": status, "criteria": crit, "pooled": pooled, "report_only": report}


# =============================================================================================
# TL statistics (sections 4.2, 5.3)
# =============================================================================================


def _fit_xy(
    px: np.ndarray, py: np.ndarray, ex: np.ndarray, ey: np.ndarray
) -> tuple[float | None, str | None]:
    return _safe(lambda: st.core_xy(px, py, ex, ey))


def tl_points(
    case: int,
    ion_lat: dict[str, np.ndarray],
    edges: dict[str, np.ndarray],
    topas_parts: list[dict[str, np.ndarray]],
    gate_sigma: list[float | None],
    depths: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Per (case, depth) point of TL (a). ``ion_lat['x'|'y']``: (B, n, D) combined batches;
    ``topas_parts``: K = 4 dicts of (n, D) run totals; ``gate_sigma``: pooled TOPAS-EM sigma_core
    per depth (gating "iff sigma_core,TOPAS-EM >= 3 h_x"; None = undeterminable = gated and
    invalid)."""
    hx = io.case_info(case)["hx_mm"]
    ex, ey = edges["x"], edges["y"]
    pts = []
    for k, d in enumerate(depths):
        ix, iy = ion_lat["x"][:, :, k], ion_lat["y"][:, :, k]
        si, ei = _fit_xy(ix.mean(axis=0), iy.mean(axis=0), ex, ey)
        si_b = [_fit_xy(ix[b], iy[b], ex, ey)[0] for b in range(ix.shape[0])]
        tx = [p["x"][:, k] for p in topas_parts]
        ty = [p["y"][:, k] for p in topas_parts]
        stt, et = _fit_xy(sum(tx), sum(ty), ex, ey)
        st_p = [_fit_xy(a, b, ex, ey)[0] for a, b in zip(tx, ty, strict=True)]
        g = gate_sigma[k]
        gated = g is None or g >= 3.0 * hx
        pt: dict[str, Any] = {
            **d,
            "gated": gated,
            "gate_sigma_topas_em": g,
            "sigma_ionmc": si,
            "sigma_topas": stt,
        }
        valid = (
            si is not None
            and stt is not None
            and None not in si_b
            and None not in st_p
            and g is not None
        )
        if valid:
            r = si / stt
            s_i = st.batch_se([b / stt for b in si_b])
            s_t = st.topas_se([si / p for p in st_p])
            res = st.band(r, math.hypot(s_i, s_t), *TL_BAND)
            res.update(sigma_i=s_i, sigma_t=s_t)
            pt["_batches"], pt["_parts"] = si_b, st_p
        else:
            why = ei or et or ("gate sigma undefined" if g is None else "a batch/part fit failed")
            res = st.band(math.nan, math.nan, *TL_BAND)
            res["note"] = f"invalid point: {why}"
        if not gated:
            res["verdict"] = "REPORT-ONLY"
        pt["a"] = res
        for name, prof, e in (
            ("ionmc", (ix.mean(axis=0), iy.mean(axis=0)), (ex, ey)),
            ("topas", (sum(tx), sum(ty)), (ex, ey)),
        ):
            mom, _ = _safe(
                lambda prof=prof, e=e: 0.5
                * (st.second_moment_sigma(prof[0], e[0]) + st.second_moment_sigma(prof[1], e[1]))
            )
            pt[f"moment20_{name}"] = mom
            sc = si if name == "ionmc" else stt
            pt[f"halo_{name}"] = (
                _safe(
                    lambda prof=prof, e=e, sc=sc: 0.5
                    * (st.halo_fraction(prof[0], e[0], sc) + st.halo_fraction(prof[1], e[1], sc))
                )[0]
                if sc
                else None
            )
        pts.append(pt)
    return pts


def _mean_r(pts: list[dict[str, Any]]) -> tuple[float, float, float]:
    """<r> over points with (sigma_i, sigma_T) from the per-batch / per-part substitution (D7)."""
    r = float(np.mean([p["sigma_ionmc"] / p["sigma_topas"] for p in pts]))
    nb = len(pts[0]["_batches"])
    rb = [float(np.mean([p["_batches"][b] / p["sigma_topas"] for p in pts])) for b in range(nb)]
    rp = [
        float(np.mean([p["sigma_ionmc"] / p["_parts"][j] for p in pts])) for j in range(st.K_TOPAS)
    ]
    return r, st.batch_se(rb), st.topas_se(rp)


def d_ext(low: list[dict[str, Any]], up: list[dict[str, Any]]) -> dict[str, Any]:
    """Section 5.3 (b): D_ext = <r>_upper/<r>_low - 1 over the gated depths in {0.75, 0.9} common to
    both cases; PASS iff |D_ext| + 1.645 sigma <= 12 %; NOT EVALUABLE (recorded, not a pass) when
    the low-end case has no gated point there (or none is common)."""
    lo = {p["f"]: p for p in low if p["f"] in DEXT_FRACTIONS and p["gated"]}
    hi = {p["f"]: p for p in up if p["f"] in DEXT_FRACTIONS and p["gated"]}
    common = sorted(set(lo) & set(hi))
    if not lo or not common:
        return {
            "verdict": "NOT EVALUABLE",
            "depths": common,
            "note": "no common gated depth in {0.75, 0.9}",
        }
    sel_l, sel_u = [lo[f] for f in common], [hi[f] for f in common]
    if any("_batches" not in p for p in sel_l + sel_u):
        out = st.tost(math.nan, 0.0, math.nan, math.nan, DEXT_TOL)
        return out | {"depths": common, "note": "invalid: a point of D_ext is invalid"}
    rl, sl_i, sl_t = _mean_r(sel_l)
    ru, su_i, su_t = _mean_r(sel_u)
    d = ru / rl - 1.0
    s_i = (1.0 + d) * math.hypot(su_i / ru, sl_i / rl)
    s_t = (1.0 + d) * math.hypot(su_t / ru, sl_t / rl)
    return st.tost(d, 0.0, s_i, s_t, DEXT_TOL) | {"depths": common, "r_low": rl, "r_upper": ru}


def tl_statistics(points: list[dict[str, Any]], dext: dict[str, Any] | None) -> dict[str, Any]:
    """Row status of TL-EM / TL-FULL from its points (criterion (a)) and, for an upper-end case,
    criterion (b). FAIL if any gating verdict fails; NOT EVALUABLE (not a pass) if (a) has no
    gated point or (b) is not evaluable; else PASS."""
    gated = [p["a"] for p in points if p["gated"]]
    verdicts = [g["verdict"] for g in gated]
    a = "NOT EVALUABLE" if not gated else ("PASS" if all(v == "PASS" for v in verdicts) else "FAIL")
    b = dext["verdict"] if dext is not None else None
    allv = [a] + ([b] if b else [])
    status = "FAIL" if "FAIL" in allv else ("NOT EVALUABLE" if "NOT EVALUABLE" in allv else "PASS")
    clean = [{k: v for k, v in p.items() if not k.startswith("_")} for p in points]
    return {"status": status, "criterion_a": a, "points": clean, "criterion_b": dext}


# =============================================================================================
# grading glue and CLI
# =============================================================================================


class Context:
    def __init__(
        self,
        repo: Path,
        plan: dict[str, Any],
        review: dict[str, Any],
        index: list[dict[str, Any]],
        records: Path,
        outputs: Path,
        seal: dict[str, Any],
        samples_root: Path,
        trust: dict[str, Any],
    ) -> None:
        self.trust = trust
        self.reader: io.CommittedReader | None = None
        self.repo, self.plan, self.review, self.index = repo, plan, review, index
        self.records, self.outputs, self.seal, self.samples_root = (
            records,
            outputs,
            seal,
            samples_root,
        )
        self.runs: list[io.TopasRun] = []
        self.chunks_used: list[dict[str, Any]] = []

    def topas(self, case: int, config: str) -> list[io.TopasRun]:
        runs = [
            io.load_topas_run(
                self.repo,
                self.records,
                self.outputs,
                case,
                config,
                p,
                rid,
                self.trust,
                self.reader,
            )
            for p, rid in io.runs_for(self.index, case, config)
        ]
        self.runs.extend(runs)
        return runs

    def chunks(self, case: int) -> list[dict[str, Any]]:
        ch = io.case_chunks(self.seal, case)
        for c in ch:
            entry = {
                "case": case,
                "k": c["k"],
                "run_id": c["run_id"],
                "npz_sha256": c["npz_sha256"],
            }
            if entry not in self.chunks_used:
                self.chunks_used.append(entry)
        return ch


def run_chronology(
    repo: Path,
    reader: io.CommittedReader | None,
    run: io.TopasRun,
    plan: dict[str, Any],
    review_records: Path,
) -> dict[str, Any]:
    """R2-TOPAS-4 chronology of one loaded run: ``check_run_validity`` with the seal note and
    review record pinned at the run's signed ``code_sha`` (which must equal HEAD's seal note)."""
    pinned = run.pinned  # seal note and review record at the run's signed code_sha
    if pinned["plan"] != plan:
        raise io.RecordError(
            f"{run.run_id}: the seal note at code_sha {run.code_sha} differs from HEAD's"
        )
    v = check_run_validity(
        repo, pinned["plan"], pinned["review"], run.chronology_record(), run.run_id
    )
    return {
        "run_id": run.run_id, "config": run.config, "case": run.case, "part": run.part,
        "code_sha": run.code_sha, "review_finished_at": pinned["review"]["finished_at"],
        "review_sha256_at_code_sha": pinned["review_sha256"],
        "review_sha256_at_head": reader.reads.get(
            f"{reader.rel(review_records)}/review.json"
        ) if reader else None,
        "plan_sha256_at_code_sha": pinned["plan_sha256"],
        "manifest_bundle_set_digest": pinned["manifest_bundle_set_digest"], **v,
    }  # fmt: skip


def _validity(ctx: Context) -> list[dict[str, Any]]:
    return [run_chronology(ctx.repo, ctx.reader, r, ctx.plan, ctx.records) for r in ctx.runs]


def _r80_em_topas(ctx: Context, case: int) -> tuple[float, list[io.TopasRun]]:
    runs = ctx.topas(case, "idd-em")
    parts = [io.topas_idd(r) for r in runs]
    hz = parts[0][1]
    pooled = sum(p[0] for p in parts) / float(sum(r.histories for r in runs)) / hz
    return st.depth_quantities(pooled, hz, io.case_info(case)["r_plan_mm"])["r80"], runs


def _lateral(
    ctx: Context, case: int, row: str, r80_em: float
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    info = io.case_info(case)
    length = info["length_mm"]
    slab_edges = np.arange(round(length) + 1, dtype=np.float64)
    depths = []
    for f in DEPTH_FRACTIONS:
        z = f * r80_em
        k = int(np.searchsorted(slab_edges, z, side="right") - 1)
        if not 0 <= k < slab_edges.size - 1:
            raise io.RecordError(f"depth {z} mm outside the lateral slabs")
        depths.append(
            {
                "f": f,
                "z_mm": z,
                "slab": k,
                "slab_mm": [float(slab_edges[k]), float(slab_edges[k + 1])],
            }
        )
    slabs = [d["slab"] for d in depths]

    def topas_profiles(
        runs: list[io.TopasRun], sx: str, sy: str
    ) -> tuple[list[dict[str, np.ndarray]], dict[str, np.ndarray]]:
        parts, edges = [], {}
        for r in runs:
            vx, ex, wz = io.topas_lateral(r, sx, 0)
            vy, ey, wz2 = io.topas_lateral(r, sy, 1)
            if abs(wz - 1.0) > 1e-9 or abs(wz2 - 1.0) > 1e-9 or vx.shape[1] != slab_edges.size - 1:
                raise io.RecordError(f"{r.run_id}: lateral depth bins are not 1 mm over L")
            parts.append({"x": vx[:, slabs], "y": vy[:, slabs]})
            edges = {"x": ex, "y": ey}
        return parts, edges

    em_parts, edges = topas_profiles(ctx.topas(case, "lat-em"), "lat_x", "lat_y")
    gate = [
        _fit_xy(
            sum(p["x"][:, k] for p in em_parts),
            sum(p["y"][:, k] for p in em_parts),
            edges["x"],
            edges["y"],
        )[0]
        for k in range(len(slabs))
    ]
    extra: dict[str, Any] = {}
    if row == "TL-FULL":
        full_runs = ctx.topas(case, "lat-full")
        parts, _ = topas_profiles(full_runs, "lat_x_primary", "lat_y_primary")
        allp, _ = topas_profiles(full_runs, "lat_x", "lat_y")
        sig_all = [
            _fit_xy(
                sum(p["x"][:, k] for p in allp),
                sum(p["y"][:, k] for p in allp),
                edges["x"],
                edges["y"],
            )[0]
            for k in range(len(slabs))
        ]
        extra["sigma_topas_full_all_particle"] = sig_all
    else:
        parts = em_parts
    chunks = ctx.chunks(case)
    gx = io.load_ionmc(ctx.samples_root, case, chunks, "latx", slabs)
    gy = io.load_ionmc(ctx.samples_root, case, chunks, "laty", slabs)
    io.same_edges(gx.edges[0], edges["x"], "lateral x")
    io.same_edges(gy.edges[1], edges["y"], "lateral y")
    io.same_edges(gx.edges[2], slab_edges, "lateral depth slabs")
    io.same_edges(gy.edges[2], slab_edges, "lateral depth slabs")
    pts = tl_points(case, {"x": gx.batches, "y": gy.batches}, edges, parts, gate, depths)
    if row == "TL-FULL":
        extra["sigma_full_primary_over_em"] = [
            (p["sigma_topas"] / g) if p["sigma_topas"] and g else None
            for p, g in zip(pts, gate, strict=True)
        ]
    return pts, extra


def evaluator_repo() -> Path:
    """The git toplevel of the checkout that contains this evaluator (the only repository the
    production CLI grades from; Codex REVIEW-a50374b1)."""
    here = Path(__file__).resolve()
    proc = subprocess.run(  # noqa: S603 - fixed argv
        ["git", "-C", str(here.parent), "rev-parse", "--show-toplevel"],
        capture_output=True, text=True, check=False,
    )  # fmt: skip
    if proc.returncode != 0:
        raise io.RecordError(f"the evaluator is not inside a git checkout: {proc.stderr.strip()}")
    top = Path(proc.stdout.strip()).resolve()
    if not here.is_relative_to(top):
        raise io.RecordError("the evaluator does not live inside its git toplevel")
    return top


def grade(
    records: Path | str,
    seal_path: Path | str,
    case: int,
    row: str,
    outputs: Path | None = None,
    samples_root: Path | None = None,
    *,
    repo: Path | None = None,
) -> dict[str, Any]:
    """Grade one (case, row). Raises on any missing/malformed record (caller: ERROR).

    Every repository input (acceptance-plan seal note, review record, run index, signed run
    records, artifact listings, sample seal, trust file, bundles) is read as committed at the
    evaluator repository's HEAD (bundles at the run's code_sha) and a differing working copy fails
    closed. ``repo`` is a test-only override (the CLI never exposes it): production grades use
    :func:`evaluator_repo`."""
    if row not in ROWS:
        raise io.RecordError(f"unknown row {row!r}; rows are {ROWS}")
    repo = Path(repo).resolve() if repo is not None else evaluator_repo()
    reader = io.CommittedReader(repo)
    info = io.case_info(case)
    plan = io.load_plan(reader)
    review = io.load_review(reader, records, plan)
    seal, seal_digest = io.load_seal(reader, seal_path, plan["seal_sha"])
    trust, trust_info = io.committed_trust(reader)
    ctx = Context(
        repo,
        plan,
        review,
        io.load_index(reader, records),
        Path(records) if Path(records).is_absolute() else repo / records,
        outputs or repo / ".ionmc-cache" / "reference-runs",
        seal,
        samples_root or repo,
        trust,
    )
    ctx.reader = reader
    report: dict[str, Any] = {
        "schema": REPORT_SCHEMA,
        "row": row,
        "case": case,
        "case_info": info,
        "trust_file": trust_info,
        "evaluator": {
            "head_sha": reader.head,
            "committed_inputs_clean": True,
            "committed_inputs": reader.reads,  # filled while grading (path -> sha256)
        },
        "held_out": info["set"] == "held-out",
        "plan": plan,
        "review": {k: review[k] for k in ("review_id", "head_sha", "status", "finished_at")},
        "seal_record": {"path": io.SEAL_RECORD_PATH, "sha256": seal_digest},
        "implementation_decisions": DECISIONS,
        "not_computed": NOT_COMPUTED,
    }
    # load every TOPAS record the grade needs, then check the chronology BEFORE reading numbers
    cases = [case]
    lowc, _std, upc = io.ION_ENDS[info["ion"]]
    tl_b = row.startswith("TL") and case == upc
    if tl_b:
        cases.append(lowc)
    configs = {
        "TD-EM": ["idd-em"],
        "TD-FULL": ["idd-full"],
        "TL-EM": ["idd-em", "lat-em"],
        "TL-FULL": ["idd-em", "lat-em", "lat-full"],
    }[row]
    for c in cases:
        for cfg in configs:
            ctx.topas(c, cfg)
    if row == "TD-FULL" and any(e["case"] == case and e["config"] == "idd-em" for e in ctx.index):
        ctx.topas(case, "idd-em")
    validity = _validity(ctx)
    report["run_validity"] = validity
    if not all(v["valid"] for v in validity):
        bad = [v["run_id"] for v in validity if not v["valid"]]
        return report | {
            "status": "FAIL",
            "category": "blinding control",
            "rerun_allowed": False,
            "reason": f"blinding control: run(s) {bad} fail check_run_validity; "
            "the TD/TL rows of this case are a FAIL and are never re-graded",
            "topas_runs": [r.summary() for r in ctx.runs],
            "ionmc_chunks": [],
        }
    # per-part run validity and the replicate gate (section 6 "Run validity (each TOPAS part)")
    part_checks = [
        {"run_id": r.run_id, "case": r.case, "config": r.config, "part": r.part,
         **io.run_checks(r)}
        for r in ctx.runs
    ]  # fmt: skip
    report["module_markers_unverifiable"] = sorted(
        {m for x in part_checks for m in x["module_markers_unverifiable"]}
    )
    groups: dict[tuple[int, str], list[io.TopasRun]] = {}
    for r in ctx.runs:
        groups.setdefault((r.case, r.config), []).append(r)
    gates = [
        {"case": c, "config": cfg, "problems": io.replicate_gate(rs)}
        for (c, cfg), rs in sorted(groups.items())
    ]
    report["part_validity"] = part_checks
    report["replicate_gate"] = gates
    if any(x["problems"] for x in part_checks + gates):
        return report | {
            "status": "FAIL",
            "category": "run validity",
            "rerun_allowed": True,
            "reason": "run validity (section 6, each TOPAS part): an invalid part is re-run at a "
            "fresh re-run base declared by amendment, the invalid part preserved; not graded",
            "topas_runs": [r.summary() for r in ctx.runs],
            "ionmc_chunks": [],
        }
    ctx.runs = []  # reload while reading outputs, so the report lists the files actually read
    if row.startswith("TD"):
        cfg = "idd-em" if row == "TD-EM" else "idd-full"
        runs = ctx.topas(case, cfg)
        parts = [io.topas_idd(r) for r in runs]
        hz = parts[0][1]
        if any(abs(p[1] - hz) > 1e-12 for p in parts) or abs(hz - info["hz_mm"]) > 1e-9:
            raise io.RecordError(
                f"IDD bin width {hz} mm differs from the case h_z {info['hz_mm']} mm"
            )
        grid = io.load_ionmc(ctx.samples_root, case, ctx.chunks(case), "idd")
        io.same_edges(grid.edges[2], hz * np.arange(parts[0][0].size + 1), "IDD depth")
        io.same_edges(grid.edges[0], [-io.FIELD_MM / 2, io.FIELD_MM / 2], "IDD x extent")
        io.same_edges(grid.edges[1], [-io.FIELD_MM / 2, io.FIELD_MM / 2], "IDD y extent")
        em_parts = prim = None
        if row == "TD-FULL":
            prim = [
                (io.topas_idd(r, "idd_primary_dose.csv")[0], r.histories)
                for r in runs
                if "idd_primary_dose.csv" in r.expected_outputs
            ] or None
            if any(e["case"] == case and e["config"] == "idd-em" for e in ctx.index):
                em_parts = [(io.topas_idd(r)[0], r.histories) for r in ctx.topas(case, "idd-em")]
        res = td_statistics(
            row,
            case,
            grid.batches,
            [(p[0], r.histories) for p, r in zip(parts, runs, strict=True)],
            hz,
            em_parts,
            prim,
        )
    else:
        r80_em, _ = _r80_em_topas(ctx, case)
        pts, extra = _lateral(ctx, case, row, r80_em)
        dext = None
        if tl_b:
            r80_low, _ = _r80_em_topas(ctx, lowc)
            low_pts, _ = _lateral(ctx, lowc, row, r80_low)
            dext = d_ext(low_pts, pts)
            dext["low_case"], dext["upper_case"] = lowc, case
        res = tl_statistics(pts, dext)
        res["r80_topas_em_mm"] = r80_em
        res["report_only"] = extra
        if not tl_b:
            res["criterion_b_note"] = (
                "graded with the upper-end case"
                if case in (lowc, _std)
                else "no D_ext for the held-out set (section 5.4)"
            )
    return (
        report
        | res
        | {
            "topas_runs": [r.summary() for r in ctx.runs],
            "ionmc_chunks": ctx.chunks_used,
        }
    )


CHECK_SCHEMA = "ionmc-v3-008-r2-check-run-1"
SMOKE_RECORDS = "validation/results/transport/v8r2-topas/smoke-runs"
VERDICTS = ("VALID", "RUN-VALIDITY FAIL", "BLINDING FAIL", "ERROR")


def check_run(
    run_id: str,
    bundle: str,
    records_dir: Path | str | None = None,
    outputs: Path | str | None = None,
    *,
    repo: Path | None = None,
) -> dict[str, Any]:
    """Run validity of ONE signed TOPAS run (evidence part or smoke), nothing graded.

    The same steps as :func:`grade` for each run it uses, through the same loaders:
    ``io.load_signed_run`` (signature against the committed trust file at HEAD, schema 3, run id
    derived from the configuration id, request binding incl. ``request_sha256`` and the bundle
    committed at ``code_sha``, frozen-manifest binding of ``case.json`` / ``input.txt``,
    stdout/stderr sha256 against the listing and the signed record, ``output_files_sha256``
    reproduced from the listing's ``work/`` entries), :func:`run_chronology` (seal note and
    ``<DEFAULT_RECORDS>/review.json`` at the run's ``code_sha``; ``code_sha`` an ancestor of HEAD)
    and ``io.run_checks`` for the configuration of the bundle name (``io.parse_bundle``). The
    scorer outputs (``work/``) are NOT read: only ``request.json``, ``stdout.txt`` and
    ``stderr.txt`` are read from ``<outputs>/<run_id>/`` (each hashed against the listing); the
    outputs are bound through the signed ``output_files_sha256`` recomputed from the listing's
    hashes. The replicate gate needs a group of parts and is not applicable to a single run.
    The seed is checked against the section-7 formula for an evidence part; a smoke's seed is
    bound through its committed bundle only.

    ``records_dir`` holds ``<run_id>/record.json`` and ``<run_id>/artifacts.json`` (default
    ``<DEFAULT_RECORDS>/runs`` for an evidence part, :data:`SMOKE_RECORDS` for a smoke); like every
    repository input it is read as committed at HEAD. Verdict: VALID, BLINDING FAIL (chronology),
    RUN-VALIDITY FAIL (``run_checks`` problems), ERROR (any record gap; raised to the caller)."""
    repo = Path(repo).resolve() if repo is not None else evaluator_repo()
    reader = io.CommittedReader(repo)
    info = io.parse_bundle(bundle)
    if records_dir is None:
        records_dir = (
            Path(io.DEFAULT_RECORDS) / "runs" if info["kind"] == "part" else Path(SMOKE_RECORDS)
        )
    records_dir = Path(records_dir) if Path(records_dir).is_absolute() else repo / records_dir
    review_records = repo / io.DEFAULT_RECORDS
    plan = io.load_plan(reader)
    review = io.load_review(reader, review_records, plan)
    trust, trust_info = io.committed_trust(reader)
    seeds = None
    if info["kind"] == "part":
        seeds = [io.topas_seed(info["case"], info["part"])]
    run = io.load_signed_run(
        repo, records_dir / run_id,
        Path(outputs) if outputs is not None else repo / ".ionmc-cache" / "reference-runs",
        run_id, bundle, info["config"], trust, reader, review_records=review_records,
        case=info["case"], part=info["part"], expected_seeds=seeds,
    )  # fmt: skip
    chron = run_chronology(repo, reader, run, plan, review_records)
    checks = io.run_checks(run)
    if not chron["valid"]:
        verdict = "BLINDING FAIL"
    elif checks["problems"]:
        verdict = "RUN-VALIDITY FAIL"
    else:
        verdict = "VALID"
    return {
        "schema": CHECK_SCHEMA,
        "run_id": run_id,
        "bundle": bundle,
        "bundle_kind": info["kind"],
        "config": info["config"],
        "case": info["case"],
        "part": info["part"],
        "verdict": verdict,
        "checks": {
            "record_binding": "PASS",
            "chronology": "PASS" if chron["valid"] else "FAIL",
            "chronology_ancestry": chron["ancestry"],
            "chronology_review_before_start": chron["review_before_start"],
            "run_checks": "FAIL" if checks["problems"] else "PASS",
            "run_checks_problems": checks["problems"],
            "module_markers_unverifiable": checks["module_markers_unverifiable"],
            "replicate_gate": "not applicable (single run)",
        },
        "code_sha": run.code_sha,
        "configuration_id": run.configuration_id,
        "reference_started": run.reference_started,
        "finished_at": run.finished_at,
        "seconds": run.seconds,
        "histories": run.histories,
        "seeds": run.case_doc.get("seeds"),
        "review": {
            "review_id": review["review_id"],
            "finished_at_at_code_sha": chron["review_finished_at"],
            "review_sha256_at_code_sha": chron["review_sha256_at_code_sha"],
            "review_sha256_at_head": chron["review_sha256_at_head"],
            "plan_sha256_at_code_sha": chron["plan_sha256_at_code_sha"],
        },
        "plan": plan,
        "manifest_bundle_set_digest": chron["manifest_bundle_set_digest"],
        "files_sha256": dict(sorted(run.read.items())),
        "outputs_note": "scorer outputs (work/) not read; output_files_sha256 recomputed from the "
        "artifact listing's work/ hashes and sizes",
        "trust_file": trust_info,
        "evaluator": {
            "head_sha": reader.head,
            "committed_inputs": dict(sorted(reader.reads.items())),
        },
    }


def _json_safe(v: Any) -> Any:
    if isinstance(v, dict):
        return {k: _json_safe(x) for k, x in v.items()}
    if isinstance(v, list | tuple):
        return [_json_safe(x) for x in v]
    if isinstance(v, np.generic):
        v = v.item()
    if isinstance(v, float) and not math.isfinite(v):
        return None
    if isinstance(v, np.ndarray):
        raise TypeError("raw arrays are never written to the report")
    return v


def _head(repo: Path) -> dict[str, Any]:
    try:
        sha = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()  # noqa: S603,S607
    except (OSError, subprocess.CalledProcessError):
        sha = None
    return {"evaluator_head": sha, "generated_at": datetime.now(UTC).isoformat()}


_FAIL_CLOSED = (
    ChronologyError, io.RecordError, st.StatisticsError, OSError, ValueError, KeyError, TypeError,
)  # fmt: skip


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="V3-008 r2 TOPAS-agreement evaluator")
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("grade")
    g.add_argument(
        "--records",
        type=Path,
        default=Path(io.DEFAULT_RECORDS),
        help=f"committed evidence directory inside the evaluator repository "
        f"(default {io.DEFAULT_RECORDS})",
    )
    g.add_argument("--seal", type=Path, default=Path(io.SEAL_RECORD_PATH))
    g.add_argument("--case", required=True, type=int)
    g.add_argument("--row", required=True, choices=ROWS)
    g.add_argument("--outputs", type=Path, default=None)
    g.add_argument("--samples-root", type=Path, default=None)
    g.add_argument("--out", type=Path, default=None)
    c = sub.add_parser("check-run", help="run validity of one signed run (evidence or smoke)")
    c.add_argument("--run-id", required=True)
    c.add_argument("--bundle", required=True, help="bundle directory name")
    c.add_argument(
        "--records-dir",
        type=Path,
        default=None,
        help=f"committed directory holding <run_id>/record.json and artifacts.json (default "
        f"{io.DEFAULT_RECORDS}/runs for an evidence part, {SMOKE_RECORDS} for a smoke)",
    )
    c.add_argument("--outputs", type=Path, default=None)
    c.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    if a.cmd == "check-run":
        return _main_check_run(a)
    try:
        rep = grade(a.records, a.seal, a.case, a.row, a.outputs, a.samples_root)
        code = {"PASS": 0}.get(rep["status"], 1)
    except _FAIL_CLOSED as exc:
        rep = {
            "schema": REPORT_SCHEMA,
            "row": a.row,
            "case": a.case,
            "status": "ERROR",
            "error": f"{type(exc).__name__}: {exc}",
            "note": "fail closed: nothing graded; not a pass",
        }
        print(f"r2_evaluate: ERROR (fail closed): {exc}", file=sys.stderr)
        code = 2
    _write_report(rep, a.out)
    return code


def _main_check_run(a: argparse.Namespace) -> int:
    try:
        rep = check_run(a.run_id, a.bundle, a.records_dir, a.outputs)
        code = 0 if rep["verdict"] == "VALID" else 1
    except _FAIL_CLOSED as exc:
        rep = {
            "schema": CHECK_SCHEMA,
            "run_id": a.run_id,
            "bundle": a.bundle,
            "verdict": "ERROR",
            "error": f"{type(exc).__name__}: {exc}",
            "note": "fail closed: run validity not established",
        }
        print(f"r2_evaluate check-run: ERROR (fail closed): {exc}", file=sys.stderr)
        code = 2
    _write_report(rep, a.out)
    return code


def _write_report(rep: dict[str, Any], out: Path | None) -> None:
    try:
        own = evaluator_repo()
    except io.RecordError:
        own = Path(__file__).resolve().parent
    rep = _json_safe(rep | _head(own))
    text = json.dumps(rep, indent=1, sort_keys=True, allow_nan=False)
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)


if __name__ == "__main__":
    sys.exit(main())
