"""V3-008 revision 2: run validity of the signed smoke re-runs (plan amendment 2026-11-09, Codex
REVIEW-8429c50e: the re-runs "must pass ``r2_io.run_checks`` and the R2-TOPAS-4 chronology against
the seal review").

Usage::

    python validation/scripts/ion_em/r2_smoke_reruns.py --ledger <ledger.jsonl> [--ledger ...]
        [--outputs <dir>] [--out validation/results/transport/v8r2-topas/smoke-reruns.json]
        [--confirmation-set-complete]

``--ledger`` is repeatable (e.g. the R2-SMOKE-LP ledger plus the re-runs made at the current
head); a run id listed twice is refused. ``--confirmation-set-complete`` asks the summary to state
``confirmation_set_complete: true``; it is set only when every run of the given ledgers is VALID,
otherwise the summary records ``false`` with the reason and the exit status is 1. The summary
records the real evaluator HEAD (the checkout that contains this script; every run's
``check_run`` must have read its committed inputs at that same HEAD).

Each ledger has one JSON object per line (``case`` = bundle directory name, ``seed``, ``run_id``,
``code_sha``, ``started_at``, ``exec_real_s``, ``previous_unsigned.run_id``); it is the
orchestrator's execution ledger and is NOT evidence: every verdict comes from
``r2_evaluate.check_run`` on the signed record committed under
``validation/results/transport/v8r2-topas/smoke-runs/<run_id>/`` (the ledger's ``code_sha`` and
``started_at`` are compared with the signed record and a difference is an ERROR). Smokes are
execution records, never evidence. Exit 0 iff every run is VALID.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r2_evaluate as ev  # noqa: E402

SCHEMA = "ionmc-v3-008-r2-smoke-reruns-1"
DEFAULT_OUT = "validation/results/transport/v8r2-topas/smoke-reruns.json"


def check_entry(entry: dict[str, Any], outputs: Path | None) -> dict[str, Any]:
    row: dict[str, Any] = {
        "case": entry["case"], "seed": entry["seed"], "run_id": entry["run_id"],
        "exec_real_s": entry.get("exec_real_s"),
        "previous_unsigned_run_id": (entry.get("previous_unsigned") or {}).get("run_id"),
    }  # fmt: skip
    try:
        rep = ev.check_run(entry["run_id"], entry["case"], None, outputs)
        for key in ("code_sha", "started_at"):
            signed = rep[key] if key == "code_sha" else rep["reference_started"]
            if entry.get(key) != signed:
                raise ev.io.RecordError(f"ledger {key} {entry.get(key)!r} != signed {signed!r}")
        if rep["seeds"] != [entry["seed"]]:
            raise ev.io.RecordError(f"ledger seed {entry['seed']} != bundle seeds {rep['seeds']}")
    except ev._FAIL_CLOSED as exc:
        return row | {"verdict": "ERROR", "error": f"{type(exc).__name__}: {exc}"}
    return row | {
        "code_sha": rep["code_sha"], "started_at": rep["reference_started"],
        "seconds": rep["seconds"], "histories": rep["histories"], "config": rep["config"],
        "verdict": rep["verdict"], "checks": rep["checks"], "review": rep["review"],
        "trust_file_sha256": rep["trust_file"]["sha256"],
        "evaluator_head": rep["evaluator"]["head_sha"],
    }  # fmt: skip


def read_ledgers(paths: list[Path]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """(entries in ledger order, provenance ``[{"path", "sha256", "n_runs"}]``); a run id that
    occurs twice (within or across ledgers) raises RecordError."""
    entries: list[dict[str, Any]] = []
    prov = []
    for path in paths:
        blob = path.read_bytes()
        rows = [json.loads(ln) for ln in blob.decode("utf-8").splitlines() if ln.strip()]
        prov.append({"path": str(path), "sha256": hashlib.sha256(blob).hexdigest(),
                     "n_runs": len(rows)})  # fmt: skip
        entries += rows
    ids = [e.get("run_id") for e in entries]
    dup = sorted({i for i in ids if ids.count(i) > 1}, key=str)
    if dup:
        raise ev.io.RecordError(f"run id(s) listed more than once across the ledgers: {dup}")
    return entries, prov


def confirmation(rows: list[dict[str, Any]], requested: bool) -> dict[str, Any]:
    """``confirmation_set_complete`` is true only if requested AND every run is VALID."""
    bad = [f"{r['run_id']} ({r['case']}): {r['verdict']}" for r in rows if r["verdict"] != "VALID"]
    if not requested:
        reason = "not requested (--confirmation-set-complete not given)"
    elif not rows:
        reason = "refused: no runs"
    elif bad:
        reason = f"refused: {len(bad)} run(s) not VALID: {bad}"
    else:
        reason = f"all {len(rows)} runs of the given ledgers are VALID"
    complete = requested and bool(rows) and not bad
    return {"confirmation_set_complete": complete, "confirmation_reason": reason}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ledger", type=Path, action="append", required=True)
    ap.add_argument("--outputs", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=Path(DEFAULT_OUT))
    ap.add_argument("--confirmation-set-complete", action="store_true")
    a = ap.parse_args(argv)
    head = ev.io.CommittedReader(ev.evaluator_repo()).head
    entries, ledgers = read_ledgers(a.ledger)
    rows = [check_entry(e, a.outputs) for e in entries]
    reviews = {(r["review"]["review_id"], r["review"]["finished_at_at_code_sha"]) for r in rows
               if "review" in r}  # fmt: skip
    trusts = {r["trust_file_sha256"] for r in rows if "trust_file_sha256" in r}
    for r in rows:
        if r.get("evaluator_head", head) != head:  # HEAD moved while checking: fail closed
            r.update(verdict="ERROR", error=f"checked at HEAD {r['evaluator_head']}, not {head}")
        for key in ("review", "trust_file_sha256", "evaluator_head"):
            r.pop(key, None)
    conf = confirmation(rows, a.confirmation_set_complete)
    doc = {
        "schema": SCHEMA,
        "task_id": "V3-008",
        "purpose": "run validity of signed revision-2 smoke runs (plan amendment 2026-11-09, "
        "REVIEW-8429c50e): r2_io.run_checks + R2-TOPAS-4 chronology against the seal review; "
        "execution records only, never evidence",
        "checker": "validation/scripts/ion_em/r2_evaluate.py check-run (r2_evaluate.check_run)",
        "ledgers": ledgers,
        "scorer_outputs_read": False,
        "outputs_note": "only request.json, stdout.txt, stderr.txt were materialized and read; the "
        "scorer outputs are bound through the signed output_files_sha256 recomputed from the "
        "artifact listing's work/ hashes",
        "seal_review": [{"review_id": i, "finished_at": f} for i, f in sorted(reviews)],
        "trust_file": {"path": ev.io.DEFAULT_TRUST, "sha256": sorted(trusts)},
        "evaluator_head": head,
        "n_runs": len(rows),
        "all_valid": bool(rows) and all(r["verdict"] == "VALID" for r in rows),
        **conf,
        "runs": rows,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    for r in rows:
        print(f"{r['run_id']} {r['case']}: {r['verdict']}")
    print(f"confirmation_set_complete: {conf['confirmation_set_complete']} "
          f"({conf['confirmation_reason']})")  # fmt: skip
    ok = doc["all_valid"] and (conf["confirmation_set_complete"] or not a.confirmation_set_complete)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
