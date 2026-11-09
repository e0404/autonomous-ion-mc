"""V3-008 revision 2: run validity of the signed smoke re-runs (plan amendment 2026-11-09, Codex
REVIEW-8429c50e: the re-runs "must pass ``r2_io.run_checks`` and the R2-TOPAS-4 chronology against
the seal review").

Usage::

    python validation/scripts/ion_em/r2_smoke_reruns.py --ledger <ledger-signed.jsonl>
        [--outputs <dir>] [--out validation/results/transport/v8r2-topas/smoke-reruns.json]

The ledger has one JSON object per line (``case`` = bundle directory name, ``seed``, ``run_id``,
``code_sha``, ``started_at``, ``exec_real_s``, ``previous_unsigned.run_id``); it is the
orchestrator's execution ledger and is NOT evidence: every verdict comes from
``r2_evaluate.check_run`` on the signed record committed under
``validation/results/transport/v8r2-topas/smoke-runs/<run_id>/`` (the ledger's ``code_sha`` and
``started_at`` are compared with the signed record and a difference is an ERROR). Smokes are
execution records, never evidence. Exit 0 iff every run is VALID.
"""

from __future__ import annotations

import argparse
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


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ledger", type=Path, required=True)
    ap.add_argument("--outputs", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=Path(DEFAULT_OUT))
    a = ap.parse_args(argv)
    entries = [json.loads(ln) for ln in a.ledger.read_text().splitlines() if ln.strip()]
    rows = [check_entry(e, a.outputs) for e in entries]
    reviews = {(r["review"]["review_id"], r["review"]["finished_at_at_code_sha"]) for r in rows
               if "review" in r}  # fmt: skip
    trusts = {r["trust_file_sha256"] for r in rows if "trust_file_sha256" in r}
    heads = {r["evaluator_head"] for r in rows if "evaluator_head" in r}
    for r in rows:
        for key in ("review", "trust_file_sha256", "evaluator_head"):
            r.pop(key, None)
    doc = {
        "schema": SCHEMA,
        "task_id": "V3-008",
        "purpose": "run validity of the 14 signed smoke re-runs (plan amendment 2026-11-09, "
        "REVIEW-8429c50e): r2_io.run_checks + R2-TOPAS-4 chronology against the seal review; "
        "execution records only, never evidence",
        "checker": "validation/scripts/ion_em/r2_evaluate.py check-run (r2_evaluate.check_run)",
        "scorer_outputs_read": False,
        "outputs_note": "only request.json, stdout.txt, stderr.txt were materialized and read; the "
        "scorer outputs are bound through the signed output_files_sha256 recomputed from the "
        "artifact listing's work/ hashes",
        "seal_review": [{"review_id": i, "finished_at": f} for i, f in sorted(reviews)],
        "trust_file": {"path": ev.io.DEFAULT_TRUST, "sha256": sorted(trusts)},
        "evaluator_head": sorted(heads),
        "n_runs": len(rows),
        "all_valid": bool(rows) and all(r["verdict"] == "VALID" for r in rows),
        "runs": rows,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    for r in rows:
        print(f"{r['run_id']} {r['case']}: {r['verdict']}")
    return 0 if doc["all_valid"] else 1


if __name__ == "__main__":
    sys.exit(main())
