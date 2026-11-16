"""Production consumer of signed codex-worker (REVIEW-...) review records.

The review service (infrastructure/experiment_v3/review.py, V3-INFRA-006)
writes ``<state>/reviews/<REVIEW-id>/record.json`` (schema_version 1) when a
review finishes, with an Ed25519 attestation (signer "codex-worker"). The
record is the dict inspect_codex_review returns: review_id, task_id, head_sha,
base_sha, status, started_at, finished_at (both the service's own UTC
timestamps), model, provider, tool, experiment_id, and for a finished review
with a structured report also report_sha256 and the parsed report (a failed
review carries ``error`` instead).

This module loads such a record and verifies it: signature by the
"codex-worker" key, schema, status/verdict consistency and timestamps, and
optionally the identity the caller expects (task id, head SHA, base SHA) and
the archived report.json bytes. Status values follow
infrastructure/experiment_v3/attest.py: verified, unverified (no attestation,
or the codex-worker key is not in the trust map, e.g. before its key is
published), invalid (anything wrong). Whether "unverified" is tolerated is the
caller's decision; chronology gates should fail closed on anything but
"verified".

CLI: ``python -m infrastructure.validation.review_records <record.json>``
prints ``{"status", "reason", "review_id", ...}`` and exits 0 only when the
record is verified.
"""

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

SIGNER = "codex-worker"
SCHEMA_VERSION = 1
REVIEW_ID_PATTERN = re.compile(r"REVIEW-[0-9a-f]{32}")
SHA_PATTERN = re.compile(r"[0-9a-f]{40}")
DIGEST_PATTERN = re.compile(r"[0-9a-f]{64}")
STATUSES = ("passed", "changes_required", "failed")
REQUIRED_STRINGS = ("task_id", "model", "provider", "tool", "experiment_id")


class ReviewRecordMissing(Exception):
    """No readable signed record exists at the given path or review id."""


def _repo_import():
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def _attest():
    _repo_import()
    from infrastructure.experiment_v3 import attest

    return attest


REVIEW_STATE = None  # tests monkeypatch; None means common.STATE


def reviews_root(state=None) -> Path:
    """``<state>/reviews``; default state as the review service uses it."""
    if state is None:
        state = REVIEW_STATE
    if state is None:
        _repo_import()
        from infrastructure.experiment_v3.common import STATE

        state = STATE
    return Path(state) / "reviews"


def load_review_record(path) -> dict:
    """Load a signed review ``record.json`` (raises ReviewRecordMissing)."""
    path = Path(path)
    if not path.is_file():
        raise ReviewRecordMissing(f"review record {path}: not found")
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ReviewRecordMissing(f"review record {path} unreadable: {exc}") from exc
    if not isinstance(record, dict):
        raise ReviewRecordMissing(f"review record {path} is not an object")
    return record


def review_record_path(review_id: str, state=None) -> Path:
    if not isinstance(review_id, str) or not REVIEW_ID_PATTERN.fullmatch(review_id):
        raise ReviewRecordMissing(f"invalid review id: {review_id!r}")
    return reviews_root(state) / review_id / "record.json"


def _utc_time(value):
    """Timezone-aware datetime from an ISO-8601 string, else None."""
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return None
    return parsed if parsed.tzinfo is not None else None


def current_experiment_id() -> str:
    """The v3 experiment identifier (infrastructure/experiment_v3/paths.py).

    The review service signs ``experiment_id`` from the same
    ``.ionmc-condition.json`` that defines this default.
    """
    _repo_import()
    from infrastructure.experiment_v3.paths import EXPERIMENT_ID

    return EXPERIMENT_ID


def _report_passes(report, head_sha):
    """review.validate_report (raises ValueError on a malformed report)."""
    _repo_import()
    from infrastructure.experiment_v3.review import validate_report

    return validate_report(report, head_sha)


def schema_problem(record: dict) -> str | None:
    """Reason the signed payload is malformed or inconsistent, else None."""
    if record.get("schema_version") != SCHEMA_VERSION:
        return f"schema_version is {record.get('schema_version')!r}, expected {SCHEMA_VERSION}"
    if not REVIEW_ID_PATTERN.fullmatch(str(record.get("review_id"))):
        return f"review_id is malformed: {record.get('review_id')!r}"
    for field in ("head_sha", "base_sha"):
        if not SHA_PATTERN.fullmatch(str(record.get(field))):
            return f"{field} is not a 40-hex commit SHA: {record.get(field)!r}"
    for field in REQUIRED_STRINGS:
        value = record.get(field)
        if not isinstance(value, str) or not value.strip():
            return f"{field} is not a non-blank string"
    status = record.get("status")
    if status not in STATUSES:
        return f"status is {status!r}, expected one of {list(STATUSES)}"
    started = _utc_time(record.get("started_at"))
    if started is None:
        return (
            "started_at missing or not a timezone-aware ISO-8601 time: "
            f"{record.get('started_at')!r}"
        )
    finished = _utc_time(record.get("finished_at"))
    if finished is None or not started < finished:
        return (
            f"finished_at {record.get('finished_at')!r} is missing, malformed or "
            f"not later than started_at {record.get('started_at')!r}"
        )
    if status == "failed":
        if not isinstance(record.get("error"), str) or not record["error"].strip():
            return "failed review without an error"
        if "report" in record or "report_sha256" in record:
            return "failed review carries a report"
        return None
    if not DIGEST_PATTERN.fullmatch(str(record.get("report_sha256"))):
        return "report_sha256 missing or not 64 hex digits"
    try:
        passes = _report_passes(record.get("report"), record["head_sha"])
    except ValueError as exc:
        return f"report invalid: {exc}"
    if passes != (status == "passed"):
        return (
            f"status {status!r} inconsistent with report verdict "
            f"{record['report'].get('verdict')!r} and its findings"
        )
    return None


def verify_review_record(
    record: dict,
    trusted: dict | None = None,
    *,
    expected_task_id: str | None = None,
    expected_head_sha: str | None = None,
    expected_base_sha: str | None = None,
    expected_experiment_id: str | None = None,
    report_path=None,
) -> tuple[str, str]:
    """Return (status, reason): verified | unverified | invalid.

    ``trusted`` is ``{signer: {key_id: public_key_hex}}`` (default: the
    committed trust file). Signature/trust problems keep their attest
    classification (an unsigned record or an untrusted key is "unverified").
    A malformed or inconsistent payload (schema_problem), an identity
    differing from an expected value (replay protection), or a report.json at
    ``report_path`` whose sha256 differs from the signed report_sha256 makes
    the record "invalid" even when the signature is genuine. The signed
    ``experiment_id`` must always equal ``expected_experiment_id`` (default:
    the current v3 experiment, current_experiment_id), so a genuine record of
    another experiment is never accepted.
    """
    if not isinstance(record, dict):
        return "invalid", "review record is not an object"
    attest = _attest()
    if trusted is None:
        trusted = attest.load_trusted()
    status, reason = attest.classify_record(record, trusted, SIGNER)
    if status == "invalid":
        return status, reason
    problem = schema_problem(record)
    if problem:
        return "invalid", problem
    expected = {
        "task_id": expected_task_id,
        "head_sha": expected_head_sha,
        "base_sha": expected_base_sha,
        "experiment_id": expected_experiment_id or current_experiment_id(),
    }
    for field, want in expected.items():
        if want is not None and record.get(field) != want:
            return "invalid", (
                f"signed identity mismatch: {field} is {record.get(field)!r}, expected {want!r}"
            )
    if report_path is not None and record.get("report_sha256"):
        try:
            actual = hashlib.sha256(Path(report_path).read_bytes()).hexdigest()
        except OSError as exc:
            return "invalid", f"report.json unreadable: {exc}"
        if actual != record["report_sha256"]:
            return "invalid", "report.json does not match the signed report_sha256"
    return status, reason


def review_evidence(task_id: str, head_sha: str, trusted: dict, state=None) -> list[dict]:
    """Signed review records of ``task_id`` at ``head_sha`` under the review state.

    Each entry: {review_id, status, reason, review_status, ok}; ``ok`` requires
    a verified record (trusted codex-worker key) bound to this task and SHA
    whose review passed. Used by the local-validation gate as key-publication
    evidence for a codex-worker key (V3-INFRA-006).
    """
    out = []
    root = reviews_root(state)
    for path in sorted(root.glob("REVIEW-*/record.json")):
        try:
            record = load_review_record(path)
        except ReviewRecordMissing as exc:
            out.append(
                {
                    "review_id": path.parent.name,
                    "status": "missing",
                    "reason": str(exc),
                    "review_status": None,
                    "ok": False,
                }
            )
            continue
        if record.get("task_id") != task_id or record.get("head_sha") != head_sha:
            continue
        status, reason = verify_review_record(
            record,
            trusted,
            expected_task_id=task_id,
            expected_head_sha=head_sha,
            report_path=path.parent / "report.json",
        )
        out.append(
            {
                "review_id": path.parent.name,
                "status": status,
                "reason": reason,
                "review_status": record.get("status"),
                "ok": status == "verified" and record.get("status") == "passed",
            }
        )
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Verify a signed Codex review record.")
    parser.add_argument("record", type=Path, help="path to a review record.json")
    parser.add_argument("--trust-file", type=Path, default=None)
    parser.add_argument("--task-id", default=None)
    parser.add_argument("--head-sha", default=None)
    parser.add_argument("--base-sha", default=None)
    parser.add_argument("--experiment-id", default=None)
    args = parser.parse_args(argv)
    try:
        record = load_review_record(args.record)
    except ReviewRecordMissing as exc:
        print(json.dumps({"status": "invalid", "reason": str(exc)}))
        return 1
    trusted = _attest().load_trusted(args.trust_file)
    report = args.record.parent / "report.json"
    status, reason = verify_review_record(
        record,
        trusted,
        expected_task_id=args.task_id,
        expected_head_sha=args.head_sha,
        expected_base_sha=args.base_sha,
        expected_experiment_id=args.experiment_id,
        report_path=report if report.is_file() else None,
    )
    raw = record.get("attestation")
    block: dict = raw if isinstance(raw, dict) else {}
    print(
        json.dumps(
            {
                "status": status,
                "reason": reason,
                "review_id": record.get("review_id"),
                "task_id": record.get("task_id"),
                "head_sha": record.get("head_sha"),
                "review_status": record.get("status"),
                "started_at": record.get("started_at"),
                "finished_at": record.get("finished_at"),
                "key_id": block.get("key_id"),
                "signed_at_utc": block.get("signed_at_utc"),
            },
            indent=2,
        )
    )
    return 0 if status == "verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
