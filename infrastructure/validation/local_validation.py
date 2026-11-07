#!/usr/bin/env python3

import argparse
import json
import os
import re
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path


WORKTREE_ROOT = Path(__file__).resolve().parents[2].parent / (Path(__file__).resolve().parents[2].name + "-worktrees")

VALIDATION_ROOT = Path(
    os.environ.get(
        "IONMC_VALIDATION_DIR",
        str(Path.home() / ".local/share/ionmc-experiment/validation"),
    )
).expanduser()


def git_text(*args: str, cwd: Path) -> str:
    return subprocess.check_output(
        ["git", *args],
        cwd=cwd,
        text=True,
        stderr=subprocess.PIPE,
    ).strip()


def normalize_task_id(task_id: str) -> str:
    value = task_id.strip().upper()

    if not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{1,63}", value):
        raise ValueError(f"Invalid task ID: {task_id}")

    return value


def worktree_for(task_id: str) -> Path:
    return WORKTREE_ROOT / normalize_task_id(task_id).lower()


def worktree_state(task_id: str) -> dict:
    path = worktree_for(task_id)

    if not path.exists():
        raise RuntimeError(f"Task worktree does not exist: {path}")

    status = git_text("status", "--porcelain", cwd=path)

    return {
        "path": path,
        "branch": git_text("branch", "--show-current", cwd=path),
        "sha": git_text("rev-parse", "HEAD", cwd=path),
        "dirty": bool(status),
        "porcelain_status": status,
    }


SIGNER = "validation-manager"

# Trust map consulted by the gate. Deliberately not overridable by environment
# variable (tests monkeypatch this attribute).
TRUSTED_KEYS_PATH = Path(__file__).resolve().parents[1] / "keys" / "trusted_attestation_keys.json"

BOOTSTRAP_NOTE = "unverified (no trusted key published yet)"


def _attest():
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from infrastructure.experiment_v3 import attest

    return attest


def load_trust() -> dict:
    return _attest().load_trusted(TRUSTED_KEYS_PATH)


def enforcing(trusted: dict) -> bool:
    """True once the trust map lists at least one key for the signer."""
    return bool(trusted.get(SIGNER))


def classify(rec: dict, task_id: str, sha: str, trusted: dict) -> tuple[str, str]:
    """Return (attestation_status, reason): verified | unverified | invalid.

    Besides the signature, the signed identity must match the requested task
    and SHA, so a genuine record copied into another task/SHA directory is
    invalid.
    """
    status, reason = _attest().classify_record(rec, trusted, SIGNER)
    if status == "invalid":
        return status, reason
    if rec.get("task_id") != normalize_task_id(task_id) or rec.get("git_sha") != sha:
        return "invalid", (
            f"signed identity mismatch: record is for task {rec.get('task_id')!r} "
            f"sha {rec.get('git_sha')!r}, requested {normalize_task_id(task_id)!r} {sha!r}"
        )
    return status, reason


def records_for(task_id: str, sha: str) -> list[dict]:
    """Load records, each annotated with derived attestation_status/_reason.

    The derived fields are added after verification and are not part of the
    signed payload; they are never read from disk (any such field in a file is
    overwritten).
    """
    directory = VALIDATION_ROOT / normalize_task_id(task_id) / sha

    if not directory.exists():
        return []

    trusted = load_trust()
    records = []

    for path in sorted(directory.glob("*.json")):
        try:
            rec = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(rec, dict):
                continue
        except Exception:
            continue
        rec.pop("attestation_status", None)
        rec.pop("attestation_reason", None)
        status, reason = classify(rec, task_id, sha, trusted)
        rec["attestation_status"] = status
        rec["attestation_reason"] = reason
        records.append(rec)

    return records


def evaluate(records: list[dict], dirty: bool) -> dict:
    """Gate decision over annotated records.

    Enforcement rule (fail closed once trust exists):
    - If the trust map lists at least one key for "validation-manager", only
      records with attestation_status "verified" count as passing; unverified
      and invalid records are rejected with an explicit reason.
    - "invalid" records (attestation present but wrong, or signed identity not
      matching the requested task/SHA) are always rejected, even in bootstrap.
    - Bootstrap exception: while the trust map has no key for the signer
      (before the operator publishes it) unverified records are tolerated and
      the summary reports them as "unverified (no trusted key published yet)".
      The exception ends as soon as a key is published.
    - Every record with status "failed" counts as failing regardless of its
      attestation (a forged failure can only block, never release, the gate).
    """
    enforce = enforcing(load_trust())
    passing, rejected = [], []

    for r in records:
        att = r.get("attestation_status", "unverified")
        if att == "invalid" or (enforce and att != "verified"):
            rejected.append(
                {
                    "validation_id": r.get("validation_id"),
                    "attestation_status": att,
                    "reason": r.get("attestation_reason"),
                }
            )
        elif r.get("status") == "passed":
            passing.append(r)

    failing = [r for r in records if r.get("status") == "failed"]
    unverified = [r for r in records if r.get("attestation_status") == "unverified"]

    summary = {
        "enforcement": "enforcing" if enforce else "bootstrap",
        "verified_records": sum(
            1 for r in records if r.get("attestation_status") == "verified"
        ),
        "unverified_records": len(unverified),
        "invalid_records": sum(
            1 for r in records if r.get("attestation_status") == "invalid"
        ),
        "rejected_records": rejected,
        "passing_records": len(passing),
        "failing_records": len(failing),
        "gate_passed": not dirty and len(passing) > 0 and len(failing) == 0,
    }
    if not enforce and unverified:
        summary["attestation_note"] = (
            f"{len(unverified)} record(s) {BOOTSTRAP_NOTE}; tolerated only "
            "until the operator publishes the validation-manager key"
        )
    return summary


def record(
    task_id: str,
    suite: str,
    status: str,
    summary: str,
    details_json: str,
):
    state = worktree_state(task_id)

    if state["dirty"]:
        raise RuntimeError(
            "Refusing to record validation for a dirty worktree. "
            "Commit the exact state being validated first."
        )

    try:
        details = json.loads(details_json)
    except json.JSONDecodeError as exc:
        raise ValueError("--details-json must contain valid JSON") from exc

    now = datetime.now(timezone.utc)

    validation_id = (
        f"VAL-{now.strftime('%Y%m%d-%H%M%S')}-"
        f"{uuid.uuid4().hex[:6].upper()}"
    )

    record = {
        "schema_version": 1,
        "validation_id": validation_id,
        "task_id": normalize_task_id(task_id),
        "recorded_at_utc": now.isoformat(),
        "suite": suite,
        "status": status,
        "summary": summary,
        "details": details,
        "git_branch": state["branch"],
        "git_sha": state["sha"],
        "worktree": str(state["path"]),
    }

    # V3-INFRA-002: Ed25519 attestation by the service-held key (additive field).
    record = _attest().sign_with_service_key(record, SIGNER)

    directory = VALIDATION_ROOT / normalize_task_id(task_id) / state["sha"]
    directory.mkdir(parents=True, exist_ok=True)

    path = directory / f"{validation_id}.json"

    path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "status": "recorded",
                "validation": record,
                "record_file": str(path),
            },
            indent=2,
        )
    )


def inspect(task_id: str):
    state = worktree_state(task_id)
    records = records_for(task_id, state["sha"])
    verdict = evaluate(records, state["dirty"])

    result = {
        "status": "present" if records else "missing",
        "task_id": normalize_task_id(task_id),
        "git_branch": state["branch"],
        "git_sha": state["sha"],
        "worktree_dirty": state["dirty"],
        "records": records,
        **verdict,
    }

    print(json.dumps(result, indent=2))


def check(task_id: str):
    state = worktree_state(task_id)
    records = records_for(task_id, state["sha"])
    verdict = evaluate(records, state["dirty"])

    result = {
        "status": "passed" if verdict["gate_passed"] else "failed",
        "task_id": normalize_task_id(task_id),
        "git_sha": state["sha"],
        "worktree_dirty": state["dirty"],
        **{k: v for k, v in verdict.items() if k != "gate_passed"},
    }

    print(json.dumps(result, indent=2))

    if not verdict["gate_passed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    p_record = sub.add_parser("record")
    p_record.add_argument("--task-id", required=True)
    p_record.add_argument("--suite", required=True)
    p_record.add_argument(
        "--status",
        choices=["passed", "failed"],
        required=True,
    )
    p_record.add_argument("--summary", required=True)
    p_record.add_argument("--details-json", default="{}")

    p_inspect = sub.add_parser("inspect")
    p_inspect.add_argument("--task-id", required=True)

    p_check = sub.add_parser("check")
    p_check.add_argument("--task-id", required=True)

    args = parser.parse_args()

    if args.command == "record":
        record(
            args.task_id,
            args.suite,
            args.status,
            args.summary,
            args.details_json,
        )

    elif args.command == "inspect":
        inspect(args.task_id)

    elif args.command == "check":
        check(args.task_id)


if __name__ == "__main__":
    main()
