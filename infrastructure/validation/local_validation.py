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

# Tasks for which bootstrap tolerance (unverified records counted while the
# signer has no trusted key) applies: only the migration task that introduces
# attestations. Every other task fails closed whatever the trust map state.
BOOTSTRAP_TASKS = frozenset({"V3-INFRA-002"})

ATTESTATION_REQUIRED = (
    "attestation required: no trusted validation-manager key published "
    "(operator runbook) or record not verified"
)
HOST_RUN_SIGNER = "host-runner"


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
        rec.pop("host_run_attestation", None)
        status, reason = classify(rec, task_id, sha, trusted)
        rec["attestation_status"] = status
        rec["attestation_reason"] = reason
        rec["host_run_attestation"] = check_host_runs(rec, task_id, sha, trusted)
        records.append(rec)

    return records


def bootstrap_allowed(task_id: str, trusted: dict, signer: str = SIGNER) -> bool:
    """Bootstrap tolerance: signer has no trusted key AND task is a migration task."""
    return not trusted.get(signer) and normalize_task_id(task_id) in BOOTSTRAP_TASKS


def host_run_ids(rec: dict) -> list[str]:
    """Host run ids referenced by a validation record's details."""
    details = rec.get("details")
    if not isinstance(details, dict):
        return []
    ids = []
    value = details.get("host_run_id")
    if isinstance(value, str) and value.strip():
        ids.append(value.strip())
    many = details.get("host_run_ids")
    if isinstance(many, list):
        ids.extend(v.strip() for v in many if isinstance(v, str) and v.strip())
    return ids


def check_host_runs(rec: dict, task_id: str, sha: str, trusted: dict) -> list[dict]:
    """Verify every signed run record a validation record references.

    Each entry: {run_id, status, reason, ok}. ``ok`` is False for a missing,
    invalid or (when not tolerated) unverified run record, a task/SHA that
    differs from the validation record's, a non-zero exit code or a timeout.
    """
    from infrastructure.validation import host_run_records as hrr

    tolerated = bootstrap_allowed(task_id, trusted, HOST_RUN_SIGNER)
    out = []
    for run_id in host_run_ids(rec):
        entry = {"run_id": run_id, "status": "missing", "reason": "", "ok": False}
        out.append(entry)
        try:
            run = hrr.load_signed_run_record(run_id, normalize_task_id(task_id))
        except hrr.RunRecordMissing as exc:
            entry["reason"] = f"host run record missing: {exc}"
            continue
        status, reason = hrr.verify_run_record(
            run,
            rec.get("task_id", ""),
            rec.get("git_sha", ""),
            run_id,
            trusted,
        )
        entry["status"], entry["reason"] = status, reason
        if status == "invalid":
            continue
        if status != "verified" and not tolerated:
            entry["reason"] = f"host run record {status} (attestation required): {reason}"
            continue
        if run.get("exit_code") != 0 or run.get("timed_out") is not False:
            entry["reason"] = (
                f"host run did not succeed: exit_code={run.get('exit_code')!r} "
                f"timed_out={run.get('timed_out')!r}"
            )
            continue
        entry["ok"] = True
    return out


def evaluate(records: list[dict], dirty: bool, task_id: str) -> dict:
    """Gate decision over annotated records for ``task_id``.

    Enforcement rule (fail closed):
    - Only records with attestation_status "verified" count as passing, except
      in bootstrap. Unverified and invalid records are rejected with a reason.
    - Bootstrap (unverified records counted, noted as "unverified (no trusted
      key published yet)") applies ONLY when task_id is in BOOTSTRAP_TASKS
      (the migration task V3-INFRA-002) AND the trust map has no key for
      "validation-manager". For every other task unverified records are
      rejected whatever the trust map state, with ATTESTATION_REQUIRED; all
      other tasks' gates are therefore blocked until the operator publishes the
      keys.
    - "invalid" records (present but wrong, or signed identity not matching the
      requested task/SHA) are always rejected, even in bootstrap.
    - A passed record that references host runs (details.host_run_id) is
      rejected unless every referenced signed run record is present, bound to
      the same task/SHA, verified (or bootstrap-tolerated, see check_host_runs)
      and succeeded (exit_code 0, not timed out).
    - Every record with status "failed" counts as failing regardless of its
      attestation (a forged failure can only block, never release, the gate).
    """
    trusted = load_trust()
    enforce = enforcing(trusted)
    tolerate = bootstrap_allowed(task_id, trusted)
    passing, rejected = [], []

    for r in records:
        att = r.get("attestation_status", "unverified")
        runs = r.get("host_run_attestation") or []
        bad_run = next((x for x in runs if not x.get("ok")), None)
        if att == "invalid" or (att != "verified" and not tolerate):
            reason = r.get("attestation_reason")
            if att == "unverified" and not enforce:
                reason = ATTESTATION_REQUIRED
            rejected.append(
                {
                    "validation_id": r.get("validation_id"),
                    "attestation_status": att,
                    "reason": reason,
                }
            )
        elif r.get("status") == "passed" and bad_run is not None:
            rejected.append(
                {
                    "validation_id": r.get("validation_id"),
                    "attestation_status": att,
                    "reason": f"host run {bad_run['run_id']}: {bad_run['reason']}",
                }
            )
        elif r.get("status") == "passed":
            passing.append(r)

    failing = [r for r in records if r.get("status") == "failed"]
    unverified = [r for r in records if r.get("attestation_status") == "unverified"]

    summary = {
        "enforcement": "enforcing" if enforce else ("bootstrap" if tolerate else "blocked"),
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
    if tolerate and unverified:
        summary["attestation_note"] = (
            f"{len(unverified)} record(s) {BOOTSTRAP_NOTE}; tolerated only for "
            f"bootstrap task(s) {sorted(BOOTSTRAP_TASKS)} until the operator "
            "publishes the validation-manager key"
        )
    if not enforce and not tolerate:
        summary["attestation_required"] = ATTESTATION_REQUIRED
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
    verdict = evaluate(records, state["dirty"], task_id)

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
    verdict = evaluate(records, state["dirty"], task_id)

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
