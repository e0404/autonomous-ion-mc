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


LEGACY_NO_BINDING = "host run referenced without output_manifest_sha256 binding"
MANIFEST_PATTERN = re.compile(r"[0-9a-f]{64}")
_NO_CLAIM = object()


def host_run_refs(rec: dict) -> list[tuple[str, object]]:
    """Host runs referenced by a validation record as (run_id, claim).

    Convention: ``details.host_runs = [{"run_id": "RUN-...",
    "output_manifest_sha256": "<64 hex or null>"}]``. ``claim`` is the stated
    digest (str or None) when the key is present, else the ``_NO_CLAIM``
    sentinel. The legacy ``host_run_id`` / ``host_run_ids`` forms never carry a
    claim.
    """
    details = rec.get("details")
    if not isinstance(details, dict):
        return []
    refs: list[tuple[str, object]] = []
    runs = details.get("host_runs")
    if isinstance(runs, list):
        for item in runs:
            if isinstance(item, dict) and isinstance(item.get("run_id"), str):
                run_id = item["run_id"].strip()
                if run_id:
                    refs.append((run_id, item.get("output_manifest_sha256", _NO_CLAIM)))
    value = details.get("host_run_id")
    if isinstance(value, str) and value.strip():
        refs.append((value.strip(), _NO_CLAIM))
    many = details.get("host_run_ids")
    if isinstance(many, list):
        refs.extend((v.strip(), _NO_CLAIM) for v in many if isinstance(v, str) and v.strip())
    return refs


def host_run_ids(rec: dict) -> list[str]:
    """Host run ids referenced by a validation record's details."""
    return [run_id for run_id, _ in host_run_refs(rec)]


def run_success_problem(run: dict) -> str | None:
    """Reason a signed run record cannot support validation, else None."""
    checks = [
        ("succeeded", run.get("succeeded") is True),
        ("exit_code", run.get("exit_code") == 0),
        ("timed_out", run.get("timed_out") is False),
        (
            "committed_state_unchanged",
            "committed_state_unchanged" not in run or run["committed_state_unchanged"] is True,
        ),
        ("dirty_after", "dirty_after" not in run or run["dirty_after"] is False),
        ("output_publication_error", run.get("output_publication_error") is None),
    ]
    bad = [f"{name}={run.get(name)!r}" for name, ok in checks if not ok]
    if bad:
        return "host run did not succeed: " + ", ".join(bad)
    return None


def check_host_runs(rec: dict, task_id: str, sha: str, trusted: dict) -> list[dict]:
    """Verify every signed run record a validation record references.

    Each entry: {run_id, status, reason, ok}. ``ok`` is False for a missing,
    invalid or (when not tolerated) unverified run record, a task/SHA that
    differs from the validation record's, an output manifest digest that
    differs from the validation record's claim (null included; a signed record
    without the field is invalid), a legacy
    reference without a manifest claim (when not tolerated), or a run that did
    not succeed (see run_success_problem).
    """
    from infrastructure.validation import host_run_records as hrr

    tolerated = bootstrap_allowed(task_id, trusted, HOST_RUN_SIGNER)
    details = rec.get("details") if isinstance(rec.get("details"), dict) else {}
    lookup_missing = details.get("host_run_manifest_lookup") == "missing"
    out = []
    for run_id, claim in host_run_refs(rec):
        entry = {"run_id": run_id, "status": "missing", "reason": "", "ok": False}
        out.append(entry)
        if claim is not _NO_CLAIM and not (
            claim is None or (isinstance(claim, str) and MANIFEST_PATTERN.fullmatch(claim))
        ):
            entry["reason"] = "output_manifest_sha256 claim is not 64 hex digits or null"
            continue
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
            expected_output_manifest_sha256=claim if isinstance(claim, str) else None,
        )
        entry["status"], entry["reason"] = status, reason
        if status == "invalid":
            continue
        if status != "verified" and not tolerated:
            entry["reason"] = f"host run record {status} (attestation required): {reason}"
            continue
        if "output_manifest_sha256" not in run:
            entry["status"] = "invalid"
            entry["reason"] = (
                "signed run record has no output_manifest_sha256 field "
                "(absence is not a null manifest)"
            )
            continue
        if claim is None and run["output_manifest_sha256"] is not None:
            entry["status"] = "invalid"
            entry["reason"] = (
                "output_manifest_sha256 mismatch: claimed null, signed run record has "
                f"{run['output_manifest_sha256']!r}"
            )
            continue
        if claim is _NO_CLAIM and (not tolerated or lookup_missing):
            entry["reason"] = LEGACY_NO_BINDING + (
                " (run record was missing at record time)" if lookup_missing else ""
            )
            continue
        problem = run_success_problem(run)
        if problem:
            entry["reason"] = problem
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
    - A passed record that references host runs (details.host_runs, each with
      run_id and output_manifest_sha256) is rejected unless every referenced
      signed run record is present, bound to the same task/SHA and manifest
      digest, verified (or bootstrap-tolerated, see check_host_runs) and
      succeeded (succeeded true, exit_code 0, not timed out, state unchanged,
      no publication error). Legacy host_run_id(s) references carry no
      manifest binding and are rejected outside bootstrap.
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


def fill_host_run_manifests(details, task_id: str):
    """Fill output_manifest_sha256 of details.host_runs entries at record time.

    The digest comes from the signed run record under the protected host-runs
    root, so the caller cannot mis-state it. If a run record is missing the
    entry is left unbound and details.host_run_manifest_lookup = "missing" is
    set; the gate rejects such a record later.
    """
    if not isinstance(details, dict) or not isinstance(details.get("host_runs"), list):
        return details
    from infrastructure.validation import host_run_records as hrr

    for item in details["host_runs"]:
        if not isinstance(item, dict) or "output_manifest_sha256" in item:
            continue
        run_id = item.get("run_id")
        try:
            run = hrr.load_signed_run_record(run_id, task_id)
        except hrr.RunRecordMissing:
            details["host_run_manifest_lookup"] = "missing"
            continue
        item["output_manifest_sha256"] = run.get("output_manifest_sha256")
    return details


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

    details = fill_host_run_manifests(details, normalize_task_id(task_id))

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
