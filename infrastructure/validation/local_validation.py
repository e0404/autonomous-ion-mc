#!/usr/bin/env python3

import argparse
import json
import os
import re
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path


WORKTREE_ROOT = Path.home() / "aiprojects" / "ion-mc-worktrees"

VALIDATION_ROOT = Path(
    os.environ.get(
        "IONMC_VALIDATION_DIR",
        str(Path.home() / ".local/share/ionmc-experiment/validation"),
    )
)


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


def records_for(task_id: str, sha: str) -> list[dict]:
    directory = VALIDATION_ROOT / normalize_task_id(task_id) / sha

    if not directory.exists():
        return []

    records = []

    for path in sorted(directory.glob("*.json")):
        try:
            records.append(json.loads(path.read_text(encoding="utf-8")))
        except Exception:
            continue

    return records


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

    passing = [r for r in records if r.get("status") == "passed"]
    failing = [r for r in records if r.get("status") == "failed"]

    result = {
        "status": "present" if records else "missing",
        "task_id": normalize_task_id(task_id),
        "git_branch": state["branch"],
        "git_sha": state["sha"],
        "worktree_dirty": state["dirty"],
        "records": records,
        "passing_records": len(passing),
        "failing_records": len(failing),
        "gate_passed": (
            not state["dirty"]
            and len(passing) > 0
            and len(failing) == 0
        ),
    }

    print(json.dumps(result, indent=2))


def check(task_id: str):
    state = worktree_state(task_id)
    records = records_for(task_id, state["sha"])

    passing = [r for r in records if r.get("status") == "passed"]
    failing = [r for r in records if r.get("status") == "failed"]

    gate_passed = (
        not state["dirty"]
        and len(passing) > 0
        and len(failing) == 0
    )

    result = {
        "status": "passed" if gate_passed else "failed",
        "task_id": normalize_task_id(task_id),
        "git_sha": state["sha"],
        "worktree_dirty": state["dirty"],
        "passing_records": len(passing),
        "failing_records": len(failing),
    }

    print(json.dumps(result, indent=2))

    if not gate_passed:
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
