#!/usr/bin/env python3

import json
import os
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


EXPERIMENT_ID = os.environ.get("IONMC_EXPERIMENT_ID", "unknown")

TELEMETRY_DIR = Path(
    os.environ.get(
        "IONMC_TELEMETRY_DIR",
        str(Path.home() / ".local/share/ionmc-experiment/telemetry"),
    )
)

OUTPUT = TELEMETRY_DIR / "experiment-events.jsonl"


def git(*args):
    try:
        return subprocess.check_output(
            ["git", *args],
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=2,
        ).strip()
    except Exception:
        return None


def git_metadata():
    root = git("rev-parse", "--show-toplevel")

    if root is None:
        return {
            "git_root": None,
            "git_branch": None,
            "git_sha": None,
            "git_dirty": None,
        }

    status = git("status", "--porcelain")

    return {
        "git_root": root,
        "git_branch": git("branch", "--show-current"),
        "git_sha": git("rev-parse", "HEAD"),
        "git_dirty": bool(status) if status is not None else None,
    }


def selected_payload(data):
    event = data.get("hook_event_name")

    fields_by_event = {
        "SessionStart": [
            "source",
            "model",
            "agent_type",
            "session_title",
        ],
        "SessionEnd": [
            "reason",
        ],
        "TaskCreated": [
            "task_id",
            "task_subject",
            "teammate_name",
        ],
        "TaskCompleted": [
            "task_id",
            "task_subject",
            "teammate_name",
        ],
        "SubagentStart": [
            "agent_id",
            "agent_type",
        ],
        "SubagentStop": [
            "agent_id",
            "agent_type",
        ],
        "Notification": [
            "notification_type",
            "title",
        ],
        "Stop": [
            "stop_hook_active",
        ],
        "StopFailure": [
            "error",
            "error_details",
        ],
        "WorktreeCreate": [
            "name",
        ],
        "WorktreeRemove": [
            "worktree_path",
        ],
        "ConfigChange": [
            "source",
            "file_path",
        ],
        "PreModelSwitch": [
            "from_model",
            "to_model",
        ],
        "PostModelSwitch": [
            "from_model",
            "to_model",
        ],
    }

    return {
        key: data.get(key)
        for key in fields_by_event.get(event, [])
        if key in data
    }


def main():
    try:
        data = json.load(sys.stdin)
    except Exception as exc:
        print(f"ionmc telemetry hook: invalid input: {exc}", file=sys.stderr)
        return 0

    TELEMETRY_DIR.mkdir(parents=True, exist_ok=True)

    record = {
        "schema_version": 1,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "experiment_id": EXPERIMENT_ID,
        "event": data.get("hook_event_name"),
        "session_id": data.get("session_id"),
        "permission_mode": data.get("permission_mode"),
        "cwd": data.get("cwd"),
        "hostname": socket.gethostname(),
        **git_metadata(),
        "data": selected_payload(data),
    }

    with OUTPUT.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, separators=(",", ":"), ensure_ascii=False))
        f.write("\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
