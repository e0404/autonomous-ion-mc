"""Small standard-library primitives shared by the v3 infrastructure."""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from infrastructure.experiment_v3.paths import ROOT as ROOT
from infrastructure.experiment_v3.paths import SHARE_ROOT

STATE = Path(os.environ.get("IONMC_V3_STATE", str(SHARE_ROOT))).expanduser()


def now():
    return datetime.now(UTC).isoformat()


def digest(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def file_hash(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def exact_state(root):
    sha = git(root, "rev-parse", "HEAD")
    if git(root, "status", "--porcelain", "--untracked-files=all"):
        raise ValueError("A clean committed checkout is required")
    return sha


def safe_path(root, relative):
    root = Path(root).resolve()
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError("Expected a nonempty relative path without parent traversal")
    result = root / path
    if result.is_symlink() or not result.resolve().is_relative_to(root):
        raise ValueError("Path escapes its allowed root")
    return result


def inventory(root, *, allow_internal_links=False):
    root = Path(root)
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink() and not (
            allow_internal_links
            and path.resolve().is_relative_to(root.resolve())
            and path.resolve().is_file()
        ):
            raise ValueError(
                f"Symlinks are not allowed in case/output bundles: {path.name}"
            )
        if path.is_file():
            files[path.relative_to(root).as_posix()] = {
                "sha256": file_hash(path),
                "bytes": path.stat().st_size,
            }
        elif not path.is_dir():
            raise ValueError("Special files are not allowed")
    return files


EVENTS = {
    "reference_started",
    "reference_finished",
    "reference_failed",
    "data_acquired",
    "data_role_assigned",
    "validation_strategy_changed",
    "internal_evidence_rejected",
    "scientific_failure",
    "reference_contradiction",
    "model_changed",
    "tolerance_changed",
    "performance_failure",
    "workflow_failure",
    "permission_prompt",
    "notification_result",
    "release_evaluated",
    "release_promoted",
    "intervention_requested",
    "intervention_resolved",
    "targets_frozen",
    "agent_routed",
}


def event(name, *, sha=None, task_id=None, details=None, state=None):
    if name not in EVENTS:
        raise ValueError(f"Unknown v3 event: {name}")
    record = {
        "schema_version": 2,
        "experiment_id": "experiment-v3",
        "at": now(),
        "event": name,
        "sha": sha,
        "task_id": task_id,
        "details": details or {},
    }
    # Pass scientific record IDs, never environments, prompts or credentials.
    event_root = (
        Path(state)
        if state is not None
        else Path(
            os.environ.get("IONMC_V3_EVENT_DIR", str(STATE / "telemetry"))
        ).expanduser()
    )
    path = event_root / "events.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o600)
    try:
        os.write(fd, (json.dumps(record, allow_nan=False) + "\n").encode())
    finally:
        os.close(fd)
    return record


def identifier(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", value):
        raise ValueError("Invalid record identifier")
    return value
