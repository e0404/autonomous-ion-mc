"""Bounded, deduplicated operator alerts when Claude stops at a usage limit."""

import fcntl
import hashlib
import json
import os
import stat
import time
import uuid
from pathlib import Path

from infrastructure.experiment_v3.common import STATE, event
from infrastructure.experiment_v3.paths import EXPERIMENT_ID

SUCCESS_COOLDOWN = 15 * 60
FAILURE_COOLDOWN = 60


def limit_kind(payload):
    # Only emit fixed categories; never forward raw API messages or reset strings.
    text = str(payload.get("last_assistant_message", ""))[:4000].lower()
    for phrase, kind in (
        ("fable limit", "fable"),
        ("opus limit", "opus"),
        ("sonnet limit", "sonnet"),
        ("weekly limit", "weekly"),
        ("session limit", "session"),
    ):
        if phrase in text:
            return kind
    return "unspecified"


def handle(payload, *, sender, state=None):
    root = (
        Path(state)
        if state is not None
        else Path(
            os.environ.get("IONMC_V3_EVENT_DIR", str(STATE / "telemetry"))
        ).expanduser()
    )
    session = hashlib.sha256(
        str(payload.get("session_id") or "unknown").encode()
    ).hexdigest()[:16]
    key = hashlib.sha256(f"{EXPERIMENT_ID}:{session}".encode()).hexdigest()
    record_id = "LIMIT-" + uuid.uuid4().hex[:16]
    details = {
        "record_id": record_id,
        "session": session,
        "error": "rate_limit",
        "limit_kind": limit_kind(payload),
    }
    event("runtime_limit", details=details, state=root)
    directory = root / "notification-cooldowns"
    directory.mkdir(parents=True, exist_ok=True)
    # Anchor access to this directory and reject links to other host files.
    parent = os.open(directory, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        fd = os.open(
            key + ".json", os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600, dir_fd=parent
        )
    finally:
        os.close(parent)
    with os.fdopen(fd, "r+") as stream:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise ValueError(
                "Notification state must be a regular file with no hard-link aliases"
            )
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            status = "suppressed"
        else:
            now = time.time()
            try:
                previous = json.load(stream)
                remaining = float(previous["retry_after"]) - now
            except (ValueError, TypeError, KeyError):
                remaining = 0
            if 0 < remaining <= SUCCESS_COOLDOWN:
                status = "suppressed"
            else:
                # Reserve the retry window before network I/O, including crashes.
                def save(delay):
                    stream.seek(0)
                    json.dump({"retry_after": now + delay}, stream)
                    stream.truncate()
                    stream.flush()

                save(FAILURE_COOLDOWN)
                try:
                    result = sender(
                        record_id,
                        "runtime",
                        "usage-limit",
                        f"IonMC {EXPERIMENT_ID}: Claude stopped at a rate/usage limit. "
                        f"Session: {session}. Limit category: {details['limit_kind']}. "
                        "Inspect Claude for quota/reset details "
                        "and resume when available. "
                        "This alert does not restart Claude or switch models.",
                    )
                    status = result.get("status", "failed")
                    if status not in ("sent", "not-configured", "failed"):
                        status = "failed"
                except Exception:
                    status = "failed"
                save(SUCCESS_COOLDOWN if status == "sent" else FAILURE_COOLDOWN)
    event(
        "notification_result",
        details={"record_id": record_id, "session": session, "status": status},
        state=root,
    )
    return status
