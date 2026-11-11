"""Read-only archive access and byte-preserving host-run record transfer.

Only record.json is exported. Signatures are neither created nor interpreted;
consumers must verify the complete record against their committed trust root.
Directory descriptors keep writes confined even if a task renames cache paths.
"""

import hashlib
import json
import os
import re
import stat
import uuid
from contextlib import contextmanager
from pathlib import Path

from infrastructure.experiment_v3.paths import ROOT, SHARE_ROOT
from infrastructure.experiment_v3.worker_boundary import validate_task

RUN_PATTERN = re.compile(r"RUN-\d{8}T\d{6}Z-[0-9a-f]{8}")
MAX_RECORD_BYTES = 4 * 1024 * 1024
MAX_BATCH_BYTES = 64 * 1024 * 1024
DIR_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
FILE_FLAGS = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC


def _task(task_id):
    if not isinstance(task_id, str):
        raise ValueError("Expected a V3 task ID")
    value = task_id.strip().upper()
    if not re.fullmatch(r"V3-[A-Z0-9][A-Z0-9._-]{0,60}", value):
        raise ValueError("Expected a safe V3 task ID")
    return value


def _run(run_id):
    if not isinstance(run_id, str) or not RUN_PATTERN.fullmatch(run_id):
        raise ValueError("Expected a host RUN-<timestamp>-<suffix> identifier")
    return run_id


@contextmanager
def _directory(name, parent=None, *, create=False):
    if create:
        try:
            os.mkdir(name, mode=0o700, dir_fd=parent)
        except FileExistsError:
            pass
    fd = os.open(name, DIR_FLAGS, dir_fd=parent)
    try:
        yield fd
    finally:
        os.close(fd)


@contextmanager
def _archive(task_id, state):
    with (
        _directory(Path(state)) as state_fd,
        _directory("host-runs", state_fd) as runs_fd,
        _directory(task_id, runs_fd) as task_fd,
    ):
        yield task_fd


def _read(parent, name):
    fd = os.open(name, FILE_FLAGS, dir_fd=parent)
    with os.fdopen(fd, "rb") as stream:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode) or info.st_size > MAX_RECORD_BYTES:
            raise ValueError("Record must be a regular file of at most 4 MiB")
        data = stream.read(MAX_RECORD_BYTES + 1)
    if len(data) > MAX_RECORD_BYTES:
        raise ValueError("Record exceeds 4 MiB")
    return data


def _record(task_fd, task_id, run_id):
    with _directory(run_id, task_fd) as run_fd:
        data = _read(run_fd, "record.json")
    record = json.loads(data)
    if (
        not isinstance(record, dict)
        or record.get("task_id") != task_id
        or record.get("run_id") != run_id
    ):
        raise ValueError("Archived record task/run identity mismatch")
    return data, record


def _metadata(run_id, data, record):
    sha = record.get("sha")
    if sha is not None and (
        not isinstance(sha, str) or not re.fullmatch(r"[0-9a-f]{40}", sha)
    ):
        raise ValueError("Malformed archived code SHA")
    return {
        "run_id": run_id,
        "record_available": True,
        "sha": sha,
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        # Presence is deliberately not a cryptographic verification result.
        "attestation_present": isinstance(record.get("attestation"), dict),
    }


def list_records(task_id, offset=0, limit=50, *, state=SHARE_ROOT):
    """Paginate archived runs, including runs without a persisted record.json.

    Stable timestamp order; new runs append in normal operation. No stdout,
    stderr, argv, private keys or execution snapshots are returned.
    """
    task_id = _task(task_id)
    if (
        type(offset) is not int
        or type(limit) is not int
        or offset < 0
        or not 1 <= limit <= 100
    ):
        raise ValueError("Use offset >= 0 and limit 1..100")
    with _archive(task_id, state) as task_fd:
        names = sorted(n for n in os.listdir(task_fd) if RUN_PATTERN.fullmatch(n))
        records = []
        for run_id in names[offset : offset + limit]:
            try:
                data, record = _record(task_fd, task_id, run_id)
            except FileNotFoundError:
                records.append({"run_id": run_id, "record_available": False})
            else:
                records.append(_metadata(run_id, data, record))
    return {
        "task_id": task_id,
        "total_runs": len(names),
        "next_offset": offset + limit if offset + limit < len(names) else None,
        "records": records,
    }


def _publish(destination_fd, name, data):
    """Atomic publication without overwriting, following links or changing modes."""
    try:
        existing = _read(destination_fd, name)
    except FileNotFoundError:
        pass
    else:
        if existing != data:
            raise ValueError(f"Conflicting destination record: {name}")
        return "existing"
    staging = ".record-" + uuid.uuid4().hex
    fd = os.open(
        staging,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC,
        0o600,
        dir_fd=destination_fd,
    )
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(
                staging,
                name,
                src_dir_fd=destination_fd,
                dst_dir_fd=destination_fd,
                follow_symlinks=False,
            )
        except FileExistsError:
            if _read(destination_fd, name) != data:
                raise ValueError(f"Conflicting destination record: {name}") from None
            status = "existing"
        else:
            status = "copied"
        if _read(destination_fd, name) != data:
            raise ValueError(f"Record changed during transfer: {name}")
        return status
    finally:
        os.unlink(staging, dir_fd=destination_fd)


def materialize_records(task_id, run_ids, *, state=SHARE_ROOT, root=ROOT):
    """Copy 1..100 records (max 64 MiB) into the same registered task's cache.

    Dirty worktrees are allowed: this retrieves existing evidence, not execution
    or validation. Each file is published atomically; retrying a partially
    completed batch is safe. Unsigned records remain unsigned, without fallback
    to result.json or access to other archive files.
    """
    task_id = _task(task_id)
    if not isinstance(run_ids, list) or not 1 <= len(run_ids) <= 100:
        raise ValueError("Select 1..100 unique run IDs")
    selected = [_run(run_id) for run_id in run_ids]
    if len(set(selected)) != len(selected):
        raise ValueError("Select unique run IDs")
    root = Path(root).resolve()
    worktree = root.parent / (root.name + "-worktrees") / task_id.lower()
    _, worktree = validate_task(task_id, worktree, root=root)
    batch = []
    size = 0
    with _archive(task_id, state) as task_fd:
        for run_id in selected:
            data, record = _record(task_fd, task_id, run_id)
            size += len(data)
            if size > MAX_BATCH_BYTES:
                raise ValueError("Batch exceeds 64 MiB; select fewer run IDs")
            batch.append((run_id, data, _metadata(run_id, data, record)))
    destination = worktree / ".ionmc-cache" / "host-run-records"
    files = []
    with (
        _directory(worktree) as work_fd,
        _directory(".ionmc-cache", work_fd, create=True) as cache_fd,
        _directory("host-run-records", cache_fd, create=True) as destination_fd,
    ):
        for run_id, data, metadata in batch:
            name = run_id + ".json"
            status = _publish(destination_fd, name, data)
            files.append(
                {**metadata, "path": str(destination / name), "status": status}
            )
    return {"task_id": task_id, "files": files, "bytes": size}
