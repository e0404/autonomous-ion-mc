"""Hash-verified reference output transfer without bulk data in model context."""

import json
import os
import shutil
import tempfile
from pathlib import Path

from infrastructure.experiment_v3.common import (
    STATE,
    file_hash,
    identifier,
    safe_path,
    write_json,
)


def manifest(run_id, *, state=STATE):
    archive = Path(state) / "references" / identifier(run_id)
    record = json.loads((archive / "result.json").read_text())
    return archive, record, record["artifact_files"]


def archived_file(archive, name):
    relative = Path(name)
    if relative.is_absolute() or ".." in relative.parts or not relative.parts:
        raise ValueError("Invalid archive path")
    # Some engines emit internal aliases (FRED Dose.mhd). Resolve only within
    # this immutable archive, then copy bytes rather than exporting symlinks.
    source = (archive / relative).resolve(strict=True)
    if not source.is_relative_to(archive.resolve()) or not source.is_file():
        raise ValueError("Artifact escapes archive or is not a regular file")
    return source


def list_artifacts(run_id, *, offset=0, limit=50, state=STATE):
    if (
        type(offset) is not int
        or type(limit) is not int
        or offset < 0
        or not 1 <= limit <= 100
    ):
        raise ValueError("Use offset >= 0 and limit 1..100")
    _, record, files = manifest(run_id, state=state)
    names = sorted(files)
    return {
        "run_id": run_id,
        "sha": record["code_sha"],
        "total_files": len(names),
        "next_offset": offset + limit if offset + limit < len(names) else None,
        "files": [
            {"path": name, **files[name]} for name in names[offset : offset + limit]
        ],
    }


def materialize(run_id, worktree, paths=None, *, state=STATE):
    archive, record, files = manifest(run_id, state=state)
    selected = sorted(files) if paths is None else paths
    if not selected or len(selected) > 10000 or len(set(selected)) != len(selected):
        raise ValueError("Select 1..10000 unique archived paths")
    # The fixed destination prevents arbitrary writes. Parent symlinks are rejected.
    relative = ".ionmc-cache/reference-runs/" + identifier(run_id)
    destination = safe_path(worktree, relative)
    for parent in [destination, *destination.parents]:
        if parent == Path(worktree).resolve():
            break
        if parent.is_symlink():
            raise ValueError("Symlink destination is forbidden")
    source_files = {}
    for name in selected:
        if name not in files:
            raise ValueError("Unknown artifact path; use list_reference_artifacts")
        source = archived_file(archive, name)
        if not source.is_file() or file_hash(source) != files[name]["sha256"]:
            raise ValueError("Archived artifact hash mismatch")
        source_files[name] = source
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".transfer-", dir=destination.parent))
    try:
        for name, source in source_files.items():
            target = safe_path(staging, name)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            if file_hash(target) != files[name]["sha256"]:
                raise ValueError("Artifact changed during transfer")
        provenance = {
            "run_id": run_id,
            "source_sha": record["code_sha"],
            "engine": record["engine"],
            "files": {n: files[n] for n in selected},
        }
        write_json(staging / "transfer-manifest.json", provenance)
        if destination.exists():
            existing = destination / "transfer-manifest.json"
            if (
                not existing.is_file()
                or json.loads(existing.read_text()) != provenance
                or any(
                    file_hash(safe_path(destination, n)) != files[n]["sha256"]
                    for n in selected
                )
            ):
                raise ValueError(
                    "Destination exists with different content; no files overwritten"
                )
        else:
            os.rename(staging, destination)
        return {
            "run_id": run_id,
            "path": str(destination),
            "manifest": str(destination / "transfer-manifest.json"),
            "files": len(selected),
            "bytes": sum(files[n]["bytes"] for n in selected),
            "source_sha": record["code_sha"],
        }
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def read_text(run_id, path, *, offset=0, limit=2000, state=STATE):
    if (
        type(offset) is not int
        or type(limit) is not int
        or offset < 0
        or not 1 <= limit <= 4000
    ):
        raise ValueError(
            "Text previews are limited to 4000 bytes; materialize files for analysis"
        )
    archive, _, files = manifest(run_id, state=state)
    if path not in files:
        raise ValueError("Unknown artifact path; use list_reference_artifacts")
    source = archived_file(archive, path)
    if file_hash(source) != files[path]["sha256"]:
        raise ValueError("Archived artifact hash mismatch")
    with source.open("rb") as stream:
        stream.seek(offset)
        data = stream.read(limit)
    try:
        text = data.decode("utf-8")
        if "\x00" in text:
            raise UnicodeError()
    except UnicodeError as exc:
        raise ValueError(
            "Binary artifact; use materialize_reference_artifacts"
        ) from exc
    return {
        "text": text,
        "offset": offset,
        "bytes": len(data),
        "total_bytes": source.stat().st_size,
        "truncated": offset + len(data) < source.stat().st_size,
        "sha256": files[path]["sha256"],
    }
