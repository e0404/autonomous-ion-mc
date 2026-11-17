"""Protected validation snapshots from a recorded Git object, not live source."""

import json
import os
import subprocess
from pathlib import Path

from infrastructure.experiment_v3.host_inputs import (
    DIRECTORY_FLAGS,
    InputStore,
    sync_dir,
)
from infrastructure.experiment_v3.storage_guard import StorageGuard


def create(
    root,
    worktree,
    sha,
    branch,
    destination,
    *,
    include_inputs=True,
    input_store=None,
    guard=None,
):
    """Fetch exactly sha; capture ignored inputs through a protected shared store.

    Source edits/restoration during execution cannot affect this checkout. Git
    metadata is independent, so source ref/index updates cannot change provenance.
    Ignored .ionmc-cache data is an input snapshot, never committed source evidence.
    """
    root, worktree, destination = map(Path, (root, worktree, destination))
    if include_inputs:
        guard = guard or StorageGuard(destination.parent)
        guard.require()
    env = dict(os.environ, GIT_CONFIG_NOSYSTEM="1", GIT_CONFIG_GLOBAL="/dev/null")

    def run(*args):
        return subprocess.check_output(
            ["git", "-c", "core.hooksPath=/dev/null", *args],
            env=env,
            stderr=subprocess.PIPE,
            text=True,
        ).strip()

    run("-c", "init.templateDir=", "init", str(destination))
    run("-C", str(destination), "fetch", "--no-tags", str(root / ".git"), sha)
    run("-C", str(destination), "checkout", "-B", branch, sha)
    if run("-C", str(destination), "rev-parse", "HEAD") != sha:
        raise ValueError("Snapshot SHA mismatch")
    # Git symlinks are unnecessary for this runtime and complicate nested mounts.
    if any(p.is_symlink() for p in destination.rglob("*")):
        raise ValueError("Validation snapshots must not contain committed symlinks")
    if not include_inputs:
        return destination
    entries = []
    fd = os.open(worktree, DIRECTORY_FLAGS)
    try:
        try:
            inputs = os.open(".ionmc-cache", DIRECTORY_FLAGS, dir_fd=fd)
        except FileNotFoundError:
            inputs = None
        if inputs is not None:
            try:
                with InputStore(
                    input_store or destination.parent / "input-objects", guard
                ) as store:
                    entries = store.capture_tree(inputs, destination / ".ionmc-cache")
            finally:
                os.close(inputs)
    finally:
        os.close(fd)
    manifest = {"schema_version": 1, "sha": sha, "entries": entries}
    manifest_path = destination.parent / "input-manifest.json"
    with manifest_path.open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, sort_keys=True, indent=2)
        stream.write("\n")
        stream.flush()
        os.fchmod(stream.fileno(), 0o444)
        os.fsync(stream.fileno())
    sync_dir(destination.parent)
    return destination
