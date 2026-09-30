"""Protected validation snapshots from a recorded Git object, not live source."""

import os
import shutil
import stat
import subprocess
from pathlib import Path

DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW


def copy_regular_tree(source_fd, destination):
    """Copy data through anchored descriptors; never follow input symlinks."""
    destination.mkdir()
    for name in os.listdir(source_fd):
        info = os.stat(name, dir_fd=source_fd, follow_symlinks=False)
        if stat.S_ISDIR(info.st_mode):
            child = os.open(name, DIRECTORY_FLAGS, dir_fd=source_fd)
            try:
                copy_regular_tree(child, destination / name)
            finally:
                os.close(child)
        elif stat.S_ISREG(info.st_mode):
            fd = os.open(
                name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=source_fd
            )
            with os.fdopen(fd, "rb") as incoming:
                if not stat.S_ISREG(os.fstat(incoming.fileno()).st_mode):
                    raise ValueError("Snapshot input changed to a special file")
                with (destination / name).open("xb") as outgoing:
                    shutil.copyfileobj(incoming, outgoing)
        else:
            raise ValueError("Snapshot inputs must not contain links or special files")


def create(root, worktree, sha, branch, destination):
    """Fetch exactly sha into a private object database; copy only ignored inputs.

    Source edits/restoration during execution cannot affect this checkout. Git
    metadata is independent, so source ref/index updates cannot change provenance.
    Ignored .ionmc-cache data is an input snapshot, never committed source evidence.
    """
    root, worktree, destination = map(Path, (root, worktree, destination))
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
    fd = os.open(worktree, DIRECTORY_FLAGS)
    try:
        try:
            inputs = os.open(".ionmc-cache", DIRECTORY_FLAGS, dir_fd=fd)
        except FileNotFoundError:
            inputs = None
        if inputs is not None:
            try:
                copy_regular_tree(inputs, destination / ".ionmc-cache")
            finally:
                os.close(inputs)
    finally:
        os.close(fd)
    return destination
