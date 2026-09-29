"""Exact source identity of the running ``ionmc`` code.

Two situations are distinguished explicitly rather than guessed:

* running from a Git checkout (development, task worktrees, host-runner
  validation): the commit SHA and dirty state are read live from Git;
* running from an installed wheel: the identity stamped at build time by
  ``hatch_build.py`` is reported, including whether the build tree was dirty.

If neither is available the result says so (``source == "unknown"``) instead of
pretending a clean known SHA.
"""

from __future__ import annotations

import json
import subprocess
from importlib import resources
from pathlib import Path
from typing import TypedDict


class CodeIdentity(TypedDict):
    version: str
    git_sha: str | None
    git_dirty: bool | None
    source: str


def _repository_root(start: Path) -> Path | None:
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    return None


def _git(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True, stderr=subprocess.DEVNULL
    ).strip()


def _from_git(root: Path) -> CodeIdentity | None:
    try:
        sha = _git(root, "rev-parse", "HEAD")
        dirty = bool(_git(root, "status", "--porcelain"))
    except (OSError, subprocess.CalledProcessError):
        return None
    from ionmc import __version__

    return {"version": __version__, "git_sha": sha, "git_dirty": dirty, "source": "git"}


def _from_build_info() -> CodeIdentity | None:
    try:
        text = resources.files("ionmc").joinpath("_build_info.json").read_text()
    except (FileNotFoundError, OSError):
        return None
    info = json.loads(text)
    return {
        "version": str(info.get("version")),
        "git_sha": info.get("git_sha"),
        "git_dirty": info.get("git_dirty"),
        "source": "wheel:" + str(info.get("source")),
    }


def code_identity() -> CodeIdentity:
    """Return the exact code identity, preferring a live Git checkout."""
    root = _repository_root(Path(__file__).resolve().parent)
    if root is not None and (identity := _from_git(root)) is not None:
        return identity
    if (identity := _from_build_info()) is not None:
        return identity
    from ionmc import __version__

    return {
        "version": __version__,
        "git_sha": None,
        "git_dirty": None,
        "source": "unknown",
    }
