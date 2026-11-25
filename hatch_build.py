"""Hatchling custom build hook: stamp the code identity into the wheel (decision 0045 section 3).

Active for the wheel target only. Computes the Git tree object id of ``src/ionmc`` (a copy of the
routine in ``ionmc.provenance``: the hook must not import the package under build) and ONE
read-only ``git rev-parse HEAD HEAD:src/ionmc`` with every ``GIT_*`` variable removed (status
"dirty" when the tree differs from the commit's), writes ``_build_info.json`` into a temporary
directory and maps it to ``ionmc/_build_info.json`` with ``force_include``. Nothing is
ever written under ``src/ionmc``. Without a reachable repository (for example a build from an
sdist) the stamp carries ``git_commit = null`` and ``status = "unknown"``; that is never an error.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import stat
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

_STAMP = "_build_info.json"
_HEX40 = re.compile(r"^[0-9a-f]{40}$")


def _blob_oid(data: bytes) -> bytes:
    return hashlib.sha1(b"blob %d\0" % len(data) + data).digest()


def _tree_oid(directory: Path, top: bool) -> bytes | None:
    entries: list[tuple[bytes, bytes, bytes, bytes]] = []
    for child in directory.iterdir():
        name = child.name
        if name == "__pycache__" or name.endswith(".pyc") or (top and name == _STAMP):
            continue
        raw = os.fsencode(name)
        st = child.lstat()
        if stat.S_ISLNK(st.st_mode):
            entries.append((raw, b"120000", raw, _blob_oid(os.fsencode(os.readlink(child)))))
        elif stat.S_ISDIR(st.st_mode):
            sub = _tree_oid(child, False)
            if sub is not None:
                entries.append((raw + b"/", b"40000", raw, sub))
        elif stat.S_ISREG(st.st_mode):
            mode = b"100755" if st.st_mode & stat.S_IXUSR else b"100644"
            entries.append((raw, mode, raw, _blob_oid(child.read_bytes())))
    if not entries:
        return None
    entries.sort(key=lambda e: e[0])
    body = b"".join(mode + b" " + name + b"\0" + oid for _, mode, name, oid in entries)
    return hashlib.sha1(b"tree %d\0" % len(body) + body).digest()


def _git_head(root: Path) -> tuple[str, str | None] | None:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env["LC_ALL"] = "C"
    try:
        proc = subprocess.run(  # noqa: S603 - fixed argument list, no shell
            ["git", "rev-parse", "HEAD", "HEAD:src/ionmc"],  # noqa: S607
            cwd=root, env=env, capture_output=True, text=True, timeout=5.0, check=False,
            stdin=subprocess.DEVNULL,
        )  # fmt: skip
    except (OSError, subprocess.SubprocessError):
        return None
    lines = proc.stdout.split()
    if proc.returncode != 0 or not lines or not _HEX40.match(lines[0]):
        return None
    tree = lines[1] if len(lines) > 1 and _HEX40.match(lines[1]) else None
    return lines[0], tree


class StampBuildHook(BuildHookInterface):
    PLUGIN_NAME = "custom"

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        if self.target_name != "wheel":
            return
        root = Path(self.root)
        package = root / "src" / "ionmc"
        oid = _tree_oid(package, True) if package.is_dir() else None
        head = _git_head(root)
        commit = None if head is None else head[0]
        if commit is None or oid is None:
            status = "unknown"
        else:  # a tree that differs from the commit's package tree is a dirty build
            status = "ok" if head is not None and head[1] == oid.hex() else "dirty"
        stamp = {
            "git_commit": commit,
            "package_tree_oid": None if oid is None else oid.hex(),
            "status": status,
        }
        self._tmp = tempfile.TemporaryDirectory(prefix="ionmc-stamp-")
        path = Path(self._tmp.name) / _STAMP
        path.write_text(json.dumps(stamp, sort_keys=True) + "\n", encoding="utf-8")
        build_data["force_include"][str(path)] = f"ionmc/{_STAMP}"

    def finalize(self, version: str, build_data: dict[str, Any], artifact_path: str) -> None:
        tmp = getattr(self, "_tmp", None)
        if tmp is not None:
            tmp.cleanup()
