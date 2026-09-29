"""Stamp the exact source identity into built wheels.

Persisted scientific results must recover the code SHA and dirty state outside
the producing process (V2-OUTPUT). When ``ionmc`` runs from a Git checkout the
identity is read live; an installed wheel instead carries this build-time
record. A wheel built from a dirty tree is marked dirty rather than clean.
"""

from __future__ import annotations

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


def git_identity(root: Path) -> dict:
    def run(*args: str) -> str:
        return subprocess.check_output(
            ["git", "-C", str(root), *args], text=True, stderr=subprocess.DEVNULL
        ).strip()

    try:
        sha = run("rev-parse", "HEAD")
        dirty = bool(run("status", "--porcelain"))
        return {"git_sha": sha, "git_dirty": dirty, "source": "build-time-git"}
    except (OSError, subprocess.CalledProcessError):
        return {"git_sha": None, "git_dirty": None, "source": "unknown"}


class BuildInfoHook(BuildHookInterface):
    PLUGIN_NAME = "custom"

    def initialize(self, version: str, build_data: dict) -> None:
        root = Path(self.root)
        target = root / "src" / "ionmc" / "_build_info.json"
        info = git_identity(root)
        if info["source"] == "unknown" and target.is_file():
            # Building a wheel from an sdist: keep the identity stamped when the
            # sdist itself was built from Git instead of reporting "unknown".
            previous = json.loads(target.read_text())
            if previous.get("git_sha"):
                info = {k: previous.get(k) for k in ("git_sha", "git_dirty", "source")}
                info["source"] = "sdist:" + str(previous["source"])
        info |= {
            "version": self.metadata.version,
            "built_at_utc": datetime.now(UTC).isoformat(),
        }
        target.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")
        archive_path = (
            "ionmc/_build_info.json"
            if self.target_name == "wheel"
            else "src/ionmc/_build_info.json"
        )
        build_data.setdefault("force_include", {})[str(target)] = archive_path
