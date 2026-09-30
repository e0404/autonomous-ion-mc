"""Require the upstream bubblewrap fix for CVE-2026-87766."""

import re
import shutil
import subprocess

from infrastructure.experiment_v3.common import STATE


def inspect():
    dedicated = STATE / "runtime/bin/bwrap"
    executable = str(dedicated) if dedicated.is_file() else shutil.which("bwrap")
    if not executable:
        return {"ready": False, "version": None, "executable": None}
    result = subprocess.run(
        [executable, "--version"],
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )
    match = re.search(r"bubblewrap (\d+)\.(\d+)\.(\d+)", result.stdout)
    version = tuple(map(int, match.groups())) if match else ()
    return {
        "ready": result.returncode == 0 and version >= (0, 12, 0),
        "version": ".".join(map(str, version)),
        "executable": executable,
    }


def require():
    status = inspect()
    if not status["ready"]:
        raise RuntimeError("bubblewrap >= 0.12.0 is required (CVE-2026-87766)")
    return status["executable"]
