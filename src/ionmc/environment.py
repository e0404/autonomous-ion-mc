"""Runtime environment description for provenance and diagnostics."""

from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy

from ionmc._version import __version__

_BUILD_INFO = Path(__file__).with_name("_build_info.json")


def _read_build_info() -> dict[str, Any]:
    """Return the tracked build info (git SHA and dirty flag) if present, else empty."""
    try:
        data = json.loads(_BUILD_INFO.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _describe_warp() -> dict[str, Any]:
    """Describe Warp and its devices; import and init are lazy and failure-tolerant."""
    try:
        import warp

        warp.config.log_level = warp.LOG_WARNING
        warp.init()
        devices: list[dict[str, Any]] = []
        for device in warp.get_devices():
            devices.append(
                {
                    "alias": device.alias,
                    "is_cuda": bool(device.is_cuda),
                    "name": getattr(device, "name", None),
                    "arch": getattr(device, "arch", None) if device.is_cuda else None,
                }
            )
        return {"version": getattr(warp, "__version__", None), "devices": devices, "error": None}
    except Exception as exc:  # noqa: BLE001 - diagnostics must never raise
        return {"version": None, "devices": [], "error": f"{type(exc).__name__}: {exc}"}


def describe_environment() -> dict[str, Any]:
    """Describe the software environment as a JSON-serialisable dictionary.

    Git information is read from the tracked ``_build_info.json`` file if it exists;
    git is never executed.
    """
    build = _read_build_info()
    warp_info = _describe_warp()
    return {
        "ionmc_version": __version__,
        "python_version": platform.python_version(),
        "numpy_version": numpy.__version__,
        "warp_version": warp_info["version"],
        "warp_error": warp_info["error"],
        "warp_devices": warp_info["devices"],
        "git_sha": build.get("git_sha"),
        "git_dirty": build.get("git_dirty"),
    }
