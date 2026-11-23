"""Runtime environment description for provenance and diagnostics."""

from __future__ import annotations

import platform
from typing import Any

import numpy

from ionmc._version import __version__
from ionmc.provenance import code_identity


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

    Code identity comes from :func:`ionmc.provenance.code_identity` (decision 0045 section 3): the
    build stamp ``ionmc/_build_info.json`` of an installed wheel, else exactly one read-only,
    timeout-bounded ``git`` query with every inherited ``GIT_*`` variable removed, else unknown.
    That query is the only subprocess the library runs; it never uses the network.
    """
    code = code_identity()
    warp_info = _describe_warp()
    return {
        "ionmc_version": __version__,
        "python_version": platform.python_version(),
        "numpy_version": numpy.__version__,
        "warp_version": warp_info["version"],
        "warp_error": warp_info["error"],
        "warp_devices": warp_info["devices"],
        "git_sha": code["git_commit"],
        "git_dirty": code["git_dirty"],
    }
