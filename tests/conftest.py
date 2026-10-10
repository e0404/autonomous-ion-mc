"""Test configuration.

The repository has no installable Python packages yet -
`infrastructure/tasks/`, `infrastructure/validation/`, and
`infrastructure/diagnostics/` are plain script directories without
`__init__.py`. To keep that convention, each such directory that has
tests gets its own explicit `sys.path` entry here rather than relying
on package installation or implicit rootdir insertion.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DIAGNOSTICS_DIR = REPO_ROOT / "infrastructure" / "diagnostics"

if str(DIAGNOSTICS_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTICS_DIR))


import pytest  # noqa: E402


@pytest.fixture(autouse=True)
def _clear_nuclear_device_cache() -> object:
    """Drop the per-process nuclear device cache around every test.

    The packed nuclear device is cached per process (``ionmc.transport.nuclear_device``). Tests
    that monkeypatch ``pack_nuclear``, the table rows or kernel limits must not see a device built
    by an earlier test, nor leave theirs behind (only if the module is already imported)."""
    nd = sys.modules.get("ionmc.transport.nuclear_device")
    if nd is not None:
        nd.clear_nuclear_device_cache()
    yield
    nd = sys.modules.get("ionmc.transport.nuclear_device")
    if nd is not None:
        nd.clear_nuclear_device_cache()
