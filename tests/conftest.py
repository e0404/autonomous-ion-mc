"""Test configuration.

The repository has no installable Python packages yet -
`infrastructure/tasks/`, `infrastructure/validation/`, and
`infrastructure/diagnostics/` are plain script directories without
`__init__.py`. To keep that convention, each such directory that has
tests gets its own explicit `sys.path` entry here rather than relying
on package installation or implicit rootdir insertion.

The `ionmc` package itself (under `src/ionmc`) is installed normally (see
`pyproject.toml`), so tests for it are imported as `import ionmc...` with
no `sys.path` manipulation needed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DIAGNOSTICS_DIR = REPO_ROOT / "infrastructure" / "diagnostics"

if str(DIAGNOSTICS_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTICS_DIR))


def pytest_configure(config: pytest.Config) -> None:
    """Register custom markers (mirrors `[tool.pytest.ini_options].markers`).

    Explicit registration here keeps marker discovery working even when
    tests are collected without the project's `pyproject.toml` (for
    example via an ad hoc `pytest` invocation from a different rootdir).
    """
    config.addinivalue_line("markers", "warp: requires warp-lang")
    config.addinivalue_line("markers", "cuda: requires a CUDA device")
    config.addinivalue_line("markers", "network: requires network access")


@pytest.fixture(scope="session")
def warp_module():
    """Return the imported `warp` module, skipping the test if unavailable."""
    return pytest.importorskip("warp")


@pytest.fixture
def cuda_available(warp_module) -> bool:
    """Skip the test if no CUDA device is visible to `warp`."""
    if not warp_module.get_cuda_devices():
        pytest.skip("no CUDA device available")
    return True
