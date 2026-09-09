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
