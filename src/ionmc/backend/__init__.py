"""Execution backends: how one physics source runs on reference Python and Warp.

See :mod:`ionmc.backend.mathlib` for the shared-source mechanism and
:mod:`ionmc.backend.reference` for loading a module as a pure-Python
reference copy. Decision ``0005`` records the design.
"""

from __future__ import annotations

from ionmc.backend.mathlib import HAVE_WARP, MathNamespace, current, func, warp_module

__all__ = ["HAVE_WARP", "MathNamespace", "current", "func", "warp_module"]
