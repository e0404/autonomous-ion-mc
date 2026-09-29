"""Load the shared step physics once per math namespace.

``physics("warp")`` returns the module whose ``MATH`` is Warp's builtins (for
kernels); ``physics("python")`` returns the same source bound to Python's
``math`` (for the float64 reference backend). Both come from
``step_physics.py``; nothing is duplicated.
"""

from __future__ import annotations

import importlib.util
import sys
from functools import cache
from pathlib import Path
from types import ModuleType

_SOURCE = Path(__file__).with_name("step_physics.py")


@cache
def physics(variant: str) -> ModuleType:
    mathlib: ModuleType
    if variant == "warp":
        from ionmc.transport import mathlib_warp

        mathlib = mathlib_warp
    elif variant == "python":
        from ionmc.transport import mathlib_py

        mathlib = mathlib_py
    else:
        raise ValueError(f"unknown physics variant {variant!r}")
    name = f"ionmc.transport._step_physics_{variant}"
    spec = importlib.util.spec_from_file_location(name, _SOURCE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    module.MATH = mathlib  # type: ignore[attr-defined]
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
