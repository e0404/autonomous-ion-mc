"""Helper for precision-generic ``@wp.func`` factories.

Warp registers a function in its module by name and argument types. Two factory calls
(float32 and float64) produce functions with identical argument types in some cases (for
example a function taking only ``wp.uint32``), which would overwrite each other. The
decorator below gives every generated function a unique name suffix before registration.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import warp as wp

SUPPORTED_REALS = ("float32", "float64")


def check_real(real: Any) -> str:
    """Return ``'float32'`` or ``'float64'`` for a Warp real type, else raise ValueError."""
    name = getattr(real, "__name__", None)
    if real is wp.float32 or name == "float32":
        return "float32"
    if real is wp.float64 or name == "float64":
        return "float64"
    raise ValueError(f"precision must be wp.float32 or wp.float64, got {real!r}")


def named_func(suffix: str) -> Callable[[Callable[..., Any]], Any]:
    """Decorator factory: register a ``@wp.func`` named ``<name>_<suffix>``."""

    def decorate(fn: Callable[..., Any]) -> Any:
        fn.__name__ = f"{fn.__name__}_{suffix}"
        fn.__qualname__ = fn.__name__
        return wp.func(fn)

    return decorate
