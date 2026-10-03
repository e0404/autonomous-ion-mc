"""Small validation helpers shared by the configuration dataclasses (fail closed)."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from ionmc.errors import UnsupportedCombinationError


def fail(message: str) -> UnsupportedCombinationError:
    """Return the engine's configuration error for ``message``."""
    return UnsupportedCombinationError(message)


def real(name: str, value: Any, *, positive: bool = False, nonnegative: bool = False) -> float:
    """Return ``value`` as a finite float, optionally > 0 or >= 0; ``bool`` is rejected."""
    if isinstance(value, bool) or not isinstance(value, int | float | np.integer | np.floating):
        raise fail(f"{name} must be a real number, got {value!r}")
    v = float(value)
    if not math.isfinite(v):
        raise fail(f"{name} must be finite, got {value!r}")
    if positive and not v > 0.0:
        raise fail(f"{name} must be > 0, got {value!r}")
    if nonnegative and v < 0.0:
        raise fail(f"{name} must be >= 0, got {value!r}")
    return v


def integer(name: str, value: Any, *, minimum: int | None = None) -> int:
    """Return ``value`` as an int (``bool`` and non-integral floats are rejected)."""
    if isinstance(value, bool) or not isinstance(value, int | np.integer):
        raise fail(f"{name} must be an integer, got {value!r}")
    v = int(value)
    if minimum is not None and v < minimum:
        raise fail(f"{name} must be >= {minimum}, got {value!r}")
    return v


def triple(name: str, value: Any, *, positive: bool = False) -> tuple[float, float, float]:
    """Return three finite floats (optionally all > 0)."""
    if not isinstance(value, tuple | list | np.ndarray) or len(value) != 3:
        raise fail(f"{name} must be a sequence of three numbers, got {value!r}")
    return (
        real(f"{name}[0]", value[0], positive=positive),
        real(f"{name}[1]", value[1], positive=positive),
        real(f"{name}[2]", value[2], positive=positive),
    )


def int_triple(name: str, value: Any) -> tuple[int, int, int]:
    """Return three positive integers."""
    if not isinstance(value, tuple | list | np.ndarray) or len(value) != 3:
        raise fail(f"{name} must be a sequence of three integers, got {value!r}")
    return (
        integer(f"{name}[0]", value[0], minimum=1),
        integer(f"{name}[1]", value[1], minimum=1),
        integer(f"{name}[2]", value[2], minimum=1),
    )


def choice(name: str, value: Any, allowed: Sequence[str]) -> str:
    """Return ``value`` if it is one of ``allowed``, else raise."""
    if value not in allowed:
        raise fail(f"{name} must be one of {tuple(allowed)}, got {value!r}")
    return str(value)
