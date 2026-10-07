"""Pure-Python stand-in for the ``warp`` names used inside the shared ``@wp.func`` factories.

The Python reference backend must not execute Warp at Python scope (decision 0037/0039: Warp
1.17 Python-scope builtin dispatch segfaulted intermittently in spawned workers). The factories
in ``funcs.py``, ``physics/em.py`` and ``physics/kinematics.py`` are therefore re-executed from
their own source text with ``wp`` replaced by this module (see ``ionmc._wpfunc.python_twin``).
Only the names those sources use are provided; arithmetic is IEEE double precision. Warp's
C-like edge cases are kept: ``log(0) = -inf``, ``log`` and ``sqrt`` of a negative number are
NaN (Python's ``math`` raises instead).
"""

from __future__ import annotations

import math
from types import SimpleNamespace
from typing import Any

_INF = math.inf
_NAN = math.nan


class Vec3:
    """Three-component double vector with the operations the shared functions use."""

    __slots__ = ("x", "y", "z")

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x, self.y, self.z = float(x), float(y), float(z)

    def __getitem__(self, i: int) -> float:
        return (self.x, self.y, self.z)[i]

    def __iter__(self) -> Any:
        return iter((self.x, self.y, self.z))

    def __len__(self) -> int:
        return 3

    def __add__(self, o: Vec3) -> Vec3:
        return Vec3(self.x + o.x, self.y + o.y, self.z + o.z)

    def __sub__(self, o: Vec3) -> Vec3:
        return Vec3(self.x - o.x, self.y - o.y, self.z - o.z)

    def __mul__(self, s: float) -> Vec3:
        return Vec3(self.x * s, self.y * s, self.z * s)

    __rmul__ = __mul__

    def __neg__(self) -> Vec3:
        return Vec3(-self.x, -self.y, -self.z)

    def __repr__(self) -> str:
        return f"Vec3({self.x!r}, {self.y!r}, {self.z!r})"


def _vector(length: int = 3, dtype: Any = float) -> type[Vec3]:
    if length != 3:
        raise ValueError("only 3-vectors are provided")
    return Vec3


def _log(x: float) -> float:
    if x > 0.0:
        return math.log(x)
    return -_INF if x == 0.0 else _NAN


def _log10(x: float) -> float:
    if x > 0.0:
        return math.log10(x)
    return -_INF if x == 0.0 else _NAN


def _sqrt(x: float) -> float:
    return math.sqrt(x) if x >= 0.0 else _NAN


def _pow(x: float, y: float) -> float:
    try:
        return math.pow(x, y)
    except (ValueError, OverflowError, ZeroDivisionError):
        return _NAN


def _exp(x: float) -> float:
    try:
        return math.exp(x)
    except OverflowError:
        return _INF


def _floor(x: float) -> float:
    """Warp ``floor`` returns a float (non-finite values pass through)."""
    return float(math.floor(x)) if math.isfinite(x) else x


def _normalize(v: Vec3) -> Vec3:
    n = math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
    return Vec3(v.x / n, v.y / n, v.z / n) if n > 0.0 else Vec3(0.0, 0.0, 0.0)


wp = SimpleNamespace(
    exp=_exp,
    log=_log,
    log10=_log10,
    sqrt=_sqrt,
    pow=_pow,
    cos=math.cos,
    sin=math.sin,
    floor=_floor,
    abs=abs,
    min=min,
    max=max,
    normalize=_normalize,
    constant=lambda x: x,
    float64=float,
    types=SimpleNamespace(vector=_vector),
)
"""Namespace substituted for ``warp`` when a factory is re-executed as a pure-Python twin."""
