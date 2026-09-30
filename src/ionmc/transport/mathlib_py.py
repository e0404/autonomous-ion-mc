"""Math namespace with Python/NumPy float64 semantics (reference variant of the step physics).

Function names and semantics mirror the Warp builtins used by the step
physics: ``floor`` returns a float, ``randf``/``randn`` draw from a
``numpy.random.Generator`` passed as the "state".
"""

from __future__ import annotations

import math

import numpy as np

sqrt = math.sqrt
log = math.log
exp = math.exp
sin = math.sin
cos = math.cos
atan2 = math.atan2
abs = builtins_abs = abs


def floor(x: float) -> float:
    return float(math.floor(x))


def min(a: float, b: float) -> float:
    return a if a < b else b


def max(a: float, b: float) -> float:
    return a if a > b else b


def randf(state: np.random.Generator) -> float:
    return float(state.random())


def randn(state: np.random.Generator) -> float:
    return float(state.standard_normal())


def add3(array, i: int, j: int, k: int, value: float) -> float:
    """In-place accumulate into a 3-D NumPy array (mirrors ``wp.atomic_add``)."""
    old = array[i, j, k]
    array[i, j, k] = old + value
    return old
