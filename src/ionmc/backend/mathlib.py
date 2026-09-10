"""Shared-source mechanism: one physics source, several execution bindings.

Physics functions in :mod:`ionmc.physics` are written once against a small
*math namespace* ``m`` (``m.sqrt``, ``m.log``, ``m.where`` ...) and decorated
with ``@func``. Which concrete implementation ``m`` denotes is decided when the
physics module is **loaded**:

``warp`` binding (default when ``warp`` is importable)
    ``m.<name>`` is the Warp builtin (``wp.sqrt`` ...) and ``@func`` is
    ``wp.func``. The functions can be called from Warp kernels on the CPU and
    CUDA devices. Because Warp's own ``Function.__call__`` falls back to the
    undecorated Python function when a ``@wp.func`` is called from Python
    scope (verified on Warp 1.17.0, host run ``RUN-20260910T004325Z-472e71df``),
    the very same objects are also callable from Python - but Warp's builtins
    evaluate in float32 in Python scope, so this is not the reference path.

``python`` binding
    ``m.<name>`` is the standard-library ``math`` function (float64) and
    ``@func`` is the identity. This is the **reference execution path**: the
    identical source text, executed by CPython in double precision, with no
    dependency on Warp being installed.

``numpy`` binding
    ``m.<name>`` is the numpy ufunc; the same source evaluates elementwise on
    arrays. Used for fast reference calculations and validation scripts.

The binding is selected through :func:`bind` while a module is being imported
(see :mod:`ionmc.backend.reference`), so the same file can exist several times
in ``sys.modules`` under different names with different bindings.

Only functions whose semantics agree across Warp, ``math`` and numpy are
exposed. In particular the physics code must not rely on Python's ``**`` or
``%`` on negative operands, and must use ``m.where`` instead of data-dependent
``if`` so that the numpy binding stays vectorised.
"""

from __future__ import annotations

import contextlib
import math
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from typing import Any

try:  # pragma: no cover - exercised only where warp is installed
    import warp as wp

    HAVE_WARP = True
except ImportError:  # pragma: no cover
    wp = None
    HAVE_WARP = False


def warp_module() -> Any:
    """Return the imported ``warp`` module or raise ``ImportError``."""
    if wp is None:
        raise ImportError(
            "warp-lang is not installed; only the python/numpy bindings work"
        )
    return wp


@dataclass(frozen=True)
class MathNamespace:
    """The math functions available to shared physics source.

    Attribute names follow Warp's builtin names so that the ``warp`` binding
    is a direct attribute copy.
    """

    name: str
    sqrt: Callable[..., Any]
    log: Callable[..., Any]
    log10: Callable[..., Any]
    exp: Callable[..., Any]
    pow: Callable[..., Any]
    abs: Callable[..., Any]
    min: Callable[..., Any]
    max: Callable[..., Any]
    where: Callable[..., Any]
    floor: Callable[..., Any]
    sin: Callable[..., Any]
    cos: Callable[..., Any]
    atan2: Callable[..., Any]
    acos: Callable[..., Any]
    #: Type constructor turning a Python literal into the working scalar type.
    scalar: Callable[..., Any]
    pi: float


def _python_where(cond: bool, a: Any, b: Any) -> Any:
    return a if cond else b


PYTHON_MATH = MathNamespace(
    name="python",
    sqrt=math.sqrt,
    log=math.log,
    log10=math.log10,
    exp=math.exp,
    pow=math.pow,
    abs=abs,
    min=min,
    max=max,
    where=_python_where,
    floor=math.floor,
    sin=math.sin,
    cos=math.cos,
    atan2=math.atan2,
    acos=math.acos,
    scalar=float,
    pi=math.pi,
)


def numpy_math() -> MathNamespace:
    """Math namespace bound to numpy ufuncs (imported lazily)."""
    import numpy as np

    return MathNamespace(
        name="numpy",
        sqrt=np.sqrt,
        log=np.log,
        log10=np.log10,
        exp=np.exp,
        pow=np.power,
        abs=np.abs,
        min=np.minimum,
        max=np.maximum,
        where=np.where,
        floor=np.floor,
        sin=np.sin,
        cos=np.cos,
        atan2=np.arctan2,
        acos=np.arccos,
        scalar=np.float64,
        pi=math.pi,
    )


def warp_math() -> MathNamespace:
    """Math namespace bound to Warp builtins (requires warp)."""
    w = warp_module()
    return MathNamespace(
        name="warp",
        sqrt=w.sqrt,
        log=w.log,
        log10=w.log10,
        exp=w.exp,
        pow=w.pow,
        abs=w.abs,
        min=w.min,
        max=w.max,
        where=w.where,
        floor=w.floor,
        sin=w.sin,
        cos=w.cos,
        atan2=w.atan2,
        acos=w.acos,
        scalar=float,
        pi=math.pi,
    )


# ---------------------------------------------------------------------------
# Binding selection
# ---------------------------------------------------------------------------

_DEFAULT_BINDING = "warp" if HAVE_WARP else "python"
_active_binding: list[str] = [_DEFAULT_BINDING]


def default_binding() -> str:
    """The binding used when no :func:`bind` context is active."""
    return _DEFAULT_BINDING


@contextlib.contextmanager
def bind(binding: str) -> Iterator[None]:
    """Select the math binding for physics modules imported inside the block."""
    if binding not in ("warp", "python", "numpy"):
        raise ValueError(f"unknown binding {binding!r}")
    _active_binding.append(binding)
    try:
        yield
    finally:
        _active_binding.pop()


def current() -> MathNamespace:
    """The math namespace for the currently active binding."""
    binding = _active_binding[-1]
    if binding == "warp":
        return warp_math()
    if binding == "numpy":
        return numpy_math()
    return PYTHON_MATH


def _identity_decorator(f: Callable[..., Any]) -> Callable[..., Any]:
    return f


def func(f: Callable[..., Any]) -> Any:
    """``wp.func`` under the warp binding, identity otherwise.

    Under the ``warp`` binding the returned object is a Warp ``Function`` that
    is usable from kernels and, thanks to Warp's Python-scope fallback, also
    directly callable; under the other bindings it is the plain function.
    """
    if _active_binding[-1] == "warp":
        return warp_module().func(f)
    return _identity_decorator(f)
