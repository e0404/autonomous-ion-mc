"""Helper for precision-generic ``@wp.func`` factories.

Warp registers a function in its module by name and argument types. Two factory calls
(float32 and float64) produce functions with identical argument types in some cases (for
example a function taking only ``wp.uint32``), which would overwrite each other. The
decorator below gives every generated function a unique name suffix before registration.
"""

from __future__ import annotations

import ast
import functools
import inspect
import sys
from collections.abc import Callable
from types import SimpleNamespace
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


def _py_named_func(suffix: str) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Twin of :func:`named_func`: the plain Python function (same name suffix), no Warp."""

    def decorate(fn: Callable[..., Any]) -> Callable[..., Any]:
        fn.__name__ = f"{fn.__name__}_{suffix}"
        fn.__qualname__ = fn.__name__
        return fn

    return decorate


@functools.cache
def python_twin(factory: Callable[..., SimpleNamespace]) -> SimpleNamespace:
    """Pure-Python (float64) twin of a ``make_*(real)`` ``@wp.func`` factory.

    The source text of the factory (``inspect.getsource`` of the function under
    ``functools.cache``) is compiled again in its module's namespace with ``wp`` replaced by
    :data:`ionmc._pyshim.wp`, ``named_func`` by the plain-function decorator, ``check_real`` by
    ``float64`` and every other ``make_*`` factory by its own twin, and called with
    ``real = float``.
    The returned namespace has the same attributes as the Warp one; the functions evaluate the
    identical expressions in the identical order, with no Warp call at Python scope. The Warp
    functions stay the objects used by kernels, so one source text defines both."""
    from ionmc._pyshim import wp as shim

    fn = getattr(factory, "__wrapped__", factory)
    tree = ast.parse(inspect.getsource(fn))
    node = tree.body[0]
    assert isinstance(node, ast.FunctionDef)
    node.decorator_list = []
    module = sys.modules[fn.__module__]
    ns: dict[str, Any] = dict(vars(module))
    ns.update(wp=shim, named_func=_py_named_func, check_real=lambda real: "float64")
    ns.update(getattr(module, "PYTHON_TWIN_OVERRIDES", {}))  # python replacements of module names
    for key, val in vars(module).items():
        if key.startswith("make_") and hasattr(val, "__wrapped__") and val is not factory:
            ns[key] = functools.partial(lambda f, real: python_twin(f), val)
    exec(compile(tree, inspect.getsourcefile(fn) or "<twin>", "exec"), ns)
    twin: SimpleNamespace = ns[node.name](float)
    return twin
