"""Load a shared-source physics module under a non-Warp binding.

:func:`load_bound_module` executes the *same source file* as an ``ionmc.physics``
module a second time, under a different name and with the requested math
binding active, and returns the resulting module object. The module's
functions are then plain Python callables computing in float64 (``python``
binding) or numpy ufunc compositions (``numpy`` binding). This is the
reference execution path required by ``EXPERIMENT.md``: independent of Warp,
yet "minimally transformed" - the transformation is the binding of ``m``.
"""

from __future__ import annotations

import contextlib
import importlib
import importlib.util
import sys
from collections.abc import Iterator, Sequence
from types import ModuleType

from ionmc.backend import mathlib


@contextlib.contextmanager
def _rebound_dependencies(binding: str, names: Sequence[str]) -> Iterator[None]:
    """Temporarily install same-binding copies of shared-source dependencies.

    A shared-source physics module that ``from ionmc.physics.X import f`` calls
    another shared-source module ``X`` must, in a non-default binding, pick up
    the *same-binding* ``X`` rather than the default (Warp) one. For the
    duration of the block, ``sys.modules[X]`` points at ``X`` loaded under
    ``binding``; the originals are restored afterwards. Reentrant within one
    thread; the loads are cached under their ``@binding`` aliases.
    """
    saved: dict[str, ModuleType | None] = {}
    try:
        for name in names:
            bound = load_bound_module(name, binding)
            saved[name] = sys.modules.get(name)
            sys.modules[name] = bound
        yield
    finally:
        for name, original in saved.items():
            if original is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = original


def load_bound_module(
    qualified_name: str,
    binding: str,
    rebind_dependencies: Sequence[str] = (),
) -> ModuleType:
    """Return ``qualified_name`` executed under ``binding``.

    Parameters
    ----------
    qualified_name:
        Dotted name of a shared-source module, e.g. ``"ionmc.physics.stopping"``.
    binding:
        ``"python"``, ``"numpy"`` or ``"warp"``.
    rebind_dependencies:
        Other shared-source modules that ``qualified_name`` imports functions
        from and calls; each is loaded under ``binding`` and installed under
        its real name while ``qualified_name`` is executed, so the cross-module
        calls stay within one binding (e.g. transport calling the tabulated
        stopping power). The dependencies are part of the cache key.

    The returned module is registered in ``sys.modules`` as
    ``"<qualified_name>@<binding>"`` (with a dependency suffix when rebinding)
    and cached, so repeated calls are cheap. The ``warp`` binding simply
    returns the normally imported module when Warp is the default binding.
    """
    if binding == "warp" and mathlib.default_binding() == "warp":
        return importlib.import_module(qualified_name)

    dep_key = ("+" + ",".join(rebind_dependencies)) if rebind_dependencies else ""
    alias = f"{qualified_name}@{binding}{dep_key}"
    cached = sys.modules.get(alias)
    if cached is not None:
        return cached

    spec = importlib.util.find_spec(qualified_name)
    if spec is None or spec.origin is None:
        raise ImportError(f"cannot locate source of {qualified_name!r}")

    new_spec = importlib.util.spec_from_file_location(alias, spec.origin)
    if new_spec is None or new_spec.loader is None:
        raise ImportError(f"cannot build a loader for {qualified_name!r}")
    module = importlib.util.module_from_spec(new_spec)
    module.__package__ = qualified_name.rpartition(".")[0]
    sys.modules[alias] = module
    try:
        with mathlib.bind(binding), _rebound_dependencies(binding, rebind_dependencies):
            new_spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(alias, None)
        raise
    return module


def python_reference(qualified_name: str) -> ModuleType:
    """Shortcut for ``load_bound_module(name, "python")``."""
    return load_bound_module(qualified_name, "python")


def numpy_reference(qualified_name: str) -> ModuleType:
    """Shortcut for ``load_bound_module(name, "numpy")``."""
    return load_bound_module(qualified_name, "numpy")
