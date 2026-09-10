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

import importlib
import importlib.util
import sys
from types import ModuleType

from ionmc.backend import mathlib


def load_bound_module(qualified_name: str, binding: str) -> ModuleType:
    """Return ``qualified_name`` executed under ``binding``.

    Parameters
    ----------
    qualified_name:
        Dotted name of a shared-source module, e.g. ``"ionmc.physics.stopping"``.
    binding:
        ``"python"``, ``"numpy"`` or ``"warp"``.

    The returned module is registered in ``sys.modules`` as
    ``"<qualified_name>@<binding>"`` and cached, so repeated calls are cheap.
    The ``warp`` binding simply returns the normally imported module when
    Warp is the default binding.
    """
    if binding == "warp" and mathlib.default_binding() == "warp":
        return importlib.import_module(qualified_name)

    alias = f"{qualified_name}@{binding}"
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
        with mathlib.bind(binding):
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
