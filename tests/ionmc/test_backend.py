"""Tests of the shared-source execution mechanism."""

from __future__ import annotations

import math
import sys

import numpy as np
import pytest

from ionmc.backend import mathlib, reference


def test_python_namespace_semantics() -> None:
    m = mathlib.PYTHON_MATH
    assert m.where(True, 1.0, 2.0) == 1.0
    assert m.where(False, 1.0, 2.0) == 2.0
    assert m.max(1.0, 2.0) == 2.0 and m.min(1.0, 2.0) == 1.0
    assert m.sqrt(4.0) == 2.0 and m.log(math.e) == pytest.approx(1.0)
    assert m.pow(2.0, 3.0) == 8.0


def test_numpy_namespace_is_vectorised() -> None:
    m = mathlib.numpy_math()
    x = np.array([0.5, 2.0])
    np.testing.assert_array_equal(m.where(x > 1.0, x, -x), np.array([-0.5, 2.0]))
    np.testing.assert_array_equal(m.max(x, 1.0), np.array([1.0, 2.0]))


def test_bind_selects_current_namespace() -> None:
    default = mathlib.default_binding()
    with mathlib.bind("python"):
        assert mathlib.current().name == "python"
        assert mathlib.func(lambda x: x) is not None
        with mathlib.bind("numpy"):
            assert mathlib.current().name == "numpy"
        assert mathlib.current().name == "python"
    assert mathlib._active_binding[-1] == default
    with pytest.raises(ValueError):
        with mathlib.bind("fortran"):
            pass


def test_load_bound_module_creates_independent_copies() -> None:
    py = reference.load_bound_module("ionmc.physics.stopping", "python")
    npy = reference.load_bound_module("ionmc.physics.stopping", "numpy")
    assert py is not npy
    assert py.m.name == "python" and npy.m.name == "numpy"
    assert "ionmc.physics.stopping@python" in sys.modules
    # cached
    assert reference.python_reference("ionmc.physics.stopping") is py
    # same source, same numbers
    assert py.beta_squared(100.0, 938.272) == npy.beta_squared(100.0, 938.272)
    with pytest.raises(ImportError):
        reference.load_bound_module("ionmc.does_not_exist", "python")


def test_default_binding_matches_warp_availability() -> None:
    assert mathlib.default_binding() == ("warp" if mathlib.HAVE_WARP else "python")
    if not mathlib.HAVE_WARP:
        with pytest.raises(ImportError):
            mathlib.warp_module()
