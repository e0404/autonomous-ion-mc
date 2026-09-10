"""Tests for the pure-Python mirror of Warp's RNG.

The reference values were dumped from a Warp 1.17.0 kernel executing
``wp.rand_init(42, i)`` followed by one ``wp.randf``, ``wp.randu`` and
``wp.randn`` on the experiment workstation (Warp CPU and CUDA gave identical
states, floats and ints; host run ``RUN-20260910T004325Z-472e71df``, task
DEV-002). They pin the mirror to the kernel behaviour bit-for-bit.
"""

from __future__ import annotations

import math

import pytest

from ionmc import rng

KERNEL_STATES = [
    2041448361,
    2311944962,
    1095600430,
    3208355652,
    1945815008,
    1281180028,
    2094792678,
    3115480925,
]
KERNEL_FLOATS = [
    0.8659751415252686,
    0.013862192630767822,
    0.4421435594558716,
    0.10313504934310913,
    0.7510697245597839,
    0.8886032104492188,
    0.4199472665786743,
    0.49427831172943115,
]
KERNEL_INTS = [
    798497746,
    3339178658,
    4037432846,
    1145322783,
    1826705985,
    2281956139,
    4518472,
    1413447399,
]
KERNEL_NORMALS = [
    0.7695993781089783,
    -0.05680134519934654,
    -1.0362950563430786,
    -0.4782019555568695,
    -0.3645292818546295,
    -0.6014437675476074,
    0.18708495795726776,
    -0.8180838823318481,
]


def test_rand_init_matches_kernel_states() -> None:
    assert [rng.rand_init(42, i) for i in range(8)] == KERNEL_STATES


@pytest.mark.parametrize("i", range(8))
def test_randf_randu_randn_sequence_matches_kernel(i: int) -> None:
    # The kernel drew randf, then randu, then randn from the same state.
    state = rng.RandomState(42, i)
    assert state.randf() == KERNEL_FLOATS[i]
    assert state.randu() == KERNEL_INTS[i]
    # randn is evaluated in float32 by Warp: agreement to ~1e-6 absolute.
    assert math.isclose(state.randn(), KERNEL_NORMALS[i], abs_tol=2.0e-6)


def test_randf_range_and_resolution() -> None:
    state = rng.RandomState(7)
    values = [state.randf() for _ in range(10_000)]
    assert all(0.0 <= v < 1.0 for v in values)
    # 24 random bits: every value is a multiple of 2**-24.
    assert all((v * 16777216.0).is_integer() for v in values)


def test_randf_scaled_and_randi_bounds() -> None:
    state = rng.RandomState(3, 5)
    for _ in range(1000):
        v = state.randf(-2.0, 3.0)
        assert -2.0 <= v < 3.0
    for _ in range(1000):
        k = state.randi(10, 20)
        assert 10 <= k < 20


def test_streams_are_deterministic_and_distinct() -> None:
    a = [rng.RandomState(11, 4).randf() for _ in range(3)]
    b = [rng.RandomState(11, 4).randf() for _ in range(3)]
    c = [rng.RandomState(11, 5).randf() for _ in range(3)]
    assert a == b
    assert a != c


def test_rand_init_single_argument_form() -> None:
    assert rng.rand_init(42) == rng.rand_pcg(42)
    assert rng.RandomState(42).state == rng.rand_pcg(42)
