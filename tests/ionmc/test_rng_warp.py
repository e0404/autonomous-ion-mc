"""Live check of the pure-Python RNG mirror against a Warp kernel.

``tests/ionmc/test_rng.py`` pins the mirror to values dumped on the experiment
workstation; this module re-derives them from the installed Warp on every
run where Warp is available (CPU device), so a Warp upgrade that changes the
generator or the operand order of ``randn`` is caught automatically
(decision 0005, *Expected tradeoffs*).
"""

from __future__ import annotations

import math

import pytest

from ionmc import rng

wp = pytest.importorskip("warp")

pytestmark = pytest.mark.warp

N_STREAMS = 64
SEED = 20260910


@wp.kernel
def _dump_kernel(
    seed: int,
    states: wp.array(dtype=wp.uint32),  # type: ignore[valid-type]
    floats: wp.array(dtype=float),  # type: ignore[valid-type]
    ints: wp.array(dtype=wp.uint32),  # type: ignore[valid-type]
    normals: wp.array(dtype=float),  # type: ignore[valid-type]
    scaled: wp.array(dtype=float),  # type: ignore[valid-type]
    bounded: wp.array(dtype=int),  # type: ignore[valid-type]
):
    i = wp.tid()
    state = wp.rand_init(seed, i)
    states[i] = state
    floats[i] = wp.randf(state)
    ints[i] = wp.randu(state)
    normals[i] = wp.randn(state)
    scaled[i] = wp.randf(state, -2.0, 3.0)
    bounded[i] = wp.randi(state, 10, 20)


@pytest.fixture(scope="module")
def kernel_dump(warp_module):
    wp.init()
    device = "cpu"
    arrays = {
        "states": wp.zeros(N_STREAMS, dtype=wp.uint32, device=device),
        "floats": wp.zeros(N_STREAMS, dtype=float, device=device),
        "ints": wp.zeros(N_STREAMS, dtype=wp.uint32, device=device),
        "normals": wp.zeros(N_STREAMS, dtype=float, device=device),
        "scaled": wp.zeros(N_STREAMS, dtype=float, device=device),
        "bounded": wp.zeros(N_STREAMS, dtype=int, device=device),
    }
    wp.launch(
        _dump_kernel,
        dim=N_STREAMS,
        inputs=[SEED, *arrays.values()],
        device=device,
    )
    wp.synchronize_device(device)
    return {k: v.numpy() for k, v in arrays.items()}


def test_states_and_integers_are_bit_exact(kernel_dump) -> None:
    mirror_states = [rng.rand_init(SEED, i) for i in range(N_STREAMS)]
    assert mirror_states == [int(v) for v in kernel_dump["states"]]
    for i in range(N_STREAMS):
        state = rng.RandomState(SEED, i)
        state.randf()
        assert state.randu() == int(kernel_dump["ints"][i])


def test_uniforms_are_exact_and_normals_float32_close(kernel_dump) -> None:
    for i in range(N_STREAMS):
        state = rng.RandomState(SEED, i)
        assert state.randf() == float(kernel_dump["floats"][i])
        state.randu()
        assert math.isclose(
            state.randn(), float(kernel_dump["normals"][i]), abs_tol=2e-6
        )
        # Warp scales in float32 ((high - low) * u + low); the mirror in float64.
        assert math.isclose(
            state.randf(-2.0, 3.0), float(kernel_dump["scaled"][i]), abs_tol=1e-6
        )
        assert state.randi(10, 20) == int(kernel_dump["bounded"][i])
