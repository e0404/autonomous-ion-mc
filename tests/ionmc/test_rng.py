"""U6: Philox4x32-10 known answers, kernel/Python identity, uniform mapping, counters."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest
import warp as wp

from ionmc.config import DiagnosticsOptions, SimulationConfig
from ionmc.errors import CounterOverflowError
from ionmc.rng.philox import (
    MAX_GENEALOGY_ID,
    PURPOSE_SOURCE,
    PURPOSE_TRANSPORT,
    child_genealogy_id,
    draw_block,
    key_from_seed,
    make_counter,
    make_philox,
    philox4x32_10_py,
    u01_py,
)
from ionmc.simulation import Simulation

wp.config.log_level = wp.LOG_WARNING

PH32 = make_philox(wp.float32)
PH64 = make_philox(wp.float64)

# Random123 known-answer vectors for philox4x32 with 10 rounds (kat_vectors).
KAT = [
    ((0, 0, 0, 0), (0, 0), (0x6627E8D5, 0xE169C58D, 0xBC57AC4C, 0x9B00DBD8)),
    (
        (0xFFFFFFFF,) * 4,
        (0xFFFFFFFF, 0xFFFFFFFF),
        (0x408F276D, 0x41C83B0E, 0xA20BC7C6, 0x6D5451FD),
    ),
    (
        (0x243F6A88, 0x85A308D3, 0x13198A2E, 0x03707344),
        (0xA4093822, 0x299F31D0),
        (0xD16CFE09, 0x94FDCCEB, 0x5001E420, 0x24126EA1),
    ),
]


@wp.kernel(module="unique")
def _block_kernel(
    cs: wp.array(dtype=wp.vec4ui),  # type: ignore[valid-type]
    ks: wp.array(dtype=wp.vec2ui),  # type: ignore[valid-type]
    out: wp.array(dtype=wp.vec4ui),  # type: ignore[valid-type]
    u32: wp.array2d(dtype=wp.float32),  # type: ignore[valid-type]
    u64: wp.array2d(dtype=wp.float64),  # type: ignore[valid-type]
) -> None:
    i = wp.tid()
    r = PH64.philox4x32_10(cs[i], ks[i])
    out[i] = r
    for j in range(4):
        u32[i, j] = PH32.u01(r[j])
        u64[i, j] = PH64.u01(r[j])


def _run_blocks(
    counters: np.ndarray, keys: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = counters.shape[0]
    out = wp.zeros(n, dtype=wp.vec4ui, device="cpu")
    u32 = wp.zeros((n, 4), dtype=wp.float32, device="cpu")
    u64 = wp.zeros((n, 4), dtype=wp.float64, device="cpu")
    wp.launch(
        _block_kernel,
        dim=n,
        inputs=[
            wp.array(counters, dtype=wp.vec4ui, device="cpu"),
            wp.array(keys, dtype=wp.vec2ui, device="cpu"),
            out,
            u32,
            u64,
        ],
        device="cpu",
    )
    return out.numpy(), u32.numpy(), u64.numpy()


def test_u6a_random123_known_answers_python_and_kernel() -> None:
    counters = np.array([k[0] for k in KAT], dtype=np.uint32)
    keys = np.array([k[1] for k in KAT], dtype=np.uint32)
    kernel, _, _ = _run_blocks(counters, keys)
    for i, (c, k, expected) in enumerate(KAT):
        assert philox4x32_10_py(c, k) == expected
        assert tuple(int(x) for x in kernel[i]) == expected


def test_u6b_kernel_and_python_are_bit_identical_on_20000_random_blocks() -> None:
    rng = np.random.default_rng(12345)
    n = 20000
    counters = rng.integers(0, 2**32, (n, 4), dtype=np.uint64).astype(np.uint32)
    keys = rng.integers(0, 2**32, (n, 2), dtype=np.uint64).astype(np.uint32)
    kernel, u32, u64 = _run_blocks(counters, keys)
    expected = np.array(
        [
            philox4x32_10_py(tuple(int(x) for x in counters[i]), tuple(int(x) for x in keys[i]))
            for i in range(n)
        ],
        dtype=np.uint32,
    )
    assert np.array_equal(kernel, expected)
    # the uniform maps agree with the Python mapping bit for bit
    w = expected.astype(np.uint64)
    ref64 = ((w >> np.uint64(8)).astype(np.float64) + 0.5) * 2.0**-24
    ref32 = (((w >> np.uint64(9)).astype(np.float32)) + np.float32(0.5)) * np.float32(2.0**-23)
    assert np.array_equal(u64, ref64)
    assert np.array_equal(u32, ref32)
    assert u01_py(int(expected[0, 0]), "float64") == ref64[0, 0]
    assert u01_py(int(expected[0, 0]), "float32") == float(ref32[0, 0])


def test_u6c_u01_strictly_inside_unit_interval_in_both_precisions() -> None:
    words = np.array([0, 0xFF, 0x1FF, 0x7FFFFFFF, 0xFFFFFF00, 0xFFFFFFFF], dtype=np.uint32)
    for w in words:
        for precision in ("float32", "float64"):
            u = u01_py(int(w), precision)
            assert 0.0 < u < 1.0
    # kernel check including the word that rounds to exactly 1.0 with the 24-bit float32 form
    n = len(words)
    counters = np.zeros((n, 4), dtype=np.uint32)
    keys = np.zeros((n, 2), dtype=np.uint32)
    _, u32, u64 = _run_blocks(counters, keys)
    assert np.all((u32 > 0.0) & (u32 < 1.0)) and np.all((u64 > 0.0) & (u64 < 1.0))
    # the regression of decision 0039: the 24-bit form in float32 would give exactly 1.0
    assert np.float32((0xFFFFFFFF >> 8) + 0.5) * np.float32(2.0**-24) == np.float32(1.0)
    assert u01_py(0xFFFFFFFF, "float32") < 1.0
    assert u01_py(0xFFFFFFFF, "float64") < 1.0
    with pytest.raises(ValueError):
        u01_py(2**32, "float64")
    with pytest.raises(ValueError):
        u01_py(0, "float16")


def test_u6d_counter_streams_are_disjoint_and_partition_invariant(
    make_config: Callable[..., SimulationConfig],
) -> None:
    n = 200
    diag = DiagnosticsOptions(trace_histories=n)
    used: dict[int, set[tuple[int, int, int, int]]] = {}
    block_counts: dict[tuple[int, int], np.ndarray] = {}
    for seed, n_batches in ((7, 2), (8, 2), (7, 4)):
        cfg = make_config(
            energy=4.0,
            n=n,
            n_batches=n_batches,
            seed=seed,
            diagnostics=diag,
            physics_kwargs={"max_energy_loss_fraction": 0.2},
        )
        tr = Simulation(cfg).run().diagnostics["trace"]
        last = np.zeros(n, dtype=int)  # Philox blocks drawn by each history
        for h, b in zip(tr["history"].astype(int), tr["blocks"].astype(int), strict=True):
            last[h] = max(last[h], b)
        block_counts[(seed, n_batches)] = last
        if n_batches == 2:
            counters: list[tuple[int, int, int, int]] = []
            for h in range(n):
                counters.append(make_counter(h, 0, 0, PURPOSE_SOURCE))
                counters.extend(make_counter(h, 0, k, PURPOSE_TRANSPORT) for k in range(last[h]))
            assert len(set(counters)) == len(counters)  # no counter is used twice
            assert {c[3] for c in counters} == {PURPOSE_SOURCE, PURPOSE_TRANSPORT}
            used[seed] = set(counters)
    # the same histories use the same counters whatever the batch partition
    assert np.array_equal(block_counts[(7, 2)], block_counts[(7, 4)])
    # different seeds: different keys, hence different streams for the same counters
    assert key_from_seed(7) != key_from_seed(8)
    assert draw_block(key_from_seed(7), 0, 0, 0, 0) != draw_block(key_from_seed(8), 0, 0, 0, 0)
    # the full 128-bit outputs of all used counters under one key are pairwise distinct
    outs = {philox4x32_10_py(c, key_from_seed(7)) for c in used[7]}
    assert len(outs) == len(used[7])


def test_u6e_counter_bounds_raise() -> None:
    assert make_counter(2**32 - 1, MAX_GENEALOGY_ID - 1, 2**32 - 1, 2) == (
        2**32 - 1,
        MAX_GENEALOGY_ID - 1,
        2**32 - 1,
        2,
    )
    for bad in (
        dict(history=2**32, genealogy_id=0, block=0, purpose=0),
        dict(history=-1, genealogy_id=0, block=0, purpose=0),
        dict(history=0, genealogy_id=MAX_GENEALOGY_ID, block=0, purpose=0),
        dict(history=0, genealogy_id=0, block=2**32, purpose=0),
        dict(history=0, genealogy_id=0, block=0, purpose=3),
    ):
        with pytest.raises(CounterOverflowError):
            make_counter(**bad)  # type: ignore[arg-type]
    for bad_seed in (-1, 2**64, 1.5, True):
        with pytest.raises(CounterOverflowError):
            key_from_seed(bad_seed)  # type: ignore[arg-type]
    assert key_from_seed(2**64 - 1) == (0xFFFFFFFF, 0xFFFFFFFF)
    assert key_from_seed(0x1_0000_0002) == (2, 1)
    # genealogy encoding: id = parent + b * 32**g, 31 children, 6 generations
    assert child_genealogy_id(0, 0, 1) == 1
    assert child_genealogy_id(0, 0, 31) == 31
    assert child_genealogy_id(5, 1, 3) == 5 + 3 * 32
    last = 0
    for g in range(6):
        last = child_genealogy_id(last, g, 31)
    assert last == 32**6 - 1 == MAX_GENEALOGY_ID - 1
    with pytest.raises(CounterOverflowError):
        child_genealogy_id(0, 0, 32)  # a 32nd child
    with pytest.raises(CounterOverflowError):
        child_genealogy_id(0, 0, 0)
    with pytest.raises(CounterOverflowError):
        child_genealogy_id(last, 6, 1)  # a 7th generation
    # all ids of a full tree are distinct (disjointness of the encoding)
    ids = {0}
    frontier = [(0, 0)]
    for _g in range(3):
        nxt = []
        for pid, gen in frontier:
            for b in (1, 2, 31):
                cid = child_genealogy_id(pid, gen, b)
                assert cid not in ids
                ids.add(cid)
                nxt.append((cid, gen + 1))
        frontier = nxt
