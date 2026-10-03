"""T13: partition invariance (workers, chunks, history ranges) and the fail-closed pool.

Per-history tally rows are reduced exactly, so counters and tallies are bit-identical for any
partition of the histories; the deposit grids differ only by the order of float additions.
"""

from __future__ import annotations

import math
import multiprocessing
from collections.abc import Callable
from fractions import Fraction

import numpy as np
import pytest

from ionmc.config import DiagnosticsOptions, SimulationConfig, validate
from ionmc.errors import TransportWorkerError, UnsupportedCombinationError
from ionmc.simulation import Simulation
from ionmc.transport.parity import compare_partition
from ionmc.transport.pool import history_ranges, run_pool
from ionmc.transport.reference import run_reference_range
from ionmc.transport.tally import (
    PartialTransport,
    exact_components,
    merge_partials,
    rows_to_partial,
)

MakeConfig = Callable[..., SimulationConfig]


def test_exact_components_equal_the_exact_sum() -> None:
    rng = np.random.default_rng(7)
    x = rng.normal(size=5000) * np.exp(rng.uniform(-30, 30, 5000))
    x = np.concatenate([x, -x[:2500]])  # heavy cancellation
    exact = sum(Fraction(float(v)) for v in x)
    comps = exact_components(x)
    assert sum(Fraction(c) for c in comps) == exact
    assert math.fsum(comps) == float(exact)
    # splitting into chunks in any order gives the same correctly rounded sum
    parts = np.array_split(rng.permutation(x), 7)
    split = math.fsum(c for p in parts for c in exact_components(p))
    assert split == float(exact)
    assert exact_components(np.zeros(5)) == []
    with pytest.raises(ValueError, match="finite"):
        exact_components(np.array([1.0, np.nan]))


def test_history_ranges_tile() -> None:
    for n, w in ((10, 3), (7, 7), (100, 8), (1, 1)):
        r = history_ranges(n, w)
        assert r[0][0] == 0 and r[-1][1] == n and len(r) == w
        assert all(a[1] == b[0] for a, b in zip(r, r[1:], strict=False))
        assert all(h1 > h0 for h0, h1 in r)


def _partial(h0: int, h1: int, n_grids: int = 1) -> PartialTransport:
    n = h1 - h0
    return rows_to_partial(
        h0,
        h1,
        np.ones((n, 6 + n_grids)),
        np.zeros((n, 7), dtype=np.int32),
        [np.zeros((2, 3)) for _ in range(n_grids)],
        None,
    )


def test_merge_partials_fails_closed_on_gap_overlap_or_shortfall() -> None:
    merge_partials([_partial(0, 2), _partial(2, 5)], 5, 1)  # complete: fine
    with pytest.raises(ValueError, match="gap or overlap"):
        merge_partials([_partial(0, 2), _partial(3, 5)], 5, 1)
    with pytest.raises(ValueError, match="gap or overlap"):
        merge_partials([_partial(0, 3), _partial(2, 5)], 5, 1)
    with pytest.raises(ValueError, match="expected"):
        merge_partials([_partial(0, 2), _partial(2, 4)], 5, 1)
    with pytest.raises(ValueError, match="wrong number"):
        merge_partials([_partial(0, 5, n_grids=2)], 5, 1)


def test_python_reference_history_ranges_are_partition_independent(
    make_config: MakeConfig,
) -> None:
    """T13 (python): any split of the histories into ranges gives identical tallies and
    counters and a deposit grid equal to 1e-12 (the same float additions in another order)."""
    cfg = make_config(energy=12.0, n=12, n_batches=3, energy_sigma=0.2, lateral_sigma=0.5, seed=3)
    eff = validate(cfg)
    whole = merge_partials([run_reference_range(eff, 0, 12)], 12, 1)
    split = merge_partials(
        [run_reference_range(eff, a, b) for a, b in ((0, 5), (5, 6), (6, 12))], 12, 1
    )
    assert split.tallies == whole.tallies and split.outside_mev == whole.outside_mev
    assert split.counters == whole.counters
    np.testing.assert_allclose(split.edep_mev[0], whole.edep_mev[0], rtol=1e-12, atol=0.0)


def _compare(a: object, b: object) -> dict:  # type: ignore[type-arg]
    verdict = compare_partition(a, b)  # type: ignore[arg-type]
    return verdict


@pytest.mark.parametrize("precision", ["float32", "float64"])
def test_t13_warp_cpu_one_versus_three_workers(make_config: MakeConfig, precision: str) -> None:
    """T13 (CI): 1 versus 3 worker processes on warp-cpu: counters and tallies identical, deposit
    grid within 1e-5 (float32) or 1e-12 (float64); every worker loaded the cached kernel."""

    def run(workers: int):  # type: ignore[no-untyped-def]
        cfg = make_config(
            energy=60.0,
            n=720,
            n_batches=6,
            seed=11,
            energy_sigma=0.3,
            lateral_sigma=1.0,
            backend="warp-cpu",
            precision=precision,
            run_kwargs={"cpu_workers": workers, "worker_timeout_s": 300.0},
        )
        return Simulation(cfg).run()

    one, three = run(1), run(3)
    verdict = _compare(one, three)
    assert verdict["tallies_identical"] and verdict["counters_identical"]
    assert verdict["pass"], verdict["deposit"]
    parts = three.transport_report["partials"]
    assert [p["worker"] for p in parts] == [0, 1, 2]
    assert len({p["pid"] for p in parts}) == 3
    # the parent compiled the kernel before spawning, so workers only load it from the cache
    assert all(p["compile_s"] < 0.5 for p in parts), [p["compile_s"] for p in parts]
    assert three.device.endswith("x 3 processes")
    assert three.energy_balance.relative_residual < (1e-5 if precision == "float32" else 1e-12)


def test_t13_warp_cpu_chunk_sizes(make_config: MakeConfig) -> None:
    """Chunking a launch (2**10 versus a single chunk) changes nothing but the order of float
    adds in the deposit grid."""

    def run(chunk: int):  # type: ignore[no-untyped-def]
        cfg = make_config(
            energy=40.0,
            n=3000,
            n_batches=6,
            seed=5,
            backend="warp-cpu",
            precision="float32",
            run_kwargs={"chunk_histories": chunk},
        )
        return Simulation(cfg).run()

    small, large = run(2**10), run(2**18)
    assert small.transport_report["partials"][0]["n_chunks"] == 3
    assert large.transport_report["partials"][0]["n_chunks"] == 1
    verdict = _compare(small, large)
    assert verdict["tallies_identical"] and verdict["counters_identical"] and verdict["pass"]


def test_t13_python_one_versus_three_workers(make_config: MakeConfig) -> None:
    """The python backend also splits over spawned workers (T12 needs it): identical tallies,
    deposit grid to 1e-12, diagnostics concatenated in history order."""

    def run(workers: int):  # type: ignore[no-untyped-def]
        cfg = make_config(
            energy=15.0,
            n=12,
            n_batches=3,
            seed=21,
            lateral_sigma=0.5,
            diagnostics=DiagnosticsOptions(track_end_positions=True, trace_histories=5),
            run_kwargs={"cpu_workers": workers, "worker_timeout_s": 300.0},
        )
        return Simulation(cfg).run()

    one, three = run(1), run(3)
    verdict = _compare(one, three)
    assert verdict["tallies_identical"] and verdict["counters_identical"] and verdict["pass"]
    for key in ("end_position_mm", "end_code", "end_energy_mev"):
        assert np.array_equal(one.diagnostics[key], three.diagnostics[key]), key
    for col, v in one.diagnostics["trace"].items():
        assert np.array_equal(v, three.diagnostics["trace"][col]), col


@pytest.mark.parametrize("fault", ["raise", "exit", "hang"])
def test_pool_fails_closed_without_partial_result(make_config: MakeConfig, fault: str) -> None:
    """A worker that raises, dies or hangs terminates the whole run: no result, no stray
    processes (the hanging worker is killed at the timeout)."""
    cfg = make_config(
        energy=10.0,
        n=6,
        n_batches=2,
        backend="warp-cpu",
        precision="float32",
        run_kwargs={"cpu_workers": 3, "worker_timeout_s": 12.0},
    )
    eff = validate(cfg)
    with pytest.raises(TransportWorkerError) as err:
        run_pool(eff, _fault=fault)
    text = str(err.value)
    if fault == "raise":
        assert "injected worker failure" in text and "worker 1" in text
    elif fault == "exit":
        assert "exit code 3" in text or "exited with code 3" in text
    else:
        assert "worker_timeout_s" in text
    assert multiprocessing.active_children() == []


def test_run_options_chunk_workers_and_trace_rules(make_config: MakeConfig) -> None:
    with pytest.raises(UnsupportedCombinationError, match="power of two"):
        make_config(run_kwargs={"chunk_histories": 3000})
    with pytest.raises(UnsupportedCombinationError, match=">= 1024"):
        make_config(run_kwargs={"chunk_histories": 512})
    with pytest.raises(UnsupportedCombinationError, match="cpu_workers"):
        Simulation(make_config(n=2, run_kwargs={"cpu_workers": 4}))
    with pytest.raises(UnsupportedCombinationError, match="cpu_workers"):
        Simulation(make_config(n=4, run_kwargs={"cpu_workers": 1000}))
    eff = validate(make_config(n=4, run_kwargs={"cpu_workers": 2, "chunk_histories": 2048}))
    assert eff.summary()["chunk_histories"] == 2048 and eff.summary()["cpu_workers"] == 2
