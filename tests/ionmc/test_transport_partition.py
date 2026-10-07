"""T13: partition invariance (workers, chunks, history ranges) and the fail-closed pool.

Per-history tally rows are reduced exactly, so counters and tallies are bit-identical for any
partition of the histories; the deposit grids are int64 fixed-point quanta (integer atomic
adds) and are bit-identical as well.
"""

from __future__ import annotations

import math
import multiprocessing
import os
from collections.abc import Callable
from fractions import Fraction

import numpy as np
import pytest

from ionmc.config import DiagnosticsOptions, SimulationConfig, validate
from ionmc.errors import TransportWorkerError, UnsupportedCombinationError
from ionmc.simulation import Simulation
from ionmc.transport.parity import compare_partition, format_partition_verdict
from ionmc.transport.pool import history_ranges, run_pool
from ionmc.transport.reference import run_reference_range
from ionmc.transport.tally import (
    PartialTransport,
    exact_components,
    merge_partials,
    rows_to_partial,
)


def _workers(n: int) -> int:
    """Worker processes used only for speed: 1 in the single-process diagnostic mode (histories
    and seeds unchanged)."""
    return 1 if os.environ.get("IONMC_SINGLE_PROCESS") == "1" else n


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
        np.ones((n, 6 + 2 * n_grids)),
        np.zeros((n, 9), dtype=np.int32),
        [np.zeros((2, 3), dtype=np.int64) for _ in range(n_grids)],
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
    """T13 (python): any split of the histories into ranges gives bit-identical tallies, counters
    and deposit grid (exact summation, int64 fixed-point grids). A failure reports the seed, the
    partition and the first differing field or voxel.

    OPEN DETERMINISM OBSERVATION (2026-10-04): once, during development of V3-003B, this test
    failed in a full-suite run: the "whole" run contained an escaped history (about 10 MeV in the
    ``escaped`` tally) that the split run did not. It was not reproduced afterwards: 1,920
    histories compared bitwise twice, 30 reruns of this test alone, 12 reruns after the scoring
    tests, 3 cold Warp caches and 192,000 Python-scope Warp calls under heap churn were all
    identical. The cause is unknown; together with an intermittent SIGABRT inside Python-scope
    Warp evaluation seen on the GPU host it is treated as an open observation about the
    reliability of the Python-scope ``@wp.func`` path. The repeatability steps of the LV and HR
    suites (``t-r1``) and ``test_repeatability_two_runs_identical`` watch for it."""
    seed, partition = 3, ((0, 5), (5, 6), (6, 12))
    cfg = make_config(
        energy=12.0, n=12, n_batches=3, energy_sigma=0.2, lateral_sigma=0.5, seed=seed
    )
    eff = validate(cfg)
    whole_part = run_reference_range(eff, 0, 12)
    split_parts = [run_reference_range(eff, a, b) for a, b in partition]
    whole = merge_partials([whole_part], 12, 1)
    split = merge_partials(split_parts, 12, 1)
    # per-history rows are independent of the partition: report the first difference exactly
    where = f"seed={seed} partition={partition}"
    for k in whole.tallies:
        assert split.tallies[k] == whole.tallies[k], (
            f"{where}: tally {k}: whole={whole.tallies[k]!r} split={split.tallies[k]!r}; "
            f"whole parts components={[p.tally_components for p in [whole_part]]} "
            f"split components={[p.tally_components for p in split_parts]}"
        )
    assert split.outside_mev == whole.outside_mev, where
    assert split.counters == whole.counters, f"{where}: {split.counters} vs {whole.counters}"
    diff = np.argwhere(split.edep_mev[0] != whole.edep_mev[0])
    assert diff.size == 0, (
        f"{where}: first differing voxel (batch, voxel) {diff[0].tolist()}: "
        f"whole={whole.edep_mev[0][tuple(diff[0])]!r} split={split.edep_mev[0][tuple(diff[0])]!r}"
    )


def _compare(a: object, b: object) -> dict:  # type: ignore[type-arg]
    verdict = compare_partition(a, b)  # type: ignore[arg-type]
    return verdict


@pytest.mark.multiprocess
@pytest.mark.parametrize("precision", ["float32", "float64"])
def test_t13_warp_cpu_one_versus_three_workers(make_config: MakeConfig, precision: str) -> None:
    """T13 (CI): 1 versus 3 worker processes on warp-cpu: counters, tallies and int64 deposit
    grids bit-identical (both precisions); every worker loaded the cached kernel."""

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
    assert verdict["pass"], format_partition_verdict(verdict)
    parts = three.transport_report["partials"]
    assert [p["worker"] for p in parts] == [0, 1, 2]
    assert len({p["pid"] for p in parts}) == 3
    # the parent compiled the kernel before spawning, so workers only load it from the cache
    assert all(p["compile_s"] < 0.5 for p in parts), [p["compile_s"] for p in parts]
    assert three.device.endswith("x 3 processes")
    # the pooled compile time includes the parent's compile before the workers were spawned
    assert all("parent_compile_s" in p for p in parts)
    assert three.timings["compile"] >= max(p["parent_compile_s"] for p in parts)
    assert three.energy_balance.relative_residual < (1e-5 if precision == "float32" else 1e-12)


def test_t13_warp_cpu_chunk_sizes(make_config: MakeConfig) -> None:
    """Chunking a launch (2**10 versus a single chunk) changes nothing: counters, tallies and the
    int64 deposit grid are bit-identical."""

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


@pytest.mark.multiprocess
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


@pytest.mark.multiprocess
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


def test_repeatability_two_runs_identical(make_config: MakeConfig) -> None:
    """CI-size repeatability (python 2 x 40 histories at 40 MeV, 2 workers): the
    same configuration and seed give bit-identical tallies, counters, deposit grids and
    per-history end states; an altered seed is detected (the comparison has power)."""
    from ionmc.config import DiagnosticsOptions
    from ionmc.simulation import Simulation
    from ionmc.transport.parity import compare_runs_bitwise

    def run(seed: int):  # type: ignore[no-untyped-def]
        cfg = make_config(
            energy=40.0, n=40, n_batches=20, seed=seed, lateral_sigma=0.5, energy_sigma=0.3,
            diagnostics=DiagnosticsOptions(track_end_positions=True),
            run_kwargs={"cpu_workers": _workers(2), "worker_timeout_s": 600.0},
        )  # fmt: skip
        return Simulation(cfg).run()

    a, b, c = run(9), run(9), run(10)
    same = compare_runs_bitwise(a, b)
    assert same["identical"], same["first_difference"]
    other = compare_runs_bitwise(a, c)
    assert not other["identical"] and other["first_difference"] is not None
