"""Tests of the reproducible benchmarking harness (decision 0030).

These cover the harness primitives (provenance, timing statistics, the scientific
digest, cross-backend agreement, JSON serialisation) on the reference path only, so
they are fast and need neither Warp nor a GPU. The actual performance numbers are
never asserted (hardware-dependent); the physics digest and the harness contract are.
"""

from __future__ import annotations

import io

import numpy as np
import pytest

from ionmc import benchmarking as bm


def test_capture_provenance_has_required_keys() -> None:
    prov = bm.capture_provenance()
    for key in (
        "ionmc_version",
        "python",
        "numpy",
        "platform",
        "cpu_count",
        "warp_version",
        "warp_devices",
    ):
        assert key in prov
    assert isinstance(prov["warp_devices"], list)


def test_measure_statistics_and_throughput() -> None:
    counter = {"n": 0}

    def fn() -> None:
        counter["n"] += 1

    t = bm.measure(fn, repeats=5, warmup=2, work_units=1000.0, backend="reference")
    # warmup + repeats calls, and only repeats are timed
    assert counter["n"] == 7
    assert len(t.seconds) == 5
    assert t.seconds_min <= t.seconds_median <= max(t.seconds)
    assert t.throughput_per_s == pytest.approx(1000.0 / t.seconds_median)
    assert t.backend == "reference"


def test_measure_rejects_bad_arguments() -> None:
    with pytest.raises(ValueError):
        bm.measure(lambda: None, repeats=0)
    with pytest.raises(ValueError):
        bm.measure(lambda: None, warmup=-1)


def test_measure_runs_sync_callback() -> None:
    calls = {"sync": 0}

    def sync() -> None:
        calls["sync"] += 1

    bm.measure(lambda: None, repeats=3, warmup=1, sync=sync)
    assert calls["sync"] == 4  # once per warmup + per timed repeat


def test_array_digest_is_deterministic_and_sensitive() -> None:
    a = np.linspace(0.0, 1.0, 50)
    assert bm.array_digest(a) == bm.array_digest(a.copy())
    b = a.copy()
    b[10] += 1e-3  # above the 6-decimal rounding floor
    assert bm.array_digest(a) != bm.array_digest(b)
    # a sub-rounding perturbation does not change the digest (stability)
    c = a.copy()
    c[10] += 1e-9
    assert bm.array_digest(a) == bm.array_digest(c)


def test_relative_agreement() -> None:
    a = np.array([0.0, 1.0, 2.0, 1.0])
    assert bm.relative_agreement(a, a) == {
        "total_rel_diff": 0.0,
        "cumulative_rel_diff": 0.0,
        "max_bin_rel_diff": 0.0,
    }
    b = a.copy()
    b[2] += 0.2  # peak is 2.0, total is 4.0
    agree = bm.relative_agreement(a, b)
    assert agree["max_bin_rel_diff"] == pytest.approx(0.1)
    assert agree["total_rel_diff"] == pytest.approx(0.2 / 4.0)
    # cumulative: cumsum diff is 0.2 from bin 2 onward -> max 0.2, /total 4.0
    assert agree["cumulative_rel_diff"] == pytest.approx(0.2 / 4.0)


def test_agreement_catches_uniform_scale_divergence() -> None:
    """A uniform-scale divergence (a units / normalisation / accumulation-constant
    bug) must fail the gate on arrays that carry absolute magnitude. This is the
    guarantee the benchmark physics gate relies on: comparing per-history curves,
    not unit-normalised shapes (which would hide any constant factor)."""
    a = np.array([1.0, 3.0, 5.0, 2.0])
    scaled = a * 1.10  # 10 % more total energy, identical shape
    agree = bm.relative_agreement(a, scaled)
    assert agree["total_rel_diff"] == pytest.approx(0.10)
    assert agree["cumulative_rel_diff"] > 0.05  # would be ~0 on unit-normalised curves
    # sanity: normalising away the scale hides it -> shows why we must not do that
    norm = bm.relative_agreement(a / a.sum(), scaled / scaled.sum())
    assert norm["total_rel_diff"] < 1e-12


def test_make_sync_cuda_branch(monkeypatch) -> None:
    """For a non-CPU Warp device, make_sync returns a callable that synchronises
    that device (the honest-GPU-timing path), even without a real GPU present."""
    if not bm.mathlib.HAVE_WARP:
        pytest.skip("Warp not installed")
    calls: list[str] = []

    class _FakeWarp:
        def synchronize_device(self, device: str) -> None:
            calls.append(device)

    monkeypatch.setattr(bm.mathlib, "warp_module", lambda: _FakeWarp())
    sync = bm.make_sync("warp", "cuda:0")
    assert sync is not None
    sync()
    assert calls == ["cuda:0"]


def test_backend_label_and_sync_selection() -> None:
    assert bm.backend_label("python", None) == "reference"
    assert bm.backend_label("warp", "cpu") == "warp:cpu"
    assert bm.backend_label("warp", "cuda:0") == "warp:cuda:0"
    # reference and Warp CPU are synchronous -> no sync callback
    assert bm.make_sync("python", None) is None
    assert bm.make_sync("warp", "cpu") is None


def test_dump_report_writes_valid_json_with_numpy() -> None:
    report = bm.new_report("unit")
    report["value"] = np.float32(1.5)
    report["array"] = np.arange(3)
    buf = io.StringIO()
    bm.dump_report(report, buf)
    import json

    parsed = json.loads(buf.getvalue())
    assert parsed["benchmark"] == "unit"
    assert parsed["value"] == 1.5
    assert parsed["array"] == [0, 1, 2]
    assert parsed["schema_version"] == 1


def test_end_to_end_reference_benchmark(pstar_cache_root) -> None:
    """A tiny real depth-dose benchmark on the reference path produces a positive
    throughput and a stable digest (the driver's own building blocks)."""
    from ionmc.data import MCSQUARE_PSTAR_WATER, cache
    from ionmc.data.stopping_tables import load_stopping_table
    from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine
    from ionmc.transport.geometry import WaterSlab

    proton = load_stopping_table(
        cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root)
    )
    grid = DepthDoseGrid(250.0, 500)
    eng = TransportEngine(
        proton, WaterSlab(250.0), grid, straggling=False, nuclear=False
    )
    src = PencilBeamSource(150.0)

    r1 = eng.run(src, 5, seed=1, path="python").edep_mev
    r2 = eng.run(src, 5, seed=1, path="python").edep_mev
    # deterministic: identical digests across runs
    assert bm.array_digest(r1 / r1.sum()) == bm.array_digest(r2 / r2.sum())

    t = bm.measure(
        lambda: eng.run(src, 5, seed=1, path="python"),
        repeats=2,
        warmup=1,
        work_units=5.0,
        backend="reference",
    )
    assert t.throughput_per_s is not None and t.throughput_per_s > 0.0
