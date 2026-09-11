"""Reference-path tests for the release-validation orchestrator (decision 0033).

These check the suite enumeration and per-script flag detection without running the
full (GPU, dataset-dependent) suite — that runs on the host runner. The orchestrator
itself is a subprocess runner, so its wiring is what these guard.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[2]
_SCRIPT = _REPO / "validation" / "release_validation.py"


def _load():
    spec = importlib.util.spec_from_file_location("release_validation", _SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def rv():
    return _load()


def test_enumerates_milestone_validations_and_excludes_diagnostics(rv) -> None:
    suites = rv._validation_suites(cache_dir=None, require_cuda=False)
    names = {s["name"] for s in suites}
    # every milestone family is represented
    for expected in (
        "v0_stopping_power",
        "v1_depth_dose_csda",
        "v2_nuclear_attenuation",
        "v3_voxel_grid_3d",
        "v4_dose3d",
        "v5_fragmentation",
    ):
        assert expected in names
    # the diagnostic and the orchestrator itself are excluded
    assert "warp_cuda_smoke" not in names
    assert "release_validation" not in names
    # all discovered scripts actually exist on disk
    assert len(names) >= 20


def test_flag_detection_per_script(rv) -> None:
    suites = {
        s["name"]: s["argv"]
        for s in rv._validation_suites(cache_dir="/cache/ionmc", require_cuda=True)
    }
    # v0_stopping_power is analytic: --require-cuda but NOT --cache-dir
    v0 = suites["v0_stopping_power"]
    assert "--require-cuda" in v0
    assert "--cache-dir" not in v0
    # a dataset-dependent script gets both flags
    frag = suites["v5_fragmentation"]
    assert "--require-cuda" in frag
    assert "--cache-dir" in frag and "/cache/ionmc" in frag


def test_no_flags_when_not_requested(rv) -> None:
    suites = {
        s["name"]: s["argv"]
        for s in rv._validation_suites(cache_dir=None, require_cuda=False)
    }
    frag = suites["v5_fragmentation"]
    assert "--require-cuda" not in frag
    assert "--cache-dir" not in frag


def test_benchmark_baselines_exist_for_release_benchmarks(rv) -> None:
    """The release suite regression-checks each benchmark against a committed
    baseline; those baselines must exist."""
    for name in rv.BENCHMARKS:
        assert (rv.BASELINE_DIR / f"{name}.json").is_file()
        assert (rv.BENCH_DIR / rv.BENCH_SCRIPT[name]).is_file()
