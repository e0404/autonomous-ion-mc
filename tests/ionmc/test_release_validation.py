"""Reference-path tests for the release-validation orchestrator (decision 0033).

These check the suite enumeration and per-script flag detection without running the
full (GPU, dataset-dependent) suite — that runs on the host runner. The orchestrator
itself is a subprocess runner, so its wiring is what these guard.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
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


def test_aggregation_all_pass(rv, monkeypatch) -> None:
    """With every suite returning 0 (and benchmarks' regression 0), the suite is
    release-ready and main() returns 0."""
    monkeypatch.setattr(rv, "_run", lambda argv, timeout: (0, '{"ok": true}', ""))
    monkeypatch.setattr(rv, "_regression_check", lambda b, s, t: (0, ""))
    monkeypatch.setattr("sys.argv", ["release_validation"])

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = rv.main()
    report = json.loads(buf.getvalue())
    assert rc == 0
    assert report["release_ready"] is True
    assert report["failed"] == []
    assert report["n_passed"] == report["n_suites"]


def test_aggregation_one_validation_fails_fails_all(rv, monkeypatch) -> None:
    """A single milestone validation returning 3 makes the whole suite not-ready
    (exit 3), and the failing suite is named with a diagnostic tail."""

    def fake_run(argv, timeout):
        # fail exactly one validation script by name
        if "v1_depth_dose_csda" in " ".join(argv):
            return 3, '{"all_gates_passed": false}', "gate failure detail"
        return 0, '{"ok": true}', ""

    monkeypatch.setattr(rv, "_run", fake_run)
    monkeypatch.setattr(rv, "_regression_check", lambda b, s, t: (0, ""))
    monkeypatch.setattr("sys.argv", ["release_validation"])

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = rv.main()
    report = json.loads(buf.getvalue())
    assert rc == 3
    assert report["release_ready"] is False
    assert "v1_depth_dose_csda" in report["failed"]
    entry = next(r for r in report["results"] if r["name"] == "v1_depth_dose_csda")
    assert entry["diagnostic"]  # a tail was captured for diagnosis


def test_benchmark_requires_gate_and_regression(rv, monkeypatch) -> None:
    """A benchmark passes its physics gate (rc 0) but the regression check fails:
    the benchmark suite must count as failed (fail closed)."""
    monkeypatch.setattr(rv, "_run", lambda argv, timeout: (0, '{"report": 1}', ""))
    monkeypatch.setattr(
        rv, "_regression_check", lambda b, s, t: (3, "physics regression")
    )
    monkeypatch.setattr("sys.argv", ["release_validation", "--only", "bench:"])

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = rv.main()
    report = json.loads(buf.getvalue())
    assert rc == 3
    for entry in report["results"]:
        assert entry["kind"] == "benchmark"
        assert entry["returncode"] == 0
        assert entry["regression_passed"] is False
        assert entry["passed"] is False


def test_benchmark_empty_report_fails_closed(rv, monkeypatch) -> None:
    """A benchmark that exits 0 but produces no stdout report must fail the
    regression check (not be silently skipped)."""
    rc, err = rv._regression_check("dose3d", "   ", timeout=5)
    assert rc != 0
    assert "no report" in err
