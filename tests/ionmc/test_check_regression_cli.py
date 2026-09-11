"""Exit-code contract tests for the regression-check CLI (decision 0032).

The CLI's exit code is the machine-readable V6 gate, so its branches are tested
directly: 0 on a physics match, 3 on a physics regression, 4 on a bad/missing input,
and the ``--emit-baseline`` round trip. Also checks that the *committed* baselines
under ``benchmarks/baselines/`` are well-formed and self-consistent. Reference path
only (no Warp/GPU).
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from ionmc import benchmarking as bm

_REPO = Path(__file__).resolve().parents[2]
_CLI = _REPO / "benchmarks" / "check_regression.py"
_BASELINES = _REPO / "benchmarks" / "baselines"


def _load_cli():
    spec = importlib.util.spec_from_file_location("check_regression", _CLI)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def cli():
    return _load_cli()


def _report(digest="abc123", integral=150.0, peak=None):
    r = {
        "benchmark": "dose3d",
        "config": {"energy_mev": 150.0},
        "provenance": {"machine": "x86_64", "git_sha": "deadbeef"},
        "digest": {
            "reference": {
                "integral_per_history_mev": integral,
                "shape_digest": digest,
            }
        },
    }
    if peak is not None:
        r["peak_throughput_per_s"] = peak
    return r


def _write(tmp_path, name, obj) -> str:
    p = tmp_path / name
    p.write_text(json.dumps(obj), encoding="utf-8")
    return str(p)


def _run(cli, monkeypatch, argv):
    monkeypatch.setattr("sys.argv", ["check_regression", *argv])
    return cli.main()


def test_exit0_on_physics_match(cli, tmp_path, monkeypatch) -> None:
    base = _write(
        tmp_path, "b.json", bm.make_baseline(_report(peak={"warp:cuda:0": 900.0}))
    )
    rep = _write(tmp_path, "r.json", _report(peak={"warp:cuda:0": 1200.0}))
    assert _run(cli, monkeypatch, ["--report", rep, "--baseline", base]) == 0


def test_exit3_on_physics_regression(cli, tmp_path, monkeypatch) -> None:
    base = _write(tmp_path, "b.json", bm.make_baseline(_report(digest="abc123")))
    rep = _write(tmp_path, "r.json", _report(digest="CHANGED"))
    assert _run(cli, monkeypatch, ["--report", rep, "--baseline", base]) == 3


def test_exit4_on_missing_report_file(cli, tmp_path, monkeypatch) -> None:
    base = _write(tmp_path, "b.json", bm.make_baseline(_report()))
    missing = str(tmp_path / "nope.json")
    assert _run(cli, monkeypatch, ["--report", missing, "--baseline", base]) == 4


def test_exit4_on_error_report_without_digest(cli, tmp_path, monkeypatch) -> None:
    """A driver's dataset-uncached error report is valid JSON but not a benchmark
    report; the CLI must reject it cleanly with exit 4, not crash."""
    err = _write(
        tmp_path, "err.json", {"schema_version": 1, "error": "dataset X not cached"}
    )
    base = _write(tmp_path, "b.json", bm.make_baseline(_report()))
    assert _run(cli, monkeypatch, ["--report", err, "--baseline", base]) == 4


def test_exit4_when_baseline_missing_flag(cli, tmp_path, monkeypatch) -> None:
    rep = _write(tmp_path, "r.json", _report())
    assert _run(cli, monkeypatch, ["--report", rep]) == 4


def test_exit4_on_malformed_baseline(cli, tmp_path, monkeypatch) -> None:
    rep = _write(tmp_path, "r.json", _report())
    bad = _write(tmp_path, "b.json", {"not": "a baseline"})
    assert _run(cli, monkeypatch, ["--report", rep, "--baseline", bad]) == 4


def test_emit_baseline_roundtrip(cli, tmp_path, monkeypatch, capsys) -> None:
    rep = _write(tmp_path, "r.json", _report(peak={"warp:cuda:0": 900.0}))
    assert _run(cli, monkeypatch, ["--report", rep, "--emit-baseline"]) == 0
    baseline = json.loads(capsys.readouterr().out)
    assert baseline["physics"]["benchmark"] == "dose3d"
    assert baseline["throughput_by_backend"] == {"warp:cuda:0": 900.0}
    assert baseline["recorded_on"]["git_sha"] == "deadbeef"


@pytest.mark.parametrize("name", ["depth_dose_csda", "dose3d"])
def test_committed_baselines_are_well_formed_and_self_consistent(name) -> None:
    """Each committed baseline is well-formed, and a report reconstructed from its own
    fingerprint checks green — guarding against baseline drift/corruption."""
    baseline = json.loads((_BASELINES / f"{name}.json").read_text(encoding="utf-8"))
    phys = baseline["physics"]
    assert phys["benchmark"] == name
    assert phys["reference_digest"]["shape_digest"]
    assert baseline["throughput_by_backend"]
    assert "recorded_on" in baseline
    # a run identical to the baseline's own physics must pass the gate
    synthetic = {
        "benchmark": phys["benchmark"],
        "config": phys["config"],
        "digest": {"reference": phys["reference_digest"]},
        "peak_throughput_per_s": baseline["throughput_by_backend"],
    }
    result = bm.compare_to_baseline(synthetic, baseline)
    assert result["physics_ok"] is True
    assert result["reference_digest_match"] is True
