"""Regression tests for the fail-closed behaviour of the measurement runner scripts."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "validation" / "scripts" / "warp-architecture"


def run_runner(tmp_path: Path, **env: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", str(SCRIPTS / "run_all.sh"), str(tmp_path / "out")],
        env={"PATH": "/usr/bin:/bin", "PYTHON": sys.executable, **env},
        text=True,
        capture_output=True,
        check=False,
    )


@pytest.mark.parametrize("variable,value", [("REPEATS", "0"), ("REPEATS", "x"), ("STEP_TIMEOUT", "0"), ("STEP_TIMEOUT", "-5")])
def test_runner_rejects_invalid_controls(tmp_path: Path, variable: str, value: str) -> None:
    result = run_runner(tmp_path, **{variable: value})
    assert result.returncode != 0
    assert "must be a positive integer" in result.stderr
    assert not (tmp_path / "out" / "environment.txt").exists()


def test_runner_refuses_existing_results_directory(tmp_path: Path) -> None:
    (tmp_path / "out").mkdir()
    result = run_runner(tmp_path)
    assert result.returncode != 0
    assert "refusing to reuse" in result.stderr


def write_archive(root: Path, sha: str, names: list[str], manifest: list[str]) -> None:
    root.mkdir()
    (root / "environment.txt").write_text(f"git_sha={sha}\n")
    (root / "manifest.txt").write_text("\n".join(manifest) + "\n")
    for name in names:
        (root / f"{name}.txt").write_text(f"# command: x\n# git_sha: {sha}\n# started_utc: t\n# step_timeout_s: 1\nresult\n# exit=0\n")


def summarize(root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPTS / "summarize.py"), str(root)], text=True, capture_output=True, check=False)


def test_summarizer_requires_complete_manifest(tmp_path: Path) -> None:
    write_archive(tmp_path / "ok", "a" * 40, ["01-r1-step", "08-other"], ["01-r1-step", "08-other"])
    assert summarize(tmp_path / "ok").returncode == 0
    write_archive(tmp_path / "missing", "a" * 40, ["01-r1-step"], ["01-r1-step", "08-other"])
    result = summarize(tmp_path / "missing")
    assert result.returncode != 0 and "missing=['08-other']" in result.stderr


def test_summarizer_rejects_foreign_sha(tmp_path: Path) -> None:
    root = tmp_path / "mixed"
    write_archive(root, "a" * 40, ["01-r1-step"], ["01-r1-step"])
    (root / "01-r1-step.txt").write_text("# command: x\n# git_sha: " + "b" * 40 + "\n# started_utc: t\n# step_timeout_s: 1\n# exit=0\n")
    result = summarize(root)
    assert result.returncode != 0 and "does not match environment.txt" in result.stderr
