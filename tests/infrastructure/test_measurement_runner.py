"""Regression tests for the fail-closed behaviour of the measurement runner scripts."""

from __future__ import annotations

import re
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


def write_archive(root: Path, sha: str, repeats: int, names: list[str], manifest: list[str]) -> None:
    root.mkdir()
    (root / "environment.txt").write_text(f"git_sha={sha}\nrepeats={repeats}\n")
    (root / "manifest.txt").write_text("\n".join(manifest) + "\n")
    for name in names:
        (root / f"{name}.txt").write_text(f"# command: x\n# git_sha: {sha}\n# started_utc: t\n# step_timeout_s: 1\nresult\n# exit=0\n")


def summarize(root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPTS / "summarize.py"), str(root)], text=True, capture_output=True, check=False)


def load_summarizer():
    import importlib.util

    spec = importlib.util.spec_from_file_location("summarize", SCRIPTS / "summarize.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_suite_definition_matches_runner() -> None:
    """The summarizer's fixed suite must name exactly the steps run_all.sh executes."""
    module = load_summarizer()
    text = (SCRIPTS / "run_all.sh").read_text()
    runner_repeated = set(re.findall(r"^run_repeated (\d\d) (\S+)", text, re.M))
    runner_single = set(re.findall(r"^run (\d\d) (\S+)", text, re.M))
    assert {f"{n}-{s}" for n, s in runner_repeated} == set(module.REPEATED_STEPS)
    assert {f"{n}-{s}" for n, s in runner_single} == set(module.SINGLE_STEPS)


def complete_names(repeats: int) -> list[str]:
    return sorted(load_summarizer().expected_names(repeats))


def test_summarizer_accepts_complete_archive(tmp_path: Path) -> None:
    names = complete_names(1)
    write_archive(tmp_path / "ok", "a" * 40, 1, names, names)
    result = summarize(tmp_path / "ok")
    assert result.returncode == 0, result.stderr
    assert "## Statistics" in result.stdout


@pytest.mark.parametrize("defect", ["missing_file", "truncated_manifest", "duplicate_manifest", "empty_manifest", "wrong_repeats"])
def test_summarizer_rejects_incomplete_archive(tmp_path: Path, defect: str) -> None:
    names = complete_names(2)
    files, manifest, repeats = list(names), list(names), 2
    if defect == "missing_file":
        files = files[:-1]
    if defect == "truncated_manifest":
        files = files[:-1]
        manifest = manifest[:-1]
    if defect == "duplicate_manifest":
        manifest = manifest + [manifest[0]]
    if defect == "empty_manifest":
        manifest = []
    if defect == "wrong_repeats":
        repeats = 3
    write_archive(tmp_path / defect, "a" * 40, repeats, files, manifest)
    result = summarize(tmp_path / defect)
    assert result.returncode != 0, defect
    assert "not match" in result.stderr or "empty or contains duplicates" in result.stderr


def test_summarizer_rejects_foreign_sha(tmp_path: Path) -> None:
    names = complete_names(1)
    root = tmp_path / "mixed"
    write_archive(root, "a" * 40, 1, names, names)
    (root / f"{names[0]}.txt").write_text("# command: x\n# git_sha: " + "b" * 40 + "\n# started_utc: t\n# step_timeout_s: 1\n# exit=0\n")
    result = summarize(root)
    assert result.returncode != 0 and "does not match environment.txt" in result.stderr
