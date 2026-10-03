"""Fail-closed behaviour of the LV/HR runner and its summariser (no transport is run)."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"
SHA = "a" * 40


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPTS / "run_suite.py"), *args],
        capture_output=True,
        text=True,
        cwd=REPO,
        timeout=300,
    )


def _head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=REPO, check=True
    ).stdout.strip()


def _load_summarize():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("summarize", SCRIPTS / "summarize.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_runner_refuses_existing_output_directory(tmp_path: Path) -> None:
    out = REPO / "validation" / "generated" / "transport" / "test-existing"
    out.mkdir(parents=True, exist_ok=True)
    try:
        r = _run("--suite", "lv", "--out", str(out), "--expected-sha", _head(), "--workers", "2")
        assert r.returncode != 0 and "refusing to reuse" in r.stderr
        assert list(out.iterdir()) == []  # nothing was written
    finally:
        out.rmdir()


def test_runner_requires_generated_directory_and_matching_sha(tmp_path: Path) -> None:
    r = _run("--suite", "lv", "--out", str(tmp_path / "x"), "--expected-sha", _head())
    assert r.returncode != 0 and "validation/generated" in r.stderr
    out = REPO / "validation" / "generated" / "transport" / "test-sha"
    assert not out.exists()
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", SHA)
    assert r.returncode != 0 and "differs from git rev-parse HEAD" in r.stderr
    assert not out.exists()
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", "abc")
    assert r.returncode != 0 and "40 lowercase hex" in r.stderr and not out.exists()
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", _head(), "--workers", "1")
    assert r.returncode != 0 and not out.exists()
    r = _run("--suite", "lv", "--out", str(out))
    assert r.returncode != 0  # --expected-sha is mandatory


def _archive(d: Path, *, exit_code: int = 0, doc: dict | None = None, sha: str = SHA) -> None:  # type: ignore[type-arg]
    d.mkdir()
    (d / "environment.txt").write_text(f"git_sha={SHA}\n")
    (d / "manifest.txt").write_text("01-a\n")
    body = f"# command: x\n# git_sha: {sha}\n# started_utc: now\n# step_timeout_s: 1\n"
    if doc is not None:
        body += "#JSON-BEGIN\n" + json.dumps(doc) + "\n#JSON-END\n"
    (d / "01-a.txt").write_text(body + f"\n# exit={exit_code}\n")


def test_summarize_verdicts(tmp_path: Path) -> None:
    summ = _load_summarize()
    good = tmp_path / "good"
    _archive(good, doc={"pass": True, "reduced": False})
    assert summ.main([str(good), "--expected-sha", SHA]) == 0
    s = json.loads((good / "summary.json").read_text())
    assert s["pass"] and s["conformant"]

    reduced = tmp_path / "reduced"
    _archive(reduced, doc={"pass": True, "reduced": True, "histories": 5, "frozen_histories": 9})
    assert summ.main([str(reduced), "--expected-sha", SHA]) == 0
    assert not json.loads((reduced / "summary.json").read_text())["conformant"]

    for name, kw in {
        "failed": {"exit_code": 1, "doc": {"pass": True}},
        "failed_verdict": {"doc": {"pass": False}},
        "other_sha": {"sha": "b" * 40},
    }.items():
        d = tmp_path / name
        _archive(d, **kw)
        assert summ.main([str(d), "--expected-sha", SHA]) == 1, name

    extra = tmp_path / "extra"
    _archive(extra)
    (extra / "02-b.txt").write_text("stray")
    assert summ.main([str(extra), "--expected-sha", SHA]) == 1
    missing = tmp_path / "missing"
    _archive(missing)
    (missing / "01-a.txt").unlink()
    assert summ.main([str(missing), "--expected-sha", SHA]) == 1
    trailer = tmp_path / "trailer"
    _archive(trailer)
    t = trailer / "01-a.txt"
    t.write_text(t.read_text().replace("# exit=0", ""))
    assert summ.main([str(trailer), "--expected-sha", SHA]) == 1


@pytest.mark.parametrize("suite", ["lv", "hr"])
def test_suite_manifests_are_fixed(suite: str) -> None:
    spec = importlib.util.spec_from_file_location("run_suite", SCRIPTS / "run_suite.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    steps = mod.suite_steps(suite, 4, 1.0)
    names = [s[0] for s in steps]
    assert names == sorted(names) and len(set(names)) == len(names)
    if suite == "hr":
        assert all(s[2].get("IONMC_REQUIRE_CUDA") == "1" for s in steps)
        assert "cuda32" in " ".join(" ".join(s[1]) for s in steps)
    else:
        assert any("t14" in n for n in names) and any("t1-" in n for n in names)
