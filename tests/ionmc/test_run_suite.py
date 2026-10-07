"""Fail-closed behaviour of the LV/HR runner and its summariser (subset labelling, combining,
source attestation, dirty trees); steps are only run for the git/snapshot end-to-end checks."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"
SHA = "a" * 40
WORKERS = "1" if os.environ.get("IONMC_SINGLE_PROCESS") == "1" else "2"


def _load(name: str) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run(*args: str, root: Path = REPO) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(root / "validation/scripts/transport/run_suite.py"), *args],
        capture_output=True,
        text=True,
        cwd=root,
        timeout=600,
        env=_git_env(root),
    )


def _head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=REPO, check=True
    ).stdout.strip()


def test_runner_refuses_existing_output_directory() -> None:
    out = REPO / "validation" / "generated" / "transport" / "test-existing"
    out.mkdir(parents=True, exist_ok=True)
    try:
        r = _run(
            "--suite", "lv", "--out", str(out), "--expected-sha", _head(), "--workers", WORKERS
        )
        assert r.returncode != 0 and "refusing to reuse" in r.stderr
        assert list(out.iterdir()) == []
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
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", _head(), "--workers", "0")
    assert r.returncode != 0 and "--workers >= 1" in r.stderr and not out.exists()
    assert _run("--suite", "lv", "--out", str(out)).returncode != 0  # SHA is mandatory


def test_only_selects_steps_or_fails_before_creating_anything() -> None:
    out = REPO / "validation" / "generated" / "transport" / "test-only"
    assert not out.exists()
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", _head(), "--workers", WORKERS,
             "--only", "99")  # fmt: skip
    assert r.returncode != 0 and "does not select" in r.stderr and not out.exists()


@pytest.mark.parametrize("suite", ["lv", "hr"])
def test_suite_manifests_are_fixed(suite: str) -> None:
    mod = _load("run_suite")
    steps = mod.suite_steps(suite, 4, 1.0)
    names = [s[0] for s in steps]
    assert names == sorted(names) and len(set(names)) == len(names)
    assert names == mod.full_step_names(suite, 2)
    text = " ".join(" ".join(s[1]) for s in steps)
    for sub in ("t12-python-sample", "t12-accelerated-samples", "t12-compare"):
        assert sub in text
    py = [n for n in names if "t12-python-sample" in n]
    assert len(py) == 2 and py[0].endswith("1of2") and py[1].endswith("2of2")
    three = [n for n in mod.full_step_names(suite, 3) if "t12-python-sample" in n]
    assert [n.split("-")[-1] for n in three] == ["1of3", "2of3", "3of3"]
    # the comparison comes after all sample steps
    assert names.index([n for n in names if "t12-compare" in n][0]) > max(
        names.index(n) for n in names if "sample" in n
    )
    if suite == "hr":
        assert all(s[2].get("IONMC_REQUIRE_CUDA") == "1" for s in steps if "cuda" in s[0])
        assert "cuda32" in text
    else:
        assert any("t14" in n for n in names) and any("t1-" in n for n in names)


# -- archives built by hand ------------------------------------------------------------------
def _env(*, dirty: str = "no", source: str = "git", hashes: dict[str, str] | None = None) -> str:
    lines = [
        "suite=lv", f"git_sha={SHA}", f"sha_source={source}", f"tree_dirty={dirty}",
        "scale=1.0", "seed_base=20391004", "python_parts=2", "only=", "source_hashes:",
    ]  # fmt: skip
    lines += [f"  {h}  {p}" for p, h in (hashes or {"src/ionmc/a.py": "0" * 64}).items()]
    return "\n".join(lines) + "\n"


def _archive(d: Path, names: list[str], *, doc: dict | None = None, exit_code: int = 0,  # type: ignore[type-arg]
             env: str | None = None, sha: str = SHA, identity: dict | None = None) -> None:  # type: ignore[type-arg]  # fmt: skip
    """A hand-made archive: every non-pytest step prints a result document naming itself, the
    suite and the SHA (``identity`` overrides those fields; ``doc=None`` prints no document)."""
    summ = _load("summarize")
    d.mkdir()
    (d / "environment.txt").write_text(env if env is not None else _env())
    (d / "manifest.txt").write_text("".join(f"{n}\n" for n in names))
    for n in names:
        body = f"# command: x\n# git_sha: {sha}\n# started_utc: now\n# step_timeout_s: 1\n"
        tag = summ.expected_tag(n)
        if doc is not None and tag is not None:
            full = {**doc, "step": tag, "suite": "lv", "git_sha": SHA, "seed_base": 20391004,
                    **(identity or {})}  # fmt: skip
            body += "#JSON-BEGIN\n" + json.dumps(full) + "\n#JSON-END\n"
        (d / f"{n}.txt").write_text(body + f"\n# exit={exit_code}\n")


def _full() -> list[str]:
    return _load("run_suite").full_step_names("lv", 2)


def test_summarize_verdicts_subset_and_source(tmp_path: Path) -> None:
    summ = _load("summarize")
    full = _full()
    good = tmp_path / "good"
    _archive(good, full, doc={"pass": True, "reduced": False})
    assert summ.main([str(good), "--expected-sha", SHA]) == 0
    s = json.loads((good / "summary.json").read_text())
    assert s["pass"] and not s["subset"] and s["conformant"] and s["source_ok"]

    subset = tmp_path / "subset"  # --only: never conformant
    _archive(subset, full[:2], doc={"pass": True})
    assert summ.main([str(subset), "--expected-sha", SHA]) == 0
    s = json.loads((subset / "summary.json").read_text())
    assert s["subset"] and not s["conformant"]

    reduced = tmp_path / "reduced"
    _archive(reduced, full, doc={"pass": True, "reduced": True})
    assert summ.main([str(reduced), "--expected-sha", SHA]) == 0
    assert not json.loads((reduced / "summary.json").read_text())["conformant"]

    dirty = tmp_path / "dirty"  # unknown dirty state without attestation: not conformant
    _archive(dirty, full, doc={"pass": True}, env=_env(dirty="unknown", source="declared"))
    assert summ.main([str(dirty), "--expected-sha", SHA]) == 0
    s = json.loads((dirty / "summary.json").read_text())
    assert s["pass"] and not s["source_ok"] and not s["conformant"]

    for name, kw in {
        "no_document": {"doc": None},
        "wrong_step": {"doc": {"pass": True}, "identity": {"step": "t2"}},
        "wrong_suite": {"doc": {"pass": True}, "identity": {"suite": "hr"}},
        "wrong_sha_in_document": {"doc": {"pass": True}, "identity": {"git_sha": "e" * 40}},
        "failed": {"exit_code": 1, "doc": {"pass": True}},
        "failed_verdict": {"doc": {"pass": False}},
        "other_sha": {"sha": "b" * 40},
    }.items():
        d = tmp_path / name
        _archive(d, full, **kw)
        assert summ.main([str(d), "--expected-sha", SHA]) == 1, name
    extra = tmp_path / "extra"
    _archive(extra, full)
    (extra / "99-stray.txt").write_text("stray")
    assert summ.main([str(extra), "--expected-sha", SHA]) == 1
    missing = tmp_path / "missing"
    _archive(missing, full)
    (missing / f"{full[0]}.txt").unlink()
    assert summ.main([str(missing), "--expected-sha", SHA]) == 1
    trailer = tmp_path / "trailer"
    _archive(trailer, full)
    t = trailer / f"{full[0]}.txt"
    t.write_text(t.read_text().replace("# exit=0", ""))
    assert summ.main([str(trailer), "--expected-sha", SHA]) == 1


def test_combine_requires_exactly_the_full_suite(tmp_path: Path) -> None:
    summ = _load("summarize")
    full = _full()
    a, b = full[:6], full[6:]
    _archive(tmp_path / "a", a, doc={"pass": True})
    _archive(tmp_path / "b", b, doc={"pass": True})
    out = tmp_path / "combined.json"
    assert summ.main(["--combine", str(tmp_path / "a"), str(tmp_path / "b"), "--expected-sha", SHA,
                      "--out", str(out)]) == 0  # fmt: skip
    s = json.loads(out.read_text())
    assert s["pass"] and s["complete"] and s["conformant"] and not s["missing_steps"]

    # a missing step: not complete, not conformant (and the verdict fails)
    _archive(tmp_path / "c", b[:-1], doc={"pass": True})
    out2 = tmp_path / "combined2.json"
    assert summ.main(["--combine", str(tmp_path / "a"), str(tmp_path / "c"), "--expected-sha", SHA,
                      "--out", str(out2)]) == 1  # fmt: skip
    s2 = json.loads(out2.read_text())
    assert not s2["pass"] and not s2["conformant"] and s2["missing_steps"] == [b[-1]]

    # a duplicated step is an error
    _archive(tmp_path / "d", b + [a[0]], doc={"pass": True})
    out3 = tmp_path / "combined3.json"
    assert summ.main(["--combine", str(tmp_path / "a"), str(tmp_path / "d"), "--expected-sha", SHA,
                      "--out", str(out3)]) == 1  # fmt: skip
    assert any("appears in" in p for p in json.loads(out3.read_text())["problems"])

    # different source hashes: the archives are not parts of one run
    _archive(tmp_path / "e", b, doc={"pass": True}, env=_env(hashes={"src/ionmc/a.py": "1" * 64}))
    out4 = tmp_path / "combined4.json"
    assert summ.main(["--combine", str(tmp_path / "a"), str(tmp_path / "e"), "--expected-sha", SHA,
                      "--out", str(out4)]) == 1  # fmt: skip
    # an existing output file is never overwritten
    with pytest.raises(SystemExit):
        summ.main(["--combine", str(tmp_path / "a"), "--expected-sha", SHA, "--out", str(out)])


# -- git and snapshot end to end -----------------------------------------------------------------
def _git_env(root: Path) -> dict[str, str]:
    """Environment for git in sandboxes without a home directory or identity."""
    return dict(os.environ, HOME=str(root.parent), GIT_CONFIG_GLOBAL="/dev/null",
                GIT_CONFIG_NOSYSTEM="1")  # fmt: skip


def _make_repo(root: Path) -> str:
    """A throw-away git repository (in tmp) holding the scripts and sources under test."""
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc", "generated")
    for rel in ("src/ionmc", "tests/ionmc", "validation/scripts/transport", "benchmarks/transport"):
        shutil.copytree(REPO / rel, root / rel, ignore=ignore)
    for rel in (
        "pyproject.toml",
        "uv.lock",
        "validation/plans/v3-003-acceptance.md",
        "validation/plans/v3-004-acceptance.md",
        "tests/data/synthetic_lookup.json",
    ):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(REPO / rel, root / rel)
    (root / ".gitignore").write_text("validation/generated/\nbenchmarks/generated/\n")
    env = _git_env(root)
    base = ["git", "-c", "user.name=t", "-c", "user.email=t@t", "-C", str(root)]
    try:
        for cmd in (["init", "-q"], ["add", "-A"], ["commit", "-q", "-m", "snapshot"]):
            subprocess.run([*base, *cmd], check=True, env=env, capture_output=True)
        return subprocess.run([*base, "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True, env=env).stdout.strip()  # fmt: skip
    except (subprocess.CalledProcessError, OSError) as exc:
        pytest.skip(f"cannot create a throw-away git repository in this sandbox: {exc}")


def test_dirty_tree_is_refused_and_clean_tree_runs(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    sha = _make_repo(repo)
    out = repo / "validation" / "generated" / "transport" / "run"
    args = ("--suite", "lv", "--out", str(out), "--expected-sha", sha, "--workers", WORKERS,
            "--only", "03")  # fmt: skip
    r = _run(*args, root=repo)
    assert r.returncode == 0, r.stderr[-1500:] + r.stdout[-1500:]
    s = json.loads((out / "summary.json").read_text())
    assert s["tree_dirty"] == "no" and s["source_ok"] and s["subset"] and not s["conformant"]
    env = (out / "environment.txt").read_text()
    assert "source_hashes:" in env and "src/ionmc/config.py" in env
    # dirty: refused before anything is written
    (repo / "src/ionmc/config.py").write_text(
        (repo / "src/ionmc/config.py").read_text() + "\n# x\n"
    )
    out2 = repo / "validation" / "generated" / "transport" / "run2"
    r = _run(*args[:3], str(out2), *args[4:], root=repo)
    assert r.returncode != 0 and "dirty working tree" in r.stderr and not out2.exists()


def test_snapshot_without_git_is_attested_later(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    sha = _make_repo(repo)
    snap = tmp_path / "snap"
    shutil.copytree(repo, snap, ignore=shutil.ignore_patterns(".git", "generated"))
    out = snap / "validation" / "generated" / "transport" / "run"
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", sha, "--workers", WORKERS,
             "--only", "03", root=snap)  # fmt: skip
    assert r.returncode == 0, r.stderr[-1500:] + r.stdout[-1500:]
    s = json.loads((out / "summary.json").read_text())
    assert s["sha_source"].startswith("declared") and s["tree_dirty"] == "unknown"
    assert not s["source_ok"]
    # attest where git exists: the archive's hashes equal the blobs of the commit
    shutil.copytree(out, tmp_path / "archive")
    summ = subprocess.run(
        [sys.executable, str(repo / "validation/scripts/transport/summarize.py"),
         str(tmp_path / "archive"), "--expected-sha", sha, "--attest-sha", sha],
        capture_output=True, text=True, env=_git_env(repo),
    )  # fmt: skip
    assert summ.returncode == 0, summ.stderr
    a = json.loads((tmp_path / "archive" / "summary.json").read_text())
    assert a["attestation"]["valid"] and a["source_ok"]
    # a tampered source in the snapshot is detected
    snap2 = tmp_path / "snap2"
    shutil.copytree(repo, snap2, ignore=shutil.ignore_patterns(".git", "generated"))
    (snap2 / "src/ionmc/config.py").write_text((snap2 / "src/ionmc/config.py").read_text() + "#t\n")
    out2 = snap2 / "validation" / "generated" / "transport" / "run"
    r = _run("--suite", "lv", "--out", str(out2), "--expected-sha", sha, "--workers", WORKERS,
             "--only", "03", root=snap2)  # fmt: skip
    assert r.returncode == 0
    shutil.copytree(out2, tmp_path / "archive2")
    subprocess.run(
        [sys.executable, str(repo / "validation/scripts/transport/summarize.py"),
         str(tmp_path / "archive2"), "--expected-sha", sha, "--attest-sha", sha],
        capture_output=True, text=True,
    )  # fmt: skip
    b = json.loads((tmp_path / "archive2" / "summary.json").read_text())
    assert not b["attestation"]["valid"] and not b["source_ok"]
    assert any("content differs: src/ionmc/config.py" in m for m in b["attestation"]["mismatches"])

    # every execution-defining file is covered: an altered test file, fixture, lock file or
    # plan also fails the attestation
    for rel in ("tests/ionmc/conftest.py", "uv.lock", "validation/plans/v3-003-acceptance.md",
                "pyproject.toml"):  # fmt: skip
        snap3 = tmp_path / f"snap-{Path(rel).name}"
        shutil.copytree(repo, snap3, ignore=shutil.ignore_patterns(".git", "generated"))
        (snap3 / rel).write_text((snap3 / rel).read_text() + "\n#t\n")
        out3 = snap3 / "validation" / "generated" / "transport" / "run"
        r = _run("--suite", "lv", "--out", str(out3), "--expected-sha", sha, "--workers", WORKERS,
                 "--only", "03", root=snap3)  # fmt: skip
        assert r.returncode == 0, r.stderr[-800:]
        arch = tmp_path / f"arch-{Path(rel).name}"
        shutil.copytree(out3, arch)
        subprocess.run([sys.executable, str(repo / "validation/scripts/transport/summarize.py"),
                        str(arch), "--expected-sha", sha, "--attest-sha", sha],
                       capture_output=True, text=True, env=_git_env(repo))  # fmt: skip
        c = json.loads((arch / "summary.json").read_text())
        assert not c["attestation"]["valid"], rel
        assert any(f"content differs: {rel}" in m for m in c["attestation"]["mismatches"]), rel


def test_t12_split_steps_roundtrip_and_tamper_detection(tmp_path: Path) -> None:
    """The python sample in two history ranges, the accelerated samples and the comparison, each
    a separate step with its own output. The comparison trusts nothing it is given: it rebuilds the
    expected configuration itself, verifies every producer archive completely, recomputes the
    source hashes of the executing tree, and refuses tampered files or step archives, tampered
    metadata, a wrong seed or SHA, other source hashes (a stored attestation is ignored) and
    missing parts. (Tiny histories: the verdict itself is not asserted, only the plumbing.)"""
    import hashlib

    import numpy as np

    sha = "c" * 40
    run_suite = _load("run_suite")
    current = {
        str(f.relative_to(run_suite.REPO)): run_suite.sha256(f) for f in run_suite.source_files()
    }
    env = dict(os.environ, IONMC_RUN_SHA=sha, IONMC_RUN_SUITE="lv", PYTHONPATH=str(REPO / "src"))
    base = ["--energy", "70", "--lateral-bin", "0.5", "--half-width", "8", "--scale", "0.01"]
    steps = str(SCRIPTS / "steps.py")
    own, prod = tmp_path / "own", tmp_path / "prod"
    names = ["01-t12-python-sample-1of2", "02-t12-python-sample-2of2", "03-t12-accelerated-samples"]
    for d in (own, prod):
        d.mkdir()
        (d / "environment.txt").write_text(_env(hashes=current).replace(SHA, sha))
    (prod / "manifest.txt").write_text("".join(f"{n}\n" for n in names))

    def step(
        name: str, *args: str, environ: dict | None = None
    ) -> subprocess.CompletedProcess[str]:  # type: ignore[type-arg]
        return subprocess.run([sys.executable, steps, name, *base, *args], capture_output=True,
                              text=True, env=environ or env, timeout=900)  # fmt: skip

    def archive(fname: str, proc: subprocess.CompletedProcess[str]) -> None:
        body = f"# command: x\n# git_sha: {sha}\n# started_utc: n\n# step_timeout_s: 1\n"
        (prod / f"{fname}.txt").write_text(body + proc.stdout + f"\n# exit={proc.returncode}\n")

    samples = str(prod / "samples")
    for i in (1, 2):
        p = step(
            "t12-python-sample", "--part", f"{i}/2", "--out-dir", samples, "--workers", WORKERS
        )
        assert p.returncode == 0, p.stdout[-2000:] + p.stderr[-2000:]
        archive(names[i - 1], p)
    p = step("t12-accelerated-samples", "--samples", "cpu32", "--out-dir", samples)
    assert p.returncode == 0, p.stderr[-2000:]
    archive(names[2], p)
    args = ("t12-compare", "--pairs", "python:cpu32", "--dirs", str(own), str(prod))
    c = step(*args)
    doc = json.loads(c.stdout.split("#JSON-BEGIN")[1].split("#JSON-END")[0])
    assert doc["step"] == "t12-compare" and doc["suite"] == "lv" and doc["git_sha"] == sha
    assert doc["reduced"] and doc["samples"]["python"]["parts"] == 2
    assert doc["samples"]["python"]["reduced"] and doc["samples"]["cpu32"]["reduced"]
    assert doc["samples"]["python"]["frozen_histories"] == 4000
    assert doc["samples"]["cpu32"]["frozen_histories"] == 1_000_000
    assert doc["t12"]["python_vs_cpu32"]["permutation"]["n_perm"] == 2000

    def refused(proc: subprocess.CompletedProcess[str], text: str) -> None:
        assert proc.returncode != 0 and text in (proc.stderr + proc.stdout), proc.stderr[-800:]

    refused(step(*args, environ=dict(env, IONMC_RUN_SHA="d" * 40)), "SHA")  # wrong SHA
    refused(step(*args, "--seed", "99"), "differs from expected")  # wrong seed
    refused(step(*args, "--energy", "71"), "differs from expected")  # other configuration
    # a tampered producer step archive (its trailer says the step failed): the archive is refused
    good = (prod / f"{names[0]}.txt").read_text()
    (prod / f"{names[0]}.txt").write_text(good.replace("# exit=0", "# exit=1"))
    refused(step(*args), "does not verify")
    (prod / f"{names[0]}.txt").write_text(good)
    assert "#JSON-BEGIN" in step(*args).stdout
    # other source hashes in the producer's archive are refused - also when a stored attestation
    # claims validity (it is ignored; the comparison recomputes it against the commit)
    other = dict(current, **{next(iter(current)): "1" * 64})
    (prod / "environment.txt").write_text(_env(hashes=other).replace(SHA, sha))
    refused(step(*args), "source hashes differ")
    (prod / "summary.json").write_text(
        json.dumps({"attestation": {"valid": True, "attested_sha": sha}})
    )
    refused(step(*args), "source hashes differ")
    (prod / "summary.json").unlink()
    (prod / "environment.txt").write_text(_env(hashes=current).replace(SHA, sha))
    # this run's own archive must match the tree that is executing
    (own / "environment.txt").write_text(_env(hashes=other).replace(SHA, sha))
    refused(step(*args), "differ from the executing tree")
    (own / "environment.txt").write_text(_env(hashes=current).replace(SHA, sha))
    # tampered metadata with a consistently updated sha256 in the producer's archive
    f = prod / "samples" / "t12-python-part-1-of-2.npz"
    data = dict(np.load(f))
    meta = json.loads(str(data["meta"]))
    meta["energy_mev"] = 71.0
    data["meta"] = np.array(json.dumps(meta, sort_keys=True))
    old = f.read_bytes()
    np.savez(f, **data)
    new_digest = hashlib.sha256(f.read_bytes()).hexdigest()
    arch = prod / f"{names[0]}.txt"
    arch.write_text(arch.read_text().replace(hashlib.sha256(old).hexdigest(), new_digest))
    refused(step(*args), "differs from expected")
    # tampered sample file: sha256 mismatch
    f.write_bytes(f.read_bytes() + b"x")
    refused(step(*args), "sha256 mismatch")
    # a missing part: the sample no longer tiles the history range
    f.unlink()
    assert step(*args).returncode != 0


def test_step_timeout_kills_the_whole_process_group(tmp_path: Path) -> None:
    """A step that spawns a long-sleeping child and exceeds its timeout is archived with exit 124
    and leaves no descendant alive (the step runs in its own process group)."""
    run_suite = _load("run_suite")
    script = (
        "import subprocess, sys, time\n"
        "child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(300)'])\n"
        "print('CHILD', child.pid, flush=True)\n"
        "time.sleep(300)\n"
    )
    out = tmp_path / "step.txt"
    with out.open("w") as fh:
        code = run_suite.run_step(
            [sys.executable, "-c", script], cwd=tmp_path, env=dict(os.environ), stdout=fh,
            timeout=3.0, grace=2.0,
        )  # fmt: skip
    assert code == 124
    pid = int(next(ln for ln in out.read_text().splitlines() if ln.startswith("CHILD")).split()[1])
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)  # no descendant remains


def test_step_that_finishes_leaves_no_orphans_and_keeps_its_exit_code(tmp_path: Path) -> None:
    run_suite = _load("run_suite")
    script = (
        "import subprocess, sys\n"
        "child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(300)'])\n"
        "print('CHILD', child.pid, flush=True)\n"
        "sys.exit(3)\n"
    )
    out = tmp_path / "step.txt"
    with out.open("w") as fh:
        code = run_suite.run_step(
            [sys.executable, "-c", script], cwd=tmp_path, env=dict(os.environ), stdout=fh,
            timeout=30.0, grace=2.0,
        )  # fmt: skip
    assert code == 3
    pid = int(next(ln for ln in out.read_text().splitlines() if ln.startswith("CHILD")).split()[1])
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)


def test_suites_forward_an_inner_timeout_below_the_step_timeout() -> None:
    run_suite = _load("run_suite")
    steps = run_suite.suite_steps("lv", 4, 1.0, step_timeout=1000)
    scripted = [s for s in steps if "steps.py" in " ".join(s[1])]
    assert scripted and all("--timeout" in s[1] for s in scripted)
    for s in scripted:
        assert float(s[1][s[1].index("--timeout") + 1]) < 1000
    assert all("--timeout" not in s[1] for s in steps if "pytest" in s[0])


def test_seed_base_is_forwarded_recorded_and_checked(tmp_path: Path) -> None:
    """``--seed-base`` reaches every scripted step, is part of the archive identity (so subsets made
    with different bases cannot be combined) and every step document must carry the archive's
    base."""
    run_suite, summ = _load("run_suite"), _load("summarize")
    steps = run_suite.suite_steps("lv", 4, 1.0, seed_base=20391004)
    scripted = [s for s in steps if "steps.py" in " ".join(s[1])]
    assert scripted and all(s[1][s[1].index("--seed-base") + 1] == "20391004" for s in scripted)
    full = _full()
    env = _env()
    good = tmp_path / "good"
    _archive(good, full, doc={"pass": True}, env=env, identity={"seed_base": 20391004})
    assert summ.main([str(good), "--expected-sha", SHA]) == 0
    bad = tmp_path / "bad"  # a step that ran with another base than the archive records
    _archive(bad, full, doc={"pass": True}, env=env, identity={"seed_base": 20271004})
    assert summ.main([str(bad), "--expected-sha", SHA]) == 1
    # archives with different seed bases are not parts of one run
    a, b = full[:6], full[6:]
    _archive(tmp_path / "a", a, doc={"pass": True}, env=env, identity={"seed_base": 20391004})
    other = _env().replace("seed_base=20391004", "seed_base=20271004")
    _archive(tmp_path / "b", b, doc={"pass": True}, env=other, identity={"seed_base": 20271004})
    out = tmp_path / "combined.json"
    assert summ.main(["--combine", str(tmp_path / "a"), str(tmp_path / "b"),
                      "--expected-sha", SHA, "--out", str(out)]) == 1  # fmt: skip
    assert any("differ" in p for p in json.loads(out.read_text())["problems"])


def test_only_the_qualification_seed_base_can_be_conformant(tmp_path: Path) -> None:
    """The default base is the qualification base 20391004; an archive made with a consumed base
    (20271004) or without a recorded base verifies but is never conformant, with the reason
    recorded."""
    run_suite, summ = _load("run_suite"), _load("summarize")
    assert run_suite.DEFAULT_SEED_BASE == 20391004 == summ.QUALIFICATION_SEED_BASE
    for suite_args in (run_suite.suite_steps("lv", 4, 1.0), run_suite.suite_steps("hr", 4, 1.0)):
        scripted = [s for s in suite_args if "steps.py" in " ".join(s[1])]
        assert all(s[1][s[1].index("--seed-base") + 1] == "20391004" for s in scripted)
    full = _full()
    ok = tmp_path / "ok"
    _archive(ok, full, doc={"pass": True})
    assert summ.main([str(ok), "--expected-sha", SHA]) == 0
    s = json.loads((ok / "summary.json").read_text())
    assert s["conformant"] and s["non_conformant_reasons"] == []
    rehearsal = tmp_path / "rehearsal"
    _archive(rehearsal, full, doc={"pass": True},
             env=_env().replace("seed_base=20391004", "seed_base=20271004"),
             identity={"seed_base": 20271004})  # fmt: skip
    assert summ.main([str(rehearsal), "--expected-sha", SHA]) == 0  # it verifies ...
    r = json.loads((rehearsal / "summary.json").read_text())
    assert r["pass"] and not r["conformant"]  # ... but cannot qualify
    assert any("rehearsal" in x for x in r["non_conformant_reasons"])
    missing = tmp_path / "missing"
    _archive(missing, full, doc={"pass": True}, env=_env().replace("seed_base=20391004\n", ""))
    assert summ.main([str(missing), "--expected-sha", SHA]) == 0
    m = json.loads((missing / "summary.json").read_text())
    assert not m["conformant"] and "not recorded" in m["non_conformant_reasons"][0]


@pytest.mark.parametrize(
    "base",
    [20261004, 20271004, 20281004, 20291004, 20301004, 20311004, 20321004, 20331004, 20341004],
)
def test_non_qualification_bases_are_recorded_as_such(tmp_path: Path, base: int) -> None:
    """20261004 (rehearsal), 20271004 (T9 investigation), 20281004, 20291004, 20301004, 20311004,
    20321004 and 20331004 (first to sixth qualification attempts) verify but never qualify."""
    summ = _load("summarize")
    d = tmp_path / "x"
    env = _env().replace("seed_base=20391004", f"seed_base={base}")
    _archive(d, _full(), doc={"pass": True}, env=env, identity={"seed_base": base})
    assert summ.main([str(d), "--expected-sha", SHA]) == 0
    s = json.loads((d / "summary.json").read_text())
    assert not s["conformant"] and s["seed_base"] == str(base)
    assert "non-qualification" in s["non_conformant_reasons"][0]


@pytest.mark.parametrize("suite", ["lv", "hr"])
def test_t12_compare_includes_the_float64_control_pair(suite: str) -> None:
    """python versus warp-cpu float64 (T1 bit-identical code paths) is the same-algorithm control
    pair of the T12 comparison in both suites."""
    run_suite = _load("run_suite")
    cmp_args = next(s[1] for s in run_suite.suite_steps(suite, 4, 1.0) if "t12-compare" in s[1])
    pairs = cmp_args[cmp_args.index("--pairs") + 1].split(",")
    assert "python:cpu64" in pairs and "python:cpu32" in pairs and "cpu32:cpu64" in pairs


# -- single-process diagnostic mode ------------------------------------------------------------
DIAG_ENV_LINES = "execution_mode=single-process-diagnostic\nworkers=1\n"
DEFERRED = "multiprocessing-specific check deferred in single-process diagnostic mode"


def _diag_archive(
    d: Path, *, standard: bool = False, defer: tuple[str, ...] = ("04-t13-workers",)
) -> None:
    full = _full()
    env = _env() if standard else _env() + DIAG_ENV_LINES
    _archive(d, full, doc={"pass": True, "reduced": False}, env=env)
    for n in defer:
        (d / f"{n}.txt").write_text(
            f"# command: x\n# git_sha: {SHA}\n# started_utc: now\n# step_timeout_s: 1\n"
            f"# status: deferred\n# reason: {DEFERRED}\n\n# exit=0\n"
        )


def test_single_process_mode_defers_only_the_worker_partition_step(tmp_path: Path) -> None:
    mod = _load("run_suite")
    assert mod.deferred_step_names("lv") == ["04-t13-workers"]
    assert mod.deferred_step_names("hr") == []
    assert mod.SINGLE_PROCESS_ENV["IONMC_SINGLE_PROCESS"] == "1"
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        assert mod.SINGLE_PROCESS_ENV[k] == "1"
    lv = mod.suite_steps("lv", 1, 1.0, single_process=True)
    assert lv[0][1][-2:] == ["-m", "not multiprocess"]
    assert "not multiprocess" not in mod.suite_steps("lv", 4, 1.0)[0][1]
    hr = mod.suite_steps("hr", 1, 1.0, single_process=True)
    assert hr[0][1][-1] == "cuda and not multiprocess"
    # histories and seeds are those of the standard suite; only the worker count differs
    std = mod.suite_steps("lv", 4, 1.0, seed_base=20391004)
    one = mod.suite_steps("lv", 1, 1.0, seed_base=20391004)
    for (n1, c1, _), (n2, c2, _) in zip(std[1:], one[1:], strict=True):
        assert n1 == n2
        norm = [x for i, x in enumerate(c1) if not (c1[i - 1] == "--workers" or x == "--workers")]
        norm2 = [x for i, x in enumerate(c2) if not (c2[i - 1] == "--workers" or x == "--workers")]
        assert norm == norm2, n1


def test_diagnostic_archive_passes_but_is_never_conformant(tmp_path: Path) -> None:
    summ = _load("summarize")
    d = tmp_path / "diag"
    _diag_archive(d)
    assert summ.main([str(d), "--expected-sha", SHA]) == 0
    s = json.loads((d / "summary.json").read_text())
    assert s["pass"] and not s["subset"] and not s["conformant"]
    assert s["deferred_steps"] == ["04-t13-workers"]
    assert s["execution_mode"] == "single-process-diagnostic"
    assert "deferred multiprocessing checks" in s["non_conformant_reasons"]
    assert s["steps"]["04-t13-workers"]["status"] == "deferred"
    # deferral is only legitimate in a diagnostic archive and only for the listed steps
    std = tmp_path / "std"
    _diag_archive(std, standard=True)
    assert summ.main([str(std), "--expected-sha", SHA]) == 1
    other = tmp_path / "other"
    _diag_archive(other, defer=("04-t13-workers", "02-t1-trace-parity-256x150MeV"))
    assert summ.main([str(other), "--expected-sha", SHA]) == 1
    # an executed failing step is still a failure
    bad = tmp_path / "bad"
    _diag_archive(bad)
    (bad / "02-t1-trace-parity-256x150MeV.txt").write_text(
        f"# command: x\n# git_sha: {SHA}\n# started_utc: now\n# step_timeout_s: 1\n\n# exit=1\n"
    )
    assert summ.main([str(bad), "--expected-sha", SHA]) == 1


def test_combine_carries_the_deferred_list(tmp_path: Path) -> None:
    summ = _load("summarize")
    full = _full()
    a, b = tmp_path / "a", tmp_path / "b"
    _archive(a, full[:5], doc={"pass": True}, env=_env() + DIAG_ENV_LINES)
    (a / "04-t13-workers.txt").write_text(
        f"# command: x\n# git_sha: {SHA}\n# started_utc: now\n# step_timeout_s: 1\n"
        f"# status: deferred\n# reason: {DEFERRED}\n\n# exit=0\n"
    )
    _archive(b, full[5:], doc={"pass": True}, env=_env() + DIAG_ENV_LINES)
    out = tmp_path / "combined.json"
    assert summ.main(["--combine", str(a), str(b), "--expected-sha", SHA, "--out", str(out)]) == 0
    c = json.loads(out.read_text())
    assert c["pass"] and c["deferred_steps"] == ["04-t13-workers"] and not c["conformant"]
    assert "deferred multiprocessing checks" in c["non_conformant_reasons"]


def test_runner_single_process_end_to_end(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    sha = _make_repo(repo)
    out = repo / "validation" / "generated" / "transport" / "run"
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", sha, "--workers", "1",
             "--only", "03", "04", root=repo)  # fmt: skip
    assert r.returncode == 0, r.stderr[-1500:] + r.stdout[-1500:]
    s = json.loads((out / "summary.json").read_text())
    assert s["pass"] and s["execution_mode"] == "single-process-diagnostic" and not s["conformant"]
    assert s["deferred_steps"] == ["04-t13-workers"]
    text = (out / "04-t13-workers.txt").read_text()
    assert "# status: deferred" in text and DEFERRED in text
    env = (out / "environment.txt").read_text()
    assert "execution_mode=single-process-diagnostic" in env and "OMP_NUM_THREADS=1" in env
    assert "IONMC_SINGLE_PROCESS=1" in env


# -- V3-004 suites (lv4, hr4) --------------------------------------------------------------------
QUAL4 = 20401004


def _env4(suite: str = "lv4", base: int = QUAL4) -> str:
    return (
        _env()
        .replace("suite=lv", f"suite={suite}")
        .replace("seed_base=20391004", f"seed_base={base}")
    )


def _archive4(d: Path, names: list[str], suite: str = "lv4", base: int = QUAL4,
              env: str | None = None, doc: dict | None = None) -> None:  # type: ignore[type-arg]  # fmt: skip
    summ = _load("summarize")
    d.mkdir()
    (d / "environment.txt").write_text(env if env is not None else _env4(suite, base))
    (d / "manifest.txt").write_text("".join(f"{n}\n" for n in names))
    for n in names:
        body = f"# command: x\n# git_sha: {SHA}\n# started_utc: now\n# step_timeout_s: 1\n"
        tag = summ.expected_tag(n)
        if tag is not None:
            full = {"pass": True, "reduced": False, "step": tag, "suite": suite, "git_sha": SHA,
                    "seed_base": base, **(doc or {})}  # fmt: skip
            body += "#JSON-BEGIN\n" + json.dumps(full) + "\n#JSON-END\n"
        (d / f"{n}.txt").write_text(body + "\n# exit=0\n")


def test_v4_suite_manifests_seed_base_and_hashed_files() -> None:
    mod = _load("run_suite")
    lv, hr = mod.full_step_names("lv4", 2), mod.full_step_names("hr4", 2)
    assert [n.split("-", 1)[1] for n in lv] == [
        "pytest-scoring-warp-cpu", "a16-qualified-path-regression",
        "a11-lv-python-vs-warp-cpu-256x150MeV", "a15-chunks-cpu", "a15-workers",
        "a7-step-independence", "a8-offline-let", "a13-let-profile-exploratory",
        "a9-part-1of2", "a9-part-2of2", "a9-compare",
    ]  # fmt: skip
    assert [n.split("-", 1)[1] for n in hr] == [
        "pytest-cuda-scoring", "a15-chunks-cuda", "a11-hr-channel-parity"
    ]  # fmt: skip
    assert names_sorted(lv) and names_sorted(hr)
    for suite in ("lv4", "hr4"):
        assert mod.DEFAULT_SEED_BASES[suite] == QUAL4
        steps = mod.suite_steps(suite, 4, 0.5)
        scripted = [s for s in steps if "steps_v4.py" in " ".join(s[1])]
        assert scripted
        for s in scripted:
            assert s[1][s[1].index("--seed-base") + 1] == str(QUAL4)
            assert "--timeout" in s[1] and "--scale" in s[1] or "a16" in s[0]
    # the V3-003 suites keep their bases and file sets
    assert mod.DEFAULT_SEED_BASES["lv"] == 20391004 == mod.DEFAULT_SEED_BASE
    assert mod.source_file_list("lv") == mod.SOURCE_FILES
    assert "validation/plans/v3-004-acceptance.md" in mod.source_file_list("lv4")
    assert "tests/data/synthetic_lookup.json" in mod.source_file_list("hr4")
    rel = {str(f.relative_to(mod.REPO)) for f in mod.source_files("lv4")}
    assert {"src/ionmc/let_offline.py", "validation/scripts/transport/steps_v4.py",
            "validation/plans/v3-004-acceptance.md",
            "tests/data/synthetic_lookup.json"} <= rel  # fmt: skip
    assert "validation/plans/v3-004-acceptance.md" not in {
        str(f.relative_to(mod.REPO)) for f in mod.source_files("lv")
    }
    assert all(s[2].get("IONMC_REQUIRE_CUDA") == "1" for s in mod.suite_steps("hr4", 4, 1.0)
               if "cuda" in s[0] or "hr-" in s[0])  # fmt: skip


def names_sorted(names: list[str]) -> bool:
    return names == sorted(names) and len(set(names)) == len(names)


def test_v4_timeouts_deferred_step_and_single_process_markers() -> None:
    mod = _load("run_suite")
    assert mod.deferred_step_names("lv4") == ["05-a15-workers"]
    assert mod.deferred_step_names("hr4") == []
    one = mod.suite_steps("lv4", 1, 1.0, single_process=True)
    assert one[0][1][-2:] == ["-m", "not multiprocess"]
    assert mod.suite_steps("hr4", 1, 1.0, single_process=True)[0][1][-1] == (
        "cuda and not multiprocess"
    )
    assert mod.step_timeout_s("lv4", "09-a9-part-1of2", 1500) == 3300
    assert mod.step_timeout_s("lv4", "09-a9-part-1of2", 5000) == 5000
    assert mod.step_timeout_s("hr4", "03-a11-hr-channel-parity", 1500) == 3600
    assert mod.step_timeout_s("lv4", "04-a15-chunks-cpu", 1500) == 1500
    std = mod.suite_steps("lv4", 4, 1.0)
    for (n1, c1, _), (n2, c2, _) in zip(std[1:], one[1:], strict=True):
        assert n1 == n2

        def strip(c: list[str]) -> list[str]:
            return [x for i, x in enumerate(c) if not (c[i - 1] == "--workers" or x == "--workers")]

        assert strip(c1) == strip(c2) or "a15-workers" in n1 or "a9" in n1, n1  # fmt: skip


def test_v4_summary_conformance_deferred_and_combine(tmp_path: Path) -> None:
    summ, mod = _load("summarize"), _load("run_suite")
    full = mod.full_step_names("lv4", 2)
    ok = tmp_path / "ok"
    _archive4(ok, full)
    assert summ.main([str(ok), "--expected-sha", SHA]) == 0
    s = json.loads((ok / "summary.json").read_text())
    assert s["pass"] and s["conformant"] and not s["subset"], s["non_conformant_reasons"]
    reh = tmp_path / "reh"  # the rehearsal base verifies but never qualifies
    _archive4(reh, full, base=20351004)
    assert summ.main([str(reh), "--expected-sha", SHA]) == 0
    r = json.loads((reh / "summary.json").read_text())
    assert r["pass"] and not r["conformant"]
    assert any("rehearsal" in x and "non-qualification" in x for x in r["non_conformant_reasons"])
    used = tmp_path / "used"  # the consumed base 20361004 (amendment 4) never qualifies
    _archive4(used, full, base=20361004)
    assert summ.main([str(used), "--expected-sha", SHA]) == 0
    u = json.loads((used / "summary.json").read_text())
    assert u["pass"] and not u["conformant"]
    assert any("consumed" in x for x in u["non_conformant_reasons"])
    old = (
        tmp_path / "old-qual"
    )  # 20381004: the V3-004 base of the previous head, consumed (V3-003D)
    _archive4(old, full, base=20381004)
    assert summ.main([str(old), "--expected-sha", SHA]) == 0
    o = json.loads((old / "summary.json").read_text())
    assert o["pass"] and not o["conformant"]
    # single-process diagnostic archive: the workers step may be deferred, never conformant
    d = tmp_path / "diag"
    _archive4(d, full, env=_env4() + DIAG_ENV_LINES)
    (d / "05-a15-workers.txt").write_text(
        f"# command: x\n# git_sha: {SHA}\n# started_utc: now\n# step_timeout_s: 1\n"
        f"# status: deferred\n# reason: {DEFERRED}\n\n# exit=0\n"
    )
    assert summ.main([str(d), "--expected-sha", SHA]) == 0
    sd = json.loads((d / "summary.json").read_text())
    assert sd["deferred_steps"] == ["05-a15-workers"] and not sd["conformant"]
    bad = tmp_path / "bad"  # another step may not be deferred
    _archive4(bad, full, env=_env4() + DIAG_ENV_LINES)
    (bad / "06-a7-step-independence.txt").write_text((d / "05-a15-workers.txt").read_text())
    assert summ.main([str(bad), "--expected-sha", SHA]) == 1
    # a document of the wrong tag or seed base fails
    wrong = tmp_path / "wrong"
    _archive4(wrong, full, doc={"seed_base": 1})
    assert summ.main([str(wrong), "--expected-sha", SHA]) == 1
    # combine the two halves
    a, b = tmp_path / "a", tmp_path / "b"
    _archive4(a, full[:6])
    _archive4(b, full[6:])
    out = tmp_path / "combined.json"
    assert summ.main(["--combine", str(a), str(b), "--expected-sha", SHA, "--out", str(out)]) == 0
    c = json.loads(out.read_text())
    assert c["pass"] and c["complete"] and c["conformant"] and not c["missing_steps"]


def test_v4_runner_single_process_smoke_end_to_end(tmp_path: Path) -> None:
    """The cheapest steps of lv4 through the real runner (a8 at scale 0.01 and the deferred
    workers step), single-process, rehearsal base."""
    repo = tmp_path / "repo"
    sha = _make_repo(repo)
    out = repo / "validation" / "generated" / "transport" / "run"
    r = _run("--suite", "lv4", "--out", str(out), "--expected-sha", sha, "--workers", "1",
             "--scale", "0.01", "--seed-base", "20351004", "--only", "05", "07",
             root=repo)  # fmt: skip
    assert r.returncode == 0, r.stderr[-1500:] + r.stdout[-1500:]
    s = json.loads((out / "summary.json").read_text())
    assert s["pass"] and s["suite"] == "lv4" and not s["conformant"] and s["subset"]
    assert s["deferred_steps"] == ["05-a15-workers"]
    assert s["seed_base"] == "20351004"
    assert s["steps"]["07-a8-offline-let"]["reduced"] is True
    assert "validation/plans/v3-004-acceptance.md" in (out / "environment.txt").read_text()


def test_a9_parts_roundtrip_and_tamper_detection(tmp_path: Path) -> None:
    """The two A9 parts and the comparison (tiny: 4 seeds of 1000 histories; the verdict itself is
    not asserted, only the plumbing): the comparison verifies the producer archive, the file hashes,
    the metadata and the deterministic seed split, and refuses tampering."""
    sha = "c" * 40
    run_suite = _load("run_suite")
    current = {
        str(f.relative_to(run_suite.REPO)): run_suite.sha256(f)
        for f in run_suite.source_files("lv4")
    }
    env = dict(os.environ, IONMC_RUN_SHA=sha, IONMC_RUN_SUITE="lv4", PYTHONPATH=str(REPO / "src"))
    steps = str(SCRIPTS / "steps_v4.py")
    own, prod = tmp_path / "own", tmp_path / "prod"
    names = ["01-a9-part-1of2", "02-a9-part-2of2"]
    for d in (own, prod):
        d.mkdir()
        text = _env(hashes=current).replace(SHA, sha).replace("suite=lv", "suite=lv4")
        (d / "environment.txt").write_text(text.replace("seed_base=20391004", "seed_base=20351004"))
    (prod / "manifest.txt").write_text("".join(f"{n}\n" for n in names))
    common = ["--scale", "0.01", "--seeds", "4", "--seed-base", "20351004"]

    def step(*args: str, environ: dict | None = None) -> subprocess.CompletedProcess[str]:  # type: ignore[type-arg]
        return subprocess.run([sys.executable, steps, *args, *common], capture_output=True,
                              text=True, env=environ or env, timeout=900)  # fmt: skip

    samples = str(prod / "samples")
    for i in (1, 2):
        p = step("a9-part", "--part", f"{i}/2", "--out-dir", samples)
        assert p.returncode == 0, p.stdout[-1500:] + p.stderr[-1500:]
        body = f"# command: x\n# git_sha: {sha}\n# started_utc: n\n# step_timeout_s: 1\n"
        (prod / f"{names[i - 1]}.txt").write_text(body + p.stdout + "\n# exit=0\n")
    args = ("a9-compare", "--dirs", str(own), str(prod))
    c = step(*args)
    assert "#JSON-BEGIN" in c.stdout, c.stderr[-1500:]
    doc = json.loads(c.stdout.split("#JSON-BEGIN")[1].split("#JSON-END")[0])
    assert doc["step"] == "a9-compare" and doc["seeds"] == 4 and doc["reduced"]
    assert doc["criterion"]["z_std_band"] == [0.87, 1.13] and set(doc["coverage"]) == {
        "plateau_0.5R", "bragg_peak", "distal_80pct"}  # fmt: skip

    def refused(proc: subprocess.CompletedProcess[str], text: str) -> None:
        assert proc.returncode != 0 and text in (proc.stderr + proc.stdout), proc.stderr[-800:]

    refused(step(*args, environ=dict(env, IONMC_RUN_SHA="d" * 40)), "SHA")
    f = prod / "samples" / "a9-part-1-of-2.npz"
    good = f.read_bytes()
    f.write_bytes(good + b"x")
    refused(step(*args), "sha256 mismatch")
    f.write_bytes(good)
    (prod / "samples" / "a9-part-2-of-2.npz").unlink()
    assert step(*args).returncode != 0  # a missing part


def test_a9_amended_criterion_synthetic() -> None:
    """Plan footnote 5: z sd in [0.87, 1.13], |mean z| < 0.18, mean coverage in [0.93, 0.97]."""
    mod = _load("steps_v4")
    rng = np.random.default_rng(7)
    z = rng.standard_normal(200)
    z = (z - z.mean()) / z.std(ddof=1)  # exactly calibrated
    assert mod.a9_depth_pass(z)["calibrated"]
    low = mod.a9_depth_pass(z * 1.18)  # sigma 15 % too low
    assert not low["z_std_ok"] and not low["calibrated"]
    assert mod.a9_depth_pass(z * 1.119)["calibrated"]  # the observed peak value at 20371004
    biased = mod.a9_depth_pass(z + 0.3)
    assert biased["z_std_ok"] and not biased["z_mean_ok"] and not biased["calibrated"]
    assert mod.a9_cover_all_pass(0.9435) and mod.a9_cover_all_pass(0.9464)
    assert not mod.a9_cover_all_pass(0.92) and not mod.a9_cover_all_pass(0.98)


# -- A16 intended-change exception (review finding: bound to V3-003D and verified at run time) ---
def _a16_cmd(mod: ModuleType) -> list[str]:
    steps = mod.suite_steps("lv4", 4, 1.0)
    return next(s[1] for s in steps if s[0].endswith("a16-qualified-path-regression"))


def _v3003d_record() -> dict:  # type: ignore[type-arg]
    """The V3-003D A16 intended-change record as the historical block of the V3-003D plan states
    it (``A16_INTENDED_CHANGE`` of ``run_suite.py`` is ``None`` since V3-005A): the fixture of the
    gate and binding tests, which exercise the machinery a future task re-introduces."""
    st = _load("steps_v4")
    text = (REPO / st.PLAN_FILE).read_text()
    block = st.plan_block(text)
    inner = block[len(st.PLAN_BEGIN) : -len(st.PLAN_END)]
    rec = json.loads(inner[inner.index("{") : inner.rindex("}") + 1])
    rec["plan_block_sha256"] = hashlib.sha256(block.encode()).hexdigest()
    return rec  # type: ignore[no-any-return]


def test_a16_has_no_recorded_exception_and_the_baseline_is_the_v3_003d_merge() -> None:
    run, st = _load("run_suite"), _load("steps_v4")
    assert run.A16_INTENDED_CHANGE is None  # V3-003D record deleted (amendment 6 of V3-004)
    assert re.fullmatch(r"[0-9a-f]{40}", st.A16_BASELINE)  # a full develop merge SHA
    assert st.A16_BASELINE.startswith("f3a1dd62")  # V3-003D merge commit
    assert not st.A16_BASELINE.startswith(_v3003d_record()["baseline"])  # advanced past a524f209
    cmd = _a16_cmd(run)  # the suite runs the gated regression mode, without a record
    assert cmd[cmd.index("--mode") + 1] == "regression"
    assert "--intended-change-record" not in cmd


def test_a16_suite_passes_a_recorded_exception_to_the_intended_change_mode() -> None:
    mod = _load("run_suite")
    rec = _v3003d_record()
    ident = ("task", "baseline", "identity_field", "baseline_value", "new_value")
    assert {k: rec[k] for k in ident} == {
        "task": "V3-003D", "baseline": "a524f209", "identity_field": "range_construction",
        "baseline_value": None, "new_value": "exact-loglog-quadrature-v1",
    }  # fmt: skip
    mod.A16_INTENDED_CHANGE = rec  # a future task records its exception
    cmd = _a16_cmd(mod)
    assert cmd[cmd.index("--mode") + 1] == "intended-change"
    assert json.loads(cmd[cmd.index("--intended-change-record") + 1]) == rec


def test_a16_intended_change_verification_is_fail_closed() -> None:
    st = _load("steps_v4")
    rec = {**_v3003d_record(), "baseline": st.A16_BASELINE[:8]}  # a record naming this baseline
    new = rec["new_value"]
    ok = st.verify_intended_change(rec, st.A16_BASELINE, None, new)
    assert ok["verified_baseline_identity"] is None and ok["verified_current_identity"] == new
    with pytest.raises(SystemExit, match="needs the A16_INTENDED_CHANGE record"):
        st.verify_intended_change(None, st.A16_BASELINE, None, new)  # no record
    with pytest.raises(SystemExit, match="already carries"):
        st.verify_intended_change(rec, st.A16_BASELINE, new, new)  # baseline merged: delete record
    with pytest.raises(SystemExit, match="record expects"):
        st.verify_intended_change(rec, st.A16_BASELINE, None, None)  # no change in the tree
    with pytest.raises(SystemExit, match="record expects"):
        st.verify_intended_change(rec, st.A16_BASELINE, None, "other-construction")
    with pytest.raises(SystemExit, match="names baseline"):
        st.verify_intended_change({**rec, "baseline": "deadbeef"}, st.A16_BASELINE, None, new)
    with pytest.raises(SystemExit, match="names baseline"):  # the V3-003D record is out of date
        st.verify_intended_change(_v3003d_record(), st.A16_BASELINE, None, new)
    with pytest.raises(SystemExit, match="exactly the keys"):
        st.verify_intended_change({"task": "V3-003D"}, st.A16_BASELINE, None, new)


def test_a16_step_refuses_intended_change_without_record() -> None:
    st = _load("steps_v4")
    ns = argparse.Namespace(mode="intended-change", intended_change_record=None)
    with pytest.raises(SystemExit, match="needs the A16_INTENDED_CHANGE record"):
        st.step_a16(ns)


def test_lv_pytest_step_requires_the_nist_cache() -> None:
    mod = _load("run_suite")
    step = next(s for s in mod.suite_steps("lv", 4, 1.0) if "pytest-warp-cpu" in s[0])
    assert "tests/ionmc/test_range_quadrature.py" in step[1]
    assert step[2]["IONMC_REQUIRE_NIST"] == "1" and step[2]["IONMC_CACHE_DIR"].endswith(
        ".ionmc-cache/ionmc-data"
    )


# -- A16 intended-change: binding to the plan block and the gated comparison ---------------------
def _synthetic_raw(spec: str = "t13:warp-cpu:float64") -> dict[str, np.ndarray]:
    """Small physically valid raw A16 values with the field layout of ``a16_digest.py``."""
    rng = np.random.default_rng(3)
    z = np.arange(40)
    prof = np.exp(-0.5 * ((z - 30) / 6.0) ** 2) + 0.2
    grid = (prof[None, None, None, :] * (1.0 + 0.1 * rng.random((4, 3, 3, 40)))).astype(np.float64)
    n, nt = 20, 10
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    pos = np.column_stack(
        [rng.uniform(-20, 20, n), rng.uniform(-20, 20, n), rng.uniform(100, 150, n)]
    )
    e_end = rng.uniform(1.9, 2.0, n)
    raw = {
        "grid.dose.batch_energy_mev": grid,
        "energy_balance.in_grid_mev": np.array([grid.sum()]),
        "energy_balance.step_deposit_mev": np.array(grid.sum() * 0.99),
        "energy_balance.cutoff_mev": np.array(grid.sum() * 0.01),
        "energy_balance.quantization_mev": np.array([1e-7]),
        "diagnostics.end_position_mm": pos,
        "diagnostics.end_direction": u,
        "diagnostics.end_energy_mev": e_end,
        "diagnostics.trace_end_energy_mev": e_end[:4].copy(),
        "counters.n_steps": np.array(1234),
        "valid": np.array(True),
        "diagnostics.trace.step": np.arange(1.0, nt + 1),
        "diagnostics.trace.history": np.repeat(np.arange(2.0), nt // 2),
        "diagnostics.trace.ix": np.zeros(nt),
        "diagnostics.trace.iy": np.zeros(nt),
        "diagnostics.trace.iz": np.arange(nt, dtype=float),
        "diagnostics.trace.reason": np.ones(nt),
        "diagnostics.trace.attempts": np.ones(nt),
        "diagnostics.trace.blocks": np.arange(2.0, 2.0 * nt + 1, 2.0),
        "diagnostics.trace.energy_mev": np.linspace(90.0, 20.0, nt),
        "diagnostics.trace.x_mm": rng.uniform(-3, 3, nt),
        "diagnostics.trace.y_mm": rng.uniform(-3, 3, nt),
        "diagnostics.trace.z_mm": np.linspace(1.0, 90.0, nt),
        "diagnostics.trace.deposit_mev": rng.uniform(0.5, 1.5, nt),
        "diagnostics.trace.step_mm": np.full(nt, 1.0),
        "diagnostics.trace.ux": np.full(nt, 0.01),
        "diagnostics.trace.uy": np.full(nt, 0.02),
        "diagnostics.trace.uz": np.full(nt, 0.99),
    }
    return {f"{spec}|{k}": v for k, v in raw.items()}


T13, T1 = "t13:warp-cpu:float64", "t1:warp-cpu:float64"


def _perturbed(
    base: dict[str, np.ndarray], spec: str = T13, **edits: object
) -> dict[str, np.ndarray]:
    """The intended change in miniature (small end-depth and energy-balance shifts), then ``edits``
    (``field__with__dots`` -> array or function of the field)."""
    cur = {k: v.copy() for k, v in base.items()}
    cur[f"{spec}|diagnostics.end_position_mm"][:, 2] -= 0.003
    for f, fac in (("step_deposit_mev", 1 + 1e-6), ("cutoff_mev", 1 - 1e-6)):
        cur[f"{spec}|energy_balance.{f}"] = cur[f"{spec}|energy_balance.{f}"] * fac
    for k, v in edits.items():
        key = f"{spec}|" + k.replace("__", ".")
        cur[key] = v(cur[key]) if callable(v) else np.asarray(v)  # type: ignore[operator]
    return cur


def _violations(rec: dict, base: dict, cur: dict, spec: str) -> list[str]:  # type: ignore[type-arg]
    res = _load("steps_v4").a16_gate(rec, base, cur)
    return list(res["specs"][spec]["violations"]) if not res["ok"] else []


def _with(fn):  # type: ignore[no-untyped-def]
    def edit(a: np.ndarray) -> np.ndarray:
        a = a.copy()
        fn(a)
        return a

    return edit


def test_a16_gate_accepts_the_intended_change_and_rejects_everything_else() -> None:
    st = _load("steps_v4")
    rec = _v3003d_record()
    base = _synthetic_raw()
    ok = st.a16_gate(rec, base, _perturbed(base))
    assert ok["ok"], ok
    # a non-allowlisted field (a counter) differs
    bad = _violations(rec, base, _perturbed(base, counters__n_steps=1235), T13)
    assert any("without an allowlisted bound" in v for v in bad)

    # an allowed field leaves its bound (a scoring bug: 5 % more deposit in one layer)
    def bump(g: np.ndarray) -> None:
        g[..., 30] *= 1.05

    assert _violations(rec, base, _perturbed(base, grid__dose__batch_energy_mev=_with(bump)), T13)

    # a 1 mm end-depth shift (wrong kernel branch) exceeds the end-position bound
    def shift(p: np.ndarray) -> None:
        p[:, 2] += 1.0

    assert _violations(rec, base, _perturbed(base, diagnostics__end_position_mm=_with(shift)), T13)
    # the discrete trace row count changing in t13 is a violation
    cur = _perturbed(base)
    cur[f"{T13}|diagnostics.trace.step"] = np.arange(11.0)
    assert not st.a16_gate(rec, base, cur)["ok"]
    # a field that exists only in one tree
    cur = _perturbed(base)
    cur[f"{T13}|counters.extra"] = np.array(1)
    assert not st.a16_gate(rec, base, cur)["ok"]
    # a changed discrete trace column (the t13 digests must stay identical there)
    assert _violations(rec, base, _perturbed(base, diagnostics__trace__reason=lambda a: a + 1), T13)


def test_a16_gate_t1_bounds_every_allowlisted_field() -> None:
    rec = _v3003d_record()
    base = _synthetic_raw("t1:warp-cpu:float64")
    assert _violations(rec, base, _perturbed(base, T1), T1) == []
    # every differing field of the t1 allowlist is covered by existing bound keys
    for fam in ("t1", "t13"):
        for f, keys in rec["allowed_differing_fields"][fam].items():
            assert keys and all(k in rec["bounds"][fam] for k in keys), (fam, f)

    def edit(**kw: object) -> list[str]:
        return _violations(rec, base, _perturbed(base, T1, **kw), T1)

    def moved(p: np.ndarray) -> None:
        p[:15, 0] += 5.0  # 75 % of the histories move by 5 mm

    assert any("moved_gt1mm" in v for v in edit(diagnostics__end_position_mm=_with(moved)))
    assert any(
        "outside" in v
        for v in edit(diagnostics__end_position_mm=_with(lambda p: p.__setitem__((0, 0), 50.0)))
    )
    assert any(
        "end_energy_mean_abs" in v for v in edit(diagnostics__end_energy_mev=lambda e: e - 0.05)
    )
    assert any("over_cut" in v for v in edit(diagnostics__end_energy_mev=lambda e: e + 0.5))
    assert any("end_direction_mean_abs" in v for v in edit(
        diagnostics__end_direction=lambda d: d[::-1].copy()))  # fmt: skip
    assert any("norm_err" in v for v in edit(diagnostics__end_direction=lambda d: d * 1.001))
    assert any(
        "trace_energy_over_e0" in v for v in edit(diagnostics__trace__energy_mev=lambda e: e + 20)
    )
    assert any(
        "trace_step_over_max" in v for v in edit(diagnostics__trace__step_mm=lambda s: s + 1.5)
    )
    assert any(
        "deposit_row_max" in v for v in edit(diagnostics__trace__deposit_mev=lambda d: d * 4)
    )
    assert any(
        "trace_position_outside" in v for v in edit(diagnostics__trace__x_mm=lambda x: x + 40)
    )
    assert any("direction_over_unit" in v for v in edit(diagnostics__trace__uz=lambda u: u + 0.5))
    assert any("attempts_max" in v for v in edit(diagnostics__trace__attempts=lambda a: a + 9))
    assert any("discrete_out_of_range" in v for v in edit(diagnostics__trace__iz=lambda a: a + 30))
    assert any(
        "discrete_out_of_range" in v for v in edit(diagnostics__trace__reason=lambda a: a + 4)
    )
    assert any(
        "trace_rows_rel" in v or "field" in v
        for v in edit(diagnostics__trace__energy_mev=lambda e: np.append(e, np.arange(8.0)))
    )
    assert any("in_grid_rel" in v for v in edit(energy_balance__in_grid_mev=lambda a: a * 1.001))

    def spike(g: np.ndarray) -> None:
        g[..., 10] *= 8

    assert any("profile" in v or "depth-dose" in v for v in edit(
        grid__dose__batch_energy_mev=_with(spike)))  # fmt: skip


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
@pytest.mark.parametrize("spec", [T13, T1])
def test_a16_gate_rejects_non_finite_values_in_every_family(spec: str, bad: float) -> None:
    """NaN/Inf must never satisfy a bound or pass a digest comparison (fail closed)."""
    rec = _v3003d_record()
    fam = spec.split(":")[0]
    base = _synthetic_raw(spec)
    assert _violations(rec, base, _perturbed(base, spec), spec) == []
    for field in rec["allowed_differing_fields"][fam]:
        key = f"{spec}|{field}"
        assert key in base, field

        def poison(a: np.ndarray) -> np.ndarray:
            a = np.array(a, dtype=np.float64)
            a.flat[0] = bad
            return a

        for tree in ("current", "baseline", "both"):
            b, c = dict(base), _perturbed(base, spec)
            if tree in ("current", "both"):
                c[key] = poison(c[key])
            if tree in ("baseline", "both"):
                b[key] = poison(b[key])
            v = _violations(rec, b, c, spec)
            assert any("non-finite" in x for x in v), (field, tree, v)
    # a digest-compared (non-allowlisted) field, stable NaN bytes in both trees
    b, c = dict(base), _perturbed(base, spec)
    b[f"{spec}|counters.n_steps"] = c[f"{spec}|counters.n_steps"] = np.array(bad)
    assert any("non-finite" in x for x in _violations(rec, b, c, spec))


def test_a16_caps_are_nan_safe() -> None:
    st = _load("steps_v4")
    assert st._over(np.array([1.0, np.nan]), 5.0) == float("inf")
    assert st._over(np.array([1.0, 7.0]), 5.0) == 2.0 and st._over(np.array([1.0]), 5.0) == 0.0
    assert st._outside(np.array([[0.0, np.nan, 0.0]]), [-1] * 3, [1] * 3) == float("inf")
    mixed = {"a": np.array([1.0]), "b": np.array([np.inf]), "c": np.array(True)}
    mixed["d"] = np.array([np.nan, 1.0])
    assert st.nonfinite_fields(mixed) == ["b", "d"]


_RECORD_TEXT = (
    "head = 1\nA16_INTENDED_CHANGE: dict[str, Any] | None = {\n"
    '    "bounds": {"a": 1},\n    "source_digest": "1111",\n'
    '    "plan_block_sha256": "2222",\n}\n"""doc"""\ntail = 2\n'
)


def test_a16_normalization_masks_exactly_the_record_literals() -> None:
    run = _load("run_suite")
    rec = _v3003d_record()
    text = _RECORD_TEXT
    lits = ['"source_digest": "1111"', '"plan_block_sha256": "2222"']
    # decoys with the same shape before and after the record are not masked
    decoy = '"source_digest": "deadbeef"\n"plan_block_sha256": "cafe"\n'
    doctored = decoy + text + decoy
    out = run.a16_normalize(doctored, "record")
    masked = [f'"{k}": "MASKED"' for k in ("source_digest", "plan_block_sha256")]
    assert sum(out.count(m) for m in masked) == 2 and out.count(decoy) == 2
    # exactly the two literals changed, at their positions; every other byte is identical
    expect = doctored
    for lit, m in zip(lits, masked, strict=True):
        expect = expect.replace(lit, m)
    assert out == expect
    # the plan: one literal (the source digest) inside the block, decoys outside untouched
    plan = (REPO / run.A16_PLAN_FILE).read_text()
    plit = f'"source_digest": "{rec["source_digest"]}"'
    assert plan.count(plit) == 1
    pd = run.a16_normalize(decoy + plan + decoy, "plan")
    assert pd.count('"source_digest": "MASKED"') == 1 and pd.count(decoy) == 2
    assert pd == (decoy + plan + decoy).replace(plit, '"source_digest": "MASKED"')
    # a record or block with a missing or extra literal fails closed
    with pytest.raises(SystemExit, match="expected 2"):
        run.a16_normalize(text.replace(lits[0], '"source_digest_x": "0"'), "record")
    with pytest.raises(SystemExit, match="expected 1"):
        run.a16_normalize(plan.replace(plit, plit + ',\n "plan_block_sha256": "0"'), "plan")
    with pytest.raises(SystemExit, match="not closed"):
        run.a16_normalize(text.replace("}\n", ""), "record")
    # no record: nothing to mask
    assert run.a16_normalize("no record here", "record") == "no record here"


def test_a16_normalization_handles_the_none_form() -> None:
    """With ``A16_INTENDED_CHANGE = None`` (the normal state) nothing is masked, even if a later
    brace-closed literal or look-alike digest literals exist (they stay hashed)."""
    run = _load("run_suite")
    none_form = "A16_INTENDED_CHANGE: dict[str, Any] | None = None\n"
    later = 'X = {\n    "source_digest": "1111",\n    "plan_block_sha256": "2222",\n}\n'
    for text in (
        "head = 1\n" + none_form + '"""doc"""\n',
        "head = 1\n" + none_form + later,
        "head = 1\n" + none_form + '"""doc"""\n}\n' + later,
    ):
        assert run.a16_normalize(text, "record") == text
    # the checked-in file is in the None form: its text is hashed unchanged
    rs = (REPO / run.A16_RUN_SUITE_FILE).read_text()
    assert run.A16_RECORD_START in rs and "\n" + none_form in rs
    assert run.a16_normalize(rs, "record") == rs
    # re-introducing a record in the documented form is masked again
    recorded = rs.replace(none_form, _RECORD_TEXT.split("\n", 1)[1].split('"""doc"""')[0], 1)
    norm = run.a16_normalize(recorded, "record")
    assert recorded != rs and norm.count('"source_digest": "MASKED"') == 1
    assert norm.count('"plan_block_sha256": "MASKED"') == 1


def test_a16_source_digest_binds_the_source_tree() -> None:
    run, st = _load("run_suite"), _load("steps_v4")
    assert (run.A16_PLAN_BEGIN, run.A16_PLAN_END) == (st.PLAN_BEGIN, st.PLAN_END)
    digest = run.a16_source_digest()
    assert re.fullmatch(r"[0-9a-f]{64}", digest)
    entries = run.a16_source_entries()
    assert run.a16_digest_of(entries) == digest
    rs, plan = run.A16_RUN_SUITE_FILE, run.A16_PLAN_FILE
    assert rs in entries and plan in entries

    def changed(rel: str, old: str, new: str) -> str:
        e = dict(entries)
        text = e[rel].decode()
        assert old in text, old
        e[rel] = text.replace(old, new, 1).encode()
        return run.a16_digest_of(e)

    # run_suite.py (the gate bounds' home when a record exists, and the digest function) and the
    # plan block are hashed; only the self-referential literals are masked
    assert changed(rs, "def a16_normalize", "def a16_normalise") != digest  # the digest function
    assert changed(plan, '"voxel_rel_max": 0.037', '"voxel_rel_max": 0.5') != digest  # plan block
    raw = (REPO / plan).read_text()
    lit = f'"source_digest": "{_v3003d_record()["source_digest"]}"'
    assert raw.count(lit) == 1  # the historical block's literal is masked, not hashed
    e = dict(entries)
    e[plan] = run.a16_normalize(raw.replace(lit, '"source_digest": "0"'), "plan").encode()
    assert run.a16_digest_of(e) == digest
    assert changed(rs, '"MASKED"', '"MASKED"') == digest
    assert st.verify_source_digest({"source_digest": digest}, digest) == digest
    with pytest.raises(SystemExit, match="exact source state"):
        st.verify_source_digest({"source_digest": digest}, "0" * 64)
    cmd = [sys.executable, str(SCRIPTS / "run_suite.py"), "--print-a16-source-digest"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    assert out.stdout.strip() == digest


def test_a16_plan_binding_is_verified() -> None:
    st = _load("steps_v4")
    rec = _v3003d_record()
    text = (REPO / st.PLAN_FILE).read_text()
    assert st.verify_plan_binding(rec, text) == rec["plan_block_sha256"]
    with pytest.raises(SystemExit, match="hashes to"):
        st.verify_plan_binding(rec, text.replace("in_grid_rel", "in_grid_rell", 1))
    with pytest.raises(SystemExit, match="hashes to"):
        st.verify_plan_binding(
            {**rec, "bounds": rec["bounds"]} | {"plan_block_sha256": "0" * 64}, text
        )
    changed = {
        **rec,
        "bounds": {**rec["bounds"], "t13": {**rec["bounds"]["t13"], "voxel_rel_max": 0.5}},
    }
    with pytest.raises(SystemExit, match="does not state the record"):
        st.verify_plan_binding(changed, text)
    with pytest.raises(SystemExit, match="exactly one delimited"):
        st.verify_plan_binding(rec, "no block here")


def test_v3_005_acceptance_plan_is_frozen_and_consistent() -> None:
    plan = (REPO / "validation" / "plans" / "v3-005-acceptance.md").read_text(encoding="utf-8")
    rows = re.findall(r"^\| ([A-Za-z0-9-]+) \|", plan, flags=re.MULTILINE)
    for row in ("P1", "P2", "P3", "P4", "P5", "N1", "V1", "V1b", "V2", "V2-probe", "V2b", "V3"):
        assert row in rows, row
    for row in ("V4", "V5", "V6", "V7", "V8", "V9", "R1", "X1", "C1", "D6", "E1"):
        assert row in rows, row
    d6 = next(line for line in plan.splitlines() if line.startswith("| D6 |"))
    assert "Tier 1" in d6 and "Tier 2" in d6 and "G(150) ≤ 0.05" in d6 and "D ≤ 1e-3" in d6
    assert "20421004" in plan and "9448d5e" in plan and "0041" in plan
    r1 = next(line for line in plan.splitlines() if line.startswith("| R1 |"))
    assert "`A16_INTENDED_CHANGE` stays `None`" in r1
    assert "DEFERRED" in plan and "partition invariance of nuclear runs" in plan
    assert (REPO / "decisions" / "0041-proton-nuclear-interactions.md").is_file()


# -- V3-005A slice-A suite (lv5) -------------------------------------------------------------------
QUAL5 = 20421004
PLAN5 = REPO / "validation" / "plans" / "v3-005-acceptance.md"
LV5_STEPS = [
    "lv5-throughput", "n1-v1-v1b-d6",
    "v2-100-s0", "v2-100-s1", "v2-100-s2", "v2-150-s0", "v2-150-s1", "v2-150-s2",
    "v2-200-s0", "v2-200-s1", "v2-200-s2", "v2-combine",
    "v2-probe-s05-s0", "v2-probe-s05-s1", "v2-probe-fe-s0", "v2-probe-fe-s1", "v2-probe-fe-s2",
    "v2-probe-combine", "v3-lv", "v4-v4b", "x1", "e1", "r1-a16-t1-regression",
    "v3-workers-partition",
]  # fmt: skip
# plan row -> lv5 step that evaluates it; the rest is CI, slice B, or covered by a pytest tier
LV5_ROW_STEP = {
    "N1": "n1-v1-v1b-d6", "V1": "n1-v1-v1b-d6", "V1b": "n1-v1-v1b-d6", "D6": "n1-v1-v1b-d6",
    "V2": "v2-combine", "V2-probe": "v2-probe-combine", "V3": "v3-lv", "V4": "v4-v4b",
    "V4b": "v4-v4b", "X1": "x1", "E1": "e1", "R1": "r1-a16-t1-regression",
}  # fmt: skip
LV5_NOT_IN_SUITE = {"P1", "P2", "P3", "P4", "P5", "V9", "C1"}  # CI tier (tests/ionmc)
LV5_SLICE_B = {"V2b", "V5", "V6", "V7", "V8"}


def _plan_rows() -> list[str]:
    rows = re.findall(r"^\| ([A-Z][A-Za-z0-9-]*) \|", PLAN5.read_text(), flags=re.M)
    return [r for r in rows if r not in ("Tier", "Row(s)", "#", "CI", "LV", "HR")]


def test_lv5_manifest_seeds_hashed_files_and_deferred() -> None:
    mod = _load("run_suite")
    names = mod.full_step_names("lv5", 2)
    assert [n.split("-", 1)[1] for n in names] == LV5_STEPS
    assert names_sorted(names)
    assert mod.DEFAULT_SEED_BASES["lv5"] == QUAL5 == mod.V5_QUALIFICATION_SEED_BASE
    for suite in ("lv", "lv4", "hr4"):  # the earlier suites are unchanged
        assert mod.DEFAULT_SEED_BASES[suite] != QUAL5
    steps = mod.suite_steps("lv5", 1, 0.5)
    scripted = [s for s in steps if "steps_v5.py" in " ".join(s[1])]
    assert len(scripted) == len(steps)  # every lv5 step is a steps_v5 step (no CUDA, no hr5)
    for name, cmd, env in steps:
        assert cmd[cmd.index("--seed-base") + 1] == str(QUAL5), name
        assert "--timeout" in cmd and env["IONMC_REQUIRE_DATA"] == "1", name
        assert "cuda" not in " ".join(cmd).lower() and "warp" not in cmd[2:3], name
        if name.split("-", 1)[1] != "lv5-throughput" and "combine" not in name:
            assert "--scale" in cmd or "workers-partition" in name or "n1-" in name, name
    assert mod.deferred_step_names("lv5") == ["24-v3-workers-partition"]
    for e in (100, 150, 200):  # the long steps get the 3300 s floor
        assert mod.step_timeout_s("lv5", f"03-v2-{e}-s0", 1500) == 3300
    assert mod.step_timeout_s("lv5", "03-v2-100-s0", 5000) == 5000
    assert mod.step_timeout_s("lv5", "12-v2-combine", 1500) == 1500
    got = {str(f.relative_to(mod.REPO)) for f in mod.source_files("lv5")}
    assert {
        "validation/plans/v3-005-acceptance.md",
        "decisions/0041-proton-nuclear-interactions.md",
        "validation/scripts/transport/steps_v5.py",
        "validation/scripts/transport/nuclear_checks.py",
        "validation/plans/v3-004-acceptance.md",
        "src/ionmc/nuclear/events.py",
        "src/ionmc/nuclear/tables.py",
    } <= got  # fmt: skip
    assert "validation/plans/v3-005-acceptance.md" not in mod.source_file_list("lv4")
    # lv5 always runs one worker (plan, Execution under the single-process directive)
    assert all("--workers" not in c or c[c.index("--workers") + 1] == "1" for _, c, _ in steps)


def test_lv5_plan_rows_are_all_accounted_for_and_seeds_match_the_plan() -> None:
    mod, v5 = _load("run_suite"), _load("steps_v5")
    plan = PLAN5.read_text()
    rows = set(_plan_rows())
    assert {"N1", "V1", "V1b", "V2", "V2-probe", "V3", "V4", "V4b", "X1", "E1", "D6", "R1"} <= rows
    accounted = set(LV5_ROW_STEP) | LV5_NOT_IN_SUITE | LV5_SLICE_B | {"V3-CI"}
    assert rows - accounted == set(), (
        rows - accounted
    )  # a new plan row needs a step or a declaration
    suite = {n.split("-", 1)[1] for n in mod.full_step_names("lv5", 2)}
    assert set(LV5_ROW_STEP.values()) <= suite
    # the r_index table and the base of the plan equal the declarations of steps_v5.py
    text = re.search(r"`r_index` is fixed in `steps_v5.py` \(([^)]*)\)", plan)
    assert text
    pairs = {k.lower(): int(v) for k, v in re.findall(r"([A-Za-z0-9-]+): (\d+)", text.group(1))}
    assert pairs == {"v2-100": 1, "v2-150": 2, "v2-200": 3, "v2-probe-s05": 4, "v2-probe-fe": 5,
                     "v3-lv": 6, "x1": 7, "e1": 8, "v4": 9}  # fmt: skip
    assert pairs == v5.R_INDEX
    assert (v5.R_INDEX["v2-probe-s05"], v5.R_INDEX["v2-probe-fe"], v5.R_INDEX["v4"]) == (4, 5, 9)
    assert "20421004" in plan and v5.QUALIFICATION_SEED_BASE == QUAL5
    assert v5.REHEARSAL_SEED_BASE == 20431004
    nc = _load("nuclear_checks")
    assert nc.V4_BASE_SEED == QUAL5 + 1000 * v5.R_INDEX["v4"]
    # frozen history counts of the plan
    assert (v5.V2_SHARDS * v5.V2_SHARD_N, v5.V3_HISTORIES, v5.X1_N, v5.E1_N, v5.V4_EVENTS) == (
        240_000, {150.0: 20_000, 250.0: 10_000}, 20_000, 100_000, 100_000)  # fmt: skip
    assert v5.PROBES["s05"][0] * v5.PROBES["s05"][1] == 100_000
    assert v5.PROBES["fe"][0] * v5.PROBES["fe"][1] >= 100_000
    assert (v5.V2_TOL, v5.PROBE_TOL, v5.PROBE_SIGMA_MAX, v5.X1_Z) == (0.003, 0.002, 7e-4, 3.0)


def test_lv5_summary_tags_seed_blockers_and_step_documents() -> None:
    summ, mod = _load("summarize"), _load("run_suite")
    v5src = (SCRIPTS / "steps_v5.py").read_text()
    printed = set(re.findall(r'"step": "([a-z0-9-]+)"', v5src))
    tags = {summ.expected_tag(n) for n in mod.full_step_names("lv5", 2)}
    assert tags == printed, tags ^ printed  # every step prints the document its name expects
    assert summ.seed_blockers(QUAL5, "lv5") == []
    assert summ.seed_blockers(20431004, "lv5") and summ.seed_blockers(20401004, "lv5")
    assert summ.seed_blockers(QUAL5, "lv4")  # the lv5 base does not qualify the V3-004 suite
    assert summ.V5_QUALIFICATION_SEED_BASE == QUAL5


def test_v5_partials_are_hash_verified_and_bound(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Shard partials carry ``content_sha256`` and the run/table bindings; an altered value or a
    different table id is rejected, valid partials load (decision 0041, finding 5)."""
    v5 = _load("steps_v5")
    monkeypatch.setenv("IONMC_RUN_SHA", SHA)
    monkeypatch.setenv("IONMC_RUN_SUITE", "lv5")
    monkeypatch.setattr(
        v5, "table_record", lambda: {"table_id": v5.TABLE_ID, "npz_sha256": "c" * 64,
                                     "json_sha256": "d" * 64}
    )  # fmt: skip
    a = argparse.Namespace(out_dir=str(tmp_path), scale=1.0, dirs=[str(tmp_path)])
    v5.write_partial(a, "p-s0", {"row": "p", "shard": 0, "ratio": [1.0, 0.5], "valid": True})
    docs = v5.load_partials(a, ["p-s0.json"])
    assert docs[0]["run_sha"] == SHA and docs[0]["table_id"] == v5.TABLE_ID
    assert docs[0]["content_sha256"] == v5.content_digest(docs[0])
    path = tmp_path / "p-s0.json"
    good = json.loads(path.read_text())
    bad = {**good, "ratio": [1.0, 0.5000001]}  # one altered value, hash left in place
    path.write_text(json.dumps(bad))
    with pytest.raises(SystemExit, match="content_sha256"):
        v5.load_partials(a, ["p-s0.json"])
    other = {**good, "table_id": "e" * 64}  # re-sealed, but bound to a different table
    other["content_sha256"] = v5.content_digest(other)
    path.write_text(json.dumps(other))
    with pytest.raises(SystemExit, match="table_id"):
        v5.load_partials(a, ["p-s0.json"])
    path.write_text(json.dumps(good))
    assert v5.load_partials(a, ["p-s0.json"])[0]["ratio"] == [1.0, 0.5]


def test_v5_imported_partials_need_the_attested_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Imported partials are accepted only with their exact digest in ``--partials-manifest``
    (the PARTIAL stdout lines of the shard steps); a value altered together with a recomputed
    ``content_sha256`` or a partial absent from the manifest is refused (Codex finding 3)."""
    v5, rs = _load("steps_v5"), _load("run_suite")
    monkeypatch.setenv("IONMC_RUN_SHA", SHA)
    monkeypatch.setenv("IONMC_RUN_SUITE", "lv5")
    monkeypatch.setattr(
        v5, "table_record", lambda: {"table_id": v5.TABLE_ID, "npz_sha256": "c" * 64,
                                     "json_sha256": "d" * 64}
    )  # fmt: skip
    shard_dir, cur = tmp_path / "shard", tmp_path / "cur"
    cur.mkdir()
    sa = argparse.Namespace(out_dir=str(shard_dir), scale=1.0, dirs=[str(shard_dir)])
    v5.write_partial(sa, "p-s0", {"row": "p", "shard": 0, "ratio": [1.0, 0.5], "valid": True})
    line = capsys.readouterr().out.strip()
    path = shard_dir / "p-s0.json"
    good = json.loads(path.read_text())
    assert line == f"PARTIAL p-s0.json {good['content_sha256']}"
    man = tmp_path / "manifest.json"
    a = argparse.Namespace(out_dir=str(cur), scale=1.0, dirs=[str(cur), str(shard_dir)],
                           partials_manifest=str(man))  # fmt: skip
    with pytest.raises(SystemExit, match="partials-manifest"):  # imported dirs, no manifest
        v5.load_partials(argparse.Namespace(**{**vars(a), "partials_manifest": None}),
                         ["p-s0.json"])  # fmt: skip
    man.write_text(json.dumps({}))
    with pytest.raises(SystemExit, match="not in the partials manifest"):
        v5.load_partials(a, ["p-s0.json"])
    man.write_text(json.dumps({"p-s0.json": good["content_sha256"]}))
    assert v5.load_partials(a, ["p-s0.json"])[0]["ratio"] == [1.0, 0.5]
    forged = {**good, "ratio": [1.0, 0.9]}  # altered AND resealed: the file is self-consistent
    forged["content_sha256"] = v5.content_digest(forged)
    path.write_text(json.dumps(forged))
    with pytest.raises(SystemExit, match="differs from the manifest"):
        v5.load_partials(a, ["p-s0.json"])
    # run_suite refuses imported directories of lv5 without a manifest and forwards it otherwise
    with pytest.raises(SystemExit, match="partials-manifest"):
        rs.suite_steps("lv5", 1, 1.0, out=cur, import_dirs=[str(shard_dir)])
    steps = rs.suite_steps("lv5", 1, 1.0, out=cur, import_dirs=[str(shard_dir)],
                           partials_manifest=str(man))  # fmt: skip
    comb = [s[1] for s in steps if "combine" in s[0]]
    assert comb and all(str(man) in c for c in comb)
