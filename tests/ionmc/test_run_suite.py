"""Fail-closed behaviour of the LV/HR runner and its summariser (subset labelling, combining,
source attestation, dirty trees); steps are only run for the git/snapshot end-to-end checks."""

from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"
SHA = "a" * 40


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
        r = _run("--suite", "lv", "--out", str(out), "--expected-sha", _head(), "--workers", "2")
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
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", _head(), "--workers", "1")
    assert r.returncode != 0 and not out.exists()
    assert _run("--suite", "lv", "--out", str(out)).returncode != 0  # SHA is mandatory


def test_only_selects_steps_or_fails_before_creating_anything() -> None:
    out = REPO / "validation" / "generated" / "transport" / "test-only"
    assert not out.exists()
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", _head(), "--workers", "2",
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
        "scale=1.0", "seed_base=20331004", "python_parts=2", "only=", "source_hashes:",
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
            full = {**doc, "step": tag, "suite": "lv", "git_sha": SHA, "seed_base": 20331004,
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
    for rel in ("pyproject.toml", "uv.lock", "validation/plans/v3-003-acceptance.md"):
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
    args = ("--suite", "lv", "--out", str(out), "--expected-sha", sha, "--workers", "2",
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
    r = _run("--suite", "lv", "--out", str(out), "--expected-sha", sha, "--workers", "2",
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
    r = _run("--suite", "lv", "--out", str(out2), "--expected-sha", sha, "--workers", "2",
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
        r = _run("--suite", "lv", "--out", str(out3), "--expected-sha", sha, "--workers", "2",
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
        p = step("t12-python-sample", "--part", f"{i}/2", "--out-dir", samples, "--workers", "2")
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
    steps = run_suite.suite_steps("lv", 4, 1.0, seed_base=20331004)
    scripted = [s for s in steps if "steps.py" in " ".join(s[1])]
    assert scripted and all(s[1][s[1].index("--seed-base") + 1] == "20331004" for s in scripted)
    full = _full()
    env = _env()
    good = tmp_path / "good"
    _archive(good, full, doc={"pass": True}, env=env, identity={"seed_base": 20331004})
    assert summ.main([str(good), "--expected-sha", SHA]) == 0
    bad = tmp_path / "bad"  # a step that ran with another base than the archive records
    _archive(bad, full, doc={"pass": True}, env=env, identity={"seed_base": 20271004})
    assert summ.main([str(bad), "--expected-sha", SHA]) == 1
    # archives with different seed bases are not parts of one run
    a, b = full[:6], full[6:]
    _archive(tmp_path / "a", a, doc={"pass": True}, env=env, identity={"seed_base": 20331004})
    other = _env().replace("seed_base=20331004", "seed_base=20271004")
    _archive(tmp_path / "b", b, doc={"pass": True}, env=other, identity={"seed_base": 20271004})
    out = tmp_path / "combined.json"
    assert summ.main(["--combine", str(tmp_path / "a"), str(tmp_path / "b"),
                      "--expected-sha", SHA, "--out", str(out)]) == 1  # fmt: skip
    assert any("differ" in p for p in json.loads(out.read_text())["problems"])


def test_only_the_qualification_seed_base_can_be_conformant(tmp_path: Path) -> None:
    """The default base is the qualification base 20331004; an archive made with a consumed base
    (20271004) or without a recorded base verifies but is never conformant, with the reason
    recorded."""
    run_suite, summ = _load("run_suite"), _load("summarize")
    assert run_suite.DEFAULT_SEED_BASE == 20331004 == summ.QUALIFICATION_SEED_BASE
    for suite_args in (run_suite.suite_steps("lv", 4, 1.0), run_suite.suite_steps("hr", 4, 1.0)):
        scripted = [s for s in suite_args if "steps.py" in " ".join(s[1])]
        assert all(s[1][s[1].index("--seed-base") + 1] == "20331004" for s in scripted)
    full = _full()
    ok = tmp_path / "ok"
    _archive(ok, full, doc={"pass": True})
    assert summ.main([str(ok), "--expected-sha", SHA]) == 0
    s = json.loads((ok / "summary.json").read_text())
    assert s["conformant"] and s["non_conformant_reasons"] == []
    rehearsal = tmp_path / "rehearsal"
    _archive(rehearsal, full, doc={"pass": True},
             env=_env().replace("seed_base=20331004", "seed_base=20271004"),
             identity={"seed_base": 20271004})  # fmt: skip
    assert summ.main([str(rehearsal), "--expected-sha", SHA]) == 0  # it verifies ...
    r = json.loads((rehearsal / "summary.json").read_text())
    assert r["pass"] and not r["conformant"]  # ... but cannot qualify
    assert any("rehearsal" in x for x in r["non_conformant_reasons"])
    missing = tmp_path / "missing"
    _archive(missing, full, doc={"pass": True}, env=_env().replace("seed_base=20331004\n", ""))
    assert summ.main([str(missing), "--expected-sha", SHA]) == 0
    m = json.loads((missing / "summary.json").read_text())
    assert not m["conformant"] and "not recorded" in m["non_conformant_reasons"][0]


@pytest.mark.parametrize(
    "base", [20261004, 20271004, 20281004, 20291004, 20301004, 20311004, 20321004]
)
def test_non_qualification_bases_are_recorded_as_such(tmp_path: Path, base: int) -> None:
    """20261004 (rehearsal), 20271004 (T9 investigation), 20281004, 20291004, 20301004, 20311004 and
    20321004 (first to fifth qualification attempts) verify but never qualify."""
    summ = _load("summarize")
    d = tmp_path / "x"
    env = _env().replace("seed_base=20331004", f"seed_base={base}")
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
