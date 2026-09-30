import subprocess
from pathlib import Path

import pytest

from infrastructure.experiment_v3 import worker_boundary
from infrastructure.host_runner import host_runner


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


@pytest.fixture
def registered(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    git(root, "init", "-b", "v3/develop")
    (root / "README.md").write_text("tracked source")
    (root / ".gitignore").write_text("validation/generated/\nbenchmarks/generated/\n")
    git(root, "add", ".")
    git(
        root,
        "-c",
        "user.name=Test",
        "-c",
        "user.email=test@example.invalid",
        "commit",
        "-m",
        "baseline",
    )
    wt = tmp_path / "repo-worktrees/v3-test"
    git(root, "worktree", "add", "-b", "task/v3-test-worker", str(wt))
    return root, wt


def test_worker_binds_to_registered_task(registered, tmp_path):
    root, wt = registered
    assert worker_boundary.validate_task("v3-test", wt, root=root) == ("V3-TEST", wt)
    for task in ["../escape", "/tmp/escape", "V2-001", "V3-A/../../escape"]:
        with pytest.raises(ValueError):
            worker_boundary.validate_task(task, wt, root=root)
    for path in [root, tmp_path / "archived-v2", wt.parent / "other"]:
        with pytest.raises(ValueError):
            worker_boundary.validate_task("V3-TEST", path, root=root)
    git(wt, "switch", "-c", "task/v3-other")
    with pytest.raises(ValueError, match="disagree"):
        worker_boundary.validate_task("V3-TEST", wt, root=root)


def test_host_mounts_source_read_only_and_only_ignored_outputs_rw(
    registered, monkeypatch, tmp_path
):
    _, wt = registered
    from infrastructure.experiment_v3 import sandbox_runtime

    monkeypatch.setattr(sandbox_runtime, "require", lambda: "/usr/bin/bwrap")
    exists = Path.exists
    monkeypatch.setattr(Path, "exists", lambda p: str(p) == "/dev/dxg" or exists(p))
    monkeypatch.setattr(host_runner, "HOST_VENV", wt)
    output_root = tmp_path / "protected-run/outputs"
    cmd = host_runner.build_bwrap_command(wt, ["/bin/true"], output_root)
    mounts = [
        (cmd[i], cmd[i + 1], cmd[i + 2])
        for i in range(len(cmd) - 2)
        if cmd[i] in ("--bind", "--ro-bind")
    ]
    assert ("--ro-bind", str(wt), "/workspace") in mounts
    assert ("--bind", str(wt), "/workspace") not in mounts
    assert (
        "--bind",
        str(output_root / "validation/generated"),
        "/workspace/validation/generated",
    ) in mounts
    assert not git(wt, "status", "--porcelain")
    (wt / "validation/generated/forbidden.txt").write_text("tracked")
    git(wt, "add", "-f", "validation/generated/forbidden.txt")
    with pytest.raises(ValueError, match="committed"):
        host_runner.build_bwrap_command(wt, ["/bin/true"], tmp_path / "other-output")


def test_host_run_rejects_source_mutation_even_on_exit_zero(
    registered, monkeypatch, tmp_path
):
    root, wt = registered
    monkeypatch.setattr(host_runner, "REPO", root)
    sha = git(wt, "rev-parse", "HEAD")
    monkeypatch.setattr(
        host_runner,
        "inspect_worktree",
        lambda _: {"path": wt, "branch": "task/v3-test-worker", "sha": sha},
    )
    monkeypatch.setattr(host_runner, "RUN_ROOT", tmp_path / "runs")
    monkeypatch.setattr(host_runner, "CACHE_ROOT", tmp_path / "cache")
    monkeypatch.setattr(
        host_runner, "build_bwrap_command", lambda *a: ["simulated-validation"]
    )
    original = host_runner.subprocess.run

    def execute(cmd, **kwargs):
        if cmd == ["simulated-validation"]:
            (wt / "README.md").write_text("concurrent mutation")
            return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")
        return original(cmd, **kwargs)

    monkeypatch.setattr(host_runner.subprocess, "run", execute)
    result = host_runner.run_validation("V3-TEST", ["/bin/true"], 10)
    assert result["exit_code"] == 0
    assert result["dirty_after"] and not result["committed_state_unchanged"]
    assert not result["succeeded"]


def test_partial_release_retry_requires_successful_asset_upload(tmp_path, monkeypatch):
    import json

    from infrastructure.experiment_v3 import promote

    (tmp_path / ".ionmc-condition.json").write_text(
        json.dumps({"integration_branch": "v3/develop", "release_branch": "v3/main"})
    )
    monkeypatch.setattr(
        promote,
        "evaluate",
        lambda *a: {
            "release_ready": True,
            "code_sha": "source",
            "report_sha256": "hash",
        },
    )

    def git_result(root, *args):
        if args[0] == "branch":
            return "v3/develop"
        if args[0] == "tag":
            return "ionmc-v3-1.0.0"
        if args[-1].endswith("^{tree}"):
            return "same-tree"
        return "source" if args[-1] == "origin/v3/develop" else "release"

    monkeypatch.setattr(promote, "git", git_result)
    monkeypatch.setattr(
        promote.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0)
    )
    uploaded = []
    fail_upload = True

    def gh_result(root, *args):
        if args[:2] == ("pr", "list"):
            return json.dumps(
                [{"headRefOid": "source", "state": "MERGED", "number": 1}]
            )
        if args[:2] == ("pr", "view"):
            return json.dumps(
                {
                    "headRefOid": "source",
                    "state": "MERGED",
                    "statusCheckRollup": [
                        {"name": n, "conclusion": "SUCCESS"}
                        for n in ("pre-commit", "tests", "docs")
                    ],
                }
            )
        assert args[:2] == ("release", "upload")
        uploaded.append(args)
        if fail_upload:
            raise RuntimeError("upload failed")
        return ""

    monkeypatch.setattr(promote, "gh", gh_result)
    promoted = []
    monkeypatch.setattr(promote, "event", lambda *a, **k: promoted.append(k))
    args = (
        tmp_path,
        "plan",
        "freeze",
        tmp_path / "report.json",
        tmp_path,
        "ionmc-v3-1.0.0",
    )
    with pytest.raises(RuntimeError, match="upload failed"):
        promote.promote(*args)
    assert not promoted
    fail_upload = False
    assert promote.promote(*args)["promoted"]
    assert len(uploaded) == 2 and len(promoted) == 1


def test_host_rejects_archived_symlink_and_unregistered_clone(registered, monkeypatch):
    root, wt = registered
    monkeypatch.setattr(host_runner, "REPO", root)
    monkeypatch.setattr(host_runner, "WORKTREE_ROOT", wt.parent)
    assert host_runner.inspect_worktree("V3-TEST")["path"] == wt
    for name in ["v2-001", "v3-archive"]:
        (wt.parent / name).symlink_to(wt, target_is_directory=True)
        with pytest.raises(ValueError):
            host_runner.inspect_worktree(name)
    other = wt.parent / "v3-unregistered"
    subprocess.run(
        ["git", "clone", str(root), str(other)], check=True, capture_output=True
    )
    with pytest.raises(ValueError, match="registered"):
        host_runner.inspect_worktree("V3-UNREGISTERED")


def test_output_publication_preserves_hardlinked_source_and_rejects_links(
    registered, tmp_path
):
    import os

    _, wt = registered
    generated = wt / "validation/generated"
    generated.mkdir(parents=True)
    (wt / "benchmarks/generated").mkdir(parents=True)
    os.link(wt / "README.md", generated / "prior-source-alias")
    run = tmp_path / "RUN-test"
    fresh = run / "outputs/validation/generated"
    fresh.mkdir(parents=True)
    (fresh / "README.md").write_text("new output")
    paths = host_runner.publish_outputs(wt, run)
    assert paths == [str(generated / "RUN-test")]
    assert (generated / "RUN-test/README.md").read_text() == "new output"
    assert (wt / "README.md").read_text() == "tracked source"
    run2 = tmp_path / "RUN-link"
    unsafe = run2 / "outputs/validation/generated"
    unsafe.mkdir(parents=True)
    (unsafe / "escape").symlink_to(wt / "README.md")
    with pytest.raises(ValueError, match="without links"):
        host_runner.publish_outputs(wt, run2)
    assert (wt / "README.md").read_text() == "tracked source"


def test_output_publication_rejects_destination_symlink(registered, tmp_path):
    _, wt = registered
    (wt / "validation").mkdir()
    (wt / "validation/generated").symlink_to(tmp_path, target_is_directory=True)
    run = tmp_path / "RUN-destination"
    (run / "outputs/validation/generated").mkdir(parents=True)
    with pytest.raises(OSError):
        host_runner.publish_outputs(wt, run)
    assert not (tmp_path / "RUN-destination/RUN-destination").exists()


def test_snapshot_source_ignores_mutation_restored_before_exit(registered, tmp_path):
    from infrastructure.experiment_v3.host_snapshot import create

    root, wt = registered
    sha = git(wt, "rev-parse", "HEAD")
    original = (wt / "README.md").read_text()
    snapshot = create(root, wt, sha, "task/v3-test-worker", tmp_path / "snapshot")
    (wt / "README.md").write_text("different source during execution")
    observed = (snapshot / "README.md").read_text()
    (wt / "README.md").write_text(original)
    assert not git(wt, "status", "--porcelain")
    assert observed == original
    assert git(snapshot, "rev-parse", "HEAD") == sha
    assert (snapshot / ".git").is_dir()
    assert (snapshot / "README.md").stat().st_ino != (wt / "README.md").stat().st_ino


def test_snapshot_rejects_symlink_input_directory(registered, tmp_path):
    from infrastructure.experiment_v3.host_snapshot import create

    root, wt = registered
    (wt / ".ionmc-cache").symlink_to(tmp_path)
    with pytest.raises(OSError):
        create(
            root,
            wt,
            git(wt, "rev-parse", "HEAD"),
            "task/v3-test-worker",
            tmp_path / "snapshot",
        )


def test_sandbox_version_gate_rejects_old_runtime(monkeypatch):
    from infrastructure.experiment_v3 import sandbox_runtime

    monkeypatch.setattr(
        sandbox_runtime, "inspect", lambda: {"ready": False, "version": "0.9.0"}
    )
    with pytest.raises(RuntimeError, match="0.12.0"):
        sandbox_runtime.require()


def test_real_sandbox_snapshot_provenance_and_restore_race(
    registered, tmp_path, monkeypatch
):
    import os

    if os.environ.get("IONMC_TEST_HOST_SANDBOX") != "1":
        pytest.skip("operator host sandbox test; set IONMC_TEST_HOST_SANDBOX=1")
    root, wt = registered
    monkeypatch.setattr(host_runner, "REPO", root)
    monkeypatch.setattr(host_runner, "WORKTREE_ROOT", wt.parent)
    monkeypatch.setattr(host_runner, "RUN_ROOT", tmp_path / "protected/runs")
    expected = git(wt, "rev-parse", "HEAD")
    original_run = subprocess.run
    source = wt / "README.md"
    original = source.read_bytes()

    def restore_race(cmd, **kwargs):
        if cmd and Path(cmd[0]).name == "bwrap":
            source.write_text("concurrently changed live source")
            try:
                return original_run(cmd, **kwargs)
            finally:
                source.write_bytes(original)
        return original_run(cmd, **kwargs)

    monkeypatch.setattr(host_runner.subprocess, "run", restore_race)
    code = "\n".join(
        [
            "import subprocess; from pathlib import Path",
            "def git(*a): return subprocess.check_output(['git',*a],text=True).strip()",
            f"assert git('rev-parse','HEAD') == {expected!r}",
            "assert git('branch','--show-current') == 'task/v3-test-worker'",
            "assert not git('status','--porcelain')",
            "assert Path('README.md').read_text() == 'tracked source'",
            "Path('validation/generated/result').write_text('passed')",
            "print('snapshot provenance and restored live mutation passed')",
        ]
    )
    result = host_runner.run_validation("V3-TEST", ["python", "-c", code], 30)
    assert result["succeeded"], result
    assert result["committed_state_unchanged"]
    assert source.read_bytes() == original
    assert (Path(result["output_paths"][0]) / "result").read_text() == "passed"
