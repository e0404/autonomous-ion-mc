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
    registered, monkeypatch
):
    _, wt = registered
    monkeypatch.setattr(host_runner.shutil, "which", lambda _: "/usr/bin/bwrap")
    exists = Path.exists
    monkeypatch.setattr(Path, "exists", lambda p: str(p) == "/dev/dxg" or exists(p))
    monkeypatch.setattr(host_runner, "HOST_VENV", wt)
    cmd = host_runner.build_bwrap_command(wt, ["/bin/true"])
    mounts = [
        (cmd[i], cmd[i + 1], cmd[i + 2])
        for i in range(len(cmd) - 2)
        if cmd[i] in ("--bind", "--ro-bind")
    ]
    assert ("--ro-bind", str(wt), "/workspace") in mounts
    assert ("--bind", str(wt), "/workspace") not in mounts
    assert (
        "--bind",
        str(wt / "validation/generated"),
        "/workspace/validation/generated",
    ) in mounts
    assert not git(wt, "status", "--porcelain")
    (wt / "validation/generated/forbidden.txt").write_text("tracked")
    git(wt, "add", "-f", "validation/generated/forbidden.txt")
    with pytest.raises(ValueError, match="committed"):
        host_runner.build_bwrap_command(wt, ["/bin/true"])


def test_host_run_rejects_source_mutation_even_on_exit_zero(
    registered, monkeypatch, tmp_path
):
    _, wt = registered
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
