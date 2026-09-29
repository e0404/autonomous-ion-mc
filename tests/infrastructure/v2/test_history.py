import json
import runpy
import subprocess
from pathlib import Path

import pytest

from infrastructure.experiment_v2 import history


def git(root, *args):
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True, stderr=subprocess.PIPE
    ).strip()


def commit(root, message):
    git(root, "add", ".")
    git(
        root,
        "-c",
        "user.name=Test",
        "-c",
        "user.email=test@example.invalid",
        "commit",
        "-m",
        message,
    )
    return git(root, "rev-parse", "HEAD")


@pytest.fixture
def isolated(tmp_path):
    remote = tmp_path / "remote"
    remote.mkdir()
    git(remote, "init", "-b", "v2/main")
    (remote / "baseline").write_text("infrastructure only")
    commit(remote, "baseline")
    git(remote, "checkout", "-b", "develop")
    (remote / "science").write_text("excluded v1 implementation")
    excluded = commit(remote, "v1 scientific result")
    git(remote, "tag", "experiment-v1-final")
    git(remote, "checkout", "-b", "v2/develop", "v2/main")
    (remote / ".ionmc-condition.json").write_text(
        json.dumps(
            {
                "integration_branch": "v2/develop",
                "release_branch": "v2/main",
                "history_access": {
                    "mode": "v2-only",
                    "excluded_commit_ids": [excluded],
                },
            }
        )
    )
    commit(remote, "v2 infrastructure")
    clone = tmp_path / "clone"
    # --no-local exercises transfer of reachable objects, not local hardlinks.
    git(
        tmp_path,
        "clone",
        "--no-local",
        "--single-branch",
        "--no-tags",
        "--branch",
        "v2/develop",
        str(remote),
        str(clone),
    )
    git(
        clone,
        "config",
        "--add",
        "remote.origin.fetch",
        "+refs/heads/v2/main:refs/remotes/origin/v2/main",
    )
    git(clone, "fetch", "origin")
    return clone, remote, excluded


def test_independent_clone_allows_task_and_release_workflow(isolated, tmp_path):
    clone, _, excluded = isolated
    assert history.inspect(clone)["ready"]
    assert (
        subprocess.run(
            ["git", "-C", str(clone), "cat-file", "-e", excluded],
            stderr=subprocess.DEVNULL,
        ).returncode
        != 0
    )
    task = tmp_path / "task"
    git(clone, "worktree", "add", "-b", "task/v2-test", str(task))
    (task / "new-code").write_text("new v2 work")
    commit(task, "v2 task")
    git(task, "push", "origin", "HEAD:refs/heads/task/v2-test")
    assert history.inspect(clone)["ready"]
    # Controlled promotion fetches both branches; common ancestry stays intact.
    git(clone, "fetch", "origin", "v2/develop", "v2/main")
    git(clone, "merge-base", "--is-ancestor", "origin/v2/main", "v2/develop")
    assert history.inspect(clone)["ready"]
    assert not history.inspect(task)["checks"]["independent_git_directory"]


def test_shared_v1_worktree_fails(isolated, tmp_path):
    _, remote, _ = isolated
    worktree = tmp_path / "shared"
    git(remote, "worktree", "add", "--detach", str(worktree), "v2/develop")
    result = history.inspect(worktree)
    assert not result["ready"]
    assert not result["checks"]["independent_git_directory"]


@pytest.mark.parametrize(
    "mutation,failed_check",
    [
        (
            [
                "config",
                "--add",
                "remote.origin.fetch",
                "+refs/heads/*:refs/remotes/origin/*",
            ],
            "restricted_fetch",
        ),
        (["config", "remote.origin.tagOpt", "--tags"], "no_automatic_tags"),
        (["remote", "add", "other", "https://example.invalid/other"], "only_origin"),
        (["tag", "experiment-v1-fake"], "only_v2_refs"),
    ],
)
def test_broadened_access_fails(isolated, mutation, failed_check):
    clone, _, _ = isolated
    git(clone, *mutation)
    result = history.inspect(clone)
    assert not result["ready"]
    assert not result["checks"][failed_check]


def test_excluded_fetch_fails_even_without_a_branch(isolated):
    clone, _, excluded = isolated
    git(clone, "fetch", "origin", excluded)
    result = history.inspect(clone)
    assert result["checks"]["only_v2_refs"]
    assert not result["checks"]["excluded_objects_absent"]
    assert not result["ready"]


def test_object_alternates_fail(isolated):
    clone, remote, _ = isolated
    (clone / ".git/objects/info/alternates").write_text(str(remote / ".git/objects"))
    result = history.inspect(clone)
    assert not result["ready"]
    assert not result["checks"]["no_shared_objects"]


def test_missing_policy_fails(isolated):
    clone, _, _ = isolated
    path = clone / ".ionmc-condition.json"
    condition = json.loads(path.read_text())
    del condition["history_access"]
    path.write_text(json.dumps(condition))
    assert not history.inspect(clone)["ready"]


def test_project_uses_separate_v2_telemetry():
    root = Path(__file__).resolve().parents[3]
    for name in (".claude/settings.json", "experiment/v2/claude-settings.json"):
        settings = json.loads((root / name).read_text())
        assert settings["env"]["IONMC_TELEMETRY_DIR"].endswith("/telemetry/v2")
        assert (
            "~/.local/share/ionmc-experiment/telemetry"
            not in (settings["sandbox"]["filesystem"]["allowRead"])
        )


@pytest.mark.parametrize(
    "source,variable,constant",
    [
        (".claude/hooks/log_event.py", "IONMC_TELEMETRY_DIR", "TELEMETRY_DIR"),
        (
            "infrastructure/telemetry/record_event.py",
            "IONMC_TELEMETRY_DIR",
            "TELEMETRY_DIR",
        ),
        ("infrastructure/codex/codex_worker.py", "IONMC_CODEX_RAW_DIR", "RAW_ROOT"),
        (
            "infrastructure/validation/local_validation.py",
            "IONMC_VALIDATION_DIR",
            "VALIDATION_ROOT",
        ),
    ],
)
def test_runtime_expands_project_log_paths(monkeypatch, source, variable, constant):
    root = Path(__file__).resolve().parents[3]
    monkeypatch.setenv(variable, "~/.local/share/ionmc-experiment/test-v2")
    namespace = runpy.run_path(str(root / source))
    assert namespace[constant] == Path.home() / ".local/share/ionmc-experiment/test-v2"


def test_preflight_requires_history_boundary(isolated, tmp_path, monkeypatch):
    from infrastructure.experiment_v2 import preflight, remote
    from infrastructure.host_runner import host_runner

    clone, _, _ = isolated
    path = clone / ".ionmc-condition.json"
    condition = json.loads(path.read_text())
    condition.update(
        experiment_id="experiment-v2", scientific_implementation_inherited=False
    )
    path.write_text(json.dumps(condition))
    commit(clone, "test condition")
    runtime = tmp_path / "runtime"
    (runtime / "bin").mkdir(parents=True)
    (runtime / "bin/python").touch()
    monkeypatch.setattr(host_runner, "HOST_VENV", runtime)
    original_exists = Path.exists
    monkeypatch.setattr(
        Path, "exists", lambda p: str(p) == "/dev/dxg" or original_exists(p)
    )
    monkeypatch.setattr(preflight.shutil, "which", lambda _: "/usr/bin/test")
    monkeypatch.setattr(preflight, "load_engine", lambda *a: {"version": "test"})
    monkeypatch.setattr(
        preflight,
        "run_reference",
        lambda *a, **k: {"succeeded": True, "run_id": "test"},
    )
    monkeypatch.setattr(preflight, "send_notification", lambda *a: {"status": "sent"})
    monkeypatch.setattr(remote, "inspect", lambda _: {"ready": True})
    assert preflight.check(clone, smoke=True, notify=True, state=tmp_path)[
        "ready_for_unattended_launch"
    ]
    git(
        clone,
        "config",
        "--add",
        "remote.origin.fetch",
        "+refs/heads/*:refs/remotes/origin/*",
    )
    result = preflight.check(clone, smoke=True, notify=True, state=tmp_path)
    assert not result["history_isolation"]["ready"]
    assert not result["ready_for_unattended_launch"]
