"""Archived record recovery, signature preservation and filesystem boundaries."""

import importlib.util
import json
import os
import subprocess
from pathlib import Path

import pytest

from infrastructure.experiment_v3 import attest
from infrastructure.experiment_v3 import host_records as hr

TASK = "V3-008"
RUN = "RUN-20261111T032108Z-0f2c091b"
RUN2 = "RUN-20261111T072739Z-2ca81843"
SHA = "a" * 40


@pytest.fixture
def setup(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()

    def git(*args):
        subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)

    git("init")
    git(
        "-c",
        "user.name=Fixture",
        "-c",
        "user.email=fixture@example.invalid",
        "commit",
        "--allow-empty",
        "-m",
        "fixture",
    )
    worktree = tmp_path / "repo-worktrees/v3-008"
    git("worktree", "add", "-b", "task/v3-008-test", str(worktree))
    state = tmp_path / "state"
    archive = state / "host-runs" / TASK / RUN
    archive.mkdir(parents=True)
    key = attest.keypair_from_seed(bytes(range(32)))
    record = attest.sign_record(
        {
            "task_id": TASK,
            "run_id": RUN,
            "sha": SHA,
            "stdout": "complete output\n",
            "stderr": "",
            "output_manifest_sha256": "b" * 64,
        },
        key,
        "host-runner",
    )
    data = (json.dumps(record, indent=3) + "\n").encode()
    (archive / "record.json").write_bytes(data)
    return root, worktree, state, archive, key, data


def materialize(setup, ids=None):
    root, _, state, *_ = setup
    return hr.materialize_records(TASK, ids or [RUN], root=root, state=state)


def test_signed_bytes_survive_recovery_dirty_worktree_and_retry(setup):
    _, worktree, state, archive, key, data = setup
    (worktree / "unfinished.txt").write_text("ongoing science")
    listing = hr.list_records("v3-008", state=state)
    assert listing["records"][0]["attestation_present"] is True
    result = materialize(setup)
    path = Path(result["files"][0]["path"])
    assert path == worktree / ".ionmc-cache/host-run-records" / (RUN + ".json")
    assert path.read_bytes() == data == (archive / "record.json").read_bytes()
    attest.verify_record(
        json.loads(path.read_bytes()),
        {"host-runner": {key.key_id: key.public_key.hex()}},
        "host-runner",
    )
    assert result["files"][0]["sha256"] == listing["records"][0]["sha256"]
    assert materialize(setup)["files"][0]["status"] == "existing"
    assert (worktree / "unfinished.txt").read_text() == "ongoing science"
    assert set(path.parent.iterdir()) == {path}


def test_missing_unsigned_and_pagination(setup):
    _, _, state, archive, _, _ = setup
    (archive.parent / RUN2).mkdir()
    (archive.parent / RUN2 / "result.json").write_text('{"not": "a signed record"}')
    (archive / "record.json").write_text(
        json.dumps({"task_id": TASK, "run_id": RUN, "sha": SHA})
    )
    page = hr.list_records(TASK, limit=1, state=state)
    assert page["total_runs"] == 2 and page["next_offset"] == 1
    assert page["records"][0]["attestation_present"] is False
    assert hr.list_records(TASK, offset=1, state=state)["records"] == [
        {"run_id": RUN2, "record_available": False}
    ]
    assert materialize(setup)["files"][0]["attestation_present"] is False
    with pytest.raises(FileNotFoundError):
        materialize(setup, [RUN2])


@pytest.mark.parametrize(
    "offset,limit", [(-1, 1), (0, 0), (0, 101), (True, 1), (0, 1.5)]
)
def test_bad_pagination(setup, offset, limit):
    with pytest.raises(ValueError):
        hr.list_records(TASK, offset, limit, state=setup[2])


@pytest.mark.parametrize(
    "task", ["../V3-008", "V1-008", "V2-008", "/V3-008", "V3-008/x"]
)
def test_task_path_and_history_restrictions(setup, task):
    with pytest.raises(ValueError):
        hr.list_records(task, state=setup[2])
    with pytest.raises(ValueError):
        hr.materialize_records(task, [RUN], root=setup[0], state=setup[2])


@pytest.mark.parametrize(
    "ids", [[], [RUN, RUN], ["../record.json"], [RUN + "/x"], [RUN] * 101, [None]]
)
def test_invalid_selection(setup, ids):
    with pytest.raises(ValueError):
        hr.materialize_records(TASK, ids, root=setup[0], state=setup[2])
    assert not (setup[1] / ".ionmc-cache").exists()


@pytest.mark.parametrize("field,value", [("task_id", "V3-OTHER"), ("run_id", RUN2)])
def test_identity_mismatch(setup, field, value):
    path = setup[3] / "record.json"
    record = json.loads(path.read_bytes())
    record[field] = value
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="identity"):
        materialize(setup)


def test_conflict_is_not_overwritten(setup):
    path = Path(materialize(setup)["files"][0]["path"])
    path.write_text("conflict")
    with pytest.raises(ValueError, match="Conflicting"):
        materialize(setup)
    assert path.read_text() == "conflict"


@pytest.mark.parametrize(
    "level",
    [
        ".ionmc-cache",
        ".ionmc-cache/host-run-records",
        ".ionmc-cache/host-run-records/" + RUN + ".json",
    ],
)
def test_destination_symlinks_rejected(setup, tmp_path, level):
    link = setup[1] / level
    link.parent.mkdir(parents=True, exist_ok=True)
    outside = tmp_path / "outside"
    if level.endswith(".json"):
        outside.write_bytes(setup[5])  # Even an identical symlink must be rejected.
    else:
        outside.mkdir()
    link.symlink_to(outside)
    with pytest.raises(OSError):
        materialize(setup)
    if outside.is_dir():
        assert list(outside.iterdir()) == []
    else:
        assert outside.read_bytes() == setup[5]


@pytest.mark.parametrize("level", ["record", "run", "task", "store"])
def test_archive_symlinks_rejected(setup, tmp_path, level):
    archive = setup[3]
    path = {
        "record": archive / "record.json",
        "run": archive,
        "task": archive.parent,
        "store": archive.parent.parent,
    }[level]
    moved = tmp_path / "moved"
    path.rename(moved)
    path.symlink_to(moved)
    with pytest.raises(OSError):
        materialize(setup)


def test_fifo_does_not_block(setup):
    path = setup[3] / "record.json"
    path.unlink()
    os.mkfifo(path)
    with pytest.raises(ValueError, match="regular"):
        materialize(setup)


def test_bounds_and_missing_batch_dont_publish(setup, monkeypatch):
    monkeypatch.setattr(hr, "MAX_RECORD_BYTES", 5)
    with pytest.raises(ValueError):
        materialize(setup)
    monkeypatch.setattr(hr, "MAX_RECORD_BYTES", 4 * 1024 * 1024)
    monkeypatch.setattr(hr, "MAX_BATCH_BYTES", 5)
    with pytest.raises(ValueError, match="Batch"):
        materialize(setup)
    monkeypatch.setattr(hr, "MAX_BATCH_BYTES", 64 * 1024 * 1024)
    with pytest.raises(FileNotFoundError):
        materialize(setup, [RUN, RUN2])
    assert not (setup[1] / ".ionmc-cache").exists()


def test_unregistered_worktree_rejected(setup):
    root, worktree, *_ = setup
    subprocess.run(
        ["git", "-C", str(root), "worktree", "remove", str(worktree)], check=True
    )
    worktree.mkdir()
    with pytest.raises(ValueError, match="registered"):
        materialize(setup)


def test_destination_parent_swap_stays_confined(setup, tmp_path, monkeypatch):
    original = hr._publish
    outside = tmp_path / "outside"
    outside.mkdir()

    def swap(fd, name, data):
        cache = setup[1] / ".ionmc-cache"
        cache.rename(setup[1] / "old-cache")
        cache.symlink_to(outside)
        return original(fd, name, data)

    monkeypatch.setattr(hr, "_publish", swap)
    materialize(setup)
    assert list(outside.iterdir()) == []
    assert (
        setup[1] / "old-cache/host-run-records" / (RUN + ".json")
    ).read_bytes() == setup[5]


def test_mcp_registration_and_dispatch(setup, monkeypatch):
    pytest.importorskip("mcp.server.mcpserver")
    path = Path(__file__).resolve().parents[2] / "infrastructure/mcp/host_runner_mcp.py"
    spec = importlib.util.spec_from_file_location("host_record_mcp_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(
        module, "list_records", lambda *args: hr.list_records(*args, state=setup[2])
    )
    monkeypatch.setattr(
        module,
        "materialize_records",
        lambda *args: hr.materialize_records(*args, state=setup[2], root=setup[0]),
    )
    assert module.list_host_run_records(TASK)["total_runs"] == 1
    assert (
        module.materialize_host_run_records(TASK, [RUN])["files"][0]["status"]
        == "copied"
    )
    import asyncio

    names = {tool.name for tool in asyncio.run(module.mcp.list_tools())}
    assert {
        "run_host_validation",
        "list_host_run_records",
        "materialize_host_run_records",
    } <= names
