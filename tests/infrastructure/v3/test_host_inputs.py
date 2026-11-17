"""Input storage regressions: reuse, isolation, provenance and admission failures."""

import hashlib
import json
import os
from types import SimpleNamespace

import pytest

from infrastructure.experiment_v3 import host_inputs as hi
from infrastructure.experiment_v3 import storage_guard as sg


@pytest.fixture
def guard(tmp_path, monkeypatch):
    monkeypatch.setattr(sg.platform, "release", lambda: "Linux")
    monkeypatch.delenv("IONMC_HOST_BACKING_VOLUME", raising=False)
    (tmp_path / "host-storage.json").write_text('{"reserve_bytes": 1}')
    return sg.StorageGuard(tmp_path)


def capture(source, destination, store, guard):
    fd = os.open(source, hi.DIRECTORY_FLAGS)
    try:
        with hi.InputStore(store, guard) as objects:
            return objects.capture_tree(fd, destination)
    finally:
        os.close(fd)


def test_repeated_inputs_share_protected_inodes_not_live_files(tmp_path, guard):
    source = tmp_path / "source"
    source.mkdir()
    (source / "empty").mkdir()
    (source / "a").write_bytes(b"reference data" * 1000)
    (source / "b").write_bytes((source / "a").read_bytes())
    store = tmp_path / "store"
    first = capture(source, tmp_path / "one", store, guard)
    second = capture(source, tmp_path / "two", store, guard)
    assert first == second
    files = [tmp_path / n for n in ("one/a", "one/b", "two/a", "two/b")]
    assert len({p.stat().st_ino for p in files}) == 1
    assert files[0].stat().st_ino != (source / "a").stat().st_ino
    assert files[0].stat().st_mode & 0o222 == 0
    assert (tmp_path / "two/empty").is_dir()
    assert len(list(store.iterdir())) == 2  # object and lock; no staging garbage
    (source / "a").write_bytes(b"updated live input")
    assert files[0].read_bytes() == b"reference data" * 1000
    third = capture(source, tmp_path / "three", store, guard)
    assert third != first
    assert (tmp_path / "three/a").stat().st_ino != files[0].stat().st_ino
    assert (tmp_path / "three/b").stat().st_ino == files[0].stat().st_ino


def test_execute_permissions_are_distinct_objects(tmp_path, guard):
    source = tmp_path / "source"
    source.mkdir()
    for name, mode in (("plain", 0o644), ("executable", 0o755)):
        (source / name).write_text("same")
        (source / name).chmod(mode)
    capture(source, tmp_path / "snapshot", tmp_path / "store", guard)
    plain, executable = [
        tmp_path / "snapshot" / name for name in ("plain", "executable")
    ]
    assert plain.stat().st_ino != executable.stat().st_ino
    assert plain.stat().st_mode & 0o777 == 0o444
    assert executable.stat().st_mode & 0o777 == 0o555


@pytest.mark.parametrize("kind", ["symlink", "directory-link", "fifo"])
def test_unsafe_inputs_rejected_without_following(tmp_path, guard, kind):
    source = tmp_path / "source"
    source.mkdir()
    if kind == "fifo":
        os.mkfifo(source / "bad")
    else:
        (source / "bad").symlink_to(tmp_path if kind == "directory-link" else __file__)
    with pytest.raises(ValueError, match="links or special"):
        capture(source, tmp_path / "snapshot", tmp_path / "store", guard)


@pytest.mark.parametrize("damage", ["bytes", "writable", "symlink", "fifo"])
def test_corrupt_existing_object_is_not_used_or_overwritten(tmp_path, guard, damage):
    source = tmp_path / "source"
    source.mkdir()
    (source / "a").write_bytes(b"abc")
    store = tmp_path / "store"
    capture(source, tmp_path / "first", store, guard)
    target = next(p for p in store.iterdir() if p.name != "store.lock")
    if damage in ("symlink", "fifo"):
        target.unlink()
        if damage == "symlink":
            target.symlink_to(source / "a")
        else:
            os.mkfifo(target)
    else:
        target.chmod(0o644)
        if damage == "bytes":
            target.write_bytes(b"xyz")
            target.chmod(0o444)
    with pytest.raises((ValueError, OSError)):
        capture(source, tmp_path / "second", store, guard)
    assert not (tmp_path / "second/a").exists()


@pytest.mark.parametrize("stage", ["hash", "copy", "directory"])
def test_source_mutation_fails_closed(tmp_path, guard, monkeypatch, stage):
    source = tmp_path / "source"
    source.mkdir()
    (source / "a").write_bytes(b"abc")
    original = hi.digest_stream

    def change(stream, outgoing=None, **kwargs):
        result = original(stream, outgoing, **kwargs)
        if stage == "directory":
            (source / "added").touch()
        elif (stage == "copy") == (outgoing is not None):
            before = (source / "a").stat()
            (source / "a").write_bytes(b"xyz")
            # Restoring mtime and file size must not hide a ctime change.
            os.utime(source / "a", ns=(before.st_atime_ns, before.st_mtime_ns))
        return result

    monkeypatch.setattr(hi, "digest_stream", change)
    with pytest.raises(ValueError, match="changed"):
        capture(source, tmp_path / "snapshot", tmp_path / "store", guard)
    assert not list((tmp_path / "store").glob("staging-*"))


def test_failed_copy_is_retriable_without_publishing_partial_object(
    tmp_path, guard, monkeypatch
):
    source = tmp_path / "source"
    source.mkdir()
    (source / "a").write_bytes(b"abc")
    original = hi.digest_stream

    def fail(stream, outgoing=None, **kwargs):
        if outgoing is not None:
            outgoing.write(b"a")
            raise OSError("simulated write failure")
        return original(stream, **kwargs)

    monkeypatch.setattr(hi, "digest_stream", fail)
    with pytest.raises(OSError, match="simulated"):
        capture(source, tmp_path / "failed", tmp_path / "store", guard)
    assert [p.name for p in (tmp_path / "store").iterdir()] == ["store.lock"]
    monkeypatch.setattr(hi, "digest_stream", original)
    capture(source, tmp_path / "retry", tmp_path / "store", guard)
    assert (tmp_path / "retry/a").read_bytes() == b"abc"


def test_guard_accounts_for_new_bytes_and_backing_volume(tmp_path, monkeypatch, guard):
    backing = tmp_path / "backing"
    backing.mkdir()
    (tmp_path / "host-storage.json").write_text(
        json.dumps(
            {
                "reserve_bytes": 10,
                "backing_volume": str(backing),
            }
        )
    )
    monkeypatch.setattr(
        sg.shutil,
        "disk_usage",
        lambda p: SimpleNamespace(free=12 if p == backing else 100),
    )
    guard = sg.StorageGuard(tmp_path)
    guard.require(2)
    with pytest.raises(RuntimeError, match="insufficient free space"):
        guard.require(3)
    monkeypatch.setattr(
        sg.shutil,
        "disk_usage",
        lambda p: SimpleNamespace(free=9 if p == tmp_path else 100),
    )
    with pytest.raises(RuntimeError, match="insufficient free space"):
        guard.require()


def test_low_space_stops_before_copy_and_reuse_needs_no_payload_space(
    tmp_path, guard, monkeypatch
):
    source = tmp_path / "source"
    source.mkdir()
    (source / "a").write_bytes(b"abc")
    store = tmp_path / "store"
    monkeypatch.setattr(sg.shutil, "disk_usage", lambda p: SimpleNamespace(free=2))
    with pytest.raises(RuntimeError, match="insufficient free space"):
        capture(source, tmp_path / "failed", store, guard)
    assert not (tmp_path / "failed/a").exists()
    monkeypatch.setattr(sg.shutil, "disk_usage", lambda p: SimpleNamespace(free=100))
    capture(source, tmp_path / "first", store, guard)
    monkeypatch.setattr(sg.shutil, "disk_usage", lambda p: SimpleNamespace(free=2))
    capture(source, tmp_path / "second", store, guard)
    assert (tmp_path / "first/a").stat().st_ino == (tmp_path / "second/a").stat().st_ino


def test_wsl_requires_explicit_backing_volume(tmp_path, guard, monkeypatch):
    monkeypatch.setattr(sg.platform, "release", lambda: "microsoft-standard-WSL2")
    with pytest.raises(RuntimeError, match="backing-volume check is unconfigured"):
        sg.StorageGuard(tmp_path)
    monkeypatch.setenv("IONMC_HOST_BACKING_VOLUME", str(tmp_path))
    sg.StorageGuard(tmp_path).require()
    monkeypatch.setenv("IONMC_HOST_BACKING_VOLUME", str(tmp_path / "missing"))
    with pytest.raises(ValueError, match="existing absolute"):
        sg.StorageGuard(tmp_path)


def test_manifest_binds_captured_bytes_to_exact_commit(tmp_path, monkeypatch, guard):
    import subprocess

    from infrastructure.experiment_v3.host_snapshot import create

    root = tmp_path / "repo"
    root.mkdir()

    def git(*args):
        return subprocess.check_output(
            ["git", "-C", str(root), *args], text=True
        ).strip()

    git("init", "-b", "task/test")
    (root / "source.txt").write_text("committed")
    git("add", ".")
    git(
        "-c",
        "user.name=Fixture",
        "-c",
        "user.email=fixture@example.invalid",
        "commit",
        "-m",
        "fixture",
    )
    sha = git("rev-parse", "HEAD")
    cache = root / ".ionmc-cache"
    cache.mkdir()
    (cache / "data").write_bytes(b"data")
    for run in ("run1", "run2"):
        parent = tmp_path / run
        parent.mkdir()
        snapshot = create(
            root,
            root,
            sha,
            "task/test",
            parent / "snapshot",
            input_store=tmp_path / "objects",
            guard=guard,
        )
        assert (snapshot / "source.txt").read_text() == "committed"
        manifest = json.loads((parent / "input-manifest.json").read_text())
        assert manifest["sha"] == sha
        entry = manifest["entries"][1]
        assert entry["path"] == ".ionmc-cache/data"
        assert entry["sha256"] == hashlib.sha256(b"data").hexdigest()
        assert (snapshot / entry["path"]).read_bytes() == b"data"
    assert (tmp_path / "run1/snapshot/.ionmc-cache/data").stat().st_ino == (
        tmp_path / "run2/snapshot/.ionmc-cache/data"
    ).stat().st_ino


def test_concurrent_captures_publish_one_object(tmp_path, guard):
    from concurrent.futures import ThreadPoolExecutor

    source = tmp_path / "source"
    source.mkdir()
    (source / "data").write_bytes(b"data" * 10000)
    store = tmp_path / "objects"
    with ThreadPoolExecutor(max_workers=2) as workers:
        results = list(
            workers.map(
                lambda name: capture(source, tmp_path / name, store, guard),
                ("first", "second"),
            )
        )
    assert results[0] == results[1]
    assert len(list(store.iterdir())) == 2
    assert (tmp_path / "first/data").stat().st_ino == (
        tmp_path / "second/data"
    ).stat().st_ino


def test_admission_failure_does_not_create_run_or_launch_process(
    tmp_path, monkeypatch, guard
):
    from infrastructure.host_runner import host_runner

    monkeypatch.setattr(host_runner, "RUN_ROOT", tmp_path / "host-runs")
    monkeypatch.setattr(host_runner, "CACHE_ROOT", tmp_path / "cache")
    monkeypatch.setattr(
        host_runner,
        "inspect_worktree",
        lambda _: {
            "sha": "a" * 40,
            "path": tmp_path,
            "branch": "task/test",
        },
    )
    monkeypatch.setattr(sg.shutil, "disk_usage", lambda p: SimpleNamespace(free=0))
    monkeypatch.setattr(
        host_runner, "create_run_dir", lambda _: pytest.fail("created run")
    )
    monkeypatch.setattr(
        host_runner.subprocess, "run", lambda *a, **k: pytest.fail("launched process")
    )
    with pytest.raises(RuntimeError, match="insufficient free space"):
        host_runner.run_validation("V3-TEST", ["true"], 30)


def test_abrupt_exit_staging_is_recovered_under_lock(tmp_path, guard):
    import subprocess
    import sys

    source = tmp_path / "source"
    source.mkdir()
    (source / "data").write_bytes(b"data" * 1000)
    store = tmp_path / "objects"
    code = """
import os, sys
from pathlib import Path
from types import SimpleNamespace
from infrastructure.experiment_v3 import host_inputs as hi
source, store, destination = map(Path, sys.argv[1:])
original = hi.digest_stream
def crash(stream, outgoing=None, **kwargs):
    if outgoing is not None:
        outgoing.write(b"partial staging bytes")
        outgoing.flush()
        os._exit(77)
    return original(stream, **kwargs)
hi.digest_stream = crash
fd = os.open(source, hi.DIRECTORY_FLAGS)
with hi.InputStore(store, SimpleNamespace(require=lambda *a: None)) as objects:
    objects.capture_tree(fd, destination)
"""
    proc = subprocess.run(
        [
            sys.executable,
            "-B",
            "-c",
            code,
            str(source),
            str(store),
            str(tmp_path / "failed"),
        ],
        check=False,
    )
    assert proc.returncode == 77
    assert len(list(store.glob("staging-*"))) == 1
    capture(source, tmp_path / "retry", store, guard)
    assert not list(store.glob("staging-*"))
    assert (tmp_path / "retry/data").read_bytes() == b"data" * 1000
    assert len(list(store.iterdir())) == 2
