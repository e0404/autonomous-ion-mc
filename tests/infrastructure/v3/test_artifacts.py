import json

import pytest

from infrastructure.experiment_v3 import artifacts
from infrastructure.experiment_v3.common import file_hash, write_json


@pytest.fixture
def bundle(tmp_path):
    archive = tmp_path / "state/references/REF-test"
    output = archive / "work/dose.bin"
    output.parent.mkdir(parents=True)
    output.write_bytes(b"\x00\x01" * 100000)
    log = archive / "stdout.txt"
    log.write_text("finished\n")
    files = {
        str(p.relative_to(archive)): {"bytes": p.stat().st_size, "sha256": file_hash(p)}
        for p in (output, log)
    }
    write_json(
        archive / "result.json",
        {"code_sha": "abc", "engine": "topas", "artifact_files": files},
    )
    wt = tmp_path / "task"
    wt.mkdir()
    return tmp_path / "state", archive, wt


def test_transfer_large_binary_without_inline_content(bundle):
    state, archive, wt = bundle
    listing = artifacts.list_artifacts("REF-test", limit=1, state=state)
    assert listing["total_files"] == 2 and listing["next_offset"] == 1
    result = artifacts.materialize("REF-test", wt, state=state)
    assert len(json.dumps(result)) < 1000
    assert result["bytes"] == 200009
    assert file_hash(
        wt / ".ionmc-cache/reference-runs/REF-test/work/dose.bin"
    ) == file_hash(archive / "work/dose.bin")
    assert artifacts.materialize("REF-test", wt, state=state) == result
    assert (
        artifacts.read_text("REF-test", "stdout.txt", state=state)["text"]
        == "finished\n"
    )
    with pytest.raises(ValueError, match="Binary"):
        artifacts.read_text("REF-test", "work/dose.bin", state=state)
    with pytest.raises(ValueError):
        artifacts.read_text("REF-test", "stdout.txt", limit=16000, state=state)


def test_corruption_unknown_paths_and_symlinks_rejected(bundle, tmp_path):
    state, archive, wt = bundle
    with pytest.raises(ValueError):
        artifacts.materialize("REF-test", wt, ["../secret"], state=state)
    (archive / "work/dose.bin").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="hash"):
        artifacts.materialize("REF-test", wt, state=state)
    outside = tmp_path / "outside"
    outside.mkdir()
    (wt / ".ionmc-cache").symlink_to(outside)
    with pytest.raises(ValueError):
        artifacts.materialize("REF-test", wt, ["stdout.txt"], state=state)
    assert not list(outside.iterdir())


def test_engine_internal_alias_is_copied_as_bytes(bundle, tmp_path):
    state, archive, wt = bundle
    alias = archive / "work/alias.bin"
    alias.symlink_to("dose.bin")
    record = json.loads((archive / "result.json").read_text())
    record["artifact_files"]["work/alias.bin"] = record["artifact_files"][
        "work/dose.bin"
    ]
    write_json(archive / "result.json", record)
    artifacts.materialize("REF-test", wt, ["work/alias.bin"], state=state)
    target = wt / ".ionmc-cache/reference-runs/REF-test/work/alias.bin"
    assert target.is_file() and not target.is_symlink()
    alias.unlink()
    outside = tmp_path / "private"
    outside.write_bytes(b"secret")
    alias.symlink_to(outside)
    with pytest.raises(ValueError, match="escapes"):
        artifacts.read_text("REF-test", "work/alias.bin", state=state)
