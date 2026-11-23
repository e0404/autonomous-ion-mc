# ruff: noqa: E501
"""Tests of validation/scripts/stopping/env_digest.py (V3-013 sealed environment digests)."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "validation/scripts/stopping/env_digest.py"


def _load():  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location("_t_env_digest", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ED = _load()


def test_tree_digest_is_over_sorted_path_hash_pairs() -> None:
    es = ED._env_seal()
    a = es.tree_digest({"b": "2" * 64, "a": "1" * 64})
    assert a == es.tree_digest({"a": "1" * 64, "b": "2" * 64})
    assert a != es.tree_digest({"a": "1" * 64, "b": "3" * 64})
    assert a != es.tree_digest({"a": "1" * 64, "c": "2" * 64})  # the path is part of the digest


def test_directory_files_excludes_caches_and_hashes_every_file(tmp_path: Path) -> None:
    es = ED._env_seal()
    (tmp_path / "pkg/__pycache__").mkdir(parents=True)
    (tmp_path / "pkg/__pycache__/x.pyc").write_bytes(b"c")
    (tmp_path / "pkg/m.py").write_text("x = 1\n")
    (tmp_path / "pkg/data.bin").write_bytes(b"\x00\x01")
    (tmp_path / "site-packages").mkdir()
    (tmp_path / "site-packages/y.py").write_text("y\n")
    files = es.directory_files(tmp_path)
    assert sorted(files) == ["pkg/data.bin", "pkg/m.py"]


def test_distribution_files_cover_record_and_libs() -> None:
    es = ED._env_seal()
    files = es.distribution_files("numpy")
    assert any(k.endswith("RECORD") for k in files)
    assert any(k.startswith("numpy.libs/") for k in files) or not (
        list((ROOT / ".venv").glob("lib/python*/site-packages/numpy.libs"))
    )
    assert len(files) > 100
    with pytest.raises(es.EnvironmentSealError, match="not installed"):
        es.distribution_files("no-such-distribution-xyz")


def test_differences_name_the_component() -> None:
    a = {"schema": "s", "digests": {"numpy": "1", "warp": "2"}, "versions": {"numpy": "x"},
         "file_counts": {"numpy": 3}, "sys_version": "v", "platform": "p"}  # fmt: skip
    assert ED.differences(a, dict(a)) == []
    b = json.loads(json.dumps(a))
    b["digests"]["warp"] = "3"
    assert ED.differences(a, b) == ["digests.warp"]
    b["platform"] = "q"
    assert ED.differences(a, b) == ["digests.warp", "platform"]
    assert "digests" in ED.differences(a, {})


def test_warp_kernel_cache_attestation(tmp_path: Path) -> None:
    (tmp_path / "wp_mod_1").mkdir()
    (tmp_path / "wp_mod_1/k.cpp").write_text("src\n")
    (tmp_path / "wp_mod_1/k.so").write_bytes(b"\x7fELF")
    att = ED.warp_kernel_cache_attestation(tmp_path)
    assert att["count"] == 2 and sorted(att["files"]) == ["wp_mod_1/k.cpp", "wp_mod_1/k.so"]
    (tmp_path / "wp_mod_1/k.so").write_bytes(b"\x7fELG")
    assert ED.warp_kernel_cache_attestation(tmp_path)["digest"] != att["digest"]
    with pytest.raises(FileNotFoundError):
        ED.warp_kernel_cache_attestation(tmp_path / "absent")


def test_environment_digests_with_a_temporary_uv_lock(tmp_path: Path) -> None:
    (tmp_path / "uv.lock").write_text("lock\n")
    rec = ED.environment_digests(tmp_path)
    assert set(rec["digests"]) == set(ED.COMPONENTS) and rec["schema"] == ED.SCHEMA
    assert all(rec["file_counts"][k] >= 1 for k in ED.COMPONENTS)
    (tmp_path / "uv.lock").write_text("lock2\n")
    assert ED.differences(rec, ED.environment_digests(tmp_path)) == ["digests.uv_lock"]
    assert str(Path.home()) not in json.dumps(rec)  # no absolute path leaks


def test_write_check_round_trip_and_tamper(tmp_path: Path) -> None:
    out = tmp_path / "env.json"
    cmd = [sys.executable, "-I", "-B", str(SCRIPT)]
    w = subprocess.run([*cmd, "--write", str(out)], capture_output=True, text=True, check=False)
    assert w.returncode == 0, w.stderr
    ok = subprocess.run([*cmd, "--check", str(out)], capture_output=True, text=True, check=False)
    assert ok.returncode == 0 and "matches" in ok.stdout
    doc = json.loads(out.read_text())
    doc["digests"]["warp"] = "0" * 64
    out.write_text(json.dumps(doc))
    bad = subprocess.run([*cmd, "--check", str(out)], capture_output=True, text=True, check=False)
    assert bad.returncode == 1 and "digests.warp" in bad.stderr
    out.write_text("not json")
    assert (
        subprocess.run([*cmd, "--check", str(out)], capture_output=True, check=False).returncode
        == 2
    )
