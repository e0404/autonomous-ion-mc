"""Tests for the external-data cache."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from ionmc.data import (
    DataCache,
    DatasetRecord,
    DatasetSpec,
    DatasetUnavailableError,
    IntegrityError,
    default_cache_root,
)

PAYLOAD = b"hello dataset\n"
SPEC = DatasetSpec(
    dataset_id="test/ds",
    version="v1",
    url="https://example.invalid/x",
    method="POST",
    post_fields={"a": "b"},
    license_basis="test",
    citation="cite",
)


class FakeNet:
    def __init__(self, payload: bytes = PAYLOAD) -> None:
        self.payload = payload
        self.calls = 0

    def __call__(self, spec: DatasetSpec) -> bytes:
        self.calls += 1
        return self.payload


def test_fetch_stores_object_and_record(tmp_path: Path) -> None:
    net = FakeNet()
    cache = DataCache(tmp_path, download=net)
    data, rec = cache.fetch(SPEC)
    digest = hashlib.sha256(PAYLOAD).hexdigest()
    assert data == PAYLOAD
    assert rec.source == "network"
    assert rec.sha256 == digest
    assert rec.bytes == len(PAYLOAD)
    assert (tmp_path / "objects" / digest).read_bytes() == PAYLOAD
    assert (tmp_path / "datasets" / "test__ds" / "v1.json").is_file()
    assert cache.verify(SPEC)
    assert [r.dataset_id for r in cache.list_records()] == ["test/ds"]


def test_second_fetch_uses_cache(tmp_path: Path) -> None:
    net = FakeNet()
    cache = DataCache(tmp_path, download=net)
    cache.fetch(SPEC)
    data, rec = cache.fetch(SPEC)
    assert net.calls == 1
    assert rec.source == "cache"
    assert data == PAYLOAD


def test_offline_behaviour(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    net = FakeNet()
    cache = DataCache(tmp_path, download=net)
    with pytest.raises(DatasetUnavailableError, match="test/ds"):
        cache.fetch(SPEC, offline=True)
    monkeypatch.setenv("IONMC_OFFLINE", "1")
    with pytest.raises(DatasetUnavailableError):
        cache.fetch(SPEC)
    monkeypatch.delenv("IONMC_OFFLINE")
    cache.fetch(SPEC)
    monkeypatch.setenv("IONMC_OFFLINE", "1")
    _, rec = cache.fetch(SPEC)
    assert rec.source == "cache"
    assert net.calls == 1


def test_pinned_hash_mismatch_stores_nothing(tmp_path: Path) -> None:
    spec = DatasetSpec(
        dataset_id="test/ds", version="v1", url="u", expected_sha256="0" * 64
    )
    cache = DataCache(tmp_path, download=FakeNet())
    with pytest.raises(IntegrityError):
        cache.fetch(spec)
    assert not (tmp_path / "objects").exists()
    assert cache.record(spec) is None
    assert cache.list_records() == []


def test_pinned_hash_match(tmp_path: Path) -> None:
    spec = DatasetSpec(
        dataset_id="test/ds",
        version="v1",
        url="u",
        expected_sha256=hashlib.sha256(PAYLOAD).hexdigest(),
    )
    _, rec = DataCache(tmp_path, download=FakeNet()).fetch(spec)
    assert rec.expected_sha256 == rec.sha256


def test_corrupted_object_is_refetched(tmp_path: Path) -> None:
    net = FakeNet()
    cache = DataCache(tmp_path, download=net)
    _, rec = cache.fetch(SPEC)
    (tmp_path / "objects" / rec.sha256).write_bytes(b"corrupt")
    assert not cache.verify(SPEC)
    data, rec2 = cache.fetch(SPEC)
    assert data == PAYLOAD
    assert rec2.source == "network"
    assert net.calls == 2
    assert cache.verify(SPEC)


def test_data_dir_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("IONMC_DATA_DIR", str(tmp_path / "d"))
    assert default_cache_root() == tmp_path / "d"
    assert DataCache().root == tmp_path / "d"
    monkeypatch.delenv("IONMC_DATA_DIR")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "x"))
    assert default_cache_root() == tmp_path / "x" / "ionmc" / "data"
    monkeypatch.delenv("XDG_CACHE_HOME")
    assert default_cache_root() == Path.home() / ".cache" / "ionmc" / "data"


def test_record_json_round_trip(tmp_path: Path) -> None:
    _, rec = DataCache(tmp_path, download=FakeNet()).fetch(SPEC)
    rec.transformations.append("parsed")
    assert DatasetRecord.from_json(rec.to_json()) == rec
