"""Content-addressed local cache for externally acquired datasets.

Layout under the cache root::

    objects/<sha256>                          exact downloaded bytes
    datasets/<dataset_id with / -> __>/<version>.json   provenance records

No dataset is stored in Git; data are downloaded on first use, verified by
SHA-256 and reused offline afterwards.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import urllib.parse
import urllib.request
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

TIMEOUT_S = 60.0


class DataError(Exception):
    """Base class for external-data errors."""


class DatasetUnavailableError(DataError):
    """The dataset is not cached and cannot be downloaded (offline mode)."""


class IntegrityError(DataError):
    """Downloaded or cached bytes do not match the expected SHA-256."""


def default_cache_root() -> Path:
    """Return the cache root: ``IONMC_DATA_DIR``, XDG cache, or ``~/.cache``."""
    env = os.environ.get("IONMC_DATA_DIR")
    if env:
        return Path(env).expanduser()
    xdg = os.environ.get("XDG_CACHE_HOME")
    if xdg:
        return Path(xdg) / "ionmc" / "data"
    return Path.home() / ".cache" / "ionmc" / "data"


@dataclass(frozen=True)
class DatasetSpec:
    """Immutable description of one upstream dataset and how to obtain it.

    ``version`` identifies the upstream evaluation/query and is part of the
    cache key; ``expected_sha256`` optionally pins the exact bytes.
    """

    dataset_id: str
    version: str
    url: str
    method: Literal["GET", "POST"] = "GET"
    post_fields: dict[str, str] | None = None
    license_basis: str = ""
    citation: str = ""
    expected_sha256: str | None = None
    description: str = ""


@dataclass
class DatasetRecord:
    """JSON-serializable provenance record of a cached dataset."""

    dataset_id: str
    version: str
    url: str
    method: str
    post_fields: dict[str, str] | None
    sha256: str
    bytes: int
    retrieved_at_utc: str
    license_basis: str
    citation: str
    expected_sha256: str | None = None
    transformations: list[str] = field(default_factory=list)
    source: str = "network"

    def to_json(self) -> str:
        """Serialize to indented JSON."""
        return json.dumps(asdict(self), indent=2, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> DatasetRecord:
        """Deserialize from :meth:`to_json` output."""
        return cls(**json.loads(text))


def _urllib_download(spec: DatasetSpec) -> bytes:
    from ionmc import __version__

    headers = {"User-Agent": f"ionmc-data/{__version__}"}
    data = None
    if spec.method == "POST":
        data = urllib.parse.urlencode(spec.post_fields or {}).encode("ascii")
    req = urllib.request.Request(
        spec.url, data=data, headers=headers, method=spec.method
    )
    with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp:
        return bytes(resp.read())


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".tmp-")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


class DataCache:
    """Local content-addressed cache with provenance records."""

    def __init__(
        self,
        root: Path | None = None,
        download: Callable[[DatasetSpec], bytes] | None = None,
    ) -> None:
        self.root = Path(root) if root is not None else default_cache_root()
        self._download = download if download is not None else _urllib_download

    def _object_path(self, sha256: str) -> Path:
        return self.root / "objects" / sha256

    def _record_path(self, spec: DatasetSpec) -> Path:
        name = spec.dataset_id.replace("/", "__")
        return self.root / "datasets" / name / f"{spec.version}.json"

    def record(self, spec: DatasetSpec) -> DatasetRecord | None:
        """Return the stored provenance record, or ``None`` if absent."""
        path = self._record_path(spec)
        if not path.is_file():
            return None
        try:
            return DatasetRecord.from_json(path.read_text(encoding="utf-8"))
        except (ValueError, TypeError):
            return None

    def _cached_bytes(self, spec: DatasetSpec) -> tuple[bytes, DatasetRecord] | None:
        rec = self.record(spec)
        if rec is None:
            return None
        obj = self._object_path(rec.sha256)
        if not obj.is_file():
            return None
        data = obj.read_bytes()
        if _sha256(data) != rec.sha256:
            return None
        if spec.expected_sha256 and rec.sha256 != spec.expected_sha256:
            return None
        return data, rec

    def verify(self, spec: DatasetSpec) -> bool:
        """True if record and object exist and the object re-hashes correctly."""
        return self._cached_bytes(spec) is not None

    def fetch(
        self, spec: DatasetSpec, *, offline: bool = False
    ) -> tuple[bytes, DatasetRecord]:
        """Return dataset bytes and record, downloading only if not cached."""
        cached = self._cached_bytes(spec)
        if cached is not None:
            data, rec = cached
            return data, DatasetRecord(**{**asdict(rec), "source": "cache"})
        if offline or os.environ.get("IONMC_OFFLINE") == "1":
            raise DatasetUnavailableError(
                f"Dataset {spec.dataset_id!r} version {spec.version!r} is not in "
                f"the local cache ({self.root}) and offline mode is active. "
                f"Unset IONMC_OFFLINE / offline=False to download it from "
                f"{spec.url}, or copy a populated cache to {self.root}."
            )
        data = self._download(spec)
        digest = _sha256(data)
        if spec.expected_sha256 and digest != spec.expected_sha256:
            raise IntegrityError(
                f"Dataset {spec.dataset_id!r} version {spec.version!r} from "
                f"{spec.url}: SHA-256 {digest} does not match pinned "
                f"{spec.expected_sha256}. Nothing was stored."
            )
        rec = DatasetRecord(
            dataset_id=spec.dataset_id,
            version=spec.version,
            url=spec.url,
            method=spec.method,
            post_fields=dict(spec.post_fields) if spec.post_fields else None,
            sha256=digest,
            bytes=len(data),
            retrieved_at_utc=datetime.now(UTC).isoformat(timespec="seconds"),
            license_basis=spec.license_basis,
            citation=spec.citation,
            expected_sha256=spec.expected_sha256,
            source="network",
        )
        _atomic_write(self._object_path(digest), data)
        _atomic_write(self._record_path(spec), rec.to_json().encode("utf-8"))
        return data, rec

    def list_records(self) -> list[DatasetRecord]:
        """Return all provenance records in the cache, sorted by id/version."""
        out: list[DatasetRecord] = []
        for path in sorted((self.root / "datasets").glob("*/*.json")):
            try:
                out.append(DatasetRecord.from_json(path.read_text(encoding="utf-8")))
            except (ValueError, TypeError):
                continue
        return out
