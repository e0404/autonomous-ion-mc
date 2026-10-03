"""Content-addressed cache for external datasets.

Layout under the cache directory::

    objects/<sha256>              raw dataset bytes, named by their SHA-256
    manifests/<dataset_id>.json   provenance record (see :func:`write_manifest`)
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any

CACHE_ENV = "IONMC_CACHE_DIR"
MANIFEST_FIELDS = (
    "dataset_id",
    "version",
    "url",
    "method",
    "post_body",
    "sha256",
    "bytes",
    "retrieved_at",
    "license",
    "citation",
)


class IntegrityError(RuntimeError):
    """Raised when data do not match their pinned SHA-256 or byte count."""


def resolve_cache_dir(explicit: str | os.PathLike[str] | None = None) -> Path:
    """Return the cache directory (no units).

    Precedence: ``explicit`` argument, then environment variable ``IONMC_CACHE_DIR``,
    then ``~/.cache/ionmc``. The directory is not created here.
    """
    if explicit is not None:
        return Path(explicit).expanduser()
    env = os.environ.get(CACHE_ENV)
    if env:
        return Path(env).expanduser()
    return Path.home() / ".cache" / "ionmc"


def sha256_bytes(data: bytes) -> str:
    """Return the hexadecimal SHA-256 digest of ``data`` (bytes in, hex string out)."""
    return hashlib.sha256(data).hexdigest()


def object_path(sha256: str, cache_dir: Path) -> Path:
    """Return the path of the cached object with the given SHA-256."""
    return cache_dir / "objects" / sha256


def manifest_path(dataset_id: str, cache_dir: Path) -> Path:
    """Return the path of the manifest of ``dataset_id``."""
    return cache_dir / "manifests" / f"{dataset_id}.json"


def write_manifest(manifest: dict[str, Any], cache_dir: Path) -> Path:
    """Write a provenance manifest (a dict with the keys in ``MANIFEST_FIELDS``).

    ``bytes`` is a byte count and ``retrieved_at`` a UTC ISO-8601 timestamp.
    """
    missing = [k for k in MANIFEST_FIELDS if k not in manifest]
    if missing:
        raise ValueError(f"manifest lacks fields: {missing}")
    path = manifest_path(str(manifest["dataset_id"]), cache_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_manifest(dataset_id: str, cache_dir: Path) -> dict[str, Any] | None:
    """Return the manifest of ``dataset_id`` or ``None`` if it is not cached."""
    path = manifest_path(dataset_id, cache_dir)
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise IntegrityError(f"manifest {path} is not a JSON object")
    return data


def store_object(data: bytes, cache_dir: Path) -> Path:
    """Store ``data`` under its SHA-256 name (atomic rename) and return the path."""
    digest = sha256_bytes(data)
    path = object_path(digest, cache_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    tmp.replace(path)
    return path


def verify(dataset_id: str, cache_dir: Path | None = None) -> Path:
    """Re-hash the cached object of ``dataset_id`` and return its path.

    Raises ``FileNotFoundError`` if manifest or object is missing and
    :class:`IntegrityError` if the stored bytes no longer match the manifest hash
    or the registry-pinned hash.
    """
    from ionmc.data.registry import DATASETS

    cdir = resolve_cache_dir(cache_dir)
    manifest = read_manifest(dataset_id, cdir)
    if manifest is None:
        raise FileNotFoundError(f"dataset {dataset_id!r} is not cached in {cdir}")
    digest = str(manifest["sha256"])
    path = object_path(digest, cdir)
    if not path.is_file():
        raise FileNotFoundError(f"cached object {path} is missing")
    actual = sha256_bytes(path.read_bytes())
    if actual != digest:
        raise IntegrityError(f"{dataset_id}: object hash {actual} != manifest hash {digest}")
    pinned = DATASETS.get(dataset_id)
    if pinned is not None and pinned.sha256 != digest:
        raise IntegrityError(f"{dataset_id}: manifest hash {digest} != registry {pinned.sha256}")
    return path
