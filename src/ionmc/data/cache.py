"""Dataset specifications, cache location, acquisition and integrity checks."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import urllib.request
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

#: Environment variable overriding the cache directory.
CACHE_ENV_VAR = "IONMC_CACHE_DIR"

#: Schema version of the manifest files written next to cached datasets.
MANIFEST_SCHEMA_VERSION = 1


class DatasetNotCached(FileNotFoundError):
    """The requested dataset is not in the local cache (offline use phase)."""


class IntegrityError(ValueError):
    """A downloaded or cached file does not match the specification checksum."""


@dataclass(frozen=True)
class DatasetSpec:
    """Immutable description of one external data file.

    Parameters
    ----------
    name:
        Cache-safe identifier (letters, digits, ``-``, ``_``, ``.``).
    version:
        Immutable version identifier of the source, e.g. a Git commit SHA or
        a release tag. Two specs with the same name but different versions
        are cached side by side.
    url:
        Download location. The URL must resolve to exactly the bytes whose
        checksum is ``sha256``; for Git-hosted files use a commit-pinned raw
        URL, never a branch name.
    sha256:
        Hex SHA-256 of the file contents.
    filename:
        Name of the cached file.
    license:
        License under which the file may be used and redistributed.
    provenance:
        Free-form record of what the data are and where they originate.
    """

    name: str
    version: str
    url: str
    sha256: str
    filename: str
    license: str
    provenance: dict[str, str]

    def __post_init__(self) -> None:
        for field_name in ("name", "version"):
            value = getattr(self, field_name)
            if not value or not all(c.isalnum() or c in "-_." for c in value):
                raise ValueError(f"{field_name}={value!r} is not cache-safe")
        if len(self.sha256) != 64 or any(
            c not in "0123456789abcdef" for c in self.sha256.lower()
        ):
            raise ValueError("sha256 must be a 64-character hex digest")
        if "/" in self.filename or self.filename in ("", ".", ".."):
            raise ValueError(f"unsafe filename {self.filename!r}")

    @property
    def key(self) -> str:
        """Relative cache path of this dataset's directory."""
        return f"{self.name}/{self.version}"


def cache_dir(override: str | os.PathLike[str] | None = None) -> Path:
    """Resolve the dataset cache directory.

    Precedence: explicit ``override`` argument, then the ``IONMC_CACHE_DIR``
    environment variable, then ``$XDG_CACHE_HOME/ionmc``, then
    ``~/.cache/ionmc``. The directory is *not* created here.
    """
    if override is not None:
        return Path(override).expanduser()
    env = os.environ.get(CACHE_ENV_VAR)
    if env:
        return Path(env).expanduser()
    xdg = os.environ.get("XDG_CACHE_HOME")
    if xdg:
        return Path(xdg).expanduser() / "ionmc"
    return Path.home() / ".cache" / "ionmc"


def dataset_dir(spec: DatasetSpec, root: str | os.PathLike[str] | None = None) -> Path:
    """Directory holding the cached file and manifest of ``spec``."""
    return cache_dir(root) / "datasets" / spec.name / spec.version


def cached_file_path(
    spec: DatasetSpec, root: str | os.PathLike[str] | None = None
) -> Path:
    """Path where the cached data file of ``spec`` lives (may not exist)."""
    return dataset_dir(spec, root) / spec.filename


def manifest_path(
    spec: DatasetSpec, root: str | os.PathLike[str] | None = None
) -> Path:
    """Path of the JSON manifest written when ``spec`` was acquired."""
    return dataset_dir(spec, root) / "manifest.json"


def sha256_of_file(path: str | os.PathLike[str]) -> str:
    """Hex SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_cached(spec: DatasetSpec, root: str | os.PathLike[str] | None = None) -> bool:
    """True if the data file and its manifest are present in the cache."""
    return (
        cached_file_path(spec, root).is_file() and manifest_path(spec, root).is_file()
    )


def read_manifest(
    spec: DatasetSpec, root: str | os.PathLike[str] | None = None
) -> dict[str, Any]:
    """Return the manifest of a cached dataset (raises ``DatasetNotCached``)."""
    path = manifest_path(spec, root)
    if not path.is_file():
        raise DatasetNotCached(_not_cached_message(spec, root))
    with open(path, encoding="utf-8") as handle:
        result: dict[str, Any] = json.load(handle)
        return result


def load_path(
    spec: DatasetSpec,
    root: str | os.PathLike[str] | None = None,
    verify: bool = True,
) -> Path:
    """Offline *use* phase: return the cached file of ``spec``.

    Raises :class:`DatasetNotCached` with an actionable message if the file
    is absent, and :class:`IntegrityError` if ``verify`` is set and the cached
    bytes do not match ``spec.sha256`` (a corrupted or tampered cache entry).
    No network access is attempted.
    """
    path = cached_file_path(spec, root)
    if not path.is_file():
        raise DatasetNotCached(_not_cached_message(spec, root))
    if verify:
        actual = sha256_of_file(path)
        if actual != spec.sha256.lower():
            raise IntegrityError(
                f"cached file {path} has sha256 {actual}, expected {spec.sha256}; "
                "delete it and re-acquire"
            )
    return path


def acquire(
    spec: DatasetSpec,
    root: str | os.PathLike[str] | None = None,
    force: bool = False,
    timeout_seconds: float = 120.0,
) -> Path:
    """*Acquire* phase: download, verify and cache ``spec`` (needs network).

    The download goes to a temporary file in the target directory, its
    SHA-256 is checked against the spec, and only then is it renamed into
    place (atomic on POSIX). A manifest with provenance, license, checksum,
    source URL and retrieval time is written next to it. An already cached,
    intact file is returned without network access unless ``force`` is set.
    """
    target = cached_file_path(spec, root)
    if not force and target.is_file() and manifest_path(spec, root).is_file():
        if sha256_of_file(target) == spec.sha256.lower():
            return target
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(spec.url, headers={"User-Agent": "ionmc-data/1"})
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        payload = response.read()
    actual = hashlib.sha256(payload).hexdigest()
    if actual != spec.sha256.lower():
        raise IntegrityError(
            f"downloaded {spec.url}: sha256 {actual} does not match the "
            f"specification {spec.sha256}; the source may have changed or the "
            "download was corrupted. Nothing was written to the cache."
        )
    fd, tmp_name = tempfile.mkstemp(dir=target.parent, prefix=".partial-")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload)
        os.replace(tmp_name, target)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise
    manifest = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "spec": asdict(spec),
        "retrieved_at_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "size_bytes": len(payload),
        "sha256": actual,
    }
    with open(manifest_path(spec, root), "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return target


def _not_cached_message(spec: DatasetSpec, root: str | os.PathLike[str] | None) -> str:
    where = dataset_dir(spec, root)
    return (
        f"dataset {spec.name!r} version {spec.version!r} is not in the local "
        f"cache ({where}). Acquire it with network access, e.g.\n"
        f"    python -m ionmc.data acquire {spec.name}\n"
        f"or set {CACHE_ENV_VAR} to a cache that already contains it."
    )
