"""External physics-data layer: versioned, checksummed, locally cached datasets.

The layer separates two phases (decision ``0007``):

* **acquire** (needs network): :func:`acquire` downloads a dataset described
  by an immutable :class:`DatasetSpec`, verifies its SHA-256 checksum and
  writes it atomically into the local cache together with a JSON manifest
  recording provenance and retrieval time;
* **use** (offline): :func:`load_path` returns the cached file for a spec, or
  raises :class:`DatasetNotCached` with the exact command that acquires it.

The cache location is configurable (:func:`cache_dir`). Datasets are
identified by name *and* immutable version (a Git commit or release tag), and
the checksum is part of the spec, so a cached file is only ever accepted if
it is byte-identical to what the spec describes.
"""

from __future__ import annotations

from ionmc.data.cache import (
    DatasetNotCached,
    DatasetSpec,
    IntegrityError,
    acquire,
    cache_dir,
    is_cached,
    load_path,
    manifest_path,
    read_manifest,
    sha256_of_file,
)
from ionmc.data.registry import DATASETS, MCSQUARE_G4_WATER, MCSQUARE_PSTAR_WATER

__all__ = [
    "DATASETS",
    "MCSQUARE_G4_WATER",
    "MCSQUARE_PSTAR_WATER",
    "DatasetNotCached",
    "DatasetSpec",
    "IntegrityError",
    "acquire",
    "cache_dir",
    "is_cached",
    "load_path",
    "manifest_path",
    "read_manifest",
    "sha256_of_file",
]
