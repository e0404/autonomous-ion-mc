"""Download external datasets into the verified cache (https only, hash pinned)."""

from __future__ import annotations

import os
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from ionmc.data.cache import (
    IntegrityError,
    read_manifest,
    resolve_cache_dir,
    sha256_bytes,
    store_object,
    verify,
    write_manifest,
)
from ionmc.data.registry import DATASETS, Dataset

TIMEOUT_S = 30.0
OFFLINE_ENV = "IONMC_OFFLINE"


class OfflineError(RuntimeError):
    """Raised when a dataset is not cached and network access is disabled."""


class _HttpsOnlyRedirect(urllib.request.HTTPRedirectHandler):
    """Follow redirects only to https URLs."""

    def redirect_request(  # type: ignore[no-untyped-def]
        self, req, fp, code, msg, headers, newurl
    ):
        if not newurl.lower().startswith("https://"):
            raise urllib.error.URLError(f"refusing non-https redirect to {newurl!r}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _download(dataset: Dataset) -> bytes:
    """Download the dataset bytes (https only, 30 s timeout)."""
    if not dataset.url.lower().startswith("https://"):
        raise ValueError(f"only https URLs are allowed: {dataset.url}")
    body = dataset.post_body.encode("ascii") if dataset.post_body is not None else None
    headers = {"User-Agent": "ionmc-data/1"}
    if body is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    request = urllib.request.Request(dataset.url, data=body, headers=headers, method=dataset.method)
    opener = urllib.request.build_opener(_HttpsOnlyRedirect)
    with opener.open(request, timeout=TIMEOUT_S) as response:  # noqa: S310
        data: bytes = response.read()
    return data


def is_offline(offline: bool = False) -> bool:
    """True if ``offline`` is set or the environment variable ``IONMC_OFFLINE`` is ``1``."""
    return offline or os.environ.get(OFFLINE_ENV) == "1"


def fetch(
    dataset_id: str,
    cache_dir: str | os.PathLike[str] | None = None,
    offline: bool = False,
) -> Path:
    """Return the path of the verified cached dataset, downloading it if needed.

    A cached object is re-hashed (no network). Otherwise, unless offline, the dataset
    is downloaded, its SHA-256 and size are checked against the registry before anything
    is stored (:class:`IntegrityError` on mismatch), and a manifest is written.
    Raises :class:`OfflineError` when offline and not cached, ``KeyError`` for unknown ids.
    """
    if dataset_id not in DATASETS:
        raise KeyError(f"unknown dataset {dataset_id!r}; known: {sorted(DATASETS)}")
    dataset = DATASETS[dataset_id]
    cdir = resolve_cache_dir(cache_dir)
    manifest = read_manifest(dataset_id, cdir)
    if manifest is not None:
        return verify(dataset_id, cdir)
    if is_offline(offline):
        raise OfflineError(f"dataset {dataset_id!r} is not cached in {cdir} and offline mode is on")
    data = _download(dataset)
    digest = sha256_bytes(data)
    if digest != dataset.sha256 or len(data) != dataset.bytes:
        raise IntegrityError(
            f"{dataset_id}: downloaded sha256 {digest} ({len(data)} bytes) does not match "
            f"pinned {dataset.sha256} ({dataset.bytes} bytes); nothing stored"
        )
    path = store_object(data, cdir)
    record: dict[str, Any] = {
        "dataset_id": dataset.id,
        "version": dataset.version,
        "url": dataset.url,
        "method": dataset.method,
        "post_body": dataset.post_body,
        "sha256": digest,
        "bytes": len(data),
        "retrieved_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "license": dataset.license,
        "citation": dataset.citation,
    }
    write_manifest(record, cdir)
    return path
