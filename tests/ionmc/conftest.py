"""Fixtures shared by the ionmc tests: cached external datasets."""

from __future__ import annotations

import urllib.error

import pytest

from ionmc.data import MCSQUARE_G4_WATER, MCSQUARE_PSTAR_WATER, cache
from ionmc.data.cache import DatasetSpec


def _cached_or_acquired(spec: DatasetSpec, tmp_root):
    """Return a cache root holding ``spec``: the configured cache if it already
    has it, otherwise a session-temporary cache filled over the network
    (skipping the test when the network is unavailable)."""
    if cache.is_cached(spec):
        return None  # default cache
    try:
        cache.acquire(spec, tmp_root)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        pytest.skip(f"dataset {spec.name} not cached and network unavailable: {exc}")
    return tmp_root


@pytest.fixture(scope="session")
def dataset_cache_root(tmp_path_factory):
    """Cache root (or None for the default cache) containing the water tables."""
    tmp_root = tmp_path_factory.mktemp("ionmc-data-cache")
    roots = {
        spec.name: _cached_or_acquired(spec, tmp_root)
        for spec in (MCSQUARE_PSTAR_WATER, MCSQUARE_G4_WATER)
    }
    return roots


@pytest.fixture(scope="session")
def pstar_cache_root(dataset_cache_root):
    return dataset_cache_root[MCSQUARE_PSTAR_WATER.name]


@pytest.fixture(scope="session")
def g4_cache_root(dataset_cache_root):
    return dataset_cache_root[MCSQUARE_G4_WATER.name]
