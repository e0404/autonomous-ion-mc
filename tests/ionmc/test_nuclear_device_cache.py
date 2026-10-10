"""The packed nuclear device is built once per process and key (V3-005B C37).

The V7 shard crash of the 2026-10-08 host runs: ``pack_nuclear`` re-run for every ``run_range``
call, 18000 times per shard.

Data-backed (``IONMC_CACHE_DIR``; skipped without a built nuclear table, a failure with
``IONMC_REQUIRE_DATA=1``), warp-cpu float64 / float32, 64 histories."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import ionmc.transport.nuclear_device as nd
from ionmc.simulation import Simulation
from ionmc.transport.warp_driver import run_warp_range
from tests.ionmc.test_nuclear_transport_loop import _config, _table_id

N = 64
REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def tid() -> str:
    return _table_id()


@pytest.fixture(autouse=True)
def _fresh_cache() -> Any:
    nd.clear_nuclear_device_cache()
    yield
    nd.clear_nuclear_device_cache()


def _eff(tid: str, precision: str = "float64", nuclear: bool = True) -> Any:
    cfg = _config(tid, energy=150.0, n=N, track_end=False)
    phys = cfg.physics if nuclear else replace(cfg.physics, nuclear=False, nuclear_table_id=None)
    cfg = replace(cfg, physics=phys, run=replace(cfg.run, backend="warp-cpu", precision=precision))
    return Simulation(cfg).effective


def _digest(part: Any) -> str:
    h = hashlib.sha256()
    for g in part.edep:
        h.update(np.ascontiguousarray(g).tobytes())
    for c in part.tally_components:
        h.update(np.ascontiguousarray(c).tobytes())
    h.update(np.asarray(part.counter_sums, dtype=np.int64).tobytes())
    if part.channel_acc is not None:
        h.update(np.ascontiguousarray(part.channel_acc).tobytes())
    return h.hexdigest()


def test_successive_runs_are_identical_and_pack_once(
    tid: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = {"pack": 0, "ame": 0}
    real_pack, real_ame = nd.pack_nuclear, nd.cached_ame2020

    def counting_pack(*a: Any, **k: Any) -> Any:
        calls["pack"] += 1
        return real_pack(*a, **k)

    def counting_ame() -> Any:
        calls["ame"] += 1
        return real_ame()

    monkeypatch.setattr(nd, "pack_nuclear", counting_pack)
    monkeypatch.setattr(nd, "cached_ame2020", counting_ame)
    eff = _eff(tid)
    p1 = run_warp_range(eff, 0, N, "cpu")
    p2 = run_warp_range(eff, 0, N, "cpu")
    p3 = run_warp_range(eff, 0, N // 2, "cpu")  # a different history range reuses the same device
    assert calls == {"pack": 1, "ame": 1}  # packed once, AME2020 parsed once
    assert _digest(p1) == _digest(p2)
    assert (
        p1.meta["nuclear_device_sha256"]
        == p2.meta["nuclear_device_sha256"]
        == p3.meta["nuclear_device_sha256"]
    )
    assert list(p1.counter_sums) == list(p2.counter_sums) and sum(p1.counter_sums) == 0
    assert len(nd._DEVICE_CACHE) == 1
    # the cached device is the object of the first call and its bytes are unchanged by the runs
    (dev,) = nd._DEVICE_CACHE.values()
    assert dev.device_sha256() == dev.sha256 == p1.meta["nuclear_device_sha256"]


def test_key_and_invalidation_equal_a_fresh_build(tid: str) -> None:
    eff64, eff32 = _eff(tid, "float64"), _eff(tid, "float32")
    t = eff64.nuclear.table
    k64 = nd.nuclear_device_key(t, eff64.geometry.materials, nd.DEFAULT_F_E, nd.wp.float64, "cpu")
    k32 = nd.nuclear_device_key(t, eff64.geometry.materials, nd.DEFAULT_F_E, nd.wp.float32, "cpu")
    assert k64 != k32 and k64[0] == t.table_id and k64[1] == t.info["npz_sha256"]
    assert k64[2] == tuple(m.name for m in eff64.geometry.materials) and k64[3] == nd.DEFAULT_F_E
    assert k64[4] == "float64" and k32[4] == "float32" and k64[5] == "cpu"
    first = run_warp_range(eff64, 0, N, "cpu")
    other = run_warp_range(eff32, 0, N, "cpu")  # a second key: both kept (at most a handful)
    assert len(nd._DEVICE_CACHE) == 2
    assert first.meta["nuclear_device_sha256"] != other.meta["nuclear_device_sha256"]
    nd.clear_nuclear_device_cache()  # invalidation: the next run packs again
    again = run_warp_range(eff64, 0, N, "cpu")
    assert _digest(again) == _digest(first)
    assert again.meta["nuclear_device_sha256"] == first.meta["nuclear_device_sha256"]
    # the cache never grows beyond a handful of entries
    for i in range(nd._MAX_CACHED_DEVICES + 3):
        nd.cached_nuclear_device(
            t,
            eff64.geometry.materials,
            real=nd.wp.float64,
            device="cpu",
            f_e=nd.DEFAULT_F_E + 1e-3 * (i + 1),
        )
    assert len(nd._DEVICE_CACHE) <= nd._MAX_CACHED_DEVICES


def test_cached_run_equals_a_fresh_process_run(tid: str) -> None:
    eff = _eff(tid)
    mine = run_warp_range(eff, 0, N, "cpu")
    mine2 = run_warp_range(eff, 0, N, "cpu")  # served from the cache
    code = (
        "import hashlib, json, sys\n"
        "import numpy as np\n"
        "from dataclasses import replace\n"
        "from ionmc.simulation import Simulation\n"
        "from ionmc.transport.warp_driver import run_warp_range\n"
        "from tests.ionmc.test_nuclear_transport_loop import _config\n"
        f"cfg = _config({tid!r}, energy=150.0, n={N}, track_end=False)\n"
        "cfg = replace(cfg, run=replace(cfg.run, backend='warp-cpu', precision='float64'))\n"
        "eff = Simulation(cfg).effective\n"
        f"p = run_warp_range(eff, 0, {N}, 'cpu')\n"
        "from tests.ionmc.test_nuclear_device_cache import _digest\n"
        "print(json.dumps({'digest': _digest(p), 'sha': p.meta['nuclear_device_sha256']}))\n"
    )
    env = {
        **os.environ,
        "PYTHONPATH": f"{REPO / 'src'}{os.pathsep}{REPO}",
        "PYTHONFAULTHANDLER": "1",
    }
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, env=env, cwd=REPO, timeout=900
    )
    assert out.returncode == 0, out.stderr[-2000:]
    fresh = json.loads(out.stdout.strip().splitlines()[-1])
    assert fresh["digest"] == _digest(mine) == _digest(mine2)
    assert fresh["sha"] == mine.meta["nuclear_device_sha256"]


def test_em_only_path_does_not_touch_the_nuclear_cache(tid: str) -> None:
    """The EM-only transport (A16 nuclear=False digest path) never builds or uses a nuclear
    device."""
    eff = _eff(tid, nuclear=False)
    assert eff.nuclear is None
    part = run_warp_range(eff, 0, N, "cpu")
    assert nd._DEVICE_CACHE == {} and nd._AME_CACHE == {}
    assert "nuclear_device_sha256" not in part.meta
