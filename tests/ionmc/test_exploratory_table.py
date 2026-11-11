"""V3-005C C6 (E-shape): the forced-isotropic exploratory elastic table (plan Amendment 15
(b)6). Built into the cache at run time from the pinned sources (never committed); the loader
accepts it only when the configuration names it AND sets ``allow_exploratory_table``; no
lv5c/hr5c evidence producer can use it. Data-backed (skip without the cache, fail under
``IONMC_REQUIRE_DATA=1``)."""

from __future__ import annotations

import json
import os
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ionmc.config import PhysicsOptions
from ionmc.data import cache
from ionmc.errors import UnsupportedCombinationError
from ionmc.nuclear import elastic_build as EB
from ionmc.nuclear.elastic_tables import ElasticTable, ElasticTableUnqualifiedError
from ionmc.simulation import Simulation
from tests.ionmc.test_elastic_transport import _cfg
from tests.ionmc.test_elastic_transport import tid as tid  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def built() -> tuple[Any, Any]:
    cdir = cache.resolve_cache_dir(None)
    try:
        exp = EB.build_elastic_proton(cdir, exploratory_forced_isotropic=True)
    except FileNotFoundError:
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail("data for the elastic table are not cached")
        pytest.skip("data not cached")
    return exp, ElasticTable.load(cdir, exp.table_id, allow_exploratory=True)


def _pin() -> str:
    return str(json.loads((ROOT / "src/ionmc/data/elastic_table_pin.json").read_text())["table_id"])


def test_id_is_distinct_deterministic_and_the_pinned_id_is_unchanged(built: Any) -> None:
    exp, tab = built
    assert (
        exp.table_id != _pin()
        and tab.exploratory
        and tab.info["exploratory"]["exploratory_forced_isotropic"] is True
    )
    again = EB.build_elastic_proton(
        cache.resolve_cache_dir(None), exploratory_forced_isotropic=True
    )
    assert again.table_id == exp.table_id
    pinned = ElasticTable.load(None, _pin())
    assert not pinned.exploratory and "exploratory" not in pinned.info


def test_loader_refuses_without_the_flag(built: Any) -> None:
    exp, _ = built
    with pytest.raises(ElasticTableUnqualifiedError, match="exploratory"):
        ElasticTable.load(None, exp.table_id)


def test_only_the_nodes_below_the_domain_differ(built: Any) -> None:
    _, tab = built
    pin = ElasticTable.load(None, _pin())
    a, b = tab.arrays, pin.arrays
    assert np.array_equal(a["grid_e_mev"], b["grid_e_mev"]) and np.array_equal(
        a["target_e_min_mev"], b["target_e_min_mev"]
    )
    assert np.array_equal(a["sigma_barn"][0], b["sigma_barn"][0]) and np.array_equal(
        a["edges_mu"][0], b["edges_mu"][0]
    )
    n_q = a["edges_mu"].shape[2] - 1
    for it in range(1, a["sigma_barn"].shape[0]):
        k_min = int(np.searchsorted(a["grid_e_mev"], a["target_e_min_mev"][it]))
        assert np.array_equal(a["sigma_barn"][it, k_min:], b["sigma_barn"][it, k_min:])
        assert np.array_equal(a["edges_mu"][it, k_min:], b["edges_mu"][it, k_min:])
        assert np.all(a["valid"][it] == 1) and np.all(b["valid"][it, :k_min] == 0)
        assert np.all(
            a["edges_mu"][it, :k_min] == np.linspace(-1.0, 1.0, n_q + 1)
        )  # isotropic in mu_CM
    assert np.any(
        a["sigma_barn"][3, : np.searchsorted(a["grid_e_mev"], a["target_e_min_mev"][3])] > 0
    )  # BGG sigma


def test_configuration_rules_and_the_marking_of_summary_and_result(built: Any, tid: str) -> None:  # noqa: F811
    exp, _ = built
    cfg = _cfg(tid, 100.0, 4, elastic_table_id=exp.table_id)
    with pytest.raises(UnsupportedCombinationError):
        _ = Simulation(cfg).effective  # named but not allowed
    ok = replace(cfg, physics=replace(cfg.physics, allow_exploratory_table=True))
    sim = Simulation(ok)
    assert sim.effective.summary()["elastic"]["exploratory"] is True
    res = sim.run()
    assert res.diagnostics["exploratory"] is True
    normal = Simulation(_cfg(tid, 100.0, 4))
    assert "exploratory" not in normal.effective.summary()["elastic"]
    assert "exploratory" not in normal.run().diagnostics
    with pytest.raises(UnsupportedCombinationError):  # never selected implicitly
        PhysicsOptions(
            stopping=cfg.physics.stopping, nuclear=True, elastic=True, allow_exploratory_table=True
        )
    with pytest.raises(UnsupportedCombinationError):
        replace(ok.physics, elastic_table_id=None)


def test_evidence_builders_and_partials_refuse_exploratory_tables(built: Any) -> None:
    import importlib
    import sys

    sys.path.insert(0, str(ROOT / "validation" / "scripts" / "transport"))
    c = importlib.import_module("steps_v5c")
    exp, _ = built
    import steps_v5 as v5

    from ionmc import materials as M

    geo, grid, _, _ = v5.depth_box(M.WATER, 100.0)
    cfg = c.nuc_config_c(energy=100.0, n=10, seed=1, geometry=geo, grid=grid)
    assert c.require_evidence_config(cfg) is cfg
    bad_id = replace(cfg, physics=replace(cfg.physics, elastic_table_id=exp.table_id))
    with pytest.raises(SystemExit, match="not the pinned"):
        c.require_evidence_config(bad_id)
    bad_flag = replace(bad_id, physics=replace(bad_id.physics, allow_exploratory_table=True))
    with pytest.raises(SystemExit, match="exploratory"):
        c.require_evidence_config(bad_flag)
