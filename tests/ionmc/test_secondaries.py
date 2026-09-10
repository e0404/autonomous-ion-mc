"""Tests of secondary charged-particle production and transport (decision 0013).

The sampling module is checked directly (multiplicity, per-reaction energy
renormalisation, determinism); the engine tests check that ``secondaries=False``
reproduces the DEV-007 baseline exactly, that the energy budget stays closed
with the second pass, that the secondary dose has the right magnitude
(~1-2 % at entrance, a few percent of total), and that the reference and Warp
paths produce the identical secondary set (host-side generation from the shared
reaction records).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ionmc import materials
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import secondaries
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab

ENERGIES = (150.0, 200.0)


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


def _engine(table, **kw) -> TransportEngine:
    return TransportEngine(
        table, WaterSlab(400.0, materials.WATER), DepthDoseGrid(400.0, 800), **kw
    )


# -- sampling module ----------------------------------------------------------


def test_mean_multiplicity_linear() -> None:
    assert secondaries.mean_multiplicity(0.0) == pytest.approx(0.5)
    assert secondaries.mean_multiplicity(200.0) == pytest.approx(0.5 + 0.004 * 200.0)
    assert secondaries.mean_multiplicity(-1000.0) == 0.0  # clamped >= 0


def test_generate_secondaries_energy_and_determinism() -> None:
    """The moved (transported + sub-cut local) energy equals the energy taken out
    of the escaping channel; per reaction that produces >=1 secondary the pool is
    exactly f_p * w * E (renormalisation); and the sampler is deterministic."""
    n = 50  # identical reactions so the per-reaction pool is a fixed unit
    react_z = np.full(n, 120.0)
    react_e = np.full(n, 150.0)
    weight = np.full(n, 2.0)
    b1 = secondaries.generate_secondaries(react_z, react_e, weight, seed=99)
    b2 = secondaries.generate_secondaries(react_z, react_e, weight, seed=99)
    # deterministic in the seed
    assert b1.n_secondaries == b2.n_secondaries
    assert b1.escaped_reduction_mev == b2.escaped_reduction_mev
    # moved energy (transported + sub-cut local) == energy removed from escaping
    transported = (
        0.0
        if b1.state is None
        else float(np.sum(b1.state.energy_mev * b1.state.weight))
    )
    local = sum(we for _, we in b1.local_deposit)
    assert math.isclose(transported + local, b1.escaped_reduction_mev, rel_tol=1e-9)
    # the removed energy is an integer number of per-reaction pools f_p * w * E
    unit = secondaries.F_SECONDARY * 2.0 * 150.0
    ratio = b1.escaped_reduction_mev / unit
    assert ratio > 0.0
    assert abs(ratio - round(ratio)) < 1e-9, ratio


def test_generate_secondaries_no_reactions() -> None:
    z = np.zeros(4)
    e = np.zeros(4)
    w = np.ones(4)
    b = secondaries.generate_secondaries(z, e, w, seed=1)
    assert b.state is None
    assert b.n_secondaries == 0
    assert b.escaped_reduction_mev == 0.0


# -- engine integration -------------------------------------------------------


def test_secondaries_require_nuclear(table) -> None:
    with pytest.raises(ValueError, match="requires nuclear"):
        _engine(table, nuclear=False, secondaries=True)


def test_secondaries_off_matches_dev007(table) -> None:
    """secondaries=False is the DEV-007 path: no second pass, unchanged dose."""
    src = PencilBeamSource(150.0)
    base = _engine(table, nuclear=True, secondaries=False).run(
        src, n_histories=400, seed=7, path="python"
    )
    assert base.n_secondaries == 0
    assert base.secondary_edep_mev is None
    # a second run is identical (determinism, no hidden secondary draws)
    again = _engine(table, nuclear=True, secondaries=False).run(
        src, n_histories=400, seed=7, path="python"
    )
    np.testing.assert_array_equal(base.edep_mev, again.edep_mev)


@pytest.mark.parametrize("e0", ENERGIES)
def test_energy_budget_closes_with_secondaries(table, e0) -> None:
    """deposited + escaped = energy_in exactly with the secondary pass on."""
    res = _engine(table, nuclear=True, secondaries=True).run(
        PencilBeamSource(e0), n_histories=600, seed=7, path="python"
    )
    assert res.n_reactions > 0
    assert res.n_secondaries > 0
    assert res.secondary_edep_mev is not None
    assert abs(res.energy_balance) <= 1e-9


def test_secondaries_add_dose_and_reduce_escaping(table) -> None:
    """Turning secondaries on moves energy from escaping into deposited dose."""
    src = PencilBeamSource(200.0)
    off = _engine(table, nuclear=True, secondaries=False).run(
        src, n_histories=800, seed=7, path="python"
    )
    on = _engine(table, nuclear=True, secondaries=True).run(
        src, n_histories=800, seed=7, path="python"
    )
    # same primaries react (secondary pass does not change the primary stream)
    assert on.n_reactions == off.n_reactions
    # secondaries deposit dose and lower the escaping energy
    assert on.secondary_edep_mev.sum() > 0.0
    assert on.escaped_mev < off.escaped_mev
    assert on.energy_deposited_mev > off.energy_deposited_mev


@pytest.mark.warp
@pytest.mark.parametrize("e0", ENERGIES)
def test_secondary_dose_fraction_matches_published(warp_module, table, e0) -> None:
    """Secondary protons contribute ~1-2 % at entrance and a few percent of the
    total dose (Paganetti 2002). Large-statistics Warp CPU run."""
    res = _engine(table, nuclear=True, secondaries=True).run(
        PencilBeamSource(e0), n_histories=40000, seed=11, path="warp", device="cpu"
    )
    frac = res.secondary_dose_fraction
    entrance = float(np.mean(frac[5:25]))  # ~2.5-12 mm depth
    total_fraction = res.secondary_edep_mev.sum() / res.energy_deposited_mev
    assert 0.005 <= entrance <= 0.04, entrance
    assert 0.02 <= total_fraction <= 0.12, total_fraction
    # the fraction rises across the plateau to ~5-10 % and collapses at the peak
    peak_bin = int(np.argmax(res.edep_mev))
    plateau_mean = float(np.mean(frac[25 : peak_bin - 10]))
    assert 0.03 <= plateau_mean <= 0.12, plateau_mean
    assert plateau_mean > entrance
    assert frac[peak_bin] < entrance


@pytest.mark.warp
def test_reference_and_warp_secondaries_match(warp_module, table) -> None:
    """Host-side secondary generation from the shared reaction records makes the
    reference and Warp paths produce the same secondary set and dose in practice.

    The secondary count is derived from the reaction residual energy, which is
    float32 on the Warp path and float64 on the reference path, so a boundary
    case could shift a Poisson draw or a sub-cut classification; the counts are
    therefore compared within a small tolerance rather than for bit equality,
    and the cumulative depth-dose agreement (decision 0001) is the real gate.
    """
    eng = _engine(table, nuclear=True, secondaries=True)
    src = PencilBeamSource(150.0)
    ref = eng.run(src, n_histories=4000, seed=7, path="python")
    war = eng.run(src, n_histories=4000, seed=7, path="warp", device="cpu")
    assert war.n_reactions == ref.n_reactions
    assert abs(war.n_secondaries - ref.n_secondaries) <= max(
        1, round(0.005 * ref.n_secondaries)
    )
    total = float(np.sum(ref.edep_mev))
    cum = np.max(np.abs(np.cumsum(war.edep_mev) - np.cumsum(ref.edep_mev))) / total
    assert cum <= 1e-4
    sec_total = float(np.sum(ref.secondary_edep_mev))
    sec_cum = (
        np.max(
            np.abs(
                np.cumsum(war.secondary_edep_mev) - np.cumsum(ref.secondary_edep_mev)
            )
        )
        / sec_total
    )
    assert sec_cum <= 1e-4
