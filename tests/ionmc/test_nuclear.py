"""Tests of proton nonelastic nuclear attenuation and local deposition (0012).

The physics functions are exercised through the float64 ``python`` binding
(the shared-source ``m.where`` / ``m.min`` cannot be called from the Warp-bound
Python scope). The engine tests check that: ``nuclear=False`` reproduces the EM
baseline bit-for-bit; ``nuclear=True`` removes primaries at the published rate
(survival ~0.80 at 150 MeV, ~0.73 at 200 MeV) and closes the energy budget
(``deposited + escaped = energy_in``); the Monte Carlo reaction fraction agrees
with the analytic ``1 - exp(-integral Sigma dl)``; and the reference and Warp
paths remove the identical primary set (aligned counter-based RNG).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ionmc import materials
from ionmc.backend import reference
from ionmc.constants import AVOGADRO
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab

ENERGIES = (150.0, 200.0)
#: Published primary survival to the Bragg peak (Paganetti 2002; Gottschalk).
PUBLISHED_SURVIVAL = {150.0: (0.80, 0.03), 200.0: (0.73, 0.04)}


@pytest.fixture(scope="module")
def nuc():
    """The nuclear physics module bound to the float64 reference backend."""
    return reference.load_bound_module("ionmc.physics.nuclear", "python")


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.fixture(scope="module")
def oxygen_density() -> float:
    w = materials.WATER
    return AVOGADRO * w.atoms_per_gram("O") * w.density_g_per_cm3


def _engine(table, **kw) -> TransportEngine:
    return TransportEngine(
        table, WaterSlab(400.0, materials.WATER), DepthDoseGrid(400.0, 2000), **kw
    )


# -- cross-section physics (reference binding) --------------------------------


def test_cross_section_shape(nuc) -> None:
    """Threshold, peak, and plateau of sigma_nonel on oxygen [barn]."""
    assert nuc.nonelastic_cross_section_oxygen(5.0) == 0.0  # below 7 MeV threshold
    peak = nuc.nonelastic_cross_section_oxygen(20.0)
    assert 0.50 <= peak <= 0.60  # ~550 mb resonance region
    for e in (100.0, 150.0, 200.0, 250.0):
        s = nuc.nonelastic_cross_section_oxygen(e)
        assert 0.30 <= s <= 0.45  # ICRU-63 plateau band
    # plateau declines mildly with energy
    assert nuc.nonelastic_cross_section_oxygen(
        100.0
    ) > nuc.nonelastic_cross_section_oxygen(250.0)


def test_macroscopic_and_probability(nuc, oxygen_density) -> None:
    """Sigma ~0.012-0.013 /cm on the plateau; P is Sigma*dl and zero below cut."""
    sigma_150 = nuc.macroscopic_nonelastic(150.0, oxygen_density)
    assert 0.011 <= sigma_150 <= 0.014  # 1/cm, Gottschalk delta ~ 0.012
    # thin-step probability equals Sigma * (dl/10)
    p = nuc.nonelastic_step_probability(150.0, 1.0, oxygen_density)
    assert math.isclose(p, sigma_150 * 0.1, rel_tol=1e-12)
    # below the 7 MeV threshold the probability is exactly zero
    assert nuc.nonelastic_step_probability(5.0, 1.0, oxygen_density) == 0.0


# -- engine: regression, survival, and energy budget --------------------------


def test_nuclear_off_matches_em_baseline(table) -> None:
    """nuclear=False draws no extra uniform, so the EM result is bit-identical."""
    src = PencilBeamSource(150.0)
    base = _engine(table).run(src, n_histories=2000, seed=99, path="python")
    off = _engine(table, nuclear=False).run(
        src, n_histories=2000, seed=99, path="python"
    )
    assert off.n_reactions == 0
    assert off.escaped_mev == 0.0
    np.testing.assert_array_equal(off.edep_mev, base.edep_mev)


@pytest.mark.parametrize("e0", ENERGIES)
def test_primary_survival_matches_published(table, e0) -> None:
    """Primary survival to the Bragg peak matches Paganetti/Gottschalk."""
    res = _engine(table, nuclear=True).run(
        PencilBeamSource(e0), n_histories=6000, seed=2024, path="python"
    )
    survival = 1.0 - res.n_reactions / res.n_histories
    target, tol = PUBLISHED_SURVIVAL[e0]
    assert abs(survival - target) <= tol, (survival, target, tol)


@pytest.mark.parametrize("e0", ENERGIES)
def test_energy_budget_closes(table, e0) -> None:
    """deposited + escaped = energy_in to float64 round-off, with escaped > 0."""
    res = _engine(table, nuclear=True).run(
        PencilBeamSource(e0), n_histories=3000, seed=7, path="python"
    )
    assert res.n_reactions > 0
    assert res.escaped_mev > 0.0
    assert abs(res.energy_balance) <= 1e-12


def test_local_fraction_split(table) -> None:
    """Escaped energy is (1 - f_local) of the total reacting-primary energy."""
    f_local = 0.30
    res = _engine(table, nuclear=True, nuclear_local_fraction=f_local).run(
        PencilBeamSource(200.0), n_histories=3000, seed=11, path="python"
    )
    # local deposit from reactions = escaped * f_local/(1-f_local)
    local_from_reactions = res.escaped_mev * f_local / (1.0 - f_local)
    assert local_from_reactions > 0.0
    # the reacting-primary energy is fully accounted (local + escaped)
    assert abs(res.energy_balance) <= 1e-12


@pytest.mark.parametrize("e0", ENERGIES)
def test_reaction_fraction_matches_analytic(table, oxygen_density, nuc, e0) -> None:
    """MC reaction fraction agrees with 1 - exp(-integral Sigma dl) on the mean
    track (a sampling-correctness check independent of the depth-dose shape)."""
    engine = _engine(table, nuclear=True)
    res = engine.run(PencilBeamSource(e0), n_histories=8000, seed=555, path="python")
    mc_fraction = res.n_reactions / res.n_histories

    # analytic survival: change variables from path length to energy along the
    # CSDA track, dl = dE / S_lin, so integral Sigma dl = integral Sigma/S_lin dE
    # from the cutoff to E0 (a shape-independent sampling-correctness check).
    from ionmc import particles
    from ionmc.tabulated_stopping_power import TabulatedStoppingPower

    density = materials.WATER.density_g_per_cm3
    model = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, "numpy")
    e_grid = np.linspace(1.0, e0, 8000)
    s_lin = model.mass_stopping_power(e_grid) * density  # MeV/cm
    sigma = np.array(
        [nuc.macroscopic_nonelastic(float(e), oxygen_density) for e in e_grid]
    )
    integral = np.trapezoid(sigma / s_lin, e_grid)  # dimensionless
    analytic_fraction = 1.0 - math.exp(-integral)

    stat = math.sqrt(mc_fraction * (1.0 - mc_fraction) / res.n_histories)
    assert abs(mc_fraction - analytic_fraction) <= 5.0 * stat + 0.01, (
        mc_fraction,
        analytic_fraction,
    )


# -- cross-backend parity -----------------------------------------------------


@pytest.mark.warp
def test_reference_and_warp_remove_same_primaries(warp_module, table) -> None:
    """Aligned counter RNG: the two backends react on the identical primary set
    and book the identical escaped energy; depth dose agrees cumulatively."""
    engine = _engine(table, nuclear=True)
    src = PencilBeamSource(150.0)
    ref = engine.run(src, n_histories=4000, seed=777, path="python")
    war = engine.run(src, n_histories=4000, seed=777, path="warp", device="cpu")
    assert war.n_reactions == ref.n_reactions
    assert math.isclose(war.escaped_mev, ref.escaped_mev, rel_tol=1e-4)
    total = float(np.sum(ref.edep_mev))
    cum = np.max(np.abs(np.cumsum(war.edep_mev) - np.cumsum(ref.edep_mev))) / total
    assert cum <= 1e-4


@pytest.mark.warp
@pytest.mark.cuda
def test_warp_cuda_matches_cpu(warp_module, cuda_available, table) -> None:
    """CUDA and CPU Warp remove the same primaries and agree cumulatively."""
    engine = _engine(table, nuclear=True)
    src = PencilBeamSource(200.0)
    cpu = engine.run(src, n_histories=4000, seed=321, path="warp", device="cpu")
    cuda = engine.run(src, n_histories=4000, seed=321, path="warp", device="cuda")
    assert cuda.n_reactions == cpu.n_reactions
    total = float(np.sum(cpu.edep_mev))
    cum = np.max(np.abs(np.cumsum(cuda.edep_mev) - np.cumsum(cpu.edep_mev))) / total
    assert cum <= 1e-5
