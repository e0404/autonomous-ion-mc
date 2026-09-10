"""Tests of multiple Coulomb scattering and lateral spread (decision 0011).

The primary metric is the lateral spread sigma_x(z), validated against the
Fermi-Eyges oracle built from the same scattering power (a tight
sampling/hinge-correctness check) and against published values. Large-N runs
use the Warp CPU path; the reference path is checked at small N.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials, particles
from ionmc.constants import PROTON_MASS_MEV
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import transport as tphys
from ionmc.physics.fermi_eyges import lateral_sigma_x_mm
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import (
    DepthDoseGrid,
    DepthLateralGrid,
    PencilBeamSource,
    TransportEngine,
    WaterSlab,
)

ENERGIES = (150.0, 200.0)
PUBLISHED_SIGMA_08R = {150.0: 2.4, 200.0: 3.9}  # mm, Gottschalk / Fermi-Eyges
X0 = materials.WATER.radiation_length_g_per_cm2


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.fixture(scope="module")
def engine(table):
    return TransportEngine(
        table, WaterSlab(400.0, materials.WATER), DepthDoseGrid(400.0, 10)
    )


@pytest.fixture(scope="module")
def model(table):
    return TabulatedStoppingPower(
        table, materials.WATER, particles.PROTON, path="numpy"
    )


def _range_mm(model, e0: float) -> float:
    return float(model.csda_range(e0)[0]) * 10.0


def _grid(range_mm: float) -> DepthLateralGrid:
    return DepthLateralGrid(
        depth_mm=range_mm * 1.05,
        n_depth=int(range_mm * 1.05 / 0.5),
        half_width_mm=25.0,
        n_lateral=500,
    )


def test_theta0_and_pv_match_highland_hand_values() -> None:
    assert tphys.momentum_times_velocity(100.0, PROTON_MASS_MEV) == pytest.approx(
        190.369, rel=1e-4
    )
    assert tphys.momentum_times_velocity(200.0, PROTON_MASS_MEV) == pytest.approx(
        364.859, rel=1e-4
    )
    # theta0 over 1 mm water, bracket = 1 (scattering-power form)
    th100 = tphys.highland_theta0(
        100.0,
        PROTON_MASS_MEV,
        1.0,
        1.0,
        1.0,
        materials.WATER.radiation_length_g_per_cm2,
    )
    assert th100 * 1000.0 == pytest.approx(3.761, rel=1e-3)


def test_scattering_requires_radiation_length(table) -> None:
    # a material with no radiation length (default 0.0) cannot scatter
    bare = materials.Material(
        "bare_water",
        1.0,
        {"H": 0.111894, "O": 0.888106},
        materials.MeanExcitationEnergy(75.0, "x"),
    )
    assert bare.radiation_length_g_per_cm2 == 0.0
    eng = TransportEngine(table, WaterSlab(400.0, bare), DepthDoseGrid(400.0, 10))
    with pytest.raises(ValueError):
        eng.run_scattering(PencilBeamSource(150.0), _grid(150.0), 1, path="python")
    # and with_mean_excitation_energy now PRESERVES the radiation length
    icru90 = materials.WATER.with_mean_excitation_energy(78.0, "ICRU90")
    assert (
        icru90.radiation_length_g_per_cm2 == materials.WATER.radiation_length_g_per_cm2
    )


@pytest.mark.warp
@pytest.mark.parametrize("e0", ENERGIES)
def test_lateral_sigma_matches_fermi_eyges(warp_module, engine, model, e0) -> None:
    r = _range_mm(model, e0)
    grid = _grid(r)
    res = engine.run_scattering(PencilBeamSource(e0), grid, 30000, seed=42, path="warp")
    depths = np.array([0.5 * r, 0.8 * r])
    fe = lateral_sigma_x_mm(
        model, e0, depths, materials.WATER.radiation_length_g_per_cm2
    )
    mc = np.array([res.sigma_x_at_depth(d) for d in depths])
    assert np.all(np.abs(mc / fe - 1.0) <= 0.03), (e0, mc, fe)
    # secondary: sigma_x(0.8R) vs published
    assert abs(mc[1] / PUBLISHED_SIGMA_08R[e0] - 1.0) <= 0.08, (e0, mc[1])
    assert abs(res.energy_balance) <= 5e-5


def test_reference_scattering_conserves_energy_and_spreads(engine, model) -> None:
    r = _range_mm(model, 150.0)
    res = engine.run_scattering(
        PencilBeamSource(150.0), _grid(r), 64, seed=5, path="python"
    )
    assert abs(res.energy_balance) <= 1e-9
    # lateral spread grows with depth
    s_shallow = res.sigma_x_at_depth(0.3 * r)
    s_deep = res.sigma_x_at_depth(0.8 * r)
    assert 0.0 < s_shallow < s_deep


@pytest.mark.warp
def test_depth_dose_marginal_has_a_bragg_peak(warp_module, engine, model) -> None:
    r = _range_mm(model, 150.0)
    grid = _grid(r)
    res = engine.run_scattering(
        PencilBeamSource(150.0), grid, 20000, seed=3, path="warp"
    )
    dd = res.depth_dose_mev
    centers = grid.depth_centers_mm
    peak_bin = int(np.argmax(dd))
    entrance = dd[centers < 20.0].mean()
    assert centers[peak_bin] > 0.9 * r  # peak near end of range
    assert dd[peak_bin] > 3.0 * entrance
    # detour: projected peak is at or slightly before the CSDA range (<1%)
    assert centers[peak_bin] <= r * 1.001


@pytest.mark.warp
def test_cross_backend_sigma_x_consistent(warp_module, engine, model) -> None:
    # reference (small N) and Warp (same seed/streams) agree on sigma_x
    r = _range_mm(model, 150.0)
    grid = _grid(r)
    ref = engine.run_scattering(
        PencilBeamSource(150.0), grid, 400, seed=11, path="python"
    )
    warp = engine.run_scattering(
        PencilBeamSource(150.0), grid, 400, seed=11, path="warp"
    )
    depth = 0.8 * r
    assert abs(ref.sigma_x_at_depth(depth) - warp.sigma_x_at_depth(depth)) <= 0.1  # mm


@pytest.mark.cuda
def test_warp_cuda_scattering_matches_cpu(warp_module, cuda_available, engine, model):
    r = _range_mm(model, 150.0)
    grid = _grid(r)
    cpu = engine.run_scattering(
        PencilBeamSource(150.0), grid, 20000, seed=9, path="warp", device="cpu"
    )
    cuda = engine.run_scattering(
        PencilBeamSource(150.0), grid, 20000, seed=9, path="warp", device="cuda:0"
    )
    d = 0.8 * r
    assert abs(cpu.sigma_x_at_depth(d) - cuda.sigma_x_at_depth(d)) <= 0.02  # mm


@pytest.mark.warp
@pytest.mark.parametrize("e0", ENERGIES)
def test_detour_factor_small(warp_module, engine, model, e0) -> None:
    # multiple scattering shortens the mean PROJECTED stopping depth relative to
    # the CSDA path range by the detour factor (< 0.1 %); it must not be larger.
    r = _range_mm(model, e0)
    grid = _grid(r)
    res = engine.run_scattering(PencilBeamSource(e0), grid, 20000, seed=8, path="warp")
    assert res.n_stopped > 19000
    detour = res.range_mean_mm / r - 1.0
    assert -0.003 < detour <= 1e-4, (e0, res.range_mean_mm, r, detour)
