"""Tests of proton energy-loss straggling and the Bragg peak (decision 0010).

The clean physics metric is the range straggling ``sigma_R`` (std of the
stopping depths), which isolates straggling from the depth-dose peak shape. It
is compared with the model's own analytic Bohr integral (a sampling-correctness
check) and with Bortfeld's empirical law (a physics check). Large-statistics
runs use the Warp CPU path (fast); the reference path is checked at small N.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ionmc import materials, particles
from ionmc.constants import PROTON_MASS_MEV
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import transport as transport_physics
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab

ENERGIES = (100.0, 150.0, 200.0)


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.fixture(scope="module")
def engine(table):
    return TransportEngine(
        table,
        WaterSlab(400.0, materials.WATER),
        DepthDoseGrid(400.0, 2000),
        straggling=True,
    )


def _analytic_sigma_r_mm(table, e0: float, floor_mev: float = 2.0) -> float:
    """Bohr range straggling sqrt(integral (dOmega^2/dx)/S^3 dE), in mm (rho=1).

    The integral starts at ``floor_mev`` (the MC's straggling floor), so the
    analytic reference and the Monte Carlo accumulate straggling over the same
    energy range.
    """
    tab = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, path="numpy")
    e = np.linspace(floor_mev, e0, 4000)
    tau = e / PROTON_MASS_MEV
    gamma = 1.0 + tau
    beta2 = tau * (tau + 2.0) / (gamma * gamma)
    f_rel = (1.0 - 0.5 * beta2) / (1.0 - beta2)
    za = materials.WATER.electrons_per_gram_ratio
    d_omega2_dx = (
        transport_physics.BOHR_K_MEV2_CM2_PER_MOL * za * 1.0 * f_rel
    )  # MeV^2/cm
    s_lin = tab.mass_stopping_power(e) * 1.0  # MeV/cm
    sigma_r_cm = math.sqrt(np.trapezoid(d_omega2_dx / s_lin**3, e))
    return sigma_r_cm * 10.0


def _csda_range_mm(table, e0: float) -> float:
    tab = TabulatedStoppingPower(
        table, materials.WATER, particles.PROTON, path="python"
    )
    return float(tab.csda_range(e0)[0]) * 10.0


def _bortfeld_sigma_mm(range_mm: float) -> float:
    return 0.012 * (range_mm / 10.0) ** 0.935 * 10.0


def test_bohr_sigma_matches_hand_value() -> None:
    # dOmega^2/dx at 100 MeV in water: 0.0969 MeV^2/cm (research table).
    sig_1cm = transport_physics.bohr_straggling_sigma(
        100.0, PROTON_MASS_MEV, 10.0, 1.0, 0.5551, 1.0
    )  # 10 mm = 1 cm
    assert sig_1cm**2 == pytest.approx(0.0969, rel=0.01)


@pytest.mark.warp
@pytest.mark.parametrize("e0", ENERGIES)
def test_mean_range_is_unbiased(warp_module, engine, table, e0) -> None:
    res = engine.run(PencilBeamSource(e0), n_histories=20000, seed=1, path="warp")
    assert res.n_stopped > 19000
    assert abs(res.range_mean_mm / _csda_range_mm(table, e0) - 1.0) <= 1e-3


@pytest.mark.warp
@pytest.mark.parametrize("e0", ENERGIES)
def test_range_straggling_matches_analytic_and_bortfeld(warp_module, engine, table, e0):
    res = engine.run(PencilBeamSource(e0), n_histories=30000, seed=7, path="warp")
    sigma_mc = res.range_sigma_mm
    sigma_ana = _analytic_sigma_r_mm(table, e0)
    sigma_bort = _bortfeld_sigma_mm(_csda_range_mm(table, e0))
    assert abs(sigma_mc / sigma_ana - 1.0) <= 0.03, (e0, sigma_mc, sigma_ana)
    assert abs(sigma_mc / sigma_bort - 1.0) <= 0.10, (e0, sigma_mc, sigma_bort)
    # sigma_R / R in the 0.9-1.2 % band
    frac = sigma_mc / res.range_mean_mm
    assert 0.009 <= frac <= 0.012, (e0, frac)


@pytest.mark.warp
def test_bragg_peak_exists_and_conserves_energy(warp_module, engine) -> None:
    res = engine.run(PencilBeamSource(150.0), n_histories=20000, seed=3, path="warp")
    centers = res.grid.centers_mm
    peak_bin = int(np.argmax(res.edep_mev))
    entrance = res.edep_mev[centers < 20.0].mean()
    # a Bragg peak: the maximum is deep and well above the entrance plateau
    assert centers[peak_bin] > 140.0
    assert res.edep_mev[peak_bin] > 3.0 * entrance
    assert abs(res.energy_balance) <= 5e-5


def test_reference_path_straggles_and_conserves_energy(engine) -> None:
    # small N on the float64 reference path: energy is conserved exactly and the
    # stopping depths have non-zero spread (straggling is active)
    res = engine.run(PencilBeamSource(150.0), n_histories=64, seed=5, path="python")
    assert abs(res.energy_balance) <= 1e-9
    assert res.range_sigma_mm > 0.5  # mm, straggling present


@pytest.mark.warp
def test_batched_uncertainty_and_statistical_consistency(warp_module, engine) -> None:
    a = engine.run_batched(
        PencilBeamSource(150.0), n_histories=20000, n_batches=10, seed=100, path="warp"
    )
    b = engine.run_batched(
        PencilBeamSource(150.0), n_histories=20000, n_batches=10, seed=500, path="warp"
    )
    assert a.mean_edep_mev.shape == (a.grid.n_bins,)
    # standard error is positive where there is dose
    peak = a.mean_edep_mev.max()
    mask = a.mean_edep_mev > 0.01 * peak
    assert np.all(a.standard_error_mev[mask] > 0.0)
    # two independent estimates are statistically consistent (decision 0010)
    combined = np.sqrt(a.standard_error_mev**2 + b.standard_error_mev**2)
    t = (a.mean_edep_mev[mask] - b.mean_edep_mev[mask]) / combined[mask]
    assert 0.5 <= math.sqrt(np.mean(t**2)) <= 1.6, math.sqrt(np.mean(t**2))
    assert np.max(np.abs(t)) <= 5.0


@pytest.mark.warp
def test_cross_backend_cumulative_with_straggling(warp_module, engine) -> None:
    # shared RNG streams -> reference and Warp agree far better than statistics
    ref = engine.run(PencilBeamSource(150.0), n_histories=200, seed=11, path="python")
    warp = engine.run(PencilBeamSource(150.0), n_histories=200, seed=11, path="warp")
    total = np.sum(ref.edep_mev)
    cum = np.max(np.abs(np.cumsum(warp.edep_mev) - np.cumsum(ref.edep_mev))) / total
    assert cum <= 1e-4, cum


@pytest.mark.cuda
def test_warp_cuda_matches_cpu_with_straggling(warp_module, cuda_available, engine):
    cpu = engine.run(
        PencilBeamSource(150.0), n_histories=20000, seed=9, path="warp", device="cpu"
    )
    cuda = engine.run(
        PencilBeamSource(150.0), n_histories=20000, seed=9, path="warp", device="cuda:0"
    )
    total = np.sum(cpu.edep_mev)
    cum = np.max(np.abs(np.cumsum(cuda.edep_mev) - np.cumsum(cpu.edep_mev))) / total
    assert cum <= 1e-5, cum
