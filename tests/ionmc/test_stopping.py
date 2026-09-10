"""Tests of the analytical stopping-power model on the reference paths.

Acceptance tolerances are those fixed in decision 0006 *before* the
comparison was run: 1.0 percent against the PSTAR-derived water table for
proton energies of 10 MeV and above. Below 10 MeV the analytic layer is
documented as degraded; the bounds asserted there are regression guards, not
accuracy claims.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ionmc import materials, particles
from ionmc.backend import reference
from ionmc.reference_data.pstar_water import PSTAR_WATER_STOPPING_POWER
from ionmc.stopping_power import (
    AnalyticStoppingPower,
    Corrections,
    build_parameters,
)

TOL_ABOVE_10_MEV = 0.010  # decision 0006
REGRESSION_GUARD_BELOW_10_MEV = {1.0: 0.10, 2.0: 0.05, 5.0: 0.02}


@pytest.fixture(scope="module")
def water_proton_python() -> AnalyticStoppingPower:
    return AnalyticStoppingPower(materials.WATER, particles.PROTON, path="python")


@pytest.fixture(scope="module")
def water_proton_numpy() -> AnalyticStoppingPower:
    return AnalyticStoppingPower(materials.WATER, particles.PROTON, path="numpy")


def test_matches_pstar_within_one_percent_above_10_mev(water_proton_python) -> None:
    energies = sorted(e for e in PSTAR_WATER_STOPPING_POWER if e >= 10.0)
    s = water_proton_python.mass_stopping_power(energies)
    for e, value in zip(energies, s, strict=True):
        ref = PSTAR_WATER_STOPPING_POWER[e]
        assert abs(value / ref - 1.0) <= TOL_ABOVE_10_MEV, (e, value, ref)


def test_regression_guard_below_10_mev(water_proton_python) -> None:
    for e, bound in REGRESSION_GUARD_BELOW_10_MEV.items():
        value = float(water_proton_python.mass_stopping_power(e)[0])
        ref = PSTAR_WATER_STOPPING_POWER[e]
        assert abs(value / ref - 1.0) <= bound, (e, value, ref)


def test_numpy_path_is_bitwise_identical_to_python_path(
    water_proton_python, water_proton_numpy
) -> None:
    # The grid extends to 1500 MeV so that the density-effect ``mid`` branch
    # (x >= x0, i.e. protons above about 705 MeV in water), the only place
    # where ``m.pow`` is evaluated, is exercised as well.
    energies = np.concatenate(
        [np.array(sorted(PSTAR_WATER_STOPPING_POWER)), np.linspace(500.0, 1500.0, 21)]
    )
    a = water_proton_python.mass_stopping_power(energies)
    b = water_proton_numpy.mass_stopping_power(energies)
    np.testing.assert_allclose(b, a, rtol=1e-14, atol=0.0)
    ra = water_proton_python.csda_range([50.0, 100.0])
    rb = water_proton_numpy.csda_range([50.0, 100.0])
    np.testing.assert_allclose(rb, ra, rtol=1e-14, atol=0.0)


def test_csda_ranges_of_protons_in_water(water_proton_numpy) -> None:
    # PSTAR CSDA ranges (g/cm^2) as recalled from the NIST tables by two
    # independent research passes; not fetched in-session, hence treated as
    # informational with a 0.5 percent guard (decision 0006). The analytic
    # integration starts at 1 MeV and neglects ~0.0025 g/cm^2 below it.
    recalled = {100.0: 7.718, 150.0: 15.77, 200.0: 25.96, 250.0: 37.94}
    r = water_proton_numpy.csda_range(sorted(recalled))
    for value, (e, ref) in zip(r, sorted(recalled.items()), strict=True):
        assert abs(value / ref - 1.0) <= 0.005, (e, value, ref)
    # Range difference over the model's validity window: tighter check.
    assert abs((r[2] - r[0]) / (25.96 - 7.718) - 1.0) <= 0.002


def test_csda_range_quadrature_converged(water_proton_numpy) -> None:
    r200 = water_proton_numpy.csda_range([200.0], n_steps=200)[0]
    r2000 = water_proton_numpy.csda_range([200.0], n_steps=2000)[0]
    assert abs(r200 / r2000 - 1.0) < 1e-6


def test_csda_range_rejects_bad_arguments(water_proton_numpy) -> None:
    with pytest.raises(ValueError):
        water_proton_numpy.csda_range([100.0], n_steps=201)
    with pytest.raises(ValueError):
        water_proton_numpy.csda_range([0.5])


def test_correction_magnitudes_are_physical() -> None:
    e = np.array([10.0, 100.0, 250.0])
    full = AnalyticStoppingPower(materials.WATER, particles.PROTON, path="numpy")
    s_full = full.mass_stopping_power(e)

    def without(**kw):
        return AnalyticStoppingPower(
            materials.WATER, particles.PROTON, Corrections(**kw), path="numpy"
        ).mass_stopping_power(e)

    # Shell correction reduces S; ~0.8% at 10 MeV, <0.1% at 100 MeV.
    d_shell = s_full / without(shell=False) - 1.0
    assert -0.012 < d_shell[0] < -0.005
    assert abs(d_shell[1]) < 0.001
    # Barkas increases S for positive projectiles; ~0.3% at 10 MeV.
    d_barkas = s_full / without(barkas=False) - 1.0
    assert 0.001 < d_barkas[0] < 0.005
    assert 0.0 < d_barkas[2] < 0.001
    # Bloch reduces S; below 0.1% for protons at 10 MeV.
    d_bloch = s_full / without(bloch=False) - 1.0
    assert -0.001 < d_bloch[0] < 0.0
    # Density effect is zero for water below ~900 MeV.
    d_density = s_full / without(density_effect=False) - 1.0
    assert np.all(np.abs(d_density) < 1e-12)
    s_1gev = full.mass_stopping_power(1000.0)
    s_1gev_no_delta = AnalyticStoppingPower(
        materials.WATER,
        particles.PROTON,
        Corrections(density_effect=False),
        path="numpy",
    ).mass_stopping_power(1000.0)
    assert -0.003 < float((s_1gev / s_1gev_no_delta - 1.0)[0]) < 0.0


def test_icru90_i_value_lowers_stopping_power_by_about_half_a_percent() -> None:
    e = np.array([20.0, 100.0, 200.0])
    s75 = AnalyticStoppingPower(materials.WATER_ICRU49, particles.PROTON, path="numpy")
    s78 = AnalyticStoppingPower(materials.WATER_ICRU90, particles.PROTON, path="numpy")
    ratio = s78.mass_stopping_power(e) / s75.mass_stopping_power(e) - 1.0
    assert np.all(ratio < -0.004) and np.all(ratio > -0.008)


def test_bloch_series_matches_polynomial_expansion_for_small_y() -> None:
    mod = reference.python_reference("ionmc.physics.stopping")
    beta2 = 0.1
    y2 = 1.0 * mod.ALPHA_SQUARED / beta2
    series = mod.bloch_correction(1.0, beta2)
    poly = -y2 * (1.20206 - y2 * (1.042 - 0.855 * y2 + 0.343 * y2 * y2))
    assert math.isclose(series, poly, rel_tol=2e-3)


def test_barkas_function_interpolates_table_and_decays_beyond() -> None:
    from ionmc.physics.barkas_table import BARKAS_F, BARKAS_W

    mod = reference.python_reference("ionmc.physics.stopping")
    n = len(BARKAS_W)
    for w, f in zip(BARKAS_W, BARKAS_F, strict=True):
        value = mod.barkas_function(w, BARKAS_W, BARKAS_F, n)
        assert math.isclose(value, f, rel_tol=1e-12)
    mid = mod.barkas_function(0.15, BARKAS_W, BARKAS_F, n)
    assert min(BARKAS_F[8], BARKAS_F[9]) < mid < max(BARKAS_F[8], BARKAS_F[9])
    assert math.isclose(mod.barkas_function(20.0, BARKAS_W, BARKAS_F, n), 0.0025 * 0.5)
    assert mod.barkas_function(0.001, BARKAS_W, BARKAS_F, n) == BARKAS_F[0]


def test_max_energy_transfer_limits() -> None:
    mod = reference.python_reference("ionmc.physics.stopping")
    m_p = particles.PROTON.rest_energy_mev
    # Non-relativistic limit: T_max ~ 4 (m_e/M) T
    t = 1.0
    assert math.isclose(
        mod.max_energy_transfer(t, m_p),
        4.0 * mod.ELECTRON_MASS_MEV / m_p * t,
        rel_tol=2e-3,
    )
    # 1000 MeV: mass-ratio terms matter at the 0.35 percent level.
    t = 1000.0
    gamma = 1.0 + t / m_p
    bg2 = gamma * gamma - 1.0
    naive = 2.0 * mod.ELECTRON_MASS_MEV * bg2
    exact = mod.max_energy_transfer(t, m_p)
    assert 0.002 < 1.0 - exact / naive < 0.003


def test_parameters_for_water() -> None:
    p = build_parameters(materials.WATER, particles.PROTON)
    assert p.n_elem == 2 and p.n_table == 47
    assert math.isclose(sum(p.elem_f), 1.0, rel_tol=1e-12)
    assert p.elem_z == (1.0, 8.0) and p.elem_b == (1.8, 1.8)
    assert math.isclose(p.i_mev, 75.0e-6)
    assert p.use_density == 1.0
    p_no = build_parameters(materials.WATER_ICRU90, particles.PROTON)
    assert p_no.use_density == 0.0  # no density-effect parameters for 78 eV


def test_provenance_records_i_value_and_corrections(water_proton_python) -> None:
    prov = water_proton_python.provenance()
    assert prov["mean_excitation_energy_ev"] == 75.0
    assert "ICRU Report 49" in prov["mean_excitation_energy_source"]
    assert prov["corrections"] == {
        "density_effect": True,
        "shell": True,
        "barkas": True,
        "bloch": True,
    }
    assert prov["path"] == "python"


def test_unknown_path_rejected() -> None:
    with pytest.raises(ValueError):
        AnalyticStoppingPower(materials.WATER, particles.PROTON, path="fortran")


def test_kinematics_are_cancellation_free_in_float32() -> None:
    """Regression guard for the float32 error found in decision 0006.

    ``(beta gamma)^2`` written as ``gamma^2 - 1`` loses bits for low-energy
    protons in float32; the shared source uses ``tau (tau + 2)``. The numpy
    binding evaluates the shared source in whatever dtype it is given, so it
    serves as a float32 probe without Warp.
    """
    npy = reference.numpy_reference("ionmc.physics.stopping")
    py = reference.python_reference("ionmc.physics.stopping")
    mass = particles.PROTON.rest_energy_mev
    energies = np.linspace(2.0, 10.0, 201)
    t32, m32 = energies.astype(np.float32), np.float32(mass)
    exact = np.array([py.beta_gamma_squared(float(t), mass) for t in energies])
    shared = npy.beta_gamma_squared(t32, m32).astype(np.float64)
    naive = ((np.float32(1.0) + t32 / m32) ** 2 - np.float32(1.0)).astype(np.float64)
    assert shared.dtype == np.float64 and t32.dtype == np.float32
    assert np.max(np.abs(shared / exact - 1.0)) < 2.0e-7
    assert np.max(np.abs(naive / exact - 1.0)) > 1.0e-6, "naive form lost no bits?"
    beta2_exact = np.array([py.beta_squared(float(t), mass) for t in energies])
    beta2_shared = npy.beta_squared(t32, m32).astype(np.float64)
    assert np.max(np.abs(beta2_shared / beta2_exact - 1.0)) < 3.0e-7


def test_projectile_charge_scaling_of_corrections() -> None:
    """The ion generalisation: Barkas is odd in z, Bloch depends on z^2/beta^2."""
    py = reference.python_reference("ionmc.physics.stopping")
    p1 = build_parameters(materials.WATER, particles.PROTON)
    args = list(p1.as_call_args())
    del args[2]  # za_ratio is not an argument of stopping_number
    charge_index = 1
    t = 100.0
    contributions = {}
    for z in (1.0, 2.0, 6.0):
        a = list(args)
        a[charge_index] = z
        on = py.stopping_number(t, *a)
        a[-2] = 0.0  # use_barkas off
        off = py.stopping_number(t, *a)
        contributions[z] = on - off
    assert contributions[2.0] == pytest.approx(2.0 * contributions[1.0], rel=1e-12)
    assert contributions[6.0] == pytest.approx(6.0 * contributions[1.0], rel=1e-12)
    beta2 = py.beta_squared(t, particles.PROTON.rest_energy_mev)
    assert py.bloch_correction(2.0, beta2) == pytest.approx(
        py.bloch_correction(1.0, beta2 / 4.0), rel=1e-12
    )
    # Alpha at the same velocity as a 100 MeV proton: S scales as z^2 up to
    # the projectile-mass terms of T_max (about 4e-5 here).
    same_velocity = (
        t * particles.ALPHA.rest_energy_mev / particles.PROTON.rest_energy_mev
    )
    s_alpha = AnalyticStoppingPower(
        materials.WATER, particles.ALPHA, Corrections(barkas=False, bloch=False)
    ).mass_stopping_power(same_velocity)[0]
    s_p = AnalyticStoppingPower(
        materials.WATER, particles.PROTON, Corrections(barkas=False, bloch=False)
    ).mass_stopping_power(t)[0]
    # T_max grows with projectile mass (smaller m_e/M terms), so the ratio is
    # slightly above one.
    assert 1.0 < s_alpha / (4.0 * s_p) < 1.001
