"""Analytic stopping layer and table construction versus ICRU 90 water tables.

These are ion-specific tabulated evaluation checks of the *analytic* Bethe
layer in its documented domain (above the blend window), using the shipped
ICRU 90 tables (protons, alphas, carbon ions in liquid water, I = 78 eV) as
the reference. Below the blend window the same tables are construction data
for water, so agreement there is by construction and is not asserted as
evidence.
"""

import numpy as np
import pytest

from ionmc.data import icru90
from ionmc.materials import Material, get_material
from ionmc.physics import stopping
from ionmc.physics.tables import (
    bethe_extrapolated_source,
    build_stopping_table,
    icru90_water_source,
    water_shape_source,
)
from ionmc.species import get_species

WATER = get_material("water")


@pytest.mark.parametrize(
    "species,t_min,tolerance",
    [("proton", 16.0, 0.0025), ("he4", 16.0, 0.0030), ("c12", 32.0, 0.0060)],
)
def test_bethe_layer_matches_icru90_above_blend_window(species, t_min, tolerance):
    table = icru90.water_table(species)
    t = table.energy_per_nucleon_mev
    mask = (t >= t_min) & (t <= 500.0)
    ref = table.electronic[mask]
    mine = stopping.mass_stopping_power(get_species(species), WATER, t[mask])
    rel = mine / ref - 1.0
    assert np.max(np.abs(rel)) < tolerance, (
        species,
        t[mask][np.argmax(np.abs(rel))],
        rel,
    )


def test_bethe_layer_is_inaccurate_below_its_domain_for_carbon():
    """Documented limitation: the analytic layer is >1 % low for carbon < 10 MeV/u."""
    table = icru90.water_table("c12")
    t = 5.0
    ref = table.electronic_at(t)
    mine = stopping.mass_stopping_power(get_species("c12"), WATER, np.array([t]))[0]
    assert mine / ref - 1.0 < -0.02


def test_shell_correction_lowers_low_energy_stopping():
    p = get_species("proton")
    t = np.array([10.0, 100.0])
    bare = stopping.mass_stopping_power(p, WATER, t, corrections=False)
    corrected = stopping.mass_stopping_power(p, WATER, t)
    assert corrected[0] < bare[0]
    assert abs(corrected[1] / bare[1] - 1.0) < 0.005


def test_effective_charge_limits():
    beta = np.array([1e-4, 0.01, 0.3, 0.9])
    z6 = stopping.effective_charge(6, beta)
    assert z6[0] < 0.1 and abs(z6[-1] - 6.0) < 1e-6
    assert np.all(np.diff(z6) > 0)


def test_i_value_dependence_matches_ICRU49_vs_ICRU90_offset():
    """Going from I = 78 eV to 75 eV raises stopping by ~0.4-0.5 % at 100 MeV."""
    w75 = Material("water75", 1.0, {"H": 0.111894, "O": 0.888106}, 75.0)
    p = get_species("proton")
    ratio = stopping.mass_stopping_power(p, w75, 100.0) / stopping.mass_stopping_power(
        p, WATER, 100.0
    )
    assert 1.003 < float(ratio) < 1.007


def test_water_table_is_continuous_monotonic_and_invertible():
    tbl = build_stopping_table("proton", WATER, icru90_water_source)
    s = tbl.electronic_mev_cm2_g
    r = tbl.csda_range_g_cm2
    assert np.all(np.isfinite(s)) and np.all(s > 0)
    assert np.all(np.diff(r) > 0)
    # smooth: log-log slope bounded and without kinks above 1 MeV/u (a blend
    # discontinuity would appear as a spike in the slope differences)
    t = tbl.t_mev_per_u
    slope = np.diff(np.log(s)) / np.diff(np.log(t))
    above = t[:-1] > 1.0
    assert np.abs(slope[above]).max() < 1.0
    assert np.abs(np.diff(slope[above])).max() < 0.05
    e = tbl.energy_at_range(tbl.range_at(150.0))
    assert e == pytest.approx(150.0, rel=1e-4)
    assert tbl.provenance["low_energy_source"].startswith("ICRU 90")


@pytest.mark.parametrize(
    "species,t,expected_rel",
    [
        ("proton", 100.0, 7.759),
        ("proton", 200.0, 26.09),
        ("he4", 100.0, None),
        ("c12", 100.0, 2.597),
        ("c12", 290.0, None),
    ],
)
def test_csda_ranges_match_icru90(species, t, expected_rel):
    """CSDA range from the constructed table versus the ICRU 90 table (g/cm²)."""
    tbl = build_stopping_table(species, WATER, icru90_water_source)
    ref = icru90.water_table(species).csda_range_at(t)
    mine = float(tbl.range_at(t))
    assert mine == pytest.approx(float(ref), rel=0.004), (species, t, mine, ref)
    if expected_rel is not None:
        assert float(ref) == pytest.approx(expected_rel, rel=2e-3)


def test_carbon_290_mev_per_u_range_in_water_is_about_163_mm():
    tbl = build_stopping_table("c12", WATER, icru90_water_source)
    assert 16.0 < float(tbl.range_at(290.0)) < 16.6  # g/cm² (water density 1) -> mm/10


def test_other_ions_use_scaled_proton_table_and_fallbacks():
    tbl = build_stopping_table("o16", WATER, icru90_water_source)
    assert "effective-charge" in tbl.provenance["low_energy_scaling"]
    lung = get_material("lung_inflated")
    s, prov = water_shape_source(get_species("proton"), lung, np.array([1.0, 5.0]))
    assert "water shape" in prov["low_energy_source"] and np.all(s > 0)
    s2, prov2 = bethe_extrapolated_source(
        get_species("proton"), lung, np.geomspace(1e-3, 100, 50)
    )
    assert prov2["low_energy_source"] == "bethe-extrapolated" and np.all(
        np.isfinite(s2)
    )


def test_stopping_power_scales_with_electron_density_and_i_value():
    p = get_species("proton")
    bone = get_material("bone_cortical")
    t = np.array([100.0])
    ratio = stopping.mass_stopping_power(p, bone, t) / stopping.mass_stopping_power(
        p, WATER, t
    )
    # bone: lower Z/A (0.515 vs 0.555) and higher I -> ~7-9 % lower mass stopping power
    assert 0.88 < float(ratio[0]) < 0.94
