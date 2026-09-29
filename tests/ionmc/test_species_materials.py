import math

import numpy as np
import pytest

from ionmc.materials import ELEMENTS, Material, get_material, list_materials
from ionmc.species import PROTON_MASS_MEV, SPECIES, beta_gamma, get_species


def test_species_masses_and_aliases():
    p = get_species("proton")
    assert p.z == 1 and p.a == 1 and p.mass_mev == PROTON_MASS_MEV
    assert get_species("alpha") is get_species("he4")
    c = get_species("carbon")
    assert c.z == 6 and c.a == 12
    assert abs(c.mass_mev - 11174.86) < 0.5  # 12 u minus 6 electron masses
    assert abs(get_species("o16").mass_mev - 14895.08) < 0.5
    with pytest.raises(ValueError):
        get_species("unobtainium")
    for s in SPECIES.values():
        assert s.mass_mev / s.a == pytest.approx(931.5, rel=0.01)


def test_beta_gamma_per_nucleon():
    beta, gamma = beta_gamma(get_species("proton"), 150.0)
    assert gamma == pytest.approx(1.0 + 150.0 / PROTON_MASS_MEV)
    assert beta == pytest.approx(0.5066, abs=2e-4)
    # equal velocity for equal energy per nucleon (up to the mass-per-nucleon spread)
    beta_c, _ = beta_gamma(get_species("c12"), 150.0)
    assert beta_c == pytest.approx(beta, rel=5e-3)  # mass/nucleon of C-12 < m_p


def test_water_material_and_bragg_additivity():
    water = get_material("water")
    assert water.mean_excitation_ev == 78.0
    assert water.electrons_per_gram == pytest.approx(0.5551, rel=1e-3)
    # Bragg additivity of elemental I values gives ~75 eV for H2O (ICRU 37 style)
    assert 65.0 < water.bragg_additivity_i_ev() < 80.0
    assert get_material("G4_WATER") is water


def test_all_materials_are_consistent():
    for name in list_materials():
        m = get_material(name)
        assert abs(sum(m.composition.values()) - 1.0) < 1e-3
        assert m.density_g_cm3 > 0 and m.mean_excitation_ev > 0
        assert set(m.composition) <= set(ELEMENTS)


def test_invalid_material_rejected():
    with pytest.raises(ValueError):
        Material("bad", 1.0, {"H": 0.5, "O": 0.4}, 78.0)
    with pytest.raises(ValueError):
        Material("bad", 1.0, {"Xx": 1.0}, 78.0)
    with pytest.raises(ValueError):
        get_material("kryptonite")


def test_density_override_keeps_composition():
    lung = get_material("lung_inflated").with_density(0.5)
    assert lung.density_g_cm3 == 0.5
    assert lung.composition == get_material("lung_inflated").composition
    assert math.isclose(lung.mean_excitation_ev, 75.3)
    assert isinstance(np.asarray(lung.electrons_per_gram), np.ndarray)
