"""Tests for unit conventions, constants and material definitions."""

from __future__ import annotations

import math

import pytest

from ionmc import constants, materials, units


def test_bethe_prefactor_matches_pdg_value() -> None:
    # PDG: K = 4 pi N_A r_e^2 m_e c^2 = 0.307075 MeV mol^-1 cm^2
    assert math.isclose(constants.BETHE_K_MEV_CM2_PER_MOL, 0.307075, rel_tol=2e-6)


def test_unit_conversions_round_trip() -> None:
    assert units.ev_to_mev(units.mev_to_ev(1.234)) == pytest.approx(1.234)
    assert units.mm_to_cm(units.cm_to_mm(2.5)) == pytest.approx(2.5)
    # 1 MeV cm^2/g in water (1 g/cm^3) is 0.1 MeV/mm
    assert units.mass_stopping_power_to_linear(1.0, 1.0) == pytest.approx(0.1)
    # 7.718 g/cm^2 in water is 77.18 mm
    assert units.mass_range_to_linear(7.718, 1.0) == pytest.approx(77.18)
    assert units.kinetic_energy_per_nucleon(4800.0, 12) == pytest.approx(400.0)
    assert units.kinetic_energy_from_per_nucleon(400.0, 12) == pytest.approx(4800.0)


def test_water_definition_matches_pstar_material_276() -> None:
    w = materials.WATER
    assert w.density_g_per_cm3 == 1.0
    assert w.mass_fractions == {"H": 0.111894, "O": 0.888106}
    assert w.mean_excitation_energy.value_ev == 75.0
    # <Z/A> of water: PSTAR quotes 0.55509
    assert w.electrons_per_gram_ratio == pytest.approx(0.55509, rel=2e-4)
    assert materials.WATER_ICRU90.mean_excitation_energy.value_ev == 78.0
    assert materials.WATER_ICRU90.mass_fractions == w.mass_fractions


def test_material_validation() -> None:
    with pytest.raises(ValueError):
        materials.Material(
            "bad", 1.0, {"H": 0.5, "O": 0.4}, materials.MeanExcitationEnergy(75.0, "x")
        )
    with pytest.raises(ValueError):
        materials.Material(
            "bad", 1.0, {"Xx": 1.0}, materials.MeanExcitationEnergy(75.0, "x")
        )
    with pytest.raises(ValueError):
        materials.Material(
            "bad", -1.0, {"H": 1.0}, materials.MeanExcitationEnergy(75.0, "x")
        )


def test_electron_density_of_water() -> None:
    # 3.34e23 electrons per cm^3
    assert materials.WATER.electron_density_per_cm3 == pytest.approx(3.343e23, rel=1e-3)
