"""Tests for elements, materials and projectiles."""

import math

import pytest

from ionmc.materials import (
    AIR,
    MATERIALS,
    WATER,
    Material,
    fractions_from_atom_counts,
    water,
)
from ionmc.physics.projectiles import ALPHA, CARBON12, PROJECTILES, PROTON


def test_water_z_over_a() -> None:
    assert WATER.z_over_a == pytest.approx(0.555087, abs=1e-5)
    assert sum(WATER.mass_fractions.values()) == pytest.approx(1.0)
    assert WATER.mean_excitation_eV == 78.0
    assert water(75.0).mean_excitation_eV == 75.0
    assert WATER.sternheimer is not None and WATER.sternheimer.x0 == 0.2400


def test_electron_density_of_water() -> None:
    # 3.343e23 electrons per gram of water at Z/A = 0.55509 mol/g.
    assert WATER.electron_density_cm3 == pytest.approx(3.343e23, rel=1e-3)


@pytest.mark.parametrize("material", list(MATERIALS.values()), ids=list(MATERIALS))
def test_predefined_materials_are_consistent(material: Material) -> None:
    assert sum(material.mass_fractions.values()) == pytest.approx(1.0, abs=1e-5)
    assert 0.4 < material.z_over_a < 0.6
    # Bragg additivity is within 20 % of the tabulated I of the same material.
    assert material.mean_excitation_eV / math.exp(material.ln_I_bragg) < 1.25
    assert material.ln_I == pytest.approx(math.log(material.mean_excitation_eV))


def test_bragg_additivity_when_I_missing() -> None:
    m = Material("w", 1.0, fractions_from_atom_counts({"H": 2, "O": 1}))
    assert m.mean_excitation_eV == pytest.approx(math.exp(m.ln_I_bragg))
    assert 65.0 < m.mean_excitation_eV < 72.0


def test_bad_fractions_rejected() -> None:
    with pytest.raises(ValueError):
        Material("bad", 1.0, {"H": 0.5, "O": 0.4})


def test_with_I_keeps_composition() -> None:
    m = AIR.with_I(90.0)
    assert m.mean_excitation_eV == 90.0 and m.density_g_cm3 == AIR.density_g_cm3


def test_projectiles() -> None:
    assert PROTON.mass_mev == 938.27208816
    assert set(PROJECTILES) == {
        "proton", "deuteron", "triton", "helium3", "alpha", "carbon12", "oxygen16"
    }  # fmt: skip
    assert ALPHA.mass_mev == pytest.approx(3727.379, abs=1e-3)
    # carbon-12 nucleus: 12 u minus six electron masses
    assert CARBON12.mass_mev == pytest.approx(12 * 931.49410242 - 6 * 0.51099895, abs=1e-6)
    assert CARBON12.z == 6 and CARBON12.a == 12
