"""Transported particle species.

Kinetic energies of ions are expressed per nucleon (MeV/u) throughout the
physics layer because electronic stopping depends on velocity; total kinetic
energy in MeV is used at API boundaries where explicitly stated. Masses are
atomic-mass-unit based nuclear masses (electrons excluded) in MeV/c², from the
2020 atomic mass evaluation (AME2020) values rounded to the keV level.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

AMU_MEV: Final[float] = 931.49410242  # MeV per atomic mass unit (CODATA 2018)
ELECTRON_MASS_MEV: Final[float] = 0.51099895
PROTON_MASS_MEV: Final[float] = 938.27208816


@dataclass(frozen=True)
class Species:
    """A charged nuclear species transported by the code."""

    name: str
    z: int  # nuclear charge number
    a: int  # mass number
    mass_mev: float  # nuclear rest mass, MeV/c^2

    @property
    def mass_per_nucleon_mev(self) -> float:
        return self.mass_mev / self.a

    def total_kinetic_mev(self, energy_per_nucleon: float) -> float:
        return energy_per_nucleon * self.a

    def energy_per_nucleon_mev(self, total_kinetic_mev: float) -> float:
        return total_kinetic_mev / self.a


def _nuclear_mass(atomic_mass_u: float, z: int) -> float:
    # Atomic mass minus electron masses (binding of electrons neglected, < 1 keV)
    return atomic_mass_u * AMU_MEV - z * ELECTRON_MASS_MEV


# Atomic masses in u (AME2020); the proton uses the CODATA value directly.
_TABLE: dict[str, tuple[int, int, float]] = {
    "proton": (1, 1, 1.007276466621),  # nuclear mass in u
    "deuteron": (1, 2, 2.014101778),
    "triton": (1, 3, 3.016049281),
    "he3": (2, 3, 3.016029322),
    "he4": (2, 4, 4.002603254),
    "li6": (3, 6, 6.015122887),
    "li7": (3, 7, 7.016003434),
    "be7": (4, 7, 7.016928717),
    "be9": (4, 9, 9.012183065),
    "be10": (4, 10, 10.013534695),
    "b10": (5, 10, 10.012936862),
    "b11": (5, 11, 11.009305167),
    "c10": (6, 10, 10.016853218),
    "c11": (6, 11, 11.011432597),
    "c12": (6, 12, 12.0),
    "n13": (7, 13, 13.005738609),
    "n14": (7, 14, 14.003074004),
    "n15": (7, 15, 15.000108899),
    "o14": (8, 14, 14.008596359),
    "o15": (8, 15, 15.003065618),
    "o16": (8, 16, 15.994914619),
}

SPECIES: dict[str, Species] = {}
for _name, (_z, _a, _mass_u) in _TABLE.items():
    if _name == "proton":
        _mass = PROTON_MASS_MEV
    else:
        _mass = _nuclear_mass(_mass_u, _z)
    SPECIES[_name] = Species(_name, _z, _a, _mass)

ALIASES: dict[str, str] = {
    "p": "proton",
    "h1": "proton",
    "d": "deuteron",
    "h2": "deuteron",
    "t": "triton",
    "h3": "triton",
    "alpha": "he4",
    "helium": "he4",
    "helium-4": "he4",
    "carbon": "c12",
    "carbon-12": "c12",
    "oxygen": "o16",
    "oxygen-16": "o16",
}

PRIMARY_SPECIES: tuple[str, ...] = ("proton", "he4", "c12", "o16")


def get_species(name: str | Species) -> Species:
    """Resolve a species by canonical name or alias (case-insensitive)."""
    if isinstance(name, Species):
        return name
    key = name.strip().lower()
    key = ALIASES.get(key, key)
    try:
        return SPECIES[key]
    except KeyError:
        raise ValueError(
            f"Unknown species {name!r}; known: {sorted(SPECIES)} "
            f"(aliases: {sorted(ALIASES)})"
        ) from None


def beta_gamma(species: Species, energy_per_nucleon: float) -> tuple[float, float]:
    """Relativistic beta and gamma for the given kinetic energy per nucleon."""
    gamma = 1.0 + energy_per_nucleon * species.a / species.mass_mev
    beta = (1.0 - 1.0 / (gamma * gamma)) ** 0.5
    return beta, gamma
