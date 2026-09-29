"""Materials: elemental composition, density and mean excitation energy.

Compositions (mass fractions) follow ICRU Report 44/46 tissue substitutes as
tabulated in the Geant4 NIST material database (public, derived from ICRU),
densities in g/cm³. Mean excitation energies: elements from ICRU 37/49; water,
air and graphite from ICRU Report 90 (2016); compounds without a recommended
value use Bragg additivity of the elemental values. The values are explicit
data, not fitted to any Monte Carlo code.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

# Element data: Z, standard atomic weight (g/mol), mean excitation energy I (eV)
# I values: ICRU Report 37 (1984) / ICRU 49 elemental table.
ELEMENTS: dict[str, tuple[int, float, float]] = {
    "H": (1, 1.008, 19.2),
    "C": (6, 12.011, 78.0),  # ICRU 90 graphite value; amorphous carbon 81 (ICRU 37)
    "N": (7, 14.007, 82.0),
    "O": (8, 15.999, 95.0),
    "F": (9, 18.998, 115.0),
    "Na": (11, 22.990, 149.0),
    "Mg": (12, 24.305, 156.0),
    "Al": (13, 26.982, 166.0),
    "Si": (14, 28.085, 173.0),
    "P": (15, 30.974, 173.0),
    "S": (16, 32.06, 180.0),
    "Cl": (17, 35.45, 174.0),
    "Ar": (18, 39.948, 188.0),
    "K": (19, 39.098, 190.0),
    "Ca": (20, 40.078, 191.0),
    "Ti": (22, 47.867, 233.0),
    "Fe": (26, 55.845, 286.0),
    "Cu": (29, 63.546, 322.0),
    "Zn": (30, 65.38, 330.0),
    "Pb": (82, 207.2, 823.0),
}


@dataclass(frozen=True)
class Material:
    """A homogeneous material with mass-fraction composition."""

    name: str
    density_g_cm3: float
    composition: dict[str, float]  # element symbol -> mass fraction (sums to 1)
    mean_excitation_ev: float
    source: str = ""
    aliases: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        total = sum(self.composition.values())
        if abs(total - 1.0) > 1e-3:
            raise ValueError(f"{self.name}: mass fractions sum to {total:.4f}, not 1")
        unknown = set(self.composition) - set(ELEMENTS)
        if unknown:
            raise ValueError(f"{self.name}: unknown elements {sorted(unknown)}")
        if self.density_g_cm3 <= 0 or self.mean_excitation_ev <= 0:
            raise ValueError(f"{self.name}: density and I must be positive")

    @property
    def electrons_per_gram(self) -> float:
        """Z/A weighted electron density in mol electrons per gram."""
        return sum(
            w * ELEMENTS[e][0] / ELEMENTS[e][1] for e, w in self.composition.items()
        )

    @property
    def z_over_a(self) -> float:
        return self.electrons_per_gram

    def bragg_additivity_i_ev(self) -> float:
        """Mean excitation energy from Bragg additivity of the elemental values."""
        num = sum(
            w * ELEMENTS[e][0] / ELEMENTS[e][1] * math.log(ELEMENTS[e][2])
            for e, w in self.composition.items()
        )
        return math.exp(num / self.electrons_per_gram)

    def with_density(self, density_g_cm3: float, name: str | None = None) -> Material:
        return Material(
            name or f"{self.name}@{density_g_cm3:g}",
            density_g_cm3,
            dict(self.composition),
            self.mean_excitation_ev,
            self.source + " (density overridden)",
        )


def _m(name, rho, comp, i_ev=None, source="", aliases=()):
    mat = Material(name, rho, comp, i_ev or 1.0, source, tuple(aliases))
    if i_ev is None:
        mat = Material(
            name,
            rho,
            comp,
            mat.bragg_additivity_i_ev(),
            source + "; I by Bragg additivity",
            tuple(aliases),
        )
    return mat


# Compositions: Geant4 NIST database (ICRU 44/46 tissue substitutes), mass fractions.
MATERIALS: dict[str, Material] = {}
for _mat in [
    _m(
        "water",
        1.0,
        {"H": 0.111894, "O": 0.888106},
        78.0,
        "ICRU 90 (I = 78 eV); composition H2O",
        ("g4_water", "liquid_water"),
    ),
    _m(
        "air",
        0.00120479,
        {"C": 0.000124, "N": 0.755267, "O": 0.231781, "Ar": 0.012827},
        85.7,
        "ICRU 90 (I = 85.7 eV); dry air, sea level",
        ("g4_air",),
    ),
    _m("graphite", 2.21, {"C": 1.0}, 81.0, "ICRU 90 (I = 81 eV)"),
    _m(
        "pmma",
        1.19,
        {"H": 0.080538, "C": 0.599848, "O": 0.319614},
        74.0,
        "ICRU 37 (I = 74 eV); C5H8O2",
        ("lucite", "plexiglass"),
    ),
    _m(
        "polyethylene",
        0.94,
        {"H": 0.143711, "C": 0.856289},
        57.4,
        "ICRU 37 (I = 57.4 eV)",
    ),
    _m(
        "adipose_tissue",
        0.95,
        {
            "H": 0.114,
            "C": 0.598,
            "N": 0.007,
            "O": 0.278,
            "Na": 0.001,
            "S": 0.001,
            "Cl": 0.001,
        },
        63.2,
        "ICRP adipose tissue (Geant4 G4_ADIPOSE_TISSUE_ICRP), I = 63.2 eV",
    ),
    _m(
        "muscle_skeletal",
        1.05,
        {
            "H": 0.102,
            "C": 0.143,
            "N": 0.034,
            "O": 0.710,
            "Na": 0.001,
            "P": 0.002,
            "S": 0.003,
            "Cl": 0.001,
            "K": 0.004,
        },
        75.3,
        "ICRP skeletal muscle (G4_MUSCLE_SKELETAL_ICRP), I = 75.3 eV",
    ),
    _m(
        "soft_tissue",
        1.03,
        {
            "H": 0.105,
            "C": 0.256,
            "N": 0.027,
            "O": 0.602,
            "Na": 0.001,
            "P": 0.002,
            "S": 0.003,
            "Cl": 0.002,
            "K": 0.003,
        },
        72.3,
        "ICRP soft tissue (G4_TISSUE_SOFT_ICRP), I = 72.3 eV",
    ),
    _m(
        "lung_inflated",
        0.26,
        {
            "H": 0.103,
            "C": 0.105,
            "N": 0.031,
            "O": 0.749,
            "Na": 0.002,
            "P": 0.002,
            "S": 0.003,
            "Cl": 0.003,
            "K": 0.002,
        },
        75.3,
        "ICRP lung, inflated (G4_LUNG_ICRP composition; density 0.26), I = 75.3 eV",
    ),
    _m(
        "bone_compact",
        1.85,
        {
            "H": 0.063984,
            "C": 0.278,
            "N": 0.027,
            "O": 0.410016,
            "Mg": 0.002,
            "P": 0.07,
            "S": 0.002,
            "Ca": 0.147,
        },
        91.9,
        "ICRU compact bone (G4_BONE_COMPACT_ICRU), I = 91.9 eV",
    ),
    _m(
        "bone_cortical",
        1.92,
        {
            "H": 0.034,
            "C": 0.155,
            "N": 0.042,
            "O": 0.435,
            "Na": 0.001,
            "Mg": 0.002,
            "P": 0.103,
            "S": 0.003,
            "Ca": 0.225,
        },
        110.0,
        "ICRP cortical bone (G4_BONE_CORTICAL_ICRP), I = 110 eV",
    ),
    _m("aluminium", 2.699, {"Al": 1.0}, 166.0, "ICRU 37", ("aluminum", "al")),
    _m("lead", 11.35, {"Pb": 1.0}, 823.0, "ICRU 37", ("pb",)),
    _m("titanium", 4.54, {"Ti": 1.0}, 233.0, "ICRU 37", ("ti",)),
]:
    MATERIALS[_mat.name] = _mat
    for _alias in _mat.aliases:
        MATERIALS[_alias] = _mat


def get_material(name: str | Material) -> Material:
    if isinstance(name, Material):
        return name
    key = name.strip().lower()
    try:
        return MATERIALS[key]
    except KeyError:
        raise ValueError(
            f"Unknown material {name!r}; known: {list_materials()}"
        ) from None


def list_materials() -> list[str]:
    return sorted({mat.name for mat in MATERIALS.values()})
