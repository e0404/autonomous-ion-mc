"""Elements and materials for stopping-power calculations.

Mean excitation energies of the elements are the values tabulated by NIST for the
ESTAR/PSTAR/ASTAR programs (ICRU Report 37 (1984) and ICRU Report 49 (1993)): H 19.2,
C 81.0 (graphite-like carbon as used in the NIST element table), N 82.0, O 95.0,
Na 149, Mg 156, P 173, S 180, Cl 174, Ar 188, K 190, Ca 191 eV. Atomic weights are
standard (IUPAC) values as used in the Geant4 element builder, ``G4NistElementBuilder``.

Compositions, densities and mean excitation energies of the predefined materials are
taken from the Geant4 NIST material builder, ``G4NistMaterialBuilder.cc`` (Geant4 v11.4.2,
https://github.com/Geant4/geant4/blob/v11.4.2/source/materials/src/G4NistMaterialBuilder.cc),
which for the ICRP tissues follows ICRU Report 46 (1992) / ICRP Publication 23 (1975)
as stated in that file's change log, and ICRU Report 44 (1989) for ICRU tissues.
Sternheimer parameters: R.M. Sternheimer, M.J. Berger, S.M. Seltzer, At. Data Nucl.
Data Tables 30 (1984) 261.

Units: densities in g/cm3, mean excitation energies in eV, A in g/mol.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np

N_A = 6.02214076e23  # 1/mol (exact, SI 2019)


@dataclass(frozen=True)
class Element:
    """A chemical element: symbol, atomic number, molar mass [g/mol], I [eV]."""

    symbol: str
    Z: int
    A_g_mol: float
    I_eV: float


ELEMENTS: dict[str, Element] = {
    e.symbol: e
    for e in (
        Element("H", 1, 1.00794, 19.2),
        Element("C", 6, 12.0107, 81.0),
        Element("N", 7, 14.0067, 82.0),
        Element("O", 8, 15.9994, 95.0),
        Element("Na", 11, 22.98977, 149.0),
        Element("Mg", 12, 24.305, 156.0),
        Element("P", 15, 30.97376, 173.0),
        Element("S", 16, 32.065, 180.0),
        Element("Cl", 17, 35.453, 174.0),
        Element("Ar", 18, 39.948, 188.0),
        Element("K", 19, 39.0983, 190.0),
        Element("Ca", 20, 40.078, 191.0),
    )
}
_BY_Z = {e.Z: e for e in ELEMENTS.values()}


@dataclass(frozen=True)
class SternheimerParameters:
    """Sternheimer density-effect parameters (dimensionless): x0, x1, Cbar, a, m."""

    x0: float
    x1: float
    cbar: float
    a: float
    m: float

    def __post_init__(self) -> None:
        for name in ("x0", "x1", "cbar", "a", "m"):
            if not math.isfinite(getattr(self, name)):
                raise ValueError(f"Sternheimer parameter {name} is not finite")


@dataclass(frozen=True)
class Material:
    """A homogeneous material.

    ``mass_fractions`` maps element symbol to mass fraction (normalised to sum 1 within
    1e-6). ``I_eV`` is the mean excitation energy [eV]; if None it follows from Bragg
    additivity (:attr:`ln_I_bragg`). ``sternheimer`` may be None (density effect = 0).
    """

    name: str
    density_g_cm3: float
    mass_fractions: Mapping[str, float]
    I_eV: float | None = None
    sternheimer: SternheimerParameters | None = None
    source: str = ""
    _fractions: tuple[tuple[Element, float], ...] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        fractions = dict(self.mass_fractions)
        if not fractions:
            raise ValueError(f"{self.name}: mass fractions are empty")
        for symbol, w in fractions.items():
            if symbol not in ELEMENTS:
                raise ValueError(f"{self.name}: unknown element symbol {symbol!r}")
            if not (math.isfinite(w) and w > 0.0):
                raise ValueError(f"{self.name}: mass fraction of {symbol} must be finite and > 0")
        total = math.fsum(fractions.values())
        if abs(total - 1.0) > 1e-6:
            raise ValueError(f"{self.name}: mass fractions sum to {total}, not 1 (tolerance 1e-6)")
        if not (math.isfinite(self.density_g_cm3) and self.density_g_cm3 > 0.0):
            raise ValueError(f"{self.name}: density must be finite and positive")
        if self.I_eV is not None and not (math.isfinite(self.I_eV) and self.I_eV > 0.0):
            raise ValueError(f"{self.name}: I_eV must be finite and positive")
        fr = tuple((ELEMENTS[s], w / total) for s, w in fractions.items())
        object.__setattr__(self, "mass_fractions", MappingProxyType(fractions))
        object.__setattr__(self, "_fractions", fr)

    @property
    def z_over_a(self) -> float:
        """Mean Z/A [mol/g] = sum_i w_i Z_i / A_i."""
        return sum(w * e.Z / e.A_g_mol for e, w in self._fractions)

    @property
    def electron_density_cm3(self) -> float:
        """Electron density [1/cm3] = N_A * (Z/A) * density."""
        return N_A * self.z_over_a * self.density_g_cm3

    @property
    def ln_I_bragg(self) -> float:
        """Bragg-additivity ln I (I in eV): sum_i w_i (Z/A)_i ln I_i / (Z/A)."""
        s = sum(w * e.Z / e.A_g_mol * math.log(e.I_eV) for e, w in self._fractions)
        return s / self.z_over_a

    @property
    def mean_excitation_eV(self) -> float:
        """Mean excitation energy [eV]: ``I_eV`` if set, else the Bragg-additivity value."""
        return self.I_eV if self.I_eV is not None else math.exp(self.ln_I_bragg)

    @property
    def ln_I(self) -> float:
        """Natural log of the mean excitation energy in eV."""
        return math.log(self.mean_excitation_eV)

    def with_I(self, I_eV: float, name: str | None = None) -> Material:
        """Return a copy with a different mean excitation energy [eV]."""
        return Material(
            name or f"{self.name} (I={I_eV:g} eV)",
            self.density_g_cm3,
            dict(self.mass_fractions),
            I_eV,
            self.sternheimer,
            self.source,
        )


def fractions_from_atom_counts(counts: dict[str, int]) -> dict[str, float]:
    """Convert a chemical formula (symbol -> atoms) to mass fractions."""
    mass = {s: n * ELEMENTS[s].A_g_mol for s, n in counts.items()}
    total = float(np.sum(list(mass.values())))
    return {s: m / total for s, m in mass.items()}


_G4 = "G4NistMaterialBuilder v11.4.2: "
_STERNHEIMER_WATER = SternheimerParameters(x0=0.2400, x1=2.8004, cbar=3.5017, a=0.09116, m=3.4773)


def water(I_eV: float = 78.0) -> Material:
    """Liquid water (1.0 g/cm3, H2O) with mean excitation energy ``I_eV`` [eV].

    The default of 78 eV is the ICRU Report 90 (2016) recommendation (also used by
    G4_WATER); NIST PSTAR/ASTAR use 75 eV (ICRU Report 49).
    """
    return Material(
        "water",
        1.0,
        fractions_from_atom_counts({"H": 2, "O": 1}),
        I_eV,
        _STERNHEIMER_WATER,
        _G4 + "G4_WATER; I and Sternheimer parameters as cited in module docstring",
    )


WATER = water()
AIR = Material(
    "air",
    0.00120479,
    # Geant4 weights sum to 0.999999; normalised here to sum to 1.
    {
        k: v / 0.999999
        for k, v in {"C": 0.000124, "N": 0.755267, "O": 0.231781, "Ar": 0.012827}.items()
    },
    85.7,
    None,
    _G4 + "G4_AIR (dry air, near sea level)",
)
PMMA = Material(
    "pmma",
    1.19,
    fractions_from_atom_counts({"C": 5, "H": 8, "O": 2}),
    74.0,
    None,
    _G4 + "G4_PLEXIGLASS (C5H8O2)",
)
ADIPOSE_TISSUE_ICRP = Material(
    "adipose_tissue_icrp",
    0.95,
    {"H": 0.114, "C": 0.598, "N": 0.007, "O": 0.278, "Na": 0.001, "S": 0.001, "Cl": 0.001},
    63.2,
    None,
    _G4 + "G4_ADIPOSE_TISSUE_ICRP (ICRU 46 / ICRP 23)",
)
MUSCLE_SKELETAL_ICRP = Material(
    "muscle_skeletal_icrp",
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
    None,
    _G4 + "G4_MUSCLE_SKELETAL_ICRP (ICRU 46 / ICRP 23)",
)
BONE_COMPACT_ICRU = Material(
    "bone_compact_icru",
    1.85,
    {
        "H": 0.064,
        "C": 0.278,
        "N": 0.027,
        "O": 0.410,
        "Mg": 0.002,
        "P": 0.07,
        "S": 0.002,
        "Ca": 0.147,
    },
    91.9,
    None,
    _G4 + "G4_BONE_COMPACT_ICRU (ICRU 44)",
)
LUNG_ICRP = Material(
    "lung_icrp",
    1.04,
    {
        "H": 0.105,
        "C": 0.083,
        "N": 0.023,
        "O": 0.779,
        "Na": 0.002,
        "P": 0.001,
        "S": 0.002,
        "Cl": 0.003,
        "K": 0.002,
    },
    75.3,
    None,
    _G4 + "G4_LUNG_ICRP (ICRU 46 / ICRP 23)",
)

MATERIALS: dict[str, Material] = {
    m.name: m
    for m in (
        WATER,
        AIR,
        PMMA,
        ADIPOSE_TISSUE_ICRP,
        MUSCLE_SKELETAL_ICRP,
        BONE_COMPACT_ICRU,
        LUNG_ICRP,
    )
}
