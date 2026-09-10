"""Material definitions for stopping-power and transport calculations.

A :class:`Material` carries what the electronic stopping-power model needs:
mass density, elemental composition by mass fraction, a mean excitation
energy with explicit provenance, and (optionally) Sternheimer density-effect
parameters. The mean excitation energy is a *first-class, provenance-tagged
parameter* rather than a hidden constant because the two authoritative
recommendations for liquid water differ by about 4 percent (75 eV in
ICRU 37/49 and the NIST PSTAR tables; 78 +- 2 eV in ICRU 90), which changes
the stopping power by roughly one percent at therapeutic energies. Results
must record which value was used (decision ``0006``).

Only the elements needed for water and common tissue substitutes are listed
here; the table grows with the material library in Stage 3 of the roadmap.
Element atomic weights are standard atomic weights (IUPAC 2021 abridged);
element mean excitation energies are the ICRU 37 values as used by NIST
(ESTAR/PSTAR) and Geant4's NIST material builder. They enter only the shell
correction of compounds (see :mod:`ionmc.physics.stopping`); the compound
``I`` is never derived from them for liquid water.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Element:
    """A chemical element.

    ``atomic_weight`` is in g/mol (numerically the mass number in u);
    ``mean_excitation_energy_ev`` is the ICRU 37 elemental ``I`` in eV.
    """

    symbol: str
    atomic_number: int
    atomic_weight: float
    mean_excitation_energy_ev: float


#: Elements by symbol.
ELEMENTS: dict[str, Element] = {
    e.symbol: e
    for e in (
        Element("H", 1, 1.008, 19.2),
        Element("C", 6, 12.011, 81.0),
        Element("N", 7, 14.007, 82.0),
        Element("O", 8, 15.999, 95.0),
        Element("Na", 11, 22.990, 149.0),
        Element("Mg", 12, 24.305, 156.0),
        Element("P", 15, 30.974, 173.0),
        Element("S", 16, 32.06, 180.0),
        Element("Cl", 17, 35.45, 174.0),
        Element("Ar", 18, 39.95, 188.0),
        Element("K", 19, 39.098, 190.0),
        Element("Ca", 20, 40.078, 191.0),
    )
}


@dataclass(frozen=True)
class MeanExcitationEnergy:
    """Mean excitation energy ``I`` in eV together with its provenance."""

    value_ev: float
    source: str


@dataclass(frozen=True)
class DensityEffectParameters:
    """Sternheimer density-effect parameters.

    The parameterisation is that of Sternheimer, Berger and Seltzer, At. Data
    Nucl. Data Tables 30 (1984) 261, with ``x = log10(beta * gamma)``:

    * ``x >= x1``:       ``delta = 2 ln(10) x - cbar``
    * ``x0 <= x < x1``:  ``delta = 2 ln(10) x - cbar + a (x1 - x)^m``
    * ``x < x0``:        ``delta = delta0 * 10^(2 (x - x0))``

    ``delta0`` is zero for insulators. The parameter set is tied to the
    mean excitation energy it was derived with (``i_value_ev``).
    """

    cbar: float
    x0: float
    x1: float
    a: float
    m: float
    delta0: float
    i_value_ev: float
    source: str


@dataclass(frozen=True)
class Material:
    """A homogeneous material.

    Parameters
    ----------
    name:
        Human-readable identifier.
    density_g_per_cm3:
        Mass density in g/cm^3.
    mass_fractions:
        Mapping element symbol -> mass fraction; the fractions must sum to one
        within ``1e-6``.
    mean_excitation_energy:
        Mean excitation energy with provenance.
    density_effect:
        Sternheimer parameters, or ``None`` if unknown (the density effect is
        then omitted, which is accurate to better than 0.1 percent for
        condensed media below about 900 MeV/u but not for gases).
    source:
        Provenance of the composition and density.
    """

    name: str
    density_g_per_cm3: float
    mass_fractions: dict[str, float]
    mean_excitation_energy: MeanExcitationEnergy
    density_effect: DensityEffectParameters | None = None
    source: str = ""
    _elements: tuple[Element, ...] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        total = sum(self.mass_fractions.values())
        if abs(total - 1.0) > 1.0e-6:
            raise ValueError(
                f"mass fractions of {self.name!r} sum to {total!r}, expected 1"
            )
        unknown = [s for s in self.mass_fractions if s not in ELEMENTS]
        if unknown:
            raise ValueError(f"unknown element symbol(s) {unknown!r} in {self.name!r}")
        if self.density_g_per_cm3 <= 0.0:
            raise ValueError("density must be positive")
        if self.density_effect is not None and (
            abs(self.density_effect.i_value_ev - self.mean_excitation_energy.value_ev)
            > 1.0e-9
        ):
            raise ValueError(
                "density-effect parameters were derived for "
                f"I = {self.density_effect.i_value_ev} eV but the material uses "
                f"I = {self.mean_excitation_energy.value_ev} eV"
            )
        object.__setattr__(
            self,
            "_elements",
            tuple(ELEMENTS[s] for s in self.mass_fractions),
        )

    @property
    def elements(self) -> tuple[Element, ...]:
        """Elements in the order of ``mass_fractions``."""
        return self._elements

    def atoms_per_gram(self, symbol: str) -> float:
        """Number of atoms of ``symbol`` per gram, in units of mol/g."""
        return self.mass_fractions[symbol] / ELEMENTS[symbol].atomic_weight

    @property
    def electrons_per_gram_ratio(self) -> float:
        """Effective ``<Z/A>`` in mol/g: sum_i w_i Z_i / A_i.

        This is the quantity that multiplies the Bethe prefactor for a
        compound (Bragg additivity of the electron density).
        """
        return sum(
            self.atoms_per_gram(s) * ELEMENTS[s].atomic_number
            for s in self.mass_fractions
        )

    def electron_fraction(self, symbol: str) -> float:
        """Fraction of the material's electrons belonging to ``symbol``."""
        e = ELEMENTS[symbol]
        electrons = self.atoms_per_gram(symbol) * e.atomic_number
        return electrons / self.electrons_per_gram_ratio

    @property
    def electron_density_per_cm3(self) -> float:
        """Electron number density in 1/cm^3."""
        from ionmc.constants import AVOGADRO

        return AVOGADRO * self.electrons_per_gram_ratio * self.density_g_per_cm3

    def with_mean_excitation_energy(
        self, value_ev: float, source: str, name: str | None = None
    ) -> Material:
        """Return a copy with a different mean excitation energy.

        Density-effect parameters are dropped because they are tied to the
        ``I`` value they were fitted with.
        """
        return Material(
            name=name if name is not None else self.name,
            density_g_per_cm3=self.density_g_per_cm3,
            mass_fractions=dict(self.mass_fractions),
            mean_excitation_energy=MeanExcitationEnergy(value_ev, source),
            density_effect=None,
            source=self.source,
        )


_SBS_1984 = (
    "Sternheimer, Berger & Seltzer, At. Data Nucl. Data Tables 30 (1984) 261, "
    "as tabulated in Geant4 G4DensityEffectData"
)

#: Liquid water as defined by NIST PSTAR (material 276): density 1.0 g/cm^3,
#: mass fractions H 0.111894 / O 0.888106, I = 75.0 eV (ICRU 37/49), with the
#: Sternheimer-Berger-Seltzer density-effect parameters (Geant4 ``G4_WATER``
#: entry: Cbar 3.5017, x0 0.2400, x1 2.8004, a 0.09116, m 3.4773, delta0 0).
WATER_ICRU49: Material = Material(
    name="water_liquid",
    density_g_per_cm3=1.0,
    mass_fractions={"H": 0.111894, "O": 0.888106},
    mean_excitation_energy=MeanExcitationEnergy(
        75.0, "ICRU Report 49 (1993) / NIST PSTAR material 276"
    ),
    density_effect=DensityEffectParameters(
        cbar=3.5017,
        x0=0.2400,
        x1=2.8004,
        a=0.09116,
        m=3.4773,
        delta0=0.0,
        i_value_ev=75.0,
        source=_SBS_1984,
    ),
    source="NIST PSTAR material 276 composition and density",
)

#: Liquid water with the ICRU Report 90 (2016) mean excitation energy. No
#: density-effect parameters are attached (the SBS set was fitted for 75 eV);
#: the effect is below 0.1 percent for protons under 900 MeV.
WATER_ICRU90: Material = WATER_ICRU49.with_mean_excitation_energy(
    78.0, "ICRU Report 90 (2016), I = 78 +- 2 eV", name="water_liquid_icru90"
)

#: Default water definition. The ICRU 49 value is the default because the
#: reference tables used for validation (PSTAR-derived) are built on it;
#: decision 0006 records the policy.
WATER: Material = WATER_ICRU49
