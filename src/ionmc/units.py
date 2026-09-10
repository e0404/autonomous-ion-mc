"""Internal unit system and conventions of ionmc.

The conventions below are fixed by decision ``0004`` (see ``decisions/``) and
apply to every public and internal API unless a function's docstring says
otherwise explicitly. Quantities never carry implicit units: names of
variables and parameters that are not in the internal unit carry a suffix
(``_cm``, ``_ev``, ``_g_per_cm2``).

Internal units
--------------

=====================  ==================  =============================
Quantity               Unit                Notes
=====================  ==================  =============================
kinetic energy         MeV                 total kinetic energy of the particle,
                                           not per nucleon (helpers convert)
rest energy / mass     MeV                 m c^2
length, position       mm                  treatment-planning convention
time                   not used
mass density           g/cm^3              as in CT calibration and ICRU tables
mass stopping power    MeV cm^2/g          as in ICRU 49 / PSTAR tables
linear stopping power  MeV/mm              mass stopping power * density / 10
mass range (CSDA)      g/cm^2              as in ICRU 49 / PSTAR tables
linear range           mm                  mass range / density * 10
mean excitation energy eV                  as tabulated by ICRU
angle                  rad
charge                 e (elementary)
=====================  ==================  =============================

Coordinate system
-----------------

Right-handed Cartesian coordinates in mm. Geometry objects define their own
origin; the beam frame and any patient/DICOM frame are related by explicit
rigid transforms documented where they are introduced (Stage 3 of the
roadmap). No axis is assumed to be "the beam axis" by the physics layer.

Helper functions convert between the internal units and the units that
external data sources and users commonly employ.
"""

from __future__ import annotations

from ionmc.constants import ATOMIC_MASS_UNIT_MEV, EV_PER_MEV, MM_PER_CM


def mev_to_ev(energy_mev: float) -> float:
    """Convert an energy from MeV to eV."""
    return energy_mev * EV_PER_MEV


def ev_to_mev(energy_ev: float) -> float:
    """Convert an energy from eV to MeV."""
    return energy_ev / EV_PER_MEV


def cm_to_mm(length_cm: float) -> float:
    """Convert a length from cm to mm."""
    return length_cm * MM_PER_CM


def mm_to_cm(length_mm: float) -> float:
    """Convert a length from mm to cm."""
    return length_mm / MM_PER_CM


def mass_stopping_power_to_linear(
    mass_stopping_power_mev_cm2_per_g: float, density_g_per_cm3: float
) -> float:
    """Convert a mass stopping power [MeV cm^2/g] to a linear one [MeV/mm]."""
    return mass_stopping_power_mev_cm2_per_g * density_g_per_cm3 / MM_PER_CM


def mass_range_to_linear(
    mass_range_g_per_cm2: float, density_g_per_cm3: float
) -> float:
    """Convert a mass range [g/cm^2] to a linear range [mm]."""
    return mass_range_g_per_cm2 / density_g_per_cm3 * MM_PER_CM


def kinetic_energy_per_nucleon(kinetic_energy_mev: float, mass_number: int) -> float:
    """Kinetic energy per nucleon [MeV/u] from total kinetic energy [MeV]."""
    return kinetic_energy_mev / mass_number


def kinetic_energy_from_per_nucleon(
    kinetic_energy_mev_per_u: float, mass_number: int
) -> float:
    """Total kinetic energy [MeV] from kinetic energy per nucleon [MeV/u]."""
    return kinetic_energy_mev_per_u * mass_number


def rest_energy_from_mass_number(mass_number: int) -> float:
    """Approximate rest energy [MeV] of a nucleus with the given mass number.

    Uses ``A * m_u c^2``; exact isotopic masses are supplied by the particle
    definitions where the difference matters.
    """
    return mass_number * ATOMIC_MASS_UNIT_MEV
