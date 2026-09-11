"""Projectile definitions.

The architecture must not assume proton-only transport; every projectile is a
:class:`Particle` with charge, rest energy and mass number. Rest energies are
CODATA 2018 / AME 2020 values in MeV.
"""

from __future__ import annotations

from dataclasses import dataclass

from ionmc.constants import PROTON_MASS_MEV


@dataclass(frozen=True)
class Particle:
    """A charged projectile.

    ``charge`` is in units of the elementary charge (the *ion* charge used in
    the Bethe formula, i.e. the nuclear charge for fully stripped ions);
    ``rest_energy_mev`` is ``m c^2``; ``mass_number`` is ``A``.
    """

    name: str
    charge: float
    rest_energy_mev: float
    mass_number: int


PROTON: Particle = Particle("proton", 1.0, PROTON_MASS_MEV, 1)
#: Alpha particle (helium-4 nucleus), rest energy 3727.379 MeV.
ALPHA: Particle = Particle("alpha", 2.0, 3727.3794066, 4)
#: Carbon-12 nucleus, rest energy 11174.86 MeV (12 u minus 6 electrons' mass,
#: binding energies neglected at the 1e-6 level).
CARBON_12: Particle = Particle("carbon-12", 6.0, 11174.8628, 12)
#: Oxygen-16 nucleus, rest energy 14895.08 MeV.
OXYGEN_16: Particle = Particle("oxygen-16", 8.0, 14895.0796, 16)
#: Boron-11 nucleus, rest energy ~10252.5 MeV (atomic mass 11.009305 u minus 5
#: electrons). Used as the representative "heavy" carbon fragment (Z=3-5 lumped)
#: in the bounded fragmentation model (decision 0029).
BORON_11: Particle = Particle("boron-11", 5.0, 10252.5476, 11)

PARTICLES: dict[str, Particle] = {
    p.name: p for p in (PROTON, ALPHA, CARBON_12, OXYGEN_16, BORON_11)
}
