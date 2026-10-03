"""Projectile definitions: charge, mass number and nuclear mass.

Nuclear masses [MeV/c2]. Proton, deuteron, triton and alpha are the CODATA 2018 values
(Tiesinga et al., Rev. Mod. Phys. 93 (2021) 025010; proton 938.27208816, electron
0.51099895 MeV/c2, atomic mass unit 931.49410242 MeV/c2). Helium-3, carbon-12 and oxygen-16
are computed from the AME2020 atomic masses (Wang et al., Chin. Phys. C 45 (2021) 030003;
3.01602932008, 12 (exact) and 15.99491461957 u) minus Z electron masses, neglecting the
electron binding energy.
"""

from __future__ import annotations

from dataclasses import dataclass

ELECTRON_MASS_MEV = 0.51099895000
AMU_MEV = 931.49410242


@dataclass(frozen=True)
class Projectile:
    """Ion: charge number ``z``, mass number ``a``, nuclear mass ``mass_mev`` [MeV/c2]."""

    name: str
    symbol: str
    z: int
    a: int
    mass_mev: float


def _nucleus(name: str, symbol: str, z: int, a: int, atomic_mass_u: float) -> Projectile:
    mass = atomic_mass_u * AMU_MEV - z * ELECTRON_MASS_MEV
    return Projectile(name, symbol, z, a, mass)


PROTON = Projectile("proton", "p", 1, 1, 938.27208816)
DEUTERON = Projectile("deuteron", "d", 1, 2, 1875.61294257)
TRITON = Projectile("triton", "t", 1, 3, 2808.92113298)
HELIUM3 = _nucleus("helium3", "He3", 2, 3, 3.01602932008)
ALPHA = Projectile("alpha", "He4", 2, 4, 3727.3794066)
CARBON12 = _nucleus("carbon12", "C12", 6, 12, 12.0)
OXYGEN16 = _nucleus("oxygen16", "O16", 8, 16, 15.99491461957)

PROJECTILES: dict[str, Projectile] = {
    p.name: p for p in (PROTON, DEUTERON, TRITON, HELIUM3, ALPHA, CARBON12, OXYGEN16)
}
