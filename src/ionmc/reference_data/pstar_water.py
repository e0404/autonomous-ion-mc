"""PSTAR proton stopping powers in liquid water: a committed reference subset.

Values are the NIST PSTAR (ICRU 49, ``I`` = 75 eV) *total* mass stopping
power of protons in liquid water in MeV cm^2/g, transcribed from the file
``Materials/Water/PSTAR_Stop_Pow.dat`` of the Apache-2.0 licensed MCsquare
repository (UCLouvain; https://gitlab.com/openmcsquare/MCsquare, GitHub
mirror https://github.com/e0404/MCsquare) at commit
``211eefe6eaf2b8572d196d17f546f35ffb0ae0cf`` (master on 2026-09-10). Two
independent fetches of that file agreed on every value below. Only the
native PSTAR grid points (four significant figures) are kept; the MCsquare
file also contains interpolated 8-digit values at intermediate energies and a
non-physical entry at 0 MeV, both excluded.

Above 1 MeV the nuclear stopping power of protons in water is below 1e-4 of
the electronic one, so these totals equal the electronic stopping power at
the four quoted significant figures and may be compared directly with an
electronic-only model (decision ``0006``).

Attribution: NIST PSTAR (Standard Reference Database 124, Berger et al.),
redistributed through MCsquare, Universite catholique de Louvain, Apache-2.0.
The 20 numbers below are a de-minimis subset kept for regression testing.
"""

from __future__ import annotations

#: Kinetic energy [MeV] -> total mass stopping power [MeV cm^2/g].
PSTAR_WATER_STOPPING_POWER: dict[float, float] = {
    1.0: 260.8,
    2.0: 158.6,
    5.0: 79.11,
    10.0: 45.67,
    15.0: 32.92,
    20.0: 26.07,
    30.0: 18.76,
    40.0: 14.88,
    50.0: 12.45,
    60.0: 10.78,
    70.0: 9.559,
    80.0: 8.625,
    90.0: 7.888,
    100.0: 7.289,
    150.0: 5.445,
    200.0: 4.492,
    250.0: 3.911,
    300.0: 3.520,
    350.0: 3.241,
    400.0: 3.032,
}

#: Provenance record attached to every comparison that uses this table.
PROVENANCE: dict[str, str] = {
    "dataset": "NIST PSTAR, protons in liquid water (material 276), I = 75 eV",
    "via": "MCsquare Materials/Water/PSTAR_Stop_Pow.dat",
    "repository": "https://github.com/e0404/MCsquare",
    "commit": "211eefe6eaf2b8572d196d17f546f35ffb0ae0cf",
    "license": "Apache-2.0 (UCLouvain); underlying NIST SRD 124 data",
    "quantity": "total mass stopping power, MeV cm^2/g",
    "retrieved": "2026-09-10",
}
