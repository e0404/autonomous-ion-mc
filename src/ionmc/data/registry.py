"""Registered external datasets (immutable versions with checksums).

Every entry pins a Git commit of the source repository and the SHA-256 of the
exact file, so acquisition is reproducible and tamper-evident. Licensing is
recorded per dataset; only data that may legally be redistributed by the
project are registered here (decision ``0007``).
"""

from __future__ import annotations

from ionmc.data.cache import DatasetSpec

_MCSQUARE_COMMIT = "211eefe6eaf2b8572d196d17f546f35ffb0ae0cf"
_MCSQUARE_RAW = f"https://raw.githubusercontent.com/e0404/MCsquare/{_MCSQUARE_COMMIT}"
_MCSQUARE_LICENSE = (
    "Apache-2.0 (MCsquare, Universite catholique de Louvain; attribution to "
    "UCLouvain required by the MCsquare LICENSE file)"
)

#: NIST PSTAR total mass stopping power of protons in liquid water on a
#: 0.5 MeV grid, 0-400 MeV, as redistributed by MCsquare. Rows at the native
#: PSTAR grid energies are the 4-significant-figure PSTAR values (I = 75 eV);
#: the other rows were interpolated by the MCsquare authors; the 0 MeV row is
#: non-physical and dropped by the parser.
MCSQUARE_PSTAR_WATER = DatasetSpec(
    name="mcsquare-pstar-water",
    version=_MCSQUARE_COMMIT,
    url=f"{_MCSQUARE_RAW}/Materials/Water/PSTAR_Stop_Pow.dat",
    sha256="96cd80538ee823f19136a21dceb314d5833111969a507f0697e9289131462c00",
    filename="PSTAR_Stop_Pow.dat",
    license=_MCSQUARE_LICENSE,
    provenance={
        "dataset": "NIST PSTAR (SRD 124) total mass stopping power, protons in "
        "liquid water (material 276), I = 75 eV, MeV cm^2/g",
        "repository": "https://github.com/e0404/MCsquare",
        "path": "Materials/Water/PSTAR_Stop_Pow.dat",
        "commit": _MCSQUARE_COMMIT,
        "grid": "0.0-400.0 MeV in 0.5 MeV steps (801 rows)",
        "note": "underlying NIST SRD 124 data; redistributed via MCsquare",
    },
)

#: Geant4-derived proton stopping power in water on the same grid (MCsquare's
#: alternative table, I = 78 eV physics), used as an independent second data
#: source for the configurable-physics-source requirement.
MCSQUARE_G4_WATER = DatasetSpec(
    name="mcsquare-g4-water",
    version=_MCSQUARE_COMMIT,
    url=f"{_MCSQUARE_RAW}/Materials/Water/G4_Stop_Pow.dat",
    sha256="183f33834b8dcd78642b506c858654b11cc4c9c019d843692e8e41df45635a6b",
    filename="G4_Stop_Pow.dat",
    license=_MCSQUARE_LICENSE,
    provenance={
        "dataset": "Geant4-derived total mass stopping power, protons in liquid "
        "water, MeV cm^2/g (MCsquare G4_Stop_Pow.dat)",
        "repository": "https://github.com/e0404/MCsquare",
        "path": "Materials/Water/G4_Stop_Pow.dat",
        "commit": _MCSQUARE_COMMIT,
        "grid": "0.0-400.0 MeV in 0.5 MeV steps (801 rows)",
    },
)

#: All registered datasets by name.
DATASETS: dict[str, DatasetSpec] = {
    spec.name: spec for spec in (MCSQUARE_PSTAR_WATER, MCSQUARE_G4_WATER)
}
