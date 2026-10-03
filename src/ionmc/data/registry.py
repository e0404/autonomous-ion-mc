"""Registry of external datasets with pinned SHA-256 hashes.

A hash mismatch on download is an integrity failure: the data are rejected.
NIST tables are SRD 124 and must not be committed to Git; they live in the cache.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Dataset:
    """Description of one external dataset (``bytes`` is the exact size in bytes)."""

    id: str
    version: str
    url: str
    method: str
    post_body: str | None
    sha256: str
    bytes: int
    license: str
    citation: str
    parser: str
    description: str


_STAR_URL = "https://physics.nist.gov/cgi-bin/Star/apdata.pl"
_STAR_BODY = (
    "prog={prog}&matno=276&ShowDefault=on&NumofEnergies=0&character=space"
    "&electronic=on&nuclear=on&total=on&csda=on&project=on&detour=on"
)
_STAR_LICENSE = (
    "NIST SRD 124; copyright claimed, attribution required (https://www.nist.gov/open/license)"
)
_STAR_CITATION = (
    "M.J. Berger, J.S. Coursey, M.A. Zucker, J. Chang (2005), ESTAR, PSTAR, and ASTAR: "
    "Computer Programs for Calculating Stopping-Power and Range Tables for Electrons, "
    "Protons, and Helium Ions (version 2.0.1), NIST, https://doi.org/10.18434/T4NC7P"
)

DATASETS: dict[str, Dataset] = {
    d.id: d
    for d in (
        Dataset(
            id="nist-pstar-water-2005",
            version="PSTAR 2.0.1 (2005), liquid water (matno 276)",
            url=_STAR_URL,
            method="POST",
            post_body=_STAR_BODY.format(prog="PSTAR"),
            sha256="f5fc6482235f1e779bb9d221014ed0a343cf28dec88cc1bde19267a38daead18",
            bytes=9308,
            license=_STAR_LICENSE,
            citation=_STAR_CITATION,
            parser="nist_star",
            description="NIST PSTAR proton stopping powers and ranges in liquid water (I = 75 eV).",
        ),
        Dataset(
            id="nist-astar-water-2005",
            version="ASTAR 2.0.1 (2005), liquid water (matno 276)",
            url=_STAR_URL,
            method="POST",
            post_body=_STAR_BODY.format(prog="ASTAR"),
            sha256="dd4a2e504e19beca1f893ca1e2a7d75629bcabc04e002d8e62ef176f3b69b46f",
            bytes=8568,
            license=_STAR_LICENSE,
            citation=_STAR_CITATION,
            parser="nist_star",
            description="NIST ASTAR alpha-particle stopping powers and ranges in liquid water.",
        ),
        Dataset(
            id="geant4-icru90-stopping-11.4.2",
            version="Geant4 v11.4.2 G4ICRU90StoppingData.cc",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/materials/src/G4ICRU90StoppingData.cc"
            ),
            method="GET",
            post_body=None,
            sha256="afd1479497a605d9675c1fe2a3870ddb3f305110fb2d76e13bb600fa3ab281cd",
            bytes=9639,
            license="Geant4 Software License",
            citation=(
                "Geant4 Collaboration: S. Agostinelli et al., Nucl. Instrum. Meth. A 506 (2003) "
                "250; J. Allison et al., Nucl. Instrum. Meth. A 835 (2016) 186. The file embeds "
                "ICRU Report 90 (2016) proton and alpha stopping-power data."
            ),
            parser="icru90",
            description="ICRU 90 proton/alpha electronic stopping arrays (I_water = 78 eV).",
        ),
    )
}
