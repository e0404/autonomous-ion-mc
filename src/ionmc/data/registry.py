"""Registry of external datasets with pinned SHA-256 hashes.

A hash mismatch on download or import is an integrity failure: the data are rejected.
NIST tables are SRD 124 and must not be committed to Git; they live in the cache. The same holds
for the nuclear data of decision 0041 (ENDF/B-VIII.0 protons, AME2020, EXFOR entries, geant-val):
only the registry entries and the EXFOR manifest (``exfor_manifest.json``, identifiers and hashes,
no values) are committed.

Every dataset carries an evidence ``role`` (``experiment/v3/TELEMETRY.md``): ``construction``,
``calibration``, ``evaluation`` or ``exploratory``. Parser names are identifiers only; this module
imports no parser.
"""

from __future__ import annotations

from dataclasses import dataclass

ROLES = ("construction", "calibration", "evaluation", "exploratory")

_ENDF_LICENSE = "No licence text; free public download from the NNDC (BNL)"
_EXFOR_LICENSE = "CC BY 4.0 (EXFOR, IAEA Nuclear Data Section)"
_EXFOR_CITATION = (
    "{experiment}; N. Otuka et al., Nucl. Data Sheets 120 (2014) 272 (EXFOR). EXFOR master entry "
    "{entry}"
)


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
    role: str = "exploratory"

    def __post_init__(self) -> None:
        if self.role not in ROLES:
            raise ValueError(f"{self.id}: role {self.role!r} is not one of {ROLES}")


_STAR_URL = "https://physics.nist.gov/cgi-bin/Star/apdata.pl"
_STAR_BODY = (
    "prog={prog}&matno=276&ShowDefault=on&NumofEnergies=0&character=space"
    "&electronic=on&nuclear=on&total=on&csda=on&project=on&detour=on"
)
_STAR_LICENSE = (
    "NIST SRD 124: copyright secured under 15 U.S.C. 290e, all rights reserved "
    "(https://www.nist.gov/open/license); use-only, downloaded by each user, never redistributed"
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
            role="evaluation",
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
            role="evaluation",
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
            role="evaluation",
        ),
        Dataset(
            id="endf-b8.0-protons",
            version="ENDF/B-VIII.0 proton sublibrary (ENDF-B-VIII.0_protons.zip)",
            url="https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_protons.zip",
            method="GET",
            post_body=None,
            sha256="27bcafb89cf0444c53c6b9f3dd17618c62e2e6b3694f70d31c4f502399f103f7",
            bytes=14112244,
            license=_ENDF_LICENSE,
            citation=(
                "D.A. Brown et al., Nucl. Data Sheets 148 (2018) 1 (ENDF/B-VIII.0); "
                "M.B. Chadwick et al., Nucl. Sci. Eng. 131 (1999) 293 (LA150 evaluations)"
            ),
            parser="endf6",
            description=(
                "ENDF/B-VIII.0 proton sublibrary (LA150 evaluations up to 150 MeV): MF3/MT5 "
                "non-elastic cross sections and MF6/MT5 product yields and spectra."
            ),
            role="construction",
        ),
        Dataset(
            id="ame2020-mass",
            version="AME2020 mass table (mass_1.mas20.txt)",
            url="https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt",
            method="GET",
            post_body=None,
            sha256="e8599c6d7f724fac91934e59f1b9de8fb8f63e820f4b39456b790665ed2a3307",
            bytes=472648,
            license="No licence text; free public download from the AMDC (IAEA NDS)",
            citation="M. Wang et al., Chin. Phys. C 45 (2021) 030003 (AME2020)",
            parser="ame2020",
            description="AME2020 atomic masses and mass excesses (mass_1.mas20.txt).",
            role="construction",
        ),
        Dataset(
            id="exfor-d0356",
            version="EXFOR master entry D0356",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/d/d0356.txt",
            method="GET",
            post_body=None,
            sha256="2ef17fb10aaaeb5096dc92f9c3175f51f2252163eb4a83aa15669027371747bf",
            bytes=12798,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="A. Auce et al., Phys. Rev. C 71 (2005) 064606", entry="D0356"
            ),
            parser="exfor",
            description=(
                "Proton reaction cross sections on 12C, 40Ca, 58Ni, 90Zr, 208Pb (published 2005); "
                "evaluation data for the informative row V1b."
            ),
            role="evaluation",
        ),
        Dataset(
            id="exfor-c1862",
            version="EXFOR master entry C1862",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c1862.txt",
            method="GET",
            post_body=None,
            sha256="157aca7e1b8f6cb6814a5fa99fa749823721af6829cbf84712a56fce159ae4ef",
            bytes=18063,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="I. Slaus et al., Phys. Rev. C 12 (1975) 1093", entry="C1862"
            ),
            parser="exfor",
            description=(
                "Proton non-elastic cross sections on 12C, 9Be, 16O, Si (published 1975, before "
                "1997: shares lineage with LA150); exploratory, report-only."
            ),
            role="exploratory",
        ),
        Dataset(
            id="geant-val-exfor-inelastic-7",
            version="geant-val experiment curves of article inspireId -7",
            url="https://geant-val.cern.ch/api/getExpPlotsByInspireId?inspire_id=-7",
            method="GET",
            post_body=None,
            sha256="fa7ac90fd259f728e5948c71b7a3636bec7b7d9756abeb60837c7083dc2a4c49",
            bytes=12233,
            license="No licence stated",
            citation=(
                "Geant4 validation portal (geant-val), EXFOR-derived experimental inelastic "
                "cross-section curves (article inspireId -7)"
            ),
            parser="geant_val_json",
            description=(
                "geant-val JSON export of EXFOR-derived proton inelastic cross-section curves; "
                "exploratory, report-only; acquired through the provenance tool and staged "
                "with `ionmc data import`."
            ),
            role="exploratory",
        ),
    )
}


def identify(content_sha256: str) -> Dataset | None:
    """Return the registered dataset whose pinned SHA-256 equals ``content_sha256``, or None."""
    for ds in DATASETS.values():
        if ds.sha256 == content_sha256:
            return ds
    return None
