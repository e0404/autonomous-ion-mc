"""Registry of external datasets with pinned SHA-256 hashes.

A hash mismatch on download or import is an integrity failure: the data are rejected.
NIST tables are SRD 124 and must not be committed to Git; they live in the cache. The same holds
for the nuclear data of decision 0041 (ENDF/B-VIII.0 protons, AME2020, EXFOR entries, geant-val):
only the registry entries and the EXFOR manifest (``exfor_manifest.json``, identifiers and hashes,
no values) are committed.

Every dataset carries an evidence ``role`` (``experiment/v3/TELEMETRY.md``): ``construction``,
``calibration``, ``evaluation`` or ``exploratory``. Parser names are identifiers only; this module
imports no parser.

V3-005C (acceptance Amendment 14 (l), step C1b): the elastic-scattering acquisitions are registered
here (EXFOR entries, the arXiv:1409.1938 e-print, six Geant4 11.4.2 source files). ``role`` is the
role recorded by the acquisition tool (the tool accepts only the four roles above); the finer
label of the amendment (e.g. "related-model evidence (shared/related lineage)", "report-only",
"evaluation with shared beam lineage") is in ``lineage`` and in the description, and where the two
differ the label decides how a row is judged. Data finding (decision 0041, rejected alternative
(i)): above about 20 MeV the O-16 MT2 data of the LA150 proton evaluation (MF3 and the MF6 LAW=5
tables) are a numerical copy of C-12; the elastic builder asserts this fail-closed and never uses
the O-16 MT2 as a construction input.
"""

from __future__ import annotations

from dataclasses import dataclass

ROLES = ("construction", "calibration", "evaluation", "exploratory")

_ENDF_LICENSE = "No licence text; free public download from the NNDC (BNL)"
_EXFOR_LICENSE = "CC BY 4.0 (EXFOR, IAEA Nuclear Data Section)"
_GEANT4_CITATION = (
    "Geant4 Collaboration: S. Agostinelli et al., Nucl. Instrum. Meth. A 506 (2003) 250; "
    "J. Allison et al., Nucl. Instrum. Meth. A 835 (2016) 186 (release 11.4.2 source file)"
)
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
    lineage: str = ""

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
        Dataset(
            id="exfor-c0063",
            version="EXFOR master entry C0063",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c0063.txt",
            method="GET",
            post_body=None,
            sha256="717c570c9fa3c75b64f1fcc67bea7a618001539de471296cbf4c81677d7cbea8",
            bytes=55809,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="J.J. Kelly et al., Phys. Rev. C 39 (1989) 1222", entry="C0063"
            ),
            parser="exfor",
            description=(
                "p+16O elastic angular distribution at 135 MeV (V10-A angular dataset). Tool role"
                " evaluation; Amendment 14 label: related-model evidence (shared/related "
                "lineage); pre-1997."
            ),
            role="evaluation",
            lineage="related-model evidence (shared/related lineage); pre-1997",
        ),
        Dataset(
            id="exfor-c0141",
            version="EXFOR master entry C0141",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c0141.txt",
            method="GET",
            post_body=None,
            sha256="1e4db380e24f124311b6969abcdd0a0b24e44442c7a29c2121b56cb96582c3f0",
            bytes=40743,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="J.J. Kelly et al., Phys. Rev. C 41 (1990) 2504", entry="C0141"
            ),
            parser="exfor",
            description=(
                "p+16O elastic angular distribution at about 180 MeV (V10-A angular dataset). "
                "Tool role evaluation; Amendment 14 label: related-model evidence (shared/related"
                " lineage); pre-1997."
            ),
            role="evaluation",
            lineage="related-model evidence (shared/related lineage); pre-1997",
        ),
        Dataset(
            id="exfor-c0148",
            version="EXFOR master entry C0148",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c0148.txt",
            method="GET",
            post_body=None,
            sha256="5ae92bb492d77471163e08b7e223beb62a9b346b038ff929d72857074fde36fe",
            bytes=223155,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="H. Seifert et al., Phys. Rev. C 47 (1993) 1615", entry="C0148"
            ),
            parser="exfor",
            description=(
                "p+16O and p+40Ca elastic angular distributions at about 200 MeV (V10-A angular "
                "dataset). Tool role evaluation; Amendment 14 label: related-model evidence "
                "(shared/related lineage); pre-1997."
            ),
            role="evaluation",
            lineage="related-model evidence (shared/related lineage); pre-1997",
        ),
        Dataset(
            id="exfor-c0550",
            version="EXFOR master entry C0550",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c0550.txt",
            method="GET",
            post_body=None,
            sha256="2c824fd44ec4e419b9846581f08d78c19bf020fc9794415d593768a8d5e57336",
            bytes=12798,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="C.W. Glover et al., Phys. Rev. C 31 (1985) 1", entry="C0550"
            ),
            parser="exfor",
            description=(
                "p+16O elastic angular distribution at 200 MeV (V10-A angular dataset). Tool role"
                " evaluation; Amendment 14 label: related-model evidence (shared/related "
                "lineage); pre-1997."
            ),
            role="evaluation",
            lineage="related-model evidence (shared/related lineage); pre-1997",
        ),
        Dataset(
            id="exfor-c0057",
            version="EXFOR master entry C0057",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c0057.txt",
            method="GET",
            post_body=None,
            sha256="3d02ae236159ec44fc0185e9148a092d506ebfc83961fd6f4c39abc8999748f0",
            bytes=20817,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="H.O. Meyer et al., Phys. Rev. C 23 (1981) 616", entry="C0057"
            ),
            parser="exfor",
            description=(
                "p+12C elastic angular distribution at 200 MeV (V10-A angular dataset). Tool role"
                " evaluation; Amendment 14 label: related-model evidence (shared/related "
                "lineage); pre-1997."
            ),
            role="evaluation",
            lineage="related-model evidence (shared/related lineage); pre-1997",
        ),
        Dataset(
            id="exfor-c0055",
            version="EXFOR master entry C0055",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c0055.txt",
            method="GET",
            post_body=None,
            sha256="4512b9670f9ec0c6e0747e9b506c3604bc90e4ce0155c04455bb52bd30ae70d6",
            bytes=20736,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="H.O. Meyer et al., Phys. Rev. C 27 (1983) 459", entry="C0055"
            ),
            parser="exfor",
            description=(
                "p+12C elastic angular distributions at 120-200 MeV (V10-A angular dataset). Tool"
                " role evaluation; Amendment 14 label: related-model evidence (shared/related "
                "lineage); pre-1997."
            ),
            role="evaluation",
            lineage="related-model evidence (shared/related lineage); pre-1997",
        ),
        Dataset(
            id="exfor-c1420",
            version="EXFOR master entry C1420",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c1420.txt",
            method="GET",
            post_body=None,
            sha256="13364995fbd8a97fabc530ff2fc8f7e8068fa0bc49c5fdf8e3f8a85d0527b6dd",
            bytes=52002,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="H.O. Meyer et al., Phys. Rev. C 37 (1988) 544", entry="C1420"
            ),
            parser="exfor",
            description=(
                "p+12C elastic angular distribution at 250 MeV (V10-A angular dataset). Tool role"
                " evaluation; Amendment 14 label: related-model evidence (shared/related "
                "lineage); pre-1997."
            ),
            role="evaluation",
            lineage="related-model evidence (shared/related lineage); pre-1997",
        ),
        Dataset(
            id="exfor-o0553",
            version="EXFOR master entry O0553",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/o/o0553.txt",
            method="GET",
            post_body=None,
            sha256="01c1259f9e6b13f92e6d765a50e1ebe21a95fb438da15f45ee3e6bcdc9319336",
            bytes=41553,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="C. Rolland et al., Nucl. Phys. 80 (1966) 625", entry="O0553"
            ),
            parser="exfor",
            description=(
                "p+12C and p+40Ca elastic angular distributions at 75 and 150 MeV (V10-A angular "
                "dataset). Tool role evaluation; Amendment 14 label: related-model evidence "
                "(shared/related lineage); pre-1997; no 16O subentry (12C, 40Ca, 209Bi, 140Ce)."
            ),
            role="evaluation",
            lineage=(
                "related-model evidence (shared/related lineage); pre-1997; no 16O "
                "subentry (12C, 40Ca, 209Bi, 140Ce)"
            ),
        ),
        Dataset(
            id="exfor-o0274",
            version="EXFOR master entry O0274",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/o/o0274.txt",
            method="GET",
            post_body=None,
            sha256="3db34b417938f35219903cd12f07f7dd74aa22ed9dce67d142e7b826d3131cab",
            bytes=42444,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="K. Strauch, F. Titus, Phys. Rev. 103 (1956) 200", entry="O0274"
            ),
            parser="exfor",
            description=(
                "p+12C elastic angular distribution at 96 MeV (V10-A angular dataset). Tool role "
                "evaluation; Amendment 14 label: related-model evidence (shared/related lineage);"
                " pre-1997."
            ),
            role="evaluation",
            lineage="related-model evidence (shared/related lineage); pre-1997",
        ),
        Dataset(
            id="exfor-13753",
            version="EXFOR master entry 13753",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/1/13753.txt",
            method="GET",
            post_body=None,
            sha256="73e308154d37c5d9f7bdb98af216de197f8f3ef4ed63bfa7a2b850d6798c1c8f",
            bytes=1147365,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="W.P. Abfalterer et al., Phys. Rev. C 63 (2001) 044608", entry="13753"
            ),
            parser="exfor",
            description=(
                "n+C total cross section 5.3-559 MeV (post-1997 level proxy for the elastic "
                "level, V10-A (c) and D9-BGG). Tool role evaluation; Amendment 14 label: "
                "evaluation (post-1997); gating C level proxy for V10-A (c) (sigma_tot(n+C) minus"
                " sigma_nonel(p+C))."
            ),
            role="evaluation",
            lineage=(
                "evaluation (post-1997); gating C level proxy for V10-A (c) "
                "(sigma_tot(n+C) minus sigma_nonel(p+C))"
            ),
        ),
        Dataset(
            id="exfor-13569",
            version="EXFOR master entry 13569",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/1/13569.txt",
            method="GET",
            post_body=None,
            sha256="82efd34a7f0c0fd3338ed626689ded512a7af45d638071a7dfdd0390593de558",
            bytes=558576,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="R.W. Finlay et al., Phys. Rev. C 47 (1993) 237", entry="13569"
            ),
            parser="exfor",
            description=(
                "n+C and n+O total cross sections (pre-1997; report-only O level proxy). Tool "
                "role exploratory; Amendment 14 label: report-only (pre-1997); O level proxy with"
                " O0579-004."
            ),
            role="exploratory",
            lineage="report-only (pre-1997); O level proxy with O0579-004",
        ),
        Dataset(
            id="exfor-o0579",
            version="EXFOR master entry O0579",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/o/o0579.txt",
            method="GET",
            post_body=None,
            sha256="815bf4f2a65cd68f01020b06ae1707220882238447609bcae6cead9621805562",
            bytes=16767,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="A. Ingemarsson et al., Nucl. Phys. A 653 (1999) 341", entry="O0579"
            ),
            parser="exfor",
            description=(
                "p+12C and p+16O reaction cross sections at 65.5 MeV (post-1997 level proxy). "
                "Tool role evaluation; Amendment 14 label: evaluation (post-1997); C sigma_R "
                "gates the C level (V10-A (c)); O sigma_R is report-only."
            ),
            role="evaluation",
            lineage=(
                "evaluation (post-1997); C sigma_R gates the C level (V10-A "
                "(c)); O sigma_R is report-only"
            ),
        ),
        Dataset(
            id="exfor-o1226",
            version="EXFOR master entry O1226",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/o/o1226.txt",
            method="GET",
            post_body=None,
            sha256="55c760598165cf2a6d89783f40407b9fd0515b2811bdc4f47a2722c6f6383c95",
            bytes=9558,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="M. Mahjour-Shafiei et al., Phys. Rev. C 70 (2004) 024004", entry="O1226"
            ),
            parser="exfor",
            description=(
                "p-p elastic angular distribution at 190 MeV (report-only). Tool role "
                "exploratory; Amendment 14 label: report-only (p-p dsigma/dOmega, V10-pp-c)."
            ),
            role="exploratory",
            lineage="report-only (p-p dsigma/dOmega, V10-pp-c)",
        ),
        Dataset(
            id="exfor-o0145",
            version="EXFOR master entry O0145",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/o/o0145.txt",
            method="GET",
            post_body=None,
            sha256="b839286518f064febdb98131ef2922c9b4d54f5c444ee0857a6f26329f18ec8c",
            bytes=28431,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="M. Avan et al., Phys. Rev. C 30 (1984) 521", entry="O0145"
            ),
            parser="exfor",
            description=(
                "p-p elastic angular distribution at about 200 MeV (report-only). Tool role "
                "exploratory; Amendment 14 label: report-only (p-p dsigma/dOmega, V10-pp-c)."
            ),
            role="exploratory",
            lineage="report-only (p-p dsigma/dOmega, V10-pp-c)",
        ),
        Dataset(
            id="exfor-c2637",
            version="EXFOR master entry C2637",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c2637.txt",
            method="GET",
            post_body=None,
            sha256="3c09bafee1b340ee69c740e798eebe9d5f60b1dd857b51941c86655975ee90fd",
            bytes=4131,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="O. Chamberlain et al., Phys. Rev. 93 (1954) 1424", entry="C2637"
            ),
            parser="exfor",
            description=(
                "p-p elastic angular distribution at 225 MeV (report-only). Tool role "
                "exploratory; Amendment 14 label: report-only (p-p dsigma/dOmega, V10-pp-c)."
            ),
            role="exploratory",
            lineage="report-only (p-p dsigma/dOmega, V10-pp-c)",
        ),
        Dataset(
            id="exfor-c2606",
            version="EXFOR master entry C2606",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/c/c2606.txt",
            method="GET",
            post_body=None,
            sha256="c9a6d4991d375c33937a0629f33091375b78da30a16709fc7e074780200d3d6f",
            bytes=12555,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="O. Chamberlain et al., Phys. Rev. 83 (1951) 923", entry="C2606"
            ),
            parser="exfor",
            description=(
                "p-p elastic angular distributions at 163-250 MeV (report-only). Tool role "
                "exploratory; Amendment 14 label: report-only (p-p dsigma/dOmega, V10-pp-c)."
            ),
            role="exploratory",
            lineage="report-only (p-p dsigma/dOmega, V10-pp-c)",
        ),
        Dataset(
            id="exfor-o2057",
            version="EXFOR master entry O2057",
            url="https://www-nds.iaea.org/nrdc/exfor-master/entry/o/o2057.txt",
            method="GET",
            post_body=None,
            sha256="5691502f3e9635abd1ebd44e43490bf03a7295de06d6ced43e86e0fa3c7c611a",
            bytes=2997,
            license=_EXFOR_LICENSE,
            citation=_EXFOR_CITATION.format(
                experiment="J.M. Cassels, Proc. Phys. Soc. A 69 (1956) 495", entry="O2057"
            ),
            parser="exfor",
            description=(
                "p-p elastic angular distribution at 147 MeV (report-only). Tool role "
                "exploratory; Amendment 14 label: report-only (p-p dsigma/dOmega, V10-pp-c)."
            ),
            role="exploratory",
            lineage="report-only (p-p dsigma/dOmega, V10-pp-c)",
        ),
        Dataset(
            id="geant4-g4barashenkovdata-hh-11.4.2",
            version="Geant4 v11.4.2 G4BarashenkovData.hh",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/processes/hadronic/cross_sections/include/G4BarashenkovData.hh"
            ),
            method="GET",
            post_body=None,
            sha256="0d0102fec187f53ae187ea550181fb1f528913ec496f3a546c2983f1a017034c",
            bytes=19864,
            license="Geant4 Software License",
            citation=_GEANT4_CITATION,
            parser="geant4-source",
            description=(
                "Barashenkov nucleon-nucleus total and reaction cross-section arrays "
                "(sigma_tot(n+A), sigma_inel(p+A)); numbers hand-transcribed into "
                "src/ionmc/physics/elastic. Tool role construction (sigma_el(p+A), S(E) for p-p "
                "above 150 MeV) and the "
                "D9-BGG lineage record; Geant4-derived, numbers hand-transcribed, no code copied."
            ),
            role="construction",
            lineage="construction input (Geant4-derived; shared lineage with TOPAS sigma_el)",
        ),
        Dataset(
            id="geant4-g4bggnucleonelasticxs-cc-11.4.2",
            version="Geant4 v11.4.2 G4BGGNucleonElasticXS.cc",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/processes/hadronic/cross_sections/src/G4BGGNucleonElasticXS.cc"
            ),
            method="GET",
            post_body=None,
            sha256="2e7fb6a0de45ba2a09abb6ef5085513f5acca6c2d6fb0f9e4d463fff5ca41d9c",
            bytes=9165,
            license="Geant4 Software License",
            citation=_GEANT4_CITATION,
            parser="geant4-source",
            description=(
                "BGG nucleon-nucleus elastic cross section (Barashenkov 14 MeV-91 GeV, Coulomb-"
                "factor rule below 14 MeV, 1.0115 x p-p formula for Z=1); formulas re-implemented"
                " by hand. Tool role construction (sigma_el(p+A), S(E) for p-p above 150 MeV) and"
                " the "
                "D9-BGG lineage record; Geant4-derived, numbers hand-transcribed, no code copied."
            ),
            role="construction",
            lineage="construction input (Geant4-derived; shared lineage with TOPAS sigma_el)",
        ),
        Dataset(
            id="geant4-g4nucleonnuclearcrosssection-cc-11.4.2",
            version="Geant4 v11.4.2 G4NucleonNuclearCrossSection.cc",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/processes/hadronic/cross_sections/src/G4NucleonNuclearCrossSection.cc"
            ),
            method="GET",
            post_body=None,
            sha256="99df4ddd49990ee38aeaa00c8588d0b0dfe87c7cd03f864388ff05a0e5a0c392",
            bytes=3967,
            license="Geant4 Software License",
            citation=_GEANT4_CITATION,
            parser="geant4-source",
            description=(
                "Nucleon-nucleus cross sections via the Barashenkov component (sigma_el = "
                "max(sigma_tot - sigma_inel, 0)). Tool role construction (sigma_el(p+A), S(E) for"
                " p-p above 150 MeV) and the "
                "D9-BGG lineage record; Geant4-derived, numbers hand-transcribed, no code copied."
            ),
            role="construction",
            lineage="construction input (Geant4-derived; shared lineage with TOPAS sigma_el)",
        ),
        Dataset(
            id="geant4-g4componentbarnucleonnucleusxsc-cc-11.4.2",
            version="Geant4 v11.4.2 G4ComponentBarNucleonNucleusXsc.cc",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/processes/hadronic/cross_sections/src/G4ComponentBarNucleonNucleusXsc.cc"
            ),
            method="GET",
            post_body=None,
            sha256="a54b06d8eff1db8c18f361e19d4a6606bf105658d211b269fcb87a098dc2ab0c",
            bytes=9721,
            license="Geant4 Software License",
            citation=_GEANT4_CITATION,
            parser="geant4-source",
            description=(
                "Barashenkov component: tabulated Z list, energy grids, linear interpolation in E"
                " and A^(2/3) interpolation in Z. Tool role construction (sigma_el(p+A), S(E) for"
                " p-p above 150 MeV) and the "
                "D9-BGG lineage record; Geant4-derived, numbers hand-transcribed, no code copied."
            ),
            role="construction",
            lineage="construction input (Geant4-derived; shared lineage with TOPAS sigma_el)",
        ),
        Dataset(
            id="geant4-g4hadronnucleonxsc-cc-11.4.2",
            version="Geant4 v11.4.2 G4HadronNucleonXsc.cc",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/processes/hadronic/cross_sections/src/G4HadronNucleonXsc.cc"
            ),
            method="GET",
            post_body=None,
            sha256="5c22c4697d114df24a2be6a020882620f5da9c717c87f4e6fabf5d3ddb755505",
            bytes=39970,
            license="Geant4 Software License",
            citation=_GEANT4_CITATION,
            parser="geant4-source",
            description=(
                "Hadron-nucleon cross sections (p-p elastic 23 + 50 sqrt(ln(0.73/p)^7) mb below "
                "0.73 GeV/c); formula re-implemented by hand. Tool role construction "
                "(sigma_el(p+A), S(E) for p-p above 150 MeV) and the "
                "D9-BGG lineage record; Geant4-derived, numbers hand-transcribed, no code copied."
            ),
            role="construction",
            lineage="construction input (Geant4-derived; shared lineage with TOPAS sigma_el)",
        ),
        Dataset(
            id="geant4-g4hadronelasticphysics-cc-11.4.2",
            version="Geant4 v11.4.2 G4HadronElasticPhysics.cc",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/physics_lists/constructors/hadron_elastic/src/G4HadronElasticPhysics.cc"
            ),
            method="GET",
            post_body=None,
            sha256="aef75b9b96a01c7a1768d392dad19ab9f5d5a8c25ee011af967c06ecd0c4ded8",
            bytes=9311,
            license="Geant4 Software License",
            citation=_GEANT4_CITATION,
            parser="geant4-source",
            description=(
                "Hadron elastic physics constructor: proton elastic = G4BGGNucleonElasticXS + "
                "G4ChipsElasticModel (the lineage record for TOPAS g4h-elastic). Tool role "
                "construction (sigma_el(p+A), S(E) for p-p above 150 MeV) and the "
                "D9-BGG lineage record; Geant4-derived, numbers hand-transcribed, no code copied."
            ),
            role="construction",
            lineage="construction input (Geant4-derived; shared lineage with TOPAS sigma_el)",
        ),
        Dataset(
            id="arxiv-1409-1938-source",
            version="arXiv:1409.1938 e-print (source tarball, https://arxiv.org/src/1409.1938)",
            url="https://arxiv.org/src/1409.1938",
            method="GET",
            post_body=None,
            sha256="92cdaffbc4237e5178e8a83b3ceee6e15f708c2f4ebbd16aff63f249dc990ac2",
            bytes=681055,
            license=(
                "arXiv non-exclusive distribution licence unless stated by the authors "
                "(licence field of the abs page not verified)"
            ),
            citation=(
                "B. Gottschalk et al., arXiv:1409.1938 (2014), e-print source files; table "
                "tbl:dmlg "
                "(measured depth-dose) and Appendix D (model-dependent fit)"
            ),
            parser="latex-source",
            description=(
                "Gottschalk et al. e-print source with the measured log10 dose table (46 depths x "
                "10 radii, 177 MeV beam) for row V6. Tool role evaluation; Amendment 14 label: "
                "evaluation with shared beam lineage (calibration-type inputs: p1 and p5-p7 come "
                "from the paper's fit to the same table)."
            ),
            role="evaluation",
            lineage="evaluation with shared beam lineage (calibration-type inputs)",
        ),
        Dataset(
            id="geant4-g4nuclearradii-cc-11.4.2",
            version="Geant4 v11.4.2 G4NuclearRadii.cc",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/processes/hadronic/util/src/G4NuclearRadii.cc"
            ),
            method="GET",
            post_body=None,
            sha256="53cae34b31621536d2b0a782e19ce649866aa03dfdd67ba56791e2b051d3b18a",
            bytes=8251,
            license="Geant4 Software License",
            citation=_GEANT4_CITATION,
            parser="geant4-source",
            description=(
                "Nuclear radii (RadiusCB with the r0[Z] table) and the Coulomb factor of the BGG "
                "rule below 14 MeV (acquisition 1d8cf519263344b4b8848b2ae869ad3a). Tool role "
                "construction; Geant4-derived, numbers hand-transcribed, "
                "no code copied."
            ),
            role="construction",
            lineage="construction input (Geant4-derived; shared lineage with TOPAS sigma_el)",
        ),
        Dataset(
            id="geant4-g4isotopelist-hh-11.4.2",
            version="Geant4 v11.4.2 G4IsotopeList.hh",
            url=(
                "https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/"
                "source/processes/hadronic/util/include/G4IsotopeList.hh"
            ),
            method="GET",
            post_body=None,
            sha256="d343b8271391891259e09150b6f9d27d0b3a146051b1b17d307bc49f7d439f40",
            bytes=10794,
            license="Geant4 Software License",
            citation=_GEANT4_CITATION,
            parser="geant4-source",
            description=(
                "Effective atomic masses aeff[Z] used by the Barashenkov A^(2/3) interpolation in"
                " Z (acquisition 5c2506b13afe44b1b6db09fd39a767f3). Tool role construction; "
                "Geant4-derived, numbers hand-transcribed, "
                "no code copied."
            ),
            role="construction",
            lineage="construction input (Geant4-derived; shared lineage with TOPAS sigma_el)",
        ),
        Dataset(
            id="pdg-rpp2022-pp-elastic",
            version="PDG rpp2022 pp_elastic.dat",
            url="https://pdg.lbl.gov/2022/hadronic-xsections/rpp2022-pp_elastic.dat",
            method="GET",
            post_body=None,
            sha256="41268d804bd1fc3f165c88a36f055df5f65cfad3c7acdccdf06c3dd442b73dfd",
            bytes=17101,
            license="No licence text; free public download from the PDG (LBNL)",
            citation=(
                "R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01"
            ),
            parser="pdg-xsec",
            description=(
                "PDG compilation of the p-p elastic cross section (acquisition "
                "b742ae4345934c188eb664de5bed355d). Tool role construction: it feeds the "
                "above-150 MeV scaling ratio S(E). It is also the report-only V10 comparison "
                "(shared lineage with S(E))."
            ),
            role="construction",
            lineage="construction input for S(E) and report-only V10 comparison (shared lineage)",
        ),
    )
}


def identify(content_sha256: str) -> Dataset | None:
    """Return the registered dataset whose pinned SHA-256 equals ``content_sha256``, or None."""
    for ds in DATASETS.values():
        if ds.sha256 == content_sha256:
            return ds
    return None
