"""NIST PSTAR/ASTAR stopping-power tables (NIST SRD 124, ICRU Report 49)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ionmc.data.cache import DataCache, DatasetRecord, DatasetSpec

ENDPOINT = "https://physics.nist.gov/cgi-bin/Star/apdata.pl"
VERSION = "ICRU49"
LICENSE_BASIS = (
    "NIST Standard Reference Database 124, publicly accessible US government "
    "data; cached locally, not redistributed"
)
CITATION = (
    "M.J. Berger, J.S. Coursey, M.A. Zucker, J. Chang, NIST Standard Reference "
    "Database 124: Stopping-Power and Range Tables for Electrons, Protons, and "
    "Helium Ions (2017), https://dx.doi.org/10.18434/T4NC7P"
)

# Only codes known with certainty are listed. Additional material codes must
# be verified against https://physics.nist.gov/PhysRefData/Star/Text/PSTAR-t.html
# before being added (e.g. skeletal muscle is deliberately omitted).
MATERIAL_NUMBERS: dict[str, str] = {
    # Verified 2026-09-30 by querying apdata.pl and reading the returned title line.
    "hydrogen": "001",
    "carbon_amorphous": "006",
    "nitrogen": "007",
    "oxygen": "008",
    "aluminum": "013",
    "silicon": "014",
    "argon": "018",
    "titanium": "022",
    "iron": "026",
    "copper": "029",
    "lead": "082",
    "a150_tissue_plastic": "099",
    "adipose_tissue_icrp": "103",
    "air_dry": "104",
    "bone_compact_icru": "119",
    "bone_cortical_icrp": "120",
    "polyethylene": "221",
    "pmma": "223",
    "polystyrene": "226",
    "water_liquid": "276",
    "water_vapor": "277",
    "graphite": "906",
}
# Additional codes must be verified against
# https://physics.nist.gov/PhysRefData/Star/Text/PSTAR-t.html before use; calcium
# and several ICRP tissues (muscle, soft tissue, lung) have no PSTAR entry.

# ionmc material name -> NIST STAR material key (None: no STAR table)
_IONMC_TO_STAR: dict[str, str] = {
    "water": "water_liquid",
    "air": "air_dry",
    "graphite": "graphite",
    "pmma": "pmma",
    "polyethylene": "polyethylene",
    "adipose_tissue": "adipose_tissue_icrp",
    "bone_compact": "bone_compact_icru",
    "bone_cortical": "bone_cortical_icrp",
    "aluminium": "aluminum",
    "lead": "lead",
    "titanium": "titanium",
}


def material_key(material: object) -> str | None:
    """NIST STAR material key for an ``ionmc.materials.Material`` (or its name)."""
    name = getattr(material, "name", material)
    return _IONMC_TO_STAR.get(str(name))


_NAMES = ("PSTAR", "ASTAR")


def spec(program: str, material: str) -> DatasetSpec:
    """Build the :class:`DatasetSpec` for ``program`` ("PSTAR"/"ASTAR")."""
    prog = program.upper()
    if prog not in _NAMES:
        raise ValueError(f"program must be PSTAR or ASTAR, got {program!r}")
    if material not in MATERIAL_NUMBERS:
        raise KeyError(
            f"unknown material {material!r}; known: {sorted(MATERIAL_NUMBERS)}"
        )
    fields = {
        "prog": prog,
        "matno": MATERIAL_NUMBERS[material],
        "ShowDefault": "on",
        "NumofEnergies": "0",
        "character": "space",
        "electronic": "on",
        "nuclear": "on",
        "total": "on",
        "csda": "on",
        "project": "on",
        "detour": "on",
    }
    return DatasetSpec(
        dataset_id=f"nist-{prog.lower()}/{material}",
        version=VERSION,
        url=ENDPOINT,
        method="POST",
        post_fields=fields,
        license_basis=LICENSE_BASIS,
        citation=CITATION,
        description=f"NIST {prog} stopping powers and ranges, {material}",
    )


@dataclass
class StarTable:
    """Parsed PSTAR/ASTAR table.

    Energies are in MeV (total kinetic energy of the ion); stopping powers in
    MeV cm2/g; ranges in g/cm2; detour factor is dimensionless.
    """

    program: str
    material: str
    energy_mev: NDArray[np.float64]
    electronic: NDArray[np.float64]
    nuclear: NDArray[np.float64]
    total: NDArray[np.float64]
    csda_range_g_cm2: NDArray[np.float64]
    projected_range_g_cm2: NDArray[np.float64]
    detour_factor: NDArray[np.float64]

    def energy_per_nucleon_mev(self) -> NDArray[np.float64]:
        """Kinetic energy per nucleon [MeV/u] (ASTAR total energy divided by 4)."""
        if self.program.upper() == "ASTAR":
            return self.energy_mev / 4.0
        return self.energy_mev

    def electronic_at(self, t_per_nucleon: float | NDArray[np.float64]) -> object:
        """Electronic stopping power [MeV cm2/g] by log-log interpolation.

        ``t_per_nucleon`` in MeV/u; values outside the table raise ValueError.
        """
        t = np.asarray(t_per_nucleon, dtype=np.float64)
        x = self.energy_per_nucleon_mev()
        if np.any(t < x[0]) or np.any(t > x[-1]):
            raise ValueError(f"energy outside table range [{x[0]:g}, {x[-1]:g}] MeV/u")
        y = np.exp(np.interp(np.log(t), np.log(x), np.log(self.electronic)))
        return float(y) if y.ndim == 0 else y


def parse_star_table(
    text: str, program: str = "PSTAR", material: str = ""
) -> StarTable:
    """Parse the plain-text table returned by the NIST STAR CGI.

    A data line has 7 whitespace-separated tokens whose first parses as float;
    all other (header) lines are skipped.
    """
    rows: list[list[float]] = []
    for line in text.splitlines():
        tok = line.split()
        if len(tok) != 7:
            continue
        try:
            rows.append([float(t) for t in tok])
        except ValueError:
            continue
    if not rows:
        raise ValueError("no data lines found in STAR table text")
    arr = np.asarray(rows, dtype=np.float64)
    return StarTable(
        program=program.upper(),
        material=material,
        energy_mev=arr[:, 0].copy(),
        electronic=arr[:, 1].copy(),
        nuclear=arr[:, 2].copy(),
        total=arr[:, 3].copy(),
        csda_range_g_cm2=arr[:, 4].copy(),
        projected_range_g_cm2=arr[:, 5].copy(),
        detour_factor=arr[:, 6].copy(),
    )


def load_star_table(
    cache: DataCache, program: str, material: str, *, offline: bool = False
) -> tuple[StarTable, DatasetRecord]:
    """Fetch (cache-first) and parse a STAR table with its provenance record."""
    data, rec = cache.fetch(spec(program, material), offline=offline)
    table = parse_star_table(data.decode("utf-8"), program, material)
    return table, rec
