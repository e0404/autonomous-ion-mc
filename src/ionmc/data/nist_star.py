"""Parser for NIST PSTAR/ASTAR plain-text tables (SRD 124).

The text has a title line, a material line, a three-line column header and then
whitespace-separated rows of seven numbers: kinetic energy [MeV], electronic, nuclear and
total stopping power [MeV cm2/g], CSDA range and projected range [g/cm2] and the
detour factor (dimensionless). For ASTAR the energy is the *total* alpha kinetic energy.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

N_COLUMNS = 7


@dataclass(frozen=True)
class StarTable:
    """A NIST STAR table; all arrays are float64 with one entry per energy row.

    Units: ``energy_mev`` [MeV, total kinetic energy of the projectile]; ``s_electronic``,
    ``s_nuclear``, ``s_total`` [MeV cm2/g]; ``csda_range``, ``projected_range`` [g/cm2];
    ``detour`` dimensionless.
    """

    program: str
    material: str
    energy_mev: NDArray[np.float64]
    s_electronic: NDArray[np.float64]
    s_nuclear: NDArray[np.float64]
    s_total: NDArray[np.float64]
    csda_range: NDArray[np.float64]
    projected_range: NDArray[np.float64]
    detour: NDArray[np.float64]


def parse_star_text(text: str) -> StarTable:
    """Parse the text of a PSTAR/ASTAR response into a :class:`StarTable`.

    Raises ``ValueError`` if no rows are found, a row does not have seven numeric columns,
    or the energies are not strictly increasing.
    """
    lines = [ln.rstrip() for ln in text.splitlines()]
    nonblank = [ln for ln in lines if ln.strip()]
    if len(nonblank) < 2:
        raise ValueError("STAR text too short")
    program = nonblank[0].split(":")[0].strip()
    material = nonblank[1].strip()
    rows: list[list[float]] = []
    for ln in lines:
        parts = ln.split()
        if not parts:
            continue
        try:
            values = [float(p) for p in parts]
        except ValueError:
            if rows:
                raise ValueError(f"non-numeric line after data start: {ln!r}") from None
            continue  # header
        if len(values) != N_COLUMNS:
            raise ValueError(f"expected {N_COLUMNS} columns, got {len(values)}: {ln!r}")
        rows.append(values)
    if not rows:
        raise ValueError("no data rows found")
    arr = np.asarray(rows, dtype=np.float64)
    if np.any(np.diff(arr[:, 0]) <= 0.0):
        raise ValueError("energies are not strictly increasing")
    return StarTable(
        program=program,
        material=material,
        energy_mev=arr[:, 0].copy(),
        s_electronic=arr[:, 1].copy(),
        s_nuclear=arr[:, 2].copy(),
        s_total=arr[:, 3].copy(),
        csda_range=arr[:, 4].copy(),
        projected_range=arr[:, 5].copy(),
        detour=arr[:, 6].copy(),
    )


def load_star_table(path: str | Path) -> StarTable:
    """Read and parse a cached PSTAR/ASTAR file."""
    return parse_star_text(Path(path).read_text(encoding="ascii"))
