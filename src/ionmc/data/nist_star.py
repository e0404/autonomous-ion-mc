"""Parser for NIST PSTAR/ASTAR plain-text tables (SRD 124).

The text has a title line, a material line, a three-line column header and then
whitespace-separated rows of seven numbers: kinetic energy [MeV], electronic, nuclear and
total stopping power [MeV cm2/g], CSDA range and projected range [g/cm2] and the
detour factor (dimensionless). For ASTAR the energy is the *total* alpha kinetic energy.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from ionmc.data.cache import provenance_of, sha256_bytes

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
    content_sha256: str = ""
    dataset_id: str | None = None
    version: str | None = None
    sha256: str | None = None
    retrieved_at: str | None = None

    def provenance(self) -> dict[str, str | None]:
        """Dataset provenance: dataset_id, version, pinned sha256, content_sha256, retrieved_at."""
        return {
            "dataset_id": self.dataset_id,
            "version": self.version,
            "sha256": self.sha256,
            "content_sha256": self.content_sha256,
            "retrieved_at": self.retrieved_at,
        }


def parse_star_text(text: str) -> StarTable:
    """Parse the text of a PSTAR/ASTAR response into a :class:`StarTable`.

    The result records ``content_sha256`` of the UTF-8 text but no dataset identity (use
    :func:`load_star_table` for verified provenance).
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
    if not np.all(np.isfinite(arr)):
        raise ValueError("table contains NaN or infinity")
    positive = arr[:, [0, 1, 3, 4]]
    if np.any(positive <= 0.0) or np.any(arr[:, [2, 5, 6]] < 0.0):
        raise ValueError(
            "energy, electronic/total stopping power and CSDA range must be positive; "
            "other columns non-negative"
        )
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
        content_sha256=sha256_bytes(text.encode("utf-8")),
    )


def load_star_table(path: str | Path, *, allow_unverified: bool = False) -> StarTable:
    """Read, verify and parse a cached PSTAR/ASTAR file.

    The bytes must hash to the pinned SHA-256 of a registered dataset (``IntegrityError``
    otherwise); the table then carries ``dataset_id``, ``version``, ``sha256``,
    ``content_sha256`` and ``retrieved_at``. With ``allow_unverified=True`` any file is
    parsed, ``dataset_id`` is None and only ``content_sha256`` is recorded.
    """
    path = Path(path)
    raw = path.read_bytes()
    prov = provenance_of(path, raw, allow_unverified)
    table = parse_star_text(raw.decode("ascii"))
    return replace(
        table,
        content_sha256=str(prov["content_sha256"]),
        dataset_id=prov["dataset_id"],
        version=prov["version"],
        sha256=prov["sha256"],
        retrieved_at=prov["retrieved_at"],
    )
