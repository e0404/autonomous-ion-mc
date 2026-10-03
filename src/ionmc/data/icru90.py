"""Extract ICRU 90 liquid-water electronic stopping arrays from Geant4 source.

``G4ICRU90StoppingData.cc`` (Geant4 v11.4.2) declares ``T0_proton`` (57 values) and
``T0_alpha`` (49 values) energy grids in MeV and, per material, the arrays
``e{i}_proton`` / ``e{i}_alpha`` of electronic stopping power in MeV cm2/g (``AddData``
multiplies by ``MeV cm2/g``). The material index ``i`` follows ``nameNIST_ICRU90 =
{"G4_AIR", "G4_WATER", "G4_GRAPHITE"}``, so water is index 1. The proton grid is the
proton kinetic energy. The alpha grid ends at 1000 MeV, the same span as NIST ASTAR (total
kinetic energy of the alpha particle, up to 250 MeV/u); it is treated as total alpha kinetic
energy, which ``compare_nist.py`` checks against ASTAR on the shared energies.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from ionmc.data.cache import provenance_of, sha256_bytes

WATER_INDEX = 1
_MATERIAL_ORDER = ("G4_AIR", "G4_WATER", "G4_GRAPHITE")


@dataclass(frozen=True)
class Icru90Water:
    """ICRU 90 liquid-water electronic stopping power arrays (float64).

    ``*_energy_mev`` are total projectile kinetic energies [MeV]; ``*_stopping`` are
    mass electronic stopping powers [MeV cm2/g].
    """

    proton_energy_mev: NDArray[np.float64]
    proton_stopping: NDArray[np.float64]
    alpha_energy_mev: NDArray[np.float64]
    alpha_stopping: NDArray[np.float64]
    content_sha256: str = ""
    dataset_id: str | None = None
    version: str | None = None
    sha256: str | None = None
    retrieved_at: str | None = None


def _array(source: str, name: str) -> NDArray[np.float64]:
    """Return the numbers of the C++ array ``name[...] = { ... };`` as float64."""
    match = re.search(rf"{name}\[\d+\]\s*=\s*\{{([^}}]*)\}}", source)
    if match is None:
        raise ValueError(f"array {name} not found in source")
    tokens = [t.strip().rstrip("f") for t in match.group(1).split(",") if t.strip()]
    return np.asarray([float(t) for t in tokens], dtype=np.float64)


def parse_icru90_source(source: str, material: str = "G4_WATER") -> Icru90Water:
    """Parse the Geant4 source text and return the arrays of ``material``.

    Raises ``ValueError`` if arrays are missing or their lengths differ from the grids.
    """
    index = _MATERIAL_ORDER.index(material)
    t_p = _array(source, "T0_proton")
    t_a = _array(source, "T0_alpha")
    s_p = _array(source, f"e{index}_proton")
    s_a = _array(source, f"e{index}_alpha")
    if s_p.shape != t_p.shape or s_a.shape != t_a.shape:
        raise ValueError("stopping arrays and energy grids differ in length")
    for arr in (t_p, t_a, s_p, s_a):
        if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
            raise ValueError("energies and stopping powers must be finite and positive")
    if np.any(np.diff(t_p) <= 0.0) or np.any(np.diff(t_a) <= 0.0):
        raise ValueError("energy grids are not strictly increasing")
    return Icru90Water(t_p, s_p, t_a, s_a, content_sha256=sha256_bytes(source.encode("utf-8")))


def load_icru90_water(path: str | Path, *, allow_unverified: bool = False) -> Icru90Water:
    """Read, verify and parse the cached Geant4 source file; return the liquid-water arrays.

    The bytes must hash to the pinned SHA-256 of a registered dataset (``IntegrityError``
    otherwise) and the result carries ``dataset_id``, ``version``, ``sha256``,
    ``content_sha256`` and ``retrieved_at``; ``allow_unverified=True`` parses any file
    and records only ``content_sha256``.
    """
    path = Path(path)
    raw = path.read_bytes()
    prov = provenance_of(path, raw, allow_unverified)
    data = parse_icru90_source(raw.decode("utf-8"))
    return replace(
        data,
        content_sha256=str(prov["content_sha256"]),
        dataset_id=prov["dataset_id"],
        version=prov["version"],
        sha256=prov["sha256"],
        retrieved_at=prov["retrieved_at"],
    )
