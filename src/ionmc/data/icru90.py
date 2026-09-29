"""ICRU Report 90 stopping-power tables for liquid water (protons, alphas, carbon).

The values are transcribed from the NIST addendum *Update to ESTAR, PSTAR,
and ASTAR Databases* (tables A.9, A.12 and A.15; mean excitation energy
78 eV), a public NIST document. They are shipped with the package as a
small key-data file (``icru90_water.json``, ≈ 8 kB) because they are the
reference low-energy stopping data for the default water tables and must be
available offline; the file header records the source URL, the PDF SHA-256
and the extraction method.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources

import numpy as np

_SPECIES_KEYS = {"proton": "proton", "he4": "he4", "c12": "c12"}


@dataclass(frozen=True)
class Icru90Table:
    species: str
    mass_number: int
    energy_mev: np.ndarray  # total kinetic energy
    electronic: np.ndarray  # MeV cm²/g
    nuclear: np.ndarray
    total: np.ndarray
    csda_range_g_cm2: np.ndarray
    provenance: dict

    @property
    def energy_per_nucleon_mev(self) -> np.ndarray:
        return self.energy_mev / self.mass_number

    def electronic_at(self, t_per_nucleon: np.ndarray | float) -> np.ndarray:
        """Log-log interpolated electronic mass stopping power (clamped outside)."""
        lt = np.log(
            np.clip(
                np.asarray(t_per_nucleon, dtype=np.float64),
                self.energy_per_nucleon_mev[0],
                self.energy_per_nucleon_mev[-1],
            )
        )
        return np.exp(
            np.interp(lt, np.log(self.energy_per_nucleon_mev), np.log(self.electronic))
        )

    def csda_range_at(self, t_per_nucleon: np.ndarray | float) -> np.ndarray:
        lt = np.log(
            np.clip(
                np.asarray(t_per_nucleon, dtype=np.float64),
                self.energy_per_nucleon_mev[0],
                self.energy_per_nucleon_mev[-1],
            )
        )
        return np.exp(
            np.interp(
                lt, np.log(self.energy_per_nucleon_mev), np.log(self.csda_range_g_cm2)
            )
        )


def _load() -> dict:
    text = resources.files("ionmc.data").joinpath("icru90_water.json").read_text()
    return json.loads(text)


def available_species() -> list[str]:
    return sorted(_load()["tables"])


def provenance() -> dict:
    data = _load()
    return {k: v for k, v in data.items() if k != "tables"}


def water_table(species: str) -> Icru90Table:
    """Return the ICRU 90 liquid-water table for 'proton', 'he4' or 'c12'."""
    data = _load()
    key = _SPECIES_KEYS.get(species.lower())
    if key is None or key not in data["tables"]:
        raise ValueError(
            f"No ICRU 90 water table for species {species!r}; "
            f"available: {available_species()}"
        )
    rows = np.array(data["tables"][key]["rows"], dtype=np.float64)
    return Icru90Table(
        key,
        int(data["tables"][key]["mass_number"]),
        rows[:, 0],
        rows[:, 1],
        rows[:, 2],
        rows[:, 3],
        rows[:, 4],
        provenance(),
    )
