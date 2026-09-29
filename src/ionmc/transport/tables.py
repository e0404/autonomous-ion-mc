"""Pack per-(species, material) physics tables into dense arrays for transport.

The transport step needs, on one common logarithmic grid of kinetic energy
per nucleon t:

* electronic mass stopping power S(t) [MeV cm²/g],
* CSDA range R(t) [g/cm²] as ln R, inverted exactly by the step physics,
* per material: Z/A (for straggling), scattering length ρX_S [g/cm²],
* per species: charge z, mass number A, nuclear mass [MeV].

Everything is float64 here; backends cast as needed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ionmc.materials import Material
from ionmc.physics import stopping
from ionmc.physics.scattering import scattering_length_g_cm2
from ionmc.physics.tables import (
    GRID_POINTS,
    GRID_T_MAX,
    GRID_T_MIN,
    LowEnergySource,
    StoppingTable,
    build_stopping_table,
)
from ionmc.species import SPECIES, Species, get_species


@dataclass(frozen=True)
class TableSet:
    species: tuple[Species, ...]
    materials: tuple[Material, ...]
    t_grid: np.ndarray  # [nt] MeV/u, log-spaced
    stopping: np.ndarray  # [ns, nm, nt] MeV cm²/g
    csda_range: np.ndarray  # [ns, nm, nt] g/cm²
    log_range: (
        np.ndarray
    )  # [ns, nm, nt] ln(csda_range), inverted exactly by the step physics
    z_over_a: np.ndarray  # [nm]
    rho_x_s: np.ndarray  # [nm] g/cm²
    charge: np.ndarray  # [ns] int
    mass_number: np.ndarray  # [ns] int
    mass_mev: np.ndarray  # [ns]
    provenance: dict[str, Any]

    @property
    def t_log_min(self) -> float:
        return float(np.log(self.t_grid[0]))

    @property
    def t_inv_dlog(self) -> float:
        return float((self.t_grid.size - 1) / np.log(self.t_grid[-1] / self.t_grid[0]))

    def species_index(self, species: Species | str) -> int:
        sp = get_species(species)
        return [s.name for s in self.species].index(sp.name)

    def material_index(self, material: Material) -> int:
        return [m.name for m in self.materials].index(material.name)


def build_table_set(
    species: list[Species | str],
    materials: list[Material],
    low_energy: LowEnergySource | None,
) -> TableSet:
    sp = tuple(get_species(s) for s in species)
    mats = tuple(materials)
    t = np.geomspace(GRID_T_MIN, GRID_T_MAX, GRID_POINTS)
    nt = t.size
    stop = np.zeros((len(sp), len(mats), nt))
    rng = np.zeros_like(stop)
    prov: dict[str, Any] = {"tables": {}}
    for i, s in enumerate(sp):
        for j, m in enumerate(mats):
            tbl: StoppingTable = build_stopping_table(s, m, low_energy, grid=t)
            stop[i, j] = tbl.electronic_mev_cm2_g
            rng[i, j] = tbl.csda_range_g_cm2
            prov["tables"][f"{s.name}/{m.name}"] = tbl.provenance
    return TableSet(
        sp,
        mats,
        t,
        stop,
        rng,
        np.log(rng),
        np.array([m.z_over_a for m in mats]),
        np.array([scattering_length_g_cm2(m) for m in mats]),
        np.array([s.z for s in sp], dtype=np.int32),
        np.array([s.a for s in sp], dtype=np.int32),
        np.array([s.mass_mev for s in sp]),
        prov,
    )


def all_species() -> list[str]:
    return sorted(SPECIES)


def effective_charge_factor(species: Species, t_mev_per_u: float) -> float:
    beta, _ = stopping.kinematics(species, np.array([t_mev_per_u]))
    return float((stopping.effective_charge(species.z, beta) / species.z)[0])
