"""Bounded carbon-12 nuclear fragmentation and the distal dose tail (decision 0029).

Carbon ions fragment strongly; the lighter charged fragments (H, He, Li-B) have a
longer range than the primary carbon at the same velocity and deposit a dose tail
beyond the Bragg peak. This module produces that tail with a deterministic model
that reuses the existing multi-ion CSDA transport (decisions 0027/0028) and needs
no transport-kernel change:

* the primary carbon survival is ``S(z) = exp(-Sigma*z)`` with a constant reaction
  cross-section ``sigma_R ~ 1.4 barn`` per water molecule;
* the reactions in each depth bin emit, per representative species (proton, alpha,
  boron-11), a forward fragment at the same velocity ``E_f = A_f*(E_C/12)``,
  transported by that species' z^2-scaled table into the depth-dose grid.

The total depth dose is the attenuated primary plus the summed fragment doses.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ionmc.constants import AVOGADRO
from ionmc.data.stopping_tables import (
    StoppingTable,
    scale_ion_stopping_table,
)
from ionmc.particles import ALPHA, BORON_11, CARBON_12, PROTON, Particle
from ionmc.physics.nuclear import BARN_TO_CM2
from ionmc.transport.depth_dose import DepthDoseGrid
from ionmc.transport.engine import TransportEngine
from ionmc.transport.geometry import WaterSlab
from ionmc.transport.source import PencilBeamSource
from ionmc.transport.state import ParticleState, Species, Status

#: Carbon-12 total reaction cross-section in water [barn per water molecule],
#: energy-independent over the therapeutic range (decision 0029).
CARBON_SIGMA_R_BARN: float = 1.4
#: Molar mass of water [g/mol].
WATER_G_PER_MOL: float = 18.01528
#: Representative charged fragments and their per-reaction multiplicities.
FRAGMENT_SPECIES: tuple[tuple[Particle, float], ...] = (
    (PROTON, 2.0),
    (ALPHA, 0.7),
    (BORON_11, 0.35),
)


def macroscopic_carbon_reaction_per_cm(density_g_per_cm3: float = 1.0) -> float:
    """Macroscopic carbon reaction rate ``Sigma`` [1/cm] in water (decision 0029)."""
    n_mol = AVOGADRO * density_g_per_cm3 / WATER_G_PER_MOL
    return n_mol * CARBON_SIGMA_R_BARN * BARN_TO_CM2


@dataclass
class FragmentationResult:
    """Carbon depth dose decomposed into the attenuated primary and the fragment
    tail (decision 0029). All energy arrays are [MeV per depth bin]."""

    grid: DepthDoseGrid
    primary_edep_mev: np.ndarray  # attenuated primary (survivors)
    fragment_edep_mev: np.ndarray  # summed over fragment species
    per_species_edep_mev: dict[str, np.ndarray] = field(default_factory=dict)
    energy_in_mev: float = 0.0
    escaped_mev: float = 0.0
    primary_survival_at_peak: float = float("nan")

    @property
    def total_edep_mev(self) -> np.ndarray:
        return self.primary_edep_mev + self.fragment_edep_mev

    def tail_to_peak(self, distal_mm: float = 10.0) -> float:
        """Ratio of the fragment-tail dose a distance ``distal_mm`` beyond the total
        Bragg peak to the peak dose (the fragment-tail magnitude metric)."""
        total = self.total_edep_mev
        centers = self.grid.centers_mm
        kpk = int(total.argmax())
        peak = float(total[kpk])
        if peak <= 0.0:
            return float("nan")
        z_tail = centers[kpk] + distal_mm
        k_tail = int(np.argmin(np.abs(centers - z_tail)))
        return float(total[k_tail] / peak)


def carbon_fragmentation_depth_dose(
    proton_table: StoppingTable,
    slab: WaterSlab,
    grid: DepthDoseGrid,
    source: PencilBeamSource,
    n_histories: int = 1,
    seed: int = 12345,
    path: str = "python",
    device: str = "cpu",
    straggling: bool = True,
) -> FragmentationResult:
    """Carbon depth dose with the distal fragment tail (decision 0029).

    ``source`` is the carbon primary beam (its ``energy_mev`` is the total carbon
    energy, e.g. 3480 MeV = 290 MeV/u). ``proton_table`` is the base proton stopping
    table from which the carbon and fragment tables are z^2-scaled. Returns the
    attenuated primary dose, the summed and per-species fragment doses, and the
    energy bookkeeping.
    """
    carbon_table = scale_ion_stopping_table(proton_table, CARBON_12)
    eng = TransportEngine(
        carbon_table,
        slab,
        grid,
        particle=CARBON_12,
        straggling=straggling,
        nuclear=False,
    )
    primary = eng.run(source, n_histories, seed=seed, path=path, device=device)
    d0 = primary.edep_mev
    energy_in = primary.energy_in_mev

    # -- carbon survival and per-bin reaction weights ------------------------
    sigma = macroscopic_carbon_reaction_per_cm(1.0)  # homogeneous water
    centers_cm = grid.centers_mm / 10.0
    edges_cm = grid.edges_mm / 10.0
    s_center = np.exp(-sigma * centers_cm)
    s_edge = np.exp(-sigma * edges_cm)
    d_primary = d0 * s_center  # attenuated primary (survivors deposit dose)
    # weighted number of reactions in each bin (proportional to the primary weight)
    total_weight = float(np.sum(source.sample(n_histories, seed).weight))
    react_per_bin = total_weight * (s_edge[:-1] - s_edge[1:])

    # -- residual carbon energy at each production depth (range inversion) ----
    r_total = float(
        np.interp(
            source.energy_mev,
            carbon_table.energy_mev,
            carbon_table.csda_range_g_per_cm2,
        )
    )
    resid_range = r_total - centers_cm  # water density 1 g/cm^3 -> cm == g/cm^2
    e_carbon = np.where(
        resid_range > carbon_table.csda_range_g_per_cm2[0],
        np.interp(
            resid_range,
            carbon_table.csda_range_g_per_cm2,
            carbon_table.energy_mev,
        ),
        0.0,
    )

    # -- per-species fragment transport --------------------------------------
    fragment_total = grid.empty()
    per_species: dict[str, np.ndarray] = {}
    a_c = float(CARBON_12.mass_number)
    for species, mult in FRAGMENT_SPECIES:
        # production bins with a defined residual energy and a reaction weight
        active = (e_carbon > 0.0) & (react_per_bin > 0.0)
        if not np.any(active):
            per_species[species.name] = grid.empty()
            continue
        a_f = float(species.mass_number)
        z_prod = grid.centers_mm[active]
        e_frag = a_f * e_carbon[active] / a_c  # same velocity: E_f = A_f*(E_C/12)
        w_frag = mult * react_per_bin[active]
        # a fragment "history" per production bin (deterministic weighted state)
        n_frag = int(z_prod.shape[0])
        state = ParticleState.allocate(n_frag)
        state.position_mm[:, 2] = z_prod
        state.direction[:] = np.array([0.0, 0.0, 1.0])
        state.energy_mev[:] = e_frag
        state.weight[:] = w_frag
        state.species[:] = int(Species.PROTON)
        state.status[:] = int(Status.ALIVE)
        state.rng_state[:] = np.arange(1, n_frag + 1, dtype=np.uint32)
        frag_table = scale_ion_stopping_table(proton_table, species)
        frag_eng = TransportEngine(
            frag_table,
            slab,
            grid,
            particle=species,
            straggling=straggling,
            nuclear=False,
        )
        s_edep, _, _, _, _, _, _, _ = frag_eng._transport(state, path, device)
        per_species[species.name] = s_edep
        fragment_total = fragment_total + s_edep

    total = d_primary + fragment_total
    escaped = energy_in - float(np.sum(total))
    survival_at_peak = float(s_center[int(d0.argmax())])
    return FragmentationResult(
        grid=grid,
        primary_edep_mev=d_primary,
        fragment_edep_mev=fragment_total,
        per_species_edep_mev=per_species,
        energy_in_mev=energy_in,
        escaped_mev=escaped,
        primary_survival_at_peak=survival_at_peak,
    )
