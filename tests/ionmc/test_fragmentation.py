"""Tests of the bounded carbon-12 fragmentation dose tail (decision 0029).

``carbon_fragmentation_depth_dose`` attenuates the primary carbon by a constant
reaction cross-section and emits forward, same-velocity fragments (proton, alpha,
boron-11) that produce the distal dose tail beyond the Bragg peak. The tests cover
the cross-section, the tail magnitude (via the resolution-robust integrated
distal-dose fraction) and reach, primary survival, the fragment energy budget
(reconciled independently, not tautologically), and reference/Warp CPU agreement.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.constants import AVOGADRO
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table, scale_ion_stopping_table
from ionmc.fragmentation import (
    CARBON_SIGMA_R_BARN,
    FRAGMENT_SPECIES,
    carbon_fragmentation_depth_dose,
    macroscopic_carbon_reaction_per_cm,
)
from ionmc.particles import CARBON_12
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine
from ionmc.transport.geometry import WaterSlab

CARBON_RANGE_MM = 163.0  # 290 MeV/u carbon range in water


@pytest.fixture(scope="module")
def proton_table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.fixture(scope="module")
def frag(proton_table):
    grid = DepthDoseGrid(550.0, 275)  # 0-55 cm, 2 mm bins
    return carbon_fragmentation_depth_dose(
        proton_table,
        WaterSlab(550.0),
        grid,
        PencilBeamSource(3480.0),
        n_histories=1,
        seed=1,
        path="python",
        straggling=False,
    )


def test_macroscopic_cross_section() -> None:
    """Sigma = n_mol * sigma_R, ~0.047 /cm for 1.4 barn per water molecule."""
    n_mol = AVOGADRO * 1.0 / 18.01528
    expected = n_mol * CARBON_SIGMA_R_BARN * 1.0e-24
    assert macroscopic_carbon_reaction_per_cm(1.0) == pytest.approx(expected)
    assert 0.04 < macroscopic_carbon_reaction_per_cm(1.0) < 0.05


def test_fragment_tail_present(frag) -> None:
    """The integrated distal-dose fraction (the resolution-robust tail metric) is a
    substantial fraction of the deposited dose: > 10 % lands more than 10 mm distal
    to the Bragg peak, tapering with margin. Canonical ~16 %."""
    assert 0.10 <= frag.distal_dose_fraction(10.0) <= 0.25
    # monotone decrease with distance from the peak (a real tail, not a plateau)
    assert (
        frag.distal_dose_fraction(5.0)
        > frag.distal_dose_fraction(10.0)
        > frag.distal_dose_fraction(20.0)
        > 0.10
    )


def test_tail_metric_is_resolution_robust(proton_table) -> None:
    """The gate metric (integrated distal-dose fraction) is invariant to the bin
    width, unlike the single-bin ``tail_to_peak`` point ratio (which is why the
    latter is a diagnostic, not the gate)."""
    fracs = []
    for n_bins in (138, 275, 550, 1100):  # 4, 2, 1, 0.5 mm bins
        r = carbon_fragmentation_depth_dose(
            proton_table,
            WaterSlab(550.0),
            DepthDoseGrid(550.0, n_bins),
            PencilBeamSource(3480.0),
            n_histories=1,
            seed=1,
            path="python",
            straggling=False,
        )
        fracs.append(r.distal_dose_fraction(10.0))
    assert max(fracs) - min(fracs) < 0.01  # ~0.164 across an 8x resolution range


def test_primary_only_has_no_tail(proton_table) -> None:
    """A primary-only carbon run (no fragmentation) has ~0 dose beyond the peak, so
    the distal-dose-fraction metric isolates the tail as a fragmentation effect."""
    grid = DepthDoseGrid(550.0, 275)
    carbon = scale_ion_stopping_table(proton_table, CARBON_12)
    eng = TransportEngine(
        carbon,
        WaterSlab(550.0),
        grid,
        particle=CARBON_12,
        straggling=False,
        nuclear=False,
    )
    r = eng.run(PencilBeamSource(3480.0), n_histories=1, seed=1, path="python")
    dd = r.edep_mev
    centers = grid.centers_mm
    kpk = int(dd.argmax())
    primary_distal = float(dd[centers > centers[kpk] + 10.0].sum() / dd.sum())
    assert primary_distal < 1e-3


def test_primary_survival(frag) -> None:
    """The surviving primary fraction at the peak matches exp(-Sigma*R) ~ 0.47."""
    sigma = macroscopic_carbon_reaction_per_cm(1.0)
    expected = float(np.exp(-sigma * CARBON_RANGE_MM / 10.0))
    assert frag.primary_survival_at_peak == pytest.approx(expected, rel=0.05)
    assert 0.4 <= frag.primary_survival_at_peak <= 0.55


def test_tail_reach(frag) -> None:
    """Fragment dose stays > 1 % of the peak out to >= 1.5x the carbon range."""
    total = frag.total_edep_mev
    centers = frag.grid.centers_mm
    kpk = int(total.argmax())
    significant = np.where(total > 0.01 * total[kpk])[0]
    reach_mm = centers[significant[-1]]
    assert reach_mm >= 1.5 * CARBON_RANGE_MM


def test_fragment_energy_budget_reconciles(frag, proton_table) -> None:
    """The fragment energy budget is checked three ways, not by construction:

    1. the injected fragment KE, reconstructed *independently* from the reaction
       weights, multiplicities and residual carbon energy, matches the value the
       transport recorded;
    2. that injected KE is deposited in full (the grid contains the fragments), so
       ``fragment_edep.sum() == fragment_energy_injected``;
    3. only then does ``deposited + escaped == E0`` (with escaped >= 0) close.
    """
    # (1) independent reconstruction of Sum_species Sum_bins w_frag * E_frag
    grid = frag.grid
    carbon = scale_ion_stopping_table(proton_table, CARBON_12)
    sigma = macroscopic_carbon_reaction_per_cm(1.0)
    s_edge = np.exp(-sigma * grid.edges_mm / 10.0)
    react = s_edge[:-1] - s_edge[1:]  # total weight ~ 1 for a unit-weight history
    r_tot = float(np.interp(3480.0, carbon.energy_mev, carbon.csda_range_g_per_cm2))
    resid = r_tot - grid.centers_mm / 10.0
    e_carbon = np.where(
        resid > carbon.csda_range_g_per_cm2[0],
        np.interp(resid, carbon.csda_range_g_per_cm2, carbon.energy_mev),
        0.0,
    )
    n_a_over_12 = sum(m * sp.mass_number for sp, m in FRAGMENT_SPECIES) / 12.0
    injected_independent = float(np.sum(react * e_carbon * n_a_over_12))
    assert frag.fragment_energy_injected_mev == pytest.approx(
        injected_independent, rel=1e-6
    )
    # (2) the injected KE is deposited in full (fragments stay inside the grid)
    assert frag.fragment_edep_mev.sum() == pytest.approx(
        frag.fragment_energy_injected_mev, rel=1e-4
    )
    # (3) overall balance, escaped non-negative (target fragments / neutrons escape)
    deposited = float(frag.total_edep_mev.sum())
    assert deposited + frag.escaped_mev == pytest.approx(frag.energy_in_mev, rel=1e-9)
    assert frag.escaped_mev >= 0.0
    assert 0.0 < frag.fragment_edep_mev.sum() < frag.energy_in_mev


@pytest.mark.warp
def test_warp_cpu_fragmentation_matches_reference(warp_module, proton_table) -> None:
    """Deterministic total depth dose (primary + fragments) agrees reference vs
    Warp CPU (the primary and each fragment species use the validated CSDA
    kernels)."""
    grid = DepthDoseGrid(550.0, 275)
    slab = WaterSlab(550.0)
    src = PencilBeamSource(3480.0)
    ref = carbon_fragmentation_depth_dose(
        proton_table,
        slab,
        grid,
        src,
        n_histories=1,
        seed=1,
        path="python",
        straggling=False,
    )
    cpu = carbon_fragmentation_depth_dose(
        proton_table,
        slab,
        grid,
        src,
        n_histories=1,
        seed=1,
        path="warp",
        device="cpu",
        straggling=False,
    )
    rt = ref.total_edep_mev
    ct = cpu.total_edep_mev
    denom = float(rt.sum())
    assert abs(float(ct.sum()) - denom) / denom < 5e-5
    assert np.max(np.abs(ct - rt)) / float(rt.max()) < 5e-3
