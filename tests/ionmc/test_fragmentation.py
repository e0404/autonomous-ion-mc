"""Tests of the bounded carbon-12 fragmentation dose tail (decision 0029).

``carbon_fragmentation_depth_dose`` attenuates the primary carbon by a constant
reaction cross-section and emits forward, same-velocity fragments (proton, alpha,
boron-11) that produce the distal dose tail beyond the Bragg peak. The tests cover
the cross-section, the tail magnitude and reach, primary survival, energy
bookkeeping, and reference/Warp CPU agreement. Realistic 2 mm bins are used (the
tail-to-peak ratio depends on the peak's binning; measured carbon curves have a
straggling-broadened peak at a few-mm resolution).
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.constants import AVOGADRO
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table, scale_ion_stopping_table
from ionmc.fragmentation import (
    CARBON_SIGMA_R_BARN,
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
    """The distal fragment tail is 8-20 % of the Bragg peak (canonical ~15 %)."""
    for distal in (5.0, 10.0, 20.0):
        assert 0.08 <= frag.tail_to_peak(distal) <= 0.20


def test_primary_only_has_no_tail(proton_table) -> None:
    """A primary-only carbon run (no fragmentation) has ~0 dose beyond the peak,
    isolating the tail as a fragmentation effect."""
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
    kpk = int(dd.argmax())
    # a bin ~10 mm beyond the primary peak has negligible primary dose...
    k_beyond = kpk + 5
    assert dd[k_beyond] / dd[kpk] < 1e-3
    # ...while the fragmentation run has a real tail there
    assert frag_tail_ratio(proton_table) > 0.05


def frag_tail_ratio(proton_table) -> float:
    grid = DepthDoseGrid(550.0, 275)
    r = carbon_fragmentation_depth_dose(
        proton_table,
        WaterSlab(550.0),
        grid,
        PencilBeamSource(3480.0),
        n_histories=1,
        seed=1,
        path="python",
        straggling=False,
    )
    return r.tail_to_peak(10.0)


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


def test_energy_conservation(frag) -> None:
    """total deposited + escaped == energy in, with escaped >= 0 (the fragments
    carry ~0.72 of each reaction's energy; the rest escapes)."""
    total = float(frag.total_edep_mev.sum())
    assert total + frag.escaped_mev == pytest.approx(frag.energy_in_mev, rel=1e-9)
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
