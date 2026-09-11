"""Tests of carbon-12 / oxygen-16 primary transport via z-squared scaling
(decision 0028).

`scale_ion_stopping_table` is species-generic, so carbon and oxygen are
transported exactly as helium (decision 0027): build the scaled table, run with
the corresponding particle. The tests confirm the scaling accuracy against the
independent Bethe model (Z-dependent), the range identity, and the carbon Bragg
range; nuclear fragmentation (the distal fragment tail) is deferred.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.constants import PROTON_MASS_MEV
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table, scale_ion_stopping_table
from ionmc.materials import WATER
from ionmc.particles import CARBON_12, OXYGEN_16
from ionmc.stopping_power import AnalyticStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine
from ionmc.transport.geometry import WaterSlab


@pytest.fixture(scope="module")
def proton_table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.mark.parametrize(("particle", "tol"), [(CARBON_12, 0.02), (OXYGEN_16, 0.03)])
def test_carbon_oxygen_stopping_vs_bethe(proton_table, particle, tol: float) -> None:
    """The scaled table matches the independent Bethe model within the ion-specific
    tolerance (the z^2 scaling omits Barkas/Bloch terms that grow with z)."""
    tab = scale_ion_stopping_table(proton_table, particle)
    bethe = AnalyticStoppingPower(WATER, particle, path="numpy")
    # top energy stays inside the scaled table's domain: the proton table ceils at
    # 400 MeV, so the ion table reaches only 400*(m_p/m_ion)^-1 ~ 397 MeV/u for
    # carbon; 380 MeV/u avoids np.interp silently clamping (decision 0028).
    for e_per_u in (10.0, 50.0, 100.0, 150.0, 250.0, 380.0):
        e = e_per_u * particle.mass_number
        s_tab = float(np.interp(e, tab.energy_mev, tab.stopping_mev_cm2_per_g))
        s_bethe = float(bethe.mass_stopping_power(e)[0])
        assert s_tab == pytest.approx(s_bethe, rel=tol)


@pytest.mark.parametrize("particle", [CARBON_12, OXYGEN_16])
def test_ion_range_identity(proton_table, particle) -> None:
    """R_ion(E) = (m_ion/(z^2 m_p)) * R_p(E*m_p/m_ion) exactly by construction."""
    tab = scale_ion_stopping_table(proton_table, particle)
    mr = particle.rest_energy_mev / PROTON_MASS_MEV
    pref = particle.rest_energy_mev / (particle.charge**2 * PROTON_MASS_MEV)
    for e_per_u in (100.0, 290.0):
        e = e_per_u * particle.mass_number
        r_ion = float(np.interp(e, tab.energy_mev, tab.csda_range_g_per_cm2))
        r_p = float(
            np.interp(
                e / mr, proton_table.energy_mev, proton_table.csda_range_g_per_cm2
            )
        )
        assert r_ion == pytest.approx(pref * r_p, rel=1e-6)


def test_carbon_bragg_range(proton_table) -> None:
    """A 290 MeV/u carbon beam (3480 MeV) stops at ~163 mm in water (its clinical
    range), conserving energy. This is the PRIMARY Bragg curve; the distal fragment
    tail (nuclear fragmentation) is deferred (decision 0028)."""
    tab = scale_ion_stopping_table(proton_table, CARBON_12)
    grid = DepthDoseGrid(250.0, 500)
    eng = TransportEngine(
        tab,
        WaterSlab(250.0),
        grid,
        particle=CARBON_12,
        straggling=False,
        nuclear=False,
    )
    r = eng.run(PencilBeamSource(3480.0), n_histories=1, seed=1, path="python")
    kpk = int(r.edep_mev.argmax())
    assert 159.0 <= grid.centers_mm[kpk] <= 167.0  # ~163 mm
    assert abs(r.energy_balance) < 1e-9


@pytest.mark.warp
def test_warp_cpu_carbon_matches_reference(warp_module, proton_table) -> None:
    """Deterministic carbon depth dose: Warp CPU matches the reference."""
    tab = scale_ion_stopping_table(proton_table, CARBON_12)
    grid = DepthDoseGrid(250.0, 250)
    eng = TransportEngine(
        tab,
        WaterSlab(250.0),
        grid,
        particle=CARBON_12,
        straggling=False,
        nuclear=False,
    )
    src = PencilBeamSource(3480.0)
    ref = eng.run(src, n_histories=1, seed=2, path="python")
    cpu = eng.run(src, n_histories=1, seed=2, path="warp", device="cpu")
    denom = float(ref.edep_mev.sum())
    # heavier-ion float32 budget: carbon (z^2=36, 3480 MeV, ~500 steps) accumulates
    # ~1.4e-5 relative float32 rounding in the total, more than the proton/helium
    # ~6e-6 -- benign summation over the larger step count and total energy.
    assert abs(float(cpu.edep_mev.sum()) - denom) / denom < 5e-5
    assert (
        np.max(np.abs(cpu.edep_mev - ref.edep_mev)) / float(ref.edep_mev.max()) < 5e-3
    )
