"""Tests of helium-4 transport via equal-velocity z-squared scaling (decision 0027).

``scale_ion_stopping_table`` builds a helium mass-stopping table from the proton
table (``E_ion = E_p*m_ion/m_p``, ``S_ion = z^2*S_p``); the transport kernels run
unchanged with that table and the ALPHA particle. The tests cover the scaling
construction, the range identity, agreement with the independent analytic Bethe
model, the 600 MeV helium Bragg range, and reference/Warp CPU agreement.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.constants import PROTON_MASS_MEV
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table, scale_ion_stopping_table
from ionmc.materials import WATER
from ionmc.particles import ALPHA
from ionmc.stopping_power import AnalyticStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine
from ionmc.transport.geometry import WaterSlab


@pytest.fixture(scope="module")
def proton_table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.fixture(scope="module")
def helium_table(proton_table):
    return scale_ion_stopping_table(proton_table, ALPHA)


def test_scale_ion_table_z2(proton_table, helium_table) -> None:
    mr = ALPHA.rest_energy_mev / PROTON_MASS_MEV  # m_He/m_p ~ 3.9726
    # energy grid scaled by the mass ratio, stopping by z^2 = 4
    assert np.allclose(helium_table.energy_mev, proton_table.energy_mev * mr)
    assert np.allclose(
        helium_table.stopping_mev_cm2_per_g,
        4.0 * proton_table.stopping_mev_cm2_per_g,
    )
    assert helium_table.provenance["ion"] == "alpha"
    assert helium_table.provenance["ion_charge"] == 2.0


def test_helium_range_identity(proton_table, helium_table) -> None:
    """R_He(E) = (m_He/(z^2 m_p)) * R_p(E*m_p/m_He) exactly by construction."""
    mr = ALPHA.rest_energy_mev / PROTON_MASS_MEV
    prefactor = ALPHA.rest_energy_mev / (4.0 * PROTON_MASS_MEV)
    for e_he in (200.0, 400.0, 600.0):
        r_he = float(
            np.interp(e_he, helium_table.energy_mev, helium_table.csda_range_g_per_cm2)
        )
        r_p = float(
            np.interp(
                e_he / mr,
                proton_table.energy_mev,
                proton_table.csda_range_g_per_cm2,
            )
        )
        assert r_he == pytest.approx(prefactor * r_p, rel=1e-6)
    # 600 MeV (150 MeV/u) helium ranges to ~15.86 g/cm^2 (== 150 MeV proton depth)
    r600 = float(
        np.interp(600.0, helium_table.energy_mev, helium_table.csda_range_g_per_cm2)
    )
    assert r600 == pytest.approx(15.86, abs=0.1)


@pytest.mark.parametrize("e_per_u", [10.0, 50.0, 100.0, 150.0, 250.0])
def test_helium_stopping_vs_bethe(helium_table, e_per_u: float) -> None:
    """The scaled-PSTAR helium table matches the independent analytic Bethe model
    (which reproduces PSTAR to <0.1% for protons) to ~1.5% (decision 0027)."""
    e_he = e_per_u * 4.0
    s_table = float(
        np.interp(e_he, helium_table.energy_mev, helium_table.stopping_mev_cm2_per_g)
    )
    bethe = AnalyticStoppingPower(WATER, ALPHA, path="numpy")
    s_bethe = float(bethe.mass_stopping_power(e_he)[0])
    assert s_table == pytest.approx(s_bethe, rel=0.015)


def test_helium_bragg_range(helium_table) -> None:
    """A 600 MeV (150 MeV/u) helium beam stops at ~158 mm in water, conserving
    energy (the Bragg-curve target, == a 150 MeV proton range)."""
    grid = DepthDoseGrid(200.0, 400)
    eng = TransportEngine(
        helium_table,
        WaterSlab(200.0),
        grid,
        particle=ALPHA,
        straggling=False,
        nuclear=False,
    )
    r = eng.run(PencilBeamSource(600.0), n_histories=1, seed=1, path="python")
    dd = r.edep_mev
    centers = grid.centers_mm
    kpk = int(dd.argmax())
    assert 155.0 <= centers[kpk] <= 161.0  # Bragg peak near 158 mm
    assert abs(r.energy_balance) < 1e-9  # contained beam conserves energy


@pytest.mark.warp
def test_warp_cpu_helium_matches_reference(warp_module, helium_table) -> None:
    """Deterministic helium depth dose: Warp CPU matches the reference (the
    species-agnostic kernels transport helium consistently)."""
    grid = DepthDoseGrid(200.0, 200)
    eng = TransportEngine(
        helium_table,
        WaterSlab(200.0),
        grid,
        particle=ALPHA,
        straggling=False,
        nuclear=False,
    )
    src = PencilBeamSource(600.0)
    ref = eng.run(src, n_histories=1, seed=2, path="python")
    cpu = eng.run(src, n_histories=1, seed=2, path="warp", device="cpu")
    denom = float(ref.edep_mev.sum())
    assert abs(float(cpu.edep_mev.sum()) - denom) / denom < 1e-5
    assert (
        np.max(np.abs(cpu.edep_mev - ref.edep_mev)) / float(ref.edep_mev.max()) < 5e-3
    )
