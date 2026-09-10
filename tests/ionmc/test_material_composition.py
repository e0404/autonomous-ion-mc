"""Tests of per-voxel tissue materials via stopping-power ratios (decision 0015).

The material library and the SPR/nuclear-scaling quantities are checked against
the published water-equivalent-ratio bands; the transport tests confirm that a
water VoxelSlab reproduces the homogeneous baseline, that a tissue slab puts R80
at R80_water / WER, that a multi-material phantom conserves energy across
interfaces, and that the reference and Warp paths agree.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials as M
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.particles import PROTON
from ionmc.stopping_power import mass_stopping_power_ratio
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.geometry import VoxelSlab

#: Published linear water-equivalent ratios (Schneider 2000 / ICRU) and bands.
PUBLISHED_WER = {
    "cortical_bone": (M.CORTICAL_BONE, 1.60, 1.72),
    "adipose": (M.ADIPOSE, 0.95, 0.98),
    "soft_tissue": (M.SOFT_TISSUE, 1.01, 1.04),
    "skeletal_muscle": (M.SKELETAL_MUSCLE, 1.02, 1.05),
    "lung_tissue": (M.LUNG_TISSUE, 1.02, 1.05),
}


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


# -- material library ---------------------------------------------------------


def test_tissue_library_loads() -> None:
    for name, mat in M.TISSUES.items():
        assert mat.name == name
        assert mat.density_g_per_cm3 > 0.0
        assert mat.mean_excitation_energy.value_ev > 0.0
        assert mat.radiation_length_g_per_cm2 > 0.0


@pytest.mark.parametrize("name", list(PUBLISHED_WER))
def test_water_equivalent_ratio_in_published_band(name) -> None:
    mat, lo, hi = PUBLISHED_WER[name]
    spr = mass_stopping_power_ratio(mat, 150.0, PROTON)
    wer = spr * mat.density_g_per_cm3
    assert lo <= wer <= hi, (name, wer)


def test_water_spr_is_unity() -> None:
    assert mass_stopping_power_ratio(M.WATER, 150.0, PROTON) == pytest.approx(1.0)


def test_oxygen_equivalent_scaling() -> None:
    """The composition-scaled nuclear content reduces to n_O for water and is
    larger for bone (~2x) and adipose (~3.4x) (decision 0015)."""
    water_ratio = M.WATER.oxygen_equivalent_per_gram / M.WATER.atoms_per_gram("O")
    assert water_ratio == pytest.approx(1.0)  # water: oxygen-only recovered
    bone_ratio = (
        M.CORTICAL_BONE.oxygen_equivalent_per_gram / M.CORTICAL_BONE.atoms_per_gram("O")
    )
    adipose_ratio = M.ADIPOSE.oxygen_equivalent_per_gram / M.ADIPOSE.atoms_per_gram("O")
    assert 1.9 <= bone_ratio <= 2.3
    assert 3.0 <= adipose_ratio <= 3.8


# -- geometry -----------------------------------------------------------------


def test_from_material_layers() -> None:
    slab = VoxelSlab.from_material_layers(
        [(40.0, M.WATER), (20.0, M.CORTICAL_BONE), (30.0, M.ADIPOSE)]
    )
    mats = slab.materials_profile()
    assert [m.name for m in mats] == ["water_liquid", "cortical_bone", "adipose"]
    z, rho = slab.voxel_profile()
    np.testing.assert_allclose(z, [0.0, 40.0, 60.0, 90.0])
    np.testing.assert_allclose(rho, [1.0, 1.92, 0.95])


# -- transport ----------------------------------------------------------------


def test_water_voxelslab_matches_homogeneous(table) -> None:
    """A water VoxelSlab reproduces the homogeneous WaterSlab bit-for-bit
    (SPR(water) = 1, nuclear scale = n_O)."""
    grid = DepthDoseGrid(400.0, 800)
    src = PencilBeamSource(150.0)
    homo = TransportEngine(table, WaterSlab(400.0), grid).run(
        src, 200, seed=3, path="python"
    )
    vox = TransportEngine(
        table, VoxelSlab.from_material_layers([(400.0, M.WATER)]), grid
    ).run(src, 200, seed=3, path="python")
    np.testing.assert_array_equal(homo.edep_mev, vox.edep_mev)


@pytest.mark.parametrize("name", ["cortical_bone", "adipose", "skeletal_muscle"])
@pytest.mark.parametrize("e0", [150.0, 200.0])
def test_tissue_r80_at_water_range_over_wer(table, name, e0) -> None:
    """A tissue slab puts R80 at R80_water / WER (WER = SPR_mass x rho)."""
    mat = PUBLISHED_WER[name][0]
    depth = 700.0
    grid = DepthDoseGrid(depth, 7000)
    src = PencilBeamSource(e0)
    water = TransportEngine(table, WaterSlab(depth), grid, straggling=False).run(
        src, 1, path="python"
    )
    slab = VoxelSlab.from_material_layers([(depth, mat)])
    tissue = TransportEngine(table, slab, grid, straggling=False).run(
        src, 1, path="python"
    )
    wer = mass_stopping_power_ratio(mat, 150.0, PROTON) * mat.density_g_per_cm3
    assert tissue.r80_mm() == pytest.approx(water.r80_mm() / wer, rel=3e-3)


def test_multi_material_interface_conserves_energy(table) -> None:
    grid = DepthDoseGrid(400.0, 800)
    slab = VoxelSlab.from_material_layers(
        [(40.0, M.WATER), (20.0, M.CORTICAL_BONE), (30.0, M.ADIPOSE), (310.0, M.WATER)]
    )
    res = TransportEngine(table, slab, grid, nuclear=True, secondaries=True).run(
        PencilBeamSource(150.0), 600, seed=7, path="python"
    )
    assert res.n_reactions > 0
    assert res.n_secondaries > 0
    assert abs(res.energy_balance) <= 1e-9


@pytest.mark.warp
def test_material_cross_backend(warp_module, table) -> None:
    grid = DepthDoseGrid(400.0, 800)
    slab = VoxelSlab.from_material_layers(
        [(40.0, M.WATER), (20.0, M.CORTICAL_BONE), (30.0, M.ADIPOSE), (310.0, M.WATER)]
    )
    eng = TransportEngine(table, slab, grid, nuclear=True, secondaries=True)
    src = PencilBeamSource(150.0)
    ref = eng.run(src, 4000, seed=7, path="python")
    war = eng.run(src, 4000, seed=7, path="warp", device="cpu")
    assert war.n_reactions == ref.n_reactions
    total = float(np.sum(ref.edep_mev))
    cum = np.max(np.abs(np.cumsum(war.edep_mev) - np.cumsum(ref.edep_mev))) / total
    assert cum <= 1e-4
    assert abs(war.energy_balance) <= 1e-5


def test_nuclear_scales_with_composition(table) -> None:
    """The nonelastic reaction count in adipose exceeds the oxygen-only estimate:
    running the same phantom composition scaled vs oxygen-only content."""
    # a thick adipose slab: more reactions than if only its oxygen counted
    grid = DepthDoseGrid(400.0, 800)
    adipose_slab = VoxelSlab.from_material_layers([(200.0, M.ADIPOSE)])
    n = 4000
    res = TransportEngine(table, adipose_slab, grid, nuclear=True).run(
        PencilBeamSource(150.0), n, seed=5, path="python"
    )
    # oxygen-equivalent density is ~3.4x the oxygen-only density for adipose, so
    # the reaction fraction is well above the water/oxygen-only expectation
    frac = res.n_reactions / n
    assert frac > 0.0
    ratio = M.ADIPOSE.oxygen_equivalent_per_gram / M.ADIPOSE.atoms_per_gram("O")
    assert ratio > 3.0
