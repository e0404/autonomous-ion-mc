"""Tests of 1-D voxelized density-heterogeneous transport (decision 0014).

The geometry (``VoxelSlab``) is checked directly (layer boundaries, water-
equivalent thickness, validation); the transport tests confirm that a uniform
voxel slab reproduces the homogeneous ``WaterSlab`` bit-for-bit, that the R80
depth scales as ``1/rho`` (water-equivalent thickness), that a dense layer
shifts the peak by its extra water-equivalent thickness, that the energy budget
still closes with nuclear and secondary transport in a heterogeneous phantom,
and that the reference and Warp paths agree across a density interface.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.geometry import VoxelSlab

ENERGIES = (150.0, 200.0)


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


# -- geometry -----------------------------------------------------------------


def test_from_layers_boundaries_and_wet() -> None:
    slab = VoxelSlab.from_layers([(50.0, 1.0), (20.0, 1.85), (30.0, 1.0)])
    z, rho = slab.voxel_profile()
    np.testing.assert_allclose(z, [0.0, 50.0, 70.0, 100.0])
    np.testing.assert_allclose(rho, [1.0, 1.85, 1.0])
    assert slab.depth_mm == 100.0
    # WET to 100 mm = 50*1 + 20*1.85 + 30*1 = 117 mm
    assert slab.water_equivalent_depth_mm(100.0) == pytest.approx(117.0)
    # WET partway into the dense layer (60 mm) = 50 + 10*1.85 = 68.5
    assert slab.water_equivalent_depth_mm(60.0) == pytest.approx(68.5)


def test_voxelslab_validation() -> None:
    with pytest.raises(ValueError):
        VoxelSlab(np.array([0.0, 1.0]), np.array([1.0, 2.0]))  # mismatched lengths
    with pytest.raises(ValueError):
        VoxelSlab(np.array([0.0, 2.0, 1.0]), np.array([1.0, 1.0]))  # not increasing
    with pytest.raises(ValueError):
        VoxelSlab(np.array([1.0, 2.0]), np.array([1.0]))  # first boundary not 0
    with pytest.raises(ValueError):
        VoxelSlab(np.array([0.0, 1.0]), np.array([-1.0]))  # non-positive density


# -- transport ----------------------------------------------------------------


def test_uniform_voxels_match_homogeneous(table) -> None:
    """A uniform voxel slab (density 1.0) reproduces WaterSlab bit-for-bit: the
    merge collapses equal-density voxels so the step limiting is a no-op."""
    grid = DepthDoseGrid(400.0, 800)
    src = PencilBeamSource(150.0)
    homo = TransportEngine(table, WaterSlab(400.0, materials.WATER), grid).run(
        src, 200, seed=3, path="python"
    )
    vox = TransportEngine(table, VoxelSlab.uniform(400.0, 1.0, 200), grid).run(
        src, 200, seed=3, path="python"
    )
    np.testing.assert_array_equal(homo.edep_mev, vox.edep_mev)


@pytest.mark.parametrize("rho", (0.5, 1.2))
@pytest.mark.parametrize("e0", ENERGIES)
def test_r80_scales_with_inverse_density(table, e0, rho) -> None:
    """Water-equivalent thickness: a uniform slab of density rho puts R80 at
    R80_water / rho (CSDA range scales exactly as 1/rho)."""
    # deep enough that the 1/rho-scaled range still stops inside the phantom
    depth = 700.0
    grid = DepthDoseGrid(depth, 7000)  # 0.1 mm bins to resolve the shift
    src = PencilBeamSource(e0)
    water = TransportEngine(table, WaterSlab(depth), grid, straggling=False).run(
        src, 1, path="python"
    )
    dense = TransportEngine(
        table, VoxelSlab.uniform(depth, rho, 1), grid, straggling=False
    ).run(src, 1, path="python")
    assert dense.r80_mm() == pytest.approx(water.r80_mm() / rho, rel=3e-3)


@pytest.mark.parametrize("e0", ENERGIES)
def test_dense_layer_shifts_peak_by_wet(table, e0) -> None:
    """A dense layer shifts the peak proximally by its extra water-equivalent
    thickness (thickness x (rho - 1))."""
    grid = DepthDoseGrid(400.0, 4000)
    src = PencilBeamSource(e0)
    water = TransportEngine(table, WaterSlab(400.0), grid, straggling=False).run(
        src, 1, path="python"
    )
    thickness, rho = 20.0, 1.85
    layered = VoxelSlab.from_layers(
        [(40.0, 1.0), (thickness, rho), (360.0 - thickness, 1.0)]
    )
    lay = TransportEngine(table, layered, grid, straggling=False).run(
        src, 1, path="python"
    )
    shift = water.r80_mm() - lay.r80_mm()
    assert shift == pytest.approx(thickness * (rho - 1.0), abs=0.5)


def test_energy_conserved_in_heterogeneous_phantom(table) -> None:
    """Energy budget closes with nuclear + secondaries across a density step."""
    grid = DepthDoseGrid(400.0, 800)
    layered = VoxelSlab.from_layers([(50.0, 1.0), (20.0, 1.85), (330.0, 1.0)])
    res = TransportEngine(table, layered, grid, nuclear=True, secondaries=True).run(
        PencilBeamSource(150.0), 600, seed=7, path="python"
    )
    assert res.n_reactions > 0
    assert res.n_secondaries > 0
    assert abs(res.energy_balance) <= 1e-9


def test_scattering_rejects_heterogeneous_slab(table) -> None:
    """The homogeneous-only 3-D scattering path rejects a heterogeneous VoxelSlab
    rather than silently using the front-voxel density (decision 0014)."""
    from ionmc.transport.depth_dose import DepthLateralGrid

    layered = VoxelSlab.from_layers([(40.0, 1.0), (20.0, 1.85), (340.0, 1.0)])
    eng = TransportEngine(table, layered, DepthDoseGrid(400.0, 800))
    lat = DepthLateralGrid(
        depth_mm=400.0, n_depth=800, half_width_mm=25.0, n_lateral=100
    )
    with pytest.raises(NotImplementedError, match="homogeneous"):
        eng.run_scattering(PencilBeamSource(150.0), lat, 1, path="python")


@pytest.mark.warp
def test_heterogeneous_cross_backend(warp_module, table) -> None:
    """Reference and Warp CPU agree across a density interface (with nuclear and
    secondary transport on)."""
    grid = DepthDoseGrid(400.0, 800)
    layered = VoxelSlab.from_layers([(50.0, 1.0), (20.0, 1.85), (330.0, 1.0)])
    eng = TransportEngine(table, layered, grid, nuclear=True, secondaries=True)
    src = PencilBeamSource(150.0)
    ref = eng.run(src, 4000, seed=7, path="python")
    war = eng.run(src, 4000, seed=7, path="warp", device="cpu")
    assert war.n_reactions == ref.n_reactions
    total = float(np.sum(ref.edep_mev))
    cum = np.max(np.abs(np.cumsum(war.edep_mev) - np.cumsum(ref.edep_mev))) / total
    assert cum <= 1e-4
    assert abs(war.energy_balance) <= 1e-5
