"""Tests of density-heterogeneous 3-D multiple-scattering transport (0016).

A ``VoxelSlab`` of water at varying density drives the scattering path's local
density lookup. The tests confirm homogeneous equivalence (bit-exact), the
lateral spread scaling with density and across a layered interface (validated
against the Fermi-Eyges oracle, including its piecewise-density variant), energy
conservation, the non-water rejection, and cross-backend agreement.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials, particles
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import fermi_eyges as fe
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.depth_dose import DepthLateralGrid
from ionmc.transport.geometry import VoxelSlab

X0 = materials.WATER.radiation_length_g_per_cm2


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.fixture(scope="module")
def model(table):
    return TabulatedStoppingPower(table, materials.WATER, particles.PROTON, "numpy")


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(
        depth_mm=depth, n_depth=int(depth * 2), half_width_mm=25.0, n_lateral=600
    )


def test_uniform_voxels_match_homogeneous_scattering(table) -> None:
    """A uniform water VoxelSlab reproduces the homogeneous WaterSlab scattering
    bit-for-bit (the merge collapses it to one voxel)."""
    src = PencilBeamSource(150.0)
    lat = _lat(250.0)
    homo = TransportEngine(table, WaterSlab(250.0), DepthDoseGrid(250.0, 10))
    vox = TransportEngine(
        table, VoxelSlab.uniform(250.0, 1.0, 100), DepthDoseGrid(250.0, 10)
    )
    a = homo.run_scattering(src, lat, 400, seed=3, path="python")
    b = vox.run_scattering(src, lat, 400, seed=3, path="python")
    np.testing.assert_array_equal(a.edep_zx_mev, b.edep_zx_mev)


def test_scattering_supports_non_water_material(table) -> None:
    """The 3-D scattering path now supports non-water materials (decision 0017):
    a bone slab runs and conserves energy rather than raising."""
    eng = TransportEngine(
        table,
        VoxelSlab.from_material_layers([(100.0, materials.CORTICAL_BONE)]),
        DepthDoseGrid(250.0, 10),
    )
    res = eng.run_scattering(
        PencilBeamSource(150.0), _lat(250.0), 100, seed=1, path="python"
    )
    assert abs(res.energy_balance) <= 1e-9


@pytest.mark.warp
def test_uniform_density_sigma_x_matches_fermi_eyges(warp_module, table, model) -> None:
    """A uniform slab of density rho reproduces the Fermi-Eyges sigma_x at that
    density (lateral spread scales with density)."""
    rho = 1.2
    lat = _lat(140.0)
    r = TransportEngine(
        table, VoxelSlab.uniform(140.0, rho, 1), DepthDoseGrid(140.0, 10)
    ).run_scattering(
        PencilBeamSource(150.0), lat, 40000, seed=11, path="warp", device="cpu"
    )
    r0 = float(model.csda_range(150.0)[0]) * 10.0
    for frac in (0.5, 0.8):
        z = frac * r0 / rho
        mc = r.sigma_x_at_depth(z)
        oracle = float(
            fe.lateral_sigma_x_mm(
                model, 150.0, np.array([z]), X0, density_g_per_cm3=rho
            )[0]
        )
        assert abs(mc / oracle - 1.0) <= 0.03, (z, mc, oracle)


@pytest.mark.warp
def test_layered_sigma_x_matches_heterogeneous_fermi_eyges(
    warp_module, table, model
) -> None:
    """A water/dense/water phantom reproduces the piecewise-density Fermi-Eyges
    sigma_x (the interface produces the correct kink)."""
    layered = VoxelSlab.from_layers([(60.0, 1.0), (30.0, 1.4), (160.0, 1.0)])
    zb, dens = layered.voxel_profile()
    lat = _lat(250.0)
    r = TransportEngine(table, layered, DepthDoseGrid(250.0, 10)).run_scattering(
        PencilBeamSource(150.0), lat, 40000, seed=11, path="warp", device="cpu"
    )
    for z in (60.0, 90.0, 120.0):
        mc = r.sigma_x_at_depth(z)
        oracle = float(
            fe.lateral_sigma_x_heterogeneous_mm(
                model, 150.0, np.array([z]), X0, zb, dens
            )[0]
        )
        assert abs(mc / oracle - 1.0) <= 0.03, (z, mc, oracle)


def test_layered_scattering_conserves_energy(table) -> None:
    layered = VoxelSlab.from_layers([(60.0, 1.0), (30.0, 1.4), (160.0, 1.0)])
    r = TransportEngine(table, layered, DepthDoseGrid(250.0, 10)).run_scattering(
        PencilBeamSource(150.0), _lat(250.0), 400, seed=7, path="python"
    )
    assert abs(r.energy_balance) <= 1e-9


@pytest.mark.warp
def test_layered_scattering_cross_backend(warp_module, table) -> None:
    layered = VoxelSlab.from_layers([(60.0, 1.0), (30.0, 1.4), (160.0, 1.0)])
    eng = TransportEngine(table, layered, DepthDoseGrid(250.0, 10))
    src = PencilBeamSource(150.0)
    lat = _lat(250.0)
    ref = eng.run_scattering(src, lat, 4000, seed=7, path="python")
    war = eng.run_scattering(src, lat, 4000, seed=7, path="warp", device="cpu")
    # sigma_x agrees tightly; the depth dose agrees within the 3-D float32 budget
    for z in (60.0, 90.0, 120.0):
        assert abs(ref.sigma_x_at_depth(z) - war.sigma_x_at_depth(z)) <= 0.05
    dd_ref = ref.depth_dose_mev
    dd_war = war.depth_dose_mev
    total = float(np.sum(dd_ref))
    cum = np.max(np.abs(np.cumsum(dd_war) - np.cumsum(dd_ref))) / total
    assert cum <= 5e-4
