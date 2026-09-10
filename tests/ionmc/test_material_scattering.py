"""Tests of non-water materials on the 3-D scattering path (decision 0017).

The scattering path now uses the per-voxel water-equivalent density for energy
loss and the per-voxel physical density and material radiation length for the
multiple scattering. The tests confirm the material lateral spread against the
material-aware Fermi-Eyges oracle (homogeneous bone and a bone interface), the
water-equivalent range, energy conservation, and cross-backend agreement.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials as M
from ionmc import particles
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.physics import fermi_eyges as fe
from ionmc.stopping_power import mass_stopping_power_ratio
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine
from ionmc.transport.depth_dose import DepthLateralGrid
from ionmc.transport.geometry import VoxelSlab


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


@pytest.fixture(scope="module")
def model(table):
    return TabulatedStoppingPower(table, M.WATER, particles.PROTON, "numpy")


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(
        depth_mm=depth, n_depth=int(depth * 2), half_width_mm=25.0, n_lateral=600
    )


def _spr(mat) -> float:
    return mass_stopping_power_ratio(mat, 150.0, particles.PROTON)


def _material_arrays(slab: VoxelSlab):
    zb, rho = slab.voxel_profile()
    mats = slab.materials_profile()
    we = np.array([_spr(m) * rho[i] for i, m in enumerate(mats)])
    phys = np.asarray(rho, dtype=np.float64)
    rl = np.array([m.radiation_length_g_per_cm2 for m in mats])
    return zb, we, phys, rl


def test_scattering_requires_radiation_length(table) -> None:
    """A material without a radiation length is rejected on the scattering path."""
    no_x0 = M.CORTICAL_BONE.__class__(
        name="bone_no_x0",
        density_g_per_cm3=1.92,
        mass_fractions=dict(M.CORTICAL_BONE.mass_fractions),
        mean_excitation_energy=M.CORTICAL_BONE.mean_excitation_energy,
    )  # radiation_length defaults to 0.0
    eng = TransportEngine(
        table,
        VoxelSlab.from_material_layers([(100.0, no_x0)]),
        DepthDoseGrid(200.0, 10),
    )
    with pytest.raises(ValueError, match="radiation length"):
        eng.run_scattering(PencilBeamSource(150.0), _lat(200.0), 1, path="python")


@pytest.mark.warp
def test_bone_sigma_x_matches_material_oracle(warp_module, table, model) -> None:
    """A homogeneous bone slab reproduces the material-aware Fermi-Eyges sigma_x
    (physical density and bone X0) and puts R80 at R_water / WER."""
    bone = VoxelSlab.from_material_layers([(120.0, M.CORTICAL_BONE)])
    zb, we, phys, rl = _material_arrays(bone)
    r = TransportEngine(table, bone, DepthDoseGrid(120.0, 10)).run_scattering(
        PencilBeamSource(150.0), _lat(120.0), 40000, seed=11, path="warp", device="cpu"
    )
    wer = _spr(M.CORTICAL_BONE) * M.CORTICAL_BONE.density_g_per_cm3
    r0 = float(model.csda_range(150.0)[0]) * 10.0
    for frac in (0.4, 0.6):
        z = frac * r0 / wer
        mc = r.sigma_x_at_depth(z)
        oracle = float(
            fe.lateral_sigma_x_material_mm(
                model, 150.0, np.array([z]), zb, we, phys, rl
            )[0]
        )
        assert abs(mc / oracle - 1.0) <= 0.03, (z, mc, oracle)


@pytest.mark.warp
def test_material_interface_sigma_x(warp_module, table, model) -> None:
    """A water/bone/water phantom reproduces the material-aware Fermi-Eyges
    sigma_x across the interface."""
    layered = VoxelSlab.from_material_layers(
        [(40.0, M.WATER), (20.0, M.CORTICAL_BONE), (190.0, M.WATER)]
    )
    zb, we, phys, rl = _material_arrays(layered)
    r = TransportEngine(table, layered, DepthDoseGrid(250.0, 10)).run_scattering(
        PencilBeamSource(150.0), _lat(250.0), 40000, seed=11, path="warp", device="cpu"
    )
    for z in (40.0, 60.0, 90.0):
        mc = r.sigma_x_at_depth(z)
        oracle = float(
            fe.lateral_sigma_x_material_mm(
                model, 150.0, np.array([z]), zb, we, phys, rl
            )[0]
        )
        assert abs(mc / oracle - 1.0) <= 0.03, (z, mc, oracle)


def test_material_interface_conserves_energy(table) -> None:
    layered = VoxelSlab.from_material_layers(
        [(40.0, M.WATER), (20.0, M.CORTICAL_BONE), (190.0, M.ADIPOSE)]
    )
    r = TransportEngine(table, layered, DepthDoseGrid(250.0, 10)).run_scattering(
        PencilBeamSource(150.0), _lat(250.0), 400, seed=7, path="python"
    )
    assert abs(r.energy_balance) <= 1e-9


@pytest.mark.warp
def test_material_scattering_cross_backend(warp_module, table) -> None:
    layered = VoxelSlab.from_material_layers(
        [(40.0, M.WATER), (20.0, M.CORTICAL_BONE), (190.0, M.WATER)]
    )
    eng = TransportEngine(table, layered, DepthDoseGrid(250.0, 10))
    src = PencilBeamSource(150.0)
    lat = _lat(250.0)
    ref = eng.run_scattering(src, lat, 4000, seed=7, path="python")
    war = eng.run_scattering(src, lat, 4000, seed=7, path="warp", device="cpu")
    for z in (40.0, 60.0, 90.0):
        assert abs(ref.sigma_x_at_depth(z) - war.sigma_x_at_depth(z)) <= 0.05
    dd_ref = ref.depth_dose_mev
    cum = np.max(np.abs(np.cumsum(war.depth_dose_mev) - np.cumsum(dd_ref))) / np.sum(
        dd_ref
    )
    assert cum <= 5e-4
