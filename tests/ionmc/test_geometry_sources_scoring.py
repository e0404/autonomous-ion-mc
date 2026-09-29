import numpy as np
import pytest

from ionmc.geometry import homogeneous_box, slab_phantom
from ionmc.scoring import ScoringGrid, Tally
from ionmc.sources import PencilBeam, orthonormal_frame


def test_homogeneous_box_conventions():
    geo = homogeneous_box((120.0, 120.0, 300.0), 2.0, "water")
    assert geo.shape == (60, 60, 150)
    assert geo.origin_mm == (-60.0, -60.0, 0.0)
    assert geo.edges_mm(2)[0] == 0.0 and geo.edges_mm(2)[-1] == 300.0
    assert geo.voxel_volume_cm3 == pytest.approx(0.008)
    assert geo.mass_g().sum() == pytest.approx(120 * 120 * 300 * 1e-3)
    with pytest.raises(ValueError):
        homogeneous_box((121.0, 120.0, 300.0), 2.0)


def test_slab_phantom_materials_and_density():
    geo = slab_phantom(
        (20.0, 20.0, 100.0),
        2.0,
        [(20.0, 40.0, "bone_cortical", None), (60.0, 80.0, "lung_inflated", 0.3)],
    )
    names = [m.name for m in geo.materials]
    assert names == ["water", "bone_cortical", "lung_inflated"]
    zc = geo.centers_mm(2)
    assert np.all(geo.material_index[0, 0, (zc > 20) & (zc < 40)] == 1)
    assert np.all(geo.density_g_cm3[0, 0, (zc > 60) & (zc < 80)] == np.float32(0.3))
    assert np.all(geo.density_g_cm3[0, 0, zc < 20] == np.float32(1.0))


def test_pencil_beam_sampling_statistics():
    beam = PencilBeam(
        "proton",
        150.0,
        (1.0, 2.0, -5.0),
        (0.0, 0.0, 1.0),
        sigma_mm=3.0,
        sigma_energy_fraction=0.01,
        sigma_angle_rad=0.002,
    )
    rng = np.random.default_rng(0)
    s = beam.sample(20000, rng)
    assert s["position_mm"][:, 2].std() == 0.0
    assert s["position_mm"][:, 0].std() == pytest.approx(3.0, rel=0.05)
    assert s["position_mm"][:, 0].mean() == pytest.approx(1.0, abs=0.1)
    assert np.allclose(np.linalg.norm(s["direction"], axis=1), 1.0)
    assert s["direction"][:, 0].std() == pytest.approx(0.002, rel=0.05)
    assert s["energy_mev_per_u"].std() / 150.0 == pytest.approx(0.01, rel=0.05)
    u, v, w = orthonormal_frame((0.3, -0.4, 0.5))
    assert abs(np.dot(u, w)) < 1e-12 and abs(np.dot(u, v)) < 1e-12
    oblique = PencilBeam("c12", 290.0, direction=(1.0, 1.0, 0.0))
    assert oblique.direction == pytest.approx((2**-0.5, 2**-0.5, 0.0))
    assert oblique.describe()["energy_total_mev"] == pytest.approx(3480.0)
    with pytest.raises(ValueError):
        PencilBeam("proton", -1.0)


def test_scoring_grid_mass_matches_geometry_for_coarse_and_shifted_grids():
    geo = homogeneous_box((20.0, 20.0, 40.0), 2.0, "water")
    fine = ScoringGrid.from_geometry(geo)
    assert fine.voxel_mass_g(geo).sum() == pytest.approx(geo.mass_g().sum())
    coarse = ScoringGrid.coarse(geo, 4.0)
    assert coarse.shape == (5, 5, 10)
    assert coarse.voxel_mass_g(geo).sum() == pytest.approx(geo.mass_g().sum())
    assert coarse.voxel_mass_g(geo)[0, 0, 0] == pytest.approx(0.064)
    shifted = fine.shifted((1.0, 0.0, 0.0))  # half a voxel: last column half outside
    m = shifted.voxel_mass_g(geo)
    assert m[-1, 0, 0] == pytest.approx(0.004)
    assert m[0, 0, 0] == pytest.approx(0.008)


def test_tally_batch_statistics():
    t = Tally("e", (2,), 4)
    for value in (1.0, 2.0, 3.0, 4.0):
        t.current[:] = value * 10
        t.close_batch(10)
    assert np.allclose(t.mean(), 2.5)
    assert np.allclose(t.standard_error(), np.sqrt(5.0 / 3.0 / 4.0))
    t1 = Tally("e", (2,), 1)
    t1.current[:] = 1.0
    t1.close_batch(1)
    assert np.all(np.isnan(t1.standard_error()))
