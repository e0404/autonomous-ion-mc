"""Tests of arbitrary beam incidence via beam-frame transport (decision 0018).

The scattering path transports in a canonical beam frame (origin at the entry
point, +z' along the beam) and looks up the voxel by the material coordinate
``u = normal . position``. These tests confirm:

* the axis-aligned default (+z beam, +z slab normal) is unchanged;
* a rigid rotation of the whole scene (beam + slab normal) leaves the
  frame-invariant observables unchanged -- bit-exact for a deterministic
  straight ray (scattering off) and within statistics with scattering on;
* oblique incidence on a lab-fixed slab lands the peak at the expected
  water-equivalent depth (geometric R80 scales as R80_normal / cos(theta));
* the pure depth-dose path rejects a non-+z beam.

All checks use the reference Python path (no GPU); cross-backend parity is a
host-runner validation gate (``validation/v3_arbitrary_incidence.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import (
    DepthDoseGrid,
    PencilBeamSource,
    TransportEngine,
    WaterSlab,
)
from ionmc.transport.depth_dose import DepthLateralGrid
from ionmc.transport.engine import _beam_frame, _transverse_frame
from ionmc.transport.geometry import VoxelSlab


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(
        depth_mm=depth, n_depth=int(depth * 2), half_width_mm=25.0, n_lateral=600
    )


def _rot_x(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def _rot_y(deg: float) -> np.ndarray:
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _r80_mm(res) -> float:
    """Distal 80 %-of-max depth of the beam-frame integral depth dose [mm]."""
    dd = res.depth_dose_mev
    centers = res.grid.depth_centers_mm
    peak = float(np.max(dd))
    i_peak = int(np.argmax(dd))
    level = 0.8 * peak
    for i in range(i_peak, len(dd) - 1):
        if dd[i] >= level >= dd[i + 1]:
            f = (dd[i] - level) / (dd[i] - dd[i + 1])
            return float(centers[i] + f * (centers[i + 1] - centers[i]))
    return float(centers[i_peak])


def test_beam_frame_material_coordinate_identity() -> None:
    """The beam-frame coefficients satisfy the defining identity
    ``u0 + m.p_beam == normal.(p0 + R.p_beam)`` for arbitrary direction, normal,
    entry point and beam-frame position (falsifies a wrong m0/m1/m2/u0, which the
    physics gates -- all with normal==direction or on-axis rays -- cannot)."""
    rng = np.random.default_rng(20260913)
    for _ in range(200):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        p0 = rng.normal(size=3) * 40.0
        e1x, e1y, e1z, e2x, e2y, e2z = _transverse_frame(*d)
        rot = np.array([[e1x, e2x, d[0]], [e1y, e2y, d[1]], [e1z, e2z, d[2]]])
        m0, m1, m2, u0 = _beam_frame(d, n, p0)
        for _ in range(5):
            pb = rng.normal(size=3) * 30.0
            lhs = u0 + m0 * pb[0] + m1 * pb[1] + m2 * pb[2]
            rhs = float(n @ (p0 + rot @ pb))
            assert abs(lhs - rhs) < 1e-9


def test_oblique_heterogeneous_interior_crossing(table) -> None:
    """A beam tilted by theta through a lab-fixed +z *heterogeneous* slab crosses
    the same water-equivalent path -- layer by layer -- as a normal beam through
    the same layers each stretched by 1/cos(theta). Deposited energy and the
    beam-depth profile agree to round-off (deterministic, scattering off). This
    exercises interior oblique voxel-boundary crossings, not just the final escape
    (decision 0018)."""
    theta = 25.0
    cos_t = float(np.cos(np.radians(theta)))
    layers = [(20.0, 1.0), (15.0, 1.7), (25.0, 1.0)]  # interior density boundaries
    stretched = [(t / cos_t, rho) for t, rho in layers]
    grid = DepthDoseGrid(300.0, 10)
    lat = _lat(300.0)
    kw = dict(straggling=False, scattering=False)
    normal_eff = TransportEngine(
        table, VoxelSlab.from_layers(stretched), grid, **kw
    ).run_scattering(PencilBeamSource(150.0), lat, 1, seed=3, path="python")
    d = _rot_y(theta) @ np.array([0.0, 0.0, 1.0])
    oblique = TransportEngine(
        table, VoxelSlab.from_layers(layers), grid, **kw
    ).run_scattering(
        PencilBeamSource(150.0, direction=tuple(d)), lat, 1, seed=3, path="python"
    )
    assert 0.0 < oblique.energy_deposited_mev < 150.0  # escapes the back face
    e_ref = normal_eff.energy_deposited_mev
    assert abs(oblique.energy_deposited_mev - e_ref) / e_ref < 1e-6
    denom = np.sum(normal_eff.depth_dose_mev)
    cum = np.max(
        np.abs(np.cumsum(oblique.depth_dose_mev) - np.cumsum(normal_eff.depth_dose_mev))
    )
    assert cum / denom < 1e-9


def test_source_direction_is_normalised() -> None:
    src = PencilBeamSource(150.0, direction=(0.0, 0.0, 3.0))
    np.testing.assert_allclose(src.direction_hat, [0.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="non-zero"):
        PencilBeamSource(150.0, direction=(0.0, 0.0, 0.0))


def test_slab_normal_is_normalised() -> None:
    slab = WaterSlab(200.0, normal=(0.0, 0.0, 5.0))
    np.testing.assert_allclose(slab.normal_hat, [0.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="non-zero"):
        VoxelSlab.uniform(200.0, 1.0, 1, normal=(0.0, 0.0, 0.0))


def test_depth_dose_path_rejects_oblique_beam(table) -> None:
    """The longitudinal depth-dose ``run`` accepts only +z beams (decision 0018)."""
    eng = TransportEngine(table, WaterSlab(250.0), DepthDoseGrid(250.0, 250))
    with pytest.raises(ValueError, match="arbitrary beam incidence"):
        eng.run(PencilBeamSource(150.0, direction=(0.2, 0.0, 1.0)), 1, path="python")


def test_scattering_toggle_removes_lateral_spread(table) -> None:
    """With scattering off the pencil beam stays on axis (sigma_x ~ 0); with it
    on the lateral spread is finite."""
    src = PencilBeamSource(150.0)
    lat = _lat(250.0)
    off = TransportEngine(
        table, WaterSlab(250.0), DepthDoseGrid(250.0, 10), scattering=False
    ).run_scattering(src, lat, 1, seed=1, path="python")
    on = TransportEngine(
        table, WaterSlab(250.0), DepthDoseGrid(250.0, 10)
    ).run_scattering(src, lat, 200, seed=1, path="python")
    sig_off = off.sigma_x_at_depth(120.0)
    sig_on = on.sigma_x_at_depth(120.0)
    assert sig_off < 1e-9
    assert sig_on > 1.0


def test_rotated_scene_equivalence_deterministic(table) -> None:
    """Scattering off + straggling off: a rigid rotation of beam+slab leaves the
    beam-frame depth dose bit-for-bit invariant (the frame transform is the only
    difference and the ray is deterministic)."""
    lat = _lat(250.0)
    axis = TransportEngine(
        table,
        WaterSlab(250.0),
        DepthDoseGrid(250.0, 10),
        straggling=False,
        scattering=False,
    ).run_scattering(PencilBeamSource(150.0), lat, 1, seed=5, path="python")
    r = _rot_y(35.0) @ _rot_x(20.0)
    d = r @ np.array([0.0, 0.0, 1.0])
    rot = TransportEngine(
        table,
        WaterSlab(250.0, normal=tuple(d)),
        DepthDoseGrid(250.0, 10),
        straggling=False,
        scattering=False,
    ).run_scattering(
        PencilBeamSource(150.0, direction=tuple(d)), lat, 1, seed=5, path="python"
    )
    denom = axis.energy_deposited_mev
    assert abs(rot.energy_deposited_mev - denom) / denom < 1e-12
    assert np.max(np.abs(rot.depth_dose_mev - axis.depth_dose_mev)) / denom < 1e-9


def test_rotated_scene_equivalence_statistical(table) -> None:
    """Scattering on: a rigid rotation of beam+slab conserves total energy to
    round-off and reproduces R80 within statistics (only the azimuthal gauge of
    the reconstructed frame differs, so histories agree only in distribution)."""
    lat = _lat(250.0)
    axis = TransportEngine(
        table, WaterSlab(250.0), DepthDoseGrid(250.0, 10)
    ).run_scattering(PencilBeamSource(150.0), lat, 400, seed=9, path="python")
    r = _rot_y(28.0) @ _rot_x(-17.0)
    d = r @ np.array([0.0, 0.0, 1.0])
    rot = TransportEngine(
        table, WaterSlab(250.0, normal=tuple(d)), DepthDoseGrid(250.0, 10)
    ).run_scattering(
        PencilBeamSource(150.0, direction=tuple(d)), lat, 400, seed=9, path="python"
    )
    denom = axis.energy_deposited_mev
    assert abs(rot.energy_deposited_mev - denom) / denom < 1e-9
    assert abs(_r80_mm(rot) - _r80_mm(axis)) < 1.0


def test_oblique_incidence_traverses_sec_theta_wet(table) -> None:
    """A beam tilted by theta through a lab-fixed +z slab of thickness D crosses
    the same water-equivalent path as a normal beam through D / cos(theta): both
    escape the back face having deposited the same energy (deterministic,
    scattering off). This exercises the material-coordinate plane traversal
    (u = normal . position), which a stopping beam in a thick slab would hide."""
    theta = 30.0
    cos_t = float(np.cos(np.radians(theta)))
    depth = 80.0  # thin enough that a 150 MeV proton escapes the back face
    lat = _lat(300.0)
    eng_kw = dict(straggling=False, scattering=False)
    normal_eff = TransportEngine(
        table, WaterSlab(depth / cos_t), DepthDoseGrid(300.0, 10), **eng_kw
    ).run_scattering(PencilBeamSource(150.0), lat, 1, seed=3, path="python")
    d = _rot_y(theta) @ np.array([0.0, 0.0, 1.0])
    oblique = TransportEngine(
        table, WaterSlab(depth), DepthDoseGrid(300.0, 10), **eng_kw
    ).run_scattering(
        PencilBeamSource(150.0, direction=tuple(d)), lat, 1, seed=3, path="python"
    )
    # both escape with residual energy; the deposited energy (same WET path) must
    # agree to the step-discretisation level
    assert 0.0 < oblique.energy_deposited_mev < 150.0
    e_ref = normal_eff.energy_deposited_mev
    assert abs(oblique.energy_deposited_mev - e_ref) / e_ref < 1e-3


def _rotated(theta_y: float, theta_x: float) -> tuple:
    d = _rot_y(theta_y) @ _rot_x(theta_x) @ np.array([0.0, 0.0, 1.0])
    return tuple(d)


@pytest.mark.warp
def test_warp_cpu_matches_reference_rotated(warp_module, table) -> None:
    """For a rigidly rotated scene the Warp CPU path reproduces the reference
    lateral spread (float32 vs float64 round-off; decision 0001)."""
    d = _rotated(28.0, -17.0)
    lat = _lat(250.0)
    eng = TransportEngine(table, WaterSlab(250.0, normal=d), DepthDoseGrid(250.0, 10))
    src = PencilBeamSource(150.0, direction=d)
    ref = eng.run_scattering(src, lat, 400, seed=9, path="python")
    warp = eng.run_scattering(src, lat, 400, seed=9, path="warp", device="cpu")
    assert abs(ref.sigma_x_at_depth(120.0) - warp.sigma_x_at_depth(120.0)) <= 0.1


@pytest.mark.warp
def test_warp_deterministic_rotation_equivalence(warp_module, table) -> None:
    """Scattering off: the Warp CPU frame transform makes a rigidly rotated scene
    deposit the same total energy as the axis-aligned run (float32 round-off)."""
    lat = _lat(250.0)
    kw = dict(straggling=False, scattering=False)
    axis = TransportEngine(
        table, WaterSlab(250.0), DepthDoseGrid(250.0, 10), **kw
    ).run_scattering(PencilBeamSource(150.0), lat, 1, seed=5, path="warp", device="cpu")
    d = _rotated(35.0, 20.0)
    rot = TransportEngine(
        table, WaterSlab(250.0, normal=d), DepthDoseGrid(250.0, 10), **kw
    ).run_scattering(
        PencilBeamSource(150.0, direction=d), lat, 1, seed=5, path="warp", device="cpu"
    )
    denom = axis.energy_deposited_mev
    assert abs(rot.energy_deposited_mev - denom) / denom < 1e-5


@pytest.mark.cuda
def test_warp_cuda_matches_cpu_rotated(warp_module, cuda_available, table) -> None:
    """A rotated-scene run agrees on Warp CPU and CUDA (deterministic-kernel
    parity, decision 0001)."""
    d = _rotated(28.0, -17.0)
    lat = _lat(250.0)
    eng = TransportEngine(table, WaterSlab(250.0, normal=d), DepthDoseGrid(250.0, 10))
    src = PencilBeamSource(150.0, direction=d)
    cpu = eng.run_scattering(src, lat, 20000, seed=9, path="warp", device="cpu")
    cuda = eng.run_scattering(src, lat, 20000, seed=9, path="warp", device="cuda:0")
    assert abs(cpu.sigma_x_at_depth(120.0) - cuda.sigma_x_at_depth(120.0)) <= 0.02
