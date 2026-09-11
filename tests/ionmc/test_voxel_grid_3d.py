"""Tests of the 3-D voxel grid with ray/voxel DDA traversal (decision 0020).

``VoxelGrid3D`` gives the scattering path a true Nx*Ny*Nz voxel geometry with an
Amanatides-Woo DDA. The tests confirm: geometry bookkeeping; a single-column grid
reproduces the 1-D ``VoxelSlab`` (the reduction safety anchor); a homogeneous box
reproduces ``WaterSlab``; an oblique beam accumulates the correct water-equivalent
path (vs an independent Siddon oracle) and exits the box sides; energy is
conserved for a contained beam; and the reference and Warp CPU paths agree.
All checks use the reference Python path except those marked ``warp``/``cuda``.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc import materials, particles
from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.tabulated_stopping_power import TabulatedStoppingPower
from ionmc.transport import (
    DepthDoseGrid,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
    WaterSlab,
)
from ionmc.transport.depth_dose import DepthLateralGrid
from ionmc.transport.geometry import VoxelSlab


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(depth, int(depth * 2), 30.0, 400)


def _range_mm(table, e0: float) -> float:
    model = TabulatedStoppingPower(table, materials.WATER, particles.PROTON, "numpy")
    return float(model.csda_range(e0)[0]) * 10.0


def _siddon_wet_mm(grid: VoxelGrid3D, p0, d, path_mm: float) -> float:
    """Independent water-equivalent path integral ``integral rho dl`` along the ray
    from ``p0`` in unit direction ``d`` for ``path_mm`` (a plain per-mm sampler,
    algorithmically distinct from the transport DDA)."""
    o = np.asarray(grid.origin_mm)
    h = np.asarray(grid.spacing_mm)
    n = np.array([grid.nx, grid.ny, grid.nz])
    rho = grid.density_g_per_cm3
    p0 = np.asarray(p0, dtype=float)
    d = np.asarray(d, dtype=float)
    d = d / np.linalg.norm(d)
    ds = 1.0e-3
    wet = 0.0
    s = 0.5 * ds
    while s < path_mm:
        x = p0 + s * d
        idx = np.floor((x - o) / h).astype(int)
        if np.all(idx >= 0) and np.all(idx < n):
            wet += float(rho[idx[0], idx[1], idx[2]]) * ds
        s += ds
    return wet


def test_grid_geometry_basics() -> None:
    g = VoxelGrid3D.uniform((2, 3, 4), (1.0, 1.0, 2.0), 1.0, origin_mm=(-1.0, 0.0, 0.0))
    assert (g.nx, g.ny, g.nz) == (2, 3, 4)
    assert g.n_voxels == 24
    assert g.bbox_hi_mm == (1.0, 3.0, 8.0)
    d = np.arange(24, dtype=float).reshape(2, 3, 4) + 1.0
    gg = VoxelGrid3D(d, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    assert gg.density_flat()[(1 * gg.ny + 2) * gg.nz + 3] == d[1, 2, 3]
    with pytest.raises(ValueError, match="positive"):
        VoxelGrid3D.uniform((1, 1, 1), (1.0, 1.0, 0.0), 1.0)
    with pytest.raises(ValueError, match="positive"):
        VoxelGrid3D(np.zeros((1, 1, 1)), (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))


def test_from_voxel_slab_requires_uniform_z() -> None:
    good = VoxelSlab(
        z_boundaries_mm=np.arange(0, 51, 10.0),
        density_g_per_cm3=np.array([1.0, 1.5, 1.0, 0.8, 1.0]),
    )
    g = VoxelGrid3D.from_voxel_slab(good, transverse_extent_mm=40.0)
    assert (g.nx, g.ny, g.nz) == (1, 1, 5)
    np.testing.assert_allclose(g.density_g_per_cm3[0, 0, :], [1.0, 1.5, 1.0, 0.8, 1.0])
    bad = VoxelSlab.from_layers([(20.0, 1.0), (10.0, 1.6)])  # non-uniform spacing
    with pytest.raises(ValueError, match="uniform z-spacing"):
        VoxelGrid3D.from_voxel_slab(bad, transverse_extent_mm=40.0)


def test_reduction_to_voxel_slab(table) -> None:
    """A single-column grid reproduces the 1-D VoxelSlab depth dose (deterministic
    straight ray) to round-off -- the strict-generalization safety anchor."""
    z = np.arange(0, 301, 10.0)
    rho = 1.0 + 0.3 * np.cos(np.arange(len(z) - 1))  # distinct adjacent, all > 0
    slab = VoxelSlab(z_boundaries_mm=z, density_g_per_cm3=rho)
    grid = VoxelGrid3D.from_voxel_slab(slab, transverse_extent_mm=60.0)
    lat = _lat(300.0)
    src = PencilBeamSource(150.0)
    kw = dict(straggling=False, scattering=False)
    a = TransportEngine(table, slab, DepthDoseGrid(300.0, 10), **kw).run_scattering(
        src, lat, 1, seed=4, path="python"
    )
    b = TransportEngine(table, grid, DepthDoseGrid(300.0, 10), **kw).run_scattering(
        src, lat, 1, seed=4, path="python"
    )
    denom = np.sum(a.depth_dose_mev)
    cum = np.max(np.abs(np.cumsum(b.depth_dose_mev) - np.cumsum(a.depth_dose_mev)))
    assert cum / denom < 1e-9
    assert abs(b.energy_deposited_mev - a.energy_deposited_mev) / denom < 1e-12


def _r80(res) -> float:
    dd = res.depth_dose_mev
    c = res.grid.depth_centers_mm
    peak = float(dd.max())
    i = int(dd.argmax())
    lvl = 0.8 * peak
    for k in range(i, len(dd) - 1):
        if dd[k] >= lvl >= dd[k + 1]:
            f = (dd[k] - lvl) / (dd[k] - dd[k + 1])
            return float(c[k] + f * (c[k + 1] - c[k]))
    return float(c[i])


def test_homogeneous_box_matches_water_slab(table) -> None:
    """A homogeneous water box reproduces the WaterSlab depth dose to the step-
    partition discretization level (the box clips at every 1 mm z-face while the
    merged WaterSlab takes physics-limited steps, so they agree at the midpoint-
    rule O(step) level, not bit-for-bit; deterministic)."""
    lat = _lat(300.0)
    src = PencilBeamSource(150.0)
    kw = dict(straggling=False, scattering=False)
    box = VoxelGrid3D.uniform(
        (80, 80, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-40.0, -40.0, 0.0)
    )
    b = TransportEngine(table, box, DepthDoseGrid(300.0, 10), **kw).run_scattering(
        src, lat, 1, seed=4, path="python"
    )
    s = TransportEngine(
        table, WaterSlab(300.0), DepthDoseGrid(300.0, 10), **kw
    ).run_scattering(src, lat, 1, seed=4, path="python")
    denom = np.sum(s.depth_dose_mev)
    cum = np.max(np.abs(np.cumsum(b.depth_dose_mev) - np.cumsum(s.depth_dose_mev)))
    assert abs(b.energy_deposited_mev - s.energy_deposited_mev) / denom < 1e-4
    assert abs(_r80(b) - _r80(s)) < 1e-2
    assert cum / denom < 3e-3


def test_oblique_uniform_box_wet_range(table) -> None:
    """An oblique beam through a uniform rho box accumulates WET = rho * path, so
    its beam-depth range is R_water/rho independent of the tilt (deterministic).
    Exercises 3-D traversal on all axes."""
    rho = 1.2
    r_water = _range_mm(table, 150.0)
    lat = _lat(320.0)
    kw = dict(straggling=False, scattering=False)
    box = VoxelGrid3D.uniform(
        (200, 200, 320), (1.0, 1.0, 1.0), rho, origin_mm=(-100.0, -100.0, 0.0)
    )
    eng = TransportEngine(table, box, DepthDoseGrid(320.0, 10), **kw)
    for theta in (0.0, 15.0):
        a = np.radians(theta)
        d = (np.sin(a), 0.0, np.cos(a))
        res = eng.run_scattering(
            PencilBeamSource(150.0, direction=d), lat, 1, seed=3, path="python"
        )
        assert res.n_stopped == 1  # contained in the wide box
        assert abs(res.range_mean_mm - r_water / rho) / (r_water / rho) < 0.02


def test_oblique_heterogeneous_wet_vs_siddon(table) -> None:
    """An oblique beam through a grid with an off-axis dense insert stops at the
    beam-depth where the independent Siddon WET integral reaches the water CSDA
    range (deterministic; a genuine 3-D interior-crossing physics check)."""
    r_water = _range_mm(table, 150.0)
    # water box with a dense slab insert spanning x>0 for a mid-z band
    dens = np.ones((160, 160, 300))
    dens[80:, :, 80:180] = 1.7  # dense quadrant/band (off the entry axis)
    box = VoxelGrid3D(dens, origin_mm=(-80.0, -80.0, 0.0), spacing_mm=(1.0, 1.0, 1.0))
    lat = _lat(300.0)
    a = np.radians(18.0)
    d = (np.sin(a), 0.0, np.cos(a))
    res = TransportEngine(
        table, box, DepthDoseGrid(300.0, 10), straggling=False, scattering=False
    ).run_scattering(
        PencilBeamSource(150.0, direction=d), lat, 1, seed=3, path="python"
    )
    assert res.n_stopped == 1
    # Siddon: the path length at which the independent WET integral reaches r_water
    p0 = (0.0, 0.0, 0.0)
    lo, hi = 0.0, 300.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _siddon_wet_mm(box, p0, d, mid) < r_water:
            lo = mid
        else:
            hi = mid
    siddon_range = 0.5 * (lo + hi)
    assert abs(res.range_mean_mm - siddon_range) / siddon_range < 1e-2


def test_oblique_beam_exits_box_side(table) -> None:
    """A steeply tilted beam through a narrow box exits a side face (escaped, not
    stopped), depositing less than its full energy -- the DDA handles lateral
    escape."""
    lat = _lat(300.0)
    box = VoxelGrid3D.uniform(
        (60, 60, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-30.0, -30.0, 0.0)
    )
    a = np.radians(25.0)
    d = (np.sin(a), 0.0, np.cos(a))
    res = TransportEngine(
        table, box, DepthDoseGrid(300.0, 10), straggling=False, scattering=False
    ).run_scattering(
        PencilBeamSource(150.0, direction=d), lat, 1, seed=3, path="python"
    )
    assert res.n_stopped == 0
    assert 0.0 < res.energy_deposited_mev < 150.0


def test_energy_conservation_contained(table) -> None:
    """A beam contained in a wide box deposits all its energy (scattering on)."""
    lat = _lat(300.0)
    box = VoxelGrid3D.uniform(
        (120, 120, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-60.0, -60.0, 0.0)
    )
    res = TransportEngine(table, box, DepthDoseGrid(300.0, 10)).run_scattering(
        PencilBeamSource(150.0), lat, 200, seed=7, path="python"
    )
    assert abs(res.energy_balance) <= 1e-9


def test_depth_dose_path_rejects_grid(table) -> None:
    box = VoxelGrid3D.uniform((2, 2, 10), (1.0, 1.0, 10.0), 1.0)
    eng = TransportEngine(table, box, DepthDoseGrid(100.0, 10))
    with pytest.raises(ValueError, match="VoxelGrid3D"):
        eng.run(PencilBeamSource(150.0), 1, path="python")


@pytest.mark.warp
def test_warp_cpu_reduction(warp_module, table) -> None:
    """The Warp CPU grid kernel reproduces the reference grid path on a
    deterministic +z homogeneous box (float32 round-off)."""
    lat = _lat(300.0)
    box = VoxelGrid3D.uniform(
        (60, 60, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-30.0, -30.0, 0.0)
    )
    eng = TransportEngine(
        table, box, DepthDoseGrid(300.0, 10), straggling=False, scattering=False
    )
    src = PencilBeamSource(150.0)
    ref = eng.run_scattering(src, lat, 1, seed=4, path="python")
    cpu = eng.run_scattering(src, lat, 1, seed=4, path="warp", device="cpu")
    denom = np.sum(ref.depth_dose_mev)
    cum = np.max(np.abs(np.cumsum(cpu.depth_dose_mev) - np.cumsum(ref.depth_dose_mev)))
    assert cum / denom < 5e-4


@pytest.mark.warp
def test_warp_cpu_matches_reference_scattering(warp_module, table) -> None:
    """Warp CPU vs reference on the grid path with scattering on: lateral spread
    agrees (float32 round-off, decision 0001)."""
    lat = _lat(250.0)
    box = VoxelGrid3D.uniform(
        (100, 100, 250), (1.0, 1.0, 1.0), 1.0, origin_mm=(-50.0, -50.0, 0.0)
    )
    eng = TransportEngine(table, box, DepthDoseGrid(250.0, 10))
    src = PencilBeamSource(150.0)
    ref = eng.run_scattering(src, lat, 400, seed=11, path="python")
    cpu = eng.run_scattering(src, lat, 400, seed=11, path="warp", device="cpu")
    assert abs(ref.sigma_x_at_depth(120.0) - cpu.sigma_x_at_depth(120.0)) <= 0.1


@pytest.mark.cuda
def test_warp_cuda_matches_cpu(warp_module, cuda_available, table) -> None:
    """A grid run agrees on Warp CPU and CUDA (deterministic-kernel parity)."""
    lat = _lat(250.0)
    box = VoxelGrid3D.uniform(
        (100, 100, 250), (1.0, 1.0, 1.0), 1.0, origin_mm=(-50.0, -50.0, 0.0)
    )
    eng = TransportEngine(table, box, DepthDoseGrid(250.0, 10))
    src = PencilBeamSource(150.0)
    cpu = eng.run_scattering(src, lat, 20000, seed=9, path="warp", device="cpu")
    cuda = eng.run_scattering(src, lat, 20000, seed=9, path="warp", device="cuda:0")
    assert abs(cpu.sigma_x_at_depth(120.0) - cuda.sigma_x_at_depth(120.0)) <= 0.02
