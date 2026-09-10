"""Tests of scoring grids decoupled from the transport grid (decision 0019).

The scoring grid's resolution AND alignment (depth origin, lateral centre) are
independent of the transport voxel grid. The integral dose is invariant under
scoring-grid changes as long as the grid contains the dose; the depth-dose curve
agrees across resolutions; and lateral spread is invariant under a lateral shift.
All checks use the reference Python path (no GPU); cross-backend parity is a
host-runner validation gate (``validation/v3_scoring_grid.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import DepthDoseGrid, PencilBeamSource, TransportEngine, WaterSlab
from ionmc.transport.depth_dose import DepthLateralGrid


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


def test_grid_origin_shifts_edges_and_centers() -> None:
    base = DepthDoseGrid(100.0, 10)
    shifted = DepthDoseGrid(100.0, 10, origin_mm=15.0)
    np.testing.assert_allclose(shifted.edges_mm, base.edges_mm + 15.0)
    np.testing.assert_allclose(shifted.centers_mm, base.centers_mm + 15.0)
    lat = DepthLateralGrid(
        100.0, 10, 20.0, 40, depth_origin_mm=5.0, lateral_center_mm=3.0
    )
    np.testing.assert_allclose(lat.depth_centers_mm[0], 5.0 + 0.5 * lat.depth_bin_mm)
    np.testing.assert_allclose(lat.lateral_lo_mm, 3.0 - 20.0)
    assert abs(lat.lateral_centers_mm.mean() - 3.0) < 1e-9


def _run(table, grid, **eng_kw):
    return TransportEngine(
        table, WaterSlab(300.0), DepthDoseGrid(300.0, 10), **eng_kw
    ).run_scattering(PencilBeamSource(150.0), grid, 1, seed=4, path="python")


def test_integral_dose_invariant_under_resolution(table) -> None:
    """Deterministic straight ray (scattering off): the total deposited energy
    and the cumulative depth dose are invariant when the depth resolution changes
    (the deposition splits energy by depth overlap)."""
    kw = dict(straggling=False, scattering=False)
    coarse = _run(table, DepthLateralGrid(300.0, 150, 25.0, 200), **kw)
    fine = _run(table, DepthLateralGrid(300.0, 600, 25.0, 200), **kw)
    assert coarse.energy_deposited_mev > 100.0
    rel = abs(coarse.energy_deposited_mev - fine.energy_deposited_mev) / (
        fine.energy_deposited_mev
    )
    assert rel < 1e-12
    # the fine grid (600) nests exactly 4:1 in the coarse (150) over the same
    # extent, so summing each group of 4 fine depth bins must equal the coarse
    # bin (energy per depth interval is conserved by the overlap deposition)
    rebinned = fine.depth_dose_mev.reshape(150, 4).sum(axis=1)
    denom = np.sum(coarse.depth_dose_mev)
    assert np.max(np.abs(rebinned - coarse.depth_dose_mev)) / denom < 1e-12


def test_integral_dose_invariant_under_alignment(table) -> None:
    """The total deposited energy is invariant under a depth-origin and
    lateral-centre shift of the scoring grid, as long as it still contains the
    dose (deterministic straight ray)."""
    kw = dict(straggling=False, scattering=False)
    base = _run(table, DepthLateralGrid(300.0, 300, 25.0, 200), **kw)
    shifted = _run(
        table,
        DepthLateralGrid(
            320.0, 320, 30.0, 240, depth_origin_mm=-20.0, lateral_center_mm=4.0
        ),
        **kw,
    )
    e0 = base.energy_deposited_mev
    assert abs(shifted.energy_deposited_mev - e0) / e0 < 1e-12


def test_partial_coverage_drops_outside_dose(table) -> None:
    """A scoring grid that starts past the entrance captures strictly less energy
    (dose upstream of the origin is dropped), confirming the origin truly shifts
    the scored window."""
    kw = dict(straggling=False, scattering=False)
    full = _run(table, DepthLateralGrid(300.0, 300, 25.0, 200), **kw)
    downstream = _run(
        table, DepthLateralGrid(200.0, 200, 25.0, 200, depth_origin_mm=100.0), **kw
    )
    assert 0.0 < downstream.energy_deposited_mev < full.energy_deposited_mev


def test_lateral_shift_preserves_sigma_x(table) -> None:
    """With scattering on, the lateral spread sigma_x(z) is invariant under a
    lateral-centre shift of the grid (the centres move with it), as long as the
    grid still contains the lateral distribution."""
    eng = TransportEngine(table, WaterSlab(250.0), DepthDoseGrid(250.0, 10))
    src = PencilBeamSource(150.0)
    centred = eng.run_scattering(
        src, DepthLateralGrid(250.0, 500, 30.0, 600), 400, seed=7, path="python"
    )
    shifted = eng.run_scattering(
        src,
        DepthLateralGrid(250.0, 500, 30.0, 600, lateral_center_mm=5.0),
        400,
        seed=7,
        path="python",
    )
    for z in (90.0, 130.0):
        assert abs(centred.sigma_x_at_depth(z) - shifted.sigma_x_at_depth(z)) < 1e-9


def test_lateral_window_shift_selects_region(table) -> None:
    """The lateral centre truly shifts the scored window (discriminates the
    lateral origin, which the translation-invariant sigma_x cannot): with
    scattering off the beam stays on the x=0 axis, so a window centred on the beam
    captures the dose while one shifted off the beam captures none."""
    kw = dict(straggling=False, scattering=False)
    centred = _run(table, DepthLateralGrid(300.0, 300, 3.0, 60), **kw)  # x in [-3, 3]
    off = _run(
        table, DepthLateralGrid(300.0, 300, 3.0, 60, lateral_center_mm=10.0), **kw
    )  # x in [7, 13], excludes the beam
    assert centred.energy_deposited_mev > 100.0
    assert off.energy_deposited_mev == 0.0


def test_lateral_shift_reports_beam_centre(table) -> None:
    """On a laterally shifted grid the energy-weighted mean lateral position is the
    beam's true axis (x=0), not the grid centre -- so the deposition honours the
    lateral origin rather than ignoring it (deterministic straight ray)."""
    kw = dict(straggling=False, scattering=False)
    shifted = _run(
        table, DepthLateralGrid(300.0, 300, 30.0, 600, lateral_center_mm=8.0), **kw
    )
    w = shifted.edep_zx_mev.sum(axis=0)
    centers = shifted.grid.lateral_centers_mm
    mean_x = float((w * centers).sum() / w.sum())
    assert abs(mean_x) < 0.2  # beam is at x=0, far from the grid centre (8 mm)


@pytest.mark.warp
def test_warp_cpu_origin_partial_coverage(warp_module, table) -> None:
    """The Warp CPU scattering kernel honours the depth origin: a grid starting
    past the entrance captures strictly less than one covering the full dose
    (discriminating -- if the origin were ignored both would capture the same)."""
    eng = TransportEngine(
        table,
        WaterSlab(300.0),
        DepthDoseGrid(300.0, 10),
        straggling=False,
        scattering=False,
    )
    src = PencilBeamSource(150.0)
    full = eng.run_scattering(
        src,
        DepthLateralGrid(300.0, 300, 25.0, 200),
        1,
        seed=4,
        path="warp",
        device="cpu",
    )
    downstream = eng.run_scattering(
        src,
        DepthLateralGrid(200.0, 200, 25.0, 200, depth_origin_mm=100.0),
        1,
        seed=4,
        path="warp",
        device="cpu",
    )
    assert 0.0 < downstream.energy_deposited_mev < full.energy_deposited_mev


@pytest.mark.cuda
def test_warp_cuda_matches_cpu_shifted_grid(warp_module, cuda_available, table) -> None:
    """A shifted scoring grid agrees on Warp CPU and CUDA (deterministic-kernel
    parity, decision 0001)."""
    eng = TransportEngine(table, WaterSlab(250.0), DepthDoseGrid(250.0, 10))
    src = PencilBeamSource(150.0)
    grid = DepthLateralGrid(
        250.0, 500, 30.0, 600, depth_origin_mm=-10.0, lateral_center_mm=4.0
    )
    cpu = eng.run_scattering(src, grid, 20000, seed=9, path="warp", device="cpu")
    cuda = eng.run_scattering(src, grid, 20000, seed=9, path="warp", device="cuda:0")
    assert abs(cpu.sigma_x_at_depth(120.0) - cuda.sigma_x_at_depth(120.0)) <= 0.02
