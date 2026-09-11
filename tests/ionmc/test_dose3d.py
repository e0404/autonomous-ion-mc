"""Tests of lab-frame 3-D dose scoring on the voxel-grid path (decision 0021).

``DoseGrid3D`` accumulates deposited energy per lab voxel; the 3-D voxel-grid
transport path deposits each step at its lab midpoint. The tests confirm the
grid bookkeeping, energy conservation into the dose grid, depth-dose consistency
with the beam-frame scoring, grid-independence, the non-grid rejection, and
reference/Warp CPU agreement.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.transport import (
    DepthDoseGrid,
    DoseGrid3D,
    PencilBeamSource,
    TransportEngine,
    VoxelGrid3D,
    WaterSlab,
)
from ionmc.transport.depth_dose import DepthLateralGrid


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(depth, int(depth * 2), 30.0, 400)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (120, 120, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-60.0, -60.0, 0.0)
    )


def _dose(nzv: int = 150, sp: float = 2.0) -> DoseGrid3D:
    n = int(120 / (sp if sp else 2.0))
    return DoseGrid3D(
        shape=(n, n, nzv),
        origin_mm=(-60.0, -60.0, 0.0),
        spacing_mm=(sp, sp, 300.0 / nzv),
    )


def test_dose_grid_basics() -> None:
    g = DoseGrid3D(shape=(2, 3, 4), spacing_mm=(1.0, 2.0, 2.0))
    assert (g.nx, g.ny, g.nz) == (2, 3, 4)
    assert g.n_voxels == 24
    assert abs(g.voxel_volume_cm3 - 4.0e-3) < 1e-12  # 1*2*2 mm^3 -> cm^3
    assert g.empty().shape == (2, 3, 4)
    d = g.dose_gy(np.ones((2, 3, 4)), density_g_per_cm3=1.0)
    assert np.all(d > 0.0)
    with pytest.raises(ValueError, match="positive"):
        DoseGrid3D(shape=(1, 1, 1), spacing_mm=(1.0, 1.0, 0.0))
    with pytest.raises(ValueError, match="at least one voxel"):
        DoseGrid3D(shape=(0, 1, 1))


def test_dose_conserves_energy_deterministic(table) -> None:
    """A deterministic contained beam deposits all its energy into the dose grid
    (dose sum == deposited energy, exactly)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(300.0, 10), straggling=False, scattering=False
    )
    res = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(300.0),
        1,
        seed=4,
        path="python",
        dose_grid=_dose(),
    )
    assert res.dose3d_mev is not None
    assert abs(float(res.dose3d_mev.sum()) - res.energy_deposited_mev) < 1e-9


def test_dose_conserves_energy_scattering(table) -> None:
    """A contained scattering beam conserves energy into the dose grid."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(300.0, 10))
    res = eng.run_scattering(
        PencilBeamSource(150.0),
        _lat(300.0),
        200,
        seed=7,
        path="python",
        dose_grid=_dose(),
    )
    assert abs(float(res.dose3d_mev.sum()) - res.energy_deposited_mev) < 1e-9
    assert abs(res.energy_balance) <= 1e-9


def _r80_from_profile(centers: np.ndarray, prof: np.ndarray) -> float:
    peak = float(prof.max())
    i = int(prof.argmax())
    lvl = 0.8 * peak
    for k in range(i, len(prof) - 1):
        if prof[k] >= lvl >= prof[k + 1]:
            f = (prof[k] - lvl) / (prof[k] - prof[k + 1])
            return float(centers[k] + f * (centers[k + 1] - centers[k]))
    return float(centers[i])


def test_dose_depth_profile_matches_beam_frame(table) -> None:
    """The 3-D dose projected onto the beam axis (z-marginal) reproduces the
    beam-frame depth dose R80 to within a dose voxel (deterministic +z beam)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(300.0, 10), straggling=False, scattering=False
    )
    dose = _dose(nzv=300, sp=2.0)  # 1 mm z-voxels
    res = eng.run_scattering(
        PencilBeamSource(150.0), _lat(300.0), 1, seed=4, path="python", dose_grid=dose
    )
    zc, prof = dose.axis_marginals(res.dose3d_mev)
    dd = res.depth_dose_mev
    ddz = res.grid.depth_centers_mm
    r80_3d = _r80_from_profile(zc, prof)
    r80_bf = _r80_from_profile(ddz, dd)
    assert abs(r80_3d - r80_bf) <= 1.0  # within a 1 mm dose voxel


def test_dose_grid_independence(table) -> None:
    """The total 3-D dose is invariant to the dose grid's resolution and alignment
    (it only rebins the same deposits), for a contained beam."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(300.0, 10), straggling=False, scattering=False
    )
    src = PencilBeamSource(150.0)
    lat = _lat(300.0)
    a = eng.run_scattering(
        src, lat, 1, seed=4, path="python", dose_grid=_dose(150, 2.0)
    )
    b = eng.run_scattering(
        src,
        lat,
        1,
        seed=4,
        path="python",
        dose_grid=DoseGrid3D(
            shape=(30, 30, 100),
            origin_mm=(-45.0, -45.0, -10.0),
            spacing_mm=(3.0, 3.0, 3.0),
        ),
    )
    ta, tb = float(a.dose3d_mev.sum()), float(b.dose3d_mev.sum())
    assert abs(ta - tb) / ta < 1e-12


def test_dose_grid_rejected_on_non_grid_geometry(table) -> None:
    eng = TransportEngine(table, WaterSlab(250.0), DepthDoseGrid(250.0, 10))
    with pytest.raises(ValueError, match="VoxelGrid3D transport path"):
        eng.run_scattering(
            PencilBeamSource(150.0),
            _lat(250.0),
            1,
            seed=1,
            path="python",
            dose_grid=_dose(),
        )


@pytest.mark.warp
def test_warp_cpu_dose_matches_reference_deterministic(warp_module, table) -> None:
    """Deterministic +z: the Warp CPU dose grid matches the reference (float32
    point-deposition round-off)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(300.0, 10), straggling=False, scattering=False
    )
    src = PencilBeamSource(150.0)
    lat = _lat(300.0)
    ref = eng.run_scattering(src, lat, 1, seed=4, path="python", dose_grid=_dose())
    cpu = eng.run_scattering(
        src, lat, 1, seed=4, path="warp", device="cpu", dose_grid=_dose()
    )
    assert (
        abs(float(cpu.dose3d_mev.sum()) - float(ref.dose3d_mev.sum()))
        / float(ref.dose3d_mev.sum())
        < 1e-5
    )
    assert (
        np.max(np.abs(cpu.dose3d_mev - ref.dose3d_mev)) / float(ref.dose3d_mev.max())
        < 1e-3
    )


@pytest.mark.cuda
def test_warp_cuda_dose_matches_cpu(warp_module, cuda_available, table) -> None:
    eng = TransportEngine(table, _box(), DepthDoseGrid(300.0, 10))
    src = PencilBeamSource(150.0)
    lat = _lat(300.0)
    cpu = eng.run_scattering(
        src, lat, 20000, seed=9, path="warp", device="cpu", dose_grid=_dose()
    )
    cuda = eng.run_scattering(
        src, lat, 20000, seed=9, path="warp", device="cuda:0", dose_grid=_dose()
    )
    denom = float(cpu.dose3d_mev.sum())
    assert abs(float(cuda.dose3d_mev.sum()) - denom) / denom < 1e-4
