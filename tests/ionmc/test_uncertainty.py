"""Tests of batch-based statistical uncertainty for 3-D dose (decision 0024).

``run_scattering_batched`` runs independent history batches into a ``DoseGrid3D``
and reports the per-voxel mean dose and the standard error of the mean (and, with
``score_let``, the LET_d mean and SEM). The tests cover the batch-count guard,
energy conservation of the mean, uncertainty sanity (positive in dose, zero
outside), the optional LET fields, and the ``1/sqrt(N)`` scaling law.
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
)
from ionmc.transport.depth_dose import DepthLateralGrid


@pytest.fixture(scope="module")
def table(pstar_cache_root):
    from ionmc.data import cache

    return load_stopping_table(cache.load_path(MCSQUARE_PSTAR_WATER, pstar_cache_root))


def _lat(depth: float) -> DepthLateralGrid:
    return DepthLateralGrid(depth, int(depth * 2), 30.0, 60)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (21, 21, 400), (1.0, 1.0, 0.5), 1.0, origin_mm=(-10.5, -10.5, 0.0)
    )


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(21, 21, 100), origin_mm=(-10.5, -10.5, 0.0), spacing_mm=(1.0, 1.0, 2.0)
    )


def test_batch_count_guard(table) -> None:
    eng = TransportEngine(table, _box(), DepthDoseGrid(200.0, 10))
    with pytest.raises(ValueError, match="at least two batches"):
        eng.run_scattering_batched(
            PencilBeamSource(150.0),
            _lat(200.0),
            _dose(),
            40,
            n_batches=1,
            path="python",
        )


def test_batch_mean_energy_conservation(table) -> None:
    """The batch-mean total dose equals the per-batch input energy for a contained
    beam (each batch conserves energy; so does their mean)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    r = eng.run_scattering_batched(
        PencilBeamSource(150.0),
        _lat(200.0),
        _dose(),
        40,
        n_batches=4,
        seed=3,
        path="python",
    )
    expected = r.histories_per_batch * 150.0
    assert abs(float(r.mean_dose3d_mev.sum()) - expected) / expected < 1e-9


def test_batch_uncertainty_sanity(table) -> None:
    """SEM is positive in high-dose voxels (stochastic transport), zero where there
    is no dose, and the high-dose mean relative uncertainty is a small fraction."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(200.0, 10))  # scattering on
    r = eng.run_scattering_batched(
        PencilBeamSource(150.0),
        _lat(200.0),
        _dose(),
        400,
        n_batches=8,
        seed=5,
        path="python",
    )
    assert r.mean_dose3d_mev.shape == (21, 21, 100)
    assert r.standard_error_mev.shape == (21, 21, 100)
    hi = r.mean_dose3d_mev > 0.5 * r.mean_dose3d_mev.max()
    assert np.all(r.standard_error_mev[hi] > 0.0)
    assert np.all(r.standard_error_mev[r.mean_dose3d_mev == 0.0] == 0.0)
    rse = r.relative_standard_error
    assert np.all(rse >= 0.0)
    mru = r.mean_relative_uncertainty(0.5)
    assert 0.0 < mru < 1.0


def test_batch_let_fields(table) -> None:
    """``score_let`` adds LET_d mean and SEM grids; otherwise they are None."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(200.0, 10), straggling=False, scattering=False
    )
    no_let = eng.run_scattering_batched(
        PencilBeamSource(150.0),
        _lat(200.0),
        _dose(),
        40,
        n_batches=4,
        seed=3,
        path="python",
    )
    assert no_let.mean_let_d_kev_um is None
    assert no_let.standard_error_let_kev_um is None
    with_let = eng.run_scattering_batched(
        PencilBeamSource(150.0),
        _lat(200.0),
        _dose(),
        40,
        n_batches=4,
        seed=3,
        path="python",
        score_let=True,
    )
    assert with_let.mean_let_d_kev_um.shape == (21, 21, 100)
    assert with_let.standard_error_let_kev_um.shape == (21, 21, 100)
    # entrance LET_d ~ 0.5 keV/um near the beam axis at shallow depth
    assert with_let.mean_let_d_kev_um.max() > 0.0


@pytest.mark.warp
def test_se_scaling_one_over_sqrt_n(warp_module, table) -> None:
    """The standard error of the mean scales as 1/sqrt(N): quadrupling the
    histories halves the high-dose mean relative SEM (decision 0024). Uses a shared
    high-dose mask so both runs average over the identical voxel set."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(200.0, 10))  # stochastic
    src = PencilBeamSource(150.0)
    lat, dose = _lat(200.0), _dose()
    rN = eng.run_scattering_batched(
        src, lat, dose, 40000, n_batches=20, seed=100, path="warp", device="cpu"
    )
    r4 = eng.run_scattering_batched(
        src, lat, dose, 160000, n_batches=20, seed=100, path="warp", device="cpu"
    )
    mask = r4.mean_dose3d_mev > 0.2 * r4.mean_dose3d_mev.max()
    uN = float(np.mean(rN.relative_standard_error[mask]))
    u4 = float(np.mean(r4.relative_standard_error[mask]))
    ratio = u4 / uN
    assert 0.35 <= ratio <= 0.71  # 1/sqrt(4) = 0.5 within a generous stat band
