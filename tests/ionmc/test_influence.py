"""Tests of beamlet-resolved scoring and sparse influence matrices (decision 0022).

Several pencil-beam beamlets are transported (separately and batched together);
the sum of the per-beamlet doses must equal the batched broad-field dose (the
first V4 gate), the sparse matrix must round-trip and threshold correctly, and
the reference and Warp CPU paths must agree.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.influence import SparseInfluenceMatrix, assemble_influence_matrix
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


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (120, 120, 250), (1.0, 1.0, 1.0), 1.0, origin_mm=(-60.0, -60.0, 0.0)
    )


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(40, 40, 125), origin_mm=(-60.0, -60.0, 0.0), spacing_mm=(3.0, 3.0, 2.0)
    )


def _lat() -> DepthLateralGrid:
    return DepthLateralGrid(250.0, 500, 30.0, 400)


def _beamlets() -> list[PencilBeamSource]:
    return [
        PencilBeamSource(150.0, position_mm=(x, 0.0, 0.0), beamlet=b)
        for b, x in enumerate([-20.0, -7.0, 7.0, 20.0])
    ]


def test_state_concatenation_batches_beamlets(table) -> None:
    """run_scattering_multi transports the concatenation of the per-beamlet states
    (n_histories per beamlet)."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    srcs = _beamlets()
    res = eng.run_scattering_multi(srcs, _lat(), n_histories=50, seed=1, path="python")
    assert res.n_histories == 50 * len(srcs)
    assert abs(res.energy_in_mev - 50 * len(srcs) * 150.0) < 1e-9


def test_beamlet_sum_equals_broadfield(table) -> None:
    """The sum of the per-beamlet influence rows equals the batched broad-field
    dose to round-off (same seed-partitioned histories; the first V4 gate)."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    srcs = _beamlets()
    dose = _dose()
    m = assemble_influence_matrix(
        eng, srcs, _lat(), dose, n_histories=100, seed=1000, path="python"
    )
    combined = eng.run_scattering_multi(
        srcs, _lat(), n_histories=100, seed=1000, path="python", dose_grid=dose
    )
    tot = m.total_dose()
    cb = combined.dose3d_mev
    denom = float(cb.sum())
    assert np.max(np.abs(tot - cb)) / denom < 1e-9
    assert abs(denom - 100 * len(srcs) * 150.0) < 1e-6  # contained, energy conserved


def test_influence_matrix_is_sparse_and_thresholds(table) -> None:
    """The matrix is sparse; thresholding drops low-dose voxels while keeping
    almost all of the energy, and to_dense reproduces the retained rows."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    srcs = _beamlets()
    dose = _dose()
    full = assemble_influence_matrix(
        eng, srcs, _lat(), dose, n_histories=100, seed=1000, path="python"
    )
    thr = assemble_influence_matrix(
        eng,
        srcs,
        _lat(),
        dose,
        n_histories=100,
        seed=1000,
        path="python",
        threshold_frac=0.01,
    )
    assert full.nnz < full.n_beamlets * full.n_voxels  # genuinely sparse
    assert thr.nnz <= full.nnz
    assert float(thr.total_dose().sum()) / float(full.total_dose().sum()) > 0.99
    dense = full.to_dense()
    assert dense.shape == (full.n_beamlets, full.n_voxels)
    np.testing.assert_array_equal(dense[0], full.beamlet_dose_flat(0))


def test_influence_matrix_save_load(table, tmp_path) -> None:
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    m = assemble_influence_matrix(
        eng, _beamlets(), _lat(), _dose(), n_histories=50, seed=7, path="python"
    )
    p = str(tmp_path / "influence.npz")
    m.save(p)
    loaded = SparseInfluenceMatrix.load(p)
    assert loaded.grid_shape == m.grid_shape
    np.testing.assert_array_equal(loaded.indptr, m.indptr)
    np.testing.assert_array_equal(loaded.indices, m.indices)
    np.testing.assert_array_equal(loaded.data, m.data)
    np.testing.assert_array_equal(loaded.beamlet_ids, m.beamlet_ids)


def test_run_scattering_multi_requires_sources(table) -> None:
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    with pytest.raises(ValueError, match="at least one beamlet"):
        eng.run_scattering_multi([], _lat(), n_histories=1, path="python")


@pytest.mark.warp
def test_warp_cpu_beamlet_sum_equals_broadfield(warp_module, table) -> None:
    """On Warp CPU the sum of per-beamlet doses matches the batched broad-field
    dose within the DoseGrid3D float32 budget."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    srcs = _beamlets()
    dose = _dose()
    m = assemble_influence_matrix(
        eng, srcs, _lat(), dose, n_histories=2000, seed=5, path="warp", device="cpu"
    )
    combined = eng.run_scattering_multi(
        srcs,
        _lat(),
        n_histories=2000,
        seed=5,
        path="warp",
        device="cpu",
        dose_grid=dose,
    )
    tot = m.total_dose()
    cb = combined.dose3d_mev
    denom = float(cb.sum())
    assert abs(float(tot.sum()) - denom) / denom < 1e-5
    assert np.max(np.abs(tot - cb)) / float(cb.max()) < 5e-3


@pytest.mark.cuda
def test_warp_cuda_broadfield_matches_cpu(warp_module, cuda_available, table) -> None:
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    srcs = _beamlets()
    dose = _dose()
    cpu = eng.run_scattering_multi(
        srcs,
        _lat(),
        n_histories=20000,
        seed=9,
        path="warp",
        device="cpu",
        dose_grid=dose,
    )
    cuda = eng.run_scattering_multi(
        srcs,
        _lat(),
        n_histories=20000,
        seed=9,
        path="warp",
        device="cuda:0",
        dose_grid=dose,
    )
    denom = float(cpu.dose3d_mev.sum())
    assert abs(float(cuda.dose3d_mev.sum()) - denom) / denom < 1e-4
