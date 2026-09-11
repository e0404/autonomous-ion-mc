"""Tests of beamlet-resolved (planning-aware) uncertainty (decision 0025).

``assemble_influence_matrix_batched`` transports each beamlet over independent
history batches and stores the per-voxel mean dose plus its standard error
(``data_sigma``) in the sparse influence matrix. The tests cover the sigma/data
alignment and accessors, energy conservation of the summed means, the
missing-sigma guard on the exact single-run matrix, save/load round-tripping,
and the per-beamlet ``1/sqrt(N)`` scaling.
"""

from __future__ import annotations

import numpy as np
import pytest

from ionmc.data import MCSQUARE_PSTAR_WATER
from ionmc.data.stopping_tables import load_stopping_table
from ionmc.influence import (
    assemble_influence_matrix,
    assemble_influence_matrix_batched,
)
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


def _lat() -> DepthLateralGrid:
    return DepthLateralGrid(250.0, 500, 30.0, 400)


def _box() -> VoxelGrid3D:
    return VoxelGrid3D.uniform(
        (60, 60, 300), (1.0, 1.0, 1.0), 1.0, origin_mm=(-30.0, -30.0, 0.0)
    )


def _dose() -> DoseGrid3D:
    return DoseGrid3D(
        shape=(30, 30, 125), origin_mm=(-30.0, -30.0, 0.0), spacing_mm=(2.0, 2.0, 2.0)
    )


def _beamlets() -> list[PencilBeamSource]:
    return [
        PencilBeamSource(150.0, position_mm=(x, 0.0, 0.0), beamlet=b)
        for b, x in enumerate([-8.0, 0.0, 8.0])
    ]


def test_batched_matrix_alignment_and_accessors(table) -> None:
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    m = assemble_influence_matrix_batched(
        eng, _beamlets(), _lat(), _dose(), 200, n_batches=4, seed=1000, path="python"
    )
    assert m.data_sigma is not None
    assert m.data_sigma.shape == m.data.shape  # aligned per nnz
    assert np.all(m.data_sigma >= 0.0)
    # beamlet_sigma_flat is zero off the beamlet's own voxels
    s0 = m.beamlet_sigma_flat(0)
    lo, hi = int(m.indptr[0]), int(m.indptr[1])
    off = np.ones(m.n_voxels, dtype=bool)
    off[m.indices[lo:hi]] = False
    assert np.all(s0[off] == 0.0)
    # total_sigma = sqrt(sum of per-beamlet variances) has the grid shape
    assert m.total_sigma().shape == (30, 30, 125)
    # per-beamlet quality metric is a small positive fraction
    for b in range(m.n_beamlets):
        assert 0.0 < m.beamlet_relative_uncertainty(b, 0.5) < 1.0


def test_summed_beamlet_means_conserve_energy(table) -> None:
    """The sum of the per-beamlet mean doses equals the per-batch input energy of
    the contained beamlets (each batched mean conserves energy)."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(250.0, 10), straggling=False, scattering=False
    )
    m = assemble_influence_matrix_batched(
        eng, _beamlets(), _lat(), _dose(), 40, n_batches=4, seed=1000, path="python"
    )
    per_batch = 40 // 4
    expected = len(_beamlets()) * per_batch * 150.0
    assert abs(float(m.total_dose().sum()) - expected) / expected < 1e-9


def test_single_run_matrix_has_no_sigma(table) -> None:
    """The exact single-run matrix (decision 0022) carries no uncertainty; the
    per-beamlet uncertainty accessors raise."""
    eng = TransportEngine(
        table, _box(), DepthDoseGrid(250.0, 10), straggling=False, scattering=False
    )
    m = assemble_influence_matrix(
        eng, _beamlets(), _lat(), _dose(), 4, seed=1000, path="python"
    )
    assert m.data_sigma is None
    with pytest.raises(ValueError, match="no per-beamlet uncertainty"):
        m.beamlet_sigma_flat(0)
    with pytest.raises(ValueError, match="no per-beamlet uncertainty"):
        m.total_sigma()
    with pytest.raises(ValueError, match="no per-beamlet uncertainty"):
        m.beamlet_relative_uncertainty(0)


def test_save_load_round_trips_sigma(table, tmp_path) -> None:
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    m = assemble_influence_matrix_batched(
        eng, _beamlets(), _lat(), _dose(), 200, n_batches=4, seed=1000, path="python"
    )
    p = str(tmp_path / "influence_unc.npz")
    m.save(p)
    from ionmc.influence import SparseInfluenceMatrix

    loaded = SparseInfluenceMatrix.load(p)
    assert loaded.data_sigma is not None
    assert np.array_equal(loaded.data_sigma, m.data_sigma)
    assert np.array_equal(loaded.data, m.data)


@pytest.mark.warp
def test_per_beamlet_se_scaling(warp_module, table) -> None:
    """A representative beamlet's high-dose mean relative SEM scales as 1/sqrt(N)
    (decision 0025), inheriting the decision-0024 estimator per beamlet."""
    eng = TransportEngine(table, _box(), DepthDoseGrid(250.0, 10))
    srcs = _beamlets()
    lat, dose = _lat(), _dose()
    mN = assemble_influence_matrix_batched(
        eng, srcs, lat, dose, 8000, n_batches=16, seed=100, path="warp", device="cpu"
    )
    m4 = assemble_influence_matrix_batched(
        eng, srcs, lat, dose, 32000, n_batches=16, seed=100, path="warp", device="cpu"
    )
    uN = mN.beamlet_relative_uncertainty(1, 0.2)
    u4 = m4.beamlet_relative_uncertainty(1, 0.2)
    assert 0.35 <= u4 / uN <= 0.71
