"""Sparse dose-influence matrices for treatment planning (decision ``0022``).

A dose-influence matrix has one row per beamlet and one column per dose voxel; row
``b`` is beamlet ``b``'s dose distribution. It is stored sparse (CSR-style
``indptr``/``indices``/``data``) because a pencil beam illuminates only a small
fraction of the patient. This module assembles the matrix by transporting each
beamlet with a :class:`~ionmc.transport.DoseGrid3D` scorer (decision ``0021``) and
thresholding its dose, and exports it in a documented ``.npz`` format usable by
inverse planning.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ionmc.transport import (
    DepthLateralGrid,
    DoseGrid3D,
    PencilBeamSource,
    TransportEngine,
)


@dataclass(frozen=True)
class SparseInfluenceMatrix:
    """A ``(n_beamlets, n_voxels)`` dose-influence matrix in CSR form [MeV].

    ``indptr`` has ``n_beamlets + 1`` entries; row ``b`` occupies
    ``indices[indptr[b]:indptr[b+1]]`` (flat C-order voxel indices,
    ``flat = (i*ny + j)*nz + k``) with the matching ``data`` doses. ``grid_shape``
    is the dose grid's ``(nx, ny, nz)``.
    """

    grid_shape: tuple[int, int, int]
    beamlet_ids: np.ndarray  # (n_beamlets,) int
    indptr: np.ndarray  # (n_beamlets + 1,) int
    indices: np.ndarray  # (nnz,) int flat voxel index
    data: np.ndarray  # (nnz,) float [MeV]

    @property
    def n_beamlets(self) -> int:
        return int(self.indptr.shape[0] - 1)

    @property
    def n_voxels(self) -> int:
        nx, ny, nz = self.grid_shape
        return int(nx * ny * nz)

    @property
    def nnz(self) -> int:
        return int(self.data.shape[0])

    def beamlet_dose_flat(self, row: int) -> np.ndarray:
        """Dense ``(n_voxels,)`` dose [MeV] for beamlet at ``row``."""
        out = np.zeros(self.n_voxels, dtype=np.float64)
        lo, hi = int(self.indptr[row]), int(self.indptr[row + 1])
        out[self.indices[lo:hi]] = self.data[lo:hi]
        return out

    def total_dose(self) -> np.ndarray:
        """Sum over all beamlets, as a ``(nx, ny, nz)`` dose [MeV] -- the
        broad-field dose the matrix represents."""
        tot = np.zeros(self.n_voxels, dtype=np.float64)
        np.add.at(tot, self.indices, self.data)
        return tot.reshape(self.grid_shape)

    def to_dense(self) -> np.ndarray:
        """The full ``(n_beamlets, n_voxels)`` dense matrix [MeV]."""
        dense = np.zeros((self.n_beamlets, self.n_voxels), dtype=np.float64)
        for b in range(self.n_beamlets):
            lo, hi = int(self.indptr[b]), int(self.indptr[b + 1])
            dense[b, self.indices[lo:hi]] = self.data[lo:hi]
        return dense

    def save(self, path: str) -> None:
        """Write the matrix to ``path`` as an ``.npz`` archive."""
        np.savez(
            path,
            grid_shape=np.asarray(self.grid_shape, dtype=np.int64),
            beamlet_ids=self.beamlet_ids,
            indptr=self.indptr,
            indices=self.indices,
            data=self.data,
        )

    @classmethod
    def load(cls, path: str) -> SparseInfluenceMatrix:
        """Read a matrix written by :meth:`save`."""
        with np.load(path) as z:
            shape = tuple(int(v) for v in z["grid_shape"])
            return cls(
                grid_shape=(shape[0], shape[1], shape[2]),
                beamlet_ids=z["beamlet_ids"],
                indptr=z["indptr"],
                indices=z["indices"],
                data=z["data"],
            )


def assemble_influence_matrix(
    engine: TransportEngine,
    sources: list[PencilBeamSource],
    grid: DepthLateralGrid,
    dose_grid: DoseGrid3D,
    n_histories: int = 1,
    seed: int = 12345,
    path: str = "warp",
    device: str = "cpu",
    threshold_frac: float = 0.0,
) -> SparseInfluenceMatrix:
    """Transport each beamlet with a 3-D dose scorer and assemble the sparse
    influence matrix (decision 0022). Beamlet ``i`` is seeded ``seed + i`` (the
    same partition :meth:`TransportEngine.run_scattering_multi` uses), so summing
    the beamlet rows reproduces the batched broad-field dose exactly. Voxels below
    ``threshold_frac`` of a beamlet's peak dose are dropped from that row."""
    if not sources:
        raise ValueError("need at least one beamlet source")
    indptr = [0]
    idx_parts: list[np.ndarray] = []
    data_parts: list[np.ndarray] = []
    beamlet_ids: list[int] = []
    for i, src in enumerate(sources):
        res = engine.run_scattering(
            src,
            grid,
            n_histories,
            seed=seed + i,
            path=path,
            device=device,
            dose_grid=dose_grid,
        )
        assert res.dose3d_mev is not None
        d = res.dose3d_mev.ravel()
        thr = threshold_frac * float(d.max()) if threshold_frac > 0.0 else 0.0
        nz = np.nonzero(d > thr)[0]
        idx_parts.append(nz.astype(np.int64))
        data_parts.append(d[nz].astype(np.float64))
        beamlet_ids.append(int(src.beamlet))
        indptr.append(indptr[-1] + int(nz.size))
    empty_i = np.empty(0, dtype=np.int64)
    empty_d = np.empty(0, dtype=np.float64)
    return SparseInfluenceMatrix(
        grid_shape=(dose_grid.nx, dose_grid.ny, dose_grid.nz),
        beamlet_ids=np.asarray(beamlet_ids, dtype=np.int64),
        indptr=np.asarray(indptr, dtype=np.int64),
        indices=np.concatenate(idx_parts) if idx_parts else empty_i,
        data=np.concatenate(data_parts) if data_parts else empty_d,
    )
