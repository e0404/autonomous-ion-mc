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
    #: Optional per-nnz standard error of the mean dose [MeV], aligned with
    #: ``data`` (decision 0025). None for the exact single-run matrix (0022).
    data_sigma: np.ndarray | None = None

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

    def beamlet_sigma_flat(self, row: int) -> np.ndarray:
        """Dense ``(n_voxels,)`` per-voxel dose SEM [MeV] for beamlet ``row`` (0
        off-row and where no uncertainty was estimated). Requires ``data_sigma``."""
        if self.data_sigma is None:
            raise ValueError("this matrix carries no per-beamlet uncertainty")
        out = np.zeros(self.n_voxels, dtype=np.float64)
        lo, hi = int(self.indptr[row]), int(self.indptr[row + 1])
        out[self.indices[lo:hi]] = self.data_sigma[lo:hi]
        return out

    def total_sigma(self) -> np.ndarray:
        """Broad-field per-voxel dose SEM ``sqrt(Sum_b sigma_b^2)`` [MeV] as a
        ``(nx, ny, nz)`` grid: the beamlets are independent batched estimates, so
        their variances add. Requires ``data_sigma``."""
        if self.data_sigma is None:
            raise ValueError("this matrix carries no per-beamlet uncertainty")
        var = np.zeros(self.n_voxels, dtype=np.float64)
        np.add.at(var, self.indices, self.data_sigma**2)
        return np.sqrt(var).reshape(self.grid_shape)

    def beamlet_relative_uncertainty(
        self, row: int, min_dose_frac: float = 0.5
    ) -> float:
        """Mean relative SEM over the beamlet's voxels whose dose exceeds
        ``min_dose_frac`` of that beamlet's peak dose -- the per-beamlet MC quality
        metric (decision 0025). nan if no voxel clears the threshold."""
        if self.data_sigma is None:
            raise ValueError("this matrix carries no per-beamlet uncertainty")
        lo, hi = int(self.indptr[row]), int(self.indptr[row + 1])
        d = self.data[lo:hi]
        s = self.data_sigma[lo:hi]
        if d.size == 0:
            return float("nan")
        mask = d > max(0.0, min_dose_frac) * float(d.max())
        if not np.any(mask):
            return float("nan")
        return float(np.mean(s[mask] / d[mask]))

    def to_dense(self) -> np.ndarray:
        """The full ``(n_beamlets, n_voxels)`` dense matrix [MeV]."""
        dense = np.zeros((self.n_beamlets, self.n_voxels), dtype=np.float64)
        for b in range(self.n_beamlets):
            lo, hi = int(self.indptr[b]), int(self.indptr[b + 1])
            dense[b, self.indices[lo:hi]] = self.data[lo:hi]
        return dense

    def save(self, path: str) -> None:
        """Write the matrix to ``path`` as an ``.npz`` archive (including
        ``data_sigma`` when present)."""
        shape = np.asarray(self.grid_shape, dtype=np.int64)
        if self.data_sigma is not None:
            np.savez(
                path,
                grid_shape=shape,
                beamlet_ids=self.beamlet_ids,
                indptr=self.indptr,
                indices=self.indices,
                data=self.data,
                data_sigma=self.data_sigma,
            )
        else:
            np.savez(
                path,
                grid_shape=shape,
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
                data_sigma=z["data_sigma"] if "data_sigma" in z else None,
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


def assemble_influence_matrix_batched(
    engine: TransportEngine,
    sources: list[PencilBeamSource],
    grid: DepthLateralGrid,
    dose_grid: DoseGrid3D,
    n_histories: int,
    n_batches: int = 10,
    seed: int = 12345,
    path: str = "warp",
    device: str = "cpu",
    threshold_frac: float = 0.0,
) -> SparseInfluenceMatrix:
    """Assemble a sparse influence matrix carrying a per-beamlet statistical
    uncertainty (decision 0025). Each beamlet ``i`` is transported with
    ``run_scattering_batched`` over ``n_batches`` independent history batches,
    seeded from ``seed + i*n_batches`` (so per-batch seeds never collide across
    beamlets); the stored row is the per-voxel mean dose and ``data_sigma`` its
    standard error of the mean. Voxels below ``threshold_frac`` of the beamlet's
    peak *mean* dose are dropped, and ``data_sigma`` is kept for exactly the
    retained voxels so it stays aligned with ``data``."""
    if not sources:
        raise ValueError("need at least one beamlet source")
    indptr = [0]
    idx_parts: list[np.ndarray] = []
    data_parts: list[np.ndarray] = []
    sigma_parts: list[np.ndarray] = []
    beamlet_ids: list[int] = []
    for i, src in enumerate(sources):
        res = engine.run_scattering_batched(
            src,
            grid,
            dose_grid,
            n_histories,
            n_batches=n_batches,
            seed=seed + i * n_batches,
            path=path,
            device=device,
        )
        mean = res.mean_dose3d_mev.ravel()
        sem = res.standard_error_mev.ravel()
        thr = threshold_frac * float(mean.max()) if threshold_frac > 0.0 else 0.0
        nz = np.nonzero(mean > thr)[0]
        idx_parts.append(nz.astype(np.int64))
        data_parts.append(mean[nz].astype(np.float64))
        sigma_parts.append(sem[nz].astype(np.float64))
        beamlet_ids.append(int(src.beamlet))
        indptr.append(indptr[-1] + int(nz.size))
    empty = np.empty(0, dtype=np.float64)
    return SparseInfluenceMatrix(
        grid_shape=(dose_grid.nx, dose_grid.ny, dose_grid.nz),
        beamlet_ids=np.asarray(beamlet_ids, dtype=np.int64),
        indptr=np.asarray(indptr, dtype=np.int64),
        indices=np.concatenate(idx_parts) if idx_parts else np.empty(0, np.int64),
        data=np.concatenate(data_parts) if data_parts else empty,
        data_sigma=np.concatenate(sigma_parts) if sigma_parts else empty,
    )
