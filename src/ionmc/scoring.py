"""Scoring grids, voxel masses and the batch estimator (float64 numpy).

A scoring grid is an independent regular grid (same conventions as
:class:`ionmc.geometry.VoxelGeometry`; it need not align with the transport grid). Energy is
deposited at the midpoint of each step into the voxel containing it. The mass of a scoring
voxel is the exact overlap integral of the piecewise-constant density of the geometry over
the voxel (the geometry is vacuum outside, so partially covered voxels have a correspondingly
smaller mass). Dose uses 1 MeV/g = 1.602176634e-10 Gy.

Batch estimator: histories are assigned to batches by ``history_index mod n_batches``; with
per-batch energy sums ``E_b`` [MeV] and ``N/B`` histories per batch, ``x_b = E_b / (N/B)``
is the per-primary estimate of batch ``b``; the mean is ``mean(x_b)`` and the variance of
the mean ``sum (x_b - mean)^2 / (B (B - 1))``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ionmc._validate import fail, int_triple, triple
from ionmc.geometry import VoxelGeometry

MEV_PER_G_TO_GY = 1.602176634e-10
MAX_SCORING_GRIDS = 4


@dataclass(frozen=True)
class ScoringGrid:
    """A regular scoring grid: ``origin_mm`` (corner of voxel (0,0,0)), ``spacing_mm`` (> 0),
    ``shape`` ``(nx, ny, nz)``, C-ordered ``[ix, iy, iz]``; ``name`` labels the result."""

    origin_mm: tuple[float, float, float]
    spacing_mm: tuple[float, float, float]
    shape: tuple[int, int, int]
    name: str = "dose"

    def __post_init__(self) -> None:
        object.__setattr__(self, "origin_mm", triple("origin_mm", self.origin_mm))
        object.__setattr__(self, "spacing_mm", triple("spacing_mm", self.spacing_mm, positive=True))
        object.__setattr__(self, "shape", int_triple("shape", self.shape))
        if not isinstance(self.name, str) or not self.name:
            raise fail("name must be a non-empty string")

    @property
    def n_voxels(self) -> int:
        """Number of voxels."""
        return self.shape[0] * self.shape[1] * self.shape[2]

    def edges_mm(self, axis: int) -> NDArray[np.float64]:
        """Voxel boundary coordinates along ``axis`` [mm], ``shape[axis] + 1`` values."""
        return self.origin_mm[axis] + self.spacing_mm[axis] * np.arange(
            self.shape[axis] + 1, dtype=np.float64
        )


def overlap_matrix(
    edges_a: NDArray[np.float64], edges_b: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Overlap length [mm] ``O[i, j]`` of interval ``i`` of ``edges_a`` and ``j`` of ``edges_b``."""
    lo = np.maximum(edges_a[:-1, None], edges_b[None, :-1])
    hi = np.minimum(edges_a[1:, None], edges_b[None, 1:])
    out: NDArray[np.float64] = np.clip(hi - lo, 0.0, None)
    return out


def voxel_mass_g(grid: ScoringGrid, geometry: VoxelGeometry) -> NDArray[np.float64]:
    """Mass [g] of every scoring voxel, shape ``grid.shape``: the exact overlap integral of the
    geometry density [g/cm3] over the voxel (volumes converted from mm3 to cm3)."""
    ox = overlap_matrix(grid.edges_mm(0), _geometry_edges(geometry, 0))
    oy = overlap_matrix(grid.edges_mm(1), _geometry_edges(geometry, 1))
    oz = overlap_matrix(grid.edges_mm(2), _geometry_edges(geometry, 2))
    rho = geometry.densities_g_cm3()
    mass = np.einsum("ia,jb,kc,abc->ijk", ox, oy, oz, rho, optimize=True)
    out: NDArray[np.float64] = mass * 1.0e-3
    return out


def _geometry_edges(geometry: VoxelGeometry, axis: int) -> NDArray[np.float64]:
    return geometry.origin_mm[axis] + geometry.spacing_mm[axis] * np.arange(
        geometry.shape[axis] + 1, dtype=np.float64
    )


@dataclass(frozen=True, eq=False)
class BatchStatistics:
    """Batch-method estimate per voxel: ``mean`` and ``variance_of_mean`` (per primary,
    same units as the input divided by histories per batch), ``n_nonzero`` batches."""

    mean: NDArray[np.float64]
    variance_of_mean: NDArray[np.float64]
    n_nonzero: NDArray[np.int32]


def reduce_batches(batch_sums: NDArray[np.float64], histories_per_batch: int) -> BatchStatistics:
    """Reduce per-batch sums ``[B, n_voxels]`` (float64) to per-primary mean and variance."""
    s = np.asarray(batch_sums, dtype=np.float64)
    if s.ndim != 2 or s.shape[0] < 2:
        raise ValueError("batch sums must have shape [B >= 2, n_voxels]")
    if histories_per_batch < 1:
        raise ValueError("histories_per_batch must be >= 1")
    if not np.all(np.isfinite(s)):
        raise ValueError("batch sums must be finite")
    b = s.shape[0]
    x = s / float(histories_per_batch)
    mean = x.mean(axis=0)
    var = ((x - mean) ** 2).sum(axis=0) / (b * (b - 1))
    return BatchStatistics(mean, var, (s > 0.0).sum(axis=0).astype(np.int32))
