"""Scoring grids and tallies.

A :class:`ScoringGrid` is an axis-aligned Cartesian grid (mm) that may differ
in resolution and alignment from the transport geometry (V1-MUST-011). Energy
deposited along a transport step is distributed over the scoring voxels the
step crosses, proportionally to path length (uniform loss rate along the
step), so grid refinement or shifts do not bias the integral. Dose (Gy per
primary) = energy (MeV) / mass (g) × 1.602176634e-10 J/MeV × 1000 g/kg.

Tallies accumulate per batch: sums and sums of squares over batches give the
mean and the standard error of the mean (V1-MUST-030).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ionmc.geometry import VoxelGeometry, vec3, vec3i

MEV_TO_J = 1.602176634e-13
GY_PER_MEV_PER_G = MEV_TO_J * 1000.0  # Gy for 1 MeV deposited in 1 g


@dataclass(frozen=True)
class ScoringGrid:
    origin_mm: tuple[float, float, float]
    spacing_mm: tuple[float, float, float]
    shape: tuple[int, int, int]

    @classmethod
    def from_geometry(cls, geometry: VoxelGeometry) -> ScoringGrid:
        return cls(geometry.origin_mm, geometry.spacing_mm, geometry.shape)

    @classmethod
    def coarse(
        cls, geometry: VoxelGeometry, spacing_mm: float | tuple[float, float, float]
    ) -> ScoringGrid:
        sp = vec3(
            (spacing_mm,) * 3 if isinstance(spacing_mm, int | float) else spacing_mm
        )
        shape = vec3i(
            round(geometry.shape[a] * geometry.spacing_mm[a] / sp[a]) for a in range(3)
        )
        return cls(geometry.origin_mm, sp, shape)

    def shifted(self, delta_mm: tuple[float, float, float]) -> ScoringGrid:
        return ScoringGrid(
            vec3(o + d for o, d in zip(self.origin_mm, delta_mm, strict=True)),
            self.spacing_mm,
            self.shape,
        )

    def edges_mm(self, axis: int) -> np.ndarray:
        return self.origin_mm[axis] + self.spacing_mm[axis] * np.arange(
            self.shape[axis] + 1
        )

    def centers_mm(self, axis: int) -> np.ndarray:
        return self.origin_mm[axis] + self.spacing_mm[axis] * (
            np.arange(self.shape[axis]) + 0.5
        )

    @property
    def voxel_volume_cm3(self) -> float:
        return float(np.prod(self.spacing_mm)) * 1e-3

    def voxel_mass_g(self, geometry: VoxelGeometry) -> np.ndarray:
        """Mass of each scoring voxel from the transport density (exact box overlaps)."""
        rho = geometry.density_g_cm3.astype(np.float64)
        mats = []
        for a in range(3):
            te = geometry.edges_mm(a)
            se = self.edges_mm(a)
            lo = np.maximum(se[:-1, None], te[None, :-1])
            hi = np.minimum(se[1:, None], te[None, 1:])
            mats.append(
                np.clip(hi - lo, 0.0, None)
            )  # [n_score, n_transport] overlap length (mm)
        # sum_ijk Ox[a,i] Oy[b,j] Oz[c,k] rho[i,j,k]  (mm^3 g/cm^3 -> g via 1e-3)
        m = np.einsum("ai,ijk->ajk", mats[0], rho)
        m = np.einsum("bj,ajk->abk", mats[1], m)
        m = np.einsum("ck,abk->abc", mats[2], m)
        return m * 1e-3

    def describe(self) -> dict:
        return {
            "origin_mm": list(self.origin_mm),
            "spacing_mm": list(self.spacing_mm),
            "shape": list(self.shape),
            "units": {"dose": "Gy per primary", "energy": "MeV per primary"},
        }


@dataclass
class Tally:
    """Batch-wise accumulator of a per-voxel quantity."""

    name: str
    shape: tuple[int, ...]
    batches: int
    sum: np.ndarray = field(init=False)
    sum_sq: np.ndarray = field(init=False)
    current: np.ndarray = field(init=False)
    completed_batches: int = 0

    def __post_init__(self) -> None:
        self.sum = np.zeros(self.shape, dtype=np.float64)
        self.sum_sq = np.zeros(self.shape, dtype=np.float64)
        self.current = np.zeros(self.shape, dtype=np.float64)

    def close_batch(self, histories_in_batch: int) -> None:
        per_primary = self.current / float(histories_in_batch)
        self.sum += per_primary
        self.sum_sq += per_primary * per_primary
        self.current[...] = 0.0
        self.completed_batches += 1

    def mean(self) -> np.ndarray:
        return self.sum / max(self.completed_batches, 1)

    def standard_error(self) -> np.ndarray:
        """Standard error of the batch mean; NaN where undefined (< 2 batches)."""
        n = self.completed_batches
        if n < 2:
            return np.full(self.shape, np.nan)
        mean = self.sum / n
        var = (self.sum_sq / n - mean * mean) * n / (n - 1)
        return np.sqrt(np.maximum(var, 0.0) / n)
