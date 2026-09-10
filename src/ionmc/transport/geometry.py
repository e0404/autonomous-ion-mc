"""Geometry for proton depth-dose transport (decisions ``0009``, ``0014``).

A :class:`WaterSlab` is a semi-infinite homogeneous medium entered at ``z = 0``
along +z, with a maximum depth beyond which a history has escaped. A
:class:`VoxelSlab` (decision ``0014``) generalises it to a 1-D stack of voxels
along the beam axis, each with its own mass density, so the transport becomes
density-heterogeneous (water-equivalent-thickness behaviour, layered phantoms).
Both expose :meth:`voxel_profile`, which the engine uses to drive per-step
density; a homogeneous slab is a single voxel spanning the depth, so the
heterogeneous path reproduces the homogeneous one exactly.

The scoring grid is separate (:mod:`ionmc.transport.depth_dose`).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ionmc.materials import WATER, Material


@dataclass(frozen=True)
class WaterSlab:
    """A homogeneous slab of ``material`` from ``z = 0`` to ``z = depth_mm``."""

    depth_mm: float = 400.0
    material: Material = WATER

    def __post_init__(self) -> None:
        if self.depth_mm <= 0.0:
            raise ValueError("slab depth must be positive")

    @property
    def density_g_per_cm3(self) -> float:
        return self.material.density_g_per_cm3

    def voxel_profile(self) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(z_boundaries_mm, density_g_per_cm3)`` as a single voxel.

        The boundaries are ``[0, depth_mm]`` and the density is the slab density,
        so the voxelized transport path reproduces the homogeneous one exactly.
        """
        return (
            np.array([0.0, self.depth_mm], dtype=np.float64),
            np.array([self.density_g_per_cm3], dtype=np.float64),
        )


@dataclass(frozen=True)
class VoxelSlab:
    """A 1-D stack of voxels along +z, each with its own mass density.

    Construct from explicit boundary/density arrays, or with :meth:`from_layers`
    from a list of ``(thickness_mm, density_g_per_cm3)`` layers (exact interface
    positions), or with :meth:`uniform` for a single-density slab discretised
    into ``n`` voxels. All voxels share ``material`` (its composition is used for
    the density-scaled physics; per-voxel *composition* is a later task).
    """

    z_boundaries_mm: np.ndarray
    density_g_per_cm3: np.ndarray
    material: Material = WATER
    _z: np.ndarray = field(init=False, repr=False, compare=False)
    _rho: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        z = np.ascontiguousarray(self.z_boundaries_mm, dtype=np.float64)
        rho = np.ascontiguousarray(self.density_g_per_cm3, dtype=np.float64)
        if z.ndim != 1 or rho.ndim != 1:
            raise ValueError("boundaries and densities must be 1-D")
        if z.shape[0] != rho.shape[0] + 1:
            raise ValueError("z_boundaries must have one more entry than densities")
        if z.shape[0] < 2:
            raise ValueError("need at least one voxel")
        if not np.all(np.diff(z) > 0.0):
            raise ValueError("z_boundaries must be strictly increasing")
        if z[0] != 0.0:
            raise ValueError("the first boundary must be z = 0")
        if not np.all(rho > 0.0):
            raise ValueError("all voxel densities must be positive")
        object.__setattr__(self, "_z", z)
        object.__setattr__(self, "_rho", rho)

    @classmethod
    def from_layers(
        cls, layers: list[tuple[float, float]], material: Material = WATER
    ) -> VoxelSlab:
        """Build from ``(thickness_mm, density_g_per_cm3)`` layers (front to back)."""
        if not layers:
            raise ValueError("need at least one layer")
        thicknesses = [t for t, _ in layers]
        densities = [d for _, d in layers]
        if any(t <= 0.0 for t in thicknesses):
            raise ValueError("layer thicknesses must be positive")
        boundaries = np.concatenate([[0.0], np.cumsum(thicknesses)])
        return cls(
            z_boundaries_mm=boundaries,
            density_g_per_cm3=np.asarray(densities, dtype=np.float64),
            material=material,
        )

    @classmethod
    def uniform(
        cls,
        depth_mm: float,
        density: float,
        n_voxels: int,
        material: Material = WATER,
    ) -> VoxelSlab:
        """A single-density slab discretised into ``n_voxels`` equal voxels."""
        if depth_mm <= 0.0 or n_voxels < 1:
            raise ValueError("depth must be positive and n_voxels >= 1")
        boundaries = np.linspace(0.0, depth_mm, n_voxels + 1)
        densities = np.full(n_voxels, float(density), dtype=np.float64)
        return cls(
            z_boundaries_mm=boundaries,
            density_g_per_cm3=densities,
            material=material,
        )

    @property
    def depth_mm(self) -> float:
        return float(self._z[-1])

    @property
    def n_voxels(self) -> int:
        return int(self._rho.shape[0])

    def water_equivalent_depth_mm(self, z_mm: float) -> float:
        """Integrated water-equivalent thickness ``integral rho dl`` from 0 to
        ``z_mm`` [mm] (WET; the depth at which the peak lands scales with this)."""
        wet = 0.0
        z = self._z
        rho = self._rho
        for v in range(self.n_voxels):
            lo = z[v]
            hi = z[v + 1]
            if z_mm <= lo:
                break
            seg = min(z_mm, hi) - lo
            wet += rho[v] * seg
        return float(wet)

    def voxel_profile(self) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(z_boundaries_mm, density_g_per_cm3)``."""
        return self._z, self._rho
