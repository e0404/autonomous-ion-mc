"""Transport geometry: a voxel grid of materials (and optional densities) or a box phantom.

Conventions (decision 0037): right-handed Cartesian world coordinates in mm; ``origin_mm`` is
the corner of voxel ``(0, 0, 0)``; voxel ``(i, j, k)`` spans ``origin + (i, j, k) * spacing``
to ``origin + (i + 1, j + 1, k + 1) * spacing``; arrays are C-ordered ``[ix, iy, iz]``.
The world is the union of the voxels; everything outside is vacuum (particles that leave do
not return). Vacuum *voxels* (density <= 0) are not supported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ionmc._frozen import freeze_array
from ionmc._validate import fail, int_triple, real, triple
from ionmc.materials import Material


@dataclass(frozen=True, eq=False)
class VoxelGeometry:
    """Regular voxel grid of materials.

    ``material_index`` is an integer array of shape ``shape`` indexing ``materials``.
    ``density_g_cm3`` (optional, same shape) overrides the nominal density of each voxel's
    material [g/cm3]; it must be finite and positive everywhere.
    """

    origin_mm: tuple[float, float, float]
    spacing_mm: tuple[float, float, float]
    shape: tuple[int, int, int]
    materials: tuple[Material, ...]
    material_index: NDArray[np.int32]
    density_g_cm3: NDArray[np.float64] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "origin_mm", triple("origin_mm", self.origin_mm))
        object.__setattr__(self, "spacing_mm", triple("spacing_mm", self.spacing_mm, positive=True))
        object.__setattr__(self, "shape", int_triple("shape", self.shape))
        if not isinstance(self.materials, tuple) or not self.materials:
            raise fail("materials must be a non-empty tuple of Material")
        for m in self.materials:
            if not isinstance(m, Material):
                raise fail(f"materials must contain Material objects, got {m!r}")
        raw = np.asarray(self.material_index)
        if raw.shape != self.shape:
            raise fail(f"material_index shape {raw.shape} does not match shape {self.shape}")
        if raw.dtype.kind not in "iu":
            raise fail("material_index must have an integer dtype")
        if raw.size and (raw.min() < 0 or raw.max() >= len(self.materials)):
            raise fail(
                f"material_index values must lie in [0, {len(self.materials)}), "
                f"found [{raw.min()}, {raw.max()}]"
            )
        object.__setattr__(self, "material_index", freeze_array(raw, np.int32, "material_index"))
        if self.density_g_cm3 is not None:
            d = np.array(self.density_g_cm3, dtype=np.float64)
            if d.shape != self.shape:
                raise fail(f"density_g_cm3 shape {d.shape} does not match shape {self.shape}")
            if not np.all(np.isfinite(d)):
                raise fail("density_g_cm3 must be finite")
            if np.any(d <= 0.0):
                raise fail("density_g_cm3 must be > 0 everywhere (vacuum voxels are unsupported)")
            object.__setattr__(self, "density_g_cm3", freeze_array(d, np.float64, "density_g_cm3"))

    @property
    def lower_mm(self) -> tuple[float, float, float]:
        """Lower corner of the world box [mm]."""
        return self.origin_mm

    @property
    def upper_mm(self) -> tuple[float, float, float]:
        """Upper corner of the world box [mm]."""
        o, s, n = self.origin_mm, self.spacing_mm, self.shape
        return (o[0] + s[0] * n[0], o[1] + s[1] * n[1], o[2] + s[2] * n[2])

    @property
    def n_voxels(self) -> int:
        """Number of voxels."""
        return self.shape[0] * self.shape[1] * self.shape[2]

    def densities_g_cm3(self) -> NDArray[np.float64]:
        """Density of every voxel [g/cm3], shape ``shape`` (nominal where not overridden)."""
        if self.density_g_cm3 is not None:
            return self.density_g_cm3
        nominal = np.array([m.density_g_cm3 for m in self.materials], dtype=np.float64)
        out: NDArray[np.float64] = nominal[self.material_index]
        return out


@dataclass(frozen=True)
class BoxPhantom:
    """A homogeneous box: the world is the box, vacuum outside; a 1x1x1 voxel geometry."""

    lower_mm: tuple[float, float, float]
    size_mm: tuple[float, float, float]
    material: Material

    def __post_init__(self) -> None:
        object.__setattr__(self, "lower_mm", triple("lower_mm", self.lower_mm))
        object.__setattr__(self, "size_mm", triple("size_mm", self.size_mm, positive=True))
        if not isinstance(self.material, Material):
            raise fail(f"material must be a Material, got {self.material!r}")
        real("density", self.material.density_g_cm3, positive=True)

    def to_geometry(self) -> VoxelGeometry:
        """The equivalent one-voxel :class:`VoxelGeometry`."""
        return VoxelGeometry(
            origin_mm=self.lower_mm,
            spacing_mm=self.size_mm,
            shape=(1, 1, 1),
            materials=(self.material,),
            material_index=np.zeros((1, 1, 1), dtype=np.int32),
        )
