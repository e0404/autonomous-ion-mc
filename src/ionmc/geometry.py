"""Transport geometry: voxelized phantoms with material and density per voxel.

Conventions (V1-MUST-037/038): right-handed Cartesian coordinates in
millimetres; voxel ``(i, j, k)`` spans
``origin + (i, j, k) * spacing`` to ``origin + (i+1, j+1, k+1) * spacing`` along
``(x, y, z)``. Arrays are indexed ``[i, j, k]`` (x fastest in the physical
sense, C-order in memory with k contiguous). A homogeneous box is a voxel
geometry with one voxel per axis unless a finer transport grid is requested.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ionmc.materials import Material, get_material

Vec3 = tuple[float, float, float]


def vec3(values) -> Vec3:
    a, b, c = (float(v) for v in values)
    return (a, b, c)


def vec3i(values) -> tuple[int, int, int]:
    a, b, c = (int(v) for v in values)
    return (a, b, c)


@dataclass(frozen=True)
class VoxelGeometry:
    origin_mm: tuple[float, float, float]
    spacing_mm: tuple[float, float, float]
    material_index: np.ndarray  # int32 [nx, ny, nz], index into ``materials``
    density_g_cm3: np.ndarray  # float32 [nx, ny, nz]
    materials: tuple[Material, ...]
    name: str = "voxel-geometry"
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if (
            self.material_index.ndim != 3
            or self.density_g_cm3.shape != self.material_index.shape
        ):
            raise ValueError(
                "material_index and density_g_cm3 must be 3-D arrays of equal shape"
            )
        if self.material_index.min() < 0 or self.material_index.max() >= len(
            self.materials
        ):
            raise ValueError("material_index out of range of the materials list")
        if np.any(self.density_g_cm3 <= 0) or not np.all(
            np.isfinite(self.density_g_cm3)
        ):
            raise ValueError("densities must be positive and finite")
        if any(s <= 0 for s in self.spacing_mm):
            raise ValueError("voxel spacing must be positive")

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(n) for n in self.material_index.shape)  # type: ignore[return-value]

    @property
    def extent_mm(
        self,
    ) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
        return tuple(
            (self.origin_mm[a], self.origin_mm[a] + self.shape[a] * self.spacing_mm[a])
            for a in range(3)
        )  # type: ignore[return-value]

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

    def mass_g(self) -> np.ndarray:
        return self.density_g_cm3.astype(np.float64) * self.voxel_volume_cm3

    def describe(self) -> dict:
        return {
            "name": self.name,
            "origin_mm": list(self.origin_mm),
            "spacing_mm": list(self.spacing_mm),
            "shape": list(self.shape),
            "materials": [m.name for m in self.materials],
            "density_range_g_cm3": [
                float(self.density_g_cm3.min()),
                float(self.density_g_cm3.max()),
            ],
            "coordinate_convention": "right-handed Cartesian, mm; voxel (i,j,k) spans origin+(i,j,k)*spacing to origin+(i+1,j+1,k+1)*spacing",
            **self.metadata,
        }


def homogeneous_box(
    size_mm: tuple[float, float, float],
    spacing_mm: tuple[float, float, float] | float,
    material: Material | str = "water",
    origin_mm: tuple[float, float, float] | None = None,
    density_g_cm3: float | None = None,
) -> VoxelGeometry:
    """A box of one material. Default origin centres x and y and starts z at 0."""
    mat = get_material(material)
    if isinstance(spacing_mm, int | float):
        spacing = (float(spacing_mm),) * 3
    else:
        spacing = vec3(spacing_mm)
    shape = vec3i(round(size_mm[a] / spacing[a]) for a in range(3))
    for a in range(3):
        if abs(shape[a] * spacing[a] - size_mm[a]) > 1e-6:
            raise ValueError(
                f"size {size_mm[a]} mm is not a multiple of spacing {spacing[a]} mm on axis {a}"
            )
    if origin_mm is None:
        origin_mm = (-size_mm[0] / 2.0, -size_mm[1] / 2.0, 0.0)
    rho = mat.density_g_cm3 if density_g_cm3 is None else density_g_cm3
    return VoxelGeometry(
        vec3(origin_mm),
        spacing,
        np.zeros(shape, dtype=np.int32),
        np.full(shape, rho, dtype=np.float32),
        (mat,),
        name=f"homogeneous-{mat.name}",
    )


def slab_phantom(
    size_mm: tuple[float, float, float],
    spacing_mm: float,
    slabs: list[tuple[float, float, Material | str, float | None]],
    background: Material | str = "water",
) -> VoxelGeometry:
    """Water box with slabs along z: (z_start_mm, z_end_mm, material, density or None)."""
    geo = homogeneous_box(size_mm, spacing_mm, background)
    mats = [geo.materials[0]]
    index = geo.material_index.copy()
    rho = geo.density_g_cm3.copy()
    zc = geo.centers_mm(2)
    for z0, z1, material, density in slabs:
        mat = get_material(material)
        if mat not in mats:
            mats.append(mat)
        m = mats.index(mat)
        sel = (zc >= z0) & (zc < z1)
        index[:, :, sel] = m
        rho[:, :, sel] = mat.density_g_cm3 if density is None else density
    return VoxelGeometry(
        geo.origin_mm, geo.spacing_mm, index, rho, tuple(mats), name="slab-phantom"
    )
