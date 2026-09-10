"""Homogeneous geometry for Stage-1 transport (decision ``0009``).

A :class:`WaterSlab` is a semi-infinite homogeneous medium entered at ``z = 0``
along +z, with a maximum depth beyond which a history has escaped. It carries
the material (for the mass density and, later, material-dependent physics). The
scoring grid is separate (:mod:`ionmc.transport.depth_dose`).
"""

from __future__ import annotations

from dataclasses import dataclass

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
