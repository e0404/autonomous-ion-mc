"""Particle sources."""

from __future__ import annotations

import math
from dataclasses import dataclass

from ionmc._validate import fail, real, triple
from ionmc.physics.projectiles import Projectile


@dataclass(frozen=True)
class PencilBeamSource:
    """A pencil beam: one position, one direction, Gaussian energy and lateral spread.

    ``position_mm`` is the point of emission [mm]; the beam travels in vacuum until it enters
    the geometry. ``direction`` is any nonzero vector; the normalised value is recorded in
    the effective configuration. ``kinetic_energy_mev`` is the mean total kinetic energy of
    the ion [MeV]; ``energy_sigma_mev`` its Gaussian standard deviation [MeV];
    ``lateral_sigma_mm`` the standard deviation of the position offset in each of the two
    directions perpendicular to the beam [mm]. No axis is privileged.
    """

    projectile: Projectile
    position_mm: tuple[float, float, float]
    direction: tuple[float, float, float]
    kinetic_energy_mev: float
    energy_sigma_mev: float = 0.0
    lateral_sigma_mm: float = 0.0

    def __post_init__(self) -> None:
        if not isinstance(self.projectile, Projectile):
            raise fail(f"projectile must be a Projectile, got {self.projectile!r}")
        object.__setattr__(self, "position_mm", triple("position_mm", self.position_mm))
        d = triple("direction", self.direction)
        if math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) <= 1e-12:
            raise fail("direction must be a nonzero vector")
        object.__setattr__(self, "direction", d)
        real("kinetic_energy_mev", self.kinetic_energy_mev, positive=True)
        real("energy_sigma_mev", self.energy_sigma_mev, nonnegative=True)
        real("lateral_sigma_mm", self.lateral_sigma_mm, nonnegative=True)

    def unit_direction(self) -> tuple[float, float, float]:
        """The normalised direction (float64)."""
        d = self.direction
        n = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        return (d[0] / n, d[1] / n, d[2] / n)
