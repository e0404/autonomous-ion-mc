"""Particle sources."""

from __future__ import annotations

import math
from dataclasses import dataclass

from ionmc._validate import fail, real, triple
from ionmc.physics.projectiles import Projectile

MAX_ANGULAR_SIGMA_MRAD = 100.0
"""Small-angle limit of the angular spread (fail closed above)."""


@dataclass(frozen=True)
class PencilBeamSource:
    """A pencil beam: one position, one direction, Gaussian energy and lateral spread.

    ``position_mm`` is the point of emission [mm]; the beam travels in vacuum until it enters
    the geometry. ``direction`` is any nonzero vector; the normalised value is recorded in
    the effective configuration. ``kinetic_energy_mev`` is the mean total kinetic energy of
    the ion [MeV]; ``energy_sigma_mev`` its Gaussian standard deviation [MeV];
    ``lateral_sigma_mm`` the standard deviation of the position offset in each of the two
    directions perpendicular to the beam [mm]. No axis is privileged.

    ``angular_sigma_mrad`` (default 0: no extra random number is drawn, every existing
    configuration is bit-identical) is the standard deviation of the projected angle in each
    transverse plane [mrad]; ``x_theta_correlation`` the coefficient rho of the bivariate Gaussian
    (x, theta) per plane: ``theta = sigma_theta (rho g1 + sqrt(1 - rho^2) g2)`` with ``g1`` the
    normal of the plane's lateral offset (purpose 1 block 0) and ``g2`` a normal of purpose 1
    block 1. The direction is ``normalise(d + tan(theta_1) e1 + tan(theta_2) e2)``, i.e.
    (tan thx, tan thy, 1) normalised for ``d = z``; positive theta tilts toward +e.
    """

    projectile: Projectile
    position_mm: tuple[float, float, float]
    direction: tuple[float, float, float]
    kinetic_energy_mev: float
    energy_sigma_mev: float = 0.0
    lateral_sigma_mm: float = 0.0
    angular_sigma_mrad: float = 0.0
    x_theta_correlation: float = 0.0

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
        real("angular_sigma_mrad", self.angular_sigma_mrad, nonnegative=True)
        real("x_theta_correlation", self.x_theta_correlation)
        if self.angular_sigma_mrad > MAX_ANGULAR_SIGMA_MRAD:
            raise fail(
                f"angular_sigma_mrad ({self.angular_sigma_mrad}) exceeds the small-angle limit "
                f"{MAX_ANGULAR_SIGMA_MRAD} mrad"
            )
        if not -1.0 < self.x_theta_correlation < 1.0:
            raise fail("x_theta_correlation must satisfy |rho| < 1")
        if self.x_theta_correlation != 0.0 and not (
            self.angular_sigma_mrad > 0.0 and self.lateral_sigma_mm > 0.0
        ):
            raise fail(
                "a nonzero x_theta_correlation needs angular_sigma_mrad > 0 and "
                "lateral_sigma_mm > 0"
            )

    def unit_direction(self) -> tuple[float, float, float]:
        """The normalised direction (float64)."""
        d = self.direction
        n = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        return (d[0] / n, d[1] / n, d[2] / n)
