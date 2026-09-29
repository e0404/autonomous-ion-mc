"""Particle sources.

A :class:`PencilBeam` emits one species from a point or Gaussian spot with a
given direction (unit vector, any orientation), optional Gaussian angular
divergence and optional Gaussian energy spread. Energies are kinetic energy
per nucleon (MeV/u); for protons this equals the total kinetic energy.
Sampling is done in NumPy (reference path) or in kernels from the same
parameters; both use the same definitions below.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ionmc.species import Species, get_species


def _unit(v: tuple[float, float, float]) -> tuple[float, float, float]:
    a = np.asarray(v, dtype=np.float64)
    n = float(np.linalg.norm(a))
    if n == 0.0:
        raise ValueError("direction must be non-zero")
    return tuple(float(x) for x in a / n)  # type: ignore[return-value]


def orthonormal_frame(
    direction: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (u, v, w) with w = direction; u, v span the transverse plane."""
    w = np.asarray(_unit(direction))
    helper = np.array([1.0, 0.0, 0.0]) if abs(w[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(w, helper)
    u /= np.linalg.norm(u)
    v = np.cross(w, u)
    return u, v, w


@dataclass(frozen=True)
class PencilBeam:
    species: Species | str
    energy_mev_per_u: float
    position_mm: tuple[float, float, float] = (0.0, 0.0, -1.0)
    direction: tuple[float, float, float] = (0.0, 0.0, 1.0)
    sigma_mm: float = 0.0  # Gaussian spot sigma (isotropic in the transverse plane)
    sigma_energy_fraction: float = 0.0  # Gaussian relative energy spread (sigma/E)
    sigma_angle_rad: float = 0.0  # Gaussian angular divergence per transverse axis
    name: str = "pencil-beam"

    def __post_init__(self) -> None:
        object.__setattr__(self, "species", get_species(self.species))
        object.__setattr__(self, "direction", _unit(self.direction))
        if self.energy_mev_per_u <= 0:
            raise ValueError("energy must be positive")
        if min(self.sigma_mm, self.sigma_energy_fraction, self.sigma_angle_rad) < 0:
            raise ValueError("spread parameters must be non-negative")

    def sample(self, n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
        """Sample n primaries: positions (mm), unit directions, energies (MeV/u)."""
        u, v, w = orthonormal_frame(self.direction)
        pos = np.tile(np.asarray(self.position_mm, dtype=np.float64), (n, 1))
        if self.sigma_mm > 0:
            a = rng.normal(0.0, self.sigma_mm, n)
            b = rng.normal(0.0, self.sigma_mm, n)
            pos += a[:, None] * u + b[:, None] * v
        dirs = np.tile(w, (n, 1))
        if self.sigma_angle_rad > 0:
            ta = rng.normal(0.0, self.sigma_angle_rad, n)
            tb = rng.normal(0.0, self.sigma_angle_rad, n)
            dirs = dirs + ta[:, None] * u + tb[:, None] * v
            dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        e = np.full(n, self.energy_mev_per_u, dtype=np.float64)
        if self.sigma_energy_fraction > 0:
            e = e * (1.0 + rng.normal(0.0, self.sigma_energy_fraction, n))
            e = np.maximum(e, 1e-3)
        return {"position_mm": pos, "direction": dirs, "energy_mev_per_u": e}

    def describe(self) -> dict:
        sp = get_species(self.species)
        return {
            "type": "pencil-beam",
            "name": self.name,
            "species": sp.name,
            "energy_mev_per_u": self.energy_mev_per_u,
            "energy_total_mev": self.energy_mev_per_u * sp.a,
            "position_mm": list(self.position_mm),
            "direction": list(self.direction),
            "sigma_mm": self.sigma_mm,
            "sigma_energy_fraction": self.sigma_energy_fraction,
            "sigma_angle_rad": self.sigma_angle_rad,
        }
