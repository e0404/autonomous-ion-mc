"""Particle-state layout for proton transport (decision ``0009``).

Structure-of-arrays over the history batch. Species and beamlet identity and a
per-history RNG state are present from the first version so that secondaries
(Stage 2) and beamlet-resolved influence matrices (Stage 4) do not force a
state-layout rewrite. Positions are in mm, energy in MeV; arrays are float64
here (the reference layout) and copied to float32 for the Warp device.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np


class Species(IntEnum):
    """Transported particle species (integer ids carried in the state)."""

    PROTON = 0


class Status(IntEnum):
    """Per-history status."""

    ALIVE = 0
    STOPPED = 1  # reached the energy cutoff
    ESCAPED = 2  # left the geometry
    REACTED = 4  # removed by a nonelastic nuclear reaction (3 marks truncation)


@dataclass
class ParticleState:
    """Structure-of-arrays batch of ``n`` histories.

    Attributes are ``(n, 3)`` for vectors and ``(n,)`` for scalars.
    """

    position_mm: np.ndarray
    direction: np.ndarray
    energy_mev: np.ndarray
    weight: np.ndarray
    species: np.ndarray
    beamlet: np.ndarray
    rng_state: np.ndarray
    status: np.ndarray

    @property
    def size(self) -> int:
        return int(self.energy_mev.shape[0])

    @classmethod
    def allocate(cls, n: int) -> ParticleState:
        """Allocate an ``n``-history state (all fields zeroed / proton / alive)."""
        return cls(
            position_mm=np.zeros((n, 3), dtype=np.float64),
            direction=np.tile(np.array([0.0, 0.0, 1.0]), (n, 1)),
            energy_mev=np.zeros(n, dtype=np.float64),
            weight=np.ones(n, dtype=np.float64),
            species=np.full(n, int(Species.PROTON), dtype=np.int32),
            beamlet=np.zeros(n, dtype=np.int32),
            rng_state=np.zeros(n, dtype=np.uint32),
            status=np.full(n, int(Status.ALIVE), dtype=np.int32),
        )
