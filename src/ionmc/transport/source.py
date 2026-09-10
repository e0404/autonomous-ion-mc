"""Beam sources for Stage-1 transport (decision ``0009``).

A :class:`PencilBeamSource` produces a monoenergetic proton pencil beam
entering a slab at ``z = 0`` along +z from a fixed lateral position. All
histories are identical (deterministic transport in DEV-004); the batch/beamlet
machinery is present so that stochastic sources and multiple beamlets slot in
without an API change.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ionmc.rng import rand_init
from ionmc.transport.state import ParticleState, Species, Status


@dataclass(frozen=True)
class PencilBeamSource:
    """Monoenergetic proton pencil beam along +z."""

    energy_mev: float
    position_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    beamlet: int = 0

    def __post_init__(self) -> None:
        if self.energy_mev <= 0.0:
            raise ValueError("beam energy must be positive")

    def sample(self, n_histories: int, seed: int) -> ParticleState:
        """Emit ``n_histories`` identical histories with per-history RNG streams."""
        if n_histories <= 0:
            raise ValueError("n_histories must be positive")
        state = ParticleState.allocate(n_histories)
        state.position_mm[:] = np.asarray(self.position_mm, dtype=np.float64)
        state.direction[:] = np.array([0.0, 0.0, 1.0])
        state.energy_mev[:] = self.energy_mev
        state.species[:] = int(Species.PROTON)
        state.beamlet[:] = self.beamlet
        state.status[:] = int(Status.ALIVE)
        # counter-based per-history streams (decision 0005): state = rand_init(seed, i)
        state.rng_state[:] = np.array(
            [rand_init(seed, i) for i in range(n_histories)], dtype=np.uint32
        )
        return state
