"""Beam sources for Stage-1 transport (decisions ``0009``, ``0018``).

A :class:`PencilBeamSource` produces a monoenergetic proton pencil beam entering
the geometry at a fixed ``position_mm`` along a unit ``direction``. The default
direction is +z (the historical axis-aligned beam); an arbitrary direction gives
oblique/arbitrary beam incidence, transported in a canonical beam frame by the
scattering path (decision ``0018``). All histories are identical (deterministic
transport in DEV-004); the batch/beamlet machinery is present so that stochastic
sources and multiple beamlets slot in without an API change.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ionmc.rng import rand_init
from ionmc.transport.state import ParticleState, Species, Status


@dataclass(frozen=True)
class PencilBeamSource:
    """Monoenergetic proton pencil beam from ``position_mm`` along ``direction``.

    ``direction`` need not be normalised; it is normalised on use. It defaults to
    +z so existing axis-aligned callers are unchanged. The pure depth-dose path
    (:meth:`TransportEngine.run`) is longitudinal and requires +z; arbitrary
    directions are supported by the 3-D scattering path (decision ``0018``).
    """

    energy_mev: float
    position_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    direction: tuple[float, float, float] = (0.0, 0.0, 1.0)
    beamlet: int = 0

    def __post_init__(self) -> None:
        if self.energy_mev <= 0.0:
            raise ValueError("beam energy must be positive")
        norm = float(np.linalg.norm(np.asarray(self.direction, dtype=np.float64)))
        if norm <= 0.0:
            raise ValueError("beam direction must be a non-zero vector")

    @property
    def direction_hat(self) -> np.ndarray:
        """The unit beam direction."""
        d = np.asarray(self.direction, dtype=np.float64)
        return d / np.linalg.norm(d)

    def sample(self, n_histories: int, seed: int) -> ParticleState:
        """Emit ``n_histories`` identical histories with per-history RNG streams."""
        if n_histories <= 0:
            raise ValueError("n_histories must be positive")
        state = ParticleState.allocate(n_histories)
        state.position_mm[:] = np.asarray(self.position_mm, dtype=np.float64)
        state.direction[:] = self.direction_hat
        state.energy_mev[:] = self.energy_mev
        state.species[:] = int(Species.PROTON)
        state.beamlet[:] = self.beamlet
        state.status[:] = int(Status.ALIVE)
        # counter-based per-history streams (decision 0005): state = rand_init(seed, i)
        state.rng_state[:] = np.array(
            [rand_init(seed, i) for i in range(n_histories)], dtype=np.uint32
        )
        return state
