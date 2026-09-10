"""Integral depth-dose scoring grid (decision ``0009``).

A uniform grid of ``n_bins`` depth bins over ``[0, depth_mm]`` accumulating the
energy deposited per bin (MeV). For an on-axis pencil beam with no lateral
transport this integral depth dose (energy per depth interval) is the
physically meaningful 1-D quantity; volumetric dose needs a lateral extent
introduced with multiple scattering (a later task). The grid is decoupled from
the transport geometry, as ``REQUIREMENTS.md`` requires.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DepthDoseGrid:
    """Uniform depth bins over ``[0, depth_mm]``."""

    depth_mm: float
    n_bins: int

    def __post_init__(self) -> None:
        if self.depth_mm <= 0.0 or self.n_bins <= 0:
            raise ValueError("depth and bin count must be positive")

    @property
    def bin_width_mm(self) -> float:
        return self.depth_mm / self.n_bins

    @property
    def edges_mm(self) -> np.ndarray:
        return np.linspace(0.0, self.depth_mm, self.n_bins + 1)

    @property
    def centers_mm(self) -> np.ndarray:
        e = self.edges_mm
        return 0.5 * (e[:-1] + e[1:])

    def empty(self) -> np.ndarray:
        """A zeroed energy-deposition array (MeV per bin)."""
        return np.zeros(self.n_bins, dtype=np.float64)
