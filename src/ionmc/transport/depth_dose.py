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
    """Uniform depth bins over ``[origin_mm, origin_mm + depth_mm]``.

    ``origin_mm`` (default 0) shifts the grid along the beam axis independently
    of the transport geometry, so the scoring grid's resolution *and* alignment
    are decoupled from the transport grid (decision ``0019``); the default
    reproduces the historical grid anchored at ``z = 0``.
    """

    depth_mm: float
    n_bins: int
    origin_mm: float = 0.0

    def __post_init__(self) -> None:
        if self.depth_mm <= 0.0 or self.n_bins <= 0:
            raise ValueError("depth and bin count must be positive")

    @property
    def bin_width_mm(self) -> float:
        return self.depth_mm / self.n_bins

    @property
    def edges_mm(self) -> np.ndarray:
        return np.linspace(
            self.origin_mm, self.origin_mm + self.depth_mm, self.n_bins + 1
        )

    @property
    def centers_mm(self) -> np.ndarray:
        e = self.edges_mm
        return 0.5 * (e[:-1] + e[1:])

    def empty(self) -> np.ndarray:
        """A zeroed energy-deposition array (MeV per bin)."""
        return np.zeros(self.n_bins, dtype=np.float64)


@dataclass(frozen=True)
class DepthLateralGrid:
    """2-D scoring grid: depth ``z`` in ``[depth_origin_mm, depth_origin_mm +
    depth_mm]`` x lateral ``x`` in ``[lateral_center_mm +/- half_width_mm]``
    (decisions 0011, 0019).

    Energy is accumulated per (depth, lateral-x) bin, summed over the second
    transverse axis ``y`` (a marginal projection). The depth marginal (sum over
    ``x``) is the integral depth dose; the second moment in ``x`` per depth
    slice gives the lateral spread ``sigma_x(z)``. ``depth_origin_mm`` and
    ``lateral_center_mm`` (both default 0) shift the grid independently of the
    transport geometry, decoupling scoring-grid alignment from the transport
    grid; the defaults reproduce the historical grid anchored at ``z = 0`` and
    centred at ``x = 0``.
    """

    depth_mm: float
    n_depth: int
    half_width_mm: float
    n_lateral: int
    depth_origin_mm: float = 0.0
    lateral_center_mm: float = 0.0

    def __post_init__(self) -> None:
        if self.depth_mm <= 0.0 or self.n_depth <= 0:
            raise ValueError("depth and depth-bin count must be positive")
        if self.half_width_mm <= 0.0 or self.n_lateral <= 0:
            raise ValueError("half width and lateral-bin count must be positive")

    @property
    def depth_bin_mm(self) -> float:
        return self.depth_mm / self.n_depth

    @property
    def lateral_bin_mm(self) -> float:
        return 2.0 * self.half_width_mm / self.n_lateral

    @property
    def lateral_lo_mm(self) -> float:
        """Lower lateral edge ``lateral_center_mm - half_width_mm``."""
        return self.lateral_center_mm - self.half_width_mm

    @property
    def depth_centers_mm(self) -> np.ndarray:
        edges = np.linspace(
            self.depth_origin_mm, self.depth_origin_mm + self.depth_mm, self.n_depth + 1
        )
        return 0.5 * (edges[:-1] + edges[1:])

    @property
    def lateral_centers_mm(self) -> np.ndarray:
        edges = np.linspace(
            self.lateral_center_mm - self.half_width_mm,
            self.lateral_center_mm + self.half_width_mm,
            self.n_lateral + 1,
        )
        return 0.5 * (edges[:-1] + edges[1:])

    def empty(self) -> np.ndarray:
        """A zeroed ``(n_depth, n_lateral)`` energy array [MeV]."""
        return np.zeros((self.n_depth, self.n_lateral), dtype=np.float64)

    def depth_dose(self, edep_zx: np.ndarray) -> np.ndarray:
        """Integral depth dose [MeV per depth bin] = sum over the lateral axis."""
        return edep_zx.sum(axis=1)

    def sigma_x_mm(self, edep_zx: np.ndarray) -> np.ndarray:
        """Lateral RMS ``sigma_x`` [mm] per depth slice (nan where no dose).

        ``sigma_x^2 = sum_x x^2 E / sum_x E - (sum_x x E / sum_x E)^2`` using the
        lateral bin centres; symmetric beams have mean ~0.
        """
        x = self.lateral_centers_mm
        w = edep_zx
        total = w.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = (w * x).sum(axis=1) / total
            mean_sq = (w * x * x).sum(axis=1) / total
            var = mean_sq - mean * mean
        sigma = np.sqrt(np.clip(var, 0.0, None))
        sigma[total <= 0.0] = np.nan
        return sigma
