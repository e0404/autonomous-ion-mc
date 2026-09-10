"""Parsing and preparation of tabulated stopping-power files.

:func:`parse_two_column_table` reads the MCsquare ``*_Stop_Pow.dat`` layout
(two whitespace-separated columns: kinetic energy in MeV, mass stopping power
in MeV cm^2/g, one row per line, no header) and :func:`prepare_stopping_table`
turns it into the :class:`StoppingTable` used by the tabulated physics layer
(decision ``0008``):

* rows with non-positive energy are dropped (the MCsquare files carry a
  non-physical ``0 MeV`` row), energies must be strictly increasing and the
  values positive and finite;
* the **monotone cubic Hermite (PCHIP) slopes** of Fritsch & Carlson (SIAM J.
  Numer. Anal. 17 (1980) 235) with the end-point rule of Moler
  (*Numerical Computing with MATLAB*), which is what MATLAB ``pchip`` and
  SciPy ``PchipInterpolator`` implement, are computed at every grid point;
* the **cumulative CSDA range** from the first grid point to every grid
  point is obtained by 4-point Gauss-Legendre quadrature of ``1/S`` over the
  Hermite interpolant of each segment, the same rule the shared-source
  partial-segment integral uses, so the two agree at grid points to
  rounding (relative quadrature error below 2e-7 on the steepest 0.5 MeV
  segment of the water table and far smaller elsewhere).

Preparation is deterministic and depends only on the file contents, so it is
repeated on every load rather than cached.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

#: Gauss-Legendre points used for the cumulative range (matches the kernel).
RANGE_QUADRATURE_POINTS = 4


class TableFormatError(ValueError):
    """The file does not have the expected two-column numeric layout."""


def parse_two_column_table(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Read ``energy  value`` rows; returns two float64 arrays in file order."""
    energies: list[float] = []
    values: list[float] = []
    with open(path, encoding="utf-8") as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise TableFormatError(
                    f"{path}:{lineno}: expected 2 columns, got {len(parts)}"
                )
            try:
                energies.append(float(parts[0]))
                values.append(float(parts[1]))
            except ValueError as exc:
                raise TableFormatError(f"{path}:{lineno}: non-numeric entry") from exc
    if not energies:
        raise TableFormatError(f"{path}: no data rows")
    return np.asarray(energies, dtype=np.float64), np.asarray(values, dtype=np.float64)


@dataclass(frozen=True)
class StoppingTable:
    """A prepared stopping-power table for one material and projectile.

    Attributes
    ----------
    energy_mev, stopping_mev_cm2_per_g:
        The tabulated grid (strictly increasing, positive) and values.
    slope:
        PCHIP derivative ``dS/dE`` [MeV cm^2/g per MeV] at every grid point.
    csda_range_g_per_cm2:
        Cumulative CSDA range from the first grid point to each grid point.
        The range below the first grid point is *not* included
        (``range_floor_energy_mev``); decision ``0008`` quantifies it.
    provenance:
        Free-form metadata (dataset spec, checksum, retrieval time ...).
    """

    energy_mev: np.ndarray
    stopping_mev_cm2_per_g: np.ndarray
    slope: np.ndarray
    csda_range_g_per_cm2: np.ndarray
    provenance: dict[str, Any]

    @property
    def size(self) -> int:
        return int(self.energy_mev.shape[0])

    @property
    def range_floor_energy_mev(self) -> float:
        return float(self.energy_mev[0])

    @property
    def bisection_steps(self) -> int:
        """Iterations of the fixed-count bisection that locate any energy."""
        return max(1, math.ceil(math.log2(self.size)))


def pchip_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Fritsch-Carlson monotone cubic Hermite slopes with Moler's end rule.

    Interior points use the weighted harmonic mean of the neighbouring secant
    slopes (zero where the secants change sign or vanish); the end points use
    the shape-preserving three-point formula.
    """
    h = np.diff(x)
    delta = np.diff(y) / h
    d = np.zeros_like(y)
    if y.shape[0] == 2:
        d[:] = delta[0]
        return d
    # interior
    d0, d1 = delta[:-1], delta[1:]
    h0, h1 = h[:-1], h[1:]
    w1 = 2.0 * h1 + h0
    w2 = h1 + 2.0 * h0
    same_sign = (d0 * d1) > 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        harmonic = (w1 + w2) / (w1 / d0 + w2 / d1)
    d[1:-1] = np.where(same_sign, harmonic, 0.0)
    # end points
    d[0] = _edge_slope(h[0], h[1], delta[0], delta[1])
    d[-1] = _edge_slope(h[-1], h[-2], delta[-1], delta[-2])
    return d


def _edge_slope(h0: float, h1: float, m0: float, m1: float) -> float:
    d = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    if np.sign(d) != np.sign(m0):
        return 0.0
    if np.sign(m0) != np.sign(m1) and abs(d) > 3.0 * abs(m0):
        return 3.0 * m0
    return float(d)


def hermite_eval(t: Any, h: Any, y0: Any, y1: Any, d0: Any, d1: Any) -> np.ndarray:
    """Cubic Hermite segment value at normalised position ``t`` in [0, 1]."""
    t2 = t * t
    t3 = t2 * t
    return (
        (2.0 * t3 - 3.0 * t2 + 1.0) * y0
        + (t3 - 2.0 * t2 + t) * h * d0
        + (-2.0 * t3 + 3.0 * t2) * y1
        + (t3 - t2) * h * d1
    )


def segment_range_increments(
    x: np.ndarray, y: np.ndarray, d: np.ndarray, n_points: int = RANGE_QUADRATURE_POINTS
) -> np.ndarray:
    """``integral dE / S(E)`` over every segment of the Hermite interpolant."""
    nodes, weights = np.polynomial.legendre.leggauss(n_points)
    h = np.diff(x)
    total = np.zeros_like(h)
    for xi, wi in zip(nodes, weights, strict=True):
        t = 0.5 * (xi + 1.0)
        s = hermite_eval(t, h, y[:-1], y[1:], d[:-1], d[1:])
        total += wi / s
    return 0.5 * h * total


def prepare_stopping_table(
    energy_mev: np.ndarray,
    stopping_mev_cm2_per_g: np.ndarray,
    provenance: dict[str, Any] | None = None,
    drop_nonpositive_energy: bool = True,
) -> StoppingTable:
    """Validate a raw table and precompute PCHIP slopes and cumulative ranges."""
    e = np.asarray(energy_mev, dtype=np.float64)
    s = np.asarray(stopping_mev_cm2_per_g, dtype=np.float64)
    if e.shape != s.shape or e.ndim != 1:
        raise TableFormatError(
            "energy and stopping-power columns must be 1-D and equal"
        )
    if drop_nonpositive_energy:
        keep = e > 0.0
        e, s = e[keep], s[keep]
    if e.shape[0] < 2:
        raise TableFormatError(
            "a stopping table needs at least two positive-energy rows"
        )
    if np.any(np.diff(e) <= 0.0):
        raise TableFormatError("energies must be strictly increasing")
    if np.any(s <= 0.0) or not np.all(np.isfinite(s)):
        raise TableFormatError("stopping powers must be positive and finite")
    d = pchip_slopes(e, s)
    increments = segment_range_increments(e, s, d)
    cumulative = np.concatenate([[0.0], np.cumsum(increments)])
    return StoppingTable(
        energy_mev=e,
        stopping_mev_cm2_per_g=s,
        slope=d,
        csda_range_g_per_cm2=cumulative,
        provenance=dict(provenance or {}),
    )


def load_stopping_table(
    path: str | Path, provenance: dict[str, Any] | None = None
) -> StoppingTable:
    """Parse and prepare a two-column stopping-power file."""
    e, s = parse_two_column_table(path)
    return prepare_stopping_table(e, s, provenance)
