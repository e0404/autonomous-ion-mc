"""Offline LET from fluence spectra (row A8 of ``validation/plans/v3-004-acceptance.md``).

A ``fluence_spectrum`` tally (decision 0040) stores, per voxel, the path length ``Phi_b`` that was
travelled with a piece-mean energy per nucleon in energy bin ``b`` (plus an underflow and an
overflow bin). With the unrestricted water stopping power ``S_w`` evaluated at the centre of each
bin, the track- and dose-averaged LET follow without a transport run::

    LET_t = sum_b Phi_b S_b / sum_b Phi_b        LET_d = sum_b Phi_b S_b^2 / sum_b Phi_b S_b

This is the post-processing route of a stored spectrum; it differs from the online value only by
the discretisation of the energy axis (a second-order term in the bin width, because ``S`` is
sampled at the centre of a bin) and by the ramp curvature inside a piece. The centre of a log axis
is the geometric mean of the edges, that of a linear axis the arithmetic mean. The underflow and
overflow weights are not part of the sums: they are returned so that a caller can require them to
be negligible.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True, eq=False)
class OfflineLET:
    """Result of :func:`let_from_spectrum`: ``let_t`` and ``let_d`` (NaN where the in-range weight
    is zero), the per-point in-range weight and the relative weight outside the edge range."""

    let_t: NDArray[np.float64]
    let_d: NDArray[np.float64]
    weight: NDArray[np.float64]
    outside_fraction: NDArray[np.float64]


def log_edges(e_min: float, bins_per_decade: int, n_decades: float) -> NDArray[np.float64]:
    """Uniform logarithmic bin edges ``e_min 10^(i / bins_per_decade)``, ``i = 0 .. n`` with
    ``n = round(bins_per_decade n_decades)`` bins."""
    if not (e_min > 0.0 and bins_per_decade >= 1 and n_decades > 0.0):
        raise ValueError("need e_min > 0, bins_per_decade >= 1 and n_decades > 0")
    n = round(bins_per_decade * n_decades)
    return e_min * 10.0 ** (np.arange(n + 1, dtype=np.float64) / bins_per_decade)


def bin_centres(edges: Sequence[float], *, log_axis: bool = True) -> NDArray[np.float64]:
    """Centres of the bins of ``edges`` (geometric mean for a log axis, else arithmetic)."""
    e = np.asarray(edges, dtype=np.float64)
    if e.ndim != 1 or e.size < 2 or np.any(np.diff(e) <= 0.0) or e[0] <= 0.0:
        raise ValueError("edges must be >= 2 strictly increasing positive values")
    return np.sqrt(e[:-1] * e[1:]) if log_axis else 0.5 * (e[:-1] + e[1:])


def let_from_spectrum(
    spectrum: NDArray[np.float64],
    edges_mev_per_u: Sequence[float],
    s_of_energy: Callable[[float], float],
    *,
    log_axis: bool = True,
    a_nucleon: int = 1,
) -> OfflineLET:
    """``LET_t`` and ``LET_d`` from a fluence spectrum.

    ``spectrum`` has the bins on its last axis, ``len(edges) + 1`` entries (underflow, the bins,
    overflow), any leading shape (voxels). ``s_of_energy`` maps a kinetic energy [MeV] of the
    projectile to the water stopping power (for example ``tables.s_water``), evaluated at
    ``a_nucleon`` times the bin centre."""
    sp = np.asarray(spectrum, dtype=np.float64)
    nb = len(edges_mev_per_u) - 1
    if sp.shape[-1] != nb + 2:
        raise ValueError(f"spectrum needs {nb + 2} entries on its last axis, got {sp.shape[-1]}")
    if not np.all(np.isfinite(sp)) or np.any(sp < 0.0):
        raise ValueError("spectrum weights must be finite and non-negative")
    centres = bin_centres(edges_mev_per_u, log_axis=log_axis)
    s = np.array([s_of_energy(float(a_nucleon * c)) for c in centres], dtype=np.float64)
    if not np.all(np.isfinite(s)) or np.any(s <= 0.0):
        raise ValueError("stopping power must be positive and finite at every bin centre")
    inside = sp[..., 1:-1]
    outside = sp[..., 0] + sp[..., -1]
    w = inside.sum(axis=-1)
    m1 = inside @ s
    m2 = inside @ (s * s)
    with np.errstate(divide="ignore", invalid="ignore"):
        let_t = np.where(w > 0.0, m1 / np.where(w > 0.0, w, 1.0), math.nan)
        let_d = np.where(m1 > 0.0, m2 / np.where(m1 > 0.0, m1, 1.0), math.nan)
        total = w + outside
        out_frac = np.where(total > 0.0, outside / np.where(total > 0.0, total, 1.0), math.nan)
    return OfflineLET(let_t, let_d, w, out_frac)
