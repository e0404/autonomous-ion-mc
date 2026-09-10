"""Tabulated stopping power and CSDA range by interpolation (shared source).

The functions take a prepared table (:class:`ionmc.data.stopping_tables.
StoppingTable`) as flat arrays and evaluate

* the mass stopping power by **monotone cubic Hermite (PCHIP)
  interpolation** between grid points, with the slopes precomputed at
  preparation time; the interpolant reproduces the grid values exactly, is
  shape preserving and needs no transcendental functions;
* the CSDA range as the precomputed cumulative range at the lower grid point
  plus a 4-point Gauss-Legendre integral of ``1/S`` over the partial segment
  (decision ``0008``).

Locating the segment uses a **fixed-count bisection**: ``n_steps =
ceil(log2(n))`` iterations of ``lo/hi`` updates written with ``m.where``, so
the loop has a data-independent trip count and the same source runs as a Warp
function (per thread), as pure Python (float64 reference) and vectorised
with numpy (decision ``0005`` shared-source rules).

Energies outside the table evaluate the boundary segment's cubic
(extrapolation); callers decide whether that is acceptable
(``TabulatedStoppingPower`` raises outside the table by default).

Units: energies in MeV, stopping power in MeV cm^2/g, range in g/cm^2.
"""

from __future__ import annotations

from typing import Any

from ionmc.backend import mathlib

m = mathlib.current()
func = mathlib.func

#: 4-point Gauss-Legendre nodes (on [-1, 1]) and weights for the partial segment.
GL4_X0: float = -0.8611363115940526
GL4_X1: float = -0.3399810435848563
GL4_X2: float = 0.3399810435848563
GL4_X3: float = 0.8611363115940526
GL4_W0: float = 0.3478548451374538
GL4_W1: float = 0.6521451548625461
GL4_W2: float = 0.6521451548625461
GL4_W3: float = 0.3478548451374538


@func
def locate_segment(e: float, table_e: Any, n: int, n_steps: int) -> int:
    """Index ``i`` in ``[0, n - 2]`` with ``table_e[i] <= e < table_e[i+1]``.

    Bisection with ``n_steps`` iterations (``2^n_steps >= n``); energies below
    the table give 0 and energies above give ``n - 2``.
    """
    # ``lo`` and ``hi`` are mutated inside the loop, so they must be dynamic
    # variables for Warp: initialising with ``int(...)`` marks them mutable
    # (Warp refuses to mutate a constant-initialised variable inside a dynamic
    # loop; the same shared-source rule as decision 0006). ``noqa: UP018`` keeps
    # ruff from rewriting ``int(0)`` back to the constant ``0``.
    lo = int(0)  # noqa: UP018, RUF046 - int() marks the variable mutable for Warp
    hi = int(n - 1)
    for _ in range(n_steps):
        mid = (lo + hi) // 2
        go_right = table_e[mid] <= e
        lo = m.where(go_right, mid, lo)
        hi = m.where(go_right, hi, mid)
    return m.min(lo, n - 2)


@func
def hermite_segment(
    t: float, h: float, y0: float, y1: float, d0: float, d1: float
) -> float:
    """Cubic Hermite value at normalised position ``t`` of a segment of width ``h``."""
    t2 = t * t
    t3 = t2 * t
    return (
        (2.0 * t3 - 3.0 * t2 + 1.0) * y0
        + (t3 - 2.0 * t2 + t) * h * d0
        + (-2.0 * t3 + 3.0 * t2) * y1
        + (t3 - t2) * h * d1
    )


@func
def tabulated_mass_stopping_power(
    kinetic_energy: float,
    table_e: Any,
    table_s: Any,
    table_d: Any,
    n: int,
    n_steps: int,
) -> float:
    """Mass stopping power [MeV cm^2/g] at ``kinetic_energy`` [MeV]."""
    i = locate_segment(kinetic_energy, table_e, n, n_steps)
    e0 = table_e[i]
    h = table_e[i + 1] - e0
    t = (kinetic_energy - e0) / h
    return hermite_segment(t, h, table_s[i], table_s[i + 1], table_d[i], table_d[i + 1])


@func
def tabulated_csda_range(
    kinetic_energy: float,
    table_e: Any,
    table_s: Any,
    table_d: Any,
    table_range: Any,
    n: int,
    n_steps: int,
) -> float:
    """CSDA range [g/cm^2] from the first grid energy up to ``kinetic_energy``.

    ``R(E) = R_i + integral_{E_i}^{E} dE'/S(E')`` with ``S`` the Hermite
    interpolant of segment ``i``, by 4-point Gauss-Legendre quadrature.
    """
    i = locate_segment(kinetic_energy, table_e, n, n_steps)
    e0 = table_e[i]
    h = table_e[i + 1] - e0
    y0 = table_s[i]
    y1 = table_s[i + 1]
    d0 = table_d[i]
    d1 = table_d[i + 1]
    width = kinetic_energy - e0
    half = 0.5 * width
    center = e0 + half
    total = float(0.0)  # noqa: UP018 - dynamic variable for Warp
    total += GL4_W0 / hermite_segment(
        (center + half * GL4_X0 - e0) / h, h, y0, y1, d0, d1
    )
    total += GL4_W1 / hermite_segment(
        (center + half * GL4_X1 - e0) / h, h, y0, y1, d0, d1
    )
    total += GL4_W2 / hermite_segment(
        (center + half * GL4_X2 - e0) / h, h, y0, y1, d0, d1
    )
    total += GL4_W3 / hermite_segment(
        (center + half * GL4_X3 - e0) / h, h, y0, y1, d0, d1
    )
    return table_range[i] + half * total
