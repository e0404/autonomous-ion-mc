"""Shared-source physics of the continuous-slowing-down transport step.

These functions are written against the math namespace ``m`` and decorated with
``@func`` (decision ``0005``), so the identical text runs as a Warp function per
thread, as float64 Python (the reference oracle) and vectorised with numpy.
They contain the *physics* of one transport step; the stepping loop, atomic
scoring and thread indexing live in the Warp kernel and in the mirrored Python
driver (decision ``0009``).

The stopping power is looked up from a prepared table
(:mod:`ionmc.physics.tabulated`) passed in as flat arrays. Units: energy in
MeV, length in mm, mass density in g/cm^3, mass stopping power in MeV cm^2/g,
linear stopping power in MeV/mm (decision ``0004``).
"""

from __future__ import annotations

from typing import Any

from ionmc.backend import mathlib
from ionmc.physics.tabulated import tabulated_mass_stopping_power

m = mathlib.current()
func = mathlib.func

#: mm per cm; linear stopping power [MeV/mm] = mass SP [MeV cm^2/g] * rho / 10.
MM_PER_CM: float = 10.0


@func
def linear_stopping_power(
    energy: float,
    density: float,
    table_e: Any,
    table_s: Any,
    table_d: Any,
    n: int,
    n_steps: int,
) -> float:
    """Linear stopping power ``S_lin = S_mass * rho / 10`` [MeV/mm]."""
    s_mass = tabulated_mass_stopping_power(
        energy, table_e, table_s, table_d, n, n_steps
    )
    return s_mass * density / MM_PER_CM


@func
def energy_loss_step_length(
    energy: float,
    max_fraction: float,
    max_step_mm: float,
    density: float,
    table_e: Any,
    table_s: Any,
    table_d: Any,
    n: int,
    n_steps: int,
) -> float:
    """Step length [mm] that loses at most ``max_fraction`` of ``energy``.

    ``dl = max_fraction * E / S_lin(E)``, capped at ``max_step_mm``. The caller
    additionally limits the step to the next depth-bin boundary and the
    geometry exit (decision ``0009``).
    """
    s_lin = linear_stopping_power(
        energy, density, table_e, table_s, table_d, n, n_steps
    )
    dl_energy = max_fraction * energy / s_lin
    return m.min(dl_energy, max_step_mm)


@func
def midpoint_energy_loss(
    energy: float,
    step_mm: float,
    density: float,
    table_e: Any,
    table_s: Any,
    table_d: Any,
    n: int,
    n_steps: int,
) -> float:
    """Energy lost [MeV] over a step of ``step_mm`` by the midpoint (RK2) rule.

    ``S`` is evaluated at the energy after half the Euler loss, so the estimate
    is second-order accurate in the step length. The result is clamped to
    ``energy`` so a step never removes more than the particle has.
    """
    s0 = linear_stopping_power(energy, density, table_e, table_s, table_d, n, n_steps)
    e_mid = m.max(energy - 0.5 * s0 * step_mm, 0.0)
    s_mid = linear_stopping_power(e_mid, density, table_e, table_s, table_d, n, n_steps)
    return m.min(s_mid * step_mm, energy)
