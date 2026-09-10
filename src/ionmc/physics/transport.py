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
from ionmc.constants import BETHE_K_MEV_CM2_PER_MOL, ELECTRON_MASS_MEV
from ionmc.physics.tabulated import tabulated_mass_stopping_power

m = mathlib.current()
func = mathlib.func

#: mm per cm; linear stopping power [MeV/mm] = mass SP [MeV cm^2/g] * rho / 10.
MM_PER_CM: float = 10.0

#: Bohr energy-straggling constant ``K m_e c^2 = 4 pi r_e^2 (m_e c^2)^2 N_A``
#: [MeV^2 cm^2/mol]; the Bohr variance per unit mass thickness is this times
#: ``(Z/A) z^2 f(beta)`` (decision 0010).
BOHR_K_MEV2_CM2_PER_MOL: float = BETHE_K_MEV_CM2_PER_MOL * ELECTRON_MASS_MEV


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


@func
def _beta_squared(energy: float, rest_energy: float) -> float:
    """``beta^2`` from kinetic and rest energy, cancellation-free (decision 0006)."""
    tau = energy / rest_energy
    gamma = 1.0 + tau
    return tau * (tau + 2.0) / (gamma * gamma)


@func
def bohr_straggling_sigma(
    energy: float,
    rest_energy: float,
    step_mm: float,
    density: float,
    za_ratio: float,
    charge: float,
) -> float:
    """Standard deviation [MeV] of the Bohr energy-loss straggling over a step.

    Bohr variance per unit path length ``dOmega^2/dx =
    K m_e c^2 (Z/A) rho z^2 f(beta)`` with ``K m_e c^2`` = ``BOHR_K``
    [MeV^2 cm^2/mol] and the relativistic factor
    ``f(beta) = (1 - beta^2/2)/(1 - beta^2)`` (decision 0010). ``step_mm`` is
    converted to cm; the result is ``sqrt(dOmega^2/dx * step)``.
    """
    beta2 = _beta_squared(energy, rest_energy)
    f_rel = (1.0 - 0.5 * beta2) / (1.0 - beta2)
    variance = (
        BOHR_K_MEV2_CM2_PER_MOL
        * za_ratio
        * density
        * (step_mm / MM_PER_CM)
        * charge
        * charge
        * f_rel
    )
    return m.sqrt(m.max(variance, 0.0))


@func
def straggled_energy_loss(
    mean_loss: float, sigma: float, variate: float, energy: float
) -> float:
    """Energy lost [MeV] over a step with Gaussian energy-loss straggling.

    ``mean_loss + variate * sigma``, clamped to ``[0, energy]`` (a step never
    removes more than the particle has, nor adds energy). ``variate`` is a
    standard-normal sample drawn in the execution layer (decision 0005: the RNG
    is not called from shared source; the variate is passed in).
    """
    loss = mean_loss + variate * sigma
    return m.min(m.max(loss, 0.0), energy)


#: Highland/Lynch-Dahl leading constant [MeV] and log coefficient (decision 0011).
HIGHLAND_CONSTANT_MEV: float = 13.6
HIGHLAND_LOG_COEFF: float = 0.038


@func
def momentum_times_velocity(energy: float, rest_energy: float) -> float:
    """``beta c p = p v = T (T + 2 M) / (T + M)`` [MeV] (cancellation-free)."""
    return energy * (energy + 2.0 * rest_energy) / (energy + rest_energy)


@func
def highland_theta0(
    energy: float,
    rest_energy: float,
    charge: float,
    step_mm: float,
    density: float,
    radiation_length_g_per_cm2: float,
) -> float:
    """Projected RMS multiple-scattering angle [rad] over a step (decision 0011).

    ``theta0 = (13.6 MeV / (beta c p)) z sqrt(x / X0)`` with the mass thickness
    ``x = step_mm/10 * density`` [g/cm^2] and ``X0`` the radiation length
    [g/cm^2]. This is the **projected** (one-plane) RMS angle; the space-angle
    variance is ``2 theta0^2``. The Lynch-Dahl logarithmic bracket is omitted
    (the scattering-power form): it is not additive across the ~1 mm steps used
    here and its per-step value is pathological (negative log); the bracket-1
    form matches published lateral spreads within a few percent and makes the
    per-step model and the Fermi-Eyges oracle use the identical scattering power
    (decision 0011).
    """
    pv = momentum_times_velocity(energy, rest_energy)
    x_over_x0 = (step_mm / MM_PER_CM) * density / radiation_length_g_per_cm2
    return HIGHLAND_CONSTANT_MEV / pv * charge * m.sqrt(m.max(x_over_x0, 0.0))
