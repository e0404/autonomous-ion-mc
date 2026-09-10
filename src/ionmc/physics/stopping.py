"""Analytical electronic stopping power of charged particles (shared source).

This module implements the Bethe formula with the standard corrections, in the
form used by ICRU Report 49 / NIST PSTAR above their theory threshold and by
the PDG review "Passage of particles through matter":

    -(1/rho) dE/dx = K z^2 (Z/A) / beta^2 * L(beta)

    L = 1/2 ln(2 m_e c^2 beta^2 gamma^2 T_max / I^2) - beta^2 - delta/2 - C/Z
        + z L1 + L2

with ``K = 4 pi N_A r_e^2 m_e c^2 = 0.307075 MeV cm^2/mol``, ``T_max`` the
maximum energy transfer to a free electron including the projectile-mass
terms, ``delta`` the Sternheimer density-effect correction, ``C/Z`` the shell
correction, ``z L1`` the Barkas (z^3) correction and ``L2`` the Bloch (z^4)
correction (``L2`` includes its z^2 factor here).

Model choices and their sources (decision ``0006``):

* **Density effect**: Sternheimer-Berger-Seltzer parameterisation with the
  per-material constants of ``ionmc.materials``.
* **Shell correction**: the Barkas-Berger empirical formula
  ``C(I, eta) = A(eta) I^2 + B(eta) I^3`` with ``eta = beta gamma`` and ``I``
  in eV (the formula also appears in Leo, *Techniques for Nuclear and
  Particle Physics Experiments*, eq. 2.33; the coefficients used here are
  transcribed from Geant4 ``G4IonisParamElm.cc``, which is the authoritative
  source for their exact values), applied per target element with the
  ICRU 37 elemental ``I`` and combined by electron-weighted Bragg additivity
  ``sum_k n_k C_k / N_e`` (which is what Geant4's ``ShellCorrectionSTD``
  evaluates). The formula is valid for ``eta >= 0.13``; below the velocity
  of an 8 MeV proton it is frozen and tapered logarithmically to zero at the
  velocity of a 2 MeV proton. This is a velocity-based variant of the Geant4
  convention, whose thresholds are ``8 MeV / M`` and ``2 MeV / m_p`` in
  ``tau`` (identical for protons, and within 0.1 percent of the taper factor
  at 5 MeV; for heavier ions Geant4's thresholds move to lower velocity).
  The analytic layer is therefore only claimed accurate for proton energies
  >= 10 MeV (decision ``0006`` fixes the acceptance tolerances).
* **Barkas correction**: Ashley-Ritchie-Brandt ``L1 = 1.29 F(b / sqrt(X)) /
  (sqrt(Z X) X)``, ``X = beta^2 / (alpha^2 Z)``, with the tabulated ``F`` and
  the element parameter ``b`` of :mod:`ionmc.physics.barkas_table`, combined
  over target elements with electron weights (Bragg additivity per electron;
  Geant4 weights by atoms instead, which for water gives an ``L1`` smaller
  by a factor 1.6-1.8, i.e. 0.07 percent of ``S`` at 10 MeV, see decision
  ``0006``).
* **Bloch correction**: exact series ``L2 = -y^2 sum_j 1 / (j (j^2 + y^2))``,
  ``y = z alpha / beta`` (Geant4 ``BlochCorrection``), summed to 16 terms
  with an integral tail estimate.

All functions are *shared source*: they use only the math namespace ``m``
and plain ``float`` scalars (Warp: float32) or the small tables passed in as
arrays, so the same text executes as Warp functions in kernels, as pure
Python in float64 (the reference path) or vectorised with numpy. Data
dependent branching is expressed with ``m.where``. The Warp instantiation is
single precision; the float64 reference is the Python binding (decision
``0005``).

Units: energies in MeV (``I`` is passed in MeV), masses as rest energies in
MeV, stopping power in MeV cm^2/g.
"""

from __future__ import annotations

from typing import Any

from ionmc.backend import mathlib
from ionmc.constants import (
    BETHE_K_MEV_CM2_PER_MOL,
    ELECTRON_MASS_MEV,
    FINE_STRUCTURE_CONSTANT,
    PROTON_MASS_MEV,
)
from ionmc.physics.barkas_table import BARKAS_SCALE

m = mathlib.current()
func = mathlib.func

#: 2 ln(10), the slope of the density-effect asymptote.
TWO_LN10: float = 4.605170185988092

#: alpha^2
ALPHA_SQUARED: float = FINE_STRUCTURE_CONSTANT * FINE_STRUCTURE_CONSTANT


def _beta2_gamma2_of_proton(kinetic_energy_mev: float) -> float:
    tau = kinetic_energy_mev / PROTON_MASS_MEV
    return tau * (tau + 2.0)


#: (beta gamma)^2 of an 8 MeV proton: below this velocity the shell correction
#: formula is outside its fitted range and is frozen (Geant4 convention).
SHELL_BG2_FREEZE: float = _beta2_gamma2_of_proton(8.0)

#: (beta gamma)^2 of a 2 MeV proton: the frozen shell correction is tapered
#: logarithmically to zero here (Geant4 convention).
SHELL_BG2_ZERO: float = _beta2_gamma2_of_proton(2.0)

#: Number of explicit terms of the Bloch series.
BLOCH_TERMS: int = 16


@func
def lorentz_gamma(kinetic_energy: float, rest_energy: float) -> float:
    """Lorentz factor from kinetic energy and rest energy (same units)."""
    return 1.0 + kinetic_energy / rest_energy


@func
def beta_gamma_squared(kinetic_energy: float, rest_energy: float) -> float:
    """``(beta gamma)^2 = tau (tau + 2)`` with ``tau = T / (M c^2)``.

    This form has no cancellation, unlike ``gamma^2 - 1``, which loses about
    five digits in float32 for a 5 MeV proton (``gamma - 1 ~ 5e-3``).
    """
    tau = kinetic_energy / rest_energy
    return tau * (tau + 2.0)


@func
def beta_squared(kinetic_energy: float, rest_energy: float) -> float:
    """``beta^2`` from kinetic energy and rest energy (same units)."""
    gamma = lorentz_gamma(kinetic_energy, rest_energy)
    return beta_gamma_squared(kinetic_energy, rest_energy) / (gamma * gamma)


@func
def max_energy_transfer(kinetic_energy: float, rest_energy: float) -> float:
    """Maximum kinetic energy transferable to a free electron, in MeV.

    ``T_max = 2 m_e c^2 beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2)``.
    """
    gamma = lorentz_gamma(kinetic_energy, rest_energy)
    bg2 = beta_gamma_squared(kinetic_energy, rest_energy)
    ratio = ELECTRON_MASS_MEV / rest_energy
    return 2.0 * ELECTRON_MASS_MEV * bg2 / (1.0 + 2.0 * gamma * ratio + ratio * ratio)


@func
def density_effect_delta(
    x: float,
    cbar: float,
    x0: float,
    x1: float,
    a: float,
    mexp: float,
    delta0: float,
) -> float:
    """Sternheimer density-effect correction ``delta`` at ``x = log10(beta gamma)``."""
    high = TWO_LN10 * x - cbar
    mid = high + a * m.pow(m.max(x1 - x, 0.0), mexp)
    low = delta0 * m.pow(10.0, 2.0 * (x - x0))
    return m.where(x >= x1, high, m.where(x >= x0, mid, low))


@func
def shell_correction(bg2: float, shell_m2: float, shell_m3: float) -> float:
    """Shell correction per electron, ``sum_k n_k C_k(I_k, eta) / N_e``.

    ``shell_m2 = sum_k n_k I_k^2 / N_e`` and ``shell_m3 = sum_k n_k I_k^3 / N_e``
    (``I`` in eV) are material constants; ``bg2 = (beta gamma)^2``.
    """
    inv = 1.0 / m.max(bg2, SHELL_BG2_FREEZE)
    inv2 = inv * inv
    inv3 = inv2 * inv
    a_eta = (0.422377 * inv + 0.0304043 * inv2 - 0.00038106 * inv3) * 1.0e-6
    b_eta = (3.858019 * inv - 0.1667989 * inv2 + 0.00157955 * inv3) * 1.0e-9
    c = a_eta * shell_m2 + b_eta * shell_m3
    taper = m.log(m.max(bg2, SHELL_BG2_ZERO) / SHELL_BG2_ZERO) / m.log(
        SHELL_BG2_FREEZE / SHELL_BG2_ZERO
    )
    return c * m.where(bg2 >= SHELL_BG2_FREEZE, 1.0, taper)


@func
def bloch_correction(charge: float, beta2: float) -> float:
    """Bloch correction ``L2 = -y^2 sum_j 1/(j (j^2 + y^2))``, ``y = z alpha/beta``."""
    y2 = charge * charge * ALPHA_SQUARED / beta2
    # ``float(0.0)`` declares a mutable (dynamic) variable for Warp's code
    # generator; a bare literal would be a constant that loops cannot update.
    total = float(0.0)  # noqa: UP018 - dynamic variable for Warp
    for j in range(1, BLOCH_TERMS + 1):
        fj = float(j)
        total += 1.0 / (fj * (fj * fj + y2))
    # tail: sum_{j > N} 1/j^3 ~ integral from N + 1/2
    edge = float(BLOCH_TERMS) + 0.5
    total += 0.5 / (edge * edge)
    return -y2 * total


@func
def barkas_function(w: float, table_w: Any, table_f: Any, n_table: int) -> float:
    """Piecewise-linear interpolation of the tabulated ARB function ``F(W)``.

    Below the first tabulated ``W`` the first value is used; above the last,
    the value decays as ``F(W_max) W_max / W`` (Geant4 convention).
    """
    value = table_f[0]
    for i in range(n_table - 1):
        w0 = table_w[i]
        w1 = table_w[i + 1]
        t = (w - w0) / (w1 - w0)
        candidate = table_f[i] + t * (table_f[i + 1] - table_f[i])
        value = m.where(w >= w0, m.where(w < w1, candidate, value), value)
    wmax = table_w[n_table - 1]
    return m.where(w >= wmax, table_f[n_table - 1] * wmax / w, value)


@func
def barkas_l1(
    beta2: float,
    target_z: float,
    b: float,
    table_w: Any,
    table_f: Any,
    n_table: int,
) -> float:
    """Barkas term ``L1`` per unit projectile charge for one target element."""
    x = beta2 / (ALPHA_SQUARED * target_z)
    w = b / m.sqrt(x)
    return (
        BARKAS_SCALE
        * barkas_function(w, table_w, table_f, n_table)
        / (m.sqrt(target_z * x) * x)
    )


@func
def stopping_number(
    kinetic_energy: float,
    rest_energy: float,
    charge: float,
    i_mev: float,
    cbar: float,
    x0: float,
    x1: float,
    a: float,
    mexp: float,
    delta0: float,
    shell_m2: float,
    shell_m3: float,
    elem_z: Any,
    elem_b: Any,
    elem_f: Any,
    n_elem: int,
    table_w: Any,
    table_f: Any,
    n_table: int,
    use_density: float,
    use_shell: float,
    use_barkas: float,
    use_bloch: float,
) -> float:
    """Stopping number ``L`` of the Bethe formula with selectable corrections.

    The ``use_*`` switches are 0.0 or 1.0 so that the same source serves the
    configurable-fidelity requirement without branching.
    """
    gamma = lorentz_gamma(kinetic_energy, rest_energy)
    bg2 = beta_gamma_squared(kinetic_energy, rest_energy)
    beta2 = bg2 / (gamma * gamma)
    tmax = max_energy_transfer(kinetic_energy, rest_energy)
    l0 = 0.5 * m.log(2.0 * ELECTRON_MASS_MEV * bg2 * tmax / (i_mev * i_mev)) - beta2
    x = 0.5 * m.log10(bg2)
    delta = density_effect_delta(x, cbar, x0, x1, a, mexp, delta0)
    shell = shell_correction(bg2, shell_m2, shell_m3)
    barkas = float(0.0)  # noqa: UP018 - dynamic variable for Warp
    for k in range(n_elem):
        barkas += elem_f[k] * barkas_l1(
            beta2, elem_z[k], elem_b[k], table_w, table_f, n_table
        )
    bloch = bloch_correction(charge, beta2)
    return (
        l0
        - 0.5 * delta * use_density
        - shell * use_shell
        + charge * barkas * use_barkas
        + bloch * use_bloch
    )


@func
def mass_stopping_power(
    kinetic_energy: float,
    rest_energy: float,
    charge: float,
    za_ratio: float,
    i_mev: float,
    cbar: float,
    x0: float,
    x1: float,
    a: float,
    mexp: float,
    delta0: float,
    shell_m2: float,
    shell_m3: float,
    elem_z: Any,
    elem_b: Any,
    elem_f: Any,
    n_elem: int,
    table_w: Any,
    table_f: Any,
    n_table: int,
    use_density: float,
    use_shell: float,
    use_barkas: float,
    use_bloch: float,
) -> float:
    """Electronic mass stopping power in MeV cm^2/g."""
    beta2 = beta_squared(kinetic_energy, rest_energy)
    number = stopping_number(
        kinetic_energy,
        rest_energy,
        charge,
        i_mev,
        cbar,
        x0,
        x1,
        a,
        mexp,
        delta0,
        shell_m2,
        shell_m3,
        elem_z,
        elem_b,
        elem_f,
        n_elem,
        table_w,
        table_f,
        n_table,
        use_density,
        use_shell,
        use_barkas,
        use_bloch,
    )
    return BETHE_K_MEV_CM2_PER_MOL * charge * charge * za_ratio / beta2 * number


@func
def csda_range_increment(
    energy_low: float,
    energy_high: float,
    n_steps: int,
    rest_energy: float,
    charge: float,
    za_ratio: float,
    i_mev: float,
    cbar: float,
    x0: float,
    x1: float,
    a: float,
    mexp: float,
    delta0: float,
    shell_m2: float,
    shell_m3: float,
    elem_z: Any,
    elem_b: Any,
    elem_f: Any,
    n_elem: int,
    table_w: Any,
    table_f: Any,
    n_table: int,
    use_density: float,
    use_shell: float,
    use_barkas: float,
    use_bloch: float,
) -> float:
    """``integral_{E_low}^{E_high} dE / S(E)`` in g/cm^2 by composite Simpson.

    The integration variable is ``u = ln E`` (``dE/S = E/S du``), on
    ``n_steps`` (even) uniform sub-intervals. With ``n_steps = 200`` the
    quadrature error is below 1e-6 relative for 1-1000 MeV protons.
    """
    u0 = m.log(energy_low)
    u1 = m.log(energy_high)
    h = (u1 - u0) / float(n_steps)
    total = float(0.0)  # noqa: UP018 - dynamic variable for Warp
    for i in range(n_steps + 1):
        u = u0 + h * float(i)
        e = m.exp(u)
        s = mass_stopping_power(
            e,
            rest_energy,
            charge,
            za_ratio,
            i_mev,
            cbar,
            x0,
            x1,
            a,
            mexp,
            delta0,
            shell_m2,
            shell_m3,
            elem_z,
            elem_b,
            elem_f,
            n_elem,
            table_w,
            table_f,
            n_table,
            use_density,
            use_shell,
            use_barkas,
            use_bloch,
        )
        weight = 2.0 + 2.0 * float(i % 2)
        weight = m.where(i == 0, 1.0, m.where(i == n_steps, 1.0, weight))
        total += weight * e / s
    return total * h / 3.0
