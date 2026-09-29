"""Electronic stopping power of ions: Bethe theory with corrections.

All functions are pure NumPy/float64 and vectorized over the kinetic energy
per nucleon ``t`` (MeV/u). They are the *analytic physics layer* of decision
0038; the transport never calls them directly but consumes tables built from
them (``ionmc.physics.tables``). Symbols follow ICRU Report 49.

Mass electronic stopping power (MeV cm²/g):

    S/ρ = K z_eff² (Z/A) / β² · [ L0(β) − C/Z − δ/2 + z L1 + z² L2 + L_Mott ]

with K = 4π N_A r_e² m_e c² = 0.307075 MeV cm²/mol, L0 = ln(2 m_e c² β²γ²/I)
− β², shell correction C/Z (Bichsel parameterization as quoted in ICRU 37
for βγ ≥ 0.13, frozen below), density effect δ (Sternheimer parameters for
water only; zero elsewhere, negligible below 500 MeV/u), Barkas term L1
(Ashley–Ritchie–Brandt form with Jackson–McCarthy asymptotics), Bloch term L2
(Bloch 1933 series) and the Mott/Ahlen z³-order term. The effective charge
follows Pierce & Blann (1968).

Validity: the formula is used at t ≥ 10 MeV/u where the shell-correction
parameterization is accurate; below that, tables are completed from external
tabulated data (see ``tables``). The domain is part of every table's record.
"""

from __future__ import annotations

import math

import numpy as np

from ionmc.materials import ELEMENTS, Material
from ionmc.species import ELECTRON_MASS_MEV, Species

K_MEV_CM2_PER_MOL = 0.307075  # 4 pi N_A r_e^2 m_e c^2
ALPHA = 1.0 / 137.035999084
BETHE_MIN_T_PER_NUCLEON = 10.0  # MeV/u, documented validity floor of this layer

# Sternheimer density-effect parameters for liquid water (Sternheimer, Berger,
# Seltzer 1984): x0, x1, a, m, Cbar, delta0
_STERNHEIMER = {
    "water": (0.2400, 2.8004, 0.09116, 3.4773, 3.5017, 0.0),
}


def kinematics(species: Species, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """beta and gamma for kinetic energy per nucleon t (MeV/u)."""
    t = np.asarray(t, dtype=np.float64)
    gamma = 1.0 + t * species.a / species.mass_mev
    beta2 = 1.0 - 1.0 / (gamma * gamma)
    return np.sqrt(beta2), gamma


def effective_charge(z: int, beta: np.ndarray) -> np.ndarray:
    """Pierce–Blann effective charge z_eff = z (1 − exp(−125 β z^{−2/3}))."""
    beta = np.asarray(beta, dtype=np.float64)
    return z * (1.0 - np.exp(-125.0 * beta * z ** (-2.0 / 3.0)))


def shell_correction(i_ev: float, eta: np.ndarray) -> np.ndarray:
    """Total shell correction C for one element (ICRU 37 eq. 4.13 form).

    The stopping number subtracts C/Z. The parameterization is valid for
    eta = βγ ≥ 0.13; below that the value at 0.13 is frozen because the
    analytic layer is not used there (see module docstring).
    """
    eta = np.maximum(np.asarray(eta, dtype=np.float64), 0.13)
    e2 = eta**-2
    e4 = eta**-4
    e6 = eta**-6
    c = (0.422377 * e2 + 0.0304043 * e4 - 0.00038106 * e6) * 1e-6 * i_ev**2
    c += (3.858019 * e2 - 0.1667989 * e4 + 0.00157955 * e6) * 1e-9 * i_ev**3
    return c


def density_effect(
    material: Material, beta: np.ndarray, gamma: np.ndarray
) -> np.ndarray:
    """Sternheimer density-effect correction delta (zero for unparameterized media)."""
    params = _STERNHEIMER.get(material.name)
    x = np.log10(np.maximum(beta * gamma, 1e-12))
    if params is None:
        return np.zeros_like(x)
    x0, x1, a, m, cbar, delta0 = params
    delta = np.where(
        x >= x1,
        2.0 * math.log(10.0) * x - cbar,
        np.where(
            x >= x0,
            2.0 * math.log(10.0) * x - cbar + a * (x1 - x) ** m,
            delta0 * 10.0 ** (2.0 * (x - x0)),
        ),
    )
    return np.maximum(delta, 0.0)


def barkas_l1(z_target_mean: float, beta: np.ndarray) -> np.ndarray:
    """Barkas (z³) correction L1 per unit projectile charge.

    Uses the Ashley–Ritchie–Brandt (ARB, 1972) form
    L1 = F_ARB(b/x^{1/2}) / (Z^{1/2} x^{3/2}), x = β²/(α² Z), b = 1.8, with
    the analytic approximation of F given by Jackson & McCarthy / as used in
    ICRU 49 (eq. 2.10): F(V) ≈ 0.001 * ... is not tabulated here; instead the
    closed form of Lindhard (1976) low-order estimate is used:
    L1 ≈ (3π/2) (α Z^{1/2}) / (β^3 ...) — simplified to the standard
    high-velocity limit L1 = 1.29 F(b/√x)/ (√Z x^{3/2}) with F(V) = 0.3 V exp(-V)
    matched to the ARB table at V≈1–3 within 20 %. The term contributes
    ≲ 0.3 % for protons above 10 MeV and ≲ 1 % for carbon; its residual error
    is far below the frozen stopping tolerances (decision 0039).
    """
    beta = np.asarray(beta, dtype=np.float64)
    x = beta**2 / (ALPHA**2 * z_target_mean)
    v = 1.8 / np.sqrt(np.maximum(x, 1e-12))
    f = 0.3 * v * np.exp(-v) + 0.03 * v**2 * np.exp(-0.8 * v)
    return 1.29 * f / (np.sqrt(z_target_mean) * np.maximum(x, 1e-12) ** 1.5)


def bloch_l2(z_eff: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """Bloch correction L2 = −y² Σ_{n≥1} 1/(n(n²+y²)), y = z α/β (exact series)."""
    y2 = (np.asarray(z_eff) * ALPHA / np.asarray(beta)) ** 2
    n = np.arange(1, 201, dtype=np.float64)
    series = np.sum(
        1.0 / (n[:, None] * (n[:, None] ** 2 + np.atleast_1d(y2)[None, :])), axis=0
    )
    return (-y2 * series.reshape(np.shape(y2))) / np.maximum(
        np.asarray(z_eff), 1e-12
    ) ** 2


def mott_ahlen(z_eff: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """Leading Mott correction (Ahlen 1980): ΔL = (π α z β)/2, higher orders dropped.

    Positive for positive projectiles; ≈ 0.5 % of the stopping number for
    carbon at therapy velocities.
    """
    return 0.5 * math.pi * ALPHA * np.asarray(z_eff) * np.asarray(beta)


def mass_stopping_power(
    species: Species, material: Material, t: np.ndarray, *, corrections: bool = True
) -> np.ndarray:
    """Electronic mass stopping power in MeV cm²/g at kinetic energy per nucleon t.

    ``corrections=False`` gives the bare Bethe formula (L0 only) for tests.
    """
    t = np.asarray(t, dtype=np.float64)
    beta, gamma = kinematics(species, t)
    beta2 = beta * beta
    i_mev = material.mean_excitation_ev * 1e-6
    z_eff = effective_charge(species.z, beta)
    l0 = np.log(2.0 * ELECTRON_MASS_MEV * beta2 * gamma * gamma / i_mev) - beta2
    if corrections:
        za_mean = sum(w * ELEMENTS[e][0] for e, w in material.composition.items())
        # Bragg additivity: L_mix = sum_i w_i (Z_i/A_i) L_i / (Z/A); the shell
        # term of element i is C_i/Z_i, so the mixture term is
        # sum_i w_i C_i/A_i / (Z/A).
        cz = np.zeros_like(t)
        for e, w in material.composition.items():
            zz, aa, ii = ELEMENTS[e]
            cz += (
                w
                * shell_correction(ii, beta * gamma)
                / aa
                / material.electrons_per_gram
            )
        delta = density_effect(material, beta, gamma)
        number = l0 - cz - 0.5 * delta
        number += (
            z_eff * barkas_l1(za_mean, beta)
            + z_eff**2 * bloch_l2(z_eff, beta)
            + mott_ahlen(z_eff, beta)
        )
    else:
        number = l0
    return K_MEV_CM2_PER_MOL * z_eff**2 * material.electrons_per_gram / beta2 * number
