"""Lindhard-Sorensen (LS) stopping correction for a point nucleus, and the residual
LS - (Bloch + first-order Mott) omitted by the ionmc analytic Bethe model (V3-008, check P3).

All quantities are corrections to the stopping number L in the convention used by ionmc
(``src/ionmc/physics/stopping.py``)::

    S/rho = K q^2 (Z/A) beta^-2 [ L0 + q L1 + Bloch + Mott + ... ]

i.e. the "4 pi r_e^2" convention in which the first-order Mott term is ``(pi/2) alpha beta z``.
(Geant4 uses the "2 pi r_e^2" convention, in which every Delta L is twice as large.)

Equations
---------
[LS96] J. Lindhard, A. H. Sorensen, Phys. Rev. A 53, 2443 (1996).
[KHV]  P. B. Kats, K. V. Halenka, O. O. Voskresenskaya, arXiv:2008.04819v2, eqs. (9)-(12)
       (eq. (10) is the point-nucleus LS correction "as represented in [LS96]").
[MARS] I. L. Rakhno, N. V. Mokhov, S. I. Striganov, arXiv:physics/0505043, eq. (LS) (same sum,
       written term by term).

Point-nucleus LS correction ([KHV] eq. 10; [MARS] eq. LS)::

    dL_LS = beta^2/2 + sum_{k>=1} [ k/eta^2 (k-1)/(2k-1) sin^2(d_k - d_{k-1})
                                  + k/eta^2 (k+1)/(2k+1) sin^2(d_{-k} - d_{-k-1})
                                  + k / ((4k^2-1)(gamma^2 k^2 + eta^2)) - 1/k ]

with eta = alpha z / beta and the Dirac-Coulomb phase shifts ([KHV] eqs. 9-10)::

    d_k = phi_k - arg Gamma(rho_k + 1 + i eta) + (pi/2)(l - rho_k),
    exp(2 i phi_k) = (k - i xi) / (rho_k - i eta),   xi = eta/gamma,
    rho_k = sqrt(k^2 - (alpha z)^2),   l = k (k > 0),  l = -k - 1 (k < 0).

Only sin^2 of phase differences enters, so every phase is needed modulo pi; d_0 is never
needed (its coefficient k-1 vanishes at k = 1).

Non-relativistic limit (alpha z -> 0 at fixed eta, gamma -> 1): phi_k = 0 (k>0), d_k = -arg
Gamma(l+1+i eta) for both signs of k, sin^2 of both differences = eta^2/(k^2+eta^2), the
k-coefficients sum to 1 and the summand becomes k/(k^2+eta^2) - 1/k: the Bloch correction
([KHV] eq. 11)  dL_B = -eta^2 sum_k 1/(k (k^2 + eta^2)).  (Checked algebraically and in the tests.)

First-order Mott term (Ahlen, Rev. Mod. Phys. 52, 121 (1980); [KHV] eq. 5; Geant4
G4EmCorrections::MottCorrection in the 2x convention):  dL_M1 = (pi/2) alpha beta z.

Numerics
--------
* Phases in numpy ``longdouble`` (x86-64 80-bit, eps 1.1e-19). ``l - rho_k`` is evaluated as
  ``(alpha z)^2/(|k| + rho_k)`` (minus 1 for k < 0) to avoid cancellation at large k.
* Complex log-Gamma: Stirling series with 8 Bernoulli terms after upward recurrence to
  Re(w) >= 16 (truncation error < 1e-21 relative to the last kept term's scale).
* The summand decays as c2/k^2 (observed numerically: c2 equals the first-order Mott term to
  the precision of the fit; the partial-wave sum of the Mott interference term converges slowly,
  as is known for Mott series). The sum is taken explicitly to K and the tail
  sum_{k>K} t_k is added from a fit t_k = c2/k^2 + c3/k^3 + c4/k^4 to the summand at
  k = K, K/2, K/4, with sum_{k>K} k^-p from the Euler-Maclaurin formula.
* Convergence criterion: K is doubled from ``k_start`` until the tail-corrected sums at K and
  K/2 agree to ``tol`` (default 1e-11 absolute in Delta L); otherwise RuntimeError.

Finite nuclear size is neglected: it modifies the phase shifts only when gamma m_e c R / hbar is
of order 1/2 ([MARS]); for 16O (R ~ 3 fm) at <= 430 MeV/u (gamma <= 1.47) that ratio is
~0.011, i.e. gamma would have to be ~60 times larger.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# Constants identical to ionmc (src/ionmc/physics/stopping.py, projectiles.py).
ALPHA_FS = 1.0 / 137.035999084
ELECTRON_MASS_MEV = 0.51099895000
AMU_MEV = 931.49410242

#: (mass number, nuclear mass in MeV) exactly as in ionmc ``physics/projectiles.py``.
ION_MASS: dict[int, tuple[int, float]] = {
    1: (1, 938.27208816),
    2: (4, 3727.3794066),
    6: (12, 12.0 * AMU_MEV - 6 * ELECTRON_MASS_MEV),
    8: (16, 15.99491461957 * AMU_MEV - 8 * ELECTRON_MASS_MEV),
}

_LD = np.longdouble
_CLD = np.clongdouble
# Stirling coefficients B_2n / (2n (2n-1)), n = 1..8.
_STIRLING = ((1, 12), (-1, 360), (1, 1260), (-1, 1680), (1, 1188), (-691, 360360), (1, 156),
             (-3617, 122400))  # fmt: skip
_HALF_LN_2PI = _LD("0.918938533204672741780329736405617639861")


def clgamma(z: np.ndarray) -> np.ndarray:
    """Complex log-Gamma (clongdouble) for Re z > 0, imaginary part correct modulo 2 pi.

    Recurrence ln Gamma(z) = ln Gamma(z + n) - sum_{j<n} ln(z + j) up to Re >= 16, then the
    Stirling series ln Gamma(w) = (w - 1/2) ln w - w + ln(2 pi)/2 + sum_n B_2n/(2n(2n-1) w^(2n-1)).
    """
    w = np.array(z, dtype=_CLD, copy=True, ndmin=1)
    if np.any(w.real <= 0):
        raise ValueError("clgamma needs Re z > 0")
    acc = np.zeros_like(w)
    for _ in range(64):
        m = w.real < 16
        if not m.any():
            break
        acc[m] += np.log(w[m])
        w[m] += 1
    lw = np.log(w)
    r = (w - _LD(0.5)) * lw - w + _HALF_LN_2PI
    inv = 1 / w
    inv_w2 = inv * inv
    p = inv
    for a, b in _STIRLING:
        r = r + (_LD(a) / _LD(b)) * p
        p = p * inv_w2
    return r - acc


def kinematics(z: int, e_mev_per_u: float) -> tuple[float, float]:
    """Return (beta, gamma) for nucleus ``z`` (1, 2, 6, 8; ionmc masses) at E [MeV/u]."""
    a, m = ION_MASS[z]
    gamma = 1.0 + e_mev_per_u * a / m
    beta = math.sqrt(1.0 - 1.0 / (gamma * gamma))
    return beta, gamma


def ls_summand(alpha_z: float, beta: float, k_max: int) -> tuple[np.ndarray, float, float]:
    """Summand t_k (k = 1..k_max, longdouble) of the LS sum, plus beta^2/2 and eta.

    ``alpha_z`` = alpha * (projectile charge); may be non-integer charge (tests, effective q).
    """
    az = _LD(alpha_z)
    b = _LD(beta)
    if not (0 < b < 1) or not (0 < az < 1):
        raise ValueError("need 0 < beta < 1 and 0 < alpha z < 1")
    g = 1 / np.sqrt(1 - b * b)
    eta = az / b
    xi = eta / g
    k = np.arange(1, k_max + 2, dtype=_LD)  # 1 .. K+1 (|k|)
    rho = np.sqrt(k * k - az * az)
    arg_gamma = clgamma(rho + 1 + 1j * eta).imag  # same for +k and -k
    small = az * az / (k + rho)  # = |k| - rho_k, cancellation-free
    half_pi = _LD(np.pi) / 2
    # phi_k = (1/2) [arg(k - i xi) - arg(rho_k - i eta)]   (mod pi)
    phi_p = (np.arctan2(-xi, k) - np.arctan2(-eta, rho)) / 2
    phi_m = (np.arctan2(-xi, -k) - np.arctan2(-eta, rho)) / 2
    d_p = phi_p - arg_gamma + half_pi * small  # l - rho = k - rho
    d_m = phi_m - arg_gamma + half_pi * (small - 1)  # l - rho = (|k| - 1) - rho
    kk = k[:k_max]
    diff_p = np.zeros(k_max, dtype=_LD)
    diff_p[1:] = d_p[1:k_max] - d_p[: k_max - 1]  # d_k - d_{k-1}, k >= 2 (k = 1 has weight 0)
    diff_m = d_m[:k_max] - d_m[1 : k_max + 1]  # d_{-k} - d_{-k-1}
    t = (
        kk / (eta * eta)
        * ((kk - 1) / (2 * kk - 1) * np.sin(diff_p) ** 2
           + (kk + 1) / (2 * kk + 1) * np.sin(diff_m) ** 2)
        + kk / ((4 * kk * kk - 1) * (g * g * kk * kk + eta * eta))
        - 1 / kk
    )  # fmt: skip
    return t, float(b * b / 2), float(eta)


def _tail_power(p: int, k: int) -> float:
    """sum_{n > k} n^-p (p >= 2) by Euler-Maclaurin (4 terms; error O(k^(-p-5)))."""
    kk = float(k)
    return float(
        kk ** (1 - p) / (p - 1)
        - 0.5 * kk**-p
        + p * kk ** (-p - 1) / 12.0
        - p * (p + 1) * (p + 2) * kk ** (-p - 3) / 720.0
    )


def _sum_with_tail(t: np.ndarray, k: int) -> tuple[float, float, float]:
    """Explicit sum of t_1..t_k plus fitted tail. Returns (total, tail, c2)."""
    ks = np.array([k, k // 2, k // 4], dtype=np.float64)
    vals = np.array([float(t[int(x) - 1]) for x in ks])
    a = np.stack([ks**-2, ks**-3, ks**-4], axis=1)
    c2, c3, c4 = np.linalg.solve(a, vals)
    tail = c2 * _tail_power(2, k) + c3 * _tail_power(3, k) + c4 * _tail_power(4, k)
    return float(np.sum(t[:k], dtype=_LD)) + tail, tail, c2


@dataclass(frozen=True)
class LSResult:
    """Converged LS correction and its convergence record."""

    delta_l: float
    k_used: int
    change_vs_half_k: float
    tail: float
    tail_coeff_c2: float
    eta: float


def delta_l_ls(
    alpha_z: float, beta: float, *, k_start: int = 2000, k_limit: int = 256000, tol: float = 1e-11
) -> LSResult:
    """Point-nucleus Lindhard-Sorensen Delta L ([KHV] eq. 10), converged in k (see module doc)."""
    k = k_start
    t, b2, eta = ls_summand(alpha_z, beta, k)
    prev, _, _ = _sum_with_tail(t, k // 2)
    while True:
        cur, tail, c2 = _sum_with_tail(t, k)
        if abs(cur - prev) <= tol:
            return LSResult(float(cur + b2), k, float(cur - prev), float(tail), float(c2), float(eta))
        if 2 * k > k_limit:
            raise RuntimeError(f"LS sum not converged: |change| = {abs(cur - prev):.3e} at K={k}")
        prev = cur
        k *= 2
        t, _, _ = ls_summand(alpha_z, beta, k)


def bloch_ionmc(beta2: float, q: float) -> float:
    """Bloch term exactly as ionmc ``bloch_term``: -y^2 sum_{n<=2000} 1/(n(n^2+y^2)) + tail.

    y = q alpha / beta; tail approximation 0.5/(2000.5)^2 as in the repository code.
    """
    y2 = (q * ALPHA_FS) ** 2 / beta2
    n = np.arange(1, 2001, dtype=np.float64)
    s = float(np.sum(1.0 / (n * (n * n + y2)))) + 0.5 / (2000 + 0.5) ** 2
    return float(-y2 * s)


def bloch_exact(eta: float) -> float:
    """Converged Bloch term -eta^2 sum_k 1/(k(k^2+eta^2)) (reference for the tests)."""
    n = np.arange(1, 200001, dtype=np.float64)
    s = float(np.sum(1.0 / (n * (n * n + eta * eta)))) + 0.5 / 200000.5**2
    return float(-eta * eta * s)


def mott1_ionmc(beta: float, q: float) -> float:
    """First-order Mott term exactly as ionmc: 0.5 pi alpha beta q (Ahlen 1980)."""
    return 0.5 * math.pi * ALPHA_FS * beta * q


def effective_charge_ionmc(z: int, beta: float) -> float:
    """ionmc Barkas/Pierce-Blann effective charge z(1 - exp(-125 beta z^(-2/3))), z >= 2."""
    if z < 2:
        return float(z)
    return z * (1.0 - math.exp(-125.0 * beta * z ** (-2.0 / 3.0)))


def residual(z: int, e_mev_per_u: float, *, q: float | None = None, **ls_kwargs: float) -> float:
    """Delta L_LS - (Bloch + first-order Mott) for nucleus ``z`` at E [MeV/u].

    ``q``: charge used in LS, Bloch and Mott (default: bare charge z). ionmc applies Bloch and
    Mott with its effective charge; pass ``q=effective_charge_ionmc(z, beta)`` for that variant.
    """
    beta, _ = kinematics(z, e_mev_per_u)
    qq = float(z) if q is None else float(q)
    ls = delta_l_ls(ALPHA_FS * qq, beta, **ls_kwargs).delta_l  # type: ignore[arg-type]
    return float(ls - bloch_ionmc(beta * beta, qq) - mott1_ionmc(beta, qq))
