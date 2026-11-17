"""Lindhard-Sorensen (LS) stopping-number correction for a point nucleus (V3-008, option ``ls``).

Port of the validated reference ``validation/scripts/ion_em/lindhard_sorensen.py`` (numpy only,
float64 inputs, ``longdouble`` internally). The correction is returned in the ionmc convention
of ``physics/stopping.py`` (the "4 pi r_e^2" convention, in which the first-order Mott term is
``(pi/2) alpha beta z``)::

    S/rho = K q^2 (Z/A) beta^-2 [ L0 + q L1 + dL ],  dL = dL_LS(q_ls, beta)

Equations: Kats, Halenka, Voskresenskaya, arXiv:2008.04819v2, eqs. (9)-(10) and Rakhno, Mokhov,
Striganov, arXiv:physics/0505043 (Lindhard and Sorensen, Phys. Rev. A 53, 2443 (1996), was not
read directly)::

    dL_LS = beta^2/2 + sum_{k>=1} [ k/eta^2 (k-1)/(2k-1) sin^2(d_k - d_{k-1})
                                  + k/eta^2 (k+1)/(2k+1) sin^2(d_{-k} - d_{-k-1})
                                  + k / ((4k^2-1)(gamma^2 k^2 + eta^2)) - 1/k ]
    eta = alpha z / beta,  d_k = phi_k - arg Gamma(rho_k + 1 + i eta) + (pi/2)(l - rho_k),
    exp(2 i phi_k) = (k - i eta/gamma)/(rho_k - i eta),  rho_k = sqrt(k^2 - (alpha z)^2),
    l = k (k > 0), l = -k - 1 (k < 0).

The sum contains the Bloch term exactly and the first-order Mott term as its 1/k^2 tail; the
point-nucleus LS correction therefore stands in for ``bloch_term + 0.5 pi alpha beta q``.
Finite nuclear size is neglected (irrelevant for gamma <= 1.5, z <= 8).

Numerics (identical to the reference): phases in ``longdouble``; complex log-Gamma by upward
recurrence to Re >= 16 and an 8-term Stirling series; ``l - rho_k`` as ``(alpha z)^2/(|k|+rho_k)``;
the sum is taken to K = 2000 and the tail ``sum_{k>K}`` is added from a fit
``c2/k^2 + c3/k^3 + c4/k^4`` to the summand at K, K/2, K/4 (Euler-Maclaurin power sums). K is
doubled until the tail-corrected sums at K and K/2 agree to 1e-11 (absolute in dL), otherwise
RuntimeError. Energies are processed in batches (one 2-D array per batch); a point that does not
converge at K = 2000 is redone alone with doubling K.

Cost (single thread, measured on the development host): about 1.4 ms per point, i.e. about 1.4 s
for a 1000-point grid at z = 8 (below the 10 s budget, so no caching or reduced K is used).

Physics status. The option ``BetheOptions.ls`` of ``physics/stopping.py`` exists because the
"mixed" charge convention (LS at the bare charge z, Bloch and Mott at the effective charge q; the
Geant4 indexing) met the pre-declared P3 trigger of decision 0043 sec. 1, while the ionmc
convention (every correction at q) did not. Which convention is closer to the stopping of a
partially stripped ion is unsettled (it needs a screened-ion LS calculation or data). The default
``off`` is the orchestrator's amendment decision. See decision 0043 sec. 1 and
``validation/plans/v3-008/p3-lindhard-sorensen.md``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

ALPHA_FS = 1.0 / 137.035999084

_LD = np.longdouble
_CLD = np.clongdouble
# Stirling coefficients B_2n / (2n (2n-1)), n = 1..8.
_STIRLING = ((1, 12), (-1, 360), (1, 1260), (-1, 1680), (1, 1188), (-691, 360360), (1, 156),
             (-3617, 122400))  # fmt: skip
_HALF_LN_2PI = _LD("0.918938533204672741780329736405617639861")

K_START = 2000
K_LIMIT = 256000
TOL = 1.0e-11
_BATCH = 16


def clgamma(z: ArrayLike) -> NDArray[np.clongdouble]:
    """Complex log-Gamma (clongdouble) for Re z > 0, imaginary part correct modulo 2 pi."""
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
    out: NDArray[np.clongdouble] = r - acc
    return out


def _summand(
    az: NDArray[np.longdouble], b: NDArray[np.longdouble], k_max: int
) -> NDArray[np.longdouble]:
    """Summand t_k, k = 1..k_max, shape (n, k_max), for per-row ``alpha z`` and ``beta``."""
    az = az[:, None]
    b = b[:, None]
    g = 1 / np.sqrt(1 - b * b)
    eta = az / b
    xi = eta / g
    k = np.arange(1, k_max + 2, dtype=_LD)[None, :]
    rho = np.sqrt(k * k - az * az)
    arg_gamma = clgamma(rho + 1 + 1j * eta).imag
    small = az * az / (k + rho)
    half_pi = _LD(np.pi) / 2
    phi_p = (np.arctan2(-xi, k) - np.arctan2(-eta, rho)) / 2
    phi_m = (np.arctan2(-xi, -k) - np.arctan2(-eta, rho)) / 2
    d_p = phi_p - arg_gamma + half_pi * small
    d_m = phi_m - arg_gamma + half_pi * (small - 1)
    kk = k[:, :k_max]
    diff_p = np.zeros((az.shape[0], k_max), dtype=_LD)
    diff_p[:, 1:] = d_p[:, 1:k_max] - d_p[:, : k_max - 1]
    diff_m = d_m[:, :k_max] - d_m[:, 1 : k_max + 1]
    t: NDArray[np.longdouble] = (
        kk / (eta * eta)
        * ((kk - 1) / (2 * kk - 1) * np.sin(diff_p) ** 2
           + (kk + 1) / (2 * kk + 1) * np.sin(diff_m) ** 2)
        + kk / ((4 * kk * kk - 1) * (g * g * kk * kk + eta * eta))
        - 1 / kk
    )  # fmt: skip
    return t


def _tail_power(p: int, k: int) -> float:
    """sum_{n > k} n^-p (p >= 2) by Euler-Maclaurin (4 terms)."""
    kk = float(k)
    return float(
        kk ** (1 - p) / (p - 1)
        - 0.5 * kk**-p
        + p * kk ** (-p - 1) / 12.0
        - p * (p + 1) * (p + 2) * kk ** (-p - 3) / 720.0
    )


def _sum_with_tail(t: NDArray[np.longdouble], k: int) -> NDArray[np.float64]:
    """Row-wise explicit sum of t_1..t_k plus the fitted tail (float64, one value per row)."""
    ks = np.array([k, k // 2, k // 4], dtype=np.float64)
    vals = np.stack([t[:, int(x) - 1].astype(np.float64) for x in ks], axis=1)  # (n, 3)
    a = np.stack([ks**-2, ks**-3, ks**-4], axis=1)
    coef = np.linalg.solve(a[None, :, :], vals[:, :, None])[:, :, 0]  # (n, 3)
    tail = (
        coef[:, 0] * _tail_power(2, k)
        + coef[:, 1] * _tail_power(3, k)
        + coef[:, 2] * _tail_power(4, k)
    )
    head = np.sum(t[:, :k], axis=1, dtype=_LD).astype(np.float64)
    return head + tail


def _converge_one(az: float, beta: float) -> float:
    k = K_START
    azl, bl = np.array([az], dtype=_LD), np.array([beta], dtype=_LD)
    t = _summand(azl, bl, k)
    prev = _sum_with_tail(t, k // 2)[0]
    while True:
        cur = _sum_with_tail(t, k)[0]
        if abs(cur - prev) <= TOL:
            return float(cur + beta * beta / 2)
        if 2 * k > K_LIMIT:
            raise RuntimeError(f"LS sum not converged: |change| = {abs(cur - prev):.3e} at K={k}")
        prev = cur
        k *= 2
        t = _summand(azl, bl, k)


def delta_l_ls(charge: ArrayLike, beta: ArrayLike) -> NDArray[np.float64]:
    """Point-nucleus Lindhard-Sorensen ``dL_LS`` (ionmc convention) for ``alpha*charge`` and beta.

    ``charge`` (projectile charge in units of e, not necessarily an integer) and ``beta`` are
    broadcast to one dimension; float64 result of the same length. Requires
    ``0 < beta < 1`` and ``0 < alpha*charge < 1``.
    """
    q, b = np.broadcast_arrays(
        np.atleast_1d(np.asarray(charge, dtype=np.float64)),
        np.atleast_1d(np.asarray(beta, dtype=np.float64)),
    )
    q = q.reshape(-1)
    b = b.reshape(-1)
    az = ALPHA_FS * q
    if not (np.all((b > 0) & (b < 1)) and np.all((az > 0) & (az < 1))):
        raise ValueError("need 0 < beta < 1 and 0 < alpha z < 1")
    out = np.empty(q.size, dtype=np.float64)
    for s in range(0, q.size, _BATCH):
        sl = slice(s, s + _BATCH)
        azl, bl = az[sl].astype(_LD), b[sl].astype(_LD)
        t = _summand(azl, bl, K_START)
        prev = _sum_with_tail(t, K_START // 2)
        cur = _sum_with_tail(t, K_START)
        ok = np.abs(cur - prev) <= TOL
        res = cur + (b[sl] * b[sl] / 2)
        for j in np.flatnonzero(~ok):
            res[j] = _converge_one(float(az[sl][j]), float(b[sl][j]))
        out[sl] = res
    return out
