"""Tripathi light-system proton-nucleus reaction cross section (build-time, numpy only).

Source: R. K. Tripathi, F. A. Cucinotta, J. W. Wilson, "Universal parameterization of absorption
cross sections - light systems", NASA/TP-1999-209726, equations (3)-(15) and Table 2; the
formulas are written from the paper (the Geant4 transcription was read only to cross-check
constants). Proton projectile (``n(p)+X`` row of the paper: ``T1 = 23``)::

    sigma = pi r0^2 (A_P^1/3 + A_T^1/3 + dE)^2 (1 - R_c B / E_cm) X_m,      r0 = 1.1 fm
    B = 1.44 Z_P Z_T / R [MeV],   R = r_P + r_T + 1.2 (A_P^1/3 + A_T^1/3) / E_cm^1/3 [fm]
    r_i = 1.29 r_rms,i
    dE = 1.85 S + 0.16 S / E_cm^1/3 - C_E + 0.91 (A_T - 2 Z_T) Z_P / (A_T A_P)
    S = A_P^1/3 A_T^1/3 / (A_P^1/3 + A_T^1/3)
    C_E = D [1 - exp(-E / T1)] - 0.292 exp(-E / 792) cos(0.229 E^0.453)
    D = 1.85 + 0.16 / (1 + exp((500 - E) / 200))
    X_m = 1 - X1 exp(-E / (X1 S_L)),  X1 = 2.83 - 3.1e-2 A_T + 1.7e-4 A_T^2
    S_L = 1.2 + 1.6 [1 - exp(-E / 15)]

``E`` is the lab energy per nucleon [MeV] (= T for a proton), ``E_cm`` the total kinetic energy in
the p + target centre of mass [MeV] (as in the Geant4 transcription; it enters only the Coulomb
term and a weak ``E_cm^-1/3`` term). ``R_c = 1`` for targets heavier than Li. Approximations
declared: nuclear mass ``A u - Z m_e`` (no binding); rounded electron-scattering rms charge radii
(de Vries et al., ADNDT 36 (1987) 495) for the supported nuclides (a 5 % radius change alters sigma
by < 0.5 % at 100 MeV, tested); NOT cross-checked against the paper PDF beyond the constants above.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

U_MEV = 931.49410242
M_E_MEV = 0.51099895
M_P_MEV = 938.27208816
R_RMS_FM = {1: 0.875, 12: 2.472, 14: 2.540, 16: 2.730, 27: 3.061, 31: 3.187, 40: 3.478}
"""Rounded rms charge radii [fm] by mass number (proton 0.875 fm)."""
T1_PROTON = 23.0


def sigma_tripathi_light_p(
    e_mev: ArrayLike, z_t: int, a_t: int, radius_scale: float = 1.0
) -> NDArray[np.float64]:
    """Proton + (Z_t, A_t) reaction cross section [barn] at proton kinetic energy ``e_mev`` [MeV]
    (module docstring). Raises ``ValueError`` for a target without a tabulated radius and for
    non-positive energies; negative values of the formula are clamped to 0."""
    if a_t not in R_RMS_FM:
        raise ValueError(f"no rms radius for A_t={a_t}; supported {sorted(R_RMS_FM)}")
    e = np.asarray(e_mev, dtype=np.float64)
    if np.any(e <= 0.0):
        raise ValueError("energies must be positive")
    m_t = a_t * U_MEV - z_t * M_E_MEV
    e_lab = e + M_P_MEV
    s2 = M_P_MEV**2 + m_t**2 + 2.0 * e_lab * m_t
    e_cm = np.sqrt(s2) - m_t - M_P_MEV
    e_cm13 = np.cbrt(e_cm)
    ap13, at13 = 1.0, float(a_t) ** (1.0 / 3.0)
    s = ap13 * at13 / (ap13 + at13)
    d = 1.85 + 0.16 / (1.0 + np.exp((500.0 - e) / 200.0))
    c_e = d * (1.0 - np.exp(-e / T1_PROTON)) - 0.292 * np.exp(-e / 792.0) * np.cos(0.229 * e**0.453)
    de = 1.85 * s + 0.16 * s / e_cm13 - c_e + 0.91 * (a_t - 2 * z_t) * 1.0 / (a_t * 1.0)
    x1 = 2.83 - 3.1e-2 * a_t + 1.7e-4 * a_t**2
    s_l = 1.2 + 1.6 * (1.0 - np.exp(-e / 15.0))
    x_m = 1.0 - x1 * np.exp(-e / (x1 * s_l))
    radius = 1.29 * radius_scale * (R_RMS_FM[1] + R_RMS_FM[a_t]) + 1.2 * (ap13 + at13) / e_cm13
    b = 1.44 * 1.0 * z_t / radius
    xr = 1.1 * (ap13 + at13 + de)
    sigma_fm2 = np.pi * xr * xr * (1.0 - b / e_cm) * x_m
    return np.asarray(np.maximum(sigma_fm2 * 0.01, 0.0), dtype=np.float64)  # 1 fm^2 = 0.01 b


def extension_factor(
    e_mev: ArrayLike, z_t: int, a_t: int, e_anchor: float = 150.0
) -> NDArray[np.float64]:
    """``sigma_TL(E) / sigma_TL(E_anchor)``; exactly 1.0 at ``E == E_anchor`` (decision 0041)."""
    e = np.asarray(e_mev, dtype=np.float64)
    ratio = sigma_tripathi_light_p(e, z_t, a_t) / float(sigma_tripathi_light_p(e_anchor, z_t, a_t))
    return np.where(e == e_anchor, 1.0, ratio)
