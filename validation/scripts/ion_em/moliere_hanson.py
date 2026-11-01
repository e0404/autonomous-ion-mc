"""Moliere-Hanson multiple-scattering reference for protons and light ions (V3-008, check P1).

Pure functions, numpy only (the optional report driver ``main`` imports ``ionmc`` lazily for the
range-energy relation). Every quantity is the PROJECTED Gaussian width theta_0 [rad] of the whole
slab (Hanson width), the convention of ``ionmc.physics.em`` (projected scattering power) and of
the column theta_Hanson of Gottschalk, Med. Phys. 37 (2010) 352 (arXiv:0908.1413).

Theory (Moliere 1948; Bethe, Phys. Rev. 89 (1953) 1256; thick-slab and compound rules as stated
by Gottschalk et al., NIM B 74 (1993) 467 and Gottschalk arXiv:1804.00022 eqs. chicsq-Hanson):

* characteristic angle (space angle), Bethe's Z(Z+1) form, compound by mass fractions w_i::

      chi_c^2 = K_C z^2 sum_i w_i Z_i (Z_i + 1) / A_i  *  int_0^t dt' / (pv)^2
      K_C = 4 pi N_A r_e^2 (m_e c^2)^2 = 0.156915 MeV^2 cm^2 / g

* screening angle per element (Moliere; Gottschalk 2018 eqs. chiasq, c1, c2)::

      chi_a,i^2 = chi_0,i^2 (1.13 + 3.76 alpha_i'^2),
      chi_0,i = alpha m_e c^2 Z_i^(1/3) / (0.885 pc),
      alpha_i' = z Z_i alpha / beta

  compound: ln chi_a^2 = sum_i W_i ln chi_a,i^2 with W_i proportional to w_i Z_i (Z_i + 1) / A_i;
  thick slab: ln chi_a^2 averaged along the slab with weight 1/(pv)^2 (Moliere's prescription);
* b = ln(chi_c^2 / (1.167 chi_a^2)), 1.167 = exp(2 C_E - 1); B - ln B = b;
* theta_0 = chi_c sqrt(B - 1.2) / sqrt(2)   (Hanson et al. 1951; Gottschalk 2010 eq. thetaHanson).

Electron treatment (P1b; ``electrons=`` switch, default ``"fano"``):

* ``"fano"`` (the reference): Moliere's Z^2 form (Z^2 in chi_c^2 and in the compound weights W_i)
  plus Fano's heavy-particle correction, Fano, Phys. Rev. 93 (1954) 117, abstract (b): "For
  incident heavy particles, leave Z^2 unaltered but increase b by
  Z^-1 {ln[1130 Z^-4/3 (c^2/v^2 - 1)^-1] - u_in - v^2/(2c^2)}". Gottschalk uses this variant for
  theta_Hanson (Makarova, Gottschalk, Sauerwein, arXiv:1610.01279 sec. 2: "Moliere's rather than
  Bethe's form of the theory (Z^2 rather than Z(Z+1)) and using the Fano correction"). -u_in = 3.6
  for hydrogen (Fano) and 5.0 otherwise ("as suggested by Gottschalk", MICE, arXiv:2209.10251).
  Compound/thick-slab rule (derived here, not sourced): delta_b = sum_i w_i Z_i/A_i D_i /
  sum_i w_i Z_i^2/A_i, with D_i averaged along the slab with weight 1/(pv)^2.
* ``"bethe"`` (P1 record): Bethe's Z(Z+1) in chi_c^2 and W_i, no Fano term.

Point nucleus.

Candidate (d) correction (frozen formulas, see ``f_scr``, ``n_omega``, ``c_omega``)::

      f_dM,ion = f_dM(pv*, p1v1*) + c_Omega * ln( (M / m_p) / F_scr(beta, z; material) )
      ln F_scr = sum_i W_i [ ln(1.13 + 3.76 (z Z_i alpha/beta)^2)
                             - ln(1.13 + 3.76 (Z_i alpha/beta)^2) ]
      c_Omega  = (c_2 + c_3 L_1) / (ln 10 * n_Omega),   L_1 = lg(1 - (pv/p1v1)^2)
      n_Omega  = d ln R_p / d ln(pv_p) - d ln beta^2 / d ln(pv_p)
                 (proton, at the velocity-scaled energy E* = E m_p / M)

Units: MeV, g/cm^2, rad.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

import numpy as np
from numpy.typing import NDArray

# --- constants (sources) -----------------------------------------------------------------------
ALPHA = 1.0 / 137.035999084  # CODATA 2018 fine-structure constant
R_E_CM = 2.8179403262e-13  # CODATA 2018 classical electron radius [cm]
ME_MEV = 0.51099895000  # CODATA 2018 electron rest energy [MeV]
N_A = 6.02214076e23  # SI 2019 exact [1/mol]
M_P_MEV = 938.27208816  # CODATA 2018 proton rest energy [MeV]
EULER_GAMMA = 0.5772156649015329
K_CHI_C = 4.0 * math.pi * N_A * R_E_CM**2 * ME_MEV**2
"""4 pi N_A r_e^2 (m_e c^2)^2 [MeV^2 cm^2/g] = 0.156915 (Gottschalk 1993; 2018 c3)."""
CHI_0_MEV = ALPHA * ME_MEV / 0.885
"""alpha m_e c^2 / 0.885 [MeV]: chi_0 = CHI_0_MEV Z^(1/3) / pc (Thomas-Fermi 0.885 a_0 Z^-1/3)."""
SCREEN_A, SCREEN_B = 1.13, 3.76
"""Moliere screening factor 1.13 + 3.76 alpha'^2 (Bethe 1953; Gottschalk 2018 eq. chiasq)."""
OMEGA_FACTOR = math.exp(2.0 * EULER_GAMMA - 1.0)
"""1.167 = exp(2 C_E - 1) in b = ln(chi_c^2 / (1.167 chi_a^2))."""
HANSON_SHIFT = 1.2
"""theta_0 = chi_c sqrt(B - 1.2)/sqrt(2) (Hanson et al., Phys. Rev. 84 (1951) 634)."""
E_S_MEV = 15.0
U_IN_NEG = {1: 3.6}
"""-u_in of Fano's correction: 3.6 for hydrogen (Fano 1954), else U_IN_NEG_DEFAULT."""
U_IN_NEG_DEFAULT = 5.0
"""-u_in = 5.0 for Z > 1 (Gottschalk's value as quoted by MICE, arXiv:2209.10251)."""
FANO_CONST = 1130.0
ELECTRONS = ("fano", "bethe")
FDM = (0.5244, 0.1975, 0.2320, -0.0098)
"""Gottschalk 2010 f_dM coefficients, identical to ionmc.physics.em.FDM_COEFFICIENTS."""

Element = tuple[int, float, float]
"""(Z, A [g/mol], mass fraction w)."""

WATER: tuple[Element, ...] = (
    (1, 1.00794, 2 * 1.00794 / (2 * 1.00794 + 15.9994)),
    (8, 15.9994, 15.9994 / (2 * 1.00794 + 15.9994)),
)
"""H2O with the atomic weights of ``ionmc.materials`` (H 1.00794, O 15.9994)."""


# --- kinematics --------------------------------------------------------------------------------
def pv_mev(t_mev: NDArray[np.float64] | float, m_mev: float) -> NDArray[np.float64]:
    """p v = T (T + 2M) / (T + M) [MeV]."""
    t = np.asarray(t_mev, dtype=np.float64)
    return t * (t + 2.0 * m_mev) / (t + m_mev)


def pc_mev(t_mev: NDArray[np.float64] | float, m_mev: float) -> NDArray[np.float64]:
    """p c = sqrt(T (T + 2M)) [MeV]."""
    t = np.asarray(t_mev, dtype=np.float64)
    return np.sqrt(t * (t + 2.0 * m_mev))


def beta2(t_mev: NDArray[np.float64] | float, m_mev: float) -> NDArray[np.float64]:
    """beta^2 = T (T + 2M) / (T + M)^2."""
    t = np.asarray(t_mev, dtype=np.float64)
    return t * (t + 2.0 * m_mev) / (t + m_mev) ** 2


# --- Moliere quantities ------------------------------------------------------------------------
def zz1_sum(elements: Sequence[Element]) -> float:
    """sum_i w_i Z_i (Z_i + 1) / A_i [mol/g]."""
    return math.fsum(w * zt * (zt + 1) / a for zt, a, w in elements)


def z2_sum(elements: Sequence[Element]) -> float:
    """sum_i w_i Z_i^2 / A_i [mol/g] (Moliere's Z^2 form, used with the Fano correction)."""
    return math.fsum(w * zt * zt / a for zt, a, w in elements)


def _check_electrons(electrons: str) -> None:
    if electrons not in ELECTRONS:
        raise ValueError(f"electrons must be one of {ELECTRONS}, got {electrons!r}")


def screening_weights(elements: Sequence[Element], electrons: str = "fano") -> NDArray[np.float64]:
    """W_i normalised: w_i Z_i^2/A_i (``"fano"``) or w_i Z_i (Z_i + 1)/A_i (``"bethe"``)."""
    _check_electrons(electrons)
    extra = 0 if electrons == "fano" else 1
    raw = np.array([w * zt * (zt + extra) / a for zt, a, w in elements], dtype=np.float64)
    out: NDArray[np.float64] = raw / float(raw.sum())
    return out


def ln_screening_factor(z: float, zt: int, b2: NDArray[np.float64]) -> NDArray[np.float64]:
    """ln(1.13 + 3.76 alpha'^2), alpha' = z Z alpha / beta."""
    return np.log(SCREEN_A + SCREEN_B * (z * zt * ALPHA) ** 2 / b2)


def ln_chi_a2(
    z: float,
    elements: Sequence[Element],
    t_mev: NDArray[np.float64],
    m_mev: float,
    electrons: str = "fano",
) -> NDArray[np.float64]:
    """Compound ln chi_a^2 at energies ``t_mev`` (weights from :func:`screening_weights`)."""
    pc = pc_mev(t_mev, m_mev)
    b2 = beta2(t_mev, m_mev)
    out = np.zeros_like(pc)
    for wgt, (zt, _a, _w) in zip(screening_weights(elements, electrons), elements, strict=True):
        out += wgt * (
            2.0 * np.log(CHI_0_MEV * zt ** (1.0 / 3.0) / pc) + ln_screening_factor(z, zt, b2)
        )
    return out


def fano_d(zt: int, b2: NDArray[np.float64]) -> NDArray[np.float64]:
    """Fano's D = ln[1130 Z^-4/3 beta^2/(1 - beta^2)] - u_in - beta^2/2 (Fano 1954, (b))."""
    neg_u = U_IN_NEG.get(zt, U_IN_NEG_DEFAULT)
    b2a = np.asarray(b2, dtype=np.float64)
    out: NDArray[np.float64] = (
        np.log(FANO_CONST * zt ** (-4.0 / 3.0) * b2a / (1.0 - b2a)) + neg_u - 0.5 * b2a
    )
    return out


def fano_delta_b(elements: Sequence[Element], b2: NDArray[np.float64]) -> NDArray[np.float64]:
    """Increase of Moliere's b: sum_i w_i Z_i/A_i D_i / sum_i w_i Z_i^2/A_i (= D/Z for an
    element; the compound weighting electrons ~ Z, nucleus ~ Z^2 is derived, not sourced)."""
    num = np.zeros_like(np.asarray(b2, dtype=np.float64))
    for zt, a, w in elements:
        num += w * zt / a * fano_d(zt, b2)
    return num / z2_sum(elements)


def solve_b(b: float) -> float:
    """Root B > 1 of B - ln B = b (Newton from B = b + ln b); requires b > 1."""
    if not b > 1.0:
        raise ValueError(f"Moliere b = {b} <= 1: slab too thin for Moliere theory")
    x = b + math.log(b)
    for _ in range(60):
        step = (x - math.log(x) - b) / (1.0 - 1.0 / x)
        x -= step
        if abs(step) < 1e-15 * x:
            break
    return x


def _gauss_nodes(t_g_cm2: float, n: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Gauss-Legendre nodes on [0, t] with t' = t (1 - (1-s)^2) (clusters nodes at the slab end,
    where 1/(pv)^2 grows fastest)."""
    gx, gw = np.polynomial.legendre.leggauss(n)
    s = 0.5 * (gx + 1.0)
    nodes = t_g_cm2 * (1.0 - (1.0 - s) ** 2)
    weights = t_g_cm2 * 2.0 * (1.0 - s) * 0.5 * gw
    return nodes, weights


def moliere_slab(
    z: float,
    m_mev: float,
    elements: Sequence[Element],
    t_g_cm2: float,
    energy_at: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    n: int = 400,
    electrons: str = "fano",
) -> dict[str, float]:
    """Moliere quantities of a homogeneous slab of mass thickness ``t_g_cm2``.

    ``energy_at(t')`` returns the CSDA kinetic energy [MeV] of the projectile (charge ``z``,
    rest energy ``m_mev``) at mass depth t' [g/cm^2]. Returns chi_c^2 [rad^2], ln chi_a^2
    (slab average), b, B and theta_0 [rad] (projected Hanson width).
    """
    if not t_g_cm2 > 0.0:
        raise ValueError("slab thickness must be positive")
    nodes, weights = _gauss_nodes(t_g_cm2, n)
    energies = energy_at(nodes)
    inv_pv2 = 1.0 / pv_mev(energies, m_mev) ** 2
    integral = float(np.sum(weights * inv_pv2))
    _check_electrons(electrons)
    zsum = z2_sum(elements) if electrons == "fano" else zz1_sum(elements)
    chi_c2 = K_CHI_C * z * z * zsum * integral
    lna = (
        float(np.sum(weights * inv_pv2 * ln_chi_a2(z, elements, energies, m_mev, electrons)))
        / integral
    )
    delta_b = 0.0
    if electrons == "fano":
        db = fano_delta_b(elements, beta2(energies, m_mev))
        delta_b = float(np.sum(weights * inv_pv2 * db)) / integral
    b = math.log(chi_c2) - lna - math.log(OMEGA_FACTOR) + delta_b
    big_b = solve_b(b)
    theta0 = math.sqrt(chi_c2 * (big_b - HANSON_SHIFT) / 2.0)
    return {
        "chi_c2": chi_c2,
        "ln_chi_a2": lna,
        "delta_b": delta_b,
        "b": b,
        "B": big_b,
        "theta0": theta0,
    }


def theta0_hanson(
    z: float,
    m_mev: float,
    elements: Sequence[Element],
    t_g_cm2: float,
    energy_at: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    n: int = 400,
    electrons: str = "fano",
) -> float:
    """Projected Moliere-Hanson width theta_0 [rad] of the slab (see :func:`moliere_slab`)."""
    return moliere_slab(z, m_mev, elements, t_g_cm2, energy_at, n, electrons)["theta0"]


# --- f_dM and the MCS candidates (numpy re-statement of ionmc.physics.em) ------------------------
def inverse_scattering_length(elements: Sequence[Element]) -> float:
    """1/X_S [cm^2/g] (Gottschalk 2010 eq. XS, natural log, Bragg-additive)."""
    pref = ALPHA * N_A * R_E_CM**2
    return math.fsum(
        w * pref * zt * zt / a * (2.0 * math.log(33219.0 * (a * zt) ** (-1.0 / 3.0)) - 1.0)
        for zt, a, w in elements
    )


def fdm(pv: NDArray[np.float64], p1v1: float) -> NDArray[np.float64]:
    """Gottschalk's f_dM(pv, p1v1), clamped at 0 (as ``ionmc.physics.em.fdm``)."""
    om = 1.0 - (pv / p1v1) ** 2
    out = np.zeros_like(pv)
    ok = om > 0.0
    l1 = np.log10(om[ok])
    l2 = np.log10(pv[ok])
    out[ok] = FDM[0] + FDM[1] * l1 + FDM[2] * l2 + FDM[3] * l1 * l2
    return np.maximum(out, 0.0)


def f_scr(
    z: float, b2: NDArray[np.float64], elements: Sequence[Element], electrons: str = "fano"
) -> NDArray[np.float64]:
    """F_scr(beta, z; material): ratio of the Moliere screening factors of the ion and of a proton
    at the same velocity, geometric mean with the reference's weights W_i (frozen formula,
    candidate (d)). Fano's delta_b depends on beta only and cancels at equal velocity."""
    out = np.zeros_like(np.asarray(b2, dtype=np.float64))
    for wgt, (zt, _a, _w) in zip(screening_weights(elements, electrons), elements, strict=True):
        out += wgt * (ln_screening_factor(z, zt, b2) - ln_screening_factor(1.0, zt, b2))
    return np.exp(out)


def n_omega(
    e_p_mev: NDArray[np.float64],
    range_p: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    stopping_p: Callable[[NDArray[np.float64]], NDArray[np.float64]],
) -> NDArray[np.float64]:
    """n_Omega = d ln(R_p / beta^2) / d ln(pv) of a proton at kinetic energy ``e_p_mev``
    (analytic: d ln R/d ln E = E / (S R); kinematic derivatives exact). ``range_p`` [g/cm^2] and
    ``stopping_p`` [MeV cm^2/g] are the proton tables of the same material."""
    e = np.asarray(e_p_mev, dtype=np.float64)
    m = M_P_MEV
    dlnr = e / (stopping_p(e) * range_p(e))
    dlnpv = 1.0 + e / (e + 2.0 * m) - e / (e + m)
    dlnb2 = 1.0 + e / (e + 2.0 * m) - 2.0 * e / (e + m)
    return (dlnr - dlnb2) / dlnpv


def c_omega(
    pv_star: NDArray[np.float64], p1v1_star: float, n_om: NDArray[np.float64]
) -> NDArray[np.float64]:
    """c_Omega = (c_2 + c_3 L_1) / (ln 10 n_Omega): sensitivity of f_dM to ln Omega at fixed x/R
    (frozen formula, candidate (d))."""
    om = np.clip(1.0 - (pv_star / p1v1_star) ** 2, 1e-300, None)
    l1 = np.log10(om)
    return (FDM[2] + FDM[3] * l1) / (math.log(10.0) * n_om)


def theta0_candidate(
    candidate: str,
    z: float,
    m_mev: float,
    elements: Sequence[Element],
    t_g_cm2: float,
    e0_mev: float,
    energy_at: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    proton_range: Callable[[NDArray[np.float64]], NDArray[np.float64]] | None = None,
    proton_stopping: Callable[[NDArray[np.float64]], NDArray[np.float64]] | None = None,
    n: int = 400,
    electrons: str = "fano",
) -> float:
    """theta_0 [rad] of the slab from integrating a candidate scattering power along the CSDA path.

    ``a``: f_dM(pv*, p1v1*), pv* = pv m_p/M;
    ``b``: f_dM(pv, p1v1) with the ion's own pv;
    ``d``: ``a`` plus c_Omega ln((M/m_p)/F_scr) (needs the proton tables);
    ``e`` (informative, Gottschalk arXiv:1204.4470 eq. fdMi): f_dM at pv f_I,
    f_I = (z sqrt(M/m_p))^-1.018;
    ``c`` (slab proxy, NOT the differential Kanematsu form): Lynch-Dahl/PDG Highland with z^2 in
    the logarithm (see ``theta0_highland_z2``).
    """
    if candidate == "c":
        return theta0_highland_z2(z, m_mev, elements, t_g_cm2, e0_mev, energy_at, n)
    nodes, weights = _gauss_nodes(t_g_cm2, n)
    e = energy_at(nodes)
    pv = pv_mev(e, m_mev)
    p1v1 = float(pv_mev(e0_mev, m_mev))
    ratio = M_P_MEV / m_mev
    if candidate == "a":
        f = fdm(pv * ratio, p1v1 * ratio)
    elif candidate == "b":
        f = fdm(pv, p1v1)
    elif candidate == "e":
        fi = (z * math.sqrt(m_mev / M_P_MEV)) ** -1.018
        f = fdm(pv * fi, p1v1 * fi)
    elif candidate == "d":
        if proton_range is None or proton_stopping is None:
            raise ValueError("candidate d needs the proton range and stopping functions")
        e_star = e * ratio
        n_om = n_omega(e_star, proton_range, proton_stopping)
        corr = c_omega(pv * ratio, p1v1 * ratio, n_om) * np.log(
            (m_mev / M_P_MEV) / f_scr(z, beta2(e, m_mev), elements, electrons)
        )
        f = np.maximum(fdm(pv * ratio, p1v1 * ratio) + corr, 0.0)
    else:
        raise ValueError(f"unknown candidate {candidate!r}")
    power = f * (z * E_S_MEV / pv) ** 2 * inverse_scattering_length(elements)
    return math.sqrt(float(np.sum(weights * power)))


_TSAI = {1: (5.31, 6.144), 2: (4.79, 5.621), 3: (4.74, 5.805), 4: (4.71, 5.924)}


def radiation_length(elements: Sequence[Element]) -> float:
    """X_0 [g/cm^2], Tsai's formula as given by the PDG (as ``ionmc.transport.mcs_checks``)."""
    inv = 0.0
    for zt, a, w in elements:
        az = ALPHA * zt
        coul = (
            az
            * az
            * (1 / (1 + az * az) + 0.20206 - 0.0369 * az**2 + 0.0083 * az**4 - 0.002 * az**6)
        )
        lrad, lradp = _TSAI.get(
            zt, (math.log(184.15 * zt ** (-1 / 3)), math.log(1194.0 * zt ** (-2 / 3)))
        )
        inv += w * 4 * ALPHA * R_E_CM**2 * N_A / a * (zt * zt * (lrad - coul) + zt * lradp)
    return 1.0 / inv


def theta0_highland_z2(
    z: float,
    m_mev: float,
    elements: Sequence[Element],
    t_g_cm2: float,
    e0_mev: float,
    energy_at: Callable[[NDArray[np.float64]], NDArray[np.float64]],
    n: int = 400,
) -> float:
    """Candidate (c) slab proxy: PDG/Lynch-Dahl Highland with z^2 in the logarithm, thick-slab
    form with int dt/(X_0 (pv)^2) (Gottschalk 1993 generalisation) and the incident beta_1."""
    x0 = radiation_length(elements)
    nodes, weights = _gauss_nodes(t_g_cm2, n)
    pv = pv_mev(energy_at(nodes), m_mev)
    integral = float(np.sum(weights / (x0 * pv**2)))
    b1 = float(beta2(e0_mev, m_mev))
    return 13.6 * z * math.sqrt(integral) * (1.0 + 0.038 * math.log(t_g_cm2 * z * z / (x0 * b1)))


# --- report driver (imports ionmc lazily) --------------------------------------------------------
def _ionmc_path(
    projectile_name: str, e_per_u: float
) -> tuple[object, float, Callable[..., NDArray[np.float64]]]:
    from ionmc.materials import WATER as W  # noqa: PLC0415
    from ionmc.physics.projectiles import PROJECTILES  # noqa: PLC0415
    from ionmc.physics.stopping import BetheStoppingSource  # noqa: PLC0415

    proj = PROJECTILES[projectile_name]
    table = BetheStoppingSource().table(W, proj)
    r1 = float(table.range_at(e_per_u))

    def energy_at(t: NDArray[np.float64]) -> NDArray[np.float64]:
        rem = np.maximum(r1 - t, float(table.csda_range_g_cm2[0]))
        return np.asarray(table.energy_from_range(rem)) * proj.a

    return proj, r1, energy_at


def ion_table(
    fractions: Sequence[float] = (0.1, 0.5, 0.9), electrons: str = "fano"
) -> list[dict[str, float | str]]:
    """Report-only: theta_candidate / theta_MH for He, C, O in water (I11 planning grid)."""
    from ionmc.materials import WATER as W  # noqa: PLC0415
    from ionmc.physics.projectiles import PROTON  # noqa: PLC0415
    from ionmc.physics.stopping import BetheStoppingSource  # noqa: PLC0415

    ptab = BetheStoppingSource().table(W, PROTON)

    def rp(e: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.asarray(ptab.range_at(np.clip(e, ptab.energy_per_u[0], ptab.energy_per_u[-1])))

    def sp(e: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.asarray(ptab.stopping_at(np.clip(e, ptab.energy_per_u[0], ptab.energy_per_u[-1])))

    rows: list[dict[str, float | str]] = []
    for name in ("alpha", "carbon12", "oxygen16"):
        for e_u in (50.0, 100.0, 150.0, 200.0, 300.0, 430.0):
            proj, r1, energy_at = _ionmc_path(name, e_u)
            for frac in fractions:
                t = frac * r1
                pz, pm = float(proj.z), proj.mass_mev  # type: ignore[attr-defined]
                ref = theta0_hanson(pz, pm, WATER, t, energy_at, electrons=electrons)
                row: dict[str, float | str] = {"ion": name, "e_u": e_u, "x_over_r": frac}
                for cand in ("a", "b", "c", "d", "e"):
                    th = theta0_candidate(
                        cand,
                        pz,
                        pm,
                        WATER,
                        t,
                        e_u * proj.a,  # type: ignore[attr-defined]
                        energy_at,
                        rp,
                        sp,
                        electrons=electrons,
                    )
                    row[cand] = th / ref - 1.0
                rows.append(row)
    return rows


def main() -> None:
    for r in ion_table():
        print(
            f"{r['ion']:9s} {r['e_u']:6.0f} {r['x_over_r']:4.2f} "
            + " ".join(f"{c}={100 * float(r[c]):+6.2f}%" for c in ("a", "b", "c", "d", "e"))
        )


if __name__ == "__main__":
    main()
