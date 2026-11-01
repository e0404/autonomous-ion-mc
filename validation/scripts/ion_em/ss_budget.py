"""V5-SS branch-B budget: single-scattering widths of the Geant4 ion path in 1 mm water (V3-008).

Pre-freeze planning calculation for plan row V5-SS and ``validation/plans/v3-008/
p10-single-scattering-availability.md`` section 7.5.1 (budget, discrimination, cost) and
section 7.8 (this module). Numpy only; no ionmc, TOPAS or Geant4 output is read; nothing is
computed or written at import. Every number it prints is SEEN pre-freeze, a planning value and
not qualification evidence. It is the committed, tested form of the planner's session-scratch
scripts ``ss_budget.py`` (e50e3450...), ``ss_check.py`` (ff8be412...), ``ss_cost.py``
(2e378188...) and ``ss_signal.py`` (1dc8b5eb...), with the same constants and numerics.

Method (P10 section 7.5.1)
--------------------------
* Small-angle screened-Rutherford single scattering. For a species i with characteristic angle
  chi_c,i^2 (space-angle convention of Moliere/Bethe) and screening angle a_i, the exponent of the
  characteristic function of the projected angle (Moliere; Bethe, Phys. Rev. 89 (1953) 1256,
  eqs. for Omega(eta)) is::

      Omega(eta) = sum_i 2 chi_c,i^2 int_0^theta_max,i theta (1 - J0(eta theta))
                                                   / (theta^2 + a_i^2)^2 dtheta

  For theta_max = infinity the integral is closed form (Gradshteyn-Ryzhik 6.532.4)::

      Omega_i(eta) = chi_c,i^2 (1 - eta a_i K1(eta a_i)) / a_i^2

  A finite theta_max (electron scattering) subtracts eta^2 G(eta theta_max) with
  G(S) = int_S^inf (1 - J0(s)) / s^3 ds (the a_i^2 term is negligible there, a_i << theta_max).
  Small eta: Omega -> (chi_c^2 eta^2 / 4) (-ln(a^2 eta^2 / 4) + 1 - 2 C_E), Moliere's form with
  chi_a = a; large eta: Omega -> chi_c^2 / a^2, the mean number of collisions.
* Projected density f(x) = (1/pi) int_0^inf exp(-Omega(eta)) cos(eta x) deta. In a window
  |x| <= W: P(W) = (1/pi) int e^-Omega 2 sin(eta W)/eta deta and the second moment
  M2(W) = (1/pi) int e^-Omega g(eta) deta with
  g = 2 (W^2 sin(eta W)/eta + 2 W cos(eta W)/eta^2 - 2 sin(eta W)/eta^3) (series for small eta W).
  The windowed rms is sqrt(M2 / P). W = 2.5 theta_H.
* theta_H is the Moliere-Hanson width with Moliere's Z^2 and Fano's Delta b
  (``moliere_hanson.py`` conventions, constant energy in the 1 mm slab): b = ln chi_c^2 - ln
  chi_a^2 - (2 C_E - 1) + Delta b_Fano, B - ln B = b, theta_H = chi_c sqrt((B - 1.2) / 2). It
  stands in for ionmc's Gaussian; its windowed rms is that of a Gaussian of sigma theta_H.

Model terms (Geant4 v11.4.2 sources, pinned in P10 sections 7.1 and 7.4.1; role lineage)
-----------------------------------------------------------------------------------------
Tags and sha256 prefixes: [CS] ``G4CoulombScattering.cc`` db354d80, [ICS]
``G4IonCoulombScatteringModel.cc`` ef9b31d4, [IXS] ``G4IonCoulombCrossSection.cc`` 87989b93,
[WOKVI] ``G4WentzelOKandVIxSection.cc`` a64c94b1, [EMP] ``G4EmParameters.cc`` aaa143c1.

* Model per particle: [CS] L120-148, alpha and GenericIon -> G4IonCoulombScatteringModel, the
  proton -> G4eCoulombScatteringModel + G4WentzelOKandVIxSection.
* Ion path, nuclear only ([ICS] L125-148; [IXS] L186-219, no electron term); Rutherford factor
  Z^2 z^2 per atom ([IXS] L80-81, L190); universal screening length
  x = a0 (z^0.23 + Z^0.23), a0 = m_e c^2 / 0.88534, ScreenRSquare = alpha^2 x^2 ([IXS] L68,
  L167-182); Born factor 1.13 + corr 3.76 (z Z alpha / beta)^2 with corr = (10 pi Z z alpha)^0.04,
  also for hydrogen, no cap ([IXS] L153-155); no form factor, no spin term. Here: a_i^2 =
  (x/p)^2 (1.13 + corr 3.76 alpha'^2). The CM kinematics ([IXS] L107-135) leave beta unchanged
  and do not matter at these angles (P10 section 7.4.1).
* Proton path ([WOKVI] InitialiseA, SetupTarget, ComputeMaxElectronScattering): Thomas-Fermi
  a0 Z^(1/3) with the factor 1 + exp(-0.001 Z^2) for Z >= 2; screenZ = min(1.13 Z, 1.13 + 3.76
  alpha'^2) for Z > 1, no Born term for Z = 1; electron scattering with weight Z and the angle
  limit 1 - cos theta = T_max m_e / p^2. Used only for the superseded branch-A design.
* The mapping's physics ("Moliere-Fano"): one effective species with chi_c^2 (Moliere Z^2) and
  a_eff^2 = chi_a^2 exp(-Delta b_Fano), so that Omega reproduces b. Variants used for the bias
  decomposition: ``TF_nocorr`` (ion path with the Thomas-Fermi radius and corr = 1, no
  electrons), ``TF_corr`` (with corr), ``U_corr`` (the ion path itself).

Units: MeV, g/cm^2, rad; percentages as stated.
"""

from __future__ import annotations

import functools
import math
from collections.abc import Callable, Sequence

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
# (chi_c^2 [rad^2], screening angle a [rad], upper angle limit [rad] or None)
Species = tuple[float, float, float | None]
Case = tuple[int, float]  # (z, kinetic energy per nucleon [MeV/u])

ALPHA = 1.0 / 137.035999084  # CODATA 2018
ME_MEV = 0.51099895  # CODATA 2018 electron rest energy
M_P_MEV = 938.27208816  # CODATA 2018 proton rest energy
R_E_CM = 2.8179403262e-13  # CODATA 2018 classical electron radius
N_A = 6.02214076e23  # SI 2019, exact
EULER_GAMMA = 0.5772156649015329
K_CHI_C = 4.0 * math.pi * N_A * R_E_CM**2 * ME_MEV**2  # MeV^2 cm^2 / g per (Z^2/A)
# Nuclear rest energies [MeV] (rounded; the values of the scratch calculation).
M_NUC_MEV = {1: M_P_MEV, 2: 3727.379, 6: 11174.863, 8: 14895.080}
A_NUC = {1: 1, 2: 4, 6: 12, 8: 16}
_W_H = 2.0 * 1.00794 / (2.0 * 1.00794 + 15.9994)
WATER: tuple[tuple[int, float, float], ...] = (
    (1, 1.00794, _W_H),
    (8, 15.9994, 1.0 - _W_H),
)  # (Z, A [g/mol], mass fraction)
T_SLAB_G_CM2 = 0.1  # 1 mm water
MOLIERE_TF = 0.885  # Moliere chi_0 constant (moliere_hanson.py)
A0_G4 = ALPHA * ME_MEV / 0.88534  # [IXS] L68 / [WOKVI] InitialiseA, times alpha
FANO_U_IN_NEG = {1: 3.6}  # -u_in (Fano 1954; 5.0 otherwise, Gottschalk)
FANO_U_IN_NEG_DEFAULT = 5.0
HANSON_SHIFT = 1.2

WINDOW_K = 2.5  # window half-width in units of theta_H (plan row V5-SS)
ETA_MAX_UNITS = 40.0  # eta grid upper end in units of 1/theta_H
N_ETA = 200_001

DESIGN_CASES: tuple[Case, ...] = (
    (2, 150.0),
    (6, 150.0),
    (8, 150.0),
    (2, 250.0),
    (6, 250.0),
    (8, 250.0),
    (6, 430.0),
    (8, 430.0),
)
GATED_PAIRS: tuple[tuple[Case, Case], ...] = (
    ((2, 150.0), (6, 150.0)),
    ((2, 150.0), (8, 150.0)),
    ((2, 250.0), (6, 250.0)),
    ((2, 250.0), (8, 250.0)),
    ((6, 430.0), (8, 430.0)),
)
INFO_PAIRS: tuple[tuple[Case, Case], ...] = (
    ((6, 290.0), (8, 290.0)),
    ((6, 350.0), (8, 350.0)),
)
# Superseded r6 ion/proton design (branch A): the ions; protons at equal velocity.
BRANCH_A_IONS: tuple[Case, ...] = (
    (2, 150.0),
    (6, 290.0),
    (8, 350.0),
    (2, 250.0),
    (6, 430.0),
    (8, 430.0),
)
# All ion cases evaluated by the scratch calculation (for the range statements).
ALL_ION_CASES: tuple[Case, ...] = (
    (2, 150.0),
    (6, 290.0),
    (8, 350.0),
    (2, 250.0),
    (6, 430.0),
    (8, 430.0),
    (6, 150.0),
    (8, 290.0),
    (6, 350.0),
    (6, 250.0),
    (8, 150.0),
    (8, 250.0),
    (2, 290.0),
)

# Tolerance budget of V5-SS (P10 section 7.5.1; percent)
TOLERANCE_PCT = 3.0
SIGMA_MAX_PCT = 1.0
BIAS_BOUND_PCT = 0.9  # computed |bias| <= 0.69 plus <= 0.2 pp for the approximations
MODEL_BOUND_PCT = 2.0  # I11-rel gate
STAT_SIGMA_PCT = (0.7, 1.0)
# Scorer-free 1 mm G4_WATER g4em-standard_SS smoke, GenericIon(6,12) at 290 MeV/u, 100 histories,
# Ts/MaxStepNumber 1e8, no killed track: wall time per primary as relayed by the orchestrator
# (REF-35df82fc..., case carbon-water-290mevu-ss-slab1mm-smoke at ce4d3d6a). Execution only.
MEASURED_C290_S_PER_PRIMARY = 0.023


# --- kinematics and Moliere-Hanson width --------------------------------------------------------


def kinematics(z: int, t_per_u: float) -> tuple[float, float, float, float]:
    """Return (M, p, beta^2, p v) [MeV] for a nucleus of charge z (1, 2, 6, 8) at T/A."""
    m = M_NUC_MEV[z]
    t = t_per_u * A_NUC[z]
    p = math.sqrt(t * (t + 2.0 * m))
    b2 = p * p / (t + m) ** 2
    return m, p, b2, p * math.sqrt(b2)


def proton_energy_equal_velocity(z: int, t_per_u: float) -> float:
    """Proton kinetic energy [MeV] at the velocity of ion z at T/A: (gamma - 1) m_p."""
    gamma = 1.0 + t_per_u * A_NUC[z] / M_NUC_MEV[z]
    return (gamma - 1.0) * M_P_MEV


def solve_big_b(b: float) -> float:
    """Root B of B - ln B = b (Newton from b + ln b)."""
    x = b + math.log(b)
    for _ in range(80):
        x -= (x - math.log(x) - b) / (1.0 - 1.0 / x)
    return x


def moliere_fano(z: int, p: float, b2: float, pv: float, fano: bool = True) -> dict[str, float]:
    """Moliere collision number and Hanson width in the 1 mm water slab at constant energy.

    chi_c^2 = K_C z^2 sum_i w_i Z_i^2/A_i t/(pv)^2 (Moliere Z^2 form); ln chi_a^2 =
    sum_i W_i ln[chi_0,i^2 (1.13 + 3.76 (z Z_i alpha)^2/beta^2)], chi_0,i = alpha m_e Z_i^(1/3) /
    (0.885 p), W_i proportional to w_i Z_i^2/A_i; Delta b_Fano = sum_i w_i Z_i/A_i D_i /
    sum_i w_i Z_i^2/A_i, D_i = ln(1130 Z_i^-4/3 beta^2/(1-beta^2)) - u_in - beta^2/2 (Fano,
    Phys. Rev. 93 (1954) 117); b = ln chi_c^2 - ln chi_a^2 - (2 C_E - 1) + Delta b;
    theta_H = chi_c sqrt((B - 1.2)/2). Returns chi_c2, lna (ln chi_a^2), db, b, B, th.
    """
    zs = sum(w * zt * zt / a for zt, a, w in WATER)
    chi_c2 = K_CHI_C * z * z * zs * T_SLAB_G_CM2 / pv**2
    lna = 0.0
    for zt, a, w in WATER:
        wt = w * zt * zt / a / zs
        lna += wt * (
            2.0 * math.log(ALPHA * ME_MEV / MOLIERE_TF * zt ** (1.0 / 3.0) / p)
            + math.log(1.13 + 3.76 * (z * zt * ALPHA) ** 2 / b2)
        )
    db = 0.0
    if fano:
        num = 0.0
        for zt, a, w in WATER:
            neg_u = FANO_U_IN_NEG.get(zt, FANO_U_IN_NEG_DEFAULT)
            d = math.log(1130.0 * zt ** (-4.0 / 3.0) * b2 / (1.0 - b2)) + neg_u - 0.5 * b2
            num += w * zt / a * d
        db = num / zs
    b = math.log(chi_c2) - lna - (2.0 * EULER_GAMMA - 1.0) + db
    big_b = solve_big_b(b)
    return {
        "chi_c2": chi_c2,
        "lna": lna,
        "db": db,
        "b": b,
        "B": big_b,
        "th": math.sqrt(chi_c2 * (big_b - HANSON_SHIFT) / 2.0),
    }


# --- special functions (numpy only) -------------------------------------------------------------


@functools.cache
def _digamma_int() -> tuple[float, ...]:
    """psi(n) for n = 0 (unused placeholder 0.0), 1, ..., 89."""
    psi = [0.0, -EULER_GAMMA]
    for n in range(2, 90):
        psi.append(psi[-1] + 1.0 / (n - 1))
    return tuple(psi)


def one_minus_xk1(x: FloatArray | Sequence[float]) -> FloatArray:
    """1 - x K1(x) for x > 0.

    x < 2: the ascending series K1(x) = (1/x) + ln(x/2) I1(x) - (x/4) sum_k (psi(k+1) +
    psi(k+2)) (x^2/4)^k / (k! (k+1)!) (Abramowitz-Stegun 9.6.11), i.e.
    1 - x K1 = -x ln(x/2) I1(x) + (x^2/4) sum_k (psi(k+1)+psi(k+2)) (x^2/4)^k/(k!(k+1)!);
    x >= 2: K1(x) = int_0^inf exp(-x cosh t) cosh t dt (A-S 9.6.24), trapezoid on [0, 12].
    """
    psi = _digamma_int()
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    s = x < 2.0
    xs = x[s]
    q = xs * xs / 4.0
    i1 = np.zeros_like(xs)
    sm = np.zeros_like(xs)
    term = np.ones_like(xs)
    for k in range(40):
        if k > 0:
            term = term * q / (k * (k + 1))
        i1 += term
        sm += (psi[k + 1] + psi[k + 2]) * term
    i1 *= xs / 2.0
    out[s] = -xs * np.log(xs / 2.0) * i1 + q * sm
    xl = x[~s]
    if xl.size:
        tt = np.linspace(0.0, 12.0, 4001)
        k1 = np.trapezoid(np.exp(-np.outer(xl, np.cosh(tt))) * np.cosh(tt), tt, axis=1)
        out[~s] = 1.0 - xl * k1
    return out


def bessel_j0(s: FloatArray | Sequence[float], chunk: int = 1000) -> FloatArray:
    """J0(s) = (1/pi) int_0^pi cos(s sin tau) dtau (A-S 9.1.18), trapezoid with 2049 nodes."""
    s = np.asarray(s, dtype=float)
    tau = np.linspace(0.0, math.pi, 2049)
    st = np.sin(tau)
    out = np.empty_like(s)
    for i in range(0, s.size, chunk):
        out[i : i + chunk] = np.trapezoid(np.cos(np.outer(s[i : i + chunk], st)), tau, axis=1)
    return out / math.pi


@functools.cache
def _g_table() -> tuple[FloatArray, FloatArray]:
    """G(S) = int_S^inf (1 - J0(s))/s^3 ds on a log grid 1e-5..1e3 (tail 1/(2 S^2) beyond)."""
    sg = np.logspace(-5.0, 3.0, 6001)
    f = (1.0 - bessel_j0(sg)) / sg**3
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(sg))))
    gg = (cum[-1] - cum) + 1.0 / (2.0 * sg[-1] ** 2)
    return sg, gg


def g_tail(s: FloatArray) -> FloatArray:
    """G(S) by log interpolation; G(S) ~ G(S0) + (1/4) ln(S0/S) below the grid."""
    sg, gg = _g_table()
    s = np.asarray(s, dtype=float)
    out = np.interp(np.log(s), np.log(sg), gg)
    small = s < sg[0]
    out[small] = gg[0] + 0.25 * np.log(sg[0] / s[small])
    return out


def omega(eta: FloatArray, species: Sequence[Species]) -> FloatArray:
    """Omega(eta) = sum_i chi_c,i^2 (1 - eta a_i K1(eta a_i))/a_i^2 - [eta^2 G(eta theta_max)]."""
    tot = np.zeros_like(eta)
    for c2, a, thm in species:
        full = one_minus_xk1(eta * a) / (2.0 * a * a)
        if thm is not None:
            full = full - eta**2 * g_tail(eta * thm)
        tot += 2.0 * c2 * full
    return tot


# --- windowed rms -------------------------------------------------------------------------------


def windowed_rms_char(
    char: Callable[[FloatArray], FloatArray],
    window: float,
    th_scale: float,
    eta_max_units: float = ETA_MAX_UNITS,
    n: int = N_ETA,
) -> tuple[float, float]:
    """(rms, P) of the projected density with characteristic function char(eta), |x| <= window."""
    eta = np.linspace(1e-9, eta_max_units / th_scale, n)
    e = char(eta)
    s, c = np.sin(eta * window), np.cos(eta * window)
    p0 = float(np.trapezoid(e * 2.0 * s / eta, eta)) / math.pi
    u = eta * window
    w = window
    g = np.where(
        u > 0.05,
        2.0 * (w * w * s / eta + 2.0 * w * c / eta**2 - 2.0 * s / eta**3),
        2.0 * w**3 * (1.0 / 3.0 - u * u / 10.0 + u**4 / 168.0 - u**6 / 6480.0),
    )
    m2 = float(np.trapezoid(e * g, eta)) / math.pi
    return math.sqrt(m2 / p0), p0


def windowed_rms(
    species: Sequence[Species],
    window: float,
    th_scale: float,
    eta_max_units: float = ETA_MAX_UNITS,
    n: int = N_ETA,
) -> tuple[float, float]:
    """(rms, P) in |x| <= window for exp(-Omega(eta)) of the given species."""
    return windowed_rms_char(
        lambda eta: np.exp(-omega(eta, species)), window, th_scale, eta_max_units, n
    )


def gauss_windowed(sigma: float, window: float) -> float:
    """rms of a Gaussian of width sigma truncated to |x| <= window (closed form)."""
    k = window / sigma
    phi = math.exp(-k * k / 2.0) / math.sqrt(2.0 * math.pi)
    prob = math.erf(k / math.sqrt(2.0))
    return sigma * math.sqrt(1.0 - 2.0 * k * phi / prob)


# --- species --------------------------------------------------------------------------------------


def species_ion_path(
    z: int,
    p: float,
    b2: float,
    radius: str = "universal",
    corr_on: bool = True,
    electrons: bool = False,
    m: float | None = None,
) -> list[Species]:
    """Geant4 ion path ([ICS], [IXS]) per element of water; variants for the decomposition.

    radius "universal": x = a0 (z^0.23 + Z^0.23) ([IXS] L167-182); "tf": a0 Z^(1/3).
    corr_on: corr = (10 pi Z z alpha)^0.04 ([IXS] L153-155), else 1. electrons=True adds the
    hypothetical electron species (weight Z, limit 2 m_e / M) that the ion path does NOT have.
    """
    sp: list[Species] = []
    for zt, a, w in WATER:
        c2 = K_CHI_C * z * z * (w * zt * zt / a) * T_SLAB_G_CM2 / (p * p * b2)
        x = A0_G4 * (z**0.23 + zt**0.23) if radius == "universal" else A0_G4 * zt ** (1.0 / 3.0)
        corr = (5.0 * 2.0 * math.pi * zt * z * ALPHA) ** 0.04 if corr_on else 1.0
        a2 = (x / p) ** 2 * (1.13 + corr * 3.76 * (zt * z * ALPHA) ** 2 / b2)
        sp.append((c2, math.sqrt(a2), None))
        if electrons:
            if m is None:
                raise ValueError("electrons=True needs the projectile mass m")
            ce2 = K_CHI_C * z * z * (w * zt / a) * T_SLAB_G_CM2 / (p * p * b2)
            sp.append((ce2, math.sqrt(a2), 2.0 * ME_MEV / m))
    return sp


def species_proton_path(
    p: float, b2: float, m: float = M_P_MEV, electrons: bool = True
) -> tuple[list[Species], float]:
    """Geant4 proton path ([WOKVI]) per element of water and the electron angle limit.

    Screening: a0^2 Z^(2/3) (1 + exp(-0.001 Z^2) for Z >= 2) times min(1.13 Z, 1.13 + 3.76
    (Z alpha)^2/beta^2) for Z > 1 and 1 for Z = 1; electrons with weight Z up to
    theta_e = sqrt(2 T_max m_e)/p (1 - cos theta = T_max m_e / p^2).
    """
    sp: list[Species] = []
    g = 1.0 / math.sqrt(1.0 - b2)
    tmax = 2.0 * ME_MEV * b2 * g * g / (1.0 + 2.0 * g * ME_MEV / m + (ME_MEV / m) ** 2)
    the = math.sqrt(2.0 * tmax * ME_MEV) / p
    for zt, a, w in WATER:
        c2 = K_CHI_C * (w * zt * zt / a) * T_SLAB_G_CM2 / (p * p * b2)
        x2 = A0_G4**2 * zt ** (2.0 / 3.0) * (1.0 if zt == 1 else 1.0 + math.exp(-0.001 * zt * zt))
        if zt == 1:
            a2 = x2 / p**2
        else:
            a2 = min(1.13 * zt, 1.13 + 3.76 * (zt * ALPHA) ** 2 / b2) * x2 / p**2
        sp.append((c2, math.sqrt(a2), None))
        if electrons:
            ce2 = K_CHI_C * (w * zt / a) * T_SLAB_G_CM2 / (p * p * b2)
            sp.append((ce2, math.sqrt(a2), the))
    return sp, the


def species_moliere(
    z: int, p: float, b2: float, pv: float, fano: bool = True
) -> tuple[list[Species], dict[str, float]]:
    """One effective species with chi_c^2 and a_eff^2 = chi_a^2 exp(-Delta b), i.e. Moliere b."""
    mf = moliere_fano(z, p, b2, pv, fano=fano)
    a_eff = math.sqrt(math.exp(mf["lna"] - mf["db"]))
    return [(mf["chi_c2"], a_eff, None)], mf


# --- per-case widths and pair decomposition -------------------------------------------------------


def case_widths(
    z: int,
    t_per_u: float,
    k_window: float = WINDOW_K,
    eta_max_units: float = ETA_MAX_UNITS,
    n: int = N_ETA,
    variants: bool = True,
) -> dict[str, float]:
    """Windowed rms widths [rad] of one case (window k_window theta_H).

    Keys: b, B, th (theta_H), db, gauss (the Gaussian of theta_H), U_corr (the Geant4 ion path);
    with variants: MF (Moliere-Fano), TF_nocorr, TF_corr. For z = 1 also P_ref (proton path)
    and P_ref_noe (without its electrons), and the_p.
    """
    m, p, b2, pv = kinematics(z, t_per_u)
    mf = moliere_fano(z, p, b2, pv)
    w = k_window * mf["th"]
    th = mf["th"]
    out = {"z": float(z), "tu": t_per_u, "b": mf["b"], "B": mf["B"], "th": th, "db": mf["db"]}
    out["gauss"] = gauss_windowed(th, w)

    def wr(sp: Sequence[Species]) -> float:
        return windowed_rms(sp, w, th, eta_max_units, n)[0]

    out["U_corr"] = wr(species_ion_path(z, p, b2, "universal", True))
    if variants:
        out["MF"] = wr(species_moliere(z, p, b2, pv, fano=True)[0])
        out["TF_nocorr"] = wr(species_ion_path(z, p, b2, "tf", False))
        out["TF_corr"] = wr(species_ion_path(z, p, b2, "tf", True))
    if z == 1:
        sp, the = species_proton_path(p, b2)
        out["P_ref"] = wr(sp)
        out["P_ref_noe"] = wr(species_proton_path(p, b2, electrons=False)[0])
        out["the_p"] = the
    return out


def electron_variant_width(z: int, t_per_u: float) -> float:
    """Windowed rms of the ion path WITH hypothetical electron scattering (diagnostic only)."""
    m, p, b2, pv = kinematics(z, t_per_u)
    mf = moliere_fano(z, p, b2, pv)
    sp = species_ion_path(z, p, b2, "universal", True, electrons=True, m=m)
    return windowed_rms(sp, WINDOW_K * mf["th"], mf["th"])[0]


def pair_terms(ra: dict[str, float], rb: dict[str, float]) -> dict[str, float]:
    """Branch-B pair A/B at a common velocity: signal, bias and its telescoping decomposition.

    model = gauss_A/gauss_B (ionmc stand-in), ref = U_corr_A/U_corr_B (Geant4 ion path).
    bias = model/ref - 1 = (1+estimator)(1+electrons)(1+corr)(1+radius) - 1 with
    estimator = model/(MF ratio), electrons = (MF ratio)/(TF_nocorr ratio),
    corr = (TF_nocorr ratio)/(TF_corr ratio), radius = (TF_corr ratio)/(U_corr ratio).
    All values are fractions (not percent).
    """

    def r(key: str) -> float:
        return ra[key] / rb[key]

    model, ref = r("gauss"), r("U_corr")
    return {
        "signal_ref": ref - 1.0,
        "signal_model": ra["th"] / rb["th"] - 1.0,
        "bias": model / ref - 1.0,
        "estimator": model / r("MF") - 1.0,
        "electrons": r("MF") / r("TF_nocorr") - 1.0,
        "corr": r("TF_nocorr") / r("TF_corr") - 1.0,
        "radius": r("TF_corr") / r("U_corr") - 1.0,
    }


def branch_b_table(
    pairs: Sequence[tuple[Case, Case]] = GATED_PAIRS,
) -> list[tuple[Case, Case, dict[str, float]]]:
    """pair_terms for each pair (widths computed once per case)."""
    cache: dict[Case, dict[str, float]] = {}
    rows = []
    for a, b in pairs:
        for c in (a, b):
            if c not in cache:
                cache[c] = case_widths(*c)
        rows.append((a, b, pair_terms(cache[a], cache[b])))
    return rows


def branch_a_table() -> list[tuple[Case, dict[str, float]]]:
    """Superseded ion/proton design: ratio_MHgauss/ratio_ref - 1 at equal velocity."""
    rows = []
    for z, tu in BRANCH_A_IONS:
        ri = case_widths(z, tu)
        tp = proton_energy_equal_velocity(z, tu)
        rp = case_widths(1, tp)
        model = ri["gauss"] / rp["gauss"]
        ref = ri["U_corr"] / rp["P_ref"]
        rows.append(
            (
                (z, tu),
                {
                    "total": model / ref - 1.0,
                    "ion_ref_vs_mf": ri["U_corr"] / ri["MF"] - 1.0,
                    "p_ref_vs_mf": rp["P_ref"] / rp["MF"] - 1.0,
                    "p_ref_noe_vs_mf": rp["P_ref_noe"] / rp["MF"] - 1.0,
                    "t_p": tp,
                },
            )
        )
    return rows


# --- discrimination (Hanson-ratio alternatives) ---------------------------------------------------


def hanson_width_variant(
    z: int, t_per_u: float, born: bool = True, z_b: int | None = None
) -> float:
    """theta_H with the Born term removed (born=False) or with b evaluated for charge z_b.

    chi_c keeps the charge z; only the collision number b uses z_b (z_b = 1: a b without
    z-dependence); the Fano Delta b is that of the full model.
    """
    m, p, b2, pv = kinematics(z, t_per_u)
    zz = z if z_b is None else z_b
    zs = sum(w * zt * zt / a for zt, a, w in WATER)
    chi_c2 = K_CHI_C * z * z * zs * T_SLAB_G_CM2 / pv**2
    chi_c2_b = K_CHI_C * zz * zz * zs * T_SLAB_G_CM2 / pv**2
    lna = 0.0
    for zt, a, w in WATER:
        wt = w * zt * zt / a / zs
        f = 1.13 + (3.76 * (zz * zt * ALPHA) ** 2 / b2 if born else 0.0)
        lna += wt * (
            2.0 * math.log(ALPHA * ME_MEV / MOLIERE_TF * zt ** (1.0 / 3.0) / p) + math.log(f)
        )
    db = moliere_fano(z, p, b2, pv)["db"]
    b = math.log(chi_c2_b) - lna - (2.0 * EULER_GAMMA - 1.0) + db
    big_b = solve_big_b(b)
    return math.sqrt(chi_c2 * (big_b - HANSON_SHIFT) / 2.0)


def discrimination(a: Case, b: Case) -> dict[str, float]:
    """Hanson ratio A/B - 1 and the changes of the ratio for "no Born" and "b without z"."""
    full = hanson_width_variant(*a) / hanson_width_variant(*b)
    no_born = hanson_width_variant(*a, born=False) / hanson_width_variant(*b, born=False)
    no_z = hanson_width_variant(*a, z_b=1) / hanson_width_variant(*b, z_b=1)
    return {"full": full - 1.0, "no_born": no_born / full - 1.0, "b_without_z": no_z / full - 1.0}


# --- cost -----------------------------------------------------------------------------------------


def events_per_slab(z: int, t_per_u: float) -> float:
    """Mean number of elastic events per primary in the 1 mm slab, ion path: sum chi_c^2/a^2."""
    m, p, b2, pv = kinematics(z, t_per_u)
    return sum(c2 / a**2 for c2, a, _ in species_ion_path(z, p, b2))


def moliere_omega_nofano(z: int, t_per_u: float) -> float:
    """Moliere's Omega_0 = exp(b) without the Fano term (for the 0.43-0.59 comparison)."""
    m, p, b2, pv = kinematics(z, t_per_u)
    return math.exp(moliere_fano(z, p, b2, pv, fano=False)["b"])


def cost_table(
    us_per_event: float = 2.0, primaries: float = 2e4, cases: Sequence[Case] = DESIGN_CASES
) -> tuple[list[tuple[Case, float, float, float]], float]:
    """Rows (case, events/primary, Omega_0, seconds) and the total [s] at us_per_event."""
    rows = []
    tot = 0.0
    for c in cases:
        n_ev = events_per_slab(*c)
        sec = n_ev * primaries * us_per_event * 1e-6
        tot += sec
        rows.append((c, n_ev, moliere_omega_nofano(*c), sec))
    return rows, tot


def cost_from_measured(
    s_per_primary: float,
    measured_case: Case = (6, 290.0),
    primaries: float = 2e4,
    cases: Sequence[Case] = DESIGN_CASES,
) -> tuple[list[tuple[Case, float]], float, float]:
    """Scale a measured wall time per primary of one case by the event counts of the others.

    Returns (rows (case, seconds), total seconds, implied microseconds per event). The scaling
    assumes the cost is proportional to the number of elastic events (an upper bound if part
    of the measured time is fixed per-run overhead).
    """
    n_ref = events_per_slab(*measured_case)
    rows = [(c, s_per_primary * events_per_slab(*c) / n_ref * primaries) for c in cases]
    return rows, sum(s for _, s in rows), s_per_primary / n_ref * 1e6


# --- checks (numerical convergence, window anchor, probability) -----------------------------------


def pair_bias(
    a: Case,
    b: Case,
    k_window: float = WINDOW_K,
    eta_max_units: float = ETA_MAX_UNITS,
    n: int = N_ETA,
) -> float:
    """Pair bias model/ref - 1 [fraction] with a changed window or eta grid."""
    ra = case_widths(*a, k_window=k_window, eta_max_units=eta_max_units, n=n, variants=False)
    rb = case_widths(*b, k_window=k_window, eta_max_units=eta_max_units, n=n, variants=False)
    return (ra["gauss"] / rb["gauss"]) / (ra["U_corr"] / rb["U_corr"]) - 1.0


def window_probability(z: int, t_per_u: float, k_window: float) -> float:
    """P(|x| < k theta_H) of the Moliere-Fano projected distribution."""
    m, p, b2, pv = kinematics(z, t_per_u)
    sp, mf = species_moliere(z, p, b2, pv)
    return windowed_rms(sp, k_window * mf["th"], mf["th"])[1]


# --- budget and gate ------------------------------------------------------------------------------


def budget(
    bias: float = BIAS_BOUND_PCT,
    model: float = MODEL_BOUND_PCT,
    stat: Sequence[float] = STAT_SIGMA_PCT,
) -> dict[str, float]:
    """Tolerance budget [percent]: quadrature (at the largest sigma) and conservative range.

    quadrature = sqrt(bias^2 + model^2 + sigma_max^2); conservative = bias + sqrt(model^2 +
    sigma^2) for sigma in stat (bias added linearly).
    """
    smax = max(stat)
    cons = [bias + math.hypot(model, s) for s in stat]
    return {
        "quadrature": math.sqrt(bias**2 + model**2 + smax**2),
        "conservative_min": min(cons),
        "conservative_max": max(cons),
    }


def gate(
    deviation_pct: float,
    sigma_pct: float,
    tol: float = TOLERANCE_PCT,
    sigma_max: float = SIGMA_MAX_PCT,
) -> str:
    """V5-SS rule per gated pair: INCONCLUSIVE (counts as FAIL) if sigma > 1 %, else PASS iff
    |ratio_model/ratio_ref - 1| <= 3 %."""
    if sigma_pct > sigma_max:
        return "INCONCLUSIVE"
    return "PASS" if abs(deviation_pct) <= tol else "FAIL"


# --- report ---------------------------------------------------------------------------------------


def _rng(vals: Sequence[float]) -> str:
    return f"{100 * min(vals):+.2f} to {100 * max(vals):+.2f} %"


def main() -> None:
    """Print the aggregates of P10 section 7.5.1 (SEEN pre-freeze planning values)."""
    print("V5-SS branch-B budget (validation/scripts/ion_em/ss_budget.py); SEEN pre-freeze")
    allc = {c: case_widths(*c) for c in ALL_ION_CASES}
    print("\n(a) absolute windowed widths, 13 ion cases")
    ele = [r["TF_nocorr"] / r["MF"] - 1 for r in allc.values()]
    ele8 = [allc[c]["TF_nocorr"] / allc[c]["MF"] - 1 for c in DESIGN_CASES]
    rad = [r["U_corr"] / r["TF_corr"] - 1 for r in allc.values()]
    cor = [r["TF_corr"] / r["TF_nocorr"] - 1 for r in allc.values()]
    mfg = [r["MF"] / r["gauss"] - 1 for r in allc.values()]
    print(f" no electrons (TF_nocorr/MF - 1): {_rng(ele)}; the 8 design cases {_rng(ele8)}")
    print(f" corr^0.04 (TF_corr/TF_nocorr - 1): {_rng(cor)}")
    print(f" universal radius (U_corr/TF_corr - 1): {_rng(rad)}")
    print(f" Moliere_w/Gaussian_w - 1: {_rng(mfg)}")
    fac = [allc[c]["gauss"] / allc[c]["U_corr"] - 1 for c in DESIGN_CASES]
    print(f" expected reference factor theta_Gauss,w/theta_SS,w - 1 (8 cases): {_rng(fac)}")

    print("\n(b) superseded branch A (ion/proton at equal velocity): ratio_model/ratio_ref - 1")
    for (z, tu), r in branch_a_table():
        print(
            f" z={z} {tu:g} MeV/u (T_p {r['t_p']:.2f} MeV): {100 * r['total']:+.2f} % "
            f"[ion path vs MF {100 * r['ion_ref_vs_mf']:+.2f} %, proton path vs MF "
            f"{100 * r['p_ref_vs_mf']:+.2f} %, without its electrons "
            f"{100 * r['p_ref_noe_vs_mf']:+.2f} %]"
        )

    print("\n(c) branch B, gated pairs: signal, bias and decomposition [%]")
    rows = branch_b_table(GATED_PAIRS + INFO_PAIRS)
    for a, b, t in rows:
        tag = "gated" if (a, b) in GATED_PAIRS else "info"
        print(
            f" {a[0]}/{b[0]} {a[1]:g} ({tag}): signal {100 * t['signal_ref']:+.2f}, "
            f"bias {100 * t['bias']:+.2f} [radius {100 * t['radius']:+.2f} / corr "
            f"{100 * t['corr']:+.2f} / electrons {100 * t['electrons']:+.2f} / estimator "
            f"{100 * t['estimator']:+.2f}]"
        )
    gated = [t for a, b, t in rows if (a, b) in GATED_PAIRS]
    print(f" gated bias range {_rng([t['bias'] for t in gated])}")
    print(f" He-pair signal range {_rng([t['signal_ref'] for t in gated[:4]])}")
    for a, b in GATED_PAIRS[:2]:
        print(f" G4-electrons-if-added change {a[0]}/{b[0]} {a[1]:g}: ", end="")
        ue = electron_variant_width(*a) / electron_variant_width(*b)
        u = allc[a]["U_corr"] / allc[b]["U_corr"]
        print(f"{100 * (u / ue - 1):+.2f} %")

    print("\n(d) discrimination (Hanson ratio changes) [%]")
    for a, b in GATED_PAIRS + INFO_PAIRS:
        d = discrimination(a, b)
        print(
            f" {a[0]}/{b[0]} {a[1]:g}: Hanson ratio-1 {100 * d['full']:+.2f}, no-Born "
            f"{100 * d['no_born']:+.2f}, b without z {100 * d['b_without_z']:+.2f}"
        )

    print("\n(e) checks")
    for a, b in (GATED_PAIRS[0], GATED_PAIRS[4]):
        base = pair_bias(a, b)
        fine = pair_bias(a, b, eta_max_units=60.0, n=600_001)
        lo, hi = pair_bias(a, b, k_window=2.375), pair_bias(a, b, k_window=2.625)
        print(
            f" {a[0]}/{b[0]} {a[1]:g}: base {100 * base:+.3f}, fine grid {100 * fine:+.3f}, "
            f"window x0.95 {100 * lo:+.3f}, x1.05 {100 * hi:+.3f} [%]"
        )
    print(f" P(|x| < 50 theta_H), C 290: {window_probability(6, 290.0, 50.0):.5f}")

    print("\n(f) cost (8 design cases x 2e4 primaries)")
    ctab, tot = cost_table(2.0)
    for c, n_ev, om, sec in ctab:
        print(
            f" z={c[0]} {c[1]:g}: events/primary {n_ev:.3g}, Omega_0 {om:.3g}, "
            f"ratio {n_ev / om:.2f}, {sec:.0f} s at 2 us/event"
        )
    print(f" total {tot:.0f} s = {tot / 3600:.2f} h at 2 us/event; {3 * tot / 3600:.2f} h at 6 us")
    print(f" events/primary C 290: {events_per_slab(6, 290.0):.4g}")
    mrows, mtot, us_ev = cost_from_measured(MEASURED_C290_S_PER_PRIMARY)
    print(
        f" scaled from the 1 mm C 290 smoke ({1e3 * MEASURED_C290_S_PER_PRIMARY:g} ms/primary, "
        f"{us_ev:.2f} us/event): per case {min(s for _, s in mrows):.0f}-"
        f"{max(s for _, s in mrows):.0f} s, total {mtot:.0f} s = {mtot / 3600:.2f} h"
    )

    bud = budget()
    print(
        f"\n(g) budget: quadrature {bud['quadrature']:.2f} %, conservative "
        f"{bud['conservative_min']:.2f}-{bud['conservative_max']:.2f} %; tolerance "
        f"{TOLERANCE_PCT} %, sigma <= {SIGMA_MAX_PCT} %"
    )


if __name__ == "__main__":
    main()
