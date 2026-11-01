"""ENDF-6 MF6 LAW=5 (charged-particle elastic scattering) parser and reconstruction (slice C).

Acceptance Amendment 14 (b), implemented literally for H-1: "<= 150 MeV: the LA150 Hale LAW=5
LTP=1 (LIDP=1) reconstruction (ENDF-102 eqs 6.9/6.10/6.14, b/sr CM, factor 2 pi on integration).
The tables interpolate the ratio sigma_e/sigma_c and never P_NI. The transported density is
NI = sigma_e - sigma_c over theta_CM >= 16.26 deg (|mu_CM| <= 0.96, the NJOY `umin` convention)."

Conventions (ENDF-102 section 6.2.7, NJOY2016 ``acecpe``/``coul``; all b/sr in the CM frame)
-------------------------------------------------------------------------------------------
* Records: ``[MAT,6,MT / SPI, 0, LIDP, 0, NR, NE / E_int]TAB2`` then per energy
  ``[MAT,6,MT / 0, E, LTP, 0, NW, NL / A]LIST``. LTP=1: ``LIDP=1`` NW = 3NL+3 (``b_0..b_NL``,
  then ``Re a_l, Im a_l``, l = 0..NL); ``LIDP=0`` NW = 4NL+3 (``b_0..b_2NL``, then the a_l).
  LTP=12 (also 14, 15 in the format; only 12 is accepted here): NW = 2NL, ``NL`` pairs
  ``(mu_i, P_NI(mu_i))`` lin-lin. Any other LTP, any LAW other than 5, a TAB2 interpolation other
  than lin-lin or an inconsistent NW raises :class:`UnsupportedEndfError` (fail closed).
* Coulomb (eqs 6.9/6.10, non-relativistic, no screening): ``sigma_cd = eta^2 / (k^2 (1-mu)^2)``;
  identical particles ``sigma_ci = 2 eta^2 / (k^2 (1-mu^2)) [(1+mu^2)/(1-mu^2) +
  (-1)^{2s} cos(eta ln((1+mu)/(1-mu))) / (2s+1)]`` (s = 1/2: coefficient -1/2); ``k`` and
  ``eta`` of eqs 6.11/6.12 (``k`` in b^{-1/2}).
* LTP=1, identical particles (eq 6.14): ``sigma_ei = sigma_ci - 2 eta/(1-mu^2) Re{ sum_l [(1+mu)
  e^{i eta ln((1-mu)/2)} + (-1)^l (1-mu) e^{i eta ln((1+mu)/2)}] (2l+1)/2 a_l P_l(mu) }
  + sum_{l=0}^{NL} (4l+1)/2 b_l P_{2l}(mu)`` (phase sign exp(+i eta ln ...) as NJOY ``coul``);
  distinguishable (eq 6.13) ``sigma_ed = sigma_cd - 2 eta/(1-mu) Re{e^{i eta ln((1-mu)/2)}
  sum_l (2l+1)/2 a_l P_l} + sum_{l=0}^{2NL} (2l+1)/2 b_l P_l``.
* LTP=12 (eq 6.19): ``sigma_NI(mu) = sigma_MF3/MT2(E) P_NI(mu)`` over ``[mu_min, mu_max]``;
  ``sigma_e = sigma_c + sigma_NI``.
* **2 pi and the half sphere.** All densities are per steradian, so an integrated cross section is
  ``2 pi int dmu``. For identical particles (H-1) each event is counted ONCE: the integral runs
  over the half sphere ``0 <= mu <= mu_cut`` (``theta_CM`` from the cut to 90 deg), which is the
  ENDF identical-particle range (``mu_min = 0``). The transported density is symmetric in mu.
* Energy interpolation: the ratio ``R(mu) = sigma_e / sigma_c`` of the two bracketing ENDF nodes
  (each at its own k, eta) is interpolated linearly in E (TAB2 INT = 2), and
  ``NI(mu, E) = sigma_c(mu, E) (R(mu, E) - 1)`` with the Coulomb density at E; P_NI is never
  interpolated.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ionmc.data.endf6 import (
    EndfMaterial,
    EndfSection,
    UnsupportedEndfError,
    _cont_from_line,
    _Reader,
)

AMU_EV = 931.49410242e6
MN_AMU = 1.00866491595
HBARC_EV_FM = 197.3269804e6
ALPHA = 1.0 / 137.035999084
MU_CUT_PP = 0.96
"""|mu_CM| bound of the transported p-p density: theta_CM >= 16.26 deg (NJOY ``umin``)."""


@dataclass(frozen=True)
class Law5Energy:
    """One LIST record: incident energy [eV], ``ltp``, ``nl`` and the coefficient array."""

    energy_ev: float
    ltp: int
    nl: int
    data: NDArray[np.float64]


@dataclass(frozen=True)
class Law5Section:
    """MF6 LAW=5 section of one MT: header values, ``spi`` (spin), ``lidp`` and the records."""

    za: float
    awr: float
    lct: int
    zap: int
    awp: float
    spi: float
    lidp: int
    records: tuple[Law5Energy, ...]

    @property
    def energies_ev(self) -> NDArray[np.float64]:
        return np.array([r.energy_ev for r in self.records])


def _expected_nw(ltp: int, nl: int, lidp: int) -> int:
    if ltp == 1:
        return 3 * nl + 3 if lidp == 1 else 4 * nl + 3
    if ltp == 12:
        return 2 * nl
    raise UnsupportedEndfError(f"MF6 LAW=5 LTP={ltp} is not supported (only 1 and 12)")


def parse_law5(section: EndfSection) -> Law5Section:
    """Parse an MF6 section whose single product has LAW=5 (module docstring); fail closed."""
    if section.mf != 6:
        raise UnsupportedEndfError("not an MF6 section")
    r = _Reader(section.lines)
    head = r.cont()
    if head.n1 != 1:
        raise UnsupportedEndfError(f"MF6 LAW=5 section with NK={head.n1} products")
    phead, _yield = r.tab1()
    if phead.l2 != 5:
        raise UnsupportedEndfError(f"MF6 product LAW={phead.l2}, expected 5")
    thead, _nbt, interp = r.tab2()
    lidp, ne = thead.l1, thead.n2
    if lidp not in (0, 1):
        raise UnsupportedEndfError(f"MF6 LAW=5 LIDP={lidp} is not supported")
    if np.any(interp != 2):
        raise UnsupportedEndfError("LAW=5 energy interpolation other than lin-lin (INT=2)")
    recs: list[Law5Energy] = []
    for _ in range(ne):
        lh = r.cont()
        body = r.list_body(lh)
        ltp, nw, nl = lh.l1, lh.n1, lh.n2
        if nw != _expected_nw(ltp, nl, lidp):
            raise UnsupportedEndfError(f"LAW=5 LTP={ltp} NW={nw} inconsistent with NL={nl}")
        recs.append(Law5Energy(lh.c2, ltp, nl, body))
    if not r.exhausted:
        raise UnsupportedEndfError("unexpected records after the LAW=5 distribution")
    energies = [x.energy_ev for x in recs]
    if energies != sorted(energies) or len(set(energies)) != len(energies):
        raise UnsupportedEndfError("LAW=5 incident energies must be strictly increasing")
    return Law5Section(head.c1, head.c2, head.l2, int(round(phead.c1)), phead.c2, thead.c1, lidp,
                       tuple(recs))  # fmt: skip


def projectile_awi(material: EndfMaterial) -> float:
    """AWI (projectile mass in neutron masses) of the MF1/MT451 header (third record, C1)."""
    return float(_cont_from_line(material.sections[(1, 451)].lines[2]).c1)


def coulomb_params(e_ev: float, awr: float, awi: float, z1: int, z2: int) -> tuple[float, float]:
    """``(k [b^{-1/2}], eta)`` of ENDF-102 eqs 6.11/6.12 (non-relativistic) for a projectile of
    ``awi`` neutron masses and charge ``z1`` on a target of ``awr`` neutron masses and ``z2``."""
    m1 = awi * MN_AMU
    a = awr / awi
    k = a / (1.0 + a) * math.sqrt(2.0 * m1 * AMU_EV * e_ev) / HBARC_EV_FM * 10.0
    eta = z1 * z2 * math.sqrt(ALPHA**2 * m1 * AMU_EV / (2.0 * e_ev))
    return k, eta


def sigma_c(
    mu: NDArray[np.float64], k: float, eta: float, lidp: int, spi: float
) -> NDArray[np.float64]:
    """Coulomb density [b/sr], eq 6.9 (``lidp`` 0) or 6.10 (``lidp`` 1)."""
    mu = np.asarray(mu, dtype=np.float64)
    if lidp == 0:
        return eta**2 / k**2 / (1.0 - mu) ** 2
    sign = (-1.0) ** round(2 * spi)
    return (2.0 * eta**2 / k**2 / (1.0 - mu**2)) * (
        (1.0 + mu**2) / (1.0 - mu**2)
        + sign * np.cos(eta * np.log((1.0 + mu) / (1.0 - mu))) / (2.0 * spi + 1.0)
    )


def legendre(x: NDArray[np.float64], n: int) -> list[NDArray[np.float64]]:
    """``[P_0(x), ..., P_n(x)]`` by the three-term recurrence."""
    p = [np.ones_like(x), x]
    for lval in range(1, n):
        p.append(((2 * lval + 1) * x * p[lval] - lval * p[lval - 1]) / (lval + 1))
    return p[: n + 1]


def sigma_e_ltp1(
    mu: NDArray[np.float64], rec: Law5Energy, lidp: int, k: float, eta: float, spi: float
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """``(sigma_e, nuclear term, interference term)`` [b/sr] of an LTP=1 record (eqs 6.13/6.14)."""
    if rec.ltp != 1:
        raise UnsupportedEndfError("sigma_e_ltp1 needs an LTP=1 record")
    mu = np.asarray(mu, dtype=np.float64)
    nl, a = rec.nl, rec.data
    nb = nl + 1 if lidp == 1 else 2 * nl + 1
    b, ca = a[:nb], a[nb : nb + 2 * (nl + 1)]
    al = ca[0::2] + 1j * ca[1::2]
    p = legendre(mu, 2 * nl + 1)
    e1 = np.exp(1j * eta * np.log((1.0 - mu) / 2.0))
    if lidp == 1:
        e2 = np.exp(1j * eta * np.log((1.0 + mu) / 2.0))
        nuc = sum((4 * lval + 1) / 2 * b[lval] * p[2 * lval] for lval in range(nl + 1))
        s = sum(
            ((1.0 + mu) * e1 + (-1.0) ** lval * (1.0 - mu) * e2)
            * (2 * lval + 1)
            / 2
            * al[lval]
            * p[lval]
            for lval in range(nl + 1)
        )
        intf = -2.0 * eta / (1.0 - mu**2) * np.real(s)
    else:
        nuc = sum((2 * lval + 1) / 2 * b[lval] * p[lval] for lval in range(2 * nl + 1))
        s = e1 * sum((2 * lval + 1) / 2 * al[lval] * p[lval] for lval in range(nl + 1))
        intf = -2.0 * eta / (1.0 - mu) * np.real(s)
    return sigma_c(mu, k, eta, lidp, spi) + intf + nuc, np.asarray(nuc), np.asarray(intf)


def ratio_ltp1(
    sec: Law5Section, e_ev: float, mu: NDArray[np.float64], awi: float, z1: int, z2: int
) -> NDArray[np.float64]:
    """``R(mu, E) = sigma_e / sigma_c`` for an LTP=1 section at ``e_ev``: the ratio of the two
    bracketing nodes interpolated linearly in E (exactly the node ratio at a node). Outside the
    tabulated range raises ``ValueError``."""
    en = sec.energies_ev
    if not en[0] <= e_ev <= en[-1]:
        raise ValueError(f"E={e_ev} eV outside the LAW=5 range [{en[0]}, {en[-1]}]")
    j = int(np.clip(np.searchsorted(en, e_ev, side="right") - 1, 0, en.size - 2))
    out = []
    for jj in (j, j + 1):
        rec = sec.records[jj]
        k, eta = coulomb_params(rec.energy_ev, sec.awr, awi, z1, z2)
        se, _, _ = sigma_e_ltp1(mu, rec, sec.lidp, k, eta, sec.spi)
        out.append(se / sigma_c(mu, k, eta, sec.lidp, sec.spi))
    f = (e_ev - en[j]) / (en[j + 1] - en[j])
    return np.asarray((1.0 - f) * out[0] + f * out[1])


def ni_density_ltp1(
    sec: Law5Section, e_ev: float, mu: NDArray[np.float64], awi: float, z1: int, z2: int
) -> NDArray[np.float64]:
    """Transported density ``NI = sigma_c (R - 1)`` [b/sr] at ``e_ev`` (ratio interpolation)."""
    k, eta = coulomb_params(e_ev, sec.awr, awi, z1, z2)
    r = ratio_ltp1(sec, e_ev, mu, awi, z1, z2)
    return sigma_c(mu, k, eta, sec.lidp, sec.spi) * (r - 1.0)


def simpson(y: NDArray[np.float64], x: NDArray[np.float64]) -> float:
    """Composite Simpson rule on a uniform grid with an odd number of points."""
    h = (x[-1] - x[0]) / (x.size - 1)
    return float(h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum()))


def ni_cross_section_b(
    sec: Law5Section, e_ev: float, mu_cut: float, awi: float, z1: int, z2: int, n: int = 20001
) -> float:
    """``2 pi int_0^{mu_cut} NI dmu`` [b]: the half-sphere (each event once) NI cross section of an
    identical-particle LTP=1 section above ``theta_CM = arccos(mu_cut)``."""
    mu = np.linspace(0.0, mu_cut, n)
    return 2.0 * math.pi * simpson(ni_density_ltp1(sec, e_ev, mu, awi, z1, z2), mu)


def ltp12_cross_sections_b(
    sec: Law5Section, index: int, sigma_mf3_b: float, awi: float, z1: int, z2: int
) -> dict[str, float]:
    """For an LTP=12 node at ``sec.records[index]``: ``{"ni_b", "coulomb_b", "total_b",
    "mu_max", "int_p"}`` over ``mu in [-1, mu_max]`` with ``2 pi`` (b); ``ni_b = 2 pi sigma_MF3
    int P dmu`` (P lin-lin) and ``coulomb_b`` is the analytic integral of eq 6.9 (``lidp`` 0)."""
    rec = sec.records[index]
    if rec.ltp != 12 or sec.lidp != 0:
        raise UnsupportedEndfError("ltp12_cross_sections_b needs LTP=12 and LIDP=0")
    mu, p = rec.data[0::2], rec.data[1::2]
    int_p = float(np.sum(0.5 * (p[1:] + p[:-1]) * np.diff(mu)))
    k, eta = coulomb_params(rec.energy_ev, sec.awr, awi, z1, z2)
    mu_max = float(mu[-1])
    coul = 2.0 * math.pi * eta**2 / k**2 * (1.0 / (1.0 - mu_max) - 0.5)
    ni = 2.0 * math.pi * sigma_mf3_b * int_p
    return {"ni_b": ni, "coulomb_b": coul, "total_b": ni + coul, "mu_max": mu_max, "int_p": int_p}
