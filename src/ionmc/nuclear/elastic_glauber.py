"""M2 = optical-limit Glauber shape for the proton elastic table (V3-005D plan section 4.1,
Amendment 20 (b') 1-5). Exploratory model; numpy and stdlib only.

Physics (fm, MeV, hbar c = 197.327 MeV fm; 1 mb = 0.1 fm^2)
-----------------------------------------------------------
* NN profile (per channel ch = pp, pn): ``Gamma(b) = (1 - i alpha) sigma / (4 pi beta)
  exp(-b^2 / (2 beta))`` with sigma, alpha from the GWU SAID SM16 SGT and ALFA tables (linear in
  T_lab, no extrapolation) and ``beta = (1 + alpha^2) sigma / (16 pi)`` (Horiuchi et al. Eq. 10,
  applied per channel, so that sigma_el = sigma_tot). The plan-literal variant (Horiuchi Table I,
  isospin mean, held constant below its lowest energy) is :class:`HoriuchiBeta`.
* Eikonal phase (optical limit, point-nucleon densities rho_p = Z T, rho_n = N T with the same
  normalised thickness function T, ``int T d^2b = 1``)::

      chi(b) = sum_ch n_ch (1/2) sigma_ch (i + alpha_ch) C_ch(b),   C_ch = T (x) f_ch,
      f_ch(b) = exp(-b^2/(2 beta_ch)) / (2 pi beta_ch)        (n_pp = Z, n_pn = N)

  and ``F(q) = i k int b db J0(q b) [1 - exp(i chi(b))]``, ``dsigma/dOmega_CM = |F|^2``,
  ``q = 2 k sin(theta/2)``, ``k = p_CM / hbar c`` of the p + A system (the table's ``p_cm_mev``).
  For a sharp black disk this is ``i k R J1(qR)/q``, i.e. the pinned builder's shape.
* Densities: HO ``rho ~ (1 + alpha r^2/a^2) exp(-r^2/a^2)`` with charge form factor
  ``F_ch = [1 - alpha q^2 a^2 / (2 (2 + 3 alpha))] exp(-q^2 a^2/4)`` (rms^2 = a^2 (6 + 15 alpha) /
  (4 + 6 alpha)); 2pF ``rho ~ 1 / (1 + exp((r - c)/z))``. Proton-size unfolding in form-factor
  space, ``F_pt = F_ch / exp(-q^2 r_p^2/6)``. For HO this is again exactly a HO form with
  ``a'^2 = a^2 - 2 r_p^2/3`` and ``alpha' = 4 c1 / (a'^2 - 6 c1)``, ``c1 = alpha a^2 / (2 (2 + 3
  alpha))`` (rms^2_pt = rms^2_ch - r_p^2 exactly) and the Gaussian-profile convolution is analytic.
  By the projection-slice theorem the 2D Hankel transform of T(b) is F_pt(q), so for 2pF the
  folding ``T (x) f`` is the numerical Hankel quadrature
  ``(1/2 pi) int q dq J0(q b) F_pt(q) exp(-beta q^2 / 2)`` cut at ``q_cut`` (the unfolded integrand
  is exponentially small there; its size at the cut is recorded as ``tail``).
* Coulomb: nuclear amplitude only, exactly as the pinned black-disk shape (Amendment 14 (b)):
  Coulomb stays in the MCS model, no Coulomb-nuclear interference (plan 4.1 item 5).

M2' (``M2Config.pauli``, CLI ``--model m2prime``): the same construction with the Bertulani & De
Conti (PRC 81, 064603) Pauli-blocked in-medium sigma_NN evaluated at the local density in the
eikonal phase (:func:`pauli_factor`, :class:`PauliThickness`); alpha stays the free SAID value and
beta the unitarity rule of the free cross section. ``m2`` itself is unchanged.

Also here: the p-p level re-anchoring above 150 MeV (SAID DSG minus the relativistic Mott cross
section, Amendment 20 (b')4 as revised by the coordinator) and the SAID text readers.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

HBARC = 197.3269804
MP = 938.27208816
ALPHA_EM = 1.0 / 137.035999084
R_PROTON_FM = 0.8414
MB_TO_FM2 = 0.1
DATA_DIR = Path(__file__).resolve().parents[1] / "data"
INPUTS_JSON = DATA_DIR / "elastic_m2_inputs.json"
HORIUCHI_JSON = DATA_DIR / "elastic_m2_horiuchi_table.json"
F64 = NDArray[np.float64]
CF64 = NDArray[np.complex128]


class M2Error(RuntimeError):
    """Fail-closed error of the M2 construction (plan 4.2 trigger 1 or 2)."""


# ---------------------------------------------------------------------------------------------
# Bessel J0 (Abramowitz & Stegun 9.4.1 / 9.4.3, |error| < 5e-8), numpy only
# ---------------------------------------------------------------------------------------------
def _horner(coeffs: tuple[float, ...], t: F64) -> F64:
    """``sum_k coeffs[k] t^k`` (Horner); a loop, not a nested expression (mypy/numpy overload
    resolution is exponential in the nesting depth of a one-line polynomial)."""
    acc: F64 = np.full_like(t, coeffs[-1])
    for c in reversed(coeffs[:-1]):
        acc = acc * t + c
    return acc


_J0_SMALL = (1.0, -2.2499997, 1.2656208, -0.3163866, 0.0444479, -0.0039444, 0.00021)
_J0_F0 = (0.79788456, -0.00000077, -0.00552740, -0.00009512, 0.00137237, -0.00072805, 0.00014476)
_J0_TH = (
    -0.78539816,
    -0.04166397,
    -0.00003954,
    0.00262573,
    -0.00054125,
    -0.00029333,
    0.00013558,
)  # theta0 = x + sum_k c_k / x^k with c_0 the constant (k = 0 term is not divided)


def bessel_j0(x: F64) -> F64:
    """J0 by Abramowitz & Stegun 9.4.1 (|x| <= 3) and 9.4.3 (|x| > 3), |error| < 5e-8."""
    xa = np.abs(np.asarray(x, dtype=np.float64))
    small = xa <= 3.0
    xs = np.where(small, 3.0, xa)  # safe argument of the large-x branch
    ps = _horner(_J0_SMALL, (np.where(small, xa, 3.0) / 3.0) ** 2)
    y = 3.0 / xs
    f0 = _horner(_J0_F0, y)
    th = xs + _horner(_J0_TH, y)
    pl = f0 * np.cos(th) / np.sqrt(xs)
    return np.asarray(np.where(small, ps, pl))


# ---------------------------------------------------------------------------------------------
# SAID text readers
# ---------------------------------------------------------------------------------------------
_NUM = r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?"
_ROW3 = re.compile(rf"^\s*({_NUM})\s+({_NUM})\s+({_NUM})\s*$")
_ROW5 = re.compile(rf"^\s*({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s+({_NUM})\s*$")
_OBS = re.compile(r"^\s*(PP|NP)\s+(SGT|ALFA)\b")


@dataclass(frozen=True)
class SaidSeries:
    """One observable of one channel versus T_lab (prediction block only)."""

    channel: str
    observable: str
    t_mev: F64
    value: F64
    headers: tuple[str, ...]
    sha256: tuple[str, ...]

    def at(self, t: float) -> float:
        """Linear interpolation in T_lab; outside the tabulated range is an error (no
        extrapolation)."""
        if not self.t_mev[0] - 1e-9 <= t <= self.t_mev[-1] + 1e-9:
            raise M2Error(
                f"{self.channel} {self.observable}: T_lab {t} outside "
                f"[{self.t_mev[0]}, {self.t_mev[-1]}]"
            )
        return float(np.interp(t, self.t_mev, self.value))


def parse_said_prediction(text: str) -> tuple[str, str, str, F64, F64]:
    """``(header, channel, observable, T_lab, value)`` of the SAID forward-observable output.

    Format: a header line (``SM16 0-3.0 GEV ...``), a line ``PP SGT at ACM = 0.00 <date>``, a line
    ``TLAB Obs,Err`` and rows ``T_lab value error``. The prediction block ends at the first line
    that is not a three-number row (experimental data blocks that follow are ignored)."""
    lines = [ln.rstrip() for ln in text.splitlines()]
    header, i = "", 0
    for i, ln in enumerate(lines):  # noqa: B007
        if ln.strip() and not ln.lstrip().startswith(("<", "&")):
            header = ln.strip()
            break
    obs = None
    for j in range(i + 1, len(lines)):
        obs = _OBS.match(lines[j])
        if obs:
            start = j + 1
            break
    if obs is None:
        raise M2Error("no 'PP|NP SGT|ALFA' line in the SAID output")
    if "TLAB" in lines[start]:
        start += 1
    t_l: list[float] = []
    v_l: list[float] = []
    for ln in lines[start:]:
        m = _ROW3.match(ln)
        if not m:
            break
        t_l.append(float(m.group(1)))
        v_l.append(float(m.group(2)))
    if len(t_l) < 2:
        raise M2Error("fewer than two prediction rows in the SAID output")
    return header, obs.group(1), obs.group(2), np.asarray(t_l), np.asarray(v_l)


def read_said_series(paths: list[Path]) -> SaidSeries:
    """Merge several output files of the same channel/observable (e.g. 5-30 step 1 and 5-300 step
    5); rows at equal T_lab must agree to 1e-9 relative. The sha256 of every file read is kept."""
    ts: dict[float, float] = {}
    chan = obs_name = None
    heads, shas = [], []
    for p in paths:
        raw = Path(p).read_bytes()
        shas.append(hashlib.sha256(raw).hexdigest())
        h, c, o, t, v = parse_said_prediction(raw.decode("ascii", "replace"))
        if chan is None:
            chan, obs_name = c, o
        elif (c, o) != (chan, obs_name):
            raise M2Error(f"{p}: {c} {o} mixed with {chan} {obs_name}")
        heads.append(h)
        for ti, vi in zip(t, v, strict=True):
            key = round(float(ti), 6)
            if key in ts and abs(ts[key] - vi) > 1e-9 * max(abs(vi), 1e-30) + 1e-12:
                raise M2Error(f"{p}: inconsistent duplicate row at T_lab = {ti}")
            ts[key] = float(vi)
    assert chan is not None and obs_name is not None
    order = sorted(ts)
    return SaidSeries(chan, obs_name, np.asarray(order), np.asarray([ts[k] for k in order]),
                      tuple(heads), tuple(shas))  # fmt: skip


def parse_said_dsg(text: str) -> tuple[F64, F64]:
    """``(A_cm [deg], DSG_cm)`` of a SAID DSG output segment (columns ``A(cm) DSG(cm) Err A(lab)
    DSG(lab)``); every five-number row is taken. Units are not stated by SAID (assumed mb/sr and
    verified in :func:`dsg_anchor`)."""
    a, d = [], []
    for ln in text.splitlines():
        m = _ROW5.match(ln)
        if m:
            a.append(float(m.group(1)))
            d.append(float(m.group(2)))
    if not a:
        raise M2Error("no DSG rows found")
    return np.asarray(a), np.asarray(d)


def merge_dsg_segments(texts: list[str]) -> tuple[F64, F64]:
    """Merge the three 1-degree segments (1-60, 61-120, 121-179) of one energy; the angles must be
    exactly 1..179 degrees in the union."""
    al, dl = [], []
    for t in texts:
        a, d = parse_said_dsg(t)
        al.append(a)
        dl.append(d)
    a, idx = np.unique(np.concatenate(al), return_index=True)
    d = np.concatenate(dl)[idx]
    if a.size != 179 or abs(a[0] - 1.0) > 1e-9 or abs(a[-1] - 179.0) > 1e-9:
        raise M2Error(f"DSG union has {a.size} angles ({a[0]}..{a[-1]}), expected 1..179")
    return a, d


# ---------------------------------------------------------------------------------------------
# NN inputs
# ---------------------------------------------------------------------------------------------
@dataclass
class HoriuchiBeta:
    """Plan-literal beta_NN (sensitivity ``--beta horiuchi-table``): the isospin-mean Table I of
    Horiuchi et al. (PRC 75, 044607), linear in T_lab, held constant below the lowest and above the
    highest tabulated energy. Values come from ``elastic_m2_horiuchi_table.json``."""

    energy_mev: F64
    beta_fm2: F64
    provenance: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def load(cls, path: Path = HORIUCHI_JSON) -> HoriuchiBeta:
        doc = json.loads(Path(path).read_text(encoding="utf-8"))
        tab = doc["table"]
        e, b = np.asarray(tab["energy_mev"], float), np.asarray(tab["beta_fm2"], float)
        if e.size < 2 or e.size != b.size or np.any(np.diff(e) <= 0) or np.any(b <= 0):
            raise M2Error(f"{path}: Horiuchi Table I values are not transcribed (or invalid)")
        return cls(e, b, {k: v for k, v in doc.items() if k != "table"})

    def beta(self, t: float) -> float:
        return float(np.interp(t, self.energy_mev, self.beta_fm2))


@dataclass
class NNInputs:
    """Per-channel sigma_NN, alpha_NN, beta_NN versus T_lab."""

    series: dict[tuple[str, str], SaidSeries]
    horiuchi: HoriuchiBeta | None = None

    @classmethod
    def from_files(
        cls, files: dict[str, list[Path]], horiuchi: HoriuchiBeta | None = None
    ) -> NNInputs:
        """``files`` keys ``pp_sgt``, ``pp_alfa``, ``np_sgt``, ``np_alfa`` -> list of files."""
        out: dict[tuple[str, str], SaidSeries] = {}
        for key, paths in files.items():
            s = read_said_series([Path(p) for p in paths])
            if (s.channel.lower(), s.observable.lower()) != tuple(key.split("_")):
                raise M2Error(f"{key}: files hold {s.channel} {s.observable}")
            out[(s.channel, s.observable)] = s
        for ch in ("PP", "NP"):
            for ob in ("SGT", "ALFA"):
                if (ch, ob) not in out:
                    raise M2Error(f"missing SAID series {ch} {ob}")
        return cls(out, horiuchi)

    @property
    def t_min_mev(self) -> float:
        """Lowest T_lab covered by every series."""
        return max(float(s.t_mev[0]) for s in self.series.values())

    def sigma_fm2(self, ch: str, t: float) -> float:
        return MB_TO_FM2 * self.series[(ch, "SGT")].at(t)

    def alpha(self, ch: str, t: float) -> float:
        return self.series[(ch, "ALFA")].at(t)

    def beta_fm2(self, ch: str, t: float) -> float:
        """``(1 + alpha^2) sigma / (16 pi)`` per channel (or the Horiuchi table value)."""
        if self.horiuchi is not None:
            return self.horiuchi.beta(t)
        a = self.alpha(ch, t)
        return (1.0 + a * a) * self.sigma_fm2(ch, t) / (16.0 * math.pi)

    def provenance(self) -> dict[str, Any]:
        return {f"{c}_{o}": {"headers": sorted(set(s.headers)), "sha256": list(s.sha256),
                             "t_range_mev": [float(s.t_mev[0]), float(s.t_mev[-1])],
                             "n_rows": int(s.t_mev.size)}
                for (c, o), s in sorted(self.series.items())}  # fmt: skip


# ---------------------------------------------------------------------------------------------
# densities
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class HODensity:
    """Harmonic-oscillator form ``rho ~ (1 + alpha r^2/a^2) exp(-r^2/a^2)`` (a in fm)."""

    alpha: float
    a_fm: float

    def coeff(self) -> float:
        """``c1`` of ``F(q) = (1 - c1 q^2) exp(-q^2 a^2/4)``."""
        return self.alpha * self.a_fm**2 / (2.0 * (2.0 + 3.0 * self.alpha))

    def rms2(self) -> float:
        return self.a_fm**2 * (6.0 + 15.0 * self.alpha) / (4.0 + 6.0 * self.alpha)

    def form_factor(self, q: F64) -> F64:
        q2 = np.asarray(q) ** 2
        return np.asarray((1.0 - self.coeff() * q2) * np.exp(-q2 * self.a_fm**2 / 4.0))

    def unfolded(self, r_p: float = R_PROTON_FM) -> HODensity:
        """Point density: ``F_pt = F_ch / exp(-q^2 r_p^2/6)`` is again of the HO form."""
        c1 = self.coeff()
        a2 = self.a_fm**2 - 2.0 * r_p**2 / 3.0
        if a2 <= 6.0 * c1:
            raise M2Error("proton-size unfolding of the HO density is not of HO form (a'^2<=6c1)")
        return HODensity(4.0 * c1 / (a2 - 6.0 * c1), math.sqrt(a2))

    def density(self, r: F64) -> F64:
        """Normalised ``int rho d^3r = 1``."""
        a, al = self.a_fm, self.alpha
        n = 1.0 / (math.pi**1.5 * a**3 * (1.0 + 1.5 * al))
        return np.asarray(n * (1.0 + al * (r / a) ** 2) * np.exp(-((r / a) ** 2)))

    def gauss_folded_thickness(self, b: F64, beta: float) -> F64:
        """``T (x) f`` [fm^-2] for the Gaussian profile ``f = exp(-b^2/2beta)/(2 pi beta)``:
        ``(1/4 pi u) exp(-b^2/4u) [1 + c1 (b^2/(4u^2) - 1/u)]``, ``u = a^2/4 + beta/2``
        (``beta = 0`` gives T itself)."""
        u = self.a_fm**2 / 4.0 + beta / 2.0
        b2 = np.asarray(b) ** 2
        c1 = self.coeff()
        return np.asarray(np.exp(-b2 / (4.0 * u)) / (4.0 * math.pi * u)
                          * (1.0 + c1 * (b2 / (4.0 * u * u) - 1.0 / u)))  # fmt: skip


@dataclass(frozen=True)
class FermiDensity:
    """Two-parameter Fermi ``rho ~ 1/(1 + exp((r - c)/z))`` (fm); with ``w != 0`` the
    three-parameter Fermi ``rho ~ (1 + w r^2/c^2)/(1 + exp((r - c)/z))``."""

    c_fm: float
    z_fm: float
    w: float = 0.0

    def _rho(self, r: F64) -> F64:
        return np.asarray((1.0 + self.w * (r / self.c_fm) ** 2)
                          / (1.0 + np.exp((r - self.c_fm) / self.z_fm)))  # fmt: skip

    def _r(self, dr: float = 0.01) -> F64:
        return np.arange(0.0, self.c_fm + 40.0 * self.z_fm, dr)

    def rms2(self) -> float:
        r = self._r()
        w = self._rho(r)
        return float(np.sum(w * r**4) / np.sum(w * r**2))

    def density(self, r: F64) -> F64:
        """Normalised ``int rho d^3r = 1`` (numerical normalisation on the internal grid)."""
        rg = self._r()
        norm = float(np.trapezoid(4.0 * np.pi * rg**2 * self._rho(rg), rg))
        return np.asarray(self._rho(np.asarray(r, dtype=np.float64)) / norm)

    def form_factor(self, q: F64, dr: float = 0.01) -> F64:
        r = self._r(dr)
        w = self._rho(r) * r  # rho r
        q = np.asarray(q, dtype=np.float64)
        out = np.empty_like(q)
        norm = float(np.sum(w * r))  # int rho r^2 dr
        for i0 in range(0, q.size, 256):
            qq = q[i0 : i0 + 256, None]
            safe = np.where(qq < 1e-12, 1.0, qq)
            kern = np.where(qq < 1e-12, r[None, :], np.sin(safe * r[None, :]) / safe)
            ff = np.sum(w[None, :] * kern,
                        axis=1)  # fmt: skip
            out[i0 : i0 + 256] = ff / norm
        return out


def scale_density(d: HODensity | FermiDensity, factor: float) -> HODensity | FermiDensity:
    """Surrogate mapping of decision 0041: all lengths scaled by ``(A/A_ref)^(1/3)``."""
    if isinstance(d, HODensity):
        return HODensity(d.alpha, d.a_fm * factor)
    return FermiDensity(d.c_fm * factor, d.z_fm * factor, d.w)


def rms_matched_point(dens: FermiDensity, r_p: float = R_PROTON_FM) -> FermiDensity:
    """Fermi point density of the same z and w as ``dens`` whose size is solved so that
    ``rms_pt^2 = rms_ch^2 - r_p^2`` (the ``rms-matched`` unfolding of :class:`Thickness`)."""
    target = math.sqrt(dens.rms2() - r_p**2)
    c = solve_size(lambda x: FermiDensity(x, dens.z_fm, dens.w), target, 0.8 * dens.c_fm,
                   1.2 * dens.c_fm)  # fmt: skip
    return FermiDensity(c, dens.z_fm, dens.w)


class Thickness:
    """Point-nucleon thickness function folded with the NN Gaussian, on a fixed b grid.

    ``conv(beta)`` returns ``C(b) = (T (x) f_beta)(b)`` [fm^-2]. HO: analytic. 2pF: Hankel
    quadrature of ``F_pt(q) exp(-beta q^2/2)`` with ``q`` up to ``q_cut`` (trapezoid, ``dq``)."""

    def __init__(self, dens: HODensity | FermiDensity, b: F64, r_p: float = R_PROTON_FM,
                 q_cut: float = 24.0, dq: float = 0.02, unfold: str = "formfactor",
                 dr: float = 0.01) -> None:  # fmt: skip
        """``unfold`` (Fermi-type densities only): ``formfactor`` = the plan's
        ``F_pt = F_ch / exp(-q^2 r_p^2/6)`` (exact, but the integrand ``F_pt exp(-beta q^2/2)``
        diverges at large q whenever ``beta < r_p^2/3 = 0.236 fm^2``, i.e. above ~125 MeV; the
        ``tail`` attribute records the integrand at ``q_cut``); ``rms-matched`` = a Fermi point
        density of the same z and w whose size is solved so that ``rms_pt^2 = rms_ch^2 - r_p^2``
        (well posed; NOT the plan's rule, exploratory only)."""
        self.dens, self.b, self.r_p, self.unfold = dens, b, r_p, unfold
        self.dq = dq
        self.tail = 0.0
        if isinstance(dens, HODensity):
            self.point: HODensity | None = dens.unfolded(r_p)
        else:
            self.point = None
            self.q = np.arange(0.0, q_cut + 0.5 * dq, dq)
            if unfold == "formfactor":
                f_pt = dens.form_factor(self.q, dr) * np.exp(self.q**2 * r_p**2 / 6.0)
            elif unfold == "rms-matched":
                target = math.sqrt(dens.rms2() - r_p**2)
                c = solve_size(
                    lambda x: FermiDensity(x, dens.z_fm, dens.w),
                    target,
                    0.8 * dens.c_fm,
                    1.2 * dens.c_fm,
                )
                self.point_fermi = FermiDensity(c, dens.z_fm, dens.w)
                f_pt = self.point_fermi.form_factor(self.q, dr)
            else:
                raise M2Error(f"unknown unfold mode {unfold!r}")
            self.f_pt = f_pt
            w = np.full(self.q.size, dq)
            w[0] = w[-1] = 0.5 * dq
            self.kern = bessel_j0(np.outer(b, self.q)) * (self.q * w / (2.0 * math.pi))

    def conv(self, beta: float) -> F64:
        if self.point is not None:
            return self.point.gauss_folded_thickness(self.b, beta)
        env = self.f_pt * np.exp(-beta * self.q**2 / 2.0)
        self.tail = max(self.tail, float(np.max(np.abs(env[-max(2, self.q.size // 20) :]))))
        # Euler-Maclaurin end-point term of the q trapezoid: the integrand q env J0(q b) has the
        # slope env(0) at q = 0, so the sum is low by dq^2 env(0) / (24 pi) for every b
        return np.asarray(self.kern @ env + self.dq**2 * env[0] / (24.0 * math.pi))


RHO0_FM3 = 0.17
"""Saturation density rho_0 of the Bertulani & De Conti medium modification [fm^-3]."""


PAULI_VARIANTS = ("clementel-villi", "bdc-pn")
HB2M_MEV_FM2 = HBARC**2 / (2.0 * MP)
"""``hbar^2 / 2 m_N`` [MeV fm^2] (``m_N`` = proton mass)."""


def fermi_energy_mev(rho_tot: F64) -> F64:
    """``E_F = (hbar^2 / 2 m_N) (3 pi^2 rho_tot / 2)^(2/3)`` [MeV] of symmetric nuclear matter of
    total (point-nucleon) density ``rho_tot`` [fm^-3] (``k_F = (3 pi^2 rho_tot/2)^(1/3)``)."""
    kf = (1.5 * math.pi**2 * np.maximum(np.asarray(rho_tot, dtype=np.float64), 0.0)) ** (1.0 / 3.0)
    return np.asarray(HB2M_MEV_FM2 * kf * kf)


def pauli_factor(
    e_mev: float, rho_p: F64, rho_n: F64, variant: str, scale: float = 1.0
) -> F64:  # fmt: skip
    """In-medium (Pauli-blocked) reduction ``sigma_eff(E, rho) / sigma_free(E)`` of a bare nucleon
    of laboratory kinetic energy ``e_mev`` [MeV] on a local nucleon density ``rho_p``, ``rho_n``
    [fm^-3]. The reading is an explicit REQUIRED parameter ``variant``:

    ``clementel-villi`` (the variant of record): flux-weighted Pauli blocking of a bare nucleon
    scattering on a nucleon Fermi gas (symmetric matter), in the local-density approximation.
    The closed form below was verified by an independent Monte Carlo derivation (session record
    M2PRIME.md, "R-C verification"); no printed open source was located, so no equation
    citation is given. With ``x = E_F / E``, ``E_F = (hbar c)^2 k_F^2 / (2 m_p c^2)``,
    ``k_F = (3 pi^2 rho_tot / 2)^(1/3)``, ``rho_tot = rho_p + rho_n`` the local point-nucleon
    density [fm^-3] (``E_F(0.17) = 38.392 MeV``), and ``E`` the bare laboratory kinetic energy::

        E >= 2 E_F :  1 - 7 x / 5
        E_F < E < 2 E_F :  1 - 7 x / 5 + (2 x / 5) (2 - E/E_F)^(5/2)
        E <= E_F :  0

    (both branches equal 0.3 at E = 2 E_F, so the factor is continuous; the same factor is used
    for the pp and the pn channel; ``rho_p``/``rho_n`` enter only through the sum).

    ``bdc-pn`` (report-only sensitivity): Bertulani & De Conti, PRC 81, 064603 (2010),
    arXiv:1004.2096, Eq. (14) with ``rho1 = rho_p``, ``rho2 = rho_n`` (local densities)::

        [1 + 1.892 (|rho1 - rho2| / (rho~ rho0))^2.75]^-1 *
        { 1 - 37.02 rho~^(2/3) / E        if E > 46.27 rho~^(2/3)
          E / (231.38 rho~^(2/3))         otherwise }

    ``rho~ = (rho1 + rho2) / rho0``, ``rho0 = 0.17 fm^-3`` (zero asymmetry term for N = Z; the
    published constants leave a 6e-5 step between the branches).

    ``scale`` multiplies the densities (a test hook: ``scale -> 0`` recovers the free cross
    section, i.e. M2). The single-Fermi-sphere reading of Eq. (14) (rho1 = 0) is deliberately NOT
    implemented (no free limit)."""
    if variant not in PAULI_VARIANTS:
        raise M2Error(f"unknown Pauli variant {variant!r}; expected one of {PAULI_VARIANTS}")
    rp = np.asarray(rho_p, dtype=np.float64)
    rn = np.asarray(rho_n, dtype=np.float64)
    tot = np.maximum(rp + rn, 1e-300)
    if variant == "clementel-villi":
        ef = np.asarray(fermi_energy_mev(scale * tot))
        x = ef / e_mev
        hi = 1.0 - 1.4 * x
        arg = np.clip(2.0 - 1.0 / np.where(x > 0.0, x, 1.0), 0.0, None)
        lo = hi + 0.4 * x * arg**2.5
        return np.asarray(np.where(e_mev <= ef, 0.0, np.where(e_mev >= 2.0 * ef, hi, lo)))
    rt23 = np.maximum(scale * tot / RHO0_FM3, 1e-300) ** (2.0 / 3.0)
    asym = 1.0 / (1.0 + 1.892 * (np.abs(rp - rn) / tot) ** 2.75)
    base = np.where(e_mev > 46.27 * rt23, 1.0 - 37.02 * rt23 / e_mev, e_mev / (231.38 * rt23))
    return np.asarray(asym * base)


class PauliThickness:
    """M2' (Pauli-blocked) point-nucleon thickness, folded with the NN Gaussian, on a fixed b grid.

    The phase is the z integral of the local-density profile with ``sigma_NN`` evaluated at the
    local density: ``chi = sum_ch n_ch (1/2) sigma_ch (i + alpha_ch) (T_eff (x) f_ch)`` with the
    effective density ``rho_eff(r) = rho_pt(r) P(E, Z rho_pt(r), N rho_pt(r))`` (``P`` =
    :func:`pauli_factor`, not renormalised; ``alpha`` and ``beta`` are the free-NN values). By the
    projection-slice theorem the Hankel transform of ``T_eff`` is ``F_eff(q) = int rho_eff j0(q r)
    d^3r``, so the folding is the same Hankel quadrature as :class:`Thickness`. The point density
    is the unfolded HO density, or the rms-matched Fermi density (no form-factor unfolding here:
    a position-space point density is needed)."""

    def __init__(self, dens: HODensity | FermiDensity, b: F64, z: int, n: int, *, variant: str,
                 r_p: float = R_PROTON_FM, q_cut: float = 24.0, dq: float = 0.02,
                 dr: float = 0.01, pauli_scale: float = 1.0) -> None:  # fmt: skip
        self.b, self.z, self.n, self.scale, self.variant = b, z, n, pauli_scale, variant
        self.dr, self.dq = dr, dq
        self.tail = 0.0
        self.pauli_peak = 1.0
        self.point: HODensity | FermiDensity = (
            dens.unfolded(r_p) if isinstance(dens, HODensity) else rms_matched_point(dens, r_p)
        )
        self.r = np.arange(0.0, 16.0 * dens.a_fm if isinstance(dens, HODensity)
                           else dens.c_fm + 40.0 * dens.z_fm, dr)  # fmt: skip
        self.rho_pt = self.point.density(self.r)
        self.q = np.arange(0.0, q_cut + 0.5 * dq, dq)
        w = np.full(self.q.size, dq)
        w[0] = w[-1] = 0.5 * dq
        self.kern = bessel_j0(np.outer(b, self.q)) * (self.q * w / (2.0 * math.pi))
        self._e = float("nan")
        self.f_eff: F64 = np.zeros(0)

    def set_energy(self, e_mev: float) -> None:
        """Compute ``F_eff(q)`` for the laboratory energy ``e_mev``."""
        pf = pauli_factor(
            e_mev, self.z * self.rho_pt, self.n * self.rho_pt, self.variant, self.scale
        )
        self.pauli_peak = float(pf[0])  # the density peaks at r = 0 for the HO/Fermi forms used
        w = self.rho_pt * pf * self.r
        f = np.empty(self.q.size)
        for i0 in range(0, self.q.size, 128):
            qq = self.q[i0 : i0 + 128, None]
            safe = np.where(qq < 1e-12, 1.0, qq)
            kern = np.where(qq < 1e-12, self.r[None, :], np.sin(safe * self.r[None, :]) / safe)
            f[i0 : i0 + 128] = 4.0 * math.pi * np.sum(w[None, :] * kern, axis=1) * self.dr
        self.f_eff = f
        self._e = e_mev

    def conv(self, beta: float) -> F64:
        env = self.f_eff * np.exp(-beta * self.q**2 / 2.0)
        self.tail = max(self.tail, float(np.max(np.abs(env[-max(2, self.q.size // 20) :]))))
        return np.asarray(self.kern @ env + self.dq**2 * env[0] / (24.0 * math.pi))  # see Thickness


def b_grid(dens: HODensity | FermiDensity, db: float = 0.025) -> F64:
    """Uniform b grid ``0..b_max`` with ``b_max`` where the folded thickness is < ~1e-11."""
    if isinstance(dens, HODensity):
        bmax = 6.5 * math.sqrt(dens.a_fm**2 + 2.0)
    else:
        bmax = dens.c_fm + 26.0 * dens.z_fm
    return np.arange(0.0, bmax + 0.5 * db, db)


# ---------------------------------------------------------------------------------------------
# eikonal phase, amplitude, shape
# ---------------------------------------------------------------------------------------------
def eikonal_chi(th: Thickness | PauliThickness, nn: NNInputs, t_mev: float, z: int, n: int) -> CF64:
    """``chi(b)`` (complex) of the optical limit, summed over the pp (Z protons) and pn (N
    neutrons) channels at the NN laboratory energy ``t_mev``."""
    chi = np.zeros(th.b.size, dtype=np.complex128)
    for ch, cnt in (("PP", z), ("NP", n)):
        if cnt == 0:
            continue
        sig, al, be = nn.sigma_fm2(ch, t_mev), nn.alpha(ch, t_mev), nn.beta_fm2(ch, t_mev)
        chi += cnt * 0.5 * sig * (1j + al) * th.conv(be)
    return chi


def trapz_weights(b: F64) -> F64:
    w = np.full(b.size, float(b[1] - b[0]))
    w[0] = w[-1] = 0.5 * w[0]
    return w


def eikonal_amplitude(b: F64, chi: CF64, k_fm: float, s: F64) -> NDArray[np.complex128]:
    """``F(q) = i k int b db J0(q b) [1 - exp(i chi)]`` at ``q = 2 k s`` [fm]."""
    g = (1.0 - np.exp(1j * chi)) * b * trapz_weights(b)
    out = np.empty(s.size, dtype=np.complex128)
    for i0 in range(0, s.size, 1024):
        q = 2.0 * k_fm * s[i0 : i0 + 1024]
        out[i0 : i0 + 1024] = 1j * k_fm * (bessel_j0(np.outer(q, b)) @ g)
    return out


@dataclass(frozen=True)
class ShapeResult:
    s: F64
    density: F64  # |F|^2 [fm^2/sr] on the s grid
    sigma_el_fm2: float  # int |F|^2 dOmega over the full sphere
    sigma_tot_ot_fm2: float  # (4 pi/k) Im F(0)
    sigma_tot_b_fm2: float  # 4 pi int b db Re(1 - e^{i chi}) (direct)
    sigma_react_fm2: float  # 2 pi int b db (1 - |e^{i chi}|^2)


def m2_shape(b: F64, chi: CF64, k_fm: float, n_s: int) -> ShapeResult:
    s = np.linspace(0.0, 1.0, n_s)
    f = eikonal_amplitude(b, chi, k_fm, s)
    dens = np.abs(f) ** 2
    g = 4.0 * s * dens  # per ds; dmu = 4 s ds (up to sign), dOmega = 2 pi dmu
    sig_el = 2.0 * math.pi * float(np.sum(0.5 * (g[1:] + g[:-1]) * np.diff(s)))
    w = trapz_weights(b) * b
    e = np.exp(1j * chi)
    return ShapeResult(
        s, dens, sig_el, 4.0 * math.pi / k_fm * float(f[0].imag),
        4.0 * math.pi * float(np.sum(w * (1.0 - e).real)),
        2.0 * math.pi * float(np.sum(w * (1.0 - np.abs(e) ** 2))),
    )  # fmt: skip


def quantiles_mu(s: F64, dens: F64, n_q: int) -> tuple[F64, float]:
    """Edges (n_q + 1, ascending mu_CM from -1 to 1) of the ``n_q`` equiprobable bins of the
    density ``dens(s)`` (per solid angle) with ``mu = 1 - 2 s^2``: the same convention and
    machinery as the pinned builder's ``disk_quantiles`` (CDF = trapezoid sum of ``4 s dens`` on
    the s grid, inverted linearly in s). Also returns ``|cumulative - 1|`` after normalisation
    (the normalisation check of the discrete CDF)."""
    g = 4.0 * s * dens
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(s))))
    cum /= cum[-1]
    target = 1.0 - np.arange(n_q + 1) / n_q
    s_q = np.interp(target, cum, s)
    mu = 1.0 - 2.0 * s_q * s_q
    mu[0], mu[-1] = -1.0, 1.0
    return np.maximum.accumulate(mu), float(abs(cum[-1] - 1.0))


def helm_q_h(dens: HODensity, r_p: float = R_PROTON_FM) -> dict[str, float]:
    """Plan 4.2 trigger 4: ``q_H = 4.493 / R_H`` of the Helm form factor whose ``R`` is matched to
    the q^4 cumulant of the point density (``kappa = <r^4> - (5/3) <r^2>^2``, invariant under
    Gaussian folding; Helm: ``kappa = -6 R^4 / 35``) and whose surface width follows from
    ``<r^2> = 3 R^2/5 + 3 s^2`` (equivalently ``s_pt^2 = s_ch^2 - r_p^2/3``). HO densities only."""
    p = dens.unfolded(r_p)
    r = np.linspace(0.0, 16.0 * p.a_fm, 200001)
    w = 4.0 * np.pi * r * r * p.density(r)
    r2 = float(np.trapezoid(w * r**2, r))
    kap = float(np.trapezoid(w * r**4, r)) - 5.0 / 3.0 * r2 * r2
    r_h = (-35.0 * kap / 6.0) ** 0.25
    return {
        "r_h_fm": r_h,
        "s_pt_fm": math.sqrt((r2 - 0.6 * r_h * r_h) / 3.0),
        "q_h_fm_inv": 4.493 / r_h,
    }


def first_minimum_q(s: F64, dens: F64, k_fm: float) -> float:
    """q [fm^-1] of the first local minimum of the density (plan 4.2 trigger 4); nan if none."""
    d = dens
    for i in range(1, d.size - 1):
        if d[i] < d[i - 1] and d[i] <= d[i + 1]:
            return float(2.0 * k_fm * s[i])
    return float("nan")


# ---------------------------------------------------------------------------------------------
# p-p level re-anchoring above 150 MeV
# ---------------------------------------------------------------------------------------------
def mott_pp_endf_mb_sr(t_lab_mev: float, theta_cm_deg: F64, awi: float) -> F64:
    """Identical spin-1/2 Coulomb (Mott) density of ENDF-102 eq 6.10 in mb/sr, the very function
    (:func:`ionmc.nuclear.law5.sigma_c`, ``coulomb_params`` with the MF1 ``awi``) that the LAW=5
    reconstruction subtracts to form NI::

        (eta / 2k)^2 [1/sin^4 + 1/cos^4 - cos(eta ln tan^2(theta/2)) / (sin^2 cos^2)]
        = 2 eta^2/(k^2 (1 - mu^2)) [(1 + mu^2)/(1 - mu^2) - cos(eta ln((1+mu)/(1-mu)))/2]

    (``theta`` = theta_CM; the interference coefficient is ``(-1)^{2s}/(2s+1) = -1/2``; the sign of
    the cos term in the amendment's printed 'half' form is that of this form), ``k`` and ``eta`` of
    ENDF-102 eqs 6.11/6.12 (non-relativistic, ``eta = alpha/v_rel``)."""
    from ionmc.nuclear import law5

    k, eta = law5.coulomb_params(t_lab_mev * 1e6, awi, awi, 1, 1)
    mu = np.cos(np.radians(np.asarray(theta_cm_deg, dtype=np.float64)))
    return np.asarray(1.0e3 * law5.sigma_c(mu, k, eta, 1, 0.5))


def mott_pp_moller_mb_sr(t_lab_mev: float, theta_cm_deg: F64) -> F64:
    """Cross-reference only: relativistic Born (Moller) spin-averaged identical-fermion Coulomb
    cross section, CM, mb/sr: ``alpha^2 (hbar c)^2 (2E^2 - m^2)^2 / (4 E^2 p^4) [4/sin^4 - 3/sin^2
    + (p^2/(2E^2 - m^2))^2 (1 + 4/sin^2)]`` with ``E = sqrt(s)/2`` (no Coulomb phase). Its
    non-relativistic limit is :func:`mott_pp_endf_mb_sr` at ``theta = 90 deg`` (phase = 1)."""
    s_mand = 2.0 * MP * MP + 2.0 * MP * (t_lab_mev + MP)
    e = 0.5 * math.sqrt(s_mand)
    p2 = e * e - MP * MP
    sin2 = np.sin(np.radians(np.asarray(theta_cm_deg, dtype=np.float64))) ** 2
    pref = (ALPHA_EM * HBARC) ** 2 * (2 * e * e - MP * MP) ** 2 / (4.0 * e * e * p2 * p2)
    br = 4.0 / sin2**2 - 3.0 / sin2 + (p2 / (2 * e * e - MP * MP)) ** 2 * (1.0 + 4.0 / sin2)
    return np.asarray(pref * br * 10.0)


def dsg_anchor(
    t_mev: float, a_cm_deg: F64, dsg_mb_sr: F64, sigma_sgt_mb: float, awi: float,
    mu_cut: float = 0.96,
) -> dict[str, Any]:  # fmt: skip
    """``sigma_anchor = 2 pi int_0^{mu_cut} [DSG_SAID - Mott_pp] dmu`` [mb] (half sphere, each
    event once = the table's P6 item (3) convention, so it is like-for-like with the table partial
    and with SGT) from the 1-degree DSG scan (mb/sr assumed). ``DSG - Mott`` is interpolated
    linearly in theta_CM in its logarithm (angle grid 0.02 deg); if it is not positive at a node
    the energy is unusable (``usable`` False, ``sigma_anchor_mb`` nan)."""
    diff = dsg_mb_sr - mott_pp_endf_mb_sr(t_mev, a_cm_deg, awi)
    out: dict[str, Any] = {"t_mev": t_mev, "sigma_sgt_mb": sigma_sgt_mb, "usable": bool(
        np.all(diff > 0.0)), "min_diff_mb_sr": float(diff.min())}  # fmt: skip
    if not out["usable"]:
        out.update(sigma_anchor_mb=float("nan"), ratio_to_sgt=float("nan"))
        return out
    th_cut = math.degrees(math.acos(mu_cut))
    th = np.linspace(th_cut, 90.0, int((90.0 - th_cut) / 0.02) + 1)
    rad = np.radians(th)
    f = np.exp(np.interp(th, a_cm_deg, np.log(diff))) * np.sin(rad)
    sig = 2.0 * math.pi * float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(rad)))
    out.update(sigma_anchor_mb=sig, ratio_to_sgt=sig / sigma_sgt_mb)
    return out


CONTINUITY_TOL = 0.03
SGT_RATIO_BAND = (0.80, 1.05)
RAMP_FROM_MEV = 125.0
ANCHOR_FROM_MEV = 150.0


@dataclass
class PPReanchor:
    """Level re-anchoring of the p-p partial cross section (Amendment 20 (b')4, revision 3).

    ``rows``: per scan energy the :func:`dsg_anchor` record plus ``sigma_table_mb`` (the table's
    partial at that energy) and ``ratio_to_table``. Blocking checks (assertions with recorded
    values): units/frame ``ratio_to_sgt`` in [0.80, 1.05] at every scan energy and usable DSG;
    continuity ``|ratio_to_table - 1| <= 3 %`` at 100 and 125 MeV. If one fails, ``enabled`` is
    False and ``reason`` says which (the table level is then left unchanged)."""

    rows: list[dict[str, Any]]
    enabled: bool
    reason: str

    @classmethod
    def build(cls, rows: list[dict[str, Any]]) -> PPReanchor:
        reasons = []
        for r in rows:
            if not r["usable"]:
                reasons.append(f"DSG - Mott not positive at {r['t_mev']:g} MeV")
            elif not SGT_RATIO_BAND[0] <= r["ratio_to_sgt"] <= SGT_RATIO_BAND[1]:
                reasons.append(f"sigma_anchor/SGT = {r['ratio_to_sgt']:.4f} at {r['t_mev']:g} MeV")
        for r in rows:
            if r["t_mev"] in (100.0, 125.0) and r["usable"]:
                if abs(r["ratio_to_table"] - 1.0) > CONTINUITY_TOL:
                    reasons.append(
                        f"continuity: sigma_anchor/sigma_table = {r['ratio_to_table']:.4f} "
                        f"at {r['t_mev']:g} MeV"
                    )
        return cls(rows, not reasons, "; ".join(reasons))

    def _row(self, t: float) -> dict[str, Any]:
        return next(r for r in self.rows if r["t_mev"] == t)

    def sigma_mb(self, e_mev: float, sigma_table_mb: float) -> float:
        """New partial [mb] at a node: unchanged below 125 MeV (and when disabled); the table
        value times a factor ramping linearly from 1 (125 MeV) to ``r(150)`` below 150 MeV; the
        SAID anchor linear in T_lab between the scan energies for E >= 150 MeV."""
        if not self.enabled or e_mev <= RAMP_FROM_MEV:
            return sigma_table_mb
        if e_mev < ANCHOR_FROM_MEV:
            r150 = self._row(150.0)["ratio_to_table"]
            ramp = (e_mev - RAMP_FROM_MEV) / (ANCHOR_FROM_MEV - RAMP_FROM_MEV)
            return float(sigma_table_mb * (1.0 + (r150 - 1.0) * ramp))
        t = [r["t_mev"] for r in self.rows]
        return float(np.interp(e_mev, t, [r["sigma_anchor_mb"] for r in self.rows]))


# ---------------------------------------------------------------------------------------------
# parameter provenance JSON
# ---------------------------------------------------------------------------------------------
def load_inputs(path: Path = INPUTS_JSON) -> dict[str, Any]:
    return dict(json.loads(Path(path).read_text(encoding="utf-8")))


def density_from_entry(entry: dict[str, Any]) -> HODensity | FermiDensity:
    if entry["form"] == "HO":
        return HODensity(float(entry["alpha"]), float(entry["a_fm"]))
    if entry["form"] in ("2pF", "3pF"):
        return FermiDensity(float(entry["c_fm"]), float(entry["z_fm"]), float(entry.get("w", 0.0)))
    raise M2Error(f"unsupported density form {entry['form']!r}")


def solve_size(make: Any, rms_target: float, lo: float, hi: float) -> float:
    """Bisection for the size parameter ``x`` with ``sqrt(make(x).rms2()) = rms_target``."""
    flo = math.sqrt(make(lo).rms2()) - rms_target
    fhi = math.sqrt(make(hi).rms2()) - rms_target
    if flo * fhi > 0:
        raise M2Error("size parameter not bracketed")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fm = math.sqrt(make(mid).rms2()) - rms_target
        if fm * flo > 0:
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def target_density(
    inputs: dict[str, Any], name: str, shape_override: dict[str, float] | None = None
) -> tuple[HODensity | FermiDensity, dict[str, Any]]:
    """Charge density of target ``name`` and its provenance entry.

    Tabulated targets (``form`` with parameters) are taken as stored. A ``surrogate`` entry
    (IAEA-radius-constrained, coordinator decision) takes the shape parameter from the nearest
    tabulated nucleus (HO: ``alpha`` = mean of the entries in ``alpha_from``; 2pF: ``z`` of
    ``d_from``) and solves the size parameter (HO ``a``, 2pF ``c``) numerically so that the charge
    rms radius equals the IAEA value in ``iaea_rms_charge_radius_fm`` exactly.
    ``shape_override[name]`` replaces the shape parameter (diffuseness-sensitivity option)."""
    ent = inputs["densities"].get(name)
    if ent is None or ent.get("status") == "pending":
        raise M2Error(f"{name}: no density parameters in {INPUTS_JSON.name} (pending acquisition)")
    if not ent.get("surrogate"):
        return density_from_entry(ent), ent
    rms = float(inputs["iaea_rms_charge_radius_fm"]["values"][name])
    ov = (shape_override or {}).get(name)
    if ent["form"] == "HO":
        al = (
            ov
            if ov is not None
            else float(np.mean([inputs["densities"][n]["alpha"] for n in ent["alpha_from"]]))
        )
        a = solve_size(lambda x: HODensity(al, x), rms, 0.5, 5.0)
        return HODensity(al, a), ent
    z = ov if ov is not None else float(inputs["densities"][ent["d_from"]]["z_fm"])
    c = solve_size(lambda x: FermiDensity(x, z), rms, 1.0, 8.0)
    return FermiDensity(c, z), ent


def rms_check(dens: HODensity | FermiDensity, name: str, inputs: dict[str, Any]) -> dict[str, Any]:
    """rms charge radius of the density against the IAEA value read from the provenance JSON
    (tolerance ``rms_tolerance_fm``); returns the recorded values (no hard-coded reference)."""
    ref = inputs["iaea_rms_charge_radius_fm"]["values"].get(name)
    if ref is None:
        raise M2Error(f"{name}: no IAEA rms radius in the provenance JSON")
    tol = float(inputs["iaea_rms_charge_radius_fm"]["tolerance_fm"])
    rms = math.sqrt(dens.rms2())
    return {"rms_model_fm": rms, "rms_iaea_fm": float(ref), "delta_fm": rms - float(ref),
            "tolerance_fm": tol, "pass": bool(abs(rms - float(ref)) <= tol)}  # fmt: skip


# ---------------------------------------------------------------------------------------------
# build driver
# ---------------------------------------------------------------------------------------------
BETA_MODES = ("unitarity", "horiuchi-table")
PP_LEVEL_MODES = ("none", "said-anchored-ge150")
NORM_TOL = 1.0e-12
CONVERGENCE_TOL = 1.0e-3
TAIL_TOL = 1.0e-6
RMS_TOL_SURROGATE_NOTE = "surrogate: identically satisfied by construction, not a pass"


@dataclass(frozen=True)
class M2Config:
    """Options of the M2 build (kept out of ``ElasticBuildOptions``, whose canonical JSON is part
    of the pinned table id)."""

    input_dir: Path
    dsg_files: dict[float, list[Path]] | None = None
    beta: str = "unitarity"  # --beta
    pp_level: str = "none"  # --pp-level
    unfold: str = "rms-matched"  # declared rule for 2pF/3pF; "formfactor" diverges, see Thickness
    n_s: int = 4097
    db_fm: float = 0.025
    q_cut: float = 24.0
    shape_override: dict[str, float] | None = None
    inputs_json: Path = INPUTS_JSON
    horiuchi_json: Path = HORIUCHI_JSON
    allow_unpinned_inputs: bool = False
    pauli: bool = False  # M2': Bertulani & De Conti Pauli-blocked local-density phase
    pauli_variant: str | None = None  # REQUIRED with pauli: clementel-villi | bdc-pn
    pauli_emin_mev: float | None = None  # REQUIRED with pauli: nodes below keep the black disk
    pauli_scale: float = 1.0
    measured_sigma_tot: Path | None = None  # report-only comparison table (JSON), never embedded

    def variant(self) -> str:
        if self.pauli:
            return "M2prime-c" if self.beta == "horiuchi-table" else "M2prime"
        if self.beta == "horiuchi-table":
            return "V-c" if self.pp_level == "none" else "V-c+anchored"
        return "V-a" if self.pp_level != "none" else "V-b"

    def __post_init__(self) -> None:
        if self.pauli and self.pauli_variant not in PAULI_VARIANTS:
            raise M2Error(f"pauli needs an explicit pauli_variant in {PAULI_VARIANTS}")
        if self.pauli and self.pauli_emin_mev is None:
            raise M2Error("pauli needs an explicit pauli_emin_mev (no default)")
        if self.pp_level == "unchanged":
            object.__setattr__(self, "pp_level", "none")
        if self.beta not in BETA_MODES or self.pp_level not in PP_LEVEL_MODES:
            raise M2Error("unknown --beta / --pp-level mode")
        if self.pp_level != "none" and not self.dsg_files:
            raise M2Error("--pp-level said-anchored-ge150 needs the SAID DSG scans (dsg_files)")


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def classify_said_dir(input_dir: Path) -> dict[str, list[Path]]:
    """SGT/ALFA prediction files of ``input_dir`` grouped by their content (file names are free)."""
    groups: dict[str, list[Path]] = {}
    for p in sorted(Path(input_dir).iterdir()):
        if not p.is_file():
            continue
        try:
            _, ch, ob, _, _ = parse_said_prediction(p.read_bytes().decode("ascii", "replace"))
        except M2Error:
            continue
        groups.setdefault(f"{ch.lower()}_{ob.lower()}", []).append(p)
    return groups


class M2Model:
    """The M2 shape for a list of targets ``(name, z, a)`` (see the module docstring)."""

    def __init__(self, cfg: M2Config, targets: list[tuple[str, int, int]]) -> None:
        self.cfg = cfg
        self.inputs = load_inputs(cfg.inputs_json)
        self.r_p = float(self.inputs["r_proton_fm"]["value"])
        hor = HoriuchiBeta.load(cfg.horiuchi_json) if cfg.beta == "horiuchi-table" else None
        self.nn = NNInputs.from_files(classify_said_dir(cfg.input_dir), hor)
        self.pinned_said = {v for v in self.inputs["said"].get("file_sha256", {}).values() if v}
        self.pinned_dsg = {v for v in self.inputs["pp_dsg"].get("file_sha256", {}).values() if v}
        self.unpinned = not self.pinned_said
        if self.unpinned and not cfg.allow_unpinned_inputs:
            raise M2Error(
                "said.file_sha256 in the inputs JSON is empty: the SAID files are not pinned "
                "(use allow_unpinned_inputs / --allow-unpinned-inputs for an exploratory build)"
            )
        read = {h for s in self.nn.series.values() for h in s.sha256}
        if self.pinned_said and not read <= self.pinned_said:
            raise M2Error(f"SAID files with unpinned sha256: {sorted(read - self.pinned_said)}")
        self.targets = {n: (z, a) for n, z, a in targets}
        self.dens: dict[str, HODensity | FermiDensity] = {}
        self.entries: dict[str, dict[str, Any]] = {}
        self.thick: dict[str, Thickness | PauliThickness] = {}
        self.rms: dict[str, dict[str, Any]] = {}
        for name, _, _ in targets:
            d, ent = target_density(self.inputs, name, cfg.shape_override)
            self.dens[name], self.entries[name] = d, ent
            rec = rms_check(d, name, self.inputs)
            if ent.get("surrogate"):
                rec["pass"] = None
                rec["note"] = RMS_TOL_SURROGATE_NOTE
            elif not rec["pass"]:
                raise M2Error(f"{name}: rms radius check failed: {rec}")
            self.rms[name] = rec
            self.thick[name] = self._make_thickness(name, cfg.db_fm)

    def _make_thickness(
        self, name: str, db_fm: float, *, dq: float = 0.02, q_cut: float | None = None,
        dr: float = 0.01,
    ) -> Thickness | PauliThickness:  # fmt: skip
        d, cfg = self.dens[name], self.cfg
        qc = cfg.q_cut if q_cut is None else q_cut
        if cfg.pauli:
            if cfg.unfold == "formfactor" and isinstance(d, FermiDensity):
                raise M2Error("pauli (M2') needs a position-space point density: use rms-matched")
            z, a = self.targets[name]
            assert cfg.pauli_variant is not None
            return PauliThickness(d, b_grid(d, db_fm), z, a - z, variant=cfg.pauli_variant,
                                  r_p=self.r_p, q_cut=qc, dq=dq, dr=dr,
                                  pauli_scale=cfg.pauli_scale)  # fmt: skip
        return Thickness(d, b_grid(d, db_fm), self.r_p, qc, dq, cfg.unfold, dr)

    def check_dsg_pinned(self, shas: list[str]) -> bool:
        """Fail closed on a DSG file whose hash is not pinned; returns True if the DSG hashes
        are unpinned (empty JSON list) and that was allowed."""
        if not self.pinned_dsg:
            if not self.cfg.allow_unpinned_inputs:
                raise M2Error("pp_dsg.file_sha256 in the inputs JSON is empty (unpinned DSG)")
            self.unpinned = True
            return True
        if not set(shas) <= self.pinned_dsg:
            raise M2Error(f"DSG files with unpinned sha256: {sorted(set(shas) - self.pinned_dsg)}")
        return False

    def _chi(self, name: str, th: Thickness | PauliThickness, e_mev: float) -> CF64:
        z, a = self.targets[name]
        if isinstance(th, PauliThickness):
            th.set_energy(e_mev)
        chi = eikonal_chi(th, self.nn, e_mev, z, a - z)
        if th.tail > TAIL_TOL:
            raise M2Error(
                f"{name} at {e_mev:g} MeV: the form-factor unfolding integrand is {th.tail:.3g} at "
                f"q_cut (diverges for beta < r_p^2/3); use unfold='rms-matched' (exploratory)"
            )
        return chi

    def node(self, name: str, e_mev: float, pcm_mev: float, n_q: int) -> tuple[F64, dict[str, Any]]:
        th = self.thick[name]
        k = pcm_mev / HBARC
        b = th.b
        chi = self._chi(name, th, e_mev)
        edge_val = abs(1.0 - complex(np.exp(1j * chi[-1])))
        res = m2_shape(b, chi, k, self.cfg.n_s)
        edges, dev = quantiles_mu(res.s, res.density, n_q)
        g = 4.0 * res.s * res.density
        norm = float(np.sum(0.5 * (g[1:] + g[:-1]) * np.diff(res.s)))
        re_norm = float(np.sum(0.5 * ((g / norm)[1:] + (g / norm)[:-1]) * np.diff(res.s)))
        rec = {"e_mev": e_mev, "sigma_el_m2_mb": 10.0 * res.sigma_el_fm2,
               "sigma_tot_m2_mb": 10.0 * res.sigma_tot_ot_fm2,
               "sigma_react_m2_mb": 10.0 * res.sigma_react_fm2,
               "norm_deviation": max(dev, abs(re_norm - 1.0)),
               "min_density_fm2_sr": float(res.density.min()),
               "unitarity_residual": abs(res.sigma_tot_ot_fm2 - res.sigma_el_fm2
                                         - res.sigma_react_fm2) / res.sigma_tot_ot_fm2,
               "optical_theorem_residual": abs(res.sigma_tot_ot_fm2 - res.sigma_tot_b_fm2)
                                           / res.sigma_tot_ot_fm2,
               "q_first_min_fm_inv": first_minimum_q(res.s, res.density, k),
               "b_edge_residual": edge_val}  # fmt: skip
        if isinstance(th, PauliThickness):
            rec["pauli_factor_peak"] = th.pauli_peak
        return edges, rec

    def convergence(self, name: str, e_mev: float, pcm_mev: float) -> dict[str, Any]:
        """Step-halving record of the integration steps (plan 4.2 trigger 2, tolerance 1e-3), each
        component reported separately: ``b_step`` (impact-parameter step halved), ``s_grid`` (the
        angular grid of the amplitude doubled), ``radial_step`` (radial step of the form-factor /
        density transform halved), ``q_step`` (Hankel q step
        halved) and ``q_cutoff`` (cutoff x 1.5).
        The metric of each is the largest of the relative changes of sigma_el and sigma_tot and of
        the density at the common angles relative to its forward value. The radial and q
        components do not exist for the analytic HO thickness (``None``). ``max`` combines them."""
        k = pcm_mev / HBARC
        th1 = self.thick[name]
        n_s = self.cfg.n_s
        r1 = m2_shape(th1.b, self._chi(name, th1, e_mev), k, n_s)
        analytic = isinstance(th1, Thickness) and th1.point is not None  # HO: analytic fold

        def metric(th: Thickness | PauliThickness, n: int = n_s) -> float:
            r2 = m2_shape(th.b, self._chi(name, th, e_mev), k, n)
            d2 = r2.density[:: (n - 1) // (n_s - 1)]
            return max(abs(r2.sigma_el_fm2 / r1.sigma_el_fm2 - 1.0),
                       abs(r2.sigma_tot_ot_fm2 / r1.sigma_tot_ot_fm2 - 1.0),
                       float(np.max(np.abs(d2 - r1.density)) / r1.density[0]))  # fmt: skip

        db, qc = self.cfg.db_fm, self.cfg.q_cut
        comp: dict[str, float | None] = {
            "b_step": metric(self._make_thickness(name, db / 2)),
            "s_grid": metric(th1, 2 * (n_s - 1) + 1),
            "radial_step": None if analytic else metric(self._make_thickness(name, db, dr=0.005)),
            "q_step": None if analytic else metric(self._make_thickness(name, db, dq=0.01)),
            "q_cutoff": None if analytic
            else metric(self._make_thickness(name, db, q_cut=1.5 * qc)),
        }  # fmt: skip
        mx = max(v for v in comp.values() if v is not None)
        return {"components": comp, "max": mx, "tol": CONVERGENCE_TOL,
                "passed": mx <= CONVERGENCE_TOL, "e_mev": e_mev}  # fmt: skip

    def provenance(self) -> dict[str, Any]:
        inp = self.inputs
        dens = {}
        for n, ent in self.entries.items():
            d = self.dens[n]
            dens[n] = {"form": ent["form"], "surrogate": bool(ent.get("surrogate")),
                       "parameters": {k: v for k, v in vars(d).items()},
                       "sources": ent.get("sources"), "rms_check": self.rms[n]}  # fmt: skip
        return {"said": inp["said"] | {"read": self.nn.provenance()},
                "iaea": inp["iaea_rms_charge_radius_fm"], "r_proton": inp["r_proton_fm"],
                "densities": dens, "inputs_json_sha256": _sha256(self.cfg.inputs_json),
                "inputs_unpinned": self.unpinned,
                "horiuchi_json_sha256": _sha256(self.cfg.horiuchi_json)
                if self.cfg.beta == "horiuchi-table" else None}  # fmt: skip


def black_disk_limit_check(n_s: int = 2049) -> dict[str, float]:
    """Sharp-disk limit: ``chi = 50 i`` inside R, 0 outside gives ``F = i k R J1(qR)/q``; the max
    deviation of the normalised eikonal density from the analytic ``|2 J1(qR)/(qR)|^2`` shape
    (relative to the forward value) at ``R = 4 fm``, ``k = 2.5 fm^-1`` (``qR`` up to 20) is
    returned. Tolerance stated by the builder: 1e-3."""
    from ionmc.nuclear.elastic_build import disk_amplitude_sq

    r_d, k = 4.0, 2.5
    b = np.arange(0.0, 4.0 + 0.005 * 3, 0.005)
    b = b[b <= r_d + 1e-12]
    chi = np.full(b.size, 50.0j)
    s = np.linspace(0.0, 1.0, n_s)
    f = eikonal_amplitude(b, chi, k, s)
    d = np.abs(f) ** 2
    ref = disk_amplitude_sq(2.0 * k * s * r_d)
    return {"max_rel_dev": float(np.max(np.abs(d / d[0] - ref))), "tolerance": 1e-3}


def load_dsg_dir(dsg_dir: Path) -> dict[float, list[Path]]:
    """``{T_lab: [seg1, seg2, seg3]}`` from ``dsg_dir/dsg_manifest.json`` (``{"100": ["f1", "f2",
    "f3"], ...}``, file names relative to the directory): the SAID DSG outputs carry no T_lab."""
    d = Path(dsg_dir)
    man = json.loads((d / "dsg_manifest.json").read_text(encoding="utf-8"))
    return {float(t): [d / n for n in names] for t, names in man.items()}


def load_measured_sigma_tot(path: Path | None) -> dict[str, dict[float, float]]:
    """Caller-supplied measured total cross sections for the report-only comparison: JSON
    ``{"C-12": {"100": 492.5, ...}, "O-16": {...}}`` (T_lab [MeV] -> sigma_tot [mb]). ``None`` gives
    an empty table (comparison skipped). No measured value is embedded in the source."""
    if path is None:
        return {}
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    return {t: {float(e): float(v) for e, v in rows.items()} for t, rows in doc.items()}


def sigma_tot_comparison(
    measured: dict[str, dict[float, float]], nodes: dict[str, list[dict[str, Any]]]
) -> dict[str, Any]:  # fmt: skip
    """Report-only: model sigma_tot (linear in E between nodes) against the supplied table;
    energies outside the model nodes are reported as ``None``."""
    out: dict[str, Any] = {}
    for tgt, rows in measured.items():
        recs = nodes.get(tgt)
        if not recs:
            out[tgt] = "target not built"
            continue
        e = [r["e_mev"] for r in recs]
        sg = [r["sigma_tot_m2_mb"] for r in recs]
        out[tgt] = {}
        for en, meas in sorted(rows.items()):
            if e[0] <= en <= e[-1]:
                mod = float(np.interp(en, e, sg))
                out[tgt][f"{en:g}_mev"] = {"measured_mb": meas, "model_mb": mod,
                                           "model_over_measured": mod / meas}  # fmt: skip
            else:
                out[tgt][f"{en:g}_mev"] = None
    return out
