"""Deterministic and sampled producers of the V3-005C rows V10, V10-A, X-ENDF and D9-BGG (plan
Amendment 14 (d), (e); decision 0041 slice C).

Pure functions on the data layer and on the built elastic table (``ionmc.nuclear.*``); the steps of
``steps_v5d.py`` call them. This module never contains a measured value: data come from the
hash-verified cache (``IONMC_CACHE_DIR``) and are read at run time. The functions below are added in
stages; this first stage holds the shared loaders."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.data import cache, endf6, exfor
from ionmc.data.ame import load_ame2020
from ionmc.nuclear import bgg, law5
from ionmc.nuclear import events as ev
from ionmc.nuclear.build import TARGETS
from ionmc.nuclear.elastic_build import (
    XENDF_ENERGIES_MEV,
    XENDF_TARGETS,
    XENDF_THETA_DEG,
    disk_density_mb_sr,
)
from ionmc.nuclear.elastic_tables import ElasticTable


def h1_law5(cache_dir: str | Path | None = None) -> tuple[Any, float]:
    """LA150 H-1 MF6/MT2 LAW=5 section and the projectile AWI (as the elastic table builder)."""
    cdir = cache.resolve_cache_dir(cache_dir)
    zip_path = cache.verify("endf-b8.0-protons", cdir)
    h1 = endf6.parse_endf(endf6.read_member(zip_path, "ENDF-B-VIII.0_protons/p-001_H_001.endf"))
    return law5.parse_law5(h1.sections[(6, 2)]), law5.projectile_awi(h1)


def target_masses_mev(cache_dir: str | Path | None = None) -> dict[str, float]:
    """AME2020 nuclear mass ``m_t`` [MeV] of every p + A table target (the builder's ``m_t_mev``)."""
    cdir = cache.resolve_cache_dir(cache_dir)
    ame = load_ame2020(cache.verify("ame2020-mass", cdir).read_text(encoding="ascii"))
    return {t.name: float(ev.build_event_model(ame, t.z, t.a).m_t_mev) for t in TARGETS}


def sigma_pp_direct_mb(sec: Any, awi: float, e_mev: float, anchor_mev: float = 150.0,
                       n_mu: int = 4001) -> float:
    """Direct (builder-formula) evaluation of the transported p-p sigma [mb]: the LA150 nuclear-plus-
    interference integral over 0 <= mu_CM <= 0.96 up to the 150 MeV anchor, the anchor value times
    ``pp_bgg(E) / pp_bgg(150)`` above it."""
    def ni(e: float) -> float:
        mu = np.linspace(0.0, law5.MU_CUT_PP, n_mu)
        d = law5.ni_density_ltp1(sec, e * 1e6, mu, awi, 1, 1)
        return float(2.0 * np.pi * law5.simpson(d, mu) * 1.0e3)

    if e_mev <= anchor_mev:
        return ni(e_mev)
    return ni(anchor_mev) * bgg.pp_bgg_mb(e_mev) / bgg.pp_bgg_mb(anchor_mev)


# =============================================================================================
# kinematics, Rutherford and the model value
# =============================================================================================



MP_MEV = bgg.M_PROTON_MEV
ALPHA_HBARC_MEV_FM = bgg.FINE_STRUCTURE * bgg.HBARC_MEV_FM
FM2_TO_MB = 10.0
V10A_ANGLES_DEG = (10.0, 15.0, 20.0)
V10A_BRACKET_DEG = 3.0
V10A_MIN_CUT = 0.8
V10A_TOL_A, V10A_TOL_B, V10A_TOL_C = 0.30, 0.20, 0.12
V10A_COVERAGE = {"a_points": 12, "o16_low_mev": (130.0, 140.0), "o16_high_mev": (175.0, 205.0), "c12_datasets": 1,
                 "b_datasets": 3, "c_points": 4}  # fmt: skip


def two_body_frame(e_lab_mev: float, m_t_mev: float) -> dict[str, float]:
    """CM-frame quantities of p + target at proton kinetic energy ``e_lab_mev`` (relativistic)."""
    e_tot = e_lab_mev + MP_MEV
    p_lab = math.sqrt(e_lab_mev * e_lab_mev + 2.0 * e_lab_mev * MP_MEV)
    w = math.sqrt(MP_MEV**2 + m_t_mev**2 + 2.0 * e_tot * m_t_mev)
    p_cm = p_lab * m_t_mev / w
    e1 = math.sqrt(MP_MEV**2 + p_cm**2)
    e2 = math.sqrt(m_t_mev**2 + p_cm**2)
    return {"w": w, "p_cm": p_cm, "e1_cm": e1, "e2_cm": e2, "beta_cm": p_lab / (e_tot + m_t_mev),
            "gamma_cm": (e_tot + m_t_mev) / w, "beta_star": p_cm / e1}  # fmt: skip


def rutherford_cm_mb_sr(theta_cm_deg: NDArray[np.float64], e_lab_mev: float, m_t_mev: float, z_t: int,
                        z_p: int = 1) -> NDArray[np.float64]:
    """Unscreened relativistic point-Rutherford ``d sigma / d Omega_CM`` [mb/sr] of two spin-0 point charges:
    ``(Z1 Z2 alpha hbar c)^2 (E1* E2*)^2 / (4 W^2 (p* c)^4 sin^4(theta/2))`` (Born, equals the classical
    Rutherford formula in the non-relativistic limit; the factor 10 converts fm^2 to mb)."""
    f = two_body_frame(e_lab_mev, m_t_mev)
    s = np.sin(np.radians(np.asarray(theta_cm_deg, dtype=float)) / 2.0)
    coef = (z_p * z_t * ALPHA_HBARC_MEV_FM) ** 2 * (f["e1_cm"] * f["e2_cm"]) ** 2 / (4.0 * f["w"] ** 2 * f["p_cm"] ** 4)
    return np.asarray(coef * FM2_TO_MB / s**4)


def lab_to_cm(theta_lab_deg: NDArray[np.float64], e_lab_mev: float, m_t_mev: float
              ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Proton laboratory angle -> CM angle [deg] and the factor ``d Omega_L / d Omega_CM`` (so that
    ``d sigma/d Omega_CM = d sigma/d Omega_L x factor``) for elastic scattering off a heavier target
    (single-valued: the target is heavier than the proton, ``beta_cm / beta* < 1``)."""
    f = two_body_frame(e_lab_mev, m_t_mev)
    r = f["beta_cm"] / f["beta_star"]
    if r >= 1.0:
        raise ValueError("lab_to_cm needs a target heavier than the projectile (beta_cm < beta*)")

    def lab_of(tc: float) -> float:
        return math.atan2(math.sin(tc), f["gamma_cm"] * (math.cos(tc) + r))

    out, jac = [], []
    for tl in np.radians(np.asarray(theta_lab_deg, dtype=float)):
        lo, hi = 0.0, math.pi
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if lab_of(mid) < tl else (lo, mid)
        tc = 0.5 * (lo + hi)
        h = 1e-6
        dtl = (lab_of(tc + h) - lab_of(tc - h)) / (2 * h)
        out.append(math.degrees(tc))
        jac.append(math.sin(tl) * dtl / math.sin(tc))
    return np.array(out), np.array(jac)


class ModelTable:
    """The pinned elastic table with the model ``d sigma / d Omega_CM`` of the p + A shape."""

    def __init__(self, table_id: str, cache_dir: str | Path | None = None) -> None:
        self.tab = ElasticTable.load(cache_dir, table_id)
        self.masses = target_masses_mev(cache_dir)

    def index(self, name: str) -> int:
        return self.tab.target_names.index(name)

    def sigma_el_mb(self, name: str, e_mev: float) -> float:
        return float(self.tab.sigma_barn(self.index(name), e_mev) * 1.0e3)

    def nuclear_mb_sr(self, name: str, e_mev: float, theta_cm_deg: NDArray[np.float64]) -> NDArray[np.float64]:
        it = self.index(name)
        a = self.tab.arrays
        grid = np.asarray(a["grid_e_mev"], dtype=float)
        r_fm = float(np.interp(e_mev, grid, a["radius_fm"][it]))
        pcm = float(np.interp(e_mev, grid, a["p_cm_mev"][it]))
        return disk_density_mb_sr(self.sigma_el_mb(name, e_mev) * 1.0e-3, r_fm, pcm,
                                  np.asarray(theta_cm_deg, dtype=float))  # fmt: skip

    def total_mb_sr(self, name: str, e_mev: float, theta_cm_deg: NDArray[np.float64]
                    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Nuclear plus Rutherford (incoherent, the transport's content) and the Rutherford part."""
        z = int(self.tab.arrays["target_z"][self.index(name)])
        ru = rutherford_cm_mb_sr(theta_cm_deg, e_mev, self.masses[name], z)
        return self.nuclear_mb_sr(name, e_mev, theta_cm_deg) + ru, ru


# =============================================================================================
# V10-A: p + A nuclear elastic vs measured data (plan Amendment 14 (e), items 1-8)
# =============================================================================================


@dataclass(frozen=True)
class AngularSet:
    """One gating angular dataset: EXFOR entry, subentry, target, label and the frame of the data."""

    dataset_id: str  # registry id, e.g. exfor-c0063
    subentry: str
    target: str
    lab_frame: bool = False
    energies_from_column: bool = False  # c0055: every energy of its EN column in 120-200 MeV


GATING_SETS = (
    AngularSet("exfor-c0063", "C0063002", "O-16"),  # Kelly 1989, 135 MeV
    AngularSet("exfor-c0141", "C0141002", "O-16", lab_frame=True),  # Kelly 1990, 179.9 MeV
    AngularSet("exfor-c0148", "C0148002", "O-16"),  # Seifert 1993, 201.4 MeV
    AngularSet("exfor-c0550", "C0550002", "O-16"),  # Glover 1985, 200 MeV
    AngularSet("exfor-o0274", "O0274009", "C-12"),  # Strauch & Titus 1956, 96 MeV
    AngularSet("exfor-o0553", "O0553008", "C-12"),  # Rolland 1966, 75 MeV
    AngularSet("exfor-o0553", "O0553002", "C-12"),  # Rolland 1966, 150 MeV
    AngularSet("exfor-c0055", "C0055002", "C-12", energies_from_column=True),  # Meyer 1983
    AngularSet("exfor-c0057", "C0057002", "C-12"),  # Meyer 1981, 200 MeV
    AngularSet("exfor-c1420", "C1420002", "C-12"),  # Meyer 1988, 250 MeV
)
V10A_ENERGY_WINDOW_C0055 = (120.0, 200.0)
_ENTRY_CACHE: dict[tuple[str, str], Any] = {}


def parse_entry(dataset_id: str, cache_dir: str | Path | None = None) -> Any:
    cdir = cache.resolve_cache_dir(cache_dir)
    key = (str(cdir), dataset_id)
    if key not in _ENTRY_CACHE:
        _ENTRY_CACHE[key] = exfor.parse_entry(cache.verify(dataset_id, cdir).read_text(encoding="ascii"))
    return _ENTRY_CACHE[key]


def entry_energy_mev(entry: Any, subentry_id: str) -> float | None:
    """Incident energy [MeV] of a subentry: its own COMMON ``EN`` if present, else the COMMON ``EN`` of the
    entry's first subentry (EXFOR convention: a COMMON of subentry 001 applies to all subentries)."""
    subs = {s.subentry_id: s for s in entry.subentries}
    for s in (subs[subentry_id], entry.subentries[0]):
        if s.common is not None and "EN" in s.common.heads:
            i = s.common.heads.index("EN")
            v = s.common.values[i]
            if v is not None:
                return float(exfor.energy_total_mev(np.array([v]), s.common.units[i], 1)[0])
    return None


def read_angular(spec: AngularSet, cache_dir: str | Path | None = None) -> list[dict[str, Any]]:
    """Measured ``d sigma/d Omega`` of a gating dataset, one record per incident energy: ``energy_mev``,
    ``theta_lab``/``theta_cm`` [deg] (CM after the laboratory conversion of item 3), ``ds_cm_mb_sr`` and its
    relative error (a fraction; absent errors are NaN). Rows with non-finite or non-positive data are
    dropped and counted."""
    entry = parse_entry(spec.dataset_id, cache_dir)
    sub = next(s for s in entry.subentries if s.subentry_id == spec.subentry)
    d, h = sub.data, sub.data.heads
    ang_i = next(i for i, x in enumerate(h) if x.startswith("ANG") and not x.startswith("ANG-ERR"))
    dat_i = next(i for i, x in enumerate(h) if x.startswith("DATA") and not x.startswith("DATA-ERR"))
    err_i = next((i for i, x in enumerate(h) if x in ("DATA-ERR", "ERR-S", "ERR-T")), None)
    ang, ds = d.column(ang_i), d.column(dat_i)
    if d.units[dat_i] != "MB/SR":
        raise ValueError(f"{spec.subentry}: DATA unit {d.units[dat_i]!r} is not MB/SR")
    rel = np.full(ang.shape, np.nan)
    if err_i is not None:
        e = d.column(err_i)
        if d.units[err_i] == "PER-CENT":
            rel = e / 100.0
        elif d.units[err_i] == "MB/SR":
            rel = e / ds
    en_col = d.column(h.index("EN")) if "EN" in h else None
    if spec.energies_from_column:
        energies = sorted({float(x) for x in en_col if V10A_ENERGY_WINDOW_C0055[0] <= x <= V10A_ENERGY_WINDOW_C0055[1]})
    else:
        e0 = entry_energy_mev(entry, spec.subentry)
        if e0 is None:
            raise ValueError(f"{spec.subentry}: no incident energy in COMMON")
        energies = [e0]
    out = []
    for en in energies:
        sel = np.isfinite(ang) & np.isfinite(ds) & (ds > 0.0)
        if en_col is not None:
            sel &= en_col == en
        a, y, r = ang[sel], ds[sel], rel[sel]
        order = np.argsort(a, kind="stable")
        a, y, r = a[order], y[order], r[order]
        if spec.lab_frame:
            tc, jac = lab_to_cm(a, en, target_masses_mev()[spec.target])
            y, a = y * jac, tc
            order = np.argsort(a, kind="stable")
            a, y, r = a[order], y[order], r[order]
        out.append({"dataset": f"{spec.subentry}@{en:g}", "entry": spec.dataset_id, "subentry": spec.subentry,
                    "target": spec.target, "energy_mev": en, "theta_cm": a, "ds_cm_mb_sr": y, "rel_err": r,
                    "n_points": int(a.size), "lab_frame_converted": spec.lab_frame})  # fmt: skip
    return out


def first_minimum_deg(theta: NDArray[np.float64], y: NDArray[np.float64]) -> float | None:
    """theta_min,1: the first measured point above 10 deg that is lower than both of its neighbours."""
    for i in range(1, theta.size - 1):
        if theta[i] > 10.0 and y[i] < y[i - 1] and y[i] < y[i + 1]:
            return float(theta[i])
    return None


def loglin(theta: NDArray[np.float64], y: NDArray[np.float64], x: float,
           bracket_deg: float = V10A_BRACKET_DEG) -> float | None:
    """Log-linear interpolation in theta of the two data points that bracket ``x``, both within
    ``bracket_deg`` of it; ``None`` otherwise (the point is not evaluated; item 5)."""
    for i in range(theta.size - 1):
        if theta[i] <= x <= theta[i + 1] and theta[i + 1] > theta[i]:
            if x - theta[i] <= bracket_deg and theta[i + 1] - x <= bracket_deg:
                w = (x - theta[i]) / (theta[i + 1] - theta[i])
                return float(math.exp((1 - w) * math.log(y[i]) + w * math.log(y[i + 1])))
            return None
    return None


def model_first_minimum_deg(model: ModelTable, name: str, e_mev: float) -> float | None:
    th = np.linspace(10.0, 170.0, 1601)
    m, _ = model.total_mb_sr(name, e_mev, th)
    for i in range(1, th.size - 1):
        if m[i] < m[i - 1] and m[i] < m[i + 1]:
            return float(th[i])
    return None


def trapezoid_2pi_sin(theta_deg: NDArray[np.float64], y: NDArray[np.float64]) -> float:
    t = np.radians(theta_deg)
    f = 2.0 * math.pi * np.sin(t) * y
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(t)))


def evaluate_angular(model: ModelTable, rec: dict[str, Any]) -> dict[str, Any]:
    """(a) and (b) of one dataset record (items 4-6)."""
    th, y, name, e = rec["theta_cm"], rec["ds_cm_mb_sr"], rec["target"], rec["energy_mev"]
    tmin = first_minimum_deg(th, y)
    cut_src = "measured"
    if tmin is None:
        tmin_model = model_first_minimum_deg(model, name, e)
        cut_src, cut = "model", (V10A_MIN_CUT * tmin_model if tmin_model else math.inf)
    else:
        cut = V10A_MIN_CUT * tmin
    points = []
    for ang in V10A_ANGLES_DEG:
        row: dict[str, Any] = {"theta_deg": ang}
        if ang >= cut:
            row.update(evaluated=False, reason=f"theta >= {V10A_MIN_CUT} theta_min,1 ({cut_src}: {cut:.2f} deg)")
        else:
            d = loglin(th, y, ang)
            if d is None:
                row.update(evaluated=False, reason="no bracketing data pair within 3 deg")
            else:
                tot, ru = model.total_mb_sr(name, e, np.array([ang]))
                row.update(evaluated=True, data_mb_sr=d, model_mb_sr=float(tot[0]), rutherford_mb_sr=float(ru[0]),
                           rutherford_share=float(ru[0] / tot[0]), ratio_minus_1=float(tot[0] / d - 1.0),
                           pass_=bool(abs(tot[0] / d - 1.0) <= V10A_TOL_A))  # fmt: skip
        points.append(row)
    b: dict[str, Any] = {"evaluated": False}
    if tmin is None:
        b["reason"] = "no measured minimum: (b) is not evaluated (item 5)"
    else:
        d10 = loglin(th, y, 10.0)
        if d10 is None:
            b["reason"] = "no bracketing data pair around 10 deg"
        else:
            sel = (th > 10.0) & (th <= tmin)
            nodes = np.concatenate(([10.0], th[sel]))
            ydat = np.concatenate(([d10], y[sel]))
            tot, _ = model.total_mb_sr(name, e, nodes)
            i_d, i_m = trapezoid_2pi_sin(nodes, ydat), trapezoid_2pi_sin(nodes, tot)
            b = {"evaluated": True, "nodes": int(nodes.size), "theta_min_deg": tmin, "data_mb": i_d, "model_mb": i_m,
                 "ratio_minus_1": i_m / i_d - 1.0, "pass_": bool(abs(i_m / i_d - 1.0) <= V10A_TOL_B)}  # fmt: skip
    return {"dataset": rec["dataset"], "target": name, "energy_mev": e, "n_points": rec["n_points"],
            "theta_min_deg": tmin, "cut_deg": cut if math.isfinite(cut) else None, "cut_source": cut_src,
            "a": points, "b": b}  # fmt: skip


def _xs_column(entry: Any, subentry_id: str, unit_to_mb: bool = True) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    sub = next(s for s in entry.subentries if s.subentry_id == subentry_id)
    h = sub.data.heads
    en = exfor.energy_total_mev(sub.data.column(h.index("EN")), sub.data.units[h.index("EN")], 1)
    xs = exfor.xs_to_mb(sub.data.column(h.index("DATA")), sub.data.units[h.index("DATA")])
    ei = next((i for i, x in enumerate(h) if x in ("ERR-S", "ERR-T", "DATA-ERR")), None)
    er = exfor.xs_to_mb(sub.data.column(ei), sub.data.units[ei]) if ei is not None else 0.05 * xs
    er = np.where(np.isfinite(er) & (er > 0), er, 0.05 * xs)
    ok = np.isfinite(en) & np.isfinite(xs) & (xs > 0)
    o = np.argsort(en[ok], kind="stable")
    return en[ok][o], xs[ok][o], er[ok][o]


def neutron_total_at(entry: Any, subentry_id: str, e_mev: float) -> tuple[float, float] | None:
    """sigma_tot(n+A) [mb] linearly interpolated in E between the two neighbouring data points, with its
    propagated error (item 2); ``None`` outside the measured range."""
    en, xs, er = _xs_column(entry, subentry_id)
    j = int(np.searchsorted(en, e_mev))
    if j <= 0 or j >= en.size:
        return None
    w = (e_mev - en[j - 1]) / (en[j] - en[j - 1])
    return float((1 - w) * xs[j - 1] + w * xs[j]), float(math.hypot((1 - w) * er[j - 1], w * er[j]))


def level_points(model: ModelTable, cache_dir: str | Path | None = None) -> dict[str, Any]:
    """(c): sigma_el of the table vs sigma_tot(n+C) - sigma_nonel(p+C) at every Auce 2005 D0356-002 energy and at
    Ingemarsson O0579-003 (65.5 MeV); the oxygen proxy (O0579-004 minus Finlay 13569-005) is report-only."""
    ntot_c = parse_entry("exfor-13753", cache_dir)
    gating, report = [], []
    for ds_id, sub, target, tot_entry, tot_sub, gate in (
        ("exfor-d0356", "D0356002", "C-12", "exfor-13753", "13753007", True),
        ("exfor-o0579", "O0579003", "C-12", "exfor-13753", "13753007", True),
        ("exfor-o0579", "O0579004", "O-16", "exfor-13569", "13569005", False),
    ):
        pn = parse_entry(ds_id, cache_dir)
        en, xs, er = _xs_column(pn, sub)
        tot = parse_entry(tot_entry, cache_dir) if tot_entry != "exfor-13753" else ntot_c
        for e, s_non, e_non in zip(en, xs, er, strict=True):
            row: dict[str, Any] = {"set": ds_id, "subentry": sub, "target": target, "energy_mev": float(e)}
            nt = neutron_total_at(tot, tot_sub, float(e))
            if nt is None or not 1.0 <= e <= 250.0:
                row.update(evaluated=False, reason="outside the neutron-total data range or the table domain")
            else:
                proxy, perr = nt[0] - s_non, math.hypot(nt[1], e_non)
                mod = model.sigma_el_mb(target, float(e))
                row.update(evaluated=True, sigma_tot_n_mb=nt[0], sigma_nonel_p_mb=float(s_non), proxy_mb=proxy,
                           proxy_err_mb=perr, model_sigma_el_mb=mod, ratio_minus_1=mod / proxy - 1.0,
                           pass_=bool(abs(mod / proxy - 1.0) <= V10A_TOL_C))  # fmt: skip
            (gating if gate else report).append(row)
    return {"gating": gating, "report_only": report}


def v10_a(model: ModelTable, cache_dir: str | Path | None = None) -> dict[str, Any]:
    """Row V10-A from the committed code and the cached data (deterministic). ``status`` is ``pass``, ``fail`` or
    ``not executed (evidence missing)`` (item 7: never a pass without the minimum coverage)."""
    datasets = []
    for spec in GATING_SETS:
        for rec in read_angular(spec, cache_dir):
            datasets.append(evaluate_angular(model, rec))
    lp = level_points(model, cache_dir)
    a_eval = [(d, p) for d in datasets for p in d["a"] if p["evaluated"]]
    b_eval = [d for d in datasets if d["b"]["evaluated"]]
    c_eval = [r for r in lp["gating"] if r["evaluated"]]
    lo, hi = V10A_COVERAGE["o16_low_mev"], V10A_COVERAGE["o16_high_mev"]
    cov = {
        "a_points": len(a_eval), "a_points_required": V10A_COVERAGE["a_points"],
        "o16_dataset_130_140": any(d["target"] == "O-16" and lo[0] <= d["energy_mev"] <= lo[1] and any(p["evaluated"] for p in d["a"]) for d in datasets),
        "o16_dataset_175_205": any(d["target"] == "O-16" and hi[0] <= d["energy_mev"] <= hi[1] and any(p["evaluated"] for p in d["a"]) for d in datasets),
        "c12_datasets": sum(1 for d in datasets if d["target"] == "C-12" and any(p["evaluated"] for p in d["a"])),
        "b_datasets": len(b_eval), "b_datasets_required": V10A_COVERAGE["b_datasets"],
        "c_points": len(c_eval), "c_points_required": V10A_COVERAGE["c_points"],
    }  # fmt: skip
    met = bool(cov["a_points"] >= V10A_COVERAGE["a_points"] and cov["o16_dataset_130_140"] and cov["o16_dataset_175_205"]
               and cov["c12_datasets"] >= V10A_COVERAGE["c12_datasets"] and cov["b_datasets"] >= V10A_COVERAGE["b_datasets"]
               and cov["c_points"] >= V10A_COVERAGE["c_points"])  # fmt: skip
    a_ok = all(p["pass_"] for _, p in a_eval)
    b_ok = all(d["b"]["pass_"] for d in b_eval)
    c_ok = all(r["pass_"] for r in c_eval)
    failures = {
        "a": [{"dataset": d["dataset"], "target": d["target"], "energy_mev": d["energy_mev"], "theta_deg": p["theta_deg"],
               "data_mb_sr": p["data_mb_sr"], "model_mb_sr": p["model_mb_sr"], "model_over_data_minus_1": p["ratio_minus_1"],
               "rutherford_share": p["rutherford_share"]} for d, p in a_eval if not p["pass_"]],
        "b": [{"dataset": d["dataset"], "target": d["target"], "energy_mev": d["energy_mev"], "theta_min_deg": d["b"]["theta_min_deg"],
               "data_mb": d["b"]["data_mb"], "model_mb": d["b"]["model_mb"], "model_over_data_minus_1": d["b"]["ratio_minus_1"]}
              for d in b_eval if not d["b"]["pass_"]],
        "c": [{k: r[k] for k in ("set", "subentry", "target", "energy_mev", "sigma_tot_n_mb", "sigma_nonel_p_mb", "proxy_mb",
                                 "proxy_err_mb", "model_sigma_el_mb", "ratio_minus_1")} for r in c_eval if not r["pass_"]],
    }  # fmt: skip
    po_fail = bool(any(f["target"] == "O-16" for f in (*failures["a"], *failures["b"])))
    status = "not executed (evidence missing)" if not met else ("pass" if a_ok and b_ok and c_ok else "fail")
    return {"status": status, "failures": failures, "a_or_b_fails_for_pO": po_fail, "coverage": cov, "coverage_met": met, "a_pass": a_ok, "b_pass": b_ok, "c_pass": c_ok,
            "tolerances": {"a": V10A_TOL_A, "b": V10A_TOL_B, "c": V10A_TOL_C}, "datasets": datasets,
            "level_points": lp, "pass": bool(met and a_ok and b_ok and c_ok)}  # fmt: skip


# =============================================================================================
# sampled cross-checks (r18, r19): histograms of table samples against the model and the data
# =============================================================================================
def sampled_dsigma_mb_sr(model: ModelTable, name: str, e_mev: float, edges_deg: NDArray[np.float64], seed: int,
                         n: int) -> dict[str, Any]:
    """``n`` mu_CM samples of the table at (``name``, ``e_mev``) histogrammed in the angle bins ``edges_deg``:
    the nuclear ``d sigma/d Omega`` [mb/sr] of the sampled histogram, its statistical error, the analytic model
    value of the same bins (bin average) and the z of the counts."""
    t = model.index(name)
    rng = np.random.default_rng(seed)
    mu = model.tab.sample_mu(t, e_mev, rng.random(n))
    th = np.degrees(np.arccos(np.clip(mu, -1.0, 1.0)))
    cnt = np.histogram(th, bins=edges_deg)[0].astype(float)
    sig = model.sigma_el_mb(name, e_mev)
    d_cos = np.cos(np.radians(edges_deg[:-1])) - np.cos(np.radians(edges_deg[1:]))
    samp = cnt / n * sig / (2.0 * math.pi * d_cos)
    err = np.sqrt(np.maximum(cnt, 1.0)) / n * sig / (2.0 * math.pi * d_cos)
    ana = np.empty_like(samp)
    for i in range(samp.size):
        fine = np.linspace(edges_deg[i], edges_deg[i + 1], 41)
        dens = model.nuclear_mb_sr(name, e_mev, fine)
        w = np.sin(np.radians(fine))
        ana[i] = np.trapezoid(dens * w, np.radians(fine)) / np.trapezoid(w, np.radians(fine))
    p = ana * 2.0 * math.pi * d_cos / sig
    z = (cnt - n * p) / np.sqrt(np.maximum(n * p * (1 - p), 1e-300))
    return {"seed": seed, "events": n, "counts": cnt.tolist(), "sampled_mb_sr": samp.tolist(),
            "sampled_err_mb_sr": err.tolist(), "model_mb_sr": ana.tolist(), "z_counts": z.tolist(),
            "max_abs_z": float(np.max(np.abs(z)))}  # fmt: skip


def v10a_xcheck(model: ModelTable, seeds: list[int], n: int, cache_dir: str | Path | None = None) -> dict[str, Any]:
    """Report-only sampled cross-check of V10-A: per gating dataset energy, ``n`` sampled events histogrammed in bins
    centred on the measured angles (data angles up to 60 deg), compared with the data and the model."""
    out, k = [], 0
    for spec in GATING_SETS:
        for rec in read_angular(spec, cache_dir):
            th = rec["theta_cm"][rec["theta_cm"] <= 60.0]
            if th.size < 3:
                continue
            mid = 0.5 * (th[1:] + th[:-1])
            edges = np.concatenate(([max(th[0] - 0.5 * (th[1] - th[0]), 0.1)], mid, [th[-1] + 0.5 * (th[-1] - th[-2])]))
            res = sampled_dsigma_mb_sr(model, rec["target"], rec["energy_mev"], edges, seeds[k], n)
            y = rec["ds_cm_mb_sr"][: th.size]
            tot, _ = model.total_mb_sr(rec["target"], rec["energy_mev"], th)
            res.update(dataset=rec["dataset"], target=rec["target"], energy_mev=rec["energy_mev"],
                       theta_deg=th.tolist(), data_mb_sr=y.tolist(),
                       sampled_over_data=(np.asarray(res["sampled_mb_sr"]) / y).tolist(),
                       model_total_over_data=(tot / y).tolist())  # fmt: skip
            out.append(res)
            k += 1
    return {"datasets": out, "events_per_dataset": n, "report_only": True,
            "note": "the sampled histogram is the nuclear part only (Coulomb is the MCS)"}  # fmt: skip


# =============================================================================================
# X-ENDF (plan Amendment 14 (d): model dsigma/dOmega_CM at 20-40 deg vs LA150 NI)
# =============================================================================================
def x_endf(model: ModelTable, cache_dir: str | Path | None = None) -> dict[str, Any]:
    """Model nuclear ``d sigma/d Omega_CM`` at 20-40 deg against the LA150 nuclear-plus-interference
    density for C-12, N-14, Ca-40 at 50 / 100 / 150 MeV (LTP=12 ratio interpolation, ``ltp12_ni_density_ratio``),
    recomputed here and compared with the arrays the table builder stored. Report only."""
    from ionmc.nuclear.elastic_build import ltp12_ni_density_ratio

    cdir = cache.resolve_cache_dir(cache_dir)
    zip_path = cache.verify("endf-b8.0-protons", cdir)
    th = np.asarray(XENDF_THETA_DEG)
    mu = np.cos(np.radians(th))
    stems = {t.name: (t.member, t.z) for t in TARGETS}
    rows: dict[str, Any] = {}
    for name in XENDF_TARGETS:
        stem, z = stems[name]
        mat = endf6.parse_endf(endf6.read_member(zip_path, f"ENDF-B-VIII.0_protons/{stem}.endf"))
        sec, awi = law5.parse_law5(mat.sections[(6, 2)]), law5.projectile_awi(mat)
        per_e: dict[str, Any] = {}
        for ie, e in enumerate(XENDF_ENERGIES_MEV):
            idx = [r.energy_ev for r in sec.records].index(e * 1.0e6)
            sig = float(mat.cross_section(2).interpolate(e * 1.0e6)[0])
            la = np.asarray(ltp12_ni_density_ratio(sec, idx, sig, mu, awi, z)) * 1.0e3  # b/sr -> mb/sr
            mod = model.nuclear_mb_sr(name, e, th)
            it = XENDF_TARGETS.index(name)
            stored_la = np.asarray(model.tab.arrays["xendf_la150_ni_mb_sr"][it, ie])
            stored_mod = np.asarray(model.tab.arrays["xendf_model_mb_sr"][it, ie])
            per_e[f"{e:g}"] = {
                "theta_deg": th.tolist(), "la150_ni_mb_sr": la.tolist(), "model_mb_sr": mod.tolist(),
                "ratio_model_over_la150": (mod / la).tolist(),
                "builder_arrays_agree": {"la150_max_rel": float(np.max(np.abs(la / stored_la - 1.0))),
                                         "model_max_rel": float(np.max(np.abs(mod / stored_mod - 1.0)))},
                "note": "not like-for-like: CNI is present in ENDF and absent in the model"}  # fmt: skip
        rows[name] = per_e
    return {"targets": rows, "report_only": True,
            "o16": "excluded: its MT2 is a copy of C-12 (decision 0041 data finding)"}  # fmt: skip


def x_endf_xcheck(model: ModelTable, seeds: list[int], n: int) -> dict[str, Any]:
    """r19 sampled cross-check: ``n`` table samples per (target, energy) histogrammed in the X-ENDF angle bins."""
    th = np.asarray(XENDF_THETA_DEG)
    edges = np.concatenate(([th[0] - 1.25], 0.5 * (th[1:] + th[:-1]), [th[-1] + 1.25]))
    out, k = {}, 0
    for name in XENDF_TARGETS:
        for e in XENDF_ENERGIES_MEV:
            out[f"{name}@{e:g}"] = sampled_dsigma_mb_sr(model, name, e, edges, seeds[k], n)
            k += 1
    return {"cases": out, "events_per_case": n, "angle_edges_deg": edges.tolist(), "report_only": True}


# =============================================================================================
# V10 (p-p level above Hale's range, report only) and D9-BGG
# =============================================================================================
V10_P_RANGE_GEV = (0.556, 0.729)
V10_FLAG_TOL = 0.08


def parse_pdg_pp(text: str) -> list[dict[str, Any]]:
    """Rows of the PDG rpp2022 ``pp_elastic.dat`` (fixed width text): ``p_lab`` [GeV/c], sigma and the errors [mb]
    and the reference text. Lines that do not start with ``<int> <float>`` are header text and are skipped."""
    rows = []
    for line in text.splitlines():
        tok = line.split()
        if len(tok) < 10 or not tok[0].isdigit():
            continue
        try:
            nums = [float(x) for x in tok[1:9]]
        except ValueError:
            continue
        rows.append({"index": int(tok[0]), "p_lab_gev": nums[0], "p_min": nums[1], "p_max": nums[2],
                     "sigma_mb": nums[3], "stat_plus": nums[4], "stat_minus": nums[5], "sys_plus_pct": nums[6],
                     "sys_minus_pct": nums[7], "reference": " ".join(tok[9:])})  # fmt: skip
    return rows


def v10_pp_level(model: ModelTable, cache_dir: str | Path | None = None) -> dict[str, Any]:
    """Row V10 (report only): the transported sigma_pp(E) at every PDG ``pp_elastic`` point with
    0.556 < p_lab <= 0.729 GeV/c; per-point ratios and the |model/data - 1| <= 0.08 flag, no verdict."""
    cdir = cache.resolve_cache_dir(cache_dir)
    rows = parse_pdg_pp(cache.verify("pdg-rpp2022-pp-elastic", cdir).read_text(encoding="ascii", errors="replace"))
    pts = []
    for r in rows:
        p = r["p_lab_gev"]
        if not V10_P_RANGE_GEV[0] < p <= V10_P_RANGE_GEV[1]:
            continue
        t = (math.sqrt(p * p * 1.0e6 + MP_MEV**2) - MP_MEV)
        mod = model.sigma_el_mb("H-1", t)
        err = math.hypot(max(r["stat_plus"], r["stat_minus"]), r["sigma_mb"] * max(r["sys_plus_pct"], r["sys_minus_pct"]) / 100.0)
        ratio = mod / r["sigma_mb"]
        pts.append({**r, "t_mev": t, "model_mb": mod, "data_err_mb": err, "ratio": ratio,
                    "within_8_percent": bool(abs(ratio - 1.0) <= V10_FLAG_TOL)})  # fmt: skip
    return {"points": pts, "n_points": len(pts), "p_lab_range_gev": list(V10_P_RANGE_GEV),
            "t_range_mev": [pts[0]["t_mev"], pts[-1]["t_mev"]] if pts else None,
            "n_within_8_percent": sum(1 for p in pts if p["within_8_percent"]), "flag_tolerance": V10_FLAG_TOL,
            "verdict": "none (report only)",
            "v10_pp_b": "SAID SP07 is unreachable from this environment: dropped and recorded (plan V10-pp-b)",
            "v10_pp_c": "report-only: inventory of the EXFOR p-p entries (shape comparison not computed here)",
            "evidence_gap": None if pts else (
                f"the acquired PDG pp_elastic file has {len(rows)} points with p_lab >= "
                f"{min(r['p_lab_gev'] for r in rows):g} GeV/c and none in the plan window "
                f"{V10_P_RANGE_GEV[0]} < p_lab <= {V10_P_RANGE_GEV[1]} GeV/c (T = 152.4-249.9 MeV)"),
            "disclosure": "the 160 MeV point was seen in research; the row was not blind at registration"}  # fmt: skip


def v10_pp_c_inventory(cache_dir: str | Path | None = None) -> dict[str, Any]:
    out = {}
    for ds in ("exfor-o1226", "exfor-o0145", "exfor-c2637", "exfor-c2606", "exfor-o2057"):
        try:
            ent = parse_entry(ds, cache_dir)
        except exfor.ExforError as exc:  # the P1 parser does not read every number format of these entries
            out[ds] = {"parse_error": str(exc)}
            continue
        subs = []
        for s in ent.subentries:
            if s.data is None or not s.reactions:
                continue
            subs.append({"subentry": s.subentry_id, "reaction": s.reactions[0].text[:70], "heads": list(s.data.heads),
                         "n_rows": len(s.data.rows)})  # fmt: skip
        out[ds] = subs
    return out


def d9_bgg(model: ModelTable, cache_dir: str | Path | None = None) -> dict[str, Any]:
    """Row D9-BGG (report only): the acquired Barashenkov sigma_inel(p+A) and sigma_el(p+A) arrays against the
    post-1997 EXFOR set in the V1b / Amendment 7 (i) windows (``nuclear_d9.WINDOWS``). sigma_inel: D0356 (evaluation
    role) with the pre-1997 and geant-val sets report-only; sigma_el: the level proxy of V10-A (c)
    (sigma_tot(n+A) - sigma_nonel(p+A)) at the post-1997 energies."""
    import nuclear_d9 as d9

    cdir = cache.resolve_cache_dir(cache_dir)
    zmap = {"C-12": 6, "N-14": 7, "O-16": 8, "Ca-40": 20}

    def inel(t: str) -> Any:
        return lambda e: np.array([bgg.barashenkov_mb(float(x), zmap[t])[1] for x in e])

    data = d9.load_data(cdir)
    sets: dict[str, Any] = {}
    for sid, role in d9.SETS.items():
        per_t = {}
        for t in zmap:
            pts = data.get(sid, {}).get(t)
            per_t[t] = {"data": "none"} if pts is None or pts.size == 0 else {
                "data": "present", "n_points": int(pts.shape[0]), "windows": d9.evaluate_model(pts, inel(t))}
        sets[sid] = {"role": role, "targets": per_t}
    lp = level_points(model, cdir)
    proxy: dict[str, Any] = {}
    for t in ("C-12", "O-16"):
        rows = [r for r in (*lp["gating"], *lp["report_only"]) if r["target"] == t and r["evaluated"]]
        if not rows:
            proxy[t] = {"data": "none"}
            continue
        pts = np.array([[r["energy_mev"], r["proxy_mb"], abs(r["proxy_err_mb"])] for r in rows])
        proxy[t] = {"data": "present", "n_points": len(rows),
                    "windows": d9.evaluate_model(pts, lambda e, t=t: np.array(
                        [bgg.barashenkov_elastic_mb(float(x), zmap[t]) for x in e]))}  # fmt: skip
    return {"sigma_inel": {"evaluation_set": d9.EVALUATION_SET, "sets": sets},
            "sigma_el_level_proxy": {"model": "Barashenkov sigma_tot(n+A) - sigma_inel(p+A) (>= 14 MeV)",
                                     "targets": proxy}, "windows": [list(w) for w in d9.WINDOWS],
            "report_only": True, "lineage": "BGG is the Geant4 lineage of the TOPAS reference (shared lineage)"}  # fmt: skip
