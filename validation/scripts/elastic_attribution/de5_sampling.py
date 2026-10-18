"""D-E5, table side (exploratory, report-only; Amendment 19 (d)): the pinned elastic table sampled in the D-E4 frame.

Cases 0-9 = rows 14-23 of the plan (section 3.2): water (p+O-16 and p+H-1) at 30 / 50 / 100 / 150 / 200 MeV (cases 0-4)
and graphite (p+C-12, 2.0 g/cm3) at the same energies (cases 5-9). Per case ``N_EVENTS`` = 1e6 elastic events are sampled
at the **fixed beam energy** (thin target: no slowing down, no multiple scattering, no EM), seed
``SeedSequence(20991000 + case)`` (Amendment 19 (f)). The target of every event is chosen with probability proportional to
``n_t sigma_t(E)`` of the case material; ``mu_CM`` comes from the P7 sampler (``ElasticTable.sample_mu``, the reference
inverse-CDF sampler of the P7 row, ``steps_v5d.p7_case``); the lab frame is the exact two-body kinematics of
``ionmc.nuclear.elastic_kin.two_body``. p-p: the table density already stops at ``theta_CM = 16.26 deg`` (``|mu| <= 0.96``,
asserted); the Geant4 side applies the same cut; both protons are scored (the faster one is the primary, the slower the
secondary). Kinetic-energy filter of the D-E4 scorer: records with ``E > E0 - 1e-4 MeV`` (undeflected primaries) are dropped
(plan 3.2); the fraction removed from the p+A classes is reported.

Signatures (plan section 2.2, H2 / H3 rows and the descriptive characterisation), per class (p+O, p-p, p+C):
- rate per primary ``n sigma L`` (L = 1 cm at 100-200 MeV, 2 mm at 30 / 50 MeV; N = 4e6 / 2e7 primaries) and its expected
  count; p-p rate = the rate above 16.26 deg CM;
- ``P(theta_lab > theta_w)`` for w = 2, 3 mm with ``R(E)(1 - cos theta_w) = w``, R = the Bethe CSDA range of water at the
  beam energy (the water-geometry thresholds, used for graphite too); ``theta_med`` (lab median of the primary); the tail moment
  ``T = <(1 - cos theta_lab) 1[theta_lab > 3 theta_med]>`` (all of the primary);
- p-p: ``<T_slow>`` and the exceedance of the fast proton (the same ``P`` entries);
- descriptive (report-only): the first local minimum above 3 deg of the 1-deg CM binned density (the frozen estimator; it also
  picks noise dips, so a variant needing a 3-sigma Poisson dip is reported next to it; the q range uses the frozen minimum), the large-q log-slope
  d ln(d sigma/d q^2)/d ln q above that minimum (12 q bins to theta_CM = 90 deg, bins with >= 30 counts), and a least-squares
  3-D Helm (R_eff, s) fit to the sampled density (free amplitude) - the Helm fit is a fit to the table, not to data;
- classification of the D-E4 scheme (:func:`classify`, applied to the table events and later to the Geant4 phase space):
  relative residual of the primary's energy to the elastic kinematic curves of p+H and p+A; the smaller residual wins; a
  *conflict* track (both below 1e-3, i.e. theta_lab below the crossover of about 1.8 deg) is kept: p-p iff the event has a secondary
  proton, else p+A. The raw conflict fraction is a diagnostic; the per-class rate after the rule is reported. Tracks off both curves
  are counted and **excluded** before the signatures (:func:`select_tracks`, shared with the Geant4 parser); a p-p event is accepted when its
  fast proton lies on the p+H curve (addendum (c)2; the table side is p-p single scattering, so every table p-p event qualifies); the Geant4 ratio hook
  compares ``rate_per_primary_after_rule`` (KE filter, exclusion and classification applied on both sides).
- the first-minimum, q-slope and Helm descriptors are recorded, not interpreted.

Geant4 ratios are made later by :func:`signature_ratios` from a Geant4 per-class summary of the same keys (D-E4); no Geant4
number is used here. Run: ``python -m validation.scripts.elastic_attribution.de5_sampling`` (about 1 minute).
"""

from __future__ import annotations

import math
import time
from pathlib import Path
from typing import Any

import numpy as np

from . import common as c

ENERGIES = (30, 50, 100, 150, 200)
N_EVENTS = 1_000_000
SEED_BASE = 20991000  # plus the case index 0-9 (Amendment 19 (f))
KE_CUT_MEV = 1.0e-4
PP_CUT_DEG = 16.26
SELECTION_RULE = "fast-proton-on-curve-v2"  # addendum (c)2: a p-p event is accepted when its fast proton lies on the p+H curve
GRAPHITE_RHO = 2.0
HBARC = 197.3269804
SIGNATURES = (  # keys of a class summary that the Geant4 hook compares
    "rate_per_primary_after_rule", "P_theta_lab_gt_theta_2mm", "P_theta_lab_gt_theta_3mm", "theta_lab_median_deg", "T_tail_moment",
    "mean_t_other_mev", "first_minimum_deg", "first_minimum_3sigma_deg", "qslope",
)  # fmt: skip
_CL = "cases/*/classes/*"
GIT_ALLOW = (
    "selection_rule", "events_per_case", "seed_base", "table_id", "code/*", "code/scripts_sha256/*", "ke_cut_mev", "classification_rule", "descriptor_note", "thin_target",
    *(f"cases/*/{k}" for k in ("case", "row", "material", "energy_mev", "length_cm", "primaries", "seed", "events_sampled",
                               "range_water_mm", "theta_2mm_deg", "theta_3mm_deg", "total_rate_per_primary",
                               "pp_share_of_elastic_events")),
    *(f"{_CL}/{k}" for k in (
        "P_theta_lab_gt_theta_2mm", "P_theta_lab_gt_theta_3mm", "T_tail_moment", "counts_in_first_minimum_bin", "event_share",
        "events", "expected_events_in_primaries", "first_minimum_3sigma_deg", "first_minimum_deg", "ke_cut_fraction_removed",
        "mean_t_other_mev", "qslope", "qslope_bins", "rate_per_primary", "rate_per_primary_after_ke_cut",
        "rate_per_primary_after_rule", "off_curve_fraction_excluded", "theta_cm_median_deg",
        "theta_cm_min_deg", "theta_lab_mean_1_minus_cos", "theta_lab_median_deg", "theta_slow_lab_median_deg")),
    f"{_CL}/classification/*", f"{_CL}/helm/*",
)  # fmt: skip
"""Git allowlist of ``de5.json`` (syntax: :mod:`common`): per-primary rates, class shares, angle statistics, exceedance
probabilities, tail moments, mean secondary energies, classification diagnostics and descriptors. Not kept: ``sigma_el_mb``
(absolute table cross sections), every histogram (8-bin and full), wall time."""

GIT_FILES = {"de5.json": (GIT_ALLOW, ())}
HIST_INTERIOR_DEG = [2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 60.0]  # 8 lab-angle bins (7 interior edges)
HIST_INTERIOR_T = [0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8]  # primary T1 / E0, 8 bins


def _build_cases(seed_base: int = SEED_BASE) -> tuple[dict[str, Any], ...]:
    out = []
    for i, (mat, e) in enumerate([(m, e) for m in ("water", "graphite") for e in ENERGIES]):
        thin = e <= 50
        row = 14 + i
        out.append({
            "case": i, "row": row, "material": mat, "energy_mev": e, "length_cm": 0.2 if thin else 1.0,
            "primaries": 20_000_000 if thin else 4_000_000, "seed": seed_base + i,
            "target": "O-16" if mat == "water" else "C-12",  # the p+A class of the material (water also has H-1)
            "topas_material": "G4_WATER" if mat == "water" else "G4_C", "topas_seed_base": 20270000 + 1000 * row,
            # frozen D-E4 design (plan 3.2): split into 3 equal parts at 30 / 50 MeV, one part otherwise; the TOPAS input is
            # bound to these constants by de4_parse.identify_case
            "part_histories": (6_666_667, 6_666_667, 6_666_666) if thin else (4_000_000,),
            "target_hl_mm": (2.0 if thin else 10.0) / 2.0, "target_rmax_mm": 10.0, "sphere_rmax_mm": 500.0,
            "topas_modules": ("g4h-elastic_HP",),
        })  # fmt: skip
    return tuple(out)


CASES = _build_cases()
"""The frozen D-E4 / D-E5 case matrix (plan 3.2, Amendment 19 (c)): case index 0-9 = rows 14-23; water rows 14-18, graphite rows
19-23, 30 / 50 / 100 / 150 / 200 MeV; thickness 2 mm at 30 / 50 MeV and 1 cm otherwise; primaries 2e7 / 4e6. The parser
(:mod:`de4_parse`) identifies a materialized TOPAS case against this table."""


def case_table(seed_base: int = SEED_BASE) -> list[dict[str, Any]]:
    """Cases 0-9 = rows 14-23 (a copy of :data:`CASES`; ``seed_base`` replaces the sampler seeds only: 20991300 for the M2' diagnostic D4,
    20991200 for the report-only M2 arm; the case matrix and the TOPAS identity are unchanged)."""
    if seed_base == SEED_BASE:
        return [dict(cs) for cs in CASES]
    return [dict(cs) for cs in _build_cases(seed_base)]


def select_tracks(
    cos_lab: np.ndarray,
    t1: np.ndarray,
    e0: float,
    m_h: float,
    m_a: float,
    has_secondary_proton: Any,
    ke_cut: float = KE_CUT_MEV,
) -> dict[str, np.ndarray]:
    """The frozen D-E4 track selection, shared by the table side and the Geant4 parser: the kinetic-energy filter removes tracks with
    ``E > E0 - ke_cut`` (undeflected primaries); of the remaining tracks those off both two-body curves (:func:`classify`) are
    counted and **excluded**. ``keep`` = neither removed nor off-curve; ``is_pp`` and ``conflict`` as :func:`classify`."""
    is_pp, conflict, off = classify(cos_lab, t1, e0, m_h, m_a, has_secondary_proton)
    ke_removed = np.asarray(t1, dtype=float) > e0 - ke_cut
    off = off & ~ke_removed
    return {
        "is_pp": is_pp,
        "conflict": conflict,
        "off": off,
        "ke_removed": ke_removed,
        "keep": ~ke_removed & ~off,
    }


def kin_curve_t(e0: float, m_t: float, cos_theta: np.ndarray) -> np.ndarray:
    """Kinetic energy of an elastically scattered proton at lab angle ``theta`` (forward root, exact relativistic):
    ``(A^2 - B^2) E1^2 - 2 A C E1 + C^2 + B^2 m^2 = 0`` with A = E0 + M, B = p0 cos(theta), C = m^2 + E0 M (E0 total energy)."""
    from ionmc.nuclear.bgg import M_PROTON_MEV as m

    et = e0 + m
    p0 = math.sqrt(et * et - m * m)
    a, b, cc = et + m_t, p0 * cos_theta, m * m + et * m_t
    disc = np.sqrt(np.maximum(cc * cc - m * m * (a * a - b * b), 0.0))
    return (a * cc + b * disc) / (a * a - b * b) - m


def classify(
    cos_lab: np.ndarray,
    t1: np.ndarray,
    e0: float,
    m_h: float,
    m_a: float,
    has_secondary_proton: Any,
    tol: float = 1.0e-3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """D-E4 per-track classification rule (table side now, Geant4 phase space later). Residuals of the primary's kinetic energy
    ``t1`` at its lab angle to the elastic kinematic curves of p+H (mass ``m_h``) and p+A (``m_a``), relative to ``t1``. The class
    with the smaller residual is assigned; a *conflict* track (both residuals < ``tol``; the curves lie within 1e-3 of each other up
    to ``theta_lab`` of about 1.8 deg) is **not** excluded: it is p-p iff the event has a secondary proton (``has_secondary_proton``,
    bool or array), otherwise p+A. Returns ``(is_pp, conflict, off_curve)``; ``off_curve`` = the smaller residual exceeds ``tol``
    (such tracks are classified by the smaller residual and counted)."""
    t1 = np.asarray(t1, dtype=float)
    rh = np.abs(t1 - kin_curve_t(e0, m_h, cos_lab)) / np.maximum(t1, 1e-9)
    ra = np.abs(t1 - kin_curve_t(e0, m_a, cos_lab)) / np.maximum(t1, 1e-9)
    conflict = (rh < tol) & (ra < tol)
    is_pp = np.where(conflict, np.asarray(has_secondary_proton, dtype=bool), rh < ra)
    return is_pp, conflict, np.minimum(rh, ra) > tol


def curve_crossover_deg(e0: float, m_h: float, m_a: float, tol: float = 1.0e-3) -> float:
    """Lab angle below which the p+H and p+A curves lie within ``tol`` (relative) of each other."""
    th = np.arange(0.0, 20.0, 0.01)
    c = np.cos(np.radians(th))
    d = np.abs(kin_curve_t(e0, m_h, c) - kin_curve_t(e0, m_a, c)) / kin_curve_t(e0, m_h, c)
    bad = np.nonzero(d >= tol)[0]
    return float(th[bad[0]]) if bad.size else float("nan")


def first_local_min_deg(
    th: np.ndarray, f: np.ndarray, above: float
) -> tuple[float | None, float | None]:
    for i in range(1, th.size - 1):
        if th[i] > above and f[i] < f[i - 1] and f[i] < f[i + 1]:
            return float(th[i]), float(f[i])
    return None, None


def _j1s(x: np.ndarray) -> np.ndarray:
    """Spherical Bessel j1(x) = (sin x - x cos x) / x^2."""
    x = np.asarray(x, dtype=float)
    small = x < 1e-3
    xs = np.where(small, 1.0, x)
    return np.where(small, x / 3.0, (np.sin(xs) - xs * np.cos(xs)) / xs**2)


def helm_fit(q: np.ndarray, dens: np.ndarray) -> dict[str, float]:
    """Report-only least-squares fit of ``|3 j1(qR)/(qR) exp(-q^2 s^2 / 2)|^2`` (free amplitude) to ``ln dens(q)``."""
    best = (math.inf, 0.0, 0.0)
    for r in np.arange(2.0, 6.0001, 0.05):
        for s in np.arange(0.0, 1.5001, 0.05):
            m = np.log(
                np.maximum((3.0 * _j1s(q * r) / (q * r) * np.exp(-0.5 * (q * s) ** 2)) ** 2, 1e-300)
            )
            d = np.log(dens) - m
            chi = float(np.sum((d - d.mean()) ** 2))
            if chi < best[0]:
                best = (chi, float(r), float(s))
    return {
        "R_eff_fm": best[1],
        "s_fm": best[2],
        "rms_log_residual": math.sqrt(best[0] / q.size),
        "n_bins": int(q.size),
    }


def shape_descriptors(mu: np.ndarray, pcm: float) -> dict[str, Any]:
    """First minimum (1-deg CM bins above 3 deg), large-q log slope and Helm fit of a p+A CM sample."""
    th = np.degrees(np.arccos(np.clip(mu, -1.0, 1.0)))
    edges = np.arange(0.0, 181.0, 1.0)
    cnt, _ = np.histogram(th, bins=edges)
    dens = cnt / (np.cos(np.radians(edges[:-1])) - np.cos(np.radians(edges[1:])))
    cen = 0.5 * (edges[1:] + edges[:-1])
    tmin, _ = first_local_min_deg(cen, dens.astype(float), 3.0)
    # the frozen estimator picks any noise dip; the 3-sigma variant needs a Poisson dip of >= 3 sigma against both neighbours
    sig_min = None
    for i in range(1, cen.size - 1):
        if cen[i] > 3.0 and cnt[i] > 0 and dens[i] < dens[i - 1] and dens[i] < dens[i + 1]:
            w = np.cos(np.radians(edges[i])) - np.cos(np.radians(edges[i + 1]))
            err = math.sqrt(cnt[i]) / w
            dips = [
                (dens[j] - dens[i])
                / math.hypot(
                    err,
                    math.sqrt(max(cnt[j], 1))
                    / (np.cos(np.radians(edges[j])) - np.cos(np.radians(edges[j + 1]))),
                )
                for j in (i - 1, i + 1)
            ]
            if min(dips) >= 3.0:
                sig_min = float(cen[i])
                break
    out: dict[str, Any] = {"first_minimum_deg": tmin, "first_minimum_3sigma_deg": sig_min,
                           "counts_in_first_minimum_bin": None if tmin is None else int(cnt[int(tmin)]),
                           "qslope": None, "qslope_bins": 0, "helm": None}  # fmt: skip
    if tmin is None:
        return out
    k = pcm / HBARC
    q_lo, q_hi = 2 * k * math.sin(math.radians(tmin) / 2), 2 * k * math.sin(math.radians(90.0) / 2)
    qs = 2 * k * np.sin(np.radians(th) / 2)
    qe = np.linspace(q_lo, q_hi, 13)
    qc, _ = np.histogram(qs, bins=qe)
    dq = qe[1:] ** 2 - qe[:-1] ** 2
    qcen = 0.5 * (qe[1:] + qe[:-1])
    ok = qc >= 30
    if ok.sum() >= 3:
        slope = np.polyfit(np.log(qcen[ok]), np.log(qc[ok] / dq[ok]), 1)[0]
        out["qslope"], out["qslope_bins"] = float(slope), int(ok.sum())
    qh = 2 * k * np.sin(np.radians(np.linspace(3.0, 90.0, 61)) / 2)
    qe2 = np.linspace(qh[0], qh[-1], 25)
    c2, _ = np.histogram(qs, bins=qe2)
    d2 = c2 / (qe2[1:] ** 2 - qe2[:-1] ** 2)
    ok2 = c2 >= 50
    if ok2.sum() >= 6:
        out["helm"] = helm_fit(0.5 * (qe2[1:] + qe2[:-1])[ok2], d2[ok2])
    return out


def sample_case(tab: Any, water: c.Water, case: dict[str, Any], seed_base: int = SEED_BASE) -> dict[str, Any]:
    """The table sampled for one case; see the module docstring."""
    from ionmc.nuclear.bgg import M_PROTON_MEV as MP
    from ionmc.nuclear.elastic_kin import two_body

    e = float(case["energy_mev"])
    names = tuple(tab.target_names)
    assert seed_base <= case["seed"] < seed_base + 10 and 20990000 <= case["seed"] < 21000000
    if case["material"] == "water":
        dens = {n: t["n_per_cm3"] * t["scale"] for n, t in water.targets.items()}
    else:
        el = tab.info["elements"]["C"]
        dens = {
            str(el["target"]): c.N_AVOGADRO
            * GRAPHITE_RHO
            / float(el["a_g_mol"])
            * float(el["sigma_scale"])
        }
    cls = sorted(dens, key=names.index)
    sig = {n: float(tab.sigma_barn(names.index(n), e)) * 1e-24 for n in cls}  # cm2
    ns = {n: dens[n] * sig[n] for n in cls}  # 1/cm
    tot = sum(ns.values())
    rng = np.random.default_rng(np.random.SeedSequence(case["seed"]))
    pick = rng.choice(len(cls), size=N_EVENTS, p=[ns[n] / tot for n in cls])
    r_w = float(water.range_mm(e))
    th_w = {w: math.degrees(math.acos(1.0 - w / r_w)) for w in (2.0, 3.0)}
    out: dict[str, Any] = {"case": case["case"], "row": case["row"], "material": case["material"], "energy_mev": e,
                           "length_cm": case["length_cm"], "primaries": case["primaries"], "seed": case["seed"],
                           "events_sampled": N_EVENTS, "range_water_mm": r_w, "theta_2mm_deg": th_w[2.0], "theta_3mm_deg": th_w[3.0],
                           "sigma_el_mb": {n: sig[n] * 1e27 for n in cls}, "classes": {}}  # fmt: skip
    for k, name in enumerate(cls):
        it = names.index(name)
        n_k = int((pick == k).sum())
        m_t = float(tab.arrays["target_mass_mev"][it])
        mu = tab.sample_mu(it, e, rng.random(n_k))
        pp = name == "H-1"
        if pp:
            assert np.abs(mu).max() <= 0.96 + 1e-12  # theta_CM >= 16.26 deg
        tb = two_body(e, MP, m_t, mu)
        if pp:
            fast = tb.t1 >= tb.t2
            cos2 = tb.pz2 / np.hypot(tb.px2, tb.pz2)
            t1, cos1 = np.where(fast, tb.t1, tb.t2), np.where(fast, tb.cos_lab1, cos2)
            t_other = np.where(fast, tb.t2, tb.t1)
            cos_other = np.where(fast, cos2, tb.cos_lab1)
        else:
            t1, cos1, t_other = tb.t1, tb.cos_lab1, tb.t2
        mh = float(tab.arrays["target_mass_mev"][names.index("H-1")])
        m_a = (
            m_t
            if not pp
            else float(
                tab.arrays["target_mass_mev"][
                    names.index(next((n for n in cls if n != "H-1"), "O-16"))
                ]
            )
        )
        sel = select_tracks(
            cos1, t1, e, mh, m_a, has_secondary_proton=pp
        )  # KE filter, classification, off-curve exclusion
        keep = sel["keep"]
        removed = float(sel["ke_removed"].mean())
        surv = ~sel["ke_removed"]
        off_frac = float(sel["off"][surv].mean()) if surv.any() else 0.0
        th = np.degrees(np.arccos(np.clip(cos1, -1, 1)))[keep]
        t1k, cos1k, tok, mu_k = (
            t1[keep],
            cos1[keep],
            t_other[keep],
            mu[keep],
        )  # the selected events (KE filter, exclusion)
        th_med = float(np.median(th))
        rate = dens[name] * sig[name] * case["length_cm"]
        h_th, _ = np.histogram(th, bins=[0.0, *HIST_INTERIOR_DEG, 180.0])
        h_t, _ = np.histogram(t1k / e, bins=[0.0, *HIST_INTERIOR_T, 1.0 + 1e-9])
        cls_doc: dict[str, Any] = {
            "events": n_k, "event_share": n_k / N_EVENTS, "ke_cut_fraction_removed": removed,
            "rate_per_primary": rate, "rate_per_primary_after_ke_cut": rate * (1.0 - removed),
            "off_curve_fraction_excluded": off_frac,
            "expected_events_in_primaries": rate * case["primaries"],
            "P_theta_lab_gt_theta_2mm": float(np.mean(th > th_w[2.0])), "P_theta_lab_gt_theta_3mm": float(np.mean(th > th_w[3.0])),
            "theta_lab_median_deg": th_med, "theta_lab_mean_1_minus_cos": float(np.mean(1.0 - cos1k)),
            "T_tail_moment": float(np.mean((1.0 - cos1k) * (th > 3.0 * th_med))), "mean_t_other_mev": float(tok.mean()),
            "theta_cm_median_deg": float(np.median(np.degrees(np.arccos(np.clip(np.abs(mu_k) if pp else mu_k, -1, 1))))),
            "hist_theta_lab_counts_8bins": h_th.tolist(), "hist_theta_lab_interior_edges_deg": HIST_INTERIOR_DEG,
            "hist_t1_over_e0_counts_8bins": h_t.tolist(), "hist_t1_over_e0_interior_edges": HIST_INTERIOR_T,
            "hist_theta_lab_1deg_full": np.histogram(th, bins=np.arange(0.0, 181.0))[0].tolist(),
        }  # fmt: skip
        if pp:
            cls_doc["theta_cm_min_deg"] = float(np.degrees(np.arccos(np.abs(mu_k).max())))
            cls_doc["theta_slow_lab_median_deg"] = float(
                np.median(np.degrees(np.arccos(np.clip(cos_other, -1, 1))))
            )
            cls_doc["first_minimum_deg"] = cls_doc["first_minimum_3sigma_deg"] = cls_doc[
                "qslope"
            ] = None
        else:
            sd = shape_descriptors(
                mu_k, float(np.interp(e, tab.arrays["grid_e_mev"], tab.arrays["p_cm_mev"][it]))
            )
            cls_doc.update(sd)
        # D-E4 classification rule on the kept table events
        is_pp, conflict = sel["is_pp"][keep], sel["conflict"][surv]
        assigned_own = is_pp if pp else ~is_pp
        own_frac = float(assigned_own.mean())
        cls_doc["rate_per_primary_after_rule"] = (
            rate * (1.0 - removed) * (1.0 - off_frac) * own_frac
        )
        cls_doc["classification"] = {
            "conflict_fraction_raw": float(conflict.mean()), "assigned_to_own_class_fraction": own_frac,
            "misclassified_fraction": 1.0 - own_frac, "off_curve_fraction": off_frac,
            "rate_per_primary_after_rule": cls_doc["rate_per_primary_after_rule"],
            "curve_crossover_deg": curve_crossover_deg(e, mh, m_a),
        }  # fmt: skip
        out["classes"][name] = cls_doc
    out["total_rate_per_primary"] = sum(v["rate_per_primary"] for v in out["classes"].values())
    if "H-1" in out["classes"]:
        out["pp_share_of_elastic_events"] = out["classes"]["H-1"]["event_share"]
    return out


def signature_ratios(
    table_doc: dict[str, Any], geant4: dict[str, dict[str, dict[str, float]]]
) -> dict[str, Any]:
    """Hook for the Geant4 comparison (plan D-E5): ``geant4[str(case)][class][signature]`` (the keys of :data:`SIGNATURES` of a
    D-E4 per-class summary) -> ``ionmc / Geant4`` per case, class and signature. A missing case, class or signature key on either
    side raises ``KeyError`` (nothing is silently dropped) and a Geant4 case summary whose ``_selection_rule`` is not :data:`SELECTION_RULE`
    raises ``ValueError``; a value that is ``None`` or a zero Geant4 value gives ``None``."""
    out: dict[str, Any] = {}
    for key, case in table_doc["cases"].items():
        if key not in geant4:
            raise KeyError(f"Geant4 summary has no case {key}")
        if geant4[key].get("_selection_rule") != SELECTION_RULE:
            raise ValueError(
                f"Geant4 summary of case {key} was not made under the selection rule {SELECTION_RULE!r}"
            )
        out[key] = {}
        for name, cls in case["classes"].items():
            if name not in geant4[key]:
                raise KeyError(f"Geant4 summary of case {key} has no class {name}")
            g = geant4[key][name]
            row = {}
            for s in SIGNATURES:
                if s not in cls:
                    raise KeyError(f"table record of case {key}, class {name} has no signature {s}")
                if s not in g:
                    raise KeyError(
                        f"Geant4 summary of case {key}, class {name} has no signature {s}"
                    )
                a, b = cls[s], g[s]
                row[s] = None if a is None or not b else float(a) / float(b)
            out[key][name] = row
    return out


D4_SEED_BASES = {"m2prime": 20991300, "m2": 20991200}
"""Amendment 21 (e): D4 seeds of the model under test (M2', 20991300 + case) and of the report-only M2 arm (20991200 + case)."""


def run(table_id: str | None = None, cache_dir: Path | None = None, seed_base: int = SEED_BASE, tag: str | None = None) -> None:
    """The D-E5 table-side record. ``table_id`` (in the scratch ``cache_dir``) samples a variant / builder table with ``seed_base``; the record
    then goes to ``de5-<tag>/`` (the committed pinned-table ``de5/`` record is never overwritten) and carries ``table_id``."""
    t0 = time.perf_counter()
    tab = c.load_table(table_id, cache_dir)
    water = c.Water(tab)
    cases = {str(cs["case"]): sample_case(tab, water, cs, seed_base) for cs in case_table(seed_base)}
    doc = {"banner": c.BANNER, "selection_rule": SELECTION_RULE, "events_per_case": N_EVENTS, "seed_base": seed_base, "ke_cut_mev": KE_CUT_MEV, "cases": cases,
           "descriptor_note": "first-minimum, q-slope and Helm descriptors are recorded, not interpreted (frozen 1-deg estimator picks noise dips in a 256-bin piecewise-uniform table; the 3-sigma variant misses the J1 zero)",
        "classification_rule": "conflict tracks (both residuals < 1e-3, theta_lab below the ~1.8 deg crossover) are p-p iff the event has a secondary proton, else p+A",
        "thin_target": "fixed beam energy, no slowing down, no multiple scattering (as the D-E4 EM-free cases)",
           "wall_s": time.perf_counter() - t0}  # fmt: skip
    if table_id is not None:
        doc["table_id"] = table_id
        doc["code"] = c.code_provenance()
    c.write_outputs(
        c.OUT / (f"de5-{tag}" if tag else "de5"), "de5.json", doc, "de5-report.md", _header(table_id, cache_dir, seed_base), render_report, GIT_ALLOW
    )


def _header(table_id: str | None = None, cache_dir: Path | None = None, seed_base: int = SEED_BASE) -> str:
    cmd = "python -m validation.scripts.elastic_attribution.de5_sampling" + (f" --table-id {table_id}" if table_id else "")
    return c.header(
        "D-E5 (table side): the pinned elastic table sampled in the D-E4 frame",
        cmd,
        c.table_files(table_id, cache_dir),
        [
            f"1e6 events per case from the P7 sampler (`ElasticTable.sample_mu`), seeds SeedSequence({seed_base} + case), thin target "
            "at the beam energy, exact two-body kinematics, p-p cut at 16.26 deg CM, D-E4 kinetic-energy filter. Signature "
            "definitions: module docstring (plan 2.2). Thresholds theta_w use the Bethe range of water."
        ],
    )


def render_report(doc: dict[str, Any]) -> str:
    out = ["", "## Per case and class", "",
           "| case (row) | material | E | class | rate/primary | share | KE cut removed | theta_med,lab | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV | first min (deg; 3-sigma) | q-slope | Helm R_eff / s |",
           "|" + "---|" * 15]  # fmt: skip
    for k, cs in doc["cases"].items():
        for name, v in cs["classes"].items():
            h = v.get("helm")
            out.append(
                f"| {k} ({cs['row']}) | {cs['material']} | {cs['energy_mev']:g} | {name} | {v['rate_per_primary']:.4g} | "
                f"{v['event_share']:.3f} | {v['ke_cut_fraction_removed']:.2e} | {v['theta_lab_median_deg']:.2f} | {v['P_theta_lab_gt_theta_2mm']:.4f} | "
                f"{v['P_theta_lab_gt_theta_3mm']:.4f} | {v['T_tail_moment']:.4f} | {v['mean_t_other_mev']:.3f} | "
                f"{'-' if v.get('first_minimum_deg') is None else format(v['first_minimum_deg'], '.1f')}; {'-' if v.get('first_minimum_3sigma_deg') is None else format(v['first_minimum_3sigma_deg'], '.1f')} | "
                f"{'-' if v.get('qslope') is None else format(v['qslope'], '.2f')} | "
                f"{'-' if not h else format(h['R_eff_fm'], '.2f') + ' / ' + format(h['s_fm'], '.2f')} |"
            )
    out += [
        "",
        "## Thresholds, p-p cut and the classification rule",
        "",
        "Conflict tracks (both kinematic residuals < 1e-3; the p+H and p+A curves agree to 1e-3 below the crossover angle listed) "
        "are not excluded: p-p iff the event has a secondary proton, else p+A.",
        "",
        "| case | theta_2mm | theta_3mm | class | theta_CM min (p-p) | crossover (deg) | raw conflict | misclassified | off curve | rate/primary after rule |",
        "|" + "---|" * 10,
    ]
    for k, cs in doc["cases"].items():
        for name, v in cs["classes"].items():
            cl = v["classification"]
            out.append(
                f"| {k} | {cs['theta_2mm_deg']:.2f} | {cs['theta_3mm_deg']:.2f} | {name} | "
                f"{'-' if 'theta_cm_min_deg' not in v else format(v['theta_cm_min_deg'], '.2f')} | {cl['curve_crossover_deg']:.2f} | {cl['conflict_fraction_raw']:.2e} | "
                f"{cl['misclassified_fraction']:.2e} | {cl['off_curve_fraction']:.2e} | {cl['rate_per_primary_after_rule']:.4g} |"
            )
    out += [
        "",
        "The first-minimum, q-slope and Helm descriptors are recorded, not interpreted: the frozen 1-degree estimator picks noise dips "
        "in a 256-bin piecewise-uniform table, and the 3-sigma variant misses the J1 zero.",
        "",
        "Absolute cross sections and angle / energy histograms are not committed (decision 0041); they are in the generated record. "
        "Geant4 ratios (ionmc / Geant4) are produced later by `signature_ratios` from the D-E4 summary.",
    ]
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> None:
    import argparse

    ap = argparse.ArgumentParser(description="D-E5 table side (default: the pinned table, seeds 20991000 + case)")
    ap.add_argument("--table-id", default=None, help="a variant / builder table id installed in --scratch-cache")
    ap.add_argument("--scratch-cache", type=Path, default=None)
    ap.add_argument("--seed-base", type=int, default=SEED_BASE, help="m2prime D4: 20991300; m2 (report-only): 20991200")
    ap.add_argument("--tag", default=None, help="output directory suffix (de5-<tag>); required with --table-id")
    a = ap.parse_args(argv)
    if a.table_id and not a.tag:
        raise SystemExit("--table-id needs --tag (the committed de5/ record is never overwritten by a variant)")
    if a.table_id and a.seed_base == SEED_BASE:
        raise SystemExit(f"a variant table needs its own seed base, one of {sorted(D4_SEED_BASES.values())}")
    run(a.table_id, a.scratch_cache, a.seed_base, a.tag)


if __name__ == "__main__":
    main()
