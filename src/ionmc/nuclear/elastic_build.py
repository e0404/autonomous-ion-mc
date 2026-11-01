"""Builder of the derived proton elastic-scattering table (decision 0041 slice C, V3-005C).

``build_elastic_proton(cache_dir, options)`` reads the hash-pinned sources from the cache (the
ENDF/B-VIII.0 proton sublibrary = LA150, AME2020 and the Geant4 11.4.2 Barashenkov/BGG files, whose
hand-transcribed numbers are compared with the pinned files) and writes
``<cache>/derived/elastic-proton-<id>.npz`` and ``.json``. Single process, numpy and stdlib only,
deterministic (the id is the sha256 of the canonical sidecar, as for the nuclear table).

Targets (axis 0): ``H-1`` followed by the seven LA150 targets of the nuclear table (same order, same
``ELEMENT_TARGET`` map and ``(A_el/A_ref)^(2/3)`` scaling for surrogate elements; hydrogen is
its own
target). Grid: the union grid convention of the nuclear table (1 MeV, the ENDF H-1 LAW=5 nodes <=
150 MeV, a uniform ln E grid with >= 50 points per decade) plus every Barashenkov node and 14
MeV, cut
at 250 MeV, plus adaptive midpoint refinement below 14 MeV where the BGG Coulomb-factor rule is
nonlinear (``sigma_el`` lin-lin between nodes reproduces the direct evaluation to 5e-4).

p + A (C, N, O, Al, Si, P, Ca), 1-250 MeV
    ``sigma_el(E)``: BGG/Barashenkov (:mod:`ionmc.nuclear.bgg`). Shape (Amendment 14 (b)): "the
    nuclear-only black-disk form dsigma/dOmega_CM proportional to |2 J1(qR)/(qR)|^2,
    q = 2 p_CM sin(theta_CM/2). R comes from pi (R + lambdabar)^2 = sigma_nonel(E) of the
    transport's
    own table (LA150 MT5 <= 150 MeV, Tripathi-shape extension above), so no elastic data enter the
    shape." ``lambdabar = hbar/p_CM``. The inverse CDF of ``mu_CM`` over [-1, 1] is stored as
    the 257
    edges of 256 equiprobable bins per node; R, ``sigma_nonel`` and the inversion residual are
    recorded per node. A node with ``R <= 0`` or ``sigma_nonel = 0`` has no diffraction radius:
    below
    ``e_min_shape`` (the lowest node from which all nodes are valid) ``sigma_el`` is set to 0
    and the
    value is recorded (analogous to ``e_min_pp``; the proton range there is below 1 mm).

H-1
    <= 150 MeV the LA150 Hale LAW=5 LTP=1 reconstruction (:mod:`ionmc.nuclear.law5`): transported
    density NI = sigma_e - sigma_c over ``|mu_CM| <= 0.96`` (ratio interpolation in E), ``sigma`` =
    2 pi int_0^0.96 NI dmu in barn (HALF sphere: identical protons, each event once), the edges
    cover
    the symmetric range [-0.96, 0.96] of the same density. ``e_min_pp`` = the lowest grid node from
    which every node and node midpoint has a non-negative density over |mu| <= 0.96 (decision 0041,
    ratified); ``sigma = 0`` below it. The frozen row P6 FAILS for H-1 below ``e_min_pp``: this is
    kept as a recorded failure (``p6_negative_density_found``, all negative nodes/midpoints with
    their values), the revised-domain check (``no_negative_density_in_domain``) is separate.
    150-250 MeV: "sigma(E) = sigma_NI(150) S(E)/S(150), with S = 1.0115 (23 + 50 sqrt(ln(0.73/p)^7))
    mb (Geant4 ``G4HadronNucleonXsc`` p-p formula, p_lab in GeV/c; no PDG data enter). The shape is
    held fixed in mu_CM."

X-ENDF (report-only) and the O-16 finding
    LA150 C-12/N-14/Ca-40 LTP=12 NI densities (ratio interpolation in mu) at 20-40 deg for
    50/100/150
    MeV are stored for the cross-check row. The builder asserts, fail closed, that the O-16 MT2 data
    are a numerical copy of C-12 above 24 MeV (decision 0041 data finding), that the LA150 O-16
    MF3/MT5 sigma_nonel is NOT a copy of C-12's (``o16_mt5_copy_check``: the R rule and the
    nuclear table use MT5) and never uses LAW=5 data of any target except H-1 as a construction
    input.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.data import cache, endf6
from ionmc.data.ame import load_ame2020
from ionmc.data.registry import DATASETS
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import ELEMENTS
from ionmc.nuclear import bgg, law5
from ionmc.nuclear import elastic_glauber as glauber
from ionmc.nuclear import elastic_m3 as m3mod
from ionmc.nuclear import events as ev
from ionmc.nuclear.build import (
    TARGETS,
    BuildError,
    SpeciesTables,
    TargetTables,
    element_rows,
    table_id,
    union_grid,
    write_npz_deterministic,
)

BUILDER_VERSION = "ionmc-elastic-proton-builder-5"
MF6_MIN_INFORMATIVE = 10
SCHEMA = "ionmc-elastic-proton-table-2"
SOURCE_IDS = (
    "endf-b8.0-protons",
    "ame2020-mass",
    "geant4-g4barashenkovdata-hh-11.4.2",
    "geant4-g4bggnucleonelasticxs-cc-11.4.2",
    "geant4-g4nucleonnuclearcrosssection-cc-11.4.2",
    "geant4-g4componentbarnucleonnucleusxsc-cc-11.4.2",
    "geant4-g4hadronnucleonxsc-cc-11.4.2",
    "geant4-g4nuclearradii-cc-11.4.2",
    "geant4-g4isotopelist-hh-11.4.2",
)
"""Construction sources (Amendment 15 (c)): S(E) uses the Geant4 p-p formula only; the PDG p-p
compilation is the report-only V10 comparison and is not a source of this table."""
E_MAX_MEV = 250.0
E_ANCHOR_MEV = 150.0
N_Q = 256
HBARC = 197.3269804
MP = 938.27208816
PP_TARGET = "H-1"
TARGET_NAMES = (PP_TARGET, *(t.name for t in TARGETS))
QUALIFICATION_FIELDS = (
    "no_negative_density_in_domain",
    "o16_mt5_sigma_not_c12_copy",
    "o16_mt5_spectra_independent",
    "shape_normalised",
)
E_MIN_PP_LIMIT_MEV = 15.0
"""The builder fails if E_min,pp exceeds this (Amendment 15 (a)3)."""
E_MIN_SHAPE_LIMIT_MEV = 10.0
"""The builder fails if any target's ``e_min_shape`` exceeds this (declared low-energy limitation
of p + A elastic scattering, Amendment 15 (b)5)."""
N_H_WATER_CM3 = 6.69e22
N_O_WATER_CM3 = 3.34e22
MODEL_REVISIONS = (
    "Amendment 15: P6 negativity of H-1 below E_min,pp is a recorded failure of the "
    "frozen row, not a pass; the table domain is revised to [E_min,pp, 250] MeV (p-p) and "
    "[e_min_shape(target), 250] MeV (p + A); S(E) uses the Geant4 p-p formula only (PDG data "
    "are report-only); the o16_mt5_copy_check is implemented. Amendment 16: the frozen MF6 "
    "per-product rule stays recorded as failed; the flag o16_mt5_not_c12_copy is replaced by "
    "o16_mt5_sigma_not_c12_copy and o16_mt5_spectra_independent (KS test of the product spectra)."
)
XENDF_TARGETS = ("C-12", "N-14", "Ca-40")
XENDF_ENERGIES_MEV = (50.0, 100.0, 150.0)
XENDF_THETA_DEG = tuple(float(x) for x in np.arange(20.0, 40.0001, 2.5))
REFINE_TOL = 5.0e-4
MAX_REFINE_ROUNDS = 10


@dataclass(frozen=True)
class ElasticBuildOptions:
    """Builder options (the canonical JSON enters the table id)."""

    points_per_decade: int = 50
    n_quantiles: int = N_Q
    n_shape: int = 16385
    n_mu_check: int = 4001

    def canonical(self) -> str:
        return json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))


@dataclass(frozen=True)
class ElasticBuildResult:
    npz_path: Path
    json_path: Path
    table_id: str
    info: dict[str, Any]


def element_rows_elastic() -> dict[str, dict[str, Any]]:
    """``element_rows`` of the nuclear table plus hydrogen (own target ``H-1``, scale 1)."""
    rows = element_rows(TARGETS)
    rows["H"] = {
        "target": PP_TARGET, "surrogate": False, "a_g_mol": ELEMENTS["H"].A_g_mol,
        "a_ref": 1, "sigma_scale": 1.0,
    }  # fmt: skip
    return rows


# ---------------------------------------------------------------------------------------------
# Bessel black-disk shape
# ---------------------------------------------------------------------------------------------
def bessel_j1(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """J1(x) = (1/2 pi) int_0^2pi cos(t - x sin t) dt by the periodic trapezoid rule (exponentially
    convergent for N >> x; numpy only, no scipy)."""
    xa = np.asarray(x, dtype=np.float64)
    n = int(2 * float(xa.max(initial=0.0))) + 64
    t = 2.0 * np.pi * np.arange(n) / n
    out = np.zeros_like(xa)
    for ti in t:
        out += np.cos(ti - xa * math.sin(ti))
    return out / n


def disk_amplitude_sq(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """``|2 J1(x)/x|^2`` (1 at x = 0)."""
    xa = np.asarray(x, dtype=np.float64)
    safe = np.where(xa < 1e-8, 1.0, xa)
    a = np.where(xa < 1e-8, 1.0, 2.0 * bessel_j1(safe) / safe)
    return np.asarray(a * a)


def kinematics_cm(t_mev: float, m_t_mev: float) -> tuple[float, float]:
    """``(p_CM [MeV/c], sqrt(s) [MeV])`` of a proton of kinetic energy ``t_mev`` on ``m_t_mev``."""
    e1 = t_mev + MP
    p = math.sqrt(e1 * e1 - MP * MP)
    rs = math.sqrt(MP * MP + m_t_mev**2 + 2.0 * e1 * m_t_mev)
    return p * m_t_mev / rs, rs


def disk_radius_fm(sigma_nonel_b: float, t_mev: float, m_t_mev: float) -> tuple[float, float]:
    """``(R, lambdabar)`` [fm]: ``pi (R + lambdabar)^2 = sigma_nonel`` (1 b = 100 fm^2), ``lambdabar
    = hbar / p_CM``."""
    pcm, _ = kinematics_cm(t_mev, m_t_mev)
    lam = HBARC / pcm
    return math.sqrt(100.0 * sigma_nonel_b / math.pi) - lam, lam


def disk_quantiles(
    r_fm: float, pcm_mev: float, n_s: int, n_q: int
) -> tuple[NDArray[np.float64], float]:
    """Edges (n_q + 1, ascending mu_CM from -1 to 1) of the ``n_q`` equiprobable bins of
    ``dsigma/dOmega proportional to |2 J1(qR)/(qR)|^2`` and the relative deviation of the
    cumulative from 1 (normalisation check). In s = sin(theta/2): mu = 1 - 2 s^2, q = (2 p/hbar) s,
    density per ds = 4 s |.|^2; the CDF is the trapezoid sum on ``n_s`` nodes, inverted linearly."""
    s = np.linspace(0.0, 1.0, n_s)
    g = 4.0 * s * disk_amplitude_sq(2.0 * pcm_mev / HBARC * s * r_fm)
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(s))))
    cum /= cum[-1]
    target = 1.0 - np.arange(n_q + 1) / n_q  # F(s) = 1 - i/n: ascending mu
    s_q = np.interp(target, cum, s)
    mu = 1.0 - 2.0 * s_q * s_q
    mu[0], mu[-1] = -1.0, 1.0
    return np.maximum.accumulate(mu), float(abs(cum[-1] - 1.0))


def symmetric_quantiles(
    mu: NDArray[np.float64], dens: NDArray[np.float64], n_q: int
) -> NDArray[np.float64]:
    """Edges of ``n_q`` equiprobable bins of a density tabulated on the ascending ``mu`` grid (CDF
    by trapezoid, inverted linearly)."""
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(mu))))
    cum /= cum[-1]
    edges = np.interp(np.arange(n_q + 1) / n_q, cum, mu)
    edges[0], edges[-1] = mu[0], mu[-1]
    return np.maximum.accumulate(edges)


# ---------------------------------------------------------------------------------------------
# transcription check and O-16 finding
# ---------------------------------------------------------------------------------------------
def _array(text: str, name: str) -> list[float]:
    m = re.search(
        r"G4double (?:G4NuclearRadii::)?" + name + r"\s*\[\d*\]\s*=\s*\{(.*?)\};", text, re.S
    )
    if m is None:
        raise BuildError(f"array {name!r} not found in a pinned Geant4 source")
    body = re.sub(r"//.*", "", m.group(1))
    return [float(x) for x in re.split(r"[\s,]+", body) if x]


def verify_transcription(cdir: Path) -> dict[str, bool]:
    """Compare the arrays hand-transcribed in :mod:`ionmc.nuclear.bgg` with the pinned Geant4
    files; ``BuildError`` on any difference."""
    bar = cache.verify("geant4-g4barashenkovdata-hh-11.4.2", cdir).read_text(encoding="utf-8")
    iso = cache.verify("geant4-g4isotopelist-hh-11.4.2", cdir).read_text(encoding="utf-8")
    rad = cache.verify("geant4-g4nuclearradii-cc-11.4.2", cdir).read_text(encoding="utf-8")
    ok = {
        "e1": list(bgg._E1_GEV) == _array(bar, "e1"),
        "e2": list(bgg._E2_GEV) == _array(bar, "e2"),
        "e3": list(bgg._E3_GEV) == _array(bar, "e3"),
    }
    for el, z in (("c", 6), ("n", 7), ("o", 8), ("al", 13), ("si", 14), ("ca", 20)):
        ok[f"{el}_m_t"] = list(bgg.BARASHENKOV_M_T[z]) == _array(bar, f"{el}_m_t")
        ok[f"{el}_p_in"] = list(bgg.BARASHENKOV_P_IN[z]) == _array(bar, f"{el}_p_in")
    aeff, r0 = _array(iso, "aeff"), _array(rad, "r0")
    ok["aeff"] = all(aeff[z] == v for z, v in bgg.AEFF.items())
    ok["r0"] = all(r0[z] == v for z, v in bgg.R0_CB.items())
    bad = [k for k, v in ok.items() if not v]
    if bad:
        raise BuildError(
            f"transcribed Barashenkov/BGG numbers differ from the pinned sources: {bad}"
        )
    return ok


def o16_copy_finding(c12: endf6.EndfMaterial, o16: endf6.EndfMaterial) -> dict[str, Any]:
    """Compare MF3/MT2 and the MF6 LAW=5 LTP=12 P_NI tables of O-16 with C-12 at every incident
    energy >= 24 MeV of the common LAW=5 grid (decision 0041 data finding: numerical copy). Returns
    the finding; raises :class:`BuildError` when the data are NOT a copy (the decision no longer
    describes the data and must be revisited before any table is built)."""
    sc, so = law5.parse_law5(c12.sections[(6, 2)]), law5.parse_law5(o16.sections[(6, 2)])
    tc, to = c12.cross_section(2), o16.cross_section(2)
    ec, eo = sc.energies_ev, so.energies_ev
    if not np.array_equal(ec, eo):
        raise BuildError("O-16 and C-12 LAW=5 energy grids differ: not the documented copy")
    max_mf3 = max_p = 0.0
    n_nodes = 0
    first = None
    for i, e in enumerate(ec):
        if e < 24.0e6:
            continue
        rc, ro = sc.records[i], so.records[i]
        vc = float(tc.interpolate(e)[0])
        vo = float(to.interpolate(e)[0])
        max_mf3 = max(max_mf3, abs(vc - vo) / max(abs(vo), 1e-30))
        pc, po = rc.data[1::2], ro.data[1::2]
        if rc.data.size != ro.data.size or not np.array_equal(rc.data[0::2], ro.data[0::2]):
            raise BuildError("O-16 and C-12 mu grids differ: not the documented copy")
        max_p = max(max_p, float(np.max(np.abs(pc - po))) / float(np.max(np.abs(pc))))
        n_nodes += 1
        first = e * 1e-6 if first is None else first
    out = {
        "o16_mt2_is_c12_copy": bool(max_mf3 <= 1e-4 and max_p <= 1e-4),
        "from_mev": first, "nodes_compared": n_nodes,
        "max_relative_difference_mf3": max_mf3, "max_difference_over_max_p": max_p,
        "tolerance": 1e-4,
        "consequence": "LA150 O-16 MT2 is never a construction input (C1-ext); p+O elastic uses "
        "BGG sigma_el and the black-disk shape",
    }  # fmt: skip
    if not out["o16_mt2_is_c12_copy"]:
        raise BuildError(f"O-16 MT2 is not a C-12 copy any more: {out}")
    return out


REL_DIFF_THRESHOLD = 1.0e-3
MIN_SHARED_NODES = 10
E_COMPARE_MAX_EV = 150.0e6


def _first_positive_ev(tab: endf6.Tab1) -> float:
    pos = tab.x[tab.y > 0.0]
    if pos.size == 0:
        raise BuildError("MF3/MT5 sigma_nonel is nowhere positive")
    return float(pos[0])


def _mt5_compare_grid(tc: endf6.Tab1, to: endf6.Tab1) -> tuple[NDArray[np.float64], str]:
    """Incident energies [eV] of the O-16 / C-12 MT5 comparison (Amendment 15 (d)1): the shared
    nodes in ``[max threshold, 150 MeV]``; with fewer than 10 of them, the union of both node sets
    in the same interval (linear interpolation, INT=2, evaluated by the callers)."""
    lo = max(_first_positive_ev(tc), _first_positive_ev(to))
    hi = min(E_COMPARE_MAX_EV, float(tc.x[-1]), float(to.x[-1]))
    shared = np.intersect1d(tc.x, to.x)
    shared = shared[(shared >= lo) & (shared <= hi)]
    if shared.size >= MIN_SHARED_NODES:
        return shared, "shared"
    union = np.union1d(tc.x, to.x)
    union = union[(union >= lo) & (union <= hi)]
    if union.size < 2:
        raise BuildError("O-16 / C-12 MF3/MT5: fewer than two comparable energies")
    return union, "union_linear_interpolation"


def _mf6_yields(
    mat: endf6.EndfMaterial, e_ev: NDArray[np.float64]
) -> dict[int, NDArray[np.float64]]:
    """MF6/MT5 yield of each product (by ZAP) at ``e_ev`` (clipped to each yield's range;
    products with the same ZAP are summed)."""
    out: dict[int, NDArray[np.float64]] = {}
    for prod in mat.products(5).products:
        yt = prod.yield_
        v = np.array([float(yt.interpolate(float(min(max(e, yt.x[0]), yt.x[-1])))[0])
                      for e in e_ev])  # fmt: skip
        out[prod.zap] = out.get(prod.zap, 0.0) + v
    return out


_PRODUCT_NAMES = {1: "n", 1001: "p", 1002: "d", 1003: "t", 2003: "he3", 2004: "alpha", 0: "gamma"}


def o16_mt5_copy_check(c12: endf6.EndfMaterial, o16: endf6.EndfMaterial) -> dict[str, Any]:
    """O-16 MF3/MT5 and MF6/MT5 are not a C-12 copy (Amendment 15 (d)1, both parts required).

    Compared energies: the shared MF3/MT5 nodes in ``[max threshold, 150 MeV]``, or (fewer than 10
    shared) the union grid with linear interpolation. ``o16_mt5_sigma_not_c12_copy`` (the sigma
    half of Amendment 15 (d)1 as frozen) is true only if ``|s_C - s_O| / max`` > 1e-3 at > 50 % of
    the energies AND the coefficient of variation of ``s_O / s_C`` > 1e-3 (a scaled copy has a
    constant ratio); it is the only flag set here. The MF6 per-product rule as frozen (all compared
    energies as denominator) is recorded as ``mf6_rule_as_frozen`` with its pass flag (it FAILS on
    LA150 for two sparse products and no longer qualifies anything, Amendment 16 item 1); the
    informative-node revision, written after that outcome was observed, is recorded as
    ``mf6_rule_revised`` and is report-only. A missing product counts as a zero yield. Parser or
    data errors propagate as :class:`BuildError` (never swallowed), and so does a failed check
    (fail closed)."""
    try:
        tc, to = c12.cross_section(5), o16.cross_section(5)
        e_ev, grid_kind = _mt5_compare_grid(tc, to)
        vc = np.array([float(tc.interpolate(float(e))[0]) for e in e_ev])
        vo = np.array([float(to.interpolate(float(e))[0]) for e in e_ev])
        yc, yo = _mf6_yields(c12, e_ev), _mf6_yields(o16, e_ev)
    except BuildError:
        raise
    except (endf6.EndfError, KeyError, ValueError) as exc:
        raise BuildError(f"O-16 / C-12 MT5 copy check cannot be evaluated: {exc!r}") from exc
    den = np.maximum(vc, vo)
    keep = den > 0.0
    rel = np.abs(vc[keep] - vo[keep]) / den[keep]
    frac = float(np.mean(rel > REL_DIFF_THRESHOLD))
    ratio = vo[keep] / np.where(vc[keep] > 0.0, vc[keep], np.nan)
    ratio = ratio[np.isfinite(ratio)]
    cv = float(np.std(ratio) / np.mean(ratio)) if ratio.size > 1 else 0.0
    sigma_ok = bool(frac > 0.5 and cv > REL_DIFF_THRESHOLD)
    frozen_prod: dict[str, Any] = {}
    revised_prod: dict[str, Any] = {}
    for zap in sorted(set(yc) | set(yo)):
        a = yc.get(zap, np.zeros(e_ev.size))
        b = yo.get(zap, np.zeros(e_ev.size))
        name = _PRODUCT_NAMES.get(zap, f"zap{zap}")
        d6 = np.maximum(a, b)
        if not np.any(d6 > 0.0):
            frozen_prod[name] = revised_prod[name] = None  # zero in both: nothing to compare
            continue
        # relative difference at every compared energy; both zero -> 0 (not differing), one zero
        # and the other nonzero -> 1.0
        inform = d6 > 0.0
        r6 = np.zeros(e_ev.size)
        r6[inform] = np.abs(a - b)[inform] / d6[inform]
        differs = r6 > REL_DIFF_THRESHOLD
        n_all, n_inf, n_diff = int(e_ev.size), int(np.sum(inform)), int(np.sum(differs))
        f_all, f_inf = n_diff / n_all, n_diff / n_inf
        frozen_prod[name] = {
            "zap": zap, "numerator": n_diff, "denominator": n_all, "fraction": f_all,
            "nodes_both_zero": n_all - n_inf, "differs_from_c12": bool(f_all > 0.5),
        }  # fmt: skip
        gating = n_inf >= MF6_MIN_INFORMATIVE
        revised_prod[name] = {
            "zap": zap, "informative_nodes": n_inf, "numerator": n_diff, "fraction": f_inf,
            "median_relative_difference": float(np.median(r6[inform])),
            "max_relative_difference": float(r6[inform].max()),
            "min_relative_difference": float(r6[inform].min()),
            "would_gate_under_revised_rule": gating,
            "differs_from_c12": bool(f_inf > 0.5),
        }  # fmt: skip
    frozen_fail = sorted(
        k for k, v in frozen_prod.items() if v is not None and not v["differs_from_c12"]
    )
    mf6_frozen = {
        "rule": "all compared energies are the denominator (both-zero nodes count as not "
        "differing); every nonzero product must exceed 50 %",
        "products": frozen_prod,
        "passes": not frozen_fail,
        "failing_products": frozen_fail,
    }
    mf6_revised = {
        "role": "report-only (post-outcome analysis, Amendment 16 item 1)",
        "rule": "denominator = informative compared energies (at least one material nonzero); a "
        f"product gates only with >= {MF6_MIN_INFORMATIVE} informative energies, else report-only; "
        "every gating product must exceed 50 %",
        "min_informative": MF6_MIN_INFORMATIVE, "products": revised_prod,
        "passes": all(v["differs_from_c12"] for v in revised_prod.values()
                      if v is not None and v["would_gate_under_revised_rule"]),
    }  # fmt: skip
    out: dict[str, Any] = {
        "compared_grid": grid_kind, "n_compared": int(e_ev.size),
        "compared_from_mev": float(e_ev[0] * 1e-6), "compared_to_mev": float(e_ev[-1] * 1e-6),
        "nodes_compared": int(rel.size), "fraction_rel_diff_gt_1e-3": frac,
        "max_relative_difference": float(rel.max()), "min_relative_difference": float(rel.min()),
        "median_relative_difference": float(np.median(rel)),
        "threshold": REL_DIFF_THRESHOLD, "required_fraction_gt": 0.5,
        "ratio_o16_over_c12_cv": cv, "required_cv_gt": REL_DIFF_THRESHOLD,
        "o16_mt5_sigma_not_c12_copy": sigma_ok,
        "mf6_rule_as_frozen": mf6_frozen, "mf6_rule_revised": mf6_revised,
    }  # fmt: skip
    if not sigma_ok:
        raise BuildError(f"LA150 O-16 MF3/MT5 sigma looks like a copy of C-12 (fail closed): {out}")
    return out


SPECTRA_PRODUCTS = (("n", 1), ("p", 1001), ("alpha", 2004))
SPECTRA_ENERGIES_MEV = (50.0, 100.0, 150.0)
SPECTRA_D_MIN = 0.02
SPECTRA_MIN_COMBINATIONS = 8
SPECTRA_ZERO_FLOOR = 1.0e-12
"""A KS distance below this is floating-point noise of the normalisation and counts as exactly 0
(a 1e-6 rescaling of a spectrum changes the normalised cumulative at the 1e-16 level)."""


def _spectrum_cdf(
    prod: endf6.Product, sp: SpeciesTables, e_mev: float
) -> tuple[NDArray[np.float64], NDArray[np.float64], str]:
    """Normalised angle-integrated E' cumulative distribution of one MF6/MT5 product at the
    incident energy ``e_mev``: ``(E' knots [MeV], cumulative at the knots, kind)``.

    LAW=1, LEP=1 (histogram pdf b_0 on [E'_j, E'_j+1)), exactly as the nuclear table builder reads
    it (:class:`SpeciesTables`). At a tabulated incident energy the cumulative is the exact
    piecewise-linear integral of the histogram; between nodes the builder's quantile
    interpolation of the 65 equiprobable edges is used (cumulative i/64 at edge i). A product
    outside the tabulated incident range fails closed."""
    en = sp.energies
    if not (en[0] <= e_mev <= en[-1]):
        raise BuildError(f"incident energy {e_mev} MeV outside the MF6 range of ZAP={sp.zap}")
    node = np.nonzero(np.abs(en - e_mev) <= 1e-9 * e_mev)[0]
    if node.size:
        d = prod.distributions[int(node[0])]
        x = d.rows[:, 0] * 1.0e-6
        w = d.rows[:-1, 1] * np.diff(x)
        total = float(w.sum())
        if not total > 0.0:
            raise BuildError(f"ZAP={sp.zap} has no positive spectrum weight at {e_mev} MeV")
        return x, np.concatenate(([0.0], np.cumsum(w) / total)), "tabulated_node_histogram"
    edges = sp.rows_at(e_mev)[0]
    return edges, np.arange(edges.size) / (edges.size - 1.0), "quantile_interpolation_between_nodes"


def o16_mt5_spectra_check(c12: endf6.EndfMaterial, o16: endf6.EndfMaterial) -> dict[str, Any]:
    """O-16 MF6/MT5 product energy spectra are not a C-12 copy (Amendment 16 item 2, X-MT5-SPEC).

    For n, p, alpha at 50, 100 and 150 MeV the normalised angle-integrated secondary-energy
    distributions of both nuclei are compared by the Kolmogorov-Smirnov distance D = max |F_O -
    F_C| on the union E' grid. ``o16_mt5_spectra_independent`` is true iff D > 0.02 at >= 8 of the
    9 combinations AND D > 0 at all 9 (D below :data:`SPECTRA_ZERO_FLOOR` counts as 0). Fail closed
    (:class:`BuildError`) if a product or incident energy cannot be reconstructed or the flag is
    false."""
    try:
        secs = {"C-12": c12.products(5), "O-16": o16.products(5)}
        tabs: dict[str, dict[int, tuple[endf6.Product, SpeciesTables]]] = {}
        for name, sec in secs.items():
            by_zap = {p.zap: p for p in sec.products}
            tabs[name] = {}
            for pname, zap in SPECTRA_PRODUCTS:
                prod = by_zap.get(zap)
                if prod is None or prod.law != 1:
                    raise BuildError(f"{name}: no MF6/MT5 LAW=1 product {pname} (ZAP={zap})")
                tabs[name][zap] = (prod, SpeciesTables(prod))
        rows: list[dict[str, Any]] = []
        for pname, zap in SPECTRA_PRODUCTS:
            for e in SPECTRA_ENERGIES_MEV:
                xc, fc, kc = _spectrum_cdf(*tabs["C-12"][zap], e)
                xo, fo, ko = _spectrum_cdf(*tabs["O-16"][zap], e)
                grid = np.union1d(xc, xo)
                d = float(np.max(np.abs(np.interp(grid, xo, fo) - np.interp(grid, xc, fc))))
                rows.append({
                    "product": pname, "zap": zap, "incident_mev": e, "D": d,
                    "D_effective": 0.0 if d < SPECTRA_ZERO_FLOOR else d,
                    "n_union_grid": int(grid.size),
                    "grid_from_mev": float(grid[0]), "grid_to_mev": float(grid[-1]),
                    "c12_reconstruction": kc, "o16_reconstruction": ko,
                })  # fmt: skip
    except BuildError:
        raise
    except (endf6.EndfError, KeyError, ValueError, AttributeError) as exc:
        raise BuildError(f"O-16 / C-12 MT5 spectra cannot be evaluated: {exc!r}") from exc
    d_eff = np.array([r["D_effective"] for r in rows])
    n_above = int(np.sum(d_eff > SPECTRA_D_MIN))
    independent = bool(n_above >= SPECTRA_MIN_COMBINATIONS and np.all(d_eff > 0.0))
    out: dict[str, Any] = {
        "law": "MF6/MT5 LAW=1 (continuum energy-angle), LEP=1 histogram, angle-integrated (b_0), "
        "read by SpeciesTables as in the nuclear table builder; between incident nodes the "
        "builder's quantile interpolation is used",
        "statistic": "Kolmogorov-Smirnov distance on the union E' grid",
        "rule": f"D > {SPECTRA_D_MIN} at >= {SPECTRA_MIN_COMBINATIONS} of 9 combinations and "
        "D > 0 at all 9",
        "zero_floor": SPECTRA_ZERO_FLOOR, "combinations": rows,
        "D_values": {f"{r['product']}@{r['incident_mev']:g}": r["D"] for r in rows},
        "n_above_threshold": n_above, "n_zero": int(np.sum(d_eff == 0.0)),
        "o16_mt5_spectra_independent": independent,
    }  # fmt: skip
    if not independent:
        raise BuildError(
            f"LA150 O-16 MF6/MT5 product spectra look like a copy of C-12 (fail closed): {out}"
        )
    return out


def check_domain_limit(name: str, e_min_shape_mev: float) -> None:
    """``BuildError`` if a p + A target's ``e_min_shape`` exceeds :data:`E_MIN_SHAPE_LIMIT_MEV`."""
    if e_min_shape_mev > E_MIN_SHAPE_LIMIT_MEV:
        raise BuildError(
            f"{name}: e_min_shape = {e_min_shape_mev:.4g} MeV exceeds the declared limitation "
            f"bound {E_MIN_SHAPE_LIMIT_MEV} MeV"
        )


def residual_range_g_cm2(e_mev: float) -> float:
    """CSDA range [g/cm2] of a proton of ``e_mev`` in liquid water (project Bethe table)."""
    return float(_water_tables().range_g_cm2(0, e_mev))


def _water_tables() -> Any:
    from ionmc.materials import WATER
    from ionmc.physics.projectiles import PROTON
    from ionmc.physics.stopping import BetheStoppingSource
    from ionmc.transport.tables import TransportTables

    return TransportTables.from_stopping_tables([BetheStoppingSource().table(WATER, PROTON)])


def omitted_events_per_history(
    e_cut_mev: float, n_cm3: float, sigma_mb: Callable[[float], float], e_lo_mev: float = 1.0
) -> float:
    """Expected number of omitted elastic events per history of a proton that slows down through
    ``[e_lo, e_cut]`` in liquid water: ``n int sigma(E) dx`` with ``dx = dE / S(E)`` (project Bethe
    stopping table, density 1 g/cm3). ``sigma`` is never evaluated or held constant below ``e_lo``
    (1 MeV, the lowest node of the table); the residual range ``R(e_lo)`` is reported separately
    by the caller. Nuclear attrition is ignored (<< 1)."""
    tab = _water_tables()
    e = np.linspace(e_lo_mev, e_cut_mev, 4001)
    f = np.array([sigma_mb(float(x)) * 1e-27 / tab.stopping_mass(0, float(x)) for x in e])
    integral = float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(e)))
    return float(n_cm3 * integral)


PP_OMIT_SEMANTICS = (
    "total variation of the signed nuclear-plus-interference correction to Rutherford scattering "
    "over the half sphere above the cut, integrated along the residual path; an event-equivalent "
    "magnitude, not an expected number of physical events; the weighted form is the "
    "recoil-energy-weighted total variation, not energy transferred; limitations: the signed "
    "correction can cancel, the physical sigma_pp below E_min,pp is not available from the "
    "evaluation"
)


def pp_omitted_ni_correction(
    sec: law5.Law5Section,
    awi: float,
    grid: NDArray[np.float64],
    k_pp: int,
    mu_h: NDArray[np.float64],
    e_lo_mev: float = 1.0,
) -> dict[str, Any]:
    """DATA-BASED total-variation diagnostic of the omitted p-p nuclear-plus-interference (NI)
    correction to Coulomb scattering below E_min,pp (Amendment 15 (a)2; see ``PP_OMIT_SEMANTICS``;
    not a model bound and not an event count).

    At each grid node ``E`` in ``[e_lo, E_min,pp]`` the Hale LAW=5 reconstruction gives
    ``M(E) = 2 pi int_0^{0.96} |rho_NI(mu, E)| dmu`` (absolute value, so cancellations do not hide
    magnitude; half sphere, each event once) and the energy-weighted ``W(E) = 2 pi int |rho_NI|
    T (1 - mu)/2 dmu`` (kinetic energy handed to the recoil at CM angle theta_CM, equal masses).
    Then the total variation ``N = n_H int M (dx/dE) dE`` and ``E_omit = n_H int W (dx/dE) dE``
    over the path of a proton slowing down from E_min,pp to ``e_lo`` (``dx/dE = 1/S``, project
    water Bethe table, trapezoid on the grid nodes). Below ``e_lo`` M is not extrapolated: the
    residual range ``R(e_lo)`` is reported separately."""
    tab = _water_tables()
    nodes = grid[(grid >= e_lo_mev) & (grid <= float(grid[k_pp]))]
    m = np.zeros(nodes.size)
    w = np.zeros(nodes.size)
    for i, e in enumerate(nodes):
        d = np.abs(law5.ni_density_ltp1(sec, float(e) * 1e6, mu_h, awi, 1, 1))
        m[i] = 2.0 * math.pi * law5.simpson(d, mu_h)  # barn
        w[i] = 2.0 * math.pi * law5.simpson(d * float(e) * (1.0 - mu_h) / 2.0, mu_h)  # barn MeV
    inv_s = np.array([1.0 / tab.stopping_mass(0, float(e)) for e in nodes])  # g/cm2/MeV
    de = np.diff(nodes)
    n_h = N_H_WATER_CM3

    def trap(f: NDArray[np.float64]) -> float:
        return float(np.sum(0.5 * (f[1:] + f[:-1]) * de))

    return {
        "total_variation_per_history": n_h * 1e-24 * trap(m * inv_s),
        "weighted_total_variation_mev": n_h * 1e-24 * trap(w * inv_s),
        "nodes_mev": nodes.tolist(), "m_barn": m.tolist(), "w_barn_mev": w.tolist(),
        "max_mean_recoil_fraction_of_t": float(np.max(w / (nodes * np.maximum(m, 1e-300)))),
        "residual_range_below_1mev_g_cm2": float(tab.range_g_cm2(0, e_lo_mev)),
        "semantics": PP_OMIT_SEMANTICS,
        "n_h_cm3": n_h, "e_cut_mev": float(grid[k_pp]), "e_lo_mev": e_lo_mev,
        "method": (
            "data-based, not a model estimate: below E_min,pp the transport omits the nuclear-plus-"
            "interference correction to Coulomb scattering. From the Hale LAW=5 LTP=1 "
            "reconstruction at every grid node of [1 MeV, E_min,pp]: M(E) = 2 pi int_0^0.96 "
            "|rho_NI| dmu (absolute value; half sphere, each event once); N = n_H int M (dx/dE) dE "
            "over the residual path of a proton from E_min,pp to 1 MeV (project water Bethe "
            "table, trapezoid on the grid nodes, n_H = 6.69e22 cm^-3); the energy is the same "
            "integral with M replaced by 2 pi int |rho_NI| T (1-mu)/2 dmu. M is not extrapolated "
            "below 1 MeV: the residual range there is reported separately."
        ),
    }  # fmt: skip


def reject_endf_mt2_construction(target: str) -> None:
    """C1-ext: the LA150 LAW=5 data of any target other than H-1 may not be a construction input
    (O-16 is a C-12 copy; C-12/N-14/Ca-40 are report-only cross-checks)."""
    if target != PP_TARGET:
        raise UnsupportedCombinationError(
            f"LA150 MT2 of {target} cannot be a construction input (decision 0041: O-16 is a "
            "copy of C-12, nuclear-plus-interference double counts the Coulomb scattering)"
        )


def ltp12_ni_density_ratio(
    sec: law5.Law5Section,
    index: int,
    sigma_mf3_b: float,
    mu: NDArray[np.float64],
    awi: float,
    z: int,
) -> NDArray[np.float64]:
    """X-ENDF: NI density [b/sr] of an LTP=12 node with R(mu) = 1 + sigma_MF3 P(mu)/sigma_c
    interpolated linearly in mu (ratio, never P_NI) and NI = sigma_c (R - 1)."""
    rec = sec.records[index]
    k, eta = law5.coulomb_params(rec.energy_ev, sec.awr, awi, 1, z)
    mu_t, p = rec.data[0::2], rec.data[1::2]
    ratio_t = sigma_mf3_b * p / law5.sigma_c(mu_t, k, eta, 0, 0.0)
    ratio = np.interp(mu, mu_t, ratio_t)
    return np.asarray(law5.sigma_c(mu, k, eta, 0, 0.0) * ratio)


def _noop(_: str) -> None:
    return None


# ---------------------------------------------------------------------------------------------
# grid
# ---------------------------------------------------------------------------------------------
def barashenkov_nodes_mev() -> NDArray[np.float64]:
    e = np.unique(
        np.concatenate([np.asarray(v, dtype=np.float64) * bgg.GEV_MEV for v in
                        (bgg._E1_GEV, bgg._E2_GEV, bgg._E3_GEV)])
    )  # fmt: skip
    return e[(e >= 1.0) & (e <= E_MAX_MEV)]


def _merge(grid: NDArray[np.float64], extra: NDArray[np.float64]) -> NDArray[np.float64]:
    allv = np.unique(np.concatenate((grid, extra)))
    keep = np.ones(allv.size, dtype=bool)
    keep[1:] = np.diff(allv) > 1e-12 * allv[1:]
    return allv[keep]


def build_grid_elastic(
    h1_nodes_mev: NDArray[np.float64], points_per_decade: int
) -> NDArray[np.float64]:
    """Union grid (module docstring): nuclear convention + Barashenkov nodes, cut at 250 MeV."""
    bar = barashenkov_nodes_mev()
    grid, _ = union_grid(
        np.concatenate((h1_nodes_mev, bar[bar <= E_ANCHOR_MEV])), points_per_decade
    )
    grid = _merge(grid[grid < E_MAX_MEV], bar)
    grid = grid[grid <= E_MAX_MEV]
    if grid[-1] != E_MAX_MEV:
        raise BuildError("the grid does not end at 250 MeV")
    return grid


def refine_grid(
    grid: NDArray[np.float64], f: Callable[[float], list[float]], lo: float, hi: float
) -> tuple[NDArray[np.float64], int]:
    """Insert midpoints in the intervals inside ``[lo, hi]`` MeV until the lin-lin interpolation of
    every component of ``f(E)`` reproduces the direct value at the midpoint to ``REFINE_TOL``
    (relative to the larger of the value and 1e-3 of ``f(hi)``). Returns ``(grid, nodes_added)``."""
    ref = np.maximum(np.asarray(f(hi)), 1e-30)
    n0 = grid.size
    for _ in range(MAX_REFINE_ROUNDS):
        mids = []
        for a, b in zip(grid[:-1], grid[1:], strict=True):
            if a < lo or b > hi:
                continue
            m = 0.5 * (a + b)
            fm, fa, fb = np.asarray(f(m)), np.asarray(f(float(a))), np.asarray(f(float(b)))
            err = np.abs(0.5 * (fa + fb) - fm) / np.maximum(fm, 1e-3 * ref)
            if np.any(err > REFINE_TOL):
                mids.append(m)
        if not mids:
            break
        grid = _merge(grid, np.asarray(mids))
    return grid, int(grid.size - n0)


H1_REFINE_FROM_MEV = 13.0
"""Lower bound of the H-1 sigma_NI refinement (just above the first node with a non-negative
density; the positivity itself is decided independently by ``e_min_pp``)."""


def disk_density_mb_sr(
    sigma_el_b: float, r_fm: float, pcm_mev: float, theta_deg: NDArray[np.float64], n_s: int = 65537
) -> NDArray[np.float64]:
    """Model ``dsigma/dOmega_CM`` [mb/sr] of the black-disk shape normalised to ``sigma_el_b``:
    ``sigma_el A(theta) / (2 pi int_{-1}^{1} A dmu)`` with ``A = |2 J1(qR)/(qR)|^2``."""
    s = np.linspace(0.0, 1.0, n_s)
    g = 4.0 * s * disk_amplitude_sq(2.0 * pcm_mev / HBARC * s * r_fm)
    norm = float(np.sum(0.5 * (g[1:] + g[:-1]) * np.diff(s)))  # int A dmu over [-1, 1]
    sin_half = np.sin(np.radians(theta_deg) / 2.0)
    a = disk_amplitude_sq(2.0 * pcm_mev / HBARC * sin_half * r_fm)
    return np.asarray(1.0e3 * sigma_el_b * a / (2.0 * math.pi * norm))


def _shape_norm_deviation(r_fm: float, pcm_mev: float, n_s: int) -> float:
    """Relative difference of the shape integral at ``n_s`` and at 4 (n_s - 1) + 1 nodes."""
    out = []
    for n in (n_s, 4 * (n_s - 1) + 1):
        s = np.linspace(0.0, 1.0, n)
        g = 4.0 * s * disk_amplitude_sq(2.0 * pcm_mev / HBARC * s * r_fm)
        out.append(float(np.sum(0.5 * (g[1:] + g[:-1]) * np.diff(s))))
    return abs(out[0] / out[1] - 1.0)


def pp_scan(
    grid: NDArray[np.float64],
    min_node: NDArray[np.float64],
    min_mid: NDArray[np.float64],
    i150: int,
) -> tuple[int, list[dict[str, Any]]]:
    """Pass 1 of P6/P6-D: ``(k_pp, negative points)``. ``k_pp`` is the lowest node from which every
    node and node midpoint up to ``grid[i150]`` has a non-negative density; the points list holds
    EVERY negative node/midpoint (the recorded failure of the frozen P6 item (4)). ``BuildError``
    if ``E_min,pp`` exceeds ``i150`` or :data:`E_MIN_PP_LIMIT_MEV`."""
    ok = (min_node[: i150 + 1] >= 0.0) & (np.append(min_mid[:i150], np.inf) >= 0.0)
    bad = np.nonzero(~ok)[0]
    k_pp = int(bad.max()) + 1 if bad.size else 0
    if k_pp > i150:
        raise BuildError("no energy up to 150 MeV has a non-negative p-p density")
    if float(grid[k_pp]) > E_MIN_PP_LIMIT_MEV:
        raise BuildError(
            f"E_min,pp = {float(grid[k_pp]):.4g} MeV exceeds the bound {E_MIN_PP_LIMIT_MEV} MeV"
        )
    pts: list[dict[str, Any]] = [
        {"kind": "node", "e_mev": float(grid[k]), "min_density_mb_sr": 1e3 * float(min_node[k])}
        for k in range(i150 + 1) if min_node[k] < 0.0
    ] + [
        {"kind": "midpoint", "e_mev": 0.5 * float(grid[k] + grid[k + 1]),
         "min_density_mb_sr": 1e3 * float(min_mid[k])}
        for k in range(i150) if min_mid[k] < 0.0
    ]  # fmt: skip
    pts.sort(key=lambda r: float(r["e_mev"]))
    return k_pp, pts


def verify_pp_domain(
    grid: NDArray[np.float64], k_pp: int, i150: int, check: Callable[[float], float]
) -> dict[str, Any]:
    """Pass 2 of P6-D (1), separate from the pass that fixed ``k_pp`` (the density is re-evaluated
    by ``check(e_mev) -> min density over |mu| <= 0.96`` on a finer mu grid): every node and node
    midpoint in ``[E_min,pp, 150]`` MeV must be non-negative, else ``BuildError``."""
    pts = [float(grid[k]) for k in range(k_pp, i150 + 1)]
    pts += [0.5 * float(grid[k] + grid[k + 1]) for k in range(k_pp, i150)]
    worst = min(float(check(e)) for e in pts)
    if worst < 0.0:
        raise BuildError(
            f"P6-D: a negative p-p density ({1e3 * worst:.3g} mb/sr) at or above E_min,pp = "
            f"{float(grid[k_pp]):.4g} MeV"
        )
    return {"n_points_checked": len(pts), "min_density_mb_sr": 1e3 * worst}


def _nan_to_none(x: Any) -> Any:
    if isinstance(x, float) and not math.isfinite(x):
        return None
    if isinstance(x, dict):
        return {k: _nan_to_none(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_nan_to_none(v) for v in x]
    return x


def _pp_reanchor(
    m2: glauber.M2Config, model: glauber.M2Model, grid: NDArray[np.float64],
    sigma: NDArray[np.float64], sig_h: NDArray[np.float64], i150: int, k_pp: int,
    s_ref: float, awi: float,
) -> dict[str, Any]:  # fmt: skip
    """Amendment 20 (b')4 (revision 3): SAID DSG - Mott anchor of the H-1 level; modifies
    ``sigma[0]`` in place when the blocking checks pass; returns the record (old/new level ratio at
    every node >= 100 MeV, the scan-energy rows, the checks)."""
    assert m2.dsg_files is not None
    sgt = model.nn.series[("PP", "SGT")]
    e_lo = grid[k_pp : i150 + 1]
    rows = []
    for t in sorted(m2.dsg_files):
        a, d = glauber.merge_dsg_segments(
            [Path(f).read_bytes().decode("ascii", "replace") for f in m2.dsg_files[t]]
        )
        row = glauber.dsg_anchor(t, a, d, sgt.at(t), awi)
        row["files_sha256"] = [glauber._sha256(Path(f)) for f in m2.dsg_files[t]]
        model.check_dsg_pinned(row["files_sha256"])
        row["sigma_table_mb"] = (
            1e3 * float(np.interp(t, e_lo, sigma[0, k_pp : i150 + 1]))
            if t <= E_ANCHOR_MEV
            else 1e3 * sig_h[i150] * bgg.pp_bgg_mb(t) / s_ref
        )
        row["ratio_to_table"] = row.get("sigma_anchor_mb", float("nan")) / row["sigma_table_mb"]
        rows.append(row)
    rea = glauber.PPReanchor.build(rows)
    old = 1e3 * sigma[0].copy()
    if rea.enabled:
        for k in range(k_pp, grid.size):
            if grid[k] > glauber.RAMP_FROM_MEV:
                sigma[0, k] = 1e-3 * rea.sigma_mb(float(grid[k]), old[k])
    nodes = [
        {"e_mev": float(grid[k]), "sigma_old_mb": float(old[k]),
         "sigma_new_mb": float(1e3 * sigma[0, k]), "ratio_new_over_old":
         float(sigma[0, k] * 1e3 / old[k])}
        for k in range(k_pp, grid.size) if grid[k] >= 100.0
    ]  # fmt: skip
    return {"enabled": rea.enabled, "reason": rea.reason,
            "pp_level_reanchoring": (
                "applied" if rea.enabled else f"not applied ({rea.reason})"
            ),
            "scan_rows": rows, "nodes": nodes,
            "checks": {"continuity_tol": glauber.CONTINUITY_TOL,
                       "sgt_ratio_band": list(glauber.SGT_RATIO_BAND)},
            "mott": "ENDF-102 eq 6.10 (law5.sigma_c), identical spin-1/2, half sphere 2 pi "
            "int_0^0.96 [DSG - Mott] dmu; ln(DSG - Mott) linear in theta_CM"}  # fmt: skip


def _pauli_record(
    m2: glauber.M2Config, nodes: dict[str, list[dict[str, Any]]],
    pinned: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:  # fmt: skip
    """M2' block: the formula pinned in :func:`glauber.pauli_factor` and the Pauli factor at the
    peak density per target at 100/150/200 MeV."""
    peak = {
        n: {f"{e:g}_mev": min(v, key=lambda r: abs(r["e_mev"] - e))["pauli_factor_peak"]
            for e in (100.0, 150.0, 200.0)}
        for n, v in nodes.items() if v
    }  # fmt: skip
    return {"pauli": {
        "reference": "Clementel & Villi 1955 / Hussein, Rego & Bertulani 1991 (variant of record); "
        "Bertulani & De Conti, PRC 81, 064603 (2010), arXiv:1004.2096 (bdc-pn)",
        "rho0_fm3": glauber.RHO0_FM3, "pauli_scale": m2.pauli_scale,
        "pauli_emin_mev": m2.pauli_emin_mev,
        "nodes_applied_per_target": {n: len(v) for n, v in nodes.items()},
        "nodes_below_emin_black_disk_per_target": {n: len(v) for n, v in pinned.items()},
        "nodes_below_emin_black_disk": pinned,
        "pauli_variant": m2.pauli_variant,
        "reading": {"clementel-villi": "Fermi-gas LDA, E_F(rho_tot), variant of record",
                    "bdc-pn": "Eq. 14 with rho1 = local rho_p, rho2 = local rho_n; "
                    "report-only sensitivity"}.get(str(m2.pauli_variant)),
        "formula": "see glauber.pauli_factor: clementel-villi 1 - 7E_F/(5E) (E >= 2E_F), "
        "low branch 1 - 7E_F/(5E) + (2/5)(E_F/E)(2 - E/E_F)^(5/2), E_F = (hbar^2/2m)(3 pi^2 "
        "rho_tot/2)^(2/3); bdc-pn: Bertulani & De Conti Eq. 14 (arXiv:1004.2096)",
        "alpha": "free SAID", "beta": "unitarity rule of the free sigma",
        "pauli_factor_at_peak_density": peak}}  # fmt: skip


def _trigger2_record(nodes: dict[str, list[dict[str, Any]]], conv: float) -> dict[str, Any]:
    """Plan 4.2 trigger 2 / Amendment 21: step-halving convergence, every component separately
    (b step, s grid, radial step, q step, q cutoff), the maximum per component over all evaluated
    nodes, and the rows at the candidate energies (nodes nearest 100/150/200 MeV)."""
    comps = ("b_step", "s_grid", "radial_step", "q_step", "q_cutoff")
    recs = [(n, r["convergence"]) for n, v in nodes.items() for r in v if "convergence" in r]
    cmax = {
        c: max((cv["components"][c] for _, cv in recs if cv["components"][c] is not None),
               default=None)
        for c in comps
    }  # fmt: skip
    cand = [{"target": n, "e_mev": cv["e_mev"], "components": cv["components"], "max": cv["max"]}
            for n, cv in recs if cv.get("candidate_energy")]  # fmt: skip
    return {"tol": glauber.CONVERGENCE_TOL, "max": conv, "passed": conv <= glauber.CONVERGENCE_TOL,
            "component_max": cmax, "n_evaluated": len(recs), "candidate_energies": cand,
            "note": "radial/q components are None for analytic HO thickness"}  # fmt: skip


def m2_record(
    m2: glauber.M2Config, model: glauber.M2Model, nodes: dict[str, list[dict[str, Any]]],
    conv: float, norm: float, pp_rec: dict[str, Any] | None,
    pinned: dict[str, list[dict[str, Any]]] | None = None,
) -> dict[str, Any]:  # fmt: skip
    """Provenance and construction checks of an M2 table (plan 4.2), recorded and asserted."""
    if norm > glauber.NORM_TOL:
        raise BuildError(f"M2 normalisation deviation {norm:.3g} exceeds {glauber.NORM_TOL}")
    if conv > glauber.CONVERGENCE_TOL:
        raise BuildError(f"M2 step-halving convergence {conv:.3g} exceeds 1e-3 (plan 4.2 trig. 2)")
    bd = glauber.black_disk_limit_check()
    if bd["max_rel_dev"] > bd["tolerance"]:
        raise BuildError(f"M2 black-disk limit check failed: {bd}")
    worst_unit = max(r["unitarity_residual"] for v in nodes.values() for r in v)
    worst_ot = max(r["optical_theorem_residual"] for v in nodes.values() for r in v)
    min_dens = min(r["min_density_fm2_sr"] for v in nodes.values() for r in v)
    if min_dens < 0.0:
        raise BuildError("M2 density is negative")
    level = {}
    for name, v in nodes.items():
        level[name] = {
            f"{e:g}_mev": float(
                np.interp(e, [r["e_mev"] for r in v], [r["m2_over_bgg"] for r in v])
            )
            for e in (100.0, 150.0, 200.0)
        }
    t3 = {n: [{"e_mev": float(e.split("_")[0]), "ratio": ratio,
               "out_of_band": not 0.7 <= ratio <= 1.4} for e, ratio in level[n].items()]
          for n in ("C-12", "O-16")}  # fmt: skip
    trigger3 = {"band": [0.7, 1.4], "targets": t3,
                "fired": any(sum(r["out_of_band"] for r in v) >= 2 for v in t3.values()),
                "note": "plan 4.2 trigger 3: >= 2 of 100/150/200 MeV outside the band for O or C; "
                "recorded here, NOT acted upon by the builder (orchestrator decision)"}  # fmt: skip
    trigger4: dict[str, Any] = {"band": [0.6, 1.0]}
    for n in ("C-12", "O-16"):
        dn = model.dens[n]
        assert isinstance(dn, glauber.HODensity)  # C-12, O-16 are HO
        qh = glauber.helm_q_h(dn)
        rows = {}
        for e in (150.0, 200.0):
            q = min(nodes[n], key=lambda r: abs(r["e_mev"] - e))["q_first_min_fm_inv"]
            ok = isinstance(q, float) and math.isfinite(q)
            rows[f"{e:g}_mev"] = {
                "q_min_m2_fm_inv": q if ok else "not evaluable (no true minimum)",
                "ratio_to_q_h": q / qh["q_h_fm_inv"] if ok else None,
                "in_band": (0.6 <= q / qh["q_h_fm_inv"] <= 1.0) if ok else None,
            }
        trigger4[n] = {"helm": qh, "rows": rows}  # fmt: skip
    return {
        "level_sanity_trigger_3": trigger3, "radius_sanity_trigger_4": trigger4,
        "unfold_note": "rms-matched is the declared rule for 2pF/3pF; 'formfactor' diverges when "
        "beta < r_p^2/3 = 0.236 fm^2 (E > ~125 MeV)",
        "variant": m2.variant(), "beta": m2.beta, "pp_level": m2.pp_level, "unfold": m2.unfold,
        "options": {"n_s": m2.n_s, "db_fm": m2.db_fm, "q_cut": m2.q_cut,
                    "shape_override": m2.shape_override},
        "provenance": model.provenance(),
        "trigger_2_convergence": _trigger2_record(nodes, conv),
        "checks": {"normalisation_max_deviation": norm, "normalisation_tol": glauber.NORM_TOL,
                   "convergence_max": conv, "convergence_tol": glauber.CONVERGENCE_TOL,
                   "unitarity_max_residual": worst_unit, "optical_theorem_max_residual": worst_ot,
                   "min_density_fm2_sr": min_dens, "black_disk_limit": bd,
                   "rms_radius": {n: model.rms[n] for n in model.rms},
                   "m2_over_bgg_sigma_el_at_100_150_200": level,
                   "unfolding_tail_max": {n: model.thick[n].tail for n in model.thick}},
        "first_minimum_q_150_mev": {
            n: min(v, key=lambda r: abs(r["e_mev"] - 150.0))["q_first_min_fm_inv"]
            for n, v in nodes.items()
        },
        "targets": {n: {"nodes": v} for n, v in nodes.items()},
        "pp_reanchor": pp_rec,
        "note": "Coulomb stays in the MCS model; nuclear-only amplitude; no Coulomb-nuclear "
        "interference (plan 4.1 item 5); sigma_el = pinned BGG",
        **(_pauli_record(m2, nodes, pinned or {}) if m2.pauli else {}),
        **(
            {"sigma_tot_comparison_report_only": glauber.sigma_tot_comparison(
                glauber.load_measured_sigma_tot(m2.measured_sigma_tot), nodes)}
            if m2.measured_sigma_tot is not None else {}
        ),
    }  # fmt: skip


def build_elastic_proton(
    cache_dir: str | Path | None = None,
    options: ElasticBuildOptions | None = None,
    log: Callable[[str], None] = _noop,
    *,
    exploratory_forced_isotropic: bool = False,
    m2: glauber.M2Config | None = None,
    m3: m3mod.M3Config | None = None,
    out_dir: Path | None = None,
) -> ElasticBuildResult:
    """Build the table (module docstring). ``BuildError`` on any fail-closed rule.

    ``exploratory_forced_isotropic`` (E-shape only, plan Amendment 15 (b)6; never evidence): the
    nodes of every p + A target BELOW its declared domain ``e_min_shape`` carry the BGG
    ``sigma_el`` (valid) with an isotropic ``mu_CM`` distribution (edges uniform in [-1, 1]);
    every node inside the domain, every ``sigma_el`` there, the grid and H-1 are those of the
    model table. The table JSON carries ``exploratory`` (so its id differs and the loader
    refuses it unless the configuration names it and sets ``allow_exploratory_table``).

    ``m2`` (V3-005D, exploratory, never evidence): the p + A angular shapes are the optical-limit
    Glauber shapes of :mod:`ionmc.nuclear.elastic_glauber` (``sigma_el`` stays the pinned BGG
    value; the domains ``e_min_shape``, the grid and the H-1 shape are those of the model table);
    ``m2.pp_level = "said-anchored-ge150"`` additionally re-anchors the H-1 level at >= 125 MeV
    (Amendment 20 (b')4). The table carries ``exploratory`` with the variant name and ``m2`` with
    the complete input provenance and the construction checks. ``m3`` (V3-005D fallback model,
    exploratory, never evidence): the C-12 and O-16 shapes are the ENDF/B-VIII.0 neutron MF4/MT2
    angular distributions at the same lab energy (:mod:`ionmc.nuclear.elastic_m3`); all other
    targets keep the black-disk shape (recorded); the evaluations end at 150 MeV (verified
    2026-10-23) and ``m3.above_emax`` ("error" | "hold-last" | "black-disk") decides the nodes
    above it.
    ``m2`` and ``m3`` are mutually exclusive. ``out_dir`` overrides the output
    directory (default ``<cache>/derived``); scratch builds go to a git-ignored directory."""
    t_start = time.perf_counter()
    opt = options or ElasticBuildOptions()
    cdir = cache.resolve_cache_dir(cache_dir)
    zip_path = cache.verify("endf-b8.0-protons", cdir)
    ame_tab = load_ame2020(cache.verify("ame2020-mass", cdir).read_text(encoding="ascii"))
    for sid in SOURCE_IDS:
        cache.verify(sid, cdir)
    transcription = verify_transcription(cdir)
    n_q = opt.n_quantiles

    def member(stem: str) -> endf6.EndfMaterial:
        return endf6.parse_endf(endf6.read_member(zip_path, f"ENDF-B-VIII.0_protons/{stem}.endf"))

    h1 = member("p-001_H_001")
    sec, awi = law5.parse_law5(h1.sections[(6, 2)]), law5.projectile_awi(h1)
    if sec.lidp != 1 or {r.ltp for r in sec.records} != {1}:
        raise BuildError("H-1 MF6/MT2 is not LTP=1 LIDP=1")
    mats = {t.name: member(t.member) for t in TARGETS}
    models = {t.name: ev.build_event_model(ame_tab, t.z, t.a) for t in TARGETS}
    tts = {t.name: TargetTables(t, mats[t.name], models[t.name]) for t in TARGETS}
    m2model = glauber.M2Model(m2, [(t.name, t.z, t.a) for t in TARGETS]) if m2 is not None else None
    if m2 is not None and m3 is not None:
        raise BuildError("m2 and m3 are mutually exclusive")
    m3model = m3mod.M3Model(m3) if m3 is not None else None
    m3_nodes: dict[str, list[dict[str, Any]]] = {t.name: [] for t in TARGETS}
    m3_above = {t.name: 0 for t in TARGETS}
    m3_mid = {t.name: 0 for t in TARGETS}
    m3_norm = 0.0
    m2_nodes: dict[str, list[dict[str, Any]]] = {t.name: [] for t in TARGETS}
    m2_conv = 0.0
    m2_below = {t.name: 0 for t in TARGETS}
    m2_pinned: dict[str, list[dict[str, Any]]] = {t.name: [] for t in TARGETS}
    m2_norm = 0.0
    o16 = o16_copy_finding(mats["C-12"], mats["O-16"])
    log(f"O-16 MT2 is a C-12 copy from {o16['from_mev']} MeV: recorded")
    o16_mt5 = o16_mt5_copy_check(mats["C-12"], mats["O-16"])
    log(f"O-16 MT5 differs from C-12 at {o16_mt5['fraction_rel_diff_gt_1e-3']:.0%} of the nodes")
    o16_spec = o16_mt5_spectra_check(mats["C-12"], mats["O-16"])
    log(f"O-16 MT5 spectra KS D: {o16_spec['D_values']}")

    e_h1 = sec.energies_ev * 1e-6
    grid = build_grid_elastic(e_h1[(e_h1 >= 1.0) & (e_h1 <= E_ANCHOR_MEV)], opt.points_per_decade)

    def direct_mb(e: float) -> list[float]:
        return [bgg.bgg_elastic_mb(e, t.z, models[t.name].m_t_mev) for t in TARGETS]

    grid, n_refined = refine_grid(grid, direct_mb, 0.0, bgg.BGG_LOW_ENERGY_MEV)
    mu_r = np.linspace(0.0, law5.MU_CUT_PP, opt.n_mu_check)

    def sigma_pp_nb(e: float) -> list[float]:
        d = law5.ni_density_ltp1(sec, e * 1e6, mu_r, awi, 1, 1)
        return [2.0 * math.pi * law5.simpson(d, mu_r)]

    grid, n_refined_pp = refine_grid(grid, sigma_pp_nb, H1_REFINE_FROM_MEV, E_ANCHOR_MEV)
    n_grid, n_t = grid.size, len(TARGET_NAMES)
    cand_k = {int(np.argmin(np.abs(grid - e))) for e in (100.0, 150.0, 200.0)}
    i150 = int(np.searchsorted(grid, E_ANCHOR_MEV))
    log(f"grid {n_grid} nodes ({n_refined} below 14 MeV, {n_refined_pp} for H-1)")

    sigma = np.zeros((n_t, n_grid))
    edges = np.zeros((n_t, n_grid, n_q + 1))
    r_fm = np.zeros((n_t, n_grid))
    sig_nonel = np.zeros((n_t, n_grid))
    lam_bar = np.zeros((n_t, n_grid))
    pcm_arr = np.zeros((n_t, n_grid))
    valid = np.zeros((n_t, n_grid), dtype=np.int8)
    resid = np.zeros((n_t, n_grid))
    info_t: list[dict[str, Any]] = []

    # ---- H-1 -------------------------------------------------------------------------------
    mu_h = np.linspace(0.0, law5.MU_CUT_PP, opt.n_mu_check)
    mu_full = np.concatenate((-mu_h[::-1], mu_h[1:]))
    sig_h = np.zeros(n_grid)
    row_h = np.zeros((n_grid, n_q + 1))
    min_node = np.full(n_grid, np.inf)
    min_mid = np.full(n_grid, np.inf)
    for k in range(i150 + 1):
        e_ev = float(grid[k]) * 1e6
        d = law5.ni_density_ltp1(sec, e_ev, mu_h, awi, 1, 1)
        min_node[k] = float(d.min())
        sig_h[k] = 2.0 * math.pi * law5.simpson(d, mu_h)
        dfull = np.concatenate((d[::-1], d[1:]))
        row_h[k] = symmetric_quantiles(mu_full, np.maximum(dfull, 0.0), n_q)
        if k < i150:
            dm = law5.ni_density_ltp1(sec, 0.5 * (grid[k] + grid[k + 1]) * 1e6, mu_h, awi, 1, 1)
            min_mid[k] = float(dm.min())
    k_pp, p6_points = pp_scan(grid, min_node, min_mid, i150)  # pass 1; P6 failure diagnostic
    p6_found = bool(p6_points)  # recorded: the frozen P6 item (4) fails below E_min,pp
    mu_fine = np.linspace(0.0, law5.MU_CUT_PP, 2 * opt.n_mu_check - 1)
    pp_check = verify_pp_domain(  # pass 2 (P6-D (1)): independent finer-mu re-evaluation
        grid, k_pp, i150,
        lambda e: float(law5.ni_density_ltp1(sec, e * 1e6, mu_fine, awi, 1, 1).min()),
    )  # fmt: skip
    in_domain_ok = True  # verify_pp_domain raised otherwise
    s_ref = bgg.pp_bgg_mb(E_ANCHOR_MEV)
    for k in range(k_pp, n_grid):
        if k <= i150:
            sigma[0, k] = sig_h[k]
            edges[0, k] = row_h[k]
        else:
            sigma[0, k] = sig_h[i150] * bgg.pp_bgg_mb(float(grid[k])) / s_ref
            edges[0, k] = row_h[i150]
    edges[0, :k_pp] = row_h[k_pp] if k_pp <= i150 else row_h[i150]
    valid[0, k_pp:] = 1
    for k in range(n_grid):
        pcm_arr[0, k] = kinematics_cm(float(grid[k]), MP)[0]
    h_mid_err = 0.0
    for k in range(k_pp, i150):
        e_mid = 0.5 * (grid[k] + grid[k + 1])
        direct = (
            2.0
            * math.pi
            * law5.simpson(law5.ni_density_ltp1(sec, e_mid * 1e6, mu_h, awi, 1, 1), mu_h)
        )
        h_mid_err = max(h_mid_err, abs(0.5 * (sigma[0, k] + sigma[0, k + 1]) - direct) / direct)
    info_t.append(
        {"name": PP_TARGET, "z": 1, "a": 1, "kind": "law5-ltp1", "e_min_pp_mev": float(grid[k_pp]),
         "e_min_pp_index": k_pp, "negative_density_below_e_min_pp": p6_found,
         "p6_negative_density_found": p6_found, "first_negative_nodes": p6_points,
         "e_min_pp": float(grid[k_pp]), "no_negative_density_in_domain": in_domain_ok,
         "p6d_pass2": pp_check,
         "sigma_150_mev_barn": float(sigma[0, i150]), "s_150_mb": s_ref,
         "max_midpoint_interpolation_error": h_mid_err,
         "half_sphere": "sigma counts each event once (0 <= mu <= 0.96, 2 pi); edges cover the "
         "symmetric range [-0.96, 0.96] of the same density"}
    )  # fmt: skip
    log(f"H-1: E_min_pp = {grid[k_pp]:.4g} MeV, sigma(150) = {sigma[0, i150] * 1e3:.3f} mb")
    pp_rec: dict[str, Any] | None = None
    hale_curve = [  # diagnostic (report only): the table's own Hale-NI partial sigma_pp(E)
        {"e_mev": float(grid[k]), "sigma_pp_partial_mb": 1e3 * float(sigma[0, k])}
        for k in range(k_pp, n_grid)
    ]  # fmt: skip
    if m2 is not None and m2.pp_level != "none" and m2model is not None:
        pp_rec = _pp_reanchor(m2, m2model, grid, sigma, sig_h, i150, k_pp, s_ref, awi)

    # ---- p + A -----------------------------------------------------------------------------
    sigma_10: dict[str, float] = {}
    max_dev = 0.0
    for it, spec in enumerate(TARGETS, start=1):
        tt, model = tts[spec.name], models[spec.name]
        s_nonel = tt.sigma_barn(grid)
        sig_nonel[it] = s_nonel
        for k in range(n_grid):
            r_fm[it, k], lam_bar[it, k] = disk_radius_fm(
                float(s_nonel[k]), float(grid[k]), model.m_t_mev
            )
            pcm_arr[it, k] = kinematics_cm(float(grid[k]), model.m_t_mev)[0]
        good = (s_nonel > 0.0) & (r_fm[it] > 0.0)
        bad_idx = np.nonzero(~good)[0]
        k_min = int(bad_idx.max()) + 1 if bad_idx.size else 0
        if k_min >= n_grid:
            raise BuildError(f"{spec.name}: no node with a positive disk radius")
        check_domain_limit(spec.name, float(grid[k_min]))
        sig_mb = np.array([bgg.bgg_elastic_mb(float(e), spec.z, model.m_t_mev) for e in grid])
        sigma[it, k_min:] = sig_mb[k_min:] * 1e-3
        valid[it, k_min:] = 1
        for k in range(k_min, n_grid):
            edges[it, k], _ = disk_quantiles(
                float(r_fm[it, k]), float(pcm_arr[it, k]), opt.n_shape, n_q
            )
            resid[it, k] = abs(
                math.pi * (r_fm[it, k] + lam_bar[it, k]) ** 2 / (100.0 * s_nonel[k]) - 1.0
            )
            if k % 11 == 0:
                max_dev = max(
                    max_dev,
                    _shape_norm_deviation(float(r_fm[it, k]), float(pcm_arr[it, k]), opt.n_shape),
                )
        if m2model is not None:
            for k in range(k_min, n_grid):
                if grid[k] < m2model.nn.t_min_mev:  # below the SAID range: black disk kept
                    m2_below[spec.name] += 1
                    continue
                if m2 is not None and m2.pauli and grid[k] < (m2.pauli_emin_mev or 0.0):
                    # below pauli_emin: the pinned black-disk node exactly (edges untouched)
                    m2_pinned[spec.name].append({"e_mev": float(grid[k]), "pauli_applied": False})
                    continue
                edges[it, k], rec = m2model.node(
                    spec.name, float(grid[k]), float(pcm_arr[it, k]), n_q
                )
                rec["sigma_el_bgg_mb"] = 1e3 * float(sigma[it, k])
                rec["m2_over_bgg"] = rec["sigma_el_m2_mb"] / rec["sigma_el_bgg_mb"]
                m2_nodes[spec.name].append(rec)
                m2_norm = max(m2_norm, rec["norm_deviation"])
                if (k - k_min) % 25 == 0 or k in cand_k:  # + the plan's candidate energies
                    rec["convergence"] = m2model.convergence(
                        spec.name, float(grid[k]), float(pcm_arr[it, k])
                    )
                    rec["convergence"]["candidate_energy"] = k in cand_k
                    m2_conv = max(m2_conv, rec["convergence"]["max"])
        if m3model is not None and m3 is not None and m3model.has(spec.name):
            emax = m3model.emax_mev(spec.name)
            pcm_top = kinematics_cm(emax, model.m_t_mev)[0]
            for k in range(k_min, n_grid):
                e = float(grid[k])
                if e > emax * (1.0 + 1e-12):
                    m3_above[spec.name] += 1
                    if m3.above_emax == "error":
                        raise BuildError(
                            f"{spec.name}: node {e:.4g} MeV is above the ENDF evaluation range "
                            f"({emax:g} MeV); choose m3.above_emax = hold-last | black-disk"
                        )
                    if m3.above_emax == "black-disk":
                        continue
                edges[it, k], rec = m3model.node(spec.name, e, n_q, float(pcm_arr[it, k]), pcm_top)
                m3_nodes[spec.name].append(rec)
                m3_norm = max(m3_norm, rec["norm_deviation"])
            # positivity at the energy midpoints of the domain as well (P6-D style, no clipping)
            for k in range(k_min, n_grid - 1):
                if m3.above_emax == "black-disk" and grid[k] > emax * (1.0 + 1e-12):
                    continue
                m3model.check_midpoint(spec.name, 0.5 * float(grid[k] + grid[k + 1]))
                m3_mid[spec.name] += 1
        edges[it, :k_min] = edges[it, k_min]
        sigma_10[spec.name] = bgg.bgg_elastic_mb(10.0, spec.z, model.m_t_mev)
        info_t.append(
            {"name": spec.name, "z": spec.z, "a": spec.a, "kind": "bgg-black-disk",
             "e_min_shape_mev": float(grid[k_min]), "e_min_shape_index": k_min,
             "sigma_el_10_mev_mb": sigma_10[spec.name],
             "sigma_el_14_mev_mb": bgg.bgg_elastic_mb(14.0, spec.z, model.m_t_mev),
             "sigma_el_100_mev_mb": bgg.bgg_elastic_mb(100.0, spec.z, model.m_t_mev),
             "sigma_el_150_mev_mb": bgg.bgg_elastic_mb(150.0, spec.z, model.m_t_mev),
             "radius_fm_100_150_200_mev": [
                 float(np.interp(e, grid, r_fm[it])) for e in (100.0, 150.0, 200.0)
             ],
             "max_inversion_residual": float(resid[it].max())}
        )  # fmt: skip
        log(
            f"{spec.name}: e_min_shape {grid[k_min]:.4g} MeV, "
            f"R(150) = {np.interp(150.0, grid, r_fm[it]):.3f} fm"
        )

    # ---- table vs direct sigma_el at every node midpoint ------------------------------------
    max_mid = 0.0
    for it, spec in enumerate(TARGETS, start=1):
        k_min = info_t[it]["e_min_shape_index"]
        for k in range(max(k_min, 0), n_grid - 1):
            m = 0.5 * (grid[k] + grid[k + 1])
            d_b = bgg.bgg_elastic_mb(float(m), spec.z, models[spec.name].m_t_mev) * 1e-3
            if d_b > 1e-3 * sigma[it, -1]:
                max_mid = max(max_mid, abs(0.5 * (sigma[it, k] + sigma[it, k + 1]) - d_b) / d_b)

    # ---- X-ENDF cross-check tables (report-only) ---------------------------------------------
    th = np.asarray(XENDF_THETA_DEG)
    mu_x = np.cos(np.radians(th))
    xe_la150 = np.zeros((len(XENDF_TARGETS), len(XENDF_ENERGIES_MEV), th.size))
    xe_model = np.zeros_like(xe_la150)
    for ix, name in enumerate(XENDF_TARGETS):
        spec = next(t for t in TARGETS if t.name == name)
        s5, tab = law5.parse_law5(mats[name].sections[(6, 2)]), mats[name].cross_section(2)
        for ie, e in enumerate(XENDF_ENERGIES_MEV):
            idx = [r.energy_ev for r in s5.records].index(e * 1e6)
            sg = float(tab.interpolate(e * 1e6)[0])
            xe_la150[ix, ie] = 1e3 * ltp12_ni_density_ratio(
                s5, idx, sg, mu_x, law5.projectile_awi(mats[name]), spec.z
            )
            it = 1 + [t.name for t in TARGETS].index(name)
            xe_model[ix, ie] = disk_density_mb_sr(
                float(np.interp(e, grid, sigma[it])), float(np.interp(e, grid, r_fm[it])),
                float(np.interp(e, grid, pcm_arr[it])), th,
            )  # fmt: skip

    if m3model is not None:
        max_dev = m3_norm
    if m2model is not None:
        max_dev = max(m2_norm, 0.0) if m2_conv <= glauber.CONVERGENCE_TOL else float("inf")
    alt_kind: str | None = None
    if m2model is not None and m2 is not None:
        alt_kind = f"m2prime-{m2.pauli_variant}" if m2.pauli else "m2-glauber"
    for it, spec in enumerate(TARGETS, start=1):
        if alt_kind is not None:
            info_t[it]["kind"] = alt_kind
        elif m3model is not None and m3model.has(spec.name):
            info_t[it]["kind"] = "m3-endf-mirror"
    alt_variant = m2model is not None or m3model is not None
    if alt_variant:  # the cross-check array is the black-disk formula: omitted, with a reason
        xe_model[:] = np.nan
    qual = np.array(
        [
            bool(k_pp <= i150 and in_domain_ok),
            bool(o16_mt5["o16_mt5_sigma_not_c12_copy"]),
            bool(o16_spec["o16_mt5_spectra_independent"]),
            bool(max_dev <= (glauber.NORM_TOL if (m2model or m3model) else 1e-6)),
        ],
        dtype=np.int8,
    )
    e_min_t = np.array([float(grid[k_pp])] + [info_t[i]["e_min_shape_mev"] for i in range(1, n_t)])
    # NI cross section above the cut and the expected omitted events below the domain
    sigma_ni_mb = {
        f"{e:g}_mev": 1e3 * 2.0 * math.pi * law5.simpson(
            law5.ni_density_ltp1(sec, e * 1e6, mu_h, awi, 1, 1), mu_h
        )
        for e in (15.0, 20.0)
    }  # fmt: skip
    i_o = [t.name for t in TARGETS].index("O-16") + 1
    spec_o = next(t for t in TARGETS if t.name == "O-16")
    e_o = float(e_min_t[i_o])
    basis = (
        "per history of a proton slowing down from the cut energy to 1 MeV in liquid water "
        "(Bethe stopping table, 1 g/cm3), nuclear attrition ignored"
    )
    pp_omit = pp_omitted_ni_correction(sec, awi, grid, k_pp, mu_h)
    pp_omit["residual_range_at_cut_g_cm2"] = residual_range_g_cm2(float(grid[k_pp]))
    pp_omit["basis"] = basis
    pa_omitted = {
        "basis": basis, "target": "O-16", "e_cut_mev": e_o, "n_o_cm3": N_O_WATER_CM3,
        "sigma": "BGG sigma_el over [1 MeV, e_min_shape]; not evaluated below 1 MeV",
        "events_per_history": omitted_events_per_history(
            e_o, N_O_WATER_CM3, lambda e: bgg.bgg_elastic_mb(e, spec_o.z, models["O-16"].m_t_mev)
        ),
        "residual_range_g_cm2": residual_range_g_cm2(e_o),
        "residual_range_below_1mev_g_cm2": pp_omit["residual_range_below_1mev_g_cm2"],
        "note": "BGG sigma_el is defined down to 1 MeV; the events of the final residual range "
        "below 1 MeV are not estimated (not extrapolated)",
    }  # fmt: skip
    elastic_domain = {PP_TARGET: [float(grid[k_pp]), E_MAX_MEV]}
    elastic_domain.update(
        {t.name: [float(e_min_t[i + 1]), E_MAX_MEV] for i, t in enumerate(TARGETS)}
    )
    for i, t in enumerate(TARGETS, start=1):
        info_t[i]["sigma_bgg_below_domain_max_mb"] = max(
            bgg.bgg_elastic_mb(float(e), t.z, models[t.name].m_t_mev)
            for e in np.linspace(1.0, float(e_min_t[i]), 200)
        )
    info_t[0]["sigma_bgg_below_domain_max_mb"] = None  # H-1: not BGG
    if exploratory_forced_isotropic:
        iso = np.linspace(-1.0, 1.0, n_q + 1)
        for it, spec in enumerate(TARGETS, start=1):
            k_min = int(info_t[it]["e_min_shape_index"])
            for k in range(k_min):
                sigma[it, k] = (
                    bgg.bgg_elastic_mb(float(grid[k]), spec.z, models[spec.name].m_t_mev) * 1e-3
                )
                valid[it, k] = 1
                edges[it, k] = iso
    arr: dict[str, NDArray[Any]] = {
        "grid_e_mev": grid, "sigma_barn": sigma, "edges_mu": edges, "radius_fm": r_fm,
        "sigma_nonel_barn": sig_nonel, "lambda_bar_fm": lam_bar, "p_cm_mev": pcm_arr,
        "valid": valid, "inversion_residual": resid, "qualification": qual,
        "target_e_min_mev": e_min_t,
        "xendf_theta_deg": th, "xendf_energies_mev": np.asarray(XENDF_ENERGIES_MEV),
        "xendf_la150_ni_mb_sr": xe_la150, "xendf_model_mb_sr": xe_model,
        "target_z": np.array([1] + [t.z for t in TARGETS], dtype=np.int64),
        "target_a": np.array([1] + [t.a for t in TARGETS], dtype=np.int64),
        "target_mass_mev": np.array([MP] + [models[t.name].m_t_mev for t in TARGETS]),
    }  # fmt: skip
    out_dir = Path(out_dir) if out_dir is not None else cdir / "derived"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = out_dir / "elastic-proton-building.npz"
    npz_sha = write_npz_deterministic(tmp, arr)
    sig_sha = {
        t.name: hashlib.sha256(np.ascontiguousarray(sig_nonel[i + 1]).tobytes()).hexdigest()
        for i, t in enumerate(TARGETS)
    }
    info: dict[str, Any] = {
        "schema": SCHEMA, "builder_version": BUILDER_VERSION, "options": asdict(opt),
        "sources": {
            s: {"sha256": DATASETS[s].sha256, "version": DATASETS[s].version} for s in SOURCE_IDS
        },
        "npz_sha256": npz_sha, "npz_arrays": {k: list(v.shape) for k, v in sorted(arr.items())},
        "targets": info_t, "target_names": list(TARGET_NAMES), "elements": element_rows_elastic(),
        "grid": {"e_min_mev": float(grid[0]), "e_max_mev": float(grid[-1]), "n_points": int(n_grid),
                 "anchor_index": i150, "nodes_added_by_refinement_below_14_mev": n_refined,
            "nodes_added_by_h1_refinement": n_refined_pp,
                 "refine_tolerance": REFINE_TOL, "interpolation": "lin-lin in E between nodes",
                 "max_midpoint_sigma_el_interpolation_error": max_mid},
        "edges": {"n_quantiles": n_q, "mu_range": "[-1, 1] for p+A; [-0.96, 0.96] for H-1",
                  "meaning": "n_q + 1 ascending mu_CM edges of n_q equiprobable bins"},
        "units": {"sigma_barn": "barn per nucleus [target, node]; H-1 half sphere",
                  "radius_fm": "black-disk R [fm]: pi (R + lambdabar)^2 = sigma_nonel",
                  "sigma_nonel_barn": "the transport's own LA150 MT5 (+ Tripathi shape > 150 MeV)"},
        "sigma_nonel_sha256": sig_sha,
        "convention": {"identical_particles": "half sphere, 2 pi, once per event, ratio interp.",
                       "pp_cut_mu": law5.MU_CUT_PP,
                       "pp_above_150_mev": "sigma_NI(150) S(E)/S(150), shape fixed in mu_CM"},
        "sigma_el_10_mev_mb": sigma_10,
        "o16_finding": o16, "o16_mt2_copy_recorded": o16["o16_mt2_is_c12_copy"],
        "o16_mt5_finding": o16_mt5, "o16_mt5_spectra_finding": o16_spec,
        "transcription_verified": transcription,
        "elastic_domain": elastic_domain, "model_revisions": MODEL_REVISIONS,
        "negative_density_below_e_min_pp": p6_found,
        "p6": {"negative_density_below_e_min_pp": p6_found, "first_negative_nodes": p6_points,
               "e_min_pp": float(grid[k_pp]), "no_negative_density_in_domain": in_domain_ok,
               "note": "the frozen row P6 item (4) fails for H-1 below e_min_pp; recorded as a "
               "failure, not converted into a pass by the revised domain (row P6-D)"},
        "sigma_ni_above_cut_mb": sigma_ni_mb,
        "pp_omitted_ni_correction_total_variation_per_history_150mev": pp_omit[
            "total_variation_per_history"
        ],
        "pp_omitted_ni_correction_weighted_total_variation_mev": pp_omit[
            "weighted_total_variation_mev"
        ],
        "pp_omitted_ni_correction": pp_omit,
        "pa_omitted_events_per_history_150mev": pa_omitted,
        "xendf": {"targets": list(XENDF_TARGETS), "energies_mev": list(XENDF_ENERGIES_MEV),
                  "report_only": True,
                  "note": "LA150 NI contains Coulomb-nuclear interference, the model does not"},
        "shape_normalisation_max_deviation": max_dev,
        "qualification_fields": list(QUALIFICATION_FIELDS),
    }  # fmt: skip
    if m2model is not None and m2 is not None:
        info["exploratory"] = {
            "model": "m2prime-pauli-blocked" if m2.pauli else "m2-optical-glauber",
            "variant": m2.variant(),
            "pauli_variant": m2.pauli_variant if m2.pauli else None,
            "pauli_emin_mev": m2.pauli_emin_mev if m2.pauli else None,
            "beta": m2.beta,
            "pp_level": m2.pp_level, "unfold": m2.unfold,
            "inputs_unpinned": bool(m2model.unpinned),
            "purpose": "V3-005D diagnostic (plan 4.1 / Amendment 20), report-only, never evidence",
            "change": "p + A mu_CM shapes = optical-limit Glauber (sigma_el, domains, grid, H-1 "
            "shape as the model table)" + ("; H-1 level re-anchored >= 125 MeV"
                                           if m2.pp_level != "none" else ""),
        }  # fmt: skip
        rec_m2 = m2_record(m2, m2model, m2_nodes, m2_conv, m2_norm, pp_rec, m2_pinned)
        info["m2"] = _nan_to_none(rec_m2)
        info["m2"]["nodes_below_said_range_black_disk_kept"] = {
            "t_min_mev": m2model.nn.t_min_mev,
            "count_per_target": m2_below,
            "note": "SAID tables start at 5 MeV and are not extrapolated; below it the pinned "
            "black-disk shape is kept (e_min_shape of N-14, Al-27, Si-28, P-31 is < 5 MeV)",
        }
        info["m2"]["hale_ni_partial_pp_curve_before_reanchoring"] = hale_curve
    if alt_variant:
        info["xendf"]["model_omitted_reason"] = (
            "xendf_model_mb_sr is the black-disk formula and does not describe the shape of this "
            "alternative-model table (kind per target); the array is NaN"
        )
    if m3model is not None and m3 is not None:
        info["exploratory"] = {
            "model": "m3-endf-mirror", "above_emax": m3.above_emax,
            "inputs_unpinned": bool(m3model.unpinned),
            "endf": {n: {"sha256": f["sha256"], "mat": f["mat"]} for n, f in m3model.files.items()},
            "purpose": "V3-005D fallback shape (plan row (ii'), isospin mirror), never evidence",
            "change": "C-12 and O-16 mu_CM shapes = ENDF/B-VIII.0 neutron MF4/MT2 at the same lab "
            "energy; sigma_el, domains, grid, H-1 and all other targets as the model table",
        }  # fmt: skip
        info["m3"] = {
            "provenance": m3model.provenance(),
            "targets": {
                t.name: (
                    {"shape": "endf-mf4-mt2", "n_nodes_endf": len(m3_nodes[t.name]),
                     "n_nodes_above_emax": m3_above[t.name],
                     "above_emax_policy": m3.above_emax,
                     "n_midpoints_checked_positive": m3_mid[t.name],
                     "n_nodes_q_scaled": sum(bool(r["q_scaled"]) for r in m3_nodes[t.name]),
                     "e_top_mev": m3model.emax_mev(t.name),
                     "coverage_limitation": m3model.emax_mev(t.name) < 150.0,
                     "min_density": min((r["min_f"] for r in m3_nodes[t.name]), default=None),
                     "max_norm_deviation_of_endf_density": max(
                         (abs(r["norm_int_f_dmu"] - 1.0) for r in m3_nodes[t.name]), default=None
                     ),
                     "max_first_moment_1_minus_mu": max(
                         (r["first_moment_1_minus_mu"] for r in m3_nodes[t.name]), default=None
                     ),
                     "gates": "MF1 identity, copy check (181 pts), normalisation <= 1e-6, "
                     "positivity at nodes and midpoints (no clipping), monotone CDF, finite first "
                     "moment, E_top recorded; any failure raises (no table)"}
                    if m3model.has(t.name) else
                    {"shape": "black-disk (pinned; not named in the plan row (ii'))"}
                )
                for t in TARGETS
            },
            "note": "ENDF/B-VIII.0 neutron MF4/MT2 evaluations verified to extend to 150 MeV "
            "(C-12 tabulated 20-150 MeV, O-16 30-150 MeV; record validation/results/"
            "reference_data/endf-b-viii0-neutrons-c12-o16, 2026-10-23); no data above "
            "150 MeV; sigma_el stays BGG",
        }  # fmt: skip
    if exploratory_forced_isotropic:
        info["exploratory"] = {
            "exploratory_forced_isotropic": True,
            "purpose": "E-shape (Amendment 15 (b)6), report-only, never evidence",
            "change": "p + A nodes below e_min_shape: BGG sigma_el (valid) with mu_CM isotropic on "
            "[-1, 1]; everything inside the declared domain, the grid and H-1 as the model table",
        }
    tid = table_id(info)
    info["table_id"] = tid
    npz_path, json_path = (
        out_dir / f"elastic-proton-{tid}.npz",
        out_dir / f"elastic-proton-{tid}.json",
    )
    tmp.replace(npz_path)
    json_path.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    out = dict(info)
    out["timing_s"] = time.perf_counter() - t_start
    return ElasticBuildResult(npz_path, json_path, tid, out)
