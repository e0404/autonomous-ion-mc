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
from ionmc.nuclear import events as ev
from ionmc.nuclear.build import (
    TARGETS,
    BuildError,
    TargetTables,
    element_rows,
    table_id,
    union_grid,
    write_npz_deterministic,
)

BUILDER_VERSION = "ionmc-elastic-proton-builder-4"
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
QUALIFICATION_FIELDS = ("no_negative_density_in_domain", "o16_mt5_not_c12_copy", "shape_normalised")
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
    "are report-only); the o16_mt5_copy_check is implemented."
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
    shared) the union grid with linear interpolation. ``o16_mt5_not_c12_copy`` is true only if
    (a) sigma: ``|s_C - s_O| / max`` > 1e-3 at > 50 % of the energies AND the coefficient of
    variation of ``s_O / s_C`` > 1e-3 (a scaled copy has a constant ratio), and (b) the revised MF6
    rule (specification revision before any transport result): per product of MF6/MT5 the
    denominator is the number of INFORMATIVE compared energies (at least one material nonzero; the
    relative difference is 1.0 where exactly one yield is zero); a product with at least
    ``MF6_MIN_INFORMATIVE`` informative energies gates and must differ by > 1e-3 at > 50 % of them,
    the others are report-only. The rule as frozen (all compared energies as denominator) is
    recorded as ``mf6_rule_as_frozen`` with its pass flag. A missing product counts as a zero
    yield. Parser or data errors propagate as :class:`BuildError` (never swallowed), and so does a
    failed check (fail closed)."""
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
            "role": "gating" if gating else "report-only",
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
        "rule": "denominator = informative compared energies (at least one material nonzero); a "
        f"product gates only with >= {MF6_MIN_INFORMATIVE} informative energies, else report-only; "
        "every gating product must exceed 50 %",
        "min_informative": MF6_MIN_INFORMATIVE, "products": revised_prod,
        "passes": all(v["differs_from_c12"] for v in revised_prod.values()
                      if v is not None and v["role"] == "gating"),
    }  # fmt: skip
    yields_ok = bool(mf6_revised["passes"])
    out: dict[str, Any] = {
        "compared_grid": grid_kind, "n_compared": int(e_ev.size),
        "compared_from_mev": float(e_ev[0] * 1e-6), "compared_to_mev": float(e_ev[-1] * 1e-6),
        "nodes_compared": int(rel.size), "fraction_rel_diff_gt_1e-3": frac,
        "max_relative_difference": float(rel.max()), "min_relative_difference": float(rel.min()),
        "median_relative_difference": float(np.median(rel)),
        "threshold": REL_DIFF_THRESHOLD, "required_fraction_gt": 0.5,
        "ratio_o16_over_c12_cv": cv, "required_cv_gt": REL_DIFF_THRESHOLD,
        "sigma_not_copy": sigma_ok, "mf6_yields_not_copy": yields_ok,
        "mf6_rule_as_frozen": mf6_frozen, "mf6_rule_revised": mf6_revised,
        "o16_mt5_not_c12_copy": bool(sigma_ok and yields_ok),
    }  # fmt: skip
    if not out["o16_mt5_not_c12_copy"]:
        raise BuildError(
            f"LA150 O-16 MF3/MT5 or MF6/MT5 looks like a copy of C-12 (fail closed): {out}"
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


def build_elastic_proton(
    cache_dir: str | Path | None = None,
    options: ElasticBuildOptions | None = None,
    log: Callable[[str], None] = _noop,
) -> ElasticBuildResult:
    """Build the table (module docstring). ``BuildError`` on any fail-closed rule."""
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
    o16 = o16_copy_finding(mats["C-12"], mats["O-16"])
    log(f"O-16 MT2 is a C-12 copy from {o16['from_mev']} MeV: recorded")
    o16_mt5 = o16_mt5_copy_check(mats["C-12"], mats["O-16"])
    log(f"O-16 MT5 differs from C-12 at {o16_mt5['fraction_rel_diff_gt_1e-3']:.0%} of the nodes")

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

    qual = np.array(
        [
            bool(k_pp <= i150 and in_domain_ok),
            bool(o16_mt5["o16_mt5_not_c12_copy"]),
            bool(max_dev <= 1e-6),
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
    out_dir = cdir / "derived"
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
        "o16_mt5_finding": o16_mt5, "transcription_verified": transcription,
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
