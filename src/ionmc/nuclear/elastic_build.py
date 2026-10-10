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
    ratified); ``sigma = 0`` below it, the value and the first-negative diagnostics are recorded.
    150-250 MeV: "sigma(E) = sigma_NI(150) S(E)/S(150), with S = 1.0115 (23 + 50 sqrt(ln(0.73/p)^7))
    mb (BGG/PDG systematics, p_lab in GeV/c). The shape is held fixed in mu_CM."

X-ENDF (report-only) and the O-16 finding
    LA150 C-12/N-14/Ca-40 LTP=12 NI densities (ratio interpolation in mu) at 20-40 deg for
    50/100/150
    MeV are stored for the cross-check row. The builder asserts, fail closed, that the O-16 MT2 data
    are a numerical copy of C-12 above 24 MeV (decision 0041 data finding) and never uses LAW=5 data
    of any target except H-1 as a construction input.
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

BUILDER_VERSION = "ionmc-elastic-proton-builder-1"
SCHEMA = "ionmc-elastic-proton-table-1"
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
    "pdg-rpp2022-pp-elastic",
)
E_MAX_MEV = 250.0
E_ANCHOR_MEV = 150.0
N_Q = 256
HBARC = 197.3269804
MP = 938.27208816
PP_TARGET = "H-1"
TARGET_NAMES = (PP_TARGET, *(t.name for t in TARGETS))
QUALIFICATION_FIELDS = ("no_negative_density", "o16_copy_asserted", "shape_normalised")
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
    ok = (min_node >= 0.0) & (min_mid >= 0.0)
    ok[i150 + 1 :] = True
    bad = np.nonzero(~ok)[0]
    k_pp = int(bad.max()) + 1 if bad.size else 0
    if k_pp > i150:
        raise BuildError("no energy up to 150 MeV has a non-negative p-p density")
    first_neg = [
        {"e_mev": float(grid[k]), "min_node_mb_sr": 1e3 * float(min_node[k]),
         "min_midpoint_mb_sr": 1e3 * float(min_mid[k])}
        for k in bad[-12:]
    ]  # fmt: skip
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
         "e_min_pp_index": k_pp, "first_negative_nodes": first_neg,
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
        [bool(k_pp <= i150), bool(o16["o16_mt2_is_c12_copy"]), bool(max_dev <= 1e-6)], dtype=np.int8
    )
    arr: dict[str, NDArray[Any]] = {
        "grid_e_mev": grid, "sigma_barn": sigma, "edges_mu": edges, "radius_fm": r_fm,
        "sigma_nonel_barn": sig_nonel, "lambda_bar_fm": lam_bar, "p_cm_mev": pcm_arr,
        "valid": valid, "inversion_residual": resid, "qualification": qual,
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
        "o16_finding": o16, "transcription_verified": transcription,
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
