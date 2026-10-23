"""D-E9 / B0 (exploratory, report-only, never evidence): an INDEPENDENT reconstruction of the pinned table's p-p level, 100-150 MeV.

The check reads the raw ENDF/B-VIII.0 (LA150) H-1 file of the cache (the zip member ``p-001_H_001.endf`` of the dataset the pinned build
uses; its sha256 is read from the pinned table's sidecar) and implements the LAW=5, LTP=1, LIDP=1 reconstruction anew. It does NOT import
``elastic_build``, ``ionmc.data.endf6`` or any LAW=5 code of the package (independence rule of the D-E9 plan); it reads only the pinned npz /
json (node energies and the table's own node values) and the cache objects. The ENDF-6 formats manual (ENDF-102) definitions used:

* MF6 / MT2 record layout, Sec. 6.2.1 (CONT, TAB1, TAB2, LIST) and LAW=5 (Sec. 6.2.7): ``[MAT,6,2/ SPI,0,LIDP,0,NR,NE/ E_int]TAB2`` then per
  energy ``[MAT,6,2/ 0,E,LTP,0,NW,NL/ A]LIST``; LTP=1, LIDP=1: ``NW = 3 NL + 3`` (``b_0..b_NL``, then ``Re a_l, Im a_l``, ``l = 0..NL``).
* eq. 6.10 (identical particles, spin s = 1/2): ``sigma_ci = 2 eta^2/(k^2 (1-mu^2)) [ (1+mu^2)/(1-mu^2) - cos(eta ln((1+mu)/(1-mu)))/2 ]``.
* eq. 6.14 (LTP=1, identical particles): ``sigma_ei = sigma_ci - 2 eta/(1-mu^2) Re{ sum_l [ (1+mu) exp(i eta ln((1-mu)/2)) +
  (-1)^l (1-mu) exp(i eta ln((1+mu)/2)) ] (2l+1)/2 a_l P_l(mu) } + sum_l (4l+1)/2 b_l P_2l(mu)``.
* eqs. 6.11 / 6.12 (non-relativistic): ``k = A/(1+A) sqrt(2 m_1 E)/hbar`` (b^-1/2), ``eta = z_1 z_2 alpha sqrt(m_1 c^2/(2E))``.
* The transported density is ``NI = sigma_e - sigma_c``; the table interpolates the RATIO ``sigma_e/sigma_c`` of the two bracketing energy
  nodes linearly in E (TAB2 INT = 2) and multiplies with the Coulomb density at E. The cross section is ``2 pi int_0^0.96 NI dmu`` (barn;
  half sphere, each identical-particle event counted once; ``|mu_CM| <= 0.96`` = theta_CM >= 16.26 deg).

Numerics differ from the builder on purpose: Gauss-Legendre quadrature with 1024 nodes on [0, 0.96] (the builder: a 4001-point uniform Simpson grid, ``n_mu_check``; its P6-D re-evaluation uses 8001; a planner note fixes the plan text, which spoke of a "four times finer" grid). The quadrature error is not assumed: ``quadrature_convergence`` records GL 1024 against GL 2048, Simpson 4001 and Simpson 20001 at the 100 / 125 / 150 MeV nodes (acceptance, pre-declared: relative difference < 1e-6, negligible against the 1e-3 match tolerance),
Legendre sums by ``numpy.polynomial.legendre.legval``, an own fixed-field ENDF reader.

Planning note (plan, B0): the plan names a 1 MeV grid over 100-200 MeV, but the raw LA150 H-1 evaluation ends at 150 MeV, so the match and the dip test stay on
the raw-data domain 100-150 MeV. Above 150 MeV the builder carries the 150 MeV level by the table's own S(E) ratio; ``carry_up_check`` reconstructs the table level
above 150 MeV independently from the reconstructed 150 MeV level x ``bgg.pp_bgg_mb(E) / bgg.pp_bgg_mb(150)`` (``ionmc.nuclear.bgg`` only: the Geant4 p-p formula,
not LAW=5 code) and compares with the pinned nodes (ratios; match within 1e-3), with the 1 MeV grid 150-200 MeV recorded as ratios to the 150 MeV level.

Definitions fixed before computing (plan section 2, B0):
* "match": relative difference reconstruction / table - 1 within 1e-3 at every table node and every node midpoint (linear interpolation of the
  table node values, the builder's own midpoint convention) in [100, 150] MeV;
* "dip persists": an interior minimum of the reconstructed level on the 1 MeV grid 100-150 MeV that lies more than 1 % below BOTH endpoint values;
* reading: match AND dip -> "evaluation feature at the LA150 top nodes"; otherwise -> "builder defect route". No other reading is offered.

Only ratios, booleans and node energies are recorded: no absolute cross section. ``leak_violations`` checks the projection against every
absolute value computed here.

Run: ``python -m validation.scripts.elastic_attribution.de9_pp_level``.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import zipfile
from pathlib import Path
from typing import Any

import numpy as np

from . import common as c

MEMBER = "ENDF-B-VIII.0_protons/p-001_H_001.endf"
SOURCE_KEY = "endf-b8.0-protons"
MU_CUT = 0.96
N_GAUSS = 1024
TOL_MATCH = 1.0e-3
DIP_MARGIN = 0.01
E_LO, E_HI = 100.0, 150.0
NAMED_NODES = (100.0, 125.0, 130.0, 140.0, 150.0)
AMU_EV, MN_AMU, HBARC_EV_FM, ALPHA = 931.49410242e6, 1.00866491595, 197.3269804e6, 1.0 / 137.035999084
STEP_DIR = "de9"

GIT_ALLOW = (
    "banner", "disclosure", "named_nodes/*/*", "nodes/*/*", "midpoints/*/*", "match/*", "dip/*", "levels_over_125mev/*/*", "reading",
    "carry_up_check/**", "quadrature_convergence/**", "provenance/*",
)  # fmt: skip
GIT_FILES = {"de9-b0.json": (GIT_ALLOW, ())}
PRE_RUN_LOOKS = (
    "Before the frozen run the independent parser was validated on the real LA150 H-1 file at the 100 MeV node only (reconstruction / table - 1 = 2e-12 relative); "
    "no other node, no midpoint, no dip test and no verdict was computed before it. The reading rules were declared in the plan before that look. The quadrature convergence was also looked at that one node (not a level comparison)."
)
DISCLOSURE = "D-E9 B0: exploratory; leads identified after the D-E8 STOP; not a re-evaluation. Independent reconstruction of the pinned p-p level; ratios only."


# ---------------------------------------------------------------------------------------------------------------------------------
# own ENDF-6 reader (fixed 11-character fields)
# ---------------------------------------------------------------------------------------------------------------------------------
def endf_float(field: str) -> float:
    """``1.001000+3`` -> 1001.0 (ENDF omits the exponent letter); blank -> 0."""
    s = field.strip()
    if not s:
        return 0.0
    for i in range(len(s) - 1, 0, -1):
        if s[i] in "+-" and s[i - 1] not in "eEdD+-":
            return float(s[:i] + "e" + s[i:])
    return float(s.replace("D", "e"))


def section_lines(text: str, mf: int, mt: int) -> list[str]:
    """The data lines (columns 1-66) of section (MF, MT) in file order."""
    out = []
    for ln in text.splitlines():
        if len(ln) >= 75 and ln[70:72].strip() == str(mf) and ln[72:75].strip() == str(mt):
            out.append(ln[:66])
    return out


class _Cursor:
    def __init__(self, lines: list[str]) -> None:
        self.lines, self.i = lines, 0

    def fields(self) -> list[str]:
        ln = self.lines[self.i]
        self.i += 1
        return [ln[j * 11 : (j + 1) * 11] for j in range(6)]

    def head(self) -> tuple[float, float, int, int, int, int]:
        f = self.fields()
        return endf_float(f[0]), endf_float(f[1]), int(endf_float(f[2])), int(endf_float(f[3])), int(endf_float(f[4])), int(endf_float(f[5]))

    def values(self, n: int) -> list[float]:
        out: list[float] = []
        while len(out) < n:
            out += [endf_float(x) for x in self.fields()]
        return out[:n]


def parse_h1_mt2(text: str) -> dict[str, Any]:
    """MF6 / MT2 of H-1 (one product, LAW=5): ``{awr, spi, lidp, records: [{e_ev, ltp, nl, data}]}``; AWI from MF1 / MT451. Fail closed."""
    cur = _Cursor(section_lines(text, 6, 2))
    _za, awr, _l1, _lct, nk, _ = cur.head()
    if nk != 1:
        raise ValueError(f"MF6 section with NK={nk} products")
    _zap, _awp, _lip, law, nr, np_ = cur.head()
    if law != 5:
        raise ValueError(f"LAW={law}, expected 5")
    cur.values(2 * nr)  # yield TAB1: interpolation pairs, then the points (each block starts on a new line)
    cur.values(2 * np_)
    spi, _c2, lidp, _l2, nr2, ne = cur.head()
    inter = cur.values(2 * nr2)
    if any(int(v) != 2 for v in inter[1::2]):
        raise ValueError("energy interpolation other than lin-lin")
    if lidp != 1:
        raise ValueError(f"LIDP={lidp}: only identical particles (LIDP=1) are implemented")
    recs = []
    for _ in range(ne):
        _c1, e_ev, ltp, _l2b, nw, nl = cur.head()
        if ltp != 1 or nw != 3 * nl + 3:
            raise ValueError(f"LTP={ltp}, NW={nw}, NL={nl}: not an LTP=1 / LIDP=1 record")
        recs.append({"e_ev": e_ev, "ltp": ltp, "nl": nl, "data": np.array(cur.values(nw))})
    if cur.i != len(cur.lines):
        raise ValueError("unexpected records after the LAW=5 distribution")
    en = [r["e_ev"] for r in recs]
    if en != sorted(en) or len(set(en)) != len(en):
        raise ValueError("incident energies must increase strictly")
    m451 = section_lines(text, 1, 451)
    awi = endf_float(m451[2][:11])
    return {"awr": awr, "awi": awi, "spi": spi, "lidp": lidp, "records": recs}


# ---------------------------------------------------------------------------------------------------------------------------------
# eqs 6.10-6.14 and the half-sphere cross section
# ---------------------------------------------------------------------------------------------------------------------------------
def coulomb_k_eta(e_ev: float, awr: float, awi: float) -> tuple[float, float]:
    """eqs 6.11 / 6.12 for z1 = z2 = 1 (proton on proton)."""
    m1 = awi * MN_AMU
    a = awr / awi
    k = a / (1.0 + a) * math.sqrt(2.0 * m1 * AMU_EV * e_ev) / HBARC_EV_FM * 10.0
    eta = math.sqrt(ALPHA**2 * m1 * AMU_EV / (2.0 * e_ev))
    return k, eta


def sigma_coulomb(mu: np.ndarray, k: float, eta: float) -> np.ndarray:
    """eq. 6.10, identical spin-1/2 particles (b/sr)."""
    return (2.0 * eta**2 / (k**2 * (1.0 - mu**2))) * ((1.0 + mu**2) / (1.0 - mu**2) - 0.5 * np.cos(eta * np.log((1.0 + mu) / (1.0 - mu))))


def sigma_e(mu: np.ndarray, rec: dict[str, Any], k: float, eta: float) -> np.ndarray:
    """eq. 6.14 (b/sr): Coulomb + interference + nuclear for one LTP=1, LIDP=1 record."""
    nl, d = rec["nl"], rec["data"]
    b, ca = d[: nl + 1], d[nl + 1 : nl + 1 + 2 * (nl + 1)]
    a = ca[0::2] + 1j * ca[1::2]
    c_l = (2.0 * np.arange(nl + 1) + 1.0) / 2.0 * a
    nuc_c = np.zeros(2 * nl + 1)
    nuc_c[0::2] = (4.0 * np.arange(nl + 1) + 1.0) / 2.0 * b
    nuclear = np.polynomial.legendre.legval(mu, nuc_c)
    sign = (-1.0) ** np.arange(nl + 1)
    e1 = np.exp(1j * eta * np.log((1.0 - mu) / 2.0))
    e2 = np.exp(1j * eta * np.log((1.0 + mu) / 2.0))
    s = (1.0 + mu) * e1 * np.polynomial.legendre.legval(mu, c_l) + (1.0 - mu) * e2 * np.polynomial.legendre.legval(mu, sign * c_l)
    interference = -2.0 * eta / (1.0 - mu**2) * np.real(s)
    return np.asarray(sigma_coulomb(mu, k, eta) + interference + nuclear)


def ratio_at(sec: dict[str, Any], e_ev: float, mu: np.ndarray) -> np.ndarray:
    """``sigma_e / sigma_c``: the ratios of the two bracketing nodes (each at its own k, eta) interpolated linearly in E."""
    recs = sec["records"]
    en = np.array([r["e_ev"] for r in recs])
    if not en[0] <= e_ev <= en[-1]:
        raise ValueError(f"E={e_ev} eV outside the tabulated range")
    j = int(np.clip(np.searchsorted(en, e_ev, side="right") - 1, 0, en.size - 2))
    out = []
    for rec in (recs[j], recs[j + 1]):
        k, eta = coulomb_k_eta(rec["e_ev"], sec["awr"], sec["awi"])
        out.append(sigma_e(mu, rec, k, eta) / sigma_coulomb(mu, k, eta))
    f = (e_ev - en[j]) / (en[j + 1] - en[j])
    return np.asarray((1.0 - f) * out[0] + f * out[1])


def sigma_ni_barn(sec: dict[str, Any], e_mev: float, n: int | None = None) -> float:
    """``2 pi int_0^0.96 (sigma_e - sigma_c) dmu`` in barn by Gauss-Legendre quadrature (half sphere)."""
    x, w = np.polynomial.legendre.leggauss(n or N_GAUSS)
    mu, wt = 0.5 * MU_CUT * (x + 1.0), 0.5 * MU_CUT * w
    k, eta = coulomb_k_eta(e_mev * 1e6, sec["awr"], sec["awi"])
    ni = sigma_coulomb(mu, k, eta) * (ratio_at(sec, e_mev * 1e6, mu) - 1.0)
    return float(2.0 * math.pi * np.sum(wt * ni))


# ---------------------------------------------------------------------------------------------------------------------------------
# the check
# ---------------------------------------------------------------------------------------------------------------------------------
def dip_test(energies: np.ndarray, level: np.ndarray, margin: float = DIP_MARGIN) -> dict[str, Any]:
    """Interior minimum more than ``margin`` below both endpoint values (ratios only)."""
    i = int(np.argmin(level[1:-1])) + 1
    lo, hi = float(level[i] / level[0] - 1.0), float(level[i] / level[-1] - 1.0)
    return {"persists": bool(level[i] < (1.0 - margin) * level[0] and level[i] < (1.0 - margin) * level[-1]), "min_energy_mev": float(energies[i]),
            "min_over_lower_endpoint_minus_1": lo, "min_over_upper_endpoint_minus_1": hi, "margin": margin}  # fmt: skip


def compare(sec: dict[str, Any], grid_e: np.ndarray, table_sigma_barn: np.ndarray) -> dict[str, Any]:
    """Reconstruction against the table at the nodes and node midpoints of [100, 150] MeV, the dip test on the 1 MeV grid, the reading."""
    idx = [k for k, e in enumerate(grid_e) if E_LO - 1e-9 <= e <= E_HI + 1e-9]
    nodes: dict[str, Any] = {}
    mids: dict[str, Any] = {}
    worst = 0.0
    absolutes: list[float] = []
    for pos, k in enumerate(idx):
        e = float(grid_e[k])
        rec = sigma_ni_barn(sec, e)
        absolutes += [rec, float(table_sigma_barn[k])]
        r = rec / float(table_sigma_barn[k])
        nodes[f"{e:g}"] = {"rec_over_table": r, "within": bool(abs(r - 1.0) <= TOL_MATCH)}
        worst = max(worst, abs(r - 1.0))
        if pos + 1 < len(idx):
            k2 = idx[pos + 1]
            em = 0.5 * (e + float(grid_e[k2]))
            recm = sigma_ni_barn(sec, em)
            tm = 0.5 * (float(table_sigma_barn[k]) + float(table_sigma_barn[k2]))
            absolutes += [recm, tm]
            rm = recm / tm
            mids[f"{em:g}"] = {"rec_over_table": rm, "within": bool(abs(rm - 1.0) <= TOL_MATCH)}
            worst = max(worst, abs(rm - 1.0))
    e1 = np.arange(E_LO, E_HI + 0.5, 1.0)
    lev = np.array([sigma_ni_barn(sec, float(e)) for e in e1])
    absolutes += [float(x) for x in lev]
    dip = dip_test(e1, lev)
    tab = {float(grid_e[k]): float(table_sigma_barn[k]) for k in idx}
    ref_rec, ref_tab = float(lev[int(np.argmin(np.abs(e1 - 125.0)))]), tab.get(125.0)
    named = {f"{e:g}": {"rec_over_table": nodes[f"{e:g}"]["rec_over_table"], "within": nodes[f"{e:g}"]["within"]} for e in NAMED_NODES if f"{e:g}" in nodes}
    levels = {f"{e:g}": {"rec": float(lev[int(np.argmin(np.abs(e1 - e)))]) / ref_rec, "table": (tab[e] / ref_tab) if ref_tab and e in tab else None} for e in NAMED_NODES}
    all_within = bool(worst <= TOL_MATCH)
    reading = "evaluation feature at the LA150 top nodes (the independent reconstruction reproduces the step)" if (all_within and dip["persists"]) else "builder defect route (mismatch beyond 1e-3 or no dip)"
    return {"named_nodes": named, "nodes": nodes, "midpoints": mids, "match": {"max_abs_rel_diff": worst, "tolerance": TOL_MATCH, "all_within": all_within},
            "dip": dip, "levels_over_125mev": levels, "reading": reading, "_absolutes": absolutes}  # fmt: skip


def simpson_uniform(y: np.ndarray, x: np.ndarray) -> float:
    """Composite Simpson rule on a uniform grid with an odd number of points (own implementation)."""
    h = (x[-1] - x[0]) / (x.size - 1)
    return float(h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum()))


def _ni(sec: dict[str, Any], e_mev: float, mu: np.ndarray) -> np.ndarray:
    k, eta = coulomb_k_eta(e_mev * 1e6, sec["awr"], sec["awi"])
    return np.asarray(sigma_coulomb(mu, k, eta) * (ratio_at(sec, e_mev * 1e6, mu) - 1.0))


def quadrature_convergence(sec: dict[str, Any], nodes: tuple[float, ...] = (100.0, 125.0, 150.0)) -> dict[str, Any]:
    """The quadrature error of the reconstruction: GL 1024 (the production rule), GL 2048, Simpson 4001 (the builder's grid) and Simpson 20001, as relative
    differences to GL 2048 at the named nodes; acceptance (pre-declared): every difference below 1e-6. Ratios only."""
    out: dict[str, Any] = {}
    worst = 0.0
    absolutes: list[float] = []
    for e in nodes:
        def gl(n: int) -> float:
            x, w = np.polynomial.legendre.leggauss(n)
            mu, wt = 0.5 * MU_CUT * (x + 1.0), 0.5 * MU_CUT * w
            return float(2.0 * math.pi * np.sum(wt * _ni(sec, e, mu)))

        def sp(n: int) -> float:
            mu = np.linspace(0.0, MU_CUT, n)
            return float(2.0 * math.pi * simpson_uniform(_ni(sec, e, mu), mu))

        ref = gl(2048)
        vals = {"gl1024": gl(N_GAUSS), "simpson4001": sp(4001), "simpson20001": sp(20001)}
        absolutes += [ref, *vals.values()]
        rel = {f"{k}_vs_gl2048": v / ref - 1.0 for k, v in vals.items()}
        worst = max(worst, *(abs(v) for v in rel.values()))
        out[f"{e:g}"] = rel
    return {"nodes": out, "max_abs_rel_diff": worst, "acceptance": 1e-6, "accepted": bool(worst < 1e-6), "_absolutes": absolutes}


def carry_up_check(sec: dict[str, Any], grid_e: np.ndarray, table_sigma_barn: np.ndarray, level_150_barn: float | None = None) -> dict[str, Any]:
    """Carry-up rule above 150 MeV: sigma(E) = sigma_rec(150) S(E)/S(150) with S = ``bgg.pp_bgg_mb`` (the table's own rule), against the pinned nodes in
    (150, 250] MeV; the 175 / 200 / 250 MeV nodes are named. Ratios only; also the carried level on the 1 MeV grid 150-200 MeV relative to 150 MeV."""
    from ionmc.nuclear import bgg  # Geant4 p-p formula (S(E)); not LAW=5 / ENDF code

    rec150 = level_150_barn if level_150_barn is not None else sigma_ni_barn(sec, 150.0)
    s150 = bgg.pp_bgg_mb(150.0)
    nodes: dict[str, Any] = {}
    worst = 0.0
    absolutes = [rec150]
    for k, e in enumerate(grid_e):
        if 150.0 + 1e-9 < e <= 250.0 + 1e-9:
            carried = rec150 * bgg.pp_bgg_mb(float(e)) / s150
            r = carried / float(table_sigma_barn[k])
            absolutes += [carried, float(table_sigma_barn[k])]
            nodes[f"{float(e):g}"] = {"carried_over_table": r, "within": bool(abs(r - 1.0) <= TOL_MATCH)}
            worst = max(worst, abs(r - 1.0))
    named = {f"{e:g}": nodes[f"{e:g}"] for e in (175.0, 200.0, 250.0) if f"{e:g}" in nodes}
    grid = {f"{e:g}": bgg.pp_bgg_mb(float(e)) / s150 for e in np.arange(150.0, 200.5, 1.0)}
    return {"named_nodes": named, "nodes": nodes, "match": {"max_abs_rel_diff": worst, "tolerance": TOL_MATCH, "all_within": bool(nodes) and worst <= TOL_MATCH},
            "grid_150_200_over_150": grid, "rule": "sigma(E) = sigma_rec(150 MeV) x S(E)/S(150 MeV), S = the table's own bgg.pp_bgg_mb", "_absolutes": absolutes}  # fmt: skip


def leak_violations(proj: Any, absolutes: list[float]) -> list[str]:
    """Projected numbers equal (1e-9 relative) to an absolute cross section (barn or mb), and cross-section key names."""
    from .de7_m2prime_record import _numbers

    ab = [x * f for x in absolutes for f in (1.0, 1e3, 1e-3)]
    exempt = {"tolerance", "margin", "min_energy_mev"}
    bad = [f"{p}: equals an absolute cross section" for p, x in _numbers(proj) if p.rsplit("/", 1)[-1] not in exempt and x not in (0.1, 0.01, 1.0, 1e-3)
           and any(a != 0.0 and abs(x - a) <= 1e-9 * abs(a) for a in ab)]  # fmt: skip
    return bad + c.git_safety_violations(proj)


def locate_raw(cache_dir: Path) -> tuple[Path, dict[str, Any]]:
    """The cache object of the LA150 zip the pinned build used (sha256 from the pinned sidecar), verified against that hash."""
    side = json.loads((Path(cache_dir) / "derived" / f"elastic-proton-{c.TABLE_ID}.json").read_text())
    sha = side["sources"][SOURCE_KEY]["sha256"]
    p = Path(cache_dir) / "objects" / sha
    if c.sha256_file(p) != sha:
        raise SystemExit(f"{p}: sha256 differs from the pinned source {sha}")
    npz_sha = c.sha256_file(Path(cache_dir) / "derived" / f"elastic-proton-{c.TABLE_ID}.npz")
    if npz_sha != side["npz_sha256"]:
        raise SystemExit("the pinned npz sha256 differs from its sidecar")
    return p, {"source_sha256": sha, "table_id": side["table_id"], "npz_sha256": npz_sha}


def assemble(sec: dict[str, Any], grid_e: np.ndarray, table_sigma_barn: np.ndarray, ident: dict[str, Any], member_sha256: str) -> tuple[dict[str, Any], list[float]]:
    """The B0 record (full document) and the list of absolute values computed (for the leak check) from the parsed section, the pinned node values and
    the identity of the inputs (``source_sha256`` of the LA150 zip, ``table_id`` and ``npz_sha256`` of the pinned table)."""
    res = compare(sec, grid_e, table_sigma_barn)
    absolutes = res.pop("_absolutes")
    carry = carry_up_check(sec, grid_e, table_sigma_barn)
    absolutes += carry.pop("_absolutes")
    res["carry_up_check"] = carry
    conv = quadrature_convergence(sec)
    absolutes += conv.pop("_absolutes")
    res["quadrature_convergence"] = conv
    doc = {"banner": c.BANNER, "disclosure": DISCLOSURE, **res,
           "provenance": {**ident, "member": MEMBER, "member_sha256": member_sha256, "code_sha256": c.sha256_file(Path(__file__)), "quadrature": f"Gauss-Legendre {N_GAUSS}",
                          **c.repo_state(), "pre_run_looks": PRE_RUN_LOOKS}}  # fmt: skip
    return doc, absolutes


def _write(doc: dict[str, Any], absolutes: list[float], inputs: list[Path], out: Path) -> None:
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), GIT_ALLOW)
    bad = leak_violations(proj, absolutes)
    if bad:
        raise SystemExit(f"leak check failed: {bad[:6]}")
    hdr = c.header("D-E9 B0: independent reconstruction of the pinned p-p level, 100-150 MeV", "python -m validation.scripts.elastic_attribution.de9_pp_level", inputs)
    c.write_outputs(out, "de9-b0.json", doc, "de9-b0-report.md", hdr, render_report, GIT_ALLOW)


def run(cache_dir: Path | None = None, out: Path | None = None, attempts: Path | None = None) -> dict[str, Any]:
    """The whole B0 run, once only. The projected record ``<out>/git/de9-b0.json`` is the target: an existing target is refused before any work, the attempt
    record (``de8-attempts.json``) is written BEFORE the inputs are located and parsed, a failure anywhere (parsing, calculation, leak check, writing) is recorded as
    ``failed`` (two attempts without a document are the bound), and the sha256 of the published projection is recorded for the consumer to authenticate it."""
    out = Path(out or c.OUT / STEP_DIR)
    target = out / "git" / "de9-b0.json"
    att = c.begin_attempt(target, "de9_pp_level.run", {"quadrature": N_GAUSS}, attempts)
    try:
        cd = Path(cache_dir or c.CACHE_DIR)
        zp, ident = locate_raw(cd)
        raw = zipfile.ZipFile(zp).read(MEMBER)
        sec = parse_h1_mt2(raw.decode("ascii"))
        npz = np.load(cd / "derived" / f"elastic-proton-{c.TABLE_ID}.npz", allow_pickle=False)
        doc, absolutes = assemble(sec, npz["grid_e_mev"], npz["sigma_barn"][0], ident, hashlib.sha256(raw).hexdigest())
        _write(doc, absolutes, [zp], out)
    except BaseException:
        c.end_attempt(att, "failed", attempts)
        raise
    c.end_attempt(att, "published-valid", attempts, target)
    return doc


def render_report(doc: dict[str, Any]) -> str:
    m, d = doc["match"], doc["dip"]
    out = ["", doc["disclosure"], "", f"Match (relative difference <= {m['tolerance']:g} at every node and midpoint of 100-150 MeV): **{m['all_within']}** (max |rec/table - 1| = {m['max_abs_rel_diff']:.2e}).", "",
           f"Dip (interior minimum more than {d['margin']:.0%} below both endpoints): **{d['persists']}** (minimum at {d['min_energy_mev']:g} MeV; {d['min_over_lower_endpoint_minus_1']:+.4f} against 100 MeV, {d['min_over_upper_endpoint_minus_1']:+.4f} against 150 MeV).", "",
           f"Reading (pre-declared): **{doc['reading']}**.", "", "| node (MeV) | reconstruction / table | within 1e-3 |", "|---|---|---|"]  # fmt: skip
    for e, r in doc["named_nodes"].items():
        out.append(f"| {e} | {r['rec_over_table']:.6f} | {r['within']} |")
    qc = doc.get("quadrature_convergence")
    if qc:
        out += ["", f"Quadrature convergence (GL 1024 / Simpson 4001 / Simpson 20001 against GL 2048 at 100 / 125 / 150 MeV): max |relative difference| {qc['max_abs_rel_diff']:.2e}; accepted (< {qc['acceptance']:g}): **{qc['accepted']}**."]
    cu = doc.get("carry_up_check")
    if cu:
        out += ["", f"Carry-up check above 150 MeV (rule: {cu['rule']}): all nodes within {cu['match']['tolerance']:g}: **{cu['match']['all_within']}** (max |carried/table - 1| = {cu['match']['max_abs_rel_diff']:.2e}).", "",
                "| node (MeV) | carried / table | within |", "|---|---|---|"]
        out += [f"| {e} | {r['carried_over_table']:.6f} | {r['within']} |" for e, r in cu.get("named_nodes", {}).items()]
    out += ["", "Level relative to the 125 MeV node:", "", "| node (MeV) | reconstruction | table |", "|---|---|---|"]
    for e, r in doc["levels_over_125mev"].items():
        out.append(f"| {e} | {r['rec']:.4f} | {'-' if r.get('table') is None else format(r['table'], '.4f')} |")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    run()
