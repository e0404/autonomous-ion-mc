"""V3-008 P1c: attribution of the Pb residual of the Moliere-Hanson reference and the
compound/water theta_Hanson check (freeze condition FC-1), as a committed, reproducible module.

Replaces the former scratch scripts (runtime monkeypatches of the reference, unhashed inputs).
The committed reference ``moliere_hanson.py`` is used read-only; the variants that the scratch
scripts obtained by monkeypatching (Fano -u_in, Delta b rule, quadrature nodes) are re-stated
here as explicit arguments of :func:`slab_theta`, which equals ``MH.theta0_hanson`` at its
defaults (tested).

Inputs (all hash-pinned; read at run time, never from the network):

* Gottschalk, Med. Phys. 37 (2010) 352, arXiv:0908.1413 LaTeX source (sha256 67fb1478...): the
  published theta_Hanson tables ``tbl:theta0Single`` (Be/Al/Cu/Pb at 158.6 MeV, rho R1 in each
  header), ``tbl:B`` and ``tbl:modStack``. The paper states that ICRU 49 (MIXED.RET) ranges with a
  cubic spline in ln R(ln T) were used.
* Makarova, Gottschalk, Sauerwein, arXiv:1610.01279 PDF (sha256 934fd316...): Table 1 (Go93
  measured theta_0 and the LOOKUP theta_0(Hanson)); text extracted with ``mutool`` (or
  ``pdftotext``) and parsed by :func:`parse_makarova_table1`.
* NIST PSTAR (ICRU 49) CSDA ranges, the cached POST responses (public domain; only ratios are
  derived and stored): Be, Al, Cu, Pb, polystyrene. Parsed by :func:`parse_pstar`.

The ICRU 49 range-energy path is a CSDA range function R(E) turned into an energy-at-depth
function (:func:`path_from_range_function`): natural cubic spline of ln R(ln T) through
Gottschalk's (E, residual range) nodes ("S"), our table scaled by kappa(E) through the same
nodes ("K"), or a spline through the PSTAR table.

The Pb clearing is an ATTRIBUTION (the failure of the committed Pb test on ionmc Bethe ranges is
an input difference in the range table), not a pass of that test.

Run: ``python validation/scripts/ion_em/p1c_attribution.py [--out FILE]`` (aggregates only, no
PSTAR values).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tarfile
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.materials import AIR, ALUMINIUM, BERYLLIUM, COPPER, LEAD, WATER, Material
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource

REPO = Path(__file__).resolve().parents[3]

_SPEC = importlib.util.spec_from_file_location(
    "moliere_hanson", REPO / "validation/scripts/ion_em/moliere_hanson.py"
)
assert _SPEC is not None and _SPEC.loader is not None
MH = importlib.util.module_from_spec(_SPEC)
sys.modules.setdefault("moliere_hanson", MH)
_SPEC.loader.exec_module(MH)

Element = tuple[int, float, float]
Path1 = Callable[[NDArray[np.float64]], NDArray[np.float64]]

E1 = 158.6
"""Incident energy [MeV] of tbl:theta0Single."""
MP = PROTON.mass_mev

GOTTSCHALK_SHA256 = "67fb1478e51534064f0e5363b1ee4160da2a37a600f51897b56dd74ea22ba5ed"
MAKAROVA_PDF_SHA256 = "934fd31644eb97734a0f4c9a756a40d49012d3c126ef2b4391e088be7a092aa4"
PSTAR_SHA256 = {  # NIST PSTAR POST responses (plain text), keyed by matno
    "004": "0fcbb24d9b7f07af7898ac3954802d08010d22339e195ee7cd18fe2638380757",  # beryllium
    "013": "3fe6459ec4abde6052bcb63f427ea856bbc35dd8506014191145b4eb9253db5e",  # aluminum
    "029": "986e8577db37120c405f6686394c1250c3e7f1d6af7e88547b964de2e99e0063",  # copper
    "082": "11dbcd1b014d926ef3e5b8bdd2bf9e539c153b5a4c8d51a62d81c34f5b0a05a7",  # lead
    "219": "8974efe7242af76f1df4fb696a7b957c696f6ca7de711e01c3c68a2f40c04040",  # polycarbonate
    "223": "1e13828ffc3dee1e5748d2d922e5ecf2c583380903bb787f84fbe2886dc549e9",  # PMMA (lucite)
    "226": "dc65ded99e0cb9b455482404642054c628f0639fd75ed8ec5ebc0e5e68f981c1",  # polystyrene
}
PSTAR_NAMES = {
    "004": "BERYLLIUM",
    "013": "ALUMINUM",
    "029": "COPPER",
    "082": "LEAD",
    "219": "POLYCARBONATE (MAKROLON, LEXAN)",
    "223": "POLYMETHYL METHACRALATE (LUCITE, PERSPEX, PLEXIGLASS)",
    "226": "POLYSTYRENE",
}
MATERIALS = {"Be": BERYLLIUM, "Al": ALUMINIUM, "Cu": COPPER, "Pb": LEAD}
PSTAR_MATNO = {"Be": "004", "Al": "013", "Cu": "029", "Pb": "082"}


class InputUnavailable(RuntimeError):
    """A hash-pinned input is not in the cache or cannot be extracted (tests skip or fail)."""


# --------------------------------------------------------------------------------- input access
def cache_dirs() -> list[Path]:
    """Directories searched for sha256-named objects: ``$IONMC_CACHE_DIR[/objects]`` first, then
    the repository ``.ionmc-cache/reference`` and ``.ionmc-cache/ionmc-data/objects``."""
    dirs: list[Path] = []
    env = os.environ.get("IONMC_CACHE_DIR")
    if env:
        dirs += [Path(env) / "objects", Path(env)]
    dirs += [
        REPO / ".ionmc-cache" / "reference",
        REPO / ".ionmc-cache" / "ionmc-data" / "objects",
    ]
    return dirs


def cached_object(sha256: str) -> bytes:
    """Bytes of the cached object named ``sha256`` (content hash verified)."""
    for d in cache_dirs():
        p = d / sha256
        if p.is_file():
            raw = p.read_bytes()
            if hashlib.sha256(raw).hexdigest() != sha256:
                raise ValueError(f"{p} does not have sha256 {sha256}")
            return raw
    raise InputUnavailable(f"object sha256 {sha256[:8]}... not found in {[str(d) for d in cache_dirs()]}")


def gottschalk_tex() -> str:
    """LaTeX source of arXiv:0908.1413 (ScatPowerV2.tex) from the cached e-print."""
    raw = cached_object(GOTTSCHALK_SHA256)
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as tar:
        member = next(m for m in tar.getmembers() if m.name.endswith("ScatPowerV2.tex"))
        handle = tar.extractfile(member)
        assert handle is not None
        return handle.read().decode("utf-8", errors="replace")


def makarova_pdf_text() -> str:
    """Text of the cached arXiv:1610.01279 PDF via ``mutool draw -F txt`` or ``pdftotext``."""
    raw = cached_object(MAKAROVA_PDF_SHA256)
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        pdf = Path(td) / "makarova.pdf"
        pdf.write_bytes(raw)
        if shutil.which("mutool"):
            cmd = ["mutool", "draw", "-F", "txt", "-o", "-", str(pdf)]
        elif shutil.which("pdftotext"):
            cmd = ["pdftotext", str(pdf), "-"]
        else:
            raise InputUnavailable("neither mutool nor pdftotext is available to read the PDF")
        out = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=120)
    return out.stdout


# ------------------------------------------------------------------------------------- parsers
_NUM = r"(-?\d+\.\d+)"
_ROW11 = re.compile(r"^\s*" + r"\s*&\s*".join([_NUM] * 11) + r"\s*(?:\\\\)?\s*$")


def theta_table(t: str) -> dict[str, dict[str, Any]]:
    """``tbl:theta0Single``: ``{mat: {'rhoR1': .., 'rows': [(x/R1, rho x, E, theta_H mrad)]}}``."""
    start = t.index("%%% Theta0 table:")
    end = t.index("tbl:theta0Single", start)
    out: dict[str, dict[str, Any]] = {}
    mat = None
    for line in t[start:end].splitlines():
        h = re.search(r"158\.60\\,MeV p on (\w+), \$\\rho R_1=([\d.]+)\$", line)
        if h:
            mat = h.group(1)
            out[mat] = {"rhoR1": float(h.group(2)), "rows": []}
            continue
        m = _ROW11.match(line)
        if m and mat:
            v = [float(x) for x in m.groups()]
            out[mat]["rows"].append((v[0], v[1], v[2], v[3]))
    return out


MAKAROVA_MATERIALS = (
    "Beryllium", "Polystyrene", "Carbon", "Lexan", "Nylon", "Lucite", "Teflon",
    "Aluminum", "Copper", "Zinc", "Brass", "Tin", "Lead", "Uranium",
)  # fmt: skip


def parse_makarova_table1(text: str) -> dict[str, list[tuple[float, float, float, float, float]]]:
    """Table 1 of arXiv:1610.01279 from extracted PDF text: per material rows
    ``(thickness g/cm2, theta0 Go93 exptl, theta0 Hanson LOOKUP, Urban, Wentzel)`` [mrad]."""
    a = text.index("Table 1: Selected results")
    b = text.index("Table 2: Material properties", a)
    cur: str | None = None
    vals: dict[str, list[float]] = {}
    for tok in text[a:b].split():
        if tok in MAKAROVA_MATERIALS:
            cur = tok
            vals.setdefault(cur, [])
        elif cur and re.fullmatch(r"\d+\.\d+", tok):
            vals[cur].append(float(tok))
    out: dict[str, list[tuple[float, float, float, float, float]]] = {}
    for k, v in vals.items():
        if len(v) % 5:
            raise ValueError(f"Table 1 {k}: {len(v)} numbers is not a multiple of 5")
        out[k] = [tuple(v[i : i + 5]) for i in range(0, len(v), 5)]  # type: ignore[misc]
    return out


def parse_pstar(raw: bytes | str) -> tuple[str, NDArray[np.float64], NDArray[np.float64]]:
    """NIST PSTAR POST response: ``(material name, kinetic energy [MeV], CSDA range [g/cm2])``."""
    text = raw.decode("utf-8", errors="replace") if isinstance(raw, bytes) else raw
    lines = text.splitlines()
    name = next(ln.strip() for ln in lines[1:] if ln.strip())
    e: list[float] = []
    r: list[float] = []
    for ln in lines:
        v = ln.split()
        if len(v) == 7 and all(re.fullmatch(r"\d\.\d+E[+-]\d+|\d\.\d+", x) for x in v):
            e.append(float(v[0]))
            r.append(float(v[4]))
    return name, np.array(e), np.array(r)


def pstar_data(matno: str) -> tuple[str, NDArray[np.float64], NDArray[np.float64]]:
    return parse_pstar(cached_object(PSTAR_SHA256[matno]))


# --------------------------------------------------------------------------- range-energy paths
def natural_cubic_spline(x: Any, y: Any) -> Callable[[Any], NDArray[np.float64]]:
    """Natural cubic spline through (x, y); returns a vectorised evaluator."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    h = np.diff(x)
    a = np.zeros((n, n))
    rhs = np.zeros(n)
    a[0, 0] = a[-1, -1] = 1.0
    for i in range(1, n - 1):
        a[i, i - 1] = h[i - 1]
        a[i, i] = 2 * (h[i - 1] + h[i])
        a[i, i + 1] = h[i]
        rhs[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    c = np.linalg.solve(a, rhs)
    b = (y[1:] - y[:-1]) / h - h * (2 * c[:-1] + c[1:]) / 3
    d = (c[1:] - c[:-1]) / (3 * h)

    def f(xq: Any) -> NDArray[np.float64]:
        xq = np.asarray(xq, float)
        i = np.clip(np.searchsorted(x, xq) - 1, 0, n - 2)
        dx = xq - x[i]
        return np.asarray(y[i] + b[i] * dx + c[i] * dx**2 + d[i] * dx**3)

    return f


def log_spline_range(e_mev: Any, range_g_cm2: Any) -> Callable[[Any], NDArray[np.float64]]:
    """R(E) as the natural cubic spline of ln R(ln T) through the nodes (the interpolation
    Gottschalk states for ICRU 49 / MIXED.RET)."""
    s = natural_cubic_spline(np.log(e_mev), np.log(range_g_cm2))
    return lambda q: np.exp(s(np.log(np.asarray(q, float))))


def path_from_range_function(
    range_of_e: Callable[[Any], NDArray[np.float64]], e_lo: float, e_hi: float, r1: float,
    npts: int = 20000,
) -> Path1:  # fmt: skip
    """Energy at depth t from a monotone range function R(E) [g/cm2] by dense log-log inversion;
    the beam has residual range ``r1 - t``."""
    eg = np.exp(np.linspace(math.log(e_lo), math.log(e_hi), npts))
    rg = range_of_e(eg)
    assert np.all(np.diff(rg) > 0), "range function not monotone"

    def energy_at(t: Any) -> NDArray[np.float64]:
        rem = r1 - np.asarray(t, float)
        return np.asarray(np.exp(np.interp(np.log(rem), np.log(rg), np.log(eg))))

    return energy_at


def ours_path(material: Material, e0: float = E1, src: BetheStoppingSource | None = None) -> Any:
    """``(table, R1, energy_at)`` on the ionmc analytic Bethe CSDA range (the committed P1 path)."""
    src = src or BetheStoppingSource()
    table = src.table(material, PROTON)
    r1 = float(table.range_at(e0))

    def energy_at(t: Any) -> NDArray[np.float64]:
        rem = np.maximum(r1 - t, float(table.csda_range_g_cm2[0]))
        return np.asarray(table.energy_from_range(rem))

    return table, r1, energy_at


def gottschalk_nodes(rows: Sequence[Any], r1g: float, use: Sequence[float]) -> Any:
    """Gottschalk's (E, residual range) nodes sorted by energy, incl. the entrance node."""
    ne = np.array([E1] + [eg for f, _x, eg, _t in rows if f in use])
    nr = np.array([r1g] + [r1g - x for f, x, _eg, _t in rows if f in use])
    o = np.argsort(ne)
    return ne[o], nr[o]


def spline_path(rows: Sequence[Any], r1g: float, use: Sequence[float]) -> Path1:
    """Variant S: natural cubic spline of ln R(ln T) through Gottschalk's own nodes."""
    ne, nr = gottschalk_nodes(rows, r1g, use)
    return path_from_range_function(log_spline_range(ne, nr), ne[0] * 0.999, E1, r1g)


def kappa_path(table: Any, rows: Sequence[Any], r1g: float, use: Sequence[float]) -> Any:
    """Variant K: the ionmc table range times kappa(E) = R_G/R_ours interpolated (in ln E) through
    Gottschalk's nodes. Returns ``(path, ln kappa at the nodes)``."""
    ne, nr = gottschalk_nodes(rows, r1g, use)
    lk = np.log(nr / np.asarray(table.range_at(ne)))

    def rg(e: Any) -> NDArray[np.float64]:
        return np.asarray(table.range_at(e)) * np.exp(np.interp(np.log(e), np.log(ne), lk))

    return path_from_range_function(rg, ne[0] * 0.999, E1, r1g), lk


def pstar_range_function(matno: str) -> tuple[Callable[[Any], NDArray[np.float64]], int]:
    """Natural-spline ln R(ln T) through the cached PSTAR CSDA table; also its point count."""
    _name, e, r = pstar_data(matno)
    return log_spline_range(e, r), len(e)


# -------------------------------------------------------------------- Moliere-Hanson evaluation
def elements(m: Material) -> tuple[Element, ...]:
    return tuple((e.Z, e.A_g_mol, w) for e, w in m._fractions)


def _fano_d(zt: int, b2: Any, neg_u_default: float) -> NDArray[np.float64]:
    neg_u = MH.U_IN_NEG.get(zt, neg_u_default)
    b2a = np.asarray(b2, dtype=np.float64)
    out: NDArray[np.float64] = (
        np.log(MH.FANO_CONST * zt ** (-4.0 / 3.0) * b2a / (1.0 - b2a)) + neg_u - 0.5 * b2a
    )
    return out


def _uniform_nodes(t: float, n: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    gx, gw = np.polynomial.legendre.leggauss(n)
    return 0.5 * t * (gx + 1.0), 0.5 * t * gw


def rule_electron_weighted(els: Sequence[Element], b2: Any) -> NDArray[np.float64]:
    """Alternative compound Fano rule: sum (w Z/A)(D/Z) / sum w Z/A."""
    num = sum(w * z / a * MH.fano_d(z, b2) / z for z, a, w in els)
    return np.asarray(num / sum(w * z / a for z, a, w in els))


def rule_mass_weighted(els: Sequence[Element], b2: Any) -> NDArray[np.float64]:
    """Alternative compound Fano rule: sum w (D/Z)."""
    return np.asarray(sum(w * MH.fano_d(z, b2) / z for z, a, w in els))


def slab_theta(
    elems: Sequence[Element],
    t: float,
    energy_at: Path1,
    n: int = 400,
    electrons: str = "fano",
    *,
    neg_u_in_default: float | None = None,
    delta_b: Callable[[Sequence[Element], Any], Any] | None = None,
    nodes: Callable[[float, int], tuple[NDArray[np.float64], NDArray[np.float64]]] | None = None,
) -> dict[str, float]:
    """``MH.moliere_slab`` (z = 1, proton) with explicit, local variants: ``neg_u_in_default``
    (-u_in for Z > 1), ``delta_b(els, b2)`` (the compound Fano rule; a function returning zeros
    switches the Fano term off) and ``nodes(t, n)`` (the quadrature). At the defaults the result
    equals ``MH.moliere_slab`` exactly (same operations; tested)."""
    neg_u = MH.U_IN_NEG_DEFAULT if neg_u_in_default is None else neg_u_in_default
    node_x, weights = (nodes or MH._gauss_nodes)(t, n)
    energies = energy_at(node_x)
    inv_pv2 = 1.0 / MH.pv_mev(energies, MP) ** 2
    integral = float(np.sum(weights * inv_pv2))
    MH._check_electrons(electrons)
    zsum = MH.z2_sum(elems) if electrons == "fano" else MH.zz1_sum(elems)
    chi_c2 = MH.K_CHI_C * 1.0 * 1.0 * zsum * integral
    lna = (
        float(np.sum(weights * inv_pv2 * MH.ln_chi_a2(1.0, elems, energies, MP, electrons)))
        / integral
    )
    db_mean = 0.0
    if electrons == "fano":
        b2 = MH.beta2(energies, MP)
        if delta_b is None:
            num = np.zeros_like(np.asarray(b2, dtype=np.float64))
            for zt, a, w in elems:
                num += w * zt / a * _fano_d(zt, b2, neg_u)
            db = num / MH.z2_sum(elems)
        else:
            db = delta_b(elems, b2)
        db_mean = float(np.sum(weights * inv_pv2 * db)) / integral
    b = math.log(chi_c2) - lna - math.log(MH.OMEGA_FACTOR) + db_mean
    big_b = MH.solve_b(b)
    return {
        "chi_c2": chi_c2, "ln_chi_a2": lna, "delta_b": db_mean, "b": b, "B": big_b,
        "theta0": math.sqrt(chi_c2 * (big_b - MH.HANSON_SHIFT) / 2.0),
    }  # fmt: skip


def theta(elems: Sequence[Element], t: float, energy_at: Path1, **kw: Any) -> float:
    """Projected Hanson width [rad] of a slab (see :func:`slab_theta`)."""
    return slab_theta(elems, t, energy_at, **kw)["theta0"]


def dev(th_rad: float, th_h_mrad: float) -> float:
    """theta/theta_Hanson - 1 (theta in rad, theta_Hanson in mrad)."""
    return 1e3 * th_rad / th_h_mrad - 1.0


def _mx(d: dict[Any, Any]) -> float:
    return max(abs(v) for v in d.values() if v is not None)


# --------------------------------------------------------------------- (1) Pb attribution
def pb_attribution(
    name: str,
    mat: Material,
    tbl: dict[str, Any],
    quad_ns: Sequence[int] = (50, 100, 200, 400, 800, 1600, 3200),
    dense_table: bool = True,
) -> dict[str, Any]:
    """Residuals of the committed reference against theta_Hanson of ``tbl:theta0Single`` for one
    element on the ionmc path (baseline) and on the ICRU 49 node paths, plus the Fano and
    quadrature controls (all relative deviations; x = (x/R1) * R1 of the respective path)."""
    el = elements(mat)
    g = tbl[name]
    rows = g["rows"]
    r1g = g["rhoR1"]
    table, r1o, e_o = ours_path(mat)
    out: dict[str, Any] = {"rhoR1_ours_over_G_minus_1": r1o / r1g - 1.0}
    out["baseline"] = {f: dev(theta(el, f * r1o, e_o), th) for f, _x, _e, th in rows}
    out["residual_range_ours_at_EG_over_G_minus_1"] = {
        f: float(table.range_at(eg)) / (r1g - x) - 1.0 for f, x, eg, _t in rows
    }
    out["abs_thickness_on_our_path"] = {
        f: (dev(theta(el, x, e_o), th) if x < 0.995 * r1o else None) for f, x, _e, th in rows
    }
    fr_all = [f for f, *_ in rows]
    fr_mid = [f for f in fr_all if f >= 0.1]
    ek_all, lk_all = kappa_path(table, rows, r1g, fr_all)
    ek_mid, _ = kappa_path(table, rows, r1g, fr_mid)
    variants: dict[str, Any] = {}
    for vname, path in (
        ("kappa_all_nodes", ek_all),
        ("kappa_nodes_ge_0.1", ek_mid),
        ("spline_lnR_lnT_nodes_ge_0.1", spline_path(rows, r1g, fr_mid)),
        ("spline_lnR_lnT_all_nodes", spline_path(rows, r1g, fr_all)),
    ):
        variants[vname] = {f: dev(theta(el, f * r1g, path), th) for f, _x, _e, th in rows}
        variants[vname + "_node_E_resid_MeV_max"] = max(
            abs(float(path(np.array([x]))[0]) - eg) for f, x, eg, _t in rows if f >= 0.1
        )
    out["icru49_path"] = variants
    out["ln_kappa_range"] = [float(lk_all.min()), float(lk_all.max())]

    sens = {}
    for f, _x, _e, _t in rows:
        if f >= 0.5:
            h = 1e-3
            tp = theta(el, f * r1o * (1 + h), e_o)
            tm = theta(el, f * r1o * (1 - h), e_o)
            sens[f] = math.log(tp / tm) / math.log((1 + h) / (1 - h))
    out["dlntheta_dlnx_fixed_path"] = sens

    fano: dict[str, Any] = {}
    for neg_u in (3.6, 4.0, 5.0, 6.0):
        fano[f"fano_neg_u_in_{neg_u}"] = {
            f: dev(theta(el, f * r1o, e_o, neg_u_in_default=neg_u), th) for f, _x, _e, th in rows
        }
    fano["moliere_Z2_no_delta_b"] = {
        f: dev(theta(el, f * r1o, e_o, delta_b=lambda _els, b2: np.zeros_like(np.asarray(b2))), th)
        for f, _x, _e, th in rows
    }
    fano["bethe_ZZ1"] = {
        f: dev(theta(el, f * r1o, e_o, electrons="bethe"), th) for f, _x, _e, th in rows
    }
    fano["delta_b_at_EG"] = {
        f: float(MH.fano_delta_b(el, MH.beta2(np.array([eg]), MH.M_P_MEV))[0])
        for f, _x, eg, _t in rows
    }
    fano["slab_b_B"] = {
        f: {k: MH.moliere_slab(1.0, MH.M_P_MEV, el, f * r1o, e_o)[k] for k in ("b", "B", "delta_b")}
        for f, _x, _e, _t in rows
    }
    out["electrons"] = fano

    quad: dict[str, Any] = {}
    for n in quad_ns:
        quad[f"n{n}"] = {f: dev(theta(el, f * r1o, e_o, n=n), th) for f, _x, _e, th in rows}
    quad["uniform_GL_4000"] = {
        f: dev(theta(el, f * r1o, e_o, n=4000, nodes=_uniform_nodes), th)
        for f, _x, _e, th in rows
    }
    if dense_table:
        dense = BetheStoppingSource(points_per_decade=1000)
        _t2, r1d, e_d = ours_path(mat, src=dense)
        quad["table_1000_per_decade"] = {
            f: dev(theta(el, f * r1d, e_d), th) for f, _x, _e, th in rows
        }
    out["numerics"] = quad
    return out


def pstar_element_check(name: str, mat: Material, tbl: dict[str, Any]) -> dict[str, Any]:
    """Element on the PSTAR (ICRU 49) CSDA path (natural ln R(ln T) spline)."""
    rp, npts = pstar_range_function(PSTAR_MATNO[name])
    g = tbl[name]
    r1p = float(rp(E1))
    el = elements(mat)
    path = path_from_range_function(rp, 5.0, E1, r1p)
    _t, r1o, _e = ours_path(mat)
    return {
        "pstar_points": npts,
        "R1_pstar_over_G_minus_1": r1p / g["rhoR1"] - 1.0,
        "R1_ours_over_pstar_minus_1": r1o / r1p - 1.0,
        "resid_range_pstar_at_EG_over_G_minus_1": {
            f: float(rp(eg)) / (g["rhoR1"] - x) - 1.0 for f, x, eg, _ in g["rows"]
        },
        "theta_on_pstar_path_over_thetaH_minus_1": {
            f: dev(theta(el, f * r1p, path), th) for f, _x, _e2, th in g["rows"]
        },
    }


# ------------------------------------------------------------------- (2) compound / water check
ZA = {
    "H": (1, 1.00794), "C": (6, 12.0107), "N": (7, 14.0067), "O": (8, 15.9994),
    "F": (9, 18.9984032), "Be": (4, 9.012182), "Pb": (82, 207.2),
}  # fmt: skip


def els_from_fractions(frac: dict[str, float]) -> tuple[Element, ...]:
    tot = sum(frac.values())
    return tuple((ZA[s][0], ZA[s][1], w / tot) for s, w in frac.items())


COMP: dict[str, tuple[dict[str, float], float]] = {
    "Polystyrene": ({"H": 0.077418, "C": 0.922582}, 68.7),  # G4_POLYSTYRENE (C8H8)
    "Lexan": ({"C": 0.741, "O": 0.185, "H": 0.074}, 73.1),  # Go93 composition (arXiv:1610.01279)
    "Lexan_C16H14O3": ({"H": 0.055491, "C": 0.755751, "O": 0.188758}, 73.1),  # stoichiometric
    "Lucite": ({"C": 0.600, "O": 0.320, "H": 0.081}, 74.0),
    "Nylon": ({"C": 0.549, "O": 0.244, "N": 0.107, "H": 0.100}, 63.9),
    "Carbon": ({"C": 1.0}, 78.0),
    "Beryllium": ({"Be": 1.0}, 63.7),
    "Lead": ({"Pb": 1.0}, 823.0),
}
"""Mass fractions and ICRU 49 mean excitation energies [eV] of the Go93 targets (compositions of
arXiv:1610.01279 Table 3 / ``tab:materials``; Nylon ranges in LOOKUP are Janni 82, here the same
Bethe model)."""
TEFLON = {"C": 0.240183, "F": 0.759817}
"""C2F4 mass fractions (Teflon points are thin: constant energy)."""


def make_material(name: str, frac: dict[str, float], i_ev: float) -> Material:
    tot = sum(frac.values())
    return Material(name, 1.0, {k: v / tot for k, v in frac.items()}, i_ev, None, "p1c")


def _const(t: Any) -> NDArray[np.float64]:
    return np.full_like(np.asarray(t, float), E1)


def _path_of(tab: Any, r1: float, scale: float = 1.0) -> Path1:
    lo = float(tab.csda_range_g_cm2[0])
    return lambda t: np.asarray(
        tab.energy_from_range(np.maximum(r1 - np.asarray(t, float) / scale, lo))
    )


def compound_check(
    makarova: dict[str, list[tuple[float, ...]]], g_tex: str, src: BetheStoppingSource
) -> dict[str, Any]:
    """Compounds against the LOOKUP theta_Hanson of arXiv:1610.01279 Table 1, the tbl:B b check
    (water/air/Be/Cu/Pb) and tbl:modStack rows 1 and 12 of arXiv:0908.1413."""
    res: dict[str, Any] = {}
    mak: dict[str, list[dict[str, Any]]] = {}
    for name, rows in makarova.items():
        keys = [name] + (["Lexan_C16H14O3"] if name == "Lexan" else [])
        for key in keys:
            if key not in COMP and name != "Teflon":
                continue
            for thick, th_exp, th_h, _u, _w in rows:
                entry: dict[str, Any] = {"t_g_cm2": thick}
                if name == "Teflon":
                    if thick > 0.1:
                        continue
                    el = els_from_fractions(TEFLON)
                    eat: Path1 = _const
                    entry["path"] = "constant 158.6 MeV (energy loss < 0.3 MeV)"
                else:
                    frac, i_ev = COMP[key]
                    el = els_from_fractions(frac)
                    tab = src.table(make_material(key, frac, i_ev), PROTON)
                    r1 = float(tab.range_at(E1))
                    eat = _path_of(tab, r1)
                    entry["x_over_R1_ours"] = thick / r1
                    sens = {}
                    for fi in (0.9, 1.1):
                        tb2 = src.table(make_material(key, frac, i_ev * fi), PROTON)
                        r2 = float(tb2.range_at(E1))
                        sens[f"I_x{fi}"] = {
                            "R1_rel": r2 / r1 - 1.0,
                            "theta_rel": theta(el, thick, _path_of(tb2, r2)) / theta(el, thick, eat)
                            - 1.0,
                        }
                    entry["range_sensitivity"] = sens
                th = theta(el, thick, eat)
                entry["ours_over_theta_hanson_LOOKUP_minus_1"] = dev(th, th_h)
                entry["ours_over_measured_Go93_minus_1"] = dev(th, th_exp)
                entry["LOOKUP_hanson_over_measured_minus_1"] = th_h / th_exp - 1.0
                entry["bethe_ZZ1_over_theta_hanson_minus_1"] = dev(
                    theta(el, thick, eat, electrons="bethe"), th_h
                )
                s = MH.moliere_slab(1.0, MP, el, thick, eat)
                entry["B"] = s["B"]
                entry["delta_b"] = s["delta_b"]
                mak.setdefault(key, []).append(entry)
    res["makarova_tab_results"] = mak

    tb = g_tex[
        g_tex.index("\\begin{center}\\begin{tabular}{crrcrr}") : g_tex.index("\\label{tbl:B}")
    ]
    tblb: dict[str, dict[str, float]] = {}
    for line in tb.splitlines():
        m = re.match(r"^(\w+)&([\d.]+)&([\d.]+)&([\d.]+)(\\,cm)?&([\d.]+)&([\d.]+)", line)
        if m:
            rho, x = float(m.group(2)), float(m.group(4))
            t_gcm2 = rho * (x if m.group(5) else x * 1e-4)
            tblb[m.group(1)] = {
                "rho": rho, "rhoR1": float(m.group(3)), "t_g_cm2": t_gcm2,
                "xR1e6": float(m.group(6)), "thetaH_mrad": float(m.group(7)),
            }  # fmt: skip
    el_b = {
        "water": elements(WATER), "air": elements(AIR), "Be": els_from_fractions({"Be": 1.0}),
        "Cu": ((29, 63.546, 1.0),), "Pb": elements(LEAD),
    }  # fmt: skip
    outb: dict[str, Any] = {}
    b_target = 4.5 - math.log(4.5)
    for k, v in tblb.items():
        el = el_b[k]
        s = MH.moliere_slab(1.0, MP, el, v["t_g_cm2"], _const)
        t45 = v["t_g_cm2"] * math.exp(b_target - s["b"])
        s45 = MH.moliere_slab(1.0, MP, el, t45, _const)
        try:
            sb_b = MH.moliere_slab(1.0, MP, el, v["t_g_cm2"], _const, electrons="bethe")["b"]
        except ValueError:
            sb_b = float("nan")  # Bethe variant b <= 1 at this thickness
        outb[k] = {
            "B_ours_at_G_thickness": s["B"],
            "b_ours_minus_b_G": s["b"] - b_target,
            "t_B4.5_ours_over_G_minus_1": t45 / v["t_g_cm2"] - 1.0,
            "check_B_at_t45": s45["B"],
            "theta_ours_over_G_minus_1": dev(s["theta0"], v["thetaH_mrad"]),
            "theta_rounding_half_interval_rel": 0.0005 / v["thetaH_mrad"],
            "thickness_rounding_half_interval_rel": 0.005
            / (v["t_g_cm2"] / v["rho"] * (1 if k == "air" else 1e4)),
            "bethe_b_minus_b_G": sb_b - b_target,
            "delta_b": s["delta_b"],
        }
    res["gottschalk2010_tblB"] = outb

    ms = g_tex[g_tex.index("\\label{tbl:modStack}") - 4000 : g_tex.index("\\label{tbl:modStack}")]
    prow = [
        [float(x) for x in re.split(r"\s*&\s*", ln.strip().rstrip("\\").strip())]
        for ln in ms.splitlines()
        if re.match(r"^\s*\d+&\s*[\d.]+&", ln)
    ]
    mrows = {int(r[0]): r for r in prow}
    mod: dict[str, Any] = {}
    for idx, key, frac, i_ev, thick, e_out in (
        (1, "Pb", {"Pb": 1.0}, 823.0, mrows[1][1], mrows[1][3]),
        (12, "Lexan", COMP["Lexan"][0], 73.1, mrows[12][2], mrows[12][4]),
        (12, "Lexan_C16H14O3", COMP["Lexan_C16H14O3"][0], 73.1, mrows[12][2], mrows[12][4]),
    ):
        th_h = mrows[idx][6]
        el = els_from_fractions(frac)
        tab = src.table(make_material(key, frac, i_ev), PROTON)
        r1 = float(tab.range_at(230.0))
        eat = lambda t, tb=tab, rr=r1: np.asarray(  # noqa: E731
            tb.energy_from_range(rr - np.asarray(t, float))
        )
        e_out_ours = float(eat(np.array([thick]))[0])
        r_out = float(tab.range_at(e_out))
        kappa = thick / (r1 - r_out)  # constant kappa pinning the exit energy to the table
        eat_k = _path_of(tab, r1, kappa)
        mod[f"row{idx}_{key}"] = {
            "x_over_R1_ours": thick / r1,
            "E_out_ours_minus_table_MeV": e_out_ours - e_out,
            "ours_over_thetaH_minus_1": dev(theta(el, thick, eat), th_h),
            "kappa_pinned_minus_1": kappa - 1.0,
            "ours_pinned_over_thetaH_minus_1": dev(theta(el, thick, eat_k), th_h),
        }
    res["gottschalk2010_modStack"] = mod
    return res


def lookup_value(
    makarova: dict[str, list[tuple[float, ...]]], material: str, thickness: str
) -> tuple[float, float]:
    """``(theta_exptl, theta_Hanson LOOKUP)`` [mrad] of the Table 1 row with that thickness."""
    for row in makarova[material]:
        if f"{row[0]:g}" == f"{float(thickness):g}":
            return row[1], row[2]
    raise KeyError((material, thickness))


def thick_on_pstar_path(makarova: dict[str, list[tuple[float, ...]]]) -> dict[str, Any]:
    """Thick Table 1 points (Pb, polystyrene, Be) on the PSTAR path vs LOOKUP."""
    out: dict[str, Any] = {}
    poly = ((1, 1.00794, 0.077418), (6, 12.0107, 0.922582))
    for label, matno, el, t in (
        ("Lead", "082", elements(LEAD), "31.566"),
        ("Polystyrene", "226", poly, "15.751"),
        ("Beryllium", "004", elements(BERYLLIUM), "20.313"),
    ):
        rp, _n = pstar_range_function(matno)
        r1p = float(rp(E1))
        path = path_from_range_function(rp, 5.0, E1, r1p)
        _exp, th_h = lookup_value(makarova, label, t)
        out[label] = {
            "x_over_R1_pstar": float(t) / r1p,
            "ours_on_pstar_path_over_LOOKUP_minus_1": dev(theta(el, float(t), path), th_h),
        }
    return out


def alt_compound_rules(makarova: dict[str, list[tuple[float, ...]]]) -> dict[str, Any]:
    """Thin-slab (constant energy) theta vs LOOKUP for the committed compound Fano rule and the
    electron-/mass-weighted alternatives (discrimination power of the comparison)."""
    alt: dict[str, Any] = {}
    luc = (
        (6, 12.0107, 0.600 / 1.001), (8, 15.9994, 0.320 / 1.001), (1, 1.00794, 0.081 / 1.001),
    )  # fmt: skip
    poly = ((1, 1.00794, 0.077418), (6, 12.0107, 0.922582))
    for label, el, mat, t in (
        ("Polystyrene_0.347", poly, "Polystyrene", "0.347"),
        ("Lucite_0.366", luc, "Lucite", "0.366"),
    ):
        th_h = lookup_value(makarova, mat, t)[1]
        d = {}
        for rname, fn in (
            ("committed", None),
            ("electron_weighted", rule_electron_weighted),
            ("mass_weighted", rule_mass_weighted),
        ):
            d[rname] = dev(theta(el, float(t), _const, delta_b=fn), th_h)
        alt[label] = d
    return alt


def water_b_shift_effect() -> dict[str, Any]:
    """Relative theta change of a b shift of -0.040 (tbl:B water finding) on water slabs."""
    wel = elements(WATER)
    eff: dict[str, Any] = {}
    for e0, f in ((150.0, 0.1), (150.0, 0.5), (150.0, 0.9)):
        _t, r1, eat = ours_path(WATER, e0)
        s = MH.moliere_slab(1.0, MP, wel, f * r1, eat)
        bb = MH.solve_b(s["b"] - 0.040)
        eff[f"water_{e0:.0f}MeV_xR{f}"] = {
            "B": s["B"],
            "theta_rel_change": math.sqrt((bb - 1.2) / (s["B"] - 1.2)) - 1.0,
        }
    return eff


# ------------------------------------------------------------------------------- aggregates
def compute_aggregates() -> dict[str, Any]:
    """All aggregates (ratios/differences; no NIST value) as the former ``p1c-aggregates.json``."""
    tbl = theta_table(gottschalk_tex())
    makarova = parse_makarova_table1(makarova_pdf_text())
    src = BetheStoppingSource()
    pb = {n: pb_attribution(n, m, tbl) for n, m in MATERIALS.items()}
    ps = {n: pstar_element_check(n, m, tbl) for n, m in MATERIALS.items()}
    cc = compound_check(makarova, gottschalk_tex(), src)
    agg: dict[str, Any] = {
        "units": "relative (fraction); multiply by 100 for %",
        "pb_attribution": {},
        "controls_Be_Al_Cu": {},
    }
    for m, o in pb.items():
        base = o["baseline"]
        num = o["numerics"]
        e = o["electrons"]
        row = {
            "baseline_dev": base,
            "baseline_max_abs": _mx(base),
            "R1_ours_over_G_minus_1": o["rhoR1_ours_over_G_minus_1"],
            "R1_ours_over_PSTAR_minus_1": ps[m]["R1_ours_over_pstar_minus_1"],
            "R1_PSTAR_over_G_minus_1": ps[m]["R1_pstar_over_G_minus_1"],
            "resid_range_ours_at_EG_over_G_minus_1": o["residual_range_ours_at_EG_over_G_minus_1"],
            "dlntheta_dlnx_fixed_path": o["dlntheta_dlnx_fixed_path"],
            "G_absolute_thickness_on_our_path_dev": o["abs_thickness_on_our_path"],
            "icru49_path_dev": {
                k: v for k, v in o["icru49_path"].items() if not k.endswith("_max")
            },
            "icru49_path_max_abs_over_variants": max(
                _mx(v) for k, v in o["icru49_path"].items() if not k.endswith("_max")
            ),
            "pstar_path_dev": ps[m]["theta_on_pstar_path_over_thetaH_minus_1"],
            "pstar_path_max_abs": _mx(ps[m]["theta_on_pstar_path_over_thetaH_minus_1"]),
            "fano_variants_minus_baseline": {
                k: {f: v[f] - base[f] for f in base}
                for k, v in e.items()
                if k not in ("delta_b_at_EG", "slab_b_B")
            },
            "slab_delta_b_range": [
                min(x["delta_b"] for x in e["slab_b_B"].values()),
                max(x["delta_b"] for x in e["slab_b_B"].values()),
            ],
            "numerics_max_abs_change_vs_baseline": max(
                abs(num[k][f] - base[f]) for k in num for f in base
            ),
        }
        (agg["pb_attribution"] if m == "Pb" else agg["controls_Be_Al_Cu"])[m] = row
    alt = alt_compound_rules(makarova)
    agg["compound"] = {
        "makarova_LOOKUP": {
            k: [
                {
                    "t_g_cm2": e["t_g_cm2"],
                    "x_over_R1_ours": e.get("x_over_R1_ours"),
                    "ours_over_LOOKUP_minus_1": e["ours_over_theta_hanson_LOOKUP_minus_1"],
                    "ours_over_Go93_measured_minus_1": e["ours_over_measured_Go93_minus_1"],
                    "LOOKUP_over_Go93_measured_minus_1": e["LOOKUP_hanson_over_measured_minus_1"],
                    "bethe_ZZ1_over_LOOKUP_minus_1": e["bethe_ZZ1_over_theta_hanson_minus_1"],
                    "range_sensitivity_I_pm10pct": e.get("range_sensitivity"),
                }
                for e in v
            ]
            for k, v in cc["makarova_tab_results"].items()
        },
        "makarova_thick_on_pstar_path": thick_on_pstar_path(makarova),
        "gottschalk2010_tblB": cc["gottschalk2010_tblB"],
        "gottschalk2010_modStack": cc["gottschalk2010_modStack"],
        "alt_compound_rules_minus_committed": {
            k: {r: v[r] - v["committed"] for r in v if r != "committed"} for k, v in alt.items()
        },
        "water_b_shift_minus_0.040_theta_effect": water_b_shift_effect(),
    }
    return agg


def print_summary(agg: dict[str, Any]) -> None:
    """Aggregate tables (percent / ratios only)."""
    pc = lambda x: f"{100 * x:+7.3f}%"  # noqa: E731
    pb = agg["pb_attribution"]["Pb"]
    print("Pb theta/theta_Hanson - 1 by x/R1 (percent)")
    print("  path                         " + " ".join(f"{f:>8}" for f in pb["baseline_dev"]))
    for label, d in (
        ("ionmc Bethe ranges (baseline)", pb["baseline_dev"]),
        ("Gottschalk-node spline path  ", pb["icru49_path_dev"]["spline_lnR_lnT_nodes_ge_0.1"]),
        ("kappa(E)-scaled ionmc table  ", pb["icru49_path_dev"]["kappa_nodes_ge_0.1"]),
        ("PSTAR spline path            ", pb["pstar_path_dev"]),
    ):
        print(f"  {label} " + " ".join(f"{100 * v:+8.3f}" for v in d.values()))
    print(
        f"Pb: max|dev| nodes-spline {pc(max(_mx(v) for k, v in pb['icru49_path_dev'].items() if k.startswith('spline')))},"
        f" PSTAR {pc(pb['pstar_path_max_abs'])}, baseline max {pc(pb['baseline_max_abs'])}"
    )
    print(
        f"Pb: R1(ours)/R1(G)-1 {pc(pb['R1_ours_over_G_minus_1'])}, R1(PSTAR)/R1(G)-1 "
        f"{pc(pb['R1_PSTAR_over_G_minus_1'])}, numerics max change "
        f"{pb['numerics_max_abs_change_vs_baseline']:.2e}"
    )
    for m, r in agg["controls_Be_Al_Cu"].items():
        print(
            f"{m}: R1(ours)/R1(G)-1 {pc(r['R1_ours_over_G_minus_1'])}, max|dev| ICRU49 paths "
            f"{pc(r['icru49_path_max_abs_over_variants'])}, PSTAR {pc(r['pstar_path_max_abs'])}"
        )
    print("Compounds, ours/LOOKUP - 1 (percent) [Bethe Z(Z+1) variant]")
    for k, v in agg["compound"]["makarova_LOOKUP"].items():
        print(
            f"  {k:15s} "
            + "; ".join(
                f"{e['t_g_cm2']:g} g/cm2: {pc(e['ours_over_LOOKUP_minus_1'])} "
                f"[{pc(e['bethe_ZZ1_over_LOOKUP_minus_1'])}]"
                for e in v
            )
        )
    print("alternative compound Fano rules minus committed:", agg["compound"]["alt_compound_rules_minus_committed"])
    print("modStack:", {k: pc(v["ours_over_thetaH_minus_1"]) for k, v in agg["compound"]["gottschalk2010_modStack"].items()})
    print(
        "tbl:B b_ours - b_G:",
        {k: round(v["b_ours_minus_b_G"], 4) for k, v in agg["compound"]["gottschalk2010_tblB"].items()},
    )



def _json_safe(obj: Any) -> Any:
    """Recursively convert numpy scalars/arrays to Python values and NaN/inf to ``None`` (standard JSON)."""
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_json_safe(v) for v in obj.tolist()]
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return None if (f != f or f in (float("inf"), float("-inf"))) else f
    if isinstance(obj, (np.integer, int)) and not isinstance(obj, bool):
        return int(obj)
    return obj

def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--out", type=Path, default=REPO / "validation/plans/v3-008/p1c-aggregates.json",
        help="aggregates JSON (ratios only; no PSTAR values)",
    )  # fmt: skip
    args = ap.parse_args(argv)
    agg = compute_aggregates()
    print_summary(agg)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    agg["_null_means"] = (
        "undefined value (NaN or infinity in the computation, e.g. bethe_b_minus_b_G where the "
        "configuration has no Bethe-path b); serialised as null so the file is standard JSON"
    )
    args.out.write_text(json.dumps(_json_safe(agg), indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
