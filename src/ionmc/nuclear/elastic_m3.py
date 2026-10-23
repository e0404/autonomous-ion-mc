"""M3 = (ii') shape for the proton elastic table (V3-005D plan section 4, Amendment 20 (b)3).

The p + A angular shape is the angular distribution of the elastic scattering of a NEUTRON of the
same laboratory kinetic energy on the same nucleus (isospin mirror; the Coulomb-free nuclear
amplitude of p + A is approximated by that of n + A, no Coulomb-nuclear interference, no
isospin-breaking). It is taken from the ENDF/B-VIII.0 neutron evaluation, MF4/MT2 (elastic angular
distributions). Only the targets named in the plan row (C-12, O-16) get this shape; every other
target keeps the pinned black-disk shape (recorded). ``sigma_el`` stays Barashenkov/BGG.

ENDF-6 MF4 layout read here (ENDF-6 Formats Manual ch. 4): HEAD ``[ZA, AWR, 0, LTT, 0, 0]``,
CONT ``[0, AWR, LI, LCT, 0, 0]`` (LI = 1 isotropic; LCT = 1 laboratory, 2 centre-of-mass frame),
then for LTT = 1 a TAB2 (NE energies) and NE LIST records ``[T, E, LT, 0, NL, 0 / a_1 .. a_NL]``
(``f(mu) = sum_{l=0}^{NL} (2l+1)/2 a_l P_l(mu)``, ``a_0 = 1``); LTT = 2 a TAB2 and NE TAB1 records
``[T, E, LT, 0, NR, NP / mu, f(mu)]``; LTT = 3 the Legendre block (low energies) followed by the
tabulated block (high energies). Energy interpolation: INT of the TAB2 (1 histogram or 2 lin-lin;
others fail closed), applied to ``f`` at fixed mu (for Legendre data to the coefficients, which is
identical for lin-lin). Temperature-dependent records (LT != 0) fail closed. A laboratory-frame
distribution is converted exactly to the CM frame for elastic scattering with the mass ratio
``A = AWR``: ``mu_lab = (1 + A mu) / sqrt(A^2 + 2 A mu + 1)``, ``f_cm = f_lab dmu_lab/dmu``,
``dmu_lab/dmu = A^2 (A + mu) / (A^2 + 2 A mu + 1)^{3/2}``.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ionmc.data import endf6
from ionmc.nuclear import law5
from ionmc.nuclear.elastic_glauber import DATA_DIR, quantiles_mu

F64 = NDArray[np.float64]
M3_INPUTS_JSON = DATA_DIR / "elastic_m3_inputs.json"
ABOVE_EMAX_POLICIES = ("error", "hold-last", "black-disk")
COPY_CHECK_E_MIN_EV = 20.0e6
COPY_CHECK_N_MU = 181
COPY_CHECK_REL = 1.0e-3
GATE_N_MU = 40001


class M3Error(RuntimeError):
    """Fail-closed error of the M3 construction."""


@dataclass
class Mf4Angular:
    """Parsed MF4/MT2 angular distribution of one material."""

    awr: float
    ltt: int
    li: int
    lct: int
    leg_e: F64  # eV
    leg_a: list[F64]  # a_1..a_NL per energy
    leg_int: int
    tab_e: F64
    tab: list[endf6.Tab1]
    tab_int: int

    @property
    def emax_ev(self) -> float:
        if self.li == 1:  # isotropic at all energies (no records)
            return math.inf
        parts = [self.leg_e, self.tab_e]
        return float(max(p[-1] for p in parts if p.size))

    @property
    def emin_ev(self) -> float:
        if self.li == 1:
            return 0.0
        parts = [self.leg_e, self.tab_e]
        return float(min(p[0] for p in parts if p.size))

    def record_energies_ev(self) -> F64:
        """Sorted unique incident energies [eV] of all records (Legendre and tabulated)."""
        return np.unique(np.concatenate([self.leg_e, self.tab_e]))

    def _legendre_f(self, mu: F64, a: F64) -> F64:
        p = law5.legendre(mu, a.size)
        out = 0.5 * p[0]
        for lval in range(1, a.size + 1):
            out = out + 0.5 * (2 * lval + 1) * a[lval - 1] * p[lval]
        return np.asarray(out)

    def f_native(self, mu: F64, e_ev: float) -> F64:
        """``f(mu)`` (normalised ``int f dmu = 1``) in the evaluation's own frame at ``e_ev``."""
        mu = np.asarray(mu, dtype=np.float64)
        if self.li == 1:
            return np.full_like(mu, 0.5)
        if not self.emin_ev - 1e-9 <= e_ev <= self.emax_ev * (1 + 1e-12):
            raise M3Error(f"E = {e_ev} eV outside the MF4 range [{self.emin_ev}, {self.emax_ev}]")
        use_tab = self.tab_e.size > 0 and (self.leg_e.size == 0 or e_ev >= self.tab_e[0])
        if use_tab:
            e, interp = self.tab_e, self.tab_int
            j = int(np.clip(np.searchsorted(e, e_ev, side="right") - 1, 0, e.size - 2))
            if e.size == 1:
                return np.asarray(self.tab[0].interpolate(mu))
            f0 = self.tab[j].interpolate(mu)
            f1 = self.tab[j + 1].interpolate(mu)
        else:
            e, interp = self.leg_e, self.leg_int
            j = int(np.clip(np.searchsorted(e, e_ev, side="right") - 1, 0, e.size - 2))
            n = max(self.leg_a[j].size, self.leg_a[j + 1].size)
            pad = [np.pad(self.leg_a[i], (0, n - self.leg_a[i].size)) for i in (j, j + 1)]
            f0, f1 = self._legendre_f(mu, pad[0]), self._legendre_f(mu, pad[1])
        frac = 0.0 if interp == 1 else (e_ev - e[j]) / (e[j + 1] - e[j])
        return np.asarray((1.0 - frac) * f0 + frac * f1)

    def f_cm(self, mu_cm: F64, e_ev: float) -> F64:
        """``f`` per unit CM cosine (CM data are returned as stored; LCT = 1 is converted)."""
        mu_cm = np.asarray(mu_cm, dtype=np.float64)
        if self.lct == 2 or self.li == 1:
            return self.f_native(mu_cm, e_ev)
        a = self.awr
        s = np.sqrt(a * a + 2.0 * a * mu_cm + 1.0)
        mu_lab = np.clip((1.0 + a * mu_cm) / s, -1.0, 1.0)
        return np.asarray(self.f_native(mu_lab, e_ev) * a * a * (a + mu_cm) / s**3)


def _check_int(nbt: NDArray[np.int64], law: NDArray[np.int64]) -> int:
    if nbt.size != 1 or int(law[0]) not in (1, 2):
        raise endf6.UnsupportedEndfError(f"MF4 energy interpolation {law.tolist()} not supported")
    return int(law[0])


def parse_mf4(section: endf6.EndfSection) -> Mf4Angular:
    """MF4 section (module docstring). Raises :class:`endf6.UnsupportedEndfError` (fail closed)
    for anything outside the supported subset."""
    if section.mf != 4:
        raise endf6.EndfError("not an MF4 section")
    rd = endf6._Reader(section.lines)
    return _parse_mf4(rd, rd.cont())


def _parse_mf4(rd: endf6._Reader, head: endf6.Cont) -> Mf4Angular:
    ltt = head.l2  # HEAD: ZA, AWR, 0, LTT, 0, 0
    awr = head.c2
    c2 = rd.cont()  # 0, AWR, LI, LCT, 0, 0
    li, lct = c2.l1, c2.l2
    if lct not in (1, 2):
        raise endf6.UnsupportedEndfError(f"LCT = {lct}")
    empty: F64 = np.zeros(0)
    if li == 1:
        return Mf4Angular(awr, ltt, li, lct, empty, [], 2, empty, [], 2)
    if li != 0 or ltt not in (1, 2, 3):
        raise endf6.UnsupportedEndfError(f"LI = {li}, LTT = {ltt}")
    leg_e: list[float] = []
    leg_a: list[F64] = []
    leg_int = 2
    tab_e: list[float] = []
    tab: list[endf6.Tab1] = []
    tab_int = 2

    def read_leg() -> None:
        nonlocal leg_int
        h, nbt, law = rd.tab2()
        leg_int = _check_int(nbt, law)
        for _ in range(h.n2):
            lh = rd.cont()
            if lh.l1 != 0:
                raise endf6.UnsupportedEndfError("temperature-dependent MF4 record (LT != 0)")
            leg_e.append(lh.c2)
            leg_a.append(rd.list_body(lh))

    def read_tab() -> None:
        nonlocal tab_int
        h, nbt, law = rd.tab2()
        tab_int = _check_int(nbt, law)
        for _ in range(h.n2):
            th = rd.cont()
            if th.l1 != 0:
                raise endf6.UnsupportedEndfError("temperature-dependent MF4 record (LT != 0)")
            tab_e.append(th.c2)
            tab.append(rd.tab1_body(th))

    if ltt in (1, 3):
        read_leg()
    if ltt in (2, 3):
        read_tab()
    return Mf4Angular(awr, ltt, li, lct, np.asarray(leg_e), leg_a, leg_int,
                      np.asarray(tab_e), tab, tab_int)  # fmt: skip


def mf1_header(mat: endf6.EndfMaterial) -> dict[str, Any]:
    """The MF1/MT451 header fields used by the copy check."""
    r = endf6._Reader(mat.sections[(1, 451)].lines)
    h1 = r.cont()
    h2 = r.cont()
    h3 = r.cont()
    h4 = r.cont()
    text = [r.text() for _ in range(h4.n1)]
    return {"za": h1.c1, "awr": h1.c2, "nlib": h1.n1, "nmod": h1.n2,
            "elis": h2.c1, "nfor": h2.n2, "awi": h3.c1, "emax_ev": h3.c2, "lrel": h3.l1,
            "nsub": h3.n1, "nver": h3.n2, "text": text, "mat": mat.mat}  # fmt: skip


def mf1_copy_check(name: str, mat: endf6.EndfMaterial, expected: dict[str, Any]) -> dict[str, Any]:
    """MF1/MT451 header check of an acquired file against the ENDF/B-VIII.0 neutron catalogue
    entries stored in ``elastic_m3_inputs.json`` (``expected``): ZA, MAT, library (NLIB = 0,
    ENDF/B), version (NVER = 8), release (LREL), sublibrary (NSUB = 10, incident neutron) and
    AWI = 1 (neutron projectile). The description text is recorded verbatim, not matched.
    Raises :class:`M3Error` on any mismatch (the file is not the catalogue material)."""
    h = mf1_header(mat)
    checks = {
        "za": h["za"] == expected["za"], "mat": h["mat"] == expected["mat"],
        "nlib": h["nlib"] == expected["nlib"], "nver": h["nver"] == expected["nver"],
        "lrel": h["lrel"] == expected["lrel"], "nsub": h["nsub"] == expected["nsub"],
        "awi": abs(h["awi"] - 1.0) < 1e-9,
    }  # fmt: skip
    if not all(checks.values()):
        raise M3Error(f"{name}: MF1 does not match the catalogue entry: {checks} {h}")
    return {"checks": checks, "header": h, "expected": expected}


@dataclass(frozen=True)
class M3Config:
    """Options of the M3 build (kept out of ``ElasticBuildOptions``)."""

    input_dir: Path
    above_emax: str = "error"
    n_s: int = 8193
    inputs_json: Path = M3_INPUTS_JSON
    allow_unpinned_inputs: bool = False

    def __post_init__(self) -> None:
        if self.above_emax not in ABOVE_EMAX_POLICIES:
            raise M3Error(f"above_emax must be one of {ABOVE_EMAX_POLICIES}")


class M3Model:
    """n + A MF4/MT2 shapes for the targets of the plan row (read from ``input_dir``: every
    parseable ENDF file whose ZA matches an entry of the inputs JSON is used)."""

    def __init__(self, cfg: M3Config) -> None:
        self.cfg = cfg
        self.inputs = json.loads(Path(cfg.inputs_json).read_text(encoding="utf-8"))
        self.mf4: dict[str, Mf4Angular] = {}
        self.files: dict[str, dict[str, Any]] = {}
        by_za = {int(v["za"]): n for n, v in self.inputs["targets"].items()}
        pinned = {v for v in self.inputs.get("file_sha256", {}).values() if v}
        self.unpinned = not pinned
        if self.unpinned and not cfg.allow_unpinned_inputs:
            raise M3Error("file_sha256 in elastic_m3_inputs.json is empty (unpinned ENDF files)")
        for p in sorted(Path(cfg.input_dir).iterdir()):
            if not p.is_file():
                continue
            raw = p.read_bytes()
            try:
                mat = endf6.parse_endf(raw.decode("ascii", "replace"))
            except endf6.EndfError:
                continue
            name = by_za.get(int(round(mat.za)))
            if name is None:
                continue
            sha = hashlib.sha256(raw).hexdigest()
            if pinned and sha not in pinned:
                raise M3Error(f"{p.name}: sha256 {sha} is not pinned")
            if (4, 2) not in mat.sections:
                raise M3Error(f"{p.name}: no MF4/MT2")
            exp = {**self.inputs["catalogue_common"], **self.inputs["targets"][name]}
            chk = mf1_copy_check(name, mat, exp)
            self.mf4[name] = parse_mf4(mat.sections[(4, 2)])
            self.files[name] = {"file": p.name, "sha256": sha, "mat": mat.mat, "mf1_check": chk}
        missing = set(by_za.values()) - set(self.mf4)
        if missing:
            raise M3Error(f"no ENDF file found for {sorted(missing)}")
        self.copy_check = self._cross_copy_check()

    def _cross_copy_check(self) -> dict[str, Any]:
        """Pre-registered copy check (acceptance (b'')4 item 2, Amendment 15 (d)1 transposed):
        the O-16 and C-12 MF4/MT2 angular distributions (Legendre or tabulated, converted to CM)
        are compared at EVERY shared incident energy above 20 MeV on 181 equally spaced mu_CM
        points. NOT a copy iff the relative difference ``|a - b| / max(|a|, |b|)`` exceeds 1e-3
        at > 50 % of the points at > 50 % of the shared energies; otherwise ``M3Error``. Fail
        closed: if O-16 and C-12 are not both present, or share no incident energy above 20 MeV,
        the criterion cannot be established and ``M3Error`` is raised (no table)."""
        o, c = self.mf4.get("O-16"), self.mf4.get("C-12")
        if o is None or c is None:
            raise M3Error("copy check cannot be established: O-16 and C-12 are not both present")
        shared = [
            float(e) for e in o.record_energies_ev()
            if e > COPY_CHECK_E_MIN_EV and np.any(np.isclose(c.record_energies_ev(), e, rtol=1e-12,
                                                             atol=0.0))
        ]  # fmt: skip
        if not shared:
            raise M3Error(
                "copy check cannot be established: O-16 and C-12 share no record energy above "
                "20 MeV"
            )
        mu = np.linspace(-1.0, 1.0, COPY_CHECK_N_MU)
        n_diff = 0
        for e in shared:
            a, b = o.f_cm(mu, e), c.f_cm(mu, e)
            rel = np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)
            rel = np.where((a == 0.0) & (b == 0.0), 0.0, rel)
            n_diff += int(np.mean(rel > COPY_CHECK_REL) > 0.5)
        frac = n_diff / len(shared)
        if frac <= 0.5:
            raise M3Error(
                f"O-16 MF4 distributions are a C-12 copy above 20 MeV: they differ by > 1e-3 at "
                f"> 50 % of the points at only {n_diff}/{len(shared)} shared energies"
            )
        return {"evaluated": True, "n_mu_points": COPY_CHECK_N_MU, "shared_energies": len(shared),
                "energies_that_differ": n_diff, "fraction": frac}  # fmt: skip

    def has(self, name: str) -> bool:
        return name in self.mf4

    def emax_mev(self, name: str) -> float:
        return self.mf4[name].emax_ev * 1e-6

    def _gates(self, name: str, e_ev: float) -> dict[str, Any]:
        """Verification gates on the CM distribution at ``e_ev`` (fail closed, no clipping):
        positivity (no value below -1e-12 of the maximum), normalisation ``|int f dmu - 1| <=
        1e-6`` (Simpson on 40001 mu points), monotone CDF and a finite first moment <1 - mu>."""
        mu = np.linspace(-1.0, 1.0, GATE_N_MU)
        f = self.mf4[name].f_cm(mu, e_ev)
        if not np.all(np.isfinite(f)):
            raise M3Error(f"{name} at {e_ev * 1e-6:.6g} MeV: non-finite angular distribution")
        fmax = float(np.max(np.abs(f)))
        if float(f.min()) < -1e-12 * fmax:
            raise M3Error(
                f"{name} at {e_ev * 1e-6:.6g} MeV: negative angular density {f.min():.3g}"
            )
        h = mu[1] - mu[0]
        simpson = h / 3.0 * (f[0] + f[-1] + 4.0 * f[1:-1:2].sum() + 2.0 * f[2:-1:2].sum())
        if abs(simpson - 1.0) > 1e-6:
            raise M3Error(f"{name} at {e_ev * 1e-6:.6g} MeV: normalisation {simpson!r} (> 1e-6)")
        cdf = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * h)))
        if np.any(np.diff(cdf) < -1e-15):
            raise M3Error(f"{name} at {e_ev * 1e-6:.6g} MeV: CDF is not monotone")
        mom = float(np.sum(0.5 * ((1.0 - mu) * f)[1:] * h + 0.5 * ((1.0 - mu) * f)[:-1] * h))
        if not math.isfinite(mom) or not 0.0 <= mom <= 2.0:
            raise M3Error(f"{name} at {e_ev * 1e-6:.6g} MeV: first moment {mom!r}")
        return {"norm_int_f_dmu": float(simpson), "min_f": float(f.min()),
                "first_moment_1_minus_mu": mom, "cdf_monotone": True}  # fmt: skip

    def check_midpoint(self, name: str, e_mev: float) -> None:
        """Positivity (and finiteness) at an energy midpoint of the pinned grid; above E_top the
        distribution is that of E_top. ``M3Error`` on a negative value (no clipping)."""
        m = self.mf4[name]
        self._gates(name, min(e_mev * 1e6, m.emax_ev))

    def node(
        self, name: str, e_mev: float, n_q: int, pcm_mev: float | None = None,
        pcm_top_mev: float | None = None,
    ) -> tuple[F64, dict[str, Any]]:  # fmt: skip
        """256-bin edges of the p + A shape at ``e_mev`` after all gates. At or below E_top the
        n + A distribution at ``e_mev`` itself. Above E_top (policy ``hold-last``; the caller
        must pass the p + A CM momenta ``pcm_mev`` at the node and ``pcm_top_mev`` at E_top): the
        E_top distribution is HELD FIXED IN MOMENTUM TRANSFER q = 2 p sin(theta/2) and rescaled to
        the node's p_CM: its equiprobable edges in ``s' = sin(theta_top/2)`` map to
        ``s = s' p_top / p`` (the distribution of q is unchanged, the CDF in q is identical). The
        node's q range exceeds the E_top kinematic limit 2 p_top, so the density is zero for
        ``s > p_top / p`` (q beyond the E_top limit); the normalisation is unchanged (1)."""
        m = self.mf4[name]
        above = e_mev * 1e6 > m.emax_ev * (1.0 + 1e-12)
        e_use = m.emax_ev if above else e_mev * 1e6
        gates = self._gates(name, e_use)
        s = np.linspace(0.0, 1.0, self.cfg.n_s)
        d = m.f_cm(1.0 - 2.0 * s * s, e_use)
        edges, dev = quantiles_mu(s, d, n_q)
        rec: dict[str, Any] = {"e_mev": e_mev, "e_used_mev": e_use * 1e-6, "q_scaled": False,
                               "norm_deviation": dev, **gates}  # fmt: skip
        if above:
            if pcm_mev is None or pcm_top_mev is None or pcm_mev < pcm_top_mev:
                raise M3Error("hold-last needs p_CM at the node and at E_top (node >= E_top)")
            ratio = pcm_top_mev / pcm_mev
            s_edges = np.sqrt(np.maximum((1.0 - edges) / 2.0, 0.0)) * ratio
            edges = np.maximum.accumulate(1.0 - 2.0 * s_edges * s_edges)
            edges[-1] = 1.0  # edges[0] = 1 - 2 (p_top/p)^2 > -1: zero density below it
            rec.update(q_scaled=True, p_ratio_top_over_node=ratio, s_max=ratio)
        return edges, rec

    def provenance(self) -> dict[str, Any]:
        digest = hashlib.sha256(Path(self.cfg.inputs_json).read_bytes()).hexdigest()
        return {"files": self.files, "cross_copy_check": self.copy_check,
                "inputs_json_sha256": digest, "inputs_unpinned": self.unpinned,
                "above_emax": self.cfg.above_emax,
                "emax_mev": {n: self.emax_mev(n) for n in self.mf4},
                "frames": {n: ("CM" if m.lct == 2 else "lab (converted)")
                           for n, m in self.mf4.items()},
                "ltt": {n: m.ltt for n, m in self.mf4.items()}}  # fmt: skip
