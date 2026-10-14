"""D-E0 (exploratory, report-only; Amendment 19): deterministic kernels of the elastic channel.

Definitions are those of plan section 2.1 (``validation/plans/v3-005d-elastic-attribution.md``), evaluated on
the pinned elastic table with the project (Bethe, as V5/V11) CSDA range of water and the exact two-body kinematics:

- depth deficit ``Delta(E, mu) = R(E) - R(T1) cos(theta_lab)`` of the primary (for p-p the faster proton);
- ``K_w = sum_t int_0^{0.98 R0} n_t sigma_t(E(z)) P(Delta > w) dz``, w = 1, 2, 3, 5 mm (no straggling);
- ``K*`` with the archived EM-only ionmc IDD as kernel: ``int n sigma <1 - IDD_EM(z_pk -+ Delta)/IDD_EM(z_pk)> dz``
  (the plan formula writes ``z_pk - Delta`` ("as frozen", ``K_star_literal``; a SIGN ERROR, see the dated note in the
  report: a range shortened by Delta puts the proton's dose at ``z_pk`` where the unshifted curve is at ``z_pk + Delta``,
  ``K_star_shifted``, the physically meaningful kernel; the plan's prose describes it). Both are computed; no conclusion
  is drawn from the as-frozen value);
- ``G`` plateau gain (operational form, plan wording under-specified): over the window 20-60 mm,
  ``G = G_path + G_recoil + G_pp_slow``; ``G_path`` = sum over scatterings upstream of the depth z' of
  ``n sigma <(1/cos theta) 1[R(T1) cos theta >= z'-z] - 1>`` weighted by the stopping power S(z') and normalised to
  ``int S dz'`` over the window; ``G_recoil`` = local recoil energy of p+A, ``G_pp_slow`` = the slow proton of
  p-p deposited locally (upper estimate), both ``int n sigma <T> dz / int S dz'``;
- ``theta_w(E)``: ``R(E)(1 - cos theta_w) = w`` (w = 2, 3 mm); tail moment ``T = <(1 - cos theta_lab) 1[theta_lab > 3 theta_med]>``.

Run: ``python -m validation.scripts.elastic_attribution.de0_kernels``.
"""

from __future__ import annotations

import difflib
import json
import math
from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import NDArray

from . import common as c
from .de1_shapes import FAMILIES, Shapes, table_bin_density_first_minimum_deg

ENERGIES = (50, 100, 150, 200)
W_MM = (1.0, 2.0, 3.0, 5.0)
N_U = 2048
U = (np.arange(N_U) + 0.5) / N_U
N_Z = 300
R_FRACTION = 0.98
WINDOW_MM = (20.0, 60.0)
SIGMA_R_MM = {150: 1.9, 200: 3.1}  # plan 2.1: Gaussian K* form (report-only)
PLANNING_R_O_FM = {100: 2.58, 150: 2.60, 200: 2.64}  # decision 0041 D-1 planning values
PLANNING_THETA_RMS_O = {100: 21.0, 150: 17.0, 200: 14.0}
EXPECTATION = {150: 0.985, 200: 0.986}  # Amendment 14 (f)
PLANNING_GAIN = {150: 0.012, 200: 0.014}  # Amendment 14 (c) V11 row
REVIEW_CRUDE = {  # reviewer's removal.py table (review.md B.1), 150 MeV, p+O only: moment, removal 2 mm, 3 mm
    "sharp": (0.0305, 0.119, 0.091), "s0.5": (0.0182, 0.107, 0.077), "s0.9": (0.0123, 0.091, 0.061),
    "Rx1.2": (0.0250, 0.096, 0.069), "gauss": (0.0196, 0.125, 0.096)}  # fmt: skip


def rng_mm_factory(water: c.Water) -> Callable[[Any], Any]:
    return water.range_mm


def primary_event(
    sh: Shapes, name: str, kind: str, e_mev: float, u: NDArray[np.float64]
) -> dict[str, NDArray[np.float64]]:
    """Primary proton after one elastic event at ``e_mev``: lab kinetic energy, lab cosine, the energy given to
    the other body (p+A recoil, p-p slower proton) and ``mu_CM``. p-p: the faster proton continues."""
    from ionmc.nuclear.bgg import M_PROTON_MEV
    from ionmc.nuclear.elastic_kin import two_body

    it = sh.index(name)
    m_t = float(sh.a["target_mass_mev"][it])
    mu = sh.mu_quantiles("table" if name != "O-16" else kind, name, e_mev, u)
    tb = two_body(e_mev, M_PROTON_MEV, m_t, mu)
    if name != "H-1":
        return {"mu": mu, "t1": tb.t1, "cos": tb.cos_lab1, "t_other": tb.t2}
    cos2 = tb.pz2 / np.hypot(tb.px2, tb.pz2)
    first = tb.t1 >= tb.t2
    return {"mu": mu, "t1": np.where(first, tb.t1, tb.t2), "cos": np.where(first, tb.cos_lab1, cos2),
            "t_other": np.where(first, tb.t2, tb.t1)}  # fmt: skip


def deficit_mm(
    e_mev: float, ev: dict[str, NDArray[np.float64]], rng: Callable[[Any], Any]
) -> NDArray[np.float64]:
    """``Delta = R(E) - R(T1) cos(theta_lab)`` [mm] (a backward proton has ``Delta > R``)."""
    return np.asarray(float(rng(e_mev)) - np.asarray(rng(ev["t1"])) * ev["cos"])


def theta_w_deg(e_mev: float, w_mm: float, rng: Callable[[Any], Any]) -> float:
    return math.degrees(math.acos(1.0 - w_mm / float(rng(e_mev))))


def pointwise(sh: Shapes, water: c.Water, name: str, kind: str, e_mev: float) -> dict[str, Any]:
    """Shape moments of one interaction at ``e_mev``: CM and lab quantities, threshold-angle exceedance, tail moment."""
    ev = primary_event(sh, name, kind, e_mev, U)
    mu, cos = ev["mu"], ev["cos"]
    th_lab = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))
    th_med = float(np.median(th_lab))
    d = deficit_mm(e_mev, ev, water.range_mm)
    out: dict[str, Any] = {
        "mean_1_minus_mu_cm": float(np.mean(1.0 - mu)), "theta_rms_cm_deg": float(np.degrees(np.sqrt(np.mean(np.arccos(mu) ** 2)))),
        "theta_med_lab_deg": th_med, "mean_1_minus_cos_lab": float(np.mean(1.0 - cos)),
        "T_tail_moment": float(np.mean((1.0 - cos) * (th_lab > 3.0 * th_med))), "mean_t_other_mev": float(np.mean(ev["t_other"])),
        "range_mm": float(water.range_mm(e_mev)),
    }  # fmt: skip
    for w in (2.0, 3.0):
        tw = theta_w_deg(e_mev, w, water.range_mm)
        out[f"theta_{w:g}mm_deg"] = tw
        out[f"P_theta_lab_gt_theta_{w:g}mm"] = float(np.mean(th_lab > tw))
    for w in W_MM:
        out[f"P_deficit_gt_{w:g}mm"] = float(np.mean(d > w))
    return out


def sigma_nsigma(sh: Shapes, water: c.Water | None, name: str, e_mev: float) -> float:
    """``n sigma`` [1/cm] of table target ``name`` in water (``n`` = 1 per cm3 if the target is not in water)."""
    sig = float(sh.tab.sigma_barn(sh.index(name), e_mev)) * 1.0e-24
    if water is not None and name in water.targets:
        t = water.targets[name]
        return t["n_per_cm3"] * t["scale"] * sig
    return sig


def idd_curve(energy: int) -> tuple[Callable[[Any], Any], float]:
    z, idd, met = c.idd_mean(energy, "emonly")
    zpk = float(met["peak_depth_mm"])
    fpk = float(np.interp(zpk, z, idd))
    return (lambda x: np.interp(np.asarray(x, dtype=float), z, idd, left=0.0, right=0.0) / fpk), zpk


def kstar_node(f: Callable[[Any], Any], zpk: float, d: NDArray[np.float64]) -> tuple[float, float]:
    """Node average of ``1 - IDD(z_pk -+ Delta)/IDD(z_pk)``: ``(as frozen, shifted)``. A proton whose range is shortened by
    Delta has the depth-dose curve ``IDD(z + Delta)``, so its contribution at ``z_pk`` is ``IDD(z_pk + Delta)`` (shifted,
    the physically meaningful form). The plan formula ``IDD(z_pk - Delta)`` (as frozen) samples the proximal rise: a sign
    error of the plan formula (its prose describes the shifted form); see the dated note in the D-E0 report."""
    return float(np.mean(1.0 - f(zpk - d))), float(np.mean(1.0 - f(zpk + d)))


def kernels(
    sh: Shapes,
    water: c.Water,
    e0: int,
    kind: str,
    targets: tuple[str, ...],
    with_kstar: bool,
    n_z: int = N_Z,
) -> dict[str, Any]:
    """``P_el``, ``K_w`` (all ``W_MM``) and, if the IDD exists, ``K*`` (literal, shifted, Gaussian) per target and total.
    ``targets`` not in water (``C-12``) get ``n = 1`` per cm3 along the water CSDA path (hypothetical, per unit n)."""
    r0 = float(water.range_mm(e0))
    zc = (np.arange(n_z) + 0.5) / n_z * R_FRACTION * r0
    dz_cm = R_FRACTION * r0 / n_z / 10.0
    en = water.energy_at_depth(float(e0), zc)
    idd = idd_curve(e0) if with_kstar else None
    res: dict[str, Any] = {}
    for name in targets:
        p_el, kw = 0.0, dict.fromkeys(W_MM, 0.0)
        ks = {"literal": 0.0, "shifted": 0.0, "gauss": 0.0}
        for z_i, e in zip(zc, en, strict=True):
            e = float(e)
            lo = float(sh.tab.arrays["target_e_min_mev"][sh.index(name)])
            if e < max(lo, 1.0):
                continue
            ns = sigma_nsigma(sh, water, name, e) * dz_cm
            if ns <= 0.0:
                continue
            ev = primary_event(sh, name, kind, e, U)
            d = deficit_mm(e, ev, water.range_mm)
            p_el += ns
            for w in W_MM:
                kw[w] += ns * float(np.mean(d > w))
            if idd is not None:
                f, zpk = idd
                k_lit, k_sh = kstar_node(f, zpk, d)
                ks["literal"] += ns * k_lit
                ks["shifted"] += ns * k_sh
                ks["gauss"] += ns * float(
                    np.mean(1.0 - np.exp(-(d**2) / (2.0 * SIGMA_R_MM[e0] ** 2)))
                )
        res[name] = {"P_el": p_el, **{f"K_{w:g}mm": kw[w] for w in W_MM}}
        if idd is not None:
            res[name].update({f"K_star_{k}": v for k, v in ks.items()})
    return res


def plateau_gain(
    sh: Shapes, water: c.Water, e0: int, kind: str, dz_mm: float = 1.0
) -> dict[str, float]:
    """``G`` of the module docstring over the window 20-60 mm (O and H-1 of water)."""
    z1, z2 = WINDOW_MM
    zc = np.arange(dz_mm / 2.0, z2, dz_mm)
    zp = zc[zc > z1]
    e_c, e_p = water.energy_at_depth(float(e0), zc), water.energy_at_depth(float(e0), zp)
    s_p = water.stopping_mev_per_cm(e_p)
    denom = float(np.sum(s_p) * dz_mm / 10.0)
    g_path = g_rec = g_slow = g_stop = 0.0
    for name in ("O-16", "H-1"):
        for i, (zi, e) in enumerate(zip(zc, e_c, strict=True)):
            ns = sigma_nsigma(sh, water, name, float(e)) * dz_mm / 10.0
            if ns <= 0.0:
                continue
            ev = primary_event(sh, name, kind, float(e), U)
            cos, r1 = ev["cos"], np.asarray(water.range_mm(ev["t1"]))
            later = zp > zi
            if later.any():
                dist = zp[later][:, None] - zi
                reach = (cos[None, :] > 0.05) & (r1[None, :] * cos[None, :] >= dist)
                net = np.mean(
                    np.where(reach, 1.0 / np.maximum(cos[None, :], 0.05), 0.0) - 1.0, axis=1
                )
                g_path += ns * float(np.sum(s_p[later] * net) * dz_mm / 10.0)
                # with the stopping power of the continuing proton at z' (residual range R(T1) - path)
                res = np.maximum(r1[None, :] - dist / np.maximum(cos[None, :], 0.05), 0.05)
                e_after = np.asarray(water.stop.energy_from_range(res * water.rho / 10.0))
                s_after = water.stopping_mev_per_cm(e_after)
                term = (
                    np.where(reach, s_after / np.maximum(cos[None, :], 0.05), 0.0)
                    - s_p[later][:, None]
                )
                g_stop += ns * float(np.sum(np.mean(term, axis=1)) * dz_mm / 10.0)
            if zi > z1:
                if name == "O-16":
                    g_rec += ns * float(np.mean(ev["t_other"]))
                else:
                    g_slow += ns * float(np.mean(ev["t_other"]))
    g_path, g_rec, g_slow, g_stop = g_path / denom, g_rec / denom, g_slow / denom, g_stop / denom
    return {"G_path": g_path, "G_recoil": g_rec, "G_pp_slow_local_upper": g_slow, "G_total": g_path + g_rec + g_slow,
            "G_path_plus_recoil": g_path + g_rec, "G_path_with_post_event_stopping": g_stop,
            "G_total_with_post_event_stopping": g_stop + g_rec + g_slow,
            "post_event_stopping_term": g_stop - g_path}  # fmt: skip


def f_pred(k_star: float, g: float) -> float:
    """Kernel prediction of F = (peak/plateau)_EM+el / (peak/plateau)_EM: ``(1 - K*) / (1 + G)``."""
    return (1.0 - k_star) / (1.0 + g)


# -- the reviewer's crude estimate, recomputed with its own definitions ---------------------------------
def crude_reviewer_estimate(sh: Shapes, e0: int, kind: str) -> dict[str, float]:
    """Re-implementation of ``review-scripts/removal.py`` (Bragg-Kleeman R = 0.0022 E^1.77 cm, straight line, p+O only,
    step threshold, no straggling, no p-p, nearest-grid radius and sigma, non-relativistic recoil) for shape ``kind``
    (``sharp`` = analytic disk at the nearest-grid radius)."""
    from ionmc.nuclear.bgg import M_PROTON_MEV as mp

    it = sh.index("O-16")
    grid = np.asarray(sh.a["grid_e_mev"], dtype=float)
    alpha, p = 0.0022, 1.77
    rng = lambda e: alpha * np.maximum(e, 0.0) ** p  # noqa: E731
    n_o = 3.343e22
    mu = np.linspace(-1.0, 1.0, 4001)
    m_t = float(sh.a["target_mass_mev"][it])
    r0 = float(rng(e0))
    zs = np.linspace(0.0, r0 * R_FRACTION, 120)
    acc_p = acc_tr = 0.0
    accw = dict.fromkeys(W_MM, 0.0)
    for j in range(len(zs) - 1):
        z, dz = 0.5 * (zs[j] + zs[j + 1]), zs[j + 1] - zs[j]
        e = ((r0 - z) / alpha) ** (1.0 / p)
        if e < 7:
            continue
        i = int(np.argmin(abs(grid - e)))
        sig = float(sh.a["sigma_barn"][it, i]) * 1e-24
        if sig <= 0:
            continue
        k = float(sh.a["p_cm_mev"][it, i]) / 197.327
        r = float(sh.a["radius_fm"][it, i])
        q = 2.0 * k * np.sqrt((1.0 - mu) / 2.0)
        w = _crude_weight(kind, q, r)
        w = w / np.trapezoid(w, mu)
        thl = np.arctan2(np.sqrt(1.0 - mu**2), mu + mp / m_t)
        t_rec = 2.0 * mp * m_t / (mp + m_t) ** 2 * e * (1.0 - mu)
        cosl = np.cos(thl)
        defic = (rng(e) - rng(e - t_rec) * np.clip(cosl, 0.0, None)) * 10.0
        d_p = n_o * sig * dz
        acc_p += d_p
        acc_tr += d_p * np.trapezoid(w * (1.0 - cosl), mu)
        for x in W_MM:
            accw[x] += d_p * np.trapezoid(w * (defic > x), mu)
    return {"P_el": acc_p, "sum_1_minus_cos": acc_tr, **{f"removal_{x:g}mm": accw[x] for x in W_MM}}


def _crude_weight(kind: str, q: NDArray[np.float64], r: float) -> NDArray[np.float64]:
    from .de1_shapes import family_weight

    return family_weight("disk" if kind == "sharp" else kind, np.maximum(q, 1e-12), r)


# -- H6 case-file diff --------------------------------------------------------------------------------
def h6_case_diff() -> dict[str, Any]:
    base = c.REPO / "validation" / "reference_cases" / "topas"
    out: dict[str, Any] = {"pairs": []}
    ok_all = True
    for e in (150, 200):
        for seed in (1, 2, 3):
            a, b = (
                f"proton-water-{e}mev-idd-r20-emelastic-seed{seed}",
                f"proton-water-{e}mev-idd-r20-emonly-seed{seed}",
            )
            rec: dict[str, Any] = {"emelastic": a, "emonly": b, "files": {}}
            for fn in ("input.txt", "case.json"):
                la = (base / a / fn).read_text().splitlines()
                lb = (base / b / fn).read_text().splitlines()
                diff = [
                    ln
                    for ln in difflib.unified_diff(lb, la, lineterm="", n=0)
                    if ln[:1] in "+-" and ln[:3] not in ("+++", "---")
                ]
                rec["files"][fn] = diff
            allowed = ("# Native OpenTOPAS", "i:Ts/Seed", "sv:Ph/Default/Modules")
            unexpected = [
                ln for ln in rec["files"]["input.txt"] if not ln[1:].lstrip().startswith(allowed)
            ]
            rec["input_txt_unexpected_lines"] = unexpected
            rec["input_txt_only_modules_seed_header"] = not unexpected
            ok_all &= not unexpected
            out["pairs"].append(rec)
    out["all_input_txt_differ_only_in_header_seed_modules"] = ok_all
    out["files_hashed"] = {rel_: c.sha256_file(c.REPO / rel_) for rel_ in _h6_files(base)}
    return out


def _h6_files(base: Any) -> list[str]:
    fs = []
    for e in (150, 200):
        for mode in ("emelastic", "emonly"):
            for seed in (1, 2, 3):
                for fn in ("input.txt", "case.json"):
                    fs.append(c.rel(base / f"proton-water-{e}mev-idd-r20-{mode}-seed{seed}" / fn))
    return fs


# -- orchestration -------------------------------------------------------------------------------------
def equivalent_shape(
    sh: Shapes, water: c.Water, e0: int, target_k: float, param: str
) -> float | None:
    """Report-only: the surface ``s`` (fm) or radius factor ``Rx`` for which ``K*`` (shifted form; O + H-1) equals
    ``target_k`` (the table value scaled by the observed TOPAS/ionmc peak-removal ratio). TOPAS-seen diagnostic;
    never a model parameter."""

    def k_of(x: float) -> float:
        kind = f"s{x}" if param == "s" else f"Rx{x}"
        k = kernels(sh, water, e0, kind, ("O-16", "H-1"), True, n_z=120)
        return float(sum(v["K_star_shifted"] for v in k.values()))

    lo, hi = (0.0, 1.5) if param == "s" else (1.0, 1.6)
    klo, khi = k_of(lo), k_of(hi)
    if not (min(klo, khi) <= target_k <= max(klo, khi)):
        return None
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        km = k_of(mid)
        if (km - target_k) * (klo - target_k) > 0:
            lo, klo = mid, km
        else:
            hi = mid
    return 0.5 * (lo + hi)


def run() -> None:
    tab = c.load_table()
    sh, water = Shapes(tab), c.Water(tab)
    ref = c.v11_reference()
    kinds = ("table", *FAMILIES)
    point: dict[str, Any] = {}
    for e in ENERGIES:
        point[str(e)] = {}
        for name in ("O-16", "C-12", "H-1"):
            ks = (
                ("table",)
                if name != "O-16" and name != "C-12"
                else (("table", "disk") if name == "C-12" else kinds)
            )
            point[str(e)][name] = {k: pointwise(sh, water, name, k, float(e)) for k in ks}
        point[str(e)]["O-16"]["table"]["first_minimum_bin_density_deg"] = (
            table_bin_density_first_minimum_deg(sh, "O-16", float(e))
        )
        point[str(e)]["C-12"]["table"]["first_minimum_bin_density_deg"] = (
            table_bin_density_first_minimum_deg(sh, "C-12", float(e))
        )
    kern: dict[str, Any] = {}
    for e in ENERGIES:
        kern[str(e)] = {}
        have_idd = e in (150, 200)
        for k in ("table", *FAMILIES) if have_idd or e in (50, 100) else ("table",):
            if e in (50, 100) and k not in ("table", "disk"):
                continue
            kk = kernels(sh, water, e, k, ("O-16", "H-1"), have_idd)
            if k == "table":
                kk["C-12 (hypothetical, n=1 per cm3 along the water path)"] = kernels(
                    sh, water, e, "table", ("C-12",), have_idd
                )["C-12"]
            if e >= 100:
                kk["G"] = plateau_gain(sh, water, e, k)
            kern[str(e)][k] = kk
    recon = _reconcile(sh, water, ref, kern)
    quarantine(sh, water, ref, recon)
    crude = {
        str(e): {
            k: crude_reviewer_estimate(sh, e, k)
            for k in ("sharp", "s0.5", "s0.9", "Rx1.2", "gauss")
        }
        for e in (150, 200)
    }
    h6 = h6_case_diff()
    doc = {"banner": c.BANNER, "table_id": c.TABLE_ID, "energies_mev": list(ENERGIES), "w_mm": list(W_MM),
           "pointwise": point, "kernels": kern, "reconciliation": recon, "crude_recomputed": crude,
           "crude_cited_review_150": {k: list(v) for k, v in REVIEW_CRUDE.items()}, "h6": h6,
           "v11_reference": {str(k): v for k, v in ref.items()}}  # fmt: skip
    inputs = [
        *c.table_files(),
        c.RUN_V11,
        c.RUN_V10A,
        *(c.V11_PARTIALS[k] for k in sorted(c.V11_PARTIALS)),
    ]
    inputs += [c.REPO / r for r in sorted(h6["files_hashed"])]
    c.write_outputs(c.OUT / "de0", "kernels.json", doc, "de0-report.md", _report(doc, inputs))


def _total(kern: dict[str, Any], key: str) -> float:
    return float(sum(kern[t][key] for t in ("O-16", "H-1")))


def theta_rms_s(sh: Shapes, water: c.Water, e: int, target_deg: float) -> float:
    """Surface ``s`` (fm) of the ``s<value>`` family at the table radius whose p+O ``theta_rms,CM`` equals the planning
    value of decision 0041 (17 / 14 deg): a planning-moment match, no V11 or TOPAS number involved."""
    lo, hi = 0.0, 1.5
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if pointwise(sh, water, "O-16", f"s{mid:.6f}", float(e))["theta_rms_cm_deg"] > target_deg:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _kg(sh: Shapes, water: c.Water, e: int, kind: str) -> tuple[float, float, float]:
    k = kernels(sh, water, e, kind, ("O-16", "H-1"), True)
    g = plateau_gain(sh, water, e, kind)["G_total"]
    ks = _total(k, "K_star_shifted")
    return ks, g, f_pred(ks, g)


def _reconcile(
    sh: Shapes, water: c.Water, ref: dict[int, Any], kern: dict[str, Any]
) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for e in (150, 200):
        r = ref[e]
        kt = kern[str(e)]["table"]
        k_lit, k_sh, g = _total(kt, "K_star_literal"), _total(kt, "K_star_shifted"), kt["G"]
        obs_removal_ionmc = 1.0 - r["F_ionmc"] * (1.0 + r["gain_ionmc"])
        obs_removal_topas = 1.0 - r["F_topas"] * (1.0 + r["gain_topas"])
        f_plan = EXPECTATION[e] * r["F_topas"]
        r_tab, plan_r = sh.r_pcm("O-16", float(e))[0], PLANNING_R_O_FM[e]
        s_plan = theta_rms_s(sh, water, e, PLANNING_THETA_RMS_O[e])
        f_tab = f_pred(k_sh, g["G_total"])
        radius = _kg(sh, water, e, f"Rx{plan_r / r_tab:.6f}")
        shape = _kg(sh, water, e, f"s{s_plan:.6f}")
        gap_obs = f_plan - r["F_ionmc"]
        out[str(e)] = {
            "observed": {"F_ionmc": r["F_ionmc"], "F_topas": r["F_topas"], "gain_ionmc": r["gain_ionmc"], "gain_topas": r["gain_topas"],
                         "ratio_F_ionmc_over_topas": r["F_ionmc"] / r["F_topas"],
                         "peak_removal_ionmc_obs": obs_removal_ionmc, "peak_removal_topas_obs": obs_removal_topas},
            "planning_expectation": {"ratio": EXPECTATION[e], "implied_F_ionmc": f_plan,
                                     "implied_peak_removal_with_planning_gain": 1.0 - f_plan * (1.0 + PLANNING_GAIN[e]),
                                     "planning_gain_used": PLANNING_GAIN[e], "gap_to_observed_F_ionmc": gap_obs},
            "kernel_table": {"K_star_as_frozen_sign_error": k_lit, "K_star_shifted": k_sh, "G_total": g["G_total"],
                             "G_total_with_post_event_stopping": g["G_total_with_post_event_stopping"],
                             "F_pred_shifted": f_tab, "F_pred_over_F_obs": f_tab / r["F_ionmc"],
                             "peak_removal_pred_over_obs": k_sh / obs_removal_ionmc},
            "planning_radius_disk": {"R_O_fm": plan_r, "R_table_fm": r_tab, "K_star_shifted": radius[0], "G": radius[1], "F_pred": radius[2],
                                     "dF_vs_table": radius[2] - f_tab},
            "planning_theta_rms_shape": {"theta_rms_cm_deg_target": PLANNING_THETA_RMS_O[e], "s_fm_matching": s_plan,
                                         "K_star_shifted": shape[0], "G": shape[1], "F_pred": shape[2], "dF_vs_table": shape[2] - f_tab},
        }  # fmt: skip
        d = out[str(e)]
        d["miss_explained"] = {
            "radius_fraction_of_gap_to_observed": d["planning_radius_disk"]["dF_vs_table"] / gap_obs,
            "shape_fraction_of_gap_to_observed": d["planning_theta_rms_shape"]["dF_vs_table"] / gap_obs,
            "shape_fraction_of_gap_to_kernel_table": d["planning_theta_rms_shape"]["dF_vs_table"] / (f_plan - f_tab),
            "gap_to_kernel_table": f_plan - f_tab,
        }  # fmt: skip
    return out


def quarantine(sh: Shapes, water: c.Water, ref: dict[int, Any], recon: dict[str, Any]) -> None:
    """Data-seen (V11/TOPAS) equivalent-shape search: written only to ``de0/quarantine-data-seen.json``."""
    doc: dict[str, Any] = {
        "banner": "QUARANTINED: fitted to V11/TOPAS numbers, not pre-registered. Must never be cited, "
        "used in a report, a decision or a selection. " + c.BANNER
    }
    for e in (150, 200):
        r = recon[str(e)]
        target = (
            r["kernel_table"]["K_star_shifted"]
            * r["observed"]["peak_removal_topas_obs"]
            / r["observed"]["peak_removal_ionmc_obs"]
        )
        doc[str(e)] = {
            "s_fm": equivalent_shape(sh, water, e, target, "s"),
            "radius_factor": equivalent_shape(sh, water, e, target, "Rx"),
        }
    (c.OUT / "de0").mkdir(parents=True, exist_ok=True)
    (c.OUT / "de0" / "quarantine-data-seen.json").write_text(
        json.dumps(doc, indent=1, sort_keys=True) + "\n"
    )


def _g(x: Any, f: str = ".4g") -> str:
    return "-" if x is None else format(x, f)


def _report(doc: dict[str, Any], inputs: list[Any]) -> str:
    cmd = "python -m validation.scripts.elastic_attribution.de0_kernels"
    k, p, rc = doc["kernels"], doc["pointwise"], doc["reconciliation"]
    out = [
        c.header(
            "D-E0: deterministic kernels of the elastic channel",
            cmd,
            inputs,
            [
                "Definitions: module docstring of `de0_kernels.py` (plan section 2.1). Water, Bethe CSDA range (the source of the V5/V11 runs), "
                "straight-line path, exact two-body kinematics, no straggling, 2048 equiprobable quantile points per event, 300 depth nodes. "
                "The plan's G wording is under-specified; the operational form is in the docstring."
            ],
        )
    ]
    out += [
        "",
        "## Dated note (2026-10-14): sign of the frozen K* formula",
        "",
        "The plan (section 2.1, frozen by Amendment 19 (a)) writes K* = int n sigma <1 - IDD_EM(z_pk - Delta)/IDD_EM(z_pk)> dz. This formula has a "
        "sign error: a proton whose range is shortened by Delta has the depth-dose curve IDD_EM(z + Delta), so its contribution at z_pk is "
        "IDD_EM(z_pk + Delta) (distal fall-off); IDD_EM(z_pk - Delta) is the proximal rise, i.e. a range *extended* by Delta. The plan's own prose "
        "('moves that proton's contribution ... to z_pk - Delta') describes the shifted form. The as-frozen form is computed and labelled "
        "'as frozen (sign error)'; **no conclusion is drawn from it**. The shifted form is the physically meaningful kernel "
        "(F predicted 0.853 / 0.832 against 0.857 / 0.844 observed in V11; peak removal about 10 % high). In D-E3 both forms are to be reported and "
        "validated against the shifted one. The frozen plan is not edited. The kernel is adequate for ranking and sensitivity to about 10 %, not for a 0.010 closure.",
        "",
        "G caveat: G (operational form, see the docstring) omits the higher stopping power of the proton that continues after a p-p event "
        "(E_fast = E - T_slow); that term is positive. The column `G_total incl. post-p-p stopping` adds the post-event stopping power of every "
        "continuing proton at z' (all events, exact residual range); its size is in `post_event_stopping_term` of kernels.json. "
        "G is still the plan's incomplete descriptive quantity and enters only the descriptive peak-versus-gain ordering.",
        "",
        "## Kernels in water (p+O and p-p; per primary)",
        "",
        "| E | shape | P_el(O) | P_el(pp) | K_1mm | K_2mm | K_3mm | K_5mm | K* as frozen (sign error; no conclusion) | K* shifted (meaningful) | K*_gauss (shifted form) | G_total | G_path | G_recoil | G_pp_slow(up.) | G_total incl. post-p-p stopping |",
        "|" + "---|" * 16,
    ]
    for e in ENERGIES:
        for kind, kk in k[str(e)].items():
            g = kk.get("G", {})
            tot = (
                lambda key: sum(kk[t].get(key, 0.0) for t in ("O-16", "H-1"))
                if key in kk["O-16"]
                else None
            )  # noqa: E731
            out.append(
                f"| {e} | {kind} | {_g(kk['O-16']['P_el'])} | {_g(kk['H-1']['P_el'])} | "
                + " | ".join(_g(tot(f"K_{w:g}mm")) for w in W_MM)
                + f" | {_g(tot('K_star_literal'))} | {_g(tot('K_star_shifted'))} | {_g(tot('K_star_gauss'))} | {_g(g.get('G_total'))} | {_g(g.get('G_path'))} | {_g(g.get('G_recoil'))} | {_g(g.get('G_pp_slow_local_upper'))} | {_g(g.get('G_total_with_post_event_stopping'))} |"
            )
    out += [
        "",
        "Per target (table shape): K_w of p+O alone and p-p alone, and the conditional removal probability K_w/P_el.",
        "",
        "| E | target | K_1mm | K_2mm | K_3mm | K_5mm | K_2mm/P_el | K_3mm/P_el |",
        "|" + "---|" * 8,
    ]
    for e in ENERGIES:
        kt = k[str(e)]["table"]
        for t in ("O-16", "H-1", "C-12 (hypothetical, n=1 per cm3 along the water path)"):
            d = kt[t]
            cols = ["-"] * 4 if t.startswith("C-12") else [_g(d[f"K_{w:g}mm"]) for w in W_MM]
            out.append(
                f"| {e} | {t[:4]} | "
                + " | ".join(cols)
                + f" | {_g(d['K_2mm'] / d['P_el'])} | {_g(d['K_3mm'] / d['P_el'])} |"
            )
    out += [
        "",
        "## Pointwise shape moments at the entrance energy (table; per interaction)",
        "",
        "| E | target | <1-mu_CM> | theta_rms,CM | theta_med,lab | theta_2mm | P(th>th_2mm) | theta_3mm | P(th>th_3mm) | T (3 th_med) | <T_other> MeV | P(D>2mm) | P(D>3mm) | first min (bin density) |",
        "|" + "---|" * 14,
    ]
    for e in ENERGIES:
        for name in ("O-16", "C-12", "H-1"):
            d = p[str(e)][name]["table"]
            out.append(
                f"| {e} | {name} | {_g(d['mean_1_minus_mu_cm'])} | {_g(d['theta_rms_cm_deg'], '.1f')} | {_g(d['theta_med_lab_deg'], '.2f')} | {_g(d['theta_2mm_deg'], '.2f')} | {_g(d['P_theta_lab_gt_theta_2mm'])} | {_g(d['theta_3mm_deg'], '.2f')} | {_g(d['P_theta_lab_gt_theta_3mm'])} | {_g(d['T_tail_moment'])} | {_g(d['mean_t_other_mev'])} | {_g(d['P_deficit_gt_2mm'])} | {_g(d['P_deficit_gt_3mm'])} | {_g(d.get('first_minimum_bin_density_deg'), '.1f')} |"
            )
    out += [
        "",
        "## Pointwise: p+O shape families",
        "",
        "| E | shape | <1-mu_CM> | theta_rms,CM | theta_med,lab | P(th>th_2mm) | P(th>th_3mm) | T | P(D>2mm) | P(D>3mm) |",
        "|" + "---|" * 10,
    ]
    for e in ENERGIES:
        for kind, d in p[str(e)]["O-16"].items():
            out.append(
                f"| {e} | {kind} | {_g(d['mean_1_minus_mu_cm'])} | {_g(d['theta_rms_cm_deg'], '.1f')} | {_g(d['theta_med_lab_deg'], '.2f')} | {_g(d['P_theta_lab_gt_theta_2mm'])} | {_g(d['P_theta_lab_gt_theta_3mm'])} | {_g(d['T_tail_moment'])} | {_g(d['P_deficit_gt_2mm'])} | {_g(d['P_deficit_gt_3mm'])} |"
            )
    out += ["", "## Reconciliation of the planning expectation (Amendment 14 (f): 0.985 / 0.986) with the observed ratio", "",
            "Observed V5 ratio ionmc/TOPAS peak/plateau: 0.9634 / 0.9511 (step 22, -3.66 % / -4.89 %); V11 F ratio below. "
            "Kernel prediction F = (1 - K*_shifted)/(1 + G), O + p-p, table shape. The as-frozen K* is shown for the record only (sign error).", "",
            "| quantity | 150 MeV | 200 MeV |", "|---|---|---|"]  # fmt: skip
    rows = [
        ("F ionmc (V11, archived)", lambda r: r["observed"]["F_ionmc"]),
        ("F TOPAS (V11 rows 10/11, archived)", lambda r: r["observed"]["F_topas"]),
        ("F ionmc / F TOPAS", lambda r: r["observed"]["ratio_F_ionmc_over_topas"]),
        ("peak removal ionmc = 1 - F(1+gain)", lambda r: r["observed"]["peak_removal_ionmc_obs"]),
        ("peak removal TOPAS = 1 - F(1+gain)", lambda r: r["observed"]["peak_removal_topas_obs"]),
        ("planning expectation ratio", lambda r: r["planning_expectation"]["ratio"]),
        (
            "planning-implied F ionmc (= ratio x F TOPAS)",
            lambda r: r["planning_expectation"]["implied_F_ionmc"],
        ),
        (
            "gap: planning-implied F - observed F ionmc",
            lambda r: r["planning_expectation"]["gap_to_observed_F_ionmc"],
        ),
        (
            "K* as frozen (sign error), table (no conclusion)",
            lambda r: r["kernel_table"]["K_star_as_frozen_sign_error"],
        ),
        ("K* shifted, table (O + pp)", lambda r: r["kernel_table"]["K_star_shifted"]),
        (
            "G total, table / incl. post-p-p stopping",
            lambda r: f"{r['kernel_table']['G_total']:.4f} / {r['kernel_table']['G_total_with_post_event_stopping']:.4f}",
        ),
        ("gain ionmc observed", lambda r: r["observed"]["gain_ionmc"]),
        ("F predicted, table (shifted K*)", lambda r: r["kernel_table"]["F_pred_shifted"]),
        ("F predicted / F observed ionmc", lambda r: r["kernel_table"]["F_pred_over_F_obs"]),
        (
            "peak removal predicted / observed",
            lambda r: r["kernel_table"]["peak_removal_pred_over_obs"],
        ),
        (
            "R_O table (fm) / planning R_O (fm)",
            lambda r: f"{r['planning_radius_disk']['R_table_fm']:.3f} / {r['planning_radius_disk']['R_O_fm']}",
        ),
        (
            "RADIUS: F predicted at planning R_O (disk)",
            lambda r: r["planning_radius_disk"]["F_pred"],
        ),
        ("RADIUS: dF vs table", lambda r: r["planning_radius_disk"]["dF_vs_table"]),
        (
            "SHAPE: s (fm) giving the planning theta_rms,CM (17 / 14 deg)",
            lambda r: r["planning_theta_rms_shape"]["s_fm_matching"],
        ),
        ("SHAPE: F predicted at that s", lambda r: r["planning_theta_rms_shape"]["F_pred"]),
        ("SHAPE: dF vs table", lambda r: r["planning_theta_rms_shape"]["dF_vs_table"]),
        (
            "radius share of the gap (planning F - observed F)",
            lambda r: r["miss_explained"]["radius_fraction_of_gap_to_observed"],
        ),
        (
            "shape share of the gap (planning F - observed F)",
            lambda r: r["miss_explained"]["shape_fraction_of_gap_to_observed"],
        ),
        (
            "shape share of the gap to the kernel's own table prediction",
            lambda r: r["miss_explained"]["shape_fraction_of_gap_to_kernel_table"],
        ),
    ]
    for lab, fn in rows:
        out.append(
            f"| {lab} | "
            + " | ".join(
                _g(fn(rc[str(e)])) if not isinstance(fn(rc[str(e)]), str) else fn(rc[str(e)])
                for e in (150, 200)
            )
            + " |"
        )
    out += [
        "",
        "## Reviewer's crude removal estimate (150 / 200 MeV, p+O only): cited against recomputed",
        "",
        "Cited (review.md B.1, 150 MeV): moment, removal 2 mm, 3 mm. Recomputed with the same definitions "
        "(`crude_reviewer_estimate`: Bragg-Kleeman, straight line, step threshold) and with the plan definitions "
        "(project Bethe range, exact kinematics, table-interpolated radius, K_w of the O channel).",
        "",
        "| E | shape | cited moment / 2 mm / 3 mm | crude recomputed | plan K_w O-only 2 mm / 3 mm | relative change vs sharp: crude 2mm / 3mm / moment | plan 2mm / 3mm |",
        "|" + "---|" * 7,
    ]
    crude = doc["crude_recomputed"]
    for e in (150, 200):
        base_c, base_p = crude[str(e)]["sharp"], k[str(e)]["disk"]["O-16"]
        for lab, kd in (
            ("sharp", "disk"),
            ("s0.5", "s0.5"),
            ("s0.9", "s0.9"),
            ("Rx1.2", "Rx1.2"),
            ("gauss", "gauss"),
        ):
            cc, pp = crude[str(e)][lab], k[str(e)][kd]["O-16"]
            cited = REVIEW_CRUDE[lab] if e == 150 else None
            out.append(
                f"| {e} | {lab} | {'-' if cited is None else '/'.join(map(str, cited))} | {cc['sum_1_minus_cos']:.4f} / {cc['removal_2mm']:.4f} / {cc['removal_3mm']:.4f} | {pp['K_2mm']:.4f} / {pp['K_3mm']:.4f} | "
                f"{cc['removal_2mm'] / base_c['removal_2mm'] - 1:+.1%} / {cc['removal_3mm'] / base_c['removal_3mm'] - 1:+.1%} / {cc['sum_1_minus_cos'] / base_c['sum_1_minus_cos'] - 1:+.1%} | "
                f"{pp['K_2mm'] / base_p['K_2mm'] - 1:+.1%} / {pp['K_3mm'] / base_p['K_3mm'] - 1:+.1%} |"
            )
    h6 = doc["h6"]
    out += [
        "",
        "## H6: TOPAS rows 10/11 (emelastic) against the frozen EM-only cases (emonly)",
        "",
        f"`input.txt` differs only in the header comment, the seed and the module list for all 6 pairs: **{h6['all_input_txt_differ_only_in_header_seed_modules']}**.",
        "",
        "Differing lines of the seed-1 pairs (emelastic `+` against emonly `-`):",
        "",
    ]
    for rec in h6["pairs"]:
        if rec["emelastic"].endswith("seed1"):
            out += [f"- {rec['emelastic']} vs {rec['emonly']}:"] + [
                f"  - `{ln[:150]}`" for ln in rec["files"]["input.txt"]
            ]
    out += [
        "",
        "`case.json` differs in the descriptive text, row/mode/modules and the seed (not scored by the engine). Hashed inputs are listed above.",
    ]
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    run()
