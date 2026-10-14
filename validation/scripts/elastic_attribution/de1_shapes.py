"""D-E1 (exploratory, report-only; Amendment 19): the pinned p+A shape and the Amendment 19 (b) shape
families against the V10-A data, plus the shape-family machinery used by D-E0.

Per V10-A dataset: model against data ``d sigma / d Omega_CM`` at the data angles, the data-integrated
partial ``sigma_tr`` over the measured range (10 deg to the last data angle), the same for the model and
for the families at their a-priori parameters (no fit), the location of the measured minima against the
table's first minimum, and a one-parameter R diagnostic fit that is **data-seen** and never used for a
model parameter.

Run: ``python -m validation.scripts.elastic_attribution.de1_shapes``.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from numpy.typing import NDArray

from . import common as c

FAMILIES = ("disk", "Rx1.2", "s0.5", "s0.9", "gauss")
"""Analytic p+A shapes (Amendment 19 (b)): ``disk`` = black disk at the table radius R(E) (the form of the
pinned table before its 256-quantile tabulation), ``Rx1.2`` = disk at 1.2 R, ``s0.5`` / ``s0.9`` = disk times
``exp(-q^2 s^2)`` (s in fm), ``gauss`` = ``exp(-q^2 R^2 / 4)``. Normalised over the full sphere; sigma unchanged."""
HBARC = 197.3269804
J1_FIRST_ZERO = 3.8317059702075125
N_S = 16385

_J1_X = np.linspace(0.0, 120.0, 60001)
_J1_Y: NDArray[np.float64] | None = None


def _j1(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """J1 by interpolation in a table of the repository's ``bessel_j1`` (error < 1e-6)."""
    global _J1_Y
    if _J1_Y is None:
        from ionmc.nuclear.elastic_build import bessel_j1

        _J1_Y = np.asarray(bessel_j1(_J1_X))
    return np.interp(np.asarray(x, dtype=float), _J1_X, _J1_Y)


def disk_a(x: NDArray[np.float64]) -> NDArray[np.float64]:
    """``|2 J1(x)/x|^2`` (1 at x = 0)."""
    xa = np.asarray(x, dtype=float)
    safe = np.where(xa < 1e-8, 1.0, xa)
    a = np.where(xa < 1e-8, 1.0, 2.0 * _j1(safe) / safe)
    return np.asarray(a * a)


def family_weight(kind: str, q_fm: NDArray[np.float64], r_fm: float) -> NDArray[np.float64]:
    """Unnormalised angular weight of shape family ``kind`` at momentum transfer ``q_fm`` [1/fm]. Besides the five
    names of :data:`FAMILIES`, ``s<value>`` (surface ``s`` in fm) and ``Rx<factor>`` are accepted (used for the
    report-only equivalent-shape search of D-E0)."""
    if kind == "disk":
        return disk_a(q_fm * r_fm)
    if kind.startswith("Rx"):
        return disk_a(q_fm * float(kind[2:]) * r_fm)
    if kind.startswith("s"):
        s = float(kind[1:])
        return np.asarray(disk_a(q_fm * r_fm) * np.exp(-((q_fm * s) ** 2)))
    if kind == "gauss":
        return np.asarray(np.exp(-((q_fm * r_fm) ** 2) / 4.0))
    raise ValueError(kind)


class Shapes:
    """The pinned table plus the analytic families, at the table radius and CM momentum."""

    def __init__(self, table: Any) -> None:
        self.tab = table
        self.a = table.arrays
        self.names = table.target_names
        self._s = np.linspace(0.0, 1.0, N_S)

    def index(self, name: str) -> int:
        return int(self.names.index(name))

    def r_pcm(self, name: str, e_mev: float) -> tuple[float, float]:
        it = self.index(name)
        g = np.asarray(self.a["grid_e_mev"], dtype=float)
        return (
            float(np.interp(e_mev, g, self.a["radius_fm"][it])),
            float(np.interp(e_mev, g, self.a["p_cm_mev"][it])),
        )

    def _w_s(
        self, kind: str, name: str, e_mev: float, r_override: float | None = None
    ) -> NDArray[np.float64]:
        r, pcm = self.r_pcm(name, e_mev)
        if r_override is not None:
            r = r_override
        return family_weight(kind, 2.0 * pcm / HBARC * self._s, r)

    def norm(self, w_s: NDArray[np.float64]) -> float:
        """``int w dmu`` over [-1, 1] (``dmu = 4 s ds``)."""
        g = 4.0 * self._s * w_s
        return float(np.sum(0.5 * (g[1:] + g[:-1]) * np.diff(self._s)))

    def density_mb_sr(
        self,
        kind: str,
        name: str,
        e_mev: float,
        theta_deg: NDArray[np.float64],
        r_override: float | None = None,
    ) -> NDArray[np.float64]:
        """Nuclear ``d sigma / d Omega_CM`` [mb/sr] of family ``kind`` normalised to the table sigma_el."""
        r, pcm = self.r_pcm(name, e_mev)
        if r_override is not None:
            r = r_override
        w_s = self._w_s(kind, name, e_mev, r)
        sig_mb = float(self.tab.sigma_barn(self.index(name), e_mev)) * 1.0e3
        s = np.sin(np.radians(np.asarray(theta_deg, dtype=float)) / 2.0)
        w = family_weight(kind, 2.0 * pcm / HBARC * s, r)
        return np.asarray(sig_mb * w / (2.0 * math.pi * self.norm(w_s)))

    def mu_quantiles(
        self, kind: str, name: str, e_mev: float, u: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """``mu_CM`` at the cumulative fractions ``u`` (ascending ``mu``) of ``kind`` (``table`` = the pinned
        256-quantile table, otherwise the analytic family)."""
        if kind == "table":
            return np.asarray(self.tab.sample_mu(self.index(name), e_mev, u))
        w_s = self._w_s(kind, name, e_mev)
        g = 4.0 * self._s * w_s
        f = np.concatenate(([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(self._s))))
        f /= f[-1]  # fraction of events with s' < s (forward ones)
        s_u = np.interp(1.0 - np.asarray(u, dtype=float), f, self._s)
        return np.asarray(1.0 - 2.0 * s_u * s_u)


def first_local_min_deg(
    theta: NDArray[np.float64], f: NDArray[np.float64], above: float = 10.0
) -> float | None:
    for i in range(1, theta.size - 1):
        if theta[i] > above and f[i] < f[i - 1] and f[i] < f[i + 1]:
            return float(theta[i])
    return None


def _trapz(x: NDArray[np.float64], y: NDArray[np.float64]) -> float:
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def partial_moments(
    theta_deg: NDArray[np.float64], f_mb_sr: NDArray[np.float64]
) -> dict[str, float]:
    """Trapezoid (in theta) partial ``sigma`` and ``sigma_tr = 2 pi int (1 - cos) f sin`` [mb] over the nodes."""
    t = np.radians(theta_deg)
    return {"sigma_mb": 2.0 * math.pi * _trapz(t, f_mb_sr * np.sin(t)),
            "sigma_tr_mb": 2.0 * math.pi * _trapz(t, f_mb_sr * np.sin(t) * (1.0 - np.cos(t)))}  # fmt: skip


def merge_duplicates(
    th: NDArray[np.float64], y: NDArray[np.float64], rel: NDArray[np.float64], tol: float = 0.3
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Merge points whose angles differ by less than ``tol`` deg (the archived sets repeat angles from two runs):
    angle mean, value mean, relative error of the mean (1/sqrt(n) of the mean relative error; NaN if absent)."""
    t_o, y_o, r_o = [], [], []
    i = 0
    while i < th.size:
        j = i
        while j + 1 < th.size and th[j + 1] - th[i] < tol:
            j += 1
        n = j - i + 1
        t_o.append(float(np.mean(th[i : j + 1])))
        y_o.append(float(np.mean(y[i : j + 1])))
        r_o.append(float(np.mean(rel[i : j + 1])) / math.sqrt(n))
        i = j + 1
    return np.array(t_o), np.array(y_o), np.array(r_o)


def measured_minima(
    th: NDArray[np.float64],
    y: NDArray[np.float64],
    rel: NDArray[np.float64],
    k_cm: float,
    significant: bool,
) -> list[dict[str, Any]]:
    """Measured local minima above 10 deg. ``significant=False``: the V10-A rule (a point lower than both of its
    neighbours; duplicates not merged). ``significant=True``: duplicates merged, and the point must lie at least 10 %
    below both neighbours and by at least 3 combined errors below the lower neighbour. Per minimum: the neighbour
    ratios, the depth in sigma and the disk radius ``3.8317 / q`` the position implies."""
    if significant:
        th, y, rel = merge_duplicates(th, y, rel)
    out = []
    for i in range(1, th.size - 1):
        if th[i] > 10.0 and y[i] < y[i - 1] and y[i] < y[i + 1]:
            q = 2.0 * k_cm * math.sin(math.radians(th[i]) / 2.0)
            lo = i - 1 if y[i - 1] <= y[i + 1] else i + 1
            err_i, err_n = rel[i] * y[i], rel[lo] * y[lo]
            sig = (
                (y[lo] - y[i]) / math.hypot(err_i, err_n)
                if np.isfinite(err_i) and np.isfinite(err_n)
                else None
            )
            ratios = [float(y[i - 1] / y[i]), float(y[i + 1] / y[i])]
            if significant and not (min(ratios) >= 1.1 and sig is not None and sig >= 3.0):
                continue
            out.append({"theta_cm_deg": float(th[i]), "data_mb_sr": float(y[i]), "neighbour_ratio": ratios,
                        "drop_in_sigma": None if sig is None else float(sig),
                        "implied_disk_radius_fm": J1_FIRST_ZERO / q})  # fmt: skip
    return out


def frozen_criterion_minima(
    th: NDArray[np.float64], y: NDArray[np.float64], rel: NDArray[np.float64], ru_share: Any
) -> list[dict[str, Any]]:
    """The frozen D-E1 diffraction-minimum criterion (plan section 5 text): the dip is present in at least 2 adjacent
    points (distinct angles, duplicates merged; both lie below the two flanking points), the lowest point is deeper
    than 2 combined sigma_stat below the lower flank, and the model Rutherford share at it is below 20 %."""
    t, v, r = merge_duplicates(th, y, rel)
    out = []
    for i in range(1, t.size - 2):
        if t[i] <= 10.0:
            continue
        lo_i = i if v[i] <= v[i + 1] else i + 1
        flank = min(v[i - 1], v[i + 2])
        if max(v[i], v[i + 1]) < flank:
            fl = i - 1 if v[i - 1] <= v[i + 2] else i + 2
            err = (
                math.hypot(r[lo_i] * v[lo_i], r[fl] * v[fl])
                if np.isfinite(r[lo_i]) and np.isfinite(r[fl])
                else float("nan")
            )
            sig = (flank - v[lo_i]) / err if np.isfinite(err) and err > 0 else None
            share = float(ru_share(t[lo_i]))
            out.append({"theta_cm_deg": float(t[lo_i]), "depth_sigma": sig, "rutherford_share": share,
                        "met": bool(sig is not None and sig > 2.0 and share < 0.2)})  # fmt: skip
    return out


def table_binned_first_minimum_deg(
    sh: Shapes, name: str, e_mev: float, bin_deg: float = 1.0, above: float = 3.0
) -> float | None:
    """The Amendment 19 section 2.2 estimator on the pinned table: first local minimum, above ``above`` deg, of the
    CM density binned in ``bin_deg`` bins (bin centre), from the table's own quantiles (2^16 points)."""
    u = (np.arange(1 << 16) + 0.5) / (1 << 16)
    th = np.degrees(np.arccos(np.clip(sh.mu_quantiles("table", name, e_mev, u), -1.0, 1.0)))
    edges = np.arange(0.0, 180.0 + bin_deg, bin_deg)
    cnt, _ = np.histogram(th, bins=edges)
    dens = cnt / (np.cos(np.radians(edges[:-1])) - np.cos(np.radians(edges[1:])))
    cen = 0.5 * (edges[1:] + edges[:-1])
    return first_local_min_deg(cen, dens.astype(float), above)


def table_bin_density_first_minimum_deg(
    sh: Shapes, name: str, e_mev: float, above: float = 10.0
) -> float | None:
    """First local minimum above ``above`` deg of the density of the table's equiprobable bins (1 / (n_q d mu) at the
    bin centre angle; the estimator behind the reviewer's 35.3 / 29.2 deg for p+O at 135 / 179.9 MeV)."""
    ed = sh.tab.edges_at(sh.index(name), e_mev)
    th = np.degrees(np.arccos(0.5 * (ed[1:] + ed[:-1])))
    dens = 1.0 / (ed.size - 1) / np.diff(ed)
    o = np.argsort(th)
    return first_local_min_deg(th[o], dens[o], above)


def run() -> None:
    v10 = c.import_transport_module("nuclear_v10")
    tab = c.load_table()
    sh = Shapes(tab)
    model = v10.ModelTable(c.TABLE_ID, c.CACHE_DIR)
    datasets: list[dict[str, Any]] = []
    for spec in v10.GATING_SETS:
        for rec in v10.read_angular(spec, c.CACHE_DIR):
            datasets.append(_one_dataset(v10, sh, model, rec))
    doc = {"banner": c.BANNER, "table_id": c.TABLE_ID, "datasets": datasets,
           "families": list(FAMILIES), "note": "R fit is data-seen and is never a model parameter"}  # fmt: skip
    inputs = [*c.table_files(), c.RUN_V10A]
    c.write_outputs(c.OUT / "de1", "per_dataset.json", doc, "de1-report.md", _report(doc, inputs))


def _one_dataset(v10: Any, sh: Shapes, model: Any, rec: dict[str, Any]) -> dict[str, Any]:
    th, y, rel_err = rec["theta_cm"], rec["ds_cm_mb_sr"], rec["rel_err"]
    name, e = rec["target"], float(rec["energy_mev"])
    r_tab, pcm = sh.r_pcm(name, e)
    k = pcm / HBARC
    tot, ru = model.total_mb_sr(name, e, th)
    fam = {kind: sh.density_mb_sr(kind, name, e, th) + ru for kind in FAMILIES}
    # measured minima against the model minima
    mins = measured_minima(th, y, rel_err, k, significant=False)
    sig_mins = measured_minima(th, y, rel_err, k, significant=True)
    fine_ru = np.linspace(5.0, 170.0, 3301)
    ru_f = model.total_mb_sr(name, e, fine_ru)
    frozen = frozen_criterion_minima(
        th, y, rel_err, lambda a: np.interp(a, fine_ru, ru_f[1] / ru_f[0])
    )
    tmin_data = v10.first_minimum_deg(th, y)
    fine = np.linspace(10.0, 170.0, 3201)
    nuc_fine = model.nuclear_mb_sr(name, e, fine)
    tot_fine, ru_fine = model.total_mb_sr(name, e, fine)
    tmin_nuc, tmin_tot = first_local_min_deg(fine, nuc_fine), first_local_min_deg(fine, tot_fine)
    th_zero = math.degrees(2.0 * math.asin(min(1.0, J1_FIRST_ZERO / (2.0 * k * r_tab))))
    d_min: dict[str, Any] = {"measured_minima": mins, "significant_measured_minima": sig_mins, "frozen_criterion_minima": frozen,
                             "frozen_criterion_met_at_v10a_minimum": any(f["met"] and tmin_data is not None and abs(f["theta_cm_deg"] - tmin_data) < 0.5 for f in frozen),
                             "v10a_cut_deg_08_theta_min": None if tmin_data is None else 0.8 * tmin_data,
                             "table_first_minimum_binned_1deg": table_binned_first_minimum_deg(sh, name, e),
                             "table_first_minimum_bin_density": table_bin_density_first_minimum_deg(sh, name, e), "first_measured_minimum_v10a_rule": tmin_data,
                             "table_first_minimum_nuclear_deg": tmin_nuc, "table_first_minimum_total_deg": tmin_tot,
                             "table_J1_first_zero_deg": th_zero, "table_radius_fm": r_tab}  # fmt: skip
    if tmin_data is not None:
        i = int(np.argmin(abs(fine - tmin_data)))
        d_min["rutherford_share_at_measured_minimum"] = float(ru_fine[i] / tot_fine[i])
        d_min["measured_to_table_minimum_ratio"] = tmin_data / tmin_tot if tmin_tot else None
        qs = [2.0 * k * math.sin(math.radians(m["theta_cm_deg"]) / 2.0) for m in mins]
        if len(qs) >= 2:
            d_min["q2_over_q1_measured_minima"] = qs[1] / qs[0]  # disk: 7.0156/3.8317 = 1.831
        d_min["data_first_minimum_is_significant"] = any(
            abs(m["theta_cm_deg"] - tmin_data) < 0.3 for m in sig_mins
        )
    # partial integrals from 10 deg (log-linear data value at 10 deg if bracketed) to a window end
    d10 = v10.loglin(th, y, 10.0)
    sig_el_mb = model.sigma_el_mb(name, e)
    part: dict[str, Any] = {}
    for wname, hi_req in (("full", float(th.max())), ("10-30", 30.0), ("10-40", 40.0)):
        hi = min(hi_req, float(th.max()))
        sel = (th > 10.0) & (th <= hi)
        nodes = np.concatenate(([10.0], th[sel])) if d10 is not None else th[sel]
        ydat = np.concatenate(([d10], y[sel])) if d10 is not None else y[sel]
        if nodes.size < 3:
            continue
        tot_n, ru_n = model.total_mb_sr(name, e, nodes)
        blk: dict[str, Any] = {"theta_range_deg": [float(nodes[0]), float(nodes[-1])], "n_nodes": int(nodes.size),
                               "window_truncated_by_data": hi < hi_req,
                               "data": partial_moments(nodes, ydat), "model": partial_moments(nodes, tot_n)}  # fmt: skip
        fine = np.linspace(nodes[0], nodes[-1], 3001)
        blk["model_fine_grid"] = partial_moments(fine, model.total_mb_sr(name, e, fine)[0])
        blk["rutherford_only"] = partial_moments(nodes, ru_n)
        for kind in FAMILIES:
            blk[kind] = partial_moments(nodes, sh.density_mb_sr(kind, name, e, nodes) + ru_n)
        for key in ("model", *FAMILIES):
            blk[key]["sigma_tr_ratio_to_data"] = (
                blk[key]["sigma_tr_mb"] / blk["data"]["sigma_tr_mb"]
            )
            blk[key]["sigma_ratio_to_data"] = blk[key]["sigma_mb"] / blk["data"]["sigma_mb"]
        blk["data_partial_sigma_over_table_sigma_el"] = blk["data"]["sigma_mb"] / sig_el_mb
        part[wname] = blk
    # data-seen one-parameter diagnostic: R of the black-disk shape (sigma_el, Rutherford fixed)
    fits = {}
    windows = {"10-30deg": (10.0, 30.0), "10-40deg": (10.0, 40.0),
               "to_V10A_first_minimum": (10.0, tmin_data if tmin_data is not None else 30.0)}  # fmt: skip
    for wname, (lo, hi) in windows.items():
        m = (th >= lo) & (th <= hi)
        if m.sum() < 3:
            fits[wname] = None
            continue
        rs = np.arange(1.0, 8.0001, 0.01)
        ru_m = ru[m]
        chi = np.array(
            [
                float(
                    np.sum(
                        (
                            np.log(
                                sh.density_mb_sr("disk", name, e, th[m], r_override=float(r)) + ru_m
                            )
                            - np.log(y[m])
                        )
                        ** 2
                    )
                )
                for r in rs
            ]
        )
        j = int(np.argmin(chi))
        fits[wname] = {"R_fit_fm": float(rs[j]), "R_table_fm": r_tab, "ratio": float(rs[j] / r_tab),
                       "rms_log_residual": math.sqrt(chi[j] / m.sum()), "n_points": int(m.sum()),
                       "rms_log_residual_at_table_R": math.sqrt(float(np.sum((np.log(tot[m]) - np.log(y[m])) ** 2)) / m.sum()),
                       "label": "data-seen diagnostic; never a model parameter"}  # fmt: skip
    return {"dataset": rec["dataset"], "target": name, "energy_mev": e, "n_points": rec["n_points"],
            "theta_range_deg": [float(th.min()), float(th.max())], "theta_cm_deg": th.tolist(), "data_mb_sr": y.tolist(), "rel_err": [None if not np.isfinite(x) else float(x) for x in rel_err],
            "model_total_mb_sr": tot.tolist(), "model_ratio_minus_1": (tot / y - 1.0).tolist(),
            "rutherford_mb_sr": ru.tolist(), "family_total_mb_sr": {kk: vv.tolist() for kk, vv in fam.items()},
            "minima": d_min, "partial": part, "R_fit": fits}  # fmt: skip


def _fmt(x: Any, f: str = ".3g") -> str:
    return "-" if x is None else format(x, f)


def _report(doc: dict[str, Any], inputs: list[Any]) -> str:
    ds = doc["datasets"]
    cmd = "python -m validation.scripts.elastic_attribution.de1_shapes"
    out = [
        c.header(
            "D-E1: pinned p+A shape and shape families against the V10-A data",
            cmd,
            inputs,
            [
                "Model = pinned table (analytic black disk of the table radius, plus unscreened Rutherford, incoherent; "
                "the V10-A model value). Families: Amendment 19 (b) at their a-priori parameters (no fit). "
                "Partial integrals: trapezoid over the data nodes from 10 deg to the window end (10-30, 10-40 deg, and the last data angle), nuclear + Rutherford. "
                "The R fit is **data-seen**, a diagnostic, never a model parameter."
            ],
        )
    ]
    out += ["", "## Minima: measured against the table", "",
            "| dataset | target | E | theta_min,1 (V10-A rule) | significant? | table first min: bin-density / 1-deg bins / analytic disk (total) / J1 zero | Ru share at data min | R_table | R implied by V10-A min (3.8317/q) | significant minima (deg, merged duplicates, >=10 % and >=3 sigma below both neighbours) |",
            "|---|---|---|---|---|---|---|---|---|---|"]  # fmt: skip
    for d in ds:
        m = d["minima"]
        sig_list = (
            ", ".join(f"{x['theta_cm_deg']:.1f}" for x in m["significant_measured_minima"])
            or "none"
        )
        rim = (
            ", ".join(f"{x['implied_disk_radius_fm']:.2f}" for x in m["measured_minima"][:1]) or "-"
        )
        out.append(
            f"| {d['dataset']} | {d['target']} | {d['energy_mev']:g} | {_fmt(m['first_measured_minimum_v10a_rule'], '.2f')} | "
            f"{m.get('data_first_minimum_is_significant', '-')} | "
            f"{_fmt(m['table_first_minimum_bin_density'], '.1f')} / {_fmt(m['table_first_minimum_binned_1deg'], '.1f')} / {_fmt(m['table_first_minimum_total_deg'], '.1f')} / {m['table_J1_first_zero_deg']:.1f} | "
            f"{_fmt(m.get('rutherford_share_at_measured_minimum'), '.3f')} | {m['table_radius_fm']:.2f} | {rim} | {sig_list} |"
        )
    out += ["", "**Frozen criterion** (plan: dip in >= 2 adjacent distinct-angle points, deeper than 2 sigma_stat, model Rutherford share < 20 %), "
            "applied to every dataset; the `significant` column above is my stricter (>= 10 % and >= 3 sigma) cross-check.", "",
            "| dataset | V10-A rule theta_min,1 | drop vs lower neighbour (sigma, duplicates not merged) | V10-A cut 0.8 theta_min | frozen criterion met at that minimum | frozen-criterion minima anywhere (deg; depth sigma; Ru share) |",
            "|---|---|---|---|---|---|"]  # fmt: skip
    for d in ds:
        m = d["minima"]
        f0 = m["measured_minima"][0] if m["measured_minima"] else None
        lst = (
            "; ".join(
                f"{x['theta_cm_deg']:.1f} ({_fmt(x['depth_sigma'], '.1f')}; {x['rutherford_share']:.2f}; {'met' if x['met'] else 'not met'})"
                for x in m["frozen_criterion_minima"]
            )
            or "none"
        )
        out.append(
            f"| {d['dataset']} | {_fmt(m['first_measured_minimum_v10a_rule'], '.2f')} | {_fmt(f0 and f0['drop_in_sigma'], '.2f')} | {_fmt(m['v10a_cut_deg_08_theta_min'], '.2f')} | "
            f"{m['frozen_criterion_met_at_v10a_minimum'] if f0 else '-'} | {lst} |"
        )
    out += [
        "",
        "**Finding.** The V10-A minima of C0063 @ 135 MeV (17.17 deg) and C0141 @ 179.9 MeV (14.99 deg) are duplicate-angle artefacts "
        "(a 0.74 sigma / 0.48 sigma wiggle between two runs at the same angle), not diffraction minima: the data fall monotonically over 10-38 deg. "
        "The recorded V10-A evaluation used them to set the (a) cut 0.8 theta_min,1 = 13.74 deg (C0063; 11.99 deg for C0141), which excluded the 15 and 20 deg (a) points of C0063 and C0141 "
        '(`03-v10-a.txt`: "theta >= 0.8 theta_min,1 (measured: 13.74 deg)"). This is a defect of the V10-A minimum rule (duplicate angles not merged), recorded here as a finding; '
        "the V10-A verdict is unchanged and nothing is re-evaluated. The statement 'first minima about 2x too wide / R about 4.9 fm' rested on these minima and is withdrawn.",
        "",
    ]
    for wname, title in (
        ("10-30", "10-30 deg"),
        ("10-40", "10-40 deg"),
        ("full", "10 deg to the last data angle"),
    ):
        out += ["", f"## Partial sigma_tr, window {title} (model / data; families at a-priori parameters)", "",
                "| dataset | range (deg) | data sigma_tr (mb) | model | " + " | ".join(FAMILIES) + " | model sigma/data sigma | data partial sigma / table sigma_el | Ru-only / data |",
                "|---|---|---|---|" + "---|" * len(FAMILIES) + "---|---|---|"]  # fmt: skip
        for d in ds:
            p = d["partial"].get(wname)
            if p is None:
                out.append(
                    f"| {d['dataset']} | not covered | - | - |" + " - |" * (len(FAMILIES) + 3)
                )
                continue
            flag = " (window cut by data)" if p["window_truncated_by_data"] else ""
            out.append(
                f"| {d['dataset']} | {p['theta_range_deg'][0]:.0f}-{p['theta_range_deg'][1]:.0f}{flag} | {p['data']['sigma_tr_mb']:.3g} | "
                f"{p['model']['sigma_tr_ratio_to_data']:.3f} | "
                + " | ".join(f"{p[k]['sigma_tr_ratio_to_data']:.3f}" for k in FAMILIES)
                + f" | {p['model']['sigma_ratio_to_data']:.3f} | {p['data_partial_sigma_over_table_sigma_el']:.2f} | {p['rutherford_only']['sigma_tr_mb'] / p['data']['sigma_tr_mb']:.3f} |"
            )
    out += [
        "",
        "The `full` window is dominated by the model's high-angle diffraction lobes (model 100-600x above the data at 60-90 deg for the O sets), "
        "so the 10-30 / 10-40 deg windows are the like-for-like comparison. Datasets with `data partial sigma / table sigma_el` above 1 "
        "(the data integral alone exceeds the table sigma_el) are not credible against the table: **O0553008 @ 75 MeV (p+C, Rolland 1966)** "
        "has a 10 deg data value about 8.4x the model; its EXFOR subentry gives MB/SR with a 5 % digital error and no other factor, "
        "so the cause is not an ingestion unit error that can be identified here (it is a source/transcription or frame issue of that subentry); "
        "it is reported but excluded from summaries and nothing was changed.",
        "",
    ]
    out += ["", "## Model / data - 1 at the V10-A angles (nuclear + Rutherford)", "",
            "| dataset | 10 deg | 15 deg | 20 deg | 30 deg |", "|---|---|---|---|---|"]  # fmt: skip
    for d in ds:
        th = np.array(d["theta_cm_deg"])
        r = np.array(d["model_ratio_minus_1"])
        cells = []
        for ang in (10.0, 15.0, 20.0, 30.0):
            j = int(np.argmin(abs(th - ang)))
            cells.append(f"{r[j]:+.2f} @{th[j]:.1f}" if abs(th[j] - ang) <= 1.5 else "-")
        out.append(f"| {d['dataset']} | " + " | ".join(cells) + " |")
    out += ["", "## One-parameter R diagnostic (DATA-SEEN; never a model parameter)", "",
            "Black-disk shape with sigma_el and the Rutherford term fixed; least squares in log(dsigma) on the stated CM window "
            "(R scanned 1-8 fm in 0.01 fm steps). The windows 10-30 and 10-40 deg are fixed a priori; `to_V10A_first_minimum` "
            "ends at the V10-A rule minimum and is only as meaningful as that minimum.", "",
            "| dataset | window | R_fit (fm) | R_table (fm) | R_fit/R_table | rms log resid at R_fit | at R_table |", "|---|---|---|---|---|---|---|"]  # fmt: skip
    for d in ds:
        for w, f in d["R_fit"].items():
            if f is None:
                out.append(f"| {d['dataset']} | {w} | - | - | - | - | - |")
            else:
                out.append(
                    f"| {d['dataset']} | {w} | {f['R_fit_fm']:.2f} | {f['R_table_fm']:.2f} | {f['ratio']:.2f} | "
                    f"{f['rms_log_residual']:.3f} | {f['rms_log_residual_at_table_R']:.3f} |"
                )
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    run()
