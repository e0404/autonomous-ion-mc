"""D9 (V3-005B, C17): offline, report-only evaluation of alternative non-elastic cross-section
models against the post-1997 EXFOR set (decision 0041 limitation 1; plan Amendment 7 (i)).

Models: ``la150`` (the rows of the nuclear table, i.e. what the runtime uses), ``tripathi``
(Tripathi light-system, ``ionmc.physics.tripathi``) and ``bgg`` (Geant4 BGG/Barashenkov), which is
**not evaluable offline**: its parameterisation (Barashenkov tables / Glauber-Gribov constants) is
not in the repository or the cache, and none is transcribed here.

Data roles (kept separate): ``exfor-d0356`` (Auce et al. 2005, post-1997) is the evaluation role
and the only set the decision rule uses; ``exfor-c1862`` (pre-1997) and ``geant-val`` (EXFOR-derived
curves) are report-only. Per element and window [20,40), [40,70), [70,110), [110,160),
[160,250] MeV: weighted mean ratio model/data, its error with the PDG scale factor, ``chi2/ndf``
of the scatter about the weighted mean (ndf = n - 1) and ``chi2`` of the ratios against unity
(ndf = n). Pre-registered rule (plan section E): an alternative is *preferred* for an element iff
``|1 - r_LA150| - |1 - r_alt| > sigma_comb`` in every window with evaluation-role data, with
``sigma_comb = sqrt(err_LA150^2 + err_alt^2)`` (quadrature of the two window errors; the plan does
not define the correlation, this is the literal reading). No model is adopted here.

Usage: ``python nuclear_d9.py --table-id ID [--cache-dir DIR] [--out-json F] [--out-md F]``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from nuclear_checks import WINDOWS, weighted_ratio  # noqa: E402

TARGETS = {"C-12": (6, 12), "N-14": (7, 14), "O-16": (8, 16), "Ca-40": (20, 40)}
_EXFOR_TARGET = {"6-C-12": "C-12", "7-N-14": "N-14", "8-O-16": "O-16", "20-CA-40": "Ca-40"}
_GV_TARGET = {"C": "C-12", "N": "N-14", "O": "O-16", "Ca": "Ca-40"}
SETS = {
    "exfor-d0356": "evaluation (post-1997)",
    "exfor-c1862": "report-only (pre-1997)",
    "geant-val-exfor-inelastic-7": "report-only (EXFOR-derived curves)",
}
EVALUATION_SET = "exfor-d0356"
ALTERNATIVES = ("tripathi", "bgg")
BGG_NOT_EVALUABLE = {
    "evaluable": False,
    "reason": "Geant4 BGG/Barashenkov parameterisation not present in the repository or cache",
    "acquisition_needed": {
        "what": "G4BGGNucleonInelasticXS / G4ComponentBarNucleonNucleusXsc data tables "
        "(Barashenkov p+A sigma_inel and the Glauber-Gribov constants)",
        "url": "https://github.com/Geant4/geant4 (source/processes/hadronic/cross_sections)",
        "licence": "Geant4 Software License; shares lineage with TOPAS QGSP_BIC_HP, so "
        "report-only role (adopting it would remove the V5 sigma independence)",
        "role": "report-only alternative sigma, offline comparison",
    },
}


def window_key(lo: float, hi: float) -> str:
    return f"{lo:g}-{hi:g}"


def in_window(e: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Half-open window [lo, hi); the last window (hi = 250) is closed."""
    return (e >= lo) & ((e < hi) if hi < 250.0 else (e <= hi))


def window_stats(ratio: np.ndarray, sigma: np.ndarray) -> dict[str, Any] | None:
    """Weighted mean ratio, error (PDG scale), chi2/ndf of the scatter and chi2 against unity."""
    r = weighted_ratio(ratio, sigma)
    if r is None:
        return None
    mean, err, scale = r
    w = 1.0 / sigma**2
    chi2 = float(np.sum(w * (ratio - mean) ** 2))
    ndf = int(ratio.size - 1)
    chi2_unity = float(np.sum(w * (ratio - 1.0) ** 2))
    return {
        "ratio": mean,
        "err": err,
        "pdg_scale": scale,
        "n": int(ratio.size),
        "chi2_scatter": chi2,
        "ndf": ndf,
        "chi2_per_ndf": chi2 / ndf if ndf > 0 else None,
        "chi2_vs_unity": chi2_unity,
        "chi2_vs_unity_per_n": chi2_unity / ratio.size,
    }


def evaluate_model(
    points: np.ndarray, model_mb: Callable[[np.ndarray], np.ndarray]
) -> dict[str, Any]:
    """Per-window statistics of ``model/data`` for rows ``(E, sigma_mb, err_mb)``."""
    out: dict[str, Any] = {}
    for lo, hi in WINDOWS:
        m = in_window(points[:, 0], lo, hi)
        if not m.any():
            out[window_key(lo, hi)] = None
            continue
        e, xs, er = points[m, 0], points[m, 1], points[m, 2]
        ratio = model_mb(e) / xs
        out[window_key(lo, hi)] = window_stats(ratio, ratio * er / xs)
    return out


def preferred(la: dict[str, Any], alt: dict[str, Any]) -> dict[str, Any]:
    """Pre-registered rule: alt beats LA150 by more than sigma_comb in every window with data."""
    rows, ok = {}, True
    for k, a in la.items():
        b = alt.get(k)
        if a is None or b is None:
            continue
        gain = abs(1.0 - a["ratio"]) - abs(1.0 - b["ratio"])
        comb = math.hypot(a["err"], b["err"])
        rows[k] = {"gain": gain, "sigma_comb": comb, "pass": bool(gain > comb)}
        ok &= gain > comb
    return {"windows": rows, "preferred": bool(rows) and ok, "n_windows": len(rows)}


def _exfor_points(text: str) -> dict[str, np.ndarray]:
    from ionmc.data import exfor

    entry = exfor.parse_entry(text)
    pts: dict[str, list[tuple[float, float, float]]] = {}
    for sub in entry.subentries:
        if sub.data is None or not sub.reactions:
            continue
        txt = sub.reactions[0].text
        tgt = next((v for k, v in _EXFOR_TARGET.items() if f"{k}(P,NON)" in txt), None)
        h = sub.data.heads
        if tgt is None or "EN" not in h or "DATA" not in h:
            continue
        en, xs = sub.data.column(h.index("EN")), sub.data.column(h.index("DATA"))
        er_i = [i for i, k in enumerate(h) if k.startswith("ERR") or k == "DATA-ERR"]
        er = sub.data.column(er_i[0]) if er_i else 0.05 * xs
        er = np.where(np.isfinite(er) & (er > 0), er, 0.05 * xs)
        en = exfor.energy_total_mev(en, sub.data.units[h.index("EN")], 1)
        xs = exfor.xs_to_mb(xs, sub.data.units[h.index("DATA")])
        er = exfor.xs_to_mb(er, sub.data.units[er_i[0]]) if er_i else er
        pts.setdefault(tgt, []).extend(zip(en, xs, er, strict=True))
    return {t: _clean(np.array(r)) for t, r in pts.items()}


def _clean(a: np.ndarray) -> np.ndarray:
    return a[np.all(np.isfinite(a), axis=1) & (a[:, 1] > 0)]


def _geant_val_points(text: str) -> dict[str, np.ndarray]:
    from ionmc.data import geant_val

    pts: dict[str, list[np.ndarray]] = {}
    for c in geant_val.parse_geant_val(text):
        name = _GV_TARGET.get(c.target.replace("12", "").replace("16", ""))
        if c.beam.lower() != "proton" or name is None or c.x.size == 0:
            continue
        f = 1.0 if "mb" in c.y_axis.lower() else 1.0e3
        stat = np.nan_to_num(np.asarray(c.y_stat_plus, dtype=float))
        sys_ = np.nan_to_num(np.asarray(c.y_sys_plus, dtype=float))
        err = np.hypot(stat, sys_)
        err = np.where(err > 0, err, 0.05 * c.y) * f
        pts.setdefault(name, []).append(np.column_stack((c.x, c.y * f, err)))
    return {t: _clean(np.vstack(v)) for t, v in pts.items()}


def load_data(cdir: Path) -> dict[str, dict[str, np.ndarray]]:
    """``{set id: {target: rows (E MeV, sigma mb, err mb)}}`` for the three data sets."""
    from ionmc.data import cache

    out = {d: _exfor_points(cache.verify(d, cdir).read_text(encoding="ascii"))
           for d in ("exfor-d0356", "exfor-c1862")}  # fmt: skip
    gv = "geant-val-exfor-inelastic-7"
    out[gv] = _geant_val_points(cache.verify(gv, cdir).read_text(encoding="utf-8"))
    return out


def run_d9(
    data: dict[str, dict[str, np.ndarray]],
    models: dict[str, Callable[[str], Callable[[np.ndarray], np.ndarray]]],
) -> dict[str, Any]:
    """Evaluate every model on every set and target; ``models[name](target)`` gives mb(E)."""
    out: dict[str, Any] = {
        "windows": [list(w) for w in WINDOWS],
        "roles": SETS,
        "evaluation_set": EVALUATION_SET,
        "bgg": BGG_NOT_EVALUABLE,
        "sets": {},
        "decision": {},
    }
    for sid, role in SETS.items():
        res: dict[str, Any] = {"role": role, "targets": {}}
        for tgt in TARGETS:
            pts = data.get(sid, {}).get(tgt)
            if pts is None or pts.size == 0:
                res["targets"][tgt] = {"data": "none"}
                continue
            res["targets"][tgt] = {
                "data": "present",
                "n_points": int(pts.shape[0]),
                "e_range_mev": [float(pts[:, 0].min()), float(pts[:, 0].max())],
                "models": {name: evaluate_model(pts, f(tgt)) for name, f in models.items()},
            }
        out["sets"][sid] = res
    for sid in SETS:
        dec: dict[str, Any] = {}
        for tgt, rec in out["sets"][sid]["targets"].items():
            if rec["data"] != "present":
                dec[tgt] = {"status": "no data"}
                continue
            dec[tgt] = {
                alt: preferred(rec["models"]["la150"], rec["models"][alt])
                for alt in ALTERNATIVES
                if alt in rec["models"]
            }
        out["decision"][sid] = dec
    return out


def _fmt(s: dict[str, Any] | None) -> str:
    if s is None:
        return "gap"
    c = "-" if s["chi2_per_ndf"] is None else f"{s['chi2_per_ndf']:.2f}"
    u = s["chi2_vs_unity_per_n"]
    return f"{s['ratio']:.3f} ± {s['err']:.3f} (n={s['n']}, χ²/ndf {c}, χ²₁/n {u:.2f})"


def markdown(res: dict[str, Any]) -> str:
    keys = [window_key(lo, hi) for lo, hi in WINDOWS]
    lines = ["# D9 report: alternative sigma_nonel vs EXFOR (report-only, no model adopted)", ""]
    for sid, s in res["sets"].items():
        lines += [f"## {sid} - {s['role']}", ""]
        for tgt, rec in s["targets"].items():
            if rec["data"] != "present":
                lines += [f"p+{tgt}: no data in this set.", ""]
                continue
            lo, hi = rec["e_range_mev"]
            lines += [f"### p+{tgt} ({rec['n_points']} points, {lo:.1f}-{hi:.1f} MeV)", ""]
            lines += ["| model | " + " | ".join(keys) + " |", "|---" * (len(keys) + 1) + "|"]
            for name, wins in rec["models"].items():
                lines.append(f"| {name} | " + " | ".join(_fmt(wins[k]) for k in keys) + " |")
            lines.append("| bgg | " + " | ".join("not evaluable" for _ in keys) + " |")
            lines.append("")
            for alt, d in res["decision"][sid][tgt].items():
                lines.append(f"- rule, {alt} vs la150: preferred = {d['preferred']} "
                             f"({d['n_windows']} windows with data)")  # fmt: skip
            lines.append("")
    b = res["bgg"]["acquisition_needed"]
    lines += [
        "## BGG",
        "",
        f"Not evaluable offline: {res['bgg']['reason']}.",
        f"Needed: {b['what']}; URL {b['url']}; licence {b['licence']}; role {b['role']}.",
    ]
    return "\n".join(lines) + "\n"


def run(cache_dir: str | Path | None, table_id: str) -> dict[str, Any]:
    from ionmc.data import cache
    from ionmc.nuclear.tables import NuclearTable
    from ionmc.physics.tripathi import sigma_tripathi_light_p

    cdir = cache.resolve_cache_dir(cache_dir)
    tab = NuclearTable.load(cdir, table_id)
    grid = tab.arrays["grid_e_mev"]

    def la150(t: str) -> Callable[[np.ndarray], np.ndarray]:
        row = tab.arrays["sigma_barn"][tab.target_names.index(t)]
        return lambda e: 1e3 * np.interp(e, grid, row)

    def tripathi(t: str) -> Callable[[np.ndarray], np.ndarray]:
        z, a = TARGETS[t]
        return lambda e: 1e3 * sigma_tripathi_light_p(e, z, a)

    res = run_d9(load_data(cdir), {"la150": la150, "tripathi": tripathi})
    res["table_id"] = table_id
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table-id", required=True)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--out-json", default=None)
    ap.add_argument("--out-md", default=None)
    args = ap.parse_args()
    res = run(args.cache_dir, args.table_id)
    if args.out_json:
        Path(args.out_json).write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    md = markdown(res)
    if args.out_md:
        Path(args.out_md).write_text(md, encoding="utf-8")
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
