"""Deterministic nuclear-table checks N1, V1 and V1b of V3-005A (decision 0041).

N1   runtime Sigma_mass (table grid, linear in ln E) vs the independent numpy composition from the
     native ENDF TAB1 interpolation, at every ENDF node and node midpoint, 1-150 MeV, water and the
     ICRU/ICRP tissues of ``ionmc.materials``: rel <= 1e-3.
V1   per-element sigma (table) vs MF3/MT5 at every ENDF node: rel <= 1e-3; extension continuity
     sigma(150+)/sigma(150-) - 1 within 1e-6.
V1b  informative: weighted mean ratio table/EXFOR per window with the PDG scale factor
     (D0356 gating-class evidence; C1862 and geant-val report-only).

Usage: ``python nuclear_checks.py --table-id ID [--cache-dir DIR] [--out summary.json]``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from ionmc import materials as M
from ionmc.data import cache, endf6, exfor, geant_val
from ionmc.nuclear import build as B
from ionmc.nuclear.tables import NuclearTable

WINDOWS = ((20.0, 40.0), (40.0, 70.0), (70.0, 110.0), (110.0, 160.0), (160.0, 250.0))
MATERIALS = (
    M.WATER,
    M.ADIPOSE_TISSUE_ICRP,
    M.MUSCLE_SKELETAL_ICRP,
    M.BONE_COMPACT_ICRU,
    M.LUNG_ICRP,
)


def _targets(cdir: Path) -> dict[str, B.TargetTables]:
    from ionmc.data.ame import load_ame2020

    zp = cache.verify("endf-b8.0-protons", cdir)
    ame_tab = load_ame2020(cache.verify("ame2020-mass", cdir).read_text(encoding="ascii"))
    out = {}
    for spec in B.TARGETS:
        mat = endf6.parse_endf(endf6.read_member(zp, f"ENDF-B-VIII.0_protons/{spec.member}.endf"))
        out[spec.name] = B.TargetTables(spec, mat, B.ev.build_event_model(ame_tab, spec.z, spec.a))
    return out


def weighted_ratio(ratio: np.ndarray, sigma: np.ndarray) -> tuple[float, float, float] | None:
    """Weighted mean, its error and the PDG scale factor of ``ratio`` with errors ``sigma``."""
    if ratio.size == 0:
        return None
    w = 1.0 / sigma**2
    mean = float(np.sum(w * ratio) / np.sum(w))
    err = float(1.0 / math.sqrt(np.sum(w)))
    chi2 = float(np.sum(w * (ratio - mean) ** 2))
    scale = max(1.0, math.sqrt(chi2 / (ratio.size - 1))) if ratio.size > 1 else 1.0
    return mean, err * scale, scale


def run_checks(cache_dir: str | Path | None, table_id: str) -> dict[str, Any]:
    cdir = cache.resolve_cache_dir(cache_dir)
    tab = NuclearTable.load(cdir, table_id)
    tts = _targets(cdir)
    grid = tab.arrays["grid_e_mev"]
    ln_grid = np.log(grid)
    elements = tab.info["elements"]
    # ---- V1: per element at every ENDF node ------------------------------------------------
    v1: dict[str, Any] = {}
    all_nodes: list[np.ndarray] = []
    for i, name in enumerate(tab.target_names):
        x = tts[name].sigma_tab.x * 1e-6
        y = tts[name].sigma_tab.y
        sel = (x >= 1.0) & (x <= 150.0)
        x, y = x[sel], y[sel]
        all_nodes.append(x)
        run = np.interp(np.log(x), ln_grid, tab.arrays["sigma_barn"][i])
        pos = y > 0.0
        rel = np.abs(run[pos] / y[pos] - 1.0)
        strong = y[pos] >= 0.01 * y.max()
        v1[name] = {
            "max_rel_all_positive_nodes": float(rel.max(initial=0.0)),
            "max_rel": float(rel[strong].max(initial=0.0)),
            "max_rel_above_10_mev": float(rel[x[pos] >= 10.0].max(initial=0.0)),
            "at_mev": float(x[pos][strong][int(np.argmax(rel[strong]))]) if strong.any() else None,
            "max_abs_barn_zero_nodes": float(np.abs(run[~pos]).max(initial=0.0)),
        }
    scale = {}
    for sym, row in elements.items():
        i = tab.target_names.index(row["target"])
        s150 = float(np.interp(math.log(150.0), ln_grid, tab.arrays["sigma_barn"][i]))
        scale[sym] = row["sigma_scale"] * s150
    cont = {}
    for name, tt in tts.items():
        lo = float(tt.sigma_barn(np.array([150.0]))[0])
        hi = float(tt.sigma_barn(np.array([150.0 * (1.0 + 1e-9)]))[0])
        cont[name] = hi / lo - 1.0
    # ---- N1: compositions at nodes and midpoints -------------------------------------------
    nodes = np.unique(np.concatenate(all_nodes))
    pts = np.unique(np.concatenate((nodes, 0.5 * (nodes[:-1] + nodes[1:]))))
    pts = pts[(pts >= 1.0) & (pts <= 150.0)]
    sig_native = {n: tts[n].sigma_native_barn(pts) for n in tts}
    n1: dict[str, Any] = {}
    for mat in MATERIALS:
        rows = tab.material_rows(mat)
        run = np.interp(np.log(pts), ln_grid, rows.sigma_mass_cm2_g)
        ref = B.material_sigma_mass(mat, sig_native, elements)
        pos = ref > 0.0
        rel = np.abs(run[pos] / ref[pos] - 1.0)
        strong = ref[pos] >= 0.01 * ref.max()
        n1[mat.name] = {
            "max_rel_all_positive_points": float(rel.max(initial=0.0)),
            "max_rel": float(rel[strong].max(initial=0.0)),
            "max_rel_above_10_mev": float(rel[pts[pos] >= 10.0].max(initial=0.0)),
            "at_mev": float(pts[pos][strong][int(np.argmax(rel[strong]))])
            if strong.any()
            else None,
            "n_points": int(pos.sum()),
            "majorant_ge_sigma": bool(np.all(rows.sigma_hat_window >= rows.sigma_mass_cm2_g)),
        }
    # ---- V1b -------------------------------------------------------------------------------
    v1b = _v1b(cdir, tab, tts)
    summary = {
        "table_id": table_id,
        "N1": {
            "tolerance": 1e-3,
            "materials": n1,
            "max_rel": max(v["max_rel"] for v in n1.values()),
        },
        "V1": {
            "tolerance": 1e-3,
            "targets": v1,
            "max_rel": max(v["max_rel"] for v in v1.values()),
            "extension_continuity": cont,
            "extension_continuity_max_abs": max(abs(v) for v in cont.values()),
        },
        "V1b": v1b,
    }
    summary["N1"]["pass"] = summary["N1"]["max_rel"] <= 1e-3
    summary["V1"]["pass"] = (
        summary["V1"]["max_rel"] <= 1e-3 and summary["V1"]["extension_continuity_max_abs"] <= 1e-6
    )
    return summary


def _v1b(cdir: Path, tab: NuclearTable, tts: dict[str, B.TargetTables]) -> dict[str, Any]:
    ln_grid = np.log(tab.arrays["grid_e_mev"])

    def table_mb(name: str, e: np.ndarray) -> np.ndarray:
        i = tab.target_names.index(name)
        return 1e3 * np.interp(np.log(e), ln_grid, tab.arrays["sigma_barn"][i])

    out: dict[str, Any] = {"windows": [list(w) for w in WINDOWS], "sets": {}}
    specs = {"exfor-d0356": "evaluation (post-1997)", "exfor-c1862": "report-only (pre-1997)"}
    for did, role in specs.items():
        entry = exfor.parse_entry(cache.verify(did, cdir).read_text(encoding="ascii"))
        res: dict[str, Any] = {"role": role}
        pts: dict[str, list[tuple[float, float, float]]] = {}
        for sub in entry.subentries:
            if sub.data is None or not sub.reactions:
                continue
            txt = sub.reactions[0].text
            tgt = "C-12" if "6-C-12(P,NON)" in txt else "O-16" if "8-O-16(P,NON)" in txt else None
            tgt = "Ca-40" if "20-CA-40(P,NON)" in txt else tgt
            h = sub.data.heads
            if tgt is None or "EN" not in h or "DATA" not in h:
                continue
            en, xs = sub.data.column(h.index("EN")), sub.data.column(h.index("DATA"))
            er_i = [i for i, k in enumerate(h) if k.startswith("ERR")]
            er = sub.data.column(er_i[0]) if er_i else 0.05 * xs
            er = np.where(np.isfinite(er) & (er > 0), er, 0.05 * xs)
            en = exfor.energy_total_mev(en, sub.data.units[h.index("EN")], 1)
            xs = exfor.xs_to_mb(xs, sub.data.units[h.index("DATA")])
            er = exfor.xs_to_mb(er, sub.data.units[er_i[0]]) if er_i else er
            pts.setdefault(tgt, []).extend(zip(en, xs, er, strict=True))
        for tgt, rows in pts.items():
            a = np.array(rows)
            a = a[np.all(np.isfinite(a), axis=1) & (a[:, 1] > 0)]
            if a.size == 0:
                continue
            ratio = table_mb(tgt, a[:, 0]) / a[:, 1]
            sig = ratio * a[:, 2] / a[:, 1]
            wins = {}
            for lo, hi in WINDOWS:
                m = (a[:, 0] >= lo) & (a[:, 0] < hi if hi < 250 else a[:, 0] <= hi)
                r = weighted_ratio(ratio[m], sig[m]) if m.any() else None
                wins[f"{lo:g}-{hi:g}"] = (
                    None
                    if r is None
                    else {"ratio": r[0], "err": r[1], "pdg_scale": r[2], "n": int(m.sum())}
                )
            res[tgt] = wins
        out["sets"][did] = res
    gv_path = cache.verify("geant-val-exfor-inelastic-7", cdir)
    gv = {}
    for c in geant_val.parse_geant_val(gv_path.read_text(encoding="utf-8")):
        if c.beam.lower() != "proton":
            continue
        name = {"C": "C-12", "O": "O-16"}.get(c.target.replace("12", "").replace("16", ""))
        if name is None or c.x.size == 0:
            continue
        y_mb = c.y * (1.0 if "mb" in c.y_axis.lower() else 1.0e3)
        ratio = table_mb(name, np.clip(c.x, 1.0, 250.0)) / y_mb
        gv[f"{c.target}#{c.record_id}"] = {
            "energies_mev": c.x.tolist(),
            "ratio": ratio.tolist(),
        }
    out["sets"]["geant-val-exfor-inelastic-7"] = {"role": "report-only", "curves": gv}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table-id", required=True)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    summary = run_checks(args.cache_dir, args.table_id)
    text = json.dumps(summary, indent=2, sort_keys=True)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if summary["N1"]["pass"] and summary["V1"]["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
