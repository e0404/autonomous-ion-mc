"""Deterministic nuclear-table checks N1, V1, V1b, V4, V4b and D6 of V3-005A (decision 0041).

N1   runtime Sigma_mass (union grid, lin-lin in E) vs the independent numpy composition from the
     native ENDF TAB1 interpolation, at every ENDF node and node midpoint, 1-150 MeV, water and
     the ICRU/ICRP tissues of ``ionmc.materials``: rel <= 1e-3.
V1   per-element sigma (table) vs MF3/MT5 at every ENDF node: rel <= 1e-3; extension continuity
     sigma(150+)/sigma(150-) - 1 within 1e-6.
V1b  informative: weighted mean ratio table/EXFOR per window with the PDG scale factor
     (D0356 gating-class evidence; C1862 and geant-val report-only).
V4   event sampler vs ENDF, C-12 and O-16 at 100 and 150 MeV, 1e5 events: yields within 1 % + 3
     sem, sum y <E'> + recoil within 5 %.
V4b  regression guard: <E'> n, p, d within 1 %, alpha within 6 %, |mean Delta_lab| <= 0.12
     E_avail, exact P_accept >= 0.99, ledger closure <= 1e-9 MeV on every event.
D6   alpha local-deposition numbers of the table JSON with the tier and B1-ceiling booleans, and
     the capacity bound ``transport_energy_bound_mev`` with B_L in water (rho = 1 g/cm3).

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


def light_masses(ame: dict[tuple[int, int], Any]) -> list[float]:
    """Nuclear masses [MeV] of n, p, d, alpha and gamma (0) read directly from the AME2020 table."""
    from ionmc.data.ame import nuclear_mass_mev

    return [nuclear_mass_mev(ame, z, a) for z, a in ((0, 1), (1, 1), (1, 2), (2, 4))] + [0.0]


def residual_mass(ame: dict[tuple[int, int], Any], z_r: int, a_r: int) -> float:
    """AME2020 mass of the residual; the total break-up (0, 0) is the empty residual of mass 0."""
    from ionmc.data.ame import nuclear_mass_mev

    return 0.0 if (z_r, a_r) == (0, 0) else nuclear_mass_mev(ame, int(z_r), int(a_r))


def recompute_binding_events(
    counts: Any, z_r: Any, a_r: Any, z_t: int, a_t: int, ame: dict[tuple[int, int], Any]
) -> Any:
    """Binding ``sum m_out + M_r - m_p - M_t`` of every event from its light-product counts (n, p,
    d, alpha, gamma) and residual (Z_r, A_r), with all masses read from ``ame`` (the sampler's own
    ``binding_mev`` is not used)."""
    from ionmc.data.ame import nuclear_mass_mev

    m = light_masses(ame)
    m_p, m_t = m[1], nuclear_mass_mev(ame, z_t, a_t)
    m_r = np.array([residual_mass(ame, int(z), int(a)) for z, a in zip(z_r, a_r, strict=True)])
    return np.asarray(counts, dtype=np.float64) @ np.array(m) + m_r - m_p - m_t


def recompute_binding_total(
    nuc_diag: dict[str, Any], targets: list[dict[str, Any]], ame: dict[tuple[int, int], Any]
) -> float:
    """The run's summed binding [MeV] recomputed from the nuclear diagnostics block of the Result
    (per-target event, light-product and residual counts) and AME2020 masses."""
    from ionmc.data.ame import nuclear_mass_mev

    m = light_masses(ame)
    total = 0.0
    for tgt, rec in nuc_diag.items():
        info = targets[int(tgt)]
        m_t = nuclear_mass_mev(ame, int(info["z"]), int(info["a"]))
        total += sum(rec["light"][k] * m[i] for i, k in enumerate(("n", "p", "d", "a", "g")))
        for key, cnt in rec["residual"].items():
            z_r, a_r = (int(x) for x in key.split(","))
            total += cnt * residual_mass(ame, z_r, a_r)
        total -= rec["events"] * (m[1] + m_t)
    return float(total)


def _targets(cdir: Path) -> dict[str, B.TargetTables]:
    from ionmc.data.ame import load_ame2020

    zp = cache.verify("endf-b8.0-protons", cdir)
    ame_tab = load_ame2020(cache.verify("ame2020-mass", cdir).read_text(encoding="ascii"))
    out = {}
    for spec in B.TARGETS:
        mat = endf6.parse_endf(endf6.read_member(zp, f"ENDF-B-VIII.0_protons/{spec.member}.endf"))
        model = B.ev.build_event_model(ame_tab, spec.z, spec.a)
        out[spec.name] = B.TargetTables(spec, mat, model)
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
        run = np.interp(x, grid, tab.arrays["sigma_barn"][i])
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
        s150 = float(np.interp(150.0, grid, tab.arrays["sigma_barn"][i]))
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
        run = np.interp(pts, grid, rows.sigma_mass_cm2_g)
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
    v4 = v4_checks(tab, tts)
    summary = {
        "table_id": table_id,
        "N1": {
            "tolerance": 1e-3,
            "grid_points": int(grid.size),
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
        "V4": v4,
        "D6": d6_report(tab),
    }
    summary["N1"]["pass"] = summary["N1"]["max_rel"] <= 1e-3
    summary["V1"]["pass"] = (
        summary["V1"]["max_rel"] <= 1e-3 and summary["V1"]["extension_continuity_max_abs"] <= 1e-6
    )
    return summary


V4_CASES = (("C-12", 100.0), ("C-12", 150.0), ("O-16", 100.0), ("O-16", 150.0))
V4_EVENTS = 100_000
V4_BASE_SEED = 20421004 + 1000 * 9  # plan r_index 9 (V4), shard 0


def v4_checks(tab: NuclearTable, tts: dict[str, B.TargetTables]) -> dict[str, Any]:
    """V4 and V4b on the loaded table (module docstring)."""
    grid = tab.arrays["grid_e_mev"]
    out: dict[str, Any] = {"events": V4_EVENTS, "cases": {}}
    ok_v4 = ok_v4b = True
    for ic, (name, e) in enumerate(V4_CASES):
        it = tab.target_names.index(name)
        model = tts[name].model
        rows = B.ev.interp_rows(
            grid, tab.arrays["lam"][it], tab.arrays["edges_mev"][it], tab.arrays["r_pre"][it],
            tab.arrays["recoil_t_cm_mev"][it], e,
        )  # fmt: skip
        endf = tts[name].rows_at(e)
        src = B.ev.CounterUniforms(V4_BASE_SEED + ic)
        diag = B.diagnostics_block(model, rows, endf, e, V4_EVENTS, V4_BASE_SEED + ic)
        # ledger closure on every event
        batch = B.ev.sample_events(model, rows, e, V4_EVENTS, src, want_particles=True)
        assert batch.particle_event is not None and batch.particle_lab is not None
        assert batch.particle_species is not None
        t_lab = batch.particle_lab[:, 0] - model.species_mass_mev[batch.particle_species]
        sum_t = np.bincount(batch.particle_event, weights=t_lab, minlength=V4_EVENTS)
        acc = batch.accepted
        clos = np.abs(
            e - sum_t[acc] - batch.recoil_t_mev[acc] - batch.binding_mev[acc]
            - batch.imbalance_mev[acc]
        )  # fmt: skip
        ratio = np.array(diag["yield_ratio"])
        tol = 0.01 + np.array(diag["yield_ratio_3sem"])
        v4_yield = bool(np.all(np.abs(ratio - 1.0) <= tol))
        v4_sum = abs(diag["sum_y_e_prime_ratio"] - 1.0) <= 0.05
        er = np.array(diag["mean_e_prime_ratio"])
        sem_e = np.full(5, np.nan)
        for sp in range(5):
            sel_sp = batch.particle_species == sp
            if sel_sp.sum() > 1 and endf["mean_ecm"][sp] > 0:
                sem_e[sp] = (
                    batch.e_cm_mev[sel_sp].std() / math.sqrt(sel_sp.sum()) / endf["mean_ecm"][sp]
                )
        lim = np.array([0.01, 0.01, 0.01, 0.06, np.inf]) + 3.0 * np.nan_to_num(sem_e)
        e_pass = np.abs(er - 1.0) <= lim
        v4b_e = bool(np.all(e_pass[:4]))
        v4b_d = abs(diag["delta_lab_mean_mev"]) <= 0.12 * diag["e_avail_mev"]
        v4b_p = diag["p_accept_exact"] >= 0.99
        v4b_c = float(clos.max()) <= 1e-9
        diag.update(
            {
                "delta_lab_over_e_avail": diag["delta_lab_mean_mev"] / diag["e_avail_mev"],
                "max_ledger_closure_mev": float(clos.max()),
                "mean_e_prime_ratio_sem": sem_e.tolist(),
                "mean_e_prime_ratio_pass": e_pass.tolist(),
                "v4_yields_pass": v4_yield,
                "v4_sum_pass": bool(v4_sum),
                "v4b_pass": bool(v4b_e and v4b_d and v4b_p and v4b_c),
            }
        )
        ok_v4 &= v4_yield and bool(v4_sum)
        ok_v4b &= bool(v4b_e and v4b_d and v4b_p and v4b_c)
        out["cases"][f"{name}@{e:g}"] = diag
    out["V4_pass"] = bool(ok_v4)
    out["V4b_pass"] = bool(ok_v4b)
    return out


def d6_report(tab: NuclearTable) -> dict[str, Any]:
    """D6 numbers, tier and ceiling booleans of the table JSON, the capacity bound and B_L."""
    from ionmc.physics.projectiles import PROTON
    from ionmc.physics.stopping import BetheStoppingSource

    g = tab.info["gate_d6"]
    terms = tab.info["transport_path_bound_terms"]
    stop = BetheStoppingSource().table(M.WATER, PROTON)

    def path_mm(e_hi: float) -> float:
        """CSDA path [mm] in water (rho = 1 g/cm3) from the proton stopping table."""
        return float(10.0 * stop.range_at(e_hi))

    # the deuteron term uses the proton range (no deuteron table yet)
    b_l = 1.25 * (path_mm(250.0) + sum(t["n_max"] * path_mm(t["t_lab_max_mev"]) for t in terms.values()))
    return {
        "numbers": {k: dict(v) for k, v in g["numbers"].items()},
        "tier1_pass": g["tier1_pass"],
        "tier2_pass": g["tier2_pass"],
        "ceiling_pass": g["ceiling_pass"],
        "ceiling_pass_energy_weighted_range": g["ceiling_pass_energy_weighted_range"],
        "transport_energy_bound_mev": float(tab.info["transport_energy_bound_mev"]),
        "transport_path_bound_terms": {k: dict(v) for k, v in terms.items()},
        "path_water_250_mev_mm": path_mm(250.0),
        "path_terms_water_mm": {k: t["n_max"] * path_mm(t["t_lab_max_mev"]) for k, t in terms.items()},
        "b_l_water_rho1_mm": b_l,
        "b_l_note": "deuteron term evaluated with the proton range",
    }


def _v1b(cdir: Path, tab: NuclearTable, tts: dict[str, B.TargetTables]) -> dict[str, Any]:
    grid = tab.arrays["grid_e_mev"]

    def table_mb(name: str, e: np.ndarray) -> np.ndarray:
        i = tab.target_names.index(name)
        return 1e3 * np.interp(e, grid, tab.arrays["sigma_barn"][i])

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
