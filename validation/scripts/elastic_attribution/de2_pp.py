"""D-E2 (exploratory, report-only; Amendment 19): p-p audit table.

Per energy (20, 50, 100, 150, 200 MeV): sigma_pp of the pinned table (Hale NI, theta_CM >= 16.26 deg; above 150 MeV
anchored times the BGG ratio), the builder-formula value, the BGG formula ``pp_bgg_mb`` (a different quantity: total
elastic p-p), PDG rpp2022 points where the acquired file has any; ``sigma_ET = sigma * <T_slow>`` (energy-transfer
cross section); the slow-proton energy and the angle distribution of the continuing (faster) proton; and the
probability of a p-p event per primary along the water CSDA path.

Run: ``python -m validation.scripts.elastic_attribution.de2_pp``.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from . import common as c
from .de0_kernels import N_Z, R_FRACTION, U, deficit_mm, primary_event, theta_w_deg
from .de1_shapes import Shapes

ENERGIES = (20, 50, 100, 150, 200)
_PE = "per_energy/*"
GIT_ALLOW = (
    "energies_mev", "pdg_covers_20_200_mev", "pdg_min_t_mev_exceeds_300", "pdg_n_points",
    *(f"{_PE}/{k}" for k in (
        "P_deficit_gt_2mm", "P_deficit_gt_3mm", "P_theta_fast_gt_theta_2mm", "P_theta_fast_gt_theta_3mm",
        "mean_1_minus_cos_fast_lab", "mean_t_slow_mev", "mean_t_slow_over_T", "n_pdg_points_within_15pct_in_T",
        "ratio_table_over_bgg")),
    f"{_PE}/theta_fast_cm_deg/*", f"{_PE}/theta_fast_lab_deg/*", "per_primary_water_path/*/*",
)  # fmt: skip
GIT_ALLOW_LONG = ("energies_mev",)
GIT_FILES = {"pp_audit.json": (GIT_ALLOW, GIT_ALLOW_LONG)}
"""Git allowlist of ``pp_audit.json`` (syntax: :mod:`common`): the table / BGG ratio, mean slow-proton energy, angle statistics,
per-primary p-p probabilities. Not kept: the absolute sigma_pp (table, builder formula, BGG), sigma_ET, the PDG points."""


def run() -> None:
    from ionmc.data import cache
    from ionmc.nuclear import bgg

    v10 = c.import_transport_module("nuclear_v10")
    tab = c.load_table()
    sh, water = Shapes(tab), c.Water(tab)
    sec, awi = v10.h1_law5(c.CACHE_DIR)
    pdg_path = cache.verify("pdg-rpp2022-pp-elastic", c.CACHE_DIR)
    pdg = v10.parse_pdg_pp(pdg_path.read_text(encoding="ascii", errors="replace"))
    pdg_t = [
        {
            **r,
            "t_mev": math.sqrt(r["p_lab_gev"] ** 2 * 1e6 + bgg.M_PROTON_MEV**2) - bgg.M_PROTON_MEV,
        }
        for r in pdg
    ]
    rows: dict[str, Any] = {}
    for e in ENERGIES:
        ev = primary_event(sh, "H-1", "table", float(e), U)
        sig = float(tab.sigma_barn(sh.index("H-1"), e)) * 1e3
        t_slow = ev["t_other"]
        th_cm = np.degrees(
            np.arccos(np.clip(np.abs(ev["mu"]), 0, 1))
        )  # the faster proton has cos_CM >= 0
        th_lab = np.degrees(np.arccos(np.clip(ev["cos"], -1, 1)))
        d = deficit_mm(float(e), ev, water.range_mm)
        near = [r for r in pdg_t if abs(r["t_mev"] - e) / e < 0.15]
        rows[str(e)] = {
            "sigma_pp_table_mb": sig, "sigma_pp_builder_formula_mb": float(v10.sigma_pp_direct_mb(sec, awi, float(e))),
            "sigma_pp_bgg_formula_mb": float(bgg.pp_bgg_mb(float(e))),
            "ratio_table_over_bgg": sig / float(bgg.pp_bgg_mb(float(e))),
            "pdg_points_within_15pct_in_T": [{"t_mev": r["t_mev"], "sigma_mb": r["sigma_mb"]} for r in near],
            "n_pdg_points_within_15pct_in_T": len(near),
            "mean_t_slow_mev": float(np.mean(t_slow)), "mean_t_slow_over_T": float(np.mean(t_slow)) / e,
            "sigma_ET_mb_mev": sig * float(np.mean(t_slow)),
            "theta_fast_lab_deg": {f"p{q}": float(np.percentile(th_lab, q)) for q in (10, 50, 90)},
            "theta_fast_cm_deg": {f"p{q}": float(np.percentile(th_cm, q)) for q in (10, 50, 90)},
            "mean_1_minus_cos_fast_lab": float(np.mean(1 - ev["cos"])),
            "P_theta_fast_gt_theta_2mm": float(np.mean(th_lab > theta_w_deg(e, 2.0, water.range_mm))),
            "P_theta_fast_gt_theta_3mm": float(np.mean(th_lab > theta_w_deg(e, 3.0, water.range_mm))),
            "P_deficit_gt_2mm": float(np.mean(d > 2.0)), "P_deficit_gt_3mm": float(np.mean(d > 3.0)),
        }  # fmt: skip
    # probability of p-p per primary and transferred energy along the water path
    path: dict[str, Any] = {}
    for e0 in (100, 150, 200):
        r0 = float(water.range_mm(e0))
        zc = (np.arange(N_Z) + 0.5) / N_Z * R_FRACTION * r0
        dz = R_FRACTION * r0 / N_Z / 10.0
        p_t = p_b = p_t_hi = et = 0.0
        for en in water.energy_at_depth(float(e0), zc):
            en = float(en)
            nh = water.targets["H-1"]["n_per_cm3"] * water.targets["H-1"]["scale"]
            s_t = float(tab.sigma_barn(sh.index("H-1"), en)) * 1e-24
            p_t += nh * s_t * dz
            p_t_hi += nh * s_t * dz * (en > 50.0)
            if en < 275.0:
                p_b += nh * float(bgg.pp_bgg_mb(en)) * 1e-27 * dz * (en >= 12.5)
            if s_t > 0:
                et += (
                    nh
                    * s_t
                    * dz
                    * float(np.mean(primary_event(sh, "H-1", "table", en, U)["t_other"]))
                )
        path[str(e0)] = {"P_pp_table": p_t, "P_pp_table_above_50MeV": p_t_hi, "P_pp_bgg_formula_E_ge_12.5": p_b,
                         "transferred_energy_mev_per_primary": et}  # fmt: skip
    pdg_low = sorted(pdg_t, key=lambda r: r["t_mev"])[:5]
    doc = {"banner": c.BANNER, "energies_mev": list(ENERGIES), "per_energy": rows, "per_primary_water_path": path,
           "pdg_lowest_points": [{"t_mev": r["t_mev"], "p_lab_gev": r["p_lab_gev"], "sigma_mb": r["sigma_mb"]} for r in pdg_low],
           "pdg_covers_20_200_mev": bool(min(r["t_mev"] for r in pdg_t) <= 200.0),
           "pdg_min_t_mev_exceeds_300": bool(min(r["t_mev"] for r in pdg_t) > 300.0),
           "pdg_n_points": len(pdg_t),
           "note": "BGG formula = total elastic p-p; the table sigma is the NI part above 16.26 deg CM (not the same quantity)"}  # fmt: skip
    inputs = [*c.table_files(), pdg_path]
    c.write_outputs(
        c.OUT / "de2",
        "pp_audit.json",
        doc,
        "de2-report.md",
        _header(inputs),
        render_report,
        GIT_ALLOW,
        GIT_ALLOW_LONG,
    )


def _header(inputs: list[Any]) -> str:
    cmd = "python -m validation.scripts.elastic_attribution.de2_pp"
    return c.header(
        "D-E2: p-p audit",
        cmd,
        inputs,
        [
            "Table sigma = Hale NI above 16.26 deg CM (above 150 MeV anchor x BGG ratio); BGG formula = total elastic, a different quantity "
            "(absolute values are not committed: decision 0041). Angles: the continuing (faster) proton. No download was made; PDG = the cached rpp2022 file."
        ],
    )


def render_report(doc: dict[str, Any]) -> str:
    """Markdown body, a pure function of the (projected) JSON record."""
    out: list[str] = []
    out += [
        "",
        "| E | table / BGG formula | <T_slow> MeV | <T_slow>/E | th_fast lab p10/p50/p90 | th_fast CM p50 | P(th>th_2mm) | P(th>th_3mm) | PDG within 15 % |",
        "|" + "---|" * 9,
    ]
    for e, r in sorted(doc["per_energy"].items(), key=lambda kv: float(kv[0])):
        t = r["theta_fast_lab_deg"]
        pdg = str(r["n_pdg_points_within_15pct_in_T"])
        out.append(
            f"| {e} | {r['ratio_table_over_bgg']:.3f} | "
            f"{r['mean_t_slow_mev']:.2f} | {r['mean_t_slow_over_T']:.3f} | {t['p10']:.1f} / {t['p50']:.1f} / {t['p90']:.1f} | {r['theta_fast_cm_deg']['p50']:.1f} | "
            f"{r['P_theta_fast_gt_theta_2mm']:.3f} | {r['P_theta_fast_gt_theta_3mm']:.3f} | {pdg} |"
        )
    out += [
        "",
        "## p-p events per primary along the water CSDA path",
        "",
        "| E0 | P(pp) table | P(pp) table, E>50 MeV | P(pp) BGG formula (E >= 12.5) | transferred energy (MeV/primary) |",
        "|---|---|---|---|---|",
    ]
    for e, r in sorted(doc["per_primary_water_path"].items(), key=lambda kv: float(kv[0])):
        out.append(
            f"| {e} | {r['P_pp_table']:.4f} | {r['P_pp_table_above_50MeV']:.4f} | {r['P_pp_bgg_formula_E_ge_12.5']:.4f} | {r['transferred_energy_mev_per_primary']:.3f} |"
        )
    out += [
        "",
        f"PDG file: {doc['pdg_n_points']} points; covers 20-200 MeV: {doc['pdg_covers_20_200_mev']}; "
        f"lowest point above 300 MeV: {doc['pdg_min_t_mev_exceeds_300']} (values and energies not reproduced, decision 0041).",
    ]
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    run()
