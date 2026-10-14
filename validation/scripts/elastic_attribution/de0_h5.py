"""D-E0 / H5 (exploratory, report-only; Amendment 19): EM-only peak/plateau, ionmc against TOPAS.

Plan hypothesis table row H5 (MCS overlap): "EM-only peak/plateau ionmc against TOPAS from the archived V11 EM-only inputs";
excluded if the EM-only peak/plateau agree within +-2 %; H5 counts only if the difference exceeds 2 % **and** is in the
direction that amplifies the V11 depression difference.

Direction: F = (peak/plateau)_EM+el / (peak/plateau)_EM, so a larger ionmc EM-only peak/plateau than TOPAS (ratio > 1) lowers
F_ionmc relative to F_TOPAS (amplifies Delta_pk); a smaller one reduces it.

Reused definitions (the V11 compare step; nothing re-derived): ionmc partials ``compare_idd_v5.load_ionmc_v11`` (l.769) and
``jackknife_metrics`` (l.186, delete-one batch jackknife); TOPAS EM-only groups of the three frozen seeds
``load_topas_mode_groups`` (l.642) with the runs selected by ``steps_v5c._v11_reference_runs`` (l.670) and the committed cases
of ``steps_v5c._cases_dir`` (l.665), metrics ``seed_metrics`` (l.200; mean of the per-seed metrics, SD / sqrt(n));
``curve_metrics`` (l.168): plateau = mean of the bins with depth in ``PLATEAU_MM`` = 20-60 mm (l.92; bin centres
(i + 0.5) * 0.5 mm), peak = maximum bin, R80 = ``rm.r80`` of the peak-normalised curve; ratio ``ratio_est`` (l.670, delta method).

Run: ``python -m validation.scripts.elastic_attribution.de0_h5``.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

from . import common as c

GIT_ALLOW = ("tolerance", "verdict/*", "per_energy/*/*", "per_energy/*/*/*", "per_energy/*/*/*/*")
"""Git allowlist of ``h5.json`` (syntax: :mod:`common`): EM-only ratios and metrics (no nuclear data), verdicts, flags."""
GIT_FILES = {"h5.json": (GIT_ALLOW, ())}
TOL = 0.02
STEP_DIR = "de0h5"
REF_DIR = c.REPO / ".ionmc-cache" / "reference-runs"


def _cmp() -> Any:
    path = c.REPO / "validation" / "scripts" / "reference" / "compare_idd_v5.py"
    spec = importlib.util.spec_from_file_location("compare_idd_v5", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def compute() -> tuple[dict[str, Any], list[Path]]:
    cmp = _cmp()
    v5c = c.import_transport_module("steps_v5c")
    ion = {e: cmp.load_ionmc_v11(c.V11_PARTIALS[(e, "emonly")], e, "emonly") for e in (150, 200)}
    nz = {e: ion[e][0].shape[1] for e in ion}
    groups = cmp.load_topas_mode_groups(v5c._v11_reference_runs(REF_DIR), nz, v5c._cases_dir())
    ref = c.v11_reference()
    out: dict[str, Any] = {}
    inputs: list[Path] = [c.RUN_V11, *(c.V11_PARTIALS[(e, "emonly")] for e in (150, 200))]
    for e in (150, 200):
        g = groups[(e, "emonly")]
        im, tm = cmp.jackknife_metrics(ion[e][0]), cmp.seed_metrics(g["curves"])
        pp = cmp.ratio_est(im["peak_over_plateau"], tm["peak_over_plateau"])
        pl = cmp.ratio_est(im["plateau"], tm["plateau"])
        pk = cmp.ratio_est(cmp.Est(im["peak_over_plateau"].value * im["plateau"].value, 0.0, 1.0),
                           cmp.Est(tm["peak_over_plateau"].value * tm["plateau"].value, 0.0, 1.0))  # fmt: skip
        r80 = cmp.difference(im["r80_mm"], tm["r80_mm"], relative=False)
        d = pp.value - 1.0
        out[str(e)] = {
            "ionmc_over_topas_peak_over_plateau": cmp._ci90(pp), "ionmc_over_topas_plateau": cmp._ci90(pl),
            "ionmc_over_topas_peak_value_ratio_no_se": pk.value, "r80_difference_mm_ionmc_minus_topas": cmp._ci90(r80),
            "ionmc": {k: {"value": im[k].value, "se": im[k].se} for k in ("peak_over_plateau", "plateau", "r80_mm")},
            "topas": {k: {"value": tm[k].value, "se": tm[k].se, "seeds": len(g["seeds"])} for k in ("peak_over_plateau", "plateau", "r80_mm")},
            "within_2pct_point": bool(abs(d) <= TOL),
            "ci90_inside_2pct": bool(pp.value - cmp.t95(pp.df) * pp.se > 1.0 - TOL and pp.value + cmp.t95(pp.df) * pp.se < 1.0 + TOL),
            "direction": "amplifies Delta_pk (ionmc EM-only peak/plateau above TOPAS)" if d > 0 else "reduces Delta_pk (ionmc EM-only peak/plateau below TOPAS)",
            "amplifies_delta_pk": bool(d > 0),
            "h5_excluded": bool(abs(d) <= TOL),
            "h5_counts": bool(abs(d) > TOL and d > 0),
            "implied_shift_of_F_ratio": 1.0 / pp.value - 1.0,
            "delta_pk_F_ionmc_over_F_topas_minus_1": ref[e]["F_ionmc"] / ref[e]["F_topas"] - 1.0,
            "topas_runs": sorted(g["runs"]),
        }  # fmt: skip
        inputs += [
            REF_DIR / r / "work" / "idd_dose.csv"
            for r in sorted(g["runs"])
            if (REF_DIR / r / "work" / "idd_dose.csv").is_file()
        ]
    doc = {"banner": c.BANNER, "tolerance": TOL, "per_energy": out,
           "verdict": {e: ("H5 excluded (EM-only within +-2 %)" if o["h5_excluded"] else
                           "H5 counts (>2 % and amplifying)" if o["h5_counts"] else "outside +-2 % but not amplifying: H5 does not count")
                       for e, o in out.items()}}  # fmt: skip
    return doc, inputs


def run() -> None:
    doc, inputs = compute()
    cmd = "python -m validation.scripts.elastic_attribution.de0_h5"
    hdr = c.header("D-E0 / H5: EM-only peak/plateau, ionmc against TOPAS (archived V11 EM-only inputs)", cmd, inputs,
                   ["Plateau window 20-60 mm and all metric definitions are those of the V11 compare step (module docstring). "
                    "Peak/plateau and plateau: ionmc batch-jackknife, TOPAS mean of the 3 frozen seeds. 90 % intervals by the "
                    "compare step's Welch df. Only ratios and metrics are recorded."])  # fmt: skip
    c.write_outputs(
        c.OUT / STEP_DIR, "h5.json", doc, f"{STEP_DIR}-report.md", hdr, render_report, GIT_ALLOW
    )


def render_report(doc: dict[str, Any]) -> str:
    out = ["", "| E | peak/plateau ionmc / TOPAS (90 % CI) | plateau ratio | R80 ionmc - TOPAS (mm) | within 2 % (point / CI) | direction | verdict |",
           "|---|---|---|---|---|---|---|"]  # fmt: skip
    for e, o in doc["per_energy"].items():
        p, q, r = (
            o["ionmc_over_topas_peak_over_plateau"],
            o["ionmc_over_topas_plateau"],
            o["r80_difference_mm_ionmc_minus_topas"],
        )
        out.append(
            f"| {e} | {p['value']:.4f} [{p['ci90'][0]:.4f}, {p['ci90'][1]:.4f}] | {q['value']:.4f} [{q['ci90'][0]:.4f}, {q['ci90'][1]:.4f}] | "
            f"{r['value']:+.3f} [{r['ci90'][0]:+.3f}, {r['ci90'][1]:+.3f}] | {o['within_2pct_point']} / {o['ci90_inside_2pct']} | {o['direction']} | {doc['verdict'][e]} |"
        )
    out += [
        "",
        "| E | EM-only contribution to F_ionmc/F_TOPAS - 1 if uncorrected (= 1/ratio - 1) | V11 F_ionmc/F_TOPAS - 1 (Delta_pk, archived) | ionmc peak/plateau (EM-only) | TOPAS peak/plateau (EM-only) |",
        "|---|---|---|---|---|",
    ]
    for e, o in doc["per_energy"].items():
        out.append(
            f"| {e} | {o['implied_shift_of_F_ratio']:+.4f} | {o['delta_pk_F_ionmc_over_F_topas_minus_1']:+.4f} | "
            f"{o['ionmc']['peak_over_plateau']['value']:.4f} +- {o['ionmc']['peak_over_plateau']['se']:.4f} | "
            f"{o['topas']['peak_over_plateau']['value']:.4f} +- {o['topas']['peak_over_plateau']['se']:.4f} |"
        )
    out += [
        "",
        f"Rule (plan H5): excluded if the EM-only peak/plateau agree within +-{doc['tolerance']:.0%}; H5 counts only if the difference "
        "exceeds 2 % and amplifies Delta_pk. F is a ratio, so a common EM difference cancels to first order.",
    ]
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    run()
