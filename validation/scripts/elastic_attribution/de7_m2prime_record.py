"""D-E7 (exploratory, report-only; Amendments 20 (b'') and 21 (a), (d), (f)): Git-safe projection of the official M2 and M2' builder records.

Inputs: the builder record JSONs of ``python -m ionmc data build elastic-proton --model m2|m2prime`` (``elastic-proton-<table id>.json``;
fields ``exploratory``, ``table_id``, ``m2`` with ``targets/<T>/nodes``, ``checks``, ``trigger_2_convergence``,
``level_sanity_trigger_3``, ``radius_sanity_trigger_4``, ``pauli``, ``provenance``, ``pp_reanchor``,
``sigma_tot_comparison_report_only``) and the output of the researcher's scratch script
``validation/plans/v3-005d-m2prime/scripts/m2prime2.py`` (captured text, or run here with ``IONMC_CACHE`` set).

Projected per target and candidate energy (100 / 150 / 200 MeV; node values interpolated linearly in E as the builder does for the
trigger): sigma_el / Barashenkov (``m2_over_bgg``, the committed BGG/Barashenkov level), sigma_R / Barashenkov (sigma_inel(p + A) of
``bgg.barashenkov_mb``; **withheld for M2**, Amendment 21 (f)), the sigma_tot comparison as pass/fail against the +-10 % band only,
q_min and q_min / q_H, the Pauli factor at peak density; trigger 2 (component maxima, pass), trigger 3 (>= 2 of 3 in [0.7, 1.4] per
target, the margins), trigger 4, the Pauli energy domain and node counts, the pins (SAID / DSG file sha256, inputs JSON sha256), the
table id and the exploratory block, and the builder-vs-scratch reproduction record as RELATIVE differences with the 2 % rule.

Nothing absolute is recorded: no cross section of any variant (mb, barn, fm^2), no measured value, no % deviation from measurement.
``leak_check`` asserts that no projected number equals (to 1e-9 relative) any absolute cross section read from the records or from the
scratch output; it runs on every generation and in the test.

Scratch resolution: ``m2prime2.py`` prints cross sections to 0.1 mb and sigma_el/Barashenkov to 3 digits, which bounds the comparison
resolution (about 1e-3 relative); the scratch Barashenkov sigma_el are the rounded literals of the script.

Run: ``python -m validation.scripts.elastic_attribution.de7_m2prime_record --m2 DIR_OR_JSON --m2prime DIR_OR_JSON [--scratch-text FILE]``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

from . import common as c

STEP_DIR = "de7"
ENERGIES = (100, 150, 200)
BAND3 = (0.7, 1.4)
SIGMA_TOT_BAND = 0.10
REPRO_TOL = 0.02
TARGET_Z = {"C-12": 6, "N-14": 7, "O-16": 8, "Al-27": 13, "Si-28": 14, "P-31": 15, "Ca-40": 20}
SCRATCH_SCRIPT = c.REPO / "validation" / "plans" / "v3-005d-m2prime" / "scripts" / "m2prime2.py"
SCRATCH_LABEL = "R-C Clementel-Villi LDA"
SCRATCH_RE = re.compile(
    r"^(?P<t>C-12|O-16) (?P<e>\d+) " + re.escape(SCRATCH_LABEL) + r"\s+sTot\s+(?P<tot>[\d.]+) \([^)]*\) "
    r"sR\s+(?P<sr>[\d.]+) sEl\s+(?P<sel>[\d.]+) el/Bar (?P<bar>[\d.]+) qmin (?P<q>[\d.]+)",
    re.M,
)
ABS_NODE_KEYS = ("sigma_el_m2_mb", "sigma_tot_m2_mb", "sigma_react_m2_mb", "sigma_el_bgg_mb")
LEAK_CONSTANTS = (0.1, 0.7, 1.4, 0.6, 1.0, 0.02, 0.001, 80.0, 100.0, 150.0, 200.0)  # frozen plan constants and energy labels
LEAK_EXEMPT_KEYS = ("band", "e_mev", "energy_mev", "tolerance", "tol", "pauli_emin_mev", "emin")  # declared configuration, not results
LEAK_FRAGMENTS = ("measured", "model_over", "percent_dev", "deviation_from")

GIT_ALLOW = (
    "variants/*/identity/*",
    "variants/*/exploratory/*",
    "variants/*/withheld/*",
    "variants/*/cells/**",
    "variants/*/sigma_tot_band/**",
    "variants/*/trigger2/**",
    "variants/*/trigger3/**",
    "variants/*/trigger4/**",
    "variants/*/pauli/**",
    "variants/*/pins/**",
    "reproduction/**",
    "leak_check/*",
)
"""Git allowlist of ``de7.json`` (syntax: :mod:`common`): ratios, pass/fail flags, q_min, factors, counts, hashes, relative differences."""
GIT_FILES = {"de7.json": (GIT_ALLOW, ())}


# ---------------------------------------------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------------------------------------------
def load_builder_record(path: Path) -> dict[str, Any]:
    """The builder record JSON (a file, or a directory holding exactly one ``elastic-proton-*.json``)."""
    p = Path(path)
    if p.is_dir():
        files = sorted(p.glob("elastic-proton-*.json"))
        if len(files) != 1:
            raise SystemExit(f"{p}: expected exactly one elastic-proton-*.json, found {len(files)}")
        p = files[0]
    return json.loads(p.read_text())  # type: ignore[no-any-return]


def _finite(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def interp(nodes: list[dict[str, Any]], key: str, energy: float) -> float | None:
    """``key`` linear in E between the nodes (None outside the node range or for a non-finite value)."""
    rows = sorted((r for r in nodes if _finite(r.get(key))), key=lambda r: r["e_mev"])
    if not rows or not rows[0]["e_mev"] <= energy <= rows[-1]["e_mev"]:
        return None
    lo = max((r for r in rows if r["e_mev"] <= energy), key=lambda r: r["e_mev"])
    hi = min((r for r in rows if r["e_mev"] >= energy), key=lambda r: r["e_mev"])
    if hi["e_mev"] == lo["e_mev"]:
        return float(lo[key])
    w = (energy - lo["e_mev"]) / (hi["e_mev"] - lo["e_mev"])
    return float(lo[key] + w * (hi[key] - lo[key]))


def nearest(nodes: list[dict[str, Any]], key: str, energy: float) -> Any:
    """``key`` of the node nearest to ``energy`` (None if absent or non-finite)."""
    if not nodes:
        return None
    v = min(nodes, key=lambda r: abs(r["e_mev"] - energy)).get(key)
    return v if _finite(v) else None


def barashenkov_inelastic_mb(target: str, energy: float) -> float | None:
    from ionmc.nuclear import bgg

    z = TARGET_Z.get(target)
    if z is None:
        return None
    try:
        return float(bgg.barashenkov_mb(energy, z)[1])
    except bgg.BggError:
        return None


def ratio(a: float | None, b: float | None) -> float | None:
    return None if a is None or b is None or b == 0.0 else a / b


# ---------------------------------------------------------------------------------------------------------------------------------
# projection of one builder record
# ---------------------------------------------------------------------------------------------------------------------------------
def project_variant(info: dict[str, Any], *, withhold_sigma_r: bool) -> dict[str, Any]:
    m2 = info["m2"]
    targets = {t: v["nodes"] for t, v in m2["targets"].items()}
    t4 = m2.get("radius_sanity_trigger_4", {})
    comp = m2.get("sigma_tot_comparison_report_only") or {}
    cells: dict[str, Any] = {}
    band: dict[str, Any] = {}
    for t, nodes in targets.items():
        q_h = (t4.get(t) or {}).get("helm", {}).get("q_h_fm_inv")
        cells[t] = {}
        flags: dict[str, bool | None] = {}
        for e in ENERGIES:
            q = nearest(nodes, "q_first_min_fm_inv", e)
            cell: dict[str, Any] = {
                "el_over_barashenkov": interp(nodes, "m2_over_bgg", e),
                "q_min_fm_inv": q,
                "q_min_over_q_h": ratio(q, q_h),
                "pauli_factor_at_peak": nearest(nodes, "pauli_factor_peak", e),
            }
            if not withhold_sigma_r:
                cell["react_over_barashenkov"] = ratio(interp(nodes, "sigma_react_m2_mb", e), barashenkov_inelastic_mb(t, e))
            cells[t][f"{e}_mev"] = cell
            row = (comp.get(t) or {}).get(f"{e}_mev") if isinstance(comp.get(t), dict) else None
            mom = row.get("model_over_measured") if isinstance(row, dict) else None
            flags[f"{e}_mev"] = None if not _finite(mom) else bool(abs(mom - 1.0) <= SIGMA_TOT_BAND)
        known = [v for v in flags.values() if v is not None]
        band[t] = {"pass_per_cell": flags, "band": SIGMA_TOT_BAND,
                   "pass_at_least_2_of_3": (sum(known) >= 2) if len(known) == 3 else None,
                   "pass_all_3": all(known) if len(known) == 3 else None}  # fmt: skip
    return {
        "identity": {"table_id": info.get("table_id"), "npz_sha256": info.get("npz_sha256")},
        "exploratory": {k: (info.get("exploratory") or {}).get(k) for k in
                        ("model", "variant", "pauli_variant", "pauli_emin_mev", "beta", "pp_level", "unfold", "inputs_unpinned")},  # fmt: skip
        "withheld": ({"sigma_react_over_barashenkov": "withheld for M2 (Amendment 21 (f)): the ratios with the Amendment 20 M2 percentages would recover the measured sigma_tot"}
                     if withhold_sigma_r else {}),  # fmt: skip
        "cells": cells,
        "sigma_tot_band": band,
        "trigger2": _trigger2(m2.get("trigger_2_convergence") or {}),
        "trigger3": _trigger3(cells, m2),
        "trigger4": _trigger4(t4),
        "pauli": _pauli(m2),
        "pins": _pins(m2),
    }


def _trigger2(t2: dict[str, Any]) -> dict[str, Any]:
    return {
        "tolerance": t2.get("tol"), "max": t2.get("max"), "passed": t2.get("passed"), "n_evaluated": t2.get("n_evaluated"),
        "component_max": t2.get("component_max"),
        "candidate_energies": [{"target": r["target"], "e_mev": r["e_mev"], "components": r["components"], "max": r["max"]}
                               for r in t2.get("candidate_energies", [])],
    }  # fmt: skip


def _trigger3(cells: dict[str, Any], m2: dict[str, Any]) -> dict[str, Any]:
    lo, hi = BAND3
    rec = (m2.get("checks") or {}).get("m2_over_bgg_sigma_el_at_100_150_200", {})
    out: dict[str, Any] = {"band": list(BAND3), "targets": {}}
    for t in ("C-12", "O-16"):
        rows = []
        for e in ENERGIES:
            r = cells[t][f"{e}_mev"]["el_over_barashenkov"]
            in_band = None if r is None else bool(lo <= r <= hi)
            recorded = (rec.get(t) or {}).get(f"{e}_mev")
            rows.append({
                "e_mev": e, "ratio": r, "in_band": in_band,
                "margin_to_nearest_bound": None if r is None else min(r - lo, hi - r),
                "excess_over_upper_bound_rel": None if r is None else r / hi - 1.0,
                "matches_builder_record": None if r is None or not _finite(recorded) else bool(abs(r - recorded) <= 1e-9),
            })  # fmt: skip
        n_in = sum(1 for r in rows if r["in_band"])
        out["targets"][t] = {"rows": rows, "n_in_band": n_in, "passes": bool(n_in >= 2 and all(r["in_band"] is not None for r in rows))}
    out["passes"] = all(v["passes"] for v in out["targets"].values())
    out["rule"] = "plan 4.2 trigger 3 (Amendment 21 (d)): sigma_el / Barashenkov in [0.7, 1.4] at >= 2 of 100 / 150 / 200 MeV, for C-12 and for O-16"
    return out


def _trigger4(t4: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {"band": t4.get("band"), "targets": {}}
    for t in ("C-12", "O-16"):
        rec = t4.get(t)
        if not rec:
            continue
        rows = {k: {"q_min_fm_inv": v.get("q_min_m2_fm_inv"), "ratio_to_q_h": v.get("ratio_to_q_h"), "in_band": v.get("in_band")}
                for k, v in rec["rows"].items()}  # fmt: skip
        used = next((k for k in ("150_mev", "200_mev") if rows.get(k, {}).get("in_band") is not None), None)
        out["targets"][t] = {"q_h_fm_inv": rec["helm"].get("q_h_fm_inv"), "rows": rows, "verdict_cell": used,
                             "passes": None if used is None else bool(rows[used]["in_band"]),
                             "note": "150 MeV, or 200 MeV where 150 MeV has no true minimum (Amendment 21 (d))"}  # fmt: skip
    ev = [v["passes"] for v in out["targets"].values()]
    out["passes"] = bool(ev) and all(x is True for x in ev)
    return out


def _pauli(m2: dict[str, Any]) -> dict[str, Any]:
    p = m2.get("pauli")
    if not p:
        return {}
    return {
        "pauli_variant": p.get("pauli_variant"), "pauli_scale": p.get("pauli_scale"), "pauli_emin_mev": p.get("pauli_emin_mev"),
        "rho0_fm3": p.get("rho0_fm3"),
        "nodes_applied_per_target": p.get("nodes_applied_per_target"),
        "nodes_below_emin_black_disk_per_target": p.get("nodes_below_emin_black_disk_per_target"),
        "nodes_below_said_range_black_disk_kept": (m2.get("nodes_below_said_range_black_disk_kept") or {}).get("count_per_target"),
        "factor_at_peak_density_recorded": p.get("pauli_factor_at_peak_density"),
    }  # fmt: skip


def _pins(m2: dict[str, Any]) -> dict[str, Any]:
    prov = m2.get("provenance") or {}
    read = (prov.get("said") or {}).get("read") or {}
    dsg = sorted({h for r in ((m2.get("pp_reanchor") or {}).get("scan_rows") or []) for h in r.get("files_sha256", [])})
    return {
        "said_files_sha256": {k: list(v.get("sha256", [])) for k, v in sorted(read.items())},
        "dsg_files_sha256": dsg,
        "inputs_json_sha256": prov.get("inputs_json_sha256"),
        "horiuchi_json_sha256": prov.get("horiuchi_json_sha256"),
        "inputs_unpinned": prov.get("inputs_unpinned"),
    }  # fmt: skip


# ---------------------------------------------------------------------------------------------------------------------------------
# builder-vs-scratch reproduction (relative differences only)
# ---------------------------------------------------------------------------------------------------------------------------------
def parse_scratch(text: str) -> dict[str, dict[int, dict[str, float]]]:
    """The R-C rows of ``m2prime2.py`` output: ``{target: {E: {tot, sr, sel, bar, q}}}`` (absolute values; kept in memory only)."""
    out: dict[str, dict[int, dict[str, float]]] = {}
    for m in SCRATCH_RE.finditer(text):
        out.setdefault(m["t"], {})[int(m["e"])] = {k: float(m[k]) for k in ("tot", "sr", "sel", "bar", "q")}
    return out


def run_scratch(cache_dir: Path) -> str:
    env = {**os.environ, "IONMC_CACHE": str(cache_dir), "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    res = subprocess.run([sys.executable, str(SCRATCH_SCRIPT)], env=env, capture_output=True, text=True, check=True, cwd=str(c.REPO))
    return res.stdout


def reproduction(info: dict[str, Any], scratch: dict[str, dict[int, dict[str, float]]]) -> dict[str, Any]:
    targets = {t: v["nodes"] for t, v in info["m2"]["targets"].items()}
    cells: dict[str, Any] = {}
    worst = 0.0
    complete = True
    for t in ("C-12", "O-16"):
        cells[t] = {}
        for e in ENERGIES:
            s = scratch.get(t, {}).get(e)
            nodes = targets[t]
            b = {"tot": interp(nodes, "sigma_tot_m2_mb", e), "sr": interp(nodes, "sigma_react_m2_mb", e),
                 "sel": interp(nodes, "sigma_el_m2_mb", e), "bar": interp(nodes, "m2_over_bgg", e),
                 "q": nearest(nodes, "q_first_min_fm_inv", e)}  # fmt: skip
            if s is None or any(b[k] is None for k in ("tot", "sr", "sel", "bar")):
                complete = False
                cells[t][f"{e}_mev"] = {"available": False}
                continue
            d = {k: b[k] / s[k] - 1.0 for k in ("tot", "sr", "sel", "bar")}  # builder over scratch, relative
            d["q"] = None if b["q"] is None else b["q"] / s["q"] - 1.0  # q_min may be absent (no true minimum)
            x = max(abs(d["tot"]), abs(d["sr"]), abs(d["sel"]))
            worst = max(worst, x)
            cells[t][f"{e}_mev"] = {
                "available": True,
                "rel_diff_sigma_tot": d["tot"], "rel_diff_sigma_react": d["sr"], "rel_diff_sigma_el": d["sel"],
                "rel_diff_el_over_barashenkov": d["bar"], "rel_diff_q_min": d["q"],
                "max_abs_rel_diff_cross_sections": x, "within_2pct": bool(x <= REPRO_TOL),
            }  # fmt: skip
    return {
        "scratch": f"{SCRATCH_SCRIPT.relative_to(c.REPO)} (reading R-C), as printed",
        "cells": cells, "complete": complete, "max_abs_rel_diff_cross_sections": worst,
        "tolerance": REPRO_TOL,
        "rule_2pct_verdict": ("all six cells within 2 % in sigma_tot, sigma_R, sigma_el" if complete and worst <= REPRO_TOL
                              else "INCOMPLETE: a cell is missing" if not complete
                              else "above 2 % in at least one cell: independent defect review required before the trigger outcome is recorded"),  # fmt: skip
        "note": "relative differences only (builder / scratch - 1); the 2 % rule is evaluated on the absolute values at run time",
    }  # fmt: skip


# ---------------------------------------------------------------------------------------------------------------------------------
# leak check
# ---------------------------------------------------------------------------------------------------------------------------------
def collect_absolutes(infos: list[dict[str, Any]], scratch: dict[str, dict[int, dict[str, float]]]) -> list[float]:
    """Every absolute cross section of the records and the scratch output, in mb, barn and fm^2."""
    vals: list[float] = []
    for info in infos:
        for v in info["m2"]["targets"].values():
            for r in v["nodes"]:
                for k in ABS_NODE_KEYS:
                    if _finite(r.get(k)):
                        vals.append(float(r[k]))
    for rows in scratch.values():
        for r in rows.values():
            vals += [r["tot"], r["sr"], r["sel"]]
    return [x * f for x in vals for f in (1.0, 1e-3, 0.1)]


def _numbers(doc: Any, path: str = "") -> list[tuple[str, float]]:
    out: list[tuple[str, float]] = []
    if isinstance(doc, dict):
        for k, v in doc.items():
            out += _numbers(v, f"{path}/{k}")
    elif isinstance(doc, list):
        for i, v in enumerate(doc):
            out += _numbers(v, f"{path}/{i}")
    elif _finite(doc) and not isinstance(doc, int):
        out.append((path, float(doc)))
    return out


def leak_violations(proj: Any, absolutes: list[float]) -> list[str]:
    """Paths of the projection that hold an absolute cross section (to 1e-9 relative), a cross-section key name, or a measured-value /
    percentage key."""
    bad = [f"{p}: equals an absolute cross section" for p, x in _numbers(proj)
           if p.rsplit("/", 1)[-1] not in LEAK_EXEMPT_KEYS and p.split("/")[-2:-1] != ["band"] and x not in LEAK_CONSTANTS
           and any(a != 0.0 and abs(x - a) <= 1e-9 * abs(a) for a in absolutes)]  # fmt: skip
    bad += c.git_safety_violations(proj)

    def keys(d: Any, path: str = "") -> list[str]:
        out: list[str] = []
        if isinstance(d, dict):
            for k, v in d.items():
                if any(f in str(k).lower() for f in LEAK_FRAGMENTS):
                    out.append(f"{path}/{k}: measured-value or percentage key")
                out += keys(v, f"{path}/{k}")
        elif isinstance(d, list):
            for i, v in enumerate(d):
                out += keys(v, f"{path}/{i}")
        return out

    return bad + keys(proj)


# ---------------------------------------------------------------------------------------------------------------------------------
# document, report, run
# ---------------------------------------------------------------------------------------------------------------------------------
def build_document(m2_info: dict[str, Any], m2p_info: dict[str, Any], scratch_text: str | None, rb_info: dict[str, Any] | None = None) -> dict[str, Any]:
    scratch = parse_scratch(scratch_text) if scratch_text else {}
    doc: dict[str, Any] = {
        "banner": c.BANNER,
        "variants": {"M2": project_variant(m2_info, withhold_sigma_r=True), "M2prime": project_variant(m2p_info, withhold_sigma_r=False),
                     **({"RB": project_variant(rb_info, withhold_sigma_r=False)} if rb_info is not None else {})},  # fmt: skip
        "reproduction": reproduction(m2p_info, scratch) if scratch else {"complete": False, "rule_2pct_verdict": "INCOMPLETE: no scratch output given"},
    }  # fmt: skip
    absolutes = collect_absolutes([m2_info, m2p_info, *([rb_info] if rb_info is not None else [])], scratch)
    proj = c.git_projection(json.loads(json.dumps(doc, default=float)), GIT_ALLOW)
    bad = leak_violations(proj, absolutes)
    if bad:
        raise SystemExit(f"leak check failed: {bad[:8]}")
    doc["leak_check"] = {"passed": True, "absolute_values_checked": len(absolutes)}
    return doc


def _f(x: Any, nd: int = 3) -> str:
    return "-" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def _yn(x: Any) -> str:
    return "-" if x is None else ("pass" if x else "fail")


def render_report(doc: dict[str, Any]) -> str:
    out = ["", "Builder records of the official M2 and M2' constructions (Amendments 20 (b''), 21). Ratios, pass/fail flags, counts and hashes only.", ""]
    for name, v in doc["variants"].items():
        ex, idt = v["exploratory"], v["identity"]
        out += [f"## {name}", "", f"Table id `{idt.get('table_id')}`; model `{ex.get('model')}`, variant `{ex.get('variant')}`, Pauli variant `{ex.get('pauli_variant')}`, "
                f"Pauli E_min {_f(ex.get('pauli_emin_mev'), 1)} MeV, beta `{ex.get('beta')}`, pp level `{ex.get('pp_level')}`, unfold `{ex.get('unfold')}`, "
                f"inputs unpinned: {ex.get('inputs_unpinned')}.", ""]  # fmt: skip
        if v.get("withheld"):
            out += [f"sigma_R / Barashenkov: {next(iter(v['withheld'].values()))}", ""]
        react = not v.get("withheld")
        out += ["| target | E (MeV) | sigma_el / Barashenkov | " + ("sigma_R / Barashenkov | " if react else "") + "sigma_tot in +-10 % band | q_min (1/fm) | q_min / q_H | Pauli factor (peak) |",
                "|---|---|---|" + ("---|" if react else "") + "---|---|---|---|"]  # fmt: skip
        for t, ce in v["cells"].items():
            for k, cell in ce.items():
                e = k.split("_")[0]
                out.append(f"| {t} | {e} | {_f(cell['el_over_barashenkov'])} | " + (f"{_f(cell.get('react_over_barashenkov'))} | " if react else "")
                           + f"{_yn(v['sigma_tot_band'][t]['pass_per_cell'][k])} | {_f(cell['q_min_fm_inv'])} | {_f(cell['q_min_over_q_h'])} | {_f(cell['pauli_factor_at_peak'])} |")  # fmt: skip
        out += ["", "sigma_tot against the measured values, pass/fail per cell (band +-10 %):", ""]
        for t, b in v["sigma_tot_band"].items():
            out.append(f"- {t}: at least 2 of 3: {_yn(b['pass_at_least_2_of_3'])}; all 3: {_yn(b['pass_all_3'])}")
        t2 = v["trigger2"]
        out += ["", f"Trigger 2: max {_f(t2.get('max'), 6)} against {t2.get('tolerance')}, passed: {t2.get('passed')}, nodes evaluated {t2.get('n_evaluated')}. "
                "Component maxima: " + ", ".join(f"{k} {_f(x, 6)}" for k, x in (t2.get("component_max") or {}).items()) + ".", ""]  # fmt: skip
        t3 = v["trigger3"]
        out += [f"Trigger 3 (>= 2 of 3 in [{t3['band'][0]}, {t3['band'][1]}], C-12 and O-16): **{'passes' if t3['passes'] else 'fails'}**.", "",
                "| target | E (MeV) | ratio | in band | margin to nearest bound | excess over upper bound |", "|---|---|---|---|---|---|"]  # fmt: skip
        for t, tr in t3["targets"].items():
            for r in tr["rows"]:
                out.append(f"| {t} | {r['e_mev']} | {_f(r['ratio'], 4)} | {_yn(r['in_band'])} | {_f(r['margin_to_nearest_bound'], 4)} | {_f(r['excess_over_upper_bound_rel'], 4)} |")
            out.append(f"| {t} | verdict | {tr['n_in_band']} of 3 in band | {_yn(tr['passes'])} | | |")
        t4 = v["trigger4"]
        out += ["", f"Trigger 4 (q_min / q_H in {t4.get('band')}): **{_yn(t4['passes'])}**.", ""]
        for t, tr in t4["targets"].items():
            out.append(f"- {t}: q_H {_f(tr['q_h_fm_inv'])} 1/fm; " + "; ".join(
                f"{k.split('_')[0]} MeV q_min {_f(r['q_min_fm_inv'])}, q_min/q_H {_f(r['ratio_to_q_h'])}, in band {_yn(r['in_band'])}" for k, r in tr["rows"].items())
                + f"; verdict cell {tr['verdict_cell']}: {_yn(tr['passes'])}")  # fmt: skip
        p = v.get("pauli")
        if p:
            out += ["", f"Pauli: variant `{p['pauli_variant']}`, scale {p['pauli_scale']}, E_min {p['pauli_emin_mev']} MeV, rho_0 {p['rho0_fm3']} 1/fm^3. "
                    f"Nodes with the factor applied per target: {p['nodes_applied_per_target']}; black-disk nodes below E_min: {p['nodes_below_emin_black_disk_per_target']}; "
                    f"nodes below the SAID range (black disk kept): {p['nodes_below_said_range_black_disk_kept']}."]  # fmt: skip
        pins = v["pins"]
        out += ["", "Pins:", f"- inputs JSON sha256 `{pins.get('inputs_json_sha256')}`; Horiuchi JSON sha256 `{pins.get('horiuchi_json_sha256')}`; inputs unpinned: {pins.get('inputs_unpinned')}"]
        for k, hs in pins.get("said_files_sha256", {}).items():
            out.append(f"- SAID {k}: " + ", ".join(f"`{h}`" for h in hs))
        out += ["- DSG files: " + (", ".join(f"`{h}`" for h in pins.get("dsg_files_sha256", [])) or "none"), ""]
    r = doc["reproduction"]
    out += ["## Builder-vs-scratch reproduction (M2', six trigger cells; relative differences, builder / scratch - 1)", ""]
    if r.get("cells"):
        out += ["| target | E (MeV) | sigma_tot | sigma_R | sigma_el | sigma_el / Barashenkov | q_min | max (tot, R, el) | within 2 % |", "|---|---|---|---|---|---|---|---|---|"]
        for t, ce in r["cells"].items():
            for k, x in ce.items():
                if not x.get("available"):
                    out.append(f"| {t} | {k.split('_')[0]} | n/a | | | | | | |")
                    continue
                out.append(f"| {t} | {k.split('_')[0]} | {x['rel_diff_sigma_tot']:+.4f} | {x['rel_diff_sigma_react']:+.4f} | {x['rel_diff_sigma_el']:+.4f} | "
                           f"{x['rel_diff_el_over_barashenkov']:+.4f} | {'-' if x['rel_diff_q_min'] is None else format(x['rel_diff_q_min'], '+.4f')} | {x['max_abs_rel_diff_cross_sections']:.4f} | {_yn(x['within_2pct'])} |")  # fmt: skip
        out.append("")
    out += [f"2 % rule: **{r['rule_2pct_verdict']}**.", "", f"Leak check (no absolute cross section in the projection): passed {doc['leak_check']['passed']}, "
            f"{doc['leak_check']['absolute_values_checked']} absolute values compared."]  # fmt: skip
    return "\n".join(out) + "\n"


def run(m2: Path, m2prime: Path, scratch_text: Path | None, out: Path | None = None, cache_dir: Path | None = None, rb: Path | None = None) -> dict[str, Any]:
    text = Path(scratch_text).read_text() if scratch_text else run_scratch(Path(cache_dir or c.REPO / ".ionmc-cache"))
    doc = build_document(load_builder_record(m2), load_builder_record(m2prime), text, load_builder_record(rb) if rb else None)
    inputs = [p for d in (m2, m2prime, *([rb] if rb else [])) for p in ([Path(d)] if Path(d).is_file() else sorted(Path(d).glob("elastic-proton-*.json")))]
    inputs += [Path(scratch_text)] if scratch_text else []
    hdr = c.header("D-E7: official M2 / M2' builder records, Git-safe projection",
                   "python -m validation.scripts.elastic_attribution.de7_m2prime_record --m2 DIR --m2prime DIR --scratch-text FILE", inputs,
                   ["Relative differences against the researcher's scratch numbers only; no absolute cross section, no measured value, no % deviation. "
                    "sigma_R / Barashenkov is withheld for M2 (Amendment 21 (f))."])  # fmt: skip
    c.write_outputs(Path(out or c.OUT / STEP_DIR), "de7.json", doc, "de7-report.md", hdr, render_report, GIT_ALLOW)
    return doc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--m2", type=Path, required=True)
    ap.add_argument("--m2prime", type=Path, required=True)
    ap.add_argument("--rb", type=Path, default=None, help="R-B (bdc-pn, report-only) builder record directory or json")
    ap.add_argument("--scratch-text", type=Path, default=None, help="captured stdout of m2prime2.py (default: run it)")
    ap.add_argument("--cache-dir", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    run(a.m2, a.m2prime, a.scratch_text, a.out, a.cache_dir, a.rb)


if __name__ == "__main__":
    main()
