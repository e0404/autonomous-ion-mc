"""D-E4 analysis (exploratory, report-only; Amendment 19 (c), plan section 3.2): TOPAS phase space -> per-class summary.

Input: the TOPAS binary phase space of an EM-free ``g4h-elastic_HP`` thin-target case (``ps.header`` + ``ps.phsp``; protons
leaving the 500 mm vacuum sphere), one or several parts of one case (``--parts DIR1 DIR2 DIR3``, concatenated). The case is
identified from the materialized inputs against the frozen table :data:`de5_sampling.CASES` (case index 0-9 = rows 14-23): target
material (``G4_WATER`` -> p+O-16 and p-p; ``G4_C`` -> p+C-12), thickness (2 x ``Ge/Target/HL``), beam energy and the sum of the
parts' histories are read from ``input.txt`` and must match one entry, and the ``ps.header`` maximum proton energy must equal E0. The record layout is parsed from
``ps.header`` (no hard-coded offsets; the byte count per record is asserted to equal the sum of the field sizes and the file size
to equal ``records x bytes``).

Analysis (plan 3.2):
- exit angle ``theta_lab = arccos(W)``, ``W = +-sqrt(1 - U^2 - V^2)`` (the sign flag of the header);
- kinetic-energy filter: primary records (parent ID 0) with ``E > E0 - 1e-4 MeV`` are the undeflected primaries and are removed
  (OpenTOPAS has no current-KE filter, so they are recorded and removed here); the removed fraction is reported;
- "has a secondary proton" of an event: the event has a proton record with **parent ID != 0** (the p-p partner; the number of events
  with more than one proton record is reported as a cross-check);
- :func:`de5_sampling.select_tracks` (the same function as the table side): KE filter, classification (smaller kinematic-curve
  residual to p+H / p+target wins; a conflict track, both residuals < 1e-3, is p-p iff the event has a secondary proton, else p+A)
  and the frozen **exclusion of tracks off both curves** (counted, fraction reported). p-p events enter the H-1 sample only if
  the event has exactly two proton records and both lie on the p+H curve at E0 (others are counted as ``pp_events_off_curve``);
- p-p events (event with a secondary proton): the faster of the two protons is the "fast proton" (the primary of the ionmc
  convention), the slower the "slow/secondary proton"; the class ``H-1`` signatures use the fast proton, ``mean_t_other_mev`` the
  slow proton's energy; p+A (``O-16``): ``mean_t_other_mev`` = ``E0 - T1`` (the recoil, energy conservation);
- per-class rates per primary ``N_class / N_primaries`` with the Poisson standard error;
- the signature values with the same constants and functions as :mod:`de5_sampling` (thresholds ``theta_w`` from the Bethe range of
  water, ``theta_med``, tail moment ``T`` over 3 theta_med, the descriptors ``shape_descriptors`` on the CM angle reconstructed
  from ``theta_lab`` and the energy with ``nuclear_v10.lab_to_cm``; recorded, not interpreted).

``geant4_summary`` has the shape ``{str(case): {class: {signature: value}}}`` that
``de5_sampling.signature_ratios(table_doc, geant4_summary)`` takes. Parts of one case are concatenated (counts add, primaries add).

Outputs (``--out``, default ``validation/generated/elastic_attribution/de4``): ``de4-case<N>.json`` per case through
``common.write_outputs`` (full record + Git projection with its Markdown), and ``--combine`` writes ``de4.json`` / ``de4-report.md``
from every ``de4-case*.json`` of the directory.

Run: ``python -m validation.scripts.elastic_attribution.de4_parse --parts DIR [DIR ...]`` and ``... --combine``.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Any

import numpy as np

from . import common as c
from . import de5_sampling as de5

NAMES = {
    "Position X": ("X", "f4"), "Position Y": ("Y", "f4"), "Position Z": ("Z", "f4"), "Direction Cosine X": ("U", "f4"),
    "Direction Cosine Y": ("V", "f4"), "Energy": ("E", "f4"), "Weight": ("weight", "f4"), "Particle Type": ("pdg", "i4"),
    "Flag to tell if Third Direction Cosine is Negative": ("w_negative", "u1"),
    "Flag to tell if this is the First Scored Particle from this History": ("first_in_history", "u1"),
    "Event ID": ("event", "i4"), "Track ID": ("track", "i4"), "Parent ID": ("parent", "i4"),
}  # fmt: skip
SIZE = {"f4": 4, "i4": 4, "u1": 1}
PROTON = 2212
_CL = "classes/*"
GIT_ALLOW_CASE = (
    "case", "row", "material", "target", "energy_mev", "length_cm", "primaries", "records", "secondary_definition",
    "pp_events_all", "pp_events_off_curve", "pp_events_after_cm_cut", "pp_cm_cut_deg", "pp_events_removed_by_cm_cut",
    "theta_2mm_deg", "theta_3mm_deg", "ke_filter/*", "classification/*", "secondary_on_p+H_curve/*", "off_curve/*",
    "identification/*", "geant4_summary/*/*/*",
    *(f"{_CL}/{k}" for k in (
        "events", "rate_per_primary", "rate_per_primary_se", "rate_per_primary_after_rule", "theta_lab_median_deg",
        "P_theta_lab_gt_theta_2mm", "P_theta_lab_gt_theta_3mm", "T_tail_moment", "theta_lab_mean_1_minus_cos",
        "mean_t_other_mev", "first_minimum_deg", "first_minimum_3sigma_deg", "qslope", "qslope_bins",
        "counts_in_first_minimum_bin", "slow_proton_theta_lab_median_deg")),
    f"{_CL}/helm/*",
)  # fmt: skip
GIT_ALLOW = tuple(f"cases/*/{p}" for p in GIT_ALLOW_CASE)  # the combined record de4.json
"""Git allowlists (syntax: :mod:`common`): counts, per-primary rates, angle statistics, exceedances, identification of the case.
Not kept: histograms, the part directories (host paths)."""
GIT_FILES = {
    "de4-case*.json": (GIT_ALLOW_CASE, ()),
    "de4.json": (GIT_ALLOW + ("geant4_summary/*/*/*",), ()),
}


def parse_header(path: Path) -> dict[str, Any]:
    """Record layout, number of histories and of scored particles from a TOPAS binary ``.header``."""
    text = Path(path).read_text()
    fields = []
    for line in text.splitlines():
        m = re.match(r"^(f4|i4|b1):\s*(.*?)\s*(\[.*\])?$", line)
        if m:
            key = next((k for k in NAMES if m.group(2).startswith(k)), None)
            if key is None:
                raise ValueError(f"unknown header field {line!r}")
            fields.append(NAMES[key])
    nb = int(re.search(r"Number of Bytes per Particle:\s*(\d+)", text).group(1))  # type: ignore[union-attr]
    hist = int(re.search(r"Number of Original Histories:\s*(\d+)", text).group(1))  # type: ignore[union-attr]
    scored = int(re.search(r"Number of Scored Particles:\s*(\d+)", text).group(1))  # type: ignore[union-attr]
    assert sum(SIZE[t] for _, t in fields) == nb, "header byte count != sum of field sizes"
    return {"fields": fields, "bytes": nb, "histories": hist, "scored": scored}


def read_phsp(directory: Path) -> tuple[np.ndarray, dict[str, Any]]:
    """Structured array of the records of ``directory/ps.phsp`` (layout from ``ps.header``)."""
    d = Path(directory)
    h = parse_header(d / "ps.header")
    dt = np.dtype([(n, "<" + t) for n, t in h["fields"]])
    assert dt.itemsize == h["bytes"]
    raw = (d / "ps.phsp").read_bytes()
    assert len(raw) == h["scored"] * h["bytes"], "file size != scored particles x bytes per record"
    return np.frombuffer(raw, dtype=dt), h


def beam_energy_mev(directory: Path) -> float:
    return float(parse_input(directory)["energy_mev"])


def parse_input(directory: Path) -> dict[str, Any]:
    """Target material, thickness, beam energy and histories of a materialized TOPAS case (``input.txt``)."""
    text = (Path(directory) / "input.txt").read_text()

    def grab(pat: str) -> str:
        m = re.search(pat, text, re.M)
        assert m, f"{directory}: input.txt has no match for {pat!r}"
        return m.group(1)

    hl = re.search(r"^d:Ge/Target/HL\s*=\s*([0-9.eE+-]+)\s*(mm|cm)", text, re.M)
    assert hl, f"{directory}: no Ge/Target/HL"
    length_cm = 2.0 * float(hl.group(1)) * (0.1 if hl.group(2) == "mm" else 1.0)
    return {"material": grab(r'^s:Ge/Target/Material\s*=\s*"([^"]+)"'),
            "length_cm": length_cm, "energy_mev": float(grab(r"^d:So/Beam/BeamEnergy\s*=\s*([0-9.eE+-]+)\s*MeV")),
            "histories": int(grab(r"^i:So/Beam/NumberOfHistoriesInRun\s*=\s*(\d+)"))}  # fmt: skip


def header_max_ke(directory: Path) -> float:
    m = re.search(
        r"Maximum Kinetic Energy of proton:\s*([0-9.eE+-]+)\s*MeV",
        (Path(directory) / "ps.header").read_text(),
    )
    assert m, f"{directory}: no maximum proton energy in ps.header"
    return float(m.group(1))


def identify_case(
    dirs: list[Path], case: int | None = None
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Match the parts against :data:`de5_sampling.CASES`; all parts must agree, their primaries must add up to the case and the
    phase-space header maximum energy must equal E0. Returns ``(case entry, identification record)``."""
    infos = [parse_input(d) for d in dirs]
    key = {(i["material"], i["energy_mev"], round(i["length_cm"], 6)) for i in infos}
    assert len(key) == 1, f"parts differ in material / energy / thickness: {sorted(key)}"
    mat, e0, length = next(iter(key))
    hits = [
        cs
        for cs in de5.CASES
        if cs["topas_material"] == mat
        and cs["energy_mev"] == e0
        and abs(cs["length_cm"] - length) < 1e-6
    ]
    assert len(hits) == 1, f"no unique case for {mat}, {e0} MeV, {length} cm"
    cs = hits[0]
    if case is not None:
        assert case == cs["case"], (
            f"--case {case} but the inputs are case {cs['case']} (row {cs['row']})"
        )
    total = sum(i["histories"] for i in infos)
    assert total == cs["primaries"], (
        f"parts hold {total} primaries, case {cs['case']} needs {cs['primaries']}"
    )
    for d in dirs:
        assert abs(header_max_ke(d) - e0) <= 1e-4 * e0, (
            f"{d}: ps.header maximum energy {header_max_ke(d)} != E0 {e0}"
        )
    ident = {"material_topas": mat, "length_cm": length, "energy_mev": e0, "parts": len(dirs), "primaries_total": total,
             "header_max_ke_equals_e0": True}  # fmt: skip
    return cs, ident


def analyse(dirs: list[Path], case: int | None = None) -> dict[str, Any]:
    from ionmc.nuclear.bgg import M_PROTON_MEV as MP

    nv10 = c.import_transport_module("nuclear_v10")
    tab = c.load_table()
    water = c.Water(tab)
    names = tuple(tab.target_names)
    cs, ident = identify_case(dirs, case)
    case, e0, target = cs["case"], float(cs["energy_mev"]), cs["target"]
    has_h = cs["material"] == "water"
    n_prim, recs = 0, []
    for k, d in enumerate(dirs):
        r, h = read_phsp(d)
        assert (r["pdg"] == PROTON).all(), "non-proton record"
        n_prim += h["histories"]
        r = r.copy()
        r["event"] = r["event"] + k * (1 << 24)  # keep events of different parts apart
        recs.append(r)
    r = np.concatenate(recs)
    cos = np.sqrt(
        np.maximum(1.0 - r["U"].astype(float) ** 2 - r["V"].astype(float) ** 2, 0.0)
    ) * np.where(r["w_negative"] == 1, -1.0, 1.0)
    e = r["E"].astype(float)
    prim = r["parent"] == 0
    sec = ~prim
    ev_all, ev_cnt = np.unique(r["event"], return_counts=True)
    ev_sec = np.unique(r["event"][sec])
    ev_has_sec = np.isin(r["event"], ev_sec)
    mh = float(tab.arrays["target_mass_mev"][names.index("H-1")])
    ma = float(tab.arrays["target_mass_mev"][names.index(target)])
    # primaries: KE filter, classification and the frozen off-curve exclusion (the same function as the table side)
    pi = np.nonzero(prim)[0]
    sel = de5.select_tracks(cos[pi], e[pi], e0, mh, ma, ev_has_sec[pi])
    ke_removed, off, conflict, is_pp = sel["ke_removed"], sel["off"], sel["conflict"], sel["is_pp"]
    surv = ~ke_removed
    pa_keep = sel["keep"] & ~is_pp & ~ev_has_sec[pi]
    pa_e, pa_cos = e[pi][pa_keep], cos[pi][pa_keep]
    # p-p events (water only): exactly two proton records, both on the p+H curve at E0, then the CM cut on the fast proton
    order = np.argsort(r["event"], kind="stable")
    ev_sorted = r["event"][order]
    fast_e, fast_cos, slow_e, slow_cos = [], [], [], []
    n_multi = n_off_pp = 0
    for ev in ev_sec:
        lo, hi = np.searchsorted(ev_sorted, [ev, ev + 1])
        idx = order[lo:hi]
        if idx.size != 2:
            n_multi += 1
            continue
        idx = idx[np.argsort(-e[idx])]
        res = np.abs(e[idx] - de5.kin_curve_t(e0, mh, cos[idx])) / np.maximum(e[idx], 1e-9)
        if (res >= 1.0e-3).any():
            n_off_pp += 1
            continue
        (
            fast_e.append(e[idx[0]]),
            fast_cos.append(cos[idx[0]]),
            slow_e.append(e[idx[1]]),
            slow_cos.append(cos[idx[1]]),
        )
    fast_e, fast_cos, slow_e, slow_cos = map(np.asarray, (fast_e, fast_cos, slow_e, slow_cos))
    gam = nv10.two_body_frame(e0, MP)["gamma_cm"]
    th_cm_fast = 2.0 * np.degrees(np.arctan(gam * np.tan(np.arccos(np.clip(fast_cos, -1.0, 1.0)))))
    keep_cm = th_cm_fast >= de5.PP_CUT_DEG  # 16.26 deg CM cut (applied to Geant4 as to the table)
    n_pp_ok = int(fast_e.size)
    fast_e, fast_cos, slow_e, slow_cos = (
        fast_e[keep_cm],
        fast_cos[keep_cm],
        slow_e[keep_cm],
        slow_cos[keep_cm],
    )
    r_w = float(water.range_mm(e0))
    th_w = {x: math.degrees(math.acos(1.0 - x / r_w)) for x in (2.0, 3.0)}
    pcm = float(
        np.interp(e0, tab.arrays["grid_e_mev"], tab.arrays["p_cm_mev"][names.index(target)])
    )

    def cls_doc(
        t1: np.ndarray, cs_: np.ndarray, t_other: np.ndarray, cm_mass: float | None
    ) -> dict[str, Any]:
        count = int(t1.size)
        th = np.degrees(np.arccos(np.clip(cs_, -1.0, 1.0)))
        med = float(np.median(th)) if th.size else float("nan")
        nan = float("nan")
        doc: dict[str, Any] = {
            "events": count, "rate_per_primary": count / n_prim, "rate_per_primary_after_rule": count / n_prim,
            "rate_per_primary_se": math.sqrt(max(count, 1)) / n_prim,
            "theta_lab_median_deg": med, "P_theta_lab_gt_theta_2mm": float(np.mean(th > th_w[2.0])) if th.size else nan,
            "P_theta_lab_gt_theta_3mm": float(np.mean(th > th_w[3.0])) if th.size else nan,
            "T_tail_moment": float(np.mean((1.0 - cs_) * (th > 3.0 * med))) if th.size else nan,
            "theta_lab_mean_1_minus_cos": float(np.mean(1.0 - cs_)) if th.size else nan,
            "mean_t_other_mev": float(np.mean(t_other)) if th.size else nan,
            "hist_theta_lab_counts_8bins": np.histogram(th, bins=[0.0, *de5.HIST_INTERIOR_DEG, 180.0])[0].tolist(),
            "hist_theta_lab_interior_edges_deg": de5.HIST_INTERIOR_DEG,
            "hist_t1_over_e0_counts_8bins": np.histogram(t1 / e0, bins=[0.0, *de5.HIST_INTERIOR_T, 1.0 + 1e-9])[0].tolist(),
            "hist_t1_over_e0_interior_edges": de5.HIST_INTERIOR_T,
            "hist_theta_lab_1deg_full": np.histogram(th, bins=np.arange(0.0, 181.0))[0].tolist(),
            "first_minimum_deg": None, "first_minimum_3sigma_deg": None, "qslope": None,
        }  # fmt: skip
        if cm_mass is not None and th.size > 50:
            tcm, _ = nv10.lab_to_cm(th, e0, cm_mass)
            doc.update(de5.shape_descriptors(np.cos(np.radians(tcm)), pcm))
        return doc

    classes = {target: cls_doc(pa_e, pa_cos, e0 - pa_e, ma)}
    if has_h:
        classes["H-1"] = cls_doc(fast_e, fast_cos, slow_e, None)
        classes["H-1"]["slow_proton_theta_lab_median_deg"] = (
            float(np.median(np.degrees(np.arccos(np.clip(slow_cos, -1, 1)))))
            if slow_cos.size else float("nan"))  # fmt: skip
    nsurv = int(surv.sum())
    pool = surv & ~ev_has_sec[pi]  # deflected primaries without a secondary proton: the p+A pool
    n_off_pa = int((off & ~ev_has_sec[pi]).sum())
    cl = {
        "primary_records": int(prim.sum()), "deflected_primaries": nsurv, "off_curve_tracks_excluded": int(off.sum()),
        "off_curve_pa_pool_excluded": n_off_pa,
        "off_curve_fraction": float(n_off_pa / pool.sum()) if pool.any() else float("nan"),
        "assigned_pp": int((is_pp & surv).sum()), "assigned_pA": int((~is_pp & surv).sum()),
        "conflict_tracks_raw": int((conflict & surv).sum()),
        "conflict_fraction_raw": float((conflict & surv).sum() / nsurv) if nsurv else float("nan"),
        "conflict_assigned_pp_by_secondary": int((conflict & surv & is_pp).sum()),
        "pA_events_with_secondary_proton": int((~is_pp & surv & ev_has_sec[pi]).sum()),
        "curve_crossover_deg": de5.curve_crossover_deg(e0, mh, ma),
    }  # fmt: skip
    sec_res = np.abs(e[sec] - de5.kin_curve_t(e0, mh, cos[sec])) / np.maximum(e[sec], 1e-9)
    ev_multi_n = int((ev_cnt > 1).sum())
    doc = {
        "banner": c.BANNER, "case": case, "row": cs["row"], "material": cs["material"], "target": target,
        "energy_mev": e0, "length_cm": cs["length_cm"], "primaries": n_prim, "parts": [str(d) for d in dirs],
        "records": int(r.size), "identification": ident,
        "ke_filter": {"threshold_mev": e0 - 1.0e-4, "removed_primaries": int(ke_removed.sum()),
                      "removed_fraction": float(ke_removed.sum() / n_prim), "primary_records": int(prim.sum()),
                      "secondary_records_parent_ne_0": int(sec.sum())},  # fmt: skip
        "secondary_definition": "event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = "
                                f"{ev_multi_n}); a p-p event needs exactly two proton records",
        "off_curve": {"primaries_excluded": n_off_pa, "primaries_fraction": cl["off_curve_fraction"],
                      "pp_events_excluded": n_off_pp, "pp_events_not_two_records": n_multi},  # fmt: skip
        "pp_events_all": int(ev_sec.size), "pp_events_off_curve": n_off_pp + n_multi, "pp_events_after_cm_cut": int(fast_e.size),
        "pp_cm_cut_deg": de5.PP_CUT_DEG, "pp_events_removed_by_cm_cut": n_pp_ok - int(fast_e.size),
        "classification": cl, "classes": classes,
        "secondary_on_p+H_curve": {"max_rel_residual": float(sec_res.max()) if sec_res.size else float("nan"),
                                   "fraction_below_1e-3": float(np.mean(sec_res < 1e-3)) if sec_res.size else float("nan")},
        "theta_2mm_deg": th_w[2.0], "theta_3mm_deg": th_w[3.0],
    }  # fmt: skip
    doc["geant4_summary"] = {
        str(case): {k: {s: v.get(s) for s in de5.SIGNATURES} for k, v in classes.items()}
    }
    return doc


def render_case(doc: dict[str, Any]) -> str:
    k, cl, oc = doc["ke_filter"], doc["classification"], doc["off_curve"]
    out = ["", f"Case {doc['case']} (row {doc['row']}: {doc['material']}, {doc['energy_mev']:g} MeV, {doc['length_cm']:g} cm, "
           f"{doc['primaries']} primaries, {doc['records']} proton records, target {doc['target']}).", "",
           f"KE filter (E > E0 - 1e-4 MeV): {k['removed_primaries']} primaries removed (fraction {k['removed_fraction']:.5f}); "
           f"{k['primary_records']} primary records, {k['secondary_records_parent_ne_0']} secondary proton records (parent != 0). "
           f"Secondary definition: {doc['secondary_definition']}.", "",
           "| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |", "|---|---|---|---|---|---|---|---|"]  # fmt: skip
    for name, v in doc["classes"].items():
        out.append(
            f"| {name} | {v['events']} | {v['rate_per_primary']:.4g} +- {v['rate_per_primary_se']:.2g} | {v['theta_lab_median_deg']:.2f} | "
            f"{v['P_theta_lab_gt_theta_2mm']:.4f} | {v['P_theta_lab_gt_theta_3mm']:.4f} | {v['T_tail_moment']:.4f} | {v['mean_t_other_mev']:.3f} |"
        )
    out += ["", f"p-p events: {doc['pp_events_all']} with a secondary proton; {oc['pp_events_excluded'] + oc['pp_events_not_two_records']} excluded "
            f"(off the p+H curve or not exactly two proton records); {doc['pp_events_removed_by_cm_cut']} removed by the {doc['pp_cm_cut_deg']} deg CM cut; "
            f"{doc['pp_events_after_cm_cut']} kept (class H-1).", "",
            f"Classification: {cl['deflected_primaries']} deflected primaries -> {cl['assigned_pA']} p+A, {cl['assigned_pp']} p-p; "
            f"off both kinematic curves: {cl['off_curve_pa_pool_excluded']} p+A-pool tracks excluded (fraction {cl['off_curve_fraction']:.3e}; "
            f"{cl['off_curve_tracks_excluded']} including primaries of p-p events); raw conflict tracks "
            f"{cl['conflict_tracks_raw']} ({cl['conflict_fraction_raw']:.3e}; crossover {cl['curve_crossover_deg']:.2f} deg), of which "
            f"{cl['conflict_assigned_pp_by_secondary']} p-p by the secondary rule. Secondary protons on the p+H curve: "
            f"max residual {doc['secondary_on_p+H_curve']['max_rel_residual']:.2e}, fraction < 1e-3 {doc['secondary_on_p+H_curve']['fraction_below_1e-3']:.4f}.", "",
            "The first-minimum, q-slope and Helm descriptors are recorded, not interpreted."]  # fmt: skip
    return "\n".join(out) + "\n"


def render_report(doc: dict[str, Any]) -> str:
    """Markdown body of the combined record ``de4.json`` (one section per case)."""
    out: list[str] = []
    for n in sorted(doc["cases"], key=int):
        out += ["", f"## Case {n}", render_case(doc["cases"][n])]
    return "\n".join(out) + "\n"


def _cmd(args: str) -> str:
    return f"python -m validation.scripts.elastic_attribution.de4_parse {args}"


def run(dirs: list[Path], outdir: Path, case: int | None = None) -> dict[str, Any]:
    doc = analyse(dirs, case)
    inputs = [
        p for d in dirs for p in (d / "ps.header", d / "ps.phsp", d / "input.txt") if p.is_file()
    ]
    name = f"de4-case{doc['case']}"
    c.write_outputs(outdir, f"{name}.json", doc, f"{name}-report.md",
                    c.header(f"D-E4 analysis, case {doc['case']}: TOPAS phase space -> per-class summary",
                             _cmd("--parts " + " ".join(rel_path(d) for d in dirs)), inputs),
                    render_case, GIT_ALLOW_CASE)  # fmt: skip
    return doc


def rel_path(p: Path) -> str:
    return c.rel(Path(p))


def combine(outdir: Path) -> dict[str, Any]:
    """``de4.json`` / ``de4-report.md`` from every ``de4-case*.json`` (full records) of ``outdir``."""
    import json

    files = sorted(p for p in Path(outdir).glob("de4-case*.json"))
    cases = {str(json.loads(p.read_text())["case"]): json.loads(p.read_text()) for p in files}
    doc = {"banner": c.BANNER, "cases": cases,
           "geant4_summary": {k: v["geant4_summary"][k] for k, v in cases.items()}}  # fmt: skip
    allow = GIT_ALLOW + ("geant4_summary/*/*/*",)
    c.write_outputs(outdir, "de4.json", doc, "de4-report.md",
                    c.header("D-E4 analysis: combined per-case summaries", _cmd("--combine"), files), render_report, allow)  # fmt: skip
    return doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument(
        "--parts",
        type=Path,
        nargs="+",
        help="work directories (ps.header, ps.phsp, input.txt) of the parts of ONE case",
    )
    ap.add_argument(
        "--case",
        type=int,
        default=None,
        help="optional: asserted against the case identified from the inputs",
    )
    ap.add_argument(
        "--combine",
        action="store_true",
        help="write de4.json / de4-report.md from the de4-case*.json in --out",
    )
    ap.add_argument("--out", type=Path, default=c.OUT / "de4")
    a = ap.parse_args(argv)
    if a.combine:
        combine(a.out)
    elif a.parts:
        run(a.parts, a.out, a.case)
    else:
        ap.error("give --parts DIR [DIR ...] or --combine")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
