"""D-E4 analysis (exploratory, report-only; Amendment 19 (c), plan section 3.2): TOPAS phase space -> per-class summary.

Input: the TOPAS binary phase space of an EM-free ``g4h-elastic_HP`` thin-target case (``ps.header`` + ``ps.phsp``; protons
leaving the 500 mm vacuum sphere), one or several parts of one case (concatenated). The record layout is parsed from
``ps.header`` (no hard-coded offsets; the byte count per record is asserted to equal the sum of the field sizes and the file size
to equal ``records x bytes``).

Analysis (plan 3.2):
- exit angle ``theta_lab = arccos(W)``, ``W = +-sqrt(1 - U^2 - V^2)`` (the sign flag of the header);
- kinetic-energy filter: primary records (parent ID 0) with ``E > E0 - 1e-4 MeV`` are the undeflected primaries and are removed
  (OpenTOPAS has no current-KE filter, so they are recorded and removed here); the removed fraction is reported;
- "has a secondary proton" of an event: the event has a proton record with **parent ID != 0** (the p-p partner; the number of events
  with more than one proton record is reported as a cross-check);
- classification of every deflected primary with :func:`de5_sampling.classify` (smaller kinematic-curve residual to p+H / p+O wins;
  a conflict track, both residuals < 1e-3, is p-p iff the event has a secondary proton, else p+A);
- p-p events (event with a secondary proton): the faster of the two protons is the "fast proton" (the primary of the ionmc
  convention), the slower the "slow/secondary proton"; the class ``H-1`` signatures use the fast proton, ``mean_t_other_mev`` the
  slow proton's energy; p+A (``O-16``): ``mean_t_other_mev`` = ``E0 - T1`` (the recoil, energy conservation);
- per-class rates per primary ``N_class / N_primaries`` with the Poisson standard error;
- the signature values with the same constants and functions as :mod:`de5_sampling` (thresholds ``theta_w`` from the Bethe range of
  water, ``theta_med``, tail moment ``T`` over 3 theta_med, the descriptors ``shape_descriptors`` on the CM angle reconstructed
  from ``theta_lab`` and the energy with ``nuclear_v10.lab_to_cm``; recorded, not interpreted).

``geant4_summary`` has the shape ``{str(case): {class: {signature: value}}}`` that
``de5_sampling.signature_ratios(table_doc, geant4_summary)`` takes. Parts of one case are concatenated by passing all of them
to :func:`analyse` (counts add, primaries add).

Run: ``python -m validation.scripts.elastic_attribution.de4_parse --case 3 --ps DIR [DIR ...]``.
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
GIT_ALLOW = (
    "case", "material", "energy_mev", "length_cm", "primaries", "records", "secondary_definition", "pp_events_all",
    "pp_events_after_cm_cut", "pp_cm_cut_deg", "pp_events_removed_by_cm_cut", "theta_2mm_deg", "theta_3mm_deg",
    "ke_filter/*", "classification/*", "secondary_on_p+H_curve/*", "geant4_summary/*/*/*",
    *(f"{_CL}/{k}" for k in (
        "events", "rate_per_primary", "rate_per_primary_se", "theta_lab_median_deg", "P_theta_lab_gt_theta_2mm",
        "P_theta_lab_gt_theta_3mm", "T_tail_moment", "theta_lab_mean_1_minus_cos", "mean_t_other_mev", "first_minimum_deg",
        "first_minimum_3sigma_deg", "qslope", "qslope_bins", "counts_in_first_minimum_bin", "slow_proton_theta_lab_median_deg")),
    f"{_CL}/helm/*",
)  # fmt: skip
"""Git allowlist of a D-E4 summary (syntax: :mod:`common`): counts, per-primary rates, angle statistics, exceedances. Not kept:
histograms, the part directories (host paths)."""
GIT_FILES = {"de4*.json": (GIT_ALLOW, ())}


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
    m = re.search(
        r"^d:So/Beam/BeamEnergy\s*=\s*([0-9.eE+-]+)\s*MeV",
        (Path(directory) / "input.txt").read_text(),
        re.M,
    )
    assert m, "no BeamEnergy in input.txt"
    return float(m.group(1))


def analyse(
    dirs: list[Path], case: int = 3, material: str = "water", length_cm: float = 1.0
) -> dict[str, Any]:
    from ionmc.nuclear.bgg import M_PROTON_MEV as MP

    nv10 = c.import_transport_module("nuclear_v10")
    tab = c.load_table()
    water = c.Water(tab)
    names = tuple(tab.target_names)
    e0s = {beam_energy_mev(d) for d in dirs}
    assert len(e0s) == 1, f"parts differ in beam energy: {e0s}"
    e0 = e0s.pop()
    n_prim, recs = 0, []
    for k, d in enumerate(dirs):
        r, h = read_phsp(d)
        assert (r["pdg"] == PROTON).all(), "non-proton record"
        n_prim += h["histories"]
        r = r.copy()
        r["event"] = r["event"] + k * (1 << 24)  # keep events of different parts apart
        recs.append(r)
    r = np.concatenate(recs)
    w = np.sqrt(
        np.maximum(1.0 - r["U"].astype(float) ** 2 - r["V"].astype(float) ** 2, 0.0)
    ) * np.where(r["w_negative"] == 1, -1.0, 1.0)
    cos = w
    e = r["E"].astype(float)
    prim = r["parent"] == 0
    sec = ~prim
    undeflected = prim & (e > e0 - 1.0e-4)
    ev_has_sec = np.isin(r["event"], r["event"][sec])
    ev_multi = np.isin(
        r["event"],
        *[
            np.unique(r["event"], return_counts=True)[0][
                np.unique(r["event"], return_counts=True)[1] > 1
            ]
        ],
    )
    mh = float(tab.arrays["target_mass_mev"][names.index("H-1")])
    ma = float(tab.arrays["target_mass_mev"][names.index("O-16")])
    # deflected primaries
    dp = prim & ~undeflected
    is_pp, conflict, off = de5.classify(cos[dp], e[dp], e0, mh, ma, ev_has_sec[dp])
    # p-p events: every event with a secondary proton; fast / slow of the two protons (primary + first secondary)
    pp_events = np.unique(r["event"][sec])
    fast_e, fast_cos, slow_e, slow_cos = [], [], [], []
    order = np.argsort(r["event"], kind="stable")
    ev_sorted = r["event"][order]
    for ev in pp_events:
        lo, hi = np.searchsorted(ev_sorted, [ev, ev + 1])
        idx = order[lo:hi]
        idx = idx[np.argsort(-e[idx])]
        fast_e.append(e[idx[0]])
        fast_cos.append(cos[idx[0]])
        slow_e.append(e[idx[1]])
        slow_cos.append(cos[idx[1]])
    fast_e, fast_cos, slow_e, slow_cos = map(np.asarray, (fast_e, fast_cos, slow_e, slow_cos))
    # 16.26 deg CM cut on p-p (the table density stops there; applied to Geant4 too): equal masses, the fast proton has
    # tan(theta_lab) = tan(theta_cm / 2) / gamma_cm, so theta_cm = 2 arctan(gamma_cm tan(theta_lab))
    gam = nv10.two_body_frame(e0, MP)["gamma_cm"]
    th_cm_fast = 2.0 * np.degrees(np.arctan(gam * np.tan(np.arccos(np.clip(fast_cos, -1.0, 1.0)))))
    keep_cm = th_cm_fast >= de5.PP_CUT_DEG
    n_pp_all = int(fast_e.size)
    fast_e, fast_cos, slow_e, slow_cos = (
        fast_e[keep_cm],
        fast_cos[keep_cm],
        slow_e[keep_cm],
        slow_cos[keep_cm],
    )
    pa_mask = ~is_pp & ~ev_has_sec[dp]
    pa_e, pa_cos = e[dp][pa_mask], cos[dp][pa_mask]
    r_w = float(water.range_mm(e0))
    th_w = {x: math.degrees(math.acos(1.0 - x / r_w)) for x in (2.0, 3.0)}
    pcm = float(
        np.interp(e0, tab.arrays["grid_e_mev"], tab.arrays["p_cm_mev"][names.index("O-16")])
    )

    def cls_doc(
        t1: np.ndarray, cs: np.ndarray, t_other: np.ndarray, count: int, cm_mass: float | None
    ) -> dict[str, Any]:
        th = np.degrees(np.arccos(np.clip(cs, -1.0, 1.0)))
        med = float(np.median(th)) if th.size else float("nan")
        doc: dict[str, Any] = {
            "events": int(count), "rate_per_primary": count / n_prim, "rate_per_primary_se": math.sqrt(max(count, 1)) / n_prim,
            "theta_lab_median_deg": med, "P_theta_lab_gt_theta_2mm": float(np.mean(th > th_w[2.0])) if th.size else float("nan"),
            "P_theta_lab_gt_theta_3mm": float(np.mean(th > th_w[3.0])) if th.size else float("nan"),
            "T_tail_moment": float(np.mean((1.0 - cs) * (th > 3.0 * med))) if th.size else float("nan"),
            "theta_lab_mean_1_minus_cos": float(np.mean(1.0 - cs)) if th.size else float("nan"),
            "mean_t_other_mev": float(np.mean(t_other)) if th.size else float("nan"),
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

    classes = {
        "O-16": cls_doc(pa_e, pa_cos, e0 - pa_e, int(pa_mask.sum()), ma),
        "H-1": cls_doc(fast_e, fast_cos, slow_e, int(fast_e.size), None),
    }
    classes["H-1"]["slow_proton_theta_lab_median_deg"] = (
        float(np.median(np.degrees(np.arccos(np.clip(slow_cos, -1, 1)))))
        if slow_cos.size
        else float("nan")
    )
    cl = {
        "deflected_primaries": int(dp.sum()), "assigned_pp": int(is_pp.sum()), "assigned_pA": int((~is_pp).sum()),
        "conflict_tracks_raw": int(conflict.sum()), "conflict_fraction_raw": float(conflict.mean()) if conflict.size else float("nan"),
        "off_curve_tracks": int(off.sum()), "off_curve_fraction": float(off.mean()) if off.size else float("nan"),
        "conflict_assigned_pp_by_secondary": int((conflict & is_pp).sum()),
        "pA_events_with_secondary_proton": int((~is_pp & ev_has_sec[dp]).sum()),
        "pp_events_primary_not_assigned_pp": int((ev_has_sec[dp] & ~is_pp).sum()),
        "curve_crossover_deg": de5.curve_crossover_deg(e0, mh, ma),
    }  # fmt: skip
    sec_cs, sec_e = cos[sec], e[sec]
    sec_res = np.abs(sec_e - de5.kin_curve_t(e0, mh, sec_cs)) / np.maximum(sec_e, 1e-9)
    doc = {
        "banner": c.BANNER, "case": case, "material": material, "energy_mev": e0, "length_cm": length_cm,
        "primaries": n_prim, "parts": [str(d) for d in dirs], "records": int(r.size),
        "ke_filter": {"threshold_mev": e0 - 1.0e-4, "removed_primaries": int(undeflected.sum()), "removed_fraction": float(undeflected.sum() / n_prim),
                      "primary_records": int(prim.sum()), "secondary_records_parent_ne_0": int(sec.sum())},  # fmt: skip
        "secondary_definition": "event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = "
                                f"{int(np.unique(r['event'][ev_multi]).size)})",
        "pp_events_all": n_pp_all, "pp_events_after_cm_cut": int(fast_e.size), "pp_cm_cut_deg": de5.PP_CUT_DEG,
        "pp_events_removed_by_cm_cut": n_pp_all - int(fast_e.size), "classification": cl, "classes": classes,
        "secondary_on_p+H_curve": {"max_rel_residual": float(sec_res.max()) if sec_res.size else float("nan"),
                                   "fraction_below_1e-3": float(np.mean(sec_res < 1e-3)) if sec_res.size else float("nan")},
        "theta_2mm_deg": th_w[2.0], "theta_3mm_deg": th_w[3.0],
    }  # fmt: skip
    doc["geant4_summary"] = {
        str(case): {k: {s: v.get(s) for s in de5.SIGNATURES} for k, v in classes.items()}
    }
    return doc


def render_report(doc: dict[str, Any]) -> str:
    k, cl = doc["ke_filter"], doc["classification"]
    out = ["", f"Case {doc['case']} ({doc['material']}, {doc['energy_mev']:g} MeV, {doc['primaries']} primaries, {doc['records']} proton records).", "",
           f"KE filter (E > E0 - 1e-4 MeV): {k['removed_primaries']} primaries removed (fraction {k['removed_fraction']:.5f}); "
           f"{k['primary_records']} primary records, {k['secondary_records_parent_ne_0']} secondary proton records (parent != 0). "
           f"Secondary definition: {doc['secondary_definition']}.", "",
           "| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |", "|---|---|---|---|---|---|---|---|"]  # fmt: skip
    for name, v in doc["classes"].items():
        out.append(
            f"| {name} | {v['events']} | {v['rate_per_primary']:.4g} +- {v['rate_per_primary_se']:.2g} | {v['theta_lab_median_deg']:.2f} | "
            f"{v['P_theta_lab_gt_theta_2mm']:.4f} | {v['P_theta_lab_gt_theta_3mm']:.4f} | {v['T_tail_moment']:.4f} | {v['mean_t_other_mev']:.3f} |"
        )
    out += ["", f"p-p events: {doc['pp_events_all']} with a secondary proton; {doc['pp_events_removed_by_cm_cut']} removed by the {doc['pp_cm_cut_deg']} deg CM cut "
            f"(theta_cm from the fast proton), {doc['pp_events_after_cm_cut']} kept (class H-1).", "", f"Classification: {cl['deflected_primaries']} deflected primaries -> {cl['assigned_pA']} p+A, {cl['assigned_pp']} p-p; raw conflict tracks "
            f"{cl['conflict_tracks_raw']} ({cl['conflict_fraction_raw']:.3e}; crossover {cl['curve_crossover_deg']:.2f} deg), of which {cl['conflict_assigned_pp_by_secondary']} "
            f"p-p by the secondary rule; off-curve tracks {cl['off_curve_tracks']} ({cl['off_curve_fraction']:.3e}). Secondary protons on the p+H curve: "
            f"max residual {doc['secondary_on_p+H_curve']['max_rel_residual']:.2e}, fraction < 1e-3 {doc['secondary_on_p+H_curve']['fraction_below_1e-3']:.4f}.", "",
           "The first-minimum, q-slope and Helm descriptors are recorded, not interpreted."]  # fmt: skip
    return "\n".join(out) + "\n"


def run(case: int, dirs: list[Path], outdir: Path, name: str) -> dict[str, Any]:
    doc = analyse(dirs, case)
    inputs = [
        p for d in dirs for p in (d / "ps.header", d / "ps.phsp", d / "input.txt") if p.is_file()
    ]
    cmd = (
        f"python -m validation.scripts.elastic_attribution.de4_parse --case {case} --ps "
        + " ".join(str(d) for d in dirs)
    )
    c.write_outputs(
        outdir,
        f"{name}.json",
        doc,
        f"{name}.md",
        c.header("D-E4 analysis: TOPAS phase space -> per-class summary", cmd, inputs),
        render_report,
        GIT_ALLOW,
    )
    return doc


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--case", type=int, required=True)
    ap.add_argument(
        "--ps",
        type=Path,
        nargs="+",
        required=True,
        help="work directories holding ps.header / ps.phsp / input.txt",
    )
    ap.add_argument("--out", type=Path, default=c.OUT / "de4")
    ap.add_argument("--name", default="de4")
    a = ap.parse_args(argv)
    run(a.case, a.ps, a.out, a.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
