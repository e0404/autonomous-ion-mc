"""V6-F sizing table (revision 2, after review.md) from the C8 rehearsal partials. Planning evidence only;
never a verdict.

Inputs (read-only):
  * the pinned arXiv:1409.1938 e-print (tarball sha256 92cdaffb..., member halo.tex sha256 e5a270e7...)
    from the worktree data cache. ``tbl:dmlg`` is parsed with the repository parser
    ``gottschalk_v6.parse_dmlg``; ``tbl:parms`` with ``parse_parms``, which asserts the frozen p1, p2, p5-p7.
  * the three 1e6-history rehearsal partials (rehearsal seed base 20505000, head c847cc75). The seeds
    actually used are base + 1000 r_index + k:
      v6-s0               RUN-20261011T030353Z-0b25db97  seed 20517000 (r_index 12)
      v6-sens-nominal-s0  RUN-20261011T030829Z-537a9857  seed 20525000 (r_index 20)
      v6-sens-st-hi-s0    RUN-20261011T031141Z-88a536e0  seed 20525000 (r_index 20, sigma_theta x 1.2)

Blindness rule (enforced by construction):
  * Measured table. Used only for (a) the measured fraction of the r = 0 maximum, (b) the MD-fit aura
    share D_au/D_meas, where D_au comes from tbl:parms p23-p25 (no model output), and (c) the measured
    denominator of the V6-sens-F screen precision (sigma of the paired DIFFERENCE / D_meas).
  * Model partials. Used only for precision quantities: the relative SE of the gating sample v6-s0,
    non-empty batch counts, and the paired SE of the difference D_var - D_nom.
  * Never computed: a model mean on its own, a model/measured ratio, a residual, or a sensitivity
    ratio D_var/D_nom.
  * The ratio-SE helper ``gottschalk_v6.ratio_paired`` is no longer called. Revision 1 called it in
    memory and discarded the ratio.
  * The revision-1 column of model-normalised ratio SE was removed. Placed next to the
    measured-normalised SE of this revision it would reveal D_nom/D_meas, so the revision-1 table must
    never be combined with this one.
  * Remaining weak leak, disclosed: rel_se_1e6 (model-normalised, sample v6-s0) and
    sens_sigma_diff_over_meas_1e6 (measured-normalised, a different sample) are both per point. Their
    quotient involves D_model/D_meas times an unknown, point-dependent factor sigma_diff/se. That factor
    depends on the CRN correlation and is not a residual.

Outputs: se_table.json next to this script. Usage: python3 -I se_table.py <worktree>
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
import json
import math
import sys
import tarfile
from pathlib import Path

import numpy as np

TARBALL_SHA = "92cdaffbc4237e5178e8a83b3ceee6e15f708c2f4ebbd16aff63f249dc990ac2"
HALO_SHA = "e5a270e79c5a6d0e91f6c06936892a34d691c22b72f28ebfe6f0af7f8a5023d6"
RUNS = {
    "v6-s0": ("RUN-20261011T030353Z-0b25db97", "v6-s0.json",
              "615fcf31480c3adc24f871c1eee2f9ec666da1b8eece1df2db7930838fc462a3", 20517000),
    "sens-nominal": ("RUN-20261011T030829Z-537a9857", "v6-sens-nominal-s0.json",
                     "b6b16abb464b530e74d7c59070a2da2e9c86d4fbea037fb28a49160159ecc3f6", 20525000),
    "sens-st-hi": ("RUN-20261011T031141Z-88a536e0", "v6-sens-st-hi-s0.json",
                   "48f0f9f4c3d72b26161361f431029741e6315eb3e3f2d9e44bdda9db9d8c93c7", 20525000),
}
SHARD_N = 1_000_000
SEC_PER_SHARD_CPU = 184.0  # Amendment 14a (a)2: 184 s per 1e6 on warp-cpu f64 (5437 hist/s)
SE_HALO, SE_R0 = 0.03, 0.005  # frozen sizing targets of Amendment 14 (g)
SIGMA_SCREEN_F = 0.0375  # V6-sens-F sizing: 2 sigma'/D_meas <= 0.075, half of the 0.15 screen tolerance
SAFETY = 1.25  # SE-of-SE inflation for sizing the gating sample (20 batches)
Z_F_CM = 21.6  # deepest data depth <= R0 + sigma_sem = 21.614 cm (design 4)
HALO_FLOOR = 1e-4  # frozen sizing-subset floor (Amendment 14 (g), 17 (c)6)
AURA_MAX = 0.05  # design D7: points with D_au/D_meas > 0.05 excluded
N_F, K_F = 16, 8  # declared shards (design 5)
FLOORS = (1e-1, 1e-2, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 0.0)


def load_module(worktree: Path):
    path = worktree / "validation" / "scripts" / "reference" / "gottschalk_v6.py"
    spec = importlib.util.spec_from_file_location("gottschalk_v6", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gottschalk_v6"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_tex(worktree: Path) -> str:
    raw = (worktree / ".ionmc-cache" / "ionmc-data" / "objects" / TARBALL_SHA).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == TARBALL_SHA, "tarball hash mismatch"
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tar:
        tex = tar.extractfile("halo.tex").read()
    assert hashlib.sha256(tex).hexdigest() == HALO_SHA, "halo.tex hash mismatch"
    return tex.decode("utf-8", errors="replace")


def load_partial(worktree: Path, key: str) -> np.ndarray:
    run, name, content, seed = RUNS[key]
    p = worktree / "validation" / "generated" / run / "transport" / "reh-c847cc7-v6" / "samples" / name
    d = json.loads(p.read_text())
    assert d["content_sha256"] == content and d["valid"] and d["n"] == SHARD_N and d["n_batches"] == 20
    assert d["seed"] == seed and d["seed_base"] == 20505000
    return np.asarray(d["dose_batches"], dtype=np.float64)  # [20, 46, 10], MeV/g per primary


def rel_se_only(b: np.ndarray) -> np.ndarray:
    """Relative SE of the batch mean; the mean itself is not returned."""
    m = b.mean(axis=0)
    se = b.std(axis=0, ddof=1) / math.sqrt(b.shape[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(m > 0, se / m, np.inf)


def paired_diff_se(var: np.ndarray, nom: np.ndarray) -> np.ndarray:
    """Paired batch SE of mean(var) - mean(nom) (absolute, MeV/g/p); no ratio and no mean returned."""
    d = var - nom
    return d.std(axis=0, ddof=1) / math.sqrt(d.shape[0])


def aura_md(r_cm: float, p: dict[str, float]) -> float:
    """Gottschalk MD-fit aura dose (halo.tex l.784, 797-798): p23 E_2D(r_n, r_au), r_au = p24 (1 + p25 r_n),
    r_n = r / (R0/3), E_2D(r, r_e) = exp(-r/r_e) / (2 pi r_e^2); z-independent; MeV/g/p."""
    rn = r_cm / (p["p1"] / 3.0)
    rau = p["p24"] * (1.0 + p["p25"] * rn)
    return p["p23"] * math.exp(-rn / rau) / (2.0 * math.pi * rau * rau)


def shards_needed(x: float, target: float, safety: float = 1.0) -> int | None:
    if not math.isfinite(x):
        return None
    return max(1, math.ceil((safety * x / target) ** 2))


def sign_depths(sel: list[dict]) -> int:
    return len({p["z_cm"] for p in sel if p["r_cm"] <= 3} & {p["z_cm"] for p in sel if p["r_cm"] >= 4})


def main() -> None:
    worktree = Path(sys.argv[1]).resolve()
    out = Path(__file__).resolve().parent / "se_table.json"
    g = load_module(worktree)
    tex = load_tex(worktree)
    z, r, logd, n_missing = g.parse_dmlg(tex)
    parms = g.parse_parms(tex)
    dmeas = np.power(10.0, logd)
    r0 = int(np.where(r == 0)[0][0])
    dmax = float(np.nanmax(dmeas[:, r0]))
    jmax = int(np.nanargmax(dmeas[:, r0]))

    b = load_partial(worktree, "v6-s0")
    rse = rel_se_only(b)
    nonempty = (b > 0).sum(axis=0)
    bn, bv = load_partial(worktree, "sens-nominal"), load_partial(worktree, "sens-st-hi")
    sdiff = paired_diff_se(bv, bn)
    ne_sens = np.minimum((bn > 0).sum(axis=0), (bv > 0).sum(axis=0))

    points = []
    for j in range(len(z)):
        for i in range(len(r)):
            missing = not np.isfinite(dmeas[j, i])
            rs = float(rse[j, i])
            dau = aura_md(float(r[i]), parms)
            sdm = None if missing else float(sdiff[j, i] / dmeas[j, i])
            points.append({
                "z_cm": float(z[j]), "r_cm": float(r[i]), "group": "r0" if r[i] == 0 else "halo",
                "missing": missing,
                "measured_fraction_of_r0_max": None if missing else float(dmeas[j, i] / dmax),
                "aura_md_MeV_per_g": dau,
                "aura_share_of_measured": None if missing else float(dau / dmeas[j, i]),
                "rel_se_1e6": rs if math.isfinite(rs) else None,
                "batches_nonempty_of_20": int(nonempty[j, i]),
                "shards_for_3pct": shards_needed(rs, SE_HALO),
                "shards_for_0p5pct": shards_needed(rs, SE_R0),
                "shards_for_3pct_safety": shards_needed(rs, SE_HALO, SAFETY),
                "shards_for_0p5pct_safety": shards_needed(rs, SE_R0, SAFETY),
                "sens_sigma_diff_over_meas_1e6_st_hi": sdm,
                "sens_batches_nonempty_min_of_20": int(ne_sens[j, i]),
                "sens_k_for_2sigma_0p075": None if sdm is None else shards_needed(sdm, SIGMA_SCREEN_F),
            })

    def worst(sel, key):
        vals = [p[key] for p in sel]
        if any(v is None for v in vals):
            return {"max": None, "n_undefined": sum(v is None for v in vals)}
        return {"max": max(vals) if vals else 0, "n_undefined": 0}

    def argworst(sel, key):
        fin = [p for p in sel if p[key] is not None]
        if not fin:
            return None
        w = max(fin, key=lambda p: p[key])
        return {"z_cm": w["z_cm"], "r_cm": w["r_cm"], key: w[key]}

    present = [p for p in points if not p["missing"]]
    halo = [p for p in present if p["group"] == "halo"]
    r0pts = [p for p in present if p["group"] == "r0"]
    floors = []
    for f in FLOORS:
        hs = [p for p in halo if p["measured_fraction_of_r0_max"] >= f]
        floors.append({"floor": f, "halo_points": len(hs), "halo_shards_3pct": worst(hs, "shards_for_3pct"),
                       "halo_worst_point_3pct": argworst(hs, "shards_for_3pct"),
                       "sens_k_2sigma_0p075": worst(hs, "sens_k_for_2sigma_0p075"), "sign_test_depths": sign_depths(hs)})
    zcuts = []
    for zc in (19.0, 20.0, 21.0, 21.4, 21.6, 21.8, 22.0, 22.2):
        rs_ = [p for p in r0pts if p["z_cm"] <= zc + 1e-9]
        zcuts.append({"z_max_cm": zc, "r0_points": len(rs_), "r0_shards_0p5pct": worst(rs_, "shards_for_0p5pct"),
                      "r0_shards_0p5pct_safety": worst(rs_, "shards_for_0p5pct_safety")})

    # ---- V6-F sets (design.md, revision 2)
    f_r0 = [p for p in r0pts if p["z_cm"] <= Z_F_CM + 1e-9]
    h_f = [p for p in halo if p["z_cm"] <= Z_F_CM + 1e-9 and p["measured_fraction_of_r0_max"] >= HALO_FLOOR]
    h_fa = [p for p in h_f if p["aura_share_of_measured"] <= AURA_MAX]  # H_F after D7
    by_r = lambda sel: {str(int(rr)): sum(1 for p in sel if p["r_cm"] == rr) for rr in (1, 2, 3, 4, 5, 6, 7, 8, 10)}
    aura_stats = {}
    for rr in (1, 2, 3, 4, 5, 6):
        s = sorted(p["aura_share_of_measured"] for p in h_f if p["r_cm"] == rr)
        aura_stats[str(rr)] = {"n": len(s), "median": s[len(s) // 2] if s else None, "max": s[-1] if s else None,
                               "n_above_0p05": sum(x > AURA_MAX for x in s), "n_above_0p10": sum(x > 0.10 for x in s)}
    sd1 = max(p["sens_sigma_diff_over_meas_1e6_st_hi"] for p in h_fa)
    k_rule = math.ceil((sd1 / SIGMA_SCREEN_F) ** 2)
    n_rule = max(worst(f_r0, "shards_for_0p5pct_safety")["max"], worst(h_fa, "shards_for_3pct_safety")["max"])
    v6f = {
        "z_F_cm": Z_F_CM, "halo_floor": HALO_FLOOR, "aura_max_share": AURA_MAX,
        "r0_points": len(f_r0),
        "H_F_before_aura_points": len(h_f), "H_F_before_aura_by_r": by_r(h_f), "H_F_before_aura_sign_depths": sign_depths(h_f),
        "aura_share_by_r_in_H_F_before_exclusion": aura_stats,
        "halo_points": len(h_fa), "halo_points_by_r": by_r(h_fa), "sign_test_depths": sign_depths(h_fa),
        "aura_excluded_points": len(h_f) - len(h_fa), "aura_excluded_by_r": by_r([p for p in h_f if p not in h_fa]),
        "r0_shards_0p5pct": worst(f_r0, "shards_for_0p5pct"), "r0_shards_0p5pct_safety": worst(f_r0, "shards_for_0p5pct_safety"),
        "halo_shards_3pct": worst(h_fa, "shards_for_3pct"), "halo_shards_3pct_safety": worst(h_fa, "shards_for_3pct_safety"),
        "halo_worst_3pct": argworst(h_fa, "shards_for_3pct"), "r0_worst_0p5pct": argworst(f_r0, "shards_for_0p5pct"),
        "N_F_rule_shards": n_rule, "N_F_declared_shards": N_F,
        "expected_worst_rel_se_at_N_F": {"halo": max(p["rel_se_1e6"] for p in h_fa) / math.sqrt(N_F),
                                         "r0": max(p["rel_se_1e6"] for p in f_r0) / math.sqrt(N_F)},
        "screen_sigma_diff_over_meas_1e6_max": sd1,
        "screen_worst": argworst(h_fa, "sens_sigma_diff_over_meas_1e6_st_hi"),
        "k_F_rule_2sigma_0p075": k_rule, "k_F_declared": K_F,
        "screen_2sigma_at_k_F_max": 2 * sd1 / math.sqrt(K_F),
        "screen_n_points_2sigma_gt_0p075_at_k_F": sum(2 * p["sens_sigma_diff_over_meas_1e6_st_hi"] / math.sqrt(K_F) > 0.075 for p in h_fa),
        "screen_k_rule_on_H_F_before_aura": math.ceil((max(p["sens_sigma_diff_over_meas_1e6_st_hi"] for p in h_f) / SIGMA_SCREEN_F) ** 2),
        "cpu_planning_seconds": {"gating": N_F * SEC_PER_SHARD_CPU, "sens": 7 * K_F * SEC_PER_SHARD_CPU,
                                 "total_h": (N_F + 7 * K_F) * SEC_PER_SHARD_CPU / 3600},
    }
    doc = {
        "revision": 2,
        "purpose": "V6-F sizing (planning evidence; rehearsal seeds 20517000 / 20525000 at base 20505000 are never evidence)",
        "blindness": __doc__.split("Blindness rule (enforced by construction):")[1].split("Outputs:")[0].strip(),
        "inputs": {k: {"run": v[0], "file": v[1], "content_sha256": v[2], "seed": v[3]} for k, v in RUNS.items()},
        "data": {"tarball_sha256": TARBALL_SHA, "halo_tex_sha256": HALO_SHA, "n_depths": len(z), "n_radii": len(r),
                 "n_missing": n_missing, "n_present": len(present), "n_r0_present": len(r0pts), "n_halo_present": len(halo),
                 "r0_max_measured_MeV_per_g": dmax, "r0_max_depth_cm": float(z[jmax]),
                 "parms_used": {k: parms[k] for k in ("p1", "p2", "p23", "p24", "p25")}},
        "constants": {"shard_histories": SHARD_N, "sec_per_shard_cpu_f64": SEC_PER_SHARD_CPU, "se_target_halo": SE_HALO,
                      "se_target_r0": SE_R0, "sigma_screen_F": SIGMA_SCREEN_F, "safety_factor_on_rse": SAFETY,
                      "scaling": "shards = ceil((safety * x_1e6 / target)^2), SE ~ N^-1/2; null = zero model dose in 1e6",
                      "screen_proxy": "paired SE of D_st-hi - D_nominal (CRN) over D_meas; the only rehearsed pair"},
        "v6f": v6f, "floors": floors, "r0_depth_cuts": zcuts, "points": points,
    }
    out.write_text(json.dumps(doc, indent=1) + "\n")
    print(json.dumps(v6f, indent=1))


if __name__ == "__main__":
    main()
