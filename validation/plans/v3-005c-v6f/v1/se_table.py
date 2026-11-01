"""QUARANTINED REVISION-1 SCRIPT (kept for provenance only; superseded by ../se_table.py, revision 2).
Its per-point precision columns and the model-normalised sigma_Delta were removed (review re-check R2-I1).
Never commit or archive any table it produced next to the revision-2 table.

V6-F sizing table from the C8 rehearsal partials (planning evidence only; never a verdict).

Inputs (read-only):
  * the pinned arXiv:1409.1938 e-print (sha256 92cdaffb...; member halo.tex sha256 e5a270e7...) from the
    worktree data cache, parsed with the repository parser ``gottschalk_v6.parse_dmlg``;
  * the three 1e6-history rehearsal partials at base 20505000 (head c847cc75):
      v6-s0               RUN-20261011T030353Z-0b25db97  (r_index 12, seed 20517000, nominal)
      v6-sens-nominal-s0  RUN-20261011T030829Z-537a9857  (r_index 20, seed 20525000, nominal)
      v6-sens-st-hi-s0    RUN-20261011T031141Z-88a536e0  (r_index 20, seed 20525000, sigma_theta x 1.2)

Blindness rule of this script (stated because the design must stay free of post-hoc tuning):
  * the measured table is used ONLY for the point's measured fraction of the r = 0 maximum (point
    selection by a measured-dose floor, the same use as the frozen 1e-4 sizing subset);
  * the model partials are used ONLY for precision quantities (relative SE, non-empty batch counts, the
    paired CRN ratio SE). No model mean dose, no model/measured ratio, no residual and no sensitivity
    ratio |D_var/D_nom - 1| is computed, printed or written.

Outputs: se_table.json next to this script.

Usage: python3 -I se_table.py <worktree>
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
              "615fcf31480c3adc24f871c1eee2f9ec666da1b8eece1df2db7930838fc462a3"),
    "sens-nominal": ("RUN-20261011T030829Z-537a9857", "v6-sens-nominal-s0.json",
                     "b6b16abb464b530e74d7c59070a2da2e9c86d4fbea037fb28a49160159ecc3f6"),
    "sens-st-hi": ("RUN-20261011T031141Z-88a536e0", "v6-sens-st-hi-s0.json",
                   "48f0f9f4c3d72b26161361f431029741e6315eb3e3f2d9e44bdda9db9d8c93c7"),
}
SHARD_N = 1_000_000
SEC_PER_SHARD_CPU = 184.0  # Amendment 14a (a)2: 184 s per 1e6 on warp-cpu f64 (5437 hist/s)
SE_HALO, SE_R0 = 0.03, 0.005  # the frozen sizing targets of Amendment 14 (g)
SIGMA_DELTA_FROZEN = 0.075  # 2 sigma_Delta < 0.15 alone (the Amendment 14a (a)3 sizing reading)
SIGMA_DELTA_F = 0.025  # V6-sens-F sizing target: 2 sigma_Delta = 0.05 leaves 0.10 of the 0.15 for |R - 1|
SAFETY = 1.25  # SE-of-SE inflation for sizing (20 batches; see design.md section 5)
FLOORS = (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 0.0)


def load_module(worktree: Path):
    path = worktree / "validation" / "scripts" / "reference" / "gottschalk_v6.py"
    spec = importlib.util.spec_from_file_location("gottschalk_v6", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gottschalk_v6"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_table(worktree: Path, g):
    raw = (worktree / ".ionmc-cache" / "ionmc-data" / "objects" / TARBALL_SHA).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == TARBALL_SHA, "tarball hash mismatch"
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tar:
        tex = tar.extractfile("halo.tex").read()
    assert hashlib.sha256(tex).hexdigest() == HALO_SHA, "halo.tex hash mismatch"
    return g.parse_dmlg(tex.decode("utf-8", errors="replace"))


def load_partial(worktree: Path, key: str) -> np.ndarray:
    run, name, content = RUNS[key]
    p = worktree / "validation" / "generated" / run / "transport" / "reh-c847cc7-v6" / "samples" / name
    d = json.loads(p.read_text())
    assert d["content_sha256"] == content and d["valid"] and d["n"] == SHARD_N and d["n_batches"] == 20
    return np.asarray(d["dose_batches"], dtype=np.float64)  # [20, 46, 10], MeV/g per primary


def shards_needed(rse: float, target: float, safety: float = 1.0) -> int | None:
    if not math.isfinite(rse):
        return None
    return max(1, math.ceil((safety * rse / target) ** 2))


def main() -> None:
    worktree = Path(sys.argv[1]).resolve()
    out = Path(__file__).resolve().parent / "se_table.json"
    g = load_module(worktree)
    z, r, logd, n_missing = load_table(worktree, g)
    dmeas = np.power(10.0, logd)
    r0 = int(np.where(r == 0)[0][0])
    dmax = float(np.nanmax(dmeas[:, r0]))
    jmax = int(np.nanargmax(dmeas[:, r0]))

    b = load_partial(worktree, "v6-s0")
    mean, se = g.batch_mean_se(b)
    with np.errstate(divide="ignore", invalid="ignore"):
        rse = np.where(mean > 0, se / mean, np.inf)
    nonempty = (b > 0).sum(axis=0)

    # QUARANTINE (revision 2, R2-I1): the model-normalised paired sigma_Delta (ratio_paired) is no longer
    # computed; next to the revision-2 measured-normalised sigma it would give a model/measured map.
    sd = np.full(rse.shape, np.inf)
    ne_sens = np.zeros(rse.shape, dtype=int)

    points = []
    for j in range(len(z)):
        for i in range(len(r)):
            missing = not np.isfinite(dmeas[j, i])
            frac = None if missing else float(dmeas[j, i] / dmax)
            rs = float(rse[j, i])
            sdi = float(sd[j, i])
            points.append({
                "z_cm": float(z[j]), "r_cm": float(r[i]), "group": "r0" if r[i] == 0 else "halo",
                "missing": missing, "measured_fraction_of_r0_max": frac,
                "rel_se_1e6": rs if math.isfinite(rs) else None,
                "batches_nonempty_of_20": int(nonempty[j, i]),
                "shards_for_3pct": shards_needed(rs, SE_HALO),
                "shards_for_0p5pct": shards_needed(rs, SE_R0),
                "shards_for_3pct_safety": shards_needed(rs, SE_HALO, SAFETY),
                "shards_for_0p5pct_safety": shards_needed(rs, SE_R0, SAFETY),
                "sens_sigma_delta_1e6_st_hi": sdi if math.isfinite(sdi) else None,
                "sens_batches_nonempty_min_of_20": int(ne_sens[j, i]),
                "sens_k_for_sigma_0p075": shards_needed(sdi, SIGMA_DELTA_FROZEN),
                "sens_k_for_sigma_0p025": shards_needed(sdi, SIGMA_DELTA_F),
                "sens_k_for_sigma_0p025_safety": shards_needed(sdi, SIGMA_DELTA_F, SAFETY),
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
        rs_ = [p for p in r0pts if p["measured_fraction_of_r0_max"] >= f]
        floors.append({
            "floor": f, "halo_points": len(hs), "r0_points": len(rs_),
            "halo_shards_3pct": worst(hs, "shards_for_3pct"),
            "halo_shards_3pct_safety": worst(hs, "shards_for_3pct_safety"),
            "halo_worst_point_3pct": argworst(hs, "shards_for_3pct"),
            "r0_shards_0p5pct": worst(rs_, "shards_for_0p5pct"),
            "r0_shards_0p5pct_safety": worst(rs_, "shards_for_0p5pct_safety"),
            "r0_worst_point_0p5pct": argworst(rs_, "shards_for_0p5pct"),
            "sens_k_0p075": worst(hs, "sens_k_for_sigma_0p075"),
            "sens_k_0p025": worst(hs, "sens_k_for_sigma_0p025"),
            "sens_k_0p025_safety": worst(hs, "sens_k_for_sigma_0p025_safety"),
            "sens_worst_point_0p025": argworst(hs, "sens_k_for_sigma_0p025"),
            "halo_min_batches_nonempty": min((p["batches_nonempty_of_20"] for p in hs), default=None),
            "halo_depths_with_both_groups": len({p["z_cm"] for p in hs if p["r_cm"] <= 3}
                                                & {p["z_cm"] for p in hs if p["r_cm"] >= 4}),
        })
    # r = 0 depth cut (z <= zc), the alternative r = 0 selection
    zcuts = []
    for zc in (19.0, 20.0, 21.0, 21.4, 21.6, 21.8, 22.0, 22.2):
        rs_ = [p for p in r0pts if p["z_cm"] <= zc + 1e-9]
        zcuts.append({"z_max_cm": zc, "r0_points": len(rs_), "r0_shards_0p5pct": worst(rs_, "shards_for_0p5pct"),
                      "r0_shards_0p5pct_safety": worst(rs_, "shards_for_0p5pct_safety"),
                      "min_measured_fraction": min(p["measured_fraction_of_r0_max"] for p in rs_)})

    # ---- V6-F (design.md): z <= Z_F at every radius; halo additionally measured >= 1e-4 of the r = 0 max
    zf = 21.6
    f_r0 = [p for p in r0pts if p["z_cm"] <= zf + 1e-9]
    f_h = [p for p in halo if p["z_cm"] <= zf + 1e-9 and p["measured_fraction_of_r0_max"] >= 1e-4]
    sd1 = max(p["sens_sigma_delta_1e6_st_hi"] for p in f_h)
    n_f = max(worst(f_r0, "shards_for_0p5pct_safety")["max"], worst(f_h, "shards_for_3pct_safety")["max"])
    k_f = math.ceil((sd1 / 0.0375) ** 2)
    excl_r0 = [p for p in r0pts if p not in f_r0]
    excl_h = [p for p in halo if p not in f_h]
    v6f = {
        "z_F_cm": zf, "halo_floor": 1e-4,
        "r0_points": len(f_r0), "halo_points": len(f_h),
        "halo_points_by_r": {str(rr): sum(1 for p in f_h if p["r_cm"] == rr) for rr in (1, 2, 3, 4, 5, 6, 7, 8, 10)},
        "sign_test_depths": len({p["z_cm"] for p in f_h if p["r_cm"] <= 3} & {p["z_cm"] for p in f_h if p["r_cm"] >= 4}),
        "r0_shards_0p5pct": worst(f_r0, "shards_for_0p5pct"), "r0_shards_0p5pct_safety": worst(f_r0, "shards_for_0p5pct_safety"),
        "halo_shards_3pct": worst(f_h, "shards_for_3pct"), "halo_shards_3pct_safety": worst(f_h, "shards_for_3pct_safety"),
        "halo_worst_3pct": argworst(f_h, "shards_for_3pct"), "r0_worst_0p5pct": argworst(f_r0, "shards_for_0p5pct"),
        "N_F_rule_shards": n_f, "N_F_declared_shards": 16,
        "sens_sigma_delta_1e6_max": sd1, "sens_worst": argworst(f_h, "sens_k_for_sigma_0p025"),
        "k_F_rule_sigma_0p0375": k_f, "k_F_declared": 8,
        "sens_2sigma_at_k8_max": 2 * sd1 / math.sqrt(8),
        "excluded_r0_points": [p["z_cm"] for p in excl_r0],
        "excluded_halo_points": len(excl_h),
        "excluded_halo_beyond_zF": sum(1 for p in excl_h if p["z_cm"] > zf + 1e-9),
        "excluded_halo_below_floor_within_zF": sum(1 for p in excl_h if p["z_cm"] <= zf + 1e-9),
        "excluded_halo_zero_model_dose_1e6": sum(1 for p in excl_h if p["rel_se_1e6"] is None),
        "cpu_seconds": {"gating": 16 * SEC_PER_SHARD_CPU, "sens_7x8": 56 * SEC_PER_SHARD_CPU,
                        "total": 72 * SEC_PER_SHARD_CPU, "total_h": 72 * SEC_PER_SHARD_CPU / 3600},
        "cuda_seconds_planning_10x": 16 * SEC_PER_SHARD_CPU / 10, "cuda_seconds_worst_1x": 16 * SEC_PER_SHARD_CPU,
    }
    doc = {
        "v6f": v6f,
        "purpose": "V6-F sizing (planning evidence; rehearsal base 20505000 is never evidence)",
        "blindness": "measured table used only for measured_fraction_of_r0_max; model partials only for precision "
                     "quantities; no model mean, residual or sensitivity ratio computed",
        "inputs": {k: {"run": v[0], "file": v[1], "content_sha256": v[2]} for k, v in RUNS.items()},
        "data": {"tarball_sha256": TARBALL_SHA, "halo_tex_sha256": HALO_SHA, "n_depths": len(z),
                 "n_radii": len(r), "n_missing": n_missing, "n_present": len(present), "n_r0_present": len(r0pts),
                 "n_halo_present": len(halo), "r0_max_measured_MeV_per_g": dmax, "r0_max_depth_cm": float(z[jmax])},
        "constants": {"shard_histories": SHARD_N, "sec_per_shard_cpu_f64": SEC_PER_SHARD_CPU,
                      "se_target_halo": SE_HALO, "se_target_r0": SE_R0, "sigma_delta_frozen_reading": SIGMA_DELTA_FROZEN,
                      "sigma_delta_F": SIGMA_DELTA_F, "safety_factor_on_rse": SAFETY,
                      "scaling": "shards = ceil((safety * rse_1e6 / target)^2), SE ~ N^-1/2; null = zero model dose in "
                                 "1e6 (undefined)",
                      "sens_proxy": "sigma_Delta from the only rehearsed variant pair (st-hi vs nominal, CRN)"},
        "floors": floors, "r0_depth_cuts": zcuts, "points": points,
    }
    keep = ("z_cm", "r_cm", "group", "missing", "measured_fraction_of_r0_max")
    doc["points"] = [{k: p[k] for k in keep} for p in points]  # R2-I1: no per-point precision column
    doc["quarantine"] = "revision-1 script, neutralised in revision 2: no per-point precision columns, no sigma_Delta"
    out.write_text(json.dumps(doc, indent=1) + "\n")
    # short console summary (precision quantities only)
    print(f"present {len(present)} (r0 {len(r0pts)}, halo {len(halo)}), missing {n_missing}; r0 max at z={z[jmax]}")
    for f in floors:
        print(f"floor {f['floor']:>7g}: halo {f['halo_points']:3d} r0 {f['r0_points']:2d} | halo3% "
              f"{f['halo_shards_3pct']} saf {f['halo_shards_3pct_safety']['max']} | r0 0.5% {f['r0_shards_0p5pct']} "
              f"saf {f['r0_shards_0p5pct_safety']['max']} | sens.075 {f['sens_k_0p075']['max']} sens.025 "
              f"{f['sens_k_0p025']} saf {f['sens_k_0p025_safety']['max']} | minNE {f['halo_min_batches_nonempty']} "
              f"| signdepths {f['halo_depths_with_both_groups']} | worst {f['halo_worst_point_3pct']}")
    for c in zcuts:
        print("zcut", c)
    print("V6-F", json.dumps(v6f, indent=1))


if __name__ == "__main__":
    main()
