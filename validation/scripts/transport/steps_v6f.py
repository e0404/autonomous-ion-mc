"""Row V6-F / V6-sens-F (Gottschalk 2015 halo, reduced scope): producers, the V6-F timing step and the comparator
(plan Amendment 18 (b), Amendment 20 (e)-(g); full text ``validation/plans/v3-005c-v6f-followup.md`` revision 2).

Usage::

    python validation/scripts/transport/steps_v6f.py <step> [--shard K] [--variant NAME] [--scale F] ...

Steps: ``v6f-shard`` (warp-cpu float64, r_index 7, k = 0..15), ``v6f-sens-shard`` (warp-cpu float64, r_index 8,
k = 0..7, ``--variant`` nominal or one of the six beam variants, common random numbers), ``v6f-cuda-shard``
(warp-cuda float32, r_index 9, k = 0..15), ``v6f-compare`` (consumer) and ``v6f-timing`` (one 1e6-history CUDA
float32 shard of the V6-F beam at the rehearsal seed 20514000: timing only, never evidence).

The producers are the existing ``steps_v6._produce`` (same configuration, scoring and sphere average; only the
row label, the step name, the seeds and the kinds differ). Seeds are ``B + 1000 * r_index + k`` with the V6-F base
``B = 20551004`` (Amendment 20 (f)2): 20558004-20558019, 20559004-20559011, 20560004-20560019. The configuration
(elastic on through the D-1 table, E0, beam, sphere average) is the one of ``steps_v6.v6_config``; V6-F certifies
that configuration only. Nothing here is a V6 step: V6 / V6-sens stay recorded as failed (Amendment 18 (c)).

The comparator (``evaluate_v6f``, pure) follows v6f-followup items 3, 4, 6, 7 and 9:

* z_F = deepest data depth <= R0 + sigma_sem = p1 (1 + p2); P_F,0 = non-missing r = 0 points with z <= z_F, each
  gating at 2 %; the deeper r = 0 points are report-only (D1);
* H_F = non-missing r >= 1 cm points with z <= z_F, measured dose >= 1e-4 of the measured r = 0 maximum and
  aura share D_au / D_meas <= 0.05, D_au = p23 E_2D(r_n, r_au) (D2, D7);
* H_g = points of H_F with |D_var - D_nom| / D_meas + 2 sigma' / D_meas < 0.15 for all six beam variants (D3);
* coverage minimum (D8): H_g has >= 40 points, >= 3 radii and >= 12 sign-test depths after excluding ties; else the
  halo is "not established (beam-sensitive)", never a pass;
* halo on H_g: 15 % per point, >= 90 % of points, median |residual| <= 0.10, sign-test p > 0.01; r = 0: 2 %;
* the CUDA sample is judged on the same H_g; V6-F fails if the backends differ on r = 0, halo, overall or coverage
  (Amendment 17 (c)8); per-point |z| > 4 is listed.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

import steps as base
import steps_v5 as v5
import steps_v5b as v5b
import steps_v5c as v5c
import steps_v6 as s6
from ionmc.simulation import Simulation

V6F_BASE = 20551004  # V6-F row base (Amendment 20 (f)2)
V6F_R_CPU, V6F_R_SENS, V6F_R_CUDA = 7, 8, 9
V6F_CPU_SHARDS, V6F_SENS_SHARDS, V6F_CUDA_SHARDS = 16, 8, 16
V6F_TIMING_SEED = 20514000  # rehearsal base 20505000, r_index 9, k = 0 (Amendment 20 (f)6)
V6F_BEAM_VARIANTS = ("sx-lo", "sx-hi", "st-lo", "st-hi", "rho-lo", "rho-hi")
V6F_SENS_RUNS = ("nominal", *V6F_BEAM_VARIANTS)
V6F_AURA_MAX = 0.05
V6F_FLOOR = 1.0e-4  # of the measured r = 0 maximum
V6F_MIN_POINTS, V6F_MIN_RADII, V6F_MIN_DEPTHS = 40, 3, 12
V6F_SCREEN_TOL, V6F_SCREEN_NSIGMA = 0.15, 2.0
V6F_R0_TARGET_REL_SE, V6F_HALO_TARGET_REL_SE = 0.005, 0.03
V6F_OLDEST_REHEARSAL_MAX = 20525999  # largest rehearsal seed actually used (route 4.4)
V6F_ORIGINAL_V6_BASES = (20493004, 20494004, 20501004)  # never-run original V6 seeds, base + k
V6F_STEP_LIMIT_PLANNED_S = 1200.0
V6F_STEP_LIMIT_HARD_S = 1500.0
V6F_ROW = "V6-F"
NOT_ESTABLISHED = "not established (beam-sensitive)"


# -- names, seeds -------------------------------------------------------------------------------
def cpu_names() -> tuple[str, ...]:
    return tuple(f"v6f-s{k}.json" for k in range(V6F_CPU_SHARDS))


def sens_names() -> tuple[str, ...]:
    return tuple(f"v6f-sens-{v}-s{k}.json" for v in V6F_SENS_RUNS for k in range(V6F_SENS_SHARDS))


def cuda_names() -> tuple[str, ...]:
    return tuple(f"v6f-cuda-s{k}.json" for k in range(V6F_CUDA_SHARDS))


def seeds_of(base_seed: int) -> dict[str, list[int]]:
    """Seeds of a V6-F execution at base ``base_seed``: B + 1000 r_index + k."""
    return {"cpu": [base_seed + 1000 * V6F_R_CPU + k for k in range(V6F_CPU_SHARDS)],
            "sens": [base_seed + 1000 * V6F_R_SENS + k for k in range(V6F_SENS_SHARDS)],
            "cuda": [base_seed + 1000 * V6F_R_CUDA + k for k in range(V6F_CUDA_SHARDS)]}  # fmt: skip


EVIDENCE_BASES_USED = (20481004, 20491004)  # lv5c / hr5c qualification bases (dd976c5a): evidence seeds observed


def consumed_evidence_seeds() -> set[int]:
    """Every evidence seed used up to now (v6f-followup 7 check (a): "``consumed_seed_set()`` and every evidence
    seed used up to then"): ``steps_v5c.consumed_seed_set()`` plus the conservative envelope ``base + 1000 r + k``
    (r_index 0..20, k < 1000) of the lv5c and hr5c bases, whose seeds were observed at dd976c5a."""
    out = set(v5c.consumed_seed_set())
    for b in EVIDENCE_BASES_USED:
        out |= {b + 1000 * r + k for r in range(21) for k in range(1000)}
    return out


def check_producer_seed(seed: int, *, consumed: set[int] | None = None) -> None:
    """Fail closed before a shard is run: the seed must not be a consumed evidence seed (check (a), including the
    lv5c / hr5c envelopes), must exceed the largest rehearsal seed used (check (b); the V6-F seeds are
    B + 7000 + ... > 20525999) and must not be one of the never-run original V6 seeds (check (c))."""
    used = consumed_evidence_seeds() if consumed is None else consumed
    original = {b + k for b in V6F_ORIGINAL_V6_BASES for k in range(1000)}
    if seed in used:
        raise SystemExit(f"V6-F seed {seed} is a consumed evidence seed (check (a); fail closed)")
    if seed in original:
        raise SystemExit(f"V6-F seed {seed} is a never-run original V6 seed (check (c); fail closed)")
    if seed <= V6F_OLDEST_REHEARSAL_MAX:
        raise SystemExit(f"V6-F seed {seed} is not above the largest rehearsal seed used (check (b); fail closed)")


def _other_rows_envelope() -> set[int]:
    """Seeds of the other rows of the executing task (check (d)): every lv5d seed (base 20531004, r_index 0..20,
    k < 1000) and the hr5d rows that draw seeds (hr 20541004, r_index 1, 2, 3)."""
    out = {20531004 + 1000 * r + k for r in range(21) for k in range(1000)}
    out |= {20541004 + 1000 * r + k for r in (1, 2, 3) for k in range(1000)}
    return out


def check_v6f_seeds(
    cpu: list[int], sens: list[int], cuda: list[int], *, rehearsal_used_max: int = V6F_OLDEST_REHEARSAL_MAX,
    consumed: set[int] | None = None, other_rows: set[int] | None = None,
) -> dict[str, bool]:
    """Checks (a)-(e) of v6f-followup 7 for the three samples; fails closed (``SystemExit``) on any violation.
    (a) disjoint from ``consumed_seed_set()``; (b) every seed above the largest rehearsal seed used; (c) disjoint
    from the never-run original V6 seeds ``20493004 / 20494004 / 20501004 + k`` (k < 1000); (d) disjoint from the
    other rows of the executing task; (e) CPU, sensitivity and CUDA samples pairwise disjoint and each free of
    repeats."""
    consumed = consumed_evidence_seeds() if consumed is None else consumed
    other = _other_rows_envelope() if other_rows is None else other_rows
    original = {b + k for b in V6F_ORIGINAL_V6_BASES for k in range(1000)}
    groups = {"cpu": cpu, "sens": sens, "cuda": cuda}
    every = [s for g in groups.values() for s in g]
    bad: list[str] = []
    if any(len(set(g)) != len(g) for g in groups.values()):
        bad.append("(e) repeated seed inside a sample")
    for (n1, g1), (n2, g2) in ((("cpu", cpu), ("sens", sens)), (("cpu", cpu), ("cuda", cuda)),
                               (("sens", sens), ("cuda", cuda))):  # fmt: skip
        if set(g1) & set(g2):
            bad.append(f"(e) {n1} and {n2} samples share a seed")
    if set(every) & consumed:
        bad.append("(a) overlaps consumed_seed_set()")
    if min(every) <= rehearsal_used_max:
        bad.append(f"(b) seed {min(every)} is not above the largest rehearsal seed {rehearsal_used_max}")
    if set(every) & original:
        bad.append("(c) overlaps the never-run original V6 seeds")
    if set(every) & other:
        bad.append("(d) overlaps another row of the executing task")
    if bad:
        raise SystemExit("V6-F seed checks failed: " + "; ".join(bad))
    return {"a": True, "b": True, "c": True, "d": True, "e": True}


# -- producers ----------------------------------------------------------------------------------
def step_v6f_shard(a: argparse.Namespace) -> int:
    if not 0 <= a.shard < V6F_CPU_SHARDS:
        raise SystemExit(f"v6f-shard: --shard outside 0..{V6F_CPU_SHARDS - 1}")
    check_producer_seed(s6.seed_of(V6F_R_CPU, a.shard))
    return s6._produce(a, f"v6f-s{a.shard}", "v6f", "nominal", V6F_R_CPU, "warp-cpu", "float64",
                       row=V6F_ROW, step="v6f-shard")  # fmt: skip


def step_v6f_sens_shard(a: argparse.Namespace) -> int:
    if a.variant not in V6F_SENS_RUNS:
        raise SystemExit(f"v6f-sens-shard: --variant must be one of {V6F_SENS_RUNS}")
    if not 0 <= a.shard < V6F_SENS_SHARDS:
        raise SystemExit(f"v6f-sens-shard: --shard outside 0..{V6F_SENS_SHARDS - 1}")
    check_producer_seed(s6.seed_of(V6F_R_SENS, a.shard))
    return s6._produce(a, f"v6f-sens-{a.variant}-s{a.shard}", "v6f-sens", a.variant, V6F_R_SENS,
                       "warp-cpu", "float64", row="V6-sens-F", step="v6f-sens-shard")  # fmt: skip


def step_v6f_cuda_shard(a: argparse.Namespace) -> int:
    if not 0 <= a.shard < V6F_CUDA_SHARDS:
        raise SystemExit(f"v6f-cuda-shard: --shard outside 0..{V6F_CUDA_SHARDS - 1}")
    check_producer_seed(s6.seed_of(V6F_R_CUDA, a.shard))
    backend, precision = s6.cuda_target()
    return s6._produce(a, f"v6f-cuda-s{a.shard}", "v6f-cuda", "nominal", V6F_R_CUDA, backend, precision,
                       row=V6F_ROW, step="v6f-cuda-shard")  # fmt: skip


# -- timing and shard plans ---------------------------------------------------------------------
def cuda_shard_plan(rate_hist_per_s: float, total_histories: int = V6F_CUDA_SHARDS * s6.V6_SHARD_N,
                    *, planned_s: float = V6F_STEP_LIMIT_PLANNED_S,
                    hard_s: float = V6F_STEP_LIMIT_HARD_S) -> dict[str, Any]:  # fmt: skip
    """CUDA shard sizing from the measured timing rate: shards needed so that one step plans at most ``planned_s``
    (1200 s) and the per-shard time at the declared 1e6-history shards against ``hard_s`` (1500 s)."""
    if not rate_hist_per_s > 0.0:
        raise ValueError("rate must be positive")
    need = max(1, math.ceil(total_histories / rate_hist_per_s / planned_s))
    declared = V6F_CUDA_SHARDS
    per = total_histories / rate_hist_per_s / declared
    return {"rate_hist_per_s": rate_hist_per_s, "total_histories": total_histories, "planned_s_total": total_histories / rate_hist_per_s,
            "declared_shards": declared, "planned_s_per_declared_shard": per, "shards_needed_for_planned_limit": need,
            "planned_limit_s": planned_s, "hard_limit_s": hard_s, "ok": bool(declared >= need and per <= planned_s)}  # fmt: skip


def v6f_plans() -> dict[str, tuple[int, str, int | None]]:
    """``name -> (histories, rate key, declared shards)`` of the V6-F CPU families, in the format of
    ``steps_v5d.declared_plans`` for ``shard_plan_c`` (rate key ``v5``; planning 184 s x 1.672 ~ 308 s per shard).
    Not registered in ``declared_plans`` here (suite constants are not touched): the throughput step of the
    executing suite adds these entries."""
    n = s6.V6_SHARD_N
    return {
        f"v6f gating ({V6F_CPU_SHARDS} shards x 1e6)": (V6F_CPU_SHARDS * n, "v5", V6F_CPU_SHARDS),
        f"v6f-sens ({len(V6F_SENS_RUNS)} runs x {V6F_SENS_SHARDS} shards x 1e6)":
            (len(V6F_SENS_RUNS) * V6F_SENS_SHARDS * n, "v5", len(V6F_SENS_RUNS) * V6F_SENS_SHARDS),
    }  # fmt: skip


def step_v6f_timing(a: argparse.Namespace) -> int:
    """One 1e6-history CUDA float32 shard of the V6-F beam at seed 20514000; timing only, no partial for any
    consumer. Fixes the CUDA shard sizes (<= 1200 s planned, <= 1500 s hard) and the budget line."""
    backend, precision = s6.cuda_target()
    n = v5b.scaled(s6.V6_SHARD_N, a.scale, 2000, s6.V6_BATCHES)
    cfg = s6.v6_config("nominal", n, V6F_TIMING_SEED, backend, precision, a.timeout)
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    ok = bool(v5.clean(res))
    doc = {"step": "v6f-timing", "row": "none (timing only; never evidence)", "seed": V6F_TIMING_SEED,
           "backend": backend, "precision": precision, "n": n, "wall_s": wall, "hist_per_s": n / wall,
           "effective": {"backend": res.effective_config.backend, "precision": res.effective_config.precision,
                         "device": res.device},
           "plan": cuda_shard_plan(n / wall), "counters_clean": ok, "reduced": n < s6.V6_SHARD_N,
           "pass": bool(ok)}  # fmt: skip
    return v5b.finish5b(doc, s6.V6_SHARD_N, n, n < s6.V6_SHARD_N)


# -- comparator (pure) --------------------------------------------------------------------------
def aura_dose(r_cm: float, p: dict[str, float]) -> float:
    """MD-fit aura dose of the source, ``p23 E_2D(r_n, r_au)``, ``r_au = p24 (1 + p25 r_n)``, ``r_n = r / (R0/3)``,
    ``E_2D(r, r_e) = exp(-r / r_e) / (2 pi r_e^2)``; z-independent, MeV/g per primary (v6f-followup 2)."""
    rn = r_cm / (p["p1"] / 3.0)
    rau = p["p24"] * (1.0 + p["p25"] * rn)
    return p["p23"] * math.exp(-rn / rau) / (2.0 * math.pi * rau * rau)


def _binom_p(n: int, k: int) -> float:
    if n == 0:
        return 1.0
    m = min(k, n - k)
    return float(min(1.0, 2.0 * sum(math.comb(n, i) for i in range(m + 1)) / 2.0**n))


def point_sets(z: NDArray[np.float64], r: NDArray[np.float64], d_meas: NDArray[np.float64],
               parms: dict[str, float], *, z_f: float | None = None,
               aura_max: float | None = V6F_AURA_MAX) -> dict[str, Any]:  # fmt: skip
    """P_F,0 and H_F from the measured table and the frozen / aura parameters only (never model output).
    ``z_f=None``: the deepest data depth <= p1 (1 + p2). ``aura_max=None``: no aura exclusion (contrary reading a).
    Every present point gets its aura share in ``aura_share``."""
    z_lim = parms["p1"] * (1.0 + parms["p2"])
    if z_f is None:
        z_f = max(float(x) for x in z if float(x) <= z_lim + 1e-9)
    present = np.isfinite(d_meas)
    i0 = int(np.where(r == 0)[0][0])
    r0_max = float(np.nanmax(d_meas[:, i0]))
    p0, p0_report = [], []
    for j in range(len(z)):
        if present[j, i0]:
            (p0 if float(z[j]) <= z_f + 1e-9 else p0_report).append((j, i0))
    hf: list[tuple[int, int]] = []
    excluded: dict[str, list[tuple[int, int]]] = {"beyond_z_f": [], "below_floor": [], "aura": []}
    aura_share: dict[tuple[int, int], float] = {}
    for j in range(len(z)):
        for i in range(len(r)):
            if r[i] < 1 or not present[j, i]:
                continue
            aura_share[(j, i)] = aura_dose(float(r[i]), parms) / float(d_meas[j, i])
            if float(z[j]) > z_f + 1e-9:
                excluded["beyond_z_f"].append((j, i))
            elif float(d_meas[j, i]) < V6F_FLOOR * r0_max:
                excluded["below_floor"].append((j, i))
            elif aura_max is not None and aura_share[(j, i)] > aura_max:
                excluded["aura"].append((j, i))
            else:
                hf.append((j, i))
    return {"z_f": z_f, "p_f0": p0, "r0_report_only": p0_report, "h_f": hf, "excluded": excluded,
            "aura_share": aura_share, "r0_max_measured": r0_max}  # fmt: skip


def beam_screen(h_f: list[tuple[int, int]], d_meas: NDArray[np.float64], nominal: NDArray[np.float64],
                sens: dict[str, NDArray[np.float64]], *, n_sigma: float = V6F_SCREEN_NSIGMA,
                tol: float = V6F_SCREEN_TOL, model_normalised: bool = False) -> dict[str, Any]:  # fmt: skip
    """H_g: points with ``|D_var - D_nom| / D_meas + n_sigma sigma' / D_meas < tol`` for every one of the six beam
    variants; sigma' is the paired batch SE of ``D_var - D_nom`` (batches on axis 0 of ``[B, n_z, n_r]``). A
    non-finite metric excludes the point. ``model_normalised`` is contrary reading (d): ``|D_var / D_nom - 1| +
    n_sigma sigma_delta`` (paired ratio SE, ``gottschalk_v6.ratio_paired``)."""
    if set(sens) != set(V6F_BEAM_VARIANTS):
        raise ValueError(f"beam screen needs exactly the six variants {V6F_BEAM_VARIANTS}")
    nb = nominal.shape[0]
    metric = {p: 0.0 for p in h_f}
    per_variant: dict[str, dict[tuple[int, int], float]] = {}
    for name in V6F_BEAM_VARIANTS:
        var = sens[name]
        if var.shape != nominal.shape:
            raise ValueError("sensitivity batches do not pair with the nominal")
        diff = var - nominal
        se = diff.std(axis=0, ddof=1) / math.sqrt(nb)
        mv: dict[tuple[int, int], float] = {}
        if model_normalised:
            from_ratio = np.abs(var.mean(axis=0) / nominal.mean(axis=0) - 1.0)
            ym = nominal.mean(axis=0)
            rr = var.mean(axis=0) / ym
            se_r = np.sqrt(((var - rr[None] * nominal) ** 2).sum(axis=0) / (nb * (nb - 1))) / np.abs(ym)
        for (j, i) in h_f:
            if model_normalised:
                m = float(from_ratio[j, i] + n_sigma * se_r[j, i])
            else:
                m = float((abs(diff[:, j, i].mean()) + n_sigma * se[j, i]) / d_meas[j, i])
            m = m if math.isfinite(m) else math.inf
            mv[(j, i)] = m
            metric[(j, i)] = max(metric[(j, i)], m)
        per_variant[name] = mv
    h_g = [p for p in h_f if metric[p] < tol]
    return {"h_g": h_g, "metric_max": metric, "per_variant": per_variant, "tolerance": tol, "n_sigma": n_sigma,
            "n_screened_out": len(h_f) - len(h_g)}  # fmt: skip


def _residuals(d_meas: NDArray[np.float64], batches: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    mean = batches.mean(axis=0)
    se = batches.std(axis=0, ddof=1) / math.sqrt(batches.shape[0])
    return mean / d_meas - 1.0, se


def sign_test(h_g: list[tuple[int, int]], z: NDArray[np.float64], r: NDArray[np.float64],
              resid: NDArray[np.float64]) -> dict[str, Any]:  # fmt: skip
    """No r-trend: per depth with H_g points in both the 1-3 cm and the >= 4 cm group, the sign of
    ``median(>= 4 cm) - median(1-3 cm)``; ties are excluded and counted; exact two-sided binomial p."""
    n_pos = n_neg = n_tie = 0
    for j in range(len(z)):
        lo = [resid[jj, i] for (jj, i) in h_g if jj == j and 1 <= r[i] <= 3]
        hi = [resid[jj, i] for (jj, i) in h_g if jj == j and r[i] >= 4]
        if not lo or not hi:
            continue
        d = float(np.median(hi) - np.median(lo))
        if d > 0:
            n_pos += 1
        elif d < 0:
            n_neg += 1
        else:
            n_tie += 1
    used = n_pos + n_neg
    p = _binom_p(used, n_pos)
    return {"n_depths_after_ties": used, "n_depths_with_ties": used + n_tie, "n_pos": n_pos, "n_neg": n_neg,
            "n_ties_excluded": n_tie, "p": p, "pass": bool(p > 0.01)}  # fmt: skip


def coverage_minimum(h_g: list[tuple[int, int]], r: NDArray[np.float64], sign: dict[str, Any]) -> dict[str, Any]:
    radii = sorted({float(r[i]) for (_, i) in h_g})
    ok = bool(len(h_g) >= V6F_MIN_POINTS and len(radii) >= V6F_MIN_RADII
              and sign["n_depths_after_ties"] >= V6F_MIN_DEPTHS)  # fmt: skip
    return {"n_points": len(h_g), "radii_cm": radii, "n_radii": len(radii),
            "sign_test_depths_after_ties": sign["n_depths_after_ties"],
            "minimum": {"points": V6F_MIN_POINTS, "radii": V6F_MIN_RADII, "depths": V6F_MIN_DEPTHS}, "ok": ok}  # fmt: skip


def evaluate_backend(z: NDArray[np.float64], r: NDArray[np.float64], d_meas: NDArray[np.float64],
                     sets: dict[str, Any], h_g: list[tuple[int, int]],
                     batches: NDArray[np.float64]) -> dict[str, Any]:  # fmt: skip
    """Verdict components of one backend's sample on the fixed point sets (the same ``h_g`` for every backend)."""
    resid, se = _residuals(d_meas, batches)
    mean = batches.mean(axis=0)

    def pt(j: int, i: int) -> dict[str, Any]:
        return {"z_cm": float(z[j]), "r_cm": float(r[i]), "D_meas": float(d_meas[j, i]), "D_ionmc": float(mean[j, i]),
                "residual": float(resid[j, i]), "se": float(se[j, i]),
                "rel_se": float(se[j, i] / mean[j, i]) if mean[j, i] > 0 else None,
                "aura_share": sets["aura_share"].get((j, i))}  # fmt: skip

    r0_pts = [{**pt(j, i), "pass": bool(abs(resid[j, i]) <= 0.02),
               "se_flag": bool(mean[j, i] > 0 and se[j, i] / mean[j, i] > 2 * V6F_R0_TARGET_REL_SE)}
              for (j, i) in sets["p_f0"]]  # fmt: skip
    r0_report = [pt(j, i) for (j, i) in sets["r0_report_only"]]
    r0_pass = bool(r0_pts) and all(p["pass"] for p in r0_pts)
    halo_pts = [{**pt(j, i), "pass": bool(abs(resid[j, i]) <= 0.15),
                 "se_flag": bool(mean[j, i] > 0 and se[j, i] / mean[j, i] > 2 * V6F_HALO_TARGET_REL_SE)}
                for (j, i) in h_g]  # fmt: skip
    absres = np.array([abs(p["residual"]) for p in halo_pts])
    frac = float(np.mean(absres <= 0.15)) if absres.size else 0.0
    med = float(np.median(absres)) if absres.size else math.inf
    sign = sign_test(h_g, z, r, resid)
    cov = coverage_minimum(h_g, r, sign)
    halo_clauses = bool(halo_pts and frac >= 0.90 and med <= 0.10 and sign["pass"])
    halo_status = ("pass" if halo_clauses else "fail") if cov["ok"] else NOT_ESTABLISHED
    return {"r0": {"points": r0_pts, "pass": r0_pass, "n_points": len(r0_pts), "report_only": r0_report},
            "halo": {"points": halo_pts, "frac_pass": frac, "median_abs": med, "sign_test": sign,
                     "clauses_pass": halo_clauses, "status": halo_status},
            "coverage": cov, "overall_pass": bool(r0_pass and halo_status == "pass")}  # fmt: skip


def backend_agreement(cpu: dict[str, Any], cuda: dict[str, Any]) -> dict[str, Any]:
    """Amendment 17 (c)8: V6-F fails if the two backends differ on r = 0, halo, overall or coverage. The
    per-point paired z (CPU vs CUDA, from the batch SEs) is listed for every point; |z| > 4 is a finding."""
    comps = {"r0": (cpu["r0"]["pass"], cuda["r0"]["pass"]),
             "halo": (cpu["halo"]["status"], cuda["halo"]["status"]),
             "overall": (cpu["overall_pass"], cuda["overall_pass"]),
             "coverage": (cpu["coverage"]["ok"], cuda["coverage"]["ok"])}  # fmt: skip
    disagree = sorted(k for k, (x, y) in comps.items() if x != y)
    zs: list[dict[str, Any]] = []
    for grp in ("r0", "halo"):
        for pc, pg in zip(cpu[grp]["points"], cuda[grp]["points"], strict=True):
            den = math.hypot(pc["se"], pg["se"])
            zv = (pc["D_ionmc"] - pg["D_ionmc"]) / den if den > 0 else math.nan
            zs.append({"group": grp, "z_cm": pc["z_cm"], "r_cm": pc["r_cm"], "z": zv})
    findings = [q for q in zs if math.isfinite(q["z"]) and abs(q["z"]) > 4.0]
    return {"components": {k: list(v) for k, v in comps.items()}, "disagreement": disagree,
            "agree": not disagree, "paired_z": zs, "findings_abs_z_gt_4": findings}  # fmt: skip


def evaluate_v6f(z: NDArray[np.float64], r: NDArray[np.float64], d_meas: NDArray[np.float64],
                 parms: dict[str, float], cpu_batches: NDArray[np.float64], cuda_batches: NDArray[np.float64],
                 sens_nominal: NDArray[np.float64], sens: dict[str, NDArray[np.float64]]) -> dict[str, Any]:  # fmt: skip
    """The V6-F verdict. ``cpu_batches`` / ``cuda_batches`` ``[B, n_z, n_r]`` are the gating samples,
    ``sens_nominal`` and ``sens`` (six beam variants) the common-random-number sensitivity runs. Passes iff the CPU
    r = 0 and halo components pass (halo only if the coverage minimum holds), the CUDA sample passes on the same
    sets, and the backends agree. Contrary readings (a)-(e) of v6f-followup 9 are reported, never gating."""
    sets = point_sets(z, r, d_meas, parms)
    scr = beam_screen(sets["h_f"], d_meas, sens_nominal, sens)
    cpu = evaluate_backend(z, r, d_meas, sets, scr["h_g"], cpu_batches)
    cuda = evaluate_backend(z, r, d_meas, sets, scr["h_g"], cuda_batches)
    agree = backend_agreement(cpu, cuda)
    passed = bool(cpu["overall_pass"] and cuda["overall_pass"] and agree["agree"])
    if cpu["halo"]["status"] == NOT_ESTABLISHED or cuda["halo"]["status"] == NOT_ESTABLISHED:
        verdict = NOT_ESTABLISHED if agree["agree"] and cpu["r0"]["pass"] else "fail"
    else:
        verdict = "pass" if passed else "fail"
    if verdict != "pass":
        passed = False
    contrary = _contrary(z, r, d_meas, parms, cpu_batches, sens_nominal, sens, sets, scr)
    return {"sets": {"z_f": sets["z_f"], "n_p_f0": len(sets["p_f0"]), "n_r0_report_only": len(sets["r0_report_only"]),
                     "n_h_f": len(sets["h_f"]), "n_h_g": len(scr["h_g"]),
                     "excluded_counts": {k: len(v) for k, v in sets["excluded"].items()},
                     "aura_share_points": [{"z_cm": float(z[j]), "r_cm": float(r[i]), "aura_share": s}
                                           for (j, i), s in sorted(sets["aura_share"].items())]},
            "beam_screen": {"n_screened_out": scr["n_screened_out"], "tolerance": scr["tolerance"],
                            "n_sigma": scr["n_sigma"],
                            "per_point_max_metric": [{"z_cm": float(z[j]), "r_cm": float(r[i]), "metric": m}
                                                     for (j, i), m in sorted(scr["metric_max"].items())]},
            "cpu": cpu, "cuda": cuda, "backend_agreement": agree, "verdict": verdict, "pass": passed,
            "contrary_readings": contrary, "monitor_note": s6.MONITOR_NOTE}  # fmt: skip


def _contrary(z: NDArray[np.float64], r: NDArray[np.float64], d_meas: NDArray[np.float64], parms: dict[str, float],
              cpu_batches: NDArray[np.float64], nominal: NDArray[np.float64],
              sens: dict[str, NDArray[np.float64]], sets: dict[str, Any], scr: dict[str, Any]) -> dict[str, Any]:  # fmt: skip
    """Readings (a)-(e) of v6f-followup 9 on the CPU sample (report only)."""

    def halo_of(h: list[tuple[int, int]]) -> dict[str, Any]:
        res = evaluate_backend(z, r, d_meas, sets, h, cpu_batches)
        return {"n_points": len(h), "status": res["halo"]["status"], "frac_pass": res["halo"]["frac_pass"],
                "median_abs": res["halo"]["median_abs"], "sign_p": res["halo"]["sign_test"]["p"]}

    out: dict[str, Any] = {}
    s_a = point_sets(z, r, d_meas, parms, aura_max=None)
    out["a_no_aura_exclusion"] = halo_of(beam_screen(s_a["h_f"], d_meas, nominal, sens)["h_g"])
    zs = sorted(float(x) for x in z)
    deeper = [x for x in zs if x > sets["z_f"]]
    if deeper:
        s_b = point_sets(z, r, d_meas, parms, z_f=deeper[0])
        r0b = evaluate_backend(z, r, d_meas, s_b, [], cpu_batches)["r0"]
        out["b_z_f_next_depth"] = {"z_f": deeper[0], "r0_pass": r0b["pass"],
                                   **halo_of(beam_screen(s_b["h_f"], d_meas, nominal, sens)["h_g"])}  # fmt: skip
    out["c_one_sigma_screen"] = halo_of(beam_screen(sets["h_f"], d_meas, nominal, sens, n_sigma=1.0)["h_g"])
    out["d_model_normalised_screen"] = halo_of(
        beam_screen(sets["h_f"], d_meas, nominal, sens, model_normalised=True)["h_g"])
    screened = [p for p in sets["h_f"] if p not in set(scr["h_g"])]
    out["e_screened_out_points"] = halo_of(screened) if screened else {"n_points": 0}
    return out


# -- compare step -------------------------------------------------------------------------------
def step_v6f_compare(a: argparse.Namespace) -> int:
    g = s6.gv6()
    check_v6f_seeds_base = a.seed
    s6.check_base(check_v6f_seeds_base)
    doc: dict[str, Any] = {"step": "v6f-compare", "row": V6F_ROW, "lineage": g.LINEAGE}
    reasons: list[str] = []
    try:
        cpu_docs = s6._load_group(a, cpu_names(), V6F_R_CPU, "v6f", "nominal", "warp-cpu", "float64", V6F_CPU_SHARDS)
        cuda_docs = s6._load_group(a, cuda_names(), V6F_R_CUDA, "v6f-cuda", "nominal", *s6.cuda_target(),
                                   V6F_CUDA_SHARDS)  # fmt: skip
        sens_docs = {v: s6._load_group(a, tuple(f"v6f-sens-{v}-s{k}.json" for k in range(V6F_SENS_SHARDS)),
                                       V6F_R_SENS, "v6f-sens", v, "warp-cpu", "float64", V6F_SENS_SHARDS)
                     for v in V6F_SENS_RUNS}  # fmt: skip
        common = s6._consistent(cpu_docs)
        _, wdigest = s6.weights()
        if common["weights_digest"] != wdigest or common["beam"] != s6.beam_table()["nominal"]:
            raise SystemExit("shard weights digest or beam differ from the current code")
        for v, ds in sens_docs.items():
            if s6._consistent(ds)["beam"] != s6.beam_table()[v]:
                raise SystemExit(f"v6f-sens {v}: beam differs from the current code")
        if s6._consistent(cuda_docs)["beam"] != common["beam"]:
            raise SystemExit("CUDA shards use another beam than the CPU shards")
        sens_seed = {v: [d["seed"] for d in ds] for v, ds in sens_docs.items()}
        if len({tuple(s) for s in sens_seed.values()}) != 1:
            raise SystemExit("v6f-sens runs do not share seeds (common random numbers)")
        checks = check_v6f_seeds([d["seed"] for d in cpu_docs], next(iter(sens_seed.values())),
                                 [d["seed"] for d in cuda_docs])  # fmt: skip
        z, r, meas, data = s6._data(g)
        tex, _ = g.load_halo_tex()
        parms = g.parse_parms(tex)
        stack = s6._stack
        ver = evaluate_v6f(z, r, meas, parms, stack(cpu_docs), stack(cuda_docs), stack(sens_docs["nominal"]),
                           {v: stack(sens_docs[v]) for v in V6F_BEAM_VARIANTS})  # fmt: skip
        doc.update(data=data, beam=s6._beam_doc(g, common), n_histories_cpu=sum(d["n"] for d in cpu_docs),
                   n_histories_cuda=sum(d["n"] for d in cuda_docs), seed_checks=checks,
                   seeds={"cpu": [d["seed"] for d in cpu_docs], "cuda": [d["seed"] for d in cuda_docs],
                          "sens": next(iter(sens_seed.values()))}, **ver)  # fmt: skip
        if ver["verdict"] != "pass":
            reasons.append(f"verdict: {ver['verdict']}")
        if ver["backend_agreement"]["disagreement"]:
            reasons.append("backend disagreement on " + ", ".join(ver["backend_agreement"]["disagreement"]))
    except (ValueError, KeyError) as exc:
        doc.update(error=f"{type(exc).__name__}: {exc}", verdict="fail")
        doc["pass"] = False
        reasons.append(str(exc))
    doc["pass"] = bool(doc.get("pass", False))
    doc["reasons"] = reasons
    doc["attestation"] = v5.attestation_block(a)
    n_cpu = V6F_CPU_SHARDS * s6.V6_SHARD_N
    return v5b.finish5b(doc, n_cpu, n_cpu, False)


STEPS = {
    "v6f-shard": step_v6f_shard, "v6f-sens-shard": step_v6f_sens_shard, "v6f-cuda-shard": step_v6f_cuda_shard,
    "v6f-compare": step_v6f_compare, "v6f-timing": step_v6f_timing,
}  # fmt: skip


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--workers", default="1", help="accepted for uniformity; one process")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--variant", default="nominal", help="v6f-sens-shard: nominal, sx-lo, ... rho-hi")
    ap.add_argument("--out-dir", default="samples")
    ap.add_argument("--dirs", nargs="+", default=["."])
    ap.add_argument("--partials-manifest", default=None)
    ap.add_argument("--seed-base", "--seed", dest="seed", type=int, default=V6F_BASE)
    ap.add_argument("--timeout", type=float, default=None)
    args = ap.parse_args(argv)
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    s6.check_base(args.seed)
    base.SEED_BASE = args.seed
    v5c.v4.base.SEED_BASE = args.seed
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
