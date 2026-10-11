"""Validation steps of the V3-005 acceptance plan, slice C, part 2 (suites ``lv5c`` and ``hr5c``).

Usage::

    python validation/scripts/transport/steps_v5d.py <step> [--scale F] [--seed-base N] ...

``steps_v5c.py`` holds the configuration builders, V7-R and V11; this module holds the remaining lv5c
producers and consumers (plan Amendments 14 and 17): ``lv5c-throughput`` (r1), ``r1-d2-regression``
(r3), ``v8-lv`` (r2), ``p7-sampler`` (r15), the V5 re-run ``v5-ionmc`` / ``v5-compare`` /
``v5-attribution`` (r4-7), V2b-el, V7-scan-el, V7-shift-el, E1-B, E-shape, the deterministic rows
(``nuclear_v10.py``) and the hr5c steps (V8-stat with V8-nuc, V7 f32/f64 on CUDA). The frozen steps of
``steps_v5.py`` / ``steps_v5b.py`` are imported, never edited; every elastic-on configuration is built by
``steps_v5c.nuc_config_c`` / ``wcfg_c``.

Seeds: ``seed = base + 1000 * r_index + k`` (Amendment 14 (j)); ``R_INDEX_C`` / ``R_INDEX_HRC`` below.
Every step that combines samples checks that their seed sets are pairwise disjoint and disjoint from the
consumed evidence seed sets (``steps_v5c.consumed_seed_set``); E-shape uses the absolute rehearsal-family
seeds 20509500 + k (Amendment 17 (b)3).
"""

from __future__ import annotations

import argparse
import contextlib
import json
import math
import sys
import time
from collections.abc import Iterator
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

import run_suite
import steps as base
import steps_v4 as v4
import steps_v5 as v5
import steps_v5b as v5b
import steps_v5c as c
from ionmc import materials as M
from ionmc.config import DiagnosticsOptions, SimulationConfig
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Simulation

HERE = Path(__file__).resolve().parent
LV5C_SEED_BASE = c.LV5C_SEED_BASE
HR5C_SEED_BASE = 20491004
R_INDEX_C = {"throughput": 1, "v8-lv": 2, "r1-d2": 3, "v5-150-on": 4, "v5-150-off": 5, "v5-200-on": 6,
             "v5-200-off": 7, "v2b-el": 8, "v7-scan-el": 9, "v7-shift-el": 10, "v6": 12, "e1-b": 13,
             "v11": 14, "p7": 15, "v7-r": 16, "v7-r-diag": 17, "v10-a-xcheck": 18, "x-endf": 19,
             "v6-sens": 20}  # fmt: skip
R_INDEX_HRC = {"v8-stat": 1, "v7-f32": 2, "v6-cuda": 3}
E_SHAPE_BASE = 20509500  # absolute seeds of the rehearsal family, Amendment 15 (b)6 / 17 (b)3
REHEARSAL_RANGE = (20505000, 20509999)
E1B_SHARDS: int | None = 3
"""Number of ``e1-b-s{k}`` shards (>= 2, total >= 1e5; Amendment 14a: 3); ``None`` makes the step ``e1-b-unsized``
fail closed (Amendment 17 (b)6). ``run_suite.E1B_SHARDS`` must equal it (a test)."""
E1B_N = 100_000
E1B_BATCHES = 20


def e1b_shard_n(shards: int) -> int:
    """Histories per ``e1-b`` shard: ``ceil(1e5 / shards)`` rounded up to a multiple of the 20 batches (the
    configuration requires ``n_histories % n_batches == 0``): 33340 for 3 shards (total 100020 >= 1e5; the
    split 33334/33333/33333 is not a multiple of 20)."""
    return -(-(-(-E1B_N // shards)) // E1B_BATCHES) * E1B_BATCHES


V8NUC_PYTHON_SHARDS, V8NUC_PYTHON_N, V8NUC_BATCHES = 4, 30_000, 20  # 4 x 20 batches x 1500 (Amendment 14 (i))
V8NUC_WARP_N, V8NUC_WARP_BATCHES = 1_000_000, 100
V8_SAMPLES_C = {  # name -> (backend, precision, shard index k, histories, batches), Amendment 17 (b)2
    **{f"python-s{k}": ("python", "float64", k, V8NUC_PYTHON_N, V8NUC_BATCHES) for k in range(4)},
    "cpu64": ("warp-cpu", "float64", 4, V8NUC_WARP_N, V8NUC_WARP_BATCHES),
    "cuda32": ("warp-cuda", "float32", 5, V8NUC_WARP_N, V8NUC_WARP_BATCHES),
    "cuda64": ("warp-cuda", "float64", 6, V8NUC_WARP_N, V8NUC_WARP_BATCHES),
}  # fmt: skip
V8_PAIRS_C = (("python", "cpu64"), ("cpu64", "cuda32"), ("cpu64", "cuda64"), ("cuda32", "cuda64"))

NUC_LOCAL = v5b.NUC_LOCAL
REC = TallyRequest("recoil", "dose", "edep", species=("elastic_recoil",))
"""The recoil depth channel of the V5 "on" producers (Amendment 17 (b)5): the energy of the elastic p + A
recoil (pseudo-species ``elastic_recoil``) per depth bin; tally-neutral for the grid energy."""

# R1-D2 (C) pre-registered rule (Amendment 17 (b)1)
R1D2_Z_MIN = 3.0
R1D2_DECLARED = {"peak_over_plateau": -1, "plateau": +1, "lateral_rms": +1}  # sign of (elastic - no elastic)
R1D2_TOTAL_TOL = 1.0e-3  # |delta total deposit| / E0
R1D2_LATERAL_HALF_MM, R1D2_LATERAL_BIN_MM, R1D2_SLAB = 50.0, 1.0, (0.45, 0.55)

# attribution ladder (Amendment 14 (f), Amendment 17 (b)4)
STOP_LOW, STOP_ABS, RUNG2_LOW = -0.02, 0.05, 0.02
EXPECTATION = {  # plan (f), planning estimates, not criteria
    "150": {"peak_over_plateau_ratio": 0.985, "range": [0.95, 1.03], "delta_idd_percent": [-2.5, [-8, 3]]},
    "200": {"peak_over_plateau_ratio": 0.986, "range": [0.90, 0.99], "note": "at delta = 5 mm"},
    "basis": "straight-line CSDA estimate with Barashenkov sigma_el and P(p-O) 0.21 / 0.24, +-4 % on the "
    "factor; the TOPAS X-case route gives 0.998 / 0.996; plateau, R80 and total deposit stay passing",
}  # fmt: skip


# -- seeds ---------------------------------------------------------------------------------------
def seed_c(row: str, k: int = 0) -> int:
    return base.SEED_BASE + 1000 * R_INDEX_C[row] + k


def seed_hrc(row: str, k: int = 0) -> int:
    return base.SEED_BASE + 1000 * R_INDEX_HRC[row] + k


def seed_guard(named: dict[str, list[int]], *, allow_rehearsal: bool = False) -> dict[str, Any]:
    """Fail closed (``SystemExit``) unless the seed lists of ``named`` are pairwise disjoint and disjoint
    from every consumed evidence seed (plan l.380). Seeds of the rehearsal family 20505000-20509999 are
    refused unless ``allow_rehearsal`` (E-shape). Returns the record stored in the step document."""
    names = list(named)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            if set(named[a]) & set(named[b]):
                raise SystemExit(f"seed sets {a!r} and {b!r} intersect (fail closed)")
        if len(set(named[a])) != len(named[a]):
            raise SystemExit(f"seed list {a!r} repeats a seed")
    allseeds = {s for v in named.values() for s in v}
    bad = sorted(allseeds & c.consumed_seed_set())
    if bad:
        raise SystemExit(f"seeds {bad} belong to a consumed evidence seed set (fail closed)")
    reh = sorted(s for s in allseeds if REHEARSAL_RANGE[0] <= s <= REHEARSAL_RANGE[1])
    if reh and not allow_rehearsal:
        raise SystemExit(f"seeds {reh} are rehearsal-family seeds, never evidence (fail closed)")
    return {"seed_lists": {k: list(v) for k, v in named.items()}, "disjoint": True,
            "disjoint_from_consumed": True, "consumed_bases": list(c.CONSUMED_SEED_BASES_C)}  # fmt: skip


def finish(doc: dict[str, Any], frozen: Any, used: Any, reduced: bool) -> int:
    return v5b.finish5b(doc, frozen, used, reduced)


def _cmp() -> Any:
    import importlib.util

    path = v5b.REPO / "validation" / "scripts" / "reference" / "compare_idd_v5.py"
    spec = importlib.util.spec_from_file_location("compare_idd_v5", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def status_partial(a: argparse.Namespace, name: str, doc: dict[str, Any]) -> str:
    """Sealed status partial ``<name>.json`` (a step cannot read another step's status; the consumers
    ``v5-attribution`` load these through the protected-manifest path)."""
    return c.write_partial_c(a, name, {"row": name, "status_of": doc["step"], **doc})


@contextlib.contextmanager
def scaled_elastic_rows(factor: float) -> Iterator[None]:
    """Multiply every Sigma_el (and its majorants and cumulative rows) of ``ElasticTable.material_rows``."""
    from ionmc.nuclear.elastic_tables import ElasticTable

    orig = ElasticTable.material_rows

    def scaled_fn(self: Any, material: Any, f_e: float = 0.02) -> Any:
        m = orig(self, material, f_e)
        return replace(m, sigma_mass_cm2_g=m.sigma_mass_cm2_g * factor,
                       sigma_hat_window=m.sigma_hat_window * factor,
                       sigma_hat_end=m.sigma_hat_end * factor,
                       cum_sigma_mass_cm2_g=m.cum_sigma_mass_cm2_g * factor)  # fmt: skip

    ElasticTable.material_rows = scaled_fn  # type: ignore[method-assign,assignment]
    try:
        yield
    finally:
        ElasticTable.material_rows = orig  # type: ignore[method-assign]


def elastic_off(cfg: SimulationConfig) -> SimulationConfig:
    """``cfg`` with the elastic channel off (the lv5b physics, ``elastic=False``, no table id)."""
    return replace(cfg, physics=replace(cfg.physics, elastic=False, elastic_only=False,
                                        elastic_table_id=None))


# -- lv5c-throughput (r1) --------------------------------------------------------------------------
def measure_c(backend: str, precision: str, e: float, n: int, **kw: Any) -> dict[str, Any]:
    geo, grid, _, _ = v5.depth_box(M.WATER, e)
    if backend != "python":
        Simulation(c.wcfg_c(backend, precision, energy=e, n=200, seed=1, geometry=geo, grid=grid,
                            n_batches=2)).run()  # fmt: skip
    cfg = c.wcfg_c(backend, precision, energy=e, n=n, seed=1, geometry=geo, grid=grid, n_batches=2, **kw)
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    return v5b.rate_of(res, n, time.perf_counter() - t0)


def block_rate(n: int) -> float:
    """hist/s of ``batch_estimates`` (500-history blocks, the V7 / V8 / E1-style producers)."""
    geo, grid = v5b.coarse_depth(v5b.V7_BINS)
    nb = max(2, n // v5b.V7_BATCH_N)
    t0 = time.perf_counter()
    c_est = v5b.batch_estimates(c.wcfg_c("warp-cpu", "float64", energy=150.0, n=nb * v5b.V7_BATCH_N,
                                         seed=1, geometry=geo, grid=grid, n_batches=nb,
                                         tallies=v5b.V7_TALLIES))  # fmt: skip
    return float(nb * v5b.V7_BATCH_N / (time.perf_counter() - t0)) if c_est else 0.0


def declared_plans() -> dict[str, tuple[int, str, int | None]]:
    """``name -> (histories, rate key, declared shards)`` of every lv5c / hr5c step that must fit
    ``STEP_LIMIT_S / MARGIN`` (``None`` shards: not yet declared, reported as needed)."""
    v6 = run_suite.V6_SHARDS
    out: dict[str, tuple[int, str, int | None]] = {
        "v7r (16 shards x 4.5e6, block rate)": (16 * 4_500_000, "blocks", 16),
        "v2b-el (both variants, 1e6 each)": (v5b.V2B_N, "v2b", v5b.V2B_SHARDS),
        "v7-scan-el (1e4 + 1e5 + 1e6, plus f32 1e5)": (sum(v5b.V7_SCAN_N) + v5b.V7_F32_N, "blocks", 1),
        "v7-shift-el (3 layouts)": (3 * v5b.V7_SHIFT_N, "blocks", 1),
        "v5-ionmc (per step)": (v5b.V5_N, "v5", 1),
        "v11-ionmc (per step)": (v5b.V5_N, "v5", 1),
        "v8-nuc python sample (4 x 3e4, 150 MeV)": (V8NUC_PYTHON_SHARDS * V8NUC_PYTHON_N, "python", V8NUC_PYTHON_SHARDS),
        "e1-b (python 150 MeV, 1e5)": (E1B_N, "python", E1B_SHARDS),
        "e-shape (python 150 MeV, 2 x 2e4)": (2 * E_SHAPE_N, "python", 1),
    }  # fmt: skip
    if v6 is not None:
        out[f"v6 ({v6} shards x 1e6)"] = (v6 * 1_000_000, "v5", v6)
    return out


def shard_plan_c(items: dict[str, tuple[int, str, int | None]], rates: dict[str, float],
                 limit_s: float, margin: float) -> dict[str, Any]:
    """``name -> (histories, rate key, declared shards)``: planned time per shard at the measured rate, the shards
    needed so that one step takes at most ``limit_s / margin`` and whether the declared count suffices. An
    undeclared count (``None``) is reported with the shards needed and ``ok = None`` (fixed in Amendment 14a)."""
    out = {}
    for name, (n, key, declared) in items.items():
        r = rates[key]
        need = max(1, math.ceil(n / r / (limit_s / margin)))
        d = declared
        out[name] = {"histories": n, "rate_key": key, "rate_hist_per_s": r, "planned_s_total": n / r,
                     "planned_s_per_declared_shard": None if d is None else n / r / d, "declared_shards": d,
                     "shards_needed_with_margin": need, "histories_per_step_with_margin": int(limit_s / margin * r),
                     "limit_s": limit_s, "ok": None if d is None else bool(d >= need)}  # fmt: skip
    return out


def plan_table(rates: dict[str, float], plans: dict[str, tuple[int, str, int | None]]) -> dict[str, Any]:
    """The shard table against ``run_suite.HOST_STEP_WALL_LIMIT_S`` (1500 s; Amendment 14 (k)) with the margin of
    ``steps_v5.MARGIN``; the timeout floors (3300 / 3600) are kill guards and play no role here."""
    return shard_plan_c(plans, rates, float(run_suite.HOST_STEP_WALL_LIMIT_S), v5.MARGIN)


E_SHAPE_N = 20_000


def step_throughput(a: argparse.Namespace) -> int:
    npy = v5b.scaled(v5b.THROUGHPUT_PY_N, a.scale, 20)
    nw = v5b.scaled(v5b.THROUGHPUT_WARP_N, a.scale, 500)
    out: dict[str, Any] = {}
    for e in v5.THROUGHPUT_ENERGIES:
        out[f"python-f64@{e:g}"] = measure_c("python", "float64", e, npy, timeout=a.timeout)
    for prec in ("float64", "float32"):
        out[f"warp-cpu-{prec}@150"] = measure_c("warp-cpu", prec, 150.0, nw)
    geo, grid, _, _ = v5.depth_box(M.WATER, 150.0)
    for tag, (smax, frac) in zip("ab", v5b.V2B_VARIANTS, strict=True):
        cfg = c.wcfg_c("warp-cpu", "float64", energy=150.0, n=nw, seed=1, geometry=geo, grid=grid,
                       mcs=False, straggling=False, max_step=smax, frac=frac, tallies=v5.PRIMARY)  # fmt: skip
        t0 = time.perf_counter()
        res = Simulation(cfg).run()
        out[f"v2b-{tag}"] = v5b.rate_of(res, nw, time.perf_counter() - t0)
    for e in v5b.V5_ENERGIES:
        geo5, grid5, _, _ = v5b.v5_geometry(e)
        cfg = c.wcfg_c("warp-cpu", "float64", energy=e, n=nw, seed=1, geometry=geo5, grid=grid5)
        t0 = time.perf_counter()
        res = Simulation(cfg).run()
        out[f"v5-geometry@{e:g}"] = v5b.rate_of(res, nw, time.perf_counter() - t0)
    rates = {"python": out["python-f64@150"]["hist_per_s"],  # E1-B, E-shape and V8-nuc are 150 MeV python runs
             "python@250": out["python-f64@250"]["hist_per_s"],  # reported alongside, used by no plan
             "v5": min(out[f"v5-geometry@{e:g}"]["hist_per_s"] for e in v5b.V5_ENERGIES),
             "blocks": block_rate(v5b.scaled(20 * v5b.V7_BATCH_N, a.scale, 2 * v5b.V7_BATCH_N,
                                             v5b.V7_BATCH_N))}  # fmt: skip
    rates["v2b"] = 1.0 / (1.0 / out["v2b-a"]["hist_per_s"] + 1.0 / out["v2b-b"]["hist_per_s"])
    table = plan_table(rates, declared_plans())
    ok = all(v["counters_clean"] for v in out.values())
    fits = all(v["ok"] for v in table.values() if v["ok"] is not None)
    doc = {"step": "lv5c-throughput", "table": v5.table_record(), **c.elastic_table_record(),
           "seed": seed_c("throughput"), "throughput": out, "rates_hist_per_s": rates,
           "step_limit_s": run_suite.HOST_STEP_WALL_LIMIT_S, "margin": v5.MARGIN,
           "effective_s_per_step": run_suite.HOST_STEP_WALL_LIMIT_S / v5.MARGIN,
           "timeout_floors_are_kill_guards": {"lv5c": 3300, "hr5c": 3600}, "shard_table": table,
           "shard_table_ok": bool(fits),
           "undeclared": sorted(k for k, v in table.items() if v["ok"] is None),
           "row": "none (throughput yields no verdict; lv5c r_index 1 seeds equal the consumed lv5b "
           "V7-replicate seeds, plan Amendment 14 (j))",
           "pass": bool(ok and (a.scale < 1.0 or fits))}  # fmt: skip
    return finish(doc, v5b.THROUGHPUT_WARP_N, nw, nw < v5b.THROUGHPUT_WARP_N)


# -- r1-d2-regression (r3) -------------------------------------------------------------------------
R1D2_PLAN_ROW = (
    "R1-D2 (plan Amendment 14 (d), l.206): (A) `nuclear=False`: the A16 digest is bit-identical to f3a1dd62 "
    "and the intended-change set is empty ... (B) `nuclear=True` with Sigma_el = 0 (internal switch): "
    "bit-identical to the V3-005B nuclear-on outputs (python primary-only digest `3edfce9c...` = b84fdf38, "
    "and the C13 secondary-proton record) ... (C) The record for nuclear-on outputs, with directions "
    "declared now: peak/plateau down, plateau up, lateral spread up, |delta total deposit| <= 0.1 % of E0 "
    "| (A, B) bit-identical; (C) record present, reviewed and matching every declared direction"
)


def paired_z(d: NDArray[np.float64]) -> dict[str, Any]:
    """Paired batch-mean difference, its standard error (ddof = 1), z and the 90 % interval."""
    d = np.asarray(d, dtype=float)
    n = d.size
    mean = float(d.mean())
    se = float(d.std(ddof=1) / math.sqrt(n))
    q = v5b.student_t_quantile(0.95, n - 1)
    return {"mean": mean, "se": se, "z": mean / se if se > 0 else math.inf if mean else 0.0,
            "n_batches": n, "t95": q, "ci90": [mean - q * se, mean + q * se]}  # fmt: skip


def lateral_rms_batches(batch_energy: NDArray[np.float64], half_mm: float, bin_mm: float) -> NDArray[np.float64]:
    """Dose-weighted rms lateral radius [mm] per batch of a ``[B, nx * ny]`` lateral energy map of one depth
    slab (equal voxel masses, so the energy weights equal the dose weights)."""
    n = int(round(2 * half_mm / bin_mm))
    w = np.asarray(batch_energy, dtype=float).reshape(-1, n, n)
    x = -half_mm + (np.arange(n) + 0.5) * bin_mm
    r2 = x[:, None] ** 2 + x[None, :] ** 2
    return np.sqrt((w * r2).sum(axis=(1, 2)) / w.sum(axis=(1, 2)))


def r1d2_judgement(paired: dict[str, dict[str, Any]], total: dict[str, Any], e0: float) -> dict[str, Any]:
    """Direction match per declared output: the paired mean difference has the declared sign and |z| >= 3;
    the total-deposit point criterion |delta| / E0 <= 1e-3 with its 90 % interval (reported)."""
    out: dict[str, Any] = {}
    for name, sign in R1D2_DECLARED.items():
        p = paired[name]
        out[name] = {**p, "declared_sign": sign,
                     "matches": bool(p["mean"] * sign > 0 and abs(p["z"]) >= R1D2_Z_MIN)}  # fmt: skip
    frac = {"mean": total["mean"] / e0, "ci90": [x / e0 for x in total["ci90"]], "se": total["se"] / e0}
    out["total_deposit"] = {**frac, "limit": R1D2_TOTAL_TOL,
                            "matches": bool(abs(frac["mean"]) <= R1D2_TOTAL_TOL)}  # fmt: skip
    out["all_match"] = bool(all(v["matches"] for v in out.values()))
    return out


def r1d2_measure(e: float, k: int, n: int, timeout: float | None) -> dict[str, Any]:
    """Paired CRN runs (warp-cpu float64, seed ``r3 + k`` for both variants) at ``e`` MeV."""
    from ionmc.reference import metrics as rm

    geo, grid, nz, _ = v5b.v5_geometry(e)
    seed = seed_c("r1-d2", k)
    idd: dict[bool, NDArray[np.float64]] = {}
    clean = True
    for el in (False, True):
        cfg = c.wcfg_c("warp-cpu", "float64", energy=e, n=n, seed=seed, geometry=geo, grid=grid,
                       n_batches=v5b.V5_BATCHES, timeout=timeout)  # fmt: skip
        res = Simulation(cfg if el else elastic_off(cfg)).run()
        clean &= bool(v5.clean(res) and res.energy_balance.relative_residual <= 1e-12)
        bin_g = M.WATER.density_g_cm3 * v5b.V5_DZ_MM / 10.0
        idd[el] = res.grid("dose").batch_energy_mev.reshape(v5b.V5_BATCHES, nz) / bin_g
    depth = (np.arange(nz) + 0.5) * v5b.V5_DZ_MM
    plateau = (depth >= 20.0) & (depth <= 60.0)
    bin_g = M.WATER.density_g_cm3 * v5b.V5_DZ_MM / 10.0
    m = {el: {"peak_over_plateau": idd[el].max(axis=1) / idd[el][:, plateau].mean(axis=1),
              "plateau": idd[el][:, plateau].mean(axis=1), "total": idd[el].sum(axis=1) * bin_g}
         for el in idd}  # fmt: skip
    r80 = float(rm.r80(depth, rm.normalize_to_peak(idd[False].mean(axis=0))))
    # lateral run: dose-weighted rms radius in [0.45, 0.55] R80 of the elastic-off run (1 mm grid, 50 mm)
    z0, dz = R1D2_SLAB[0] * r80, (R1D2_SLAB[1] - R1D2_SLAB[0]) * r80
    nl = int(round(2 * R1D2_LATERAL_HALF_MM / R1D2_LATERAL_BIN_MM))
    lgrid = ScoringGrid((-R1D2_LATERAL_HALF_MM, -R1D2_LATERAL_HALF_MM, z0),
                        (R1D2_LATERAL_BIN_MM, R1D2_LATERAL_BIN_MM, dz), (nl, nl, 1), name="dose")  # fmt: skip
    rms: dict[bool, NDArray[np.float64]] = {}
    for el in (False, True):
        cfg = c.wcfg_c("warp-cpu", "float64", energy=e, n=n, seed=seed_c("r1-d2", k + 2), geometry=geo,
                       grid=lgrid, n_batches=v5b.V5_BATCHES, timeout=timeout)  # fmt: skip
        res = Simulation(cfg if el else elastic_off(cfg)).run()
        clean &= bool(v5.clean(res) and res.energy_balance.relative_residual <= 1e-12)
        rms[el] = lateral_rms_batches(res.grid("dose").batch_energy_mev, R1D2_LATERAL_HALF_MM,
                                      R1D2_LATERAL_BIN_MM)
    paired = {"peak_over_plateau": paired_z(m[True]["peak_over_plateau"] - m[False]["peak_over_plateau"]),
              "plateau": paired_z(m[True]["plateau"] - m[False]["plateau"]),
              "lateral_rms": paired_z(rms[True] - rms[False])}  # fmt: skip
    total = paired_z(m[True]["total"] - m[False]["total"])
    j = r1d2_judgement(paired, total, e)
    return {"energy": e, "seeds": {"idd": seed, "lateral": seed_c("r1-d2", k + 2)}, "r80_off_mm": r80,
            "slab_mm": [z0, z0 + dz], "means": {"elastic": {"peak_over_plateau": float(m[True]["peak_over_plateau"].mean()),
            "plateau": float(m[True]["plateau"].mean()), "lateral_rms_mm": float(rms[True].mean())},
            "no_elastic": {"peak_over_plateau": float(m[False]["peak_over_plateau"].mean()),
            "plateau": float(m[False]["plateau"].mean()), "lateral_rms_mm": float(rms[False].mean())}},
            "judgement": j, "valid": bool(clean)}  # fmt: skip


def step_r1_d2(a: argparse.Namespace) -> int:
    a16 = v5._captured(v5.step_r1, argparse.Namespace(scale=a.scale, timeout=a.timeout))
    # (B) Sigma_el = 0 is bit-identical to elastic=False, and the elastic=False python nuclear-on digest
    # with the primary-only toggle equals the lv5b record (3edfce9c = b84fdf38)
    t0 = time.perf_counter()
    off = r1_nuc_digest_c(False)
    with scaled_elastic_rows(0.0):
        zero = r1_nuc_digest_c(True)
    import ionmc.transport.reference as reference

    old = reference.SECONDARY_NUCLEAR
    reference.SECONDARY_NUCLEAR = False
    try:
        toggle = r1_nuc_digest_c(False)
    finally:
        reference.SECONDARY_NUCLEAR = old
    wall_b = time.perf_counter() - t0
    b_ok = bool(zero["sha256"] == off["sha256"] and toggle["sha256"] == v5b.R1_NUC_DIGEST
                and off["valid"] and zero["valid"] and toggle["valid"])  # fmt: skip
    n = v5b.scaled(v5b.V5_N, a.scale, 2000, v5b.V5_BATCHES)
    runs = [r1d2_measure(e, k, n, a.timeout) for k, e in enumerate(v5b.V5_ENERGIES)]
    c_ok = bool(all(r["judgement"]["all_match"] and r["valid"] for r in runs))
    seeds = {f"{r['energy']:g}": list(r["seeds"].values()) for r in runs}
    guard = seed_guard({f"r1-d2-{e}": v for e, v in seeds.items()})
    a_ok = bool(a16["pass"] and run_suite.A16_INTENDED_CHANGE is None)
    doc = {
        "step": "r1-d2", "plan_row": R1D2_PLAN_ROW, "table": v5.table_record(), **c.elastic_table_record(),
        "A": {"pass": a_ok, "a16_baseline": a16["a16_baseline"], "digest_sha256": a16["digest_sha256"],
              "intended_change_set": a16["intended_change_set"]},
        "B": {"pass": b_ok, "sigma_el_zero_digest": zero["sha256"], "elastic_false_digest": off["sha256"],
              "primary_only_toggle_digest": toggle["sha256"], "recorded": v5b.R1_NUC_DIGEST,
              "baseline_commit": v5b.R1_NUC_BASELINE, "wall_s": wall_b,
              "digest_record": "r1_nuc_digest_c reproduces 3edfce9c... = b84fdf38. The frozen lv5b document "
              "also hashes EnergyBalance.elastic_mev (added by C3), so steps_v5b.r1_nuc_digest no longer "
              "reproduces 3edfce9c at HEAD. Field-by-field comparison of the python nuclear-on run "
              "(150 MeV, 400 histories, seed 20351004, SECONDARY_NUCLEAR off) between the source tree of "
              "7aae5bb and this tree: batch_energy, channels, counters, valid, nuclear diagnostics and every "
              "balance field are identical; the only difference is the new key elastic_mev ({} here). "
              "Elastic-only additions (elastic_mev when empty or zero, the appended elastic diagnostic "
              "counters) are left out of the compared document",
              "u3_slot": "free (C3 found the candidate-block slot free; the branch 'slot not free' of the "
              "row is not taken)"},
        "C": {"pass": c_ok, "rule": "paired common-random-number runs, warp-cpu float64, 150 and 200 MeV, "
              "elastic on versus off on the V5 geometry (peak/plateau, plateau, total deposit) and on a "
              "1 mm lateral grid of half-width 50 mm (dose-weighted rms radius in [0.45, 0.55] R80 of the "
              "elastic-off run); a direction matches iff the paired batch-mean difference has the declared "
              "sign with |z| >= 3 (20 batches); |delta total deposit| <= 0.1 % of E0 (90 % interval "
              "reported). The comparison with the lv5b V5 partials at 7aae5bb1 (independent seeds) is "
              "report-only and is made in the dossier, not in this step (Amendment 17 (b)1)",
              "declared_directions": {"peak_over_plateau": "down", "plateau": "up", "lateral_rms": "up"},
              "runs": runs, "expected_magnitudes": {"peak_percent": [-11, -15], "plateau_percent": [1.2, 1.4]}},
        "frozen_suite_pin": "steps_v5.nuc_config carries elastic=False (C3b); lv5c/hr5c use "
        "steps_v5c.nuc_config_c (elastic on, table " + c.ELASTIC_TABLE_ID[:8] + ")",
        "seeds": guard, "pass": bool(a_ok and b_ok and c_ok)}  # fmt: skip
    doc["status_partial"] = status_partial(a, "r1-d2-status", {"step": "r1-d2", "pass": doc["pass"],
                                           "reduced": a.scale < 1.0})
    return finish(doc, 256 + 4 * v5b.V5_N, n * 4 + a16["histories"], bool(a16["reduced"]) or n < v5b.V5_N)


def r1_nuc_digest_c(elastic: bool) -> dict[str, Any]:
    """The python nuclear-on digest of ``v5b.r1_nuc_digest`` (150 MeV, ``R1_NUC_N`` histories, seed 20351004),
    with the configuration of ``steps_v5.nuc_config`` (``elastic=False``) or ``nuc_config_c`` (``elastic=True``).
    ``EnergyBalance.elastic_mev`` (added by C3; empty or all zero unless elastic events occur) is left out
    of the hashed document when it is empty or zero: the lv5b record 3edfce9c... predates the field, and with it dropped the
    document of an elastic=False run is the b84fdf38 one byte for byte (checked against the tree of 7aae5bb)."""
    import dataclasses
    import hashlib

    geo, grid, _, _ = v5.depth_box(M.WATER, 150.0)
    build = c.nuc_config_c if elastic else v5.nuc_config
    cfg = build(energy=150.0, n=v5b.R1_NUC_N, seed=v5b.R1_NUC_SEED, geometry=geo, grid=grid, n_batches=2,
                tallies=(*v5.PRIMARY, *v5.EDEP))  # fmt: skip
    res = Simulation(cfg).run()
    arrays = {"batch_energy": res.grid("dose").batch_energy_mev,
              "channels": res.channel_raw.acc if res.channel_raw is not None else np.zeros(0)}  # fmt: skip
    bal = dataclasses.asdict(res.energy_balance)
    if all(v == 0 for v in bal.get("elastic_mev", {}).values()):  # absent, empty or all zero (Sigma_el = 0)
        bal.pop("elastic_mev", None)
    doc = {"arrays": {k: hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest()
                      for k, v in arrays.items()},
           "counters": res.counters.as_dict(), "valid": bool(res.valid), "balance": bal,
           "nuclear": res.diagnostics.get("nuclear", {})}  # fmt: skip
    from ionmc.transport.tally import ELASTIC_COUNTER_NAMES

    # the elastic additions are not part of the identity claim: the diagnostic counters appended by the
    # elastic channel (they are zero or diagnostics only) are left out
    doc["counters"] = {k: v for k, v in doc["counters"].items() if k not in ELASTIC_COUNTER_NAMES}
    text = json.dumps(doc, sort_keys=True, default=base._json)
    return {"sha256": hashlib.sha256(text.encode()).hexdigest(), "valid": bool(res.valid), "doc": doc}


# -- v8-lv (r2) ------------------------------------------------------------------------------------
def v8_pair_c(e: float, n: int, seed: int, timeout: float | None) -> dict[str, Any]:
    """``v5b.v8_pair`` with elastic on and the nuclear_local / elastic_recoil channels compared."""
    from ionmc.transport.reference import run_reference_range
    from ionmc.transport.warp_driver import run_warp_range

    geo, grid, _, _ = v5.depth_box(M.WATER, e)
    effs = {}
    for backend in ("python", "warp-cpu"):
        cfg = c.wcfg_c(backend, "float64", energy=e, n=n, seed=seed, geometry=geo, grid=grid, n_batches=2,
                       tallies=(NUC_LOCAL, REC), timeout=timeout)  # fmt: skip
        cfg = replace(cfg, diagnostics=DiagnosticsOptions(track_end_positions=True, trace_histories=n))
        effs[backend] = Simulation(cfg).effective
    t0 = time.perf_counter()
    py = run_reference_range(effs["python"], 0, n)
    t_py = time.perf_counter() - t0
    t0 = time.perf_counter()
    wp_ = run_warp_range(effs["warp-cpu"], 0, n, "cpu")
    t_wp = time.perf_counter() - t0
    m = v5b.trace_compare(py, wp_)
    chan = float(np.max(np.abs(np.asarray(py.channel_acc, dtype=np.int64)
                               - np.asarray(wp_.channel_acc, dtype=np.int64)))) * v5b.QUANTUM_MEV
    from ionmc.transport.tally import tally_column_index

    def count(part: Any, name: str) -> float:
        return math.fsum(part.tally_components[tally_column_index(len(part.tally_components), True, name)])

    m.update(python_s=t_py, warp_cpu_s=t_wp, energy_mev=e, histories=n, seed=seed,
             channels_max_abs_mev=chan, elastic_pp_events=count(py, "elastic_pp_events"),
             elastic_pa_events=count(py, "elastic_pa_events"),
             elastic_recoil_local_mev=count(py, "elastic_recoil_local"),
             nuclear_device_sha256=wp_.meta.get("nuclear_device_sha256"))  # fmt: skip
    m["n_elastic_events"] = m["elastic_pp_events"] + m["elastic_pa_events"]
    m["max_continuous"] = max(m["max_continuous"], chan)
    return m


def step_v8_lv(a: argparse.Namespace) -> int:
    k = v5b.scaled(v5b.V8_K, a.scale, 8, 2)
    nd = v5b.scaled(v5b.V8_DENSE_N, a.scale, 16, 2)
    ci = v8_pair_c(150.0, k, seed_c("v8-lv", 0), a.timeout)
    with v5b.scaled_rows(v5b.V8_DENSE_FACTOR), scaled_elastic_rows(v5b.V8_DENSE_FACTOR):
        dense = v8_pair_c(100.0, nd, seed_c("v8-lv", 1), a.timeout)
    need = v5b.V8_DENSE_MIN_EVENTS if a.scale == 1.0 else max(1, nd // 2)
    dense["events_required"] = need

    def row_ok(m: dict[str, Any]) -> bool:
        return bool(m.get("discrete_identical") and m.get("max_continuous", 1.0) <= v5b.V8_TOL)

    ok = row_ok(ci) and row_ok(dense) and dense.get("n_events", 0) >= need and dense["n_elastic_events"] > 0
    guard = seed_guard({"v8-lv": [seed_c("v8-lv", 0), seed_c("v8-lv", 1)]})
    doc = {"step": "v8-lv", "table": v5.table_record(), **c.elastic_table_record(),
           "tolerance": v5b.V8_TOL, "ci": ci, "nuclear_dense": {**dense, "sigma_factor": v5b.V8_DENSE_FACTOR},
           "counters": "python and warp counter vectors (incl. elastic_below_domain / pp_below_domain) are "
           "compared inside trace_compare (discrete_identical)", "seeds": guard, "pass": bool(ok)}  # fmt: skip
    doc["status_partial"] = status_partial(a, "v8-lv-status", {"step": "v8-lv", "pass": bool(ok),
                                           "reduced": a.scale < 1.0})
    return finish(doc, v5b.V8_K + v5b.V8_DENSE_N, k + nd, a.scale < 1.0)


# -- p7-sampler (r15) ------------------------------------------------------------------------------
P7_CASES = tuple([("H-1", e) for e in (20.0, 50.0, 100.0, 150.0, 200.0)]
                 + [(t, e) for t in ("C-12", "O-16") for e in (20.0, 50.0, 100.0, 150.0, 200.0, 250.0)])
P7_N, P7_BINS, P7_P_MIN, P7_SIGMA_RTOL = 100_000, 64, 0.001, 1.0e-3
P7_PLAN_ROW = (
    "P7 (plan Amendment 14 (d), l.197): Sampler vs the built tables, offline: H-1 at 20 / 50 / 100 / 150 / "
    "200 MeV; C-12 and O-16 at 20 / 50 / 100 / 150 / 200 / 250 MeV; 1e5 events per case; mu_CM histogram "
    "in >= 50 equiprobable bins of the table CDF; sigma_el table vs direct builder evaluation at every "
    "node and node midpoint | chi2 p > 0.001 in every case; sigma relative error <= 1e-3 at every point"
)


def p7_sigma_check(tab: Any) -> dict[str, Any]:
    """sigma_el of the table against the direct builder evaluation at every node and node midpoint of every
    target, inside the declared domain (below it sigma must be exactly 0). p + A: the BGG elastic cross
    section with the AME2020 target mass; H-1: the LA150 reconstruction up to the 150 MeV anchor and the
    anchor times the BGG p-p ratio above it (``nuclear_v10.sigma_pp_direct_mb``)."""
    import nuclear_v10 as v10
    from ionmc.nuclear import bgg
    from ionmc.nuclear.build import TARGETS

    grid = np.asarray(tab.arrays["grid_e_mev"], dtype=float)
    dom = tab.elastic_domain()
    masses = v10.target_masses_mev()
    sec, awi = v10.h1_law5()
    z = {t.name: t.z for t in TARGETS}
    out: dict[str, Any] = {}
    worst = 0.0
    for name in tab.target_names:
        it = tab.target_names.index(name)
        e_min = dom[name][0]
        sig = np.asarray(tab.arrays["sigma_barn"][it], dtype=float) * 1.0e3  # mb

        def direct(e: float, name: str = name) -> float:
            if name == "H-1":
                return v10.sigma_pp_direct_mb(sec, awi, e)
            return float(bgg.bgg_elastic_mb(e, z[name], masses[name]))

        pts = [(float(e), float(s)) for e, s in zip(grid, sig, strict=True)]
        # midpoints entirely inside the domain are evaluated; entirely below it must be zero; the one that
        # straddles the domain edge has no direct value and is skipped
        mids = [(float(0.5 * (grid[k] + grid[k + 1])), float(0.5 * (sig[k] + sig[k + 1])))
                for k in range(grid.size - 1) if grid[k] >= e_min - 1e-12 or grid[k + 1] < e_min - 1e-12]
        worst_t, n_pts, below_ok = 0.0, 0, True
        for items in (pts, mids):
            for e, s in items:
                if e < e_min - 1e-12:
                    below_ok &= bool(s == 0.0)
                    continue
                d = direct(e)
                worst_t = max(worst_t, abs(s / d - 1.0)) if d > 0 else math.inf
                n_pts += 1
        out[name] = {"max_rel_error": worst_t, "points": n_pts, "zero_below_domain": below_ok,
                     "domain_mev": list(dom[name]), "pass": bool(worst_t <= P7_SIGMA_RTOL and below_ok)}
        worst = max(worst, worst_t)
    return {"targets": out, "max_rel_error": worst, "rtol": P7_SIGMA_RTOL,
            "pass": bool(all(v["pass"] for v in out.values()))}  # fmt: skip


def p7_case(tab: Any, name: str, e: float, seed: int, n: int) -> dict[str, Any]:
    from ionmc.nuclear.elastic_tables import chi2_equiprobable

    t = tab.target_names.index(name)
    rng = np.random.default_rng(seed)
    mu = tab.sample_mu(t, e, rng.random(n))
    chi2, p, dof = chi2_equiprobable(mu, tab.edges_at(t, e), P7_BINS)
    return {"target": name, "energy_mev": e, "seed": seed, "events": n, "bins": P7_BINS, "dof": dof,
            "chi2": chi2, "p": p, "pass": bool(p > P7_P_MIN)}  # fmt: skip


def step_p7(a: argparse.Namespace) -> int:
    from ionmc.nuclear.elastic_tables import ElasticTable

    tab = ElasticTable.load(None, c.ELASTIC_TABLE_ID)
    n = v4.scaled(P7_N, a.scale, 5000)
    cases = [p7_case(tab, nm, e, seed_c("p7", k), n) for k, (nm, e) in enumerate(P7_CASES)]
    sig = p7_sigma_check(tab)
    guard = seed_guard({"p7": [seed_c("p7", k) for k in range(len(P7_CASES))]})
    ok = bool(all(x["pass"] for x in cases) and sig["pass"])
    doc = {"step": "p7-sampler", "plan_row": P7_PLAN_ROW, "table": v5.table_record(),
           **c.elastic_table_record(), "cases": cases, "n_cases": len(cases),
           "sampling": "numpy default_rng(seed = 20496004 + case index) uniforms through ElasticTable.sample_mu "
           "(the reference inverse CDF; the Warp samplers are P7-tested in CI)",
           "sigma_el": sig, "family_error_note": "with 17 cases a correct sampler fails at least one with "
           "probability 1.7 % (plan)", "seeds": guard, "pass": ok}  # fmt: skip
    doc["status_partial"] = status_partial(a, "p7-status", {"step": "p7-sampler", "pass": ok,
                                           "reduced": n < P7_N})
    return finish(doc, P7_N, n, n < P7_N)


# -- V5 re-run (r4-7) ------------------------------------------------------------------------------
def v5_seed(energy: float, nuclear: bool) -> int:
    return seed_c(f"v5-{energy:g}-{'on' if nuclear else 'off'}")


def step_v5_ionmc(a: argparse.Namespace) -> int:
    """V5 producer of lv5c: warp-cpu float64, 1e5 histories in 20 batches, V5 geometry. "on" = ``nuc_config_c``
    (nuclear + elastic) with the recoil depth channel; "off" = ``nuclear=False`` (the lv5b physics)."""
    from ionmc.reference import metrics as rm

    e, nuclear = a.energy, a.nuclear == "on"
    row = f"v5-{e:g}-{a.nuclear}"
    n = v5b.scaled(v5b.V5_N, a.scale, 2000, v5b.V5_BATCHES)
    geo, grid, nz, r_mm = v5b.v5_geometry(e)
    seed = v5_seed(e, nuclear)
    cfg = c.wcfg_c("warp-cpu", "float64", energy=e, n=n, seed=seed, geometry=geo, grid=grid,
                   nuclear=nuclear, n_batches=v5b.V5_BATCHES, tallies=(REC,) if nuclear else (),
                   timeout=a.timeout)  # fmt: skip
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    g = res.grid("dose")
    rho = M.WATER.density_g_cm3
    bin_g = rho * v5b.V5_DZ_MM / 10.0
    idd_b = g.batch_energy_mev.reshape(v5b.V5_BATCHES, nz) / bin_g
    idd = idd_b.mean(axis=0)
    depth = (np.arange(nz) + 0.5) * v5b.V5_DZ_MM
    curve = rm.normalize_to_peak(idd)
    plateau = (depth >= 20.0) & (depth <= 60.0)
    tot_b = idd_b.sum(axis=1) * bin_g
    met = {"peak_depth_mm": rm.peak_depth(depth, curve), "r80_mm": rm.r80(depth, curve),
           "plateau_mean_20_60_mm": float(idd[plateau].mean()),
           "peak_over_plateau": float(idd.max() / idd[plateau].mean()),
           "total_in_grid_mev_per_primary": float(tot_b.mean()),
           "total_in_grid_sem": float(tot_b.std(ddof=1) / math.sqrt(v5b.V5_BATCHES)),
           "in_grid_over_e0": float(tot_b.mean() / e)}  # fmt: skip
    ok = bool(v5.clean(res) and res.energy_balance.relative_residual <= 1e-12)
    extra: dict[str, Any] = {}
    if nuclear:
        rec_b = v5b.tally_batches(res, "recoil") / bin_g  # MeV/(g/cm^2)/primary, same units as idd
        extra = {"recoil_idd_batches": rec_b.tolist(), "recoil_unit": "MeV/(g/cm^2)/primary",
                 "recoil_total_mev_per_primary": float((rec_b.sum(axis=1) * bin_g).mean()),
                 "elastic_recoil_local_mev": float(res.energy_balance.nuclear_mev.get("elastic_recoil_local", 0.0))}
    part = c.write_partial_c(a, row, {
        "row": row, "energy": e, "nuclear": nuclear, "seed": seed, "n": n, "n_batches": v5b.V5_BATCHES,
        "bin_mm": v5b.V5_DZ_MM, "lateral_half_mm": v5b.V5_HALF_MM, "density_g_cm3": rho,
        "unit": "MeV/(g/cm^2)/primary", "idd_batches": idd_b.tolist(), "metrics": met, "valid": ok,
        "reduced": n < v5b.V5_N, "nuclear_mev": dict(res.energy_balance.nuclear_mev) if nuclear else {},
        **extra})  # fmt: skip
    doc = {"step": "v5-ionmc", "row": row, "table": v5.table_record() if nuclear else None,
           **c.elastic_table_record(), "seed": seed, "seeds": seed_guard({row: [seed]}),
           "range_csda_mm": r_mm, "wall_s": wall, "hist_per_s": n / wall, "metrics": met,
           "relative_residual": float(res.energy_balance.relative_residual),
           "counters": res.counters.as_dict(), "partial": part,
           "comparison": "pending: step v5-compare / v5-attribution", "pass": ok}  # fmt: skip
    return finish(doc, v5b.V5_N, n, n < v5b.V5_N)


V5_NAMES = tuple(f"v5-{e}-{t}.json" for e in (150, 200) for t in ("on", "off"))
V5_COMPARE_PLAN = (
    "V5 re-run (plan Amendment 14 (f), l.254): V5 is re-run at lv5c r_index 4-7 with the frozen tolerances, "
    "frozen reference runs and frozen comparator. TOPAS gates and MCsquare is report-only"
)


def _v5_inputs(a: argparse.Namespace, extra: tuple[str, ...] = ()) -> tuple[list[dict[str, Any]], Any, Path, Path]:
    cmp_ = _cmp()
    parts = [p for p, _ in c.load_partials_c(a, list(V5_NAMES) + list(extra))]
    ref_dir = Path(a.reference_dir) if a.reference_dir else v5b.REPO / ".ionmc-cache" / "reference-runs"
    return parts, cmp_, ref_dir, c._cases_dir()


def v5_verdict(a: argparse.Namespace) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """The frozen-comparator V5 verdict from the four lv5c partials (``SystemExit`` / ``IddError`` on a lineage
    or input error); also checks that the four seed sets are disjoint and unconsumed."""
    import tempfile

    parts, cmp_, ref_dir, cases_dir = _v5_inputs(a)
    seeds = {p["row"]: [p["seed"]] for p in parts}
    for p in parts:
        if p["seed"] != v5_seed(p["energy"], p["nuclear"]) or p.get("seed_base") != base.SEED_BASE:
            raise SystemExit(f"{p['row']}: seed {p['seed']} differs from the r_index rule")
    seed_guard(seeds)
    topas, mc = v5b._v5_reference_runs(ref_dir)
    with tempfile.TemporaryDirectory() as tmp:
        for nm, p in zip(V5_NAMES, parts, strict=True):
            (Path(tmp) / nm).write_text(json.dumps(p, sort_keys=True))
        verdict = cmp_.build_verdict(Path(tmp), topas, mc, cases_dir)
    verdict.update(c.v11_bound_cases_identity(verdict, cases_dir))
    return verdict, parts


def step_v5_compare(a: argparse.Namespace) -> int:
    doc: dict[str, Any] = {"step": "v5-compare", "plan_rule": V5_COMPARE_PLAN}
    reduced = False
    try:
        verdict, parts = v5_verdict(a)
        reduced = any(p["reduced"] for p in parts)
        doc.update(verdict=verdict, pass_gating=bool(verdict["pass"]), cases_in_source_identity=True, error=None,
                   **c.elastic_table_record())  # fmt: skip
        doc["pass"] = bool(verdict["pass"] and not reduced)
    except (_cmp().IddError, SystemExit, ValueError, OSError) as exc:
        doc.update(cases_in_source_identity=False, verdict=None, error=f"{type(exc).__name__}: {exc}",
                   **{"pass": False})  # fmt: skip
    doc["attestation"] = v5.attestation_block(a)  # after ALL loads: imported partials are listed (never omitted)
    return finish(doc, v5b.V5_N, v5b.V5_N, reduced)


# -- attribution ladder (Amendment 14 (f), Amendment 17 (b)4/(b)5) ---------------------------------
ATTRIBUTION_PLAN = (
    "Attribution ladder (plan Amendment 14 (f), l.260-270): It applies only if V5 fails. An attribution is "
    "recorded as an attribution, and V5 stays failed. Rung 0 (implementation): P6, P7, P5-ext, R1-D2 and "
    "V8-LV must pass. Rung 1 (elastic-model difference): V11 at the failing energy; if V11 fails in the "
    "direction of the V5 residual, the failure is attributed to the p+A angular shape. Rung 2 (V11 passes, "
    "ionmc above TOPAS by 2-5 % on peak/plateau): the non-elastic factor F_ne, attributed iff the 90 % "
    "interval of F_ne,ionmc/F_ne,TOPAS - 1 lies beyond 0.02 on the side of the residual. Rung 3 (nuclear-"
    "recoil local deposition), only if rungs 1-2 exclude: attributed iff the recorded recoil energy "
    "deposited in [peak - 10 mm, R80] is >= 50 % of the residual's energy equivalent; otherwise "
    "unexplained. Stop rule: if peak/plateau fails low (ionmc below TOPAS by more than 2 %) or by more "
    "than 5 % in either direction, the task stops: a recorded scientific failure, escalated to the "
    "orchestrator, with no tolerance change. (Amendment 17 (b)4: the document is written, the step exits "
    "non-zero, the remaining steps still run, no fresh base is taken.)"
)
RUNG3_FRACTION = 0.5
CI_RUNG0 = ("P6", "P6-D", "P5-ext")  # CI rows, evaluated by the merge gate, not inside the archive


def attribution_ladder(
    v5: dict[str, dict[str, Any]], rung0: dict[str, dict[str, Any]], v11: dict[str, dict[str, Any]],
    fne: dict[str, dict[str, Any]], rung3: dict[str, dict[str, Any]],
) -> dict[str, Any]:  # fmt: skip
    """Pure ladder. Per energy key: ``v5[e] = {"pass": bool, "pp_diff": float, "failing": [criteria]}``;
    ``rung0[name] = {"status": "pass"|"fail"|"ci"}`` (``"ci"``: CI gate, not evaluated in the archive);
    ``v11[e] = {"pass": bool, "F_diff": float}``; ``fne[e] = {"computed": bool, "attributed_if_residual_above",
    "attributed_if_residual_below", "ratio_minus_1": {...}}``; ``rung3[e] = {"computed": bool, "e_res_mev",
    "e_rec_mev", ...}``. The residual side is the sign of the V5 peak/plateau difference (ionmc/TOPAS - 1)."""
    v5_fails = any(not x["pass"] for x in v5.values())
    doc: dict[str, Any] = {"applies": bool(v5_fails), "plan_rule": ATTRIBUTION_PLAN,
                           "expectation_preregistered": EXPECTATION, "energies": {}}  # fmt: skip
    doc["v5_stays_failed"] = bool(v5_fails)
    stop = []
    for e, x in v5.items():
        d = x["pp_diff"]
        if v5_fails and (d < STOP_LOW or abs(d) > STOP_ABS):
            stop.append({"energy": e, "pp_diff": d,
                         "reason": "fails low (< -2 %)" if d < STOP_LOW else "by more than 5 %"})
    doc["stop_rule_triggered"] = bool(stop)
    doc["stop_rule"] = stop
    r0_fail = sorted(k for k, v in rung0.items() if v["status"] == "fail")
    doc["rung0"] = {"items": rung0, "failed": r0_fail,
                    "ci_gates_not_evaluated_here": sorted(k for k, v in rung0.items() if v["status"] == "ci"),
                    "pass": not r0_fail}  # fmt: skip
    if not v5_fails:
        doc["attribution"] = "not applicable (V5 passes)"
    else:
        for e, x in v5.items():
            if x["pass"]:
                doc["energies"][e] = {"attribution": "not applicable (V5 passes at this energy)"}
                continue
            side = 1 if x["pp_diff"] > 0 else -1 if x["pp_diff"] < 0 else 0
            rec: dict[str, Any] = {"failing": x["failing"], "pp_diff": x["pp_diff"], "residual_side": side}
            if r0_fail:
                rec.update(attribution="implementation defect (rung 0)", failed=r0_fail,
                           action="fix, review, one fresh base")
            elif "peak_over_plateau" not in x["failing"]:
                rec.update(attribution="unexplained: V5 fails on " + ", ".join(x["failing"]) +
                           ", not on peak/plateau; rungs 1-3 address the peak/plateau residual")
            else:
                v = v11[e]
                v11_dir = 1 if v["F_diff"] > 0 else -1
                rec["rung1"] = {"v11_pass": v["pass"], "v11_F_diff": v["F_diff"],
                                "v11_fails_in_residual_direction": bool(not v["pass"] and v11_dir == side)}  # fmt: skip
                if rec["rung1"]["v11_fails_in_residual_direction"]:
                    rec.update(attribution="p+A angular shape (rung 1; radius rule of D-1)",
                               note="recorded scientific failure; (ii') only through its trigger in (e)")
                else:
                    r2_applies = bool(v["pass"] and x["pp_diff"] > RUNG2_LOW)
                    rec["rung2"] = {"applies": r2_applies}
                    if r2_applies:
                        f = fne[e]
                        hit = bool(f.get("computed") and f["attributed_if_residual_above"])
                        rec["rung2"].update(fne=f, attributed=hit)
                    if r2_applies and rec["rung2"]["attributed"]:
                        rec.update(attribution="inclusive-yield energy excess / alpha local deposition (rung 2)")
                    else:
                        r3 = rung3[e]
                        ok3 = bool(r3.get("computed") and r3["e_res_mev"] > 0
                                   and r3["e_rec_mev"] >= RUNG3_FRACTION * r3["e_res_mev"])
                        rec["rung3"] = {**r3, "fraction_required": RUNG3_FRACTION, "attributed": ok3}
                        rec.update(attribution="nuclear-recoil local deposition (rung 3)" if ok3 else "unexplained")
            doc["energies"][e] = rec
        att = {r["attribution"] for r in doc["energies"].values()}
        doc["attribution"] = sorted(att)
    if doc["stop_rule_triggered"]:
        doc["message"] = ("recorded scientific failure, escalate to the orchestrator, no tolerance change; "
                          "the remaining steps still run (Amendment 17 (b)4), no fresh base for V5")
    doc["pass"] = not doc["stop_rule_triggered"]  # the step exits non-zero iff the stop rule fired
    return doc


def window_energy(depth: NDArray[np.float64], lo: float, hi: float, curve: NDArray[np.float64],
                  bin_g: float) -> float:
    """Energy [MeV per primary] of an IDD curve [MeV/(g/cm^2)/primary] over the depth window [lo, hi]."""
    sel = (depth >= lo) & (depth <= hi)
    return float(curve[sel].sum() * bin_g)


def step_v5_attribution(a: argparse.Namespace) -> int:
    import tempfile

    from ionmc.reference import metrics as rm

    doc: dict[str, Any] = {"step": "v5-attribution", "report_only": True}
    cmp_ = _cmp()
    try:
        verdict, parts = v5_verdict(a)
        stat = {k: c.load_partials_c(a, [f"{k}-status.json"])[0][0] for k in ("p7", "r1-d2", "v8-lv")}
        rung0: dict[str, dict[str, Any]] = {
            "P7": {"status": "pass" if stat["p7"]["pass"] else "fail"},
            "R1-D2": {"status": "pass" if stat["r1-d2"]["pass"] else "fail"},
            "V8-LV": {"status": "pass" if stat["v8-lv"]["pass"] else "fail"},
            **{k: {"status": "ci"} for k in CI_RUNG0}}  # fmt: skip
        v11_parts = [p for p, _ in c.load_partials_c(a, list(c.V11_NAMES))]
        c.v11_seed_check(v11_parts)
        ref_dir = Path(a.reference_dir) if a.reference_dir else v5b.REPO / ".ionmc-cache" / "reference-runs"
        cases_dir = c._cases_dir()
        with tempfile.TemporaryDirectory() as tmp:
            for p in [*v11_parts, *parts]:
                (Path(tmp) / f"{p['row']}.json").write_text(json.dumps(p, sort_keys=True))
            v11v = cmp_.build_v11_verdict(Path(tmp), c._v11_reference_runs(ref_dir), cases_dir)
        topas, _ = v5b._v5_reference_runs(ref_dir)
        by = {p["row"]: p for p in parts}
        v5in: dict[str, Any] = {}
        v11in: dict[str, Any] = {}
        fnein: dict[str, Any] = {}
        r3in: dict[str, Any] = {}
        nzs = {e: np.asarray(by[f"v5-{e}-on"]["idd_batches"]).shape[1] for e in (150, 200)}
        groups = cmp_.load_engine_groups(topas, "topas", nzs, cases_dir)
        for e in (150, 200):
            rows = verdict["engines"]["topas"]["energies"][str(e)]
            fail = sorted(k for k in cmp_.TOLERANCES if not rows[k]["pass"])
            v5in[str(e)] = {"pass": bool(rows["_all_pass"]), "pp_diff": float(rows["peak_over_plateau"]["diff"]),
                            "failing": fail}  # fmt: skip
            v11in[str(e)] = {"pass": bool(v11v["energies"][str(e)]["pass"]),
                             "F_diff": float(v11v["energies"][str(e)]["F"]["diff"])}  # fmt: skip
            fnein[str(e)] = v11v["f_ne"][str(e)]
            on = by[f"v5-{e}-on"]
            idd_on = np.asarray(on["idd_batches"]).mean(axis=0)
            rec_on = np.asarray(on["recoil_idd_batches"]).mean(axis=0)
            top = groups[(e, True)]["curves"].mean(axis=0)
            depth = (np.arange(idd_on.size) + 0.5) * v5b.V5_DZ_MM
            bin_g = float(on["density_g_cm3"]) * v5b.V5_DZ_MM / 10.0
            r80 = float(rm.r80(depth, rm.normalize_to_peak(top)))
            peak = float(rm.peak_depth(depth, rm.normalize_to_peak(top)))
            lo = peak - 10.0
            e_res = window_energy(depth, lo, r80, idd_on - top, bin_g)
            e_rec = window_energy(depth, lo, r80, rec_on, bin_g)
            r3in[str(e)] = {"computed": True, "window_mm": [lo, r80], "peak_topas_mm": peak,
                            "e_res_mev": e_res, "e_rec_mev": e_rec, "e_rec_over_e_res": e_rec / e_res if e_res else None,
                            "residual_energy_definition": "E_res = sum over the window [TOPAS peak - 10 mm, TOPAS R80] of "
                            "(IDD_ionmc,on - IDD_TOPAS,full) x bin mass thickness, per primary [MeV]; definition "
                            "pre-registered in Amendment 14a (plan gives none); rung 3 is considered only if "
                            "rungs 1-2 exclude",
                            "recoil_definition": "sum over the window of the elastic-recoil depth channel x bin mass thickness"}  # fmt: skip
        lad = attribution_ladder(v5in, rung0, v11in, fnein, r3in)
        doc.update(lad, v5_pass=bool(verdict["pass"]), v11_pass=bool(v11v["pass"]),
                   rung0_sources={k: v["content_sha256"] for k, v in stat.items()},
                   **c.elastic_table_record())  # fmt: skip
    except (cmp_.IddError, SystemExit, ValueError, OSError, KeyError) as exc:
        doc.update(error=f"{type(exc).__name__}: {exc}", stop_rule_triggered=False, **{"pass": False})
    doc["attestation"] = v5.attestation_block(a)  # after ALL loads: imported partials are listed (never omitted)
    return finish(doc, v5b.V5_N, v5b.V5_N, False)


# -- elastic-on versions of the frozen lv5b producers (config via nuc_config_c) ----------------------
def v7_config_c(backend: str, precision: str, n: int, nb: int, seed: int, *, grid: ScoringGrid | None = None,
                geo: Any = None, timeout: float | None = None) -> SimulationConfig:
    g0, gr0, _, _ = v5.depth_box(M.WATER, 150.0)
    return c.wcfg_c(backend, precision, energy=150.0, n=n, seed=seed, geometry=geo or g0, grid=grid or gr0,
                    n_batches=nb, tallies=v5b.V7_TALLIES, timeout=timeout)  # fmt: skip


def est_ok(est: dict[str, Any]) -> bool:
    """Validity of ``batch_estimates`` output with elastic on: every counter but the two diagnostics is zero."""
    return bool(c.invalid_counters(est) == 0 and est["batch_assignment_ok"])


# -- V2b-el (r8) -----------------------------------------------------------------------------------
V2B_EL_PLAN = ("V2b-el (plan Amendment 14 (d)): Frozen row V2b re-executed with the elastic channel on "
               "(warp-cpu float64, s_max 0.1 vs 1.0 mm, 1e6 per variant, common random numbers) | As V2b "
               "(frozen, unchanged) | falsification, LV, r_index 8")


def step_v2b_el_shard(a: argparse.Namespace) -> int:
    k = a.shard
    n = v5b.scaled(v5b.V2B_N // v5b.V2B_SHARDS, a.scale, 2000, v5b.V2B_BATCHES)
    geo, grid, _, _ = v5.depth_box(M.WATER, 150.0)
    ratios, okc, wall = [], True, 0.0
    for smax, frac in v5b.V2B_VARIANTS:
        cfg = c.wcfg_c("warp-cpu", "float64", energy=150.0, n=n, seed=seed_c("v2b-el", k), geometry=geo,
                       grid=grid, mcs=False, straggling=False, max_step=smax, frac=frac,
                       tallies=v5.PRIMARY, n_batches=v5b.V2B_BATCHES)  # fmt: skip
        t0 = time.perf_counter()
        res = Simulation(cfg).run()
        wall += time.perf_counter() - t0
        okc &= v5.clean(res)
        ratios.append(v5.primary_ratio(res))
    path = c.write_partial_c(a, f"v2b-el-s{k}", {
        "row": "v2b-el", "shard": k, "seed": seed_c("v2b-el", k), "n": n, "valid": bool(okc),
        "diff": (ratios[0] - ratios[1]).tolist(), "s_b": ratios[1].tolist(),
        "reduced": n < v5b.V2B_N // v5b.V2B_SHARDS})  # fmt: skip
    doc = {"step": "v2b-el-shard", "row": "v2b-el", "shard": k, "seed": seed_c("v2b-el", k),
           "table": v5.table_record(), **c.elastic_table_record(), "wall_s": wall, "hist_per_s": 2 * n / wall,
           "partial": path, "pass": bool(okc)}  # fmt: skip
    return finish(doc, v5b.V2B_N // v5b.V2B_SHARDS, n, n < v5b.V2B_N // v5b.V2B_SHARDS)


def step_v2b_el_combine(a: argparse.Namespace) -> int:
    parts = [p for p, _ in c.load_partials_c(a, [f"v2b-el-s{k}.json" for k in range(v5b.V2B_SHARDS)])]
    for k, p in enumerate(parts):
        if p["seed"] != seed_c("v2b-el", k) or p["shard"] != k or not p["valid"]:
            raise SystemExit(f"v2b-el shard {k}: seed/shard mismatch or invalid result")
    guard = seed_guard({"v2b-el": [p["seed"] for p in parts]})
    diff = np.concatenate([np.array(p["diff"]) for p in parts])
    n = sum(p["n"] for p in parts)
    nz = diff.shape[1]
    bins = v5.depth_bins(150.0, nz)
    mean = diff.mean(axis=0)[bins]
    sem = diff.std(axis=0, ddof=1)[bins] / math.sqrt(diff.shape[0])
    conclusive = bool(sem[-1] <= v5.PROBE_SIGMA_MAX)
    ok = bool(conclusive and np.all(np.abs(mean) <= v5.PROBE_TOL))
    doc = {"step": "v2b-el-combine", "plan_row": V2B_EL_PLAN, "table": v5.table_record(),
           **c.elastic_table_record(), "seeds": guard,
           "variants": {"a": {"s_max_mm": v5b.V2B_VARIANTS[0][0], "f_e": v5b.V2B_VARIANTS[0][1]},
                        "b": {"s_max_mm": v5b.V2B_VARIANTS[1][0], "f_e": v5b.V2B_VARIANTS[1][1]}},
           "histories_per_variant": n, "batches": int(diff.shape[0]), "depths_mm": [float(b) for b in bins],
           "diff": mean.tolist(), "sem": sem.tolist(), "max_abs_diff": float(np.abs(mean).max()),
           "sigma_delta_deepest": float(sem[-1]), "sigma_max": v5.PROBE_SIGMA_MAX, "conclusive": conclusive,
           "tolerance": v5.PROBE_TOL, "pass": ok, "attestation": v5.attestation_block(a)}  # fmt: skip
    return finish(doc, v5b.V2B_N, n, any(p["reduced"] for p in parts))


# -- V7-scan-el (r9), V7-shift-el (r10) ------------------------------------------------------------
def f32_vs_f64_c(backend: str, n: int, nb: int, seed_f64: int, seed_f32: int, timeout: Any) -> dict[str, Any]:
    """``v5b.f32_vs_f64`` with the elastic configuration and the elastic-aware validity rule."""
    from ionmc.transport import parity

    runs, wall = {}, {}
    for prec, sd in (("float64", seed_f64), ("float32", seed_f32)):
        t0 = time.perf_counter()
        runs[prec] = v5b.batch_estimates(v7_config_c(backend, prec, n, nb, sd, timeout=timeout))
        wall[prec] = time.perf_counter() - t0
    ok = bool(all(est_ok(r) for r in runs.values()))
    out: dict[str, Any] = {}
    for name in v5b.ESTIMATORS:
        x, y = (v5b.scalar_estimators(runs[p])[name] for p in ("float64", "float32"))
        z = parity.scalar_z(x, y)
        out[name] = {**z, "pass": bool(abs(z["z"]) < v5b.V7_F32_Z)}
        ok &= out[name]["pass"]
    return {"backend": backend, "histories": n, "batches": nb, "seeds": [seed_f64, seed_f32], "wall_s": wall,
            "z_limit": v5b.V7_F32_Z, "estimators": out, "pass": bool(ok)}  # fmt: skip


def step_v7_scan_el(a: argparse.Namespace) -> int:
    ns = [v5b.scaled(n, a.scale, 1000 * (j + 1), v5b.V7_BATCHES) for j, n in enumerate(v5b.V7_SCAN_N)]
    rows: dict[str, Any] = {}
    series: dict[str, list[float]] = {k: [] for k in v5b.ESTIMATORS}
    ok = True
    for j, n in enumerate(ns):
        t0 = time.perf_counter()
        est = v5b.batch_estimates(v7_config_c("warp-cpu", "float64", n, v5b.V7_BATCHES, seed_c("v7-scan-el", j),
                                              timeout=a.timeout))  # fmt: skip
        row: dict[str, Any] = {"histories": n, "seed": seed_c("v7-scan-el", j), "wall_s": time.perf_counter() - t0,
                               "counters_invalid": c.invalid_counters(est), "counters_by_name": c.counters_by_name(est)}  # fmt: skip
        ok &= est_ok(est)
        for name, x in v5b.scalar_estimators(est).items():
            mean, sem, r = v5b.rel_se(x)
            row[name] = {"mean": mean, "sem": sem, "relative_se": r}
            series[name].append(r)
        rows[f"N={n}"] = row
    slopes = {}
    for name, rs in series.items():
        slope = float(np.polyfit(np.log(ns), np.log(rs), 1)[0])
        s_ok = bool(abs(slope - v5b.V7_SLOPE) <= v5b.V7_SLOPE_TOL)
        slopes[name] = {"slope": slope, "limit": [v5b.V7_SLOPE - v5b.V7_SLOPE_TOL, v5b.V7_SLOPE + v5b.V7_SLOPE_TOL],
                        "pass": s_ok}  # fmt: skip
        ok &= s_ok
    nf = v5b.scaled(v5b.V7_F32_N, a.scale, 2000, v5b.V7_BATCHES)
    f32 = f32_vs_f64_c("warp-cpu", nf, v5b.V7_BATCHES, seed_c("v7-scan-el", 3), seed_c("v7-scan-el", 4), a.timeout)
    ok &= f32["pass"]
    guard = seed_guard({"v7-scan-el": [seed_c("v7-scan-el", j) for j in range(5)]})
    doc = {"step": "v7-scan-el", "table": v5.table_record(), **c.elastic_table_record(), "seeds": guard,
           "batches": v5b.V7_BATCHES, "runs": rows, "slopes": slopes, "f32_vs_f64": f32, "pass": bool(ok)}  # fmt: skip
    return finish(doc, sum(v5b.V7_SCAN_N) + v5b.V7_F32_N, sum(ns) + nf, a.scale < 1.0)


def step_v7_shift_el(a: argparse.Namespace) -> int:
    n = v5b.scaled(v5b.V7_SHIFT_N, a.scale, 2000, v5b.V7_SHIFT_BATCHES)
    _, _, nz, _ = v5.depth_box(M.WATER, 150.0)
    tot: dict[str, NDArray[np.float64]] = {}
    ok = True
    for kind in ("base", "shifted", "refined"):
        geo, grid = v5b.lateral_layout(kind, nz)
        est = v5b.batch_estimates(v7_config_c("warp-cpu", "float64", n, v5b.V7_SHIFT_BATCHES,
                                              seed_c("v7-shift-el"), grid=grid, geo=geo, timeout=a.timeout))  # fmt: skip
        ok &= est_ok(est)
        tot[kind] = est["sec_p"].sum(axis=1)
    out = {}
    for kind, tol in (("shifted", v5b.V7_SHIFT_RTOL), ("refined", v5b.V7_REFINE_RTOL)):
        d = tot[kind] - tot["base"]
        rel = float(abs(d.mean()) / tot["base"].mean())
        out[kind] = {"tolerance": tol, "relative_mean_difference": rel,
                     "max_abs_batch_difference_mev": float(np.abs(d).max()),
                     "base_total": float(tot["base"].mean()), "pass": bool(rel <= tol)}  # fmt: skip
        ok &= out[kind]["pass"]
    doc = {"step": "v7-shift-el", "table": v5.table_record(), **c.elastic_table_record(),
           "seeds": seed_guard({"v7-shift-el": [seed_c("v7-shift-el")]}), "seed": seed_c("v7-shift-el"),
           "histories": n, "batches": v5b.V7_SHIFT_BATCHES, "layouts": out,
           "estimator": "whole-grid secondary-proton edep per primary (common random numbers)",
           "pass": bool(ok)}  # fmt: skip
    return finish(doc, v5b.V7_SHIFT_N, n, n < v5b.V7_SHIFT_N)


# -- E1-B (r13) and E-shape ------------------------------------------------------------------------
def e1b_shards() -> int:
    if E1B_SHARDS is None:
        raise SystemExit("e1-b: the shard count is fixed in Amendment 14a (E1B_SHARDS is None); fail closed")
    if E1B_SHARDS < 2 or E1B_SHARDS * e1b_shard_n(E1B_SHARDS) < E1B_N:
        raise SystemExit("e1-b: need >= 2 shards with a total of >= 1e5 histories")
    return E1B_SHARDS


def step_e1_b_unsized(a: argparse.Namespace) -> int:
    doc = {"step": "e1-b-unsized", "error": "E1B_SHARDS is None: the number of e1-b shards is fixed in "
           "Amendment 14a from the throughput step (>= 2 shards, total >= 1e5 histories)", "pass": False}
    return finish(doc, 0, 0, False)


E1B_GEOMETRY = {"half_mm": 60.0, "depth_mm": 300.0, "bin_mm": 0.5}


def step_e1_b_shard(a: argparse.Namespace) -> int:
    from ionmc.reference import metrics as rm

    ns = e1b_shards()
    k = a.shard
    if not 0 <= k < ns:
        raise SystemExit("e1-b-shard: --shard outside 0..E1B_SHARDS-1")
    nfull = e1b_shard_n(ns)
    nb = E1B_BATCHES
    n = v5b.scaled(nfull, a.scale, 2000, nb)
    g = E1B_GEOMETRY
    geo = v5b.BoxPhantom((-g["half_mm"],) * 2 + (0.0,), (2 * g["half_mm"], 2 * g["half_mm"], g["depth_mm"]), M.WATER)
    grid = ScoringGrid((-g["half_mm"], -g["half_mm"], 0.0), (2 * g["half_mm"], 2 * g["half_mm"], g["bin_mm"]),
                       (1, 1, int(g["depth_mm"] / g["bin_mm"])), name="dose")  # fmt: skip
    seed = seed_c("e1-b", k)
    cfg = c.nuc_config_c(energy=150.0, n=n, seed=seed, geometry=geo, grid=grid, n_batches=nb, timeout=a.timeout)
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    idd_b = res.grid("dose").batch_energy_mev.reshape(nb, -1) / n * nb  # MeV per history, per batch
    depth = (np.arange(idd_b.shape[1]) + 0.5) * g["bin_mm"]
    idd = idd_b.mean(axis=0)
    curve = rm.normalize_to_peak(idd)
    met = {"peak_depth_mm": rm.peak_depth(depth, curve), "r80_mm": rm.r80(depth, curve), "r90_mm": rm.r90(depth, curve)}
    ok = bool(v5.clean(res))
    part = c.write_partial_c(a, f"e1-b-s{k}", {"row": "e1-b", "shard": k, "seed": seed, "n": n, "n_batches": nb,
                             "bin_mm": g["bin_mm"], "idd_mev_per_history_batches": idd_b.tolist(),
                             "metrics": met, "valid": ok, "reduced": n < nfull})  # fmt: skip
    doc = {"step": "e1-b-shard", "shard": k, "seed": seed, "table": v5.table_record(), **c.elastic_table_record(),
           "wall_s": wall, "hist_per_s": n / wall, "metrics": met, "partial": part, "pass": ok}  # fmt: skip
    return finish(doc, nfull, n, n < nfull)


def step_e1_b_compare(a: argparse.Namespace) -> int:
    """Report-only E1-B document: the shards combined (pooled batch curves) with the metrics; the comparison
    with TOPAS and MCsquare is made against the materialized V3-010B runs by the dossier (exploratory, as E1)."""
    ns = e1b_shards()
    parts = [p for p, _ in c.load_partials_c(a, [f"e1-b-s{k}.json" for k in range(ns)])]
    for k, p in enumerate(parts):
        if p["seed"] != seed_c("e1-b", k) or p["shard"] != k or not p["valid"]:
            raise SystemExit(f"e1-b shard {k}: seed/shard mismatch or invalid result")
    guard = seed_guard({"e1-b": [p["seed"] for p in parts]})
    from ionmc.reference import metrics as rm

    batches = np.concatenate([np.asarray(p["idd_mev_per_history_batches"]) for p in parts])
    depth = (np.arange(batches.shape[1]) + 0.5) * parts[0]["bin_mm"]
    idd = batches.mean(axis=0)
    curve = rm.normalize_to_peak(idd)
    met = {"peak_depth_mm": rm.peak_depth(depth, curve), "r80_mm": rm.r80(depth, curve),
           "r90_mm": rm.r90(depth, curve), "plateau_over_peak": float(curve[(depth > 20) & (depth < 60)].mean())}  # fmt: skip
    n = sum(p["n"] for p in parts)
    doc = {"step": "e1-b-compare", "status": "exploratory: report only (plan row E1-B)", "histories": n,
           "table": v5.table_record(), **c.elastic_table_record(), "seeds": guard, "metrics": met,
           "idd_mev_per_history": idd.tolist(), "attestation": v5.attestation_block(a), "pass": True}  # fmt: skip
    return finish(doc, E1B_N, n, any(p["reduced"] for p in parts))


E_SHAPE_PLAN = ("E-shape (plan Amendment 15 (b)6 and the addendum table; Amendment 17 (b)3): exploratory low-energy shape "
                "sensitivity: python, 150 MeV, water, 2e4 histories per variant, sigma_el below E_min,shape off "
                "(the model) vs forced isotropic in mu_CM with BGG sigma_el (test-only builder option); seeds of "
                "the rehearsal family, base 20509500, shard k -> 20509500 + k | Report-only, never evidence")


E_SHAPE_WINDOW_MM = 10.0


def exploratory_table_id() -> str:
    """Build (deterministically, from the pinned sources) the forced-isotropic E-shape table into the cache and
    return its id (never committed; the loader refuses it unless a configuration names it explicitly)."""
    from ionmc.data import cache
    from ionmc.nuclear import elastic_build as eb

    return str(eb.build_elastic_proton(cache.resolve_cache_dir(None), exploratory_forced_isotropic=True).table_id)


def e_shape_compare(idd_a: NDArray[np.float64], idd_b: NDArray[np.float64], depth: NDArray[np.float64],
                    window_mm: float = E_SHAPE_WINDOW_MM) -> dict[str, Any]:
    """Amendment 15 (b)6 comparison of two variants run with the same seeds (batch ``b`` of one is paired with batch
    ``b`` of the other): ``delta IDD(z) = B - A`` in the last ``window_mm`` before the model's R80 with the batch
    standard error, and the R80 / R90 shift (mean of the paired per-batch differences and its standard error)."""
    from ionmc.reference import metrics as rm

    nb = idd_a.shape[0]
    mean_a = idd_a.mean(axis=0)
    r80 = float(rm.r80(depth, rm.normalize_to_peak(mean_a)))
    sel = (depth >= r80 - window_mm) & (depth <= r80)
    d = idd_b - idd_a
    se = d.std(axis=0, ddof=1) / math.sqrt(nb)

    def shift(fn: Any) -> dict[str, float]:
        a = np.array([fn(depth, rm.normalize_to_peak(x)) for x in idd_a])
        b = np.array([fn(depth, rm.normalize_to_peak(x)) for x in idd_b])
        dd = b - a
        return {"model_mm": float(fn(depth, rm.normalize_to_peak(mean_a))),
                "forced_isotropic_mm": float(fn(depth, rm.normalize_to_peak(idd_b.mean(axis=0)))),
                "shift_mm": float(dd.mean()), "se_mm": float(dd.std(ddof=1) / math.sqrt(nb))}  # fmt: skip

    return {"window_mm": [r80 - window_mm, r80], "depth_mm": depth[sel].tolist(),
            "delta_idd": d.mean(axis=0)[sel].tolist(), "delta_idd_se": se[sel].tolist(),
            "delta_idd_relative_to_model": (d.mean(axis=0)[sel] / mean_a[sel]).tolist(),
            "delta_idd_window_mean": float(d.mean(axis=0)[sel].mean()),
            "r80": shift(rm.r80), "r90": shift(rm.r90)}  # fmt: skip


def step_e_shape(a: argparse.Namespace) -> int:
    """E-shape (Amendment 15 (b)6): python, 150 MeV, water, ``E_SHAPE_N`` histories per variant with the same seed
    20509500 + k: the model (sigma_el below E_min,shape off) and the forced-isotropic exploratory table
    (``allow_exploratory_table``, built into the cache at step time). Report-only, never evidence."""
    n = v5b.scaled(E_SHAPE_N, a.scale, 2000, 20)
    seed = E_SHAPE_BASE + a.shard
    guard = seed_guard({"e-shape": [seed]}, allow_rehearsal=True)
    geo, grid, nz, _ = v5.depth_box(M.WATER, 150.0, dz=0.5, half_mm=40.0)
    cfg_a = c.nuc_config_c(energy=150.0, n=n, seed=seed, geometry=geo, grid=grid, n_batches=20, timeout=a.timeout)
    tid_b = exploratory_table_id()
    cfg_b = replace(cfg_a, physics=replace(cfg_a.physics, elastic_table_id=tid_b, allow_exploratory_table=True))
    out: dict[str, Any] = {}
    idd: dict[str, NDArray[np.float64]] = {}
    wall = {}
    for name, cfg in (("model", cfg_a), ("forced_isotropic", cfg_b)):
        t0 = time.perf_counter()
        res = Simulation(cfg).run()
        wall[name] = time.perf_counter() - t0
        idd[name] = res.grid("dose").batch_energy_mev.reshape(20, nz) / (n / 20)  # MeV per history, per batch
        cnt = res.counters.as_dict()
        out[name] = {"valid": bool(v5.clean(res)), "exploratory": bool(res.diagnostics.get("exploratory", False)),
                     "elastic_below_domain": int(cnt.get("elastic_below_domain", 0)),
                     "pp_below_domain": int(cnt.get("pp_below_domain", 0)),
                     "table_id": res.effective_config.summary()["elastic"]["table_id"]}  # fmt: skip
    depth = (np.arange(nz) + 0.5) * 0.5
    cmpd = e_shape_compare(idd["model"], idd["forced_isotropic"], depth)
    ok = bool(out["model"]["valid"] and out["forced_isotropic"]["valid"] and not out["model"]["exploratory"]
              and out["forced_isotropic"]["exploratory"])  # fmt: skip
    doc = {"step": "e-shape", "plan_row": E_SHAPE_PLAN, "verdict": "executed", "report_only": True,
           "tag": "exploratory", "never_evidence": True, "seed": seed, "seeds": guard, "histories_per_variant": n,
           "variants": out, "comparison": cmpd, "exploratory_table_id": tid_b,
           "reading": "forced isotropic in mu_CM with BGG sigma_el below E_min,shape (p + A), the model inside the "
           "declared domain unchanged; window = the last 10 mm before the model's R80; the batches of the two "
           "variants are paired (same seed). With about 45 expected affected events in 2e4 histories the "
           "comparison is an upper-bound check; it changes no row and no tolerance (Amendment 15 (b)6).",
           **c.elastic_table_record(), "table": v5.table_record(), "wall_s": wall, "histories": n,
           "pass": ok}  # fmt: skip
    return finish(doc, E_SHAPE_N, n, n < E_SHAPE_N)


# -- hr5c: V8-stat with V8-nuc, V7 f32/f64 on CUDA --------------------------------------------------
V8NUC_PLAN = (
    "V8-nuc (plan Amendment 14 (d), (i)): python:cpu64 `nuc_local` T12 parity with the (i) sample, as part of the "
    "full hr5c V8 statistical parity | T12 conclusive and passing under its unchanged rule. An inconclusive "
    "verdict is not a pass, and no larger sample is taken afterwards | parity, HR, hr5c r_index 1. (i): 4 shards "
    "x 20 batches x 1500 histories = 1.2e5 (hr5c r_index 1, shards k = 0..3); the cpu64 side is unchanged "
    "(100 x 1e4)"
)


def v8_config_c(name: str, a: argparse.Namespace) -> tuple[SimulationConfig, int, int]:
    backend, prec, k, n_full, nb = V8_SAMPLES_C[name]
    n = v5b.scaled(n_full, a.scale, nb * 100 if backend != "python" else nb * 100, nb)
    geo, grid, _, _ = v5.depth_box(M.WATER, 150.0)
    cfg = c.wcfg_c(backend, prec, energy=150.0, n=n, seed=seed_hrc("v8-stat", k), geometry=geo, grid=grid,
                   n_batches=nb, tallies=v5b.V7_TALLIES, timeout=a.timeout)  # fmt: skip
    return cfg, n, n_full


def step_v8_stat_sample(a: argparse.Namespace) -> int:
    name = a.sample
    if name not in V8_SAMPLES_C:
        raise SystemExit(f"v8-stat-sample: unknown sample {name!r}")
    cfg, n, n_full = v8_config_c(name, a)
    t0 = time.perf_counter()
    est = v5b.batch_estimates(cfg)
    wall = time.perf_counter() - t0
    ok = est_ok(est)
    path = c.write_partial_c(a, f"v8-stat-{name}", {
        "row": "v8-stat", "sample": name, "backend": cfg.run.backend, "precision": cfg.run.precision,
        "seed": cfg.run.seed, "n": n, "n_batches": cfg.run.n_batches, "valid": ok,
        "idd": est["idd"].tolist(), "sec_p": est["sec_p"].tolist(), "nuc_local": est["nuc_local_dose"].tolist(),
        "nuclear_local": est["nuclear_local"].tolist(), "escaped_neutral": est["escaped_neutral"].tolist(),
        "counters_by_name": c.counters_by_name(est), "reduced": n < n_full})  # fmt: skip
    doc = {"step": "v8-stat-sample", "sample": name, "seed": cfg.run.seed, "table": v5.table_record(),
           **c.elastic_table_record(), "wall_s": wall, "hist_per_s": n / wall, "partial": path, "pass": ok}  # fmt: skip
    return finish(doc, n_full, n, n < n_full)


def v8_groups(docs: dict[str, dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    return {"python": [docs[f"python-s{k}"] for k in range(V8NUC_PYTHON_SHARDS)], "cpu64": [docs["cpu64"]],
            "cuda32": [docs["cuda32"]], "cuda64": [docs["cuda64"]]}


def v8_nuc_verdict(t12: dict[str, Any]) -> dict[str, Any]:
    """Row V8-nuc: the ``nuc_local`` profile of the python versus cpu64 T12 verdict is conclusive and passing; an
    inconclusive verdict is not a pass and no larger sample is taken."""
    rec = t12["python_vs_cpu64"]["arrays"]["nuc_local"]
    verdict = rec.get("verdict")
    return {"verdict": verdict, "conclusive": verdict in ("pass", "fail"), "pass": bool(verdict == "pass"),
            "singleton_mass_fraction": rec.get("singleton_mass_fraction"), "n_groups": rec.get("n_groups"),
            "p_value": rec.get("p_value"), "p_value_max_z": rec.get("p_value_max_z"),
            "rule": "T12 unchanged; inconclusive is not a pass; no larger sample afterwards"}  # fmt: skip


def step_v8_stat_compare(a: argparse.Namespace) -> int:
    from ionmc.transport import parity

    names = list(V8_SAMPLES_C)
    docs = {nm: c.load_partials_c(a, [f"v8-stat-{nm}.json"])[0][0] for nm in names}
    for nm, p in docs.items():
        if p["seed"] != seed_hrc("v8-stat", V8_SAMPLES_C[nm][2]) or not p["valid"]:
            raise SystemExit(f"v8-stat {nm}: seed mismatch or invalid result")
    guard = seed_guard({nm: [p["seed"]] for nm, p in docs.items()})
    groups = v8_groups(docs)
    obs = {g: v5b.observables(ps) for g, ps in groups.items()}
    out, ok = {}, True
    for x, y in V8_PAIRS_C:
        v = parity.t12_compare(obs[x], obs[y])
        out[f"{x}_vs_{y}"] = v
        ok &= bool(v["pass"])
    nuc = v8_nuc_verdict(out)
    n = {g: sum(p["n"] for p in ps) for g, ps in groups.items()}
    frozen = {"python": V8NUC_PYTHON_SHARDS * V8NUC_PYTHON_N, "cpu64": V8NUC_WARP_N, "cuda32": V8NUC_WARP_N,
              "cuda64": V8NUC_WARP_N}  # fmt: skip
    doc = {"step": "v8-stat-compare", "plan_row": V8NUC_PLAN, "table": v5.table_record(),
           **c.elastic_table_record(), "t12": out, "histories": n, "frozen_histories": frozen,
           "seeds": guard, "v8_nuc": nuc, "reduced": any(p["reduced"] for p in docs.values()),
           "pass": bool(ok and nuc["pass"]), "attestation": v5.attestation_block(a)}  # fmt: skip
    return base.emit(doc)


def step_v7_f32_cuda(a: argparse.Namespace) -> int:
    n = v5b.scaled(v5b.HR5_F32_N, a.scale, 2000, v5b.V7_BATCHES)
    out = f32_vs_f64_c("warp-cuda", n, v5b.V7_BATCHES, seed_hrc("v7-f32", 0), seed_hrc("v7-f32", 1), a.timeout)
    guard = seed_guard({"v7-f32": [seed_hrc("v7-f32", 0), seed_hrc("v7-f32", 1)]})
    doc = {"step": "v7-f32", "table": v5.table_record(), **c.elastic_table_record(), "result": out,
           "seeds": guard, "pass": out["pass"]}  # fmt: skip
    return finish(doc, v5b.HR5_F32_N, n, n < v5b.HR5_F32_N)


# -- deterministic rows and their r18 / r19 sampled cross-checks (nuclear_v10.py) --------------------
PLAN_FILE = "validation/plans/v3-005-acceptance.md"
PLAN_QUOTES = {  # row -> (plan line, first sentence(s) of the criterion as written)
    "V10": (201, "p-p level above Hale's range: the transported sigma_pp(E) = sigma_NI,theta>=16.26deg(150 MeV)*S(E)/S(150) at "
            "**every** PDG rpp2022 `pp_elastic` point with 0.556 < p_lab <= 0.729 GeV/c (T = 152.4-249.9 MeV) | **Report "
            "only.** The per-point ratios and the flag |model/data - 1| <= 0.08 are recorded; there is no verdict"),
    "V10-A": (202, "p+A nuclear elastic vs measured data (entries and execution in (e)): (a) dsigma/dOmega_CM at theta_CM = 10, 15, "
              "20 deg for every gating angular dataset; (b) the partial integral from 10 deg to the data's first minimum; (c) the "
              "sigma_el level vs sigma_tot(n+C) - sigma_nonel(p+C) at the post-1997 energies | (a) every evaluated point has "
              "|model/data - 1| <= 0.30; (b) every evaluated dataset has |model/data - 1| <= 0.20; (c) every gating point has "
              "|model/proxy - 1| <= 0.12; plus the minimum coverage of (e). Otherwise the row fails, or it is 'not executed' when "
              "the coverage is missing"),
    "V10-A-execution": (218, "(e) items 1-8: gating datasets, level points, frames, model value (nuclear Bessel + unscreened "
                        "relativistic point-Rutherford, added incoherently), measured value (log-linear between the bracketing "
                        "points within 3 deg; theta >= 0.8 theta_min,1 not evaluated), partial integral (trapezoid), minimum coverage "
                        "(item 7, l.238: >= 12 (a) points; >= 1 p+16O dataset at 130-140 MeV and one at 175-205 MeV; >= 1 p+12C; "
                        ">= 3 (b) datasets; >= 4 (c) points), r_index-18 cross-check report-only (item 8)"),
    "X-ENDF": (203, "Model nuclear dsigma/dOmega_CM at 20-40 deg vs LA150 NI (ratio interpolation) for C-12, N-14, Ca-40 at 50 / 100 / "
               "150 MeV | Report only | related evaluation, LV, r_index 19 (sampled cross-check)"),
    "D9-BGG": (208, "The acquired Barashenkov sigma_inel(p+A) and sigma_el(p+A) arrays vs the post-1997 EXFOR set, windows of V1b / "
               "Amendment 7(i) | Report only | related model (TOPAS lineage), LV, deterministic"),
}  # fmt: skip
V10A_XCHECK_N, XENDF_XCHECK_N = 1_000_000, 1_000_000


def _model() -> Any:
    import nuclear_v10 as v10

    return v10.ModelTable(c.ELASTIC_TABLE_ID)


def _plan_doc(row: str) -> dict[str, Any]:
    line, text = PLAN_QUOTES[row]
    return {"plan_file": PLAN_FILE, "plan_line": line, "plan_text": text}


def step_v10_pp_level(a: argparse.Namespace) -> int:
    import nuclear_v10 as v10

    m = _model()
    res = v10.v10_pp_level(m)
    res["v10_pp_c_inventory"] = v10.v10_pp_c_inventory()
    doc = {"step": "v10-pp-level", "row": "V10", **_plan_doc("V10"), **c.elastic_table_record(), "result": res,
           "verdict": "none (report-only; executed, no points in the frozen window)", "evidence_gap": res["evidence_gap"], "pass": True}  # fmt: skip
    return base.emit(doc)


FALLBACK_RULE = ("plan (e)9 (l.241): the neutron shape replaces the Bessel shape only if V10-A (a) or (b) fails for p+O, "
                 "V10-A (c) passes and P7 passes; never because of V5 or V11 alone")


def fallback_trigger(a_or_b_fails_for_po: bool, c_passes: bool, p7_passes: bool | str, *,
                     provisional: bool) -> dict[str, Any]:
    """The (e)9 fallback trigger. ``p7_passes`` is ``True`` / ``False`` or ``"ci-only"`` (not established here);
    the trigger fires only with all three conditions established."""
    fires = bool(a_or_b_fails_for_po and c_passes and p7_passes is True)
    if not c_passes:
        reason = ("(c) fails, so the pre-registered neutron-shape fallback is not triggered; recorded scientific "
                  "failure, new decision required (orchestrator)")
    elif not a_or_b_fails_for_po:
        reason = "(a)/(b) do not fail for p+O: the fallback is not triggered"
    elif p7_passes is not True:
        reason = ("P7 fails: the fallback is not triggered" if p7_passes is False
                  else "P7 not established here (provisional): evaluated by the step v10-a-trigger")
    else:
        reason = "all three conditions of (e)9 hold: the (ii') fallback fires (exploratory results; new decision, next base)"
    return {"rule": FALLBACK_RULE, "a_or_b_fails_for_pO": bool(a_or_b_fails_for_po), "c_passes": bool(c_passes),
            "p7_passes": p7_passes, "fires": fires, "reason": reason, "provisional": provisional}  # fmt: skip


def step_v10_a(a: argparse.Namespace) -> int:
    import nuclear_v10 as v10

    res = v10.v10_a(_model())
    trigger = fallback_trigger(res["a_or_b_fails_for_pO"], bool(res["c_pass"]), "ci-only", provisional=True)
    doc = {"step": "v10-a", "row": "V10-A", **_plan_doc("V10-A"), "execution_rules": _plan_doc("V10-A-execution"),
           **c.elastic_table_record(), "result": res, "status": res["status"],
           "fallback_ii_prime_trigger": trigger, "failures": res["failures"],
           "trigger_final": "evaluated by the step v10-a-trigger (after p7-sampler)",
           "pass": bool(res["pass"])}  # fmt: skip
    doc["status_partial"] = status_partial(a, "v10-a-status", {"step": "v10-a", "pass": bool(res["pass"]),
                                           "status": res["status"], "a_or_b_fails_for_pO": res["a_or_b_fails_for_pO"],
                                           "c_pass": bool(res["c_pass"])})
    return base.emit(doc)


def step_v10_a_trigger(a: argparse.Namespace) -> int:
    """Final evaluation of the (e)9 fallback trigger from the sealed v10-a status partial and the p7 status partial
    (fails closed when either is missing); report only: ``fires`` is information for the orchestrator."""
    doc: dict[str, Any] = {"step": "v10-a-trigger", "row": "V10-A", "report_only": True, **_plan_doc("V10-A-execution")}
    try:
        va = c.load_partials_c(a, ["v10-a-status.json"])[0][0]
        p7 = c.load_partials_c(a, ["p7-status.json"])[0][0]
        if va.get("reduced") or p7.get("reduced"):
            doc["reduced_inputs"] = True
        trig = fallback_trigger(bool(va["a_or_b_fails_for_pO"]), bool(va["c_pass"]), bool(p7["pass"]), provisional=False)
        doc.update(fallback_ii_prime_trigger=trig, v10_a_status=va["status"], v10_a_pass=bool(va["pass"]),
                   inputs={"v10-a-status": va["content_sha256"], "p7-status": p7["content_sha256"]},
                   **c.elastic_table_record(), error=None, **{"pass": True})  # fmt: skip
    except (SystemExit, KeyError, OSError, ValueError) as exc:
        doc.update(error=f"{type(exc).__name__}: {exc}", **{"pass": False})
    doc["attestation"] = v5.attestation_block(a)  # after ALL loads
    return base.emit(doc)


def step_v10_a_xcheck(a: argparse.Namespace) -> int:
    import nuclear_v10 as v10

    n = v4.scaled(V10A_XCHECK_N, a.scale, 20_000)
    seeds = [seed_c("v10-a-xcheck", k) for k in range(40)]
    guard = seed_guard({"v10-a-xcheck": seeds})
    res = v10.v10a_xcheck(_model(), seeds, n)
    res["seeds_used"] = seeds[: len(res["datasets"])]
    doc = {"step": "v10-a-xcheck", "row": "V10-A", "report_only": True, **_plan_doc("V10-A-execution"),
           **c.elastic_table_record(), "seeds": guard, "result": res, "pass": True, "reduced": n < V10A_XCHECK_N,
           "histories": n, "frozen_histories": V10A_XCHECK_N}  # fmt: skip
    return base.emit(doc)


def step_x_endf(a: argparse.Namespace) -> int:
    import nuclear_v10 as v10

    n = v4.scaled(XENDF_XCHECK_N, a.scale, 20_000)
    seeds = [seed_c("x-endf", k) for k in range(9)]
    guard = seed_guard({"x-endf": seeds})
    m = _model()
    res = v10.x_endf(m)
    res["sampled_cross_check"] = v10.x_endf_xcheck(m, seeds, n)
    doc = {"step": "x-endf", "row": "X-ENDF", "report_only": True, **_plan_doc("X-ENDF"), **c.elastic_table_record(),
           "seeds": guard, "result": res, "pass": True, "reduced": n < XENDF_XCHECK_N, "histories": n,
           "frozen_histories": XENDF_XCHECK_N}  # fmt: skip
    return base.emit(doc)


def step_d9_bgg(a: argparse.Namespace) -> int:
    import nuclear_v10 as v10

    doc = {"step": "d9-bgg", "row": "D9-BGG", "report_only": True, **_plan_doc("D9-BGG"),
           "sigma_el_part": "uses the V10-A (c) level proxy " + PLAN_QUOTES["V10-A"][1][PLAN_QUOTES["V10-A"][1].index("(c)"):][:150],
           **c.elastic_table_record(), "result": v10.d9_bgg(_model()), "pass": True}  # fmt: skip
    return base.emit(doc)


STEPS = {
    "lv5c-throughput": step_throughput,
    "r1-d2": step_r1_d2,
    "v8-lv": step_v8_lv,
    "p7-sampler": step_p7,
    "v5-ionmc": step_v5_ionmc,
    "v5-compare": step_v5_compare,
    "v5-attribution": step_v5_attribution,
    "v2b-el-shard": step_v2b_el_shard,
    "v2b-el-combine": step_v2b_el_combine,
    "v7-scan-el": step_v7_scan_el,
    "v7-shift-el": step_v7_shift_el,
    "e1-b-shard": step_e1_b_shard,
    "e1-b-compare": step_e1_b_compare,
    "e1-b-unsized": step_e1_b_unsized,
    "e-shape": step_e_shape,
    "v8-stat-sample": step_v8_stat_sample,
    "v8-stat-compare": step_v8_stat_compare,
    "v7-f32": step_v7_f32_cuda,
    "v10-pp-level": step_v10_pp_level,
    "v10-a": step_v10_a,
    "v10-a-trigger": step_v10_a_trigger,
    "v10-a-xcheck": step_v10_a_xcheck,
    "x-endf": step_x_endf,
    "d9-bgg": step_d9_bgg,
}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--workers", default="1", help="accepted for uniformity; one process")
    ap.add_argument("--energy", type=float, default=150.0, help="v5-ionmc: 150 or 200")
    ap.add_argument("--nuclear", choices=("on", "off"), default="on", help="v5-ionmc")
    ap.add_argument("--shard", type=int, default=0, help="sharded steps: shard index")
    ap.add_argument("--sample", default="cpu64", help="v8-stat-sample: name of V8_SAMPLES_C")
    ap.add_argument("--reference-dir", default=None, help="v5-compare / v5-attribution: dir of REF-* runs")
    ap.add_argument("--out-dir", default="samples", help="partial files")
    ap.add_argument("--dirs", nargs="+", default=["."], help="consumer steps: archive directories")
    ap.add_argument("--partials-manifest", default=None, help="consumer steps: manifest of imports")
    ap.add_argument("--seed-base", "--seed", dest="seed", type=int, default=LV5C_SEED_BASE,
                    help="lv5c base 20481004, hr5c base 20491004, rehearsals 20505000-20509499")  # fmt: skip
    ap.add_argument("--timeout", type=float, default=None)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    if args.step == "v5-ionmc" and args.energy not in v5b.V5_ENERGIES:
        raise SystemExit("v5-ionmc: --energy 150 or 200")
    if args.seed in c.CONSUMED_SEED_BASES_C:
        raise SystemExit(f"seed base {args.seed} is a consumed evidence base (lv5b/hr5/b40d8121)")
    base.SEED_BASE = args.seed
    v4.base.SEED_BASE = args.seed
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
