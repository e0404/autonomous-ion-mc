"""Validation steps of the V3-005 acceptance plan, slice B (suites ``lv5b`` and ``hr5``; decision 0041).

Usage::

    python validation/scripts/transport/steps_v5b.py <step> [--scale F] [--seed-base N] ...

Rows (``validation/plans/v3-005-acceptance.md``, Amendments 6 and 7): V8 (trajectory parity python vs
warp-cpu float64 with nuclear on; statistical parity python / warp-cpu / CUDA), R1 for the nuclear
branch, V5 (the ionmc side: absolute IDD at r = 20 cm; the comparison with TOPAS and MCsquare is
``validation/scripts/reference/compare_idd_v5.py`` on the written partials; the step ``v5-compare``
runs it, enforced, after the four ``v5-ionmc`` steps and consumes the partials only, no seed), V2b and V7, plus the
throughput measurement and its shard table. V6, E1-B and the p-p elastic rows belong to V3-005C
(``R_INDEX`` 12 and 13 stay reserved). Every step prints one JSON document between ``#JSON-BEGIN``
and ``#JSON-END`` (the conventions of ``steps_v5.py``); ``--scale`` < 1 gives a labelled,
non-conformant run; the rehearsal base is ``REHEARSAL_SEED_BASE`` (2046xxxx, never evidence).

Seeds (Amendment 6): ``seed = base + 1000 * r_index + shard`` with the lv5b base 20471004 (Amendment 11; 20441004 consumed) and the
hr5 base 20451004; ``R_INDEX`` (lv5b): throughput 1, v8-lv 2, r1-nuc 3, v5-150-on 4, v5-150-off 5,
v5-200-on 6, v5-200-off 7, v2b 8, v7-scan 9, v7-shift 10, v7-rep 11 (replicate j: shard j, the
1e6 reference: shard ``V7_REPLICATES``); ``R_INDEX_HR`` (hr5): v8-stat 1, v7-f32 2. The python and
warp samples of v8-stat use the shards of ``V8_SAMPLES``.

Definitions made here where the plan leaves them open (all recorded in the documents):

* V8-LV: 256 histories at 150 MeV (water box, all physics on) and a nuclear-dense set (every
  Sigma x 40, 100 MeV, ``V8_DENSE_N`` histories, at least 50 events). The python reference and the
  warp-cpu float64 kernel must agree in the trace (steps, voxels, reasons, RNG blocks), the end
  state, the nuclear trace (event ids, species, genealogy, counts, Z_r, A_r, attempts: identical;
  energies, directions, positions: <= 1e-10) and every counter; tallies and grids <= 1e-10.
* R1-nuc: the A16 ``nuclear=False`` regression of ``steps_v5.step_r1`` and a fixed python
  nuclear-on run (``R1_NUC_N`` histories, seed 20351004, 150 MeV) whose digest must equal
  ``R1_NUC_DIGEST``, measured with the source tree of commit b84fdf38 (the lv5 archive commit; the
  archive stores only verdicts, so the reference values are recomputed from that tree and recorded
  here). After C13 (secondary-proton interactions) the comparison becomes an intended-change record
  (``R1_NUC_INTENDED_CHANGE``).
* V5: warp-cpu float64, nuclear on and off, 1e5 histories in 20 batches, 400 x 400 mm water box
  (r = 20 cm), 0.5 mm depth bins; the partial holds the per-batch IDD in MeV/(g/cm^2)/primary.
* V2b: warp-cpu float64, MCS and straggling off, nuclear on, s_max 0.1 mm against 1.0 mm with
  common random numbers (one seed per shard for both variants); criterion and validity as V2-probe
  (|mean difference| <= 0.002 at every 10 mm depth, standard error at the deepest depth <= 7e-4).
* V7: the estimators are the per-batch secondary-proton dose (edep of generation > 0 protons),
  ``nuclear_local`` and the escaped neutral energy (neutron + gamma); the per-batch values come from
  the exact batch partials (``run_range`` over the batch history ranges, which equals the batch
  decomposition of the full run). The N-scan uses ``V7_BATCHES`` = 100 batches per N so that the
  standard error of the slope of ln(SE) against ln N (about 0.02) resolves the frozen +-0.05.
  Coverage: ``V7_REPLICATES`` replicates of ``V7_REP_N`` histories, ``V7_BINS`` depth bins of the
  secondary-proton dose profile, 1 sigma intervals against the mean of a ``V7_REF_N`` reference.
  Pass rule (Amendment 12, Codex REVIEW-d3f216ca; ``v7_tost_verdict``): replicate-level TOST. With
  ``f_j`` = covered bins / bins used of replicate j, the 90 % two-sided t interval of mean(f_j)
  (``t_{0.95,R-1}`` times sd/sqrt(R)) must lie inside [``V7_COV_LOW``, ``V7_COV_HIGH``] = [0.640, 0.700],
  with at least ``V7_MIN_INTERVALS`` intervals and 10 bins used. The pooled-point gate 0.68 +- 0.03 is
  only reported (``legacy_point_gate_pass``).
  Grid shift and refinement (decision of this module, no frozen tolerance): the whole-grid secondary
  proton dose of a lateral half-voxel shift agrees to ``V7_SHIFT_RTOL`` (edge strips) and of a 2x
  lateral refinement to ``V7_REFINE_RTOL`` (quantisation). f32 vs f64 on the same backend: paired
  z of every estimator within 3.
* hr5 V8-stat: T12-style samples (``parity.t12_compare``) of the IDD, the secondary-proton and
  ``nuclear_local`` profiles and the scalars total deposit, R80, ``nuclear_local`` and escaped neutral
  energy for python (sharded) against warp-cpu float64, CUDA float32 and CUDA float64.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import hashlib
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
from ionmc import materials as M
from ionmc.config import DiagnosticsOptions, SimulationConfig
from ionmc.geometry import BoxPhantom
from ionmc.nuclear.tables import NuclearTable
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Result, Simulation
from ionmc.transport import parity
from ionmc.transport.run import run_range
from ionmc.transport.tally import NUCLEAR_TALLY_NAMES, QUANTUM_MEV

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FORMAT = 1
QUALIFICATION_SEED_BASE = 20471004  # lv5b (Amendment 11; 20441004 consumed by the failed V7 replicate step)
HR5_SEED_BASE = 20451004  # hr5 (Amendment 6)
REHEARSAL_SEED_BASE = 20461004  # V3-005B rehearsals, never evidence
R_INDEX = {"throughput": 1, "v8-lv": 2, "r1-nuc": 3, "v5-150-on": 4, "v5-150-off": 5,
           "v5-200-on": 6, "v5-200-off": 7, "v2b": 8, "v7-scan": 9, "v7-shift": 10,
           "v7-rep": 11}  # fmt: skip
R_INDEX_RESERVED = {"v6": 12, "e1-b": 13}  # V3-005C
R_INDEX_HR = {"v8-stat": 1, "v7-f32": 2}
R_INDEX_HR_RESERVED = {"v6-cuda": 3}  # V3-005C
STEP_LIMIT_S, MARGIN = v5.STEP_LIMIT_S, v5.MARGIN

# -- declared history counts ---------------------------------------------------------------------
THROUGHPUT_PY_N, THROUGHPUT_WARP_N = 200, 10_000
V8_K, V8_DENSE_N, V8_DENSE_FACTOR, V8_DENSE_MIN_EVENTS, V8_TOL = 256, 96, 40.0, 50, 1e-10
R1_NUC_N, R1_NUC_SEED = 400, 20351004
R1_NUC_BASELINE = "b84fdf389e8bae3b3820b73497086d071c3325ad"
R1_NUC_DIGEST = "3edfce9ce5af0e1008d3012440dce09c22e86a308d80fe4b2c926882dc208749"
NUCLEAR_INTENDED_CHANGE: dict[str, Any] = {
    "task": "V3-005B-C13",
    "change": "secondary protons (generation >= 1, species 0) undergo non-elastic interactions on "
              "the python and warp backends (decision 0041 section 3; deuterons still none)",
    "baseline": R1_NUC_BASELINE,
    "baseline_digest": R1_NUC_DIGEST,
    "regression_toggle": "SECONDARY_NUCLEAR = False (ionmc.transport.reference and "
                         "ionmc.transport.kernels_nuclear): primary-only histories, which must "
                         "reproduce baseline_digest bitwise",
    "may_change": ["nuclear event count per history", "secondary-proton dose and fluence",
                   "nuclear tallies (nuclear_local, escaped neutral energy, binding, alpha_local)",
                   "nuclear diagnostics block", "energy_balance nuclear fields"],
    "must_not_change": ["the nuclear=False digest (A16 c862edf7...)",
                        "primary-only histories (toggle off) vs baseline_digest",
                        "per-particle RNG streams and slot layout of the primary",
                        "energy balance closure (<= 1e-12) and zero fail-closed counters"],
    "expected": {  # plan Amendment 7: secondary-proton events per history and their kinetic energy
        "150": {"events_per_history": 0.0037, "energy_share_of_e0": 0.0012},
        "200": {"events_per_history": 0.0096, "energy_share_of_e0": 0.0032},
    },
    "measured_by": "tests/ionmc/test_nuclear_secondary.py::test_secondary_event_rate_against_the_plan",
}
"""The intended-change record of the nuclear-on outputs after C13 (A16-style, decision 0039 section
A16): the R1-nuc step compares the nuclear-on digest with the lv5 archive baseline only through
this record: with ``SECONDARY_NUCLEAR`` off the digest must still equal ``R1_NUC_DIGEST`` (gated),
with it on the digest differs by the secondary interactions and the run must be valid."""
R1_NUC_INTENDED_CHANGE: dict[str, Any] | None = NUCLEAR_INTENDED_CHANGE
V5_N, V5_BATCHES, V5_DZ_MM, V5_HALF_MM = 100_000, 20, 0.5, 200.0
V5_ENERGIES = (150.0, 200.0)
V2B_N, V2B_BATCHES = 1_000_000, 20  # per variant in total
V2B_SHARDS = run_suite.V2B_SHARDS  # defined in run_suite.py (host interpreter has no ionmc)
V2B_VARIANTS = ((0.1, 0.02), (1.0, 0.02))  # (s_max mm, f_E): a, b
V7_SCAN_N, V7_BATCHES, V7_SLOPE, V7_SLOPE_TOL = (10_000, 100_000, 1_000_000), 100, -0.5, 0.05
V7_F32_N, V7_F32_Z = 100_000, 3.0
V7_SHIFT_N, V7_SHIFT_BATCHES, V7_SHIFT_RTOL, V7_REFINE_RTOL = 100_000, 20, 1e-3, 1e-6
# Amendment 12 (Codex REVIEW-d3f216ca): 300 replicates (30 -> 90 in Amendment 11); the reference shard
# index V7_REPLICATES = 300 stays inside the row's 1000-seed block.
V7_REPLICATES, V7_REP_N, V7_REF_N, V7_BINS, V7_REP_BATCHES = 300, 10_000, 1_000_000, 12, 20
# Replicate-level TOST (two one-sided tests, alpha = 0.05 each, i.e. a 90 % two-sided t interval of the
# mean per-replicate coverage) against the equivalence region [V7_COV_LOW, V7_COV_HIGH]: centre 0.670 =
# P(|t_19| <= 1), the nominal coverage of one-sigma intervals from the SEM of 20 batches, margin +-0.03.
V7_COV_LOW, V7_COV_HIGH, V7_TOST_ALPHA = 0.640, 0.700, 0.05
V7_MIN_INTERVALS = 300
# former point-estimate gate (pooled coverage within 0.68 +- 0.03): reported only, no longer decides ``pass``
V7_LEGACY_TARGET, V7_LEGACY_TOL = 0.68, 0.03
HR5_PYTHON_SHARDS, HR5_PYTHON_N, HR5_PYTHON_BATCHES = 2, 12_000, 20
HR5_WARP_N, HR5_WARP_BATCHES, HR5_F32_N = 1_000_000, 100, 1_000_000
V8_SAMPLES = {  # name -> (backend, precision, shard index, histories, batches)
    "python-s0": ("python", "float64", 0, HR5_PYTHON_N, HR5_PYTHON_BATCHES),
    "python-s1": ("python", "float64", 1, HR5_PYTHON_N, HR5_PYTHON_BATCHES),
    "cpu64": ("warp-cpu", "float64", 2, HR5_WARP_N, HR5_WARP_BATCHES),
    "cuda32": ("warp-cuda", "float32", 3, HR5_WARP_N, HR5_WARP_BATCHES),
    "cuda64": ("warp-cuda", "float64", 4, HR5_WARP_N, HR5_WARP_BATCHES),
}
V8_PAIRS = (("python", "cpu64"), ("cpu64", "cuda32"), ("cpu64", "cuda64"), ("cuda32", "cuda64"))
# measured warp-cpu / python rates of this host (hist/s, single process, nuclear on, 150 MeV): the
# declared shard counts must keep every step at or below STEP_LIMIT_S / MARGIN at these rates;
# the throughput step re-measures them and fails if a declared count no longer fits
DECLARED_RATES = {"python": 63.8, "warp-cpu-f64": 10185.0, "v2b-a": 1500.0, "v2b-b": 7000.0}


def seed_of(row: str, shard: int = 0) -> int:
    return base.SEED_BASE + 1000 * R_INDEX[row] + shard


def seed_hr(row: str, shard: int = 0) -> int:
    return base.SEED_BASE + 1000 * R_INDEX_HR[row] + shard


def finish5b(doc: dict[str, Any], frozen: Any, used: Any, reduced: bool) -> int:
    return v5.finish5(doc, frozen, used, reduced)


def wcfg(backend: str, precision: str, **kw: Any) -> SimulationConfig:
    """``steps_v5.nuc_config`` (python float64) re-targeted to ``backend``/``precision``."""
    cfg = v5.nuc_config(**kw)
    return replace(cfg, run=replace(cfg.run, backend=backend, precision=precision))


def scaled(n: int, scale: float, floor: int, multiple: int = 1) -> int:
    return v4.scaled(n, scale, floor, multiple)


def rate_of(res: Result, n: int, wall: float) -> dict[str, Any]:
    return {"histories": n, "wall_s": wall, "hist_per_s": n / wall,
            "transport_s": res.timings.get("transport"), "counters_clean": v5.clean(res)}  # fmt: skip


# -- per-batch estimators ------------------------------------------------------------------------
SEC_P = TallyRequest("sec_p", "dose", "edep", species=("proton",), generation="secondary")
NUC_LOCAL = TallyRequest("nuc_local", "dose", "edep", species=("nuclear_local",))
V7_TALLIES = (SEC_P, NUC_LOCAL)


def channel_slice(eff: Any, name: str) -> tuple[int, int, float]:
    plan = eff.channels
    assert plan is not None
    c = plan.channels[{q.name: q.numerator for q in plan.quantities}[name]]
    return c.offset, c.size, float(c.quantum)


def batch_estimates(cfg: SimulationConfig) -> dict[str, Any]:
    """Per-block, per-primary estimators of ``cfg`` (nuclear on, tallies ``sec_p`` and ``nuc_local``
    on grid 0): ``sec_p`` / ``nuc_local_dose`` and ``idd`` profiles ``[B, voxels]``, the scalars
    ``nuclear_local``, ``escaped_neutral`` ``[B]`` and the summed counters. The ``B = n_batches``
    batches are the CONTIGUOUS history blocks ``[b N/B, (b+1) N/B)`` (the kernels assign the
    internal batch index ``h mod B``, which interleaves, and the exact tally columns exist only per
    partial result, so the blocks are run as partial results of a two-batch configuration whose rows are summed); histories
    are independent, so the blocks are valid batches. The trajectories do not depend on the batch
    structure (counter-based RNG per history), so the sum over blocks equals the full run."""
    cfg1 = replace(cfg, run=replace(cfg.run, n_batches=2))
    eff = Simulation(cfg1).effective
    n, nb = cfg.run.n_histories, cfg.run.n_batches
    if n % nb:
        raise SystemExit("batch_estimates: n_histories must be a multiple of n_batches")
    hpb = n // nb
    off_s, size_s, q_s = channel_slice(eff, "sec_p")
    off_l, size_l, q_l = channel_slice(eff, "nuc_local")
    k0 = len(NUCLEAR_TALLY_NAMES)
    ix = {name: i for i, name in enumerate(NUCLEAR_TALLY_NAMES)}
    out: dict[str, list[Any]] = {k: [] for k in
                                 ("sec_p", "nuc_local_dose", "idd", "nuclear_local", "escaped_neutral")}
    counters = 0
    for b in range(nb):
        part = run_range(eff, b * hpb, (b + 1) * hpb)
        comps = part.tally_components
        col = {name: math.fsum(comps[len(comps) - k0 + i]) for name, i in ix.items()}
        acc = part.channel_acc
        assert acc is not None
        out["sec_p"].append(acc[:, off_s : off_s + size_s].sum(axis=0) * q_s / hpb)
        out["nuc_local_dose"].append(acc[:, off_l : off_l + size_l].sum(axis=0) * q_l / hpb)
        out["idd"].append(part.edep[0].sum(axis=0) * QUANTUM_MEV / hpb)
        out["nuclear_local"].append(col["nuclear_local"] / hpb)
        out["escaped_neutral"].append(
            (col["nuclear_escaped_neutron"] + col["nuclear_escaped_gamma"]) / hpb)  # fmt: skip
        counters += int(sum(part.counter_sums))
    res = {k: np.array(v) for k, v in out.items()}
    res["counters_sum"] = counters
    res["batch_assignment_ok"] = True  # blocks by construction (kept for the document schema)
    return res


def rel_se(x: NDArray[np.float64]) -> tuple[float, float, float]:
    """Mean, standard error of the mean and relative standard error of batch values."""
    mean, sem = float(x.mean()), float(x.std(ddof=1) / math.sqrt(x.size))
    return mean, sem, (sem / mean if mean != 0.0 else math.inf)


def v7_config(backend: str, precision: str, n: int, nb: int, seed: int, *,
              grid: ScoringGrid | None = None, geo: Any = None,
              timeout: float | None = None) -> SimulationConfig:  # fmt: skip
    g0, gr0, _, _ = v5.depth_box(M.WATER, 150.0)
    return wcfg(backend, precision, energy=150.0, n=n, seed=seed, geometry=geo or g0,
                grid=grid or gr0, n_batches=nb, tallies=V7_TALLIES, timeout=timeout)  # fmt: skip


# -- throughput and the shard table -------------------------------------------------------------
def measure(backend: str, precision: str, e: float, n: int, **kw: Any) -> dict[str, Any]:
    geo, grid, _, _ = v5.depth_box(M.WATER, e)
    if backend != "python":  # compile once, outside the timed run
        Simulation(wcfg(backend, precision, energy=e, n=200, seed=1, geometry=geo, grid=grid,
                        n_batches=2)).run()  # fmt: skip
    cfg = wcfg(backend, precision, energy=e, n=n, seed=1, geometry=geo, grid=grid, n_batches=2, **kw)
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    return rate_of(res, n, time.perf_counter() - t0)


def shard_plan(items: dict[str, tuple[int, str, int]], rates: dict[str, float]) -> dict[str, Any]:
    """``name -> (histories, rate key, declared shards)``: planned time per shard at the measured
    rate, the shards needed for ``STEP_LIMIT_S / MARGIN`` and whether the declared count suffices."""
    out = {}
    for name, (n, key, declared) in items.items():
        r = rates[key]
        need = max(1, math.ceil(n / r / (STEP_LIMIT_S / MARGIN)))
        out[name] = {"histories": n, "rate_hist_per_s": r, "planned_s_total": n / r,
                     "planned_s_per_declared_shard": n / r / declared, "declared_shards": declared,
                     "shards_needed_with_margin": need, "histories_per_3300s_step_with_25pct_margin":
                     int(STEP_LIMIT_S / MARGIN * r), "ok": bool(declared >= need)}  # fmt: skip
    return out


def step_throughput(a: argparse.Namespace) -> int:
    npy = scaled(THROUGHPUT_PY_N, a.scale, 20)
    nw = scaled(THROUGHPUT_WARP_N, a.scale, 500)
    out: dict[str, Any] = {}
    for e in v5.THROUGHPUT_ENERGIES:
        out[f"python-f64@{e:g}"] = measure("python", "float64", e, npy, timeout=a.timeout)
    for prec in ("float64", "float32"):
        for e in v5.THROUGHPUT_ENERGIES:
            out[f"warp-cpu-{prec}@{e:g}"] = measure("warp-cpu", prec, e, nw)
    geo, grid, _, _ = v5.depth_box(M.WATER, 150.0)
    for tag, (smax, frac) in zip("ab", V2B_VARIANTS, strict=True):  # V2b configuration, nuclear on
        cfg = wcfg("warp-cpu", "float64", energy=150.0, n=nw, seed=1, geometry=geo, grid=grid,
                   mcs=False, straggling=False, max_step=smax, frac=frac, tallies=v5.PRIMARY)  # fmt: skip
        t0 = time.perf_counter()
        res = Simulation(cfg).run()
        out[f"v2b-{tag}-smax{smax:g}"] = rate_of(res, nw, time.perf_counter() - t0)
    for e in V5_ENERGIES:  # V5 geometry (r = 20 cm, 0.5 mm bins), nuclear on
        geo5, grid5, _, _ = v5_geometry(e)
        cfg = wcfg("warp-cpu", "float64", energy=e, n=nw, seed=1, geometry=geo5, grid=grid5)
        t0 = time.perf_counter()
        res = Simulation(cfg).run()
        out[f"v5-geometry@{e:g}"] = rate_of(res, nw, time.perf_counter() - t0)
    rates = {
        "python": min(out[f"python-f64@{e:g}"]["hist_per_s"] for e in v5.THROUGHPUT_ENERGIES),
        "warp-cpu-f64": min(out[f"warp-cpu-float64@{e:g}"]["hist_per_s"] for e in v5.THROUGHPUT_ENERGIES),
        "v2b-a": out["v2b-a-smax0.1"]["hist_per_s"], "v2b-b": out["v2b-b-smax1"]["hist_per_s"],
        "v5": min(out[f"v5-geometry@{e:g}"]["hist_per_s"] for e in V5_ENERGIES),
    }  # fmt: skip
    v2b_rate = 1.0 / (1.0 / rates["v2b-a"] + 1.0 / rates["v2b-b"])  # both variants per history
    table = shard_plan({
        "v5-ionmc (per step)": (V5_N, "v5", 1),
        "v2b (both variants)": (V2B_N, "v2b", V2B_SHARDS),
        "v7-scan (1e4 + 1e5 + 1e6, f64, plus f32 1e5)": (sum(V7_SCAN_N) + V7_F32_N, "warp-cpu-f64", 1),
        "v7-shift (3 layouts)": (3 * V7_SHIFT_N, "warp-cpu-f64", 1),
        "v7-rep (300 x 1e4 + 1e6 reference)": (V7_REPLICATES * V7_REP_N + V7_REF_N, "warp-cpu-f64", 1),
        "hr5 python sample (2.4e4)": (HR5_PYTHON_SHARDS * HR5_PYTHON_N, "python", HR5_PYTHON_SHARDS),
    }, {**rates, "v2b": v2b_rate})  # fmt: skip
    ok = all(v["counters_clean"] for v in out.values())
    doc = {"step": "lv5b-throughput", "table": v5.table_record(), "throughput": out,
           "rates_hist_per_s": {**rates, "v2b_combined": v2b_rate}, "step_limit_s": STEP_LIMIT_S,
           "margin": MARGIN, "shard_table": table,
           "shard_table_ok": bool(all(v["ok"] for v in table.values())),
           "declared_rates_reference": DECLARED_RATES,
           "pass": bool(ok and (a.scale < 1.0 or all(v["ok"] for v in table.values())))}  # fmt: skip
    return finish5b(doc, THROUGHPUT_WARP_N, nw, nw < THROUGHPUT_WARP_N)


# -- V8-LV ---------------------------------------------------------------------------------------
@contextlib.contextmanager
def scaled_rows(factor: float) -> Iterator[None]:
    """Multiply every Sigma (and majorant) of ``NuclearTable.material_rows`` by ``factor``."""
    orig = NuclearTable.material_rows

    def scaled_fn(self: NuclearTable, material: Any, f_e: float = 0.02) -> Any:
        m = orig(self, material, f_e)
        return replace(m, sigma_mass_cm2_g=m.sigma_mass_cm2_g * factor,
                       sigma_hat_window=m.sigma_hat_window * factor,
                       sigma_hat_end=m.sigma_hat_end * factor,
                       cum_sigma_mass_cm2_g=m.cum_sigma_mass_cm2_g * factor)  # fmt: skip

    NuclearTable.material_rows = scaled_fn  # type: ignore[method-assign,assignment]
    try:
        yield
    finally:
        NuclearTable.material_rows = orig  # type: ignore[method-assign]


def trace_compare(py: Any, wp_: Any) -> dict[str, Any]:
    """Parity of a python and a warp partial (the discrete columns identical, continuous maxima)."""
    out: dict[str, Any] = {}
    dp, dw = py.diagnostics, wp_.diagnostics
    ident = bool(np.array_equal(dp.trace_int, dw.trace_int) and np.array_equal(dp.end_code, dw.end_code))
    out["trace"] = float(np.max(np.abs(dp.trace_float - dw.trace_float), initial=0.0))
    out["end"] = max(float(np.max(np.abs(dp.end_position_mm - dw.end_position_mm), initial=0.0)),
                     float(np.max(np.abs(dp.end_energy_mev - dw.end_energy_mev), initial=0.0)))
    tp, tw = py.meta["nuclear_trace"], wp_.meta["nuclear_trace"]
    for key in ("events", "secondaries"):
        x, y = tp[key], tw[key]
        if x.shape != y.shape:
            return {"discrete_identical": False, "shape_mismatch": key,
                    "shapes": [list(x.shape), list(y.shape)]}  # fmt: skip
        disc = list(range(11)) if key == "events" else [0, 1, 2, 3, 4]
        cont = [11] if key == "events" else list(range(5, 12))
        ident &= bool(np.array_equal(x[:, disc], y[:, disc]))
        out[key] = float(np.max(np.abs(x[:, cont] - y[:, cont]), initial=0.0))
        out[f"n_{key}"] = int(x.shape[0])
    sp = [math.fsum(c) for c in py.tally_components]
    sw = [math.fsum(c) for c in wp_.tally_components]
    ident &= list(py.counter_sums) == list(wp_.counter_sums)
    out["tallies"] = max(abs(x - y) for x, y in zip(sp, sw, strict=True))
    out["edep"] = max(float(np.max(np.abs(x - y))) for x, y in zip(py.edep, wp_.edep, strict=True)) * QUANTUM_MEV
    out["discrete_identical"] = bool(ident)
    cont_keys = ("trace", "end", "events", "secondaries", "tallies", "edep")
    out["max_continuous"] = max(out[k] for k in cont_keys)
    out["counters_python"], out["counters_warp"] = list(py.counter_sums), list(wp_.counter_sums)
    return out


def v8_pair(e: float, n: int, seed: int, a: argparse.Namespace) -> dict[str, Any]:
    from ionmc.transport.reference import run_reference_range
    from ionmc.transport.warp_driver import run_warp_range

    geo, grid, _, _ = v5.depth_box(M.WATER, e)
    effs = {}
    for backend in ("python", "warp-cpu"):
        cfg = wcfg(backend, "float64", energy=e, n=n, seed=seed, geometry=geo, grid=grid, n_batches=2,
                   timeout=a.timeout)  # fmt: skip
        cfg = replace(cfg, diagnostics=DiagnosticsOptions(track_end_positions=True, trace_histories=n))
        effs[backend] = Simulation(cfg).effective
    t0 = time.perf_counter()
    py = run_reference_range(effs["python"], 0, n)
    t_py = time.perf_counter() - t0
    t0 = time.perf_counter()
    wp_ = run_warp_range(effs["warp-cpu"], 0, n, "cpu")
    t_wp = time.perf_counter() - t0
    m = trace_compare(py, wp_)
    m.update(python_s=t_py, warp_cpu_s=t_wp, energy_mev=e, histories=n, seed=seed,
             nuclear_device_sha256=wp_.meta.get("nuclear_device_sha256"))  # fmt: skip
    return m


def step_v8_lv(a: argparse.Namespace) -> int:
    k = scaled(V8_K, a.scale, 8)
    nd = scaled(V8_DENSE_N, a.scale, 16)
    ci = v8_pair(150.0, k, seed_of("v8-lv", 0), a)
    with scaled_rows(V8_DENSE_FACTOR):
        dense = v8_pair(100.0, nd, seed_of("v8-lv", 1), a)
    need = V8_DENSE_MIN_EVENTS if a.scale == 1.0 else max(1, nd // 2)
    dense["events_required"] = need

    def row_ok(m: dict[str, Any]) -> bool:
        return bool(m.get("discrete_identical") and m.get("max_continuous", 1.0) <= V8_TOL)

    ok = row_ok(ci) and row_ok(dense) and dense.get("n_events", 0) >= need
    doc = {"step": "v8-lv", "table": v5.table_record(), "tolerance": V8_TOL, "ci": ci,
           "nuclear_dense": {**dense, "sigma_factor": V8_DENSE_FACTOR}, "pass": bool(ok)}  # fmt: skip
    return finish5b(doc, V8_K + V8_DENSE_N, k + nd, a.scale < 1.0)


# -- R1-nuc --------------------------------------------------------------------------------------
def r1_nuc_digest() -> dict[str, Any]:
    """Digest of the fixed python nuclear-on run (public API of the lv5 commit only): the per-batch
    grid, the channel accumulators, counters, every energy-balance field, the nuclear event counts."""
    geo, grid, _, _ = v5.depth_box(M.WATER, 150.0)
    cfg = v5.nuc_config(energy=150.0, n=R1_NUC_N, seed=R1_NUC_SEED, geometry=geo, grid=grid,
                        n_batches=2, tallies=(*v5.PRIMARY, *v5.EDEP))  # fmt: skip
    res = Simulation(cfg).run()
    arrays = {"batch_energy": res.grid("dose").batch_energy_mev,
              "channels": res.channel_raw.acc if res.channel_raw is not None else np.zeros(0)}  # fmt: skip
    doc = {
        "arrays": {k: hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest()
                   for k, v in arrays.items()},
        "counters": res.counters.as_dict(), "valid": bool(res.valid),
        "balance": dataclasses.asdict(res.energy_balance),
        "nuclear": res.diagnostics.get("nuclear", {}),
    }  # fmt: skip
    text = json.dumps(doc, sort_keys=True, default=base._json)
    return {"sha256": hashlib.sha256(text.encode()).hexdigest(), "valid": bool(res.valid),
            "summary": {"counters": doc["counters"], "events": sum(
                r["events"] for r in doc["nuclear"].values())}}  # fmt: skip


def step_r1_nuc(a: argparse.Namespace) -> int:
    import run_suite

    a16 = v5._captured(v5.step_r1, argparse.Namespace(scale=a.scale, timeout=a.timeout))
    t0 = time.perf_counter()
    dig = r1_nuc_digest()
    wall = time.perf_counter() - t0
    recorded = R1_NUC_DIGEST
    match = dig["sha256"] == recorded
    toggle = None
    if R1_NUC_INTENDED_CHANGE is None:
        nuc_ok = bool(dig["valid"] and match)
    else:  # intended change: the primary-only toggle reproduces the baseline, the run is valid
        import ionmc.transport.reference as reference

        old = reference.SECONDARY_NUCLEAR
        reference.SECONDARY_NUCLEAR = False
        try:
            toggle = r1_nuc_digest()
        finally:
            reference.SECONDARY_NUCLEAR = old
        nuc_ok = bool(dig["valid"] and toggle["valid"] and toggle["sha256"] == recorded)
    doc = {
        "step": "r1-nuc", "a16_regression": {"pass": a16["pass"], "baseline": a16["a16_baseline"],
                                             "digest_sha256": a16["digest_sha256"],
                                             "intended_change_set": a16["intended_change_set"],
                                             "t1": a16["t1"]},
        "nuclear_on_python": {"histories": R1_NUC_N, "seed": R1_NUC_SEED, "digest": dig["sha256"],
                              "recorded": recorded, "baseline_commit": R1_NUC_BASELINE,
                              "matches_recorded": match, "summary": dig["summary"], "wall_s": wall,
                              "intended_change": R1_NUC_INTENDED_CHANGE,
                              "primary_only_toggle_digest": None if toggle is None
                              else toggle["sha256"]},
        "a16_intended_change_none": run_suite.A16_INTENDED_CHANGE is None,
        "pass": bool(a16["pass"] and nuc_ok and run_suite.A16_INTENDED_CHANGE is None),
    }  # fmt: skip
    return finish5b(doc, R1_NUC_N + 256, R1_NUC_N + a16["histories"], bool(a16["reduced"]))


# -- V5 (ionmc side) -----------------------------------------------------------------------------
def v5_geometry(e: float) -> tuple[BoxPhantom, ScoringGrid, int, float]:
    return v5.depth_box(M.WATER, e, dz=V5_DZ_MM, half_mm=V5_HALF_MM)


def step_v5_ionmc(a: argparse.Namespace) -> int:
    from ionmc.reference import metrics as rm

    e, nuclear = a.energy, a.nuclear == "on"
    row = f"v5-{e:g}-{a.nuclear}"
    n = scaled(V5_N, a.scale, 2000, V5_BATCHES)
    geo, grid, nz, r_mm = v5_geometry(e)
    cfg = wcfg("warp-cpu", "float64", energy=e, n=n, seed=seed_of(row), geometry=geo, grid=grid,
               nuclear=nuclear, n_batches=V5_BATCHES)  # fmt: skip
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    g = res.grid("dose")
    rho = M.WATER.density_g_cm3
    bin_g_cm2 = rho * V5_DZ_MM / 10.0  # mass thickness of a bin of the full-width scoring grid
    idd_b = g.batch_energy_mev.reshape(V5_BATCHES, nz) / bin_g_cm2  # MeV/(g/cm^2)/primary
    idd = idd_b.mean(axis=0)
    depth = (np.arange(nz) + 0.5) * V5_DZ_MM
    curve = rm.normalize_to_peak(idd)
    plateau = (depth >= 20.0) & (depth <= 60.0)
    tot_b = idd_b.sum(axis=1) * bin_g_cm2  # MeV per primary deposited in the grid
    met = {"peak_depth_mm": rm.peak_depth(depth, curve), "r80_mm": rm.r80(depth, curve),
           "plateau_mean_20_60_mm": float(idd[plateau].mean()),
           "peak_over_plateau": float(idd.max() / idd[plateau].mean()),
           "total_in_grid_mev_per_primary": float(tot_b.mean()),
           "total_in_grid_sem": float(tot_b.std(ddof=1) / math.sqrt(V5_BATCHES)),
           "in_grid_over_e0": float(tot_b.mean() / e)}  # fmt: skip
    ok = bool(v5.clean(res) and res.energy_balance.relative_residual <= 1e-12)
    part = v5.write_partial(a, row, {
        "row": row, "energy": e, "nuclear": nuclear, "seed": seed_of(row), "n": n,
        "n_batches": V5_BATCHES, "bin_mm": V5_DZ_MM, "lateral_half_mm": V5_HALF_MM,
        "density_g_cm3": rho, "unit": "MeV/(g/cm^2)/primary", "idd_batches": idd_b.tolist(),
        "metrics": met, "valid": ok, "reduced": n < V5_N,
        "nuclear_mev": dict(res.energy_balance.nuclear_mev) if nuclear else {},
    })  # fmt: skip
    doc = {"step": "v5-ionmc", "row": row, "table": v5.table_record() if nuclear else None,
           "seed": seed_of(row), "range_csda_mm": r_mm, "wall_s": wall, "hist_per_s": n / wall,
           "metrics": met, "relative_residual": float(res.energy_balance.relative_residual),
           "counters": res.counters.as_dict(), "partial": part,
           "comparison": "pending: validation/scripts/reference/compare_idd_v5.py (C15)",
           "pass": ok}  # fmt: skip
    return finish5b(doc, V5_N, n, n < V5_N)


def _v5_reference_runs(ref_dir: Path) -> tuple[list[Path], list[Path]]:
    """TOPAS and MCsquare run directories of ``ref_dir/REF-*`` whose case.json carries a V5 block
    (other materialized runs are not V5 evidence); an unreadable directory stops the step."""
    from ionmc.reference.runs import load_run

    topas: list[Path] = []
    mc: list[Path] = []
    for d in sorted(ref_dir.glob("REF-*")):
        run = load_run(d)
        v5blk = run.case.get("v5")
        if isinstance(v5blk, dict) and v5blk.get("row") == "V5":
            {"topas": topas, "mcsquare": mc}.setdefault(run.engine, []).append(d)
    return topas, mc


def _bound_case_names(node: Any) -> set[str]:
    """Every ``bound_case`` string (``<engine>/<case>``) anywhere in a comparator verdict."""
    if isinstance(node, dict):
        found = {node["bound_case"]} if isinstance(node.get("bound_case"), str) else set()
        return found.union(*(_bound_case_names(v) for v in node.values()))
    if isinstance(node, list):
        return set().union(*(_bound_case_names(v) for v in node))
    return set()


def bound_cases_source_identity(verdict: dict[str, Any], cases_dir: Path) -> dict[str, Any]:
    """Record that the cases the V5 runs were bound to are part of the suite source identity
    (review fe367d22): every file of every bound case directory must be in
    ``run_suite.source_file_list('lv5b')`` (hashed into ``environment.txt`` and attested), else
    ``SystemExit``."""
    import run_suite

    listed = set(run_suite.source_file_list("lv5b"))
    paths: list[str] = []
    for name in sorted(_bound_case_names(verdict)):
        d = cases_dir / name
        rels = sorted(f"{run_suite.V5_CASE_ROOT}/{name}/{f.relative_to(d).as_posix()}"
                      for f in d.rglob("*") if f.is_file())  # fmt: skip
        if not rels:
            raise SystemExit(f"bound case {name}: no files under {d}")
        absent = [r for r in rels if r not in listed]
        if absent:
            raise SystemExit(f"bound case {name}: {absent} not in the suite source identity")
        paths += rels
    if not paths:
        raise SystemExit("v5-compare: the verdict records no bound case")
    return {"cases_in_source_identity": True, "bound_case_paths": paths}


def step_v5_compare(a: argparse.Namespace) -> int:
    """Row V5 verdict (C19 F3): the four verified ``v5-ionmc`` partials against the materialized
    TOPAS and MCsquare runs (``--reference-dir``, default ``<repo>/.ionmc-cache/reference-runs``;
    ``/workspace/...`` in the host snapshot) through ``compare_idd_v5.build_verdict``. The step
    consumes partials only (no seed of its own). Every run must be byte-identical to a committed case
    of ``<repo>/validation/reference_cases`` (C20 G1; the bound file hashes are in the verdict). ``pass`` is the gating (TOPAS) verdict; any
    lineage or input error is a failed document with the reason, never a silent skip."""
    import importlib.util
    import tempfile

    path = REPO / "validation" / "scripts" / "reference" / "compare_idd_v5.py"
    spec = importlib.util.spec_from_file_location("compare_idd_v5", path)
    assert spec and spec.loader
    cmp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cmp)
    names = [f"v5-{e}-{t}.json" for e in (150, 200) for t in ("on", "off")]
    parts = v5.load_partials(a, names)  # hash- and binding-verified, attested
    ref_dir = Path(a.reference_dir) if a.reference_dir else REPO / ".ionmc-cache" / "reference-runs"
    cases_dir = REPO / "validation" / "reference_cases"  # frozen cases of this snapshot (C20 G1)
    doc: dict[str, Any] = {"step": "v5-compare", "reference_dir": str(ref_dir),
                           "cases_dir": str(cases_dir),
                           "attestation": v5.attestation_block(a)}  # fmt: skip
    reduced = any(p["reduced"] for p in parts)
    try:
        topas, mc = _v5_reference_runs(ref_dir)
        with tempfile.TemporaryDirectory() as tmp:
            for nm, p in zip(names, parts, strict=True):
                (Path(tmp) / nm).write_text(json.dumps(p, sort_keys=True))
            verdict = cmp.build_verdict(Path(tmp), topas, mc, cases_dir)
        doc.update(bound_cases_source_identity(verdict, cases_dir))  # SystemExit: failed document
        doc.update(verdict=verdict, pass_gating=bool(verdict["pass"]), error=None)
        doc["pass"] = bool(verdict["pass"] and not reduced)
    except (cmp.IddError, SystemExit, ValueError, OSError) as exc:
        doc.update(cases_in_source_identity=False, verdict=None, error=f"{type(exc).__name__}: {exc}", **{"pass": False})
    return finish5b(doc, V5_N, min(p["n"] for p in parts), reduced)


# -- V2b -----------------------------------------------------------------------------------------
def step_v2b_shard(a: argparse.Namespace) -> int:
    k = a.shard
    n = scaled(V2B_N // V2B_SHARDS, a.scale, 2000, V2B_BATCHES)
    geo, grid, nz, _ = v5.depth_box(M.WATER, 150.0)
    ratios, okc, wall = [], True, 0.0
    for smax, frac in V2B_VARIANTS:
        cfg = wcfg("warp-cpu", "float64", energy=150.0, n=n, seed=seed_of("v2b", k), geometry=geo,
                   grid=grid, mcs=False, straggling=False, max_step=smax, frac=frac,
                   tallies=v5.PRIMARY, n_batches=V2B_BATCHES)  # fmt: skip
        t0 = time.perf_counter()
        res = Simulation(cfg).run()
        wall += time.perf_counter() - t0
        okc &= v5.clean(res)
        ratios.append(v5.primary_ratio(res))
    path = v5.write_partial(a, f"v2b-s{k}", {
        "row": "v2b", "shard": k, "seed": seed_of("v2b", k), "n": n, "valid": bool(okc),
        "diff": (ratios[0] - ratios[1]).tolist(), "s_b": ratios[1].tolist(),
        "reduced": n < V2B_N // V2B_SHARDS})  # fmt: skip
    doc = {"step": "v2b-shard", "row": "v2b", "shard": k, "seed": seed_of("v2b", k),
           "table": v5.table_record(), "wall_s": wall, "hist_per_s": 2 * n / wall,
           "partial": path, "pass": bool(okc)}  # fmt: skip
    return finish5b(doc, V2B_N // V2B_SHARDS, n, n < V2B_N // V2B_SHARDS)


def step_v2b_combine(a: argparse.Namespace) -> int:
    parts = v5.load_partials(a, [f"v2b-s{k}.json" for k in range(V2B_SHARDS)])
    for k, p in enumerate(parts):
        if p["seed"] != seed_of("v2b", k) or p["shard"] != k or not p["valid"]:
            raise SystemExit(f"v2b shard {k}: seed/shard mismatch or invalid result")
    diff = np.concatenate([np.array(p["diff"]) for p in parts])
    n = sum(p["n"] for p in parts)
    nz = diff.shape[1]
    bins = v5.depth_bins(150.0, nz)
    mean = diff.mean(axis=0)[bins]
    sem = diff.std(axis=0, ddof=1)[bins] / math.sqrt(diff.shape[0])
    conclusive = bool(sem[-1] <= v5.PROBE_SIGMA_MAX)
    ok = bool(conclusive and np.all(np.abs(mean) <= v5.PROBE_TOL))
    doc = {"step": "v2b-combine", "table": v5.table_record(),
           "variants": {"a": {"s_max_mm": V2B_VARIANTS[0][0], "f_e": V2B_VARIANTS[0][1]},
                        "b": {"s_max_mm": V2B_VARIANTS[1][0], "f_e": V2B_VARIANTS[1][1]}},
           "histories_per_variant": n, "batches": int(diff.shape[0]),
           "depths_mm": [float(b) for b in bins], "diff": mean.tolist(), "sem": sem.tolist(),
           "max_abs_diff": float(np.abs(mean).max()), "sigma_delta_deepest": float(sem[-1]),
           "sigma_max": v5.PROBE_SIGMA_MAX, "conclusive": conclusive, "tolerance": v5.PROBE_TOL,
           "pass": ok, "attestation": v5.attestation_block(a)}  # fmt: skip
    return finish5b(doc, V2B_N, n, any(p["reduced"] for p in parts))


# -- V7 ------------------------------------------------------------------------------------------
ESTIMATORS = ("sec_p_dose", "nuclear_local", "escaped_neutral")


def scalar_estimators(est: dict[str, Any]) -> dict[str, NDArray[np.float64]]:
    return {"sec_p_dose": est["sec_p"].sum(axis=1), "nuclear_local": est["nuclear_local"],
            "escaped_neutral": est["escaped_neutral"]}  # fmt: skip


def f32_vs_f64(backend: str, n: int, nb: int, seed_f64: int, seed_f32: int, timeout: Any) -> dict[str, Any]:
    """Paired z of every V7 estimator between float32 and float64 on ``backend`` (independent
    seeds: the trajectories diverge chaotically, so the samples are independent)."""
    runs, wall = {}, {}
    for prec, sd in (("float64", seed_f64), ("float32", seed_f32)):
        t0 = time.perf_counter()
        runs[prec] = batch_estimates(v7_config(backend, prec, n, nb, sd, timeout=timeout))
        wall[prec] = time.perf_counter() - t0
    out, ok = {}, bool(all(r["counters_sum"] == 0 and r["batch_assignment_ok"] for r in runs.values()))
    for name in ESTIMATORS:
        x, y = (scalar_estimators(runs[p])[name] for p in ("float64", "float32"))
        z = parity.scalar_z(x, y)
        out[name] = {**z, "pass": bool(abs(z["z"]) < V7_F32_Z)}
        ok &= out[name]["pass"]
    return {"backend": backend, "histories": n, "batches": nb, "seeds": [seed_f64, seed_f32],
            "wall_s": wall, "z_limit": V7_F32_Z, "estimators": out, "pass": bool(ok)}  # fmt: skip


def step_v7_scan(a: argparse.Namespace) -> int:
    ns = [scaled(n, a.scale, 1000 * (j + 1), V7_BATCHES) for j, n in enumerate(V7_SCAN_N)]
    rows: dict[str, Any] = {}
    series: dict[str, list[float]] = {k: [] for k in ESTIMATORS}
    ok = True
    for j, n in enumerate(ns):
        t0 = time.perf_counter()
        est = batch_estimates(v7_config("warp-cpu", "float64", n, V7_BATCHES, seed_of("v7-scan", j),
                                        timeout=a.timeout))  # fmt: skip
        row: dict[str, Any] = {"histories": n, "seed": seed_of("v7-scan", j),
                               "wall_s": time.perf_counter() - t0, "counters_sum": est["counters_sum"],
                               "batch_assignment_ok": est["batch_assignment_ok"]}  # fmt: skip
        ok &= est["counters_sum"] == 0 and est["batch_assignment_ok"]
        for name, x in scalar_estimators(est).items():
            mean, sem, r = rel_se(x)
            row[name] = {"mean": mean, "sem": sem, "relative_se": r}
            series[name].append(r)
        rows[f"N={n}"] = row
    slopes = {}
    for name, rs in series.items():
        slope = float(np.polyfit(np.log(ns), np.log(rs), 1)[0])
        s_ok = bool(abs(slope - V7_SLOPE) <= V7_SLOPE_TOL)
        slopes[name] = {"slope": slope, "limit": [V7_SLOPE - V7_SLOPE_TOL, V7_SLOPE + V7_SLOPE_TOL],
                        "pass": s_ok}  # fmt: skip
        ok &= s_ok
    nf = scaled(V7_F32_N, a.scale, 2000, V7_BATCHES)
    f32 = f32_vs_f64("warp-cpu", nf, V7_BATCHES, seed_of("v7-scan", 3), seed_of("v7-scan", 4), a.timeout)
    ok &= f32["pass"]
    doc = {"step": "v7-scan", "table": v5.table_record(), "batches": V7_BATCHES, "runs": rows,
           "slopes": slopes, "f32_vs_f64": f32, "pass": bool(ok)}  # fmt: skip
    return finish5b(doc, sum(V7_SCAN_N) + V7_F32_N, sum(ns) + nf, a.scale < 1.0)


def lateral_layout(kind: str, nz: int) -> tuple[BoxPhantom, ScoringGrid]:
    """Geometry (half width 45 mm) and the scoring grid ``dose`` of a layout: ``base`` 4 mm lateral
    voxels over +-40 mm, ``shifted`` the same shifted by half a voxel, ``refined`` 2 mm voxels."""
    geo = BoxPhantom((-45.0, -45.0, 0.0), (90.0, 90.0, float(nz)), M.WATER)
    sp, o = {"base": (4.0, -40.0), "shifted": (4.0, -38.0), "refined": (2.0, -40.0)}[kind]
    n = int(round(80.0 / sp))
    return geo, ScoringGrid((o, o, 0.0), (sp, sp, 1.0), (n, n, nz), name="dose")


def step_v7_shift(a: argparse.Namespace) -> int:
    n = scaled(V7_SHIFT_N, a.scale, 2000, V7_SHIFT_BATCHES)
    _, _, nz, _ = v5.depth_box(M.WATER, 150.0)
    tot: dict[str, NDArray[np.float64]] = {}
    ok = True
    for kind in ("base", "shifted", "refined"):
        geo, grid = lateral_layout(kind, nz)
        est = batch_estimates(v7_config("warp-cpu", "float64", n, V7_SHIFT_BATCHES,
                                        seed_of("v7-shift"), grid=grid, geo=geo, timeout=a.timeout))  # fmt: skip
        ok &= est["counters_sum"] == 0 and est["batch_assignment_ok"]
        tot[kind] = est["sec_p"].sum(axis=1)
    out = {}
    for kind, tol in (("shifted", V7_SHIFT_RTOL), ("refined", V7_REFINE_RTOL)):
        d = tot[kind] - tot["base"]
        rel = float(abs(d.mean()) / tot["base"].mean())
        out[kind] = {"tolerance": tol, "relative_mean_difference": rel,
                     "max_abs_batch_difference_mev": float(np.abs(d).max()),
                     "base_total": float(tot["base"].mean()), "pass": bool(rel <= tol)}  # fmt: skip
        ok &= out[kind]["pass"]
    doc = {"step": "v7-shift", "table": v5.table_record(), "seed": seed_of("v7-shift"),
           "histories": n, "batches": V7_SHIFT_BATCHES, "layouts": out,
           "estimator": "whole-grid secondary-proton edep per primary (common random numbers)",
           "pass": bool(ok)}  # fmt: skip
    return finish5b(doc, V7_SHIFT_N, n, n < V7_SHIFT_N)


def coarse_depth(nbins: int) -> tuple[BoxPhantom, ScoringGrid]:
    geo, _, nz, r_mm = v5.depth_box(M.WATER, 150.0)
    dzc = 1.1 * r_mm / nbins
    return geo, ScoringGrid((-40.0, -40.0, 0.0), (80.0, 80.0, dzc), (1, 1, nbins), name="dose")


def tally_batches(res: Result, name: str) -> NDArray[np.float64]:
    plan = res.effective_config.channels
    assert plan is not None
    return v4.batch_values(res, {q.name: q.numerator for q in plan.quantities}[name])


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction of the regularized incomplete beta function (modified Lentz)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 500):
        m2 = 2 * m
        for aa in (m * (b - m) * x / ((qam + m2) * (a + m2)), -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))):
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) > tiny else tiny)
            c = 1.0 + aa / c
            c = c if abs(c) > tiny else tiny
            h *= d * c
        if abs(d * c - 1.0) < 1e-15 and m > 1:
            break
    return h


def betainc_reg(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta function I_x(a, b) (continued fraction; no scipy)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lbt) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbt) * _betacf(b, a, 1.0 - x) / b


def student_abs_prob(df: int, h: float) -> float:
    """P(|t_df| <= h) = 1 - I_{df/(df+h^2)}(df/2, 1/2) for the Student t distribution."""
    return 1.0 - betainc_reg(df / 2.0, 0.5, df / (df + h * h))


def student_t_quantile(p: float, df: int) -> float:
    """Quantile of Student t with ``df`` degrees of freedom (``0.5 < p < 1``), by bisection on the
    exact CDF ``student_abs_prob``: P(t <= q) = (1 + P(|t| <= q)) / 2."""
    lo, hi = 0.0, 1.0
    while (1.0 + student_abs_prob(df, hi)) / 2.0 < p:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (1.0 + student_abs_prob(df, mid)) / 2.0 < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def coverage_cluster_stats(covered_j: NDArray[np.int64], bins_used: int, n_batches: int) -> dict[str, Any]:
    """Cluster-aware coverage statistics: the replicates are the independent clusters (the bins of
    one replicate come from the same histories). ``f_j = covered_j / bins_used``; the SE of the pooled
    coverage (their mean) is ``sd(f_j, ddof=1) / sqrt(R)`` with a 90 % t interval on R - 1 df.
    ``nominal_coverage`` is P(|t_{nb-1}| <= 1): the coverage a one-sigma interval built from the sample
    SEM of ``n_batches`` batches has for a normal estimator."""
    f = np.asarray(covered_j, dtype=float) / bins_used
    r = f.size
    se = float(f.std(ddof=1) / math.sqrt(r))
    half = student_t_quantile(0.95, r - 1) * se
    return {"per_replicate_covered": [int(c) for c in covered_j], "per_replicate_fraction": f.tolist(),
            "coverage_cluster_se": se, "t_quantile_90": student_t_quantile(0.95, r - 1),
            "coverage_ci90_cluster": [float(f.mean() - half), float(f.mean() + half)],
            "nominal_coverage": student_abs_prob(n_batches - 1, 1.0)}  # fmt: skip


def cluster_doc(covered_j: NDArray[np.int64], used: NDArray[np.bool_], ref: NDArray[np.float64],
                ref_sem: NDArray[np.float64], rep_means: NDArray[np.float64], nb: int) -> dict[str, Any]:
    """Cluster-level record of the v7-rep coverage (Amendment 11, review db4474cd). The reference
    is one finite sample shared by every replicate, so its variance is a common-mode term:
    ``reference_variance_ratio`` is the mean over used bins of ``SEM_ref^2`` over the variance of the
    replicate means (expected ~ n_rep / n_ref)."""
    out = coverage_cluster_stats(covered_j, int(used.sum()), nb)
    var_rep = rep_means[:, used].var(axis=0, ddof=1)
    ratio = ref_sem[used] ** 2 / var_rep
    out.update({"reference_profile_sem": ref_sem.tolist(),
                "reference_relative_se": np.where(used, ref_sem / np.where(used, ref, 1.0), 0.0).tolist(),
                "reference_variance_ratio_per_bin": ratio.tolist(),
                "reference_variance_ratio": float(ratio.mean())})  # fmt: skip
    return out


def v7_tost_verdict(covered_j: NDArray[np.int64], bins_used: int, intervals: int,
                    low: float, high: float, alpha: float) -> dict[str, Any]:
    """Replicate-level TOST of the v7-rep coverage (Amendment 12, Codex REVIEW-d3f216ca): with
    ``f_j = covered_j / bins_used`` the 1 - 2 alpha two-sided t interval of ``mean(f_j)`` (replicates
    are the independent clusters) must lie entirely inside ``[low, high]``; also requires
    ``intervals >= V7_MIN_INTERVALS`` and ``bins_used >= 10``. Pure function of its arguments."""
    f = np.asarray(covered_j, dtype=float) / bins_used
    r = f.size
    m = float(f.mean())
    se = float(f.std(ddof=1) / math.sqrt(r))
    t = student_t_quantile(1.0 - alpha, r - 1)
    ci = [m - t * se, m + t * se]
    reasons = []
    if ci[0] < low:
        reasons.append(f"CI lower bound {ci[0]:.4f} below {low}")
    if ci[1] > high:
        reasons.append(f"CI upper bound {ci[1]:.4f} above {high}")
    if intervals < V7_MIN_INTERVALS:
        reasons.append(f"intervals {intervals} < {V7_MIN_INTERVALS}")
    if bins_used < 10:
        reasons.append(f"bins_used {bins_used} < 10")
    return {"mean": m, "se": se, "t": t, "ci": ci, "low": low, "high": high, "alpha": alpha,
            "replicates": r, "pass": not reasons, "reasons": reasons}  # fmt: skip


def step_v7_rep(a: argparse.Namespace) -> int:
    geo, grid = coarse_depth(V7_BINS)
    n_rep = scaled(V7_REP_N, a.scale, 2000, V7_REP_BATCHES)
    n_ref = scaled(V7_REF_N, a.scale, 20_000, V7_REP_BATCHES)

    def run(n: int, k: int) -> NDArray[np.float64]:
        cfg = wcfg("warp-cpu", "float64", energy=150.0, n=n, seed=seed_of("v7-rep", k), geometry=geo,
                   grid=grid, n_batches=V7_REP_BATCHES, tallies=(SEC_P,), timeout=a.timeout)  # fmt: skip
        res = Simulation(cfg).run()
        if not v5.clean(res):
            raise SystemExit(f"v7-rep shard {k}: invalid result")
        return tally_batches(res, "sec_p")

    t0 = time.perf_counter()
    ref_b = run(n_ref, V7_REPLICATES)
    ref = ref_b.mean(axis=0)
    ref_sem = ref_b.std(axis=0, ddof=1) / math.sqrt(ref_b.shape[0])
    used = ref > 0.0
    covered = intervals = 0
    per_bin = np.zeros(V7_BINS)
    covered_j = np.zeros(V7_REPLICATES, dtype=np.int64)
    rep_means = np.zeros((V7_REPLICATES, V7_BINS))
    for j in range(V7_REPLICATES):
        b = run(n_rep, j)
        mean, sem = b.mean(axis=0), b.std(axis=0, ddof=1) / math.sqrt(b.shape[0])
        hit = (np.abs(mean - ref) <= sem) & used
        covered += int(hit.sum())
        intervals += int(used.sum())
        per_bin += hit
        covered_j[j] = int(hit.sum())
        rep_means[j] = mean
    cov = covered / intervals if intervals else math.nan
    tost = v7_tost_verdict(covered_j, int(used.sum()), intervals, V7_COV_LOW, V7_COV_HIGH, V7_TOST_ALPHA)
    legacy = bool(intervals >= V7_MIN_INTERVALS and int(used.sum()) >= 10
                  and abs(cov - V7_LEGACY_TARGET) <= V7_LEGACY_TOL)  # fmt: skip
    doc = {"step": "v7-rep", "table": v5.table_record(), "replicates": V7_REPLICATES,
           "histories_per_replicate": n_rep, "reference_histories": n_ref, "bins": V7_BINS,
           "bins_used": int(used.sum()), "intervals": intervals, "covered": covered,
           "coverage": cov, "target": V7_LEGACY_TARGET, "tolerance": V7_LEGACY_TOL,
           "tost": tost, "legacy_point_gate_pass": legacy,
           "per_bin_covered": per_bin.tolist(), "reference_profile": ref.tolist(),
           **cluster_doc(covered_j, used, ref, ref_sem, rep_means, V7_REP_BATCHES),
           "wall_s": time.perf_counter() - t0, "pass": bool(tost["pass"])}  # fmt: skip
    return finish5b(doc, V7_REPLICATES * V7_REP_N + V7_REF_N, V7_REPLICATES * n_rep + n_ref,
                    n_rep < V7_REP_N)


# -- hr5 -----------------------------------------------------------------------------------------
def v8_config(name: str, a: argparse.Namespace) -> tuple[SimulationConfig, int, int]:
    backend, prec, k, n_full, nb = V8_SAMPLES[name]
    n = scaled(n_full, a.scale, nb * 100, nb)
    geo, grid, _, _ = v5.depth_box(M.WATER, 150.0)
    cfg = wcfg(backend, prec, energy=150.0, n=n, seed=seed_hr("v8-stat", k), geometry=geo, grid=grid,
               n_batches=nb, tallies=V7_TALLIES, timeout=a.timeout)  # fmt: skip
    return cfg, n, n_full


def step_v8_stat_sample(a: argparse.Namespace) -> int:
    name = a.sample
    cfg, n, n_full = v8_config(name, a)
    t0 = time.perf_counter()
    est = batch_estimates(cfg)
    wall = time.perf_counter() - t0
    ok = bool(est["counters_sum"] == 0 and est["batch_assignment_ok"])
    path = v5.write_partial(a, f"v8-stat-{name}", {
        "row": "v8-stat", "sample": name, "backend": cfg.run.backend, "precision": cfg.run.precision,
        "seed": cfg.run.seed, "n": n, "n_batches": cfg.run.n_batches, "valid": ok,
        "idd": est["idd"].tolist(), "sec_p": est["sec_p"].tolist(),
        "nuc_local": est["nuc_local_dose"].tolist(), "nuclear_local": est["nuclear_local"].tolist(),
        "escaped_neutral": est["escaped_neutral"].tolist(), "reduced": n < n_full,
    })  # fmt: skip
    doc = {"step": "v8-stat-sample", "sample": name, "seed": cfg.run.seed, "table": v5.table_record(),
           "wall_s": wall, "hist_per_s": n / wall, "partial": path, "pass": ok}  # fmt: skip
    return finish5b(doc, n_full, n, n < n_full)


def observables(parts: list[dict[str, Any]]) -> parity.T12Observables:
    cat = {k: np.concatenate([np.array(p[k]) for p in parts]) for k in
           ("idd", "sec_p", "nuc_local", "nuclear_local", "escaped_neutral")}
    idd = cat["idd"]
    nb = idd.shape[0]
    arrays = {"idd": idd, "sec_p": cat["sec_p"], "nuc_local": cat["nuc_local"]}
    scalars = {parity.TOTAL_DEPOSIT: idd.sum(axis=1),
               "r80_mm": np.array([parity._r80(idd[i], V8_DZ_MM) for i in range(nb)]),
               "nuclear_local_mev": cat["nuclear_local"], "escaped_neutral_mev": cat["escaped_neutral"]}  # fmt: skip
    hpb = parts[0]["n"] // parts[0]["n_batches"]
    return parity.T12Observables(arrays, scalars, nb, hpb, parts[0]["precision"])


V8_DZ_MM = 1.0


def step_v8_stat_compare(a: argparse.Namespace) -> int:
    names = ["python-s0", "python-s1", "cpu64", "cuda32", "cuda64"]
    docs = {nm: v5.load_partials(a, [f"v8-stat-{nm}.json"])[0] for nm in names}
    for nm, p in docs.items():
        if p["seed"] != seed_hr("v8-stat", V8_SAMPLES[nm][2]) or not p["valid"]:
            raise SystemExit(f"v8-stat {nm}: seed mismatch or invalid result")
    groups = {"python": [docs["python-s0"], docs["python-s1"]], "cpu64": [docs["cpu64"]],
              "cuda32": [docs["cuda32"]], "cuda64": [docs["cuda64"]]}
    obs = {g: observables(ps) for g, ps in groups.items()}
    out, ok = {}, True
    for x, y in V8_PAIRS:
        v = parity.t12_compare(obs[x], obs[y])
        out[f"{x}_vs_{y}"] = v
        ok &= bool(v["pass"])
    n = {g: sum(p["n"] for p in ps) for g, ps in groups.items()}
    frozen = {"python": HR5_PYTHON_SHARDS * HR5_PYTHON_N, "cpu64": HR5_WARP_N, "cuda32": HR5_WARP_N,
              "cuda64": HR5_WARP_N}  # fmt: skip
    doc = {"step": "v8-stat-compare", "table": v5.table_record(), "t12": out, "histories": n,
           "frozen_histories": frozen, "reduced": any(p["reduced"] for p in docs.values()),
           "pass": bool(ok), "attestation": v5.attestation_block(a)}  # fmt: skip
    return base.emit(doc)


def step_v7_f32_cuda(a: argparse.Namespace) -> int:
    n = scaled(HR5_F32_N, a.scale, 2000, V7_BATCHES)
    out = f32_vs_f64("warp-cuda", n, V7_BATCHES, seed_hr("v7-f32", 0), seed_hr("v7-f32", 1), a.timeout)
    doc = {"step": "v7-f32", "table": v5.table_record(), "result": out, "pass": out["pass"]}
    return finish5b(doc, HR5_F32_N, n, n < HR5_F32_N)


STEPS = {
    "lv5b-throughput": step_throughput,
    "v8-lv": step_v8_lv,
    "r1-nuc": step_r1_nuc,
    "v5-ionmc": step_v5_ionmc,
    "v5-compare": step_v5_compare,
    "v2b-shard": step_v2b_shard,
    "v2b-combine": step_v2b_combine,
    "v7-scan": step_v7_scan,
    "v7-shift": step_v7_shift,
    "v7-rep": step_v7_rep,
    "v8-stat-sample": step_v8_stat_sample,
    "v8-stat-compare": step_v8_stat_compare,
    "v7-f32": step_v7_f32_cuda,
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--workers", default="1", help="accepted for uniformity; slice B runs one process")
    ap.add_argument("--energy", type=float, default=150.0, help="v5-ionmc: 150 or 200")
    ap.add_argument("--nuclear", choices=("on", "off"), default="on", help="v5-ionmc")
    ap.add_argument("--shard", type=int, default=0, help="v2b-shard: shard index")
    ap.add_argument("--sample", choices=sorted(V8_SAMPLES), default="cpu64", help="v8-stat-sample")
    ap.add_argument("--out-dir", default="samples", help="partial files")
    ap.add_argument("--dirs", nargs="+", default=["."], help="combine steps: archive directories")
    ap.add_argument("--reference-dir", default=None, help="v5-compare: dir of REF-* runs")
    ap.add_argument("--partials-manifest", default=None, help="combine steps: manifest of imports")
    ap.add_argument("--seed-base", "--seed", dest="seed", type=int, default=QUALIFICATION_SEED_BASE,
                    help="lv5b base 20471004, hr5 base 20451004, rehearsals 2046xxxx")  # fmt: skip
    ap.add_argument("--timeout", type=float, default=None)
    args = ap.parse_args(argv)
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    if args.step == "v5-ionmc" and args.energy not in V5_ENERGIES:
        raise SystemExit("v5-ionmc: --energy 150 or 200")
    if args.step == "v2b-shard" and not 0 <= args.shard < V2B_SHARDS:
        raise SystemExit("v2b-shard: --shard outside 0..V2B_SHARDS-1")
    base.SEED_BASE = args.seed
    v4.base.SEED_BASE = args.seed
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
