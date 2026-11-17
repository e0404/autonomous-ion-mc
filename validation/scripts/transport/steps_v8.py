"""V3-008 revision 2: ionmc ion-transport sample steps for the TOPAS-agreement rows (suite ``lv8``).

Usage::

    python validation/scripts/transport/steps_v8.py td-tl-sample --case N [--k K] [--scale F]
                                                                  [--out-dir DIR]

REGISTERED, NOT RUN. The twelve cases of ``validation/plans/v3-008/r2-topas-target.md`` section 3
(1-9 main, 10-12 held-out) are sampled once each by ionmc (warp-cpu, float64, 1e6 histories, 40
batches, all EM physics on, default options: ``ls = "off"``, ``E_cut = 2 MeV/u``) with three scoring
grids on the TOPAS bin edges (section 4): the depth dose ``idd`` (1 x 1 x L/h_z over the full 120 x 120 mm
field), and the projected lateral profiles ``latx`` (120/h_x x 1 x L) and ``laty`` (1 x 120/h_x x L)
in 1 mm depth bins over the full field. The ionmc sample does not depend on any TOPAS output (sample
seal, plan section "Order").

Seeds (plan section C, r2-topas-target section 7): ``seed = 20581004 + 100 r_index + k`` with
``r_index = 29 + case`` (cases 1-12 give r_index 30-41), ``k = 0`` for the single run and
``k = 1, 2, ...`` for chunks of a run split to stay below 1500 s per host call (chunks are
independent runs; their per-batch grids are combined by the evaluator, not here). The held-out set
after a model amendment uses the re-run base 20621004 with the same r_index (``--seed-base``).

The step writes ``ionmc-case<N>-k<K>.npz`` (per-grid batch energy arrays [MeV] and the grid edges; ``np.savez_compressed``, lossless)
and prints one JSON document (``#JSON-BEGIN`` / ``#JSON-END``) with the sha256 of the npz, the
effective-configuration summary, the counters, validity and the energy balance. Every invalid
result or nonzero counter fails the step. No TD/TL statistic is computed here.

V5-SS ionmc thin-slab samples (row V5-SS, branch B; standalone host calls, not an lv8 step)::

    python validation/scripts/transport/steps_v8.py ss-sample --case N [--k K] [--scale F]
                                                    [--seed-base 20631004] [--out-dir DIR]

writes ``ionmc-ss-case<N>-k<K>.npz`` (per-primary exit angles) and one step document; see
:func:`step_ss_sample` and the "V5-SS execution note (2026-11-13)" of the acceptance plan.

Verification rows of plan section A (suites ``lv8`` and ``hr8`` of ``run_suite.py``)::

    steps_v8.py {i6-lv,i8,i13,t1-ion,r1-ion,r2,t13-ion} [--point P] [--case-index I] [--scale F]
                [--seed-base N]
    steps_v8.py shard --row ROW --unit U --shard J [--out-dir DIR] [--scale F] [--seed-base N]
    steps_v8.py combine --row ROW --dirs OUT [IMPORTED ...] [--partials-manifest M] [--scale F]

Implemented: I6 (LV part), I8 + I8b (one step per standard point), I13 (one step per case),
T1-ion LV, R1-ion, R2 LV (lv8) and T13-ion (hr8); and, sharded under amendment 4c (2e5 warp-cpu
histories per configuration; ``shard`` steps write sealed partials, one ``combine`` step per row
applies the criterion): I9, I12, I14, I15, T12-ion LV (lv8) and the T12-ion CUDA parts (hr8;
CPU f64 parts re-run inside hr8, see ``T12_SAMPLES``). Measured in the
sandbox (2026-11-11, C 290 MeV/u, all EM on, 1 mm steps, one process): warp-cpu float64 about
120 histories/s (8.5 ms per history), python about 0.4-0.5 s per history - some 60 times below
the planning throughput of plan section E (7e3/s), so every 1e6-history row needs about 2.4 h of
warp-cpu time and must be sharded into >= 7 host calls (a budget decision for the orchestrator).
Seeds: ``base + 100 * ROW_R_INDEX[row] + k`` (``row_seed``; k documented in ``ROW_R_INDEX``);
every document carries the seed, ``reduced`` and ``pass``; any invalid run or non-zero counter
fails the step. Definitions made here where the plan leaves them open: I8 uses the end depth of
the transported ion (where E reaches E_cut = 2 MeV/u, before the continuous residual) and the
quadrature from E0 down to E_cut; I8b compares sigma_R/R(transport) with the proton at the same
CSDA range (quadrature, lower limit 2 MeV) divided by sqrt(M/m_p), and for carbon with the band
0.3-0.4 % widened by +-15 %; I13 uses 0.1 mm IDD bins of a 1.1 R water box and ``parity.r80_of``
on the mean profile; R1-ion runs python 100 and warp-cpu float32/float64 10^4 histories twice;
T13-ion runs 1e5 histories per chunk size and precision.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

import steps as base
from ionmc.config import DiagnosticsOptions, PhysicsOptions, RunOptions, SimulationConfig
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, PROTON, Projectile
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport.tables import TransportTables

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FORMAT = 1
LV8_BASE = 20581004
RERUN_BASE_1 = 20621004
R_INDEX_OFFSET = 29  # r_index = 29 + case: cases 1..12 -> 30..41
N_HISTORIES = 1_000_000
N_BATCHES = 40
FIELD_MM = 120.0
SOURCE_Z_MM = -0.001  # as the TOPAS bundles (zero-size pencil on the axis from z = -0.001 mm)
MEMORY_BUDGET_BYTES = 2**32  # the 3 grids x 40 batches of the longest cases exceed 2 GiB (192 MB x ...)

IONS: dict[str, Projectile] = {"he": ALPHA, "carbon": CARBON12, "oxygen": OXYGEN16}
# case -> (ion, MeV/u, box length L mm, IDD bin h_z mm, lateral bin h_x mm, set); r2-topas-target section 3
CASES: dict[int, tuple[str, float, float, float, float, str]] = {
    1: ("he", 50.0, 50.0, 0.02, 0.01, "main"),
    2: ("he", 150.0, 250.0, 0.1, 0.02, "main"),
    3: ("he", 250.0, 500.0, 0.1, 0.02, "main"),
    4: ("carbon", 80.0, 50.0, 0.02, 0.01, "main"),
    5: ("carbon", 290.0, 250.0, 0.1, 0.02, "main"),
    6: ("carbon", 430.0, 450.0, 0.1, 0.02, "main"),
    7: ("oxygen", 100.0, 50.0, 0.02, 0.01, "main"),
    8: ("oxygen", 350.0, 250.0, 0.1, 0.02, "main"),
    9: ("oxygen", 430.0, 350.0, 0.1, 0.02, "main"),
    10: ("he", 100.0, 150.0, 0.05, 0.02, "held-out"),
    11: ("carbon", 150.0, 100.0, 0.05, 0.02, "held-out"),
    12: ("oxygen", 200.0, 100.0, 0.05, 0.02, "held-out"),
}
K_MAX = 99


def r_index_of(case: int) -> int:
    """The lv8 row index of a case (30-41)."""
    if case not in CASES:
        raise ValueError(f"unknown case {case}; cases are 1..12")
    return R_INDEX_OFFSET + case


def seed_of(case: int, k: int = 0, seed_base: int = LV8_BASE) -> int:
    """``seed_base + 100 r_index + k`` (``0 <= k <= 99``)."""
    if not 0 <= k <= K_MAX:
        raise ValueError(f"chunk index k must be in 0..{K_MAX}")
    return seed_base + 100 * r_index_of(case) + k


def grids_of(case: int) -> tuple[ScoringGrid, ScoringGrid, ScoringGrid]:
    """The three scoring grids (``idd``, ``latx``, ``laty``) on the TOPAS bin edges."""
    _ion, _eu, length, hz, hx, _set = CASES[case]
    half = FIELD_MM / 2.0
    nz_idd = round(length / hz)
    nx_lat = round(FIELD_MM / hx)
    nz_lat = round(length)
    return (
        ScoringGrid((-half, -half, 0.0), (FIELD_MM, FIELD_MM, hz), (1, 1, nz_idd), name="idd"),
        ScoringGrid((-half, -half, 0.0), (hx, FIELD_MM, 1.0), (nx_lat, 1, nz_lat), name="latx"),
        ScoringGrid((-half, -half, 0.0), (FIELD_MM, hx, 1.0), (1, nx_lat, nz_lat), name="laty"),
    )


def sample_config(
    case: int, k: int = 0, *, n_histories: int = N_HISTORIES, n_batches: int = N_BATCHES,
    seed_base: int = LV8_BASE, backend: str = "warp-cpu",
) -> SimulationConfig:  # fmt: skip
    """The ionmc sample configuration of ``case`` (defaults: the frozen run; the arguments allow a
    labelled reduced run)."""
    ion, e_u, length, _hz, _hx, _set = CASES[case]
    proj = IONS[ion]
    return SimulationConfig(
        source=PencilBeamSource(
            proj, (0.0, 0.0, SOURCE_Z_MM), (0.0, 0.0, 1.0), e_u * proj.a
        ),
        geometry=BoxPhantom(
            (-FIELD_MM / 2.0, -FIELD_MM / 2.0, 0.0), (FIELD_MM, FIELD_MM, length), WATER
        ),
        scoring=grids_of(case),
        physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource()),
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision="float64",
            seed=seed_of(case, k, seed_base),
            n_histories=n_histories,
            n_batches=n_batches,
            memory_budget_bytes=MEMORY_BUDGET_BYTES,
        ),
    )


def step_td_tl_sample(a: argparse.Namespace) -> int:
    """Run one sample (the sample-seal artefact; no TD/TL statistic)."""
    n = max(N_BATCHES, int(round(N_HISTORIES * a.scale / N_BATCHES)) * N_BATCHES)
    cfg = sample_config(a.case, a.k, n_histories=n, seed_base=a.seed_base)
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, Any] = {}
    for g, grid in zip(res.grids, cfg.scoring, strict=True):
        arrays[f"{g.name}_batch_energy_mev"] = np.asarray(g.batch_energy_mev)
        for ax in range(3):
            arrays[f"{g.name}_edges{ax}_mm"] = grid.edges_mm(ax)
    npz = out / f"ionmc-case{a.case}-k{a.k}.npz"
    np.savez_compressed(npz, **arrays)  # lateral grids are mostly zero: ~10-20x smaller on disk
    digest = hashlib.sha256(npz.read_bytes()).hexdigest()
    counters = res.counters.as_dict()
    doc = {
        "format": FORMAT,
        "row": "td-tl-sample",
        "case": a.case,
        "r_index": r_index_of(a.case),
        "k": a.k,
        "seed": cfg.run.seed,
        "histories": n,
        "frozen_histories": N_HISTORIES,
        "reduced": n != N_HISTORIES,
        "backend": cfg.run.backend,
        "valid": bool(res.valid),
        "counters": counters,
        "npz": npz.name,
        "npz_sha256": digest,
        "wall_s": wall,
        "effective_config": res.effective_config.summary(),
        "pass": bool(res.valid and not any(counters.values())),
    }
    return base.emit(doc)


# == V5-SS ionmc thin-slab samples (plan row V5-SS, branch B; "V5-SS execution note 2026-11-13") ==
# Seeds ``SS_SEED_BASE + 100 * SS_R_INDEX[case] + k`` (k = chunk index, 0 for the single 4e5 run).
# SS_SEED_BASE = 20631004 is the declared re-run base ``RERUN_BASE_2`` of plan section C (family
# ``rerun2``, r_index <= 9, never drawn before), used here for the FIRST V5-SS run by the dated
# section-C note of 2026-11-13. These r_index values are NOT entries of ``ROW_R_INDEX``: that table
# indexes the lv8 base 20581004, and the two families are disjoint (20631004-20632003 vs
# 20581004-20587103; ``tests/ionmc/test_v3_008_seed_families.py``).
SS_SEED_BASE = 20631004
SS_CASES: tuple[int, ...] = (22, 25, 26, 27, 28, 29, 30, 31)
SS_R_INDEX: dict[int, int] = {case: i for i, case in enumerate(SS_CASES)}
SS_HISTORIES = 400_000
SS_SLAB_MM = 1.0
SS_WINDOW_K = 2.5  # = ss_budget.WINDOW_K (the evaluator imports it from there)
SS_Z_TOL_MM = 1e-9
# case -> (ion, MeV/u); plan section C dated note 2026-11-01 (branch B case numbers)
SS_CASE_TABLE: dict[int, tuple[str, float]] = {
    22: ("he", 150.0), 25: ("carbon", 150.0), 26: ("oxygen", 150.0), 27: ("carbon", 250.0),
    28: ("he", 250.0), 29: ("carbon", 430.0), 30: ("oxygen", 430.0), 31: ("oxygen", 250.0),
}  # fmt: skip


def ss_seed(case: int, k: int = 0, seed_base: int = SS_SEED_BASE) -> int:
    """``seed_base + 100 SS_R_INDEX[case] + k`` (``0 <= k <= 99``)."""
    if case not in SS_R_INDEX:
        raise ValueError(f"unknown V5-SS case {case}; cases are {SS_CASES}")
    if not 0 <= k <= K_MAX:
        raise ValueError(f"chunk index k must be in 0..{K_MAX}")
    return seed_base + 100 * SS_R_INDEX[case] + k


def ss_config(
    case: int, k: int = 0, *, n_histories: int = SS_HISTORIES, seed_base: int = SS_SEED_BASE,
    backend: str = "warp-cpu",
) -> SimulationConfig:  # fmt: skip
    """The V5-SS ionmc configuration of ``case``: the default all-EM physics of
    :func:`sample_config` (MCS on, straggling on, nuclear off, Bethe stopping, ``ls = "off"``,
    default step controls), warp-cpu float64, 40 batches, a pencil beam from z = -0.001 mm on the
    axis into the 120 x 120 x 1 mm water slab (the TOPAS slab), one dummy 1-bin grid (content
    unused) and the escape records."""
    ion, e_u = SS_CASE_TABLE[case]
    proj = IONS[ion]
    half = FIELD_MM / 2.0
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, SOURCE_Z_MM), (0.0, 0.0, 1.0), e_u * proj.a),
        geometry=BoxPhantom((-half, -half, 0.0), (FIELD_MM, FIELD_MM, SS_SLAB_MM), WATER),
        scoring=(
            ScoringGrid(
                (-half, -half, 0.0), (FIELD_MM, FIELD_MM, SS_SLAB_MM), (1, 1, 1), name="slab"
            ),
        ),
        physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource()),
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision="float64",
            seed=ss_seed(case, k, seed_base),
            n_histories=n_histories,
            n_batches=N_BATCHES,
            memory_budget_bytes=MEMORY_BUDGET_BYTES,
        ),
        diagnostics=DiagnosticsOptions(escape_records=True),
    )


def pooled_rms(tx: np.ndarray, ty: np.ndarray, window: float | None = None) -> tuple[float, int]:
    """Pooled projected rms about zero of x and y, optionally restricted to ``|theta| <= window``:
    sqrt(sum[tx^2 1(|tx|<=W) + ty^2 1(|ty|<=W)] / sum[1(|tx|<=W) + 1(|ty|<=W)]). Returns
    (rms, number of projections used)."""
    t = np.concatenate([np.asarray(tx, np.float64), np.asarray(ty, np.float64)])
    if window is not None:
        t = t[np.abs(t) <= window]
    if t.size == 0:
        raise ValueError("no projection inside the window")
    return float(np.sqrt(np.mean(t * t))), int(t.size)


def step_ss_sample(a: argparse.Namespace) -> int:
    """One V5-SS ionmc thin-slab sample (the SS-seal artefact; no V5-SS statistic is graded).

    Exit angles of every primary at the far face: theta_x = atan2(u_x, u_z), theta_y =
    atan2(u_y, u_z) from the escape records. No secondary track exists under this configuration
    (nuclear off; ionmc transports no delta rays or other secondaries without the nuclear
    channels), so the diagnostics hold exactly one end record per history; the step verifies it
    (escape histories == 0..n-1) instead of assuming it. Valid only if ``res.valid``, every
    counter is 0, every history escaped, every exit z is within 1e-9 mm of 1 mm and u_z > 0
    (exit through the far face); otherwise ``pass`` is false and the exit status 1."""
    n = max(N_BATCHES, int(round(SS_HISTORIES * a.scale / N_BATCHES)) * N_BATCHES)
    cfg = ss_config(a.case, a.k, n_histories=n, seed_base=a.seed_base)
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    d = res.diagnostics
    hist = np.asarray(d["escape_history"], dtype=np.int64)
    pos = np.asarray(d["escape_position_mm"], dtype=np.float64)
    u = np.asarray(d["escape_direction"], dtype=np.float64)
    problems = []
    if hist.size != n or not np.array_equal(hist, np.arange(n)):
        problems.append(f"escape records {hist.size} != histories {n} (or not one per history)")
    if u.size and not np.all(u[:, 2] > 0.0):
        problems.append("a primary left with u_z <= 0")
    if pos.size and not np.all(np.abs(pos[:, 2] - SS_SLAB_MM) <= SS_Z_TOL_MM):
        problems.append("an exit z differs from the far face by more than 1e-9 mm")
    tx, ty = np.arctan2(u[:, 0], u[:, 2]), np.arctan2(u[:, 1], u[:, 2])
    anchor, _ = pooled_rms(tx, ty) if hist.size else (float("nan"), 0)
    window = SS_WINDOW_K * anchor
    theta_w, n_in = pooled_rms(tx, ty, window) if hist.size else (float("nan"), 0)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    npz = out / f"ionmc-ss-case{a.case}-k{a.k}.npz"
    np.savez_compressed(npz, theta_x=tx, theta_y=ty, z_exit_mm=pos[:, 2] if pos.size else pos,
                        escape_history=hist, theta_rms_untruncated=np.float64(anchor))  # fmt: skip
    counters = res.counters.as_dict()
    ok = bool(res.valid and not any(counters.values()) and not problems)
    doc = {
        "format": FORMAT, "row": "v5ss", "case": a.case, "r_index": SS_R_INDEX[a.case],
        "k": a.k, "seed": cfg.run.seed, "histories": n, "frozen_histories": SS_HISTORIES,
        "reduced": n != SS_HISTORIES, "backend": cfg.run.backend,
        "precision": cfg.run.precision, "valid": bool(res.valid), "counters": counters,
        "npz": npz.name, "npz_sha256": hashlib.sha256(npz.read_bytes()).hexdigest(),
        "wall_s": wall, "effective_config": res.effective_config.summary(),
        "n_escaped": int(hist.size), "theta_anchor_rad": anchor, "window_k": SS_WINDOW_K,
        "theta_w_rad": theta_w, "n_in_window": n_in, "problems": problems, "pass": ok,
    }  # fmt: skip
    return base.emit(doc)


# == lv8 / hr8 verification rows (plan section A) =================================================
# Seeds: ``seed_base + 100 * ROW_R_INDEX[row] + k`` (plan section C; r_index <= 60, k < 100). The
# td-tl-sample rows hold r_index 30-41; the verification rows below use 42-60 (unused: 52-60
# under lv8; hr8 has its own base 20591004 and uses 54-55). ``k`` is documented per step.
HR8_BASE = 20591004
ROW_R_INDEX: dict[str, int] = {
    "i8": 42,  # k = standard-point index 0 He 150, 1 C 290, 2 O 350
    "i13": 43,  # k = case index 0..8 of I13_CASES
    "t1-ion": 48,  # k = 0 C 290, 1 O 350 (python and warp-cpu share the seed: trace parity)
    "r1-ion": 50,  # k = spec index; both repeats of a spec use the same seed
    "i6-lv": 51,  # k = case index 0..8 (deterministic run; the seed is recorded only)
    "t13-ion": 55,  # hr8; k = 0 float32, 1 float64 (both chunk sizes share the seed)
    # sharded rows (amendment 4c): k = position of the shard in ``shard_plan(row)`` (units in
    # order, shards within a unit; every shard its own seed); ``combine`` runs no transport
    "i9": 44,  # units he150, c290, o350 x 2 shards: k 0-5
    "i12": 45,  # theta-s01..s5 (10, 4, 2, 2 shards: k 0-17), sigma-s01..s5 (k 18-35)
    "i14": 46,  # s01-f0005 16, s01-f002 10, s05-f0005 8, s05-f002 4, s1-f0005 8, s1-f002 2: k 0-47
    "i15": 47,  # theta v1-sh0/sh1, v2-sh0/sh1 (2 each: k 0-7); sigma v1 (4 each), v2 (2): k 8-19
    "t12-ion-lv": 49,  # python 2 shards (k 0-1), cpu64 2 shards (k 2-3)
    "t12-ion-hr": 54,  # hr8: c290-cpu64 k 0-1, c290-cuda32 2, c290-cuda64 3, he150-cpu64 4-5,
    #                    he150-cuda32 6, o350-cpu64 7-8, o350-cuda32 9
}
RESERVED_R_INDEX: dict[str, int] = {}
"""r_index declared for rows not implemented yet (none left: I9 44, I12 45, I14 46, I15 47 and
T12-ion 49 / 54 are in ``ROW_R_INDEX`` since the amendment-4c commit)."""
STANDARD: tuple[tuple[str, float], ...] = (("he", 150.0), ("carbon", 290.0), ("oxygen", 350.0))
I13_CASES: tuple[tuple[str, float], ...] = (
    ("he", 50.0), ("he", 150.0), ("he", 250.0), ("carbon", 80.0), ("carbon", 290.0),
    ("carbon", 430.0), ("oxygen", 100.0), ("oxygen", 350.0), ("oxygen", 430.0),
)  # fmt: skip
I8_N, I13_N, I8_TOL, I8B_FINDING, I13_TOL = 100_000, 100_000, 0.05, 0.15, 0.002
I8B_LITERATURE_C = (0.003, 0.004)  # sigma_R/R(C) ~ 0.3-0.4 % (Schardt et al. RMP 2010; TO VERIFY)
T1_K = 256
R1_SPECS = ("python:float64:100", "warp-cpu:float64:10000", "warp-cpu:float32:10000")
T13_N = 100_000


def row_seed(row: str, k: int = 0, seed_base: int = LV8_BASE) -> int:
    """``seed_base + 100 * ROW_R_INDEX[row] + k`` (``0 <= k <= 99``)."""
    if not 0 <= k <= K_MAX:
        raise ValueError(f"chunk index k must be in 0..{K_MAX}")
    return seed_base + 100 * ROW_R_INDEX[row] + k


def scaled(n: int, scale: float, batches: int = 20) -> int:
    return n if scale == 1.0 else max(batches, int(n * scale) // batches * batches)


def ion_tables(proj: Projectile) -> TransportTables:
    return TransportTables.from_stopping_tables(
        [BetheStoppingSource().table(WATER, proj)], projectile=proj
    )


def csda_mm(proj: Projectile, e_u: float) -> float:
    """R: CSDA range of the project table at the source energy (plan notation) [mm of water]."""
    return ion_tables(proj).range_g_cm2(0, e_u * proj.a) * 10.0 / WATER.density_g_cm3


def ion_config(
    proj: Projectile, e_u: float, *, geometry: Any, scoring: tuple[ScoringGrid, ...],
    backend: str = "warp-cpu", precision: str = "float64", seed: int, n: int, n_batches: int = 20,
    mcs: bool = True, straggling: bool = True, max_step: float = 1.0,
    diag: DiagnosticsOptions | None = None, chunk: int | None = None, frac: float = 0.02,
) -> SimulationConfig:  # fmt: skip
    """Ion pencil beam on the axis at z = 0 (the source API of ``sample_config``), Bethe stopping
    with the default options, every other physics option at its default (in particular the
    end-of-range model ``"default"`` = ``continuous_residual_v1`` for ions)."""
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), e_u * proj.a),
        geometry=geometry,
        scoring=scoring,
        physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource(), straggling=straggling,
                               multiple_scattering=mcs, max_step_mm=max_step,
                               max_energy_loss_fraction=frac),  # fmt: skip
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision=precision,  # type: ignore[arg-type]
            seed=seed, n_histories=n, n_batches=n_batches,
            memory_budget_bytes=MEMORY_BUDGET_BYTES,
            **({"chunk_histories": chunk} if chunk else {}),
        ),  # fmt: skip
        diagnostics=diag or DiagnosticsOptions(),
    )


def water_box(depth: float, dz: float = 1.0, half: float = 30.0) -> tuple[BoxPhantom, ScoringGrid]:
    nz = int(math.ceil(depth / dz))
    geo = BoxPhantom((-half, -half, 0.0), (2 * half, 2 * half, nz * dz), WATER)
    grid = ScoringGrid((-half, -half, 0.0), (2 * half, 2 * half, dz), (1, 1, nz), name="idd")
    return geo, grid


def finish8(doc: dict[str, Any], frozen: Any, used: Any, reduced: bool) -> int:
    return base.emit({"format": FORMAT, **doc, "reduced": bool(reduced), "frozen_histories": frozen,
                      "histories": used})  # fmt: skip


def clean(res: Any) -> bool:
    return bool(res.valid and not res.counters.any_nonzero)


# -- I8 / I8b ---------------------------------------------------------------------------------
def sigma_range_quadrature(proj: Projectile, e0: float, e_lo: float) -> tuple[float, float]:
    """``(sigma_R, R(e0) - R(e_lo))`` [mm of water] from the numpy quadrature of I4b,
    sigma_R^2 = int (d sigma^2 / dx) / S^3 dE, with the implemented Bohr variance (python twin of
    ``em.bohr_variance``, bare charge z as in transport) and the Bethe stopping WITH the default
    corrections (I8: "I4b quadrature with corrections on"); total-energy axis, Simpson in ln E."""
    from ionmc._wpfunc import python_twin
    from ionmc.physics.em import make_em
    from ionmc.physics.stopping import bethe_mass_stopping

    o = BetheStoppingSource().options
    em = python_twin(make_em)
    n = 4001
    u = np.linspace(math.log(e_lo), math.log(e0), n)
    e = np.exp(u)
    s = bethe_mass_stopping(
        e / proj.a, proj, WATER, shell=o.shell, barkas=o.barkas, bloch=o.bloch, mott=o.mott,
        effective_charge=o.effective_charge,
    )  # fmt: skip
    var_t = np.array([em.bohr_variance(float(x), proj.mass_mev, float(proj.z), WATER.z_over_a,
                                       1.0, 10.0) for x in e])  # fmt: skip
    w = np.full(n, 2.0)
    w[1::2] = 4.0
    w[0] = w[-1] = 1.0
    h = u[1] - u[0]
    rng = float(h / 3.0 * np.sum(w * e / s))
    sig2 = float(h / 3.0 * np.sum(w * e * var_t / s**3))
    return math.sqrt(sig2) * 10.0 / WATER.density_g_cm3, rng * 10.0 / WATER.density_g_cm3


def proton_energy_at_range(r_mm: float) -> float:
    """Proton kinetic energy [MeV] whose CSDA range (project table) is ``r_mm`` (bisection)."""
    lo, hi = 2.0, 499.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if base.r_csda_mm(mid) < r_mm else (lo, mid)
    return 0.5 * (lo + hi)


def step_i8(a: argparse.Namespace) -> int:
    """I8 (gating) and I8b (report-only) at one standard point (``--point`` 0 He 150, 1 C 290,
    2 O 350): sigma of the end depth (the depth where E reaches E_cut, MCS off, straggling on)
    of ``I8_N`` warp-cpu float64 histories against the quadrature from E0 down to E_cut."""
    ion, e_u = STANDARD[a.point]
    proj = IONS[ion]
    n = scaled(I8_N, a.scale)
    r = csda_mm(proj, e_u)
    geo, grid = water_box(1.1 * r)
    seed = row_seed("i8", a.point, a.seed_base)
    cfg = ion_config(proj, e_u, geometry=geo, scoring=(grid,), seed=seed,
                     n=n, mcs=False, diag=DiagnosticsOptions(track_end_positions=True))  # fmt: skip
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    z = np.asarray(res.diagnostics["end_position_mm"])[:, 2]
    codes = np.asarray(res.diagnostics["end_code"])
    sig = float(np.std(z, ddof=1))
    e_cut = res.effective_config.e_cut_mev
    sq, rq = sigma_range_quadrature(proj, e_u * proj.a, e_cut)
    dev = sig / sq - 1.0
    ok = bool(clean(res) and np.all(codes == 0) and abs(dev) <= I8_TOL)
    # I8b (report): proton at the same range, scaled by 1/sqrt(M/m_p); literature band for C
    e_p = proton_energy_at_range(r)
    sp, rp = sigma_range_quadrature(PROTON, e_p, 2.0)
    scaled_p = sp / rp / math.sqrt(proj.mass_mev / PROTON.mass_mev)
    rel = sig / r
    i8b: dict[str, Any] = {
        "sigma_R_over_R_transport": rel, "sigma_R_over_R_quadrature": sq / r,
        "proton_energy_same_range_mev": e_p, "proton_sigma_R_over_R_scaled": scaled_p,
        "ratio_vs_proton_scaled": rel / scaled_p,
        "finding_vs_proton": abs(rel / scaled_p - 1.0) > I8B_FINDING,
    }  # fmt: skip
    if ion == "carbon":
        lo, hi = I8B_LITERATURE_C
        i8b["literature_band_C"] = [lo, hi]
        i8b["finding_vs_literature"] = not (lo * (1 - I8B_FINDING) <= rel <= hi * (1 + I8B_FINDING))
    doc = {
        "step": "i8", "ion": ion, "e_mev_per_u": e_u, "seed": cfg.run.seed, "csda_range_mm": r,
        "e_cut_mev": e_cut, "mean_end_depth_mm": float(z.mean()), "sigma_end_depth_mm": sig,
        "sigma_quadrature_mm": sq, "quadrature_range_mm": rq, "relative_deviation": dev,
        "tolerance": I8_TOL, "end_codes_all_cutoff": bool(np.all(codes == 0)),
        "counters": res.counters.as_dict(), "valid": bool(res.valid), "wall_s": wall,
        "i8b_report_only": i8b, "pass": ok,
    }  # fmt: skip
    return finish8(doc, I8_N, n, n < I8_N)


# -- I13 --------------------------------------------------------------------------------------
def step_i13(a: argparse.Namespace) -> int:
    """I13 at one case (``--case-index`` into ``I13_CASES``): R80 of the IDD (0.1 mm bins, MCS
    off, straggling on, ``I13_N`` warp-cpu float64 histories) against R: |R80/R - 1| <= 0.2 %."""
    ion, e_u = I13_CASES[a.case_index]
    proj = IONS[ion]
    n = scaled(I13_N, a.scale)
    r = csda_mm(proj, e_u)
    geo, grid = water_box(1.1 * r, dz=0.1)
    cfg = ion_config(proj, e_u, geometry=geo, scoring=(grid,), n=n, mcs=False,
                     seed=row_seed("i13", a.case_index, a.seed_base))  # fmt: skip
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    from ionmc.transport.parity import r80_of

    idd = np.asarray(res.grid("idd").batch_energy_mev).reshape(cfg.run.n_batches, -1)
    r80 = r80_of(idd.mean(axis=0), 0.1)
    r80_b = [r80_of(b, 0.1) for b in idd]
    dev = r80 / r - 1.0
    doc = {
        "step": "i13", "ion": ion, "e_mev_per_u": e_u, "seed": cfg.run.seed, "csda_range_mm": r,
        "r80_mm": r80, "r80_batch_sem_mm": float(np.std(r80_b, ddof=1) / math.sqrt(len(r80_b))),
        "relative_deviation": dev, "tolerance": I13_TOL, "counters": res.counters.as_dict(),
        "valid": bool(res.valid), "wall_s": wall,
        "pass": bool(clean(res) and math.isfinite(r80) and abs(dev) <= I13_TOL),
    }  # fmt: skip
    return finish8(doc, I13_N, n, n < I13_N)


# -- I6 (LV part) -----------------------------------------------------------------------------
def step_i6(a: argparse.Namespace) -> int:
    """I6 LV: deterministic CSDA end depth (straggling and MCS off) at the standard points and
    ends, python and warp-cpu float64, 4 histories each:
    R(E0) - R(E_cut) - 1e-4 R <= z_end <= R(E0) + 1e-4 R (the CI test's criterion)."""
    out: dict[str, Any] = {}
    ok = True
    for i, (ion, e_u) in enumerate(I13_CASES):
        proj = IONS[ion]
        r = csda_mm(proj, e_u)
        geo, grid = water_box(1.1 * r, dz=2.0)
        for backend in ("python", "warp-cpu"):
            cfg = ion_config(proj, e_u, geometry=geo, scoring=(grid,), backend=backend, n=4,
                             n_batches=2, mcs=False, straggling=False, max_step=2.0,
                             seed=row_seed("i6-lv", i, a.seed_base),
                             diag=DiagnosticsOptions(track_end_positions=True))  # fmt: skip
            res = Simulation(cfg).run()
            t = res.effective_config.tables
            r0 = t.range_g_cm2(0, e_u * proj.a) * 10.0
            rc = t.range_g_cm2(0, res.effective_config.e_cut_mev) * 10.0
            z = np.asarray(res.diagnostics["end_position_mm"])[:, 2]
            lo, hi = r0 - rc - 1e-4 * r0, r0 + 1e-4 * r0
            p = bool(clean(res) and np.all((z >= lo) & (z <= hi))
                     and np.all(np.asarray(res.diagnostics["end_code"]) == 0))  # fmt: skip
            out[f"{ion}-{e_u:g}-{backend}"] = {
                "z_end_mm": z.tolist(), "lower": lo, "upper": hi,
                "residual_range_at_e_cut_mm": rc, "pass": p,
            }  # fmt: skip
            ok &= p
    return finish8({"step": "i6-lv", "cases": out, "pass": ok}, None, None, False)


# -- T1-ion LV ----------------------------------------------------------------------------------
def step_t1_ion(a: argparse.Namespace) -> int:
    """T1-ion LV: K = 256 histories each of C 290 and O 350 MeV/u, all EM on, python against
    warp-cpu float64 on the V3-003 T1 geometry (12 x 12 x 40 voxels of 5 mm, 2 mm steps, 2 mm
    scoring grid); discrete trace columns identical, continuous within rtol = atol = 1e-10
    (``parity.compare_traces``)."""
    from ionmc.geometry import VoxelGeometry
    from ionmc.transport.parity import compare_traces

    k = max(4, int(T1_K * a.scale))
    shape = (12, 12, 40)
    geo = VoxelGeometry((-30.0, -30.0, 0.0), (5.0, 5.0, 5.0), shape, (WATER,),
                        np.zeros(shape, dtype=np.int32))  # fmt: skip
    grid = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 100), name="dose")
    diag = DiagnosticsOptions(track_end_positions=True, trace_histories=k)
    out: dict[str, Any] = {}
    ok = True
    for i, (ion, e_u) in enumerate((("carbon", 290.0), ("oxygen", 350.0))):
        runs = {}
        for backend in ("python", "warp-cpu"):
            cfg = ion_config(IONS[ion], e_u, geometry=geo, scoring=(grid,), backend=backend,
                             seed=row_seed("t1-ion", i, a.seed_base), n=k, n_batches=2,
                             max_step=2.0, diag=diag)  # fmt: skip
            t0 = time.perf_counter()
            runs[backend] = Simulation(cfg).run()
            out.setdefault(ion, {})[f"{backend}_seconds"] = time.perf_counter() - t0
        v = compare_traces(runs["python"].diagnostics, runs["warp-cpu"].diagnostics)
        v["pass"] = bool(v["pass"] and all(clean(r) for r in runs.values()))
        out[ion].update(t1=v, seed=row_seed("t1-ion", i, a.seed_base))
        ok &= v["pass"]
    return finish8({"step": "t1-ion", "k": k, "runs": out, "pass": ok}, T1_K, k, k < T1_K)


# -- R1-ion ---------------------------------------------------------------------------------------
def step_r1_ion(a: argparse.Namespace) -> int:
    """R1-ion: C 290 MeV/u, all EM on, each spec of ``R1_SPECS`` run twice with the same
    configuration and seed; tallies, counters, deposit grids and end states bit-identical
    (``parity.compare_runs_bitwise``; the V3-003 repeatability step with the ion source)."""
    from ionmc.transport.parity import compare_runs_bitwise

    proj, e_u = CARBON12, 290.0
    geo, _ = water_box(1.1 * csda_mm(proj, e_u))
    nz = int(1.1 * csda_mm(proj, e_u) / 2)
    grid = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, nz), name="dose")
    out: dict[str, Any] = {}
    ok = True
    used, frozen = [], []
    for i, spec in enumerate(R1_SPECS):
        backend, prec, n_s = spec.split(":")
        n = scaled(int(n_s), a.scale)
        used.append(n)
        frozen.append(int(n_s))
        cfg = ion_config(proj, e_u, geometry=geo, scoring=(grid,), backend=backend, precision=prec,
                         seed=row_seed("r1-ion", i, a.seed_base), n=n, max_step=2.0,
                         diag=DiagnosticsOptions(track_end_positions=True))  # fmt: skip
        runs = [Simulation(cfg).run() for _ in range(2)]
        v = compare_runs_bitwise(runs[0], runs[1])
        v.update(valid=all(clean(r) for r in runs), n_histories=n, seed=cfg.run.seed)
        out[spec] = v
        ok &= bool(v["identical"] and v["valid"])
    return finish8({"step": "r1-ion", "r1": out, "pass": ok}, frozen, used, used != frozen)


# -- R2 LV ---------------------------------------------------------------------------------------
def _captured(fn: Any, ns: argparse.Namespace) -> tuple[int, dict[str, Any]]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = fn(ns)
    text = buf.getvalue()
    return rc, json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])


def step_r2(a: argparse.Namespace) -> int:
    """R2 LV (proton regression; any difference = STOP): (a) the A16 nuclear-off digests of the
    tree under test against ``steps_v4.A16_BASELINE`` (``steps_v4.step_a16 --mode regression``,
    the gate of the lv5 ``r1`` step), ``run_suite.A16_INTENDED_CHANGE`` must stay None; (b) the T1
    proton trace fixture ``tests/ionmc/data/trace_baseline_20mev.npz``: its frozen guard
    ``test_default_trace_is_unchanged_from_the_stored_baseline`` (python and warp-cpu float64) is
    executed, and the exact (bitwise) equality of every column of the python trace is recorded."""
    import importlib.util

    import run_suite
    import steps_v4 as v4

    ns = argparse.Namespace(mode="regression", intended_change_record=None)
    rc, a16 = _captured(v4.step_a16, ns)
    spec = importlib.util.spec_from_file_location(
        "v8_trace_guard", REPO / "tests" / "ionmc" / "test_transport_diagnostic.py")  # fmt: skip
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    guard: dict[str, Any] = {"test": "test_default_trace_is_unchanged_from_the_stored_baseline"}
    try:
        mod.test_default_trace_is_unchanged_from_the_stored_baseline()
        guard["pass"] = True
    except AssertionError as exc:
        guard.update(error=str(exc)[:2000], **{"pass": False})
    fixture = np.load(mod.BASELINE)
    tr = Simulation(mod._config("python", legacy=True)).run().diagnostics["trace"]
    from ionmc.transport.tally import TRACE_COLUMNS

    bitwise = {c: bool(np.array_equal(tr[c], fixture[c])) for c in TRACE_COLUMNS}
    doc = {
        "step": "r2", "a16": {"pass": a16.get("pass"), "exit": rc, "mode": a16.get("mode"),
                              "baseline_ref": a16.get("baseline_ref"),
                              "baseline_digests_sha256": a16.get("baseline_digests_sha256")},
        "a16_baseline_expected": v4.A16_BASELINE,
        "intended_change_set": run_suite.A16_INTENDED_CHANGE,
        "trace_fixture_guard": guard, "trace_fixture_python_bitwise": bitwise,
        "trace_fixture_python_bitwise_all": all(bitwise.values()),
    }  # fmt: skip
    doc["pass"] = bool(rc == 0 and a16.get("pass") and a16.get("baseline_ref") == v4.A16_BASELINE
                       and run_suite.A16_INTENDED_CHANGE is None and guard["pass"])  # fmt: skip
    return finish8(doc, None, None, False)


# -- T13-ion (hr8, CUDA) ---------------------------------------------------------------------------
def step_t13_ion(a: argparse.Namespace) -> int:
    """T13-ion: C 290 MeV/u on warp-cuda, chunks 2^10 vs 2^18 histories, float32 and float64,
    same seed per precision: counters, tallies and deposit grids bit-identical within each
    precision (``parity.compare_partition``)."""
    from ionmc.transport.parity import compare_partition

    proj, e_u = CARBON12, 290.0
    r = csda_mm(proj, e_u)
    geo, _ = water_box(1.1 * r)
    nz = int(1.1 * r / 2)
    grid = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, nz), name="dose")
    n = scaled(T13_N, a.scale)
    out: dict[str, Any] = {}
    ok = True
    for i, prec in enumerate(("float32", "float64")):
        runs = [Simulation(ion_config(proj, e_u, geometry=geo, scoring=(grid,), backend="warp-cuda",
                                      precision=prec, seed=row_seed("t13-ion", i, a.seed_base),
                                      n=n, chunk=c)).run() for c in (2**10, 2**18)]  # fmt: skip
        v = compare_partition(runs[0], runs[1])
        v["valid"] = all(clean(x) for x in runs)
        out[f"chunks_1024_vs_262144_{prec}"] = v
        ok &= bool(v["pass"] and v["valid"])
    return finish8({"step": "t13-ion", "t13": out, "pass": ok}, T13_N, n, n < T13_N)


# == sharded LV/HR rows: I9, I12, I14, I15, T12-ion (amendment 4c) ===============================
# Amendment 4c (host budget correction): 2e5 warp-cpu histories per configuration, run as shards
# of <= 1500 s (separate host steps, ``shard``) whose hash-sealed partial documents are combined
# by one ``combine`` step per row (the lv5 v2-shard / v2-combine pattern; ``steps_v5`` digest,
# manifest and attestation helpers). Shard sizes follow the throughput measured in the sandbox
# (2026-11-11, C 290 MeV/u, all EM on, warp-cpu float64, one process, histories/s):
# s_max/f_E 1/0.02: 111, 1/0.005: 35, 0.5/0.02: 78, 0.5/0.005: 34, 0.1/0.02: 23, 0.1/0.005: 18
# (He 150 ~150, O 350 ~100 at 1/0.02); every shard is planned at <= ~1000 s. ``k`` of a shard is
# its position in ``shard_plan(row)`` (units in order, shards within a unit).
SHARD_BATCHES = 20
C290_SHARD_N: dict[tuple[float, float], int] = {
    (0.1, 0.005): 12_500, (0.1, 0.02): 20_000, (0.5, 0.005): 25_000, (0.5, 0.02): 50_000,
    (1.0, 0.005): 25_000, (1.0, 0.02): 100_000, (5.0, 0.02): 100_000,
}  # fmt: skip
"""Histories per shard of a C 290 configuration (s_max mm, f_E) at <= ~1000 s planned."""
C290_RATE: dict[tuple[float, float], float] = {
    (0.1, 0.005): 18.0, (0.1, 0.02): 23.0, (0.5, 0.005): 34.0, (0.5, 0.02): 78.0,
    (1.0, 0.005): 35.0, (1.0, 0.02): 111.0, (5.0, 0.02): 111.0,
}  # fmt: skip
"""Measured sandbox throughput [histories/s] (5 mm assumed no faster than 1 mm)."""
ROW_N = 200_000  # amendment 4c: warp-cpu histories per configuration
I9_FRACS, I9_BIN_MM, I9_TOL = (0.25, 0.5, 0.75, 0.9), 0.1, 0.02
I12_SMAX, I12_THETA_TOL, I12_SIGMA_TOL = (0.1, 0.5, 1.0, 5.0), 0.005, 0.01
I14_CONFIGS = ((0.1, 0.005), (0.1, 0.02), (0.5, 0.005), (0.5, 0.02), (1.0, 0.005), (1.0, 0.02))
I14_REFERENCE, I14_R80_TOL, I14_BATCHES = (1.0, 0.02), 0.05, 40
I15_CONFIGS = ((1.0, False), (1.0, True), (2.0, False), (2.0, True))  # (voxel mm, half shift)
T12_HALF_WIDTH_MM, T12_LATERAL_BIN_MM = 20.0, 0.2
T12_MAX_STEP_MM = 1.0
"""Step limit of the T12-ion samples: the lv8 default of every ion row (1 mm), NOT the proton
T12 default (``max_step = lateral bin = 0.2 mm``); the comparison is between backends with the
same configuration (backend parity), see the report of this implementation."""


def _sname(s: float) -> str:
    return {0.1: "s01", 0.5: "s05", 1.0: "s1", 5.0: "s5"}[s]


def _fname(f: float) -> str:
    return {0.005: "f0005", 0.02: "f002"}[f]


# unit -> (ion, MeV/u, backend, precision, frozen total histories, shards, batches per shard)
T12_SAMPLES: dict[str, dict[str, tuple[str, float, str, str, int, int, int]]] = {
    "t12-ion-lv": {
        "python": ("carbon", 290.0, "python", "float64", 4000, 2, 20),
        "cpu64": ("carbon", 290.0, "warp-cpu", "float64", ROW_N, 2, 20),
    },
    "t12-ion-hr": {  # the CPU parts run inside hr8 (seed base 20591004), not imported from lv8
        "c290-cpu64": ("carbon", 290.0, "warp-cpu", "float64", ROW_N, 2, 20),
        "c290-cuda32": ("carbon", 290.0, "warp-cuda", "float32", 1_000_000, 1, 100),
        "c290-cuda64": ("carbon", 290.0, "warp-cuda", "float64", 1_000_000, 1, 100),
        "he150-cpu64": ("he", 150.0, "warp-cpu", "float64", ROW_N, 2, 20),
        "he150-cuda32": ("he", 150.0, "warp-cuda", "float32", 1_000_000, 1, 100),
        "o350-cpu64": ("oxygen", 350.0, "warp-cpu", "float64", ROW_N, 2, 20),
        "o350-cuda32": ("oxygen", 350.0, "warp-cuda", "float32", 1_000_000, 1, 100),
    },
}
T12_PAIRS: dict[str, tuple[tuple[str, str], ...]] = {
    "t12-ion-lv": (("python", "cpu64"),),
    "t12-ion-hr": (("c290-cpu64", "c290-cuda32"), ("c290-cuda32", "c290-cuda64"),
                   ("he150-cpu64", "he150-cuda32"), ("o350-cpu64", "o350-cuda32")),
}  # fmt: skip
"""Compared pairs; ``c290-cuda32`` vs ``c290-cuda64`` is the f32-vs-f64 probe."""
STANDARD_UNITS = {"he150": ("he", 150.0), "c290": ("carbon", 290.0), "o350": ("oxygen", 350.0)}


def shard_plan(row: str) -> list[tuple[str, int, int, int, dict[str, Any]]]:
    """``(unit, shards, histories per shard, batches per shard, meta)`` of a sharded row, in the
    order that defines ``k`` (``row_seed(row, k)``)."""
    if row == "i9":
        return [(u, 2, ROW_N // 2, SHARD_BATCHES, {"point": u}) for u in STANDARD_UNITS]
    if row == "i12":
        return [(f"{obs}-{_sname(s)}", ROW_N // C290_SHARD_N[s, 0.02], C290_SHARD_N[s, 0.02],
                 SHARD_BATCHES, {"obs": obs, "smax": s})
                for obs in ("theta", "sigma") for s in I12_SMAX]  # fmt: skip
    if row == "i14":
        return [(f"{_sname(s)}-{_fname(f)}", ROW_N // C290_SHARD_N[s, f], C290_SHARD_N[s, f],
                 SHARD_BATCHES, {"smax": s, "frac": f}) for s, f in I14_CONFIGS]  # fmt: skip
    if row == "i15":  # 1 mm voxels add boundary steps: their sigma shards are halved
        return [(f"{obs}-v{v:g}-sh{int(sh)}", n_sh, ROW_N // n_sh, SHARD_BATCHES,
                 {"obs": obs, "voxel": v, "shift": sh})
                for obs in ("theta", "sigma") for v, sh in I15_CONFIGS
                for n_sh in [4 if (obs == "sigma" and v == 1.0) else 2]]  # fmt: skip
    if row in T12_SAMPLES:
        return [(u, sh, n // sh, b, {"sample": u}) for u, (_i, _e, _b, _p, n, sh, b)
                in T12_SAMPLES[row].items()]  # fmt: skip
    raise SystemExit(f"unknown sharded row {row!r}")


def shard_k(row: str, unit: str, shard: int) -> int:
    k = 0
    for u, n_sh, _n, _b, _m in shard_plan(row):
        if u == unit:
            if not 0 <= shard < n_sh:
                raise SystemExit(f"{row} {unit}: shard must be in 0..{n_sh - 1}")
            return k + shard
        k += n_sh
    raise SystemExit(f"{row}: unknown unit {unit!r}")


def planned_shard_seconds(row: str, unit: str) -> float:
    """Planned wall time [s] of one shard at frozen scale (sandbox throughput above; python
    0.5 s per history; CUDA unknown before the first hr8 step, reported as 0)."""
    _u, _sh, n, _b, m = next(x for x in shard_plan(row) if x[0] == unit)
    if row == "i9":
        return n / {"he150": 150.0, "c290": 111.0, "o350": 100.0}[m["point"]]
    if row in ("i12", "i15"):
        rate = C290_RATE[m.get("smax", 1.0), 0.02] * (0.7 if m.get("voxel") == 1.0 else 1.0)
        return n / rate * (0.55 if m["obs"] == "theta" else 1.0)
    if row == "i14":
        return n / C290_RATE[m["smax"], m["frac"]]
    ion, _e, backend, _p, _t, _s, _b = T12_SAMPLES[row][unit]
    if backend == "warp-cuda":
        return 0.0
    if backend == "python":
        return 0.5 * n
    return n / {"he": 150.0, "carbon": 111.0, "oxygen": 100.0}[ion]


def ion_path(proj: Projectile, e_u: float) -> Any:
    """The CSDA path of ``proj`` in water (``ionmc.transport.ion_mcs_checks.IonPath``, the path of
    the I10 CI test): the implemented T_ion per mm in absolute scale (``power``; the shared
    ``mcs_variance`` with the material 1/X_S = ``inverse_scattering_length_cm2_per_g(WATER)`` and
    rho = ``WATER.density_g_cm3``, the values the transport uses) and its quadrature
    ``theta2_quadrature``; also the integrand of the Fermi-Eyges A2(z)
    (``steps.fermi_eyges_a2``)."""
    from ionmc.transport.ion_mcs_checks import IonPath

    return IonPath(proj, e_u)


def _lateral_run(proj: Projectile, e_u: float, fracs: tuple[float, ...], bin_mm: float, *,
                 voxel: float, shift: bool, smax: float, n: int, seed: int, batches: int
                 ) -> tuple[Any, dict[str, Any]]:  # fmt: skip
    """Deposits in 1 mm slabs at ``floor(f R)`` with ``bin_mm`` lateral bins (MCS on, straggling
    off; scoring independent of the transport voxels; the T7/T8/T14 observable with the ion)."""
    path = ion_path(proj, e_u)
    r = path.r1_mm
    zs = {f: math.floor(f * r / base.SLAB_MM) * base.SLAB_MM for f in fracs}
    geo, _ = base._slab_geometry(max(zs.values()) + base.SLAB_MM + 3.0, voxel, shift, shift)
    w = math.ceil(8.0 * math.sqrt(base.fermi_eyges_a2(path, 0.9 * r)) / bin_mm) * bin_mm
    nl = int(round(2 * w / bin_mm))
    grids = tuple(ScoringGrid((-w, -w, zs[f]), (bin_mm, bin_mm, base.SLAB_MM), (nl, nl, 1),
                              name=f"slab{int(round(100 * f))}") for f in fracs)  # fmt: skip
    cfg = ion_config(proj, e_u, geometry=geo, scoring=grids, seed=seed, n=n, n_batches=batches,
                     straggling=False, max_step=smax)  # fmt: skip
    res = Simulation(cfg).run()
    slabs = {}
    for f, g in zip(fracs, grids, strict=True):
        e = np.asarray(res.grid(g.name).energy_mev)[:, :, 0]
        slabs[str(f)] = {"px": e.sum(axis=1).tolist(), "py": e.sum(axis=0).tolist(), "x0_mm": -w,
                         "bin_mm": bin_mm, "z0_mm": zs[f]}  # fmt: skip
    return res, {"slabs": slabs, "csda_range_mm": r}


def _theta_run(proj: Projectile, e_u: float, *, voxel: float, shift: bool, smax: float, n: int,
               seed: int, batches: int) -> tuple[Any, dict[str, Any]]:  # fmt: skip
    """Exit angles of the ions leaving a 0.5 R water slab (exit plane exact for every voxel
    size and shift; MCS on, straggling off), from the escape records (the T8/T14 observable)."""
    x = 0.5 * csda_mm(proj, e_u)
    geo, z_end = base._slab_geometry(x, voxel, shift, shift, exit_clip=True)
    h = base.HALF_WORLD_MM
    grid = ScoringGrid((-h, -h, 0.0), (2 * h, 2 * h, max(x, 1.0)), (1, 1, 1), name="slab")
    cfg = ion_config(proj, e_u, geometry=geo, scoring=(grid,), seed=seed, n=n, n_batches=batches,
                     straggling=False, max_step=smax,
                     diag=DiagnosticsOptions(escape_records=True))  # fmt: skip
    res = Simulation(cfg).run()
    u = np.asarray(res.diagnostics["escape_direction"])
    tx, ty = u[:, 0] / u[:, 2], u[:, 1] / u[:, 2]
    return res, {"sum_theta2": float(np.sum(tx * tx + ty * ty)), "n_escaped": int(len(u)),
                 "z_end_mm": z_end, "thickness_mm": x}  # fmt: skip


def t12_ion_layout(proj: Projectile, e_u: float) -> tuple[Any, tuple[ScoringGrid, ...], Any]:
    """Water box (1.3 R deep) and the frozen T12 scoring grids (``parity.t12_scoring_grids``:
    1 mm IDD, lateral slabs at z/R 0.25, 0.5, 0.9 with 0.2 mm bins) at the ion's CSDA range."""
    from ionmc.transport.parity import t12_scoring_grids

    r = csda_mm(proj, e_u)
    grids, layout = t12_scoring_grids(r, half_width_mm=T12_HALF_WIDTH_MM,
                                      lateral_bin_mm=T12_LATERAL_BIN_MM)  # fmt: skip
    m = T12_HALF_WIDTH_MM + 10.0
    return BoxPhantom((-m, -m, 0.0), (2 * m, 2 * m, 1.3 * r), WATER), grids, layout


def _run_unit(row: str, meta: dict[str, Any], n: int, seed: int, b: int
              ) -> tuple[Any, dict[str, Any]]:  # fmt: skip
    if row == "i9":
        ion, e_u = STANDARD_UNITS[meta["point"]]
        return _lateral_run(IONS[ion], e_u, I9_FRACS, I9_BIN_MM, voxel=-1.0, shift=False,
                            smax=1.0, n=n, seed=seed, batches=b)  # fmt: skip
    if row in ("i12", "i15"):
        v, sh = meta.get("voxel", -1.0), meta.get("shift", False)
        smax = meta["smax"] if row == "i12" else min(1.0, v)
        if meta["obs"] == "theta":
            return _theta_run(CARBON12, 290.0, voxel=v, shift=sh, smax=smax, n=n, seed=seed,
                              batches=b)  # fmt: skip
        fr = (0.9,) if row == "i12" else (0.5, 0.9)
        return _lateral_run(CARBON12, 290.0, fr, base.LATERAL_BIN_MM, voxel=v, shift=sh,
                            smax=smax, n=n, seed=seed, batches=b)  # fmt: skip
    if row == "i14":
        r = csda_mm(CARBON12, 290.0)
        depth = 1.3 * r
        geo = BoxPhantom((-50.0, -50.0, 0.0), (100.0, 100.0, depth), WATER)
        grid = ScoringGrid((-40.0, -40.0, 0.0), (80.0, 80.0, 1.0), (1, 1, int(math.ceil(depth))),
                           name="idd")  # fmt: skip
        cfg = ion_config(CARBON12, 290.0, geometry=geo, scoring=(grid,), seed=seed, n=n,
                         n_batches=b, max_step=meta["smax"], frac=meta["frac"])  # fmt: skip
        res = Simulation(cfg).run()
        idd = np.asarray(res.grid("idd").batch_energy_mev).reshape(b, -1)
        return res, {"idd_batches": idd.tolist(), "csda_range_mm": r}
    ion, e_u, backend, prec, _t, _s, _b = T12_SAMPLES[row][meta["sample"]]
    from ionmc.transport.parity import t12_observables

    geo, grids, layout = t12_ion_layout(IONS[ion], e_u)
    cfg = ion_config(IONS[ion], e_u, geometry=geo, scoring=grids, backend=backend, precision=prec,
                     seed=seed, n=n, n_batches=b, max_step=T12_MAX_STEP_MM)  # fmt: skip
    res = Simulation(cfg).run()
    obs = t12_observables(res, layout)
    return res, {"arrays": {k: v.tolist() for k, v in obs.arrays.items()},
                 "scalars": {k: v.tolist() for k, v in obs.scalars.items()},
                 "n_batches": obs.n_batches, "histories_per_batch": obs.histories_per_batch,
                 "precision": obs.precision, "backend": backend}  # fmt: skip


def partial_bindings(a: argparse.Namespace) -> dict[str, Any]:
    """What a v8 partial is bound to: the run (SHA, suite), the seed base and the scale."""
    import os

    sha = os.environ.get("IONMC_RUN_SHA", "unknown")
    return {"format": FORMAT, "git_sha": sha, "run_sha": sha,
            "suite": os.environ.get("IONMC_RUN_SUITE", "unknown"), "seed_base": a.seed_base,
            "scale": a.scale}  # fmt: skip


def write_partial8(a: argparse.Namespace, name: str, doc: dict[str, Any]) -> str:
    """Hash-sealed partial (``steps_v5.content_digest``) plus the ``PARTIAL <name> <sha>`` stdout
    line that ``run_suite.echo_partial_lines`` forwards (the protected record of the shard)."""
    import steps_v5 as v5

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    doc = {**partial_bindings(a), **doc}
    doc["content_sha256"] = v5.content_digest(doc)
    path = out / f"{name}.json"
    path.write_text(json.dumps(doc, sort_keys=True, default=base._json))
    print(f"PARTIAL {path.name} {doc['content_sha256']}", flush=True)
    return str(path)


def load_partials8(a: argparse.Namespace, names: list[str]) -> list[dict[str, Any]]:
    """``steps_v5.load_partials`` with the v8 bindings: each partial must exist exactly once in
    ``--dirs`` (or their ``samples``), match its ``content_sha256`` and the run bindings; an
    imported one (not in ``dirs[0]``) must be in ``--partials-manifest`` with its digest and a
    host run id. Any defect stops the step (exit != 0)."""
    import steps_v5 as v5

    imported = len(a.dirs) > 1
    manifest: dict[str, dict[str, str]] = {}
    if imported:
        manifest, msha = v5.read_manifest(getattr(a, "partials_manifest", None))
        a.attested_manifest_path, a.attested_manifest_sha256 = str(a.partials_manifest), msha
    a.attested_partials = getattr(a, "attested_partials", [])
    out = []
    for nm in names:
        hits = [(i, p) for i, d in enumerate(a.dirs)
                for p in (Path(d) / nm, Path(d) / "samples" / nm) if p.is_file()]  # fmt: skip
        if len(hits) != 1:
            raise SystemExit(f"partial {nm}: found {len(hits)} copies in {a.dirs} (need exactly 1)")
        idx, hit = hits[0]
        doc = json.loads(hit.read_text())
        digest = v5.content_digest(doc)
        if doc.get("content_sha256") != digest:
            raise SystemExit(f"partial {nm}: content_sha256 does not match the document")
        rec: dict[str, Any] = {"name": nm, "content_sha256": digest, "origin": "current-run",
                               "host_run_id": None}  # fmt: skip
        if idx != 0:
            if nm not in manifest or manifest[nm]["content_sha256"] != digest:
                raise SystemExit(f"partial {nm}: imported but not in the manifest with its digest")
            rec.update(origin="imported", host_run_id=manifest[nm]["host_run_id"])
        for key, want in partial_bindings(a).items():
            if doc.get(key) != want:
                raise SystemExit(f"partial {nm}: {key} is {doc.get(key)!r}, expected {want!r}")
        a.attested_partials.append(rec)
        out.append(doc)
    return out


def step_shard(a: argparse.Namespace) -> int:
    """One shard of a sharded row (``--row``, ``--unit``, ``--shard``): runs the unit's
    configuration with ``seed = row_seed(row, k)`` and writes the sealed partial
    ``<row>-<unit>-s<shard>.json``. Invalid results and nonzero counters fail the step."""
    k = shard_k(a.row, a.unit, a.shard)
    _u, _sh, n_nom, b, meta = next(x for x in shard_plan(a.row) if x[0] == a.unit)
    n = scaled(n_nom, a.scale, b)
    seed = row_seed(a.row, k, a.seed_base)
    t0 = time.perf_counter()
    res, payload = _run_unit(a.row, meta, n, seed, b)
    wall = time.perf_counter() - t0
    ok = clean(res)
    path = write_partial8(a, f"{a.row}-{a.unit}-s{a.shard}", {
        "row": a.row, "unit": a.unit, "shard": a.shard, "k": k, "seed": seed, "n": n,
        "n_nominal": n_nom, "batches": b, "valid": ok, "counters": res.counters.as_dict(),
        **payload,
    })  # fmt: skip
    doc = {"step": f"{a.row}-shard", "row": a.row, "unit": a.unit, "shard": a.shard, "k": k,
           "seed": seed, "partial": path, "wall_s": wall, "hist_per_s": n / wall,
           "planned_wall_s_frozen": planned_shard_seconds(a.row, a.unit),
           "counters": res.counters.as_dict(), "valid": bool(res.valid), "pass": ok}  # fmt: skip
    return finish8(doc, n_nom, n, n < n_nom)


def _unit_parts(a: argparse.Namespace) -> dict[str, list[dict[str, Any]]]:
    """Every shard partial of ``a.row`` by unit, verified (seed, k, n, validity); fail closed."""
    out: dict[str, list[dict[str, Any]]] = {}
    for unit, n_sh, n_nom, b, _m in shard_plan(a.row):
        parts = load_partials8(a, [f"{a.row}-{unit}-s{j}.json" for j in range(n_sh)])
        for j, p in enumerate(parts):
            k = shard_k(a.row, unit, j)
            if (p.get("row"), p.get("unit"), p.get("shard"), p.get("k")) != (a.row, unit, j, k) or (
                p.get("seed") != row_seed(a.row, k, a.seed_base)
                or p.get("n") != scaled(n_nom, a.scale, b) or p.get("valid") is not True
            ):  # fmt: skip
                raise SystemExit(f"{a.row} {unit} shard {j}: identity, seed, count or validity")
        out[unit] = parts
    return out


def _pooled_sigma2(parts: list[dict[str, Any]], frac: str) -> float:
    """Sheppard-corrected projected variance [mm^2] of the n-weighted pooled deposit profiles."""
    n = np.array([p["n"] for p in parts], dtype=float)
    px = sum(w * np.array(p["slabs"][frac]["px"]) for w, p in zip(n, parts, strict=True)) / n.sum()
    py = sum(w * np.array(p["slabs"][frac]["py"]) for w, p in zip(n, parts, strict=True)) / n.sum()
    s = parts[0]["slabs"][frac]
    xc = s["x0_mm"] + (np.arange(len(px)) + 0.5) * s["bin_mm"]
    v = 0.5 * (base._profile_var(px, xc) + base._profile_var(py, xc))
    return float(v - s["bin_mm"] ** 2 / 12.0)


def _pooled_theta(parts: list[dict[str, Any]]) -> dict[str, Any]:
    ne = sum(p["n_escaped"] for p in parts)
    return {"theta_rms": math.sqrt(sum(p["sum_theta2"] for p in parts) / (2.0 * ne)),
            "n_escaped": ne, "histories": sum(p["n"] for p in parts),
            "z_end_mm": parts[0]["z_end_mm"]}  # fmt: skip


def combine_i9(units: dict[str, list[dict[str, Any]]]) -> tuple[dict[str, Any], bool]:
    """I9: Sheppard-corrected lateral variance at z/R in I9_FRACS against the Fermi-Eyges A2 at
    the slab centre with the implemented T_ion: |var / A2 - 1| <= 2 % (the row's literal
    'lateral variance within 2 % of A2'; sigma ratios reported). A2 integrates ``ion_path``'s
    T_ion with the material 1/X_S of water (``inverse_scattering_length_cm2_per_g``)."""
    out, ok = {}, True
    for unit, parts in units.items():
        ion, e_u = STANDARD_UNITS[unit]
        path = ion_path(IONS[ion], e_u)
        rows = {}
        for f in I9_FRACS:
            v = _pooled_sigma2(parts, str(f))
            zc = parts[0]["slabs"][str(f)]["z0_mm"] + 0.5 * base.SLAB_MM
            a2 = base.fermi_eyges_a2(path, zc)
            p = bool(abs(v / a2 - 1.0) <= I9_TOL)
            rows[str(f)] = {"variance_mm2": v, "fermi_eyges_a2_mm2": a2, "z_centre_mm": zc,
                            "variance_over_a2": v / a2, "sigma_over_sqrt_a2": math.sqrt(v / a2),
                            "pass": p}  # fmt: skip
            ok &= p
        out[unit] = {"slabs": rows, "histories": sum(p["n"] for p in parts), "tolerance": I9_TOL}
    return out, ok


def combine_i12(units: dict[str, list[dict[str, Any]]]) -> tuple[dict[str, Any], bool]:
    """I12: pairwise theta_rms spread <= 0.5 %, sigma(0.9 R) spread <= 1 %, every theta_rms
    within 0.5 % of the I10 quadrature (``theta2_quadrature`` at the exact exit plane; T_ion with
    the material 1/X_S of water, ``inverse_scattering_length_cm2_per_g``)."""
    path = ion_path(CARBON12, 290.0)
    th = {s: _pooled_theta(units[f"theta-{_sname(s)}"]) for s in I12_SMAX}
    sg = {s: math.sqrt(_pooled_sigma2(units[f"sigma-{_sname(s)}"], "0.9")) for s in I12_SMAX}
    quad = math.sqrt(path.theta2_quadrature(th[I12_SMAX[0]]["z_end_mm"]))
    t = [v["theta_rms"] for v in th.values()]
    spread_t, spread_s = max(t) / min(t) - 1.0, max(sg.values()) / min(sg.values()) - 1.0
    dev_q = max(abs(x / quad - 1.0) for x in t)
    ok = bool(spread_t <= I12_THETA_TOL and spread_s <= I12_SIGMA_TOL and dev_q <= I12_THETA_TOL)
    return {"theta": {str(s): v for s, v in th.items()},
            "sigma_0.9R_mm": {str(s): v for s, v in sg.items()},
            "theta_quadrature": quad, "theta_pairwise_spread": spread_t,
            "sigma_pairwise_spread": spread_s, "theta_max_dev_from_quadrature": dev_q,
            "tolerances": {"theta": I12_THETA_TOL, "sigma": I12_SIGMA_TOL}}, ok  # fmt: skip


def rebatch(x: np.ndarray, b: int) -> np.ndarray:
    """Means of consecutive groups of equal-size batches (per-primary batch means stay exact)."""
    if x.shape[0] % b:
        raise SystemExit(f"cannot rebatch {x.shape[0]} batches into {b}")
    return x.reshape(b, x.shape[0] // b, *x.shape[1:]).mean(axis=1)


def combine_i14(units: dict[str, list[dict[str, Any]]]) -> tuple[dict[str, Any], bool]:
    """I14 (T9 pattern, ``steps.step_t9``): every configuration's pooled IDD (rebatched to 40
    batches) against the reference (s_max 1 mm, f_E 0.02) by ``parity.t12_compare`` (grouped
    chi-square AND the max-T statistic of footnote 17, inconclusive = fail); |Delta R80| <= 0.05 mm
    over ALL pairs of configurations (R80 = mean of the batch R80s, as in T9)."""
    from ionmc.transport import parity

    r = units[f"{_sname(1.0)}-{_fname(0.02)}"][0]["csda_range_mm"]
    layout = parity.T12Layout("idd", (), r, 1.0, 80.0)
    obs = {}
    for unit, parts in units.items():
        idd = rebatch(np.concatenate([np.array(p["idd_batches"]) for p in parts]), I14_BATCHES)
        hpb = sum(p["n"] for p in parts) // I14_BATCHES
        obs[unit] = parity.t12_observables_from_grids({"idd": idd[:, None, None, :]}, (), layout,
                                                      hpb)  # fmt: skip
    ref = f"{_sname(I14_REFERENCE[0])}-{_fname(I14_REFERENCE[1])}"
    r80 = {u: float(o.scalars["r80_mm"].mean()) for u, o in obs.items()}
    dr = max(r80.values()) - min(r80.values())
    out: dict[str, Any] = {"r80_mm": r80, "max_pairwise_delta_r80_mm": dr,
                           "tolerance_mm": I14_R80_TOL, "reference": ref,
                           "comparisons": {}}  # fmt: skip
    ok = bool(dr <= I14_R80_TOL)
    for unit, o in obs.items():
        if unit != ref:
            v = parity.t12_compare(obs[ref], o)["arrays"]["idd"]
            out["comparisons"][unit] = v
            ok &= bool(v["pass"])
    return out, ok


def combine_i15(units: dict[str, list[dict[str, Any]]]) -> tuple[dict[str, Any], bool]:
    """I15 (reduced T14): pairwise theta_rms spread <= 0.5 % and sigma spread <= 1 % at 0.5 R
    and 0.9 R over the four voxel configurations; the ratios to the I10 quadrature and to the
    Fermi-Eyges sigma are informative (both from ``ion_path``: T_ion with the material 1/X_S of
    water, ``inverse_scattering_length_cm2_per_g``)."""
    path = ion_path(CARBON12, 290.0)
    rows: dict[str, Any] = {}
    for v, sh in I15_CONFIGS:
        key = f"v{v:g}-sh{int(sh)}"
        th = _pooled_theta(units[f"theta-{key}"])
        quad = math.sqrt(path.theta2_quadrature(th["z_end_mm"]))
        th["theta_over_quadrature"] = th["theta_rms"] / quad
        sp = units[f"sigma-{key}"]
        sg = {f: math.sqrt(_pooled_sigma2(sp, f)) for f in ("0.5", "0.9")}
        zc = {f: sp[0]["slabs"][f]["z0_mm"] + 0.5 * base.SLAB_MM for f in sg}
        fe = {f: math.sqrt(base.fermi_eyges_a2(path, zc[f])) for f in sg}
        rows[key] = {"theta": th, "sigma_mm": sg,
                     "sigma_over_fermi_eyges": {f: sg[f] / fe[f] for f in sg}}  # fmt: skip
    t = [x["theta"]["theta_rms"] for x in rows.values()]
    sp5 = [x["sigma_mm"]["0.5"] for x in rows.values()]
    sp9 = [x["sigma_mm"]["0.9"] for x in rows.values()]
    res = {"theta_pairwise_spread": max(t) / min(t) - 1.0,
           "sigma_pairwise_spread_0.5R": max(sp5) / min(sp5) - 1.0,
           "sigma_pairwise_spread_0.9R": max(sp9) / min(sp9) - 1.0}  # fmt: skip
    worst_sigma = max(res["sigma_pairwise_spread_0.5R"], res["sigma_pairwise_spread_0.9R"])
    ok = bool(res["theta_pairwise_spread"] <= I12_THETA_TOL and worst_sigma <= I12_SIGMA_TOL)
    return {"runs": rows, **res}, ok


def combine_t12(row: str, units: dict[str, list[dict[str, Any]]]) -> tuple[dict[str, Any], bool]:
    """T12-ion: the shards of each sample concatenated batch-wise (equal histories per batch
    required) and every pair of ``T12_PAIRS`` compared by ``parity.t12_compare`` (the frozen
    T12 rules, footnotes 9-22; inconclusive = fail)."""
    from ionmc.transport import parity

    obs = {}
    for unit, parts in units.items():
        hpb = {p["histories_per_batch"] for p in parts}
        if len(hpb) != 1:
            raise SystemExit(f"{row} {unit}: shards differ in histories per batch")
        cat = {f: {k: np.concatenate([np.array(p[f][k]) for p in parts]) for k in parts[0][f]}
               for f in ("arrays", "scalars")}  # fmt: skip
        obs[unit] = parity.T12Observables(cat["arrays"], cat["scalars"],
                                          sum(p["n_batches"] for p in parts), hpb.pop(),
                                          parts[0]["precision"])  # fmt: skip
    out, ok = {}, True
    for x, y in T12_PAIRS[row]:
        v = parity.t12_compare(obs[x], obs[y])
        out[f"{x}_vs_{y}"] = v
        ok &= bool(v["pass"])
    return out, ok


def step_combine(a: argparse.Namespace) -> int:
    """Combine step of a sharded row: loads every shard partial (fail closed on a missing,
    altered, foreign, invalid or nonzero-counter shard) and applies the row's criterion."""
    import steps_v5 as v5

    units = _unit_parts(a)
    if a.row in T12_SAMPLES:
        out, ok = combine_t12(a.row, units)
    else:
        out, ok = {"i9": combine_i9, "i12": combine_i12, "i14": combine_i14,
                   "i15": combine_i15}[a.row](units)  # fmt: skip
    frozen = sum(n_sh * n for _u, n_sh, n, _b, _m in shard_plan(a.row))
    used = sum(p["n"] for parts in units.values() for p in parts)
    doc = {"step": f"{a.row}-combine", "row": a.row, a.row.replace("-", "_"): out,
           "shards": {u: [p["seed"] for p in ps] for u, ps in units.items()},
           "attestation": v5.attestation_block(a), "pass": bool(ok)}  # fmt: skip
    return finish8(doc, frozen, used, used < frozen)


SHARDED_ROWS = ("i9", "i12", "i14", "i15", "t12-ion-lv", "t12-ion-hr")


STEPS = {
    "td-tl-sample": step_td_tl_sample,
    "ss-sample": step_ss_sample,
    "i6-lv": step_i6,
    "i8": step_i8,
    "i13": step_i13,
    "t1-ion": step_t1_ion,
    "r1-ion": step_r1_ion,
    "r2": step_r2,
    "t13-ion": step_t13_ion,
    "shard": step_shard,
    "combine": step_combine,
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument(
        "--case",
        type=int,
        default=None,
        help="td-tl-sample: case 1..12; ss-sample: a V5-SS case of SS_CASES",
    )
    ap.add_argument("--k", type=int, default=0, help="td-tl-sample / ss-sample: chunk 0..99")
    ap.add_argument("--point", type=int, default=0, help="i8: 0 He 150, 1 C 290, 2 O 350")
    ap.add_argument("--case-index", type=int, default=0, help="i13: index 0..8 of I13_CASES")
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--seed-base", type=int, default=None,
                    help="lv8 base 20581004 (re-run base 20621004); hr8 steps default to 20591004")
    ap.add_argument("--timeout", type=float, default=None, help="accepted for run_suite uniformity")
    ap.add_argument("--out-dir", default="samples")
    ap.add_argument("--row", choices=SHARDED_ROWS, default=None, help="shard/combine: the row")
    ap.add_argument("--unit", default=None, help="shard: the unit of shard_plan(row)")
    ap.add_argument("--shard", type=int, default=0, help="shard: index within the unit")
    ap.add_argument("--dirs", nargs="+", default=None,
                    help="combine: the current output directory first, then imported archives")
    ap.add_argument("--partials-manifest", default=None,
                    help="combine with imported directories: steps_v5.read_manifest schema")
    args = ap.parse_args(argv)
    if args.seed_base is None and args.step == "ss-sample":
        args.seed_base = SS_SEED_BASE
    if args.seed_base is None:
        hr = args.step == "t13-ion" or args.row == "t12-ion-hr"
        args.seed_base = HR8_BASE if hr else LV8_BASE
    if args.step in ("shard", "combine") and args.row is None:
        raise SystemExit("shard/combine need --row")
    if args.step == "shard":
        shard_k(args.row, args.unit, args.shard)  # range check
    if args.step == "combine" and not args.dirs:
        raise SystemExit("combine needs --dirs")
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    if args.step == "td-tl-sample":
        if args.case not in CASES:
            raise SystemExit("--case must be in 1..12")
        seed_of(args.case, args.k, args.seed_base)  # range check
    if args.step == "ss-sample":
        if args.case not in SS_R_INDEX:
            raise SystemExit(f"--case must be one of {SS_CASES}")
        ss_seed(args.case, args.k, args.seed_base)  # range check
    if not 0 <= args.point < len(STANDARD) or not 0 <= args.case_index < len(I13_CASES):
        raise SystemExit("--point in 0..2, --case-index in 0..8")
    base.SEED_BASE = args.seed_base
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
