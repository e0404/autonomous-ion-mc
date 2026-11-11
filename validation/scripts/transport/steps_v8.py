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

Verification rows of plan section A (suites ``lv8`` and ``hr8`` of ``run_suite.py``)::

    steps_v8.py {i6-lv,i8,i13,t1-ion,r1-ion,r2,t13-ion} [--point P] [--case-index I] [--scale F]
                [--seed-base N]

Implemented: I6 (LV part), I8 + I8b (one step per standard point), I13 (one step per case),
T1-ion LV, R1-ion, R2 LV (lv8) and T13-ion (hr8). NOT implemented in this revision: I9, I12, I14,
I15 and T12-ion (LV and HR); their r_index are reserved in ``RESERVED_R_INDEX``. Measured in the
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
}
RESERVED_R_INDEX = {"i9": 44, "i12": 45, "i14": 46, "i15": 47, "t12-ion-lv": 49, "t12-ion-hr": 54}
"""r_index declared for the rows NOT implemented in this revision of ``steps_v8.py`` (I9, I12,
I14, I15, T12-ion): kept free so that their later implementation uses exactly these indices."""
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
    diag: DiagnosticsOptions | None = None, chunk: int | None = None,
) -> SimulationConfig:  # fmt: skip
    """Ion pencil beam on the axis at z = 0 (the source API of ``sample_config``), Bethe stopping
    with the default options, every other physics option at its default (in particular the
    end-of-range model ``"default"`` = ``continuous_residual_v1`` for ions)."""
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), e_u * proj.a),
        geometry=geometry,
        scoring=scoring,
        physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource(), straggling=straggling,
                               multiple_scattering=mcs, max_step_mm=max_step),  # fmt: skip
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


STEPS = {
    "td-tl-sample": step_td_tl_sample,
    "i6-lv": step_i6,
    "i8": step_i8,
    "i13": step_i13,
    "t1-ion": step_t1_ion,
    "r1-ion": step_r1_ion,
    "r2": step_r2,
    "t13-ion": step_t13_ion,
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--case", type=int, default=None, help="td-tl-sample: case 1..12")
    ap.add_argument("--k", type=int, default=0, help="td-tl-sample: chunk 0..99")
    ap.add_argument("--point", type=int, default=0, help="i8: 0 He 150, 1 C 290, 2 O 350")
    ap.add_argument("--case-index", type=int, default=0, help="i13: index 0..8 of I13_CASES")
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--seed-base", type=int, default=None,
                    help="lv8 base 20581004 (re-run base 20621004); hr8 steps default to 20591004")
    ap.add_argument("--timeout", type=float, default=None, help="accepted for run_suite uniformity")
    ap.add_argument("--out-dir", default="samples")
    args = ap.parse_args(argv)
    if args.seed_base is None:
        args.seed_base = HR8_BASE if args.step == "t13-ion" else LV8_BASE
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    if args.step == "td-tl-sample":
        if args.case not in CASES:
            raise SystemExit("--case must be in 1..12")
        seed_of(args.case, args.k, args.seed_base)  # range check
    if not 0 <= args.point < len(STANDARD) or not 0 <= args.case_index < len(I13_CASES):
        raise SystemExit("--point in 0..2, --case-index in 0..8")
    base.SEED_BASE = args.seed_base
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
