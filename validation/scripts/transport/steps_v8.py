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

The step writes ``ionmc-case<N>-k<K>.npz`` (per-grid batch energy arrays [MeV] and the grid edges)
and prints one JSON document (``#JSON-BEGIN`` / ``#JSON-END``) with the sha256 of the npz, the
effective-configuration summary, the counters, validity and the energy balance. Every invalid
result or nonzero counter fails the step. No TD/TL statistic is computed here.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

import steps as base
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, Projectile
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource

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
    np.savez(npz, **arrays)
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


STEPS = {"td-tl-sample": step_td_tl_sample}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--case", type=int, required=True, help="case 1..12 (r2-topas-target section 3)")
    ap.add_argument("--k", type=int, default=0, help="chunk index 0..99 (seed offset)")
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--seed-base", type=int, default=LV8_BASE, help="lv8 base (re-run base 20621004)")
    ap.add_argument("--out-dir", default="samples")
    args = ap.parse_args(argv)
    if args.case not in CASES:
        raise SystemExit("--case must be in 1..12")
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    seed_of(args.case, args.k, args.seed_base)  # range check
    base.SEED_BASE = args.seed_base
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
