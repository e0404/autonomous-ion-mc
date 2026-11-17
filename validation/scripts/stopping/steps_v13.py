"""V3-013 seal material: ionmc sample definitions and model-only predictions (suite ``lv13``).

Usage::

    python validation/scripts/stopping/steps_v13.py vc-sample  --case N --model M [--k K] [--scale F] [--out-dir DIR]
    python validation/scripts/stopping/steps_v13.py vc2d-sample --case N --model M --backend {python,warp-cpu}
                                                                [--k K] [--out-dir DIR]
    python validation/scripts/stopping/steps_v13.py predict [--out PATH]

REGISTERED, NOT RUN (plan ``validation/plans/v3-013-acceptance.md`` sections A.6, C and D steps 3-4).
Cases 1-3 are He 75 / 125 / 200 MeV/u in water; models are ``icru90-water-v1`` and the default
``bethe-ionmc-v1``.

``vc-sample`` is the V-C1 / V-C2 / V-C1-A sample of one case and model: warp-cpu float64, 1e6
histories in 40 batches, IDD scoring grid 1 x 1 x L/h_z over the full 120 x 120 mm field on the edges
of the TOPAS bundle of ``gen_v13_bundles.py``, MCS and straggling on, nuclear off, the default
end-of-range model, requested water I = 78 eV, ``stopping_model = M``. Seeds (plan section C, family
``ion13``): ``seed = 20612004 + 100 r_index + k`` with ``r_index = case - 1`` (0-2); the seed does not
depend on the model, so the pair of a case runs with identical seeds (common random numbers, V-C2).
It writes ``ionmc-vc-case<N>-<model>-k<K>.npz`` (per-batch E(z) [MeV] and the grid edges) and one JSON
step document (``#JSON-BEGIN`` / ``#JSON-END``): the effective-configuration summary including the
``stopping`` identity block, counters, energy balance, seeds, histories, the environment digests of
``env_digest.py`` and ``warp_modules``.

``vc2d-sample`` is the deterministic V-C2-D run: MCS and straggling off, 2000 histories, seed of
``r_index`` 15 (``20612004 + 1500 + k``). It records per history the end depth ``z_end`` and the end
energy ``E_end`` and the telescoped end point ``Z_m = mean(z_end + R_m(E_end))`` with ``R_m`` from the
run's own transport tables, and ``R_m(E0)``.

Warp backends start with a fresh, empty kernel cache: ``WARP_CACHE_PATH`` is set to a new directory
inside the output directory BEFORE ``warp`` is imported (nothing of ``ionmc`` is imported before
that), a non-empty directory refuses the start, and at the end the sha256 of every file of the cache
is recorded as ``warp_modules`` (plan section D step 3).

``predict`` computes model-only quantities (no TOPAS data) and writes
``validation/results/stopping/v3-013/vc-predictions.json``: ``Delta R_pred`` and ``Delta R_pred,D``,
``R_CSDA`` of both models at the three energies, both informativeness preconditions, the resulting
gating map, the identity pins of both models (spec / content / table digests), the ``src/ionmc`` tree
hash and the sha256 of this file.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FORMAT = 1
SUITE = "lv13"
PREDICTIONS_SCHEMA = "ionmc-v13-vc-predictions-1"
DEFAULT_PREDICTIONS = REPO / "validation" / "results" / "stopping" / "v3-013" / "vc-predictions.json"
SELF_REL = "validation/scripts/stopping/steps_v13.py"

# -- seeds (plan section C, family ion13) --
ION13_BASE = 20612004  # seed = base + 100 r_index + k
K_MAX = 99
R_INDEX_VC2D = 15

# -- run definitions (plan section A.6) --
N_HISTORIES = 1_000_000
N_BATCHES = 40
N_HISTORIES_D = 2000
N_BATCHES_D = 20
FIELD_MM = 120.0
SOURCE_Z_MM = -0.001  # as the TOPAS bundles (zero-size pencil on the axis from z = -0.001 mm)
MEMORY_BUDGET_BYTES = 2**32
MODEL_DEFAULT = "bethe-ionmc-v1"
MODEL_ICRU90 = "icru90-water-v1"
MODELS = (MODEL_ICRU90, MODEL_DEFAULT)
A_NUCLEONS = 4  # He-4
# case -> (MeV/u, box length L mm, IDD bin h_z mm, R_plan mm of the default model, rounded to 0.1 mm);
# the sizing table of plan section A.6 and the bundles of gen_v13_bundles.py
CASES: dict[int, tuple[float, float, float, float]] = {
    1: (75.0, 100.0, 0.05, 46.7),
    2: (125.0, 200.0, 0.05, 115.7),
    3: (200.0, 350.0, 0.1, 261.9),
}

# -- V-C2 / V-C2-D literals (plan section A.6) --
TOL_FLOOR_MM = 0.02  # V-C2 PASS: |dR80 - dR_pred| + 1.645 sigma <= max(0.02 mm, 0.10 |dR_pred|)
TOL_FRACTION = 0.10
Z_TOST = 1.645
SE_R80_1E5_MM = 0.01  # R80 standard error of the V3-008 TD-EM He budget (1e5 histories)
N_BUDGET = 1e5
SIGMA_EXP_MM = math.sqrt(2.0) * SE_R80_1E5_MM * math.sqrt(N_BUDGET / N_HISTORIES)
VC2_PRECONDITION_MM = 2.0 * (TOL_FLOOR_MM + Z_TOST * SIGMA_EXP_MM)
VC2D_INFORMATIVE_FRACTION = 5e-4  # |dR_pred,D| >= 5e-4 R_CSDA,default(E0)
VC2D_ROUND_TRIP = 1e-4  # PASS (i) and (ii): <= 1e-4 R

HAVE_ENV_DIGEST_NOTE = "env_digest.py not present at this commit"


# -- pure definitions (no ionmc, no warp) --
def r_index_of(case: int) -> int:
    """The ion13 r_index of the V-C1/V-C2 sample of ``case`` (0-2)."""
    if case not in CASES:
        raise ValueError(f"unknown case {case}; cases are 1..3")
    return case - 1


def seed_of(case: int, k: int = 0) -> int:
    """``20612004 + 100 r_index + k``; identical for both models of a case."""
    if not 0 <= k <= K_MAX:
        raise ValueError(f"chunk index k must be in 0..{K_MAX}")
    return ION13_BASE + 100 * r_index_of(case) + k


def d_seed_of(case: int, k: int = 0) -> int:
    """The V-C2-D seed ``20612004 + 100 * 15 + k`` (r_index 15; the case enters through the energy)."""
    if case not in CASES:
        raise ValueError(f"unknown case {case}; cases are 1..3")
    if not 0 <= k <= K_MAX:
        raise ValueError(f"chunk index k must be in 0..{K_MAX}")
    return ION13_BASE + 100 * R_INDEX_VC2D + k


def grid_spec(case: int) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[int, int, int]]:
    """``(origin, spacing, shape)`` [mm] of the IDD grid of ``case`` on the TOPAS bin edges."""
    _eu, length, hz, _r = CASES[case]
    half = FIELD_MM / 2.0
    return (-half, -half, 0.0), (FIELD_MM, FIELD_MM, hz), (1, 1, round(length / hz))


def idd_edges_z_mm(case: int) -> np.ndarray:
    """Depth edges of the IDD grid [mm] (the same arithmetic as ``ScoringGrid.edges_mm``)."""
    origin, spacing, shape = grid_spec(case)
    return origin[2] + spacing[2] * np.arange(shape[2] + 1, dtype=np.float64)


def vc2_tolerance_mm(delta_r_pred_mm: float) -> float:
    """``max(0.02 mm, 0.10 |dR_pred|)`` (V-C2 tolerance)."""
    return max(TOL_FLOOR_MM, TOL_FRACTION * abs(delta_r_pred_mm))


def vc2_informative(delta_r_pred_mm: float) -> bool:
    """``|dR_pred| >= 2 (0.02 mm + 1.645 sigma_exp)``, sigma_exp = sqrt 2 * 0.01 mm * sqrt(1e5 / 1e6)."""
    return abs(delta_r_pred_mm) >= VC2_PRECONDITION_MM


def vc2d_informative(delta_r_pred_d_mm: float, r_csda_default_mm: float) -> bool:
    """``|dR_pred,D| >= 5e-4 R_CSDA,default(E0)``."""
    return abs(delta_r_pred_d_mm) >= VC2D_INFORMATIVE_FRACTION * r_csda_default_mm


def gaussian_average(f: Any, e0: float, sigma: float, n: int = 33) -> float:
    """``E[f(E)]`` for ``E ~ N(e0, sigma^2)`` (Gauss-Hermite; ``f(e0)`` for ``sigma = 0``)."""
    if sigma <= 0.0:
        return float(f(e0))
    x, w = np.polynomial.hermite.hermgauss(n)
    return float(sum(wi * f(e0 + math.sqrt(2.0) * sigma * xi) for xi, wi in zip(x, w)) / math.sqrt(math.pi))


# -- output helpers --
def _json_default(o: Any) -> Any:
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def emit(doc: dict[str, Any]) -> int:
    """Print the step document and return the exit status (0 iff ``pass``)."""
    doc = {
        **doc,
        "suite": os.environ.get("IONMC_RUN_SUITE", SUITE),
        "git_sha": os.environ.get("IONMC_RUN_SHA", "unknown"),
    }
    print("#JSON-BEGIN")
    print(json.dumps(doc, indent=1, sort_keys=True, default=_json_default))
    print("#JSON-END")
    return 0 if doc.get("pass") else 1


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# -- environment digests and the fresh Warp kernel cache --
def environment_digests() -> tuple[dict[str, Any] | None, str | None]:
    """The digests of ``env_digest.py`` (lazy; ``(None, note)`` if the module does not exist yet)."""
    path = HERE / "env_digest.py"
    if not path.is_file():
        return None, HAVE_ENV_DIGEST_NOTE
    spec = importlib.util.spec_from_file_location("_v13_env_digest", path)
    if spec is None or spec.loader is None:
        return None, "env_digest.py could not be loaded"
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for name in ("environment_digests", "compute_digests", "compute_environment_digests", "compute"):
        fn = getattr(mod, name, None)
        if callable(fn):
            return fn(), None
    return None, "env_digest.py provides none of the expected entry points"


def prepare_warp_cache(out_dir: Path, tag: str) -> Path:
    """A new empty kernel-cache directory under ``out_dir``, set as ``WARP_CACHE_PATH`` before
    ``warp`` is imported; refuses to start if it is non-empty or ``warp`` is already imported."""
    d = out_dir / f"warp-cache-{tag}"
    if d.exists() and any(d.iterdir()):
        raise SystemExit(f"refusing to start: the Warp kernel cache directory {d} is not empty")
    if "warp" in sys.modules:
        raise SystemExit("warp was imported before the fresh kernel cache was set (run-validity failure)")
    d.mkdir(parents=True, exist_ok=True)
    os.environ["WARP_CACHE_PATH"] = str(d.resolve())
    return d


def warp_modules(cache_dir: Path) -> dict[str, str]:
    """sha256 of every file of the kernel cache (relative path -> digest)."""
    return {
        p.relative_to(cache_dir).as_posix(): _sha256_file(p)
        for p in sorted(cache_dir.rglob("*"))
        if p.is_file()
    }


# -- configuration (ionmc imported lazily) --
def build_config(
    case: int, model: str, *, n_histories: int, n_batches: int, backend: str, seed: int,
    deterministic: bool = False, track_end: bool = False,
) -> Any:  # fmt: skip
    """The ionmc configuration of ``case`` and ``model`` (He, water I = 78 eV, nuclear off, default
    end-of-range model). ``deterministic`` switches MCS and straggling off (V-C2-D)."""
    from ionmc.config import DiagnosticsOptions, PhysicsOptions, RunOptions, SimulationConfig
    from ionmc.geometry import BoxPhantom
    from ionmc.materials import WATER
    from ionmc.physics.projectiles import ALPHA
    from ionmc.scoring import ScoringGrid
    from ionmc.sources import PencilBeamSource

    if model not in MODELS:
        raise ValueError(f"model must be one of {MODELS}")
    e_u, length, _hz, _r = CASES[case]
    origin, spacing, shape = grid_spec(case)
    half = FIELD_MM / 2.0
    off: dict[str, Any] = {"multiple_scattering": False, "straggling": False} if deterministic else {}
    return SimulationConfig(
        source=PencilBeamSource(ALPHA, (0.0, 0.0, SOURCE_Z_MM), (0.0, 0.0, 1.0), e_u * ALPHA.a),
        geometry=BoxPhantom((-half, -half, 0.0), (FIELD_MM, FIELD_MM, length), WATER),
        scoring=(ScoringGrid(origin, spacing, shape, name="idd"),),
        physics=PhysicsOptions(nuclear=False, stopping_model=model, **off),
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision="float64",
            seed=seed,
            n_histories=n_histories,
            n_batches=n_batches,
            memory_budget_bytes=MEMORY_BUDGET_BYTES,
        ),
        diagnostics=DiagnosticsOptions(track_end_positions=track_end),
    )


def _balance(res: Any) -> dict[str, Any]:
    b = res.energy_balance
    d = {f.name: getattr(b, f.name) for f in dataclasses.fields(b)}
    d["closure_residual_mev"] = b.closure_residual_mev
    d["relative_residual"] = b.relative_residual
    return d


def _run_header(case: int, model: str, k: int, seed: int, r_index: int) -> dict[str, Any]:
    digests, note = environment_digests()
    return {
        "format": FORMAT,
        "case": case,
        "ion": "he",
        "e_mev_per_u": CASES[case][0],
        "model": model,
        "r_index": r_index,
        "k": k,
        "seed": seed,
        "seed_family": {"family": "ion13", "base": ION13_BASE, "r_index": r_index, "k": k},
        "environment_digests": digests,
        "environment_digests_note": note,
    }


# -- steps --
def step_vc_sample(a: argparse.Namespace) -> int:
    """Run one V-C1/V-C2 sample (a sample-seal artefact; no V-C statistic is computed here)."""
    n = max(N_BATCHES, int(round(N_HISTORIES * a.scale / N_BATCHES)) * N_BATCHES)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache_dir = prepare_warp_cache(out, f"vc-case{a.case}-{a.model}-k{a.k}")  # before any warp import
    header = _run_header(a.case, a.model, a.k, seed_of(a.case, a.k), r_index_of(a.case))
    from ionmc.simulation import Simulation

    cfg = build_config(a.case, a.model, n_histories=n, n_batches=N_BATCHES, backend="warp-cpu", seed=header["seed"])
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    grid = res.grid("idd")
    arrays: dict[str, Any] = {"idd_batch_energy_mev": np.asarray(grid.batch_energy_mev)}
    for ax in range(3):
        arrays[f"idd_edges{ax}_mm"] = cfg.scoring[0].edges_mm(ax)
    npz = out / f"ionmc-vc-case{a.case}-{a.model}-k{a.k}.npz"
    np.savez_compressed(npz, **arrays)
    counters = res.counters.as_dict()
    summary = res.effective_config.summary()
    doc = {
        **header,
        "row": "vc-sample",
        "histories": n,
        "frozen_histories": N_HISTORIES,
        "reduced": n != N_HISTORIES,
        "n_batches": N_BATCHES,
        "backend": cfg.run.backend,
        "precision": cfg.run.precision,
        "valid": bool(res.valid),
        "counters": counters,
        "energy_balance": _balance(res),
        "npz": npz.name,
        "npz_sha256": _sha256_file(npz),
        "wall_s": wall,
        "effective_config": summary,
        "stopping_identity": summary["stopping"],
        "warp_modules": warp_modules(cache_dir),
        "pass": bool(res.valid and not any(counters.values())),
    }
    return emit(doc)


def step_vc2d_sample(a: argparse.Namespace) -> int:
    """Run the deterministic V-C2-D sample (MCS and straggling off) and telescope the end points."""
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cache_dir = prepare_warp_cache(out, f"vc2d-case{a.case}-{a.model}-{a.backend}-k{a.k}") if a.backend != "python" else None
    header = _run_header(a.case, a.model, a.k, d_seed_of(a.case, a.k), R_INDEX_VC2D)
    from ionmc.materials import WATER
    from ionmc.simulation import Simulation

    cfg = build_config(
        a.case, a.model, n_histories=N_HISTORIES_D, n_batches=N_BATCHES_D, backend=a.backend,
        seed=header["seed"], deterministic=True, track_end=True,
    )  # fmt: skip
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    z_end = np.asarray(res.diagnostics["end_position_mm"])[:, 2]
    e_end = np.asarray(res.diagnostics["end_energy_mev"])
    codes = np.asarray(res.diagnostics["end_code"])
    tables = res.effective_config.tables  # the run's own tables
    to_mm = 10.0 / WATER.density_g_cm3
    r_end = np.array([tables.range_g_cm2(0, float(e)) * to_mm for e in e_end])
    e0 = float(cfg.source.kinetic_energy_mev)
    r_e0 = float(tables.range_g_cm2(0, e0)) * to_mm
    z_m = float(np.mean(z_end + r_end))
    e_cut = float(res.effective_config.e_cut_mev)
    out.mkdir(parents=True, exist_ok=True)
    npz = out / f"ionmc-vc2d-case{a.case}-{a.model}-{a.backend}-k{a.k}.npz"
    np.savez_compressed(npz, z_end_mm=z_end, e_end_mev=e_end, r_end_mm=r_end, end_code=codes)
    counters = res.counters.as_dict()
    clean = bool(res.valid and not any(counters.values()))
    all_cutoff = bool(np.all(codes == 0))
    below_cut = bool(np.all(e_end <= e_cut * (1.0 + 1e-12)))
    summary = res.effective_config.summary()
    doc = {
        **header,
        "row": "vc2d-sample",
        "histories": N_HISTORIES_D,
        "n_batches": N_BATCHES_D,
        "backend": cfg.run.backend,
        "precision": cfg.run.precision,
        "valid": bool(res.valid),
        "counters": counters,
        "energy_balance": _balance(res),
        "e0_mev": e0,
        "e_cut_mev": e_cut,
        "r_e0_mm": r_e0,
        "z_m_mm": z_m,
        "z_m_minus_r_e0_relative": (z_m - r_e0) / r_e0,
        "z_end_mm": {"min": float(z_end.min()), "max": float(z_end.max()), "mean": float(z_end.mean())},
        "e_end_mev": {"min": float(e_end.min()), "max": float(e_end.max())},
        "end_codes_all_cutoff": all_cutoff,
        "e_end_all_at_or_below_cut": below_cut,
        "npz": npz.name,
        "npz_sha256": _sha256_file(npz),
        "wall_s": wall,
        "effective_config": summary,
        "stopping_identity": summary["stopping"],
        "warp_modules": None if cache_dir is None else warp_modules(cache_dir),
        "pass": bool(clean and all_cutoff and below_cut),
    }
    return emit(doc)


# -- predictions (model-only; no TOPAS data) --
def _git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True, text=True).stdout.strip()


def compute_predictions(cache_dir: str | None = None) -> dict[str, Any]:
    """The content of ``vc-predictions.json`` (no file is written). Needs the ICRU 90 cache object."""
    from ionmc.config import validate
    from ionmc.materials import WATER
    from ionmc.physics.projectiles import ALPHA
    from ionmc.physics.stopping import BetheStoppingSource
    from ionmc.transport.tables import TransportTables

    if cache_dir is not None:
        os.environ["IONMC_CACHE_DIR"] = cache_dir
    to_mm = 10.0 / WATER.density_g_cm3
    plan_tables = TransportTables.from_stopping_tables([BetheStoppingSource().table(WATER, ALPHA)], projectile=ALPHA)
    eff = {
        (m, c): validate(build_config(c, m, n_histories=N_HISTORIES, n_batches=N_BATCHES, backend="warp-cpu", seed=seed_of(c)))
        for m in MODELS
        for c in CASES
    }
    energies = []
    for c, (e_u, length, hz, r_plan_table) in CASES.items():
        e0 = e_u * ALPHA.a
        r_plan = round(float(plan_tables.range_g_cm2(0, e0)) * to_mm, 1)
        if r_plan != r_plan_table:
            raise SystemExit(f"case {c}: R_plan {r_plan} mm differs from the plan sizing table {r_plan_table} mm")
        t_def, t_icru = eff[MODEL_DEFAULT, c].tables, eff[MODEL_ICRU90, c].tables
        r_def = float(t_def.range_g_cm2(0, e0)) * to_mm
        r_icru = float(t_icru.range_g_cm2(0, e0)) * to_mm
        sigma_e = 0.0  # the TD-EM source is monoenergetic (BeamEnergySpread 0)
        d_d = r_icru - r_def
        d = gaussian_average(
            lambda e, ti=t_icru, td=t_def: (float(ti.range_g_cm2(0, e)) - float(td.range_g_cm2(0, e))) * to_mm,
            e0,
            sigma_e,
        )
        vc2_inf, vc2d_inf = vc2_informative(d), vc2d_informative(d_d, r_def)
        energies.append(
            {
                "case": c,
                "e_mev_per_u": e_u,
                "e0_total_mev": e0,
                "sigma_e_mev": sigma_e,
                "r_plan_mm": r_plan,
                "box_length_mm": length,
                "idd_bin_mm": hz,
                "r_csda_default_mm": r_def,
                "r_csda_icru90_mm": r_icru,
                "delta_r_pred_d_mm": d_d,
                "delta_r_pred_mm": d,
                "vc2_tolerance_mm": vc2_tolerance_mm(d),
                "vc2_precondition_holds": vc2_inf,
                "vc2_gate": "gating" if vc2_inf else "non-informative (gate for this energy is V-C2-D alone)",
                "vc2d_threshold_mm": VC2D_INFORMATIVE_FRACTION * r_def,
                "vc2d_informative": vc2d_inf,
                "vc2d_gate": "gating" if vc2d_inf else "substitute: CC-T (both backends) + CC-ID + PER-1 for icru90-water-v1",
                "vc2d_round_trip_tolerance_mm": {m: VC2D_ROUND_TRIP * r for m, r in ((MODEL_DEFAULT, r_def), (MODEL_ICRU90, r_icru))},
                "vc2d_difference_tolerance_mm": VC2D_ROUND_TRIP * r_def,
            }
        )
    pins: dict[str, Any] = {}
    blocks: dict[str, Any] = {}
    for m in MODELS:
        per_case = {c: eff[m, c].summary() for c in CASES}
        first = per_case[1]
        pins[m] = {
            "model_id_effective": first["stopping"]["model_id_effective"],
            "version": first["stopping"]["version"],
            "spec_sha256": first["stopping"]["spec_sha256"],
            "content_sha256": first["stopping"]["content_sha256"],
            "tables_sha256": {str(c): s["tables"]["sha256"] for c, s in per_case.items()},
            "tables_materials": {str(c): s["tables"]["materials"] for c, s in per_case.items()},
        }
        blocks[m] = {str(c): s["stopping"] for c, s in per_case.items()}
    tree = _git("rev-parse", "HEAD:src/ionmc")
    dirty = _git("status", "--porcelain", "--", "src/ionmc")
    if dirty:
        raise SystemExit(f"src/ionmc has uncommitted changes; the tree hash would not describe the run:\n{dirty}")
    return {
        "schema": PREDICTIONS_SCHEMA,
        "plan": "validation/plans/v3-013-acceptance.md section A.6 (V-C2, V-C2-D)",
        "note": "model-only quantities (ranges and differences in mm of water); no TOPAS data, no stopping array",
        "constants": {
            "se_r80_1e5_mm": SE_R80_1E5_MM,
            "n_budget": N_BUDGET,
            "n_histories": N_HISTORIES,
            "sigma_exp_mm": SIGMA_EXP_MM,
            "vc2_tolerance_floor_mm": TOL_FLOOR_MM,
            "vc2_tolerance_fraction": TOL_FRACTION,
            "tost_z": Z_TOST,
            "vc2_precondition_mm": VC2_PRECONDITION_MM,
            "vc2d_informative_fraction": VC2D_INFORMATIVE_FRACTION,
            "vc2d_round_trip": VC2D_ROUND_TRIP,
        },
        "energies": energies,
        "models": pins,
        "identity_blocks": blocks,
        "src_ionmc_tree": tree,
        "script": {"path": SELF_REL, "sha256": _sha256_file(Path(__file__).resolve())},
    }


def step_predict(a: argparse.Namespace) -> int:
    doc = compute_predictions()
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out} sha256 {_sha256_file(out)}")
    for e in doc["energies"]:
        print(
            f"case {e['case']} He {e['e_mev_per_u']:g} MeV/u: R_def {e['r_csda_default_mm']:.4f} mm, "
            f"dR_pred,D {e['delta_r_pred_d_mm']:+.5f} mm, V-C2 {e['vc2_gate']}, V-C2-D {e['vc2d_gate']}"
        )
    return 0


STEPS = {"vc-sample": step_vc_sample, "vc2d-sample": step_vc2d_sample, "predict": step_predict}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--case", type=int, default=None, help="1..3 = He 75 / 125 / 200 MeV/u")
    ap.add_argument("--model", choices=MODELS, default=None)
    ap.add_argument("--backend", choices=("python", "warp-cpu"), default=None, help="vc2d-sample only")
    ap.add_argument("--k", type=int, default=0, help="chunk 0..99")
    ap.add_argument("--scale", type=float, default=1.0, help="vc-sample: history-count factor (reduced run)")
    ap.add_argument("--out-dir", default="samples")
    ap.add_argument("--out", default=str(DEFAULT_PREDICTIONS), help="predict: output file")
    ap.add_argument("--timeout", type=float, default=None, help="accepted for run_suite uniformity")
    a = ap.parse_args(argv)
    if a.step != "predict":
        if a.case not in CASES:
            raise SystemExit("--case must be 1, 2 or 3")
        if a.model is None:
            raise SystemExit("--model is required")
        seed_of(a.case, a.k)  # range check
    if a.step == "vc2d-sample" and a.backend is None:
        raise SystemExit("vc2d-sample needs --backend")
    if not 0.0 < a.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    return STEPS[a.step](a)


if __name__ == "__main__":
    sys.exit(main())
