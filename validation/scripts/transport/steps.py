"""Validation steps of the V3-003 transport acceptance plan (T1, T2, T8-T10, T12-T14).

Usage::

    python validation/scripts/transport/steps.py <step> [--workers N] [--n HISTORIES] ...

Each step prints one JSON document (aggregates only) between the lines ``#JSON-BEGIN`` and
``#JSON-END`` (Warp writes its own messages to stdout) and exits 0 iff every criterion of the
step passes. ``run_suite.py`` archives the output verbatim. Criteria are those frozen in
``validation/plans/v3-003-acceptance.md``; the history counts are the frozen ones unless ``--n``
reduces them, in which case the document says ``"reduced": true`` and a pass is not a
conformant result. Everything uses the offline analytic Bethe stopping source (I = 78 eV), the
production float32 Warp CPU backend unless stated, and the helpers of
``ionmc.transport.parity`` and ``ionmc.transport.mcs_checks`` that the tests use.

Seeds: every statistical step derives all its seeds deterministically from one
``--seed-base`` (default 20331004, the qualification base; 20261004 was the rehearsal,
20271004 the T9 investigation, 20281004, 20291004, 20301004, 20311004 and 20321004 the first to
fifth consumed qualification attempts):
  T12 samples ``base + 1000 k`` (k = 1 python, 2 cpu32, 3 cpu64, 4 cuda32);
  T9 ``base + i`` (i-th case); T8 ``base + i`` (theta) and ``base + 10 + i`` (sigma);
  T14 ``base + i`` (theta), ``base + 100 + i`` (sigma), ``base + 1000`` (control);
  T10 ``base + i`` (i-th orientation); T1, T13 and the repeatability step ``base``.
The base is recorded in the step documents and in the metadata of every T12 sample, and the
comparison refuses samples made with another base.

Observables (see ``docs/architecture/transport.md``, section Validation runner): T8 and T14 take
the lateral sigma from the deposited energy in fixed 0.2 mm lateral bins and 1 mm slabs at
z/R = 0.5 and 0.9 (scoring grids independent of the transport voxels; Fermi-Eyges A2 comparison)
and the exit ``theta_rms`` of a 0.5 R1 slab from the escape records (T14 ends the world at the
slab thickness with ``VoxelGeometry.z_exit_mm``); T9 uses 1 mm depth-dose bins; T10 projects
1 mm path grids onto the beam axis (``projected_idd``, 0.1 mm bins) for R80 and uses a
whole-world grid for the total energy. Scoring steps may be longer than the scoring voxels: the
track-length scoring walks each leg over the voxels it crosses up to a piece bound computed at
validation (there is no step-length rule). The T14 control is the "truncate-first" diagnostic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
)
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Result, Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport import mcs_checks as mc
from ionmc.transport import parity
from ionmc.transport.tables import TransportTables

STOP = BetheStoppingSource()
TABLES = TransportTables.from_stopping_tables([STOP.table(WATER, PROTON)])


def r_csda_mm(e_mev: float) -> float:
    """CSDA range of the project table in water [mm]."""
    return TABLES.range_g_cm2(0, e_mev) * 10.0 / WATER.density_g_cm3


STEP_MEMORY_BUDGET = 2**34
"""Memory budget [bytes] of the deposit accumulators that the validation steps declare (16 GiB:
host RAM headroom for the per-worker private copies); recorded in the step documents. A step
never requests more worker processes than this budget allows (``plan_workers``)."""
RUN_PLANS: list[dict[str, Any]] = []


def plan_workers(
    requested: int,
    backend: str,
    n_batches: int,
    scoring: tuple[ScoringGrid, ...],
    n_histories: int,
    budget: int = STEP_MEMORY_BUDGET,
) -> dict[str, Any]:
    """Worker processes to use: CUDA never uses CPU workers; CPU runs use at most ``requested``
    workers, at most one per history and at most ``budget // per-worker accumulator bytes``
    (``n_batches * voxels * 8`` bytes per worker), never fewer than one."""
    per_worker = n_batches * sum(g.n_voxels for g in scoring) * 8
    if backend == "warp-cuda":
        used = 1
    else:
        used = max(1, min(requested, n_histories, budget // max(per_worker, 1)))
    return {
        "backend": backend,
        "workers_requested": requested,
        "workers_used": used,
        "accumulator_bytes_per_worker": per_worker,
        "memory_budget_bytes": budget,
    }


def run_cfg(
    *,
    energy: float,
    geometry: Any,
    scoring: tuple[ScoringGrid, ...],
    backend: str = "warp-cpu",
    precision: str = "float32",
    seed: int = 1,
    n: int,
    n_batches: int = 20,
    workers: int = 1,
    mcs: bool = True,
    straggling: bool = True,
    max_step: float = 1.0,
    frac: float = 0.02,
    position: tuple[float, float, float] = (0.0, 0.0, 0.0),
    direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
    diag: DiagnosticsOptions | None = None,
    timeout: float | None = None,
    chunk: int | None = None,
    trunc_diag: bool = False,
    memory_budget: int | None = None,
    straggling_model: str = "bohr_gamma_v1",
) -> Result:
    kw = {"chunk_histories": chunk} if chunk else {}
    plan = plan_workers(
        workers, backend, n_batches, scoring, n, memory_budget or STEP_MEMORY_BUDGET
    )
    RUN_PLANS.append(plan)
    workers = plan["workers_used"]
    kw["memory_budget_bytes"] = plan["memory_budget_bytes"]
    cfg = SimulationConfig(
        source=PencilBeamSource(PROTON, position, direction, energy),
        geometry=geometry,
        scoring=scoring,
        physics=PhysicsOptions(
            nuclear=False,
            stopping=STOP,
            straggling=straggling,
            multiple_scattering=mcs,
            max_step_mm=max_step,
            max_energy_loss_fraction=frac,
            truncated_hinge_diagnostic=trunc_diag,
            straggling_model=straggling_model,
        ),
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision=precision,  # type: ignore[arg-type]
            seed=seed,
            n_histories=n,
            n_batches=n_batches,
            cpu_workers=workers,
            worker_timeout_s=timeout,
            **kw,
        ),
        diagnostics=diag or DiagnosticsOptions(),
    )
    return Simulation(cfg).run()


def finish(doc: dict[str, Any], frozen_n: int | None, n: int | None) -> int:
    doc["reduced"] = bool(frozen_n is not None and n is not None and n < frozen_n)
    doc["frozen_histories"] = frozen_n
    doc["histories"] = n
    return emit(doc)


# the qualification base; 20261004 rehearsal, 20271004 T9 investigation, 20281004, 20291004,
# 20301004, 20311004 and 20321004 consumed by the first to fifth qualification attempts (all
# non-qualification)
DEFAULT_SEED_BASE = 20331004
SEED_BASE = DEFAULT_SEED_BASE


def emit(doc: dict[str, Any]) -> int:
    """Print the result document (with the run identity) and return the exit status."""
    if RUN_PLANS:
        doc = {
            **doc,
            "worker_plan": {
                "runs": len(RUN_PLANS),
                "workers_requested": max(p["workers_requested"] for p in RUN_PLANS),
                "workers_used_min": min(p["workers_used"] for p in RUN_PLANS),
                "workers_used_max": max(p["workers_used"] for p in RUN_PLANS),
                "max_accumulator_bytes_per_worker": max(
                    p["accumulator_bytes_per_worker"] for p in RUN_PLANS
                ),
                "memory_budget_bytes": STEP_MEMORY_BUDGET,
            },
        }
    doc = {
        **doc,
        "seed_base": SEED_BASE,
        "suite": os.environ.get("IONMC_RUN_SUITE", "unknown"),
        "git_sha": os.environ.get("IONMC_RUN_SHA", "unknown"),
    }
    print("#JSON-BEGIN")
    print(json.dumps(doc, indent=1, sort_keys=True, default=_json))
    print("#JSON-END")
    return 0 if doc.get("pass") else 1


def _json(o: Any) -> Any:
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


# -- T1 ------------------------------------------------------------------------------------------
def step_t1(a: argparse.Namespace) -> int:
    k, energy = a.k, a.energy
    shape = (12, 12, 40)
    geo = VoxelGeometry(
        (-30.0, -30.0, 0.0), (5.0, 5.0, 5.0), shape, (WATER,), np.zeros(shape, dtype=np.int32)
    )
    grid = (ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 100), name="dose"),)
    diag = DiagnosticsOptions(track_end_positions=True, trace_histories=k)
    common = dict(
        energy=energy,
        geometry=geo,
        scoring=grid,
        seed=a.seed,
        n=k,
        n_batches=2,
        max_step=2.0,
        diag=diag,
        precision="float64",
    )
    t0 = time.perf_counter()
    ref = run_cfg(backend="python", workers=a.workers, timeout=a.timeout, **common)
    t_py = time.perf_counter() - t0
    t0 = time.perf_counter()
    wrp = run_cfg(backend="warp-cpu", **common)
    t_wp = time.perf_counter() - t0
    v = parity.compare_traces(ref.diagnostics, wrp.diagnostics)
    v["python_seconds"], v["warp_cpu_seconds"] = t_py, t_wp
    v["counters_python"] = ref.counters.as_dict()
    v["counters_warp"] = wrp.counters.as_dict()
    v["pass"] = bool(v["pass"] and not ref.counters.any_nonzero and not wrp.counters.any_nonzero)
    doc = {
        "step": "t1",
        "energy_mev": energy,
        "k": k,
        "workers": a.workers,
        "t1": v,
        "pass": v["pass"],
    }
    return finish(doc, 256, k)


# -- T2 ------------------------------------------------------------------------------------------
def step_t2(a: argparse.Namespace) -> int:
    backend = a.backend
    out = {}
    ok = True
    for e in (100.0, 150.0, 200.0):
        depth = 1.1 * r_csda_mm(e)
        res = run_cfg(
            energy=e,
            geometry=BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, depth), WATER),
            scoring=(ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, int(depth / 2))),),
            backend=backend,
            n=4,
            n_batches=2,
            mcs=False,
            straggling=False,
            max_step=2.0,
            diag=DiagnosticsOptions(track_end_positions=True),
        )
        r0 = r_csda_mm(e)
        r_cut = TABLES.range_g_cm2(0, 2.0) * 10.0 / WATER.density_g_cm3
        z = res.diagnostics["end_position_mm"][:, 2]
        lo, hi = r0 - r_cut - 1e-4 * r0, r0 + 1e-4 * r0
        passed = bool(np.all((z >= lo) & (z <= hi)) and res.valid)
        out[str(e)] = {"z_end": z.tolist(), "lower": lo, "upper": hi, "pass": passed}
        ok &= passed
    return finish({"step": "t2", "backend": backend, "t2": out, "pass": ok}, None, None)


# -- repeatability (python-scope Warp functions and the kernels must be deterministic) -----------
def step_repeat(a: argparse.Namespace) -> int:
    """Run each ``backend:precision:histories`` of ``--runs`` twice (same configuration and seed,
    150 MeV, all physics on) and require bit-identical tallies, counters, deposit grids and
    per-history end states."""
    depth = 1.1 * r_csda_mm(150.0)
    geo = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, depth), WATER)
    grid = (ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, int(depth / 2))),)
    out: dict[str, Any] = {}
    ok = True
    n_min = None
    for spec in a.runs.split(","):
        backend, prec, n_s = spec.split(":")
        n = int(n_s)
        n_min = n if n_min is None else min(n_min, n)
        runs = []
        for _ in range(2):
            runs.append(
                run_cfg(
                    energy=150.0,
                    geometry=geo,
                    scoring=grid,
                    backend=backend,
                    precision=prec,
                    seed=a.seed,
                    n=n,
                    n_batches=20,
                    workers=a.workers if backend in ("python", "warp-cpu") else 1,
                    timeout=a.timeout,
                    max_step=2.0,
                    diag=DiagnosticsOptions(track_end_positions=True),
                )
            )
        v = parity.compare_runs_bitwise(runs[0], runs[1])
        v["valid"] = bool(runs[0].valid and runs[1].valid)
        v["n_histories"] = n
        out[spec] = v
        ok &= bool(v["identical"] and v["valid"])
    frozen = 100_000 if "cuda" in a.runs else None
    return finish({"step": "t-r1", "t_r1": out, "pass": ok}, frozen, n_min)


# -- T13 -----------------------------------------------------------------------------------------
def step_t13(a: argparse.Namespace) -> int:
    e = a.energy
    depth = 1.1 * r_csda_mm(e)
    geo = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, depth), WATER)
    grid = (ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, int(depth / 2))),)
    out: dict[str, Any] = {}
    ok = True
    n = a.n
    base = dict(energy=e, geometry=geo, scoring=grid, seed=a.seed, n=n, n_batches=20)
    if a.mode == "workers":
        if a.workers < 2:
            raise SystemExit("t13 workers mode needs --workers >= 2")
        for prec in ("float32", "float64"):
            one = run_cfg(precision=prec, workers=1, **base)
            many = run_cfg(precision=prec, workers=a.workers, timeout=a.timeout, **base)
            v = parity.compare_partition(one, many)
            out[f"workers_1_vs_{a.workers}_{prec}"] = v
            ok &= v["pass"]
    else:  # chunk sizes on the selected backend (warp-cuda in the HR suite)
        small = run_cfg(backend=a.backend, precision="float32", chunk=2**10, **base)
        large = run_cfg(backend=a.backend, precision="float32", chunk=2**18, **base)
        v = parity.compare_partition(small, large)
        v["chunk_seconds_small"] = max(
            p["chunk_seconds"][0] for p in small.transport_report["partials"]
        )
        v["register_count"] = large.transport_report["partials"][0].get("register_count")
        out[f"chunks_1024_vs_262144_{a.backend}"] = v
        ok &= v["pass"]
    return finish({"step": "t13", "mode": a.mode, "t13": out, "pass": ok}, None, n)


# -- T12 (split into sample steps and a comparison step) -----------------------------------------
# name -> (backend, precision, frozen history count, batches, seed index)
SAMPLES = {
    "python": ("python", "float64", 4000, 40, 1),
    "cpu32": ("warp-cpu", "float32", 1_000_000, 100, 2),
    "cpu64": ("warp-cpu", "float64", 1_000_000, 100, 3),
    "cuda32": ("warp-cuda", "float32", 1_000_000, 100, 4),
}
SAMPLE_FORMAT = 1


def sample_histories(name: str, scale: float) -> int:
    """History count of a sample: its frozen count times ``scale`` (a multiple of its batches)."""
    _, _, n, b, _ = SAMPLES[name]
    return n if scale == 1.0 else max(b, int(n * scale) // b * b)


def sample_config(a: argparse.Namespace, name: str, workers: int) -> tuple[Any, parity.T12Layout]:
    """Configuration of sample ``name`` with the worker count planned from the memory budget
    (CUDA samples always use one CPU worker); the requested count is not changed by the plan."""
    backend, prec, _, b, k = SAMPLES[name]
    n = sample_histories(name, a.scale)
    kwargs = dict(
        energy_mev=a.energy,
        backend=backend,
        precision=prec,
        seed=a.seed + 1000 * k,  # a distinct seed per sample
        n_histories=n,
        n_batches=b,
        lateral_bin_mm=a.lateral_bin,
        half_width_mm=a.half_width,
        timeout_s=a.timeout,
    )
    cfg, layout = parity.t12_config(workers=1, **kwargs)
    plan = plan_workers(workers, backend, b, cfg.scoring, n)
    return parity.t12_config(
        workers=plan["workers_used"], memory_budget_bytes=plan["memory_budget_bytes"], **kwargs
    )


def config_fingerprint(eff_summary: dict[str, Any]) -> str:
    """sha256 of the effective-configuration summary without the parallelism settings."""
    s = dict(eff_summary)
    s.pop("cpu_workers", None)
    s.pop("chunk_histories", None)
    return hashlib.sha256(json.dumps(s, sort_keys=True, default=str).encode()).hexdigest()


def expected_sample(a: argparse.Namespace, name: str) -> dict[str, Any]:
    """The metadata a sample of this comparison must carry, from the configuration built here
    (never from the sample file): the effective-configuration fingerprint (a CUDA sample is
    validated as if a device existed, so the comparison runs on any host), seed, history count,
    batches, backend, precision, energy, scoring geometry parameters and the frozen count."""
    import ionmc.config as ionmc_config
    from ionmc.config import validate

    cfg, _ = sample_config(a, name, 1)
    real = ionmc_config.cuda_available
    ionmc_config.cuda_available = lambda: True  # type: ignore[assignment]
    try:
        fingerprint = config_fingerprint(validate(cfg).summary())
    finally:
        ionmc_config.cuda_available = real  # type: ignore[assignment]
    return {
        "format": SAMPLE_FORMAT,
        "name": name,
        "backend": cfg.run.backend,
        "precision": cfg.run.precision,
        "n_total": cfg.run.n_histories,
        "n_batches": cfg.run.n_batches,
        "seed": cfg.run.seed,
        "seed_base": a.seed,
        "energy_mev": a.energy,
        "lateral_bin_mm": a.lateral_bin,
        "half_width_mm": a.half_width,
        "frozen_histories": SAMPLES[name][2],
        "config_fingerprint": fingerprint,
    }


def current_source_hashes() -> dict[str, str]:
    """sha256 of every execution-defining file of the tree this step runs from (recomputed, not
    read from any archive)."""
    import run_suite

    return {
        str(f.relative_to(run_suite.REPO)): run_suite.sha256(f) for f in run_suite.source_files()
    }


def verify_archives(dirs: list[Path], expected_sha: str) -> None:
    """Fail closed unless every archive that supplies samples is trustworthy: nothing stored is
    taken on trust.

    * this run's own archive (``dirs[0]``, still being written) must carry the expected SHA and
      the ``source_hashes`` that are recomputed from the tree now executing;
    * every other archive must verify completely (``summarize.verify``: manifest, step headers,
      trailers and result documents; any invalid step fails it), carry the expected SHA, and have
      ``source_hashes`` byte-equal to the recomputed ones - or, when git can read the repository,
      hashes that match the blobs of the expected commit (the attestation is recomputed here,
      a stored ``attestation`` is ignored).
    """
    import summarize

    def env_of(d: Path) -> dict[str, Any]:
        f = d / "environment.txt"
        if not f.exists():
            raise SystemExit(f"{d} has no environment.txt: sample source not verifiable")
        return summarize.parse_env(f.read_text())

    current = current_source_hashes()
    own = env_of(dirs[0])
    if own.get("git_sha") != expected_sha:
        raise SystemExit(f"{dirs[0]}: environment SHA {own.get('git_sha')} != {expected_sha}")
    if own["source_hashes"] != current:
        raise SystemExit(f"{dirs[0]}: recorded source hashes differ from the executing tree")
    for d in dirs[1:]:
        summary = summarize.verify(d, expected_sha)
        if not summary["pass"]:
            bad = {n: s["problems"] for n, s in summary["steps"].items() if not s["pass"]}
            raise SystemExit(f"{d}: archive does not verify: {summary['problems']} {bad}")
        env = env_of(d)
        if env.get("git_sha") != expected_sha:
            raise SystemExit(f"{d}: archive SHA {env.get('git_sha')} != {expected_sha}")
        if env["source_hashes"] == current:
            continue
        if not summarize.attest(env, expected_sha)["valid"]:
            raise SystemExit(
                f"{d}: source hashes differ from the executing tree and from the commit"
            )


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_sample_part(
    path: Path, name: str, h0: int, h1: int, a: argparse.Namespace, workers: int
) -> dict[str, Any]:
    """Transport the histories ``[h0, h1)`` of sample ``name`` and write the exact per-batch
    results (int64 deposit quanta, tally expansions, counters) to ``path`` (npz)."""
    from ionmc.config import validate
    from ionmc.transport.pool import run_pool
    from ionmc.transport.run import run_range
    from ionmc.transport.tally import concat_partials

    cfg, _ = sample_config(a, name, workers)
    workers_used = cfg.run.cpu_workers
    eff = validate(cfg)
    t0 = time.perf_counter()
    if workers_used > 1 and cfg.run.backend != "warp-cuda":
        parts = run_pool(eff, h_range=(h0, h1))
    else:
        parts = [run_range(eff, h0, h1)]
    seconds = time.perf_counter() - t0
    part = concat_partials(parts)
    flat = np.array([c for col in part.tally_components for c in col], dtype=np.float64)
    lens = np.array([len(col) for col in part.tally_components], dtype=np.int64)
    meta = {
        "format": SAMPLE_FORMAT,
        "name": name,
        "backend": cfg.run.backend,
        "precision": cfg.run.precision,
        "n_total": cfg.run.n_histories,
        "n_batches": cfg.run.n_batches,
        "seed": cfg.run.seed,
        "seed_base": a.seed,
        "h0": h0,
        "h1": h1,
        "energy_mev": a.energy,
        "lateral_bin_mm": a.lateral_bin,
        "half_width_mm": a.half_width,
        "frozen_histories": SAMPLES[name][2],
        "git_sha": os.environ.get("IONMC_RUN_SHA", "unknown"),
        "config_fingerprint": config_fingerprint(eff.summary()),
        "workers_requested": workers,
        "workers_used": workers_used,
        "memory_budget_bytes": cfg.run.memory_budget_bytes,
        "seconds": seconds,
        "device": parts[0].meta.get("device"),
    }
    arrays = {f"edep_{i}": e for i, e in enumerate(part.edep)}
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        meta=np.array(json.dumps(meta, sort_keys=True)),
        comps_flat=flat,
        comps_len=lens,
        counters=np.array(part.counter_sums, dtype=np.int64),
        **arrays,
    )
    return {"file": path.name, "sha256": file_sha256(path), "meta": meta}


def _sample_doc(entries: list[dict[str, Any]], a: argparse.Namespace) -> dict[str, Any]:
    ok = all(sum(e["counters"]) == 0 for e in entries)
    return {"samples": entries, "pass": ok, "format": SAMPLE_FORMAT}


def _write_parts(a: argparse.Namespace, names: list[str], part: tuple[int, int]) -> int:
    out_dir = Path(a.out_dir)
    entries = []
    for name in names:
        n_total = sample_histories(name, a.scale)
        i, nparts = part
        h0, h1 = (i - 1) * n_total // nparts, i * n_total // nparts
        fname = f"t12-{name}-part-{i}-of-{nparts}.npz"
        e = save_sample_part(out_dir / fname, name, h0, h1, a, a.workers)
        data = np.load(out_dir / fname)
        e["counters"] = [int(x) for x in data["counters"]]
        e["sample_file"] = f"samples/{fname}"
        entries.append(e)
    doc = _sample_doc(entries, a)
    doc["step"] = "t12-sample"
    reduced = any(e["meta"]["n_total"] < e["meta"]["frozen_histories"] for e in entries)
    doc["reduced"] = reduced
    return emit(doc)


def step_t12_python_sample(a: argparse.Namespace) -> int:
    i, n = (int(x) for x in a.part.split("/"))
    if not 1 <= i <= n:
        raise SystemExit("--part must be i/n with 1 <= i <= n")
    return _write_parts(a, ["python"], (i, n))


def step_t12_accelerated(a: argparse.Namespace) -> int:
    return _write_parts(a, a.samples.split(","), (1, 1))


def _producer_docs(dirs: list[Path]) -> dict[str, tuple[Path, dict[str, Any], str]]:
    """Map sample file name -> (archive step file, document entry, header sha) of every sample
    step archived in ``dirs``."""
    found: dict[str, tuple[Path, dict[str, Any], str]] = {}
    for d in dirs:
        for txt in sorted(d.glob("[0-9][0-9]-*.txt")):
            text = txt.read_text()
            if "#JSON-BEGIN" not in text:
                continue
            try:
                doc = json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])
            except json.JSONDecodeError:
                continue
            if doc.get("step") != "t12-sample":
                continue
            sha = next(
                (ln.split(": ", 1)[1] for ln in text.splitlines() if ln.startswith("# git_sha: ")),
                "",
            )
            for e in doc["samples"]:
                found[e["file"]] = (txt, e, sha)
    return found


def load_sample(
    name: str, dirs: list[Path], expected_sha: str, a: argparse.Namespace
) -> tuple[Any, dict[str, Any]]:
    """Load, hash-verify and merge all parts of sample ``name`` from ``dirs``: every part file
    must be recorded by a producer step archived with the expected SHA, with the recorded sha256;
    the parts must share one configuration and tile ``[0, n_total)`` exactly."""
    from ionmc.transport.tally import PartialTransport, merge_partials

    verify_archives(dirs, expected_sha)
    expected = expected_sample(a, name)
    producers = _producer_docs(dirs)
    parts = []
    metas = []
    for fname, (_txt, entry, sha) in sorted(producers.items()):
        if entry["meta"]["name"] != name:
            continue
        path = next((d / "samples" / fname for d in dirs if (d / "samples" / fname).exists()), None)
        if path is None:
            raise SystemExit(f"sample file {fname} recorded but not found in {dirs}")
        if file_sha256(path) != entry["sha256"]:
            raise SystemExit(f"sha256 mismatch for {fname}")
        if sha != expected_sha or entry["meta"]["git_sha"] != expected_sha:
            raise SystemExit(f"{fname} was produced at SHA {sha}/{entry['meta']['git_sha']}")
        data = np.load(path)
        meta = json.loads(str(data["meta"]))
        lens = data["comps_len"]
        flat = data["comps_flat"]
        comps, pos = [], 0
        for n in lens:
            comps.append([float(x) for x in flat[pos : pos + int(n)]])
            pos += int(n)
        edep = [
            data[k]
            for k in sorted(
                (k for k in data.files if k.startswith("edep_")), key=lambda k: int(k.split("_")[1])
            )
        ]
        parts.append(
            PartialTransport(
                meta["h0"], meta["h1"], comps, [int(x) for x in data["counters"]], edep, None, {}
            )
        )
        metas.append(meta)
    if not parts:
        raise SystemExit(f"no parts of sample {name} found")
    for m in metas:  # every part against the configuration built here
        for key, want in expected.items():
            if m.get(key) != want:
                raise SystemExit(
                    f"sample {name}: metadata {key}={m.get(key)!r} differs from expected {want!r}"
                )
        if m["git_sha"] != expected_sha:
            raise SystemExit(f"sample {name}: produced at SHA {m['git_sha']}")
    n_total = metas[0]["n_total"]
    cfg, layout = sample_config(a, name, 1)
    raw = merge_partials(parts, n_total, len(cfg.scoring))  # fails closed on gaps/overlaps
    info = {
        "name": name,
        "backend": metas[0]["backend"],
        "precision": metas[0]["precision"],
        "n_histories": n_total,
        "frozen_histories": SAMPLES[name][2],
        "reduced": n_total < SAMPLES[name][2],
        "n_batches": metas[0]["n_batches"],
        "seed": metas[0]["seed"],
        "parts": len(parts),
        "files": sorted(f for f, v in producers.items() if v[1]["meta"]["name"] == name),
        "counters": raw.counters,
        "valid": not any(raw.counters.values()),
    }
    hpb = n_total // metas[0]["n_batches"]
    grids = {
        g.name: raw.edep_mev[i].reshape((metas[0]["n_batches"], *g.shape)) / hpb
        for i, g in enumerate(cfg.scoring)
    }
    obs = parity.t12_observables_from_grids(grids, cfg.scoring, layout, hpb, metas[0]["precision"])
    return obs, info


def step_t12_compare(a: argparse.Namespace) -> int:
    dirs = [Path(d) for d in a.dirs]
    pairs = [p.split(":") for p in a.pairs.split(",")]
    names = sorted({x for p in pairs for x in p})
    expected = os.environ.get("IONMC_RUN_SHA", "")
    obs, info = {}, {}
    for name in names:
        obs[name], info[name] = load_sample(name, dirs, expected, a)
    out = {}
    ok = all(i["valid"] for i in info.values())
    for x, y in pairs:
        v = parity.t12_compare(obs[x], obs[y])
        out[f"{x}_vs_{y}"] = v
        ok &= v["pass"]
    reduced = any(i["reduced"] for i in info.values())  # each sample against its own frozen count
    doc = {
        "step": "t12-compare",
        "energy_mev": a.energy,
        "samples": info,
        "t12": out,
        "pass": ok,
        "lateral_bin_mm": a.lateral_bin,
        "reduced": reduced,
        "frozen_histories": {n: i["frozen_histories"] for n, i in info.items()},
        "histories": {n: i["n_histories"] for n, i in info.items()},
    }
    return emit(doc)


# -- slab observables (T8, T14) ------------------------------------------------------------------
LATERAL_BIN_MM = 0.2
SLAB_MM = 1.0
HALF_WORLD_MM = 60.0


def fermi_eyges_a2(path: mc.ProtonPath, z_mm: float) -> float:
    """A2(z) = integral of (z - z')^2 T_dM(z') dz' over [0, z] (projected, mm^2)."""
    gx, gw = np.polynomial.legendre.leggauss(200)
    s = 0.5 * (gx + 1.0)
    pts = z_mm * s**4
    jac = 4.0 * z_mm * s**3 * 0.5 * gw
    return float(np.sum(jac * (z_mm - pts) ** 2 * path.power(pts)))


def _thickness(frac: float) -> float:
    """The actual depth ``frac * R1`` of the slab [mm] (no rounding: the exit plane is exact)."""
    return frac * r_csda_mm(150.0)


def _slab_geometry(
    depth: float, voxel: float, shift_z: bool, shift_xy: bool, exit_clip: bool = False
) -> tuple[Any, float]:
    """Water world of the given depth from the source plane z = 0: one voxel (``voxel <= 0``) or a
    grid of cubic voxels, optionally shifted by half a voxel along the beam and laterally.
    With ``exit_clip`` the far face of the world is the plane ``z = depth`` for every voxel size
    and shift (the grid overhangs it; ``VoxelGeometry.z_exit_mm``), so the physical slab does not
    depend on the grid. Returns the geometry and the z of the far face of the world."""
    half = HALF_WORLD_MM
    if voxel <= 0.0:
        return BoxPhantom((-half, -half, 0.0), (2 * half, 2 * half, depth), WATER), depth
    z0 = -0.5 * voxel if shift_z else 0.0
    if exit_clip:
        nz = int(math.ceil((depth - z0) / voxel - 1e-9))
    else:
        nz = int(round(depth / voxel)) + (1 if shift_z else 0)
    nxy = int(round(2 * half / voxel)) + (1 if shift_xy else 0)
    xy0 = -half - (0.5 * voxel if shift_xy else 0.0)
    shape = (nxy, nxy, nz)
    geo = VoxelGeometry(
        (xy0, xy0, z0),
        (voxel,) * 3,
        shape,
        (WATER,),
        np.zeros(shape, dtype=np.int32),
        z_exit_mm=depth if exit_clip else None,
    )
    return geo, (depth if exit_clip else z0 + nz * voxel)


def slab_theta(
    *,
    thickness: float,
    voxel: float,
    shift_z: bool,
    shift_xy: bool,
    smax: float,
    n: int,
    workers: int,
    seed: int,
    timeout: float | None,
    trunc_diag: bool = False,
) -> dict[str, float]:
    """Exit-angle observable: projected rms angle of protons leaving a water slab of the given
    thickness (MCS on, straggling off), from the escape records."""
    geo, z_end = _slab_geometry(thickness, voxel, shift_z, shift_xy, exit_clip=True)
    half = HALF_WORLD_MM
    res = run_cfg(
        energy=150.0,
        geometry=geo,
        scoring=(
            ScoringGrid((-half, -half, 0.0), (2 * half, 2 * half, max(thickness, 1.0)), (1, 1, 1)),
        ),
        seed=seed,
        n=n,
        n_batches=20,
        workers=workers,
        timeout=timeout,
        mcs=True,
        straggling=False,
        max_step=smax,
        diag=DiagnosticsOptions(escape_records=True, track_end_positions=trunc_diag),
        trunc_diag=trunc_diag,
    )
    u = res.diagnostics["escape_direction"]
    tx, ty = u[:, 0] / u[:, 2], u[:, 1] / u[:, 2]
    return {
        "theta_rms": math.sqrt(float(np.mean(tx * tx + ty * ty)) / 2.0),
        "z_end_mm": z_end,
        "n_escaped": len(res.diagnostics["escape_history"]),
        "valid": float(res.valid),
        "control_residual_max": (
            float(res.diagnostics["control_residual"].max()) if trunc_diag else 0.0
        ),
        "control_residual_fraction_above_tolerance": (
            float((res.diagnostics["control_residual"] > 1e-3).mean()) if trunc_diag else 0.0
        ),
        "control_residual_quantiles_50_90_99": (
            [float(x) for x in np.quantile(res.diagnostics["control_residual"], [0.5, 0.9, 0.99])]
            if trunc_diag
            else [0.0, 0.0, 0.0]
        ),
        # truncate-first control: vector sum of the snap displacements of each history [mm]
        "control_displacement_rms_mm": (
            float(np.sqrt(np.mean(np.sum(res.diagnostics["control_displacement"] ** 2, axis=1))))
            if trunc_diag
            else 0.0
        ),
        "control_displacement_quantiles_50_90_99_mm": (
            [
                float(x)
                for x in np.quantile(
                    np.linalg.norm(res.diagnostics["control_displacement"], axis=1),
                    [0.5, 0.9, 0.99],
                )
            ]
            if trunc_diag
            else [0.0, 0.0, 0.0]
        ),
        "sigma_exit_mm": float(
            math.sqrt(
                float(
                    np.mean(
                        res.diagnostics["escape_position_mm"][:, 0] ** 2
                        + res.diagnostics["escape_position_mm"][:, 1] ** 2
                    )
                )
                / 2.0
            )
        ),
    }


def _profile_var(p: np.ndarray, xc: np.ndarray) -> float:
    m = float((p * xc).sum() / p.sum())
    return float((p * (xc - m) ** 2).sum() / p.sum())


def sigma_half_width(path: mc.ProtonPath) -> float:
    """Half width of the lateral scoring grids: 8 sigma of the widest slab (z/R = 0.9) in units of
    the lateral bin."""
    sig = math.sqrt(fermi_eyges_a2(path, 0.9 * path.r1_mm))
    return math.ceil(8.0 * sig / LATERAL_BIN_MM) * LATERAL_BIN_MM


def slab_sigma(
    *,
    fracs: tuple[float, ...],
    voxel: float,
    shift_z: bool,
    shift_xy: bool,
    smax: float,
    n: int,
    workers: int,
    seed: int,
    timeout: float | None,
    path: mc.ProtonPath,
) -> dict[str, dict[str, float]]:
    """The frozen observable of T7/T8/T14: lateral sigma of the deposited energy in fixed
    0.2 mm lateral bins and 1 mm slabs at z/R in ``fracs`` (scoring grids independent of the
    transport voxels; water, MCS on, straggling off, 150 MeV), Sheppard-corrected
    (``sigma^2 - bin^2 / 12``), with the Fermi-Eyges prediction at the slab centre."""
    zs = {f: math.floor(f * path.r1_mm / SLAB_MM) * SLAB_MM for f in fracs}
    depth = max(zs.values()) + SLAB_MM + 3.0
    geo, _ = _slab_geometry(depth, voxel, shift_z, shift_xy)
    w = sigma_half_width(path)
    nl = int(round(2 * w / LATERAL_BIN_MM))
    grids = tuple(
        ScoringGrid(
            (-w, -w, zs[f]),
            (LATERAL_BIN_MM, LATERAL_BIN_MM, SLAB_MM),
            (nl, nl, 1),
            name=f"slab{int(round(100 * f))}",
        )
        for f in fracs
    )
    res = run_cfg(
        energy=150.0,
        geometry=geo,
        scoring=grids,
        seed=seed,
        n=n,
        n_batches=10,
        workers=workers,
        timeout=timeout,
        mcs=True,
        straggling=False,
        max_step=smax,
    )
    out: dict[str, dict[str, float]] = {}
    for f, g in zip(fracs, grids, strict=True):
        e = np.asarray(res.grid(g.name).energy_mev)[:, :, 0]  # (nx, ny)
        xc = g.origin_mm[0] + (np.arange(nl) + 0.5) * LATERAL_BIN_MM
        px, py = e.sum(axis=1), e.sum(axis=0)

        sigma2 = 0.5 * (_profile_var(px, xc) + _profile_var(py, xc)) - LATERAL_BIN_MM**2 / 12.0
        zc = zs[f] + 0.5 * SLAB_MM
        a2 = fermi_eyges_a2(path, zc)
        out[str(f)] = {
            "sigma_mm": math.sqrt(sigma2),
            "z_centre_mm": zc,
            "fermi_eyges_sigma_mm": math.sqrt(a2),
            "sigma_over_fermi_eyges": math.sqrt(sigma2 / a2),
            "valid": float(res.valid),
        }
    return out


def step_t8(a: argparse.Namespace) -> int:
    path = mc.ProtonPath(TABLES, WATER, 150.0)
    x = _thickness(0.5)
    quad_theta = math.sqrt(path.theta2_quadrature(x))
    steps = (0.1, 0.5, 1.0, 5.0)
    timings: dict[str, float] = {}
    theta, sigma = {}, {}
    for i, sm in enumerate(steps):
        t0 = time.perf_counter()
        theta[str(sm)] = slab_theta(
            thickness=x,
            voxel=-1.0,
            shift_z=False,
            shift_xy=False,
            smax=sm,
            n=a.n,
            workers=a.workers,
            seed=a.seed + i,
            timeout=a.timeout,
        )
        timings[f"theta_smax{sm}_s"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        sigma[str(sm)] = slab_sigma(
            fracs=(0.9,),
            voxel=-1.0,
            shift_z=False,
            shift_xy=False,
            smax=sm,
            n=a.n,
            workers=a.workers,
            seed=a.seed + 10 + i,
            timeout=a.timeout,
            path=path,
        )["0.9"]
        timings[f"sigma_smax{sm}_s"] = time.perf_counter() - t0
    th = [theta[k]["theta_rms"] for k in theta]
    sg = [sigma[k]["sigma_mm"] for k in sigma]
    spread_t, spread_s = max(th) / min(th) - 1.0, max(sg) / min(sg) - 1.0
    vs_quad = max(abs(v / quad_theta - 1.0) for v in th)
    ok_t = bool(
        spread_t <= 0.005 and vs_quad <= 0.005 and all(r["valid"] == 1.0 for r in theta.values())
    )
    ok_s = bool(spread_s <= 0.01 and all(r["valid"] == 1.0 for r in sigma.values()))
    out = {
        "theta_0.5R1": {
            "thickness_mm": x,
            "runs": theta,
            "pairwise_spread": spread_t,
            "tolerance_pairwise": 0.005,
            "quadrature": quad_theta,
            "max_dev_from_quadrature": vs_quad,
            "pass": ok_t,
        },
        "sigma_0.9R_deposit": {
            "runs": sigma,
            "pairwise_spread": spread_s,
            "tolerance_pairwise": 0.01,
            "pass": ok_s,
        },
    }
    return finish(
        {"step": "t8", "t8": out, "timings": timings, "pass": ok_t and ok_s}, 1_000_000, a.n
    )


def step_t14(a: argparse.Namespace) -> int:
    path = mc.ProtonPath(TABLES, WATER, 150.0)
    x05 = _thickness(0.5)
    rows: dict[str, Any] = {}
    timings: dict[str, float] = {}
    i = 0
    for voxel in (0.5, 1.0, 2.0, 5.0):
        for sz, sxy in ((False, False), (True, True)):
            key = f"voxel{voxel}_shift{int(sz)}"
            smax = min(1.0, voxel)
            i += 1
            t0 = time.perf_counter()
            th = slab_theta(
                thickness=x05,
                voxel=voxel,
                shift_z=sz,
                shift_xy=sxy,
                smax=smax,
                n=a.n,
                workers=a.workers,
                seed=a.seed + i,
                timeout=a.timeout,
            )
            th["theta_over_quadrature"] = th["theta_rms"] / math.sqrt(
                path.theta2_quadrature(th["z_end_mm"])
            )
            sg = slab_sigma(
                fracs=(0.5, 0.9),
                voxel=voxel,
                shift_z=sz,
                shift_xy=sxy,
                smax=smax,
                n=a.n,
                workers=a.workers,
                seed=a.seed + 100 + i,
                timeout=a.timeout,
                path=path,
            )
            timings[key + "_s"] = time.perf_counter() - t0
            rows[key] = {"theta": th, "sigma": sg}
    # frozen negative control: the hinge angle sampled for the truncated length (diagnostic
    # switch) must change theta_rms by < 0.5 % at 1 mm voxels (same seed as the default run)
    t0 = time.perf_counter()
    kw = dict(
        thickness=x05,
        voxel=1.0,
        shift_z=False,
        shift_xy=False,
        smax=1.0,
        n=a.n,
        workers=a.workers,
        seed=a.seed + 1000,
        timeout=a.timeout,
    )
    default_1mm = slab_theta(**kw)
    control_1mm = slab_theta(**kw, trunc_diag=True)
    timings["negative_control_s"] = time.perf_counter() - t0
    control = {
        "theta_rms_default": default_1mm["theta_rms"],
        "theta_rms_control": control_1mm["theta_rms"],
        "relative_change": control_1mm["theta_rms"] / default_1mm["theta_rms"] - 1.0,
        "bound": 0.005,
        "valid": bool(default_1mm["valid"] == 1.0 and control_1mm["valid"] == 1.0),
        # truncate-first control: the snap moves only the position (never the direction); the
        # positional bias is the vector sum of a history's snap displacements, whose RMS over the
        # histories must be below 1 % of the lateral sigma at the exit plane (ENFORCED); the
        # quantiles of |sum| and of the per-step displacement/step are informative
        "snap_displacement": {
            "rms_of_vector_sum_mm": control_1mm["control_displacement_rms_mm"],
            "sigma_exit_mm": control_1mm["sigma_exit_mm"],
            "relative_to_sigma": control_1mm["control_displacement_rms_mm"]
            / control_1mm["sigma_exit_mm"],
            "bound": 0.01,
            "quantiles_50_90_99_of_abs_sum_mm": control_1mm[
                "control_displacement_quantiles_50_90_99_mm"
            ],
            "per_step_over_step_informative": {
                "max_over_histories": control_1mm["control_residual_max"],
                "fraction_of_histories_above_1e-3": control_1mm[
                    "control_residual_fraction_above_tolerance"
                ],
                "quantiles_50_90_99": control_1mm["control_residual_quantiles_50_90_99"],
            },
        },
    }
    control["pass"] = bool(
        abs(control["relative_change"]) < 0.005
        and control["valid"]
        and control["snap_displacement"]["relative_to_sigma"]
        < control["snap_displacement"]["bound"]
    )
    # the frozen observable: the RAW exit theta_rms of the same 0.5 R1 slab (the world ends at
    # z = thickness for every voxel size and shift), compared pairwise; the ratio to the U5
    # quadrature is the second frozen bound (within 0.5 %)
    th_raw = [rows[k]["theta"]["theta_rms"] for k in rows]
    ths = [rows[k]["theta"]["theta_over_quadrature"] for k in rows]
    s5 = [rows[k]["sigma"]["0.5"]["sigma_over_fermi_eyges"] for k in rows]
    s9 = [rows[k]["sigma"]["0.9"]["sigma_over_fermi_eyges"] for k in rows]
    abs_s5 = [rows[k]["sigma"]["0.5"]["sigma_mm"] for k in rows]
    abs_s9 = [rows[k]["sigma"]["0.9"]["sigma_mm"] for k in rows]
    res = {
        "theta_pairwise_spread": max(th_raw) / min(th_raw) - 1.0,
        "theta_max_dev_from_quadrature": max(abs(t - 1.0) for t in ths),
        "sigma_pairwise_spread_0.5R": max(abs_s5) / min(abs_s5) - 1.0,
        "sigma_pairwise_spread_0.9R": max(abs_s9) / min(abs_s9) - 1.0,
        "sigma_max_dev_from_fermi_eyges": max(abs(t - 1.0) for t in s5 + s9),
    }
    valid = all(
        r["theta"]["valid"] == 1.0 and all(v["valid"] == 1.0 for v in r["sigma"].values())
        for r in rows.values()
    )
    res["pass"] = bool(
        res["theta_pairwise_spread"] <= 0.005
        and res["theta_max_dev_from_quadrature"] <= 0.005
        and max(res["sigma_pairwise_spread_0.5R"], res["sigma_pairwise_spread_0.9R"]) <= 0.01
        and res["sigma_max_dev_from_fermi_eyges"] <= 0.02
        and valid
        and control["pass"]
    )
    doc = {
        "step": "t14",
        "runs": rows,
        "t14": res,
        "pass": res["pass"],
        "negative_control": control,
        "timings": timings,
    }
    return finish(doc, 1_000_000, a.n)


# -- T9 ------------------------------------------------------------------------------------------
def step_t9(a: argparse.Namespace) -> int:
    """T9; ``--physics`` selects a diagnostic variant (``no-straggling``, ``no-mcs``; not part of
    the suites) and ``--straggling-model`` the sampler (default ``bohr_gamma_v1``; the old
    ``bohr_gauss_clamped_gamma_v1`` for the comparison)."""
    mcs, straggling = a.physics != "no-mcs", a.physics != "no-straggling"
    e = 150.0
    depth = 1.3 * r_csda_mm(e)
    nz = int(math.ceil(depth))
    geo = BoxPhantom((-50.0, -50.0, 0.0), (100.0, 100.0, depth), WATER)
    grid = (ScoringGrid((-40.0, -40.0, 0.0), (80.0, 80.0, 1.0), (1, 1, nz), name="idd"),)
    layout = parity.T12Layout("idd", (), r_csda_mm(e), 1.0, 80.0)
    cases = {
        "s1.0_f0.02": (1.0, 0.02),
        "s0.5_f0.02": (0.5, 0.02),
        "s0.1_f0.02": (0.1, 0.02),
        "s1.0_f0.005": (1.0, 0.005),
    }
    obs = {}
    for i, (name, (s, f)) in enumerate(cases.items()):
        res = run_cfg(
            energy=e,
            geometry=geo,
            scoring=grid,
            seed=a.seed + i,
            n=a.n,
            n_batches=20,
            workers=a.workers,
            timeout=a.timeout,
            max_step=s,
            frac=f,
            mcs=mcs,
            straggling=straggling,
            straggling_model=a.straggling_model,
        )
        obs[name] = parity.t12_observables(res, layout)
    ref = obs["s1.0_f0.02"]
    out = {}
    ok = True
    for name, o in obs.items():
        if name == "s1.0_f0.02":
            continue
        v = parity.t12_compare(ref, o)
        dr = abs(v["scalars"]["r80_mm"]["a"] - v["scalars"]["r80_mm"]["b"])
        idd = v["arrays"]["idd"]
        passed = bool(dr <= 0.1 and idd["pass"])  # grouped profile verdict
        # where the profiles differ: bins (1 mm = depth in mm) with |z| > 3 and the sign of
        # z = (reference - test) / se for bins above 1 % of the maximum
        ma, mt = ref.arrays["idd"].mean(axis=0), o.arrays["idd"].mean(axis=0)
        se = np.sqrt(
            ref.arrays["idd"].var(axis=0, ddof=1) / ref.n_batches
            + o.arrays["idd"].var(axis=0, ddof=1) / o.n_batches
        )
        sel = 0.5 * (ma + mt) > 0.01 * (0.5 * (ma + mt)).max()
        zz = np.where(sel & (se > 0), (ma - mt) / np.where(se > 0, se, 1.0), 0.0)
        far = np.nonzero(np.abs(zz) > 3.0)[0]
        out[name] = {
            "delta_r80_mm": dr,
            "idd_chi2": idd,
            "z_bins_above_3": {"depth_mm": far.tolist(), "z": [float(zz[i]) for i in far]},
            "pass": passed,
        }
        ok &= passed
    doc = {
        "step": "t9",
        "t9": out,
        "pass": ok,
        "variant": a.physics,
        "straggling_model": a.straggling_model,
    }
    return finish(doc, None, a.n)


# -- T10 -----------------------------------------------------------------------------------------
DIRECTIONS = {
    "+x": (1.0, 0.0, 0.0),
    "+y": (0.0, 1.0, 0.0),
    "+z": (0.0, 0.0, 1.0),
    "-z": (0.0, 0.0, -1.0),
    "(1,1,0)": (1.0, 1.0, 0.0),
    "(1,1,1)": (1.0, 1.0, 1.0),
}
T10_BATCHES = 10


def _t10_grid(start: np.ndarray, u: np.ndarray, r: float) -> ScoringGrid:
    """1 mm scoring grid around the path of the beam (axis-aligned bounding box plus 15 mm)."""
    end = start + r * u
    lo = np.floor(np.minimum(start, end) - 15.0)
    hi = np.ceil(np.maximum(start, end) + 15.0)
    shape = tuple(int(x) for x in (hi - lo))
    return ScoringGrid(tuple(float(x) for x in lo), (1.0, 1.0, 1.0), shape)  # type: ignore[arg-type]


def step_t10(a: argparse.Namespace) -> int:
    e = 150.0
    r = r_csda_mm(e)
    geo = BoxPhantom((-100.0, -100.0, -100.0), (200.0, 200.0, 200.0), WATER)
    nb = int(math.ceil(1.3 * r / 0.1))  # 0.1 mm bins of the projected depth
    # one voxel covering the whole world: every deposit, so its per-batch value is the total energy
    # deposited (the path box of the depth-dose omits large-angle tails beyond its margin)
    world = ScoringGrid((-100.0, -100.0, -100.0), (200.0, 200.0, 200.0), (1, 1, 1), name="world")
    obs: dict[str, Any] = {}
    for i, (name, d) in enumerate(DIRECTIONS.items()):
        u = np.array(d) / np.linalg.norm(d)
        start = -0.5 * r * u
        grid = _t10_grid(start, u, r)
        res = run_cfg(
            energy=e,
            geometry=geo,
            scoring=(grid, world),  # path box for the depth-dose, whole world for the energy
            seed=a.seed + i,
            n=a.n,
            n_batches=T10_BATCHES,
            workers=a.workers,
            timeout=a.timeout,
            position=tuple(float(x) for x in start),
            direction=d,
            diag=DiagnosticsOptions(track_end_positions=True),
        )
        be = np.asarray(res.grids[0].batch_energy_mev)  # (B, nx, ny, nz) per primary
        idd = parity.projected_idd(be, grid, tuple(float(x) for x in start), d, n_bins=nb)
        idd_1mm = idd[:, : nb // 10 * 10].reshape(T10_BATCHES, -1, 10).sum(axis=2)
        pos = res.diagnostics["end_position_mm"] - start
        t_end = pos @ u
        obs[name] = {
            "r80_mm": parity.r80_of(idd.mean(axis=0), 0.1),
            "idd_1mm": idd_1mm,
            "precision": res.precision,
            "total_b": np.asarray(res.grid("world").batch_energy_mev).reshape(T10_BATCHES),
            "in_path_box_mev": float(be.sum() / T10_BATCHES),
            "mean_end_depth_mm": float(t_end.mean()),
            "counters": res.counters.as_dict(),
            "valid": res.valid,
        }
    ref = obs["+z"]
    out = {}
    ok = True
    for name, o in obs.items():
        if name == "+z":
            continue
        v = parity.t12_compare(
            parity.T12Observables({"idd_1mm": ref["idd_1mm"]}, {}, T10_BATCHES, a.n // T10_BATCHES),
            parity.T12Observables({"idd_1mm": o["idd_1mm"]}, {}, T10_BATCHES, a.n // T10_BATCHES),
        )
        d_r80 = abs(o["r80_mm"] - ref["r80_mm"])
        # the total energy is fixed by energy conservation (no variance): the opt-in 1e-9 shortcut
        # the total energy is fixed by energy conservation: the rule of plan footnote 19 (the
        # deterministic T4 bound of the less precise sample plus 3.5 combined standard errors)
        e_z = parity.deterministic_scalar_verdict(
            ref["total_b"], o["total_b"], ref["precision"], o["precision"]
        )
        oblique = name.startswith("(")
        chi_ok = bool(v["arrays"]["idd_1mm"]["pass"])  # grouped profile verdict, inconclusive fails
        r80_ok = d_r80 <= 0.3
        energy_ok = bool(e_z["pass"])
        # frozen: axis permutations -> chi-square p > 0.001 (obliques: the chi-square is only
        # informative); obliques -> |dR80| <= 0.3 mm; all: total energy by the deterministic
        # rule of footnote 19, no counters
        crit = r80_ok if oblique else chi_ok
        passed = bool(crit and energy_ok and not any(o["counters"].values()) and o["valid"])
        out[name] = {
            "r80_mm": o["r80_mm"],
            "delta_r80_mm": d_r80,
            "r80_ok": r80_ok,
            "idd_chi2": v["arrays"]["idd_1mm"],
            "permutation": v["permutation"],
            "total_energy": e_z,
            "mean_end_depth_mm_informative": o["mean_end_depth_mm"],
            "counters": o["counters"],
            "pass": passed,
        }
        ok &= passed
    doc = {"step": "t10", "reference_r80_mm": ref["r80_mm"], "t10": out, "pass": ok}
    return finish(doc, 1_000_000, a.n)


STEPS = {
    "t-r1": step_repeat,
    "t12-python-sample": step_t12_python_sample,
    "t12-accelerated-samples": step_t12_accelerated,
    "t12-compare": step_t12_compare,
    "t1": step_t1,
    "t2": step_t2,
    "t13": step_t13,
    "t8": step_t8,
    "t9": step_t9,
    "t10": step_t10,
    "t14": step_t14,
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--workers", default="1", help="integer or 'auto' (all cores)")
    ap.add_argument("--runs", default="python:float64:400", help="t-r1: backend:precision:n,...")
    ap.add_argument(
        "--physics",
        choices=("default", "no-straggling", "no-mcs"),
        default="default",
        help="t9 diagnostic variant (not part of the suites)",
    )
    ap.add_argument(
        "--straggling-model",
        default="bohr_gamma_v1",
        help="t9: straggling sampler (bohr_gauss_clamped_gamma_v1 or bohr_gamma_v1)",
    )
    ap.add_argument("--part", default="1/1", help="t12-python-sample: history range i/n")
    ap.add_argument("--samples", default="cpu32", help="t12-accelerated-samples: names")
    ap.add_argument("--out-dir", default="samples", help="directory of the sample files")
    ap.add_argument("--dirs", nargs="+", default=["."], help="t12-compare: archive directories")
    ap.add_argument("--n", type=int, default=1_000_000)
    ap.add_argument("--k", type=int, default=256)
    ap.add_argument("--energy", type=float, default=150.0)
    ap.add_argument(
        "--seed-base",
        "--seed",
        dest="seed",
        type=int,
        default=DEFAULT_SEED_BASE,
        help="base of every seed of the step (derivation: see the module docstring)",
    )
    ap.add_argument("--timeout", type=float, default=None)
    ap.add_argument("--backend", default="warp-cpu")
    ap.add_argument("--mode", choices=("workers", "chunks"), default="workers")
    ap.add_argument("--pairs", default="python:cpu32,python:cpu64,cpu32:cpu64")
    ap.add_argument("--scale", type=float, default=1.0, help="T12 history-count factor")
    ap.add_argument("--lateral-bin", type=float, default=0.2)
    ap.add_argument("--half-width", type=float, default=20.0)
    args = ap.parse_args(argv)
    global SEED_BASE
    SEED_BASE = args.seed
    args.workers = (os.cpu_count() or 1) if args.workers == "auto" else int(args.workers)
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
