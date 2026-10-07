"""Validation steps of the V3-005 acceptance plan, slice A (suite ``lv5``; decision 0041).

Usage::

    python validation/scripts/transport/steps_v5.py <step> [--scale F] [--seed-base N] ...

Rows (``validation/plans/v3-005-acceptance.md``): N1, V1, V1b, D6 (deterministic, the table JSON and
``nuclear_checks.py``), V2 (primary survival, 3 energies x 3 shards of 8e4, pooled by ``v2-combine``),
V2-probe (common-random-number step-size probes, pooled by ``v2-probe-combine``), V3-LV (amended
energy balance), V4 and V4b (offline event sampler), X1 (paired nuclear on/off), E1 (exploratory IDD),
R1 (``nuclear=False`` regression: the A16 digest at baseline f3a1dd62 with the intended-change set
``None`` and the T1 trace) and the throughput measurement. All statistical steps run the python
backend with ``nuclear=True`` (decision 0041) in one process; the 1-vs-N worker partition check
``v3-workers-partition`` exists but is deferred by ``run_suite.py`` (plan, "Deferred multiprocess
checks"). Every step prints one JSON document between ``#JSON-BEGIN`` and ``#JSON-END`` (the
conventions of ``steps.py`` and ``steps_v4.py``): ``format``, ``frozen_histories``, ``histories``,
``reduced`` (``--scale`` < 1 gives a labelled, non-conformant run), ``pass`` and the criterion values;
every nuclear document carries the table id and the sha256 of the table files. Every invalid result
or non-zero fail-closed counter fails the step.

Seeds (plan, Seeds): ``seed = base + 1000 * r_index + shard`` with ``--seed-base`` the lv5 base
(qualification 20421004; the rehearsal family is 2043xxxx) and shards counted from 0. ``R_INDEX``:
V2-100 1, V2-150 2, V2-200 3, V2-probe-s05 4, V2-probe-fE 5, V3-LV 6, X1 7, E1 8, V4 9. V3-LV uses
its six runs as shards 0..5 (material-major: water, tissue, bone; 150 then 250 MeV); X1 runs on and
off with the same seed (shard 0); V4 and V4b use ``nuclear_checks.V4_BASE_SEED`` patched to
``base + 9000`` (case ``i`` adds ``i``). N1, V1, V1b, D6 are deterministic and R1 uses the fixed A16
digest seed (20351004) and the T1 seed ``--seed-base``.

Definitions made here where the plan leaves them open (all recorded in the documents):

* V2: the primary-survival profile is the ``generation="primary"`` proton fluence in 1 mm depth bins
  of a water box (MCS and straggling off, s_max 1 mm, f_E 0.02) divided, per batch, by its first bin.
  The reference is exp(-int Sigma/S dE) along the deterministic CSDA path (same tables), averaged
  over each 1 mm bin and divided by the first bin's average. The row is the pooled TOST over every
  batch of the 3 shards with margin 0.003 and alpha = 0.05 per depth (z = 1.645), at every 10 mm
  depth (bin index = depth / 1 mm) up to the CSDA range minus 5 mm.
* V2-probe: the same physics as V2 (MCS and straggling off) at 150 MeV; a shard runs both variants
  with one seed; the criterion is the pooled batch mean of ``S_a - S_b`` (each normalised to its
  first bin); valid only if the pooled standard error at the deepest depth is <= 7e-4.
* V3-LV: homogeneous water, ``MUSCLE_SKELETAL_ICRP`` (tissue) and ``BONE_COMPACT_ICRU`` (bone)
  boxes of depth 1.1 R(E) of the Bethe table; the relative residuals of the amended balance and of
  the per-grid identity are <= 1e-12 and every counter is 0. The ``nuclear_binding`` tally is recomputed independently from the
  reference's per-target event, light-product and residual (Z_r, A_r) counts with AME2020 masses
  read from the cached AME text: ``|recomputed - tally| <= 1e-9 |tally| + 1e-9`` MeV (gating).
* X1: peak = the pooled-off-profile peak bin, plateau = mean over bin centres 20-60 mm; the paired
  statistic is the batch difference of ``peak / plateau`` (on - off), z = -mean / sem, pass iff the
  difference is negative and z > 3.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import math
import os
import sys
import time
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np
from numpy.typing import NDArray

import nuclear_checks as nc
import steps as base
import steps_v4 as v4
from ionmc import materials as M
from ionmc.config import DiagnosticsOptions, PhysicsOptions, RunOptions, SimulationConfig
from ionmc.data import cache
from ionmc.geometry import BoxPhantom
from ionmc.nuclear.tables import NuclearTable
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.simulation import Result, Simulation
from ionmc.sources import PencilBeamSource

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FORMAT = 1
QUALIFICATION_SEED_BASE = 20421004
REHEARSAL_SEED_BASE = 20431004
TABLE_ID = "3bcf146e38dd2b5581bd1ff245c127a7d059e6789421a5d3c6760974f2784504"
A16_R1_DIGEST = "c862edf799dcb83542e5071219b5703c7665f8f792f6bbc49ab32918418b1f8f"
A16_R1_BASELINE = "f3a1dd62ea2f57a3f4c07999935f317b3c044871"
R_INDEX = {"v2-100": 1, "v2-150": 2, "v2-200": 3, "v2-probe-s05": 4, "v2-probe-fe": 5,
           "v3-lv": 6, "x1": 7, "e1": 8, "v4": 9}  # fmt: skip
STEP_LIMIT_S = 3300
MARGIN = 1.25  # planned step time <= STEP_LIMIT_S / MARGIN

V2_ENERGIES = (100.0, 150.0, 200.0)
V2_SHARDS, V2_SHARD_N = 3, 80_000
V2_TOL, V2_ALPHA, V2_BATCHES = 0.003, 0.05, 20
PROBE_TOL, PROBE_SIGMA_MAX = 0.002, 7e-4
PROBES = {  # name -> (shards, histories per shard, (a: s_max mm, f_E), (b: s_max mm, f_E))
    "s05": (2, 50_000, (0.5, 0.02), (1.0, 0.02)),
    "fe": (3, 34_000, (1.0, 0.005), (1.0, 0.02)),
}
V3_HISTORIES = {150.0: 20_000, 250.0: 10_000}
V3_MATERIALS = (("water", M.WATER), ("tissue", M.MUSCLE_SKELETAL_ICRP), ("bone", M.BONE_COMPACT_ICRU))
X1_N, X1_BATCHES, X1_Z = 20_000, 20, 3.0
X1_PLATEAU_MM = (20.0, 60.0)
E1_N, E1_BATCHES = 100_000, 20
V4_EVENTS = 100_000
THROUGHPUT_N, THROUGHPUT_ENERGIES = 200, (150.0, 250.0)
PARTITION_N, PARTITION_WORKERS = 3000, 3


def seed_of(row: str, shard: int = 0) -> int:
    return base.SEED_BASE + 1000 * R_INDEX[row] + shard


def finish5(doc: dict[str, Any], frozen: Any, used: Any, reduced: bool) -> int:
    doc = {"format": FORMAT, **doc}
    return v4.finish4(doc, frozen, used, reduced)


def table_record() -> dict[str, Any]:
    """Id and file hashes of the nuclear table every nuclear step loads."""
    cdir = cache.resolve_cache_dir(None)
    out: dict[str, Any] = {"table_id": TABLE_ID}
    for ext in ("json", "npz"):
        p = cdir / "derived" / f"nuclear-proton-{TABLE_ID}.{ext}"
        out[f"{ext}_sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def nuc_config(
    *, energy: float, n: int, seed: int, geometry: Any, grid: ScoringGrid, nuclear: bool = True,
    mcs: bool = True, straggling: bool = True, max_step: float = 1.0, frac: float = 0.02,
    n_batches: int = 20, tallies: tuple[TallyRequest, ...] = (), workers: int = 1,
    timeout: float | None = None,
) -> SimulationConfig:  # fmt: skip
    """Python backend, float64, all EM physics as given; ``nuclear`` selects decision 0041."""
    return SimulationConfig(
        source=PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), energy),
        geometry=geometry,
        scoring=(grid,),
        tallies=tallies,
        physics=PhysicsOptions(
            nuclear=nuclear, nuclear_table_id=TABLE_ID if nuclear else None,
            stopping=BetheStoppingSource(), straggling=straggling, multiple_scattering=mcs,
            max_step_mm=max_step, max_energy_loss_fraction=frac,
        ),
        run=RunOptions(
            backend="python", precision="float64", seed=seed, n_histories=n,
            n_batches=n_batches, cpu_workers=workers, worker_timeout_s=timeout,
        ),
        diagnostics=DiagnosticsOptions(),
    )  # fmt: skip


def clean(res: Result) -> bool:
    return bool(res.valid and not res.counters.any_nonzero)


def depth_box(
    mat: Any, energy: float, dz: float = 1.0, half_mm: float = 40.0
) -> tuple[BoxPhantom, ScoringGrid, int, float]:
    """Box of depth 1.1 R (Bethe table of ``mat``) and a ``(1, 1, nz)`` depth grid named ``dose``."""
    r_mm = 10.0 * float(BetheStoppingSource().table(mat, PROTON).range_at(energy)) / mat.density_g_cm3
    nz = int(math.ceil(1.1 * r_mm / dz))
    geo = BoxPhantom((-half_mm, -half_mm, 0.0), (2 * half_mm, 2 * half_mm, nz * dz), mat)
    grid = ScoringGrid((-half_mm, -half_mm, 0.0), (2 * half_mm, 2 * half_mm, dz), (1, 1, nz),
                       name="dose")  # fmt: skip
    return geo, grid, nz, r_mm


PRIMARY = (TallyRequest("fl_prim", "dose", "fluence", species=("proton",), generation="primary"),)


EDEP = (TallyRequest("edep", "dose", "edep"),)


def edep_batches(res: Result) -> NDArray[np.float64]:
    """Per-batch energy deposit per history in the depth bins, ``[B, nz]`` (tally ``edep``)."""
    plan = res.effective_config.channels
    assert plan is not None
    return v4.batch_values(res, {q.name: q.numerator for q in plan.quantities}["edep"])


def primary_ratio(res: Result) -> NDArray[np.float64]:
    """Per-batch primary fluence profile divided by its first bin, ``[B, nz]``."""
    plan = res.effective_config.channels
    assert plan is not None
    idx = {q.name: q.numerator for q in plan.quantities}
    f = v4.batch_values(res, idx["fl_prim"])
    return f / f[:, :1]


def reference_survival(sim: Simulation, energy: float, nz: int, dz: float) -> NDArray[np.float64]:
    """Bin-averaged exp(-int Sigma/S dE) along the deterministic CSDA path, divided by bin 0."""
    eff = sim.effective
    assert eff.nuclear is not None
    rows, tab = eff.nuclear.rows[0], eff.tables
    e = np.geomspace(energy, eff.requested.physics.e_cut_mev, 40001)  # descending
    f = np.array([rows.sigma_at(x) / tab.stopping_mass(0, x) for x in e])
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * -np.diff(e))))
    r = np.array([tab.range_g_cm2(0, x) for x in e])
    z = (r[0] - r) * 10.0 / M.WATER.density_g_cm3  # mm along the path, increasing
    sub = (np.arange(nz * 20) + 0.5) * dz / 20.0
    s = np.where(sub <= z[-1], np.exp(-np.interp(sub, z, cum)), 0.0).reshape(nz, 20).mean(axis=1)
    return s / s[0]


def depth_bins(energy: float, nz: int, dz: float = 1.0) -> list[int]:
    """Bin indices of every 10 mm depth up to the CSDA range minus 5 mm."""
    r = base.r_csda_mm(energy)
    return [int(z / dz) for z in np.arange(10.0, r - 5.0 + 1e-9, 10.0) if int(z / dz) < nz]


def content_digest(doc: dict[str, Any]) -> str:
    """sha256 of the canonical JSON (sorted keys, compact separators) of ``doc`` without its
    ``content_sha256`` field; the document is normalised through one JSON round trip first, so
    the digest of a written file equals the digest recomputed from the parsed file."""
    norm = json.loads(json.dumps(doc, default=base._json))
    norm.pop("content_sha256", None)
    text = json.dumps(norm, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()


def bindings(a: argparse.Namespace) -> dict[str, Any]:
    """What a partial is bound to: the run (SHA, suite), the nuclear table and the seed base."""
    rec = table_record()
    return {"format": FORMAT, "git_sha": os.environ.get("IONMC_RUN_SHA", "unknown"),
            "run_sha": os.environ.get("IONMC_RUN_SHA", "unknown"),
            "suite": os.environ.get("IONMC_RUN_SUITE", "unknown"), "table_id": rec["table_id"],
            "table_npz_sha256": rec["npz_sha256"], "seed_base": base.SEED_BASE,
            "scale": a.scale}  # fmt: skip


def write_partial(a: argparse.Namespace, name: str, doc: dict[str, Any]) -> str:
    """Write a hash-sealed partial: the run bindings plus ``content_sha256`` over the rest."""
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    doc = {**bindings(a), **doc}
    doc["content_sha256"] = content_digest(doc)
    path = out / f"{name}.json"
    path.write_text(json.dumps(doc, sort_keys=True, default=base._json))
    # the host-runner record retains stdout: the orchestrator writes the partials manifest of the
    # combine steps from these lines (name -> content_sha256), an attestation the archive's own
    # files cannot forge
    print(f"PARTIAL {path.name} {doc['content_sha256']}", flush=True)
    return str(path)


def scaled_n(n: int, scale: float) -> int:
    return v4.scaled(n, scale, 2000, V2_BATCHES)


# -- throughput ----------------------------------------------------------------------------------
def step_throughput(a: argparse.Namespace) -> int:
    """Measured python nuclear-on throughput (1 process, all physics on, 200 histories incl. setup)."""
    out: dict[str, Any] = {}
    ok = True
    for e in THROUGHPUT_ENERGIES:
        geo, grid, _, _ = depth_box(M.WATER, e)
        t0 = time.perf_counter()
        res = Simulation(nuc_config(energy=e, n=THROUGHPUT_N, seed=a.seed, geometry=geo, grid=grid,
                                    n_batches=2, timeout=a.timeout)).run()  # fmt: skip
        wall = time.perf_counter() - t0
        rate = THROUGHPUT_N / wall
        ok &= clean(res)
        out[f"{e:g}"] = {
            "histories": THROUGHPUT_N, "wall_s": wall, "hist_per_s": rate,
            "histories_per_3300s_step_with_25pct_margin": int(STEP_LIMIT_S / MARGIN * rate),
            "counters_clean": clean(res),
        }  # fmt: skip
    rate_min = min(v["hist_per_s"] for v in out.values())
    plan = {
        f"v2-{e:g}-shard": {"n": V2_SHARD_N,
                            "planned_s": V2_SHARD_N / out[f"{150.0 if e <= 150 else 250.0:g}"]["hist_per_s"]}
        for e in V2_ENERGIES
    }  # fmt: skip
    doc = {"step": "lv5-throughput", "table": table_record(), "throughput": out,
           "rate_min_hist_per_s": rate_min, "step_limit_s": STEP_LIMIT_S, "margin": MARGIN,
           "v2_shard_plan_at_measured_rate": plan, "pass": bool(ok)}  # fmt: skip
    return finish5(doc, THROUGHPUT_N, THROUGHPUT_N, False)


# -- N1, V1, V1b, D6 ------------------------------------------------------------------------------
def step_n1(a: argparse.Namespace) -> int:
    nc.v4_checks = lambda *_a, **_k: {"skipped": "V4 and V4b run in the step v4-v4b"}  # type: ignore[assignment]
    s = nc.run_checks(None, TABLE_ID)
    d6 = s["D6"]
    ok = bool(s["N1"]["pass"] and s["V1"]["pass"] and d6["ceiling_pass"])
    doc = {
        "step": "n1", "table": table_record(),
        "N1": {k: s["N1"][k] for k in ("tolerance", "max_rel", "grid_points", "pass")},
        "N1_materials_max_rel": {k: v["max_rel"] for k, v in s["N1"]["materials"].items()},
        "V1": {k: s["V1"][k] for k in ("tolerance", "max_rel", "extension_continuity",
                                       "extension_continuity_max_abs", "pass")},
        "V1b": {"gating": False, **s["V1b"]},
        "D6": {**d6, "pass": bool(d6["ceiling_pass"]),
               "rule": "amendment 2: both tiers false is accepted iff the ceiling holds"},
        "pass": ok,
    }  # fmt: skip
    return finish5(doc, None, None, False)


# -- V2 shards and pooling ------------------------------------------------------------------------
def step_v2_shard(a: argparse.Namespace) -> int:
    e, k = a.energy, a.shard
    row = f"v2-{e:g}"
    n = scaled_n(V2_SHARD_N, a.scale)
    geo, grid, nz, _ = depth_box(M.WATER, e)
    cfg = nuc_config(energy=e, n=n, seed=seed_of(row, k), geometry=geo, grid=grid, mcs=False,
                     straggling=False, tallies=PRIMARY, timeout=a.timeout)  # fmt: skip
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    ratio = primary_ratio(res)
    path = write_partial(a, f"{row}-s{k}", {
        "row": row, "energy": e, "shard": k, "seed": seed_of(row, k), "n": n, "valid": clean(res),
        "ratio": ratio.tolist(), "reduced": n < V2_SHARD_N,
    })  # fmt: skip
    doc = {"step": "v2-shard", "row": row, "shard": k, "seed": seed_of(row, k), "table": table_record(),
           "wall_s": wall, "hist_per_s": n / wall, "partial": path, "counters": res.counters.as_dict(),
           "pass": clean(res)}  # fmt: skip
    return finish5(doc, V2_SHARD_N, n, n < V2_SHARD_N)


def read_manifest(path: str | None) -> dict[str, str]:
    """The partials manifest: a JSON object mapping a partial file name to its expected
    ``content_sha256`` (written by the orchestrator from the ``PARTIAL`` lines of the host-runner
    records of the shard steps)."""
    if path is None:
        raise SystemExit("--import-dirs requires --partials-manifest (name -> content_sha256)")
    doc = json.loads(Path(path).read_text())
    if not isinstance(doc, dict) or not all(
        isinstance(k, str) and isinstance(v, str) for k, v in doc.items()
    ):
        raise SystemExit(f"partials manifest {path}: need a JSON object of name -> sha256 string")
    return doc


def load_partials(a: argparse.Namespace, names: list[str]) -> list[dict[str, Any]]:
    """Partial files ``names`` of ``--dirs`` (their ``samples`` directories included), verified
    against their ``content_sha256`` (canonical JSON of the document without that field) and the
    run bindings (SHA, suite, table id and table file hash, seed base, scale, format; the same
    ones for every partial because each must equal the current run); a missing, altered or
    foreign file stops the step. ``a.dirs[0]`` is the current output directory, whose partials
    this run produced; a partial found in any other (imported) directory must be listed in the
    ``--partials-manifest`` with exactly its recomputed digest, because a ``content_sha256``
    stored in the file can be recomputed by whoever alters it. Without imported directories only
    partials of the current directory are accepted. As for the A9 samples of lv4 there is no
    foreign-SHA override."""
    imported = len(a.dirs) > 1
    manifest = read_manifest(getattr(a, "partials_manifest", None)) if imported else {}
    out = []
    for nm in names:
        hits = [(i, p) for i, d in enumerate(a.dirs)
                for p in (Path(d) / nm, Path(d) / "samples" / nm) if p.is_file()]  # fmt: skip
        if len(hits) != 1:
            raise SystemExit(f"partial {nm}: found {len(hits)} copies in {a.dirs} (need exactly 1)")
        idx, hit = hits[0]
        doc = json.loads(hit.read_text())
        digest = content_digest(doc)
        if doc.get("content_sha256") != digest:
            raise SystemExit(f"partial {nm}: content_sha256 does not match the document")
        if idx != 0:
            if nm not in manifest:
                raise SystemExit(f"partial {nm}: imported but not in the partials manifest")
            if manifest[nm] != digest:
                raise SystemExit(
                    f"partial {nm}: digest {digest} differs from the manifest value {manifest[nm]}"
                )
        for key, want in bindings(a).items():
            if doc.get(key) != want:
                raise SystemExit(f"partial {nm}: {key} is {doc.get(key)!r}, expected {want!r}")
        out.append(doc)
    return out


def tost(diff: NDArray[np.float64], tol: float) -> dict[str, Any]:
    """Pooled two-one-sided test of ``|mean(diff)| < tol`` with the batch standard error."""
    b = diff.shape[0]
    mean = diff.mean(axis=0)
    sem = diff.std(axis=0, ddof=1) / math.sqrt(b)
    z = NormalDist().inv_cdf(1.0 - V2_ALPHA)
    with np.errstate(divide="ignore", invalid="ignore"):
        zlo, zhi = (mean + tol) / sem, (tol - mean) / sem
    ok = (sem > 0) & (zlo >= z) & (zhi >= z)
    return {"mean": mean, "sem": sem, "z_lower": zlo, "z_upper": zhi, "ok": ok}


def step_v2_combine(a: argparse.Namespace) -> int:
    out: dict[str, Any] = {}
    ok, used, reduced = True, 0, False
    for e in V2_ENERGIES:
        row = f"v2-{e:g}"
        parts = load_partials(a, [f"{row}-s{k}.json" for k in range(V2_SHARDS)])
        for k, p in enumerate(parts):
            if p["seed"] != seed_of(row, k) or p["shard"] != k or not p["valid"]:
                raise SystemExit(f"{row} shard {k}: seed/shard mismatch or invalid result")
        ratio = np.concatenate([np.array(p["ratio"]) for p in parts])
        n = sum(p["n"] for p in parts)
        used += n
        reduced |= any(p["reduced"] for p in parts)
        geo, grid, nz, _ = depth_box(M.WATER, e)
        sim = Simulation(nuc_config(energy=e, n=2000, seed=0, geometry=geo, grid=grid, mcs=False,
                                    straggling=False, tallies=PRIMARY))  # fmt: skip
        ref = reference_survival(sim, e, nz, 1.0)
        bins = depth_bins(e, nz)
        t = tost(ratio[:, bins] - ref[bins], V2_TOL)
        row_ok = bool(t["ok"].all())
        ok &= row_ok
        out[row] = {
            "histories": n, "batches": int(ratio.shape[0]), "depths_mm": [float(b) for b in bins],
            "s_measured": ratio.mean(axis=0)[bins].tolist(), "s_reference": ref[bins].tolist(),
            "diff": t["mean"].tolist(), "sem": t["sem"].tolist(),
            "z_lower_min": float(t["z_lower"].min()), "z_upper_min": float(t["z_upper"].min()),
            "max_abs_diff": float(np.abs(t["mean"]).max()), "tolerance": V2_TOL,
            "alpha": V2_ALPHA, "pass": row_ok,
        }  # fmt: skip
    doc = {"step": "v2-combine", "table": table_record(), "v2": out, "pass": bool(ok)}
    return finish5(doc, 3 * V2_SHARDS * V2_SHARD_N, used, reduced)


# -- V2-probe -------------------------------------------------------------------------------------
def step_probe_shard(a: argparse.Namespace) -> int:
    name, k = a.probe, a.shard
    shards, n_full, var_a, var_b = PROBES[name]
    n = scaled_n(n_full, a.scale)
    row = f"v2-probe-{name}"
    geo, grid, nz, _ = depth_box(M.WATER, 150.0)
    ratios, okc, wall = [], True, 0.0
    for smax, frac in (var_a, var_b):
        cfg = nuc_config(energy=150.0, n=n, seed=seed_of(row, k), geometry=geo, grid=grid, mcs=False,
                         straggling=False, max_step=smax, frac=frac, tallies=PRIMARY,
                         timeout=a.timeout)  # fmt: skip
        t0 = time.perf_counter()
        res = Simulation(cfg).run()
        wall += time.perf_counter() - t0
        okc &= clean(res)
        ratios.append(primary_ratio(res))
    path = write_partial(a, f"{row}-s{k}", {
        "row": row, "shard": k, "seed": seed_of(row, k), "n": n, "valid": okc,
        "diff": (ratios[0] - ratios[1]).tolist(), "s_b": ratios[1].tolist(),
        "reduced": n < n_full,
    })  # fmt: skip
    doc = {"step": "v2-probe-shard", "row": row, "shard": k, "seed": seed_of(row, k),
           "table": table_record(), "wall_s": wall, "hist_per_s": 2 * n / wall, "partial": path,
           "pass": bool(okc)}  # fmt: skip
    return finish5(doc, n_full, n, n < n_full)


def step_probe_combine(a: argparse.Namespace) -> int:
    out: dict[str, Any] = {}
    ok, used, reduced, frozen = True, 0, False, 0
    for name, (shards, n_full, var_a, var_b) in PROBES.items():
        row = f"v2-probe-{name}"
        parts = load_partials(a, [f"{row}-s{k}.json" for k in range(shards)])
        for k, p in enumerate(parts):
            if p["seed"] != seed_of(row, k) or p["shard"] != k or not p["valid"]:
                raise SystemExit(f"{row} shard {k}: seed/shard mismatch or invalid result")
        diff = np.concatenate([np.array(p["diff"]) for p in parts])
        n = sum(p["n"] for p in parts)
        used += n
        frozen += shards * n_full
        reduced |= any(p["reduced"] for p in parts)
        nz = diff.shape[1]
        bins = depth_bins(150.0, nz)
        mean = diff.mean(axis=0)[bins]
        sem = diff.std(axis=0, ddof=1)[bins] / math.sqrt(diff.shape[0])
        conclusive = bool(sem[-1] <= PROBE_SIGMA_MAX)
        row_ok = bool(conclusive and np.all(np.abs(mean) <= PROBE_TOL))
        ok &= row_ok
        out[row] = {
            "variants": {"a": {"s_max_mm": var_a[0], "f_e": var_a[1]},
                         "b": {"s_max_mm": var_b[0], "f_e": var_b[1]}},
            "histories": n, "batches": int(diff.shape[0]), "depths_mm": [float(b) for b in bins],
            "diff": mean.tolist(), "sem": sem.tolist(), "max_abs_diff": float(np.abs(mean).max()),
            "sigma_delta_deepest": float(sem[-1]), "sigma_max": PROBE_SIGMA_MAX,
            "conclusive": conclusive, "tolerance": PROBE_TOL, "pass": row_ok,
        }  # fmt: skip
    doc = {"step": "v2-probe-combine", "table": table_record(), "probes": out, "pass": bool(ok)}
    return finish5(doc, frozen, used, reduced)


# -- V3-LV ----------------------------------------------------------------------------------------
BINDING_RTOL, BINDING_ATOL_MEV = 1e-9, 1e-9


def binding_recompute(res: Result, tab_info: dict[str, Any], ame: Any) -> dict[str, Any]:
    """Independent check of the ``nuclear_binding`` tally: the binding of the run recomputed from
    the per-target event, light-product and residual counts of the reference (``res.diagnostics
    ["nuclear"]``) and AME2020 masses read from the cached AME text (not from the table arrays);
    ``|recomputed - tally| <= 1e-9 |tally| + 1e-9 MeV``."""
    tally = float(res.energy_balance.nuclear_mev["nuclear_binding"])
    nuc = res.diagnostics.get("nuclear", {})
    got = nc.recompute_binding_total(nuc, tab_info["targets"], ame)
    diff = abs(got - tally)
    return {"tally_mev": tally, "recomputed_mev": got, "abs_diff_mev": diff,
            "events": int(sum(r["events"] for r in nuc.values())),
            "pass": bool(diff <= BINDING_RTOL * abs(tally) + BINDING_ATOL_MEV)}  # fmt: skip


def step_v3_lv(a: argparse.Namespace) -> int:
    from ionmc.data.ame import load_ame2020

    cdir = cache.resolve_cache_dir(None)
    ame = load_ame2020(cache.verify("ame2020-mass", cdir).read_text(encoding="ascii"))
    tab_info = NuclearTable.load(cdir, TABLE_ID).info
    out: dict[str, Any] = {}
    ok, used, frozen, reduced, j = True, 0, 0, False, 0
    for mname, mat in V3_MATERIALS:
        for e, n_full in V3_HISTORIES.items():
            n = v4.scaled(n_full, a.scale, 2000, 2)
            geo, grid, _, _ = depth_box(mat, e)
            cfg = nuc_config(energy=e, n=n, seed=seed_of("v3-lv", j), geometry=geo, grid=grid,
                             n_batches=2, timeout=a.timeout)  # fmt: skip
            t0 = time.perf_counter()
            res = Simulation(cfg).run()
            wall = time.perf_counter() - t0
            b = res.energy_balance
            rel, grel = float(b.relative_residual), float(b.grid_relative_residual(0))
            bind = binding_recompute(res, tab_info, ame)
            row_ok = bool(clean(res) and rel <= 1e-12 and grel <= 1e-12 and bind["pass"])
            ok &= row_ok
            used, frozen, reduced = used + n, frozen + n_full, reduced or n < n_full
            out[f"{mname}@{e:g}"] = {
                "seed": seed_of("v3-lv", j), "histories": n, "wall_s": wall,
                "relative_residual": rel, "grid_relative_residual": grel,
                "initial_mev": float(b.initial_mev), "nuclear_mev": dict(b.nuclear_mev),
                "nuclear_imbalance_over_initial": float(b.nuclear_mev.get("nuclear_imbalance", 0.0))
                / float(b.initial_mev),
                "truncated_mev": float(b.truncated_mev), "unaccounted_mev": float(b.unaccounted_mev),
                "binding_recompute": bind,
                "valid": bool(res.valid), "counters": res.counters.as_dict(), "pass": row_ok,
            }  # fmt: skip
            j += 1
    doc = {"step": "v3-lv", "table": table_record(), "runs": out, "tolerance": 1e-12,
           "binding_recompute": {"rtol": BINDING_RTOL, "atol_mev": BINDING_ATOL_MEV,
                                 "evaluated": True, "gating": True},
           "pass": bool(ok)}  # fmt: skip
    return finish5(doc, frozen, used, reduced)


# -- V4, V4b --------------------------------------------------------------------------------------
def step_v4(a: argparse.Namespace) -> int:
    nc.V4_EVENTS = v4.scaled(V4_EVENTS, a.scale, 2000)
    nc.V4_BASE_SEED = seed_of("v4")
    cdir = cache.resolve_cache_dir(None)
    tab = NuclearTable.load(cdir, TABLE_ID)
    res = nc.v4_checks(tab, nc._targets(cdir))
    cases = {}
    for k, c in res["cases"].items():
        cases[k] = {
            "yield_ratio": c["yield_ratio"], "yield_ratio_3sem": c["yield_ratio_3sem"],
            "sum_y_e_prime_ratio": c["sum_y_e_prime_ratio"], "mean_e_prime_ratio": c["mean_e_prime_ratio"],
            "mean_e_prime_ratio_sem": c["mean_e_prime_ratio_sem"],
            "delta_lab_over_e_avail": c["delta_lab_over_e_avail"], "p_accept_exact": c["p_accept_exact"],
            "max_ledger_closure_mev": c["max_ledger_closure_mev"], "v4_yields_pass": c["v4_yields_pass"],
            "v4_sum_pass": c["v4_sum_pass"], "v4b_pass": c["v4b_pass"],
        }  # fmt: skip
    ok = bool(res["V4_pass"] and res["V4b_pass"])
    doc = {"step": "v4-v4b", "table": table_record(), "events_per_case": nc.V4_EVENTS,
           "seed": seed_of("v4"), "V4_pass": res["V4_pass"], "V4b_pass": res["V4b_pass"],
           "cases": cases, "full_cases": res["cases"], "pass": ok}  # fmt: skip
    return finish5(doc, V4_EVENTS, nc.V4_EVENTS, nc.V4_EVENTS < V4_EVENTS)


# -- X1 -------------------------------------------------------------------------------------------
def step_x1(a: argparse.Namespace) -> int:
    n = v4.scaled(X1_N, a.scale, 2000, X1_BATCHES)
    geo, grid, nz, _ = depth_box(M.WATER, 150.0)
    runs, wall = {}, {}
    for tag, nuclear in (("on", True), ("off", False)):
        cfg = nuc_config(energy=150.0, n=n, seed=seed_of("x1"), geometry=geo, grid=grid,
                         nuclear=nuclear, n_batches=X1_BATCHES, tallies=EDEP, timeout=a.timeout)  # fmt: skip
        t0 = time.perf_counter()
        runs[tag] = Simulation(cfg).run()
        wall[tag] = time.perf_counter() - t0
    z = (np.arange(nz) + 0.5) * 1.0
    plateau = (z >= X1_PLATEAU_MM[0]) & (z <= X1_PLATEAU_MM[1])
    idd = {t: edep_batches(r) for t, r in runs.items()}
    peak = int(np.argmax(idd["off"].mean(axis=0)))
    ratio = {t: v[:, peak] / v[:, plateau].mean(axis=1) for t, v in idd.items()}
    d = ratio["on"] - ratio["off"]
    sem = float(d.std(ddof=1) / math.sqrt(d.size))
    zst = float(-d.mean() / sem) if sem > 0 else math.nan
    ok = bool(all(clean(r) for r in runs.values()) and d.mean() < 0 and zst > X1_Z)
    doc = {
        "step": "x1", "table": table_record(), "seed": seed_of("x1"), "peak_bin": peak,
        "plateau_mm": list(X1_PLATEAU_MM), "ratio_on": float(ratio["on"].mean()),
        "ratio_off": float(ratio["off"].mean()), "paired_diff_mean": float(d.mean()),
        "paired_diff_sem": sem, "paired_z": zst, "z_required": X1_Z, "wall_s": wall,
        "counters_on": runs["on"].counters.as_dict(), "valid": {t: bool(r.valid) for t, r in runs.items()},
        "pass": ok,
    }  # fmt: skip
    return finish5(doc, 2 * X1_N, 2 * n, n < X1_N)


# -- E1 (exploratory, report only) ----------------------------------------------------------------
def step_e1(a: argparse.Namespace) -> int:
    from ionmc.reference import metrics as rm

    n = v4.scaled(E1_N, a.scale, 2000, E1_BATCHES)
    dz = 0.5
    geo = BoxPhantom((-60.0, -60.0, 0.0), (120.0, 120.0, 300.0), M.WATER)
    grid = ScoringGrid((-60.0, -60.0, 0.0), (120.0, 120.0, dz), (1, 1, int(300.0 / dz)), name="dose")
    cfg = nuc_config(energy=150.0, n=n, seed=seed_of("e1"), geometry=geo, grid=grid,
                     n_batches=E1_BATCHES, timeout=a.timeout)  # fmt: skip
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    idd = res.grid("dose").energy_mev.reshape(-1) / n
    depth = (np.arange(idd.size) + 0.5) * dz
    curve = rm.normalize_to_peak(idd)
    met = {"peak_depth_mm": rm.peak_depth(depth, curve), "r80_mm": rm.r80(depth, curve),
           "r90_mm": rm.r90(depth, curve),
           "falloff_80_20_mm": rm.distal_falloff_80_20(depth, curve),
           "plateau_over_peak": float(curve[(depth > 20) & (depth < 60)].mean())}  # fmt: skip
    doc = {
        "step": "e1", "status": "exploratory: report only; compare with the V3-010B TOPAS and "
        "MCsquare runs (compare_depth_dose.py, materialized by the orchestrator)",
        "table": table_record(), "seed": seed_of("e1"), "bin_width_mm": dz, "wall_s": wall,
        "metrics": met, "idd_mev_per_history": idd.tolist(), "valid": clean(res),
        "nuclear_mev": dict(res.energy_balance.nuclear_mev), "pass": clean(res),
    }  # fmt: skip
    return finish5(doc, E1_N, n, n < E1_N)


# -- R1 -------------------------------------------------------------------------------------------
def _captured(fn: Any, ns: argparse.Namespace) -> dict[str, Any]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        fn(ns)
    text = buf.getvalue()
    return json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])


def step_r1(a: argparse.Namespace) -> int:
    import run_suite

    a16 = _captured(v4.step_a16, argparse.Namespace(mode="regression", intended_change_record=None))
    k = max(4, int(256 * a.scale))
    t1 = _captured(base.step_t1, argparse.Namespace(k=k, energy=150.0, seed=base.SEED_BASE,
                                                    workers=1, timeout=a.timeout))  # fmt: skip
    digest = a16.get("baseline_digests_sha256")
    doc = {
        "step": "r1", "a16_baseline": a16["baseline_ref"], "a16_baseline_expected": A16_R1_BASELINE,
        "intended_change_set": run_suite.A16_INTENDED_CHANGE, "digest_seed": 20351004,
        "digest_sha256": digest, "digest_recorded": A16_R1_DIGEST,
        "digest_matches_recorded": digest == A16_R1_DIGEST,
        "digest_gate": "informational: sha256 of the sorted-key JSON of the baseline digests; the "
        "gating comparison is the field-by-field identity of the tree under test with the baseline",
        "a16": {"pass": a16["pass"], "mode": a16["mode"],
                "specs": {s: {"tally_neutral": v["tally_neutrality"]["identical"],
                              "identical_to_baseline": v["no_tallies_vs_baseline"]["identical"]}
                          for s, v in a16["specs"].items()}},
        "t1": {"pass": t1["pass"], "k": t1["histories"], "verdict": t1["t1"].get("pass")},
    }  # fmt: skip
    doc["pass"] = bool(
        a16["pass"] and t1["pass"] and a16["baseline_ref"] == A16_R1_BASELINE
        and run_suite.A16_INTENDED_CHANGE is None
    )
    return finish5(doc, 256, k, k < 256)


# -- deferred: 1 versus N workers ----------------------------------------------------------------
def step_partition(a: argparse.Namespace) -> int:
    """1-vs-N worker partition invariance of nuclear runs (deferred in the single-process mode)."""
    if os.environ.get("IONMC_SINGLE_PROCESS") == "1":
        raise SystemExit("v3-workers-partition needs worker processes: deferred (IONMC_SINGLE_PROCESS=1)")
    geo, grid, _, _ = depth_box(M.WATER, 150.0)
    res = {}
    for w in (1, PARTITION_WORKERS):
        cfg = nuc_config(energy=150.0, n=PARTITION_N, seed=seed_of("v3-lv", 7), geometry=geo,
                         grid=grid, n_batches=6, workers=w, timeout=a.timeout)  # fmt: skip
        res[w] = Simulation(cfg).run()
    same = bool(np.array_equal(res[1].grid("dose").batch_energy_mev,
                               res[PARTITION_WORKERS].grid("dose").batch_energy_mev)
                and res[1].counters.as_dict() == res[PARTITION_WORKERS].counters.as_dict()
                and dict(res[1].energy_balance.nuclear_mev)
                == dict(res[PARTITION_WORKERS].energy_balance.nuclear_mev))  # fmt: skip
    doc = {"step": "v3-workers", "table": table_record(), "workers": [1, PARTITION_WORKERS],
           "identical": same, "pass": bool(same and all(clean(r) for r in res.values()))}  # fmt: skip
    return finish5(doc, PARTITION_N, PARTITION_N, False)


STEPS = {
    "lv5-throughput": step_throughput,
    "n1": step_n1,
    "v2-shard": step_v2_shard,
    "v2-combine": step_v2_combine,
    "v2-probe-shard": step_probe_shard,
    "v2-probe-combine": step_probe_combine,
    "v3-lv": step_v3_lv,
    "v4-v4b": step_v4,
    "x1": step_x1,
    "e1": step_e1,
    "r1": step_r1,
    "v3-workers": step_partition,
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--workers", default="1", help="accepted for uniformity; lv5 runs one process")
    ap.add_argument("--energy", type=float, default=150.0, help="v2-shard: 100, 150 or 200")
    ap.add_argument("--shard", type=int, default=0, help="v2-shard, v2-probe-shard: shard index")
    ap.add_argument("--probe", choices=sorted(PROBES), default="s05", help="v2-probe-shard")
    ap.add_argument("--out-dir", default="samples", help="partial files of the shard steps")
    ap.add_argument("--dirs", nargs="+", default=["."], help="combine steps: archive directories")
    ap.add_argument("--partials-manifest", default=None, help="combine steps: JSON name -> "
                    "content_sha256 of the imported partials (required with imported --dirs)")
    ap.add_argument("--seed-base", "--seed", dest="seed", type=int, default=QUALIFICATION_SEED_BASE,
                    help="lv5 base (qualification 20421004; rehearsals 2043xxxx)")  # fmt: skip
    ap.add_argument("--timeout", type=float, default=None)
    args = ap.parse_args(argv)
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    if args.step == "v2-shard" and (args.energy not in V2_ENERGIES or not 0 <= args.shard < V2_SHARDS):
        raise SystemExit("v2-shard: --energy in 100/150/200, --shard in 0..2")
    if args.step == "v2-probe-shard" and not 0 <= args.shard < PROBES[args.probe][0]:
        raise SystemExit("v2-probe-shard: --shard outside the shards of the probe")
    base.SEED_BASE = args.seed
    v4.base.SEED_BASE = args.seed
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
