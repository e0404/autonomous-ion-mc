"""Validation steps of the V3-005 acceptance plan, slice C (suite ``lv5c``; decision 0041 slice C).

Usage::

    python validation/scripts/transport/steps_v5c.py <step> [--scale F] [--seed-base N] ...

Configuration builders of the elastic-on suites (C6 brief section 1) and the steps of row V7-R
(plan Amendment 14 (h), Amendment 17 (a)): ``v7r-shard`` (one of 16 shards of 4.5e6 histories with the
block-sum sidecar), ``v7r-ref``, ``v7r-combine`` (the gates of ``v7r.py``) and ``v7r-diag``
(exploratory); row V11 (plan Amendment 14 (c)): ``v11-ionmc`` (EM+elastic and EM-only producers,
r_index 14, common random numbers per energy) and ``v11-compare``. The frozen steps of ``steps_v5.py`` / ``steps_v5b.py`` are imported, never edited, and
build their configurations with ``elastic=False``; this module re-enables the elastic channel
explicitly with the table id pinned in ``src/ionmc/data/elastic_table_pin.json``.

Block-sum sidecar (``v7r-s{k}-escaped_neutral.npz``)
----------------------------------------------------
Written next to the JSON partial ``v7r-s{k}.json`` by ``np.savez`` (uncompressed; numpy fixes the zip
timestamps, so equal arrays give equal bytes). It holds the one array ``escaped_neutral_block_sums``,
shape ``[450, 20]`` (replicate, block), float64, C order, little endian: the true SUM of the escaped
neutral energy (neutron + gamma, ``nuclear_escaped_neutron + nuclear_escaped_gamma``) over the 500
histories of each block in MeV (``units: "MeV per 500-history block"``, Amendment 17 (a)1). The JSON
partial carries ``sidecar: {file, sha256 (of the file bytes), array_sha256 (of the C-order
little-endian float64 bytes), array, shape, dtype, units, batch_histories, blocks_per_replicate}``,
so its ``content_sha256`` binds the sidecar. The JSON replicate ``mean`` is the per-primary value
``mean(sums) / 500`` (asserted to 1e-14 relative at write and at load).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

import steps as base
import steps_v4 as v4
import steps_v5 as v5
import steps_v5b as v5b
import v7r
from ionmc.config import SimulationConfig
from ionmc.data import cache

FORMAT = 1
LV5C_SEED_BASE = 20481004  # lv5c qualification base (Amendment 14 (j))
ELASTIC_TABLE_ID = "1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591"
"""The derived elastic table of the pin ``src/ionmc/data/elastic_table_pin.json`` (a test asserts it)."""
V7R_SHARDS = v7r.V7R_SHARDS  # 16; run_suite.V7R_SHARDS (host interpreter has no ionmc) is asserted equal
V7R_ESTIMATORS = (("sec_p", "profile"), ("nuclear_local", "profile"), ("escaped_neutral", "scalar"))
SIDECAR_ARRAY = "escaped_neutral_block_sums"
SIDECAR_ESTIMATOR = "escaped_neutral"
CONSUMED_SEED_BASES_C = (20471004, 20451004, 20441004)  # lv5b, hr5, b40d8121 (evidence, consumed)
CONSUMED_R_INDICES = range(0, 13)  # r_index 0..12 of the consumed bases (lv5b uses 1..11, hr5 1..3)
CONSUMED_SHARDS = range(0, 33)


# -- configuration builders (C6 brief section 1) ------------------------------------------------
def nuc_config_c(**kw: Any) -> SimulationConfig:
    """``steps_v5.nuc_config`` with the hadronic elastic channel ON and the pinned table id. Do not
    use ``steps_v5.nuc_config`` directly for lv5c/hr5c: it is pinned ``elastic=False`` (frozen suites)."""
    cfg = v5.nuc_config(**kw)
    return replace(
        cfg,
        physics=replace(
            cfg.physics, elastic=True, elastic_only=False, elastic_table_id=ELASTIC_TABLE_ID
        ),
    )


def nuc_config_c_elastic_only(**kw: Any) -> SimulationConfig:
    """As ``nuc_config_c`` with the non-elastic channel off (EM + elastic only, row V11)."""
    cfg = nuc_config_c(**kw)
    return replace(cfg, physics=replace(cfg.physics, elastic_only=True))


def wcfg_c(backend: str, precision: str, **kw: Any) -> SimulationConfig:
    """``nuc_config_c`` re-targeted to ``backend`` / ``precision`` (as ``steps_v5b.wcfg``)."""
    cfg = nuc_config_c(**kw)
    return replace(cfg, run=replace(cfg.run, backend=backend, precision=precision))


def elastic_table_record() -> dict[str, Any]:
    """Id and npz hash of the pinned elastic table, read from the cache files."""
    p = cache.resolve_cache_dir(None) / "derived" / f"elastic-proton-{ELASTIC_TABLE_ID}.npz"
    return {"elastic_table_id": ELASTIC_TABLE_ID,
            "elastic_table_npz_sha256": hashlib.sha256(p.read_bytes()).hexdigest()}  # fmt: skip


def bindings_c(a: argparse.Namespace) -> dict[str, Any]:
    """``steps_v5.bindings`` plus the elastic table id and hash and ``elastic: true``."""
    return {**v5.bindings(a), **elastic_table_record(), "elastic": True}


def write_partial_c(a: argparse.Namespace, name: str, doc: dict[str, Any]) -> str:
    """``steps_v5.write_partial`` with the lv5c bindings inside the sealed document."""
    return v5.write_partial(a, name, {**elastic_table_record(), "elastic": True, **doc})


def _find_partial(a: argparse.Namespace, name: str) -> Path:
    hits = [p for d in a.dirs for p in (Path(d) / name, Path(d) / "samples" / name) if p.is_file()]
    if len(hits) != 1:
        raise SystemExit(f"partial {name}: found {len(hits)} copies in {a.dirs} (need exactly 1)")
    return hits[0]


def load_partials_c(a: argparse.Namespace, names: list[str]) -> list[tuple[dict[str, Any], Path]]:
    """``steps_v5.load_partials`` (digest, manifest and run bindings) plus the lv5c bindings, returning
    ``(document, path)`` pairs (the sidecar must be in the directory of its JSON hit)."""
    want = bindings_c(a)
    out = []
    for nm in names:
        (doc,) = v5.load_partials(a, [nm])
        for key in ("elastic", "elastic_table_id", "elastic_table_npz_sha256"):
            if doc.get(key) != want[key]:
                raise SystemExit(f"partial {nm}: {key} is {doc.get(key)!r}, expected {want[key]!r}")
        out.append((doc, _find_partial(a, nm)))
    return out


# -- block-sum sidecar ---------------------------------------------------------------------------
def array_sha256(arr: NDArray[np.float64]) -> str:
    """sha256 of the C-order little-endian float64 bytes of ``arr``."""
    return hashlib.sha256(np.ascontiguousarray(arr, dtype="<f8").tobytes()).hexdigest()


def write_v7r_sidecar(path: Path, sums: NDArray[np.float64]) -> dict[str, Any]:
    """Write the sidecar ``path`` (the array ``escaped_neutral_block_sums`` ``[replicates, 20]``) and
    return its ``sidecar`` record for the JSON partial."""
    arr = np.ascontiguousarray(sums, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != v7r.V7R_BLOCKS_PER_REP or not np.all(np.isfinite(arr)):
        raise SystemExit("v7r sidecar: need a finite [replicates, 20] float64 array")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **{SIDECAR_ARRAY: arr})
    return {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "array_sha256": array_sha256(arr), "array": SIDECAR_ARRAY, "shape": list(arr.shape),
            "dtype": "float64", "units": v7r.V7R_UNITS, "batch_histories": v7r.V7R_BATCH_N,
            "blocks_per_replicate": v7r.V7R_BLOCKS_PER_REP}  # fmt: skip


def load_v7r_sidecar(doc: dict[str, Any], partial_path: Path) -> NDArray[np.float64]:
    """The block sums bound by ``doc["sidecar"]``: file in the directory of ``partial_path``, file and
    array sha256, shape, dtype, units and finiteness verified, and the JSON replicate mean equal to
    ``mean(sums) / batch_histories`` (1e-14 relative); fails closed (``SystemExit``)."""
    rec = doc.get("sidecar")
    if not isinstance(rec, dict):
        raise SystemExit(f"{partial_path.name}: no sidecar record")
    need = ("file", "sha256", "array_sha256", "array", "shape", "dtype", "units", "batch_histories",
            "blocks_per_replicate")  # fmt: skip
    if any(k not in rec for k in need):
        raise SystemExit(f"{partial_path.name}: incomplete sidecar record")
    if rec["file"] != Path(rec["file"]).name:
        raise SystemExit(f"{partial_path.name}: sidecar file must be a plain file name")
    p = partial_path.parent / rec["file"]
    if not p.is_file():
        raise SystemExit(f"{partial_path.name}: sidecar {rec['file']} not found next to the partial")
    raw = p.read_bytes()
    if hashlib.sha256(raw).hexdigest() != rec["sha256"]:
        raise SystemExit(f"{rec['file']}: file sha256 differs from the partial")
    try:
        import io

        with np.load(io.BytesIO(raw), allow_pickle=False) as z:
            if list(z.files) != [rec["array"]] or rec["array"] != SIDECAR_ARRAY:
                raise SystemExit(f"{rec['file']}: unexpected array set {list(z.files)}")
            arr = z[SIDECAR_ARRAY]
    except SystemExit:
        raise
    except Exception as exc:  # corrupt npz, pickled object arrays, ...
        raise SystemExit(f"{rec['file']}: unreadable sidecar: {exc!r}") from exc
    reps = int(doc.get("replicates", -1))
    if arr.dtype != np.dtype("float64"):
        raise SystemExit(f"{rec['file']}: dtype {arr.dtype}, expected float64")
    if (list(arr.shape) != list(rec["shape"]) or arr.shape != (reps, v7r.V7R_BLOCKS_PER_REP)
            or rec["dtype"] != "float64" or rec["units"] != v7r.V7R_UNITS
            or rec["batch_histories"] != v7r.V7R_BATCH_N
            or doc.get("batch_histories") != v7r.V7R_BATCH_N
            or rec["blocks_per_replicate"] != v7r.V7R_BLOCKS_PER_REP):  # fmt: skip
        raise SystemExit(f"{rec['file']}: shape/dtype/units differ from the partial")
    if not np.all(np.isfinite(arr)):
        raise SystemExit(f"{rec['file']}: non-finite block sums")
    if array_sha256(arr) != rec["array_sha256"]:
        raise SystemExit(f"{rec['file']}: array sha256 differs from the partial")
    est = doc.get("estimators", {}).get(SIDECAR_ESTIMATOR)
    if est is None:
        raise SystemExit(f"{partial_path.name}: no {SIDECAR_ESTIMATOR} estimator")
    mean = np.array(est["mean"], float).reshape(-1)
    sem = np.array(est["sem"], float).reshape(-1)
    per_primary = arr / v7r.V7R_BATCH_N
    ok = (np.allclose(per_primary.mean(axis=1), mean, rtol=1e-14, atol=0.0)
          and np.allclose(per_primary.std(axis=1, ddof=1) / math.sqrt(arr.shape[1]), sem,
                          rtol=1e-12, atol=0.0))  # fmt: skip
    if not ok:
        raise SystemExit(f"{rec['file']}: JSON replicate mean/sem inconsistent with the block sums")
    return arr.astype(np.float64)


# -- shard and reference producers --------------------------------------------------------------
def _seed(k: int) -> int:
    return base.SEED_BASE + 1000 * v7r.V7R_R_INDEX + k


def v7r_config(n: int, nb: int, seed: int, *, timeout: float | None = None) -> SimulationConfig:
    """warp-cpu float64, 150 MeV, 12-bin coarse depth geometry, elastic on (``nuc_config_c``)."""
    geo, grid = v5b.coarse_depth(v5b.V7_BINS)
    return wcfg_c("warp-cpu", "float64", energy=150.0, n=n, seed=seed, geometry=geo, grid=grid,
                  n_batches=nb, tallies=v5b.V7_TALLIES, timeout=timeout)  # fmt: skip


def counter_names() -> tuple[str, ...]:
    """Names of the counter vector of an elastic-on run, in the order of ``merge_partials``."""
    from ionmc.transport.tally import COUNTER_NAMES, ELASTIC_COUNTER_NAMES, NUCLEAR_COUNTER_NAMES

    return COUNTER_NAMES + NUCLEAR_COUNTER_NAMES + ELASTIC_COUNTER_NAMES


def counters_by_name(est: dict[str, Any]) -> dict[str, int]:
    """Raw per-counter sums of ``batch_estimates`` output, by name (recorded in the partial)."""
    vec = np.asarray(est["counter_vector"], dtype=np.int64)
    names = counter_names()
    if vec.size != len(names):
        raise SystemExit(f"v7r: counter vector of {vec.size} entries, expected {len(names)}")
    return {nm: int(v) for nm, v in zip(names, vec, strict=True)}


def invalid_counters(est: dict[str, Any]) -> int:
    """Sum of the validity counters: every counter BY NAME except the diagnostics that
    ``ionmc.simulation.ElasticTransportCounters.any_nonzero`` ignores (``elastic_below_domain``,
    ``pp_below_domain``; Amendment 15 (b)2: they never invalidate a result)."""
    from ionmc.simulation import ElasticTransportCounters

    ignore = set(ElasticTransportCounters.DIAGNOSTIC)
    return sum(v for k, v in counters_by_name(est).items() if k not in ignore)


def step_v7r_shard(a: argparse.Namespace) -> int:
    k = a.shard
    n = v4.scaled(v7r.V7R_SHARD_N, a.scale, 2 * v5b.V7_REP_N, v5b.V7_REP_N)
    nb = n // v7r.V7R_BATCH_N
    t0 = time.perf_counter()
    est = v5b.batch_estimates(v7r_config(n, nb, _seed(k), timeout=a.timeout))
    wall = time.perf_counter() - t0
    ok = bool(invalid_counters(est) == 0 and est["batch_assignment_ok"])
    stats = {nm: v5b.v7_replicate_stats(x) for nm, x in v5b.v7_batch_arrays(est).items()}
    reps = nb // v7r.V7R_BLOCKS_PER_REP
    sums = np.asarray(est["block_sums"], float).reshape(reps, v7r.V7R_BLOCKS_PER_REP)
    side = write_v7r_sidecar(Path(a.out_dir) / f"v7r-s{k}-{SIDECAR_ESTIMATOR}.npz", sums)
    path = write_partial_c(a, f"v7r-s{k}", {
        "row": "v7-r", "kind": "shard", "shard": k, "seed": _seed(k), "n": n, "n_batches": nb,
        "batch_histories": v7r.V7R_BATCH_N, "blocks_per_replicate": v7r.V7R_BLOCKS_PER_REP,
        "replicates": reps, "bins": v5b.V7_BINS, "valid": ok, "counters_sum": est["counters_sum"], "counters_invalid": invalid_counters(est), "counters_by_name": counters_by_name(est),
        "wall_s": wall, "estimators": {nm: {"mean": m.tolist(), "sem": s.tolist()}
                                       for nm, (m, s) in stats.items()},
        "sidecar": side, "reduced": n < v7r.V7R_SHARD_N})  # fmt: skip
    load_ok = load_v7r_sidecar(json.loads(Path(path).read_text()), Path(path)).shape == sums.shape
    doc = {"step": "v7r-shard", "shard": k, "seed": _seed(k), "table": v5.table_record(),
           **elastic_table_record(), "histories": n, "batches": nb, "replicates": reps,
           "wall_s": wall, "hist_per_s": n / wall, "counters_sum": est["counters_sum"], "counters_invalid": invalid_counters(est), "counters_by_name": counters_by_name(est),
           "partial": path, "sidecar": side, "pass": bool(ok and load_ok)}  # fmt: skip
    return v5b.finish5b(doc, v7r.V7R_SHARD_N, n, n < v7r.V7R_SHARD_N)


def step_v7r_ref(a: argparse.Namespace) -> int:
    n = v4.scaled(v7r.V7R_REF_N, a.scale, 20_000, v7r.V7R_REF_BLOCKS)
    seed = _seed(v7r.V7R_REF_INDEX)
    t0 = time.perf_counter()
    est = v5b.batch_estimates(v7r_config(n, v7r.V7R_REF_BLOCKS, seed, timeout=a.timeout))
    wall = time.perf_counter() - t0
    ok = bool(invalid_counters(est) == 0 and est["batch_assignment_ok"])
    ref = {}
    for nm, x in v5b.v7_batch_arrays(est).items():
        x = np.asarray(x, float).reshape(v7r.V7R_REF_BLOCKS, -1)
        ref[nm] = {"mean": x.mean(axis=0).tolist(),
                   "sem": (x.std(axis=0, ddof=1) / math.sqrt(v7r.V7R_REF_BLOCKS)).tolist(),
                   "blocks": x.tolist()}  # fmt: skip
    path = write_partial_c(a, "v7r-ref", {
        "row": "v7-r", "kind": "reference", "seed": seed, "n": n, "n_batches": v7r.V7R_REF_BLOCKS,
        "bins": v5b.V7_BINS, "valid": ok, "counters_sum": est["counters_sum"], "counters_invalid": invalid_counters(est), "counters_by_name": counters_by_name(est), "wall_s": wall,
        "estimators": ref, "reduced": n < v7r.V7R_REF_N})  # fmt: skip
    doc = {"step": "v7r-ref", "seed": seed, "table": v5.table_record(), **elastic_table_record(),
           "histories": n, "batches": v7r.V7R_REF_BLOCKS, "wall_s": wall, "hist_per_s": n / wall,
           "counters_sum": est["counters_sum"], "counters_invalid": invalid_counters(est), "counters_by_name": counters_by_name(est), "partial": path, "pass": ok}  # fmt: skip
    return v5b.finish5b(doc, v7r.V7R_REF_N, n, n < v7r.V7R_REF_N)


# -- combine -------------------------------------------------------------------------------------
def consumed_seed_set() -> set[int]:
    """Seeds of the consumed evidence bases (``base + 1000 r + k`` over the r_index and shard ranges
    ever used by lv5b / hr5 / b40d8121); every combine step checks disjointness (plan l.380)."""
    return {b + 1000 * r + k for b in CONSUMED_SEED_BASES_C for r in CONSUMED_R_INDICES
            for k in CONSUMED_SHARDS}  # fmt: skip


def v7r_seed_list(seed_base: int, n_shards: int = V7R_SHARDS) -> dict[str, list[int]]:
    return {"shards": [seed_base + 1000 * v7r.V7R_R_INDEX + k for k in range(n_shards)],
            "reference": [seed_base + 1000 * v7r.V7R_R_INDEX + v7r.V7R_REF_INDEX]}  # fmt: skip


def v7r_combine(
    shards: list[tuple[dict[str, Any], Path]], ref: dict[str, Any], seed_base: int, *,
    n_shards: int = V7R_SHARDS, n_boot: int = v5b.V7_BOOT_N,
) -> dict[str, Any]:
    """Document of row V7-R from the 16 shard partials (with sidecars) and the reference partial.
    Pure given its inputs; ``n_shards`` and ``n_boot`` exist for the CI-scale dry run only (the step
    uses the defaults). Fails closed (``SystemExit``) on a missing/invalid item, unequal replicate
    counts, a seed mismatch or a seed set intersecting a consumed evidence seed set."""
    if len(shards) != n_shards:
        raise SystemExit(f"v7r: need {n_shards} shard partials, got {len(shards)}")
    seeds = v7r_seed_list(seed_base, n_shards)
    for k, (p, _) in enumerate(shards):
        if p.get("seed") != seeds["shards"][k] or p.get("shard") != k or not p.get("valid"):
            raise SystemExit(f"v7r shard {k}: seed/shard mismatch or invalid result")
        if p.get("seed_base") != seed_base:
            raise SystemExit(f"v7r shard {k}: seed base differs from the run")
    if ref.get("seed") != seeds["reference"][0] or not ref.get("valid"):
        raise SystemExit("v7r reference: seed mismatch or invalid result")
    all_seeds = seeds["shards"] + seeds["reference"]
    if len(set(all_seeds)) != len(all_seeds):
        raise SystemExit("v7r: the seed lists of the shards and the reference intersect")
    bad = sorted(set(all_seeds) & consumed_seed_set())
    if bad:
        raise SystemExit(f"v7r: seeds {bad} belong to a consumed evidence seed set (fail closed)")
    reps = {p["replicates"] for p, _ in shards}
    if len(reps) != 1:
        raise SystemExit("v7r: all shards must hold equal replicate counts (pairing shard s with s + 8)")
    rps = reps.pop()
    docs = [p for p, _ in shards]
    sums = np.concatenate([load_v7r_sidecar(p, path) for p, path in shards])
    ref_mean_en = np.array(ref["estimators"][SIDECAR_ESTIMATOR]["mean"], float)
    ref_blocks_en = np.array(ref["estimators"][SIDECAR_ESTIMATOR]["blocks"], float)
    estimators: dict[str, Any] = {}
    for name, kind in V7R_ESTIMATORS:
        try:
            means = np.concatenate([np.array(p["estimators"][name]["mean"], float) for p in docs])
            sems = np.concatenate([np.array(p["estimators"][name]["sem"], float) for p in docs])
            rmean = np.array(ref["estimators"][name]["mean"], float)
            rsem = np.array(ref["estimators"][name]["sem"], float)
            rblocks = np.array(ref["estimators"][name]["blocks"], float)
        except (KeyError, TypeError, ValueError) as exc:
            raise SystemExit(f"v7r: estimator {name} missing or malformed in a partial: {exc!r}") from exc
        if (means.ndim != 2 or means.shape != sems.shape or rmean.shape != (means.shape[1],)
                or rsem.shape != rmean.shape or not np.all(np.isfinite([means, sems]))
                or rblocks.shape != (v7r.V7R_REF_BLOCKS, means.shape[1])
                or not np.all(np.isfinite(rblocks)) or not np.all(np.isfinite([rmean, rsem]))):  # fmt: skip
            raise SystemExit(f"v7r: estimator {name}: inconsistent shapes or non-finite values")
        verdict = v5b.v7_estimator_verdict(kind, means, sems, rmean, rsem, rblocks, n_ref=int(ref["n"]),
                                           n_rep=v5b.V7_REP_N, boot_seed=int(ref["seed"]),
                                           n_boot=n_boot)  # fmt: skip
        if name != SIDECAR_ESTIMATOR:
            estimators[name] = verdict
            continue
        spread = ref_blocks_en.std(axis=0, ddof=1)
        ref_used = bool(ref_mean_en[0] > 0.0 and spread[0] > 0.0)
        gate = v7r.escaped_neutral_verdict(sums, seed_base, boot_seed=int(ref["seed"]),
                                           ref_used=ref_used, reps_per_shard=rps, n_boot=n_boot)  # fmt: skip
        gate["amendment13_rule_report_only"] = verdict  # never used in ``pass``
        gate["amendment13_rule_report_only_note"] = (
            "the 20-block Student-t rule of Amendment 13; it cannot turn the Amendment 13 failure "
            "into a pass (Amendment 14 (h))")
        gate["reference_bin_used"] = ref_used
        estimators[name] = gate
    gates = bool(all(e["pass"] for e in estimators.values()))
    return {
        "replicates": int(rps * n_shards), "estimators": estimators, "units": v7r.V7R_UNITS,
        "histories_per_replicate": v5b.V7_REP_N, "reference_histories": ref["n"],
        "shard_histories": [p["n"] for p in docs], "bins": v5b.V7_BINS, "seed_base": seed_base,
        "seeds": seeds, "seeds_disjoint_from_consumed": True,
        "consumed_seed_bases": list(CONSUMED_SEED_BASES_C), "rng": v7r.rng_record(),
        "rule": {"region": list(v7r.V7R_REGION), "B": v7r.V7R_B, "alpha": v7r.V7R_ALPHA,
                 "ranks": [v7r.V7R_RANK_LO, v7r.V7R_RANK_HI], "quantile_rule": v7r.V7R_QUANTILE_RULE,
                 "seed_rule": v7r.SEED_RULE, "alpha_tost": v5b.V7_ALPHA_TOST,
                 "pairs": int(rps * n_shards) // 2,
                 "pairing": "(shard s, replicate j) with (shard s + 8, replicate j)",
                 "heldout": "shards 12-15", "evaluation": "shards 0-11", "min_pairs": v5b.V7_MIN_INTERVALS,
                 "profile_estimators": "Amendment 13 two-gate rule, frozen functions",
                 "escaped_neutral": "bootstrap-t two-gate rule; Amendment 13 rule report-only"},
        "gates_pass": gates,
        "row_pass_condition": (
            "V7-R passes iff every pytest-v7r-calibration[-s{k}] step passed at the evidence SHA "
            "(identical hashed-source digest, Amendment 17 (a)7) AND gates_pass; the conjunction is "
            "made by summarize.py (a step cannot read the pytest step status)"),
        "pass": gates,  # the gates alone; summarize.py adds the calibration conjunction
    }  # fmt: skip


def step_v7r_combine(a: argparse.Namespace) -> int:
    names = [f"v7r-s{k}.json" for k in range(V7R_SHARDS)] + ["v7r-ref.json"]
    parts = load_partials_c(a, names)
    doc = {"step": "v7r-combine", "table": v5.table_record(), **elastic_table_record(),
           **v7r_combine(parts[:-1], parts[-1][0], base.SEED_BASE),
           "attestation": v5.attestation_block(a)}  # fmt: skip
    n = sum(p["n"] for p, _ in parts)
    return v5b.finish5b(doc, V7R_SHARDS * v7r.V7R_SHARD_N + v7r.V7R_REF_N, n,
                        any(p["reduced"] for p, _ in parts))  # fmt: skip


# -- exploratory diagnostic ----------------------------------------------------------------------
def v7r_diag(shards: list[tuple[dict[str, Any], Path]]) -> dict[str, Any]:
    """Per-replicate excess kurtosis of the block sums for b in {5, 10, 20} blocks per replicate (exact
    aggregation of the recorded sums). No pass or fail. Sub-block sums (b = 40, 80) are not emitted
    by the shards (Amendment 17 (a)8): recorded as not run."""
    sums = np.concatenate([load_v7r_sidecar(p, path) for p, path in shards])
    out: dict[str, Any] = {}
    for b in v7r.V7R_AGG_B:
        agg = v7r.aggregate_blocks(sums, b)
        k = np.array([v7r.excess_kurtosis(row) for row in agg])
        fin = k[np.isfinite(k)]
        out[str(b)] = {
            "blocks_per_replicate": b, "histories_per_block": v7r.V7R_BATCH_N * (20 // b),
            "replicates": int(agg.shape[0]), "kurtosis_finite": int(fin.size),
            "excess_kurtosis_mean": float(fin.mean()) if fin.size else math.nan,
            "excess_kurtosis_median": float(np.median(fin)) if fin.size else math.nan,
            "excess_kurtosis_q": [float(q) for q in np.quantile(fin, [0.05, 0.25, 0.75, 0.95])]
            if fin.size else [],
            "pooled_excess_kurtosis_of_blocks": v7r.excess_kurtosis(agg)}  # fmt: skip
    return {"diagnostic": "excess kurtosis of the escaped-neutral block sums", "per_b": out,
            "sub_blocks_b40_b80": "not emitted (Amendment 17 (a)8); decided in Amendment 14a",
            "report_only": True}  # fmt: skip


def step_v7r_diag(a: argparse.Namespace) -> int:
    names = [f"v7r-s{k}.json" for k in range(V7R_SHARDS)]
    parts = load_partials_c(a, names)
    doc = {"step": "v7r-diag", **v7r_diag(parts), "pass": True,
           "attestation": v5.attestation_block(a)}  # fmt: skip
    n = sum(p["n"] for p, _ in parts)
    return v5b.finish5b(doc, V7R_SHARDS * v7r.V7R_SHARD_N, n, any(p["reduced"] for p, _ in parts))


# -- V11 (plan Amendment 14 (c) row V11; Amendment 14 (j) r_index 14) -----------------------------
V11_R_INDEX = 14
V11_K = {150: 0, 200: 1}  # shard index k per energy; emel and emonly of an energy share the seed (CRN)
V11_MODES = ("emel", "emonly")
V11_NAMES = tuple(f"v11-{e}-{m}.json" for e in (150, 200) for m in V11_MODES)
V11_STEP_NAMES = tuple(f"v11-ionmc-{e}-{m}" for e in (150, 200) for m in V11_MODES)


def v11_seed(energy: float) -> int:
    """``base + 1000 * 14 + k`` with k = 0 (150 MeV) or 1 (200 MeV); EM+elastic and EM-only share it."""
    return base.SEED_BASE + 1000 * V11_R_INDEX + V11_K[int(energy)]


def v11_config(energy: float, mode: str, n: int, seed: int, geometry: Any, grid: Any,
               timeout: float | None = None) -> SimulationConfig:
    """warp-cpu float64, V5 geometry, 20 batches. ``emel``: ``nuc_config_c_elastic_only`` (non-elastic
    channel off, elastic on); ``emonly``: ``nuc_config_c(nuclear=False)``, the same EM model as the
    frozen lv5b V5 "off" producer (the elastic switches have no effect without ``nuclear=True``)."""
    if mode not in V11_MODES:
        raise SystemExit(f"v11: mode {mode!r} not in {V11_MODES}")
    build = nuc_config_c_elastic_only if mode == "emel" else nuc_config_c
    cfg = build(energy=energy, n=n, seed=seed, geometry=geometry, grid=grid,
                nuclear=mode == "emel", n_batches=v5b.V5_BATCHES, timeout=timeout)  # fmt: skip
    return replace(cfg, run=replace(cfg.run, backend="warp-cpu", precision="float64"))


def step_v11_ionmc(a: argparse.Namespace) -> int:
    from ionmc.reference import metrics as rm

    e, mode = a.energy, a.mode
    row = f"v11-{e:g}-{mode}"
    seed = v11_seed(e)
    n = v5b.scaled(v5b.V5_N, a.scale, 2000, v5b.V5_BATCHES)
    geo, grid, nz, r_mm = v5b.v5_geometry(e)
    cfg = v11_config(e, mode, n, seed, geo, grid, timeout=a.timeout)
    t0 = time.perf_counter()
    res = v5b.Simulation(cfg).run()
    wall = time.perf_counter() - t0
    g = res.grid("dose")
    rho = v5b.M.WATER.density_g_cm3
    bin_g_cm2 = rho * v5b.V5_DZ_MM / 10.0
    idd_b = g.batch_energy_mev.reshape(v5b.V5_BATCHES, nz) / bin_g_cm2  # MeV/(g/cm^2)/primary
    idd = idd_b.mean(axis=0)
    depth = (np.arange(nz) + 0.5) * v5b.V5_DZ_MM
    curve = rm.normalize_to_peak(idd)
    plateau = (depth >= 20.0) & (depth <= 60.0)
    tot_b = idd_b.sum(axis=1) * bin_g_cm2
    met = {"peak_depth_mm": rm.peak_depth(depth, curve), "r80_mm": rm.r80(depth, curve),
           "plateau_mean_20_60_mm": float(idd[plateau].mean()),
           "peak_over_plateau": float(idd.max() / idd[plateau].mean()),
           "total_in_grid_mev_per_primary": float(tot_b.mean()),
           "in_grid_over_e0": float(tot_b.mean() / e)}  # fmt: skip
    ok = bool(v5.clean(res) and res.energy_balance.relative_residual <= 1e-12)
    part = write_partial_c(a, row, {
        "row": row, "mode": mode, "energy": e, "nuclear": mode == "emel", "elastic_only": mode == "emel",
        "seed": seed, "n": n, "n_batches": v5b.V5_BATCHES, "bin_mm": v5b.V5_DZ_MM,
        "lateral_half_mm": v5b.V5_HALF_MM, "density_g_cm3": rho, "unit": "MeV/(g/cm^2)/primary",
        "idd_batches": idd_b.tolist(), "metrics": met, "valid": ok, "reduced": n < v5b.V5_N})  # fmt: skip
    doc = {"step": "v11-ionmc", "row": row, "mode": mode, "seed": seed, "range_csda_mm": r_mm,
           "table": v5.table_record(), **elastic_table_record(), "wall_s": wall,
           "hist_per_s": n / wall, "metrics": met,
           "relative_residual": float(res.energy_balance.relative_residual),
           "counters": res.counters.as_dict(), "partial": part,
           "comparison": "pending: step v11-compare", "pass": ok}  # fmt: skip
    return v5b.finish5b(doc, v5b.V5_N, n, n < v5b.V5_N)


def v11_seed_check(parts: list[dict[str, Any]]) -> dict[str, Any]:
    """Seed rules of row V11 (fail closed, ``SystemExit``): each energy's emel and emonly partials carry
    the shared seed ``base + 14000 + k``, the two energies differ, the set is disjoint from every
    consumed evidence seed, and every partial carries the run's seed base."""
    by = {p["row"]: p for p in parts}
    sd = parts[0].get("seed_base") if parts else None
    seeds: dict[str, int] = {}
    for e in (150, 200):
        want = int(sd) + 1000 * V11_R_INDEX + V11_K[e] if sd is not None else -1
        for m in V11_MODES:
            p = by.get(f"v11-{e}-{m}")
            if p is None or p.get("seed") != want or p.get("seed_base") != sd or not p.get("valid"):
                raise SystemExit(f"v11 {e} {m}: missing, invalid or seed != {want}")
            seeds[f"{e}-{m}"] = want
    if seeds["150-emel"] == seeds["200-emel"]:
        raise SystemExit("v11: the two energies share a seed")
    bad = sorted(set(seeds.values()) & consumed_seed_set())
    if bad:
        raise SystemExit(f"v11: seeds {bad} belong to a consumed evidence seed set (fail closed)")
    return {"seeds": seeds, "seed_base": sd, "r_index": V11_R_INDEX, "common_random_numbers": True,
            "seeds_disjoint_from_consumed": True}  # fmt: skip


def _partial_present(a: argparse.Namespace, name: str) -> bool:
    return any(p.is_file() for d in a.dirs for p in (Path(d) / name, Path(d) / "samples" / name))


def v11_bound_cases_identity(verdict: dict[str, Any], cases_dir: Path) -> dict[str, Any]:
    """Every file of every bound case directory is in ``run_suite.source_file_list('lv5c')``."""
    import run_suite

    listed = set(run_suite.source_file_list("lv5c"))
    paths: list[str] = []
    for name in sorted(v5b._bound_case_names(verdict)):
        d = cases_dir / name
        rels = sorted(f"{run_suite.V5_CASE_ROOT}/{name}/{f.relative_to(d).as_posix()}"
                      for f in d.rglob("*") if f.is_file())  # fmt: skip
        absent = [r for r in rels if r not in listed] if rels else ["<no files>"]
        if absent:
            raise SystemExit(f"bound case {name}: {absent} not in the lv5c source identity")
        paths += rels
    return {"cases_in_source_identity": True, "bound_case_paths": paths}


def _cases_dir() -> Path:
    """The frozen committed cases of this snapshot (the only directory the V11 runs may be bound to)."""
    return v5b.REPO / "validation" / "reference_cases"


def _v11_reference_runs(ref_dir: Path) -> list[Path]:
    """TOPAS runs of ``ref_dir/REF-*`` whose case.json carries a V5 (full, emonly), V11 (emelastic) or
    X-elastic-factor (noelastic) block; other materialized runs are not V11 evidence."""
    from ionmc.reference.runs import load_run

    out = []
    for d in sorted(ref_dir.glob("REF-*")):
        run = load_run(d)
        blk = run.case.get("v5")
        if run.engine == "topas" and isinstance(blk, dict) and blk.get("row") in (
                "V5", "V11", "X-elastic-factor"):  # fmt: skip
            out.append(d)
    return out


def step_v11_compare(a: argparse.Namespace) -> int:
    """Row V11 verdict from the four ``v11-*`` partials (and the optional lv5c ``v5-{e}-on`` partials
    for F_ne) against the materialized TOPAS runs through ``compare_idd_v5.build_v11_verdict``.
    Consumes partials only. Lineage or input errors give a failed document, never a silent skip."""
    import importlib.util
    import tempfile

    path = v5b.REPO / "validation" / "scripts" / "reference" / "compare_idd_v5.py"
    spec = importlib.util.spec_from_file_location("compare_idd_v5", path)
    assert spec and spec.loader
    cmp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cmp)
    parts = [p for p, _ in load_partials_c(a, list(V11_NAMES))]
    on_names = [f"v5-{e}-on.json" for e in (150, 200) if _partial_present(a, f"v5-{e}-on.json")]
    on_parts = [p for p, _ in load_partials_c(a, on_names)] if on_names else []
    seeds = v11_seed_check(parts)
    for p in on_parts:
        if p["seed"] in set(seeds["seeds"].values()):
            raise SystemExit(f"v11: V5-on seed {p['seed']} equals a V11 seed (F_ne needs independence)")
    ref_dir = Path(a.reference_dir) if a.reference_dir else v5b.REPO / ".ionmc-cache" / "reference-runs"
    cases_dir = _cases_dir()
    doc: dict[str, Any] = {"step": "v11-compare", "reference_dir": str(ref_dir),
                           "cases_dir": str(cases_dir), **seeds,
                           "v5_on_partials_used": on_names,
                           "attestation": v5.attestation_block(a)}  # fmt: skip
    reduced = any(p["reduced"] for p in [*parts, *on_parts])
    try:
        topas = _v11_reference_runs(ref_dir)
        with tempfile.TemporaryDirectory() as tmp:
            for p in [*parts, *on_parts]:
                (Path(tmp) / f"{p['row']}.json").write_text(json.dumps(p, sort_keys=True))
            verdict = cmp.build_v11_verdict(Path(tmp), topas, cases_dir)
        doc.update(v11_bound_cases_identity(verdict, cases_dir))
        doc.update(verdict=verdict, plan_rule=verdict["plan_rule"], pass_gating=bool(verdict["pass"]),
                   error=None)  # fmt: skip
        doc["pass"] = bool(verdict["pass"] and not reduced)
    except (cmp.IddError, SystemExit, ValueError, OSError) as exc:
        doc.update(cases_in_source_identity=False, verdict=None, plan_rule=cmp.V11_PLAN_RULE,
                   error=f"{type(exc).__name__}: {exc}", **{"pass": False})  # fmt: skip
    return v5b.finish5b(doc, v5b.V5_N, min(p["n"] for p in parts), reduced)


STEPS = {
    "v7r-shard": step_v7r_shard,
    "v7r-ref": step_v7r_ref,
    "v7r-combine": step_v7r_combine,
    "v7r-diag": step_v7r_diag,
    "v11-ionmc": step_v11_ionmc,
    "v11-compare": step_v11_compare,
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--workers", default="1", help="accepted for uniformity; one process")
    ap.add_argument("--shard", type=int, default=0, help="v7r-shard: shard index 0..15")
    ap.add_argument("--energy", type=float, default=150.0, help="v11-ionmc: 150 or 200")
    ap.add_argument("--mode", choices=V11_MODES, default="emel", help="v11-ionmc: emel or emonly")
    ap.add_argument("--reference-dir", default=None, help="v11-compare: dir of REF-* runs")
    ap.add_argument("--out-dir", default="samples", help="partial files")
    ap.add_argument("--dirs", nargs="+", default=["."], help="combine steps: archive directories")
    ap.add_argument("--partials-manifest", default=None, help="combine steps: manifest of imports")
    ap.add_argument("--seed-base", "--seed", dest="seed", type=int, default=LV5C_SEED_BASE,
                    help="lv5c base 20481004, rehearsals 20505000-20509499")  # fmt: skip
    ap.add_argument("--timeout", type=float, default=None)
    args = ap.parse_args(argv)
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    if args.step == "v7r-shard" and not 0 <= args.shard < V7R_SHARDS:
        raise SystemExit("v7r-shard: --shard outside 0..15")
    if args.step == "v11-ionmc" and int(args.energy) not in V11_K:
        raise SystemExit("v11-ionmc: --energy 150 or 200")
    base.SEED_BASE = args.seed
    v4.base.SEED_BASE = args.seed
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
