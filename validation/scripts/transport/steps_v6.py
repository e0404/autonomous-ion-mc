"""Validation steps of row V6 / V6-sens (Gottschalk 2015 177 MeV depth-dose and halo; plan Amendment 14 (g),
Amendment 17 (c); suites ``lv5c`` and ``hr5c``).

Usage::

    python validation/scripts/transport/steps_v6.py <step> [--shard K] [--variant NAME] [--scale F] ...

Steps: ``v6-shard`` (lv5c, warp-cpu float64, r_index 12), ``v6-sens-shard`` (lv5c, r_index 20, the nominal
and the eight variants with common random numbers), ``v6-cuda-shard`` (hr5c, warp-cuda float32,
r_index 3), ``v6-compare`` / ``v6-cuda-compare`` (consumers, no seed), ``v6-e0`` (the deterministic E0
record) and ``v6-unsized`` (fails closed until Amendment 14a fixes ``V6_SHARDS``).

Producers never read the measured dose table: the beam parameters are the plan-frozen ``tbl:parms``
values (``gottschalk_v6.PARMS_FROZEN``, verified against the e-print by ``v6-e0`` and ``v6-compare``), the
data depths are coordinates (``DEPTHS_CM``), and E0 comes from the range table. The dose table is read
only by the comparator inside ``step_v6_compare`` / ``step_v6_cuda_compare`` (isolation test).

Observable and scoring: an axisymmetric 0.5 mm x 0.5 mm annular (r, z) tally of the water phantom
(``ScoringGrid(kind="annular")``, 210 x 520 cells, axis x = y = 0), dose = edep / (rho V) in MeV/g per
primary, averaged over the 0.4 cm diameter sphere at each data (r, z) with the recorded overlap weights.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

import steps as base
import steps_v5 as v5
import steps_v5b as v5b
import steps_v5c as v5c
from ionmc import materials as M
from ionmc.config import SimulationConfig
from ionmc.geometry import BoxPhantom
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
V6_R_INDEX = 12  # lv5c warp-cpu float64
V6_SENS_R_INDEX = 20  # lv5c, nominal + eight variants, common random numbers
V6_CUDA_R_INDEX = 3  # hr5c warp-cuda float32
V6_SHARD_N = 1_000_000
V6_BATCHES = 20
V6_CUDA_CHUNK = 2**14
V6_SHARDS: int | None = None
"""N_V6, the number of 1e6-history shards: ``None`` until Amendment 14a fixes it after the rehearsal (the
suite then fails closed when V6 steps are requested; ``run_suite.V6_SHARDS`` must equal it)."""
LV5C_SEED_BASE = 20481004
HR5C_SEED_BASE = 20491004
FORMAT = 1
GRID_NAME = "rz"
MONITOR_NOTE = "the +-3 % monitor uncertainty of the measurement is reported, not added to any tolerance"
ANNOTATION = {"sigma_sem_cm": 0.293, "d80_range_shift_percent_per_0p1mm": 1.2,
              "d50_range_shift_percent_per_0p1mm": 2.8, "plateau_percent_per_0p1mm_below_0p95_R0": 0.2}  # fmt: skip
SEED_BASES_REFUSED = v5c.CONSUMED_SEED_BASES_C


def gv6() -> Any:
    """The comparator module (``validation/scripts/reference/gottschalk_v6.py``), loaded by path."""
    if "gottschalk_v6" in sys.modules:
        return sys.modules["gottschalk_v6"]
    path = REPO / "validation" / "scripts" / "reference" / "gottschalk_v6.py"
    spec = importlib.util.spec_from_file_location("gottschalk_v6", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gottschalk_v6"] = mod
    spec.loader.exec_module(mod)
    return mod


# -- beam, geometry, grid -----------------------------------------------------------------------
def e0_functions() -> dict[str, Any]:
    """The range-table callables of the E0 procedure (the same call as ``steps_v5.depth_box``)."""
    t = BetheStoppingSource().table(M.WATER, PROTON)
    return {"rho": M.WATER.density_g_cm3, "range": lambda e: float(t.range_at(e)),
            "inverse": lambda r: float(t.energy_from_range(r)),
            "stopping": lambda e: float(t.stopping_at(e)), "table": t}  # fmt: skip


def e0_for_r0(r0_cm: float) -> tuple[float, float]:
    """``(E0, sigma_E)`` for a range ``R0`` [cm]: bisection of the water range table and
    ``sigma_E = 0.008 R0 rho S_mass(E0)`` (the plan's conversion with dR/dE = 1/S of the same table)."""
    g, f = gv6(), e0_functions()
    e0, _, _ = g.solve_e0(f["range"], r0_cm * f["rho"])
    return e0, g.SIGMA_R_FRACTION * r0_cm * f["rho"] * f["stopping"](e0)


def beam_table() -> dict[str, dict[str, float]]:
    return gv6().beam_parameters(gv6().PARMS_FROZEN, e0_for_r0)


def v6_geometry() -> tuple[BoxPhantom, ScoringGrid]:
    g = gv6()
    geo = BoxPhantom((-200.0, -200.0, 0.0), (400.0, 400.0, 260.0), M.WATER)
    grid = ScoringGrid((0.0, 0.0, 0.0), (g.GRID_DR_MM, g.GRID_DR_MM, g.GRID_DZ_MM),
                       (g.GRID_NR, 1, g.GRID_NZ), name=GRID_NAME, kind="annular")  # fmt: skip
    return geo, grid


def grid_record(grid: ScoringGrid) -> dict[str, Any]:
    return {"kind": grid.kind, "origin_mm": list(grid.origin_mm), "spacing_mm": list(grid.spacing_mm),
            "shape": list(grid.shape), "phantom": "water box [-200,200]^2 x [0,260] mm, beam at (0,0,0) along +z"}  # fmt: skip


def v6_config(variant: str, n: int, seed: int, backend: str, precision: str,
              timeout: float | None = None) -> SimulationConfig:
    bp = beam_table()[variant]
    geo, grid = v6_geometry()
    cfg = v5c.wcfg_c(backend, precision, energy=bp["E0_MeV"], n=n, seed=seed, geometry=geo, grid=grid,
                     n_batches=V6_BATCHES, timeout=timeout)  # fmt: skip
    src = PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), bp["E0_MeV"],
                           energy_sigma_mev=bp["sigma_E_MeV"], lateral_sigma_mm=bp["sigma_x_mm"],
                           angular_sigma_mrad=bp["sigma_theta_mrad"], x_theta_correlation=bp["rho_xtheta"])  # fmt: skip
    cfg = replace(cfg, source=src)
    if backend == "warp-cuda":
        cfg = replace(cfg, run=replace(cfg.run, chunk_histories=V6_CUDA_CHUNK))
    return cfg


_WEIGHTS: dict[str, Any] = {}


def weights() -> tuple[dict[tuple[int, int], Any], str]:
    if not _WEIGHTS:
        g = gv6()
        w = g.all_sphere_weights(g.DEPTHS_CM, g.RADII_CM)
        _WEIGHTS.update(w=w, digest=g.weights_digest(w))
    return _WEIGHTS["w"], _WEIGHTS["digest"]


# -- seeds and names ----------------------------------------------------------------------------
def shard_names(n: int | None = None) -> tuple[str, ...]:
    n = require_shards(n)
    return tuple(f"v6-s{k}.json" for k in range(n))


def sens_names(n: int | None = None) -> tuple[str, ...]:
    n = require_shards(n)
    return tuple(f"v6-sens-{v}-s{k}.json" for v in gv6().SENS_NAMES for k in range(n))


def cuda_names(n: int | None = None) -> tuple[str, ...]:
    n = require_shards(n)
    return tuple(f"v6-cuda-s{k}.json" for k in range(n))


def require_shards(n: int | None = None) -> int:
    n = V6_SHARDS if n is None else n
    if n is None:
        raise SystemExit("V6_SHARDS is not fixed: N_V6 comes from the C8 rehearsal and Amendment 14a (fail closed)")
    return n


def seed_of(r_index: int, k: int) -> int:
    return base.SEED_BASE + 1000 * r_index + k


def check_base(seed_base: int) -> None:
    if seed_base in SEED_BASES_REFUSED:
        raise SystemExit(f"seed base {seed_base} is a consumed evidence base (fail closed)")


# -- producers ----------------------------------------------------------------------------------
def _produce(a: argparse.Namespace, name: str, kind: str, variant: str, r_index: int, backend: str,
             precision: str) -> int:
    g = gv6()
    k = a.shard
    seed = seed_of(r_index, k)
    n = v5b.scaled(V6_SHARD_N, a.scale, 2000, V6_BATCHES)
    bp = beam_table()[variant]
    cfg = v6_config(variant, n, seed, backend, precision, a.timeout)
    t0 = time.perf_counter()
    res = Simulation(cfg).run()
    wall = time.perf_counter() - t0
    grid = res.grid(GRID_NAME)
    nb = V6_BATCHES
    dose = (grid.batch_energy_mev.reshape(nb, g.GRID_NR, g.GRID_NZ)
            / grid.mass_g.reshape(1, g.GRID_NR, g.GRID_NZ))  # MeV/g per primary  # fmt: skip
    w, digest = weights()
    avg = g.sphere_average_batches(dose, w, len(g.DEPTHS_CM), len(g.RADII_CM))
    eb = res.energy_balance
    ok = bool(v5.clean(res) and eb.relative_residual <= 1e-12 and np.all(np.isfinite(avg)))
    part = v5c.write_partial_c(a, name, {
        "row": "V6-sens" if kind == "sens" else "V6", "name": name, "kind": kind, "variant": variant,
        "shard": k, "seed": seed, "r_index": r_index, "n": n, "n_batches": nb, "backend": backend,
        "precision": precision, "beam": bp, "grid": grid_record(v6_geometry()[1]),
        "weights_digest": digest, "sphere": {"diameter_mm": g.SPHERE_DIAMETER_MM, "subsamples": g.SPHERE_SUBSAMPLES},
        "depths_cm": list(g.DEPTHS_CM), "radii_cm": list(g.RADII_CM), "unit": "MeV/g per primary",
        "dose_batches": avg.tolist(), "valid": ok, "reduced": n < V6_SHARD_N,
        "outside_over_initial": float(eb.outside_mev[0] / eb.initial_mev),
        "convention": "(tan thx, tan thy, 1) normalised; theta = sigma_theta (rho g1 + sqrt(1 - rho^2) g2)",
    })  # fmt: skip
    doc = {"step": {"v6": "v6-shard", "sens": "v6-sens-shard", "cuda": "v6-cuda-shard"}[kind],
           "name": name.removesuffix(".json"), "row": "V6", "kind": kind, "variant": variant, "seed": seed,
           "r_index": r_index, "backend": backend, "precision": precision, "beam": bp, "wall_s": wall,
           "hist_per_s": n / wall, "weights_digest": digest, "counters": res.counters.as_dict(),
           "relative_residual": float(eb.relative_residual), "partial": part, "pass": ok,
           "comparison": "pending: v6-compare"}  # fmt: skip
    return v5b.finish5b(doc, V6_SHARD_N, n, n < V6_SHARD_N)


def step_v6_shard(a: argparse.Namespace) -> int:
    return _produce(a, f"v6-s{a.shard}", "v6", "nominal", V6_R_INDEX, "warp-cpu", "float64")


def step_v6_sens_shard(a: argparse.Namespace) -> int:
    if a.variant not in gv6().SENS_NAMES:
        raise SystemExit(f"v6-sens-shard: --variant must be one of {gv6().SENS_NAMES}")
    return _produce(a, f"v6-sens-{a.variant}-s{a.shard}", "sens", a.variant, V6_SENS_R_INDEX,
                    "warp-cpu", "float64")  # fmt: skip


def cuda_target() -> tuple[str, str]:
    """``(backend, precision)`` of the hr5c producers: warp-cuda float32. A test-only override
    (``IONMC_V6_CUDA_TEST_BACKEND``, e.g. warp-cpu) emulates the sens-free CUDA path on a host without a GPU; it
    is refused unless running under pytest."""
    over = os.environ.get("IONMC_V6_CUDA_TEST_BACKEND")
    if over:
        if "PYTEST_CURRENT_TEST" not in os.environ:
            raise SystemExit("IONMC_V6_CUDA_TEST_BACKEND is a test-only override (refused outside pytest)")
        return over, "float64"  # float64 so that the emulation can be made bit-identical to the lv5c sample
    return "warp-cuda", "float32"


def step_v6_cuda_shard(a: argparse.Namespace) -> int:
    backend, precision = cuda_target()
    return _produce(a, f"v6-cuda-s{a.shard}", "cuda", "nominal", V6_CUDA_R_INDEX, backend, precision)


def step_v6_unsized(a: argparse.Namespace) -> int:
    require_shards()
    return 0


# -- E0 record ----------------------------------------------------------------------------------
def step_v6_e0(a: argparse.Namespace) -> int:
    """Deterministic E0 record (plan (g)); writes ``e0_record.json`` into ``--out-dir`` and fails
    closed if a cross-check exceeds 1e-6 relative or the e-print ``tbl:parms`` differs from the
    frozen values."""
    import json

    g = gv6()
    f = e0_functions()
    geo, grid = v6_geometry()
    cfg = v6_config("nominal", 2000, seed_of(V6_R_INDEX, 0), "warp-cpu", "float64")
    eff = Simulation(cfg).effective
    tabs = eff.tables
    ident = {"effective_transport_table": "ionmc.transport.tables.TransportTables of the effective "
                                          "nuc_config_c configuration (water = material 0), built from the Bethe stopping table",
             "transport_tables_sha256": tabs.sha256, "stopping_source": cfg.physics.stopping.name,
             "range_construction": "exact-loglog-quadrature-v1", "water_I_eV": M.WATER.mean_excitation_eV,
             "primary_range_call": "BetheStoppingSource().table(M.WATER, PROTON).range_at (as steps_v5.depth_box)"}  # fmt: skip
    rec = g.e0_record(g.PARMS_FROZEN["p1"], f["rho"], f["range"], f["inverse"], f["stopping"],
                      lambda e: float(tabs.range_g_cm2(0, e)), ident, code_sha=base.os.environ.get("IONMC_RUN_SHA"))
    record: dict[str, Any] = {}
    try:
        tex, record = g.load_halo_tex()
        parms = g.parse_parms(tex)
        rec["tbl_parms_match_frozen"] = all(parms[k] == v for k, v in g.PARMS_FROZEN.items())
        rec["beam"] = g.beam_numbers(parms)
    except (g.V6DataError, ValueError) as exc:
        rec["tbl_parms_match_frozen"] = False
        rec["data_error"] = f"{type(exc).__name__}: {exc}"
        record = getattr(exc, "record", {})
    rec["data"] = record
    rec["beam_table"] = beam_table()
    rec["pass"] = bool(rec["pass"] and rec["tbl_parms_match_frozen"])
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "e0_record.json").write_text(json.dumps(rec, indent=1, sort_keys=True, default=base._json) + "\n")
    return base.emit({"step": "v6-e0", **rec})


# -- comparison ---------------------------------------------------------------------------------
def _load_group(a: argparse.Namespace, names: tuple[str, ...], r_index: int, kind: str, variant: str,
                backend: str, precision: str, n_shards: int) -> list[dict[str, Any]]:
    docs = [d for d, _ in v5c.load_partials_c(a, list(names))]
    for k, d in enumerate(docs):
        want = d.get("seed_base", -1) + 1000 * r_index + k
        if (d.get("seed") != want or d.get("shard") != k or d.get("kind") != kind or d.get("variant") != variant
                or d.get("backend") != backend or d.get("precision") != precision or not d.get("valid")):
            raise SystemExit(f"{d.get('name')}: shard binding (seed/shard/kind/variant/backend/valid) mismatch")
    if len(docs) != n_shards:
        raise SystemExit("incomplete shard set")
    return docs


def _consistent(docs: list[dict[str, Any]]) -> dict[str, Any]:
    keys = ("beam", "grid", "weights_digest", "n", "n_batches", "seed_base", "r_index", "unit", "sphere")
    first = docs[0]
    for d in docs[1:]:
        for k in keys:
            if d.get(k) != first.get(k):
                raise SystemExit(f"shards disagree on {k}: {d.get('name')}")
    return {k: first[k] for k in keys}


def _stack(docs: list[dict[str, Any]]) -> Any:
    return np.concatenate([np.array(d["dose_batches"], dtype=np.float64) for d in docs], axis=0)


def _d80(z_cm: Any, curve: Any) -> float | None:
    j = int(np.argmax(curve))
    thr = 0.8 * curve[j]
    for i in range(j, len(z_cm) - 1):
        if curve[i] >= thr > curve[i + 1]:
            t = (curve[i] - thr) / (curve[i] - curve[i + 1])
            return float(z_cm[i] + t * (z_cm[i + 1] - z_cm[i]))
    return None


def _data(g: Any) -> tuple[Any, Any, Any, dict[str, Any]]:
    tex, record = g.load_halo_tex()
    z, r, logd, nmiss = g.parse_dmlg(tex)
    parms = g.parse_parms(tex)
    meas = g.measured_dose(logd)
    return z, r, meas, {"record": record, "n_points": int(np.isfinite(meas).sum()), "n_missing": nmiss,
                        "sha256_tarball": record["sha256_tarball"], "sha256_halo_tex": record["sha256_halo_tex"],
                        "tbl_parms": {k: parms[k] for k in g.PARMS_FROZEN}}  # fmt: skip


def _beam_doc(g: Any, common: dict[str, Any]) -> dict[str, Any]:
    b = g.beam_numbers(g.PARMS_FROZEN)
    return {**common["beam"], "A1_cm_rad": b["A1_cm_rad"], "emittance_pi_mm_mrad": b["emittance_pi_mm_mrad"],
            "f": b["f"], "R0_cm": g.PARMS_FROZEN["p1"], "convention": b["convention"],
            "roles": {"tbl:dmlg": g.DATA_ROLE, "tbl:parms": g.BEAM_ROLE}}  # fmt: skip


def _compare(a: argparse.Namespace, *, cuda: bool) -> int:
    g = gv6()
    n_shards = require_shards()
    check_base(base.SEED_BASE)
    step = "v6-cuda-compare" if cuda else "v6-compare"
    doc: dict[str, Any] = {"step": step, "row": "V6", "lineage": g.LINEAGE, "attestation": v5.attestation_block(a)}
    reasons: list[str] = []
    reasons_extra: list[str] = []
    try:
        if cuda:
            docs = _load_group(a, cuda_names(), V6_CUDA_R_INDEX, "cuda", "nominal", *cuda_target(), n_shards)
        else:
            docs = _load_group(a, shard_names(), V6_R_INDEX, "v6", "nominal", "warp-cpu", "float64", n_shards)
        common = _consistent(docs)
        _, wdigest = weights()
        if common["weights_digest"] != wdigest or common["beam"] != beam_table()["nominal"]:
            raise SystemExit("shard weights digest or beam differ from the current code")
        sens = nominal = None
        sens_docs: dict[str, list[dict[str, Any]]] = {}
        if not cuda:
            for v in g.SENS_NAMES:
                sens_docs[v] = _load_group(a, tuple(f"v6-sens-{v}-s{k}.json" for k in range(n_shards)),
                                           V6_SENS_R_INDEX, "sens", v, "warp-cpu", "float64", n_shards)
                if _consistent(sens_docs[v])["beam"] != beam_table()[v]:
                    raise SystemExit(f"v6-sens {v}: beam differs from the current code")
            sens_seed = {v: [d["seed"] for d in ds] for v, ds in sens_docs.items()}
            if len({tuple(s) for s in sens_seed.values()}) != 1:
                raise SystemExit("v6-sens variants do not share seeds (common random numbers)")
            nominal = _stack(sens_docs["nominal"])
            sens = {v: _stack(ds) for v, ds in sens_docs.items() if v != "nominal"}
        own = [d["seed"] for d in docs]
        if len(set(own)) != len(own) or set(own) & set(v5c.consumed_seed_set()):
            raise SystemExit("V6 seeds overlap or are consumed (fail closed)")
        sens_seeds = sorted({d["seed"] for ds in sens_docs.values() for d in ds})
        if set(own) & set(sens_seeds):
            raise SystemExit("the V6 r12 sample shares seeds with V6-sens (must be independent)")
        all_seeds = [*own, *sens_seeds]
        z, r, meas, data = _data(g)
        batches = _stack(docs)
        lv = _optional_verdict(a) if cuda else None
        fixed = lv["status"] if lv is not None and lv.get("status") in ("gating", "report-only") else None
        ver = g.evaluate(z, r, meas, batches, sens, nominal, fixed_status=fixed)
        if not cuda:
            ver["v6_sens"]["r0_variants_report"] = g.r0_variants_report(
                z, r, sens_all(sens_docs), nominal, ver["halo"]["points"])
        mean = batches.mean(axis=0)
        r0c = int((r == 0).nonzero()[0][0])
        d80 = _d80(z, mean[:, r0c])
        doc.update(
            data=data, beam=_beam_doc(g, common), n_histories=sum(d["n"] for d in docs), n_shards=n_shards,
            n_batches_total=int(batches.shape[0]), grid=common["grid"], weights_digest=common["weights_digest"],
            r0=ver["r0"], halo=ver["halo"], v6_sens=ver["v6_sens"], status=ver["status"],
            d80_ionmc_cm=d80, R0_cm=g.PARMS_FROZEN["p1"], d80_note="ionmc's own d80 from the sphere-averaged r = 0 curve at the data depths; never used to adjust E0",
            annotation=ANNOTATION, monitor_uncertainty_note=MONITOR_NOTE, seeds={"v6": own, "all": sorted(set(all_seeds))},
            plan_rule=g.PLAN_QUOTE_STATUS,
        )  # fmt: skip
        if cuda:
            doc["lv5c_verdict"] = lv
            if lv is not None and fixed is not None:
                comb = g.combine_backends(lv, ver)
                doc["backend_combination"] = comb
                doc["status"] = lv["status"]
                doc["evaluated_status_both"] = comb["evaluated_status_both"]
            else:
                reasons_extra.append("no usable sealed lv5c verdict (v6-lv5c-verdict.json): fail closed")
        else:
            v5c_seal = v5c.write_partial_c(a, LV5C_VERDICT.removesuffix(".json"), {"row": "V6", "verdict": _compact(ver)})
            doc["verdict_partial"] = v5c_seal
        status = doc["status"]
        if cuda:
            gp = doc.get("backend_combination", {}).get("pass")
            doc["gate_pass"] = gp
        else:
            gp = ver["gate_pass"]
            doc["gate_pass"] = gp
        reasons = [*_reasons(doc), *reasons_extra]
        if cuda:  # the backend combination decides (both statuses); no sealed verdict: fail closed
            comb_pass = doc.get("backend_combination", {}).get("pass")
            doc["pass"] = bool(comb_pass) if comb_pass is not None else False
            if comb_pass is False and doc["backend_combination"]["disagreement_components"]:
                reasons.append("backend disagreement on " + ", ".join(doc["backend_combination"]["disagreement_components"]))
        else:
            doc["pass"] = bool(gp) if status == "gating" else (status == "report-only" and True)
        if status == "pending":
            doc["pass"] = False
    except g.V6DataError as exc:
        doc.update(error=f"V6DataError: {exc}", data_record=exc.record)
        doc["pass"] = False
        reasons = [f"data pin mismatch recorded: {exc}"]
    except (ValueError, KeyError) as exc:
        doc.update(error=f"{type(exc).__name__}: {exc}")
        doc["pass"] = False
        reasons = [str(exc)]
    doc["reasons"] = reasons
    return v5b.finish5b(doc, V6_SHARD_N * n_shards, V6_SHARD_N * n_shards, False)


def sens_all(sens_docs: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    return {v: _stack(ds) for v, ds in sens_docs.items() if v != "nominal"}


def _compact(ver: dict[str, Any]) -> dict[str, Any]:
    return {"r0": ver["r0"], "halo": {k: v for k, v in ver["halo"].items() if k != "points"} | {"points": ver["halo"]["points"]},
            "status": ver["status"], "gate_pass": ver["gate_pass"], "v6_sens": {"condition_holds": ver["v6_sens"]["condition_holds"]}}  # fmt: skip


LV5C_VERDICT = "v6-lv5c-verdict.json"
LV5C_VERDICT_BINDINGS = ("git_sha", "run_sha", "table_id", "table_npz_sha256", "elastic", "elastic_table_id",
                         "elastic_table_npz_sha256")  # fmt: skip


def _optional_verdict(a: argparse.Namespace) -> dict[str, Any] | None:
    """The sealed lv5c verdict of ``v6-compare`` (``None`` if absent; any defect is a SystemExit). It is written
    by another suite (``lv5c``, its own seed base), so the run bindings of ``load_partials`` (suite and seed base)
    do not apply; instead the seal, the manifest entry of an imported file, the code and table bindings, the
    producing suite and the recorded seed base are checked."""
    hits = [(i, p) for i, d in enumerate(a.dirs) for p in (Path(d) / LV5C_VERDICT, Path(d) / "samples" / LV5C_VERDICT)
            if p.is_file()]  # fmt: skip
    if not hits:
        return None
    if len(hits) != 1:
        raise SystemExit(f"{LV5C_VERDICT}: found {len(hits)} copies in {a.dirs} (need exactly 1)")
    idx, hit = hits[0]
    doc = json.loads(hit.read_text())
    digest = v5.content_digest(doc)
    if doc.get("content_sha256") != digest:
        raise SystemExit(f"{LV5C_VERDICT}: content_sha256 does not match the document")
    record: dict[str, Any] = {"name": LV5C_VERDICT, "content_sha256": digest, "origin": "current-run", "host_run_id": None}
    if idx != 0:
        manifest, msha = v5.read_manifest(getattr(a, "partials_manifest", None))
        ent = manifest.get(LV5C_VERDICT)
        if ent is None or ent["content_sha256"] != digest:
            raise SystemExit(f"{LV5C_VERDICT}: imported but not in the partials manifest with this digest")
        record.update(origin="imported", host_run_id=ent["host_run_id"])
        a.attested_manifest_path, a.attested_manifest_sha256 = str(a.partials_manifest), msha
    want = v5c.bindings_c(a)
    for key in LV5C_VERDICT_BINDINGS:
        if doc.get(key) != want[key]:
            raise SystemExit(f"{LV5C_VERDICT}: {key} is {doc.get(key)!r}, expected {want[key]!r}")
    if doc.get("suite") != "lv5c" or doc.get("seed_base") in SEED_BASES_REFUSED:
        raise SystemExit(f"{LV5C_VERDICT}: not a verdict of an unconsumed lv5c run")
    if not hasattr(a, "attested_partials"):
        a.attested_partials = []
    a.attested_partials.append(record)
    return doc["verdict"]


def _reasons(doc: dict[str, Any]) -> list[str]:
    out: list[str] = []
    if not doc["r0"]["pass"]:
        out.append(f"r = 0: {sum(not p['pass'] for p in doc['r0']['points'])} point(s) outside 2 %")
    if not doc["halo"]["pass"]:
        h = doc["halo"]
        out.append(f"halo: frac_pass {h['frac_pass']:.3f}, median {h['median_abs']:.3f}, sign p {h['sign_test']['p']:.3g}")
    if doc["status"] == "report-only":
        out.append("V6-sens condition failed: V6 as a whole is report-only (recorded)")
    if doc["status"] == "pending":
        out.append("status pending: V6-sens samples absent")
    return out


def step_v6_compare(a: argparse.Namespace) -> int:
    return _compare(a, cuda=False)


def step_v6_cuda_compare(a: argparse.Namespace) -> int:
    return _compare(a, cuda=True)


STEPS = {
    "v6-shard": step_v6_shard, "v6-sens-shard": step_v6_sens_shard, "v6-cuda-shard": step_v6_cuda_shard,
    "v6-compare": step_v6_compare, "v6-cuda-compare": step_v6_cuda_compare, "v6-e0": step_v6_e0,
    "v6-unsized": step_v6_unsized,
}  # fmt: skip


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("step", choices=sorted(STEPS))
    ap.add_argument("--scale", type=float, default=1.0, help="history-count factor (reduced run)")
    ap.add_argument("--workers", default="1", help="accepted for uniformity; one process")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--variant", default="nominal", help="v6-sens-shard: nominal, sx-lo, ... r0-hi")
    ap.add_argument("--out-dir", default="samples")
    ap.add_argument("--dirs", nargs="+", default=["."])
    ap.add_argument("--partials-manifest", default=None)
    ap.add_argument("--seed-base", "--seed", dest="seed", type=int, default=LV5C_SEED_BASE)
    ap.add_argument("--timeout", type=float, default=None)
    args = ap.parse_args(argv)
    if not 0.0 < args.scale <= 1.0:
        raise SystemExit("need 0 < --scale <= 1")
    if args.step in ("v6-shard", "v6-sens-shard", "v6-cuda-shard") and args.shard < 0:
        raise SystemExit("--shard must be >= 0")
    check_base(args.seed)
    base.SEED_BASE = args.seed
    v5c.v4.base.SEED_BASE = args.seed
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
