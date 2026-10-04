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

Deviations forced by the engine's rule ``max_step_mm <= smallest scoring spacing`` and by the
regular-grid geometry are listed in ``docs/architecture/transport.md`` (section Validation
runner): T8/T14 sigma is the standard deviation of the exit position of particles leaving a water
slab of the stated thickness (the Fermi-Eyges A2 quantity), T9 uses 1 mm IDD bins, T10 uses the
distribution of the projected track-end depth.
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
) -> Result:
    kw = {"chunk_histories": chunk} if chunk else {}
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
    backend, prec, _, b, k = SAMPLES[name]
    return parity.t12_config(
        energy_mev=a.energy,
        backend=backend,
        precision=prec,
        seed=a.seed + 1000 * k,  # a distinct seed per sample
        n_histories=sample_histories(name, a.scale),
        n_batches=b,
        workers=workers,
        lateral_bin_mm=a.lateral_bin,
        half_width_mm=a.half_width,
        timeout_s=a.timeout,
    )


def config_fingerprint(eff_summary: dict[str, Any]) -> str:
    """sha256 of the effective-configuration summary without the parallelism settings."""
    s = dict(eff_summary)
    s.pop("cpu_workers", None)
    s.pop("chunk_histories", None)
    return hashlib.sha256(json.dumps(s, sort_keys=True, default=str).encode()).hexdigest()


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
    eff = validate(cfg)
    t0 = time.perf_counter()
    if workers > 1 and cfg.run.backend != "warp-cuda":
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
        "h0": h0,
        "h1": h1,
        "energy_mev": a.energy,
        "lateral_bin_mm": a.lateral_bin,
        "half_width_mm": a.half_width,
        "frozen_histories": SAMPLES[name][2],
        "git_sha": os.environ.get("IONMC_RUN_SHA", "unknown"),
        "config_fingerprint": config_fingerprint(eff.summary()),
        "workers": workers,
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
    print("#JSON-BEGIN")
    print(json.dumps(doc, indent=1, sort_keys=True, default=_json))
    print("#JSON-END")
    return 0 if doc["pass"] else 1


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
        edep = [data[k] for k in sorted((k for k in data.files if k.startswith("edep_")),
                                        key=lambda k: int(k.split("_")[1]))]  # fmt: skip
        parts.append(
            PartialTransport(
                meta["h0"], meta["h1"], comps, [int(x) for x in data["counters"]], edep, None, {}
            )
        )
        metas.append(meta)
    if not parts:
        raise SystemExit(f"no parts of sample {name} found")
    for key in ("config_fingerprint", "seed", "n_total", "n_batches", "git_sha"):
        if len({m[key] for m in metas}) != 1:
            raise SystemExit(f"parts of sample {name} disagree on {key}")
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
    return parity.t12_observables_from_grids(grids, cfg.scoring, layout, hpb), info


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
    print("#JSON-BEGIN")
    print(json.dumps(doc, indent=1, sort_keys=True, default=_json))
    print("#JSON-END")
    return 0 if ok else 1


# -- slab observables (T8, T14) ------------------------------------------------------------------
def slab_observables(
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
    energy: float = 150.0,
    trunc_diag: bool = False,
) -> dict[str, float]:
    """Exit angle and exit-position spread of protons crossing a water slab (MCS on, straggling
    off). ``voxel`` None-like (<=0) means one voxel (box)."""
    half = 60.0
    if voxel <= 0.0:
        geo: Any = BoxPhantom((-half, -half, 0.0), (2 * half, 2 * half, thickness), WATER)
        z_end = thickness
    else:
        nz = int(round(thickness / voxel)) + (1 if shift_z else 0)
        nxy = int(round(2 * half / voxel)) + (1 if shift_xy else 0)
        z0 = -0.5 * voxel if shift_z else 0.0
        xy0 = -half - (0.5 * voxel if shift_xy else 0.0)
        shape = (nxy, nxy, nz)
        geo = VoxelGeometry(
            (xy0, xy0, z0), (voxel,) * 3, shape, (WATER,), np.zeros(shape, dtype=np.int32)
        )
        z_end = z0 + nz * voxel
    res = run_cfg(
        energy=energy,
        geometry=geo,
        scoring=(
            ScoringGrid((-half, -half, 0.0), (2 * half, 2 * half, max(thickness, smax)), (1, 1, 1)),
        ),
        seed=seed,
        n=n,
        n_batches=20,
        workers=workers,
        timeout=timeout,
        mcs=True,
        straggling=False,
        max_step=smax,
        diag=DiagnosticsOptions(escape_records=True),
        trunc_diag=trunc_diag,
    )
    d = res.diagnostics
    n_esc = len(d["escape_history"])
    u = d["escape_direction"]
    tx, ty = u[:, 0] / u[:, 2], u[:, 1] / u[:, 2]
    theta_rms = math.sqrt(float(np.mean(tx * tx + ty * ty)) / 2.0)
    p = d["escape_position_mm"]
    sigma = math.sqrt(float(np.mean(p[:, 0] ** 2 + p[:, 1] ** 2)) / 2.0)
    return {
        "theta_rms": theta_rms,
        "sigma_mm": sigma,
        "z_end_mm": z_end,
        "n_escaped": n_esc,
        "valid": float(res.valid),
    }


def fermi_eyges_a2(path: mc.ProtonPath, z_mm: float) -> float:
    """A2(z) = integral of (z - z')^2 T_dM(z') dz' over [0, z] (projected, mm^2)."""
    gx, gw = np.polynomial.legendre.leggauss(200)
    s = 0.5 * (gx + 1.0)
    pts = z_mm * s**4
    jac = 4.0 * z_mm * s**3 * 0.5 * gw
    return float(np.sum(jac * (z_mm - pts) ** 2 * path.power(pts)))


def _thickness(frac: float) -> float:
    return 5.0 * round(frac * r_csda_mm(150.0) / 5.0)


def step_t8(a: argparse.Namespace) -> int:
    path = mc.ProtonPath(TABLES, WATER, 150.0)
    out: dict[str, Any] = {}
    cases = {"theta_0.5R1": 0.5, "sigma_0.9R": 0.9}
    ok = True
    for label, frac in cases.items():
        x = _thickness(frac)
        quad_theta = math.sqrt(path.theta2_quadrature(x))
        quad_sigma = math.sqrt(fermi_eyges_a2(path, x))
        rows = {}
        for i, s in enumerate((0.1, 0.5, 1.0, 5.0)):
            r = slab_observables(
                thickness=x,
                voxel=-1.0,
                shift_z=False,
                shift_xy=False,
                smax=s,
                n=a.n,
                workers=a.workers,
                seed=a.seed + i,
                timeout=a.timeout,
            )
            rows[str(s)] = r
        key = "theta_rms" if frac == 0.5 else "sigma_mm"
        vals = [rows[k][key] for k in rows]
        spread = max(vals) / min(vals) - 1.0
        tol_pair = 0.005 if frac == 0.5 else 0.01
        ref = quad_theta if frac == 0.5 else quad_sigma
        vs_quad = max(abs(v / ref - 1.0) for v in vals) if frac == 0.5 else None
        passed = bool(
            spread <= tol_pair
            and (vs_quad is None or vs_quad <= 0.005)
            and all(r["valid"] == 1.0 for r in rows.values())
        )
        out[label] = {
            "thickness_mm": x,
            "runs": rows,
            "pairwise_spread": spread,
            "tolerance_pairwise": tol_pair,
            "quadrature": ref,
            "max_dev_from_quadrature": vs_quad,
            "pass": passed,
        }
        ok &= passed
    return finish({"step": "t8", "t8": out, "pass": ok}, 1_000_000, a.n)


def step_t14(a: argparse.Namespace) -> int:
    path = mc.ProtonPath(TABLES, WATER, 150.0)
    rows: dict[str, Any] = {}
    thick = {"0.5": _thickness(0.5), "0.9": _thickness(0.9)}
    i = 0
    for voxel in (0.5, 1.0, 2.0, 5.0):
        for sz, sxy in ((False, False), (True, True)):
            key = f"voxel{voxel}_shift{int(sz)}"
            rows[key] = {}
            for fk, x in thick.items():
                i += 1
                r = slab_observables(
                    thickness=x,
                    voxel=voxel,
                    shift_z=sz,
                    shift_xy=sxy,
                    smax=min(1.0, voxel),
                    n=a.n,
                    workers=a.workers,
                    seed=a.seed + i,
                    timeout=a.timeout,
                )
                xa = r["z_end_mm"]
                r["theta_over_quadrature"] = r["theta_rms"] / math.sqrt(path.theta2_quadrature(xa))
                r["sigma_over_fermi_eyges"] = r["sigma_mm"] / math.sqrt(fermi_eyges_a2(path, xa))
                rows[key][fk] = r
    # frozen negative control: the hinge angle sampled for the truncated length (diagnostic
    # switch) must change theta_rms by < 0.5 % at 1 mm voxels (same seed as the default run)
    x05 = thick["0.5"]
    default_1mm = slab_observables(
        thickness=x05, voxel=1.0, shift_z=False, shift_xy=False, smax=1.0, n=a.n,
        workers=a.workers, seed=a.seed + 1000, timeout=a.timeout,
    )  # fmt: skip
    control_1mm = slab_observables(
        thickness=x05, voxel=1.0, shift_z=False, shift_xy=False, smax=1.0, n=a.n,
        workers=a.workers, seed=a.seed + 1000, timeout=a.timeout, trunc_diag=True,
    )  # fmt: skip
    control = {
        "theta_rms_default": default_1mm["theta_rms"],
        "theta_rms_control": control_1mm["theta_rms"],
        "relative_change": control_1mm["theta_rms"] / default_1mm["theta_rms"] - 1.0,
        "bound": 0.005,
        "valid": bool(default_1mm["valid"] == 1.0 and control_1mm["valid"] == 1.0),
    }
    control["pass"] = bool(abs(control["relative_change"]) < 0.005 and control["valid"])
    th = [rows[k]["0.5"]["theta_over_quadrature"] for k in rows]
    s5 = [rows[k]["0.5"]["sigma_over_fermi_eyges"] for k in rows]
    s9 = [rows[k]["0.9"]["sigma_over_fermi_eyges"] for k in rows]
    res = {
        "theta_pairwise_spread": max(th) / min(th) - 1.0,
        "theta_max_dev_from_quadrature": max(abs(t - 1.0) for t in th),
        "sigma_pairwise_spread_0.5R": max(s5) / min(s5) - 1.0,
        "sigma_pairwise_spread_0.9R": max(s9) / min(s9) - 1.0,
        "sigma_max_dev_from_fermi_eyges": max(abs(t - 1.0) for t in s5 + s9),
    }
    res["pass"] = bool(
        res["theta_pairwise_spread"] <= 0.005
        and res["theta_max_dev_from_quadrature"] <= 0.005
        and max(res["sigma_pairwise_spread_0.5R"], res["sigma_pairwise_spread_0.9R"]) <= 0.01
        and res["sigma_max_dev_from_fermi_eyges"] <= 0.02
        and all(r[f]["valid"] == 1.0 for r in rows.values() for f in r)
        and control["pass"]
    )
    doc = {
        "step": "t14",
        "runs": rows,
        "t14": res,
        "pass": res["pass"],
        "negative_control": control,
    }
    return finish(doc, 1_000_000, a.n)


# -- T9 ------------------------------------------------------------------------------------------
def step_t9(a: argparse.Namespace) -> int:
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
        passed = bool(dr <= 0.1 and idd["p_value"] > parity.P_VALUE_MIN)
        out[name] = {"delta_r80_mm": dr, "idd_chi2": idd, "pass": passed}
        ok &= passed
    return finish({"step": "t9", "t9": out, "pass": ok}, None, a.n)


# -- T10 -----------------------------------------------------------------------------------------
DIRECTIONS = {
    "+x": (1.0, 0.0, 0.0),
    "+y": (0.0, 1.0, 0.0),
    "+z": (0.0, 0.0, 1.0),
    "-z": (0.0, 0.0, -1.0),
    "(1,1,0)": (1.0, 1.0, 0.0),
    "(1,1,1)": (1.0, 1.0, 1.0),
}


def step_t10(a: argparse.Namespace) -> int:
    e = 150.0
    r = r_csda_mm(e)
    geo = BoxPhantom((-100.0, -100.0, -100.0), (200.0, 200.0, 200.0), WATER)
    grid = (ScoringGrid((-100.0, -100.0, -100.0), (200.0, 200.0, 200.0), (1, 1, 1)),)
    nb = int(math.ceil(1.3 * r))
    obs = {}
    for i, (name, d) in enumerate(DIRECTIONS.items()):
        u = np.array(d) / np.linalg.norm(d)
        start = tuple(float(x) for x in (-0.5 * r * u))
        res = run_cfg(
            energy=e,
            geometry=geo,
            scoring=grid,
            seed=a.seed + i,
            n=a.n,
            n_batches=20,
            workers=a.workers,
            timeout=a.timeout,
            position=start,
            direction=d,
            diag=DiagnosticsOptions(track_end_positions=True),
        )
        pos = res.diagnostics["end_position_mm"] - np.array(start)
        t = pos @ u
        batch = np.arange(len(t)) % 20
        hist = np.zeros((20, nb))
        for b in range(20):
            hist[b] = np.histogram(t[batch == b], bins=nb, range=(0.0, float(nb)))[0]
        mean_depth = np.array([t[batch == b].mean() for b in range(20)])
        total_b = np.asarray(res.grids[0].batch_energy_mev).reshape(20)  # whole-world grid
        obs[name] = (
            parity.T12Observables(
                {"end_depth": hist}, {"mean_end_depth_mm": mean_depth}, 20, a.n // 20
            ),
            res.counters.as_dict(),
            float(t.mean()),
            total_b,
        )
    ref_name = "+z"
    out = {}
    ok = True
    for name, (o, counters, md, total_b) in obs.items():
        if name == ref_name:
            continue
        v = parity.t12_compare(obs[ref_name][0], o)
        dmean = abs(md - obs[ref_name][2])
        # frozen: total energy equal within 3 sigma (batch standard errors of both orientations)
        e_z = parity.scalar_z(obs[ref_name][3], total_b)
        oblique = name.startswith("(")
        chi_ok = v["arrays"]["end_depth"]["p_value"] > parity.P_VALUE_MIN
        crit = (dmean <= 0.3) if oblique else chi_ok
        energy_ok = abs(e_z["z"]) < 3.0
        passed = bool(crit and energy_ok and not any(counters.values()))
        out[name] = {
            "delta_mean_end_depth_mm": dmean,
            "end_depth_chi2": v["arrays"]["end_depth"],
            "permutation": v["permutation"],
            "total_energy": {**e_z, "bound": 3.0, "pass": energy_ok},
            "counters": counters,
            "pass": passed,
        }
        ok &= passed
    return finish({"step": "t10", "t10": out, "pass": ok}, None, a.n)


STEPS = {
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
    ap.add_argument("--part", default="1/1", help="t12-python-sample: history range i/n")
    ap.add_argument("--samples", default="cpu32", help="t12-accelerated-samples: names")
    ap.add_argument("--out-dir", default="samples", help="directory of the sample files")
    ap.add_argument("--dirs", nargs="+", default=["."], help="t12-compare: archive directories")
    ap.add_argument("--n", type=int, default=1_000_000)
    ap.add_argument("--k", type=int, default=256)
    ap.add_argument("--energy", type=float, default=150.0)
    ap.add_argument("--seed", type=int, default=20261004)
    ap.add_argument("--timeout", type=float, default=None)
    ap.add_argument("--backend", default="warp-cpu")
    ap.add_argument("--mode", choices=("workers", "chunks"), default="workers")
    ap.add_argument("--pairs", default="python:cpu32,cpu32:cpu64")
    ap.add_argument("--scale", type=float, default=1.0, help="T12 history-count factor")
    ap.add_argument("--lateral-bin", type=float, default=0.2)
    ap.add_argument("--half-width", type=float, default=20.0)
    args = ap.parse_args(argv)
    args.workers = (os.cpu_count() or 1) if args.workers == "auto" else int(args.workers)
    return STEPS[args.step](args)


if __name__ == "__main__":
    sys.exit(main())
