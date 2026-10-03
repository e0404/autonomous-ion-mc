"""EM-only transport benchmark harness (informative aggregates, NOT performance targets).

No performance targets exist yet (``targets: null``; they are frozen in V3-012, and the protocol
forbids targets derived from stripped physics). Every result is labelled ``"EM-only"``:
electromagnetic transport without nuclear interactions, so it is not comparable with a full-physics
run.

Workload (default): 150 MeV protons in a 120 x 120 x 300 mm water box with 1 mm transport voxels,
a 2 mm scoring grid, float32, 1e6 histories. Backends: warp-cpu with W workers (default 1 and
16) and warp-cuda. Every measurement runs in a fresh Python process (``--child``), 3 repeats; a
cold start (empty Warp kernel cache) is measured once per backend as repeat 0. Metrics: wall
time, histories/s, compile seconds, kernel/transport seconds, host reduction seconds, peak RSS.

Usage::

    python benchmarks/transport/bench_em.py --out benchmarks/generated/transport/<new>.json \
        [--histories 1000000] [--workers 1 16] [--cuda] [--repeats 3]

The output file must not exist and must lie under ``benchmarks/generated/``.
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"


def child(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(SRC))
    import numpy as np

    from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig
    from ionmc.geometry import VoxelGeometry
    from ionmc.materials import WATER
    from ionmc.physics.projectiles import PROTON
    from ionmc.physics.stopping import BetheStoppingSource
    from ionmc.scoring import ScoringGrid
    from ionmc.simulation import Simulation
    from ionmc.sources import PencilBeamSource

    shape = (120, 120, 300)
    geo = VoxelGeometry(
        (-60.0, -60.0, 0.0), (1.0, 1.0, 1.0), shape, (WATER,), np.zeros(shape, dtype=np.int32)
    )
    cfg = SimulationConfig(
        source=PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), 150.0),
        geometry=geo,
        scoring=(ScoringGrid((-60.0, -60.0, 0.0), (2.0, 2.0, 2.0), (60, 60, 150)),),
        physics=PhysicsOptions(nuclear=False, stopping=BetheStoppingSource(), max_step_mm=1.0),
        run=RunOptions(
            backend=args.backend,
            precision="float32",
            seed=1,
            n_histories=args.histories,
            n_batches=20,
            cpu_workers=args.child_workers,
        ),
    )
    t0 = time.perf_counter()
    sim = Simulation(cfg)
    res = sim.run()
    wall = time.perf_counter() - t0
    parts = res.transport_report.get("partials", [])
    out = {
        "label": "EM-only",
        "backend": args.backend,
        "workers": args.child_workers,
        "histories": args.histories,
        "wall_s": wall,
        "histories_per_s": args.histories / wall,
        "setup_s": res.timings["setup"],
        "compile_s": res.timings["compile"],
        "transport_s": res.timings["transport"],
        "host_reduce_s": res.timings["reduce"],
        "kernel_s": max((sum(p.get("chunk_seconds", [])) for p in parts), default=None),
        "peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
        "valid": res.valid,
        "counters_nonzero": res.counters.any_nonzero,
        "device": res.device,
        "register_count": parts[0].get("register_count") if parts else None,
    }
    print("#BENCH " + json.dumps(out))


def measure(backend: str, workers: int, histories: int, cache: str | None) -> dict[str, object]:
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    if cache is not None:
        env["WARP_CACHE_PATH"] = cache
    cmd = [
        sys.executable, str(Path(__file__).resolve()), "--child", "--backend", backend,
        "--child-workers", str(workers), "--histories", str(histories),
    ]  # fmt: skip
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=7200)  # noqa: S603
    lines = [ln for ln in r.stdout.splitlines() if ln.startswith("#BENCH ")]
    if r.returncode != 0 or not lines:
        raise SystemExit(f"measurement failed ({backend}, {workers}): {r.stderr[-2000:]}")
    out: dict[str, object] = json.loads(lines[-1][len("#BENCH ") :])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out")
    ap.add_argument("--histories", type=int, default=1_000_000)
    ap.add_argument("--workers", type=int, nargs="*", default=[1, 16])
    ap.add_argument("--cuda", action="store_true")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--child", action="store_true")
    ap.add_argument("--backend", default="warp-cpu")
    ap.add_argument("--child-workers", type=int, default=1)
    args = ap.parse_args()
    if args.child:
        child(args)
        return 0
    if not args.out:
        raise SystemExit("--out is required")
    out = Path(args.out).resolve()
    if not out.is_relative_to(REPO / "benchmarks" / "generated"):
        raise SystemExit("--out must lie under benchmarks/generated/")
    if out.exists():
        raise SystemExit(f"refusing to overwrite {out}")
    if args.histories < 1000 or args.repeats < 1:
        raise SystemExit("need --histories >= 1000 and --repeats >= 1")
    configs = [("warp-cpu", w) for w in args.workers] + ([("warp-cuda", 1)] if args.cuda else [])
    results = []
    for backend, workers in configs:
        with tempfile.TemporaryDirectory() as cold:  # empty kernel cache: cold start
            results.append(
                {"kind": "cold_start", **measure(backend, workers, args.histories, cold)}
            )
        for k in range(args.repeats):
            results.append(
                {"kind": f"repeat_{k + 1}", **measure(backend, workers, args.histories, None)}
            )
    sha = subprocess.run(
        ["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True
    ).stdout.strip()  # noqa: S603,S607
    doc = {
        "label": "EM-only informative aggregates; no performance targets (targets: null, V3-012)",
        "targets": None,
        "git_sha": sha or None,
        "workload": "150 MeV protons, 120x120x300 mm water, 1 mm voxels, 2 mm scoring, float32",
        "results": results,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
