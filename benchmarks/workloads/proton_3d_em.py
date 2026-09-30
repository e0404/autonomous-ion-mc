#!/usr/bin/env python3
"""Timing of the proton-3d workload envelope with the physics implemented so far.

NOT a performance-target measurement: nuclear interactions, secondaries and LET
scoring are not implemented yet, so the numbers only track the electromagnetic
transport cost on each backend. Geometry, source and scoring follow
``experiment/v2/performance.json`` (150 MeV protons, 120 x 120 x 300 mm water,
1 mm transport voxels, 2 mm scoring voxels, 10 batches).

Usage: python benchmarks/workloads/proton_3d_em.py --backend warp-cuda \
           --histories 100000 1000000 --output benchmarks/generated/proton-3d-em.json
"""

from __future__ import annotations

import argparse
import json
import platform
import resource
import time
from pathlib import Path

from ionmc.config import SimulationConfig
from ionmc.geometry import homogeneous_box
from ionmc.provenance import code_identity
from ionmc.scoring import ScoringGrid
from ionmc.simulation import build_tables, run
from ionmc.sources import PencilBeam


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--backend", default="warp-cuda")
    p.add_argument("--precision", default="float32")
    p.add_argument("--histories", nargs="+", type=int, default=[100000])
    p.add_argument("--batches", type=int, default=10)
    p.add_argument("--repeats", type=int, default=2)
    p.add_argument("--output", type=Path)
    a = p.parse_args(argv)
    t_import = time.perf_counter()
    geo = homogeneous_box((120.0, 120.0, 300.0), 1.0, "water")
    beam = PencilBeam("proton", 150.0, (0.0, 0.0, -1.0))
    grid = ScoringGrid.coarse(geo, 2.0)
    tables = build_tables(SimulationConfig(beam, geo, 1, 1), offline=True)
    setup_seconds = time.perf_counter() - t_import
    report = {
        "schema_version": 1,
        "purpose": "electromagnetic-only timing of the proton-3d envelope; not a target measurement",
        "code": code_identity(),
        "platform": platform.platform(),
        "backend": a.backend,
        "precision": a.precision,
        "physics": "energy loss, Gamma/Bohr straggling, differential-Moliere MCS; nuclear=False",
        "setup_seconds": setup_seconds,
        "runs": [],
    }
    for n in a.histories:
        for r in range(a.repeats):
            cfg = SimulationConfig(
                beam,
                geo,
                histories=n,
                batches=a.batches,
                seed=1000 + r,
                backend=a.backend,
                precision=a.precision,
                scoring=grid,
            )
            t0 = time.perf_counter()
            res = run(cfg, tables=tables, offline=True)
            wall = time.perf_counter() - t0
            tm = res.metadata["timing"]
            acc = res.metadata["energy_accounting_per_primary"]
            report["runs"].append(
                {
                    "histories": n,
                    "repeat": r,
                    "wall_seconds": wall,
                    "transport_wall_seconds": tm["wall_seconds"],
                    "kernel_seconds": tm.get("kernel_seconds"),
                    "host_seconds": tm.get("host_seconds"),
                    "first_launch_seconds": tm.get("first_launch_seconds"),
                    "histories_per_second_transport": n / tm["wall_seconds"],
                    "steps_per_history": acc["steps_per_history"],
                    "escaped_mev_per_primary": acc["escaped_mev"],
                    "max_dose_gy_per_primary": float(
                        res.dose_mean[~__import__("numpy").isnan(res.dose_mean)].max()
                    ),
                    "host_peak_rss_mib": resource.getrusage(
                        resource.RUSAGE_SELF
                    ).ru_maxrss
                    / 1024.0,
                }
            )
            print(json.dumps(report["runs"][-1]), flush=True)
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
