"""V3-013 local validation, plan row R1-ion (``validation/plans/v3-013-acceptance.md`` A.1).

Usage::

    python validation/scripts/stopping/lv13_r1ion.py --backend {python,warp-cpu}
                                                     [--histories N] [--out-dir DIR] [--cache-dir DIR]

Repeatability C 290 MeV/u in water on the V3-008 R1-ion configuration (``steps_v8.py``
``step_r1_ion``: 60 x 60 mm water box of depth ceil(1.1 R) in 1 mm slabs, 2 mm scoring grid
``dose`` of 30 x 30 x int(1.1 R / 2) voxels, 2 mm step limit, 20 batches, float64, MCS and
straggling on, end-of-range model ``default``, nuclear off, ``track_end_positions``), with the
two selection paths of the default model:

* selection A: ``stopping_model="bethe-ionmc-v1"``;
* selection B: ``stopping=BetheStoppingSource()``.

Each selection is run twice with the same seed (plan section C, family ion13:
``20612004 + 100 r_index + k``, k = 0; r_index 11 python, 12 warp-cpu). Compared bit for bit: every
array of every grid result, the energy balance, the counters, the end states of every history, the
``EffectiveConfig.summary()`` ``tables`` block (and its sha256), and the ``stopping`` identity
block, which must name ``bethe-ionmc-v1`` for both selections and may differ only in
``selected_by``. Default histories are the V3-008 R1-ion counts (python 100, warp-cpu 10 000).
Prints ``#JSON-BEGIN`` / ``#JSON-END`` and exits non-zero on any difference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

SEED_BASE = 20612004  # ion13 family (plan section C)
R_INDEX = {"python": 11, "warp-cpu": 12}
K = 0
DEFAULT_HISTORIES = {"python": 100, "warp-cpu": 10_000}  # V3-008 R1_SPECS (steps_v8.py:373)
N_BATCHES = 20
E_U = 290.0
MAX_STEP_MM = 2.0
MODEL_ID = "bethe-ionmc-v1"
RUNS = (("A1", "A"), ("A2", "A"), ("B1", "B"), ("B2", "B"))
PAIRS = (("A1", "A2"), ("B1", "B2"), ("A1", "B1"), ("A2", "B2"))
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def seed_for(backend: str) -> int:
    return SEED_BASE + 100 * R_INDEX[backend] + K


def sha256(a: Any) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def _canon(doc: Any) -> str:
    return json.dumps(doc, sort_keys=True, default=str)


def make_config(selection: str, backend: str, n: int) -> Any:
    """The V3-008 R1-ion configuration with selection ``A`` (registered id) or ``B`` (object)."""
    from ionmc.config import DiagnosticsOptions, PhysicsOptions, RunOptions, SimulationConfig
    from ionmc.geometry import BoxPhantom
    from ionmc.materials import WATER
    from ionmc.physics.projectiles import CARBON12
    from ionmc.physics.stopping import BetheStoppingSource
    from ionmc.scoring import ScoringGrid
    from ionmc.sources import PencilBeamSource
    from ionmc.transport.tables import TransportTables

    tables = TransportTables.from_stopping_tables(
        [BetheStoppingSource().table(WATER, CARBON12)], projectile=CARBON12
    )
    range_mm = tables.range_g_cm2(0, E_U * CARBON12.a) * 10.0 / WATER.density_g_cm3
    depth = 1.1 * range_mm
    nz_box = int(math.ceil(depth / 1.0))  # water_box(depth), dz = 1 mm
    geo = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, nz_box * 1.0), WATER)
    nz = int(depth / 2)
    grid = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, nz), name="dose")
    sel: dict[str, Any] = (
        {"stopping_model": MODEL_ID} if selection == "A" else {"stopping": BetheStoppingSource()}
    )
    return SimulationConfig(
        source=PencilBeamSource(CARBON12, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), E_U * CARBON12.a),
        geometry=geo,
        scoring=(grid,),
        physics=PhysicsOptions(
            nuclear=False,
            straggling=True,
            multiple_scattering=True,
            max_step_mm=MAX_STEP_MM,
            max_energy_loss_fraction=0.02,
            **sel,
        ),
        run=RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision="float64",
            seed=seed_for(backend),
            n_histories=n,
            n_batches=N_BATCHES,
            memory_budget_bytes=2**32,
        ),
        diagnostics=DiagnosticsOptions(track_end_positions=True),
    )


def fingerprint(res: Any) -> dict[str, Any]:
    """Everything compared: arrays by sha256 of their bytes (dtype and shape included)."""
    arrays: dict[str, np.ndarray] = {}
    for g in res.grids:
        for name, v in vars(g).items():
            if isinstance(v, np.ndarray):
                arrays[f"grid.{g.name}.{name}"] = v
    for k in ("end_position_mm", "end_direction", "end_energy_mev", "end_code"):
        arrays[f"end_state.{k}"] = np.asarray(res.diagnostics[k])
    eb = res.energy_balance
    summary = res.effective_config.summary()
    return {
        "arrays": {k: [str(v.dtype), list(v.shape), sha256(v)] for k, v in sorted(arrays.items())},
        "energy_balance": {k: repr(v) for k, v in sorted(vars(eb).items())},
        "counters": res.counters.as_dict(),
        "valid": bool(res.valid),
        "tables": _canon(summary["tables"]),
        "tables_sha256": summary["tables"]["sha256"],
        "stopping": summary["stopping"],
        "effective_tables_sha256": res.effective_config.tables.sha256,
    }


def _diff(a: dict[str, Any], b: dict[str, Any]) -> list[str]:
    out = [f"array {k}" for k in a["arrays"] if a["arrays"][k] != b["arrays"].get(k)]
    out += [f"array {k} missing" for k in set(b["arrays"]) ^ set(a["arrays"])]
    for key in ("energy_balance", "counters", "valid", "tables", "tables_sha256"):
        if a[key] != b[key]:
            out.append(key)
    sa = {k: v for k, v in a["stopping"].items() if k != "selected_by"}
    sb = {k: v for k, v in b["stopping"].items() if k != "selected_by"}
    if _canon(sa) != _canon(sb):
        out.append("stopping identity (other than selected_by)")
    return out


def compare_selections(backend: str, n: int) -> dict[str, Any]:
    """Run the four runs on ``backend`` with ``n`` histories and compare them bit for bit."""
    from ionmc.simulation import Simulation
    from ionmc.transport.parity import compare_runs_bitwise

    results: dict[str, Any] = {}
    wall: dict[str, float] = {}
    for label, sel in RUNS:
        t0 = time.perf_counter()
        results[label] = Simulation(make_config(sel, backend, n)).run()
        wall[label] = time.perf_counter() - t0
    fps = {k: fingerprint(r) for k, r in results.items()}
    comparisons: dict[str, Any] = {}
    for a, b in PAIRS:
        diffs = _diff(fps[a], fps[b])
        parity = compare_runs_bitwise(results[a], results[b])
        comparisons[f"{a}-vs-{b}"] = {
            "identical": not diffs and bool(parity["identical"]),
            "differences": diffs,
            "compare_runs_bitwise": parity,
        }
    identity: dict[str, Any] = {}
    for label, sel in RUNS:
        st = fps[label]["stopping"]
        identity[label] = {
            "model_id_requested": st["model_id_requested"],
            "model_id_effective": st["model_id_effective"],
            "selected_by": st["selected_by"],
            "spec_sha256": st["spec_sha256"],
        }
    expected_by = {"A": "stopping_model", "B": "stopping-object"}
    identity_ok = all(
        identity[lb]["model_id_effective"] == MODEL_ID
        and identity[lb]["model_id_requested"] == MODEL_ID
        and identity[lb]["selected_by"] == expected_by[sel]
        for lb, sel in RUNS
    ) and identity["A1"]["selected_by"] != identity["B1"]["selected_by"]
    valid = all(
        bool(r.valid) and not r.counters.any_nonzero for r in results.values()
    )
    ok = valid and identity_ok and all(c["identical"] for c in comparisons.values())
    return {
        "backend": backend,
        "histories": n,
        "n_batches": N_BATCHES,
        "seed": seed_for(backend),
        "comparisons": comparisons,
        "identity": identity,
        "identity_ok": identity_ok,
        "valid": valid,
        "tables_sha256": {k: v["tables_sha256"] for k, v in fps.items()},
        "tally_sha256": {k: {a: v[2] for a, v in f["arrays"].items()} for k, f in fps.items()},
        "wall_s": wall,
        "pass": bool(ok),
    }


def _git_sha() -> str:
    env = os.environ.get("IONMC_RUN_SHA")
    if env:
        return env
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--backend", choices=sorted(R_INDEX), required=True)
    ap.add_argument("--histories", type=int, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--cache-dir", default=None)
    a = ap.parse_args(argv)
    if a.cache_dir is not None:
        os.environ["IONMC_CACHE_DIR"] = a.cache_dir
    n = a.histories if a.histories is not None else DEFAULT_HISTORIES[a.backend]
    t0 = time.perf_counter()
    doc = compare_selections(a.backend, n)
    doc.update(
        {
            "step": "lv13-r1ion",
            "row": "R1-ion",
            "git_sha": _git_sha(),
            "reduced": n != DEFAULT_HISTORIES[a.backend],
            "default_histories": DEFAULT_HISTORIES[a.backend],
            "total_wall_s": time.perf_counter() - t0,
        }
    )
    text = json.dumps(doc, indent=1, sort_keys=True, default=str)
    if a.out_dir is not None:
        a.out_dir.mkdir(parents=True, exist_ok=True)
        (a.out_dir / f"lv13-r1ion-{a.backend}.json").write_text(text + "\n")
    print("#JSON-BEGIN")
    print(text)
    print("#JSON-END")
    return 0 if doc["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
