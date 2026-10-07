"""Digests of the qualified outputs of fixed configurations (row A16 of
``validation/plans/v3-004-acceptance.md``).

Standalone and deliberately restricted to the public API that exists at the baseline commit
the qualified-path baseline commit (``A16_BASELINE`` of ``steps_v4.py``, a524f209 at V3-003D; ``SimulationConfig`` without scoring channels), so that ``steps_v4.py a16`` can run
it twice against two source trees (the baseline extracted with ``git archive`` and the tree under
test) with ``PYTHONPATH`` set accordingly. ``--tallies all`` adds scoring channels of every
quantity kind to the configuration (only possible with the tree under test).

Usage::

    python a16_digest.py --specs t1:python:float64,t13:warp-cpu:float32 [--tallies all]

Configurations: ``t1`` (the V3-003 T1 geometry: 12 x 12 x 40 voxels of 5 mm, 2 mm scoring grid, 100
MeV, 32 histories in 2 batches, 2 mm steps, end positions and a trace of every history) and ``t13``
(the T13 layout: a 60 x 60 x 1.1 R water box, 2 mm scoring grid, 150 MeV, 20 batches, 1 mm steps,
chunks of 2^10 histories on the Warp backends; 200 histories on python, 20000 on Warp, 4 traced
histories). Everything uses the analytic Bethe source. The output is one JSON document
``{"source": <ionmc.__file__>, "digests": {spec: {field: sha256}}}``; the fields are the per-grid
batch energy arrays (exactly the int64 quanta), every field of the energy balance, the counters and
every diagnostic array (end state and trace).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

import ionmc
from ionmc.config import DiagnosticsOptions, PhysicsOptions, RunOptions, SimulationConfig
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource

REPO = Path(__file__).resolve().parents[3]
SEED = 20351004  # fixed: the digests are compared between two trees, not a statistical result


def _config(name: str, backend: str, precision: str, tallies: str) -> SimulationConfig:
    stop = BetheStoppingSource()
    if name == "t1":
        shape = (12, 12, 40)
        geo: Any = VoxelGeometry(
            (-30.0, -30.0, 0.0), (5.0, 5.0, 5.0), shape, (WATER,), np.zeros(shape, dtype=np.int32)
        )
        grid = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 100), name="dose")
        energy, n, nb, step, chunk, trace = 100.0, 32, 2, 2.0, None, 32
    elif name == "t13":
        energy = 150.0
        depth = 1.1 * 158.6
        geo = BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, depth), WATER)
        grid = ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, int(depth / 2)),
                           name="dose")  # fmt: skip
        n = 200 if backend == "python" else 20000
        nb, step, chunk, trace = 20, 1.0, (None if backend == "python" else 2**10), 4
    else:
        raise SystemExit(f"unknown configuration {name!r}")
    if precision != "float64":
        trace = 0  # the per-step trace exists in float64 only
    extra: dict[str, Any] = {} if chunk is None else {"chunk_histories": chunk}
    kw: dict[str, Any] = {}
    if tallies == "all":
        kw = _channel_requests()
    return SimulationConfig(
        source=PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), energy),
        geometry=geo,
        scoring=(grid,),
        physics=PhysicsOptions(nuclear=False, stopping=stop, max_step_mm=step),
        run=RunOptions(  # type: ignore[arg-type]
            backend=backend, precision=precision, seed=SEED, n_histories=n, n_batches=nb, **extra
        ),
        diagnostics=DiagnosticsOptions(track_end_positions=True, trace_histories=trace),
        **kw,
    )


def _channel_requests() -> dict[str, Any]:
    from ionmc.lookup import LookupTable
    from ionmc.scoring import TallyRequest

    lk = LookupTable.from_file(REPO / "tests" / "data" / "synthetic_lookup.json")
    edges = tuple(np.geomspace(2.0, 110.0, 9))
    tallies = (
        TallyRequest("ed", "dose", "edep"),
        TallyRequest("lt", "dose", "let_t"),
        TallyRequest("ld", "dose", "let_d"),
        TallyRequest("le", "dose", "let_d_eps"),
        TallyRequest("fl", "dose", "fluence"),
        TallyRequest("fe", "dose", "lookup_sum", lookup=lk.name),
        TallyRequest("sp", "dose", "fluence_spectrum", energy_edges_mev_per_u=edges),
    )
    return {"tallies": tallies, "lookups": (lk,)}


def _hash(x: Any) -> str:
    h = hashlib.sha256()
    if isinstance(x, np.ndarray):
        h.update(str(x.dtype).encode() + str(x.shape).encode() + np.ascontiguousarray(x).tobytes())
    else:
        h.update(repr(x).encode())
    return h.hexdigest()


RAW: dict[str, Any] = {}
"""Raw values of the last digested fields (``field -> array or scalar``), for ``--dump``."""


def _put(key: str, obj: Any, out: dict[str, str]) -> None:
    out[key] = _hash(obj)
    RAW[key] = obj


def _flatten(prefix: str, obj: Any, out: dict[str, str]) -> None:
    if isinstance(obj, dict):
        for k in sorted(obj, key=str):
            _flatten(f"{prefix}.{k}", obj[k], out)
    elif isinstance(obj, np.ndarray):
        _put(prefix, obj, out)
    elif isinstance(obj, tuple | list) and obj and isinstance(obj[0], np.ndarray):
        for i, v in enumerate(obj):
            _flatten(f"{prefix}[{i}]", v, out)
    else:
        _put(prefix, obj, out)


def digest(spec: str, tallies: str) -> dict[str, str]:
    name, backend, precision = spec.split(":")
    res = Simulation(_config(name, backend, precision, tallies)).run()
    out: dict[str, str] = {}
    RAW.clear()
    for g in res.grids:
        _put(f"grid.{g.name}.batch_energy_mev", np.asarray(g.batch_energy_mev), out)
    _flatten("energy_balance", asdict(res.energy_balance), out)
    _flatten("counters", res.counters.as_dict(), out)
    _flatten("diagnostics", res.diagnostics, out)
    _put("valid", res.valid, out)
    return out


def table_identity() -> dict[str, Any]:
    """Identity fields of the analytic water table of the tree under test (``None`` = absent, as
    at the V3-003D baseline, whose tables carry no ``range_construction``); verified by
    ``steps_v4.py a16 --mode intended-change``."""
    meta = BetheStoppingSource().table(WATER, PROTON).metadata
    return {"range_construction": meta.get("range_construction")}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--specs", required=True, help="name:backend:precision,...")
    ap.add_argument("--tallies", choices=("none", "all"), default="none")
    ap.add_argument("--dump", default=None, help="write the raw values to this .npz (A16 bounds)")
    args = ap.parse_args(argv)
    raw: dict[str, Any] = {}
    digests: dict[str, dict[str, str]] = {}
    for sp in args.specs.split(","):
        digests[sp] = digest(sp, args.tallies)
        for k, v in RAW.items():
            raw[f"{sp}|{k}"] = np.asarray(v) if not isinstance(v, str) else np.asarray(repr(v))
    if args.dump:
        np.savez(args.dump, **raw)
    doc = {
        "source": str(Path(ionmc.__file__).resolve()),
        "tallies": args.tallies,
        "table_identity": table_identity(),
        "digests": digests,
    }
    print("#DIGEST-BEGIN")
    print(json.dumps(doc, sort_keys=True))
    print("#DIGEST-END")
    return 0


if __name__ == "__main__":
    sys.exit(main())
