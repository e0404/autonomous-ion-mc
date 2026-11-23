"""Generator of the committed golden fixture ``golden/`` (plan V3-006 rows REP-2 and REG-4).

Run from the repository root with the project environment::

    python tests/ionmc/data/v3_006/gen_golden.py

It builds ONE synthetic, KB-scale ``ionmc.result/1`` bundle with ``ionmc.results.save`` from
hand-made arrays (a real ``EffectiveConfig`` of a tiny tilted-source configuration; no transport,
no table value, no stopping data) and stores it as ``golden/``: ``result.json``, ``report.md``,
``plots/*.svg`` and ``arrays/<member>.npy`` (one small file per catalogue array). The compressed
``arrays.npz`` is NOT committed (REG-4: no committed ``.npz``); the test assembles it
deterministically with the writer of ``ionmc.results`` (``write_npz``) from the ``.npy`` files.
The bundle is written with ``zlib_level = 0`` so that the assembled archive is stored, not
deflated, and its bytes do not depend on the zlib version. Code identity and environment are fixed
synthetic records (a snapshot), so regeneration on another machine gives the same bytes.
"""

from __future__ import annotations

import dataclasses
import json
import os
import shutil
import tempfile
import zipfile
from pathlib import Path
from typing import Any

import numpy as np

from ionmc import provenance as prov
from ionmc import results
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.scoring import QuantityResult, ScoringGrid, TallyRequest
from ionmc.simulation import EnergyBalance, Result, TransportCounters, _grid_result
from ionmc.sources import PencilBeamSource

HERE = Path(__file__).resolve().parent
OUT = HERE / "golden"
SHAPE = (2, 8, 2)
NB = 4
N_HIST = 8
TREE_OID = "0123456789abcdef0123456789abcdef01234567"
CODE = {
    "package_tree_oid": TREE_OID,
    "ionmc_version": "0.1.0.dev0",
    "git_commit": None,
    "git_dirty": None,
    "source": "unavailable",
    "status": "unknown",
    "reason": "synthetic golden fixture: not produced by a code version",
}
ENVIRONMENT: dict[str, Any] = {
    "python_version": "3.12.0",
    "platform": "synthetic-golden-fixture",
    "numpy_version": "0",
    "warp_version": None,
    "warp_devices": [],
    "warp_error": None,
}


def config() -> SimulationConfig:
    """A tilted proton pencil beam: the largest direction component is the y axis (axis 1)."""
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, -12.0, 2.0), (0.2, 0.95, 0.1), 20.0),
        # the last y layer of the grid lies outside the box: its voxel mass is 0 (dose undefined)
        BoxPhantom((-2.0, -8.0, 0.0), (4.0, 14.0, 4.0), WATER),
        (ScoringGrid((-2.0, -8.0, 0.0), (2.0, 2.0, 2.0), SHAPE, "dose"),),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1"),
        RunOptions("python", "float64", 11, N_HIST, n_batches=NB),
        tallies=(TallyRequest("let", "dose", "let_d"),),
    )


def _quantity() -> QuantityResult:
    n = int(np.prod(SHAPE))
    mean: Any = (1.0 + 0.1 * np.arange(n)).reshape(SHAPE)
    mask = np.ones(SHAPE, dtype=bool)
    mask[:, 0, :] = False  # undefined LET in the first y layer: a gap in the plot
    mean = np.where(mask, mean, np.nan)
    return QuantityResult(
        "let", "let_d", "ratio", "MeV/mm", "sum(l S^2)/sum(l S), S in water",
        mean, 0.05 * mean, mask,
        np.where(mask, 4, 0).astype(np.int32), np.where(mask, 1e-9, np.nan), (3, 3),
    )  # fmt: skip


def _pieces() -> dict[str, QuantityResult]:
    spec = (
        ("edep_excluded_from_let", "edep", "MeV", results._EXCLUDED_DEFINITION),
        ("scoring_pieces", "count", "1", results._PIECE_DEFINITIONS["scoring_pieces"]),
        ("scoring_pieces_local", "count", "1", results._PIECE_DEFINITIONS["scoring_pieces_local"]),
    )
    out = {}
    for n, q, u, d in spec:
        one = np.ones(SHAPE)
        out[n] = QuantityResult(
            n, q, "linear", u, d, one, 0.0 * one, np.ones(SHAPE, dtype=bool),
            np.full(SHAPE, 4, dtype=np.int32), 0.0 * one, (3,),
        )  # fmt: skip
    return out


def synthetic_result() -> Result:
    cfg = config()
    eff = validate(cfg)
    # synthetic per-batch sums [NB, nvox]: a Bragg-like rise along y, noise from a fixed seed
    rng = np.random.default_rng(20261123)
    y = (np.arange(SHAPE[1]) + 0.5) / SHAPE[1]
    shape_y = np.exp(-(((y - 0.7) / 0.2) ** 2)) + 0.1
    per_vox = np.broadcast_to(shape_y[None, :, None], SHAPE).reshape(-1)
    sums = rng.uniform(0.8, 1.2, (NB, per_vox.size)) * per_vox[None, :] * 4.0
    sums[:, 0] = 0.0  # one voxel with no deposit: undefined
    g = _grid_result(cfg.scoring[0], sums, eff, True)
    g = dataclasses.replace(g, quantities={"let": _quantity(), **_pieces()})
    in_grid = float(np.sum(g.energy_mev)) * N_HIST
    total = in_grid + 1.0
    bal = EnergyBalance(
        total + 0.25, total - 0.5, 0.5, 0.25, 0.0, 0.0, (in_grid,), (total - in_grid,)
    )
    env = {
        **ENVIRONMENT,
        "ionmc_version": CODE["ionmc_version"],
        "provenance_snapshot": {"code": dict(CODE), "environment": dict(ENVIRONMENT)},
    }
    return Result(
        True, cfg, eff, "python", "cpu", "float64", 11, dict(eff.rng), N_HIST, NB,
        (g,), bal, TransportCounters(), {"total": 0.1}, env, {},
    )  # fmt: skip


def build(out: Path = OUT) -> Path:
    """Write the committed golden directory ``out`` (replacing it)."""
    old_oid = prov.package_tree_oid
    prov.package_tree_oid = lambda *a, **k: TREE_OID
    cwd = Path.cwd()
    try:
        with tempfile.TemporaryDirectory() as td:
            os.chdir(td)
            bundle = results.save(synthetic_result(), "golden", zlib_level=0)  # relative path
            if out.exists():
                shutil.rmtree(out)
            out.mkdir(parents=True)
            for name in ("result.json", "report.md"):
                shutil.copyfile(bundle / name, out / name)
            if (bundle / "plots").is_dir():
                shutil.copytree(bundle / "plots", out / "plots")
            (out / "arrays").mkdir()
            with zipfile.ZipFile(bundle / "arrays.npz") as zf:
                for info in zf.infolist():
                    (out / "arrays" / info.filename).write_bytes(zf.read(info))
    finally:
        os.chdir(cwd)
        prov.package_tree_oid = old_oid
    doc: Any = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert doc["run"]["output_path"] == "golden"
    return out


if __name__ == "__main__":
    print(build())
