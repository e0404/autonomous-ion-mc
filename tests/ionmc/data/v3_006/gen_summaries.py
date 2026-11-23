"""Generate ``summaries-95937553.json`` (row REG-2 of plan v3-006).

``EffectiveConfig.summary()`` of a fixed set of configurations as produced at the V3-006 base
commit 95937553 (V3-013 integrated), before any V3-006 code. The fixture holds identities and
digests only (stopping and transport table sha256, nuclear table id), never table values. The
only normalisation is the value of every ``cache_dir`` key (a machine path), replaced by
``"<cache_dir>"``; the test applies the same function.

Run once at the base: ``IONMC_CACHE_DIR=<cache> python gen_summaries.py`` (single process). The
cases that need a cache (``icru90-water-v1``, the nuclear table) are written only when the cache
is present; the test skips a case whose cache is missing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.lookup import LookupTable
from ionmc.materials import PMMA, WATER
from ionmc.physics.projectiles import ALPHA, CARBON12, PROTON, Projectile
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.sources import PencilBeamSource

HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "summaries-95937553.json"
LOOKUP = HERE.parents[2] / "data" / "synthetic_lookup.json"
CACHE_DIR_PLACEHOLDER = "<cache_dir>"
SEED = 5


def normalise(obj: Any) -> Any:
    """Replace the value of every ``cache_dir`` key by the placeholder."""
    if isinstance(obj, dict):
        return {
            k: CACHE_DIR_PLACEHOLDER if k == "cache_dir" and v is not None else normalise(v)
            for k, v in obj.items()
        }
    if isinstance(obj, list):
        return [normalise(x) for x in obj]
    return obj


def dumps(obj: Any) -> str:
    """The text form compared byte for byte (sorted keys, two-space indent, final newline)."""
    return json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n"


def _box_config(
    projectile: Projectile,
    energy: float,
    model: str,
    backend: str,
    precision: str,
    *,
    nuclear: bool = False,
    nuclear_table_id: str | None = None,
) -> SimulationConfig:
    grid = ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 2.0), (4, 4, 20))
    return SimulationConfig(
        PencilBeamSource(projectile, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), energy),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), WATER),
        (grid,),
        PhysicsOptions(nuclear=nuclear, stopping_model=model, nuclear_table_id=nuclear_table_id),
        RunOptions(backend, precision, SEED, 8, n_batches=4),  # type: ignore[arg-type]
    )


def _slab_config() -> SimulationConfig:
    import numpy as np

    index = np.zeros((3, 3, 8), dtype=np.int32)
    index[:, :, 5:] = 1
    geometry = VoxelGeometry((-6.0, -6.0, 0.0), (4.0, 4.0, 4.0), (3, 3, 8), (WATER, PMMA), index)
    grid = ScoringGrid((-6.0, -6.0, 0.0), (4.0, 4.0, 4.0), (3, 3, 8))
    table = LookupTable.from_file(LOOKUP)
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 30.0),
        geometry,
        (grid,),
        PhysicsOptions(nuclear=False, stopping_model="bethe-ionmc-v1"),
        RunOptions("warp-cpu", "float64", SEED, 8, n_batches=4),
        tallies=(
            TallyRequest("let", "dose", "let_d"),
            TallyRequest("weighted", "dose", "lookup_sum", lookup=table.name),
        ),
        lookups=(table,),
    )


def nuclear_table_id() -> str | None:
    """The nuclear table the fixture uses: the full build with the largest id of the cache, as
    ``tests/ionmc/test_nuclear_transport_config.py`` chooses it; None without a cache."""
    from ionmc.data import cache
    from ionmc.nuclear.tables import NuclearTable, NuclearTableError

    found: list[tuple[bool, str]] = []
    for p in sorted((cache.resolve_cache_dir(None) / "derived").glob("nuclear-proton-*.json")):
        tid = p.stem.removeprefix("nuclear-proton-")
        try:
            info = NuclearTable.load(None, tid).info
        except NuclearTableError:
            continue
        if "transport_path_bound_terms" in info:
            found.append((info["options"]["diagnostic_nodes_mev"] is None, tid))
    return max(found)[1] if found else None


def build_case(name: str, params: dict[str, Any]) -> SimulationConfig:
    """The configuration of case ``name`` (``params`` carries the nuclear table id)."""
    bethe = "bethe-ionmc-v1"
    if name == "proton_box_python":
        return _box_config(PROTON, 100.0, bethe, "python", "float64")
    if name == "proton_box_warp_cpu_float32":
        return _box_config(PROTON, 100.0, bethe, "warp-cpu", "float32")
    if name == "voxel_slab_let_lookup_warp_cpu":
        return _slab_config()
    if name == "carbon_box_warp_cpu":
        return _box_config(CARBON12, 1200.0, bethe, "warp-cpu", "float64")
    if name == "nuclear_proton_python":
        return _box_config(
            PROTON,
            100.0,
            bethe,
            "python",
            "float64",
            nuclear=True,
            nuclear_table_id=params["nuclear_table_id"],
        )
    if name == "alpha_box_icru90_warp_cpu":
        return _box_config(ALPHA, 400.0, "icru90-water-v1", "warp-cpu", "float64")
    raise KeyError(name)


CASES = {
    "proton_box_python": {"requires": []},
    "proton_box_warp_cpu_float32": {"requires": []},
    "voxel_slab_let_lookup_warp_cpu": {"requires": []},
    "carbon_box_warp_cpu": {"requires": []},
    "nuclear_proton_python": {"requires": ["nuclear table"]},
    "alpha_box_icru90_warp_cpu": {"requires": ["icru90-water-v1"]},
}


def summary_of(name: str, params: dict[str, Any]) -> dict[str, Any]:
    """The normalised ``validate(config).summary()`` of case ``name``."""
    return normalise(validate(build_case(name, params)).summary())  # type: ignore[no-any-return]


def main() -> int:
    tid = nuclear_table_id()
    cases: dict[str, Any] = {}
    for name, meta in CASES.items():
        params = {"nuclear_table_id": tid} if name == "nuclear_proton_python" else {}
        if name == "nuclear_proton_python" and tid is None:
            print(f"skipped {name}: no nuclear table in the cache", file=sys.stderr)
            continue
        cases[name] = {
            "requires": meta["requires"],
            "params": params,
            "summary": summary_of(name, params),
        }
    doc = {
        "base_commit": "95937553",
        "normalisation": "the value of every cache_dir key is replaced by '<cache_dir>'",
        "cases": cases,
    }
    FIXTURE.write_text(dumps(doc), encoding="utf-8")
    print(f"wrote {FIXTURE} ({len(cases)} cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
