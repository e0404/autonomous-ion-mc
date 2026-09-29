"""JSON run configuration for the command line (``ionmc run``).

Schema (all lengths mm, energies MeV/u):

```json
{
  "source": {"species": "proton", "energy_mev_per_u": 150, "position_mm": [0, 0, -1],
             "direction": [0, 0, 1], "sigma_mm": 0, "sigma_energy_fraction": 0, "sigma_angle_rad": 0},
  "geometry": {"type": "box", "size_mm": [120, 120, 300], "spacing_mm": 1, "material": "water",
               "slabs": [[z0, z1, "bone_cortical", null]]},
  "scoring": {"spacing_mm": 2} | {"origin_mm": [...], "spacing_mm": [...], "shape": [...]},
  "histories": 10000, "batches": 10, "seed": 1, "backend": "python", "precision": "float64",
  "physics": {"straggling": true, "multiple_scattering": true, "nuclear": false, ...},
  "scorers": ["energy", "dose"]
}
```
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ionmc.config import PhysicsConfig, SimulationConfig
from ionmc.geometry import VoxelGeometry, homogeneous_box, slab_phantom, vec3, vec3i
from ionmc.scoring import ScoringGrid
from ionmc.sources import PencilBeam


def geometry_from_dict(spec: dict[str, Any]) -> VoxelGeometry:
    kind = spec.get("type", "box")
    if kind != "box":
        raise ValueError(f"unsupported geometry type {kind!r}; supported: 'box'")
    size = vec3(spec["size_mm"])
    spacing = spec.get("spacing_mm", 1.0)
    material = spec.get("material", "water")
    slabs = spec.get("slabs")
    if slabs:
        return slab_phantom(
            size,
            float(spacing),
            [
                (float(a), float(b), m, None if d is None else float(d))
                for a, b, m, d in slabs
            ],
            material,
        )
    return homogeneous_box(
        size, spacing, material, spec.get("origin_mm"), spec.get("density_g_cm3")
    )


def scoring_from_dict(
    spec: dict[str, Any] | None, geometry: VoxelGeometry
) -> ScoringGrid | None:
    if not spec:
        return None
    if "shape" in spec:
        return ScoringGrid(
            vec3(spec["origin_mm"]),
            vec3(spec["spacing_mm"]),
            vec3i(spec["shape"]),
        )
    grid = ScoringGrid.coarse(
        geometry,
        spec["spacing_mm"]
        if not isinstance(spec["spacing_mm"], list)
        else tuple(spec["spacing_mm"]),
    )
    if "shift_mm" in spec:
        grid = grid.shifted(vec3(spec["shift_mm"]))
    return grid


def config_from_dict(spec: dict[str, Any]) -> SimulationConfig:
    src = spec["source"]
    source = PencilBeam(
        src["species"],
        float(src["energy_mev_per_u"]),
        vec3(src.get("position_mm", (0.0, 0.0, -1.0))),
        vec3(src.get("direction", (0.0, 0.0, 1.0))),
        float(src.get("sigma_mm", 0.0)),
        float(src.get("sigma_energy_fraction", 0.0)),
        float(src.get("sigma_angle_rad", 0.0)),
        src.get("name", "pencil-beam"),
    )
    geometry = geometry_from_dict(spec["geometry"])
    physics = PhysicsConfig(**spec.get("physics", {}))
    return SimulationConfig(
        source,
        geometry,
        int(spec["histories"]),
        int(spec.get("batches", 10)),
        int(spec.get("seed", 12345)),
        spec.get("backend", "python"),
        spec.get("precision", "float64"),
        scoring_from_dict(spec.get("scoring"), geometry),
        physics,
        tuple(spec.get("scorers", ("energy", "dose"))),
    )


def load_config(path: str | Path) -> SimulationConfig:
    return config_from_dict(json.loads(Path(path).read_text()))
