"""Fixed pencil configurations of plan v3-007 row REG-2 (pencil API only, so the same module runs
against the 3f3736a6 tree that the golden digests were computed from)."""

from __future__ import annotations

import hashlib
import json
from typing import Any

from ionmc.config import DiagnosticsOptions, PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.configio import to_dict
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import ALPHA, PROTON
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.sources import PencilBeamSource

BETHE = "bethe-ionmc-v1"


def _cfg(species: Any = PROTON, energy: float = 20.0, sigma: float = 0.0,
         direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
         physics: dict[str, Any] | None = None, run: dict[str, Any] | None = None,
         tallies: tuple[TallyRequest, ...] = ()) -> SimulationConfig:  # fmt: skip
    ph = {"nuclear": False, "stopping_model": BETHE, "multiple_scattering": species.z < 2}
    rn = {"backend": "python", "precision": "float64", "seed": 5, "n_histories": 40,
          "n_batches": 4}  # fmt: skip
    return SimulationConfig(
        PencilBeamSource(species, (0.0, 0.0, -5.0), direction, energy, sigma),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 5.0), (4, 4, 4)),),
        PhysicsOptions(**(ph | (physics or {}))),
        RunOptions(**(rn | (run or {}))),
        DiagnosticsOptions(),
        tallies=tallies,
    )


def fixed_configs() -> dict[str, SimulationConfig]:
    let = (TallyRequest("let", "dose", "let_d"),)
    return {
        "proton-let": _cfg(tallies=let, direction=(0.0, 0.0, 2.0), sigma=0.1,
                           physics={"end_of_range_model": "point_v1"}, run={"max_steps": 1000}),
        "proton-warp": _cfg(run={"backend": "warp-cpu", "precision": "float32"}),
        "proton-no-straggling": _cfg(physics={"straggling": False, "multiple_scattering": False,
                                              "stopping_model": BETHE + "+ls-eff"}),
        "alpha-mcs-off": _cfg(ALPHA, 400.0, physics={"multiple_scattering": False}),
        "alpha-mcs-on": _cfg(ALPHA, 400.0, sigma=1.0, physics={"multiple_scattering": True}),
    }  # fmt: skip


def digests() -> dict[str, tuple[str, str]]:
    """name -> (sha256 of the canonical JSON of summary(), of to_dict())."""
    out = {}
    for name, c in fixed_configs().items():
        s = json.dumps(validate(c).summary(), sort_keys=True, separators=(",", ":"))
        d = json.dumps(to_dict(c), sort_keys=True, separators=(",", ":"))
        out[name] = (hashlib.sha256(s.encode()).hexdigest(), hashlib.sha256(d.encode()).hexdigest())
    return out
