"""The truncated-hinge diagnostic (T14 negative control): off by default and then
bit-for-bit the engine of the stored baseline trace; on, python and warp-cpu float64 still
follow the same trajectory (T1) and the trajectory differs from the default."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
)
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import VoxelGeometry
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.simulation import Simulation
from ionmc.sources import PencilBeamSource
from ionmc.transport.parity import compare_traces
from ionmc.transport.tally import TRACE_COLUMNS

BASELINE = Path(__file__).parent / "data" / "trace_baseline_20mev.npz"


def _config(backend: str = "python", *, diagnostic: bool = False) -> SimulationConfig:
    shape = (16, 16, 16)
    geo = VoxelGeometry(
        (-12.0, -12.0, 0.0), (1.5, 1.5, 1.5), shape, (WATER,), np.zeros(shape, dtype=np.int32)
    )
    return SimulationConfig(
        PencilBeamSource(PROTON, (0.0, 0.0, 0.0), (0.1, 0.05, 1.0), 20.0, 0.0, 0.3),
        geo,
        (ScoringGrid((-12.0, -12.0, 0.0), (2.0, 2.0, 2.0), (12, 12, 12)),),
        PhysicsOptions(
            nuclear=False,
            stopping=BetheStoppingSource(),
            max_step_mm=1.0,
            truncated_hinge_diagnostic=diagnostic,
        ),
        RunOptions(
            backend=backend,  # type: ignore[arg-type]
            precision="float64",
            seed=7,
            n_histories=4,
            n_batches=2,
        ),
        DiagnosticsOptions(trace_histories=4),
    )


def test_default_trace_is_unchanged_from_the_stored_baseline() -> None:
    """The baseline trace was produced before the diagnostic existed (commit 6d58480): with the
    flag off the python and warp-cpu float64 traces still reproduce it (rtol 1e-9 across
    platforms' libm; on the generating host they are bit-identical)."""
    base = np.load(BASELINE)
    for backend in ("python", "warp-cpu"):
        res = Simulation(_config(backend)).run()
        tr = res.diagnostics["trace"]
        for name in TRACE_COLUMNS:
            if name in ("history", "step", "ix", "iy", "iz", "reason", "blocks", "attempts"):
                assert np.array_equal(tr[name], base[name]), (backend, name)
            else:
                np.testing.assert_allclose(tr[name], base[name], rtol=1e-9, atol=1e-12)


def test_diagnostic_on_keeps_python_warp_parity_and_changes_trajectories() -> None:
    py = Simulation(_config("python", diagnostic=True)).run()
    wp = Simulation(_config("warp-cpu", diagnostic=True)).run()
    v = compare_traces(py.diagnostics, wp.diagnostics)
    assert v["pass"], v
    base = np.load(BASELINE)
    ux = py.diagnostics["trace"]["ux"]
    assert len(ux) != len(base["ux"]) or not np.allclose(ux, base["ux"], rtol=0, atol=1e-12)
    assert py.effective_config.summary()["physics"]["truncated_hinge_diagnostic"] is True
    default = Simulation(_config("python")).run()
    assert default.effective_config.summary()["physics"]["truncated_hinge_diagnostic"] is False


def test_diagnostic_flag_is_validated() -> None:
    cfg = _config()
    with pytest.raises(UnsupportedCombinationError, match="bool"):
        replace(cfg, physics=replace(cfg.physics, truncated_hinge_diagnostic=1))  # type: ignore[arg-type]
