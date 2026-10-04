"""Worker and memory planning of the validation steps: with ``--workers auto`` on a many-core host
no step builds an invalid configuration (the per-worker private accumulators stay within the
declared memory budget), and CUDA configurations never carry CPU workers."""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from ionmc.config import SimulationConfig, validate
from ionmc.scoring import ScoringGrid

SCRIPTS = Path(__file__).resolve().parents[2] / "validation" / "scripts" / "transport"


def _load(name: str) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _Stop(Exception):
    pass


@pytest.fixture
def steps(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    monkeypatch.setattr(os, "cpu_count", lambda: 64)  # a 64-core host
    monkeypatch.setattr("ionmc.config.cuda_available", lambda: True)  # configurations only
    return _load("steps")


def test_plan_workers_rules(steps: ModuleType) -> None:
    grid = (ScoringGrid((0, 0, 0), (1, 1, 1), (100, 100, 100)),)  # 1e6 voxels
    big = steps.plan_workers(64, "warp-cpu", 100, grid, 10**6, budget=2**31)
    assert big["accumulator_bytes_per_worker"] == 100 * 10**6 * 8
    assert big["workers_used"] == 2**31 // (100 * 10**6 * 8) == 2
    assert steps.plan_workers(64, "warp-cpu", 2, grid, 10**6, budget=2**40)["workers_used"] == 64
    assert steps.plan_workers(64, "warp-cpu", 100, grid, 10**6, budget=1)["workers_used"] == 1
    assert steps.plan_workers(64, "warp-cuda", 2, grid, 10**6)["workers_used"] == 1
    assert steps.plan_workers(64, "python", 2, grid, 10)["workers_used"] == 10  # <= histories


def test_every_t12_sample_builds_a_valid_configuration_on_64_cores(steps: ModuleType) -> None:
    import argparse

    a = argparse.Namespace(energy=150.0, lateral_bin=0.2, half_width=20.0, scale=1.0, seed=7,
                           timeout=None)  # fmt: skip
    workers = steps_workers = 64
    for name in steps.SAMPLES:
        cfg, _ = steps.sample_config(a, name, workers)
        validate(cfg)  # raises on an invalid configuration
        used = cfg.run.cpu_workers
        if name == "cuda32":
            assert used == 1 and cfg.run.backend == "warp-cuda"  # never CPU workers on CUDA
        else:
            assert 1 <= used <= steps_workers
        per = cfg.run.n_batches * sum(g.n_voxels for g in cfg.scoring) * 8
        assert used * per <= cfg.run.memory_budget_bytes


@pytest.mark.parametrize(
    "argv",
    [
        ["t1", "--k", "4"],
        ["t2"],
        ["t2", "--backend", "warp-cuda"],
        ["t13", "--mode", "workers", "--n", "2000"],
        ["t13", "--mode", "chunks", "--backend", "warp-cuda", "--n", "2000"],
        ["t-r1", "--runs", "warp-cuda:float32:2000"],
        ["t-r1", "--runs", "python:float64:400,warp-cpu:float64:2000"],
        ["t8", "--n", "2000"],
        ["t9", "--n", "2000"],
        ["t10", "--n", "2000"],
        ["t14", "--n", "2000"],
    ],
    ids=lambda a: " ".join(a),
)
def test_steps_build_valid_configurations_with_workers_auto(
    steps: ModuleType, monkeypatch: pytest.MonkeyPatch, argv: list[str]
) -> None:
    """The first configuration of every step validates with ``--workers auto`` on a mocked 64-core
    host (a step that forwarded 64 CPU workers to CUDA, or ignored the accumulator budget, would
    raise ``UnsupportedCombinationError`` here), and carries a recorded plan."""
    seen: list[SimulationConfig] = []

    class FakeSimulation:
        def __init__(self, cfg: SimulationConfig) -> None:
            validate(cfg)
            seen.append(cfg)

        def run(self) -> Any:
            raise _Stop

    monkeypatch.setattr(steps, "Simulation", FakeSimulation)
    with pytest.raises(_Stop):
        steps.main([*argv, "--workers", "auto"])
    cfg = seen[0]
    if cfg.run.backend == "warp-cuda":
        assert cfg.run.cpu_workers == 1
    assert steps.RUN_PLANS and steps.RUN_PLANS[-1]["workers_used"] == cfg.run.cpu_workers
    assert cfg.run.cpu_workers <= 64
    per = cfg.run.n_batches * sum(g.n_voxels for g in cfg.scoring) * 8
    assert cfg.run.cpu_workers * per <= cfg.run.memory_budget_bytes


def test_largest_step_grids_fit_the_budget_with_auto_workers(
    steps: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The heaviest accumulators (T10's 1 mm path grid of the oblique beam, T8/T14 slab grids of
    0.2 mm bins) are validated with 64 requested workers: the plan reduces the workers."""
    seen: list[SimulationConfig] = []

    class FakeSimulation:
        def __init__(self, cfg: SimulationConfig) -> None:
            validate(cfg)
            seen.append(cfg)

        def run(self) -> Any:
            raise _Stop

    monkeypatch.setattr(steps, "Simulation", FakeSimulation)
    from ionmc.geometry import BoxPhantom
    from ionmc.materials import WATER

    r = steps.r_csda_mm(150.0)
    u = np.ones(3) / np.sqrt(3.0)
    grid = steps._t10_grid(-0.5 * r * u, u, r)
    geo = BoxPhantom((-100.0, -100.0, -100.0), (200.0, 200.0, 200.0), WATER)
    with pytest.raises(_Stop):
        steps.run_cfg(energy=150.0, geometry=geo, scoring=(grid,), n=2000, n_batches=10,
                      workers=64)  # fmt: skip
    plan = steps.RUN_PLANS[-1]
    assert plan["accumulator_bytes_per_worker"] == 10 * grid.n_voxels * 8
    assert plan["workers_used"] == min(
        64, plan["memory_budget_bytes"] // plan["accumulator_bytes_per_worker"]
    )
    assert seen[0].run.cpu_workers == plan["workers_used"]


def test_resolve_workers_auto_uses_every_cpu(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(os, "cpu_count", lambda: 64)
    assert _load("run_suite").resolve_workers("auto") == 64
    assert _load("run_suite").resolve_workers("12") == 12
