"""Shared fixtures for the transport-engine tests (offline: analytic Bethe stopping only)."""

from __future__ import annotations

import os
from collections.abc import Callable
from typing import Any

import pytest

from ionmc.config import (
    DiagnosticsOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
)
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.materials import WATER
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping import BetheStoppingSource
from ionmc.scoring import ScoringGrid
from ionmc.sources import PencilBeamSource

SINGLE_PROCESS = os.environ.get("IONMC_SINGLE_PROCESS") == "1"


@pytest.fixture(scope="session")
def bethe() -> BetheStoppingSource:
    """Analytic Bethe source (I = 78 eV water), never any downloaded data."""
    return BetheStoppingSource()


@pytest.fixture
def make_config(bethe: BetheStoppingSource) -> Callable[..., SimulationConfig]:
    """Factory for a small configuration (python float64 by default); keyword overrides."""

    def factory(
        *,
        energy: float = 20.0,
        n: int = 2,
        n_batches: int = 2,
        seed: int = 1,
        position: tuple[float, float, float] = (0.0, 0.0, 0.0),
        direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
        energy_sigma: float = 0.0,
        lateral_sigma: float = 0.0,
        geometry: VoxelGeometry | BoxPhantom | None = None,
        scoring: tuple[ScoringGrid, ...] | None = None,
        mcs: bool = True,
        straggling: bool = True,
        max_step: float = 2.0,
        e_cut: float = 2.0,
        max_steps: int | None = None,
        diagnostics: DiagnosticsOptions | None = None,
        backend: str = "python",
        precision: str = "float64",
        allow_invalid: bool = False,
        physics_kwargs: dict[str, Any] | None = None,
        run_kwargs: dict[str, Any] | None = None,
    ) -> SimulationConfig:
        geo = geometry or BoxPhantom((-30.0, -30.0, 0.0), (60.0, 60.0, 60.0), WATER)
        grids = scoring or (ScoringGrid((-30.0, -30.0, 0.0), (2.0, 2.0, 2.0), (30, 30, 30)),)
        return SimulationConfig(
            source=PencilBeamSource(
                PROTON, position, direction, energy, energy_sigma, lateral_sigma
            ),
            geometry=geo,
            scoring=grids,
            physics=PhysicsOptions(
                nuclear=False,
                stopping=bethe,
                straggling=straggling,
                multiple_scattering=mcs,
                e_cut_mev=e_cut,
                max_step_mm=max_step,
                **(physics_kwargs or {}),
            ),
            run=RunOptions(
                backend=backend,  # type: ignore[arg-type]
                precision=precision,  # type: ignore[arg-type]
                seed=seed,
                n_histories=n,
                n_batches=n_batches,
                max_steps=max_steps,
                allow_invalid_result=allow_invalid,
                **(run_kwargs or {}),
            ),
            diagnostics=diagnostics or DiagnosticsOptions(),
        )

    return factory


@pytest.fixture(scope="session", autouse=True)
def _default_device_is_cpu() -> None:
    """Make the suite independent of the host: Warp's default device is the CPU, so a test that
    forgets an explicit ``device=`` never allocates on a GPU while its kernel runs on the CPU
    (a segfault on a CUDA host). CUDA tests name ``cuda:0`` explicitly."""
    import warp as wp

    wp.set_device("cpu")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Tests marked ``cuda`` skip without a usable CUDA device; ``IONMC_REQUIRE_CUDA=1`` turns
    that skip into a failure (the GPU host runner sets it, so a missing GPU cannot pass).
    In the single-process diagnostic mode (``IONMC_SINGLE_PROCESS=1``) tests marked
    ``multiprocess`` (they need several worker processes) are skipped."""
    if SINGLE_PROCESS:
        skip_mp = pytest.mark.skip(reason="deferred: single-process diagnostic mode")
        for item in items:
            if item.get_closest_marker("multiprocess"):
                item.add_marker(skip_mp)
    from ionmc.config import cuda_available

    if any("cuda" in item.keywords for item in items) and not cuda_available():
        required = os.environ.get("IONMC_REQUIRE_CUDA") == "1"
        for item in items:
            if "cuda" in item.keywords:
                if required:
                    item.fixturenames.insert(0, "_cuda_required_fail")
                else:
                    item.add_marker(pytest.mark.skip(reason="no CUDA device"))


@pytest.fixture
def _cuda_required_fail() -> None:
    pytest.fail("IONMC_REQUIRE_CUDA=1 but no CUDA device is available")
