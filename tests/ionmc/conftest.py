"""Shared fixtures for the transport-engine tests (offline: analytic Bethe stopping only)."""

from __future__ import annotations

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


@pytest.fixture(scope="session")
def bethe() -> BetheStoppingSource:
    """Analytic Bethe source (I = 78 eV water), never any downloaded data."""
    return BetheStoppingSource()


@pytest.fixture
def make_config(bethe: BetheStoppingSource) -> Callable[..., SimulationConfig]:
    """Factory for a small python-backend configuration; keyword overrides."""

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
                backend="python",
                precision="float64",
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
