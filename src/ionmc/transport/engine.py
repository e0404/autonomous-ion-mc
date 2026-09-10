"""Transport engine: reference Python and Warp execution of CSDA depth dose.

The engine runs a monoenergetic proton pencil beam through a homogeneous water
slab and returns the integral depth dose (energy deposited per depth bin). The
**reference Python** path is a scalar per-history loop calling the float64
shared-source step functions (:mod:`ionmc.physics.transport`); the **warp** path
launches the equivalent kernel on the CPU or CUDA device (decision ``0009``).
Both paths use the same step algorithm and the same tabulated stopping power,
so the reference path is a faithful oracle for the accelerated one.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from ionmc.backend import mathlib, reference
from ionmc.data.stopping_tables import StoppingTable
from ionmc.transport.depth_dose import DepthDoseGrid
from ionmc.transport.geometry import WaterSlab
from ionmc.transport.source import PencilBeamSource
from ionmc.transport.state import ParticleState, Status

#: Default fractional energy loss per step and the absolute step cap [mm].
DEFAULT_MAX_FRACTION: float = 0.02
DEFAULT_MAX_STEP_MM: float = 1.0
#: Default energy cutoff [MeV]; below it the residual energy is deposited locally.
DEFAULT_ENERGY_CUT_MEV: float = 0.5
#: Hard per-history step-count cap (a truncated history is counted, not dropped).
DEFAULT_MAX_STEPS: int = 100_000


@dataclass(frozen=True)
class DepthDoseResult:
    """Outcome of a depth-dose transport run."""

    edep_mev: np.ndarray  # energy deposited per depth bin [MeV]
    grid: DepthDoseGrid
    n_histories: int
    energy_in_mev: float  # total kinetic energy started
    truncated: int  # histories that hit the step cap
    path: str
    device: str | None

    @property
    def energy_deposited_mev(self) -> float:
        return float(np.sum(self.edep_mev))

    @property
    def energy_balance(self) -> float:
        """Relative energy imbalance ``(deposited - in) / in``."""
        return (self.energy_deposited_mev - self.energy_in_mev) / self.energy_in_mev

    def r80_mm(self) -> float:
        """Distal 80 %-of-maximum depth of the integral depth dose [mm]."""
        return _distal_level_depth(self.edep_mev, self.grid, 0.80)


class TransportEngine:
    """Runs CSDA depth-dose transport for one stopping-power table and geometry."""

    def __init__(
        self,
        table: StoppingTable,
        slab: WaterSlab,
        grid: DepthDoseGrid,
        max_fraction: float = DEFAULT_MAX_FRACTION,
        max_step_mm: float = DEFAULT_MAX_STEP_MM,
        energy_cut_mev: float = DEFAULT_ENERGY_CUT_MEV,
        max_steps: int = DEFAULT_MAX_STEPS,
    ) -> None:
        if not (0.0 < max_fraction < 1.0):
            raise ValueError("max_fraction must be in (0, 1)")
        self.table = table
        self.slab = slab
        self.grid = grid
        self.max_fraction = max_fraction
        self.max_step_mm = max_step_mm
        self.energy_cut_mev = max(energy_cut_mev, float(table.energy_mev[0]))
        self.max_steps = max_steps

    def run(
        self,
        source: PencilBeamSource,
        n_histories: int = 1,
        seed: int = 12345,
        path: str = "python",
        device: str = "cpu",
    ) -> DepthDoseResult:
        state = source.sample(n_histories, seed)
        energy_in = float(np.sum(state.energy_mev * state.weight))
        if path == "warp":
            edep, truncated = self._run_warp(state, device)
        elif path == "python":
            edep, truncated = self._run_reference(state)
        else:
            raise ValueError(
                f"unknown transport path {path!r} (use 'python' or 'warp')"
            )
        return DepthDoseResult(
            edep_mev=edep,
            grid=self.grid,
            n_histories=n_histories,
            energy_in_mev=energy_in,
            truncated=truncated,
            path=path,
            device=device if path == "warp" else None,
        )

    # -- reference Python path ------------------------------------------------

    def _run_reference(self, state: ParticleState) -> tuple[np.ndarray, int]:
        tp = reference.load_bound_module(
            "ionmc.physics.transport",
            "python",
            rebind_dependencies=["ionmc.physics.tabulated"],
        )
        t = self.table
        args: tuple[Any, ...] = (
            t.energy_mev,
            t.stopping_mev_cm2_per_g,
            t.slope,
            t.size,
            t.bisection_steps,
        )
        density = self.slab.density_g_per_cm3
        dz = self.grid.bin_width_mm
        geom_depth = self.slab.depth_mm  # transport is bounded by the medium...
        n_bins = self.grid.n_bins  # ...scoring only within the grid extent
        edep = self.grid.empty()
        truncated = 0
        for h in range(state.size):
            e = float(state.energy_mev[h])
            z = float(state.position_mm[h, 2])
            w = float(state.weight[h])
            step = 0
            status = Status.ALIVE
            while status == Status.ALIVE and step < self.max_steps:
                cur_bin = math.floor(z / dz)
                boundary = (cur_bin + 1) * dz
                dl_e = tp.energy_loss_step_length(
                    e, self.max_fraction, self.max_step_mm, density, *args
                )
                dl = min(dl_e, boundary - z + 1.0e-6, geom_depth - z)
                de = tp.midpoint_energy_loss(e, dl, density, *args)
                if 0 <= cur_bin < n_bins:
                    edep[cur_bin] += w * de
                z += dl
                e -= de
                step += 1
                if e <= self.energy_cut_mev:
                    dep_bin = math.floor(z / dz)
                    if 0 <= dep_bin < n_bins:
                        edep[dep_bin] += w * e
                    e = 0.0
                    status = Status.STOPPED
                if z >= geom_depth:
                    status = Status.ESCAPED
            if step >= self.max_steps and status == Status.ALIVE:
                truncated += 1
        return edep, truncated

    # -- Warp path ------------------------------------------------------------

    def _run_warp(self, state: ParticleState, device: str) -> tuple[np.ndarray, int]:
        if not mathlib.HAVE_WARP:
            raise ImportError(
                "warp-lang is not installed; the warp path is unavailable"
            )
        from ionmc.backend.warp_transport import DepthDoseKernel

        kernel = DepthDoseKernel(self.table, device)
        return kernel.run(
            energy0=state.energy_mev,
            z0=state.position_mm[:, 2],
            weight=state.weight,
            density=self.slab.density_g_per_cm3,
            max_fraction=self.max_fraction,
            max_step_mm=self.max_step_mm,
            bin_width_mm=self.grid.bin_width_mm,
            geom_depth_mm=self.slab.depth_mm,
            energy_cut_mev=self.energy_cut_mev,
            max_steps=self.max_steps,
            n_bins=self.grid.n_bins,
        )


def _distal_level_depth(edep: np.ndarray, grid: DepthDoseGrid, level: float) -> float:
    """Depth [mm] where the profile falls through ``level * max`` on the distal side.

    Linear interpolation between the two bin centres bracketing the crossing
    beyond the maximum. Returns ``nan`` if the profile is all zero.
    """
    if not np.any(edep > 0.0):
        return float("nan")
    peak = int(np.argmax(edep))
    threshold = level * edep[peak]
    centers = grid.centers_mm
    for i in range(peak, len(edep) - 1):
        if edep[i] >= threshold >= edep[i + 1]:
            e0, e1 = edep[i], edep[i + 1]
            c0, c1 = centers[i], centers[i + 1]
            if e0 == e1:
                return float(c0)
            frac = (e0 - threshold) / (e0 - e1)
            return float(c0 + frac * (c1 - c0))
    return float(centers[-1])
