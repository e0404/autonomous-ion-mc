"""Proton transport in homogeneous water (Stage 1).

Implemented so far: continuous-slowing-down longitudinal transport of a
monoenergetic proton pencil beam with integral depth-dose scoring, on the
reference Python path and the Warp CPU/CUDA path (task DEV-004, decision
``0009``); Bohr energy-loss straggling giving a realistic Bragg peak, with
batch-based statistical uncertainty (task DEV-005, decision ``0010``).
Multiple Coulomb scattering and lateral spread are a later Stage-1 task.
"""

from __future__ import annotations

from ionmc.transport.depth_dose import DepthDoseGrid, DepthLateralGrid, DoseGrid3D
from ionmc.transport.engine import (
    BatchedDepthDoseResult,
    BatchedDoseResult,
    DepthDoseResult,
    ScatteringResult,
    TransportEngine,
)
from ionmc.transport.geometry import VoxelGrid3D, WaterSlab
from ionmc.transport.source import PencilBeamSource
from ionmc.transport.state import ParticleState, Species

__all__ = [
    "BatchedDepthDoseResult",
    "BatchedDoseResult",
    "DepthDoseGrid",
    "DepthDoseResult",
    "DepthLateralGrid",
    "DoseGrid3D",
    "ParticleState",
    "PencilBeamSource",
    "ScatteringResult",
    "Species",
    "TransportEngine",
    "VoxelGrid3D",
    "WaterSlab",
]
