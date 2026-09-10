"""Proton transport in homogeneous water (Stage 1).

Implemented so far (task DEV-004, decision ``0009``): continuous-slowing-down
longitudinal transport of a monoenergetic proton pencil beam with integral
depth-dose scoring, on the reference Python path and the Warp CPU/CUDA path.
Energy-loss straggling and multiple Coulomb scattering are later Stage-1 tasks.
"""

from __future__ import annotations

from ionmc.transport.depth_dose import DepthDoseGrid
from ionmc.transport.engine import DepthDoseResult, TransportEngine
from ionmc.transport.geometry import WaterSlab
from ionmc.transport.source import PencilBeamSource
from ionmc.transport.state import ParticleState, Species

__all__ = [
    "DepthDoseGrid",
    "DepthDoseResult",
    "ParticleState",
    "PencilBeamSource",
    "Species",
    "TransportEngine",
    "WaterSlab",
]
