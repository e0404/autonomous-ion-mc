"""Simulation configuration and the fail-closed capability contract (V2-CAP).

``PhysicsConfig`` lists every physics switch the public API accepts. A
backend must either honour a requested switch or raise
:class:`UnsupportedConfigurationError` *before* producing any result; a
requested switch is never silently ignored. ``capabilities()`` reports what
each backend supports so users can discover the contract without reading
code.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace

from ionmc.geometry import VoxelGeometry
from ionmc.scoring import ScoringGrid
from ionmc.sources import PencilBeam


class UnsupportedConfigurationError(ValueError):
    """Requested physics/scoring/backend combination is not implemented."""


@dataclass(frozen=True)
class PhysicsConfig:
    energy_loss: bool = True
    straggling: bool = True
    multiple_scattering: bool = True
    nuclear: bool = False  # not implemented yet: requesting it raises
    transport_secondaries: bool = False  # not implemented yet
    max_step_mm: float = 1.0
    energy_step_fraction: float = 0.10
    cutoff_mev_per_u: float = 0.5
    straggling_model: str = "bohr-gamma"
    mcs_model: str = "differential-moliere"
    lateral_displacement: str = "random-hinge"

    def validate(self) -> None:
        if not 0.0 < self.max_step_mm <= 100.0:
            raise ValueError("max_step_mm must be in (0, 100]")
        if not 0.0 < self.energy_step_fraction <= 0.5:
            raise ValueError("energy_step_fraction must be in (0, 0.5]")
        if self.cutoff_mev_per_u <= 0.0:
            raise ValueError("cutoff must be positive")
        if self.straggling_model not in ("bohr-gamma",):
            raise UnsupportedConfigurationError(
                f"straggling model {self.straggling_model!r}"
            )
        if self.mcs_model not in ("differential-moliere",):
            raise UnsupportedConfigurationError(f"mcs model {self.mcs_model!r}")
        if self.lateral_displacement not in ("random-hinge",):
            raise UnsupportedConfigurationError(
                f"lateral displacement {self.lateral_displacement!r}"
            )


@dataclass(frozen=True)
class SimulationConfig:
    source: PencilBeam
    geometry: VoxelGeometry
    histories: int
    batches: int = 10
    seed: int = 12345
    backend: str = "python"
    precision: str = "float64"
    scoring: ScoringGrid | None = None
    physics: PhysicsConfig = field(default_factory=PhysicsConfig)
    scorers: tuple[str, ...] = ("energy", "dose")

    def validate(self) -> None:
        if self.histories < 1 or self.batches < 1 or self.batches > self.histories:
            raise ValueError("need histories >= batches >= 1")
        self.physics.validate()
        for s in self.scorers:
            if s not in ("energy", "dose"):
                raise UnsupportedConfigurationError(f"scorer {s!r} is not implemented")

    def with_defaults(self) -> SimulationConfig:
        return replace(
            self, scoring=self.scoring or ScoringGrid.from_geometry(self.geometry)
        )

    def requested(self) -> dict:
        return {
            "source": self.source.describe(),
            "geometry": self.geometry.describe(),
            "scoring": (
                self.scoring or ScoringGrid.from_geometry(self.geometry)
            ).describe(),
            "histories": self.histories,
            "batches": self.batches,
            "seed": self.seed,
            "backend": self.backend,
            "precision": self.precision,
            "physics": asdict(self.physics),
            "scorers": list(self.scorers),
        }


_WARP_CAPS = {
    "precision": ["float32", "float64"],
    "physics": {
        "energy_loss": [True],
        "straggling": [True, False],
        "multiple_scattering": [True, False],
        "nuclear": [False],
        "transport_secondaries": [False],
    },
    "scorers": ["energy", "dose"],
    "sources": ["pencil-beam"],
    "geometries": ["voxel"],
    "species": "any species in ionmc.species (electromagnetic transport only)",
    "energy_range_mev_per_u": [1.0, 1000.0],
}

BACKEND_CAPABILITIES: dict[str, dict] = {
    "warp-cpu": dict(_WARP_CAPS),
    "warp-cuda": dict(_WARP_CAPS),
    "python": {
        "precision": ["float64"],
        "physics": {
            "energy_loss": [True],
            "straggling": [True, False],
            "multiple_scattering": [True, False],
            "nuclear": [False],
            "transport_secondaries": [False],
        },
        "scorers": ["energy", "dose"],
        "sources": ["pencil-beam"],
        "geometries": ["voxel"],
        "species": "any species in ionmc.species (electromagnetic transport only)",
        "energy_range_mev_per_u": [1.0, 1000.0],
    },
}


def capabilities() -> dict:
    """Discoverable capability contract of the installed backends."""
    return BACKEND_CAPABILITIES


def check_supported(config: SimulationConfig) -> None:
    """Raise UnsupportedConfigurationError for any unsupported request (fail closed)."""
    config.validate()
    caps = BACKEND_CAPABILITIES.get(config.backend)
    if caps is None:
        raise UnsupportedConfigurationError(
            f"backend {config.backend!r} is not available; known: {sorted(BACKEND_CAPABILITIES)}"
        )
    if config.precision not in caps["precision"]:
        raise UnsupportedConfigurationError(
            f"backend {config.backend!r} does not support precision {config.precision!r}"
        )
    for key, allowed in caps["physics"].items():
        value = getattr(config.physics, key)
        if value not in allowed:
            raise UnsupportedConfigurationError(
                f"backend {config.backend!r} does not support physics.{key}={value!r} (allowed: {allowed})"
            )
    for s in config.scorers:
        if s not in caps["scorers"]:
            raise UnsupportedConfigurationError(
                f"backend {config.backend!r} does not provide scorer {s!r}"
            )
    lo, hi = caps["energy_range_mev_per_u"]
    if not lo <= config.source.energy_mev_per_u <= hi:
        raise UnsupportedConfigurationError(
            f"source energy {config.source.energy_mev_per_u} MeV/u outside supported range [{lo}, {hi}]"
        )
