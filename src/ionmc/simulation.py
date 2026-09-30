"""Public entry point: ``run(config)`` dispatches to a backend after the capability check."""

from __future__ import annotations

from dataclasses import asdict

from ionmc.config import (
    SimulationConfig,
    UnsupportedConfigurationError,
    check_supported,
)
from ionmc.data import DataCache
from ionmc.physics.tables import default_low_energy_source
from ionmc.results import SimulationResult, build_result
from ionmc.transport.tables import TableSet, build_table_set


def build_tables(
    config: SimulationConfig, cache: DataCache | None = None, *, offline: bool = False
) -> TableSet:
    species = [config.source.species]
    materials = list(config.geometry.materials)
    return build_table_set(
        species, materials, default_low_energy_source(cache, offline=offline)
    )


def run(
    config: SimulationConfig,
    *,
    tables: TableSet | None = None,
    cache: DataCache | None = None,
    offline: bool = False,
) -> SimulationResult:
    check_supported(config)
    cfg = config.with_defaults()
    if tables is None:
        tables = build_tables(cfg, cache, offline=offline)
    if cfg.backend == "python":
        from ionmc.transport.reference import run_reference

        out = run_reference(cfg, tables)
    elif cfg.backend in ("warp-cpu", "warp-cuda"):
        from ionmc.transport.warp_backend import run_warp

        out = run_warp(cfg, tables)
    else:  # pragma: no cover - guarded by check_supported
        raise UnsupportedConfigurationError(cfg.backend)
    assert cfg.scoring is not None
    mass = cfg.scoring.voxel_mass_g(cfg.geometry)
    effective = cfg.requested()
    effective["physics"] = asdict(cfg.physics)
    effective["histories_per_batch"] = out["histories_per_batch"]
    effective["rng"] = out["rng"]
    effective["tables"] = tables.provenance
    return build_result(
        config.requested(),
        effective,
        out["energy"],
        mass,
        out["accounting"].per_primary(),
        {
            "timing": {
                "wall_seconds": out["wall_seconds"],
                "backend": cfg.backend,
                "kernel_seconds": out.get("kernel_seconds"),
                "host_seconds": out.get("host_seconds"),
                "upload_seconds": out.get("upload_seconds"),
                "first_launch_seconds": out.get("first_launch_seconds"),
                "device": out.get("device"),
            }
        },
    )
