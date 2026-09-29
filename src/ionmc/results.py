"""Simulation results: in-memory container and persisted, self-describing files.

``SimulationResult.save(path)`` writes ``<path>.npz`` (arrays) and
``<path>.json`` (metadata) so a separate process can reopen the result with
``load``. The metadata includes the requested and effective configuration,
code identity, geometry/scoring descriptions with coordinates, units, table
provenance, energy accounting, histories, batches, seeds/RNG identity and
the uncertainty definition. Fields that are undefined are explicit ``null``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.provenance import code_identity
from ionmc.scoring import GY_PER_MEV_PER_G


@dataclass
class SimulationResult:
    energy_mean: np.ndarray  # MeV per primary per scoring voxel
    energy_stderr: np.ndarray  # standard error of the batch mean (NaN if undefined)
    dose_mean: np.ndarray  # Gy per primary
    dose_stderr: np.ndarray
    mass_g: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def scoring_edges_mm(self) -> list[np.ndarray]:
        s = self.metadata["scoring"]
        return [
            np.asarray(s["origin_mm"][a])
            + s["spacing_mm"][a] * np.arange(s["shape"][a] + 1)
            for a in range(3)
        ]

    def depth_dose(self, axis: int = 2) -> tuple[np.ndarray, np.ndarray]:
        """Integral (laterally summed) energy per primary along ``axis`` and its bin centres."""
        other = tuple(a for a in range(3) if a != axis)
        prof = self.energy_mean.sum(axis=other)
        s = self.metadata["scoring"]
        centers = s["origin_mm"][axis] + s["spacing_mm"][axis] * (
            np.arange(s["shape"][axis]) + 0.5
        )
        return centers, prof

    def save(self, path: str | Path) -> tuple[Path, Path]:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        npz = path.with_suffix(".npz")
        js = path.with_suffix(".json")
        np.savez_compressed(
            npz,
            energy_mean=self.energy_mean,
            energy_stderr=self.energy_stderr,
            dose_mean=self.dose_mean,
            dose_stderr=self.dose_stderr,
            mass_g=self.mass_g,
        )
        meta = dict(self.metadata)
        meta["arrays"] = {
            "file": npz.name,
            "energy_mean": {
                "units": "MeV per primary",
                "shape": list(self.energy_mean.shape),
            },
            "energy_stderr": {
                "units": "MeV per primary",
                "definition": "standard error of the mean over independent batches; NaN where undefined",
            },
            "dose_mean": {"units": "Gy per primary"},
            "dose_stderr": {"units": "Gy per primary"},
            "mass_g": {"units": "g per scoring voxel"},
        }
        js.write_text(
            json.dumps(meta, indent=2, sort_keys=True, default=_json_default) + "\n"
        )
        return npz, js

    @classmethod
    def load(cls, path: str | Path) -> SimulationResult:
        path = Path(path)
        meta = json.loads(path.with_suffix(".json").read_text())
        with np.load(path.with_suffix(".npz")) as data:
            return cls(
                data["energy_mean"],
                data["energy_stderr"],
                data["dose_mean"],
                data["dose_stderr"],
                data["mass_g"],
                meta,
            )


def _json_default(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(f"not JSON serializable: {type(obj)}")


def build_result(
    config_requested: dict,
    config_effective: dict,
    tally,
    mass_g: np.ndarray,
    accounting: dict,
    extra: dict,
) -> SimulationResult:
    energy_mean = tally.mean()
    energy_err = tally.standard_error()
    with np.errstate(divide="ignore", invalid="ignore"):
        dose = np.where(mass_g > 0, energy_mean / mass_g * GY_PER_MEV_PER_G, np.nan)
        dose_err = np.where(mass_g > 0, energy_err / mass_g * GY_PER_MEV_PER_G, np.nan)
    meta = {
        "format": "ionmc-result/1",
        "code": code_identity(),
        "requested": config_requested,
        "effective": config_effective,
        "scoring": config_effective["scoring"],
        "geometry": config_effective["geometry"],
        "energy_accounting_per_primary": accounting,
        "uncertainty": {
            "method": "batch statistics",
            "batches": config_effective["batches"],
            "definition": "standard error of the mean of per-batch per-primary values; NaN where fewer than 2 batches",
        },
        **extra,
    }
    return SimulationResult(energy_mean, energy_err, dose, dose_err, mass_g, meta)


def write_text_summary(result: SimulationResult, path: str | Path) -> Path:
    """Human-readable summary: run identity, energy accounting and the integral depth-dose."""
    path = Path(path).with_suffix(".txt")
    m = result.metadata
    z, prof = result.depth_dose(2)
    acc = m["energy_accounting_per_primary"]
    lines = [
        "IonMC simulation summary",
        f"code: {m['code']}",
        f"backend: {m['effective']['backend']}  precision: {m['effective']['precision']}",
        f"source: {m['effective']['source']}",
        f"geometry: {m['geometry']['name']} shape {m['geometry']['shape']} spacing {m['geometry']['spacing_mm']} mm",
        f"scoring: shape {m['scoring']['shape']} spacing {m['scoring']['spacing_mm']} mm origin {m['scoring']['origin_mm']} mm",
        f"histories: {m['effective']['histories']} in {m['effective']['batches']} batches, seed {m['effective']['seed']}",
        "energy per primary (MeV): "
        + ", ".join(f"{k}={v:.4f}" for k, v in acc.items() if k != "steps_per_history"),
        f"steps per history: {acc['steps_per_history']:.1f}",
        f"max dose: {np.nanmax(result.dose_mean):.4e} Gy/primary; relative standard error at max: "
        + (
            f"{float(result.dose_stderr.flat[int(np.nanargmax(result.dose_mean))] / np.nanmax(result.dose_mean)):.3f}"
            if np.isfinite(result.dose_stderr).any()
            else "undefined"
        ),
        "",
        "integral depth-energy profile (laterally summed), z centre [mm]  energy [MeV/primary]  rel. std. error",
    ]
    err = np.sqrt(np.nansum(np.nan_to_num(result.energy_stderr) ** 2, axis=(0, 1)))
    for zi, ei, si in zip(z, prof, err, strict=True):
        rel = si / ei if ei > 0 else float("nan")
        lines.append(f"{zi:10.2f}  {ei:14.6e}  {rel:8.3f}")
    path.write_text("\n".join(lines) + "\n")
    return path
