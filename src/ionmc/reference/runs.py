"""Load materialized reference runs and extract depth-dose curves per engine.

A run directory is ``.ionmc-cache/reference-runs/<run-id>`` as written by
``materialize_reference_artifacts``. All readers fail closed (``RunError``) on missing files,
missing metadata or inconsistent dimensions.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .parsers import ParseError, read_metaimage, read_topas_csv

_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0}


class RunError(ValueError):
    """Raised when a reference run is incomplete or inconsistent."""


@dataclass
class ReferenceRun:
    run_id: str
    engine: str
    source_sha: str
    run_dir: Path
    files: dict[str, dict[str, Any]]  # manifest path -> {bytes, sha256}
    case: dict[str, Any]
    request: dict[str, Any]

    @property
    def histories(self) -> int:
        return int(self.case["histories"])


@dataclass
class DepthDose:
    depth_mm: np.ndarray
    dose: np.ndarray
    unit: str
    histories: int
    files: list[str]  # run-relative paths of the files used


def _json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise RunError(f"cannot read {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise RunError(f"{path}: expected a JSON object")
    return data


def load_run(run_dir: str | Path) -> ReferenceRun:
    run_dir = Path(run_dir)
    manifest = _json(run_dir / "transfer-manifest.json")
    for key in ("run_id", "engine", "source_sha", "files"):
        if key not in manifest:
            raise RunError(f"{run_dir}: manifest lacks {key}")
    case_path = run_dir / "inputs" / "case.json"
    if not case_path.is_file():
        raise RunError(f"{run_dir}: inputs/case.json missing")
    request_path = run_dir / "request.json"
    request = _json(request_path) if request_path.is_file() else {}
    case = _json(case_path)
    if "histories" not in case:
        raise RunError(f"{case_path}: no histories")
    return ReferenceRun(
        run_id=str(manifest["run_id"]),
        engine=str(manifest["engine"]),
        source_sha=str(manifest["source_sha"]),
        run_dir=run_dir,
        files=dict(manifest["files"]),
        case=case,
        request=request,
    )


def _require(run: ReferenceRun, rel: str) -> Path:
    if rel not in run.files:
        raise RunError(f"{run.run_id}: {rel} not in transfer manifest")
    path = run.run_dir / rel
    if not path.is_file():
        raise RunError(f"{run.run_id}: {rel} not materialized")
    return path


def file_hashes(run: ReferenceRun, rels: list[str]) -> dict[str, str]:
    return {rel: str(run.files[rel]["sha256"]) for rel in rels}


def _topas(run: ReferenceRun) -> DepthDose:
    rel = "work/idd_dose.csv"
    s = read_topas_csv(_require(run, rel))
    if s.bins[0] != 1 or s.bins[1] != 1:
        raise RunError(f"{rel}: expected a 1x1xN column, got bins {s.bins}")
    if "Sum" not in s.values:
        raise RunError(f"{rel}: no Sum column")
    unit = s.bin_unit[2]
    if unit not in _TO_MM:
        raise RunError(f"{rel}: unknown length unit {unit!r}")
    width = s.bin_width[2] * _TO_MM[unit]
    dose = s.values["Sum"][0, 0, :]
    if not np.all(np.isfinite(dose)):
        raise RunError(f"{rel}: missing bins")
    depth = (np.arange(dose.size) + 0.5) * width
    return DepthDose(depth, dose, f"{s.unit} (run total)", run.histories, [rel])


def _fred(run: ReferenceRun) -> DepthDose:
    rel = "work/out/score/Phantom.Dose.mhd"
    img = read_metaimage(_require(run, rel))
    nx, ny, nz = img.dims
    if img.data.shape != (nz, ny, nx):
        raise RunError(f"{rel}: inconsistent dims")
    dose = img.data.astype(np.float64).sum(axis=(1, 2))
    # FRED writes mm; Offset is the first voxel centre (0.5 mm for 1 mm voxels at z=0).
    depth = img.offset[2] + np.arange(nz) * img.spacing[2]
    if abs(float(depth[0]) - 0.5 * img.spacing[2]) > 1e-6:
        raise RunError(f"{rel}: z offset {img.offset[2]} is not half a voxel; depth origin unclear")
    return DepthDose(depth, dose, "Gy per primary (lateral sum)", run.histories, [rel])


def _mcsquare(run: ReferenceRun) -> DepthDose:
    rel = "work/Outputs/Dose.mhd"
    _require(run, "work/Outputs/Dose.raw")
    img = read_metaimage(_require(run, rel))
    nx, ny, nz = img.dims
    if img.data.shape != (nz, ny, nx):
        raise RunError(f"{rel}: inconsistent dims")
    dose_y = img.data.astype(np.float64).sum(axis=(0, 2))  # per y index
    # gantry 0: beam enters at y_max and travels in -y, so depth = y_max - y
    dose = dose_y[::-1]
    depth = (np.arange(ny) + 0.5) * img.spacing[1]
    return DepthDose(
        depth,
        dose,
        "MCsquare native per proton (lateral sum)",
        run.histories,
        [rel, "work/Outputs/Dose.raw"],
    )


_EXTRACTORS = {"topas": _topas, "fred": _fred, "mcsquare": _mcsquare}


def depth_dose(run: ReferenceRun) -> DepthDose:
    """Integral depth-dose curve (depth in mm increasing from the entrance face)."""
    try:
        extract = _EXTRACTORS[run.engine]
    except KeyError as exc:
        raise RunError(f"no depth-dose extractor for engine {run.engine!r}") from exc
    try:
        dd = extract(run)
    except ParseError as exc:
        raise RunError(str(exc)) from exc
    if dd.depth_mm.shape != dd.dose.shape or dd.dose.size < 3:
        raise RunError(f"{run.run_id}: inconsistent depth-dose arrays")
    return dd
