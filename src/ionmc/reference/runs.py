"""Load materialized reference runs and extract depth-dose curves per engine.

A run directory is ``.ionmc-cache/reference-runs/<run-id>`` as written by
``materialize_reference_artifacts``. All readers fail closed (``RunError``) on missing files,
missing metadata or inconsistent dimensions.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .parsers import ParseError, parse_metaimage, parse_topas_csv

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
    read_digests: dict[str, str] = field(default_factory=dict)  # sha256 of the bytes read

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


def _json_bytes(blob: bytes, label: str) -> dict[str, Any]:
    try:
        data = json.loads(blob.decode("utf-8"))
    except ValueError as exc:
        raise RunError(f"cannot parse {label}: {exc}") from exc
    if not isinstance(data, dict):
        raise RunError(f"{label}: expected a JSON object")
    return data


def _verified_bytes(
    run_dir: Path, files: dict[str, dict[str, Any]], rel: str, digests: dict[str, str]
) -> bytes:
    """Read ``rel`` and check size and sha256 against the transfer manifest (fail closed)."""
    entry = files.get(rel)
    if not isinstance(entry, dict) or "sha256" not in entry or "bytes" not in entry:
        raise RunError(f"{run_dir.name}: {rel} not in transfer manifest")
    path = run_dir / rel
    if path.is_symlink() or not path.is_file():
        raise RunError(f"{run_dir.name}: {rel} not materialized as a regular file")
    blob = path.read_bytes()
    if len(blob) != entry["bytes"]:
        raise RunError(f"{run_dir.name}: {rel} size {len(blob)} != manifest {entry['bytes']}")
    digest = hashlib.sha256(blob).hexdigest()
    if digest != entry["sha256"]:
        raise RunError(f"{run_dir.name}: {rel} sha256 mismatch against transfer manifest")
    digests[rel] = digest
    return blob


def load_run(run_dir: str | Path) -> ReferenceRun:
    run_dir = Path(run_dir)
    try:
        manifest = _json_bytes((run_dir / "transfer-manifest.json").read_bytes(), "manifest")
    except OSError as exc:
        raise RunError(f"{run_dir}: transfer-manifest.json unreadable: {exc}") from exc
    for key in ("run_id", "engine", "source_sha", "files"):
        if key not in manifest:
            raise RunError(f"{run_dir}: manifest lacks {key}")
    if not isinstance(manifest["files"], dict):
        raise RunError(f"{run_dir}: manifest files must be an object")
    files = dict(manifest["files"])
    digests: dict[str, str] = {}
    case = _json_bytes(
        _verified_bytes(run_dir, files, "inputs/case.json", digests), "inputs/case.json"
    )
    if "histories" not in case:
        raise RunError(f"{run_dir}: case.json has no histories")
    request: dict[str, Any] = {}
    if "request.json" in files:
        request = _json_bytes(
            _verified_bytes(run_dir, files, "request.json", digests), "request.json"
        )
    return ReferenceRun(
        run_id=str(manifest["run_id"]),
        engine=str(manifest["engine"]),
        source_sha=str(manifest["source_sha"]),
        run_dir=run_dir,
        files=files,
        case=case,
        request=request,
        read_digests=digests,
    )


def _read(run: ReferenceRun, rel: str) -> bytes:
    return _verified_bytes(run.run_dir, run.files, rel, run.read_digests)


def file_hashes(run: ReferenceRun, rels: list[str]) -> dict[str, str]:
    """sha256 digests of the bytes actually read (verified against the manifest)."""
    missing = [r for r in rels if r not in run.read_digests]
    if missing:
        raise RunError(f"{run.run_id}: files not read/verified: {missing}")
    return {rel: run.read_digests[rel] for rel in rels}


def _mhd(run: ReferenceRun, rel: str) -> Any:
    """Parse a manifested .mhd; a separate raw file must be a plain name in the same directory."""
    directory = rel.rpartition("/")[0]

    def loader(name: str) -> bytes:
        if "/" in name or "\\" in name or name in ("", ".", ".."):
            raise RunError(f"{rel}: unsafe ElementDataFile {name!r}")
        return _read(run, f"{directory}/{name}" if directory else name)

    return parse_metaimage(_read(run, rel), rel, loader)


def _topas(run: ReferenceRun) -> DepthDose:
    rel = "work/idd_dose.csv"
    s = parse_topas_csv(_read(run, rel).decode("utf-8"), rel)
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
    img = _mhd(run, rel)
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
    img = _mhd(run, rel)
    if img.header.get("ElementDataFile") != "Dose.raw":
        raise RunError(f"{rel}: ElementDataFile must be Dose.raw")
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
