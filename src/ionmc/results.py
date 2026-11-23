"""Persisted result bundle ``ionmc.result/1`` (decision 0045 sections 1 and 2).

A bundle is a directory with ``result.json`` (strict canonical JSON), ``arrays.npz`` (deterministic
zip of ``.npy`` members), ``report.md`` and (later) ``plots/*.svg``. :func:`save` writes it
atomically from a :class:`~ionmc.simulation.Result`; :func:`load` reads it fail-closed
(:class:`ResultFormatError`) under explicit resource budgets and never reconstructs transport
internals or recomputes ``evidence_at_run``.

Scope of this module version (plan V3-006, step S2a): ``report.md`` is a minimal deterministic
placeholder generated here (S5 replaces it by the full report and adds ``plots/``); the
``effective_extra.resolutions`` list is empty with an explicit ``status`` (``reconcile`` is S4);
``provenance.code`` is built from the run's ``describe_environment()`` record (S3 replaces the
builder). No table values are ever persisted: array names are checked against an allowlist.

Strictness of ``load``: unknown/missing fields are refused at every level of the blocks this
module defines; ``requested`` is checked by ``configio.validate_dict_structure``,
``input_manifest`` by key forms, ``effective`` and ``provenance`` against the DECLARED schema of
``ionmc._result_shapes`` (every nested field typed; open engine tables listed there) and bound to
``requested``, the run block and the arrays (see ``_check_effective_binding``).

Integrity is two-layered: the sha256 of every file (``files``; depends on the zlib level) and the
``content_sha256`` of every array (``dtype.str | shape | C-order bytes``; the digest of
:func:`ionmc.configio.array_content_sha256`, independent of compression); the content digests are
the authoritative identity of the data.
"""

from __future__ import annotations

import contextlib
import getpass
import hashlib
import json
import math
import os
import platform
import re
import shutil
import socket
import struct
import tempfile
import zipfile
import zlib
from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field, fields, replace
from pathlib import Path
from typing import Any

import numpy as np

from ionmc import _result_shapes as shapes
from ionmc.config import validate
from ionmc.configio import (
    ConfigError,
    DatasetError,
    InputManifest,
    array_content_sha256,
    from_dict,
    manifest_for_api_config,
    reader_item_sha256s,
    sha256_canonical,
    to_dict,
    validate_dict_structure,
)
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import VoxelGeometry
from ionmc.lookup import LookupTable
from ionmc.physics.projectiles import PROJECTILES
from ionmc.physics.stopping_selection import validate_identity_block
from ionmc.scoring import MEV_PER_G_TO_GY, min_defined_batches
from ionmc.simulation import (
    IonTransportCounters,
    NuclearEnergyBalance,
    NuclearTransportCounters,
    Result,
    TransportCounters,
)
from ionmc.transport.channels import (
    CLASS_LOCAL,
    CLASS_STEP,
    GEN_MAX,
    _generation_range,
    _species_row,
    _spectrum_spec,
)
from ionmc.transport.tally import NUCLEAR_TALLY_NAMES

SCHEMA = "ionmc.result"
SCHEMA_VERSION = "1.0"
KNOWN_VERSIONS = ("1.0",)
KIND = "transport"

MAX_TOTAL_BYTES = 8 * 2**30
MAX_ARRAY_BYTES = 4 * 2**30
MAX_JSON_BYTES = 16 * 2**20
MAX_CENTRAL_DIRECTORY_BYTES = 2**20
MAX_JSON_TOKENS = 2_000_000
MAX_JSON_DEPTH = 64
MAX_NPY_HEADER_BYTES = 1024
ZLIB_LEVEL = 6
_EOCD_WINDOW = 65_557

GRID_CONVENTION = (
    "origin = corner of voxel (0,0,0); C order [ix, iy, iz]; lengths in mm (decision 0037)"
)
UNIT_VOCABULARY = frozenset(
    {"MeV", "Gy", "g", "mm", "mm^-2", "keV/um", "MeV/mm", "g/cm3", "1", "bool", "count"}
)
UNITS_TABLE = {
    "energy": "MeV",
    "dose": "Gy",
    "length": "mm",
    "mass": "g",
    "density": "g/cm3",
    "let": "keV/um",
    "fluence": "mm^-2",
    "counts": "count",
    "dimensionless": "1",
    "mask": "bool",
}

# array-name allowlist: (base, units, definition) of every persisted array
GRID_ARRAYS: dict[str, tuple[str, str]] = {
    "energy_mev": ("MeV", "mean deposited energy per primary history"),
    "energy_std_mev": ("MeV", "standard error of the mean of energy_mev (batch statistics)"),
    "dose_gy": ("Gy", "dose per primary history (NaN where the voxel mass is 0)"),
    "dose_std_gy": ("Gy", "standard error of the mean of dose_gy"),
    "relative_uncertainty": ("1", "std / mean (NaN where the mean is 0)"),
    "defined_mask": ("bool", "positive mean energy and positive voxel mass"),
    "mass_g": ("g", "exact voxel mass"),
    "n_nonzero_batches": ("count", "batches with a nonzero sum"),
    "sparse_mask": ("bool", "defined voxels with fewer than B/2 nonzero batches"),
    "batch_energy_mev": ("MeV", "per-batch energy estimates E_b / (N/B)"),
}
QUANTITY_ARRAYS = ("mean", "std", "defined_mask", "n_nonzero", "rounding_bound")
GEOMETRY_ARRAYS: dict[str, tuple[str, str]] = {
    "material_index": ("1", "index of the voxel material in the geometry material list"),
    "density_g_cm3": ("g/cm3", "per-voxel density override"),
}
DIAGNOSTIC_ARRAYS: dict[str, tuple[str, str]] = {
    "end_position_mm": ("mm", "track end position"),
    "end_code": ("1", "track end reason code"),
    "end_direction": ("1", "track end unit direction"),
    "control_residual": ("1", "control residual per track"),
    "control_displacement": ("1", "control displacement per track"),
    "end_energy_mev": ("MeV", "track end kinetic energy"),
    "escape_history": ("count", "history index of an escaping particle"),
    "escape_position_mm": ("mm", "escape position"),
    "escape_direction": ("1", "escape unit direction"),
    "escape_energy_mev": ("MeV", "escape kinetic energy"),
    "trace_end_history": ("count", "traced history index"),
    "trace_end_code": ("1", "end code of a traced history"),
    "trace_end_energy_mev": ("MeV", "end energy of a traced history"),
}
EFFECTIVE_TOP_KEYS = shapes.EFFECTIVE_TOP_KEYS
STAT_TEXT = {
    "normalisation": "per primary history",
    "uncertainty": (
        "batch standard error of the mean over B batches (history h in batch h mod B); "
        "ratio quantities by the delta method (decision 0040)"
    ),
    "defined_mask": "positive mean energy and positive voxel mass",
    "sparse_mask": "defined voxels with fewer than B/2 nonzero batches",
}
LOW_BATCH_COUNT = 10  # GridResult.low_batch_count is n_batches < 10 (simulation.py)
_ALLOWED_DTYPES = frozenset({"<f8", "<f4", "<i4", "<i8", "<u4", "|b1", "|u1", "|i1"})
_NAME_RE = re.compile(r"^[a-z0-9_-]{1,64}$")
_MEMBER_RE = re.compile(r"^[a-z0-9_-]{1,64}\.npy$")
_FILE_RE = re.compile(r"^(?:[a-z0-9_-]{1,64}\.[a-z0-9]{1,8}|plots/[a-z0-9_-]{1,64}\.svg)$")
_SHA_RE = re.compile(r"^[0-9a-f]{64}$")
_TOKEN_RE = re.compile(
    rb'"(?:[^"\\]++|\\.)*+"|-?[0-9][0-9.eE+\-]*+|true|false|null|[{}\[\]:,]', re.DOTALL
)


class ResultFormatError(Exception):
    """A result bundle is malformed, tampered with, unsafe or over budget (fail closed)."""


def _fail(msg: str) -> ResultFormatError:
    return ResultFormatError(msg)


# ---------------------------------------------------------------------------------------------
# safe names and canonical JSON


def safe_name(raw: str, used: set[str], maxlen: int = 20) -> str:
    """Map a user name to ``[a-z0-9_-]`` (lower case, <= ``maxlen``); collisions get an index
    suffix. ``used`` is updated."""
    base = re.sub(r"[^a-z0-9_-]+", "-", str(raw).lower()).strip("-")[:maxlen] or "x"
    name, i = base, 1
    while name in used:
        suffix = f"-{i}"
        name = base[: maxlen - len(suffix)] + suffix
        i += 1
    used.add(name)
    return name


def canonical_bytes(doc: Any) -> bytes:
    """Canonical result.json bytes: sorted keys, two-space indent, UTF-8, one final newline."""
    try:
        text = json.dumps(doc, allow_nan=False, sort_keys=True, indent=2, ensure_ascii=False)
        return text.encode("utf-8") + b"\n"
    except (ValueError, TypeError, UnicodeEncodeError) as exc:
        raise _fail(f"result document is not strict JSON: {exc}") from exc


def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _sha256_path(path: Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
            n += len(chunk)
    return h.hexdigest(), n


def _scalar(value: float | None, units: str, reason: str = "") -> dict[str, Any]:
    if value is None or not math.isfinite(value):
        return {"value": None, "status": "undefined", "reason": reason, "units": units}
    return {"value": float(value), "units": units}


# ---------------------------------------------------------------------------------------------
# catalogue (writer side)


@dataclass(frozen=True)
class _Arr:
    name: str
    data: np.ndarray[Any, Any]
    axes: tuple[str, ...]
    grid: str | None
    units: str
    definition: str
    origin: dict[str, Any]

    def entry(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "dtype": self.data.dtype.str,
            "shape": [int(n) for n in self.data.shape],
            "axes": list(self.axes),
            "grid": self.grid,
            "units": self.units,
            "definition": self.definition,
            "content_sha256": array_content_sha256(self.data),
            "origin": self.origin,
        }


_SPATIAL = ("x", "y", "z")


def _member_pattern(origin: Mapping[str, Any]) -> re.Pattern[str] | None:
    kind, arr = origin.get("kind"), origin.get("array")
    if not isinstance(arr, str):
        return None
    sn = r"[a-z0-9_-]{1,20}"
    if kind == "grid" and arr in GRID_ARRAYS:
        return re.compile(rf"grid__{sn}__{re.escape(arr)}")
    if kind == "quantity" and arr in QUANTITY_ARRAYS:
        return re.compile(rf"grid__{sn}__q__{sn}__{re.escape(arr)}")
    if kind == "geometry" and arr in GEOMETRY_ARRAYS:
        return re.compile(rf"geometry__{re.escape(arr)}")
    if kind == "diagnostic" and (arr in DIAGNOSTIC_ARRAYS or arr.startswith("trace_")):
        return re.compile(rf"diag__{re.escape(arr)}")
    return None


def _trace_columns() -> tuple[str, ...]:
    from ionmc.transport.tally import TRACE_COLUMNS

    return tuple(TRACE_COLUMNS)


def _unit_ok(units: str, origin: Mapping[str, Any]) -> bool:
    if units in UNIT_VOCABULARY:
        return True
    q = origin.get("quantity")
    return bool(
        origin.get("kind") == "quantity"
        and isinstance(q, str)
        and q.startswith("lookup_")
        and units
        and units.isprintable()
    )


def check_allowlist(entries: list[dict[str, Any]]) -> None:
    """Refuse any array name outside the allowlist (never persisted: table arrays)."""
    seen: set[str] = set()
    for e in entries:
        name, origin = e.get("name"), e.get("origin")
        if not isinstance(name, str) or not isinstance(origin, Mapping):
            raise _fail(f"catalogue entry without name or origin: {e!r}")
        pat = _member_pattern(origin)
        if pat is None or not pat.fullmatch(name):
            raise _fail(f"array {name!r} is not in the persisted-array allowlist")
        if origin.get("kind") == "diagnostic" and origin["array"].startswith("trace_"):
            if origin["array"][len("trace_") :] not in _trace_columns() and (
                origin["array"] not in DIAGNOSTIC_ARRAYS
            ):
                raise _fail(f"array {name!r} is not in the persisted-array allowlist")
        if name in seen:
            raise _fail(f"duplicate array name {name!r}")
        seen.add(name)
        if not _MEMBER_RE.fullmatch(name + ".npy"):
            raise _fail(f"unsafe array name {name!r}")
        if not _unit_ok(str(e.get("units")), origin):
            raise _fail(f"array {name!r}: units {e.get('units')!r} outside the vocabulary")


def _collect_arrays(result: Result, persist_batches: bool) -> list[_Arr]:
    out: list[_Arr] = []
    used_g: set[str] = set()
    for g in result.grids:
        gs = safe_name(g.name, used_g)
        used_q: set[str] = set()
        for base, (units, definition) in GRID_ARRAYS.items():
            if base == "batch_energy_mev" and not persist_batches:
                continue
            axes = (("batch",) if base == "batch_energy_mev" else ()) + _SPATIAL
            out.append(
                _Arr(
                    f"grid__{gs}__{base}",
                    np.ascontiguousarray(getattr(g, base)),
                    axes,
                    g.name,
                    units,
                    definition,
                    {"kind": "grid", "grid": g.name, "array": base},
                )
            )
        for qi, (qname, q) in enumerate(g.quantities.items()):
            qs = safe_name(qname, used_q)
            axes = _SPATIAL + (("energy_bin",) if np.ndim(q.mean) == 4 else ())
            for base in QUANTITY_ARRAYS:
                units = {"defined_mask": "bool", "n_nonzero": "count"}.get(base, q.units)
                out.append(
                    _Arr(
                        f"grid__{gs}__q__{qs}__{base}",
                        np.ascontiguousarray(getattr(q, base)),
                        axes,
                        g.name,
                        units,
                        f"{q.definition} [{base}]",
                        {
                            "kind": "quantity",
                            "grid": g.name,
                            "name": q.name,
                            "quantity": q.quantity,
                            "index": qi,
                            "array": base,
                        },
                    )
                )
    geom = result.effective_config.geometry
    if isinstance(geom, VoxelGeometry):
        for base in GEOMETRY_ARRAYS:
            arr = geom.material_index if base == "material_index" else geom.density_g_cm3
            if arr is None:
                continue
            units, definition = GEOMETRY_ARRAYS[base]
            out.append(
                _Arr(
                    f"geometry__{base}",
                    np.ascontiguousarray(arr),
                    _SPATIAL,
                    None,
                    units,
                    definition,
                    {"kind": "geometry", "array": base},
                )
            )
    out.extend(_diagnostic_arrays(result))
    return out


def _diagnostic_arrays(result: Result) -> list[_Arr]:
    out: list[_Arr] = []
    cols = _trace_columns()
    for key, val in result.diagnostics.items():
        if key in ("nuclear", "trace_columns"):
            continue
        if key == "trace" and isinstance(val, Mapping):
            for col, data in val.items():
                if col not in cols:
                    raise _fail(f"diagnostic trace column {col!r} is not in the allowlist")
                out.append(_diag(f"trace_{col}", data, "1", "trace column " + str(col)))
            continue
        if key not in DIAGNOSTIC_ARRAYS:
            raise _fail(f"diagnostic array {key!r} is not in the persisted-array allowlist")
        units, definition = DIAGNOSTIC_ARRAYS[key]
        out.append(_diag(key, val, units, definition))
    return out


def _diag(key: str, data: Any, units: str, definition: str) -> _Arr:
    a = np.ascontiguousarray(np.asarray(data))
    axes = ("record",) + tuple(f"c{i}" for i in range(a.ndim - 1))
    return _Arr(
        f"diag__{key}", a, axes, None, units, definition, {"kind": "diagnostic", "array": key}
    )


# ---------------------------------------------------------------------------------------------
# deterministic npz


def write_npz(path: Path, arrays: list[_Arr], zlib_level: int = ZLIB_LEVEL) -> None:
    """Write ``arrays`` (catalogue order) as a deterministic zip of ``.npy`` members."""
    with zipfile.ZipFile(
        path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=zlib_level, allowZip64=True
    ) as zf:
        for a in sorted(arrays, key=lambda x: x.name):
            zi = zipfile.ZipInfo(a.name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            zi._compresslevel = zlib_level  # type: ignore[attr-defined]
            zi.create_system = 3
            zi.external_attr = 0o600 << 16
            with zf.open(zi, "w", force_zip64=a.data.nbytes > 2**30) as fh:
                np.lib.format.write_array(fh, a.data, allow_pickle=False)


# ---------------------------------------------------------------------------------------------
# document blocks


def _strip(obj: Any, keys: frozenset[str]) -> Any:
    if isinstance(obj, Mapping):
        return {k: _strip(v, keys) for k, v in obj.items() if k not in keys}
    if isinstance(obj, list | tuple):
        return [_strip(v, keys) for v in obj]
    return obj


def _profile_axis(direction: Any) -> int:
    d = [abs(float(c)) for c in direction]
    return d.index(max(d))


def _material_def(m: Any) -> dict[str, Any]:
    st = m.sternheimer
    return {
        "name": m.name,
        "density_g_cm3": m.density_g_cm3,
        "mass_fractions": {k: float(v) for k, v in sorted(m.mass_fractions.items())},
        "I_eV_requested": m.I_eV,
        "I_eV_effective": float(math.exp(m.ln_I)),
        "sternheimer": None
        if st is None
        else {k: float(getattr(st, k)) for k in ("x0", "x1", "cbar", "a", "m")},
        "source": m.source,
    }


def _balance(result: Result) -> dict[str, Any]:
    eb = result.energy_balance
    d: dict[str, Any] = {
        f.name: (
            list(getattr(eb, f.name))
            if isinstance(getattr(eb, f.name), tuple)
            else getattr(eb, f.name)
        )
        for f in fields(eb)
        if f.name != "nuclear_mev"
    }
    d["closure_residual_mev"] = eb.closure_residual_mev
    if isinstance(eb, NuclearEnergyBalance):
        d["nuclear_mev"] = {k: float(v) for k, v in sorted(eb.nuclear_mev.items())}
    return d


def _grid_stats(result: Result) -> list[dict[str, Any]]:
    out = []
    for g in result.grids:
        m = np.asarray(g.defined_mask, dtype=bool)
        any_def = bool(m.any())
        why = "no voxel has positive mean energy and mass"
        dmax = float(np.max(g.dose_gy[m])) if any_def else None
        rmax = float(np.max(g.relative_uncertainty[m])) if any_def else None
        out.append(
            {
                "name": g.name,
                "valid": bool(g.valid),
                "low_batch_count": bool(g.low_batch_count),
                "shape": [int(n) for n in g.energy_mev.shape],
                "n_defined_voxels": int(m.sum()),
                "n_sparse_voxels": int(np.asarray(g.sparse_mask).sum()),
                "energy_total_mev": _scalar(float(np.sum(g.energy_mev)), "MeV"),
                "dose_max_gy": _scalar(dmax, "Gy", why),
                "relative_uncertainty_max": _scalar(rmax, "1", why),
            }
        )
    return out


def _data_block(result: Result) -> dict[str, Any]:
    eff = result.effective_config.summary()
    stopping = eff.get("stopping", {})
    nuc = eff.get("nuclear")
    return {
        "stopping_data_used": _strip(stopping.get("data_used", []), frozenset({"cache_dir"})),
        "stopping_content_sha256": stopping.get("content_sha256"),
        "table_sha256": eff["tables"]["sha256"],
        "table_stopping_source": eff["tables"]["stopping_source"],
        "table_materials": eff["tables"]["materials"],
        "nuclear": None if nuc is None else _strip(nuc, frozenset({"cache_dir"})),
        "lookups": [dict(lk) for lk in result.lookups],
    }


def _code_block(env: Mapping[str, Any]) -> dict[str, Any]:
    sha = env.get("git_sha")
    return {
        "ionmc_version": env.get("ionmc_version"),
        "git_commit": sha,
        "git_dirty": env.get("git_dirty"),
        "package_tree_oid": None,
        "source": "stamp" if sha else "unavailable",
        "status": "known" if sha else "unknown",
        "reason": "" if sha else "no build stamp (code identity builder is replaced in S3)",
    }


def _cache_dir_fields(stopping: Mapping[str, Any]) -> tuple[Any, str]:
    cd = stopping.get("cache_dir")
    if cd is None:
        return None, "not-used"
    if os.environ.get("IONMC_CACHE_DIR") == cd:
        return cd, "IONMC_CACHE_DIR"
    if str(Path.home() / ".cache" / "ionmc") == cd:
        return cd, "default"
    return cd, "explicit"


def build_document(
    result: Result,
    manifest: InputManifest | None,
    arrays: list[_Arr],
    persist_batches: bool,
    output_path: str,
) -> dict[str, Any]:
    """The result.json document without the ``files`` block."""
    eff = result.effective_config
    summary = eff.summary()
    cfg = result.requested_config
    manifest_dict = (manifest or manifest_for_api_config(cfg)).to_dict()
    env = result.environment
    stopping = summary.get("stopping", {})
    cache_dir, cache_source = _cache_dir_fields(stopping)
    counters = result.counters.as_dict()
    ch = result.channel_raw
    channel_counters = {
        "lookup_out_of_domain": 0 if ch is None else int(ch.lookup_out_of_domain),
        "path_bound_exceeded": 0 if ch is None else int(ch.path_bound_exceeded),
    }
    ph = summary["physics"]
    geom = eff.geometry
    by_origin = {
        (a.origin["kind"], a.origin["array"]): a.name
        for a in arrays
        if a.origin["kind"] == "geometry"
    }
    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "kind": KIND,
        "valid": bool(result.valid),
        "requested": to_dict(cfg, manifest=manifest),
        "input_manifest": manifest_dict,
        "effective": summary,
        "effective_extra": {
            "materials": [_material_def(m) for m in geom.materials],
            "geometry_arrays": {
                "material_index": by_origin.get(("geometry", "material_index")),
                "density_g_cm3": by_origin.get(("geometry", "density_g_cm3")),
            },
            "grid_convention": GRID_CONVENTION,
            "profile_axis": _profile_axis(summary["unit_direction"]),
            "resolutions": [],
            "resolutions_status": "not-computed: reconcile() is added in step S4",
        },
        "provenance": {
            "code": _code_block(env),
            "environment": {
                "python_version": env.get("python_version"),
                "platform": f"{platform.system()} {platform.release()} {platform.machine()}",
                "numpy_version": env.get("numpy_version"),
                "warp_version": env.get("warp_version"),
                "warp_devices": env.get("warp_devices"),
                "warp_error": env.get("warp_error"),
            },
            "data": _data_block(result),
            "models": {
                "stopping_model_id_effective": stopping.get("model_id_effective"),
                "stopping_version": stopping.get("version"),
                "stopping_spec_sha256": stopping.get("spec_sha256"),
                "evidence_at_run": stopping.get("evidence_at_run"),
                "straggling_model": ph["straggling_model"] if ph["straggling"] else None,
                "mcs_model": ph["mcs_model"] if ph["multiple_scattering"] else None,
                "end_of_range_model": ph["end_of_range_model_effective"],
                "delta_electrons": ph["delta_electrons"],
                "nuclear_table_id": (summary.get("nuclear") or {}).get("table_id"),
            },
        },
        "run": {
            "backend": result.backend,
            "device": result.device,
            "precision": result.precision,
            "production": summary["production"],
            "seed": result.seed,
            "rng": dict(result.rng),
            "n_histories": result.n_histories,
            "n_batches": result.n_batches,
            "histories_per_batch": result.n_histories // result.n_batches,
            "cpu_workers": summary["cpu_workers"],
            "chunk_histories": summary["chunk_histories"],
            "timings": {"deterministic": False, "values_s": dict(result.timings)},
            "output_path": output_path,
            "cache_dir": cache_dir,
            "cache_dir_source": cache_source,
        },
        "statistics": {
            **STAT_TEXT,
            "batch_arrays": (
                {"status": "persisted"}
                if persist_batches
                else {"status": "not-persisted", "reason": "user option persist_batches=false"}
            ),
            "grids": _grid_stats(result),
        },
        "accounting": {
            "energy_balance": _balance(result),
            "counters": counters,
            "channel_counters": channel_counters,
            "valid": bool(result.valid),
            "reasons": _invalidating(counters, channel_counters),
        },
        "units": dict(UNITS_TABLE),
        "arrays": [a.entry() for a in sorted(arrays, key=lambda x: x.name)],
    }


# ---------------------------------------------------------------------------------------------
# path hygiene (set P) and report placeholder

_P_PATTERNS = tuple(
    re.compile(p)
    for p in (
        r"^/run/output_path$",
        r"^/run/cache_dir$",
        r"^/input_manifest/base_dir$",
        r"^/input_manifest/config_file/path$",
        r"^/input_manifest/files/\d+/path$",
        r"^/effective/(stopping|nuclear)(/.*)?/cache_dir$",
    )
)


def _walk_strings(obj: Any, ptr: str = "") -> Iterator[tuple[str, str]]:
    if isinstance(obj, str):
        yield ptr, obj
    elif isinstance(obj, Mapping):
        for k, v in obj.items():
            yield from _walk_strings(v, f"{ptr}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk_strings(v, f"{ptr}/{i}")


def check_path_hygiene(doc: Mapping[str, Any]) -> None:
    """Outside the pointer set P no string may hold the home path, user name or host name."""
    names: set[str] = set()
    for fn in (getpass.getuser, socket.gethostname):
        with contextlib.suppress(Exception):
            names.add(fn())
    names.discard("")
    home = str(Path.home())
    for ptr, s in _walk_strings(doc):
        if any(p.match(ptr) for p in _P_PATTERNS):
            continue
        parts = set(re.split(r"[/\\]", s)) | {s}
        if (home not in ("", "/") and home in s) or (names & parts):
            raise _fail(f"string at {ptr} holds a home path, user name or host name")


def report_placeholder(doc: Mapping[str, Any]) -> str:
    """Minimal deterministic report (replaced by the full generator in step S5)."""
    run, code = doc["run"], doc["provenance"]["code"]
    lines = ["# ionmc result report (placeholder)", ""]
    if not doc["valid"]:
        lines += [
            "**INVALID RESULT: transport-limit counters are nonzero; arrays are not "
            "physically meaningful.**",
            "",
        ]
    lines += [
        f"- schema: {doc['schema']} {doc['schema_version']} ({doc['kind']})",
        f"- code version: {code['ionmc_version']}; commit: {code['git_commit'] or 'unknown'}",
        f"- backend: {run['backend']}, precision: {run['precision']}, seed: {run['seed']}",
        f"- histories: {run['n_histories']} in {run['n_batches']} batches",
        f"- arrays: {len(doc['arrays'])}",
    ]
    for g in doc["statistics"]["grids"]:
        lines.append(
            f"- grid {g['name']}: defined voxels {g['n_defined_voxels']}, valid {g['valid']}"
        )
    return "\n".join([*lines, "", "Not a clinical claim.", ""])


# ---------------------------------------------------------------------------------------------
# save


def _honoured(result: Result, arrays: list[_Arr]) -> None:
    cfg = result.requested_config
    have = {
        (a.origin.get("grid"), a.origin.get("name"))
        for a in arrays
        if a.origin["kind"] == "quantity"
    }
    for t in cfg.tallies:
        if (t.grid, t.name) not in have:
            raise _fail(f"requested tally {t.name!r} on grid {t.grid!r} left no array")
    diag = {a.origin["array"] for a in arrays if a.origin["kind"] == "diagnostic"}
    d = cfg.diagnostics
    for flag, key in (
        (d.track_end_positions, "end_position_mm"),
        (d.escape_records, "escape_history"),
    ):
        if flag and key not in diag:
            raise _fail(f"requested diagnostic {key!r} left no array")
    if d.trace_histories > 0 and "trace_end_history" not in diag:
        raise _fail("requested trace diagnostics left no array")
    if cfg.physics.nuclear != isinstance(result.energy_balance, NuclearEnergyBalance):
        raise _fail("nuclear request and nuclear energy balance disagree")


def save(
    result: Result,
    path: str | os.PathLike[str],
    *,
    manifest: InputManifest | None = None,
    persist_batches: bool = True,
    overwrite: bool = False,
    zlib_level: int = ZLIB_LEVEL,
) -> Path:
    """Write ``result`` as a bundle directory at ``path`` (temporary sibling + rename)."""
    dest = Path(path)
    parent = dest.parent
    if not parent.is_dir():
        raise _fail(f"parent directory {parent.name!r} does not exist")
    if os.path.lexists(dest):
        if not overwrite:
            raise FileExistsError(f"{dest.name!r} exists (pass overwrite=True)")
        if dest.is_symlink() or not (dest / "result.json").is_file():
            raise _fail("refusing to overwrite something that is not a result bundle")
    if manifest is not None:
        manifest.verify(result.requested_config)
    arrays = _collect_arrays(result, persist_batches)
    doc = build_document(result, manifest, arrays, persist_batches, str(dest))
    check_allowlist(doc["arrays"])
    _honoured(result, arrays)
    check_path_hygiene(doc)
    tmp = Path(tempfile.mkdtemp(dir=parent, prefix=".ionmc-tmp-"))
    try:
        write_npz(tmp / "arrays.npz", arrays, zlib_level)
        (tmp / "report.md").write_text(report_placeholder(doc), encoding="utf-8", newline="\n")
        doc["files"] = {}
        for name in ("arrays.npz", "report.md"):
            sha, size = _sha256_path(tmp / name)
            doc["files"][name] = {"sha256": sha, "size": size}
        (tmp / "result.json").write_bytes(canonical_bytes(doc))
        os.chmod(tmp, 0o755)
        old = None
        if os.path.lexists(dest):
            old = parent / f".ionmc-old-{os.getpid()}"
            os.rename(dest, old)
        os.rename(tmp, dest)
        if old is not None:
            shutil.rmtree(old)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    return dest


# ---------------------------------------------------------------------------------------------
# load: JSON pre-scan and structure


def prescan_json(data: bytes, max_tokens: int, max_depth: int) -> None:
    """Single pass over raw JSON bytes: every token (values, opening and closing delimiters,
    commas and colons) counts against one budget and the nesting depth is tracked. Allocation
    does not grow with the content (lazy regex iteration, early exit)."""
    tokens = depth = 0
    for m in _TOKEN_RE.finditer(data):
        c = m.group()[0]
        tokens += 1
        if c in b"{[":
            depth += 1
            if depth > max_depth:
                raise _fail(f"result.json nesting depth exceeds {max_depth}")
        elif c in b"}]":
            depth -= 1
        if tokens > max_tokens:
            raise _fail(f"result.json has more than {max_tokens} JSON tokens")


def _dup_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    d: dict[str, Any] = {}
    for k, v in pairs:
        if k in d:
            raise _fail(f"duplicate JSON key {k!r}")
        d[k] = v
    return d


def _parse_float(s: str) -> float:
    v = float(s)
    if not math.isfinite(v):
        raise _fail(f"non-finite JSON number {s}")
    return v


def _parse_constant(s: str) -> Any:
    raise _fail(f"JSON constant {s} is not allowed (NaN and infinity are never written)")


_TOP = (
    "schema schema_version kind valid requested input_manifest effective effective_extra "
    "provenance run statistics accounting units arrays files"
).split()
_RUN = (
    "backend device precision production seed rng n_histories n_batches histories_per_batch "
    "cpu_workers chunk_histories timings output_path cache_dir cache_dir_source"
).split()
_STAT = "normalisation uncertainty defined_mask sparse_mask batch_arrays grids".split()
_GRIDSTAT = (
    "name valid low_batch_count shape n_defined_voxels n_sparse_voxels energy_total_mev "
    "dose_max_gy relative_uncertainty_max"
).split()
_ACC = "energy_balance counters channel_counters valid reasons".split()
_PROV = "code environment data models".split()
_EXTRA = (
    "materials geometry_arrays grid_convention profile_axis resolutions resolutions_status".split()
)
_MANIFEST = "source base_dir config_file files config_sha256".split()
_ENTRY = "name dtype shape axes grid units definition content_sha256 origin".split()
_ORIGIN_KEYS = {"kind", "array", "grid", "name", "quantity", "index"}


def _keys(
    obj: Any, required: list[str], where: str, optional: tuple[str, ...] = ()
) -> dict[str, Any]:
    if not isinstance(obj, dict):
        raise _fail(f"{where}: expected an object")
    missing = [k for k in required if k not in obj]
    if missing:
        raise _fail(f"{where}: missing required field(s) {missing}")
    unknown = sorted(set(obj) - set(required) - set(optional))
    if unknown:
        raise _fail(f"{where}: unknown field(s) {unknown}")
    return obj


_UNAVAILABLE_KEYS = {"status", "reason"}
_MANIFEST_FILE = "role path sha256 content_sha256 item_sha256".split()


def _str_or_unavailable(v: Any, where: str) -> None:
    ok = isinstance(v, str) or (
        isinstance(v, dict) and set(v) == _UNAVAILABLE_KEYS and v["status"] == "unavailable"
    )
    if not ok:
        raise _fail(f"{where}: expected a string or an unavailable record")


def _is_unavailable(v: Any) -> bool:
    return isinstance(v, dict) and set(v) == _UNAVAILABLE_KEYS and v["status"] == "unavailable"


def _check_manifest(m: dict[str, Any]) -> None:
    """Strict keys and types of ``input_manifest`` (``InputManifest.to_dict`` forms)."""
    _need(m["source"] in ("toml", "api"), f"input_manifest.source: {m['source']!r} unknown")
    toml = m["source"] == "toml"
    cf = m["config_file"]
    if toml:
        _need(isinstance(m["base_dir"], str), "input_manifest.base_dir: toml needs a path")
        _need(
            isinstance(cf, dict) and set(cf) == {"path", "sha256"}
            and isinstance(cf["path"], str) and _hex(cf["sha256"]),
            "input_manifest.config_file: expected {path, sha256}",
        )  # fmt: skip
    else:
        _need(_is_unavailable(m["base_dir"]) and _is_unavailable(cf), "api manifest: unavailable")
    _need(_hex(m["config_sha256"]), "input_manifest.config_sha256 is malformed")
    _need(isinstance(m["files"], list), "input_manifest.files: expected a list")
    for i, f in enumerate(m["files"]):
        w = f"input_manifest.files[{i}]"
        _keys(f, _MANIFEST_FILE, w)
        _need(isinstance(f["role"], str), f"{w}.role: expected a string")
        _str_or_unavailable(f["path"], f"{w}.path")
        _str_or_unavailable(f["sha256"], f"{w}.sha256")
        _need(_is_unavailable(f["sha256"]) or _hex(f["sha256"]), f"{w}.sha256: not a digest")
        _need(_is_unavailable(f["path"]) == _is_unavailable(f["sha256"]), f"{w}: path/sha256")
        _need(_is_unavailable(f["path"]) == (not toml), f"{w}: availability disagrees with source")
        for k in ("content_sha256", "item_sha256"):
            _need(_hex(f[k]), f"{w}.{k} is malformed")


_DIAGNOSTIC_COUNTERS = frozenset({"mcs_table_bracket_miss"})  # never invalidates (simulation.py)
_COUNTER_KEYS = {
    "base": tuple(f.name for f in fields(TransportCounters)),
    "ion": tuple(f.name for f in fields(IonTransportCounters)),
    "nuclear": tuple(f.name for f in fields(NuclearTransportCounters)),
}
_CHANNEL_COUNTER_KEYS = ("lookup_out_of_domain", "path_bound_exceeded")


def _invalidating(counters: Mapping[str, int], channel: Mapping[str, int]) -> list[str]:
    """Names of the nonzero counters that invalidate a result (not the diagnostic ones)."""
    both = {**counters, **channel}
    return sorted(k for k, v in both.items() if v and k not in _DIAGNOSTIC_COUNTERS)


def _run_type(doc: Mapping[str, Any]) -> str:
    eff = doc["effective"]
    if "nuclear" in eff:
        return "nuclear"
    return "ion" if eff["source"]["z"] >= 2 else "base"


def _is_int(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


_BALANCE_LISTS = ("in_grid_mev", "outside_mev", "quantization_mev")
_BALANCE_SCALARS = (
    "initial_mev step_deposit_mev cutoff_mev escaped_mev truncated_mev unaccounted_mev "
    "closure_residual_mev"
).split()


def _is_num(v: Any) -> bool:
    return isinstance(v, int | float) and not isinstance(v, bool) and math.isfinite(v)


def _check_accounting(doc: dict[str, Any]) -> None:
    acc = doc["accounting"]
    rt = _run_type(doc)
    for blk, keys in (("counters", _COUNTER_KEYS[rt]), ("channel_counters", _CHANNEL_COUNTER_KEYS)):
        c = acc[blk]
        if not isinstance(c, dict) or set(c) != set(keys):
            raise _fail(f"accounting.{blk}: keys must be exactly {sorted(keys)} ({rt} run)")
        if not all(_is_int(v) and v >= 0 for v in c.values()):
            raise _fail(f"accounting.{blk}: values must be non-negative integers")
    if not isinstance(acc["valid"], bool):
        raise _fail("accounting.valid must be a boolean")
    if not (isinstance(acc["reasons"], list) and all(isinstance(r, str) for r in acc["reasons"])):
        raise _fail("accounting.reasons: expected a list of strings")
    eb = acc["energy_balance"]
    want = {*_BALANCE_SCALARS, *_BALANCE_LISTS} | ({"nuclear_mev"} if rt == "nuclear" else set())
    if not isinstance(eb, dict) or set(eb) != want:
        raise _fail(f"accounting.energy_balance: keys must be exactly {sorted(want)}")
    if not all(_is_num(eb[k]) for k in _BALANCE_SCALARS):
        raise _fail("accounting.energy_balance: scalar fields must be finite numbers")
    for k in _BALANCE_LISTS:
        if not (isinstance(eb[k], list) and all(_is_num(x) for x in eb[k])):
            raise _fail(f"accounting.energy_balance.{k}: expected a list of numbers")
    if rt == "nuclear" and not (
        isinstance(eb["nuclear_mev"], dict) and all(_is_num(x) for x in eb["nuclear_mev"].values())
    ):
        raise _fail("accounting.energy_balance.nuclear_mev: expected numbers")


_CACHE_SOURCES = ("explicit", "IONMC_CACHE_DIR", "default", "not-used")
_MATERIAL_KEYS = (
    "name density_g_cm3 mass_fractions I_eV_requested I_eV_effective sternheimer source"
).split()
_SCALAR_FLOAT = ("value", "units")
_SCALAR_UNDEF = ("value", "status", "reason", "units")


def _need(cond: bool, msg: str) -> None:
    if not cond:
        raise _fail(msg)


def _is_str_list(v: Any) -> bool:
    return isinstance(v, list) and all(isinstance(x, str) for x in v)


def _hex(v: Any) -> bool:
    return isinstance(v, str) and bool(_SHA_RE.fullmatch(v))


def _check_scalar_record(r: Any, where: str) -> None:
    _need(isinstance(r, dict) and r.get("units") in UNIT_VOCABULARY, f"{where}: malformed record")
    if "status" in r:
        _need(set(r) == set(_SCALAR_UNDEF), f"{where}: malformed undefined record")
        _need(r["status"] == "undefined" and r["value"] is None, f"{where}: malformed undefined")
        _need(isinstance(r["reason"], str), f"{where}: reason must be a string")
    else:
        _need(set(r) == set(_SCALAR_FLOAT) and _is_num(r["value"]), f"{where}: malformed value")


def _check_run(r: dict[str, Any]) -> None:
    for k in ("backend", "device", "precision", "output_path"):
        _need(isinstance(r[k], str), f"run.{k}: expected a string")
    for k in ("seed", "n_histories", "n_batches", "histories_per_batch", "cpu_workers"):
        _need(_is_int(r[k]) and r[k] >= 0, f"run.{k}: expected a non-negative integer")
    _need(r["n_batches"] >= 1, "run.n_batches must be positive")
    _need(_is_int(r["chunk_histories"]), "run.chunk_histories: expected an integer")
    _need(isinstance(r["production"], bool), "run.production: expected a boolean")
    _need(isinstance(r["rng"], dict), "run.rng: expected an object")
    _need(r["cache_dir"] is None or isinstance(r["cache_dir"], str), "run.cache_dir: type")
    _need(r["cache_dir_source"] in _CACHE_SOURCES, "run.cache_dir_source: unknown value")
    t = r["timings"]
    _need(isinstance(t, dict) and set(t) == {"deterministic", "values_s"}, "run.timings: keys")
    _need(t["deterministic"] is False, "run.timings.deterministic must be false")
    _need(
        isinstance(t["values_s"], dict) and all(_is_num(v) for v in t["values_s"].values()),
        "run.timings.values_s: expected numbers",
    )


def _check_statistics(st: dict[str, Any]) -> None:
    for k in ("normalisation", "uncertainty", "defined_mask", "sparse_mask"):
        _need(isinstance(st[k], str), f"statistics.{k}: expected a string")
    ba = st["batch_arrays"]
    _need(isinstance(ba, dict), "statistics.batch_arrays: expected an object")
    if ba.get("status") == "persisted":
        _need(set(ba) == {"status"}, "statistics.batch_arrays: unknown field(s)")
    else:
        _need(
            set(ba) == {"status", "reason"} and ba["status"] == "not-persisted"
            and isinstance(ba["reason"], str),
            "statistics.batch_arrays: malformed record",
        )  # fmt: skip
    for i, g in enumerate(st["grids"]):
        w = f"statistics.grids[{i}]"
        _need(isinstance(g, dict), f"{w}: expected an object")
        _need(isinstance(g["name"], str), f"{w}.name: expected a string")
        for k in ("valid", "low_batch_count"):
            _need(isinstance(g[k], bool), f"{w}.{k}: expected a boolean")
        _need(
            isinstance(g["shape"], list) and all(_is_int(n) and n >= 0 for n in g["shape"]),
            f"{w}.shape: expected integers",
        )
        for k in ("n_defined_voxels", "n_sparse_voxels"):
            _need(_is_int(g[k]) and g[k] >= 0, f"{w}.{k}: expected a non-negative integer")
        for k in ("energy_total_mev", "dose_max_gy", "relative_uncertainty_max"):
            _check_scalar_record(g[k], f"{w}.{k}")


def _check_extra(x: dict[str, Any]) -> None:
    _need(isinstance(x["materials"], list), "effective_extra.materials: expected a list")
    for i, m in enumerate(x["materials"]):
        w = f"effective_extra.materials[{i}]"
        _keys(m, _MATERIAL_KEYS, w)
        _need(isinstance(m["name"], str) and isinstance(m["source"], str), f"{w}: strings")
        _need(_is_num(m["density_g_cm3"]) and _is_num(m["I_eV_effective"]), f"{w}: numbers")
        _need(m["I_eV_requested"] is None or _is_num(m["I_eV_requested"]), f"{w}.I_eV_requested")
        _need(
            isinstance(m["mass_fractions"], dict)
            and all(_is_num(v) for v in m["mass_fractions"].values()),
            f"{w}.mass_fractions: expected numbers",
        )
        _need(m["sternheimer"] is None or isinstance(m["sternheimer"], dict), f"{w}.sternheimer")
    ga = x["geometry_arrays"]
    _need(
        isinstance(ga, dict) and set(ga) == set(GEOMETRY_ARRAYS)
        and all(v is None or isinstance(v, str) for v in ga.values()),
        "effective_extra.geometry_arrays: malformed",
    )  # fmt: skip
    _need(x["grid_convention"] == GRID_CONVENTION, "effective_extra.grid_convention differs")
    _need(_is_int(x["profile_axis"]) and 0 <= x["profile_axis"] <= 2, "profile_axis: 0, 1 or 2")
    _need(
        isinstance(x["resolutions"], list) and all(isinstance(r, dict) for r in x["resolutions"]),
        "effective_extra.resolutions: expected a list of objects",
    )
    _need(isinstance(x["resolutions_status"], str), "effective_extra.resolutions_status: type")


def _schema(value: Any, spec: Any, path: str) -> None:
    try:
        shapes.check(value, spec, path)
    except shapes.SchemaError as exc:
        raise _fail(str(exc)) from exc


def _check_effective(eff: dict[str, Any]) -> None:
    """``effective`` against the declared schema of ``EffectiveConfig.summary()`` (every nested
    field typed; exact integers, finite numbers, list lengths)."""
    try:
        shapes.check_effective(eff)
    except shapes.SchemaError as exc:
        raise _fail(str(exc)) from exc


def _check_provenance(p: dict[str, Any]) -> None:
    for blk, spec in (
        ("code", shapes.PROVENANCE_CODE),
        ("environment", shapes.PROVENANCE_ENVIRONMENT),
        ("data", shapes.PROVENANCE_DATA),
        ("models", shapes.PROVENANCE_MODELS),
    ):
        _schema(p[blk], spec, f"provenance.{blk}")
    _need(p["code"]["source"] in ("stamp", "git", "unavailable"), "provenance.code.source")


def _bind(cond: bool, what: str) -> None:
    _need(cond, f"{what} disagrees with requested")


def _check_run_binding(doc: dict[str, Any]) -> None:
    """Bind run, provenance and statistics scalars to ``requested`` (the configuration) and to
    ``effective``: a field of ``run`` or ``effective`` that is changed without ``requested`` is
    refused."""
    eff, run, rq = doc["effective"], doc["run"], doc["requested"]
    rrun = rq["run"]
    for k in ("backend", "precision", "seed", "n_histories", "n_batches", "cpu_workers",
              "chunk_histories"):  # fmt: skip
        _bind(run[k] == rrun[k], f"run.{k}")
        _bind(eff[k] == rrun[k], f"effective.{k}")
    _need(run["production"] == eff["production"], "run.production disagrees with effective")
    _need(run["rng"] == eff["rng"], "run.rng disagrees with effective.rng")
    _need(run["n_histories"] >= 1, "run.n_histories must be positive")
    _need(run["n_histories"] % run["n_batches"] == 0, "n_histories is not a multiple of n_batches")
    _need(run["histories_per_batch"] == run["n_histories"] // run["n_batches"],
          "run.histories_per_batch disagrees with n_histories / n_batches")  # fmt: skip
    dev = run["device"].lower()
    _need(
        ("cuda" in dev) == (run["backend"] == "warp-cuda") or run["backend"] == "python",
        "run.device disagrees with the requested backend",
    )
    ms = rrun["max_steps"]
    _bind(eff["max_steps_origin"] == ("user" if ms is not None else "computed"),
          "effective.max_steps_origin")  # fmt: skip
    if ms is not None:
        _bind(eff["max_steps"] == ms, "effective.max_steps")
    _need(eff["max_steps"] >= 1 and eff["scoring_pieces"] >= 1, "effective step bounds")
    cd = run["cache_dir"]
    _need(cd == (eff["stopping"] or {}).get("cache_dir"), "run.cache_dir disagrees with effective")
    _need((cd is None) == (run["cache_dir_source"] == "not-used"), "run.cache_dir_source")
    p = doc["provenance"]
    st, nuc = eff["stopping"], eff.get("nuclear")
    ph = eff["physics"]
    m = p["models"]
    want_models = {
        "stopping_model_id_effective": st["model_id_effective"],
        "stopping_version": st["version"],
        "stopping_spec_sha256": st["spec_sha256"],
        "evidence_at_run": st["evidence_at_run"],
        "straggling_model": ph["straggling_model"] if ph["straggling"] else None,
        "mcs_model": ph["mcs_model"] if ph["multiple_scattering"] else None,
        "end_of_range_model": ph["end_of_range_model_effective"],
        "delta_electrons": ph["delta_electrons"],
        "nuclear_table_id": None if nuc is None else nuc["table_id"],
    }
    _need(m == want_models, "provenance.models disagrees with effective")
    d = p["data"]
    want_data = {
        "stopping_data_used": _strip(st["data_used"], frozenset({"cache_dir"})),
        "stopping_content_sha256": st["content_sha256"],
        "table_sha256": eff["tables"]["sha256"],
        "table_stopping_source": eff["tables"]["stopping_source"],
        "table_materials": eff["tables"]["materials"],
        "nuclear": None if nuc is None else _strip(nuc, frozenset({"cache_dir"})),
    }
    _need({k: d[k] for k in want_data} == want_data, "provenance.data disagrees with effective")
    _need(len(d["lookups"]) == len(rq.get("lookups", [])),
          "provenance.data.lookups disagrees with requested.lookups")  # fmt: skip
    for pl, rl in zip(d["lookups"], rq.get("lookups", []), strict=True):
        _need(pl["units"] == rl.get("units") and pl["content_sha256"] == rl["content_sha256"],
              "provenance.data.lookups disagrees with requested.lookups")  # fmt: skip
    _need(doc["effective_extra"]["profile_axis"] == _profile_axis(eff["unit_direction"]),
          "effective_extra.profile_axis disagrees with the source direction")  # fmt: skip
    _need(doc["accounting"]["valid"] == doc["valid"], "accounting.valid disagrees with valid")


_PHYSICS_BOUND = (
    "nuclear straggling straggling_model multiple_scattering mcs_model delta_electrons e_cut_mev "
    "max_step_mm max_energy_loss_fraction range_alpha range_rho_f_mm short_step_fraction "
    "truncated_hinge_diagnostic e_cut_ion_mev_per_u end_of_range_model"
).split()


def _check_effective_binding(doc: dict[str, Any]) -> None:
    """Bind the nested ``effective`` block to ``requested`` and check value semantics (positive
    bounds, finite geometry): physics options, source, geometry, materials, scoring grids and
    the compiled channel plan."""
    eff, rq = doc["effective"], doc["requested"]
    ph, rph = eff["physics"], rq["physics"]
    for k in _PHYSICS_BOUND:
        _bind(ph[k] == rph[k], f"effective.physics.{k}")
    _bind(("nuclear" in eff) == bool(rph["nuclear"]), "effective.nuclear presence")
    if "nuclear" in eff:
        _bind(ph["e_cut_deuteron_mev"] == rph["e_cut_deuteron_mev"], "physics.e_cut_deuteron_mev")
        _need(ph["e_cut_deuteron_mev"] > 0, "effective.physics.e_cut_deuteron_mev must be positive")
    smid = rph["stopping_model"]
    if isinstance(smid, str):
        _bind(eff["stopping"]["model_id_requested"] == smid, "effective.stopping model id")
    for k in ("e_cut_mev", "e_cut_ion_mev_per_u", "e_cut_effective_mev", "max_step_mm",
              "range_alpha", "range_rho_f_mm", "short_step_fraction"):  # fmt: skip
        _need(ph[k] > 0, f"effective.physics.{k} must be positive")
    _need(0 < ph["max_energy_loss_fraction"] <= 1, "max_energy_loss_fraction must be in (0, 1]")
    # source
    src, rs = eff["source"], rq["source"]
    rp = rs["projectile"]
    if isinstance(rp, str):
        reg = PROJECTILES[rp]
        want = {"projectile": rp, "z": reg.z, "a": reg.a, "mass_mev": reg.mass_mev}
    else:
        want = {
            k: rp[{"projectile": "name"}.get(k, k)] for k in ("projectile", "z", "a", "mass_mev")
        }
    for k, v in want.items():
        _bind(src[k] == v, f"effective.source.{k}")
    for k in ("position_mm", "kinetic_energy_mev", "energy_sigma_mev", "lateral_sigma_mm"):
        _bind(src[k] == rs[k], f"effective.source.{k}")
    _need(src["kinetic_energy_mev"] > 0 and src["energy_sigma_mev"] >= 0, "source energies")
    _need(src["lateral_sigma_mm"] >= 0 and src["z"] >= 1 and src["a"] >= 1, "source values")
    d, rd = eff["unit_direction"], rs["direction"]
    nrm = math.sqrt(sum(float(c) ** 2 for c in rd))
    _need(nrm > 0, "requested.source.direction is zero")
    _need(all(abs(a - float(b) / nrm) <= 1e-9 for a, b in zip(d, rd, strict=True)),
          "effective.unit_direction is not the normalised requested direction")  # fmt: skip
    # geometry and materials
    geo, rg = eff["geometry"], rq["geometry"]
    if rg["type"] == "voxel":
        want_g = {k: rg[k] for k in ("origin_mm", "spacing_mm", "shape", "materials", "z_exit_mm")}
    else:
        want_g = {"origin_mm": rg["lower_mm"], "spacing_mm": rg["size_mm"], "shape": [1, 1, 1],
                  "materials": [rg["material"]], "z_exit_mm": None}  # fmt: skip
    _bind(geo == want_g, "effective.geometry")
    names = {m["name"] for m in rq.get("materials", [])}
    _need(set(geo["materials"]) <= names and bool(geo["materials"]),
          "effective.geometry.materials are not defined in requested.materials")  # fmt: skip
    _need(all(x > 0 for x in geo["spacing_mm"]) and all(n >= 1 for n in geo["shape"]),
          "effective.geometry: spacing and shape must be positive")  # fmt: skip
    _need([m["name"] for m in doc["effective_extra"]["materials"]] == geo["materials"],
          "effective_extra.materials disagrees with effective.geometry.materials")  # fmt: skip
    for g in eff["scoring"]:
        _need(all(x > 0 for x in g["spacing_mm"]) and all(n >= 1 for n in g["shape"]),
              "effective.scoring: spacing and shape must be positive")  # fmt: skip
    _check_channel_plan(doc)
    try:
        validate_identity_block(eff["stopping"])
    except (ValueError, KeyError, TypeError, IndexError) as exc:
        raise _fail(f"effective.stopping identity block is inconsistent: {exc}") from exc
    # statistics definitions and units
    st = doc["statistics"]
    for k, v in STAT_TEXT.items():
        _need(st[k] == v, f"statistics.{k} differs from the engine definition")
    nb = doc["run"]["n_batches"]
    for g in st["grids"]:
        _need(g["valid"] == doc["valid"], "statistics.grids.valid disagrees with valid")
        _need(g["low_batch_count"] == (nb < LOW_BATCH_COUNT), "statistics low_batch_count")
        for k, u in (("energy_total_mev", "MeV"), ("dose_max_gy", "Gy"),
                     ("relative_uncertainty_max", "1")):  # fmt: skip
            _need(g[k]["units"] == u, f"statistics.grids.{k}: units must be {u!r}")
    _check_energy_structure(doc)


def _check_energy_structure(doc: dict[str, Any]) -> None:
    """Accounting identities that need no array (decision 0045 amendment A-1, item 4): per grid
    ``in_grid + quantization + outside = step_deposit + cutoff + nuclear_local`` and the closure
    ``initial - (step + cutoff + escaped + truncated + unaccounted + signed nuclear terms) =
    closure_residual_mev``; tolerance ``1e-9 max(1, |reference|)`` (sums in different float
    order); a nuclear result persists exactly the engine's nuclear key set."""
    eb, n = doc["accounting"]["energy_balance"], len(doc["effective"]["scoring"])
    for k in ("in_grid_mev", "outside_mev"):
        _need(len(eb[k]) == n, f"accounting.energy_balance.{k}: one entry per scoring grid")
    _need(len(eb["quantization_mev"]) in (0, n), "energy_balance.quantization_mev length")
    nuc = eb.get("nuclear_mev")
    nuc_dest = local = 0.0
    if "nuclear_mev" in eb:
        _need(
            isinstance(nuc, dict) and set(nuc) == set(NUCLEAR_TALLY_NAMES)
            and all(_is_num(v) for v in nuc.values()),
            f"energy_balance.nuclear_mev: keys must be exactly {sorted(NUCLEAR_TALLY_NAMES)}",
        )  # fmt: skip
        assert nuc is not None
        local = nuc["nuclear_local"]
        nuc_dest = (nuc["nuclear_local"] + nuc["nuclear_escaped_neutron"]
                    + nuc["nuclear_escaped_gamma"] + nuc["nuclear_binding"]
                    + nuc["nuclear_imbalance"])  # fmt: skip
    total = sum(eb[k] for k in ("step_deposit_mev", "cutoff_mev", "escaped_mev", "truncated_mev",
                                "unaccounted_mev")) + nuc_dest  # fmt: skip
    tol = 1e-9 * max(1.0, abs(eb["initial_mev"]))
    resid = eb["initial_mev"] - total - eb["closure_residual_mev"]
    _need(abs(resid) <= tol, "energy_balance: closure_residual_mev disagrees with identity")
    grid_total = eb["step_deposit_mev"] + eb["cutoff_mev"] + local
    gtol = 1e-9 * max(1.0, abs(grid_total))
    for i in range(n):
        quant = eb["quantization_mev"][i] if eb["quantization_mev"] else 0.0
        _need(eb["in_grid_mev"][i] >= 0 and eb["outside_mev"][i] >= 0,
              f"energy_balance: negative in-grid or outside energy of grid {i}")  # fmt: skip
        _need(abs(eb["in_grid_mev"][i] + quant + eb["outside_mev"][i] - grid_total) <= gtol,
              f"energy_balance: in_grid + quantization + outside != deposits, grid {i}"
              )  # fmt: skip


# ---- the channel plan derived from the REQUESTED tallies (no value of effective.tallies is used)

_QRULES = {  # quantity -> (numerator kind, class, denominator, units, definition)
    "fluence": ("L", CLASS_STEP, None, "mm", "sum of path lengths"),
    "let_t": ("LS", CLASS_STEP, ("L", CLASS_STEP), "MeV/mm",
              "sum(l S)/sum(l), S = unrestricted electronic stopping power in water"),
    "let_d": ("LS2", CLASS_STEP, ("LS", CLASS_STEP), "MeV/mm", "sum(l S^2)/sum(l S), S in water"),
    "let_d_eps": ("ES", CLASS_STEP, ("E", CLASS_STEP), "MeV/mm",
                  "sum(eps S)/sum(eps) over step deposits, S in water"),
}  # fmt: skip


def _derive_channel_plan(
    doc: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int, int]:
    """``(channels, quantities, total_size, n_residual_columns)`` that ``compile_channels`` builds
    from the requested tallies, grids and lookup units: the channel specs, their deduplication
    and order, the automatic channels of every used grid, the sizes, offsets and the quantity
    definitions. The numerical quanta and bounds (``k``, ``bound_per_history``) are not derived."""
    rq = doc["requested"]
    grids = rq["scoring"]
    gi_of = {g["name"]: i for i, g in enumerate(grids)}
    lnames = [x["name"] for x in doc["provenance"]["data"]["lookups"]]
    lunits = [r["units"] for r in rq.get("lookups", [])]
    specs: dict[tuple[Any, ...], dict[str, Any]] = {}

    def get(kind: str, gi: int, cls: int, species: Any, gen: tuple[int, int], li: int = -1,
            sp: Any = None, pseudo: bool = False) -> tuple[Any, ...]:  # fmt: skip
        row = _species_row(None if species is None else tuple(species), include_pseudo=pseudo)
        key = (kind, gi, cls, row, gen, li, None if sp is None else (sp.edges, sp.log))
        specs.setdefault(key, {"kind": kind, "grid": gi, "class_mask": cls, "row": row,
                               "gen": gen, "li": li, "sp": sp})  # fmt: skip
        return key

    qdefs: list[tuple[dict[str, Any], str, tuple[Any, ...], Any, str, str]] = []
    used: set[int] = set()
    both = CLASS_STEP | CLASS_LOCAL
    for t in rq.get("tallies", []):
        gi = gi_of[t["grid"]]
        used.add(gi)
        gen, q, sps = _generation_range(t["generation"]), t["quantity"], t["species"]
        li, lu, sp = -1, "", None
        if t["lookup"] is not None:
            li = lnames.index(t["lookup"])
            lu = lunits[li]
        if t["energy_edges_mev_per_u"] is not None:
            sp = _spectrum_spec(t["name"], tuple(t["energy_edges_mev_per_u"]))

        def spec(kind: str, cls: int, li: int = -1, sp: Any = None, _g: int = gi,
                 _gen: tuple[int, int] = gen, _s: Any = sps) -> tuple[Any, ...]:  # fmt: skip
            pseudo = kind == "E" and bool(cls & CLASS_LOCAL) and _s is None
            return get(kind, _g, cls, _s, _gen, li, sp, pseudo=pseudo)

        if q in ("edep", "dose"):
            qdefs.append((t, "linear", spec("E", both), None, "MeV", f"sum of deposits ({q})"))
        elif q in _QRULES:
            nk, nc, den, units, defn = _QRULES[q]
            qdefs.append((t, "ratio" if den else "linear", spec(nk, nc),
                          None if den is None else spec(*den), units, defn))  # fmt: skip
        elif q == "lookup_sum":
            qdefs.append((t, "linear", spec("FE", CLASS_STEP, li=li), None, f"MeV*[{lu}]",
                          "sum(eps f)"))  # fmt: skip
        elif q == "lookup_dose_avg":
            qdefs.append((t, "ratio", spec("FE", CLASS_STEP, li=li), spec("E", CLASS_STEP), lu,
                          "sum(eps f)/sum(eps) over step deposits"))  # fmt: skip
        elif q == "fluence_spectrum":
            qdefs.append((t, "linear", spec("FL", CLASS_STEP, sp=sp), None, "mm",
                          "sum of path lengths per energy bin per nucleon (under/overflow "
                          "included)"))  # fmt: skip
        else:
            raise _fail(f"requested tally {t['name']!r}: unknown quantity {q!r}")
    for gi in sorted(used):
        get("N", gi, CLASS_STEP, None, (0, GEN_MAX))
        get("N", gi, CLASS_LOCAL, None, (0, GEN_MAX), pseudo=True)
        get("E", gi, CLASS_LOCAL, None, (0, GEN_MAX), pseudo=True)
    ordered = sorted(specs.items(), key=lambda kv: kv[1]["grid"])
    channels: list[dict[str, Any]] = []
    index: dict[tuple[Any, ...], int] = {}
    offset = n_res = 0
    for key, s in ordered:
        shp = grids[s["grid"]]["shape"]
        nvox = int(shp[0]) * int(shp[1]) * int(shp[2])
        sp = s["sp"]
        size = nvox * (sp.n_bins + 2 if sp is not None else 1)
        index[key] = len(channels)
        channels.append({
            "kind": s["kind"], "grid": s["grid"], "class_mask": s["class_mask"],
            "species": [i for i, v in enumerate(s["row"]) if v], "generation": list(s["gen"]),
            "lookup": s["li"],
            "spectrum": None if sp is None else {"edges": list(sp.edges), "log": sp.log},
            "size": size, "offset": offset,
        })  # fmt: skip
        offset += size
        n_res += s["kind"] != "N"
    quantities = [
        {"name": t["name"], "kind": kind, "numerator": index[num],
         "denominator": None if den is None else index[den], "units": units, "definition": defn}
        for t, kind, num, den, units, defn in qdefs
    ]  # fmt: skip
    return channels, quantities, offset, n_res


def _check_channel_plan(doc: dict[str, Any]) -> None:
    """``effective.tallies`` must EQUAL the plan derived from ``requested`` (its channels except
    the numerical quanta and bounds, its quantities, sizes, residual columns, histories per batch
    and memory); without requested tallies it must be null."""
    eff, rq = doc["effective"], doc["requested"]
    tl = eff["tallies"]
    if not rq.get("tallies"):
        _need(tl is None, "effective.tallies must be null without requested tallies")
        return
    _need(tl is not None, "effective.tallies is null although tallies were requested")
    assert tl is not None
    channels, quantities, total, n_res = _derive_channel_plan(doc)
    got = tl["channels"]
    _need(len(got) == len(channels), "effective.tallies.channels differ from the requested plan")
    for c, w in zip(got, channels, strict=True):
        core = {k: v for k, v in c.items() if k not in ("k", "bound_per_history")}
        _need(core == w, "effective.tallies.channels differ from the plan from requested.tallies")
        # k may be negative (a coarse quantum is legal): only the bound is signed
        _need(c["bound_per_history"] >= 0, "channel bound is negative")
    _need(
        tl["quantities"] == quantities,
        "effective.tallies.quantities differ from the plan from requested.tallies",
    )
    run = doc["run"]
    nvox = sum(int(g["shape"][0]) * int(g["shape"][1]) * int(g["shape"][2]) for g in rq["scoring"])
    want = {
        "n_channels": len(channels), "total_size": total, "n_residual_columns": n_res,
        "histories_per_batch": run["n_histories"] // run["n_batches"],
        "memory_bytes": run["n_batches"] * (nvox + total) * 8 * run["cpu_workers"],
    }  # fmt: skip
    _need({k: tl[k] for k in want} == want, "effective.tallies sizes differ from the plan")
    _need(tl["lookups"] == doc["provenance"]["data"]["lookups"],
          "effective.tallies.lookups differ from provenance.data.lookups")  # fmt: skip


def _check_structure(doc: Any) -> None:
    _keys(doc, _TOP, "result.json")
    if doc["schema"] != SCHEMA:
        raise _fail(f"schema {doc['schema']!r} is not {SCHEMA!r}")
    if doc["schema_version"] not in KNOWN_VERSIONS:
        raise _fail(f"unknown schema_version {doc['schema_version']!r}; known: {KNOWN_VERSIONS}")
    if doc["kind"] != KIND:
        raise _fail(f"unknown kind {doc['kind']!r}")
    if not isinstance(doc["valid"], bool):
        raise _fail("valid must be a boolean")
    _check_run(_keys(doc["run"], _RUN, "run"))
    st = _keys(doc["statistics"], _STAT, "statistics")
    if not isinstance(st["grids"], list):
        raise _fail("statistics.grids: expected a list")
    for i, g in enumerate(st["grids"]):
        _keys(g, _GRIDSTAT, f"statistics.grids[{i}]")
    _check_statistics(st)
    _keys(doc["accounting"], _ACC, "accounting")
    for k in ("requested", "effective", "units", "provenance"):
        if not isinstance(doc[k], dict):
            raise _fail(f"{k}: expected an object")
    _keys(doc["provenance"], _PROV, "provenance")
    _need(doc["units"] == UNITS_TABLE, "units: differs from the unit table")
    _check_extra(_keys(doc["effective_extra"], _EXTRA, "effective_extra"))
    try:
        validate_dict_structure(doc["requested"])
    except ConfigError as exc:
        raise _fail(f"requested: {exc}") from exc
    _check_effective(doc["effective"])
    _check_provenance(doc["provenance"])
    _check_run_binding(doc)
    _check_effective_binding(doc)
    _check_accounting(doc)
    _check_manifest(_keys(doc["input_manifest"], _MANIFEST, "input_manifest"))
    if not isinstance(doc["arrays"], list):
        raise _fail("arrays: expected a list")
    for i, e in enumerate(doc["arrays"]):
        _keys(e, _ENTRY, f"arrays[{i}]")
        if not isinstance(e["origin"], dict) or not set(e["origin"]) <= _ORIGIN_KEYS:
            raise _fail(f"arrays[{i}].origin: unknown field(s)")
        for k in ("name", "dtype", "units", "definition"):
            if not isinstance(e[k], str):
                raise _fail(f"arrays[{i}].{k}: expected a string")
        if not (isinstance(e["axes"], list) and all(isinstance(a, str) for a in e["axes"])):
            raise _fail(f"arrays[{i}].axes: expected a list of strings")
        if e["grid"] is not None and not isinstance(e["grid"], str):
            raise _fail(f"arrays[{i}].grid: expected a string or null")
        for k, v in e["origin"].items():
            if not (_is_int(v) if k == "index" else isinstance(v, str)):
                raise _fail(f"arrays[{i}].origin.{k}: wrong type")
        if not (isinstance(e["shape"], list) and all(_is_int(n) for n in e["shape"])):
            raise _fail(f"arrays[{i}].shape: expected a list of integers")
        if not _hex(e["content_sha256"]):
            raise _fail(f"arrays[{i}].content_sha256 is malformed")
    files = doc["files"]
    if not isinstance(files, dict):
        raise _fail("files: expected an object")
    for name, f in files.items():
        _keys(f, ["sha256", "size"], f"files[{name!r}]")
        if not _hex(f["sha256"]):
            raise _fail(f"files[{name!r}].sha256 is malformed")
        if not _is_int(f["size"]) or f["size"] < 0:
            raise _fail(f"files[{name!r}].size is malformed")
    _check_manifest_binding(doc)


def _array_positions(doc: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    """(role, requested record) of every array position of ``requested``, in position order."""
    req = doc["requested"]
    g = req["geometry"]
    recs: list[tuple[str, dict[str, Any]]] = []
    if g["type"] == "voxel":
        recs.append(("material_index", g["material_index"]))
        if g.get("density") is not None:
            recs.append(("density", g["density"]))
    names = [x["name"] for x in doc["provenance"]["data"]["lookups"]]
    for n, r in zip(names, req.get("lookups", []), strict=True):
        recs.append((f"lookup:{n}", r))
    return recs


def _check_manifest_binding(doc: dict[str, Any]) -> None:
    """Bind ``input_manifest`` to the persisted ``requested`` document (no file is opened):
    config digest, exact unique roles for every array position (whatever the record kind),
    reference and content digests per role."""
    req, m = doc["requested"], doc["input_manifest"]
    if sha256_canonical(req) != m["config_sha256"]:
        raise _fail("input_manifest.config_sha256 does not match the requested document")
    _need(
        len(doc["provenance"]["data"]["lookups"]) == len(req.get("lookups", [])),
        "provenance.data.lookups disagrees with requested.lookups",
    )
    recs = _array_positions(doc)
    files = m["files"]
    roles = [f["role"] for f in files]
    _need(len(set(roles)) == len(roles), "input_manifest.files repeats a role")
    _need(
        sorted(roles) == sorted(r for r, _ in recs),
        f"input_manifest roles {sorted(roles)} differ from the roles implied by requested",
    )
    by_role = {f["role"]: f for f in files}
    for role, r in recs:
        f = by_role[role]
        _need(
            f["content_sha256"] == r.get("content_sha256"), f"role {role!r}: content digests differ"
        )
        if "ref" in r:
            _need(r["ref"] == role, f"requested reference {r['ref']!r} has the wrong role {role!r}")
            _need(
                f["sha256"] == r["sha256"],
                f"requested reference {role!r} disagrees with the manifest entry",
            )
        else:
            _need(
                not isinstance(f["sha256"], str),
                f"role {role!r}: expected a reference",
            )
    refs = [r["ref"] for _, r in recs if "ref" in r]
    for f in files:
        if isinstance(f["sha256"], str):
            _need(f["role"] in refs, f"manifest file {f['role']!r} is not referenced")


# ---------------------------------------------------------------------------------------------
# load: files and ZIP


_PLOT_RE = re.compile(r"^plots/[a-z0-9_-]{1,64}\.svg$")
_FIXED_FILES = ("arrays.npz", "report.md")


def _file_name_ok(name: Any) -> bool:
    return isinstance(name, str) and (name in _FIXED_FILES or bool(_PLOT_RE.fullmatch(name)))


def _bundle_file(root: Path, name: str) -> Path:
    if not _file_name_ok(name):
        raise _fail(f"unsafe file name {name!r} in the files map")
    p = root / name
    try:
        st = os.lstat(p)
    except OSError as exc:
        raise _fail(f"file {name!r} listed in the bundle is missing") from exc
    import stat as _stat

    if not _stat.S_ISREG(st.st_mode):
        raise _fail(f"{name!r} is not a regular file (symlinks are refused)")
    real_root = os.path.realpath(root)
    if os.path.commonpath([real_root, os.path.realpath(p)]) != real_root:
        raise _fail(f"{name!r} resolves outside the bundle")
    return p


def _zip_preflight(path: Path, max_entries: int, max_total: int, max_cd: int) -> None:
    size = path.stat().st_size
    if size > max_total:
        raise _fail(f"arrays.npz ({size} bytes) exceeds max_total_bytes")
    with open(path, "rb") as fh:
        tail_len = min(size, _EOCD_WINDOW)
        fh.seek(size - tail_len)
        tail = fh.read(tail_len)
        i = tail.rfind(b"PK\x05\x06")
        if i < 0 or len(tail) - i < 22:
            raise _fail("arrays.npz has no end-of-central-directory record")
        _sig, _d, _cd, n_this, n_total, cd_size, cd_off, _cl = struct.unpack(
            "<4s4H2LH", tail[i : i + 22]
        )
        if 0xFFFF in (n_this, n_total) or 0xFFFFFFFF in (cd_size, cd_off):
            eocd_pos = size - tail_len + i
            if eocd_pos < 20:
                raise _fail("arrays.npz: truncated ZIP64 locator")
            fh.seek(eocd_pos - 20)
            loc = fh.read(20)
            if loc[:4] != b"PK\x06\x07":
                raise _fail("arrays.npz: ZIP64 values without a ZIP64 locator")
            z64_off = struct.unpack("<4sLQL", loc)[2]
            if z64_off + 56 > size:
                raise _fail("arrays.npz: ZIP64 record offset outside the file")
            fh.seek(z64_off)
            rec = fh.read(56)
            if len(rec) < 56 or rec[:4] != b"PK\x06\x06":
                raise _fail("arrays.npz: malformed ZIP64 end record")
            n_total, cd_size, cd_off = struct.unpack("<QQQ", rec[32:56])
    if n_total > max_entries:
        raise _fail(
            f"arrays.npz declares {n_total} entries, more than the {max_entries} catalogued"
        )
    if cd_size > max_cd:
        raise _fail(f"arrays.npz central directory ({cd_size} bytes) exceeds the budget")
    if cd_off > size or cd_off + cd_size > size:
        raise _fail("arrays.npz central directory lies outside the file")


def _entry_nbytes(e: dict[str, Any], max_array: int) -> tuple[np.dtype[Any], tuple[int, ...], int]:
    if e["dtype"] not in _ALLOWED_DTYPES:
        raise _fail(f"array {e['name']!r}: dtype {e['dtype']!r} is not allowed")
    dt = np.dtype(e["dtype"])
    shape = tuple(e["shape"])
    if len(shape) > 4 or any(n < 0 for n in shape):
        raise _fail(f"array {e['name']!r}: invalid shape {shape}")
    n = 1
    for d in shape:
        n *= d
    nbytes = n * dt.itemsize
    if nbytes > max_array:
        raise _fail(f"array {e['name']!r}: {nbytes} bytes exceed max_array_bytes")
    return dt, shape, nbytes


def _read_member(
    zf: zipfile.ZipFile, info: zipfile.ZipInfo, e: dict[str, Any], dt: Any, shape: Any, nbytes: int
) -> np.ndarray[Any, Any]:
    name = e["name"]
    if info.file_size > nbytes + MAX_NPY_HEADER_BYTES:
        raise _fail(f"member {name!r} declares {info.file_size} bytes, over its catalogue budget")
    try:
        with zf.open(info) as fh:
            ver = np.lib.format.read_magic(fh)
            if ver == (1, 0):
                hshape, fortran, hdt = np.lib.format.read_array_header_1_0(fh, MAX_NPY_HEADER_BYTES)
            elif ver == (2, 0):
                hshape, fortran, hdt = np.lib.format.read_array_header_2_0(fh, MAX_NPY_HEADER_BYTES)
            else:
                raise _fail(f"member {name!r}: unsupported NPY version {ver}")
            if fortran or tuple(hshape) != shape or hdt.str != dt.str:
                raise _fail(f"member {name!r}: NPY header disagrees with the catalogue")
            arr = np.empty(shape, dtype=dt)
            view = memoryview(arr.reshape(-1).view(np.uint8)) if nbytes else memoryview(b"")
            got = 0
            while got < nbytes:
                chunk = fh.read(min(nbytes - got, 1 << 22))
                if not chunk:
                    raise _fail(f"member {name!r} is shorter than its catalogue size")
                view[got : got + len(chunk)] = chunk
                got += len(chunk)
            if fh.read(1):
                raise _fail(f"member {name!r} inflates beyond its declared size")
    except (zipfile.BadZipFile, zlib.error, ValueError, EOFError, OSError) as exc:
        raise _fail(f"member {name!r} is corrupt: {exc}") from exc
    return arr


# ---------------------------------------------------------------------------------------------
# load: cross-field invariants


_UNIT_OVERRIDE = {"dose": "Gy", "fluence": "mm^-2", "fluence_spectrum": "mm^-2"}
_PER_VOLUME = ("fluence", "fluence_spectrum")
_EXCLUDED_DEFINITION = "energy of cutoff and zero-length deposits, excluded from LET and lookups"
_PIECE_DEFINITIONS = {
    "scoring_pieces": "scoring pieces of class step per primary",
    "scoring_pieces_local": "scoring pieces of class local (point deposits) per primary",
}


_F8, _B1, _I4 = {"<f8"}, {"|b1"}, {"<i4", "<i8"}
_GRID_DTYPES = {
    "defined_mask": _B1, "sparse_mask": _B1, "n_nonzero_batches": _I4,
}  # fmt: skip
_QUANTITY_DTYPES = {"mean": _F8, "std": _F8, "rounding_bound": _F8, "defined_mask": _B1}


@dataclass(frozen=True)
class _Expect:
    name: str
    dtypes: set[str]
    shape: list[int] | None
    axes: list[str] | None
    units: str
    quantity: str | None = None
    index: int | None = None
    definition: str | None = None


_QPlan = tuple[str, str, str, list[int], str]  # name, quantity, units, shape, definition


def _quantity_plan(
    plan_q: Mapping[str, Mapping[str, Any]], req: Mapping[str, Any], gi: int, shape: list[int]
) -> list[_QPlan]:
    """Ordered quantities of grid ``gi`` the run must have: the requested tallies (units and
    definitions from the plan derived from ``requested``, never from ``effective.tallies``), the
    excluded-energy channel and the two piece counts that every grid with tallies has."""
    grid = req["scoring"][gi]["name"]
    out: list[_QPlan] = []
    for t in req.get("tallies", []):
        if t["grid"] != grid:
            continue
        q = t["quantity"]
        shp = list(shape)
        if q == "fluence_spectrum":
            shp.append(len(t["energy_edges_mev_per_u"]) + 1)
        pq = plan_q[t["name"]]
        defn = pq["definition"] + (" per volume" if q in _PER_VOLUME else "")
        out.append((t["name"], q, _UNIT_OVERRIDE.get(q, pq["units"]), shp, defn))
    if not out:
        return out
    out.append(("edep_excluded_from_let", "edep", "MeV", list(shape), _EXCLUDED_DEFINITION))
    for n, d in _PIECE_DEFINITIONS.items():
        out.append((n, "count", "1", list(shape), d))
    return out


def _expected_catalogue(doc: dict[str, Any], gs: dict[str, str]) -> dict[tuple[Any, ...], _Expect]:
    """Every array the bundle must hold, derived from ``requested``/``effective`` only."""
    req, eff = doc["requested"], doc["effective"]
    nb = doc["run"]["n_batches"]
    persisted = doc["statistics"]["batch_arrays"]["status"] == "persisted"
    plan_q = {q["name"]: q for q in _derive_channel_plan(doc)[1]}
    exp: dict[tuple[Any, ...], _Expect] = {}
    for gi, g in enumerate(eff["scoring"]):
        n, shape = g["name"], list(g["shape"])
        for base, (units, _d) in GRID_ARRAYS.items():
            if base == "batch_energy_mev" and not persisted:
                continue
            batch = base == "batch_energy_mev"
            exp[("grid", n, None, base)] = _Expect(
                f"grid__{gs[n]}__{base}",
                _GRID_DTYPES.get(base, _F8),
                ([nb] if batch else []) + shape,
                (["batch"] if batch else []) + list(_SPATIAL),
                units,
            )
        plan = _quantity_plan(plan_q, req, gi, shape)
        used: set[str] = set()
        for qi, (qn, q, qunits, shp, qdef) in enumerate(plan):
            qs = safe_name(qn, used)
            axes = list(_SPATIAL) + (["energy_bin"] if len(shp) == 4 else [])
            for base in QUANTITY_ARRAYS:
                qu: str = {"defined_mask": "bool", "n_nonzero": "count"}.get(base, qunits)
                dt = _I4 if base == "n_nonzero" else _QUANTITY_DTYPES[base]
                exp[("quantity", n, qn, base)] = _Expect(
                    f"grid__{gs[n]}__q__{qs}__{base}", dt, shp, axes, qu, q, qi, f"{qdef} [{base}]"
                )
    geom = eff["geometry"]
    gshape = list(geom["shape"])
    exp[("geometry", None, None, "material_index")] = _Expect(
        "geometry__material_index", _I4, gshape, list(_SPATIAL), "1"
    )
    if req["geometry"].get("density") is not None:
        exp[("geometry", None, None, "density_g_cm3")] = _Expect(
            "geometry__density_g_cm3", _F8, gshape, list(_SPATIAL), "g/cm3"
        )
    d = req.get("diagnostics", {})
    want: set[str] = set()
    if d.get("track_end_positions"):
        want |= {k for k in DIAGNOSTIC_ARRAYS if k.startswith(("end_", "control_"))}
    if d.get("escape_records"):
        want |= {k for k in DIAGNOSTIC_ARRAYS if k.startswith("escape_")}
    if d.get("trace_histories", 0) > 0:
        want |= {k for k in DIAGNOSTIC_ARRAYS if k.startswith("trace_")}
        want |= {f"trace_{c}" for c in _trace_columns()}
    for k in want:
        exp[("diagnostic", None, None, k)] = _Expect(
            f"diag__{k}", set(_ALLOWED_DTYPES), None, None, DIAGNOSTIC_ARRAYS.get(k, ("1", ""))[0]
        )
    return exp


def _grid_meta(g: Mapping[str, Any]) -> tuple[Any, ...]:
    return (g.get("name"), g.get("origin_mm"), g.get("spacing_mm"), g.get("shape"))


def _check_catalogue(doc: dict[str, Any]) -> dict[tuple[Any, Any, Any, Any], str]:
    """Compare the catalogue to the exact expectation derived from ``requested``/``effective``
    (names, dtype, shape, axes, units, quantity metadata) and require unique origins."""
    sc, rsc = doc["effective"]["scoring"], doc["requested"]["scoring"]
    if not (isinstance(sc, list) and all(isinstance(g, dict) for g in sc)):
        raise _fail("effective.scoring: expected a list of objects")
    for g in sc:
        if not (isinstance(g.get("name"), str) and isinstance(g.get("shape"), list)):
            raise _fail("effective.scoring: malformed grid entry")
        if len(g["shape"]) != 3 or not all(_is_int(n) for n in g["shape"]):
            raise _fail("effective.scoring: malformed grid shape")
    gnames = [g["name"] for g in sc]
    if len(set(gnames)) != len(gnames):
        raise _fail("effective.scoring: duplicate grid names")
    if [_grid_meta(g) for g in rsc] != [_grid_meta(g) for g in sc]:
        raise _fail("effective.scoring disagrees with requested.scoring")
    geo = doc["effective"]["geometry"]
    if not (
        isinstance(geo, dict) and isinstance(geo.get("shape"), list) and len(geo["shape"]) == 3
        and all(_is_int(n) for n in geo["shape"])
    ):  # fmt: skip
        raise _fail("effective.geometry.shape: expected three integers")
    rg = doc["requested"]["geometry"]
    if rg["type"] == "voxel" and list(rg["shape"]) != list(geo["shape"]):
        raise _fail("effective.geometry.shape disagrees with requested.geometry.shape")
    entries = doc["arrays"]
    if [e["name"] for e in entries] != sorted(e["name"] for e in entries):
        raise _fail("arrays: catalogue must be in sorted member-name order")
    by: dict[tuple[Any, Any, Any, Any], str] = {}
    qidx: dict[str, dict[str, int]] = {}
    for e in entries:
        o = e["origin"]
        key = (o.get("kind"), o.get("grid"), o.get("name"), o.get("array"))
        if key in by:  # reachable: the allowlist pattern does not bind the grid part of a name
            raise _fail(f"arrays: duplicate semantic origin {key!r}")
        by[key] = e["name"]
        if len(e["axes"]) != len(e["shape"]) and o.get("kind") != "diagnostic":
            raise _fail(f"array {e['name']!r}: axes and shape lengths differ")
        if e["grid"] != o.get("grid"):
            raise _fail(f"array {e['name']!r}: grid field and origin disagree")
        if o.get("kind") == "quantity":
            if not _is_int(o.get("index")) or not isinstance(o.get("name"), str):
                raise _fail(f"array {e['name']!r}: malformed quantity origin")
            if qidx.setdefault(o.get("grid"), {}).setdefault(o["name"], o["index"]) != o["index"]:
                raise _fail(f"quantity {o['name']!r}: inconsistent index")
    used: set[str] = set()
    gs = {n: safe_name(n, used) for n in gnames}
    exp = _expected_catalogue(doc, gs)
    if set(exp) != set(by):
        missing = sorted(set(exp) - set(by), key=str)
        extra = sorted(set(by) - set(exp), key=str)
        raise _fail(f"array inventory differs from the request: missing {missing}, extra {extra}")
    for e in entries:
        o = e["origin"]
        x = exp[(o["kind"], o.get("grid"), o.get("name"), o["array"])]
        bad = (
            e["name"] != x.name
            or e["dtype"] not in x.dtypes
            or (x.shape is not None and list(e["shape"]) != x.shape)
            or (x.axes is not None and list(e["axes"]) != x.axes)
            or (x.units is not None and e["units"] != x.units)
            or (x.quantity is not None and o.get("quantity") != x.quantity)
            or (x.index is not None and o.get("index") != x.index)
            or (x.definition is not None and e["definition"] != x.definition)
        )
        if bad:
            raise _fail(f"array {e['name']!r} does not match its semantic origin and the request")
    for g, names in qidx.items():
        if sorted(names.values()) != list(range(len(names))):
            raise _fail(f"grid {g!r}: quantity indices are not 0..n-1")
    ga = doc["effective_extra"]["geometry_arrays"]
    for base in GEOMETRY_ARRAYS:
        if ga.get(base) != by.get(("geometry", None, None, base)):
            raise _fail(f"effective_extra.geometry_arrays.{base} disagrees with the catalogue")
    return by


_RTOL_DERIVED = 1e-12  # arrays that are the engine's own elementwise transform of another
_RTOL_BATCH = 1e-9  # mean and batch standard error of the persisted batch arrays
_RATIO_QUANTITIES = ("let_t", "let_d", "let_d_eps", "lookup_dose_avg")


def _close(a: np.ndarray[Any, Any], b: np.ndarray[Any, Any], rtol: float, what: str) -> None:
    """Equal NaN pattern and equal values to ``rtol`` (relative to the largest finite value)."""
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    ok = a.shape == b.shape and bool(np.array_equal(nan_a, nan_b))
    if ok:
        x, y = a[~nan_a], b[~nan_b]
        scale = float(np.max(np.abs(y))) if y.size else 0.0
        ok = bool(np.all(np.abs(x - y) <= rtol * np.maximum(np.abs(y), 1e-3 * scale)))
    _need(ok, f"{what} disagrees with the arrays it is derived from")


def _check_array_invariants(
    doc: dict[str, Any], arrays: dict[str, np.ndarray[Any, Any]], by: Mapping[Any, str]
) -> None:
    """Array invariants of decision 0045 amendment A-1, item 5 (tolerances: derived arrays
    ``1e-12``, batch statistics ``1e-9``, in-grid energy ``1e-9`` relative): standard errors finite
    and non-negative; dose, dose standard error and relative uncertainty are the engine's
    transforms of ``energy_mev``/``energy_std_mev``/``mass_g``; ``defined_mask``,
    ``n_nonzero_batches`` and ``sparse_mask`` agree with the energy and the batch count; the batch
    arrays reproduce mean and batch standard error; quantity companions follow the decision-0040
    rules; ``in_grid_mev`` equals the energy array times the number of histories."""
    nb, nh = doc["run"]["n_batches"], doc["run"]["n_histories"]
    persisted = doc["statistics"]["batch_arrays"]["status"] == "persisted"
    eb = doc["accounting"]["energy_balance"]
    for i, g in enumerate(doc["effective"]["scoring"]):
        n = g["name"]

        def get(base: str, n: str = n) -> np.ndarray[Any, Any]:
            return arrays[by[("grid", n, None, base)]]

        en, std, mass = get("energy_mev"), get("energy_std_mev"), get("mass_g")
        dose, dstd, rel = get("dose_gy"), get("dose_std_gy"), get("relative_uncertainty")
        mask, nnz, sparse = (
            get("defined_mask").astype(bool),
            get("n_nonzero_batches"),
            get("sparse_mask").astype(bool),
        )
        _need(bool(np.all(np.isfinite(en)) and np.all(en >= 0)), f"grid {n!r}: energy_mev")
        _need(bool(np.all(np.isfinite(std)) and np.all(std >= 0)), f"grid {n!r}: energy_std_mev")
        _need(bool(np.all(np.isfinite(mass)) and np.all(mass >= 0)), f"grid {n!r}: mass_g")
        pos = mass > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            scale = np.where(pos, MEV_PER_G_TO_GY / np.where(pos, mass, 1.0), np.nan)
            _close(dose, np.where(pos, en * scale, np.nan), _RTOL_DERIVED, f"grid {n!r}: dose_gy")
            _close(dstd, np.where(pos, std * scale, np.nan), _RTOL_DERIVED,
                   f"grid {n!r}: dose_std_gy")  # fmt: skip
            _close(rel, np.where(en > 0, std / np.where(en > 0, en, 1.0), np.nan),
                   _RTOL_DERIVED, f"grid {n!r}: relative_uncertainty")  # fmt: skip
        _need(bool(np.array_equal(mask, (en > 0) & pos)), f"grid {n!r}: defined_mask")
        _need(bool(np.all((nnz >= 0) & (nnz <= nb))), f"grid {n!r}: n_nonzero_batches range")
        _need(bool(np.array_equal(nnz > 0, en > 0)), f"grid {n!r}: n_nonzero_batches vs energy")
        _need(bool(np.array_equal(sparse, mask & (nnz < 0.5 * nb))), f"grid {n!r}: sparse_mask")
        if persisted:
            x = get("batch_energy_mev").reshape(nb, -1)
            mean = x.mean(axis=0)
            sem = np.sqrt(((x - mean) ** 2).sum(axis=0) / (nb * (nb - 1)))
            _close(mean, en.reshape(-1), _RTOL_BATCH, f"grid {n!r}: batch mean")
            _close(sem, std.reshape(-1), _RTOL_BATCH, f"grid {n!r}: batch standard error")
            _need(bool(np.array_equal((x > 0).sum(axis=0), nnz.reshape(-1))),
                  f"grid {n!r}: n_nonzero_batches vs batch arrays")  # fmt: skip
        dep = float(np.sum(en)) * nh
        _need(abs(eb["in_grid_mev"][i] - dep) <= 1e-9 * max(1.0, abs(dep)),
              f"accounting.energy_balance.in_grid_mev[{i}] disagrees with the energy array of "
              f"grid {n!r}")  # fmt: skip
    ent = {e["name"]: e for e in doc["arrays"]}
    for (kind, gname, qname, base), name in by.items():
        if kind != "quantity" or base != "mean":
            continue
        quant = ent[name]["origin"]["quantity"]
        mean, qstd, qmask, qn, bound = (
            arrays[by[("quantity", gname, qname, b)]]
            for b in ("mean", "std", "defined_mask", "n_nonzero", "rounding_bound")
        )
        qmask = qmask.astype(bool)
        w = f"quantity {qname!r} of grid {gname!r}"
        _need(bool(np.all((qn >= 0) & (qn <= nb))), f"{w}: n_nonzero range")
        if quant in _RATIO_QUANTITIES:
            for a in (mean, qstd, bound):
                _need(bool(np.all(np.isnan(a[~qmask]))), f"{w}: undefined ratio values are NaN")
            _need(bool(np.all(qn[qmask] >= min_defined_batches(nb))),
                  f"{w}: defined bins need max(2, ceil(B/2)) nonzero batches")  # fmt: skip
        else:
            _need(bool(np.array_equal(qmask, mean > 0)), f"{w}: defined_mask is mean > 0")
        _need(bool(np.all(np.isfinite(mean[qmask])) and np.all(np.isfinite(qstd[qmask]))),
              f"{w}: defined values must be finite")  # fmt: skip
        _need(bool(np.all(qstd[qmask] >= 0) and np.all(bound[qmask] >= 0)),
              f"{w}: standard error and rounding bound are non-negative")  # fmt: skip
        _need(not bool(np.any(np.isinf(qstd))), f"{w}: infinite standard error")


# ---- re-derivation of ``effective`` (decision 0045 amendment A-1, items 1 and 2)

_F_MIN, _F_MAX = math.ldexp(1.0, -60), math.ldexp(1.0, 60)  # channels.FE_F_{MIN,MAX}_EXPONENT
LOOKUP_PLAN_EXCLUSIONS = (
    "effective.tallies.channels[kind == 'FE'].k",
    "effective.tallies.channels[kind == 'FE'].bound_per_history",
    "effective.tallies.bounds.f_max",
    "effective.tallies.lookups[*].content_sha256",
)
REDERIVE_EXCLUSIONS = (
    "effective.stopping.evidence_at_run (persisted verbatim, never recomputed)",
    "every cache_dir (a path of the writing machine)",
)


def _placeholder_lookup(p: Mapping[str, Any], fmax: float) -> LookupTable:
    lo, hi = p["axis_range"]
    n = int(p["n_points"])
    axis = np.geomspace(lo, hi, n) if p["axis_spacing"] == "log" else np.linspace(lo, hi, n)
    return LookupTable(
        name=p["name"], quantity=p["quantity"], units=p["units"], axis=p["axis"],
        axis_spacing=p["axis_spacing"], axis_values=axis,
        values={sp: np.full(n, fmax) for sp in p["species"]}, citation=p["citation"],
        license=p["license"], source=p["source"], synthetic=p["synthetic"],
        file_sha256=p["file_sha256"], resampling=p["resampling"],
    )  # fmt: skip


def _comparable(eff: Mapping[str, Any], lookup: bool) -> dict[str, Any]:
    e = _strip(eff, frozenset({"cache_dir"}))
    e["stopping"].pop("evidence_at_run", None)
    if lookup and e["tallies"] is not None:
        for c in e["tallies"]["channels"]:
            if c["kind"] == "FE":
                c.pop("k"), c.pop("bound_per_history")
        e["tallies"]["bounds"].pop("f_max", None)
        for lk in e["tallies"]["lookups"]:
            lk.pop("content_sha256", None)
    return e  # type: ignore[no-any-return]


def _rederive(
    doc: dict[str, Any], arrays: dict[str, np.ndarray[Any, Any]], by: Mapping[Any, str],
    rederive: bool,
) -> tuple[str, tuple[str, ...]]:  # fmt: skip
    """Rebuild the configuration from ``requested`` (geometry arrays from the persisted members,
    lookup tables as placeholders of the persisted provenance), ``validate`` it and compare its
    ``summary()`` with the persisted ``effective`` (canonical JSON), except the documented
    exclusions. Returns the binding level and the recorded limits."""
    if not rederive:
        return "typed-only", ("effective is typed and bound to requested, not re-derived",)
    rq, eff = doc["requested"], doc["effective"]
    user = rq["physics"].get("stopping")
    if isinstance(user, dict) and "status" in user:  # a user stopping object cannot be rebuilt
        return "typed-only", (
            "requested.physics.stopping is a not-reconstructible user stopping object: effective "
            "is typed and bound to requested, not re-derived",
        )
    lookups = doc["provenance"]["data"]["lookups"]
    geo = {
        role: arrays[by[("geometry", None, None, base)]]
        for role, base in (("material_index", "material_index"), ("density", "density_g_cm3"))
        if ("geometry", None, None, base) in by
    }

    def not_derivable(exc: BaseException) -> ResultFormatError:
        cached = isinstance(exc, DatasetError) or isinstance(
            exc.__cause__, FileNotFoundError | OSError
        )
        what = "dataset not cached" if cached else "the requested configuration is refused"
        return _fail(f"effective not re-derivable: {what}: {exc}")

    try:
        cfg = from_dict({k: v for k, v in rq.items() if k != "lookups"}, arrays=geo)
        if lookups:
            fmax = (eff["tallies"] or {}).get("bounds", {}).get("f_max")
            _need(
                isinstance(fmax, float) and _F_MIN <= fmax <= _F_MAX,
                "effective.tallies.bounds.f_max is outside the admissible lookup range",
            )
            cfg = replace(cfg, lookups=tuple(_placeholder_lookup(p, fmax) for p in lookups))
            for pl, p in zip(cfg.lookups, lookups, strict=True):
                got = {k: v for k, v in pl.provenance().items() if k != "content_sha256"}
                _need(got == {k: v for k, v in p.items() if k != "content_sha256"},
                      f"lookup {p['name']!r}: provenance is not a table description")  # fmt: skip
        want = validate(cfg).summary()
    except (DatasetError, UnsupportedCombinationError) as exc:
        raise not_derivable(exc) from exc
    if canonical_bytes(_comparable(want, bool(lookups))) != canonical_bytes(
        _comparable(eff, bool(lookups))
    ):
        raise _fail("effective differs from the summary re-derived from requested")
    if lookups:
        return "re-derived-except-lookup-plan", LOOKUP_PLAN_EXCLUSIONS
    return "re-derived", ()


def _check_invariants(
    doc: dict[str, Any], arrays: dict[str, np.ndarray[Any, Any]], rederive: bool = True
) -> tuple[str, tuple[str, ...]]:
    scoring = {g["name"]: tuple(g["shape"]) for g in doc["effective"].get("scoring", [])}
    nb = doc["run"]["n_batches"]
    by = _check_catalogue(doc)
    for e in doc["arrays"]:
        if e["grid"] is not None:
            if e["grid"] not in scoring:
                raise _fail(f"array {e['name']!r}: unknown grid {e['grid']!r}")
            shp, gshape = tuple(e["shape"]), scoring[e["grid"]]
            if e["axes"] and e["axes"][0] == "batch":
                if len(shp) != 4 or shp[0] != nb or shp[1:] != gshape:
                    raise _fail(f"array {e['name']!r}: shape {shp} disagrees with grid {gshape}")
            elif shp[:3] != gshape:
                raise _fail(f"array {e['name']!r}: shape {shp} disagrees with grid {gshape}")
    stats = {g["name"]: g for g in doc["statistics"]["grids"]}
    if set(stats) != set(scoring):
        raise _fail("statistics.grids disagree with the effective scoring grids")
    for gname in scoring:

        def get(base: str, gname: str = gname) -> np.ndarray[Any, Any]:
            n = by.get(("grid", gname, None, base))
            if n is None:
                raise _fail(f"grid {gname!r}: missing array {base!r}")
            return arrays[n]

        mask, mass = get("defined_mask").astype(bool), get("mass_g")
        rel, dose, dstd, en = (
            get(b) for b in ("relative_uncertainty", "dose_gy", "dose_std_gy", "energy_mev")
        )
        bad = ~mask & (mass > 0) & ~np.isnan(rel)
        if bad.any() or (mask & ~(np.isfinite(rel) & np.isfinite(dose) & (en > 0))).any():
            raise _fail(f"grid {gname!r}: undefined voxels must be NaN, never 0 (defined_mask)")
        if ((mass <= 0) & ~(np.isnan(dose) & np.isnan(dstd))).any():
            raise _fail(f"grid {gname!r}: dose must be NaN where the voxel mass is 0")
        s = stats[gname]
        for key, arr in (("dose_max_gy", dose), ("relative_uncertainty_max", rel)):
            rec = s[key]
            if mask.any():
                want = float(np.max(arr[mask]))
                if rec.get("status") is not None or rec.get("value") != want:
                    raise _fail(f"statistics.grids[{gname!r}].{key} disagrees with the arrays")
            elif rec.get("status") != "undefined" or rec.get("value") is not None:
                raise _fail(f"statistics.grids[{gname!r}].{key} must be an undefined record")
        if s["n_defined_voxels"] != int(mask.sum()):
            raise _fail(f"statistics.grids[{gname!r}].n_defined_voxels disagrees with the mask")
    for (kind, gname, qname, base), name in by.items():
        if kind == "quantity" and base == "mean":
            q = next((e for e in doc["arrays"] if e["name"] == name), None)
            if q is not None and q["origin"].get("quantity") not in (
                None,
                "dose",
                "edep",
                "fluence",
                "fluence_spectrum",
                "count",
                "lookup_sum",
            ):
                m = arrays[by[("quantity", gname, qname, "defined_mask")]].astype(bool)
                if (~m & ~np.isnan(arrays[name])).any():
                    raise _fail(f"quantity {qname!r}: ratio values must be NaN where undefined")
    for g, st in zip(doc["effective"]["scoring"], doc["statistics"]["grids"], strict=True):
        n = g["name"]
        _need(
            st["name"] == n and st["shape"] == g["shape"],
            f"statistics.grids[{n!r}] shape or name differs",
        )
        en = arrays[by[("grid", n, None, "energy_mev")]]
        _need(st["energy_total_mev"].get("value") == float(np.sum(en)),
              f"statistics.grids[{n!r}].energy_total_mev disagrees with the array")  # fmt: skip
        sp = arrays[by[("grid", n, None, "sparse_mask")]]
        _need(
            st["n_sparse_voxels"] == int(sp.sum()),
            f"statistics.grids[{n!r}].n_sparse_voxels differs",
        )
    _check_array_invariants(doc, arrays, by)
    geo = {
        role: arrays[by[("geometry", None, None, base)]]
        for role, base in (("material_index", "material_index"), ("density", "density_g_cm3"))
        if ("geometry", None, None, base) in by
    }
    items = reader_item_sha256s(doc["requested"], geo, doc["provenance"]["data"]["lookups"])
    for f in doc["input_manifest"]["files"]:
        _need(
            items.get(f["role"]) == f["item_sha256"],
            f"manifest item_sha256 of {f['role']!r} differs",
        )
    acc = doc["accounting"]
    nonzero = _invalidating(acc["counters"], acc["channel_counters"])
    if doc["valid"] != (not nonzero) or acc["valid"] != doc["valid"]:
        raise _fail("valid disagrees with the transport-limit counters")
    if acc["reasons"] != nonzero:
        raise _fail("accounting.reasons disagree with the counters")
    req = doc["requested"]
    have = {(o[1], o[2]) for o in by if o[0] == "quantity"}
    for t in req.get("tallies", []):
        if (t.get("grid"), t.get("name")) not in have:
            raise _fail(f"requested tally {t.get('name')!r} has no array in the bundle")
    want_diag = {
        "track_end_positions": "end_position_mm",
        "escape_records": "escape_history",
    }
    present = {k[3] for k in by if k[0] == "diagnostic"}
    d = req.get("diagnostics", {})
    for flag, key in want_diag.items():
        if bool(d.get(flag)) != (key in present):
            raise _fail(f"diagnostic {key!r} presence disagrees with the request")
    if (d.get("trace_histories", 0) > 0) != ("trace_end_history" in present):
        raise _fail("trace diagnostics presence disagrees with the request")
    return _rederive(doc, arrays, by, rederive)


# ---------------------------------------------------------------------------------------------
# load


@dataclass(frozen=True)
class LoadedResult:
    """A verified bundle: the JSON document and every array (read-only)."""

    path: Path
    doc: dict[str, Any]
    arrays: dict[str, np.ndarray[Any, Any]] = field(repr=False)
    effective_binding: str = "typed-only"
    binding_limits: tuple[str, ...] = ()
    rederive_exclusions: tuple[str, ...] = REDERIVE_EXCLUSIONS

    def grid_array(self, grid: str, base: str) -> np.ndarray[Any, Any]:
        for e in self.doc["arrays"]:
            o = e["origin"]
            if o["kind"] == "grid" and o.get("grid") == grid and o["array"] == base:
                return self.arrays[e["name"]]
        raise KeyError(f"no array {base!r} on grid {grid!r}")

    def dump(self) -> dict[str, Any]:
        """Canonical dump (document, content digests, file digests) for clean-process checks."""
        return {
            "doc": self.doc,
            "content_sha256": {n: array_content_sha256(a) for n, a in sorted(self.arrays.items())},
        }


def load(
    path: str | os.PathLike[str],
    *,
    max_total_bytes: int = MAX_TOTAL_BYTES,
    max_array_bytes: int = MAX_ARRAY_BYTES,
    max_json_bytes: int = MAX_JSON_BYTES,
    max_central_directory_bytes: int = MAX_CENTRAL_DIRECTORY_BYTES,
    max_json_tokens: int = MAX_JSON_TOKENS,
    max_json_depth: int = MAX_JSON_DEPTH,
    rederive: bool = True,
) -> LoadedResult:
    """Read and verify a bundle; raises :class:`ResultFormatError` on any defect (a malformed
    bundle never surfaces as a bare ``TypeError``/``KeyError``). The keyword-only budgets are
    documented in decision 0045 section 1. ``rederive`` (default True) rebuilds ``effective``
    from ``requested`` and compares (amendment A-1); ``False`` records ``effective_binding ==
    "typed-only"`` (needs no data cache; never used by ``ionmc verify``)."""
    try:
        return _load(
            path,
            rederive=rederive,
            max_total_bytes=max_total_bytes,
            max_array_bytes=max_array_bytes,
            max_json_bytes=max_json_bytes,
            max_central_directory_bytes=max_central_directory_bytes,
            max_json_tokens=max_json_tokens,
            max_json_depth=max_json_depth,
        )
    except ResultFormatError:
        raise
    except (
        KeyError, TypeError, ValueError, AttributeError, IndexError, OverflowError, OSError,
        RecursionError, ConfigError,
    ) as exc:  # fmt: skip
        raise _fail(f"malformed bundle ({type(exc).__name__}: {exc})") from exc


def _load(
    path: str | os.PathLike[str],
    *,
    max_total_bytes: int = MAX_TOTAL_BYTES,
    max_array_bytes: int = MAX_ARRAY_BYTES,
    max_json_bytes: int = MAX_JSON_BYTES,
    max_central_directory_bytes: int = MAX_CENTRAL_DIRECTORY_BYTES,
    max_json_tokens: int = MAX_JSON_TOKENS,
    max_json_depth: int = MAX_JSON_DEPTH,
    rederive: bool = True,
) -> LoadedResult:
    root = Path(path)
    if root.is_symlink() or not root.is_dir():
        raise _fail("bundle path is not a directory")
    rj = _bundle_file_json(root)
    if rj.stat().st_size > max_json_bytes:
        raise _fail("result.json exceeds max_json_bytes")
    raw = rj.read_bytes()
    prescan_json(raw, max_json_tokens, max_json_depth)
    try:
        doc = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_dup_hook,
            parse_float=_parse_float,
            parse_constant=_parse_constant,
        )
    except (ValueError, UnicodeDecodeError, RecursionError) as exc:
        if isinstance(exc, ResultFormatError):
            raise
        raise _fail(f"result.json is not valid JSON: {exc}") from exc
    _check_structure(doc)
    entries = doc["arrays"]
    check_allowlist(entries)
    files = doc["files"]
    for name in files:
        if not _file_name_ok(name):
            raise _fail(f"unsafe file name {name!r} in the files map")
    import stat as _stat

    for de in os.scandir(root):
        st = de.stat(follow_symlinks=False)
        if de.name == "result.json":
            continue
        if _stat.S_ISLNK(st.st_mode):
            raise _fail(f"{de.name!r} is a symbolic link (refused)")
        if de.name == "plots" and _stat.S_ISDIR(st.st_mode):
            for pe in os.scandir(root / "plots"):
                pst = pe.stat(follow_symlinks=False)
                if _stat.S_ISLNK(pst.st_mode) or not _stat.S_ISREG(pst.st_mode):
                    raise _fail(f"plots/{pe.name!r} is a symbolic link or not a regular file")
                if f"plots/{pe.name}" not in files:
                    raise _fail(f"unexpected entry 'plots/{pe.name}' in the bundle")
            continue
        if de.name not in files:
            raise _fail(f"unexpected entry {de.name!r} in the bundle")
    for need in _FIXED_FILES:
        if need not in files:
            raise _fail(f"files has no entry for {need}")
    for name, rec in files.items():
        p = _bundle_file(root, name)
        sha, size = _sha256_path(p)
        if size != rec["size"] or sha != rec["sha256"]:
            raise _fail(f"file {name!r} does not match its recorded sha256 and size")
    npz = root / "arrays.npz"
    _zip_preflight(npz, len(entries), max_total_bytes, max_central_directory_bytes)
    by_member = {e["name"] + ".npy": e for e in entries}
    total = 0
    sizes = {n: _entry_nbytes(e, max_array_bytes) for n, e in by_member.items()}
    total = sum(v[2] for v in sizes.values())
    if total > max_total_bytes:
        raise _fail("catalogue arrays exceed max_total_bytes")
    arrays: dict[str, np.ndarray[Any, Any]] = {}
    try:
        zf = zipfile.ZipFile(npz)
    except (zipfile.BadZipFile, OSError) as exc:
        raise _fail(f"arrays.npz is not a readable zip archive: {exc}") from exc
    with zf:
        infos = zf.infolist()
        names = [i.filename for i in infos]
        if len(set(names)) != len(names):
            raise _fail("arrays.npz has duplicate members")
        for n in names:
            if not _MEMBER_RE.fullmatch(n) or n not in by_member:
                raise _fail(f"arrays.npz member {n!r} is not in the catalogue")
        if set(names) != set(by_member):
            raise _fail("arrays.npz lacks catalogued members")
        for info in infos:
            e = by_member[info.filename]
            dt, shape, nbytes = sizes[info.filename]
            a = _read_member(zf, info, e, dt, shape, nbytes)
            if array_content_sha256(a) != e["content_sha256"]:
                raise _fail(f"array {e['name']!r}: content_sha256 mismatch")
            a.flags.writeable = False
            arrays[e["name"]] = a
    binding, limits = _check_invariants(doc, arrays, rederive)
    return LoadedResult(root, doc, arrays, binding, limits)


def _bundle_file_json(root: Path) -> Path:
    import stat as _stat

    p = root / "result.json"
    try:
        st = os.lstat(p)
    except OSError as exc:
        raise _fail("result.json is missing") from exc
    if not _stat.S_ISREG(st.st_mode):
        raise _fail("result.json is not a regular file")
    real_root = os.path.realpath(root)
    if os.path.commonpath([real_root, os.path.realpath(p)]) != real_root:
        raise _fail("result.json resolves outside the bundle")
    return p
