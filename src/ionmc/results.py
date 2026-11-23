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
``input_manifest`` by key forms, ``effective`` only at its top level (deeper levels are a verbatim
copy of ``summary()``); ``provenance.data/models/environment`` are free-form until S3 defines them.

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
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.configio import (
    ConfigError,
    InputManifest,
    array_content_sha256,
    manifest_for_api_config,
    sha256_canonical,
    to_dict,
    validate_dict_structure,
)
from ionmc.geometry import VoxelGeometry
from ionmc.simulation import (
    IonTransportCounters,
    NuclearEnergyBalance,
    NuclearTransportCounters,
    Result,
    TransportCounters,
)

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
    {"MeV", "Gy", "g", "mm", "mm^-2", "keV/um", "g/cm3", "1", "bool", "count"}
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
# Union of the top-level keys of EffectiveConfig.summary() over the REG-2 summaries fixture
# (tests/ionmc/data/v3_006/summaries-95937553.json; pinned equal by a test). Deeper levels of
# ``effective`` are a verbatim copy of summary() and are not checked key by key.
EFFECTIVE_TOP_KEYS = (
    "backend precision production seed n_histories n_batches cpu_workers chunk_histories "
    "max_steps scoring_pieces max_steps_origin unit_direction source physics mcs_mapping "
    "geometry scoring tables stopping scattering_length_g_cm2 rng tallies nuclear"
).split()
_ALLOWED_DTYPES = frozenset({"<f8", "<f4", "<i4", "<i8", "<u4", "|b1", "|u1", "|i1"})
_NAME_RE = re.compile(r"^[a-z0-9_-]{1,64}$")
_MEMBER_RE = re.compile(r"^[a-z0-9_-]{1,64}\.npy$")
_FILE_RE = re.compile(r"^(?:[a-z0-9_-]{1,64}\.[a-z0-9]{1,8}|plots/[a-z0-9_-]{1,64}\.svg)$")
_SHA_RE = re.compile(r"^[0-9a-f]{64}$")
_TOKEN_RE = re.compile(
    rb'"(?:[^"\\]++|\\.)*+"|-?[0-9][0-9.eE+\-]*+|true|false|null|[{}\[\]:]', re.DOTALL
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
            "normalisation": "per primary history",
            "uncertainty": (
                "batch standard error of the mean over B batches (history h in batch h mod B); "
                "ratio quantities by the delta method (decision 0040)"
            ),
            "defined_mask": "positive mean energy and positive voxel mass",
            "sparse_mask": "defined voxels with fewer than B/2 nonzero batches",
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
    """Single pass over raw JSON bytes: count tokens and object members, track nesting depth.
    Allocation does not grow with the content (lazy regex iteration, early exit)."""
    tokens = depth = 0
    for m in _TOKEN_RE.finditer(data):
        c = m.group()[0]
        if c in b"}]":
            depth -= 1
            continue
        tokens += 1
        if c in b"{[":
            depth += 1
            if depth > max_depth:
                raise _fail(f"result.json nesting depth exceeds {max_depth}")
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
    for k in _PROV:  # free-form until S3 defines them
        _need(isinstance(doc["provenance"][k], dict), f"provenance.{k}: expected an object")
    _need(doc["units"] == UNITS_TABLE, "units: differs from the unit table")
    _check_extra(_keys(doc["effective_extra"], _EXTRA, "effective_extra"))
    try:
        validate_dict_structure(doc["requested"])
    except ConfigError as exc:
        raise _fail(f"requested: {exc}") from exc
    unknown = sorted(set(doc["effective"]) - set(EFFECTIVE_TOP_KEYS))
    if unknown:
        raise _fail(f"effective: unknown field(s) {unknown}")
    src = doc["effective"].get("source")
    _need(isinstance(src, dict) and _is_int(src.get("z")), "effective.source.z: expected integer")
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


def _check_manifest_binding(doc: dict[str, Any]) -> None:
    """Bind ``input_manifest`` to the persisted ``requested`` document (no file is opened)."""
    req, m = doc["requested"], doc["input_manifest"]
    if sha256_canonical(req) != m["config_sha256"]:
        raise _fail("input_manifest.config_sha256 does not match the requested document")
    g = req["geometry"]
    recs: list[tuple[str | None, dict[str, Any]]] = []
    if g["type"] == "voxel":
        recs.append(("material_index", g["material_index"]))
        if g.get("density") is not None:
            recs.append(("density", g["density"]))
    recs += [(None, r) for r in req.get("lookups", [])]
    files = m["files"]
    by_role = {f["role"]: f for f in files}
    refs = []
    for want_role, r in recs:
        if "ref" not in r:
            continue
        role = r["ref"]
        ok = (
            role == want_role if want_role else isinstance(role, str) and role.startswith("lookup:")
        )
        _need(ok, f"requested reference {role!r} has the wrong role for its position")
        f = by_role.get(role)
        _need(f is not None, f"requested reference {role!r} has no manifest entry")
        assert f is not None
        _need(
            f["sha256"] == r["sha256"] and f["content_sha256"] == r["content_sha256"],
            f"requested reference {role!r} disagrees with the manifest entry",
        )
        refs.append(role)
    _need(len(set(refs)) == len(refs), "requested references repeat a manifest role")
    for f in files:
        if isinstance(f["sha256"], str):
            _need(f["role"] in refs, f"manifest file {f['role']!r} is not referenced")
    _need(len(recs) == len(files), "requested array positions and manifest files differ in number")
    _need(
        sorted(r["content_sha256"] for _, r in recs) == sorted(f["content_sha256"] for f in files),
        "requested array digests and manifest content digests differ",
    )


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


_LINEAR_UNITS = {
    "edep": "MeV",
    "dose": "Gy",
    "fluence": "mm^-2",
    "fluence_spectrum": "mm^-2",
    "let_t": "keV/um",
    "let_d": "keV/um",
    "let_d_eps": "keV/um",
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
    units: str | None  # None: any lookup unit
    quantity: str | None = None
    index: int | None = None


def _quantity_plan(
    req: Mapping[str, Any], grid: str, shape: list[int]
) -> list[tuple[str, str, str | None, list[int]]]:
    """Ordered (name, quantity, units, shape) of the quantities of ``grid`` the run must have:
    the requested tallies, then the optional excluded-energy channel and the two piece counts."""
    out: list[tuple[str, str, str | None, list[int]]] = []
    for t in req.get("tallies", []):
        if t["grid"] != grid:
            continue
        q = t["quantity"]
        shp = list(shape)
        if q == "fluence_spectrum":
            shp.append(len(t["energy_edges_mev_per_u"]) + 1)
        out.append((t["name"], q, _LINEAR_UNITS.get(q), shp))
    return out


def _expected_catalogue(doc: dict[str, Any], gs: dict[str, str]) -> dict[tuple[Any, ...], _Expect]:
    """Every array the bundle must hold, derived from ``requested``/``effective`` only."""
    req, eff = doc["requested"], doc["effective"]
    nb = doc["run"]["n_batches"]
    persisted = doc["statistics"]["batch_arrays"]["status"] == "persisted"
    exp: dict[tuple[Any, ...], _Expect] = {}
    for g in eff["scoring"]:
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
        plan = _quantity_plan(req, n, shape)
        if plan:
            tail = [
                ("scoring_pieces", "count", "1", shape),
                ("scoring_pieces_local", "count", "1", shape),
            ]
            excl: tuple[str, str, str | None, list[int]] = (
                "edep_excluded_from_let", "edep", "MeV", shape,
            )  # fmt: skip
            actual = _actual_quantity_names(doc, n)
            plan = [*plan, *([excl] if excl[0] in actual else []), *tail]
        used: set[str] = set()
        for qi, (qn, q, qunits, shp) in enumerate(plan):
            qs = safe_name(qn, used)
            axes = list(_SPATIAL) + (["energy_bin"] if len(shp) == 4 else [])
            for base in QUANTITY_ARRAYS:
                qu: str | None = {"defined_mask": "bool", "n_nonzero": "count"}.get(base, qunits)
                dt = _I4 if base == "n_nonzero" else _QUANTITY_DTYPES[base]
                exp[("quantity", n, qn, base)] = _Expect(
                    f"grid__{gs[n]}__q__{qs}__{base}", dt, shp, axes, qu, q, qi
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


def _actual_quantity_names(doc: dict[str, Any], grid: str) -> set[str]:
    return {
        e["origin"].get("name")
        for e in doc["arrays"]
        if e["origin"].get("kind") == "quantity" and e["origin"].get("grid") == grid
    }


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
    if [g["name"] for g in rsc] != gnames or [list(g["shape"]) for g in rsc] != [
        g["shape"] for g in sc
    ]:
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


def _check_invariants(doc: dict[str, Any], arrays: dict[str, np.ndarray[Any, Any]]) -> None:
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


# ---------------------------------------------------------------------------------------------
# load


@dataclass(frozen=True)
class LoadedResult:
    """A verified bundle: the JSON document and every array (read-only)."""

    path: Path
    doc: dict[str, Any]
    arrays: dict[str, np.ndarray[Any, Any]] = field(repr=False)

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
) -> LoadedResult:
    """Read and verify a bundle; raises :class:`ResultFormatError` on any defect (a malformed
    bundle never surfaces as a bare ``TypeError``/``KeyError``). The keyword-only budgets are
    documented in decision 0045 section 1."""
    try:
        return _load(
            path,
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
    _check_invariants(doc, arrays)
    return LoadedResult(root, doc, arrays)


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
