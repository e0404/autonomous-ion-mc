"""Configuration files ``ionmc.config/1``: canonical dict form, strict TOML loader, input manifest.

This module implements decision 0045 sections 6 and 7 (loader part): the canonical JSON form of a
:class:`~ionmc.config.SimulationConfig` (:func:`to_dict`), its strict inverse (:func:`from_dict`),
the TOML loader (:func:`load_toml`), the immutable :class:`InputManifest` that binds a loaded
configuration to the files it was built from, and the dataset checks of the loader
(:func:`verify_datasets`). ``SimulationConfig``, ``Simulation`` and ``Result`` are unchanged: the
manifest is a side carrier that the caller hands to ``to_dict`` and to the result writer.

Canonical form. ``to_dict`` returns plain dicts, lists, strings, numbers, booleans and ``None``;
``json.dumps(..., sort_keys=True, allow_nan=False)`` succeeds on it. It covers every field of
every configuration dataclass (a coverage check at import time fails if a field is added to a
dataclass and not to the schema here). Arrays (voxel ``material_index`` and density, lookup
tables) are never serialised by value:

* with a manifest, an array whose ``item_sha256`` equals that of a manifest file entry is emitted
  as the reference ``{"ref": role, "sha256": file hash, "content_sha256": array digest}``;
* otherwise it is emitted as ``{"status": "not-reconstructible", "reason": ..., "content_sha256":
  ...}``;
* a user ``stopping`` object is emitted as a ``not-reconstructible`` record that also carries
  ``stopping_object`` (class name, model id, spec hash) and ``reconstructible = false``.

``from_dict`` refuses every ``not-reconstructible`` record (``ConfigError`` naming the key path)
and never substitutes a default. A reference is reloaded through the manifest entry of its role:
the file is re-hashed and both the file hash and the content digest are checked.

Digests. ``content_sha256`` of an array is :func:`array_content_sha256`: sha256 over
``dtype.str``, ``|``, the comma-separated shape, ``|`` and the C-order bytes (independent of any
file format or compression). ``item_sha256`` is the sha256 of the canonical JSON of every
serialised field of the item: for a geometry array its role, dtype, shape, content digest and the
geometry fields that describe it (origin, spacing, shape, ``z_exit_mm`` and, for
``material_index``, the full material definitions); for a lookup table its complete provenance
record (name, quantity, units, axis, spacing, range, point count, species, ``file_sha256``,
``content_sha256``, citation, license, source, ``synthetic`` and ``resampling``).
``config_sha256`` of a manifest is the sha256 of the canonical JSON of ``to_dict(config,
manifest=m)``.

File formats referenced by a configuration: voxel arrays are NumPy ``.npy`` files (read with
``allow_pickle=False``), lookup tables are the JSON files of :meth:`ionmc.lookup.LookupTable.
from_file`. Every referenced file carries a mandatory sha256; relative paths resolve against the
directory of the configuration file (``..`` allowed) and are recorded as absolute paths.

No environment variable overrides a physics or run field: the loader reads the file and its
``base_dir`` only (the cache location is the one environment input, read by the stopping and
nuclear loaders exactly as for an API configuration).
"""

from __future__ import annotations

import hashlib
import io
import json
import re
import tomllib
from collections.abc import Callable, Mapping
from dataclasses import dataclass, fields, replace
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.config import (
    DiagnosticsOptions,
    InfluenceOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
)
from ionmc.data import cache
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.lookup import LookupTable
from ionmc.materials import MATERIALS, Material, SternheimerParameters
from ionmc.physics import stopping_models
from ionmc.physics.projectiles import PROJECTILES, Projectile
from ionmc.physics.stopping_models import spec_sha256
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.sources import Beamlet, BeamletSet, PencilBeamSource

SCHEMA = "ionmc.config"
SCHEMA_VERSION = "1.0"
UNAVAILABLE = {"status": "unavailable", "reason": "constructed in process"}


class ConfigError(UnsupportedCombinationError):
    """A configuration file or dict is refused (decision 0045 section 6).

    ``key_path`` names the offending key (for example ``physics.nuclear`` or ``scoring[0].name``)
    and is the first part of the message. It is an ``UnsupportedCombinationError`` (a
    ``ValueError``); the CLI maps it to exit code 2 before its base class.
    """

    def __init__(self, key_path: str, message: str) -> None:
        self.key_path = key_path
        super().__init__(f"{key_path}: {message}")


class DatasetError(Exception):
    """A cache object needed by the selected models is missing or does not match its pin.

    Deliberately not an ``UnsupportedCombinationError`` (decision 0045 section 7); the original
    error is ``__cause__``.
    """


# ---------------------------------------------------------------------------------------------
# digests


def canonical_json(obj: Any) -> str:
    """Canonical JSON text: sorted keys, no whitespace, ASCII, no NaN or infinity."""
    return json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    )


def sha256_canonical(obj: Any) -> str:
    """sha256 (hex) of :func:`canonical_json` of ``obj``."""
    return hashlib.sha256(canonical_json(obj).encode("ascii")).hexdigest()


def array_content_sha256(arr: Any) -> str:
    """Content digest of an array: sha256 over ``dtype.str``, ``|``, the comma-separated shape,
    ``|`` and the C-order bytes (compression and file format independent)."""
    a = np.ascontiguousarray(arr)
    h = hashlib.sha256()
    h.update(a.dtype.str.encode("ascii"))
    h.update(b"|")
    h.update(",".join(str(int(n)) for n in a.shape).encode("ascii"))
    h.update(b"|")
    h.update(a.tobytes(order="C"))
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------
# field schemas (key -> converter, required); coverage against the dataclasses is checked below

Conv = Callable[[Any, str], Any]


def _describe(v: Any) -> str:
    return type(v).__name__


def _bool(v: Any, kp: str) -> bool:
    if not isinstance(v, bool):
        raise ConfigError(kp, f"must be a boolean, got {_describe(v)}")
    return v


def _str(v: Any, kp: str) -> str:
    if not isinstance(v, str):
        raise ConfigError(kp, f"must be a string, got {_describe(v)}")
    return v


def _int(v: Any, kp: str) -> int:
    if isinstance(v, bool) or not isinstance(v, int):
        raise ConfigError(kp, f"must be an integer, got {_describe(v)}")
    return int(v)


def _float(v: Any, kp: str) -> float:
    if isinstance(v, bool) or not isinstance(v, int | float):
        raise ConfigError(kp, f"must be a number, got {_describe(v)}")
    x = float(v)
    if x != x or x in (float("inf"), float("-inf")):
        raise ConfigError(kp, "must be finite")
    return x


def _seq(v: Any, kp: str) -> list[Any]:
    if not isinstance(v, list | tuple):
        raise ConfigError(kp, f"must be an array, got {_describe(v)}")
    return list(v)


def _triple(conv: Conv) -> Conv:
    def f(v: Any, kp: str) -> tuple[Any, Any, Any]:
        s = _seq(v, kp)
        if len(s) != 3:
            raise ConfigError(kp, f"must have exactly 3 entries, got {len(s)}")
        return (conv(s[0], f"{kp}[0]"), conv(s[1], f"{kp}[1]"), conv(s[2], f"{kp}[2]"))

    return f


def _tuple_of(conv: Conv) -> Conv:
    def f(v: Any, kp: str) -> tuple[Any, ...]:
        return tuple(conv(x, f"{kp}[{i}]") for i, x in enumerate(_seq(v, kp)))

    return f


def _opt(conv: Conv) -> Conv:
    def f(v: Any, kp: str) -> Any:
        return None if v is None else conv(v, kp)

    return f


Schema = Mapping[str, tuple[Conv, bool]]

PHYSICS_FIELDS: Schema = {
    "nuclear": (_bool, True),
    "stopping_model": (_str, True),
    "straggling": (_bool, False),
    "straggling_model": (_str, False),
    "multiple_scattering": (_bool, False),
    "mcs_model": (_str, False),
    "delta_electrons": (_str, False),
    "e_cut_mev": (_float, False),
    "max_step_mm": (_float, False),
    "max_energy_loss_fraction": (_float, False),
    "range_alpha": (_float, False),
    "range_rho_f_mm": (_float, False),
    "short_step_fraction": (_float, False),
    "truncated_hinge_diagnostic": (_bool, False),
    "e_cut_deuteron_mev": (_float, False),
    "nuclear_table_id": (_opt(_str), False),
    "e_cut_ion_mev_per_u": (_float, False),
    "mcs_mapping_eval": (_str, False),
    "end_of_range_model": (_str, False),
}
RUN_FIELDS: Schema = {
    "backend": (_str, True),
    "precision": (_str, True),
    "seed": (_int, True),
    "n_histories": (_int, True),
    "n_batches": (_int, False),
    "cpu_workers": (_int, False),
    "max_steps": (_opt(_int), False),
    "allow_invalid_result": (_bool, False),
    "worker_timeout_s": (_opt(_float), False),
    "memory_budget_bytes": (_int, False),
    "chunk_histories": (_int, False),
}
DIAGNOSTICS_FIELDS: Schema = {
    "track_end_positions": (_bool, False),
    "escape_records": (_bool, False),
    "trace_histories": (_int, False),
}
SOURCE_FIELDS: Schema = {
    "projectile": (lambda v, kp: v, True),  # a registered name or a full definition
    "position_mm": (_triple(_float), True),
    "direction": (_triple(_float), True),
    "kinetic_energy_mev": (_float, True),
    "energy_sigma_mev": (_float, False),
    "lateral_sigma_mm": (_float, False),
}
BEAMLET_FIELDS: Schema = {
    "position_mm": (_triple(_float), True),
    "direction": (_triple(_float), True),
    "kinetic_energy_mev": (_float, True),
    "energy_sigma_mev": (_float, False),
    "lateral_sigma_mm": (_float, False),
    "n_histories": (_int, True),
    "weight": (_opt(_float), False),
    "label": (_opt(_str), False),
}
INFLUENCE_FIELDS: Schema = {
    "grid": (_str, True),
    "threshold_rel": (_float, False),
    "group_size": (_opt(_int), False),
    "persist_batches": (_bool, False),
}
SOURCE_TYPES = ("pencil", "beamlets")
PROJECTILE_FIELDS: Schema = {
    "name": (_str, True),
    "symbol": (_str, True),
    "z": (_int, True),
    "a": (_int, True),
    "mass_mev": (_float, True),
}
SCORING_FIELDS: Schema = {
    "origin_mm": (_triple(_float), True),
    "spacing_mm": (_triple(_float), True),
    "shape": (_triple(_int), True),
    "name": (_str, False),
}
TALLY_FIELDS: Schema = {
    "name": (_str, True),
    "grid": (_str, True),
    "quantity": (_str, True),
    "species": (_opt(_tuple_of(_str)), False),
    "generation": (_str, False),
    "lookup": (_opt(_str), False),
    "energy_edges_mev_per_u": (_opt(_tuple_of(_float)), False),
    "let_medium": (_str, False),
    "dose_reference": (_str, False),
}
MATERIAL_FIELDS: Schema = {
    "name": (_str, True),
    "density_g_cm3": (_float, True),
    "mass_fractions": (lambda v, kp: v, True),
    "I_eV": (_opt(_float), False),
    "sternheimer": (lambda v, kp: v, False),
    "source": (_str, False),
}
STERNHEIMER_FIELDS: Schema = {k: (_float, True) for k in ("x0", "x1", "cbar", "a", "m")}
BOX_FIELDS: Schema = {
    "type": (_str, True),
    "lower_mm": (_triple(_float), True),
    "size_mm": (_triple(_float), True),
    "material": (_str, True),
}
VOXEL_FIELDS: Schema = {
    "type": (_str, True),
    "origin_mm": (_triple(_float), True),
    "spacing_mm": (_triple(_float), True),
    "shape": (_triple(_int), True),
    "materials": (_tuple_of(_str), True),
    "z_exit_mm": (_opt(_float), False),
    "material_index": (lambda v, kp: v, True),  # array item
    "density": (lambda v, kp: v, False),  # array item or null
}
TOP_KEYS = (
    "schema",
    "schema_version",
    "source",
    "materials",
    "geometry",
    "scoring",
    "physics",
    "run",
    "diagnostics",
    "tallies",
    "lookups",
    "influence",
    "output",
)
REQUIRED_TOP = ("source", "geometry", "scoring", "physics", "run")

# schema key -> dataclass field when the names differ (decision 0045 section 6 names the density
# file key ``density``)
_RENAMED = {("VoxelGeometry", "density"): "density_g_cm3"}


def _check_schema_coverage() -> None:
    """Fail at import if a configuration dataclass has a field the schema does not know."""
    no_type = {"type"}
    pairs: list[tuple[type, set[str], set[str]]] = [
        (PhysicsOptions, set(PHYSICS_FIELDS), {"stopping"}),
        (RunOptions, set(RUN_FIELDS), set()),
        (DiagnosticsOptions, set(DIAGNOSTICS_FIELDS), set()),
        (PencilBeamSource, set(SOURCE_FIELDS), set()),
        (Beamlet, set(BEAMLET_FIELDS), set()),
        (InfluenceOptions, set(INFLUENCE_FIELDS), set()),
        (BeamletSet, {"projectile", "beamlets"}, set()),
        (Projectile, set(PROJECTILE_FIELDS), set()),
        (ScoringGrid, set(SCORING_FIELDS), set()),
        (TallyRequest, set(TALLY_FIELDS), set()),
        (Material, set(MATERIAL_FIELDS), set()),
        (SternheimerParameters, set(STERNHEIMER_FIELDS), set()),
        (BoxPhantom, set(BOX_FIELDS) - no_type, set()),
        (
            VoxelGeometry,
            {_RENAMED.get(("VoxelGeometry", k), k) for k in VOXEL_FIELDS} - no_type,
            set(),
        ),
    ]
    for cls, keys, extra in pairs:
        have = {f.name for f in fields(cls) if f.init}
        known = keys | extra
        if have != known:
            raise AssertionError(
                f"configio schema out of date for {cls.__name__}: "
                f"unknown to the schema {sorted(have - known)}, "
                f"unknown to the dataclass {sorted(known - have)}"
            )
    top = {f.name for f in fields(SimulationConfig)}
    covered = {
        "source",
        "geometry",
        "scoring",
        "physics",
        "run",
        "diagnostics",
        "tallies",
        "lookups",
        "influence",
    }
    if top != covered:
        raise AssertionError(
            f"configio schema out of date for SimulationConfig: {sorted(top ^ covered)}"
        )


_check_schema_coverage()


# ---------------------------------------------------------------------------------------------
# serialisation: config -> canonical dict


def _ser(v: Any) -> Any:
    """Plain JSON value of a field value (numpy scalars to Python, tuples to lists)."""
    if v is None or isinstance(v, bool | str):
        return v
    if isinstance(v, int | np.integer):
        return int(v)
    if isinstance(v, float | np.floating):
        return float(v)
    if isinstance(v, tuple | list | np.ndarray):
        return [_ser(x) for x in v]
    raise TypeError(f"cannot serialise {type(v).__name__}")


def _fields_dict(obj: Any, schema: Schema) -> dict[str, Any]:
    cls = type(obj).__name__
    return {k: _ser(getattr(obj, _RENAMED.get((cls, k), k))) for k in schema}


def _material_dict(m: Material) -> dict[str, Any]:
    s = m.sternheimer
    return {
        "name": m.name,
        "density_g_cm3": float(m.density_g_cm3),
        "mass_fractions": {k: float(v) for k, v in m.mass_fractions.items()},
        "I_eV": None if m.I_eV is None else float(m.I_eV),
        "sternheimer": None if s is None else _fields_dict(s, STERNHEIMER_FIELDS),
        "source": m.source,
    }


def _material_entry(m: Material) -> dict[str, Any]:
    """``{"name": n}`` for a built-in material, else the full definition."""
    return {"name": m.name} if MATERIALS.get(m.name) == m else _material_dict(m)


def _projectile_entry(p: Projectile) -> Any:
    return p.name if PROJECTILES.get(p.name) == p else _fields_dict(p, PROJECTILE_FIELDS)


def _geometry_materials(geometry: VoxelGeometry | BoxPhantom) -> tuple[Material, ...]:
    return (geometry.material,) if isinstance(geometry, BoxPhantom) else geometry.materials


def _collect_materials(config: SimulationConfig) -> list[Material]:
    seen: dict[str, Material] = {}
    for m in _geometry_materials(config.geometry):
        if m.name in seen and seen[m.name] != m:
            raise ConfigError(
                "geometry.materials",
                f"two different materials are named {m.name!r}; ionmc.config/1 refers to "
                "materials by name and cannot represent this configuration",
            )
        seen.setdefault(m.name, m)
    return list(seen.values())


def _geometry_items(geometry: VoxelGeometry) -> dict[str, dict[str, Any]]:
    """Item payloads (role -> JSON object with every serialised field) of the geometry arrays."""
    desc: dict[str, Any] = {
        "origin_mm": _ser(geometry.origin_mm),
        "spacing_mm": _ser(geometry.spacing_mm),
        "shape": _ser(geometry.shape),
        "z_exit_mm": _ser(geometry.z_exit_mm),
    }
    out: dict[str, dict[str, Any]] = {}
    for role, arr in (
        ("material_index", geometry.material_index),
        ("density", geometry.density_g_cm3),
    ):
        if arr is None:
            continue
        item: dict[str, Any] = {
            "role": role,
            "dtype": arr.dtype.str,
            "array_shape": [int(n) for n in arr.shape],
            "content_sha256": array_content_sha256(arr),
            **desc,
        }
        if role == "material_index":
            item["materials"] = [_material_dict(m) for m in geometry.materials]
        out[role] = item
    return out


def reader_item_sha256s(
    requested: Mapping[str, Any],
    geometry_arrays: Mapping[str, Any],
    lookup_provenance: list[dict[str, Any]],
    material_resolver: Callable[[str], dict[str, Any]] | None = None,
) -> dict[str, str]:
    """Role -> ``item_sha256`` recomputed from a persisted ``requested`` document, the persisted
    geometry arrays and the persisted lookup provenance records (no file is opened); equals the
    loader's value for the same configuration (see :func:`_geometry_items`). A name-only voxel
    material is expanded by ``material_resolver(name)`` (the ``_material_dict`` form), by default
    through the live ``MATERIALS`` registry; a result reader passes its frozen definitions."""
    out: dict[str, str] = {}
    g = requested["geometry"]
    if g["type"] == "voxel":
        defs = {m["name"]: m for m in requested.get("materials", [])}
        mats = []
        for name in g["materials"]:
            rec = defs[name]
            if set(rec) != {"name"}:
                mats.append(rec)
            elif material_resolver is not None:
                mats.append(material_resolver(name))
            else:
                mats.append(_material_dict(MATERIALS[name]))
        for role, arr in geometry_arrays.items():
            item: dict[str, Any] = {
                "role": role,
                "dtype": arr.dtype.str,
                "array_shape": [int(n) for n in arr.shape],
                "content_sha256": array_content_sha256(arr),
                "origin_mm": g["origin_mm"],
                "spacing_mm": g["spacing_mm"],
                "shape": g["shape"],
                "z_exit_mm": g["z_exit_mm"],
            }
            if role == "material_index":
                item["materials"] = mats
            out[role] = sha256_canonical(item)
    for prov in lookup_provenance:
        out[f"lookup:{prov['name']}"] = sha256_canonical(prov)
    return out


def _config_items(config: SimulationConfig) -> dict[str, tuple[str, str]]:
    """Role -> (content_sha256, item_sha256) of every array of the configuration."""
    out: dict[str, tuple[str, str]] = {}
    if isinstance(config.geometry, VoxelGeometry):
        for role, item in _geometry_items(config.geometry).items():
            out[role] = (item["content_sha256"], sha256_canonical(item))
    for t in config.lookups:
        out[f"lookup:{t.name}"] = (t.content_sha256, sha256_canonical(t.provenance()))
    return out


def _find_ref(
    role_hint: str | None, item_sha: str, manifest: InputManifest | None
) -> InputFile | None:
    if manifest is None:
        return None
    for f in manifest.files:
        if (
            f.sha256 is not None
            and f.item_sha256 == item_sha
            and (role_hint is None or f.role == role_hint)
        ):
            return f
    return None


def _array_record(
    content_sha: str, item_sha: str, manifest: InputManifest | None, role_hint: str | None
) -> dict[str, Any]:
    entry = _find_ref(role_hint, item_sha, manifest)
    if entry is not None:
        return {"ref": entry.role, "sha256": entry.sha256, "content_sha256": content_sha}
    return {
        "status": "not-reconstructible",
        "reason": "array exists only in memory (no input-file entry of the manifest matches it)",
        "content_sha256": content_sha,
    }


def _physics_dict(ph: PhysicsOptions) -> dict[str, Any]:
    out = _fields_dict(ph, PHYSICS_FIELDS)
    if ph.stopping is not None:
        sel = ph.selection
        out["stopping"] = {
            "status": "not-reconstructible",
            "reason": "a user stopping object cannot be rebuilt from a configuration",
            "reconstructible": False,
            "stopping_object": {
                "class": sel.source_class,
                "model_id": sel.requested_id,
                "spec_sha256": spec_sha256(sel.spec),
            },
        }
    return out


def to_dict(config: SimulationConfig, *, manifest: InputManifest | None = None) -> dict[str, Any]:
    """Canonical JSON-serialisable form of ``config`` (module docstring).

    With ``manifest``, every array or lookup table whose ``item_sha256`` equals that of a manifest
    file entry is emitted as the reference ``{"ref", "sha256", "content_sha256"}``; every other
    array is emitted as a ``not-reconstructible`` record. The result is deterministic: the same
    config and manifest give the same dict.
    """
    geom = config.geometry
    if isinstance(geom, BoxPhantom):
        geometry: dict[str, Any] = {
            "type": "box",
            "lower_mm": _ser(geom.lower_mm),
            "size_mm": _ser(geom.size_mm),
            "material": geom.material.name,
        }
    else:
        items = _geometry_items(geom)
        geometry = {
            "type": "voxel",
            "origin_mm": _ser(geom.origin_mm),
            "spacing_mm": _ser(geom.spacing_mm),
            "shape": _ser(geom.shape),
            "materials": [m.name for m in geom.materials],
            "z_exit_mm": _ser(geom.z_exit_mm),
        }
        for role in ("material_index", "density"):
            if role in items:
                it = items[role]
                geometry[role] = _array_record(
                    it["content_sha256"], sha256_canonical(it), manifest, role
                )
            else:
                geometry[role] = None
    s = config.source
    if isinstance(s, BeamletSet):
        source: dict[str, Any] = {
            "type": "beamlets",
            "projectile": _projectile_entry(s.projectile),
            "beamlets": [_fields_dict(b, BEAMLET_FIELDS) for b in s.beamlets],
        }
    else:
        source = {k: _ser(getattr(s, k)) for k in SOURCE_FIELDS if k != "projectile"} | {
            "projectile": _projectile_entry(s.projectile)
        }
    out = {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "source": source,
        "materials": [_material_entry(m) for m in _collect_materials(config)],
        "geometry": geometry,
        "scoring": [_fields_dict(g, SCORING_FIELDS) for g in config.scoring],
        "physics": _physics_dict(config.physics),
        "run": _fields_dict(config.run, RUN_FIELDS),
        "diagnostics": _fields_dict(config.diagnostics, DIAGNOSTICS_FIELDS),
        "tallies": [_fields_dict(t, TALLY_FIELDS) for t in config.tallies],
        "lookups": [
            {
                **_array_record(t.content_sha256, sha256_canonical(t.provenance()), manifest, None),
                "units": t.units,  # declared units: the reader derives the channel plan from them
            }
            for t in config.lookups
        ],
    }
    if config.influence is not None:  # absent for pencil configurations (byte-identical form)
        out["influence"] = _fields_dict(config.influence, INFLUENCE_FIELDS)
    return out


def config_sha256(config: SimulationConfig, manifest: InputManifest | None) -> str:
    """sha256 of the canonical JSON of ``to_dict(config, manifest=manifest)``."""
    return sha256_canonical(to_dict(config, manifest=manifest))


# ---------------------------------------------------------------------------------------------
# input manifest


def _sha256_file(path: str, kp: str) -> str:
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except (OSError, ValueError) as exc:
        raise ConfigError(kp, f"cannot read input file {Path(path).name!r}: {exc}") from exc


@dataclass(frozen=True)
class InputFile:
    """One referenced input file. ``path`` and ``sha256`` are None for an array of an API
    configuration (no file; serialised as an ``unavailable`` record)."""

    role: str
    path: str | None
    sha256: str | None
    content_sha256: str
    item_sha256: str
    bytes: int | None = None  # file size at load time (None for an API array)

    def to_dict(self) -> dict[str, Any]:
        return {
            "role": self.role,
            "path": dict(UNAVAILABLE) if self.path is None else self.path,
            "sha256": dict(UNAVAILABLE) if self.sha256 is None else self.sha256,
            "content_sha256": self.content_sha256,
            "item_sha256": self.item_sha256,
            "bytes": dict(UNAVAILABLE) if self.bytes is None else self.bytes,
        }


@dataclass(frozen=True)
class ConfigFileRecord:
    """The configuration file a manifest was loaded from (absolute resolved path, sha256)."""

    path: str
    sha256: str


@dataclass(frozen=True)
class InputManifest:
    """Immutable record of where a loaded configuration came from (decision 0045 section 6).

    ``source`` is ``"toml"`` or ``"api"``. ``base_dir`` and ``config_file`` are None for ``api``.
    ``files`` has one entry per array (role ``material_index``, ``density`` or
    ``lookup:<name>``). ``config_sha256`` binds the manifest to the complete configuration.
    """

    source: str
    base_dir: str | None
    config_file: ConfigFileRecord | None
    files: tuple[InputFile, ...]
    config_sha256: str

    def __post_init__(self) -> None:
        if self.source not in ("toml", "api"):
            raise ConfigError(
                "input_manifest.source", f"must be 'toml' or 'api', got {self.source!r}"
            )
        object.__setattr__(self, "files", tuple(self.files))
        roles = [f.role for f in self.files]
        if len(set(roles)) != len(roles):
            raise ConfigError("input_manifest.files", f"duplicate roles in {sorted(roles)}")
        if self.source == "toml" and (self.base_dir is None or self.config_file is None):
            raise ConfigError("input_manifest", "a toml manifest needs base_dir and config_file")

    def file(self, role: str) -> InputFile | None:
        """The entry of ``role`` or None."""
        for f in self.files:
            if f.role == role:
                return f
        return None

    def to_dict(self) -> dict[str, Any]:
        """JSON form (pointers ``/input_manifest/*`` of the result bundle)."""
        return {
            "source": self.source,
            "base_dir": dict(UNAVAILABLE) if self.base_dir is None else self.base_dir,
            "config_file": (
                dict(UNAVAILABLE)
                if self.config_file is None
                else {"path": self.config_file.path, "sha256": self.config_file.sha256}
            ),
            "files": [f.to_dict() for f in self.files],
            "config_sha256": self.config_sha256,
        }

    def verify(self, config: SimulationConfig) -> None:
        """Raise ``ConfigError`` unless this manifest still describes ``config``.

        Re-hashes every listed file, checks every ``content_sha256`` against the arrays of
        ``config`` and recomputes ``config_sha256`` from ``config`` with this manifest: a manifest
        of another configuration or a configuration changed after loading is refused. The
        configuration file is re-hashed too (``input_manifest.config_file``).
        """
        if self.config_file is not None:
            now = _sha256_file(self.config_file.path, "input_manifest.config_file")
            if now != self.config_file.sha256:
                raise ConfigError(
                    "input_manifest.config_file",
                    f"configuration file {Path(self.config_file.path).name!r} changed since it "
                    f"was loaded (sha256 {now}, recorded {self.config_file.sha256})",
                )
        items = _config_items(config)
        for i, f in enumerate(self.files):
            kp = f"input_manifest.files[{i}]({f.role})"
            if f.path is not None:
                try:
                    size = Path(f.path).stat().st_size if f.bytes is not None else None
                except OSError as exc:
                    raise ConfigError(
                        kp, f"cannot read input file {Path(f.path).name!r}: {exc}"
                    ) from exc
                if size is not None and size != f.bytes:
                    raise ConfigError(
                        kp, f"input file {Path(f.path).name!r} changed size since it was loaded"
                    )
                now = _sha256_file(f.path, kp)
                if now != f.sha256:
                    raise ConfigError(
                        kp,
                        f"input file {Path(f.path).name!r} changed since it was loaded "
                        f"(sha256 {now}, recorded {f.sha256})",
                    )
            cur = items.get(f.role)
            if cur is None:
                raise ConfigError(kp, "the configuration has no array of this role")
            if cur[0] != f.content_sha256:
                raise ConfigError(
                    kp,
                    f"content digest {cur[0]} of the configuration differs from the "
                    f"manifest's {f.content_sha256}",
                )
            if cur[1] != f.item_sha256:
                raise ConfigError(
                    kp,
                    "an item field of the configuration (provenance, geometry description) "
                    "differs from the manifest entry's item digest",
                )
            if f.path is not None and not f.role.startswith("lookup:"):
                self._check_geometry_file(f, config, kp)
            if f.path is not None and f.role.startswith("lookup:"):
                try:
                    reloaded = LookupTable.from_file(f.path, expected_sha256=f.sha256)
                except (
                    UnsupportedCombinationError,
                    OSError,
                    ValueError,
                    TypeError,
                    KeyError,
                ) as exc:
                    raise ConfigError(kp, f"cannot reload the lookup file: {exc}") from exc
                if sha256_canonical(reloaded.provenance()) != f.item_sha256:
                    raise ConfigError(
                        kp,
                        f"the lookup file {Path(f.path).name!r} no longer yields the recorded "
                        "item (a provenance field of the file changed)",
                    )
        extra = sorted(set(items) - {f.role for f in self.files})
        if extra:
            raise ConfigError(
                "input_manifest.files",
                f"the manifest lacks an entry for the array role(s) {extra} of the configuration",
            )
        if config_sha256(config, self) != self.config_sha256:
            raise ConfigError(
                "input_manifest.config_sha256",
                "the manifest does not describe this configuration (the configuration differs "
                "from the one that was loaded)",
            )

    @staticmethod
    def _check_geometry_file(f: InputFile, config: SimulationConfig, kp: str) -> None:
        """Reload the geometry array file of ``f`` and bind its content and item digests."""
        assert f.path is not None
        try:
            arr = np.load(io.BytesIO(Path(f.path).read_bytes()), allow_pickle=False)
        except (OSError, ValueError, EOFError) as exc:
            raise ConfigError(kp, f"cannot reload the array file: {exc}") from exc
        geometry = config.geometry
        if not isinstance(arr, np.ndarray) or not isinstance(geometry, VoxelGeometry):
            raise ConfigError(kp, "the file does not hold the array of a voxel geometry")
        dtype = np.int32 if f.role == "material_index" else np.float64
        try:
            arr = np.ascontiguousarray(arr, dtype=dtype)
        except (ValueError, TypeError) as exc:
            raise ConfigError(
                kp, f"the array file cannot be read as {np.dtype(dtype)}: {exc}"
            ) from exc
        item = dict(_geometry_items(geometry)[f.role])
        item["content_sha256"] = array_content_sha256(arr)
        item["dtype"] = arr.dtype.str
        item["array_shape"] = [int(n) for n in arr.shape]
        if item["content_sha256"] != f.content_sha256 or sha256_canonical(item) != f.item_sha256:
            raise ConfigError(
                kp,
                f"the array file {Path(f.path).name!r} no longer yields the recorded content "
                "and item digests",
            )


def _finish_manifest(
    source: str,
    base_dir: str | None,
    config_file: ConfigFileRecord | None,
    files: tuple[InputFile, ...],
    config: SimulationConfig,
) -> InputManifest:
    m0 = InputManifest(source, base_dir, config_file, files, "")
    return replace(m0, config_sha256=config_sha256(config, m0))


def manifest_for_api_config(config: SimulationConfig) -> InputManifest:
    """The ``api`` manifest of a configuration built in process: no files, every array listed
    with ``unavailable`` path and sha256 (so no reference is ever emitted for it)."""
    files = tuple(
        InputFile(role, None, None, content, item)
        for role, (content, item) in _config_items(config).items()
    )
    return _finish_manifest("api", None, None, files, config)


# ---------------------------------------------------------------------------------------------
# reconstruction: dict -> config


OUTPUT_FIELDS: Schema = {"path": (_opt(_str), False), "persist_batches": (_bool, False)}


@dataclass(frozen=True)
class OutputOptions:
    """The ``[output]`` table: ``path`` (absolute if given) and ``persist_batches``."""

    path: str | None = None
    persist_batches: bool = True


@dataclass(frozen=True)
class LoadedConfig:
    """Result of :func:`load_toml_full`. ``n_histories_derived`` is True when a beamlet
    configuration omitted ``[run] n_histories`` (derived: sum of the beamlet n_histories)."""

    config: SimulationConfig
    manifest: InputManifest
    output: OutputOptions
    n_histories_derived: bool = False

    @property
    def n_histories_origin(self) -> str:
        """Provenance of ``run.n_histories`` handed to ``Simulation`` (reconciliation evidence)."""
        return "sum-of-beamlets" if self.n_histories_derived else "explicit"


@dataclass(frozen=True)
class _Source:
    """A referenced file read and hash-checked for one array position."""

    path: str
    sha256: str
    raw: bytes
    ref: dict[str, Any] | None
    entry: InputFile | None


def _table(v: Any, kp: str) -> dict[str, Any]:
    if not isinstance(v, dict) or not all(isinstance(k, str) for k in v):
        raise ConfigError(kp, f"must be a table, got {_describe(v)}")
    return v


def _read(tbl: Mapping[str, Any], kp: str, schema: Schema) -> dict[str, Any]:
    """Convert the keys of ``tbl`` by ``schema`` (strict: unknown and missing keys are errors)."""
    for k in tbl:
        if k not in schema:
            raise ConfigError(f"{kp}.{k}", f"unknown key (valid keys: {', '.join(sorted(schema))})")
    out: dict[str, Any] = {}
    for name, (conv, required) in schema.items():
        if name in tbl:
            out[name] = conv(tbl[name], f"{kp}.{name}")
        elif required:
            raise ConfigError(f"{kp}.{name}", "required key is missing")
    return out


def _table_list(v: Any, kp: str) -> list[dict[str, Any]]:
    return [_table(x, f"{kp}[{i}]") for i, x in enumerate(_seq(v, kp))]


def _construct(
    cls: Callable[..., Any], kp: str, kwargs: dict[str, Any], *, dataset: bool = False
) -> Any:
    try:
        return cls(**kwargs)
    except ConfigError:
        raise
    except (UnsupportedCombinationError, ValueError, TypeError) as exc:
        # a chained error of the stopping selection comes from the dataset layer (missing or
        # unparsable cache object): not a key-path defect, keep its class for the exit-code table
        if dataset and isinstance(exc, UnsupportedCombinationError) and exc.__cause__ is not None:
            raise
        raise ConfigError(kp, str(exc)) from exc


def _material_from(v: Any, kp: str) -> Material:
    t = _table(v, kp)
    if set(t) == {"name"}:
        name = _str(t["name"], f"{kp}.name")
        m = MATERIALS.get(name)
        if m is None:
            raise ConfigError(
                f"{kp}.name",
                f"{name!r} is not a built-in material ({', '.join(sorted(MATERIALS))}) and has "
                "no inline definition (density_g_cm3, mass_fractions)",
            )
        return m
    kw = _read(t, kp, MATERIAL_FIELDS)
    fr = _table(kw["mass_fractions"], f"{kp}.mass_fractions")
    kw["mass_fractions"] = {k: _float(x, f"{kp}.mass_fractions.{k}") for k, x in fr.items()}
    if kw.get("sternheimer") is not None:
        st = _read(
            _table(kw["sternheimer"], f"{kp}.sternheimer"), f"{kp}.sternheimer", STERNHEIMER_FIELDS
        )
        kw["sternheimer"] = _construct(SternheimerParameters, f"{kp}.sternheimer", st)
    return _construct(Material, kp, kw)  # type: ignore[no-any-return]


def material_from_record(v: Any) -> Material:
    """The :class:`~ionmc.materials.Material` of a ``requested.materials`` record (``{"name"}``
    for a built-in material, else the full definition); ``ConfigError`` when invalid."""
    return _material_from(v, "materials")


def _projectile_from(v: Any, kp: str) -> Projectile:
    if isinstance(v, str):
        p = PROJECTILES.get(v)
        if p is None:
            raise ConfigError(
                kp, f"unknown projectile {v!r} (registered: {', '.join(sorted(PROJECTILES))})"
            )
        return p
    kw = _read(_table(v, kp), kp, PROJECTILE_FIELDS)
    return _construct(Projectile, kp, kw)  # type: ignore[no-any-return]


def _normalise_sha(v: Any, kp: str) -> str:
    s = _str(v, kp).lower()
    if len(s) != 64 or any(c not in "0123456789abcdef" for c in s):
        raise ConfigError(kp, "must be a 64-digit hexadecimal sha256")
    return s


def _read_source(
    item: Any, kp: str, manifest: InputManifest | None, base_dir: Path | None, ref_role: str | None
) -> _Source:
    """Read and hash-check the file behind an array item (``file`` form or ``ref`` form)."""
    t = _table(item, kp)
    if "status" in t:
        raise ConfigError(
            kp,
            f"{t.get('status')}: {t.get('reason', 'no reason recorded')}; the array cannot be "
            "reconstructed from this dict",
        )
    if "ref" in t:
        ref = _read(
            t, kp, {"ref": (_str, True), "sha256": (_str, True), "content_sha256": (_str, True)}
        )
        _normalise_sha(ref["sha256"], f"{kp}.sha256")
        _normalise_sha(ref["content_sha256"], f"{kp}.content_sha256")
        if ref_role is not None and ref["ref"] != ref_role:
            raise ConfigError(f"{kp}.ref", f"must be {ref_role!r}, got {ref['ref']!r}")
        if manifest is None:
            raise ConfigError(kp, f"references input file {ref['ref']!r} but no manifest was given")
        entry = manifest.file(ref["ref"])
        if entry is None:
            raise ConfigError(kp, f"the manifest has no file entry for role {ref['ref']!r}")
        if entry.path is None or entry.sha256 is None:
            raise ConfigError(
                kp, f"the manifest records no file for role {ref['ref']!r} (unavailable)"
            )
        if entry.sha256 != ref["sha256"]:
            raise ConfigError(f"{kp}.sha256", "differs from the sha256 of the manifest entry")
        path = entry.path
        expected = entry.sha256
        info: tuple[dict[str, Any] | None, InputFile | None] = (ref, entry)
    elif "file" in t:
        d = _read(t, kp, {"file": (_str, True), "sha256": (_str, True)})
        expected = _normalise_sha(d["sha256"], f"{kp}.sha256")
        if not d["file"]:
            raise ConfigError(f"{kp}.file", "must not be empty")
        p = Path(d["file"])
        if not p.is_absolute():
            if base_dir is None:
                raise ConfigError(f"{kp}.file", "a relative path needs a base_dir")
            p = base_dir / p
        try:
            path = str(p.resolve())
        except (OSError, ValueError, RuntimeError) as exc:
            raise ConfigError(f"{kp}.file", f"cannot resolve the path: {exc}") from exc
        info = (None, None)
    else:
        bad = next(iter(t), None)
        if bad is not None:
            raise ConfigError(f"{kp}.{bad}", "unknown key (an array is given as {file, sha256})")
        raise ConfigError(kp, "an array needs {file, sha256}")
    try:
        raw = Path(path).read_bytes()
    except (OSError, ValueError) as exc:
        raise ConfigError(f"{kp}.file", f"cannot read {Path(path).name!r}: {exc}") from exc
    got = hashlib.sha256(raw).hexdigest()
    if got != expected:
        raise ConfigError(
            f"{kp}.sha256", f"file {Path(path).name!r} has sha256 {got}, expected {expected}"
        )
    return _Source(path, got, raw, info[0], info[1])


def _check_ref_content(src: _Source, content: str, item_sha: str, kp: str) -> None:
    """Bind a reloaded item to its manifest entry: content digest and complete item digest."""
    if src.ref is None or src.entry is None:
        return
    if content != src.ref["content_sha256"] or content != src.entry.content_sha256:
        raise ConfigError(
            f"{kp}.content_sha256",
            f"the array of file {Path(src.path).name!r} has content digest {content}, "
            "which differs from the recorded one",
        )
    if item_sha != src.entry.item_sha256:
        raise ConfigError(
            kp,
            f"the item reloaded from {Path(src.path).name!r} differs in a serialised field "
            f"(item digest {item_sha}) from the manifest entry's {src.entry.item_sha256}",
        )


def _load_npy(src: _Source, kp: str) -> np.ndarray[Any, Any]:
    try:
        arr = np.load(io.BytesIO(src.raw), allow_pickle=False)
    except (ValueError, OSError, EOFError) as exc:
        raise ConfigError(
            f"{kp}.file", f"{Path(src.path).name!r} is not a readable .npy array: {exc}"
        ) from exc
    if not isinstance(arr, np.ndarray):
        raise ConfigError(f"{kp}.file", f"{Path(src.path).name!r} must hold a single .npy array")
    return arr


def _provided_array(
    item: Any, kp: str, role: str, arrays: Mapping[str, Any]
) -> np.ndarray[Any, Any]:
    """The array of position ``role`` from the in-memory provider; its content digest must equal
    the one the record (``ref`` or ``not-reconstructible`` form) states."""
    t = _table(item, kp)
    content = _normalise_sha(t.get("content_sha256"), f"{kp}.content_sha256")
    arr = arrays.get(role)
    if not isinstance(arr, np.ndarray):
        raise ConfigError(kp, f"the array provider has no array for role {role!r}")
    if array_content_sha256(arr) != content:
        raise ConfigError(
            f"{kp}.content_sha256", f"the provided array has a different content digest ({role})"
        )
    return np.array(arr)


@dataclass
class _Built:
    """Collected while reading: the files behind arrays and the output table."""

    files: list[tuple[str, str, str]]
    output: OutputOptions
    n_histories_derived: bool = False


def _geometry_from(
    v: Any,
    mats: Mapping[str, Material],
    manifest: InputManifest | None,
    base_dir: Path | None,
    built: _Built,
    arrays: Mapping[str, Any] | None = None,
) -> VoxelGeometry | BoxPhantom:
    kp = "geometry"
    t = _table(v, kp)
    kind = t.get("type")
    if kind == "box":
        kw = _read(t, kp, BOX_FIELDS)
        kw.pop("type")
        kw["material"] = _resolve_material(kw["material"], mats, f"{kp}.material")
        return _construct(BoxPhantom, kp, kw)  # type: ignore[no-any-return]
    if kind != "voxel":
        raise ConfigError(f"{kp}.type", f"must be 'box' or 'voxel', got {kind!r}")
    kw = _read(t, kp, VOXEL_FIELDS)
    kw.pop("type")
    kw["materials"] = tuple(
        _resolve_material(n, mats, f"{kp}.materials[{i}]") for i, n in enumerate(kw["materials"])
    )
    if arrays is not None:  # arrays from memory (a result bundle): no file is opened
        kw["material_index"] = _provided_array(
            kw.pop("material_index"), f"{kp}.material_index", "material_index", arrays
        )
        dens = kw.pop("density", None)
        if dens is not None:
            kw["density_g_cm3"] = _provided_array(dens, f"{kp}.density", "density", arrays)
        return _construct(VoxelGeometry, kp, kw)  # type: ignore[no-any-return]
    src_mi = _read_source(
        kw.pop("material_index"), f"{kp}.material_index", manifest, base_dir, "material_index"
    )
    kw["material_index"] = _load_npy(src_mi, f"{kp}.material_index")
    dens = kw.pop("density", None)
    src_d = None
    if dens is not None:
        src_d = _read_source(dens, f"{kp}.density", manifest, base_dir, "density")
        kw["density_g_cm3"] = _load_npy(src_d, f"{kp}.density")
    geom = _construct(VoxelGeometry, kp, kw)
    items = _geometry_items(geom)
    mi = items["material_index"]
    _check_ref_content(src_mi, mi["content_sha256"], sha256_canonical(mi), f"{kp}.material_index")
    built.files.append(("material_index", src_mi.path, src_mi.sha256))
    if src_d is not None:
        _check_ref_content(
            src_d,
            items["density"]["content_sha256"],
            sha256_canonical(items["density"]),
            f"{kp}.density",
        )
        built.files.append(("density", src_d.path, src_d.sha256))
    return geom  # type: ignore[no-any-return]


def _resolve_material(name: str, mats: Mapping[str, Material], kp: str) -> Material:
    m = mats.get(name) or MATERIALS.get(name)
    if m is None:
        raise ConfigError(
            kp, f"unknown material {name!r} (not defined in [[materials]] and not built in)"
        )
    return m


def _lookups_from(
    v: Any, manifest: InputManifest | None, base_dir: Path | None, built: _Built
) -> tuple[LookupTable, ...]:
    out: list[LookupTable] = []
    for i, item in enumerate(_seq(v, "lookups")):
        kp = f"lookups[{i}]"
        declared = None
        if isinstance(item, dict) and "units" in item:  # written by to_dict; checked below
            declared = item["units"]
            item = {k: x for k, x in item.items() if k != "units"}
        src = _read_source(item, kp, manifest, base_dir, None)
        try:
            table = LookupTable.from_file(src.path, expected_sha256=src.sha256)
        except (UnsupportedCombinationError, OSError, ValueError, TypeError, KeyError) as exc:
            raise ConfigError(kp, f"{type(exc).__name__}: {exc}") from exc
        if declared is not None and declared != table.units:
            raise ConfigError(
                f"{kp}.units", f"declares {declared!r} but the table has units {table.units!r}"
            )
        role = f"lookup:{table.name}"
        if src.ref is not None and src.ref["ref"] != role:
            raise ConfigError(
                f"{kp}.ref", f"the file holds table {table.name!r}, not {src.ref['ref']!r}"
            )
        if any(f[0] == role for f in built.files):
            raise ConfigError(kp, f"duplicate lookup table name {table.name!r}")
        _check_ref_content(src, table.content_sha256, sha256_canonical(table.provenance()), kp)
        built.files.append((role, src.path, src.sha256))
        out.append(table)
    return tuple(out)


def _check_stopping_model(kw: dict[str, Any]) -> None:
    sm = kw["stopping_model"]
    registry = stopping_models.STOPPING_MODELS
    if sm not in registry:
        raise ConfigError(
            "physics.stopping_model",
            f"unknown stopping model id {sm!r}; registered: {sorted(registry)}",
        )


def _dataset_ids(
    stopping_model: str | None,
    pins: tuple[Any, ...],
    nuclear: bool,
    nuclear_table_id: str | None,
    cache_dir: str | Path | None,
) -> tuple[str, ...]:
    """Verify every cache object the stopping pins and the nuclear table need (``DatasetError``)."""
    del stopping_model
    cdir = cache.resolve_cache_dir(cache_dir)
    ids: list[str] = []

    def verify(dataset_id: str) -> None:
        try:
            cache.verify(dataset_id, cdir)
        except (FileNotFoundError, cache.IntegrityError, OSError, ValueError) as exc:
            raise DatasetError(f"cache object {dataset_id!r} in {cdir}: {exc}") from exc
        if dataset_id not in ids:
            ids.append(dataset_id)

    for pin in pins:
        if pin.use_role == "construction-permitted":
            verify(pin.dataset_id)
    if nuclear and nuclear_table_id is not None:
        from ionmc.nuclear.tables import NuclearTable, NuclearTableError  # noqa: PLC0415

        sidecar = cdir / "derived" / f"nuclear-proton-{nuclear_table_id}.json"
        try:
            sources = sorted(json.loads(sidecar.read_text(encoding="utf-8"))["sources"])
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise DatasetError(f"nuclear table {nuclear_table_id!r} in {cdir}: {exc}") from exc
        for sid in sources:
            verify(sid)
        try:
            NuclearTable.load(cdir, nuclear_table_id)
        except (NuclearTableError, OSError, ValueError, KeyError) as exc:
            raise DatasetError(f"nuclear table {nuclear_table_id!r} in {cdir}: {exc}") from exc
        ids.append(f"nuclear-proton-{nuclear_table_id}")
    return tuple(ids)


def verify_datasets(
    config: SimulationConfig, cache_dir: str | Path | None = None
) -> tuple[str, ...]:
    """Verify every cache object the selected stopping model and, if ``physics.nuclear``, the
    nuclear table need (``cache.verify``, re-hashing the objects).

    Returns the verified dataset ids. A missing, unreadable or mismatching object raises
    :class:`DatasetError` (``__cause__`` is the original error). ``cache_dir`` None resolves as the
    stopping and nuclear loaders do (environment variable, then the default).
    """
    ph = config.physics
    return _dataset_ids(None, ph.selection.spec.data, ph.nuclear, ph.nuclear_table_id, cache_dir)


def _physics_from(v: Any) -> PhysicsOptions:
    kp = "physics"
    t = _table(v, kp)
    if "stopping" in t:
        rec = t["stopping"]
        if isinstance(rec, dict) and "status" in rec:
            raise ConfigError(
                "physics.stopping",
                f"{rec.get('status')}: {rec.get('reason', 'a user stopping object')}; select a "
                "registered model with stopping_model",
            )
        raise ConfigError(
            "physics.stopping",
            "a stopping source object cannot be given in a configuration; use "
            "physics.stopping_model (a registered id)",
        )
    kw = _read(t, kp, PHYSICS_FIELDS)
    _check_stopping_model(kw)  # before PhysicsOptions is built
    spec = stopping_models.STOPPING_MODELS[kw["stopping_model"]]
    _dataset_ids(kw["stopping_model"], spec.data, kw["nuclear"], kw.get("nuclear_table_id"), None)
    return _construct(PhysicsOptions, kp, kw, dataset=True)  # type: ignore[no-any-return]


def _output_from(v: Any, base_dir: Path | None) -> OutputOptions:
    kw = _read(
        _table(v, "output"),
        "output",
        OUTPUT_FIELDS,
    )
    path = kw.get("path")
    if path is not None:
        p = Path(path)
        if not p.is_absolute():
            if base_dir is None:
                raise ConfigError("output.path", "a relative path needs a base_dir")
            p = base_dir / p
        path = str(p.resolve())
    return OutputOptions(path, kw.get("persist_batches", True))


def _beamlets_table(tbl: Mapping[str, Any]) -> tuple[Any, list[dict[str, Any]]]:
    """Strict read of a ``type = "beamlets"`` source table: (projectile, beamlet row tables)."""
    for k in tbl:
        if k not in ("type", "projectile", "beamlets"):
            raise ConfigError(
                f"source.{k}",
                "unknown key for type 'beamlets' (valid keys: type, projectile, beamlets)",
            )
    if "projectile" not in tbl:
        raise ConfigError("source.projectile", "required key is missing")
    if "beamlets" not in tbl:
        raise ConfigError("source.beamlets", "required key is missing")
    return tbl["projectile"], _table_list(tbl["beamlets"], "source.beamlets")


def _read_beamlet(t: Mapping[str, Any], kp: str) -> dict[str, Any]:
    if "projectile" in t:
        raise ConfigError(
            f"{kp}.projectile",
            "a beamlet set has one projectile (source.projectile); mixed projectiles are refused",
        )
    return _read(t, kp, BEAMLET_FIELDS)


def _source_from(v: Any, built: _Built | None = None) -> PencilBeamSource | BeamletSet:
    tbl = _table(v, "source")
    kind = tbl.get("type", "pencil")
    if kind not in SOURCE_TYPES:
        raise ConfigError("source.type", f"must be one of {SOURCE_TYPES}, got {kind!r}")
    if kind == "beamlets":
        proj, rows = _beamlets_table(tbl)
        projectile = _projectile_from(proj, "source.projectile")
        beamlets = tuple(
            _construct(Beamlet, f"source.beamlets[{i}]", _read_beamlet(t, f"source.beamlets[{i}]"))
            for i, t in enumerate(rows)
        )
        made: BeamletSet = _construct(
            BeamletSet, "source", {"projectile": projectile, "beamlets": beamlets}
        )
        return made
    sk = _read({k: x for k, x in tbl.items() if k != "type"}, "source", SOURCE_FIELDS)
    sk["projectile"] = _projectile_from(sk["projectile"], "source.projectile")
    return _construct(PencilBeamSource, "source", sk)  # type: ignore[no-any-return]


def _from_dict(
    d: Any,
    manifest: InputManifest | None,
    base_dir: Path | None,
    arrays: Mapping[str, Any] | None = None,
) -> tuple[SimulationConfig, _Built]:
    top = _table(d, "<root>")
    for k in top:
        if k not in TOP_KEYS:
            raise ConfigError(k, f"unknown top-level key (valid keys: {', '.join(TOP_KEYS)})")
    if top.get("schema") != SCHEMA:
        raise ConfigError("schema", f"must be {SCHEMA!r}, got {top.get('schema')!r}")
    if top.get("schema_version") != SCHEMA_VERSION:
        raise ConfigError(
            "schema_version", f"must be {SCHEMA_VERSION!r}, got {top.get('schema_version')!r}"
        )
    for k in REQUIRED_TOP:
        if k not in top:
            raise ConfigError(k, "required table is missing")
    built = _Built([], OutputOptions())
    mats: dict[str, Material] = {}
    for i, m in enumerate(_seq(top.get("materials", []), "materials")):
        mat = _material_from(m, f"materials[{i}]")
        if mat.name in mats:
            raise ConfigError(f"materials[{i}].name", f"duplicate material name {mat.name!r}")
        mats[mat.name] = mat
    source = _source_from(top["source"], built)
    geometry = _geometry_from(top["geometry"], mats, manifest, base_dir, built, arrays)
    scoring = tuple(
        _construct(ScoringGrid, f"scoring[{i}]", _read(t, f"scoring[{i}]", SCORING_FIELDS))
        for i, t in enumerate(_table_list(top["scoring"], "scoring"))
    )
    physics = _physics_from(top["physics"])
    run_tbl = dict(_table(top["run"], "run"))
    if isinstance(source, BeamletSet):  # [run] n_histories is optional: derived = sum of n_j
        total = source.n_histories_total
        if "n_histories" not in run_tbl:
            run_tbl["n_histories"] = total
            built.n_histories_derived = True
        elif run_tbl["n_histories"] != total:
            raise ConfigError(
                "run.n_histories",
                f"is {run_tbl['n_histories']!r} but the beamlet n_histories sum to {total}; "
                "omit it or give the sum",
            )
    run = _construct(RunOptions, "run", _read(run_tbl, "run", RUN_FIELDS))
    diagnostics = _construct(
        DiagnosticsOptions,
        "diagnostics",
        _read(_table(top.get("diagnostics", {}), "diagnostics"), "diagnostics", DIAGNOSTICS_FIELDS),
    )
    tallies = tuple(
        _construct(TallyRequest, f"tallies[{i}]", _read(t, f"tallies[{i}]", TALLY_FIELDS))
        for i, t in enumerate(_table_list(top.get("tallies", []), "tallies"))
    )
    lookups = _lookups_from(top.get("lookups", []), manifest, base_dir, built)
    influence = None
    if "influence" in top:
        influence = _construct(
            InfluenceOptions,
            "influence",
            _read(_table(top["influence"], "influence"), "influence", INFLUENCE_FIELDS),
        )
    if "output" in top:
        built.output = _output_from(top["output"], base_dir)
    config = _construct(
        SimulationConfig,
        "<config>",
        {
            "source": source,
            "geometry": geometry,
            "scoring": scoring,
            "physics": physics,
            "run": run,
            "diagnostics": diagnostics,
            "tallies": tallies,
            "lookups": lookups,
            "influence": influence,
        },
    )
    return config, built


def from_dict(
    d: Mapping[str, Any],
    *,
    manifest: InputManifest | None = None,
    base_dir: str | Path | None = None,
    arrays: Mapping[str, Any] | None = None,
) -> SimulationConfig:
    """Strict inverse of :func:`to_dict`.

    An unknown key at any level, a wrong type or a missing required key raises ``ConfigError``
    naming the full key path. References ``{"ref", "sha256", "content_sha256"}`` are reloaded
    through ``manifest`` (file re-hashed, both digests checked); ``not-reconstructible`` records
    raise. Array items of the file form ``{"file", "sha256"}`` resolve relative paths against
    ``base_dir``. ``physics.stopping_model`` and ``physics.nuclear`` are required; an unknown
    model id raises ``ConfigError`` before ``PhysicsOptions`` is built and the cache objects of
    the selected model are verified first (``DatasetError``).

    ``arrays`` (role -> numpy array, roles ``material_index`` and ``density``) resolves the
    geometry array positions from memory instead of files (the result reader passes the persisted
    ``.npy`` members): the content digest of each array must equal the record's
    ``content_sha256``; no file is opened for the geometry. Lookup tables are never resolved from
    ``arrays`` (their values are not persisted).
    """
    bd = None if base_dir is None else Path(base_dir).resolve()
    return _from_dict(d, manifest, bd, arrays)[0]


def _check_array_record(v: Any, kp: str, required: bool, units: bool = False) -> None:
    """A serialised array position: reference, not-reconstructible record, file form or null;
    status/reason are strings, ``ref`` a string and every digest 64 lower-case hex digits. A
    lookup position (``units=True``) may carry the string ``units`` of the table beside the
    record."""
    if v is None and not required:
        return
    t = _table(v, kp)
    if units and "units" in t:
        _str(t["units"], f"{kp}.units")
        t = {k: x for k, x in t.items() if k != "units"}
    forms = (
        {"ref", "sha256", "content_sha256"},
        {"status", "reason", "content_sha256"},
        {"file", "sha256"},
    )
    if set(t) not in forms:
        raise ConfigError(kp, f"unknown or missing key(s) in array record {sorted(t)}")
    for k, x in t.items():
        if k in ("sha256", "content_sha256"):
            if not (isinstance(x, str) and re.fullmatch(r"[0-9a-f]{64}", x)):
                raise ConfigError(f"{kp}.{k}", "must be 64 lower-case hex digits")
        elif k == "status":
            if x != "not-reconstructible":
                raise ConfigError(f"{kp}.status", "must be 'not-reconstructible'")
        else:
            _str(x, f"{kp}.{k}")


def validate_dict_structure(d: Any) -> None:
    """Structure-only validation of a :func:`to_dict` document (never opens a file).

    Same key knowledge as :func:`from_dict`: an unknown or missing key, a wrong type or an
    unknown schema/version raises ``ConfigError`` naming the key path. Array positions accept the
    reference, ``not-reconstructible`` and file records; no value is resolved or constructed.
    """
    top = _table(d, "<root>")
    for k in top:
        if k not in TOP_KEYS:
            raise ConfigError(k, f"unknown top-level key (valid keys: {', '.join(TOP_KEYS)})")
    if top.get("schema") != SCHEMA or top.get("schema_version") != SCHEMA_VERSION:
        raise ConfigError("schema", f"must be {SCHEMA!r} {SCHEMA_VERSION!r}")
    for k in REQUIRED_TOP:
        if k not in top:
            raise ConfigError(k, "required table is missing")
    for i, m in enumerate(_seq(top.get("materials", []), "materials")):
        t = _table(m, f"materials[{i}]")
        if set(t) != {"name"}:
            kw = _read(t, f"materials[{i}]", MATERIAL_FIELDS)
            _table(kw["mass_fractions"], f"materials[{i}].mass_fractions")
            if kw.get("sternheimer") is not None:
                st = _table(kw["sternheimer"], f"materials[{i}].sternheimer")
                _read(st, f"materials[{i}].sternheimer", STERNHEIMER_FIELDS)
    st = _table(top["source"], "source")
    if st.get("type", "pencil") == "beamlets":
        proj, rows = _beamlets_table(st)
        for i, t in enumerate(rows):
            _read_beamlet(t, f"source.beamlets[{i}]")
        sk = {"projectile": proj}
    else:
        sk = _read({k: x for k, x in st.items() if k != "type"}, "source", SOURCE_FIELDS)
    if not isinstance(sk["projectile"], str):
        _read(_table(sk["projectile"], "source.projectile"), "source.projectile", PROJECTILE_FIELDS)
    g = _table(top["geometry"], "geometry")
    if g.get("type") == "box":
        _read(g, "geometry", BOX_FIELDS)
    elif g.get("type") == "voxel":
        _read(g, "geometry", VOXEL_FIELDS)
        _check_array_record(g["material_index"], "geometry.material_index", True)
        _check_array_record(g.get("density"), "geometry.density", False)
    else:
        raise ConfigError("geometry.type", "must be 'box' or 'voxel'")
    for i, t in enumerate(_table_list(top["scoring"], "scoring")):
        _read(t, f"scoring[{i}]", SCORING_FIELDS)
    ph = _table(top["physics"], "physics")
    schema: dict[str, Any] = dict(PHYSICS_FIELDS)
    if "stopping" in ph:  # a user stopping object: record present <=> stopping_model is null
        _check_stopping_record(ph["stopping"])
        if ph.get("stopping_model", "absent") is not None:
            raise ConfigError("physics.stopping_model", "must be null beside a stopping record")
        schema["stopping_model"] = (_opt(_str), True)
    elif ph.get("stopping_model") is None:
        raise ConfigError("physics.stopping_model", "required (a string) without a stopping record")
    _read({k: v for k, v in ph.items() if k != "stopping"}, "physics", schema)
    _read(_table(top["run"], "run"), "run", RUN_FIELDS)
    _read(_table(top.get("diagnostics", {}), "diagnostics"), "diagnostics", DIAGNOSTICS_FIELDS)
    for i, t in enumerate(_table_list(top.get("tallies", []), "tallies")):
        _read(t, f"tallies[{i}]", TALLY_FIELDS)
    for i, r in enumerate(_seq(top.get("lookups", []), "lookups")):
        _check_array_record(r, f"lookups[{i}]", True, units=True)
    if "influence" in top:
        _read(_table(top["influence"], "influence"), "influence", INFLUENCE_FIELDS)
    if "output" in top:
        _read(_table(top["output"], "output"), "output", OUTPUT_FIELDS)


def _check_stopping_record(v: Any) -> None:
    """The exact ``not-reconstructible`` record :func:`to_dict` writes for a stopping object."""
    t = _table(v, "physics.stopping")
    want = {"status", "reason", "reconstructible", "stopping_object"}
    if set(t) != want:
        raise ConfigError(
            "physics.stopping", f"keys must be exactly {sorted(want)}, got {sorted(t)}"
        )
    if t["status"] != "not-reconstructible" or t["reconstructible"] is not False:
        raise ConfigError("physics.stopping", "must be a not-reconstructible record")
    _str(t["reason"], "physics.stopping.reason")
    so = _table(t["stopping_object"], "physics.stopping.stopping_object")
    if set(so) != {"class", "model_id", "spec_sha256"}:
        raise ConfigError(
            "physics.stopping.stopping_object", "keys must be exactly class, model_id, spec_sha256"
        )
    for k, x in so.items():
        _str(x, f"physics.stopping.stopping_object.{k}")


def _parse_toml(path: str | Path) -> tuple[dict[str, Any], Path, str]:
    p = Path(path)
    try:
        resolved = p.resolve()
        raw = resolved.read_bytes()
    except (OSError, ValueError) as exc:
        raise ConfigError(
            "<file>", f"cannot read the configuration file {p.name!r}: {exc}"
        ) from exc
    try:
        doc = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ConfigError("<file>", f"{p.name!r} is not valid TOML: {exc}") from exc
    return doc, resolved, hashlib.sha256(raw).hexdigest()


def load_toml_full(path: str | Path) -> LoadedConfig:
    """Load a TOML configuration: ``(config, manifest, [output] options)``.

    ``schema = "ionmc.config"`` and ``schema_version = "1.0"`` are required at the top level; the
    tables are ``[source]``, ``[geometry]`` (``type`` ``box`` or ``voxel``), ``[[materials]]``,
    ``[[scoring]]``, ``[physics]``, ``[run]``, ``[diagnostics]``, ``[[tallies]]``, ``[[lookups]]``
    and ``[output]``. Everything is read through the strict path of :func:`from_dict`.
    """
    doc, resolved, file_sha = _parse_toml(path)
    base = resolved.parent
    config, built = _from_dict(doc, None, base)
    items = _config_items(config)
    files = tuple(
        InputFile(role, p, sha, items[role][0], items[role][1], Path(p).stat().st_size)
        for role, p, sha in built.files
    )
    manifest = _finish_manifest(
        "toml", str(base), ConfigFileRecord(str(resolved), file_sha), files, config
    )
    return LoadedConfig(config, manifest, built.output, built.n_histories_derived)


def load_toml(path: str | Path) -> tuple[SimulationConfig, InputManifest]:
    """Load a TOML configuration: ``(SimulationConfig, InputManifest)`` (see
    :func:`load_toml_full`, which also returns the ``[output]`` options)."""
    loaded = load_toml_full(path)
    return loaded.config, loaded.manifest


def load_toml_dict(path: str | Path) -> dict[str, Any]:
    """The parsed TOML document as a dict, after the same strict validation as :func:`load_toml`
    (including the file hashes and the dataset checks). Relative paths are left as written;
    pass the file's directory as ``base_dir`` to :func:`from_dict`."""
    doc, resolved, _ = _parse_toml(path)
    _from_dict(doc, None, resolved.parent)
    return doc
