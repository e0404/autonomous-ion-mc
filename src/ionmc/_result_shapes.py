"""Declared schema of ``EffectiveConfig.summary()`` and of the provenance records of a result.

Written by hand from ``ionmc.config.EffectiveConfig.summary`` (and the records it embeds:
``ChannelPlan.summary``, ``NuclearSetup.summary``, ``identity_block``, ``FdmDomainStatus.to_dict``,
``IonMapping.table_summary``, ``LookupTable.provenance``), NOT derived from observed values: a
fixture test only checks that every observed summary conforms to this schema, never the reverse.

Mini-language (a spec is a string or a tuple):

``"int"`` exact integer (a float or bool is refused), ``"num"`` finite number (int or float, not
bool), ``"str"``, ``"bool"``, ``"null"``, ``"any"`` any strict JSON value (only for records whose
interior is a verbatim copy of an engine table, documented per use);
``("or", a, b, ...)`` one of the specs; ``("list", elem, n)`` list of ``elem`` (``n`` None or the
exact length); ``("map", elem)`` object with arbitrary string keys and ``elem`` values;
``("rec", {key: spec}, optional_keys)`` object with exactly the listed keys (the optional ones may
be absent).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

Spec = Any


def opt(x: Spec) -> Spec:
    return ("or", x, "null")


def lst(x: Spec, n: int | None = None) -> Spec:
    return ("list", x, n)


def mapping(x: Spec) -> Spec:
    return ("map", x)


def rec(fields: dict[str, Spec], optional: Sequence[str] = ()) -> Spec:
    return ("rec", fields, frozenset(optional))


class SchemaError(ValueError):
    """A value does not conform to its declared spec."""


def _is_int(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _is_num(v: Any) -> bool:
    return isinstance(v, int | float) and not isinstance(v, bool) and math.isfinite(v)


def _matches(v: Any, spec: Spec, path: str) -> bool:
    try:
        check(v, spec, path)
    except SchemaError:
        return False
    return True


def check(v: Any, spec: Spec, path: str = "") -> None:
    """Raise :class:`SchemaError` (naming ``path``) unless ``v`` conforms to ``spec``."""
    if isinstance(spec, str):
        ok = {
            "int": _is_int,
            "num": _is_num,
            "str": lambda x: isinstance(x, str),
            "bool": lambda x: isinstance(x, bool),
            "null": lambda x: x is None,
            "any": lambda x: True,
        }[spec](v)
        if not ok:
            raise SchemaError(f"{path}: expected {spec}, got {type(v).__name__}")
        return
    tag = spec[0]
    if tag == "or":
        if not any(_matches(v, s, path) for s in spec[1:]):
            raise SchemaError(f"{path}: matches none of the allowed forms")
    elif tag == "list":
        if not isinstance(v, list):
            raise SchemaError(f"{path}: expected a list")
        if spec[2] is not None and len(v) != spec[2]:
            raise SchemaError(f"{path}: expected {spec[2]} elements, got {len(v)}")
        for i, x in enumerate(v):
            check(x, spec[1], f"{path}[{i}]")
    elif tag == "map":
        if not isinstance(v, dict):
            raise SchemaError(f"{path}: expected an object")
        for k, x in v.items():
            if not isinstance(k, str):
                raise SchemaError(f"{path}: non-string key")
            check(x, spec[1], f"{path}.{k}")
    else:
        fields, optional = spec[1], spec[2]
        if not isinstance(v, dict):
            raise SchemaError(f"{path}: expected an object")
        extra = sorted(set(v) - set(fields))
        if extra:
            raise SchemaError(f"{path}: unknown field(s) {extra}")
        for k, s in fields.items():
            if k not in v:
                if k not in optional:
                    raise SchemaError(f"{path}: missing required field {k!r}")
                continue
            check(v[k], s, f"{path}.{k}")


# ---------------------------------------------------------------------------------------------
# records embedded in the summary

_VEC3 = lst("num", 3)
_PIN = mapping(
    opt("str")
)  # a dataset pin: dataset_id, archive_sha256, licence, member, parser, ...
_EVIDENCE = rec(
    {
        "spec_sha256": opt("str"),
        "status": "str",
        "row": "any",  # a plan row id (str) or null today; free until decision 0044 pins it
        "plan": "str",
        "graded_at_sha": opt("str"),
    }
)
LOOKUP_PROVENANCE = rec(
    {
        "name": "str",
        "quantity": "str",
        "units": "str",
        "axis": "str",
        "axis_spacing": "str",
        "axis_range": lst("num", 2),
        "n_points": "int",
        "species": lst("str"),
        "file_sha256": opt("str"),
        "content_sha256": "str",
        "citation": opt("str"),
        "license": opt("str"),
        "source": opt("str"),
        "synthetic": "bool",
        "resampling": "any",  # free JSON record of the resampling (or null)
    }
)
_SPECTRUM = rec({"edges": lst("num"), "log": "bool"})
_CHANNEL = rec(
    {
        "kind": "str",
        "grid": "int",
        "class_mask": "int",
        "species": lst("int"),
        "generation": lst("int", 2),
        "lookup": "int",
        "spectrum": opt(_SPECTRUM),
        "size": "int",
        "offset": "int",
        "k": "int",
        "bound_per_history": "num",
    }
)
_QUANTITY = rec(
    {
        "name": "str",
        "kind": "str",
        "numerator": "int",
        "denominator": opt("int"),
        "units": "str",
        "definition": "str",
    }
)
TALLIES = rec(
    {
        "n_channels": "int",
        "total_size": "int",
        "n_residual_columns": "int",
        "histories_per_batch": "int",
        "memory_bytes": "int",
        "bounds": mapping(opt("num")),
        "lookups": lst(LOOKUP_PROVENANCE),
        "channels": lst(_CHANNEL),
        "quantities": lst(_QUANTITY),
    }
)
NUCLEAR = rec(
    {
        "table_id": "str",
        "schema": "str",
        "builder_version": "str",
        "npz_sha256": "str",
        "source_sha256": mapping("str"),
        "elements": mapping(rec({"target": "str", "surrogate": "bool", "sigma_scale": "num"})),
        "q_plus_table_mev": "num",
        "transport_energy_bound_mev": "num",
        "history_energy_bound_mev": "num",
        "recoil_t_max_mev": "num",
        "transport_path_bound_terms": "any",  # thawed table record (species -> bound terms)
        "grid_size": "int",
        "multiplicity_model": "str",
        "majorant_factor": "num",
        "majorant_window_f_e": "num",
        "e_cut_deuteron_mev": "num",
        "deuteron_tables_sha256": "str",
        "producible": lst(rec({"species": "str", "generation": "any"})),
    },
    ("cache_dir",),  # a path: stripped from provenance.data.nuclear, kept (null/str) in effective
)
_TABLE_SUMMARY = rec(
    {
        "kind": "str",
        "n_pv": "int",
        "n_r": "int",
        "ln_pv_range": lst("num", 2),
        "ln_r_range": lst("num", 2),
        "pv_lo_mev": "num",
        "p1v1_min_mev": "num",
        "p1v1_max_mev": "num",
        "n_materials": "int",
        "sha256": "str",
        "contract": "str",
        "root_residual_tol": "num",
        "root_step_tol": "num",
        "bracket_factor": "num",
        "bracket_min": "num",
        "bracket_half_width_per_material": lst("num"),
        "interp_estimate": opt("num"),
        "interp_estimate_per_material": lst(opt("num")),
        "interp_estimate_note": "str",
        "safety_factor": "num",
        "estimated_cells": "int",
        "refined_points": "int",
        "table_coverage": "num",
        "table_coverage_per_material": lst("num"),
        "table_coverage_note": "str",
    },
    (
        "contract root_residual_tol root_step_tol bracket_factor bracket_min "
        "bracket_half_width_per_material interp_estimate interp_estimate_per_material "
        "interp_estimate_note safety_factor estimated_cells refined_points table_coverage "
        "table_coverage_per_material table_coverage_note"
    ).split(),
)
_MCS_MAPPING = rec(
    {
        "mapping": "str",
        "fdm_extrapolated": "bool",
        "p1v1_eq_kinetic_energy_mev": opt("num"),
        "p1v1_eq_kinetic_energy_hi_mev": opt("num"),
        "fdm_domain_edge_mev": "num",
        "proton_table_edge_mev": "num",
        "depth_fractions": lst("num"),
        "fdm_domain": "str",
        "evaluation": "str",
        "unused": "bool",
        "unused_note": "str",
        "table": _TABLE_SUMMARY,
    },
    ("fdm_domain", "evaluation", "unused", "unused_note", "table"),
)
_STOPPING = rec(
    {
        "model_id_requested": "str",
        "model_id_effective": "str",
        "version": "str",
        "kind": "str",
        "spec_sha256": "str",
        "options": mapping("any"),  # the selection options of the model (model-specific keys)
        "interpolation": "str",
        "range_construction": "str",
        "range_start": "str",
        "domain": mapping("any"),  # per-species window and material rule (model-specific)
        "data_declared": lst(_PIN),
        "data_used": lst(
            rec(
                {
                    "pin": _PIN,
                    "role_in_this_run": "str",
                    "table": lst("str", 2),
                    "cache_dir": opt("str"),
                }
            )
        ),
        "content_sha256": opt("str"),
        "evidence_at_run": lst(_EVIDENCE),
        "decision": "str",
        "circular_rows": lst("str"),
        "selected_by": "str",
        "species_transported": lst("str"),
        "cache_dir": opt("str"),
        "domain_resolved": mapping(rec({"e_min_per_u": "num", "e_max_per_u": "num"})),
        "mcs_proton_path_stopping": "str",
    },
    ("domain_resolved", "mcs_proton_path_stopping"),
)
_PHYSICS = rec(
    {
        "nuclear": "bool",
        "straggling": "bool",
        "straggling_model": "str",
        "multiple_scattering": "bool",
        "mcs_model": "str",
        "delta_electrons": "str",
        "e_cut_mev": "num",
        "e_cut_ion_mev_per_u": "num",
        "e_cut_effective_mev": "num",
        "end_of_range_model": "str",
        "end_of_range_model_effective": "str",
        "max_step_mm": "num",
        "max_energy_loss_fraction": "num",
        "range_alpha": "num",
        "range_rho_f_mm": "num",
        "short_step_fraction": "num",
        "truncated_hinge_diagnostic": "bool",
        "scattering_E_s_mev": "num",
        "max_rejection_attempts": "int",
        "ls_requested": "str",
        "ls_effective": "str",
        "e_cut_deuteron_mev": "num",
    },
    ("e_cut_deuteron_mev",),
)

EFFECTIVE: dict[str, Spec] = {
    "backend": "str",
    "precision": "str",
    "production": "bool",
    "seed": "int",
    "n_histories": "int",
    "n_batches": "int",
    "cpu_workers": "int",
    "chunk_histories": "int",
    "max_steps": "int",
    "scoring_pieces": "int",
    "max_steps_origin": "str",
    "unit_direction": _VEC3,
    "source": rec(
        {
            "projectile": "str",
            "z": "int",
            "a": "int",
            "mass_mev": "num",
            "position_mm": _VEC3,
            "kinetic_energy_mev": "num",
            "energy_sigma_mev": "num",
            "lateral_sigma_mm": "num",
        }
    ),
    "physics": _PHYSICS,
    "mcs_mapping": _MCS_MAPPING,
    "geometry": rec(
        {
            "origin_mm": _VEC3,
            "spacing_mm": _VEC3,
            "shape": lst("int", 3),
            "materials": lst("str"),
            "z_exit_mm": opt("num"),  # null (no exit plane) or a finite number [mm]
        }
    ),
    "scoring": lst(
        rec({"name": "str", "origin_mm": _VEC3, "spacing_mm": _VEC3, "shape": lst("int", 3)})
    ),
    "tables": rec(
        {
            "sha256": "str",
            "stopping_source": "str",
            "materials": lst(mapping("any")),  # one identity record per material (engine table)
            "n_e": "int",
            "n_r": "int",
            "water_row": mapping("any"),  # identity of the LET water row (tallies only)
        },
        ("water_row",),
    ),
    "stopping": _STOPPING,
    "scattering_length_g_cm2": lst("num"),
    "rng": rec({"generator": "str", "key": "str", "counter": "str", "u01": "str"}),
    "tallies": opt(TALLIES),
    "nuclear": NUCLEAR,
}
EFFECTIVE_OPTIONAL = frozenset({"nuclear"})  # present only for nuclear=True runs
EFFECTIVE_TOP_KEYS = tuple(EFFECTIVE)

# provenance records (the writer builds them in ionmc.results)
PROVENANCE_CODE = rec(
    {
        "ionmc_version": opt("str"),
        "git_commit": opt("str"),
        "git_dirty": opt("bool"),
        "package_tree_oid": opt("str"),
        "source": "str",
        "status": "str",
        "reason": "str",
    }
)
PROVENANCE_ENVIRONMENT = rec(
    {
        "python_version": opt("str"),
        "platform": opt("str"),
        "numpy_version": opt("str"),
        "warp_version": opt("str"),
        "warp_devices": lst("any"),
        "warp_error": opt("str"),
    }
)
PROVENANCE_DATA = rec(
    {
        "stopping_data_used": lst(
            rec({"pin": _PIN, "role_in_this_run": "str", "table": lst("str", 2)})
        ),
        "stopping_content_sha256": opt("str"),
        "table_sha256": "str",
        "table_stopping_source": "str",
        "table_materials": lst(mapping("any")),
        "nuclear": opt(NUCLEAR),
        "lookups": lst(LOOKUP_PROVENANCE),
    }
)
PROVENANCE_MODELS = rec(
    {
        "stopping_model_id_effective": opt("str"),
        "stopping_version": opt("str"),
        "stopping_spec_sha256": opt("str"),
        "evidence_at_run": opt(lst(_EVIDENCE)),
        "straggling_model": opt("str"),
        "mcs_model": opt("str"),
        "end_of_range_model": opt("str"),
        "delta_electrons": opt("str"),
        "nuclear_table_id": opt("str"),
    }
)


def check_effective(eff: Any) -> None:
    """Raise :class:`SchemaError` unless ``eff`` conforms to the declared summary schema."""
    if not isinstance(eff, dict):
        raise SchemaError("effective: expected an object")
    extra = sorted(set(eff) - set(EFFECTIVE))
    if extra:
        raise SchemaError(f"effective: unknown field(s) {extra}")
    for k, spec in EFFECTIVE.items():
        if k not in eff:
            if k in EFFECTIVE_OPTIONAL:
                continue
            raise SchemaError(f"effective: missing required field {k!r}")
        check(eff[k], spec, f"effective.{k}")
