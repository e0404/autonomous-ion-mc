"""JSON shapes (keys and value types only, never values) of ``EffectiveConfig.summary()``.

Derived from the REG-2 summaries fixture (``tests/ionmc/data/v3_006/summaries-95937553.json``) by
:func:`derive`; a test regenerates the table from the fixture and compares it with ``SHAPES``.
For each top-level key: the set of JSON types seen, and for dict values the keys present in every
case (``required``), the keys seen in any case (``allowed``) and the types of each key.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


def jtype(v: Any) -> str:
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, int | float):
        return "number"
    if isinstance(v, str):
        return "str"
    if isinstance(v, list | tuple):
        return "list"
    return "dict"


def derive(summaries: list[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for key in sorted({k for s in summaries for k in s}):
        vals = [s[key] for s in summaries if key in s]
        spec: dict[str, Any] = {
            "required": len(vals) == len(summaries),
            "types": sorted({jtype(v) for v in vals}),
        }
        dicts = [v for v in vals if isinstance(v, dict)]
        if dicts:
            keys = {k for d in dicts for k in d}
            spec["keys"] = {
                k: {
                    "required": all(k in d for d in dicts),
                    "types": sorted({jtype(d[k]) for d in dicts if k in d}),
                }
                for k in sorted(keys)
            }
        out[key] = spec
    return out


SHAPES: dict[str, dict[str, Any]] = {
    "backend": {"required": True, "types": ["str"]},
    "chunk_histories": {"required": True, "types": ["number"]},
    "cpu_workers": {"required": True, "types": ["number"]},
    "geometry": {
        "keys": {
            "materials": {"required": True, "types": ["list"]},
            "origin_mm": {"required": True, "types": ["list"]},
            "shape": {"required": True, "types": ["list"]},
            "spacing_mm": {"required": True, "types": ["list"]},
            "z_exit_mm": {"required": True, "types": ["null"]},
        },
        "required": True,
        "types": ["dict"],
    },
    "max_steps": {"required": True, "types": ["number"]},
    "max_steps_origin": {"required": True, "types": ["str"]},
    "mcs_mapping": {
        "keys": {
            "depth_fractions": {"required": True, "types": ["list"]},
            "evaluation": {"required": False, "types": ["str"]},
            "fdm_domain_edge_mev": {"required": True, "types": ["number"]},
            "fdm_extrapolated": {"required": True, "types": ["bool"]},
            "mapping": {"required": True, "types": ["str"]},
            "p1v1_eq_kinetic_energy_hi_mev": {"required": True, "types": ["null", "number"]},
            "p1v1_eq_kinetic_energy_mev": {"required": True, "types": ["null", "number"]},
            "proton_table_edge_mev": {"required": True, "types": ["number"]},
            "table": {"required": False, "types": ["dict"]},
        },
        "required": True,
        "types": ["dict"],
    },
    "n_batches": {"required": True, "types": ["number"]},
    "n_histories": {"required": True, "types": ["number"]},
    "nuclear": {
        "keys": {
            "builder_version": {"required": True, "types": ["str"]},
            "deuteron_tables_sha256": {"required": True, "types": ["str"]},
            "e_cut_deuteron_mev": {"required": True, "types": ["number"]},
            "elements": {"required": True, "types": ["dict"]},
            "grid_size": {"required": True, "types": ["number"]},
            "history_energy_bound_mev": {"required": True, "types": ["number"]},
            "majorant_factor": {"required": True, "types": ["number"]},
            "majorant_window_f_e": {"required": True, "types": ["number"]},
            "multiplicity_model": {"required": True, "types": ["str"]},
            "npz_sha256": {"required": True, "types": ["str"]},
            "producible": {"required": True, "types": ["list"]},
            "q_plus_table_mev": {"required": True, "types": ["number"]},
            "recoil_t_max_mev": {"required": True, "types": ["number"]},
            "schema": {"required": True, "types": ["str"]},
            "source_sha256": {"required": True, "types": ["dict"]},
            "table_id": {"required": True, "types": ["str"]},
            "transport_energy_bound_mev": {"required": True, "types": ["number"]},
            "transport_path_bound_terms": {"required": True, "types": ["dict"]},
        },
        "required": False,
        "types": ["dict"],
    },
    "physics": {
        "keys": {
            "delta_electrons": {"required": True, "types": ["str"]},
            "e_cut_deuteron_mev": {"required": False, "types": ["number"]},
            "e_cut_effective_mev": {"required": True, "types": ["number"]},
            "e_cut_ion_mev_per_u": {"required": True, "types": ["number"]},
            "e_cut_mev": {"required": True, "types": ["number"]},
            "end_of_range_model": {"required": True, "types": ["str"]},
            "end_of_range_model_effective": {"required": True, "types": ["str"]},
            "ls_effective": {"required": True, "types": ["str"]},
            "ls_requested": {"required": True, "types": ["str"]},
            "max_energy_loss_fraction": {"required": True, "types": ["number"]},
            "max_rejection_attempts": {"required": True, "types": ["number"]},
            "max_step_mm": {"required": True, "types": ["number"]},
            "mcs_model": {"required": True, "types": ["str"]},
            "multiple_scattering": {"required": True, "types": ["bool"]},
            "nuclear": {"required": True, "types": ["bool"]},
            "range_alpha": {"required": True, "types": ["number"]},
            "range_rho_f_mm": {"required": True, "types": ["number"]},
            "scattering_E_s_mev": {"required": True, "types": ["number"]},
            "short_step_fraction": {"required": True, "types": ["number"]},
            "straggling": {"required": True, "types": ["bool"]},
            "straggling_model": {"required": True, "types": ["str"]},
            "truncated_hinge_diagnostic": {"required": True, "types": ["bool"]},
        },
        "required": True,
        "types": ["dict"],
    },
    "precision": {"required": True, "types": ["str"]},
    "production": {"required": True, "types": ["bool"]},
    "rng": {
        "keys": {
            "counter": {"required": True, "types": ["str"]},
            "generator": {"required": True, "types": ["str"]},
            "key": {"required": True, "types": ["str"]},
            "u01": {"required": True, "types": ["str"]},
        },
        "required": True,
        "types": ["dict"],
    },
    "scattering_length_g_cm2": {"required": True, "types": ["list"]},
    "scoring": {"required": True, "types": ["list"]},
    "scoring_pieces": {"required": True, "types": ["number"]},
    "seed": {"required": True, "types": ["number"]},
    "source": {
        "keys": {
            "a": {"required": True, "types": ["number"]},
            "energy_sigma_mev": {"required": True, "types": ["number"]},
            "kinetic_energy_mev": {"required": True, "types": ["number"]},
            "lateral_sigma_mm": {"required": True, "types": ["number"]},
            "mass_mev": {"required": True, "types": ["number"]},
            "position_mm": {"required": True, "types": ["list"]},
            "projectile": {"required": True, "types": ["str"]},
            "z": {"required": True, "types": ["number"]},
        },
        "required": True,
        "types": ["dict"],
    },
    "stopping": {
        "keys": {
            "cache_dir": {"required": True, "types": ["null", "str"]},
            "circular_rows": {"required": True, "types": ["list"]},
            "content_sha256": {"required": True, "types": ["null", "str"]},
            "data_declared": {"required": True, "types": ["list"]},
            "data_used": {"required": True, "types": ["list"]},
            "decision": {"required": True, "types": ["str"]},
            "domain": {"required": True, "types": ["dict"]},
            "domain_resolved": {"required": False, "types": ["dict"]},
            "evidence_at_run": {"required": True, "types": ["list"]},
            "interpolation": {"required": True, "types": ["str"]},
            "kind": {"required": True, "types": ["str"]},
            "mcs_proton_path_stopping": {"required": False, "types": ["str"]},
            "model_id_effective": {"required": True, "types": ["str"]},
            "model_id_requested": {"required": True, "types": ["str"]},
            "options": {"required": True, "types": ["dict"]},
            "range_construction": {"required": True, "types": ["str"]},
            "range_start": {"required": True, "types": ["str"]},
            "selected_by": {"required": True, "types": ["str"]},
            "spec_sha256": {"required": True, "types": ["str"]},
            "species_transported": {"required": True, "types": ["list"]},
            "version": {"required": True, "types": ["str"]},
        },
        "required": True,
        "types": ["dict"],
    },
    "tables": {
        "keys": {
            "materials": {"required": True, "types": ["list"]},
            "n_e": {"required": True, "types": ["number"]},
            "n_r": {"required": True, "types": ["number"]},
            "sha256": {"required": True, "types": ["str"]},
            "stopping_source": {"required": True, "types": ["str"]},
            "water_row": {"required": False, "types": ["dict"]},
        },
        "required": True,
        "types": ["dict"],
    },
    "tallies": {
        "keys": {
            "bounds": {"required": True, "types": ["dict"]},
            "channels": {"required": True, "types": ["list"]},
            "histories_per_batch": {"required": True, "types": ["number"]},
            "lookups": {"required": True, "types": ["list"]},
            "memory_bytes": {"required": True, "types": ["number"]},
            "n_channels": {"required": True, "types": ["number"]},
            "n_residual_columns": {"required": True, "types": ["number"]},
            "quantities": {"required": True, "types": ["list"]},
            "total_size": {"required": True, "types": ["number"]},
        },
        "required": True,
        "types": ["dict", "null"],
    },
    "unit_direction": {"required": True, "types": ["list"]},
}
