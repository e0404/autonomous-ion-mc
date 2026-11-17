"""Single selection point of the stopping source of a run (decision 0044, sections 2 to 5; S2a).

``PhysicsOptions`` calls :func:`select` once: exactly one of ``stopping`` (a source object,
identified from its instance fields by :func:`~ionmc.physics.stopping_models.spec_of`) and
``stopping_model`` (a registered id, resolved by :func:`~ionmc.physics.stopping_models.resolve`) is
given. ``validate`` derives every stopping-dependent object from that one source and then asserts,
over the table identities that came out, that they all name one effective model
(:func:`assert_single_model`). Before transport :func:`check_species_and_materials` and
:func:`check_window` fail closed on the model's declared domain (no extrapolation, no fallback
model, no automatic composite).

:func:`identity_block` builds the persisted ``EffectiveConfig.summary()["stopping"]`` block and
:func:`validate_identity_block` is its validator (used by the PER-1 round trip and, later, PER-3).

Effective id (row ROLE-3): a z = 1 projectile with a Lindhard-Sorensen variant runs the default's
tables bit for bit (the proton model never changes, decision 0043 section 1), so its effective id
is the default's; the requested id is recorded unchanged.

A ``stopping`` object of a type ``spec_of`` does not support is refused (row DOM-7) unless it
declares itself with a ``stopping_spec()`` method; it is then recorded as ``unregistered``.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ionmc._validate import fail
from ionmc.errors import UnsupportedCombinationError
from ionmc.materials import Material
from ionmc.physics import stopping_models as _models
from ionmc.physics.projectiles import Projectile
from ionmc.physics.stopping import StoppingSource
from ionmc.physics.stopping_models import (
    DEFAULT_STOPPING_MODEL_ID,
    MATERIAL_RULE_ANY,
    MATERIAL_RULE_WATER,
    UNREGISTERED_ID,
    StoppingModelSpec,
    spec_of,
    spec_sha256,
)

# Lindhard-Sorensen variants: for z = 1 the effective id is the default's (row ROLE-3).
LS_VARIANT_IDS = ("bethe-ionmc-v1+ls-eff", "bethe-ionmc-v1+ls-bare_mixed")
MCS_PROTON_PATH_STOPPING = "bethe-ionmc-v1"
BLOCK_KEYS = (
    "model_id_requested",
    "model_id_effective",
    "version",
    "kind",
    "spec_sha256",
    "options",
    "interpolation",
    "range_construction",
    "range_start",
    "domain",
    "data_declared",
    "data_used",
    "content_sha256",
    "evidence_at_run",
    "decision",
    "circular_rows",
    "selected_by",
    "species_transported",
)
_BETHE_FLAGS = ("shell", "barkas", "bloch", "mott", "effective_charge")
_REL = 1.0e-12


@dataclass(frozen=True, eq=False)
class StoppingSelection:
    """The one stopping source of a run and its identity.

    ``spec`` is the registered spec or the user-declared spec of the source (an unsupported type is
    refused, see the module docstring); ``selected_by`` is ``"stopping_model"`` or
    ``"stopping-object"``.
    """

    requested_id: str
    selected_by: str
    source: StoppingSource
    spec: StoppingModelSpec
    source_class: str

    def effective_id(self, projectile: Projectile) -> str:
        return effective_model_id(self.requested_id, projectile)

    def effective_spec(self, projectile: Projectile) -> StoppingModelSpec:
        """The spec in effect for ``projectile`` (the default's for a z = 1 LS variant)."""
        eid = self.effective_id(projectile)
        if eid != self.requested_id:
            return _models.STOPPING_MODELS[eid]
        return self.spec


def effective_model_id(requested_id: str, projectile: Projectile) -> str:
    """Effective id: the requested one, except a z = 1 projectile with an LS variant (ROLE-3)."""
    if projectile.z == 1 and requested_id in LS_VARIANT_IDS:
        return DEFAULT_STOPPING_MODEL_ID
    return requested_id


def select(stopping: StoppingSource | None, stopping_model: str | None) -> StoppingSelection:
    """The selection point: exactly one of ``stopping`` and ``stopping_model`` (row DOM-7)."""
    if (stopping is None) == (stopping_model is None):
        raise fail(
            "exactly one of PhysicsOptions.stopping (a source object) and "
            "PhysicsOptions.stopping_model (a registered id) must be given, got "
            f"stopping={'None' if stopping is None else type(stopping).__name__}, "
            f"stopping_model={stopping_model!r}"
        )
    if stopping_model is not None:
        if not isinstance(stopping_model, str):
            raise fail(f"stopping_model must be a string id, got {stopping_model!r}")
        try:
            source = _models.resolve(stopping_model)
        except (ValueError, OSError, RuntimeError) as exc:
            # unknown id, or (tabulated models) a missing / corrupt / unparsable cache object: row
            # DOM-6, fail closed with the loader's message (IntegrityError is a RuntimeError)
            raise fail(f"{type(exc).__name__}: {exc}") from exc
        return StoppingSelection(
            stopping_model,
            "stopping_model",
            source,
            _models.STOPPING_MODELS[stopping_model],
            type(source).__name__,
        )
    assert stopping is not None
    cls = type(stopping).__qualname__
    try:
        model_id, spec = spec_of(stopping)
    except (TypeError, ValueError) as exc:  # DOM-7: an unsupported source type fails closed
        raise fail(f"unsupported stopping object {type(stopping).__qualname__}: {exc}") from exc
    return StoppingSelection(model_id, "stopping-object", stopping, spec, cls)


# ---- fail-closed domain (decision 0044 section 3) ----


def _domain_text(spec: StoppingModelSpec) -> str:
    d = spec.domain
    return (
        f"species {sorted(d.species)}, material rule {d.material_rule!r}, "
        f"MeV/u windows {dict(sorted(d.e_min_per_u.items()))} to "
        f"{dict(sorted(d.e_max_per_u.items()))}"
    )


def _refuse(
    sel: StoppingSelection, species: str, material: str, energy: str
) -> UnsupportedCombinationError:
    return fail(
        f"stopping model {sel.requested_id!r} does not support species {species!r}, material "
        f"{material!r}, energy {energy}; its domain: {_domain_text(sel.spec)} (no extrapolation, "
        "no fallback model, no automatic composite)"
    )


def check_species_and_materials(
    sel: StoppingSelection, species: Iterable[Projectile], materials: Iterable[Material]
) -> None:
    """Species in the domain and the material rule for every material; else raise."""
    spec = sel.spec
    mats = tuple(materials)
    for sp in species:
        if sp.name not in spec.domain.species:
            raise _refuse(sel, sp.name, ",".join(m.name for m in mats), "any")
        for m in mats:
            if spec.domain.material_rule == MATERIAL_RULE_ANY:
                continue
            if spec.domain.material_rule == MATERIAL_RULE_WATER:
                from ionmc.physics.stopping import _require_nist_water  # noqa: PLC0415

                try:
                    _require_nist_water(m)
                except ValueError as exc:
                    raise _refuse(sel, sp.name, m.name, "any") from exc
                continue
            raise _refuse(
                sel, sp.name, m.name, f"any (unknown material rule {spec.domain.material_rule!r})"
            )


def check_window(
    sel: StoppingSelection,
    species: Projectile,
    material: str,
    what: str,
    e_lo_mev: float,
    e_hi_mev: float,
) -> None:
    """The total-energy interval ``[e_lo, e_hi]`` of ``what`` lies inside the model's per-species
    window (MeV/u windows times the mass number); else raise."""
    spec = sel.spec
    lo_u = spec.domain.e_min_per_u.get(species.name)
    hi_u = spec.domain.e_max_per_u.get(species.name)
    a = float(species.a)
    energy = f"{what} [{e_lo_mev:g}, {e_hi_mev:g}] MeV ({e_lo_mev / a:g} to {e_hi_mev / a:g} MeV/u)"
    if lo_u is None or hi_u is None:
        raise _refuse(sel, species.name, material, energy)
    resolved = getattr(sel.source, "window_per_u", None)
    if callable(resolved):  # a tabulated model: the window of the tables actually built
        r_lo, r_hi = resolved(species)
        lo_u, hi_u = max(lo_u, r_lo), min(hi_u, r_hi)
    if e_lo_mev / a < lo_u * (1.0 - _REL) or e_hi_mev / a > hi_u * (1.0 + _REL):
        raise _refuse(sel, species.name, material, energy)


# ---- single-model identity assertion ----


def _bethe_view(
    options: Mapping[str, Any], e_min_mev: float, e_max_mev: float, projectile: Projectile
) -> tuple[Any, ...]:
    ls = str(options.get("ls", "off")) if projectile.z >= 2 else "off"
    return (
        *(bool(options[k]) for k in _BETHE_FLAGS),
        ls,
        e_min_mev / projectile.a,
        e_max_mev / projectile.a,
    )


def _same_view(a: tuple[Any, ...], b: tuple[Any, ...]) -> bool:
    return a[:6] == b[:6] and all(
        math.isclose(x, y, rel_tol=_REL) for x, y in zip(a[6:], b[6:], strict=True)
    )


def table_model_id(entry: Mapping[str, Any], projectile: Projectile) -> str:
    """The model id a table identity entry (material row or water row) names.

    ``metadata["model_id"]`` if present (every model registered after 0044 carries it); else an
    analytic Bethe table is matched against the registered Bethe specs on the fields visible in
    its metadata and identity (options, effective Lindhard-Sorensen mode, grid ends); anything else
    is ``"unregistered"``. The registered models of 0044 keep their metadata byte-identical, so
    their ids are recovered here (the (metadata, projectile) -> effective id map).
    """
    meta = entry["metadata"]
    if "model_id" in meta:
        return str(meta["model_id"])
    if meta.get("source") == "nist-star":  # pre-freeze metadata, kept byte-identical (ID-NIST)
        for mid, spec in _models.STOPPING_MODELS.items():
            pins = {p.dataset_id: p.archive_sha256 for p in spec.data}
            if (
                projectile.name in spec.domain.species
                and meta.get("dataset_id") in pins
                and meta.get("sha256") == pins[meta["dataset_id"]]
            ):
                return mid
        return UNREGISTERED_ID
    if meta.get("source") != "bethe" or not all(k in meta.get("options", {}) for k in _BETHE_FLAGS):
        return UNREGISTERED_ID
    seen = _bethe_view(
        meta["options"], float(entry["e_min_mev"]), float(entry["e_max_mev"]), projectile
    )
    for mid, spec in _models.STOPPING_MODELS.items():
        if spec.kind != "analytic" or projectile.name not in spec.domain.species:
            continue
        want = _bethe_view(
            spec.options,
            spec.domain.e_min_per_u[projectile.name] * projectile.a,
            spec.domain.e_max_per_u[projectile.name] * projectile.a,
            projectile,
        )
        if _same_view(seen, want):
            return mid
    return UNREGISTERED_ID


def expected_table_model_id(sel: StoppingSelection, projectile: Projectile) -> str | None:
    """The id :func:`table_model_id` must return for the tables of ``sel`` (from the spec, by the
    same rule): the effective id of a registered selection; for an unregistered Bethe source the
    registered id with the same visible fields, else ``"unregistered"``."""
    spec = sel.effective_spec(projectile)
    if sel.requested_id != UNREGISTERED_ID:
        return sel.effective_id(projectile)
    if spec.kind != "analytic" or "barkas" not in spec.options:
        return None  # a user-declared source: only the uniformity of the ids is checked
    name = projectile.name
    entry = {
        "metadata": {"source": "bethe", "options": dict(spec.options)},
        "e_min_mev": spec.domain.e_min_per_u.get(name, math.nan) * projectile.a,
        "e_max_mev": spec.domain.e_max_per_u.get(name, math.nan) * projectile.a,
    }
    return table_model_id(entry, projectile)


def assert_single_model(
    sel: StoppingSelection,
    groups: Sequence[tuple[str, Projectile, Iterable[Mapping[str, Any] | None]]],
) -> str:
    """Every table identity of the run names the model of ``sel`` (one effective model).

    ``groups`` lists ``(label, projectile, identity entries)`` for the material tables, the LET
    water row and the deuteron tables of a nuclear run (``None`` entries are skipped). Raises
    ``UnsupportedCombinationError`` naming the label and the ids found; returns the effective id of
    the source projectile's group. A real check over the identities (rows CC-ID, CC-NEG).
    """
    found: dict[str, set[str]] = {}
    first_id = ""
    seen_ids: dict[str, set[str]] = {}
    for label, proj, entries in groups:
        want = expected_table_model_id(sel, proj)
        for entry in entries:
            if entry is None:
                continue
            got = table_model_id(entry, proj)
            if want is None:  # user-declared source: no id to compare with, ids must agree
                seen_ids.setdefault(got, set()).add(label)
                continue
            first_id = first_id or want
            if got != want:
                found.setdefault(label, set()).add(got)
    if len(seen_ids) > 1:
        raise fail(f"the tables of this run name different models: {sorted(seen_ids)}")
    if found:
        raise fail(
            f"the tables of this run do not name one stopping model: expected "
            f"{sel.requested_id!r} (effective {first_id!r}); found "
            + "; ".join(f"{k}: {sorted(v)}" for k, v in sorted(found.items()))
        )
    return first_id


# ---- persisted identity block (decision 0044 section 4) ----


def _canon(value: Any) -> Any:
    return _models._canonical(value)


def _evidence(spec: StoppingModelSpec, digest: str | None) -> list[dict[str, Any]]:
    from ionmc.physics.stopping_evidence import STOPPING_EVIDENCE  # noqa: PLC0415

    if digest is not None and digest in STOPPING_EVIDENCE:
        return [_canon(r) for r in STOPPING_EVIDENCE[digest]]
    return [
        {
            "spec_sha256": digest,
            "status": "unvalidated",
            "row": None,
            "plan": "none (a user source object, decision 0044 section 2)",
            "graded_at_sha": None,
        }
    ]


def identity_block(
    sel: StoppingSelection,
    projectile: Projectile,
    species_transported: Sequence[str],
    *,
    data_used: Sequence[Mapping[str, Any]] = (),
    content_sha256: str | None = None,
    range_construction: str | None = None,
    materials: Sequence[Material] = (),
) -> dict[str, Any]:
    """The block ``EffectiveConfig.summary()["stopping"]`` (outside the table hash).

    ``data_used`` lists, for tabulated models, only the pins read for the resolved projectile and
    material tables of the run, each as ``{"pin": ..., "role_in_this_run": "construction",
    "table": [projectile, material]}``; ``content_sha256`` is the digest of the parsed arrays. Both
    stay empty / ``None`` for the analytic models. For a source that reports what it read
    (``datasets_used(projectile)``, ``content_sha256_for(projectile)``) and that the selection
    maps to registered pins, they are derived from the source and ``materials`` (the materials
    of the tables built); an explicit argument wins. A source with ``window_per_u`` adds the
    resolved per-species window (MeV per nucleon) of the tables as ``domain_resolved``.
    """
    from ionmc.physics.stopping_evidence import STOPPING_CIRCULAR_ROWS  # noqa: PLC0415
    from ionmc.transport.tables import RANGE_CONSTRUCTION  # noqa: PLC0415

    eid = sel.effective_id(projectile)
    spec = sel.effective_spec(projectile)
    ls_requested = str(sel.spec.options.get("ls", "off"))
    ls_effective = ls_requested if projectile.z >= 2 else "off"
    digest = spec_sha256(spec)
    options: dict[str, Any] = {
        **_canon(sel.spec.options),
        "class": sel.source_class,
        "ls_requested": ls_requested,
        "ls_effective": ls_effective,
    }
    domain = _canon(spec.domain)
    source = sel.source
    used_ids = getattr(source, "datasets_used", None)
    if not data_used and callable(used_ids) and spec.data:
        data_used = [
            {
                "pin": _canon(p),
                "role_in_this_run": "construction",
                "table": [projectile.name, m.name],
            }
            for p in spec.data
            if p.dataset_id in used_ids(projectile)
            for m in materials
        ]
    digest_for = getattr(source, "content_sha256_for", None)
    if content_sha256 is None and callable(digest_for) and spec.data:
        content_sha256 = digest_for(projectile)
    fields: dict[str, Any] = {
        "version": spec.version,
        "kind": spec.kind,
        "interpolation": spec.interpolation,
        "range_start": spec.range_start,
        "data_declared": [_canon(p) for p in spec.data],
        "decision": spec.decision,
    }
    block: dict[str, Any] = {
        "model_id_requested": sel.requested_id,
        "model_id_effective": eid,
        "version": fields["version"],
        "kind": fields["kind"],
        "spec_sha256": digest,
        "options": options,
        "interpolation": fields["interpolation"],
        "range_construction": range_construction or RANGE_CONSTRUCTION,
        "range_start": fields["range_start"],
        "domain": domain,
        "data_declared": fields["data_declared"],
        "data_used": [dict(d) for d in data_used],
        "content_sha256": content_sha256,
        "evidence_at_run": _evidence(spec, digest),
        "decision": fields["decision"],
        "circular_rows": list(STOPPING_CIRCULAR_ROWS.get(digest, ())),
        "selected_by": sel.selected_by,
        "species_transported": list(species_transported),
    }
    window = getattr(source, "window_per_u", None)
    if callable(window):
        lo, hi = window(projectile)
        block["domain_resolved"] = {projectile.name: {"e_min_per_u": lo, "e_max_per_u": hi}}
    if projectile.z >= 2:
        block["mcs_proton_path_stopping"] = MCS_PROTON_PATH_STOPPING
    return block


def validate_identity_block(
    block: Mapping[str, Any], transported: Iterable[str] | None = None
) -> None:
    """Reject an inconsistent identity block (raises ``ValueError``).

    Checks the keys, that every ``data_used`` pin is one of ``data_declared``, is a construction
    use, and names a table of a projectile that was transported (``transported`` overrides the
    block's own ``species_transported``): a pin of a projectile that was not transported is
    rejected (PER-3 negative control); the requested/effective ids agree except the z = 1 LS no-op.
    """
    missing = [k for k in BLOCK_KEYS if k not in block]
    if missing:
        raise ValueError(f"stopping identity block lacks keys {missing}")
    run_species = set(block["species_transported"] if transported is None else transported)
    declared = block["data_declared"]
    for use in block["data_used"]:
        if use["pin"] not in declared:
            raise ValueError(
                f"data_used pin {use['pin'].get('dataset_id')!r} is not in data_declared"
            )
        if use.get("role_in_this_run") != "construction":
            raise ValueError("a data_used entry must have role_in_this_run == 'construction'")
        proj = use["table"][0]
        if proj not in run_species:
            raise ValueError(
                f"data_used pin {use['pin'].get('dataset_id')!r} builds a table of {proj!r}, which "
                f"was not transported in this run (transported: {sorted(run_species)})"
            )
    req, eff = block["model_id_requested"], block["model_id_effective"]
    if req != eff and not (req in LS_VARIANT_IDS and eff == DEFAULT_STOPPING_MODEL_ID):
        raise ValueError(f"model_id_requested {req!r} differs from model_id_effective {eff!r}")
