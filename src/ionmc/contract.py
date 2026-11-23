"""Capability contract (decision 0045 section 4): declared, enforced, requested, effective.

Three pieces, all pure (no I/O, no transport):

* :data:`RECONCILIATION_TABLE` and :func:`reconcile` -- the frozen requested -> effective table of
  decision 0045 section 4 item 3 as data and as the one implementation that ``Simulation``
  (after ``validate``), ``results.save`` (``effective_extra.resolutions``) and the tier-1 reader
  ``results.load`` all call, so that the three places cannot drift. Every requested field has one
  row of class ``equal``, ``derived`` (with a rule name), ``unused`` (with the reason) or
  ``requested-only``; any difference outside the table, and any requested feature without its
  effective implementation object, raises ``UnsupportedCombinationError``.
* :func:`declared_capabilities` -- the structured part of ``capabilities()``: every entry that names
  a supported or rejected combination is computed from the constants and predicates that
  ``validate`` enforces (read from ``ionmc.config`` at call time, so a patched constant changes the
  declaration and the enforcement together). Free text appears only under ``"informative"``.
* :func:`sweep` -- the finite declared matrix (species x energy-domain class x backend x precision
  x nuclear x tallies x trace x stopping model) classified against ``validate``: supported by
  ``validate``, refused before transport (``UnsupportedCombinationError`` or
  ``BackendUnavailableError``), or a third outcome, which is a failure.

Nothing here makes a claim about combinations outside the enumerated dimensions.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from ionmc import _result_schema_1_0 as frozen
from ionmc import config as _cfg
from ionmc.errors import BackendUnavailableError, UnsupportedCombinationError
from ionmc.geometry import BoxPhantom, VoxelGeometry
from ionmc.physics.projectiles import PROTON
from ionmc.physics.stopping_models import DEFAULT_STOPPING_MODEL_ID, STOPPING_MODELS
from ionmc.physics.stopping_selection import LS_VARIANT_IDS
from ionmc.sources import PencilBeamSource

EQUAL, DERIVED, UNUSED, REQUESTED_ONLY = "equal", "derived", "unused", "requested-only"
CLASSES = (EQUAL, DERIVED, UNUSED, REQUESTED_ONLY)


# ---------------------------------------------------------------------------------------------
# the frozen reconciliation table


@dataclass(frozen=True)
class Row:
    """One row of the frozen table: requested fields, effective pointers (JSON pointers into the
    effective summary; ``()`` where there is none) and the (class, rule) pairs the row can take."""

    id: str
    requested: tuple[str, ...]
    pointers: tuple[str, ...]
    outcomes: tuple[tuple[str, str], ...]
    note: str = ""


def _r(
    id: str, req: str, ptrs: str, outcomes: str, note: str = ""
) -> Row:  # compact table constructor
    return Row(
        id,
        tuple(req.split()),
        tuple(ptrs.split()),
        tuple((o.split(":")[0], o.split(":")[1]) for o in outcomes.split()),
        note,
    )


RECONCILIATION_TABLE: tuple[Row, ...] = (
    _r(
        "source",
        "source.projectile source.position_mm source.kinetic_energy_mev "
        "source.energy_sigma_mev source.lateral_sigma_mm",
        "/source/projectile /source/z /source/a /source/mass_mev /source/position_mm "
        "/source/kinetic_energy_mev /source/energy_sigma_mev /source/lateral_sigma_mm",
        "equal:equal",
    ),  # fmt: skip
    _r(
        "direction",
        "source.direction",
        "/unit_direction",
        "derived:normalise equal:equal",
        "d / |d| (sources.py); equal if already unit",
    ),  # fmt: skip
    _r(
        "geometry-box",
        "geometry",
        "/geometry/origin_mm /geometry/spacing_mm /geometry/shape "
        "/geometry/materials /geometry/z_exit_mm",
        "derived:box-to-voxel",
        "one-voxel VoxelGeometry",
    ),  # fmt: skip
    _r(
        "geometry-voxel",
        "geometry",
        "/geometry/origin_mm /geometry/spacing_mm /geometry/shape "
        "/geometry/materials /geometry/z_exit_mm",
        "equal:equal",
        "material_index and density arrays are bound by their content digests (results catalogue)",
    ),  # fmt: skip
    _r("scoring", "scoring", "/scoring", "equal:equal"),
    _r(
        "stopping-model",
        "physics.stopping_model physics.stopping",
        "/stopping/model_id_requested /stopping/model_id_effective",
        "equal:equal derived:ls-variant-z1-default derived:user-object-id",
        "exactly one of the two is required (decision 0044 section 2); absence fails at "
        "PhysicsOptions construction and is not a resolution",
    ),  # fmt: skip
    _r(
        "stopping-identity",
        "physics.nuclear tallies",
        "/stopping/species_transported /tables/materials /tables/water_row",
        "derived:deuteron-and-water-row",
        "deuteron added for nuclear runs; LET water row added when tallies are requested",
    ),  # fmt: skip
    _r(
        "ls",
        "physics.stopping.options.ls",
        "/physics/ls_requested /physics/ls_effective",
        "derived:ls-z1-off equal:equal",
    ),  # fmt: skip
    _r(
        "physics-flags",
        "physics.nuclear physics.straggling physics.multiple_scattering physics.delta_electrons",
        "/physics/nuclear /physics/straggling "
        "/physics/multiple_scattering /physics/delta_electrons",
        "equal:equal",
    ),  # fmt: skip
    _r(
        "straggling-model",
        "physics.straggling_model",
        "/physics/straggling_model",
        "equal:equal unused:straggling-off",
    ),  # fmt: skip
    _r("mcs-model", "physics.mcs_model", "/physics/mcs_model", "equal:equal unused:mcs-off"),
    _r(
        "mcs-mapping",
        "physics.mcs_mapping_eval",
        "/mcs_mapping/evaluation /mcs_mapping/unused",
        "equal:equal derived:mcs-off-root unused:proton-has-no-mapping",
    ),  # fmt: skip
    _r(
        "e-cut",
        "physics.e_cut_mev",
        "/physics/e_cut_effective_mev",
        "equal:equal unused:ion-uses-per-u",
    ),  # fmt: skip
    _r(
        "e-cut-ion",
        "physics.e_cut_ion_mev_per_u",
        "/physics/e_cut_effective_mev",
        "derived:per-u-times-A unused:proton-uses-total",
    ),  # fmt: skip
    _r(
        "nuclear-options",
        "physics.e_cut_deuteron_mev physics.nuclear_table_id",
        "/physics/e_cut_deuteron_mev /nuclear/table_id",
        "equal:equal unused:nuclear-off",
    ),  # fmt: skip
    _r(
        "end-of-range",
        "physics.end_of_range_model",
        "/physics/end_of_range_model_effective",
        "equal:equal derived:default-by-z",
    ),  # fmt: skip
    _r(
        "physics-numerics",
        "physics.max_step_mm physics.max_energy_loss_fraction "
        "physics.range_alpha physics.range_rho_f_mm physics.short_step_fraction "
        "physics.truncated_hinge_diagnostic",
        "/physics/max_step_mm "
        "/physics/max_energy_loss_fraction /physics/range_alpha /physics/range_rho_f_mm "
        "/physics/short_step_fraction /physics/truncated_hinge_diagnostic",
        "equal:equal",
    ),  # fmt: skip
    _r(
        "run",
        "run.backend run.precision run.seed run.n_histories run.n_batches run.cpu_workers",
        "/backend /precision /seed /n_histories /n_batches /cpu_workers /production",
        "equal:equal derived:float32-warp",
    ),  # fmt: skip
    _r(
        "max-steps", "run.max_steps", "/max_steps /max_steps_origin", "equal:user derived:computed"
    ),  # fmt: skip
    _r("chunk", "run.chunk_histories", "/chunk_histories", "equal:equal unused:python-backend"),
    _r(
        "run-only",
        "run.allow_invalid_result run.worker_timeout_s run.memory_budget_bytes",
        "",
        "requested-only:execution-option",
    ),  # fmt: skip
    _r(
        "diagnostics",
        "diagnostics",
        "",
        "equal:arrays-present",
        "the arrays are checked by the honoured check of save and load",
    ),  # fmt: skip
    _r("tallies", "tallies", "/tallies/quantities", "equal:in-channel-plan"),
    _r("lookups", "lookups", "/tallies/lookups", "equal:in-channel-plan"),
)
ROWS: Mapping[str, Row] = {r.id: r for r in RECONCILIATION_TABLE}


@dataclass(frozen=True)
class Resolution:
    """The recorded outcome of one table row for one configuration."""

    row: str
    cls: str
    rule: str
    requested: tuple[str, ...]
    pointers: tuple[str, ...]
    reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "row": self.row,
            "class": self.cls,
            "rule": self.rule,
            "requested": list(self.requested),
            "effective": list(self.pointers),
            "reason": self.reason,
        }


_ABSENT = object()


def pointer_get(doc: Any, pointer: str) -> Any:
    """Value at a JSON pointer ('/a/b/0') of nested dicts and lists, or ``_ABSENT``."""
    cur = doc
    for part in pointer.split("/")[1:]:
        try:
            cur = cur[int(part)] if isinstance(cur, list) else cur[part]
        except (KeyError, IndexError, TypeError, ValueError):
            return _ABSENT
    return cur


def pointer_resolves(effective: Any, pointer: str) -> bool:
    return pointer_get(effective, pointer) is not _ABSENT


def _differ(what: str, requested: Any, effective: Any) -> UnsupportedCombinationError:
    shown = "<absent>" if effective is _ABSENT else repr(effective)
    return UnsupportedCombinationError(f"{what}: requested {requested!r}, effective {shown}")


class _Ctx:
    def __init__(self, rq: Mapping[str, Any], eff: Mapping[str, Any]) -> None:
        self.rq, self.eff = rq, eff
        self.out: list[Resolution] = []

    def need(self, what: str, requested: Any, pointer: str) -> Any:
        got = pointer_get(self.eff, pointer)
        if got is _ABSENT:
            raise _differ(what, requested, got)
        return got

    def equal(self, what: str, requested: Any, pointer: str) -> None:
        got = self.need(what, requested, pointer)
        if got != requested:
            raise _differ(what, requested, got)

    def add(self, row: str, cls: str, rule: str, pointers: Sequence[str], reason: str = "") -> None:
        r = ROWS[row]
        if (cls, rule) not in r.outcomes:
            raise AssertionError(f"row {row!r} has no outcome {cls}:{rule}")  # a table bug
        self.out.append(Resolution(row, cls, rule, r.requested, tuple(pointers), reason))


def _projectile(rp: Any) -> dict[str, Any]:
    if isinstance(rp, str):
        z, a, mass = frozen.PROJECTILES[rp]
        return {"projectile": rp, "z": z, "a": a, "mass_mev": mass}
    return {k: rp[{"projectile": "name"}.get(k, k)] for k in ("projectile", "z", "a", "mass_mev")}


def _stopping_ids(rph: Mapping[str, Any]) -> tuple[str, bool]:
    if isinstance(rph.get("stopping_model"), str):
        return rph["stopping_model"], False
    obj = (rph.get("stopping") or {}).get("stopping_object")
    if isinstance(obj, Mapping) and isinstance(obj.get("model_id"), str):
        return obj["model_id"], True
    raise UnsupportedCombinationError(
        "requested stopping: neither a registered stopping_model nor a stopping object recorded"
    )


def _norm(v: Any) -> Any:
    return [_norm(x) for x in v] if isinstance(v, list | tuple) else v


def _row_source(c: _Ctx) -> None:
    rs = c.rq["source"]
    for k, v in _projectile(rs["projectile"]).items():
        c.equal(f"source.{k}", v, f"/source/{k}")
    for k in ("position_mm", "kinetic_energy_mev", "energy_sigma_mev", "lateral_sigma_mm"):
        c.equal(f"source.{k}", _norm(rs[k]), f"/source/{k}")
    c.add("source", EQUAL, "equal", ROWS["source"].pointers)
    d = [float(x) for x in rs["direction"]]
    nrm = sum(x * x for x in d) ** 0.5
    if not nrm > 0:
        raise _differ("source.direction", rs["direction"], "a zero vector")
    got = c.need("source.direction", rs["direction"], "/unit_direction")
    if not (len(got) == 3 and all(abs(a - b / nrm) <= 1e-9 for a, b in zip(got, d, strict=True))):
        raise _differ("source.direction", rs["direction"], got)
    unit = abs(nrm - 1.0) <= 1e-12
    c.add("direction", EQUAL if unit else DERIVED, "equal" if unit else "normalise",
          ["/unit_direction"])  # fmt: skip


def _row_geometry(c: _Ctx) -> None:
    rg = c.rq["geometry"]
    geo = c.need("geometry", rg, "/geometry")
    if rg["type"] == "voxel":
        want = {k: _norm(rg[k]) for k in ("origin_mm", "spacing_mm", "shape", "materials",
                                          "z_exit_mm")}  # fmt: skip
        row, cls, rule = "geometry-voxel", EQUAL, "equal"
    elif rg["type"] == "box":
        want = {"origin_mm": _norm(rg["lower_mm"]), "spacing_mm": _norm(rg["size_mm"]),
                "shape": [1, 1, 1], "materials": [rg["material"]], "z_exit_mm": None}  # fmt: skip
        row, cls, rule = "geometry-box", DERIVED, "box-to-voxel"
    else:
        raise _differ("geometry.type", rg["type"], "not a supported geometry type")
    got = {k: _norm(geo.get(k, _ABSENT)) for k in want}
    if got != want:
        raise _differ("geometry", want, got)
    c.add(row, cls, rule, ROWS[row].pointers)
    scoring = c.need("scoring", c.rq["scoring"], "/scoring")
    rs = c.rq["scoring"]
    if len(scoring) != len(rs):
        raise _differ("scoring", len(rs), len(scoring))
    for i, (g, w) in enumerate(zip(scoring, rs, strict=True)):
        for k in ("name", "origin_mm", "spacing_mm", "shape"):
            if _norm(g.get(k, _ABSENT)) != _norm(w[k]):
                raise _differ(f"scoring[{i}].{k}", w[k], g.get(k, _ABSENT))
    c.add("scoring", EQUAL, "equal", ["/scoring"])


def _row_stopping(c: _Ctx, z: int) -> None:
    rph, eff = c.rq["physics"], c.eff
    rid, is_object = _stopping_ids(rph)
    c.equal("stopping.model_id_requested", rid, "/stopping/model_id_requested")
    got = c.need("stopping.model_id_effective", rid, "/stopping/model_id_effective")
    if is_object:
        want, cls, rule = rid, DERIVED, "user-object-id"
    elif z == 1 and rid in LS_VARIANT_IDS:
        want, cls, rule = DEFAULT_STOPPING_MODEL_ID, DERIVED, "ls-variant-z1-default"
    else:
        want, cls, rule = rid, EQUAL, "equal"
    if got != want:
        raise _differ("stopping model", want, got)
    c.add("stopping-model", cls, rule, ["/stopping/model_id_requested",
                                        "/stopping/model_id_effective"])  # fmt: skip
    nuclear = bool(rph["nuclear"])
    species = [c.rq["source"]["projectile"] if isinstance(c.rq["source"]["projectile"], str)
               else c.rq["source"]["projectile"]["name"],
               *([frozen.DEUTERON_NAME] if nuclear else [])]  # fmt: skip
    c.equal("stopping.species_transported (deuteron added for nuclear runs)", species,
            "/stopping/species_transported")  # fmt: skip
    has_water = pointer_resolves(eff, "/tables/water_row")
    if has_water != bool(c.rq.get("tallies")):
        raise _differ("tallies (LET water row)", bool(c.rq.get("tallies")), has_water)
    c.add("stopping-identity", DERIVED, "deuteron-and-water-row",
          ROWS["stopping-identity"].pointers[:2]
          + (("/tables/water_row",) if has_water else ()))  # fmt: skip
    ls_req = c.need("physics.ls_requested", rid, "/physics/ls_requested")
    ls_eff = c.need("physics.ls_effective", rid, "/physics/ls_effective")
    if not is_object:
        want_ls = {"bethe-ionmc-v1+ls-eff": "eff", "bethe-ionmc-v1+ls-bare_mixed": "bare_mixed"}
        if ls_req != want_ls.get(rid, "off"):
            raise _differ("stopping option ls", want_ls.get(rid, "off"), ls_req)
    if ls_eff != (ls_req if z >= 2 else "off"):
        raise _differ("stopping option ls (rule ls-z1-off)", ls_req if z >= 2 else "off", ls_eff)
    c.add("ls", DERIVED if z == 1 else EQUAL, "ls-z1-off" if z == 1 else "equal",
          ROWS["ls"].pointers)  # fmt: skip


def _row_physics(c: _Ctx, z: int) -> None:
    rph, eff = c.rq["physics"], c.eff
    ion = z >= 2
    for k in ("nuclear", "straggling", "multiple_scattering", "delta_electrons"):
        c.equal(f"physics.{k}", rph[k], f"/physics/{k}")
    c.add("physics-flags", EQUAL, "equal", ROWS["physics-flags"].pointers)
    modes = (
        ("straggling", "straggling-model", "straggling_model", "straggling-off"),
        ("multiple_scattering", "mcs-model", "mcs_model", "mcs-off"),
    )
    for flag, row, key, off in modes:
        c.equal(f"physics.{key}", rph[key], f"/physics/{key}")
        on = bool(rph[flag])
        c.add(row, EQUAL if on else UNUSED, "equal" if on else off, [f"/physics/{key}"])
    mm = c.need("mcs_mapping", rph["mcs_mapping_eval"], "/mcs_mapping")
    if not ion:
        c.add("mcs-mapping", UNUSED, "proton-has-no-mapping", [])
    elif rph["multiple_scattering"]:
        c.equal("physics.mcs_mapping_eval", rph["mcs_mapping_eval"], "/mcs_mapping/evaluation")
        c.add("mcs-mapping", EQUAL, "equal", ["/mcs_mapping/evaluation"])
    else:
        if not (mm.get("evaluation") == "root" and mm.get("unused") is True):
            raise _differ("physics.mcs_mapping_eval (rule mcs-off-root)", "root, unused", mm)
        c.add("mcs-mapping", DERIVED, "mcs-off-root", ROWS["mcs-mapping"].pointers)
    a = c.need("source.a", rph, "/source/a")
    if ion:
        if rph["e_cut_mev"] != _cfg.PhysicsOptions.__dataclass_fields__["e_cut_mev"].default:
            raise _differ("physics.e_cut_mev (ion source)", rph["e_cut_mev"], "unused: default 2.0")
        c.equal("physics.e_cut_ion_mev_per_u x A", rph["e_cut_ion_mev_per_u"] * a,
                "/physics/e_cut_effective_mev")  # fmt: skip
        c.add("e-cut", UNUSED, "ion-uses-per-u", [])
        c.add("e-cut-ion", DERIVED, "per-u-times-A", ["/physics/e_cut_effective_mev"])
    else:
        c.equal("physics.e_cut_mev", rph["e_cut_mev"], "/physics/e_cut_effective_mev")
        c.add("e-cut", EQUAL, "equal", ["/physics/e_cut_effective_mev"])
        c.add("e-cut-ion", UNUSED, "proton-uses-total", [])
    if rph["nuclear"]:
        c.equal("physics.e_cut_deuteron_mev", rph["e_cut_deuteron_mev"],
                "/physics/e_cut_deuteron_mev")  # fmt: skip
        ptrs = ["/physics/e_cut_deuteron_mev"]
        if rph["nuclear_table_id"] is not None:
            c.equal("physics.nuclear_table_id", rph["nuclear_table_id"], "/nuclear/table_id")
            ptrs.append("/nuclear/table_id")
        elif not pointer_resolves(eff, "/nuclear"):
            raise _differ("physics.nuclear", True, _ABSENT)
        c.add("nuclear-options", EQUAL, "equal", ptrs)
    else:
        if pointer_resolves(eff, "/nuclear"):
            raise _differ("physics.nuclear", False, "a nuclear setup")
        c.add("nuclear-options", UNUSED, "nuclear-off", [])
    eor = rph["end_of_range_model"]
    want = ("continuous_residual_v1" if ion else "point_v1") if eor == "default" else eor
    c.equal("physics.end_of_range_model", want, "/physics/end_of_range_model_effective")
    c.add("end-of-range", DERIVED if eor == "default" else EQUAL,
          "default-by-z" if eor == "default" else "equal",
          ["/physics/end_of_range_model_effective"])  # fmt: skip
    for k in ("max_step_mm", "max_energy_loss_fraction", "range_alpha", "range_rho_f_mm",
              "short_step_fraction", "truncated_hinge_diagnostic"):  # fmt: skip
        c.equal(f"physics.{k}", rph[k], f"/physics/{k}")
    c.add("physics-numerics", EQUAL, "equal", ROWS["physics-numerics"].pointers)
    # the implementation objects of the requested features (a drop is never silent)
    if rph["multiple_scattering"] and not (
        pointer_resolves(eff, "/scattering_length_g_cm2") and pointer_resolves(eff, "/mcs_mapping")
    ):
        raise _differ("physics.multiple_scattering", True, "no MCS model")
    if rph["straggling"] and not isinstance(pointer_get(eff, "/physics/straggling_model"), str):
        raise _differ("physics.straggling", True, "no straggling model")


def _row_run(c: _Ctx, z: int) -> None:
    rr = c.rq["run"]
    for k in ("backend", "precision", "seed", "n_histories", "n_batches", "cpu_workers"):
        c.equal(f"run.{k}", rr[k], f"/{k}")
    c.equal("run.production (rule float32-warp)",
            rr["backend"] != "python" and rr["precision"] == "float32", "/production")  # fmt: skip
    prod = bool(c.eff["production"])
    c.add(
        "run", DERIVED if prod else EQUAL, "float32-warp" if prod else "equal", ROWS["run"].pointers
    )
    user = rr["max_steps"] is not None
    c.equal("run.max_steps_origin", "user" if user else "computed", "/max_steps_origin")
    if user:
        c.equal("run.max_steps", rr["max_steps"], "/max_steps")
    c.add("max-steps", EQUAL if user else DERIVED, "user" if user else "computed",
          ["/max_steps", "/max_steps_origin"])  # fmt: skip
    c.equal("run.chunk_histories", rr["chunk_histories"], "/chunk_histories")
    py = rr["backend"] == "python"
    c.add(
        "chunk", UNUSED if py else EQUAL, "python-backend" if py else "equal", ["/chunk_histories"]
    )
    c.add("run-only", REQUESTED_ONLY, "execution-option", [])


def _row_requests(c: _Ctx) -> None:
    c.add("diagnostics", EQUAL, "arrays-present", [])
    tallies, lookups = c.rq.get("tallies") or [], c.rq.get("lookups") or []
    plan = c.eff.get("tallies")
    if tallies or lookups:
        if not isinstance(plan, Mapping):
            raise _differ("tallies/lookups", f"{len(tallies)} tallies, {len(lookups)} lookups",
                          "no channel plan")  # fmt: skip
    plan_map: Mapping[str, Any] = plan if isinstance(plan, Mapping) else {}
    if tallies:
        names = {q.get("name") for q in plan_map["quantities"]}
        want = {t["name"] for t in tallies}
        if names != want:
            raise _differ("tallies", sorted(want), sorted(map(str, names)))
        c.add("tallies", EQUAL, "in-channel-plan", ["/tallies/quantities"])
    else:
        if plan is not None:
            raise _differ("tallies", [], "a channel plan")
        c.add("tallies", EQUAL, "in-channel-plan", [])
    if len(plan_map.get("lookups", []) if lookups else []) != len(lookups):
        raise _differ("lookups", len(lookups), len(plan_map.get("lookups", [])) if plan else 0)
    c.add("lookups", EQUAL, "in-channel-plan", ["/tallies/lookups"] if lookups else [])


def reconcile(requested: Any, effective: Any) -> tuple[Resolution, ...]:
    """Requested -> effective reconciliation (decision 0045 section 4 item 3), pure.

    ``requested`` is a ``SimulationConfig`` or its ``configio.to_dict`` form; ``effective`` an
    ``EffectiveConfig`` or its ``summary()`` dict. Returns one :class:`Resolution` per table row
    and raises ``UnsupportedCombinationError("... requested X, effective Y")`` for any difference
    outside the table and for a requested feature without its effective implementation object."""
    if isinstance(requested, _cfg.SimulationConfig):
        from ionmc.configio import to_dict

        requested = to_dict(requested)
    if isinstance(effective, _cfg.EffectiveConfig):
        effective = effective.summary()
    c = _Ctx(requested, effective)
    try:
        z = int(pointer_get(effective, "/source/z"))
    except (TypeError, ValueError) as exc:
        raise _differ("source.z", requested["source"]["projectile"], _ABSENT) from exc
    _row_source(c)
    _row_geometry(c)
    _row_stopping(c, z)
    _row_physics(c, z)
    _row_run(c, z)
    _row_requests(c)
    order = {r.id: i for i, r in enumerate(RECONCILIATION_TABLE)}
    return tuple(sorted(c.out, key=lambda r: order[r.row]))


# ---------------------------------------------------------------------------------------------
# declared capabilities (derived from the constants and predicates that validate enforces)

SPECIES = (PROTON, *_cfg.ION_PROJECTILES)
GEOMETRY_TYPES = (VoxelGeometry, BoxPhantom)
SOURCE_TYPES = (PencilBeamSource,)
DOMAIN_CLASSES = ("below", "inside", "upper-edge", "above")
OK, UNSUPPORTED, UNAVAILABLE = "supported", "UnsupportedCombinationError", "BackendUnavailableError"


def e_cut_total_mev(name: str, a: int) -> float:
    """Effective cutoff [MeV] of the species at default options (``effective_e_cut_mev``)."""
    if name == PROTON.name:
        default = _cfg.PhysicsOptions.__dataclass_fields__["e_cut_mev"].default
        if not isinstance(default, int | float):  # pragma: no cover - a dataclass invariant
            raise TypeError("PhysicsOptions.e_cut_mev has no numeric default")
        return float(default)
    return float(_cfg.E_CUT_ION_DEFAULT_MEV_PER_U * a)


def source_window_total_mev(
    species: Any, model_id: str, nuclear: bool = False
) -> tuple[float, float] | None:
    """Supported source energy window [MeV total, inclusive; E0 + 6 sigma_E must stay inside] of a
    species with a stopping model at default cutoffs; ``None`` when the model's domain excludes the
    species or the combination is refused whatever the energy (an ion with nuclear)."""
    dom = STOPPING_MODELS[model_id].domain
    if species.name not in dom.species or (nuclear and species.z >= 2):
        return None
    if nuclear and "deuteron" not in dom.species:  # secondary deuterons need a table
        return None
    lo = 2.0 * e_cut_total_mev(species.name, species.a)
    hi = float(dom.e_max_per_u[species.name]) * species.a
    if species.z >= 2:
        lo_u, hi_u = _cfg.ION_DOMAIN_MEV_PER_U[species.name]
        lo, hi = max(lo, lo_u * species.a), min(hi, hi_u * species.a)
    elif nuclear:
        hi = min(hi, float(_cfg.NUCLEAR_MAX_ENERGY_MEV))
    return lo, hi


def energy_for_class(window: tuple[float, float], domain_class: str) -> float:
    lo, hi = window
    return {"below": 0.5 * lo, "inside": 0.5 * (lo + hi), "upper-edge": hi, "above": 1.01 * hi}[
        domain_class
    ]


def run_combination_status(
    backend: str, precision: str, nuclear: bool, tallies: bool, trace: bool, workers: int,
    species: Any = PROTON,
) -> tuple[str, str]:  # fmt: skip
    """(status, rule) of a run combination at an in-domain energy: the predicates of ``validate``
    over backend x precision x {nuclear, tallies, trace, workers}."""
    if backend not in _cfg.BACKENDS or precision not in _cfg.PRECISIONS:
        return UNSUPPORTED, "unknown-backend-or-precision"
    if nuclear and species.z >= 2:
        return UNSUPPORTED, "nuclear-ion"
    if nuclear and backend != "python":
        return UNSUPPORTED, "nuclear-python-only"
    if backend == "python" and precision == "float32":
        return UNSUPPORTED, "python-float64-only"
    if backend == "warp-cuda" and workers > 1:
        return UNSUPPORTED, "cuda-single-worker"
    if workers > _cfg.MAX_CPU_WORKERS:
        return UNSUPPORTED, "max-cpu-workers"
    if trace and precision != "float64":
        return UNSUPPORTED, "trace-float64-only"
    if tallies and backend not in _cfg.CHANNEL_BACKENDS:
        return UNSUPPORTED, "tallies-need-channel-backend"
    if backend == "warp-cuda" and not _cfg.cuda_available():
        return UNAVAILABLE, "cuda-device"
    return OK, "supported"


def declared_capabilities() -> dict[str, Any]:
    """Structured declaration of what ``validate`` enforces (a section of ``capabilities()``)."""
    models = {}
    for mid in STOPPING_MODELS:
        per: dict[str, Any] = {}
        for sp in SPECIES:
            w = source_window_total_mev(sp, mid)
            wn = source_window_total_mev(sp, mid, nuclear=True)
            per[sp.name] = {
                "window_total_mev": None if w is None else {"min": w[0], "max": w[1]},
                "window_total_mev_nuclear": None if wn is None else {"min": wn[0], "max": wn[1]},
            }
        models[mid] = per
    combos = []
    for b, p, n, t, tr, wk in itertools.product(
        _cfg.BACKENDS, _cfg.PRECISIONS, (False, True), (False, True), (False, True), (1, 2)
    ):
        status, rule = run_combination_status(b, p, n, t, tr, wk)
        combos.append({"backend": b, "precision": p, "nuclear": n, "tallies": t, "trace": tr,
                       "workers_above_one": wk > 1, "status": status, "rule": rule})  # fmt: skip
    return {
        "species": [s.name for s in SPECIES],
        "source_window_by_model": models,
        "source_energy_rules": {
            "lower_total_mev": "twice-effective-cutoff",
            "sigma_fraction_max": _cfg.MAX_ENERGY_SIGMA_FRACTION,
            "nuclear_max_total_mev": _cfg.NUCLEAR_MAX_ENERGY_MEV,
            "ion_domain_mev_per_u": {k: list(v) for k, v in _cfg.ION_DOMAIN_MEV_PER_U.items()},
            "e_cut_mev_default": e_cut_total_mev(PROTON.name, 1),
            "e_cut_ion_mev_per_u_default": _cfg.E_CUT_ION_DEFAULT_MEV_PER_U,
        },
        "run_combinations": combos,
        "limits": {
            "max_cpu_workers": _cfg.MAX_CPU_WORKERS,
            "min_chunk_histories": _cfg.MIN_CHUNK_HISTORIES,
            "max_scoring_grids": _cfg.MAX_SCORING_GRIDS,
            "max_scoring_pieces": _cfg.MAX_SCORING_PIECES,
            "max_trace_buffer_bytes": _cfg.MAX_TRACE_BUFFER_BYTES,
            "channel_backends": list(_cfg.CHANNEL_BACKENDS),
        },
        "geometry_types": [t.__name__ for t in GEOMETRY_TYPES],
        "source_types": [t.__name__ for t in SOURCE_TYPES],
        "informative": {
            "scope": "the combinations above are derived from the constants and predicates of "
            "validate; no claim is made outside the enumerated dimensions",
            "tabulated_models": "swept against validate only where their cache is present",
        },
    }


# ---------------------------------------------------------------------------------------------
# the sweep: declared == enforced over the enumerated dimensions


@dataclass(frozen=True)
class Matrix:
    species: tuple[Any, ...] = SPECIES
    domain_classes: tuple[str, ...] = DOMAIN_CLASSES
    backends: tuple[str, ...] = tuple(_cfg.BACKENDS)
    precisions: tuple[str, ...] = tuple(_cfg.PRECISIONS)
    nuclear: tuple[bool, ...] = (False, True)
    tallies: tuple[bool, ...] = (False, True)
    trace: tuple[bool, ...] = (False, True)
    models: tuple[str, ...] = ("bethe-ionmc-v1", "bethe-ionmc-v1+ls-eff",
                               "bethe-ionmc-v1+ls-bare_mixed")  # fmt: skip
    workers: tuple[int, ...] = (1,)
    materials: tuple[Any, ...] = ()  # default: water
    nuclear_table_id: str | None = None

    def size(self) -> int:
        return (len(self.species) * len(self.domain_classes) * len(self.backends)
                * len(self.precisions) * len(self.nuclear) * len(self.tallies) * len(self.trace)
                * len(self.models) * len(self.workers) * max(1, len(self.materials)))  # fmt: skip


@dataclass(frozen=True)
class Outcome:
    """One swept case: what the declaration says, what ``validate`` did, and the verdict
    (``agree``, ``disagree``, ``third`` = neither success nor a refusal before transport, or
    ``not-evaluated`` = a declared-supported case that needs a cache object this tier lacks)."""

    species: str
    domain_class: str
    energy_mev: float
    backend: str
    precision: str
    nuclear: bool
    tallies: bool
    trace: bool
    model: str
    workers: int
    material: str
    declared: str
    rule: str
    actual: str
    verdict: str


def _declare(case: Mapping[str, Any]) -> tuple[str, str]:
    sp, model = case["species"], case["model"]
    window = source_window_total_mev(sp, model, case["nuclear"])
    if window is None:
        return UNSUPPORTED, "species-or-nuclear-not-in-model-domain"
    if sp.z >= 2 and case["nuclear"]:
        return UNSUPPORTED, "nuclear-ion"
    mat = case["material"]
    if (STOPPING_MODELS[model].kind == "tabulated" and mat is not None
            and mat.name != "water"):  # fmt: skip
        return UNSUPPORTED, "material-not-in-model-domain"
    status, rule = run_combination_status(case["backend"], case["precision"], case["nuclear"],
                                          case["tallies"], case["trace"], case["workers"],
                                          sp)  # fmt: skip
    if status == UNSUPPORTED:
        return status, rule
    lo, hi = window
    e = case["energy_mev"]
    if not lo <= e <= hi:
        return UNSUPPORTED, "source-energy-domain"
    return status, rule


def build_config(case: Mapping[str, Any]) -> Any:
    """The tiny ``SimulationConfig`` of one swept case (never run, only validated)."""
    from ionmc.materials import WATER
    from ionmc.scoring import ScoringGrid, TallyRequest

    sp = case["species"]
    mat = case["material"]
    return _cfg.SimulationConfig(
        PencilBeamSource(sp, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), case["energy_mev"]),
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), mat if mat is not None else WATER),
        (ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 5.0), (4, 4, 4)),),
        _cfg.PhysicsOptions(
            nuclear=case["nuclear"],
            nuclear_table_id=case["nuclear_table_id"] if case["nuclear"] else None,
            stopping_model=case["model"],
            multiple_scattering=sp.z < 2,  # the ion f_dM mapping domain is outside the dimensions
        ),
        _cfg.RunOptions(
            case["backend"], case["precision"], 1, 40, n_batches=4, cpu_workers=case["workers"]
        ),  # fmt: skip
        _cfg.DiagnosticsOptions(trace_histories=1 if case["trace"] else 0),
        tallies=(TallyRequest("let", "dose", "let_d"),) if case["tallies"] else (),
    )


def sweep(matrix: Matrix | None = None) -> list[Outcome]:
    """Classify every case of the matrix against ``validate`` (no transport)."""
    m = matrix or Matrix()
    mats = m.materials or (None,)
    out: list[Outcome] = []
    for sp, model in itertools.product(m.species, m.models):
        window = source_window_total_mev(sp, model) or source_window_total_mev(sp, model, True)
        for dc, b, p, n, t, tr, w, mat in itertools.product(
            m.domain_classes,
            m.backends,
            m.precisions,
            m.nuclear,
            m.tallies,
            m.trace,
            m.workers,
            mats,
        ):
            win = source_window_total_mev(sp, model, n)
            e = energy_for_class(win if win is not None else (window or (4.0, 100.0)), dc)
            case = {"species": sp, "model": model, "energy_mev": e, "backend": b, "precision": p,
                    "nuclear": n, "tallies": t, "trace": tr, "workers": w, "material": mat,
                    "nuclear_table_id": m.nuclear_table_id}  # fmt: skip
            declared, rule = _declare(case)
            needs_table = n and declared == OK and m.nuclear_table_id is None
            try:
                if needs_table:
                    actual, verdict = "not-evaluated", "not-evaluated"
                else:
                    _cfg.validate(build_config(case))
                    actual = OK
                    verdict = "agree" if declared == OK else "disagree"
            except (UnsupportedCombinationError, BackendUnavailableError) as exc:
                actual = type(exc).__name__
                verdict = "agree" if actual == declared else "disagree"
            except Exception as exc:  # noqa: BLE001 - a third outcome is the failure being hunted
                actual, verdict = f"third:{type(exc).__name__}", "third"
            out.append(Outcome(sp.name, dc, e, b, p, n, t, tr, model, w,
                               "water" if mat is None else mat.name, declared, rule, actual,
                               verdict))  # fmt: skip
    return out


def summarise(outcomes: Iterable[Outcome]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for o in outcomes:
        counts[o.verdict] = counts.get(o.verdict, 0) + 1
    return counts
