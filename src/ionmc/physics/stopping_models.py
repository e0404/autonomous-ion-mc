"""Registry of stopping-power models (decision 0044, section 1; V3-013 step S1).

A model is described by an immutable :class:`StoppingModelSpec` (no evidence status: that lives in
the append-only registry of :mod:`ionmc.physics.stopping_evidence`, keyed by :func:`spec_sha256`).
:func:`resolve` maps a registered ``model_id`` to the existing source class and never returns a
source of another id; :func:`spec_of` identifies a source *instance* from its fields, so a modified
source is never reported under a registered id.

Registered in step S1: the default ``bethe-ionmc-v1`` and its two Lindhard-Sorensen sensitivity
variants; in step S3 ``nist-pstar-astar-water-2005`` (its source is built in
:mod:`ionmc.physics.stopping_nist` from the cached PSTAR/ASTAR objects). ``icru90-water-v1`` (S4)
is registered by its step.

The default is guarded like ``A16_INTENDED_CHANGE``: the default id and its ``spec_sha256`` are
pinned by a test, and a change requires ``STOPPING_DEFAULT_INTENDED_CHANGE`` to name a decision.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, fields, is_dataclass, replace
from types import MappingProxyType
from typing import Any

from ionmc.data.registry import DATASETS
from ionmc.physics.stopping import (
    BetheOptions,
    BetheStoppingSource,
    NistStarStoppingSource,
    StoppingSource,
)

DEFAULT_STOPPING_MODEL_ID = "bethe-ionmc-v1"

# None in V3-013. A future change of the default id or of its spec_sha256 must replace this by a
# record that names a decision (decision 0044 section 5), mirroring A16_INTENDED_CHANGE.
STOPPING_DEFAULT_INTENDED_CHANGE: Mapping[str, str] | None = None

UNREGISTERED_ID = "unregistered"
USE_ROLES = ("construction-permitted", "evaluation-only", "lineage-only")
KINDS = ("analytic", "tabulated")

_BETHE_SPECIES = ("proton", "deuteron", "alpha", "carbon12", "oxygen16")
MATERIAL_RULE_ANY = "any ionmc.materials material (evidence: liquid water only)"
MATERIAL_RULE_WATER = "structural liquid water (composition and density, not the name)"
_BETHE_MATERIAL_RULE = MATERIAL_RULE_ANY
_BETHE_INTERPOLATION = "analytic on a log grid; range exact-loglog-quadrature-v1"
_BETHE_RANGE_START = "constant-S approximation a*E_min/S(E_min)"


def _freeze(mapping: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(dict(mapping))


@dataclass(frozen=True)
class StoppingDomain:
    """Species, material rule and per-species energy window (MeV per nucleon) of a model."""

    species: frozenset[str]
    material_rule: str
    e_min_per_u: Mapping[str, float]
    e_max_per_u: Mapping[str, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "species", frozenset(self.species))
        object.__setattr__(self, "e_min_per_u", _freeze(self.e_min_per_u))
        object.__setattr__(self, "e_max_per_u", _freeze(self.e_max_per_u))


@dataclass(frozen=True)
class DataPin:
    """A construction-data member of a tabulated model (none in S1)."""

    dataset_id: str
    archive_sha256: str
    member: str | None
    member_sha256: str | None
    parser: str
    licence: str
    use_role: str

    def __post_init__(self) -> None:
        if self.use_role not in USE_ROLES:
            raise ValueError(f"DataPin.use_role must be one of {USE_ROLES}, got {self.use_role!r}")


@dataclass(frozen=True)
class StoppingModelSpec:
    """Immutable description of a stopping model; carries no evidence status."""

    model_id: str
    version: str
    kind: str
    options: Mapping[str, Any]
    interpolation: str
    range_start: str
    data: tuple[DataPin, ...]
    domain: StoppingDomain
    decision: str

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"StoppingModelSpec.kind must be one of {KINDS}, got {self.kind!r}")
        object.__setattr__(self, "options", _freeze(self.options))
        object.__setattr__(self, "data", tuple(self.data))


def _canonical(value: Any) -> Any:
    """JSON-ready structure: dataclasses and mappings to dicts, sets to sorted lists."""
    if is_dataclass(value) and not isinstance(value, type):
        return {f.name: _canonical(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, Mapping):
        return {str(k): _canonical(v) for k, v in value.items()}
    if isinstance(value, frozenset | set):
        return sorted(_canonical(v) for v in value)
    if isinstance(value, list | tuple):
        return [_canonical(v) for v in value]
    return value


def spec_sha256(spec: StoppingModelSpec) -> str:
    """SHA-256 of the canonical JSON of ``spec`` (sorted keys, no whitespace, float ``repr``)."""
    text = json.dumps(
        _canonical(spec), sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    )
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def _bethe_options(source: BetheStoppingSource) -> dict[str, Any]:
    o = source.options
    return {
        "shell": bool(o.shell),
        "barkas": bool(o.barkas),
        "bloch": bool(o.bloch),
        "mott": bool(o.mott),
        "effective_charge": bool(o.effective_charge),
        "effective_charge_model": "barkas-125",
        "ls": str(o.ls),
        "points_per_decade": int(source.points_per_decade),
        "i_value": "material",
        "source_name": str(source.name),
    }


def _bethe_domain(source: BetheStoppingSource) -> StoppingDomain:
    lo, hi = float(source.e_min_per_u), float(source.e_max_per_u)
    return StoppingDomain(
        frozenset(_BETHE_SPECIES),
        _BETHE_MATERIAL_RULE,
        {s: lo for s in _BETHE_SPECIES},
        {s: hi for s in _BETHE_SPECIES},
    )


def _bethe_spec(
    model_id: str, version: str, source: BetheStoppingSource, decision: str
) -> StoppingModelSpec:
    return StoppingModelSpec(
        model_id,
        version,
        "analytic",
        _bethe_options(source),
        _BETHE_INTERPOLATION,
        _BETHE_RANGE_START,
        (),
        _bethe_domain(source),
        decision,
    )


def _registered_bethe(ls: str) -> BetheStoppingSource:
    return BetheStoppingSource(options=BetheOptions(ls=ls))


_LS_IDS = {
    "bethe-ionmc-v1": "off",
    "bethe-ionmc-v1+ls-eff": "eff",
    "bethe-ionmc-v1+ls-bare_mixed": "bare_mixed",
}

NIST_MODEL_ID = "nist-pstar-astar-water-2005"
PROGRAM_OF = {"proton": "PSTAR", "alpha": "ASTAR"}
DATASET_OF = {"proton": "nist-pstar-water-2005", "alpha": "nist-astar-water-2005"}
_E_MIN_PER_U = 1.0
_E_MAX_PER_U_CAP = 500.0  # the NistStarStoppingSource default; the grid ends at the table end below
_POINTS_PER_DECADE = 200
NIST_INTERPOLATION = "log-log on the NIST nodes"
NIST_RANGE_START = "NIST CSDA range at e_min (log-log interpolation)"


def nist_pins() -> tuple[DataPin, ...]:
    """The two construction pins (PSTAR, ASTAR), from the dataset registry (no table values)."""
    return tuple(
        DataPin(
            DATASET_OF[sp],
            DATASETS[DATASET_OF[sp]].sha256,
            None,
            None,
            DATASETS[DATASET_OF[sp]].parser,
            DATASETS[DATASET_OF[sp]].license,
            "construction-permitted",
        )
        for sp in ("proton", "alpha")
    )


def nist_options(source_name: str, ppd: int) -> dict[str, Any]:
    return {
        "points_per_decade": int(ppd),
        "source_name": source_name,
        "programs": dict(PROGRAM_OF),
        "grid": "e_min_per_u to min(e_max_per_u, table end); log grid",
        "i_value": "table (75 eV); requested I recorded in the table metadata",
    }


def registered_nist_spec() -> StoppingModelSpec:
    """The spec registered as ``nist-pstar-astar-water-2005`` (windows: lower end and the grid cap;
    the upper end of each program is read from its cached table at resolve time)."""
    species = frozenset(PROGRAM_OF)
    return StoppingModelSpec(
        NIST_MODEL_ID,
        "1",
        "tabulated",
        nist_options("nist-star", _POINTS_PER_DECADE),
        NIST_INTERPOLATION,
        NIST_RANGE_START,
        nist_pins(),
        StoppingDomain(
            species,
            MATERIAL_RULE_WATER,
            dict.fromkeys(species, _E_MIN_PER_U),
            dict.fromkeys(species, _E_MAX_PER_U_CAP),
        ),
        "0044",
    )


# icru90-water-v1 (S4) is added by its step.
_REGISTRY: dict[str, StoppingModelSpec] = {
    **{mid: _bethe_spec(mid, "1", _registered_bethe(ls), "0044") for mid, ls in _LS_IDS.items()},
    NIST_MODEL_ID: registered_nist_spec(),
}
STOPPING_MODELS: Mapping[str, StoppingModelSpec] = MappingProxyType(_REGISTRY)
# constructors of the registered models that are not one of the Bethe variants (empty in src; the
# test hook below is the only writer)
_FACTORIES: dict[str, Callable[[], StoppingSource]] = {}


def resolve(model_id: str, *, cache_dir: str | None = None) -> StoppingSource:
    """Return the source of the registered ``model_id``; an unknown id raises ``ValueError``.

    ``cache_dir`` is for tabulated models (the NIST model of S3, S4's ICRU 90); the analytic models
    read no data. A missing or corrupt cache object raises (``FileNotFoundError``,
    ``IntegrityError``, ``ValueError``), never a fallback.
    """
    if model_id not in STOPPING_MODELS:  # the module-level mapping, looked up at call time
        raise ValueError(
            f"unknown stopping model id {model_id!r}; registered: {sorted(STOPPING_MODELS)}"
        )
    spec = STOPPING_MODELS[model_id]
    for pin in spec.data:  # row ROLE-1: only construction-permitted data can build a source
        if pin.use_role != "construction-permitted":
            raise ValueError(
                f"stopping model {model_id!r} cannot be resolved: its data pin "
                f"{pin.dataset_id!r} has use_role {pin.use_role!r}, not 'construction-permitted'"
            )
    if model_id == NIST_MODEL_ID:  # reads the verified cache; never downloads (row DOM-6)
        from ionmc.physics.stopping_nist import load_nist_model  # noqa: PLC0415

        return load_nist_model(cache_dir)
    if model_id in _FACTORIES:
        return _FACTORIES[model_id]()
    if model_id not in _LS_IDS:
        raise ValueError(f"stopping model {model_id!r} is registered without a constructor")
    return _registered_bethe(_LS_IDS[model_id])


@contextmanager
def register_test_model(
    spec: StoppingModelSpec, factory: Callable[[], StoppingSource]
) -> Iterator[StoppingModelSpec]:
    """TEST ONLY: register ``spec`` and its constructor ``factory`` for the duration of a ``with``.

    Extends the registry (``STOPPING_MODELS``, the ``resolve`` constructor map and an
    ``unvalidated`` record in ``STOPPING_EVIDENCE``) and restores all three on exit, also on an
    exception. This is the hook of the consumer-consistency rows (plan v3-013 A.2); it is never
    used in ``src`` and a model registered through it is not part of the library. It refuses an id
    that does not start with ``"test-"``, an id already registered and a spec whose ``spec_sha256``
    equals a registered one, and a spec with data pins (a test model has none).
    """
    from ionmc.physics import stopping_evidence as ev  # noqa: PLC0415

    mid = spec.model_id
    if not mid.startswith("test-"):
        raise ValueError(f"a test model id must start with 'test-', got {mid!r}")
    if mid in _REGISTRY or mid in _FACTORIES:
        raise ValueError(f"stopping model id {mid!r} is already registered")
    if spec.data:
        raise ValueError("a test model cannot declare data pins")
    digest = spec_sha256(spec)
    if any(spec_sha256(s) == digest for s in _REGISTRY.values()):
        raise ValueError(f"spec_sha256 of {mid!r} equals a registered model's")
    _REGISTRY[mid] = spec
    _FACTORIES[mid] = factory
    ev.add_test_record(
        ev.EvidenceRecord(digest, "unvalidated", None, "test-only registration", None)
    )
    try:
        yield spec
    finally:
        _REGISTRY.pop(mid, None)
        _FACTORIES.pop(mid, None)
        ev.remove_test_record(digest)


def _nist_star_instance_spec(source: NistStarStoppingSource) -> StoppingModelSpec:
    """Spec of a directly constructed ``NistStarStoppingSource``, built from its instance fields.

    One program is not the combined model ``nist-pstar-astar-water-2005`` (that is
    :class:`~ionmc.physics.stopping_nist.NistStarCombinedSource`), so it is always reported as
    ``unregistered`` with ``data = ()`` (no registered pins, so ``data_used`` stays empty).
    """
    program = str(source.star_table.program)
    species = {"PSTAR": "proton", "ASTAR": "alpha"}.get(program)
    members = frozenset() if species is None else frozenset({species})
    lo, hi = float(source.e_min_per_u), float(source.e_max_per_u)
    return StoppingModelSpec(
        UNREGISTERED_ID,
        "1",
        "tabulated",
        {
            "points_per_decade": int(source.points_per_decade),
            "source_name": str(source.name),
            "program": program,
            "table_dataset_id": source.star_table.dataset_id,
            "i_value": "table (75 eV)",
        },
        "log-log interpolation of S onto the grid; range from the table CSDA range at e_min",
        "table CSDA range at e_min (log-log interpolation)",
        (),
        StoppingDomain(
            members,
            MATERIAL_RULE_WATER,
            {s: lo for s in members},
            {s: hi for s in members},
        ),
        "user",
    )


def _structure(spec: StoppingModelSpec) -> tuple[Any, ...]:
    """The fields that define a model apart from its name, version and decision."""
    return (
        spec.kind,
        _canonical(spec.options),
        spec.interpolation,
        spec.range_start,
        _canonical(spec.data),
        _canonical(spec.domain),
    )


def spec_of(source: StoppingSource) -> tuple[str, StoppingModelSpec]:
    """Identify ``source`` from its instance fields (decision 0044 section 2).

    An exact structural match with a registered spec (class, every option, grid bounds and
    density, start-range rule, data pins) yields ``(model_id, spec)``; any other supported source
    yields ``("unregistered", spec built from the instance)`` with ``decision = "user"``. A source
    type without support raises ``TypeError`` (fail closed). The identification is purely
    structural: the z = 1 no-op of the Lindhard-Sorensen variants is handled by the effective id
    of the selection layer, not here.
    """
    if type(source) is NistStarStoppingSource:
        return UNREGISTERED_ID, _nist_star_instance_spec(source)
    from ionmc.physics.stopping_nist import NistStarCombinedSource  # noqa: PLC0415
    from ionmc.physics.stopping_tabulated import TabulatedStoppingSource  # noqa: PLC0415

    if isinstance(source, NistStarCombinedSource):  # structural match incl. the verified pins
        candidate = source.instance_spec()
        if _structure(STOPPING_MODELS[NIST_MODEL_ID]) == _structure(candidate):
            return NIST_MODEL_ID, STOPPING_MODELS[NIST_MODEL_ID]
        return UNREGISTERED_ID, replace(candidate, model_id=UNREGISTERED_ID, decision="user")

    if isinstance(source, TabulatedStoppingSource):  # structural match, content digest included
        candidate = source.instance_spec()
        for model_id, spec in STOPPING_MODELS.items():
            if _structure(spec) == _structure(candidate):
                return model_id, spec
        return UNREGISTERED_ID, replace(candidate, model_id=UNREGISTERED_ID, decision="user")
    if type(source) is not BetheStoppingSource:
        hook = getattr(source, "stopping_spec", None)
        if callable(hook):  # the one documented hook: a source declares itself (decision "user")
            declared = hook()
            if not isinstance(declared, StoppingModelSpec):
                raise TypeError(
                    f"{type(source).__name__}.stopping_spec() must return a StoppingModelSpec, "
                    f"got {type(declared).__name__}"
                )
            if declared.data:
                raise ValueError(
                    f"{type(source).__name__}.stopping_spec() must declare data=() (a user source "
                    "has no registered data pins)"
                )
            return UNREGISTERED_ID, replace(declared, model_id=UNREGISTERED_ID, decision="user")
        raise TypeError(
            "spec_of supports BetheStoppingSource and NistStarStoppingSource only in this step, "
            f"and a source that declares itself by a stopping_spec() method; got "
            f"{type(source).__module__}.{type(source).__qualname__}"
        )
    candidate = _bethe_spec(UNREGISTERED_ID, "1", source, "user")
    for model_id, spec in STOPPING_MODELS.items():
        if _structure(spec) == _structure(candidate):
            return model_id, spec
    return UNREGISTERED_ID, candidate
