"""Stopping-model registry with the default (decision 0044 sections 1, 2, 4, 5; V3-013 step S1).

Rows ID-DEF, ID-GUARD, ID-DIAG, ID-USER (plan v3-013 A.1) and CC-HASH part 1 (A.2). The fixture
``data/stopping_identity_a9563412.json`` holds the values computed at the base commit a9563412
before any V3-013 source change (``TransportTables.identity``, its sha256 and the table metadata
of the default and of the two Lindhard-Sorensen variants).
"""

from __future__ import annotations

import ast
import json
from dataclasses import replace
from pathlib import Path
from types import MappingProxyType

import pytest

from ionmc import materials
from ionmc.physics import stopping_evidence, stopping_models
from ionmc.physics.projectiles import ALPHA, CARBON12, OXYGEN16, PROTON
from ionmc.physics.stopping import BetheOptions, BetheStoppingSource, default_source
from ionmc.physics.stopping_models import (
    DEFAULT_STOPPING_MODEL_ID,
    STOPPING_DEFAULT_INTENDED_CHANGE,
    STOPPING_MODELS,
    resolve,
    spec_of,
    spec_sha256,
)
from ionmc.transport.tables import TransportTables, thaw

REPO = Path(__file__).resolve().parents[2]
FIXTURE = Path(__file__).parent / "data" / "stopping_identity_a9563412.json"
PROJECTILES = {"proton": PROTON, "alpha": ALPHA, "carbon12": CARBON12, "oxygen16": OXYGEN16}
MATERIALS = {"WATER": materials.WATER, "PMMA": materials.PMMA}

# Pinned literals (ID-GUARD): changing either needs a STOPPING_DEFAULT_INTENDED_CHANGE record
# that names a decision (decision 0044 section 5).
PINNED_DEFAULT_ID = "bethe-ionmc-v1"
PINNED_DEFAULT_SPEC_SHA256 = "3aa2d5b987c77afb228f4c015e579be7c9935fbbf1cfdb6accda44ccd9b9cc26"
REGISTERED_S1 = ("bethe-ionmc-v1", "bethe-ionmc-v1+ls-eff", "bethe-ionmc-v1+ls-bare_mixed")


def _fixture() -> dict:
    return json.loads(FIXTURE.read_text())


def _js(value: object) -> object:
    return json.loads(json.dumps(thaw(value), default=str))


def _tables(source: BetheStoppingSource, name: str) -> TransportTables:
    return TransportTables.from_stopping_tables(
        [source.table(m, PROJECTILES[name]) for m in MATERIALS.values()],
        projectile=PROJECTILES[name],
    )


# ---- registry ----


def test_registry_holds_the_s1_models_only_and_is_immutable() -> None:
    assert tuple(STOPPING_MODELS) == (*REGISTERED_S1, "nist-pstar-astar-water-2005")  # S3
    assert isinstance(STOPPING_MODELS, MappingProxyType)
    with pytest.raises(TypeError):
        STOPPING_MODELS["x"] = STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID]  # type: ignore[index]
    spec = STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID]
    with pytest.raises(AttributeError):
        spec.version = "2"  # type: ignore[misc]
    with pytest.raises(TypeError):
        spec.options["ls"] = "eff"  # type: ignore[index]
    assert not hasattr(spec, "evidence") and not hasattr(spec, "status")


def test_resolve_unknown_id_raises_and_never_returns_another_id() -> None:
    for bad in ("nist-star", "icru90-water-v1", "", "bethe"):
        with pytest.raises(ValueError, match="unknown stopping model id"):
            resolve(bad)
    for mid in REGISTERED_S1:
        assert spec_of(resolve(mid))[0] == mid


def test_default_resolves_to_the_default_source_and_default_source_uses_it() -> None:
    assert resolve(DEFAULT_STOPPING_MODEL_ID) == BetheStoppingSource()
    assert default_source() == BetheStoppingSource()
    assert resolve("bethe-ionmc-v1+ls-eff").options == BetheOptions(ls="eff")
    assert resolve("bethe-ionmc-v1+ls-bare_mixed").options == BetheOptions(ls="bare_mixed")


def test_spec_content_of_the_registered_models() -> None:
    for mid, spec in STOPPING_MODELS.items():
        if mid not in REGISTERED_S1:  # the NIST model: tests/ionmc/test_stopping_nist_model.py
            continue
        assert spec.model_id == mid and spec.decision == "0044" and spec.version == "1"
        assert spec.kind == "analytic" and spec.data == ()
        assert spec.domain.species == {"proton", "deuteron", "alpha", "carbon12", "oxygen16"}
        assert set(spec.domain.e_min_per_u.values()) == {1.0}
        assert set(spec.domain.e_max_per_u.values()) == {500.0}
        assert spec.options["points_per_decade"] == 200
    assert STOPPING_MODELS["bethe-ionmc-v1"].options["ls"] == "off"
    assert STOPPING_MODELS["bethe-ionmc-v1+ls-eff"].options["ls"] == "eff"
    assert STOPPING_MODELS["bethe-ionmc-v1+ls-bare_mixed"].options["ls"] == "bare_mixed"
    hashes = {spec_sha256(s) for s in STOPPING_MODELS.values()}
    assert len(hashes) == len(STOPPING_MODELS)


def test_spec_sha256_is_stable_canonical_json() -> None:
    spec = STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID]
    assert spec_sha256(spec) == spec_sha256(replace(spec))
    assert spec_sha256(replace(spec, version="2")) != spec_sha256(spec)


# ---- ID-GUARD ----


def test_default_id_and_spec_hash_are_pinned_and_guarded() -> None:
    assert STOPPING_DEFAULT_INTENDED_CHANGE is None
    assert DEFAULT_STOPPING_MODEL_ID == PINNED_DEFAULT_ID
    assert spec_sha256(STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID]) == PINNED_DEFAULT_SPEC_SHA256


def _guard(default_id: str, spec_hash: str, record: object) -> None:
    """The check a future change of the default would have to pass (mirrors the A16 guard)."""
    if (default_id, spec_hash) != (PINNED_DEFAULT_ID, PINNED_DEFAULT_SPEC_SHA256):
        assert isinstance(record, dict) and str(record.get("decision", "")).strip(), (
            "a change of the default model or of its spec_sha256 needs a "
            "STOPPING_DEFAULT_INTENDED_CHANGE record naming a decision"
        )


def test_guard_fails_a_synthetic_change_without_a_decision_record() -> None:
    _guard(PINNED_DEFAULT_ID, PINNED_DEFAULT_SPEC_SHA256, None)
    for ident in (("icru90-water-v1", PINNED_DEFAULT_SPEC_SHA256), (PINNED_DEFAULT_ID, "0" * 64)):
        with pytest.raises(AssertionError, match="names a decision|naming a decision"):
            _guard(*ident, None)
        with pytest.raises(AssertionError):
            _guard(*ident, {"decision": " "})
        _guard(*ident, {"decision": "0099"})
    spec = replace(STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID], version="2")
    assert spec_sha256(spec) != PINNED_DEFAULT_SPEC_SHA256


# ---- evidence registry ----


def test_evidence_registry_is_separate_keyed_by_spec_hash_and_initial() -> None:
    ev = stopping_evidence.STOPPING_EVIDENCE
    assert isinstance(ev, MappingProxyType)
    assert set(ev) == {spec_sha256(s) for s in STOPPING_MODELS.values()}
    status = {mid: ev[spec_sha256(s)][0].status for mid, s in STOPPING_MODELS.items()}
    assert status == {
        "bethe-ionmc-v1": "as-graded:0043",
        "bethe-ionmc-v1+ls-eff": "unvalidated",
        "bethe-ionmc-v1+ls-bare_mixed": "unvalidated",
        "nist-pstar-astar-water-2005": "consistency-only",
    }
    assert all(r.spec_sha256 == h for h, rs in ev.items() for r in rs)
    assert not any("validated" == r.status for rs in ev.values() for r in rs)
    assert stopping_models.STOPPING_MODELS is STOPPING_MODELS  # no evidence import cycle


# ---- spec_of / ID-USER ----


def test_spec_of_round_trip_for_every_registered_model() -> None:
    for mid, spec in STOPPING_MODELS.items():
        if (
            mid not in REGISTERED_S1
        ):  # NIST: cache-dependent, tests/ionmc/test_stopping_nist_model.py
            continue
        got_id, got_spec = spec_of(resolve(mid))
        assert got_id == mid and got_spec == spec and got_spec.decision == "0044"


def test_id_user_modified_sources_are_unregistered_with_complete_options() -> None:
    assert spec_of(BetheStoppingSource())[0] == "bethe-ionmc-v1"
    assert spec_of(BetheStoppingSource(e_min_per_u=1))[0] == "bethe-ionmc-v1"  # 1 == 1.0
    cases = {
        "e_min": (BetheStoppingSource(e_min_per_u=1.5), None),
        "e_max": (BetheStoppingSource(e_max_per_u=20), None),
        "barkas": (BetheStoppingSource(options=BetheOptions(barkas=False)), ("barkas", False)),
        "points": (BetheStoppingSource(points_per_decade=100), ("points_per_decade", 100)),
        "name": (BetheStoppingSource(name="other"), ("source_name", "other")),
    }
    for label, (src, option) in cases.items():
        mid, spec = spec_of(src)
        assert mid == "unregistered", label
        assert spec.model_id == "unregistered" and spec.decision == "user", label
        assert spec.kind == "analytic"
        if option is not None:
            assert spec.options[option[0]] == option[1], label
        assert set(spec.options) == set(STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID].options)
        assert spec_sha256(spec) not in {spec_sha256(s) for s in STOPPING_MODELS.values()}
    assert spec_of(cases["e_min"][0])[1].domain.e_min_per_u["proton"] == 1.5
    assert spec_of(cases["e_max"][0])[1].domain.e_max_per_u["alpha"] == 20.0


def test_spec_of_fails_closed_for_an_unsupported_source_type() -> None:
    class Foreign:
        name = "foreign"

        def table(self, material, projectile):  # pragma: no cover - never called
            raise AssertionError

    with pytest.raises(TypeError, match="spec_of supports"):
        spec_of(Foreign())  # type: ignore[arg-type]


# ---- ID-DIAG ----

_DIAGNOSTIC_FILES = (
    "src/ionmc/transport/parity.py",
    "src/ionmc/transport/ion_mcs_checks.py",
    "src/ionmc/nuclear/build.py",
    "validation/scripts/transport/a16_digest.py",
)


def test_id_diag_diagnostic_sources_equal_the_default_model() -> None:
    """Every ``BetheStoppingSource`` constructed in the diagnostics and fixtures is the
    argument-free default, hence equal to ``resolve("bethe-ionmc-v1")`` (dataclass equality)."""
    for rel in _DIAGNOSTIC_FILES:
        tree = ast.parse((REPO / rel).read_text())
        calls = [
            n
            for n in ast.walk(tree)
            if isinstance(n, ast.Call)
            and isinstance(n.func, ast.Name)
            and n.func.id == "BetheStoppingSource"
        ]
        assert calls, rel
        for call in calls:
            assert not call.args and not call.keywords, (rel, call.lineno)
    assert BetheStoppingSource() == resolve("bethe-ionmc-v1")


# ---- ID-DEF and CC-HASH part 1 (values computed at a9563412) ----


def test_id_def_default_table_identity_equals_the_base_commit() -> None:
    fx = _fixture()
    assert fx["base_commit"] == "a9563412"
    src = resolve(DEFAULT_STOPPING_MODEL_ID)
    for name in PROJECTILES:
        tabs = _tables(src, name)
        assert _js(tabs.identity) == fx["identity"][name], name
        # the default's table metadata gains no model_id or any other key (decision 0044 section 4)
        for ident in _js(tabs.identity):
            assert "model_id" not in ident["metadata"]  # type: ignore[index]


@pytest.mark.host
def test_id_def_default_table_sha256_equals_the_base_commit_on_the_same_host() -> None:
    # machine dependent (platform libm ulps, decision 0044 Q4): the fixture was computed on the
    # development host; compared separately from the portable identity assertion above and run
    # only under the controlled host runner (plan row ID-DEF: sha256 in LV; Codex REVIEW-7e3c915c).
    fx = _fixture()
    src = resolve(DEFAULT_STOPPING_MODEL_ID)
    for name in PROJECTILES:
        assert _tables(src, name).sha256 == fx["sha256"][name], name


def test_cc_hash_registered_model_metadata_is_byte_identical_to_the_base_commit() -> None:
    fx = _fixture()
    labels = {
        "default": "bethe-ionmc-v1",
        "ls_eff": "bethe-ionmc-v1+ls-eff",
        "ls_bare_mixed": "bethe-ionmc-v1+ls-bare_mixed",
    }
    for key, per_model in fx["metadata"].items():
        pname, mname = key.split("/")
        for label, mid in labels.items():
            meta = resolve(mid).table(MATERIALS[mname], PROJECTILES[pname]).metadata
            assert json.dumps(_js(meta), sort_keys=True) == json.dumps(
                per_model[label], sort_keys=True
            ), (key, mid)
            assert "model_id" not in meta and "version" not in meta


def test_cc_hash_z1_ls_variants_share_the_default_metadata_by_design() -> None:
    fx = _fixture()
    for mname in MATERIALS:
        m = fx["metadata"][f"proton/{mname}"]
        assert m["ls_eff"] == m["default"] and m["ls_bare_mixed"] == m["default"]
        for pname in ("alpha", "carbon12", "oxygen16"):  # z >= 2: the LS mode is in the metadata
            z = fx["metadata"][f"{pname}/{mname}"]
            assert z["ls_eff"]["options"]["ls"] == "eff"
            assert z["ls_bare_mixed"]["options"]["ls"] == "bare_mixed"
            assert "ls" not in z["default"]["options"]
            assert z["ls_eff"] != z["default"] != z["ls_bare_mixed"]
