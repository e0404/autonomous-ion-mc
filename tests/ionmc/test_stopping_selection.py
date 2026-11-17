"""Single selection point, effective id, identity block, fail-closed domains and roles.

Rows DOM-7, ROLE-1, ROLE-2, ROLE-3, PER-1, PER-2 (plan v3-013 A.3, A.4), PER-3 negative control
(block validator), the single-model assertion (CC-ID/CC-NEG precondition), the R1-ion precondition
(``stopping_model`` and ``stopping=BetheStoppingSource()`` give bit-identical tables) and the block
keys. Decision 0044 sections 2 to 5. No reference data and no table values.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np
import pytest

import ionmc
from ionmc import config as config_module
from ionmc.config import PhysicsOptions, RunOptions, SimulationConfig, validate
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import PMMA, WATER
from ionmc.physics import exact_omega as XO
from ionmc.physics import stopping_models
from ionmc.physics.projectiles import ALPHA, CARBON12, DEUTERON, PROTON
from ionmc.physics.stopping import BetheOptions, BetheStoppingSource
from ionmc.physics.stopping_models import (
    DEFAULT_STOPPING_MODEL_ID,
    STOPPING_MODELS,
    DataPin,
    StoppingModelSpec,
    resolve,
    spec_of,
    spec_sha256,
)
from ionmc.physics.stopping_selection import (
    BLOCK_KEYS,
    assert_single_model,
    check_window,
    effective_model_id,
    select,
    table_model_id,
    validate_identity_block,
)
from ionmc.scoring import ScoringGrid
from ionmc.sources import PencilBeamSource
from ionmc.transport.tables import TransportTables, material_fingerprint, thaw
from tests.ionmc._user_sources import UserSpecMixin, user_spec

REPO = Path(__file__).resolve().parents[2]
LS_EFF = "bethe-ionmc-v1+ls-eff"
LS_BARE = "bethe-ionmc-v1+ls-bare_mixed"


def _cfg(
    proj: Any = PROTON,
    e_u: float | None = None,
    *,
    stopping: Any = None,
    stopping_model: str | None = None,
    sigma: float = 0.0,
    mcs: bool = True,
    geometry: Any = None,
    tallies: tuple[Any, ...] = (),
) -> SimulationConfig:
    e_u = e_u if e_u is not None else (20.0 if proj.z == 1 else 100.0)
    if stopping is None and stopping_model is None:
        stopping_model = DEFAULT_STOPPING_MODEL_ID
    geo = geometry or BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 40.0), WATER)
    return SimulationConfig(
        source=PencilBeamSource(proj, (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), e_u * proj.a, sigma, 0.0),
        geometry=geo,
        scoring=(ScoringGrid((-10.0, -10.0, 0.0), (2.0, 2.0, 2.0), (10, 10, 20)),),
        physics=PhysicsOptions(
            nuclear=False,
            stopping=stopping,
            stopping_model=stopping_model,
            multiple_scattering=mcs,
        ),
        run=RunOptions(backend="python", precision="float64", seed=1, n_histories=2, n_batches=2),
        tallies=tallies,
    )


# ---- DOM-7: the selection point ----


def test_dom7_unknown_id_both_neither_and_non_string() -> None:
    with pytest.raises(UnsupportedCombinationError, match="unknown stopping model id"):
        PhysicsOptions(nuclear=False, stopping_model="no-such-model")
    with pytest.raises(UnsupportedCombinationError, match="exactly one"):
        PhysicsOptions(
            nuclear=False, stopping=BetheStoppingSource(), stopping_model="bethe-ionmc-v1"
        )
    with pytest.raises(UnsupportedCombinationError, match="exactly one"):
        PhysicsOptions(nuclear=False)
    with pytest.raises(UnsupportedCombinationError):
        PhysicsOptions(nuclear=False, stopping_model=3)  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError):
        PhysicsOptions(nuclear=False, stopping=object())  # type: ignore[arg-type]
    for bad in ("nist-star", "bethe", "icru90-water"):  # not registered ids
        with pytest.raises(UnsupportedCombinationError):
            PhysicsOptions(nuclear=False, stopping_model=bad)


def test_dom7_positional_stopping_and_replace_keep_working() -> None:
    src = BetheStoppingSource()
    p = PhysicsOptions(False, src)  # positional (nuclear, stopping), as before
    assert p.stopping is src and p.stopping_source is src and p.stopping_model is None
    assert p.selection.selected_by == "stopping-object"
    q = replace(p, max_step_mm=2.5)
    assert q.stopping_source is src
    m = PhysicsOptions(nuclear=False, stopping_model=LS_EFF)
    assert m.stopping is None and m.selection.selected_by == "stopping_model"
    assert m.stopping_source.options.ls == "eff"
    with pytest.raises(UnsupportedCombinationError, match="exactly one"):
        replace(m, stopping=src)  # both given: refused; an explicit stopping_model=None is needed
    assert replace(m, stopping=src, stopping_model=None).selection.requested_id == "bethe-ionmc-v1"


def test_dom7_unsupported_stopping_object_is_refused_unless_it_declares_itself() -> None:
    class Foreign:
        name = "foreign"

        def table(self, material: Any, projectile: Any) -> Any:  # pragma: no cover
            raise AssertionError

    with pytest.raises(UnsupportedCombinationError, match="Foreign"):
        select(Foreign(), None)  # type: ignore[arg-type]
    with pytest.raises(UnsupportedCombinationError, match="Foreign"):
        PhysicsOptions(nuclear=False, stopping=Foreign())  # type: ignore[arg-type]

    class BadHook(Foreign):
        def stopping_spec(self) -> Any:
            return {"not": "a spec"}

    with pytest.raises(UnsupportedCombinationError, match="StoppingModelSpec"):
        select(BadHook(), None)  # type: ignore[arg-type]

    class Declared(Foreign, UserSpecMixin):
        pass

    sel = select(Declared(), None)  # type: ignore[arg-type]
    assert sel.requested_id == "unregistered" and sel.spec.decision == "user"
    assert sel.spec.kind == "analytic" and sel.spec.data == ()
    assert spec_of(Declared())[0] == "unregistered"  # type: ignore[arg-type]
    # a declared spec cannot claim a registered id or data pins
    claimed = replace(user_spec(Declared()), model_id="bethe-ionmc-v1")
    Declared.stopping_spec = lambda self: claimed  # type: ignore[method-assign]
    assert select(Declared(), None).requested_id == "unregistered"  # type: ignore[arg-type]


# ---- ROLE-1 ----


def _pin(role: str) -> DataPin:
    return DataPin("synthetic", "0" * 64, None, None, "none", "synthetic", role)


@pytest.mark.parametrize("role", ["evaluation-only", "lineage-only"])
def test_role1_non_construction_pin_cannot_be_resolved(
    role: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    default = STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID]
    synth = StoppingModelSpec(
        "synthetic-role-test",
        "1",
        "tabulated",
        {},
        default.interpolation,
        default.range_start,
        (_pin(role),),
        default.domain,
        "test",
    )
    patched = MappingProxyType({**STOPPING_MODELS, synth.model_id: synth})
    monkeypatch.setattr(stopping_models, "STOPPING_MODELS", patched)
    with pytest.raises(ValueError, match=role):
        resolve("synthetic-role-test")
    with pytest.raises(UnsupportedCombinationError, match=role):
        PhysicsOptions(nuclear=False, stopping_model="synthetic-role-test")
    # a construction-permitted pin passes the role gate (then the missing constructor refuses)
    ok = replace(synth, data=(_pin("construction-permitted"),))
    monkeypatch.setattr(
        stopping_models, "STOPPING_MODELS", MappingProxyType({**patched, ok.model_id: ok})
    )
    with pytest.raises(ValueError, match="without a constructor"):
        resolve("synthetic-role-test")
    with pytest.raises(ValueError):
        _pin("construction")  # the role vocabulary is closed


# ---- ROLE-2 ----


def test_role2_no_bethe_default_in_run_paths() -> None:
    with pytest.raises(UnsupportedCombinationError, match="no Bethe default"):
        XO.ExactOmegaMapping.build(CARBON12, WATER)
    with pytest.raises(UnsupportedCombinationError, match="no Bethe default"):
        XO.ExactOmegaMapping.build(CARBON12, WATER, None)
    XO.ExactOmegaMapping.build(PROTON, WATER)  # a proton has no ion table
    with pytest.raises(UnsupportedCombinationError, match="no Bethe default"):
        XO.fdm_domain_status(CARBON12, [WATER], 290.0 * 12, 0.0)
    with pytest.raises(UnsupportedCombinationError, match="no Bethe default"):
        XO.fdm_domain_status(CARBON12, [WATER], 290.0 * 12, 0.0, ion_tables={})
    only_water = {material_fingerprint(WATER): BetheStoppingSource().table(WATER, CARBON12)}
    with pytest.raises(UnsupportedCombinationError, match="no Bethe default"):  # PMMA missing
        XO.fdm_domain_status(CARBON12, [WATER, PMMA], 290.0 * 12, 0.0, ion_tables=only_water)
    XO.fdm_domain_status(CARBON12, [WATER], 290.0 * 12, 0.0, ion_tables=only_water)
    assert XO.fdm_domain_status(PROTON, [WATER], 150.0).fdm_extrapolated is False


# ---- ROLE-3 ----


def _tables(src: Any, proj: Any) -> TransportTables:
    return TransportTables.from_stopping_tables(
        [src.table(m, proj) for m in (WATER, PMMA)], projectile=proj
    )


@pytest.mark.parametrize("variant", [LS_EFF, LS_BARE])
def test_role3_z1_ls_variant_has_the_defaults_effective_id_and_tables(variant: str) -> None:
    for proj in (PROTON, DEUTERON):
        assert effective_model_id(variant, proj) == DEFAULT_STOPPING_MODEL_ID
        a = _tables(resolve(variant), proj)
        b = _tables(resolve(DEFAULT_STOPPING_MODEL_ID), proj)
        assert a.sha256 == b.sha256 and thaw(a.identity) == thaw(b.identity)
        for name in ("ln_s_mass", "ln_r_mass", "r_mass", "f_mass", "d_f", "ln_e_of_r"):
            assert np.array_equal(getattr(a, name), getattr(b, name))
    for proj in (ALPHA, CARBON12):
        assert effective_model_id(variant, proj) == variant
        assert (
            _tables(resolve(variant), proj).sha256
            != _tables(resolve(DEFAULT_STOPPING_MODEL_ID), proj).sha256
        )
    eff = validate(_cfg(PROTON, None, stopping_model=variant))
    blk = eff.summary()["stopping"]
    assert (blk["model_id_requested"], blk["model_id_effective"]) == (
        variant,
        DEFAULT_STOPPING_MODEL_ID,
    )
    assert blk["spec_sha256"] == spec_sha256(STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID])
    assert (eff.ls_requested, eff.ls_effective) == (
        variant.split("-")[-1].replace("ls-", ""),
        "off",
    )
    assert blk["options"]["ls_requested"] != "off" and blk["options"]["ls_effective"] == "off"
    ion = validate(_cfg(ALPHA, None, stopping_model=variant)).summary()["stopping"]
    assert ion["model_id_requested"] == ion["model_id_effective"] == variant
    assert ion["options"]["ls_effective"] == ion["options"]["ls_requested"] != "off"


# ---- identity block, PER-2, default byte-identity of the tables section ----


def test_block_keys_and_default_run_tables_section_unchanged() -> None:
    by_model = validate(_cfg(stopping_model="bethe-ionmc-v1")).summary()
    by_object = validate(_cfg(stopping=BetheStoppingSource())).summary()
    for blk in (by_model["stopping"], by_object["stopping"]):
        assert set(BLOCK_KEYS) <= set(blk)
        assert blk["model_id_requested"] == blk["model_id_effective"] == "bethe-ionmc-v1"
        assert blk["kind"] == "analytic" and blk["version"] == "1" and blk["decision"] == "0044"
        assert blk["data_declared"] == [] and blk["data_used"] == []
        assert blk["content_sha256"] is None and blk["circular_rows"] == []
        assert blk["species_transported"] == ["proton"]
        assert "mcs_proton_path_stopping" not in blk  # a proton run has no MCS mapping
        assert blk["evidence_at_run"][0]["status"] == "as-graded:0043"
        assert blk["options"]["i_value"] == "material"
        assert blk["range_construction"] == by_model["tables"]["materials"][0]["range_construction"]
    assert by_model["stopping"]["selected_by"] == "stopping_model"
    assert by_object["stopping"]["selected_by"] == "stopping-object"
    assert by_model["tables"] == by_object["tables"]
    json.dumps(by_model["stopping"])  # JSON-serialisable as it is
    assert set(by_model) - {"stopping"} == set(
        validate(_cfg(stopping=BetheStoppingSource(e_min_per_u=1.0))).summary()
    ) - {"stopping"}


def test_per2_ion_runs_record_the_mcs_proton_path_and_circular_rows_are_empty() -> None:
    for proj in (ALPHA, CARBON12):
        for kw in ({"stopping_model": "bethe-ionmc-v1"}, {"stopping": BetheStoppingSource()}):
            blk = validate(_cfg(proj, None, **kw)).summary()["stopping"]
            assert blk["mcs_proton_path_stopping"] == "bethe-ionmc-v1"
            assert blk["circular_rows"] == []
            assert blk["species_transported"] == [proj.name]
    off = validate(_cfg(ALPHA, None, stopping_model=LS_EFF, mcs=False)).summary()["stopping"]
    assert off["mcs_proton_path_stopping"] == "bethe-ionmc-v1"  # recorded for every ion run


def test_id_user_unregistered_sources_are_recorded_with_complete_options() -> None:
    for src in (
        BetheStoppingSource(e_max_per_u=200.0),
        BetheStoppingSource(options=BetheOptions(barkas=False)),
        BetheStoppingSource(name="other"),
    ):
        blk = validate(_cfg(stopping=src)).summary()["stopping"]
        assert blk["model_id_requested"] == blk["model_id_effective"] == "unregistered"
        assert blk["decision"] == "user" and blk["evidence_at_run"][0]["status"] == "unvalidated"
        assert blk["options"]["barkas"] is src.options.barkas
        assert blk["spec_sha256"] == spec_sha256(spec_of(src)[1])
    declared = validate(_cfg(stopping=_NamedBethe())).summary()["stopping"]
    assert declared["model_id_effective"] == "unregistered" and declared["decision"] == "user"
    assert declared["spec_sha256"] is not None and declared["kind"] == "analytic"


class _NamedBethe(UserSpecMixin):
    """A duck-typed source (not a supported type) delegating to the default Bethe tables."""

    name = "duck"

    def table(self, material: Any, projectile: Any) -> Any:
        return BetheStoppingSource().table(material, projectile)


# ---- R1-ion precondition: the two selections build bit-identical tables ----


@pytest.mark.parametrize("proj", [PROTON, ALPHA, CARBON12])
def test_r1_ion_precondition_model_and_object_selection_are_bit_identical(proj: Any) -> None:
    a = validate(_cfg(proj, None, stopping_model="bethe-ionmc-v1"))
    b = validate(_cfg(proj, None, stopping=BetheStoppingSource()))
    assert a.tables.sha256 == b.tables.sha256 and thaw(a.tables.identity) == thaw(b.tables.identity)
    for name in (
        "ln_s_mass",
        "ln_r_mass",
        "r_mass",
        "f_mass",
        "d_f",
        "ln_e_of_r",
        "e_min_mev",
        "e_max_mev",
    ):
        assert np.array_equal(getattr(a.tables, name), getattr(b.tables, name))
    assert a.max_steps == b.max_steps
    sa, sb = a.summary(), b.summary()
    assert sa["tables"] == sb["tables"] and sa["mcs_mapping"] == sb["mcs_mapping"]


# ---- single-model assertion ----


def test_single_model_assertion_is_a_real_check_over_the_identities() -> None:
    sel = select(None, "bethe-ionmc-v1")
    ident = lambda src, proj, mat=WATER: thaw(_tables_one(src, proj, mat).identity)  # noqa: E731
    good = ident(resolve("bethe-ionmc-v1"), PROTON)
    assert (
        assert_single_model(sel, [("m", PROTON, good), ("w", PROTON, [None])]) == "bethe-ionmc-v1"
    )
    for other in (
        BetheStoppingSource(options=BetheOptions(barkas=False)),
        BetheStoppingSource(e_max_per_u=200.0),
    ):
        with pytest.raises(UnsupportedCombinationError, match="do not name one stopping model"):
            assert_single_model(sel, [("m", PROTON, good), ("water", PROTON, ident(other, PROTON))])
    # ion: the LS variant names itself for z >= 2 and the default for z = 1
    sel_ls = select(None, LS_EFF)
    alpha_ls = ident(resolve(LS_EFF), ALPHA)
    assert table_model_id(alpha_ls[0], ALPHA) == LS_EFF
    assert assert_single_model(sel_ls, [("m", ALPHA, alpha_ls)]) == LS_EFF
    with pytest.raises(UnsupportedCombinationError):
        assert_single_model(sel_ls, [("m", ALPHA, ident(resolve("bethe-ionmc-v1"), ALPHA))])
    z1 = ident(resolve(LS_EFF), PROTON)
    assert table_model_id(z1[0], PROTON) == "bethe-ionmc-v1"
    assert assert_single_model(sel_ls, [("m", PROTON, z1)]) == "bethe-ionmc-v1"
    # an unregistered source: only tables of that very source pass
    mod = BetheStoppingSource(e_max_per_u=200.0)
    sel_u = select(mod, None)
    assert assert_single_model(sel_u, [("m", PROTON, ident(mod, PROTON))]) == "unregistered"
    with pytest.raises(UnsupportedCombinationError):
        assert_single_model(sel_u, [("m", PROTON, good)])


def _tables_one(src: Any, proj: Any, mat: Any) -> TransportTables:
    return TransportTables.from_stopping_tables(
        [src.table(mat, proj)], water=src.table(WATER, proj), projectile=proj
    )


def test_validate_detects_a_consumer_that_bypasses_the_selected_source(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = _cfg(PROTON, None, stopping_model="bethe-ionmc-v1")
    ion_cfg = _cfg(ALPHA, None, stopping_model=LS_EFF)
    validate(ion_cfg)
    default_tables = config_module._stopping_tables

    def default_for_ions(source: Any, materials: Any, projectile: Any = PROTON) -> Any:
        return default_tables(BetheStoppingSource(), materials, projectile)

    monkeypatch.setattr(config_module, "_stopping_tables", default_for_ions)
    with pytest.raises(UnsupportedCombinationError, match="do not name one stopping model"):
        validate(ion_cfg)
    validate(cfg)  # the patched consumer is the default's: consistent with the default selection


# ---- domain windows (fail closed, naming model, species, material, energy, domain) ----


def test_window_and_species_failures_name_model_species_material_energy_and_domain() -> None:
    sel = select(None, "bethe-ionmc-v1")
    check_window(sel, ALPHA, "water", "x", 1.0 * 4, 500.0 * 4)  # window ends are inside
    with pytest.raises(UnsupportedCombinationError) as ei:
        check_window(sel, ALPHA, "water", "source interval", 100.0, 500.0 * 4 + 1.0)
    msg = str(ei.value)
    for part in ("bethe-ionmc-v1", "alpha", "water", "source interval", "MeV/u", "species"):
        assert part in msg
    with pytest.raises(UnsupportedCombinationError, match="species 'helium3'"):
        from ionmc.physics.projectiles import Projectile

        check_window(sel, Projectile("helium3", "He3", 2, 3, 2808.39), "water", "x", 10.0, 20.0)
    with pytest.raises(UnsupportedCombinationError, match="MeV/u"):  # 501 MeV proton: above 500
        validate(_cfg(PROTON, 501.0))


# ---- PER-3 negative control (block validator) ----


def test_block_validator_rejects_a_pin_of_a_projectile_that_was_not_transported() -> None:
    blk = validate(_cfg(stopping_model="bethe-ionmc-v1")).summary()["stopping"]
    validate_identity_block(blk)
    pin = {"dataset_id": "synthetic", "use_role": "construction-permitted"}
    ok = {
        **blk,
        "data_declared": [pin],
        "data_used": [
            {"pin": pin, "role_in_this_run": "construction", "table": ["proton", "water"]}
        ],
    }
    validate_identity_block(ok)
    bad = {**ok, "data_used": [{**ok["data_used"][0], "table": ["alpha", "water"]}]}
    with pytest.raises(ValueError, match="not transported"):
        validate_identity_block(bad)
    with pytest.raises(ValueError, match="not in data_declared"):
        validate_identity_block({**ok, "data_declared": []})
    with pytest.raises(ValueError, match="construction"):
        validate_identity_block(
            {**ok, "data_used": [{**ok["data_used"][0], "role_in_this_run": "x"}]}
        )
    with pytest.raises(ValueError, match="lacks keys"):
        validate_identity_block({k: v for k, v in blk.items() if k != "decision"})
    with pytest.raises(ValueError, match="differs"):
        validate_identity_block(
            {**blk, "model_id_requested": "bethe-ionmc-v1+ls-eff", "model_id_effective": "x"}
        )


# ---- PER-1: clean-process round trip ----

_READER = r"""
import json, sys
from ionmc import materials
from ionmc.physics.stopping_models import resolve, spec_sha256
from ionmc.physics.stopping_selection import validate_identity_block
from ionmc.physics.projectiles import PROTON, ALPHA, CARBON12
from ionmc.transport.tables import TransportTables, thaw
doc = json.loads(open(sys.argv[1]).read())
blk = doc["stopping"]
validate_identity_block(blk)
proj = {"proton": PROTON, "alpha": ALPHA, "carbon12": CARBON12}[doc["source"]["projectile"]]
src = resolve(blk["model_id_effective"])
tabs = TransportTables.from_stopping_tables([src.table(materials.WATER, proj)], projectile=proj)
from ionmc.physics.stopping_models import STOPPING_MODELS, _canonical
spec = STOPPING_MODELS[blk["model_id_effective"]]
print(json.dumps({"spec_sha256": spec_sha256(spec),
                  "pins": [_canonical(p) for p in spec.data],
                  "identity": json.loads(json.dumps(thaw(tabs.identity), default=str))}))
"""


def _write_like_steps_v8(doc: dict[str, Any], path: Path) -> None:
    """The JSON writer of validation/scripts/transport/steps_v8.py (``json.dumps(doc,
    sort_keys=True, default=base._json)``; ``_json`` of steps.py maps numpy values)."""

    def _json(o: Any) -> Any:
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))

    path.write_text(json.dumps(doc, sort_keys=True, default=_json))


@pytest.mark.parametrize(
    ("proj", "model"),
    [(PROTON, "bethe-ionmc-v1"), (PROTON, LS_EFF), (ALPHA, LS_EFF), (CARBON12, "bethe-ionmc-v1")],
)
def test_per1_round_trip_in_a_clean_process(proj: Any, model: str, tmp_path: Path) -> None:
    summary = validate(_cfg(proj, None, stopping_model=model)).summary()
    doc_path = tmp_path / "step.json"
    _write_like_steps_v8(summary, doc_path)
    out = subprocess.run(  # noqa: S603
        [sys.executable, "-I", "-c", _READER, str(doc_path)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        check=False,
        timeout=300,
    )
    assert out.returncode == 0, out.stderr
    got = json.loads(out.stdout.strip().splitlines()[-1])
    blk = json.loads(doc_path.read_text())["stopping"]
    assert got["spec_sha256"] == blk["spec_sha256"]
    assert got["pins"] == blk["data_declared"]
    assert got["identity"] == json.loads(doc_path.read_text())["tables"]["materials"]


# ---- capabilities lists every registered model ----


def test_capabilities_lists_every_registered_model_with_domain_pins_and_evidence() -> None:
    cap = ionmc.capabilities()["physics"]["stopping_models"]
    assert cap["default"] == DEFAULT_STOPPING_MODEL_ID
    assert tuple(cap["models"]) == tuple(STOPPING_MODELS)
    default = cap["models"][DEFAULT_STOPPING_MODEL_ID]
    assert default["spec_sha256"] == spec_sha256(STOPPING_MODELS[DEFAULT_STOPPING_MODEL_ID])
    assert default["data_pins"] == [] and default["domain"]["species"]
    assert default["evidence"][0]["status"] == "as-graded:0043"
    assert cap["models"][LS_EFF]["evidence"][0]["status"] == "unvalidated"
    json.dumps(cap)


# ---- guards ----


def test_no_environment_override_or_fallback_chain_in_the_selection_modules() -> None:
    for rel in (
        "src/ionmc/physics/stopping_models.py",
        "src/ionmc/physics/stopping_selection.py",
    ):
        text = (REPO / rel).read_text()
        assert "environ" not in text and "getenv" not in text, rel


def test_nist_model_is_registered_and_listed_by_the_capabilities() -> None:
    # registered in S3 (cache-independent spec); its resolve needs the cached tables
    assert "nist-pstar-astar-water-2005" in STOPPING_MODELS
    cap = ionmc.capabilities()["physics"]["stopping_models"]["models"][
        "nist-pstar-astar-water-2005"
    ]
    assert [p["dataset_id"] for p in cap["data_pins"]] == [
        "nist-pstar-water-2005",
        "nist-astar-water-2005",
    ]
    assert cap["evidence"][0]["status"] == "consistency-only"
