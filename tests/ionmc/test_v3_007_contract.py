"""Rows CAP-1 (CI part), CAP-2, CAP-3, CAP-4 (partial), CAP-6, REG-2 and REG-3 of plan v3-007
(step S1): the capability contract extended by the dimensions source kind {pencil, beamlets} x
influence {off, on}. Only ``validate`` / ``reconcile`` run; no transport. Not a clinical claim.

CAP-4 is partial in S1: the dispatch controls that patch the driver ("ignore beamlet ids or the
threshold") need beamlet transport (S2/S3); here the config/contract-level honoured check is
tested (an influence request without its effective implementation raises).
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _v3_007_fixed_configs as fixed  # noqa: E402
import test_v3_006_contract as v6  # noqa: E402
import test_v3_007_beamlets as b7  # noqa: E402

from ionmc import capabilities, config, contract  # noqa: E402
from ionmc.config import InfluenceOptions, validate  # noqa: E402
from ionmc.configio import to_dict  # noqa: E402
from ionmc.errors import UnsupportedCombinationError  # noqa: E402
from ionmc.physics.projectiles import ALPHA, PROTON  # noqa: E402
from ionmc.simulation import Simulation  # noqa: E402

U = UnsupportedCombinationError
bl, cfg = b7.bl, b7.cfg


@pytest.fixture
def gate_open(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(config, "BEAMLET_TRANSPORT_AVAILABLE", True)


def _bad(out: list[contract.Outcome]) -> list[contract.Outcome]:
    return [o for o in out if o.verdict in ("disagree", "third")]


# ---- CAP-1 (CI part): the sweep over source kind x influence ------------------------------------


def test_cap1_sweep_source_kind_x_influence(gate_open: None) -> None:
    m = contract.Matrix(
        species=(PROTON,),
        models=(v6.BETHE,),
        nuclear=(False,),
        tallies=("none", "let"),
        sources=("pencil", "beamlets"),
        influence=(False, True),
    )
    assert m.size() == 1 * 4 * 3 * 2 * 1 * 2 * 2 * 1 * 1 * 1 * 2 * 2 == 384
    out = contract.sweep(m)
    assert len(out) == m.size() and _bad(out) == []
    ok = [o for o in out if o.actual == contract.OK]
    for kind in ("pencil", "beamlets"):
        assert any(o.source == kind and o.influence is False for o in ok), kind
    assert any(o.source == "beamlets" and o.influence for o in ok)
    # influence on a pencil source and influence with a tally are refused, declared as such
    assert all(
        o.actual == contract.UNSUPPORTED for o in out if o.source == "pencil" and o.influence
    )
    assert {o.rule for o in out if o.source == "pencil" and o.influence} == {
        "influence-needs-beamlets"
    } | {
        o.rule
        for o in out
        if o.source == "pencil" and o.influence and o.rule != "influence-needs-beamlets"
    }
    assert all(
        o.actual == contract.UNSUPPORTED
        for o in out
        if o.source == "beamlets" and o.tallies != "none"
    )
    assert all(
        o.actual == contract.UNSUPPORTED
        for o in out
        if o.source == "beamlets" and o.influence and o.trace
    )


def test_cap1_sweep_with_the_s1_gate_counts_every_beamlet_cell_unsupported() -> None:
    assert config.BEAMLET_TRANSPORT_AVAILABLE is False
    m = contract.Matrix(
        species=(PROTON,), models=(v6.BETHE,), domain_classes=("inside",), tallies=("none",),
        trace=(False,), sources=("beamlets",), influence=(False, True),
    )  # fmt: skip
    out = contract.sweep(m)
    assert out and _bad(out) == []
    assert {o.actual for o in out} == {contract.UNSUPPORTED}
    assert {o.rule for o in out} == {"beamlet-transport-not-available"}


def test_cap1_sweep_ions_and_nuclear_cells(
    gate_open: None, monkeypatch: pytest.MonkeyPatch
) -> None:
    tid = v6._patch_synthetic_nuclear(monkeypatch)
    m = contract.Matrix(
        species=(PROTON, ALPHA),
        models=(v6.BETHE,),
        domain_classes=("inside", "above"),
        precisions=("float64",),
        nuclear=(False, True),
        tallies=("none",),
        trace=(False,),
        sources=("pencil", "beamlets"),
        influence=(False, True),
        nuclear_table_id=tid,
    )
    out = contract.sweep(m)
    assert _bad(out) == []
    # nuclear + beamlets: accepted on python (proton), refused on the Warp backends
    nb = [o for o in out if o.source == "beamlets" and o.nuclear and o.species == "proton"
          and o.domain_class == "inside"]  # fmt: skip
    assert {(o.backend, o.actual) for o in nb} == {
        ("python", contract.OK), ("warp-cpu", contract.UNSUPPORTED),
        ("warp-cuda", contract.UNSUPPORTED),
    }  # fmt: skip


# ---- CAP-6 ---------------------------------------------------------------------------------


def test_cap6_warp_nuclear_beamlets_name_v3_005b_python_accepted(
    gate_open: None, monkeypatch: pytest.MonkeyPatch
) -> None:
    tid = v6._patch_synthetic_nuclear(monkeypatch)
    ph = {"nuclear": True, "nuclear_table_id": tid}
    validate(cfg(physics=ph))  # python + nuclear + beamlets accepted (E <= 250 MeV)
    validate(cfg(physics=ph, influence=InfluenceOptions("dose")))
    for backend, prec in (("warp-cpu", "float32"), ("warp-cuda", "float32")):
        with pytest.raises(U, match="V3-005B"):
            validate(cfg(physics=ph, run={"backend": backend, "precision": prec}))
    with pytest.raises(U, match="V3-005B"):
        validate(cfg(physics=ph, run={"backend": "warp-cpu", "precision": "float32"},
                     influence=InfluenceOptions("dose")))  # fmt: skip


# ---- CAP-2: the declaration is derived ----------------------------------------------------------


def test_cap2_beamlet_declaration_follows_the_constants(monkeypatch: pytest.MonkeyPatch) -> None:
    base = capabilities()
    assert base["sources"] == ["PencilBeamSource", "BeamletSet"]
    assert base["declared"]["source_types"] == ["PencilBeamSource", "BeamletSet"]
    d = base["declared"]["beamlets"]
    assert d["max_beamlets"] == config.MAX_BEAMLETS == 4096
    assert d["max_total_histories"] == 2**32 and d["max_cpu_workers"] == 1
    monkeypatch.setattr(config, "MAX_BEAMLETS", 7)
    got = v6._paths(base, capabilities())
    assert got == {"/declared/beamlets/max_beamlets"}
    monkeypatch.undo()
    monkeypatch.setattr(config, "BEAMLET_MAX_CPU_WORKERS", 3)
    assert v6._paths(base, capabilities()) == {"/declared/beamlets/max_cpu_workers"} | {
        p for p in v6._paths(base, capabilities()) if p.startswith("/declared/beamlets/run_comb")
    }
    monkeypatch.undo()
    monkeypatch.setattr(config, "BEAMLET_TRANSPORT_AVAILABLE", True)
    got = v6._paths(base, capabilities())
    assert "/declared/beamlets/transport_available" in got
    assert all(p.startswith("/declared/beamlets/") for p in got), got
    monkeypatch.undo()
    monkeypatch.setattr(config, "INFLUENCE_THRESHOLD_BOUNDS", (0.0, 0.5))
    assert "/declared/beamlets/influence_threshold_rel_range[1]" in v6._paths(base, capabilities())


def test_cap2_no_new_top_level_key_and_free_text_only_under_informative() -> None:
    caps = capabilities()
    pre = {"species", "species_details", "species_roles", "unsupported_species", "backends",
           "trace", "chunk_histories", "max_cpu_workers", "backend_names", "precisions",
           "physics", "geometry", "sources", "scoring", "tallies", "max_scoring_grids",
           "declared", "informative"}  # fmt: skip
    assert set(caps) == pre
    v6.test_cap2_free_text_only_under_informative()


# ---- CAP-3: reconciliation ----------------------------------------------------------------------


def _rq_eff(c: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    return to_dict(c), validate(c).summary()


def test_cap3_beamlet_rows_for_every_new_field(gate_open: None) -> None:
    assert len(contract.RECONCILIATION_TABLE) == 24  # the pencil table of 0045 is unchanged
    assert [r.id for r in contract.BEAMLET_ROWS] == [
        "beamlet-set", "beamlet-direction", "beamlet-histories", "influence", "influence-group",
    ]  # fmt: skip
    bs = (bl(-2.0), bl(0.0), bl(2.0, label="z", weight=1.0))
    on = cfg(beamlets=bs, influence=InfluenceOptions("dose"))
    rq, eff = _rq_eff(on)
    res = {r.row: r for r in contract.reconcile(on, validate(on))}
    assert contract.reconcile(on, validate(on)) == contract.reconcile(rq, eff)
    assert "source" not in res and "direction" not in res
    assert len(res) == len(contract.RECONCILIATION_TABLE) - 1 - 2 + len(contract.BEAMLET_ROWS)
    assert (res["influence-group"].cls, res["influence-group"].rule) == ("derived", "memory-budget")
    assert (res["beamlet-histories"].cls, res["beamlet-histories"].rule) == ("equal", "sum-equal")
    assert (res["beamlet-direction"].cls, res["beamlet-direction"].rule) == (
        "equal", "not-applicable")  # fmt: skip
    for r in res.values():
        assert all(contract.pointer_resolves(eff, p) for p in r.pointers), r.row
    # explicit group size: equal; no influence: unused; TOML-derived n_histories: derived
    ex = cfg(beamlets=bs, influence=InfluenceOptions("dose", group_size=2))
    assert {r.row: (r.cls, r.rule) for r in contract.reconcile(ex, validate(ex))}[
        "influence-group"] == ("equal", "equal")  # fmt: skip
    off = cfg(beamlets=bs)
    got = {r.row: (r.cls, r.rule) for r in contract.reconcile(off, validate(off))}
    assert got["influence"] == ("unused", "no-influence")
    assert got["influence-group"] == ("unused", "no-influence")
    # requested run.n_histories None (derived): the row is derived:sum-of-beamlets
    import dataclasses

    dn = dataclasses.replace(off, run=dataclasses.replace(off.run, n_histories=None))
    dv = contract.reconcile(dn, validate(dn))
    assert {r.row: (r.cls, r.rule) for r in dv}["beamlet-histories"] == (
        "derived",
        "sum-of-beamlets",
    )
    # a non-unit direction is a derived row
    sk = cfg(beamlets=(bl(), bl(1.0)))
    from ionmc.sources import Beamlet

    tilted = cfg(beamlets=(Beamlet((0.0, 0.0, -5.0), (0.0, 0.0, 2.0), 100.0, n_histories=40),
                           bl(1.0)))  # fmt: skip
    got = {r.row: (r.cls, r.rule) for r in contract.reconcile(tilted, validate(tilted))}
    assert got["beamlet-direction"] == ("derived", "normalise-per-beamlet") and sk


@pytest.mark.parametrize(
    ("pointer", "value"),
    [
        ("/source/n_beamlets", 5),
        ("/source/n_histories_total", 81),
        ("/source/envelope/e_hi_mev", 130.0),
        ("/source/beamlets_sha256", "0" * 64),
        ("/unit_direction", [0.0, 0.0, 1.0]),
        ("/unit_direction_status", "ok"),
        ("/influence/grid", "other"),
        ("/influence/threshold_rel", 0.5),
        ("/influence/persist_batches", False),
        ("/influence/n_voxels", 65),
        ("/influence/group_size_origin", "requested"),
        ("/influence/group_size", 99),
    ],
)
def test_cap3_a_patched_effective_beamlet_field_raises(
    gate_open: None, pointer: str, value: Any
) -> None:
    c = cfg(beamlets=(bl(-2.0), bl(0.0), bl(2.0)), influence=InfluenceOptions("dose"))
    rq, eff = _rq_eff(c)
    contract.reconcile(rq, eff)
    node: Any = eff
    parts = pointer.split("/")[1:]
    for p in parts[:-1]:
        node = node[p]
    node[parts[-1]] = value
    with pytest.raises(U, match="requested .*, effective"):
        contract.reconcile(rq, eff)


def test_cap3_changed_request_and_weight_edits_are_caught(gate_open: None) -> None:
    c = cfg(influence=InfluenceOptions("dose"))
    rq, eff = _rq_eff(c)
    bad = copy.deepcopy(rq)
    bad["source"]["beamlets"][0]["label"] = "renamed"  # digest covers labels
    with pytest.raises(U, match="canonical-table sha256"):
        contract.reconcile(bad, eff)
    bad = copy.deepcopy(rq)
    bad["source"]["beamlets"][1]["weight"] = 3.0  # and weights
    with pytest.raises(U, match="canonical-table sha256"):
        contract.reconcile(bad, eff)
    bad = copy.deepcopy(rq)
    bad["run"]["n_histories"] = 81
    with pytest.raises(U, match="n_histories"):
        contract.reconcile(bad, eff)
    w0 = cfg(beamlets=(bl(), bl(1.0, weight=0.5)))
    w1 = cfg(beamlets=(bl(), bl(1.0, weight=2.0)))
    assert to_dict(w0) != to_dict(w1)  # a weight-only change changes the requested block


# ---- CAP-4 (partial): honoured check at the config / contract level ----------------


def test_cap4_an_influence_request_without_its_effective_object_raises(
    gate_open: None, monkeypatch: pytest.MonkeyPatch
) -> None:
    c = cfg(influence=InfluenceOptions("dose"))
    rq, eff = _rq_eff(c)
    dropped = copy.deepcopy(eff)
    dropped["influence"] = None
    with pytest.raises(U, match="influence.grid"):
        contract.reconcile(rq, dropped)
    del dropped["influence"]
    with pytest.raises(U, match="influence.grid.*<absent>"):
        contract.reconcile(rq, dropped)
    # an effective influence block nobody requested
    off_rq, off_eff = _rq_eff(cfg())
    off_eff["influence"] = eff["influence"]
    with pytest.raises(U, match="influence"):
        contract.reconcile(off_rq, off_eff)
    # Simulation construction runs validate + reconcile: a consumer that drops the group size
    real = config.validate

    def dropping(cc: Any) -> Any:
        e = real(cc)
        object.__setattr__(e, "influence_group_origin", None)
        return e

    monkeypatch.setattr("ionmc.simulation.validate", dropping)
    with pytest.raises(U, match="influence"):
        Simulation(c)
    monkeypatch.undo()
    monkeypatch.setattr(config, "BEAMLET_TRANSPORT_AVAILABLE", True)
    Simulation(c)  # the honest effective configuration is accepted


# ---- REG-2 / REG-3 ----------------------------------------------------------------------------

GOLDEN = {  # (summary sha256, to_dict sha256) computed from the 3f3736a6 tree with the same module
    "alpha-mcs-off": ("eb19b3162cd0c5ea56182b66eb1172b6caf2f50c1e16561a8082bc75f3426223",
                      "8af0560858a1d355dbbd36a3b1ae58fb52be9e8980d26a831a5f26200864e06c"),
    "alpha-mcs-on": ("13f6a30f37027caf9288a3f95077e6eead243487b788d1ff59e51862dd827728",
                     "e48af1260bde77a0fc7f21aaed486827b695d66df1f478d96834059e048634ba"),
    "proton-let": ("0175980a16783e5a9326372694d471275a85bdba5377716f1b8343e13b65ba56",
                   "5430e6b30b5dd0d127f7c28d77b2f1317deb2cc074e078efae88f70ba41e89b4"),
    "proton-no-straggling": ("189a9c4b72b929805f0a751df91177b99f6bcbd486a8d8f42f24c3bef6770b37",
                             "efd824d21afb701de7dc5c5ab4bf78dab079553ce9b6acc719530259524b0e53"),
    "proton-warp": ("928409998d66efa7f886765dbc7f20bd87ee63b844e28b7b5b5326fdc1f16f4b",
                    "db9175d095630fceb7f2d448d726ef82f5cb720b1298127e1bf6d6af3daa36e4"),
}  # fmt: skip


def test_reg2_pencil_summary_and_request_are_byte_identical_to_3f3736a6() -> None:
    got = fixed.digests()
    assert got == GOLDEN
    for c in fixed.fixed_configs().values():
        s = validate(c).summary()
        assert "influence" not in s and "unit_direction_status" not in s
        assert "influence" not in to_dict(c) and "type" not in to_dict(c)["source"]
        assert isinstance(s["unit_direction"], list)


# Result fields and capabilities() top-level keys at commit 3f3736a6 (V3-006 merge), literal.
REQUIRED_RESULT_FIELDS = frozenset({
    "valid", "requested_config", "effective_config", "backend", "device", "precision", "seed",
    "rng", "n_histories", "n_batches", "grids", "energy_balance", "counters", "timings",
    "environment", "diagnostics", "transport_report", "channel_raw", "lookups",
})  # fmt: skip
REQUIRED_CAPABILITY_KEYS = frozenset({
    "species", "species_details", "species_roles", "unsupported_species", "backends", "trace",
    "chunk_histories", "max_cpu_workers", "backend_names", "precisions", "physics", "geometry",
    "sources", "scoring", "tallies", "max_scoring_grids", "declared", "informative",
})  # fmt: skip


def test_reg3_capabilities_and_result_fields_of_3f3736a6_are_still_present() -> None:
    caps = capabilities()
    for key in ("species", "backends", "physics", "geometry", "sources", "scoring", "tallies",
                "max_scoring_grids", "declared", "informative", "chunk_histories",
                "max_cpu_workers", "backend_names", "precisions", "trace"):  # fmt: skip
        assert key in caps, key
    d = caps["declared"]
    for key in ("species", "source_window_by_model", "source_energy_rules", "run_combinations",
                "limits", "geometry_types", "source_types", "informative"):  # fmt: skip
        assert key in d, key
    assert "PencilBeamSource" in caps["sources"] and "PencilBeamSource" in d["source_types"]
    from dataclasses import fields

    from ionmc.simulation import Result

    names = {f.name for f in fields(Result)}
    assert REQUIRED_RESULT_FIELDS <= names
    assert REQUIRED_CAPABILITY_KEYS <= set(caps)
    for dropped in ("effective_config", "grids"):  # negative control: a missing name is caught
        assert not REQUIRED_RESULT_FIELDS <= names - {dropped}
    assert not REQUIRED_CAPABILITY_KEYS <= set(caps) - {"sources"}
    assert "effective_config" in names
    assert [r.id for r in contract.RECONCILIATION_TABLE][:2] == ["source", "direction"]


# --- tier-1 binding of the derived beamlet count (plan A-1 item 1; Codex REVIEW-2b383901) ---


def _hist_doc(requested: int | None, run_n: int, eff_n: int, rows: list[int] | None,
              src_type: str = "beamlets") -> tuple[dict, dict, dict]:  # fmt: skip
    rq: dict = {"run": {"n_histories": requested}}
    if rows is not None:
        rq["source"] = {"type": src_type, "beamlets": [{"n_histories": n} for n in rows]}
    return rq, {"n_histories": run_n}, {"n_histories": eff_n}


def test_tier1_derived_histories_bind_to_the_persisted_beamlet_table() -> None:
    """A requested ``null`` binds run/effective N to the sum of the beamlet table; every
    deviation is a ``ResultFormatError`` naming the field."""
    from ionmc import results
    from ionmc.results import ResultFormatError

    results._bind_histories(*_hist_doc(None, 80, 80, [40, 40]))  # derived, consistent
    results._bind_histories(*_hist_doc(80, 80, 80, [40, 40]))  # explicit equal
    results._bind_histories(*_hist_doc(80, 80, 80, None))  # a pencil source, explicit
    bad = [
        (_hist_doc(None, 80, 80, None), "not a beamlet set"),  # null without a beamlet table
        (_hist_doc(None, 80, 80, [40, 40], "pencil"), "not a beamlet set"),
        (_hist_doc(None, 84, 84, [40, 40]), "run.n_histories"),  # run N != sum
        (_hist_doc(None, 80, 84, [40, 40]), "effective.n_histories"),  # effective N != sum
        (_hist_doc(None, 80, 80, []), "empty"),
        (_hist_doc(None, 80, 80, [40, 0]), "positive integer"),
        (_hist_doc(84, 80, 80, [40, 40]), "run.n_histories"),  # explicit but run disagrees
        (_hist_doc(80, 80, 84, [40, 40]), "effective.n_histories"),
    ]
    for (rq, run, eff), needle in bad:
        with pytest.raises(ResultFormatError, match=needle):
            results._bind_histories(rq, run, eff)


# ---- fix 4: persistence gate (S5) and source.type structural validation -------------------------


def test_beamlet_persistence_gate_is_closed_and_skip_is_conditioned() -> None:
    import inspect

    assert config.BEAMLET_PERSISTENCE_AVAILABLE is False
    src = inspect.getsource(b7.test_beam1_derived_n_histories_in_the_bundle_s2_gate)
    assert "BEAMLET_PERSISTENCE_AVAILABLE" in src  # S5 un-skips by the constant


def test_beamlet_result_bundle_entry_points_refuse_naming_s5(tmp_path: Path) -> None:
    from types import SimpleNamespace

    from ionmc import report, results

    fake = SimpleNamespace(requested_config=cfg(), effective_config=None)
    with pytest.raises(U, match=r"V3-007 S5"):
        results.build_document(fake, None, [], True, "x")  # type: ignore[arg-type]
    with pytest.raises(U, match=r"result bundles of a beamlet set are not yet available"):
        results.save(fake, tmp_path / "b")  # type: ignore[arg-type]
    assert not (tmp_path / "b").exists()
    doc = {"requested": {"source": {"type": "beamlets"}}}
    with pytest.raises(U, match=r"V3-007 S5"):
        report.render_all(doc, {})  # type: ignore[arg-type]


def test_validate_dict_structure_checks_source_type() -> None:
    """Tier 1 (``results._check_structure``) wraps this exact call for ``doc["requested"]``; a
    full minimal result document needs a transport run, so the check is exercised directly."""
    from ionmc import configio
    from ionmc.configio import ConfigError

    d = copy.deepcopy(to_dict(cfg()))
    assert d["source"]["type"] == "beamlets"
    configio.validate_dict_structure(d)
    p = copy.deepcopy(to_dict(next(iter(fixed.fixed_configs().values()))))
    assert p["source"].get("type", "pencil") == "pencil"
    configio.validate_dict_structure(p)
    p["source"]["type"] = "pencil"
    configio.validate_dict_structure(p)
    for dd in (d, p):
        for bad in ("spots", 3, None, ["beamlets"]):
            dd["source"]["type"] = bad
            with pytest.raises(ConfigError, match="source.type"):
                configio.validate_dict_structure(dd)
