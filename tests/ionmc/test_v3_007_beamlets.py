"""Rows BEAM-1, BEAM-2, BEAM-3 of plan v3-007 (step S1): the beamlet source model, its fail-closed
validation rules (each provoked separately with its positive twin), the source-consumer audit and
the offsets / batch membership. Only ``validate`` runs here; no transport. Not a clinical claim.

``config.BEAMLET_TRANSPORT_BACKENDS`` is the temporary transport gate (S2: python only, S3 adds the
Warp backends); the tests of the rules open it for every backend, one test pins the gate itself.
"""

from __future__ import annotations

import ast
import copy
import dataclasses
import re
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from ionmc import config, configio
from ionmc.config import (
    DiagnosticsOptions,
    InfluenceOptions,
    PhysicsOptions,
    RunOptions,
    SimulationConfig,
    allocation_estimate,
    validate,
)
from ionmc.configio import ConfigError
from ionmc.errors import UnsupportedCombinationError
from ionmc.geometry import BoxPhantom
from ionmc.materials import WATER
from ionmc.physics.projectiles import ALPHA, PROTON
from ionmc.scoring import ScoringGrid, TallyRequest
from ionmc.sources import (
    UNIT_DIRECTION_STATUS_BEAMLETS,
    Beamlet,
    BeamletSet,
    PencilBeamSource,
)

U = UnsupportedCombinationError
BETHE = "bethe-ionmc-v1"
GRID = ScoringGrid((-10.0, -10.0, 0.0), (5.0, 5.0, 5.0), (4, 4, 4))  # name "dose", 64 voxels


@pytest.fixture(autouse=True)
def _transport_gate_open(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> None:
    if "gate" not in request.node.name:
        monkeypatch.setattr(
            config, "BEAMLET_TRANSPORT_BACKENDS", ("python", "warp-cpu", "warp-cuda")
        )


def bl(x: float = 0.0, e: float = 100.0, n: int = 40, **kw: Any) -> Beamlet:
    return Beamlet((x, 0.0, -5.0), (0.0, 0.0, 1.0), e, n_histories=n, **kw)


def cfg(
    beamlets: tuple[Beamlet, ...] | None = None,
    projectile: Any = PROTON,
    influence: InfluenceOptions | None = None,
    run: dict[str, Any] | None = None,
    physics: dict[str, Any] | None = None,
    diagnostics: DiagnosticsOptions | None = None,
    tallies: tuple[TallyRequest, ...] = (),
    source: Any = None,
) -> SimulationConfig:
    bs = beamlets if beamlets is not None else (bl(-2.0), bl(2.0))
    src = source if source is not None else BeamletSet(projectile, bs)
    n = sum(b.n_histories for b in bs)
    rn = {"backend": "python", "precision": "float64", "seed": 5, "n_histories": n, "n_batches": 4}
    ph = {"nuclear": False, "stopping_model": BETHE, "multiple_scattering": projectile.z < 2}
    return SimulationConfig(
        src,
        BoxPhantom((-10.0, -10.0, 0.0), (20.0, 20.0, 20.0), WATER),
        (GRID,),
        PhysicsOptions(**(ph | (physics or {}))),
        RunOptions(**(rn | (run or {}))),
        diagnostics or DiagnosticsOptions(),
        tallies=tallies,
        influence=influence,
    )


def refused(c: SimulationConfig, needle: str) -> None:
    with pytest.raises(U, match=re.escape(needle)):
        validate(c)


# ---- BEAM-1: every rule with its positive twin --------------------------------------------------


def test_beam1_positive_twin_and_effective_summary() -> None:
    eff = validate(cfg())
    s = eff.summary()
    assert eff.unit_direction is None and s["unit_direction"] is None
    assert (
        s["unit_direction_status"] == UNIT_DIRECTION_STATUS_BEAMLETS == "not-applicable: beamlets"
    )
    assert s["source"]["n_beamlets"] == 2 and s["n_histories"] == 80
    assert s["influence"] is None
    with pytest.raises(ValueError, match="no single direction"):
        eff.single_unit_direction()


def test_beam1_gate_beamlet_transport_python_only_until_s3() -> None:
    """The temporary transport gate (S2): the python backend transports beamlet sets, the Warp
    backends are refused naming S3 (TODO: removed by S3); pencil sources are unaffected."""
    assert config.BEAMLET_TRANSPORT_BACKENDS == ("python",)
    validate(cfg())
    for backend in ("warp-cpu", "warp-cuda"):
        refused(cfg(run={"backend": backend}), r"not yet available on the " + backend)
        refused(cfg(run={"backend": backend}), r"V3-007 S3")
    pencil = PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 100.0)
    validate(cfg(beamlets=(bl(),), source=pencil, run={"n_histories": 40}))  # pencil unaffected


def test_beam1_beamlet_count(monkeypatch: pytest.MonkeyPatch) -> None:
    refused(cfg(beamlets=(), run={"n_histories": 40}), "between 1 and 4096 beamlets, got 0")
    validate(cfg(beamlets=(bl(),)))
    monkeypatch.setattr(config, "MAX_BEAMLETS", 2)
    refused(cfg(beamlets=(bl(-1.0), bl(0.0), bl(1.0))), "between 1 and 2 beamlets, got 3")
    validate(cfg(beamlets=(bl(-1.0), bl(0.0))))


def test_beam1_n_histories_per_beamlet() -> None:
    refused(cfg(beamlets=(bl(n=2), bl(n=40))), "n_histories (2) of beamlet 0")
    refused(cfg(beamlets=(bl(n=40), bl(n=42))), "n_histories (42) of beamlet 1")
    validate(cfg(beamlets=(bl(n=4), bl(n=40))))


def test_beam1_total_histories_limit(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(config, "MAX_TOTAL_HISTORIES", 60)
    refused(cfg(), "sum to 80, above the limit 60")
    validate(cfg(beamlets=(bl(n=28), bl(n=28))))


def test_beam1_run_n_histories_must_equal_the_sum() -> None:
    with pytest.raises(U) as ei:
        validate(cfg(run={"n_histories": 100}))
    assert "run.n_histories (100)" in str(ei.value) and "(80)" in str(ei.value)
    validate(cfg(run={"n_histories": 80}))


@pytest.mark.parametrize(
    "label",
    ["x" * 65, "has space", "a<b", "a[0]", "f(x)", "a|b", "tab\tx", "ctl\x01", "café", ""],
)
def test_beam1_label_grammar_refused(label: str) -> None:
    refused(cfg(beamlets=(bl(label=label), bl())), "does not match ^[A-Za-z0-9_.-]{1,64}$")


def test_beam1_label_positive_and_unique() -> None:
    validate(cfg(beamlets=(bl(label="a.b-c_1"), bl(label="x" * 64), bl())))
    refused(cfg(beamlets=(bl(label="same"), bl(label="same"))), "is not unique")


@pytest.mark.parametrize("w", [float("nan"), float("inf"), float("-inf"), -1.0, "1", True])
def test_beam1_weight_domain_refused(w: Any) -> None:
    refused(cfg(beamlets=(bl(), bl(weight=w))), "weight of beamlet 1")


def test_beam1_weight_none_and_zero_accepted() -> None:
    validate(cfg(beamlets=(bl(weight=None), bl(weight=0.0), bl(weight=2))))


def test_beam1_energy_domain_is_per_beamlet_not_the_mean() -> None:
    # mean energy 900 MeV is inside the alpha domain, beamlet 1 (E + 6 sigma) is above it
    ok = (bl(e=800.0, label="lo"), bl(e=900.0, label="mid"))
    validate(cfg(beamlets=ok, projectile=ALPHA))
    bad = (bl(e=800.0), bl(e=996.0, label="hi", energy_sigma_mev=8.0))
    refused(cfg(beamlets=bad, projectile=ALPHA), "alpha source beamlet 1 ('hi') energy")
    refused(cfg(beamlets=(bl(e=3.0, label="low"), bl())), "source beamlet 0 ('low') energy 3.0")


def test_beam1_sigma_rule_is_per_beamlet() -> None:
    validate(cfg(beamlets=(bl(energy_sigma_mev=0.5), bl(energy_sigma_mev=4.0))))
    refused(cfg(beamlets=(bl(), bl(energy_sigma_mev=6.0))), "energy_sigma_mev beamlet 1 (6.0)")


def test_beam1_envelope_rules_max_steps_and_capacity(monkeypatch: pytest.MonkeyPatch) -> None:
    two = (bl(e=100.0), bl(2.0, e=140.0, energy_sigma_mev=1.0))
    got = validate(cfg(beamlets=two)).max_steps
    pencil = PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 146.0)
    want = validate(cfg(beamlets=(bl(),), source=pencil, run={"n_histories": 40})).max_steps
    assert got == want  # derived from the envelope E_hi = 146
    assert two and BeamletSet(PROTON, two).envelope().e_cap_mev == 146.0
    # the fixed-point capacity check sees e_cap = max_j (E_j + 6 sigma_j)
    monkeypatch.setattr(config, "MAX_QUANTA", int(20 * 100.0 / config.QUANTUM_MEV) + 2)
    validate(cfg(beamlets=(bl(e=100.0), bl(2.0, e=100.0))))
    refused(cfg(beamlets=(bl(e=100.0), bl(2.0, e=120.0))), "exceed the capacity")


def test_beam1_cpu_workers() -> None:
    refused(cfg(run={"cpu_workers": 2}), "must be <= 1 for a beamlet set")
    validate(cfg(run={"cpu_workers": 1}))


def test_beam1_tallies_and_lookups_refused_with_a_beamlet_set() -> None:
    let = (TallyRequest("let", "dose", "let_d"),)
    refused(cfg(tallies=let), "not supported with a beamlet set")
    refused(cfg(tallies=let, influence=InfluenceOptions("dose")), "not supported with")
    validate(cfg())


def test_beam1_influence_rules() -> None:
    validate(cfg(influence=InfluenceOptions("dose")))
    pencil = PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 100.0)
    refused(
        cfg(beamlets=(bl(),), source=pencil, run={"n_histories": 40},
            influence=InfluenceOptions("dose")),
        "influence requires a BeamletSet",
    )  # fmt: skip
    refused(cfg(influence=InfluenceOptions("nope")), "influence.grid 'nope' is not a scoring grid")
    for tau in (1.0, 1.5, -0.1):
        refused(cfg(influence=InfluenceOptions("dose", threshold_rel=tau)), "threshold_rel")
    validate(cfg(influence=InfluenceOptions("dose", threshold_rel=0.99)))
    for g in (0, 3, -1):
        refused(cfg(influence=InfluenceOptions("dose", group_size=g)), "group_size")
    validate(cfg(influence=InfluenceOptions("dose", group_size=2)))
    for d in (
        DiagnosticsOptions(trace_histories=1),
        DiagnosticsOptions(track_end_positions=True),
        DiagnosticsOptions(escape_records=True),
    ):
        refused(
            cfg(influence=InfluenceOptions("dose"), diagnostics=d), "not supported with influence"
        )
    with pytest.raises(U, match="persist_batches must be a bool"):
        InfluenceOptions("dose", persist_batches="yes")  # type: ignore[arg-type]


def _estimates(c: SimulationConfig, v_inf: int = 64) -> tuple[Any, Any]:
    """The effective configuration (budget large) and a function G -> allocation_estimate."""
    eff = validate(c)
    src = c.source
    assert isinstance(src, BeamletSet)
    off = [int(x) for x in src.offsets()]

    def est(g: int) -> int:
        return allocation_estimate(
            c.run.backend, c.run.n_batches, 64, 1, eff, c.run.n_histories, c.run.chunk_histories,
            v_inf, g, off,
        )  # fmt: skip

    return eff, est


def test_beam1_group_size_memory_rule() -> None:
    bs = (bl(-2.0), bl(0.0), bl(2.0))
    base = cfg(beamlets=bs, influence=InfluenceOptions("dose"))
    eff, est = _estimates(base)
    assert eff.influence_group_size == 3 and eff.influence_group_origin == "memory-budget"
    assert est(1) < est(2) < est(3)
    two = cfg(beamlets=bs, influence=InfluenceOptions("dose"), run={"memory_budget_bytes": est(2)})
    e2 = validate(two)
    assert (e2.influence_group_size, e2.influence_group_origin) == (2, "memory-budget")
    assert e2.summary()["influence"]["group_size"] == 2
    one = cfg(beamlets=bs, influence=InfluenceOptions("dose"), run={"memory_budget_bytes": est(1)})
    assert validate(one).influence_group_size == 1
    over = {"memory_budget_bytes": est(1) - 1}
    refused(
        cfg(beamlets=bs, influence=InfluenceOptions("dose"), run=over), "at beamlet group size 1"
    )
    refused(
        cfg(beamlets=bs, influence=InfluenceOptions("dose", group_size=3),
            run={"memory_budget_bytes": est(2)}),
        "influence.group_size (3) needs",
    )  # fmt: skip
    ok = cfg(beamlets=bs, influence=InfluenceOptions("dose", group_size=2),
             run={"memory_budget_bytes": est(2)})  # fmt: skip
    assert validate(ok).influence_group_origin == "requested"


def test_beam1_group_size_scan_is_exhaustive_for_non_monotone_estimates() -> None:
    """warp-cpu, n = [2, 2, 512, 512, 2, 2], B = 2, chunk 1024: L_1 = 512, L_2 = 1024, L_3 = 516,
    so the estimate is not monotone in G; a budget admitting G = 3 but not G = 2 yields G = 3."""
    bs = tuple(bl(float(i), n=n) for i, n in enumerate((2, 2, 512, 512, 2, 2)))
    wc = {"backend": "warp-cpu", "precision": "float32", "n_batches": 2, "chunk_histories": 1024}
    base = cfg(beamlets=bs, influence=InfluenceOptions("dose"), run=wc)
    _, est = _estimates(base)
    assert est(3) < est(2) and est(1) < est(3)
    budget = est(3)  # admits G = 1 and G = 3, not G = 2 (nor G = 4, 5, 6, which are larger)
    assert [g for g in range(1, 7) if est(g) <= budget] == [1, 3]
    c = cfg(beamlets=bs, influence=InfluenceOptions("dose"),
            run=wc | {"memory_budget_bytes": budget})  # fmt: skip
    eff = validate(c)
    assert (eff.influence_group_size, eff.influence_group_origin) == (3, "memory-budget")


def test_beam1_composite_run_without_influence_has_group_one_and_budget_rule() -> None:
    eff = validate(cfg())
    assert eff.influence_group_size is None
    _, est = _estimates(cfg(), 0)
    refused(cfg(run={"memory_budget_bytes": est(1) - 1}), "at beamlet group size 1")


def test_beam1_mixed_projectiles_refused_in_the_toml_form() -> None:
    d = configio.to_dict(cfg())
    d["source"]["beamlets"][1]["projectile"] = "alpha"
    with pytest.raises(ConfigError, match="mixed projectiles"):
        configio.from_dict(d)
    configio.from_dict(configio.to_dict(cfg()))


def test_beam1_ion_mapping_failure_branch_names_the_beamlet(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bs = (bl(e=400.0, label="first"), bl(2.0, e=440.0, label="second"))
    c = cfg(beamlets=bs, projectile=ALPHA, physics={"multiple_scattering": True})
    validate(c)  # positive twin: the mapping is defined
    fake = SimpleNamespace(
        n_uncovered=1, n_band_cells=4, pv_lowest_mev=1.0, r_largest=2.0, statuses=(2,)
    )
    monkeypatch.setattr(config, "pv_eq_table_coverage", lambda *a, **k: fake)
    with pytest.raises(U) as ei:
        validate(c)
    msg = str(ei.value)
    assert "beamlet 0 ('first')" in msg and "beamlet 1 ('second')" in msg
    assert "E0 = 400 MeV" in msg and "E0 = 440 MeV" in msg
    pen = PencilBeamSource(ALPHA, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 400.0)
    cp = cfg(beamlets=(bl(),), source=pen, projectile=ALPHA, run={"n_histories": 40},
             physics={"multiple_scattering": True})  # fmt: skip
    with pytest.raises(U, match="alpha at E0 = 400 MeV is unsupported with multiple scattering"):
        validate(cp)


def test_beam1_nuclear_beamlets_on_python_and_warp_refusal_names_v3_005b() -> None:
    c = cfg(
        physics={"nuclear": True, "nuclear_table_id": "x"},
        run={"backend": "warp-cpu", "precision": "float32"},
    )
    refused(c, "V3-005B")


# ---- the config <-> TOML form -------------------------------------------------------------------


def test_beam1_dict_round_trip_and_toml_errors(tmp_path: Path) -> None:
    c = cfg(beamlets=(bl(label="a", weight=1.5), bl(2.0, energy_sigma_mev=0.5)),
            influence=InfluenceOptions("dose", threshold_rel=0.01, group_size=2,
                                       persist_batches=False))  # fmt: skip
    d = configio.to_dict(c)
    assert d["source"]["type"] == "beamlets" and d["influence"]["group_size"] == 2
    back = configio.from_dict(copy.deepcopy(d))
    assert configio.to_dict(back) == d and back == c
    bad = copy.deepcopy(d)
    bad["influence"]["nope"] = 1
    with pytest.raises(ConfigError, match="influence.nope"):
        configio.from_dict(bad)
    bad = copy.deepcopy(d)
    bad["source"]["beamlets"][0]["energy"] = 1.0
    with pytest.raises(ConfigError, match=r"source.beamlets\[0\].energy"):
        configio.from_dict(bad)
    bad = copy.deepcopy(d)
    bad["source"]["type"] = "spots"
    with pytest.raises(ConfigError, match="source.type"):
        configio.from_dict(bad)
    base = (Path(__file__).parent / "data" / "v3_006" / "proton_box.toml").read_text()
    src_tbl = base.split("[geometry]")[0].split("[source]")[0]
    rest = "[geometry]" + base.split("[geometry]")[1]
    toml = (
        src_tbl
        + '[source]\ntype = "beamlets"\nprojectile = "proton"\n'
        + "[[source.beamlets]]\nposition_mm = [0.0, 0.0, -5.0]\ndirection = [0.0, 0.0, 1.0]\n"
        + 'kinetic_energy_mev = 20.0\nn_histories = 4\nlabel = "a"\n'
        + "[[source.beamlets]]\nposition_mm = [1.0, 0.0, -5.0]\ndirection = [0.0, 0.0, 1.0]\n"
        + "kinetic_energy_mev = 21.0\nn_histories = 8\n\n"
        + '[influence]\ngrid = "dose"\nthreshold_rel = 0.0\n\n'
        + rest
    )
    toml_nohist = re.sub(r"n_histories = 8\n(?=n_batches)", "", toml)
    p = tmp_path / "b.toml"
    p.write_text(toml_nohist)
    # the PUBLIC loader carries the provenance by itself: the omitted key loads as None
    public, _manifest = configio.load_toml(p)
    assert public.run.n_histories is None
    assert configio.load_toml_full(p).config.run.n_histories is None
    assert validate(public).run.n_histories == 12  # effective N = sum of the n_j
    p.write_text(toml.replace("n_histories = 8\nn_batches", "n_histories = 12\nn_batches"))
    assert configio.load_toml(p)[0].run.n_histories == 12
    p.write_text(toml.replace("n_histories = 8\nn_batches", "n_histories = 13\nn_batches"))
    with pytest.raises(ConfigError, match="run.n_histories"):
        configio.load_toml_full(p)
    # a pencil file keeps the key required
    pencil = (Path(__file__).parent / "data" / "v3_006" / "proton_box.toml").read_text()
    p.write_text(re.sub(r"n_histories = \d+\n", "", pencil))
    with pytest.raises(ConfigError, match="n_histories"):
        configio.load_toml_full(p)

    # production path: `ionmc check` / `ionmc run` build Simulation(loaded.config)
    from ionmc.cli import _simulation_of
    from ionmc.simulation import Simulation

    def row(sim: Any) -> tuple[str, str]:
        r = {x.row: (x.cls, x.rule) for x in sim.resolutions}["beamlet-histories"]
        return r

    p.write_text(toml_nohist)
    assert row(_simulation_of(configio.load_toml_full(p))) == ("derived", "sum-of-beamlets")
    assert row(Simulation(configio.load_toml(p)[0])) == ("derived", "sum-of-beamlets")
    p.write_text(toml.replace("n_histories = 8\nn_batches", "n_histories = 12\nn_batches"))
    assert row(_simulation_of(configio.load_toml_full(p))) == ("equal", "sum-equal")
    assert row(Simulation(cfg())) == ("equal", "sum-equal")  # API-built configuration


def _derived(c: SimulationConfig) -> SimulationConfig:
    return dataclasses.replace(c, run=dataclasses.replace(c.run, n_histories=None))


def test_beam1_derived_n_histories_none_effective_and_rows() -> None:
    from ionmc import contract

    c = _derived(cfg())
    assert c.run.n_histories is None
    eff = validate(c)
    assert eff.run.n_histories == 80 and eff.summary()["n_histories"] == 80
    rows = {r.row: (r.cls, r.rule) for r in contract.reconcile(c, eff)}  # config objects
    assert rows["beamlet-histories"] == ("derived", "sum-of-beamlets")
    d = configio.to_dict(c)  # the dicts a bundle persists
    assert d["run"]["n_histories"] is None
    rows = {r.row: (r.cls, r.rule) for r in contract.reconcile(d, eff.summary())}
    assert rows["beamlet-histories"] == ("derived", "sum-of-beamlets")
    # explicit equal -> equal:sum-equal; explicit different -> refused (config and dicts)
    eq = cfg(run={"n_histories": 80})
    rows = {r.row: (r.cls, r.rule) for r in contract.reconcile(eq, validate(eq))}
    assert rows["beamlet-histories"] == ("equal", "sum-equal")
    refused(cfg(run={"n_histories": 81}), "run.n_histories (81)")
    bad = copy.deepcopy(d)
    bad["run"]["n_histories"] = 81
    with pytest.raises(U, match="n_histories"):
        contract.reconcile(bad, eff.summary())
    # Simulation(config) carries it with no side channel
    from ionmc.simulation import Simulation

    sim = Simulation(c)
    assert {r.row: r.rule for r in sim.resolutions}["beamlet-histories"] == "sum-of-beamlets"
    assert not hasattr(sim, "n_histories_origin")


def test_beam1_none_n_histories_with_a_pencil_source_is_refused() -> None:
    from ionmc.sources import PencilBeamSource

    pencil = cfg(
        source=PencilBeamSource(PROTON, (0.0, 0.0, -5.0), (0.0, 0.0, 1.0), 100.0),
        run={"n_histories": None},
    )
    refused(pencil, "run.n_histories is None")
    with pytest.raises(U, match="run.n_histories"):
        validate(pencil)


def test_beam1_n_histories_none_dict_round_trip() -> None:
    c = _derived(cfg(influence=InfluenceOptions("dose", group_size=2)))
    d = configio.to_dict(c)
    assert d["run"]["n_histories"] is None
    back = configio.from_dict(copy.deepcopy(d))
    assert back.run.n_histories is None and back == c and configio.to_dict(back) == d


def test_beam1_derived_n_histories_in_the_bundle_s2_gate(tmp_path: Path) -> None:
    """Bundle-level save/load of a derived-history run: needs transport (S2 flips the flag)."""
    if not ("python" in config.BEAMLET_TRANSPORT_BACKENDS and config.BEAMLET_PERSISTENCE_AVAILABLE):
        pytest.skip("beamlet bundles (S5) are not yet available; the flag unskips")
    import json

    from ionmc import results
    from ionmc.simulation import Simulation

    result = Simulation(_derived(cfg())).run()
    out = tmp_path / "bundle"
    results.save(result, out)
    doc = json.loads((out / "result.json").read_text())
    assert doc["requested"]["run"]["n_histories"] is None
    assert doc["effective"]["n_histories"] == 80
    rows = {r["row"]: r["rule"] for r in doc["effective_extra"]["resolutions"]}
    assert rows["beamlet-histories"] == "sum-of-beamlets"
    results.load(out)  # the reader binds requested null / effective sum


# ---- BEAM-3: offsets and batch membership -------------------------------------------------------


def test_beam3_offsets_and_batch_membership_for_unequal_n() -> None:
    bs = BeamletSet(PROTON, (bl(n=8), bl(1.0, n=16), bl(2.0, n=4)))
    s = bs.offsets()
    assert s.dtype == np.int64 and s.tolist() == [0, 8, 24, 28]
    n_b = 4
    assert bs.histories_per_batch(n_b).tolist() == [2, 4, 1]
    counts = np.zeros((3, n_b), dtype=int)
    for h in range(int(s[-1])):
        counts[bs.beamlet_index(h), h % n_b] += 1
    assert counts.tolist() == [[2] * 4, [4] * 4, [1] * 4]
    h = np.arange(28)
    assert bs.beamlet_index(h).tolist() == [0] * 8 + [1] * 16 + [2] * 4
    for bad in (-1, 28):
        with pytest.raises(ValueError, match="outside"):
            bs.beamlet_index(bad)
    with pytest.raises(ValueError, match="multiple"):
        BeamletSet(PROTON, (bl(n=6),)).histories_per_batch(4)


def test_beam3_canonical_table_and_hash() -> None:
    a = BeamletSet(PROTON, (bl(label="x"), bl(1.0, 100.0)))
    t = a.canonical_table()
    assert t.shape == (2, 12) and t.dtype == np.float64
    assert t[0, 3:6].tolist() == [0.0, 0.0, 1.0] and t[0, 10] == 0.0
    h0 = a.canonical_sha256()
    assert re.fullmatch(r"[0-9a-f]{64}", h0) and h0 == a.canonical_sha256()
    for other in (
        BeamletSet(PROTON, (bl(label="y"), a.beamlets[1])),
        BeamletSet(PROTON, (bl(label="x", weight=0.0), a.beamlets[1])),
        BeamletSet(PROTON, (a.beamlets[1], a.beamlets[0])),
    ):
        assert other.canonical_sha256() != h0
    d = Beamlet((0.0, 0.0, 0.0), (0.0, 0.0, 5.0), 100.0)
    assert d.unit_direction() == (0.0, 0.0, 1.0)
    assert a.unit_direction() is None and a.unit_directions().shape == (2, 3)
    env = a.envelope()
    assert env.e_cap_mev == env.e_hi_mev == 100.0 and env.e_lo_mev == 100.0


# ---- BEAM-2: AST source-consumer audit ---------------------------------------------------------

# Every attribute access ``<expr>.source`` on ANY expression (a name, an attribute chain, a call
# result, a subscript) is a candidate root; there is no vocabulary of base names. Each root
# occurrence is an event ``(module, scope, "<root>:<expression>")`` that the allowlist classifies
# (occurrences whose ``.source`` is not a configuration source, such as a lookup-file "source",
# are classified ``not-a-source``); the attribute reads, aliases and pass-through of a root are
# further events. Subscript reads of a stored document (``doc["requested"]["source"]``) are not
# attribute accesses and are out of scope; the TOML/dict readers of configio are inventoried
# through ``to_dict``/``from_dict`` instead.


Event = tuple[str, str, str]
Coords = dict[Event, list[tuple[int, int]]]
Key = tuple[str, str]  # (module path relative to src/ionmc, qualified scope name)
_FUNC = (ast.FunctionDef, ast.AsyncFunctionDef)


class _Package:
    """Two-phase fixpoint analysis of ``.source`` provenance over a set of modules.

    Phase 1 (``solve``) iterates until nothing new appears: per scope, the root names (``x.source``
    attribute accesses are roots by syntax; parameter and closure roots are added below), the alias
    set (``s = <source value>``, tuple unpacking, transitive), then (i) every call site that
    resolves to a function/method of the package and passes a source value by position or keyword
    makes the receiving PARAMETER a root of the callee scope, (ii) a nested ``def`` that reads an
    enclosing source name makes that name a root of the inner scope (closure), (iii) a scope that
    directly returns a source value (a root, an alias, a parameter root) is marked
    returns-source, so ``x = f(cfg)`` binds ``x`` as an alias in the caller. Phase 2 (``events``)
    emits the per-scope events. Resolvable callees: module-level and enclosing-scope functions by
    name, ``self.m``/``cls.m`` and ``Class.m`` methods, ``Class(...)`` through its ``__init__``,
    ``from .mod import f`` names (also through ``__init__`` re-exports) and ``mod.f`` on an imported
    module. Documented limits: ``obj.method(...)`` on a non-self object, inherited methods,
    dataclass-generated ``__init__``, ``getattr``/dynamic calls, external callees, ``*args``, a
    source returned inside a container, and shadowing of a root name (over-approximated: a
    reassigned root name stays a root) are NOT followed; a pass to an unresolvable callee remains a
    ``<bare>``/``<alias-passed>`` event of the caller, which the allowlist lists as such."""

    def __init__(self, texts: dict[str, str]) -> None:
        self.nodes: dict[Key, list[ast.AST]] = {}
        self.defs: dict[Key, list[ast.AST]] = {}
        self.parents: dict[Key, dict[int, ast.AST]] = {}
        self.binds: dict[Key, list[tuple[ast.AST, ast.AST]]] = {}
        self.dotted: dict[str, str] = {}
        self.imports: dict[str, dict[str, tuple[str | None, str]]] = {}
        self.roots: dict[Key, set[str]] = {}
        self.labels: dict[Key, set[str]] = {}
        self.rets: set[Key] = set()
        self.fields: set[tuple[str, str, str]] = set()  # (module, class, attribute) source fields
        self._inst: dict[Key, dict[str, Key]] = {}
        self._callee: dict[int, tuple[Key, int] | None] = {}
        for mod, text in texts.items():
            parts = list(Path(mod).with_suffix("").parts)
            self.dotted[".".join(parts[:-1] if parts[-1] == "__init__" else parts)] = mod
            self._collect(ast.parse(text), mod, [])
        for (mod, _q), nodes in self.nodes.items():
            imp = self.imports.setdefault(mod, {})
            for n in nodes:
                if isinstance(n, ast.ImportFrom):
                    tgt = self._absolute(mod, n.level, n.module)
                    for a in n.names:
                        imp[a.asname or a.name] = (tgt, a.name)
        for key, nodes in self.nodes.items():
            par: dict[int, ast.AST] = {}
            for n in nodes:
                for ch in ast.iter_child_nodes(n):
                    par[id(ch)] = n
            self.parents[key] = par
            self.binds[key] = [b for n in nodes for b in self._bindings(n)]

    # -- structure -------------------------------------------------------------------------------
    def _collect(self, node: ast.AST, mod: str, qual: list[str]) -> None:
        for ch in ast.iter_child_nodes(node):
            if isinstance(ch, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
                self.defs.setdefault((mod, ".".join([*qual, ch.name])), []).append(ch)
                self._collect(ch, mod, [*qual, ch.name])
            else:
                self.nodes.setdefault((mod, ".".join(qual) or "<module>"), []).append(ch)
                self._collect(ch, mod, qual)

    @staticmethod
    def _bindings(n: ast.AST) -> list[tuple[ast.AST, ast.AST]]:
        """(target, value) element pairs of an assignment (tuple targets paired elementwise)."""
        if not isinstance(n, ast.Assign | ast.AnnAssign | ast.AugAssign | ast.NamedExpr):
            return []
        if n.value is None:
            return []
        pairs: list[tuple[ast.AST, ast.AST]] = []
        for t in n.targets if isinstance(n, ast.Assign) else [n.target]:
            if isinstance(t, ast.Tuple) and isinstance(n.value, ast.Tuple):
                if len(t.elts) == len(n.value.elts):
                    pairs += list(zip(t.elts, n.value.elts, strict=True))
            else:
                pairs.append((t, n.value))
        return pairs

    @staticmethod
    def _absolute(mod: str, level: int, module: str | None) -> str | None:
        if level == 0:
            if module == "ionmc":
                return ""
            return module[len("ionmc.") :] if module and module.startswith("ionmc.") else None
        base = list(Path(mod).parts[:-1])
        base = base[: len(base) - (level - 1)] if level > 1 else base
        return ".".join([*base, *(module.split(".") if module else [])])

    def _is_class(self, key: Key) -> bool:
        return isinstance(self.defs[key][0], ast.ClassDef)

    def _lookup(self, mod: str, name: str, depth: int = 0) -> tuple[str, Key] | None:
        if (mod, name) in self.defs:
            return ("class" if self._is_class((mod, name)) else "func"), (mod, name)
        if depth > 5 or name not in self.imports.get(mod, {}):
            return None
        tgt, orig = self.imports[mod][name]
        if tgt is None:
            return None
        rel = self.dotted.get(tgt)
        if rel is not None:
            hit = self._lookup(rel, orig, depth + 1)
            if hit is not None:
                return hit
        sub = self.dotted.get(f"{tgt}.{orig}" if tgt else orig)
        return ("module", (sub, "")) if sub is not None else None

    def _resolve_name(self, mod: str, q: str, name: str) -> tuple[str, Key] | None:
        parts = [] if q == "<module>" else q.split(".")
        for i in range(len(parts), 0, -1):  # enclosing FUNCTION scopes (class bodies are skipped)
            pre = (mod, ".".join(parts[:i]))
            k = (mod, f"{pre[1]}.{name}")
            if pre in self.defs and isinstance(self.defs[pre][0], _FUNC) and k in self.defs:
                return ("class" if self._is_class(k) else "func"), k
        return self._lookup(mod, name)

    def _enclosing_class(self, mod: str, q: str) -> str | None:
        parts = q.split(".")
        for i in range(len(parts), 0, -1):
            pre = ".".join(parts[:i])
            if (mod, pre) in self.defs and self._is_class((mod, pre)):
                return pre
        return None

    @staticmethod
    def _decorated(d: ast.AST, name: str) -> bool:
        return any(
            isinstance(x, ast.Name) and x.id == name for x in getattr(d, "decorator_list", [])
        )

    def callee(self, call: ast.Call, key: Key) -> tuple[Key, int] | None:
        """Resolve a call to (callee scope, number of leading implicit parameters) or None."""
        if id(call) not in self._callee:
            self._callee[id(call)] = self._resolve_call(call, key)
        return self._callee[id(call)]

    def _method(self, mk: Key, bound: bool) -> tuple[Key, int] | None:
        if mk not in self.defs or self._is_class(mk):
            return None
        d = self.defs[mk][0]
        if self._decorated(d, "staticmethod"):
            return mk, 0
        return mk, 1 if bound or self._decorated(d, "classmethod") else 0

    def _resolve_call(self, call: ast.Call, key: Key) -> tuple[Key, int] | None:
        mod, q = key
        f = call.func
        hit: tuple[str, Key] | None = None
        if isinstance(f, ast.Name):
            hit = self._resolve_name(mod, q, f.id)
            attr = None
        elif isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name):
            attr = f.attr
            if f.value.id in ("self", "cls"):
                cq = self._enclosing_class(mod, q)
                return self._method((mod, f"{cq}.{attr}"), True) if cq else None
            hit = self._resolve_name(mod, q, f.value.id)
        else:
            return None
        if hit is None:
            return None
        kind, k = hit
        if attr is None:
            if kind == "func":
                return k, 0
            ik = (k[0], f"{k[1]}.__init__") if kind == "class" else None
            return (ik, 1) if ik is not None and ik in self.defs else None
        if kind == "class":
            return self._method((k[0], f"{k[1]}.{attr}"), False)
        if kind == "module" and (k[0], attr) in self.defs:
            tgt = (k[0], attr)
            if self._is_class(tgt):
                ik = (k[0], f"{attr}.__init__")
                return (ik, 1) if ik in self.defs else None
            return tgt, 0
        return None

    def _class_of_call(self, call: ast.Call, key: Key) -> Key | None:
        """The package class constructed by ``Class(...)`` / ``mod.Class(...)``, else None."""
        mod, q = key
        f = call.func
        if isinstance(f, ast.Name):
            hit = self._resolve_name(mod, q, f.id)
            return hit[1] if hit is not None and hit[0] == "class" else None
        if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name):
            hit = self._resolve_name(mod, q, f.value.id)
            if hit is not None and hit[0] == "module":
                k = (hit[1][0], f.attr)
                return k if k in self.defs and self._is_class(k) else None
        return None

    def instances(self, key: Key) -> dict[str, Key]:
        """Names of the scope (own, enclosing-function and module-level bindings) bound to an
        instance of a package class (``x = Class(...)``), plus ``self``/``cls`` of a method."""
        if key not in self._inst:
            mod, q = key
            parts = [] if q == "<module>" else q.split(".")
            scopes = [(mod, "<module>")] + [
                (mod, ".".join(parts[:i])) for i in range(1, len(parts))
            ]
            out: dict[str, Key] = {}
            for sk in [*scopes, key]:
                for t, v in self.binds.get(sk, ()):
                    if isinstance(t, ast.Name) and isinstance(v, ast.Call):
                        c = self._class_of_call(v, sk)
                        if c is not None:
                            out[t.id] = c
            cq = self._enclosing_class(mod, q)
            if cq is not None:
                out["self"] = out["cls"] = (mod, cq)
            self._inst[key] = out
        return self._inst[key]

    def _field_owner(self, n: ast.AST, key: Key) -> tuple[Key, str] | None:
        """(class-or-scope, attribute) of ``<name>.<attr>``: the class of the instance name, else
        the scope-local pseudo class of an unresolvable base name (documented limit)."""
        if not (isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name)):
            return None
        cls = self.instances(key).get(n.value.id)
        return (cls or (key[0], f"<scope>:{key[1]}:{n.value.id}")), n.attr

    def fieldread(self, n: ast.AST, key: Key) -> str | None:
        """``Class.attr`` label when ``n`` reads a stored source field, else None."""
        if isinstance(n, ast.Attribute) and isinstance(n.ctx, ast.Load):
            own = self._field_owner(n, key)
            if own is not None and (*own[0], own[1]) in self.fields:
                return f"{own[0][1]}.{own[1]}"
        return None

    def _dataclass_fields(self, ck: Key) -> list[str]:
        return [
            st.target.id
            for st in self.defs[ck][0].body  # type: ignore[attr-defined]
            if isinstance(st, ast.AnnAssign) and isinstance(st.target, ast.Name)
        ]

    # -- provenance --------------------------------------------------------------------------------
    def srcval(self, n: ast.AST, alias: set[str], key: Key) -> bool:
        if isinstance(n, ast.Attribute):
            return n.attr == "source" or self.fieldread(n, key) is not None
        if isinstance(n, ast.Name):
            return n.id in alias
        if isinstance(n, ast.IfExp):
            return self.srcval(n.body, alias, key) or self.srcval(n.orelse, alias, key)
        if isinstance(n, ast.BoolOp):
            return any(self.srcval(v, alias, key) for v in n.values)
        if isinstance(n, ast.Call):
            c = self.callee(n, key)
            return c is not None and c[0] in self.rets
        return False

    def alias_of(self, key: Key) -> set[str]:
        alias = set(self.roots.get(key, ()))
        grew = True
        while (
            grew
        ):  # transitive aliases: src2 = src, s = src2, x = f(cfg.source) with f returning it
            grew = False
            for t, v in self.binds[key]:
                if isinstance(t, ast.Name) and t.id not in alias and self.srcval(v, alias, key):
                    alias.add(t.id)
                    grew = True
        return alias

    def _params(self, d: ast.AST, skip: int) -> tuple[list[str], set[str]]:
        a = d.args  # type: ignore[attr-defined]
        pos = [x.arg for x in a.posonlyargs + a.args][skip:]
        return pos, set(pos) | {x.arg for x in a.kwonlyargs}

    def _add_root(self, key: Key, name: str, label: str) -> bool:
        grew = label not in self.labels.setdefault(key, set())
        self.labels[key].add(label)
        self.roots.setdefault(key, set()).add(name)
        return grew

    def _pass_edges(self, call: ast.Call, alias: set[str], key: Key) -> list[tuple[Key, str]]:
        c = self.callee(call, key)
        if c is None:
            return []
        ck, skip = c
        out: list[tuple[Key, str]] = []
        for d in self.defs[ck]:
            if not isinstance(d, _FUNC):
                continue
            pos, names = self._params(d, skip)
            for i, arg in enumerate(call.args):
                if isinstance(arg, ast.Starred):
                    break
                if i < len(pos) and self.srcval(arg, alias, key):
                    out.append((ck, pos[i]))
            for kw in call.keywords:
                if kw.arg in names and self.srcval(kw.value, alias, key):
                    out.append((ck, str(kw.arg)))
        return out

    def _ctor_fields(self, call: ast.Call, alias: set[str], key: Key) -> bool:
        """A source passed to the generated ``__init__`` of a dataclass-like class (no explicit
        ``__init__``) marks the receiving annotated field as a source field of that class."""
        ck = self._class_of_call(call, key)
        if ck is None or (ck[0], f"{ck[1]}.__init__") in self.defs:
            return False
        names = self._dataclass_fields(ck)
        recv = [
            names[i]
            for i, a in enumerate(call.args)
            if i < len(names) and not isinstance(a, ast.Starred) and self.srcval(a, alias, key)
        ] + [
            str(kw.arg)
            for kw in call.keywords
            if kw.arg in names and self.srcval(kw.value, alias, key)
        ]
        new = {(*ck, r) for r in recv} - self.fields
        self.fields |= new
        return bool(new)

    def solve(self) -> None:
        changed = True
        while changed:
            changed = False
            aliases = {k: self.alias_of(k) for k in self.nodes}
            for key, nodes in self.nodes.items():
                mod, q = key
                al = aliases[key]
                if key not in self.rets and any(
                    isinstance(n, ast.Return)
                    and n.value is not None
                    and self.srcval(n.value, al, key)
                    for n in nodes
                ):
                    self.rets.add(key)
                    changed = True
                for t, v in self.binds[key]:  # object-field provenance: ``obj.x = <source>``
                    own = self._field_owner(t, key)
                    if own is not None and self.srcval(v, al, key):
                        fk = (*own[0], own[1])
                        if fk not in self.fields:
                            self.fields.add(fk)
                            changed = True
                for n in nodes:
                    if isinstance(n, ast.Call):
                        for ck, p in self._pass_edges(n, al, key):
                            changed |= self._add_root(ck, p, f"<param>:{ck[1]}.{p}")
                        changed |= self._ctor_fields(n, al, key)
                pk = (mod, q.rpartition(".")[0])
                if (
                    "." in q
                    and pk in aliases
                    and aliases[pk]
                    and isinstance(self.defs.get(pk, [None])[0], _FUNC)
                    and any(isinstance(d, _FUNC) for d in self.defs.get(key, []))
                ):
                    local = {n.id for n in nodes if isinstance(n, ast.Name) and not _load(n)}
                    for d in self.defs[key]:
                        if isinstance(d, _FUNC):
                            local |= self._params(d, 0)[1]
                    for n in nodes:
                        if (
                            isinstance(n, ast.Name)
                            and _load(n)
                            and n.id in aliases[pk]
                            and n.id not in local
                        ):
                            changed |= self._add_root(key, n.id, f"<closure>:{q}.{n.id}")

    # -- events ------------------------------------------------------------------------------------
    def events(self) -> Coords:
        events: Coords = {}

        def add(n: ast.AST, ev: Event) -> None:
            events.setdefault(ev, []).append((getattr(n, "lineno", 0), getattr(n, "col_offset", 0)))

        for key, nodes in self.nodes.items():
            mod, q = key
            alias = self.alias_of(key)
            parents = self.parents[key]
            for label in self.labels.get(key, ()):
                add(self.defs[key][0], (mod, q, label))
            for n in nodes:
                par = parents.get(id(n))
                base = isinstance(par, ast.Attribute) and par.value is n
                bound = isinstance(par, ast.Assign | ast.AnnAssign | ast.AugAssign) and (
                    par.value is n
                )
                is_src = isinstance(n, ast.Attribute) and n.attr == "source"
                if isinstance(n, ast.Attribute) and self.srcval(n.value, alias, key):
                    add(n, (mod, q, n.attr))
                fld = self.fieldread(n, key)
                if fld is not None:
                    add(n, (mod, q, f"<field>:{fld}"))
                    if not base and not bound and not is_src:
                        add(n, (mod, q, "<alias-passed>"))
                if is_src:
                    add(n, (mod, q, f"<root>:{ast.unparse(n)}"))
                    if not base and not bound:
                        add(n, (mod, q, "<bare>"))
                elif (
                    (isinstance(n, ast.Name) and _load(n) or isinstance(n, ast.Call))
                    and self.srcval(n, alias, key)
                    and not base
                    and not bound
                ):
                    add(n, (mod, q, "<alias-passed>"))
        return events


def _load(n: ast.Name) -> bool:
    return isinstance(n.ctx, ast.Load)


def audit_package(texts: dict[str, str]) -> Coords:
    """Every (module, scope, event) of a ``.source`` root over the modules ``texts`` (relative path
    -> text), each mapped to the (line, column) of EVERY syntactic occurrence; interprocedural
    (call parameters, returns, closures) as documented on ``_Package``."""
    pkg = _Package(texts)
    pkg.solve()
    return pkg.events()


def audit_text(text: str, mod: str = "mod.py") -> Counter[Event]:
    """Counted multiset of the events of ``audit_coords`` (one count per syntactic occurrence)."""
    return Counter({k: len(v) for k, v in audit_coords(text, mod).items()})


def audit_coords(text: str, mod: str = "mod.py") -> Coords:
    """``audit_package`` of a single module: the root occurrence itself (``<root>:<expression
    text>``, stable under line drift), direct attribute reads, reads through aliases bound in the
    same scope (``src = cfg.source``, tuple unpacking, TRANSITIVE aliases; the alias set is a
    fixpoint), pass-through of the source or an alias in any position, by keyword, inside a
    tuple/list/set/dict literal or a return (``<bare>`` / ``<alias-passed>``), and the events that
    propagate into callee scopes (``<param>:``), closures (``<closure>:``) and returns."""
    return audit_package({mod: text})


def audit_readers(
    root: Path, overrides: dict[str, str] | None = None
) -> tuple[Counter[Event], Coords]:
    """Counted inventory and coordinates over every module of ``root`` (``overrides`` replaces the
    analysed text of a module by its relative path: used by the repository-level control)."""
    texts = {
        str(p.relative_to(root)): (overrides or {}).get(str(p.relative_to(root)), p.read_text())
        for p in sorted(root.rglob("*.py"))
    }
    coords = audit_package(texts)
    return Counter({k: len(v) for k, v in coords.items()}), coords


def inventory_diff(counts: Counter[Event], coords: Coords, allow: dict[Event, int]) -> list[str]:
    """Exact multiset comparison; a mismatch names the key and ALL its occurrence coordinates."""
    msgs = []
    for k in sorted(set(counts) | set(allow)):
        if counts.get(k, 0) != allow.get(k, 0):
            msgs.append(
                f"{k}: found {counts.get(k, 0)}, allowlisted {allow.get(k, 0)}; "
                f"(line, col) occurrences {coords.get(k, [])}"
            )
    return msgs


ENV, PER, S1, NOT = "envelope", "per-beamlet", "refused-in-S1", "not-a-source"
_PROJ = "projectile"  # one projectile per set: the same value for pencil and set
ALLOWLIST: dict[Event, tuple[str, int]] = {
    ("config.py", "EffectiveConfig.is_ion", _PROJ): (ENV, 1),
    ("config.py", "EffectiveConfig.ls_effective", _PROJ): (ENV, 1),
    ("config.py", "EffectiveConfig.stopping_identity", _PROJ): (ENV, 2),  # 2 occurrences
    ("config.py", "EffectiveConfig.summary", "<alias-passed>"): (PER, 3),  # 3 occurrences
    ("config.py", "EffectiveConfig.summary", "energy_sigma_mev"): (PER, 1),  # pencil branch only
    ("config.py", "EffectiveConfig.summary", "kinetic_energy_mev"): (PER, 1),
    ("config.py", "EffectiveConfig.summary", "lateral_sigma_mm"): (PER, 1),
    ("config.py", "EffectiveConfig.summary", "position_mm"): (PER, 1),
    ("config.py", "EffectiveConfig.summary", _PROJ): (ENV, 4),  # 4 occurrences
    ("config.py", "_beamlet_set_checks", "<alias-passed>"): (PER, 1),
    ("config.py", "_beamlet_set_checks", "beamlets"): (PER, 1),
    ("config.py", "_beamlet_set_checks", "n_beamlets"): (PER, 1),
    ("config.py", "_beamlet_set_checks", "n_histories_total"): (PER, 1),
    ("config.py", "_check_stopping_domain", "<alias-passed>"): (ENV, 1),  # _source_e_hi
    ("config.py", "_check_stopping_domain", _PROJ): (ENV, 1),
    ("config.py", "_compile_tallies", _PROJ): (ENV, 2),  # 2 occurrences
    ("config.py", "_computed_max_steps", "<alias-passed>"): (ENV, 1),  # _source_e_hi
    ("config.py", "_finish_beamlet_set", "<alias-passed>"): (PER, 1),
    ("config.py", "_finish_beamlet_set", "n_beamlets"): (PER, 1),
    ("config.py", "_finish_beamlet_set", "offsets"): (PER, 1),
    ("config.py", "_influence_checks", "<alias-passed>"): (PER, 1),
    ("config.py", "_influence_checks", "n_beamlets"): (PER, 2),  # 2 occurrences
    ("config.py", "_ion_mapping_data", "<bare>"): (
        PER,
        1,
    ),  # _source_description names the beamlets
    ("config.py", "_ion_mapping_data", _PROJ): (ENV, 1),
    ("config.py", "_nuclear_setup_checks", "<alias-passed>"): (ENV, 1),  # _source_e_hi
    ("config.py", "_nuclear_tables", _PROJ): (ENV, 1),
    ("config.py", "_water_table", _PROJ): (ENV, 1),
    ("config.py", "effective_e_cut_mev", _PROJ): (ENV, 1),
    ("config.py", "effective_end_of_range_model", _PROJ): (ENV, 1),
    ("config.py", "EffectiveConfig.run", "n_histories_total"): (PER, 1),  # derived run.n_histories
    ("config.py", "EffectiveConfig.run", "<alias-passed>"): (PER, 1),
    ("config.py", "_with_derived_histories", "n_histories_total"): (PER, 1),  # run.n_histories None
    ("config.py", "_with_derived_histories", "<bare>"): (PER, 1),
    ("config.py", "validate", "<alias-passed>"): (
        ENV,
        3,
    ),  # 3 occurrences; _source_members / envelope helpers
    ("config.py", "validate", _PROJ): (ENV, 18),  # 18 occurrences
    ("config.py", "validate", "unit_direction"): (PER, 1),  # None with status for a set
    ("configio.py", "to_dict", "<alias-passed>"): (PER, 2),  # 2 occurrences
    ("configio.py", "to_dict", "beamlets"): (PER, 1),
    ("configio.py", "to_dict", _PROJ): (ENV, 2),  # 2 occurrences
    ("results.py", "_lookup_plan", _PROJ): (ENV, 1),
    # V3-007 S1 fix 4: isinstance(result.requested_config.source, BeamletSet) persistence gate
    ("results.py", "build_document", "<bare>"): (S1, 1),  # refused until S5
    ("results.py", "save", "<bare>"): (S1, 1),
    ("transport/reference.py", "_Reference.__init__", _PROJ): (ENV, 4),  # 4 occurrences
    # V3-007 S2: the beamlet-set handle (offsets, directions) and the per-beamlet sampling reads
    ("transport/reference.py", "_Reference.__init__", "<alias-passed>"): (PER, 3),  # 3 occurrences
    # object-field provenance (S2 fix 1): `self.beamlet_set = src0 ...` stores the set; readers of
    # the stored field are per-beamlet (offsets, unit directions, beamlet count, accumulators)
    ("transport/reference.py", "_Reference.__init__", "<field>:_Reference.beamlet_set"): (PER, 3),
    ("transport/reference.py", "_Reference.__init__", "beamlets"): (PER, 1),
    ("transport/reference.py", "_Reference.__init__", "offsets"): (PER, 1),
    ("transport/reference.py", "_Reference._flush_group", "<alias-passed>"): (PER, 1),
    ("transport/reference.py", "_Reference._flush_group", "<field>:_Reference.beamlet_set"): (
        PER,
        2,
    ),
    ("transport/reference.py", "_Reference._flush_group", "n_beamlets"): (PER, 1),
    ("transport/reference.py", "_Reference.run_range", "<alias-passed>"): (
        PER,
        3,
    ),  # `bs is not None`
    ("transport/reference.py", "_Reference.run_range", "<field>:_Reference.beamlet_set"): (PER, 1),
    # a materials citation field `Material.source` stored by the dataclass constructor
    ("materials.py", "Material.with_I", "<field>:Material.source"): (NOT, 1),
    ("transport/reference.py", "_Reference._source_parameters", "<alias-passed>"): (PER, 1),
    ("transport/reference.py", "_Reference._source_parameters", "beamlets"): (PER, 1),
    ("transport/reference.py", "_Reference._source_parameters", "energy_sigma_mev"): (PER, 1),
    ("transport/reference.py", "_Reference._source_parameters", "kinetic_energy_mev"): (PER, 1),
    ("transport/reference.py", "_Reference._source_parameters", "lateral_sigma_mm"): (PER, 1),
    ("transport/reference.py", "_Reference._source_parameters", "position_mm"): (PER, 1),
    ("transport/warp_driver.py", "run_warp_range", "<alias-passed>"): (S1, 1),  # isinstance guard
    ("transport/warp_driver.py", "run_warp_range", "energy_sigma_mev"): (S1, 1),  # S3
    ("transport/warp_driver.py", "run_warp_range", "kinetic_energy_mev"): (S1, 1),
    ("transport/warp_driver.py", "run_warp_range", "lateral_sigma_mm"): (S1, 1),
    ("transport/warp_driver.py", "run_warp_range", "position_mm"): (S1, 1),
    ("transport/warp_driver.py", "run_warp_range", _PROJ): (ENV, 4),  # 4 occurrences
    # V3-007 S1 fix 7: root occurrences (every `.source` access) and tuple pass-through
    ("config.py", "EffectiveConfig.is_ion", "<root>:self.requested.source"): (ENV, 1),
    ("config.py", "EffectiveConfig.ls_effective", "<root>:self.requested.source"): (ENV, 1),
    ("config.py", "EffectiveConfig.run", "<root>:self.requested.source"): (PER, 1),
    ("config.py", "EffectiveConfig.stopping_identity", "<root>:self.requested.source"): (
        ENV,
        2,
    ),  # 2 occurrences
    ("config.py", "EffectiveConfig.summary", "<root>:self.requested.source"): (PER, 1),
    ("config.py", "PhysicsOptions.stopping_source", "<bare>"): (NOT, 1),
    ("config.py", "PhysicsOptions.stopping_source", "<root>:self.selection.source"): (NOT, 1),
    ("config.py", "SimulationConfig.__post_init__", "<bare>"): (ENV, 1),
    ("config.py", "SimulationConfig.__post_init__", "<root>:self.source"): (ENV, 1),
    ("config.py", "_beamlet_set_checks", "<bare>"): (PER, 1),
    ("config.py", "_beamlet_set_checks", "<root>:config.source"): (PER, 1),
    ("config.py", "_check_stopping_domain", "<root>:config.source"): (ENV, 1),
    ("config.py", "_compile_tallies", "<bare>"): (ENV, 1),
    ("config.py", "_compile_tallies", "<root>:config.source"): (ENV, 1),
    ("config.py", "_computed_max_steps", "<bare>"): (ENV, 1),
    ("config.py", "_computed_max_steps", "<root>:config.source"): (ENV, 1),
    ("config.py", "_finish_beamlet_set", "<bare>"): (PER, 1),
    ("config.py", "_finish_beamlet_set", "<root>:cfg.source"): (PER, 1),
    ("config.py", "_influence_checks", "<root>:config.source"): (PER, 1),
    ("config.py", "_ion_mapping_data", "<root>:config.source"): (PER, 2),  # 2 occurrences
    ("config.py", "_nuclear_setup_checks", "<bare>"): (ENV, 1),
    ("config.py", "_nuclear_setup_checks", "<root>:config.source"): (ENV, 1),
    ("config.py", "_nuclear_tables", "<root>:config.source"): (ENV, 1),
    ("config.py", "_water_table", "<root>:config.source"): (ENV, 1),
    ("config.py", "_with_derived_histories", "<root>:config.source"): (PER, 2),  # 2 occurrences
    ("config.py", "effective_e_cut_mev", "<root>:config.source"): (ENV, 1),
    ("config.py", "effective_end_of_range_model", "<root>:config.source"): (ENV, 1),
    ("config.py", "validate", "<bare>"): (ENV, 1),
    ("config.py", "validate", "<root>:config.source"): (ENV, 1),
    ("configio.py", "InputManifest.__post_init__", "<bare>"): (NOT, 3),  # 3 occurrences
    ("configio.py", "InputManifest.__post_init__", "<root>:self.source"): (NOT, 3),  # 3 occurrences
    ("configio.py", "InputManifest.to_dict", "<bare>"): (NOT, 1),
    ("configio.py", "InputManifest.to_dict", "<root>:self.source"): (NOT, 1),
    ("configio.py", "_material_dict", "<bare>"): (NOT, 1),
    ("configio.py", "_material_dict", "<root>:m.source"): (NOT, 1),
    ("configio.py", "to_dict", "<root>:config.source"): (PER, 1),
    ("lookup.py", "LookupTable.provenance", "<bare>"): (NOT, 1),
    ("lookup.py", "LookupTable.provenance", "<root>:self.source"): (NOT, 1),
    ("materials.py", "Material.with_I", "<bare>"): (NOT, 1),
    ("materials.py", "Material.with_I", "<root>:self.source"): (NOT, 1),
    ("physics/stopping.py", "BetheStoppingSource.table", "<bare>"): (NOT, 1),
    ("physics/stopping.py", "BetheStoppingSource.table", "<root>:material.source"): (NOT, 1),
    ("physics/stopping_selection.py", "assert_identity_stable", "<bare>"): (
        NOT,
        2,
    ),  # 2 occurrences
    ("physics/stopping_selection.py", "assert_identity_stable", "<root>:sel.source"): (
        NOT,
        2,
    ),  # 2 occurrences
    ("physics/stopping_selection.py", "check_window", "<bare>"): (NOT, 1),
    ("physics/stopping_selection.py", "check_window", "<root>:sel.source"): (NOT, 1),
    ("physics/stopping_selection.py", "identity_block", "<alias-passed>"): (
        NOT,
        3,
    ),  # 3 occurrences
    ("physics/stopping_selection.py", "identity_block", "<root>:sel.source"): (NOT, 1),
    ("physics/stopping_tabulated.py", "TabulatedStoppingSource.table", "<bare>"): (NOT, 1),
    (
        "physics/stopping_tabulated.py",
        "TabulatedStoppingSource.table",
        "<root>:material.source",
    ): (NOT, 1),
    ("results.py", "_lookup_plan", "<root>:cfg.source"): (ENV, 1),
    ("results.py", "_material_def", "<bare>"): (NOT, 1),
    ("results.py", "_material_def", "<root>:m.source"): (NOT, 1),
    ("results.py", "build_document", "<root>:result.requested_config.source"): (S1, 1),
    ("results.py", "save", "<root>:result.requested_config.source"): (S1, 1),
    ("transport/reference.py", "_Reference.__init__", "<root>:cfg.source"): (
        ENV,
        5,
    ),  # 5 occurrences (4 projectile reads + the beamlet-set handle)
    ("transport/reference.py", "_Reference._source_parameters", "<root>:self.cfg.source"): (PER, 1),
    ("transport/warp_driver.py", "run_warp_range", "<bare>"): (S1, 1),
    ("transport/warp_driver.py", "run_warp_range", "<root>:cfg.source"): (S1, 1),
    # ---- interprocedural fix 9: parameter / closure / return roots of the callees reached.
    # Classification: _source_members reads each member (per-beamlet); _source_e_hi,
    # _source_description and the beamlet summary use the envelope / set-level fields; the
    # physics/stopping_* scopes receive a STOPPING-model `source` (a table provider), a
    # false-positive root of the any-`.source` rule: not-a-source.
    # propagated callee EffectiveConfig._beamlet_source_summary (config.py)
    (
        "config.py",
        "EffectiveConfig._beamlet_source_summary",
        "<param>:EffectiveConfig._beamlet_source_summary.s",
    ): (ENV, 1),
    ("config.py", "EffectiveConfig._beamlet_source_summary", "canonical_sha256"): (ENV, 1),
    ("config.py", "EffectiveConfig._beamlet_source_summary", "envelope"): (ENV, 1),
    ("config.py", "EffectiveConfig._beamlet_source_summary", "n_beamlets"): (ENV, 1),
    ("config.py", "EffectiveConfig._beamlet_source_summary", "n_histories_total"): (ENV, 1),
    ("config.py", "EffectiveConfig._beamlet_source_summary", "projectile"): (ENV, 4),
    # propagated callee _source_description (config.py)
    ("config.py", "_source_description", "<alias-passed>"): (ENV, 1),
    ("config.py", "_source_description", "<param>:_source_description.src"): (ENV, 1),
    ("config.py", "_source_description", "beamlets"): (ENV, 1),
    ("config.py", "_source_description", "kinetic_energy_mev"): (ENV, 1),
    ("config.py", "_source_description", "projectile"): (ENV, 2),
    # propagated callee _source_e_hi (config.py)
    ("config.py", "_source_e_hi", "<alias-passed>"): (ENV, 1),
    ("config.py", "_source_e_hi", "<param>:_source_e_hi.src"): (ENV, 1),
    ("config.py", "_source_e_hi", "energy_sigma_mev"): (ENV, 1),
    ("config.py", "_source_e_hi", "envelope"): (ENV, 1),
    ("config.py", "_source_e_hi", "kinetic_energy_mev"): (ENV, 1),
    # propagated callee _source_members (config.py)
    ("config.py", "_source_members", "<alias-passed>"): (PER, 1),
    ("config.py", "_source_members", "<param>:_source_members.src"): (PER, 1),
    ("config.py", "_source_members", "beamlets"): (PER, 1),
    ("config.py", "_source_members", "energy_sigma_mev"): (PER, 1),
    ("config.py", "_source_members", "kinetic_energy_mev"): (PER, 1),
    # propagated callee Icru90WaterSource.instance_spec (physics/stopping_icru90.py)
    (
        "physics/stopping_icru90.py",
        "Icru90WaterSource.instance_spec",
        "<param>:Icru90WaterSource.instance_spec.self",
    ): (NOT, 1),
    ("physics/stopping_icru90.py", "Icru90WaterSource.instance_spec", "alpha"): (NOT, 7),
    ("physics/stopping_icru90.py", "Icru90WaterSource.instance_spec", "content_sha256"): (NOT, 1),
    ("physics/stopping_icru90.py", "Icru90WaterSource.instance_spec", "name"): (NOT, 1),
    ("physics/stopping_icru90.py", "Icru90WaterSource.instance_spec", "proton"): (NOT, 6),
    ("physics/stopping_icru90.py", "Icru90WaterSource.instance_spec", "start_range_sha256"): (
        NOT,
        1,
    ),
    # propagated callee _bethe_domain (physics/stopping_models.py)
    ("physics/stopping_models.py", "_bethe_domain", "<param>:_bethe_domain.source"): (NOT, 1),
    ("physics/stopping_models.py", "_bethe_domain", "e_max_per_u"): (NOT, 1),
    ("physics/stopping_models.py", "_bethe_domain", "e_min_per_u"): (NOT, 1),
    # propagated callee _bethe_exact (physics/stopping_models.py)
    ("physics/stopping_models.py", "_bethe_exact", "<alias-passed>"): (NOT, 1),
    ("physics/stopping_models.py", "_bethe_exact", "<param>:_bethe_exact.source"): (NOT, 1),
    ("physics/stopping_models.py", "_bethe_exact", "options"): (NOT, 2),
    # propagated callee _bethe_options (physics/stopping_models.py)
    ("physics/stopping_models.py", "_bethe_options", "<param>:_bethe_options.source"): (NOT, 1),
    ("physics/stopping_models.py", "_bethe_options", "name"): (NOT, 1),
    ("physics/stopping_models.py", "_bethe_options", "options"): (NOT, 1),
    ("physics/stopping_models.py", "_bethe_options", "points_per_decade"): (NOT, 1),
    # propagated callee _bethe_spec (physics/stopping_models.py)
    ("physics/stopping_models.py", "_bethe_spec", "<alias-passed>"): (NOT, 2),
    ("physics/stopping_models.py", "_bethe_spec", "<param>:_bethe_spec.source"): (NOT, 1),
    # propagated callee _exact_types (physics/stopping_models.py)
    ("physics/stopping_models.py", "_exact_types", "<alias-passed>"): (NOT, 1),
    ("physics/stopping_models.py", "_exact_types", "<param>:_exact_types.obj"): (NOT, 1),
    # propagated callee _icru_arrays_trusted (physics/stopping_models.py)
    ("physics/stopping_models.py", "_icru_arrays_trusted", "<param>:_icru_arrays_trusted.source"): (
        NOT,
        1,
    ),
    ("physics/stopping_models.py", "_icru_arrays_trusted", "alpha"): (NOT, 1),
    ("physics/stopping_models.py", "_icru_arrays_trusted", "proton"): (NOT, 1),
    # propagated callee _nist_star_instance_spec (physics/stopping_models.py)
    (
        "physics/stopping_models.py",
        "_nist_star_instance_spec",
        "<param>:_nist_star_instance_spec.source",
    ): (NOT, 1),
    ("physics/stopping_models.py", "_nist_star_instance_spec", "e_max_per_u"): (NOT, 1),
    ("physics/stopping_models.py", "_nist_star_instance_spec", "e_min_per_u"): (NOT, 1),
    ("physics/stopping_models.py", "_nist_star_instance_spec", "name"): (NOT, 1),
    ("physics/stopping_models.py", "_nist_star_instance_spec", "points_per_decade"): (NOT, 1),
    ("physics/stopping_models.py", "_nist_star_instance_spec", "star_table"): (NOT, 2),
    # propagated callee _pristine (physics/stopping_models.py)
    ("physics/stopping_models.py", "_pristine", "<alias-passed>"): (NOT, 2),
    ("physics/stopping_models.py", "_pristine", "<param>:_pristine.obj"): (NOT, 1),
    # propagated callee _tabulated_exact (physics/stopping_models.py)
    ("physics/stopping_models.py", "_tabulated_exact", "<alias-passed>"): (NOT, 1),
    ("physics/stopping_models.py", "_tabulated_exact", "<param>:_tabulated_exact.inner"): (NOT, 1),
    ("physics/stopping_models.py", "_tabulated_exact", "metadata"): (NOT, 1),
    # propagated callee spec_of (physics/stopping_models.py)
    ("physics/stopping_models.py", "spec_of", "<alias-passed>"): (NOT, 21),
    ("physics/stopping_models.py", "spec_of", "<param>:spec_of.source"): (NOT, 1),
    ("physics/stopping_models.py", "spec_of", "alpha"): (NOT, 3),
    ("physics/stopping_models.py", "spec_of", "astar"): (NOT, 3),
    ("physics/stopping_models.py", "spec_of", "cache_dir"): (NOT, 2),
    ("physics/stopping_models.py", "spec_of", "instance_spec"): (NOT, 1),
    ("physics/stopping_models.py", "spec_of", "name"): (NOT, 2),
    ("physics/stopping_models.py", "spec_of", "proton"): (NOT, 3),
    ("physics/stopping_models.py", "spec_of", "pstar"): (NOT, 3),
    # propagated callee NistStarCombinedSource.instance_spec (physics/stopping_nist.py)
    (
        "physics/stopping_nist.py",
        "NistStarCombinedSource.instance_spec",
        "<param>:NistStarCombinedSource.instance_spec.self",
    ): (NOT, 1),
    ("physics/stopping_nist.py", "NistStarCombinedSource.instance_spec", "astar"): (NOT, 6),
    ("physics/stopping_nist.py", "NistStarCombinedSource.instance_spec", "name"): (NOT, 1),
    ("physics/stopping_nist.py", "NistStarCombinedSource.instance_spec", "pstar"): (NOT, 5),
    # propagated callee _observed (physics/stopping_selection.py)
    ("physics/stopping_selection.py", "_observed", "<alias-passed>"): (NOT, 1),
    ("physics/stopping_selection.py", "_observed", "<param>:_observed.source"): (NOT, 1),
}
SRC = Path(__file__).resolve().parents[2] / "src" / "ionmc"


def test_beam2_audit_inventory_equals_the_allowlist() -> None:
    """Exact multiset equality per (module, scope, event): an unlisted key, a stale key or a changed
    count (a new occurrence of an already-listed key) fails, naming all occurrence coordinates."""
    counts, coords = audit_readers(SRC)
    allow = {k: n for k, (_c, n) in ALLOWLIST.items()}
    assert inventory_diff(counts, coords, allow) == []
    assert {c for c, _n in ALLOWLIST.values()} <= {ENV, PER, S1, NOT}
    assert all(n >= 1 for n in allow.values())


def test_beam2_inventory_size_is_pinned() -> None:
    """Decision 0046 quotes the inventory size (occurrences counted per key)."""
    assert (len(ALLOWLIST), sum(n for _c, n in ALLOWLIST.values())) == (194, 309)


def test_beam2_a_duplicate_occurrence_in_the_repository_fails_naming_the_scope() -> None:
    """Repository-level control: append a duplicate of an already-allowlisted expression to the
    same scope of one module; the key count rises by one, nothing else changes, and the failure
    message names the scope and the coordinates of every occurrence of the key."""
    rel = "config.py"
    text = (SRC / rel).read_text()
    assert "config.source" in text
    extra = "\n\ndef validate(config):\n    config.source\n"
    counts0, _ = audit_readers(SRC)
    # a second function body of the same scope name `validate`: the scope key is identical
    counts, coords = audit_readers(SRC, {rel: text + extra})
    allow = {k: n for k, (_c, n) in ALLOWLIST.items()}
    msgs = inventory_diff(counts, coords, allow)
    key = (rel, "validate", "<root>:config.source")
    assert counts[key] == counts0[key] + 1
    assert len(msgs) == 2  # the root and the bare pass-through of the same scope
    assert all("validate" in m and "occurrences [(" in m for m in msgs)
    assert any(str(key) in m for m in msgs)


def test_beam2_the_refused_in_s1_readers_cannot_see_a_beamlet_set() -> None:
    """Every reader annotated refused-in-S1 reads a pencil-only attribute; a BeamletSet has none
    of them, and ``validate`` refuses to run a set (the gate) until S2/S3 convert them."""
    b = BeamletSet(PROTON, (bl(),))
    for attr in ("kinetic_energy_mev", "energy_sigma_mev", "lateral_sigma_mm", "position_mm"):
        assert not hasattr(b, attr)
    refused_attrs = {a for (_m, _q, a), (v, _n) in ALLOWLIST.items() if v == S1}
    pencil_only = {"kinetic_energy_mev", "energy_sigma_mev", "lateral_sigma_mm", "position_mm"}
    assert {a for a in refused_attrs if not a.startswith("<root>:")} == pencil_only | {
        "<alias-passed>",
        "<bare>",
    }


def test_beam2_the_audit_finds_an_unlisted_reader(tmp_path: Path) -> None:
    (tmp_path / "mod.py").write_text(
        "def f(cfg):\n    src = cfg.source\n    return src.kinetic_energy_mev\n"
        "def g(config):\n    return h(config.source)\n"
        "def k(self):\n    a, b = self.requested.source, 1\n    return a.position_mm, b\n"
        "def other(x):\n    return x.source.kinetic_energy_mev\n"
    )
    ev, _ = audit_readers(tmp_path)
    assert set(ev) == {
        ("mod.py", "f", "<root>:cfg.source"),
        ("mod.py", "f", "kinetic_energy_mev"),
        ("mod.py", "g", "<bare>"),
        ("mod.py", "g", "<root>:config.source"),
        ("mod.py", "k", "<bare>"),
        ("mod.py", "k", "<root>:self.requested.source"),
        ("mod.py", "k", "position_mm"),
        ("mod.py", "other", "<root>:x.source"),  # any `.source` is a candidate root
        ("mod.py", "other", "kinetic_energy_mev"),
    }


def test_beam2_envelope_consumers_use_extremes_not_one_beamlet() -> None:
    """Channel/nuclear/stopping-window bounds read the envelope: the helpers return the max E+6s."""
    two = BeamletSet(PROTON, (bl(e=50.0), bl(2.0, e=140.0, energy_sigma_mev=1.0)))
    assert config._source_e_hi(two) == 146.0
    assert config._source_e_hi(PencilBeamSource(PROTON, (0, 0, 0), (0, 0, 1), 50.0, 1.0)) == 56.0


def test_beam2_every_source_attribute_access_is_a_candidate_root() -> None:
    """Synthetic negative controls of the rule: any expression may precede ``.source``; tuple, list
    and return pass-through produce an event (no suppression)."""
    cases = {
        "def f(c):\n    return c.source.projectile\n": {
            ("mod.py", "f", "<root>:c.source"),
            ("mod.py", "f", "projectile"),
        },
        "def f(result):\n    return result.effective_config.source\n": {
            ("mod.py", "f", "<root>:result.effective_config.source"),
            ("mod.py", "f", "<bare>"),
        },
        "def f():\n    return foo().source.x\n": {
            ("mod.py", "f", "<root>:foo().source"),
            ("mod.py", "f", "x"),
        },
        "def f(cfg):\n    consume((cfg.source,))\n": {
            ("mod.py", "f", "<root>:cfg.source"),
            ("mod.py", "f", "<bare>"),
        },
        "def f(cfg):\n    return cfg.source, 1\n": {
            ("mod.py", "f", "<root>:cfg.source"),
            ("mod.py", "f", "<bare>"),
        },
        "def f(cfg):\n    return [cfg.source]\n": {
            ("mod.py", "f", "<root>:cfg.source"),
            ("mod.py", "f", "<bare>"),
        },
        "def f(cfg):\n    (a, b) = (cfg.source, 1)\n    return a.position_mm\n": {
            ("mod.py", "f", "<root>:cfg.source"),
            ("mod.py", "f", "<bare>"),
            ("mod.py", "f", "position_mm"),
        },
        "def f(xs):\n    return xs[0].source.label\n": {
            ("mod.py", "f", "<root>:xs[0].source"),
            ("mod.py", "f", "label"),
        },
    }
    for text, want in cases.items():
        assert set(audit_text(text)) == want, text


def test_beam2_the_analyser_finds_every_root_alias_and_keyword_form() -> None:
    """Negative controls on synthetic text: (i) a result-rooted read, (ii) a chained alias, (iii) a
    keyword pass, plus eff/loaded roots; ``x.source`` is a root too, a subscripted document read
    is not an attribute access."""
    text = (
        "def a(result):\n    return result.requested_config.source.kinetic_energy_mev\n"
        "def b(cfg):\n    src = cfg.source\n    src2 = src\n    s = src2\n"
        "    return src2.position_mm, s.weight\n"
        "def c(cfg):\n    return f(source=cfg.source)\n"
        "def d(eff, loaded):\n    return eff.source.projectile, eff.requested.source.label,"
        " loaded.config.source.direction\n"
        "def e(self):\n    x, y = self.config.source, 0\n    z = x\n    return g(z), y\n"
        "def n(x, doc):\n    return x.source.kinetic_energy_mev, doc['requested']['source']\n"
    )
    assert set(audit_text(text)) == {
        ("mod.py", "a", "<root>:result.requested_config.source"),
        ("mod.py", "a", "kinetic_energy_mev"),
        ("mod.py", "b", "<root>:cfg.source"),
        ("mod.py", "b", "position_mm"),
        ("mod.py", "b", "weight"),
        ("mod.py", "c", "<bare>"),
        ("mod.py", "c", "<root>:cfg.source"),
        ("mod.py", "d", "<root>:eff.source"),
        ("mod.py", "d", "<root>:eff.requested.source"),
        ("mod.py", "d", "<root>:loaded.config.source"),
        ("mod.py", "d", "projectile"),
        ("mod.py", "d", "label"),
        ("mod.py", "d", "direction"),
        ("mod.py", "e", "<alias-passed>"),
        ("mod.py", "e", "<bare>"),
        ("mod.py", "e", "<root>:self.config.source"),
        ("mod.py", "n", "<root>:x.source"),  # a non-config `.source` is still a root: the
        ("mod.py", "n", "kinetic_energy_mev"),  # allowlist classifies it not-a-source
    }


def test_beam2_repeated_occurrences_are_counted_not_collapsed() -> None:
    """Synthetic negative controls of the counted multiset: a repeated root read, an alias read
    through two attributes and a duplicated pass-through each add one count per occurrence."""
    a = audit_text(
        "def f(cfg):\n    x = cfg.source.kinetic_energy_mev\n"
        "    y = cfg.source.kinetic_energy_mev\n"
    )
    assert a[("mod.py", "f", "kinetic_energy_mev")] == 2
    assert a[("mod.py", "f", "<root>:cfg.source")] == 2
    b = audit_text("def f(cfg):\n    src = cfg.source\n    return src.position_mm, src.energy\n")
    assert b[("mod.py", "f", "position_mm")] == 1 and b[("mod.py", "f", "energy")] == 1
    assert b[("mod.py", "f", "<root>:cfg.source")] == 1
    c = audit_text("def f(cfg):\n    consume(cfg.source)\n    consume(cfg.source)\n")
    assert c[("mod.py", "f", "<bare>")] == 2 and c[("mod.py", "f", "<root>:cfg.source")] == 2
    coords = audit_coords("def f(cfg):\n    consume(cfg.source)\n    consume(cfg.source)\n")
    assert coords[("mod.py", "f", "<bare>")] == [(2, 12), (3, 12)]
    # the comparison names the key and every coordinate when a count differs from the allowlist
    (msg,) = inventory_diff(
        c, coords, {k: v for k, v in c.items() if k[2] != "<bare>"} | {("mod.py", "f", "<bare>"): 1}
    )
    assert "found 2, allowlisted 1" in msg and "[(2, 12), (3, 12)]" in msg


def test_beam2_source_provenance_propagates_through_calls_returns_and_closures() -> None:
    """Interprocedural negative controls (fix 9) on synthetic text: each form makes the read in the
    CALLEE / closure / receiving alias an event, so a new singular reader there cannot hide."""
    # (a) positional pass; the callee read is an event under the callee
    a = audit_text(
        "def f(cfg):\n    return consume(cfg.source)\n"
        "def consume(src):\n    return src.kinetic_energy_mev\n"
    )
    assert a[("mod.py", "consume", "kinetic_energy_mev")] == 1
    assert a[("mod.py", "consume", "<param>:consume.src")] == 1
    assert a[("mod.py", "f", "<bare>")] == 1
    # (b) keyword pass, read through the keyword parameter
    b = audit_text(
        "def f(cfg):\n    return g(source=cfg.source)\n"
        "def g(*, source):\n    return source.position_mm\n"
    )
    assert b[("mod.py", "g", "position_mm")] == 1 and b[("mod.py", "g", "<param>:g.source")] == 1
    # (c) method pass through self.<m>; `self` is skipped in the positional mapping
    c = audit_text(
        "class K:\n    def f(self, cfg):\n        return self.h(cfg.source)\n"
        "    def h(self, s):\n        return s.energy_sigma_mev\n"
    )
    assert c[("mod.py", "K.h", "energy_sigma_mev")] == 1
    assert c[("mod.py", "K.h", "<param>:K.h.s")] == 1
    assert ("mod.py", "K.h", "<param>:K.h.self") not in c
    # (d) a closure reads an enclosing alias
    d = audit_text(
        "def outer(cfg):\n    s = cfg.source\n"
        "    def inner():\n        return s.position_mm\n    return inner\n"
    )
    assert d[("mod.py", "outer.inner", "position_mm")] == 1
    assert d[("mod.py", "outer.inner", "<closure>:outer.inner.s")] == 1
    # (e) return propagation: the receiving name is an alias in the caller
    e = audit_text(
        "def f(c):\n    return c.source\ndef g(cfg):\n    x = f(cfg)\n    return x.energy\n"
    )
    assert e[("mod.py", "g", "energy")] == 1 and e[("mod.py", "f", "<bare>")] == 1
    # (f) an unresolvable callee stays a bare pass-through of the caller; nothing propagates
    f = audit_text("def f(cfg):\n    return external(cfg.source)\n")
    assert set(f) == {("mod.py", "f", "<root>:cfg.source"), ("mod.py", "f", "<bare>")}
    # a chain: the parameter root propagates onward to a second callee (fixpoint)
    h = audit_text(
        "def f(cfg):\n    return p(cfg.source)\ndef p(a):\n    return q(a)\n"
        "def q(b):\n    return b.weight\n"
    )
    assert h[("mod.py", "q", "weight")] == 1 and h[("mod.py", "p", "<alias-passed>")] == 1


def test_beam2_provenance_crosses_modules_and_constructors() -> None:
    """Two-module synthetic package: the callee lives in another module (``from .b import``, a
    re-export through ``__init__``, ``mod.f`` on an imported module and ``Class(...)`` through its
    ``__init__``); the callee-side reads are events keyed by the callee module."""
    pkg = {
        "a.py": (
            "from .b import consume\nfrom . import b\nfrom .b import Holder\n"
            "def run(cfg):\n    consume(cfg.source)\n    b.other(cfg.source)\n"
            "    return Holder(cfg.source)\n"
        ),
        "b.py": (
            "def consume(src):\n    return src.kinetic_energy_mev\n"
            "def other(z):\n    return z.position_mm\n"
            "class Holder:\n    def __init__(self, s):\n        self.e = s.direction\n"
        ),
        "sub/__init__.py": "from ..b import consume\n",
        "sub/c.py": "from . import consume\ndef run(cfg):\n    return consume(cfg.source)\n",
    }
    ev = set(audit_package(pkg))
    assert ("b.py", "consume", "kinetic_energy_mev") in ev
    assert ("b.py", "other", "position_mm") in ev
    assert ("b.py", "Holder.__init__", "direction") in ev
    assert ("b.py", "consume", "<param>:consume.src") in ev
    assert ("b.py", "Holder.__init__", "<param>:Holder.__init__.s") in ev
    # nothing of b.py is an event without the caller passing a source
    assert not any(k[0] == "b.py" for k in audit_package({"b.py": pkg["b.py"]}))


def test_beam2_the_audit_runtime_over_the_package_is_small() -> None:
    import time

    t0 = time.perf_counter()
    audit_readers(SRC)
    assert time.perf_counter() - t0 < 5.0


def test_beam2_object_field_provenance_controls() -> None:
    """S2 fix 1 (REVIEW-b795ade2): a source stored in an object field is a source field; a later
    read through ``self`` / an instance name, an alias of it, or a pass of it is an event."""
    cls = "class Sim:\n    def __init__(self, cfg):\n        self.src = cfg.source\n"
    # (a) another method reads the stored field
    a = audit_text(cls + "    def run(self):\n        return self.src.kinetic_energy_mev\n")
    assert a[("mod.py", "Sim.run", "<field>:Sim.src")] == 1
    assert a[("mod.py", "Sim.run", "kinetic_energy_mev")] == 1
    # (b) the same through an alias of the field
    b = audit_text(cls + "    def run(self):\n        s = self.src\n        return s.position_mm\n")
    assert b[("mod.py", "Sim.run", "position_mm")] == 1
    assert b[("mod.py", "Sim.run", "<field>:Sim.src")] == 1
    # (c) a pass of the field reaches the callee parameter
    c = audit_text(
        cls + "    def run(self):\n        return self.helper(self.src)\n"
        "    def helper(self, s):\n        return s.weight\n"
    )
    assert c[("mod.py", "Sim.helper", "weight")] == 1
    assert c[("mod.py", "Sim.helper", "<param>:Sim.helper.s")] == 1
    assert c[("mod.py", "Sim.run", "<alias-passed>")] == 1
    # (d) a module-level instance read from the module and from a function
    d = audit_text(cls + "sim = Sim(cfg)\nx = sim.src.energy\ndef f():\n    return sim.src.angle\n")
    assert d[("mod.py", "<module>", "energy")] == 1 and d[("mod.py", "f", "angle")] == 1
    assert (
        d[("mod.py", "<module>", "<field>:Sim.src")] == 1
        and d[("mod.py", "f", "<field>:Sim.src")] == 1
    )
    # (e) a dataclass constructed with the source (keyword and positional) and read later
    dc = "class Holder:\n    source: object\n    n: int = 0\n"
    e = audit_text(
        dc + "def f(cfg):\n    h = Holder(source=cfg.source)\n    return h.source.direction\n"
    )
    assert e[("mod.py", "f", "direction")] == 1 and e[("mod.py", "f", "<field>:Holder.source")] == 1
    e2 = audit_text(dc + "def f(cfg):\n    h = Holder(cfg.source)\n    return h.source.direction\n")
    assert e2[("mod.py", "f", "direction")] == 1
    # an attribute store on a plain object, augmented and tuple targets
    g = audit_text("def f(cfg, o):\n    o.x, o.y = cfg.source, 1\n    return o.x.energy, o.y.z\n")
    assert g[("mod.py", "f", "energy")] == 1 and ("mod.py", "f", "z") not in g
    # (f) specificity: a stored NON-source field read later produces no event
    f = audit_text(
        "class Sim:\n    def __init__(self, cfg):\n        self.src = cfg.source\n"
        "        self.k = 5\n        self.name = cfg.name\n"
        "    def run(self):\n        return self.k.real, self.name.upper\n"
    )
    assert set(f) == {("mod.py", "Sim.__init__", "<root>:cfg.source")}
