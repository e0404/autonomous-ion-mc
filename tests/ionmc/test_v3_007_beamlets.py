"""Rows BEAM-1, BEAM-2, BEAM-3 of plan v3-007 (step S1): the beamlet source model, its fail-closed
validation rules (each provoked separately with its positive twin), the source-consumer audit and
the offsets / batch membership. Only ``validate`` runs here; no transport. Not a clinical claim.

``BEAMLET_TRANSPORT_AVAILABLE`` is the temporary S1 gate (V3-007 S2 removes it); the tests of the
rules switch it on, one test pins the gate itself.
"""

from __future__ import annotations

import ast
import copy
import re
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
        monkeypatch.setattr(config, "BEAMLET_TRANSPORT_AVAILABLE", True)


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


def test_beam1_gate_beamlet_transport_not_yet_available() -> None:
    """The temporary S1 rule: every beamlet run is refused until S2 (TODO: removed by S2)."""
    assert config.BEAMLET_TRANSPORT_AVAILABLE is False
    refused(cfg(), "beamlet transport not yet available (V3-007 S2)")
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
    loaded = configio.load_toml_full(p)
    assert loaded.n_histories_derived and loaded.config.run.n_histories == 12
    p.write_text(toml.replace("n_histories = 8\nn_batches", "n_histories = 12\nn_batches"))
    loaded = configio.load_toml_full(p)
    assert not loaded.n_histories_derived
    p.write_text(toml.replace("n_histories = 8\nn_batches", "n_histories = 13\nn_batches"))
    with pytest.raises(ConfigError, match="run.n_histories"):
        configio.load_toml_full(p)


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

ROOTS = {"config", "cfg", "self.requested", "self.config", "eff.requested", "e.requested"}


def audit_readers(root: Path) -> set[tuple[str, str, str]]:
    """Every (module, scope, attribute) read from ``<config>.source``: direct accesses, accesses
    through aliases bound in the same scope (``src = cfg.source``, tuple unpacking included) and
    pass-through of the source or an alias (``<bare>`` / ``<alias-passed>``)."""
    events: set[tuple[str, str, str]] = set()

    def is_src(n: ast.AST) -> bool:
        return isinstance(n, ast.Attribute) and n.attr == "source" and ast.unparse(n.value) in ROOTS

    def collect(node: ast.AST, qual: list[str], out: dict[str, list[ast.AST]]) -> None:
        for ch in ast.iter_child_nodes(node):
            if isinstance(ch, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
                collect(ch, [*qual, ch.name], out)
            else:
                out.setdefault(".".join(qual) or "<module>", []).append(ch)
                collect(ch, qual, out)

    for p in sorted(root.rglob("*.py")):
        mod = str(p.relative_to(root))
        scopes: dict[str, list[ast.AST]] = {}
        collect(ast.parse(p.read_text()), [], scopes)
        for q, nodes in scopes.items():
            alias: set[str] = set()
            parents: dict[int, ast.AST] = {}
            for n in nodes:
                for ch in ast.iter_child_nodes(n):
                    parents[id(ch)] = n
                if isinstance(n, ast.Assign | ast.AnnAssign | ast.NamedExpr):
                    tg = n.targets if isinstance(n, ast.Assign) else [n.target]
                    val = n.value
                    for t in tg:
                        if (
                            isinstance(t, ast.Tuple)
                            and isinstance(val, ast.Tuple)
                            and len(t.elts) == len(val.elts)
                        ):
                            for a, b in zip(t.elts, val.elts, strict=True):
                                if val is not None and is_src(b) and isinstance(a, ast.Name):
                                    alias.add(a.id)
                        elif val is not None and is_src(val) and isinstance(t, ast.Name):
                            alias.add(t.id)
            for n in nodes:
                par = parents.get(id(n))
                if isinstance(n, ast.Attribute):
                    if is_src(n.value) or (isinstance(n.value, ast.Name) and n.value.id in alias):
                        events.add((mod, q, n.attr))
                if (
                    is_src(n)
                    and not (isinstance(par, ast.Attribute) and par.value is n)
                    and not (isinstance(par, ast.Assign | ast.AnnAssign) and par.value is n)
                    and not isinstance(par, ast.Tuple)
                ):
                    events.add((mod, q, "<bare>"))
                if (
                    isinstance(n, ast.Name)
                    and n.id in alias
                    and isinstance(n.ctx, ast.Load)
                    and not (isinstance(par, ast.Attribute) and par.value is n)
                ):
                    events.add((mod, q, "<alias-passed>"))
    return events


ENV, PER, S1, NOT = "envelope", "per-beamlet", "refused-in-S1", "not-a-source"
_PROJ = "projectile"  # one projectile per set: the same value for pencil and set
ALLOWLIST: dict[tuple[str, str, str], str] = {
    ("config.py", "EffectiveConfig.is_ion", _PROJ): ENV,
    ("config.py", "EffectiveConfig.ls_effective", _PROJ): ENV,
    ("config.py", "EffectiveConfig.stopping_identity", _PROJ): ENV,
    ("config.py", "EffectiveConfig.summary", "<alias-passed>"): PER,
    ("config.py", "EffectiveConfig.summary", "energy_sigma_mev"): PER,  # pencil branch only
    ("config.py", "EffectiveConfig.summary", "kinetic_energy_mev"): PER,
    ("config.py", "EffectiveConfig.summary", "lateral_sigma_mm"): PER,
    ("config.py", "EffectiveConfig.summary", "position_mm"): PER,
    ("config.py", "EffectiveConfig.summary", _PROJ): ENV,
    ("config.py", "_beamlet_set_checks", "<alias-passed>"): PER,
    ("config.py", "_beamlet_set_checks", "beamlets"): PER,
    ("config.py", "_beamlet_set_checks", "n_beamlets"): PER,
    ("config.py", "_beamlet_set_checks", "n_histories_total"): PER,
    ("config.py", "_check_stopping_domain", "<alias-passed>"): ENV,  # _source_e_hi
    ("config.py", "_check_stopping_domain", _PROJ): ENV,
    ("config.py", "_compile_tallies", _PROJ): ENV,
    ("config.py", "_computed_max_steps", "<alias-passed>"): ENV,  # _source_e_hi
    ("config.py", "_finish_beamlet_set", "<alias-passed>"): PER,
    ("config.py", "_finish_beamlet_set", "n_beamlets"): PER,
    ("config.py", "_finish_beamlet_set", "offsets"): PER,
    ("config.py", "_influence_checks", "<alias-passed>"): PER,
    ("config.py", "_influence_checks", "n_beamlets"): PER,
    ("config.py", "_ion_mapping_data", "<bare>"): PER,  # _source_description names the beamlets
    ("config.py", "_ion_mapping_data", _PROJ): ENV,
    ("config.py", "_nuclear_setup_checks", "<alias-passed>"): ENV,  # _source_e_hi
    ("config.py", "_nuclear_tables", _PROJ): ENV,
    ("config.py", "_water_table", _PROJ): ENV,
    ("config.py", "effective_e_cut_mev", _PROJ): ENV,
    ("config.py", "effective_end_of_range_model", _PROJ): ENV,
    ("config.py", "validate", "<alias-passed>"): ENV,  # _source_members / envelope helpers
    ("config.py", "validate", _PROJ): ENV,
    ("config.py", "validate", "unit_direction"): PER,  # None with status for a set
    ("configio.py", "to_dict", "<alias-passed>"): PER,
    ("configio.py", "to_dict", "beamlets"): PER,
    ("configio.py", "to_dict", _PROJ): ENV,
    ("results.py", "_lookup_plan", _PROJ): ENV,
    ("transport/reference.py", "_Reference.__init__", _PROJ): ENV,
    ("transport/reference.py", "_Reference._history", "<alias-passed>"): S1,  # isinstance guard
    ("transport/reference.py", "_Reference._history", "energy_sigma_mev"): S1,  # S2
    ("transport/reference.py", "_Reference._history", "kinetic_energy_mev"): S1,
    ("transport/reference.py", "_Reference._history", "lateral_sigma_mm"): S1,
    ("transport/reference.py", "_Reference._history", "position_mm"): S1,
    ("transport/warp_driver.py", "run_warp_range", "<alias-passed>"): S1,  # isinstance guard
    ("transport/warp_driver.py", "run_warp_range", "energy_sigma_mev"): S1,  # S3
    ("transport/warp_driver.py", "run_warp_range", "kinetic_energy_mev"): S1,
    ("transport/warp_driver.py", "run_warp_range", "lateral_sigma_mm"): S1,
    ("transport/warp_driver.py", "run_warp_range", "position_mm"): S1,
    ("transport/warp_driver.py", "run_warp_range", _PROJ): ENV,
}
SRC = Path(__file__).resolve().parents[2] / "src" / "ionmc"


def test_beam2_audit_inventory_equals_the_allowlist() -> None:
    events = audit_readers(SRC)
    assert events - set(ALLOWLIST) == set(), sorted(events - set(ALLOWLIST))  # unlisted reader
    assert set(ALLOWLIST) - events == set(), sorted(set(ALLOWLIST) - events)  # stale entry
    assert set(ALLOWLIST.values()) <= {ENV, PER, S1, NOT}


def test_beam2_the_refused_in_s1_readers_cannot_see_a_beamlet_set() -> None:
    """Every reader annotated refused-in-S1 reads a pencil-only attribute; a BeamletSet has none
    of them, and ``validate`` refuses to run a set (the gate) until S2/S3 convert them."""
    b = BeamletSet(PROTON, (bl(),))
    for attr in ("kinetic_energy_mev", "energy_sigma_mev", "lateral_sigma_mm", "position_mm"):
        assert not hasattr(b, attr)
    refused_attrs = {a for (_m, _q, a), v in ALLOWLIST.items() if v == S1}
    pencil_only = {"kinetic_energy_mev", "energy_sigma_mev", "lateral_sigma_mm", "position_mm"}
    assert refused_attrs == pencil_only | {"<alias-passed>"}


def test_beam2_the_audit_finds_an_unlisted_reader(tmp_path: Path) -> None:
    (tmp_path / "mod.py").write_text(
        "def f(cfg):\n    src = cfg.source\n    return src.kinetic_energy_mev\n"
        "def g(config):\n    return h(config.source)\n"
        "def k(self):\n    a, b = self.requested.source, 1\n    return a.position_mm, b\n"
        "def other(x):\n    return x.source.kinetic_energy_mev\n"
    )
    ev = audit_readers(tmp_path)
    assert ev == {
        ("mod.py", "f", "kinetic_energy_mev"),
        ("mod.py", "g", "<bare>"),
        ("mod.py", "k", "position_mm"),
    }


def test_beam2_envelope_consumers_use_extremes_not_one_beamlet() -> None:
    """Channel/nuclear/stopping-window bounds read the envelope: the helpers return the max E+6s."""
    two = BeamletSet(PROTON, (bl(e=50.0), bl(2.0, e=140.0, energy_sigma_mev=1.0)))
    assert config._source_e_hi(two) == 146.0
    assert config._source_e_hi(PencilBeamSource(PROTON, (0, 0, 0), (0, 0, 1), 50.0, 1.0)) == 56.0
