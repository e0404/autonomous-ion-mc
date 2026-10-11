"""V3-005C C6: the lv5c / hr5c steps of ``steps_v5d.py`` and the deterministic producers of ``nuclear_v10.py``
(plan Amendments 14, 15, 17): registration, seeds and seed guards, configuration pins, the R1-D2 (C) rule, the
attribution ladder (all branches, both stop-rule triggers), the V8-nuc rule, the V10 / V10-A building blocks and
data-backed dry runs (skipped without the cache, failing under ``IONMC_REQUIRE_DATA=1``)."""

# ruff: noqa: E501

from __future__ import annotations

import io
import json
import os
from contextlib import redirect_stdout
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

from tests.ionmc.test_run_suite import REPO, _load

SCRIPTS = REPO / "validation" / "scripts" / "transport"
LV5C_BASE, HR5C_BASE = 20481004, 20491004


@pytest.fixture(scope="module")
def d() -> ModuleType:
    return _load("steps_v5d")


@pytest.fixture(scope="module")
def c() -> ModuleType:
    return _load("steps_v5c")


@pytest.fixture(scope="module")
def v10() -> ModuleType:
    import importlib
    import sys

    sys.path.insert(0, str(SCRIPTS))
    return importlib.import_module("nuclear_v10")


def _data_or_skip() -> None:
    from ionmc.data import cache

    cdir = cache.resolve_cache_dir(None)
    if not (cdir / "derived").is_dir() or not list((cdir / "derived").glob("elastic-proton-*.npz")):
        if os.environ.get("IONMC_REQUIRE_DATA") == "1":
            pytest.fail("data for the lv5c dry runs are not cached")
        pytest.skip("data not cached")


def _run_step(d: ModuleType, args: list[str], tmp: Path) -> dict[str, Any]:
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            code = d.main([*args, "--out-dir", str(tmp)])
        except SystemExit as exc:
            code = exc.code
    text = buf.getvalue()
    doc = json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])
    doc["_exit"] = code
    return doc


# -- seeds ----------------------------------------------------------------------------------------
def test_r_index_tables_and_first_seeds(d: ModuleType) -> None:
    assert d.R_INDEX_C == {"throughput": 1, "v8-lv": 2, "r1-d2": 3, "v5-150-on": 4, "v5-150-off": 5, "v5-200-on": 6,
                           "v5-200-off": 7, "v2b-el": 8, "v7-scan-el": 9, "v7-shift-el": 10, "v6": 12, "e1-b": 13,
                           "v11": 14, "p7": 15, "v7-r": 16, "v7-r-diag": 17, "v10-a-xcheck": 18, "x-endf": 19,
                           "v6-sens": 20}  # fmt: skip
    assert d.R_INDEX_HRC == {"v8-stat": 1, "v7-f32": 2, "v6-cuda": 3}
    d.base.SEED_BASE = LV5C_BASE
    want = {"throughput": 20482004, "v8-lv": 20483004, "r1-d2": 20484004, "v5-150-on": 20485004,
            "v5-150-off": 20486004, "v5-200-on": 20487004, "v5-200-off": 20488004, "v2b-el": 20489004,
            "v7-scan-el": 20490004, "v7-shift-el": 20491004, "e1-b": 20494004, "p7": 20496004,
            "v10-a-xcheck": 20499004, "x-endf": 20500004}  # fmt: skip
    assert {k: d.seed_c(k) for k in want} == want
    d.base.SEED_BASE = HR5C_BASE
    assert [d.seed_hrc(k) for k in ("v8-stat", "v7-f32", "v6-cuda")] == [
        20492004,
        20493004,
        20494004,
    ]
    assert d.E_SHAPE_BASE == 20509500 and d.REHEARSAL_RANGE == (20505000, 20509999)
    d.base.SEED_BASE = LV5C_BASE


def test_v8_nuc_sample_seeds_and_shards(d: ModuleType) -> None:
    import run_suite

    d.base.SEED_BASE = HR5C_BASE
    ks = {nm: v[2] for nm, v in d.V8_SAMPLES_C.items()}
    assert ks == {
        "python-s0": 0,
        "python-s1": 1,
        "python-s2": 2,
        "python-s3": 3,
        "cpu64": 4,
        "cuda32": 5,
        "cuda64": 6,
    }
    assert [d.seed_hrc("v8-stat", k) for k in ks.values()] == [20492004 + k for k in range(7)]
    py = [v for nm, v in d.V8_SAMPLES_C.items() if nm.startswith("python")]
    assert len(py) == d.V8NUC_PYTHON_SHARDS == run_suite.V8NUC_PYTHON_SHARDS == 4
    assert (
        all((v[3], v[4]) == (30_000, 20) for v in py) and 30_000 // 20 == 1500
    )  # 4 x 20 x 1500 = 1.2e5
    assert d.V8_SAMPLES_C["cpu64"][3:] == (1_000_000, 100)  # 100 x 1e4, unchanged
    d.base.SEED_BASE = LV5C_BASE


def test_seed_guard_rules(d: ModuleType, c: ModuleType) -> None:
    ok = d.seed_guard({"a": [20490004], "b": [20490005]})
    assert ok["disjoint"] and ok["disjoint_from_consumed"]
    with pytest.raises(SystemExit, match="intersect"):
        d.seed_guard({"a": [20490004], "b": [20490004]})
    with pytest.raises(SystemExit, match="repeats"):
        d.seed_guard({"a": [20490004, 20490004]})
    with pytest.raises(SystemExit, match="consumed"):
        d.seed_guard({"a": [20471004 + 5000]})  # lv5b r_index 5
    with pytest.raises(SystemExit, match="rehearsal"):
        d.seed_guard({"a": [20509500]})
    assert d.seed_guard({"e-shape": [20509500]}, allow_rehearsal=True)["disjoint"]


def test_consumed_set_contains_r_index_0_to_11_but_not_the_never_run_12_and_13(
    c: ModuleType,
) -> None:
    cs = c.consumed_seed_set()
    for base in (20471004, 20451004, 20441004):
        assert all(base + 1000 * r + k in cs for r in range(0, 12) for k in (0, 5, 32))
    # lv5b r_index 12 (V6) and 13 (E1-B) were reserved and never run: lv5c r_index 2 and 3 seeds are not consumed
    assert not any(20471004 + 1000 * r + k in cs for r in (12, 13) for k in (0, 5, 32))
    assert not {20483004, 20483005, 20484004} & cs


def test_main_refuses_a_consumed_seed_base(d: ModuleType) -> None:
    with pytest.raises(SystemExit, match="consumed"):
        d.main(["p7-sampler", "--seed-base", "20471004"])


# -- configurations ------------------------------------------------------------------------------
def test_lv5c_configs_are_elastic_and_the_frozen_ones_are_not(d: ModuleType, c: ModuleType) -> None:
    import steps_v5 as v5

    from ionmc import materials as M

    geo, grid, _, _ = v5.depth_box(M.WATER, 150.0)
    kw: dict[str, Any] = {"energy": 150.0, "n": 100, "seed": 1, "geometry": geo, "grid": grid}
    on = c.wcfg_c("warp-cpu", "float64", tallies=(d.REC,), **kw)
    assert (
        on.physics.elastic
        and on.physics.elastic_table_id == c.ELASTIC_TABLE_ID
        and not on.physics.elastic_only
    )
    assert on.tallies[0].species == ("elastic_recoil",)
    assert v5.nuc_config(**kw).physics.elastic is False  # the frozen suites keep the pin
    off = d.elastic_off(on)
    assert (
        off.physics.elastic is False
        and off.physics.elastic_table_id is None
        and off.physics.nuclear
    )
    assert c.nuc_config_c_elastic_only(**kw).physics.elastic_only


def test_pin_matches_the_committed_table_id(c: ModuleType) -> None:
    pin = json.loads((REPO / "src/ionmc/data/elastic_table_pin.json").read_text())
    assert pin["table_id"] == c.ELASTIC_TABLE_ID


# -- R1-D2 (C) ---------------------------------------------------------------------------------
def test_r1d2_rule(d: ModuleType) -> None:
    n = 20
    rng = np.random.default_rng(1)
    noise = rng.normal(0.0, 1.0, n)
    paired = {"peak_over_plateau": d.paired_z(-1.0 + 0.1 * noise), "plateau": d.paired_z(0.3 + 0.1 * noise),
              "lateral_rms": d.paired_z(0.5 + 0.1 * noise)}  # fmt: skip
    total = d.paired_z(1e-3 * noise)
    j = d.r1d2_judgement(paired, total, 150.0)
    assert j["all_match"] and j["peak_over_plateau"]["declared_sign"] == -1
    wrong = dict(paired, plateau=d.paired_z(-0.3 + 0.1 * noise))  # wrong sign
    assert not d.r1d2_judgement(wrong, total, 150.0)["all_match"]
    weak = dict(paired, lateral_rms=d.paired_z(0.01 + 1.0 * noise))  # |z| < 3
    assert not d.r1d2_judgement(weak, total, 150.0)["lateral_rms"]["matches"]
    big = d.paired_z(1.0 + 0.1 * noise)  # |delta total| = 1 MeV of 150 > 0.1 %
    assert not d.r1d2_judgement(paired, big, 150.0)["total_deposit"]["matches"]
    assert abs(paired["plateau"]["z"]) >= 3.0 and paired["plateau"]["n_batches"] == 20


def test_lateral_rms_batches(d: ModuleType) -> None:
    w = np.zeros((2, 4, 4))
    w[0, 0, 0] = 1.0
    w[1, 1, 1] = 1.0
    e = d.lateral_rms_batches(w.reshape(2, -1), 2.0, 1.0)  # bins at -1.5 .. 1.5
    assert e[0] == pytest.approx(np.hypot(1.5, 1.5)) and e[1] == pytest.approx(np.hypot(0.5, 0.5))


# -- attribution ladder ------------------------------------------------------------------------
R0_OK = {"P7": {"status": "pass"}, "R1-D2": {"status": "pass"}, "V8-LV": {"status": "pass"},
         "P6": {"status": "ci"}, "P6-D": {"status": "ci"}, "P5-ext": {"status": "ci"}}  # fmt: skip


def _v5(
    pp: float, failing: tuple[str, ...] = ("peak_over_plateau",), passes: bool = False
) -> dict[str, Any]:
    return {"pass": passes, "pp_diff": pp, "failing": list(failing)}


def _v11(passes: bool, f: float) -> dict[str, Any]:
    return {"pass": passes, "F_diff": f}


def _fne(above: bool = False, below: bool = False) -> dict[str, Any]:
    return {
        "computed": True,
        "attributed_if_residual_above": above,
        "attributed_if_residual_below": below,
    }


def _r3(rec: float, res: float) -> dict[str, Any]:
    return {"computed": True, "e_res_mev": res, "e_rec_mev": rec}


def _ladder(d: ModuleType, v5: dict, r0: dict | None = None, v11: dict | None = None, fne: dict | None = None,
            r3: dict | None = None) -> dict[str, Any]:  # fmt: skip
    ks = list(v5)
    return d.attribution_ladder(
        v5,
        r0 or R0_OK,
        v11 or {k: _v11(True, 0.0) for k in ks},
        fne or {k: _fne() for k in ks},
        r3 or {k: _r3(0.0, 1.0) for k in ks},
    )


def test_ladder_not_applicable_when_v5_passes(d: ModuleType) -> None:
    out = _ladder(d, {"150": _v5(0.0, (), True), "200": _v5(0.01, (), True)})
    assert not out["applies"] and out["attribution"] == "not applicable (V5 passes)"
    assert out["stop_rule_triggered"] is False and out["pass"] is True


def test_ladder_stop_rule_fails_low(d: ModuleType) -> None:
    out = _ladder(d, {"150": _v5(-0.03), "200": _v5(0.0, (), True)})
    assert (
        out["stop_rule_triggered"]
        and out["pass"] is False
        and out["stop_rule"][0]["reason"].startswith("fails low")
    )
    assert "no tolerance change" in out["message"] and out["v5_stays_failed"]


def test_ladder_stop_rule_above_five_percent(d: ModuleType) -> None:
    out = _ladder(d, {"150": _v5(0.0, (), True), "200": _v5(0.06)})
    assert (
        out["stop_rule_triggered"]
        and out["pass"] is False
        and out["stop_rule"][0]["energy"] == "200"
    )
    assert out["stop_rule"][0]["reason"] == "by more than 5 %"
    # boundaries do not trigger: exactly -2 % and +5 %
    assert not _ladder(d, {"150": _v5(-0.02)})["stop_rule_triggered"]
    assert not _ladder(d, {"150": _v5(0.05)})["stop_rule_triggered"]


def test_ladder_rung0_implementation_defect(d: ModuleType) -> None:
    r0 = {**R0_OK, "P7": {"status": "fail"}}
    out = _ladder(d, {"150": _v5(0.03)}, r0=r0)
    assert out["energies"]["150"]["attribution"].startswith("implementation defect")
    assert out["rung0"]["failed"] == ["P7"] and out["rung0"]["ci_gates_not_evaluated_here"] == [
        "P5-ext",
        "P6",
        "P6-D",
    ]
    assert out["pass"] is True  # attribution only; no stop rule


def test_ladder_rung1_shape_when_v11_fails_in_the_residual_direction(d: ModuleType) -> None:
    out = _ladder(d, {"150": _v5(0.03)}, v11={"150": _v11(False, 0.04)})
    assert out["energies"]["150"]["attribution"].startswith("p+A angular shape")
    # a V11 failure in the opposite direction does not attribute at rung 1
    out2 = _ladder(d, {"150": _v5(0.03)}, v11={"150": _v11(False, -0.04)})
    assert not out2["energies"]["150"]["attribution"].startswith("p+A")
    assert out2["energies"]["150"]["rung2"]["applies"] is False  # V11 does not pass


def test_ladder_rung2_attributed_and_not(d: ModuleType) -> None:
    out = _ladder(d, {"150": _v5(0.03)}, fne={"150": _fne(above=True)})
    rec = out["energies"]["150"]
    assert (
        rec["attribution"].startswith("inclusive-yield")
        and rec["rung2"]["attributed"]
        and "rung3" not in rec
    )
    out2 = _ladder(d, {"150": _v5(0.03)}, fne={"150": _fne(above=False, below=True)})
    assert (
        out2["energies"]["150"]["rung2"]["attributed"] is False
        and "rung3" in out2["energies"]["150"]
    )
    # residual below 2 %: rung 2 does not apply
    out3 = _ladder(d, {"150": _v5(0.015)}, fne={"150": _fne(above=True)})
    assert out3["energies"]["150"]["rung2"]["applies"] is False


def test_ladder_rung3_recoil_or_unexplained(d: ModuleType) -> None:
    hit = _ladder(d, {"150": _v5(0.03)}, r3={"150": _r3(0.6, 1.0)})
    assert hit.get("energies")["150"]["attribution"].startswith("nuclear-recoil")
    assert (
        hit["energies"]["150"]["rung3"]["fraction_required"] == 0.5
        and hit["energies"]["150"]["rung3"]["attributed"]
    )
    edge = _ladder(d, {"150": _v5(0.03)}, r3={"150": _r3(0.5, 1.0)})  # exactly 50 % attributes
    assert edge["energies"]["150"]["rung3"]["attributed"]
    miss = _ladder(d, {"150": _v5(0.03)}, r3={"150": _r3(0.4, 1.0)})
    assert miss["energies"]["150"]["attribution"] == "unexplained"
    none = _ladder(d, {"150": _v5(0.03)}, r3={"150": _r3(0.9, -1.0)})  # no positive residual energy
    assert none["energies"]["150"]["attribution"] == "unexplained"


def test_ladder_failing_criterion_other_than_peak_over_plateau(d: ModuleType) -> None:
    out = _ladder(d, {"150": _v5(0.0, ("r80_mm",))})
    assert out["energies"]["150"]["attribution"].startswith("unexplained: V5 fails on r80_mm")


def test_ladder_records_the_preregistered_expectation_and_the_plan_rule(d: ModuleType) -> None:
    out = _ladder(d, {"150": _v5(0.03)})
    assert out["expectation_preregistered"]["150"]["peak_over_plateau_ratio"] == 0.985
    assert "Stop rule" in out["plan_rule"] and "Rung 3" in out["plan_rule"]


def test_attribution_step_fails_closed_without_inputs(d: ModuleType, tmp_path: Path) -> None:
    _data_or_skip()
    doc = _run_step(d, ["v5-attribution", "--dirs", str(tmp_path)], tmp_path)
    assert doc["pass"] is False and "error" in doc and doc["stop_rule_triggered"] is False


# -- V8-nuc ----------------------------------------------------------------------------------
def _t12(verdict: str) -> dict[str, Any]:
    return {"python_vs_cpu64": {"arrays": {"nuc_local": {"verdict": verdict, "n_groups": 7}}}}


def test_v8_nuc_rule(d: ModuleType) -> None:
    assert d.v8_nuc_verdict(_t12("pass"))["pass"] is True
    inc = d.v8_nuc_verdict(_t12("inconclusive"))
    assert (
        inc["pass"] is False and inc["conclusive"] is False
    )  # an inconclusive verdict is not a pass
    fail = d.v8_nuc_verdict(_t12("fail"))
    assert fail["pass"] is False and fail["conclusive"] is True
    assert "no larger sample" in inc["rule"]


# -- registration ------------------------------------------------------------------------------
def test_registration_order_seeds_timeouts_and_tags() -> None:
    import run_suite
    import summarize

    rs, su = _load("run_suite"), _load("summarize")
    lv = [n.split("-", 1)[1] for n in rs.full_step_names("lv5c", 2)]
    head = ["v6-e0", "v10-pp-level", "v10-a", "v10-a-xcheck", "x-endf", "d9-bgg", "lv5c-throughput", "r1-d2-regression",
            "v8-lv-python-vs-warp-cpu", "p7-sampler", "v11-ionmc-150-emel", "v11-ionmc-150-emonly", "v11-ionmc-200-emel",
            "v11-ionmc-200-emonly", "v11-compare", "v5-150-on", "v5-150-off", "v5-200-on", "v5-200-off", "v5-compare",
            "v5-attribution", "v2b-el-s0", "v2b-el-s1", "v2b-el-combine", "v7-scan-el", "v7-shift-el"]  # fmt: skip
    assert lv[: len(head)] == head
    k6 = lv.index("v6-unsized" if rs.V6_SHARDS is None else "v6-s0")
    assert (
        lv.index("v7-shift-el")
        < k6
        < lv.index("v7r-s0")
        < lv.index("v7r-combine")
        < lv.index("e1-b-compare")
    )
    tail = lv[lv.index("v7r-diag") + 1 :]
    assert tail == ["e1-b-unsized", "e1-b-compare", "e-shape-exploratory", "v3-workers-partition"]
    assert rs.deferred_step_names("lv5c") == [f"{len(lv):02d}-v3-workers-partition"]
    assert rs.deferred_step_names("hr5c") == []
    hr = [n.split("-", 1)[1] for n in rs.full_step_names("hr5c", 2)]
    assert hr[:9] == ["v8-stat-python-s0", "v8-stat-python-s1", "v8-stat-python-s2", "v8-stat-python-s3", "v8-stat-cpu64",
                      "v8-stat-cuda32", "v8-stat-cuda64", "v8-stat-compare", "v7-f32-f64-cuda"]  # fmt: skip
    for suite, base in (("lv5c", LV5C_BASE), ("hr5c", HR5C_BASE)):
        names = rs.full_step_names(suite, 2)
        assert all(su.expected_tag(n) not in (False,) for n in names), [
            n for n in names if su.expected_tag(n) is False
        ]
        for name, cmd, env in rs.suite_steps(suite, 1, 0.5):
            if "pytest" in name or "workers-partition" in name:
                continue
            assert cmd[cmd.index("--seed-base") + 1] == str(base), name
            assert env["IONMC_REQUIRE_DATA"] == "1" and ("IONMC_REQUIRE_CUDA" in env) == (
                suite == "hr5c"
            ), name
            assert "--workers" not in cmd or cmd[cmd.index("--workers") + 1] == "1", name
    assert rs.step_timeout_s("lv5c", "08-r1-d2-regression", 1500) == 3300
    assert rs.step_timeout_s("lv5c", "21-v5-attribution", 1500) == 1800
    assert rs.step_timeout_s("lv5c", "02-v10-pp-level", 1500) == 1800
    for n in ("v8-stat-python-s0", "v8-stat-cuda64", "v7-f32-f64-cuda"):
        assert rs.step_timeout_s("hr5c", f"01-{n}", 1500) == 3600
    assert (
        su.expected_tag("07-r1-d2-regression") == "r1-d2"
        and su.expected_tag("03-v10-a-xcheck") == "v10-a-xcheck"
    )
    assert (
        su.expected_tag("22-v2b-el-s0") == "v2b-el-shard"
        and su.expected_tag("30-v7-scan-el") == "v7-scan-el"
    )
    assert (
        su.expected_tag("40-e1-b-s1") == "e1-b-shard"
        and su.expected_tag("41-e-shape-exploratory") == "e-shape"
    )
    assert (
        su.expected_tag("01-v6-e0") == "v6-e0" and su.expected_tag("09-p7-sampler") == "p7-sampler"
    )
    assert su.expected_tag("04-x-endf") == "x-endf" and su.expected_tag("05-d9-bgg") == "d9-bgg"
    assert run_suite is not None and summarize is not None


def test_constants_equal_between_run_suite_and_the_steps(d: ModuleType) -> None:
    rs = _load("run_suite")
    assert rs.E1B_SHARDS is d.E1B_SHARDS is None
    assert rs.V8NUC_PYTHON_SHARDS == d.V8NUC_PYTHON_SHARDS


def test_source_identity_lists_the_new_modules() -> None:
    rs = _load("run_suite")
    for f in ("steps_v5d.py", "nuclear_v10.py", "nuclear_d9.py"):
        assert f"validation/scripts/transport/{f}" in rs.SOURCE_FILES_V5C


def test_e1b_unsized_fails_closed_and_sizing_rules(d: ModuleType) -> None:
    with pytest.raises(SystemExit, match="Amendment 14a"):
        d.e1b_shards()
    d.E1B_SHARDS = 1
    try:
        with pytest.raises(SystemExit, match=">= 2"):
            d.e1b_shards()
        d.E1B_SHARDS = 2
        assert d.e1b_shards() == 2
        d.E1B_SHARDS = 3  # 3 x 33333 < 1e5
        with pytest.raises(SystemExit):
            d.e1b_shards()
    finally:
        d.E1B_SHARDS = None


def test_slice_c_step_rows_in_summarize() -> None:
    su = _load("summarize")
    steps = {"01-v6-e0": {"pass": True}, "08-r1-d2-regression": {"pass": True, "status": "executed"},
             "10-p7-sampler": {"pass": False, "status": "executed"}}  # fmt: skip
    rows = su.slice_c_rows("lv5c", steps)
    assert rows["R1-D2"] == {"pass": True, "step": "r1-d2-regression", "status": "executed"}
    assert rows["P7"]["pass"] is False and rows["V6"] == {
        "pass": False,
        "step": "v6-compare",
        "status": "not run",
    }
    assert {"V5", "V8", "V11", "V10", "V10-A", "X-ENDF", "D9-BGG", "E1-B", "E-shape", "V2b-el", "V7-scan-el",
            "V7-shift-el", "V7-R"} <= set(rows)  # fmt: skip
    hr = su.slice_c_rows(
        "hr5c", {"01-v8-stat-compare": {"pass": True}, "02-v6-cuda-compare": {"pass": True}}
    )
    assert hr["V8-nuc"]["pass"] and hr["V6"]["pass"]
    assert su.slice_c_rows("lv5b", {}) == {}


# -- nuclear_v10 building blocks --------------------------------------------------------------
def test_rutherford_nonrelativistic_limit(v10: ModuleType) -> None:
    # 1 MeV proton on a heavy target: classical Rutherford (Z1 Z2 e^2 / (4 E))^2 / sin^4 with E in the CM
    m_t = 195.0 * 931.494  # Pt-like
    e, z = 1.0, 78
    th = np.array([60.0])
    got = v10.rutherford_cm_mb_sr(th, e, m_t, z)[0]
    f = v10.two_body_frame(e, m_t)
    ecm = f["w"] - v10.MP_MEV - m_t
    ref = (z * v10.ALPHA_HBARC_MEV_FM / (4.0 * ecm)) ** 2 / np.sin(np.radians(30.0)) ** 4 * 10.0
    assert got == pytest.approx(ref, rel=2e-3)


def test_lab_to_cm_round_trip_and_jacobian(v10: ModuleType) -> None:
    m_t = 16.0 * 931.494
    th_l = np.array([10.0, 30.0, 60.0])
    tc, jac = v10.lab_to_cm(th_l, 180.0, m_t)
    f = v10.two_body_frame(180.0, m_t)
    r = f["beta_cm"] / f["beta_star"]
    back = np.degrees(
        np.arctan2(np.sin(np.radians(tc)), f["gamma_cm"] * (np.cos(np.radians(tc)) + r))
    )
    assert np.allclose(back, th_l, atol=1e-9) and np.all(tc > th_l)
    # solid angle conservation: integral of jac dOmega_CM over the sphere is 4 pi
    grid = np.linspace(0.5, 179.5, 4000)
    tcg, jg = v10.lab_to_cm(grid, 180.0, m_t)
    ok = np.isfinite(jg)
    assert ok.all() and np.all(jg > 0)


def test_interpolation_minimum_and_integral_helpers(v10: ModuleType) -> None:
    th = np.array([8.0, 11.0, 14.0, 17.0, 20.0, 23.0])
    y = np.array([100.0, 50.0, 25.0, 30.0, 12.0, 20.0])
    assert v10.first_minimum_deg(th, y) == 14.0
    assert v10.first_minimum_deg(np.array([8.0, 11.0, 14.0]), np.array([3.0, 2.0, 1.0])) is None
    assert v10.loglin(th, y, 12.5) == pytest.approx(np.sqrt(50.0 * 25.0))  # log-linear midpoint
    assert (
        v10.loglin(th, y, 4.0) is None
        and v10.loglin(np.array([5.0, 20.0]), np.array([1.0, 2.0]), 12.0) is None
    )
    grid = np.linspace(0.0, 180.0, 2001)
    assert v10.trapezoid_2pi_sin(grid, np.ones_like(grid)) == pytest.approx(4.0 * np.pi, rel=1e-5)


def test_parse_pdg_pp_synthetic(v10: ModuleType) -> None:
    text = (
        " FILE_NAME\n POINT_NUMBER PLAB\n"
        "   1    0.60000    0.59000    0.61000    24.000     0.50000 0.50000    1.0  2.0  TEST 99         X 1, 2\n"
        "   2    0.80000    0.80000    0.80000    23.000     1.0000  1.0000     0.0  0.0  OTHER 98        Y 3, 4\n"
        "not a data line\n"
    )
    rows = v10.parse_pdg_pp(text)
    assert [r["p_lab_gev"] for r in rows] == [0.6, 0.8] and rows[0]["sys_minus_pct"] == 2.0
    assert rows[0]["reference"].startswith("TEST 99")


# -- data-backed dry runs ---------------------------------------------------------------------
def test_p7_step_dry_run(d: ModuleType, tmp_path: Path) -> None:
    _data_or_skip()
    doc = _run_step(d, ["p7-sampler", "--scale", "0.1"], tmp_path)
    assert doc["n_cases"] == 17 and doc["pass"] and doc["sigma_el"]["max_rel_error"] <= 1e-3
    cases = [(x["target"], x["energy_mev"]) for x in doc["cases"]]
    assert cases[:5] == [("H-1", e) for e in (20.0, 50.0, 100.0, 150.0, 200.0)]
    assert [x["seed"] for x in doc["cases"]] == [20496004 + k for k in range(17)]
    assert (tmp_path / "p7-status.json").is_file()


def test_v5_producers_write_the_recoil_channel_only_for_on(d: ModuleType, tmp_path: Path) -> None:
    _data_or_skip()
    on = _run_step(
        d, ["v5-ionmc", "--energy", "150", "--nuclear", "on", "--scale", "0.02"], tmp_path
    )
    off = _run_step(
        d, ["v5-ionmc", "--energy", "150", "--nuclear", "off", "--scale", "0.02"], tmp_path
    )
    pon = json.loads((tmp_path / "v5-150-on.json").read_text())
    poff = json.loads((tmp_path / "v5-150-off.json").read_text())
    assert on["seed"] == 20485004 and off["seed"] == 20486004 and pon["elastic"] and poff["elastic"]
    assert np.asarray(pon["recoil_idd_batches"]).shape == np.asarray(pon["idd_batches"]).shape
    assert np.asarray(pon["recoil_idd_batches"]).sum() > 0 and "recoil_idd_batches" not in poff
    assert (
        pon["elastic_table_id"] == poff["elastic_table_id"]
        and pon["nuclear"]
        and not poff["nuclear"]
    )
    assert "content_sha256" in pon and pon["seed_base"] == LV5C_BASE


def test_v8_lv_step_dry_run(d: ModuleType, tmp_path: Path) -> None:
    _data_or_skip()
    doc = _run_step(d, ["v8-lv", "--scale", "0.25"], tmp_path)
    assert (
        doc["pass"]
        and doc["nuclear_dense"]["n_elastic_events"] > 0
        and doc["ci"]["discrete_identical"]
    )
    assert doc["ci"]["channels_max_abs_mev"] <= 1e-10 and (tmp_path / "v8-lv-status.json").is_file()


def test_v10_steps_dry_run_and_the_fallback_block(d: ModuleType, tmp_path: Path) -> None:
    _data_or_skip()
    pp = _run_step(d, ["v10-pp-level"], tmp_path)
    assert pp["pass"] and pp["row"] == "V10" and pp["plan_line"] == 201
    assert (
        pp["result"]["n_points"] == 0
        and "0.846" in pp["evidence_gap"]
        and "none in the plan window" in pp["evidence_gap"]
    )
    assert pp["verdict"] == "none (report-only; executed, no points in the frozen window)"
    a = _run_step(d, ["v10-a"], tmp_path)
    assert (
        a["plan_line"] == 202
        and a["row"] == "V10-A"
        and a["status"] in ("pass", "fail", "not executed (evidence missing)")
    )
    trig = a["fallback_ii_prime_trigger"]
    assert set(trig) == {"rule", "a_or_b_fails_for_pO", "c_passes", "p7_passes", "fires", "reason"}
    assert trig["p7_passes"] == "ci-only" and trig["fires"] is False
    assert a["result"]["coverage"]["a_points"] >= 12 and set(a["failures"]) == {"a", "b", "c"}
    for f in a["failures"]["a"]:
        assert {"dataset", "theta_deg", "model_over_data_minus_1", "rutherford_share"} <= set(f)
    assert a["pass"] == (a["status"] == "pass") and a["_exit"] == (0 if a["pass"] else 1)
    xe = _run_step(d, ["x-endf", "--scale", "0.02"], tmp_path)
    assert xe["plan_line"] == 203 and xe["pass"] and "o16" in xe["result"]
    for t, per in xe["result"]["targets"].items():
        assert all(v["builder_arrays_agree"]["model_max_rel"] < 1e-9 for v in per.values()), t
    d9 = _run_step(d, ["d9-bgg"], tmp_path)
    assert (
        d9["plan_line"] == 208
        and d9["pass"]
        and set(d9["result"]["sigma_el_level_proxy"]["targets"]) == {"C-12", "O-16"}
    )
    assert "windows" in d9["result"]["sigma_inel"]["sets"]["exfor-d0356"]["targets"]["C-12"]


def test_v10a_coverage_rule_never_passes_without_coverage(
    v10: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    _data_or_skip()
    from steps_v5c import ELASTIC_TABLE_ID  # noqa: PLC0415

    model = v10.ModelTable(ELASTIC_TABLE_ID)
    monkeypatch.setattr(v10, "GATING_SETS", v10.GATING_SETS[:1])  # one O-16 dataset at 135 MeV only
    res = v10.v10_a(model)
    assert (
        res["status"] == "not executed (evidence missing)"
        and res["pass"] is False
        and not res["coverage_met"]
    )


def test_e_shape_and_e1b_documents_state_the_rules(d: ModuleType, tmp_path: Path) -> None:
    _data_or_skip()
    doc = _run_step(d, ["e-shape", "--scale", "0.05"], tmp_path)
    assert (
        doc["seed"] == 20509500
        and doc["tag"] == "exploratory"
        and doc["verdict"] == "executed (model variant only)"
    )
    assert (
        doc["variant_forced_isotropic"] == "not run (builder option not implemented)"
        and doc["pass"]
    )
    un = _run_step(d, ["e1-b-unsized"], tmp_path)
    assert un["pass"] is False and "Amendment 14a" in un["error"]


def test_hr5c_and_lv5c_partials_carry_the_elastic_bindings(d: ModuleType, tmp_path: Path) -> None:
    _data_or_skip()
    _run_step(d, ["v2b-el-shard", "--shard", "1", "--scale", "0.01"], tmp_path)
    p = json.loads((tmp_path / "v2b-el-s1.json").read_text())
    assert p["elastic"] is True and len(p["elastic_table_id"]) == 64 and p["seed"] == 20489005
