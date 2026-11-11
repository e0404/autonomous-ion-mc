"""V3-008 suites ``lv8`` and ``hr8`` of ``run_suite.py`` and their steps in ``steps_v8.py``: step
lists and names, seed derivation (plan section C), the hashed source set, the summarizer tags and
a tiny-scale smoke of the cheap CPU steps (CI seeds: the steps' own derivation from the lv8 base at
``--scale`` far below 1, labelled ``reduced``; nothing here is qualification evidence)."""

from __future__ import annotations

import importlib.util
import io
import json
import sys
from contextlib import redirect_stdout
from pathlib import Path
from types import ModuleType

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "validation" / "scripts" / "transport"
RESERVED = range(20572004, 20632004)


def _load(name: str) -> ModuleType:
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


RS = _load("run_suite")
V8 = _load("steps_v8")
SUM = _load("summarize")

LV8 = [
    "01-i6-lv-deterministic-end-depth",
    "02-i8-range-straggling-he150",
    "03-i8-range-straggling-c290",
    "04-i8-range-straggling-o350",
    *[f"{i + 5:02d}-i13-r80-vs-csda-{n}" for i, n in enumerate(
        ("he50", "he150", "he250", "c80", "c290", "c430", "o100", "o350", "o430"))],
    "14-t1-ion-trace-parity-c290-o350",
    "15-r1-ion-repeatability-c290",
    "16-r2-a16-trace-fixture",
]  # fmt: skip


def test_suite_step_names_and_seed_bases() -> None:
    assert {"lv8", "hr8"} <= set(RS.SUITES)
    assert RS.full_step_names("lv8", 2) == LV8
    assert RS.full_step_names("hr8", 2) == ["01-t13-ion-chunks-cuda-c290"]
    assert RS.DEFAULT_SEED_BASES["lv8"] == 20581004 == V8.LV8_BASE
    assert RS.DEFAULT_SEED_BASES["hr8"] == 20591004 == V8.HR8_BASE
    assert RS.DEFERRED_STEPS["lv8"] == () and RS.DEFERRED_STEPS["hr8"] == ()
    assert RS.step_timeout_s("lv8", "16-r2-a16-trace-fixture", 1500) == 3300
    for name in LV8:  # every step maps to the document tag its step prints
        assert SUM.expected_tag(name) not in (None, False), name
    assert SUM.expected_tag("01-t13-ion-chunks-cuda-c290") == "t13-ion"
    assert SUM.expected_tag("14-t1-ion-trace-parity-c290-o350") == "t1-ion"
    assert SUM.expected_tag("15-r1-ion-repeatability-c290") == "r1-ion"
    assert SUM.seed_blockers(20581004, "lv8") == [] and SUM.seed_blockers(20591004, "hr8") == []
    assert SUM.seed_blockers(20601004, "lv8")


def test_argv_carries_seed_base_timeout_and_cuda_env() -> None:
    steps = RS.suite_steps("lv8", 1, 0.01, seed_base=20601004)
    for _name, cmd, _env in steps:
        assert cmd[1].endswith("steps_v8.py")
        assert cmd[cmd.index("--seed-base") + 1] == "20601004"
        assert "--timeout" in cmd
    by = {n: c for n, c, _ in steps}
    assert by["03-i8-range-straggling-c290"][2:5] == ["i8", "--point", "1"]
    assert "--scale" not in by["16-r2-a16-trace-fixture"]
    (name, cmd, env), = RS.suite_steps("hr8", 1, 1.0)  # fmt: skip
    assert env == {"IONMC_REQUIRE_CUDA": "1"} and cmd[2] == "t13-ion"
    assert cmd[cmd.index("--seed-base") + 1] == "20591004"


def test_hashed_source_set() -> None:
    v8 = RS.source_file_list("lv8")
    assert v8 == RS.source_file_list("hr8")
    assert v8[: len(RS.SOURCE_FILES_V5)] == RS.SOURCE_FILES_V5
    for f in ("validation/plans/v3-008-acceptance.md", "decisions/0043-ion-em-transport.md",
              "validation/scripts/transport/steps_v8.py"):  # fmt: skip
        assert f in v8 and (ROOT / f).is_file()


def test_row_r_index_and_seeds() -> None:
    idx = {**V8.ROW_R_INDEX, **V8.RESERVED_R_INDEX}
    assert len(set(idx.values())) == len(idx)
    assert all(42 <= r <= 60 for r in idx.values())  # 30-41 are the td-tl-sample rows
    assert not set(idx.values()) & {V8.r_index_of(c) for c in V8.CASES}
    assert V8.row_seed("i8", 1) == 20581004 + 4200 + 1
    assert V8.row_seed("t13-ion", 1, V8.HR8_BASE) == 20591004 + 5500 + 1
    seeds = {V8.row_seed(r, k, b) for r in V8.ROW_R_INDEX for k in range(100)
             for b in (V8.LV8_BASE, V8.HR8_BASE)}  # fmt: skip
    assert all(s in RESERVED for s in seeds)
    td = {V8.seed_of(c, k) for c in V8.CASES for k in range(100)}
    assert not seeds & td
    with pytest.raises(ValueError):
        V8.row_seed("i8", 100)


def _run(argv: list[str]) -> tuple[int, dict]:  # type: ignore[type-arg]
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = V8.main(argv)
    text = buf.getvalue()
    return rc, json.loads(text.split("#JSON-BEGIN", 1)[1].split("#JSON-END", 1)[0])


def test_i8_quadrature_is_finite_and_scales_like_the_proton() -> None:
    """sigma_R/R(C 290) with corrections on is close to sigma_R/R(proton at equal velocity) times
    sqrt(m_p/M) (the I4b identity holds within 1e-3 without corrections; a few % with them)."""
    proj = V8.CARBON12
    sq, rq = V8.sigma_range_quadrature(proj, 290.0 * proj.a, 2.0 * proj.a)
    assert 0.001 < sq / rq < 0.006
    r = V8.csda_mm(proj, 290.0)
    assert abs(rq / r - 1.0) < 0.01


def test_smoke_i6_and_i8_tiny_scale() -> None:
    rc, doc = _run(["i6-lv"])
    assert rc == 0 and doc["pass"] and doc["step"] == "i6-lv" and len(doc["cases"]) == 18
    rc, doc = _run(["i8", "--point", "0", "--scale", "0.002"])
    assert doc["reduced"] and doc["histories"] == 200 and doc["seed"] == 20581004 + 4200
    assert doc["counters"] and doc["valid"] and doc["end_codes_all_cutoff"]
    assert "i8b_report_only" in doc and rc in (0, 1)  # 200 histories: the 5 % bound may fail
