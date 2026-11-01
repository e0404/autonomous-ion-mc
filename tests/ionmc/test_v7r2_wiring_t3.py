"""Wiring task T3 of the lv5d/hr5d suites (plan Amendment 22 rule 2, Amendment 20 (f)6 / (f)8; Codex REVIEW-fe47f290):
the full lv5d fixture preflight (absent / inconsistent / stale / generation-exempt / valid), the measured calibration
rate gate (absent / valid / malformed record, re-split), ``--print-hashed-source-digest`` for lv5d / hr5d, the
``v7r2-s{k}`` generation path and the builder's ``--layout v7r2``. Unit tests only: no host run, no seed drawn."""

# ruff: noqa: E501

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"
DIGEST = "d" * 64
RATE_ENV = "IONMC_V7R2_CAL_RATE_RECORD"


def _load(name: str) -> ModuleType:
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def rs() -> ModuleType:
    return _load("run_suite")


@pytest.fixture(scope="module")
def s5c() -> ModuleType:
    return _load("steps_v5c")


@pytest.fixture(scope="module")
def s5d() -> ModuleType:
    return _load("steps_v5d")


@pytest.fixture(scope="module")
def cal() -> ModuleType:
    return _load("v7r2_cal_rate")


@pytest.fixture(autouse=True)
def _no_rate_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(RATE_ENV, raising=False)


# -- synthetic v7r2-layout fixture (the helpers of test_v7r2_calibration.test_v7r2_rehearsal_loader_layout) ----
def build_v7r2_fixture(d: Path, s5c: ModuleType, digest: str, *, seed_offset: int = 0, layout_shards: int = 4,
                       code_sha: str = "c0de" * 10, seed_base: int = 20505000) -> dict[str, Any]:  # fmt: skip
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(4)
    prov: dict[str, Any] = {"host_run_ids": ["RUN-X"], "hashed_source_digest": digest, "seed_base": seed_base,
                            "code_sha": code_sha}  # fmt: skip
    digests = {}
    for k in range(layout_shards):
        sums = rng.gamma(5.0, 100.0, (225, 20))
        side = s5c.write_v7r_sidecar(d / f"v7r2-s{k}-escaped_neutral.npz", sums)
        per = sums / 500.0
        doc = {"row": "v7-r2", "kind": "shard", "shard": k, "seed": seed_base + 16000 + k + seed_offset,
               "seed_base": seed_base, "n": 2_250_000, "replicates": 225, "batch_histories": 500,
               "blocks_per_replicate": 20, "valid": True, "reduced": False, "git_sha": code_sha, "sidecar": side,
               "estimators": {"escaped_neutral": {"mean": per.mean(axis=1).tolist(),
                                                  "sem": (per.std(axis=1, ddof=1) / 20**0.5).tolist()}}}  # fmt: skip
        doc["content_sha256"] = s5c.v5.content_digest(doc)
        (d / f"v7r2-s{k}.json").write_text(json.dumps(doc))
        digests[f"v7r2-s{k}.json"] = doc["content_sha256"]
    prov["partial_digests"] = digests
    (d / "PROVENANCE.json").write_text(json.dumps(prov))
    return prov


SHA = "a" * 40


@pytest.fixture(scope="module")
def cov() -> ModuleType:
    """The real seal function of the calibration records (``test_v7r_coverage._seal``, used by the host parts)."""
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(
        "test_v7r_coverage_module_t3", REPO / "tests" / "ionmc" / "test_v7r_coverage.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_artifacts(d: Path, cov: ModuleType, *, wall: float = 996.6, part: int = 0, surrogate: str = "archive_fitted",
                   sha: str = SHA, total: float | None = None, over: dict[str, Any] | None = None,
                   unseal: bool = False, log: str | None = None) -> tuple[Path, Path]:  # fmt: skip
    """A genuine sealed V7-R2 calibration part record (sealed by the real ``_seal``) and a synthetic but well-formed
    passing pytest log of that one part (a controlled host run is not needed to test the verification)."""
    d.mkdir(parents=True, exist_ok=True)
    n = 83 if part < 2 else 82
    rec: dict[str, Any] = {
        "schema": "ionmc-v7r2-calibration-record-1", "surrogate": surrogate, "part": part, "n": n,
        "reps_env": 3200, "layout": {"replicates": 10800, "reps_per_shard": 225, "parts": 39},
        "run_sha": sha, "suite": "lv5d", "wall_mc_s": wall, "wall_per_replicate_s": wall / n,
    }  # fmt: skip
    rec.update(over or {})
    sealed = cov._seal(rec)
    if unseal:
        sealed["wall_mc_s"] = sealed["wall_mc_s"] + 1.0  # tampered after sealing
    path = d / f"v7r2-calibration-{surrogate}-p{part:02d}.json"
    path.write_text(json.dumps(sealed, sort_keys=True))
    total = wall + 70.0 if total is None else total
    lg = d / "part.log"
    lg.write_text(
        log
        if log is not None
        else (
            f"collected 78 items / 77 deselected / 1 selected\nPARTIAL {path.name} {cov._seal(rec)['content_sha256']}\n"
            f".\n====== 1 passed, 77 deselected in {total:.2f}s ({int(total) // 3600}:{int(total) // 60 % 60:02d}:{int(total) % 60:02d}) ======\n"
        )
    )
    return path, lg


def write_rate(
    path: Path, cal: ModuleType, cov: ModuleType, *, s_per_rep: float = 12.0, **over: Any
) -> Path:
    """A measured-rate record built by ``build_record`` from genuine artifacts (``over``: fields replaced afterwards)."""
    part, lg = make_artifacts(path.parent / f"{path.stem}-art", cov, wall=s_per_rep * 83)
    rec = cal.build_record(part, lg, part.parent, SHA)
    rec["host"] = {"node": "slowhost", "cpu_count": 8, "platform": "linux", "python": "3.12"}
    rec.update(over)
    path.write_text(json.dumps(rec))
    return path


@pytest.fixture
def lv5d_tree(rs: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """``rs.REPO`` is an empty temporary tree and the lv5d hashed-source digest is a constant; returns the fixture dir."""
    monkeypatch.setattr(rs, "REPO", tmp_path)
    monkeypatch.setattr(rs, "hashed_source_digest", lambda suite="lv5c": DIGEST)
    return tmp_path / rs.REHEARSAL_FIXTURE_DIR


# -- (1) the full lv5d preflight -------------------------------------------------------------------------------
def test_lv5d_preflight_absent_fixture_is_refused(rs: ModuleType, lv5d_tree: Path) -> None:
    refusal = rs.rehearsal_preflight("lv5d")
    assert refusal and rs.REHEARSAL_PREFLIGHT_STEP in refusal and "absent" in refusal


def test_lv5d_preflight_valid_fixture_needs_the_rate_record(
    rs: ModuleType,
    s5c: ModuleType,
    cal: ModuleType,
    cov: ModuleType,
    lv5d_tree: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    build_v7r2_fixture(lv5d_tree, s5c, DIGEST)
    refusal = rs.rehearsal_preflight("lv5d")
    assert refusal and "calibration rate not measured" in refusal  # fixture ok, rate not measured
    monkeypatch.setenv(RATE_ENV, str(tmp_path / "missing.json"))
    assert "not ok" in (rs.rehearsal_preflight("lv5d") or "")  # a path that is no file
    monkeypatch.setenv(RATE_ENV, str(write_rate(tmp_path / "rate.json", cal, cov)))
    assert (
        rs.rehearsal_preflight("lv5d") is None
    )  # a valid synthetic v7r2-layout fixture and a valid record pass
    # a selection without the throughput step and without a V7-R2 pytest step needs no rate record
    monkeypatch.delenv(RATE_ENV)
    assert rs.rehearsal_preflight("lv5d", ["v6f-s0"]) is None
    assert "calibration rate" in (rs.rehearsal_preflight("lv5d", ["lv5d-throughput"]) or "")
    assert "calibration rate" in (
        rs.rehearsal_preflight("lv5d", ["pytest-v7r2-calibration-s3"]) or ""
    )


@pytest.mark.parametrize(
    ("how", "expect"),
    [
        ("digest", "inconsistent"),  # a corrupted partial / shard digest recorded in PROVENANCE
        ("seed", "seed is not seed_base"),
        ("layout", "need exactly 4 shards"),  # the V7-R layout (v7r-s*) is no v7r2 fixture
        ("three", "need exactly 4 shards"),
        ("tamper", "content_sha256 does not reproduce"),
    ],
)
def test_lv5d_preflight_inconsistent_fixture_is_refused(
    how: str, expect: str, rs: ModuleType, s5c: ModuleType, lv5d_tree: Path
) -> None:
    prov = build_v7r2_fixture(lv5d_tree, s5c, DIGEST, seed_offset=7 if how == "seed" else 0,
                              layout_shards=3 if how == "three" else 4)  # fmt: skip
    if how == "digest":
        prov["partial_digests"]["v7r2-s2.json"] = "0" * 64
        (lv5d_tree / "PROVENANCE.json").write_text(json.dumps(prov))
    elif how == "layout":
        for p in lv5d_tree.glob("v7r2-s*"):
            p.rename(p.with_name(p.name.replace("v7r2-s", "v7r-s")))
    elif how == "tamper":
        p = lv5d_tree / "v7r2-s1.json"
        doc = json.loads(p.read_text())
        doc["wall_s"] = 1.0
        p.write_text(json.dumps(doc))
    refusal = rs.rehearsal_preflight("lv5d")
    assert refusal and rs.REHEARSAL_PREFLIGHT_STEP in refusal
    assert expect in refusal or how == "digest" and "differs from the digest recorded" in refusal, (
        refusal
    )


def test_lv5d_preflight_stale_fixture_is_refused(
    rs: ModuleType, s5c: ModuleType, lv5d_tree: Path
) -> None:
    build_v7r2_fixture(lv5d_tree, s5c, "e" * 64)  # recorded digest != the current lv5d digest
    refusal = rs.rehearsal_preflight("lv5d")
    assert refusal and "stale" in refusal and DIGEST in refusal


def test_lv5d_preflight_uses_the_lv5d_digest(rs: ModuleType, s5c: ModuleType, lv5d_tree: Path,
                                             monkeypatch: pytest.MonkeyPatch) -> None:  # fmt: skip
    seen: list[str] = []
    monkeypatch.setattr(
        rs, "hashed_source_digest", lambda suite="lv5c": seen.append(suite) or DIGEST
    )
    build_v7r2_fixture(lv5d_tree, s5c, DIGEST)
    rs.rehearsal_preflight("lv5d", ["v6f-s0"])
    assert seen == ["lv5d"]


@pytest.mark.parametrize(
    ("only", "seed_base", "exempt"),
    [
        (["v7r2-s0"], 20505000, True),
        (["41-v7r2-s3"], 20509499, True),  # registered names as main() passes them
        (["41-v7r2-s0", "42-v7r2-s1"], 20505000, True),
        (["v7r2-s4"], 20505000, False),  # the rehearsal has the shards 0..3 only
        (["v7r2-s47"], 20505000, False),
        (["v7r2-s0", "v7r2-ref"], 20505000, False),
        (["v7r2-combine"], 20505000, False),
        (["pytest-v7r2-calibration-s0"], 20505000, False),
        (["v7r-s0"], 20505000, False),  # the V7-R shard name is not an lv5d generation step
        (["v7r2-s0"], 20531004, False),  # the lv5d evidence base
        (["v7r2-s0"], 20509500, False),
        (["v7r2-s0"], None, False),
        (None, 20505000, False),
        (["v7r2-s"], 20505000, False),
    ],
)
def test_lv5d_generation_run_is_exempt_only_for_v7r2_shards_in_the_rehearsal_family(
    only: list[str] | None, seed_base: int | None, exempt: bool, rs: ModuleType, lv5d_tree: Path
) -> None:
    assert rs.rehearsal_generation_run(only, seed_base, "lv5d") is exempt
    refusal = rs.rehearsal_preflight("lv5d", only, seed_base)
    assert (refusal is None) is exempt  # no fixture is present here
    assert rs.rehearsal_generation_run(only, seed_base, "lv5b") is False
    assert (
        rs.rehearsal_generation_run(["v7r2-s0"], 20505000, "lv5c") is False
    )  # the V7-R suite keeps its own names


def test_hr5d_and_other_suites_are_not_preflighted(rs: ModuleType, lv5d_tree: Path) -> None:
    """hr5d registers no V7-R2 pytest step and consumes no rehearsal fixture."""
    steps = rs.suite_steps("hr5d", 1, 1.0)
    assert not any(rs.is_v7r_calibration_step(n) for n, _, _ in steps)
    assert rs.rehearsal_preflight("hr5d") is None and rs.rehearsal_preflight("hr5c") is None
    assert rs.rehearsal_preflight("lv5c") and "absent" in rs.rehearsal_preflight("lv5c")  # type: ignore[operator]


def test_lv5c_preflight_is_unchanged(rs: ModuleType, lv5d_tree: Path) -> None:
    d = lv5d_tree
    d.mkdir(parents=True)
    (d / "PROVENANCE.json").write_text(json.dumps({"hashed_source_digest": DIGEST}))
    assert rs.rehearsal_preflight("lv5c") is None  # presence and digest only (archived behaviour)


def test_calibration_steps_keep_the_required_fixture_environment(rs: ModuleType) -> None:
    steps = dict((n.split("-", 1)[1], e) for n, _, e in rs.suite_steps("lv5d", 1, 1.0))
    for name in (
        "pytest-v7r2-kappa-s0",
        "pytest-v7r2-calibration-s0",
        "pytest-v7r2-calibration-s77",
    ):
        assert steps[name]["IONMC_V7R_FIXTURES"] == "required"


# -- (2) the measured-rate gate of the throughput step ---------------------------------------------------------
def test_rep_plan_table_fails_closed_without_a_measured_rate(s5d: ModuleType) -> None:
    rep = s5d.rep_plan_table("lv5d")
    assert (
        rep["ok"] is False
        and rep["reason"] == "not ok: calibration rate not measured"
        and rep["plans"] == {}
    )
    with pytest.raises(SystemExit, match="calibration rate not measured"):
        s5d.declared_plans("lv5d")
    with pytest.raises(SystemExit, match="belong to lv5d"):
        s5d.rep_plan_table("hr5d")


def test_rep_plan_table_is_recomputed_from_the_measured_rate(
    s5d: ModuleType,
    cal: ModuleType,
    cov: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rec = write_rate(tmp_path / "rate.json", cal, cov, s_per_rep=12.0)
    monkeypatch.setenv(RATE_ENV, str(rec))
    rep = s5d.rep_plan_table("lv5d")
    assert rep["ok"] and rep["rate_source"] == "measured" and rep["reason"] is None
    assert rep["rate_record"]["sha256"] == hashlib.sha256(rec.read_bytes()).hexdigest()
    assert (
        rep["rate_record"]["host"]["node"] == "slowhost"
        and rep["rate_record"]["git_sha"] == "a" * 40
    )
    cal_plan = rep["plans"]["v7r2 calibration (3200 reps x 39 parts per surrogate)"]
    assert (
        cal_plan["s_per_item"] == 12.0 and cal_plan["fixed_s"] == 70.0
    )  # the measured rate, not the 12.2 default
    assert s5d.declared_plans(
        "lv5d", cal_rate=cal.load_record(rec)
    )  # the gated table is returned with a record
    # a slower measured rate exceeds the 1100 s target of 39 parts: not ok, the re-split is stated
    slow = write_rate(tmp_path / "slow.json", cal, cov, s_per_rep=13.0)
    bad = s5d.rep_plan_table("lv5d", str(slow))
    assert bad["ok"] is False and bad["calibration_parts_needed"] > 39
    assert "re-split" in bad["reason"] and str(bad["calibration_parts_needed"]) in bad["reason"]


@pytest.mark.parametrize(
    "case",
    [
        "not json",
        "not an object",
        "schema",
        "negative rate",
        "nan rate",
        "rate vs wall",
        "no host",
        "short sha",
        "no log_dir",
        "bad part digest",
        "fixed below planning value",
        "zero reps",
    ],
)
def test_rep_plan_table_malformed_record_is_not_ok(
    case: str, s5d: ModuleType, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    p = write_rate(tmp_path / "r.json", cal, cov)
    doc = json.loads(p.read_text())
    mutate = {
        "schema": lambda d: d.update(schema="x"), "negative rate": lambda d: d.update(s_per_rep=-1.0),
        "nan rate": lambda d: d.update(s_per_rep=float("nan")), "rate vs wall": lambda d: d.update(wall_mc_s=1.0),
        "no host": lambda d: d.pop("host"), "short sha": lambda d: d.update(git_sha="abc"),
        "no log_dir": lambda d: d.update(log_dir=""), "bad part digest": lambda d: d.update(part_log_sha256="z"),
        "fixed below planning value": lambda d: d.update(fixed_s=10.0), "zero reps": lambda d: d.update(reps=0),
    }  # fmt: skip
    if case == "not json":
        p.write_text("{")
    elif case == "not an object":
        p.write_text("[1]")
    else:
        mutate[case](doc)
        p.write_text(json.dumps(doc))
    rep = s5d.rep_plan_table("lv5d", str(p))
    assert rep["ok"] is False and rep["reason"].startswith("not ok:") and rep["plans"] == {}


def test_throughput_step_gate_is_pass_false_without_the_record(
    s5d: ModuleType, cal: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The step source: the lv5d verdict needs ``rep_ok`` even at reduced scale and records the reason in the
    document (the step runs the real Simulation, so the wiring is checked on the source, the table above in values)."""
    import inspect

    src = inspect.getsource(s5d.step_throughput)
    assert '"pass": bool(ok and rep_ok and (a.scale < 1.0 or fits))' in src
    assert '"cal_rate_status"' in src and "_declared_plans(suite)" in src


def test_cal_rate_builder_derives_the_record_from_a_part_record_and_log(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    part, log = make_artifacts(tmp_path / "art", cov, wall=996.6, total=1100.6)
    out = tmp_path / "rate.json"
    args = ["--part-record", str(part), "--part-log", str(log), "--log-dir", str(part.parent),
            "--git-sha", SHA, "--out", str(out)]  # fmt: skip
    assert cal.main(args) == 0
    got = cal.load_record(out)["record"]
    assert abs(got["s_per_rep"] - 996.6 / 83) < 1e-12 and abs(got["fixed_s"] - 104.0) < 1e-6
    assert got["part_record_sha256"] == hashlib.sha256(part.read_bytes()).hexdigest()
    with pytest.raises(SystemExit, match="exists"):
        cal.main(args)
    log.write_text("2 failed in 3s\n")
    with pytest.raises(ValueError, match="summary"):
        cal.build_record(part, log, part.parent, SHA)
    out.unlink()
    with pytest.raises(SystemExit, match="not ok"):
        cal.main(args)
    assert not out.exists()


# -- Codex REVIEW-6519d4ff: the measured-rate record is only accepted with verified controlled-run artifacts ------------
def _rate_from(tmp_path: Path, cal: ModuleType, cov: ModuleType, **kw: Any) -> Path:
    """Build a record from artifacts made with ``kw`` through the unchecked path (artifacts are valid only if ``kw`` is
    empty): the record is assembled by hand from the files so that ``load_record`` is what refuses."""
    part, lg = make_artifacts(tmp_path / "art", cov, **kw)
    doc = json.loads(part.read_text())
    wall, n = doc["wall_mc_s"], doc["n"]
    total = wall + 70.0
    rec = {"schema": cal.CAL_RATE_SCHEMA, "s_per_rep": wall / n, "fixed_s": 70.0, "reps": n, "wall_mc_s": wall,
           "wall_total_s": total, "host": {"node": "h"}, "git_sha": SHA, "log_dir": str(part.parent),
           "part_record": part.name, "part_record_sha256": hashlib.sha256(part.read_bytes()).hexdigest(),
           "part_log": lg.name, "part_log_sha256": hashlib.sha256(lg.read_bytes()).hexdigest()}  # fmt: skip
    out = tmp_path / "rate.json"
    out.write_text(json.dumps(rec))
    return out


def test_genuine_artifacts_load_and_authorise_the_39_part_plan(
    s5d: ModuleType, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    got = cal.load_record(_rate_from(tmp_path, cal, cov))
    assert cal.require_validated(got) == got
    assert s5d.declared_plans("lv5d", cal_rate=got)["v7r2 reference (1e6, block rate)"]
    assert s5d.rep_plan_table("lv5d", got["path"])["ok"]


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"unseal": True}, "not sealed"),
        ({"sha": "b" * 40}, "run_sha"),
        ({"over": {"suite": "lv5c"}}, "suite"),
        ({"over": {"layout": {"replicates": 14400, "reps_per_shard": 225, "parts": 39}}}, "layout"),
        ({"over": {"schema": "ionmc-v7r-calibration-record-1"}}, "schema"),
        ({"over": {"n": 40}}, "planned size"),
        ({"over": {"surrogate": "rehearsal_block_sums"}}, "disagree"),
        ({"over": {"reps_env": 100}}, "layout / reps_env"),
        ({"over": {"run_sha": "unknown"}}, "run_sha"),
        ({"log": "====== 1 passed in 1100.60s (0:18:20) ======\n"}, "PARTIAL"),
        (
            {
                "log": "PARTIAL v7r2-calibration-archive_fitted-p00.json "
                + "0" * 64
                + "\n====== 1 passed in 1100.60s ======\n"
            },
            "PARTIAL",
        ),
        ({"log": "PARTIAL x.json " + "0" * 64 + "\n====== 1 failed in 3.00s ======\n"}, "passed"),
    ],
)
def test_load_record_refuses_bad_artifacts(
    kw: dict[str, Any], match: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    with pytest.raises(ValueError, match=match):
        cal.load_record(_rate_from(tmp_path, cal, cov, **kw))


@pytest.mark.parametrize(
    "summary",
    [
        "====== 1 failed, 1 passed in 1100.60s ======",
        "====== 2 passed in 1100.60s ======",
        "====== 1 passed, 1 skipped in 1100.60s ======",
        "====== 1 passed, 1 xfailed in 1100.60s ======",
        "====== 1 passed, 1 error in 1100.60s ======",
        "====== 1 passed in 1100.60s ======\n====== 1 passed in 1100.60s ======",
        "no summary at all",
    ],
)
def test_load_record_refuses_non_passing_logs(
    summary: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    part, _ = make_artifacts(tmp_path / "ref", cov)
    digest = json.loads(part.read_text())["content_sha256"]
    log = f"PARTIAL {part.name} {digest}\n{summary}\n"
    with pytest.raises(ValueError, match="summary|passed"):
        cal.load_record(_rate_from(tmp_path, cal, cov, log=log))
    log = (
        f"PARTIAL {part.name} {digest}\nFAILED tests/x.py::t\n====== 1 passed in 1100.60s ======\n"
    )
    with pytest.raises(ValueError, match="FAILED"):
        cal.load_record(_rate_from(tmp_path, cal, cov, log=log))


@pytest.mark.parametrize(
    ("mutate", "match"),
    [
        ("missing part", "part_record .* is not a file"),
        ("missing log", "part_log .* is not a file"),
        ("part hash", "part_record .*sha256 differs"),
        ("log hash", "part_log .*sha256 differs"),
        ("reps", "reps"),
        ("wall", "wall_mc_s"),
        ("total", "wall_total_s"),
        ("git sha", "run_sha"),
        ("dir", "is not a file"),
        ("path name", "plain file name"),
    ],
)
def test_load_record_refuses_inconsistent_record_vs_artifacts(
    mutate: str, match: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    path = _rate_from(tmp_path, cal, cov)
    rec = json.loads(path.read_text())
    art = Path(rec["log_dir"])
    if mutate == "missing part":
        (art / rec["part_record"]).unlink()
    elif mutate == "missing log":
        (art / rec["part_log"]).unlink()
    elif mutate == "part hash":
        rec["part_record_sha256"] = "0" * 64
    elif mutate == "log hash":
        rec["part_log_sha256"] = "0" * 64
    elif mutate == "reps":
        rec["reps"] = 82
        rec["s_per_rep"] = rec["wall_mc_s"] / 82
    elif mutate == "wall":
        rec["wall_mc_s"] += 100.0
        rec["s_per_rep"] = rec["wall_mc_s"] / rec["reps"]
    elif mutate == "total":
        rec["wall_total_s"] += 1.0
    elif mutate == "git sha":
        rec["git_sha"] = "b" * 40
    elif mutate == "dir":
        rec["log_dir"] = str(tmp_path / "nowhere")
    elif mutate == "path name":
        rec["part_log"] = "../art/part.log"
    path.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match=match):
        cal.load_record(path)


def test_build_record_refuses_unverifiable_inputs(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    part, lg = make_artifacts(tmp_path / "ok", cov)
    for kw, match in (({"unseal": True}, "not sealed"), ({"sha": "b" * 40}, "run_sha"),
                      ({"over": {"suite": "lv5c"}}, "suite")):  # fmt: skip
        bad, bad_log = make_artifacts(tmp_path / f"bad{len(match)}{match[:2]}", cov, **kw)
        with pytest.raises(ValueError, match=match):
            cal.build_record(bad, bad_log, bad.parent, SHA)
    with pytest.raises(ValueError, match="lie in log_dir"):
        cal.build_record(part, lg, tmp_path, SHA)
    with pytest.raises(ValueError, match="git_sha"):
        cal.build_record(part, lg, part.parent, "abc")
    plain = lg.parent / "plain.json"  # any JSON with n and wall_mc_s is no longer accepted
    plain.write_text(json.dumps({"n": 83, "wall_mc_s": 996.6}))
    with pytest.raises(ValueError, match="not sealed"):
        cal.build_record(plain, lg, lg.parent, SHA)


def test_declared_plans_refuses_a_plain_dict_and_a_changed_file(
    s5d: ModuleType, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    for fake in ({"record": "validated by v7r2_cal_rate.load_record"},
                 {"record": {}, "sha256": "0" * 64, "path": str(tmp_path / "none.json")}, 3):  # fmt: skip
        with pytest.raises(SystemExit, match="not ok"):
            s5d.declared_plans("lv5d", cal_rate=fake)  # type: ignore[arg-type]
    path = _rate_from(tmp_path, cal, cov)
    got = cal.load_record(path)
    forged = {**got, "record": {**got["record"], "s_per_rep": 1.0}}
    with pytest.raises(SystemExit, match="differs"):
        s5d.declared_plans("lv5d", cal_rate=forged)
    assert s5d.declared_plans("lv5d", cal_rate=got)
    s5d.declared_plans("lv5c", cal_rate=None)  # the frozen suites are untouched by the gate


# -- --print-hashed-source-digest, generation names, builder ----------------------------------------------------
def test_print_hashed_source_digest_accepts_lv5d_and_hr5d(
    rs: ModuleType, capsys: pytest.CaptureFixture[str]
) -> None:
    for suite in ("lv5d", "hr5d", "lv5c"):
        assert rs.main(["--print-hashed-source-digest", suite]) == 0
        assert capsys.readouterr().out.strip() == rs.hashed_source_digest(suite)
    assert rs.main(["--print-hashed-source-digest"]) == 0  # the historical form: lv5c
    assert capsys.readouterr().out.strip() == rs.hashed_source_digest("lv5c")
    assert rs.main(["--suite", "lv5d", "--print-hashed-source-digest"]) == 0
    assert capsys.readouterr().out.strip() == rs.hashed_source_digest("lv5d")
    assert rs.hashed_source_digest("lv5d") != rs.hashed_source_digest("lv5c")
    with pytest.raises(SystemExit, match="unknown suite"):
        rs.main(["--print-hashed-source-digest", "lv5e"])


def test_rehearsal_shard_count_constant_matches_the_loader(rs: ModuleType) -> None:
    assert rs.V7R2_REHEARSAL_SHARDS == _load("v7r_rehearsal").R2_SHARDS == 4
    names = [n for n, _, _ in rs.suite_steps("lv5d", 1, 1.0)]
    for k in range(
        rs.V7R2_REHEARSAL_SHARDS
    ):  # the generation steps exist in lv5d under these names
        assert any(rs.step_base_name(n) == f"v7r2-s{k}" for n in names)


@pytest.fixture(scope="module")
def builder() -> ModuleType:
    return _load("v7r_rehearsal_fixture")


@pytest.fixture
def git_repo(
    tmp_path: Path, builder: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> tuple[Path, str]:
    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / "a.txt").write_text("x\n")
    env = ["-c", "user.name=t", "-c", "user.email=t@example.org", "-c", "commit.gpgsign=false"]
    for cmd in (["init", "-q"], ["add", "a.txt"], [*env, "commit", "-q", "-m", "init"]):
        subprocess.run(["git", *cmd], cwd=repo, check=True, capture_output=True)
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo, check=True, capture_output=True, text=True
    ).stdout.strip()
    monkeypatch.setattr(builder, "REPO", repo)
    return repo, head


def _archive(tmp_path: Path, s5c: ModuleType, head: str, shards: int = 4) -> Path:
    src = tmp_path / "src"
    build_v7r2_fixture(src, s5c, "0" * 64, code_sha=head, layout_shards=shards)
    arch = tmp_path / "archive"
    for k in range(shards):
        suite = arch / f"RUN-{k}" / "transport" / "lv5d-x"
        (suite / "samples").mkdir(parents=True)
        shutil.copyfile(src / f"v7r2-s{k}.json", suite / "samples" / f"v7r2-s{k}.json")
        shutil.copyfile(
            src / f"v7r2-s{k}-escaped_neutral.npz", suite / f"v7r2-s{k}-escaped_neutral.npz"
        )
    return arch


def _argv(arch: Path, out: Path, sha: str, shards: int = 4, *extra: str) -> list[str]:
    return ["--archive", str(arch), "--run-ids", *[f"RUN-{k}" for k in range(shards)], "--code-sha", sha,
            "--seed-base", "20505000", "--out", str(out), *extra]  # fmt: skip


def test_builder_layout_v7r2_produces_a_verifying_lv5d_fixture(
    tmp_path: Path, s5c: ModuleType, rs: ModuleType, builder: ModuleType, git_repo: tuple[Path, str]
) -> None:
    _, head = git_repo
    arch = _archive(tmp_path, s5c, head)
    out = tmp_path / "fixture"
    assert builder.main(_argv(arch, out, head, 4, "--layout", "v7r2")) == 0
    prov = json.loads((out / "PROVENANCE.json").read_text())
    assert prov["hashed_source_digest"] == rs.hashed_source_digest(
        "lv5d"
    )  # the lv5d digest, not the lv5c one
    assert prov["hashed_source_digest"] != rs.hashed_source_digest("lv5c")
    assert prov["seed_base"] == 20505000 and prov["code_sha"] == head and prov["clean_tree"] is True
    assert (
        sorted(prov["partial_digests"]) == [f"v7r2-s{k}.json" for k in range(4)]
        and len(prov["created_from"]) == 8
    )
    blocks, info = _load("v7r_rehearsal").load_rehearsal(out, s5c, rs, layout="v7r2")
    assert blocks.shape == (900, 20) and info["layout"] == "v7r2" and info["shards"] == 4
    with pytest.raises(SystemExit, match="exists"):
        builder.main(_argv(arch, out, head, 4, "--layout", "v7r2"))


def test_builder_layout_v7r2_refuses_a_wrong_shard_set_and_the_default_layout_ignores_v7r2_shards(
    tmp_path: Path, s5c: ModuleType, builder: ModuleType, git_repo: tuple[Path, str]
) -> None:
    _, head = git_repo
    arch = _archive(tmp_path, s5c, head, shards=3)
    out = tmp_path / "fixture"
    with pytest.raises(SystemExit, match="verification failed.*need exactly 4 shards"):
        builder.main(_argv(arch, out, head, 3, "--layout", "v7r2"))
    assert not out.exists()
    with pytest.raises(
        SystemExit, match=r"no v7r-s\{k\}.json partial"
    ):  # --layout v7r (default) sees no v7r-s shards
        builder.main(_argv(arch, out, head, 3))
    with pytest.raises(SystemExit):  # an unknown layout is refused by argparse
        builder.main(_argv(arch, out, head, 3, "--layout", "v7r3"))
