"""Wiring task T3 of the lv5d/hr5d suites (plan Amendment 22 rule 2, Amendment 20 (f)6 / (f)8; Codex REVIEW-fe47f290):
the full lv5d fixture preflight (absent / inconsistent / stale / generation-exempt / valid), the measured calibration
rate gate (absent / valid / malformed record, re-split), ``--print-hashed-source-digest`` for lv5d / hr5d, the
``v7r2-s{k}`` generation path and the builder's ``--layout v7r2``. Unit tests only: no host run, no seed drawn."""

# ruff: noqa: E501

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import platform
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
RUN_ID = "RUN-20261101-r4b-0001"
VAL_ID = "VAL-20261101-120000-ABCDEF"


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
                   unseal: bool = False, log: str | None = None, drop: tuple[str, ...] = ()) -> tuple[Path, Path]:  # fmt: skip
    """A genuine sealed V7-R2 calibration part record (sealed by the real ``_seal``) and a synthetic but well-formed
    passing pytest log of that one part (a controlled host run is not needed to test the verification)."""
    d.mkdir(parents=True, exist_ok=True)
    n = 83 if part < 2 else 82
    rec: dict[str, Any] = {
        "schema": "ionmc-v7r2-calibration-record-1", "surrogate": surrogate, "part": part, "n": n,
        "reps_env": 3200, "layout": {"replicates": 10800, "reps_per_shard": 225, "parts": 39},
        "run_sha": sha, "suite": "lv5d", "wall_mc_s": wall, "wall_per_replicate_s": wall / n,
        "host": {"node": platform.node() or "node", "platform": platform.platform(), "machine": platform.machine(),
                 "cpu_count": os.cpu_count(), "python": platform.python_version()},
        "host_run_id": RUN_ID,
    }  # fmt: skip
    rec.update(over or {})
    for k in drop:
        rec.pop(k)
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
    rec = cal.build_record(
        part, lg, part.parent, SHA, RUN_ID, VAL_ID
    )  # host copied from the sealed record
    rec.update(over)
    path.write_text(json.dumps(rec))
    return path


def _plan_for(cal: ModuleType, rate: Path | str, *, name: str = "plan.md", **pin_over: str) -> Path:
    """A temporary acceptance plan (next to the record) pinning the record as it is NOW on disk (the genuine pin the
    orchestrator would write after verifying the protected records); ``pin_over`` replaces pin fields. Unparseable or
    incomplete records give a plan without a pin."""
    rate = Path(rate)
    plan = rate.parent / name
    lines = ["# temporary acceptance plan", ""]
    try:
        rec = json.loads(rate.read_text())
        pin = {"sha256": cal.record_digest(rec), "git_sha": rec["git_sha"], "host_run_id": rec["host_run_id"],
               "validation_id": rec["validation_id"], "host_node": rec["host"]["node"], **pin_over}  # fmt: skip
        lines.append("RATE-RECORD-OF-RECORD: " + " ".join(f"{k}={v}" for k, v in pin.items()))
    except (ValueError, KeyError, TypeError, AttributeError, OSError):
        pass
    plan.write_text("\n".join(lines) + "\n")
    return plan


def _pinned_load(cal: ModuleType, rate: Path | str) -> dict[str, Any]:
    """``load_record`` of the record against a temporary plan holding its genuine pin (``plan_path`` argument)."""
    return cal.load_record(rate, plan_path=_plan_for(cal, rate))


@pytest.fixture
def pin_default(cal: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Any:
    """Make the REPOSITORY-plan path of the consumers (environment record) resolve to a temporary plan pinning the given
    record: the module constant is patched in the test only (production code takes no plan from the environment)."""
    import importlib

    real = importlib.import_module("v7r2_cal_rate")  # the module the consumers imported by name

    def go(rate: Path) -> Path:
        plan = _plan_for(cal, rate)
        for m in (real, cal):
            monkeypatch.setattr(m, "ACCEPTANCE_PLAN", plan)
        return plan

    return go


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
    pin_default: Any,
) -> None:
    build_v7r2_fixture(lv5d_tree, s5c, DIGEST)
    refusal = rs.rehearsal_preflight("lv5d")
    assert refusal and "calibration rate not measured" in refusal  # fixture ok, rate not measured
    monkeypatch.setenv(RATE_ENV, str(tmp_path / "missing.json"))
    assert "not ok" in (rs.rehearsal_preflight("lv5d") or "")  # a path that is no file
    rate = write_rate(tmp_path / "rate.json", cal, cov)
    monkeypatch.setenv(RATE_ENV, str(rate))
    assert "not pinned" in (
        rs.rehearsal_preflight("lv5d") or ""
    )  # genuine record, repository plan has no pin
    pin_default(rate)
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
    pin_default: Any,
) -> None:
    rec = write_rate(tmp_path / "rate.json", cal, cov, s_per_rep=12.0)
    monkeypatch.setenv(RATE_ENV, str(rec))
    assert (
        "not pinned" in s5d.rep_plan_table("lv5d")["reason"]
    )  # the repository plan carries no pin
    pin_default(rec)
    rep = s5d.rep_plan_table("lv5d")
    assert rep["ok"] and rep["rate_source"] == "measured" and rep["reason"] is None
    assert rep["rate_record"]["sha256"] == hashlib.sha256(rec.read_bytes()).hexdigest()
    assert (
        rep["rate_record"]["host"]["node"] == (platform.node() or "node")
        and rep["rate_record"]["git_sha"] == "a" * 40
    )
    cal_plan = rep["plans"]["v7r2 calibration (3200 reps x 39 parts per surrogate)"]
    assert (
        cal_plan["s_per_item"] == 12.0 and cal_plan["fixed_s"] == 70.0
    )  # the measured rate, not the 12.2 default
    assert s5d.declared_plans(
        "lv5d", cal_rate=_pinned_load(cal, rec), plan_path=str(_plan_for(cal, rec))
    )  # the gated table is returned with a record
    # a slower measured rate exceeds the 1100 s target of 39 parts: not ok, the re-split is stated
    slow = write_rate(tmp_path / "slow.json", cal, cov, s_per_rep=13.0)
    bad = s5d.rep_plan_table("lv5d", str(slow), str(_plan_for(cal, slow)))
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
    rep = s5d.rep_plan_table("lv5d", str(p), str(_plan_for(cal, p)))
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
            "--git-sha", SHA, "--host-run-id", RUN_ID, "--validation-id", VAL_ID, "--out", str(out)]  # fmt: skip
    assert cal.main(args) == 0
    got = _pinned_load(cal, out)["record"]
    assert abs(got["s_per_rep"] - 996.6 / 83) < 1e-12 and abs(got["fixed_s"] - 104.0) < 1e-6
    assert got["part_record_sha256"] == hashlib.sha256(part.read_bytes()).hexdigest()
    with pytest.raises(SystemExit, match="exists"):
        cal.main(args)
    log.write_text("2 failed in 3s\n")
    with pytest.raises(ValueError, match="summary"):
        cal.build_record(part, log, part.parent, SHA, RUN_ID, VAL_ID)
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
           "wall_total_s": total, "host": dict(doc.get("host") or {"node": "x", "platform": "p", "machine": "m", "cpu_count": 1, "python": "3"}), "host_run_id": RUN_ID, "validation_id": VAL_ID,
           "git_sha": SHA, "log_dir": str(part.parent),
           "part_record": part.name, "part_record_sha256": hashlib.sha256(part.read_bytes()).hexdigest(),
           "part_log": lg.name, "part_log_sha256": hashlib.sha256(lg.read_bytes()).hexdigest()}  # fmt: skip
    out = tmp_path / "rate.json"
    out.write_text(json.dumps(rec))
    return out


def test_genuine_artifacts_load_and_authorise_the_39_part_plan(
    s5d: ModuleType, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    got = _pinned_load(cal, _rate_from(tmp_path, cal, cov))
    plan = str(_plan_for(cal, got["path"]))
    assert cal.require_validated(got, plan) == got
    assert s5d.declared_plans("lv5d", cal_rate=got, plan_path=plan)[
        "v7r2 reference (1e6, block rate)"
    ]
    assert s5d.rep_plan_table("lv5d", got["path"], plan)["ok"]


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
        _pinned_load(cal, _rate_from(tmp_path, cal, cov, **kw))


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
        _pinned_load(cal, _rate_from(tmp_path, cal, cov, log=log))
    log = (
        f"PARTIAL {part.name} {digest}\nFAILED tests/x.py::t\n====== 1 passed in 1100.60s ======\n"
    )
    with pytest.raises(ValueError, match="FAILED"):
        _pinned_load(cal, _rate_from(tmp_path, cal, cov, log=log))


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
        _pinned_load(cal, path)


def test_build_record_refuses_unverifiable_inputs(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    part, lg = make_artifacts(tmp_path / "ok", cov)
    for kw, match in (({"unseal": True}, "not sealed"), ({"sha": "b" * 40}, "run_sha"),
                      ({"over": {"suite": "lv5c"}}, "suite")):  # fmt: skip
        bad, bad_log = make_artifacts(tmp_path / f"bad{len(match)}{match[:2]}", cov, **kw)
        with pytest.raises(ValueError, match=match):
            cal.build_record(bad, bad_log, bad.parent, SHA, RUN_ID, VAL_ID)
    with pytest.raises(ValueError, match="lie in log_dir"):
        cal.build_record(part, lg, tmp_path, SHA, RUN_ID, VAL_ID)
    with pytest.raises(ValueError, match="git_sha"):
        cal.build_record(part, lg, part.parent, "abc", RUN_ID, VAL_ID)
    plain = lg.parent / "plain.json"  # any JSON with n and wall_mc_s is no longer accepted
    plain.write_text(json.dumps({"n": 83, "wall_mc_s": 996.6}))
    with pytest.raises(ValueError, match="not sealed"):
        cal.build_record(plain, lg, lg.parent, SHA, RUN_ID, VAL_ID)


def test_declared_plans_refuses_a_plain_dict_and_a_changed_file(
    s5d: ModuleType, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    path = _rate_from(tmp_path, cal, cov)
    plan = _plan_for(cal, path)
    for fake in ({"record": "validated by v7r2_cal_rate.load_record"},
                 {"record": {}, "sha256": "0" * 64, "path": str(tmp_path / "none.json")}, 3):  # fmt: skip
        with pytest.raises(SystemExit, match="not ok"):
            s5d.declared_plans("lv5d", cal_rate=fake, plan_path=str(plan))  # type: ignore[arg-type]
    got = _pinned_load(cal, path)
    forged = {**got, "record": {**got["record"], "s_per_rep": 1.0}}
    with pytest.raises(SystemExit, match="differs"):
        s5d.declared_plans("lv5d", cal_rate=forged, plan_path=str(plan))
    assert s5d.declared_plans("lv5d", cal_rate=got, plan_path=str(plan))
    s5d.declared_plans("lv5c", cal_rate=None)  # the frozen suites are untouched by the gate


# -- Codex REVIEW-323f4c12: the measured-rate record is conformant only when pinned in the acceptance plan --------
PIN_FIELDS = ("sha256", "git_sha", "host_run_id", "validation_id", "host_node")


def _pin_text(**f: str) -> str:
    return "RATE-RECORD-OF-RECORD: " + " ".join(f"{k}={f[k]}" for k in PIN_FIELDS)


def _genuine_pin(cal: ModuleType, rate: Path) -> dict[str, str]:
    rec = json.loads(rate.read_text())
    return {"sha256": cal.record_digest(rec), "git_sha": rec["git_sha"], "host_run_id": rec["host_run_id"],
            "validation_id": rec["validation_id"], "host_node": rec["host"]["node"]}  # fmt: skip


def test_pinned_record_is_conformant_and_unpinned_is_not(
    s5d: ModuleType, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    rate = _rate_from(tmp_path, cal, cov)
    plan = tmp_path / "plan.md"
    plan.write_text("# plan\nno pin here\n")
    with pytest.raises(
        ValueError, match="not pinned in the acceptance plan .*nonconformant until the orchestrator"
    ):
        cal.load_record(rate, plan_path=plan)
    plan.write_text("# plan\n" + _pin_text(**_genuine_pin(cal, rate)) + "\n")
    got = cal.load_record(rate, plan_path=plan)
    assert cal.require_validated(got, plan) == got
    assert cal.pin_line(got["record"]) == _pin_text(**_genuine_pin(cal, rate))
    assert cal.record_digest(got["record"]) == hashlib.sha256(
        json.dumps(got["record"], sort_keys=True, separators=(",", ":")).encode()).hexdigest()  # fmt: skip
    assert s5d.rep_plan_table("lv5d", str(rate), str(plan))["ok"]
    # the pin of a record that has been tampered with (any field) no longer matches
    rec = json.loads(rate.read_text())
    rec["log_dir"] = str(Path(rec["log_dir"]))  # same value: still pinned
    rate.write_text(json.dumps(rec))
    assert cal.load_record(rate, plan_path=plan)
    plan.write_text("# plan\n")  # the pin is removed
    with pytest.raises(ValueError, match="not pinned"):
        cal.require_validated(got, plan)
    with pytest.raises(SystemExit, match="not pinned"):
        s5d.declared_plans("lv5d", cal_rate=got, plan_path=str(plan))
    refused = s5d.rep_plan_table("lv5d", str(rate), str(plan))
    assert refused["ok"] is False and "not pinned" in refused["reason"] and refused["plans"] == {}


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("sha256", "0" * 64),
        ("git_sha", "b" * 40),
        ("host_run_id", "RUN-20261101-r4b-9999"),
        ("validation_id", "VAL-20261101-120000-FFFFFF"),
        ("host_node", "othernode"),
    ],
)
def test_pin_with_a_wrong_field_is_refused(
    field: str, value: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    rate = _rate_from(tmp_path, cal, cov)
    pin = {**_genuine_pin(cal, rate), field: value}
    plan = tmp_path / "plan.md"
    plan.write_text(_pin_text(**pin) + "\n")
    with pytest.raises(ValueError, match=f"does not match the record .*{field}"):
        cal.load_record(rate, plan_path=plan)


def test_last_pin_wins_and_a_second_pin_supersedes_the_first(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    rate = _rate_from(tmp_path, cal, cov)
    good = _pin_text(**_genuine_pin(cal, rate))
    other = _pin_text(**{**_genuine_pin(cal, rate), "sha256": "1" * 64})
    plan = tmp_path / "plan.md"
    plan.write_text(f"{other}\ntext\n{good}\n")
    assert cal.load_record(rate, plan_path=plan)  # the newest (last) pin names the record
    plan.write_text(f"{good}\ntext\n{other}\n")
    with pytest.raises(ValueError, match="does not match the record .*sha256"):
        cal.load_record(rate, plan_path=plan)  # an older pin does not keep the record conformant


def test_pin_grammar_is_strict_and_the_angle_bracket_form_is_no_pin(
    cal: ModuleType, tmp_path: Path
) -> None:
    sha, git, val = "a" * 64, "b" * 40, "VAL-20261101-120000-ABCDEF"
    good = f"RATE-RECORD-OF-RECORD: sha256={sha} git_sha={git} host_run_id=RUN-1234 validation_id={val} host_node=n1"
    plan = tmp_path / "p.md"
    plan.write_text(good + "\n")
    assert cal.plan_pin(plan) == {"sha256": sha, "git_sha": git, "host_run_id": "RUN-1234", "validation_id": val,
                                  "host_node": "n1"}  # fmt: skip
    for bad in (
        "RATE-RECORD-OF-RECORD: sha256=<64 hex> git_sha=<40 hex> host_run_id=<id> validation_id=<VAL-...> host_node=<node>",
        " " + good,  # not a whole line from the first column
        "see `" + good + "`",
        good.replace(sha, "A" * 64),
        good.replace(git, "b" * 39),
        good.replace(val, "VAL-1"),
        good + " extra=1",
        good.replace("RATE-RECORD-OF-RECORD:", "RATE-RECORD:"),
    ):
        plan.write_text(bad + "\n")
        assert cal.plan_pin(plan) is None, bad
    with pytest.raises(ValueError, match="unreadable"):
        cal.plan_pin(tmp_path / "missing.md")


def test_the_repository_plan_is_the_default_and_refuses_a_synthetic_record(
    s5d: ModuleType,
    cal: ModuleType,
    cov: ModuleType,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    assert cal.ACCEPTANCE_PLAN == REPO / "validation" / "plans" / "v3-005-acceptance.md"
    assert cal.ACCEPTANCE_PLAN.is_file()
    assert (
        cal.plan_pin() is None
    )  # no measurement exists: the committed plan must carry no pin (yet)
    rate = _rate_from(
        tmp_path, cal, cov
    )  # fully self-consistent synthetic artifacts, well-formed ids
    with pytest.raises(ValueError, match="not pinned"):
        cal.load_record(rate)
    monkeypatch.setenv(RATE_ENV, str(rate))
    with pytest.raises(ValueError, match="not pinned"):
        cal.record_from_env()
    assert "not pinned" in s5d.rep_plan_table("lv5d")["reason"]
    # a plan is never taken from the environment
    monkeypatch.setenv("IONMC_V7R2_ACCEPTANCE_PLAN", str(_plan_for(cal, rate)))
    with pytest.raises(ValueError, match="not pinned"):
        cal.load_record(rate)


def test_the_committed_plan_note_shows_the_grammar_without_being_a_pin(cal: ModuleType) -> None:
    text = cal.ACCEPTANCE_PLAN.read_text()
    assert "REVIEW-323f4c12" in text and "rate-record pin" in text
    assert "RATE-RECORD-OF-RECORD: sha256=<64 hex>" in text  # the grammar is shown ...
    assert cal.plan_pin() is None  # ... and is no pin


def test_cli_print_pin_line_prints_the_exact_line_and_never_writes_the_plan(
    cal: ModuleType, cov: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    part, lg = make_artifacts(tmp_path / "art", cov)
    out = tmp_path / "rate.json"
    plan_before = cal.ACCEPTANCE_PLAN.read_bytes()
    assert cal.main(["--part-record", str(part), "--part-log", str(lg), "--log-dir", str(part.parent),
                     "--git-sha", SHA, "--host-run-id", RUN_ID, "--validation-id", VAL_ID, "--out", str(out),
                     "--print-pin-line"]) == 0  # fmt: skip
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[-1] == _pin_text(**_genuine_pin(cal, out)) and lines[-1] == cal.pin_line(
        json.loads(out.read_text())
    )
    assert cal.ACCEPTANCE_PLAN.read_bytes() == plan_before
    plan = tmp_path / "plan.md"
    plan.write_text("note\n" + lines[-1] + "\n")  # copying the printed line opens the gate
    assert cal.load_record(out, plan_path=plan)
    out.unlink()
    assert cal.main(["--part-record", str(part), "--part-log", str(lg), "--log-dir", str(part.parent),
                     "--git-sha", SHA, "--host-run-id", RUN_ID, "--validation-id", VAL_ID,
                     "--out", str(out)]) == 0  # fmt: skip
    assert "RATE-RECORD-OF-RECORD" not in capsys.readouterr().out  # only on request


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


# -- Codex REVIEW-ffe77aac: derived planning fields, host identity and run attestation are bound -------------------------
def _reviewer_example(tmp_path: Path, cal: ModuleType, cov: ModuleType) -> Path:
    """wall_mc_s 996.6, reps 83, wall_total_s 1100.6: the correct fixed_s is 104 and the plan exceeds the target."""
    part, lg = make_artifacts(tmp_path / "art", cov, wall=996.6, total=1100.6)
    out = tmp_path / "rate.json"
    assert cal.main(["--part-record", str(part), "--part-log", str(lg), "--log-dir", str(part.parent),
                     "--git-sha", SHA, "--host-run-id", RUN_ID, "--validation-id", VAL_ID, "--out", str(out)]) == 0  # fmt: skip
    return out


def test_altered_fixed_s_is_refused_and_the_genuine_example_says_re_split(
    s5d: ModuleType, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    out = _reviewer_example(tmp_path, cal, cov)
    rec = json.loads(out.read_text())
    assert abs(rec["fixed_s"] - 104.0) < 1e-6
    genuine = s5d.rep_plan_table("lv5d", str(out), str(_plan_for(cal, out)))
    assert genuine["ok"] is False and "re-split" in genuine["reason"]
    rec["fixed_s"] = 70.0  # every artifact hash is preserved
    out.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="fixed_s .* differs"):
        _pinned_load(cal, out)
    refused = s5d.rep_plan_table("lv5d", str(out), str(_plan_for(cal, out)))
    assert refused["ok"] is False and refused["plans"] == {} and "fixed_s" in refused["reason"]
    rec["fixed_s"] = (
        1.0e6  # a pessimistic edit is refused as well: the value is derived, not chosen
    )
    out.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="fixed_s .* differs"):
        _pinned_load(cal, out)


def test_altered_s_per_rep_is_refused(cal: ModuleType, cov: ModuleType, tmp_path: Path) -> None:
    out = _rate_from(tmp_path, cal, cov)
    rec = json.loads(out.read_text())
    rec["s_per_rep"] *= 0.9
    out.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="s_per_rep .* differs"):
        _pinned_load(cal, out)


@pytest.mark.parametrize("field", ["node", "platform", "machine", "cpu_count", "python"])
def test_altered_host_field_is_refused(
    field: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    out = _rate_from(tmp_path, cal, cov)
    rec = json.loads(out.read_text())
    rec["host"][field] = rec["host"][field] + 1 if field == "cpu_count" else "slowhost"
    out.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="host .* differs from the sealed part record"):
        _pinned_load(cal, out)
    rec["host"].pop(field)
    out.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="host must be an object with exactly"):
        _pinned_load(cal, out)


def test_host_is_taken_from_the_sealed_record_by_the_builder(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    host = {
        "node": "slow-node",
        "platform": "Linux-x",
        "machine": "x86_64",
        "cpu_count": 4,
        "python": "3.12.1",
    }
    part, lg = make_artifacts(tmp_path / "art", cov, over={"host": host})
    assert cal.build_record(part, lg, part.parent, SHA, RUN_ID, VAL_ID)["host"] == host


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"drop": ("host",)}, "sealed host identity unusable"),
        ({"over": {"host": {"node": "n"}}}, "host"),
        (
            {
                "over": {
                    "host": {
                        "node": "",
                        "platform": "p",
                        "machine": "m",
                        "cpu_count": 1,
                        "python": "3",
                    }
                }
            },
            "node",
        ),
        ({"drop": ("host_run_id",)}, "sealed host identity unusable"),
        ({"over": {"host_run_id": "unknown"}}, "sealed host identity unusable"),
        ({"over": {"host_run_id": "short"}}, "sealed host identity unusable"),
    ],
)
def test_sealed_part_record_without_usable_host_identity_is_refused(
    kw: dict[str, Any], match: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    with pytest.raises(ValueError, match=match):
        _pinned_load(cal, _rate_from(tmp_path, cal, cov, **kw))
    part, lg = make_artifacts(tmp_path / "b", cov, **kw)
    with pytest.raises(ValueError, match=match):
        cal.build_record(part, lg, part.parent, SHA, RUN_ID, VAL_ID)


@pytest.mark.parametrize(
    ("edit", "match"),
    [
        (lambda r: r.pop("host_run_id"), "record fields differ"),
        (lambda r: r.update(host_run_id="unknown"), "host_run_id must be"),
        (lambda r: r.update(host_run_id=""), "host_run_id must be"),
        (
            lambda r: r.update(host_run_id="RUN-OTHER-00000002"),
            "host_run_id .* differs from the sealed",
        ),
        (lambda r: r.pop("validation_id"), "record fields differ"),
        (lambda r: r.update(validation_id=""), "validation_id must be"),
        (lambda r: r.update(validation_id="VAL-2026-1101-abcdef"), "validation_id must be"),
        (lambda r: r.update(validation_id="VAL-20261101-120000-abcdef"), "validation_id must be"),
        (lambda r: r.update(validation_id=None), "validation_id must be"),
        (lambda r: r.update(planning_note="x"), "undeclared"),
        (lambda r: r.update(parts_authorised=39), "undeclared"),
    ],
)
def test_rate_record_ids_and_undeclared_fields(
    edit: Any, match: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    out = _rate_from(tmp_path, cal, cov)
    assert _pinned_load(cal, out)
    rec = json.loads(out.read_text())
    edit(rec)
    out.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match=match):
        _pinned_load(cal, out)


def test_build_record_and_cli_refuse_missing_or_unknown_ids(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    part, lg = make_artifacts(tmp_path / "art", cov)
    for ids in (
        ("unknown", VAL_ID),
        ("", VAL_ID),
        (RUN_ID, ""),
        (RUN_ID, "VAL-bad"),
        (RUN_ID, "unknown"),
    ):
        with pytest.raises(ValueError, match="host_run_id|validation_id"):
            cal.build_record(part, lg, part.parent, SHA, *ids)
    base = ["--part-record", str(part), "--part-log", str(lg), "--log-dir", str(part.parent), "--git-sha", SHA,
            "--out", str(tmp_path / "r.json")]  # fmt: skip
    for extra in ([], ["--host-run-id", RUN_ID], ["--validation-id", VAL_ID]):
        with pytest.raises(SystemExit):  # argparse: both ids are required
            cal.main(base + extra)
    with pytest.raises(SystemExit, match="not ok"):
        cal.main(base + ["--host-run-id", "unknown", "--validation-id", VAL_ID])
    assert not (tmp_path / "r.json").exists()


def test_calibration_part_writer_seals_the_host_identity() -> None:
    """The part writer (``test_v7r2_calibration._host_identity``) seals exactly the block the loader requires."""
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(
        "test_v7r2_calibration_t3", REPO / "tests" / "ionmc" / "test_v7r2_calibration.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    cal = _load("v7r2_cal_rate")
    got = mod._host_identity()
    assert set(got["host"]) == set(cal.HOST_FIELDS) and got["host"]["node"] == platform.node()
    cal.validate_host(got["host"])
    assert got["host_run_id"] == os.environ.get("IONMC_HOST_RUN_ID", "unknown")
