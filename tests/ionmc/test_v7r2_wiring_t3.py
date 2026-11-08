"""Wiring task T3 of the lv5d/hr5d suites (plan Amendment 22 rule 2, Amendment 20 (f)6 / (f)8; Codex REVIEW-fe47f290):
the full lv5d fixture preflight (absent / inconsistent / stale / generation-exempt / valid), the measured calibration
rate gate (absent / valid / malformed record, re-split; the signed RUN- / VAL- attestation chain of V3-INFRA-002,
Codex REVIEW-61f4b0e3 / REVIEW-cf3f019e), ``--print-hashed-source-digest`` for lv5d / hr5d, the ``v7r2-s{k}``
generation path and the builder's ``--layout v7r2``. Unit tests only: no host run, no seed drawn. The signed record
copies are synthetic: they are signed in the test with throw-away Ed25519 keys (``infrastructure/experiment_v3/attest.py``)
listed in a temporary trust file, never with the service keys."""

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
import types
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "validation" / "scripts" / "transport"
COMMITTED_TRUST = REPO / "infrastructure" / "keys" / "trusted_attestation_keys.json"
R4B_DIR = REPO / "validation" / "results" / "transport" / "v7r2-r4b"
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


# -- throw-away attestation keys and the temporary trust file (V3-INFRA-002 record format) ------------------------
def _attest_module() -> ModuleType:
    """``infrastructure/experiment_v3/attest.py`` (the signing half; tests may import it, the script may not)."""
    spec = importlib.util.spec_from_file_location(
        "ionmc_attest_t3", REPO / "infrastructure" / "experiment_v3" / "attest.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses resolve the module by name
    spec.loader.exec_module(mod)
    return mod


ATTEST = _attest_module()
KEYS = {
    "validation-manager": ATTEST.keypair_from_seed(
        hashlib.sha256(b"t3 validation-manager").digest()
    ),
    "host-runner": ATTEST.keypair_from_seed(hashlib.sha256(b"t3 host-runner").digest()),
    "rogue": ATTEST.keypair_from_seed(hashlib.sha256(b"t3 untrusted").digest()),  # listed nowhere
}


@pytest.fixture(scope="module")
def trust_file(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A trust file listing the test keys of both signers (the rogue key is not listed)."""
    p = tmp_path_factory.mktemp("trust") / "trusted_attestation_keys.json"
    p.write_text(json.dumps({s: [{"key_id": KEYS[s].key_id, "public_key": KEYS[s].public_key.hex(),
                                  "added": "2026-11-08", "note": "test key"}]
                             for s in ("validation-manager", "host-runner")}, indent=2))  # fmt: skip
    return p


@pytest.fixture(autouse=True)
def _test_trust(cal: ModuleType, trust_file: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The default trust file of the module instances under test (``cal`` and the one the consumers import by name) is
    the temporary one; tests of the committed state pass ``trust_path`` explicitly."""
    import importlib

    for m in (cal, importlib.import_module("v7r2_cal_rate")):
        monkeypatch.setattr(m, "TRUST_FILE", trust_file)


def _sign(doc: dict[str, Any], signer: str, key: str | None = None) -> dict[str, Any]:
    return ATTEST.sign_record(doc, KEYS[key or signer], signer)


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
RUN_ID = "RUN-20261101T120000Z-0123abcd"  # the controlled host-runner run id
TOKEN = "R4B-20261101T120000Z-aaaaaaaa"  # the in-run token (IONMC_HOST_RUN_ID of the wrapper), sealed in the part record
VAL_ID = "VAL-20261101-120000-ABCDEF"
VAL_FILE = f"{VAL_ID}.json"
RUN_FILE = f"{RUN_ID}.json"
MANIFEST = (
    "e" * 64
)  # the signed output manifest digest of the run (claimed by the validation record)


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


def _artifact_lines(part: Path, lg: Path) -> list[str]:
    return [f"CAL-RATE-ARTIFACT part_record {part.name} {hashlib.sha256(part.read_bytes()).hexdigest()}",
            f"CAL-RATE-ARTIFACT part_log {lg.name} {hashlib.sha256(lg.read_bytes()).hexdigest()}"]  # fmt: skip


def make_artifacts(d: Path, cov: ModuleType, *, wall: float = 996.6, part: int = 0, surrogate: str = "archive_fitted",
                   sha: str = SHA, total: float | None = None, over: dict[str, Any] | None = None,
                   unseal: bool = False, log: str | None = None, drop: tuple[str, ...] = (), val: bool = True,
                   val_over: dict[str, Any] | None = None, val_det: dict[str, Any] | None = None,
                   val_drop: tuple[str, ...] = (), val_key: str | None = None, val_tamper: dict[str, Any] | None = None,
                   run_over: dict[str, Any] | None = None, run_drop: tuple[str, ...] = (), run_key: str | None = None,
                   run_tamper: dict[str, Any] | None = None, stdout: str | None = None,
                   argv: list[str] | None = None) -> tuple[Path, Path]:  # fmt: skip
    """A genuine sealed V7-R2 calibration part record (sealed by the real ``_seal``) and a synthetic but well-formed
    passing pytest log of that one part, and (``val``) signed copies of the host-runner record (``RUN_FILE``: token line
    and artifact lines in its stdout) and of the validation-manager record (``VAL_FILE``: passed, ``details.host_runs``
    naming the run and its manifest). ``*_over`` / ``val_det`` replace fields before signing, ``*_drop`` remove them,
    ``*_key`` signs with another key ("rogue": untrusted), ``*_tamper`` replaces fields AFTER signing."""
    d.mkdir(parents=True, exist_ok=True)
    n = 83 if part < 2 else 82
    rec: dict[str, Any] = {
        "schema": "ionmc-v7r2-calibration-record-1", "surrogate": surrogate, "part": part, "n": n,
        "reps_env": 3200, "layout": {"replicates": 10800, "reps_per_shard": 225, "parts": 39},
        "run_sha": sha, "suite": "lv5d", "wall_mc_s": wall, "wall_per_replicate_s": wall / n,
        "host": {"node": platform.node() or "node", "platform": platform.platform(), "machine": platform.machine(),
                 "cpu_count": os.cpu_count(), "python": platform.python_version()},
        "host_run_id": TOKEN,
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
    if val:
        out = stdout if stdout is not None else "\n".join(
            [f"token {TOKEN}", "kappa parts: 6 passed", *_artifact_lines(path, lg)]) + "\n"  # fmt: skip
        run: dict[str, Any] = {
            "schema_version": 1, "run_id": RUN_ID, "task_id": "V3-005B", "branch": "task/v3-005b", "sha": sha,
            "dirty_before": False, "dirty_after": False, "sha_after": sha, "committed_state_unchanged": True,
            "argv": argv if argv is not None else ["python", "-c", "chain wrapper"], "exit_code": 0,
            "timed_out": False, "stdout": out, "stderr": "", "stdout_truncated": False, "stderr_truncated": False,
            "stdout_total_chars": len(out), "output_paths": [], "output_publication_error": None, "succeeded": True,
            "stdout_sha256": hashlib.sha256(out.encode()).hexdigest(), "stderr_sha256": hashlib.sha256(b"").hexdigest(),
            "output_manifest_sha256": MANIFEST, **(run_over or {}),
        }  # fmt: skip
        for k in run_drop:
            run.pop(k)
        signed_run = {**_sign(run, "host-runner", run_key), **(run_tamper or {})}
        (d / RUN_FILE).write_text(json.dumps(signed_run, indent=2) + "\n")
        doc: dict[str, Any] = {
            "schema_version": 1, "validation_id": VAL_ID, "task_id": "V3-005B",
            "recorded_at_utc": "2026-11-01T12:00:00.000000+00:00", "suite": "host-r4b-calibration-rate-measurement",
            "status": "passed", "summary": "synthetic signed record", "git_sha": sha,
            "details": {"head_sha": sha, "host_runs": [{"run_id": RUN_ID, "output_manifest_sha256": MANIFEST}],
                        **(val_det or {})},
            **(val_over or {}),
        }  # fmt: skip
        for k in val_drop:
            doc.pop(k)
        signed_val = {**_sign(doc, "validation-manager", val_key), **(val_tamper or {})}
        (d / VAL_FILE).write_text(json.dumps(signed_val, indent=2) + "\n")
    return path, lg


def _build(
    cal: ModuleType, part: Path, lg: Path, log_dir: Path | None = None, **kw: Any
) -> dict[str, Any]:
    """``build_record`` with the signed record copies next to the part record."""
    args: dict[str, Any] = {"git_sha": SHA, "host_run_record": part.parent / RUN_FILE,
                            "validation_record": part.parent / VAL_FILE, **kw}  # fmt: skip
    return cal.build_record(part, lg, log_dir or part.parent, **args)


def _cli(part: Path, lg: Path, out: Path, *extra: str) -> list[str]:
    return ["--part-record", str(part), "--part-log", str(lg), "--log-dir", str(part.parent), "--git-sha", SHA,
            "--host-run-record", str(part.parent / RUN_FILE), "--validation-record", str(part.parent / VAL_FILE),
            "--out", str(out), *extra]  # fmt: skip


def write_rate(
    path: Path, cal: ModuleType, cov: ModuleType, *, s_per_rep: float = 12.0, **over: Any
) -> Path:
    """A measured-rate record built by ``build_record`` from genuine artifacts (``over``: fields replaced afterwards)."""
    part, lg = make_artifacts(path.parent / f"{path.stem}-art", cov, wall=s_per_rep * 83)
    rec = _build(
        cal, part, lg
    )  # host copied from the sealed record, attestation from the signed copies
    rec.update(over)
    path.write_text(json.dumps(rec))
    return path


PIN_FIELDS = ("sha256", "git_sha", "host_run_id", "validation_id", "host_run_record_sha256", "validation_record_sha256",
              "trust_file_sha256")  # fmt: skip


def _pin_of(cal: ModuleType, rec: dict[str, Any]) -> dict[str, str]:
    """The genuine pin fields of a record (the order of the pin line)."""
    return {"sha256": cal.record_digest(rec), **{k: rec[k] for k in PIN_FIELDS[1:]}}


def _plan_for(cal: ModuleType, rate: Path | str, *, name: str = "plan.md", **pin_over: str) -> Path:
    """A temporary acceptance plan (next to the record) pinning the record as it is NOW on disk (the genuine pin the
    orchestrator would write after verifying the protected records); ``pin_over`` replaces pin fields. Unparseable or
    incomplete records give a plan without a pin."""
    rate = Path(rate)
    plan = rate.parent / name
    lines = ["# temporary acceptance plan", ""]
    try:
        rec = json.loads(rate.read_text())
        pin = {**_pin_of(cal, rec), **pin_over}
        lines.append("RATE-RECORD-OF-RECORD: " + " ".join(f"{k}={pin[k]}" for k in PIN_FIELDS))
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
    msg = rs.rehearsal_preflight("lv5d") or ""
    assert (
        "not pinned" in msg or "does not match the record" in msg
    )  # genuine synthetic record, not pinned in the repository plan
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
    reason = s5d.rep_plan_table("lv5d")["reason"]
    assert (
        "not pinned" in reason or "does not match the record" in reason
    )  # the repository plan does not pin it
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
    args = _cli(part, log, out)
    assert cal.main(args) == 0
    got = _pinned_load(cal, out)["record"]
    assert abs(got["s_per_rep"] - 996.6 / 83) < 1e-12 and abs(got["fixed_s"] - 104.0) < 1e-6
    assert got["part_record_sha256"] == hashlib.sha256(part.read_bytes()).hexdigest()
    with pytest.raises(SystemExit, match="exists"):
        cal.main(args)
    log.write_text("2 failed in 3s\n")
    with pytest.raises(ValueError, match="summary"):
        _build(cal, part, log)
    out.unlink()
    with pytest.raises(SystemExit, match="not ok"):
        cal.main(args)
    assert not out.exists()


# -- Codex REVIEW-6519d4ff: the measured-rate record is only accepted with verified controlled-run artifacts ------------
def _sig(path: Path) -> dict[str, Any]:
    block = json.loads(path.read_text()).get("attestation") or {}
    return {k: block.get(k) for k in ("key_id", "payload_sha256", "signature")}


def _rate_from(tmp_path: Path, cal: ModuleType, cov: ModuleType, **kw: Any) -> Path:
    """Build a record from artifacts made with ``kw`` through the unchecked path (artifacts are valid only if ``kw`` is
    empty): the record is assembled by hand from the files so that ``load_record`` is what refuses."""
    part, lg = make_artifacts(tmp_path / "art", cov, **kw)
    d = part.parent
    doc = json.loads(part.read_text())
    wall, n = doc["wall_mc_s"], doc["n"]
    total = wall + 70.0

    def sha(p: Path) -> str:
        return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else "0" * 64

    rec = {"schema": cal.CAL_RATE_SCHEMA, "s_per_rep": wall / n, "fixed_s": 70.0, "reps": n, "wall_mc_s": wall,
           "wall_total_s": total, "host": dict(doc.get("host") or {"node": "x", "platform": "p", "machine": "m", "cpu_count": 1, "python": "3"}),
           "git_sha": SHA, "log_dir": str(d),
           "part_record": part.name, "part_record_sha256": sha(part), "part_log": lg.name, "part_log_sha256": sha(lg),
           "host_run_id": RUN_ID, "host_run_record": RUN_FILE, "host_run_record_sha256": sha(d / RUN_FILE),
           "host_run_signature": _sig(d / RUN_FILE) if (d / RUN_FILE).exists() else {},
           "validation_id": VAL_ID, "validation_record": VAL_FILE, "validation_record_sha256": sha(d / VAL_FILE),
           "validation_signature": _sig(d / VAL_FILE) if (d / VAL_FILE).exists() else {},
           "trust_file_sha256": sha(cal.TRUST_FILE)}  # fmt: skip
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
            _build(cal, bad, bad_log)
    with pytest.raises(ValueError, match="lie in log_dir"):
        _build(cal, part, lg, tmp_path)
    with pytest.raises(ValueError, match="git_sha"):
        _build(cal, part, lg, git_sha="abc")
    plain = lg.parent / "plain.json"  # any JSON with n and wall_mc_s is no longer accepted
    plain.write_text(json.dumps({"n": 83, "wall_mc_s": 996.6}))
    with pytest.raises(ValueError, match="not sealed"):
        _build(cal, plain, lg)


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
def _pin_text(**f: str) -> str:
    return "RATE-RECORD-OF-RECORD: " + " ".join(f"{k}={f[k]}" for k in PIN_FIELDS)


def _genuine_pin(cal: ModuleType, rate: Path) -> dict[str, str]:
    return _pin_of(cal, json.loads(rate.read_text()))


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
        ("host_run_id", "RUN-20261101T120000Z-ffffffff"),
        ("validation_id", "VAL-20261101-120000-FFFFFF"),
        ("host_run_record_sha256", "0" * 64),
        ("validation_record_sha256", "0" * 64),
        ("trust_file_sha256", "0" * 64),
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
    run, tok, rsha, vsha, tsha = (
        "RUN-20261101T120000Z-0123abcd",
        "R4B-20261101T120000Z-aaaaaaaa",
        "f" * 64,
        "c" * 64,
        "d" * 64,
    )
    good = (f"RATE-RECORD-OF-RECORD: sha256={sha} git_sha={git} host_run_id={run} validation_id={val} "
            f"host_run_record_sha256={rsha} validation_record_sha256={vsha} trust_file_sha256={tsha}")  # fmt: skip
    plan = tmp_path / "p.md"
    plan.write_text(good + "\n")
    assert cal.plan_pin(plan) == {"sha256": sha, "git_sha": git, "host_run_id": run, "validation_id": val,
                                  "host_run_record_sha256": rsha, "validation_record_sha256": vsha,
                                  "trust_file_sha256": tsha}  # fmt: skip
    for bad in (
        "RATE-RECORD-OF-RECORD: sha256=<64 hex> git_sha=<40 hex> host_run_id=<RUN-...> validation_id=<VAL-...> "
        "host_run_record_sha256=<64 hex> validation_record_sha256=<64 hex> trust_file_sha256=<64 hex>",
        # the grammars of 2026-11-01 and 2026-11-07 (token / host node, no signed-copy digests) are no pins
        f"RATE-RECORD-OF-RECORD: sha256={sha} git_sha={git} host_run_id={tok} validation_id={val} host_node=n1",
        f"RATE-RECORD-OF-RECORD: sha256={sha} git_sha={git} host_run_id={run} validation_id={val} "
        f"host_run_token={tok} validation_record_sha256={vsha} host_node=n1",
        good.replace(run, tok),  # an in-run token is not a runner id
        good.replace(f" host_run_record_sha256={rsha}", ""),
        good.replace(f" validation_record_sha256={vsha}", ""),
        good.replace(f" trust_file_sha256={tsha}", ""),
        good.replace(vsha, "c" * 63),
        good + " host_node=n1",
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
    rate = _rate_from(
        tmp_path, cal, cov
    )  # fully self-consistent synthetic artifacts, signed by the test keys
    refused = "not pinned|does not match the record"
    with pytest.raises(ValueError, match=refused):
        cal.load_record(rate)  # the repository plan does not pin it
    monkeypatch.setenv(RATE_ENV, str(rate))
    with pytest.raises(ValueError, match=refused):
        cal.record_from_env()
    reason = s5d.rep_plan_table("lv5d")["reason"]
    assert "not pinned" in reason or "does not match the record" in reason
    # a plan is never taken from the environment
    monkeypatch.setenv("IONMC_V7R2_ACCEPTANCE_PLAN", str(_plan_for(cal, rate)))
    with pytest.raises(ValueError, match=refused):
        cal.load_record(rate)


def test_the_committed_rate_record_state_is_consistent(cal: ModuleType) -> None:
    """The committed plan and ``v7r2-r4b/cal-rate.json`` agree under the COMMITTED trust file. Until the re-measurement
    of 2026-11-08 is pinned, the plan holds only the superseded pins of the old grammars (history, no pin) and the
    token-attested record of 2026-11-07 is refused as legacy; once a pin of the signed grammar exists, the committed
    record must load against it with the committed keys."""
    text = cal.ACCEPTANCE_PLAN.read_text()
    assert (
        "host_run_token=R4B-20261102T010633Z-b7c44f51" in text
    )  # the superseded pin lines stay as history
    pin = cal.plan_pin()
    committed = R4B_DIR / "cal-rate.json"
    if pin is None:
        rec = json.loads(committed.read_text())
        assert rec.get("schema") == "ionmc-v7r2-cal-rate-1"
        with pytest.raises(ValueError, match="legacy token-attested"):
            cal.load_record(committed, trust_path=COMMITTED_TRUST)
        return
    got = cal.load_record(committed, trust_path=COMMITTED_TRUST)
    assert cal.record_digest(got["record"]) == pin["sha256"]
    assert (
        got["record"]["trust_file_sha256"]
        == hashlib.sha256(COMMITTED_TRUST.read_bytes()).hexdigest()
    )


def test_the_default_trust_file_is_the_committed_one(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    """Unpatched, the module verifies against ``infrastructure/keys/trusted_attestation_keys.json``, which does not list
    the test keys: a record signed by them is unverified there and refused (no service key is used in tests)."""
    fresh = _load("v7r2_cal_rate")
    assert fresh.TRUST_FILE == COMMITTED_TRUST
    rate = _rate_from(tmp_path, cal, cov)
    plan = _plan_for(cal, rate)
    assert cal.load_record(rate, plan_path=plan)  # the test trust file
    with pytest.raises(ValueError, match="trust_file_sha256|unverified"):
        cal.load_record(rate, plan_path=plan, trust_path=COMMITTED_TRUST)
    part, lg = make_artifacts(tmp_path / "b", cov)
    with pytest.raises(ValueError, match="unverified record rejected"):
        _build(cal, part, lg, trust_path=COMMITTED_TRUST)


def test_the_committed_plan_note_shows_the_grammar_without_being_a_pin(cal: ModuleType) -> None:
    text = cal.ACCEPTANCE_PLAN.read_text()
    assert "REVIEW-323f4c12" in text and "rate-record pin" in text
    assert "RATE-RECORD-OF-RECORD: sha256=<64 hex>" in text  # the grammar is shown ...
    assert (
        "trust_file_sha256=<64 hex>" in text
    )  # ... including the signed grammar of 2026-11-08 ...
    # ... and the placeholder lines are no pins (the parser only accepts hex digests)
    pin = cal.plan_pin()
    assert pin is None or all("<" not in v for v in pin.values())


def test_cli_print_pin_line_prints_the_exact_line_and_never_writes_the_plan(
    cal: ModuleType, cov: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    part, lg = make_artifacts(tmp_path / "art", cov)
    out = tmp_path / "rate.json"
    plan_before = cal.ACCEPTANCE_PLAN.read_bytes()
    assert cal.main(_cli(part, lg, out, "--print-pin-line")) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[-1] == _pin_text(**_genuine_pin(cal, out)) and lines[-1] == cal.pin_line(
        json.loads(out.read_text())
    )
    assert cal.ACCEPTANCE_PLAN.read_bytes() == plan_before
    plan = tmp_path / "plan.md"
    plan.write_text("note\n" + lines[-1] + "\n")  # copying the printed line opens the gate
    assert cal.load_record(out, plan_path=plan)
    out.unlink()
    assert cal.main(_cli(part, lg, out)) == 0
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
    assert cal.main(_cli(part, lg, out)) == 0
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
    assert _build(cal, part, lg)["host"] == host


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
        _build(cal, part, lg)


@pytest.mark.parametrize(
    ("edit", "match"),
    [
        (lambda r: r.pop("host_run_id"), "record fields differ"),
        (lambda r: r.update(host_run_id="unknown"), "host_run_id must be"),
        (lambda r: r.update(host_run_id=""), "host_run_id must be"),
        (
            lambda r: r.update(host_run_id=TOKEN),
            "host_run_id must be",
        ),  # an in-run token is no runner id
        (lambda r: r.update(host_run_id="RUN-20261101-r4b-0001"), "host_run_id must be"),
        (lambda r: r.update(host_run_id="RUN-20261101T120000Z-0123ABCD"), "host_run_id must be"),
        (
            lambda r: r.update(host_run_id="RUN-20261101T120000Z-ffffffff"),
            "host_run_id of the record differs",
        ),
        (lambda r: r.pop("host_run_record"), "record fields differ"),
        (lambda r: r.pop("host_run_record_sha256"), "record fields differ"),
        (lambda r: r.update(host_run_record=""), "host_run_record must be"),
        (lambda r: r.update(host_run_record_sha256="0" * 64), "host_run_record .*sha256 differs"),
        (lambda r: r.pop("host_run_signature"), "record fields differ"),
        (lambda r: r.update(host_run_signature={}), "host_run_signature must be"),
        (
            lambda r: r["host_run_signature"].update(signature="0" * 128),
            "host_run_signature of the record differs",
        ),
        (
            lambda r: r["host_run_signature"].update(key_id="0" * 16),
            "host_run_signature of the record differs",
        ),
        (
            lambda r: r["validation_signature"].update(payload_sha256="0" * 64),
            "validation_signature of the record differs",
        ),
        (
            lambda r: r["validation_signature"].update(signature="x"),
            "validation_signature.signature is malformed",
        ),
        (lambda r: r.pop("validation_record"), "record fields differ"),
        (lambda r: r.pop("validation_record_sha256"), "record fields differ"),
        (lambda r: r.update(validation_record=""), "validation_record must be"),
        (lambda r: r.update(validation_record_sha256="0" * 63), "validation_record_sha256 must be"),
        (lambda r: r.pop("validation_id"), "record fields differ"),
        (lambda r: r.update(validation_id=""), "validation_id must be"),
        (lambda r: r.update(validation_id="VAL-2026-1101-abcdef"), "validation_id must be"),
        (lambda r: r.update(validation_id="VAL-20261101-120000-abcdef"), "validation_id must be"),
        (
            lambda r: r.update(validation_id="VAL-20261101-120000-FFFFFF"),
            "validation_id of the record differs",
        ),
        (lambda r: r.update(validation_id=None), "validation_id must be"),
        (lambda r: r.pop("trust_file_sha256"), "record fields differ"),
        (
            lambda r: r.update(trust_file_sha256="0" * 64),
            "trust_file_sha256 of the record differs from the trust file",
        ),
        (lambda r: r.update(host_run_token=TOKEN), "legacy token-attested"),
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


def test_build_record_and_cli_refuse_missing_inputs(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    part, lg = make_artifacts(tmp_path / "art", cov)
    with pytest.raises(ValueError, match="must be files"):
        _build(cal, part, lg, host_run_record=part.parent / "absent.json")
    with pytest.raises(ValueError, match="must be files"):
        _build(cal, part, lg, validation_record=part.parent / "absent.json")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    for name, key in ((RUN_FILE, "host_run_record"), (VAL_FILE, "validation_record")):
        shutil.copy(part.parent / name, elsewhere / name)
        with pytest.raises(ValueError, match="lie in log_dir"):
            _build(cal, part, lg, **{key: elsewhere / name})
    with pytest.raises(
        ValueError, match="validation record .*unreadable|host run record .*unreadable"
    ):
        _build(cal, part, lg, validation_record=lg)  # the pytest log is no signed JSON record
    out = tmp_path / "r.json"
    full = _cli(part, lg, out)
    for opt in ("--host-run-record", "--validation-record", "--git-sha", "--part-record"):
        i = full.index(opt)
        with pytest.raises(SystemExit):  # argparse: every input is required
            cal.main(full[:i] + full[i + 2 :])
    for removed in (
        "--host-run-id",
        "--host-run-token",
        "--validation-id",
    ):  # the token-era options are gone
        with pytest.raises(SystemExit):
            cal.main([*full, removed, "x"])
    with pytest.raises(
        SystemExit, match="not ok: .*unverified"
    ):  # --trust-file: the test keys are not listed there
        cal.main([*full, "--trust-file", str(COMMITTED_TRUST)])
    assert not out.exists()


# -- Codex REVIEW-61f4b0e3 / REVIEW-cf3f019e: the signed attestation chain (V3-INFRA-002) ------------------------------
def test_a_valid_signed_chain_builds_and_loads(
    cal: ModuleType, cov: ModuleType, tmp_path: Path, trust_file: Path
) -> None:
    part, lg = make_artifacts(tmp_path / "art", cov)
    rec = _build(cal, part, lg)
    d = part.parent
    assert rec["schema"] == "ionmc-v7r2-cal-rate-2" and "host_run_token" not in rec
    assert rec["host_run_id"] == RUN_ID and rec["validation_id"] == VAL_ID
    assert rec["host_run_record"] == RUN_FILE and rec["validation_record"] == VAL_FILE
    for key, name in (("host_run", RUN_FILE), ("validation", VAL_FILE)):
        block = json.loads((d / name).read_text())["attestation"]
        assert rec[f"{key}_signature"] == {
            k: block[k] for k in ("key_id", "payload_sha256", "signature")
        }
        assert rec[f"{key}_record_sha256"] == hashlib.sha256((d / name).read_bytes()).hexdigest()
    assert rec["trust_file_sha256"] == hashlib.sha256(trust_file.read_bytes()).hexdigest()
    out = tmp_path / "rate.json"
    out.write_text(json.dumps(rec))
    assert _pinned_load(cal, out)["record"] == rec
    assert _pinned_load(
        cal, _rate_from(tmp_path / "hand", cal, cov)
    )  # the hand-assembled record agrees


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        # the validation-manager record copy
        (
            {"val_tamper": {"status": "passed", "summary": "edited after signing"}},
            "validation record .*invalid attestation: payload hash mismatch",
        ),
        (
            {
                "val_tamper": {
                    "details": {
                        "head_sha": SHA,
                        "host_runs": [{"run_id": RUN_ID, "output_manifest_sha256": "f" * 64}],
                    }
                }
            },
            "validation record .*invalid attestation",
        ),
        ({"val_key": "rogue"}, "validation record .*unverified record rejected"),
        (
            {"val_key": "host-runner"},
            "validation record .*unverified record rejected",
        ),  # a key of the other signer
        (
            {"val_tamper": {"attestation": None}},
            "host_run_signature|validation_signature|unverified",
        ),
        ({"val_over": {"status": "failed"}}, "status 'failed' is not 'passed'"),
        ({"val_over": {"status": "recorded"}}, "is not 'passed'"),
        ({"val_over": {"task_id": "V3-005D"}}, "validation record .*task_id"),
        ({"val_over": {"git_sha": "b" * 40}}, "validation record .*git_sha"),
        (
            {
                "val_det": {
                    "host_runs": [
                        {
                            "run_id": "RUN-20261101T120000Z-ffffffff",
                            "output_manifest_sha256": MANIFEST,
                        }
                    ]
                }
            },
            "run_id .*is not the run 'RUN-20261101T120000Z-ffffffff'",
        ),
        (
            {
                "val_det": {
                    "host_runs": [
                        {"run_id": RUN_ID, "output_manifest_sha256": MANIFEST},
                        {
                            "run_id": "RUN-20261101T120000Z-ffffffff",
                            "output_manifest_sha256": MANIFEST,
                        },
                    ]
                }
            },
            "exactly one host run",
        ),
        ({"val_det": {"host_runs": []}}, "exactly one host run"),
        (
            {"val_det": {"host_runs": [{"run_id": RUN_ID, "output_manifest_sha256": None}]}},
            "is not 64 hex",
        ),
        ({"val_det": {"host_runs": [{"run_id": RUN_ID}]}}, "is not 64 hex"),
        (
            {"val_det": {"host_runs": [{"run_id": RUN_ID, "output_manifest_sha256": "f" * 64}]}},
            "output_manifest_sha256 is 'e+', expected 'f+'",
        ),
        ({"val_det": {"host_run_id": RUN_ID}}, "legacy details.host_run_id"),
        ({"val_det": {"host_run_manifest_lookup": "missing"}}, "host_run_manifest_lookup"),
        ({"val_over": {"details": "none"}}, "details is not an object"),
        # the host-runner record copy
        (
            {"run_tamper": {"exit_code": 0, "stdout": "token x\n"}},
            "host run record .*invalid attestation: payload hash mismatch",
        ),
        ({"run_key": "rogue"}, "host run record .*unverified record rejected"),
        ({"run_key": "validation-manager"}, "host run record .*unverified record rejected"),
        (
            {"run_over": {"output_manifest_sha256": "f" * 64}},
            "output_manifest_sha256 is 'f+', expected 'e+'",
        ),
        ({"run_over": {"output_manifest_sha256": None}}, "output_manifest_sha256 is None"),
        ({"run_over": {"task_id": "V3-005D"}}, "host run record .*task_id"),
        ({"run_over": {"sha": "b" * 40}}, "host run record .*git_sha"),
        ({"run_over": {"succeeded": False}}, "succeeded is False"),
        ({"run_over": {"exit_code": 1, "succeeded": True}}, "exit_code is 1"),
        ({"run_over": {"exit_code": None}}, "exit_code is None"),
        ({"run_over": {"timed_out": True}}, "timed_out is True"),
        ({"run_over": {"committed_state_unchanged": False}}, "committed_state_unchanged is False"),
        ({"run_over": {"dirty_after": True}}, "dirty_after is True"),
        ({"run_drop": ("dirty_after",)}, "dirty_after is None"),
        ({"run_over": {"output_publication_error": "disk full"}}, "output_publication_error"),
        ({"run_drop": ("output_publication_error",)}, "output_publication_error is '<absent>'"),
        ({"run_over": {"stdout_sha256": "0" * 64}}, "does not hash to the signed stdout_sha256"),
        # the binding of the part record, the part log and the sealed token to the signed run
        (
            {"stdout": "kappa parts: 6 passed\n"},
            "token .* appears neither in the signed argv nor as a stdout line",
        ),
        (
            {"stdout": f"token R4B-OTHER-00000002\nCAL-RATE-ARTIFACT part_record x {'0' * 64}\n"},
            "appears neither",
        ),
        ({"stdout": f"token {TOKEN}\n"}, "needs exactly one line 'CAL-RATE-ARTIFACT part_record"),
        (
            {
                "stdout": f"token {TOKEN}\nCAL-RATE-ARTIFACT part_record v7r2-calibration-archive_fitted-p00.json {'0' * 64}\n"
            },
            "CAL-RATE-ARTIFACT part_record",
        ),
        (
            {"over": {"host_run_id": "R4B-OTHER-00000002"}},
            "appears neither",
        ),  # the part ran under another token
    ],
)
def test_the_signed_chain_is_verified_on_build_and_on_every_load(
    kw: dict[str, Any], match: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    with pytest.raises(ValueError, match=match):
        _pinned_load(cal, _rate_from(tmp_path, cal, cov, **kw))
    part, lg = make_artifacts(tmp_path / "b", cov, **kw)  # the builder refuses the same chain
    with pytest.raises(ValueError, match=match):
        _build(cal, part, lg)


def test_the_token_may_be_bound_by_the_signed_argv_and_the_artifact_lines_by_a_truncated_tail(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    part, lg = make_artifacts(
        tmp_path / "a", cov, stdout="", argv=["python", "-c", f"TOKEN = '{TOKEN}'"]
    )
    out = "\n".join(["...tail of a long stream", *_artifact_lines(part, lg)]) + "\n"
    make_artifacts(tmp_path / "a", cov, stdout=out, argv=["python", "-c", f"TOKEN = '{TOKEN}'"],
                   run_over={"stdout_truncated": True, "stdout_total_chars": 10**6, "stdout_sha256": "0" * 64})  # fmt: skip
    rec = _build(
        cal, part, lg
    )  # token in the argv; a truncated inline stdout is not compared with stdout_sha256
    assert rec["host_run_id"] == RUN_ID
    assert cal.token_line(TOKEN) == f"token {TOKEN}"
    assert cal.artifact_lines(part, lg) == _artifact_lines(part, lg)


def test_a_duplicated_artifact_line_is_refused(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    part, lg = make_artifacts(tmp_path / "a", cov, stdout="")
    lines = _artifact_lines(part, lg)
    make_artifacts(
        tmp_path / "a", cov, stdout="\n".join([f"token {TOKEN}", *lines, lines[1]]) + "\n"
    )
    with pytest.raises(ValueError, match="exactly one line 'CAL-RATE-ARTIFACT part_log"):
        _build(cal, part, lg)


@pytest.mark.parametrize("which", [RUN_FILE, VAL_FILE])
@pytest.mark.parametrize("mutate", ["missing", "reformatted", "swapped"])
def test_the_record_copies_are_authenticated_by_their_digests(
    which: str, mutate: str, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    out = _rate_from(tmp_path, cal, cov)
    rec = json.loads(out.read_text())
    key = "host_run_record" if which == RUN_FILE else "validation_record"
    f = Path(rec["log_dir"]) / which
    if mutate == "missing":
        f.unlink()
        match = f"{key} .* is not a file"
    elif mutate == "reformatted":
        f.write_text(json.dumps(json.loads(f.read_text())))  # the same signed content, other bytes
        match = f"{key} .*sha256 differs"
    else:  # another validly signed record of the same chain is not the pinned copy
        signer = "host-runner" if which == RUN_FILE else "validation-manager"
        doc = json.loads(f.read_text())
        doc.pop("attestation")
        f.write_text(json.dumps(_sign(doc, signer), indent=2) + "\n")
        match = f"{key} .*sha256 differs"
    with pytest.raises(ValueError, match=match):
        _pinned_load(cal, out)


def test_a_re_signed_copy_with_a_fresh_digest_still_needs_the_recorded_signature(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    """Re-signing a copy (same payload, a new signature) and updating its file digest in the record leaves the recorded
    signature stale: refused."""
    out = _rate_from(tmp_path, cal, cov)
    rec = json.loads(out.read_text())
    f = Path(rec["log_dir"]) / VAL_FILE
    doc = json.loads(f.read_text())
    doc.pop("attestation")
    doc = {**doc, "summary": "re-signed"}
    f.write_text(json.dumps(_sign(doc, "validation-manager"), indent=2) + "\n")
    rec["validation_record_sha256"] = hashlib.sha256(f.read_bytes()).hexdigest()
    out.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="validation_signature of the record differs"):
        _pinned_load(cal, out)


def test_a_changed_trust_file_invalidates_the_record(
    cal: ModuleType, cov: ModuleType, tmp_path: Path, trust_file: Path
) -> None:
    out = _rate_from(tmp_path, cal, cov)
    plan = _plan_for(cal, out)
    other = tmp_path / "trust.json"
    other.write_text(trust_file.read_text() + "\n")  # same keys, other bytes
    with pytest.raises(
        ValueError, match="trust_file_sha256 of the record differs from the trust file"
    ):
        cal.load_record(out, plan_path=plan, trust_path=other)
    doc = json.loads(trust_file.read_text())
    doc["host-runner"] = []  # the host-runner key is no longer trusted
    other.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="host run record .*unverified record rejected"):
        cal.load_record(out, plan_path=plan, trust_path=other)
    with pytest.raises(ValueError, match="trust file .* unreadable"):
        cal.load_record(out, plan_path=plan, trust_path=tmp_path / "absent.json")


def test_a_legacy_token_attested_record_is_refused(
    cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    """The schema-1 record of 2026-11-07 (``host_run_token``, repository-authored validation copy) is refused with an
    explicit message, by the loader and by every consumer path."""
    legacy = {"schema": "ionmc-v7r2-cal-rate-1", "s_per_rep": 9.7, "fixed_s": 70.0, "reps": 83, "wall_mc_s": 805.1,
              "wall_total_s": 814.32, "host": {"node": "n", "platform": "p", "machine": "m", "cpu_count": 1, "python": "3"},
              "host_run_id": RUN_ID, "host_run_token": TOKEN, "validation_id": VAL_ID, "validation_record": VAL_FILE,
              "validation_record_sha256": "0" * 64, "git_sha": SHA, "log_dir": str(tmp_path), "part_record": "p.json",
              "part_record_sha256": "0" * 64, "part_log": "part.log", "part_log_sha256": "0" * 64}  # fmt: skip
    p = tmp_path / "legacy.json"
    p.write_text(json.dumps(legacy))
    with pytest.raises(ValueError, match="legacy token-attested calibration rate record .*refused"):
        _pinned_load(cal, p)
    legacy["schema"] = (
        "ionmc-v7r2-cal-rate-2"  # a token field under the new schema is refused the same way
    )
    p.write_text(json.dumps(legacy))
    with pytest.raises(ValueError, match="legacy token-attested"):
        _pinned_load(cal, p)


def test_a_record_with_the_untrusted_key_on_both_copies_is_refused_even_when_pinned(
    s5d: ModuleType, cal: ModuleType, cov: ModuleType, tmp_path: Path
) -> None:
    """A self-consistent chain signed by a key outside the trust file: the pin (which a forger could write) does not
    rescue it, and no consumer accepts it."""
    out = _rate_from(tmp_path, cal, cov, run_key="rogue", val_key="rogue")
    plan = _plan_for(cal, out)
    with pytest.raises(ValueError, match="unverified record rejected"):
        cal.load_record(out, plan_path=plan)
    rep = s5d.rep_plan_table("lv5d", str(out), str(plan))
    assert rep["ok"] is False and "unverified" in rep["reason"] and rep["plans"] == {}


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


# -- --cal-rate-record: the host runner passes argv only, so the option carries the rate record -----------------
def _main_harness(
    rs: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> dict[str, Any]:
    """Run ``rs.main`` with one fake step; returns what the preflight and the step subprocess saw."""
    seen: dict[str, Any] = {"preflight_env": [], "step_envs": []}
    monkeypatch.setenv(
        RATE_ENV, ""
    )  # registers the restore of the in-process variable the option sets
    monkeypatch.delenv(RATE_ENV)
    monkeypatch.setattr(rs, "ALLOWED_PARENTS", (tmp_path,))
    monkeypatch.setattr(rs, "suite_steps", lambda *a, **k: [("10-fake", ["true"], {})])
    monkeypatch.setattr(rs, "rehearsal_preflight",
                        lambda *a, **k: seen["preflight_env"].append(os.environ.get(RATE_ENV)))  # fmt: skip
    monkeypatch.setattr(rs, "resolve_sha", lambda e: ("a" * 40, "test"))
    monkeypatch.setattr(rs, "git", lambda *a: "")
    monkeypatch.setattr(rs, "source_files", lambda suite: [])
    monkeypatch.setattr(rs, "echo_partial_lines", lambda p: None)

    def fake_run_step(cmd: Any, **kw: Any) -> int:
        seen["step_envs"].append(kw["env"])
        return 0

    monkeypatch.setattr(rs, "run_step", fake_run_step)
    # the archive of the fake step is no archive: only the environment is under test
    monkeypatch.setitem(sys.modules, "summarize", types.SimpleNamespace(main=lambda argv: 0))
    return seen


def _suite_argv(out: Path, *extra: str) -> list[str]:
    return ["--suite", "lv5d", "--out", str(out), "--expected-sha", "a" * 40, "--single-process", "--workers", "1",
            *extra]  # fmt: skip


def test_cal_rate_record_option_sets_the_environment_of_preflight_and_steps(
    rs: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen = _main_harness(rs, tmp_path, monkeypatch)
    rec = tmp_path / "rate.json"
    rec.write_text('{"x": 1}')
    out = tmp_path / "run"
    assert rs.main(_suite_argv(out, "--cal-rate-record", str(rec))) == 0
    assert seen["preflight_env"] == [str(rec)]  # the in-process preflight sees it
    assert [e[RATE_ENV] for e in seen["step_envs"]] == [str(rec)]  # and every step subprocess
    env_txt = (out / "environment.txt").read_text()
    assert f"cal_rate_record={rec}\n" in env_txt
    assert f"cal_rate_record_sha256={hashlib.sha256(rec.read_bytes()).hexdigest()}\n" in env_txt
    # a repo-relative path is resolved against the repository
    seen2 = _main_harness(rs, tmp_path, monkeypatch)
    rel = "validation/results/transport/v7r2-r4b/cal-rate.json"
    assert rs.main(_suite_argv(tmp_path / "run2", "--cal-rate-record", rel)) == 0
    assert seen2["preflight_env"] == [str(REPO / rel)] and seen2["step_envs"][0][RATE_ENV] == str(
        REPO / rel
    )


def test_cal_rate_record_option_missing_file_fails_closed(
    rs: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen = _main_harness(rs, tmp_path, monkeypatch)
    out = tmp_path / "run"
    with pytest.raises(SystemExit, match="not ok: no such file"):
        rs.main(_suite_argv(out, "--cal-rate-record", str(tmp_path / "missing.json")))
    assert not out.exists() and seen["step_envs"] == []
    with pytest.raises(SystemExit, match="lv5d suite only"):
        rs.main(
            [
                "--suite",
                "lv5c",
                "--out",
                str(out),
                "--expected-sha",
                "a" * 40,
                "--cal-rate-record",
                "x",
            ]
        )


def test_without_the_option_the_environment_and_environment_txt_are_unchanged(
    rs: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen = _main_harness(rs, tmp_path, monkeypatch)
    out = tmp_path / "run"
    assert rs.main(_suite_argv(out)) == 0
    assert seen["preflight_env"] == [None] and RATE_ENV not in seen["step_envs"][0]
    assert "cal_rate_record" not in (out / "environment.txt").read_text()
