"""Host-runner attestation consumer and its wiring into the local-validation gate."""

import json

import pytest

from infrastructure.experiment_v3 import attest
from infrastructure.validation import host_run_records as hrr
from infrastructure.validation import local_validation as lv

BOOT = "V3-INFRA-002"
OTHER = "V3-X1"
SHA = "a" * 40
RUN = "RUN-20260101-000000-ABC123"
MANIFEST = "m" * 64


def _run(task=OTHER, sha=SHA, run_id=RUN, exit_code=0, timed_out=False, **over):
    rec = {
        "schema_version": 1, "run_id": run_id, "task_id": task, "sha": sha,
        "exit_code": exit_code, "timed_out": timed_out,
        "output_manifest_sha256": MANIFEST,
    }
    rec.update(over)
    return rec


def _val(task=OTHER, run_ids=(RUN,), **over):
    rec = {
        "schema_version": 1, "validation_id": "VAL-1", "task_id": task, "suite": "s",
        "status": "passed", "summary": "ok", "git_sha": SHA,
        "details": {"host_run_id": run_ids[0]} if len(run_ids) == 1 else {"host_run_ids": list(run_ids)},
    }
    rec.update(over)
    return rec


@pytest.fixture
def env(tmp_path, monkeypatch):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    monkeypatch.setattr(lv, "VALIDATION_ROOT", tmp_path / "val")
    monkeypatch.setattr(lv, "TRUSTED_KEYS_PATH", tmp_path / "trust.json")
    monkeypatch.setattr(hrr, "HOST_RUNS_ROOT", tmp_path / "runs")
    monkeypatch.setattr(
        lv, "worktree_state",
        lambda t: {"path": tmp_path, "branch": "b", "sha": SHA, "dirty": False},
    )

    def trust(vm_keys=(), hr_keys=()):
        def ent(k):
            return {"key_id": k.key_id, "public_key": k.public_key.hex()}
        (tmp_path / "trust.json").write_text(json.dumps({
            "validation-manager": [ent(k) for k in vm_keys],
            "host-runner": [ent(k) for k in hr_keys],
        }))

    def put_val(rec, task=OTHER):
        d = lv.VALIDATION_ROOT / task / SHA
        d.mkdir(parents=True, exist_ok=True)
        (d / "v.json").write_text(json.dumps(rec))

    def put_run(rec, task=OTHER, run_id=RUN):
        d = tmp_path / "runs" / task / run_id
        d.mkdir(parents=True, exist_ok=True)
        (d / "record.json").write_text(json.dumps(rec))

    trust()
    return vm, hr, trust, put_val, put_run


def _check(capsys, task):
    try:
        lv.check(task)
        code = 0
    except SystemExit as exc:
        code = exc.code
    return code, json.loads(capsys.readouterr().out)


def _trust_map(hr):
    return {"host-runner": {hr.key_id: hr.public_key.hex()}}


def test_verify_run_record_verified(env):
    _, hr, *_ = env
    signed = attest.sign_record(_run(), hr, "host-runner")
    assert hrr.verify_run_record(signed, OTHER, SHA, RUN, _trust_map(hr), MANIFEST) == (
        "verified", "ok",
    )


@pytest.mark.parametrize("kw", [
    {"expected_task_id": "V3-OTHER"}, {"expected_sha": "b" * 40},
    {"expected_run_id": "RUN-OTHER"}, {"expected_output_manifest_sha256": "n" * 64},
])
def test_verify_run_record_replay_invalid(env, kw):
    _, hr, *_ = env
    signed = attest.sign_record(_run(), hr, "host-runner")
    args = {"expected_task_id": OTHER, "expected_sha": SHA, "expected_run_id": RUN}
    args.update(kw)
    status, reason = hrr.verify_run_record(signed, trust=_trust_map(hr), **args)
    assert status == "invalid" and "identity mismatch" in reason


def test_verify_run_record_tampered_manifest_invalid(env):
    _, hr, *_ = env
    signed = attest.sign_record(_run(), hr, "host-runner")
    signed["output_manifest_sha256"] = "0" * 64
    status, _ = hrr.verify_run_record(signed, OTHER, SHA, RUN, _trust_map(hr), "0" * 64)
    assert status == "invalid"


def test_verify_run_record_unsigned_unverified(env):
    _, hr, *_ = env
    assert hrr.verify_run_record(_run(), OTHER, SHA, RUN, _trust_map(hr))[0] == "unverified"


def test_load_signed_run_record(env):
    *_, put_run = env
    put_run(_run())
    assert hrr.load_signed_run_record(RUN)["run_id"] == RUN
    assert hrr.load_signed_run_record(RUN, OTHER)["run_id"] == RUN
    with pytest.raises(hrr.RunRecordMissing):
        hrr.load_signed_run_record(RUN, "V3-OTHER")
    with pytest.raises(hrr.RunRecordMissing):
        hrr.load_signed_run_record("../etc")


def _enforcing_setup(env):
    vm, hr, trust, put_val, put_run = env
    trust([vm], [hr])
    return vm, hr, put_val, put_run


def test_gate_verified_run_passes_and_inspect_reports(env, capsys):
    vm, hr, put_val, put_run = _enforcing_setup(env)
    put_run(attest.sign_record(_run(), hr, "host-runner"))
    put_val(attest.sign_record(_val(), vm, "validation-manager"))
    lv.inspect(OTHER)
    out = json.loads(capsys.readouterr().out)
    att = out["records"][0]["host_run_attestation"]
    assert att == [{"run_id": RUN, "status": "verified", "reason": "ok", "ok": True}]
    assert out["gate_passed"] is True


@pytest.mark.parametrize("field,value", [("task_id", "V3-OTHER"), ("sha", "b" * 40)])
def test_gate_replayed_run_record_rejected(env, capsys, field, value):
    vm, hr, put_val, put_run = _enforcing_setup(env)
    put_run(attest.sign_record(_run(**{field: value}), hr, "host-runner"))
    put_val(attest.sign_record(_val(), vm, "validation-manager"))
    code, out = _check(capsys, OTHER)
    assert code == 1 and "identity mismatch" in out["rejected_records"][0]["reason"]


def test_gate_tampered_run_record_rejected(env, capsys):
    vm, hr, put_val, put_run = _enforcing_setup(env)
    signed = attest.sign_record(_run(), hr, "host-runner")
    signed["output_manifest_sha256"] = "0" * 64
    put_run(signed)
    put_val(attest.sign_record(_val(), vm, "validation-manager"))
    code, out = _check(capsys, OTHER)
    assert code == 1 and "payload hash mismatch" in out["rejected_records"][0]["reason"]


def test_gate_unsigned_run_rejected_when_enforcing(env, capsys):
    vm, hr, put_val, put_run = _enforcing_setup(env)
    put_run(_run())
    put_val(attest.sign_record(_val(), vm, "validation-manager"))
    code, out = _check(capsys, OTHER)
    assert code == 1 and "attestation required" in out["rejected_records"][0]["reason"]


def test_gate_missing_run_record_rejected(env, capsys):
    vm, hr, put_val, put_run = _enforcing_setup(env)
    put_val(attest.sign_record(_val(), vm, "validation-manager"))
    code, out = _check(capsys, OTHER)
    assert code == 1 and "missing" in out["rejected_records"][0]["reason"]


@pytest.mark.parametrize("over", [{"exit_code": 1}, {"exit_code": None, "timed_out": True}])
def test_gate_failed_run_rejected(env, capsys, over):
    vm, hr, put_val, put_run = _enforcing_setup(env)
    put_run(attest.sign_record(_run(**over), hr, "host-runner"))
    put_val(attest.sign_record(_val(), vm, "validation-manager"))
    code, out = _check(capsys, OTHER)
    assert code == 1 and "did not succeed" in out["rejected_records"][0]["reason"]


def test_gate_bootstrap_unsigned_run_tolerated_for_migration_task_only(env, capsys):
    vm, hr, trust, put_val, put_run = env  # empty trust map
    put_run(_run(task=BOOT), task=BOOT)
    put_val(_val(task=BOOT), task=BOOT)
    code, out = _check(capsys, BOOT)
    assert code == 0 and out["enforcement"] == "bootstrap"
    put_run(_run())
    put_val(_val())
    code, out = _check(capsys, OTHER)
    assert code == 1 and out["passing_records"] == 0


def test_gate_bootstrap_missing_run_still_rejected(env, capsys):
    vm, hr, trust, put_val, put_run = env
    put_val(_val(task=BOOT), task=BOOT)
    code, out = _check(capsys, BOOT)
    assert code == 1 and "missing" in out["rejected_records"][0]["reason"]
