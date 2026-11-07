"""Gate tests: local_validation check()/inspect() enforce attestations."""

import copy
import hashlib
import importlib.util
import json
import resource
import sys
from pathlib import Path

import pytest

from infrastructure.experiment_v3 import attest
from infrastructure.validation import local_validation as lv

ROOT = Path(__file__).resolve().parents[2]
TASK = "V3-X1"
SHA = "a" * 40


def _load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BOOT = "V3-INFRA-002"


def _rec(status="passed", task=TASK, **over):
    rec = {
        "schema_version": 1,
        "validation_id": "VAL-1",
        "task_id": task,
        "suite": "s",
        "status": status,
        "summary": "ok",
        "details": {},
        "git_sha": SHA,
    }
    rec.update(over)
    return rec


@pytest.fixture
def env(tmp_path, monkeypatch):
    key = attest.generate_keypair()
    monkeypatch.setattr(lv, "VALIDATION_ROOT", tmp_path / "val")
    monkeypatch.setattr(lv, "TRUSTED_KEYS_PATH", tmp_path / "trust.json")
    monkeypatch.setattr(
        lv,
        "worktree_state",
        lambda t: {"path": tmp_path, "branch": "b", "sha": SHA, "dirty": False},
    )

    def set_trust(*keys):
        entries = [{"key_id": k.key_id, "public_key": k.public_key.hex()} for k in keys]
        (tmp_path / "trust.json").write_text(
            json.dumps({"validation-manager": entries, "host-runner": []})
        )

    def put(rec, name="r1", task=TASK):
        d = lv.VALIDATION_ROOT / task / SHA
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{name}.json").write_text(json.dumps(rec))

    set_trust()
    return key, set_trust, put


def _check(capsys, task=TASK):
    try:
        lv.check(task)
        code = 0
    except SystemExit as exc:
        code = exc.code
    return code, json.loads(capsys.readouterr().out)


def _sign(key, rec=None, signer="validation-manager"):
    return attest.sign_record(rec or _rec(), key, signer)


def test_enforcing_verified_passes(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put(_sign(key))
    code, out = _check(capsys)
    assert code == 0 and out["status"] == "passed"
    assert out["enforcement"] == "enforcing" and out["verified_records"] == 1


def test_enforcing_tampered_rejected_as_invalid(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    signed = _sign(key)
    signed["status"] = "passed"
    signed["summary"] = "edited after signing"
    put(signed)
    code, out = _check(capsys)
    assert code == 1
    assert out["rejected_records"][0]["attestation_status"] == "invalid"
    assert out["passing_records"] == 0


def test_enforcing_untrusted_key_rejected(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put(_sign(attest.generate_keypair()))
    code, out = _check(capsys)
    assert code == 1
    assert out["rejected_records"][0]["attestation_status"] == "unverified"
    assert "not trusted" in out["rejected_records"][0]["reason"]


def test_enforcing_unsigned_rejected(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put(_rec())
    code, out = _check(capsys)
    assert code == 1
    assert out["rejected_records"][0]["reason"] == "record has no attestation"


@pytest.mark.parametrize("field,value", [("task_id", "V3-OTHER"), ("git_sha", "b" * 40)])
def test_wrong_signed_identity_rejected(env, capsys, field, value, tmp_path):
    key, set_trust, put = env
    set_trust(key)
    put(_sign(key, _rec(**{field: value})))  # validly signed, wrong identity
    code, out = _check(capsys)
    assert code == 1
    assert out["rejected_records"][0]["attestation_status"] == "invalid"
    assert "identity mismatch" in out["rejected_records"][0]["reason"]


def test_wrong_signer_rejected(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put(_sign(key, signer="host-runner"))
    code, out = _check(capsys)
    assert code == 1
    assert out["rejected_records"][0]["attestation_status"] == "invalid"


def test_enforcing_unverified_passed_plus_verified_failed_still_fails(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put(_rec(validation_id="VAL-2"), "r2")
    put(_sign(key, _rec("failed")), "r1")
    code, out = _check(capsys)
    assert code == 1 and out["failing_records"] == 1 and out["passing_records"] == 0


def test_forged_derived_fields_on_disk_ignored(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put({**_rec(), "attestation_status": "verified"})
    code, out = _check(capsys)
    assert code == 1


def test_bootstrap_tasks_constant():
    assert lv.BOOTSTRAP_TASKS == frozenset({"V3-INFRA-002"})


def test_bootstrap_tolerates_unsigned_for_migration_task(env, capsys):
    key, set_trust, put = env
    put(_rec(task=BOOT), task=BOOT)
    code, out = _check(capsys, BOOT)
    assert code == 0 and out["enforcement"] == "bootstrap"
    assert out["unverified_records"] == 1
    assert "unverified (no trusted key published yet)" in out["attestation_note"]


def test_bootstrap_signed_but_untrusted_tolerated_for_migration_task(env, capsys):
    key, set_trust, put = env
    put(_sign(key, _rec(task=BOOT)), task=BOOT)
    code, out = _check(capsys, BOOT)
    assert code == 0 and "attestation_note" in out


@pytest.mark.parametrize("task", [TASK, "V3-005B"])
def test_other_tasks_unverified_rejected_with_empty_trust(env, capsys, task):
    key, set_trust, put = env
    put(_rec(task=task), task=task)
    code, out = _check(capsys, task)
    assert code == 1 and out["enforcement"] == "blocked"
    assert out["passing_records"] == 0
    assert out["rejected_records"][0]["reason"] == lv.ATTESTATION_REQUIRED
    assert out["attestation_required"] == lv.ATTESTATION_REQUIRED
    assert "attestation_note" not in out


def test_other_task_signed_untrusted_rejected_with_empty_trust(env, capsys):
    key, set_trust, put = env
    put(_sign(key))
    code, out = _check(capsys)
    assert code == 1 and out["passing_records"] == 0


def test_other_task_verified_passes_once_key_published(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put(_sign(key))
    code, out = _check(capsys)
    assert code == 0 and out["enforcement"] == "enforcing"


def test_migration_task_enforcing_once_key_published(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put(_rec(task=BOOT), task=BOOT)
    code, out = _check(capsys, BOOT)
    assert code == 1 and out["enforcement"] == "enforcing"


@pytest.mark.parametrize("task", [BOOT, TASK])
def test_bootstrap_still_rejects_invalid(env, capsys, task):
    key, set_trust, put = env
    signed = _sign(key, _rec(task=task))
    signed["summary"] = "tampered"
    put(signed, task=task)
    code, out = _check(capsys, task)
    assert code == 1
    assert out["rejected_records"][0]["attestation_status"] == "invalid"


def test_inspect_reports_attestation_status(env, capsys):
    key, set_trust, put = env
    set_trust(key)
    put(_sign(key))
    lv.inspect(TASK)
    out = json.loads(capsys.readouterr().out)
    assert out["records"][0]["attestation_status"] == "verified"
    assert out["gate_passed"] is True


def test_committed_trust_map_is_bootstrap_or_loadable():
    trusted = attest.load_trusted(lv.TRUSTED_KEYS_PATH)
    assert set(trusted) == {"validation-manager", "host-runner"}


def test_require_verified_helper():
    av = _load("attest_verify_gate", "validation/scripts/transport/attest_verify.py")
    key = attest.generate_keypair()
    signed = attest.sign_record(_rec(), key, "host-runner")
    trusted = {"host-runner": {key.key_id: key.public_key.hex()}}
    assert av.require_verified(signed, "host-runner", trusted)["status"] == "verified"
    with pytest.raises(av.AttestationRequired):
        av.require_verified({**signed, "summary": "x"}, "host-runner", trusted)
    with pytest.raises(av.AttestationRequired):
        av.require_verified(_rec(), "host-runner", trusted)
    other = {"host-runner": {"k": "00"}}
    with pytest.raises(av.AttestationRequired):
        av.require_verified(signed, "host-runner", other)
    empty = {"host-runner": {}}
    # strict by default: unverified rejected even with an empty trust map
    with pytest.raises(av.AttestationRequired):
        av.require_verified(_rec(), "host-runner", empty)
    # explicit bootstrap flag tolerates unverified only while no key is trusted
    assert (
        av.require_verified(_rec(), "host-runner", empty, bootstrap=True)["status"]
        == "unverified"
    )
    with pytest.raises(av.AttestationRequired):
        av.require_verified(_rec(), "host-runner", other, bootstrap=True)
    with pytest.raises(av.AttestationRequired):
        av.require_verified({**signed, "summary": "x"}, "host-runner", empty, bootstrap=True)


def test_require_verified_identity_binding():
    av = _load("attest_verify_bind", "validation/scripts/transport/attest_verify.py")
    key = attest.generate_keypair()
    run = {"run_id": "RUN-1", "task_id": TASK, "sha": SHA, "output_manifest_sha256": "m" * 64}
    signed = attest.sign_record(run, key, "host-runner")
    trusted = {"host-runner": {key.key_id: key.public_key.hex()}}
    ok = av.require_verified(
        signed, "host-runner", trusted, expected_task_id=TASK, expected_git_sha=SHA,
        expected_run_id="RUN-1", expected_output_manifest_sha256="m" * 64,
    )
    assert ok["status"] == "verified"
    for kw in (
        {"expected_task_id": "V3-OTHER"},
        {"expected_git_sha": "b" * 40},
        {"expected_run_id": "RUN-2"},
        {"expected_output_manifest_sha256": "n" * 64},
    ):
        with pytest.raises(av.AttestationRequired, match="identity mismatch"):
            av.require_verified(signed, "host-runner", trusted, **kw)


def test_manifest_hash_streams_large_sparse_file(tmp_path):
    hr = _load("ionmc_host_runner_sparse", "infrastructure/host_runner/host_runner.py")
    outputs = tmp_path / "outputs"
    outputs.mkdir()
    big = outputs / "big.bin"
    size = 3 * 1024**3
    try:
        with open(big, "wb") as stream:
            stream.truncate(size)
    except OSError:
        size = 256 * 1024**2  # filesystem without large sparse files
        with open(big, "wb") as stream:
            stream.truncate(size)
    reference = hashlib.sha256()
    with open(big, "rb") as stream:
        while chunk := stream.read(1 << 20):
            reference.update(chunk)
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # KiB on Linux
    assert hr.file_sha256(big) == reference.hexdigest()
    expected = hashlib.sha256(
        f"big.bin\0{reference.hexdigest()}\n".encode()
    ).hexdigest()
    assert hr.output_manifest_sha256(outputs) == expected
    grown_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss - before
    assert grown_kib < 256 * 1024  # far below the file size: streamed, not read_bytes
