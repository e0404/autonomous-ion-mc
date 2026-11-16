"""Signed codex-worker review records (V3-INFRA-006).

execute signs the finished review record with the service-held codex-worker
key and stores it as record.json next to report.json; inspect returns it and
verifies it on read; infrastructure/validation/review_records.py classifies it
verified / unverified / invalid.
"""

import json
import subprocess
from pathlib import Path

import pytest

from infrastructure.experiment_v3 import attest, host_snapshot, review
from infrastructure.experiment_v3.common import file_hash, write_json
from infrastructure.experiment_v3.paths import EXPERIMENT_ID
from infrastructure.validation import review_records as rvr

SHA = "a" * 40
BASE = "b" * 40
RID = "REVIEW-" + "c" * 32
TASK = "V3-TEST"
README = Path(__file__).resolve().parents[3] / "infrastructure/keys/README.md"


@pytest.fixture
def keys(tmp_path, monkeypatch):
    monkeypatch.setattr(attest, "key_dir", lambda: tmp_path / "keys")
    return tmp_path / "keys"


def trust_of(key):
    return {"codex-worker": {key.key_id: key.public_key.hex()}}


def run_review(tmp_path, monkeypatch, *, verdict="pass", findings=(), exit_code=0):
    """execute() with the Codex process and the snapshot replaced."""
    state = tmp_path / "state"
    dest = state / "reviews" / RID
    write_json(
        dest / "job.json",
        {
            "review_id": RID,
            "task_id": TASK,
            "head_sha": SHA,
            "base_sha": BASE,
            "worktree": str(tmp_path / "wt"),
            "status": "running",
            "started_at": "2026-09-30T00:00:00+00:00",
            "model": "gpt-5.6-sol",
            "provider": "openai",
            "tool": "codex exec",
            "experiment_id": EXPERIMENT_ID,
        },
    )
    monkeypatch.setattr(review, "exact_state", lambda _: SHA)
    monkeypatch.setattr(host_snapshot, "create", lambda *a, **k: tmp_path)

    def codex(cmd, **kwargs):
        assert cmd[:2] == ["codex", "exec"]
        write_json(
            dest / "report.json",
            {
                "head_sha": SHA,
                "verdict": verdict,
                "summary": "Inspected the diff.",
                "findings": list(findings),
            },
        )
        return subprocess.CompletedProcess(cmd, exit_code)

    monkeypatch.setattr(review.subprocess, "run", codex)
    job = review.execute(RID, state=state)
    return state, dest, job


def test_execute_signs_record_next_to_report(tmp_path, monkeypatch, keys):
    state, dest, job = run_review(tmp_path, monkeypatch)
    assert job["status"] == "passed" and "record_error" not in job
    record = json.loads((dest / "record.json").read_text())
    key = attest.load_or_create_key("codex-worker")
    assert record["attestation"]["signer"] == "codex-worker"
    assert record["attestation"]["key_id"] == key.key_id
    attest.verify_record(record, trust_of(key), "codex-worker")
    expected = {k: v for k, v in job.items() if k != "worktree"}
    expected["schema_version"] = 1
    expected["report"] = json.loads((dest / "report.json").read_text())
    assert {k: v for k, v in record.items() if k != "attestation"} == expected
    assert record["report_sha256"] == file_hash(dest / "report.json")
    assert record["started_at"] < record["finished_at"]
    assert review.inspect(RID, state=state) == record
    assert rvr.verify_review_record(
        record,
        trust_of(key),
        expected_task_id=TASK,
        expected_head_sha=SHA,
        expected_base_sha=BASE,
        report_path=dest / "report.json",
    ) == ("verified", "ok")
    assert review.require_review(TASK, tmp_path, BASE, state=state) == RID


def test_changes_required_and_failed_reviews_are_signed(tmp_path, monkeypatch, keys):
    finding = {"severity": "important", "path": "a.py", "description": "fix"}
    state, dest, job = run_review(
        tmp_path, monkeypatch, verdict="changes_required", findings=[finding]
    )
    trust = trust_of(attest.load_or_create_key("codex-worker"))
    record = review.inspect(RID, state=state)
    assert record["status"] == "changes_required" and record["attestation"]
    assert rvr.verify_review_record(record, trust)[0] == "verified"
    state, dest, job = run_review(tmp_path, monkeypatch, exit_code=1)
    record = review.inspect(RID, state=state)
    assert record["status"] == "failed" and "report" not in record
    assert "Codex exited 1" in record["error"]
    assert rvr.verify_review_record(record, trust) == ("verified", "ok")


@pytest.mark.parametrize("target", ["record", "job", "report"])
def test_inspect_rejects_tampering(tmp_path, monkeypatch, keys, target):
    state, dest, _ = run_review(tmp_path, monkeypatch)
    if target == "record":
        record = json.loads((dest / "record.json").read_text())
        record["report"]["summary"] = "Tampered"
        write_json(dest / "record.json", record)
        match = "attestation invalid"
    elif target == "job":
        job = json.loads((dest / "job.json").read_text())
        job["base_sha"] = "d" * 40
        write_json(dest / "job.json", job)
        match = "does not match"
    else:
        report = json.loads((dest / "report.json").read_text())
        report["summary"] = "Tampered"
        write_json(dest / "report.json", report)
        match = "hash"
    with pytest.raises(ValueError, match=match):
        review.inspect(RID, state=state)
    with pytest.raises(ValueError):
        review.require_review(TASK, tmp_path, BASE, state=state)


def test_inspect_rejects_unsigned_or_wrong_signer_record(tmp_path, monkeypatch, keys):
    state, dest, _ = run_review(tmp_path, monkeypatch)
    record = json.loads((dest / "record.json").read_text())
    payload = {k: v for k, v in record.items() if k != "attestation"}
    write_json(dest / "record.json", payload)
    with pytest.raises(ValueError, match="attestation invalid"):
        review.inspect(RID, state=state)
    other = attest.sign_record(payload, attest.generate_keypair(), "reference-runner")
    write_json(dest / "record.json", other)
    with pytest.raises(ValueError, match="signer"):
        review.inspect(RID, state=state)


def test_review_finished_before_signing_still_inspects(tmp_path, monkeypatch):
    monkeypatch.setattr(review, "exact_state", lambda _: SHA)
    dest = tmp_path / "reviews" / RID
    write_json(
        dest / "report.json",
        {"head_sha": SHA, "verdict": "pass", "summary": "ok", "findings": []},
    )
    job = {
        "review_id": RID,
        "task_id": TASK,
        "head_sha": SHA,
        "base_sha": BASE,
        "worktree": "/x",
        "status": "passed",
        "started_at": "2026-09-30T00:00:00+00:00",
        "finished_at": "2026-09-30T00:10:00+00:00",
        "report_sha256": file_hash(dest / "report.json"),
    }
    write_json(dest / "job.json", job)
    result = review.inspect(RID, state=tmp_path)
    assert result["attestation"] is None and result["status"] == "passed"
    assert "worktree" not in result and result["report"]["verdict"] == "pass"
    assert review.require_review(TASK, tmp_path, BASE, state=tmp_path) == RID


def test_signing_failure_keeps_review_outcome(tmp_path, monkeypatch, keys):
    def broken(*_):
        raise attest.AttestationError("key unavailable")

    monkeypatch.setattr(attest, "load_or_create_key", broken)
    state, dest, job = run_review(tmp_path, monkeypatch)
    assert job["status"] == "passed" and "key unavailable" in job["record_error"]
    assert not (dest / "record.json").exists()
    assert review.inspect(RID, state=state)["attestation"] is None


def signed(record, key=None, signer="codex-worker"):
    return attest.sign_record(record, key or attest.generate_keypair(), signer)


def base_record():
    return {
        "schema_version": 1,
        "review_id": RID,
        "task_id": TASK,
        "head_sha": SHA,
        "base_sha": BASE,
        "status": "passed",
        "started_at": "2026-09-30T00:00:00+00:00",
        "finished_at": "2026-09-30T00:10:00+00:00",
        "model": "gpt-5.6-sol",
        "provider": "openai",
        "tool": "codex exec",
        "experiment_id": EXPERIMENT_ID,
        "report_sha256": "e" * 64,
        "report": {"head_sha": SHA, "verdict": "pass", "summary": "ok", "findings": []},
    }


def test_verifier_statuses():
    key = attest.generate_keypair()
    trust = trust_of(key)
    record = signed(base_record(), key)
    assert rvr.verify_review_record(record, trust) == ("verified", "ok")
    assert rvr.verify_review_record(record, {})[0] == "unverified"
    assert rvr.verify_review_record(base_record(), trust)[0] == "unverified"
    assert rvr.verify_review_record("x", trust)[0] == "invalid"
    tampered = {**record, "status": "changes_required"}
    assert rvr.verify_review_record(tampered, trust)[0] == "invalid"
    wrong = signed(base_record(), key, "reference-runner")
    assert rvr.verify_review_record(wrong, trust) == (
        "invalid",
        "signer 'reference-runner' != expected 'codex-worker'",
    )


@pytest.mark.parametrize(
    "change, reason",
    [
        ({"schema_version": 2}, "schema_version"),
        ({"review_id": "REVIEW-x"}, "review_id"),
        ({"head_sha": "short"}, "head_sha"),
        ({"model": ""}, "model"),
        ({"experiment_id": ""}, "experiment_id"),
        ({"experiment_id": None}, "experiment_id"),
        ({"status": "running"}, "status"),
        ({"started_at": "2026-09-30"}, "started_at"),
        ({"finished_at": "2026-09-30T00:00:00+00:00"}, "not later than"),
        ({"finished_at": "2026-09-29T00:00:00+00:00"}, "not later than"),
        ({"report_sha256": None}, "report_sha256"),
        ({"status": "changes_required"}, "inconsistent"),
        (
            {
                "report": {
                    "head_sha": BASE,
                    "verdict": "pass",
                    "summary": "s",
                    "findings": [],
                }
            },
            "report invalid",
        ),
        ({"status": "failed"}, "failed review"),
    ],
)
def test_verifier_rejects_inconsistent_payload(change, reason):
    key = attest.generate_keypair()
    status, why = rvr.verify_review_record(
        signed({**base_record(), **change}, key), trust_of(key)
    )
    assert status == "invalid" and reason in why


def test_verifier_rejects_identity_and_report_mismatch(tmp_path):
    key = attest.generate_keypair()
    record = signed(base_record(), key)
    trust = trust_of(key)
    for kwargs in (
        {"expected_task_id": "V3-OTHER"},
        {"expected_head_sha": BASE},
        {"expected_base_sha": SHA},
    ):
        status, why = rvr.verify_review_record(record, trust, **kwargs)
        assert status == "invalid" and "signed identity mismatch" in why
    (tmp_path / "report.json").write_text("{}")
    status, why = rvr.verify_review_record(
        record, trust, report_path=tmp_path / "report.json"
    )
    assert status == "invalid" and "report_sha256" in why


def test_cli(tmp_path, monkeypatch, keys, capsys):
    state, dest, _ = run_review(tmp_path, monkeypatch)
    key = attest.load_or_create_key("codex-worker")
    trust_path = tmp_path / "trust.json"
    trust_path.write_text(
        json.dumps(
            {
                "codex-worker": [
                    {"key_id": key.key_id, "public_key": key.public_key.hex()}
                ]
            }
        )
    )
    code = rvr.main(
        [
            str(dest / "record.json"),
            "--trust-file",
            str(trust_path),
            "--task-id",
            TASK,
            "--head-sha",
            SHA,
        ]
    )
    out = json.loads(capsys.readouterr().out)
    assert code == 0 and out["status"] == "verified" and out["key_id"] == key.key_id
    trust_path.write_text(json.dumps({"codex-worker": []}))
    code = rvr.main([str(dest / "record.json"), "--trust-file", str(trust_path)])
    out = json.loads(capsys.readouterr().out)
    assert code == 1 and out["status"] == "unverified"
    assert rvr.main([str(tmp_path / "missing.json")]) == 1


def test_codex_worker_is_an_optional_signer_and_documented():
    from infrastructure.validation import local_validation as lv

    assert "codex-worker" in attest.SIGNERS
    assert "codex-worker" not in lv.REQUIRED_SIGNERS
    committed = json.loads(lv.TRUSTED_KEYS_PATH.read_text())
    published = {s: m for s, m in attest.load_trusted().items() if m}
    assert lv.parse_trust_map(committed) == published
    text = README.read_text()
    for phrase in (
        "`codex-worker`",
        "V3-INFRA-006",
        "record.json",
        "sha256(public key)[:16]",
        "--print-public-key codex-worker",
    ):
        assert phrase in text, phrase


def test_verifier_requires_the_signed_experiment_id():
    key = attest.generate_keypair()
    trust = trust_of(key)
    missing = {k: v for k, v in base_record().items() if k != "experiment_id"}
    status, why = rvr.verify_review_record(signed(missing, key), trust)
    assert status == "invalid" and "experiment_id" in why
    other = signed({**base_record(), "experiment_id": "experiment-v2"}, key)
    status, why = rvr.verify_review_record(other, trust)
    assert status == "invalid" and "signed identity mismatch: experiment_id" in why
    assert rvr.verify_review_record(
        other, trust, expected_experiment_id="experiment-v2"
    ) == ("verified", "ok")
    assert rvr.current_experiment_id() == EXPERIMENT_ID
