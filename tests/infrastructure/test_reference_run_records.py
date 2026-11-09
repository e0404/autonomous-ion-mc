"""Signed, timestamped reference-run records (V3-INFRA-004).

End-to-end through run_reference with a fake engine, the production verifier
infrastructure/validation/reference_run_records.py, and its wiring into the
local-validation gate (details.reference_runs).
"""

import json
from pathlib import Path

import pytest

from infrastructure.experiment_v3 import artifacts, attest, reference
from infrastructure.experiment_v3.common import digest, file_hash
from infrastructure.validation import local_validation as lv
from infrastructure.validation import reference_run_records as rrr

SHA = "a" * 40
TASK = "V3-X1"


def _case(root, relative="case"):
    case = root / relative
    case.mkdir(parents=True)
    (case / "input").write_text("fixture")
    envelope = {
        k: "fixture"
        for k in (
            "physics",
            "source",
            "geometry",
            "materials",
            "cuts",
            "scorers",
            "units",
            "rationale",
        )
    }
    envelope.update(input="input", histories=1, seeds=[1], expected_outputs=["out"], arguments=[])
    (case / "case.json").write_text(json.dumps(envelope))


@pytest.fixture
def run(tmp_path, monkeypatch):
    """Run one fake-engine reference calculation; return (returned, run_dir, key)."""
    keys = tmp_path / "keys"
    monkeypatch.setattr(attest, "key_dir", lambda: keys)
    monkeypatch.setattr(reference, "require_bwrap", lambda: "/usr/bin/true")
    monkeypatch.setattr(reference, "exact_state", lambda root: SHA)
    monkeypatch.setattr(reference, "require_committed_inputs", lambda *args: None)
    monkeypatch.setattr(reference, "load_engine", lambda *args: {"engine": "topas"})
    # Fake engine: writes the expected output into the run's work directory.
    monkeypatch.setattr(
        reference,
        "command",
        lambda config, case, work, gpu=False: [
            "/bin/sh",
            "-c",
            f"echo hi; echo err >&2; printf data > '{work}/out'",
        ],
    )
    root = tmp_path / "repo"
    _case(root)
    state = tmp_path / "state"
    returned = reference.run_reference(root, "case", "topas", state=state, task_id=TASK)
    run_dir = state / "references" / returned["run_id"]
    key = attest.load_or_create_key("reference-runner")
    return returned, run_dir, key, state


def _trust(key):
    return {"reference-runner": {key.key_id: key.public_key.hex()}}


def _verify(record, key, **over):
    expected = {
        "expected_run_id": record["run_id"],
        "expected_code_sha": record["code_sha"],
        "expected_configuration_id": record["configuration_id"],
        "trust": _trust(key),
    }
    expected.update(over)
    return rrr.verify_reference_record(record, **expected)


def test_run_reference_returns_signed_timestamped_record(run):
    returned, run_dir, key, state = run
    on_disk = json.loads((run_dir / "record.json").read_text())
    assert returned["succeeded"] and returned == on_disk
    assert returned["schema_version"] == 3 and returned["task_id"] == TASK
    assert returned["code_sha"] == SHA and returned["engine"] == "topas"
    assert returned["attestation"]["signer"] == "reference-runner"
    assert returned["attestation"]["public_key"] == key.public_key.hex()
    attest.verify_record(on_disk, _trust(key), "reference-runner")
    assert returned["started_at"] <= returned["finished_at"]
    assert returned["seconds"] >= 0
    events = [json.loads(line) for line in (state / "events.jsonl").read_text().splitlines()]
    started = [e for e in events if e["event"] == "reference_started"]
    assert len(started) == 1
    assert started[0]["details"]["started_at"] == returned["started_at"]
    assert started[0]["details"]["run_id"] == returned["run_id"]
    assert returned["request_sha256"] == file_hash(run_dir / "request.json")
    assert returned["stdout_sha256"] == file_hash(run_dir / "stdout.txt")
    assert returned["stderr_sha256"] == file_hash(run_dir / "stderr.txt")
    assert returned["stdout"]["text"] == "hi\n"
    assert "output_files" not in returned and "artifact_files" not in returned


def test_result_index_lists_signed_record_and_binds_outputs(run):
    returned, run_dir, key, state = run
    result = json.loads((run_dir / "result.json").read_text())
    assert "attestation" not in result and result["signed_record"] == "record.json"
    assert {k: v for k, v in returned.items() if k != "attestation"}.items() <= (result.items())
    assert returned["output_files_sha256"] == digest(result["output_files"])
    assert returned["output_count"] == len(result["output_files"])
    assert result["output_files"]["out"]["bytes"] == 4
    assert result["artifact_files"]["record.json"]["sha256"] == file_hash(run_dir / "record.json")
    assert "result.json" not in result["artifact_files"]


def test_signed_record_is_downloadable_through_artifact_tools(run, tmp_path):
    returned, run_dir, key, state = run
    run_id = returned["run_id"]
    listing = artifacts.list_artifacts(run_id, limit=100, state=state)
    entry = next(f for f in listing["files"] if f["path"] == "record.json")
    assert entry["sha256"] == file_hash(run_dir / "record.json")
    preview = artifacts.read_text(run_id, "record.json", limit=4000, state=state)
    assert preview["sha256"] == entry["sha256"]
    worktree = tmp_path / "wt"
    worktree.mkdir()
    moved = artifacts.materialize(run_id, worktree, ["record.json"], state=state)
    copied = json.loads((Path(moved["path"]) / "record.json").read_text())
    assert copied == returned
    assert _verify(copied, key) == ("verified", "ok")


def test_verify_reference_record_and_loader(run):
    returned, run_dir, key, state = run
    loaded = rrr.load_signed_reference_record(returned["run_id"], state=state)
    assert loaded == returned
    assert _verify(loaded, key) == ("verified", "ok")
    status, reason = _verify(loaded, key, trust={})
    assert status == "unverified" and "not trusted" in reason
    with pytest.raises(rrr.ReferenceRecordMissing):
        rrr.load_signed_reference_record("REF-" + "0" * 20 + "-00000000", state=state)
    with pytest.raises(rrr.ReferenceRecordMissing):
        rrr.load_signed_reference_record("../" + returned["run_id"], state=state)


@pytest.mark.parametrize(
    "over",
    [
        {"expected_run_id": "REF-" + "0" * 20 + "-00000000"},
        {"expected_code_sha": "b" * 40},
        {"expected_configuration_id": "0" * 64},
    ],
)
@pytest.mark.parametrize("trusted", [True, False])
def test_replayed_identity_is_invalid(run, over, trusted):
    returned, _, key, _ = run
    if not trusted:
        over = {**over, "trust": {}}
    status, reason = _verify(returned, key, **over)
    assert status == "invalid" and "signed identity mismatch" in reason


def test_tampered_started_at_is_invalid(run):
    returned, _, key, _ = run
    tampered = {**returned, "started_at": "2000-01-01T00:00:00+00:00"}
    assert _verify(tampered, key)[0] == "invalid"


@pytest.mark.parametrize("value", [None, "", "yesterday", "2026-01-01T00:00:00"])
def test_missing_or_malformed_started_at_is_invalid_even_when_signed(run, value):
    returned, _, key, _ = run
    record = {k: v for k, v in returned.items() if k != "attestation"}
    if value is None:
        del record["started_at"]
    else:
        record["started_at"] = value
    resigned = attest.sign_record(record, key, "reference-runner")
    status, reason = _verify(resigned, key)
    assert status == "invalid" and "started_at" in reason
    assert _verify(resigned, key, trust={})[0] == "invalid"


def test_started_after_finished_is_invalid(run):
    returned, _, key, _ = run
    record = {k: v for k, v in returned.items() if k != "attestation"}
    record["started_at"] = "2999-01-01T00:00:00+00:00"
    status, reason = _verify(attest.sign_record(record, key, "reference-runner"), key)
    assert status == "invalid" and "earlier than started_at" in reason


def test_wrong_signer_or_schema_is_invalid(run, tmp_path):
    returned, _, key, _ = run
    record = {k: v for k, v in returned.items() if k != "attestation"}
    other = attest.sign_record(record, key, "host-runner")
    assert _verify(other, key)[0] == "invalid"
    old = attest.sign_record({**record, "schema_version": 2}, key, "reference-runner")
    assert _verify(old, key)[0] == "invalid"


# ------------------------------------------------------------ gate wiring --


@pytest.fixture
def gate(run, tmp_path, monkeypatch):
    returned, run_dir, rr, state = run
    vm = attest.generate_keypair()
    hr = attest.generate_keypair()
    monkeypatch.setattr(lv, "VALIDATION_ROOT", tmp_path / "val")
    monkeypatch.setattr(lv, "TRUSTED_KEYS_PATH", tmp_path / "trust.json")
    monkeypatch.setattr(lv, "trust_file_change", lambda s: ([], None))
    monkeypatch.setattr(rrr, "REFERENCE_STATE", state)
    monkeypatch.setattr(
        lv,
        "worktree_state",
        lambda t: {"path": tmp_path, "branch": "b", "sha": SHA, "dirty": False},
    )

    def trust(rr_keys=(rr,)):
        def ent(k):
            return {"key_id": k.key_id, "public_key": k.public_key.hex()}

        data = {"validation-manager": [ent(vm)], "host-runner": [ent(hr)]}
        if rr_keys is not None:
            data["reference-runner"] = [ent(k) for k in rr_keys]
        (tmp_path / "trust.json").write_text(json.dumps(data))

    def put_val(refs):
        rec = {
            "schema_version": 1,
            "validation_id": "VAL-1",
            "task_id": TASK,
            "suite": "s",
            "status": "passed",
            "summary": "ok",
            "git_sha": SHA,
            "details": {"reference_runs": refs},
        }
        d = lv.VALIDATION_ROOT / TASK / SHA
        d.mkdir(parents=True, exist_ok=True)
        (d / "v.json").write_text(json.dumps(attest.sign_record(rec, vm, "validation-manager")))

    trust()
    ref = {"run_id": returned["run_id"], "configuration_id": returned["configuration_id"]}
    return trust, put_val, ref


def _check(capsys):
    try:
        lv.check(TASK)
        code = 0
    except SystemExit as exc:
        code = exc.code
    return code, json.loads(capsys.readouterr().out)


def test_gate_accepts_verified_reference_run(gate, capsys):
    trust, put_val, ref = gate
    put_val([ref])
    code, out = _check(capsys)
    assert code == 0 and out["passing_records"] == 1


@pytest.mark.parametrize("rr_keys", [None, ()])
def test_gate_rejects_reference_run_while_key_unpublished(gate, capsys, rr_keys):
    trust, put_val, ref = gate
    trust(rr_keys)
    put_val([ref])
    code, out = _check(capsys)
    assert code == 1 and out["passing_records"] == 0
    assert "attestation required" in out["rejected_records"][0]["reason"]


@pytest.mark.parametrize("field", ["configuration_id", "run_id"])
def test_gate_rejects_replayed_reference_run(gate, capsys, field):
    trust, put_val, ref = gate
    bad = {
        **ref,
        field: "0" * 64 if field == "configuration_id" else "REF-" + "0" * 20 + "-00000000",
    }
    put_val([bad])
    code, out = _check(capsys)
    assert code == 1 and out["passing_records"] == 0


@pytest.mark.parametrize(
    "refs",
    ["x", [1], [{"run_id": ""}], [{"run_id": "REF-x", "configuration_id": "zz"}]],
)
def test_gate_and_record_reject_malformed_reference_declaration(gate, capsys, refs):
    trust, put_val, ref = gate
    put_val(refs)
    code, out = _check(capsys)
    assert code == 1 and "details.reference_runs" in out["rejected_records"][0]["reason"]
    assert lv.reference_run_declaration_problems({"reference_runs": refs})
    with pytest.raises(ValueError, match="malformed reference-run declaration"):
        lv.record(TASK, "s", "passed", "ok", json.dumps({"reference_runs": refs}))


def test_parse_trust_map_reference_runner_is_optional_but_strict_when_listed():
    vm, hr, rr, spare = (attest.generate_keypair() for _ in range(4))

    def ent(k):
        return {"key_id": k.key_id, "public_key": k.public_key.hex()}

    base = {"validation-manager": [ent(vm)], "host-runner": [ent(hr)]}
    expected = {
        "validation-manager": {vm.key_id: vm.public_key.hex()},
        "host-runner": {hr.key_id: hr.public_key.hex()},
    }
    assert lv.parse_trust_map(base, exactly_one=True) == expected
    assert lv.parse_trust_map({**base, "reference-runner": []}) == expected
    full = lv.parse_trust_map({**base, "reference-runner": [ent(rr)]}, exactly_one=True)
    assert full["reference-runner"] == {rr.key_id: rr.public_key.hex()}
    with pytest.raises(ValueError, match="exactly one reference-runner"):
        lv.parse_trust_map({**base, "reference-runner": [ent(rr), ent(spare)]}, exactly_one=True)
    with pytest.raises(ValueError, match="reference-runner"):
        lv.parse_trust_map({**base, "reference-runner": [{**ent(rr), "key_id": "0" * 16}]})
    with pytest.raises(ValueError, match="reference-runner"):
        lv.parse_trust_map({**base, "reference-runner": "not-a-list"})
    with pytest.raises(
        ValueError, match="reference-runner"
    ):  # explicit null is malformed (REVIEW-73ad5438)
        lv.parse_trust_map({**base, "reference-runner": None})
    with pytest.raises(ValueError, match="host-runner"):
        lv.parse_trust_map({"validation-manager": [ent(vm)], "reference-runner": [ent(rr)]})


def test_smoke_reports_keep_signed_record_unchanged(tmp_path, monkeypatch):
    """run_smokes must not mutate the signed record (REVIEW-0bf5f77d): the smoke
    verdict sits beside it and every record keeps a valid attestation."""
    from infrastructure.experiment_v3 import smoke

    monkeypatch.setattr(attest, "key_dir", lambda: tmp_path / "keys")
    monkeypatch.setattr(reference, "require_bwrap", lambda: "/usr/bin/true")
    monkeypatch.setattr(reference, "exact_state", lambda root: SHA)
    monkeypatch.setattr(reference, "require_committed_inputs", lambda *args: None)
    monkeypatch.setattr(reference, "load_engine", lambda name, reg: {"engine": name})
    # TOPAS-like text scorer output parses; the binary MHD engines fail scoring.
    monkeypatch.setattr(
        reference,
        "command",
        lambda config, case, work, gpu=False: [
            "/bin/sh",
            "-c",
            f"printf '0 0 0 1.5\\n' > '{work}/out'",
        ],
    )
    root = tmp_path / "repo"
    for engine in ("topas", "mcsquare", "fred"):
        _case(root, "experiment/v3/references/" + engine)
    report = smoke.run_smokes(root, tmp_path / "state")
    key = attest.load_or_create_key("reference-runner")
    engines = report["engines"]
    assert engines["topas"]["smoke"]["succeeded"] is True
    assert engines["topas"]["smoke"]["score_smoke"]["bins"] == 1
    for name in ("mcsquare", "fred"):
        assert engines[name]["smoke"]["succeeded"] is False
        assert engines[name]["smoke"]["score_error"]
        assert engines[name]["smoke"]["execution_succeeded"] is True
    assert report["all_execution_smokes_passed"] is False
    for name, entry in engines.items():
        record = entry["record"]
        assert record["succeeded"] is True and record["engine"] == name
        assert set(record) >= {"stdout", "stderr", "attestation", "started_at"}
        assert _verify(record, key) == ("verified", "ok")
        on_disk = json.loads(
            (tmp_path / "state/references" / record["run_id"] / "record.json").read_text()
        )
        assert on_disk == record
