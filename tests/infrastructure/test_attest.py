import importlib.util
import json
import os
import stat
import subprocess
from pathlib import Path

import pytest

from infrastructure.experiment_v3 import attest

ROOT = Path(__file__).resolve().parents[2]


def _load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


attest_verify = _load("attest_verify", "validation/scripts/transport/attest_verify.py")

# RFC 8032 section 7.1: (seed, public key, message, signature)
RFC_VECTORS = [
    (
        "9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60",
        "d75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a",
        "",
        "e5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e06522490155"
        "5fb8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b",
    ),
    (
        "4ccd089b28ff96da9db6c346ec114e0f5b8a319f35aba624da8cf6ed4fb8a6fb",
        "3d4017c3e843895a92b70aa74d1b7ebc9c982ccf2ec4968cc0cd55f12af4660c",
        "72",
        "92a009a9f0d4cab8720e820b5f642540a2b27b5416503f8fb3762223ebdb69da"
        "085ac1e43e15996e458f3613d0f11d8c387b2eaeb4302aeeb00d291612bb0c00",
    ),
    (
        "c5aa8df43f9f837bedb7442f31dcb7b166d38535076f094b85ce3a2e0b4458f7",
        "fc51cd8e6218a1a38da47ed00230f0580816ed13ba3303ac5deb911548908025",
        "af82",
        "6291d657deec24024827e69c3abe01a30ce548a284743a445e3680d7db5ac3ac"
        "18ff9b538d16f290ae67f760984dc6594a7c15e9716ed28dc027beceea1ec40a",
    ),
]


@pytest.mark.parametrize("seed,pub,msg,sig", RFC_VECTORS)
def test_rfc8032_vectors(seed, pub, msg, sig):
    seed, msg = bytes.fromhex(seed), bytes.fromhex(msg)
    derived = attest.public_key_from_seed(seed)
    assert derived.hex() == pub
    assert attest.sign(msg, seed).hex() == sig
    assert attest.verify(msg, bytes.fromhex(sig), derived)


def test_round_trip_and_tamper():
    key = attest.generate_keypair()
    sig = attest.sign(b"message", key.seed)
    assert attest.verify(b"message", sig, key.public_key)
    assert not attest.verify(b"messagf", sig, key.public_key)
    assert not attest.verify(b"message", sig[:-1] + bytes([sig[-1] ^ 1]), key.public_key)
    assert not attest.verify(b"message", sig, attest.generate_keypair().public_key)
    assert not attest.verify(b"message", b"short", key.public_key)
    assert not attest.verify(b"message", sig, b"\x00" * 31)


def test_canonicalisation_stable():
    a = {"b": 1, "a": {"y": 2, "x": "é"}}
    b = {"a": {"x": "é", "y": 2}, "b": 1, "attestation": {"anything": 1}}
    assert attest.canonical_payload(a) == attest.canonical_payload(b)
    assert attest.canonical_payload(a) == b'{"a":{"x":"\\u00e9","y":2},"b":1}'


def test_key_persistence_and_mode(tmp_path):
    key = attest.load_or_create_key("host-runner", tmp_path)
    again = attest.load_or_create_key("host-runner", tmp_path)
    assert key == again
    assert stat.S_IMODE(os.stat(tmp_path / "host-runner.key").st_mode) == 0o600
    assert (tmp_path / "host-runner.pub").read_text().strip() == key.public_key.hex()
    os.chmod(tmp_path / "host-runner.key", 0o644)
    with pytest.raises(attest.AttestationError):
        attest.load_or_create_key("host-runner", tmp_path)
    with pytest.raises(attest.AttestationError):
        attest.load_or_create_key("nobody", tmp_path)


def _signed(signer="validation-manager"):
    key = attest.generate_keypair()
    record = {"id": "X", "n": 1.5, "text": "ünï"}
    return key, record, attest.sign_record(record, key, signer)


def test_verify_record_trust_and_tamper():
    key, record, signed = _signed()
    assert "attestation" not in record
    block = signed["attestation"]
    assert set(block) == {"schema", "alg", "signer", "key_id", "public_key",
                          "signed_at_utc", "payload_sha256", "signature"}
    trusted = {"validation-manager": {key.key_id: key.public_key.hex()}}
    ok = attest.verify_record(signed, trusted, expected_signer="validation-manager")
    assert ok.key_id == key.key_id
    with pytest.raises(attest.AttestationError, match="not trusted"):
        attest.verify_record(signed, {})
    with pytest.raises(attest.AttestationError, match="not trusted"):
        attest.verify_record(signed, {"host-runner": trusted["validation-manager"]})
    for field, value in (("signer", "host-runner"), ("signed_at_utc", "2000-01-01T00:00:00+00:00")):
        bad = json.loads(json.dumps(signed))
        bad["attestation"][field] = value
        with pytest.raises(attest.AttestationError, match="invalid signature"):
            attest.verify_record(bad, {"host-runner": trusted["validation-manager"], **trusted})
        both = {"host-runner": trusted["validation-manager"], **trusted}
        assert attest_verify.verify_record(bad, bad["attestation"]["signer"], both)["status"] == "invalid"
    with pytest.raises(attest.AttestationError, match="signer"):
        attest.verify_record(signed, trusted, expected_signer="host-runner")
    with pytest.raises(attest.AttestationError, match="no attestation"):
        attest.verify_record(record, trusted)
    tampered = {**signed, "n": 2}
    with pytest.raises(attest.AttestationError, match="payload hash"):
        attest.verify_record(tampered, trusted)
    # attacker recomputes the hash but cannot forge the signature
    import hashlib
    forged = json.loads(json.dumps(tampered))
    forged["attestation"]["payload_sha256"] = hashlib.sha256(
        attest.canonical_payload(forged)).hexdigest()
    with pytest.raises(attest.AttestationError, match="invalid signature"):
        attest.verify_record(forged, trusted)
    # attacker substitutes own key: not trusted
    other = attest.generate_keypair()
    resigned = attest.sign_record(tampered, other, "validation-manager")
    with pytest.raises(attest.AttestationError, match="not trusted"):
        attest.verify_record(resigned, trusted)


def test_repo_verifier_agrees():
    key, record, signed = _signed("host-runner")
    trusted = {"host-runner": {key.key_id: key.public_key.hex()}}
    assert attest_verify.verify_record(signed, "host-runner", trusted)["status"] == "verified"
    assert attest_verify.verify_record(signed, "host-runner", {"host-runner": {}})["status"] == "unverified"
    assert attest_verify.verify_record(record, "host-runner", trusted)["status"] == "unverified"
    assert attest_verify.verify_record(signed, "validation-manager", trusted)["status"] == "invalid"
    assert attest_verify.verify_record({**signed, "n": 9}, "host-runner", trusted)["status"] == "invalid"
    assert attest.verify_record(signed, trusted, "host-runner").key_id == key.key_id
    sig = bytes.fromhex(signed["attestation"]["signature"])
    assert attest_verify.verify(b"m", attest.sign(b"m", key.seed), key.public_key)
    assert not attest_verify.verify(b"m", sig, key.public_key)


def test_committed_trust_file_is_empty_or_strictly_published():
    """Two permitted states: empty (trust on first use, before the key-publication
    commit) or exactly one valid, non-small-order key per signer (after it)."""
    from infrastructure.validation import local_validation as lv

    trusted = attest_verify.load_trusted()
    assert attest.load_trusted() == trusted
    # Required signers are always listed; an optional signer (reference-runner,
    # V3-INFRA-004) may be absent until its key-publication commit.
    assert set(lv.REQUIRED_SIGNERS) <= set(trusted) <= set(attest.SIGNERS)
    data = json.loads(lv.TRUSTED_KEYS_PATH.read_text())
    if all(not trusted[s] for s in lv.REQUIRED_SIGNERS):
        assert {k: v for k, v in data.items() if v} == {}
    else:
        published = {s: m for s, m in trusted.items() if m}
        assert lv.parse_trust_map(data, exactly_one=True) == published


def test_cli_print_public_key(tmp_path, monkeypatch):
    monkeypatch.setattr(attest, "key_dir", lambda: tmp_path)
    attest.main(["--print-public-key", "host-runner"])
    assert (tmp_path / "host-runner.pub").exists()


def test_local_validation_signs_record(tmp_path, monkeypatch, capsys):
    from infrastructure.validation import local_validation as lv

    keys = tmp_path / "keys"
    monkeypatch.setattr(attest, "key_dir", lambda: keys)
    monkeypatch.setattr(lv, "VALIDATION_ROOT", tmp_path / "val")
    monkeypatch.setattr(lv, "worktree_state", lambda t: {
        "path": tmp_path, "branch": "task/x", "sha": "a" * 40, "dirty": False})
    lv.record("V3-X1", "suite", "passed", "ok", '{"k": 1}')
    out = json.loads(capsys.readouterr().out)
    returned = out["validation"]
    on_disk = json.loads(Path(out["record_file"]).read_text())
    assert returned == on_disk
    pub = (keys / "validation-manager.pub").read_text().strip()
    trusted = {returned["attestation"]["key_id"]: pub}
    attest.verify_record(on_disk, {"validation-manager": trusted}, expected_signer="validation-manager")
    monkeypatch.setattr(lv, "TRUSTED_KEYS_PATH", tmp_path / "none.json")
    (tmp_path / "none.json").write_text('{"validation-manager": [], "host-runner": []}')
    loaded = lv.records_for("V3-X1", "a" * 40)[0]
    assert loaded.pop("attestation_status") == "unverified"
    loaded.pop("attestation_reason")
    assert loaded.pop("host_run_attestation") == []
    assert loaded.pop("reference_run_attestation") == []
    assert loaded == on_disk


def test_host_runner_signs_record(tmp_path, monkeypatch):
    hr = _load("ionmc_host_runner_attest", "infrastructure/host_runner/host_runner.py")
    keys = tmp_path / "keys"
    monkeypatch.setattr(attest, "key_dir", lambda: keys)
    monkeypatch.setattr(hr, "RUN_ROOT", tmp_path / "runs")
    monkeypatch.setattr(hr, "CACHE_ROOT", tmp_path / "cache")
    wt = tmp_path / "wt"
    wt.mkdir()
    monkeypatch.setattr(hr, "inspect_worktree", lambda t: {"path": wt, "branch": "task/x", "sha": "b" * 40})
    monkeypatch.setattr(hr, "git_text", lambda cwd, *a: "b" * 40 if a[0] == "rev-parse" else "")
    import infrastructure.experiment_v3.host_snapshot as snap
    monkeypatch.setattr(snap, "create", lambda *a: tmp_path / "snap")
    monkeypatch.setattr(hr, "build_bwrap_command", lambda *a: ["sh", "-c", "echo hi; echo err >&2"])
    result = hr.run_validation("V3-X1", ["true"], 30)
    run_dir = Path(result["run_directory"])
    on_disk = json.loads((run_dir / "record.json").read_text())
    assert on_disk == json.loads(json.dumps(result))
    import hashlib
    assert result["stdout_sha256"] == hashlib.sha256((run_dir / "stdout.txt").read_bytes()).hexdigest()
    assert result["stderr_sha256"] == hashlib.sha256((run_dir / "stderr.txt").read_bytes()).hexdigest()
    assert result["output_manifest_sha256"] is None
    for field in ("run_id", "task_id", "sha", "argv", "started_at", "finished_at",
                  "duration_seconds", "exit_code", "timed_out"):
        assert field in result
    pub = (keys / "host-runner.pub").read_text().strip()
    attest.verify_record(
        on_disk, {"host-runner": {result["attestation"]["key_id"]: pub}}, expected_signer="host-runner"
    )


def test_output_manifest_hash(tmp_path):
    hr = _load("ionmc_host_runner_manifest", "infrastructure/host_runner/host_runner.py")
    (tmp_path / "d").mkdir()
    (tmp_path / "d" / "a.txt").write_text("x")
    first = hr.output_manifest_sha256(tmp_path)
    (tmp_path / "d" / "a.txt").write_text("y")
    assert first and first != hr.output_manifest_sha256(tmp_path)
    assert hr.output_manifest_sha256(tmp_path / "missing") is None


def test_mcp_modules_import():
    pytest.importorskip("mcp")
    import infrastructure.mcp.host_runner_mcp  # noqa: F401
    import infrastructure.mcp.validation_mcp  # noqa: F401


IDENTITY_ENC = (1).to_bytes(32, "little")  # encoding of the identity point (0, 1)
# The other small-order points: y = 0 (order 4), y = -1 (order 2), and the order-8 points.
_SMALL_Y = [0, attest._P - 1, (1 + attest._I) * attest._inv(1 - attest._I) % attest._P]


def _small_order_encodings():
    out = [IDENTITY_ENC]
    for y in _SMALL_Y:
        for sign in (0, 1):
            enc = (y | (sign << 255)).to_bytes(32, "little")
            if attest._decompress(enc) is not None:
                out.append(enc)
    return out


@pytest.mark.parametrize("verifier", [attest, attest_verify])
def test_small_order_public_key_cannot_forge(verifier):
    """Identity (and every small-order point) as public key + R with S = 0 satisfies the
    bare equation for any message; both verifiers must reject it."""
    forged = IDENTITY_ENC + (0).to_bytes(32, "little")
    assert not verifier.verify(b"any message", forged, IDENTITY_ENC)
    for enc in _small_order_encodings():
        assert verifier._small_order(verifier._decompress(enc))
        assert not verifier.verify(b"m", enc + (0).to_bytes(32, "little"), enc)
        # small-order public key with an otherwise well-formed signature
        key = attest.generate_keypair()
        sig = attest.sign(b"m", key.seed)
        assert not verifier.verify(b"m", sig, enc)
        # small-order R with an honest public key
        assert not verifier.verify(b"m", enc + sig[32:], key.public_key)
    key = attest.generate_keypair()
    assert verifier.verify(b"m", attest.sign(b"m", key.seed), key.public_key)
    assert not verifier._small_order(verifier._decompress(key.public_key))


def test_trust_map_rejects_small_order_key():
    from infrastructure.validation import local_validation as lv

    hr = attest.generate_keypair()
    for enc in _small_order_encodings():
        data = {
            "validation-manager": [{"key_id": attest.key_id_of(enc), "public_key": enc.hex()}],
            "host-runner": [{"key_id": hr.key_id, "public_key": hr.public_key.hex()}],
        }
        with pytest.raises(ValueError, match="small-order"):
            lv.parse_trust_map(data)


# --- require_verified: migration-only bootstrap and the output-manifest expectation ---------------
MIGRATION = "V3-INFRA-002"


def _run_record(task=MIGRATION, manifest="m" * 64, drop_manifest=False):
    rec = {"run_id": "RUN-1", "task_id": task, "sha": "a" * 40, "output_manifest_sha256": manifest}
    if drop_manifest:
        del rec["output_manifest_sha256"]
    return rec


def test_bootstrap_tasks_mirror_local_validation():
    """The stdlib-only verifier mirrors (does not import) the gate's BOOTSTRAP_TASKS."""
    from infrastructure.validation import local_validation as lv

    assert attest_verify.BOOTSTRAP_TASKS == lv.BOOTSTRAP_TASKS == frozenset({MIGRATION})
    src = (ROOT / "validation/scripts/transport/attest_verify.py").read_text(encoding="utf-8")
    assert "import infrastructure" not in src and "from infrastructure" not in src


def test_bootstrap_requires_expected_task_id():
    with pytest.raises(attest_verify.AttestationRequired, match="reserved for the migration"):
        attest_verify.require_verified(_run_record(), "host-runner", {}, bootstrap=True)


def test_bootstrap_rejected_for_non_migration_task():
    rec = _run_record(task="V3-008")
    for trust in ({}, {"host-runner": {}}):
        with pytest.raises(attest_verify.AttestationRequired, match="reserved for the migration"):
            attest_verify.require_verified(
                rec, "host-runner", trust, expected_task_id="V3-008", bootstrap=True
            )
    # without bootstrap the unsigned record is rejected as well
    with pytest.raises(attest_verify.AttestationRequired, match="unverified record rejected"):
        attest_verify.require_verified(rec, "host-runner", {}, expected_task_id="V3-008")


def test_bootstrap_migration_task_with_empty_trust_map_passes():
    for trust in ({}, {"host-runner": {}}):
        out = attest_verify.require_verified(
            _run_record(), "host-runner", trust, expected_task_id=MIGRATION, bootstrap=True
        )
        assert out["status"] == "unverified"
    # the task binding still applies in bootstrap
    with pytest.raises(attest_verify.AttestationRequired, match="identity mismatch"):
        attest_verify.require_verified(
            _run_record(task="V3-008"), "host-runner", {}, expected_task_id=MIGRATION,
            bootstrap=True,
        )  # fmt: skip


def test_bootstrap_migration_task_with_trusted_key_requires_verification():
    key = attest.generate_keypair()
    other = attest.generate_keypair()
    trust = {"host-runner": {key.key_id: key.public_key.hex()}}
    kw = {"expected_task_id": MIGRATION, "bootstrap": True}
    with pytest.raises(attest_verify.AttestationRequired, match="unverified record rejected"):
        attest_verify.require_verified(_run_record(), "host-runner", trust, **kw)
    untrusted = attest.sign_record(_run_record(), other, "host-runner")
    with pytest.raises(attest_verify.AttestationRequired, match="not trusted"):
        attest_verify.require_verified(untrusted, "host-runner", trust, **kw)
    signed = attest.sign_record(_run_record(), key, "host-runner")
    assert attest_verify.require_verified(signed, "host-runner", trust, **kw)["status"] == "verified"


def _signed_trust(rec):
    key = attest.generate_keypair()
    return attest.sign_record(rec, key, "host-runner"), {"host-runner": {key.key_id: key.public_key.hex()}}


def test_output_manifest_expectation_omitted():
    """OMIT (the default) checks nothing: null, hex and a missing field all pass."""
    assert repr(attest_verify.OMIT) == "OMIT"
    for rec in (_run_record(manifest=None), _run_record(), _run_record(drop_manifest=True)):
        signed, trust = _signed_trust(rec)
        assert attest_verify.require_verified(signed, "host-runner", trust)["status"] == "verified"
        out = attest_verify.require_verified(
            signed, "host-runner", trust, expected_output_manifest_sha256=attest_verify.OMIT
        )
        assert out["status"] == "verified"


def test_output_manifest_explicit_none_matches_signed_null():
    signed, trust = _signed_trust(_run_record(manifest=None))
    out = attest_verify.require_verified(
        signed, "host-runner", trust, expected_output_manifest_sha256=None
    )
    assert out["status"] == "verified"


def test_output_manifest_explicit_none_vs_non_null_raises():
    signed, trust = _signed_trust(_run_record(manifest="m" * 64))
    with pytest.raises(attest_verify.AttestationRequired, match="output_manifest_sha256 is 'm"):
        attest_verify.require_verified(
            signed, "host-runner", trust, expected_output_manifest_sha256=None
        )


def test_output_manifest_explicit_hex_vs_null_raises():
    signed, trust = _signed_trust(_run_record(manifest=None))
    with pytest.raises(attest_verify.AttestationRequired, match="output_manifest_sha256 is None"):
        attest_verify.require_verified(
            signed, "host-runner", trust, expected_output_manifest_sha256="m" * 64
        )
    signed, trust = _signed_trust(_run_record(manifest="m" * 64))
    assert (
        attest_verify.require_verified(
            signed, "host-runner", trust, expected_output_manifest_sha256="m" * 64
        )["status"]
        == "verified"
    )


def test_output_manifest_missing_field_raises_when_expected():
    signed, trust = _signed_trust(_run_record(drop_manifest=True))
    for want in (None, "m" * 64):
        with pytest.raises(attest_verify.AttestationRequired, match="is missing"):
            attest_verify.require_verified(
                signed, "host-runner", trust, expected_output_manifest_sha256=want
            )
