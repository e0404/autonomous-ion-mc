"""Key-publication mode of the local-validation gate (trust bootstrap after merge).

End-to-end against a real git repository: an integration branch with an empty
trust file and a task branch whose only change is the published trust file.
"""

import json
import subprocess

import pytest

from infrastructure.experiment_v3 import attest
from infrastructure.validation import host_run_records as hrr
from infrastructure.validation import local_validation as lv
from infrastructure.validation import reference_run_records as rrr

TASK = "V3-INFRA-003"
RUN = "RUN-20260101-000000-KEY001"
MANIFEST = "cd" * 32
BASE = "v3/develop"
GIT_ENV = {
    "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@x", "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@x", "HOME": "/nonexistent", "GIT_CONFIG_NOSYSTEM": "1",
}


def git(*args, cwd):
    return subprocess.run(
        ["git", *args], cwd=cwd, text=True, capture_output=True, check=True, env=GIT_ENV
    ).stdout.strip()


def entry(key, **extra):
    return {"key_id": key.key_id, "public_key": key.public_key.hex(), **extra}


def trust_file(vm_keys=(), hr_keys=(), rr_keys=None):
    data = {
        "validation-manager": [entry(k) for k in vm_keys],
        "host-runner": [entry(k) for k in hr_keys],
    }
    if rr_keys is not None:  # optional signer (V3-INFRA-004)
        data["reference-runner"] = [entry(k) for k in rr_keys]
    return json.dumps(data, indent=2) + "\n"


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """Git repo: base branch with empty trust file; task branch created by the test."""
    r = tmp_path / "repo"
    r.mkdir()
    git("init", "-q", "-b", BASE, cwd=r)
    key = r / lv.KEY_FILE_REL
    key.parent.mkdir(parents=True)
    key.write_text(trust_file())
    (r / "other.txt").write_text("base\n")
    git("add", "-A", cwd=r)
    git("commit", "-q", "-m", "base", cwd=r)
    git("checkout", "-q", "-b", "task/key-pub", cwd=r)

    integration_trust = tmp_path / "integration_trust.json"
    integration_trust.write_text(trust_file())
    monkeypatch.setattr(lv, "TRUSTED_KEYS_PATH", integration_trust)
    monkeypatch.setattr(lv, "VALIDATION_ROOT", tmp_path / "val")
    monkeypatch.setattr(hrr, "HOST_RUNS_ROOT", tmp_path / "runs")
    monkeypatch.setattr(rrr, "REFERENCE_STATE", tmp_path / "state")
    monkeypatch.setattr(lv, "base_ref", lambda: BASE)
    monkeypatch.setattr(
        lv, "worktree_state",
        lambda t: {"path": r, "branch": "task/key-pub",
                   "sha": git("rev-parse", "HEAD", cwd=r), "dirty": False},
    )
    return r


def commit_branch(r, text, extra_file=None):
    (r / lv.KEY_FILE_REL).write_text(text)
    if extra_file:
        (r / extra_file).write_text("changed\n")
    git("add", "-A", cwd=r)
    git("commit", "-q", "--allow-empty", "-m", "publish keys", cwd=r)
    return git("rev-parse", "HEAD", cwd=r)


CONFIG = "ef" * 32
REF = "REF-" + CONFIG[:20] + "-0123abcd"


def put_reference_run(tmp_path, sha, rr):
    rec = {
        "schema_version": 3, "run_id": REF, "configuration_id": CONFIG,
        "code_sha": sha, "task_id": TASK, "engine": "topas",
        "started_at": "2026-01-01T00:00:00+00:00",
        "finished_at": "2026-01-01T00:01:00+00:00", "succeeded": True,
    }
    d = tmp_path / "state" / "references" / REF
    d.mkdir(parents=True, exist_ok=True)
    (d / "record.json").write_text(json.dumps(attest.sign_record(rec, rr, "reference-runner")))


def put_records(
    tmp_path, sha, vm, hr, sign_val=True, sign_run=True, with_run=True, with_ref=False
):
    run = {
        "schema_version": 1, "run_id": RUN, "task_id": TASK, "sha": sha,
        "exit_code": 0, "timed_out": False, "succeeded": True,
        "committed_state_unchanged": True, "dirty_after": False,
        "output_publication_error": None, "output_manifest_sha256": MANIFEST,
    }
    if sign_run:
        run = attest.sign_record(run, hr, "host-runner")
    d = tmp_path / "runs" / TASK / RUN
    d.mkdir(parents=True, exist_ok=True)
    (d / "record.json").write_text(json.dumps(run))
    val = {
        "schema_version": 1, "validation_id": "VAL-1", "task_id": TASK, "suite": "s",
        "status": "passed", "summary": "ok", "git_sha": sha,
        "details": {"host_runs": [{"run_id": RUN, "output_manifest_sha256": MANIFEST}]},
    }
    if not with_run:
        val["details"] = {}
    if with_ref:
        val["details"]["reference_runs"] = [{"run_id": REF, "configuration_id": CONFIG}]
    if sign_val:
        val = attest.sign_record(val, vm, "validation-manager")
    d = tmp_path / "val" / TASK / sha
    d.mkdir(parents=True, exist_ok=True)
    (d / "v.json").write_text(json.dumps(val))


def check(capsys):
    try:
        lv.check(TASK)
        code = 0
    except SystemExit as exc:
        code = exc.code
    return code, json.loads(capsys.readouterr().out)


def test_pure_key_publication_branch_verifies_against_its_own_trust_file(repo, tmp_path, capsys):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 0 and out["status"] == "passed"
    assert out["trust_source"] == "key-publication"
    assert out["enforcement"] == "enforcing" and out["verified_records"] == 1
    assert out["key_publication_evidence"] == {"validation-manager": True, "host-runner": True}


def test_key_publication_requires_host_run_evidence_for_the_host_runner_key(
    repo, tmp_path, capsys
):
    """A verified validation record alone would let an arbitrary host-runner key through."""
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, vm, hr, with_run=False)
    code, out = check(capsys)
    assert code == 1 and out["status"] == "failed"
    assert out["trust_source"] == "key-publication" and out["passing_records"] == 1
    assert out["key_publication_evidence"] == {"validation-manager": True, "host-runner": False}
    assert "host-runner" in out["key_publication_problem"]


@pytest.mark.parametrize("extra", ["validation-manager", "host-runner"])
def test_key_publication_rejects_more_than_one_key_per_signer(repo, tmp_path, capsys, extra):
    vm, hr, spare = (attest.generate_keypair() for _ in range(3))
    vms, hrs = ([vm, spare], [hr]) if extra == "validation-manager" else ([vm], [hr, spare])
    sha = commit_branch(repo, trust_file(vms, hrs))
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 1 and out["trust_source"] == "integration" and out["passing_records"] == 0


def test_key_publication_mode_requires_every_record_verified(repo, tmp_path, capsys):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, vm, hr, sign_val=False)
    code, out = check(capsys)
    assert code == 1 and out["trust_source"] == "key-publication"
    assert out["passing_records"] == 0 and out["rejected_records"]


def test_key_publication_mode_requires_host_run_verified(repo, tmp_path, capsys):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, vm, hr, sign_run=False)
    code, out = check(capsys)
    assert code == 1 and "attestation required" in out["rejected_records"][0]["reason"]


def test_record_signed_by_key_not_in_branch_trust_is_rejected(repo, tmp_path, capsys):
    vm, hr, rogue = (attest.generate_keypair() for _ in range(3))
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, rogue, hr)
    code, out = check(capsys)
    assert code == 1 and out["passing_records"] == 0


def test_branch_with_any_other_change_is_not_key_publication(repo, tmp_path, capsys):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [hr]), extra_file="other.txt")
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 1 and out["trust_source"] == "integration"
    assert out["enforcement"] == "blocked"
    assert lv.ATTESTATION_REQUIRED in out["rejected_records"][0]["reason"]
    assert lv.TRUST_CHANGE_RULE in out["trust_problem"]


def test_branch_with_new_file_is_not_key_publication(repo, tmp_path, capsys):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [hr]), extra_file="new.txt")
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 1 and out["trust_source"] == "integration"


@pytest.mark.parametrize("text", [
    trust_file([], []),
    "not json",
    json.dumps({"validation-manager": [], "host-runner": []}),
    "VM_ONLY",
    "BAD_KEY_ID",
    "SHORT_KEY",
])
def test_malformed_or_incomplete_branch_trust_file_disables_mode(repo, tmp_path, capsys, text):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    if text == "VM_ONLY":
        text = trust_file([vm], [])
    elif text == "BAD_KEY_ID":
        data = json.loads(trust_file([vm], [hr]))
        data["host-runner"][0]["key_id"] = "0" * 16
        text = json.dumps(data)
    elif text == "SHORT_KEY":
        data = json.loads(trust_file([vm], [hr]))
        data["validation-manager"][0]["public_key"] = "ab"
        text = json.dumps(data)
    sha = commit_branch(repo, text)
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 1 and out["trust_source"] == "integration"
    assert out["passing_records"] == 0


def test_rotation_pure_branch_is_evaluated_against_its_own_keys(repo, tmp_path, capsys):
    """Integration already trusts keys; a pure branch replacing them (rotation) must
    evidence the NEW keys."""
    vm, hr, old_vm, old_hr = (attest.generate_keypair() for _ in range(4))
    lv.TRUSTED_KEYS_PATH.write_text(trust_file([old_vm], [old_hr]))
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 0 and out["trust_source"] == "key-publication"
    assert out["key_publication_evidence"] == {"validation-manager": True, "host-runner": True}


def test_rotation_signed_by_the_old_keys_is_rejected(repo, tmp_path, capsys):
    vm, hr, old_vm, old_hr = (attest.generate_keypair() for _ in range(4))
    lv.TRUSTED_KEYS_PATH.write_text(trust_file([old_vm], [old_hr]))
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, old_vm, old_hr)  # the branch publishes keys it cannot evidence
    code, out = check(capsys)
    assert code == 1 and out["trust_source"] == "key-publication"
    assert out["passing_records"] == 0


def test_trust_file_change_with_other_changes_fails_even_when_enforcing(repo, tmp_path, capsys):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    lv.TRUSTED_KEYS_PATH.write_text(trust_file([vm], [hr]))
    rogue = attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [rogue]), extra_file="other.txt")
    put_records(tmp_path, sha, vm, hr)  # valid under the integration keys
    code, out = check(capsys)
    assert code == 1 and out["status"] == "failed"
    assert out["trust_source"] == "integration" and out["passing_records"] == 1
    assert lv.TRUST_CHANGE_RULE in out["trust_problem"]


def test_appending_a_second_key_fails_even_when_enforcing(repo, tmp_path, capsys):
    vm, hr, extra = (attest.generate_keypair() for _ in range(3))
    lv.TRUSTED_KEYS_PATH.write_text(trust_file([vm], [hr]))
    sha = commit_branch(repo, trust_file([vm], [hr, extra]))
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 1 and out["trust_source"] == "integration"
    assert lv.TRUST_CHANGE_RULE in out["trust_problem"]


def test_branch_not_touching_trust_file_uses_integration_keys(repo, tmp_path, capsys):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    lv.TRUSTED_KEYS_PATH.write_text(trust_file([vm], [hr]))
    (repo / "other.txt").write_text("work\n")
    git("add", "-A", cwd=repo)
    git("commit", "-q", "-m", "work", cwd=repo)
    sha = git("rev-parse", "HEAD", cwd=repo)
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 0 and out["trust_source"] == "integration"
    assert "trust_problem" not in out


def test_no_bootstrap_tolerance_for_migration_task_in_key_publication_mode(
    repo, tmp_path, capsys, monkeypatch
):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, vm, hr, sign_val=False)
    # Same records filed under the migration task id.
    src = tmp_path / "val" / TASK
    (tmp_path / "val" / "V3-INFRA-002").symlink_to(src)
    (tmp_path / "runs" / "V3-INFRA-002").symlink_to(tmp_path / "runs" / TASK)
    for p in (tmp_path / "runs" / TASK / RUN / "record.json", src / sha / "v.json"):
        rec = json.loads(p.read_text())
        rec["task_id"] = "V3-INFRA-002"
        p.write_text(json.dumps(rec))
    try:
        lv.check("V3-INFRA-002")
        code = 0
    except SystemExit as exc:
        code = exc.code
    out = json.loads(capsys.readouterr().out)
    assert code == 1 and out["trust_source"] == "key-publication"
    assert out["enforcement"] == "enforcing" and out["passing_records"] == 0


def test_parse_trust_map_round_trip():
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    parsed = lv.parse_trust_map(json.loads(trust_file([vm], [hr])))
    assert parsed == {
        "validation-manager": {vm.key_id: vm.public_key.hex()},
        "host-runner": {hr.key_id: hr.public_key.hex()},
    }


def test_undeterminable_branch_diff_fails_closed(repo, tmp_path, capsys, monkeypatch):
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    lv.TRUSTED_KEYS_PATH.write_text(trust_file([vm], [hr]))
    sha = commit_branch(repo, trust_file([vm], [hr]))
    put_records(tmp_path, sha, vm, hr)
    monkeypatch.setattr(lv, "base_ref", lambda: "refs/heads/does-not-exist")
    code, out = check(capsys)
    assert code == 1 and "cannot determine the branch changes" in out["trust_problem"]


def test_reference_runner_key_publication_requires_reference_run_evidence(
    repo, tmp_path, capsys
):
    """Publishing the optional reference-runner key (V3-INFRA-004) follows the same
    rule: exactly one key, evidenced by a verified reference run of this SHA."""
    vm, hr, rr = (attest.generate_keypair() for _ in range(3))
    lv.TRUSTED_KEYS_PATH.write_text(trust_file([vm], [hr]))
    sha = commit_branch(repo, trust_file([vm], [hr], [rr]))
    put_reference_run(tmp_path, sha, rr)
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 1 and out["trust_source"] == "key-publication"
    assert out["key_publication_evidence"] == {
        "validation-manager": True, "host-runner": True, "reference-runner": False,
    }
    assert "reference-runner" in out["key_publication_problem"]
    put_records(tmp_path, sha, vm, hr, with_ref=True)
    code, out = check(capsys)
    assert code == 0 and out["status"] == "passed"
    assert out["key_publication_evidence"]["reference-runner"] is True


def test_reference_runner_key_publication_rejects_unevidenced_rogue_key(
    repo, tmp_path, capsys
):
    vm, hr, rr, rogue = (attest.generate_keypair() for _ in range(4))
    lv.TRUSTED_KEYS_PATH.write_text(trust_file([vm], [hr]))
    sha = commit_branch(repo, trust_file([vm], [hr], [rogue]))
    put_reference_run(tmp_path, sha, rr)  # signed by the service key, not the rogue one
    put_records(tmp_path, sha, vm, hr, with_ref=True)
    code, out = check(capsys)
    assert code == 1 and out["passing_records"] == 0


def test_unpublished_reference_runner_does_not_affect_existing_signers(
    repo, tmp_path, capsys
):
    """An absent or empty reference-runner list keeps the INFRA-002/003 behaviour."""
    vm, hr = attest.generate_keypair(), attest.generate_keypair()
    sha = commit_branch(repo, trust_file([vm], [hr], []))
    put_records(tmp_path, sha, vm, hr)
    code, out = check(capsys)
    assert code == 0 and out["trust_source"] == "key-publication"
    assert out["key_publication_evidence"] == {"validation-manager": True, "host-runner": True}
