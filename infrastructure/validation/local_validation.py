#!/usr/bin/env python3

import argparse
import json
import os
import re
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path


WORKTREE_ROOT = Path(__file__).resolve().parents[2].parent / (Path(__file__).resolve().parents[2].name + "-worktrees")

VALIDATION_ROOT = Path(
    os.environ.get(
        "IONMC_VALIDATION_DIR",
        str(Path.home() / ".local/share/ionmc-experiment/validation"),
    )
).expanduser()


def git_text(*args: str, cwd: Path) -> str:
    return subprocess.check_output(
        ["git", *args],
        cwd=cwd,
        text=True,
        stderr=subprocess.PIPE,
    ).strip()


def normalize_task_id(task_id: str) -> str:
    value = task_id.strip().upper()

    if not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{1,63}", value):
        raise ValueError(f"Invalid task ID: {task_id}")

    return value


def worktree_for(task_id: str) -> Path:
    return WORKTREE_ROOT / normalize_task_id(task_id).lower()


def worktree_state(task_id: str) -> dict:
    path = worktree_for(task_id)

    if not path.exists():
        raise RuntimeError(f"Task worktree does not exist: {path}")

    status = git_text("status", "--porcelain", cwd=path)

    return {
        "path": path,
        "branch": git_text("branch", "--show-current", cwd=path),
        "sha": git_text("rev-parse", "HEAD", cwd=path),
        "dirty": bool(status),
        "porcelain_status": status,
    }


SIGNER = "validation-manager"

# Trust map consulted by the gate. Deliberately not overridable by environment
# variable (tests monkeypatch this attribute).
TRUSTED_KEYS_PATH = Path(__file__).resolve().parents[1] / "keys" / "trusted_attestation_keys.json"

BOOTSTRAP_NOTE = "unverified (no trusted key published yet)"

# Tasks for which bootstrap tolerance (unverified records counted while the
# signer has no trusted key) applies: only the migration task that introduces
# attestations. Every other task fails closed whatever the trust map state.
BOOTSTRAP_TASKS = frozenset({"V3-INFRA-002"})

ATTESTATION_REQUIRED = (
    "attestation required: no trusted validation-manager key published "
    "(operator runbook) or record not verified"
)
HOST_RUN_SIGNER = "host-runner"


def _attest():
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from infrastructure.experiment_v3 import attest

    return attest


def load_trust() -> dict:
    return _attest().load_trusted(TRUSTED_KEYS_PATH)


def enforcing(trusted: dict) -> bool:
    """True once the trust map lists at least one key for the signer."""
    return bool(trusted.get(SIGNER))


# Key-publication mode (resolves the trust bootstrap after the migration merges
# and governs every later change of the trust file, e.g. key rotation).
# The gate runs from the integration checkout, whose trust file is empty until
# a key-publication commit merges; that commit could never pass an empty-trust
# gate. Therefore, while the integration trust map has no validation-manager
# key, a task branch whose ONLY change against the integration base is the
# trust file itself, and whose trust file lists exactly one well-formed key
# for every signer, is evaluated against ITS OWN trust file, with no bootstrap
# tolerance (every record must verify) and with signed evidence required for
# every published key (evaluate). Any other change on the branch, or a
# malformed trust file, disables the mode and the task stays blocked.
KEY_FILE_REL = "infrastructure/keys/trusted_attestation_keys.json"
CONDITION_PATH = Path(__file__).resolve().parents[2] / ".ionmc-condition.json"


def base_ref() -> str:
    """Remote-tracking ref of the integration branch (monkeypatched in tests)."""
    branch = "develop"
    if CONDITION_PATH.exists():
        branch = json.loads(CONDITION_PATH.read_text(encoding="utf-8")).get(
            "integration_branch", branch
        )
    return f"origin/{branch}"


def parse_trust_map(data, *, exactly_one: bool = False) -> dict:
    """Strict parse of a trust file: {signer: {key_id: public_key_hex}}.

    Every signer in SIGNERS must be present with at least one entry (exactly
    one with ``exactly_one``, the key-publication rule: a published key that
    is not evidenced cannot hide behind an evidenced one); each entry must
    carry a key_id equal to the public key's derived id, and the key must be a
    valid curve point that is not of small order. Raises ValueError.
    """
    attest = _attest()
    if not isinstance(data, dict):
        raise ValueError("trust file is not an object")
    out: dict = {}
    for signer in attest.SIGNERS:
        entries = data.get(signer)
        if not isinstance(entries, list) or not entries:
            raise ValueError(f"trust file has no {signer} key")
        if exactly_one and len(entries) != 1:
            raise ValueError(f"trust file must publish exactly one {signer} key")
        out[signer] = {}
        for e in entries:
            if not isinstance(e, dict) or not isinstance(e.get("public_key"), str):
                raise ValueError(f"malformed {signer} entry")
            try:
                pub = bytes.fromhex(e["public_key"])
            except ValueError as exc:
                raise ValueError(f"malformed {signer} public key") from exc
            if len(pub) != 32 or e.get("key_id") != attest.key_id_of(pub):
                raise ValueError(f"{signer} entry key_id does not match public key")
            point = attest._decompress(pub)
            if point is None or attest._small_order(point):
                raise ValueError(f"{signer} public key is not a valid non-small-order point")
            out[signer][e["key_id"]] = e["public_key"]
    return out


TRUST_CHANGE_RULE = (
    "trust-file changes are allowed only in a pure key-publication branch "
    "(the only change against the integration base is "
    f"{KEY_FILE_REL}, exactly one evidenced key per signer)"
)


def trust_file_change(state: dict) -> tuple[list[str], str | None]:
    """(changed paths of the branch vs its integration merge base, error)."""
    path, sha = state["path"], state["sha"]
    try:
        base = git_text("merge-base", base_ref(), sha, cwd=path)
        changed = git_text("diff", "--name-only", base, sha, cwd=path).splitlines()
    except (subprocess.CalledProcessError, OSError) as exc:
        return [], f"cannot determine the branch changes: {exc}"
    return sorted(c for c in changed if c), None


def key_publication_trust(state: dict) -> tuple[dict | None, str]:
    """Trust map of a pure key-publication branch, else (None, reason)."""
    changed, err = trust_file_change(state)
    if err:
        return None, f"key-publication check failed: {err}"
    if changed != [KEY_FILE_REL]:
        return None, f"branch changes {changed!r}, not exactly [{KEY_FILE_REL!r}]"
    try:
        text = git_text("show", f"{state['sha']}:{KEY_FILE_REL}", cwd=state["path"])
        trusted = parse_trust_map(json.loads(text), exactly_one=True)
    except (subprocess.CalledProcessError, ValueError, OSError) as exc:
        return None, f"key-publication check failed: {exc}"
    return trusted, "key-publication branch: trust map taken from the task SHA"


def effective_trust(state: dict) -> tuple[dict, str, str | None]:
    """(trust map, trust_source, problem).

    A branch that does not touch the trust file is evaluated against the
    integration trust file ("integration"). A branch that touches it is
    evaluated against ITS OWN trust file ("key-publication") only if it is a
    pure key-publication branch (first publication or rotation alike: the only
    change is the trust file, exactly one valid key per signer); evaluate()
    then demands signed evidence for every published key. Any other change to
    the trust file yields a ``problem`` (TRUST_CHANGE_RULE) and the gate fails,
    whatever keys the integration file already trusts.
    """
    trusted = load_trust()
    changed, err = trust_file_change(state)
    if err:
        return trusted, "integration", err
    if KEY_FILE_REL not in changed:
        return trusted, "integration", None
    branch_trust, reason = key_publication_trust(state)
    if branch_trust is not None:
        return branch_trust, "key-publication", None
    return trusted, "integration", f"{TRUST_CHANGE_RULE}; {reason}"


def classify(rec: dict, task_id: str, sha: str, trusted: dict) -> tuple[str, str]:
    """Return (attestation_status, reason): verified | unverified | invalid.

    Besides the signature, the signed identity must match the requested task
    and SHA, so a genuine record copied into another task/SHA directory is
    invalid.
    """
    status, reason = _attest().classify_record(rec, trusted, SIGNER)
    if status == "invalid":
        return status, reason
    if rec.get("task_id") != normalize_task_id(task_id) or rec.get("git_sha") != sha:
        return "invalid", (
            f"signed identity mismatch: record is for task {rec.get('task_id')!r} "
            f"sha {rec.get('git_sha')!r}, requested {normalize_task_id(task_id)!r} {sha!r}"
        )
    return status, reason


def records_for(task_id: str, sha: str, trusted: dict | None = None) -> list[dict]:
    """Load records, each annotated with derived attestation_status/_reason.

    The derived fields are added after verification and are not part of the
    signed payload; they are never read from disk (any such field in a file is
    overwritten). ``trusted`` defaults to the integration trust file.
    """
    directory = VALIDATION_ROOT / normalize_task_id(task_id) / sha

    if not directory.exists():
        return []

    if trusted is None:
        trusted = load_trust()
    records = []

    for path in sorted(directory.glob("*.json")):
        try:
            rec = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(rec, dict):
                continue
        except Exception:
            continue
        rec.pop("attestation_status", None)
        rec.pop("attestation_reason", None)
        rec.pop("host_run_attestation", None)
        status, reason = classify(rec, task_id, sha, trusted)
        rec["attestation_status"] = status
        rec["attestation_reason"] = reason
        rec["host_run_attestation"] = check_host_runs(rec, task_id, sha, trusted)
        records.append(rec)

    return records


def bootstrap_allowed(task_id: str, trusted: dict, signer: str = SIGNER) -> bool:
    """Bootstrap tolerance: signer has no trusted key AND task is a migration task."""
    return not trusted.get(signer) and normalize_task_id(task_id) in BOOTSTRAP_TASKS


LEGACY_NO_BINDING = "host run referenced without output_manifest_sha256 binding"
MANIFEST_PATTERN = re.compile(r"[0-9a-f]{64}")
_NO_CLAIM = object()


def host_run_declaration_problems(details) -> list[str]:
    """Reasons the reserved host-run fields of a details dict are malformed.

    Reserved fields: ``host_runs`` (list of ``{"run_id": non-blank str,
    "output_manifest_sha256": 64 hex or null (optional)}``), ``host_run_id``
    (non-blank str) and ``host_run_ids`` (list of non-blank str). Any present
    field with another shape is a problem; a malformed declaration must never
    be silently ignored (fail closed).
    """
    problems: list[str] = []
    if not isinstance(details, dict):
        return problems
    if "host_runs" in details:
        runs = details["host_runs"]
        if not isinstance(runs, list):
            problems.append("details.host_runs is not a list")
        else:
            for i, item in enumerate(runs):
                if not isinstance(item, dict):
                    problems.append(f"details.host_runs[{i}] is not an object")
                    continue
                run_id = item.get("run_id")
                if not isinstance(run_id, str) or not run_id.strip():
                    problems.append(f"details.host_runs[{i}].run_id is not a non-blank string")
                if "output_manifest_sha256" in item:
                    claim = item["output_manifest_sha256"]
                    if not (
                        claim is None
                        or (isinstance(claim, str) and MANIFEST_PATTERN.fullmatch(claim))
                    ):
                        problems.append(
                            f"details.host_runs[{i}].output_manifest_sha256 "
                            "is not 64 hex digits or null"
                        )
    if "host_run_id" in details:
        value = details["host_run_id"]
        if not isinstance(value, str) or not value.strip():
            problems.append("details.host_run_id is not a non-blank string")
    if "host_run_ids" in details:
        many = details["host_run_ids"]
        if not isinstance(many, list) or not all(
            isinstance(v, str) and v.strip() for v in many
        ):
            problems.append("details.host_run_ids is not a list of non-blank strings")
    return problems


def host_run_refs(rec: dict) -> list[tuple[str, object]]:
    """Host runs referenced by a validation record as (run_id, claim).

    Convention: ``details.host_runs = [{"run_id": "RUN-...",
    "output_manifest_sha256": "<64 hex or null>"}]``. ``claim`` is the stated
    digest (str or None) when the key is present, else the ``_NO_CLAIM``
    sentinel. The legacy ``host_run_id`` / ``host_run_ids`` forms never carry a
    claim. Malformed entries are skipped here; callers must consult
    host_run_declaration_problems so that they are rejected, not ignored.
    """
    details = rec.get("details")
    if not isinstance(details, dict):
        return []
    refs: list[tuple[str, object]] = []
    runs = details.get("host_runs")
    if isinstance(runs, list):
        for item in runs:
            if isinstance(item, dict) and isinstance(item.get("run_id"), str):
                run_id = item["run_id"].strip()
                if run_id:
                    refs.append((run_id, item.get("output_manifest_sha256", _NO_CLAIM)))
    value = details.get("host_run_id")
    if isinstance(value, str) and value.strip():
        refs.append((value.strip(), _NO_CLAIM))
    many = details.get("host_run_ids")
    if isinstance(many, list):
        refs.extend((v.strip(), _NO_CLAIM) for v in many if isinstance(v, str) and v.strip())
    return refs


def host_run_ids(rec: dict) -> list[str]:
    """Host run ids referenced by a validation record's details."""
    return [run_id for run_id, _ in host_run_refs(rec)]


def run_success_problem(run: dict) -> str | None:
    """Reason a signed run record cannot support validation, else None."""
    checks = [
        ("succeeded", run.get("succeeded") is True),
        ("exit_code", run.get("exit_code") == 0),
        ("timed_out", run.get("timed_out") is False),
        (
            "committed_state_unchanged",
            "committed_state_unchanged" not in run or run["committed_state_unchanged"] is True,
        ),
        ("dirty_after", "dirty_after" not in run or run["dirty_after"] is False),
        ("output_publication_error", run.get("output_publication_error") is None),
    ]
    bad = [f"{name}={run.get(name)!r}" for name, ok in checks if not ok]
    if bad:
        return "host run did not succeed: " + ", ".join(bad)
    return None


def check_host_runs(rec: dict, task_id: str, sha: str, trusted: dict) -> list[dict]:
    """Verify every signed run record a validation record references.

    Each entry: {run_id, status, reason, ok}. ``ok`` is False for a missing,
    invalid or (when not tolerated) unverified run record, a task/SHA that
    differs from the validation record's, an output manifest digest that
    differs from the validation record's claim (null included; a signed record
    without the field is invalid), a malformed reserved host-run field (one
    "<malformed>" entry per problem, never ignored), a legacy
    reference without a manifest claim (when not tolerated), or a run that did
    not succeed (see run_success_problem).
    """
    from infrastructure.validation import host_run_records as hrr

    tolerated = bootstrap_allowed(task_id, trusted, HOST_RUN_SIGNER)
    details = rec.get("details") if isinstance(rec.get("details"), dict) else {}
    lookup_missing = details.get("host_run_manifest_lookup") == "missing"
    out = [
        {"run_id": "<malformed>", "status": "malformed", "reason": problem, "ok": False}
        for problem in host_run_declaration_problems(details)
    ]
    for run_id, claim in host_run_refs(rec):
        entry = {"run_id": run_id, "status": "missing", "reason": "", "ok": False}
        out.append(entry)
        if claim is not _NO_CLAIM and not (
            claim is None or (isinstance(claim, str) and MANIFEST_PATTERN.fullmatch(claim))
        ):
            entry["reason"] = "output_manifest_sha256 claim is not 64 hex digits or null"
            continue
        try:
            run = hrr.load_signed_run_record(run_id, normalize_task_id(task_id))
        except hrr.RunRecordMissing as exc:
            entry["reason"] = f"host run record missing: {exc}"
            continue
        status, reason = hrr.verify_run_record(
            run,
            rec.get("task_id", ""),
            rec.get("git_sha", ""),
            run_id,
            trusted,
            expected_output_manifest_sha256=claim if isinstance(claim, str) else None,
        )
        entry["status"], entry["reason"] = status, reason
        if status == "invalid":
            continue
        if status != "verified" and not tolerated:
            entry["reason"] = f"host run record {status} (attestation required): {reason}"
            continue
        if "output_manifest_sha256" not in run:
            entry["status"] = "invalid"
            entry["reason"] = (
                "signed run record has no output_manifest_sha256 field "
                "(absence is not a null manifest)"
            )
            continue
        if claim is None and run["output_manifest_sha256"] is not None:
            entry["status"] = "invalid"
            entry["reason"] = (
                "output_manifest_sha256 mismatch: claimed null, signed run record has "
                f"{run['output_manifest_sha256']!r}"
            )
            continue
        if claim is _NO_CLAIM and (not tolerated or lookup_missing):
            entry["reason"] = LEGACY_NO_BINDING + (
                " (run record was missing at record time)" if lookup_missing else ""
            )
            continue
        problem = run_success_problem(run)
        if problem:
            entry["reason"] = problem
            continue
        entry["ok"] = True
    return out


def evaluate(
    records: list[dict],
    dirty: bool,
    task_id: str,
    trusted: dict | None = None,
    trust_source: str = "integration",
    trust_problem: str | None = None,
) -> dict:
    """Gate decision over annotated records for ``task_id``.

    Enforcement rule (fail closed):
    - Only records with attestation_status "verified" count as passing, except
      in bootstrap. Unverified and invalid records are rejected with a reason.
    - Bootstrap (unverified records counted, noted as "unverified (no trusted
      key published yet)") applies ONLY when task_id is in BOOTSTRAP_TASKS
      (the migration task V3-INFRA-002) AND the trust map has no key for
      "validation-manager". For every other task unverified records are
      rejected whatever the trust map state, with ATTESTATION_REQUIRED; all
      other tasks' gates are therefore blocked until the operator publishes the
      keys.
    - "invalid" records (present but wrong, or signed identity not matching the
      requested task/SHA) are always rejected, even in bootstrap.
    - A passed record that references host runs (details.host_runs, each with
      run_id and output_manifest_sha256) is rejected unless every referenced
      signed run record is present, bound to the same task/SHA and manifest
      digest, verified (or bootstrap-tolerated, see check_host_runs) and
      succeeded (succeeded true, exit_code 0, not timed out, state unchanged,
      no publication error). Legacy host_run_id(s) references carry no
      manifest binding and are rejected outside bootstrap.
    - Every record with status "failed" counts as failing regardless of its
      attestation (a forged failure can only block, never release, the gate).
    - A branch that changes the trust file in any way other than a pure
      key-publication branch fails (trust_problem, TRUST_CHANGE_RULE), also
      when the integration file already trusts a key.
    - Key-publication mode (trust_source "key-publication", see
      key_publication_trust): the trust map is the task branch's own file
      (exactly one key per signer); enforcement is then "enforcing" with no
      bootstrap tolerance, and the gate additionally requires signed evidence
      for every published key: a verified passing validation record
      (validation-manager) referencing a verified successful host run
      (host-runner), both bound to this task and SHA.
    """
    if trusted is None:
        trusted = load_trust()
    enforce = enforcing(trusted)
    tolerate = bootstrap_allowed(task_id, trusted)
    passing, rejected = [], []

    for r in records:
        att = r.get("attestation_status", "unverified")
        runs = r.get("host_run_attestation") or []
        bad_run = next((x for x in runs if not x.get("ok")), None)
        if att == "invalid" or (att != "verified" and not tolerate):
            reason = r.get("attestation_reason")
            if att == "unverified" and not enforce:
                reason = ATTESTATION_REQUIRED
            rejected.append(
                {
                    "validation_id": r.get("validation_id"),
                    "attestation_status": att,
                    "reason": reason,
                }
            )
        elif r.get("status") == "passed" and bad_run is not None:
            rejected.append(
                {
                    "validation_id": r.get("validation_id"),
                    "attestation_status": att,
                    "reason": f"host run {bad_run['run_id']}: {bad_run['reason']}",
                }
            )
        elif r.get("status") == "passed":
            passing.append(r)

    failing = [r for r in records if r.get("status") == "failed"]
    unverified = [r for r in records if r.get("attestation_status") == "unverified"]

    summary = {
        "enforcement": "enforcing" if enforce else ("bootstrap" if tolerate else "blocked"),
        "trust_source": trust_source,
        "verified_records": sum(
            1 for r in records if r.get("attestation_status") == "verified"
        ),
        "unverified_records": len(unverified),
        "invalid_records": sum(
            1 for r in records if r.get("attestation_status") == "invalid"
        ),
        "rejected_records": rejected,
        "passing_records": len(passing),
        "failing_records": len(failing),
        "gate_passed": not dirty and len(passing) > 0 and len(failing) == 0,
    }
    if trust_problem:
        summary["trust_problem"] = trust_problem
        summary["gate_passed"] = False
    if trust_source == "key-publication":
        # Every published key must be evidenced by a protected signed record for
        # THIS task and SHA: the validation-manager key by a verified passing
        # validation record (already required above), the host-runner key by a
        # verified, successful host run referenced from a passing record.
        evidence = {
            SIGNER: bool(passing),
            HOST_RUN_SIGNER: any(
                x.get("status") == "verified" and x.get("ok")
                for r in passing
                for x in (r.get("host_run_attestation") or [])
            ),
        }
        summary["key_publication_evidence"] = evidence
        if not all(evidence.values()):
            missing = [k for k, v in evidence.items() if not v]
            summary["key_publication_problem"] = (
                "key-publication requires signed evidence for every published "
                f"key; none verified for {missing}"
            )
            summary["gate_passed"] = False
    if tolerate and unverified:
        summary["attestation_note"] = (
            f"{len(unverified)} record(s) {BOOTSTRAP_NOTE}; tolerated only for "
            f"bootstrap task(s) {sorted(BOOTSTRAP_TASKS)} until the operator "
            "publishes the validation-manager key"
        )
    if not enforce and not tolerate:
        summary["attestation_required"] = ATTESTATION_REQUIRED
    return summary


def fill_host_run_manifests(details, task_id: str):
    """Fill output_manifest_sha256 of details.host_runs entries at record time.

    The digest comes from the signed run record under the protected host-runs
    root, so the caller cannot mis-state it. If a run record is missing the
    entry is left unbound and details.host_run_manifest_lookup = "missing" is
    set; the gate rejects such a record later.
    """
    if not isinstance(details, dict) or not isinstance(details.get("host_runs"), list):
        return details
    from infrastructure.validation import host_run_records as hrr

    for item in details["host_runs"]:
        if not isinstance(item, dict) or "output_manifest_sha256" in item:
            continue
        run_id = item.get("run_id")
        try:
            run = hrr.load_signed_run_record(run_id, task_id)
        except hrr.RunRecordMissing:
            details["host_run_manifest_lookup"] = "missing"
            continue
        item["output_manifest_sha256"] = run.get("output_manifest_sha256")
    return details


def record(
    task_id: str,
    suite: str,
    status: str,
    summary: str,
    details_json: str,
):
    state = worktree_state(task_id)

    if state["dirty"]:
        raise RuntimeError(
            "Refusing to record validation for a dirty worktree. "
            "Commit the exact state being validated first."
        )

    try:
        details = json.loads(details_json)
    except json.JSONDecodeError as exc:
        raise ValueError("--details-json must contain valid JSON") from exc

    problems = host_run_declaration_problems(details)
    if problems:
        raise ValueError("malformed host-run declaration: " + "; ".join(problems))
    details = fill_host_run_manifests(details, normalize_task_id(task_id))

    now = datetime.now(timezone.utc)

    validation_id = (
        f"VAL-{now.strftime('%Y%m%d-%H%M%S')}-"
        f"{uuid.uuid4().hex[:6].upper()}"
    )

    record = {
        "schema_version": 1,
        "validation_id": validation_id,
        "task_id": normalize_task_id(task_id),
        "recorded_at_utc": now.isoformat(),
        "suite": suite,
        "status": status,
        "summary": summary,
        "details": details,
        "git_branch": state["branch"],
        "git_sha": state["sha"],
        "worktree": str(state["path"]),
    }

    # V3-INFRA-002: Ed25519 attestation by the service-held key (additive field).
    record = _attest().sign_with_service_key(record, SIGNER)

    directory = VALIDATION_ROOT / normalize_task_id(task_id) / state["sha"]
    directory.mkdir(parents=True, exist_ok=True)

    path = directory / f"{validation_id}.json"

    path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "status": "recorded",
                "validation": record,
                "record_file": str(path),
            },
            indent=2,
        )
    )


def inspect(task_id: str):
    state = worktree_state(task_id)
    trusted, source, problem = effective_trust(state)
    records = records_for(task_id, state["sha"], trusted)
    verdict = evaluate(records, state["dirty"], task_id, trusted, source, problem)

    result = {
        "status": "present" if records else "missing",
        "task_id": normalize_task_id(task_id),
        "git_branch": state["branch"],
        "git_sha": state["sha"],
        "worktree_dirty": state["dirty"],
        "records": records,
        **verdict,
    }

    print(json.dumps(result, indent=2))


def check(task_id: str):
    state = worktree_state(task_id)
    trusted, source, problem = effective_trust(state)
    records = records_for(task_id, state["sha"], trusted)
    verdict = evaluate(records, state["dirty"], task_id, trusted, source, problem)

    result = {
        "status": "passed" if verdict["gate_passed"] else "failed",
        "task_id": normalize_task_id(task_id),
        "git_sha": state["sha"],
        "worktree_dirty": state["dirty"],
        **{k: v for k, v in verdict.items() if k != "gate_passed"},
    }

    print(json.dumps(result, indent=2))

    if not verdict["gate_passed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    p_record = sub.add_parser("record")
    p_record.add_argument("--task-id", required=True)
    p_record.add_argument("--suite", required=True)
    p_record.add_argument(
        "--status",
        choices=["passed", "failed"],
        required=True,
    )
    p_record.add_argument("--summary", required=True)
    p_record.add_argument("--details-json", default="{}")

    p_inspect = sub.add_parser("inspect")
    p_inspect.add_argument("--task-id", required=True)

    p_check = sub.add_parser("check")
    p_check.add_argument("--task-id", required=True)

    args = parser.parse_args()

    if args.command == "record":
        record(
            args.task_id,
            args.suite,
            args.status,
            args.summary,
            args.details_json,
        )

    elif args.command == "inspect":
        inspect(args.task_id)

    elif args.command == "check":
        check(args.task_id)


if __name__ == "__main__":
    main()
