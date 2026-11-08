"""Production consumer of signed host-runner (RUN-...) records.

The host runner writes ``<host-runs root>/<task>/<RUN-id>/record.json`` with an
Ed25519 attestation (signer "host-runner"). This module loads such a record by
run id and verifies it *bound to the identity the caller expects*: task id,
git SHA and run id (and optionally the output-manifest digest). A genuine,
validly signed record replayed for another task, SHA or run id is "invalid".

Status values follow infrastructure/experiment_v3/attest.py: verified,
unverified (no attestation or key not trusted), invalid (anything wrong).
Whether "unverified" is tolerated (bootstrap) is the caller's decision; see
infrastructure/validation/local_validation.py.
"""

import json
import re
import sys
from pathlib import Path

SIGNER = "host-runner"
RUN_ID_PATTERN = re.compile(r"RUN-[A-Za-z0-9][A-Za-z0-9._-]{0,127}")


class RunRecordMissing(Exception):
    """No (unambiguous) signed run record exists for the run id."""


def _attest():
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from infrastructure.experiment_v3 import attest

    return attest


def host_runs_root() -> Path:
    """Protected host-runs root (the same directory the host runner writes)."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from infrastructure.experiment_v3.paths import SHARE_ROOT

    return SHARE_ROOT / "host-runs"


HOST_RUNS_ROOT = None  # tests monkeypatch; None means host_runs_root()


def load_signed_run_record(host_run_id: str, task_id: str | None = None) -> dict:
    """Load ``record.json`` of run ``host_run_id`` (raises RunRecordMissing).

    With ``task_id`` only that task's directory is consulted; without it the
    run id must be found under exactly one task directory.
    """
    if not isinstance(host_run_id, str) or not RUN_ID_PATTERN.fullmatch(host_run_id):
        raise RunRecordMissing(f"invalid host run id: {host_run_id!r}")
    root = HOST_RUNS_ROOT if HOST_RUNS_ROOT is not None else host_runs_root()
    if task_id is not None:
        candidates = [Path(root) / task_id.strip().upper() / host_run_id / "record.json"]
    else:
        candidates = sorted(Path(root).glob(f"*/{host_run_id}/record.json"))
    found = [path for path in candidates if path.is_file()]
    if len(found) != 1:
        raise RunRecordMissing(
            f"host run record {host_run_id}: {'not found' if not found else 'ambiguous'}"
        )
    try:
        record = json.loads(found[0].read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RunRecordMissing(f"host run record {host_run_id} unreadable: {exc}") from exc
    if not isinstance(record, dict):
        raise RunRecordMissing(f"host run record {host_run_id} is not an object")
    return record


def verify_run_record(
    record: dict,
    expected_task_id: str,
    expected_sha: str,
    expected_run_id: str,
    trust: dict,
    expected_output_manifest_sha256: str | None = None,
) -> tuple[str, str]:
    """Return (status, reason): attestation plus identity binding.

    Signature/trust problems keep their attest classification. Any mismatch of
    the signed task_id, sha (run-record field name), run_id or, when given,
    output_manifest_sha256 makes the record "invalid" even when the signature
    is genuine (replay protection).
    """
    status, reason = _attest().classify_record(record, trust, SIGNER)
    if status == "invalid":
        return status, reason
    expected = {
        "task_id": expected_task_id.strip().upper(),
        "sha": expected_sha,
        "run_id": expected_run_id,
    }
    if expected_output_manifest_sha256 is not None:
        expected["output_manifest_sha256"] = expected_output_manifest_sha256
    for field, want in expected.items():
        if record.get(field) != want:
            return "invalid", (
                f"signed identity mismatch: {field} is {record.get(field)!r}, "
                f"expected {want!r}"
            )
    return status, reason
