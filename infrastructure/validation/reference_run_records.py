"""Production consumer of signed reference-runner (REF-...) records.

The reference runner (infrastructure/experiment_v3/reference.py) writes
``<state>/references/<REF-id>/record.json`` (schema_version 3) with an Ed25519
attestation (signer "reference-runner"). This module loads such a record by
run id and verifies it *bound to the identity the caller expects*: run id,
exact code SHA and configuration id. A genuine, validly signed record replayed
for another run, SHA or configuration is "invalid". The signed ``started_at``
(the controlled start time of the run, V3-INFRA-004) is required: a record
without a well-formed ``started_at`` not later than ``finished_at`` is
"invalid".

Status values follow infrastructure/experiment_v3/attest.py: verified,
unverified (no attestation, or the reference-runner key is not in the trust
map, e.g. before its key is published), invalid (anything wrong). Whether
"unverified" is tolerated is the caller's decision; chronology gates should
fail closed on anything but "verified".
"""

import json
import re
import sys
from datetime import datetime
from pathlib import Path

SIGNER = "reference-runner"
SCHEMA_VERSION = 3
RUN_ID_PATTERN = re.compile(r"REF-([0-9a-f]{20})-[0-9a-f]{8}")
CONFIGURATION_ID_PATTERN = re.compile(r"[0-9a-f]{64}")


class ReferenceRecordMissing(Exception):
    """No readable signed record exists for the reference run id."""


def _attest():
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from infrastructure.experiment_v3 import attest

    return attest


REFERENCE_STATE = None  # tests monkeypatch; None means common.STATE


def references_root(state=None) -> Path:
    """``<state>/references``; default state as the reference runner uses it."""
    if state is None:
        state = REFERENCE_STATE
    if state is None:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        from infrastructure.experiment_v3.common import STATE

        state = STATE
    return Path(state) / "references"


def load_signed_reference_record(run_id: str, state=None) -> dict:
    """Load ``record.json`` of reference run ``run_id`` (raises ReferenceRecordMissing)."""
    if not isinstance(run_id, str) or not RUN_ID_PATTERN.fullmatch(run_id):
        raise ReferenceRecordMissing(f"invalid reference run id: {run_id!r}")
    path = references_root(state) / run_id / "record.json"
    if not path.is_file():
        raise ReferenceRecordMissing(f"reference run record {run_id}: not found")
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ReferenceRecordMissing(f"reference run record {run_id} unreadable: {exc}") from exc
    if not isinstance(record, dict):
        raise ReferenceRecordMissing(f"reference run record {run_id} is not an object")
    return record


def _utc_time(value):
    """Timezone-aware datetime from an ISO-8601 string, else None."""
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return None
    return parsed if parsed.tzinfo is not None else None


def verify_reference_record(
    record: dict,
    expected_run_id: str,
    expected_code_sha: str,
    expected_configuration_id: str,
    trust: dict,
) -> tuple[str, str]:
    """Return (status, reason): attestation plus identity and start-time binding.

    Signature/trust problems keep their attest classification. Any mismatch of
    the signed run_id, code_sha or configuration_id (replay protection), a
    run_id not derived from the configuration id, a schema_version other than
    3, or a missing/malformed ``started_at`` (or one later than
    ``finished_at``) makes the record "invalid" even when the signature is
    genuine; otherwise the attestation status (verified/unverified) is returned.
    """
    if not isinstance(record, dict):
        return "invalid", "reference run record is not an object"
    status, reason = _attest().classify_record(record, trust, SIGNER)
    if status == "invalid":
        return status, reason
    if record.get("schema_version") != SCHEMA_VERSION:
        return "invalid", (
            f"schema_version is {record.get('schema_version')!r}, expected {SCHEMA_VERSION}"
        )
    expected = {
        "run_id": expected_run_id,
        "code_sha": expected_code_sha,
        "configuration_id": expected_configuration_id,
    }
    for field, want in expected.items():
        if record.get(field) != want:
            return "invalid", (
                f"signed identity mismatch: {field} is {record.get(field)!r}, expected {want!r}"
            )
    match = RUN_ID_PATTERN.fullmatch(str(record["run_id"]))
    config = str(record["configuration_id"])
    if (
        match is None
        or not CONFIGURATION_ID_PATTERN.fullmatch(config)
        or match.group(1) != config[:20]
    ):
        return "invalid", "run_id is not derived from configuration_id"
    started = _utc_time(record.get("started_at"))
    if started is None:
        return "invalid", (
            f"started_at missing or not a timezone-aware ISO-8601 time: "
            f"{record.get('started_at')!r}"
        )
    finished = _utc_time(record.get("finished_at"))
    if finished is None or finished < started:
        return "invalid", (
            f"finished_at {record.get('finished_at')!r} is missing, malformed or "
            f"earlier than started_at {record.get('started_at')!r}"
        )
    return status, reason
