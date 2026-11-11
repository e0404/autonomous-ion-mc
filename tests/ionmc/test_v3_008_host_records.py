"""V3-008 revision 2: inventory of the committed signed host-run records of the sample chunks.

``validation/results/transport/v8r2-samples/host-records/<run_id>.json`` holds the complete signed
host-runner record of every sample chunk run of the two seals (``sample-seal.json``,
``sample-seal-heldout-r2stop3.json``); the evaluator (``r2_io.verify_host_records``) binds each
graded chunk to its record. Cache-free: reads the committed JSON and trust file only (no npz, no
published run directory). The records were requested from the operator
(IR-20261111-101840-EDCAE6); until they are committed the completeness test is a strict xfail.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SAMPLES = REPO / "validation/results/transport/v8r2-samples"
HOST_RECORDS = SAMPLES / "host-records"
SEALS = (SAMPLES / "sample-seal.json", SAMPLES / "sample-seal-heldout-r2stop3.json")
TRUST = REPO / "infrastructure/keys/trusted_attestation_keys.json"


def _attest():
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from infrastructure.experiment_v3 import attest

    return attest


def _runs() -> dict[str, dict[str, str]]:
    """run_id -> {implementation_sha, output_manifest_sha256} of both seals' ``host_runs``."""
    runs: dict[str, dict[str, str]] = {}
    for path in SEALS:
        seal = json.loads(path.read_text(encoding="utf-8"))
        for h in seal["host_runs"]:
            assert h["run_id"] not in runs, f"{h['run_id']} listed twice"
            runs[h["run_id"]] = {
                "implementation_sha": seal["implementation_sha"],
                "output_manifest_sha256": h["output_manifest_sha256"],
            }
    return runs


def test_seal_host_runs_count() -> None:
    assert len(_runs()) == 76 + 13


def test_committed_host_records_verify_and_bind() -> None:
    """Every file present must be ``<run_id>.json`` of a sealed run, signed by the trusted
    host-runner key, and carry that run id, implementation SHA and output-manifest digest."""
    runs = _runs()
    files = sorted(HOST_RECORDS.iterdir()) if HOST_RECORDS.is_dir() else []
    trusted = _attest().load_trusted(TRUST)
    for path in files:
        assert path.is_file() and not path.is_symlink(), path
        assert path.suffix == ".json" and path.stem in runs, f"{path.name}: not a sealed run"
        rec = json.loads(path.read_text(encoding="utf-8"))
        status, why = _attest().classify_record(rec, trusted, "host-runner")
        assert status == "verified", f"{path.name}: {status} ({why})"
        want = runs[path.stem]
        assert rec["run_id"] == path.stem
        assert rec["sha"] == want["implementation_sha"]
        assert rec["output_manifest_sha256"] == want["output_manifest_sha256"]


@pytest.mark.xfail(
    strict=True,
    reason="host records pending IR-20261111-101840-EDCAE6; remove this marker in the commit that "
    "adds them",
)
def test_every_sealed_run_has_a_host_record() -> None:
    present = {p.stem for p in HOST_RECORDS.glob("*.json")} if HOST_RECORDS.is_dir() else set()
    missing = sorted(set(_runs()) - present)
    assert not missing, f"{len(missing)} sealed runs without a committed host record"
