"""V3-008 revision 2: the committed run-validity summary of the 14 signed smoke re-runs
(``validation/results/transport/v8r2-topas/smoke-reruns.json``, written by
``validation/scripts/ion_em/r2_smoke_reruns.py`` through ``r2_evaluate.py check-run``) and the
committed signed records under ``smoke-runs/<run_id>/`` (pure-stdlib attestation check)."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from infrastructure.experiment_v3 import attest

ROOT = Path(__file__).resolve().parents[2]
TOPAS = ROOT / "validation" / "results" / "transport" / "v8r2-topas"
SUMMARY = TOPAS / "smoke-reruns.json"
SMOKE_RUNS = TOPAS / "smoke-runs"
TRUST = ROOT / "infrastructure" / "keys" / "trusted_attestation_keys.json"


def _summary() -> dict[str, Any]:
    if not SUMMARY.is_file():
        pytest.skip(f"{SUMMARY.relative_to(ROOT)} is not committed (smoke re-runs not checked yet)")
    return json.loads(SUMMARY.read_text(encoding="utf-8"))


def test_summary_covers_the_14_preregistered_reruns() -> None:
    doc = _summary()
    seeds = {20600411, 20600412, 20600415, 20600416, 20600419, 20600420,  # EM timing
             20600413, 20600417, 20600421,  # IDD-FULL timing
             20600424, 20600425, 20600426,  # R2-SIZE-1 LAT-FULL
             20600423, 20600427}  # step-length  # fmt: skip
    assert doc["n_runs"] == len(doc["runs"]) == 14
    assert {r["seed"] for r in doc["runs"]} == seeds
    assert len({r["run_id"] for r in doc["runs"]}) == 14
    assert doc["trust_file"]["sha256"] and len(doc["seal_review"]) == 1


def test_every_smoke_rerun_verdict_is_recorded_consistently() -> None:
    """The summary is truthful: ``all_valid`` equals "every verdict VALID"; a non-VALID run names
    its failing checks (the STOP of the 2026-11-09 amendment is recorded, never hidden). The
    confirmation requirement itself (every run of the confirmation set VALID) is asserted only
    once the summary declares ``confirmation_set_complete`` (after the R2-SMOKE-LP runs)."""
    doc = _summary()
    verdicts = {r["run_id"]: r["verdict"] for r in doc["runs"]}
    assert doc["all_valid"] is all(v == "VALID" for v in verdicts.values())
    for r in doc["runs"]:
        if r["verdict"] != "VALID":
            assert r["checks"].get("run_checks_problems") or r["verdict"] == "ERROR", r["run_id"]
    if doc.get("confirmation_set_complete"):
        bad = {k: v for k, v in verdicts.items() if v != "VALID"}
        assert not bad, f"confirmation set not VALID: {bad}"
        assert doc["all_valid"] is True


def test_committed_smoke_records_verify_against_the_trust_file() -> None:
    doc = _summary()
    trusted = attest.load_trusted(TRUST)
    for r in doc["runs"]:
        rdir = SMOKE_RUNS / r["run_id"]
        blob = (rdir / "record.json").read_bytes()
        record = json.loads(blob)
        status, reason = attest.classify_record(record, trusted, "reference-runner")
        assert status == "verified", (r["run_id"], status, reason)
        assert record["run_id"] == r["run_id"] and record["schema_version"] == 3
        pages = json.loads((rdir / "artifacts.json").read_text(encoding="utf-8"))
        listed = {e["path"]: e["sha256"] for p in pages for e in p["files"]}
        assert all(p["run_id"] == r["run_id"] for p in pages)
        assert listed["record.json"] == hashlib.sha256(blob).hexdigest()
        assert listed["stdout.txt"] == record["stdout_sha256"]
        assert listed["stderr.txt"] == record["stderr_sha256"]
        assert listed["request.json"] == record["request_sha256"]
