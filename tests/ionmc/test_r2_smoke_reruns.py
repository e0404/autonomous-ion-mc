"""V3-008 revision 2: the committed run-validity summaries of the signed smoke runs (the 14
re-runs ``validation/results/transport/v8r2-topas/smoke-reruns.json`` and the pre-registered
confirmation set ``smoke-confirmation.json``, written by
``validation/scripts/ion_em/r2_smoke_reruns.py`` through ``r2_evaluate.py check-run``) and the
committed signed records under ``smoke-runs/<run_id>/`` (pure-stdlib attestation check)."""

from __future__ import annotations

import hashlib
import json
import re
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
        _verify_record(r["run_id"], trusted)


def _verify_record(run_id: str, trusted: dict[str, Any]) -> None:
    """``smoke-runs/<run_id>/record.json`` verifies against the trust file and is bound to its
    committed artifact listing."""
    rdir = SMOKE_RUNS / run_id
    blob = (rdir / "record.json").read_bytes()
    record = json.loads(blob)
    status, reason = attest.classify_record(record, trusted, "reference-runner")
    assert status == "verified", (run_id, status, reason)
    assert record["run_id"] == run_id and record["schema_version"] == 3
    pages = json.loads((rdir / "artifacts.json").read_text(encoding="utf-8"))
    listed = {e["path"]: e["sha256"] for p in pages for e in p["files"]}
    assert all(p["run_id"] == run_id for p in pages)
    assert listed["record.json"] == hashlib.sha256(blob).hexdigest()
    assert listed["stdout.txt"] == record["stdout_sha256"]
    assert listed["stderr.txt"] == record["stderr_sha256"]
    assert listed["request.json"] == record["request_sha256"]


def _reruns_module() -> Any:
    import importlib.util
    import sys

    path = ROOT / "validation" / "scripts" / "ion_em" / "r2_smoke_reruns.py"
    spec = importlib.util.spec_from_file_location("r2_smoke_reruns_under_test", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _rows(m: Any) -> list[dict[str, Any]]:
    return [{"run_id": f"REF-{i:02d}", "case": b, "seed": seed, "verdict": "VALID"}
            for i, (b, seed) in enumerate(sorted(m.CONFIRMATION_SET))]  # fmt: skip


def test_confirmation_set_is_the_preregistered_one() -> None:
    m = _reruns_module()
    seeds = sorted(seed for _, seed in m.CONFIRMATION_SET)
    assert seeds == list(range(20600424, 20600438)) and len(m.CONFIRMATION_SET) == 14
    lp = {b for b, seed in m.CONFIRMATION_SET if seed >= 20600428}
    assert len(lp) == 10 and all("-lp-" in b for b in lp)


def test_confirmation_exact_set_is_complete() -> None:
    m = _reruns_module()
    rows = _rows(m)
    assert m.confirmation(rows, True)["confirmation_set_complete"] is True
    assert m.confirmation(rows, False)["confirmation_set_complete"] is False  # not requested


def test_confirmation_refuses_missing_duplicate_unexpected_invalid_and_arbitrary() -> None:
    m = _reruns_module()
    rows = _rows(m)
    res = m.confirmation(rows[1:], True)  # one missing
    assert res["confirmation_set_complete"] is False
    assert "missing" in res["confirmation_reason"] and rows[0]["case"] in res["confirmation_reason"]
    res = m.confirmation([*rows, dict(rows[3], run_id="REF-dup")], True)  # duplicate
    assert res["confirmation_set_complete"] is False and "duplicate" in res["confirmation_reason"]
    extra = {"run_id": "REF-x", "case": "he-water-250mevu-r2-idd-em-timing-smoke",
             "seed": 20600411, "verdict": "VALID"}  # fmt: skip
    res = m.confirmation([*rows, extra], True)  # unexpected (valid) run
    assert res["confirmation_set_complete"] is False and "unexpected" in res["confirmation_reason"]
    assert "20600411" in res["confirmation_reason"]
    wrong_seed = [dict(rows[0], seed=20600499), *rows[1:]]  # right bundle, wrong seed
    res = m.confirmation(wrong_seed, True)
    assert res["confirmation_set_complete"] is False
    assert "missing" in res["confirmation_reason"] and "unexpected" in res["confirmation_reason"]
    res = m.confirmation([dict(rows[0], verdict="RUN-VALIDITY FAIL"), *rows[1:]], True)
    assert res["confirmation_set_complete"] is False and "not VALID" in res["confirmation_reason"]
    arbitrary = [{"run_id": f"REF-{i}", "case": "c", "seed": i, "verdict": "VALID"}
                 for i in range(3)]  # fmt: skip
    res = m.confirmation(arbitrary, True)
    assert res["confirmation_set_complete"] is False and res["confirmation_reason"].startswith(
        "refused"
    )
    assert m.confirmation([], True)["confirmation_set_complete"] is False


CONFIRMATION = TOPAS / "smoke-confirmation.json"


def _confirmation_doc() -> dict[str, Any]:
    if not CONFIRMATION.is_file():
        pytest.skip(f"{CONFIRMATION.relative_to(ROOT)} is not committed")
    return json.loads(CONFIRMATION.read_text(encoding="utf-8"))


def test_committed_confirmation_summary_is_the_complete_valid_set() -> None:
    doc = _confirmation_doc()
    m = _reruns_module()
    pairs = [(r["case"], r["seed"]) for r in doc["runs"]]
    assert doc["n_runs"] == len(doc["runs"]) == 14
    assert len(set(pairs)) == 14 and set(pairs) == m.CONFIRMATION_SET
    assert len({r["run_id"] for r in doc["runs"]}) == 14
    assert all(r["verdict"] == "VALID" for r in doc["runs"]), [
        (r["run_id"], r["verdict"]) for r in doc["runs"] if r["verdict"] != "VALID"
    ]
    assert doc["all_valid"] is True and doc["confirmation_set_complete"] is True
    assert re.fullmatch(r"[0-9a-f]{40}", doc["evaluator_head"])
    assert m.confirmation(doc["runs"], True)["confirmation_set_complete"] is True


def test_committed_confirmation_records_verify_against_the_trust_file() -> None:
    doc = _confirmation_doc()
    trusted = attest.load_trusted(TRUST)
    for r in doc["runs"]:
        _verify_record(r["run_id"], trusted)
        record = json.loads((SMOKE_RUNS / r["run_id"] / "record.json").read_bytes())
        assert record["code_sha"] == r["code_sha"] and record["started_at"] == r["started_at"]


def test_several_ledgers_are_combined_and_duplicates_refused(tmp_path: Path) -> None:
    m = _reruns_module()
    a, b = tmp_path / "a.jsonl", tmp_path / "b.jsonl"
    a.write_text(json.dumps({"run_id": "REF-1"}) + "\n" + json.dumps({"run_id": "REF-2"}) + "\n")
    b.write_text(json.dumps({"run_id": "REF-3"}) + "\n")
    entries, prov = m.read_ledgers([a, b])
    assert [e["run_id"] for e in entries] == ["REF-1", "REF-2", "REF-3"]
    assert [p["n_runs"] for p in prov] == [2, 1] and all(len(p["sha256"]) == 64 for p in prov)
    b.write_text(json.dumps({"run_id": "REF-2"}) + "\n")
    with pytest.raises(m.ev.io.RecordError, match="REF-2"):
        m.read_ledgers([a, b])
