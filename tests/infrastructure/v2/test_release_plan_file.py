"""The committed release plan must satisfy the qualification schema."""

import json
from pathlib import Path

from infrastructure.experiment_v2.release import validate_plan

ROOT = Path(__file__).resolve().parents[3]


def test_committed_release_plan_is_frozen_and_complete():
    plan = json.loads((ROOT / "validation/release-plan.json").read_text())
    ids = {
        r["id"]
        for r in json.loads(
            (ROOT / "experiment/v2/requirements-index.json").read_text()
        )["requirements"]
    }
    suites = validate_plan(plan, ids)
    assert len(suites) >= 18
    assert plan["status"] == "frozen"
    for suite in suites:
        assert "DRAFT" not in suite["acceptance"]
    for workload, targets in plan["performance_targets"].items():
        assert targets["histories_per_second"]["direction"] == "min", workload
        assert targets["wall_seconds"]["direction"] == "max", workload


def test_ledger_covers_every_requirement_with_a_plan_suite():
    plan = json.loads((ROOT / "validation/release-plan.json").read_text())
    ledger = json.loads((ROOT / "validation/requirement-ledger.json").read_text())
    suite_ids = {s["id"] for s in plan["suites"]}
    index_ids = {
        r["id"]
        for r in json.loads(
            (ROOT / "experiment/v2/requirements-index.json").read_text()
        )["requirements"]
    }
    ledger_ids = {r["id"] for r in ledger["requirements"]}
    assert ledger_ids == index_ids
    for item in ledger["requirements"]:
        assert item["suites"] and set(item["suites"]) <= suite_ids, item["id"]
        assert item["status"] in ledger["status_vocabulary"]
