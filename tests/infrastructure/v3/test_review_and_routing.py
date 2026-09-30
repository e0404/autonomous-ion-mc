import pytest

from infrastructure.experiment_v3 import delegation, review
from infrastructure.experiment_v3.common import file_hash, write_json


def test_model_defaults_and_escalation():
    for role, model in [
        ("general-purpose", "sonnet"),
        ("physics-researcher", "opus"),
        ("Plan", "opus"),
        ("implementation-worker", "sonnet"),
    ]:
        result = delegation.route(
            {
                "tool_name": "Agent",
                "tool_input": {"subagent_type": role, "prompt": "Do the task"},
            }
        )
        assert result["hookSpecificOutput"]["updatedInput"]["model"] == model
        assert "permissionDecision" not in result["hookSpecificOutput"]
    request = {
        "tool_name": "Agent",
        "tool_input": {"model": "fable", "prompt": "Do the task"},
    }
    assert (
        delegation.route(request)["hookSpecificOutput"]["permissionDecision"] == "deny"
    )
    request["tool_input"]["prompt"] += (
        "\nTOP_TIER_JUSTIFICATION: Two Opus attempts found contradictory "
        "derivations; resolve the specific disputed boundary condition."
    )
    assert delegation.route(request)["_routing"]["top_tier_justification"]
    request["tool_input"]["model"] = "unknown-model"
    assert (
        delegation.route(request)["hookSpecificOutput"]["permissionDecision"] == "deny"
    )


def test_review_gate_missing_stale_failed_and_tampered(tmp_path, monkeypatch):
    monkeypatch.setattr(review, "exact_state", lambda _: "head")
    with pytest.raises(ValueError, match="required"):
        review.require_review("V3-001", tmp_path, "base", state=tmp_path)
    dest = tmp_path / "reviews/REVIEW-test"
    report = {
        "head_sha": "head",
        "verdict": "pass",
        "summary": "Inspected actual diff.",
        "findings": [],
    }
    write_json(dest / "report.json", report)
    job = {
        "review_id": "REVIEW-test",
        "task_id": "V3-001",
        "head_sha": "head",
        "base_sha": "base",
        "status": "passed",
        "started_at": "2026-09-30",
        "report_sha256": file_hash(dest / "report.json"),
    }
    write_json(dest / "job.json", job)
    assert (
        review.require_review("V3-001", tmp_path, "base", state=tmp_path)
        == "REVIEW-test"
    )
    with pytest.raises(ValueError):
        review.require_review("V3-001", tmp_path, "new-base", state=tmp_path)
    monkeypatch.setattr(review, "exact_state", lambda _: "new-head")
    with pytest.raises(ValueError):
        review.require_review("V3-001", tmp_path, "base", state=tmp_path)
    monkeypatch.setattr(review, "exact_state", lambda _: "head")
    job["status"] = "failed"
    write_json(dest / "job.json", job)
    with pytest.raises(ValueError):
        review.require_review("V3-001", tmp_path, "base", state=tmp_path)
    job["status"] = "passed"
    write_json(dest / "job.json", job)
    report["summary"] = "Tampered"
    write_json(dest / "report.json", report)
    with pytest.raises(ValueError, match="hash"):
        review.require_review("V3-001", tmp_path, "base", state=tmp_path)


def test_review_rejects_unresolved_or_malformed_findings():
    report = {
        "head_sha": "head",
        "verdict": "pass",
        "summary": "Review",
        "findings": [
            {"severity": "important", "path": "file.py", "description": "Bug"}
        ],
    }
    assert not review.validate_report(report, "head")
    report["findings"][0]["severity"] = "unrecognized"
    with pytest.raises(ValueError):
        review.validate_report(report, "head")


def test_merge_refuses_without_review(monkeypatch, tmp_path):
    from infrastructure.tasks import task_integration as integration

    monkeypatch.setattr(
        integration, "ensure_worktree", lambda _: (tmp_path, "task/v3-test")
    )
    monkeypatch.setattr(integration, "ensure_clean", lambda _: None)
    monkeypatch.setattr(integration, "require_local_validation", lambda _: None)
    monkeypatch.setattr(integration, "git", lambda *a, **k: None)
    monkeypatch.setattr(integration, "git_text", lambda *a, **k: "base")

    def refused(*a, **k):
        raise ValueError("review missing")

    monkeypatch.setattr(review, "require_review", refused)
    monkeypatch.setattr(
        integration, "gh", lambda *a, **k: pytest.fail("must not merge")
    )
    with pytest.raises(ValueError, match="review missing"):
        integration.merge_pr("V3-001")
