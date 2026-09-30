#!/usr/bin/env python3
"""Release service keeps direct GitHub CLI access outside the autonomous shell."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mcp.server.mcpserver import MCPServer

from infrastructure.experiment_v3.common import ROOT, STATE, safe_path
from infrastructure.experiment_v3.promote import promote
from infrastructure.experiment_v3.release import evaluate

mcp = MCPServer("ionmc-v3-release")


@mcp.tool()
def run_release_validation(argv: list[str], timeout_seconds: int = 600) -> dict:
    """Execute exact integration-SHA qualification in the existing host sandbox.

    Store reports under ignored validation/generated. A successful process is
    execution evidence only; qualification independently checks its artifacts.
    """
    from infrastructure.host_runner.host_runner import run_validation

    return run_validation("V3-RELEASE", argv, timeout_seconds, release_root=ROOT)


@mcp.tool()
def stage_release_evidence(report_path: str, artifact_path: str) -> dict:
    """Archive release evidence from ignored integration-checkout output paths."""
    import json
    import shutil

    from infrastructure.experiment_v3.common import exact_state, file_hash, write_json

    if not (
        report_path.startswith("validation/generated/")
        and (
            artifact_path == "validation/generated"
            or artifact_path.startswith("validation/generated/")
        )
    ):
        raise ValueError("Stage only ignored validation/generated outputs")
    sha = exact_state(ROOT)
    report_file = safe_path(ROOT, report_path)
    artifact_root = safe_path(ROOT, artifact_path)
    report = json.loads(report_file.read_text())
    if report.get("code_sha") != sha or report.get("dirty") is not False:
        raise ValueError("Evidence must describe the exact clean integration SHA")
    report_id = file_hash(report_file)
    destination = STATE / "release-evidence" / report_id
    refs = [a for suite in report.get("suites", []) for a in suite.get("artifacts", [])]
    refs += [
        a
        for workload in report.get("performance", {}).values()
        for a in workload.get("measurement_artifacts", [])
    ]
    for artifact in refs:
        source = safe_path(artifact_root, artifact["path"])
        if file_hash(source) != artifact["sha256"]:
            raise ValueError("Artifact checksum mismatch")
        target = safe_path(destination / "artifacts", artifact["path"])
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    write_json(destination / "report.json", report)
    return {
        "report_path": str((destination / "report.json").relative_to(STATE)),
        "artifact_path": str((destination / "artifacts").relative_to(STATE)),
        "source_report_sha256": report_id,
    }


@mcp.tool()
def qualify_release(
    freeze_sha: str,
    report_path: str,
    artifact_path: str,
    plan_path: str = "validation/release-plan.json",
) -> dict:
    """Evaluate full exact-SHA evidence; paths are relative to v2 state."""
    return evaluate(
        ROOT,
        plan_path,
        freeze_sha,
        safe_path(STATE, report_path),
        safe_path(STATE, artifact_path),
    )


@mcp.tool()
def promote_release(
    freeze_sha: str,
    report_path: str,
    artifact_path: str,
    tag: str,
    plan_path: str = "validation/release-plan.json",
) -> dict:
    """Requalify, create/merge release PR after exact-head CI, tag and publish release.

    Restricted to the experiment's configured remote and v2 permanent branches.
    Failures are unfinished release work; rerun after resolving the actual failure.
    """
    return promote(
        ROOT,
        plan_path,
        freeze_sha,
        safe_path(STATE, report_path),
        safe_path(STATE, artifact_path),
        tag,
    )


if __name__ == "__main__":
    mcp.run()
