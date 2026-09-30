#!/usr/bin/env python3
"""Trusted v2 reference/research service. Run outside the agent shell sandbox."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mcp.server.mcpserver import MCPServer

from infrastructure.experiment_v3.artifacts import (
    list_artifacts,
    materialize,
    read_text,
)
from infrastructure.experiment_v3.common import ROOT, STATE, event, exact_state
from infrastructure.experiment_v3.data import acquire, assign_role
from infrastructure.experiment_v3.reference import run_reference
from infrastructure.host_runner.host_runner import inspect_worktree

mcp = MCPServer("ionmc-v3-research")


@mcp.tool()
def list_reference_engines() -> dict:
    """Report provisioned versions. A listed engine is not scientifically validated."""
    path = STATE / "engines.json"
    if not path.exists():
        return {"status": "not-provisioned", "engines": []}
    data = json.loads(path.read_text())
    return {
        "engines": {
            name: {"version": c["version"], "source": c["source"]}
            for name, c in data["engines"].items()
        }
    }


@mcp.tool()
def run_reference_calculation(
    task_id: str,
    case_path: str,
    engine: str,
    timeout_seconds: int = 600,
    gpu: bool = False,
) -> dict:
    """Run a committed native reference case, preserving raw outputs and exact SHA.

    case_path is relative to the clean task worktree, containing case.json and
    engine-native inputs. No arbitrary host paths or registry overrides accepted.
    Execution success is not physical validation.
    """
    worktree = (
        {"path": ROOT, "sha": exact_state(ROOT)}
        if task_id == "V3-RELEASE"
        else inspect_worktree(task_id)
    )
    return run_reference(
        worktree["path"],
        case_path,
        engine,
        timeout=timeout_seconds,
        gpu=gpu,
        task_id=task_id,
    )


@mcp.tool()
def read_reference_artifact(
    run_id: str, path: str, offset: int = 0, limit: int = 2000
) -> dict:
    """Preview UTF-8 text (max 4000 bytes).

    Use listing and direct materialization for full/binary files.
    Never reconstruct base64 in conversation.
    """
    return read_text(run_id, path, offset=offset, limit=limit)


@mcp.tool()
def list_reference_artifacts(run_id: str, offset: int = 0, limit: int = 50) -> dict:
    """List exact paths, sizes and hashes; paginate instead of guessing names."""
    return list_artifacts(run_id, offset=offset, limit=limit)


@mcp.tool()
def materialize_reference_artifacts(
    run_id: str, task_id: str, paths: list[str] | None = None
) -> dict:
    """Copy verified outputs directly into the ignored task cache.

    Return only paths/counts. Analyze locally and return compact summaries.
    """
    return materialize(run_id, inspect_worktree(task_id)["path"], paths)


@mcp.tool()
def acquire_reference_data(
    url: str,
    license_basis: str,
    citation: str,
    role: str,
    rationale: str,
    expected_sha256: str | None = None,
    post_body: str | None = None,
) -> dict:
    """Acquire public HTTPS scientific data. No credentials or private URLs.

    post_body supports public form endpoints such as NIST; never send credentials,
    personal or institutional data. Roles disclose construction/calibration reuse.
    """
    return acquire(
        url,
        license_basis,
        citation,
        role,
        rationale,
        sha=exact_state(ROOT),
        expected_sha256=expected_sha256,
        post_body=post_body,
    )


@mcp.tool()
def assign_reference_data_role(content_sha256: str, role: str, rationale: str) -> dict:
    """Append a dataset-use role without erasing earlier uses."""
    return assign_role(content_sha256, role, rationale, sha=exact_state(ROOT))


@mcp.tool()
def materialize_reference_data(content_sha256: str, task_id: str) -> dict:
    """Copy verified reference data to the ignored task cache."""
    import shutil

    from infrastructure.experiment_v3.common import file_hash

    if len(content_sha256) != 64 or any(
        c not in "0123456789abcdef" for c in content_sha256
    ):
        raise ValueError("Invalid data hash")
    worktree = inspect_worktree(task_id)
    source = STATE / "data/objects" / content_sha256
    if file_hash(source) != content_sha256:
        raise ValueError("Cached content changed")
    destination = worktree["path"] / ".ionmc-cache/reference" / content_sha256
    from infrastructure.experiment_v3.common import safe_path

    safe_path(worktree["path"], ".ionmc-cache/reference/" + content_sha256)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    return {
        "path": str(destination),
        "sha256": content_sha256,
        "bytes": destination.stat().st_size,
    }


@mcp.tool()
def record_scientific_event(
    event_name: str,
    task_id: str,
    rationale: str,
    artifact_ids: list[str],
    previous_decision: str = "",
    new_decision: str = "",
) -> dict:
    """Record a scientific failure, contradiction or consequential strategy change."""
    allowed = {
        "validation_strategy_changed",
        "internal_evidence_rejected",
        "scientific_failure",
        "reference_contradiction",
        "model_changed",
        "tolerance_changed",
        "performance_failure",
        "workflow_failure",
    }
    if event_name not in allowed or not rationale.strip() or not artifact_ids:
        raise ValueError(
            "Scientific event, rationale and artifact identifiers required"
        )
    return event(
        event_name,
        sha=inspect_worktree(task_id)["sha"],
        task_id=task_id,
        details={
            "rationale": rationale,
            "artifact_ids": artifact_ids,
            "previous_decision": previous_decision,
            "new_decision": new_decision,
        },
    )


if __name__ == "__main__":
    mcp.run()
