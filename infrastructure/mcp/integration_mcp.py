#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path

from mcp.server.mcpserver import MCPServer


REPO_ROOT = Path(__file__).resolve().parents[2]

INTEGRATION_TOOL = (
    REPO_ROOT
    / "infrastructure"
    / "tasks"
    / "task_integration.py"
)

mcp = MCPServer("ionmc-task-integration")


def run_integration(*args: str) -> dict:
    proc = subprocess.run(
        [sys.executable, str(INTEGRATION_TOOL), *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    stream = proc.stdout if proc.stdout.strip() else proc.stderr

    try:
        result = json.loads(stream)
    except json.JSONDecodeError:
        result = {
            "status": "failed",
            "exit_code": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
        }

    if proc.returncode != 0:
        result["exit_code"] = proc.returncode

    return result


@mcp.tool()
def push_task_branch(task_id: str) -> dict:
    """Push a clean committed task branch to origin."""
    return run_integration(
        "push",
        "--task-id",
        task_id,
    )


@mcp.tool()
def create_task_pull_request(
    task_id: str,
    title: str,
    body: str,
) -> dict:
    """Create the task branch pull request into develop."""
    return run_integration(
        "create-pr",
        "--task-id",
        task_id,
        "--title",
        title,
        "--body",
        body,
    )


@mcp.tool()
def inspect_task_pull_request(task_id: str) -> dict:
    """Inspect the task's open pull request and local worktree state."""
    return run_integration(
        "inspect-pr",
        "--task-id",
        task_id,
    )


@mcp.tool()
def merge_task_pull_request(task_id: str) -> dict:
    """
    Squash-merge the task PR into develop and fast-forward the
    primary local develop checkout.

    The task branch is preserved.
    """
    return run_integration(
        "merge-pr",
        "--task-id",
        task_id,
    )


if __name__ == "__main__":
    mcp.run()
