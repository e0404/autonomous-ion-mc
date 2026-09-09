#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path

from mcp.server.mcpserver import MCPServer


REPO_ROOT = Path(__file__).resolve().parents[2]

VALIDATION_TOOL = (
    REPO_ROOT
    / "infrastructure"
    / "validation"
    / "local_validation.py"
)

mcp = MCPServer("ionmc-validation-manager")


def run_validation(*args: str) -> dict:
    proc = subprocess.run(
        [sys.executable, str(VALIDATION_TOOL), *args],
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
def record_local_validation(
    task_id: str,
    suite: str,
    status: str,
    summary: str,
    details_json: str = "{}",
) -> dict:
    """
    Record the result of a local validation suite for the exact committed
    task-branch SHA.

    The caller must truthfully report validation actually performed.
    Do not fabricate passing validation.
    """
    return run_validation(
        "record",
        "--task-id",
        task_id,
        "--suite",
        suite,
        "--status",
        status,
        "--summary",
        summary,
        "--details-json",
        details_json,
    )


@mcp.tool()
def inspect_local_validation(task_id: str) -> dict:
    """Inspect local validation records for the current task-branch SHA."""
    return run_validation(
        "inspect",
        "--task-id",
        task_id,
    )


if __name__ == "__main__":
    mcp.run()
