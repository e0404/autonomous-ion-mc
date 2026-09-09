#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path

from mcp.server.mcpserver import MCPServer


REPO_ROOT = Path(__file__).resolve().parents[2]

REQUEST_TOOL = (
    REPO_ROOT
    / "infrastructure"
    / "interventions"
    / "request_intervention.py"
)

mcp = MCPServer("ionmc-human-intervention")


@mcp.tool()
def request_human_intervention(
    task_id: str,
    category: str,
    summary: str,
    reason: str,
    requested_input: str,
    depends_on: str = "",
) -> dict:
    """
    Request intervention from the designated experiment operator.

    This tool may only be used for intervention categories permitted by
    EXPERIMENT.md. Scientific or technical uncertainty alone is not a
    valid reason for human intervention.
    """

    allowed_categories = {
        "inaccessible-information",
        "external-reference-calculation",
        "external-data",
        "authorization",
        "experiment-boundary",
    }

    if category not in allowed_categories:
        return {
            "status": "failed",
            "error": f"Invalid intervention category: {category}",
        }

    cmd = [
        sys.executable,
        str(REQUEST_TOOL),
        "--task-id",
        task_id,
        "--category",
        category,
        "--summary",
        summary,
        "--reason",
        reason,
        "--requested-input",
        requested_input,
    ]

    if depends_on:
        cmd += ["--depends-on", depends_on]

    proc = subprocess.run(
        cmd,
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


if __name__ == "__main__":
    mcp.run()
