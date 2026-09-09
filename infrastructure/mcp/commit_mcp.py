#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path

from mcp.server.mcpserver import MCPServer


REPO_ROOT = Path(__file__).resolve().parents[2]

COMMIT_TOOL = (
    REPO_ROOT
    / "infrastructure"
    / "tasks"
    / "task_commit.py"
)

mcp = MCPServer("ionmc-task-commit")


@mcp.tool()
def commit_task_changes(
    task_id: str,
    message: str,
) -> dict:
    """
    Stage and commit the current task-worktree changes using the fixed
    autonomous experiment Git identity.

    Use only after reviewing the worktree changes and running appropriate
    tests before the commit.
    """

    proc = subprocess.run(
        [
            sys.executable,
            str(COMMIT_TOOL),
            "--task-id",
            task_id,
            "--message",
            message,
        ],
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
