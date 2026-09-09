#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path

from mcp.server.mcpserver import MCPServer


REPO_ROOT = Path(__file__).resolve().parents[2]
TASKCTL = REPO_ROOT / "infrastructure" / "tasks" / "task_worktree.py"

mcp = MCPServer("ionmc-task-manager")


def run_taskctl(*args: str) -> dict:
    proc = subprocess.run(
        [sys.executable, str(TASKCTL), *args],
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

    if proc.returncode != 0 and "exit_code" not in result:
        result["exit_code"] = proc.returncode

    return result


@mcp.tool()
def create_task_worktree(task_id: str, description: str = "") -> dict:
    """Create a new isolated task branch/worktree from synchronized develop."""
    args = ["create", "--task-id", task_id]
    if description:
        args += ["--description", description]
    return run_taskctl(*args)


@mcp.tool()
def inspect_task_worktree(task_id: str) -> dict:
    """Inspect task worktree branch, SHA, and cleanliness."""
    return run_taskctl("inspect", "--task-id", task_id)


@mcp.tool()
def retire_task_worktree(task_id: str) -> dict:
    """
    Remove a clean task worktree while preserving its branch.
    Dirty worktrees are refused.
    """
    return run_taskctl("retire", "--task-id", task_id)


if __name__ == "__main__":
    mcp.run()
