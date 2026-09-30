#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mcp.server.mcpserver import MCPServer
from infrastructure.experiment_v3.review import start, inspect
from infrastructure.experiment_v3.worker_boundary import validate_task


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKER = REPO_ROOT / "infrastructure" / "codex" / "codex_worker.py"

mcp = MCPServer("ionmc-codex-worker")


@mcp.tool()
def codex_worker(
    task_id: str,
    worktree: str,
    prompt: str,
    model: str | None = None,
) -> dict:
    """
    Run a controlled Codex worker in an existing isolated Git worktree.

    Args:
        task_id: Experiment/task identifier.
        worktree: Absolute path to an existing Git worktree.
        prompt: Full instruction for the Codex worker.
        model: Optional explicit Codex model.
    """

    task_id, wt = validate_task(task_id, worktree)

    cmd = [
        sys.executable,
        str(WORKER),
        "--worktree",
        str(wt),
        "--task-id",
        task_id,
        "--prompt",
        prompt,
    ]

    if model:
        cmd += ["--model", model]

    proc = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    if not proc.stdout.strip():
        return {
            "status": "failed",
            "exit_code": proc.returncode,
            "stderr": proc.stderr,
        }

    try:
        result = json.loads(proc.stdout)
    except json.JSONDecodeError:
        return {
            "status": "failed",
            "exit_code": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
        }

    if proc.stderr:
        result["stderr"] = proc.stderr

    return result


@mcp.tool()
def start_codex_review(task_id: str) -> dict:
    """Start mandatory independent review of a clean exact-SHA task. Poll inspect_codex_review; fix findings before merge. Generic codex_worker calls do not satisfy this gate."""
    return start(task_id)


@mcp.tool()
def inspect_codex_review(review_id: str) -> dict:
    """Read protected review status and structured findings, without raw token logs."""
    return inspect(review_id)


if __name__ == "__main__":
    mcp.run()
