#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path

from mcp.server.mcpserver import MCPServer


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

    wt = Path(worktree).resolve()

    allowed_root = Path.home() / "aiprojects"

    try:
        wt.relative_to(allowed_root)
    except ValueError as exc:
        raise ValueError(
            f"Worktree must be below {allowed_root}"
        ) from exc

    if not (wt / ".git").exists():
        raise ValueError(f"Not a Git worktree: {wt}")

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


if __name__ == "__main__":
    mcp.run()
