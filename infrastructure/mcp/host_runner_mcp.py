#!/usr/bin/env python3

import json
import subprocess
import sys
from pathlib import Path

from mcp.server.mcpserver import MCPServer


REPO_ROOT = Path(__file__).resolve().parents[2]
RUNNER = (
    REPO_ROOT
    / "infrastructure"
    / "host_runner"
    / "host_runner.py"
)

mcp = MCPServer("ionmc-host-runner")


@mcp.tool()
def run_host_validation(
    task_id: str,
    argv: list[str],
    timeout_seconds: int = 600,
) -> dict:
    """
    Execute a validation command for an exact clean task-worktree SHA
    inside the controlled local host/GPU sandbox.

    The command receives GPU access but no network, Docker socket,
    user home, credentials, or arbitrary host filesystem access.

    A successful process exit is execution evidence only. It does not
    itself certify scientific validity; use the validation-manager
    separately for that.
    """

    command = [
        sys.executable,
        str(RUNNER),
        "--task-id",
        task_id,
        "--timeout-seconds",
        str(timeout_seconds),
        "--",
        *argv,
    ]

    proc = subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )

    text = proc.stdout.strip()

    try:
        result = json.loads(text)
    except json.JSONDecodeError:
        result = {
            "status": "failed",
            "exit_code": proc.returncode,
            "stdout": proc.stdout,
            "stderr": proc.stderr,
        }

    if proc.returncode != 0:
        result["runner_exit_code"] = proc.returncode

    return result


if __name__ == "__main__":
    mcp.run()
