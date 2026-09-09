#!/usr/bin/env python3

import argparse
import json
import re
import subprocess
from pathlib import Path


WORKTREE_ROOT = Path.home() / "aiprojects" / "ion-mc-worktrees"

AGENT_NAME = "Autonomous IonMC Agent"
AGENT_EMAIL = "autonomous-ionmc-agent@users.noreply.github.com"


def git(
    *args: str,
    cwd: Path,
    check: bool = True,
) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", *args],
        cwd=cwd,
        text=True,
        capture_output=True,
        check=check,
    )


def git_text(*args: str, cwd: Path) -> str:
    return git(*args, cwd=cwd).stdout.strip()


def normalize_task_id(task_id: str) -> str:
    value = task_id.strip().upper()

    if not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{1,63}", value):
        raise ValueError(f"Invalid task ID: {task_id}")

    return value


def worktree_for(task_id: str) -> Path:
    return WORKTREE_ROOT / normalize_task_id(task_id).lower()


def ensure_task_worktree(task_id: str) -> tuple[Path, str]:
    path = worktree_for(task_id)

    if not path.exists():
        raise RuntimeError(f"Task worktree does not exist: {path}")

    branch = git_text("branch", "--show-current", cwd=path)

    if not branch.startswith(("task/", "fix/", "experiment/")):
        raise RuntimeError(
            f"Refusing to commit on unexpected branch: {branch}"
        )

    return path, branch


def commit(task_id: str, message: str) -> None:
    path, branch = ensure_task_worktree(task_id)

    status_before = git_text("status", "--porcelain", cwd=path)

    if not status_before:
        raise RuntimeError("Task worktree has no changes to commit")

    # Stage the complete task worktree state. Repository .gitignore rules
    # are expected to exclude per-session/runtime artifacts.
    git("add", "-A", cwd=path)

    staged = git_text("diff", "--cached", "--name-status", cwd=path)

    if not staged:
        raise RuntimeError("No staged changes remain after git add -A")

    proc = git(
        "-c",
        f"user.name={AGENT_NAME}",
        "-c",
        f"user.email={AGENT_EMAIL}",
        "commit",
        "-m",
        message,
        cwd=path,
    )

    sha = git_text("rev-parse", "HEAD", cwd=path)

    result = {
        "status": "committed",
        "task_id": normalize_task_id(task_id),
        "branch": branch,
        "worktree": str(path),
        "commit_sha": sha,
        "author_name": AGENT_NAME,
        "author_email": AGENT_EMAIL,
        "files": staged.splitlines(),
        "git_output": proc.stdout.strip(),
    }

    print(json.dumps(result, indent=2))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", required=True)
    parser.add_argument("--message", required=True)
    args = parser.parse_args()

    try:
        commit(args.task_id, args.message)
    except Exception as exc:
        print(
            json.dumps(
                {
                    "status": "failed",
                    "error": str(exc),
                },
                indent=2,
            )
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
