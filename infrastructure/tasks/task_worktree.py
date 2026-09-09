#!/usr/bin/env python3

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path


REPO = Path.home() / "aiprojects" / "ion-mc"
WORKTREE_ROOT = Path.home() / "aiprojects" / "ion-mc-worktrees"


def git(*args: str, cwd: Path = REPO, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", *args],
        cwd=cwd,
        text=True,
        capture_output=True,
        check=check,
    )


def git_text(*args: str, cwd: Path = REPO) -> str:
    return git(*args, cwd=cwd).stdout.strip()


def slugify(value: str) -> str:
    value = value.lower().strip()
    value = re.sub(r"[^a-z0-9]+", "-", value)
    return value.strip("-")[:60]


def normalize_task_id(task_id: str) -> str:
    task_id = task_id.strip().upper()
    if not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{1,63}", task_id):
        raise ValueError(
            "task_id must contain only A-Z, 0-9, '.', '_' or '-', "
            "must start alphanumerically, and be 2-64 characters long"
        )
    return task_id


def branch_for(task_id: str, description: str | None = None) -> str:
    base = task_id.lower()
    if description:
        slug = slugify(description)
        if slug:
            return f"task/{base}-{slug}"
    return f"task/{base}"


def path_for(task_id: str) -> Path:
    return WORKTREE_ROOT / task_id.lower()


def ensure_repo_ready():
    branch = git_text("branch", "--show-current")
    if branch != "develop":
        raise RuntimeError(
            f"Main checkout must be on develop, currently on {branch!r}"
        )

    status = git_text("status", "--porcelain")
    if status:
        raise RuntimeError("Main develop checkout is not clean")

    git("fetch", "origin", "develop")

    local = git_text("rev-parse", "develop")
    remote = git_text("rev-parse", "origin/develop")

    if local != remote:
        raise RuntimeError(
            "Local develop is not identical to origin/develop. "
            "Synchronize it before creating a task."
        )


def worktree_entries():
    text = git_text("worktree", "list", "--porcelain")
    entries = []
    current = {}

    for line in text.splitlines():
        if not line:
            if current:
                entries.append(current)
                current = {}
            continue

        key, _, value = line.partition(" ")
        current[key] = value

    if current:
        entries.append(current)

    return entries


def find_worktree(path: Path):
    resolved = str(path.resolve())
    for entry in worktree_entries():
        if entry.get("worktree") == resolved:
            return entry
    return None


def create(task_id: str, description: str | None):
    task_id = normalize_task_id(task_id)
    ensure_repo_ready()

    branch = branch_for(task_id, description)
    path = path_for(task_id)

    if path.exists():
        raise RuntimeError(f"Task worktree path already exists: {path}")

    if git("show-ref", "--verify", "--quiet", f"refs/heads/{branch}", check=False).returncode == 0:
        raise RuntimeError(f"Local branch already exists: {branch}")

    WORKTREE_ROOT.mkdir(parents=True, exist_ok=True)

    git("worktree", "add", "-b", branch, str(path), "develop")

    result = {
        "status": "created",
        "task_id": task_id,
        "branch": branch,
        "worktree": str(path),
        "base_branch": "develop",
        "base_sha": git_text("rev-parse", "develop"),
        "head_sha": git_text("rev-parse", "HEAD", cwd=path),
    }

    print(json.dumps(result, indent=2))


def inspect(task_id: str):
    task_id = normalize_task_id(task_id)
    path = path_for(task_id)

    if not path.exists():
        raise RuntimeError(f"No worktree directory exists for {task_id}")

    entry = find_worktree(path)
    if entry is None:
        raise RuntimeError(f"Directory is not a registered Git worktree: {path}")

    branch = git_text("branch", "--show-current", cwd=path)
    status = git_text("status", "--porcelain", cwd=path)

    result = {
        "status": "present",
        "task_id": task_id,
        "branch": branch,
        "worktree": str(path),
        "head_sha": git_text("rev-parse", "HEAD", cwd=path),
        "dirty": bool(status),
        "porcelain_status": status,
    }

    print(json.dumps(result, indent=2))


def retire(task_id: str):
    task_id = normalize_task_id(task_id)
    path = path_for(task_id)

    if not path.exists():
        raise RuntimeError(f"No worktree directory exists for {task_id}")

    entry = find_worktree(path)
    if entry is None:
        raise RuntimeError(f"Directory is not a registered Git worktree: {path}")

    status = git_text("status", "--porcelain", cwd=path)
    if status:
        raise RuntimeError(
            "Refusing to retire dirty worktree:\n" + status
        )

    branch = git_text("branch", "--show-current", cwd=path)
    head_sha = git_text("rev-parse", "HEAD", cwd=path)

    git("worktree", "remove", str(path))

    result = {
        "status": "retired",
        "task_id": task_id,
        "branch": branch,
        "head_sha": head_sha,
        "worktree_removed": str(path),
        "branch_preserved": True,
    }

    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    p_create = sub.add_parser("create")
    p_create.add_argument("--task-id", required=True)
    p_create.add_argument("--description")

    p_inspect = sub.add_parser("inspect")
    p_inspect.add_argument("--task-id", required=True)

    p_retire = sub.add_parser("retire")
    p_retire.add_argument("--task-id", required=True)

    args = parser.parse_args()

    try:
        if args.command == "create":
            create(args.task_id, args.description)
        elif args.command == "inspect":
            inspect(args.task_id)
        elif args.command == "retire":
            retire(args.task_id)
    except Exception as exc:
        print(
            json.dumps(
                {
                    "status": "failed",
                    "error": str(exc),
                },
                indent=2,
            ),
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
