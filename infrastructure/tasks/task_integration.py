#!/usr/bin/env python3

import argparse
import json
import re
import subprocess
from pathlib import Path


REPO = Path.home() / "aiprojects" / "ion-mc"
WORKTREE_ROOT = Path.home() / "aiprojects" / "ion-mc-worktrees"


def run(
    command: list[str],
    *,
    cwd: Path = REPO,
    check: bool = True,
) -> subprocess.CompletedProcess:
    return subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=check,
    )


def git(*args: str, cwd: Path = REPO, check: bool = True):
    return run(["git", *args], cwd=cwd, check=check)


def git_text(*args: str, cwd: Path = REPO) -> str:
    return git(*args, cwd=cwd).stdout.strip()


def gh(*args: str, cwd: Path = REPO, check: bool = True):
    return run(["gh", *args], cwd=cwd, check=check)


def normalize_task_id(task_id: str) -> str:
    value = task_id.strip().upper()

    if not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{1,63}", value):
        raise ValueError(f"Invalid task ID: {task_id}")

    return value


def worktree_for(task_id: str) -> Path:
    return WORKTREE_ROOT / normalize_task_id(task_id).lower()


def ensure_worktree(task_id: str) -> tuple[Path, str]:
    path = worktree_for(task_id)

    if not path.exists():
        raise RuntimeError(f"Task worktree does not exist: {path}")

    branch = git_text("branch", "--show-current", cwd=path)

    if not branch.startswith(("task/", "fix/", "experiment/")):
        raise RuntimeError(
            f"Refusing integration from unexpected branch: {branch}"
        )

    return path, branch


def ensure_clean(path: Path):
    status = git_text("status", "--porcelain", cwd=path)

    if status:
        raise RuntimeError(
            "Task worktree must be clean before integration:\n" + status
        )


def repository_name() -> str:
    remote = git_text("remote", "get-url", "origin")

    match = re.fullmatch(
        r"https://github\.com/([^/]+)/([^/]+?)(?:\.git)?",
        remote,
    )

    if not match:
        raise RuntimeError(
            f"Unsupported GitHub origin URL: {remote}"
        )

    return f"{match.group(1)}/{match.group(2)}"


def commits_ahead(path: Path) -> int:
    output = git_text(
        "rev-list",
        "--count",
        "origin/develop..HEAD",
        cwd=path,
    )
    return int(output)


def find_open_pr(branch: str) -> dict | None:
    repo = repository_name()

    proc = gh(
        "pr",
        "list",
        "--repo",
        repo,
        "--head",
        branch,
        "--base",
        "develop",
        "--state",
        "open",
        "--json",
        "number,title,url,state,mergeStateStatus,headRefName,baseRefName",
    )

    prs = json.loads(proc.stdout)

    if not prs:
        return None

    if len(prs) != 1:
        raise RuntimeError(
            f"Expected exactly one open PR for {branch}, found {len(prs)}"
        )

    return prs[0]


def push_task(task_id: str):
    path, branch = ensure_worktree(task_id)
    ensure_clean(path)

    git("fetch", "origin", "develop", cwd=path)

    ahead = commits_ahead(path)

    if ahead <= 0:
        raise RuntimeError(
            "Task branch has no commits ahead of origin/develop"
        )

    proc = git(
        "push",
        "-u",
        "origin",
        branch,
        cwd=path,
    )

    result = {
        "status": "pushed",
        "task_id": normalize_task_id(task_id),
        "branch": branch,
        "worktree": str(path),
        "head_sha": git_text("rev-parse", "HEAD", cwd=path),
        "commits_ahead_of_develop": ahead,
        "git_stdout": proc.stdout.strip(),
        "git_stderr": proc.stderr.strip(),
    }

    print(json.dumps(result, indent=2))


def create_pr(task_id: str, title: str, body: str):
    path, branch = ensure_worktree(task_id)
    ensure_clean(path)

    existing = find_open_pr(branch)

    if existing is not None:
        raise RuntimeError(
            f"Open PR already exists: #{existing['number']} {existing['url']}"
        )

    repo = repository_name()

    proc = gh(
        "pr",
        "create",
        "--repo",
        repo,
        "--base",
        "develop",
        "--head",
        branch,
        "--title",
        title,
        "--body",
        body,
        cwd=path,
    )

    pr = find_open_pr(branch)

    if pr is None:
        raise RuntimeError(
            "gh pr create succeeded but the new PR could not be found"
        )

    result = {
        "status": "created",
        "task_id": normalize_task_id(task_id),
        "branch": branch,
        "pr": pr,
        "gh_output": proc.stdout.strip(),
    }

    print(json.dumps(result, indent=2))


def inspect_pr(task_id: str):
    path, branch = ensure_worktree(task_id)

    pr = find_open_pr(branch)

    result = {
        "status": "open" if pr is not None else "not-found",
        "task_id": normalize_task_id(task_id),
        "branch": branch,
        "worktree": str(path),
        "worktree_dirty": bool(
            git_text("status", "--porcelain", cwd=path)
        ),
        "head_sha": git_text("rev-parse", "HEAD", cwd=path),
        "pr": pr,
    }

    print(json.dumps(result, indent=2))


def sync_develop():
    branch = git_text("branch", "--show-current", cwd=REPO)

    if branch != "develop":
        raise RuntimeError(
            f"Primary checkout must be on develop, currently {branch}"
        )

    status = git_text("status", "--porcelain", cwd=REPO)

    if status:
        raise RuntimeError(
            "Primary develop checkout is dirty; refusing synchronization"
        )

    git("fetch", "origin", "develop", cwd=REPO)
    git("merge", "--ff-only", "origin/develop", cwd=REPO)


def merge_pr(task_id: str):
    path, branch = ensure_worktree(task_id)
    ensure_clean(path)

    pr = find_open_pr(branch)

    if pr is None:
        raise RuntimeError(
            f"No open develop PR found for branch {branch}"
        )

    repo = repository_name()

    proc = gh(
        "pr",
        "merge",
        str(pr["number"]),
        "--repo",
        repo,
        "--squash",
        "--delete-branch=false",
        cwd=path,
    )

    sync_develop()

    result = {
        "status": "merged",
        "task_id": normalize_task_id(task_id),
        "branch": branch,
        "pr_number": pr["number"],
        "pr_url": pr["url"],
        "develop_sha": git_text("rev-parse", "develop", cwd=REPO),
        "branch_preserved": True,
        "gh_output": proc.stdout.strip(),
    }

    print(json.dumps(result, indent=2))


def main() -> int:
    parser = argparse.ArgumentParser()

    sub = parser.add_subparsers(
        dest="command",
        required=True,
    )

    p_push = sub.add_parser("push")
    p_push.add_argument("--task-id", required=True)

    p_create = sub.add_parser("create-pr")
    p_create.add_argument("--task-id", required=True)
    p_create.add_argument("--title", required=True)
    p_create.add_argument("--body", required=True)

    p_inspect = sub.add_parser("inspect-pr")
    p_inspect.add_argument("--task-id", required=True)

    p_merge = sub.add_parser("merge-pr")
    p_merge.add_argument("--task-id", required=True)

    args = parser.parse_args()

    try:
        if args.command == "push":
            push_task(args.task_id)

        elif args.command == "create-pr":
            create_pr(
                args.task_id,
                args.title,
                args.body,
            )

        elif args.command == "inspect-pr":
            inspect_pr(args.task_id)

        elif args.command == "merge-pr":
            merge_pr(args.task_id)

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
