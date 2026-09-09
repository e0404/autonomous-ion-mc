#!/usr/bin/env python3

import sys
import argparse
import json
import re
import subprocess
from pathlib import Path


REPO = Path.home() / "aiprojects" / "ion-mc"
WORKTREE_ROOT = Path.home() / "aiprojects" / "ion-mc-worktrees"

VALIDATION_TOOL = (
    Path(__file__).resolve().parents[2]
    / "infrastructure"
    / "validation"
    / "local_validation.py"
)

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


def inspect_ci(task_id: str):
    path, branch = ensure_worktree(task_id)

    pr = find_open_pr(branch)

    if pr is None:
        raise RuntimeError(
            f"No open develop PR found for branch {branch}"
        )

    repo = repository_name()
    head_sha = git_text("rev-parse", "HEAD", cwd=path)

    remote_sha_proc = git(
        "rev-parse",
        f"origin/{branch}",
        cwd=path,
        check=False,
    )

    remote_sha = (
        remote_sha_proc.stdout.strip()
        if remote_sha_proc.returncode == 0
        else None
    )

    if remote_sha != head_sha:
        raise RuntimeError(
            "Task branch HEAD does not match pushed origin branch; "
            "CI state would not describe the exact local task SHA"
        )

    checks_proc = gh(
        "pr",
        "checks",
        str(pr["number"]),
        "--repo",
        repo,
        "--json",
        "name,state,bucket,link,workflow",
        cwd=path,
        check=False,
    )

    if checks_proc.returncode not in (0, 1, 8):
        raise RuntimeError(
            "Failed to inspect PR checks:\n"
            + (checks_proc.stderr.strip() or checks_proc.stdout.strip())
        )

    try:
        checks = json.loads(checks_proc.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"Could not parse gh pr checks output: {exc}"
        ) from exc

    runs_proc = gh(
        "run",
        "list",
        "--repo",
        repo,
        "--branch",
        branch,
        "--event",
        "pull_request",
        "--limit",
        "20",
        "--json",
        (
            "databaseId,name,workflowName,event,status,"
            "conclusion,headSha,url,createdAt,updatedAt"
        ),
        cwd=path,
    )

    runs = json.loads(runs_proc.stdout)

    matching_runs = [
        run
        for run in runs
        if run.get("headSha") == head_sha
    ]

    detailed_runs = []

    for run_info in matching_runs:
        run_id = run_info["databaseId"]

        view_proc = gh(
            "run",
            "view",
            str(run_id),
            "--repo",
            repo,
            "--json",
            "jobs",
            cwd=path,
        )

        view = json.loads(view_proc.stdout)

        detailed_runs.append(
            {
                **run_info,
                "jobs": view.get("jobs", []),
            }
        )

    result = {
        "status": "observed",
        "task_id": normalize_task_id(task_id),
        "branch": branch,
        "head_sha": head_sha,
        "pr_number": pr["number"],
        "pr_url": pr["url"],
        "checks": checks,
        "workflow_runs": detailed_runs,
    }

    print(json.dumps(result, indent=2))


def ci_failure_logs(
    task_id: str,
    max_chars: int = 20000,
):
    if max_chars < 1000 or max_chars > 100000:
        raise ValueError(
            "max_chars must be between 1000 and 100000"
        )

    path, branch = ensure_worktree(task_id)

    pr = find_open_pr(branch)

    if pr is None:
        raise RuntimeError(
            f"No open develop PR found for branch {branch}"
        )

    repo = repository_name()
    head_sha = git_text("rev-parse", "HEAD", cwd=path)

    remote_sha = git_text(
        "rev-parse",
        f"origin/{branch}",
        cwd=path,
    )

    if remote_sha != head_sha:
        raise RuntimeError(
            "Task branch HEAD does not match pushed origin branch"
        )

    runs_proc = gh(
        "run",
        "list",
        "--repo",
        repo,
        "--branch",
        branch,
        "--event",
        "pull_request",
        "--limit",
        "20",
        "--json",
        "databaseId,workflowName,status,conclusion,headSha,url",
        cwd=path,
    )

    runs = json.loads(runs_proc.stdout)

    matching_runs = [
        run
        for run in runs
        if run.get("headSha") == head_sha
    ]

    failures = []

    for run_info in matching_runs:
        run_id = run_info["databaseId"]

        view_proc = gh(
            "run",
            "view",
            str(run_id),
            "--repo",
            repo,
            "--json",
            "jobs",
            cwd=path,
        )

        jobs = json.loads(view_proc.stdout).get("jobs", [])

        for job in jobs:
            if job.get("conclusion") not in {
                "failure",
                "cancelled",
                "timed_out",
                "action_required",
            }:
                continue

            job_id = job.get("databaseId")

            if job_id is None:
                continue

            log_proc = gh(
                "run",
                "view",
                "--repo",
                repo,
                "--job",
                str(job_id),
                "--log-failed",
                cwd=path,
                check=False,
            )

            raw_log = log_proc.stdout or log_proc.stderr
            truncated = len(raw_log) > max_chars

            if truncated:
                # Failure is normally near the end of the job log.
                log = raw_log[-max_chars:]
            else:
                log = raw_log

            failures.append(
                {
                    "workflow": run_info.get("workflowName"),
                    "run_id": run_id,
                    "run_url": run_info.get("url"),
                    "job_id": job_id,
                    "job_name": job.get("name"),
                    "job_conclusion": job.get("conclusion"),
                    "steps": job.get("steps", []),
                    "log": log,
                    "log_truncated": truncated,
                }
            )

    result = {
        "status": "observed",
        "task_id": normalize_task_id(task_id),
        "branch": branch,
        "head_sha": head_sha,
        "pr_number": pr["number"],
        "failures": failures,
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

def require_local_validation(task_id: str):
    proc = subprocess.run(
        [
            sys.executable,
            str(VALIDATION_TOOL),
            "check",
            "--task-id",
            normalize_task_id(task_id),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
    )

    if proc.returncode != 0:
        detail = proc.stdout.strip() or proc.stderr.strip()
        raise RuntimeError(
            "Local validation gate has not passed for the exact "
            "task-branch SHA:\n" + detail
        )

def merge_pr(task_id: str):
    path, branch = ensure_worktree(task_id)
    ensure_clean(path)

    require_local_validation(task_id)

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

    p_ci = sub.add_parser("inspect-ci")
    p_ci.add_argument("--task-id", required=True)

    p_ci_logs = sub.add_parser("ci-failure-logs")
    p_ci_logs.add_argument("--task-id", required=True)
    p_ci_logs.add_argument(
        "--max-chars",
        type=int,
        default=20000,
    )

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

        elif args.command == "inspect-ci":
            inspect_ci(args.task_id)

        elif args.command == "ci-failure-logs":
            ci_failure_logs(
                args.task_id,
                max_chars=args.max_chars,
            )

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
