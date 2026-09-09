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

    # ------------------------------------------------------------------
    # Check runs for the exact pushed SHA
    # ------------------------------------------------------------------

    checks_proc = gh(
        "api",
        "-H",
        "Accept: application/vnd.github+json",
        f"repos/{repo}/commits/{head_sha}/check-runs",
        cwd=path,
    )

    try:
        checks_payload = json.loads(checks_proc.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            "Could not parse GitHub check-runs response:\n"
            f"stdout: {checks_proc.stdout!r}\n"
            f"stderr: {checks_proc.stderr!r}"
        ) from exc

    checks = [
        {
            "id": check.get("id"),
            "name": check.get("name"),
            "status": check.get("status"),
            "conclusion": check.get("conclusion"),
            "started_at": check.get("started_at"),
            "completed_at": check.get("completed_at"),
            "details_url": check.get("details_url"),
            "app": (
                check.get("app", {}).get("name")
                if check.get("app")
                else None
            ),
        }
        for check in checks_payload.get("check_runs", [])
    ]

    # ------------------------------------------------------------------
    # Workflow runs for the exact pushed SHA
    # ------------------------------------------------------------------

    runs_proc = gh(
        "api",
        "--method",
        "GET",
        "-H",
        "Accept: application/vnd.github+json",
        f"repos/{repo}/actions/runs",
        "-f",
        f"head_sha={head_sha}",
        "-f",
        "event=pull_request",
        "-f",
        "per_page=100",
        cwd=path,
    )

    try:
        runs_payload = json.loads(runs_proc.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            "Could not parse GitHub workflow-runs response:\n"
            f"stdout: {runs_proc.stdout!r}\n"
            f"stderr: {runs_proc.stderr!r}"
        ) from exc

    matching_runs = runs_payload.get("workflow_runs", [])

    detailed_runs = []

    for run_info in matching_runs:
        run_id = run_info["id"]

        jobs_proc = gh(
            "api",
            "--method",
            "GET",
            "-H",
            "Accept: application/vnd.github+json",
            f"repos/{repo}/actions/runs/{run_id}/jobs",
            "-f",
            "per_page=100",
            cwd=path,
        )

        try:
            jobs_payload = json.loads(jobs_proc.stdout)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"Could not parse jobs for workflow run {run_id}:\n"
                f"stdout: {jobs_proc.stdout!r}\n"
                f"stderr: {jobs_proc.stderr!r}"
            ) from exc

        normalized_jobs = []

        for job in jobs_payload.get("jobs", []):
            normalized_jobs.append(
                {
                    "id": job.get("id"),
                    "name": job.get("name"),
                    "status": job.get("status"),
                    "conclusion": job.get("conclusion"),
                    "started_at": job.get("started_at"),
                    "completed_at": job.get("completed_at"),
                    "html_url": job.get("html_url"),
                    "steps": [
                        {
                            "name": step.get("name"),
                            "status": step.get("status"),
                            "conclusion": step.get("conclusion"),
                            "number": step.get("number"),
                            "started_at": step.get("started_at"),
                            "completed_at": step.get("completed_at"),
                        }
                        for step in job.get("steps", [])
                    ],
                }
            )

        detailed_runs.append(
            {
                "id": run_info.get("id"),
                "name": run_info.get("name"),
                "workflow_name": run_info.get("name"),
                "event": run_info.get("event"),
                "status": run_info.get("status"),
                "conclusion": run_info.get("conclusion"),
                "head_sha": run_info.get("head_sha"),
                "html_url": run_info.get("html_url"),
                "created_at": run_info.get("created_at"),
                "updated_at": run_info.get("updated_at"),
                "jobs": normalized_jobs,
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
            "Task branch HEAD does not match pushed origin branch"
        )

    # ------------------------------------------------------------------
    # Workflow runs for exact pushed SHA
    # ------------------------------------------------------------------

    runs_proc = gh(
        "api",
        "--method",
        "GET",
        "-H",
        "Accept: application/vnd.github+json",
        f"repos/{repo}/actions/runs",
        "-f",
        f"head_sha={head_sha}",
        "-f",
        "event=pull_request",
        "-f",
        "per_page=100",
        cwd=path,
    )

    try:
        runs_payload = json.loads(runs_proc.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            "Could not parse GitHub workflow-runs response:\n"
            f"stdout: {runs_proc.stdout!r}\n"
            f"stderr: {runs_proc.stderr!r}"
        ) from exc

    matching_runs = runs_payload.get("workflow_runs", [])
    failures = []

    for run_info in matching_runs:
        run_id = run_info["id"]

        jobs_proc = gh(
            "api",
            "--method",
            "GET",
            "-H",
            "Accept: application/vnd.github+json",
            f"repos/{repo}/actions/runs/{run_id}/jobs",
            "-f",
            "per_page=100",
            cwd=path,
        )

        try:
            jobs_payload = json.loads(jobs_proc.stdout)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"Could not parse jobs for workflow run {run_id}:\n"
                f"stdout: {jobs_proc.stdout!r}\n"
                f"stderr: {jobs_proc.stderr!r}"
            ) from exc

        for job in jobs_payload.get("jobs", []):
            if job.get("conclusion") not in {
                "failure",
                "cancelled",
                "timed_out",
                "action_required",
            }:
                continue

            job_id = job.get("id")

            if job_id is None:
                continue

            log_proc = gh(
                "api",
                f"repos/{repo}/actions/jobs/{job_id}/logs",
                cwd=path,
                check=False,
            )

            raw_log = log_proc.stdout or log_proc.stderr
            truncated = len(raw_log) > max_chars

            if truncated:
                log = raw_log[-max_chars:]
            else:
                log = raw_log

            failures.append(
                {
                    "workflow": run_info.get("name"),
                    "run_id": run_id,
                    "run_url": run_info.get("html_url"),
                    "job_id": job_id,
                    "job_name": job.get("name"),
                    "job_conclusion": job.get("conclusion"),
                    "steps": [
                        {
                            "name": step.get("name"),
                            "status": step.get("status"),
                            "conclusion": step.get("conclusion"),
                            "number": step.get("number"),
                            "started_at": step.get("started_at"),
                            "completed_at": step.get("completed_at"),
                        }
                        for step in job.get("steps", [])
                    ],
                    "log": log,
                    "log_truncated": truncated,
                    "log_command_exit_code": log_proc.returncode,
                }
            )

    result = {
        "status": "observed",
        "task_id": normalize_task_id(task_id),
        "branch": branch,
        "head_sha": head_sha,
        "pr_number": pr["number"],
        "pr_url": pr["url"],
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
