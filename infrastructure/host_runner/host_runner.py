#!/usr/bin/env python3

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


WORKTREE_ROOT = Path.home() / "aiprojects" / "ion-mc-worktrees"
RUN_ROOT = Path.home() / ".local" / "share" / "ionmc-experiment" / "host-runs"
CACHE_ROOT = Path.home() / ".cache" / "ionmc-experiment" / "host-runner"
HOST_VENV = CACHE_ROOT / "venv"

MAX_TIMEOUT_SECONDS = 7200

SAFE_ENV = {
    "HOME": "/tmp/home",
    "XDG_CACHE_HOME": "/cache",
    "UV_CACHE_DIR": "/cache/uv",
    "WARP_CACHE_PATH": "/cache/warp",
    "NUMBA_CACHE_DIR": "/cache/numba",
    "CUDA_CACHE_PATH": "/cache/cuda",
    "PYTHONUNBUFFERED": "1",
    "VIRTUAL_ENV": "/runtime",
    "PATH": (
        "/runtime/bin"
        "/usr/lib/wsl/lib:"
        "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
    ),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def git_text(cwd: Path, *args: str) -> str:
    proc = subprocess.run(
        ["git", *args],
        cwd=cwd,
        text=True,
        capture_output=True,
        check=True,
    )
    return proc.stdout.strip()


def normalize_task_id(task_id: str) -> str:
    value = task_id.strip().upper()
    if not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{1,63}", value):
        raise ValueError(f"invalid task ID: {task_id!r}")
    return value


def task_worktree(task_id: str) -> Path:
    return WORKTREE_ROOT / normalize_task_id(task_id).lower()


def inspect_worktree(task_id: str) -> dict:
    path = task_worktree(task_id)

    if not path.is_dir():
        raise RuntimeError(f"task worktree does not exist: {path}")

    branch = git_text(path, "branch", "--show-current")
    sha = git_text(path, "rev-parse", "HEAD")
    status = git_text(path, "status", "--porcelain")

    if not branch.startswith("task/"):
        raise RuntimeError(f"unexpected task branch: {branch!r}")

    if status:
        raise RuntimeError(
            "host validation requires a clean committed task worktree"
        )

    return {
        "path": path,
        "branch": branch,
        "sha": sha,
    }


def ensure_runtime_dirs() -> None:
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)

    os.chmod(RUN_ROOT, 0o700)
    os.chmod(CACHE_ROOT, 0o700)


def create_run_dir(task_id: str) -> tuple[str, Path]:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    suffix = os.urandom(4).hex()
    run_id = f"RUN-{stamp}-{suffix}"

    path = RUN_ROOT / normalize_task_id(task_id) / run_id
    path.mkdir(parents=True, exist_ok=False)
    os.chmod(path, 0o700)

    return run_id, path


def add_optional_ro_bind(args: list[str], path: str) -> None:
    if Path(path).exists():
        args.extend(["--ro-bind", path, path])


def build_bwrap_command(worktree: Path, argv: list[str]) -> list[str]:
    if not argv:
        raise ValueError("argv must contain at least one argument")

    bwrap = shutil.which("bwrap")
    if not bwrap:
        raise RuntimeError("bubblewrap executable not found")

    if not Path("/dev/dxg").exists():
        raise RuntimeError("/dev/dxg is not available")

    args = [
        bwrap,
        "--die-with-parent",
        "--new-session",
        "--unshare-all",

        "--ro-bind", "/usr", "/usr",
        "--ro-bind", "/bin", "/bin",
        "--ro-bind", "/lib", "/lib",
    ]

    add_optional_ro_bind(args, "/lib64")
    add_optional_ro_bind(args, "/etc/ld.so.cache")
    add_optional_ro_bind(args, "/etc/ssl")
    add_optional_ro_bind(args, "/etc/ca-certificates")

    if not HOST_VENV.is_dir():
        raise RuntimeError(f"host runner virtual environment not found: {HOST_VENV}")

    args.extend([
        "--ro-bind", str(HOST_VENV), "/runtime",
    ])

    args.extend([
        "--ro-bind", "/sys", "/sys",
        "--proc", "/proc",
        "--dev", "/dev",
        "--dev-bind", "/dev/dxg", "/dev/dxg",

        "--bind", str(worktree), "/workspace",
        "--bind", str(CACHE_ROOT), "/cache",

        "--tmpfs", "/tmp",
        "--dir", "/tmp/home",

        "--chdir", "/workspace",

        "--clearenv",
    ])

    for key, value in SAFE_ENV.items():
        args.extend(["--setenv", key, value])

    args.extend(["--", *argv])
    return args


def write_json(path: Path, value: dict) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def run_validation(
    task_id: str,
    argv: list[str],
    timeout_seconds: int,
) -> dict:
    if timeout_seconds < 1 or timeout_seconds > MAX_TIMEOUT_SECONDS:
        raise ValueError(
            f"timeout_seconds must be between 1 and {MAX_TIMEOUT_SECONDS}"
        )

    ensure_runtime_dirs()
    worktree = inspect_worktree(task_id)

    run_id, run_dir = create_run_dir(task_id)

    request = {
        "schema_version": 1,
        "run_id": run_id,
        "task_id": normalize_task_id(task_id),
        "worktree": str(worktree["path"]),
        "branch": worktree["branch"],
        "sha": worktree["sha"],
        "argv": argv,
        "timeout_seconds": timeout_seconds,
        "requested_at": utc_now(),
    }

    write_json(run_dir / "request.json", request)

    bwrap_command = build_bwrap_command(worktree["path"], argv)

    started_wall = utc_now()
    started_mono = time.monotonic()

    timed_out = False

    try:
        proc = subprocess.run(
            bwrap_command,
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
        exit_code = proc.returncode
        stdout = proc.stdout
        stderr = proc.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        exit_code = None
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""

        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")

    duration = time.monotonic() - started_mono
    finished_wall = utc_now()

    (run_dir / "stdout.txt").write_text(stdout, encoding="utf-8")
    (run_dir / "stderr.txt").write_text(stderr, encoding="utf-8")

    result = {
        "schema_version": 1,
        "run_id": run_id,
        "task_id": normalize_task_id(task_id),
        "worktree": str(worktree["path"]),
        "branch": worktree["branch"],
        "sha": worktree["sha"],
        "dirty_before": False,
        "argv": argv,
        "started_at": started_wall,
        "finished_at": finished_wall,
        "duration_seconds": duration,
        "exit_code": exit_code,
        "timed_out": timed_out,
        "stdout_file": str(run_dir / "stdout.txt"),
        "stderr_file": str(run_dir / "stderr.txt"),
        "run_directory": str(run_dir),
        "succeeded": (not timed_out and exit_code == 0),
    }

    write_json(run_dir / "result.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-id", required=True)
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=600,
    )
    parser.add_argument(
        "argv",
        nargs=argparse.REMAINDER,
        help="command argv, preferably after --",
    )

    args = parser.parse_args()

    argv = args.argv
    if argv and argv[0] == "--":
        argv = argv[1:]

    try:
        result = run_validation(
            task_id=args.task_id,
            argv=argv,
            timeout_seconds=args.timeout_seconds,
        )
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

    print(json.dumps(result, indent=2))
    return 0 if result["succeeded"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
