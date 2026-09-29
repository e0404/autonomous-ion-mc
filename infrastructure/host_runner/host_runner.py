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


WORKTREE_ROOT = Path(__file__).resolve().parents[2].parent / (Path(__file__).resolve().parents[2].name + "-worktrees")
IS_V2 = (Path(__file__).resolve().parents[2] / ".ionmc-condition.json").exists()
if IS_V2:
    # The runner is also invoked as a script by the trusted MCP server.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from infrastructure.experiment_v2.paths import CACHE_ROOT as EXPERIMENT_CACHE
    from infrastructure.experiment_v2.paths import SHARE_ROOT

    RUN_ROOT = SHARE_ROOT / "host-runs"
    CACHE_ROOT = EXPERIMENT_CACHE / "host-runner"
else:
    RUN_ROOT = Path.home() / ".local/share/ionmc-experiment/host-runs"
    CACHE_ROOT = Path.home() / ".cache/ionmc-experiment/host-runner"
HOST_VENV = CACHE_ROOT / "venv"
SANDBOX_VENV = str(HOST_VENV)

MAX_TIMEOUT_SECONDS = 7200
MAX_INLINE_OUTPUT_CHARS = 65536

SAFE_ENV = {
    "HOME": "/tmp/home",
    "XDG_CACHE_HOME": "/cache",
    "UV_CACHE_DIR": "/cache/uv",
    "WARP_CACHE_PATH": "/cache/warp",
    "NUMBA_CACHE_DIR": "/cache/numba",
    "CUDA_CACHE_PATH": "/cache/cuda",
    "PYTHONUNBUFFERED": "1",
    "VIRTUAL_ENV": str(HOST_VENV),
    "PATH": (
        f"{HOST_VENV}/bin:"
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
        "git_common_dir": git_common_dir(path),
    }


def git_common_dir(worktree: Path) -> Path:
    """Absolute Git common directory backing a (possibly linked) worktree."""
    return Path(
        git_text(worktree, "rev-parse", "--path-format=absolute", "--git-common-dir")
    ).resolve()


def git_dir_binds(worktree: Path, common_dir: Path | None) -> list[str]:
    """Read-only bind of a linked worktree's Git metadata, at its host path.

    A linked worktree's ``.git`` file points at ``<main>/.git/worktrees/<id>``
    outside the workspace mount. Without that directory, code inside the
    sandbox cannot resolve its own commit SHA or dirty state, so persisted
    scientific results would lose their provenance. The primary checkout keeps
    ``.git`` inside the workspace and needs no extra mount.
    """
    if common_dir is None:
        return []
    common_dir = Path(common_dir).resolve()
    if common_dir.is_relative_to(Path(worktree).resolve()):
        return []
    return ["--ro-bind", str(common_dir), str(common_dir)]


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


def build_bwrap_command(
    worktree: Path, argv: list[str], git_common: Path | None = None
) -> list[str]:
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
    # Host identification for reproducibility records (no secrets inside).
    add_optional_ro_bind(args, "/etc/os-release")

    if not HOST_VENV.is_dir():
        raise RuntimeError(f"host runner virtual environment not found: {HOST_VENV}")

    args.extend([
        "--dir", "/home",
        "--dir", str(Path.home()),
        "--dir", str(Path.home() / ".cache"),
        "--dir", str(Path.home() / ".cache" / "ionmc-experiment"),
        "--dir", str(CACHE_ROOT.parent),
        "--dir", str(CACHE_ROOT),
        "--ro-bind", str(HOST_VENV), str(HOST_VENV),
    ])

    args.extend([
        "--ro-bind", "/sys", "/sys",
        "--proc", "/proc",
        "--dev", "/dev",
        "--dev-bind", "/dev/dxg", "/dev/dxg",

        "--bind", str(worktree), "/workspace",
        *git_dir_binds(worktree, git_common),
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


def bounded_tail(text: str, max_chars: int) -> tuple[str, bool, int]:
    total_chars = len(text)

    if total_chars <= max_chars:
        return text, False, total_chars

    return text[-max_chars:], True, total_chars


def run_validation(
    task_id: str,
    argv: list[str],
    timeout_seconds: int,
    *,
    release_root: Path | None = None,
) -> dict:
    if timeout_seconds < 1 or timeout_seconds > MAX_TIMEOUT_SECONDS:
        raise ValueError(
            f"timeout_seconds must be between 1 and {MAX_TIMEOUT_SECONDS}"
        )

    ensure_runtime_dirs()
    if release_root is None:
        worktree = inspect_worktree(task_id)
    else:
        condition = json.loads((release_root / ".ionmc-condition.json").read_text())
        if (condition.get("experiment_id") != "experiment-v2"
                or git_text(release_root, "branch", "--show-current") != condition["integration_branch"]
                or git_text(release_root, "status", "--porcelain")):
            raise RuntimeError("Release validation requires a clean v2 integration checkout")
        worktree = {"path": release_root, "branch": condition["integration_branch"],
                    "sha": git_text(release_root, "rev-parse", "HEAD"),
                    "git_common_dir": git_common_dir(release_root)}


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

    bwrap_command = build_bwrap_command(
        worktree["path"], argv, worktree.get("git_common_dir")
    )

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

    stdout_inline, stdout_truncated, stdout_total_chars = bounded_tail(
        stdout,
        MAX_INLINE_OUTPUT_CHARS,
    )

    stderr_inline, stderr_truncated, stderr_total_chars = bounded_tail(
        stderr,
        MAX_INLINE_OUTPUT_CHARS,
    )

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
        "stdout": stdout_inline,
        "stderr": stderr_inline,
        "stdout_truncated": stdout_truncated,
        "stderr_truncated": stderr_truncated,
        "stdout_total_chars": stdout_total_chars,
        "stderr_total_chars": stderr_total_chars,
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
