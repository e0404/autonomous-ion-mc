#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path


RAW_ROOT = Path(
    os.environ.get(
        "IONMC_CODEX_RAW_DIR",
        str(Path.home() / ".local/share/ionmc-experiment/raw/codex"),
    )
)

EXPERIMENT_ID = os.environ.get("IONMC_EXPERIMENT_ID", "experiment-v1")


def run_git(cwd: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args],
        cwd=cwd,
        text=True,
        stderr=subprocess.PIPE,
    ).strip()


def git_info(cwd: Path) -> dict:
    return {
        "branch": run_git(cwd, "branch", "--show-current"),
        "sha": run_git(cwd, "rev-parse", "HEAD"),
        "status": run_git(cwd, "status", "--porcelain"),
    }


def parse_codex_jsonl(path: Path) -> dict:
    thread_id = None
    final_message = None
    usage = None
    event_count = 0

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue

            event_count += 1
            event = json.loads(line)

            if event.get("type") == "thread.started":
                thread_id = event.get("thread_id")

            elif event.get("type") == "item.completed":
                item = event.get("item", {})
                if item.get("type") == "agent_message":
                    final_message = item.get("text")

            elif event.get("type") == "turn.completed":
                usage = event.get("usage")

    return {
        "thread_id": thread_id,
        "final_message": final_message,
        "usage": usage,
        "event_count": event_count,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Controlled Codex worker for the autonomous ion-MC experiment."
    )

    parser.add_argument(
        "--worktree",
        type=Path,
        required=True,
        help="Existing isolated Git worktree in which Codex executes.",
    )

    parser.add_argument(
        "--task-id",
        required=True,
        help="Experiment task identifier.",
    )

    parser.add_argument(
        "--role",
        default="external-codex-worker",
        help="Semantic worker role recorded in telemetry.",
    )

    parser.add_argument(
        "--model",
        default=None,
        help="Optional explicit Codex model.",
    )

    parser.add_argument(
        "--prompt-file",
        type=Path,
        help="Read the worker prompt from this UTF-8 file.",
    )

    parser.add_argument(
        "--prompt",
        help="Worker prompt supplied directly on the command line.",
    )

    args = parser.parse_args()

    if bool(args.prompt) == bool(args.prompt_file):
        parser.error("Specify exactly one of --prompt or --prompt-file.")

    worktree = args.worktree.resolve()

    if not (worktree / ".git").exists():
        print(f"Not a Git worktree: {worktree}", file=sys.stderr)
        return 2

    prompt = (
        args.prompt
        if args.prompt is not None
        else args.prompt_file.read_text(encoding="utf-8")
    )

    RAW_ROOT.mkdir(parents=True, exist_ok=True)

    run_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "-"
        + uuid.uuid4().hex[:8]
    )

    run_dir = RAW_ROOT / args.task_id / run_id
    run_dir.mkdir(parents=True, exist_ok=False)

    prompt_path = run_dir / "prompt.txt"
    events_path = run_dir / "events.jsonl"
    metadata_path = run_dir / "metadata.json"
    final_path = run_dir / "final-message.txt"

    prompt_path.write_text(prompt, encoding="utf-8")

    before = git_info(worktree)

    command = [
        "codex",
        "exec",
        "--json",
        "--sandbox",
        "workspace-write",
    ]

    if args.model:
        command += ["--model", args.model]

    command.append(prompt)

    started_at = datetime.now(timezone.utc)
    started_monotonic = time.monotonic()

    with events_path.open("w", encoding="utf-8") as output:
        proc = subprocess.run(
            command,
            cwd=worktree,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
            env=os.environ.copy(),
        )

    ended_at = datetime.now(timezone.utc)
    elapsed_seconds = time.monotonic() - started_monotonic

    parsed = parse_codex_jsonl(events_path)
    after = git_info(worktree)

    if parsed["final_message"] is not None:
        final_path.write_text(parsed["final_message"] + "\n", encoding="utf-8")

    metadata = {
        "schema_version": 1,
        "experiment_id": EXPERIMENT_ID,
        "run_id": run_id,
        "task_id": args.task_id,
        "role": args.role,
        "provider": "openai",
        "tool": "codex-cli",
        "worktree": str(worktree),
        "requested_model": args.model,
        "started_at_utc": started_at.isoformat(),
        "ended_at_utc": ended_at.isoformat(),
        "elapsed_seconds": elapsed_seconds,
        "exit_code": proc.returncode,
        "stderr": proc.stderr,
        "git_before": before,
        "git_after": after,
        "codex": parsed,
        "files": {
            "prompt": str(prompt_path),
            "events": str(events_path),
            "final_message": (
                str(final_path)
                if parsed["final_message"] is not None
                else None
            ),
        },
    }

    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    # stdout is intentionally minimal because this is what the calling
    # orchestrator will eventually receive.
    if proc.returncode != 0:
        print(
            json.dumps(
                {
                    "status": "failed",
                    "task_id": args.task_id,
                    "run_id": run_id,
                    "exit_code": proc.returncode,
                    "metadata": str(metadata_path),
                }
            )
        )
        return proc.returncode

    print(
        json.dumps(
            {
                "status": "completed",
                "task_id": args.task_id,
                "run_id": run_id,
                "thread_id": parsed["thread_id"],
                "usage": parsed["usage"],
                "git_before": before["sha"],
                "git_after": after["sha"],
                "dirty": bool(after["status"]),
                "final_message": parsed["final_message"],
                "metadata": str(metadata_path),
            },
            ensure_ascii=False,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
