#!/usr/bin/env python3

import argparse
import json
import os
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path


EXPERIMENT_ID = os.environ.get("IONMC_EXPERIMENT_ID", "experiment-v1")

TELEMETRY_DIR = Path(
    os.environ.get(
        "IONMC_TELEMETRY_DIR",
        str(Path.home() / ".local/share/ionmc-experiment/telemetry"),
    )
)

OUTPUT = TELEMETRY_DIR / "experiment-events.jsonl"


def git(*args: str) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *args],
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=2,
        ).strip()
    except Exception:
        return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--event", required=True)
    parser.add_argument("--task-id")
    parser.add_argument("--data-json", default="{}")
    args = parser.parse_args()

    data = json.loads(args.data_json)

    record = {
        "schema_version": 1,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "experiment_id": EXPERIMENT_ID,
        "event": args.event,
        "task_id": args.task_id,
        "hostname": socket.gethostname(),
        "cwd": str(Path.cwd()),
        "git_branch": git("branch", "--show-current"),
        "git_sha": git("rev-parse", "HEAD"),
        "data": data,
    }

    TELEMETRY_DIR.mkdir(parents=True, exist_ok=True)

    with OUTPUT.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")))
        f.write("\n")

    print(json.dumps(record, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
