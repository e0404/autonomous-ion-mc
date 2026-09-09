#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import urllib.error
import urllib.request
import uuid
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REQUEST_ROOT = ROOT / "reference_requests" / "interventions"
EVENT_RECORDER = ROOT / "infrastructure" / "telemetry" / "record_event.py"


def send_notification(
    request_id: str,
    task_id: str,
    category: str,
    summary: str,
) -> dict:
    url = os.environ.get("IONMC_NTFY_URL")

    if not url:
        return {
            "status": "not-configured",
            "detail": "IONMC_NTFY_URL is not set",
        }

    message = (
        f"{summary}\n\n"
        f"Request: {request_id}\n"
        f"Task: {task_id}\n"
        f"Category: {category}"
    )

    req = urllib.request.Request(
        url,
        data=message.encode("utf-8"),
        method="POST",
        headers={
            "Title": "IonMC: human input required",
            "Priority": "high",
            "Tags": "warning,robot",
            "Content-Type": "text/plain; charset=utf-8",
        },
    )

    try:
        with urllib.request.urlopen(req, timeout=15) as response:
            return {
                "status": "sent",
                "http_status": response.status,
            }
    except (urllib.error.URLError, TimeoutError) as exc:
        return {
            "status": "failed",
            "error": str(exc),
        }


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("--task-id", required=True)

    parser.add_argument(
        "--category",
        required=True,
        choices=[
            "inaccessible-information",
            "external-reference-calculation",
            "external-data",
            "authorization",
            "experiment-boundary",
        ],
    )

    parser.add_argument("--summary", required=True)
    parser.add_argument("--reason", required=True)
    parser.add_argument("--requested-input", required=True)
    parser.add_argument("--depends-on", default="")

    args = parser.parse_args()

    now = datetime.now(timezone.utc)

    request_id = (
        f"IR-{now.strftime('%Y%m%d-%H%M%S')}-"
        f"{uuid.uuid4().hex[:6].upper()}"
    )

    request_dir = REQUEST_ROOT / request_id
    request_dir.mkdir(parents=True, exist_ok=False)

    notification = send_notification(
        request_id=request_id,
        task_id=args.task_id,
        category=args.category,
        summary=args.summary,
    )

    payload = {
        "schema_version": 1,
        "request_id": request_id,
        "created_at_utc": now.isoformat(),
        "task_id": args.task_id,
        "category": args.category,
        "summary": args.summary,
        "reason_autonomous_resolution_not_possible": args.reason,
        "requested_input": args.requested_input,
        "decisions_or_tasks_blocked": args.depends_on or None,
        "status": "open",
        "notification": notification,
    }

    (request_dir / "request.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    readme = f"""# Human Intervention Request {request_id}

**Task:** {args.task_id}
**Category:** {args.category}
**Status:** open

## Summary

{args.summary}

## Why autonomous resolution is not possible

{args.reason}

## Exact requested input or action

{args.requested_input}

## Decision or task depending on this request

{args.depends_on or "Not specified."}
"""

    (request_dir / "README.md").write_text(
        readme,
        encoding="utf-8",
    )

    subprocess.run(
        [
            str(EVENT_RECORDER),
            "--event",
            "human_intervention_requested",
            "--task-id",
            args.task_id,
            "--data-json",
            json.dumps(
                {
                    "request_id": request_id,
                    "category": args.category,
                    "summary": args.summary,
                    "notification_status": notification["status"],
                }
            ),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    print(
        json.dumps(
            {
                "status": "created",
                "request_id": request_id,
                "request_directory": str(request_dir),
                "notification": notification,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
