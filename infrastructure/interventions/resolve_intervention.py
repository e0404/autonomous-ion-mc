#!/usr/bin/env python3

import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REQUEST_ROOT = ROOT / "reference_requests" / "interventions"
EVENT_RECORDER = ROOT / "infrastructure" / "telemetry" / "record_event.py"


def atomic_write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    tmp.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request-id", required=True)
    parser.add_argument("--response", required=True)
    parser.add_argument("--operator-note", default="")
    args = parser.parse_args()

    request_dir = REQUEST_ROOT / args.request_id
    request_file = request_dir / "request.json"
    response_file = request_dir / "response.json"

    if not request_file.is_file():
        raise SystemExit(f"Unknown intervention request: {args.request_id}")

    request = json.loads(request_file.read_text(encoding="utf-8"))

    if request.get("status") != "open":
        raise SystemExit(
            f"Request {args.request_id} is not open "
            f"(status={request.get('status')!r})"
        )

    if response_file.exists():
        raise SystemExit(
            f"Response already exists for request {args.request_id}"
        )

    now = datetime.now(timezone.utc)

    response = {
        "schema_version": 1,
        "request_id": args.request_id,
        "resolved_at_utc": now.isoformat(),
        "task_id": request.get("task_id"),
        "response": args.response,
        "operator_note": args.operator_note or None,
    }

    atomic_write_json(response_file, response)

    request["status"] = "resolved"
    request["resolved_at_utc"] = now.isoformat()
    request["response_file"] = "response.json"

    atomic_write_json(request_file, request)

    subprocess.run(
        [
            str(EVENT_RECORDER),
            "--event",
            "human_intervention_resolved",
            "--task-id",
            request.get("task_id") or "",
            "--data-json",
            json.dumps(
                {
                    "request_id": args.request_id,
                    "category": request.get("category"),
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
                "status": "resolved",
                "request_id": args.request_id,
                "request_directory": str(request_dir),
                "response_file": str(response_file),
            },
            indent=2,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
