"""Permission notifications never grant approval or resolve interventions."""

from __future__ import annotations

import json
import sys
import uuid

from infrastructure.experiment_v3.common import event
from infrastructure.interventions.request_intervention import send_notification


def handle(payload, *, sender=send_notification, state=None):
    hook = payload.get("hook_event_name")
    if hook != "PermissionRequest" and not (
        hook == "Notification"
        and payload.get("notification_type") == "permission_prompt"
    ):
        return None
    record_id = "APPROVAL-" + uuid.uuid4().hex[:16]
    event(
        "permission_prompt", details={"record_id": record_id, "hook": hook}, state=state
    )
    result = sender(
        record_id,
        "runtime",
        "permission-approval",
        "IonMC runtime is waiting for explicit permission. "
        "Inspect the approval UI; this notification grants no "
        "permission.",
    )
    # No raw hook payload, URL, tool input or exception text is retained.
    event(
        "notification_result",
        details={"record_id": record_id, "status": result.get("status", "failed")},
        state=state,
    )
    return result.get("status", "failed")


def main():
    try:
        status = handle(json.load(sys.stdin))
        if status not in (None, "sent"):
            print(
                "IonMC permission notification was not delivered; "
                "inspect local notification configuration.",
                file=sys.stderr,
            )
    except Exception:
        print(
            "IonMC permission notification failed; explicit approval still required.",
            file=sys.stderr,
        )
    # No approval output: Claude retains normal permission handling.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
