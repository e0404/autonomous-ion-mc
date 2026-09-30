"""Bridge authoritative intervention records into v2 telemetry without their bodies."""

import sys

from infrastructure.experiment_v3.common import ROOT, event, git

if __name__ == "__main__":
    event(
        sys.argv[1],
        sha=git(ROOT, "rev-parse", "HEAD"),
        task_id=sys.argv[3],
        details={"request_id": sys.argv[2]},
    )
