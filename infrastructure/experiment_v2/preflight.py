"""Fail-closed unattended-launch checks; no scientific implementation is started."""

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

from infrastructure.experiment_v2.common import (
    ROOT,
    STATE,
    exact_state,
    write_json,
)
from infrastructure.experiment_v2.reference import load_engine, run_reference
from infrastructure.interventions.request_intervention import send_notification


def check(root=ROOT, *, smoke=False, notify=False, state=STATE):
    root, state = Path(root), Path(state)
    results = {}
    try:
        results["clean_sha"] = exact_state(root)
        condition = json.loads((root / ".ionmc-condition.json").read_text())
        results["fresh_condition"] = (
            condition["experiment_id"] == "experiment-v2"
            and condition["scientific_implementation_inherited"] is False
        )
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError):
        results["fresh_condition"] = False
    results["tools"] = {
        name: shutil.which(name) is not None
        for name in ("bwrap", "git", "gh", "uv", "claude")
    }
    results["gpu_device"] = Path("/dev/dxg").exists()
    try:
        from infrastructure.host_runner.host_runner import HOST_VENV

        results["host_runtime"] = (HOST_VENV / "bin/python").exists()
    except ImportError:
        results["host_runtime"] = False
    results["engines"] = {}
    for engine in ("topas", "mcsquare", "fred"):
        try:
            config = load_engine(engine, state / "engines.json")
            results["engines"][engine] = {
                "available": True,
                "version": config["version"],
            }
            if smoke:
                run = run_reference(
                    root,
                    "experiment/v2/references/" + engine,
                    engine,
                    registry=state / "engines.json",
                    state=state,
                    timeout=120,
                )
                results["engines"][engine]["smoke"] = run["succeeded"]
                results["engines"][engine]["run_id"] = run["run_id"]
        except (ValueError, OSError, KeyError) as exc:
            results["engines"][engine] = {
                "available": False,
                "error_type": type(exc).__name__,
            }
    try:
        from infrastructure.experiment_v2.remote import inspect

        results["remote"] = inspect(root)
    except (ValueError, OSError, KeyError, IndexError, subprocess.CalledProcessError):
        results["remote"] = {"ready": False}
    results["notification"] = {
        "configured": bool(os.environ.get("IONMC_NTFY_URL")),
        "delivery_tested": False,
    }
    if notify:
        delivery = send_notification(
            "V2-PREFLIGHT",
            "setup",
            "notification-test",
            "IonMC v2 operator notification preflight. No approval "
            "or scientific intervention requested.",
        )
        results["notification"].update(delivery_tested=True, status=delivery["status"])
    # This preflight deliberately never labels an untested notification/engine as ready.
    results["ready_for_unattended_launch"] = bool(
        results.get("fresh_condition")
        and all(results["tools"].values())
        and results["gpu_device"]
        and results["host_runtime"]
        and results["remote"]["ready"]
        and all(
            e.get("available") and e.get("smoke") for e in results["engines"].values()
        )
        and results["notification"].get("status") == "sent"
    )
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--notify", action="store_true")
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    result = check(smoke=a.smoke, notify=a.notify)
    if a.output:
        write_json(a.output, result)
    print(json.dumps(result, indent=2))
    return 0 if result["ready_for_unattended_launch"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
