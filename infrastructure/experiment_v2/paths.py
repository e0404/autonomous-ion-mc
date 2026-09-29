"""Experiment-ID-based storage roots shared by preflight and host execution."""

import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def storage_roots(experiment_id):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", experiment_id):
        raise ValueError("Experiment ID must be a single safe directory name")
    return (
        Path.home() / ".cache/ionmc-experiment" / experiment_id,
        Path.home() / ".local/share/ionmc-experiment" / experiment_id,
    )


condition = ROOT / ".ionmc-condition.json"
default_id = (
    json.loads(condition.read_text())["experiment_id"]
    if condition.exists()
    else "experiment-v2"
)
EXPERIMENT_ID = os.environ.get("IONMC_EXPERIMENT_ID", default_id)
CACHE_ROOT, SHARE_ROOT = storage_roots(EXPERIMENT_ID)
