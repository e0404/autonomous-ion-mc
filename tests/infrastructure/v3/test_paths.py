import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from infrastructure.experiment_v3.paths import storage_roots


@pytest.mark.parametrize("experiment_id", ["experiment-v3", "experiment-v4"])
def test_storage_is_partitioned_by_experiment(monkeypatch, tmp_path, experiment_id):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    cache, state = storage_roots(experiment_id)
    assert cache == tmp_path / ".cache/ionmc-experiment" / experiment_id
    assert state == tmp_path / ".local/share/ionmc-experiment" / experiment_id


@pytest.mark.parametrize("experiment_id", ["", "../v1", "/tmp/v2", "a/b", "a b"])
def test_storage_rejects_directory_escape(experiment_id):
    with pytest.raises(ValueError):
        storage_roots(experiment_id)


def test_project_settings_and_runner_agree_on_storage():
    root = Path(__file__).resolve().parents[3]
    cfg = json.loads((root / ".claude/settings.json").read_text())
    env = dict(os.environ, **cfg["env"])
    env.pop("IONMC_V3_STATE", None)
    result = subprocess.check_output(
        [
            sys.executable,
            "-c",
            "import json; from infrastructure.experiment_v3.common import STATE; "
            "from infrastructure.host_runner.host_runner import CACHE_ROOT,RUN_ROOT; "
            "print(json.dumps([str(STATE),str(CACHE_ROOT),str(RUN_ROOT)]))",
        ],
        cwd=root,
        env=env,
        text=True,
    )
    state, cache, runs = map(Path, json.loads(result))
    expected_cache, expected_state = storage_roots(cfg["env"]["IONMC_EXPERIMENT_ID"])
    assert state == expected_state
    assert cache == expected_cache / "host-runner"
    assert runs == state / "host-runs"
    for key, suffix in [
        ("IONMC_TELEMETRY_DIR", "telemetry"),
        ("IONMC_V3_EVENT_DIR", "telemetry"),
        ("IONMC_CODEX_RAW_DIR", "raw/codex"),
        ("IONMC_VALIDATION_DIR", "validation"),
    ]:
        assert Path(cfg["env"][key]).expanduser() == state / suffix
    for key in ["allowRead", "allowWrite"]:
        allowed = cfg["sandbox"]["filesystem"][key]
        assert "~/.cache/ionmc-experiment" not in allowed
        assert str(expected_cache).replace(str(Path.home()), "~", 1) in allowed


def test_state_override_expands_home():
    root = Path(__file__).resolve().parents[3]
    env = dict(os.environ, IONMC_V3_STATE="~/ionmc-test-state")
    result = subprocess.check_output(
        [
            sys.executable,
            "-c",
            "from infrastructure.experiment_v3.common import STATE; print(STATE)",
        ],
        cwd=root,
        env=env,
        text=True,
    ).strip()
    assert Path(result) == Path.home() / "ionmc-test-state"
