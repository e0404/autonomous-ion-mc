"""Tests for environment description."""

import json

import ionmc
from ionmc.environment import describe_environment


def test_describe_environment_fields() -> None:
    env = describe_environment()
    assert env["ionmc_version"] == ionmc.__version__
    assert env["python_version"]
    assert env["numpy_version"].startswith("2")
    assert "git_sha" in env and "git_dirty" in env
    json.dumps(env)


def test_warp_cpu_device_listed() -> None:
    env = describe_environment()
    assert env["warp_error"] is None
    assert env["warp_version"]
    cpu = [d for d in env["warp_devices"] if not d["is_cuda"]]
    assert cpu and cpu[0]["alias"] == "cpu"
