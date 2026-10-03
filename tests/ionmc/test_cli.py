"""Tests for the ionmc command-line interface."""

import json

import pytest

import ionmc
from ionmc.cli import main


def test_version_subcommand(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["version"]) == 0
    assert capsys.readouterr().out.strip() == ionmc.__version__


def test_version_flag(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--version"]) == 0
    assert ionmc.__version__ in capsys.readouterr().out


def test_info_prints_json(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["info"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["ionmc_version"] == ionmc.__version__


def test_no_command_is_usage_error() -> None:
    assert main([]) == 2


def test_unknown_command_is_usage_error() -> None:
    assert main(["bogus"]) == 2
