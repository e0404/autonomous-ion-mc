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


def test_notices_packaged_and_printed(capsys: pytest.CaptureFixture[str]) -> None:
    from importlib.resources import files

    text = (files("ionmc") / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    assert "This product includes software developed by Members of the Geant4 Collaboration" in text
    assert "Copyright protection on this compilation of data has been secured" in text
    assert "Geant4 Software License" in text
    assert main(["notices"]) == 0
    assert "Geant4 Software License" in capsys.readouterr().out
