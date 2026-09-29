import json
import subprocess
import sys
from pathlib import Path

import ionmc
from ionmc import cli, provenance


def test_version_is_exposed():
    assert ionmc.__version__ == "0.1.0.dev0"


def test_identity_from_git_checkout_matches_git():
    identity = ionmc.code_identity()
    root = provenance._repository_root(Path(ionmc.__file__).resolve().parent)
    assert root is not None
    sha = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    assert identity["source"] == "git"
    assert identity["git_sha"] == sha
    assert isinstance(identity["git_dirty"], bool)


def test_unknown_identity_is_explicit(monkeypatch):
    monkeypatch.setattr(provenance, "_repository_root", lambda _start: None)
    monkeypatch.setattr(provenance, "_from_build_info", lambda: None)
    identity = provenance.code_identity()
    assert identity["source"] == "unknown"
    assert identity["git_sha"] is None and identity["git_dirty"] is None


def test_build_info_is_reported_when_git_is_unavailable(monkeypatch):
    monkeypatch.setattr(provenance, "_repository_root", lambda _start: None)
    monkeypatch.setattr(
        provenance,
        "_from_build_info",
        lambda: {
            "version": "0.1.0.dev0",
            "git_sha": "a" * 40,
            "git_dirty": True,
            "source": "wheel:build-time-git",
        },
    )
    identity = provenance.code_identity()
    assert identity["git_dirty"] is True and identity["git_sha"] == "a" * 40


def test_cli_version_prints_json(capsys):
    assert cli.main(["version"]) == 0
    out = json.loads(capsys.readouterr().out)
    assert set(out) == {"version", "git_sha", "git_dirty", "source"}


def test_module_entry_point_runs():
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            "from ionmc.cli import main; raise SystemExit(main(['version']))",
        ],
        capture_output=True,
        text=True,
        cwd=Path(ionmc.__file__).resolve().parents[2],
        env={"PYTHONPATH": str(Path(ionmc.__file__).resolve().parents[1])},
    )
    assert proc.returncode == 0, proc.stderr
    assert "git_sha" in proc.stdout
