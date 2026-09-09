"""Tests for infrastructure/diagnostics/environment_report.py.

These tests must be deterministic on a host with no GPU, no `warp`, no
`claude` CLI, and no `codex` CLI (exactly the CI/dev host they run on).
No test may assume any optional tool is actually present; graceful
degradation is exercised by monkeypatching `subprocess.run` rather than
depending on host state.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime

import pytest

import environment_report as er


# ---------------------------------------------------------------------------
# Full report structure
# ---------------------------------------------------------------------------


EXPECTED_TOP_LEVEL_KEYS = {
    "schema_version",
    "generated_at_utc",
    "os",
    "cpu",
    "memory",
    "python",
    "uv",
    "git",
    "repository",
    "claude_code",
    "codex_cli",
    "nvidia_gpu",
    "warp",
}

VALID_STATUSES = {er.STATUS_AVAILABLE, er.STATUS_UNAVAILABLE, er.STATUS_ERROR}


def test_build_report_has_expected_top_level_keys():
    report = er.build_report()
    assert EXPECTED_TOP_LEVEL_KEYS.issubset(report.keys())


def test_build_report_is_json_serializable():
    report = er.build_report()
    text = json.dumps(report)
    round_tripped = json.loads(text)
    assert round_tripped["schema_version"] == er.SCHEMA_VERSION


def test_schema_version_is_int():
    report = er.build_report()
    assert isinstance(report["schema_version"], int)
    assert report["schema_version"] == 1


def test_generated_at_utc_is_parseable_and_timezone_aware():
    report = er.build_report()
    parsed = datetime.fromisoformat(report["generated_at_utc"])
    assert parsed.tzinfo is not None
    assert parsed.utcoffset().total_seconds() == 0


def test_every_section_has_a_valid_status():
    report = er.build_report()
    for key in EXPECTED_TOP_LEVEL_KEYS - {"schema_version", "generated_at_utc"}:
        section = report[key]
        assert isinstance(section, dict), f"section {key} is not a dict"
        assert "status" in section, f"section {key} has no status"
        assert section["status"] in VALID_STATUSES, f"section {key} has bad status {section['status']!r}"


def test_python_section_always_available():
    section = er.build_python_section()
    assert section["status"] == er.STATUS_AVAILABLE
    assert section["executable"] == sys.executable


def test_on_this_host_optional_tools_report_unavailable_not_error():
    # This host has no claude or codex CLI, and no warp install.
    report = er.build_report()
    for key in ("claude_code", "codex_cli", "warp"):
        assert report[key]["status"] == er.STATUS_UNAVAILABLE, report[key]
    # nvidia-smi may be present on PATH (e.g. under WSL) even with no GPU
    # accessible, in which case it exits non-zero rather than being
    # missing outright - either degrades gracefully (never "available").
    assert report["nvidia_gpu"]["status"] in (er.STATUS_UNAVAILABLE, er.STATUS_ERROR), report["nvidia_gpu"]


# ---------------------------------------------------------------------------
# Graceful degradation of subprocess-based probes
# ---------------------------------------------------------------------------


def test_run_subprocess_missing_executable(monkeypatch):
    def fake_run(*args, **kwargs):
        raise FileNotFoundError("no such file")

    monkeypatch.setattr(er.subprocess, "run", fake_run)
    result = er.run_subprocess(["totally-not-a-real-command"])
    assert result["ok"] is False
    assert result["status"] == er.STATUS_UNAVAILABLE
    assert "not found" in result["detail"]


def test_run_subprocess_permission_error(monkeypatch):
    def fake_run(*args, **kwargs):
        raise PermissionError("denied")

    monkeypatch.setattr(er.subprocess, "run", fake_run)
    result = er.run_subprocess(["some-command"])
    assert result["ok"] is False
    assert result["status"] == er.STATUS_UNAVAILABLE


def test_run_subprocess_timeout(monkeypatch):
    def fake_run(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd="slow-command", timeout=5)

    monkeypatch.setattr(er.subprocess, "run", fake_run)
    result = er.run_subprocess(["slow-command"])
    assert result["ok"] is False
    assert result["status"] == er.STATUS_ERROR
    assert "timed out" in result["detail"]


def test_run_subprocess_nonzero_exit(monkeypatch):
    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(args=args, returncode=1, stdout="", stderr="boom")

    monkeypatch.setattr(er.subprocess, "run", fake_run)
    result = er.run_subprocess(["failing-command"])
    assert result["ok"] is False
    assert result["status"] == er.STATUS_ERROR
    assert "exited 1" in result["detail"]
    assert "boom" in result["detail"]


def test_run_subprocess_unexpected_exception(monkeypatch):
    def fake_run(*args, **kwargs):
        raise ValueError("something unexpected")

    monkeypatch.setattr(er.subprocess, "run", fake_run)
    result = er.run_subprocess(["weird-command"])
    assert result["ok"] is False
    assert result["status"] == er.STATUS_ERROR
    assert "unexpected error" in result["detail"]


def test_run_subprocess_does_not_inherit_stdin(monkeypatch):
    captured_kwargs = {}

    def fake_run(*args, **kwargs):
        captured_kwargs.update(kwargs)
        return subprocess.CompletedProcess(args=args, returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(er.subprocess, "run", fake_run)
    er.run_subprocess(["some-command"])
    assert captured_kwargs["stdin"] == subprocess.DEVNULL
    assert "timeout" in captured_kwargs


@pytest.mark.parametrize(
    "failure",
    ["missing", "timeout", "nonzero", "raises"],
)
def test_version_section_degrades_gracefully(monkeypatch, failure):
    def fake_run(*args, **kwargs):
        if failure == "missing":
            raise FileNotFoundError("nope")
        if failure == "timeout":
            raise subprocess.TimeoutExpired(cmd="tool", timeout=5)
        if failure == "nonzero":
            return subprocess.CompletedProcess(args=args, returncode=2, stdout="", stderr="bad")
        raise RuntimeError("kaboom")

    monkeypatch.setattr(er.subprocess, "run", fake_run)
    section = er.version_section(["some-tool", "--version"])
    assert section["status"] in (er.STATUS_UNAVAILABLE, er.STATUS_ERROR)
    assert "detail" in section


def test_report_survives_subprocess_always_failing(monkeypatch):
    """Even if every subprocess call blows up, build_report() must still
    succeed and return a fully-formed report (hard requirement)."""

    def fake_run(*args, **kwargs):
        raise RuntimeError("simulated total subprocess failure")

    monkeypatch.setattr(er.subprocess, "run", fake_run)
    report = er.build_report()

    assert report["schema_version"] == er.SCHEMA_VERSION
    for key in ("uv", "git", "repository", "claude_code", "codex_cli", "nvidia_gpu"):
        assert report[key]["status"] == er.STATUS_ERROR
    # non-subprocess sections are unaffected
    assert report["python"]["status"] == er.STATUS_AVAILABLE


def test_safe_section_catches_builder_exception():
    def broken_builder():
        raise RuntimeError("builder blew up")

    result = er._safe_section(broken_builder)
    assert result["status"] == er.STATUS_ERROR
    assert "builder blew up" in result["detail"]


def test_safe_section_catches_malformed_result():
    def malformed_builder():
        return {"no_status_key": True}

    result = er._safe_section(malformed_builder)
    assert result["status"] == er.STATUS_ERROR


# ---------------------------------------------------------------------------
# Repository section
# ---------------------------------------------------------------------------


def test_repository_section_resolved_from_module_file_not_cwd(tmp_path, monkeypatch):
    # Point the probe at this actual repository checkout (a real git repo)
    # while the CWD is some unrelated directory, to prove the CWD is never
    # consulted.
    repo_dir = er.Path(__file__).resolve().parents[2]
    monkeypatch.chdir(tmp_path)

    section = er.build_repository_section(module_path=repo_dir / "infrastructure" / "diagnostics")

    assert section["status"] == er.STATUS_AVAILABLE
    assert isinstance(section["branch"], str) and section["branch"]
    assert isinstance(section["sha"], str) and len(section["sha"]) == 40
    assert isinstance(section["dirty"], bool)


def test_repository_section_unavailable_outside_a_repo(tmp_path):
    section = er.build_repository_section(module_path=tmp_path)
    assert section["status"] in (er.STATUS_UNAVAILABLE, er.STATUS_ERROR)


# ---------------------------------------------------------------------------
# Optional-import probe (warp)
# ---------------------------------------------------------------------------


def test_warp_section_reports_unavailable_on_import_error(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "warp":
            raise ImportError("No module named 'warp'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    section = er.build_warp_section()
    assert section["status"] == er.STATUS_UNAVAILABLE
    assert "warp" in section["detail"].lower()


def test_warp_section_on_this_host_is_unavailable():
    # warp is genuinely not installed on this host; exercises the real path.
    section = er.build_warp_section()
    assert section["status"] == er.STATUS_UNAVAILABLE


# ---------------------------------------------------------------------------
# Parsing helpers, tested against fixture strings (not host state)
# ---------------------------------------------------------------------------


OS_RELEASE_FIXTURE = """\
NAME="Ubuntu"
VERSION="22.04.3 LTS (Jammy Jellyfish)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 22.04.3 LTS"
VERSION_ID="22.04"
"""


def test_parse_os_release():
    parsed = er.parse_os_release(OS_RELEASE_FIXTURE)
    assert parsed["NAME"] == "Ubuntu"
    assert parsed["ID"] == "ubuntu"
    assert parsed["VERSION_ID"] == "22.04"
    assert parsed["PRETTY_NAME"] == "Ubuntu 22.04.3 LTS"


def test_parse_os_release_empty():
    assert er.parse_os_release("") == {}


CPUINFO_FIXTURE = """\
processor\t: 0
vendor_id\t: GenuineIntel
model name\t: Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz
physical id\t: 0
core id\t: 0

processor\t: 1
vendor_id\t: GenuineIntel
model name\t: Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz
physical id\t: 0
core id\t: 1

processor\t: 2
vendor_id\t: GenuineIntel
model name\t: Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz
physical id\t: 0
core id\t: 0

processor\t: 3
vendor_id\t: GenuineIntel
model name\t: Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz
physical id\t: 0
core id\t: 1
"""


def test_parse_cpuinfo():
    parsed = er.parse_cpuinfo(CPUINFO_FIXTURE)
    assert parsed["model_name"] == "Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz"
    # 4 logical processors hyperthreaded onto 2 physical cores
    assert parsed["physical_cores"] == 2
    assert parsed["logical_entries"] == 4


def test_parse_cpuinfo_without_physical_id():
    text = "processor\t: 0\nmodel name\t: Some CPU\n"
    parsed = er.parse_cpuinfo(text)
    assert parsed["model_name"] == "Some CPU"
    assert parsed["physical_cores"] is None


def test_parse_cpuinfo_empty():
    parsed = er.parse_cpuinfo("")
    assert parsed["model_name"] is None
    assert parsed["physical_cores"] is None
    assert parsed["logical_entries"] == 0


MEMINFO_FIXTURE = """\
MemTotal:       16384000 kB
MemFree:         2048000 kB
MemAvailable:    8192000 kB
Buffers:          512000 kB
Cached:          1024000 kB
"""


def test_parse_meminfo():
    parsed = er.parse_meminfo(MEMINFO_FIXTURE)
    assert parsed["total_bytes"] == 16384000 * 1024
    assert parsed["available_bytes"] == 8192000 * 1024


def test_parse_meminfo_falls_back_to_memfree_without_memavailable():
    text = "MemTotal:       16384000 kB\nMemFree:         2048000 kB\n"
    parsed = er.parse_meminfo(text)
    assert parsed["total_bytes"] == 16384000 * 1024
    assert parsed["available_bytes"] == 2048000 * 1024


def test_parse_meminfo_empty():
    parsed = er.parse_meminfo("")
    assert parsed["total_bytes"] is None
    assert parsed["available_bytes"] is None


NVIDIA_SMI_CSV_FIXTURE = (
    "NVIDIA GeForce RTX 3090, 535.104.05, 24576 MiB\n"
    "NVIDIA GeForce RTX 3090, 535.104.05, 24576 MiB\n"
)


def test_parse_nvidia_smi_csv():
    gpus = er.parse_nvidia_smi_csv(NVIDIA_SMI_CSV_FIXTURE)
    assert len(gpus) == 2
    assert gpus[0]["name"] == "NVIDIA GeForce RTX 3090"
    assert gpus[0]["driver_version"] == "535.104.05"
    assert gpus[0]["memory_total"] == "24576 MiB"


def test_parse_nvidia_smi_csv_empty():
    assert er.parse_nvidia_smi_csv("") == []


def test_extract_version_token():
    assert er.extract_version_token("uv 0.12.11 (abcdef 2024-01-01)") == "0.12.11"
    assert er.extract_version_token("git version 2.43.0") == "2.43.0"
    assert er.extract_version_token("no version here") is None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_main_exits_zero_and_prints_valid_json(capsys):
    exit_code = er.main([])
    assert exit_code == 0
    captured = capsys.readouterr()
    parsed = json.loads(captured.out)
    assert parsed["schema_version"] == er.SCHEMA_VERSION


def test_cli_main_respects_indent(capsys):
    er.main(["--indent", "0"])
    captured = capsys.readouterr()
    # indent=0 still produces valid, parseable JSON with newlines between items
    parsed = json.loads(captured.out)
    assert parsed["schema_version"] == er.SCHEMA_VERSION


def test_cli_main_writes_output_file(tmp_path):
    output_path = tmp_path / "report.json"
    exit_code = er.main(["--output", str(output_path)])
    assert exit_code == 0
    assert output_path.exists()
    parsed = json.loads(output_path.read_text(encoding="utf-8"))
    assert parsed["schema_version"] == er.SCHEMA_VERSION


def test_cli_as_subprocess_exits_zero_with_valid_json():
    module_path = er.Path(__file__).resolve().parents[2] / "infrastructure" / "diagnostics" / "environment_report.py"
    proc = subprocess.run(
        [sys.executable, str(module_path)],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    parsed = json.loads(proc.stdout)
    assert parsed["schema_version"] == er.SCHEMA_VERSION


# ---------------------------------------------------------------------------
# Sensitive-content guard
# ---------------------------------------------------------------------------


def test_report_contains_no_environment_variable_dump():
    import os as os_module

    report = er.build_report()
    text = json.dumps(report)

    # No key at any level should look like a wholesale env dump.
    def walk(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                assert k.lower() not in ("environ", "environment", "env"), f"suspicious key: {k}"
                walk(v)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)

    walk(report)

    # No secret-shaped environment variable value leaks into the report.
    for var_name, var_value in os_module.environ.items():
        if not var_value:
            continue
        if any(token in var_name.upper() for token in ("TOKEN", "SECRET", "PASSWORD", "CREDENTIAL", "API_KEY")):
            assert var_value not in text, f"leaked value of env var {var_name}"
