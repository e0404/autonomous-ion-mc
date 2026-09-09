"""Tests for infrastructure/diagnostics/repo_capability_report.py.

Capability builders are tested against synthetic fake repository roots
built under `tmp_path`, not against real repo state, so these tests are
deterministic regardless of what infrastructure this checkout happens to
have. The only tests that touch the real checked-out repository assert
true invariants of this checkout (governing docs, this module's own
files) - see the "real repo" section near the end.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime

import pytest

import repo_capability_report as rcr


VALID_STATUSES = {rcr.STATUS_PRESENT, rcr.STATUS_PARTIAL, rcr.STATUS_ABSENT, rcr.STATUS_ERROR}


# ---------------------------------------------------------------------------
# Fake repo builder
# ---------------------------------------------------------------------------


def _make_fake_repo(tmp_path, files: dict[str, str] | None = None, mcp_servers: list[str] | None = None):
    """Build a synthetic (non-Git) fake repository root under tmp_path.

    `files` maps repo-relative path -> content (empty content still
    counts as an empty file for "non-empty" checks). `mcp_servers`, if
    given, writes a well-formed `.mcp.json` declaring those server names.
    """
    root = tmp_path / "fake_repo"
    root.mkdir()
    for relpath, content in (files or {}).items():
        path = root / relpath
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    if mcp_servers is not None:
        mcp_config = {"mcpServers": {name: {"type": "stdio", "command": "uv", "args": []} for name in mcp_servers}}
        (root / ".mcp.json").write_text(json.dumps(mcp_config), encoding="utf-8")

    return root


ALL_CAPABILITY_FILES = {
    "github_ci_workflow": [".github/workflows/ci.yml"],
    "pre_commit_config": [".pre-commit-config.yaml"],
    "local_exact_sha_validation": [
        "infrastructure/validation/local_validation.py",
        "infrastructure/mcp/validation_mcp.py",
    ],
    "task_worktree_manager": ["infrastructure/tasks/task_worktree.py", "infrastructure/mcp/task_mcp.py"],
    "autonomous_commit_helper": ["infrastructure/tasks/task_commit.py", "infrastructure/mcp/commit_mcp.py"],
    "host_gpu_validation_runner": [
        "infrastructure/host_runner/host_runner.py",
        "infrastructure/mcp/host_runner_mcp.py",
    ],
    "codex_worker": ["infrastructure/codex/codex_worker.py", "infrastructure/mcp/codex_mcp.py"],
    "github_ci_observability": [
        "infrastructure/tasks/task_integration.py",
        "infrastructure/mcp/integration_mcp.py",
    ],
}

MCP_SERVER_NAMES = {
    "local_exact_sha_validation": "validation-manager",
    "task_worktree_manager": "task-manager",
    "autonomous_commit_helper": "task-commit",
    "host_gpu_validation_runner": "host-runner",
    "codex_worker": "codex-worker",
    "github_ci_observability": "task-integration",
}


def _fully_present_files() -> dict[str, str]:
    files = {}
    for paths in ALL_CAPABILITY_FILES.values():
        for p in paths:
            files[p] = "content\n"
    for gov_doc in rcr._GOVERNING_DOCS:
        files[gov_doc] = "content\n"
    return files


# ---------------------------------------------------------------------------
# build_report shape / vocabulary invariants
# ---------------------------------------------------------------------------


EXPECTED_TOP_LEVEL_KEYS = {
    "schema_version",
    "generated_at_utc",
    "repository",
    "github_ci_workflow",
    "pre_commit_config",
    "local_exact_sha_validation",
    "task_worktree_manager",
    "autonomous_commit_helper",
    "host_gpu_validation_runner",
    "codex_worker",
    "github_ci_observability",
    "claude_native_subagents",
    "governing_docs",
    "summary",
}


def test_build_report_has_expected_top_level_keys():
    report = rcr.build_report()
    assert EXPECTED_TOP_LEVEL_KEYS.issubset(report.keys())


def test_build_report_is_json_serializable():
    report = rcr.build_report()
    text = json.dumps(report)
    round_tripped = json.loads(text)
    assert round_tripped["schema_version"] == rcr.SCHEMA_VERSION


def test_schema_version_is_int():
    report = rcr.build_report()
    assert isinstance(report["schema_version"], int)
    assert report["schema_version"] == 1


def test_generated_at_utc_is_parseable_and_timezone_aware():
    report = rcr.build_report()
    parsed = datetime.fromisoformat(report["generated_at_utc"])
    assert parsed.tzinfo is not None
    assert parsed.utcoffset().total_seconds() == 0


def test_every_section_has_a_valid_status():
    report = rcr.build_report()
    for key in EXPECTED_TOP_LEVEL_KEYS - {"schema_version", "generated_at_utc", "summary"}:
        section = report[key]
        assert isinstance(section, dict), f"section {key} is not a dict"
        assert "status" in section, f"section {key} has no status"
        assert section["status"] in VALID_STATUSES, f"section {key} has bad status {section['status']!r}"
        if section["status"] != rcr.STATUS_PRESENT:
            assert "detail" in section, f"non-present section {key} missing detail"


def test_summary_shape():
    report = rcr.build_report()
    summary = report["summary"]
    assert isinstance(summary["all_present"], bool)
    assert isinstance(summary["capabilities_not_present"], list)


def test_safe_call_catches_builder_exception():
    def broken_builder():
        raise RuntimeError("builder blew up")

    result = rcr._safe_call(broken_builder)
    assert result["status"] == rcr.STATUS_ERROR
    assert "builder blew up" in result["detail"]


def test_safe_call_catches_malformed_result():
    def malformed_builder():
        return {"no_status_key": True}

    result = rcr._safe_call(malformed_builder)
    assert result["status"] == rcr.STATUS_ERROR


# ---------------------------------------------------------------------------
# repository-root resolution / repository section
# ---------------------------------------------------------------------------


def test_resolve_repo_root_from_module_location_not_cwd(tmp_path, monkeypatch):
    real_repo_root = rcr.Path(__file__).resolve().parents[2]
    monkeypatch.chdir(tmp_path)

    resolved = rcr.resolve_repo_root(module_path=real_repo_root / "infrastructure" / "diagnostics")
    assert resolved == real_repo_root


def test_resolve_repo_root_none_outside_a_repo(tmp_path):
    assert rcr.resolve_repo_root(module_path=tmp_path) is None


def test_repository_section_resolved_from_module_file_not_cwd(tmp_path, monkeypatch):
    repo_dir = rcr.Path(__file__).resolve().parents[2]
    monkeypatch.chdir(tmp_path)

    section = rcr.build_repository_section(module_path=repo_dir / "infrastructure" / "diagnostics")

    assert section["status"] == rcr.STATUS_PRESENT
    assert isinstance(section["branch"], str)
    assert isinstance(section["sha"], str) and len(section["sha"]) == 40
    assert isinstance(section["dirty"], bool)


def test_repository_section_degrades_outside_a_repo(tmp_path):
    section = rcr.build_repository_section(module_path=tmp_path)
    assert section["status"] in (rcr.STATUS_ABSENT, rcr.STATUS_ERROR)
    assert "detail" in section


def test_build_report_degrades_gracefully_when_repo_root_unresolvable(tmp_path):
    report = rcr.build_report(module_path=tmp_path)
    assert report["schema_version"] == rcr.SCHEMA_VERSION
    for key in ALL_CAPABILITY_FILES:
        assert report[key]["status"] == rcr.STATUS_ERROR
    assert report["summary"]["all_present"] is False


# ---------------------------------------------------------------------------
# load_mcp_config
# ---------------------------------------------------------------------------


def test_load_mcp_config_missing_file(tmp_path):
    result = rcr.load_mcp_config(tmp_path)
    assert result["ok"] is False
    assert "detail" in result


def test_load_mcp_config_malformed_json(tmp_path):
    (tmp_path / ".mcp.json").write_text("{not valid json", encoding="utf-8")
    result = rcr.load_mcp_config(tmp_path)
    assert result["ok"] is False
    assert "not valid JSON" in result["detail"]


def test_load_mcp_config_wrong_shape(tmp_path):
    (tmp_path / ".mcp.json").write_text(json.dumps({"nope": {}}), encoding="utf-8")
    result = rcr.load_mcp_config(tmp_path)
    assert result["ok"] is False
    assert "mcpServers" in result["detail"]


def test_load_mcp_config_ok(tmp_path):
    (tmp_path / ".mcp.json").write_text(
        json.dumps({"mcpServers": {"task-manager": {}, "codex-worker": {}}}), encoding="utf-8"
    )
    result = rcr.load_mcp_config(tmp_path)
    assert result["ok"] is True
    assert result["servers"] == {"task-manager", "codex-worker"}


# ---------------------------------------------------------------------------
# Simple file-presence capabilities
# ---------------------------------------------------------------------------


def test_github_ci_workflow_present(tmp_path):
    root = _make_fake_repo(tmp_path, files={".github/workflows/ci.yml": "name: ci\n"})
    section = rcr.build_github_ci_workflow_section(root)
    assert section["status"] == rcr.STATUS_PRESENT


def test_github_ci_workflow_absent_when_missing(tmp_path):
    root = _make_fake_repo(tmp_path)
    section = rcr.build_github_ci_workflow_section(root)
    assert section["status"] == rcr.STATUS_ABSENT
    assert "detail" in section


def test_github_ci_workflow_absent_when_empty(tmp_path):
    root = _make_fake_repo(tmp_path, files={".github/workflows/ci.yml": ""})
    section = rcr.build_github_ci_workflow_section(root)
    assert section["status"] == rcr.STATUS_ABSENT


def test_pre_commit_config_present(tmp_path):
    root = _make_fake_repo(tmp_path, files={".pre-commit-config.yaml": "repos: []\n"})
    section = rcr.build_pre_commit_config_section(root)
    assert section["status"] == rcr.STATUS_PRESENT


def test_pre_commit_config_absent(tmp_path):
    root = _make_fake_repo(tmp_path)
    section = rcr.build_pre_commit_config_section(root)
    assert section["status"] == rcr.STATUS_ABSENT


# ---------------------------------------------------------------------------
# MCP-backed capabilities (parametrized over all six)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("capability_key", list(ALL_CAPABILITY_FILES.keys() - {"github_ci_workflow", "pre_commit_config"}))
def test_mcp_backed_capability_present(tmp_path, capability_key):
    required_files = ALL_CAPABILITY_FILES[capability_key]
    server_name = MCP_SERVER_NAMES[capability_key]
    root = _make_fake_repo(
        tmp_path,
        files={p: "content\n" for p in required_files},
        mcp_servers=[server_name],
    )
    mcp_config = rcr.load_mcp_config(root)
    builder = getattr(rcr, f"build_{capability_key}_section")
    section = builder(root, mcp_config)
    assert section["status"] == rcr.STATUS_PRESENT
    assert section["mcp_registered"] is True
    assert all(section["files"].values())


@pytest.mark.parametrize("capability_key", list(ALL_CAPABILITY_FILES.keys() - {"github_ci_workflow", "pre_commit_config"}))
def test_mcp_backed_capability_absent(tmp_path, capability_key):
    root = _make_fake_repo(tmp_path, mcp_servers=[])
    mcp_config = rcr.load_mcp_config(root)
    builder = getattr(rcr, f"build_{capability_key}_section")
    section = builder(root, mcp_config)
    assert section["status"] == rcr.STATUS_ABSENT
    assert "detail" in section


@pytest.mark.parametrize("capability_key", list(ALL_CAPABILITY_FILES.keys() - {"github_ci_workflow", "pre_commit_config"}))
def test_mcp_backed_capability_partial_missing_one_file(tmp_path, capability_key):
    required_files = ALL_CAPABILITY_FILES[capability_key]
    server_name = MCP_SERVER_NAMES[capability_key]
    # Only write the first required file, not the second.
    root = _make_fake_repo(
        tmp_path,
        files={required_files[0]: "content\n"},
        mcp_servers=[server_name],
    )
    mcp_config = rcr.load_mcp_config(root)
    builder = getattr(rcr, f"build_{capability_key}_section")
    section = builder(root, mcp_config)
    assert section["status"] == rcr.STATUS_PARTIAL
    assert "detail" in section


@pytest.mark.parametrize("capability_key", list(ALL_CAPABILITY_FILES.keys() - {"github_ci_workflow", "pre_commit_config"}))
def test_mcp_backed_capability_partial_missing_mcp_registration(tmp_path, capability_key):
    required_files = ALL_CAPABILITY_FILES[capability_key]
    root = _make_fake_repo(
        tmp_path,
        files={p: "content\n" for p in required_files},
        mcp_servers=["some-other-server"],
    )
    mcp_config = rcr.load_mcp_config(root)
    builder = getattr(rcr, f"build_{capability_key}_section")
    section = builder(root, mcp_config)
    assert section["status"] == rcr.STATUS_PARTIAL
    assert section["mcp_registered"] is False
    assert "detail" in section


@pytest.mark.parametrize("capability_key", list(ALL_CAPABILITY_FILES.keys() - {"github_ci_workflow", "pre_commit_config"}))
def test_mcp_backed_capability_partial_when_mcp_json_malformed(tmp_path, capability_key):
    required_files = ALL_CAPABILITY_FILES[capability_key]
    root = _make_fake_repo(tmp_path, files={p: "content\n" for p in required_files})
    (root / ".mcp.json").write_text("{not valid json", encoding="utf-8")
    mcp_config = rcr.load_mcp_config(root)
    assert mcp_config["ok"] is False

    builder = getattr(rcr, f"build_{capability_key}_section")
    section = builder(root, mcp_config)
    # Files are present but MCP confirmation is impossible -> not "present".
    assert section["status"] != rcr.STATUS_PRESENT
    assert section["mcp_registered"] is None
    assert "detail" in section


@pytest.mark.parametrize("capability_key", list(ALL_CAPABILITY_FILES.keys() - {"github_ci_workflow", "pre_commit_config"}))
def test_mcp_backed_capability_when_mcp_json_missing_entirely(tmp_path, capability_key):
    required_files = ALL_CAPABILITY_FILES[capability_key]
    root = _make_fake_repo(tmp_path, files={p: "content\n" for p in required_files})
    mcp_config = rcr.load_mcp_config(root)
    assert mcp_config["ok"] is False

    builder = getattr(rcr, f"build_{capability_key}_section")
    section = builder(root, mcp_config)
    assert section["status"] != rcr.STATUS_PRESENT
    assert section["mcp_registered"] is None


# ---------------------------------------------------------------------------
# claude_native_subagents
# ---------------------------------------------------------------------------


GOOD_FRONTMATTER = """---
name: implementation-worker
description: does things
---

Body text.
"""

MALFORMED_NO_CLOSING_FENCE = """---
name: broken-agent
description: no closing fence
"""

MALFORMED_NO_NAME_KEY = """---
description: no name key here
---

Body.
"""

NOT_FRONTMATTER_AT_ALL = "just a regular markdown file\nwith no frontmatter\n"


def test_parse_agent_frontmatter_name_good():
    assert rcr.parse_agent_frontmatter_name(GOOD_FRONTMATTER) == "implementation-worker"


def test_parse_agent_frontmatter_name_missing_closing_fence():
    assert rcr.parse_agent_frontmatter_name(MALFORMED_NO_CLOSING_FENCE) is None


def test_parse_agent_frontmatter_name_missing_name_key():
    assert rcr.parse_agent_frontmatter_name(MALFORMED_NO_NAME_KEY) is None


def test_parse_agent_frontmatter_name_no_frontmatter():
    assert rcr.parse_agent_frontmatter_name(NOT_FRONTMATTER_AT_ALL) is None


def test_parse_agent_frontmatter_name_empty_string():
    assert rcr.parse_agent_frontmatter_name("") is None


def test_parse_agent_frontmatter_name_quoted_value():
    text = '---\nname: "quoted-name"\n---\n'
    assert rcr.parse_agent_frontmatter_name(text) == "quoted-name"


def test_claude_native_subagents_present(tmp_path):
    root = _make_fake_repo(tmp_path, files={".claude/agents/foo.md": GOOD_FRONTMATTER})
    section = rcr.build_claude_native_subagents_section(root)
    assert section["status"] == rcr.STATUS_PRESENT
    assert section["agent_names"] == ["implementation-worker"]


def test_claude_native_subagents_absent_no_dir(tmp_path):
    root = _make_fake_repo(tmp_path)
    section = rcr.build_claude_native_subagents_section(root)
    assert section["status"] == rcr.STATUS_ABSENT


def test_claude_native_subagents_absent_only_malformed_files(tmp_path):
    root = _make_fake_repo(
        tmp_path,
        files={
            ".claude/agents/broken.md": MALFORMED_NO_CLOSING_FENCE,
            ".claude/agents/noname.md": MALFORMED_NO_NAME_KEY,
        },
    )
    section = rcr.build_claude_native_subagents_section(root)
    assert section["status"] == rcr.STATUS_ABSENT
    assert section["agent_names"] == []
    assert len(section["skipped_files"]) == 2


def test_claude_native_subagents_mixed_good_and_malformed(tmp_path):
    root = _make_fake_repo(
        tmp_path,
        files={
            ".claude/agents/good.md": GOOD_FRONTMATTER,
            ".claude/agents/broken.md": MALFORMED_NO_CLOSING_FENCE,
        },
    )
    section = rcr.build_claude_native_subagents_section(root)
    assert section["status"] == rcr.STATUS_PRESENT
    assert section["agent_names"] == ["implementation-worker"]
    assert len(section["skipped_files"]) == 1


# ---------------------------------------------------------------------------
# governing_docs
# ---------------------------------------------------------------------------


def test_governing_docs_all_present(tmp_path):
    root = _make_fake_repo(tmp_path, files={name: "content\n" for name in rcr._GOVERNING_DOCS})
    section = rcr.build_governing_docs_section(root)
    assert section["status"] == rcr.STATUS_PRESENT
    assert section["all_present"] is True


def test_governing_docs_all_absent(tmp_path):
    root = _make_fake_repo(tmp_path)
    section = rcr.build_governing_docs_section(root)
    assert section["status"] == rcr.STATUS_ABSENT
    assert section["all_present"] is False


def test_governing_docs_partial(tmp_path):
    root = _make_fake_repo(tmp_path, files={"EXPERIMENT.md": "content\n"})
    section = rcr.build_governing_docs_section(root)
    assert section["status"] == rcr.STATUS_PARTIAL
    assert section["all_present"] is False
    assert section["docs"]["EXPERIMENT.md"] is True
    assert section["docs"]["REQUIREMENTS.md"] is False


# ---------------------------------------------------------------------------
# Full synthetic fake-repo report (via build_report(module_path=...))
# ---------------------------------------------------------------------------


def test_build_report_all_present_on_fully_populated_fake_repo(tmp_path, monkeypatch):
    files = _fully_present_files()
    all_servers = list(MCP_SERVER_NAMES.values())
    root = _make_fake_repo(tmp_path, files=files, mcp_servers=all_servers)
    (root / ".claude" / "agents").mkdir(parents=True, exist_ok=True)
    (root / ".claude" / "agents" / "foo.md").write_text(GOOD_FRONTMATTER, encoding="utf-8")

    # resolve_repo_root/build_repository_section use `git`; the fake repo
    # is not a real Git repository, so stub repository resolution to point
    # capability checks at the fake root while leaving repository-section
    # git plumbing untouched (tested separately above).
    monkeypatch.setattr(rcr, "resolve_repo_root", lambda module_path=None, timeout=rcr.SUBPROCESS_TIMEOUT_SECONDS: root)
    monkeypatch.setattr(
        rcr,
        "build_repository_section",
        lambda module_path=None, timeout=rcr.SUBPROCESS_TIMEOUT_SECONDS: {
            "status": rcr.STATUS_PRESENT,
            "branch": "fake",
            "sha": "0" * 40,
            "dirty": False,
        },
    )

    report = rcr.build_report(module_path=root)
    assert report["summary"]["all_present"] is True
    assert report["summary"]["capabilities_not_present"] == []
    for key in ALL_CAPABILITY_FILES:
        assert report[key]["status"] == rcr.STATUS_PRESENT
    assert report["claude_native_subagents"]["status"] == rcr.STATUS_PRESENT
    assert report["governing_docs"]["status"] == rcr.STATUS_PRESENT


def test_build_report_summary_lists_missing_on_empty_fake_repo(tmp_path, monkeypatch):
    root = _make_fake_repo(tmp_path)
    monkeypatch.setattr(rcr, "resolve_repo_root", lambda module_path=None, timeout=rcr.SUBPROCESS_TIMEOUT_SECONDS: root)
    monkeypatch.setattr(
        rcr,
        "build_repository_section",
        lambda module_path=None, timeout=rcr.SUBPROCESS_TIMEOUT_SECONDS: {
            "status": rcr.STATUS_PRESENT,
            "branch": "fake",
            "sha": "0" * 40,
            "dirty": False,
        },
    )

    report = rcr.build_report(module_path=root)
    assert report["summary"]["all_present"] is False
    assert set(report["summary"]["capabilities_not_present"]) == set(ALL_CAPABILITY_FILES.keys()) | {
        "claude_native_subagents",
        "governing_docs",
    }


# ---------------------------------------------------------------------------
# Real-repo light integration check (true invariants only)
# ---------------------------------------------------------------------------


def test_build_report_against_real_repo_governing_docs_and_self_present():
    report = rcr.build_report()
    assert report["governing_docs"]["status"] == rcr.STATUS_PRESENT
    assert report["repository"]["status"] == rcr.STATUS_PRESENT
    # This module and its docs exist in this exact checkout - a true
    # invariant, not incidental host state.
    assert report["local_exact_sha_validation"]["files"]["infrastructure/mcp/validation_mcp.py"] is True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_cli_main_exits_zero_and_prints_valid_json(capsys):
    exit_code = rcr.main([])
    assert exit_code == 0
    captured = capsys.readouterr()
    parsed = json.loads(captured.out)
    assert parsed["schema_version"] == rcr.SCHEMA_VERSION


def test_cli_main_respects_indent(capsys):
    rcr.main(["--indent", "0"])
    captured = capsys.readouterr()
    parsed = json.loads(captured.out)
    assert parsed["schema_version"] == rcr.SCHEMA_VERSION


def test_cli_main_writes_output_file(tmp_path):
    output_path = tmp_path / "report.json"
    exit_code = rcr.main(["--output", str(output_path)])
    assert exit_code == 0
    assert output_path.exists()
    parsed = json.loads(output_path.read_text(encoding="utf-8"))
    assert parsed["schema_version"] == rcr.SCHEMA_VERSION


def test_cli_main_exits_one_with_structured_error_when_output_unwritable(tmp_path, capsys):
    unwritable_path = tmp_path / "no-such-directory" / "report.json"
    exit_code = rcr.main(["--output", str(unwritable_path)])
    assert exit_code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    error_payload = json.loads(captured.err)
    assert error_payload["status"] == rcr.STATUS_ERROR
    assert "detail" in error_payload
    assert not unwritable_path.exists()


def test_cli_main_bad_argument_exits_two():
    with pytest.raises(SystemExit) as exc_info:
        rcr.main(["--not-a-real-flag"])
    assert exc_info.value.code == 2


def test_cli_subprocess_exits_one_no_traceback_when_output_unwritable(tmp_path):
    unwritable_path = tmp_path / "no-such-directory" / "report.json"
    module_path = rcr.Path(__file__).resolve().parents[2] / "infrastructure" / "diagnostics" / "repo_capability_report.py"
    proc = subprocess.run(
        [sys.executable, str(module_path), "--output", str(unwritable_path)],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
        text=True,
    )
    assert proc.returncode == 1
    assert "Traceback" not in proc.stderr
    error_payload = json.loads(proc.stderr)
    assert error_payload["status"] == rcr.STATUS_ERROR


def test_cli_as_subprocess_exits_zero_with_valid_json():
    module_path = rcr.Path(__file__).resolve().parents[2] / "infrastructure" / "diagnostics" / "repo_capability_report.py"
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
    assert parsed["schema_version"] == rcr.SCHEMA_VERSION


# ---------------------------------------------------------------------------
# Sensitive-content guard
# ---------------------------------------------------------------------------


def test_report_contains_no_environment_variable_dump():
    import os as os_module

    report = rcr.build_report()
    text = json.dumps(report)

    def walk(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                assert k.lower() not in ("environ", "environment", "env"), f"suspicious key: {k}"
                walk(v)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)

    walk(report)

    for var_name, var_value in os_module.environ.items():
        if not var_value:
            continue
        if any(token in var_name.upper() for token in ("TOKEN", "SECRET", "PASSWORD", "CREDENTIAL", "API_KEY")):
            assert var_value not in text, f"leaked value of env var {var_name}"
