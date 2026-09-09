#!/usr/bin/env python3
"""Repository-capability diagnostics utility.

This is the final infrastructure-acceptance "dry run" (DRYRUN-003). It
verifies, as a machine-readable report, that the autonomous-development
lifecycle infrastructure introduced by DRYRUN-001
(`infrastructure/diagnostics/environment_report.py`) and DRYRUN-002
(`infrastructure/diagnostics/warp_backend_report.py`) — and by the
autonomy-lifecycle features built alongside them (task worktrees, commit
helper, host/GPU validation runner, Codex worker, GitHub CI
observability, native subagents) — is actually present and wired up in
this repository checkout. This utility itself concerns only
infrastructure; it contains no Monte Carlo/physics content.

Design goals:

* Standard library only. No third-party imports.
* **Deterministic apart from the repository's current branch/SHA.** Every
  capability check is a pure filesystem/repo-content check (a file
  exists, is non-empty, and — where relevant — a specific marker
  string/JSON key is present in a config file). Nothing shells out to a
  tool whose behavior could vary by host (no `--version` probes, no
  network access). The only environment-derived, host-varying fields
  anywhere in the report are `repository.branch`/`repository.sha`/
  `repository.dirty`, plus `generated_at_utc` (timestamp-varying by
  nature, not host-varying).
* Every probe/check function is independently guarded so it can never
  raise; the top-level report is always produced (mirrors
  `_safe_section` in `environment_report.py`).
* Every path checked is anchored at the resolved repository root, never
  at the process's current working directory.

Status vocabulary used by every capability section's "status" field
(deliberately distinct from `environment_report.py`'s host-probe
vocabulary, since this module reports *capability presence*, not *host
probe outcomes*):

* "present"  - every file/marker required for this capability was found
               (non-empty files, and any required config-declared marker
               present). The section's data fields are populated.
* "partial"  - some, but not all, of the required files/markers for this
               capability were found. A "detail" field explains what is
               missing.
* "absent"   - none of the required files/markers for this capability
               were found.
* "error"    - checking the capability failed unexpectedly (e.g.
               `.mcp.json` exists but fails to parse as JSON, or a
               filesystem read raised unexpectedly for a reason other
               than "file does not exist"). A "detail" field explains
               why.

Every section carries a "status" field from this vocabulary; every
section whose status is not "present" also carries a "detail" string.

No secrets: this utility never reads environment variables, credentials,
or remote URLs, and never dumps `os.environ`. It only reads a small,
fixed set of repository-relative file paths and `.mcp.json`'s already
version-controlled, non-secret server-registration structure.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

SCHEMA_VERSION = 1

SUBPROCESS_TIMEOUT_SECONDS = 5.0

STATUS_PRESENT = "present"
STATUS_PARTIAL = "partial"
STATUS_ABSENT = "absent"
STATUS_ERROR = "error"

VALID_STATUSES = {STATUS_PRESENT, STATUS_PARTIAL, STATUS_ABSENT, STATUS_ERROR}


# ---------------------------------------------------------------------------
# Subprocess helper (duplicated from environment_report.py deliberately -
# these diagnostic modules are plain scripts, not a shared package)
# ---------------------------------------------------------------------------


def run_subprocess(
    args: list[str],
    cwd: Path | str | None = None,
    timeout: float = SUBPROCESS_TIMEOUT_SECONDS,
) -> dict[str, Any]:
    """Run a subprocess defensively.

    Never inherits stdin, always has a timeout, and translates every
    failure mode (missing executable, permission error, timeout,
    non-zero exit, any other unexpected exception) into a structured
    result rather than letting an exception propagate.

    Returns a dict with ``ok`` (bool). On success it also carries
    ``stdout``; on failure it carries ``status`` (one of STATUS_ABSENT /
    STATUS_ERROR) and ``detail``.
    """
    try:
        proc = subprocess.run(
            list(args),
            cwd=cwd,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            text=True,
        )
    except FileNotFoundError:
        return {
            "ok": False,
            "status": STATUS_ABSENT,
            "detail": f"executable not found: {args[0]}",
        }
    except PermissionError as exc:
        return {
            "ok": False,
            "status": STATUS_ERROR,
            "detail": f"permission denied running {args[0]}: {exc}",
        }
    except subprocess.TimeoutExpired:
        return {
            "ok": False,
            "status": STATUS_ERROR,
            "detail": f"{args[0]} timed out after {timeout}s",
        }
    except Exception as exc:  # noqa: BLE001 - probes must never raise
        return {
            "ok": False,
            "status": STATUS_ERROR,
            "detail": f"unexpected error running {args[0]}: {exc}",
        }

    if proc.returncode != 0:
        stderr = (proc.stderr or "").strip()
        return {
            "ok": False,
            "status": STATUS_ERROR,
            "detail": f"{args[0]} exited {proc.returncode}: {stderr[:500]}",
        }

    return {"ok": True, "stdout": (proc.stdout or "").strip()}


# ---------------------------------------------------------------------------
# Repository-root resolution
# ---------------------------------------------------------------------------


def resolve_repo_root(
    module_path: Path | str | None = None,
    timeout: float = SUBPROCESS_TIMEOUT_SECONDS,
) -> Path | None:
    """Resolve the repository root containing this module.

    Resolved from the file location of this module (``module_path``
    defaults to ``__file__``), never from the process's current working
    directory, via ``git rev-parse --show-toplevel`` run with ``cwd`` set
    to that location. Returns ``None`` if the location is not inside a
    Git working tree or the probe otherwise fails.
    """
    module_dir = Path(module_path if module_path is not None else __file__).resolve()
    if module_dir.is_file():
        module_dir = module_dir.parent

    result = run_subprocess(["git", "rev-parse", "--show-toplevel"], cwd=module_dir, timeout=timeout)
    if not result["ok"]:
        return None
    return Path(result["stdout"]).resolve()


def build_repository_section(
    module_path: Path | str | None = None,
    timeout: float = SUBPROCESS_TIMEOUT_SECONDS,
) -> dict[str, Any]:
    """Report the branch/SHA/dirty state of the repo containing this module.

    Mirrors `environment_report.build_repository_section`: resolved from
    the module's own file location, never the process CWD.
    """
    module_dir = Path(module_path if module_path is not None else __file__).resolve()
    if module_dir.is_file():
        module_dir = module_dir.parent

    branch_result = run_subprocess(["git", "branch", "--show-current"], cwd=module_dir, timeout=timeout)
    if not branch_result["ok"]:
        return {"status": branch_result["status"], "detail": branch_result["detail"]}

    sha_result = run_subprocess(["git", "rev-parse", "HEAD"], cwd=module_dir, timeout=timeout)
    if not sha_result["ok"]:
        return {"status": sha_result["status"], "detail": sha_result["detail"]}

    status_result = run_subprocess(["git", "status", "--porcelain"], cwd=module_dir, timeout=timeout)
    if not status_result["ok"]:
        return {"status": status_result["status"], "detail": status_result["detail"]}

    return {
        "status": STATUS_PRESENT,
        "branch": branch_result["stdout"],
        "sha": sha_result["stdout"],
        "dirty": bool(status_result["stdout"]),
    }


# ---------------------------------------------------------------------------
# Small filesystem helpers
# ---------------------------------------------------------------------------


def _is_nonempty_file(path: Path) -> bool:
    """True if `path` exists, is a regular file, and has non-zero size."""
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def load_mcp_config(repo_root: Path) -> dict[str, Any]:
    """Read and parse `.mcp.json` at the repository root.

    Returns a dict with ``ok`` (bool). On success it also carries
    ``servers`` (the set of declared `mcpServers` names). On failure it
    carries ``detail`` explaining why (missing file, unreadable, not
    valid JSON, or malformed shape).
    """
    mcp_path = repo_root / ".mcp.json"
    if not mcp_path.is_file():
        return {"ok": False, "detail": f"{mcp_path} does not exist"}

    try:
        text = mcp_path.read_text(encoding="utf-8")
    except OSError as exc:
        return {"ok": False, "detail": f"failed to read {mcp_path}: {exc}"}

    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        return {"ok": False, "detail": f"{mcp_path} is not valid JSON: {exc}"}

    servers = data.get("mcpServers") if isinstance(data, dict) else None
    if not isinstance(servers, dict):
        return {"ok": False, "detail": f"{mcp_path} has no 'mcpServers' object"}

    return {"ok": True, "servers": set(servers.keys())}


# ---------------------------------------------------------------------------
# Generic "N files + optional MCP marker" capability check
# ---------------------------------------------------------------------------


def _check_files_and_mcp_marker(
    repo_root: Path,
    required_relpaths: list[str],
    mcp_server_name: str | None,
    mcp_config: dict[str, Any] | None,
) -> dict[str, Any]:
    """Shared logic for capabilities defined as "these files exist and are
    non-empty, AND .mcp.json declares this server name".

    `mcp_config` is the pre-loaded result of `load_mcp_config()` (shared
    across capabilities so a single missing/malformed `.mcp.json` is read
    only once), or `None` if this capability has no MCP-registration
    requirement.
    """
    files: dict[str, bool] = {}
    for relpath in required_relpaths:
        files[relpath] = _is_nonempty_file(repo_root / relpath)

    missing_files = [relpath for relpath, present in files.items() if not present]

    mcp_detail: str | None = None
    mcp_registered: bool | None = None
    if mcp_server_name is not None:
        if mcp_config is None or not mcp_config.get("ok"):
            mcp_registered = None
            mcp_detail = (
                mcp_config.get("detail", ".mcp.json could not be read") if mcp_config else ".mcp.json could not be read"
            )
        else:
            mcp_registered = mcp_server_name in mcp_config["servers"]
            if not mcp_registered:
                mcp_detail = f"'{mcp_server_name}' not declared in .mcp.json mcpServers"

    detail_parts: list[str] = []
    if missing_files:
        detail_parts.append(f"missing/empty files: {sorted(missing_files)}")
    if mcp_detail is not None:
        detail_parts.append(mcp_detail)

    result: dict[str, Any] = {
        "files": files,
    }
    if mcp_server_name is not None:
        result["mcp_server_name"] = mcp_server_name
        result["mcp_registered"] = mcp_registered

    all_files_present = not missing_files
    mcp_ok = mcp_server_name is None or mcp_registered is True

    if all_files_present and mcp_ok:
        result["status"] = STATUS_PRESENT
    elif not files or all(not present for present in files.values()):
        if mcp_server_name is not None and mcp_registered:
            # Files entirely missing but MCP registration present: still partial.
            result["status"] = STATUS_PARTIAL
            result["detail"] = "; ".join(detail_parts) or "capability not fully present"
        else:
            result["status"] = STATUS_ABSENT
            result["detail"] = "; ".join(detail_parts) or "no required files present"
    else:
        result["status"] = STATUS_PARTIAL
        result["detail"] = "; ".join(detail_parts) or "capability not fully present"

    return result


# ---------------------------------------------------------------------------
# Capability section builders
# ---------------------------------------------------------------------------


def build_github_ci_workflow_section(repo_root: Path) -> dict[str, Any]:
    path = repo_root / ".github" / "workflows" / "ci.yml"
    present = _is_nonempty_file(path)
    if present:
        return {"status": STATUS_PRESENT, "files": {".github/workflows/ci.yml": True}}
    return {
        "status": STATUS_ABSENT,
        "detail": f"{path} does not exist or is empty",
        "files": {".github/workflows/ci.yml": False},
    }


def build_pre_commit_config_section(repo_root: Path) -> dict[str, Any]:
    path = repo_root / ".pre-commit-config.yaml"
    present = _is_nonempty_file(path)
    if present:
        return {"status": STATUS_PRESENT, "files": {".pre-commit-config.yaml": True}}
    return {
        "status": STATUS_ABSENT,
        "detail": f"{path} does not exist or is empty",
        "files": {".pre-commit-config.yaml": False},
    }


def build_local_exact_sha_validation_section(repo_root: Path, mcp_config: dict[str, Any] | None) -> dict[str, Any]:
    return _check_files_and_mcp_marker(
        repo_root,
        ["infrastructure/validation/local_validation.py", "infrastructure/mcp/validation_mcp.py"],
        "validation-manager",
        mcp_config,
    )


def build_task_worktree_manager_section(repo_root: Path, mcp_config: dict[str, Any] | None) -> dict[str, Any]:
    return _check_files_and_mcp_marker(
        repo_root,
        ["infrastructure/tasks/task_worktree.py", "infrastructure/mcp/task_mcp.py"],
        "task-manager",
        mcp_config,
    )


def build_autonomous_commit_helper_section(repo_root: Path, mcp_config: dict[str, Any] | None) -> dict[str, Any]:
    return _check_files_and_mcp_marker(
        repo_root,
        ["infrastructure/tasks/task_commit.py", "infrastructure/mcp/commit_mcp.py"],
        "task-commit",
        mcp_config,
    )


def build_host_gpu_validation_runner_section(repo_root: Path, mcp_config: dict[str, Any] | None) -> dict[str, Any]:
    return _check_files_and_mcp_marker(
        repo_root,
        ["infrastructure/host_runner/host_runner.py", "infrastructure/mcp/host_runner_mcp.py"],
        "host-runner",
        mcp_config,
    )


def build_codex_worker_section(repo_root: Path, mcp_config: dict[str, Any] | None) -> dict[str, Any]:
    return _check_files_and_mcp_marker(
        repo_root,
        ["infrastructure/codex/codex_worker.py", "infrastructure/mcp/codex_mcp.py"],
        "codex-worker",
        mcp_config,
    )


def build_github_ci_observability_section(repo_root: Path, mcp_config: dict[str, Any] | None) -> dict[str, Any]:
    return _check_files_and_mcp_marker(
        repo_root,
        ["infrastructure/tasks/task_integration.py", "infrastructure/mcp/integration_mcp.py"],
        "task-integration",
        mcp_config,
    )


def parse_agent_frontmatter_name(text: str) -> str | None:
    """Extract the `name:` value from a `.claude/agents/*.md` file's YAML
    frontmatter.

    Deliberately minimal/defensive: does not depend on a YAML library.
    Frontmatter is the block between the first two lines that are exactly
    `---`; within it, each `key: value` line is split on the first `:`.
    Returns `None` (never raises) if the file has no well-formed
    frontmatter block or no `name` key within it.
    """
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return None

    end_index = None
    for index in range(1, len(lines)):
        if lines[index].strip() == "---":
            end_index = index
            break
    if end_index is None:
        return None

    for line in lines[1:end_index]:
        if ":" not in line:
            continue
        key, _, value = line.partition(":")
        if key.strip() == "name":
            return value.strip().strip("\"'") or None

    return None


def build_claude_native_subagents_section(repo_root: Path) -> dict[str, Any]:
    agents_dir = repo_root / ".claude" / "agents"
    if not agents_dir.is_dir():
        return {
            "status": STATUS_ABSENT,
            "detail": f"{agents_dir} does not exist",
            "agent_names": [],
            "skipped_files": [],
        }

    try:
        candidate_paths = sorted(agents_dir.glob("*.md"))
    except OSError as exc:
        return {"status": STATUS_ERROR, "detail": f"failed to list {agents_dir}: {exc}"}

    agent_names: list[str] = []
    skipped_files: list[str] = []
    for path in candidate_paths:
        try:
            text = path.read_text(encoding="utf-8")
        except OSError as exc:
            skipped_files.append(f"{path.name} (read error: {exc})")
            continue
        try:
            name = parse_agent_frontmatter_name(text)
        except Exception as exc:  # noqa: BLE001 - parser must never raise
            skipped_files.append(f"{path.name} (parse error: {exc})")
            continue
        if name:
            agent_names.append(name)
        else:
            skipped_files.append(f"{path.name} (no frontmatter 'name' key)")

    if agent_names:
        status = STATUS_PRESENT
        detail = None
    elif candidate_paths:
        status = STATUS_ABSENT
        detail = "found *.md files under .claude/agents/ but none had a parseable 'name' frontmatter key"
    else:
        status = STATUS_ABSENT
        detail = f"no *.md files found under {agents_dir}"

    result: dict[str, Any] = {
        "status": status,
        "agent_names": sorted(agent_names),
        "skipped_files": skipped_files,
    }
    if detail is not None:
        result["detail"] = detail
    return result


_GOVERNING_DOCS = ["EXPERIMENT.md", "REQUIREMENTS.md", "AGENTS.md", "CLAUDE.md"]


def build_governing_docs_section(repo_root: Path) -> dict[str, Any]:
    docs: dict[str, bool] = {name: _is_nonempty_file(repo_root / name) for name in _GOVERNING_DOCS}
    missing = [name for name, present in docs.items() if not present]

    result: dict[str, Any] = {"docs": docs, "all_present": not missing}
    if not missing:
        result["status"] = STATUS_PRESENT
    elif len(missing) == len(_GOVERNING_DOCS):
        result["status"] = STATUS_ABSENT
        result["detail"] = f"missing/empty: {missing}"
    else:
        result["status"] = STATUS_PARTIAL
        result["detail"] = f"missing/empty: {missing}"
    return result


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

# Capability keys that require .mcp.json + a set of files, in the shared
# (repo_root, mcp_config) -> section shape.
_MCP_BACKED_BUILDERS: dict[str, Callable[[Path, dict[str, Any] | None], dict[str, Any]]] = {
    "local_exact_sha_validation": build_local_exact_sha_validation_section,
    "task_worktree_manager": build_task_worktree_manager_section,
    "autonomous_commit_helper": build_autonomous_commit_helper_section,
    "host_gpu_validation_runner": build_host_gpu_validation_runner_section,
    "codex_worker": build_codex_worker_section,
    "github_ci_observability": build_github_ci_observability_section,
}

# Capability keys that only need the repo root.
_SIMPLE_BUILDERS: dict[str, Callable[[Path], dict[str, Any]]] = {
    "github_ci_workflow": build_github_ci_workflow_section,
    "pre_commit_config": build_pre_commit_config_section,
    "claude_native_subagents": build_claude_native_subagents_section,
    "governing_docs": build_governing_docs_section,
}

_CAPABILITY_KEYS = list(_SIMPLE_BUILDERS.keys()) + list(_MCP_BACKED_BUILDERS.keys())


def _safe_call(builder: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    """Outer safety net around a section builder call.

    Every builder above already guards its own failure modes, but this
    ensures a report can never fail to be produced even if a builder is
    later changed to raise unexpectedly.
    """
    try:
        result = builder()
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": f"unexpected error: {exc}"}

    if not isinstance(result, dict) or "status" not in result:
        return {"status": STATUS_ERROR, "detail": "capability check returned a malformed result"}

    return result


def build_report(module_path: Path | str | None = None) -> dict[str, Any]:
    """Build the full repository-capability diagnostics report."""
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }

    report["repository"] = _safe_call(lambda: build_repository_section(module_path=module_path))

    repo_root = resolve_repo_root(module_path=module_path)
    if repo_root is None:
        # Every capability check is anchored at the repo root; if it can't
        # be resolved at all, every capability section reports "error"
        # (a genuinely unexpected condition for a module shipped inside
        # this repository) rather than the report failing to be produced.
        error_section = {
            "status": STATUS_ERROR,
            "detail": "could not resolve repository root from this module's location",
        }
        for key in _CAPABILITY_KEYS:
            report[key] = dict(error_section)
    else:
        mcp_config = load_mcp_config(repo_root)

        for key, simple_builder in _SIMPLE_BUILDERS.items():
            report[key] = _safe_call(lambda b=simple_builder: b(repo_root))

        for key, mcp_builder in _MCP_BACKED_BUILDERS.items():
            report[key] = _safe_call(lambda b=mcp_builder: b(repo_root, mcp_config))

    not_present = sorted(
        key
        for key in _CAPABILITY_KEYS
        if report[key].get("status") != STATUS_PRESENT
    )
    report["summary"] = {
        "all_present": not not_present,
        "capabilities_not_present": not_present,
    }

    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="repo_capability_report",
        description="Report repository lifecycle-infrastructure capability diagnostics as JSON.",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation level (default: 2).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Write the report to this path instead of stdout.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point.

    Exit code policy:

    * 0 - the report was produced (even if some capabilities are
      "absent"/"partial"/"error") and successfully delivered to stdout or
      --output.
    * 1 - the report was produced but could not be written to the
      --output path (e.g. an unwritable/nonexistent parent directory).
      A structured JSON error object is printed to stderr; no traceback
      is allowed to escape.
    * 2 - argparse usage error (e.g. an unrecognized argument), raised
      by argparse itself before report generation.
    """
    args = build_arg_parser().parse_args(argv)

    report = build_report()
    text = json.dumps(report, indent=args.indent, sort_keys=True) + "\n"

    if args.output is not None:
        try:
            args.output.write_text(text, encoding="utf-8")
        except OSError as exc:
            error_payload = {
                "status": STATUS_ERROR,
                "detail": f"failed to write report to {args.output}: {exc}",
            }
            sys.stderr.write(json.dumps(error_payload) + "\n")
            return 1
    else:
        sys.stdout.write(text)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
