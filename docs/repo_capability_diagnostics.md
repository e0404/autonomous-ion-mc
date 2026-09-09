# Repository Capability Diagnostics Utility

`infrastructure/diagnostics/repo_capability_report.py` is the final
infrastructure-acceptance "dry run" (DRYRUN-003). It reports, as a single
machine-readable JSON document, whether the autonomous-development
lifecycle infrastructure introduced by DRYRUN-001
(`infrastructure/diagnostics/environment_report.py`) and DRYRUN-002
(`infrastructure/diagnostics/warp_backend_report.py`) — and the
autonomy-lifecycle features built alongside them (task worktrees, the
autonomous commit helper, the host/GPU validation runner, the Codex
worker, GitHub CI observability, and native Claude subagents) — is
actually present and wired up in this repository checkout. It is an
infrastructure-audit tool only: it contains no Monte Carlo/physics
content and does not run or validate any of the tools it checks for.

It is standard-library only.

## Determinism

Unlike `environment_report.py` (which reports live host state) this
utility is **deterministic apart from the repository's current
branch/SHA**. Every capability check is a pure filesystem/repo-content
check: does a file exist, is it non-empty, and — where relevant — does a
specific marker (a JSON key in `.mcp.json`, a `name:` key in an agent's
YAML frontmatter) appear. Nothing shells out to a tool whose behavior
could vary by host (no `--version` probes), and no network access is
used or required. The only environment-derived, host-varying fields
anywhere in the report are `repository.branch` / `repository.sha` /
`repository.dirty`, plus `generated_at_utc` (which varies by
generation time, not by host).

## What it reports

Each top-level key other than `schema_version`, `generated_at_utc`, and
`summary` is a "capability section". Every section is a JSON object with
at least a `status` field (see [Status vocabulary](#status-vocabulary)).

| Key | Capability | Presence requires |
|---|---|---|
| `repository` | (not a capability — repo identity) | `branch`, `sha`, `dirty` for the Git repository containing this module, resolved from the module's own file location, never the process CWD. |
| `github_ci_workflow` | Lightweight GitHub CI | `.github/workflows/ci.yml` exists and is non-empty. |
| `pre_commit_config` | Pre-commit hooks | `.pre-commit-config.yaml` exists and is non-empty. |
| `local_exact_sha_validation` | Local exact-SHA validation gate | `infrastructure/validation/local_validation.py` and `infrastructure/mcp/validation_mcp.py` both non-empty, and `.mcp.json` declares an `mcpServers` entry named `validation-manager`. |
| `task_worktree_manager` | Isolated task worktrees | `infrastructure/tasks/task_worktree.py` and `infrastructure/mcp/task_mcp.py` both non-empty, and `.mcp.json` declares `task-manager`. |
| `autonomous_commit_helper` | Fixed-identity autonomous commits | `infrastructure/tasks/task_commit.py` and `infrastructure/mcp/commit_mcp.py` both non-empty, and `.mcp.json` declares `task-commit`. |
| `host_gpu_validation_runner` | Sandboxed host/GPU execution | `infrastructure/host_runner/host_runner.py` and `infrastructure/mcp/host_runner_mcp.py` both non-empty, and `.mcp.json` declares `host-runner`. |
| `codex_worker` | Codex worker delegation | `infrastructure/codex/codex_worker.py` and `infrastructure/mcp/codex_mcp.py` both non-empty, and `.mcp.json` declares `codex-worker`. |
| `github_ci_observability` | GitHub CI/PR inspection tools | `infrastructure/tasks/task_integration.py` and `infrastructure/mcp/integration_mcp.py` both non-empty, and `.mcp.json` declares `task-integration`. |
| `claude_native_subagents` | Native Claude subagents | `.claude/agents/` exists and contains at least one `*.md` file whose YAML frontmatter declares a `name:` key. Discovered names are listed in `agent_names`. |
| `governing_docs` | Governing documents | `EXPERIMENT.md`, `REQUIREMENTS.md`, `AGENTS.md`, `CLAUDE.md` all exist and are non-empty at the repository root (reported individually in `docs`, plus an overall `all_present`). |

For the six capabilities backed by `.mcp.json` (`local_exact_sha_validation`
through `github_ci_observability`), each section independently reports
`files` (which of the required files were found) and `mcp_registered`
(whether the named server is declared). `.mcp.json` is read once per
report and shared across all six checks; a missing, unreadable, or
malformed `.mcp.json` is treated as a hard fact affecting every one of
these capabilities — each affected section reports `mcp_registered: null`
plus a `detail` explaining why the MCP-registration part could not be
confirmed, rather than the report failing to be produced.

## Status vocabulary

This is a different domain from `environment_report.py`'s host-probe
vocabulary (`available`/`unavailable`/`error`): this module reports
*capability presence*, not *probe outcomes*. Every section carries a
`status` field drawn from:

- **`present`** — every file/marker required for this capability was
  found. The section's data fields are populated and no `detail` is
  required (though informational fields such as `files` are still
  included).
- **`partial`** — some, but not all, of the required files/markers for
  this capability were found (e.g. one of two required files is missing,
  or the files exist but the `.mcp.json` registration is missing). A
  `detail` field explains what is missing.
- **`absent`** — none of the required files/markers for this capability
  were found.
- **`error`** — checking the capability failed unexpectedly (e.g.
  `.mcp.json` exists but is not valid JSON, or a filesystem read raised
  for a reason other than "does not exist"). A `detail` field explains
  why.

Every section whose status is not `present` carries a `detail` string.

## Summary field

`summary` is a top-level (non-"section") object:

```jsonc
"summary": {
  "all_present": true,
  "capabilities_not_present": []
}
```

- `all_present` — `true` iff every capability section (i.e. every
  top-level key other than `schema_version`, `generated_at_utc`,
  `repository`, and `summary` itself) has `status == "present"`.
- `capabilities_not_present` — the sorted list of capability keys whose
  status is not `present`, so the report is scannable without
  inspecting every section individually.

`repository` is excluded from this rollup: it reports repo identity
(branch/SHA/dirty), not a capability, and is not part of the
present/absent audit.

## JSON schema outline

```jsonc
{
  "schema_version": 1,
  "generated_at_utc": "2026-09-10T00:00:00+00:00",
  "repository": { "status": "present", "branch": "...", "sha": "...", "dirty": false },
  "github_ci_workflow": { "status": "present", "files": { ".github/workflows/ci.yml": true } },
  "pre_commit_config": { "status": "present", "files": { ".pre-commit-config.yaml": true } },
  "local_exact_sha_validation": {
    "status": "present",
    "files": { "infrastructure/validation/local_validation.py": true, "infrastructure/mcp/validation_mcp.py": true },
    "mcp_server_name": "validation-manager",
    "mcp_registered": true
  },
  "task_worktree_manager": { "...": "same shape, mcp_server_name: task-manager" },
  "autonomous_commit_helper": { "...": "same shape, mcp_server_name: task-commit" },
  "host_gpu_validation_runner": { "...": "same shape, mcp_server_name: host-runner" },
  "codex_worker": { "...": "same shape, mcp_server_name: codex-worker" },
  "github_ci_observability": { "...": "same shape, mcp_server_name: task-integration" },
  "claude_native_subagents": { "status": "present", "agent_names": ["implementation-worker", "..."], "skipped_files": [] },
  "governing_docs": { "status": "present", "docs": { "EXPERIMENT.md": true, "REQUIREMENTS.md": true, "AGENTS.md": true, "CLAUDE.md": true }, "all_present": true },
  "summary": { "all_present": true, "capabilities_not_present": [] }
}
```

## No secrets

This utility touches no secrets at all. It never reads environment
variables, credentials, or remote URLs, and never dumps `os.environ`. It
only reads a small, fixed set of repository-relative file paths, plus
`.mcp.json`'s already version-controlled, non-secret server-registration
structure (server names and launch commands, not credentials).

## Running it

```bash
# Pretty JSON to stdout (default indent: 2)
python3 infrastructure/diagnostics/repo_capability_report.py

# Custom indent
python3 infrastructure/diagnostics/repo_capability_report.py --indent 4

# Write to a file instead of stdout
python3 infrastructure/diagnostics/repo_capability_report.py --output /path/to/report.json
```

### Exit code policy

- **`0`** — the report was produced (even if some capabilities are
  `absent`/`partial`/`error`) and successfully delivered to stdout or
  `--output`.
- **`1`** — the report was produced but could not be written to the
  `--output` path (e.g. a nonexistent parent directory or a permission
  error). A structured JSON object of the form
  `{"status": "error", "detail": "..."}` is printed to stderr; no
  traceback is ever printed.
- **`2`** — an argparse usage error (e.g. an unrecognized argument or bad
  `--indent` value), raised by argparse itself before any report is
  generated.

## Running its tests

The repository does not install pytest system-wide. Run tests with:

```bash
uv run --with pytest --no-project python -m pytest tests -q
```

`tests/conftest.py` inserts `infrastructure/diagnostics/` onto `sys.path`
explicitly, matching the repository convention that `infrastructure/`
subdirectories are plain script directories without `__init__.py`, rather
than installable packages.

Tests are deterministic and do not depend on incidental host/repo state:
each capability builder is exercised against a synthetic fake repository
root built under `tmp_path`, covering present/absent/partial/malformed-
`.mcp.json` cases. The only tests that touch the real checked-out
repository assert true invariants of this checkout (e.g. that the
governing docs and this utility's own files exist), not incidental
host state.
