# Environment Diagnostics Utility

`infrastructure/diagnostics/environment_report.py` reports runtime/environment
information relevant to reproducing or debugging an experiment run, as a
single machine-readable JSON document.

It is standard-library only. It never requires `warp`, `nvidia-smi`,
`claude`, `codex`, or any other optional tool to be present — every probe
degrades gracefully and the report is always produced.

## What it reports

Each top-level key is a "section". Every section is a JSON object with at
least a `status` field (see [Status vocabulary](#status-vocabulary)).

| Key | Contents (where available) |
|---|---|
| `schema_version` | Integer report schema version (currently `1`). Not a section - always present, no `status`. |
| `generated_at_utc` | ISO-8601, timezone-aware UTC timestamp of report generation. Not a section. |
| `os` | `platform`, `release`, `version`, `machine` (from `platform`), plus `distro` parsed from `/etc/os-release` when present. |
| `cpu` | `architecture`, `logical_cores` (`os.cpu_count()`), `model_name` and `physical_cores` parsed from `/proc/cpuinfo` on Linux. |
| `memory` | `total_bytes`, `available_bytes` parsed from `/proc/meminfo` on Linux. |
| `python` | `version`, `version_info`, `implementation`, `executable`, `prefix`. Always `available`. |
| `uv` | `version`, `raw` from `uv --version`. |
| `git` | `version`, `raw` from `git --version`. |
| `repository` | `branch`, `sha`, `dirty` for the Git repository containing this module, resolved from the module's own file location - **not** the process's current working directory. |
| `claude_code` | `version`, `raw` from `claude --version`. |
| `codex_cli` | `version`, `raw` from `codex --version`. |
| `nvidia_gpu` | `gpus`: list of `{name, driver_version, memory_total}` from `nvidia-smi --query-gpu=... --format=csv,noheader`. |
| `warp` | `version`, `cuda_available` (bool), `cuda_devices` (list of device names), from an optional `import warp`. |

## Status vocabulary

Every section carries a `status` field using a fixed, small vocabulary:

- **`available`** - the probe succeeded; the section's data fields are populated.
- **`unavailable`** - the underlying tool/file/module is simply absent (executable
  not on `PATH`, optional Python module not importable, platform-specific file
  such as `/proc/cpuinfo` or `/etc/os-release` does not exist). This is an
  expected condition, not a failure, and is common on machines without a GPU,
  without `warp`, or without the `claude`/`codex` CLIs installed.
- **`error`** - the probe was attempted but failed unexpectedly (command timed
  out, exited non-zero, produced unparseable output, or an unexpected
  exception was raised). A `detail` field explains why.

A section with `status != "available"` always carries a `detail` string and
omits or nulls out its data fields. No probe failure ever propagates out of
`build_report()` - a report is always produced, and the CLI always exits `0`
on a successful report, even when every optional probe is `unavailable` or
`error`.

## JSON schema outline

```jsonc
{
  "schema_version": 1,
  "generated_at_utc": "2026-09-09T18:32:07.629422+00:00",
  "os": { "status": "available", "platform": "Linux", "release": "...", "version": "...", "machine": "x86_64", "distro": { "...": "..." } },
  "cpu": { "status": "available", "architecture": "x86_64", "logical_cores": 32, "model_name": "...", "physical_cores": 16 },
  "memory": { "status": "available", "total_bytes": 33516343296, "available_bytes": 28874686464 },
  "python": { "status": "available", "version": "...", "version_info": [3, 12, 3, "final", 0], "implementation": "CPython", "executable": "...", "prefix": "..." },
  "uv": { "status": "available", "version": "0.12.11", "raw": "uv 0.12.11 (...)" },
  "git": { "status": "available", "version": "2.43.0", "raw": "git version 2.43.0" },
  "repository": { "status": "available", "branch": "...", "sha": "...", "dirty": false },
  "claude_code": { "status": "unavailable", "detail": "executable not found: claude" },
  "codex_cli": { "status": "unavailable", "detail": "executable not found: codex" },
  "nvidia_gpu": { "status": "unavailable" | "error", "detail": "...", "gpus": [ { "name": "...", "driver_version": "...", "memory_total": "..." } ] },
  "warp": { "status": "unavailable" | "available", "detail": "...", "version": "...", "cuda_available": true, "cuda_devices": ["..."] }
}
```

## No secrets, no host-identifying personal data

The report never includes environment variables, user names, home paths of
the invoking user, tokens, or remote URLs with credentials. It does not read
or dump `os.environ`. `python.executable` / `python.prefix` are the only
fields that may incidentally contain a filesystem path (they are explicitly
required diagnostic fields); no field is a raw environment or account dump.

## Running it

```bash
# Pretty JSON to stdout (default indent: 2)
python3 infrastructure/diagnostics/environment_report.py

# Custom indent
python3 infrastructure/diagnostics/environment_report.py --indent 4

# Write to a file instead of stdout
python3 infrastructure/diagnostics/environment_report.py --output /path/to/report.json
```

Exit code is always `0` on a successful report (including when every
optional probe is `unavailable`/`error`).

## Running its tests

The repository does not install pytest system-wide. Run tests with:

```bash
uv run --with pytest --no-project python -m pytest tests -q
```

`tests/conftest.py` inserts `infrastructure/diagnostics/` onto `sys.path`
explicitly, matching the repository convention that `infrastructure/`
subdirectories (e.g. `infrastructure/tasks/`, `infrastructure/validation/`)
are plain script directories without `__init__.py`, rather than installable
packages.

Tests are deterministic and do not assume any optional tool (GPU, `warp`,
`claude`, `codex`) is present on the host. Graceful-degradation paths
(missing executable, timeout, non-zero exit, unexpected exception) are
exercised by monkeypatching `subprocess.run`, not by depending on host state.
