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
| `os` | `platform`, `release`, `version`, `machine` (from `platform`), plus `distro` parsed from `/etc/os-release` when present, else `null`. Always `available` (see below). |
| `cpu` | `architecture`, `logical_cores` (`os.cpu_count()`), plus `model_name` and `physical_cores` parsed from `/proc/cpuinfo` when present, else `null`. Always `available` (see below). |
| `memory` | `total_bytes`, `available_bytes` parsed from `/proc/meminfo` on Linux. |
| `python` | `version`, `version_info`, `implementation`, `executable`, `prefix`. Always `available`. |
| `uv` | `version`, `raw` from `uv --version`. |
| `git` | `version`, `raw` from `git --version`. |
| `repository` | `branch`, `sha`, `dirty` for the Git repository containing this module, resolved from the module's own file location - **not** the process's current working directory. |
| `claude_code` | `version`, `raw` from `claude --version`. |
| `codex_cli` | `version`, `raw` from `codex --version`. |
| `nvidia_gpu` | `gpus`: list of `{name, driver_version, memory_total}` from `nvidia-smi --query-gpu=... --format=csv,noheader`. An empty result (driver present, no GPUs enumerated) is reported as `unavailable` with `gpus: []`, not as `available` with an empty list. |
| `warp` | `version`, `cuda_available` (bool), `cuda_devices` (list of device names), from an optional `import warp`. |

`/etc/os-release` parsing (`parse_os_release`) is best-effort: it strips
matching surrounding quotes from each value but does not interpret
backslash escapes within them.

### `os` and `cpu` are always `available`

Unlike the CLI-tool probes, `os` and `cpu` are always reported as
`available`, even on a non-Linux host or one without `/proc/cpuinfo` or
`/etc/os-release`: `platform.system()`/`platform.machine()`/`os.cpu_count()`
always succeed, so there is always real data to report. In that case the
Linux-specific fields (`os.distro`, `cpu.model_name`, `cpu.physical_cores`)
are simply `null` rather than the section itself being `unavailable`.

## Status vocabulary

Every section carries a `status` field using a fixed, small vocabulary:

- **`available`** - the probe succeeded; the section's data fields are populated.
- **`unavailable`** - the underlying tool/file/module is simply absent (executable
  not on `PATH`, optional Python module not importable) or produced no data to
  report (e.g. `nvidia-smi` ran but enumerated zero GPUs). This is an expected
  condition, not a failure, and is common on machines without a GPU, without
  `warp`, or without the `claude`/`codex` CLIs installed. It is distinct from
  a tool being *present but failing*, which is `error` (see below) - for
  example, a `PermissionError` running an executable is `error`, not
  `unavailable`, because the executable was found and running it was actually
  attempted.
- **`error`** - the probe was attempted but failed unexpectedly (permission
  denied, command timed out, exited non-zero, produced unparseable or
  malformed output, or an unexpected exception was raised). A `detail` field
  explains why.

A section with `status != "available"` always carries a `detail` string and
omits or nulls out its other data fields. No probe failure ever propagates
out of `build_report()` - a report is always produced. `os`, `cpu`, and
`python` are always `available` (see below); every other section can be
`unavailable` or `error` independently of the others.

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
  "nvidia_gpu": { "status": "unavailable" | "error" | "available", "detail": "...", "gpus": [ { "name": "...", "driver_version": "...", "memory_total": "..." } ] },
  "warp": { "status": "unavailable" | "available", "detail": "...", "version": "...", "cuda_available": true, "cuda_devices": ["..."] }
}
```

## No secrets, no host-identifying personal data

The report never includes environment variables, tokens, or remote URLs with
credentials, and it does not read or dump `os.environ`. `python.executable`
and `python.prefix` are the two exceptions to "no user-specific filesystem
paths": they are explicitly required diagnostic fields, and on a typical
install they genuinely can contain a path under the invoking user's home
directory (e.g. a venv or user-site install under `/home/<user>/...`). No
other field contains a filesystem path, a user name, or any other
account-identifying data.

## Running it

```bash
# Pretty JSON to stdout (default indent: 2)
python3 infrastructure/diagnostics/environment_report.py

# Custom indent
python3 infrastructure/diagnostics/environment_report.py --indent 4

# Write to a file instead of stdout
python3 infrastructure/diagnostics/environment_report.py --output /path/to/report.json
```

### Exit code policy

- **`0`** - the report was produced (even if every optional probe is
  `unavailable`/`error`) and successfully delivered to stdout or `--output`.
- **`1`** - the report was produced but could not be written to the
  `--output` path (e.g. a nonexistent parent directory or a permission
  error). A structured JSON object of the form
  `{"status": "error", "detail": "..."}` is printed to stderr; no traceback
  is ever printed.
- **`2`** - an argparse usage error (e.g. an unrecognized argument or bad
  `--indent` value), raised by argparse itself before any report is
  generated.

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
