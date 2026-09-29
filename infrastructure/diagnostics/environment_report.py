#!/usr/bin/env python3
"""Experiment-environment diagnostics utility.

Collects runtime/environment information relevant to reproducing or
debugging experiment runs (OS, CPU, memory, Python, common CLI tool
versions, the enclosing Git repository state, GPU/driver info, and
optional NVIDIA Warp availability) and reports it as a single
machine-readable JSON document.

Design goals (see also docs/environment_diagnostics.md):

* Standard library only. `warp` is probed via an optional import and is
  never a hard dependency of this module.
* Every probe is independently guarded: a missing or broken optional
  tool degrades to a section with a non-"available" status rather than
  raising. The top-level report is always produced; the CLI exits 0
  whenever a report was produced and delivered, even when every
  optional probe is unavailable (see the CLI exit-code policy on
  `main()` and in the docs for the other exit codes).
* No secrets or host-identifying personal data: no environment-variable
  dumps, user names, tokens, or remote URLs are ever included. The only
  exception is `python.executable`/`python.prefix`, which are required
  diagnostic fields that may incidentally contain a path under the
  invoking user's home directory.

Status vocabulary used by every section's "status" field:

* "available"   - the probe succeeded; the section's data fields are populated.
* "unavailable" - the underlying tool/file/module is simply absent (e.g.
                  the executable is not on PATH, or an optional Python
                  module is not importable) or produced no data (e.g.
                  nvidia-smi ran but enumerated zero GPUs). This is an
                  expected condition, not a failure.
* "error"       - the probe was attempted but failed unexpectedly (e.g.
                  permission denied, the command timed out, exited
                  non-zero, or its output could not be parsed or was
                  malformed). A "detail" field explains why.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

SCHEMA_VERSION = 1

# A few seconds is enough for any of the CLI probes below; keeps a hung
# or misbehaving tool from ever blocking report generation.
SUBPROCESS_TIMEOUT_SECONDS = 5.0

STATUS_AVAILABLE = "available"
STATUS_UNAVAILABLE = "unavailable"
STATUS_ERROR = "error"

_VERSION_TOKEN_RE = re.compile(r"\d+(?:\.\d+){1,3}")


# ---------------------------------------------------------------------------
# Subprocess helper
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
    ``stdout``; on failure it carries ``status`` (one of
    STATUS_UNAVAILABLE / STATUS_ERROR) and ``detail``.
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
            "status": STATUS_UNAVAILABLE,
            "detail": f"executable not found: {args[0]}",
        }
    except PermissionError as exc:
        # The executable exists but could not be run: this was attempted
        # and failed, not simply absent, so it is an "error" rather than
        # "unavailable".
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


def extract_version_token(text: str) -> str | None:
    """Pull the first dotted-numeric version token out of free-form text."""
    match = _VERSION_TOKEN_RE.search(text)
    return match.group(0) if match else None


def version_section(
    args: list[str],
    parse: Callable[[str], dict[str, Any]] | None = None,
    timeout: float = SUBPROCESS_TIMEOUT_SECONDS,
) -> dict[str, Any]:
    """Build a section dict for a "<tool> --version"-style probe."""
    result = run_subprocess(args, timeout=timeout)
    if not result["ok"]:
        return {"status": result["status"], "detail": result["detail"]}

    stdout = result["stdout"]
    try:
        extra = parse(stdout) if parse else {"version": extract_version_token(stdout), "raw": stdout}
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": f"failed to parse output: {exc}"}

    return {"status": STATUS_AVAILABLE, **extra}


def _parse_labeled_version(stdout: str) -> dict[str, Any]:
    version = extract_version_token(stdout)
    if version is None:
        # Unparseable output must surface as an "error" section, not as
        # "available" with a null version - version_section() converts
        # this exception into a STATUS_ERROR result.
        raise ValueError(f"could not find a version number in output: {stdout!r}")
    return {"version": version, "raw": stdout}


# ---------------------------------------------------------------------------
# Parsing helpers (pure functions, unit-testable against fixture strings)
# ---------------------------------------------------------------------------


def parse_os_release(text: str) -> dict[str, str]:
    """Parse the contents of an /etc/os-release-style file.

    Best-effort: this strips matching surrounding single/double quotes
    from each value (the common case for this file format), but it does
    NOT interpret shell-style backslash escapes within quoted values.
    """
    data: dict[str, str] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
            value = value[1:-1]
        data[key.strip()] = value
    return data


def parse_cpuinfo(text: str) -> dict[str, Any]:
    """Parse the contents of Linux /proc/cpuinfo.

    Returns the CPU model name (from the first entry that has one) and
    the number of distinct (physical id, core id) pairs as the physical
    core count, when that information is present.
    """
    entries: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            if current:
                entries.append(current)
                current = {}
            continue
        if ":" in line:
            key, _, value = line.partition(":")
            current[key.strip()] = value.strip()
    if current:
        entries.append(current)

    model_name = None
    physical_core_ids: set[tuple[str, str]] = set()
    for entry in entries:
        if model_name is None and entry.get("model name"):
            model_name = entry["model name"]
        if "physical id" in entry and "core id" in entry:
            physical_core_ids.add((entry["physical id"], entry["core id"]))

    return {
        "model_name": model_name,
        "physical_cores": len(physical_core_ids) if physical_core_ids else None,
        "logical_entries": len(entries),
    }


_MEMINFO_UNIT_MULTIPLIERS = {"b": 1, "kb": 1024, "mb": 1024**2, "gb": 1024**3}


def parse_meminfo(text: str) -> dict[str, Any]:
    """Parse the contents of Linux /proc/meminfo into byte counts."""
    values: dict[str, int] = {}
    for line in text.splitlines():
        if ":" not in line:
            continue
        key, _, rest = line.partition(":")
        parts = rest.strip().split()
        if not parts:
            continue
        try:
            amount = int(parts[0])
        except ValueError:
            continue
        unit = parts[1].lower() if len(parts) > 1 else "kb"
        multiplier = _MEMINFO_UNIT_MULTIPLIERS.get(unit, 1024)
        values[key.strip()] = amount * multiplier

    return {
        "total_bytes": values.get("MemTotal"),
        "available_bytes": values.get("MemAvailable", values.get("MemFree")),
    }


def parse_nvidia_smi_csv(text: str) -> list[dict[str, str]]:
    """Parse `nvidia-smi --query-gpu=... --format=csv,noheader` output.

    Uses `csv.reader` so quoted fields (e.g. a GPU name containing a
    comma) are handled correctly. Each non-blank row must contain
    exactly 3 non-empty fields (name, driver_version, memory_total);
    a row that doesn't is malformed and raises ValueError so the
    caller can surface it as an "error" section rather than silently
    padding it with nulls.
    """
    gpus: list[dict[str, str]] = []
    for row in csv.reader(text.splitlines()):
        fields = [field.strip() for field in row]
        if not fields or all(field == "" for field in fields):
            continue
        if len(fields) != 3 or any(field == "" for field in fields):
            raise ValueError(f"expected 3 non-empty CSV fields (name, driver_version, memory.total), got: {row!r}")
        name, driver_version, memory_total = fields
        gpus.append(
            {
                "name": name,
                "driver_version": driver_version,
                "memory_total": memory_total,
            }
        )
    return gpus


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------


def build_os_section() -> dict[str, Any]:
    try:
        info: dict[str, Any] = {
            "status": STATUS_AVAILABLE,
            "platform": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
            "distro": None,
        }
        os_release_path = Path("/etc/os-release")
        if os_release_path.exists():
            try:
                info["distro"] = parse_os_release(os_release_path.read_text(encoding="utf-8"))
            except Exception as exc:  # noqa: BLE001
                info["distro"] = None
                info["distro_error"] = str(exc)
        return info
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": str(exc)}


def build_cpu_section() -> dict[str, Any]:
    try:
        info: dict[str, Any] = {
            "status": STATUS_AVAILABLE,
            "architecture": platform.machine(),
            "logical_cores": os.cpu_count(),
            "model_name": None,
            "physical_cores": None,
        }
        cpuinfo_path = Path("/proc/cpuinfo")
        if cpuinfo_path.exists():
            try:
                parsed = parse_cpuinfo(cpuinfo_path.read_text(encoding="utf-8"))
                info["model_name"] = parsed.get("model_name")
                info["physical_cores"] = parsed.get("physical_cores")
            except Exception as exc:  # noqa: BLE001
                info["cpuinfo_error"] = str(exc)
        return info
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": str(exc)}


def build_memory_section() -> dict[str, Any]:
    meminfo_path = Path("/proc/meminfo")
    if not meminfo_path.exists():
        return {
            "status": STATUS_UNAVAILABLE,
            "detail": "/proc/meminfo not present (non-Linux platform?)",
        }
    try:
        parsed = parse_meminfo(meminfo_path.read_text(encoding="utf-8"))
        return {"status": STATUS_AVAILABLE, **parsed}
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": str(exc)}


def build_python_section() -> dict[str, Any]:
    try:
        return {
            "status": STATUS_AVAILABLE,
            "version": sys.version,
            "version_info": list(sys.version_info),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
            "prefix": sys.prefix,
        }
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": str(exc)}


def build_repository_section(
    module_path: Path | str | None = None,
    timeout: float = SUBPROCESS_TIMEOUT_SECONDS,
) -> dict[str, Any]:
    """Report the branch/SHA/dirty state of the repo containing this module.

    Resolved from the file location of this module (``module_path``
    defaults to ``__file__``), never from the process's current working
    directory, so the result is correct regardless of where the CLI is
    invoked from.
    """
    module_dir = Path(module_path if module_path is not None else __file__).resolve().parent

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
        "status": STATUS_AVAILABLE,
        "branch": branch_result["stdout"],
        "sha": sha_result["stdout"],
        "dirty": bool(status_result["stdout"]),
    }


def build_uv_section() -> dict[str, Any]:
    return version_section(["uv", "--version"], parse=_parse_labeled_version)


def build_git_section() -> dict[str, Any]:
    return version_section(["git", "--version"], parse=_parse_labeled_version)


def build_claude_code_section() -> dict[str, Any]:
    return version_section(["claude", "--version"], parse=_parse_labeled_version)


def build_codex_cli_section() -> dict[str, Any]:
    return version_section(["codex", "--version"], parse=_parse_labeled_version)


def build_nvidia_gpu_section() -> dict[str, Any]:
    result = run_subprocess(
        [
            "nvidia-smi",
            "--query-gpu=name,driver_version,memory.total",
            "--format=csv,noheader",
        ]
    )
    if not result["ok"]:
        return {"status": result["status"], "detail": result["detail"]}

    try:
        gpus = parse_nvidia_smi_csv(result["stdout"])
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": f"failed to parse nvidia-smi output: {exc}"}

    if not gpus:
        return {
            "status": STATUS_UNAVAILABLE,
            "detail": "nvidia-smi ran successfully but reported no GPUs",
            "gpus": [],
        }

    return {"status": STATUS_AVAILABLE, "gpus": gpus}


@contextlib.contextmanager
def _silence_native_output():
    """Redirect file descriptors 1 and 2 to /dev/null for the duration.

    `warp.init()` prints an initialisation banner to stdout and, on hosts
    without a CUDA driver, error text to stderr from native code. Both
    would corrupt the JSON this tool emits on those streams, so the first
    Warp call is made with the descriptors redirected. Python-level
    `sys.stdout`/`sys.stderr` are flushed first so buffered report text is
    not lost.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        for fd in (devnull, *saved):
            os.close(fd)


def _quiet_warp_init(warp_module) -> None:
    """Initialise Warp without letting it write to stdout/stderr.

    Warp >= 1.17 suppresses its banner and module-load notices when
    `config.log_level` is raised to `LOG_WARNING`; the older `config.quiet`
    flag is deprecated there and itself prints a warning to stderr, so it is
    only used as a fallback when `log_level` does not exist.
    """
    try:
        if hasattr(warp_module.config, "log_level"):
            warp_module.config.log_level = warp_module.LOG_WARNING
        else:
            warp_module.config.quiet = True
    except Exception:  # noqa: BLE001
        pass
    with _silence_native_output():
        warp_module.init()


def build_warp_section() -> dict[str, Any]:
    """Probe for NVIDIA Warp without ever making it a hard dependency.

    Everything here is best-effort and defensively guarded: an
    unavailable/broken `warp` install must degrade to a section with a
    non-"available" status, never propagate, and never trigger an
    eager call that could plausibly crash the interpreter (e.g. no
    unconditional `warp.init()`).
    """
    try:
        import warp  # type: ignore  # noqa: PLC0415
    except ImportError as exc:
        return {"status": STATUS_UNAVAILABLE, "detail": f"warp is not importable: {exc}"}
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": f"unexpected error importing warp: {exc}"}

    try:
        _quiet_warp_init(warp)
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": f"warp.init() failed: {exc}"}

    info: dict[str, Any] = {"status": STATUS_AVAILABLE, "version": None, "cuda_available": None, "cuda_devices": None}

    try:
        info["version"] = getattr(warp, "__version__", None)
    except Exception as exc:  # noqa: BLE001
        info["version_error"] = str(exc)

    try:
        device_count = 0
        if hasattr(warp, "get_cuda_device_count"):
            device_count = int(warp.get_cuda_device_count())
        info["cuda_available"] = device_count > 0

        names: list[str] = []
        if device_count and hasattr(warp, "get_cuda_device"):
            for index in range(device_count):
                try:
                    device = warp.get_cuda_device(index)
                    names.append(str(getattr(device, "name", device)))
                except Exception as exc:  # noqa: BLE001
                    names.append(f"<error reading device {index}: {exc}>")
        info["cuda_devices"] = names
    except Exception as exc:  # noqa: BLE001
        info["cuda_available"] = None
        info["cuda_devices"] = None
        info["cuda_probe_error"] = str(exc)

    return info


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------

_SECTION_BUILDERS: dict[str, Callable[[], dict[str, Any]]] = {
    "os": build_os_section,
    "cpu": build_cpu_section,
    "memory": build_memory_section,
    "python": build_python_section,
    "uv": build_uv_section,
    "git": build_git_section,
    "repository": build_repository_section,
    "claude_code": build_claude_code_section,
    "codex_cli": build_codex_cli_section,
    "nvidia_gpu": build_nvidia_gpu_section,
    "warp": build_warp_section,
}


def _safe_section(builder: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    """Outer safety net around a section builder.

    Every builder above already guards its own failure modes, but this
    ensures a report can never fail to be produced even if a builder is
    later changed to raise unexpectedly.
    """
    try:
        result = builder()
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": f"unexpected error: {exc}"}

    if not isinstance(result, dict) or "status" not in result:
        return {"status": STATUS_ERROR, "detail": "probe returned a malformed result"}

    return result


def build_report() -> dict[str, Any]:
    """Build the full environment diagnostics report."""
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    for name, builder in _SECTION_BUILDERS.items():
        report[name] = _safe_section(builder)
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="environment_report",
        description="Report experiment environment diagnostics as JSON.",
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

    * 0 - the report was produced (even if every optional probe is
      "unavailable"/"error") and successfully delivered to stdout or
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
