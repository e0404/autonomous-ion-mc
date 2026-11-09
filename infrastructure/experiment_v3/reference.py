"""Controlled native reference-engine execution; success is NOT physics validation.

Archive of one run: ``<state>/references/<run_id>/`` with ``request.json``,
``inputs/``, ``work/`` (inputs plus native outputs), ``stdout.txt``,
``stderr.txt``, ``record.json`` (signed) and ``result.json`` (archive index).

``record.json`` -- signed reference-run record, schema_version 3
(V3-INFRA-004). Ed25519 attestation by the service-held ``reference-runner``
key (infrastructure/experiment_v3/attest.py); verified by
infrastructure/validation/reference_run_records.py. The whole record is signed,
including the bounded stdout/stderr tails. It is also what run_reference (and
the MCP tool run_reference_calculation) returns. Fields:

- schema_version (3), run_id, configuration_id (sha256 of the canonical
  request), code_sha (exact clean commit), task_id (None from the CLI), engine;
- request_sha256: sha256 of request.json exactly as written (binds engine
  runtime, case, input-file hashes, timeout, gpu, runner/sandbox hashes);
- started_at: UTC ISO-8601 time taken immediately before the controlled
  ``reference_started`` event (whose details carry the same value) and before
  the engine process is launched; finished_at: UTC time after the process ended
  and the outputs were inventoried; seconds: monotonic duration between them;
- exit_code, timed_out, error, missing_outputs, succeeded,
  scientifically_validated (always False), metadata_status;
- output_count and output_files_sha256 = common.digest(output_files), i.e.
  sha256 of the canonical JSON (sorted keys, compact separators) of the
  ``work/`` inventory {path: {sha256, bytes, ...}} listed in result.json;
- stdout, stderr: tails ({text, bytes, truncated}, at most 2000 bytes each);
  stdout_sha256, stderr_sha256: hashes of the full archived logs;
- attestation: see attest.py.

``result.json`` -- archive index (unsigned; carries no attestation block):
every record field except ``attestation``, plus ``output_files`` (inventory of
``work/``), ``artifact_files`` (inventory of the run directory, computed after
record.json was written, so record.json is listed with its sha256 and can be
listed, previewed and materialized by infrastructure/experiment_v3/artifacts.py)
and ``signed_record`` = "record.json". Only record.json is evidence; the
signature, not the index, establishes its integrity.
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import signal
import subprocess
import time
import uuid
from pathlib import Path

from infrastructure.experiment_v3 import attest
from infrastructure.experiment_v3.common import (
    STATE,
    digest,
    event,
    exact_state,
    file_hash,
    git,
    identifier,
    inventory,
    now,
    safe_path,
    write_json,
)
from infrastructure.experiment_v3.sandbox_runtime import require as require_bwrap

ENGINES = {"topas", "mcsquare", "fred"}
SCHEMA_VERSION = 3
SIGNER = "reference-runner"
INLINE_BYTES = 16000
MAX_SECONDS = 7200


def runtime_inventory(root):
    """Hash files and in-tree symlinks, including libraries and physics data."""
    root = Path(root).resolve()
    entries = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            if not path.resolve().is_relative_to(root):
                raise ValueError(f"Runtime symlink escapes installation: {path}")
            entries[str(path.relative_to(root))] = {"link": os.readlink(path)}
        elif path.is_file():
            entries[str(path.relative_to(root))] = {
                "sha256": file_hash(path),
                "mode": path.stat().st_mode & 0o777,
            }
    if not entries:
        raise ValueError("Empty runtime directory")
    return digest(entries)


def load_engine(name, registry):
    if name not in ENGINES:
        raise ValueError("Unsupported engine")
    config = json.loads(Path(registry).read_text())["engines"][name]
    if not config.get("version") or not config.get("source"):
        raise ValueError("Engine version and acquisition source are required")
    for mount in config["mounts"]:
        destination = mount["destination"]
        if (
            not destination.startswith("/opt/reference/")
            or ".." in Path(destination).parts
        ):
            raise ValueError("Runtime mounts must be below /opt/reference")
        if runtime_inventory(mount["source"]) != mount["sha256"]:
            raise ValueError(
                "Runtime contents changed; reprovision and record the change"
            )
    if not config["executable"].startswith("/opt/reference/"):
        raise ValueError("Engine executable must be in the registered runtime")
    return config


def validate_case(case):
    for key in (
        "physics",
        "source",
        "geometry",
        "materials",
        "cuts",
        "histories",
        "seeds",
        "scorers",
        "units",
        "rationale",
        "expected_outputs",
    ):
        if not case.get(key):
            raise ValueError(f"Case metadata missing: {key}")
    if type(case["histories"]) is not int or case["histories"] < 1:
        raise ValueError("histories must be a positive integer")
    if not isinstance(case["seeds"], list) or not all(
        type(x) is int and x > 0 for x in case["seeds"]
    ):
        raise ValueError("Explicit nonzero integer seeds are required")
    if not isinstance(case.get("arguments", []), list) or not all(
        isinstance(x, str) for x in case.get("arguments", [])
    ):
        raise ValueError("arguments must be an argv list")
    if not isinstance(case["expected_outputs"], list):
        raise ValueError("expected_outputs must be a list")
    safe_path(Path("/case"), case["input"])
    for output in case["expected_outputs"]:
        safe_path(Path("/case"), output)


def command(config, case, work, gpu=False):
    bwrap = require_bwrap()
    args = [bwrap, "--die-with-parent", "--new-session", "--unshare-all"]
    for path in (
        "/usr",
        "/bin",
        "/lib",
        "/lib64",
        "/etc/ld.so.cache",
        "/etc/alternatives",
    ):
        if Path(path).exists():
            args += ["--ro-bind", path, path]
    args += [
        "--proc",
        "/proc",
        "--dev",
        "/dev",
        "--tmpfs",
        "/tmp",
        "--dir",
        "/tmp/home",
    ]
    for mount in config["mounts"]:
        args += ["--ro-bind", mount["source"], mount["destination"]]
    if gpu:
        if not Path("/dev/dxg").exists():
            raise RuntimeError("This runner supports WSL GPU access via /dev/dxg only")
        args += ["--dev-bind", "/dev/dxg", "/dev/dxg", "--ro-bind", "/sys", "/sys"]
    args += ["--bind", str(work), "/work", "--chdir", "/work", "--clearenv"]
    env = {
        "HOME": "/tmp/home",
        "PATH": "/usr/bin:/bin",
        "LANG": "C.UTF-8",
        "TMPDIR": "/tmp",
        "OMP_NUM_THREADS": "1",
        **config.get("environment", {}),
    }
    for key, value in env.items():
        args += ["--setenv", key, value]
    argv = [config["executable"]]
    if config.get("driver"):
        argv.append(config["driver"])
    if config["engine"] == "fred":
        argv += ["-f", case["input"]]
    else:
        argv.append(case["input"])
    return [*args, "--", *argv, *case.get("arguments", [])]


def limits():
    # Bound raw/log files and core dumps; wall timeout kills the process group.
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_FSIZE, (512 * 1024**2, 512 * 1024**2))


def tail(path):
    with path.open("rb") as stream:
        size = path.stat().st_size
        stream.seek(max(0, size - 2000))
        return {
            "text": stream.read().decode(errors="replace"),
            "bytes": size,
            "truncated": size > 2000,
        }


def require_committed_inputs(root, case_path, files):
    tracked = set(git(root, "ls-files", "--", str(case_path)).splitlines())
    expected = {(Path(case_path) / name).as_posix() for name in files}
    if not expected or not expected <= tracked:
        raise ValueError("Every reference input must be tracked in the clean commit")


def run_reference(
    root,
    case_path,
    engine,
    *,
    registry=None,
    state=None,
    timeout=600,
    gpu=False,
    task_id=None,
):
    if type(timeout) is not int or not 1 <= timeout <= MAX_SECONDS:
        raise ValueError("timeout must be 1..7200 seconds")
    root, state = Path(root), Path(state or STATE)
    sha = exact_state(root)
    case_dir = safe_path(root, case_path)
    files = inventory(case_dir)
    require_committed_inputs(root, case_path, files)
    case = json.loads((case_dir / "case.json").read_text())
    validate_case(case)
    if not safe_path(case_dir, case["input"]).is_file():
        raise ValueError("Native input does not exist")
    config = load_engine(engine, registry or STATE / "engines.json")
    if config["engine"] != engine:
        raise ValueError("Registry engine mismatch")
    request = {
        "schema_version": 2,
        "code_sha": sha,
        "dirty": False,
        "engine": config,
        "case": case,
        "input_files": files,
        "gpu": gpu,
        "timeout_seconds": timeout,
        "runner_source_sha256": file_hash(Path(__file__)),
        "os_release_sha256": file_hash("/etc/os-release"),
        "sandbox_binary_sha256": (file_hash(require_bwrap())),
    }
    config_id = digest(request)
    run_id = "REF-" + config_id[:20] + "-" + uuid.uuid4().hex[:8]
    run_dir = state / "references" / run_id
    run_dir.mkdir(parents=True)
    shutil.copytree(case_dir, run_dir / "inputs")
    shutil.copytree(run_dir / "inputs", run_dir / "work")
    if inventory(run_dir / "inputs") != files:
        raise ValueError("Inputs changed while snapshotting")
    for output in case["expected_outputs"]:
        if safe_path(run_dir / "work", output).exists():
            raise ValueError("Expected output already exists in input bundle")
    write_json(run_dir / "request.json", request)
    request_sha256 = file_hash(run_dir / "request.json")
    # V3-INFRA-004: the controlled start time, taken at the reference_started
    # event and signed in the record (chronology gates must not derive it from
    # finished_at - seconds).
    started_at, started = now(), time.monotonic()
    event(
        "reference_started",
        sha=sha,
        task_id=task_id,
        details={
            "run_id": run_id,
            "engine": engine,
            "configuration_id": config_id,
            "started_at": started_at,
        },
        state=state,
    )
    timed_out, error, returncode = False, None, None
    with (
        (run_dir / "stdout.txt").open("wb") as out,
        (run_dir / "stderr.txt").open("wb") as err,
    ):
        try:
            proc = subprocess.Popen(
                command(config, case, run_dir / "work", gpu),
                stdout=out,
                stderr=err,
                start_new_session=True,
                preexec_fn=limits,
                env={"PATH": "/usr/bin:/bin"},
            )
            try:
                returncode = proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
        except (OSError, RuntimeError) as exc:
            error = type(exc).__name__ + ": " + str(exc)
    try:
        outputs = inventory(run_dir / "work", allow_internal_links=True)
    except ValueError as exc:
        outputs, error = {}, str(exc)
    missing = [
        name
        for name in case["expected_outputs"]
        if outputs.get(name, {}).get("bytes", 0) == 0
    ]
    record = {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "configuration_id": config_id,
        "code_sha": sha,
        "task_id": task_id,
        "engine": engine,
        "request_sha256": request_sha256,
        "started_at": started_at,
        "finished_at": now(),
        "seconds": time.monotonic() - started,
        "exit_code": returncode,
        "timed_out": timed_out,
        "error": error,
        "missing_outputs": missing,
        "succeeded": returncode == 0 and not timed_out and not error and not missing,
        "scientifically_validated": False,
        "metadata_status": "declared; inspect native logs/outputs",
        "output_count": len(outputs),
        "output_files_sha256": digest(outputs),
        "stdout": tail(run_dir / "stdout.txt"),
        "stderr": tail(run_dir / "stderr.txt"),
        "stdout_sha256": file_hash(run_dir / "stdout.txt"),
        "stderr_sha256": file_hash(run_dir / "stderr.txt"),
    }
    record = attest.sign_with_service_key(record, SIGNER)
    write_json(run_dir / "record.json", record)
    # The index is written after record.json so that the artifact inventory
    # lists the signed record (with its sha256) for the artifact tools.
    result = {k: v for k, v in record.items() if k != "attestation"} | {
        "output_files": outputs,
        "artifact_files": inventory(run_dir, allow_internal_links=True),
        "signed_record": "record.json",
    }
    write_json(run_dir / "result.json", result)
    event(
        "reference_finished" if result["succeeded"] else "reference_failed",
        sha=sha,
        task_id=task_id,
        details={"run_id": run_id, "engine": engine, "succeeded": result["succeeded"]},
        state=state,
    )
    # The signed record is bounded (inventories only by digest/count); the full
    # file inventories are in result.json.
    return record


def read_artifact(run_id, relative, *, state=None, offset=0, limit=16000):
    identifier(run_id)
    if not 0 <= offset or not 1 <= limit <= INLINE_BYTES:
        raise ValueError("Invalid bounded read")
    path = safe_path(Path(state or STATE) / "references" / run_id, relative)
    with path.open("rb") as stream:
        stream.seek(offset)
        data = stream.read(limit)
    import base64

    return {
        "base64": base64.b64encode(data).decode(),
        "offset": offset,
        "bytes": len(data),
        "total_bytes": path.stat().st_size,
        "sha256": file_hash(path),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path.cwd())
    p.add_argument("--case", required=True)
    p.add_argument("--engine", choices=sorted(ENGINES), required=True)
    p.add_argument("--registry", type=Path)
    p.add_argument("--state", type=Path)
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--gpu", action="store_true")
    a = p.parse_args()
    result = run_reference(
        a.root,
        a.case,
        a.engine,
        registry=a.registry,
        state=a.state,
        timeout=a.timeout,
        gpu=a.gpu,
    )
    print(json.dumps(result, indent=2))
    return 0 if result["succeeded"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
