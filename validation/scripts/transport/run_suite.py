"""Fail-closed runner of the V3-003 local-validation (LV) and host-runner (HR) suites.

Usage (argv only, no shell; the host runner executes exactly this)::

    python validation/scripts/transport/run_suite.py --suite {lv,hr} \
        --out validation/generated/transport/<new-dir> --expected-sha <40 hex> \
        [--workers N] [--step-timeout SECONDS] [--scale F]

Rules (the style of ``validation/scripts/warp-architecture/run_all.sh``):

* the output directory must not exist and must lie under ``validation/generated/`` or
  ``benchmarks/generated/`` of this repository (both are git-ignored); it is never reused;
* ``--expected-sha`` is mandatory (40 hex digits). If ``git`` can read this repository it must
  equal ``git rev-parse HEAD``; if the tree is a snapshot without ``.git`` the declared SHA is
  recorded as ``sha_source=declared``. A mismatch or a malformed SHA exits non-zero before any step;
* ``environment.txt`` records versions, hardware, the SHA, the dirty state and the sha256 of every
  script and transport source file;
* every step archives stdout+stderr in ``NN-name.txt`` between a header (command, SHA, start
  time, timeout) and an ``# exit=`` trailer; a timed-out step is killed and archived with
  ``exit=124``; ``manifest.txt`` lists the fixed step names of the suite;
* ``summarize.py`` verifies the manifest, the headers and the SHA of every file, then writes
  ``summary.json`` with one verdict per step; any failure, missing file or mismatch exits non-zero.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ALLOWED_PARENTS = (REPO / "validation" / "generated", REPO / "benchmarks" / "generated")
STEPS = HERE / "steps.py"
PY = sys.executable


def pytest_cmd(*targets: str, marker: str | None = None) -> list[str]:
    cmd = [PY, "-m", "pytest", "-q", "-p", "no:cacheprovider", "--no-header", *targets]
    return cmd + (["-m", marker] if marker else [])


def suite_steps(
    suite: str, workers: int, scale: float
) -> list[tuple[str, list[str], dict[str, str]]]:
    """``(name, argv, extra environment)`` of the steps of a suite, in execution order."""

    def n(x: int) -> str:
        return str(max(1000, int(x * scale)))

    st = [PY, str(STEPS)]
    w = ["--workers", str(workers)]
    if suite == "lv":
        return [
            (
                "01-pytest-warp-cpu-t2-t4-t11-c1-t13",
                pytest_cmd(
                    "tests/ionmc/test_transport_warp.py",
                    "tests/ionmc/test_transport_partition.py",
                    "tests/ionmc/test_config_validation.py",
                ),
                {},
            ),
            (
                "02-t1-trace-parity-256x150MeV",
                [*st, "t1", "--k", str(max(4, int(256 * scale))), "--energy", "150", *w],
                {},
            ),
            ("03-t2-deterministic-csda", [*st, "t2", "--backend", "warp-cpu"], {}),
            ("04-t13-workers", [*st, "t13", "--mode", "workers", "--n", n(200_000), *w], {}),
            (
                "05-t13-chunks-cpu",
                [*st, "t13", "--mode", "chunks", "--backend", "warp-cpu", "--n", n(200_000)],
                {},
            ),
            (
                "06-t12-statistical-parity",
                [*st, "t12", "--pairs", "python:cpu32,cpu32:cpu64", "--scale", str(scale), *w],
                {},
            ),
            ("07-t8-mcs-step-independence", [*st, "t8", "--n", n(1_000_000), *w], {}),
            ("08-t9-idd-step-independence", [*st, "t9", "--n", n(1_000_000), *w], {}),
            ("09-t10-rotation-invariance", [*st, "t10", "--n", n(1_000_000), *w], {}),
            ("10-t14-voxel-boundary-bias", [*st, "t14", "--n", n(1_000_000), *w], {}),
        ]
    if suite == "hr":
        return [
            (
                "01-pytest-cuda",
                pytest_cmd("tests/ionmc/test_transport_cuda.py", marker="cuda"),
                {"IONMC_REQUIRE_CUDA": "1"},
            ),
            (
                "02-t2-deterministic-csda-cuda",
                [*st, "t2", "--backend", "warp-cuda"],
                {"IONMC_REQUIRE_CUDA": "1"},
            ),
            (
                "03-t13-chunks-cuda",
                [*st, "t13", "--mode", "chunks", "--backend", "warp-cuda", "--n", n(1_000_000)],
                {"IONMC_REQUIRE_CUDA": "1"},
            ),
            (
                "04-t12-statistical-parity",
                [
                    *st,
                    "t12",
                    "--pairs",
                    "python:cpu32,cpu32:cuda32,cpu32:cpu64",
                    "--scale",
                    str(scale),
                    *w,
                ],
                {"IONMC_REQUIRE_CUDA": "1"},
            ),
        ]
    raise SystemExit(f"unknown suite {suite!r}")


def git(*args: str) -> str | None:
    try:
        r = subprocess.run(
            ["git", "-C", str(REPO), *args], capture_output=True, text=True, timeout=60
        )  # noqa: S603,S607
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout.strip() if r.returncode == 0 else None


def resolve_sha(expected: str) -> tuple[str, str]:
    if not re.fullmatch(r"[0-9a-f]{40}", expected):
        raise SystemExit(f"--expected-sha must be 40 lowercase hex digits, got {expected!r}")
    head = git("rev-parse", "HEAD")
    if head is None:
        return expected, "declared (no git repository readable)"
    if head != expected:
        raise SystemExit(f"--expected-sha {expected} differs from git rev-parse HEAD {head}")
    return head, "git"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def environment_text(sha: str, source: str, args: argparse.Namespace, workers: int) -> str:
    import numpy
    import warp

    cpu = "unknown"
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    porcelain = git("status", "--porcelain")
    dirty = "unknown" if porcelain is None else ("yes" if porcelain else "no")
    lines = [
        f"date_utc={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
        f"suite={args.suite}",
        f"python={platform.python_version()}",
        f"warp={warp.__version__}",
        f"numpy={numpy.__version__}",
        f"cpu={cpu}",
        f"logical_cpus={os.cpu_count()}",
        f"kernel={platform.release()}",
        f"git_sha={sha}",
        f"sha_source={source}",
        f"tree_dirty={dirty}",
        f"workers={workers}",
        f"step_timeout_s={args.step_timeout}",
        f"scale={args.scale}",
        "script_sha256:",
    ]
    files = sorted(HERE.glob("*.py")) + sorted((REPO / "src/ionmc/transport").glob("*.py"))
    files += [REPO / "src/ionmc/config.py", REPO / "src/ionmc/simulation.py"]
    lines += [f"  {sha256(f)}  {f.relative_to(REPO)}" for f in files]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--suite", choices=("lv", "hr"), required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--expected-sha", required=True)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--step-timeout", type=int, default=1500)
    ap.add_argument(
        "--only",
        nargs="+",
        metavar="STEP",
        help="run only these steps (full name or two-digit prefix) into this output directory; "
        "the manifest lists just them (split a long suite over several output directories)",
    )
    ap.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="factor on all history counts (< 1 gives a non-conformant, labelled run)",
    )
    args = ap.parse_args(argv)
    if args.workers < 2 or args.step_timeout < 1 or not 0.0 < args.scale <= 1.0:
        raise SystemExit("need --workers >= 2, --step-timeout >= 1 and 0 < --scale <= 1")
    out = Path(args.out).resolve()
    if not any(out.is_relative_to(p) for p in ALLOWED_PARENTS):
        raise SystemExit(
            f"--out must lie under validation/generated/ or benchmarks/generated/: {out}"
        )
    if out.exists():
        raise SystemExit(f"refusing to reuse an existing results directory: {out}")
    sha, source = resolve_sha(args.expected_sha)
    steps = suite_steps(args.suite, args.workers, args.scale)
    if args.only:
        keep = [x for x in steps if any(x[0] == o or x[0].startswith(f"{o}-") for o in args.only)]
        if len(keep) != len(set(args.only)):
            raise SystemExit(
                f"--only {args.only} does not select exactly those steps of {args.suite}"
            )
        steps = keep
    out.mkdir(parents=True)
    (out / "environment.txt").write_text(environment_text(sha, source, args, args.workers))
    (out / "manifest.txt").write_text("".join(f"{name}\n" for name, _, _ in steps))
    env_base = dict(os.environ)
    env_base.setdefault("WARP_CACHE_PATH", str(out / "warp-cache"))
    env_base["PYTHONDONTWRITEBYTECODE"] = "1"
    env_base["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env_base.get("PYTHONPATH", "")
    failures = 0
    for name, cmd, extra in steps:
        path = out / f"{name}.txt"
        started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with path.open("w") as fh:
            fh.write(f"# command: {' '.join(cmd)}\n# git_sha: {sha}\n# started_utc: {started}\n")
            fh.write(f"# step_timeout_s: {args.step_timeout}\n")
            fh.flush()
            try:
                r = subprocess.run(
                    cmd,
                    cwd=REPO,
                    env={**env_base, **extra},
                    stdout=fh,
                    stderr=subprocess.STDOUT,
                    timeout=args.step_timeout,
                )  # noqa: S603
                code = r.returncode
            except subprocess.TimeoutExpired:
                code = 124
            fh.write(f"\n# exit={code}\n")
        print(f"== {name}: exit={code}", flush=True)
        failures += code != 0
    sys.path.insert(0, str(HERE))
    import summarize

    rc = summarize.main([str(out), "--expected-sha", sha])
    print(f"done: {out} (failed_steps={failures})")
    return 1 if failures or rc else 0


if __name__ == "__main__":
    sys.exit(main())
