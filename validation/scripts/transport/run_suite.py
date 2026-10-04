"""Fail-closed runner of the V3-003 local-validation (LV) and host-runner (HR) suites.

Usage (argv only, no shell; the host runner executes exactly this)::

    python validation/scripts/transport/run_suite.py --suite {lv,hr} \
        --out validation/generated/transport/<new-dir> --expected-sha <40 hex> \
        [--workers N|auto] [--step-timeout SECONDS] [--scale F] [--python-parts N] \
        [--only STEP ...] [--import-dirs DIR ...]

Rules (the style of ``validation/scripts/warp-architecture/run_all.sh``):

* the output directory must not exist and must lie under ``validation/generated/`` or
  ``benchmarks/generated/`` of this repository (both are git-ignored); it is never reused;
* ``--expected-sha`` is mandatory (40 hex digits). If ``git`` can read this repository it must
  equal ``git rev-parse HEAD`` **and the working tree must be clean** (a dirty tree is refused
  before any step); if the tree is a snapshot without ``.git`` the declared SHA is recorded as
  ``sha_source=declared`` and later attested with ``summarize.py --attest-sha``;
* ``environment.txt`` records versions, hardware, the SHA, the dirty state and ``source_hashes``
  (sha256 of every file under ``src/ionmc``, ``tests/ionmc``, ``validation/scripts/transport`` and
  ``benchmarks/transport`` and of ``pyproject.toml``, ``uv.lock`` and the acceptance plan), in both
  the git and the snapshot case;
* every step archives stdout+stderr in ``NN-name.txt`` between a header (command, SHA, start
  time, timeout) and an ``# exit=`` trailer; a timed-out step is killed and archived with
  ``exit=124``; ``manifest.txt`` lists the step names of this run;
* ``summarize.py`` verifies the manifest, the headers and the SHA of every file and writes
  ``summary.json`` with one verdict per step; a run restricted with ``--only`` is a ``subset``
  and never ``conformant``; ``summarize.py --combine`` joins the subsets of one suite;
* the T12 statistics are split into steps with their own outputs: the python sample (in
  ``--python-parts`` history ranges), the accelerated samples and the comparison, which loads
  the saved samples (hash-verified; ``--import-dirs`` names archives of other runs holding them);
* any failure exits non-zero.
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
SOURCE_PREFIXES = (
    "src/ionmc",
    "tests/ionmc",
    "validation/scripts/transport",
    "benchmarks/transport",
)
SOURCE_FILES = ("pyproject.toml", "uv.lock", "validation/plans/v3-003-acceptance.md")
"""Every tracked file that defines what is executed and judged (code, tests and fixtures, the
project definition, the lock file and the frozen acceptance plan)."""
DEFAULT_PYTHON_PARTS = 2


def pytest_cmd(*targets: str, marker: str | None = None) -> list[str]:
    cmd = [PY, "-m", "pytest", "-q", "-p", "no:cacheprovider", "--no-header", *targets]
    return cmd + (["-m", marker] if marker else [])


def suite_steps(
    suite: str,
    workers: int,
    scale: float,
    python_parts: int = DEFAULT_PYTHON_PARTS,
    out: Path | None = None,
    import_dirs: list[str] | None = None,
) -> list[tuple[str, list[str], dict[str, str]]]:
    """``(name, argv, extra environment)`` of the steps of a suite, in execution order (the
    names depend only on ``suite`` and ``python_parts``)."""

    def n(x: int) -> str:
        return str(max(1000, int(x * scale)))

    out_dir = (out or Path("OUT")) / "samples"
    dirs = [str(out or Path("OUT")), *(import_dirs or [])]
    st = [PY, str(STEPS)]
    w = ["--workers", str(workers)]
    sc = ["--scale", str(scale)]
    cuda = {"IONMC_REQUIRE_CUDA": "1"} if suite == "hr" else {}
    steps: list[tuple[str, list[str], dict[str, str]]] = []

    def add(name: str, cmd: list[str], env: dict[str, str] | None = None) -> None:
        steps.append((f"{len(steps) + 1:02d}-{name}", cmd, env or {}))

    if suite == "lv":
        add(
            "pytest-warp-cpu-t2-t4-t11-c1-t13",
            pytest_cmd(
                "tests/ionmc/test_transport_warp.py",
                "tests/ionmc/test_transport_partition.py",
                "tests/ionmc/test_config_validation.py",
            ),
        )
        add(
            "t1-trace-parity-256x150MeV",
            [*st, "t1", "--k", str(max(4, int(256 * scale))), "--energy", "150", *w],
        )
        add("t2-deterministic-csda", [*st, "t2", "--backend", "warp-cpu"])
        add("t13-workers", [*st, "t13", "--mode", "workers", "--n", n(200_000), *w])
        add(
            "t13-chunks-cpu",
            [*st, "t13", "--mode", "chunks", "--backend", "warp-cpu", "--n", n(200_000)],
        )
        add("t-r1-python-repeatability", [*st, "t-r1", "--runs", "python:float64:400", *w])
        acc, pairs = "cpu32,cpu64", "python:cpu32,cpu32:cpu64"
    elif suite == "hr":
        add("pytest-cuda", pytest_cmd("tests/ionmc/test_transport_cuda.py", marker="cuda"), cuda)
        add("t2-deterministic-csda-cuda", [*st, "t2", "--backend", "warp-cuda"], cuda)
        add(
            "t13-chunks-cuda",
            [*st, "t13", "--mode", "chunks", "--backend", "warp-cuda", "--n", n(1_000_000)],
            cuda,
        )
        add(
            "t-r1-warp-repeatability",
            [*st, "t-r1", "--runs", "warp-cpu:float64:100000,warp-cuda:float32:100000", *w],
            cuda,
        )
        acc, pairs = "cpu32,cpu64,cuda32", "python:cpu32,cpu32:cuda32,cpu32:cpu64"
    else:
        raise SystemExit(f"unknown suite {suite!r}")
    for i in range(1, python_parts + 1):
        add(
            f"t12-python-sample-{i}of{python_parts}",
            [
                *st,
                "t12-python-sample",
                "--part",
                f"{i}/{python_parts}",
                "--out-dir",
                str(out_dir),
                *sc,
                *w,
            ],
            cuda,
        )
    add(
        "t12-accelerated-samples",
        [*st, "t12-accelerated-samples", "--samples", acc, "--out-dir", str(out_dir), *sc, *w],
        cuda,
    )
    add("t12-compare", [*st, "t12-compare", "--pairs", pairs, "--dirs", *dirs, *sc], cuda)
    if suite == "lv":
        add("t8-mcs-step-independence", [*st, "t8", "--n", n(1_000_000), *w])
        add("t9-idd-step-independence", [*st, "t9", "--n", n(1_000_000), *w])
        add("t10-rotation-invariance", [*st, "t10", "--n", n(1_000_000), *w])
        add("t14-voxel-boundary-bias", [*st, "t14", "--n", n(1_000_000), *w])
    return steps


def full_step_names(suite: str, python_parts: int) -> list[str]:
    """Names of the complete suite (used to decide ``subset`` and to combine subsets)."""
    return [s[0] for s in suite_steps(suite, 2, 1.0, python_parts)]


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


def source_files() -> list[Path]:
    """Every source file whose hash identifies the code under test."""
    files: list[Path] = [REPO / f for f in SOURCE_FILES if (REPO / f).is_file()]
    for prefix in SOURCE_PREFIXES:
        for p in sorted((REPO / prefix).rglob("*")):
            parts = set(p.relative_to(REPO).parts)
            if p.is_file() and not parts & {"__pycache__", "generated"} and p.suffix != ".pyc":
                files.append(p)
    return files


def environment_text(
    sha: str, source: str, dirty: str, args: argparse.Namespace, workers: int
) -> str:
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
        f"python_parts={args.python_parts}",
        f"only={','.join(args.only) if args.only else ''}",
        "source_hashes:",
    ]
    lines += [f"  {sha256(f)}  {f.relative_to(REPO)}" for f in source_files()]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--suite", choices=("lv", "hr"), required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--expected-sha", required=True)
    ap.add_argument("--workers", default="auto", help="integer >= 2 or 'auto' (all cores)")
    ap.add_argument("--step-timeout", type=int, default=1500)
    ap.add_argument("--python-parts", type=int, default=DEFAULT_PYTHON_PARTS)
    ap.add_argument(
        "--only",
        nargs="+",
        metavar="STEP",
        help="run only these steps (full name or two-digit prefix) into this output directory; "
        "the manifest lists just them and the summary is a non-conformant subset",
    )
    ap.add_argument(
        "--import-dirs",
        nargs="*",
        default=[],
        help="archives of other runs of this suite and SHA that hold T12 sample files",
    )
    ap.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="factor on all history counts (< 1 gives a non-conformant, labelled run)",
    )
    args = ap.parse_args(argv)
    workers = (os.cpu_count() or 1) if args.workers == "auto" else int(args.workers)
    if workers < 2 or args.step_timeout < 1 or not 0.0 < args.scale <= 1.0:
        raise SystemExit("need --workers >= 2, --step-timeout >= 1 and 0 < --scale <= 1")
    if not 1 <= args.python_parts <= 64:
        raise SystemExit("--python-parts must be in [1, 64]")
    out = Path(args.out).resolve()
    if not any(out.is_relative_to(p) for p in ALLOWED_PARENTS):
        raise SystemExit(
            f"--out must lie under validation/generated/ or benchmarks/generated/: {out}"
        )
    if out.exists():
        raise SystemExit(f"refusing to reuse an existing results directory: {out}")
    sha, source = resolve_sha(args.expected_sha)
    steps = suite_steps(args.suite, workers, args.scale, args.python_parts, out, args.import_dirs)
    if args.only:
        keep = [x for x in steps if any(x[0] == o or x[0].startswith(f"{o}-") for o in args.only)]
        if len(keep) != len(set(args.only)):
            raise SystemExit(
                f"--only {args.only} does not select exactly those steps of {args.suite}"
            )
        steps = keep
    porcelain = git("status", "--porcelain")
    if porcelain:
        raise SystemExit(
            "refusing to run on a dirty working tree (commit or stash first):\n" + porcelain[:2000]
        )
    dirty = "unknown" if porcelain is None else "no"
    out.mkdir(parents=True)
    (out / "environment.txt").write_text(environment_text(sha, source, dirty, args, workers))
    (out / "manifest.txt").write_text("".join(f"{name}\n" for name, _, _ in steps))
    env_base = dict(os.environ)
    env_base.setdefault("WARP_CACHE_PATH", str(out / "warp-cache"))
    env_base["PYTHONDONTWRITEBYTECODE"] = "1"
    env_base["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env_base.get("PYTHONPATH", "")
    env_base["IONMC_RUN_SHA"] = sha
    env_base["IONMC_RUN_SUITE"] = args.suite
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
