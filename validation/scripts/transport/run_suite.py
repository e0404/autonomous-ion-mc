"""Fail-closed runner of the V3-003 (``lv``, ``hr``) and V3-004 (``lv4``, ``hr4``) local-validation
(LV) and host-runner (HR) suites. The V3-004 suites use ``steps_v4.py``, the default seed base
20381004 and, in the hashed set, the V3-004 acceptance plan and the synthetic lookup fixture.

Usage (argv only, no shell; the host runner executes exactly this)::

    python validation/scripts/transport/run_suite.py --suite {lv,hr,lv4,hr4} \
        --out validation/generated/transport/<new-dir> --expected-sha <40 hex> \
        [--workers N|auto] [--step-timeout SECONDS] [--scale F] [--python-parts N] \
        [--only STEP ...] [--import-dirs DIR ...] [--seed-base INT]

The qualification command needs no seed flag: the default ``--seed-base`` is the qualification base
20341004 (a run with the rehearsal base 20261004 or the consumed bases 20271004, 20281004,
20291004, 20301004, 20311004, 20321004 and 20331004 is archived but never conformant)::

    python validation/scripts/transport/run_suite.py --suite hr --expected-sha <sha> --out <new-dir>

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
* ``--workers 1`` (or ``--single-process``) is the single-process diagnostic mode: every step runs
  with one worker and single-threaded numerics (``SINGLE_PROCESS_ENV``), the steps whose purpose is
  multiprocessing (``DEFERRED_STEPS``) are archived as ``deferred`` and not run, histories, seeds
  and criteria are unchanged and the archive is never ``conformant``;
* any failure exits non-zero.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import re
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ALLOWED_PARENTS = (REPO / "validation" / "generated", REPO / "benchmarks" / "generated")
STEPS = HERE / "steps.py"
STEPS_V4 = HERE / "steps_v4.py"
SUITES = ("lv", "hr", "lv4", "hr4")
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
SOURCE_FILES_V4 = (
    *SOURCE_FILES,
    "validation/plans/v3-004-acceptance.md",
    "tests/data/synthetic_lookup.json",
)
"""The hashed set of the suites ``lv4`` and ``hr4``: the V3-003 set plus the V3-004 acceptance plan
and the synthetic lookup fixture (outside ``tests/ionmc``)."""


def source_file_list(suite: str | None = None) -> tuple[str, ...]:
    """Individually hashed files of ``suite`` (default: the suite named by ``IONMC_RUN_SUITE``)."""
    suite = suite or os.environ.get("IONMC_RUN_SUITE", "")
    return SOURCE_FILES_V4 if suite in ("lv4", "hr4") else SOURCE_FILES

DEFAULT_PYTHON_PARTS = 2
QUALIFICATION_SEED_BASE = 20341004
V4_QUALIFICATION_SEED_BASE = 20381004
V4_CONSUMED_SEED_BASES = (20361004, 20371004)  # plan amendments 4 and 5
V4_REHEARSAL_SEED_BASE = 20351004
DEFAULT_SEED_BASES = {
    "lv": QUALIFICATION_SEED_BASE,
    "hr": QUALIFICATION_SEED_BASE,
    "lv4": V4_QUALIFICATION_SEED_BASE,
    "hr4": V4_QUALIFICATION_SEED_BASE,
}
REHEARSAL_SEED_BASE = 20261004
CONSUMED_SEED_BASE = 20271004  # used by the T9 investigation: not a qualification base
DEFAULT_SEED_BASE = QUALIFICATION_SEED_BASE
EXECUTION_STANDARD = "standard"
EXECUTION_SINGLE_PROCESS = "single-process-diagnostic"
SINGLE_PROCESS_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "IONMC_SINGLE_PROCESS": "1",
}
"""Environment of every step in the single-process diagnostic mode (``--workers 1``)."""
DEFERRED_REASON = (
    "multiprocessing-specific check deferred in single-process diagnostic mode "
    "(operator directive 2026-10-07)"
)
DEFERRED_STEPS = {"lv": ("t13-workers",), "hr": (), "lv4": ("a15-workers",), "hr4": ()}
"""Steps whose purpose is the worker-partition invariance (1 versus N worker processes): they are
not executed in the single-process diagnostic mode and recorded as ``deferred``. Steps that merely
use workers for speed run with one worker, with unchanged histories, seeds and criteria."""


STEP_TIMEOUT_FLOOR_S = {
    "lv4": {"a9-part-1of2": 3300, "a9-part-2of2": 3300, "a7-step-independence": 1800},
    "hr4": {"a11-hr-channel-parity": 3600},
}
"""Minimum step timeout [s] of long steps (the ``--step-timeout`` default is 1500 s); the effective
timeout is the larger of the two and is recorded in the step header."""


def step_timeout_s(suite: str, name: str, default: int) -> int:
    floor = STEP_TIMEOUT_FLOOR_S.get(suite, {}).get(name.split("-", 1)[1], 0)
    return max(default, floor)


def deferred_step_names(suite: str, python_parts: int = DEFAULT_PYTHON_PARTS) -> list[str]:
    """Full names of the steps deferred in the single-process diagnostic mode."""
    return [
        n
        for n in full_step_names(suite, python_parts)
        if n.split("-", 1)[1] in DEFERRED_STEPS[suite]
    ]


def pytest_cmd(*targets: str, marker: str | None = None) -> list[str]:
    """``marker`` is a pytest ``-m`` expression."""
    cmd = [PY, "-m", "pytest", "-q", "-p", "no:cacheprovider", "--no-header", *targets]
    return cmd + (["-m", marker] if marker else [])


def suite_steps(
    suite: str,
    workers: int,
    scale: float,
    python_parts: int = DEFAULT_PYTHON_PARTS,
    out: Path | None = None,
    import_dirs: list[str] | None = None,
    step_timeout: int = 1500,
    seed_base: int | None = None,
    single_process: bool = False,
) -> list[tuple[str, list[str], dict[str, str]]]:
    """``(name, argv, extra environment)`` of the steps of a suite, in execution order (the
    names depend only on ``suite`` and ``python_parts``). ``single_process`` deselects the
    ``multiprocess`` tests of the pytest steps (the deferred steps are chosen by the caller)."""

    if seed_base is None:
        seed_base = DEFAULT_SEED_BASES.get(suite, DEFAULT_SEED_BASE)

    def n(x: int) -> str:
        return str(max(1000, int(x * scale)))

    out_dir = (out or Path("OUT")) / "samples"
    dirs = [str(out or Path("OUT")), *(import_dirs or [])]
    st = [PY, str(STEPS)]
    w = ["--workers", str(workers)]
    sc = ["--scale", str(scale)]
    cuda = {"IONMC_REQUIRE_CUDA": "1"} if suite in ("hr", "hr4") else {}
    s4 = [PY, str(STEPS_V4)]
    steps: list[tuple[str, list[str], dict[str, str]]] = []

    def add(name: str, cmd: list[str], env: dict[str, str] | None = None) -> None:
        full = f"{len(steps) + 1:02d}-{name}"
        if cmd[:2] in ([PY, str(STEPS)], [PY, str(STEPS_V4)]):
            eff = step_timeout_s(suite, full, step_timeout)
            inner = ["--timeout", str(max(30, int(0.9 * eff)))]  # the pool cleans up first
            cmd = [*cmd, *inner, "--seed-base", str(seed_base)]  # every seed derives from the base
        steps.append((full, cmd, env or {}))

    if suite in ("lv4", "hr4"):
        _suite_steps_v4(suite, add, s4, w, sc, cuda, out_dir, dirs, single_process)
        return steps

    if suite == "lv":
        add(
            "pytest-warp-cpu-t2-t4-t11-c1-t13",
            pytest_cmd(
                "tests/ionmc/test_transport_warp.py",
                "tests/ionmc/test_transport_partition.py",
                "tests/ionmc/test_config_validation.py",
                marker="not multiprocess" if single_process else None,
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
        acc, pairs = "cpu32,cpu64", "python:cpu32,python:cpu64,cpu32:cpu64"
    elif suite == "hr":
        add(
            "pytest-cuda",
            pytest_cmd(
                "tests/ionmc/test_transport_cuda.py",
                marker="cuda and not multiprocess" if single_process else "cuda",
            ),
            cuda,
        )
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
        acc, pairs = "cpu32,cpu64,cuda32", "python:cpu32,python:cpu64,cpu32:cuda32,cpu32:cpu64"
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


def _suite_steps_v4(suite, add, s4, w, sc, cuda, out_dir, dirs, single_process):  # type: ignore[no-untyped-def]
    """Steps of the V3-004 suites (rows A7, A8, A9, A11, A13, A15, A16 of
    ``validation/plans/v3-004-acceptance.md``); every statistical step takes ``--scale``."""
    if suite == "lv4":
        add(
            "pytest-scoring-warp-cpu",
            pytest_cmd(
                "tests/ionmc/test_let_offline.py",
                "tests/ionmc/test_scoring_channels.py",
                "tests/ionmc/test_scoring_transport.py",
                "tests/ionmc/test_scoring_warp.py",
                marker="not multiprocess" if single_process else None,
            ),
        )
        add("a16-qualified-path-regression", [*s4, "a16"])
        add("a11-lv-python-vs-warp-cpu-256x150MeV", [*s4, "a11-lv", *sc])
        add("a15-chunks-cpu", [*s4, "a15", "--mode", "chunks", "--backend", "warp-cpu", *sc])
        add("a15-workers", [*s4, "a15", "--mode", "workers", "--workers", "3", *sc])
        add("a7-step-independence", [*s4, "a7", *sc, *w])
        add("a8-offline-let", [*s4, "a8", *sc, *w])
        add("a13-let-profile-exploratory", [*s4, "a13", *sc, *w])
        for i in (1, 2):
            add(
                f"a9-part-{i}of2",
                [*s4, "a9-part", "--part", f"{i}/2", "--out-dir", str(out_dir), *sc, *w],
            )
        add("a9-compare", [*s4, "a9-compare", "--dirs", *dirs, *sc])
    else:
        add(
            "pytest-cuda-scoring",
            pytest_cmd(
                "tests/ionmc/test_scoring_warp.py",
                marker="cuda and not multiprocess" if single_process else "cuda",
            ),
            cuda,
        )
        add("a15-chunks-cuda", [*s4, "a15", "--mode", "chunks", "--backend", "warp-cuda", *sc],
            cuda)  # fmt: skip
        add("a11-hr-channel-parity", [*s4, "a11-hr", *sc], cuda)


KILL_GRACE_S = 10.0


def run_step(
    cmd: list[str],
    *,
    cwd: Path,
    env: dict[str, str],
    stdout: Any,
    timeout: float,
    grace: float = KILL_GRACE_S,
) -> int:
    """Run one step in its own process group and return its exit code (124 on timeout).

    The step may spawn worker processes (the multiprocessing pool of the transport engine); on a
    timeout, and in any case when the step has ended, the whole group is terminated (SIGTERM, then
    SIGKILL after ``grace`` seconds) so that no descendant outlives its step."""
    proc = subprocess.Popen(
        cmd, cwd=cwd, env=env, stdout=stdout, stderr=subprocess.STDOUT, start_new_session=True
    )  # noqa: S603
    pgid = proc.pid  # a new session: the process is the leader of its own group
    code = 124
    try:
        code = proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        pass
    finally:
        _kill_group(pgid, grace)
        proc.wait()
    return code


def _kill_group(pgid: int, grace: float) -> None:
    """Terminate every process of the group ``pgid`` (no error if none is left)."""
    try:
        os.killpg(pgid, signal.SIGTERM)
    except ProcessLookupError:
        return
    deadline = time.monotonic() + grace
    while time.monotonic() < deadline:
        try:
            os.killpg(pgid, 0)
        except ProcessLookupError:
            return
        time.sleep(0.1)
    try:
        os.killpg(pgid, signal.SIGKILL)
    except ProcessLookupError:
        return
    time.sleep(0.2)


def full_step_names(suite: str, python_parts: int) -> list[str]:
    """Names of the complete suite (used to decide ``subset`` and to combine subsets)."""
    return [s[0] for s in suite_steps(suite, 2, 1.0, python_parts)]  # names only


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


def source_files(suite: str | None = None) -> list[Path]:
    """Every source file whose hash identifies the code under test (the set of ``suite``, default
    the suite of ``IONMC_RUN_SUITE``)."""
    files: list[Path] = [REPO / f for f in source_file_list(suite) if (REPO / f).is_file()]
    for prefix in SOURCE_PREFIXES:
        for p in sorted((REPO / prefix).rglob("*")):
            parts = set(p.relative_to(REPO).parts)
            if p.is_file() and not parts & {"__pycache__", "generated"} and p.suffix != ".pyc":
                files.append(p)
    return files


def environment_text(
    sha: str, source: str, dirty: str, args: argparse.Namespace, workers: int, single: bool = False
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
        f"execution_mode={EXECUTION_SINGLE_PROCESS if single else EXECUTION_STANDARD}",
        "single_process_env="
        + (",".join(f"{k}={v}" for k, v in SINGLE_PROCESS_ENV.items()) if single else ""),
        f"step_timeout_s={args.step_timeout}",
        f"scale={args.scale}",
        f"seed_base={args.seed_base}",
        f"python_parts={args.python_parts}",
        f"only={','.join(args.only) if args.only else ''}",
        "source_hashes:",
    ]
    lines += [f"  {sha256(f)}  {f.relative_to(REPO)}" for f in source_files(args.suite)]
    return "\n".join(lines) + "\n"


def resolve_workers(value: str) -> int:
    """``auto`` is every logical CPU of this host; otherwise an integer."""
    return (os.cpu_count() or 1) if value == "auto" else int(value)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--suite", choices=SUITES, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--expected-sha", required=True)
    ap.add_argument(
        "--workers",
        default="auto",
        help="integer >= 1 or 'auto' (all cores); 1 is the single-process diagnostic mode",
    )
    ap.add_argument(
        "--single-process",
        action="store_true",
        help="single-process diagnostic mode (implied by --workers 1): one worker per sample, "
        "single-threaded numerics, the multiprocessing-specific steps are deferred (recorded, not "
        "run) and the archive is never conformant; histories, seeds and criteria are unchanged",
    )
    ap.add_argument("--step-timeout", type=int, default=1500)
    ap.add_argument("--python-parts", type=int, default=DEFAULT_PYTHON_PARTS)
    ap.add_argument(
        "--seed-base",
        type=int,
        default=None,
        help="base of all statistical seeds (default: the qualification base of the suite, "
        "20341004 for lv/hr and 20381004 for lv4/hr4 (20351004 rehearsal and 20361004, 20371004 consumed give non-conformant archives); for lv/hr the "
        "bases 20261004, 20271004, 20281004, 20291004, 20301004, 20311004, 20321004 and "
        "20331004 give non-conformant archives); "
        "recorded in the archive",
    )
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
    if args.seed_base is None:
        args.seed_base = DEFAULT_SEED_BASES[args.suite]
    workers = 1 if args.single_process else resolve_workers(args.workers)
    if workers < 1 or args.step_timeout < 1 or not 0.0 < args.scale <= 1.0:
        raise SystemExit("need --workers >= 1, --step-timeout >= 1 and 0 < --scale <= 1")
    single = workers == 1
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
    steps = suite_steps(
        args.suite,
        workers,
        args.scale,
        args.python_parts,
        out,
        args.import_dirs,
        args.step_timeout,
        args.seed_base,
        single,
    )
    deferred = set(deferred_step_names(args.suite, args.python_parts)) if single else set()
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
    (out / "environment.txt").write_text(
        environment_text(sha, source, dirty, args, workers, single)
    )
    (out / "manifest.txt").write_text("".join(f"{name}\n" for name, _, _ in steps))
    env_base = dict(os.environ)
    env_base.setdefault("WARP_CACHE_PATH", str(out / "warp-cache"))
    env_base["PYTHONDONTWRITEBYTECODE"] = "1"
    env_base["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env_base.get("PYTHONPATH", "")
    env_base["IONMC_RUN_SHA"] = sha
    env_base["IONMC_RUN_SUITE"] = args.suite
    if single:
        env_base.update(SINGLE_PROCESS_ENV)
    else:
        env_base.pop("IONMC_SINGLE_PROCESS", None)  # a standard run is never contaminated
    failures = 0
    for name, cmd, extra in steps:
        path = out / f"{name}.txt"
        started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with path.open("w") as fh:
            fh.write(f"# command: {' '.join(cmd)}\n# git_sha: {sha}\n# started_utc: {started}\n")
            eff_timeout = step_timeout_s(args.suite, name, args.step_timeout)
            fh.write(f"# step_timeout_s: {eff_timeout}\n")
            if name in deferred:
                fh.write(f"# status: deferred\n# reason: {DEFERRED_REASON}\n\n# exit=0\n")
                print(f"== {name}: deferred", flush=True)
                continue
            fh.flush()
            code = run_step(
                cmd,
                cwd=REPO,
                env={**env_base, **extra},
                stdout=fh,
                timeout=eff_timeout,
            )
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
