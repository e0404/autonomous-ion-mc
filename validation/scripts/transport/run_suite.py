"""Fail-closed runner of the V3-003 (``lv``, ``hr``) and V3-004 (``lv4``, ``hr4``) local-validation
(LV) and host-runner (HR) suites. The V3-004 suites use ``steps_v4.py``, the default seed base
20401004 and, in the hashed set, the V3-004 acceptance plan and the synthetic lookup fixture.

Usage (argv only, no shell; the host runner executes exactly this)::

    python validation/scripts/transport/run_suite.py --suite {lv,hr,lv4,hr4,lv5,lv5b,hr5} \
        --out validation/generated/transport/<new-dir> --expected-sha <40 hex> \
        [--workers N|auto] [--step-timeout SECONDS] [--scale F] [--python-parts N] \
        [--only STEP ...] [--import-dirs DIR ...] [--seed-base INT]

The qualification command needs no seed flag: the default ``--seed-base`` is the qualification base
20391004 (a run with the rehearsal base 20261004 or the consumed bases 20271004, 20281004,
20291004, 20301004, 20311004, 20321004, 20331004 and 20341004 is archived but never conformant;
V3-003D, plan amendment 24)::

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
  ``benchmarks/transport`` and of ``pyproject.toml``, ``uv.lock`` and the acceptance plan; the suites
  ``lv5b`` and ``hr5`` add the 96 files of the 24 frozen V5 reference cases, ``V5_CASE_FILES``), in
  both the git and the snapshot case;
* every step archives stdout+stderr in ``NN-name.txt`` between a header (command, SHA, start
  time, timeout) and an ``# exit=`` trailer; a timed-out step is killed and archived with
  ``exit=124``; ``manifest.txt`` lists the step names of this run;
* ``summarize.py`` verifies the manifest, the headers and the SHA of every file and writes
  ``summary.json`` with one verdict per step; a run restricted with ``--only`` is a ``subset``
  and never ``conformant``; ``summarize.py --combine`` joins the subsets of one suite;
* the T12 statistics are split into steps with their own outputs: the python sample (in
  ``--python-parts`` history ranges), the accelerated samples and the comparison, which loads
  the saved samples (hash-verified; ``--import-dirs`` names archives of other runs holding them; the lv5 shard
  partials carry a ``content_sha256`` and are bound to the run SHA, suite, table id and seed, but that digest
  is recomputable by whoever alters a file: every partial imported with ``--import-dirs`` must therefore be listed
  with its exact digest and a ``host_run_id`` in ``--partials-manifest`` (schema in
  ``steps_v5.read_manifest``). The orchestrator builds the manifest ONLY from the ``PARTIAL <name>
  <digest>`` stdout lines of the shard steps in the protected host-runner records and lists the run ids in
  record_local_validation; the code cannot verify those records. The manifest path and its sha256 are
  recorded in ``environment.txt`` (``partials_manifest``, ``partials_manifest_sha256``), the combine steps
  print an ``attestation`` block (manifest sha256, every name/digest/host_run_id used, the run SHA) that
  ``summarize.py --combine`` carries into the combined summary. These records document what was relied on;
  they do not make the archive conformant: any step that used an imported partial (origin other than the
  current output directory) is ``conformant: false`` with the reason ``imported_partials_unverified_by_code``
  (``pass`` unaffected). Under the single-process directive sharded rows are therefore combined from imported
  partials and are non-conformant by code; they become conformant only when all shards and the combine run in
  one invocation after the operator lifts the directive; the protected validation record lists the host run
  ids; without ``--import-dirs`` only partials of the current output directory are
  accepted);
* ``--workers 1`` (or ``--single-process``) is the single-process diagnostic mode: every step runs
  with one worker and single-threaded numerics (``SINGLE_PROCESS_ENV``), the steps whose purpose is
  multiprocessing (``DEFERRED_STEPS``) are archived as ``deferred`` and not run, histories, seeds
  and criteria are unchanged and the archive is never ``conformant``;
* any failure exits non-zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
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
STEPS_V5 = HERE / "steps_v5.py"
STEPS_V5B = HERE / "steps_v5b.py"
SUITES = ("lv", "hr", "lv4", "hr4", "lv5", "lv5b", "hr5")
PY = sys.executable
SOURCE_PREFIXES = (
    "src/ionmc",
    "tests/ionmc",
    "validation/scripts/transport",
    "benchmarks/transport",
)
SOURCE_FILES = (
    "pyproject.toml",
    "uv.lock",
    "validation/plans/v3-003-acceptance.md",
    "validation/plans/v3-003d-acceptance.md",
)
"""Every tracked file that defines what is executed and judged (code, tests and fixtures, the
project definition, the lock file and the frozen acceptance plan)."""
SOURCE_FILES_V4 = (
    *SOURCE_FILES,
    "validation/plans/v3-004-acceptance.md",
    "tests/data/synthetic_lookup.json",
)
"""The hashed set of the suites ``lv4`` and ``hr4``: the V3-003 set plus the V3-004 acceptance plan
and the synthetic lookup fixture (outside ``tests/ionmc``)."""


SOURCE_FILES_V5 = (
    *SOURCE_FILES_V4,
    "validation/plans/v3-005-acceptance.md",
    "decisions/0041-proton-nuclear-interactions.md",
    "validation/scripts/transport/steps_v5.py",
    "validation/scripts/transport/nuclear_checks.py",
)
"""The hashed set of the suite ``lv5`` (V3-005A, decision 0041): the V3-004 set plus the V3-005
acceptance plan, decision 0041 and, individually, the lv5 steps and the nuclear checks (the
prefix ``validation/scripts/transport`` and the nuclear package ``src/ionmc/nuclear`` are hashed in
every suite)."""


V5_CASE_ROOT = "validation/reference_cases"
V5_CASE_FILES = (
    *(
        f"{V5_CASE_ROOT}/topas/proton-water-{e}mev-idd-r20{tag}-seed{k}/{f}"
        for e in (150, 200)
        for tag in ("", "-emonly")
        for k in (1, 2, 3)
        for f in ("case.json", "input.txt")
    ),
    *(
        f"{V5_CASE_ROOT}/mcsquare/proton-water-{e}mev-idd-r20-{nuc}-seed{k}/{f}"
        for e in (150, 200)
        for nuc in ("on", "off")
        for k in (1, 2, 3)
        for f in ("case.json", "config.txt", "Plan.txt", "BDL_mono.txt", "CT.mhd", "CT.raw")
    ),
)
"""The 96 files of the 24 frozen V5 case directories (TOPAS ``proton-water-{150,200}mev-idd-r20
[-emonly]-seed{1,2,3}``, MCsquare ``proton-water-{150,200}mev-idd-r20-{on,off}-seed{1,2,3}``) that
``compare_idd_v5`` binds the V5 reference runs to byte for byte (C20). Enumerated, not a prefix: the
V3-010B cases and the exploratory X1-X4 cases are not V5 evidence. Every file must exist
(``source_files`` fails loudly otherwise)."""

SOURCE_FILES_V5B = (
    *SOURCE_FILES_V5,
    "validation/scripts/transport/steps_v5b.py",
    "validation/scripts/reference/compare_idd_v5.py",
    "validation/scripts/reference/compare_batches.py",
    *V5_CASE_FILES,
)
"""The hashed set of the suites ``lv5b`` and ``hr5`` (V3-005B): the lv5 set plus the slice-B steps
and the 96 files of the 24 frozen V5 reference cases (``V5_CASE_FILES``, so that the gitless-snapshot
attestation covers the inputs the V5 comparator binds to; the nuclear kernels, ``src/ionmc`` and the
tests are hashed through the prefixes)."""


def source_file_list(suite: str | None = None) -> tuple[str, ...]:
    """Individually hashed files of ``suite`` (default: the suite named by ``IONMC_RUN_SUITE``)."""
    suite = suite or os.environ.get("IONMC_RUN_SUITE", "")
    if suite in ("lv5b", "hr5"):
        return SOURCE_FILES_V5B
    if suite == "lv5":
        return SOURCE_FILES_V5
    return SOURCE_FILES_V4 if suite in ("lv4", "hr4") else SOURCE_FILES


DEFAULT_PYTHON_PARTS = 2
QUALIFICATION_SEED_BASE = 20391004  # V3-003D plan, amendment 24 of the V3-003 plan
CONSUMED_SEED_BASES = (
    20271004, 20281004, 20291004, 20301004, 20311004, 20321004, 20331004, 20341004,
)  # fmt: skip
V4_QUALIFICATION_SEED_BASE = 20401004  # amendment 6 of the V3-004 plan
V4_CONSUMED_SEED_BASES = (20361004, 20371004, 20381004)  # amendments 4 and 5, 20381004 by V3-003D
V4_REHEARSAL_SEED_BASE = 20351004
V5_QUALIFICATION_SEED_BASE = 20421004  # plan of V3-005, Seeds (rehearsal family 2043xxxx)
V5B_QUALIFICATION_SEED_BASE = 20471004  # amendment 11 of the V3-005 plan (lv5b); rehearsals 2046xxxx
V5B_CONSUMED_SEED_BASES = (20441004,)  # amendment 11: consumed by the failed V7 replicate step
# (host run RUN-20261008T132316Z-11c02312 at b40d8121)
HR5_QUALIFICATION_SEED_BASE = 20451004  # amendment 6 (hr5)
V3003D_REHEARSAL_FAMILY = "2041xxxx"  # rehearsals of V3-003D, never qualification evidence
DEFAULT_SEED_BASES = {
    "lv": QUALIFICATION_SEED_BASE,
    "hr": QUALIFICATION_SEED_BASE,
    "lv4": V4_QUALIFICATION_SEED_BASE,
    "hr4": V4_QUALIFICATION_SEED_BASE,
    "lv5": V5_QUALIFICATION_SEED_BASE,
    "lv5b": V5B_QUALIFICATION_SEED_BASE,
    "hr5": HR5_QUALIFICATION_SEED_BASE,
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
DEFERRED_STEPS = {
    "lv": ("t13-workers",),
    "hr": (),
    "lv4": ("a15-workers",),
    "hr4": (),
    "lv5": ("v3-workers-partition",),
    "lv5b": ("v3-workers-partition",),
    "hr5": (),
}
"""Steps whose purpose is the worker-partition invariance (1 versus N worker processes): they are
not executed in the single-process diagnostic mode and recorded as ``deferred``. Steps that merely
use workers for speed run with one worker, with unchanged histories, seeds and criteria."""


V7_REP_SHARDS = 8
"""Number of shard steps ``v7-rep-s{k}`` of the V7 replicate coverage (Amendment 13: eight simulations of 9e6
histories, 7200 replicates); defined here for the same reason as ``V2B_SHARDS``; ``steps_v5b.V7_SHARDS`` is
taken from it."""


STEP_TIMEOUT_FLOOR_S = {
    "lv4": {"a9-part-1of2": 3300, "a9-part-2of2": 3300, "a7-step-independence": 1800},
    "hr4": {"a11-hr-channel-parity": 3600},
    "lv5": {
        **{f"v2-{e}-s{k}": 3300 for e in (100, 150, 200) for k in range(3)},
        **{f"v2-probe-s05-s{k}": 3300 for k in range(2)},
        **{f"v2-probe-fe-s{k}": 3300 for k in range(3)},
        "v3-lv": 3300, "x1": 3300, "e1": 3300, "v4-v4b": 3300, "r1-a16-t1-regression": 3300,
        "n1-v1-v1b-d6": 3300, "v3-workers-partition": 3300,
    },
    "lv5b": {
        "lv5b-throughput": 3300, "v8-lv-python-vs-warp-cpu": 3300, "r1-nuc-regression": 3300,
        **{f"v5-{e}-{t}": 3300 for e in (150, 200) for t in ("on", "off")},
        "v5-compare": 1800,
        **{f"v2b-s{k}": 3300 for k in range(8)}, "v7-scan": 3300, "v7-shift": 3300,
        **{f"v7-rep-s{k}": 3300 for k in range(V7_REP_SHARDS)}, "v7-rep-ref": 3300,
        "pytest-v7-rep-calibration": 3300, "v3-workers-partition": 3300,
    },
    "hr5": {
        "v8-stat-python-s0": 3600, "v8-stat-python-s1": 3600, "v8-stat-cpu64": 3600,
        "v8-stat-cuda32": 3600, "v8-stat-cuda64": 3600, "v7-f32-f64-cuda": 3600,
    },
}
"""Minimum step timeout [s] of long steps (the ``--step-timeout`` default is 1500 s); the effective
timeout is the larger of the two and is recorded in the step header."""


V2B_SHARDS = 2
"""Number of shards of the V2b statistics of the V3-005B plan (Amendment 8). Defined here, not imported
from ``steps_v5b.py``: the host interpreter running this script has no ``ionmc`` (the steps get it via
``PYTHONPATH``), and ``steps_v5b`` imports it; ``steps_v5b.V2B_SHARDS`` is taken from this value."""


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
    partials_manifest: str | None = None,
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
    cuda = {"IONMC_REQUIRE_CUDA": "1"} if suite in ("hr", "hr4", "hr5") else {}
    s4 = [PY, str(STEPS_V4)]
    s5 = [PY, str(STEPS_V5)]
    s5b = [PY, str(STEPS_V5B)]
    steps: list[tuple[str, list[str], dict[str, str]]] = []

    def add(name: str, cmd: list[str], env: dict[str, str] | None = None) -> None:
        full = f"{len(steps) + 1:02d}-{name}"
        if cmd[:2] in ([PY, str(STEPS)], [PY, str(STEPS_V4)], [PY, str(STEPS_V5)], [PY, str(STEPS_V5B)]):
            eff = step_timeout_s(suite, full, step_timeout)
            inner = ["--timeout", str(max(30, int(0.9 * eff)))]  # the pool cleans up first
            cmd = [*cmd, *inner, "--seed-base", str(seed_base)]  # every seed derives from the base
        steps.append((full, cmd, env or {}))

    if suite == "lv5":
        if import_dirs and partials_manifest is None:
            raise SystemExit("lv5 --import-dirs requires --partials-manifest (name -> sha256)")
        pm = ["--partials-manifest", str(partials_manifest)] if import_dirs else []
        _suite_steps_v5(add, s5, sc, out_dir, [*dirs, *pm])
        return steps
    if suite in ("lv5b", "hr5"):
        if import_dirs and partials_manifest is None:
            raise SystemExit(f"{suite} --import-dirs requires --partials-manifest (name -> sha256)")
        pm = ["--partials-manifest", str(partials_manifest)] if import_dirs else []
        if suite == "lv5b":
            _suite_steps_v5b(add, s5, s5b, sc, out_dir, [*dirs, *pm])
        else:
            _suite_steps_hr5(add, s5b, sc, out_dir, [*dirs, *pm], cuda)
        return steps
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
                "tests/ionmc/test_range_quadrature.py",
                marker="not multiprocess" if single_process else None,
            ),
            NIST_REQUIRED_ENV,
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
        rec = A16_INTENDED_CHANGE
        a16 = [
            "--mode",
            "intended-change",
            "--intended-change-record",
            json.dumps(rec, sort_keys=True),
        ]
        add(
            "a16-qualified-path-regression",
            [*s4, "a16", *(a16 if rec else ["--mode", "regression"])],
        )
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


def _suite_steps_v5(add, s5, sc, out_dir, dirs):  # type: ignore[no-untyped-def]
    """Steps of the V3-005 slice-A suite ``lv5`` (rows of ``validation/plans/v3-005-acceptance.md``;
    ``steps_v5.py``); every step runs one process with ``nuclear=True`` on the python backend, the
    per-step history counts are the frozen ones (the shards of V2 and V2-probe follow from the
    measured throughput, see ``steps_v5.py``) and every statistical step takes ``--scale``."""
    env = NUCLEAR_ENV
    add("lv5-throughput", [*s5, "lv5-throughput"], env)
    add("n1-v1-v1b-d6", [*s5, "n1"], env)
    for e in (100, 150, 200):
        for k in range(3):
            add(f"v2-{e}-s{k}",
                [*s5, "v2-shard", "--energy", str(e), "--shard", str(k), "--out-dir", str(out_dir), *sc],
                env)  # fmt: skip
    add("v2-combine", [*s5, "v2-combine", "--dirs", *dirs, *sc], env)
    for name, shards in (("s05", 2), ("fe", 3)):
        for k in range(shards):
            add(f"v2-probe-{name}-s{k}",
                [*s5, "v2-probe-shard", "--probe", name, "--shard", str(k), "--out-dir", str(out_dir), *sc],
                env)  # fmt: skip
    add("v2-probe-combine", [*s5, "v2-probe-combine", "--dirs", *dirs, *sc], env)
    add("v3-lv", [*s5, "v3-lv", *sc], env)
    add("v4-v4b", [*s5, "v4-v4b", *sc], env)
    add("x1", [*s5, "x1", *sc], env)
    add("e1", [*s5, "e1", *sc], env)
    add("r1-a16-t1-regression", [*s5, "r1", *sc], env)
    add("v3-workers-partition", [*s5, "v3-workers"], env)  # deferred: runs only when the mode is lifted


def _suite_steps_v5b(add, s5, s5b, sc, out_dir, dirs):  # type: ignore[no-untyped-def]
    """Steps of the V3-005B slice-B suite ``lv5b`` (``steps_v5b.py``; rows V8-LV, R1 for the nuclear
    branch, V5 (ionmc side), V2b, V7 of ``validation/plans/v3-005-acceptance.md``): the python and
    warp-cpu float64 backends in one process, ``nuclear=True``. V6 and E1-B (r_index 12 and 13) belong
    to V3-005C (steps 12/13 reserved). The deferred step is the 1-vs-N worker partition of nuclear
    runs; any ``cpu_workers > 1`` is deferred with it."""
    shards = V2B_SHARDS
    env = NUCLEAR_ENV
    add("lv5b-throughput", [*s5b, "lv5b-throughput", *sc], env)
    add("v8-lv-python-vs-warp-cpu", [*s5b, "v8-lv", *sc], env)
    add("r1-nuc-regression", [*s5b, "r1-nuc", *sc], env)
    for e in (150, 200):
        for t in ("on", "off"):
            add(f"v5-{e}-{t}", [*s5b, "v5-ionmc", "--energy", str(e), "--nuclear", t,
                                "--out-dir", str(out_dir), *sc], env)  # fmt: skip
    # row V5 verdict: the four partials against the materialized TOPAS/MCsquare runs
    # (``/workspace/.ionmc-cache/reference-runs/REF-*`` of the host snapshot); consumes partials only
    add("v5-compare", [*s5b, "v5-compare", "--dirs", *dirs, *sc], env)
    for k in range(shards):
        add(f"v2b-s{k}", [*s5b, "v2b-shard", "--shard", str(k), "--out-dir", str(out_dir), *sc], env)
    add("v2b-combine", [*s5b, "v2b-combine", "--dirs", *dirs, *sc], env)
    add("v7-scan", [*s5b, "v7-scan", *sc], env)
    add("v7-shift", [*s5b, "v7-shift", *sc], env)
    for k in range(V7_REP_SHARDS):  # one simulation of 9e6 histories each (7200 replicates in total)
        add(f"v7-rep-s{k}", [*s5b, "v7-rep-shard", "--shard", str(k), "--out-dir", str(out_dir), *sc], env)
    add("v7-rep-ref", [*s5b, "v7-rep-ref", "--out-dir", str(out_dir), *sc], env)
    add(  # Monte Carlo calibration of the V7 rule (minutes); its pass is the pytest exit status
        "pytest-v7-rep-calibration",
        pytest_cmd("tests/ionmc/test_v7_coverage.py", marker="calibration"),
        env,
    )
    add("v7-rep-combine", [*s5b, "v7-rep-combine", "--dirs", *dirs, *sc], env)
    add("v3-workers-partition", [*s5, "v3-workers"], env)  # deferred: runs only when the mode is lifted


def _suite_steps_hr5(add, s5b, sc, out_dir, dirs, cuda):  # type: ignore[no-untyped-def]
    """Steps of the V3-005B host-runner suite ``hr5`` (rows V8 statistical parity and V7 f32/f64 on
    CUDA); every CUDA step runs in one controlling process with ``IONMC_REQUIRE_CUDA=1``. Step 03
    (V6 on CUDA) is reserved for V3-005C."""
    env = {**NUCLEAR_ENV, **cuda}
    for nm in ("python-s0", "python-s1", "cpu64", "cuda32", "cuda64"):
        add(f"v8-stat-{nm}", [*s5b, "v8-stat-sample", "--sample", nm, "--out-dir", str(out_dir), *sc], env)
    add("v8-stat-compare", [*s5b, "v8-stat-compare", "--dirs", *dirs, *sc], env)
    add("v7-f32-f64-cuda", [*s5b, "v7-f32", *sc], env)


NUCLEAR_ENV = {
    "IONMC_REQUIRE_DATA": "1",
    "IONMC_CACHE_DIR": str(REPO / ".ionmc-cache" / "ionmc-data"),
}
"""Environment of the ``lv5`` steps: the built nuclear table (decision 0041; ``steps_v5.TABLE_ID``)
and the data it derives from are read from the hash-verified cache staged by the orchestrator."""

NIST_REQUIRED_ENV = {
    "IONMC_REQUIRE_NIST": "1",
    "IONMC_CACHE_DIR": str(REPO / ".ionmc-cache" / "ionmc-data"),
}
"""Environment of the ``lv`` pytest step that runs the NIST-water case of V3-003D row D1: the
cached PSTAR table (data layer cache ``IONMC_CACHE_DIR`` with its ``manifests/`` and ``objects/``
layout, staged by the orchestrator into the git-ignored ``.ionmc-cache/ionmc-data`` of the runner
workspace from the experiment data cache) is required; a missing cache fails the test instead of
skipping it."""
A16_INTENDED_CHANGE: dict[str, Any] | None = None
"""The recorded exception to the fail-closed A16 regression of the ``lv4`` suite (plan amendment 6
of V3-004); ``None`` while there is none, which is the normal state: the step then runs
``--mode regression`` (the default), which is GATED: the no-tally digests of the tree under test
must equal those of ``A16_BASELINE`` (``steps_v4.py``) in every compared field, and tally
neutrality must hold. The V3-003D record (baseline a524f209, identity field
``range_construction``) was deleted by the first commit of V3-005A after the V3-003D merge, which
also advanced ``A16_BASELINE`` to the merge commit f3a1dd62; the record is preserved in the
delimited block of ``validation/plans/v3-003d-acceptance.md``.

A task that intentionally changes the qualified transport path re-introduces a record: it assigns a
dict literal on the line ``A16_INTENDED_CHANGE: dict[str, Any] | None = {`` (the line must start
exactly so, and the literal must end with a line ``}`` followed by a newline, because
:func:`a16_normalize` finds the record between ``A16_RECORD_START`` and ``A16_RECORD_END``) with the
keys of ``steps_v4.verify_intended_change`` (``task``, ``baseline``, ``identity_field``,
``baseline_value``, ``new_value``, ``allowed_differing_fields``, ``physical_limits``, ``bounds``,
``source_digest``, ``plan_block_sha256``), states the same record (without ``plan_block_sha256``)
in a delimited block of its own plan file (``A16_PLAN_FILE`` and ``steps_v4.PLAN_FILE`` are then
pointed at that plan), and refreshes ``source_digest`` with ``--print-a16-source-digest``. While
it is not ``None`` the ``lv4`` step runs ``--mode intended-change`` and passes the record to the
step, which (i) verifies that the plan block hashes to ``plan_block_sha256`` and states exactly
this record, (ii) verifies that the table identity ``identity_field`` of the baseline tree
(``baseline``, which must equal ``A16_BASELINE``) is ``baseline_value`` and that of the tree under
test is ``new_value``, and (iii) gates the comparison: only ``allowed_differing_fields`` may
differ, within ``bounds``; everything else must be bit-identical. It fails closed otherwise. Once
that task is merged the baseline tree carries ``new_value``, so the next task's first commit sets
this back to ``None`` and advances ``A16_BASELINE`` (the form above)."""

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


_PARTIAL_LINE = re.compile(r"^PARTIAL (\S+\.json) ([0-9a-f]{64})$")


def echo_partial_lines(archive: Path) -> None:
    """Forward every ``PARTIAL <name>.json <64 lowercase hex>`` line of the archived step output
    verbatim to this process's stdout (nothing else of the step output), so that the protected
    host-runner record, which keeps only this stdout, carries the shard partial digests.  Fail
    closed: a line starting with ``PARTIAL `` that does not match exactly, or a name repeated
    within one archive, raises SystemExit (unrelated child output must not pass as shard evidence)."""
    seen: set[str] = set()
    for line in archive.read_text(errors="replace").splitlines():
        if not line.startswith("PARTIAL "):
            continue
        match = _PARTIAL_LINE.match(line)
        if match is None:
            raise SystemExit(f"malformed PARTIAL line in {archive}: {line!r}")
        if match.group(1) in seen:
            raise SystemExit(f"duplicate PARTIAL name {match.group(1)!r} in {archive}: {line!r}")
        seen.add(match.group(1))
        print(line, flush=True)


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
    listed = source_file_list(suite)
    missing = [f for f in V5_CASE_FILES if f in listed and not (REPO / f).is_file()]
    if missing:
        raise SystemExit(
            f"{len(missing)} frozen V5 case file(s) missing from the source tree, first: "
            f"{missing[0]}; the hashed set cannot be completed")
    files: list[Path] = [REPO / f for f in listed if (REPO / f).is_file()]
    for prefix in SOURCE_PREFIXES:
        for p in sorted((REPO / prefix).rglob("*")):
            parts = set(p.relative_to(REPO).parts)
            if p.is_file() and not parts & {"__pycache__", "generated"} and p.suffix != ".pyc":
                files.append(p)
    return files


A16_PLAN_FILE = "validation/plans/v3-003d-acceptance.md"
A16_PLAN_BEGIN, A16_PLAN_END = (
    "<!-- A16-INTENDED-CHANGE-BEGIN -->",
    "<!-- A16-INTENDED-CHANGE-END -->",
)


A16_RUN_SUITE_FILE = "validation/scripts/transport/run_suite.py"
A16_SELF_REFERENCE = re.compile(r'("(?:source_digest|plan_block_sha256)": )"[^"]*"')
A16_RECORD_START, A16_RECORD_END = "\nA16_INTENDED_CHANGE: dict[", "\n}\n"
A16_EXPECTED_MASKS = {"record": 2, "plan": 1}
"""Number of self-referential literals masked per kind of file: the record in ``run_suite.py``
holds the source digest and the plan block hash, the plan block states the source digest only."""


def a16_normalize(text: str, kind: str) -> str:
    """``text`` with exactly the self-referential digest literals of the A16 record replaced by
    the placeholder ``MASKED``: ``kind`` ``"record"`` (``run_suite.py``: only inside the
    ``A16_INTENDED_CHANGE`` dict literal) or ``"plan"`` (only inside the delimited A16 block of the
    V3-003D plan). Every other byte, including look-alike text elsewhere, stays hashed. A record
    or block that is present must hold exactly the expected number of literals (fail closed). In the
    ``A16_INTENDED_CHANGE: dict[str, Any] | None = None`` form (no record) the ``run_suite.py`` text
    is returned unchanged."""
    if kind == "record":
        lo, hi = text.find(A16_RECORD_START), None
        if lo >= 0:
            eol = text.find("\n", lo + 1)
            if not text[lo + 1 : eol].rstrip().endswith("= {"):
                return text  # the ``= None`` form: no record, nothing to mask
            end = text.find(A16_RECORD_END, lo)
            if end < 0:
                raise SystemExit("A16 normalization: the A16_INTENDED_CHANGE literal is not closed")
            hi = end + len(A16_RECORD_END)
    else:
        lo = text.find(A16_PLAN_BEGIN)
        hi = text.find(A16_PLAN_END) + len(A16_PLAN_END) if lo >= 0 else None
    if lo < 0 or hi is None:
        return text  # record deleted (the next task's first commit): nothing to mask
    inner, n = A16_SELF_REFERENCE.subn(r'\1"MASKED"', text[lo:hi])
    if n != A16_EXPECTED_MASKS[kind]:
        raise SystemExit(f"A16 normalization: {n} self-referential literals in the {kind}, "
                         f"expected {A16_EXPECTED_MASKS[kind]}")  # fmt: skip
    return text[:lo] + inner + text[hi:]


def a16_source_entries(suite: str = "lv4") -> dict[str, bytes]:
    """``{relative path: bytes}`` of the hashed source set of ``suite`` as hashed for the A16
    binding: ``run_suite.py`` and the V3-003D plan in the normalized form of :func:`a16_normalize`."""
    out: dict[str, bytes] = {}
    for f in source_files(suite):
        rel = str(f.relative_to(REPO))
        raw = f.read_bytes()
        if rel in (A16_RUN_SUITE_FILE, A16_PLAN_FILE):
            kind = "record" if rel == A16_RUN_SUITE_FILE else "plan"
            raw = a16_normalize(raw.decode(), kind).encode()
        out[rel] = raw
    return out


def a16_digest_of(entries: dict[str, bytes]) -> str:
    """sha256 over the sorted ``path sha256`` lines of ``entries``."""
    lines = [f"{rel} {hashlib.sha256(b).hexdigest()}" for rel, b in entries.items()]
    return hashlib.sha256("\n".join(sorted(lines)).encode()).hexdigest()


def a16_source_digest(suite: str = "lv4") -> str:
    """The A16 source digest: :func:`a16_digest_of` of :func:`a16_source_entries`, i.e. every
    hashed source file of ``suite`` including ``run_suite.py`` (the record, its bounds and this
    function) and the V3-003D plan block, with only the two self-referential digest values masked.
    The A16 intended-change record is valid only for the source state with this digest: any later
    change to a hashed file requires a conscious refresh (or deletion) of the record."""
    return a16_digest_of(a16_source_entries(suite))


def environment_text(
    sha: str, source: str, dirty: str, args: argparse.Namespace, workers: int, single: bool = False
) -> str:
    import numpy
    import warp

    manifest_path = getattr(args, "partials_manifest", None)
    manifest_sha = ""
    if manifest_path and Path(manifest_path).is_file():
        manifest_sha = hashlib.sha256(Path(manifest_path).read_bytes()).hexdigest()
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
        f"partials_manifest={manifest_path or ''}",
        f"partials_manifest_sha256={manifest_sha}",
        "source_hashes:",
    ]
    lines += [f"  {sha256(f)}  {f.relative_to(REPO)}" for f in source_files(args.suite)]
    return "\n".join(lines) + "\n"


def resolve_workers(value: str) -> int:
    """``auto`` is every logical CPU of this host; otherwise an integer."""
    return (os.cpu_count() or 1) if value == "auto" else int(value)


def main(argv: list[str] | None = None) -> int:
    if "--print-a16-source-digest" in (argv if argv is not None else sys.argv[1:]):
        print(a16_source_digest())
        return 0
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
        "20391004 for lv/hr and 20401004 for lv4/hr4 (the 2041xxxx rehearsals, 20351004 and the consumed "
        "20361004, 20371004, 20381004 give non-conformant archives); for lv/hr the "
        "bases 20261004, 20271004, 20281004, 20291004, 20301004, 20311004, 20321004, 20331004 and "
        "20341004 give non-conformant archives); "
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
        "--partials-manifest",
        default=None,
        help="lv5 with --import-dirs: JSON file mapping each imported partial's file name to its "
        "content_sha256, written by the orchestrator from the PARTIAL lines of the host-runner "
        "records of the shard steps",
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
    if args.suite == "lv5b" and args.seed_base in V5B_CONSUMED_SEED_BASES:
        raise SystemExit(
            f"seed base {args.seed_base} is consumed for lv5b (amendment 11 of the V3-005 plan); "
            f"the qualification base is {V5B_QUALIFICATION_SEED_BASE}"
        )
    workers = 1 if args.single_process or args.suite in ("lv5", "lv5b", "hr5") else resolve_workers(args.workers)
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
        args.partials_manifest,
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
        echo_partial_lines(path)
        print(f"== {name}: exit={code}", flush=True)
        failures += code != 0
    sys.path.insert(0, str(HERE))
    import summarize

    rc = summarize.main([str(out), "--expected-sha", sha])
    print(f"done: {out} (failed_steps={failures})")
    return 1 if failures or rc else 0


if __name__ == "__main__":
    sys.exit(main())
