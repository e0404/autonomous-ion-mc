"""Release validation suite (Stage 6, decision 0033).

`EXPERIMENT.md` requires that a release from `develop` to `main` has *passed the
defined release validation suite*, with documented validation status and reproducible
benchmark information. This script **is** that suite: it runs, on the controlled host
runner, every milestone validation (V0-V5) plus the benchmark physics gates and the
performance-regression checks against the committed baselines (V6), and aggregates a
single machine-readable release-readiness report. It passes only if **every** suite
passes.

Each milestone validation is an independent CLI following the shared exit convention
(0 pass, 3 gate failure, 4 dataset missing); this orchestrator runs each as a
subprocess with only the flags it supports, so a dataset-missing (4) or a gate
failure (3) both count as NOT release-ready. Composition is discovered by globbing
``validation/v*.py`` -- a pattern that by construction excludes both the
``warp_cuda_smoke`` diagnostic (not a milestone gate) and this orchestrator itself.

Exit 0 iff release-ready (all suites pass), 3 otherwise.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from ionmc import __version__
from ionmc.backend import mathlib

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION_DIR = REPO_ROOT / "validation"
BENCH_DIR = REPO_ROOT / "benchmarks"
BASELINE_DIR = BENCH_DIR / "baselines"

#: benchmarks whose physics gate + regression check are part of the release suite
BENCHMARKS = ("depth_dose_csda", "dose3d")
BENCH_SCRIPT = {"depth_dose_csda": "bench_depth_dose.py", "dose3d": "bench_dose3d.py"}
#: how much of a failing suite's output to keep in the report, for diagnosis
DIAG_CHARS = 800


def _supports(path: Path, flag: str) -> bool:
    return flag in path.read_text(encoding="utf-8")


def _run(argv: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a suite subprocess, returning (returncode, stdout, stderr). A timeout or a
    launch failure is reported as a non-zero code so it counts as not-ready."""
    try:
        proc = subprocess.run(
            argv,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=REPO_ROOT,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return 124, "", f"timeout after {timeout}s"
    except OSError as exc:
        return 125, "", f"failed to launch: {exc}"
    return proc.returncode, proc.stdout, proc.stderr


def _validation_suites(
    cache_dir: str | None, require_cuda: bool
) -> list[dict[str, Any]]:
    suites: list[dict[str, Any]] = []
    for script in sorted(VALIDATION_DIR.glob("v*.py")):
        argv = [sys.executable, str(script)]
        if require_cuda and _supports(script, "--require-cuda"):
            argv.append("--require-cuda")
        if cache_dir and _supports(script, "--cache-dir"):
            argv += ["--cache-dir", cache_dir]
        suites.append({"name": script.stem, "kind": "validation", "argv": argv})
    return suites


def _regression_check(bname: str, stdout: str, timeout: int) -> tuple[int, str]:
    """Run the performance-regression check on a benchmark's stdout report against
    its committed baseline. Missing/empty stdout is a failure (fail closed): the V6
    baseline comparison is a required part of the release gate."""
    if not stdout.strip():
        return 126, "benchmark produced no report to regression-check"
    baseline = BASELINE_DIR / f"{bname}.json"
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False, encoding="utf-8"
    ) as fh:
        fh.write(stdout)
        report_path = Path(fh.name)
    try:
        rc, _, err = _run(
            [
                sys.executable,
                str(BENCH_DIR / "check_regression.py"),
                "--report",
                str(report_path),
                "--baseline",
                str(baseline),
            ],
            timeout,
        )
    finally:
        report_path.unlink(missing_ok=True)
    return rc, err


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-cuda", action="store_true")
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--only", default=None, help="regex; run only matching suites")
    parser.add_argument("--timeout", type=int, default=300, help="per-suite seconds")
    args = parser.parse_args()

    report: dict[str, Any] = {
        "schema_version": 1,
        "suite": "release_validation",
        "experiment_config": {
            "ionmc_version": __version__,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "warp_version": (
                str(mathlib.warp_module().config.version) if mathlib.HAVE_WARP else None
            ),
            "require_cuda": args.require_cuda,
        },
        "results": [],
    }

    suites = _validation_suites(args.cache_dir, args.require_cuda)
    for name in BENCHMARKS:
        script = BENCH_DIR / BENCH_SCRIPT[name]
        argv = [sys.executable, str(script)]
        if args.require_cuda:
            argv.append("--require-cuda")
        if args.cache_dir:
            argv += ["--cache-dir", args.cache_dir]
        suites.append({"name": f"bench:{name}", "kind": "benchmark", "argv": argv})

    if args.only:
        pat = re.compile(args.only)
        suites = [s for s in suites if pat.search(s["name"])]

    all_passed = True
    for suite in suites:
        rc, stdout, stderr = _run(suite["argv"], args.timeout)
        passed = rc == 0
        entry: dict[str, Any] = {
            "name": suite["name"],
            "kind": suite["kind"],
            "returncode": rc,
            "passed": passed,
        }
        if suite["kind"] == "benchmark":
            # the regression check is a REQUIRED part of the release gate (fail closed)
            bname = suite["name"].split(":", 1)[1]
            rrc, rerr = _regression_check(bname, stdout, args.timeout)
            entry["regression_returncode"] = rrc
            entry["regression_passed"] = rrc == 0
            passed = passed and rrc == 0
            entry["passed"] = passed
            if rrc != 0 and rerr:
                entry["regression_diagnostic"] = rerr[:DIAG_CHARS]
        if not passed:
            # keep a bounded tail of the failing suite's output for diagnosis
            tail = (stderr or stdout or "").strip()
            if tail:
                entry["diagnostic"] = tail[-DIAG_CHARS:]
        all_passed = all_passed and passed
        report["results"].append(entry)

    report["all_passed"] = all_passed
    report["n_suites"] = len(report["results"])
    report["n_passed"] = sum(1 for r in report["results"] if r["passed"])
    report["failed"] = [r["name"] for r in report["results"] if not r["passed"]]
    report["release_ready"] = all_passed
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if all_passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
