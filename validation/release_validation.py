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
failure (3) both count as NOT release-ready. The ``warp_cuda_smoke`` diagnostic is
excluded (it is not a milestone gate and has a non-standard interface).

Exit 0 iff release-ready (all suites pass), 3 otherwise.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

from ionmc import __version__
from ionmc.backend import mathlib

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION_DIR = REPO_ROOT / "validation"
BENCH_DIR = REPO_ROOT / "benchmarks"
BASELINE_DIR = BENCH_DIR / "baselines"

#: milestone validation scripts excluded from the release suite (not a gate)
EXCLUDE = {"release_validation.py", "warp_cuda_smoke.py"}
#: benchmarks whose physics gate + regression check are part of the release suite
BENCHMARKS = ("depth_dose_csda", "dose3d")
BENCH_SCRIPT = {"depth_dose_csda": "bench_depth_dose.py", "dose3d": "bench_dose3d.py"}


def _supports(path: Path, flag: str) -> bool:
    return flag in path.read_text(encoding="utf-8")


def _run(argv: list[str], timeout: int) -> tuple[int, str]:
    """Run a suite subprocess, returning (returncode, stdout). A timeout or a
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
        return 124, ""
    except OSError as exc:
        sys.stderr.write(f"failed to launch {argv}: {exc}\n")
        return 125, ""
    return proc.returncode, proc.stdout


def _validation_suites(
    cache_dir: str | None, require_cuda: bool
) -> list[dict[str, Any]]:
    suites: list[dict[str, Any]] = []
    for script in sorted(VALIDATION_DIR.glob("v*.py")):
        if script.name in EXCLUDE:
            continue
        argv = [sys.executable, str(script)]
        if require_cuda and _supports(script, "--require-cuda"):
            argv.append("--require-cuda")
        if cache_dir and _supports(script, "--cache-dir"):
            argv += ["--cache-dir", cache_dir]
        suites.append({"name": script.stem, "kind": "validation", "argv": argv})
    return suites


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
    # benchmark physics gate + regression check per benchmark
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
        rc, stdout = _run(suite["argv"], args.timeout)
        passed = rc == 0
        entry = {
            "name": suite["name"],
            "kind": suite["kind"],
            "returncode": rc,
            "passed": passed,
        }
        # for a benchmark, also run the regression check against the committed baseline
        if suite["kind"] == "benchmark" and stdout:
            bname = suite["name"].split(":", 1)[1]
            baseline = BASELINE_DIR / f"{bname}.json"
            report_path = REPO_ROOT / f".release_report_{bname}.json"
            try:
                report_path.write_text(stdout, encoding="utf-8")
                rrc, _ = _run(
                    [
                        sys.executable,
                        str(BENCH_DIR / "check_regression.py"),
                        "--report",
                        str(report_path),
                        "--baseline",
                        str(baseline),
                    ],
                    args.timeout,
                )
            finally:
                report_path.unlink(missing_ok=True)
            reg_passed = rrc == 0
            entry["regression_returncode"] = rrc
            entry["regression_passed"] = reg_passed
            passed = passed and reg_passed
            entry["passed"] = passed
        all_passed = all_passed and passed
        report["results"].append(entry)

    report["all_passed"] = all_passed
    report["n_suites"] = len(report["results"])
    report["n_passed"] = sum(1 for r in report["results"] if r["passed"])
    report["release_ready"] = all_passed
    json.dump(report, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if all_passed else 3


if __name__ == "__main__":
    raise SystemExit(main())
