"""Performance-regression check for a benchmark report (Stage 6, decision 0032).

Compares a fresh benchmark report (produced by a driver such as
``bench_depth_dose.py`` or ``bench_dose3d.py``) against a committed baseline under
``benchmarks/baselines/``. It **gates on the reference (float64) physics digest** —
the hardware-independent "unchanged scientific outcome" of milestone V6 — and
**reports per-backend throughput ratios** (current / baseline) without gating on
them (wall-clock is machine-dependent). This is the before/after guard for every
future optimisation.

Usage::

    # produce a report, then check it against the committed baseline
    python benchmarks/bench_dose3d.py --require-cuda --cache-dir /cache/ionmc > r.json
    python benchmarks/check_regression.py --report r.json \\
        --baseline benchmarks/baselines/dose3d.json

    # (re)generate a baseline from a report (run on a trusted reference machine)
    python benchmarks/check_regression.py --report r.json --emit-baseline > b.json

Exit 0 if the physics matches the baseline, 3 on a physics regression, 4 on a
missing/invalid input file.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from ionmc.benchmarking import compare_to_baseline, dump_report, make_baseline


def _load(path: str) -> dict[str, Any] | None:
    try:
        with open(path, encoding="utf-8") as fh:
            obj = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None
    return obj if isinstance(obj, dict) else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True, help="benchmark report JSON")
    parser.add_argument("--baseline", help="committed baseline JSON to compare against")
    parser.add_argument(
        "--emit-baseline",
        action="store_true",
        help="print a baseline built from --report instead of comparing",
    )
    args = parser.parse_args()

    report = _load(args.report)
    if report is None:
        sys.stderr.write(f"cannot read report: {args.report}\n")
        return 4

    if args.emit_baseline:
        dump_report(make_baseline(report), sys.stdout)
        return 0

    if not args.baseline:
        sys.stderr.write("--baseline is required unless --emit-baseline is given\n")
        return 4
    baseline = _load(args.baseline)
    if baseline is None:
        sys.stderr.write(f"cannot read baseline: {args.baseline}\n")
        return 4

    result = compare_to_baseline(report, baseline)
    dump_report(result, sys.stdout)
    return 0 if result["physics_ok"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
