"""Compare integral depth-dose metrics of materialized reference runs.

Usage: compare_depth_dose.py --runs RUN_DIR... --output OUT.json
Writes metrics (mm) per run and pairwise differences; no raw curves.

The evidence status is always "exploratory": each input is a single-seed run with no statistical
uncertainty. Batch statistics (several runs per engine with distinct seeds read from the manifested
case.json) are computed by compare_batches.py; there is deliberately no option to relabel the
status from detached
metadata. The analysis code SHA is derived from git HEAD of this repository and a dirty flag is
recorded; --code-sha, if given, must equal HEAD.
"""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from ionmc.reference.metrics import (
    distal_falloff_80_20,
    normalize_to_peak,
    peak_depth,
    r80,
    r90,
)
from ionmc.reference.runs import depth_dose, file_hashes, load_run

EXPLORATORY_PREFIX = "exploratory: single-seed runs, no statistical uncertainty"


def exploratory_status(runs: list[dict[str, Any]]) -> str:
    """Status text with the bin widths and distal widths derived from the loaded runs."""

    def fmt(key: str) -> str:
        vals = sorted({f"{r[key]:.1f}" for r in runs}, key=float)
        return "{" + ", ".join(vals) + "}"

    return (
        f"{EXPLORATORY_PREFIX}; depth bin widths {fmt('bin_width_mm')} mm; "
        f"distal 80-20 widths {fmt('falloff_80_20_mm')} mm"
    )


METRICS = ("peak_depth_mm", "r80_mm", "r90_mm", "falloff_80_20_mm")


def _git(*args: str) -> str:
    repo = Path(__file__).resolve().parent
    out = subprocess.run(
        ["git", "-C", str(repo), *args], check=True, capture_output=True, text=True
    )
    return out.stdout.strip()


def code_state(claimed: str | None) -> tuple[str, bool]:
    """Return (HEAD sha, dirty); exit non-zero if a claimed SHA differs from HEAD."""
    head = _git("rev-parse", "HEAD")
    dirty = bool(_git("status", "--porcelain"))
    if claimed is not None and claimed != head:
        raise SystemExit(f"--code-sha {claimed} does not equal repository HEAD {head}")
    return head, dirty


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", nargs="+", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument(
        "--code-sha",
        default=None,
        help="optional assertion: must equal git HEAD of this repository (HEAD is recorded)",
    )
    args = ap.parse_args(argv)

    code_sha, dirty = code_state(args.code_sha)
    runs: list[dict[str, Any]] = []
    for run_dir in args.runs:
        run = load_run(run_dir)
        dd = depth_dose(run)
        curve = normalize_to_peak(dd.dose)
        runs.append(
            {
                "run_id": run.run_id,
                "engine": run.engine,
                "source_sha": run.source_sha,
                "case_input": run.case.get("input"),
                "histories": dd.histories,
                "dose_unit": dd.unit,
                "bin_width_mm": float(dd.depth_mm[1] - dd.depth_mm[0]),
                "output_sha256": file_hashes(run, [*dd.files, "inputs/case.json"]),
                "peak_depth_mm": peak_depth(dd.depth_mm, curve),
                "r80_mm": r80(dd.depth_mm, curve),
                "r90_mm": r90(dd.depth_mm, curve),
                "falloff_80_20_mm": distal_falloff_80_20(dd.depth_mm, curve),
            }
        )
    pairs: list[dict[str, Any]] = []
    for a, b in itertools.combinations(runs, 2):
        pairs.append(
            {
                "a": a["engine"],
                "b": b["engine"],
                **{f"{m}_diff_b_minus_a": b[m] - a[m] for m in METRICS},
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    doc = {
        "evidence_status": exploratory_status(runs),
        "analysis_code_sha": code_sha,
        "analysis_code_dirty": dirty,
        "runs": runs,
        "pairwise": pairs,
    }
    args.output.write_text(json.dumps(doc, indent=2) + "\n")
    for r in runs:
        print(
            f"{r['engine']:9s} peak {r['peak_depth_mm']:7.2f}  R90 {r['r90_mm']:7.2f}  "
            f"R80 {r['r80_mm']:7.2f}  80-20 {r['falloff_80_20_mm']:5.2f}  "
            f"bin {r['bin_width_mm']:g} mm"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
