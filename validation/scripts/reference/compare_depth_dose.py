"""Compare integral depth-dose metrics of materialized reference runs.

Usage: compare_depth_dose.py --runs RUN_DIR... --output OUT.json
Writes metrics (mm) per run and pairwise differences; no raw curves.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from ionmc.reference.metrics import (
    distal_falloff_80_20,
    normalize_to_peak,
    peak_depth,
    r80,
    r90,
)
from ionmc.reference.runs import depth_dose, file_hashes, load_run

METRICS = ("peak_depth_mm", "r80_mm", "r90_mm", "falloff_80_20_mm")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", nargs="+", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    runs = []
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
                "output_sha256": file_hashes(run, dd.files),
                "peak_depth_mm": peak_depth(dd.depth_mm, curve),
                "r80_mm": r80(dd.depth_mm, curve),
                "r90_mm": r90(dd.depth_mm, curve),
                "falloff_80_20_mm": distal_falloff_80_20(dd.depth_mm, curve),
            }
        )
    pairs = []
    for a, b in itertools.combinations(runs, 2):
        pairs.append(
            {
                "a": a["engine"],
                "b": b["engine"],
                **{f"{m}_diff_b_minus_a": b[m] - a[m] for m in METRICS},
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"runs": runs, "pairwise": pairs}, indent=2) + "\n")
    for r in runs:
        print(
            f"{r['engine']:9s} peak {r['peak_depth_mm']:7.2f}  R90 {r['r90_mm']:7.2f}  "
            f"R80 {r['r80_mm']:7.2f}  80-20 {r['falloff_80_20_mm']:5.2f}  "
            f"bin {r['bin_width_mm']:g} mm"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
