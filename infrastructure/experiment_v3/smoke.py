"""Native scorer smoke checks: finite nonzero output, not physical accuracy."""

import argparse
import json
import math
import struct
from pathlib import Path

from infrastructure.experiment_v3.common import ROOT, STATE, safe_path, write_json
from infrastructure.experiment_v3.reference import run_reference


def score_summary(engine, path):
    if engine == "topas":
        values = []
        for line in path.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            columns = line.replace(",", " ").split()
            values.append(float(columns[3]))
    else:
        blob = path.read_bytes()
        marker = b"ElementDataFile = "
        start = blob.index(marker)
        end = blob.index(b"\n", start)
        name = blob[start + len(marker) : end].decode().strip()
        header = blob[:end].decode()
        if (
            "ElementType = MET_FLOAT" not in header
            or "ElementByteOrderMSB = True" in header
            or "BinaryDataByteOrderMSB = True" in header
        ):
            raise ValueError("Smoke parser expects little-endian MET_FLOAT")
        raw = (
            blob[end + 1 :]
            if name == "LOCAL"
            else safe_path(path.parent, name).read_bytes()
        )
        values = [v[0] for v in struct.iter_unpack("<f", raw)]
    if (
        not values
        or not all(math.isfinite(v) and v >= 0 for v in values)
        or not any(v > 0 for v in values)
    ):
        raise ValueError(
            "Scoring smoke produced empty, nonfinite, negative or all-zero data"
        )
    return {
        "bins": len(values),
        "nonzero_bins": sum(v > 0 for v in values),
        "sum_native_units": sum(values),
        "physical_accuracy_validated": False,
    }


def run_smokes(root=ROOT, state=STATE):
    root, state = Path(root), Path(state)
    results = {}
    for engine in ("topas", "mcsquare", "fred"):
        result = run_reference(
            root,
            "experiment/v3/references/" + engine,
            engine,
            registry=state / "engines.json",
            state=state,
            timeout=120,
        )
        if result["succeeded"]:
            case = json.loads(
                (root / "experiment/v3/references" / engine / "case.json").read_text()
            )
            output = (
                state
                / "references"
                / result["run_id"]
                / "work"
                / case["expected_outputs"][0]
            )
            try:
                result["score_smoke"] = score_summary(engine, output)
            except (ValueError, OSError, IndexError, struct.error) as exc:
                result["succeeded"] = False
                result["score_error"] = str(exc)
        results[engine] = {
            k: v for k, v in result.items() if k not in ("stdout", "stderr")
        }
    return {
        "all_execution_smokes_passed": all(r["succeeded"] for r in results.values()),
        "engines": results,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = run_smokes()
    write_json(args.output, result)
    print(json.dumps(result, indent=2))
    return 0 if result["all_execution_smokes_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
