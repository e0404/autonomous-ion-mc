"""Provenance-bound measured calibration-part rate of the V7-R2 rep-based plans (plan R-4b, Amendment 20 (f)8;
Codex REVIEW-fe47f290 finding 2). Standard library only: ``run_suite.py`` (preflight) and ``steps_v5d.py``
(throughput gate) import it without ``ionmc``.

The rep-based V7-R2 plans (``v7r.r2_rep_plans``: 3200 repetitions over 39 parts per surrogate, the kappa parts) rest on
a slow-host seconds-per-repetition rate whose default (``v7r.R2_CAL_S_PER_REP`` = 12.2 s) is EXTRAPOLATED, and the 39-part
layout has 1.6 % headroom to the 1100 s planning target. The lv5d suite therefore accepts the rep plans only from a
MEASURED rate record, a JSON file named by the environment variable ``IONMC_V7R2_CAL_RATE_RECORD``::

    {"schema": "ionmc-v7r2-cal-rate-1",
     "s_per_rep": <wall_mc_s / n of the measured part>,   "fixed_s": <per-part fixed cost, >= the 70 s default>,
     "reps": <repetitions of the measured part (n)>,       "wall_mc_s": <its Monte Carlo wall time>,
     "wall_total_s": <pytest wall time of the part>,
     "host": {"node": ..., "cpu_count": ..., "platform": ..., "python": ...},
     "git_sha": <40 hex of the measured tree>,             "log_dir": <directory with the pytest log and part record>,
     "part_record": <file name of the sealed calibration part record>, "part_record_sha256": <64 hex>,
     "part_log": <file name of the pytest log>,            "part_log_sha256": <64 hex>}

``python v7r2_cal_rate.py --part-record ... --part-log ... --log-dir ... --git-sha ... --out ...`` builds the record
from the sealed part record (``wall_mc_s``, ``n``) and the pytest log (``1 passed in X s``); the rate is the Monte
Carlo wall time per repetition and the fixed cost is ``max(70 s, wall_total - wall_mc)`` (conservative). No
measurement is made here."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import sys
from pathlib import Path
from typing import Any

CAL_RATE_ENV = "IONMC_V7R2_CAL_RATE_RECORD"
CAL_RATE_SCHEMA = "ionmc-v7r2-cal-rate-1"
MIN_FIXED_S = 70.0  # v7r.R2_CAL_FIXED_S (the planning value; the record may only be more conservative)
NOT_MEASURED = "not ok: calibration rate not measured"
_HEX40 = re.compile(r"[0-9a-f]{40}")
_HEX64 = re.compile(r"[0-9a-f]{64}")


def _number(rec: dict[str, Any], key: str, *, positive: bool) -> float:
    v = rec.get(key)
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise ValueError(f"{key} is not a finite number")
    if (v <= 0.0) if positive else (v < 0.0):
        raise ValueError(f"{key} must be {'> 0' if positive else '>= 0'}")
    return float(v)


def validate_record(rec: Any) -> dict[str, Any]:
    """The record itself or ``ValueError`` (every field of the module docstring is required and checked)."""
    if not isinstance(rec, dict):
        raise ValueError("the record is not a JSON object")
    if rec.get("schema") != CAL_RATE_SCHEMA:
        raise ValueError(f"schema is {rec.get('schema')!r}, expected {CAL_RATE_SCHEMA!r}")
    s = _number(rec, "s_per_rep", positive=True)
    fixed = _number(rec, "fixed_s", positive=False)
    wall = _number(rec, "wall_mc_s", positive=True)
    _number(rec, "wall_total_s", positive=True)
    reps = rec.get("reps")
    if isinstance(reps, bool) or not isinstance(reps, int) or reps < 1:
        raise ValueError("reps must be a positive integer")
    if abs(s - wall / reps) > 1e-9 * max(1.0, s):
        raise ValueError("s_per_rep is not wall_mc_s / reps")
    if fixed < MIN_FIXED_S:
        raise ValueError(f"fixed_s {fixed} is below the planning value {MIN_FIXED_S}")
    host = rec.get("host")
    if not (isinstance(host, dict) and isinstance(host.get("node"), str) and host["node"]):
        raise ValueError("host must be an object with a non-empty node")
    if not (isinstance(rec.get("git_sha"), str) and _HEX40.fullmatch(rec["git_sha"])):
        raise ValueError("git_sha must be 40 lowercase hex digits")
    if not (isinstance(rec.get("log_dir"), str) and rec["log_dir"]):
        raise ValueError("log_dir must be a non-empty string")
    for key in ("part_record", "part_log"):
        if not (isinstance(rec.get(key), str) and rec[key]):
            raise ValueError(f"{key} must be a non-empty file name")
        if not (isinstance(rec.get(f"{key}_sha256"), str) and _HEX64.fullmatch(rec[f"{key}_sha256"])):
            raise ValueError(f"{key}_sha256 must be 64 lowercase hex digits")
    return rec


def load_record(path: str | os.PathLike[str]) -> dict[str, Any]:
    """``{"record", "sha256", "path"}`` of the rate record file, or ``ValueError`` (absent, unreadable, malformed)."""
    p = Path(path)
    if not p.is_file():
        raise ValueError(f"{p} is not a file")
    raw = p.read_bytes()
    try:
        rec = validate_record(json.loads(raw))
    except ValueError as exc:  # JSON errors are ValueErrors as well
        raise ValueError(f"calibration rate record {p} invalid: {exc}") from exc
    return {"record": rec, "sha256": hashlib.sha256(raw).hexdigest(), "path": str(p)}


def record_from_env() -> dict[str, Any] | None:
    """:func:`load_record` of ``IONMC_V7R2_CAL_RATE_RECORD`` (``None`` if the variable is unset or empty)."""
    path = os.environ.get(CAL_RATE_ENV, "")
    return load_record(path) if path else None


def build_record(part_record: Path, part_log: Path, log_dir: Path, git_sha: str) -> dict[str, Any]:
    """The record of one measured calibration part (see the module docstring)."""
    doc = json.loads(part_record.read_text())
    n, wall = doc.get("n"), doc.get("wall_mc_s")
    if not (isinstance(n, int) and isinstance(wall, (int, float)) and n >= 1 and wall > 0):
        raise SystemExit(f"{part_record}: no n / wall_mc_s of a calibration part")
    text = part_log.read_text()
    m = re.search(r"\b1 passed\b.* in ([0-9.]+)s", text)
    if not m:
        raise SystemExit(f"{part_log}: no '1 passed ... in X s' line (the part must pass exactly once)")
    total = float(m.group(1))
    rec = {
        "schema": CAL_RATE_SCHEMA, "s_per_rep": float(wall) / n, "fixed_s": max(MIN_FIXED_S, total - float(wall)),
        "reps": n, "wall_mc_s": float(wall), "wall_total_s": total,
        "host": {"node": platform.node(), "cpu_count": os.cpu_count(), "platform": platform.platform(),
                 "python": platform.python_version()},
        "git_sha": git_sha, "log_dir": str(log_dir),
        "part_record": part_record.name, "part_record_sha256": hashlib.sha256(part_record.read_bytes()).hexdigest(),
        "part_log": part_log.name, "part_log_sha256": hashlib.sha256(part_log.read_bytes()).hexdigest(),
    }  # fmt: skip
    return validate_record(rec)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--part-record", required=True, help="sealed v7r2-calibration-<surrogate>-p<NN>.json of the part")
    ap.add_argument("--part-log", required=True, help="the logged pytest output of that single part")
    ap.add_argument("--log-dir", required=True, help="directory holding both files (recorded)")
    ap.add_argument("--git-sha", required=True, help="full SHA of the tree the part was measured at")
    ap.add_argument("--out", required=True, help="rate record to write (refuses to overwrite)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists")
    rec = build_record(Path(a.part_record), Path(a.part_log), Path(a.log_dir), a.git_sha)
    out.write_text(json.dumps(rec, sort_keys=True, indent=1) + "\n")
    print(f"{out}: s_per_rep {rec['s_per_rep']:.4f}, fixed_s {rec['fixed_s']:.1f}, "
          f"sha256 {load_record(out)['sha256']}")  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
