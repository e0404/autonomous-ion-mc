"""Batch statistics of reference runs: several seeded runs per engine (evidence-grade protocol).

Usage: compare_batches.py --runs RUN_DIR... --output OUT.json [--code-sha SHA]

Runs are grouped by the engine recorded in each run's manifest. For every engine at least
MIN_RUNS runs with DISTINCT seeds are required. A run's seed is read only from its manifested
``inputs/case.json`` (``seeds`` must hold exactly one integer) and cross-checked against the
native input that was actually run (TOPAS ``Ts/Seed``, MCsquare ``RNG_Seed``, FRED ``-rseed``),
both verified against the transfer manifest; no CLI or detached metadata can supply seeds.
Within an engine all runs must share source commit, case ``input`` name, histories
(>= MIN_HISTORIES) and depth-bin width (<= MAX_BIN_MM). Any violation exits non-zero: this
script has no exploratory fallback (use compare_depth_dose.py for single runs). Statistics:
per-engine mean, sample SD,
standard error SD/sqrt(n) and a two-sided 95 % Student-t interval for peak depth, R90, R80, the
distal 80-20 width and, for runs with a 3-D dose (TOPAS lateral case), the Sheppard-corrected
lateral sigma at z/R80 in {0.5, 0.9} for the all-particle and (if scored)
the generation-0 "primary" 3-D dose (full-field window and a fixed +-20 mm window; see
ionmc.reference.metrics). The output status "batched" means only that these conditions hold; it is
not a grade against any frozen acceptance criterion. The analysis code SHA is git HEAD with a
dirty flag; --code-sha, if given, must equal HEAD.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from ionmc.reference.metrics import (
    distal_falloff_80_20,
    lateral_sigma2_at_depth_fractions,
    normalize_to_peak,
    peak_depth,
    r80,
    r90,
)
from ionmc.reference.runs import (
    ReferenceRun,
    RunError,
    depth_dose,
    dose_3d,
    file_hashes,
    load_run,
    read_verified,
)

MIN_RUNS = 2
MIN_HISTORIES = 100_000
MAX_BIN_MM = 0.5
FRACTIONS = (0.5, 0.9)
# 3-D scorer name -> record key (all particles; generation-0 particles only)
LATERAL_SCORERS = {"dose3d": "lateral", "dose3d_primary": "lateral_primary"}
WINDOWS: dict[str, float | None] = {"full": None, "w20": 20.0}
STATUS = "batched"
STATUS_NOTE = (
    "batched: >= 2 runs per engine with distinct seeds from manifested case.json; "
    "not graded against any frozen acceptance criterion"
)
METRICS = ("peak_depth_mm", "r80_mm", "r90_mm", "falloff_80_20_mm")

# two-sided 95 % Student-t quantiles, df = 1..30
_T95 = (
    12.7062,
    4.3027,
    3.1824,
    2.7764,
    2.5706,
    2.4469,
    2.3646,
    2.3060,
    2.2622,
    2.2281,
    2.2010,
    2.1788,
    2.1604,
    2.1448,
    2.1314,
    2.1199,
    2.1098,
    2.1009,
    2.0930,
    2.0860,
    2.0796,
    2.0739,
    2.0687,
    2.0639,
    2.0595,
    2.0555,
    2.0518,
    2.0484,
    2.0452,
    2.0423,
)


class BatchError(ValueError):
    """The runs do not satisfy the batch protocol."""


def t95(df: int) -> float:
    if df < 1:
        raise BatchError("need at least 2 values for an interval")
    return _T95[df - 1] if df <= len(_T95) else 1.960


def summarize(values: list[float]) -> dict[str, Any]:
    """Mean, sample SD, standard error and 95 % t-interval of per-run values."""
    n = len(values)
    if n < MIN_RUNS:
        raise BatchError("need at least 2 runs")
    a = np.asarray(values, dtype=np.float64)
    mean = float(a.mean())
    sd = float(a.std(ddof=1))
    se = sd / math.sqrt(n)
    h = t95(n - 1) * se
    return {"n": n, "mean": mean, "sd": sd, "se": se, "ci95": [mean - h, mean + h]}


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


_NATIVE_SEED = {
    "topas": re.compile(r"^\s*i:Ts/Seed\s*=\s*(\d+)\s*$", re.M),
    "mcsquare": re.compile(r"^\s*RNG_Seed\s+(\d+)\s*$", re.M),
}


def run_seed(run: ReferenceRun) -> int:
    """Seed from the manifested case.json, cross-checked against the native input (fail closed)."""
    seeds = run.case.get("seeds")
    if not isinstance(seeds, list) or len(seeds) != 1 or type(seeds[0]) is not int or seeds[0] <= 0:
        raise BatchError(f"{run.run_id}: case.json must hold exactly one positive integer seed")
    seed: int = seeds[0]
    native = run.case.get("input")
    if not isinstance(native, str):
        raise BatchError(f"{run.run_id}: case.json has no input file name")
    rel = f"inputs/{native}"
    if rel not in run.files:
        raise BatchError(f"{run.run_id}: native input {rel} not in transfer manifest")
    try:
        text = read_verified(run, rel).decode("utf-8", errors="replace")
    except RunError as exc:
        raise BatchError(str(exc)) from exc
    if run.engine in _NATIVE_SEED:
        found = [int(m) for m in _NATIVE_SEED[run.engine].findall(text)]
        if found != [seed]:
            raise BatchError(f"{run.run_id}: native input seed {found} != case.json seed {seed}")
    elif run.engine == "fred":
        if "rseed" in text:
            raise BatchError(f"{run.run_id}: fred input must not set the seed ({rel})")
        try:
            found_seed, _ = fred_rseed(list(run.case.get("arguments", [])))
        except BatchError as exc:
            raise BatchError(f"{run.run_id}: {exc}") from exc
        if found_seed != seed:
            raise BatchError(f"{run.run_id}: -rseed {found_seed} != case.json seed {seed}")
    else:
        raise BatchError(f"{run.run_id}: unsupported engine {run.engine!r}")
    return seed


_SEED_LINE = {
    "topas": re.compile(r"^[ \t]*i:Ts/Seed[ \t]*=.*(?:\n|$)", re.M),
    "mcsquare": re.compile(r"^[ \t]*RNG_Seed[ \t].*(?:\n|$)", re.M),
}
_NATIVE_HISTORIES = {
    "topas": re.compile(r"^[ \t]*i:So/Beam/NumberOfHistoriesInRun[ \t]*=[ \t]*(\d+)[ \t]*$", re.M),
    "mcsquare": re.compile(r"^[ \t]*Num_Primaries[ \t]+(\d+)[ \t]*$", re.M),
    "fred": re.compile(r"^[ \t]*nprim[ \t]*=[ \t]*(\d+)[ \t]*$", re.M),
}
_STDOUT_HISTORIES = {
    "topas": re.compile(r"Particle source Beam: Total number of histories: (\d+)"),
    "mcsquare": re.compile(r"Nbr primaries simulated: (\d+)"),
    "fred": re.compile(r"Num of primaries to deliver: (\d+)"),
}
_IDENTITY_KEYS = (
    "engine",
    "os_release_sha256",
    "runner_source_sha256",
    "sandbox_binary_sha256",
    "gpu",
    "dirty",
)
# case.json fields allowed to differ between seed variants: the seed list and the rationale text
_CASE_SEED_FIELDS = ("seeds", "rationale")


def _sha(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def fred_rseed(args: list[Any]) -> tuple[int, list[Any]]:
    """The single well-formed ``-rseed <int>`` pair of a FRED argv and the argv without it.

    Zero or several ``-rseed`` options, a missing value or a non-positive-integer value raise
    ``BatchError``.
    """
    idx = [i for i, a in enumerate(args) if a == "-rseed"]
    if len(idx) != 1:
        raise BatchError(f"expected exactly one '-rseed' argument, found {len(idx)}")
    i = idx[0]
    if i + 1 >= len(args) or not re.fullmatch(r"[1-9][0-9]*", str(args[i + 1])):
        raise BatchError("'-rseed' must be followed by a positive integer")
    return int(args[i + 1]), args[:i] + args[i + 2 :]


def run_fingerprint(run: ReferenceRun) -> dict[str, Any]:
    """Everything that defines the executed configuration except the seed (fail closed).

    Returns the normalized-input hash, the hashes of all other manifested inputs, the engine and
    runtime identity from request.json, and the label of the replicate group. Declared histories
    are cross-checked against the native input, FRED's -nprim and the engine's own run summary.
    """
    name = run.case.get("input")
    if not isinstance(name, str):
        raise BatchError(f"{run.run_id}: case.json has no input file name")
    rel = f"inputs/{name}"
    try:
        native = read_verified(run, rel).decode("utf-8", errors="replace")
        aux = {
            r: _sha(read_verified(run, r))
            for r in sorted(run.files)
            if r.startswith("inputs/") and r not in ("inputs/case.json", rel)
        }
        stdout = read_verified(run, "stdout.txt").decode("utf-8", errors="replace")
    except RunError as exc:
        raise BatchError(str(exc)) from exc
    if not run.request:
        raise BatchError(f"{run.run_id}: request.json missing (engine identity unknown)")
    if run.engine not in _NATIVE_HISTORIES:
        raise BatchError(f"{run.run_id}: unsupported engine {run.engine!r}")
    # histories: declared == native input == run summary
    declared = run.histories
    found = [int(m) for m in _NATIVE_HISTORIES[run.engine].findall(native)]
    if found != [declared]:
        raise BatchError(f"{run.run_id}: native input histories {found} != declared {declared}")
    args = [str(a) for a in run.case.get("arguments", [])]
    if "-nprim" in args and args[args.index("-nprim") + 1 : args.index("-nprim") + 2] != [
        str(declared)
    ]:
        raise BatchError(f"{run.run_id}: -nprim argument != declared histories {declared}")
    summary = [int(m) for m in _STDOUT_HISTORIES[run.engine].findall(stdout)]
    if declared not in summary or any(v != declared for v in summary if v != 0):
        raise BatchError(f"{run.run_id}: engine run summary histories {summary} != {declared}")
    # request.json must describe exactly the manifested inputs
    if run.request.get("case") != run.case:
        raise BatchError(f"{run.run_id}: request.json case differs from inputs/case.json")
    for fname, meta in run.request.get("input_files", {}).items():
        entry = run.files.get(f"inputs/{fname}")
        if not entry or entry.get("sha256") != meta.get("sha256"):
            raise BatchError(f"{run.run_id}: request.json input {fname} != manifested file")
    if run.request.get("dirty") is not False or run.request.get("code_sha") != run.source_sha:
        raise BatchError(f"{run.run_id}: run not made from a clean commit == source_sha")
    # seed-free normalization
    if run.engine in _SEED_LINE:
        native = _SEED_LINE[run.engine].sub("", native)
    case = {k: v for k, v in run.case.items() if k not in _CASE_SEED_FIELDS}
    if run.engine == "fred":
        case["arguments"] = fred_rseed(list(case.get("arguments", [])))[1]
    identity = {k: run.request.get(k) for k in _IDENTITY_KEYS}
    if identity["engine"] is None:
        raise BatchError(f"{run.run_id}: request.json lacks the engine identity")
    parts = {
        "case_sha256": _sha(json.dumps(case, sort_keys=True).encode()),
        "native_normalized_sha256": _sha(native.encode()),
        "aux_inputs_sha256": aux,
    }
    label = run.engine
    if run.engine == "topas":
        mods = re.findall(r"^[ \t]*sv:Ph/Default/Modules[ \t]*=.*$", native, re.M)
        if len(mods) != 1:
            raise BatchError(f"{run.run_id}: expected exactly one Ph/Default/Modules line")
        if "g4h-" not in mods[0] and "g4ion" not in mods[0]:
            label = "topas-emonly"
    return {
        "group": label,
        "config_sha256": _sha(json.dumps(parts, sort_keys=True).encode()),
        "config_parts": parts,
        "identity_sha256": _sha(json.dumps(identity, sort_keys=True).encode()),
    }


def analyse_run(run: ReferenceRun) -> dict[str, Any]:
    dd = depth_dose(run)
    curve = normalize_to_peak(dd.dose)
    bin_mm = float(dd.depth_mm[1] - dd.depth_mm[0])
    rec: dict[str, Any] = {
        "run_id": run.run_id,
        "engine": run.engine,
        "source_sha": run.source_sha,
        "case_input": run.case.get("input"),
        **run_fingerprint(run),
        "seed": run_seed(run),
        "histories": dd.histories,
        "bin_width_mm": bin_mm,
        "peak_depth_mm": peak_depth(dd.depth_mm, curve),
        "r80_mm": r80(dd.depth_mm, curve),
        "r90_mm": r90(dd.depth_mm, curve),
        "falloff_80_20_mm": distal_falloff_80_20(dd.depth_mm, curve),
    }
    files = [*dd.files, "inputs/case.json", f"inputs/{run.case['input']}"]
    if run.engine == "topas" and "work/dose3d.bin" in run.files:
        for scorer, key in LATERAL_SCORERS.items():
            if scorer != "dose3d" and f"work/{scorer}.bin" not in run.files:
                continue
            d3 = dose_3d(run, scorer)
            files = sorted({*files, *d3.files})
            lat: dict[str, Any] = {}
            for wname, w in WINDOWS.items():
                res = lateral_sigma2_at_depth_fractions(
                    d3.dose, d3.bin_mm, rec["r80_mm"], FRACTIONS, window_mm=w
                )
                lat[wname] = {f"{f:g}": v for f, v in res.items()}
            rec[key] = lat
            rec["lateral_bin_mm"] = list(d3.bin_mm)
    rec["output_sha256"] = file_hashes(run, files)
    return rec


def validate_group(engine: str, recs: list[dict[str, Any]]) -> None:
    """Protocol checks for the runs of one engine (raises BatchError)."""
    if len(recs) < MIN_RUNS:
        raise BatchError(f"engine {engine}: {len(recs)} run(s); at least {MIN_RUNS} required")
    seeds = [r["seed"] for r in recs]
    if len(set(seeds)) != len(seeds):
        raise BatchError(f"engine {engine}: duplicate seeds {sorted(seeds)}")
    ids = [r["run_id"] for r in recs]
    if len(set(ids)) != len(ids):
        raise BatchError(f"engine {engine}: the same run given more than once")
    for key in ("source_sha", "case_input", "histories", "bin_width_mm", "identity_sha256"):
        vals = {r[key] for r in recs}
        if len(vals) != 1:
            raise BatchError(f"engine {engine}: runs differ in {key}: {sorted(map(str, vals))}")
    if len({r["config_sha256"] for r in recs}) != 1:
        diff = sorted(
            part
            for part in recs[0]["config_parts"]
            if len({json.dumps(r["config_parts"][part], sort_keys=True) for r in recs}) != 1
        )
        raise BatchError(
            f"engine {engine}: runs differ in the executed configuration beyond the seed: {diff}"
        )
    if recs[0]["histories"] < MIN_HISTORIES:
        raise BatchError(f"engine {engine}: histories {recs[0]['histories']} < {MIN_HISTORIES}")
    if recs[0]["bin_width_mm"] > MAX_BIN_MM + 1e-9:
        raise BatchError(f"engine {engine}: depth bin {recs[0]['bin_width_mm']} mm > {MAX_BIN_MM}")
    has_lat = {tuple(k in r for k in LATERAL_SCORERS.values()) for r in recs}
    if len(has_lat) != 1:
        raise BatchError(f"engine {engine}: only some runs have a 3-D dose")


def engine_statistics(recs: list[dict[str, Any]]) -> dict[str, Any]:
    stats: dict[str, Any] = {m: summarize([r[m] for r in recs]) for m in METRICS}
    for key in LATERAL_SCORERS.values():
        if key not in recs[0]:
            continue
        lat: dict[str, Any] = {}
        for wname in WINDOWS:
            for f in FRACTIONS:
                fk = f"{f:g}"
                s2 = [r[key][wname][fk]["sigma2_mm2"] for r in recs]
                if min(s2) <= 0:
                    raise BatchError(f"non-positive {key} variance in window {wname}, z/R80 {fk}")
                lat[f"{wname}@{fk}"] = {
                    "sigma2_mm2": summarize(s2),
                    "sigma_mm": summarize([math.sqrt(v) for v in s2]),
                    "slab_centre_over_r80": recs[0][key][wname][fk]["slab_centre_over_r80"],
                }
        stats[key] = lat
    return stats


def batch_analysis(run_dirs: list[Path]) -> dict[str, Any]:
    recs = [analyse_run(load_run(d)) for d in run_dirs]
    groups: dict[str, list[dict[str, Any]]] = {}
    for r in recs:
        groups.setdefault(r["group"], []).append(r)
    stats: dict[str, Any] = {}
    for engine, g in sorted(groups.items()):
        validate_group(engine, g)
        stats[engine] = engine_statistics(g)
    pairwise = []
    for a, b in itertools.combinations(sorted(stats), 2):
        entry: dict[str, Any] = {"a": a, "b": b}
        for m in METRICS:
            sa, sb = stats[a][m], stats[b][m]
            se = math.hypot(sa["se"], sb["se"])
            diff = sb["mean"] - sa["mean"]
            entry[m] = {"diff_b_minus_a": diff, "se": se, "z": diff / se if se > 0 else None}
        pairwise.append(entry)
    return {
        "evidence_status": STATUS,
        "status_note": STATUS_NOTE,
        "runs": recs,
        "engines": stats,
        "pairwise": pairwise,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", nargs="+", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--code-sha", default=None, help="optional assertion: must equal git HEAD")
    args = ap.parse_args(argv)
    code_sha, dirty = code_state(args.code_sha)
    try:
        doc = batch_analysis(args.runs)
    except (BatchError, RunError) as exc:
        raise SystemExit(f"batch protocol violated: {exc}") from exc
    doc = {"analysis_code_sha": code_sha, "analysis_code_dirty": dirty, **doc}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=2) + "\n")
    for engine, st in doc["engines"].items():
        line = " ".join(
            f"{m.removesuffix('_mm')} {st[m]['mean']:.3f}+-{st[m]['se']:.3f}" for m in METRICS
        )
        print(f"{engine:9s} n={st['peak_depth_mm']['n']} {line}")
        for lkey in LATERAL_SCORERS.values():
            for key, v in st.get(lkey, {}).items():
                sig = v["sigma_mm"]
                print(f"{'':9s} {lkey}[{key}] {sig['mean']:.4f}+-{sig['se']:.4f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
