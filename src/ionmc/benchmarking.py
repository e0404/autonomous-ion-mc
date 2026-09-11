"""Reproducible performance-benchmarking harness (Stage 6, decision 0030).

Milestone V6 requires *documented benchmark results with unchanged scientific
validation outcomes before and after each optimisation*. This module provides the
shared machinery for that: capturing the full run provenance (hardware, library
versions, backend/device), timing a callable with warm-up and repeats, and pinning
a reproducible **scientific digest** of the benchmarked output so a later
optimisation can be shown to leave the physics unchanged. Benchmark *drivers* live
under ``benchmarks/`` and use this library; they are deliberately kept separate from
the ``validation/`` scientific gates (a benchmark records wall-clock time, which is
hardware-dependent and never asserted against an absolute threshold).
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import statistics
import subprocess
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

from ionmc import __version__
from ionmc.backend import mathlib


def _git_sha() -> str | None:
    """Best-effort git SHA (full ``rev-parse HEAD``) of the working tree, or
    ``None`` if git is unavailable."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def capture_provenance() -> dict[str, Any]:
    """Full run provenance so a benchmark result is interpretable and comparable:
    library versions, hardware, and the available Warp backend/devices."""
    prov: dict[str, Any] = {
        "ionmc_version": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or platform.machine(),
        "cpu_count": os.cpu_count(),
        "git_sha": _git_sha(),
    }
    if mathlib.HAVE_WARP:
        wp = mathlib.warp_module()
        quiet_warp()  # suppress the init banner before it can pollute stdout JSON
        prov["warp_version"] = str(wp.config.version)
        try:
            prov["warp_devices"] = [dv.alias for dv in wp.get_devices()]
        except Exception:  # best-effort provenance
            prov["warp_devices"] = []
    else:
        prov["warp_version"] = None
        prov["warp_devices"] = []
    return prov


def quiet_warp() -> None:
    """Silence Warp's INFO-level init banner and logging so it never pollutes a
    JSON report written to stdout. Idempotent; a no-op when Warp is absent."""
    if not mathlib.HAVE_WARP:
        return
    wp = mathlib.warp_module()
    wp.config.quiet = True
    wp.config.log_level = wp.LOG_WARNING


def backend_label(path: str, device: str | None) -> str:
    """Canonical backend name for a report, e.g. ``reference``, ``warp:cpu``,
    ``warp:cuda:0``."""
    if path != "warp":
        return "reference"
    return f"warp:{device or 'cpu'}"


def make_sync(path: str, device: str | None) -> Callable[[], None] | None:
    """A device-synchronisation callback for honest GPU timing, or ``None`` when
    the backend runs synchronously (reference / Warp CPU). Warp kernel launches on
    CUDA are asynchronous, so the timed region must wait for the device."""
    if path == "warp" and mathlib.HAVE_WARP and device and device != "cpu":
        wp = mathlib.warp_module()

        def _sync() -> None:
            wp.synchronize_device(device)

        return _sync
    return None


@dataclass
class Timing:
    """Timing statistics for one benchmarked backend."""

    label: str
    backend: str
    repeats: int
    warmup: int
    seconds_min: float
    seconds_median: float
    seconds_mean: float
    seconds_stdev: float
    work_units: float | None = None
    #: primary throughput = ``work_units / seconds_median`` (e.g. histories/s)
    throughput_per_s: float | None = None
    seconds: list[float] = field(default_factory=list)


def measure(
    fn: Callable[[], Any],
    *,
    label: str = "",
    backend: str = "",
    repeats: int = 5,
    warmup: int = 1,
    sync: Callable[[], None] | None = None,
    work_units: float | None = None,
) -> Timing:
    """Time ``fn`` over ``repeats`` runs after ``warmup`` untimed runs, returning
    the timing statistics. ``sync`` (if given) is called inside the timed region
    after ``fn`` so asynchronous backends are measured honestly. ``work_units`` (a
    problem-size count such as the number of histories) yields a throughput.
    """
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    if warmup < 0:
        raise ValueError("warmup must be >= 0")
    for _ in range(warmup):
        fn()
        if sync is not None:
            sync()
    times: list[float] = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        if sync is not None:
            sync()
        times.append(time.perf_counter() - t0)
    median = statistics.median(times)
    throughput = (
        work_units / median if (work_units is not None and median > 0.0) else None
    )
    return Timing(
        label=label,
        backend=backend,
        repeats=repeats,
        warmup=warmup,
        seconds_min=min(times),
        seconds_median=median,
        seconds_mean=statistics.fmean(times),
        seconds_stdev=statistics.pstdev(times) if len(times) > 1 else 0.0,
        work_units=work_units,
        throughput_per_s=throughput,
        seconds=times,
    )


def array_digest(a: np.ndarray, decimals: int = 6) -> str:
    """A short, order-sensitive fingerprint of a numeric array, rounded to
    ``decimals`` places in float64 so it is stable across backends to that
    precision. Two runs of the same backend with the same physics share a digest;
    a changed result changes it — the basis for the V6 'unchanged scientific outcome'
    check. For a float32 (Warp) backend the 6-decimal boundary means the digest is a
    reliable regression pin only within the same backend/hardware; cross-backend
    correctness is checked with :func:`relative_agreement`, not digest equality."""
    rounded = np.round(np.asarray(a, dtype=np.float64), decimals)
    return hashlib.sha256(np.ascontiguousarray(rounded).tobytes()).hexdigest()[:16]


def relative_agreement(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    """Cross-backend agreement between two result arrays that carry their **absolute
    magnitude** (e.g. per-history depth dose), so a uniform-scale divergence is
    caught — the correctness check that ties a benchmark's speed to the validated
    physics. Do not pass unit-normalised curves: that would hide any constant-factor
    error. Returns three relative differences, all normalised by ``sum(a)``:

    * ``total_rel_diff`` — integral difference ``|sum(b)-sum(a)|/sum(a)``;
    * ``cumulative_rel_diff`` — ``max|cumsum(a)-cumsum(b)|/sum(a)``, the edge-aware
      metric used by the established depth-dose gate (decision 0009), robust to
      sub-bin range shifts at the steep distal Bragg falloff;
    * ``max_bin_rel_diff`` — worst per-bin difference normalised by the peak bin.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    denom = float(a.sum())
    peak = float(np.max(np.abs(a)))
    total_rel = abs(float(b.sum()) - denom) / denom if denom != 0.0 else float("nan")
    cum_rel = (
        float(np.max(np.abs(np.cumsum(a) - np.cumsum(b))) / denom)
        if denom != 0.0
        else float("nan")
    )
    max_rel = float(np.max(np.abs(b - a)) / peak) if peak > 0.0 else float("nan")
    return {
        "total_rel_diff": total_rel,
        "cumulative_rel_diff": cum_rel,
        "max_bin_rel_diff": max_rel,
    }


def _jsonable(o: Any) -> Any:
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")


def timing_to_dict(t: Timing) -> dict[str, Any]:
    """A JSON-ready dict for a :class:`Timing` (including the full per-repeat
    ``seconds`` list alongside the summary statistics)."""
    return asdict(t)


def new_report(name: str) -> dict[str, Any]:
    """A benchmark report seeded with the schema version and full provenance."""
    return {
        "schema_version": 1,
        "benchmark": name,
        "provenance": capture_provenance(),
        "timings": [],
    }


def dump_report(report: dict[str, Any], stream: Any) -> None:
    """Write ``report`` as indented JSON (numpy-aware) with a trailing newline."""
    json.dump(report, stream, indent=2, default=_jsonable)
    stream.write("\n")


# -- performance-regression tracking (milestone V6, decision 0032) -------------
#
# A benchmark's *physics identity* is anchored on its reference (float64) result,
# which is deterministic and hardware-independent, so it is the stable pin across
# machines and revisions. Wall-clock throughput is machine-dependent and is compared
# only informationally. A regression check therefore GATES on the reference digest
# (the "unchanged scientific outcome" V6 requirement) and REPORTS throughput deltas.

#: reference (float64) integral must match the baseline to this relative tolerance;
#: the reference path is deterministic, so this only absorbs last-bit summation order.
BASELINE_INTEGRAL_TOL = 1.0e-9


def throughput_by_backend(report: dict[str, Any]) -> dict[str, float]:
    """Per-backend peak throughput [work/s] from a benchmark report, tolerating both
    report shapes: an explicit ``peak_throughput_per_s`` (the scaling-sweep drivers)
    or the per-backend ``timings``/``scaling`` rows (max over sizes)."""
    if isinstance(report.get("peak_throughput_per_s"), dict):
        return {k: float(v) for k, v in report["peak_throughput_per_s"].items()}
    out: dict[str, float] = {}
    for row in (*report.get("timings", []), *report.get("scaling", [])):
        tp = row.get("throughput_per_s")
        if tp is not None:
            out[row["backend"]] = max(out.get(row["backend"], 0.0), float(tp))
    return out


def physics_fingerprint(report: dict[str, Any]) -> dict[str, Any]:
    """The hardware-independent physics identity of a benchmark run: its name, the
    config that defines the workload, and the reference (float64) digest. This is
    what a regression baseline pins; two runs of the same physics share it on any
    machine."""
    return {
        "benchmark": report["benchmark"],
        "config": report.get("config", {}),
        "reference_digest": report["digest"]["reference"],
    }


def make_baseline(report: dict[str, Any]) -> dict[str, Any]:
    """Build a committable regression baseline from a benchmark report: the
    hardware-independent physics fingerprint (the gate) plus a machine-tagged
    throughput snapshot (informational)."""
    prov = report.get("provenance", {})
    return {
        "schema_version": 1,
        "physics": physics_fingerprint(report),
        "throughput_by_backend": throughput_by_backend(report),
        "recorded_on": {
            "platform": prov.get("platform"),
            "machine": prov.get("machine"),
            "warp_version": prov.get("warp_version"),
            "warp_devices": prov.get("warp_devices"),
            "ionmc_version": prov.get("ionmc_version"),
            "git_sha": prov.get("git_sha"),
        },
    }


def compare_to_baseline(
    report: dict[str, Any], baseline: dict[str, Any]
) -> dict[str, Any]:
    """Compare a fresh benchmark report against a committed baseline. GATES on the
    reference physics digest (unchanged scientific outcome) and REPORTS per-backend
    throughput ratios (current / baseline) without gating on them.

    Returns ``{"physics_ok", "benchmark_match", "reference_digest_match",
    "integral_rel_diff", "throughput"}``; ``physics_ok`` is the pass/fail signal.
    """
    cur = physics_fingerprint(report)
    base = baseline["physics"]
    benchmark_match = cur["benchmark"] == base["benchmark"]
    cur_dig = cur["reference_digest"]
    base_dig = base["reference_digest"]
    digest_match = cur_dig.get("shape_digest") == base_dig.get("shape_digest")
    base_integral = float(base_dig["integral_per_history_mev"])
    integral_rel = (
        abs(float(cur_dig["integral_per_history_mev"]) - base_integral) / base_integral
        if base_integral != 0.0
        else float("nan")
    )
    physics_ok = bool(
        benchmark_match and digest_match and integral_rel <= BASELINE_INTEGRAL_TOL
    )

    cur_tp = throughput_by_backend(report)
    base_tp = baseline.get("throughput_by_backend", {})
    throughput: dict[str, dict[str, float | None]] = {}
    for backend in sorted(set(cur_tp) | set(base_tp)):
        b = base_tp.get(backend)
        c = cur_tp.get(backend)
        ratio = (c / b) if (b and c) else None
        throughput[backend] = {"baseline": b, "current": c, "ratio": ratio}

    return {
        "physics_ok": physics_ok,
        "benchmark_match": benchmark_match,
        "reference_digest_match": digest_match,
        "integral_rel_diff": integral_rel,
        "throughput": throughput,
    }
