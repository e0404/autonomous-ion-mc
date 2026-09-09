#!/usr/bin/env python3
"""NVIDIA Warp backend runtime and cross-backend parity diagnostics.

This is the project's first cross-backend (Warp CPU vs Warp CUDA)
numerical-consistency instrument (see REQUIREMENTS.md, "Cross-backend
statistical consistency", and EXPERIMENT.md's validation list: "comparison
between Warp CPU and CUDA execution"). Its numbers are treated as
validation evidence, so the report documents everything needed to
reproduce them: the exact input generation, kernels, dtypes, devices, and
tolerances used (see the "configuration" section below).

It reports, as a single machine-readable JSON document:

* `warp`      - whether NVIDIA Warp imported at all, and its version.
* `devices`   - every Warp device Warp itself enumerates, best-effort
                described (name/arch/memory).
* `parity`    - elementwise comparison of two small deterministic kernels
                (`affine`, `transcendental`), each in `float32` and
                `float64`, run on every available device and compared
                against the Warp **CPU** device as the reference.
* `configuration` - the effective run parameters (element count, seed,
                dtypes, kernels, domain, tolerances) actually used, so a
                reported result is reproducible from the report alone.

Design goals (mirrors infrastructure/diagnostics/environment_report.py -
see also docs/warp_backend_diagnostics.md):

* `warp` and `numpy` are optional imports, used lazily only inside the
  execution path (`_import_warp()`, `WarpKernelRunner`). The module
  itself, and all of its comparison/verdict logic, must import and work
  correctly with neither installed - GitHub CI installs only pytest.
* Every probe/workload is independently guarded: one failing device,
  kernel, or dtype - or `warp` being entirely absent - must never prevent
  the rest of the report from being produced.
* The numerical comparison/verdict logic is pure Python over plain
  `list[float]` sequences (`compare_sequences`, `bits_equal`, `abs_diff`,
  `rel_diff`, `select_verdict`); `numpy`/`warp` are only ever used to move
  data into and out of Warp arrays, never inside the comparison logic
  itself, so that logic is fully unit-testable without either dependency.
* Device enumeration (`build_devices_section`) and kernel execution
  (`WarpKernelRunner`) both take an explicit module/callable parameter
  rather than reaching for a global, so a fake stand-in can be injected
  by tests through the same public parameters real callers use - no
  monkeypatching of private module internals required.
* Deliberate deviation from the `environment_report.py` convention: this
  module does NOT use `from __future__ import annotations`. Warp's
  `@wp.kernel` decorator inspects a kernel function's argument
  annotations and needs them to already be real objects (e.g.
  `wp.array(dtype=wp.float32)`) at decoration time, not the postponed
  string form PEP 563 produces. Kernels are only ever constructed inside
  `_build_kernels()`, after a real `warp` module is available, so this
  does not affect whether the module imports cleanly without `warp`.
* No environment-variable dumps, user names, tokens, or remote URLs are
  ever included.

Status vocabulary (identical to environment_report.py):

* "available"   - the probe succeeded; the section's data fields are populated.
* "unavailable" - the underlying capability is simply absent (e.g. `warp`
                  is not importable, or zero devices/no CPU device were
                  enumerated). Expected, not a failure.
* "error"       - the probe was attempted but failed unexpectedly. A
                  "detail" field explains why.

Exit code policy (see also `main()`):

* 0 - the report was produced and delivered, and no requested `--require-*`
      gate failed.
* 1 - the report was produced but could not be written to `--output`.
* 2 - argparse usage error, raised before any report is generated.
* 3 - the report was produced and delivered successfully, but a requested
      `--require-cpu` / `--require-cuda` / `--require-parity` gate was not
      satisfied. Producing a report is always success by itself; gates are
      an opt-in, separate expectation on top of that.
"""

import argparse
import json
import math
import random
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

SCHEMA_VERSION = 1

STATUS_AVAILABLE = "available"
STATUS_UNAVAILABLE = "unavailable"
STATUS_ERROR = "error"

VERDICT_IDENTICAL = "identical"
VERDICT_WITHIN_TOLERANCE = "within_tolerance"
VERDICT_EXCEEDS_TOLERANCE = "exceeds_tolerance"

KERNELS = ("affine", "transcendental")
DTYPES = ("float32", "float64")

# Deterministic input domain. Kept strictly positive so `sqrt` in the
# `transcendental` kernel is always defined, and away from 0 so `exp(-x)`
# and `sin(x)` stay well-scaled (no denormals, no argument so large that
# transcendental-implementation rounding differences vanish into zero).
DOMAIN_MIN = 0.05
DOMAIN_MAX = 4.0

# Small enough to run in well under a second on CPU or GPU, large enough
# to exercise a real range of the input domain.
DEFAULT_ELEMENTS = 4096
DEFAULT_SEED = 12345

# Per-(kernel, dtype) relative-tolerance budget applied to max_rel_diff
# when deciding "within_tolerance" vs "exceeds_tolerance" (see
# select_verdict()). Chosen so that:
#
# * `affine` uses only +, -, * on values that never approach cancellation
#   catastrophically, and deliberately avoids any single fused
#   multiply-add opportunity (`(x*x - x) + 0.5` is three separate
#   operations) - Warp CPU and CUDA are expected to be BITWISE identical,
#   so its tolerance is 0.0. A 0.0 tolerance still allows a
#   "within_tolerance" verdict (rather than "exceeds_tolerance") for the
#   edge case of two values that are numerically equal but bit-distinct,
#   e.g. -0.0 vs 0.0.
# * `transcendental` invokes sqrt/exp/sin, whose CPU and CUDA
#   implementations are different, IEEE-754-compliant but not required to
#   be identically rounded. float32 transcendental intrinsics commonly
#   differ from their CPU counterparts by a handful of ULPs (single
#   precision has ~7 decimal digits, so a few ULPs is ~1e-6..1e-7
#   relative); 1e-5 keeps a comfortable margin. float64 transcendental
#   functions typically agree far more tightly between good
#   implementations (~1e-12); 1e-9 keeps a comfortable margin without
#   being so loose it would miss a real regression.
RELATIVE_TOLERANCE: dict[tuple[str, str], float] = {
    ("affine", "float32"): 0.0,
    ("affine", "float64"): 0.0,
    ("transcendental", "float32"): 1e-5,
    ("transcendental", "float64"): 1e-9,
}


# ---------------------------------------------------------------------------
# Pure comparison helpers (no numpy/warp; operate on plain list[float])
# ---------------------------------------------------------------------------


def bits_equal(a: float, b: float) -> bool:
    """True Boolean bit-pattern equality of two Python floats.

    Unlike `a == b`, this distinguishes -0.0 from 0.0 (different bit
    patterns, but `-0.0 == 0.0` is True in IEEE 754/Python). Any NaN is
    treated as bit-equal to any other NaN (both indicate "not a valid
    result" for this diagnostic's purposes); a NaN is never bit-equal to
    a non-NaN.
    """
    a_nan = math.isnan(a)
    b_nan = math.isnan(b)
    if a_nan or b_nan:
        return a_nan and b_nan
    return struct.pack("<d", a) == struct.pack("<d", b)


def abs_diff(a: float, b: float) -> float:
    """Absolute difference, safe for inf/NaN.

    Two numerically equal values (including a == b == +-inf) have a
    difference of 0.0. Otherwise, either operand being NaN makes the
    difference NaN (undefined/incomparable), never raises.
    """
    if a == b:
        return 0.0
    if math.isnan(a) or math.isnan(b):
        return math.nan
    return abs(a - b)


def rel_diff(a: float, b: float) -> float:
    """Relative difference of `b` from reference `a`, safe for inf/NaN/zero.

    Relative to `abs(a)`. When the reference is exactly 0.0, the relative
    difference is 0.0 if the values are equal (including bit-distinct
    zeros) and otherwise `inf` (no finite relative measure exists).
    """
    diff = abs_diff(a, b)
    if math.isnan(diff):
        return math.nan
    denominator = abs(a)
    if denominator == 0.0:
        return 0.0 if diff == 0.0 else math.inf
    return diff / denominator


def _reduce_max(values: list) -> float:
    """Max of a list of floats, propagating NaN (max([], any nan) -> NaN).

    Plain `max()` on a list containing NaN gives an order-dependent,
    unreliable result because NaN compares False against everything;
    this makes "any NaN present" always win, and treats an empty list as
    a difference of 0.0.
    """
    if not values:
        return 0.0
    result = 0.0
    saw_nan = False
    for value in values:
        if math.isnan(value):
            saw_nan = True
            continue
        if value > result:
            result = value
    return math.nan if saw_nan else result


def compare_sequences(reference: list, candidate: list) -> dict:
    """Elementwise-compare `candidate` against `reference`.

    Returns a dict with `n`, `bitwise_identical`, `max_abs_diff`,
    `max_rel_diff`. Raises ValueError if the sequences have different
    lengths (a structural failure, not a numerical one - callers should
    surface it as an "error" result rather than a verdict).
    """
    if len(reference) != len(candidate):
        raise ValueError(f"length mismatch: reference has {len(reference)} elements, candidate has {len(candidate)}")

    bitwise_identical = True
    abs_diffs = []
    rel_diffs = []
    for ref_value, cand_value in zip(reference, candidate):
        if not bits_equal(ref_value, cand_value):
            bitwise_identical = False
        abs_diffs.append(abs_diff(ref_value, cand_value))
        rel_diffs.append(rel_diff(ref_value, cand_value))

    return {
        "n": len(reference),
        "bitwise_identical": bitwise_identical,
        "max_abs_diff": _reduce_max(abs_diffs),
        "max_rel_diff": _reduce_max(rel_diffs),
    }


def select_verdict(bitwise_identical: bool, max_rel_diff: float, tolerance: float) -> str:
    """Choose a verdict for one (device, kernel, dtype) comparison."""
    if bitwise_identical:
        return VERDICT_IDENTICAL
    if math.isnan(max_rel_diff):
        return VERDICT_EXCEEDS_TOLERANCE
    if max_rel_diff <= tolerance:
        return VERDICT_WITHIN_TOLERANCE
    return VERDICT_EXCEEDS_TOLERANCE


def generate_inputs(n: int, seed: int) -> list:
    """Deterministic input vector, reproducible without numpy.

    Uses `random.Random(seed)` (a private, seeded PRNG instance - not the
    module-global `random` state) so this is reproducible and does not
    disturb or depend on any other code's use of the `random` module.
    """
    rng = random.Random(seed)
    return [rng.uniform(DOMAIN_MIN, DOMAIN_MAX) for _ in range(n)]


# ---------------------------------------------------------------------------
# Optional `warp` import
# ---------------------------------------------------------------------------


def _import_warp():
    """Best-effort `import warp`. Returns (module_or_None, error_section_or_None)."""
    try:
        import warp  # type: ignore  # noqa: PLC0415
    except ImportError as exc:
        return None, {"status": STATUS_UNAVAILABLE, "detail": f"warp is not importable: {exc}"}
    except Exception as exc:  # noqa: BLE001 - probes must never raise
        return None, {"status": STATUS_ERROR, "detail": f"unexpected error importing warp: {exc}"}
    return warp, None


def build_warp_version_section(warp_module=None) -> dict:
    """The `warp` section: whether Warp imported at all, and its version."""
    if warp_module is None:
        warp_module, error = _import_warp()
        if error is not None:
            return error

    info: dict[str, Any] = {"status": STATUS_AVAILABLE, "version": None}
    try:
        info["version"] = getattr(warp_module, "__version__", None)
    except Exception as exc:  # noqa: BLE001
        info["version_error"] = str(exc)
    return info


# ---------------------------------------------------------------------------
# Device enumeration
# ---------------------------------------------------------------------------


def _get_attr_maybe_call(obj, name: str):
    """`getattr(obj, name)`, calling it if it turns out to be a method.

    Warp's Device API has used both plain attributes and zero-argument
    methods for fields like `total_memory` across versions; this
    tolerates either shape, and any failure degrades to None rather than
    raising, per the module's "unexpected Warp API shape -> null, never
    raise" guarantee.
    """
    try:
        value = getattr(obj, name, None)
    except Exception:  # noqa: BLE001
        return None
    if callable(value):
        try:
            return value()
        except Exception:  # noqa: BLE001
            return None
    return value


def _describe_device(device) -> dict:
    """Best-effort description of one Warp device. Never raises."""
    try:
        alias = getattr(device, "alias", None)
        alias = str(alias) if alias is not None else str(device)
    except Exception:  # noqa: BLE001
        alias = None

    is_cpu = _get_attr_maybe_call(device, "is_cpu")
    is_cuda = _get_attr_maybe_call(device, "is_cuda")
    name = _get_attr_maybe_call(device, "name")
    arch = _get_attr_maybe_call(device, "arch")
    total_memory = _get_attr_maybe_call(device, "total_memory")
    free_memory = _get_attr_maybe_call(device, "free_memory")

    def _bool_or_none(value):
        return bool(value) if value is not None else None

    def _int_or_none(value):
        try:
            return int(value) if value is not None else None
        except Exception:  # noqa: BLE001
            return None

    def _str_or_none(value):
        try:
            return str(value) if value is not None else None
        except Exception:  # noqa: BLE001
            return None

    return {
        "alias": alias,
        "is_cpu": _bool_or_none(is_cpu),
        "is_cuda": _bool_or_none(is_cuda),
        "name": _str_or_none(name),
        "arch": _str_or_none(arch),
        "total_memory_bytes": _int_or_none(total_memory),
        "free_memory_bytes": _int_or_none(free_memory),
    }


def _enumerate_devices(warp_module) -> list:
    """List every Warp device, tolerating older/newer Warp API shapes.

    If the preferred `get_all_devices()` API is present but raises, that
    is a genuine enumeration failure and is allowed to propagate (the
    caller reports it as "error"), rather than being silently treated as
    "API not present" and falling through to the legacy fallback.
    """
    if hasattr(warp_module, "get_all_devices"):
        return list(warp_module.get_all_devices())

    devices = []
    try:
        cpu = warp_module.get_device("cpu")
        if cpu is not None:
            devices.append(cpu)
    except Exception:  # noqa: BLE001
        pass

    try:
        count = int(warp_module.get_cuda_device_count()) if hasattr(warp_module, "get_cuda_device_count") else 0
    except Exception:  # noqa: BLE001
        count = 0

    for index in range(count):
        try:
            devices.append(warp_module.get_cuda_device(index))
        except Exception:  # noqa: BLE001
            continue

    return devices


def build_devices_section(warp_module=None) -> dict:
    """The `devices` section: every Warp device Warp itself enumerates."""
    if warp_module is None:
        warp_module, error = _import_warp()
        if error is not None:
            return {"status": error["status"], "detail": error["detail"], "devices": []}

    try:
        raw_devices = _enumerate_devices(warp_module)
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": f"failed to enumerate warp devices: {exc}", "devices": []}

    described = [_describe_device(device) for device in raw_devices]
    if not described:
        return {"status": STATUS_UNAVAILABLE, "detail": "warp reported zero devices", "devices": []}
    return {"status": STATUS_AVAILABLE, "devices": described}


# ---------------------------------------------------------------------------
# Kernels and execution (only ever touched with a real `warp` module)
# ---------------------------------------------------------------------------


def _build_kernels(wp) -> dict:
    """Build the four (kernel, dtype) Warp kernels.

    Only ever called with a real `warp` module (as `wp`), after import has
    already succeeded - never at module import time. Kernel argument
    annotations reference `wp.array`/`wp.float32`/`wp.float64` directly
    (not a closed-over dtype variable) so they resolve to real Warp types
    immediately at function-definition time, which is what `@wp.kernel`
    requires (see the module docstring for why this file avoids
    `from __future__ import annotations`).
    """

    @wp.kernel
    def affine_f32(x: wp.array(dtype=wp.float32), y: wp.array(dtype=wp.float32)):
        i = wp.tid()
        xi = x[i]
        y[i] = (xi * xi - xi) + wp.float32(0.5)

    @wp.kernel
    def affine_f64(x: wp.array(dtype=wp.float64), y: wp.array(dtype=wp.float64)):
        i = wp.tid()
        xi = x[i]
        y[i] = (xi * xi - xi) + wp.float64(0.5)

    @wp.kernel
    def transcendental_f32(x: wp.array(dtype=wp.float32), y: wp.array(dtype=wp.float32)):
        i = wp.tid()
        xi = x[i]
        y[i] = wp.sqrt(xi) * wp.exp(-xi) + wp.sin(xi)

    @wp.kernel
    def transcendental_f64(x: wp.array(dtype=wp.float64), y: wp.array(dtype=wp.float64)):
        i = wp.tid()
        xi = x[i]
        y[i] = wp.sqrt(xi) * wp.exp(-xi) + wp.sin(xi)

    return {
        ("affine", "float32"): affine_f32,
        ("affine", "float64"): affine_f64,
        ("transcendental", "float32"): transcendental_f32,
        ("transcendental", "float64"): transcendental_f64,
    }


class WarpKernelRunner:
    """Runs a (kernel, dtype) pair on a given Warp device.

    Callable as `runner(kernel_name, dtype_name, device_alias, values) ->
    list[float]`. This is the seam `build_parity_section()` depends on;
    tests can substitute any callable with that signature (e.g. a plain
    Python arithmetic stand-in) instead of this class, without needing
    `warp` installed.
    """

    def __init__(self, warp_module):
        self._wp = warp_module
        self._kernels: dict | None = None

    def _kernels_dict(self) -> dict:
        if self._kernels is None:
            self._kernels = _build_kernels(self._wp)
        return self._kernels

    def __call__(self, kernel_name: str, dtype_name: str, device_alias: str, values: list) -> list:
        import numpy as np  # noqa: PLC0415

        wp = self._wp
        kernel = self._kernels_dict()[(kernel_name, dtype_name)]
        wp_dtype = getattr(wp, dtype_name)
        np_dtype = getattr(np, dtype_name)

        x_host = np.array(values, dtype=np_dtype)
        x = wp.array(x_host, dtype=wp_dtype, device=device_alias)
        y = wp.zeros(len(values), dtype=wp_dtype, device=device_alias)
        wp.launch(kernel, dim=len(values), inputs=[x, y], device=device_alias)
        if hasattr(wp, "synchronize_device"):
            wp.synchronize_device(device_alias)
        else:
            wp.synchronize()
        return y.numpy().tolist()


# ---------------------------------------------------------------------------
# Parity section
# ---------------------------------------------------------------------------


def build_parity_section(
    devices: list,
    run_kernel: Callable[[str, str, str, list], list],
    elements: int = DEFAULT_ELEMENTS,
    seed: int = DEFAULT_SEED,
) -> dict:
    """The `parity` section: cross-backend comparison against the CPU device.

    `devices` is the list of device dicts produced by
    `build_devices_section()` (or an equivalent fake for testing).
    `run_kernel` executes one (kernel, dtype) pair on one device alias -
    normally a `WarpKernelRunner`, but any compatible callable can be
    injected.
    """
    cpu_alias = next((device["alias"] for device in devices if device.get("is_cpu") and device.get("alias")), None)
    if cpu_alias is None:
        return {
            "status": STATUS_UNAVAILABLE,
            "detail": "no warp CPU device available; CPU is required as the parity reference",
            "reference_device": None,
            "results": [],
        }

    values = generate_inputs(elements, seed)
    results = []
    any_available = False

    for kernel_name in KERNELS:
        for dtype_name in DTYPES:
            tolerance = RELATIVE_TOLERANCE[(kernel_name, dtype_name)]
            try:
                reference = run_kernel(kernel_name, dtype_name, cpu_alias, values)
            except Exception as exc:  # noqa: BLE001
                detail = f"failed to compute CPU reference: {exc}"
                for device in devices:
                    results.append(
                        {
                            "device": device.get("alias"),
                            "kernel": kernel_name,
                            "dtype": dtype_name,
                            "status": STATUS_ERROR,
                            "detail": detail,
                        }
                    )
                continue

            for device in devices:
                alias = device.get("alias")
                try:
                    candidate = run_kernel(kernel_name, dtype_name, alias, values)
                    stats = compare_sequences(reference, candidate)
                    verdict = select_verdict(stats["bitwise_identical"], stats["max_rel_diff"], tolerance)
                    results.append(
                        {
                            "device": alias,
                            "kernel": kernel_name,
                            "dtype": dtype_name,
                            "status": STATUS_AVAILABLE,
                            "n": stats["n"],
                            "max_abs_diff": stats["max_abs_diff"],
                            "max_rel_diff": stats["max_rel_diff"],
                            "bitwise_identical": stats["bitwise_identical"],
                            "verdict": verdict,
                            "tolerance": tolerance,
                        }
                    )
                    any_available = True
                except Exception as exc:  # noqa: BLE001
                    results.append(
                        {
                            "device": alias,
                            "kernel": kernel_name,
                            "dtype": dtype_name,
                            "status": STATUS_ERROR,
                            "detail": str(exc),
                        }
                    )

    return {
        "status": STATUS_AVAILABLE if any_available else STATUS_ERROR,
        "reference_device": cpu_alias,
        "results": results,
    }


# ---------------------------------------------------------------------------
# Configuration section
# ---------------------------------------------------------------------------


def build_configuration_section(elements: int, seed: int) -> dict:
    """The `configuration` section: effective run parameters actually used.

    Always `available` - it documents this run's own parameters rather
    than probing anything that can fail.
    """
    return {
        "status": STATUS_AVAILABLE,
        "elements": elements,
        "seed": seed,
        "kernels": list(KERNELS),
        "dtypes": list(DTYPES),
        "reference_device_role": "cpu",
        "domain_min": DOMAIN_MIN,
        "domain_max": DOMAIN_MAX,
        "relative_tolerances": {f"{kernel}/{dtype}": tol for (kernel, dtype), tol in RELATIVE_TOLERANCE.items()},
    }


# ---------------------------------------------------------------------------
# Report assembly
# ---------------------------------------------------------------------------


def _safe_section(builder: Callable[[], dict]) -> dict:
    """Outer safety net around a section builder, mirroring environment_report.py."""
    try:
        result = builder()
    except Exception as exc:  # noqa: BLE001
        return {"status": STATUS_ERROR, "detail": f"unexpected error: {exc}"}

    if not isinstance(result, dict) or "status" not in result:
        return {"status": STATUS_ERROR, "detail": "probe returned a malformed result"}

    return result


def build_report(
    elements: int = DEFAULT_ELEMENTS,
    seed: int = DEFAULT_SEED,
    warp_module=None,
    run_kernel: Callable[[str, str, str, list], list] | None = None,
) -> dict:
    """Build the full Warp backend diagnostics report.

    `warp_module` and `run_kernel` are injection seams for tests: pass a
    fake module implementing enough of the Warp API for
    `build_devices_section`, and/or a fake callable implementing the
    `run_kernel` signature, to exercise report assembly without `warp`
    installed.
    """
    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }

    resolved_module = warp_module
    import_error: dict | None = None
    if resolved_module is None:
        resolved_module, import_error = _import_warp()

    report["warp"] = _safe_section(
        lambda: import_error if import_error is not None else build_warp_version_section(resolved_module)
    )
    report["devices"] = _safe_section(
        lambda: import_error if import_error is not None else build_devices_section(resolved_module)
    )

    devices_section = report["devices"]
    device_list = devices_section.get("devices", []) if devices_section.get("status") == STATUS_AVAILABLE else []

    if not device_list:
        detail = (
            import_error["detail"]
            if import_error is not None
            else devices_section.get("detail", "no warp devices available for parity comparison")
        )
        report["parity"] = {
            "status": STATUS_UNAVAILABLE,
            "detail": detail,
            "reference_device": None,
            "results": [],
        }
    else:
        kernel_runner = run_kernel if run_kernel is not None else WarpKernelRunner(resolved_module)
        report["parity"] = _safe_section(
            lambda: build_parity_section(device_list, kernel_runner, elements=elements, seed=seed)
        )

    report["configuration"] = build_configuration_section(elements, seed)
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="warp_backend_report",
        description="Report the NVIDIA Warp backend runtime and cross-backend numerical parity as JSON.",
    )
    parser.add_argument("--elements", type=int, default=DEFAULT_ELEMENTS, help=f"number of input elements (default: {DEFAULT_ELEMENTS}).")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help=f"deterministic input seed (default: {DEFAULT_SEED}).")
    parser.add_argument("--indent", type=int, default=2, help="JSON indentation level (default: 2).")
    parser.add_argument("--output", type=Path, default=None, help="Write the report to this path instead of stdout.")
    parser.add_argument(
        "--require-cpu",
        action="store_true",
        help="Exit non-zero if no Warp CPU device is available.",
    )
    parser.add_argument(
        "--require-cuda",
        action="store_true",
        help="Exit non-zero if no Warp CUDA device is available.",
    )
    parser.add_argument(
        "--require-parity",
        action="store_true",
        help="Exit non-zero if the parity section is not available, or any result errored or exceeded tolerance.",
    )
    return parser


def _gate_failures(report: dict, require_cpu: bool, require_cuda: bool, require_parity: bool) -> list:
    """Evaluate the requested --require-* gates against a produced report."""
    failures = []

    devices_section = report.get("devices", {})
    device_list = devices_section.get("devices", []) if devices_section.get("status") == STATUS_AVAILABLE else []

    if require_cpu and not any(device.get("is_cpu") for device in device_list):
        failures.append("no warp CPU device available")

    if require_cuda and not any(device.get("is_cuda") for device in device_list):
        failures.append("no warp CUDA device available")

    if require_parity:
        parity = report.get("parity", {})
        if parity.get("status") != STATUS_AVAILABLE:
            failures.append(f"parity section status is {parity.get('status')!r}, not {STATUS_AVAILABLE!r}")
        else:
            bad_results = [
                result
                for result in parity.get("results", [])
                if result.get("status") != STATUS_AVAILABLE or result.get("verdict") == VERDICT_EXCEEDS_TOLERANCE
            ]
            if bad_results:
                failures.append(f"{len(bad_results)} parity result(s) errored or exceeded tolerance")

    return failures


def main(argv: list | None = None) -> int:
    """CLI entry point.

    Exit code policy - see the module docstring for the full rationale:

    * 0 - report produced and delivered; no requested gate failed.
    * 1 - report produced but could not be written to --output.
    * 2 - argparse usage error, raised before any report is generated.
    * 3 - report produced and delivered, but a requested --require-* gate failed.
    """
    args = build_arg_parser().parse_args(argv)

    report = build_report(elements=args.elements, seed=args.seed)
    text = json.dumps(report, indent=args.indent, sort_keys=True) + "\n"

    if args.output is not None:
        try:
            args.output.write_text(text, encoding="utf-8")
        except OSError as exc:
            error_payload = {"status": STATUS_ERROR, "detail": f"failed to write report to {args.output}: {exc}"}
            sys.stderr.write(json.dumps(error_payload) + "\n")
            return 1
    else:
        sys.stdout.write(text)

    failures = _gate_failures(report, args.require_cpu, args.require_cuda, args.require_parity)
    if failures:
        error_payload = {"status": STATUS_ERROR, "detail": "; ".join(failures)}
        sys.stderr.write(json.dumps(error_payload) + "\n")
        return 3

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
