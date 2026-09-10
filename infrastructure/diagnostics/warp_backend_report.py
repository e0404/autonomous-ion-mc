#!/usr/bin/env python3
"""NVIDIA Warp backend runtime and cross-backend parity diagnostics.

This is the project's first cross-backend (Warp CPU vs Warp CUDA)
numerical-consistency instrument (see REQUIREMENTS.md, "Cross-backend
statistical consistency", and EXPERIMENT.md's validation list: "comparison
between Warp CPU and CUDA execution"). Its numbers are treated as
validation evidence, so the report documents everything needed to
reproduce them: the exact input generation, workloads, dtypes, devices, and
tolerances used (see the "configuration" section below).

It reports, as a single machine-readable JSON document:

* `warp`      - whether NVIDIA Warp imported at all, and its version.
* `devices`   - every Warp device Warp itself enumerates, best-effort
                described (name/arch/memory).
* `parity`    - elementwise comparison of three small deterministic
                workloads, each in `float32` and `float64`, run on every
                available device and compared against the Warp **CPU**
                device as the reference.
* `configuration` - the effective run parameters (element count, seed,
                workloads, dtypes, domain, tolerances) actually used, so a
                reported result is reproducible from the report alone.

Design goals (mirrors infrastructure/diagnostics/environment_report.py -
see also docs/warp_backend_diagnostics.md):

* `warp` and `numpy` are optional imports, used lazily only inside the
  execution path (`_import_warp()`, `WarpKernelRunner`). The module
  itself, and all of its comparison/verdict logic, must import and work
  correctly with neither installed - GitHub CI installs only pytest.
* Every probe/workload is independently guarded: one failing device,
  workload, or dtype - or `warp` being entirely absent - must never
  prevent the rest of the report from being produced.
* The numerical comparison/verdict logic is pure Python over plain
  `list[float]` sequences (`compare_sequences`, `bits_equal`, `abs_diff`,
  `rel_diff`, `normalized_diff`, `select_verdict`); `numpy`/`warp` are only
  ever used to move data into and out of Warp arrays, never inside the
  comparison logic itself, so that logic is fully unit-testable without
  either dependency.
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

Workloads and expectation classes
----------------------------------

Round 1 of this diagnostic used only two kernels, `affine` and
`transcendental`, and claimed `affine` was free of any fused-multiply-add
(FMA) contraction opportunity. That claim was FALSE and was falsified by a
real GPU run (see the module's `RTOL`/`ATOL` comment below for the exact
measured numbers): `x*x - x` IS a multiply immediately followed by a
subtract of the same operand, i.e. exactly `fma(x, x, -x)`, and the CPU and
CUDA C++ compilers backing Warp contract it differently, producing a ~1 ULP
difference that failed the round-1 0.0 tolerance. Rather than removing
that expression, it is now its own workload with an honestly-named
"contraction_sensitive" expectation class, because the contraction
difference it exposes is real, useful evidence, not a bug to hide:

* `product` (`bitwise` class): `y[i] = x[i] * x[i]` - a single multiply,
  one rounding, nothing adjacent to contract with and nothing to
  reassociate. This is the actual bitwise control.
* `fma_exposed` (`contraction_sensitive` class): `y[i] = (x[i]*x[i] -
  x[i]) + 0.5` - deliberately exposes an FMA contraction opportunity so a
  toolchain difference is measured and visible instead of silently
  passing or silently failing an unreachable bitwise expectation.
* `transcendental` (`transcendental` class): `y[i] = sqrt(x[i]) *
  exp(-x[i]) + sin(x[i])` - CPU and CUDA use different, non-identically-
  rounded transcendental implementations.

Verdict semantics
------------------

* `bitwise` class: verdict is `identical` iff bit-for-bit identical AND no
  NaN is present anywhere in either sequence; otherwise `exceeds_tolerance`.
  This class never returns `within_tolerance` - a numerically-equal but
  bit-distinct pair (e.g. -0.0 vs 0.0) is a genuine failure for a class
  whose entire point is bitwise reproducibility.
* `contraction_sensitive` / `transcendental` classes: `identical` when
  bitwise identical and NaN-free; `within_tolerance` when NaN-free and
  every element satisfies the mixed absolute+relative criterion (see
  `normalized_diff()`); otherwise `exceeds_tolerance`.
* NaN is never a pass, in any class. Any NaN in the reference or the
  candidate output forces `exceeds_tolerance`, regardless of bit patterns.
  `bits_equal()` itself still treats any NaN as bit-equal to any other NaN
  (it is a strict, self-contained bit-pattern predicate, useful on its
  own) - but `select_verdict()` checks `nan_count` FIRST, before looking
  at `bitwise_identical` at all, so a NaN can never produce `identical`.
  Each result reports `nan_count` so this is auditable.

Exit code policy (see also `main()`):

* 0 - the report was produced and delivered, and no requested `--require-*`
      gate failed.
* 1 - the report was produced but could not be written to `--output`.
* 2 - argparse usage error (including `--elements` < 1), raised before any
      report is generated.
* 3 - the report was produced and delivered successfully, but a requested
      `--require-cpu` / `--require-cuda` / `--require-parity` gate was not
      satisfied. Producing a report is always success by itself; gates are
      an opt-in, separate expectation on top of that.
"""

import argparse
import contextlib
import json
import math
import os
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

EXPECTATION_BITWISE = "bitwise"
EXPECTATION_CONTRACTION_SENSITIVE = "contraction_sensitive"
EXPECTATION_TRANSCENDENTAL = "transcendental"

WORKLOADS = ("product", "fma_exposed", "transcendental")
DTYPES = ("float32", "float64")

WORKLOAD_EXPECTATION_CLASS: dict[str, str] = {
    "product": EXPECTATION_BITWISE,
    "fma_exposed": EXPECTATION_CONTRACTION_SENSITIVE,
    "transcendental": EXPECTATION_TRANSCENDENTAL,
}

# Deterministic input domain. Kept strictly positive so `sqrt` in the
# `transcendental` workload is always defined, and away from 0 so
# `exp(-x)` and `sin(x)` stay well-scaled (no denormals). Note that
# `sqrt(x)*exp(-x) + sin(x)` still crosses zero near x ~= 3.2 within this
# domain - see the `max_rel_diff` caveat below.
DOMAIN_MIN = 0.05
DOMAIN_MAX = 4.0

# Small enough to run in well under a second on CPU or GPU, large enough
# to exercise a real range of the input domain.
DEFAULT_ELEMENTS = 4096
DEFAULT_SEED = 12345

# Per-(workload, dtype) mixed absolute+relative tolerance budget, as
# (rtol, atol) pairs, applied elementwise as:
#
#   |candidate - reference| <= atol + rtol * |reference|
#
# (see `normalized_diff()`). `product` is `None`/not applicable: it is the
# `bitwise` class, verdicts for which never consult a numeric budget.
#
# These were validated against a real GPU run (DRYRUN-002, NVIDIA RTX
# A6000, sm_86, CUDA 12.9, driver 13.1, Warp 1.17.0, 4096 elements, seed
# 12345), comparing Warp CUDA against Warp CPU:
#
#   fma_exposed/float32:      max_abs_diff 9.54e-07, max_rel_diff 1.19e-07
#   fma_exposed/float64:      max_abs_diff 1.78e-15, max_rel_diff 2.22e-16
#   transcendental/float32:   max_abs_diff 1.19e-07, max_rel_diff 7.62e-06
#   transcendental/float64:   max_abs_diff 2.22e-16, max_rel_diff 3.38e-14
#   product/{float32,float64} on CUDA and all four workloads' CPU-vs-CPU
#   self-comparisons: bitwise identical, diffs exactly 0.0.
#
# The `transcendental/float32` *relative* figure (7.62e-6) is inflated:
# `sqrt(x)*exp(-x)+sin(x)` crosses zero near x ~= 3.2 within [0.05, 4.0],
# and a reference value near zero turns an ordinary ~1 ULP absolute
# difference into a large relative one - this is exactly why the mixed
# absolute+relative criterion below (rather than a pure relative
# tolerance) is used for the pass/fail verdict; `max_rel_diff` remains
# informational only (see its docstring note).
#
# Budgets below are chosen with a comfortable margin over the measured
# figures: `fma_exposed` needs only to absorb a single FMA-contraction
# rounding step (~1 ULP), so float32 (1e-6, 1e-6) and float64
# (1e-12, 1e-14) sit roughly an order of magnitude above the measured
# diffs. `transcendental` uses the same budgets: its measured absolute
# diffs are of the same ~1 ULP order as `fma_exposed`'s (the relative
# figure is the misleading one, per above), so the same margin applies.
TOLERANCE: dict[tuple[str, str], tuple[float, float] | None] = {
    ("product", "float32"): None,
    ("product", "float64"): None,
    ("fma_exposed", "float32"): (1e-6, 1e-6),
    ("fma_exposed", "float64"): (1e-12, 1e-14),
    ("transcendental", "float32"): (1e-6, 1e-6),
    ("transcendental", "float64"): (1e-12, 1e-14),
}


# ---------------------------------------------------------------------------
# Pure comparison helpers (no numpy/warp; operate on plain list[float])
# ---------------------------------------------------------------------------


def bits_equal(a: float, b: float) -> bool:
    """True Boolean bit-pattern equality of two Python floats.

    Unlike `a == b`, this distinguishes -0.0 from 0.0 (different bit
    patterns, but `-0.0 == 0.0` is True in IEEE 754/Python). Any NaN is
    treated as bit-equal to any other NaN, and never bit-equal to a
    non-NaN.

    This is a strict, self-contained bit-pattern predicate - it does NOT,
    by itself, decide whether a NaN counts as a "pass" for a parity
    verdict. `select_verdict()` checks `nan_count` before ever consulting
    `bitwise_identical`, specifically so that NaN-vs-NaN cannot slip
    through as `identical` there, even though this function reports it as
    bit-equal.
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
    """Relative difference of `b` from reference `a` - INFORMATIONAL ONLY.

    Relative to `abs(a)`. When the reference is exactly 0.0, the relative
    difference is 0.0 if the values are equal (including bit-distinct
    zeros) and otherwise `inf`.

    This metric is reported for context but is NEVER used to decide a
    verdict (see `normalized_diff()` for the verdict-relevant metric): a
    reference value near zero can make an ordinary, small absolute
    difference look like a huge relative one. Concretely, the
    `transcendental` workload's expression crosses zero near x ~= 3.2
    within its input domain, and a real GPU run measured
    `max_rel_diff` = 7.62e-6 for `transcendental/float32` from nothing
    more than a ~1 ULP absolute difference landing near that crossing.
    """
    diff = abs_diff(a, b)
    if math.isnan(diff):
        return math.nan
    denominator = abs(a)
    if denominator == 0.0:
        return 0.0 if diff == 0.0 else math.inf
    return diff / denominator


def normalized_diff(a: float, b: float, rtol: float, atol: float) -> float:
    """Verdict-relevant metric: the mixed absolute+relative criterion as a
    single number that passes iff the result is <= 1.0.

    Defined as `|b - a| / (atol + rtol * |a|)`. If the budget
    `atol + rtol * |a|` is exactly 0.0 (only possible if both `atol` and
    `rtol` are 0, or `rtol` is 0 and `a` is 0), the result is 0.0 when the
    values are equal and `inf` otherwise - there is no tolerance budget to
    normalize against, so anything but exact equality is an infinite
    violation of it.
    """
    diff = abs_diff(a, b)
    if math.isnan(diff):
        return math.nan
    budget = atol + rtol * abs(a)
    if budget == 0.0:
        return 0.0 if diff == 0.0 else math.inf
    return diff / budget


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


def compare_sequences(reference: list, candidate: list, rtol: float | None = None, atol: float | None = None) -> dict:
    """Elementwise-compare `candidate` against `reference`.

    `rtol`/`atol` are the mixed-criterion budget (see `normalized_diff()`);
    pass `None` for the `bitwise` expectation class, where no numeric
    budget applies.

    Returns a dict with `n`, `bitwise_identical`, `nan_count`,
    `max_abs_diff`, `max_rel_diff` (informational only - see `rel_diff()`),
    and `max_normalized_diff` (`None` when `rtol`/`atol` are `None`).
    Raises ValueError if the sequences have different lengths, or if
    either is empty (a structural failure, not a numerical one - callers
    should surface it as an "error" result rather than a verdict; an
    empty comparison must never report a vacuous `identical`/
    `within_tolerance`).
    """
    if len(reference) != len(candidate):
        raise ValueError(f"length mismatch: reference has {len(reference)} elements, candidate has {len(candidate)}")
    if len(reference) == 0:
        raise ValueError("cannot compare empty sequences (would be a vacuous pass)")

    bitwise_identical = True
    nan_count = 0
    abs_diffs = []
    rel_diffs = []
    normalized_diffs = [] if rtol is not None and atol is not None else None

    for ref_value, cand_value in zip(reference, candidate):
        if not bits_equal(ref_value, cand_value):
            bitwise_identical = False
        if math.isnan(ref_value) or math.isnan(cand_value):
            nan_count += 1
        abs_diffs.append(abs_diff(ref_value, cand_value))
        rel_diffs.append(rel_diff(ref_value, cand_value))
        if normalized_diffs is not None:
            normalized_diffs.append(normalized_diff(ref_value, cand_value, rtol, atol))

    return {
        "n": len(reference),
        "bitwise_identical": bitwise_identical,
        "nan_count": nan_count,
        "max_abs_diff": _reduce_max(abs_diffs),
        "max_rel_diff": _reduce_max(rel_diffs),
        "max_normalized_diff": _reduce_max(normalized_diffs) if normalized_diffs is not None else None,
    }


def select_verdict(expectation_class: str, bitwise_identical: bool, nan_count: int, max_normalized_diff: float | None) -> str:
    """Choose a verdict for one (device, workload, dtype) comparison.

    NaN is checked FIRST, before bitwise identity: any NaN anywhere in
    either sequence forces `exceeds_tolerance`, regardless of class or bit
    patterns (see the module docstring's "Verdict semantics" section).
    """
    if nan_count:
        return VERDICT_EXCEEDS_TOLERANCE
    if bitwise_identical:
        return VERDICT_IDENTICAL
    if expectation_class == EXPECTATION_BITWISE:
        # The bitwise class never returns "within_tolerance": a
        # numerically-equal-but-bit-distinct pair is a genuine failure
        # for a class whose entire point is bitwise reproducibility.
        return VERDICT_EXCEEDS_TOLERANCE
    if max_normalized_diff is not None and not math.isnan(max_normalized_diff) and max_normalized_diff <= 1.0:
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


@contextlib.contextmanager
def _silence_native_output():
    """Redirect file descriptors 1 and 2 to /dev/null for the duration.

    `warp.init()` prints an initialisation banner to stdout and, on hosts
    without a CUDA driver, error text to stderr from native code. Both
    would corrupt the JSON this tool emits on those streams, so the first
    Warp call is made with the descriptors redirected. Python-level
    `sys.stdout`/`sys.stderr` are flushed first so buffered report text is
    not lost.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        for fd in (devnull, *saved):
            os.close(fd)


def _quiet_warp_init(warp_module) -> None:
    """Initialise Warp without letting it write to stdout/stderr.

    Warp >= 1.17 suppresses its banner and module-load notices when
    `config.log_level` is raised to `LOG_WARNING`; the older `config.quiet`
    flag is deprecated there and itself prints a warning to stderr, so it is
    only used as a fallback when `log_level` does not exist.
    """
    try:
        if hasattr(warp_module.config, "log_level"):
            warp_module.config.log_level = warp_module.LOG_WARNING
        else:
            warp_module.config.quiet = True
    except Exception:  # noqa: BLE001
        pass
    with _silence_native_output():
        warp_module.init()


def _import_warp():
    """Best-effort `import warp`. Returns (module_or_None, error_section_or_None)."""
    try:
        import warp  # type: ignore  # noqa: PLC0415
    except ImportError as exc:
        return None, {"status": STATUS_UNAVAILABLE, "detail": f"warp is not importable: {exc}"}
    except Exception as exc:  # noqa: BLE001 - probes must never raise
        return None, {"status": STATUS_ERROR, "detail": f"unexpected error importing warp: {exc}"}
    try:
        _quiet_warp_init(warp)
    except Exception as exc:  # noqa: BLE001 - probes must never raise
        return None, {"status": STATUS_ERROR, "detail": f"warp.init() failed: {exc}"}
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
    """Best-effort description of one Warp device. Never raises.

    Two normalizations are applied to avoid reporting a value that reads
    like real hardware metadata but isn't:

    * `arch` is reported as `null` for any non-CUDA device. A real
      DRYRUN-002 GPU-host run reported the Warp CPU device's `arch` as the
      string `"0"`, which is not an architecture identifier and would
      mislead a reader into thinking it was one.
    * A `total_memory`/`free_memory` value of exactly 0 is normalized to
      `null` rather than reported as a real zero-byte quantity. The same
      run reported the CPU device's memory as 0/0 bytes; Warp reports 0
      for these fields when the underlying quantity is unavailable (on
      that host, CPU memory query requires `psutil`, which was not
      installed in the host-runner environment and also emitted a
      warning on stderr about it) - `null` honestly states "unknown"
      instead of asserting a real 0-byte device.
    """
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

    is_cuda_bool = _bool_or_none(is_cuda)
    arch_value = _str_or_none(arch) if is_cuda_bool else None

    total_memory_bytes = _int_or_none(total_memory)
    if total_memory_bytes == 0:
        total_memory_bytes = None
    free_memory_bytes = _int_or_none(free_memory)
    if free_memory_bytes == 0:
        free_memory_bytes = None

    return {
        "alias": alias,
        "is_cpu": _bool_or_none(is_cpu),
        "is_cuda": is_cuda_bool,
        "name": _str_or_none(name),
        "arch": arch_value,
        "total_memory_bytes": total_memory_bytes,
        "free_memory_bytes": free_memory_bytes,
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
    """Build the six (workload, dtype) Warp kernels.

    Only ever called with a real `warp` module (as `wp`), after import has
    already succeeded - never at module import time. Kernel argument
    annotations reference `wp.array`/`wp.float32`/`wp.float64` directly
    (not a closed-over dtype variable) so they resolve to real Warp types
    immediately at function-definition time, which is what `@wp.kernel`
    requires (see the module docstring for why this file avoids
    `from __future__ import annotations`).
    """

    @wp.kernel
    def product_f32(x: wp.array(dtype=wp.float32), y: wp.array(dtype=wp.float32)):
        i = wp.tid()
        xi = x[i]
        y[i] = xi * xi

    @wp.kernel
    def product_f64(x: wp.array(dtype=wp.float64), y: wp.array(dtype=wp.float64)):
        i = wp.tid()
        xi = x[i]
        y[i] = xi * xi

    @wp.kernel
    def fma_exposed_f32(x: wp.array(dtype=wp.float32), y: wp.array(dtype=wp.float32)):
        i = wp.tid()
        xi = x[i]
        y[i] = (xi * xi - xi) + wp.float32(0.5)

    @wp.kernel
    def fma_exposed_f64(x: wp.array(dtype=wp.float64), y: wp.array(dtype=wp.float64)):
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
        ("product", "float32"): product_f32,
        ("product", "float64"): product_f64,
        ("fma_exposed", "float32"): fma_exposed_f32,
        ("fma_exposed", "float64"): fma_exposed_f64,
        ("transcendental", "float32"): transcendental_f32,
        ("transcendental", "float64"): transcendental_f64,
    }


class WarpKernelRunner:
    """Runs a (workload, dtype) pair on a given Warp device.

    Callable as `runner(workload_name, dtype_name, device_alias, values) ->
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

    def __call__(self, workload_name: str, dtype_name: str, device_alias: str, values: list) -> list:
        import numpy as np  # noqa: PLC0415

        wp = self._wp
        kernel = self._kernels_dict()[(workload_name, dtype_name)]
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
    `run_kernel` executes one (workload, dtype) pair on one device alias -
    normally a `WarpKernelRunner`, but any compatible callable can be
    injected.

    `elements` must be >= 1: an empty input vector would make a `bitwise`
    verdict of `identical` vacuously true (there would be nothing to
    disagree on), so this is rejected as `status: "error"` rather than
    silently producing a hollow pass. (The CLI additionally rejects
    `--elements < 1` at the argparse level - see `main()`.)
    """
    if elements < 1:
        return {
            "status": STATUS_ERROR,
            "detail": f"elements must be >= 1, got {elements} (an empty comparison would be a vacuous pass)",
            "reference_device": None,
            "results": [],
        }

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

    for workload_name in WORKLOADS:
        expectation_class = WORKLOAD_EXPECTATION_CLASS[workload_name]
        for dtype_name in DTYPES:
            budget = TOLERANCE[(workload_name, dtype_name)]
            rtol, atol = budget if budget is not None else (None, None)
            try:
                reference = run_kernel(workload_name, dtype_name, cpu_alias, values)
            except Exception as exc:  # noqa: BLE001
                detail = f"failed to compute CPU reference: {exc}"
                for device in devices:
                    results.append(
                        {
                            "device": device.get("alias"),
                            "workload": workload_name,
                            "expectation_class": expectation_class,
                            "dtype": dtype_name,
                            "status": STATUS_ERROR,
                            "detail": detail,
                        }
                    )
                continue

            for device in devices:
                alias = device.get("alias")
                try:
                    candidate = run_kernel(workload_name, dtype_name, alias, values)
                    stats = compare_sequences(reference, candidate, rtol=rtol, atol=atol)
                    verdict = select_verdict(
                        expectation_class,
                        stats["bitwise_identical"],
                        stats["nan_count"],
                        stats["max_normalized_diff"],
                    )
                    results.append(
                        {
                            "device": alias,
                            "workload": workload_name,
                            "expectation_class": expectation_class,
                            "dtype": dtype_name,
                            "status": STATUS_AVAILABLE,
                            "n": stats["n"],
                            "nan_count": stats["nan_count"],
                            "max_abs_diff": stats["max_abs_diff"],
                            "max_rel_diff": stats["max_rel_diff"],
                            "max_normalized_diff": stats["max_normalized_diff"],
                            "bitwise_identical": stats["bitwise_identical"],
                            "verdict": verdict,
                            "rtol": rtol,
                            "atol": atol,
                        }
                    )
                    any_available = True
                except Exception as exc:  # noqa: BLE001
                    results.append(
                        {
                            "device": alias,
                            "workload": workload_name,
                            "expectation_class": expectation_class,
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
        "workloads": list(WORKLOADS),
        "dtypes": list(DTYPES),
        "expectation_classes": dict(WORKLOAD_EXPECTATION_CLASS),
        "reference_device_role": "cpu",
        "domain_min": DOMAIN_MIN,
        "domain_max": DOMAIN_MAX,
        "tolerances": {
            f"{workload}/{dtype}": ({"rtol": budget[0], "atol": budget[1]} if budget is not None else {"rtol": None, "atol": None})
            for (workload, dtype), budget in TOLERANCE.items()
        },
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


def _positive_int(raw_value: str) -> int:
    """argparse `type=` for `--elements`: must parse as an int >= 1.

    An empty (or otherwise non-positive) element count is rejected here,
    at the CLI boundary, as an argparse usage error (exit code 2) - the
    same "vacuous pass" concern `build_parity_section()` guards at the
    library level for any caller that bypasses the CLI.
    """
    try:
        value = int(raw_value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid int value: {raw_value!r}") from exc
    if value < 1:
        raise argparse.ArgumentTypeError(f"--elements must be >= 1, got {value}")
    return value


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="warp_backend_report",
        description="Report the NVIDIA Warp backend runtime and cross-backend numerical parity as JSON.",
    )
    parser.add_argument(
        "--elements",
        type=_positive_int,
        default=DEFAULT_ELEMENTS,
        help=f"number of input elements, must be >= 1 (default: {DEFAULT_ELEMENTS}).",
    )
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
        help=(
            "Exit non-zero unless: the parity section is available; no result errored; at least one "
            "available comparison exists against a device other than the reference device (cross-backend "
            "evidence, not just CPU-vs-itself); every bitwise-class result is identical; and every "
            "tolerance-class result is identical or within_tolerance."
        ),
    )
    return parser


def _gate_failures(report: dict, require_cpu: bool, require_cuda: bool, require_parity: bool) -> list:
    """Evaluate the requested --require-* gates against a produced report.

    `--require-parity` fails unless ALL of the following hold:

    1. the parity section's status is "available";
    2. no result has status "error";
    3. at least one "available" result compares a device OTHER than the
       reference device - a CPU-only run (comparing the CPU reference to
       itself) is not cross-backend evidence, and must not satisfy this
       gate on its own;
    4. every `bitwise`-class result has verdict `identical`;
    5. every non-bitwise-class result has verdict `identical` or
       `within_tolerance`.
    """
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
            results = parity.get("results", [])
            reference_device = parity.get("reference_device")

            errored = [result for result in results if result.get("status") != STATUS_AVAILABLE]
            if errored:
                failures.append(f"{len(errored)} parity result(s) errored")

            cross_backend = [
                result
                for result in results
                if result.get("status") == STATUS_AVAILABLE and result.get("device") != reference_device
            ]
            if not cross_backend:
                failures.append(
                    "no cross-backend parity comparison available (only the reference device was compared, "
                    "or no devices at all - this is not cross-backend evidence)"
                )

            bad_bitwise = [
                result
                for result in results
                if result.get("status") == STATUS_AVAILABLE
                and result.get("expectation_class") == EXPECTATION_BITWISE
                and result.get("verdict") != VERDICT_IDENTICAL
            ]
            if bad_bitwise:
                failures.append(f"{len(bad_bitwise)} bitwise-class parity result(s) are not identical")

            bad_tolerance = [
                result
                for result in results
                if result.get("status") == STATUS_AVAILABLE
                and result.get("expectation_class") != EXPECTATION_BITWISE
                and result.get("verdict") not in (VERDICT_IDENTICAL, VERDICT_WITHIN_TOLERANCE)
            ]
            if bad_tolerance:
                failures.append(f"{len(bad_tolerance)} tolerance-class parity result(s) exceeded tolerance")

    return failures


def main(argv: list | None = None) -> int:
    """CLI entry point.

    Exit code policy - see the module docstring for the full rationale:

    * 0 - report produced and delivered; no requested gate failed.
    * 1 - report produced but could not be written to --output.
    * 2 - argparse usage error (including `--elements` < 1), raised before
      any report is generated.
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
