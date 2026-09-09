# Warp Backend Diagnostics Utility

`infrastructure/diagnostics/warp_backend_report.py` reports the NVIDIA Warp
backend runtime and its cross-backend numerical consistency (Warp CPU vs
Warp CUDA) as a single machine-readable JSON document. It is the project's
first cross-backend (Warp CPU vs Warp CUDA) consistency instrument (see
`REQUIREMENTS.md`, "Cross-backend statistical consistency", and
`EXPERIMENT.md`'s validation list: "comparison between Warp CPU and CUDA
execution"), so its numbers are treated as validation evidence: the report
documents everything needed to reproduce them.

`warp` and `numpy` are both optional imports, used lazily only inside the
execution path. The module imports cleanly, and all of its
comparison/verdict logic works correctly, with neither installed - GitHub
CI installs only pytest.

## What it reports

| Key | Contents |
|---|---|
| `schema_version` | Integer report schema version (currently `1`). Not a section - always present, no `status`. |
| `generated_at_utc` | ISO-8601, timezone-aware UTC timestamp of report generation. Not a section. |
| `warp` | `version` from an optional `import warp`. `unavailable` (not `error`) if `warp` is simply not importable. |
| `devices` | `devices`: list of `{alias, is_cpu, is_cuda, name, arch, total_memory_bytes, free_memory_bytes}` for every device Warp itself enumerates. `unavailable` if Warp enumerates zero devices. See [Device metadata normalization](#device-metadata-normalization). |
| `parity` | The cross-backend comparison described below. `unavailable` if there is no Warp CPU device to use as the reference. |
| `configuration` | The effective run parameters actually used (element count, seed, workloads, dtypes, domain, tolerances), so a reported result is reproducible from the report alone. Always `available` - it documents this run's own parameters rather than probing anything that can fail. |

## Status vocabulary

Identical to `infrastructure/diagnostics/environment_report.py` (see
`docs/environment_diagnostics.md`):

- **`available`** - the probe succeeded; the section's data fields are populated.
- **`unavailable`** - the underlying capability is simply absent (`warp` not
  importable, zero Warp devices enumerated, no Warp CPU device available as
  the parity reference). Expected, not a failure.
- **`error`** - the probe was attempted but failed unexpectedly. A `detail`
  field explains why.

Every probe and every individual (device, workload, dtype) parity
comparison is independently guarded: one failing device, workload, dtype,
or `warp` being entirely absent never prevents the rest of the report from
being produced. An unexpected Warp device-API shape (e.g. a renamed or
missing attribute) degrades the affected field to `null` via `getattr`/try
guards rather than raising.

## The parity workloads

Three workloads, each in `float32` and `float64`, are run on every device
Warp enumerates, and each device's result is compared elementwise against
the **Warp CPU** device's result for the same workload/dtype (Warp CPU is
always the reference). If no Warp CPU device is available, `parity` is
`unavailable` with a clear `detail` - there is nothing to compare against.

Each workload has an explicit **expectation class**, reported per result as
`expectation_class` and in `configuration.expectation_classes`:

| workload | expression | expectation class | why |
|---|---|---|---|
| `product` | `y[i] = x[i] * x[i]` | `bitwise` | A single multiply, one rounding, nothing adjacent to contract with and nothing to reassociate. This is the actual bitwise control - the thing that can honestly be expected to reproduce bit-for-bit across Warp CPU and Warp CUDA. |
| `fma_exposed` | `y[i] = (x[i] * x[i] - x[i]) + 0.5` | `contraction_sensitive` | Deliberately exposes a fused-multiply-add (FMA) contraction opportunity: `x*x - x` is a multiply immediately followed by a subtraction of the same operand, i.e. exactly `fma(x, x, -x)`. Different toolchains (CPU C++ vs CUDA `nvcc`) are free to contract this differently, producing a small, real, non-bug difference. This workload exists specifically to measure and report that difference rather than hide it. |
| `transcendental` | `y[i] = sqrt(x[i]) * exp(-x[i]) + sin(x[i])` | `transcendental` | CPU and CUDA use different transcendental-function implementations, which are IEEE-754-compliant but not required to be identically rounded. |

### History: why there are three workloads, not two

Round 1 of this diagnostic had only `affine` (== today's `fma_exposed`) and
`transcendental`, and claimed `affine` had "no single fused-multiply-add
opportunity" and should therefore be bitwise identical. **That claim was
false.** A real GPU run falsified it directly - see
[Measured evidence](#measured-evidence) below: `x*x - x` **is** an FMA
opportunity, and it measured a ~1 ULP difference between Warp CPU and Warp
CUDA, failing the round-1 `0.0` tolerance it was never actually able to
meet. Rather than deleting that expression, it was kept as its own
workload (`fma_exposed`) with an honest `contraction_sensitive` class - its
contraction sensitivity is real, useful validation evidence - and a new,
genuinely single-operation `product` workload was added as the actual
bitwise control.

### Inputs

The input vector is generated by `generate_inputs(n, seed)` using
`random.Random(seed).uniform(DOMAIN_MIN, DOMAIN_MAX)` - pure Python, not
numpy - so it is reproducible from `configuration.seed` and
`configuration.elements` alone, without numpy. The domain is `[0.05, 4.0]`
(`DOMAIN_MIN`/`DOMAIN_MAX`), kept strictly positive so `sqrt` in
`transcendental` is always defined. `--elements` sets `n` (default `4096`,
chosen to run in well under a second on either backend) and **must be
`>= 1`** - `0` or a negative value is rejected as an argparse usage error
(exit code `2`); an empty comparison would otherwise make the `bitwise`
class's `identical` verdict vacuously true (nothing to disagree on).
`build_parity_section()` enforces the same `>= 1` constraint at the library
level (returning `status: "error"`), for any caller that does not go
through the CLI. `--seed` sets the seed (default `12345`).

### Per-entry fields

Each entry in `parity.results` is one `(device, workload, dtype)` triple:

```jsonc
{
  "device": "cuda:0",
  "workload": "transcendental",
  "expectation_class": "transcendental",
  "dtype": "float32",
  "status": "available",
  "n": 4096,
  "nan_count": 0,
  "max_abs_diff": 1.1920928955078125e-07,
  "max_rel_diff": 7.615368005429083e-06,
  "max_normalized_diff": 0.24,
  "bitwise_identical": false,
  "verdict": "within_tolerance",
  "rtol": 1e-06,
  "atol": 1e-06
}
```

- **`max_abs_diff`** - the largest elementwise `|candidate - reference|`.
- **`max_rel_diff`** - **informational only, not verdict-relevant** (see
  [The `max_rel_diff` caveat](#the-max_rel_diff-caveat) below).
- **`max_normalized_diff`** - the verdict-relevant metric for the
  tolerance-based classes: the largest elementwise value of the mixed
  absolute+relative criterion (see [The mixed
  criterion](#the-mixed-criterion) below), which passes iff `<= 1.0`.
  `null` for the `bitwise` class (`product`), where no numeric budget
  applies.
- **`nan_count`** - number of elements where the reference or the
  candidate (or both) was `NaN`. Any nonzero `nan_count` forces
  `verdict: "exceeds_tolerance"` (see [Verdict
  semantics](#verdict-semantics)).
- **`rtol`/`atol`** - the tolerance budget actually applied (see
  [The mixed criterion](#the-mixed-criterion)); `null`/`null` for the
  `bitwise` class.

A `(device, workload, dtype)` entry that raised (e.g. a kernel launch
failure on one specific device) gets `status: "error"` and a `detail`
string instead of the numeric fields, and does not affect any other entry.

### The mixed criterion

Tolerance-class (`contraction_sensitive`, `transcendental`) verdicts use the
elementwise mixed absolute+relative criterion, applied per element as:

```
|candidate - reference| <= atol + rtol * |reference|
```

expressed as a single number via `normalized_diff()`:

```
max_normalized_diff = max_i( |c_i - r_i| / (atol + rtol * |r_i|) )
```

which passes (is verdict-relevant "within tolerance") iff
`max_normalized_diff <= 1.0`. A pure relative tolerance was deliberately
NOT used, because a reference value near zero can turn an ordinary, tiny
absolute difference into a huge relative one - see the caveat below for the
concrete measured case that motivated this.

#### Chosen `(rtol, atol)` constants

Defined as named module constants (`TOLERANCE` in `warp_backend_report.py`):

| workload | dtype | `rtol` | `atol` | rationale |
|---|---|---|---|---|
| `product` | float32 / float64 | `null` | `null` | `bitwise` class - no numeric budget applies; verdict is decided purely by bit-pattern equality. |
| `fma_exposed` | float32 | `1e-6` | `1e-6` | One FMA-contraction rounding step (~1 ULP). Measured (see below): `max_abs_diff` 9.5e-07, `max_rel_diff` 1.19e-07. The chosen budget sits comfortably above that. |
| `fma_exposed` | float64 | `1e-12` | `1e-14` | Same reasoning, double precision. Measured: `max_abs_diff` 1.78e-15, `max_rel_diff` 2.22e-16. |
| `transcendental` | float32 | `1e-6` | `1e-6` | Measured: `max_abs_diff` 1.19e-07 (same ~1 ULP order as `fma_exposed`); the `max_rel_diff` figure (7.62e-6) is inflated by a near-zero reference, see the caveat below - the absolute figure is what the budget is actually sized against. |
| `transcendental` | float64 | `1e-12` | `1e-14` | Measured: `max_abs_diff` 2.22e-16, `max_rel_diff` 3.38e-14. Same budget as `fma_exposed`/float64. |

#### Measured evidence

These constants were validated against a real GPU run (task DRYRUN-002,
NVIDIA RTX A6000, sm_86, CUDA 12.9, driver 13.1, Warp 1.17.0, 4096
elements, seed 12345), comparing Warp CUDA against Warp CPU:

| workload/dtype | `max_abs_diff` | `max_rel_diff` | bitwise? |
|---|---|---|---|
| `fma_exposed`/float32 | `9.5367431640625e-07` | `1.1920788269276779e-07` | no |
| `fma_exposed`/float64 | `1.7763568394002505e-15` | `2.2198263041838353e-16` | no |
| `transcendental`/float32 | `1.1920928955078125e-07` | `7.615368005429083e-06` | no |
| `transcendental`/float64 | `2.220446049250313e-16` | `3.384892820249319e-14` | no |
| `product` (both dtypes) | `0.0` | `0.0` | **yes** |

All four workloads' CPU-vs-CPU self-comparisons were bitwise identical with
diffs exactly `0.0`, as expected.

This is also the direct evidence behind the workload split: `fma_exposed`
was NOT bitwise identical on real CUDA hardware (confirming the FMA
contraction), while `product` was.

#### The `max_rel_diff` caveat

`max_rel_diff` is reported for context but is **never used to decide a
verdict**. The `transcendental`/float32 row above is the concrete example
of why: `sqrt(x)*exp(-x) + sin(x)` crosses zero near `x ~= 3.2` within the
`[0.05, 4.0]` input domain, and an ordinary ~1 ULP absolute difference
(`1.19e-07`, the same order as `fma_exposed`'s) landing near that crossing
inflates into a `max_rel_diff` of `7.62e-6` - about 76% of what would have
been a `1e-5` pure-relative-tolerance budget, purely from proximity to
zero, not from any real loss of precision. Do not cite `max_rel_diff` as a
precision figure; cite `max_abs_diff` and/or `max_normalized_diff`.

### Verdict semantics

- **`bitwise` class** (`product`): verdict is `identical` **iff** bit-for-bit
  identical AND no `NaN` is present anywhere in either sequence; otherwise
  `exceeds_tolerance`. This class **never** returns `within_tolerance` -
  even a numerically-equal-but-bit-distinct pair (e.g. `-0.0` vs `0.0`) is a
  genuine failure for a class whose entire point is bitwise
  reproducibility.
- **Tolerance classes** (`contraction_sensitive`, `transcendental`):
  `identical` when bitwise identical and NaN-free; `within_tolerance` when
  NaN-free and `max_normalized_diff <= 1.0`; otherwise `exceeds_tolerance`.
- **NaN is never a pass**, in any class. Any `NaN` in the reference or the
  candidate output forces `exceeds_tolerance`, regardless of bit patterns.
  This is checked before bitwise identity is even consulted: `bits_equal()`
  itself still treats any NaN as bit-equal to any other NaN (it remains a
  strict, self-contained bit-pattern predicate on its own), but
  `select_verdict()` checks `nan_count` first, specifically so a `NaN` vs
  `NaN` comparison can never slip through as `identical`. `nan_count` is
  reported on every result so this is auditable.

### Bitwise comparison

`bits_equal(a, b)` performs a true bit-pattern comparison via
`struct.pack("<d", ...)`, not `==`: it distinguishes `-0.0` from `0.0`
(`-0.0 == 0.0` is `True` in IEEE 754/Python, but their bit patterns differ),
and treats any `NaN` as bit-equal to any other `NaN` (a `NaN` is never
bit-equal to a non-`NaN`) - see [Verdict semantics](#verdict-semantics) for
why that does not let NaN pass a verdict. `abs_diff`/`rel_diff`/
`normalized_diff` are similarly guarded: numerically-equal values
(including two equal infinities) have a difference of `0.0`; a `NaN`
operand always produces a `NaN` difference, never raises.

All of this comparison/verdict logic (`bits_equal`, `abs_diff`, `rel_diff`,
`normalized_diff`, `compare_sequences`, `select_verdict`) is pure Python
over plain `list[float]` sequences - `numpy`/`warp` are only ever used to
move data into and out of Warp arrays (in `WarpKernelRunner`), never inside
the comparison logic itself.

## Device metadata normalization

Two fields are normalized rather than reported raw, to avoid a `null`-like
value being mistaken for real hardware metadata:

- **`arch`** is reported as `null` for any device where `is_cuda` is not
  `True`. A real GPU-host run of this diagnostic reported the Warp CPU
  device's `arch` as the string `"0"`, which is not an architecture
  identifier and would mislead a reader into thinking it was one; `null`
  is reported instead for any non-CUDA device, and the real value (e.g.
  `"sm_86"`) is preserved for CUDA devices.
- **`total_memory_bytes`/`free_memory_bytes`**: a value of exactly `0` is
  normalized to `null`. The same run reported the CPU device's memory as
  `0`/`0` bytes; Warp reports `0` for these fields when the underlying
  quantity is unavailable (on that host, the CPU memory query requires
  `psutil`, which was not installed in the host-runner environment and
  also emitted a warning about it on stderr) - `null` honestly states
  "unknown" rather than asserting a real zero-byte device.

## Injection seam for testing

`build_devices_section(warp_module=None)` and `build_parity_section(devices,
run_kernel, ...)` both take their Warp module / kernel-execution callable as
an explicit parameter rather than reaching for a global. `build_report(...,
warp_module=None, run_kernel=None)` exposes the same seam at the top level.
This lets tests substitute a fake object implementing just enough of the
Warp device API (or a plain-Python `run_kernel` callable) to exercise full
report assembly - including the parity comparison and verdict logic -
without `warp` (or `numpy`) installed, and without monkeypatching any
private module internals.

## Exit code policy

- **`0`** - the report was produced and delivered (to stdout or
  `--output`), and no requested `--require-*` gate failed.
- **`1`** - the report was produced but could not be written to the
  `--output` path. A structured JSON object of the form
  `{"status": "error", "detail": "..."}` is printed to stderr; no
  traceback is ever printed.
- **`2`** - an argparse usage error (including `--elements` `< 1`), raised
  by argparse itself before any report is generated.
- **`3`** - the report was produced and delivered successfully, but a
  requested gate flag's expectation was not satisfied (see below).
  Producing a diagnostic report is always success by itself; gates are an
  opt-in, separate expectation layered on top of that, consistent with
  `environment_report.py`'s exit-code policy of "producing a report is
  success."

### Gate flags

By default, no gates are applied - this is a pure diagnostic tool. Three
optional flags turn specific expectations into a non-zero exit (`3`) when
unmet, useful for a validation/CI step that wants to assert a GPU-backed
environment is actually present and consistent:

- **`--require-cpu`** - fail if no Warp CPU device is available.
- **`--require-cuda`** - fail if no Warp CUDA device is available.
- **`--require-parity`** - fail unless **all** of the following hold:
  1. the `parity` section's `status` is `available`;
  2. no result has `status: "error"`;
  3. at least one `available` result compares a device **other than** the
     reference device - a CPU-only run (the CPU reference compared only to
     itself) is not cross-backend evidence and does not satisfy this gate
     on its own, even if every individual result "passes";
  4. every `bitwise`-class result has `verdict: "identical"`;
  5. every non-`bitwise`-class result has `verdict` `"identical"` or
     `"within_tolerance"`.

## No secrets, no host-identifying personal data

The report never includes environment variables, tokens, remote URLs, or
user names, and does not read or dump `os.environ`.

## Running it

```bash
# Pretty JSON to stdout (default indent: 2)
python3 infrastructure/diagnostics/warp_backend_report.py

# Custom element count / seed (elements must be >= 1)
python3 infrastructure/diagnostics/warp_backend_report.py --elements 8192 --seed 7

# Write to a file instead of stdout
python3 infrastructure/diagnostics/warp_backend_report.py --output /path/to/report.json

# Gate on GPU + cross-backend parity for a CI/validation step
python3 infrastructure/diagnostics/warp_backend_report.py --require-cuda --require-parity
```

## Running its tests

```bash
uv run --with pytest --no-project python -m pytest tests -q
```

`tests/conftest.py` adds `infrastructure/diagnostics/` to `sys.path`, the
same convention `environment_report.py`'s tests use. Tests are deterministic
and never assume `warp`/`numpy` are present OR absent on the host: graceful
degradation is exercised via an injected fake Warp module/kernel-runner
(`warp` absent) and via `builtins.__import__` monkeypatching (simulating an
import failure), matching `environment_report.py`'s test style. A handful of
tests that need a real `warp` install call `pytest.importorskip("warp")`
first and make no assumption about which devices that installation exposes.
