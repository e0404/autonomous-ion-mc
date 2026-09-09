# 0001 — Cross-backend (Warp CPU vs CUDA) parity methodology

- Status: accepted
- Date: 2026-09-09
- Task: DRYRUN-002
- Affects: validation strategy, numerical accuracy, reproducibility

## Problem

`EXPERIMENT.md` requires that the Warp CPU and Warp CUDA execution paths
"preserve the same physical model rather than becoming independently
maintained implementations whose behavior can silently diverge", and lists
"comparison between Warp CPU and CUDA execution" as a validation method.
`REQUIREMENTS.md` additionally requires cross-backend consistency testing.

Before any transport physics exists, the project needs to decide *what
"the same result" means* across backends. Two backends executing the same
Warp kernel source will not generally produce bit-identical floating-point
output, so an equivalence criterion has to be chosen deliberately. Choosing
it later, under pressure from a failing physics comparison, invites picking
whatever tolerance makes the current test pass.

## Context

- Warp compiles the same kernel source to two different targets: a host
  C++/LLVM toolchain for the CPU device and NVRTC/nvcc for CUDA devices.
- The two targets legitimately differ in at least two ways that affect
  results without either being wrong:
  1. **Fused multiply-add contraction.** A `a * b + c` expression may be
     contracted into a single FMA with one rounding instead of two,
     independently on each target. On CUDA this is nvcc's default
     (`-fmad=true`); it is not controlled by Warp's `fast_math` module
     option, which is off by default.
  2. **Transcendental function implementations.** `sqrt`, `exp`, `sin`
     etc. are separate implementations on each target. IEEE 754 requires
     correct rounding only for `sqrt`; the others are permitted to differ
     by a small number of ULPs.
- Both effects are properties of the *toolchain*, not of the physics, and
  are not fixed by writing better transport code.

## Candidate approaches

1. **Bitwise equality only.** Simple, unambiguous, no tolerance to argue
   about. Rejected: it is not achievable for any kernel using
   transcendentals, so it would either forbid transcendentals in physics
   kernels or force the criterion to be abandoned exactly when it starts
   to matter.
2. **A single global relative tolerance.** One number for all comparisons.
   Rejected: a tolerance loose enough for `float32` transcendentals is far
   too loose to detect a genuine algorithmic divergence in `float64`
   arithmetic, where agreement should be near-exact. A single number
   silently converts the strongest available check into the weakest.
3. **Two expectation classes (arithmetic → bitwise; transcendental →
   tolerance), pure relative tolerance.** Implemented first and
   **empirically falsified** — see *Validation outcome*. Two classes were
   too few, and a pure relative tolerance proved to be the wrong measure.
4. **Three expectation classes, mixed absolute+relative criterion
   (selected).** Classify each workload by what the *toolchain* is
   permitted to do to it, and hold each class to the tightest criterion
   that class can actually meet.

## Selected approach

Approach 4, implemented as `infrastructure/diagnostics/warp_backend_report.py`
(documented in `docs/warp_backend_diagnostics.md`).

### Expectation classes

| workload | expression | class | criterion |
|---|---|---|---|
| `product` | `x * x` | `bitwise` | bitwise identical, no tolerance |
| `fma_exposed` | `(x * x - x) + 0.5` | `contraction_sensitive` | mixed abs+rel |
| `transcendental` | `sqrt(x) * exp(-x) + sin(x)` | `transcendental` | mixed abs+rel |

`product` is a single multiply: one rounding, no adjacent add or subtract
to be contracted with, and nothing to reassociate. It is the control that
proves the instrument and the two toolchains agree exactly when nothing
licenses them to differ. `fma_exposed` deliberately *retains* an FMA
opportunity so that contraction differences are measured and visible
rather than silently absorbed; it is evidence, not a defect.

### Comparison criterion

Tolerance-based classes use an elementwise **mixed** criterion

    |candidate - reference| <= atol + rtol * |reference|

rather than a pure relative tolerance, because kernel output that passes
through zero makes a pure relative measure diverge (see *Validation
outcome*). Per (workload, dtype):

| workload | dtype | rtol | atol |
|---|---|---|---|
| `product` | both | n/a | n/a |
| `fma_exposed` | `float32` | 1e-6 | 1e-6 |
| `fma_exposed` | `float64` | 1e-12 | 1e-14 |
| `transcendental` | `float32` | 1e-6 | 1e-6 |
| `transcendental` | `float64` | 1e-12 | 1e-14 |

The criterion is reported as a single verdict-relevant number,
`max_normalized_diff` = max over elements of
`|c - r| / (atol + rtol*|r|)`, which passes iff it is `<= 1`.
`max_abs_diff` and `max_rel_diff` are still reported, but `max_rel_diff`
is explicitly informational: it is inflated wherever the reference value
approaches zero and must never be cited as a precision figure.

### Other properties

- Warp **CPU is the reference device**. Every device, including the CPU
  itself, is compared against it. The CPU-vs-CPU comparison is retained
  as a self-consistency control: it must always be `identical`, and if it
  is not, the instrument — not the backend — is at fault. Because that
  control is trivially satisfiable, it does **not** count as
  cross-backend evidence: `--require-parity` additionally demands at
  least one successful comparison against a non-reference device.
- **NaN is never a pass.** Any NaN in either output forces
  `exceeds_tolerance` regardless of bit patterns, and a `nan_count` is
  reported. A NaN means the kernel produced no valid result; two invalid
  results agreeing is not equivalence.
- The `bitwise` class admits only `identical` or `exceeds_tolerance`.
  Values that are numerically equal but bit-distinct (`-0.0` vs `0.0`)
  are a failure for a class defined as bitwise.
- Inputs are generated deterministically from an explicit seed in pure
  Python (`random.Random`), reproducible without numpy and independent of
  the numpy version. An empty input vector is an error, never a vacuous
  pass.
- Being a diagnostic, producing a report is success. The criterion is
  *enforced* only when a caller opts in via `--require-cpu` /
  `--require-cuda` / `--require-parity`, which is how a future validation
  gate consumes it.
- The report embeds its own `configuration` section (elements, seed,
  workloads, classes, domain, all tolerances), so a recorded report is
  interpretable and reproducible without the code version that produced it.

## Expected tradeoffs

- Three classes must be maintained as physics kernels grow, and every new
  kernel must be classified. This is intended: classifying a kernel forces
  an explicit statement about what numerical behavior it is allowed to have.
- The tolerances are deliberately margined rather than tight (measured
  margins in *Validation outcome*). They will catch an algorithmic or
  dispatch divergence; they will not catch a change of a few ULPs in a
  library's `sin`. That is the right side to err on for a divergence
  detector, but it means these tolerances are not a precision benchmark
  and must not be cited as one.
- The criterion is elementwise and deterministic. It is **not**
  appropriate as written for stochastic per-history transport output,
  which will need a statistical criterion. This decision covers
  deterministic kernel-level parity only.
- The CPU is the reference, which means CPU results are treated as
  definitionally correct for parity purposes. Parity is an equivalence
  check, not an accuracy check: both backends agreeing says nothing about
  either being physically right. Accuracy is a separate validation axis
  (analytic and reference-data comparison).
- `fma_exposed` is expected to differ by ~1 ULP indefinitely. It is a
  monitored, documented difference, not a bug to be fixed; if it ever
  becomes bitwise identical, that is itself informative (a toolchain or
  flag change).

## Validation strategy

The instrument is validated on two levels:

1. **Logic**, by unit tests that need neither `warp` nor `numpy`
   (`tests/infrastructure/test_warp_backend_report.py`): bit-pattern
   comparison including `-0.0`/`0.0`, NaN and inf; difference
   computation; per-class verdict selection at and across the criterion
   boundary; deterministic input generation; guarded degradation;
   vacuous-pass rejection; gate exit codes.
2. **Execution**, by running the report on the experiment workstation's
   real GPU through the controlled host runner (`run_host_validation`),
   which is the only path to GPU execution from the agent sandbox.

## Validation outcome

The first implementation of this decision used approach 3 (two classes,
pure relative tolerance) and was executed on the experiment workstation
GPU (NVIDIA RTX A6000, sm_86, CUDA Toolkit 12.9, driver 13.1, Warp
1.17.0; 4096 elements, seed 12345, domain [0.05, 4.0]) at task-branch SHA
`a3b8c0aec7a1b6901c9a17f48233214aa2799c18`. It measured, against the
Warp CPU reference:

| workload | dtype | max_abs_diff | max_rel_diff | verdict |
|---|---|---|---|---|
| `affine` (`(x*x - x) + 0.5`) | `float32` | 9.5367431640625e-07 | 1.1920788269276779e-07 | exceeds_tolerance |
| `affine` | `float64` | 1.7763568394002505e-15 | 2.2198263041838353e-16 | exceeds_tolerance |
| `transcendental` | `float32` | 1.1920928955078125e-07 | 7.615368005429083e-06 | within_tolerance (tol 1e-5) |
| `transcendental` | `float64` | 2.220446049250313e-16 | 3.384892820249319e-14 | within_tolerance (tol 1e-9) |

All four CPU-vs-CPU self-comparisons were bitwise identical with diffs of
exactly 0.0, so the instrument itself is self-consistent.

Two findings, both of which changed this decision:

1. **The bitwise expectation was falsified, because the kernel was not
   what its own documentation claimed.** The arithmetic kernel was
   `(x * x - x) + 0.5`, documented as having "no single
   fused-multiply-add opportunity". That was simply wrong: `x * x - x`
   *is* `fma(x, x, -x)`. The two toolchains do not contract it
   identically, and the resulting ~1 ULP difference (1.19e-07 relative in
   `float32`, 2.22e-16 in `float64` — both almost exactly one ULP of
   their dtype) failed the 0.0 tolerance.

   The response was **not** to loosen the tolerance. The bitwise class was
   given a genuinely contraction-free control (`product` = `x * x`, a
   single rounded operation), and the FMA-exposed expression was kept as
   its own `contraction_sensitive` class with a mixed tolerance. The
   project now measures the contraction difference deliberately instead of
   asserting it away.

2. **A pure relative tolerance was the wrong measure.** The reported
   `float32` transcendental figure of 7.6e-06 is an artifact:
   `sqrt(x)*exp(-x) + sin(x)` crosses zero at x ≈ 3.2, inside the sampled
   domain, so a ~1 ULP *absolute* difference (1.19e-07) divided by a
   near-zero reference inflates into a large *relative* one. That number
   describes the location of a zero crossing, not the agreement of the
   transcendental implementations — and it sat at only ~25% margin below
   its tolerance, so the gate was fragile for a reason unrelated to what
   it was measuring.

   The criterion is therefore now mixed (`atol + rtol*|ref|`), the domain
   deliberately keeps the zero crossing (near-zero output is a normal
   numerical case, not one to dodge), and `max_normalized_diff` was added
   so the headline number is the one the verdict actually uses.

Measured margins of the selected tolerances against the data above:
`fma_exposed`/`float32` has a budget of ~1.35e-05 against a measured
9.54e-07 (≈14x); `fma_exposed`/`float64` ~1.26e-11 against 1.78e-15
(≈7000x); `transcendental`/`float32` ~2.2e-06 against 1.19e-07 (≈18x);
`transcendental`/`float64` ~1.2e-12 against 2.22e-16 (≈5000x). The
`float64` margins are deliberately generous: the measured differences are
at the 1 ULP floor, where a tight tolerance would track toolchain
versions rather than divergence.

Both findings were independently confirmed as blocking by the external
Codex reviewer, which additionally identified that NaN-vs-NaN could
report `identical` and that `--require-parity` could pass vacuously on an
empty input vector or on a CPU-only self-comparison. Those are addressed
in *Selected approach* above.

The methodological point for the experiment record: the two-class
criterion worked exactly as intended. A bitwise class with a zero
tolerance converted a silent 1-ULP toolchain difference into a visible,
diagnosable event, and it did so on the *first* real GPU execution,
before any physics depended on it. A single global tolerance (candidate
2) would have absorbed both findings without anyone noticing that a
kernel was not doing what its documentation claimed.

The revised criterion is re-validated on real hardware at the final
task-branch SHA; see the exact-SHA local validation record for this task.
