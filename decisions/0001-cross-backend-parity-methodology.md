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

- Warp compiles the same kernel source to two different targets: a C++
  host compiler for the CPU device and NVRTC/nvcc for CUDA devices.
- The two targets legitimately differ in at least two ways that affect
  results without either being wrong:
  1. **Fused multiply-add contraction.** A `a * b + c` expression may be
     contracted into a single FMA with one rounding instead of two,
     independently on each target, depending on compiler flags.
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
   Rejected: a tolerance loose enough for `float32` transcendentals
   (~1e-5) is far too loose to detect a genuine algorithmic divergence in
   `float64` arithmetic, where agreement should be near-exact. A single
   number silently converts the strongest available check into the weakest.
3. **Two expectation classes, per dtype (selected).** Classify each
   comparison by what the *toolchain* can legitimately do to it, and hold
   each class to the tightest criterion that class can actually meet:
   - arithmetic-only kernels → expected **bitwise identical** (tolerance 0);
   - transcendental-invoking kernels → expected **within a documented
     per-dtype relative tolerance**.

## Selected approach

Approach 3, implemented as `infrastructure/diagnostics/warp_backend_report.py`
(documented in `docs/warp_backend_diagnostics.md`):

- Warp **CPU is the reference device**. Every device, including the CPU
  itself, is compared elementwise against it. Comparing the CPU to itself
  is retained deliberately as a self-consistency control: it must always
  report `identical`, and if it ever does not, the instrument itself — not
  the backend — is at fault.
- Every comparison reports `max_abs_diff`, `max_rel_diff`,
  `bitwise_identical`, the tolerance applied, and a `verdict` of
  `identical` / `within_tolerance` / `exceeds_tolerance`. The numbers are
  always reported, not just the verdict, so a tightening or loosening of
  a tolerance can be argued from recorded evidence rather than re-run.
- Inputs are generated deterministically from an explicit seed in pure
  Python (`random.Random`), so an input vector is reproducible without
  numpy and independent of the numpy version.
- The report carries its own `configuration` section (elements, seed,
  kernels, dtypes, domain, all tolerances), so a recorded report is
  interpretable and reproducible without reference to the code version
  that produced it.
- Being a diagnostic, producing a report is success. The equivalence
  criterion is only *enforced* when a caller opts in via `--require-cpu`
  / `--require-cuda` / `--require-parity`, which is how a future
  validation or CI-adjacent gate consumes it.

Tolerances are named module constants, not inline literals:

| kernel | dtype | relative tolerance | basis |
|---|---|---|---|
| `affine` | `float32` | 0.0 | arithmetic only; no legitimate source of divergence |
| `affine` | `float64` | 0.0 | as above |
| `transcendental` | `float32` | 1e-5 | a few ULPs of single precision (~1e-7) with margin |
| `transcendental` | `float64` | 1e-9 | good double-precision implementations typically agree to ~1e-12 |

## Expected tradeoffs

- Two classes must be maintained as physics kernels grow, and every new
  kernel must be classified. This is intended: classifying a kernel forces
  an explicit statement about what numerical behavior it is allowed to have.
- The transcendental tolerances are deliberately margined rather than
  tight. They will catch an algorithmic or dispatch divergence; they will
  not catch a change of a few ULPs in a library's `sin`. That is the right
  side to err on for a divergence detector, but it means these tolerances
  are not a precision benchmark and must not be cited as one.
- Tolerances are relative and elementwise-max. This is appropriate for
  smooth kernel output on a well-scaled domain. It is *not* appropriate as
  written for stochastic per-history transport output or for quantities
  that legitimately pass through zero; those will need their own criterion
  (statistical, or absolute-plus-relative). This decision covers
  deterministic kernel-level parity only.
- The CPU is the reference, which means CPU results are treated as
  definitionally correct for parity purposes. Parity is an equivalence
  check, not an accuracy check: both backends agreeing says nothing about
  either being physically right. Accuracy is a separate validation axis
  (analytic and reference-data comparison).

## Validation strategy

The instrument is validated on two levels:

1. **Logic**, by unit tests that need neither `warp` nor `numpy`
   (`tests/infrastructure/test_warp_backend_report.py`): bit-pattern
   comparison including `-0.0`/`0.0`, NaN and inf; difference
   computation; verdict selection at and across the tolerance boundary;
   deterministic input generation; guarded degradation; gate exit codes.
2. **Execution**, by running the report on the experiment workstation's
   real GPU through the controlled host runner (`run_host_validation`),
   which is the only path to GPU execution from the agent sandbox.

## Validation outcome

Pending. Logic-level unit tests pass. The execution-level GPU run through
`run_host_validation` and its empirical parity numbers are recorded in the
exact-SHA local validation record for this task; this section is completed
from that run's actual output, not in advance of it.
