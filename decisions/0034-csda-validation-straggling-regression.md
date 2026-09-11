# 0034 — Fix the CSDA depth-dose validation's straggling-default regression

- Status: accepted
- Date: 2026-09-11
- Task: DEV-029
- Affects: validation V1 (`validation/v1_depth_dose_csda.py`)

## Problem

The release validation suite (decision `0033`) surfaced a failing milestone gate:
`validation/v1_depth_dose_csda.py` failed three of its gates on the host runner —
`r80_vs_csda_range` (R80 +1.5 % vs the tabulated CSDA range), `step_size_convergence`
(step drift 0.47 % ≫ 0.05 %), and `range_cross_check` (+1.4 %). All cross-backend
gates passed and reference energy conservation was ~1e-15, so the transport physics
is correct; the failure is in the *validation configuration*.

Root cause: this validation is a **CSDA** (continuous-slowing-down, deterministic)
check — its name, docstring, and gates all compare the transport range to the
analytic CSDA range. It has been unchanged since DEV-004, where it constructed
`TransportEngine(table, slab, grid)` (pure CSDA, because straggling did not yet
exist). DEV-005 then added energy-loss straggling and made **`straggling=True` the
engine default**. From DEV-005 on, this validation silently began transporting *with*
straggling, which broadens the Bragg peak and pushes R80 ~1.5 % beyond the CSDA range
and makes the step-convergence check stochastic. The unit tests in
`tests/ionmc/test_transport.py` were written to pass `straggling=False` explicitly and
so remained green; only the host validation script carried the implicit dependency on
the old default. Confirmed against the DEV-004 validation record: `straggling=False`
reproduces its R80 exactly (+0.15/+0.09/−0.05 %).

A second, smaller point surfaced once straggling was turned off: the deterministic
Bragg peak is razor-sharp, so the **Warp float32 energy_balance** at 150 MeV is
1.335e-5 — just over the validation's `1e-5` warp-balance sub-tolerance. This is the
intrinsic float32 accumulation residual over a full proton range (the same magnitude
the V4/benchmark cross-backend budgets accept), not a defect.

## Decision

1. **The CSDA validation constructs its engines with `straggling=False`.** All three
   `TransportEngine` constructions in `v1_depth_dose_csda.py` (the main depth dose and
   the two step-convergence runs) are made explicit, matching the validation's stated
   CSDA intent and the already-correct `test_transport.py` fixtures. Straggling itself
   remains validated separately by `v1_bragg_straggling.py`.

2. **The Warp energy-conservation sub-tolerance uses the project's standard float32
   budget.** `WARP_ENERGY_BALANCE_TOL = 5e-5` replaces the inline `1e-5`, aligning this
   float32 energy-conservation check with the 5e-5 float32 budget used across V4 and
   the benchmarks (the float64 reference path stays at `1e-9`). This tolerance follows
   from float32 accumulation over a full proton range, not from the observed value.

No physics or transport-kernel change: the engine is correct and its unit tests were
already green. This is a validation-configuration and tolerance-calibration fix.

## Consequences

- `v1_depth_dose_csda` passes all seven gates again on the host (CPU + CUDA), and the
  release validation suite reaches `release_ready = true` — the last gate before the
  first tagged release.
- **Lesson recorded:** a validation must be explicit about engine flags whose defaults
  can change; relying on a default (`straggling`) let a later default flip silently
  invalidate a milestone gate without any test catching it. The gap that let this
  persist — the release validation suite not being re-run after every change — is
  itself now closed by decision `0033`.
