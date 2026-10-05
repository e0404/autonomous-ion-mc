# 0039 — Proton condensed-history transport engine: models, step algorithm and backends

- Status: accepted
- Date: 2026-10-03
- Task: V3-003A (core, reference backend), V3-003B (Warp backends)
- Affects: physical accuracy, numerical accuracy, reproducibility, architecture, validation strategy, performance

## Problem

Decision 0037 fixes the architecture (shared `@wp.func` physics, three
drivers, Philox streams, per-batch accumulation). The transport engine still
needs concrete physical models for the electromagnetic part of proton
transport, a step algorithm whose numerical behaviour can be falsified, a
fail-closed configuration contract, and a test matrix frozen before results
exist. Nuclear interactions (V3-005), LET (V3-004) and persistence (V3-006)
are deliberately out of scope here.

## Context

- Research reports `docs/research/em-physics.md` and
  `docs/research/warp-architecture.md`; the implementation plan was prepared
  by a scientific-planning agent on 2026-10-03 and ratified by the lead with
  four corrections (below). Gottschalk's scattering-power paper
  (Med. Phys. 37 (2010) 352, arXiv:0908.1413) was read in its source form.
- Verified Warp 1.17 facts (decision 0037, archived measurements):
  Python-scope execution of `@wp.func` with `wp.float64` typing matches
  numpy; concurrent CPU launches from Python threads stalled once;
  single-float32 accumulators are biased at 1e6 deposits.
- V3-002 provides `StoppingTable`s on a log grid from 1 MeV/u with CSDA
  ranges integrated from the same stopping powers.

## Candidate approaches and decisions

| Topic | Candidates | Selected | Rationale |
|---|---|---|---|
| Mean energy loss | S(E)·s; range-table inversion | inverse-range telescoping `Rinv(R0) − Rinv(R0 − ρs)`; for `ρs < 10⁻² R0` a midpoint-rule branch `S(E₀ − S(E₀)·s/2)·s` (fourth amendment 2026-10-04: the former linear `S(E₀)·s` branch below `10⁻³ R0` carried a first-order step-size bias; the midpoint branch is second order and, for such steps, more consistent than the telescoping form, whose table round-trip error accumulates per step) | telescoping is exact in CSDA up to the inverse-range table round-trip error (≤ 1e-5 relative per evaluation, which accumulates with the number of steps); the midpoint branch is second order in s/R and, below 10⁻² R₀, more step-consistent than telescoping (fourth amendment); monotone; the former float32-cancellation motive for a short branch disappeared with float64 energy bookkeeping (archived precision measurement remains as history) |
| Straggling | Gaussian only; Vavilov; Gaussian/Gamma/uniform switch (`G4IonFluctuations`); Gaussian-clamped/Gamma (`bohr_gauss_clamped_gamma_v1`); Gamma for every ratio (`bohr_gamma_v1`) | **Default (third amendment 2026-10-04) `bohr_gamma_v1`:** Bohr variance sampled as Gamma with shape k = (mean/σ)² and scale σ²/mean for every ratio (Marsaglia–Tsang, `u^(1/k)` boost for k < 1; rejection bounded by 64 attempts, exhaustion increments the fail-closed counter `straggling_rejection`); every sampled loss is capped at the remaining kinetic energy. The earlier switch `bohr_gauss_clamped_gamma_v1` (Gamma for mean/σ < 3, Gaussian clamped to [0, 2·mean] above) remains selectable | the Gamma branch reproduces the mean and the Bohr variance exactly; the Gaussian branch preserves the mean and reduces the variance by at most about 0.5 % (0.13 % of samples clamped at 0 at mean/σ = 3; the moment test allows 1 %); the `G4IonFluctuations` uniform branch has variance mean²/3 whatever the Bohr value and the truncated Gaussian near mean/σ = 2 loses variance, so they are not used; sum over steps is Gaussian by the central limit theorem; delta electrons are deposited locally so single-step tails do not reach dose |
| Multiple scattering | per-step Highland; Molière; differential Highland (Kanematsu); differential Molière (Gottschalk) | **differential Molière** scattering power `T_dM = f_dM(pv, p₁v₁)(E_s/pv)²/X_S`, `E_s = 15.0 MeV`, `f_dM = 0.5244 + 0.1975 lg(1−(pv/p₁v₁)²) + 0.2320 lg(pv) − 0.0098 lg(pv) lg(1−(pv/p₁v₁)²)` (clamped ≥ 0), scattering length `1/(ρX_S) = α N_A r_e² (Z²/A){2 ln(33219 (AZ)^{-1/3}) − 1}` Bragg-additive, applied as a Gaussian polar angle with a random hinge | step-size independent by construction (per-step Highland is not; a negative-control test proves the instrument has power); ranked best against Hanson theory and measurement in the source paper; the research report's label "differential Highland" for these coefficients was wrong and is corrected |
| Lateral displacement | explicit correlated sampling; random hinge | random hinge (move `a·s`, deflect, move `(1−a)·s`, `a` uniform) | reproduces the Fermi–Eyges second moments exactly without extra draws |
| Geometry traversal | `floor(p)` with nudge; incremental DDA | incremental DDA with the voxel index as state and plane snapping on crossing | the nudge variant stalled in the archived toy kernel |
| Scoring | split steps at scoring planes; midpoint deposit; path-length-proportional deposit | **deposit apportioned along both hinge legs with a linear stopping-power ramp S(E₀) → S(E₁)** (amendments 2026-10-04, see below; initially path-length-proportional); the former `max_step ≤ min scoring spacing` rule is replaced by a computed, fail-closed bound on the number of scoring pieces per leg (second amendment 2026-10-04, see below) | simplest unbiased choice at the enforced step size; grid refinement/misalignment is a falsification probe in V3-011 |
| End of range | transport to zero; local deposition below `E_cut` | local deposition below `E_cut = 2 MeV` (protons) | residual range ≈ 0.07 mm in water; the V3-002 table floor (1 MeV/u) must stay ≤ `E_cut/2` |
| Uniform mapping | `((w>>8)+0.5)·2⁻²⁴` in both precisions | float64: 24-bit form; float32: `((w>>9)+0.5)·2⁻²³` | the 24-bit form rounds to exactly 1.0 in float32 for the top word (amends decision 0037) |
| CPU parallelism | Python threads; processes | spawned processes with private accumulators, float64 summation by the parent | archived stall of thread-concurrent CPU launches; Warp CPU atomics are plain read-modify-write |
| Nuclear flag | default off | no default; `nuclear=True` raises until V3-005 | fail-closed capability contract (V2-CAP) |

Engineering defaults (recorded in every result, probed by the falsification
tests): `E_cut` 2 MeV, maximum energy-loss fraction per step 0.02, maximum
step 1 mm, Geant4 range step function α = 0.2, ρ_f = 0.1 mm, 20 batches.

## Step algorithm (one history)

1. If `E ≤ E_cut`: deposit `E` locally (cutoff tally) and stop.
2. Fetch `S(E₀)`, `R₀`; compute the DDA distance to the next voxel plane,
   the energy-loss, range and maximum step limits; choose the step `s`.
3. Draw block A (hinge fraction, polar radius, azimuth, reserved); block A
   is always drawn so paired probes share streams.
4. Scattering (if enabled): `E_mid = Rinv(R₀ − ρs/2)`, variance
   `T_dM(E_mid)·s`, polar angle from the 2-D Gaussian, hinge at `a·s`; the
   second leg is cut at a voxel plane if needed (the particle then snaps to
   the plane and the voxel index is incremented on that axis).
5. Energy loss for the travelled length: mean by telescoping, fluctuation
   from block B (further blocks on rejection, at most 64; exceeding the
   limit is a fail-closed counter); `E₁ = E₀ − loss ≥ 0`.
6. Deposit `E₀ − E₁` into every scoring grid along both hinge legs, per
   scoring voxel by the path length inside the voxel WEIGHTED with the
   normalized linear stopping-power ramp S(E₀) → S(E₁) over the step (per-grid
   incremental DDA; outside-grid deposits are tallied separately). Until
   2026-10-04 this was a point deposit at the hinge-path midpoint, then
   proportional to path length — see the amendments "T9 at full scale" below.
7. Leaving the world tallies `E₁` as escaped. Step-count truncation or a
   stall tallies `E₁` as truncated, increments a counter and invalidates the
   result (decision 0037).

Random streams: key = 64-bit seed; counter = (history index, genealogy id,
block index, purpose) with purpose 1 for source sampling and 0 for transport.

## Expected tradeoffs

- `T_dM` depends on the initial `pv` of each particle, which is therefore
  part of the particle state; secondaries (V3-005) take their own value at
  birth. `f_dM` is clamped at the singular point `pv → p₁v₁`.
- The per-step variance uses the midpoint rule `s·T_dM(E_mid)`. On the first
  step of a particle's life `1 − (pv/p₁v₁)²` grows linearly from zero, so
  `f_dM` has an integrable logarithmic singularity that no fixed-order
  quadrature over the step removes (measured during V3-003A: midpoint and
  Simpson rules both missed the frozen U5 bound of 2e-3 by up to 1.2e-2 at
  x/R₁ = 0.05 with 1 mm steps). The birth step therefore uses the analytic (linearized) step
  average of the logarithmic term, `lg(1 − (pv(E₁)/p₁v₁)²) − 1/ln 10`, with
  `E₁` the energy at the end of the planned step, the smooth terms at `E_mid`
  and the clamp applied after averaging (`scattering_variance_birth`). With
  this U5 holds for all frozen step lengths (worst 1.5e-3 at x/R₁ = 0.9,
  5 mm steps).
- The hinge's second leg may be truncated by a voxel plane while the angle
  was sampled for the full step: an overestimate of scattering on boundary
  steps whose size is not yet measured. U5 and T8 do not probe it (they use
  uncut homogeneous steps or vary only the maximum step); the frozen T14
  (grid-size and half-voxel alignment refinement with angular and lateral
  observables, V3-003B) bounds it, and its negative control decides whether
  the default must switch to sampling for the truncated length.
- The step loop exists twice (reference loop and kernel). It is kept thin;
  trajectory-level parity between the reference and the float64 kernel is
  the guard.
- Process-based CPU parallelism costs 1–2 s of start-up per worker.

## Validation strategy (frozen before results exist)

Criteria are frozen in `validation/plans/v3-003-acceptance.md` (committed
with V3-003A before any LV/HR result is produced): shared-function parity
(decision 0001 classes), table round trip, scattering length against
Gottschalk's table (0.3 %), deterministic MCS integrator step independence
(2e-3) with a per-step-Highland negative control (≥ 5 %), RNG known-answer
and disjointness tests, trajectory parity (100 % discrete agreement, 1e-10
continuous), deterministic CSDA end depth, finite-slab escape, energy
balance (1e-5 float32, 1e-12 float64), R80 within 0.2 % of the CSDA range,
range straggling within 5 % of the Bohr quadrature and 10 % of Bortfeld's
fit, lateral variance within 2 % of Fermi–Eyges, transport step
independence (0.5 % in θ_rms), rotation invariance, DDA adversarial cases,
statistical parity across the three backends, partition invariance, and the
fail-closed configuration rules.

## Later validation outcome

- **2026-10-04 (V3-003B) — midpoint scoring falsified and replaced.** The
  frozen T9 step-independence probe exposed point-sampling aliasing of the
  midpoint deposit with the scoring-bin edges: with straggling and MCS off,
  150 MeV in water, 1 mm IDD bins, 2e5 histories on warp-cpu, the maximum
  IDD deviation from a 0.1 mm-step run was 62 % (125–140 mm) for s_max =
  1.0 mm, 32 % (plateau) for 0.33 mm, 80 % for 0.9 mm and 66 % for 2.0 mm
  steps with 2 mm bins; straggling smears but does not remove it (4.2 % at
  s_max = 1 mm), and the python reference shows the same (62.6 %). The bump at
  130–135 mm sits where the energy-loss step limit (f_E·R) first drops below
  s_max = 1 mm. The remedy is path-length-proportional apportioning of the
  step deposit along both hinge legs (exact for a uniform dE/dx within the
  step; the residual is second order in the step's energy change) with the
  new frozen T9-CI aliasing probe; the point-midpoint result is preserved
  here as contrary evidence. Also from V3-003B: the float32 per-batch
  accumulators of decision 0037 were replaced by int64 fixed-point
  accumulators (see decision 0037, amendment 2026-10-04).
- **2026-10-04 (V3-003B) — step-length rule replaced.** The rule
  `max_step_mm ≤ smallest scoring spacing` existed only to bound the
  per-leg walk of the path-length scoring (at most 8 pieces). It made the
  frozen T8 probe (s_max up to 5 mm with 0.2 mm lateral dose bins) impossible
  to run as frozen. It is replaced by a per-run computed bound
  `pieces = 3·⌈max_step/min_spacing⌉ + 4` (maximum over scoring grids; a leg meets up to ⌈·⌉+1 planes per axis including zero-length hops at a start exactly on a plane; cap
  4096 beyond which the configuration is rejected) carried in the kernel
  control structure; exceeding it at run time increments the fail-closed
  counter `scoring_pieces_overflow` and invalidates the result — the deposit
  is never silently truncated. Steps may now exceed the scoring spacing; the
  physics step limits (f_E, s_max, voxel planes) are unchanged.
- **2026-10-04 (V3-003B) — open determinism observation.** During the
  V3-003B work the exact partition-independence test of the python
  reference failed once: the undivided run contained an escaped history of
  about 10 MeV that did not occur in the two partial runs, and the failure
  could not be reproduced (30 reruns, cold Warp caches, 8 seeds × 240
  histories compared bitwise, 192k Python-scope Warp calls under heap
  churn). Separately, the GPU host runner twice aborted (SIGABRT) inside the
  Python-scope evaluation of shared Warp functions in the U1 parity test and
  passed on retry. Both point at the Python-scope `@wp.func` execution path
  as occasionally unreliable. Mitigations: a frozen repeatability probe (R1:
  identical runs must be bit-identical in tallies, counters, grids and
  per-history end states) is added to the CI/LV/HR suites, every reference
  result used as evidence must come from a run whose R1 passed, and the
  float64 Warp-CPU kernel — bit-identical to the reference in T1 — is the
  fallback reference executor if the observation recurs. The anomaly is
  preserved here as unexplained; it is not attributed to the physics.
- **2026-10-04 (V3-003B) — the Python-scope Warp evaluator is unreliable
  (reproduced).** A pool worker of the python reference backend died with
  SIGSEGV inside `inspect._bind` called from Warp `context.call_builtin`
  (chain: `types.__truediv__` ← `kinematics.tmax_mev` ← `em.bohr_variance`
  ← `reference._history`), reproduced in 1 of 10 test runs with several
  spawned workers and in a standalone two-worker repro, never in 23
  single-worker runs. Together with the three host SIGABRTs in the same path
  and the unexplained escaped history, this establishes that Warp 1.17's
  Python-scope execution of `@wp.func` builtins is subject to native
  memory corruption under multi-process use. Consequences: (a) the pool is
  fail-closed (a dead worker aborts the run, `TransportWorkerError`, no
  partial result), so no evidence in this task can come from a crashed run;
  (b) silent corruption is bounded by the frozen cross-checks that every
  qualification run must pass — T1 (trajectories bit-identical to the
  float64 kernel), R1 (bit-identical repeats) and the T12 python-vs-float64
  control; (c) the structural remedy — a reference executor that does not
  depend on Warp's Python-scope dispatch — was implemented in V3-003B on
  2026-10-05 (originally planned as V3-003C): `ionmc._wpfunc.python_twin`
  re-executes the source text of each `make_*` factory (`funcs.py`,
  `physics/em.py`, `physics/kinematics.py`) with `wp` replaced by the
  pure-Python shim `ionmc._pyshim` (math functions, a 3-vector class), so
  one source text defines the kernel functions and their twins; the
  reference backend uses only the twins, T1 (trajectories bit-identical to
  the float64 kernel, continuous columns 1e-10) still passes unchanged, and a
  test forbids any Warp `Function` call at Python scope on the reference
  path. The python backend became about two orders of magnitude faster
  (100 MeV, 200 histories, one worker, 0.5 mm steps, setup included: 0.54 to
  59 histories/s). Rule recorded in `docs/architecture/transport.md`: Warp
  functions run only inside kernels; every Python-side evaluation (reference
  backend, tests) uses the twins. The Python-scope tests (`test_formulas`,
  `test_tables_scattering`, `test_transport_scoring`, `test_transport_clip`,
  `test_straggling_gamma`, `test_shared_funcs`) were converted accordingly.
- **2026-10-04 (V3-003B) — T9 at full scale: straggling default and deposit
  apportioning changed.** The frozen T9 probe (1e6 histories, held-out seed
  base 20271004) failed at s_max = 0.1 mm for the clamped-Gaussian/Gamma
  sampler (IDD χ² 515 over 163 bins, max|z| 6.4, permutation p at the floor,
  ΔR80 0.026 mm) and still failed with a Gamma sampler for every ratio
  (χ² 447, max|z| 5.2, p 0.001) although ΔR80 fell to 0.010 mm; s_max 0.5 mm
  and f_E 0.005 passed in both cases. The Gamma sampler is adopted as the
  default (`bohr_gamma_v1`): for a common scale θ = σ²/mean the sum of the
  per-step losses is exactly Gamma-distributed with additive shape, so the
  full energy-loss distribution — not only its first two moments — is step
  independent wherever θ varies slowly along the path, which the R80 result
  confirms. The remaining deviation is deterministic and comes from
  apportioning the step deposit uniformly along the path while the stopping
  power changes by up to ≈ 1.6 % within a 2 %-energy-loss step near the
  Bragg peak (T9-CI deviations 5e-4 … 1.3e-3, above the ≈ 3e-4 sensitivity
  of the statistical T9). The deposit was therefore apportioned with a linear
  stopping-power ramp S(E₀) → S(E₁) along the hinge path (kept: residual
  second order in the curvature of S) — but the rerun at 1e6 histories still
  failed (χ² 427, max|z| 5.2; VAL-20261004-072852-16E681), so this was NOT
  the cause. The root cause, found with a calibration control (two seeds of
  the same configuration give χ² 150–160 over 163 bins, so the statistic is
  calibrated) and a parameter scan, is the short-step energy-loss branch:
  for `ρs < 10⁻³ R₀` the loss was taken as `S(E₀)·s`, which omits the
  first-order term `(s/2)·dS/dx` — about 2e-4 of the energy lost per step,
  used by 0.1 mm steps down to ≈ 115 MeV (≈ 55 mm depth) and never by 1 mm
  steps, shifting the range by ≈ 0.009 mm. With the branch threshold at
  10⁻⁹ or 10⁻⁵ the s_max 0.1 mm vs 1 mm comparison passes (χ² 163–185,
  p 0.25–0.57, ΔR80 ≤ 0.002 mm). The branch is therefore made
  second-order accurate (midpoint stopping power) and its threshold set to
  `10⁻² R₀`: below that the midpoint branch is more step-consistent than the
  telescoping form (deterministic 40 mm loss at 150 MeV, 0.1 vs 1 mm steps:
  2.5e-6 with the midpoint branch versus 3e-5 by telescoping; 0.01 vs 0.1 mm:
  2.6e-8 versus 3e-4), because the inverse-range table round-trip error of
  the telescoping accumulates per step; above it the telescoping is needed
  (the midpoint error grows with (s/R)²). The branch's original purpose
  (float32 cancellation in the telescoping) disappeared with the float64
  energy bookkeeping. The results before the change are
  preserved in the validation ledger (VAL-20261004-070915-E1A618,
  VAL-20261004-072852-16E681).

- **2026-10-06: T12 python-vs-cpu64 total-deposit failure at 9a27084 (seed
  base 20301004) and its remedy (V3-003B).** The qualification run
  (RUN-20261005T223540Z-773dc37b, 150 MeV protons in water, 0.2 mm bins) failed
  exactly one scalar of one pair: `total_deposit_mev`, python (40 x 100
  histories) 149.9999999998186 (se 1.80e-10) against warp-cpu float64
  (100 x 1e4) 149.99999999999835 (se 1.16e-11), relative difference 1.198e-12
  above the deterministic bound 1e-12. Diagnosis: not a bookkeeping bias. The
  python glue and the Warp kernel quantize identically (`floor(x/q + 1/2)`,
  q = 2^-30 MeV, residual tallied per grid, outside deposits in float64); the
  energy balance closes to <= 1e-12 for both (grid identity
  `in_grid + quantization + outside = step_deposit + cutoff`). The bare IDD
  grid sum used by `t12_observables` omits the tallied rounding residual, a
  zero-mean error of about 1.1e-8 MeV per history (per-history sd from the
  batch spread: python 1.14e-8, cpu64 1.155e-8) that differs between samples.
  The python mean deficit (-1.81e-10) is 1.0 combined standard errors, the
  cpu64 mean (-1.7e-12) is 0.14 sigma; the previous head's python sample
  (hr-90986b4-b, seed base 20291004) shows the same size of scatter
  (per-history sd 8.7e-9) with no preferred sign. The "deterministic" bound of
  1e-12 relative (1.5e-10 MeV) is below the sampling noise of a 4000-history
  python sample (1.8e-10 MeV), so it could fail by chance in either direction.
  Remedy: the degenerate-scalar rule of `t12_compare` accepts a difference up
  to `bound * scale + 3.5 * hypot(se_a, se_b)`: the deterministic precision
  bound for the systematic part plus the rounding-noise allowance (recorded as
  `noise_allowance`); a systematic offset still fails. No kernel or shared
  physics function was changed, the T1 trace stays bit-identical. A regression
  test shows that both float64 backends close the energy balance to 1e-12 and
  that the per-primary deposit (grid + quantization residual) agrees to 1e-12
  for identical seeds, and a T12 test covers the allowance.

To be appended from committed result files.
