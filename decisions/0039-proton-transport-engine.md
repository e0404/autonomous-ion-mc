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
| Mean energy loss | S(E)·s; range-table inversion | inverse-range telescoping `Rinv(R0) − Rinv(R0 − ρs)` with a linear `S·s` branch when `ρs < 10⁻³ R0` | exact in CSDA, step-size independent, monotone; the short branch avoids float32 cancellation (archived precision measurement) |
| Straggling | Gaussian only; Vavilov; Gaussian/Gamma/uniform switch (`G4IonFluctuations`); two-moment Gaussian/Gamma | Bohr variance, sampled with a two-moment-preserving rule: Gaussian clamped to [0, 2·mean] (no resampling) for mean/σ ≥ 3, otherwise Gamma with shape k = (mean/σ)² and scale σ²/mean (Marsaglia–Tsang, `u^(1/k)` boost for k < 1; rejection bounded by 64 attempts, exhaustion increments the fail-closed counter `straggling_rejection`) | the mean is exact and the variance is exact on the Gamma branch; on the Gaussian branch the clamp (P(x<0) = 0.13 % at mean/σ = 3) removes less than 1 % of the variance; the `G4IonFluctuations` uniform branch has variance mean²/3 whatever the Bohr value and the truncated Gaussian near mean/σ = 2 loses variance, so they are not used; sum over steps is Gaussian by the central limit theorem; delta electrons are deposited locally so single-step tails do not reach dose |
| Multiple scattering | per-step Highland; Molière; differential Highland (Kanematsu); differential Molière (Gottschalk) | **differential Molière** scattering power `T_dM = f_dM(pv, p₁v₁)(E_s/pv)²/X_S`, `E_s = 15.0 MeV`, `f_dM = 0.5244 + 0.1975 lg(1−(pv/p₁v₁)²) + 0.2320 lg(pv) − 0.0098 lg(pv) lg(1−(pv/p₁v₁)²)` (clamped ≥ 0), scattering length `1/(ρX_S) = α N_A r_e² (Z²/A){2 ln(33219 (AZ)^{-1/3}) − 1}` Bragg-additive, applied as a Gaussian polar angle with a random hinge | step-size independent by construction (per-step Highland is not; a negative-control test proves the instrument has power); ranked best against Hanson theory and measurement in the source paper; the research report's label "differential Highland" for these coefficients was wrong and is corrected |
| Lateral displacement | explicit correlated sampling; random hinge | random hinge (move `a·s`, deflect, move `(1−a)·s`, `a` uniform) | reproduces the Fermi–Eyges second moments exactly without extra draws |
| Geometry traversal | `floor(p)` with nudge; incremental DDA | incremental DDA with the voxel index as state and plane snapping on crossing | the nudge variant stalled in the archived toy kernel |
| Scoring | split steps at scoring planes; midpoint deposit | midpoint deposit with `max_step ≤ min scoring spacing` enforced (otherwise the configuration is rejected) | simplest unbiased choice at the enforced step size; grid refinement/misalignment is a falsification probe in V3-011 |
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
6. Deposit `E₀ − E₁` at the hinge-path midpoint into every scoring grid
   (outside-grid deposits are tallied separately).
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
  x/R₁ = 0.05 with 1 mm steps). The birth step therefore uses the exact step
  average of the logarithmic term, `lg(1 − (pv(E₁)/p₁v₁)²) − 1/ln 10`, with
  `E₁` the energy at the end of the planned step, the smooth terms at `E_mid`
  and the clamp applied after averaging (`scattering_variance_birth`). With
  this U5 holds for all frozen step lengths (worst 1.5e-3 at x/R₁ = 0.9,
  5 mm steps).
- The hinge's second leg may be truncated by a voxel plane while the angle
  was sampled for the full step: a small overestimate of scattering on
  boundary steps, measured by the step-independence probe rather than
  corrected ad hoc.
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

To be appended from committed result files.
