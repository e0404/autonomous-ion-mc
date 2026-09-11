# 0024 — Batch-based statistical uncertainty for 3-D dose (and LET_d)

- Status: accepted
- Date: 2026-09-11
- Task: DEV-019
- Affects: scoring, uncertainty quantification, validation strategy

## Problem

`REQUIREMENTS.md` makes **statistical uncertainty estimation** a MUST (number of
histories, expectation, and *variance, standard error, or another justified
measure*), with **batch-based estimation** and **planning-aware uncertainty**
(quantifiable for beamlet-/influence-matrix calculations) as SHOULDs. The V4
milestone likewise phrases the beamlet-sum gate "within statistics". The 1-D
depth-dose path already has batch uncertainty (`run_batched` /
`BatchedDepthDoseResult`, decision `0010`), but the lab-frame 3-D dose grid
(decision `0021`) and the LET_d scorer (decision `0023`) have none. This decision
adds per-voxel batch uncertainty for the 3-D scored quantities.

## Decision

1. **Independent-batch estimator (same convention as decision 0010).**
   `TransportEngine.run_scattering_batched` runs `n_batches` independent scattering
   runs of `n_histories // n_batches` histories each, batch `b` seeded `seed+1+b`
   (independent MC samples), each accumulating a `DoseGrid3D`. The per-voxel mean
   dose is the mean of the batch dose grids and the **standard error of the mean**
   is `std(batches, ddof=1) / sqrt(n_batches)`. This reuses the shared
   `_scatter_state` dispatch, so it works identically on the reference and Warp
   CPU/CUDA paths and needs no kernel change.

2. **`BatchedDoseResult`.** Carries `mean_dose3d_mev` and `standard_error_mev`
   (both `(nx,ny,nz)`), `n_batches`, `histories_per_batch`, `path`, `device`, and
   — when `score_let` is set — `mean_let_d_kev_um` and `standard_error_let_kev_um`
   (the per-voxel mean and SEM of the batch LET_d ratios). It exposes
   `relative_standard_error` (per-voxel SEM/mean where mean>0) and
   `mean_relative_uncertainty(min_dose_frac)` — the dose-threshold-restricted mean
   relative SEM, the standard Monte-Carlo quality metric (uncertainty is reported
   only where there is appreciable dose; the low-dose tail is statistically noisy).

3. **1/sqrt(N) scaling is the correctness anchor.** For a fixed batch count, the
   high-dose-region mean relative SEM scales as `1/sqrt(N_histories)`: quadrupling
   the histories halves it (within statistical tolerance). This is the pre-
   registered validation of the estimator.

4. **Scope.** Single-source 3-D dose and LET_d uncertainty. Beamlet-/influence-
   resolved (planning-aware) uncertainty — batching `run_scattering_multi` /
   `assemble_influence_matrix` — reuses this same estimator and is the deferred
   next step.

## Validation (`validation/v4_uncertainty.py`, pre-registered per decision 0001)

- **se_scaling** — the high-dose-region mean relative SEM at `4N` histories is
  ~half that at `N` (ratio in `[0.35, 0.71]`, i.e. `1/sqrt(4)` within a generous
  statistical band), confirming the `1/sqrt(N)` law. The mean is taken over a
  **shared** high-dose voxel mask (from the higher-statistics run) so the two runs
  average over the identical voxel set — otherwise per-voxel SEM-estimate noise and
  differing self-masks would blur the comparison.
- **mean_energy_conservation** — the batch-mean total dose equals the per-history
  deposited energy for a contained beam (mean over batches conserves energy).
- **uncertainty_sanity** — the relative SEM is positive in high-dose voxels and
  zero where there is no dose; SEM ≤ mean in the high-dose region.
- **cross-backend** — the reference and Warp CPU batched mean/SEM agree in the
  high-dose region (the estimator is a deterministic function of the per-batch
  dose grids, which agree per the dose scorer's cross-backend budget, decision
  `0022`).

## Consequences

- Every 3-D dose (and LET_d) result can be accompanied by a per-voxel standard
  error and a high-dose-region relative-uncertainty summary, satisfying the MUST
  statistical-uncertainty requirement for the volumetric scorers.
- **Deferred:** beamlet-/influence-resolved uncertainty (planning-aware), a
  running/online single-pass variance estimator (currently the batches are stored
  and reduced), and the statistical (gamma/uncertainty-based) cross-backend
  spatial dose comparison that this uncertainty machinery now enables
  (decision `0022`).
