# 0025 — Beamlet-resolved (planning-aware) statistical uncertainty

- Status: accepted
- Date: 2026-09-11
- Task: DEV-020
- Affects: treatment-planning capability, uncertainty quantification, validation

## Problem

`REQUIREMENTS.md` asks (SHOULD) that uncertainty be **planning-aware** —
"quantifiable for beamlet-resolved or influence-matrix calculations" — and the V4
milestone phrases the beamlet-sum gate "within statistics". DEV-017 (decision
`0022`) assembles a sparse dose-influence matrix and DEV-019 (decision `0024`)
adds batch-based per-voxel dose uncertainty for a single source. This decision
combines them: a **per-beamlet standard error** on the influence-matrix rows.

## Decision

1. **Batched per-beamlet assembly (`assemble_influence_matrix_batched`).** For
   each beamlet `i`, transport `n_batches` independent history batches with
   `run_scattering_batched` (decision `0024`) seeded from `seed + i*n_batches` (so
   the per-batch seeds `seed + i*n_batches + 1 + b` never collide across beamlets),
   giving the per-voxel **mean dose** and its **standard error of the mean**. The
   row's stored value (`data`) is the mean dose; a parallel `data_sigma` stores the
   matching per-voxel SEM. Thresholding is on the **mean** dose (`threshold_frac`
   of the beamlet's peak mean), and `data_sigma` is carried for exactly the kept
   voxels, so it stays aligned with `data`.

2. **`SparseInfluenceMatrix.data_sigma`.** An optional `(nnz,)` array parallel to
   `data` (None for the exact single-run matrix of decision `0022`). Adds
   `beamlet_sigma_flat(row)` (dense per-voxel SEM), `beamlet_relative_uncertainty(
   row, min_dose_frac)` (the beamlet's high-dose-region mean relative SEM — the
   standard planning quality metric), and `total_sigma()` (the broad-field SEM,
   `sqrt(Σ_b σ_b²)` per voxel, since the beamlets are independent). `save`/`load`
   round-trip `data_sigma` when present.

3. **The exact single-run matrix is unchanged.** `assemble_influence_matrix`
   (decision `0022`, `data_sigma=None`) keeps its exact same-seed-partition
   additivity (`Σ` beamlet rows `==` `run_scattering_multi`, DEV-017). The batched
   matrix is a separate, statistically-independent estimate — its beamlet-sum
   agreement with the broad field is now a *statistical* check (**within
   statistics**), the V4-milestone phrasing.

## Validation (`validation/v4_beamlet_uncertainty.py`, pre-registered)

- **per_beamlet_se_scaling** — a representative beamlet's high-dose mean relative
  SEM scales as `1/√N` (4× histories → ~half, ratio in `[0.35, 0.71]`, shared
  mask), inheriting decision `0024`'s estimator check per beamlet.
- **sigma_alignment** — `data_sigma` has the same length as `data`; every kept
  voxel has a defined (finite, ≥ 0) SEM; `beamlet_sigma_flat` is 0 off-row.
- **beamlet_sum_within_statistics** — the sum of the independent-batch beamlet
  mean doses agrees with an independent broad-field batched mean **within the
  combined statistical uncertainty** (per-voxel `|Δ| ≲ k·σ` in the high-dose
  region), the V4 "within statistics" gate (complementing DEV-017's exact
  same-partition round-off check).
- **cross-backend** — the **deterministic** (scattering off) per-beamlet mean dose
  agrees reference-vs-Warp-CPU (and CUDA-vs-CPU) per voxel to the float32 budget.

## Consequences

- The influence matrix can carry a per-beamlet, per-voxel dose uncertainty and a
  per-beamlet quality metric, satisfying the planning-aware-uncertainty SHOULD and
  letting the V4 beamlet-sum gate be stated with statistics.
- **Deferred:** beamlet-resolved LET_d / fluence uncertainty, a running online
  variance (batches are stored and reduced), and correlated-uncertainty
  propagation into an optimiser.
