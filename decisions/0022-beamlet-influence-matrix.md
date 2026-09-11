# 0022 — Beamlet-resolved scoring and sparse dose-influence matrices

- Status: accepted
- Date: 2026-09-11
- Task: DEV-017
- Affects: architecture, treatment-planning capability, validation strategy

## Problem

Ion-therapy treatment planning optimises the weights of thousands of **beamlets**
(pencil-beam spots), which requires a **dose-influence matrix**: one row per
beamlet, one column per dose voxel, giving each beamlet's dose distribution. This
is the core Stage-4 capability (`REQUIREMENTS.md`: *Beamlet-resolved scoring*,
*Dose influence matrices*, *Sparse influence representations*, *Batched beamlet
execution*). It builds on the lab-frame 3-D dose scorer (decision `0021`) and the
beamlet identity already carried in the particle state.

## Decision

1. **Batched multi-beamlet transport.** `TransportEngine.run_scattering_multi`
   transports a list of beamlet sources **together** in one launch: each beamlet
   `i` is sampled with seed `seed + i` (its own per-history RNG streams and its
   `beamlet` id), the states are concatenated, and the immutable physics/material
   data is **not** reinitialised. `run_scattering`/`run_scattering_multi` now share
   a common `_scatter_state` dispatch, so single- and multi-beamlet runs use the
   identical transport.

2. **Sparse influence matrix (`SparseInfluenceMatrix`).** A `(n_beamlets,
   n_voxels)` matrix in CSR form (`indptr`/`indices`/`data`, dose in MeV; flat
   voxel index `(i·ny + j)·nz + k`), because a pencil beam illuminates only a small
   fraction of the patient. It exports to a documented `.npz` and provides
   `to_dense`, `beamlet_dose_flat`, and `total_dose` (the broad-field dose).

3. **Assembly by per-beamlet scoring, thresholded (`assemble_influence_matrix`).**
   Each beamlet is transported with a `DoseGrid3D` scorer using seed `seed + i`
   (the *same* partition `run_scattering_multi` uses), its 3-D dose flattened, and
   voxels below `threshold_frac` of that beamlet's peak dropped from its row. This
   is the "chunked-dense per-beamlet buffers thresholded into a sparse matrix"
   first implementation from the roadmap; a GPU hash-table keyed by (voxel,
   beamlet) is a later throughput option.

4. **Exact broad-field additivity.** Because each beamlet's separate run and its
   contribution to `run_scattering_multi` use the identical seed-`(seed+i)`
   histories and dose is a linear sum of per-history deposits, the sum of the
   per-beamlet rows equals the batched broad-field dose to round-off — no
   statistical tolerance needed for the partitioning check.

## Validation (`validation/v4_influence.py`, pre-registered per decision 0001)

- **beamlet_sum_equals_broadfield** — `SparseInfluenceMatrix.total_dose()` (before
  thresholding) equals `run_scattering_multi`'s dose to round-off
  (`max|Δ|/ΣE ≤ 1e-9`, reference) — the first V4 gate.
- **sparse_vs_dense** — a thresholded matrix keeps ≥ 99 % of the total energy at a
  1 % peak threshold, and its `to_dense()` reproduces the retained per-beamlet
  doses exactly.
- **energy_conservation** — the summed influence dose equals the total deposited
  energy for a contained set of beamlets.
- **cross-backend** — the tight, discriminating spatial parity is **deterministic**
  (scattering/straggling off): the broad-field dose agrees **per voxel** to the
  float32 budget (total `≤ 1e-5`, per-voxel `≤ 5e-3` of peak) both
  reference-vs-Warp-CPU *and* Warp-CPU-vs-CUDA — the straight-line trajectory is
  not chaotic, so this isolates the float32 arithmetic and certifies the kernel is
  spatially identical across backends. Under **scattering on**, per-voxel dose is
  *not* a tight cross-backend metric for *any* pair: float32 CPU and float32 CUDA
  arithmetic is not bit-identical (FMA contraction, transcendental
  implementations) and, like the reference(float64)-vs-float32 case, DDA face-flip
  decorrelation amplifies it into two independent MC estimates (measured
  CPU-vs-CUDA per-voxel `~3e-2` of peak at `N=2e4`/beamlet, recorded as a
  diagnostic). Only the **total** energy stays tight (CPU-vs-CUDA `~1e-10`,
  reference-vs-float32 `~5e-9`), because dose is a linear sum of per-history
  deposits over the same seed partition — that total is the gated stochastic
  quantity. A genuine statistical (gamma/uncertainty-based) spatial comparison of
  two MC estimates is a deferred cross-cutting validation item.

## Consequences

- A sparse dose-influence matrix can be assembled and exported for inverse
  planning; batched multi-beamlet transport reuses immutable data.
- **Deferred:** a GPU (voxel, beamlet) hash-table assembly, per-beamlet scoring
  directly in the kernel (currently one transport launch per beamlet), LET/fluence
  and species-/energy-resolved scorers, and beamlet-resolved uncertainty — the
  remaining Stage-4 items.
