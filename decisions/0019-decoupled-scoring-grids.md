# 0019 — Scoring grids decoupled from the transport grid

- Status: accepted
- Date: 2026-09-11
- Task: DEV-014
- Affects: architecture, validation strategy

## Problem

`REQUIREMENTS.md` requires **changeable scoring grids** — the scoring grid must
be selectable independently of the transport grid, in resolution *and*
alignment. Milestone **V3** closes on four gates; three are met (WET equivalence,
layered interface, rotated-vs-axis-aligned, decisions `0014`–`0018`). The
remaining gate is **grid-independence of integral dose under scoring-grid
changes**. The scoring grids (`DepthDoseGrid`, `DepthLateralGrid`) were already
resolution-independent of the transport voxel grid (deposition splits each step
by depth overlap), but they were **hard-anchored** at depth `z = 0` and lateral
`x = 0`, so alignment could not vary.

## Decision

1. **Arbitrary grid alignment.** `DepthDoseGrid` gains `origin_mm` (default 0);
   its bins span `[origin_mm, origin_mm + depth_mm]`. `DepthLateralGrid` gains
   `depth_origin_mm` and `lateral_center_mm` (both default 0); its depth bins
   span `[depth_origin_mm, depth_origin_mm + depth_mm]` and its lateral bins
   `[lateral_center_mm ± half_width_mm]`. The defaults reproduce the historical
   grids bit-for-bit.

2. **Origin-aware deposition on all paths.** The reference deposition helpers
   (`_deposit_along_step`, `_deposit_zx`) and the Warp kernels
   (`csda_depth_dose_kernel`, `csda_scattering_kernel`) index the bin as
   `floor((z − z_origin) / dz)` with bin edges `z_origin + b·dz`, and the lateral
   bin as `floor((x − x_lo) / dx)` with `x_lo = lateral_center − half_width`. This
   is a strict generalisation: `floor((z − z_origin)/dz) = b` guarantees the bin
   edge `z_origin + (b+1)·dz > z` for any `z` and any sign of `b`, so the
   depth-overlap loop remains well-formed (monotone, no negative segment) for a
   grid that starts before, at, or after the entrance.

3. **Scoring grid stays a separate object.** The transport grid (the voxel
   geometry, `VoxelSlab`/`WaterSlab`) and the scoring grid (`DepthDoseGrid`/
   `DepthLateralGrid`) are already distinct arguments; DEV-014 only makes the
   scoring grid's *alignment* free, completing the decoupling. Energy deposited
   outside the scoring grid is dropped (the existing in-range guard).

## Validation (milestone V3 closure, grid-independence)

Pre-registered before running (decision `0001` methodology), in
`validation/v3_scoring_grid.py` on the host runner (CPU + CUDA):

- **integral_dose_resolution** — a straight-ray (scattering off) run scored on a
  coarse and a 4×-finer grid of the same extent: total deposited energy equal to
  round-off, and the finer per-bin dose summed 4:1 equals the coarse per-bin dose
  (`≤ 1e-12`), i.e. energy per depth interval is grid-independent.
- **integral_dose_alignment** — the same run scored on a grid shifted in depth
  origin and lateral centre (still containing the dose): total deposited energy
  invariant to round-off (`≤ 1e-12` reference, `≤ 1e-6` Warp float32).
- **lateral_origin** — a *discriminating* check of the lateral origin (which the
  translation-invariant `σ_x` cannot show): with scattering off the beam stays on
  the `x = 0` axis, so a lateral window centred on the beam captures the dose, one
  shifted off it captures **none**, and on a wide shifted grid the energy-weighted
  mean lateral position is the beam axis (`|mean_x| < 0.2 mm`), not the grid
  centre.
- **lateral_shift_sigma_x** — with scattering on, `σ_x(z)` is invariant under a
  lateral-centre shift, `≤ 1e-6 mm` (an invariance property, complementing the
  discriminating `lateral_origin` check).
- **partial_coverage** — a grid starting past the entrance captures strictly less
  energy (the depth origin truly shifts the scored window).
- **cross-backend** — reference vs Warp CPU vs CUDA on a shifted grid:
  depth-dose cumulative `≤ 5e-4`, `σ_x` `≤ 0.05 mm` (decisions `0001`/`0011`).

This closes milestone **V3** (voxelized heterogeneous geometry and materials).

## Consequences

- The scoring grid is fully decoupled from the transport grid in resolution and
  alignment; a planning-style scoring window can be placed anywhere along and
  across the beam. Beam-frame scoring (decision `0018`) composes with this: the
  origin/centre are in the scoring frame the transport already scores into.
- **Not yet (Stage-3 capability beyond the V3 gates):** a full 3-D voxel grid
  with arbitrary per-voxel material maps and robust ray/voxel (Siddon/DDA)
  traversal for true patient geometries — a follow-on task and a Stage-4
  prerequisite. The 2-D marginal scoring grid is sufficient for the V3
  grid-independence gate; a full 3-D scoring grid arrives with Stage-4
  influence-matrix scoring.
