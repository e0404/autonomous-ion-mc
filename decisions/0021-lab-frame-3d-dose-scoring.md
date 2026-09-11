# 0021 — Lab-frame 3-D dose scoring on the voxel-grid path

- Status: accepted
- Date: 2026-09-11
- Task: DEV-016
- Affects: architecture, validation strategy

## Problem

Scoring so far is a **beam-frame 2-D marginal** (`DepthLateralGrid`: depth × one
transverse axis, summed over the other). That is sufficient for depth-dose and
lateral-spread validation, but Stage 4 (treatment-planning scoring, dose-influence
matrices) needs **volumetric dose on a lab-frame 3-D grid**: dose per voxel in a
fixed patient frame, so beamlets entering at different positions/directions
accumulate into a common grid. The 3-D voxel-grid transport path (decision 0020)
already computes the lab position `x = p0 + R·q` every step, so a 3-D scorer needs
only a deposition mapping, not new geometry. This is the prerequisite deferred by
DEV-015 and the bridge to Stage 4.

## Decision

1. **`DoseGrid3D`: a lab-frame scoring grid.** An `Nx×Ny×Nz` uniform grid (origin,
   spacing) accumulating deposited energy per voxel [MeV]; independent of the
   transport grid (its own resolution/alignment). `voxel_volume_cm3` and a
   density-scaled dose conversion are exposed.

2. **Scored on the 3-D voxel-grid transport path** (`run_scattering` with a
   `VoxelGrid3D` geometry) via an optional `dose_grid` argument. When given, each
   step additionally deposits its energy `w·de` at the **lab midpoint** of the
   step (`x = p0 + R·q_mid`) into the containing dose voxel — one accumulation per
   step, which conserves energy exactly for a beam contained in the dose grid.
   The energy-cut terminal deposit likewise lands at the stop position. The
   beam-frame marginal scoring is unchanged and still produced.

3. **Point-at-midpoint deposition** (not a scoring DDA) for this first cut: the
   physics step (~1 mm) is typically ≤ a dose voxel, so point deposition into the
   containing voxel is accurate to the O(step/voxel) discretization while
   conserving energy per step. A path-length-splitting scoring DDA is a later
   refinement.

## Validation (`validation/v4_dose3d.py`, pre-registered per decision 0001)

- **energy_conservation** — a beam contained in the dose grid: the summed 3-D
  dose equals the deposited energy to round-off (reference) / ≤ 1e-5 (Warp).
- **depth_dose_consistency** — the 3-D dose projected onto the beam axis (summed
  over the transverse voxels) reproduces the beam-frame depth dose within the
  point-vs-overlap discretization (cumulative ≤ 1e-2) and R80 ≤ one dose-voxel.
- **grid_independence** — the total 3-D dose is invariant to the dose grid's
  resolution and alignment (it only rebins the same deposits).
- **cross-backend** — reference vs Warp CPU vs CUDA 3-D dose: total ≤ 1e-5, and
  the per-voxel dose agrees within the DDA float32 face-flip budget (statistical,
  decision 0020).

## Consequences

- Volumetric lab-frame dose is available on the patient-geometry path, unblocking
  Stage-4 beamlet-resolved scoring and influence matrices (the next task carries
  beamlet identity into per-beamlet `DoseGrid3D` buffers and assembles a sparse
  matrix).
- **Deferred:** a 3-D scorer on the 1-D-geometry (`WaterSlab`/`VoxelSlab`) path,
  path-length-splitting deposition, and dose (energy/mass) unit conversion driven
  by per-voxel density (the grid currently scores energy per voxel; the material
  density for dose is available from the geometry when needed).
