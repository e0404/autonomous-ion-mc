# 0020 — 3-D voxel grid with ray/voxel DDA traversal

- Status: accepted
- Date: 2026-09-11
- Task: DEV-015
- Affects: architecture, physical accuracy, validation strategy

## Problem

Through milestone V3 the geometry is a **1-D voxel stack** (`VoxelSlab`): the
scattering path looks up the voxel by a scalar material coordinate `u = n̂·p`
along a single slab normal. `REQUIREMENTS.md` (*Voxelized geometries*, *robust
ray/voxel traversal with arbitrary incidence*) requires a true 3-D voxel grid so
a beam can transport through a patient-like `Nx×Ny×Nz` geometry. The V3 gates did
not need it (they used slab/layered phantoms); this task adds the general 3-D
geometry as the remaining Stage-3 capability and a Stage-4 prerequisite.

## Decision

1. **`VoxelGrid3D`: a lab-axis-aligned 3-D grid.** Per-voxel mass density on a
   uniform per-axis spacing, min corner at `origin_mm`, a single `material`
   composition (density-only first cut; a per-voxel material map is a pure
   data-plumbing extension of the same flat arrays and is deferred). Flat index
   convention (shared by the reference driver and the Warp kernel):
   `flat = (i·Ny + j)·Nz + k`.

2. **Amanatides–Woo voxel DDA, per-step recompute** (not Siddon). The random-hinge
   step changes direction every step and the step is already physics-bounded
   (fractional energy loss), so a whole-ray Siddon integral is the wrong
   amortisation; the incremental array-free Siddon *is* a DDA anyway. The DDA is
   phrased in the existing beam-frame idiom: the single 1-D coordinate `u = n̂·p`
   and its limit `(voxel_z[v+1]−u)/mproj` generalise to **three lab-axis
   coordinates** `u_k = p0[k] + mᵏ·q` with `mᵏ` = row `k` of the beam-frame
   rotation `R = (e1|e2|d̂)` (the 1-D path already computes exactly the `z` row),
   `rateₖ = mᵏ·D_beam`, and a per-axis distance-to-next-face
   `s_k = (o_k + (i_k+step_k)·h_k − u_k)/rateₖ`. The step is clipped to the
   nearest face: `s = min(s_phys, s_x, s_y, s_z, s_box)`. Geometry state (the
   voxel indices, the per-axis limits) is **recomputed from the position and
   direction at the start of each step** rather than carried incrementally, so no
   cross-step float state can drift between the float64 reference and the float32
   kernel, and the 1-D reduction stays exact-to-round-off.

3. **The grid is not rotated.** An oblique beam through an axis-aligned grid
   already exercises full 3-D traversal (all three `rateₖ ≠ 0`); a grid rotation
   is, by rigid-motion equivalence, the beam rotation already supported by
   decision 0018, adds no traversal coverage, and is deferred.

4. **Robustness.** Axis-parallel rays (`|rateₖ| ≤ eps`) take `s_k = +∞`
   (a large sentinel in the kernel), so the index never changes and no spurious
   clip occurs — the signed 3-D case replaces the 1-D `max(mproj, 1e-6)` floor,
   which is only valid for forward motion. On-face degeneracy is resolved by a
   direction-aware nudge of `s_tie = 1e-4 mm` (far above float32 coordinate
   resolution, far below a voxel) with a fixed axis tie-priority `x<y<z`,
   identical in the reference and kernel. Escape (status 2) when a step reaches
   the grid bounding box or a re-seeded index leaves `[0, N_k)`.

5. **Parallel drivers; the 1-D path is untouched.** `run_scattering` dispatches on
   geometry type to a sibling reference driver and Warp kernel; the
   `VoxelSlab`/`WaterSlab` code is unchanged and therefore bit-exact. Beam-frame
   marginal scoring (decisions 0018/0019) is unchanged; a 3-D lab-frame scoring
   volume is deferred (the traversal already yields the lab position each step, so
   a future 3-D scorer needs only the deposition mapping). The pure depth-dose
   `run` rejects `VoxelGrid3D` (it is longitudinal, +z only).

## Validation (`validation/v3_voxel_grid_3d.py`, pre-registered per decision 0001)

- **grid_reduction_to_slab** — a single-column `(1,1,Nz)` grid built from a
  uniform-spacing `VoxelSlab` with distinct adjacent densities, deterministic
  (scattering+straggling off), reproduces the slab depth dose to round-off
  (`max|Δcumsum|/ΣE ≤ 1e-9`, R80 `≤ 1e-4 mm`). The x/y limits are `+∞` and the
  z-arithmetic mirrors the 1-D path, so the two drivers agree bit-for-bit up to
  the float reconstruction of the uniform faces.
- **homogeneous_box_vs_slab** — a homogeneous 3-D box reproduces `WaterSlab`
  (deterministic: R80 `≤ 1e-2 mm`, total E `≤ 1e-4`; stochastic: depth-dose
  cumulative `≤ 3e-3`, σ_x within SEM).
- **oblique_wet** — an oblique beam through an off-axis dense insert matches an
  independent Siddon central-ray `∫ρ dl` oracle (`≤ 1e-3` relative WET).
- **energy_conservation** — deposited + escaped − in `≤ 1e-9` (reference) /
  `≤ 1e-5` (Warp), with an escaped-energy accumulator on the grid path.
- **cross-backend** — reference vs Warp CPU vs CUDA on a scattering grid config
  (statistical parity per decision 0001). The tight metrics are σ_x `≤ 0.05 mm`
  (physics) and the **same-precision** CPU-vs-CUDA depth-dose cumulative
  `≤ 5e-4`. The reference(float64)-vs-Warp(float32) depth-dose cumulative gets a
  DDA-specific bound `≤ 3e-3`: unlike the 1-D path (which clips only at a few
  material interfaces), the DDA clips at *every* voxel face (~160 z-crossings plus
  lateral crossings per history), so float32-vs-float64 near-corner face-decision
  flips decorrelate proportionally more histories and accumulate ~1e-3 into the
  depth dose. A deterministic (scattering-off) grid run stays tight
  (`~3e-6`, no face flips), confirming this is float32 DDA divergence, not a bug.

## Provenance

Amanatides & Woo (Eurographics 1987, the voxel DDA); Siddon (Med. Phys. 12(2),
1985, the exact radiological path — used only as the offline WET oracle); Jacobs
et al. (1998, array-free incremental Siddon = a DDA). All are published
algorithms implemented from their descriptions, no third-party code, no external
data.

## Deferred

Per-voxel material map, grid rotation, a 3-D lab-frame scoring volume,
non-uniform per-axis spacing, and CT/Hounsfield ingestion.
