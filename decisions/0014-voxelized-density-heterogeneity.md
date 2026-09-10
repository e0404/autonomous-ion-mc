# 0014 — 1-D voxelized density heterogeneity and water-equivalent transport

- Status: accepted
- Date: 2026-09-10
- Task: DEV-009
- Affects: architecture, physical accuracy, validation strategy, scientific interpretation

## Problem

Stages 1 and 2 transport protons through a single homogeneous water slab
(`WaterSlab`, one mass density). Real treatment planning runs on a **voxelized,
heterogeneous** patient geometry: each voxel has its own mass density (from CT)
and material, and the dose depends on the **water-equivalent thickness (WET)**
the beam traverses, not the geometric depth. This task opens **milestone V3** by
making the depth-dose transport **density-heterogeneous** along the beam axis: a
1-D stack of voxels each with its own mass density, with the stopping power,
energy-loss straggling and nonelastic nuclear rate all responding to the local
density. It validates the two most fundamental heterogeneity behaviours: WET
equivalence (density scaling) and layered-interface behaviour.

Full 3-D voxel geometry with arbitrary beam incidence, ray/voxel traversal and
transforms, per-voxel *material* (composition, not just density) via stopping-
power ratios, and scoring grids decoupled from the transport grid are the
remaining Stage-3 content and are **deferred** to later Stage-3 tasks. This task
is deliberately confined to the existing 1-D depth-dose axis so the change is
the geometry and the density-dependence, not a transport rewrite.

## Context and evidence

- **Mass stopping power and density.** The tabulated stopping power is a *mass*
  stopping power `S/rho` [MeV cm^2/g]; the linear energy loss over a path `dl`
  is `S/rho * rho * dl`. The physics functions already take `density` as an
  argument (`energy_loss_step_length`, `midpoint_energy_loss`,
  `bohr_straggling_sigma`), so density-dependence needs no physics change, only
  a **per-voxel density** fed to each step. The CSDA range in mm therefore
  scales as `1/rho`: a slab of density `rho` shifts the Bragg peak to
  `R_water / rho`, and the water-equivalent depth is `integral rho dl`.
- **Nonelastic nuclear rate and density.** The macroscopic nonelastic cross
  section `Sigma = n_O sigma(E)` scales with the oxygen *number* density, which
  for a fixed composition scales with mass density. So `Sigma(E, voxel) =
  Sigma_water(E) * rho[voxel]` for water at density `rho`.
- **Electronic `<Z/A>` (straggling prefactor)** is a per-gram, material-intrinsic
  quantity; for a single composition it is the same in every voxel and only the
  density scales the linear straggling. (Per-voxel `<Z/A>` and stopping-power
  ratios for *different compositions* are the DEV-010 material-heterogeneity
  step.)
- **Forward transport is monotonic in depth.** The primary (and the forward-
  emitted secondaries) move in `+z` with `dl > 0` and no backscatter, so the
  voxel index a history occupies is **non-decreasing**. A step is limited to the
  current voxel boundary so its density is unambiguous, and the voxel index is
  simply advanced when the step reaches a boundary — O(1) per step, no search in
  the inner loop.

## Decision

1. **1-D voxel geometry.** A new `VoxelSlab` describes the beam axis as an
   ascending array of voxel boundaries `z_boundaries_mm` (length `n_vox + 1`,
   `[0] = 0`, `[-1] = depth`) and a per-voxel `density_g_per_cm3` (length
   `n_vox`). It is built from a list of `(thickness_mm, density)` layers (exact
   interface positions, so interface behaviour is not smeared by discretisation)
   or from a uniform grid. All voxels share one `Material` for DEV-009 (water at
   varying density); per-voxel material is the DEV-010 extension.

2. **Unified voxel representation.** The engine always transports through a voxel
   array. A homogeneous `WaterSlab` is represented as a **single voxel** spanning
   `[0, depth]` at its density, so the voxel-boundary step limiting is a no-op
   and DEV-004..008 results are reproduced **exactly** (backward compatible).

3. **Per-step density (reference and Warp).** Each step: look up the current
   voxel's density; limit the step length to `min(dl, z_boundary[voxel+1] - z)`
   (in addition to the existing energy-loss and geometry limits); take the
   energy-loss, straggling and nuclear step with that density; after advancing
   `z`, advance the voxel index while `z` has reached the next boundary. The
   nonelastic rate uses a per-voxel oxygen number density `n_O[voxel] = AVOGADRO
   * atoms_per_gram('O') * density[voxel]`.

4. **Water-equivalent scoring is unchanged.** Dose is still scored per geometric
   depth bin (decoupled scoring grids are a later Stage-3 task). The WET
   behaviour is validated through the shifted Bragg-peak position, not by
   re-binning.

5. **Deferred:** 3-D voxel geometry, arbitrary beam incidence and ray/voxel
   traversal, per-voxel *material* composition via stopping-power ratios and
   per-voxel `<Z/A>`/radiation length, nuclear removal on the 3-D scattering
   path in heterogeneous media, and scoring grids decoupled from the transport
   grid. Voxel-boundary step limiting slightly increases the straggling step
   granularity near interfaces (a physically required consequence of the density
   change); it is a no-op for a homogeneous slab.

## Consequences

- `TransportEngine` accepts either a `WaterSlab` or a `VoxelSlab`; internally it
  builds the voxel density/boundary arrays and per-voxel oxygen densities once.
  The reference driver and the Warp kernel gain the density and boundary arrays
  and a per-history voxel index.
- Backward compatibility is a validation gate: a `VoxelSlab` of uniform water at
  density 1.0 must reproduce the homogeneous `WaterSlab` depth dose bit-for-bit.
- The peak position now scales with `1/rho`; the absolute range in mm is no
  longer the CSDA range but `R_csda / rho` (the WET range).
- The 3-D scattering path (`run_scattering`) is homogeneous only (it uses the
  front-voxel density). To avoid a silently wrong result, it **rejects** a
  heterogeneous `VoxelSlab` (`n_voxels > 1`) with a `NotImplementedError`;
  density heterogeneity on the scattering path is a later Stage-3 task.

## Acceptance targets (validation `v3_density_heterogeneity.py`)

- **Homogeneous equivalence**: a uniform `VoxelSlab` (density 1.0) reproduces the
  `WaterSlab` depth dose to round-off (reference) and within the decision-`0001`
  cumulative tolerance (Warp).
- **WET / density scaling**: a uniform slab of density `rho` (e.g. 0.5, 1.2)
  shifts the R80 depth to `R80_water / rho` within a few tenths of a percent, at
  150 and 200 MeV.
- **Layered interface**: in a two-layer phantom (e.g. water / bone-density /
  water) the peak lands at the depth predicted by the integrated water-equivalent
  thickness, and energy is conserved across the interface.
- **Energy conservation**: `deposited + escaped = energy_in` holds per the DEV-007/
  008 gates in the heterogeneous geometry.
- **Cross-backend**: reference vs Warp CPU and (host) CUDA agree cumulatively
  within the decision-`0001` tolerances for a heterogeneous phantom.
