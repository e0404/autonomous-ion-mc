# 0017 — Non-water materials on the 3-D multiple-scattering path

- Status: accepted
- Date: 2026-09-10
- Task: DEV-012
- Affects: architecture, physical accuracy, validation strategy

## Problem

DEV-011 (decision `0016`) brought the 3-D multiple-Coulomb-scattering path to
1-D voxelized density heterogeneity but only for **water** (a homogeneous
`is_water_only` guard). The depth-dose path already handles per-voxel *materials*
(DEV-010, decision `0015`). This task lifts the guard and makes the scattering
path material-heterogeneous, so a patient-like phantom (bone, adipose, lung,
tissue) transports in 3-D with both the correct water-equivalent range and the
correct lateral spread. It mirrors DEV-010 for the scattering kernel.

## Context and evidence

- **Stopping vs scattering use different densities.** As on the depth-dose path,
  the energy loss uses the **water-equivalent density** `rho_we = SPR(material)
  rho_phys` (stopping-power ratio; decision `0015`) with water's `<Z/A>` for the
  straggling. But multiple Coulomb scattering scales differently: the Highland
  scattering power uses the areal thickness of the *physical* material and the
  *material* radiation length, `theta0^2 = (13.6/pv)^2 (rho_phys s / X0_mat)`. So
  a material voxel needs **two** densities on the scattering path: `rho_we` for
  the energy loss and `rho_phys` (with `X0_mat`) for the MCS. This is exactly why
  DEV-011 deferred materials.
- **The radiation lengths are already in the library.** Every tissue `Material`
  carries a Tsai radiation length (decision `0015`), validated against PDG; only
  the per-voxel wiring is new.

## Decision

1. **Per-voxel physical density and radiation length.** The engine builds, and
   merges alongside the existing per-voxel water-equivalent and oxygen-equivalent
   densities, a per-voxel **physical density** and **radiation length**. The
   voxel merge now collapses consecutive voxels only when *all* their per-voxel
   physics quantities match, keeping every array aligned.

2. **Scattering uses the physical density and material X0.** The scattering
   reference driver and Warp kernel look up, by depth: `rho_we` for the energy-
   loss step and straggling (water `<Z/A>`, water-equivalent frame), and
   `rho_phys` + `X0_mat` for the Highland scattering power. The step is limited to
   the voxel boundary in depth exactly as in DEV-011.

3. **The guard is lifted.** `run_scattering` now accepts any 1-D voxelized
   material `VoxelSlab`; it only requires every voxel material to have a known
   radiation length (`> 0`). A uniform single-material slab still collapses to
   one voxel, so water (and any homogeneous material) reproduces the prior result
   bit-for-bit.

4. **Material-aware Fermi-Eyges oracle.** The oracle gains
   `lateral_sigma_x_material_mm`, which drives the energy-vs-depth from the
   integrated water-equivalent thickness while using the local physical density
   and material radiation length for the scattering power, giving a quantitative
   `sigma_x` reference for material interfaces.

## Consequences

- The depth-dose and scattering paths now have the **same** per-voxel material
  capability; the only remaining Stage-3 items are 3-D geometry with arbitrary
  incidence and decoupled scoring grids.
- The scattering drivers carry three per-voxel density-like arrays
  (water-equivalent, physical, radiation length) plus the boundaries.

## Acceptance targets (validation `v3_material_scattering.py`)

- **Water/homogeneous equivalence**: a water `VoxelSlab` reproduces the
  `WaterSlab` scattering result bit-for-bit; a uniform bone slab reproduces a
  homogeneous run.
- **Material lateral spread**: a homogeneous bone slab reproduces the material-
  aware Fermi-Eyges `sigma_x` (physical density and bone X0) within a few percent,
  and its R80 lands at `R_water/WER`.
- **Material interface**: a water/bone/water phantom reproduces the material-aware
  Fermi-Eyges `sigma_x` across the interface within a few percent, with energy
  conserved.
- **Cross-backend**: reference vs Warp CPU and (host) CUDA agree on `sigma_x`
  (tight) and the depth dose (3-D float32 budget) across a material interface.
