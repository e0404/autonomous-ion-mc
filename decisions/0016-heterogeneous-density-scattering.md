# 0016 — Density-heterogeneous 3-D multiple-Coulomb-scattering transport

- Status: accepted
- Date: 2026-09-10
- Task: DEV-011
- Affects: architecture, physical accuracy, validation strategy

## Problem

The 3-D transport path with multiple Coulomb scattering (MCS; decision `0011`,
`run_scattering`) is homogeneous: it uses one scalar mass density and one
radiation length for the whole slab, and (decisions `0014`/`0015`) it explicitly
**rejects** a heterogeneous or non-water `VoxelSlab`. Meanwhile the 1-D depth-
dose path is fully voxelized (density DEV-009, materials DEV-010). This task
brings the 3-D scattering path to **1-D voxelized density heterogeneity** (a
`VoxelSlab` of water at varying density along the beam axis), so the lateral
spread responds to the local density and the peak lands at the correct water-
equivalent depth, matching the depth-dose path's capability. It mirrors DEV-009
for the scattering path and is validated by the lateral spread through density
interfaces.

Non-water *materials* on the scattering path (which require the physical density
and the material radiation length for the MCS, distinct from the water-equivalent
density used for stopping) are the DEV-012 step and stay guarded here.

## Context and evidence

- **MCS scales with the local density.** The Highland projected RMS angle over a
  step uses the areal thickness `x = rho * s` and the radiation length:
  `theta0^2 = (13.6/pv)^2 (rho s / X0)` (scattering-power form, decision `0011`).
  For water at varying density only `rho` varies (`X0` and `<Z/A>` are constant),
  so the per-voxel density that already drives the depth-dose step (DEV-009) also
  drives the scattering power — no new physics, only the per-voxel lookup.
- **Depth is the voxel axis.** A `VoxelSlab` is 1-D along `z`; the scattering
  history moves in 3-D with direction `d`, advancing in depth by `s d_z`. The
  voxel a history occupies is set by its depth `pz`, and (forward transport,
  `d_z > 0`) that voxel index is non-decreasing. A step is limited so `s d_z`
  does not cross a voxel boundary, mirroring DEV-009's boundary limiting for the
  depth-dose path (now projected onto `z`).
- **Water-only keeps the frame simple.** For water the water-equivalent density
  equals the physical density, so the stopping, straggling and MCS all use the
  same per-voxel density and water's `X0`/`<Z/A>`; the DEV-010 water-equivalent
  machinery is a no-op here and is reused unchanged.

## Decision

1. **Per-voxel density on the scattering path.** The reference scattering driver
   and the Warp scattering kernel look up the current voxel's density by depth
   `pz`, limit the step so `s d_z` stays within the voxel (`s <= (z_boundary -
   pz)/d_z`, in addition to the geometry limit), and advance a non-decreasing
   voxel index as boundaries are crossed. The density drives the energy-loss
   step, the straggling and the Highland scattering power identically to the
   homogeneous case; `X0` and `<Z/A>` stay water's.

2. **The guard is narrowed, not removed.** `run_scattering` now accepts a
   density-heterogeneous **water** `VoxelSlab` (any per-voxel densities) and
   still rejects a slab containing a **non-water** material (any voxel whose
   water-equivalent density differs from its physical density), with a clear
   `NotImplementedError` — non-water materials on the scattering path are DEV-012.

3. **Homogeneous equivalence is preserved.** A uniform-density water `VoxelSlab`
   (and a `WaterSlab`) collapses to a single voxel (the DEV-009 merge), so the
   scattering result is bit-for-bit identical to the pre-DEV-011 homogeneous
   path.

4. **Heterogeneous Fermi-Eyges oracle.** The Fermi-Eyges lateral-spread oracle
   (`ionmc.physics.fermi_eyges`) gains a piecewise-density variant that
   integrates the scattering power with the local density and the water-
   equivalent energy-vs-depth, so the layered-phantom lateral `sigma_x(z)` has a
   quantitative reference across a density interface.

## Consequences

- The depth-dose and scattering paths now share the same 1-D voxel-density
  capability; only per-voxel *materials* remain unequal (depth-dose has them,
  scattering does not yet).
- The scattering drivers gain the density/boundary arrays and a per-history
  voxel index, exactly as the depth-dose drivers did in DEV-009.

## Acceptance targets (validation `v3_scattering_heterogeneity.py`)

- **Homogeneous equivalence**: a uniform-density water `VoxelSlab` reproduces the
  `WaterSlab` scattering result (reference bit-exact, Warp within the decision-
  `0001` cumulative tolerance).
- **Density scaling of lateral spread**: a uniform slab of density `rho`
  reproduces the Fermi-Eyges lateral `sigma_x` computed at that density (the peak
  at `R_water/rho`, the lateral spread scaling accordingly) within a few percent.
- **Layered interface**: in a water/dense/water phantom the lateral `sigma_x(z)`
  matches the piecewise-density Fermi-Eyges oracle within a few percent, and
  energy is conserved.
- **Cross-backend**: reference vs Warp CPU and (host) CUDA agree on the depth
  dose and `sigma_x` across a density interface within the decision-`0001`
  tolerances.
