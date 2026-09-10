# 0015 — Per-voxel tissue materials via stopping-power ratios

- Status: accepted
- Date: 2026-09-10
- Task: DEV-010
- Affects: architecture, physical accuracy, validation strategy, scientific interpretation, data (material library)

## Problem

DEV-009 (decision `0014`) made the depth-dose transport density-heterogeneous,
but every voxel is the same *material* (water) at a different mass density. Real
patient geometry is heterogeneous in **composition**: cortical bone, adipose,
muscle, lung and air stop protons differently *per gram* than water, because
their electron density per gram (`<Z/A>`) and mean excitation energy `I` differ.
Treating bone as water at bone's density overestimates its stopping (bone stops
~12 % less per gram than water), so the water-equivalent thickness (WET) is
wrong. This task adds a **tissue material library** and makes the transport use
each voxel's material through a **stopping-power ratio (SPR)**, opening the
composition axis of Stage 3.

## Context and evidence

Research (Claude physics-researcher subagent, 2026-09-10; ICRU 44/49, NIST,
Schneider 2000, PDG) and an in-repo SPR computation established:

- **Tissue compositions and I-values** (ICRU-44 / NIST): cortical bone
  (rho 1.92, I 110 eV), adipose (0.95, 64.8), soft tissue (1.03, 72.3), skeletal
  muscle (1.05, 75.3), lung tissue (1.05, 75.3), air (1.205e-3, 85.7). The
  `Material` class already carries mass fractions, a provenance-tagged `I`,
  `<Z/A>` and a radiation length.
- **Stopping-power ratio.** The mass-stopping-power ratio material:water,
  `SPR_mass(E) = (S/rho)_mat / (S/rho)_wat`, is the electron-density ratio
  `<Z/A>_mat/<Z/A>_wat` times a slowly varying logarithmic `I`-value term. The
  code already has a full Bethe `AnalyticStoppingPower` bound to a `Material`, so
  the SPR is evaluated **as the ratio of the analytic mass stopping powers**
  (material / water), which carries the shell/Barkas/Bloch corrections
  consistently and cancels model error. In-repo SPR at 150 MeV: bone 0.884,
  adipose 1.019, soft tissue 0.995, muscle 0.990, lung 0.991. The **linear
  water-equivalent ratio** `WER = SPR_mass * (rho_mat/rho_wat)` is what shifts the
  Bragg peak: bone 1.70, adipose 0.968, soft tissue 1.025, muscle 1.040, lung
  1.041 — all inside the published Schneider/ICRU bands.
- **Energy dependence of the SPR** is small: soft tissues vary < 0.25 % over
  10-250 MeV, adipose ~0.9 %, cortical bone ~2.5 % (its `I` is farthest from
  water). A single SPR evaluated at a mid-therapeutic reference energy is
  accurate to < 0.25 % for soft tissue and ~1-1.5 % for bone at the entrance.
- **Straggling.** The Bohr prefactor is `<Z/A>` (electrons per gram), a
  composition-only quantity; density enters separately. Because `SPR_mass ~
  <Z/A>_mat/<Z/A>_wat` up to the small `I`-term, transporting a voxel as *water
  at its water-equivalent density* (`SPR*rho`) with water's `<Z/A>` reproduces the
  material's linear straggling to the `I`-value log term (a sub-percent effect on
  an already-small peak-width contribution).
- **Radiation length** per material (Tsai's formula; PDG): bone 26.99 g/cm^2 (the
  only tissue far from water's 36.08), adipose 41.23, soft tissue 36.78, muscle
  36.81, lung 36.53, air 36.62. Used by multiple Coulomb scattering, which is on
  the 3-D path only; populated in the library for when heterogeneous scattering
  exists.
- **Nonelastic nuclear rate.** Oxygen-only removal (decision `0012`) is wrong
  once composition varies: it underestimates the nonelastic rate by ~2x in bone
  (Ca, P, C) and ~3.4x in adipose (carbon-rich, oxygen-poor). The rate must sum
  over the significant target nuclei.

## Decision

1. **Tissue material library.** Add `CORTICAL_BONE`, `ADIPOSE`, `SOFT_TISSUE`,
   `SKELETAL_MUSCLE`, `LUNG_TISSUE`, and `AIR` to `ionmc.materials` with the
   ICRU-44/NIST compositions, densities, I-values and Tsai radiation lengths,
   each provenance-tagged. Water's `I = 75 eV` (ICRU 49 / PSTAR) is the SPR
   denominator, to stay consistent with the water table (decision `0006`).

2. **Water-equivalent density per voxel.** The engine computes each voxel's
   **water-equivalent density** `rho_we = SPR_mass(material, E_ref) * rho_phys`,
   with `SPR_mass` the ratio of the analytic Bethe mass stopping powers
   (material / water) at a reference energy `E_ref = 150 MeV`. Transport then
   reuses the DEV-009 per-voxel-density machinery **unchanged**: the water table
   times `rho_we` gives the correct material stopping, and (per the evidence
   above) the material's linear straggling to the `I`-value term. The transport
   code is therefore untouched; only the per-voxel density array the engine
   builds changes.

3. **Composition-scaled nuclear rate.** The per-voxel nuclear "oxygen-equivalent"
   number density becomes `rho_phys * N_A * sum_{i in C,N,O,P,Ca,...} (w_i/A_i)
   (A_i/A_O)^{2/3}` (geometric `A^{2/3}` scaling of the single ICRU-63 oxygen
   cross-section shape; hydrogen excluded, decision `0012`). For water this
   reduces to `n_O`, so DEV-007/008 are reproduced exactly. This removes the
   ~2x (bone) / ~3.4x (adipose) systematic error of oxygen-only removal at the
   cost of one per-material scalar.

4. **Per-voxel material geometry.** `VoxelSlab` gains an optional per-voxel
   material (from `from_material_layers([(thickness, Material), ...])`); the
   existing density-only constructors keep a single material. Consecutive voxels
   are merged only when both their water-equivalent density *and* their nuclear
   density match, so a uniform single-material slab still collapses to one voxel
   (homogeneous equivalence preserved).

5. **Scalar SPR is the first-cut approximation.** `E_ref = 150 MeV` gives a single
   SPR per material (bone ~1-1.5 % entrance error, soft tissue < 0.25 %). The
   **energy-dependent SPR** (a per-material SPR(E) lookup or per-material stopping
   table, the standard fast-MC choice) is the immediate refinement and is
   deferred with a documented accuracy note.

## Deferred (with justification)

- Energy-dependent SPR / per-material stopping tables (bone's 2.5 % energy
  dependence; soft tissues are already < 0.25 % with the scalar).
- Element-specific measured nonelastic cross sections (C, Ca) and per-element
  secondary spectra (the `A^{2/3}` scaling is a ~10-15 % interim on the nuclear
  rate, far better than the 2-3.4x error of oxygen-only).
- Per-tissue Sternheimer density-effect parameters (< 0.1 % for protons under
  900 MeV; `density_effect=None`).
- Multiple scattering on the 3-D path for **non-water** materials (and
  heterogeneous slabs). The scattering path feeds the *physical* density into the
  water stopping table with no SPR, so it is correct only for water; it now
  rejects both a heterogeneous `VoxelSlab` and a homogeneous non-water material
  (its water-equivalent density differs from the physical density) with a
  `NotImplementedError`, rather than silently returning a wrong Bragg depth.
  Radiation lengths are populated now for that future work and validated against
  PDG.

## Consequences

- A bone voxel at 1.92 g/cm^3 now transports with `rho_we ~ 1.70` (its WER), not
  1.92, so the peak lands at the correct water-equivalent depth. The nuclear
  pre-peak attenuation in bone/adipose is ~2-3x larger than under oxygen-only.
- The transport drivers and Warp kernel are unchanged from DEV-009; the material
  physics enters entirely through the per-voxel density and nuclear-density
  arrays the engine assembles once at construction.
- `_merge_voxels` merges on the transported (water-equivalent density, nuclear
  density) pair rather than physical density alone.

## Acceptance targets (validation `v3_material_composition.py`)

- **Water-equivalent ratio per material**: a slab of each tissue at its reference
  density puts R80 at `R80_water / WER` within a few tenths of a percent
  (bone ~1.70, adipose ~0.97, muscle ~1.04, lung ~1.04, at 150/200 MeV), where
  `WER` is the analytic `SPR_mass * rho`.
- **Water back-compatibility**: a `VoxelSlab` of water reproduces the DEV-009
  homogeneous result exactly (SPR(water) = 1, nuclear scale = `n_O`).
- **Material interface**: a water/bone/water phantom puts the peak at the depth
  predicted by the integrated `sum WER_i * t_i`, and energy is conserved across
  the interface with nuclear + secondary transport.
- **Nuclear composition scaling**: the per-voxel nuclear rate in bone/adipose is
  larger than oxygen-only by the composition factor (bone ~2x, adipose ~3.4x),
  and reduces to `n_O` for water.
- **Cross-backend**: reference vs Warp CPU and (host) CUDA agree cumulatively for
  a multi-material phantom within the decision-`0001` tolerances.
