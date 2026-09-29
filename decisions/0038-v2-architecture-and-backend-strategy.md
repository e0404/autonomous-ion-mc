# 0038 — V2 architecture: shared physics functions, three execution paths, planning-oriented data flow

- Status: accepted (revisable by later decisions)
- Date: 2026-09-30
- Task: V2-001
- Affects: architecture, physical accuracy, numerical accuracy, reproducibility, performance, validation strategy

## Problem

The protocol requires a Python-first package with a reference Python path,
Warp CPU and Warp CUDA paths sharing one physical model, readable physics,
GPU-oriented execution without per-particle Python loops, and planning-oriented
repeated transport (beamlets, influence matrices). The architecture must be
fixed before implementation so that later tasks compose.

## Candidate approaches

1. **Independent pure-Python physics plus separately written Warp kernels.**
   Strongest oracle independence, but two implementations of every formula
   that can silently diverge — exactly what the protocol warns against.
2. **Warp-only physics, reference path = Warp CPU.** Violates the requirement
   for a native Python reference path and offers no float64 oracle.
3. **Shared physics functions written once as `@wp.func`, executed by (a) the
   CPython interpreter from Python scope for the reference path and (b) Warp
   generic kernels specialized for float32 and float64 (selected).** The
   reference path keeps an *independent transport loop* (per-particle Python,
   `numpy.random.Generator` streams, float64 through `wp.float64` scalars),
   so kernel scheduling, indexing, compaction, atomics and RNG are checked
   independently while the physics formulas remain single-sourced.

## Selected design

- **Layers.** `ionmc.data` (external/tabulated datasets with provenance
  cache), `ionmc.physics` (analytic models and table construction; pure
  functions), `ionmc.transport` (condensed-history stepping expressed as
  shared `wp.func` functions plus backend drivers), `ionmc.geometry`
  (homogeneous and voxel phantoms with materials/density), `ionmc.sources`
  (pencil beams, beamlets, arbitrary direction), `ionmc.scoring` (dose,
  energy deposition, LET, species, fluence, external lookup tallies;
  transport grid and independent scoring grids), `ionmc.results`
  (self-describing persisted outputs), `ionmc.cli`.
- **Transport model class.** Class-II-like condensed history for charged
  ions: continuous electronic energy loss from per-(species, material) tables
  on a logarithmic kinetic-energy-per-nucleon grid, Gaussian energy-loss
  straggling with a documented validity domain, Gaussian multiple Coulomb
  scattering, discrete nuclear interactions sampled from tabulated total
  cross sections with species-specific secondary/fragment production,
  explicit transport of charged secondaries, and explicit accounting (never
  silent local deposition) of energy carried by untransported neutral
  particles or particles leaving the phantom. Detailed model choices are
  separate decisions with their own evidence.
- **Data source separation.** Transport only sees tables (stopping power,
  range, cross sections, straggling parameters) on fixed grids; whether a
  table came from an analytic formula or an external dataset is a property
  of the table's provenance record, not of the kernels.
- **Execution.** Particle-per-thread kernels over a state-of-arrays buffer;
  secondaries are pushed to a device stack and processed in later passes;
  finished particles are compacted out; a step loop runs entirely on the
  device with a bounded number of steps per launch. Immutable physics tables
  and geometry stay resident across beamlets and batches.
- **Precision.** float32 kernels by default with float64 accumulators for
  scored quantities where measurable bias is found; the same kernels are
  launchable in float64 for falsification (V2-NUM).
- **Random numbers.** Warp counter-based streams keyed by (seed, batch,
  particle index); the reference path uses independent numpy streams with the
  same keying so batches are statistically reproducible, not bitwise.
- **Uncertainty.** Batch statistics (independent batches with distinct
  streams); per-beamlet uncertainty for influence matrices.
- **Units and coordinates.** Lengths in mm, energies in MeV (kinetic energy
  per nucleon for ions where stated), densities in g/cm³, dose in Gy per
  primary; right-handed Cartesian grid with explicit voxel edge coordinates;
  beam direction given as a unit vector; all persisted with the results.
- **Outputs.** Machine-readable NumPy archives plus JSON metadata and
  plain-text/SVG human-readable summaries, written without third-party
  plotting libraries so the clean-install workflow needs only the two runtime
  dependencies.

## Expected tradeoffs

- Sharing formulas gives up an independently derived oracle for the physics
  functions themselves; that gap is closed by analytic/tabulated tests of each
  function against external evidence, not by backend agreement.
- Python-scope execution of shared functions is slow (order 10⁴–10⁵ steps/s),
  which is adequate for validation cases of 10³–10⁴ histories only.
- float32 transcendental rounding in Python scope when plain floats are used
  is avoided by always driving the reference path with `wp.float64` scalars.

## Validation strategy

Each layer is tested in isolation (analytic and tabulated expectations), the
reference path is compared statistically to both Warp backends
(NUM-BACKEND-PARITY), and physics is validated against independent evidence
per species (PHY-* suites in `validation/release-plan.json`).
