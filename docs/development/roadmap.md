# Development roadmap

- Status: living document, maintained by the lead orchestrator
- Established: 2026-09-10 (task `DEV-001`, kickoff prompt `experiment/prompts/kickoff-v1.md`)
- Governing constraints: `EXPERIMENT.md`, `REQUIREMENTS.md`
- Related decisions: `decisions/0002-staged-development-roadmap.md`,
  `decisions/0003-out-of-worktree-development-via-staging-ref.md`

This roadmap identifies the major capability stages needed to reach the
research-grade, treatment-planning-oriented ion Monte Carlo system defined by
`REQUIREMENTS.md`, the dependencies between those stages, and the validation
milestone that closes each stage. It deliberately does **not** freeze detailed
architecture or physics-model choices; those are made per stage, recorded under
`decisions/`, and revised when evidence requires it.

Everything below the *Starting state* section describes **intended** work. A
stage is only "implemented" when this document says so explicitly and the
corresponding validation record exists for a committed SHA.

## Starting state (2026-09-10, `develop` at `ce13397`)

Implemented and integrated before this roadmap:

- the autonomous-development infrastructure (task worktrees, controlled
  commit/push/PR/merge tools, exact-SHA local validation gate, controlled
  host/GPU runner, Codex worker, human-intervention channel, telemetry);
- lightweight GitHub CI (pre-commit with hygiene hooks only, pytest,
  documentation build);
- environment, Warp-backend and repository-capability diagnostics under
  `infrastructure/diagnostics/`;
- decision `0001` fixing the Warp CPU-vs-CUDA parity methodology for
  deterministic kernels.

Not yet implemented: any Monte Carlo source package (`src/` is empty), any
physics, any validation of physics, any benchmark, any formatter, linter or
type checker in pre-commit.

Host execution environment observed at kickoff (host-runner probe
`RUN-20260910T001222Z-5ec48d9a` for task `DEV-001`, plus the diagnostics
recorded in `docs/warp_backend_diagnostics.md`): NVIDIA RTX A6000 (sm_86),
CUDA 12.9, Warp 1.17.0, numpy 2.5.3, Python 3.12.3 in the controlled
host-runner environment, which has no network access; 32 logical CPU cores.

Kickoff research findings that shape the stages below (sources are cited in
decision `0002` and will be re-cited in the per-stage decisions that act on
them):

- Warp's documentation states that kernels launched on the CPU device execute
  serially. The Warp CPU path is therefore a correctness and portability path;
  throughput targets apply to CUDA.
- Warp's own source shows that a `@wp.func`-decorated function called from
  Python scope falls back to the original Python function, so physics
  functions can be shared verbatim between the Warp kernels and a Python
  driver loop. Kernel bodies themselves (thread indexing, atomics) need a
  Python re-expression, and Warp's random-number builtins likely need a
  bit-exact pure-Python mirror. Both points are to be settled empirically in
  Stage 0.
- Warp's random-number generator is a counter-seeded PCG-style hash
  (`rand_init(seed, offset)`), which supports per-history streams that are
  identical on CPU and CUDA; Warp 1.15+ offers a deterministic execution mode
  for atomic accumulation.
- Published GPU proton/ion engines (gPMC, MOQUI, FRED, the Mayo code, goCMC)
  converge on one history per thread with an in-kernel loop over a small
  fixed-size secondary stack, float32 arithmetic with batch-wise
  accumulation, and beamlet-tagged sparse scoring; a first correct
  implementation is expected in the 10^5–10^6 protons/s band on one GPU with
  1–3×10^6 protons/s as the state-of-the-art target.
- NIST PSTAR/ASTAR tables are copyrighted Standard Reference Data and not
  freely redistributable; MCsquare publishes PSTAR-derived per-material
  stopping-power tables under Apache-2.0 at immutable Git revisions. The
  PSTAR tables use a water mean excitation energy of 75 eV; ICRU 90
  recommends 78 ± 2 eV, a difference of about one percent in stopping power.

## Architectural principles fixed by the experiment constraints

These are not roadmap decisions; they are inherited constraints that every stage
must respect (see `EXPERIMENT.md`, *Implementation Constraints*):

1. **Python-first.** The package is a Python package with a native Python API;
   principal transport physics stays readable at Python source level.
2. **One physical model, three execution paths.** Reference Python, Warp CPU
   and Warp CUDA must share the physical model. The reference path is the
   correctness oracle, "independent or minimally transformed".
3. **Backend-oriented separation** of physical modeling from execution
   (*Backend-extensible design*). Checked at the end of Stage 1 by
   demonstrating that the transport physics module has no import-time or
   call-time dependency on Warp when driven by the reference path.
4. **Treatment-planning orientation** (*Ion-therapy treatment planning
   orientation*). Repeated beamlet transport and beamlet-resolved scoring
   (dose influence matrices) drive architecture and performance decisions,
   not general-purpose MC flexibility. Stage 4 is where this becomes a
   validated capability; Stages 1–3 must not preclude it (beamlet identity in
   the particle state, immutable physics data separated from per-run data).
5. **Correctness before speed.** Performance is primary, but every
   optimization that changes numerical behavior must be validated.

## Rules applied to every validation milestone

- The decision record that opens a stage fixes, **before** any comparison is
  executed: the compared quantities, the metric, the acceptance tolerance,
  the reference dataset (with license and provenance) and the estimator or
  convention definitions on both sides. Tolerances are never chosen after
  seeing the comparison result (decision `0001` explains why).
- Decision `0001` covers only deterministic kernel-level parity between Warp
  CPU and Warp CUDA. Equivalence between the reference Python path and the
  Warp paths, and statistical consistency of stochastic transport output, are
  separate criteria that must be established by their own decisions when
  first needed (Stage 0 and Stage 1 respectively).
- Every validation result is recorded for the exact committed SHA with
  `record_local_validation`; scripts that must run on the host runner depend
  only on numpy, Warp and the standard library.

## Capability stages

Stages are ordered by dependency, not strictly by time; later stages may start
tasks once their prerequisites are integrated. Each stage names the
requirements it primarily serves (by `REQUIREMENTS.md` heading), its main
dependencies, and the validation milestone that closes it.

### Stage 0 — Package foundation, execution model and physics-data layer

Serves: *Python-first implementation*, *Reference Python execution*, *NVIDIA
Warp acceleration*, *Analytical physics layer*, *External/tabulated data
layer*, *Separation of transport from data source*, *External Data
Management* (all MUSTs and the SHOULDs *Configurable cache location* and
*Automated data preparation*), *Units and Conventions*, *Explicit random
control* (seeding design), *Pre-commit framework*, *Static typing*,
*Documentation*.

Content, split into two bounded tasks:

**`DEV-002` — package, conventions, execution model, analytic stopping
power.**

- installable `ionmc` package (`src/` layout, `pyproject.toml`; CI installs
  and tests it); pre-commit extended with a formatter, a linter and a type
  checker, mirrored in CI;
- explicit internal unit system and coordinate conventions;
- the shared-source execution model: how one Python source of a physics
  function runs as reference Python and as a Warp function/kernel on CPU and
  CUDA, and how the reference path is kept a faithful oracle. This is settled
  by numerical experiments on the host runner (Python-scope calls of
  `@wp.func` functions, Python-scope behaviour of `wp.rand_init`/`wp.randf`,
  float64 atomics on multi-dimensional arrays) and must include at least one
  kernel with data-dependent control flow and per-thread state (an iterative
  range integration and a per-thread random-number loop), so that the model
  is stressed on the features transport will need rather than only on
  elementwise evaluation. The resulting decision is explicitly **provisional
  until Stage 1** re-examines it under real transport;
- the reference-Python versus Warp equivalence criterion (including how a
  float64 reference is compared with float32 kernels) recorded as a new
  decision before the first comparison is run;
- analytical electronic stopping power for protons: Bethe formula with shell,
  density, Barkas and Bloch corrections and documented low-energy handling,
  with the mean excitation energy as an explicit, provenance-tagged
  per-material parameter (75 eV reproduces ICRU 49/PSTAR, 78 eV is the
  ICRU 90 value); CSDA range by integration;
- counter-based, per-history random-number scheme decided and documented,
  including a bit-exact pure-Python mirror of the Warp generator if Python
  scope cannot execute the Warp builtins (exercised fully only in Stage 1).

**`DEV-003` — data layer and tabulated stopping power.**

- versioned, checksummed, locally cached external physics datasets with
  configurable cache location, provenance metadata, automated preparation and
  offline reuse. The data layer separates an *acquire* phase (needs network;
  writes an immutable, checksummed cache entry) from a *use* phase (offline,
  cache-only, fails with an actionable message when data are missing),
  because neither the orchestrator sandbox nor the host runner is guaranteed
  network access;
- redistribution licensing of the tabulated source resolved and recorded
  before it is wired in (the Apache-2.0 MCsquare PSTAR-derived tables are the
  leading candidate; the NIST tables themselves are not redistributed);
- tabulated proton stopping power for water exposed through the same
  interface as the analytical layer, so transport later does not care which
  layer supplies the quantity.

Contingency for the network constraint: acquisition is exercised where
network exists (GitHub CI, the Codex worker on the host). If no environment
available to the autonomous system can fetch the dataset, a small reference
table for water (tens of kilobytes, Apache-2.0 provenance) may be committed
as test/validation reference data with the justification `REQUIREMENTS.md`
requires; a full multi-material dataset is never committed.

Depends on: nothing beyond the infrastructure baseline.

Validation milestone **V0** (closed by `DEV-003`):

- analytical stopping power versus an I-value-consistent reference table
  (75 eV against PSTAR-derived data) with tolerances fixed in the `DEV-002`
  decision before the comparison; the ICRU 90 (78 eV) offset is reported as a
  separately quantified model difference, not absorbed into the tolerance;
- CSDA ranges of protons in water at least at 100, 150 and 200 MeV versus the
  reference table, again with pre-fixed tolerances;
- identical results on Warp CPU and Warp CUDA under the decision `0001`
  criterion for the deterministic evaluation kernels, and agreement of both
  with the reference Python path under the new Python-versus-Warp criterion,
  on the experiment workstation GPU.

### Stage 1 — Proton condensed-history transport in homogeneous water

Serves: *Ion transport* (protons first), *Homogeneous reference geometries*,
*Energy deposition*, *Absorbed dose* (energy deposition divided by voxel mass
from the material density), *Statistical uncertainty estimation*,
*Batch-based estimation*, *Backend consistency*, *Parallel-safe random
streams*, *Deterministic backend reproducibility*, *Cross-backend statistical
consistency*, *GPU-oriented architecture*, *Numerical convergence*,
*Configurable accuracy/performance tradeoffs*, *Backend-extensible design*
(checked here, see principle 3).

Content:

- particle state layout (structure-of-arrays for global data, register-
  resident per-history state, species and beamlet identity fields present
  from the start), history batching, and the transport loop organization
  suited to Warp: one history per thread with an in-kernel loop over a small
  fixed-size secondary stack, a hard per-history step-count limit with a
  truncation counter, and no per-particle Python loop on accelerated
  backends; a deferred queue for rare expensive branches is kept as a later
  optimisation option;
- condensed-history step with continuous energy loss (stopping-power
  integration, possibly via a precomputed fractional-energy-loss table),
  energy-loss straggling (Gaussian with Bohr variance as the baseline used by
  the clinical fast engines MCsquare, MOQUI and FRED; Vavilov/Landau as a
  configurable alternative), and multiple Coulomb scattering (Gaussian with
  Highland/Rossi–Greisen variance, Gottschalk's non-local logarithmic term
  and random-hinge placement as the baseline; single-scattering tail as a
  configurable second stage), with configurable step-size control (fractional
  energy loss per step, range floor, voxel boundary);
- energy-deposition and dose scoring on a rectilinear grid with atomics on
  accelerated backends;
- batch-based mean/variance estimation using the fixed history-to-stream
  mapping (float32 batch tallies promoted to float64 between launches), with
  Warp's deterministic execution mode exposed as an option for reproducible
  cross-backend comparisons;
- pencil-beam source with configurable energy, spot size, divergence and
  direction;
- reference Python execution of the same physics as a per-history oracle,
  plus an independent analytical/numerical check;
- the statistical cross-backend consistency criterion for stochastic output
  recorded as a decision before the first comparison.

Depends on: Stage 0.

Validation milestone **V1**: (a) energy conservation per history and per run;
(b) R80 of the depth-dose curve versus the CSDA range (expected to agree to a
few tenths of a millimetre for protons in water once straggling is included,
with the exact tolerance fixed in the stage decision) and Bragg-peak position
below the CSDA range by the expected amount; (c) Bragg-peak width versus
analytical straggling expectations; (d) lateral spread of a pencil beam versus
Fermi–Eyges theory and Gottschalk's 1993 multiple-scattering measurements;
(e) reference-Python vs Warp CPU vs Warp CUDA statistical consistency on
identical configurations, and bitwise run-to-run reproducibility of each
backend in deterministic mode with a fixed seed; (f) random-stream tests:
per-history streams reproduce independently of launch geometry and batch
size, and standard independence tests across history indices pass; (g)
convergence with step-size parameter (R80 drift and peak-height drift between
production and ten-fold finer steps).

### Stage 2 — Nuclear interactions for protons

Serves: *Validated physical behavior* for clinically relevant dose accuracy;
*Conservation checks*; *Configurable physics models*.

Content, staged:

1. primary-fluence attenuation from total inelastic cross sections (ICRU 63
   class data; roughly one percent of primaries per centimetre of water);
2. secondary-particle production and local/transported energy deposition
   following a documented fast-MC approach — the Fippel–Soukup (2004) model
   family used by MCsquare and FRED (transported secondary protons and light
   ions, locally deposited recoils and electrons, escaping neutrons and
   photons) is the leading candidate; alternatives are compared and decided
   at stage start;
3. configurable fidelity levels so the accuracy/performance tradeoff can be
   measured.

Depends on: Stage 1.

Validation milestone **V2**: energy, charge and baryon-number bookkeeping per
interaction (including explicitly accounted escaping energy); peak-to-plateau
ratios and integral depth-dose curves in water versus published measured
pristine Bragg curves (open-access commissioning datasets), publicly
available Geant4/TOPAS reference data, and, if it can be built and run in the
experiment environment, the Apache-2.0 MCsquare engine as an independently
implemented reference; if no adequate public reference exists, a
plug-and-play reference-calculation request per `EXPERIMENT.md` is prepared.
Quantities, metrics and tolerances are fixed in the stage decision.

### Stage 3 — Voxelized heterogeneous geometry and materials

Serves: *Voxelized geometries*, *Material assignment*, *Arbitrary beam
incidence*, *Changeable scoring grids*, *Separation of transport from data
source*, optionally *CT conversion*.

Content:

- voxel geometry with per-voxel material and mass density;
- material library and material-dependent physics quantities (stopping-power
  ratios, scattering power, nuclear cross sections) through the Stage 0 data
  layer;
- robust ray/voxel traversal with arbitrary incidence and beam-frame
  transforms;
- scoring grids decoupled from the transport grid (different resolution and
  alignment).

Depends on: Stages 1 and 2 (transport physics that is being made
material-dependent).

Validation milestone **V3**: water-equivalent-thickness equivalence in slab
phantoms (density scaling), interface behaviour in layered phantoms, rotated
versus axis-aligned beam equivalence, and grid-independence of integral dose
under scoring-grid changes.

### Stage 4 — Treatment-planning scoring and influence matrices

Serves: *Ion-therapy treatment planning orientation*, *Beamlet-resolved
scoring*, *Dose influence matrices*, *Sparse influence representations*,
*Batched beamlet execution*, *Planning-oriented reuse*, *Extensible scoring*,
*LET scoring*, *Track-averaged LET*, *Fluence scoring*, *Species-resolved
scoring*, *Energy-resolved scoring*, *External biological lookup data*,
*Multi-dimensional lookup dimensions*, *Planning-aware uncertainty*.

Content:

- beamlet identity carried through transport; batched multi-beamlet launches
  without reinitialising immutable physics/material data;
- sparse dose-influence matrix assembly with thresholding, exported in a
  documented format usable by inverse planning; the first implementation is
  chunked-dense per-beamlet scoring buffers thresholded into a sparse matrix
  per chunk, with a GPU hash-table keyed by (voxel, beamlet) as the candidate
  successor if memory or throughput measurements demand it;
- general scorer architecture (dose, energy deposition, fluence, dose- and
  track-averaged LET, species- and energy-resolved quantities), with the LET
  estimator definitions (restricted or unrestricted stopping power, secondary
  treatment, energy cut-offs) fixed in the stage decision;
- externally supplied lookup tables indexed by particle species, energy,
  material or tissue class and quantity, accumulated during transport
  (biological/microdosimetric quantities) without embedding one biological
  model in the engine;
- uncertainty estimation for beamlet-resolved quantities.

Depends on: Stage 3 (heterogeneous transport) for clinically meaningful use;
API and estimator design can begin after Stage 1.

Validation milestone **V4** (**achieved 2026-09-11**, develop `98251c2`): sum of
beamlet doses equals the corresponding broad-field dose within statistics (`DEV-017`
exact same-partition + `DEV-020` within-statistics); LET estimators versus published
proton LET-in-water data computed with the *same* estimator convention (or versus a
reference calculation whose convention is documented) and versus analytical limits
(`DEV-018`); statistical uncertainty for the volumetric and beamlet-resolved
scorers (`DEV-019`, `DEV-020`); lookup-table accumulation reproduces offline
post-processing on scored spectra (`DEV-021`); sparse-vs-dense matrix agreement
(`DEV-017`).

### Stage 5 — Multi-ion transport (helium, carbon, oxygen)

Serves: *Multi-ion architecture*, *Ion transport* beyond protons.

Content:

- ion stopping power (effective charge, Z² scaling with documented
  corrections) and range data through the data layer;
- ion-specific straggling and scattering parameters;
- nuclear fragmentation model appropriate for a fast treatment-planning MC,
  with transported charged fragments;
- species-resolved scoring and LET for mixed radiation fields.

Depends on: Stages 2 and 4 (the physics and scoring infrastructure it
generalises).

Validation milestone **V5**: helium and carbon Bragg curves and fragment
tails in water versus published measurements or reference calculations;
species spectra versus literature. Reference datasets, quantities and
tolerances are fixed in the stage decision.

### Stage 6 — Performance engineering and release readiness

Serves: *High-throughput execution*, *Reproducible benchmarking*, *Efficient
batching*, *Performance regression testing*, release requirements of
`EXPERIMENT.md`.

Content:

- reproducible benchmark suite under `benchmarks/` (hardware, version,
  backend, problem size, particle count, precision, configuration,
  methodology);
- profiling-driven kernel and memory-layout optimisation, divergence and
  occupancy work, precision studies (`float32` vs `float64` accumulation);
- performance regression tracking;
- first tagged release from `develop` to `main` once release validation is
  defined and passed.

Depends on: representative workloads from Stages 3–4. Micro-benchmarks and
performance instrumentation are introduced earlier as soon as Stage 1 kernels
exist, but architecture-constraining optimisation waits for representative
workloads. Stage 5 workloads are added to the benchmark suite when available
but are not a prerequisite.

Validation milestone **V6**: documented benchmark results with unchanged
scientific validation outcomes before and after each optimisation.

## Cross-cutting tracks

- **Validation (continuous).** Layered: unit → component → transport
  scenario, under the rules stated above.
- **Documentation (continuous).** Architecture, physics, validation and API
  docs updated in the same task that changes behaviour.
- **Decisions (continuous).** Consequential choices recorded under
  `decisions/` with evidence, alternatives and validation outcome.
- **Code quality (continuous).** Formatter, linter and type checker enforced
  by pre-commit and CI from `DEV-002` onward.
- **External reference calculations.** Requested through the experiment
  protocol only where public references are insufficient (expected first for
  Stage 2 nuclear models and Stage 5 fragmentation).

## Dependency summary

```
Stage 0 ──▶ Stage 1 ──▶ Stage 2 ──▶ Stage 3 ──▶ Stage 4 ──▶ Stage 5
                 │                       │          │
                 └── scorer/API design ──┼──────────┘
                 │                       │
                 └── micro-benchmarks ──▶ Stage 6 ◀─┘

Stage 6 requires Stages 3–4; Stage 5 workloads join the benchmark suite
when available but are not a prerequisite.
```

## First milestone (selected)

**Task `DEV-002` — Stage 0, first half: `ionmc` package foundation, unit
conventions, shared-source execution model, quality tooling, and analytical
proton stopping power and CSDA range on all three execution paths.** The data
layer and the tabulated stopping-power source follow immediately as
`DEV-003`, which closes V0.

Rationale (see decision `0002`): it is the smallest unit of work that is
simultaneously (a) scientifically meaningful — electronic stopping power is the
quantity every later transport step integrates, and the analytical layer can
be validated against an I-value-consistent reference table to a tolerance
fixed in advance; (b) architecturally decisive — it forces the shared-source
execution model, the unit conventions and the quality tooling to exist and be
tested, on kernels with control flow and per-thread state, before transport is
built on them; and (c) bounded — its correctness can be checked
deterministically and reviewed in one pull request.

## Status of stages and tasks

| Stage / task | Status | Evidence |
|---|---|---|
| Stage 0, `DEV-002` (package, conventions, execution model, quality tooling, analytic stopping power) | **completed** 2026-09-10 | decisions `0004`–`0006`; `validation/v0_stopping_power.py` on the host runner, CPU + CUDA, run `RUN-20260910T080455Z-e666b469` at `d23c78a`; DEV-002 local validation records; squash-merged into `develop` |
| Stage 0, `DEV-003` (data layer, tabulated stopping power; closes V0) | **completed** 2026-09-10 | decisions `0007`, `0008`; `validation/v0_tabulated_stopping_power.py` on the host runner (CPU + CUDA); DEV-003 local validation record; **milestone V0 closed** |
| Stage 1, `DEV-004` (transport architecture, deterministic CSDA depth dose) | **completed** 2026-09-10 | decision `0009`; `validation/v1_depth_dose_csda.py` on the host runner (CPU + CUDA); DEV-004 local validation record |
| Stage 1, `DEV-005` (energy-loss straggling, Bragg peak, batch uncertainty) | **completed** 2026-09-10 | decision `0010`; `validation/v1_bragg_straggling.py` on the host runner (CPU + CUDA); DEV-005 local validation record |
| Stage 1, `DEV-006` (multiple Coulomb scattering, 3-D transport; closes V1) | **completed** 2026-09-10 | decision `0011`; `validation/v1_lateral_scattering.py` on the host runner (CPU + CUDA); DEV-006 local validation record; **milestone V1 closed** |
| Stage 2, `DEV-007` (proton nonelastic nuclear attenuation, local deposition; opens V2) | **completed** 2026-09-10 | decision `0012`; `validation/v2_nuclear_attenuation.py` on the host runner (CPU + CUDA); `tests/ionmc/test_nuclear.py`; DEV-007 local validation record |
| Stage 2, `DEV-008` (secondary charged-particle transport; closes V2) | **completed** 2026-09-10 | decision `0013`; `validation/v2_secondary_transport.py` on the host runner (CPU + CUDA); `tests/ionmc/test_secondaries.py`; DEV-008 local validation record; **milestone V2 closed** |
| Stage 3, `DEV-009` (1-D voxelized density heterogeneity, WET transport; opens V3) | **completed** 2026-09-10 | decision `0014`; `validation/v3_density_heterogeneity.py` on the host runner (CPU + CUDA); `tests/ionmc/test_voxel_geometry.py`; DEV-009 local validation record |
| Stage 3, `DEV-010` (per-voxel tissue materials via stopping-power ratios) | **completed** 2026-09-10 | decision `0015`; `validation/v3_material_composition.py` on the host runner (CPU + CUDA); `tests/ionmc/test_material_composition.py`; DEV-010 local validation record |
| Stage 3, `DEV-011` (density-heterogeneous 3-D multiple-scattering transport) | **completed** 2026-09-10 | decision `0016`; `validation/v3_scattering_heterogeneity.py` on the host runner (CPU + CUDA); `tests/ionmc/test_scattering_heterogeneity.py`; DEV-011 local validation record |
| Stage 3, `DEV-012` (non-water materials on the 3-D scattering path) | **completed** 2026-09-10 | decision `0017`; `validation/v3_material_scattering.py` on the host runner (CPU + CUDA); `tests/ionmc/test_material_scattering.py`; DEV-012 local validation record; squash-merged into `develop` at `58309bf` |
| Stage 3, `DEV-013` (arbitrary beam incidence via a beam frame; rotated-vs-axis-aligned equivalence) | **completed** 2026-09-11 | decision `0018`; `validation/v3_arbitrary_incidence.py` on the host runner (CPU + CUDA); `tests/ionmc/test_arbitrary_incidence.py`; DEV-013 local validation record; squash-merged into `develop` at `ebd719f` |
| Stage 3, `DEV-014` (scoring grids decoupled from the transport grid; grid-independence; **closed V3**) | **completed** 2026-09-11 | decision `0019`; `validation/v3_scoring_grid.py` on the host runner (CPU + CUDA); `tests/ionmc/test_scoring_grid_decoupling.py`; DEV-014 local validation record; squash-merged into `develop` at `f9e3c90`; **milestone V3 closed** |
| Stage 3, `DEV-015` (3-D voxel grid with ray/voxel DDA traversal; arbitrary incidence through true voxel geometries) | **completed** 2026-09-11 | decision `0020`; `validation/v3_voxel_grid_3d.py` on the host runner (CPU + CUDA); `tests/ionmc/test_voxel_grid_3d.py`; DEV-015 local validation record; squash-merged into `develop` at `d0df873` |
| Stage 4, `DEV-016` (lab-frame 3-D dose scoring on the voxel-grid path; prerequisite for beamlet influence matrices) | **completed** 2026-09-11 | decision `0021`; `validation/v4_dose3d.py` on the host runner (CPU + CUDA); `tests/ionmc/test_dose3d.py`; DEV-016 local validation record; squash-merged into `develop` at `e8f89df` |
| Stage 4, `DEV-017` (beamlet-resolved scoring and sparse dose-influence matrices; V4 beamlet-sum gate) | **completed** 2026-09-11 | decision `0022`; `validation/v4_influence.py` on the host runner (CPU + CUDA, all 5 gates); `tests/ionmc/test_influence.py`; DEV-017 local validation record; squash-merged into `develop` at `63883c4` |
| Stage 4, `DEV-018` (dose-averaged LET (LET_d) scoring; V4 LET gate) | **completed** 2026-09-11 | decision `0023`; `validation/v4_let.py` on the host runner (CPU + CUDA, all 6 gates); `tests/ionmc/test_let.py`; DEV-018 local validation record; squash-merged into `develop` at `3b5ddbf` |
| Stage 4, `DEV-019` (batch-based statistical uncertainty for 3-D dose; V4 "within statistics") | **completed** 2026-09-11 | decision `0024`; `validation/v4_uncertainty.py` on the host runner (CPU + CUDA, all 5 gates); `tests/ionmc/test_uncertainty.py`; DEV-019 local validation record; squash-merged into `develop` at `6eb245b` |
| Stage 4, `DEV-020` (beamlet-resolved planning-aware uncertainty; V4 beamlet-sum within statistics) | **completed** 2026-09-11 | decision `0025`; `validation/v4_beamlet_uncertainty.py` on the host runner (CPU + CUDA, all 5 gates); `tests/ionmc/test_beamlet_uncertainty.py`; DEV-020 local validation record; squash-merged into `develop` at `83ec95b` |
| Stage 4, `DEV-021` (energy-resolved fluence-spectrum scoring + lookup-table accumulation; closes V4 lookup gate) | **completed** 2026-09-11 | decision `0026`; `validation/v4_fluence.py` on the host runner (CPU + CUDA, all 5 gates); `tests/ionmc/test_fluence.py`; DEV-021 local validation record; squash-merged into `develop` at `98251c2` — **milestone V4 achieved** |
| Stage 5, `DEV-022` (first non-proton ion: helium-4 via effective-charge scaling of the proton stopping power) | **completed** 2026-09-11 | decision `0027`; `validation/v5_helium.py` on the host runner (CPU + CUDA, all 5 gates); `tests/ionmc/test_helium.py`; DEV-022 local validation record; squash-merged into `develop` at `755df41` |
| Stage 5, `DEV-023` (carbon-12 and oxygen-16 primary transport; z² scaling accuracy vs Bethe) | **in progress** 2026-09-11 | decision `0028`; `validation/v5_carbon.py`; `tests/ionmc/test_carbon.py` |
| Stages 4–6 | not started | — |

## Change log

- 2026-09-10: roadmap established (`DEV-001`), revised after independent
  review (V0 criterion scope, I-value consistency, requirement coverage,
  Stage 0 split into `DEV-002`/`DEV-003`, network contingency, RNG and
  determinism validation items, dependency diagram).
- 2026-09-10: `DEV-002` completed; status table added. Decision `0003`
  (kickoff workaround) superseded by direct worktree access. The
  "provisional until Stage 1" status of decision `0005` is unchanged.
- 2026-09-10: `DEV-003` completed, closing **Stage 0** and milestone **V0**:
  external data layer (decision `0007`) and tabulated PCHIP stopping power
  (decision `0008`).
- 2026-09-10: `DEV-004` began **Stage 1**: the proton transport architecture
  and continuous-slowing-down longitudinal depth dose (decision `0009`),
  deterministic (no straggling/scattering yet).
- 2026-09-10: `DEV-005` added Bohr energy-loss straggling and the Bragg peak
  with batch-based uncertainty (decision `0010`), and decoupled the physics
  step from the scoring grid.
- 2026-09-10: `DEV-006` added multiple Coulomb scattering with 3-D transport,
  the random-hinge lateral spread and 2-D depth-lateral scoring (decision
  `0011`), **closing milestone V1** (proton transport in homogeneous water).
  Next: **Stage 2**, nuclear interactions for protons.
- 2026-09-10: `DEV-007` began **Stage 2**: proton nonelastic nuclear
  attenuation with catastrophic primary removal (analytic ICRU-63-shape oxygen
  cross section), local-vs-escaping energy split, and an audited escaping
  channel for later secondary transport (decision `0012`). Reproduces the
  published primary survival to the Bragg peak (0.81 at 150 MeV, 0.72 at
  200 MeV). **Opens milestone V2**; V2 closes once secondary charged-particle
  transport (DEV-008) consumes the escaping channel. Data-layer note: the
  ICRU-63-derived tabulated cross section is *not* vendored (licensing); an
  openly-licensed TENDL-2021 tabulated path is deferred.
- 2026-09-10: `DEV-008` added secondary charged-particle transport (decision
  `0013`), **closing milestone V2** (nuclear interactions for protons in water).
  Each nonelastic reaction is partitioned into local heavy fragments
  (`f_heavy = 0.12`), transported secondary protons (`f_p = 0.50`, sampled with
  an evaporation+cascade energy spectrum and forward emission), and truly
  escaping neutrons/gammas/binding (`f_esc = 0.38`). A two-pass architecture
  (primaries emit one reaction record each; secondaries are generated host-side
  and transported through the same proton engine) keeps the reference and Warp
  paths bit-parity. Secondary protons contribute ~1-2 % of the local dose at
  entrance and a few percent of the total, forming the nuclear plateau; the
  energy budget stays exact. Deferred (with quantitative justification): d/t/
  alpha and neutron/gamma transport, tertiary reactions, and the ICRU-63/TENDL
  tabulated double-differential path. Next: **Stage 3**, voxelized heterogeneous
  geometry and materials.
- 2026-09-10: `DEV-009` began **Stage 3**: 1-D voxelized density heterogeneity
  along the beam axis (decision `0014`). A `VoxelSlab` gives each voxel its own
  mass density; the depth-dose transport (reference and Warp) looks up the local
  density per step and limits the step to the voxel boundary, tracking a
  non-decreasing voxel index (forward transport). Consecutive equal-density
  voxels are merged, so a uniform slab reproduces the homogeneous `WaterSlab`
  bit-for-bit. **Opens milestone V3**: validated water-equivalent-thickness
  scaling (R80 at R_water/rho to ~1e-4), layered-interface peak shift (exactly
  the extra WET), energy conservation with nuclear+secondaries across an
  interface, and cross-backend agreement. Deferred to later Stage-3 tasks: 3-D
  voxel geometry and arbitrary beam incidence, per-voxel *material* composition
  (stopping-power ratios, per-voxel <Z/A>), and scoring grids decoupled from the
  transport grid.
- 2026-09-10: `DEV-010` added per-voxel **tissue materials** (decision `0015`).
  A tissue library (ICRU-44: cortical bone, adipose, soft tissue, muscle, lung,
  air) plus a stopping-power ratio: each voxel transports as water at its
  water-equivalent density `SPR(material) x rho`, computed from the analytic
  Bethe mass-stopping-power ratio, and with a composition-scaled oxygen-
  equivalent nuclear density. The depth-dose transport is unchanged from DEV-009;
  the material physics enters through the per-voxel density and nuclear arrays.
  Validated water-equivalent ratios in the published bands (bone 1.70, adipose
  0.97, muscle 1.04), R80 at R_water/WER, ~2x/3.4x nuclear scaling for bone/
  adipose, energy conservation across material interfaces, and cross-backend
  agreement. Deferred: energy-dependent SPR / per-material stopping tables
  (bone's 2.5% energy dependence; soft tissue < 0.25% with the scalar), element-
  specific nuclear cross sections, and heterogeneous-material scattering. Next:
  the remaining Stage-3 work (3-D voxel geometry with arbitrary incidence,
  decoupled scoring grids) toward milestone V3 closure.
- 2026-09-10: `DEV-011` brought the **3-D multiple-scattering path to 1-D
  voxelized density heterogeneity** (decision `0016`), reusing the DEV-009
  per-voxel-density-by-depth pattern: the scattering reference driver and Warp
  kernel look up the local water density, limit the step's depth advance to the
  voxel boundary, and advance a non-decreasing voxel index. `run_scattering` now
  accepts a density-heterogeneous water `VoxelSlab` (still rejecting non-water
  materials, whose MCS needs the physical density and material X0 — a later
  task). Validated: homogeneous equivalence bit-exact, lateral `sigma_x` vs the
  Fermi-Eyges oracle (uniform density) and a new piecewise-density Fermi-Eyges
  variant (layered interface) to < 1 %, energy conservation, and cross-backend
  agreement. This closes the density-heterogeneity half of the scattering path;
  materials on the scattering path and 3-D geometry/arbitrary incidence remain
  for V3 closure.
- 2026-09-10: `DEV-012` added **non-water materials on the 3-D scattering path**
  (decision `0017`), mirroring DEV-010 for the scattering kernel: the energy loss
  uses the per-voxel water-equivalent density and the multiple scattering the
  per-voxel physical density and material radiation length (looked up by depth).
  The `is_water_only` guard is lifted; the voxel merge now collapses on the full
  per-voxel physics tuple. Validated a homogeneous bone slab and a water/bone/
  water interface against a new material-aware Fermi-Eyges oracle (energy vs
  depth from the water-equivalent thickness, scattering power from the physical
  density and material X0) to < 1 %, energy conservation, and cross-backend
  agreement. The depth-dose and scattering paths now have the **same** per-voxel
  material capability; the remaining Stage-3 work for V3 closure is 3-D voxel
  geometry with arbitrary beam incidence and decoupled scoring grids.
- 2026-09-11: `DEV-013` added **arbitrary beam incidence via a canonical beam
  frame** (decision `0018`). The pencil beam enters at an arbitrary position and
  unit direction and the whole scene (beam + slab normal) can be rigidly rotated;
  the scattering path transports in a beam frame (origin at the entry point, +z'
  along the beam) and looks up the voxel by the material coordinate `u = normal .
  position`, so the (non-equivariant) scattering sampler is unchanged and MCS is
  automatically rotation-covariant. When `normal = direction = +z` it reduces to
  the DEV-012 path bit-for-bit. A `scattering` toggle adds a deterministic
  straight-ray mode. Validated: a rigidly rotated scene reproduces the
  axis-aligned beam-frame depth dose (bit-exact with scattering off, within
  statistics with it on), an oblique beam traverses the same water-equivalent
  path as a normal beam through `D/cos(theta)`, and reference/Warp CPU/CUDA agree.
  Remaining for **V3 closure** (DEV-014): scoring grids decoupled from the
  transport grid in resolution and alignment (the grid-independence gate).
- 2026-09-11: `DEV-014` decoupled the **scoring grid** from the transport grid in
  resolution *and* alignment (decision `0019`), **closing milestone V3**. The
  scoring grids gain an `origin_mm` (`DepthDoseGrid`) and a `depth_origin_mm` /
  `lateral_center_mm` (`DepthLateralGrid`); the reference and Warp deposition
  index bins by `floor((z - z_origin)/dz)`, defaults reproducing the prior grids
  bit-for-bit. Validated: the total deposited energy and the energy per depth
  interval are grid-independent under resolution changes (a 4x-finer grid summed
  4:1 equals the coarse, to round-off) and depth-origin/lateral-centre shifts; a
  grid past the entrance captures strictly less; lateral `sigma_x` is invariant
  under a lateral shift; and reference/Warp CPU/CUDA agree on a shifted grid.
  **Milestone V3 (voxelized heterogeneous geometry and materials) is closed.**
  Deferred (Stage-3 capability beyond the V3 gates): a full 3-D voxel grid with
  arbitrary per-voxel material maps and robust ray/voxel (Siddon/DDA) traversal
  for true patient geometries, a follow-on task and a Stage-4 prerequisite.
- 2026-09-11: `DEV-015` added a **3-D voxel grid with ray/voxel DDA traversal**
  (decision `0020`), the remaining Stage-3 *Voxelized geometries* capability
  (beyond the V3 gates, which used slab/layered phantoms). A `VoxelGrid3D`
  (lab-axis-aligned, per-voxel mass density, single material — density-only cut)
  is traversed by an Amanatides-Woo DDA phrased in the beam-frame idiom: the 1-D
  material coordinate `u = n̂·p` and its limit generalise to three lab-axis
  coordinates `u_k = p0[k] + mᵏ·q` (mᵏ = rows of the beam-frame R) with a per-axis
  distance-to-next-face; geometry state is recomputed each step, so no float state
  drifts between backends and the 1-D `VoxelSlab` is an exact sub-case (Nx=Ny=1).
  Parallel drivers keep the 1-D path untouched (bit-exact). Validated: a
  single-column grid reproduces the VoxelSlab bit-for-bit (uniform-spacing case);
  a homogeneous box
  reproduces WaterSlab at the step-partition level; an oblique beam through an
  off-axis dense insert stops where an independent Siddon WET integral reaches the
  water CSDA range; energy is conserved for a contained beam; and reference/Warp
  CPU/CUDA agree. Deferred: per-voxel material map, grid rotation, a 3-D lab-frame
  scoring volume, non-uniform spacing, CT ingestion.
- 2026-09-11: `DEV-016` began **Stage 4** with a **lab-frame 3-D dose scoring
  grid** (decision `0021`), the volumetric-dose prerequisite for beamlet-resolved
  influence matrices. A `DoseGrid3D` (lab-axis-aligned, own resolution/alignment)
  is scored on the 3-D voxel-grid transport path: each step deposits its energy at
  its lab midpoint (`x = p0 + R·q`) into the containing dose voxel, conserving
  energy exactly for a beam contained in the dose grid. Validated (reference +
  Warp CPU/CUDA): dose sum equals the deposited energy; the 3-D dose z-marginal
  reproduces the beam-frame depth-dose R80 within a dose voxel; the total dose is
  grid-independent; and the backends agree. Deferred: a 3-D scorer on the
  1-D-geometry path, path-length-splitting deposition, and density-driven dose
  units.
- 2026-09-11: `DEV-017` added **beamlet-resolved scoring and sparse dose-influence
  matrices** (decision `0022`). `TransportEngine.run_scattering_multi` transports a
  list of beamlets together in one batched launch (each seeded `seed + i`, keeping
  its `beamlet` id and per-history streams, without reinitialising physics data);
  `assemble_influence_matrix` transports each beamlet with a `DoseGrid3D` scorer
  and thresholds its dose into a CSR `SparseInfluenceMatrix` (one row per beamlet,
  exported to `.npz`). Because the per-beamlet and batched runs share the identical
  seed-partitioned histories, the sum of the beamlet rows equals the batched
  broad-field dose to round-off -- the first **V4** gate (sum of beamlet doses ==
  broad-field dose), validated at ~1e-18 (reference). A 1 % peak threshold keeps
  > 99 % of the energy while making the matrix sparse. Deferred: a GPU
  (voxel, beamlet) hash-table assembly, per-beamlet scoring inside the kernel
  (currently one launch per beamlet), LET/fluence/species-resolved scorers, and
  beamlet-resolved uncertainty. A review follow-up corrected the cross-backend
  validation story: the tight, discriminating per-voxel spatial parity is
  *deterministic* (scattering off, reference-vs-CPU and CPU-vs-CUDA to the float32
  budget), while under scattering on the per-voxel dose decorrelates for *any*
  backend pair (float32 CPU vs CUDA is not bit-identical) and only the total stays
  tight -- a genuine statistical spatial comparison is a deferred cross-cutting
  item. Squash-merged into `develop` at `63883c4`.
- 2026-09-11: `DEV-018` begins the **V4 LET gate** with **dose-averaged LET
  (LET_d) scoring** (decision `0023`), a MUST requirement and the next named V4
  validation item (LET estimators vs published proton LET-in-water data and
  analytical limits). LET_d is scored on the lab-frame 3-D grid alongside dose as
  `Σ ε·L / Σ ε`; the denominator is the existing per-voxel dose energy and the
  numerator an added `Σ ε·L` accumulator on the reference and Warp CPU/CUDA paths.
  Validated (decision 0023): thin-voxel LET_d reproduces NIST PSTAR water within
  ~0.3 %, entrance ≈0.545 keV/µm rising to >10 keV/µm distally, dose-weight
  consistent, voxel-size convergent (~0.3 % over 0.5/1/2 mm), and deterministic
  per-voxel LET_d agrees across backends. Squash-merged into `develop` at
  `3b5ddbf`.
- 2026-09-11: `DEV-019` adds **batch-based statistical uncertainty for 3-D dose**
  (decision `0024`), the MUST statistical-uncertainty requirement for the
  volumetric scorer and the V4 "within statistics" anchor.
  `run_scattering_batched` runs independent history batches (seed `seed+1+b`) into
  a `DoseGrid3D` and reports the per-voxel mean and standard error of the mean
  (`std(batches, ddof=1)/√n_batches`), plus a high-dose-region relative-uncertainty
  summary and, with `score_let`, the LET_d mean and SEM. Validated: the SEM obeys
  the `1/√N` law (4× histories → ~half the SEM over a shared high-dose mask), the
  batch mean conserves energy, and the batched-mean total agrees across backends.
  Squash-merged into `develop` at `6eb245b`.
- 2026-09-11: `DEV-020` adds **beamlet-resolved (planning-aware) uncertainty**
  (decision `0025`), the SHOULD planning-aware-uncertainty requirement.
  `assemble_influence_matrix_batched` transports each beamlet over independent
  batches (seeded `seed + i·n_batches`) and stores the per-voxel mean dose plus its
  standard error (`data_sigma`) in the sparse matrix; `SparseInfluenceMatrix` gains
  `beamlet_sigma_flat`, `total_sigma` (√Σσ² over independent beamlets),
  `beamlet_relative_uncertainty`, and `data_sigma` npz round-tripping. Validated:
  per-beamlet SEM obeys `1/√N`, and the summed per-beamlet means agree with a
  high-statistics broad-field dose **within statistics** (≈99.9 % of high-dose
  voxels within 4σ) — the V4 beamlet-sum gate in its full statistical form,
  complementing DEV-017's exact same-partition round-off check. Deferred:
  beamlet-resolved LET_d/fluence uncertainty and correlated-uncertainty
  propagation into an optimiser. Squash-merged into `develop` at `83ec95b`.
- 2026-09-11: `DEV-021` closes the last **V4** milestone gate with
  **energy-resolved fluence-spectrum scoring + lookup-table accumulation**
  (decision `0026`). A `FluenceSpectrum` histograms the track length `Σ w·s` by the
  step-mean energy `E_mid` (the track-length fluence estimator), and an on-the-fly
  accumulator `A_gate = Σ (w·s)·w_tab[bin]` reads the same precomputed bin-centre
  lookup `w_tab` and bin rule as the offline post-processing `Σ_k counts·w_tab`, so
  the two are identical to round-off (`~4e-15` reference) — the milestone's
  "lookup-table accumulation reproduces offline post-processing on scored spectra".
  With the physical lookup `w = S_lin` the accumulator recovers the total step
  energy deposited (`A/dose ≈ 0.998`; the deficit is the terminal energy-cut
  residual + binning). Scored on the reference and Warp CPU/CUDA paths (float64
  accumulators), deterministic per-bin cross-backend to the float32 budget.
  Deferred: per-voxel/angular/secondary-species spectra, restricted lookups,
  exact 1/S sub-binning, log bins, and multiple scoring regions. Squash-merged
  into `develop` at `98251c2`, **completing milestone V4** (all four V4 gates —
  beamlet-sum, LET, uncertainty, and lookup/sparse — are validated on the host
  runner).
- 2026-09-11: **Stage 5 begins.** `DEV-022` adds the first non-proton ion,
  **helium-4**, transported by reusing the proton stopping-power table via
  equal-velocity matching (`E_p = E·m_p/m_He`, the exact mass ratio, not `E/4`)
  and effective-charge `z²` scaling (decision `0027`), validated against the
  repository's *independent* analytic Bethe model (`AnalyticStoppingPower(WATER,
  ALPHA)`, ~0.1 % agreement) and the `R_He(E) = (m_He/z²m_p)·R_p(E·m_p/m_He)`
  Bragg-range identity (600 MeV He → ~15.85 g/cm² ≈ a 150 MeV proton). NIST ASTAR
  is not reachable from the sandbox; adding an ASTAR dataset as a further reference
  is a follow-up. Deferred to later Stage-5 tasks: nuclear fragmentation and the
  fragment tail, carbon/oxygen, and species-resolved scoring. Squash-merged into
  `develop` at `755df41`.
- 2026-09-11: `DEV-023` generalises primary-ion transport to **carbon-12** (z=6)
  and **oxygen-16** (z=8) — no new transport code, since `scale_ion_stopping_table`
  is species-generic — and **quantifies the z² scaling accuracy** against the
  independent Bethe model (decision `0028`): worst-case (at the 10 MeV/u
  end-of-range) helium ~0.3 %, carbon ~1.1 %, oxygen ~2.0 %, with < 0.5 % across
  the 50–400 MeV/u plateau (the error grows with the omitted Barkas ∝ z³ and Bloch
  ∝ z⁴ terms). A 290 MeV/u carbon beam ranges to ~16.3 cm (its clinical range). All
  three clinical ions (He, C, O) are now transportable, closing the *Multi-ion
  architecture* breadth. **Critical caveat:** carbon/oxygen fragment strongly, and
  the distal **fragment tail is not modelled** — this is the primary-particle Bragg
  curve only; nuclear fragmentation is the next Stage-5 task and the remaining V5
  gate.
