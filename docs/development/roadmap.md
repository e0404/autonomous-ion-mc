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

Validation milestone **V4**: sum of beamlet doses equals the corresponding
broad-field dose within statistics; LET estimators versus published proton
LET-in-water data computed with the *same* estimator convention (or versus a
reference calculation whose convention is documented) and versus analytical
limits; lookup-table accumulation reproduces offline post-processing on
scored spectra; sparse-vs-dense matrix agreement.

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
| Stage 2, `DEV-007` (proton nonelastic nuclear attenuation, local deposition; opens V2) | **in progress** 2026-09-10 | decision `0012`; `validation/v2_nuclear_attenuation.py`; `tests/ionmc/test_nuclear.py` |
| Stages 3–6 | not started | — |

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
