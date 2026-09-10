# 0002 — Staged development roadmap and first implementation milestone

- Status: accepted
- Date: 2026-09-10
- Task: DEV-001
- Affects: architecture, validation strategy, maintainability, project sequencing

## Problem

At kickoff the repository contains complete autonomous-development
infrastructure but no Monte Carlo source code. `REQUIREMENTS.md` defines a
broad MUST set (proton and ion transport, three execution paths, voxel
geometries, treatment-planning scoring, influence matrices, LET, biological
lookups, external versioned data, reproducibility, uncertainty, validation,
performance). The order in which these capabilities are built determines how
early scientific validation becomes possible, how much early code is thrown
away when later requirements arrive, and whether the architecture is shaped by
treatment-planning workloads (as the experiment demands) or by whatever was
convenient first.

Two sequencing questions must be settled: (1) the staging of capabilities and
their validation milestones; (2) the first bounded implementation milestone.

## Context

- The experiment requires one physical model executed through reference
  Python, Warp CPU and Warp CUDA, with the reference path serving as a
  correctness oracle. Whether one Python source can literally serve all three
  paths is the single most architecture-defining unknown and is cheapest to
  settle before transport code exists.
- Every later transport step integrates electronic stopping power. It is the
  physics quantity with the most authoritative public reference data
  (ICRU 49/90, NIST PSTAR) and the tightest achievable validation tolerance.
- The external-data requirements (versioning, checksums, caching, offline
  reuse, provenance) apply to the very first tabulated quantity, so the data
  layer must exist before or with it.
- Treatment-planning scoring (beamlets, influence matrices) is the primary use
  case, but it cannot be validated without transport in heterogeneous media,
  so it is a later stage even though its API design should start early.
- The host execution environment provides Warp 1.17.0 and numpy only (no
  scipy, no pytest) with no network. Validation scripts that must run there
  are therefore constrained to numpy/Warp.

## Evidence

Kickoff research (two Claude-native research agents, 2026-09-10; web sources
read, no human contacted) established the following facts that the staging
relies on. They are recorded here so later stage decisions can cite or
overturn them.

- Warp CPU kernels execute serially ("Currently, kernels launched on CPU
  devices will be executed in serial", Warp 1.17 documentation,
  https://nvidia.github.io/warp/stable/user_guide/runtime.html). The CPU path
  is therefore a correctness path and the roadmap places throughput work on
  CUDA.
- `Function.__call__` in Warp's `warp/_src/context.py` falls back to the
  original Python function when a `@wp.func` is called from Python scope
  (https://github.com/NVIDIA/warp), which makes a shared-source reference
  path plausible; kernel bodies and the RNG builtins need confirmation by
  experiment (Stage 0).
- Warp's RNG is a counter-seeded PCG-style hash (`warp/native/rand.h`,
  after Jarzynski & Olano, "Hash Functions for GPU Rendering"); Warp 1.15
  introduced a deterministic execution mode for atomic accumulation
  (https://nvidia.github.io/warp/stable/user_guide/execution_and_performance/deterministic_execution.html).
- Published GPU proton/ion engines converge on one history per thread with a
  fixed-size in-kernel secondary stack, float32 with batch-wise
  accumulation, and beamlet-tagged sparse scoring: gPMC (Jia et al., PMB
  57:7783, 2012; Qin et al., 2016), the Mayo engine (Wan Chan Tseung et al.,
  arXiv:1409.8336), MOQUI (Lee et al., PMB 2022; MIT-licensed source at
  https://github.com/mghro/moquimc), FRED (Schiavi et al., PMB 62:7482,
  2017), MCsquare (Souris et al., Med Phys 43:1700, 2016; Apache-2.0),
  goCMC (Qin et al., PMB 2017). MCsquare reports per-beamlet sparse
  influence matrices "without significant increase of the computation time".
- NIST PSTAR/ASTAR are Standard Reference Data (SRD 124) and copyrighted
  (https://www.nist.gov/srd/public-law); MCsquare redistributes PSTAR-derived
  per-material stopping-power tables under Apache-2.0 at immutable Git
  revisions (https://gitlab.com/openmcsquare/MCsquare). ICRU 90 gives
  I(water) = 78 ± 2 eV versus the 75 eV underlying ICRU 49/PSTAR; the
  difference is about 1% in stopping power at therapeutic energies.
- Clinical fast engines using Gaussian (Bohr-variance) straggling, Gaussian
  multiple scattering with random hinge, ICRU 63 / Fippel–Soukup (Med Phys
  31:2263, 2004) nuclear models and no neutron/photon transport reach >96%
  gamma(3%/3 mm) agreement against a commercial MC in a 70-patient study
  (Gadoue & Sahoo, JACMP 26(10):e70266, 2025).
- Neither the orchestrator sandbox nor the host runner had network access
  during the kickoff session; only the Codex worker and GitHub CI did. This
  motivates the acquire/use split of the data layer.

## Candidate approaches

### Staging

1. **Physics-first, homogeneous water, protons; then geometry; then planning
   scoring; then ions; then performance (selected).** Validation is possible
   at every stage against analytical results or public reference data.
2. **Geometry/scoring-framework first, physics later.** Produces a lot of
   architecture (voxel traversal, scorers, influence matrices) with nothing
   physical to validate it against; risks placeholder abstractions, which the
   kickoff prompt explicitly warns against.
3. **Breadth-first minimal end-to-end system (crude physics in voxel geometry
   with beamlet scoring) then deepen every component.** Attractive for early
   demonstrations, but a crude physics core would need to be replaced under a
   large dependent code base, and early validation would be weak everywhere.

### First milestone

- A. **Package foundation + physics-data layer + proton stopping power/range
   on all three execution paths (selected).**
- B. **Full 1D proton transport with Bragg peak in water.** Scientifically
   more visible, but it bundles the shared-source execution model, the data
   layer, RNG, condensed-history stepping, straggling, scattering and scoring
   into one review unit; a failure of the execution-model assumption would
   surface only after all of that was written.
- C. **Backend/RNG/scheduling skeleton with a toy kernel and no physics.**
   Bounded, but a placeholder: nothing scientifically validatable.

## Selected approach

Staging approach 1, recorded in `docs/development/roadmap.md` as Stages 0–6
with validation milestones V0–V6, and first milestone A, split after
independent review into two bounded tasks: `DEV-002` (package, unit
conventions, quality tooling, shared-source execution model stressed on
kernels with control flow and per-thread state, analytical stopping power and
CSDA range) and `DEV-003` (acquire/use data layer and tabulated stopping
power, closing V0). The execution-model decision taken in `DEV-002` is
explicitly provisional until Stage 1 re-examines it under real transport.

Two rules are attached to every validation milestone: acceptance quantities,
metrics, tolerances and reference datasets are fixed in the stage-opening
decision before any comparison is run; and decision `0001` is understood to
cover only deterministic Warp CPU-vs-CUDA parity, so the reference-Python
versus Warp criterion (Stage 0) and the stochastic cross-backend criterion
(Stage 1) are separate decisions.

## Rationale

- Approach 1 keeps every stage independently validatable: V0 against
  reference stopping-power tables; V1 against range/straggling/scattering
  theory and cross-backend statistics; V2 against measured Bragg curves; V3
  against water-equivalence and rotation invariance; V4 against
  broad-field/beamlet-sum consistency and published LET; V5 against ion Bragg
  curves; V6 against unchanged validation under optimisation.
- Milestone A retires the largest architectural risk (shared-source
  execution across reference Python and Warp) on a small, deterministic
  kernel where decision `0001`'s parity criterion applies directly, before
  any stochastic transport depends on it.
- Milestone A delivers real, validated physics (stopping power, CSDA range)
  rather than placeholders, and the data layer it needs is a MUST in its own
  right.
- Treatment-planning orientation is preserved by making scorer/beamlet API
  design a cross-cutting track that begins after Stage 1 and by deferring
  architecture-constraining optimisation until representative beamlet
  workloads exist (Stage 6 depends on Stages 3–4).

## Expected tradeoffs

- No dose distribution is produced until Stage 1 completes; the first
  milestone is less visible than a Bragg peak.
- Approach 1 delays heterogeneous geometry, so early transport code must be
  written with material-dependence in mind (quantities obtained through the
  data layer, never hard-coded water constants inside transport logic) to
  avoid a rewrite at Stage 3.
- Sequencing ions after the proton planning stack means ion-specific
  architecture needs (fragment species, charge states) must be anticipated in
  the particle state layout from Stage 1 to avoid a proton-only assumption,
  which `REQUIREMENTS.md` forbids.
- The analytical stopping-power layer must be validated against an
  I-value-consistent reference (75 eV against PSTAR-derived tables); the
  ICRU 90 value is a separately quantified model difference. Comparing an
  ICRU 90 model to PSTAR tables would otherwise bake a known one-percent
  offset into the tolerance.
- Automatic data acquisition cannot be exercised from the sandbox or the host
  runner while they lack network access; it is exercised in GitHub CI and via
  the Codex worker, with a documented small committed water table as the
  fallback if no available environment can fetch the dataset.
- The shared-source execution model is first tested on deterministic
  kernels. Even with control flow and per-thread state included, this is
  weaker than real transport, which is why the decision is provisional.

## Validation strategy

The roadmap itself is validated by use: each stage's validation milestone must
be satisfied for a committed SHA before the stage is declared implemented in
`docs/development/roadmap.md`. Deviations from the staging are recorded as
amendments to this decision or superseding decisions.

## Later validation outcome

Not yet available (no stage has completed).
