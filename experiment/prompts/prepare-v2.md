# Codex Task — Prepare Autonomous IonMC Experiment v2

You are working directly in the existing repository:

`https://github.com/e0404/autonomous-ion-mc`

This task is **experiment/infrastructure setup for a second autonomous-development run**, not scientific implementation of the Monte Carlo itself.

The first autonomous run (v1) has completed and has been independently audited. Preserve it as an experimental result. Do not rewrite or clean up its history in a way that obscures what happened.

## Goal

Prepare a clearly versioned **v2 autonomous-development experiment condition** that keeps the strong autonomy of v1 while addressing weaknesses exposed by the independent audit.

The key design principle is:

> The autonomous orchestrator should remain free to decide which physics models, validation strategies, reference calculations, literature sources, and external datasets are scientifically appropriate.  
> The infrastructure should give it substantially better access to independent evidence and make unsupported claims, silent simplifications, weak outputs, and vague performance goals harder.

Do **not** turn v2 into a scripted benchmark against a fully pre-defined validation dataset.

Do **not** reduce the scientific scope to “protons first, ions optional.” Heavy-ion transport, nuclear interactions, and fragmentation are central to the difficulty and scientific value of the experiment.

## 1. Preserve v1 provenance

Before modifying the experiment definition:

- inspect the repository history, tags, experiment prompts, decisions, and current branches;
- identify the exact v1 start and final/audited states;
- preserve them with clear tags or immutable references if not already present;
- do not silently alter the canonical v1 kickoff prompt, experiment protocol, or requirements after the fact.

Introduce v2 as a new experiment condition with separately versioned files or clearly versioned sections.

Prefer explicit provenance such as `experiment/v1/...` and `experiment/v2/...`, or another clean design consistent with the existing repository.

Do not duplicate large files unnecessarily if a versioned structure or manifest can reference them cleanly.

## 2. Use the independent audit as input

The v1 audit found, among other things:

- release validation could report `release_ready=true` with zero suites executed;
- 3D planning paths silently ignored requested nuclear/secondary physics;
- the thin-step stochastic straggling algorithm became strongly biased rather than converging;
- carbon fragmentation had source-position, seed, and history-scaling defects;
- heavy-ion validation was too dependent on internally related models and proton-derived data;
- complete biological lookup support and self-describing result provenance were missing;
- fine or shifted scoring grids could alias badly;
- escaping-particle energy accounting was incomplete;
- user-facing end-to-end workflows and human-readable outputs were weak;
- performance evidence did not represent a complete realistic planning workload;
- the autonomous agent never invoked the structured intervention mechanism.

Read the audit material if available in the repository/workspace. Treat it as evidence about v1, not as a fixed specification for every v2 scientific choice.

Do not simply patch each v1 bug. Improve the **experimental condition** so a fresh autonomous run has better incentives, tools, observability, and acceptance boundaries.

## 3. Preserve scientific autonomy in validation

Do **not** prescribe a fixed complete validation matrix to the orchestrator.

The v2 orchestrator should still be expected to:

- research literature autonomously;
- locate independent public datasets;
- determine appropriate validation observables;
- choose reference codes and physics settings;
- identify potential circularity in validation evidence;
- decide when TOPAS, MCsquare, FRED, analytic theory, tabulated data, or experimental measurements are suitable;
- document why a chosen reference is scientifically appropriate;
- revise its validation strategy when evidence contradicts its assumptions.

The experiment should test autonomous scientific research capability, not only autonomous code generation.

However, strengthen the protocol so that backend agreement, self-consistency, conservation identities, or comparison to a model derived from the same underlying data must not automatically be treated as sufficient independent physical validation of a major physics capability.

Require the agent to **classify the independence and provenance of its evidence** and justify sufficiency.

A useful evidence taxonomy may include concepts such as:

- measured/experimental evidence;
- authoritative ion-specific/tabulated reference data;
- independent Monte Carlo implementation;
- independent analytic/theoretical result;
- related/correlated model;
- backend parity/internal oracle;
- self-consistency/conservation.

Do not hard-code a universal ranking if a more scientifically useful representation emerges, but make validation provenance and independence explicit.

## 4. Improve access to independent references

Provide the v2 orchestrator with controlled autonomous access to reference engines.

At minimum, design/implement controlled runners or adapters for:

- TOPAS / Geant4;
- MCsquare;
- FRED, where feasible in the available environment.

The orchestrator should be able to invoke these without manual copy/paste or operator shell work.

The runners should:

- execute in controlled environments;
- not expose unrelated host credentials or data;
- record exact executable/container/version information;
- record physics configuration;
- record source, geometry, material, cuts, histories, seeds, and scorers;
- retain stdout/stderr;
- preserve raw reference outputs;
- expose bounded structured output back to the orchestrator;
- provide reproducible run IDs;
- make it possible to compare exact autonomous-code SHAs to exact reference configurations.

Prefer a common conceptual interface for reference calculations, but do not overconstrain the capabilities of TOPAS, MCsquare, or FRED merely to force them into the same lowest-common-denominator schema.

The orchestrator should be able to create its own reference input files/configurations.

If direct Geant4 access is practical in addition to TOPAS, consider enabling it, but TOPAS is sufficient as the primary flexible Geant4-based reference runner.

## 5. Improve access to external scientific data

Strengthen the data/research environment rather than pre-supplying all validation answers.

The v2 agent should be able to autonomously:

- search public literature and technical documentation;
- download legally usable reference data;
- retrieve authoritative stopping/range tables such as PSTAR/ASTAR where permitted;
- use published experimental depth-dose, fragmentation, attenuation, scattering, LET, or other datasets;
- cache exact downloaded files with hashes and provenance;
- distinguish implementation data from validation/reference data;
- record whether a dataset was used for model construction/calibration versus independent evaluation.

Do not provide a complete preselected validation corpus unless required for infrastructure testing.

The experiment should allow the agent to demonstrate that it can find appropriate independent scientific evidence itself.

## 6. Strengthen the intervention mechanism

The v1 orchestrator never requested human intervention.

Keep intervention rare, but make the rule more operational:

If the orchestrator determines that a scientifically important capability requires external evidence that it cannot autonomously obtain or execute after reasonable bounded investigation, it should use the structured intervention mechanism rather than silently substituting weaker self-consistency evidence.

Examples may include:

- inaccessible experimental datasets;
- a reference code requiring unavailable institutional data;
- a TOPAS/Geant4 reference case that cannot be executed with available infrastructure;
- external authorization or data that only the operator can supply.

Do not require intervention merely because a decision is difficult.

Also investigate whether ordinary Claude Code permission/approval prompts can trigger the same out-of-band notification channel used for experiment interventions, while remaining semantically distinct and requiring explicit approval.

This was a practical weakness in v1: unattended runs could stall on an approval prompt without operator notification.

## 7. Add explicit capability-contract requirements

A major v1 failure mode was silently accepting a requested physics configuration while executing a simpler model.

For v2, require that public simulation APIs follow a fail-closed capability contract:

- a requested physics feature must actually be active on the selected path;
- unsupported combinations must raise an explicit, actionable error;
- simulation results must record the **effective** physics configuration;
- no user-visible flag may silently become a no-op;
- the supported species / material / geometry / energy / scoring domains must be discoverable and documented.

This should apply especially to:

- nuclear interactions;
- secondary transport;
- fragmentation;
- LET;
- species-resolved scoring;
- biological lookup data;
- influence-matrix generation.

## 8. Strengthen numerical falsification requirements

Do not prescribe the exact numerical algorithms, but require adversarial numerical testing where scientifically relevant.

The orchestrator should test, as appropriate:

- decreasing transport step size;
- scoring-grid refinement;
- scoring-grid shifts relative to transport geometry;
- geometry-boundary alignment/misalignment;
- translated sources;
- sources starting outside a transport volume;
- finite-phantom escape;
- float32/float64 sensitivity;
- seed changes;
- increasing history count;
- batch-count changes;
- backend changes;
- source/beam orientation changes.

A configurable numerical parameter should either exhibit sensible convergence/stability in its documented domain or the unsupported domain should be detected/rejected.

For stochastic estimators, require tests that statistics behave sensibly with increasing histories where theoretically expected.

The purpose is not to fix v1's specific straggling implementation in advance; it is to force the v2 agent to actively attempt to falsify its own numerical methods.

## 9. Keep heavy-ion physics central

Do not redesign the experiment so that successful proton transport constitutes most of the scientific objective.

The v2 requirements should continue to make the difficult ion-therapy features central, including:

- helium/carbon/oxygen-capable architecture;
- ion-specific stopping/range physics;
- nuclear interactions;
- primary attenuation;
- projectile fragmentation;
- transported secondary ions;
- distal fragment dose;
- species-resolved behavior where scientifically relevant;
- LET / biological / microdosimetric extensibility;
- treatment-planning-oriented influence calculations.

The autonomous system may choose development order, but it should not be able to satisfy the intended scientific scope by implementing mature proton-only physics and leaving ion complexity largely untested.

## 10. Add concrete performance boundary conditions

V1 said "high throughput" but did not define sufficiently concrete application-level goals.

For v2, define representative performance/resource expectations before the scientific implementation begins.

Do not force a single arbitrary histories/s value without study.

Instead:

1. define representative treatment-planning workloads;
2. establish the available reference workstation/GPU;
3. define the metrics that must be reported;
4. allow a brief infrastructure/performance-calibration phase to set/ratify realistic targets;
5. freeze those targets before substantial optimization.

Metrics should include, where appropriate:

- wall-clock time;
- histories/s;
- cold-start/compilation time;
- host overhead versus GPU kernel time;
- GPU memory;
- host memory;
- scaling with histories;
- scaling with beamlets;
- scaling with scorer/grid size;
- full enabled physics, not only stripped-down deterministic kernels.

Include at least:

- a realistic 3D proton planning workload;
- a multi-beamlet influence-matrix workload;
- a heavy-ion / fragmentation workload.

The exact target numbers may be decided during the v2 setup/calibration stage, but once ratified they must be versioned and not silently weakened after observing implementation performance.

## 11. Add an end-to-end user/product acceptance contract

V2 must produce a usable scientific artifact, not only an internal Python library and validation scripts.

Require at least one documented end-to-end workflow from a clean installation.

The autonomous system may choose CLI, Python API, config-file workflow, or a combination, but a new scientific user must be able to:

- install the package;
- acquire/cache required data;
- define a source/geometry/scoring configuration;
- run a simulation;
- obtain dose;
- obtain LET where supported;
- obtain uncertainty information;
- obtain influence-matrix output where applicable;
- inspect human-readable results without reading test code or internal implementation details.

Require both machine-readable and human-readable output.

Examples of acceptable artifacts may include:

- HDF5 / NPZ / another self-describing result format;
- CSV profiles;
- JSON run summary;
- plots;
- HTML/Markdown run report.

Do not hard-code these exact formats if the autonomous system develops a better design, but enforce the product-level outcome.

## 12. Require self-describing persisted results

Persisted scientific outputs should contain enough provenance to be interpretable outside the producing Python process.

At minimum consider:

- code version / Git SHA;
- dirty state;
- requested and effective physics;
- particle species;
- source definition;
- geometry/material definitions;
- transport/scoring grids;
- units;
- external physics-data identifiers and hashes;
- reference/model versions;
- history count;
- batching;
- seeds and RNG identity;
- normalization;
- uncertainty definition and validity;
- escaped/truncated energy/particles where relevant;
- execution backend;
- hardware/software environment where appropriate.

Influence matrices in particular must carry enough geometry, normalization, history, and units metadata to be scientifically interpretable after reopening.

## 13. Make release validation fail closed

Redesign the experiment-level release gate so the v1 zero-suite false-positive is structurally impossible.

Requirements:

- required release-validation membership must be explicit or otherwise verifiably complete;
- zero required suites cannot pass;
- missing required validation evidence cannot pass;
- filtered/subset validation runs are diagnostic and must never claim full release readiness;
- exact Git SHA must be recorded;
- required CPU/GPU/backend evidence must be checked;
- required independent-reference evidence must be represented;
- required performance/resource gates must be represented;
- end-to-end workflow/product acceptance must be represented;
- unresolved MUST-level requirements must block release.

The release mechanism should fail closed.

If release criteria pass, the autonomous process should complete the intended promotion/tag/release workflow rather than stopping at a vague "release-ready" state.

## 14. Improve runtime approval and tool-policy robustness

Review v1 orchestration incidents and strengthen the setup where appropriate.

Known observations include:

- direct `gh pr checks` use instead of the controlled CI MCP;
- a hallucinated `/home/wahlm/...` path triggering an outside-workspace permission prompt;
- user-level `blockReadsOutsideWorkingDirectories` causing earlier sandbox friction;
- Claude Code ask/approval prompts not notifying the operator;
- milestone completion initially causing the orchestrator to stop rather than continue.

Preserve these as experiment observations.

Where technically appropriate:

- enforce controlled CI observability rather than relying only on prose;
- avoid granting unnecessary direct GitHub CLI capability if controlled tools cover the workflow;
- make long-running autonomous sessions notify the operator when blocked on a permission prompt;
- retain continuous-execution instructions so task/milestone completion is a checkpoint, not an automatic stop.

Do not overfit infrastructure to a one-off typo such as `/home/wahlm`.

## 15. Telemetry for v2

Extend experiment telemetry so later analysis can answer questions such as:

- how often external reference engines were used;
- which engine/data source was selected and why;
- how often validation strategies changed;
- when internal evidence was rejected as insufficient;
- interventions requested;
- permission prompts encountered;
- scientific failures discovered before integration;
- reference calculations that contradicted the implementation;
- model/tolerance changes after contradictory evidence;
- performance regressions/failures;
- user-workflow acceptance failures.

Where practical, preserve structured records rather than requiring later inference from chat logs.

Do not collect secrets.

## 16. Separate infrastructure/setup from the scientific v2 run

This Codex task should prepare the v2 experiment condition.

Do **not** begin autonomous scientific implementation of IonMC v2 as part of this task.

Before declaring setup complete:

- ensure the experiment protocol and requirements are versioned;
- ensure v1 remains reconstructable;
- ensure controlled reference runners work at least through smoke tests;
- ensure telemetry captures their use;
- ensure intervention/notification mechanisms work;
- ensure repository workflow/CI remains functional;
- ensure the v2 kickoff prompt can reconstruct its instructions from the repo;
- ensure no v1 scientific implementation accidentally becomes the starting implementation for v2 unless explicitly intended and documented.

Prefer a fresh scientific starting point derived from the pre-scientific-development infrastructure state plus v2 infrastructure changes, rather than repairing the completed v1 implementation.

## 17. Deliverables

Implement the v2 experiment setup in the repository with normal repository quality standards.

Produce:

1. a clear v2 experiment/protocol definition;
2. revised v2 scientific/software requirements;
3. revised agent/orchestrator operational instructions;
4. controlled TOPAS/Geant4 reference-runner integration;
5. controlled MCsquare integration;
6. controlled FRED integration if feasible;
7. external-data/research provenance improvements;
8. intervention and operator-notification improvements where technically possible;
9. capability-contract requirements;
10. numerical-falsification requirements;
11. performance-boundary framework;
12. end-to-end output/product acceptance requirements;
13. fail-closed release qualification;
14. telemetry updates;
15. an archived canonical v2 kickoff prompt;
16. documentation explaining the differences between v1 and v2.

Add automated tests/smoke tests for infrastructure where practical.

Do not claim reference runners are scientifically validated merely because smoke tests execute.

## 18. Working style

Inspect the repository before making structural decisions.

Preserve existing useful infrastructure rather than replacing it gratuitously.

Make consequential experiment-design choices explicit in version-controlled decisions or experiment documentation.

Use public documentation/research where needed to configure runner interfaces correctly.

Keep the setup reproducible.

Do not ask the human operator to choose ordinary technical design details. Make reasonable choices, document them, and proceed.

If an actual external authorization or inaccessible institutional dependency blocks setup, report it clearly.

At completion, provide:

- files changed;
- experiment-design decisions;
- v1→v2 differences;
- runner capabilities and smoke-test status;
- unresolved infrastructure limitations;
- exact recommended command/prompt to start the fresh v2 autonomous scientific run.
