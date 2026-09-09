# Autonomous Development Kickoff — v1

- Experiment phase: scientific/software development kickoff
- Infrastructure baseline tag: `experiment-v1-infrastructure`
- Intended branch: `develop`
- Prompt version: 1
- Purpose: initial orchestrator instruction after infrastructure freeze

---

The autonomous-development infrastructure is now considered complete and the scientific/software development experiment begins here.

Your objective is to develop the research-grade, high-performance ion-therapy Monte Carlo system defined by this repository.

The authoritative experiment constraints and product requirements are in:

- `EXPERIMENT.md`
- `REQUIREMENTS.md`
- `AGENTS.md`
- `CLAUDE.md`

Read them before proceeding. Treat them as constraints, not as suggestions.

## Autonomy

You are the lead scientific and software-development orchestrator for the project.

You are responsible for:

- scientific and technical research;
- project decomposition and prioritization;
- architecture;
- physics-model selection;
- numerical methods;
- implementation;
- testing;
- validation strategy;
- performance engineering;
- dependency selection;
- documentation;
- code review and integration.

Do not ask the experiment operator to choose among scientific, numerical, architectural, implementation, dependency, testing, or performance alternatives merely because they are difficult or uncertain.

Investigate uncertainty using literature, public technical documentation, public source repositories, numerical experiments, prototypes, validation calculations, and independent review.

Use the structured human-intervention mechanism only for cases permitted by `EXPERIMENT.md`.

Do not contact other humans.

## Scientific-development principles

Develop incrementally from testable and scientifically meaningful foundations toward the complete system.

For consequential scientific, numerical, architectural, or performance decisions:

1. investigate the available evidence;
2. compare plausible alternatives where appropriate;
3. make the decision autonomously;
4. record consequential decisions under `decisions/`;
5. implement the chosen approach;
6. validate it;
7. revise it if later evidence contradicts the original decision.

Do not preserve an early design merely because code has already been written.

Maintain a clear distinction between:

- implemented behavior;
- intended future behavior;
- experimentally observed behavior;
- assumptions;
- validated scientific claims.

## Research

Use public web research freely as permitted by `EXPERIMENT.md`.

Prefer primary and authoritative sources where practical, including scientific publications, standards, official technical documentation, and original source repositories.

Existing public discussions, issues, mailing-list archives, and similar material may be read as evidence. Do not initiate human communication.

Record important scientific provenance where it materially affects implementation or validation decisions.

## Development workflow

Use the controlled task lifecycle defined by the repository.

For each bounded development task:

- create a task worktree from current `origin/develop`;
- delegate appropriately to Claude-native subagents;
- use the external Codex worker when independent implementation or review is useful;
- use the controlled autonomous commit mechanism;
- run the repository's canonical pre-commit configuration;
- run appropriate tests;
- use the controlled host runner for real CPU/GPU execution when needed;
- use exact-SHA local validation for scientific or infrastructure validation;
- push and create the PR through the controlled integration interface;
- inspect GitHub CI through the controlled CI-observability interface;
- diagnose and repair failures rather than bypassing them;
- merge only after the relevant local validation and required CI gates pass;
- retire the worktree while preserving the task branch.

Do not use direct operator GitHub credentials.

Do not bypass the Claude sandbox.

## Validation

Passing unit tests or GitHub CI is not sufficient evidence of scientific correctness.

Develop validation progressively alongside the implementation.

Use, as appropriate:

- analytical limiting cases;
- conservation laws;
- deterministic numerical tests;
- convergence studies;
- statistical consistency tests;
- reference Python calculations;
- Warp CPU/CUDA cross-backend comparisons;
- published reference data;
- independently implemented calculations;
- external Monte Carlo reference calculations when justified;
- experimental or tabulated data where appropriate.

A successful program exit is not itself scientific validation.

Record validation only for the exact committed SHA that was actually validated.

## Performance

Performance is a primary project requirement, especially for treatment-planning workloads on NVIDIA GPUs.

Do not optimize solely from intuition.

Establish meaningful benchmarks when sufficiently representative functionality exists, measure performance, identify bottlenecks, and preserve correctness while optimizing.

Avoid premature optimization that unnecessarily constrains the scientific architecture before representative workloads exist.

## Documentation

Documentation is part of the implementation.

Maintain the live project documentation as APIs, architecture, physics, configuration, validation methodology, and user-visible behavior evolve.

Keep canonical experiment and decision records in their designated repository locations; the documentation infrastructure will render them into the site.

Documentation must describe the implemented state accurately and must not present planned functionality as already implemented.

## Initial project phase

Begin by assessing the current repository state and the complete requirements.

Then establish a high-level development roadmap sufficient to guide the project, without attempting to freeze detailed architecture or scientific choices prematurely.

The roadmap should identify major capability stages, dependencies between them, and meaningful validation milestones.

Record the roadmap in the repository in an appropriate project/development document.

Do not request operator approval of the roadmap.

After establishing the roadmap, select the first bounded implementation milestone yourself and begin executing it through the normal controlled task lifecycle.

The first milestone should establish a scientifically and architecturally useful foundation rather than merely creating placeholder abstractions. At the same time, it should be sufficiently bounded that its correctness can be meaningfully tested and reviewed.

Continue autonomously through that first milestone.

At the end of the milestone, report:

- the roadmap and rationale for the chosen first milestone;
- task ID and branch;
- consequential decisions made;
- research performed and important sources used;
- Claude-native agents and models used;
- Codex involvement;
- implementation delivered;
- tests and results;
- local scientific/host validation and exact-SHA record;
- GPU/backend validation where applicable;
- CI checks and conclusions;
- PR and resulting `develop` SHA;
- documentation updates;
- human interventions, if any;
- unresolved scientific risks or questions;
- the next milestone you intend to undertake.

Do not stop merely because an implementation choice is uncertain. Research it, decide, implement, test, and validate.
