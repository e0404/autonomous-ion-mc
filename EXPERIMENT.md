# Autonomous Development Experiment Protocol

## Objective

This repository is an experiment in highly autonomous scientific software development.

The goal is to develop a research-grade, high-performance Monte Carlo simulation system for ion therapy while minimizing expert human intervention. The primary application domain is inverse-planned ion radiotherapy, particularly IMPT; architectural decisions should therefore prioritize efficient repeated transport/scoring workflows relevant to treatment planning rather than general-purpose Monte Carlo flexibility. The goal is a competetive, performance-optimized ion-therapy treatment planning Monte Carlo package.

The autonomous agent system is responsible for scientific research, architecture, implementation, validation, performance optimization, dependency selection, and project management unless this protocol explicitly specifies a constraint or requires human intervention.

The quality of the autonomous decision-making process is itself an experimental outcome.

---

## Implementation Constraints

The following requirements define the intended software and are experiment constraints rather than decisions left to the autonomous system.

### Python-first implementation

The Monte Carlo system must be designed as a Python package with a native Python-facing API.

Core transport physics and interaction logic should remain readable and inspectable at the Python source level. The implementation should favor clear expression of the underlying physical algorithms rather than hiding the principal transport logic in opaque compiled extensions.

### Reference and accelerated execution

The system must support:

1. a native/reference Python execution path suitable for correctness testing, debugging, and inspection;
2. an NVIDIA Warp accelerated CPU execution path;
3. an NVIDIA Warp accelerated CUDA execution path.

The reference implementation should provide an independent or minimally transformed correctness oracle wherever practical.

CPU and CUDA acceleration through Warp should preserve the same physical model rather than becoming independently maintained implementations whose behavior can silently diverge.

### Backend-oriented architecture

The architecture should separate physical modeling from execution sufficiently that additional execution backends could plausibly be introduced in the future.

Implementing or testing additional backends is outside the scope of this experiment unless independently justified by the autonomous system for development or validation purposes.

The detailed software architecture, interfaces, abstractions, memory layout, scheduling strategy, and implementation mechanism remain decisions for the autonomous system.

### Performance and readability

High performance is a primary requirement, particularly on NVIDIA GPUs, but optimization must not make the implemented physics unnecessarily opaque or prevent independent validation.

Performance-oriented duplication or specialization is permitted when justified and validated.

---

## Fundamental Autonomy Rule

Do not ask the human developer to make a scientific, numerical, architectural, implementation, dependency, or performance decision merely because the decision is difficult or uncertain.

When uncertain:

1. investigate the problem;
2. identify plausible alternatives;
3. evaluate available evidence;
4. make a reasoned decision;
5. document consequential decisions;
6. implement the selected approach;
7. validate it;
8. revise the decision if evidence shows that it was inadequate.

Uncertainty is not, by itself, a reason for human intervention.

---

## Permitted Human Interaction

The only permitted human interaction is with the designated experiment operator through the mechanisms defined in this repository.

The autonomous system must not independently contact any other person.

Prohibited activities include, but are not limited to:

* sending email;
* sending direct messages;
* contacting software or data maintainers;
* asking questions on forums;
* posting requests for help on discussion boards;
* opening GitHub issues for the purpose of obtaining assistance;
* commenting on external issues or pull requests to request information;
* posting questions on social media;
* initiating any other person-to-person communication.

Reading publicly available websites, repositories, papers, issue trackers, forum archives, mailing-list archives, and similar sources is permitted.

Existing public human communication may be used as evidence. Creating new human communication is not permitted.

---

## Permitted Human Intervention

Human input from the designated experiment operator may be requested only when one or more of the following conditions applies.

### 1. Inaccessible information

Required information exists but cannot legally or technically be accessed by the agent system.

Examples include:

* a paper PDF behind inaccessible authentication;
* unpublished documentation;
* proprietary reference material;
* data to which only the experiment operator has access.

The request must state exactly what information is required and why.

### 2. External reference calculations

A calculation must be performed using software or infrastructure unavailable to the autonomous system.

Examples include:

* TOPAS;
* specific Geant4 configurations;
* institutional simulation infrastructure.

Reference-calculation requests must satisfy the requirements in the section **Reference Calculation Packages**.

### 3. Physical or externally held data

Required data can only be obtained from physical equipment or external systems controlled by the experiment operator.

### 4. Authorization

An irreversible or externally consequential action requires explicit authorization and is not already pre-authorized by this protocol.

Examples include:

* publication;
* expenditure;
* deployment to infrastructure outside the experimental workspace;
* modification of external projects or repositories not belonging to this experiment.

Operations on the experiment's designated GitHub repository are separately governed by the Git and Remote Repository sections below.

### 5. Ambiguous experimental boundary

Two or more interpretations of an explicit experiment constraint are materially incompatible and cannot be resolved from repository documentation.

This category must not be used for ordinary scientific uncertainty.

---

## Human Intervention Requests

Every human intervention must be recorded.

Where practical, requests should be written to:

```
reference_requests/
```

or another structured experiment log.

Each request must state:

* request identifier;
* requesting agent;
* task that triggered the request;
* information or action required;
* why autonomous resolution was not possible;
* exact requested input;
* decisions that depend on the result;
* relevant Git commit and task identifier.

Human intervention should require as little scientific or technical interpretation by the experiment operator as reasonably possible.

---

## Reference Calculation Packages

Requests for external reference calculations must be provided in a plug-and-play form.

The experiment operator should not be required to translate a scientific request into simulation inputs, decide simulation settings, manually extract requested quantities from raw output, or interpret output formats.

A reference-calculation package should, where technically possible, contain:

* a README containing exact step-by-step execution instructions;
* complete simulation input files;
* required auxiliary data;
* software/version requirements;
* expected runtime and resource requirements where known;
* reproducibility settings including random seeds where relevant;
* an executable launch script;
* an automated output collection or conversion script;
* an automated validation script that checks whether expected output files were produced;
* a manifest describing all returned files;
* units and coordinate conventions;
* requested statistical precision;
* checksums where useful;
* a clear destination for the completed result package.

The preferred interaction should be conceptually equivalent to:

```
cd reference_requests/RR-XXXX
./run.sh
./collect.sh
```

followed by returning the generated result directory.

If the external software cannot be automated completely, the request must minimize unavoidable manual actions and document them exactly.

The experiment operator should not manually parse or scientifically interpret simulation output on behalf of the autonomous system.

The autonomous system is responsible for ingesting, parsing, validating, and interpreting the returned reference data.

---

## Scientific and Technical Decisions

Consequential decisions should be documented under:

```
decisions/
```

A decision is consequential when changing it would materially affect one or more of:

* physical accuracy;
* numerical accuracy;
* reproducibility;
* computational performance;
* architecture;
* maintainability;
* dependency footprint;
* validation strategy;
* scientific interpretation.

Decision records should normally include:

* problem;
* context;
* evidence;
* candidate approaches;
* selected approach;
* rationale;
* expected tradeoffs;
* validation strategy;
* later validation outcome where available.

The autonomous system determines which decisions warrant records.

---

## Validation

Implementation completion is not sufficient.

The autonomous system is responsible for designing and executing appropriate validation.

Validation may include:

* analytic tests;
* unit tests;
* manufactured solutions;
* conservation tests;
* convergence studies;
* statistical consistency tests;
* comparison between reference Python and Warp execution;
* comparison between Warp CPU and CUDA execution;
* comparison with published data;
* comparison with independent implementations;
* comparison with requested TOPAS/Geant4 reference calculations;
* performance regression tests.

Validation methodology should be proportional to the scientific importance of the component.

A performant implementation must not be treated as scientifically correct merely because it executes successfully.

---

## Performance

Performance is a primary requirement, but correctness and quantified accuracy take precedence over unvalidated speed.

The system should explicitly investigate accuracy/performance tradeoffs where relevant.

Performance claims must be reproducible and accompanied by sufficient information about:

* hardware;
* problem size;
* software version;
* configuration;
* numerical precision;
* benchmark methodology;
* statistical workload.

Optimization decisions that materially alter numerical or physical behavior must be validated.

---

## Model and Agent Delegation

The lead orchestrator may delegate work to available subagents or external coding agents.

The orchestrator is responsible for deciding:

* whether delegation is useful;
* which agent/model class is appropriate;
* required reasoning effort;
* whether independent review is needed;
* whether parallel exploration is justified.

Model routing should not be delegated to an opaque external automatic router unless explicitly introduced as part of a later experiment revision.

---

## Code Quality and Validation

The software should be maintained as research-grade production software rather than as a sequence of experimental scripts.

### Code Quality

The repository must use automated pre-commit checks.

Pre-commit checks should cover appropriate mechanical quality requirements such as:

* formatting;
* linting;
* malformed files;
* accidental debug code;
* common repository hygiene problems.

The autonomous system may select the specific tools unless otherwise constrained.

Checks that determine whether code is acceptable for integration must also run in CI. Local Git hooks alone are not considered sufficient because they can be bypassed.

The autonomous system should establish and maintain appropriate additional quality mechanisms, which may include:

* static typing;
* documentation checks;
* unit and integration tests;
* coverage measurement;
* dependency auditing;
* reproducibility checks.

The specific implementation of these mechanisms remains an autonomous engineering decision.

### Test-supported development

Automated testing and validation must accompany implementation work.

For well-specified software behavior, test-first development is preferred where practical.

For exploratory scientific or numerical work, implementation may precede formal tests when necessary to understand the problem, but such work must not be integrated as complete until appropriate automated tests and validation have been added.

Bug fixes should normally include a regression test that fails before the fix and passes afterward.

New scientific functionality should include tests or validation appropriate to the level of physical and numerical risk.

---

## Git Workflow

Git is the authoritative record of project evolution.

### Permanent branches

Two long-lived branches are used:

```
main
develop
```

`develop` is the integration branch for ongoing development.

`main` contains releases only.

Direct implementation work must not occur on either `develop` or `main`.

### Task branches

Development work should occur on isolated task branches, preferably in isolated Git worktrees when parallel agents are active.

Suggested naming includes:

```
task/<task-id>-<description>
fix/<task-id>-<description>
experiment/<task-id>-<description>
```

Task identifiers should be traceable to orchestration or telemetry records.

### Integration into develop

Task branches must be integrated into `develop` through pull requests.

Integration uses squash merges so that each completed task produces one coherent commit on `develop`. Task branches must not be merged into `develop` unless the relevant automated tests and validation pass.

Before integration:

* required tests must pass;
* required quality checks must pass;
* relevant validation must pass;
* consequential decisions must be documented;
* unresolved implementation conflicts must be addressed.

The autonomous system may perform review using independent agents.

Human review is not required unless explicitly requested under the permitted-intervention rules.

### Preservation of task history

Task branches should not be deleted during the active experiment.

The following should be recorded for each integrated task:

* task identifier;
* branch name;
* task branch head SHA;
* pull request identifier;
* resulting squash commit SHA on `develop`;
* agents/models involved where telemetry permits.

This preserves detailed intermediate development history even though `develop` uses squash merges.

### Releases

`main` is updated only from release-ready states of `develop`.

Releases should:

* have passed the defined release validation suite;
* have documented validation status;
* have reproducible benchmark information where relevant;
* be version tagged;
* identify the corresponding experiment configuration.

The autonomous system determines when a release is scientifically and technically justified unless release authorization is explicitly reserved to the experiment operator.

### History integrity

Force pushes to shared task branches, `develop`, and `main` are prohibited except for recovery from an explicitly documented repository error.

Rewriting already-recorded experimental history is prohibited.

---

## Remote GitHub Repository

The experiment should use a designated GitHub repository as a live remote.

Operations on this repository are pre-authorized within the boundaries below and do not constitute human intervention.

The autonomous system may:

* push task branches;
* push the `develop` branch through the defined integration workflow;
* create and update pull requests;
* run and inspect CI;
* create release branches or tags where required by the workflow;
* read repository metadata and previous project history.

The autonomous system must not use the repository to solicit human assistance.

It must not:

* request reviews from external people;
* mention users to obtain help;
* open issues asking for advice;
* use Discussions to ask questions;
* contact maintainers or contributors;
* interact with unrelated repositories except through read-only research unless separately authorized.

During the active experiment, the repository should preferably remain private or otherwise configured to prevent unsolicited human participation.

Remote synchronization should occur sufficiently often that meaningful completed work and task state are preserved outside the workstation.

At minimum, task branches should be pushed at meaningful checkpoints and before integration, and integrated `develop` states should be pushed immediately.

---

## Repository Protection

The remote repository should enforce the intended workflow where technically possible.

`develop` and `main` should be protected against:

* direct pushes;
* force pushes;
* deletion.

Integration should require pull requests and passing required status checks.

A linear history should be maintained on protected integration branches.

Repository configuration that affects the experiment should itself be documented or exported where practical.

---

## Repository Discipline

Changes should be made in logically coherent commits.

Generated artifacts required for scientific reproducibility must have documented provenance.

Do not commit:

* credentials;
* access tokens;
* authentication material;
* private keys;
* unrelated personal data;
* unrelated institutional data.

Large generated or reference datasets should use an appropriate storage mechanism rather than being committed indiscriminately to Git history.

---

## External Research and Software

External web research is permitted when available.

Third-party open-source software and publicly accessible data may be investigated and integrated where scientifically and technically justified.

Public repositories may be inspected as part of normal autonomous research.

Licensing implications must be considered before introducing dependencies, copying implementation material, or redistributing data.

The provenance of externally derived algorithms, data, or substantial implementation concepts should be documented when scientifically or legally relevant.

---

## Project Scope

Subject to the explicit implementation constraints above, the autonomous system determines:

* detailed physical models;
* transport algorithms;
* numerical approaches;
* internal architecture;
* APIs below the required Python-facing interface;
* data structures;
* memory layout;
* random-number strategies;
* material and cross-section data;
* software dependencies;
* validation strategy;
* optimization path;
* execution scheduling.

The human developer may provide high-level goals and boundary conditions but should avoid steering these choices during the experiment.

---

## Experiment Changes

Material changes to:

* implementation constraints;
* autonomy rules;
* available tools;
* available models;
* permissions;
* telemetry;
* human-intervention policy;
* reference-data access;
* Git workflow;
* external connectivity;

must be version-controlled.

For major changes, create an experiment revision or Git tag so results before and after the change can be distinguished.
