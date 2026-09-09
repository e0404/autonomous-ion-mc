# Ion Therapy Monte Carlo Software Requirements

## Purpose

This document defines the scientific and software requirements for the autonomous development experiment.

`EXPERIMENT.md` defines how the autonomous-development experiment is conducted.

This document defines what the resulting Monte Carlo system should accomplish.

The autonomous system remains responsible for determining the detailed physics models, numerical algorithms, internal architecture, data structures, dependency choices, and validation strategy, subject to the constraints and priorities defined here.

---

# Requirement Levels

Requirements use the following levels.

## MUST

Required for the project to satisfy the intended scientific scope.

## SHOULD

Strongly desired unless the autonomous system identifies and documents a compelling technical or scientific reason not to implement them.

## MAY

Potentially useful extensions that are not required for successful completion of the initial experiment.

## OUT OF SCOPE

Features that should not consume significant development effort during the initial experiment unless they become necessary to satisfy a higher-priority requirement.

---

# Application Domain

## MUST — Ion-therapy treatment planning orientation

The system is intended primarily for radiation-therapy applications, with emphasis on ion-beam treatment planning and inverse-planning workflows.

Architectural and performance decisions should prioritize efficient repeated particle transport and scoring relevant to treatment planning rather than general-purpose Monte Carlo flexibility.

The primary optimization-oriented use case is IMPT-style calculation in which many individual pencil beams or beamlets must be simulated and their contributions retained separately.

---

# Language and Execution Model

## MUST — Python-first implementation

The software must be implemented as a Python package with a native Python-facing API.

Core transport physics and interaction logic should remain readable and inspectable at the Python source level.

The implementation should avoid hiding the principal physical transport model inside opaque compiled extensions unless such components are demonstrably required and independently validated.

## MUST — Reference Python execution

A native/reference Python execution path must exist for correctness testing, debugging, inspection, and validation.

Where practical, this path should function as an independent or minimally transformed correctness reference for accelerated implementations.

## MUST — NVIDIA Warp acceleration

The system must provide accelerated execution using NVIDIA Warp.

At minimum, the intended execution targets are:

- Warp CPU;
- Warp CUDA.

The same physical model should be shared between execution paths as far as technically reasonable.

## SHOULD — Backend-extensible design

The architecture should separate physical modeling from execution sufficiently that additional execution backends could plausibly be introduced in the future.

Implementing and validating such additional backends is not required during the initial experiment.

---

# Particle Scope

## MUST — Ion transport

The system must support therapeutic ion transport.

The architecture must not assume proton-only transport.

Development may proceed incrementally by particle species.

## SHOULD — Multi-ion architecture

The design should allow extension to multiple ion species such as:

- protons;
- helium ions;
- carbon ions;
- oxygen ions.

The autonomous system determines the implementation order and detailed supported species.

## MAY — Additional particle species

Additional transported particle species may be supported where scientifically useful.

The autonomous system determines which secondary species require explicit transport and which may be treated through local or condensed approximations.

---

# Geometry and Materials

## MUST — Homogeneous reference geometries

The system must support simple homogeneous geometries suitable for analytic and reference validation.

## MUST — Voxelized geometries

The system must support voxelized geometries suitable for patient-like transport calculations.

The geometry representation must support spatially varying material properties.

## MUST — Material assignment

Voxelized geometries must support material and density assignment.

## SHOULD — Arbitrary beam incidence

Transport should support arbitrary beam directions rather than assuming axis-aligned irradiation.

## MAY — CT conversion

HU-to-material and HU-to-density conversion may be implemented.

## OUT OF SCOPE — Clinical DICOM integration

Full DICOM import/export and treatment-planning-system integration are not required for the initial experiment.

---

# Treatment Planning and Influence Matrices

## MUST — Beamlet-resolved scoring

The system must support simulations in which scored quantities can be associated with the originating pencil beam or beamlet.

## MUST — Dose influence matrices

The system must support efficient construction of dose influence matrices suitable for inverse treatment planning.

Conceptually, these represent contributions from individual beamlets to scoring locations.

The autonomous system determines the exact internal and exported matrix representations.

## SHOULD — Sparse influence representations

Efficient sparse or thresholded representations should be supported where appropriate.

## SHOULD — Batched beamlet execution

The architecture should support efficient transport of many beamlets without requiring expensive reinitialization of immutable physics or material data.

## SHOULD — Planning-oriented reuse

Data that are invariant between beamlet simulations should be reusable across repeated transport calculations.

---

# Scoring Architecture

## MUST — Changeable scoring grids

The system must allow scoring in grids of different resolution or alignment than the primary transport geometry, like a reduced resolution dose grid.

## MUST — Extensible scoring

The system must provide a general scoring architecture rather than hard-coding a single dose output.

Scoring mechanisms should be extensible to additional particle-, energy-, material-, or interaction-dependent quantities.

## MUST — Absorbed dose

Absorbed dose or the quantities required to derive absorbed dose must be supported.

## MUST — Energy deposition

Energy deposition must be available for validation and dose scoring.

## MUST — LET scoring

The system must support LET-related scoring appropriate for particle therapy.

At minimum, dose-averaged LET should be supported.

The autonomous system determines the detailed estimator definitions and implementation.

## SHOULD — Track-averaged LET

Track-averaged LET should be supported where scientifically meaningful.

## SHOULD — Fluence scoring

Particle fluence should be scoreable.

## SHOULD — Species-resolved scoring

Scoring should permit separation by particle species where relevant.

## SHOULD — Energy-resolved scoring

Energy-resolved contributions should be scoreable where required for biological or physical post-processing.

---

# Biological and Microdosimetric Data

## MUST — External biological lookup data

The scoring system must support externally supplied particle- and energy-dependent biological or microdosimetric lookup data.

This should allow derived biological quantities to be accumulated during transport without embedding one specific biological model into the transport engine.

Potential uses include:

- LET-dependent quantities;
- RBE-related quantities;
- particle- and energy-dependent biological parameters;
- microdosimetric quantities;
- nanodosimetric or ionization-cluster-related quantities.

## SHOULD — Multi-dimensional lookup dimensions

Lookup data should support scientifically relevant dimensions such as:

- particle species;
- particle energy;
- material or tissue class where applicable;
- biological or microdosimetric quantity.

The exact data model remains an autonomous design decision.

---

# Physics Model and Data Architecture

## MUST — Analytical physics layer

The system must support analytical or parameterized physical models where appropriate.

## MUST — External/tabulated data layer

The system must support physics quantities obtained from external or tabulated datasets.

## MUST — Separation of transport from data source

Transport algorithms should not require extensive redesign depending on whether a physical quantity originates from:

- an analytical model;
- a parameterization;
- a tabulated reference dataset.

The detailed abstraction mechanism is left to the autonomous system.

## SHOULD — Configurable physics models

Where meaningful, the design should permit alternative physical models or data sources to be selected for validation or accuracy/performance studies.

---

# External Data Management

## MUST — No large physics datasets in Git

Large physics datasets must not be committed directly to the source repository unless there is a compelling and documented reason.

## MUST — Automatic acquisition

Required external datasets should be downloadable automatically where legally and technically possible.

## MUST — Local caching

Downloaded datasets must be cached locally to avoid unnecessary repeated downloads.

## MUST — Dataset versioning

External datasets must have explicit versions or equivalent immutable identifiers.

## MUST — Integrity verification

Downloaded data should be validated using checksums or another appropriate integrity mechanism where possible.

## MUST — Provenance

Simulation results must retain enough metadata to identify the exact external physics data used.

## MUST — Offline reuse

Once required datasets have been downloaded, normal simulations should be able to operate without network access where practical.

## SHOULD — Configurable cache location

The user should be able to control the location of the data cache.

## SHOULD — Automated data preparation

Conversion or preprocessing of downloaded data should be automated and reproducible.

---

# Random Number Generation and Reproducibility

## MUST — Explicit random control

Random-number generation must be explicitly controllable.

Simulation configuration must permit reproducible seeding.

## MUST — Parallel-safe random streams

Random-number generation must remain statistically valid under parallel CPU and GPU execution.

## SHOULD — Deterministic backend reproducibility

Repeated simulations using the same backend, configuration, and seed should be deterministic where technically practical.

## MUST — Cross-backend statistical consistency

Different execution backends are not required to produce bitwise-identical results.

They must produce statistically compatible physical results within validated tolerances.

The autonomous system determines the detailed RNG strategy.

---

# Statistical Uncertainty

## MUST — Statistical uncertainty estimation

The system must support estimation of Monte Carlo statistical uncertainty.

At minimum, it should be possible to determine:

- number of simulated histories;
- expectation estimates;
- variance, standard error, or another justified uncertainty measure.

## SHOULD — Batch-based estimation

Batch or equivalent statistically justified uncertainty estimation should be supported.

## SHOULD — Planning-aware uncertainty

Where practical, uncertainty should be quantifiable for beamlet-resolved or influence-matrix calculations.

---

# Numerical and Physical Accuracy

## MUST — Validated physical behavior

Successful execution is not sufficient evidence of physical correctness.

Physics implementations must be validated using appropriate independent references.

## MUST — Conservation checks

Relevant conservation laws should be tested where applicable.

## MUST — Numerical convergence

Numerical approximations should be subjected to convergence or sensitivity analysis where relevant.

## SHOULD — Configurable accuracy/performance tradeoffs

Where the physics permits meaningful approximation levels, these should preferably be configurable rather than permanently fixed.

The system should permit investigation of accuracy/performance tradeoffs.

---

# Performance

## MUST — High-throughput execution

The accelerated implementations must be designed for high-throughput Monte Carlo transport.

## MUST — GPU-oriented architecture

The implementation must avoid designs that fundamentally require a Python-level control loop for every transported particle or step when running accelerated backends.

## MUST — Reproducible benchmarking

Performance measurements must record enough information to reproduce the benchmark, including:

- hardware;
- software version;
- execution backend;
- problem size;
- particle count;
- numerical precision where relevant;
- configuration;
- benchmark methodology.

## SHOULD — Efficient batching

The system should support batching strategies appropriate for GPU execution and treatment-planning workloads.

## SHOULD — Performance regression testing

Important kernels and end-to-end workloads should have performance regression benchmarks.

---

# Units and Conventions

## MUST — Explicit physical units

Physical units and conventions must be documented and unambiguous.

API boundaries must not rely on undocumented implicit unit assumptions.

## MUST — Coordinate-system documentation

Spatial coordinate conventions must be explicit.

## SHOULD — Internal consistency

A consistent internal unit system should be used unless a documented reason justifies otherwise.

The autonomous system determines the specific internal unit convention.

---

# Validation Strategy

## MUST — Independent validation

Important physics components must be validated independently from their implementation wherever practical.

Potential validation sources include:

- analytic results;
- published reference data;
- independent implementations;
- reference Python execution;
- TOPAS/Geant4 calculations requested through the experiment protocol;
- conservation laws;
- statistical consistency tests.

## MUST — Backend consistency

Reference Python, Warp CPU, and Warp CUDA execution must be compared using appropriate deterministic or statistical validation criteria.

## SHOULD — Layered validation

Validation should progress from small isolated physics components to complete transport scenarios.

---

# Software Quality

## MUST — Automated testing

The project must have automated tests.

## MUST — CI integration

Tests and required quality checks must run in CI before integration.

## MUST — Pre-commit framework

The project must use automated pre-commit checks.

The autonomous system determines the specific formatting, linting, typing, and related tools.

## SHOULD — Static typing

The public API and scientifically important internal interfaces should use useful static type information where appropriate.

## SHOULD — Documentation

Scientifically significant algorithms, assumptions, units, and interfaces should be documented.

---

# Scope Control

## MUST

The initial project must prioritize:

- correct ion transport;
- Python reference execution;
- Warp CPU execution;
- Warp CUDA execution;
- treatment-planning-oriented scoring;
- dose influence matrices;
- dose and LET scoring;
- extensible biological lookup quantities;
- external versioned physics data;
- reproducibility;
- uncertainty estimation;
- scientific validation;
- high-performance execution.

## SHOULD

The project should support:

- multiple ion species;
- configurable physics models;
- sparse planning representations;
- efficient batched beamlet execution;
- extensible scoring;
- configurable physics fidelity.

## MAY

The project may include:

- CT calibration;
- phase-space input/output;
- multi-GPU execution;
- distributed computation;
- visualization utilities;
- additional execution backends;
- advanced workflow tooling.

## OUT OF SCOPE

The initial experiment should not prioritize:

- graphical user interfaces;
- full DICOM workflow integration;
- clinical TPS integration;
- regulatory or clinical certification;
- hospital deployment;
- optimization of unrelated general-purpose Monte Carlo use cases.
