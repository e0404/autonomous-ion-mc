# Agent Instructions

This repository is part of an autonomous scientific software-development experiment.

Before planning or modifying the project, read:

1. `EXPERIMENT.md` — rules governing the autonomous-development experiment.
2. `REQUIREMENTS.md` — scientific and software requirements for the Monte Carlo system.

The central autonomy rule is:

> Scientific or technical uncertainty is not a reason to ask the human developer to make the decision.

Research the issue, choose an approach, document consequential decisions, implement it, and validate it.

Human intervention is limited to the cases defined in `EXPERIMENT.md`.

Maintain:

- `decisions/` for consequential scientific and technical decisions;
- `reference_requests/` for structured requests requiring human action;
- `validation/` for validation methods and results;
- `benchmarks/` for reproducible performance evaluation.

Do not expose, copy, or commit credentials or unrelated data from the host environment.

Do not weaken the autonomy, validation, implementation, or scientific requirements without an explicit, version-controlled experiment change.