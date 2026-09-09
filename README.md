> **Status:** Active autonomous-development experiment. External contributions and unsolicited technical input are not incorporated during the experimental phase.

# Autonomous Ion Monte Carlo

Experimental development of a research-grade, high-performance Monte Carlo simulation framework for ion therapy using autonomous coding agents.

The software is intended to provide:

* a Python-first reference implementation;
* NVIDIA Warp CPU and CUDA acceleration;
* treatment-planning-oriented scoring, including dose influence matrices and LET;
* extensible physics, scoring, and biological-data interfaces;
* reproducible validation and benchmarking.

This repository is also the experimental record for studying autonomous scientific software development.

## Autonomous development infrastructure

This repository is developed as an experiment in highly autonomous scientific software engineering.

The development environment includes:

* a lead Claude Code orchestrator with project-scoped specialist subagents;
* an external Codex worker exposed through a controlled MCP interface;
* isolated Git worktrees for delegated tasks;
* structured experiment telemetry and preserved model interaction data;
* a controlled human-intervention channel for the limited cases permitted by `EXPERIMENT.md`;
* GitHub pull-request based integration, with task branches squash-merged into `develop` and releases merged from `develop` into `main`.

The autonomous system is responsible for scientific, numerical, architectural, and implementation decisions within the constraints defined by:

* `EXPERIMENT.md` — experiment protocol and autonomy rules;
* `REQUIREMENTS.md` — scientific and software requirements;
* `AGENTS.md` — operational instructions shared by coding agents.

Infrastructure supporting orchestration and experiment recording is located under `infrastructure/` and `.claude/`.

Human intervention is intentionally restricted and must use the structured intervention mechanism defined by the experiment protocol.


