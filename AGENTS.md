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

## Human intervention

When human intervention is permitted by `EXPERIMENT.md`, use the `request_human_intervention` tool rather than asking informally in the conversation.

The tool creates the authoritative intervention record and notifies the designated experiment operator.

Do not use the tool merely because a scientific or technical decision is difficult or uncertain.

When the experiment operator provides the requested information or action, record that response with `resolve_human_intervention`.

Do not fabricate, infer, summarize away, or resolve an intervention without an actual operator response. Preserve the substance of the operator's answer in the resolution record.


## Task integration workflow

Development tasks should use the controlled task lifecycle.

For a task that changes the repository:

1. create an isolated task worktree;
2. delegate or implement the work there;
3. run the required tests and local validation;
4. review consequential changes;
5. commit the completed work on the task branch;
6. push the branch with `push_task_branch`;
7. create the pull request with `create_task_pull_request`;
8. inspect the pull request and required checks;
9. squash-merge it into `develop` with `merge_task_pull_request`;
10. retire the clean task worktree.

Do not push directly to `develop` or `main`.

Task branches are preserved after integration.

### GitHub CI observability

Agents must inspect GitHub CI through the controlled task-integration tools rather than attempting to access operator GitHub credentials.

Use `inspect_task_ci` to inspect check, workflow, job, and step status for the exact pushed task-branch SHA.

If CI fails, use `get_task_ci_failure_logs` to retrieve bounded failed-job logs and diagnose the failure.

CI observations must correspond to the exact current pushed task SHA. A local commit that has not yet been pushed invalidates the relevance of earlier CI results.

Do not treat `mergeStateStatus` alone as evidence that the required CI checks passed when check-level information is available.

## Autonomous Git commits

Autonomous agents must not depend on the operator's personal Git identity or attempt to read or modify user-level Git configuration.

Use `commit_task_changes` to commit completed task-worktree changes.

Autonomous commits use the fixed experiment identity:

    Autonomous IonMC Agent <autonomous-ionmc-agent@users.noreply.github.com>

Before committing:

- inspect the worktree changes;
- ensure no unintended files are present;
- run the appropriate tests for the task.

Local scientific validation records are created only after the exact state being validated has been committed.


## Local validation gate

A task may only be merged into `develop` after appropriate local validation has been performed for the exact committed task-branch SHA.

After committing the state that was validated, record the validation result with `record_local_validation`.

A new commit invalidates prior validation for merge-gate purposes because validation records are keyed to the exact Git SHA.

Do not record a passing validation result unless the stated validation was actually performed successfully.

The required depth of validation depends on the change:

- documentation or infrastructure changes may require focused functional validation;
- software changes require appropriate automated tests;
- scientific or numerical changes require the relevant local scientific
  validation defined by `EXPERIMENT.md` and `REQUIREMENTS.md`.

GitHub CI is not a substitute for the local validation gate.

## Lightweight GitHub CI

GitHub CI is a software-quality gate, not a scientific-validation environment.

Before integration, the lightweight CI checks must pass. These checks may cover repository hygiene, syntax, unit tests, packaging, and other small deterministic tests.

Do not move full scientific validation, large Monte Carlo runs, GPU-intensive validation, reference comparisons, or performance validation into GitHub-hosted CI. Those remain local validation tasks under the exact-SHA local validation gate.

## Controlled local host/GPU execution

Agents must not bypass the Claude sandbox to obtain GPU or unrestricted host execution.

For committed task states requiring local CPU/GPU execution, use `run_host_validation`.

The runner executes the requested argv directly, without a shell, inside a second bubblewrap sandbox with:

- the task worktree mounted at `/workspace`;
- dedicated experiment caches mounted at `/cache`;
- WSL GPU access through `/dev/dxg`;
- no network access;
- no Docker socket;
- no user home directory;
- no credentials;
- no arbitrary Windows filesystem access.

Host validation requires a clean task worktree and is associated with the exact current commit SHA.

A successful command is execution evidence only. It does not constitute scientific validation by itself. The orchestrator must separately record an appropriate exact-SHA local validation result using the validation manager.

Host-runner results include bounded stdout and stderr directly in the structured MCP response so agents can interpret scientific and diagnostic output without reading the protected raw-run archive.

Complete stdout and stderr remain archived outside the task worktree. If an inline stream exceeds the configured limit, the MCP result returns the tail of the stream together with truncation and total-length metadata.

## Public development-network access

Sandboxed development tools may autonomously access public network resources needed for software and scientific development.

Common development domains configured in `.claude/settings.json` are pre-authorized so routine operations such as public Git access, dependency retrieval, and pre-commit hook installation do not require additional network approval.

These configured domains are not a strict egress allowlist. Other public network destinations may be permitted by Claude Code's normal sandbox and Auto-mode permission mechanisms.

Public web research is separately permitted by `EXPERIMENT.md`.

Agents must not bypass the sandbox or use network access to perform activities prohibited by `EXPERIMENT.md`, including initiating human communication.
