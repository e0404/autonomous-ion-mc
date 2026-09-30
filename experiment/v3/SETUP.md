# V3 setup and launch

V2 is archived incomplete. Active checkout: ~/aiprojects/ion-mc, based on the
pre-scientific v2 infrastructure commit dd33266073f84b630b3da390c9f762095c4a0870.
The archive retains uncommitted V2-006 work, remote branches, state and telemetry.
V3 uses v3/develop and v3/main in the same designated repository, with complete
infrastructure ancestry and no v1/v2 scientific commits. Do not rerun the historical
v2 bootstrap or use archived setup instructions to start v3.

The operator launcher ~/bin/claude-ionmc selects experiment-v3, enters ion-mc,
sets CLAUDE_CODE_PROJECT_DIR_NAME=experiment-v3, and retains operator notification
and OTLP endpoints. The shared Claude authentication profile remains unchanged.
Project settings set a separate auto-memory directory and Sonnet subagent default.
An Agent hook explicitly selects Opus for scientific roles and Sonnet otherwise;
Fable overrides require a recorded reason. See
[Claude subagent model selection](https://code.claude.com/docs/en/sub-agents) and
[hook input updates](https://code.claude.com/docs/en/hooks).

## Storage

- Cache: ~/.cache/ionmc-experiment/experiment-v3
- State: ~/.local/share/ionmc-experiment/experiment-v3
- Host runtime: cache/host-runner/venv; executions: state/host-runs
- Telemetry: state/telemetry; Codex workers: state/raw/codex
- Independent Codex review jobs/reports: state/reviews
- Reference engine archives: state/references; registered runtimes: state/engines.json
- Local validation: state/validation; memory: state/claude-memory

The launcher changes IONMC_V2_STATE/IONMC_V2_EVENT_DIR to
IONMC_V3_STATE/IONMC_V3_EVENT_DIR. Other IONMC_* names remain unchanged.
Do not copy old scientific data, reference outputs, reviews or memory into v3.
Dedicated execution software is reprovisioned; TOPAS/FRED installations can be
reused as immutable external runtimes. MCsquare is copied to a v3 runtime tree.
No caches or logs need to be deleted.

## Independent Codex review

The codex-worker MCP provides start_codex_review(task_id) and
inspect_codex_review(review_id). Start returns immediately; poll at reasonable
intervals, doing independent work while the reviewer runs. Codex uses read-only,
ephemeral execution and structured output, ignoring personal config and using
IONMC_CODEX_REVIEW_MODEL (default gpt-5.6-sol). Authentication uses the existing Codex
login. Generic codex_worker calls do not satisfy the review gate. The trusted
service writes reports outside the task sandbox and binds them to task SHA,
integration-base SHA and report hash. No tool accepts a caller-authored approval.
All important/blocking findings must be resolved. A new commit/base needs review.
The model's review remains fallible; the gate establishes execution and provenance.
See [Codex non-interactive execution](https://developers.openai.com/codex/noninteractive/).

## Reference outputs

Use list_reference_artifacts for a paginated manifest, then
materialize_reference_artifacts(run_id, task_id) to transfer all archived inputs,
outputs and logs to .ionmc-cache/reference-runs/<run-id>. Optional paths select
an exact subset. Hashes are verified and existing different files are not overwritten.
A transfer manifest records source SHA and hashes. These ignored files do not dirty
the committed task; /workspace exposes them to host validation. Analyze in scripts.
The response contains only paths/counts. Text preview is limited to 4000 bytes;
binary data cannot be read through the preview tool. Execution logs return short
tails; full logs stay archived. Do not move bulk data through model context.

## Start

Run the notification/engine/runtime preflight:

```bash
~/bin/claude-ionmc --preflight
```

It checks Codex executable/authentication/CLI capabilities, but does not consume a
review on every invocation. Setup validation must separately exercise a real Codex
review. Then start a new session, without --resume or --continue:

```bash
~/bin/claude-ionmc "Read experiment/prompts/kickoff-v3.md and carry out that autonomous run."
```

Preflight does not start scientific development. Permission notifications do not
authorize actions. Review/validation records are not interchangeable. All scientific
requirements and exact-SHA qualification/promotion gates remain in force.

## Validation write boundary

Host execution mounts committed /workspace and Git metadata read-only. Write
generated evidence to /workspace/validation/generated or
/workspace/benchmarks/generated. These are fresh protected directories for each
run, published afterwards under the corresponding ignored directory / <run-id>.
Use the returned output_paths to locate reports or stage release evidence; previous
outputs are not inputs to the next run. Temporary files use /tmp; compilation
caches use protected shared state/host-cache mounted at /cache. Materialized reference inputs remain read-only during
validation. The runner verifies the task SHA and clean state again after execution.
Generic Codex workers must use the registered V3 task ID and matching worktree;
archived and primary checkouts are rejected. Release retries always upload the
qualification report before claiming promotion complete.
