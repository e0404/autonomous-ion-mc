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
ephemeral execution from a protected checkout of the exact Git object (without
live ignored task inputs) and structured output, ignoring personal config and using
IONMC_CODEX_REVIEW_MODEL (default gpt-5.6-sol). Authentication uses the existing Codex
login. Generic codex_worker calls do not satisfy the review gate. The trusted
service writes reports outside the task sandbox and binds them to task SHA,
integration-base SHA and report hash. No tool accepts a caller-authored approval.
All important/blocking findings must be resolved. A new commit/base needs review.
The model's review remains fallible; the gate establishes execution and provenance.
See [Codex non-interactive execution](https://developers.openai.com/codex/noninteractive/).

## Reference outputs

Reference-data acquisition accepts downloads up to 1 GiB (1,073,741,824 bytes)
per file. The operator raised this cap for IR-20261019-000144-2C50C8 to allow
the ENDF/B-VIII.0 neutron archive. Restart the v3-research MCP server after
updating its source so the running process loads the new cap.

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

Host execution creates a protected snapshot from the recorded Git commit with an
independent Git database, then mounts it at /workspace read-only. Live worktree
edits cannot change the running source or its Git provenance. Ignored .ionmc-cache
inputs are copied as regular files; committed symlinks are rejected. Write
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

The host/reference runners and preflight require bubblewrap >= 0.12.0, including
its [CVE-2026-87766 fix](https://github.com/containers/bubblewrap/security/advisories/GHSA-pxhw-h44j-8pfx).
The dedicated binary is state/runtime/bin/bwrap; the launcher prepends that directory
to PATH for Claude as well. The system binary is unchanged.

## Usage-limit alerts

The StopFailure hook routes rate_limit errors to the existing IONMC_NTFY_URL
channel. No additional environment variables are needed. Alerts contain the
experiment and hashed session identifier, with best-effort quota category;
inspect Claude for reset details. Duplicate successful alerts are suppressed for
15 minutes per session, and failed deliveries have a 60-second retry cooldown.
These operational interruptions are logged separately from scientific interventions.
Start a fresh Claude session to load the updated hook configuration. Preflight tests
notification delivery and runtime readiness; it does not check remaining Claude
subscription quota or start the scientific experiment.

## Signed attestations (V3-INFRA-002)

Records written by the protected MCP services carry an Ed25519 `attestation`
block so that a repository-authored copy can be told apart from a real record.

- **What is signed:** the Ed25519 message is the canonical JSON (sorted keys,
  compact separators, ASCII) of `{"payload": <record without its attestation
  field>, "signer", "signed_at_utc", "key_id"}`; `payload_sha256` hashes the
  canonical payload alone. `validation-manager` signs each `VAL-...` record
  before it is written and returned; `host-runner` signs each `RUN-...` record
  (including `stdout_sha256`, `stderr_sha256` and `output_manifest_sha256`) and
  stores it as `record.json` (and `result.json`) in the run directory.
- **Block shape:** `{"schema": "ionmc-attestation-1", "alg": "ed25519", "signer",
  "key_id" (sha256(public key)[:16]), "public_key" (hex), "signed_at_utc",
  "payload_sha256", "signature" (hex)}`.
- **Keys:** created on first use under `<share root>/attestation/<signer>.key`
  (32-byte seed, hex, mode 0600) with `<signer>.pub` (public key, hex) next to
  it. Implementation is pure stdlib (`infrastructure/experiment_v3/attest.py`,
  RFC 8032, not constant time; acceptable for this single-user host).
- **Trust file:** `infrastructure/keys/trusted_attestation_keys.json`, empty on
  merge (trust on first use). A record whose key is not listed is `unverified`,
  never `verified`; a present-but-wrong attestation (altered payload, bad
  signature, wrong signer) is `invalid`.
- **Enforcement rule (local-validation gate):** `infrastructure/validation/
  local_validation.py` (`records_for`, `inspect`, `check`) verifies every record
  against the trust map with expected signer `validation-manager` and requires
  the signed `task_id` and `git_sha` to equal the requested task and SHA (else
  `invalid`). Each record gets a derived `attestation_status`
  (`verified`/`unverified`/`invalid`, plus `attestation_reason`). Only
  `verified` records count as passing; unverified and invalid ones are rejected
  with an explicit reason. `invalid` records are rejected always. **Bootstrap
  exception, migration task only:** `BOOTSTRAP_TASKS = {"V3-INFRA-002"}`; only
  for such a task, while the trust map has no `validation-manager` key, are
  unverified records tolerated ("unverified (no trusted key published yet)").
  For every other task unverified records are rejected whatever the trust map
  state ("attestation required: no trusted validation-manager key published
  (operator runbook) or record not verified"). **After this change merges,
  every other task's merge gate is blocked until the operator publishes the
  keys** (the runbook does both keys in one sitting). Records with status
  `failed` always count as failing.
- **Host-run evidence:** a validation record references host runs with
  `details.host_runs = [{"run_id": "RUN-...", "output_manifest_sha256":
  "<64 hex or null>"}]`. When `record_local_validation` receives entries with a
  `run_id` but no `output_manifest_sha256`, the service looks up the signed run
  record under the protected host-runs root and fills the digest in itself (the
  orchestrator cannot mis-state it); if the run record is missing, the record
  is still written but with `details.host_run_manifest_lookup: "missing"` and
  the gate rejects it. The gate accepts the record only if the signed run record
  `<share root>/host-runs/<task>/<RUN-id>/record.json` exists
  (`infrastructure/validation/host_run_records.py`), is `verified` by a trusted
  `host-runner` key (unverified tolerated only in bootstrap, migration task
  only; missing or invalid never), has signed `task_id`, `sha` and `run_id`
  equal to the validation record's task, SHA and the referenced id (replay
  protection), has a signed `output_manifest_sha256` equal to the claimed one
  (null must match null; mismatch rejects), and succeeded: `succeeded` true,
  `exit_code` 0, `timed_out` false, `committed_state_unchanged` true and
  `dirty_after` false (when present), and no `output_publication_error`. The
  legacy `host_run_id` / `host_run_ids` forms carry no manifest claim and are
  rejected with "host run referenced without output_manifest_sha256 binding"
  except for the bootstrap task. `inspect` reports each record's
  `host_run_attestation`.
- **Transport consumers:** `validation/scripts/transport/attest_verify.py`
  provides `require_verified(record, expected_signer, trust=None, *,
  expected_task_id, expected_git_sha, expected_run_id,
  expected_output_manifest_sha256, bootstrap=False)`. It is strict by default
  (only `verified` passes; invalid always raises), raises on any mismatch of a
  given expected identity value, and tolerates unverified records only with an
  explicit `bootstrap=True` while the signer has no trusted key. The gate above
  is the production consumer of `RUN-` records.
- **Operator runbook:** merge the PR; restart the MCP servers; run
  `python infrastructure/experiment_v3/attest.py --print-public-key validation-manager`
  and `... host-runner`; add each `key_id`/`public_key` to the trust file in a
  reviewed commit. Records created before the restart carry no attestation and
  stay `unverified`; once the `validation-manager` key is published they no
  longer count toward the gate, so re-record validation after publishing the
  key. Until both keys are published every non-migration task's gate is blocked.
