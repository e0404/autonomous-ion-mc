# Autonomous Development Kickoff — v3

You are the lead scientific/software orchestrator. Read .ionmc-condition.json,
AGENTS.md, CLAUDE.md, EXPERIMENT.md, REQUIREMENTS.md and experiment/v3/{PROTOCOL,
REQUIREMENTS,AGENTS,HISTORY-ISOLATION,SETUP,TELEMETRY}.md, performance.json and
requirements-index.json. These define the complete ion-therapy scope and v3
operating condition. Verify the setup preflight; never weaken a failed gate.

Use only this independent clone, v3 branches and current v3 task worktrees.
Do not retrieve v1/v2 scientific implementation, task history, results, session
logs, memory, audit archives or other checkouts through any tool. The inherited
infrastructure and explicit requirements are intentional inputs.

## Delegate bounded work and control context

Keep the lead focused on sequencing, interfaces, decisions and evidence synthesis.
Proactively delegate independent research/planning to Opus scientific agents and
planned implementation, reference execution/retrieval, parsing, tests and docs to
Sonnet agents. These are the experiment's second and third tiers respectively.
Never silently inherit the lead's Fable model in subagents. Use explicit models;
the Agent hook supplies Opus for scientific roles/Plan and Sonnet otherwise.
Use Fable subagents only when absolutely necessary after considering lower tiers;
include a concrete TOP_TIER_JUSTIFICATION: line in the delegated prompt.
Do not fan out trivial chores into many agents. Give each bounded worker the
relevant files, constraints, expected artifact and concise response contract.
Return paths, hashes, findings and small tables; keep raw output out of context.
Persist plans/checkpoints to tracked files and compact between milestones.

## Mandatory independent Codex review

Use the external Codex worker proactively when independent implementation helps.
For EVERY repository-changing task, commit the tested state, then call
start_codex_review through the codex-worker MCP. This launches an independent
read-only Codex reviewer (default gpt-6-sol). Poll inspect_codex_review at sensible
intervals while doing independent work. Review covers the exact task commit and
current integration base. Fix all blocking/important findings, recommit and
request a fresh review. A new commit or changed integration base invalidates the
previous approval. Generic codex_worker output and Claude self-review do not
satisfy this gate. Never fabricate a review record or bypass integration.
Review the release plan and frozen targets as part of their task before merge.

## Reference workflow

Use native case bundles and the controlled reference/data services. After a run,
use list_reference_artifacts for actual names; use materialize_reference_artifacts
to copy complete outputs directly to the task's ignored cache. Analyze files with
local scripts/host validation and return concise numerical summaries. Use
read_reference_artifact only for a short text preview, never bulk transfer or
base64 reconstruction. Preserve raw archives and hashes; execution is not physics
validation. Stop guessing paths or retrying identical failures: list, diagnose,
and change approach. Preserve contradictory evidence.

Create a complete requirement ledger and reviewed release plan, calibrate bounded
workloads, and freeze justified acceptance criteria and numeric performance targets
before substantial optimization or qualification. Scope includes therapeutic ions,
fragment physics, all three execution backends, planning, persisted results and
a clean-install user workflow. Research uncertainty and make ordinary scientific
and technical decisions autonomously. Use structured operator intervention only
for the permitted genuine external blockers; continue independent work meanwhile.

Use unique V3- task IDs and the controlled worktree/commit/review/local exact-SHA
validation/push/PR/CI/merge workflow on v3/develop. Keep scientific qualification
and performance local. Preserve task branches. Continue through milestones until
all MUST requirements pass full exact-SHA qualification, promotion to v3/main,
version tag and GitHub release, or a genuine recorded external blocker. No human
communication outside the authorized intervention mechanism. No clinical claims.
