# Operational instructions for experiment v3

Read the root protocol/requirements and this directory's PROTOCOL.md,
REQUIREMENTS.md, HISTORY-ISOLATION.md, performance.json and requirements-index.json.
The v3 kickoff is experiment/prompts/kickoff-v3.md. Setup maintenance is distinct
from the autonomous runtime. No inherited scientific MUST is relaxed.

Use v3/develop and v3/main, task IDs V3-*, and controlled MCP tools for commits,
validation, CPU/GPU execution, GitHub checks/integration, references and release.
Do not access the v1/v2 archives or scientific history. Direct gh remains denied.

For every changed task: create worktree, delegate bounded work, validate, commit,
start_codex_review, inspect_codex_review, fix/recommit/review until passed, record
exact-SHA local validation, push, open PR, inspect exact-head CI, merge, retire.
Review MUST be independent Codex execution on the exact clean committed task and
current integration base. Missing, stale, failed or unresolved reviews block
merge. Do not write protected review records or substitute self-review. Review
release criteria and plans before integration. Scientific evidence remains a
separate gate: a Codex approval is not scientific validation.

Model routing: Opus for scientific research, task planning and difficult bounded
analysis; Sonnet for implementation of agreed plans, reference execution and
retrieval, parsing, tests, docs and routine operations. General-purpose workers
default to Sonnet, never inherit Fable. The Agent hook pins missing/inherit models
according to role; top-tier Fable requires a substantive
TOP_TIER_JUSTIFICATION: line explaining why lower tiers cannot resolve the task.
Codex planned implementation defaults to gpt-6-luna; independent review uses
operator-configured gpt-6-sol. Use small parallel groups when tasks are independent.
Avoid delegation overhead for trivial commands. Preserve small checkpoints and
compact context between milestones; do not reread every document on every call.

For reference output, list_reference_artifacts supplies exact paths, sizes and
hashes. materialize_reference_artifacts copies all/selected files to
.ionmc-cache/reference-runs/<run-id> without putting contents in conversation.
Use local analysis scripts and bounded summaries. read_reference_artifact is a
UTF-8 preview capped at 4000 bytes; binary data requires materialization. Do not
reconstruct base64 strings in model output. Raw logs/evidence remain archived.

Keep a requirement ledger, freeze criteria before seeing outcomes, document
independence/calibration lineage, record contrary evidence and failures, and use
run_reference_calculation, acquire_reference_data and record_scientific_event.
Missing essential external evidence after bounded investigation requires authentic
structured intervention. A permission notification is not operator approval.
Continue through task and milestone boundaries to qualification and promotion;
never represent a partial implementation or subset of evidence as completion.
