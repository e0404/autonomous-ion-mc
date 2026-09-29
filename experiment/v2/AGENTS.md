# Operational instructions for the v2 autonomous orchestrator

Read the inherited root AGENTS.md/EXPERIMENT.md/REQUIREMENTS.md, then v2 PROTOCOL.md,
REQUIREMENTS.md, performance.json, requirements-index.json, HISTORY-ISOLATION.md and the kickoff-v2-isolated prompt.
V2 overrides take precedence only where stated. The setup maintainer is distinct
from the autonomous runtime governed here.

Work in the independent v2-only clone with .ionmc-condition.json, never the setup
branch's completed v1 package. Permanent branches are v2/develop and v2/main.
Use the existing task manager, commit, local validation, host runner and controlled
CI/integration MCPs. Use inspect_task_ci and bounded failure logs for the exact
pushed SHA; direct gh CLI is denied in the v2 agent configuration. The trusted
integration/release services retain their narrowly scoped GitHub access.

Use run_reference_calculation/list_reference_engines/read_reference_artifact for
reference execution, acquire_reference_data/assign_reference_data_role for data,
and record_scientific_event for consequential validation changes/failures. Inputs
are committed case bundles; arbitrary engine-native inputs are allowed inside the
controlled runner. Explain suitability and independence in evidence records.

Plan/freeze acceptance criteria before results. Maintain the requirement ledger,
release plan, scientific decisions and validation records continuously. Reject
insufficient internal evidence explicitly. Investigate contradictory references.
Record failed comparisons and model/tolerance changes with links to artifacts.

Missing essential external evidence after bounded investigation requires
request_human_intervention with reproducible files, attempts, limits and exact
requested operator action. Never fabricate resolution. Permission notification is
not permission approval and is not scientific intervention.

Task completion and milestones are checkpoints. Continue autonomously across them,
including research, independent review, integration and release promotion. Stop
only for successful promoted release or a genuine external blocker with a recorded
intervention; continue unrelated work if possible. Do not start with or copy v1
scientific implementation/validation. Do not fetch, browse or reconstruct v1
branches, tags, commits, PRs, logs or audit archives, including through web tools,
MCP tools or other local checkouts. Only the supplied v2 requirements and audit
lessons are experimental inputs. Do not change the restricted Git fetch settings
or attach another object database. No clinical-use claims.
