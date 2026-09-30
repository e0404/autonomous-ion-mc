# V3 evidence and telemetry

The append-only JSONL event stream is in the experiment's v3 state directory,
separate from v1 telemetry. Each event has schema/condition, UTC time, exact code
SHA when applicable, task and structured details. Reference runs automatically
emit started/finished/failed events with content configuration ID and attempt ID;
data acquisition and role changes link content hashes. Raw reference requests,
inputs, outputs, stdout and stderr are retained under references/REF-*/.

The scientific-event tool accepts these events with a required rationale and
artifact record IDs: validation_strategy_changed, internal_evidence_rejected,
scientific_failure, reference_contradiction, model_changed, tolerance_changed,
performance_failure and workflow_failure. Include old/new decision IDs and
holdout consequences where applicable. Never remove a failed attempt. This
allows analysis of changing strategies and evidence sufficiency rather than only
successful runs. The protocol requires these records; event absence cannot prove
absence of failures and the service does not infer unreported scientific meaning.

Interventions retain the existing authoritative requests and operator resolutions;
v3 events reference those records. Permission prompts and delivery status are
separate event types. Hook notifications include generic text and a random event ID,
not raw commands, prompt bodies, transcripts, environment variables or channel URLs.
Notification delivery failures are retained; no approval decision is emitted.

Dataset objects preserve exact bytes with SHA256 and acquisition metadata. Every
use records role: construction, calibration, evaluation or exploratory. Evaluation
records disclose related uses; hashes prevent renaming a calibration file to hide
its identity. Scientific independence still requires reasoned review, not merely
a different filename, engine name, or hash.

V3 records agent_routed events with selected role/model and any top-tier rationale.
Protected state/reviews holds independent Codex job status, exact task/base SHAs,
structured findings and raw JSONL usage. General Claude SubagentStart/Stop hooks
do not represent Codex execution; inspect the review records separately.
