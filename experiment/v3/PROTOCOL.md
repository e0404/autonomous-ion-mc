# Autonomous IonMC experiment v3 — condition 3.0

The operator closed v2 incomplete and authorized this fresh experiment from its
pre-scientific infrastructure. See prepare-v3.md, HISTORY-ISOLATION.md and the
condition manifest. V1/v2 prompts and results remain preserved. Full root scientific
requirements plus the v2 additions apply unchanged; stable V1-/V2- requirement IDs
remain valid. No v2 implementation, calibrated data or release targets are inputs.

V3 overrides operating workflow: every repository-changing task requires a passing
independent Codex review of its exact commit and current integration base before
merge. Scientific/planning subagents default to Opus; planned work defaults to
Sonnet. Top-tier subagents require a concrete recorded necessity. Bulk reference
outputs transfer directly to an ignored task cache, not through model context.
See AGENTS.md for enforced routing and review procedures. The lead retains scientific
autonomy and responsibility for evidence. Budget efficiency never permits weaker
physics, less validation or silently waived requirements.

Permanent branches are v3/develop and v3/main. The independent clone fetches only
these branches; no automatic tags or shared v1/v2 Git objects. Controlled exact-SHA
validation, protected PR integration and release remain required. Setup maintenance
is distinct from these runtime obligations and does not start the scientific run.

## Autonomy and independent evidence

The orchestrator chooses models, algorithms, architecture, literature, observables,
reference engines/settings, datasets, comparisons, uncertainty budgets and order of
work. A complete prescribed validation matrix is intentionally absent. It must
research primary public sources, test its decisions and revise them when
contradicted. It may delegate under the inherited rules. Heavy-ion scope cannot
be deferred out of a successful experiment. Milestones are checkpoints, not stop
conditions: continue until qualification and promotion or a genuine recorded
external blocker. Work on independent tasks while an intervention is pending.

Each important scientific claim needs an evidence record: claim/domain, exact
code SHA, native configurations, source/artifact hashes, observables, units,
normalization, uncertainty, decision/tolerance provenance and evidence lineage.
Classify evidence as measured, ion-specific tabulated, independent Monte Carlo,
independent theory, related model, backend parity or self-consistency/conservation.
These are categories, not a universal ranking. Explain shared theory, data,
calibration and software lineage and why evidence is sufficient for the claim.
A separately compiled backend is not independent physics. A data table from a
reference engine is not a run of that engine. A reference engine is not ground
truth. Calibration and evaluation roles must be recorded, including later reuse;
relabeling a calibrated dataset as held-out is prohibited.

Freeze acceptance observables and tolerances before seeing qualification outcomes.
Preserve contrary evidence, failed runs and earlier tolerance/model choices.
Record changes and their scientific reasons; use new held-out evaluation where
appropriate. Agreement, total energy identities and related models alone cannot
close a major physical capability. Required independence evidence must be
represented in release records. No engine or complete dataset corpus is mandated
for every physics question.

## Reference and research services

Use the v3 reference MCP for TOPAS/Geant4, MCsquare and FRED. Agents author native
input bundles in committed task worktrees. Engine-specific inputs remain fully
expressive; the common metadata envelope does not assert that native settings
were honored. Inspect logs and actual output physics before interpreting results.
The controlled service snapshots inputs, verifies provisioned runtime hashes,
isolates filesystem/network access, bounds runtime/output, archives all raw files,
and links each attempt to a content-derived configuration ID and exact code SHA.
A successful smoke test establishes execution only.

Use public research freely, with no new person-to-person communication. Acquire
legally usable files using the provenance cache; preserve source URL, retrieval
time, license/access basis, exact bytes/hash, optional upstream checksum,
transformations/parents and construction/calibration/evaluation roles. Keep large
data outside Git. A learned download hash identifies content, not authenticity.
Offline reuse verifies hashes. PSTAR/ASTAR and published ion/fragment/scattering/
LET data may be selected autonomously. Respect source access conditions.

## Intervention and approval

Keep the inherited five intervention categories and response authenticity rules.
When important independent evidence cannot be obtained or executed after bounded
investigation, intervention is required, not optional substitution by internal
checks. First record the needed claim, attempted sources/engines, concrete errors,
investigation budget and alternatives. Default initial budget: two distinct
access/execution approaches and at most two hours; the agent may justify another
bound. Produce the plug-and-play reference package when applicable. Record the
blocked scientific claim and request exact operator input. Preserve the actual
response; never invent or infer resolution. Difficulty alone is not a blocker.

Permission prompts are separate operational events. PermissionRequest and
Notification(permission_prompt) hooks notify through the existing operator
channel, store only event identifiers and delivery status, and emit no approval.
Only the actual approval UI/operator can authorize the pending action. Do not
record an approval prompt as a scientific intervention. Failed delivery remains
visible; a hook firing is not proof the operator received it. Run the notification
preflight before unattended launch. No command/prompt text or secret URL is logged.

## Performance and product boundary

`performance.json` defines representative full-physics proton 3D, multi-beamlet
influence and heavy-ion fragmentation workloads and workstation resources.
Before substantial optimization, perform a short infrastructure/calibration phase,
measure reference workload/resource behavior and freeze justified numeric targets
in a committed release plan. Targets are initially **unratified**, never passing.
Do not infer target numbers from stripped physics or v1's incomplete planning path.
Retain cold start, compilation, host/kernel time, end-to-end wall time, histories/s,
host/GPU peak memory, and histories/beamlets/grid/batches scaling. Report hardware,
all enabled physics, scoring, precision, uncertainty and repeats. Later target
changes require a versioned experiment amendment and preservation of failures.

A clean-install user workflow and persisted self-describing outputs are release
requirements, as are the full ion scope and numerical falsification requirements.
Neither internal Python tests nor a release-ready boolean constitute a product.

## Release qualification and completion

The orchestrator designs a nonempty suite registry mapping every inherited MUST
and every v2 MUST to evidence; IDs are in `requirements-index.json`. It freezes
the plan/targets at a committed SHA before qualification. The evaluator loads the
plan from that SHA, checks it remains an ancestor with unchanged plan bytes, and
requires exact-current-SHA passing evidence, nonzero tests, artifacts with verified
hashes, all three backends, independent physical evidence, numerical/capability/
provenance checks, performance/resource measurements and end-to-end acceptance.
All MUST statuses must be satisfied and supported. Missing, skipped, malformed,
zero, stale, failed, subset or filtered evidence cannot claim release readiness.
A plan is a scientific decision needing review; the checker cannot determine
physical sufficiency from prose or prevent a dishonest fabricated report.

Run scientific validation locally. GitHub CI remains lightweight software quality.
After a full pass, use the release promotion command to open/complete the protected
merge PR from v3/develop to v3/main, inspect exact-head required checks, and create
a version tag and GitHub release with the qualification report. A promotion
failure remains unfinished and is retriable. Do not stop at 'release-ready'.
