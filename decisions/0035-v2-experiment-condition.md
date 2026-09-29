# 0035 — Separate v2 experimental condition and evidence infrastructure

Context: v1 final/audited SHA is 3a2e8953be0034d75a29ccf2b04640353759982b.
The independent audit was read from the named workspace archive; its report/archive
hashes and historical starting references are recorded in experiment/v1/provenance.json.
The user requested experiment preparation, not a v1 physics bugfix campaign.

Decision: preserve root v1 protocol/requirements/kickoff bytes and all historical
commits. Define an additive v2 condition under experiment/v2. Reconstruct a fresh
scientific tree from the recorded pre-scientific start plus explicit infrastructure
paths. Separate v2 permanent branches protect v1 ancestry and keep existing task
integration useful with portable roots/base-branch routing.

References remain native-input engines behind a shared provenance/execution
interface. Use installed TOPAS/Geant4 and licensed FRED; use an exact public
MCsquare binary/source revision rather than silently replacing an unavailable
classic-compiler source build. Content hashes and raw logs identify what ran.
Bwrap execution excludes network/home/credentials, accepts only trusted provisioned
runtime mounts, and never falls back to unsandboxed execution for the future agent.
Runtime smoke fixtures do not define physical qualification or a complete corpus.

Scientific autonomy is preserved: evidence categories and lineage review replace
a fixed validation matrix. Heavy ions, fragment fields, falsification, capability
contracts, usable outputs and workflow acceptance become explicit release outcomes.
Workload envelopes and observed A6000 resources are fixed now; realistic numeric
targets must be ratified/frozen in a bounded early calibration phase. No fabricated
performance measurement or copied audit tolerance is treated as a passing target.

Release qualification is a separate infrastructure evaluator, not a patch to v1's
release validator. It rejects missing/empty/filtered/stale/uncorrelated evidence,
checks hashed artifacts, complete MUST membership, backends, independence and
performance/product gates. A trusted promotion service re-evaluates evidence before
remote release operations. Scientific review remains necessary; schemas cannot
prove that evidence is true or sufficient.

Permission notifications reuse the existing operator transport with generic text;
they never approve actions or resolve scientific interventions. Missing notification
configuration must remain visible in preflight. Direct gh denial is defense in
depth; actual credential separation and controlled services remain the boundary.

Validation: focused adversarial infrastructure tests, actual native reference
smokes, documentation/pre-commit build and a generated fresh-start checkout.
Measured outcomes and remaining environment limitations are in
experiment/v2/setup-validation.json. No scientific v2 implementation is begun.
