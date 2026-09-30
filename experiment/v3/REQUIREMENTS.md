# IonMC v3 requirements — inherited scientific additions (stable V2 IDs)

All original MUST requirements remain in force. Stable v1 MUST IDs and titles are
listed in requirements-index.json. These requirements specify outcomes; the agent
chooses methods, evidence and implementation order. SHOULD/MAY/out-of-scope clauses
from v1 remain unless explicitly strengthened below.

## V2-ION — MUST: independently qualified therapeutic ions

Proton, helium, carbon and oxygen-capable architecture is mandatory. Demonstrate
helium and carbon transport with ion-specific stopping/range evidence. Carbon
nuclear interactions, primary attenuation, projectile fragmentation, transported
secondary ions and distal fragment dose are central release requirements. Oxygen
must have a documented supported domain and ion-specific validation for every
claimed capability. The release cannot be proton-only. Species-resolved behavior,
LET and biological/microdosimetric extensibility must cover relevant mixed fields;
any approximations and unsupported channels require explicit domains and errors.

## V2-CAP — MUST: fail-closed capability contract

Every public simulation/dose/LET/influence path must honor requested physics or
raise an actionable unsupported-combination error before returning a result.
Nuclear, secondaries, fragmentation, LET, species scoring, biological lookup and
influence flags must not become silent no-ops. Supported species, materials,
geometries, energies, sources, scoring and backend combinations must be discoverable
and documented. Persist requested AND effective configuration. Test feature dispatch
with discriminating cases, including public planning paths, not just flag storage.

## V2-NUM — MUST: adversarial numerical falsification

Attempt to falsify each relevant approximation: decreasing step, refining/shifting
score grids, boundary alignment/misalignment, translated/outside sources, finite
phantom escape, float32/64, seed/history/batch/backend changes and beam orientation.
Explain inapplicable probes. Quantify bias and uncertainty separately. Refinement
must give sensible stability/convergence in the documented domain or reject it.
Where theoretically expected, history increases must improve statistics, including
secondary/fragment estimators. Test uncertainty validity and coverage; undefined
uncertainty must not masquerade as zero. Account separately for escaped/truncated
energy and particles where relevant; subtraction identities are not independent
conservation evidence.

## V2-EVID — MUST: physical evidence lineage and sufficiency

Important physics claims need independently justified evidence. Record the taxonomy,
provenance, shared data/model lineage, calibration/evaluation roles, observables,
normalizations, uncertainties and domains. Ion validation cannot be inferred solely
from proton data or related stopping formulas. Evidence must discriminate incorrect
physics beyond backend agreement, a plausible total integral or conservation.
Investigate contradictory references, including reference model/cut/statistical
limitations. Use the structured intervention route for genuinely inaccessible
essential evidence after bounded investigation.

## V2-BIO — MUST: external biological lookup functionality

Support externally supplied species/energy-dependent biological or microdosimetric
lookup data; include material/tissue dimensions where scientifically applicable.
Demonstrate nontrivial user-provided lookup tallies and provenance. A dot product
with the implementation's own stopping table is insufficient. No particular RBE
model or clinical certification is required.

## V2-OUTPUT — MUST: self-describing persisted scientific results

Outside the producing process, recover code version/SHA and dirty state, requested
and effective physics, species/source, geometry/material definitions, transport and
scoring grids with coordinates, units, external data identifiers/hashes, model/
reference versions, actual histories, batches, seeds/RNG identity, normalization,
uncertainty definition/validity, escape/truncation, backend and relevant environment.
Influence matrices must retain geometry, beamlet mapping, units, history and
normalization metadata. Test reopening in a separate clean process. Represent
missing/undefined quantities explicitly. Machine-readable data and human-readable
summaries/profiles/plots/reports are required; formats remain the agent's choice.

## V2-USER — MUST: complete clean-install scientific workflow

A new user must install the package, acquire/cache data, define source/geometry/
scoring, run transport and obtain dose, supported LET, uncertainty and applicable
influence output using documented public APIs/CLI/configuration. Include water,
heterogeneity, multi-beamlet and heavy-ion use where relevant. Human-readable
results must be inspectable without test code or internal implementation knowledge.
Exercise the actual documented workflow from a clean installation and offline
cache reuse; report failures as acceptance failures.

## V2-PERF — MUST: frozen application-level performance/resource gates

Ratify numeric targets for all performance.json workload families before substantial
optimization. Measure full enabled physics with meaningful accuracy/uncertainty,
wall throughput, cold start/compilation, host vs kernel time, attributable peak
host/GPU memory and scaling. Required metrics cannot be omitted because they are
hard to measure. Document measurement method, concurrency and noise. Retain failed
runs/regressions. Targets cannot silently loosen after observing implementation
performance; amendments preserve original results and rationale.

## V2-RELEASE — MUST: complete qualification and promotion

Every MUST needs explicit status and supporting suite membership. No zero-suite,
subset, stale-SHA, skipped-backend, missing-reference, missing-performance or
missing-product qualification may claim full readiness. Full CPU/reference/CUDA
local evidence and lightweight remote checks are both required. Finish the
promotion/tag/release workflow after successful qualification, preserving exact
code/evidence/condition identifiers. No unresolved MUST can be hidden by milestone
names or partial capability disclaimers.
