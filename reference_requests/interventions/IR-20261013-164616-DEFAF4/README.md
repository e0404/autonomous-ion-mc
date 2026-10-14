# Human Intervention Request IR-20261013-164616-DEFAF4

**Task:** V3-005B
**Category:** experiment-boundary
**Status:** resolved (2026-10-13T18:37:00.866542+00:00; operator response in `response.json`, recorded through `resolve_human_intervention`; the README status line was reconciled in the task branch after the response, the request tooling does not rewrite it)

## Summary

V3-005C qualification is complete and recorded at the evidence SHA dd976c5a (PR #82, head bb8c198f, Codex REVIEW-e05c771b pass, CI green at the head): the pre-registered scientific gates V5, V6, V7-R and V10-A failed and are preserved; the implementation gates (P5-ext, P6-D, P7, C1-ext, R1-D2, V8-LV, CI) and V8-nuc pass. Merge gate Amendment 14 (a)6 is not met; no passing record exists; the default is that PR #82 stays unmerged. As for V3-005B (IR-20261009-032741-C590CF, interpretation B), the question is only how the recorded failures are carried forward; nothing here proposes a relaxation, re-labelling or replacement of any gate.

## Why autonomous resolution is not possible

Facts at dd976c5a (lv5c base 20481004, hr5c base 20491004, single-process diagnostic mode; archives validation/results/transport/{lv5c,hr5c}-dd976c5-single-process.json; records lv5c VAL-20261013-163425-94C596 and hr5c VAL-20261013-163444-63C7B0, both failed, re-stating VAL-20261013-061553-2FEEC7 / VAL-20261013-061607-1BA988):
- V5: with the hadronic elastic channel the TOPAS peak/plateau discrepancy moved from +12.5/+12.3 % (lv5b at 7aae5bb1, no elastic) to −3.7/−4.9 % (150/200 MeV; tolerance 2 %); ΔIDD plateau-integral +10.3 % at 150 MeV fails (tolerance 10 %), +6.7 % at 200 MeV passes; plateau, R80 and total deposit pass; the attribution stop rule triggered (rung 1: p+A angular shape); the pre-registered expectation was a ratio of 0.985/0.986, observed 0.963/0.951.
- V11: ionmc elastic factor F = 0.857/0.844 vs TOPAS 0.892/0.887 (tolerance 0.02) and plateau gain +0.006/+0.008 (tolerance 0.005): ionmc's elastic effect is about 30 % stronger than g4h-elastic_HP, consistent with the V5 over-correction.
- V7-R: all three coverage gates pass (escaped_neutral 0.685, sec_p 0.674, nuclear_local 0.673 within [0.6527, 0.7127]), but the pre-registered calibration misses: paired gate 4.69 % at its margin on the archive-fitted surrogate (bound 4 %); on the real rehearsal block sums the paired and single gates accept 11–12 % at their margins and the joint rate exceeds 2 % at two margins; power 0.95 → "uncertainty-coverage evidence not established".
- V10-A: deterministic failure of clauses (a) at 15°/20°, (b) and (c); the (e)9 fallback trigger did not fire.
- V6/V6-sens: failed — not sized under the pre-registered rule (Amendment 14a (a)4; acknowledged by the operator in IR-20261011-031924-C205C8); the follow-up row V6-F is pre-registered and independently reviewed (Amendment 18, validation/plans/v3-005c-v6f-followup.md) for execution once in a later task.
- V8-nuc: T12 conclusive pass (p = 0.50, 160 groups). V7-scan-el: sec_p slope −0.5526 vs [−0.55, −0.45] (fail, marginal). R1-D2, V2b-el, V7-shift-el, E1-B, P7, V10, X-ENDF, D9-BGG pass; E-shape executed.
- Both archives are non-conformant by code (single-process directive, imported partials, deferred multiprocessing step). Incidents recorded in the README: a superseded partial run at 37e87ecf, a restart that zeroed four published copies, three hypervisor/WSL crashes (operator diagnostic rules applied; logs in telemetry/test-logs), one interpreter anomaly in a calibration shard (clean re-run), two wrapper defects executed directly at the same SHA and fixed in the closure, CI at 4ce33832 failing only for the missing data cache (fixed at bb8c198f).
- The PR body cannot be updated through the available tooling (create_task_pull_request refuses an existing PR; no update tool; direct gh is forbidden); the description of record is the README section "V3-005C results" and decision 0041.
Default without an operator decision: PR #82 stays unmerged; no gate is changed; any answer is recorded verbatim and changes no row verdict.

## Exact requested input or action

Choose one (or state another): (A) keep PR #82 unmerged (default) and open follow-up tasks on v3/develop from the branch state: V3-005D (p+A elastic angular shape / σ_el scale attribution and re-qualification of V5/V11 at the declared fresh bases 20511004/20521004), V3-005E (V7-R gate redesign: margin discrimination of the bootstrap-t and paired gates at ±0.03 under heavy-tailed block sums), V3-005F (V6-F execution per Amendment 18), each carrying the recorded failures in the ledger; (B) merge the branch to v3/develop as an explicitly non-qualified integration (implementation gates and CI pass; scientific gates recorded failed; no claim of validation) under a new reviewed amendment that documents the exception — the orchestrator does not recommend this and will not create a passing record; (C) another disposition. Any answer is recorded verbatim and changes no row verdict.

## Decision or task depending on this request

Not specified.

## Resolution (recorded 2026-10-13, reconciled in the task branch on 2026-10-14)

Operator response (`response.json`, `resolved_at_utc` 2026-10-13T18:37:00.866542+00:00, verbatim):

> Proceed with option A: keep PR #82 unmerged and continue dependent follow-up work from its unmerged branch state. Preserve all recorded failures and do not create an implementation-only passing record.
Prioritize elastic-model attribution and V7-R calibration repair before another full qualification run. Before expensive qualification, document and independently review the route to integration, including remaining requirement coverage, V6/V6-F treatment, and execution conformity under the existing stability restrictions.
V6-F passing does not authorize replacing V6 or merging. Bring back any concrete proposal that changes required validation scope; routine scientific investigation and implementation should continue autonomously.

Disposition A applies: PR #82 stays unmerged; follow-up work continues on the V3-005B branch from its unmerged state (commit 749f768c records the response; the planning documents under `validation/plans/v3-005d-*` and Amendment 19 follow it). All recorded failures are preserved; no implementation-only passing record is created; the proposals that change the required validation scope are collected in `validation/plans/v3-005d-integration-route.md` §4.1 (iv) for one later intervention.

Operator note stored with the response: Operator response received 2026-10-13 after CI went green at bb8c198f. Disposition A: PR #82 stays unmerged; follow-up work branches from its unmerged state; recorded failures preserved; no implementation-only passing record; priorities elastic-model attribution (V5/V11) and V7-R calibration repair; an independently reviewed integration-route document (requirement coverage, V6/V6-F treatment, execution conformity under the stability restrictions) precedes any expensive qualification; V6-F passing never replaces V6 or authorises merging; scope changes return as proposals; routine investigation continues autonomously.
