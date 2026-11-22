# Human Intervention Request IR-20261120-111407-B2596D

**Task:** V3-013
**Category:** authorization
**Status:** resolved 2026-11-22 (operator: "A: amend"; see response.json)

## Summary

Seal 14 is complete (REVIEW-827ab104 at b16e6025, pin c1cfaf9b reviewed by REVIEW-060a528a) and the ionmc samples started: 12 deterministic V-C2-D runs and the first 1e6-history sample (4411 s, all 13 pass with identity, environment and expectations equal to the seal). Checking the records against the sealed evaluator revealed a defect in the evaluator itself: warp_modules_problems() requires identical Warp kernel attestations across ALL runs of a backend, but the deterministic row (MCS and straggling off) compiles a different module set ("wp_transport_float64_diag_*", 2 modules) than the sampled row ("wp_transport_float64_plain_*" + the table kernel, 3 modules); within each row every run is identical. As sealed, the evaluator would mark all seven warp-cpu runs invalid and fail V-C1/V-C2 for a non-physical reason. The plan handles this by a prospective dated amendment with a new seal (seen samples consumed). I request authorization for that narrow amendment of the sealed evaluator; sample runs are halted until it is reviewed and re-sealed.

## Why autonomous resolution is not possible

Category 4 (authorization): the change edits the sealed evaluator after the seal, which your authorizations so far covered only for specific fixes (IR-20261118-182151-FD0238) or offered as an option (IR-20261118-201459-7CFC95, option B). Exact change requested: (a) v13_evaluate.warp_modules_problems groups runs by (backend, row) instead of backend alone, i.e. identical non-empty warp_modules across all warp-cpu runs of the same row (vc-sample; vc2d-sample) and None for every python run; the result JSON records the attestation per (backend, row); (b) because the V-C2-D seeds of k=0 have been seen (I printed the telescoped end points and their model differences while checking the runs), the amendment sets the evidence chunk to k=1 for all 18 samples (plan §C already defines k in 0..99: vc-sample seeds 20612004+100·(case−1)+1, V-C2-D seeds 20613505), i.e. the evaluator's step_expectation uses k=1 and those seeds; the 13 runs made under seal 14 stay committed as records but are not evidence; (c) tests: the rule passes with row-distinct module sets and fails on a within-row difference; seed/k expectations updated; (d) plan: dated amendment A-2 in Status, pin notes for seal 15 after its review. Nothing else changes: tolerances, acceptance criteria, sample definitions (steps_v13.py untouched), bundles, predictions. Cost: one evaluator edit + tests, manifest regeneration, two reviews, then the 18 samples again (≈ 10 h of single-process compute; the lost run is 73 min). Alternatives rejected: leaving the rule as is (false failure); running the samples first and amending afterwards (the plan requires the amendment to be reviewed before any further run, and runs at the old seal are invalid under the new one). Evidence: warp_modules_problems() on the 13 real records flags 7 of 13; per-row distinct module sets = 1 for each row. Operator note: the MCP client aborts waiting for run_host_validation after 1800 s of silence (CLAUDE_CODE_MCP_TOOL_IDLE_TIMEOUT), but the host run continues and its record appears afterwards; I poll list_host_run_records, so no action is needed unless you prefer a per-server timeout. Not a clinical claim.

## Exact requested input or action

Reply "A: amend" to authorize exactly the described evaluator amendment (group the Warp attestation rule by backend and row; evidence chunk k=1; tests; plan amendment A-2; new seal and reviews before any further sample run), or state a different disposition.

## Decision or task depending on this request

Not specified.
