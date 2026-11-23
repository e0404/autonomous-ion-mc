# Human Intervention Request IR-20261118-182151-FD0238

**Task:** V3-013
**Category:** authorization
**Status:** resolved 2026-11-18 (operator: "A: apply the evaluator cache_dir fix"; see response.json)

## Summary

A permission decision is needed for one edit of the two sealed V3-013 validation scripts (validation/scripts/stopping/steps_v13.py and v13_evaluate.py). Codex REVIEW-f5d6f67e (at aed1f072) found that the cache-directory normaliser I added rewrites the analytic model's null cache_dir to the "<run-cache-dir>" placeholder and that the evaluator then requires a real path for every sealed block, which would invalidate every default-model (bethe-ionmc-v1) evidence run. The implementation worker's edit of the two files was denied by the auto-mode permission classifier ("Undermining Oversight"). Following the rule that a subagent's permission denial must not be bypassed by the orchestrator, I have not applied the edit myself and the seal is held until you decide.

## Why autonomous resolution is not possible

Category 4 (Authorization): the denied action concerns the sealed evaluator, which is an oversight mechanism of the V3-013 plan (§D step 3), and the harness rule forbids performing on a subagent's behalf an action its permission system denied. The change itself is a correctness fix of my own earlier edit (01f40fa8), is required by an independent Codex finding, and will itself be reviewed by Codex before the seal: (a) both copies of the normaliser (steps_v13.normalise_cache_dir and the evaluator's _normalise_cache_dir) replace only non-null string cache_dir values with the placeholder and keep null as null, at the top level and in every data_used entry; (b) identity_problems() requires a real non-empty path only where the sealed value is the placeholder (data-backed models) and requires null where the sealed value is null (analytic models); (c) evaluator tests with synthetic fixtures for all four cases (bethe null passes, bethe with a path fails, icru90 with a real path passes, icru90 with null or "" fails). No other change. After it lands: regenerate vc-predictions.json on the host (steps_v13.py hash changes), regenerate the frozen manifest, record validation, push, request the Codex review of that exact SHA. Finding 2 of the same review (absolute cache paths) is already fixed and committed at d7b149b6.

## Exact requested input or action

Please choose one: (A) authorize me to apply exactly the edit described in "reason" to the two sealed files (reply "A: apply the evaluator cache_dir fix"), or (B) apply it yourself and reply "done", or (C) grant the implementation worker permission for those two paths and reply "C". No gate waiver is requested; the Codex review of the resulting SHA remains mandatory and no evidence run starts before it.

## Decision or task depending on this request

Not specified.
