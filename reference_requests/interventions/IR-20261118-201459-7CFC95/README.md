# Human Intervention Request IR-20261118-201459-7CFC95

**Task:** V3-013
**Category:** authorization
**Status:** resolved 2026-11-19 (operator: "A: prepare the PR"; see response.json)

## Summary

The frozen V3-013 ionmc samples (1e6 He histories per case and model, warp-cpu float64, one run each) cannot complete inside the controlled host runner's hard cap of 7200 s: three timing smokes at a non-evidence seed (k=99) measured 4.5 ms/history at 75 MeV/u and 8.6 ms/history at 200 MeV/u (10 000 → 93 s, 50 000 → 437 s), i.e. ≈ 75 min, ≈ 110 min and ≈ 143 min for the three cases. The seal (REVIEW-0e6ca486 at ad73d94e, pin fd30ee82 reviewed by REVIEW-7db57441) is in place and no evidence run has started. I need one authorization: (A, recommended) raise MAX_TIMEOUT_SECONDS in infrastructure/host_runner/host_runner.py (line 34, currently 7200) to 14400, which keeps the frozen plan, the seal and the single-process rule unchanged; or (B) authorize me to amend the plan and the sealed evaluator so that each 1e6 sample is produced as four chunked host runs (--scale 0.25, k = 0..3, 40 batches each, 160 paired batches in total), which requires a new seal and its two reviews.

## Why autonomous resolution is not possible

Category 4 (authorization): the runner cap is an execution limit set in the experiment's infrastructure that I must preserve ("Preserve all … execution limits"), so I will not change it on my own; and option B edits a sealed validation script, which your earlier authorization (IR-20261118-182151-FD0238) covered only for the specific cache_dir fix. Option A details: the MCP server invokes the runner script from the main clone at every call (infrastructure/mcp/host_runner_mcp.py uses REPO_ROOT), so a merged one-line change takes effect without a restart; I can prepare it as task V3-INFRA-008 through the controlled workflow (commit, Codex review, validation, PR, merge) if you say so, or you can apply it directly. Each sample run stays a single process with two threads; the six samples total ≈ 10 h of compute, run strictly one at a time. Option B is scientifically equivalent (same seed family ion13, same histories per case, common random numbers kept per chunk, delete-one-batch jackknife over 160 instead of 40 paired batches) but changes the frozen sample definition after the seal and costs ≈ 1–2 h of re-seal overhead; if chosen, the amendment would be prospective and dated, and the k=0..3 seeds have not been seen (only k=99 was used for the smokes). Measurements: RUN-20261118T200117Z-173deba7 (case 1, icru90, 10 000 histories, 52.0 s wall), RUN-20261118T200243Z-7502b1a2 (case 3, bethe, 10 000, 93.2 s), RUN-20261118T200438Z-5cb9e01e (case 3, bethe, 50 000, 437 s); all exit 0, pass, counters zero, identity blocks equal the sealed pins, environment digests equal the seal, Warp kernel attestations identical. Not a clinical claim.

## Exact requested input or action

Reply "A: prepare the PR" (I create V3-INFRA-008 raising MAX_TIMEOUT_SECONDS to 14400 through the controlled workflow and merge it after review/CI), or "A: applied" (you changed it in the main clone on v3/develop), or "B: amend" (authorizing the described edit of validation/plans/v3-013-acceptance.md and validation/scripts/stopping/v13_evaluate.py plus tests, followed by a new seal and reviews). No scientific tolerance, acceptance criterion or review gate changes are requested under either option.

## Decision or task depending on this request

Not specified.
