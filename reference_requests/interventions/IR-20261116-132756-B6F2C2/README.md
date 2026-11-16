# Human Intervention Request IR-20261116-132756-B6F2C2

**Task:** V3-INFRA-006
**Category:** authorization
**Status:** open

## Summary

One external action is needed: restart the codex-worker MCP service so that the merged V3-INFRA-006 change (signed, timestamped review records; PR #90 merged into v3/develop at 7dfab634) takes effect. After the restart the service creates its signing key on the first finished review; the orchestrator then runs one discovery review and publishes the public key through the documented V3-INFRA-003/005 procedure. This supersedes the earlier request IR-20261116-123921-312B59 (whose option (A) was withdrawn by its addendum 1): no decision beyond the restart is requested.

## Why autonomous resolution is not possible

Category 4 (Authorization: deployment to infrastructure outside the experimental workspace). The running codex-worker MCP server loads infrastructure/experiment_v3/review.py at start-up; the merged change (signing of every finished review record with a service key, signer "codex-worker"; verifier infrastructure/validation/review_records.py; key-publication evidence rule) is inactive until the server process is restarted, which only the operator can do (precedent: IR-20261107-235713-1724FF for V3-INFRA-002/003 and IR-20261109-112802-D6BAE3 for V3-INFRA-004). Why it is needed: the independent Codex review of the V3-008 amendment I3-FU-1 (nine rounds) converged on the finding that the evaluator's review-chronology pin is not authenticated by the review service; the merged change provides the service-signed record the evaluator will consume (together with the sealed launcher already on the V3-008 branch, which answers the second remaining finding). Gates met for the merge: Codex REVIEW-4a570742 (pass, one minor finding fixed in fe3d6744) and REVIEW-64ad632d (pass, no findings), local-validation records VAL-20261116-125733 and VAL-20261116-131416 (615 passed, 3 skipped), CI run 43875086986 green (docs, pre-commit, tests), squash-merged as PR #90.

## Exact requested input or action

Please restart the codex-worker MCP service (the same restart procedure as for V3-INFRA-002/003/004) and reply "restarted" (optionally with the time). Nothing else is requested; the key publication and the follow-up reviews are autonomous.

## Decision or task depending on this request

activation of signed codex-worker review records; publication of the codex-worker public key (trust-file branch); the V3-008 amendment I3-FU-1 review pins (evaluator consumes the service-signed record), hence the start of the I3 follow-up diagnostics (step 1, smokes, 58 TOPAS evidence runs)
