# Addendum 1 (2026-11-18, orchestrator)

Status while the request is open: V3-013 implementation continued without the review gate. Commits 8a3791bc
(S1 D-3 digest), dd8929dd (S2a), afe42b30 (S2b), 983eb1dd (S3), f524fe90 (S4), 5248730e (S5a), 4130d569 (S5b),
039d4aad, 6abda9eb and f1dfa55d (seal candidate: frozen manifest, host-produced predictions and environment
record) are pushed on PR #92 with exact-SHA validation records and green CI through 5248730e. Every attempt of
`start_codex_review` since the resume (about ten, the last at 12:10 UTC today) returns "Error executing tool
start_codex_review"; `inspect_codex_review` keeps working.

Consequence: the plan's seal (`validation/plans/v3-013-acceptance.md` §D step 3) is the signed Codex review of the
seal candidate's exact SHA; no ionmc sample, TOPAS smoke or evidence run may start before it. The only work left
that does not need the gate is the S6 documentation step. No waiver is requested: the runs wait for the service.
