# Human Intervention Request IR-20261003-220532-42BA5F

**Task:** V3-002R
**Category:** authorization
**Status:** resolved (2026-10-04; operator executed the remote deletion; verified absent by exact-ref query at 2026-10-03T22:12:54Z UTC; see `response.json` and decision 0038)

## Summary

Execute the authorized remote deletion of the V3-002 task branch (option B of IR-20261003-192301-6CF85A); local branch retained

## Why autonomous resolution is not possible

The operator resolved IR-20261003-192301-6CF85A with option B: the remote task branch task/v3-002-external-data-layer-materials-and-stopping-power-tables (pushed commits 7a2d817..6eb6369 contain transcribed NIST STAR rows) must not stay public, while the local branch is retained. EXPERIMENT.md pre-authorizes only pushing task branches, pushing develop through the integration workflow, pull requests, CI and release tags; deleting a remote branch is not pre-authorized, the protocol says task branches should not be deleted during the active experiment, and the autonomous system has no controlled tool for remote deletions. Execution therefore needs the operator (or an explicit one-line delegation). The autonomous system has NOT deleted anything; decision 0038 records authorization only until the remote branch is verified absent.

## Exact requested input or action

Please either (1) run, from a clone with push rights: `git push origin --delete task/v3-002-external-data-layer-materials-and-stopping-power-tables` (the local branch in /home/wahln/aiprojects/ion-mc is untouched and the squash-merged develop content 8a9ad488 is NIST-free), then reply "deleted"; or (2) reply "autonomous system may run the deletion" to delegate exactly that single command to the autonomous system. After either, the orchestrator verifies with the exact-ref query `git ls-remote --exit-code origin refs/heads/task/v3-002-external-data-layer-materials-and-stopping-power-tables` (exit code 2 = absent; never an unrestricted `--heads` enumeration, per experiment/v3/HISTORY-ISOLATION.md) and records the completion date in decision 0038.

## Decision or task depending on this request

IR-20261003-192301-6CF85A
