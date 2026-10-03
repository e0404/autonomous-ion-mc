# Human Intervention Request IR-20261003-192301-6CF85A

**Task:** V3-002
**Category:** authorization
**Status:** resolved (2026-10-04; operator response in `response.json`; option B chosen for the task branch, develop history retained; cleanup not yet performed — see decision 0038 and the follow-up request)

## Summary

Authorization to rewrite or delete pushed task-branch commits of V3-002 that contain NIST SRD 124 excerpts (documented repository error)

## Why autonomous resolution is not possible

Commits 7a2d817 through 6eb6369 on the pushed, preserved task branch task/v3-002-external-data-layer-materials-and-stopping-power-tables contain about a dozen transcribed NIST PSTAR/ASTAR rows as test fixtures and a comparison JSON with per-energy model/reference ratios from which NIST table values can be recovered. NIST's statement (https://www.nist.gov/open/license) secures copyright on SRD compilations with all rights reserved and grants no redistribution licence. The material was removed from the tree at d3945d1 before integration, so the squash merge into v3/develop and any release will not contain it; the independent Codex reviewer nevertheless flags the reachable task-branch history. EXPERIMENT.md forbids force pushes to shared task branches except for recovery from an explicitly documented repository error; decision 0038 now documents this error. Rewriting or deleting pushed history on the designated GitHub repository is an irreversible, externally consequential action that the autonomous system may not perform without explicit authorization (EXPERIMENT.md, permitted intervention 4), and the controlled integration tools provide no force-push capability.

## Exact requested input or action

Please choose one: (A) authorize the operator-performed rewrite of the task branch so that commits 7a2d817..6eb6369 are replaced by NIST-free equivalents (or the branch is reset to start at d3945d1), reporting the new branch head; (B) authorize deletion of the remote task branch after the squash merge, keeping only the local copy; or (C) decide that the private repository may retain these historical commits as an experimental record. The autonomous system proceeds with integration of the NIST-free head (the squash merge excludes the material) unless instructed otherwise; the chosen remedy will be recorded in decision 0038.

## Decision or task depending on this request

Final remedy recorded in decision 0038; whether the preserved task-branch history of V3-002 is rewritten, deleted or retained.
