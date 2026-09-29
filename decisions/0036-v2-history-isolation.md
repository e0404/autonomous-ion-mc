# 0036 — Isolate v2 from v1 implementation history on the existing remote

Date: 2026-09-29. Status: accepted by the operator before scientific kickoff.

## Problem and decision

The original v2 working tree excluded scientific implementation files, but its
Git worktree shared all v1 objects and refs. The operator requested practical
history separation while retaining the same remote and controlled workflow.

Replace that working tree with an independent full clone fetching only v2/develop
and v2/main, with no automatic tag import. The branch already starts at the
pre-scientific baseline; the v1-based setup infrastructure was copied into new v2
commits rather than merged. A shallow clone is unnecessary and would complicate
ancestry-based integration and release checks. A separate remote would provide
stronger separation but is outside the selected experimental condition.

Preserve the former worktree and v1 refs in the operator workspace. Explicitly
exclude v1 local checkout/audit paths in the v2 Claude settings. Separate runtime
logs, retain prior canonical prompts, and introduce kickoff-v2-isolated.md plus
the versioned HISTORY-ISOLATION.md amendment. Provenance hashes remain identifiers
only. Neither scientific requirements nor validation gates are weakened.

## Validation and limitations

Exercise real local Git repositories containing a divergent v1 scientific branch:
an independent narrow clone must omit its commit, while broad fetch settings,
shared worktrees, alternates, extra remotes, v1 refs and an explicit excluded-SHA
fetch must fail the check. Task creation, push and release-branch fetches must
continue to work without restoring v1 history. Run the infrastructure tests,
pre-commit and strict documentation build; record exact-SHA results through the
existing local validation gate before integration.

The remote still holds v1 and web access remains enabled. This is a policy-based
experiment boundary with checks for accidental configuration mistakes, not an
adversarial security boundary. The v2 requirements intentionally retain lessons
from the audit. Start with a new orchestrator session to avoid prior conversation
context. Installed-Claude enforcement and real notification delivery remain
separate runtime checks; unit tests do not establish those claims.
