# V3 history isolation

The active independent clone is /home/wahln/aiprojects/ion-mc. It starts at
v2 infrastructure commit dd33266073f84b630b3da390c9f762095c4a0870, before all v2
scientific work. Complete infrastructure ancestry is retained for integration;
this is not a shallow clone. No v1/v2 scientific implementation is inherited.

The operator preserved v1 at ion-mc-v1 and v2 at ion-mc-v2, with respective
*-worktrees archives, old caches/state and session histories. V2 ended incomplete
at 8a3d30ae337554e50b8712b5df3e1532d616dd89; unfinished V2-006 work and its dirty
files are preserved. Its closure record is held in the operator's v2 state.
These opaque identifiers are provenance, not permission to retrieve content.

The autonomous runtime must not fetch, inspect, reconstruct or copy earlier
scientific branches, commits, tags, PRs, documentation, outputs, logs, conversations,
memory or local archives via Git, MCP, web or filesystem. Only supplied historical
infrastructure/prompt records and scientific requirements are intentional inputs.
Do not broaden refspecs, add remotes or object alternates, or fetch old tags.
If excluded content is encountered, stop using it and record the exposure.

Fetch only v3/develop and v3/main, with remote.origin.tagOpt=--no-tags. Preserve
new V3 task branches and create only v3 release tags. The task manager queries
exact task IDs to prevent collisions without fetching other task histories.
The independent .git database must not contain excluded scientific commits.
The history guard is a practical boundary; the shared public remote is not a
network firewall. Canonical historical prompts remain unmodified.

Launch a fresh experiment-v3 session through ~/bin/claude-ionmc. Do not resume a
v1/v2 conversation or attach archived folders. Cache/state/memory/telemetry use
experiment-v3 directories; no scientific caches or outputs are copied from v2.

Codex may create local UI checkpoint refs under refs/codex/turn-diffs/. Preflight
accepts these bookkeeping refs only when their object is an exact tree already
committed in the local v3/develop ancestry. Commit, blob, unrelated tree and other
ref namespaces remain subject to the existing restrictions. This exception neither
fetches history nor changes excluded-object, remote or refspec checks.
