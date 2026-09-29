# V2 history isolation — condition 2.1

The operator approved a fresh implementation run in an independent local clone
of the existing GitHub repository. This amendment supersedes earlier instructions
to reconstruct or inspect v1 project history. Scientific requirements, autonomy,
local validation, controlled integration and release gates are unchanged.

## Permitted experimental inputs

Use the supplied v2 working tree, the pre-scientific infrastructure ancestry,
v2/develop, v2/main and task branches/PRs created for this v2 run. The supplied
requirements, preparation records and audit lessons are intentional inputs.
Commit identifiers and audit hashes in provenance files are opaque identifiers,
not permission to fetch their content. V2 remains informed by the v1 audit; this
condition isolates the prior implementation and its development/results history.

Do not fetch, inspect, reconstruct or copy v1 scientific branches, tags, commits,
PRs, releases, source, tests, results, logs, transcripts or audit archives. This
applies equally to local files, Git commands, MCP, APIs and browser/web research.
Do not broaden fetch refspecs, add remotes or object alternates, enable automatic
tag fetching, or use another checkout to circumvent the condition. Do not inspect
the same repository's legacy documentation site or published artifacts for v1
implementation details. Independent public scientific research remains allowed.
If excluded material is encountered accidentally, stop using it and record the
exposure and its likely effect without copying the material into task context.

## Clone and runtime arrangement

The integration checkout is `/home/wahln/aiprojects/ion-mc`, with its own `.git`
directory and object database. Its complete ancestry starts from the pre-scientific
baseline plus v2 infrastructure commits; it does not contain the completed v1
scientific ancestry. Keep full v2 ancestry for merge-base and release checks.
Task worktrees share only this independent v2 database.

The original checkout is archived at `/home/wahln/aiprojects/ion-mc-v1`; its
leftover task files are at `ion-mc-v1-worktrees`. The previous shared
v2 worktree is preserved at `/home/wahln/aiprojects/ion-mc-v2-before-isolation`.
Those paths, the v1 task-worktree directory and the independent audit directory
and archive are excluded in v2 Claude read permissions and sandbox settings.
Start a new Claude session in the v2 clone; do not resume a v1 or setup session.
Do not attach the excluded directories to its workspace.

At the operator's request, the independent v2 clone was moved from `ion-mc-v2`
back to `ion-mc`. It retains complete v2-only ancestry, not a shallow history.
Active tasks now use `ion-mc-worktrees`; the archived Git worktree links were
repaired without deleting or rewriting history. The launcher sets
`CLAUDE_CODE_PROJECT_DIR_NAME=experiment-v2`, and project settings explicitly use
`~/.local/share/ionmc-experiment/experiment-v2/claude-memory` for auto memory.
This avoids loading v1 memory when reusing the old checkout pathname; the shared
Claude authentication/configuration profile is preserved. Do not resume v1
conversations. See the [Claude memory documentation](https://code.claude.com/docs/en/memory).

Claude project settings select separate v2 telemetry, Codex raw-output and local
validation directories using the existing IONMC_TELEMETRY_DIR,
IONMC_CODEX_RAW_DIR and IONMC_VALIDATION_DIR variables. Their names are unchanged;
the readers now expand `~`. The general telemetry directory is no longer an
allowed sandbox read/write root. IONMC_NTFY_URL and the registered native engines
are unchanged. Old logs remain in the operator's archive.

The orchestrator chooses unique `V2-` task IDs (for example `V2-001`); the task
MCP has no automatic numbering counter. V2 creation rejects IDs without that
prefix and checks local refs plus an exact-ID remote namespace query before
creating a worktree. An existing remote ID is rejected even if its description
differs. This reads only ref metadata, not branch contents. Retain all old remote
task branches; GitHub PR numbering remains repository-wide and creates no Git
branch collision. The availability check is not an atomic distributed allocator;
concurrent orchestrators must still coordinate ID allocation.

The controlled Git services already fetch the configured v2 branches explicitly.
New task branches may be pushed and preserved normally. CI inspection and failure
logs verify the task SHA with an exact-ref `git ls-remote` query; they do not
require a local origin/task tracking ref or fetch excluded history. The clone fetch config
contains exactly:

```text
+refs/heads/v2/develop:refs/remotes/origin/v2/develop
+refs/heads/v2/main:refs/remotes/origin/v2/main
```

`remote.origin.tagOpt` is `--no-tags`. This stops automatic tag import; creating
new v2 release tags remains permitted. No v1 remote refs or tags are installed.

## Operator reproduction

Use an unused destination. These commands do not replace an existing checkout:

```bash
git clone --single-branch --branch v2/develop --no-tags \
  https://github.com/e0404/autonomous-ion-mc.git /path/to/new-v2
cd /path/to/new-v2
git config --local --add remote.origin.fetch \
  '+refs/heads/v2/main:refs/remotes/origin/v2/main'
git fetch --no-tags origin
git branch --track v2/main origin/v2/main
python3 -m infrastructure.experiment_v2.history
```

Use a network clone or `--no-local` for a local source. Do not use `--shared`,
`--reference`, a filesystem copy, or a worktree from v1. On another workstation,
update the declared absolute sandbox/task paths and provision engines normally.
Do not rewrite original prompts or move historical start/ready tags.

The history check verifies independent Git storage, no object alternates, full
v2 ancestry, restricted fetch settings, no automatic tags, expected branch
availability, permitted ref names and absence of known excluded commit objects.
Run it against the integration clone, not an individual task worktree. The main
preflight includes it as a required unattended-launch check. It does not enumerate
or prove absence of every possible v1-derived file, commit or external source.

This is a practical experimental boundary, not an inaccessible remote archive.
V1 remains on the same remote and public-network access remains available. Prompt
policy and scoped integration tools constrain use; they are not a network firewall.
Project read rules are documented in the
[Claude permission reference](https://code.claude.com/docs/en/permissions).
An operator with broader permissions can still access the preserved v1 record.

After the ordinary smoke and notification preflight passes, launch with:

```bash
claude "Read experiment/prompts/kickoff-v2-isolated.md and carry out that autonomous run."
```
