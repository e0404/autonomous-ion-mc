# 0003 — Out-of-worktree development via a staging ref (kickoff workaround)

- Status: accepted (temporary; see *Stopping condition*)
- Date: 2026-09-10
- Task: DEV-001
- Affects: reproducibility, provenance of task history, maintainability of the
  task lifecycle

## Problem

`AGENTS.md` requires development work to happen in the isolated task worktree
created by `create_task_worktree` under `~/aiprojects/ion-mc-worktrees/`,
committed there by `commit_task_changes`. During the kickoff session the
orchestrator's sandbox could not see that directory at all: the sandbox's
view of the home directory is assembled once at session start from paths
that exist at that moment, and `~/aiprojects/ion-mc-worktrees` was created
for the first time by `create_task_worktree` during the session. Neither the
shell nor the file tools could read or write the worktree, although every
controlled tool (task manager, commit, host runner, Codex worker,
integration) operated on it normally from outside the sandbox.

Independently, sandbox network egress was unavailable for the whole session
(every proxied connection was aborted, including to the pre-authorised PyPI
and GitHub hosts), so `pytest`, `pre-commit` and `zensical` could not be
installed inside the sandbox.

## Context

- The controlled lifecycle's integrity properties come from: the task branch
  being created from synchronized `develop`; the commit being made with the
  fixed autonomous identity by `commit_task_changes`; the validation record
  and host runs being keyed to the exact committed SHA; and the PR/merge tools
  acting only on that branch.
- The main checkout and the task worktree share one Git object store, and the
  main checkout *is* visible to the sandbox.
- The Codex worker executes on the host inside the worktree with `uv` and
  network access.

## Candidate approaches

1. **Request human intervention to restart the session.** Rejected: none of
   the permitted intervention categories in `EXPERIMENT.md` applies to a
   tooling inconvenience, and the problem is solvable autonomously.
2. **Let the Codex worker do all implementation from natural-language
   instructions.** Rejected as the primary mode: it would remove the
   orchestrator's direct authorship and review of scientific code and make
   Codex, rather than the orchestrator, the implementer of every task.
3. **Develop in a local clone of the task branch inside the session
   scratchpad, publish the result into the shared object store as a
   non-branch ref, and have the Codex worker materialise exactly that tree in
   the worktree (selected).**

## Selected approach

1. `git clone --branch task/<id>-<slug> /home/wahln/aiprojects/ion-mc <scratchpad>/<id>`.
2. Implement, test and review in the clone with the full tool set.
3. Snapshot: commit in the clone (staging commit, any identity), then
   `git -C /home/wahln/aiprojects/ion-mc fetch <clone> +HEAD:refs/staging/<id>`.
   This writes objects into the shared store and updates a ref that is not a
   branch and not checked out anywhere, so it is never refused and never
   touches `develop`, `main` or the task branch.
4. Codex worker, in the task worktree, runs exactly
   `git read-tree -u --reset refs/staging/<id>` (index and working tree now
   equal the staged tree; HEAD unchanged) and reports `git status`.
5. `commit_task_changes` commits with the fixed autonomous identity.
6. Integrity check from the sandbox:
   `git rev-parse task/<id>-<slug>^{tree}` must equal
   `git rev-parse refs/staging/<id>^{tree}`. Only then are validation records
   created, host runs launched, and the branch pushed.
7. Tests and quality checks that need packages unavailable in the sandbox are
   executed on the host through the Codex worker (`uv run --with pytest
   --no-project python -m pytest tests -q`, `uvx pre-commit run --all-files`,
   `uvx zensical build --clean --strict`) and independently by GitHub CI.

## Rationale

All integrity properties of the lifecycle are preserved: the task branch
still descends from synchronized `develop`, the commit is still made by the
controlled commit tool with the fixed identity, and the tree-hash equality
check proves that what was committed is byte-for-byte the tree the
orchestrator authored, tested and reviewed. The Codex worker's role is
reduced to one deterministic Git command whose effect is verified
independently.

## Expected tradeoffs

- Two extra steps per commit and one Codex invocation per commit.
- The staging commits in the clone are not part of the recorded history; the
  provenance of a task commit is the controlled commit plus this decision.
- Any unit test or quality tool that could not run in the sandbox is
  evidenced by the Codex worker's report and by CI rather than by the
  orchestrator's own execution; the validation record must say so.

## Stopping condition

This workaround is to be used only while the task worktree directory is
invisible to the orchestrator sandbox. Because `~/aiprojects/ion-mc-worktrees`
now exists on the host, later sessions are expected — but not yet verified —
to see task worktrees directly; the first action of the next session is to
check, and if direct access works this decision is marked superseded and the
normal lifecycle is used. The network workaround (host-side test execution
via the Codex worker) is likewise dropped as soon as sandbox egress works.

## Validation strategy

Per use: the tree-hash equality check in step 6, and the Codex worker's
reported `git status` showing only the intended paths. For `DEV-001` these
checks are recorded in the task's local validation record.

## Later validation outcome

To be filled in when the workaround is retired or when direct worktree access
is confirmed.
