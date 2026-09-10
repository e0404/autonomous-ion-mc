# Development

This section documents development tooling and contributor-facing workflows,
and hosts the [development roadmap](roadmap.md) (capability stages,
dependencies, validation milestones and the currently selected milestone).

## Controlled task lifecycle

`AGENTS.md` is the authoritative description of the task lifecycle. In short,
every repository change is made on a `task/<id>-<slug>` branch in an isolated
worktree under `~/aiprojects/ion-mc-worktrees/<id>`, committed with the fixed
autonomous identity, validated for its exact SHA (the merge gate refuses SHAs
without a passing local validation record), pushed, pull-requested, checked by
GitHub CI, squash-merged into `develop`, and the worktree retired with the
branch preserved.

Real CPU/GPU execution of a committed task state uses `run_host_validation`,
which runs an argv inside a second sandbox with the worktree at `/workspace`.
The host-runner Python environment currently provides `warp-lang` and `numpy`
only, without network access; scripts under `validation/` therefore depend on
nothing else.

## Test and quality commands

GitHub CI runs the following; run the same locally before integration:

```bash
uvx --from pre-commit==4.6.2 pre-commit run --all-files
uv run --with pytest --no-project python -m pytest tests -q
python3 infrastructure/docs/generate_docs.py && uvx zensical build --clean --strict
```

## Known environment caveats

- **Sandbox filesystem view is fixed at session start.** Directories created
  on the host during a session (for example the first task worktree root)
  are not visible to the orchestrator sandbox until the next session. The
  temporary workaround used in the kickoff session, its integrity argument
  and its stopping condition are recorded in decision `0003`.
- **Sandbox network egress can be unavailable.** When it is, tests and
  quality tools are executed on the host through the Codex worker and by
  GitHub CI; the validation record for the affected task states this.
