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
uvx --from pre-commit==4.6.2 pre-commit run --all-files   # hygiene, ruff, ruff-format, mypy
uv sync --extra dev && uv run python -m pytest tests -q
uv run mypy
python3 infrastructure/docs/generate_docs.py && uvx zensical build --clean --strict
```

Ruff (lint and format) and mypy apply to `src/` and `tests/ionmc/`; the
legacy scripts under `infrastructure/` and their tests under
`tests/infrastructure/` (about 120 pre-existing style violations, mostly
line length) are excluded until a separate cleanup task.
Warp tests are skipped automatically where `warp` or a CUDA device is
unavailable (markers `warp`, `cuda`).

## Environment history

- **Kickoff session (DEV-001, first DEV-002 commits).** The orchestrator
  sandbox could not see the task worktree directory (a user-level
  Claude Code setting blocked reads outside the primary working directory)
  and had no network egress. Development happened in a scratch clone whose
  tree was transferred into the worktree through a staging ref, and tests
  ran on the host through the Codex worker; decision `0003` records the
  workaround, its integrity argument and its verification.
- **From the DEV-002 resumption on 2026-09-10.** Direct worktree access and
  sandbox network egress were verified and decision `0003` was marked
  superseded. Tests, pre-commit, mypy and the documentation build run
  inside the sandbox in the task worktree; only GPU execution still goes
  through `run_host_validation`.
