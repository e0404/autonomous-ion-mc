# 0037 — V2 kickoff preflight findings and runtime-environment constraints

- Status: accepted
- Date: 2026-09-30
- Task: V2-001
- Affects: reproducibility, dependency footprint, validation strategy

## Problem

The isolated v2 run (`experiment/prompts/kickoff-v2-isolated.md`) starts with an
infrastructure preflight. The findings determine which execution paths, data
sources and dependencies the scientific implementation may rely on.

## Findings

1. **Condition and manifest.** `.ionmc-condition.json` identifies
   `experiment-v2`, `scientific_implementation_inherited: false`, integration
   branch `v2/develop`. All 71 overlay hashes match the freeze commit
   `153e2f6`; 18 of them differ from the current tree only because the later
   infrastructure commits (#60–#63) modified those files, which Git history
   records. The excluded v1 commit objects are absent from the clone.
2. **History isolation.** All checks pass except `only_v2_refs`. The single
   offender is `refs/codex/turn-diffs/checkpoints/...`, a Codex CLI turn-diff
   checkpoint that points at a *tree* object identical to the HEAD tree (108
   files, empty diff), not a commit. It is a local tool artifact, not v1
   history. An amendment tolerating tree-only Codex checkpoint refs was
   prepared but denied by the session's permission policy, so the check is
   left unchanged and the finding is recorded here for the operator. No ref
   was deleted.
3. **Sandbox view versus service view.** From the agent shell, the preflight
   also reports an unclean checkout (sandbox-mounted device-node files such as
   `.bashrc`), no `/dev/dxg`, unreadable engine registry and denied `gh`.
   These are sandbox mounts and permissions, not repository state: the
   controlled MCP services created the V2-001 worktree from a clean checkout
   and report all three reference engines registered (OpenTOPAS 4.3 / Geant4
   11.4.2, MCsquare public static binary at source revision `211eefe6`, FRED
   3.76.0). GPU execution is only available through `run_host_validation`.
4. **Host-runner runtime.** The controlled GPU sandbox uses a fixed virtual
   environment: CPython 3.12.0, `numpy` 2.5.3, `warp-lang` 1.17.0, `pytest`
   9.1.1 — nothing else, and no network. Warp reports CUDA Toolkit 12.9 and
   compiles CPU kernels in the agent sandbox.
5. **Warp capabilities verified.** `@wp.func` functions run from Python scope
   (Python floats give float32-rounded transcendental builtins; `wp.float64`
   scalars give full float64 arithmetic and builtins); generic kernels with
   `typing.Any` specialize for float32 and float64 arrays at launch;
   `wp.rand_init`/`wp.randf`/`wp.randn` are kernel-only (not callable from
   Python scope); structs, 3-D arrays, atomics and while loops work.
6. **Infrastructure tests.** 337 passed, 2 skipped in the sandbox.

## Decision

- Runtime dependencies of `ionmc` are limited to `numpy` and `warp-lang`, so
  the package runs unchanged inside the fixed host-runner environment and in a
  clean user installation. Optional extras may add plotting or IO libraries
  later, but no required scientific path may depend on them; human-readable
  plots are produced without third-party plotting libraries.
- Tests are discovered from `src/` through `pytest`'s `pythonpath` setting so
  the same command works in CI, in the sandbox and inside the host runner
  without installing the package.
- GitHub CI installs the locked environment (`uv sync --locked`) and runs only
  tests not marked `local`/`cuda`; every scientific and GPU validation runs
  locally under the exact-SHA gate.
- The stale Codex checkpoint ref is reported, not removed; the Codex worker
  may create further tool refs during the run, which are equally not history.

## Validation

`uv sync --extra dev && uv run pytest`, `pre-commit run --all-files`, the
strict documentation build, and a host-runner execution of the package tests
at the committed V2-001 SHA (recorded through the local validation gate).
