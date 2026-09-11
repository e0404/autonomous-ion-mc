# 0033 — Release validation suite and the first-release criteria

- Status: accepted
- Date: 2026-09-11
- Task: DEV-028
- Affects: release readiness (Stage 6), the `develop` → `main` release workflow

## Problem

`EXPERIMENT.md` (*Releases*) requires that a release from `develop` to `main` has
**passed the defined release validation suite**, with documented validation status,
reproducible benchmark information, a version tag, and identification of the
experiment configuration. The physics milestones **V0–V5** are complete and the
Stage-6 benchmark + regression infrastructure (**V6**) is in place, but the *release
validation suite itself was never defined* — so there was no single, reproducible,
machine-readable gate that says "this state is release-ready." This decision defines
and implements it.

## Decision

1. **The release validation suite is the union of every milestone gate.**
   `validation/release_validation.py` runs, on the controlled host runner, **every**
   milestone validation script `validation/v*.py` (V0–V5) plus the two benchmark
   physics gates and their performance-regression checks against the committed
   baselines (V6). It aggregates a single machine-readable report and exits 0 **iff
   every** suite passes. The reference-Python correctness path, the Warp CPU/CUDA
   parity, the stochastic-consistency gates and the deterministic cross-backend
   digests are therefore all re-exercised at the exact release SHA.

2. **Composition is discovered, not hand-maintained.** The orchestrator globs
   `validation/v*.py` — a pattern that by construction matches every milestone gate
   while excluding both the `warp_cuda_smoke` diagnostic (not a milestone gate,
   non-standard interface) and this orchestrator itself — and runs each as a
   subprocess with only the flags it supports (`--require-cuda`, `--cache-dir`), so new
   milestone validations are included automatically. Each milestone CLI follows the shared exit
   convention (0 pass, 3 gate failure, 4 dataset missing); a 3 **or** a 4 both count
   as not-release-ready (a release requires the datasets present and every gate
   green).

3. **Release criteria.** A `develop` state is release-ready when: the release SHA has
   passing required CI (unit tests, lint, type-check, docs — the suite is physics-only
   and does not re-run these); the release validation suite passes on the host runner
   (CPU + CUDA); the result is recorded for
   the exact SHA (`record_local_validation`); the documented validation status and the
   reproducible benchmark baselines are committed; the version is set; and the
   experiment configuration is identified (the report stamps `ionmc` / Python / NumPy
   / Warp versions and `--require-cuda`). Cutting the tag from `develop` to `main` is
   the operator-or-autonomous release step once these hold (`EXPERIMENT.md` delegates
   the release decision to the autonomous system unless reserved).

## Validation (`validation/release_validation.py`, `tests/ionmc/test_release_validation.py`)

- On the host runner (CPU + CUDA) the suite runs end to end and its
  `release_ready` flag is the gate; the aggregated per-suite report is the documented
  validation status.
- Reference-path tests guard the orchestrator's wiring — milestone-validation
  enumeration (diagnostics/itself excluded), per-script flag detection, and the
  presence of the benchmark baselines it regression-checks — in CI without a GPU.

## Consequences

- There is now a single reproducible command that certifies a `develop` state as
  release-ready, satisfying the `EXPERIMENT.md` release requirement and making the
  first tagged release a well-defined, evidence-backed step rather than a judgement
  call.
- **Deferred:** wiring the suite into CI (it needs a GPU and the datasets, so it stays
  a host-runner step); a packaging/wheel-build gate; and the actual tag-and-merge to
  `main`, which follows once the suite passes at a release SHA and the version is
  finalised.
