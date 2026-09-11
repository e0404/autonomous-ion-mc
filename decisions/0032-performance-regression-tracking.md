# 0032 — Performance-regression tracking (the milestone V6 mechanism)

- Status: accepted
- Date: 2026-09-11
- Task: DEV-027
- Affects: performance engineering (Stage 6), validation strategy (V6)

## Problem

Milestone **V6** is *documented benchmark results with unchanged scientific
validation outcomes before and after each optimisation*. Decisions `0030`/`0031`
established a benchmark harness and two physics-gated benchmarks, each of which
checks cross-backend agreement **within a single run**. What is still missing is the
**temporal** guard: a way to assert that a change (an optimisation, a refactor)
leaves a benchmark's scientific output identical to a **committed baseline**, and to
track how its performance moved. Without this there is no before/after comparison,
so no optimisation can be integrated under the V6 rule.

## Decision

1. **Anchor the physics identity on the reference (float64) result.** A benchmark's
   `digest.reference` is produced by the deterministic float64 reference path, which
   is **hardware-independent**, so it is the stable pin across machines and
   revisions. Wall-clock throughput is machine-dependent and is therefore compared
   **only informationally**. A regression check **gates on the reference digest**
   (the V6 "unchanged scientific outcome") and **reports throughput deltas**.

2. **Committed baselines.** `benchmarks/baselines/<benchmark>.json` holds the
   hardware-independent physics fingerprint (benchmark name, workload config, and the
   reference digest — the gate) plus a machine-tagged throughput snapshot (peak
   histories/s per backend, and the `recorded_on` machine/Warp/version — informative,
   never gated). Baselines are generated from a host run on the reference GPU machine
   (`benchmarks/check_regression.py --emit-baseline`) so the snapshot includes the
   CUDA throughput that the mechanism is meant to track.

3. **The check.** `ionmc.benchmarking.compare_to_baseline(report, baseline)` compares
   a fresh report to a baseline: `physics_ok` is true iff the benchmark name matches,
   the reference `shape_digest` matches, and the reference integral matches to
   `1e-9` (last-bit summation order only). It also returns per-backend throughput
   ratios (current / baseline). `benchmarks/check_regression.py` wraps it as a CLI
   that exits non-zero on a physics regression; `throughput_by_backend` tolerates both
   driver report shapes (explicit `peak_throughput_per_s`, or the `timings`/`scaling`
   rows).

## Validation (`benchmarks/check_regression.py`, `tests/ionmc/test_benchmarking.py`)

- Unit tests (reference path) cover the fingerprint extraction, baseline building,
  and the compare logic: matching physics passes with the throughput ratio computed;
  a changed reference digest, a changed integral, or a benchmark-name mismatch each
  fail `physics_ok`; a backend missing from a run yields a `None` ratio, not an error.
- On the host runner, each benchmark is re-run and checked against its committed
  baseline; the reference digest matches (the physics is unchanged at this SHA),
  confirming the mechanism end to end.

## Consequences

- Every future optimisation can now be integrated under an explicit before/after
  guard: re-run the benchmark, `check_regression` against the committed baseline, and
  the change is admissible only if the reference physics is unchanged — with the
  throughput delta documented. This unblocks the Stage 6 optimisation work.
- **Deferred:** a persisted *time series* of results (trend history beyond a single
  baseline), automatic baseline refresh in CI, per-backend throughput *regression*
  thresholds (as opposed to reporting), and wiring the check into CI (it needs a GPU
  and the datasets, so it stays a host-runner step for now).
