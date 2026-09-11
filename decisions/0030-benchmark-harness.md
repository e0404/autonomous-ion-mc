# 0030 — Reproducible benchmark harness and the first performance benchmark

- Status: accepted
- Date: 2026-09-11
- Task: DEV-025
- Affects: performance engineering (Stage 6), validation strategy (V6)

## Problem

Stage 6 requires *high-throughput execution*, *reproducible benchmarking*, *efficient
batching* and *performance-regression testing*, culminating in a first tagged release.
Validation milestone **V6** is *documented benchmark results with unchanged scientific
validation outcomes before and after each optimisation*. Before any kernel or
memory-layout optimisation can be justified, the project needs a **reproducible way to
measure performance** and to prove that a change left the physics untouched. No such
machinery exists yet (`benchmarks/` is empty). This decision establishes the harness
and the first benchmark; it does **not** optimise anything.

## Decision

1. **Harness as a library, drivers as scripts.** The reusable machinery lives in
   `ionmc.benchmarking` (importable and unit-tested like any source module); runnable
   benchmark *drivers* live under `benchmarks/` and use it, mirroring how
   `validation/` scripts use source modules. This keeps benchmarks first-class,
   testable code rather than ad-hoc scripts.

2. **Benchmarks are separate from validation gates.** A benchmark records wall-clock
   time, which is **hardware-dependent and is never asserted against an absolute
   threshold**. What a benchmark *does* gate is the physics: every driver pins a
   reproducible **scientific digest** of its output and checks cross-backend agreement
   on results that carry their **absolute magnitude** (per-history dose, never
   unit-normalised — which would hide a uniform-scale error), so a change that
   silently breaks the physics fails the benchmark. This is the concrete mechanism for
   the V6 "unchanged scientific outcome" requirement.

3. **Full provenance on every run.** `capture_provenance()` records `ionmc` version,
   Python/NumPy versions, platform/machine/processor, CPU count, git SHA, and the Warp
   version and available devices — so a recorded result is interpretable and
   comparable across machines and revisions.

4. **Honest timing.** `measure()` runs untimed warm-up iterations (excluding one-time
   Warp kernel compilation) then times `repeats` iterations with `time.perf_counter`,
   reporting min/median/mean/stdev and a throughput (`work_units / median`, e.g.
   histories/s). Because Warp CUDA launches are asynchronous, a device-synchronisation
   callback (`make_sync`) is called inside the timed region for GPU backends; the
   reference and Warp CPU backends are synchronous. Warp's INFO banner is silenced
   (`quiet_warp`) so it cannot pollute a JSON report on stdout.

5. **First benchmark — deterministic proton CSDA depth dose**
   (`benchmarks/bench_depth_dose.py`). A 150 MeV proton pencil beam in water, 0.5 mm
   bins, **straggling off** so the workload is deterministic and its output is a fixed
   digest. It times the core longitudinal transport kernel on reference, Warp CPU and
   Warp CUDA, records per-backend throughput and speedups, and gates the **per-history**
   depth dose with the established V1 metric — the edge-aware cumulative difference
   `max|cumsum(a)−cumsum(b)|/Σa` (decision `0009` /
   `validation/v1_depth_dose_csda.py`), tight to `1e-4` reference-vs-Warp and `1e-5`
   CUDA-vs-CPU (plus a `5e-3` per-bin term on CUDA-vs-CPU). The comparison is on
   per-history dose, not unit-normalised shapes, so a uniform-scale divergence is
   caught. The reference is a scalar Python oracle (~10² histories/s) timed at a small
   history count; the Warp backends are timed at a larger count. Throughput is a
   per-history rate, so the numbers stay comparable — with the caveat that the GPU is
   under-utilised at these sizes.

## Validation (`benchmarks/`, `tests/ionmc/test_benchmarking.py`)

- The harness primitives are unit-tested on the reference path (provenance keys,
  timing statistics and warm-up/repeat counts, digest determinism and sensitivity,
  cross-backend agreement arithmetic, backend/sync selection, numpy-aware JSON).
- `bench_depth_dose.py` runs on the host runner (CPU + CUDA) and reports the physics
  gate plus the timing/throughput table; the physics gate — not the timings — is the
  pass/fail criterion.

## Consequences

- There is now a reproducible, provenance-stamped way to measure performance and to
  demonstrate unchanged physics, unblocking Stage 6 optimisation work and V6 regression
  tracking. First result: Warp CPU is ~10²× the reference per-history rate on the
  deterministic depth-dose kernel (a baseline, not a target).
- **Deferred:** GPU occupancy/scaling sweeps (history-count and grid-size scans),
  benchmarks of the 3-D scattering / dose / LET / fluence / multi-ion / fragmentation
  paths, stochastic-throughput benchmarks, memory-footprint measurement, a persisted
  results series for regression tracking over time, `float32`-vs-`float64` accumulation
  precision studies, and any actual kernel or memory-layout optimisation.
