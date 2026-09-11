# Performance benchmarking

Stage 6 introduces a reproducible performance-benchmarking suite (decision
[`0030`](../generated/decisions/0030-benchmark-harness.md)). It is the basis for
validation milestone **V6** — *documented benchmark results with unchanged scientific
validation outcomes before and after each optimisation*.

## Harness and drivers

The reusable machinery lives in the source module `ionmc.benchmarking` (unit-tested
like any other source code); runnable **drivers** live under `benchmarks/` and use it,
just as `validation/` scripts use source modules.

`ionmc.benchmarking` provides:

- `capture_provenance()` — `ionmc`/Python/NumPy versions, platform, CPU count, git
  SHA, and the Warp version and available devices, stamped on every report;
- `measure(fn, ...)` — untimed warm-up (excluding one-time Warp kernel compilation)
  then timed repeats via `time.perf_counter`, reporting min/median/mean/stdev and a
  throughput (`work_units / median`, e.g. histories/s). A device-synchronisation
  callback (`make_sync`) is timed inside the region for asynchronous Warp CUDA, so the
  GPU is measured honestly;
- `array_digest()` / `relative_agreement()` — a reproducible scientific fingerprint of
  a result and the cross-backend relative difference;
- `new_report()` / `dump_report()` — a provenance-seeded report and numpy-aware JSON.

## What a benchmark asserts

A benchmark records **wall-clock time, which is hardware-dependent and is never
asserted against an absolute threshold**. What each driver *does* gate is the physics:
it pins a scientific digest of its output and checks cross-backend agreement on
**per-history** results (which retain absolute magnitude — never unit-normalised
shapes, which would hide a uniform-scale error) against the established depth-dose
budget (decision `0009` / `validation/v1_depth_dose_csda.py`). A change that silently
breaks the physics fails the benchmark; a change that only affects speed does not. This
is the concrete mechanism behind the V6 "unchanged scientific outcome" requirement.

## First benchmark

`benchmarks/bench_depth_dose.py` times the core longitudinal transport kernel — a
150 MeV proton pencil beam in water, 0.5 mm bins, energy straggling off so the workload
is deterministic — on the reference, Warp CPU and Warp CUDA backends. It records
per-backend throughput and speedups and gates the per-history depth dose with the
established edge-aware cumulative metric (1e-4 reference-vs-Warp, 1e-5 CUDA-vs-CPU).
The reference is a scalar Python oracle (~10² histories/s) timed at a
small history count; the Warp backends are timed at a larger count (throughput is a
per-history rate, so the numbers stay comparable, with the GPU under-utilised at these
sizes). Run it on the host runner:

```bash
env PYTHONPATH=/workspace/src python benchmarks/bench_depth_dose.py \
    --require-cuda --cache-dir /cache/ionmc
```

## 3-D dose benchmark + scaling sweep

`benchmarks/bench_dose3d.py` (decision `0031`) times the **representative
treatment-planning workload**: a 150 MeV proton pencil beam through a 120×120×300
(1 mm) voxel box scored into a 60×60×150 (2 mm) 3-D dose grid via ray/voxel DDA
traversal and per-step atomic dose scoring — the hot kernel (expected memory-bound),
deterministic
(scattering/straggling off). The per-history 3-D dose is gated cross-backend on the
established V4 budget (integral ≤ 1e-5, worst voxel ≤ 5e-3; decision `0021`).

Beyond the fixed-size timing it runs a **history-count scaling sweep** on the Warp
backends (CUDA `{10³,10⁴,10⁵,10⁶}`, CPU `{10³,5·10³,2·10⁴}`), reporting throughput at
each size and the per-backend peak — characterising GPU utilisation and saturation,
which the small-workload 1-D benchmark deliberately under-utilises. Run it on the host
runner:

```bash
env PYTHONPATH=/workspace/src python benchmarks/bench_dose3d.py \
    --require-cuda --cache-dir /cache/ionmc
```

## Deferred

GPU occupancy/scaling sweeps, benchmarks of the 3-D scattering / dose / LET / fluence /
multi-ion / fragmentation paths, stochastic-throughput benchmarks, memory-footprint
measurement, a persisted results series for regression tracking, `float32`-vs-`float64`
accumulation precision studies, and any actual kernel or memory-layout optimisation.
