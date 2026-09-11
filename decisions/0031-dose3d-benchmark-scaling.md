# 0031 — 3-D dose benchmark and the GPU history-count scaling sweep

- Status: accepted
- Date: 2026-09-11
- Task: DEV-026
- Affects: performance engineering (Stage 6), validation strategy (V6)

## Problem

The first benchmark (decision `0030`, `bench_depth_dose.py`) times the 1-D CSDA
depth-dose kernel, which is compute-bound and — timed at small history counts to
keep the scalar reference tractable — leaves the GPU under-utilised, so its CUDA
throughput is a floor, not a representative number. The treatment-planning-relevant
hot path is the **3-D voxel-grid dose**: ray/voxel DDA traversal writing into an
independent lab-frame dose grid by per-step atomic scatter (decisions `0020`,
`0021`). That kernel is memory-bound (atomic voxel accumulation) and is what the
throughput requirement (10⁵–10⁶ protons/s) actually refers to. It needs its own
benchmark, and the benchmark needs to be run at GPU-saturating sizes to report a
meaningful rate.

## Decision

1. **Second benchmark — deterministic 3-D dose** (`benchmarks/bench_dose3d.py`,
   reusing the `ionmc.benchmarking` harness). A 150 MeV proton pencil beam through a
   120×120×300 (1 mm) voxel box, scored into a 60×60×150 (2 mm) `DoseGrid3D`, with
   **scattering and straggling off** so the per-history 3-D dose is deterministic and
   a fixed digest.

2. **Physics gate on the per-history 3-D dose** to the established V4 float32
   budget (decision `0021` / `validation/v4_dose3d.py`): reference-vs-Warp and
   CUDA-vs-CPU each tight on the integral (`total_rel_diff ≤ 1e-5`) and the worst
   voxel (`max_bin_rel_diff ≤ 5e-3`). The cumulative (depth-ordered) metric used for
   the 1-D benchmark is not meaningful for a 3-D grid, so the integral + worst-voxel
   pair is used instead — the same pair the V4 dose gate uses. Comparisons are on
   per-history dose (absolute magnitude), never unit-normalised.

3. **History-count scaling sweep.** Each Warp backend is timed across a sweep of
   history counts — CUDA `{10³,10⁴,10⁵,10⁶}`, CPU `{10³,5·10³,2·10⁴}` (CPU runs
   serially, so it saturates far lower) — reporting throughput (histories/s) at each
   size and the per-backend peak. This characterises GPU utilisation and saturation
   and quantifies the under-utilisation caveat of the small-workload 1-D benchmark;
   the fresh dose grid allocated per run keeps the physics gate and the sweep from
   leaking accumulation between runs (its ~2 MB allocation is inside the timed region,
   a constant across backends and a legitimate part of a per-run dose calculation).

## Validation (`benchmarks/bench_dose3d.py`, `tests/ionmc/test_bench_dose3d.py`)

- On the host runner (CPU + CUDA) the physics gate passes and the sweep reports the
  throughput curve; the physics gate — not the timings — is the pass/fail criterion.
- Reference-path tests guard the driver's own logic (per-history exactness = 150 MeV,
  digest determinism, dose-grid independence between runs) in CI without a GPU.

## Consequences

- The representative memory-bound hot kernel now has a reproducible benchmark and a
  throughput-vs-size curve, giving a meaningful single-GPU protons/s figure and the
  data an optimisation would move.
- **Deferred:** scattering-on (stochastic) throughput, grid-size / voxel-count and
  beamlet-count scaling dimensions, occupancy/roofline profiling, a persisted
  regression series, and any actual kernel or memory-layout optimisation.
