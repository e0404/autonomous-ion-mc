# 0042 — Warp CPU/CUDA backends and cross-backend parity criteria

- Status: accepted
- Date: 2026-09-30
- Task: V2-005
- Affects: architecture, performance, reproducibility, validation strategy

## Problem

The accelerated paths must execute the same physical model as the reference
path (EXPERIMENT.md) without per-particle host loops (V1-MUST-035), with
parallel-safe random streams (V1-MUST-028), and be compared with statistical
criteria (V1-MUST-029/040). Decision 0001 fixed deterministic kernel-level
parity; stochastic transport output needs its own criterion.

## Decision

1. **Kernel structure.** One generic Warp kernel (`transport_kernel`)
   transports one primary per thread from birth to death in a bounded loop,
   calling the shared step physics of decision 0041 through the Warp-bound
   namespace. The same source instantiates float32 (default) and float64
   (falsification) kernels via `type(x)(literal)`.
2. **Streams.** `wp.rand_init(seed, batch_offset + history)` gives each
   history an independent counter-based stream; primaries are sampled on the
   host from `SeedSequence(seed).spawn(batches)` exactly as in the reference
   backend. Results are deterministic per (backend, precision, seed) up to
   float32 atomic ordering (verified ≤ 1e-5 relative).
3. **Scoring/accounting.** Per-batch device grids with atomic adds in the
   kernel precision, folded into float64 host tallies per batch; energy
   accounting uses float64 atomics. Boundary crossings use a precision-safe
   overshoot of 2.5e-4 mm (float32 resolution at 500 mm is 6e-5 mm); with the
   former 1e-6 mm float32 particles stalled in zero-length steps.
4. **Parity criterion (stochastic).** Between two independent runs (different
   seeds/backends), per-voxel z = (a − b)/√(σ_a² + σ_b²) over voxels above
   2 % of the maximum: std z < 1.3 (< 2.0 for small CI cases), max |z| < 4.5
   (< 5), and |mean z| < 0.6 (< 1.0 for small cases) — the mean-z bound is
   loose on purpose because depth bins of one run are correlated (every
   history feeds all bins; a first version with |mean z| < 0.25 failed at
   +0.28 with 100 000 CUDA histories although R80 agreed to 0.02 mm); tight
   checks use integral quantities: energies per primary equal within 1e-4,
   R80 within 0.15 mm (0.3 mm for small cases), plateau mean within 0.5 %,
   steps per history and cutoff deposits within 5 %. Same-seed runs across
   precisions are correlated and are only checked for consistency.
5. **Timing.** Kernel time is measured around launch + device synchronisation
   per batch; host time is the remainder of the batch loop; the first launch
   time (module load/compile) is reported separately.

## Defect found by the parity strategy (retained)

The first Warp kernel reproduced the reference bit-for-bit with straggling
and scattering disabled, but with straggling enabled it deposited 0.4 %
more in the plateau, lost 0.09 MeV per primary through the distal face and
disagreed with the deterministic plateau, while the reference agreed to
0.04 %. Isolation showed a Warp semantics defect: `wp.func` arguments are
passed by value, so the RNG state advanced inside the shared Gamma sampler
never advanced in the caller and every step of a history reused one draw
(single-thread test: two consecutive calls returned identical values).
Fix: the shared samplers return `(sample, state)` and callers rebind the
state; regression tests in `tests/ionmc/test_warp_shared_physics.py`.
After the fix the straggling-only plateau ratio Warp/deterministic is
0.9996–1.0004. The earlier CPU probe of the sampler had passed because each
thread drew once; the failure was only visible in transport, which is why
the parity suite compares full profiles, not primitives.

## Evidence (this task, 100 MeV protons in water, 1 mm transport / 1 mm depth scoring)

- Deterministic transport (no straggling/MCS): reference and Warp CPU
  float64 identical in all 100 depth bins and step counts (112.00).
- Straggling only / MCS only: plateau ratios 0.9996 and 1.0000 between the
  reference (3 000 histories) and Warp (60 000).
- Reference (2 000 histories) vs Warp CPU float32 (2 000, all physics):
  z mean +0.13, std 0.84, max |z| 2.4 over 159 bins; R80 77.55 vs 77.54 mm.
- Warp CPU float32 vs float64 (independent seeds, 60 000 each): z mean
  +0.01, std 1.08.
- CUDA parity and float64 evidence: recorded in the task validation record
  from the host-runner runs (see PR and `validation/task-history.json`).
- CUDA (host run `RUN-20260930T030324Z-664e5f1d`): 100 000 float32
  histories vs 20 000 Warp CPU: z mean +0.28, std 1.05, max |z| 2.4;
  totals equal; float32 vs float64 CUDA consistent; deterministic per seed.
- Throughput (electromagnetic physics only, `benchmarks/workloads/proton_3d_em.py`,
  host run `RUN-20260930T030345Z-6925ff2d`, 150 MeV protons in the
  120 × 120 × 300 mm envelope, 1 mm transport, 2 mm scoring, 10 batches,
  304 steps/history): RTX A6000 float32 10⁶ histories in 0.58 s transport
  wall (1.71 × 10⁶ histories/s), 0.72 s end-to-end warm; 10⁵ histories in
  0.11 s; first launch (module compile) 1.75 s; host RSS 484 MiB. Warp CPU
  single thread ≈ 1.4 × 10⁴ histories/s; Python reference ≈ 3 × 10²
  histories/s. These are *not* target measurements: nuclear interactions,
  secondaries and LET scoring are missing.

## Tradeoffs / open items

No compaction or secondary stacks yet (one primary per thread; divergence
is acceptable for electromagnetic-only protons). Float32 atomics on a
per-batch grid are order-dependent at the ulp level. Performance targets are
not claimed: nuclear physics and LET are absent (see `benchmarks/workloads`).
