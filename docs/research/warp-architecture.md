# Warp transport architecture recommendation (IonMC v3)

Author: performance specialist subagent. Status: advisory input to the lead; none of this is validated physics.

## Measurement context

All numbers in this report were obtained in the agent's sandbox on an Intel
i9-13900K (32 logical CPUs, WSL2) with Python 3.12.0, numpy 2.5.3 and
warp-lang 1.17.0 on the CPU device only (no CUDA driver visible). The toy
kernel is `validation/scripts/warp-architecture/toy_transport.py`; it is
cost-representative but **not validated physics**. Per step it does a 3D DDA in
a 120×120×300 1 mm grid with four materials, three log-spaced table lookups,
8 Philox uniforms, Gaussian straggling, Highland-type multiple Coulomb
scattering with a rotation, nuclear sampling with a secondary stack, and two
float32 atomics into a 10-batch 2 mm scoring grid, for 150 MeV protons at about
160 steps per primary.

*[Lead note, 2026-10-03: the measured values quoted in the prose of this
report are the agent's unarchived scratch measurements. The authoritative
measured values are the archived raw outputs and generated `SUMMARY.md` under
`validation/scripts/warp-architecture/results/`, produced by `run_all.sh` at a
clean committed SHA with the corrected scripts; where the two differ, the
archive supersedes the prose. The agent's claim that particle counts and
energy deposits were identical across 1/4/8/16 threads is withdrawn: the
archived configurations simulate different history sets, so partition
invariance is pending until the package repartitions the same histories and
compares them (decision 0037). The report's magnitudes (≈75k–80k histories/s
per CPU thread, ≈1.3× float32/float64 ratio, ≈10× scaling at 16 threads,
≈1.2–1.5 s cold CPU compile, 20 % RNG state reuse at 1e6×2000 draws) are
confirmed by the archive.]*

## 0. Recommendations in brief

1. Write the physics as **pure, precision-generic `@wp.func`s on scalars and vectors**. They receive already-fetched table values and already-drawn uniforms. Table fetch and RNG are thin backend adapters.
2. Use one **"whole particle life" kernel**: one thread transports one particle until it dies. Secondaries go into a **global generation queue**, which is drained by a short host loop over *generations*, not over steps.
3. Use float32 for transport state and tables in production, plus a **float64 build of the same source** for V2-NUM probes and reference parity.
4. **Do not use `wp.rand_init`/`wp.randf`.** Use a counter-based Philox4x32-10 keyed by genealogy. It is bit-identical between kernel and reference (verified).
5. Score into **float32 per-batch accumulators** and reduce in float64. A single float32 grid is measurably biased at 1e6 histories.
6. Warp CPU kernels are **serial and their atomics are non-atomic**. CPU parallelism must come from concurrent launches with private accumulators.

## 1. Shared physics across three backends (incorporating the coordinator's correction)

These Python-scope facts were verified in Warp 1.17:
- `wp.array` item indexing fails.
- `wp.rand_init` fails.
- `uint32` `>>` and `^` fail (`unsupported operand type(s) for >>: 'uint32' and 'uint32'`).
- Scalars, `vec3`, `while` loops, structs and **local `wp.types.vector(length=16)` stacks with dynamic indexing** work.
- Precision is subtle. A function annotated `float` mixes precisions in Python scope: plain arithmetic runs in float64, but builtins such as `wp.log` return float32-rounded values, and `vec3` is float32. A function typed `wp.float64` evaluates in clean float64, matching numpy bit-for-bit in the probe.

The architecture therefore has four layers:

```
physics/   pure @wp.func, generated per precision R (float32|float64) by a factory
  log_bin(x, lx0, inv_dl, nb) -> vec2(i, frac)
  interp(y0, y1, frac) -> R
  energy_loss(e, r_old_E, r_new_E, s, S_e) -> R    # telescoping + short-step branch
  highland_theta0(e, m, z, s, x0) -> R
  gauss_pair(u1, u2) -> vec2
  rotate(u, cos_t, phi) -> vec3
  plane_distance(p, u, i, dx) -> R
  nuclear_distance(u, sigma_macro) -> R
  fragment_kinematics(..., u0..u3) -> ...
adapters/
  warp:   fetch(tab: wp.array2d, m, i) (@wp.func); philox4x32_10 (@wp.func)
  python: fetch(np_tab, m, i) = R(tab[m, i]); philox_py (Python ints)
drivers/
  warp:   @wp.kernel transport_generation(...)  - thin step loop, calls physics/*
  python: def transport_particle(...)           - same loop, same RNG draw order
```

The Python reference calls the **float64 variant** of every shared function with `wp.float64`/`vec3d` arguments, so its precision is well defined. It uses numpy only for table storage.

Lookups must be composed as `log_bin` (shared) → adapter `fetch` of `y0` and `y1` → `interp` (shared). That way only the memory read differs between backends.

The driver step loop is the one piece written twice. Keep it short, with the step-level physics decisions inside shared functions. Guard it with a **trajectory-level parity test**: because Philox is bit-identical across backends and the draw order is fixed, reference(float64) and Warp CPU(float64) must agree step by step for the first K histories to about 1e-12 until a branch flips. This is far more discriminating than distribution-level parity. Cost is of order 1 ms per full step in Python (agent's scratch estimate from `pyscope_cost.py`, not archived), so the reference handles about 1e2–1e4 histories per process. Use multiprocessing for statistical parity runs.

```python
def make_physics(R):
    V3 = wp.types.vector(length=3, dtype=R)        # alias at factory scope; calling
    @wp.func                                       # wp.types.vector(3,R)(...) inside a
    def interp(y0: R, y1: R, f: R) -> R:           # kernel fails ("same type" error)
        return y0 * (R(1.0) - f) + y1 * f
    ...
    return SimpleNamespace(interp=interp, ...)
PHYS32, PHYS64 = make_physics(wp.float32), make_physics(wp.float64)
```

## 2. Particle state layout and precision

**Layout.** In a whole-life kernel the hot state (position, direction, energy, indices, RNG counter) lives in registers, so array of structs versus structure of arrays matters only for the queues. A `wp.struct` with `vec3` members works in arrays (32 B in the test). For queue records, prefer **SoA arrays** (`pos: vec3`, `dir: vec3`, `e: float32`, `w: float32`, `species: int8/int16`, `gen_id: vec2ui`, `tag: int32` = batch|beamlet). Adjacent threads read adjacent elements, which coalesces on CUDA, and arrays can be added without changing a struct ABI.

**Precision (measured):**
- **Position.** Over 300 float32 steps with MCS rotations, the final position error versus float64 was 3.5e-5 mm median and 1.7e-4 mm maximum (2000 paths). This is negligible. The real hazard is boundary logic (see §6).
- **Energy loss.** E − E_new in float32 suffers cancellation on short steps at high energy: at 150 MeV, ds = 1e-4 mm the relative error is 44%, and at 1e-2 mm it is 0.19%. DDA corner crossings produce such steps. The first toy version froze at E = 144.590164 MeV indefinitely, because forward and inverse range interpolants were not mutual inverses: E − E_new went negative and was clamped to 0.
  - Rule 1: compute ΔE as a **difference of the same inverse-range interpolant**, `Rinv(R) − Rinv(R − s)`. It is monotone, never negative, and telescopes, so per-history energy is conserved.
  - Rule 2: for s ≪ R (for example s < 1% of R) use `S(E_mid)·s`.
  - Both rules need V2-NUM step-refinement tests.
- **Accumulation.** See §4.

**Policy.** Use float32 transport state and tables, float32 atomics into per-batch grids, and float64 batch reduction and outputs. Build a float64 variant of the *same* kernels from the factory. The float64 build costs about 1.35× on CPU (archive-derived medians 77.8k vs 57.4k histories/s over three repeats, see `SUMMARY.md`); on an A6000 the FP64 ALU rate is 1/64 of FP32, so expect a much larger factor. It serves as the V2-NUM float32/float64 probe and as the trajectory-parity partner of the reference. It is not a production mode on CUDA.

## 3. Kernel organisation

**Recommendation for v1: a per-particle persistent loop (one thread runs one particle to death) and generation queues.**

```python
@wp.kernel
def transport_generation(q_in: ParticleQueue, n_in: int, q_out: ParticleQueue,
                         n_out: wp.array(dtype=wp.int32), geo: Geometry, tabs: Tables,
                         score: Scorers, ctl: Control, counters: wp.array(dtype=wp.int64)):
    t = wp.tid()
    if t >= n_in:
        return
    st = load_particle(q_in, t)                      # registers from here on
    key = wp.vec2ui(ctl.seed, ctl.stream)            # stream: batch/beamlet-independent
    ctr = wp.uint32(0); n_sec = int(0); steps = int(0)
    while st.alive:
        if steps >= ctl.max_steps:                    # bounded, reported, never silent
            tally_truncation(counters, score, st); break
        r = philox4x32_10(wp.vec4ui(st.gen_id[0], st.gen_id[1], ctr, 0), key); ctr += 1
        ...                                           # step physics via shared funcs
        if produced_secondary:
            slot = wp.atomic_add(n_out, 0, 1)
            if slot < q_out.capacity:
                store(q_out, slot, child(st, n_sec)); n_sec += 1
            else:
                wp.atomic_add(counters, OVERFLOW, wp.int64(1))   # fail closed on host
        steps += 1
```

The host loop runs `for gen in range(max_generations): launch; swap queues; read n_out`. That is about 3–10 launches per batch, not one per step. This satisfies the "no Python loop per particle or step" requirement, which a step-wise wavefront design with host-side compaction would strain: ~300+ launches per batch plus a compaction pass per step.

Why not step-wise kernels in v1?
- Event-based or wavefront designs pay off when kernels are huge and divergent, as in Geant4-on-GPU work.
- For condensed-history protons, lifetimes within a beam are coherent (similar range), and per-step branching (nuclear interaction ≈1e-3/mm) is rare.
- The whole-life loop is the simplest correct design.

**Divergence control.** Sort each generation's queue by energy with `wp.utils.radix_sort_pairs` on a quantised energy key before launch. This optimisation can wait until profiling shows the need. Query `wp.get_cuda_kernel_properties(k)` for `register_count`/`local_memory_size` and `wp.get_suggested_block_size(k)` on the host runner.

**Loop bounds.**
- Cap steps at `max_steps = ceil(k · (path_bound/min_step))`, for example 4·(nx+ny+nz) plus the range-limited count, padded generously.
- Also enforce a *progress* guard: any step of length 0 must advance a voxel index.
- Tally truncated particle count **and residual energy** into dedicated counters.
- Fail-closed rule: any step truncation, queue overflow or genealogy overflow invalidates the result in every mode; the truncated energy is tallied separately and never enters dose (V2-NUM/V2-OUTPUT). *[Lead correction, 2026-10-03: the agent's original text allowed production runs to merely report truncation; decision 0037 requires invalidation in all modes.]*

**Language limits verified in Warp 1.17.**
- `while`, `break`, `continue` and local fixed-size vectors with dynamic indexing work.
- Recursion is rejected at codegen ("maximum recursion depth exceeded"), so a secondary cascade cannot be recursive.
- There is no dynamic allocation in kernels, so all buffers are preallocated.
- `atomic_add` supports `[u]int32`, `[u]int64`, `float32`, `float64`; `atomic_or/and/xor` support integers only.
- Literal gotcha: `wp.uint32(2891336453)` is emitted as an `int32` constant (clang `-Wconstant-conversion` warning). The bits survive, but prefer `wp.constant(wp.uint64(...))` and explicit casts.

## 4. Secondaries

Options:
- **(a) Per-thread local stack:** `vecN` per component, measured to work. On CUDA a dynamically indexed local vector goes to local memory (L1-cached). 7×8×4 B = 224 B per thread, or ≈29 MB at 129k resident threads. It is fine for protons, but bounded depth forces an overflow policy and lengthens warp tails.
- **(b) Global generation queue:** recommended for all species.

Memory for the queue: a record is about 48 B (pos 12, dir 12, E 4, w 4, species 2, gen_id 8, tag 4, pad).
- Protons: about 0.2–0.5 transported secondaries per 150 MeV primary. Capacity 1e6 × 1 × 48 B ≈ 48 MB.
- Carbon at 290 MeV/u: nuclear interactions occur in about half of primaries over 30 cm, each with a few charged fragments. Capacity about 5e6 × 48 B = 240 MB, still trivial on 48 GB.

Size buffers from a calibration run multiplied by a safety factor. On overflow, re-launch the batch with doubled capacity. Never drop silently.

Determinism: the secondary's RNG identity is its **genealogy ID** (`parent_gen_id`, `child_index`), not its queue slot. Queue order from atomics is scheduling-dependent; physics results are not, apart from float summation order.

## 5. Scoring

**float32 accumulation (measured, sequential float32 adds of 0.5–8 MeV deposits into one voxel):**

| Deposits per voxel | float32 relative error | MC statistical error |
|---|---|---|
| 1e5 | 4.7e-6 | 1.6e-3 |
| 1e6 | 1.6e-5 | 5e-4 |
| 3e6 | **1.3e-3** | 2.9e-4 |
| 1e7 | **3.2e-3** | 1.6e-4 |

The hottest Bragg-peak voxel at 1e6 primaries is in the 1e6–3e6 range, so a single float32 grid is biased beyond the statistical error. The same 1e7 deposits spread over **40 batch accumulators** gave 3.5e-7.

Recommendation:
- Accumulator `[n_batches, n_vox]` float32, with batch = `history_id % n_batches`. The interleaving also spreads atomic contention.
- After transport, reduce to float64 mean and variance (batch method).
- Offer float64 atomics as a probe mode. sm_86 has native float64 `atomicAdd` in L2; its relative cost must be measured on the host runner.
- CPU is safe because each concurrent launch owns its accumulators.

**Memory** (scoring grid 60×60×150 = 540k voxels, 2.16 MB per float32 field):

| Use | Fields | 10 batches | 40 batches |
|---|---|---|---|
| proton-3d: dose, LETd numerator (Σ edep·L), LETt numerator, track-length fluence | 4 | 86 MB | 346 MB |
| carbon: per-species (Z=1..6) dose + LETd numerator | 12 | 259 MB | 1.04 GB |
| influence: 100 beamlets × dose | 100 | 2.16 GB | (8.6 GB, unnecessary) |

LETd needs two accumulators: Σ edep·L and Σ edep. For a batch estimate of the ratio, use the per-batch ratio estimator, or delta-method or jackknife on the batch sums. Document which, and flag voxels with too few batches having nonzero dose as **undefined, not 0** (V2-NUM).

**Influence workload.**
- Run all 100 beamlets × 1e4 histories in **one launch per generation**. `history_id → beamlet = id // 1e4`, and batch = `(id % 1e4) % nb`.
- Write into a dense `[nb·100, n_vox]` float32 accumulator: 2.16 GB at 10 batches, so that is the beamlet-chunk size limit.
- Sparsify on device: one pass for the per-beamlet max, a count of entries above `thr·max`, `wp.utils.array_scan` for offsets, then a fill. That gives a CSC matrix with float32 values and int32 voxel indices.
- At 1–5% occupancy that is 5k–27k nonzeros per beamlet, about 4–22 MB total.
- Reduction and sparsify passes read ≤2.2 GB, about 3 ms at the A6000's 768 GB/s, which is negligible versus transport. Transport time equals that of proton-3d (same 1e6 histories).
- Persist the threshold, the dropped-dose fraction per beamlet, the beamlet table, histories, batches and units (V2-OUTPUT).
- For thousands of beamlets, chunk beamlets so the dense buffer stays ≤ about 8 GB.

## 6. Geometry

- Transport grid: material `uint8`/`int16` plus density `float32`, 4.32 M voxels at about 22 MB.
- **Use an incremental DDA.** Voxel indices are *state*. The distance is to the index-defined plane `(i+1)·dx − p` or `i·dx − p`, clamped ≥ 0. On a boundary-limited step, increment the index of the crossed axis. Do **not** use `floor(p)` with an ε nudge along the path. The first toy did that and stalled at 1e-4 mm steps whenever a direction component was tiny, because the nudge along u does not move the coordinate across the plane. Recompute `floor(p)` only for new particles.
- The **scoring grid** is independent (own origin, spacing and shape). Compute its index from the step's midpoint, `(p_mid − origin_s)·inv_ds`.
  - v1: deposit at the midpoint *without* splitting at scoring boundaries. Steps are ≤1 mm, about half the 2 mm scoring voxel.
  - A misaligned or finer scoring grid is a V2-NUM probe: compare against an optional "split at scoring planes" mode, which is a second DDA on the scoring grid, as a falsification test.
- Arbitrary beam direction falls out of the 3D DDA. The source adapter generates (p0, u0) outside or on the phantom surface and ray-enters the grid. Test rotated beams for invariance (V2-NUM).

## 7. Tables

- Use `wp.array2d(dtype=float32)` indexed `[material, bin]`, log-spaced in E (and log R for inverse range), with 512–2048 bins.
- Precompute per material: S(E), R(E), R⁻¹, macroscopic σ_inel(E), MCS scattering constants, and species-resolved data with a leading species axis for ions.
- Cost per lookup is one `log`, two reads and one lerp. Reads of a few hundred KB of tables are L1/L2-resident on the GPU, and lanes in a warp hit nearby bins. On CPU the measured ≈83 ns per step includes three lookups.
- Precompute `lx0` and `inv_dl` per table.
- Interpolate R⁻¹ by the same rule used to build R. Test the round-trip residual as part of table qualification.

## 8. Compilation and caching

- Cold CPU compile of the toy module is about 1.2–1.4 s, with `enable_backward=False` about 0.2 s faster (archive-derived medians, see `SUMMARY.md`; the agent's scratch values 1.45/1.22 s are non-evidentiary). Disable backward for all transport modules because the kernels are not differentiated.
- **CUDA NVRTC times are unmeasured here and are a risk.** Measure them on the host runner.
- Treat physics *model selection* as runtime flags (ints in a `Control` struct) when the branches are cheap.
- Use `wp.static`/factory closures only for precision (f32/f64) and genuinely structural variants (species families). Keep the variant count ≤ about 4, and put each in its own module (`module="unique"` or separate Python modules) so a change recompiles one variant.
- The cache honours `WARP_CACHE_PATH` (`/cache/warp` on the runner) or `wp.config.kernel_cache_dir`.
- **Cold start** must use a fresh empty cache directory per measurement, for example `WARP_CACHE_PATH=/cache/warp-cold-<run-id>` or `wp.clear_kernel_cache()`. Otherwise you are measuring a warm start.

## 9. Measurement plan (APIs verified in the installed package)

- `cold_start_seconds`: process start → first result available on host, with an empty cache. Read `time.perf_counter()` at the first line of `__main__` and run the launcher wall clock externally.
- `compilation_seconds`: `t = perf_counter(); wp.load_module(module, device=dev)` (or `wp.force_load(device)`) on an empty cache.
- `wall_seconds`: the same, for the end-to-end API call including output persistence. Report warm and cold separately.
- `gpu_kernel_seconds`: either
  - `e0 = wp.Event(enable_timing=True)`, `wp.record_event(e0)` … `wp.get_event_elapsed_time(e0, e1)` (ms) around the launch sequence; or
  - `with wp.ScopedTimer("transport", synchronize=True, cuda_filter=wp.TIMING_KERNEL) as t:` then sum `t.timing_results`.
- `host_seconds`: wall time minus device-busy time. Also report host phases (setup, upload, reduce, write) with `ScopedTimer(synchronize=True)`.
- `histories_per_second`: primaries divided by wall time and by kernel time (report both).
- `gpu_peak_mib`:
  - `wp.get_mempool_used_mem_high(dev)` gives Warp allocations only. It has no reset API, so use one fresh process per workload.
  - Also sample `dev.total_memory - dev.free_memory` before and after setup and at peak. This captures context, module and local-memory reservation, but WSL2 numbers can be noisy.
  - Report both.
- `host_peak_mib`: `resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024` (KiB on Linux), in a fresh process.
- Hygiene: three or more repeats, report median and range; fixed seeds; record GPU clocks and driver; no concurrent jobs.
- Record the scaling axes from `performance.json`.

## 10. Expected throughput and main risks

**Published context.** These are rates as reported by the authors. Hardware generations differ, so they are not comparable without normalisation.
- gPMC: 1e7 protons in 6–22 s on a Tesla C2050, about 0.5–1.7 M/s ([arXiv:1409.8336](https://arxiv.org/pdf/1409.8336) summarises Jia et al. 2012).
- FRED: about 1e7 primaries/s on a single GPU ([Schiavi et al., PMB 2017](https://research.uniroma1.it/node/43406); [Acta Phys. Pol.](https://www.actaphys.uj.edu.pl/R/48/10/1625/pdf)).
- MCsquare: 1e7 protons in under a minute on many-core CPUs ([Souris et al., Med Phys 2016](https://dial.uclouvain.be/pr/boreal/en/object/boreal%3A150338)).

**Ballpark for this design (hypotheses to calibrate, not targets).**
- Full proton physics is perhaps 2–4× the toy's per-step cost, i.e. Warp CPU of order 20–40k histories/s per thread and 0.2–0.5 M/s on 16 threads (extrapolation from the archive-derived toy throughput, not a measurement).
- RTX A6000: about 1–5 M primaries/s for 150 MeV protons in water, so the 1e6 kernel time is about 0.2–1 s.
- Wall time is then dominated by Python/Warp start-up (≈1–2 s), cold NVRTC compile (expected multi-second to tens of seconds, unmeasured), data upload and output writing.
- Carbon at 290 MeV/u: 3–10× slower per primary because of fragment transport.
- Influence: about proton-3d plus less than 0.1 s of assembly.

**Main risks.**
1. RNG quality, if Warp's built-in RNG is used (see §11; mitigated by Philox).
2. float32 accumulation bias (mitigated by batches).
3. Short-step energy-loss cancellation and range-table inconsistency (§2).
4. DDA stalls (§6).
5. CUDA compile time and register/local-memory pressure of a large whole-life kernel. Measure `register_count` early.
6. Warp-tail divergence from long-lived fragments.
7. Serial Warp CPU semantics: plain read-modify-write atomics, so sharing accumulators between concurrent CPU launches corrupts results.
8. Silent queue overflow or step truncation unless counters are fail-closed.
9. The duplicated driver loop diverging between backends; mitigated by trajectory-level parity.

## 11. RNG (measured)

Warp's generator (`native/rand.h`):
- It is the Jarzynski–Olano **PCG hash** iterated on a **32-bit state**: `state = pcg(state)`.
- `rand_init(seed, i) = pcg(seed + pcg(i))`.
- `randf` returns 24-bit floats in [0, 1 − 2⁻²⁴]. It can return 0.0, so `log(u)` is −inf; use `(x>>8 + 0.5)·2⁻²⁴`.

All streams share one 2³² state space. Marking visited states in a bitmap:

| Workload | Draws | Draws revisiting another history's state |
|---|---|---|
| 1e5 histories × 1000 draws | 1e8 | 1.1% |
| 1e6 × 1000 | 1e9 | 10.7% |
| 1e6 × 2000 (≈ proton-3d) | 2e9 | **20.0%** |

Between batches seeded `s` and `s+1`, **260 of 1e6 histories are bit-identical duplicates**. This held for seed pairs (1,2), (42,43) and (1000,1001).

This violates "parallel-safe streams" (REQUIREMENTS) and would contaminate batch-variance estimates.

**Use Philox4x32-10** (Random123; it passes BigCrush):
- counter = (gen_id.hi, gen_id.lo, step counter, sub-draw);
- key = (seed, stream type).

Verified:
- Random123 known-answer vectors pass in Warp.
- The pure-Python adapter is bit-identical on 20,000 random blocks.
- Cost per uniform on one CPU thread: Philox is slightly cheaper than `wp.randf` (archive-derived medians in `SUMMARY.md`, step 13; the agent's scratch values 2.09 vs 2.42 ns are non-evidentiary).
- Results are expected to be invariant to thread count and partitioning because the generator is counter-based; this is to be tested in the package (see lead note above).

Each Philox call returns 4 uniforms. A step that needs more draws makes a second call with a different sub-draw index.
