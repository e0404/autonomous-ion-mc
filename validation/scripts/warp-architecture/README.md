# Warp architecture measurement scripts (research report `docs/research/warp-architecture.md`)

Scratch scripts used by the performance-specialist agent on 2026-10-03 to
obtain the measurements quoted in the Warp architecture research report. They
are preserved so that every number in that report can be reproduced; they are
**not** part of the `ionmc` package and the toy kernel is a cost stand-in, not
validated physics. Environment of the original measurements: Warp 1.17.0,
numpy 2.5.3, Python 3.12, Intel i9-13900K (WSL2), CPU device only. Run from
this directory with an environment containing `warp-lang` and `numpy`.

| Script | Measures | Result quoted in the report |
|---|---|---|
| `toy_transport.py` | Builds the cost-representative toy transport kernel (3D DDA, three log-table lookups, Philox draws, Gaussian straggling, Highland-type scattering, nuclear sampling, per-batch float32 atomics); `build(real)` and `tables(np_dtype)` | – |
| `run_toy.py f32|f64 N [threads]` | Throughput of the toy kernel for `N` 150 MeV primaries; with `threads > 1` runs concurrent launches with private accumulators | 75k histories/s (1 thread, f32), 55k (f64); 301k/514k/841k at 4/8/16 threads; identical totals across thread counts |
| `compile_probe.py` | Cold compile time of the toy module with and without `enable_backward` | 1.45 s / 1.22 s |
| `precision.py` | float32 sequential accumulation error vs deposits per voxel; 40-batch split; position drift over 300 float32 steps; short-step energy-loss cancellation | 1.3e-3 relative at 3e6 deposits; 3.5e-7 with 40 batches; ≤ 1.7e-4 mm drift; 44 % error for E−E_new at ds = 1e-4 mm, 150 MeV |
| `pyscope_cost.py` | Cost of calling shared `@wp.func` physics from Python scope | ≈ 300 µs per partial step |
| `probe_pyscope.py`, `pyprec.py`, `pyprec2.py` | Which constructs work in Python-scope `@wp.func` execution and their precision | array indexing, `rand_init`, `uint32` ops fail; `float`-typed functions mix precisions; `wp.float64`-typed functions match numpy |
| `recur.py`, `structt.py` | Recursion rejection at code generation; struct arrays with `vec3` members | recursion rejected; 32-byte struct works |
| `philox_lib.py` | Philox4x32-10 `@wp.func` factory imported by the toy kernel (copy of `../rng/philox_lib.py`) | – |

Raw console output of the original runs was not archived; the numbers above
are the values the agent reported. Re-running the scripts reproduces the
measurements on the same software versions (throughput depends on hardware).
Expected GPU throughput figures in the report are hypotheses, not measurements.
