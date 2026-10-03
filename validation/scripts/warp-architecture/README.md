# Warp architecture measurement scripts (research report `docs/research/warp-architecture.md`)

Scripts used by the performance-specialist agent on 2026-10-03 to obtain the
measurements quoted in the Warp architecture research report, preserved so that
every number in that report is reproducible. They are **not** part of the
`ionmc` package and the toy kernel is a cost stand-in, not validated physics.

## Reproducing the measurements

```bash
uv sync --extra dev
PYTHON="$PWD/.venv/bin/python" bash validation/scripts/warp-architecture/run_all.sh \
    validation/scripts/warp-architecture/results/<date>-<host>
```

`run_all.sh` records the exact command line of every step, bounds each step
with `timeout` (`STEP_TIMEOUT`, default 600 s) and archives stdout and stderr
verbatim as `NN-<name>.txt` next to an `environment.txt` (Python, Warp, numpy,
CPU, kernel, git SHA). `STEPS="06 07"` reruns selected steps.

## Archived run `results/2026-10-03-sandbox-cpu/`

Environment: Python 3.12.0, Warp 1.17.0, numpy 2.5.3, Intel i9-13900K (32
logical CPUs, WSL2), CPU device only, git SHA 196ee2b (scripts as committed in
this task). Raw outputs are the files in that directory; the table quotes them.

| Step / script | Exact command (see `# command:` line in the file) | Archived result |
|---|---|---|
| 01, 02 `compile_probe.py` | `compile_probe.py 1` / `compile_probe.py 0` (fresh empty `WARP_CACHE_PATH` per process) | CPU cold compile 1.35 s with `enable_backward`, 1.18 s without |
| 03 `run_toy.py f32 20000 1` | 20 000 histories, 1 thread, float32 | 79 198 histories/s; 160 steps/primary; `load_module` 1.38 s on first use |
| 04 `run_toy.py f64 20000 1` | same, float64 | 57 846 histories/s (float32 1.37× faster) |
| 05–07 `run_toy.py f32 {80000,160000,320000} {4,8,16}` | 20 000 histories per thread, concurrent launches with private accumulators | 303 968 / 588 768 / 822 201 histories/s (3.8× / 7.4× / 10.4× of step 03) |
| 08 `run_toy.py f32 20000 1` | repeat of step 03 | 78 428 histories/s (repeatability 1 %) |
| 09 `precision.py` | float32 accumulation, 40-batch split, position drift, short-step energy loss | 1.3e-3 relative error at 3e6 deposits, 3.5e-7 with 40 batches; drift median 3.5e-5 mm, max 1.7e-4 mm; 44 % error for E−E_new at ds = 1e-4 mm, 150 MeV |
| 10 `../rng/rng_seed_dupes.py` | duplicate start states between seeds s, s+1 | 260 of 1e6 for (1,2), (1000,1001), (42,43) |
| 11–13 `../rng/rng_overlap.py N L` | states revisited by other histories | 1.14 % (1e5×1000), 10.73 % (1e6×1000), 20.00 % (1e6×2000) |
| 14 `../rng/philox.py` | Random123 known-answer vectors, kernel vs Python adapter, cost | KAT pass (kernel and Python); 20 000/20 000 blocks identical; 1.96 ns/uniform Philox vs 2.31 ns `wp.randf` |

Notes on provenance and corrections:

- The agent's original scratch versions of `run_toy.py` and `compile_probe.py`
  were corrected before archiving: `run_toy.py` now refuses a history count
  that is not divisible by the thread count (the original silently simulated
  fewer histories while normalising by N) and reports `simulated=`;
  `compile_probe.py` now creates a fresh, empty kernel cache per process so
  the measured time is a cold compile. The archived numbers come from the
  corrected scripts, so they supersede the agent's originally quoted values
  (75k/55k/841k histories/s, 1.45 s/1.22 s), which were of the same magnitude.
- Energy per primary is 148.15–148.18 MeV in every run; the research report's
  statement that totals are *identical* across thread counts is not what this
  archive shows, because each configuration simulates a different history set.
  Invariance of results under partitioning is a property to be tested in the
  package implementation (decision 0037), not established by these scripts.
- The first execution of the complete suite on 2026-10-03 stalled at step 06
  (8 concurrent CPU launches) until the harness killed it after 40 min; the
  rerun of the same step completed in 0.27 s. The cause was not identified.
  Concurrent Warp CPU launches from Python threads are therefore recorded as
  a hang risk to be addressed (timeouts, process-based parallelism) when the
  CPU backend is implemented.
- `rng_overlap.py` and `philox.py` still use the deprecated `warp.config.quiet`
  flag (a deprecation warning is printed); results are unaffected.
- Expected GPU throughput figures in the research report are hypotheses, not
  measurements.
