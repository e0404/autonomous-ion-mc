# Warp architecture measurement scripts (research report `docs/research/warp-architecture.md`)

Scripts used by the performance-specialist agent on 2026-10-03 to obtain the
measurements quoted in the Warp architecture research report, preserved so that
every number in that report can be reproduced. They are **not** part of the
`ionmc` package and the toy kernel is a cost stand-in, not validated physics.

## Reproducing the measurements

```bash
uv sync --extra dev
PYTHON="$PWD/.venv/bin/python" bash validation/scripts/warp-architecture/run_all.sh \
    validation/scripts/warp-architecture/results/<date>-<host>
```

`run_all.sh` records the exact command line of every step, bounds each step
with `timeout` (`STEP_TIMEOUT`, default 600 s; a timed-out step is archived with
exit code 124), archives stdout and stderr verbatim as `NN-<name>.txt`, writes
`environment.txt` (Python, Warp, numpy, CPU, kernel, exact git SHA, dirty state
of the script directory and SHA-256 of every script) and finally generates
`SUMMARY.md` from the raw files with `summarize.py`. The script exits non-zero
if any step failed or timed out. `STEPS="06 07"` reruns selected steps.

## What each step measures

| Step | Script | Measures |
|---|---|---|
| 01, 02 | `compile_probe.py {1,0}` | CPU cold compile time of the toy module with/without `enable_backward`, using a fresh empty `WARP_CACHE_PATH` per process |
| 03, 04 | `run_toy.py {f32,f64} 20000 1` | single-thread throughput of the toy kernel (20 000 histories of 150 MeV protons), float32 vs float64 (the default kernel cache is warm here; cold compile is steps 01–02) |
| 05–07 | `run_toy.py f32 {80000,160000,320000} {4,8,16}` | throughput with concurrent CPU launches (20 000 histories per thread, private accumulators) |
| 08 | `run_toy.py f32 20000 1` | repeat of step 03 (repeatability) |
| 09 | `precision.py` | float32 accumulation error vs deposits per voxel, 40-batch split, float32 position drift over 300 steps, short-step energy-loss cancellation |
| 10 | `../rng/rng_seed_dupes.py` | histories with identical `wp.rand_init` start states between seeds s and s+1 |
| 11–13 | `../rng/rng_overlap.py N L` | fraction of `wp.randf` draws revisiting a 32-bit state already used by another history |
| 14 | `../rng/philox.py` | Philox4x32-10 Random123 known-answer vectors (kernel and Python adapter), kernel-vs-Python agreement on 20 000 random blocks, cost per uniform |

## Archived runs

The authoritative measured values are the raw files under `results/<run>/` and
the `SUMMARY.md` generated from them. Numbers quoted in the research report's
prose were the agent's unarchived scratch measurements of the same scripts and
are superseded by the archive wherever they differ.

Provenance notes:

- The agent's original scratch versions of `run_toy.py` and `compile_probe.py`
  were corrected before archiving: `run_toy.py` refuses a history count that
  is not divisible by the thread count (the original silently simulated fewer
  histories while normalising by N) and reports `simulated=`;
  `compile_probe.py` creates a fresh, empty kernel cache per process so the
  measured time is a cold compile.
- Steps 05–07 simulate a different history set per configuration, so they
  measure throughput only. Invariance of results under partitioning is a
  property of the package implementation to be tested there (decision 0037);
  these scripts do not establish it.
- During the first attempts on 2026-10-03 a concurrent-launch step (4 or 8
  Python threads issuing Warp CPU launches) stalled indefinitely once and was
  killed externally; reruns completed in well under a second. The cause was
  not identified. Concurrent Warp CPU launches from Python threads are
  therefore recorded as a hang risk for the CPU backend design (prefer
  process-based parallelism or launch timeouts).
- `rng_overlap.py` and `philox.py` use the deprecated `warp.config.quiet` flag
  (a deprecation warning is printed); results are unaffected.
- Expected GPU throughput figures in the research report are hypotheses, not
  measurements.
