# Performance calibration (V2-002)

Purpose: measure what the workstation, NVIDIA Warp and the registered
reference engines deliver *before* any IonMC physics exists, so that the
numeric performance targets frozen in `validation/release-plan.json`
(decision 0039) are derived from measurements rather than guessed or copied.
Nothing here is release evidence and no v1 timing is used.

## Contents

| File | Contents |
|---|---|
| `primitives.py` | Synthetic transport-like step kernel benchmark (Warp CPU/CUDA, float32/float64, with/without atomic scoring), cold/warm compile and init timing, memory. |
| `record-2026-09-30.json` | The calibration record: hardware, primitive measurements (host run `RUN-20260929T224600Z-098609df`), single-thread reference-engine throughput for the proton-3d envelope (TOPAS `REF-9b1409ed46826c2c82b9-ebf89c46`, MCsquare `REF-79a06b4029c303bdaae2-cf425f13`, FRED `REF-0ffa92d1284ad4b51b89-354a38ab`) and the step-count estimates used to derive targets. |

## Reproduction

```bash
# inside the controlled host runner (GPU), from a clean committed worktree:
python benchmarks/calibration/primitives.py --clear-cache --devices cpu cuda:0 \
  --particles 10000 100000 1000000 --steps 200 --repeats 3 \
  --output benchmarks/generated/primitives-cold.json
```

Reference-engine timings come from the committed case bundles under
`validation/references/calibration/` executed through
`run_reference_calculation` (single thread each).

## Method notes

* Timing is wall time around `wp.launch` + `wp.synchronize_device`, best of
  three launches after a first compile launch; the first-launch time after
  `warp.clear_kernel_cache()` is the cold compile time.
* The synthetic kernel is a ceiling, not a prediction: it has no species
  branching, no secondary stack, no voxel-boundary stepping and a single
  table. Real transport is expected to be one to two orders of magnitude
  slower per step; decision 0039 applies explicit margins.
* Warp CPU kernels run on one host thread; CUDA float64 on the RTX A6000 is
  about 25x slower than float32.
