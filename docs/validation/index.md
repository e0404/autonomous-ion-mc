# Validation

This section documents software, numerical, physical, statistical, and
cross-backend validation. The pages under this section describe the
diagnostic instruments; this page lists the scientific validation performed
so far and the rules it follows.

## Rules

- Acceptance quantities, metrics, tolerances and reference data are fixed in
  the decision that opens a stage, before the comparison is run (roadmap,
  *Rules applied to every validation milestone*).
- Decision `0001` covers deterministic Warp CPU-vs-CUDA parity by expectation
  class; decision `0005` adds the criterion for the float64 reference path
  versus the float32 Warp kernels; stochastic transport output will need its
  own criterion (Stage 1).
- Every result is recorded for the exact committed SHA with
  ``record_local_validation``. Scripts under ``validation/`` run on the
  controlled host runner and depend only on numpy, Warp and the standard
  library.

## V0 — analytical stopping power (task DEV-002)

Script: ``validation/v0_stopping_power.py``. Reference: the 20 native PSTAR
grid values for protons in liquid water kept in
``ionmc.reference_data.pstar_water`` (NIST PSTAR via the Apache-2.0 MCsquare
repository, commit ``211eefe6``; two independent fetches agreed). Warp
results below are from host run ``RUN-20260910T080455Z-e666b469`` (RTX A6000,
CUDA 12.9, Warp 1.17.0) at task SHA ``d23c78a``, recorded in the DEV-002
local validation record; the reference-path results were reproduced in the
sandbox as well.

| comparison | criterion (decision 0006) | result |
|---|---|---|
| analytic (Python, float64) vs PSTAR, 17 points from 10 to 400 MeV | ≤ 1.0 % each | max 0.075 %, all points low by 0.01–0.07 % |
| analytic vs PSTAR at 5 / 2 / 1 MeV | regression guards 2 / 5 / 10 % | +0.90 / +3.46 / +5.20 % |
| CSDA range 100–250 MeV vs recalled PSTAR values | 0.5 % guard (informational) | within 0.04 % |
| R(200) − R(100) vs recalled 18.242 g/cm² | 0.2 % | +0.02 % |
| numpy binding vs Python binding | bitwise equal | equal |
| Warp CPU / CUDA (float32) vs Python reference, 400 energies 2–400 MeV | S: rtol 1e-5; range: rtol 2e-5 | S 4.2e-7 on both devices (normalized 0.042); range 2.4e-7 (normalized 0.012) |
| Warp CPU vs CUDA | S: rtol 4e-6, atol 1e-6 (transcendental class); range: rtol 1e-5, atol 1e-6 (iterative) | S max abs 1.5e-5 MeV cm²/g = 1 float32 ULP at 2 MeV (normalized 0.043); range 9.5e-7 g/cm² (normalized 0.006) |

The ICRU 90 mean excitation energy (78 eV) is reported as a separate model
difference (−1.14 % at 1 MeV to −0.43 % at 400 MeV), not absorbed into any
tolerance.

The V0 milestone of the roadmap is closed by task ``DEV-003`` (tabulated
layer), which will replace the recalled range values by table-integrated
ones from the acquired dataset.

## Unit and regression tests

``tests/ionmc/`` covers units and constants, materials, the RNG mirror
(pinned to a kernel dump from Warp 1.17.0 on the workstation), the binding
mechanism, and the stopping-power model on the reference paths; Warp tests
are skipped where Warp or CUDA is unavailable. GitHub CI runs the suite
without a GPU.
