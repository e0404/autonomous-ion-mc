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

## V0 — tabulated stopping power and range (task DEV-003)

Script: ``validation/v0_tabulated_stopping_power.py``. Data: the MCsquare
PSTAR and Geant4 water tables through the data layer (decision `0007`),
acquired into the host-runner cache. Reference: the analytic model
(decision `0006`) and NIST PSTAR CSDA ranges for liquid water (retrieved
2026-09-10). Criteria fixed in decision `0008` before the comparison.

| comparison | criterion | result (reference paths) |
|---|---|---|
| tabulated (PCHIP) reproduces the committed 20-point PSTAR subset | rtol 1e-12 | exact |
| tabulated vs analytic, 10-400 MeV | ≤ 1.0 % | max +0.25 % (at 11 MeV) |
| tabulated CSDA range + 0.5 MeV residual vs NIST PSTAR (100/150/200/250 MeV) | ≤ 0.1 % | ≤ 0.018 % |
| analytic CSDA range + 1 MeV residual vs NIST PSTAR | ≤ 0.5 % | ≤ 0.053 % |
| numpy binding vs Python binding | bitwise equal | equal |
| Warp CPU / CUDA (float32) vs Python reference | S rtol 1e-5; range rtol 2e-5 | S 5.8e-8, range 4.5e-8 (both devices) |
| Warp CPU vs CUDA | transcendental class (rtol 4e-6, atol 1e-6) | bitwise identical (no transcendentals) |

Warp results are from host run `RUN-20260910T100641Z-80a76b4a` (RTX A6000, CUDA 12.9, Warp 1.17.0) at
task SHA `dfd9995`, recorded in the DEV-003 local validation record.

This closes milestone V0: proton stopping power in water is now available from
both an analytical model and an I-value-consistent external table, on all three
execution paths, validated against reference data with pre-fixed tolerances.

## V1 (partial) - CSDA proton depth dose (task DEV-004)

Script: ``validation/v1_depth_dose_csda.py``. Reference: the tabulated CSDA
range (task DEV-003). Criteria fixed in decision `0009` before the comparison.
This is the deterministic, no-straggling foundation of milestone V1; the
straggling/Bragg-shape and scattering/lateral parts follow in later Stage-1
tasks.

| check | criterion | result (reference paths) |
|---|---|---|
| energy conservation (reference, float64) | \|balance\| ≤ 1e-9 | ~1e-16 |
| R80 vs tabulated CSDA range, 100/150/200 MeV | ≤ 0.3 % | ≤ 0.15 % |
| step-size convergence (fixed grid, 0.02 vs 0.002) | drift ≤ 0.05 % | ≤ 3e-4 % |
| range cross-check (stopping depth vs CSDA range) | ≤ 0.2 % | ≤ 0.15 % |
| reference vs Warp CPU/CUDA, cumulative depth dose | ≤ 1e-4 of total | ≤ 2.7e-5 (both devices) |
| Warp CPU vs CUDA, cumulative depth dose | ≤ 1e-5 of total | ≤ 2.3e-8 (near bitwise) |

Warp results are from host run `RUN-20260910T135739Z-8166c2f4` (RTX A6000, CUDA 12.9, Warp 1.17.0) at
SHA `b726507`, recorded in the DEV-004 local validation record.

## V1 (straggling) - Bragg peak and range straggling (task DEV-005)

Script: ``validation/v1_bragg_straggling.py``. References: the analytic Bohr
range-straggling integral and Bortfeld's ``sigma = 0.012 R^0.935`` (decision
0010). The metric is the range straggling ``sigma_R`` (std of stopping depths),
which isolates straggling from the peak shape.

| check | criterion | result (reference paths) |
|---|---|---|
| mean range vs CSDA range, 100/150/200 MeV | within 0.1 % | <0.01 % |
| sigma_R vs analytic Bohr integral | within 3 % | <1 % |
| sigma_R vs Bortfeld 0.012 R^0.935 | within 10 % | +5-7 % (Bohr electron-binding gap) |
| sigma_R / R band | 0.9-1.2 % | 1.04-1.11 % |
| energy conservation (reference float64 / Warp float32) | 1e-9 / 5e-5 | ~1e-16 / ~1e-5 |
| Bragg peak present (peak > 3x entrance, deep) | yes | yes |
| two independent seeds statistically consistent | RMS(t) in [0.5, 1.6] | ~1.1 |
| reference vs Warp / Warp CPU vs CUDA cumulative | 1e-4 / 1e-5 | 5.8e-7 / 8.1e-8 |

Warp results are from host run `RUN-20260910T144031Z-8dd34007` (RTX A6000, CUDA 12.9, Warp 1.17.0) at
SHA `100653d`, recorded in the DEV-005 local validation record.

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
