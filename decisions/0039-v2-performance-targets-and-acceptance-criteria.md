# 0039 — Frozen v2 performance targets and physics acceptance criteria

- Status: accepted and frozen at the commit that carries this record (the
  freeze SHA is recorded in `validation/task-history.json` and the
  `targets_frozen` telemetry event); later changes require a versioned
  amendment that preserves this record and the original results.
- Date: 2026-09-30
- Task: V2-002
- Affects: performance, validation strategy, scientific interpretation

## Problem

`experiment/v2/performance.json` fixes three workload envelopes and the
workstation but leaves numeric targets `null` and unratified. The protocol
requires justified numeric targets and frozen acceptance observables and
tolerances *before* substantial optimization and before qualification
outcomes are seen, and forbids deriving targets from stripped physics or
from v1 timing.

## Evidence (calibration record `benchmarks/calibration/record-2026-09-30.json`)

Hardware: NVIDIA RTX A6000 (48 GiB, sm_86; CUDA 12.9, driver 13.4 under
WSL2), Intel i9-13900K (32 logical CPUs), 32 GiB host memory, Warp 1.17.0.

Primitive ceilings (synthetic condensed-history step with RNG, table lookup,
rotation and atomic 3-D deposit; effective particle-steps per second):

| device | precision | 10⁵ particles | 10⁶ particles |
|---|---|---|---|
| CUDA | float32 | 4.3 × 10⁹ | 2.3 × 10¹⁰ |
| CUDA | float64 | 8.3 × 10⁸ | 9.3 × 10⁸ |
| CPU (1 thread) | float32 | 3.6 × 10⁷ | – |

Cold compile of the small benchmark module: 0.84 s CUDA, 1.3 s CPU;
`warp.init` is negligible. Host RSS 434 MiB.

Reference engines, single thread, 150 MeV protons in 120×120×300 mm water,
2 mm scoring, each engine's full default physics: TOPAS/Geant4 11.4.2
(QGSP_BIC, 0.1 mm cuts, electrons transported) 1.07 × 10³ histories/s;
MCsquare 6.3 × 10⁴ histories/s; FRED (CPU) 3.0 × 10⁴ histories/s with 33 %
secondary particles.

Step estimates: a 150 MeV proton has a 157.6 mm CSDA range; with ≤ 1 mm
steps and ≤ 10 % energy loss per step plus secondaries this is ≈ 230
steps/history, so proton-3d and influence are ≈ 2.3 × 10⁸ steps; carbon at
290 MeV/u with transported fragments is estimated at ≈ 5 × 10⁸ steps.

## Candidate approaches

1. Targets at the primitive ceiling (e.g. 10⁶ histories in < 0.1 s):
   unattainable with real physics; would make the gate meaningless.
2. Targets at reference-engine CPU throughput: MCsquare already achieves
   6 × 10⁴ histories/s on one thread, so a GPU implementation merely matching
   it would not be competitive and would not distinguish a GPU-resident
   design from a host-loop design.
3. **Targets derived from the ceiling with explicit margins for physics
   complexity, checked against fast-MC literature (selected).** Real
   transport adds species/energy branching, several tables, nuclear sampling,
   secondary stacks, voxel-boundary stepping and batch launches; a factor
   50–200 below the synthetic ceiling is the plausible band, giving
   10⁵–4 × 10⁵ histories/s for protons. Published GPU proton MCs report
   10⁵–10⁶ histories/s on older hardware, so 10⁵ histories/s is a floor a
   competitive implementation must clear, not a stretch goal.

## Frozen targets (backend `warp-cuda`, full enabled physics, ≥ 10 batches)

Definitions: `wall_seconds` — end-to-end wall time of the documented public
run (configuration load, device setup, transport, statistics, persisted
outputs) with a warm kernel cache; `histories_per_second` — primary
histories divided by transport wall time (host + kernel, all batches);
`cold_start_seconds` — process start to first kernel launch with a warm
kernel cache (import, data/table load, device upload); `compilation_seconds`
— additional time of the first run after `warp.clear_kernel_cache()`;
`host_seconds` — wall time outside kernels during transport;
`gpu_kernel_seconds` — summed kernel time; `gpu_peak_mib` /
`host_peak_mib` — attributable peak device memory and peak host RSS.

| metric | proton-3d | influence | carbon-fragments |
|---|---|---|---|
| wall_seconds ≤ | 60 | 120 | 180 |
| histories_per_second ≥ | 1.0 × 10⁵ | 5.0 × 10⁴ | 2.5 × 10⁴ |
| cold_start_seconds ≤ | 30 | 30 | 30 |
| compilation_seconds ≤ | 300 | 300 | 300 |
| host_seconds ≤ | 30 | 60 | 60 |
| gpu_kernel_seconds ≤ | 20 | 40 | 80 |
| gpu_peak_mib ≤ | 4096 | 8192 | 6144 |
| host_peak_mib ≤ | 4096 | 6144 | 4096 |

Rationale per workload: proton-3d at 10⁵ histories/s is 230× below the
float32 ceiling (physics margin) and 1.6× above MCsquare's single-thread
CPU rate; influence halves the throughput floor for per-beamlet scoring,
uncertainty and sparse assembly and doubles the memory budget for 100 dense
scratch grids (216 MB) plus assembly; carbon-fragments allows ≈ 2.2× more
steps per history plus species-resolved scoring and heavier tables.
Compilation is bounded at 300 s because NVRTC compilation of a large module
with many kernels can take minutes and is cached afterwards; cold start at
30 s bounds table construction and data loading. Memory bounds keep three
concurrent workloads inside the 48 GiB device and the 32 GiB host.

Informational (recorded, not gated): the same workloads on `warp-cpu`
(single thread) and the reference Python path at reduced history counts;
the scaling series of `performance.json` (histories 10³–10⁶, beamlets
1/10/100, scoring grid 1/2/4 mm, batches 2/10/40).

Accuracy constraints attached to every performance measurement: all physics
of the workload enabled and verified from the persisted effective
configuration; batch uncertainty reported; results statistically consistent
(z-score criterion of NUM-BACKEND-PARITY) with the qualified reference-path
result of the same configuration at reduced history count.

## Frozen physics acceptance criteria

Tolerances are set from the literature record
(`validation/research/2026-09-30-*.md`), the known differences between
reference engines (Geant4 BIC/QMD fragment yields differ by 10–20 %; water
I-value lineage shifts ranges by ≈ 0.5 %) and the statistical precision
attainable with 10⁵–10⁶ histories. They are recorded in
`validation/release-plan.json` and summarized here:

- Stopping/range (all four species): tabulated stopping power within 1.5 %
  (p, He) / 3 % (C, O) of a held-out evaluation table over the documented
  domain; Monte Carlo R80 in water within max(0.5 mm, 0.5 %) (p, He) /
  max(0.7 mm, 0.7 %) (C, O) of the CSDA range of that table.
- Proton water (vs TOPAS and a non-Geant4 source): |ΔR80| ≤ 1.0 mm; distal
  80–20 % width within 0.5 mm; central-axis depth dose gamma 3 %/1 mm pass
  rate ≥ 95 % from 0 to R80 + 5 mm; lateral σ at three depths within
  max(5 %, 0.2 mm); R80 vs the non-Geant4 range reference within 1.5 mm.
- Proton nuclear: primary fraction at R80 − 10 mm within 0.03 absolute of
  TOPAS; secondary dose fraction within 30 % relative; plateau dose within
  3 %; disabling nuclear physics must change the plateau by ≥ 5 %.
- Heterogeneity: |ΔR80| ≤ 1.5 mm; 3-D gamma 3 %/2 mm ≥ 95 % vs TOPAS; a
  second engine recorded as supporting evidence.
- Helium: |ΔR80| ≤ 1.0 mm; gamma 3 %/1 mm ≥ 95 %; lateral σ within
  max(7 %, 0.3 mm); R80 vs ASTAR CSDA within 0.5 %.
- Carbon: |ΔR80| ≤ 1.0 mm; peak-to-plateau ratio within 10 %; fragment-tail
  dose (R80 + 10 … R80 + 40 mm) within 25 %; primary carbon fraction at
  R80 − 10 mm within 0.05 absolute; H + He share of tail dose within 30 %;
  measured fragment yields (if acquired) within 30 %.
- Oxygen (documented domain 100–430 MeV/u): |ΔR80| ≤ 1.5 mm; peak-to-plateau
  within 12 %; tail within 30 %.
- LET: LET_d on axis within 10 % of the TOPAS scorer of the same definition
  proximal to R80 − 5 mm and 15 % at the peak/distal region; entrance LET
  within 2 % of the analytic single-track value.
- Numerical falsification bounds, statistical coverage and backend-parity
  criteria as listed in the plan.

## Validation strategy and outcome

The plan file is validated against the qualification schema by
`tests/infrastructure/v2/test_release_plan_file.py`. Outcomes are recorded
only by the exact-SHA qualification; failures are preserved with a
`performance_failure`/`scientific_failure` event and any amendment cites
this record.
