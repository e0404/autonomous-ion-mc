# 0041 — Condensed-history transport algorithm for ions (electromagnetic core)

- Status: accepted (nuclear interactions, secondaries and LET are separate decisions)
- Date: 2026-09-30
- Task: V2-004
- Affects: physical accuracy, numerical accuracy, reproducibility, performance, validation strategy

## Problem

Choose the stepping, energy-loss, straggling and multiple-scattering
algorithms for planning-oriented ion transport such that (a) the physics is
written once and executed by the Python reference path and the Warp
kernels, (b) results are step-size and grid independent within the frozen
falsification bounds, and (c) the GPU implementation remains branch-light.

## Evidence

Research record `validation/research/2026-09-30-condensed-history-methods.md`
(fast-MC practice: range-table inversion for the mean loss, Bohr Gaussian
straggling, per-step Highland's step dependence and its remedies, Gottschalk's
differential Molière scattering power). Prototype measurements in this task:

- A separately tabulated inverse range table biased the per-step mean loss by
  0.5 % at 0.2 mm steps and produced negative losses for the µm-scale steps
  that the random hinge creates; exact inversion of the forward log-log range
  table removed both (loss error 0.0000 % at all step sizes; R80 identical at
  1 mm/10 % and 0.2 mm/2 % steps).
- With exact inversion, 100 MeV protons in water give R80 = 77.57 mm with
  straggling versus the ICRU 90 CSDA range of 77.59 mm; energy is conserved
  to 1e-9 with escaped and cutoff energies accounted separately.
- Warp's Python-scope evaluation rounds transcendental builtins to float32
  and rejects mixed float64/literal arithmetic, so a math namespace injected
  into a twice-loaded module is used instead of Python-scope Warp calls;
  `type(x)(literal)` makes generic kernels precision-agnostic.

## Decision

1. **Single-source physics.** `ionmc/transport/step_physics.py` holds every
   step formula as `@wp.func` code in the Warp subset with a `MATH`
   namespace; `ionmc.transport.shared.physics("python")` binds it to Python's
   `math` (float64 reference), `physics("warp")` to Warp builtins (CPU/CUDA,
   float32 or float64 generic kernels). Literals are written `type(x)(c)`.
2. **Step limits.** step = min(distance to the transport-voxel face + 1e-6 mm,
   `max_step_mm` (default 1 mm), path over which the CSDA energy drops by
   `energy_step_fraction` (default 10 %)). Particles enter the geometry by a
   vacuum flight to the bounding box; leaving it is counted as escaped.
3. **Mean energy loss** by exact inversion of the piecewise log-log range
   table (walk-down search from the current cell), so the mean loss over a
   step is consistent with the table for any step length.
4. **Straggling**: Gamma distribution with the CSDA mean and Bohr variance
   0.1569 z_eff² (Z/A) (1 − β²/2)/(1 − β²) ρs MeV², sampled by
   Marsaglia–Tsang; strictly positive, mean-preserving, Gaussian in the
   large-κ limit. Per-step loss spectra are not Landau/Vavilov shaped; the
   documented domain is integrated quantities (depth dose, ranges).
5. **Multiple scattering**: Gaussian polar angle with mean square
   z² f_dM (E_s/pv)² ρs/(ρX_S) from Gottschalk's differential Molière
   scattering power (E_s = 15 MeV, ρX_S from the Bragg-rule scattering
   length), uniform azimuth, applied at a uniformly random fraction of the
   step (random hinge) for unbiased lateral displacement. No single-scattering
   tail (documented limitation, to be evaluated against TOPAS/measured
   profiles).
6. **Cutoff**: below 0.5 MeV/u the residual energy is deposited locally and
   counted as `deposited_cutoff` (never silently merged).
7. **Scoring**: energy deposited along a step is distributed over the
   scoring voxels crossed, proportionally to path length, so scoring grids
   may differ from the transport grid in resolution and alignment; dose uses
   exact overlap masses.
8. **Statistics**: independent batches with `SeedSequence(seed).spawn`
   streams; standard error of the batch mean; NaN when undefined.
9. **Fail-closed contract**: `ionmc.config.check_supported` rejects any
   physics/backend/precision/scorer request the backend does not implement
   before transport starts; `nuclear=True` and secondaries raise until the
   corresponding tasks land.

## Tradeoffs

The reference path runs ≈ 4 × 10⁴ steps/s (adequate for 10²–10⁴ history
validation cases only). The Gaussian MCS core underestimates wide-angle
tails; heavy-ion MCS uses the z² scaling without an ion-specific correction.
Straggling below κ ≈ 1 per step uses a Gamma shape rather than Landau.

## Validation

`tests/ionmc/test_reference_transport.py`: energy conservation and
accounting, R80 versus ICRU 90 CSDA within the frozen tolerance, step-size
independence, grid shift/coarsening conservation, capability contract,
result reopening in a clean process, lateral σ at 30 mm depth versus
Fermi–Eyges (local), carbon range (local). Independent Monte Carlo evidence:
TOPAS electromagnetic-only case
`validation/references/cases/topas-proton-100mev-em-only` (run ID in the
validation record of this task).
