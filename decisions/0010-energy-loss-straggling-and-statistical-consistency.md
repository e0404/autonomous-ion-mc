# 0010 — Energy-loss straggling, step/deposition model, and statistical consistency

- Status: accepted
- Date: 2026-09-10
- Task: DEV-005
- Affects: physical accuracy, numerical accuracy, validation strategy, reproducibility, performance, scientific interpretation

## Problem

DEV-004 built deterministic continuous-slowing-down proton transport with a
sharp (unphysical) distal edge. This task adds **energy-loss (range)
straggling**, giving a realistic Bragg peak, and fixes, before the first
comparison: the straggling model, how a stochastic step deposits energy, the
criterion for **statistical consistency** of Monte Carlo output (the item the
roadmap defers from decision 0001 to "when first needed" — now), and the
acceptance targets for the straggled depth dose. Multiple Coulomb scattering
and lateral spread remain a later task; this task is longitudinal only.

## Context and evidence

Research (Claude physics-researcher subagent, 2026-09-10; sources cited below)
and in-repository numerical checks established:

- **Bohr energy-straggling variance** per unit path length
  `dOmega^2/dx = K m_e c^2 (Z/A) rho z^2 f(beta)`, with the Bethe constant
  `K m_e c^2 = 0.307075 x 0.5109989 = 0.156915 MeV^2 cm^2/mol` and the
  relativistic factor `f(beta) = (1 - beta^2/2)/(1 - beta^2)` (Leo §2.7;
  ICRU 49). `f` is +11 % at 100 MeV and +24 % at 200 MeV, far larger than any
  straggling tolerance, so it is **not** optional
  (https://link.springer.com/chapter/10.1007/3-540-31718-X_8).
- **Range straggling** `sigma_R^2 = integral_0^{E0} (dOmega^2/dx)/S(E)^3 dE`
  (S the linear stopping power). The in-repo integral gives sigma_R/R =
  1.11 / 1.07 / 1.04 % at 100 / 150 / 200 MeV, matching the published
  ~1.0-1.1 %.
- **Bortfeld (1997)** analytic monoenergetic straggling width
  `sigma_mono = 0.012 R^0.935` (R, sigma in cm): 0.809 / 1.585 / 2.515 mm at
  100 / 150 / 200 MeV (https://aapm.onlinelibrary.wiley.com/doi/10.1118/1.598116).
  The pure-Bohr integral runs ~5-7 % above Bortfeld because it omits the
  electron-binding (Bethe-Livingston) reduction; this sets the physics
  tolerance.
- **Per-step Gaussian is valid by the CLT**: individual steps are in the
  Vavilov regime (kappa ~ 0.1) and skewed, but the distal falloff is the sum
  of hundreds of steps, so the accumulated range straggling is Gaussian with
  variance `integral dOmega^2/dx` regardless of per-step skew, provided no
  single step dominates (step energy loss <~ 1-2 %). This is what MCsquare,
  FRED and MOQUI do.

## Selected approach

### Straggling model (shared source `ionmc.physics.transport`)

Per step of length `dl`, the mean energy loss `dE_mean` (midpoint, DEV-004)
gets a Gaussian fluctuation:

`dE = clamp(dE_mean + xi * sigma, 0, E)`, `sigma = sqrt(dOmega^2/dx * dl)`,

with `xi` a standard-normal variate drawn in the **execution layer** (the
kernel's `wp.randn`, the reference driver's `ionmc.rng` mirror; the RNG is not
called from shared source — decision 0005) from the per-history counter-based
stream. `bohr_straggling_sigma` and `straggled_energy_loss` are shared source.
Below `straggling_floor_mev` (default 2 MeV) straggling is switched off and the
step is deterministic, so the terminal-step clamp cannot bias the range; the
residual range there (<~0.1 mm) is negligible.

### Step / deposition decoupling (correction to DEV-004)

DEV-004 limited the transport step to the depth-bin boundary. With straggling
that makes `sigma` per step comparable to `dE_mean` (fine steps), and the
`[0, E]` clamp then biases the range short — measured **0.9 % short at
200 MeV**, and grid-dependent, which is unphysical. Therefore the **physics
step** is now the fractional-energy-loss step (capped by `max_step_mm`,
default 1 mm), *not* bin-limited, and its deposited energy is **distributed
across the depth bins the step spans**, proportional to the path length in
each bin. This decouples the stochastic step from the scoring grid; the mean
range is then grid-independent and unbiased (verified: mean range = CSDA range
to <0.01 % at 100-200 MeV). The scoring accumulator is **float64** on the
device (float32 tallies promoted, as decision 0005 anticipated) so energy
conservation is limited only by the float32 physics (~1e-5), not by summation.

### Statistical consistency criterion (new; complements decision 0001)

Decision 0001 covers deterministic kernel parity. For **stochastic** depth
dose, two estimates are compared by their Monte Carlo uncertainty:

1. **Independent estimates** (different seeds), each with a per-bin standard
   error `se_i` from `n_batches >= 2` batches (`TransportEngine.run_batched`):
   the standardized residual `t_i = (A_i - B_i) / sqrt(se^A_i^2 + se^B_i^2)`
   is ~N(0,1) under consistency. The gate: over bins with mean dose above 1 %
   of the peak, `RMS(t_i)` in `[0.5, 1.6]` and `max|t_i| <= 5` (a two-sided
   allowance for a few hundred bins). The integral dose and R80 agree within
   their combined standard error.
2. **Cross-backend, shared streams** (reference vs Warp, Warp CPU vs CUDA):
   because the counter-based streams are identical across backends
   (decision 0005), the histories are the *same* up to float32 rounding, so
   agreement is *stronger* than statistical and is measured, as in decision
   0009, by the cumulative depth dose: reference vs Warp `<= 1e-4` of the
   total, Warp CPU vs CUDA `<= 1e-5`. The range straggling `sigma_R` agrees
   across backends within the batch standard error.

### Acceptance criteria (fixed before the comparison), milestone V1 (straggling part)

| check | criterion |
|---|---|
| mean range vs tabulated CSDA range, 100/150/200 MeV | within 0.1 % (straggling is unbiased) |
| range straggling sigma_R vs the analytic Bohr integral | within 3 % |
| range straggling sigma_R vs Bortfeld `0.012 R^0.935` | within 10 % (documents the Bohr electron-binding gap) |
| sigma_R / R stays in 0.9-1.2 % over 100-200 MeV | yes |
| energy conservation (reference float64) | \|balance\| <= 1e-9 |
| energy conservation (Warp float32) | <= 5e-5 |
| Bragg peak-to-entrance ratio present (peak dose > entrance) | yes (a Bragg peak exists) |
| statistical consistency of two independent seeds | RMS(t) in [0.5, 1.6], max\|t\| <= 5 |
| reference vs Warp CPU/CUDA cumulative depth dose | <= 1e-4 |
| Warp CPU vs CUDA cumulative depth dose | <= 1e-5 |

## Rationale

Bohr+f(beta) Gaussian straggling is the accepted therapy-MC model and the CLT
makes the per-step Gaussian exact for the accumulated distal width, which is
the validated quantity. Validating on the **range straggling sigma_R measured
from stopping positions** (not the depth-dose edge width) isolates straggling
from the Bragg-peak shape and from the integral-vs-volumetric-dose distinction,
giving an unambiguous target against both the model's own analytic integral
(a tight sampling-correctness check) and Bortfeld's empirical law (a physics
check). The step/deposition decoupling is the standard condensed-history
structure and removes a grid-dependent bias that fine bin-limited steps
introduce once straggling is present.

## Expected tradeoffs

- Bohr straggling is ~5-7 % above the empirical (electron-binding-reduced)
  value; the +10 % tolerance covers it. A Bethe-Livingston correction could be
  added later as an empirical scale if a reference calculation demands it.
- Per-step Gaussian, not Landau/Vavilov: correct for the accumulated width by
  the CLT, not for single-step energy-loss spectra (not needed for therapy
  protons).
- Deposition splitting adds a short bounded inner loop per step; energy
  conservation in float32 is ~1e-5 (the split fractions are float32), which is
  scientifically negligible.
- Longitudinal only: no lateral spread yet, so this is still an integral depth
  dose, not a volumetric dose; multiple Coulomb scattering is the next task and
  closes milestone V1.

## Validation strategy

`tests/ionmc/test_straggling.py` (Bohr sigma formula, unbiased mean range,
sigma_R vs analytic and Bortfeld, Bragg peak present, batch uncertainty,
statistical consistency of two seeds, cross-backend cumulative on the Warp
paths) and `validation/v1_bragg_straggling.py` on the host runner for the exact
SHA (Warp CPU and CUDA), recorded below.

## Later validation outcome

DEV-005 host run `RUN-20260910T144031Z-8dd34007` (RTX A6000, CUDA 12.9, Warp 1.17.0, SHA `100653d`,
40000 histories): all eight gates passed. Mean range = CSDA range to 5e-5 at
100/150/200 MeV (straggling unbiased). Range straggling sigma_R = 0.853 / 1.687
/ 2.705 mm vs the analytic Bohr integral 0.860 / 1.694 / 2.712 mm (ratio
0.992 / 0.996 / 0.997) and Bortfeld 0.811 / 1.582 / 2.520 mm (ratio 1.052 /
1.067 / 1.073 — the expected +5-7 % electron-binding gap); sigma_R/R =
1.11 / 1.07 / 1.04 %. Energy conservation ~8e-9 (reference float64 residual in
the deposit) to ~5e-7 (Warp float32). A Bragg peak with peak/entrance ~6.8.
Two independent seeds statistically consistent: RMS(t) = 1.10, max|t| = 4.1.
Reference vs Warp CPU cumulative depth dose 5.8e-7; **Warp CPU vs CUDA 8.1e-8**,
with sigma_R agreeing to seven digits (1.68423368 vs 1.68423374 mm) — the
shared counter-based streams make the two devices sample the same histories up
to float32 rounding.
