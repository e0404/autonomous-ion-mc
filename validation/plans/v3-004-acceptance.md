# V3-004 acceptance criteria (frozen before any result)

Status: frozen with the first commit of task V3-004 (decision 0040), before any V3-004 test
result, local-validation (LV) or host-runner (HR) result exists. Later changes require a new dated
entry in the amendments section with the scientific reason. Criteria apply to the LET, fluence,
species, spectrum and lookup scoring of decision 0040 on the proton electromagnetic transport
engine of decision 0039; secondaries, nuclear interactions and persistence are covered by later
plans.

Tiers: CI = GitHub (small, Warp CPU, single process); LV = local exact-SHA validation on the CPU;
HR = controlled host runner with CUDA (numpy, warp and pytest only, no network, offline analytic
`StoppingSource`). Evidence classes follow `experiment/v3/PROTOCOL.md`. Chi-square and normal
quantiles come from `statistics.NormalDist` plus Wilson-Hilferty, without scipy.

## Notation

- `X_v`: the quantity of channel X summed over the pieces of voxel v (per primary). Channels: E
  (deposit), L (length), LS (`l S`), LS2 (`l S^2`), ES (`eps S`), FE, FL, N (piece count); `E_step`
  is the E channel of class "step" (decision 0040 section 2).
- `rho` is the float rounding term: 1e-12 for python and float64, 1e-6 for float32.
- `q_X` is the quantum of channel X (decision 0040 section 4) and `n_v` the piece count of voxel v
  from channel N (N is requested in every deterministic test). The fixed-point bound is
  `delta_X = n_v q_X / 2`.
- "Defined voxel": `ybar > 0` and `n_nonzero(y) >= max(2, ceil(B/2))` (decision 0040 section 5).

## Criteria

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| A1 | Synthetic constant-S_w water row, deterministic (MCS and straggling off), water and Al slabs | Every traversed voxel, for (X, Y) in {(LS, L), (LS2, LS), (ES, E_step)}: `\|X/Y - S\| <= rho S + (delta_X + S delta_Y)/Y`. Global closure `sum_v LS + res_LS = S (sum_v L + res_L)` within 1e-12 relative (float64 and python) and 1e-6 (float32). | theory | CI |
| A2 | Hook-level synthetic stream: 2 species x 2 generations, dyadic (l, S_mid, k, tau, eps) that are exact multiples of every quantum. Both through the Python hook and through a Warp-CPU test kernel calling the kernel-support `score_piece`. | Python and Warp int64 accumulators are bitwise equal. LET_t, LET_d, LET_d^eps and the species partials equal the closed forms exactly (0 tolerance). The species partials recombine to the all-species channels bitwise. Unproducible requests outside the test `producible` set raise. | theory | CI |
| A3 | Pencil beam, MCS and straggling off, aligned grid | Per voxel before the cutoff depth: `\|L_v/(N V_v) - 1/(dx dy)\| <= rho/(dx dy) + delta_L/(N V_v)`. Batch sd at most the same bound (hinge splits vary per history). Global `sum_v L_v + res_L = N chord` within 1e-12 (float64). | theory | CI |
| A4¹ | Deterministic CSDA 150 MeV, analytic water, scoring grid = transport grid, bins proximal to the cutoff | `LET_t dz = E(z1) - E(z2)` within 1e-4; `LET_d = int S dE / dE` (numpy quadrature on the same table) within 1e-3 | theory | CI |
| A4b¹ | As A4, but transport grid 1 mm, scoring grid 1.5 mm offset 0.25 mm (pieces split), 150 MeV | `\|LS_v/N - (E(z1) - E(z2))\| / (E(z1) - E(z2)) <= 1e-4` in every bin proximal to the cutoff. Negative control: the same run with k forced to 0 (test-only monkeypatch of the twin) exceeds 1e-4 in at least one bin. The numpy estimate gives 2.6e-4 for the control and 5e-8 for the ramp. | falsification | CI |
| A5 | All physics on, transport level (protons only exist) | Integer equality per batch and voxel on every backend: E(species = p) = E(gen = primary) = edep; `E_step` + `edep_excluded_from_let` = edep. The multi-species recombination is graded in A2. | conservation | CI |
| A6 | Cauchy-Schwarz | `LET_d >= LET_t (1 - rho) - Delta_v`, `Delta_v = LET_t (delta_LS2/LS2 + 2 delta_LS/LS + delta_L/L)`, in every defined voxel | consistency | CI |
| A7² | Step independence, s_max in {0.1, 0.5, 1} mm, 150 MeV | Plateau LET_d pairwise within 1 %; `\|LET_d - LET_d^ε\| <= 4 sigma_D + delta_v` (dose > 10 % of max, footnote 2); negative control: the trace-based `eps/l` estimator changes by at least 10 % | falsification | LV |
| A8 | LET recomputed offline from fluence spectra (200 bins/decade) and the water table | within 0.5 % of the online values | consistency | LV |
| A9³ | Ratio-error coverage, 200 seeds | fraction with `\|z\| < 1.96` in [0.91, 0.99] for LET_d at 3 depths, including the distal 80 % point | uncertainty | LV |
| A10 | Synthetic lookup `f = a + b S_w` on a uniform LET axis; constant energy table `f = 2` | `\|FE_v - a E_step,v - b ES_v\| <= rho FE_v + delta_FE + a delta_E + b delta_ES`. With `f = 2`, `FE = 2 E_step` bitwise. Recorded sha256 = file hash. | theory / provenance | CI |
| A11 | Backend parity | **CI:** python vs warp-cpu float64, T1 histories (K = 16, 100 MeV, all physics on), all channel kinds including spectrum and lookup: per batch and voxel `\|a - b\| <= 1e-10 max(\|a\|, \|b\|) + n_v q_c`; N channels and `n_nonzero` identical. **LV:** the same at K = 256, 150 MeV. **HR:** cpu-f32 vs cuda-f32, independent seeds, T12 layout: `t12_compare` with the frozen T12 rules applied to the linear depth profiles of L, LS, LS2, ES, E_step, FE (synthetic) and the depth-integrated spectrum (p > 0.001 and max-T p > 0.001; inconclusive fails). Channel totals of L, LS, LS2, ES: `\|z\| < 3.5` (`deterministic_scalar_verdict` is not used; it stays restricted to the T12 total deposit). LET_t and LET_d profiles get delta-method z, reported, non-gating. | backend parity | CI / LV / HR |
| A12 | Fail closed before transport | Raises for: unknown or unproducible species; generation "secondary" and any species other than the source projectile (no secondary transport); `let_medium` other than water; dose-to-water; lookup coverage, sha256 or species gaps; non-uniform table; non-uniform spectrum edges; negative or non-finite lookup values; a channel quantum above the precision floor; channels requested on a backend without them; memory above budget | contract | CI |
| A13 | **Exploratory, non-gating:** 150 MeV primaries, LET_d profile | Reported beside Grassberger and Paganetti 2011 (Phys. Med. Biol. 56:6677: distal maximum about 12 keV/um, primaries, other beam). Gating Monte Carlo evidence for LET is TOPAS ProtonLET (V3-010). No pass/fail. | related model | LV |
| A14 | Voxels without steps, or with `n_nonzero` below the threshold | NaN with `defined_mask` False, never 0 | consistency | CI |
| A15 | Partition invariance of all channels (mirrors T13) | `acc`, residual tallies and counters bit-identical: CPU chunk sizes 2^10 vs 2^18 (single process, CI); 1 vs 3 workers (CPU, `multiprocess`-marked, CI and LV); CUDA chunk sizes 2^10 vs 2^18 (HR) | self-consistency | CI / LV / HR |
| A16 | No regression of the qualified path | With `tallies = ()`, and for edep, tallies, counters and traces also with tallies present: python and warp-cpu float32/float64 outputs bit-identical to a524f209 at fixed seeds (T1 config and the T13 config) | regression | CI |

## Why A11-HR compares linear channels

1. The frozen T12 machinery (grouping into supported bins, batch permutation and bootstrap max-T,
   T12 amendments 17 and 22) is calibrated for additive batch sums; a merged group of ratio bins is
   not the ratio of the merged sums.
2. Per-batch ratios carry an O(1/n_b) Jensen bias that differs between samples with different hpb
   and are undefined where `y_b = 0` distally; both would produce spurious failures.
3. Equality in distribution of numerator and denominator implies ratio parity.
4. No new statistic needs calibration.

Per pair there are 7 profiles at p > 0.001, about 0.7 % family-wise false failure, pre-declared.

## Seeds

The fixed seed bases of V3-004 are distinct from every V3-003 base (20261004 ... 20341004):

- **Rehearsal base 20351004.** Full-scale rehearsal and diagnostic runs of the LV/HR steps. Its
  results are preserved as non-qualification evidence and never acceptance evidence.
- **Qualification base 20361004.** Used only for the final qualification runs of this task after the
  engine is frozen (A7, A9, A11-LV/HR, A15; A13 is not gating). It is recorded in the environment,
  the step documents and each sample's metadata, and the runner enforces conformance only for this
  base. Samples with different bases are never compared.
- **Consumption rule.** A base whose full-scale statistical results have been observed (including
  by an investigation or a failed qualification) is consumed: it is preserved as non-qualification
  evidence, it is never reused, and a new amendment re-freezes the qualification base to the next
  unused value `20371004`, `20381004`, ... before any result at the changed code is observed. CI
  tests use fixed small seeds (documented in the tests) and are not qualification evidence.
- Tolerances and frozen history counts do not change when a base is re-frozen.

## Single-process diagnostic mode

Operator directive 2026-10-07: until the host instability is understood, qualification runs execute
in single-process diagnostic mode. Every V3-004 LV/HR step accepts `--workers 1` with unchanged
histories, seeds and criteria. The only multiprocessing-specific check is the 1-vs-3-workers
sub-check of A15 (pytest marker `multiprocess`; runner
`DEFERRED_STEPS = {"lv": ("a15-workers",), "hr": ()}`); in this mode it is archived as deferred and
is not evidence. A9's 200 runs execute sequentially; its history count N is frozen from a throughput
smoke run of 10^3 histories, a configuration that is not A9's, before any A9 result exists.

## Evidence classes and scope

- A1 to A6, A10, A14: theory and conservation (deterministic bounds, never z-scores).
- A7, A9: falsification and uncertainty (statistical, LV).
- A11, A15: backend parity and self-consistency.
- A12: contract; A16: regression.
- A8: consistency (offline recomputation); A13: related model, **exploratory and non-gating**.
- No row is independent-Monte-Carlo evidence in V3-004; the gating independent comparison of LET
  is TOPAS ProtonLET (V3-010). Nothing here validates a biological model: the only lookup fixture is
  synthetic (`tests/data/synthetic_lookup.json`).

## Evidence archives (which archive satisfies which row)

| Rows | Archive (suite, qualification base 20361004) |
|---|---|
| A1-A6, A10, A12, A14, A16 (CI part), A15 (CPU chunks, 1 vs 3 workers), A11 CI | CI (GitHub), fixed small seeds; also the `pytest-scoring-warp-cpu` step of `lv4` |
| A16 against a524f209, A11-LV, A15 (LV), A7, A8, A9, A13 (non-gating) | `lv4` archive (the 1 vs 3 workers part only in a standard-mode run, deferred in single-process mode) |
| A11-HR, A15-CUDA, CUDA tests | `hr4` archive |

Rehearsal runs (base 20351004) are non-qualification evidence only.

## Amendments

¹ **2026-10-07, task V3-004, after the first A4b result at 9f52057 was observed (clarification of
rows A4 and A4b, no tolerance changed).** The frozen rows did not define `E(z)`. The first
implementation used the table CSDA energy `Rinv(R(E0) - rho z)` and A4b then failed (worst
5.3e-4 / 1.04e-3 / 1.32e-3 for a maximum step of 1 / 0.5 / 0.25 mm against 1e-4). The cause is the
trapezoid bias of the range table (`R(E)` is built by the trapezoid rule in `ln E` of `E/S` on the
200-points-per-decade grid and overestimates every increment by about 3.5e-5 relative, about 5.1 um
of the 158.6 mm range at 150 MeV; the exact analytic quadrature of the log-log interpolated `S`
gives `R(150 MeV)` = 158.625 mm against 158.630 mm tabulated), not the scoring: the engine's
short-step branch `S(E_mid) t` agrees with exact quadrature to about 3e-6, the range-inversion
branch follows the biased table, and the branch mix makes the transported `E(z)` differ from the
single-shot table inversion (-0.15 % at 1 mm, -0.41 % at 0.25 mm steps at z = 158.04 mm). Test
T2 of V3-003 (end depth within 1e-4 R, about 16 um) cannot see a 5 um bias. Since `LET_t` is the
path average of `S` of the transported particle, `LET_t dz` must equal that particle's own energy
loss. The definition is therefore clarified: `E(z)` is the **transported** particle's energy,
`E(z) = Rinv(R(E_b) - rho (z - z_b))` anchored at the last engine step boundary `(z_b, E_b)` at or
before `z`, taken from the trace of the deterministic history (`trace_histories = 1`). Under this
definition A4b gives 5.0e-5 / 4.6e-5 / 3.7e-5 for 1 / 0.5 / 0.25 mm and A4 3.8e-5 (3.7e-5 against
the table energy). The comparison with the table CSDA energy stays as a non-gating logged record.
The inconsistency of the range table is recorded as contrary evidence against the qualified
transport tables and assigned to a follow-up task, V3-003D (range table by exact quadrature of the
interpolated `S` and/or a range-carrying step); it moves the qualified T1/T9/T12 baselines by about
5 um, so it is a separate re-qualification and out of the scope of V3-004.

² **2026-10-07, task V3-004, after the first reduced A7 rehearsal at base 20351004 was observed: z = 215
at s_max = 1 mm, 3.5 / 3.9 at 0.1 / 0.5 mm, with the expected-share systematic of the plan delta
item 6 omitted from the frozen row (amendment of row A7).** The frozen row compared `LET_d` and
`LET_d^ε` within 4 sigma as if both had the same voxel expectation. They do not: the straggled
end-energy ramp shifts the expected deposit of each piece relative to the mean-loss apportionment
by `l_p g kappa (t_mid/s - 1/2) / S` (`g = -dS/dE`, `kappa = Var(eps)/s` the Bohr variance per path
length), a relative shift of at most `r = g kappa / (2 S^2)` (about 0.85e-4 to 1.9e-4 for 150 MeV
primaries). It is antisymmetric within a step, sums to zero over a whole step and survives only in
voxels that hold step starts or ends without their partners, which is the case for 1 mm steps on a
1 mm grid, where sigma of the difference is about 1e-6 relative. The criterion is amended to
`|D| <= 4 sigma_D + delta_v` with `D = LS2/LS - ES/E_step` per voxel, `sigma_D` the delta-method
standard error of `D` (batch covariance of the four channels) and `delta_v = r_v LET_d,v` evaluated
per voxel at `S_v = LET_d,v` (`E_v` from the inverse water table, `g` from its numerical derivative,
`kappa` from the engine's Bohr variance); no free constant. Derivation in the docstring of
`validation/scripts/transport/steps_v4.py`. If the observed difference exceeds the analytic bound
the amendment is not widened; the observation is reported as contrary evidence.

³ **2026-10-07, task V3-004, before any A9 result existed (calibration of row A9).** The window
[0.92, 0.98] around the nominal coverage 0.95 has a false-failure probability of about 5 % per depth
at 200 seeds (the sampling sd of a coverage fraction is 1.5 %), about 14 % over three depths. It is
replaced by the 99 % binomial interval for p = 0.95 at n = 200, [0.91, 0.99] (about 1 % per depth,
about 3 % over three depths). History count per seed (1e5 in 100 batches), 200 seeds and the
leave-one-out pooled reference are unchanged; A9 runs in two parts of 100 seeds
(`a9-part-1of2`, `a9-part-2of2`) merged by `a9-compare`.
