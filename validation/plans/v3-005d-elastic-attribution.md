# V3-005D elastic-model attribution plan (planner output; not evidence)

To be committed as `validation/plans/v3-005d-elastic-attribution.md` together with Amendment 19 (§6). The scratch scripts that produced already-seen numbers (`validation/plans/v3-005d-elastic-attribution/scripts/shape_moments.py`, the reviewer's `validation/plans/v3-005d-elastic-attribution/review-scripts/removal.py`) are committed alongside in `validation/plans/v3-005d-elastic-attribution/`, which is outside the rehearsal digest.

## Revision 4 (2026-10-14)

- Orchestrator edit (2026-10-14, after the route review): Amendment 19 (f) and the closing line name the re-qualification bases of record lv 20531004 / hr 20541004 (route document §4.4; seed test `tests/ionmc/test_v3_005d_seed_families.py`); 20511004 / 20521004 are retired unused.

This revision is text only and recomputes nothing. It answers round 3 of `validation/plans/v3-005d-elastic-attribution/review.md` (ACCEPT, conditional on these edits).
- **R3-1 (M2 buildable blind).** §4.1 now fixes, before any M2 shape is computed:
  - the NN input sources: σ_NN and α_NN from the GWU SAID NN partial-wave solution; β_NN from the Horiuchi et al. 2007 parameterisation;
  - an explicit rule below the lowest tabulated energy;
  - the de Vries 1987 density form per target for all seven targets, by a fixed preference rule;
  - the proton-size unfolding with r_p = 0.8414 fm (CODATA 2018);
  - the energy and frame conventions;
  - the lineage statements and the acquisition paths.
  If an acquisition cannot be made, M2 counts as "construction impossible", which is the declared M3 trigger.
- **R3-2.** Amendment 19 (c) now states the primary count once.
- **R3-3.** The near-identity optical-theorem trigger for M3 is replaced by input-sanity triggers that can catch wrong inputs:
  - M2's own unnormalised σ_el against the shared Barashenkov level, within a stated band;
  - M2's first-minimum position against the de Vries point-density Helm zero.
- **Confirmed numbers kept (reviewer, round 3).** The pinned table's p+O first minima are 35.3° / 29.2° CM at 135 / 179.9 MeV, against the measured 17.17° / 14.99°. A sharp disk would need R ≈ 4.9 fm. The rev-3 planner estimate of ≈ 34° / 29° is superseded by these values.


## Revision 3 (2026-10-14)

This revision answers round 2 of `validation/plans/v3-005d-elastic-attribution/review.md` (REVISE, minor; no blocking items).

- **R2-1 (pre-registration by reference).** Amendment 19 now cites the committed plan file `validation/plans/v3-005d-elastic-attribution.md` by path and sha256 (filled in at commit, in the same commit, before D-E0), not the scratch path.
- **R2-2 (H2a/H2b not identifiable).** The reviewer showed that exceedance and median ratios do not separate surface from radius (s = 0.9 fm: 1.22 / 1.24; R × 1.2: 1.17 / 1.19). Attribution is now made only at the level **"H2 (p+A shape)"**. The core/surface split is a descriptive, report-only characterisation of the Geant4 thin-target CDF (§2.2, §5 rule 4). The peak-versus-gain ordering needs a 2-paired-SE margin, otherwise it is "indeterminate". §4 no longer selects on any H2a/H2b label.
- **R2-3 (seed guard in hashed code).** Disjointness of the 2099xxxx family is asserted by a test under `tests/` only (outside the rehearsal digest). The `run_suite.py` seed-guard edit (l.1303–1313, hashed) is deferred to the re-qualification commit, which changes hashed code anyway. The D-E runs never pass through `run_suite.py`.
- **R2-4 (M1 inconsistent; selection informed by Geant4).** M1 is withdrawn as a candidate. It mixed the 2-D disk J₁ with the 3-D Helm radius and used a charge form factor as a hadronic profile. The decision is pre-committed to **M2 (optical-limit Glauber)**, with M3 = (ii′) as the fallback only on a pre-declared construction or verification failure. Selection never uses D-E5 (the Geant4 comparison). The already-seen V10-A first minima (p+16O: 17.17° at 135 MeV, C0063; 14.99° at 179.9 MeV, C0141), against the black-disk first zero of ≈ 34° / 29° CM for R_O = 2.68 fm (planner estimate), point to a core about 2× too wide. D-E1 must verify that these are true diffraction minima, not Coulomb-nuclear interference dips or single-point fluctuations. The reviewer's "p+C 42.0° / 35.9°" compares the C table with O datasets; the O values above are the like-for-like ones, and D-E1 recomputes both.
- **R2-5 (H1 cannot be established).** H1 now has an attribution path: a conditional `swap-sigma` / `swap-all` pair in D-E6, run only if the H1 "identified" signature fires. Without it H1 is reported as "not attributable by the §5 rule".
- **R2-6 (K* form).** The peak predictor is now the archived EM-only ionmc IDD used as the kernel. The Gaussian overlap form, corrected to exp(−Δ²/(2σ_R²)), is report-only.
- **R2-7 (plural scattering in D-E4).** Per-energy nσL values are stated. The target is 2 mm at 30 and 50 MeV. The off-curve fraction is reported, and the same KE cut applies to the D-E5 sampling. The finite-sample CDF error is stated.
- **R2-8 (cited removal-kernel numbers).** They are labelled crude estimates; D-E0 recomputes them and differences are expected. The script is archived with the plan.
- **Minors.** `IONMC_CACHE_DIR` and the complete scratch cache copy (§3.1).

## Revision 2 (2026-10-14)

This revision answers the independent review `validation/plans/v3-005d-elastic-attribution/review.md` (verdict REVISE). Changes by finding:

- **F1 (blocking, governing moment).** Peak depression is now modelled as a threshold-removal functional: the depth deficit R_res(1 − cos θ_lab) against the Bragg-peak width. Plateau gain is modelled as a path-length and local-deposit moment. Both are new deterministic quantities K_w and G (§2.1), computed in D-E0. The priors were re-derived. H2 is split into **H2a (tail / surface diffuseness)** and **H2b (core width / radius)**, with opposite peak-versus-gain signatures (§2.2). §5 rule 4 now discriminates them through the threshold-angle exceedance and the tail moment of the Geant4 thin-target distributions. The rev-1 claim that the 17° → 26.7° θ_rms difference explains the expectation miss is withdrawn; it is now a D-E0 question.
- **F2 (blocking, circular rule).** The core test is now a pre-declared **swap test** (D-E6): ionmc's p+A and p-p angular distributions are replaced by the empirical Geant4 thin-target distributions of D-E4 at unchanged σ, in a 2 × 2 design with an interaction term. The ionmc-only parametric variants are kept as sensitivity probes and can never "explain" Δ. Every hypothesis, H2b included, has a computable "excluded" signature (§2.2, §5).
- **F3 (major, D-E4 design).** Classification is now per track, by elastic two-body kinematics in an EM-free target scored on an enclosing vacuum sphere. Creator-process tagging is gone. Rows 16/17 (EM on, Coulomb confound) are dropped. The case matrix is fixed and covers water and graphite at 30–200 MeV. The OpenTOPAS items that must be verified before the cases are frozen are listed (§3, D-E4).
- **F4 (major, seeds).** One seed per energy is shared by EM-only and every variant (common random numbers). The family starts at 20990100. The scratch use of numpy seed 20990001 (candidates2.py, an analytic calculation, not transport) is recorded.
- **F5 (major, digest).** The rev-1 statement "no `src/` change" with the E-shape precedent was wrong: E-shape changed `src/ionmc/nuclear/elastic_build.py`, and `validation/scripts/transport/` is inside the rehearsal digest. New route: a standalone post-processing script under `validation/scripts/elastic_attribution/`, outside `SOURCE_PREFIXES`. It writes exploratory tables that are loaded through the existing `PhysicsConfig(elastic_table_id=…, allow_exploratory_table=True)` (§3.1).
- **F6 (major, pre-registration).** Amendment 19 (§6) now freezes the variants, parameters, backend and history counts, the D-E4 matrix and primaries, the classification and threshold angles, the D-E5/D-E6 sizes and seeds, the kernel definitions and the §5 rule. It is committed **before D-E0 runs**. The items the planner has already seen are declared.
- **F7 (major, §4).** The selection rule is declared data-seen (V10-A). M1's density source and R rule are fixed now (§4).
- **Minor (F8–F10).** The thin-target p-p comparison applies the 16.26° CM cut. V10-A data contain Coulomb and interference. The V5-versus-V11 equality is now stated with its uncertainty. The expectation source line is cited. 1e5 warp-cpu histories are required (the python fallback is withdrawn). Numeric and path errors are corrected: the V11 F differences are **relative**, −0.039 / −0.049; the **absolute** differences are −0.035 / −0.043. Step paths are fixed (step 03 in RUN-20261012T142343Z-081123c3; step 16 in RUN-20261012T142957Z-a8e41216; step 21 in RUN-20261012T143203Z-1b3b1ecf; step 22 in RUN-20261012T143235Z-2b818275).

Scope and frame are unchanged:
- Scope: explain why the hadronic elastic channel at the evidence SHA dd976c5a is about 30–40 % stronger than TOPAS `g4h-elastic_HP` (V11).
- Nothing here changes a frozen row, tolerance, base or verdict. Everything is **exploratory / report-only** until a re-qualification is pre-registered.
- Binding frame: operator disposition A (IR-20261013-164616-DEFAF4); Amendment 14 (a) (interpretation-B conditions), (e)9 and (f).

## 0. Executive summary

- **The residual is elastic.**
  - The V5 TOPAS peak/plateau residual is −3.66 % ± 0.41 % / −4.89 % (step 21; step 22 `pp_diff` −0.0366 / −0.0489).
  - The V11 relative difference F_ionmc/F_TOPAS − 1 is −3.90 % ± 0.29 % / −4.88 % ± 0.38 % (step 16).
  - These agree within one combined SE at both energies, so the non-elastic part is consistent with TOPAS (F_ne ≈ 1).
  - The elastic depression 1 − F is 0.143 / 0.156 in ionmc against 0.108 / 0.113 in TOPAS. Δ(1 − F) = +0.035 / +0.043, a ratio of 1.32 / 1.38.
  - The plateau gain is 0.0221 / 0.0267 against 0.0163 / 0.0184. Δgain = +0.0058 / +0.0083, a ratio of 1.36 / 1.45.
- **Leading hypothesis: the p+A angular shape (H2).** The core/radius and tail/surface parts are not separable by the declared signatures and are characterised descriptively only. Already-seen archived evidence points to the core: the measured p+16O first minima are 17.17° (135 MeV) and 14.99° (179.9 MeV), against the pinned table's first minima of 35.3° / 29.2° CM (reviewer-confirmed; a sharp disk would need R ≈ 4.9 fm). D-E1 verifies.
  - Data side: V10-A (step 03) has the model above the data by +44 … +81 % at 15° and +79 … +191 % at 20° for every ≥ 150 MeV dataset, with partial integrals of +20 … +100 %.
  - Built-table side: p+O at 150 MeV has θ_rms,CM = 26.7° against the decision-0041 planning value of 17°.
  - Peak-removal kernel (the reviewer's **crude** estimate `validation/plans/v3-005d-elastic-attribution/review-scripts/removal.py`: Bragg–Kleeman range, straight line, p+O only, step threshold, no straggling, no p-p; D-E0 recomputes it with the project range table and the IDD kernel, and differences are expected): removal at w = 2 / 3 mm changes by −10 / −15 % with a 0.5 fm surface factor, −24 / −33 % with 0.9 fm, and −19 / −24 % with R × 1.2. The Σ⟨1 − cos θ⟩ moment changes by −40 %, −60 % and −18 % respectively.
  - So the peak is about as sensitive to R as to the surface, and the exceedance and median ratios move alike for both (R2-2). The split is therefore reported, not attributed.
- **Smaller contributors a priori.**
  - H1 (σ scale): shared with TOPAS, and V10-A (c) is within −3.3 … +5.5 % at 6 of 7 C points.
  - H3 (p-p): P(pp) is about 0.037 per 150 MeV primary. Most of it falls at 12.5–30 MeV, which is irrelevant to the peak. σ_pp is +10 % above BGG only at ≥ 150 MeV.
  - H4 (recoil): 0.249 / 0.300 MeV per primary, 0.19 / 0.13 % of the plateau, at most about 0.002 of Δgain; coupled to H2 through ⟨1 − μ⟩.
  - H5 (MCS overlap): closed on provenance.
- **Plan.**
  1. Amendment 19 (pre-registration) first.
  2. Deterministic sandbox steps D-E0 to D-E2.
  3. Sandbox warp-cpu sensitivity runs D-E3 (CRN, 1e5 histories).
  4. One reference-engine set D-E4: 10 TOPAS thin-target cases, EM-free, ≤ 1500 s each, one job at a time.
  5. D-E5 table-versus-Geant4 distribution comparison.
  6. D-E6 swap test (decisive).
  7. Synthesis by the §5 rule.
  No lv5c/hr5c step and no fresh base is used.
- **Operator proposals needed (§7):**
  - the status of V10-A after any post-observation shape change (the (e)9 precedent makes it exploratory);
  - V10-A (c), which fails on the shared Barashenkov level, independent of the shape.
  Everything else is autonomous.

## 1. Evidence base (archived; values verified by the planner, the reviewer reproduced them)

| Item | Source | Value |
|---|---|---|
| V11 F | `validation/generated/RUN-20261012T142957Z-a8e41216/transport/lv5c-dd976c5-16/16-v11-compare.txt` | ionmc 0.8568 ± 0.0015 / 0.8441 ± 0.0025; TOPAS 0.8916 ± 0.0022 / 0.8874 ± 0.0024. Relative diff −0.0390, 90 % [−0.0446, −0.0334] / −0.0488 [−0.0559, −0.0417]. Absolute −0.0348 / −0.0433 |
| V11 plateau gain | same | ionmc 0.0221 ± 0.0007 / 0.0267 ± 0.0008; TOPAS 0.0163 ± 0.0004 / 0.0184 ± 0.0001 |
| TOPAS rows 12/13 (report-only) | same, `elastic_factor_rows_12_13` | F_el 0.8891 / 0.8906; plateau ratio 1.0155 / 1.0199 |
| V5 vs TOPAS | `RUN-20261012T143203Z-1b3b1ecf/transport/lv5c-dd976c5-21/21-v5-compare.txt` | peak/plateau −3.66 % / −4.89 %; ΔIDD plateau +10.35 % (fail) / +6.72 %; plateau +0.98 % / +0.62 %; R80 −0.30 / −0.32 mm; total deposit −0.30 % / −0.49 %. MCsquare (report-only) −4.78 % / −5.34 % |
| Attribution step | `RUN-20261012T143235Z-2b818275/transport/lv5c-dd976c5-22/22-v5-attribution.txt` | rung 0 pass (P7, R1-D2, V8-LV; P5-ext/P6/P6-D CI); stop rule at rung 1 |
| V10-A | `RUN-20261012T142343Z-081123c3/transport/lv5c-dd976c5-03/03-v10-a.txt`, fields `result.datasets[].a[].ratio_minus_1`, `.b.ratio_minus_1`, `result.level_points.gating[]` | (a) at 15°/20°: C0148@201.4 +0.81 / +1.91; C0550@200 +0.50 / +1.21; O0553@150 +0.55 / +1.70; C0055@160 +0.44 / +1.67; C0057@200 +0.62 / +1.81. (b): C0063@135 +0.30, C0141@179.9 +0.20, C0550 +0.66, C0055@122 +0.88, @160 +1.00, C0057 +0.86. (c): C at 81 / 100 / 119 / 141 / 158 / 180 / 65.5 MeV: +0.055 / **+0.171** / +0.020 / +0.002 / −0.025 / −0.033 / +0.030 |
| Table moments | `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff….npz` (`edges_mu`) | p+O ⟨1 − μ_CM⟩ 0.118 / 0.093 / 0.082, θ_rms 29.9 / 26.7 / 25.5° (100 / 150 / 200 MeV); p+C ⟨1 − μ_CM⟩ 0.157 / 0.115 / 0.098; σ_pp 27.3 / 26.3 / 25.6 mb (BGG formula 27.6 / 23.9 / 23.3) |
| Recoil | `.ionmc-cache/transport-samples/dd976c5/{17,19}/samples/v5-*-on.json` | 0.249 / 0.300 MeV per primary; plateau share 0.19 % / 0.13 % |
| E-shape | step 85b | R80 shift 0.0004 mm |
| Expectation | plan l.255–258 (Amendment 14 (f)); decision 0041 l.272 (planning θ_rms 21 / 17 / 14°) | 0.985 / 0.986 expected; 0.963 / 0.951 observed |

Already seen by the planner before this pre-registration, and declared: the V10-A per-point ratios, the table moments, the shape moments of `validation/plans/v3-005d-elastic-attribution/scripts/shape_moments.py` (s = 0.5 / 0.9 fm, R = 3.0 / 3.3 fm, Gaussian), the reviewer's removal-kernel numbers (`validation/plans/v3-005d-elastic-attribution/review-scripts/removal.py`), and the P(pp) estimate. Rev-1 sources are `validation/plans/v3-005d-elastic-attribution/scripts/shape_moments.py`, an inline numpy estimate (P(pp) 0.0375 ionmc / 0.0394 BGG, transfer 0.73 / 0.76 MeV; Bragg–Kleeman α = 0.0022, p = 1.77, slow proton ≈ T/4, not scripted), `validation/plans/v3-005d-elastic-attribution/scripts/binom_gate.py` and `candidates*.py` (V7-R). D-E0 recomputes every planning number in a committed script.

## 2. Model of the observables and hypotheses

### 2.1 Deterministic kernels (D-E0; definitions frozen in Amendment 19)

The straight-line CSDA path in water uses the project range table (not Bragg–Kleeman), σ_el(E) and the μ_CM distribution of the variant, with exact two-body CM→lab kinematics (target masses from the table).

- **Depth deficit** of an event at energy E: Δ(E, μ) = R(E) − R(E − T_rec(μ))·cos θ_lab(μ). T_rec is the recoil energy. For p-p the slower proton is a secondary, and the primary keeps the faster one, Δ_pp = R(E) − R(E_fast)·cos θ_fast.
- **Peak-removal kernel** K_w = Σ_targets ∫₀^{0.98R₀} n_t σ_t(E(z)) P_μ(Δ > w) dz, for w = 1, 2, 3, 5 mm. **Peak predictor K*** (rev 3): an event with deficit Δ moves that proton's contribution from the IDD shape at z_pk to z_pk − Δ. The archived EM-only ionmc IDD (V11 EM-only partials `v11-150-emonly.json`, `v11-200-emonly.json`) is used as the kernel, so K* = ∫ n σ ⟨1 − IDD_EM(z_pk − Δ)/IDD_EM(z_pk)⟩ dz. This is deterministic, needs no new run, and keeps the asymmetric Bragg peak. The Gaussian form ⟨1 − exp(−Δ²/(2σ_R²))⟩ (σ_R 1.9 / 3.1 mm) is report-only. K* only validates D-E3 and never enters §5.
- **Plateau-gain kernel** G = ∫ n σ ⟨(1/cos θ_lab − 1)·f_path + T_rec,local/E⟩ dz, evaluated over z in the plateau window 20–60 mm, plus the p-p secondary-proton deposit there.
- **Threshold angles** θ_w(E): R(E)(1 − cos θ_w) = w, for w = 2 and 3 mm. These are the angles that set K.
- **Tail moment** T(E) = ⟨(1 − cos θ_lab)·1[θ_lab > 3·θ_med(E)]⟩, where θ_med is the median deflection.

### 2.2 Hypotheses, priors and computable signatures

Notation: Δ_pk = Δ(1 − F) (observed +0.035 / +0.043), Δ_g = Δgain (+0.0058 / +0.0083). "Ratio" means ionmc/Geant4 from D-E5. "Swap" means D-E6.

| H | Statement | Prior (rev 2) | Discriminating measurement | Identified (main cause) | Excluded |
|---|---|---|---|---|---|
| H1 σ_el(p+A) scale | the ionmc Barashenkov transcription differs from Geant4's runtime BGG σ_el | low (shared lineage; P6 (5); V10-A (c) 6 of 7 within 6 %) | D-E5 rate ratio per element, from the D-E4 event count / (N·n·L), with σ_stat ≈ 1–2 %; conditional `swap-sigma` | rate ratio for O outside 1 ± 0.10 at ≥ 2 energies in the direction of Δ_pk, **and** (if run) `swap-sigma` closes ≥ 50 % of Δ_pk; without `swap-sigma`, H1 is "not attributable by the §5 rule" | rate ratio within 1 ± 0.05 at all five energies for O and C |
| H2 p+A angular shape | the black-disk shape (D-1 radius rule, sharp edge) puts too many events beyond the threshold angles θ_w and into the large-angle tail | **high** (V10-A shape failure at 15°/20°; first minima ≈ 2× too wide, to be verified in D-E1; table θ_rms 26.7° against the 17° planning value) | D-E5 ratios of P(θ_lab > θ_2mm), P(θ_lab > θ_3mm), θ_med and T; swap `pA` | swap-pA closes ≥ 50 % of Δ_pk **and** at least one of the exceedance ratios ≥ 1.15 at ≥ 3 of 5 energies or the T ratio ≥ 1.5 | exceedance, θ_med and T ratios all within 1 ± 0.10 at all energies, **and** \|swap-pA effect on F\| < 0.005 |
| H3 p-p channel | the Hale NI (θ_CM ≥ 16.26°) rate or angular / transfer distribution differs from Geant4 BGG + CHIPS p-p | low | D-E5 p-p rate above 16.26° CM (cut applied to Geant4), ⟨T_slow⟩, θ_fast exceedance; swap `pp` | rate or ⟨T_slow⟩ ratio outside 1 ± 0.10, **and** swap-pp closes ≥ 30 % of Δ_pk | ratios within 1 ± 0.05 **and** \|swap-pp effect on F\| < 0.005 at both energies |
| H4 energy bookkeeping | the locally deposited recoil, or the e_min domains, inflate the plateau or shift the peak | low | recoil channel per run (D-E3/D-E6); E-shape (done) | recoil plateau share change in swap-pA ≥ 50 % of Δ_g | recoil-explained part < 20 % of Δ_g, and E-shape R80 shift < 0.05 mm (observed 0.0004 mm) |
| H5 MCS overlap | the EM MCS contains nuclear elastic content, or the TOPAS EM/hElastic partition differs | very low (f_dM theory-calibrated; Amendment 15 (d)2) | EM-only peak/plateau ionmc against TOPAS from the archived V11 EM-only inputs (seconds) | EM-only difference > 2 % **and** in the direction that amplifies Δ_pk | EM-only peak/plateau within ±2 % (F is a ratio, so EM cancels to first order) |
| H6 comparator definition | the TOPAS rows 10/11 and V5 EM-only cases differ in more than the physics list | ≈ 0 (rows 12/13 agree with F) | D-E0 case-file diff | any non-physics difference | identical apart from the modules and seeds |

Additivity: Δ_pk = sum of the main effects + the interaction from the 2 × 2 swap (§5). It is not assumed.

**Descriptive core/surface characterisation (report-only; never an attribution or selection label).** For the Geant4 thin-target CDF and the ionmc table at each D-E4 energy:
- the core width: θ_med and the first-minimum estimator (the first local minimum of the binned dσ/dΩ_CM above 3°, with bin width fixed at 1° CM);
- the surface: the large-q log-slope of dσ/dq² over q above the first minimum;
- a report-only least-squares (R_eff, s) fit of the 3-D Helm form |3j₁(qR)/(qR)·exp(−q²s²/2)|² to the Geant4 CDF. This is a fit to Geant4, not to measured data.

The peak-versus-gain ordering in the swap (relative change of 1 − F against relative change of the plateau gain) is reported as "core-like" or "tail-like" only if the two differ by more than 2 paired SE; otherwise it is "indeterminate".

## 3. Steps (sandbox first; one compute job at a time; never concurrent with a host job or another agent's numpy job)

| Step | Where / cost | Content | Output |
|---|---|---|---|
| D-E0 | sandbox, deterministic, minutes | K_w, K*, G, θ_w, T for the built table and the D-E3 parametric variants at 50 / 100 / 150 / 200 MeV, with the project range table; reconciliation of the 0.985 / 0.986 expectation against K* (which part of the miss comes from table moments versus planning moments); H6 case-file diff; H5 EM-only check from the archived V11 inputs | kernel table (committed script) |
| D-E1 | sandbox, deterministic | per V10-A dataset: model and data at the data angles, the data-integrated partial moment from 10° to θ_max of the data, the same for the §4 candidates at their **a-priori** parameters. Data-seen; never a fit used for a model parameter | exploratory table |
| D-E2 | sandbox, deterministic | p-p: σ_pp(E) (ionmc NI, BGG formula, PDG points), T_slow distribution, θ_fast exceedance at θ_w | audit table |
| D-E3 | sandbox, warp-cpu float64, single process; **1e5 histories per run required** (the sandbox rate is measured first; if > 120 s per run, the step is run one variant at a time) | V11 harness (EM+elastic against EM-only, V5 geometry r = 20 cm, 20 batches) at 150 and 200 MeV for the parametric probes `all`, `pp-only`, `pA-only`, `pA-sigma0.8`, `pA-Rx1.2`, `pA-s0.5`, `pA-s0.9`, `pA-gauss`, plus EM-only. One seed per energy shared by all 9 runs (CRN). Reported: F, gain, recoil, R80, lateral rms (R1-D2 (C)), K* prediction against the measured F. Probes only; **never** used in §5 identification | sensitivity table; the K* model is validated or refuted against the ionmc-only F changes |
| D-E4 | host reference engine (`run_reference_calculation`, OpenTOPAS 4.3 / Geant4 11.4.2); 10 cases × ≤ 1500 s; one job at a time | Thin-target single scattering, see §3.2 | Geant4 per-class rates and angle / energy distributions |
| D-E5 | sandbox, minutes | the ionmc table sampled (P7 sampler, 1e6 events per case) in the D-E4 frame (same energies, lab conversion, 16.26° CM cut on p-p for both); ratios for every §2.2 signature | ratio table |
| D-E6 | sandbox, warp-cpu float64, 1e5 histories, CRN seeds of D-E3 | **Swap test.** Exploratory tables where the μ_CM inverse CDF of O (the only p+A target in the V5 water geometry) is replaced by the Geant4 empirical CDF of the water rows 14–18 (`swap-pA`; the graphite rows 19–23 serve H1/H2 signatures for C and the V10-A C comparison only); the H-1 CDF by the Geant4 p-p CDF above 16.26° (`swap-pp`); both (`swap-both`); σ unchanged (ionmc's). Between D-E4 energies the CDF is interpolated in the scaled variable q = 2p_CM sin(θ_CM/2) (linear in E at fixed CDF level); below 30 MeV the ionmc shape is kept (declared; K* weight below 30 MeV computed in D-E0). Runs: base, swap-pA, swap-pp, swap-both, EM-only, at 150 and 200 MeV | Δ_pk and Δ_g main effects and interaction → §5 |

Order: Amendment 19 committed (§6) → D-E0 → D-E1 → D-E2 → D-E3 → D-E4 (requested only after D-E0 to D-E3 are reported) → D-E5 → D-E6 → synthesis.

### 3.1 Exploratory-table route (no `src/` change, no rehearsal-digest change)

- A standalone script `validation/scripts/elastic_attribution/make_variant_tables.py` sits outside `SOURCE_PREFIXES` (`src/ionmc`, `tests/ionmc`, `validation/scripts/transport`, `benchmarks/transport`; run_suite.py l.95–100), so neither the suite identity nor the V7-R rehearsal digest changes. It:
  - loads the pinned table `1fd24cff…`;
  - writes variant `edges_mu` (shape probes, swaps) or scaled / zeroed `sigma_barn` rows (pp-only, pA-only, σ × 0.8);
  - adds an `exploratory` key naming the variant and the swap source (D-E4 run ids);
  - recomputes the npz sha256 and the id with `ionmc.nuclear.build.table_id`;
  - writes into a scratch cache copy selected with `IONMC_CACHE_DIR` (cache.py l.37–48), because `ElasticTable.load` is called with `cache_dir=None` (config.py l.592). The copy must hold **every** derived input the run loads: the pinned nuclear table `00e8031d…`, the transport tables, and the pinned and variant elastic tables. A preflight lists the loaded ids and compares them with the expected set.
- Runs load the table through `PhysicsConfig(elastic_table_id=<variant id>, allow_exploratory_table=True)` (config.py l.113/140, l.592; elastic_tables.py l.138–142).
- The variants carry the copied qualification flags unearned. This is disclosed in every output. The script itself asserts normalisation, the non-negative density and σ = 0 below the domain, because the exploratory key disables the loader's below-domain assertion.
- Before D-E3, a one-run check verifies that the warp backend reads the identical arrays: the base-variant table reproduces the pinned table's F bit for bit at the same seed.
- The rev-1 E-shape precedent was a `src/` change (`exploratory_forced_isotropic`, elastic_build.py l.832) and is **not** the route used here.

### 3.2 D-E4 thin-target cases (rows 14–23)

- **Geometry:** a water cylinder (rows 14–18; G4_WATER) or graphite cylinder (rows 19–23; G4_C, density 2.0 g/cm³), radius 1 cm, length **1 cm at 100–200 MeV and 2 mm at 30 and 50 MeV** (rows 14, 15, 19, 20), in a G4_Galactic world. The scoring surface is the inner surface of a vacuum sphere of radius 50 cm centred on the target, so every exiting proton is scored at its exit angle and energy.
- **Physics:** module list `g4h-elastic_HP` **only**: no EM, no inelastic, no decay. Without EM a proton changes direction or energy only through hadronic elastic scattering. Undeflected primaries keep exactly E₀.
- **Beam:** pencil, monoenergetic. E₀ = 30 / 50 / 100 / 150 / 200 MeV (rows 14–18 and 19–23 in that order). 4e6 primaries per case at 100–200 MeV and 2e7 at 30 / 50 MeV, to compensate for the 2 mm target (wall time verified by the smoke case; if > 1500 s the case is split into k = 1, 2 seeds of equal size, declared now). Expected counts: p+O events at 150 MeV ≈ 2.1e4 and p-p ≈ 6.4e3; at 200 MeV ≈ 1.4e4 and 6.2e3 (basis n·σ·L·N with the BGG σ that Geant4 uses; Poisson SE 0.7–1.3 %).
- **Seeds:** 20270000 + 1000·row + 1 (one seed per case; distinct; the 20284001 … 20293001 range is unused in `validation/reference_cases`).
- **Scoring:** a phase-space scorer on the sphere surface, protons only, with a kinetic-energy filter that drops records with E > E₀ − 1e-4 MeV. Without EM these are the undeflected primaries; the filter keeps p+A events with θ_CM ≳ 0.35° at 150 MeV, which is irrelevant for K. Recorded per proton: event id, track id, parent id, E and direction. Expected output ≤ 5e4 records per case.
- **Classification, per proton track:**
  - primary (parent id 0): compute the residuals to the two elastic kinematic curves E′(θ_lab) for p+H (two-body, m_t = m_p) and p+target (O or C). Assign the class with the smaller residual. A "conflict" class (both residuals < 1e-3 relative, which happens only at θ_lab < 0.5°) is recorded and excluded.
  - secondary protons (parent ≠ 0, same event): the p-p partner. Its presence confirms p-p.
  - Plural scattering: nσL for O in 1 cm water is ≈ 0.012 / 0.0053 / 0.0034 at 100 / 150 / 200 MeV, and for C in 1 cm graphite ≈ 0.025 / 0.0125 / 0.0084. In 2 mm targets it is ≈ 0.007 / 0.005 (O at 30 / 50 MeV) and ≈ 0.018 / 0.012 (C). The plural fraction among scattered events (≈ nσL/2) is therefore ≤ 1 % everywhere. Tracks off both kinematic curves (relative residual > 1e-3) are counted and reported per case as the "off-curve fraction" and excluded. The same KE cut and the same exclusion are applied to the D-E5 sampling.
  - The finite Geant4 sample (≈ 1–2e4 events per case) is one fixed realisation of the swapped CDF. Its effect on F is estimated with K* by resampling the Geant4 events (200 bootstrap CDFs, deterministic) and is expected to be ≲ 1e-3. It is added in quadrature to SE_r in §5.
  - p-p is reported in the CM frame with the 16.26° cut applied (both protons, identical-particle convention of Amendment 14 addendum (4)).
- **Must be verified against the OpenTOPAS 4.3 documentation before the case files are frozen (C7 practice):**
  1. the phase-space scorer parameter names for a surface on a sphere and the inclusion of event id, track id and parent id;
  2. the name of the kinetic-energy filter parameter;
  3. that a module list without any EM module is accepted by TOPAS and by the `run_reference_calculation` engine fingerprint and module allowlist;
  4. the output format (ASCII or binary) and its size;
  5. wall time ≤ 1500 s per case at 4e6 primaries (a smoke case of 1e5 primaries first, timing recorded).
  If item 3 fails, the cases fall back to `g4em-standard_opt4` + `g4h-elastic_HP` with classification restricted to θ_lab > 5°, where the single-Coulomb tail is computed and subtracted. This is declared now; the Coulomb confound is then stated.

## 4. Candidate replacement shapes (for the later decision; pre-committed now)

Data-seen declaration: the planner has seen the V10-A per-point ratios and the measured first minima. Any later use of D-E1 is therefore **data-informed with respect to V10-A**. D-E5 (Geant4, the physics of the V11 reference) and V5 / V11 / D-E3 / D-E6 outcomes are **never** selection criteria. If D-E5 is quoted in the decision at all, its Geant4-informed lineage is recorded.

| Option | Construction (fixed now) | Lineage | Status |
|---|---|---|---|
| **M2 optical-limit Glauber** (pre-committed) | Elastic amplitude F(q) = ik ∫ b db J₀(qb)[1 − exp(iχ(b))], with χ(b) = ½ σ_NN(i + α_NN) ∫ T_A(b′) f_NN(b − b′) d²b′. T_A is the thickness function of the point-nucleon density, and f_NN(b) = exp(−b²/(2β_NN))/(2πβ_NN) is the Gaussian NN profile. **All inputs, forms, conventions and rules are fixed in §4.1**, and every source is acquired and pinned in the decision before any M2 shape is computed. The shape is normalised; σ_el stays Barashenkov (D-1 level unchanged). Core and surface follow together from the density folded with the NN range; there is no free radius rule | independent of the p+A elastic data and of TOPAS; NN lineage shared with BGG (stated) | adopted if §5 attributes H2 |
| M3 = (ii′) | ENDF/B-VIII.0 neutron O-16 / C-12 MF4/MT2 (new acquisition, MF1 copy check) | related (isospin mirror); its V10-A results are exploratory ((e)9) | fallback, **only** on a pre-declared M2 failure (§4.2). Never on V5, V11, V10-A, D-E5 or D-E6 |
| Helm 3-D form (diagnostic only) | |3j₁(qR)/(qR)·exp(−q²s²/2)|² | — | report-only in D-E1 and as the Geant4 characterisation of §2.2; never a candidate |
| M4 Geant4 CHIPS | — | shared with TOPAS | rejected for gating |
| (withdrawn) M1 rev 2 | 2-D J₁ with a Helm 3-D radius; charge form factor used as a hadronic profile | — | withdrawn as inconsistent (R2-4) |

### 4.1 M2 inputs and conventions (fixed now; values are taken from the sources at acquisition, none is set by the planner)

1. **NN inputs: sources and lineage.**
   - **σ_NN and α_NN** (the ratio of the real to the imaginary forward amplitude): from the forward pp and np amplitudes and the total cross sections of the current **GWU SAID NN partial-wave solution** (Workman, Briscoe, Strakovsky et al.; https://gwdac.phys.gwu.edu). The solution id and date are recorded at acquisition.
     - Acquisition: `acquire_reference_data`, role `construction`, licence or terms recorded. The fallback is the published SP07 solution tables (Arndt, Briscoe, Strakovsky, Workman, Phys. Rev. C 76 (2007) 025209) if the live service is unreachable. If neither is reachable, M2 counts as "construction impossible" (§4.2).
     - Isospin average for the targets: σ_NN = (Z·σ_pp + N·σ_np)/A, and α_NN is weighted the same way by Im f(0) = kσ/4π.
     - Lineage: SAID fits the world NN data, which includes the pre-1998 pp data of Hale's fit (transport p-p channel). This shared lineage is stated; the p-p channel is already report-only (V10).
   - **β_NN** (the slope of the Gaussian NN profile, in fm²): from **Horiuchi, Suzuki, Abu-Ibrahim, Kohama, Phys. Rev. C 75 (2007) 044607**, the table of NN profile parameters (pp and pn), isospin-averaged as above.
     - Acquisition as above, role `construction`.
     - The energy coverage and the definition of β in that paper (profile normalisation, fm²) are verified at acquisition and recorded.
     - If it is unavailable, Ray, Phys. Rev. C 20 (1979) 1857 is used at ≥ 100 MeV only, with the low-energy rule below (recorded).
     - Lineage: NN data only, independent of the p+A elastic data and of TOPAS.
2. **Rule outside the tabulated range (the model runs at 6.6–250 MeV).**
   - σ_NN and α_NN: SAID covers 0–3 GeV, so they are evaluated at every grid node; no extrapolation.
   - β_NN: linear interpolation in T_lab between tabulated energies. Below the lowest tabulated energy E_β,min, β_NN is **held constant** at β_NN(E_β,min). Above the highest it is held at its last value (not expected within 250 MeV).
   - Rationale, stated without a claim of accuracy: below about 100 MeV the optical limit is an approximation in any case (§4.3). Holding β_NN constant is the least-assumption rule.
   - The fraction of the K* weight below E_β,min is reported from D-E0.
3. **Densities: the de Vries, De Jager & De Vries 1987 form per target** (At. Data Nucl. Data Tables 36 (1987) 495), acquired as above, role `construction`.
   - Preference rule: the **harmonic-oscillator (HO)** form for A ≤ 16 (C-12, N-14, O-16), so the thickness function T_A is analytic; the **two-parameter Fermi (2pF)** form for A ≥ 27 (Al-27, Si-28, P-31, Ca-40).
   - If the preferred form is not tabulated for a target, the next in the order 3pF → Fourier–Bessel → sum of Gaussians is used. The form actually used per target is recorded in the decision.
   - The surrogate mapping of decision 0041 (a target without its own entry uses the nearest tabulated nucleus, radius scaled by (A/A_ref)^{1/3}) applies unchanged.
4. **Point-nucleon density.**
   - The charge density is unfolded by the proton size in form-factor space: F_pt(q) = F_ch(q) / exp(−q² r_p²/6). The rms radii satisfy ⟨r²⟩_pt = ⟨r²⟩_ch − r_p².
   - r_p = **0.8414 fm** (CODATA 2018 recommended proton rms charge radius; recorded with its source). The neutron charge radius and Darwin–Foldy terms are neglected (declared).
   - ρ_n = (N/Z)·ρ_p with the same shape (exact symmetry assumed for N = Z; a declared approximation for Al-27 and P-31).
5. **Energy and frame conventions.**
   - The NN inputs are evaluated at the NN laboratory kinetic energy equal to the incident proton's T_lab (the projectile on a free nucleon at rest). σ_NN and α_NN are invariant, and β_NN is defined in the NN impact-parameter plane, which is taken as the p+A impact-parameter plane (standard optical-limit convention).
   - The p+A amplitude uses k = p_CM of the p+A system (the table's `p_cm_mev`), q = 2k sin(θ_CM/2), with F normalised so that dσ/dΩ_CM = |F|².
   - The tabulated shape is |F|² normalised to unit integral over μ_CM, then scaled to the Barashenkov σ_el (the D-1 level is unchanged). Coulomb stays in the MCS model; there is no Coulomb–nuclear interference.

### 4.2 M3 fallback triggers (replace the rev-3 optical-theorem trigger)

M3 is used only if one of these holds. Each is evaluated once, on the M2 construction, before any transport with M2.
1. **Construction impossible.** A §4.1 source cannot be acquired or does not contain the required quantity.
2. **Numerical verification fails.** The normalisation error exceeds 1e-6, or the b- and q-integration convergence (step halving) exceeds 1e-3 relative.
3. **Level sanity.** M2's own unnormalised σ_el^M2(E) = ∫|F|² dΩ differs from the shared Barashenkov σ_el by more than a factor outside **[0.7, 1.4]** at ≥ 2 of 100 / 150 / 200 MeV for O or for C.
   - The band is a planner bound, not an accuracy claim: optical-limit Glauber is expected to reproduce p+A σ_el at this energy to well within it.
   - A larger deviation indicates wrong inputs, such as units, the α/β convention, or the density normalisation.
4. **Radius sanity.** M2's first-minimum momentum transfer q_min,M2 at 150 MeV (O and C) must lie within **[0.6, 1.0] × q_H**. Here q_H = 4.493/R_H is the first zero of the 3-D Helm form factor of the de Vries point density, with R_H² = (5/3)(⟨r²⟩_pt − 3s²) and s taken from the density's surface width.
   - The NN-range folding moves the minimum to smaller q but not below about 0.6 q_H (planner bound).
   - A value outside the interval indicates wrong inputs.

Neither trigger uses V10-A, V5, V11, D-E5 or D-E6. The measured minima are reported in D-E1 (data-seen), never used to switch models.

### 4.3 Limitations stated now
Known limitation stated now: the optical-limit Glauber approximation weakens below about 100 MeV (Pauli blocking, in-medium NN amplitudes). D-E1 reports M2 against the 75 / 96 MeV C datasets, and the decision states the domain caveat. The black-disk is not kept below 100 MeV, so that the model stays a single construction.

D-E1 (data side, declared data-seen) reports for the black disk and M2: the first-minimum position against the measured one, after **verifying that each measured "first minimum" is a diffraction minimum**. Criteria: the dip is present in ≥ 2 adjacent points; it is deeper than 2 σ_stat; and it is not explained by the Rutherford-nuclear interference estimate (the Rutherford share of Amendment 14 (e)4 at that angle < 20 %). The (a) / (b) quantities are also reported. These are reported, not selection criteria. M2 is pre-committed.

## 5. Attribution established (pre-declared, computable)

Let Δ_pk(X) be the change of (1 − F) between `base` and swap X (CRN, so the SE comes from the paired batches), and likewise Δ_g(X). Then:
- the interaction is I = Δ(both) − Δ(pA) − Δ(pp);
- closure is the residual r = Δ_pk,obs − Δ(both), where Δ_pk,obs is the V11 value at dd976c5a, with its SE combined in quadrature with the swap SE;
- SE_obs = 0.0026 / 0.0034 absolute (the step-16 relative-difference SE 0.0029 / 0.0038 times F_TOPAS; recomputed in D-E0).

Attribution is **established** iff, at **both** energies:
1. `base` reproduces the archived ionmc F within 2 combined SE, and the warp table-identity check of §3.1 holds.
2. Closure: |r| ≤ max(0.010, 2·SE_r). SE_r combines SE_obs, the paired swap SE and the finite-sample CDF error of §3.2 in quadrature.
3. A named component (pA or pp) has a main effect ≥ 70 % of Δ(both), and |I| ≤ 30 % of Δ(both).
4. The D-E5 "identified" signature of §2.2 for the named hypothesis holds. A pA component is attributed as **"H2 (p+A shape)"** only. The core/surface characterisation and the peak-versus-gain ordering (2-paired-SE margin, else "indeterminate") are reported descriptively. If the H1 signature fires, D-E6 adds `swap-sigma` (per-element σ_el scaled by the D-E5 rate ratio, interpolated in E) and `swap-all` (σ + pA + pp). Rules 2–3 are then evaluated on `swap-all`, with the σ component as an additional main effect.
5. Every other hypothesis meets its "excluded" column in §2.2.

Outcomes: "partially attributed" if 2 holds and 3 or 5 fails; "unexplained" if 2 fails. D-E3 never enters this rule. The swaps are diagnostics, never a model: a Geant4-derived shape cannot become the transport model (M4 is rejected). The result is recorded in decision 0041 (new dated section) and in the results README as exploratory. V5, V11 and V10-A stay failed (Amendment 14 (a)2). An "established" attribution authorises a model-change decision (§4), never a passing record.

## 6. Pre-registration (Amendment 19, draft text; committed before D-E0 runs)

> **Amendment 19 — V3-005D elastic attribution (exploratory, report-only).** Rows D-E0 to D-E6 are exploratory and never evidence. They change no frozen row, tolerance, base or verdict; V5, V11 and V10-A at dd976c5a stay failed.
>
> (a) **Plan of record:** `validation/plans/v3-005d-elastic-attribution.md`, sha256 `<filled at commit>`, committed in the same commit as this amendment, before D-E0. Supporting already-seen scripts are in `validation/plans/v3-005d-elastic-attribution/`. **Kernels and signatures:** as §2.1–§2.2 of that file. Definitions K_w (w = 1, 2, 3, 5 mm), K*, G, θ_w (w = 2, 3 mm), T (3·θ_med). Thresholds as in the "identified" / "excluded" columns.
>
> (b) **D-E3 probes:** warp-cpu float64, single process, 1e5 histories, 20 batches, V5 geometry r = 20 cm, 150 and 200 MeV. Variants:
> - `all`;
> - `pp-only` (σ_el = 0 for every target ≠ H-1);
> - `pA-only` (σ_pp = 0);
> - `pA-sigma0.8` (σ_el × 0.8 for every target ≠ H-1);
> - `pA-Rx1.2` (black disk at 1.2·R(E));
> - `pA-s0.5` and `pA-s0.9` (sharp form × exp(−q²s²), s in fm);
> - `pA-gauss` (exp(−q²R²/4));
> - EM-only.
>
> 256 equiprobable μ_CM bins per node as in the pinned table; normalisation exact; σ unchanged except where stated.
>
> (c) **D-E4:** rows 14–23 as in §3.2 of the plan of record. Water rows 14–18 and graphite rows 19–23 at 30 / 50 / 100 / 150 / 200 MeV. Primaries per case: **2e7 at 30 and 50 MeV with a 2 mm target, and 4e6 at 100, 150 and 200 MeV with a 1 cm target.** Module `g4h-elastic_HP` only (declared fallback in §3.2); KE filter E₀ − 1e-4 MeV; kinematic-curve classification; seed 20270000 + 1000·row + 1 (k = 2 only if a case is split in two equal halves).
>
> (d) **D-E5:** 1e6 table-sampled events per case.
>
> (e) **D-E6:** swaps `swap-pA` (O CDF from water rows 14–18; graphite rows 19–23 used for D-E5 signatures only), `swap-pp` (H-1 from rows 14–18 above 16.26° CM), `swap-both`; q-scaled interpolation; ionmc shape below 30 MeV; same sizes and seeds as (b).
>
> (f) **Seeds:** investigation family 2099xxxx, starting at **20990100**. D-E3 and D-E6 use one seed per energy: **20990100** (150 MeV) and **20990101** (200 MeV), shared by every variant and by EM-only (CRN). D-E5 uses **20991000 + case index** (0–9, ordered as rows 14–23). The family is disjoint from:
> - every consumed base;
> - the declared fresh bases 20511004 / 20521004 / 20531004 / 20541004 / 20551004 and their r_index blocks (envelope base + 0 … 20999, maximum 20572003); the re-qualification bases of record are lv 20531004 / hr 20541004 (route document §4.4; 20511004 / 20521004 are retired unused, declared before any observation);
> - the rehearsal seeds used (20505000–20525999);
> - E-shape 20509500 + k;
> - the TOPAS rows.
>
> Disjointness is asserted by a test under `tests/` only (outside the rehearsal digest). The `run_suite.py` seed-guard extension (hashed, l.1303–1313) is deferred to the re-qualification commit, which changes hashed code anyway. The D-E runs are executed by `validation/scripts/elastic_attribution/` and never through `run_suite.py` or an evidence step. Disclosure: numpy seed 20990001 was used once in a planner scratch calculation (analytic V7-R binomial study, no transport, `scripts/candidates2.py`); it lies below the family start and is not used.
>
> (g) **The §5 rule** (rules 1–5, thresholds 0.010, 2·SE, 70 %, 30 %, the signature thresholds, the conditional `swap-sigma` / `swap-all`) and **the §4 pre-commitment** (M2 with the §4.1 inputs, conventions and rules; M3 only on the §4.2 triggers) are frozen with this amendment. No §4 or §5 parameter changes after D-E0 is run.
>
> (h) **Already seen and declared:** the V10-A per-point ratios and first minima, the built-table moments, the shape-moment and removal-kernel estimates of rev 1 / the review (crude; recomputed in D-E0), the reviewer's exceedance and median ratios per shape, and P(pp).
>
> (i) **Timing:** D-E0 to D-E3 and D-E5 to D-E6 are sandbox-only, run one at a time and not concurrently with any host job. D-E4 consists of 10 reference-engine calls, each ≤ 1500 s, one at a time, preceded by a 1e5-primary smoke case.
>
> A later model change goes through a new decision and a pre-registered re-qualification at the declared fresh bases of the integration route (lv 20531004, hr 20541004; V6-F, if the operator selects it, at 20551004) with **unchanged** tolerances.

## 7. Operator proposals (not autonomous)

1. **V10-A after a shape change.** Every candidate is chosen after V10-A was observed. Under the (e)9 precedent its V10-A results are exploratory, and shape gating "rests on V6 and V11". V6 is permanently failed (not sized), and V6-F is out of this task (Amendment 18 (b)5). Merge-gate entry (a)6 "V10-A" therefore cannot be met as a blind gate. Options: (i) evaluate V10-A for the new model with the frozen tolerances and a data-seen disclosure, and let it count toward (a)6; (ii) keep it exploratory and re-scope the p+A elastic claim.
2. **V10-A (c)** fails on the shared Barashenkov σ_el(p+C, 100 MeV) = 253 mb against 216 ± 21 mb (+17.1 % > 12 %), independently of the shape. This plan changes no σ_el.
3. Any change of the V11 or V5 tolerances, the ladder or stop rule, or the evidence class of V5/V11.

Autonomous: Amendment 19, D-E0 to D-E6, the shape decision under §4 (M2 pre-committed), the implementation, review, and the pre-registration of re-qualification rows with unchanged tolerances.
