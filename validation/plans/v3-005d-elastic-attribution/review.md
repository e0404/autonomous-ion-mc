# Independent review: V3-005D elastic-attribution plan (elastic-attribution-plan.md)

Reviewer: independent Claude reviewer. Repo: v3-005b worktree, HEAD 749f768c. Read-only. Scratch checks in scratchpad/v3-005d/rev/ (mom.py, removal.py, walk.py).

## A. Numerical spot-checks

| Claim | Status | Evidence |
|---|---|---|
| V11 F ionmc 0.8568/0.8441, TOPAS 0.8916/0.8874, gains 0.0221/0.0267 vs 0.0163/0.0184, rows 12/13 0.8891/0.8906 | reproduced | 16-v11-compare.txt JSON |
| "F_ionmc − F_TOPAS = −0.039/−0.049" (§0) | mislabelled | these are the relative diffs F_i/F_T − 1 (kind "relative"); the absolute diffs are −0.0348/−0.0433, consistent with Δ(1−F) = +0.035/+0.043 in §2 |
| V5 TOPAS peak/plateau −3.66/−4.89 %, ΔIDD plateau +10.35/+6.72 %, R80 −0.30/−0.32 mm, total −0.30/−0.49 %, MCsquare −4.78/−5.34 % | reproduced | 21-v5-compare.txt. The step files sit in RUN-20261012T143203Z-1b3b1ecf (21), RUN-20261012T143235Z-2b818275 (22) and RUN-20261012T142343Z-081123c3 (03), not under RUN-…a8e41216 as §1 says |
| Table moments p+O ⟨1−μ_CM⟩ 0.118/0.093/0.082, θ_rms 29.9/26.7/25.5°; p+C 0.157/0.115/0.098; H-1 σ 27.3/26.3/25.6 mb | reproduced (θ_rms 30.0/26.8/25.5) | npz edges_mu, uniform-in-μ within the 256 bins |
| Black-disk reproduction 0.092/26.6°; s = 0.5 → 0.038/16.1°; s = 0.9 → 0.025/12.8°; R = 3.3 → 0.075/24.0° | reproduced | shape_moments.py re-run |
| BGG pp 27.6/23.9/23.3 mb | reproduced | 1.0115·(23 + 50·√(ln(0.73/p))^7) |
| P(pp) 0.0375 (ionmc) | reproduced: 0.0374 | Bragg–Kleeman α = 0.0022, p = 1.77, table σ_H. Note: a large part comes from 12.5–30 MeV (σ_pp 190–90 mb), where the transfer is irrelevant to the peak. P(pp) = 0.0394 (BGG) and the 0.73/0.76 MeV transfer are not in any script; source them |
| Recoil 0.249/0.300 MeV per primary; plateau share 0.19/0.13 % | reproduced | v5-*-on.json recoil_idd_batches |
| V10-A per-point +44…+81 % (15°), +79…+191 % (20°), partial integrals | not re-verified | 03-v10-a.txt; the planner must cite the exact per-point fields |
| Expected D-E4 counts p+O ≈ 1.1e4, p-p ≈ 3.2e3 per 2e6 primaries, 1 cm water, 150 MeV | reproduced (1.06e4, 3.2e3) | n_O σ L. The quoted precisions 1.5/2.7 % do not equal 1/√N (1.0/1.8 %); state the basis |

## B. Findings

1. **[blocking] Wrong governing moment: peak depression is not proportional to σ_tr.** Elastic p+A scattering removes almost no energy. A proton leaves the Bragg peak only when its end-of-range depth deficit Δ ≈ R_res(E)(1 − cos θ_lab) (plus the small recoil loss) exceeds the peak width (σ_straggling ≈ 1.9/3.1 mm at 150/200 MeV). Peak depression therefore follows a threshold functional, ∫Σ(z) P(Δ > w) dz. That functional is dominated by the core and width of the diffraction peak, i.e. by R (H2b), not by the q⁻³ tail. The plateau gain does follow a σ_tr-like moment (path-length factor 1/cos θ − 1 plus local deposits), so it is tail-sensitive.
   - Deterministic estimate (removal.py, p+O only, 150 MeV, table σ and R(E), straight-line deficit):

     | Variant | Σ⟨1−cos⟩ | Removal at w = 2 mm | Removal at w = 3 mm |
     |---|---|---|---|
     | sharp | 0.0305 | 0.119 | 0.091 |
     | s = 0.5 | 0.0182 (−40 %) | 0.107 (−10 %) | 0.077 (−15 %) |
     | s = 0.9 | 0.0123 (−60 %) | 0.091 (−24 %) | 0.061 (−33 %) |
     | R × 1.2 | 0.0250 (−18 %) | 0.096 (−19 %) | 0.069 (−24 %) |
     | Gaussian exp(−q²R²/4) | 0.0196 (−36 %) | 0.125 (+5 %) | 0.096 (+6 %) |

     The 200 MeV values are similar.
   - Consequences:
     - (i) The H2a prior ("high") and the H2b prior ("matters mainly for the first minimum, not the tail"; "R × 1.2 moves F by ≲ 0.01") rest on σ_tr. With respect to the peak, R is as strong a lever as s.
     - (ii) The predicted H2a signature "`pA-shape-s0.5` brings F within ±0.01 of TOPAS" is not supported. My estimate gives only a −10…−15 % change in removal, against the ≈ −25…−30 % needed.
     - (iii) §5 rule 4 (σ_tr ratio ≥ 1.3) cannot tell H2a from H2b, because both raise σ_tr.
     - (iv) The §0 statement that the 17° → 26.7° θ_rms discrepancy explains the 0.985 vs 0.963 expectation miss is unproven.
   - Fix:
     - Add to D-E0/D-E1 a deterministic peak-removal kernel: the deficit distribution convolved with straggling, computed per shape.
     - Define the signatures separately for peak ratio and plateau ratio. The tail predicts a gain excess larger than the peak excess; the core/R predicts the reverse.
     - In D-E4/D-E5 compare the angular CDF at energy-dependent threshold angles (for example the angle where R_res(1 − cos θ) = 2 mm), not only ⟨1 − cos θ_lab⟩.
     - Re-derive the priors and §5 rule 4 from these quantities.

2. **[blocking] The §5 "attribution established" rule is not well-posed.**
   - Rule 2 asks that Δ(1 − F) (ionmc − TOPAS) be "explained by a sum of identified components", but D-E3 measures only ionmc-side variants, and no TOPAS-side counterpart is ever produced.
   - None of the D-E3 variants is defined as "Geant4-like". "Explaining" Δ with s0.5 or s0.9 therefore means picking whichever variant closes the gap. That is a selection on V11, which is circular and contrary to the directive.
   - Additivity of the components is assumed but never tested.
   - Rule 5 cannot be evaluated: H2b has no "excluded" signature ("—"), and the excluded signatures of H3 ("`pp-only` share consistent with D-E4 within ±0.005") and H1 refer to quantities that D-E4 does not produce (D-E4 gives thin-target rates and angles, not an F share).
   - Fix: make the core test a pre-declared swap test.
     - Replace ionmc's p+O (and C) angular distribution by the empirical Geant4 thin-target distribution from D-E4, at the same σ. Do the same for p-p.
     - Run D-E3 with each swap, in both orders, with common random numbers.
     - A component is "identified" iff its swap moves F by the declared amount. The swap is a diagnostic, never a model.
     - Give every hypothesis, including H2b, a computable excluded signature.
     - Define how the sum is evaluated, including the interaction term.

3. **[major] The D-E4 event classification will not work as written.**
   - In Geant4 hadronic elastic scattering the primary survives (same track). A creator-process field therefore never tags a deflected primary.
   - p+A recoil ions below the recoil threshold (T_rec ≈ 0.06 MeV at θ_CM = 5°, 150 MeV, O) are deposited locally and produce no track.
   - Classification is clean only in rows 14/15 (elastic physics only, no EM), where any deflection means elastic scattering. In rows 16/17 (EM on), opt4 adds WentzelVI plus single Coulomb scattering, which gives a Rutherford tail beyond 2°. The proposed H5 check ("beyond 2° must not change except by MCS convolution") is then confounded.
   - H5 has a "very low" prior and F cancels EM to first order, so rows 16/17 spend half the host budget on a non-discriminating test.
   - Fix:
     - Drop rows 16/17, or replace them with a run at another energy or with graphite.
     - Classify in rows 14/15 by deflection, plus the presence of a second proton for p-p.
     - Require the event ID and the parent/track ID.
     - Filter the phase space so that undeflected primaries are not written. 2e6 records per case, 4 energies and 3 seeds is a large output.
     - Verify against the OpenTOPAS 4.3 docs, before the case is frozen: the phase-space option names (`IncludeCreatorProcess`, `IncludeEventID`/`IncludeTrackID`, the process/angle filters); that a physics module list without any EM module is accepted by TOPAS and by the `run_reference_calculation` fingerprint/module allowlist; and the per-call ≤ 1500 s with output size.

4. **[major] The D-E3 seed plan breaks common random numbers and wastes precision.**
   - "One seed per run, 20990000 + 100·v + e" with EM-only as its own variant means EM-only and EM+el no longer share a seed. That contradicts the stated common random numbers, and the variant-to-variant differences (the decomposition target) become independent. SE(ΔF) ≈ 0.003, against a closure band of ±0.010 and a 70 % rule.
   - Fix: one seed per energy, shared by EM-only and every variant (common random numbers across variants), as in R1-D2 (C). Freeze this explicitly in Amendment 19.
   - Also: the planner's scratch candidates2.py already used numpy seed 20990001 (an RNG for an analytic calculation, not transport). Record this, or start the family at 20990100.

5. **[major] The "no `src/` change" claim is inaccurate as stated.**
   - (a) The cited E-shape precedent was a `src/` change: `exploratory_forced_isotropic` in `src/ionmc/nuclear/elastic_build.py`.
   - (b) The rehearsal digest hashes `src/ionmc`, `validation/scripts/transport` and `benchmarks/transport` by prefix, plus enumerated files (run_suite.py SOURCE_PREFIXES, hashed_source_digest). Only plans, decisions, docs, tests and README files are excluded. A script added under `validation/scripts/transport/` would therefore invalidate the V7-R rehearsal fixture (two ≈ 1260 s host shards).
   - Feasible route without a `src/` change: a standalone script outside the hashed prefixes (for example `validation/scripts/elastic_attribution/`, or under `tests/`) that:
     - post-processes the pinned table (edges_mu for the shape variants; sigma_barn zeroed or scaled for pp-only, pA-only and σ × 0.8);
     - adds an `exploratory` key;
     - recomputes npz_sha256 and the id with `ionmc.nuclear.build.table_id`;
     - loads the result through the existing `PhysicsConfig(elastic_table_id=…, allow_exploratory_table=True)` (config.py:113/140/592, elastic_tables.py:138).
   - Caveats to verify before running:
     - The copied qualification flags (shape_normalised and others) are not earned by the variants. Disclose this.
     - The exploratory key also disables the loader's below-domain assertion. Assert it in the script.
     - The warp backend must read the same arrays.
   - State this route in the plan instead of "test-only builder options".

6. **[major] Amendment 19 is not yet a sufficient pre-registration.** Besides the seeds, it must freeze:
   - the exact variant definitions and parameters (s, R factor, σ factor, swap construction), the history count and the backend;
   - the D-E4 case matrix: rows × energies × seeds. Four energies with k = 1..3 per row gives 12 cases per row, but only k = 1..3 is declared. Seeds per case must be distinct, or the sharing declared;
   - the primaries per case, the physics modules, the classification rules and the threshold angles;
   - the D-E5 sample size and seeds per case;
   - the peak-removal kernel definition (finding 1);
   - the §5 rule with computable signatures for every H (finding 2);
   - the order: D-E0 to D-E2 outputs are seen before D-E3. State that no §5 or §4 parameter may change after D-E0 to D-E2 are seen, or freeze the amendment before D-E0.

   Seed disjointness checked: 2099xxxx is clear of the consumed ionmc bases (up to 20491004 plus blocks), the rehearsal range 20505000–20525999 and the fresh bases 20511004/20521004 + ≤ 20999 (max 20542003). TOPAS rows 14–17 (20284001–20287003) are unused in reference_cases (rows 1–4, 9 and 10–13 are used).

7. **[major] The §4 selection rule is blind to V5/V11 but not to V10-A.** It selects M1 vs M2 on D-E1, which uses the V10-A data the planner has already examined (§0 quotes the model/data per point). That is acceptable only if it is declared as data-seen. Also, M1's s ≈ 0.9 fm and R must come from a named electron-scattering source (e.g. a de Vries 1987 compilation), fixed in the decision before D-E1 is computed for M1. Otherwise the choice of s is free.
   - The candidate list omits that the observed effect may be core/R-driven (finding 1). The independent M1/M2 parameters fix R too, so state which R rule M1 uses (the σ_nonel rule versus the charge radius).

8. **[minor] Comparator like-for-like.** TOPAS g4h-elastic_HP for protons is BGG σ plus G4ChipsElasticModel (no angular cut, Coulomb–nuclear interference not separated; p-p over all angles). ionmc uses Hale NI above θ_CM 16.26° and pure-nuclear black-disk p+A. The thin-target comparison must apply the 16.26° cut to Geant4 p-p before computing rates. State that the V10-A data include Coulomb and interference near the minima.

9. **[minor]** In §0, "the V5 residual is entirely elastic" needs its uncertainty (−3.66 % vs −3.90 %; −4.89 % vs −4.88 %, each with SE ≈ 0.3–0.4 %). The 0.985/0.986 expectation reconciliation should cite its source line.

10. **[minor]** The cost estimate for D-E3 assumes 6000 hist/s on the host. The sandbox rate is unmeasured, and the python fallback of 2e4 histories has SE(F) ≈ 0.0035, which is too coarse for the ±0.010 closure. Make 1e5 warp-cpu a requirement or adjust the rules.

Compliance: nothing in the plan relaxes a tolerance, re-uses an evidence seed or creates a passing record. V5, V11 and V10-A stay failed, and §7 correctly routes the scope questions to the operator. The plan satisfies the directive's intent ("elastic-model attribution before another qualification") only once findings 1–2 are fixed.

## C. Verdict

**REVISE.** There are two blocking findings (1: the governing moment; 2: the ill-posed and potentially circular attribution rule) and five major ones (3–7). No compliance violations.

---

# Round 2: re-review of Revision 2 (2026-10-14)

Scope: rev 2 checked against findings F1–F10, with code claims verified at worktree HEAD 749f768c. Additional scratch checks:
- `$TMPDIR/exc.py`: exceedance and median per shape;
- npz first-minimum location;
- 03-v10-a.txt `theta_min_deg`.

## Disposition of round-1 findings

| # | Status | Comment |
|---|---|---|
| F1 governing moment | **resolved in structure**; residual R2-2 and R2-6 | K_w, K*, G, θ_w and T are defined; H2a and H2b are split; the θ_rms claim is withdrawn |
| F2 circular rule | **resolved** | the 2×2 swap test with the Geant4 empirical CDF at unchanged σ, an interaction term, closure \|r\| ≤ max(0.010, 2·SE_r), and probes excluded from §5. SE_obs 0.0026 / 0.0034 verified (0.002885·0.8916, 0.003838·0.8874). Residual: R2-5 |
| F3 D-E4 design | **resolved**, minor residual R2-7 | EM-free, vacuum sphere, kinematic-curve classification, rows 16/17 dropped, five verification items listed, fallback declared |
| F4 seeds | **resolved** | one seed per energy (20990100 / 20990101) shared by every variant and by EM-only (CRN); D-E5 uses 20991000 + case; TOPAS 20284001–20293001 are unused in `validation/reference_cases`; the 20990001 disclosure is present |
| F5 table route | **resolved**, with one contradiction (R2-3) | Verified: SOURCE_PREFIXES = `src/ionmc`, `tests/ionmc`, `validation/scripts/transport`, `benchmarks/transport` (run_suite.py l.95–100), so `validation/scripts/elastic_attribution/` lies outside both the suite identity and the rehearsal digest. `PhysicsConfig.allow_exploratory_table` exists (config.py l.113) and requires nuclear, elastic and an explicit `elastic_table_id` (l.140–146); `ElasticTable.load(None, tid, allow_exploratory=…)` (config.py l.592; elastic_tables.py l.138–142) refuses an exploratory table without the flag. Note: load is called with cache_dir None, so a "scratch cache copy" is honoured only through `IONMC_CACHE_DIR` (cache.py l.37–48), and that copy must also hold the nuclear table and every other derived input the run loads |
| F6 Amendment 19 | **mostly resolved**; R2-1, R2-3 | sizes, rules, seeds, timing and seen items are now frozen before D-E0 |
| F7 §4 | **partially resolved**; R2-4 | the data-seen declaration and the de Vries 1987 source are present; the R rule is stated, but the M1 construction is inconsistent |
| F8–F10 | resolved | |

**Citing the reviewer's removal-kernel numbers without re-running them** is acceptable for a plan, because they are declared as already seen, are not used as thresholds, and D-E0 recomputes them in a committed script. Two conditions apply (minor, R2-8):
- (a) Label them for what they are: crude estimates (Bragg–Kleeman range, straight line, p+O only, step-function threshold, no straggling, no p-p), not a reproduction target. D-E0 will differ, because it uses the project range table and K*. Replace "reproduced in D-E0" with "recomputed in D-E0; differences expected".
- (b) If they stay cited, archive `rev/removal.py` with the plan, or drop the numbers.

## Remaining findings

**R2-1 [major] Amendment 19 freezes by reference to a scratch file.**
- Clause (a) reads "as `elastic-attribution-plan.md` rev 2 §2.1–§2.2", and (c) reads "as in §3.2". The file lives in the session scratchpad and not in the repository, so the pre-registration is not self-contained.
- Fix: commit the plan as `validation/plans/v3-005d-elastic-attribution.md` (as V6-F did with `v3-005c-v6f-followup.md`), and cite its path and sha256 in the amendment, in the same commit, before D-E0.

**R2-2 [major] The H2a/H2b discrimination signatures are confounded.**
- Exceedance and median ratios do not separate a surface change from a radius change. Black disk at 150 MeV for p+O, CM frame, θ_2mm = 9.8°:

  | Shape | P(θ > θ_2mm) | Ratio to sharp | Median θ | Ratio to sharp |
  |---|---|---|---|---|
  | sharp | 0.689 | — | 13.3° | — |
  | s = 0.5 | 0.634 | 1.09 | 12.0° | 1.11 |
  | s = 0.9 | 0.565 | 1.22 | 10.7° | 1.24 |
  | R × 1.2 | 0.589 | 1.17 | 11.2° | 1.19 |

- A pure surface (H2a) difference of s ≈ 0.9 fm therefore meets the H2b "identified" threshold (exceedance ≥ 1.20). With T ≥ 1.5 it would be labelled "H2 (both)".
- §4 selects the model on the basis of this label (M2 if H2b), so the misattribution feeds directly into the model choice.
- Fix: attribute at the level "H2 (p+A shape)" in §5. Make the H2a/H2b split descriptive, via a pre-declared two-parameter characterisation of the Geant4 thin-target CDF:
  - a first-minimum or forward-peak width estimator for the core;
  - a large-q log-slope for the surface;
  - or a declared (R_eff, s) fit of the Helm/disk form to the Geant4 CDF (a fit to Geant4, not to data; report-only).
  - Do not use the split as a gating or selection label.
  - The peak-versus-gain ordering test also needs a statistical margin (ordering beyond 2 paired SE, else "indeterminate").

**R2-3 [major] The seed-guard change contradicts "no rehearsal-digest change".**
- Amendment 19 (f) states "The seed guard refuses the family for every evidence step". The seed guards live in `validation/scripts/transport/run_suite.py` (l.1303–1313, hashed), so implementing this changes the V7-R rehearsal digest. That costs a fixture regeneration, two ≈ 1260 s host shards.
- Fix, either:
  - put the disjointness assertion in `tests/` only (excluded from the rehearsal digest) and defer the guard edit to the later re-qualification commit, which changes hashed code anyway; or
  - declare the fixture regeneration and its host cost.

**R2-4 [major] The §4 model rule is inconsistent, and it is informed by the reference engine.**
- (a) M1 combines the 2-D black-disk amplitude 2J₁(qR)/(qR) with the 3-D Helm radius R_H, defined for 3j₁(qR)/(qR) (first zero 4.49 versus 3.83). It also uses a charge form factor as a hadronic profile. The rationale "M1 cannot correct a core-width error" also contradicts the rule's own "R_H replaces the σ_nonel radius rule".
- (b) Already-seen archived evidence that the core is the dominant error: V10-A `theta_min_deg` puts the measured first minimum at 17.17° (C0063 @ 135 MeV) and 14.99° (C0141 @ 179.9 MeV). The pinned table's first p+C density minimum is at 42.0° and 35.9° CM (analytic J₁ zero 41.5° / 35.0°, R = 2.25 / 2.27 fm), a factor of about 2.4 in angle.
  - If this is confirmed (D-E1 must report it; verify that the data minimum is not a Coulomb-interference dip), no charge-radius disk (R_H(C) ≈ 2.3–2.5 fm) can move it, so M1 is inadequate a priori.
- (c) Selecting M1 versus M2 on D-E5, which compares with Geant4, the physics of the V11 reference, makes the structure choice informed by the gating reference.
- Fix:
  - Select a priori on data-side criteria in D-E1 (first-minimum position and partial moments against measured data, declared as data-seen), not on D-E5.
  - Or pre-commit to M2 (an optical-limit Glauber profile, which carries the NN range and fixes core and surface together) with M3 as fallback, and drop M1 or redefine it consistently (a Helm 3-D form factor with a declared strong-interaction radius).
  - Record the Geant4-informed lineage in the decision if D-E5 is used at all.

**R2-5 [minor] H1 cannot be "established" by the §5 rule.** The swaps keep ionmc's σ. If σ were the cause, closure (rule 2) fails and the outcome reads "unexplained".
- Fix: add a pre-declared conditional `swap-sigma` (per-element σ scaled by the D-E5 rate ratio), run only if the H1 "identified" signature fires, or state that H1 is attributed outside rule 2.

**R2-6 [minor] The K* kernel form.** exp(−Δ²/(4σ_R²)) is the overlap integral of two Gaussians. The contribution of a shifted proton to the IDD at the original peak depth goes as exp(−Δ²/(2σ_R²)) for a Gaussian peak, and the Bragg peak is asymmetric.
- Better predictor: the archived EM-only ionmc IDD as the kernel, predicted peak ratio = ⟨IDD_EM(z_pk − Δ)⟩/IDD_EM(z_pk). It is deterministic and needs no new run.
- The kernel only validates D-E3 and is not in §5, so this is minor.

**R2-7 [minor] D-E4 plural scattering is understated.** "(nσL)² < 5e-5" holds at ≥ 150 MeV, but nσL for O in 1 cm water is 0.035 at 30 MeV, and for C in graphite about 0.09. The fraction of plural events among scattered events (≈ nσL/2) is therefore 2–5 % at 30–50 MeV, off both kinematic curves.
- Fix: state the per-energy values; use a thinner target (≤ 2 mm) at 30–50 MeV or report the off-curve fraction; apply the same KE cut to the D-E5 sampling. Also state that the finite Geant4 CDF sample (~2e4 events) adds a fixed-realisation error to the swaps that CRN batch SEs do not capture; it is estimated to be small (≲ 1e-3 in F).

**R2-8 [minor]** Caveats for the cited removal-kernel numbers (see above).

## Round-2 verdict

**REVISE (minor revision). No blocking issues remain.** F1 and F2 are resolved in substance, and nothing relaxes a tolerance, reuses a seed or creates a passing record. Before Amendment 19 is committed, the four major items must be fixed:
- R2-1: commit the plan into the repo and reference it by sha;
- R2-2: the H2a/H2b split is not identifiable from the declared signatures; gate on "H2" only;
- R2-3: the seed-guard edit breaks the rehearsal digest;
- R2-4: the M1 construction is inconsistent, the core error is already evident from the archived first minima, and selection must not rest on the Geant4 comparison.

---

# Round 3: re-review of Revision 3

## Fixes verified (checked against the text and the code or data, not taken from the change list)

- **R2-1, real.**
  - Amendment 19 (a) cites `validation/plans/v3-005d-elastic-attribution.md` with its sha256 at commit, in the same commit, before D-E0.
  - The supporting scripts go to `validation/plans/v3-005d-elastic-attribution/`. That prefix is excluded from the rehearsal digest (`REHEARSAL_DIGEST_EXCLUDED`) and is not a SOURCE_PREFIX.
- **R2-2, real.**
  - §2.2 now has a single H2 row, and §5 rule 4 attributes only "H2 (p+A shape)".
  - The core/surface split is report-only. It covers θ_med, a first-minimum estimator with 1° bins, the large-q slope and a Helm (R_eff, s) fit to Geant4.
  - The peak-versus-gain ordering uses a 2-paired-SE margin and is otherwise "indeterminate".
  - §4 no longer reads the label.
- **R2-3, real.**
  - The disjointness test sits under `tests/`. `tests/` is in SOURCE_PREFIXES (the suite identity changes, which is harmless) but is excluded from the rehearsal digest.
  - The edit to `run_suite.py` l.1303–1313 is deferred.
  - The D-E runs never go through `run_suite.py`.
- **R2-4: the planner's correction is right, and my round-2 statement was wrong.**
  - `03-v10-a.txt` gives target **O-16** for C0063002@135 (θ_min 17.17°) and for C0141002@179.9 (θ_min 14.99°).
  - The pinned table's first p+O density minima are 35.3° / 29.2° CM; the analytic J₁ zeros are 34.2° / 28.8° for R = 2.66 / 2.68 fm. This matches the planner's ≈ 34° / 29°.
  - A disk would need R ≈ 4.9 fm to put the zero at the data minimum.
  - My "p+C 42.0° / 35.9°" compared the C table with O datasets. I withdraw it; the conclusion (the core is about 2× too wide, to be verified in D-E1) stands with the O values.
  - The D-E1 diffraction-minimum criteria are concrete: ≥ 2 adjacent points, deeper than 2σ, Rutherford share < 20 %.
- **R2-5, real.** A conditional `swap-sigma` / `swap-all` is added. Otherwise H1 is reported as "not attributable".
- **R2-6, real.** K* now uses the archived EM-only IDD (`v11-*-emonly.json`). The Gaussian form is corrected to 2σ² and is report-only.
- **R2-7, real.**
  - The plural-scattering numbers reproduce:
    - O in 1 cm water: 0.0118 / 0.0053 / 0.0034;
    - C in 1 cm graphite: 0.025 / 0.0125 / 0.0084;
    - 2 mm targets: O 0.0071 / 0.0055, C ≈ 0.018 / 0.012.
  - The off-curve exclusion and the same KE cut also apply to D-E5.
  - The CDF bootstrap error enters SE_r.
- **R2-8 and the cache note, real.**

## Remaining findings

**R3-1 [major, text-only] M2 is not yet specified well enough to be built blind.** The structure is correct: the optical-limit amplitude F(q) = ik∫b db J₀(qb)[1 − e^{iχ}], with iχ = −½σ_NN(1 − iα)T⊗f_NN, the Barashenkov level and a normalised shape. Before Amendment 19 is committed, three inputs are still open:
- (a) **The NN source is unnamed.** The plan says "a named NN compilation … pinned before any M2 shape". Name it now: the source of σ_pp and σ_np, α_NN and β_NN versus T_lab. It must cover the whole model domain, 6.6–250 MeV. Most classic Glauber tabulations start at about 100 MeV, so a declared extrapolation or substitute rule is needed for ≈ 6–100 MeV, where σ_np rises steeply.
- (b) **The density per target is not fixed.** Specify, per table target (C, N, O, Al, Si, P, Ca):
  - which de Vries parametrisation (HO, 2pF, 3pF or SOG) is used and its table entry;
  - the proton-size unfolding rule and the r_p value.
  - "HO / Fermi" is not unique, and N ≠ Z does not arise only for the surrogates (Cl, K and Ar map to the P-31 and Ca-40 targets).
- (c) The CM energy and frame at which the NN parameters are evaluated (T_lab per nucleon) and the k used in F(q).

Without these items the pre-commitment leaves free choices that could be made after seeing D-E1. This is fixable in text; no new review round is needed if it is added exactly.

**R3-2 [minor]** Amendment 19 (c) says "4e6 primaries each" and then "primaries 4e6 and 2e7 respectively". Delete the first.

**R3-3 [minor]** The M3 trigger "optical-theorem consistency within 5 %" is close to an identity for a self-consistent calculation. Keep it as a numerical check, but note that it cannot detect a wrong input.

## Round-3 verdict

**ACCEPT, conditional on the R3-1 text (naming the NN source with its < 100 MeV rule, the per-target density parametrisation with its unfolding, and the NN energy/frame convention) and R3-2 being inserted before Amendment 19 is committed.** No blocking issues remain. The attribution design (D-E0 to D-E6, §5) is sound and pre-registrable. Nothing relaxes a tolerance, reuses a seed or creates a passing record. The recorded V5, V11 and V10-A failures are preserved.

---

# D-E0–D-E2 review (worker implementation at 9a7b142 plus untracked files; read-only)

Scope:
- code: `validation/scripts/elastic_attribution/{common,de0_kernels,de1_shapes,de2_pp}.py`;
- outputs: `validation/results/elastic_attribution/de{0,1,2}/`;
- test: `tests/ionmc/test_elastic_attribution_de0.py`.

## 1. Kernel implementation against plan §2.1; the K* sign

- K_w, θ_w and T are implemented as defined (de0_kernels.py l.75–157): exact two-body kinematics through `ionmc.nuclear.elastic_kin`, the Bethe CSDA range, a 0.98·R₀ upper limit and 300 depth nodes. For p-p the faster proton continues.
- **The K* sign: the "shifted" form (z_pk + Δ) is the physically correct one, so the plan's formula has a sign error.**
  - A proton whose range is shortened by Δ has the depth-dose d(z) ≈ IDD_EM(z + Δ). Its contribution at the original peak depth is IDD_EM(z_pk + Δ), which lies on the distal fall-off.
  - IDD_EM(z_pk − Δ) samples the proximal rise. That is the contribution of a proton whose range was *extended* by Δ, so it systematically under-predicts removal.
  - The plan's own prose ("moves that proton's contribution from … z_pk to z_pk − Δ") describes the shifted physics correctly. Only the formula `IDD_EM(z_pk − Δ)` has the wrong sign.
  - The worker's reading of the formula is accurate, and keeping it as the definition of record while reporting the correction is the right procedure (Amendment 19 (g) freezes it).
  - **Required:**
    - Record a dated note (decision 0041 or the D-E0 record, not an edit of the frozen plan) stating that the frozen K* formula has a sign error and that K*_shifted is the physically meaningful predictor.
    - Base no conclusion on K*_literal. In particular, drop the D-E0 headline "literal K* over-predicts F 0.904/0.892", which is an artefact of the sign.
    - In D-E3, report both forms, and validate against the shifted one.
    - Add a unit test with a synthetic IDD asserting d(z_pk) = IDD(z_pk + Δ). Currently no test pins the sign.

## 2. Minima analysis and its consequences

- **The worker is right.** The V10-A rule's θ_min,1 = 17.17° (C0063 @ 135 MeV, O-16) is an artefact of duplicate angles: 17.17° gives 138.9 mb/sr and 17.24° gives 141.9 mb/sr, a 0.74σ "drop". The data fall monotonically from 570 mb/sr at 6.6° to 1.8 mb/sr at 38°. The same holds for 14.99° (C0141 @ 179.9: 143.1 vs 144.8 at the same angle; drop 0.48σ).
- Under the plan's own frozen D-E1 criteria (≥ 2 adjacent points, > 2σ, Rutherford share < 20 %), neither is a diffraction minimum.
- The worker's significance rule (≥ 10 % and ≥ 3σ below both neighbours, duplicates merged) differs from the frozen rule. **Minor:** report the frozen-rule result too; the verdict is the same.
- **Consequences:**
  - (a) **My round-2/round-3 statement "core about 2× too wide (R ≈ 4.9 fm)" is withdrawn.** It rested on a spurious minimum. The plan's §0 and §2.2 prior text ("first minima ≈ 2× too wide") needs a dated note recording that D-E1 refuted it.
  - (b) Reviewer observation, report-only: C0063 flattens at 30–34° (4.56 → 3.83 → 3.26 mb/sr, then 1.82 at 38.5°). This is consistent with a filled-in minimum near the table's 34°. The core position may be roughly right, and the discrepancy may lie in the forward fall-off and the depth of the minimum. That fits the model being +45 % / +148 % above the data at 15° / 20°.
  - (c) The §4.2 M3 trigger is unaffected, because it uses q_H from the de Vries density and not the data. The §5 signatures do not use the measured minima and are unaffected. The descriptive first-minimum estimator stays.
  - (d) **Important, record-keeping.** The recorded V10-A evaluation at dd976c5a used the spurious θ_min,1 to set the (a) cut at 0.8·θ_min = 13.74° (`03-v10-a.txt`: "theta >= 0.8 theta_min,1 (measured: 13.74 deg)"), which excluded the 15°/20° (a) points of C0063 and C0141. This is an implementation/specification defect of the V10-A minimum rule (duplicate angles not merged). Record it as a finding against the V10-A row. The recorded V10-A failure stands unchanged; nothing is re-evaluated or relabelled.

## 3. Reconciliation and kernel adequacy

- With the correct (shifted) sign, the threshold kernel reproduces ionmc's own V11 F well:
  - predicted F = (1 − K*_sh)/(1 + G) is 0.853 / 0.832, against observed 0.857 / 0.844;
  - peak removal: K*_sh 0.136 / 0.151 against observed 0.124 / 0.133, an over-prediction of about 9 % / 13 %.
  - The kernel is therefore adequate to roughly 10 % for ranking and sensitivity. It is not adequate for the 0.010 closure. That is why §5 rests on the transport swap test (D-E6) and not on the kernel, as designed.
- **The "reconciliation not reproduced" conclusion is incomplete (major for the D-E0 deliverable).**
  - The worker varied only the planning radius (+0.002), but the planning expectation used θ_rms = 17° / 14°, a different shape. The plan asks for "table moments versus planning moments".
  - The s = 0.5 family reproduces those moments (θ_rms 16.0° / 13.7°). From the worker's own table, K*_sh(s0.5) = 0.1225 / 0.1337 and G = 0.0110 / 0.0190 give F ≈ 0.868 / 0.850. Corrected for the kernel's bias on the table shape (+0.004 / +0.012), that is ≈ 0.872 / 0.862, against the planning-implied 0.878 / 0.875.
  - So the planning shape moments account for most of the miss: roughly 75 % at 150 MeV and 55 % at 200 MeV (reviewer estimate from D-E0 numbers).
  - **Fix:** add the planning-moment-equivalent shape to the reconciliation and state the remainder.
- **G under-predicts the gain** (0.0134 / 0.0205 against 0.0221 / 0.0267). The operational G omits at least the higher stopping power of the continuing proton after a p-p event (E_fast = E − T_slow). Given the plan's under-specified wording, report G as incomplete (minor). It affects only the descriptive peak-versus-gain ordering.
- **Implication for H2:**
  - The shift needed to reach TOPAS's removal (ratio 0.755 / 0.722) exceeds what even s = 0.9 (0.79 / 0.77) or R × 1.2 (0.85 / 0.85) gives.
  - p-p carries only about 22 % of K. A +10 % p-p rate (D-E2) moves K by about 0.003, which is about 10 % of Δ_pk: not negligible for the 70 % rule.
  - Independent data-side evidence for H2 remains strong: model/data +37…+193 % at 15–20° for all ≥ 135 MeV O datasets, and partial σ_tr above the data.
  - The shape hypothesis is therefore **not weakened** by D-E0. The specific "core 2× too wide" argument is gone, and the 2×2 swap test is still the only decisive step: it separates shape (pA) from p-p and measures their interaction against Geant4 at fixed σ.

## 4. Numerical spot checks

All recomputed from the archived files.

| Quantity | Result |
|---|---|
| C0063 model/data at 15.11° | 339.4/233.6 = +45 % ✓ |
| C0063 model/data at 19.39° | 186.4/75.19 = +148 % ✓ |
| C0141 model/data at 15.0° | 196.3/143.1 = +37 % ✓ |
| C0141 model/data at 20.72° | 54.85/24.57 = +123 % ✓ |
| Table first p+O minima | 35.2° / 29.2° (mine 35.3° / 29.2°) ✓ |
| σ_pp table/BGG at 150 MeV | 26.26/23.86 = 1.101 ✓ |
| P(pp) | 0.036 table vs 0.038 BGG (Bethe range; my Bragg–Kleeman estimate was 0.0374) ✓ |
| Crude removal table | reproduced exactly, and the plan-definition K_w (O-only) is within 2 % of it ✓ |
| SE of the duplicate-angle drops | 0.74σ and 0.48σ ✓ |

**Not reproduced: the headline "partial σ_tr 1.9–6.8× data".** The D-E1 table gives model/data:
- 1.48 (C0063), 2.0 (C0141), 2.77 (C0148), 4.0 (C0550), 2.2 (O0553 @ 150), 1.55 (C0055 @ 122), 3.5 (C0055 @ 160), 5.6 (C0057), 2.4 (C1420);
- but **0.33 (O0274 @ 96) and 0.004 (O0553008 @ 75 MeV)**.

The O0553008 data are ≈ 8× above both the model and plausibility at 10–16° (5800 mb/sr at 10°, against a model of 694). This looks like a unit or quantity mismatch inherited from the V10-A ingestion. Flag it and exclude it from summaries, and correct the headline to "1.5–5.6× for the ≥ 122 MeV datasets; the two low-energy C sets go the other way, one with a suspected unit error" (minor; record-keeping).

## 5. Compliance

- Every report and JSON carries the exploratory banner, and the input files are listed with sha256 (D-E2 includes the cached PDG object; no download was made).
- The R fit is labelled data-seen.
- H6: the case files differ only in the header, seed and modules (verified in the D-E0 listing).
- No row, tolerance or seed changed. The script and output paths are outside the rehearsal digest. The test is under `tests/ionmc` and is marked `requires_tables`.
- **Major (scope/compliance): the "TOPAS-equivalent surface s = 1.02 / 1.09 fm and radius factor 1.35 / 1.47" in D-E0.** These are shape parameters solved to match the V11 TOPAS peak removal. They are not in the pre-registered D-E0 definitions, and they are exactly the kind of V11-derived model numbers the operator directive keeps out of model choices. M2 is pre-committed and parameter-free, so no harm has been done yet. Remove them from the record, or move them into a clearly quarantined appendix labelled "not pre-registered; V11/Geant4-derived; must not appear in any decision or model input". Never cite them in decision 0041.

## Verdict (D-E0–D-E2)

**Accept as an exploratory record after revision of the reporting.** No blocking issue. The computations are correct, and the minima finding is right and important.

Required before the record is committed:
- (i) Record the K* sign error, lead with the shifted form, and drop conclusions drawn from the literal form.
- (ii) Complete the reconciliation with the planning-moment shape (s ≈ 0.5).
- (iii) Quarantine or remove the TOPAS-equivalent s / R × factors.
- (iv) Record the V10-A spurious-minimum defect and the refuted "2× core" prior as dated notes, with the V10-A verdict unchanged.
- (v) Correct the σ_tr headline and flag the O0553008 @ 75 dataset.

Minor items:
- report the frozen-rule minimum test;
- G incompleteness (the post-p-p stopping-power term);
- a sign unit test for K*.

---

# Round D-E3 review (D-E3 probes and the H5 EM-only check; worktree at d5c8a80 plus untracked files; read-only)

1. **Spec conformance: OK.**
   - All 8 variants and EM-only are present at 150 and 200 MeV, plus `pinned` at 150 MeV.
   - Every run: n = 1e5, 20 batches, seed 20990100 (150 MeV) or 20990101 (200 MeV), shared by every variant and by EM-only (CRN), all `valid`.
   - The runs were made one after another (telemetry 13:04–13:25Z, one process each).
   - The configuration comes from `steps_v5c.v11_config` (warp-cpu float64; emel = `nuc_config_c_elastic_only`; emonly = `nuc_config_c(nuclear=False)`) on `v5b.v5_geometry`, so V5 geometry r = 20 cm.
   - The shape variants rebuild 256 equiprobable μ_CM edges per valid node from `de1_shapes.family_weight` on a 16385-point grid, the builder's resolution. The Jacobian dμ = −4s ds is correct, and the normalisation is checked on a 4× finer grid. σ is touched only in pp-only, pA-only and pA-σ0.8.
   - Stated note: the shape variants apply to every p+A target. Only O matters in water.
   - The summary note says D-E3 never enters §5. **Conformant.**
2. **Harness reuse: OK.**
   - F and gain come from `compare_idd_v5.paired_ratio_jackknife`/`v11_energy`; the IDD, R80 and plateau from the `step_v11_ionmc` construction; recoil from `v5d.REC` (tally-neutral, Amendment 17 (b)5); lateral rms from `v5d.lateral_rms_batches` on the R1-D2 grid.
   - `all` against the archived V11 values (also 1e5 histories, same size; only the seed differs):
     - F 0.8580 ± 0.0014 vs 0.8568 ± 0.0015 (+0.6σ) and 0.8455 ± 0.0027 vs 0.8441 ± 0.0025 (+0.4σ);
     - gain 0.0221 vs 0.0221 and 0.0249 vs 0.0267 (−1.4σ).
   - Reproduced within statistics. `pinned` equals `all` bit for bit (check-base-150).
3. **Variant tables: sound.**
   - The shapes are the same `de1_shapes` families as the D-E0 kernels. D-E0 used 2048 quantile points, against 256 edges in transport; the difference is immaterial at this precision.
   - Each table has an `exploratory` key (variant, source id, purpose, flags disclosed as unearned) and a recomputed id; it loads only with `allow_exploratory_table=True`.
   - The preflight asserts that the scratch `derived/` holds exactly the expected set and logs every loaded table id.
   - The pinned npz is untouched: sha256 1028bb67… unchanged, mtime 2026-10-13.
   - The code and outputs sit outside SOURCE_PREFIXES. The modified test is under `tests/`, which is excluded from the rehearsal digest. **Rehearsal digest not touched.**
4. **Scientific reading** (contrary evidence preserved).
   - (a) **K\* validity.** pred/meas removal is 1.08–1.14 for p+A and 1.01–1.03 for pp-only. This is the D-E0 bias (+9 % / +13 %), so the kernel is validated as a ranking tool. The near-exact pp-only value fits p-p removal being almost binary.
   - (b) **Composition is multiplicative, not additive.**
     - F_pA·F_pp = 0.8986·0.9547 = 0.8579 against F_all 0.8580 (150 MeV), and 0.9044·0.9349 = 0.8455 against 0.8455 (200 MeV). In (1−F) the additive "interaction" is therefore −0.005 / −0.006, which is pure product term.
     - The gains are additive: 0.0081 + 0.0137 = 0.0218 against 0.0221; 0.0048 + 0.0198 = 0.0246 against 0.0249.
     - **Recommendation for D-E6 (§5 rule 3):** evaluate the main effects and I on −ln F (where they add), or report I relative to the product expectation. Otherwise the product term (≈ 15 % of Δ(both)) eats into the 30 % |I| budget.
   - (c) **Gap closed** toward TOPAS F (0.8916 / 0.8874; Δ = 0.0336 / 0.0419), as (F_var − F_all)/Δ:

     | Variant | 150 MeV | 200 MeV |
     |---|---|---|
     | pA-s0.9 | 87 % | 79 % |
     | pA-s0.5 | 49 % | 44 % |
     | pA-R×1.2 | 53 % | 47 % |
     | pA-σ0.8 | 55 % | 38 % |
     | pA-gauss | 4 % | 0 % |

     These are probes, not attribution.
   - (d) **Contrary to H2a as originally stated (the "q⁻³ tail" story).** pA-gauss removes the sharp-edge tail, cutting σ_tr by 36 % and the recoil by 38–42 %, yet F is unchanged (0.8592 vs 0.8580; 0.8457 vs 0.8455). Only the gain falls (0.0163 / 0.0217).
     - This confirms finding F1: F is driven by the core/width, the gain by the tail.
     - The H2 (shape) attribution stands, but any narrative that the peak deficit comes from the tail is refuted by D-E3. Record it.
   - (e) **Contrary evidence for the gain: p-p.**
     - The pp-only gain alone (0.0137 / 0.0198) is 84 % / **108 %** of TOPAS's total elastic gain (0.0163 / 0.0184).
     - p-p carries 62 % / 80 % of ionmc's gain and 26 % / 35 % of its peak removal.
     - Even pA-s0.9 leaves a gain excess at 200 MeV (0.0212 vs 0.0184).
     - **Consequence:** the gain excess cannot be closed by any p+A shape, so H3 (p-p rate +10 % vs BGG, or T_slow/deposit handling) is a live candidate for the gain and is not "low prior". The swap-pp arm of D-E6 is essential. If H3 fails its "excluded" signature, §5 rule 5 will block "established", which is the correct outcome.
   - (f) **Lateral rms.** EM+el against EM is 5.5 vs 1.54 mm (150 MeV) and 7.3 vs 2.47 mm (200 MeV); pp-only alone gives 3.6 / 6.2 mm. This rms is dominated by large-angle tails (p-p partners, wide p+A). It is report-only (R1-D2 (C) definition) and contradicts nothing. It is a reminder that the rms metric is tail-sensitive, not a core-width measure.
   - (g) **R80.** All shifts are ≤ 0.34 mm; the elastic effect on R80 is small, as expected.
5. **H5: the rule is applied correctly.** EM-only peak/plateau ionmc/TOPAS is 1.0086 [1.0048, 1.0124] and 1.0111 [1.0041, 1.0181]; both point and CI lie within ±2 %, so H5 is excluded. Definitions: V11 compare / `compare_idd_v5` (plateau 20–60 mm, ionmc jackknife, TOPAS 3-seed mean, Welch df). Direction ("amplifies Δ_pk": a sharper ionmc EM peak makes the same deficit distribution remove a larger peak fraction) is reasoned correctly and is irrelevant to the exclusion.
   - **Minor:** the column "implied shift of F ratio (1/ratio − 1) = −0.0085 / −0.0110" is mislabelled. F is a ratio, so the EM-only difference cancels to first order; 1/ratio − 1 is the shift of peak/plateau, not of F. Relabel or drop it. If wanted, estimate the second-order effect deterministically with K\*_shifted using the TOPAS EM-only IDD as kernel, instead of the ionmc one.
6. **Git record.**
   - `check-base-150.json` (305 B), `de3-summary.{json,md}` and the de0h5 `git/` pair hold scalar aggregates only: no per-batch or per-depth arrays, no experimental or compilation values.
   - The coordinator's description "ionmc simulation results only" is **not literally true for de0h5**: `h5.json` and the report carry TOPAS EM-only peak/plateau means (6.1374, 6.3900) and ratios. These are aggregate metrics of project-generated reference runs. They are not external data and are allowed under the "aggregate statistics only" rule (decision 0041; README), but describe them correctly.
   - **Important (provenance):**
     - The committed summary carries no seeds, n, batches, variant table ids, npz sha256, script sha256 or code state. The runs used **untracked** `de3_probes.py` and `make_variant_tables.py` at d5c8a80 (meta.json `git_status`). The per-run files that do hold seed, n and ids are git-ignored.
     - `check-base-150.json` asserts bit-identity but records neither the compared array digests nor the seed.
     - **Fix before committing:** add to `de3-summary.json` a provenance block per row (seed, n, n_batches, elastic_table_id, variant npz sha256, run-file sha256, telemetry run label) and globally (commit d5c8a80, the sha256 of the two scripts as run, which must equal the committed files, and the scratch `variants.json` sha). Add the idd/recoil batch-array sha256 and the seed to check-base-150.json.
7. **Minor.** An untracked TOPAS case `validation/reference_cases/topas/v3-005d-de4-smoke-water-150mev-1cm/` exists. No reference run was created (`.ionmc-cache/reference-runs` is unchanged since 10-11), so it is preparation only. Plan order: D-E4 is requested only after D-E0–D-E3 are reported, so do not submit it before this record is committed.

**Verdict: ACCEPT as exploratory record, after the item-6 provenance fields are added** (record-only, no rerun). No blocking issue. Carry 4(b) (−ln F composition for D-E6), 4(d) and 4(e) (contrary evidence) into the D-E3 record and into the H3 prior before D-E4/D-E6.

---

# Round D-E5 review (table side; plus the D-E3 provenance re-check)

## A. D-E3 provenance

Met, with one gap. Every required field is present:
- `runs` (seed, n, batches, table id, mode, wall time);
- `provenance` (HEAD d5c8a80, untracked status, backend and precision, cache directory, seeds, pinned-table sha256);
- check-base-150.json with the seed and the identical idd/recoil sha256 for both tables, plus the counters. The variant npz content is committed through the table ids, which are hashes of sidecars that contain npz_sha256.

The gap: `runtime_scripts_sha256.de3_probes.py` = bea49b84… differs from the current file (5674417c…). The committed script is therefore not the one that produced the runs. Archive the as-run copy next to the record, or record that the diff touches only the summary and reuse code. Per-run run-file sha256 would be a nice-to-have.

## B. D-E5 table side

1. **Conformance to §2.2 and the D-E4 frame: OK.**
   - Per class: rate n·σ·L, P(θ_lab > θ_w) for w = 2 and 3 mm with R the water CSDA range, θ_med, and T with the 3·θ_med cut.
   - For p-p: ⟨T_slow⟩ and the fast-proton exceedance.
   - Exact two-body lab conversion. p-p keeps the faster proton as the primary and scores both. The CM cut |μ| ≤ 0.96 is asserted on the table side, and the Geant4 side applies the same CM cut.
   - The KE cut is applied to the primary.
   - Spot checks (case 3, water 150 MeV) all match D-E0 and D-E2:
     - p-p event share 0.2488 against my nσ estimate 0.248;
     - total rate 7.07e-3 per primary;
     - p+O P(>θ_2mm) 0.689, θ_med 12.45°, T 0.0597;
     - p-p ⟨T_slow⟩ 38.7 MeV, θ_med,fast 29.6°.
   - Graphite uses the water θ_w thresholds; this is declared, as the plan states.
2. **Reuse: OK.**
   - μ_CM comes from `ElasticTable.sample_mu`, the same function `steps_v5d.p7_case` uses.
   - The lab frame comes from `elastic_kin.two_body`.
   - `kin_curve_t` is a new closed-form helper (exact relativistic forward root), used only for the D-E4 classification check, where it is needed. A unit test that `kin_curve_t(E, M, cos_lab1)` reproduces `two_body` T1 would close this.
3. **Seeds: OK.** `SeedSequence(20991000 + case)`, case 0–9 = rows 14–23, with a range assertion. `default_rng(int)` (P7) and `default_rng(SeedSequence(int))` give identical streams; this is noted only. Each case uses one stream: target choice first, then the per-class uniforms.
4. **Findings that need a dated note before D-E4 data are seen** (proposed wording; I made no edits).
   - **(a) Conflict class, important.** Both kinematic-curve residuals are < 1e-3 up to θ_lab ≈ 1.8°, not 0.5°: the p+H and p+O curves differ by ΔT/E ≈ θ²(1 − 1/16).
     - Conflicts are 0.4–1.9 % of p+O events and rise with energy. Excluding them makes the Geant4/ionmc rate and exceedance ratios depend on the shape: a narrower Geant4 core produces more conflicts. That bias is of the same order as the H1 "excluded" band of ±0.05.
     - With the 16.26° CM cut, a p-p fast proton cannot have θ_lab below about 8°, so every conflict event is p+A unless a second proton is present.
     - Proposed note: *"Amendment 19 addendum (b), dated before any D-E4 output: in D-E4/D-E5 a track whose residuals to both kinematic curves are < 1e-3 is assigned to p+A unless the same event contains a secondary proton, in which case it is p-p and the 16.26° CM cut is applied using both protons. The conflict count is reported, not excluded. The D-E5 table side applies the identical rule, and its rate is reported after the KE cut and this assignment."*
   - **(b) First-minimum estimator, minor; a note is needed because it affects the report.**
     - The frozen 1° estimator picks noise dips on the piecewise-uniform table (p+C 3.5°/5.5° at 30/50 MeV).
     - The 3σ variant finds no minimum even where the table has the real J₁ zero (p+O 150 MeV, about 34°): equiprobable bins are wide where the density is low, so the dip is smeared.
     - With 1–2e4 Geant4 events per case the frozen estimator will be noise-dominated as well.
     - Proposed note: *"The first-minimum and Helm/q-slope descriptors of §2.2 are recorded for both samples but are not interpreted. The core/surface description quotes θ_med, the exceedances, T and the 8-bin angular histograms. The frozen 1° estimator is reported next to the 3σ variant, and neither is used in any rule."*
     - The Helm fits (rms log residual 0.35–1.05; s at the 0 bound in several cases) are expected for a 2-D sharp-edge J₁ shape fitted with a 3-D Helm form. Report the fit quality and flag fits at the bounds. Descriptive only, as frozen.
   - Neither finding changes a §5 threshold or a frozen signature definition. The amendment is needed only to fix the classification rule before Geant4 data exist.
5. **Git projection: OK.** The committed `de5.json` has a longest list of 8 entries (8-bin histograms of table-sampled events). The 180-bin histogram is removed. There are no experimental, PDG or table arrays. The report is rendered from the projection.
6. **`signature_ratios`: OK.** It is defined but called nowhere; the only other references are the docstring and the report text. There are no Geant4 numbers in the code or in de5.json (no "geant" key or value).

**Verdict: ACCEPT as exploratory record.** The conflict-assignment addendum 4(a) must be committed before any D-E4 output is read. 4(b) is recommended. A has one minor gap (the as-run script copy).

---

# Round D-E4 review (de4_parse.py at d2bce93; committed and full D-E4 records; smoke REF-f8950992… re-read; read-only)

1. **Parser.**
   - (a) Frames are correct:
     - W = ±√(1−U²−V²) uses the header's sign flag (the smoke has 1 backward proton);
     - θ_lab = arccos W;
     - p+A μ_CM comes from `nuclear_v10.lab_to_cm` at E0;
     - p-p uses θ_CM = 2·arctan(γ_cm·tan θ_lab) of the fast proton, which is exact for equal masses and correct for θ_CM ≤ 90°.
     - The 16.26° CM cut is applied to the fast proton, i.e. |μ| ≤ 0.96, the same as the table side.
     - The KE filter is applied in the parser. OpenTOPAS has no current-KE filter, so the smoke records every primary. The effect is the same as the frozen filter; record this as resolved verification item 2.
     - Identification binds histories, seed (case 3: 20287001 = row 17), modules, geometry and the header max KE to the frozen CASES. The ps.phsp sha256 is in each report header, so the deleted raw files are content-addressed, and a deterministic re-run from the seeds can recreate them.
   - (b) **The "exactly two on-curve records" p-p rule is not threatened by a recoil threshold.** The smoke's minimum proton KE is 0.030 MeV, so recoil protons are produced down to tens of keV. Above the 16.26° cut, T_slow ≥ 0.02·T0 (≥ 0.6 MeV at 30 MeV) in any case.
   - (c) **Important: the off-curve p-p exclusion biases the p-p rate ratio.**
     - In the smoke, 4 of 161 secondaries are off the p+H curve. They are second interactions in the target: for example event 5078 (93.78 + 56.19 = 149.96 MeV, the secondary deflected with about 0.04 MeV lost, i.e. a p+O re-scatter). Event 30531 lost about 8.6 MeV, so it is not a simple p+O elastic re-scatter; I did not resolve it.
     - In the full cases these excluded events are 2.2 % (case 0: 556/25325) to 4.6 % (case 3: 294/6370; 3.6 % in case 4) of p-p events.
     - The table side has no plural scattering, so table/Geant4 p-p rates are biased upward by about 2–5 %. The p+A off-curve fraction (0.13–0.5 %) is small.
     - **Fix:** accept a p-p event when its fast proton is on the p+H curve (it defines θ_CM), whatever the slow proton's later history. Alternatively, add the excluded fraction back to the Geant4 rate. Report both.
     - Approximate corrected p-p ratios: 150 MeV 1.21 → about 1.16; 30 MeV 0.99 → about 0.97.
   - (d) The secondary definition (parent ID ≠ 0) is confirmed on the smoke: every secondary has parent 1. p+A events with a secondary proton (4 in case 3) are excluded, which is fine.
   - (e) The conflict rule is implemented as addendum 4(a) and shared with the table side through `select_tracks`/`classify`. The measured crossover is 1.78–1.86°. Raw conflicts are 0.3 % (30 MeV) to 2.6 % (200 MeV), and only 6–25 per case are reassigned to p-p by the secondary rule.
2. **Rates.** table/Geant4 p+A = 1.005–1.013 (O) and 1.008–1.031 (C). This is expected:
   - Geant4 BGG uses the Barashenkov elastic σ below 91 GeV, and ionmc transcribes the same arrays (acceptance plan l.163–164; the V10-A (c) note on the shared σ_el).
   - The residual 1–3 % (statistically significant at 150 MeV, SE about 0.7 %) is consistent with interpolation and A^{2/3} details plus the 0.1–0.5 % selection effects.
   - **H1 meets its "excluded" signature**: within 1 ± 0.05 at all five energies, for O and C.
   - This also means D-E4 cannot test the σ level, only shape. Say so.
3. **Scientific reading against §2.2.**
   - **H2 (shape), "identified" side.** The exceedance ratios are ≥ 1.15 at 5 of 5 energies for O and C (1.29–1.67 and 1.34–3.37 at θ_2mm), and T ≥ 1.5 everywhere (2.4–26). The §2.2 identification therefore still needs only the swap-pA ≥ 50 % condition from D-E6.
   - The table core is also wider: θ_med ratio 1.07–1.38 (O) and 1.35–1.61 (C). Both core and tail are wider, so nothing here is contrary to H2. Consistent with D-E3, the tail ratio is much larger than the core ratio.
   - **H3 (p-p), important.**
     - The rate ratio after the cut is 0.99–1.01 at ≤ 50 MeV, 1.07 at 100 MeV, and 1.21 / 1.19 at 150 / 200 MeV. The shape agrees: θ_med, P and ⟨T_slow⟩ are within ±3 %.
     - This fits D-E2. table σ_pp/BGG(total) = 0.95 / 0.99 / 1.10 / 1.10 at 50 / 100 / 150 / 200 MeV (the Hale NI(150)·S(E) rule above 150 MeV). Dividing gives Geant4's fraction of p-p events above the 16.26° CM cut: 0.94 / 0.92 / 0.91 / 0.92, i.e. a 6–9 % forward share below the cut in CHIPS. Part of that apparent forward share is the off-curve exclusion of 1(c).
     - After correcting 1(c) the 150 / 200 MeV excess is still about +15 %. That is **outside 1 ± 0.10, so the first H3 "identified" clause is met, and H3 cannot meet its "excluded" signature (1 ± 0.05)**. Under §5 rule 5, "attribution established" is then impossible unless H3 is a named component.
     - **Gap in D-E6.** swap-pp exchanges only the angular CDF at unchanged σ, and the conditional `swap-sigma` is tied to the H1 (p+A) signature. The p-p rate component is therefore measured by no arm.
     - Its size, scaling D-E3 pp-only (1−F 0.045 / 0.065, gain 0.0137 / 0.0198) by the about 15 % excess: Δ(1−F) ≈ 0.007 / 0.010 (about 20–24 % of Δ_pk) and Δgain ≈ 0.002 / 0.003 (about 35 % of Δ_g). That is the p-p gain excess flagged in D-E3 review 4(e).
     - **Required before D-E6 runs:** a dated addendum, declared as written after seeing D-E4. When the H3 rate signature fires (it has), extend `swap-sigma` to H-1 (σ_pp scaled by the corrected Geant4/table ratio, interpolated in E), and evaluate `swap-all` = σ_pp + pA + pp-shape. This is a diagnostic, never a model.
   - **H4 (recoil), note.** ⟨T_recoil⟩ table/Geant4 is 3.8–4.9 at 150–200 MeV. With ionmc's recoil plateau share of 0.19 % / 0.13 %, the recoil difference is about 0.0014 / 0.0010, i.e. about 24 % / 12 % of Δ_g. That is near the 20 % H4 "excluded" bound at 150 MeV. Evaluate it in D-E6 from the recoil channel (it is coupled to pA); do not presume exclusion.
4. **CDF extraction for D-E6.**
   - The frame and the 256 equiprobable edges are sound. The p-p symmetrisation (±μ of the fast proton on [−0.96, 0.96]) matches the table convention. Saving q-edges for the q-scaled interpolation is correct.
   - **Important: the end edges at the range limits inflate the tail.**
     - The lowest p+O bin spans [−1, 0.774] at 150 MeV although the sample minimum is μ = 0.305.
     - The table's ⟨1−μ⟩ is 0.0273 against the sample's 0.0245 (+11.5 %). This inflates the swap's gain and recoil (tail-sensitive) and slightly biases Δ_g, while K and F are nearly unaffected.
     - **Fix (declare before D-E6):** set the backward end edge to the sample extreme (μ_min), or use a declared exponential extension fitted to the last k bins. Report the −1 variant as a systematic.
   - **Finite-sample error (§3.2).** Use 200 deterministic bootstrap resamples of the selected events per case and node. Rebuild the edges, apply the q-interpolation, and evaluate K*_shifted and G. Add three items in quadrature to SE_r:
     - the bootstrap error;
     - the tail-edge systematic above;
     - a leave-one-energy-out check of the q-interpolation (predict the 150 MeV CDF from 100 and 200 MeV, and compare K* and G).

     With p+O counts of 1.4e5 → 1.3e4 and p-p counts of 2.4e4 → 5.7e3, I expect the bootstrap error on F to be ≲ 1e-3. The tail-bin and interpolation terms may dominate. Optional: two transport runs on bootstrap tables at the CRN seeds (about 30 s each) to check the K*-based estimate.
5. **Git.** The committed de4 records contain no lists at all (except the 1-entry seeds) and no histogram, CDF, edge or σ keys. They hold counts, rates per primary, angles, exceedances and moments; the reports carry run ids and input sha256. Allowlist met. Note: rate per primary with the stated n and L lets a reader recover a Geant4 σ. That is an aggregate of a project-generated run and is allowed, but state it in the README. The CDF npz and the full records stay git-ignored.

**Verdict: ACCEPT as exploratory record**, with three required items before D-E6:
- (i) the p-p off-curve rate correction, 1(c);
- (ii) the p-p σ arm addendum, 3, which means §5 cannot reach "established" without it;
- (iii) the CDF tail end-edge treatment and the declared finite-sample/interpolation error estimate, 4.

H1 is excluded, and H2's D-E5 signature holds. H3 is identified on rate at 150/200 MeV, which is contrary evidence to its "low" prior; preserve it.

---

# Round D-E6 review and the §5 verdict (records at 0637fba + regenerated generated/; plan §2.2/§5; Amendment 19 addenda (b), (c))

1. **Clauses 1–3 hold at both energies.**
   - **Clause 1.** The base F differs from the archived value by +0.0013 (150 MeV) and +0.0015 (200 MeV), against 2·√(SE_new² + SE_arch²) = 0.0041 / 0.0075. The table-identity check is bit for bit.
   - **Clause 2.**
     - 150 MeV: r = +0.0014, SE_r 0.0027.
     - 200 MeV: r = +0.0051 (1.4 SE_r), SE_r 0.0037.
     - The threshold max(0.010, 2·SE_r) = 0.010 at both energies, so the clause holds, and also on −ln F (r = +0.0016 / +0.0058).
   - **Clause 3.** pA carries 1.00 / 0.99 of Δ(both), and |I| ≈ 0. Holds.
   - **SE propagation is sound.** The paired jackknife over CRN batches is used, SE_obs = 0.0026 / 0.0034 per §5, and the CDF error enters in quadrature.
   - **The CDF error (0.00021–0.00023) is plausible within a factor of about 2.**
     - My rough binomial estimate is ≈ 0.0003–0.0004. It assumes removal probability ≈ 0.4–0.5 per node, n = 1.4e5 / 1.1e5 / 4.7e4 / 2.1e4 / 1.3e4 at 30 / 50 / 100 / 150 / 200 MeV, a K* weight spread over the nodes, and K*_O ≈ 0.07.
     - The leave-one-energy-out proxy at 200 MeV (200 is an end node, so a true leave-one-out would be an extrapolation) is 1.1e-4.
     - Even a 5× underestimate (0.001) would leave SE_r dominated by SE_obs, so neither item affects clause 2.
2. **Clause 4 holds.** The §2.2 H2 "identified" signature needs swap-pA ≥ 50 % of Δ_pk, plus either exceedance ≥ 1.15 at ≥ 3 of 5 energies or T ≥ 1.5. Observed in `de4-ratios.json` (corrected rule):
   - swap-pA closes 96 % (0.0335/0.0348) and 88 % (0.0380/0.0433) of Δ_pk;
   - the exceedance ratio at 2 mm is 1.29–1.67 at 5 of 5 energies, and at 3 mm 1.34–2.28;
   - T is 2.4–25.
3. **Clause 5 FAILS.**
   - **H1 excluded:** O 1.005–1.013 and C 1.008–1.031, all within 1 ± 0.05.
   - **H5 excluded:** within ±2 %.
   - **H6 excluded:** only the modules, the seeds and the header comment line differ. The comment is non-physics; state that.
   - **H3 NOT excluded.** The column requires "ratios within 1 ± 0.05 **and** |swap-pp effect on F| < 0.005 at both energies". The p-p rate ratio is 1.165 / 1.159 at 150 / 200 MeV, so the first conjunct fails, even though swap-pp (shape) is ≈ 0.
     - H3 is not "identified" either: swap-pp closes 0 %, and even the post-hoc σ arm closes only 6 % / 13 %, against the 30 % requirement. Its state is "neither identified nor excluded". Clause 5 requires exclusion, so it fails.
   - **H4 NOT excluded at 150 MeV.** The frozen measurement is the recoil plateau change in swap-pA (the "identified" column) against the "excluded" bound "recoil-explained part < 20 % of Δ_g". The value is 25 % (+0.00146 / 0.00585) at 150 MeV and 14 % at 200 MeV, so the exclusion fails at 150 MeV as written.
     - **Contrary reading, preserve it:** this recoil change appears only in the arms that swap the p+A shape. It is the recoil part of the H2 effect (fewer large-angle p+A events give less recoil energy), not an independent bookkeeping defect. TOPAS also deposits its sub-mm recoil ions locally.
     - The frozen §2.2 text does not provide that exemption, so the reading cannot change the verdict. It is recorded beside it.
4. **Post-hoc arms (addendum (c)3).**
   - The label rule ("closure reached only with the post-hoc p-p σ arm") applies only if clause 2 fails on the frozen 2×2. Clause 2 holds, so the extended set is reported beside the verdict and changes no label.
   - The extended results: Δ(ppS) = +0.0022 / +0.0058 (6 % / 13 % of Δ(all)), r_all −0.0007 / −0.0012, I_all ≈ 0.
   - **Consistency check.** The arm scales σ_pp by 1/1.165 only above 100 MeV (the ratio is 0.98–1.03 at 30–100 MeV, interpolated linearly). For a 150 MeV beam the path-averaged scale above 100 MeV is about −6 %; with D-E3 pp-only 1−F = 0.045 and about 70 % of the p-p removal at E > 100 MeV, the estimate is ≈ +0.002, against +0.0022 measured. At 200 MeV (about −11 % over a larger share) the estimate is ≈ +0.0057, against +0.0058.
   - This is consistent with the D-E4 rate excess and with D-E2 (σ_pp = 1.10× BGG at 150/200 MeV, plus CHIPS's 6–9 % share below 16.26°). **Correction of my D-E4 review:** my "20–24 % of Δ_pk" applied the 150 MeV excess to the whole path and overestimated the effect; the measured 6 % / 13 % stands.
5. **Interpretation**, report-only; preserve the contrary evidence.
   - **The F lever is the core width, not the tail.** The Geant4 CHIPS p+A shape closes the archived V11 peak gap to 96 % / 88 % at fixed σ. With D-E3 this places the F lever in the core width:
     - pA-gauss changed the tail and σ_tr (−36 %) but left F unchanged (0–4 % closure);
     - s0.9 closed 87 % / 79 %, R×1.2 about 50 %;
     - CHIPS is narrower in the core (θ_med ratio 1.34) and far weaker in the tail (T ratio 16–25).
   - **H2a as originally stated, rev-1 "sharp-edge q⁻³ tail inflating σ_tr causes the peak deficit", is refuted for the peak.** The tail contributes to the gain and the recoil. The attribution unit "H2 (p+A shape)" stands.
   - **Contrary evidence to preserve:**
     - (i) At 200 MeV r = +0.0051, so the frozen 2×2 leaves about 12 % of Δ_pk; the ppS arm brings it to −0.0012.
     - (ii) The gain is not closed consistently: Δ_g(pA) +0.0069 against the observed +0.0059 at 150 MeV (overshoot ≈ 1σ), and Δ_g(all) +0.0081 against +0.0059 (overshoot ≈ 2.5σ). At 200 MeV Δ_g(all) is +0.0069 against +0.0083 (undershoot).
     - (iii) The peak-versus-gain ordering for pA is "tail-like" at 150 MeV and "core-like" at 200 MeV, so the descriptive core/surface reading is energy-dependent and inconclusive.
     - (iv) H4 at 150 MeV (item 3).
6. **Swap-table construction: no material upward bias of Δ(pA) found.**
   - The 30 MeV floor (ionmc's wider shape kept below 30 MeV) biases Δ(pA) down, i.e. conservatively.
   - The KE filter drops the most forward Geant4 events, which makes the CDF slightly wider and also biases down.
   - The off-curve (plural-scatter) exclusion drops 0.1–0.5 % of slightly wider events; the bias is up but negligible.
   - Sample-extreme end edges truncate the tail; this was measured as the tail-edge systematic (dF 3–9e-5, dG ≈ −1e-4), small compared with the gain overshoot.
   - The q-scaled interpolation error (leave-one-out 1.1e-4) is small.
   - σ is unchanged, by design, and the rates agree within 3 %.
   - The p-p symmetrisation matches the table convention and is irrelevant here (Δ(pp) ≈ 0).
   - Holding the 200 MeV node fixed affects nothing, because 200 MeV is the beam energy.

## §5 verdict (wording proposed for decision 0041)

> **Partially attributed (frozen §5 rule, Amendment 19).** Clauses 1–4 hold at both energies:
> - the base reproduces the archived V11 F;
> - the closure residual of the frozen 2×2 swap is r = +0.0014 (SE_r 0.0027) at 150 MeV and +0.0051 (0.0037) at 200 MeV, against the 0.010 threshold;
> - the named component **H2 (p+A shape)** carries 100 % / 99 % of Δ(both), with a negligible interaction;
> - its D-E5 signature holds (exceedance ratio ≥ 1.15 at 5 of 5 energies; T ratio 2.4–25).
>
> **Clause 5 fails:**
> - **H3 (p-p)** does not meet its exclusion column: the Geant4/table p-p rate ratio is 1.165 / 1.159 at 150 / 200 MeV, outside 1 ± 0.05.
> - **H4 (recoil)** does not meet it at 150 MeV: the recoil part is 25 % of Δ_g, against the 20 % bound. The recoil change occurs only in the shape-swap arms; this is recorded as a contrary reading and does not change the verdict.
>
> The post-hoc p-p σ arm (addendum (c)3, declared after D-E4) contributes 6 % / 13 % of Δ(all) (r_all = −0.0007 / −0.0012). It is reported beside the verdict and changes no label.
>
> The original "sharp-edge tail" (H2a) mechanism is refuted for the peak deficit. The shape lever is the core width, while the tail affects the gain and recoil. The gain closure is inconsistent: it overshoots at 150 MeV and undershoots at 200 MeV.
>
> V5, V11 and V10-A stay failed. This is exploratory, not evidence.

**Consequence to flag (scope).** §5 says only an *established* attribution authorises a model-change decision. §4 says M2 is "adopted if §5 attributes H2". With "partially attributed" and H2 named, the two texts conflict, so the frozen record does not authorise the M2 model change by itself. Resolve it conservatively: record that the rule is not satisfied, and either raise an operator proposal (adopt M2 on the partial attribution, with H3 and H4 carried as open findings) or pre-register a follow-up that resolves H3 and H4. Do not read the conflicting texts as an authorisation.

**Verdict on the D-E6 record: ACCEPT as exploratory record**, with the §5 outcome **"partially attributed"** (clause 5: H3 at both energies, H4 at 150 MeV).

**Orchestrator note (2026-10-17):** the ratio called "Geant4/table" in the Round D-E6 §5 wording is the record's `table/Geant4` ratio (ionmc table over Geant4; `de4/de4-ratios-report.md`): the table's p-p rate is 16.5 % / 15.9 % above Geant4's at 150 / 200 MeV. Decision 0041 and Amendment 19 addendum (d) use the record's direction.

---

# Round A20 review (Amendment 20 draft; worktree HEAD 3b2a240; plan of record sha256 fd9657df… verified unchanged)

1. **Prospectivity and preservation: OK, with one wording fix (minor).**
   - The draft re-labels nothing. "Partially attributed", every failed row, the tolerances and the consumed bases are kept, and the plan of record is untouched.
   - The wording fix: rule (b) is written *after* the §5 outcome was known and is tailored to it ("clauses 1–4 hold and H2 is named"). It is therefore a post-outcome rule adopted under the operator's response, not a "prospective" reading.
   - Replace "settles, prospectively" with: "a rule written after the §5 outcome was known (Amendment 19 (d)), adopted under the operator's response of 2026-10-18; it authorises the M2 adoption for this outcome and does not re-interpret §5".
2. **M2 diagnostic (c).**
   - **(a) Major: the CONTINUE band defeats the operator's purpose** ("inexpensive candidate-model diagnostics before expensive qualification"; the V6 lesson).
     - **The F band is too wide.** A point estimate c ∈ [0.50, 0.63) or (1.37, 1.50] at 150 MeV, or c ∈ [0.50, 0.74) or (1.26, 1.50] at 200 MeV, already predicts a V11 F failure, yet still allows ≈ 44 h of host time. Exact bands with the frozen V11 rule (TOST on the 90 % interval, Welch t = 1.943 / 1.860 at df 6.5 / 8.9, σ_rel 0.002885 / 0.003838 from step 16): **63.1–136.9 %** at 150 MeV and **73.6–126.4 %** at 200 MeV. The draft's 62–138 / 73–127 % uses 1.645 with σ = 0.0030 / 0.0040; correct the arithmetic and cite the rule.
     - **Gain and V5 ΔIDD are not screened.** Gain-type quantities can predictably fail V11 and V5 even with full F closure:
       - V11 plateau gain at 200 MeV. Archived ionmc 0.0267 against TOPAS 0.0184. Even a Geant4-identical p+A shape (D-E6 swap-pA Δ_g = 0.0039) leaves 0.0044, and with t·SE ≈ 0.0014 the bound is ≈ 0.0058 > 0.005.
       - V5's ΔIDD plateau-integral clause, which failed at 150 MeV (+10.35 %, CI [4.9, 15.7] %) and is marginal at 200 MeV (CI upper 9.96 %).
     - **Fix:** make CONTINUE a set of **necessary conditions** on point estimates. Pre-register all three now:
       - (i) c_E inside the exact V11 F acceptance band at both energies;
       - (ii) the predicted V11 gain difference within ±(0.005 − t·SE) at both energies, from D5 against the archived TOPAS gain;
       - (iii) a predicted V5 ΔIDD plateau integral within ±10 %, from the archived V5 "on/off" values adjusted by the D-E3-harness plateau change of M2 against base.

       State that these conditions are necessary and never sufficient, and that V5/V11 are decided only at lv5d with fresh seeds. A screen that admits predicted failures is not "incapable of being read as acceptance evidence"; it is merely uninformative.
   - **(b) Important: disclosure.** D1 and D2 use the V11 reference (archived TOPAS F) and Geant4 thin-target ratios to decide whether M2 transport is used at all. This is not tuning, because M2 has no free parameter, but it is a stopping rule on the gating quantity.
     - Plan §4 says D-E5/V11-type outcomes are "never selection criteria". The draft must state that the CONTINUE/STOP decision is a go/no-go on the single pre-committed model, not a selection among models.
     - The lv5d V11 record and decision 0042 must say that V11 at lv5d replicates D1 on independent seeds and is not blind to D1's closure estimate.
   - **(c) STOP → no M3 is consistent with plan §4.2**, which allows M3 only on construction or verification failures and never on V5/V11/V10-A/D-E5/D-E6 outcomes.
   - **(d) Seeds.** 20991100–20991109 lie in the 2099xxxx family and are unused: the D-E5 seeds are 20991000–09, D-E3/D-E6 use 20990100/01, and the CDF bootstrap 20990200; the test will assert it. D1's deliberate reuse of 20990100/01 for CRN is declared. OK.
3. **V5 fresh base (d): OK.** There is exactly one base (lv 20531004, r_index 4–7), only after CONTINUE and with M2. Tolerances, comparator, ladder and stop rule are unchanged. The original Amendment 17 (b)4 firing and the V5 failures at b40d812, 7aae5bb and dd976c5a stay recorded, and a second base needs a new amendment.
4. **Gate set (e).**
   - **Correct:** V5, V11 (added under 2B), V7-R2 in V7-R's slot (3A), V8-nuc, and V6-F (binding, executed once). Implementation gates: P6-D, P7, P5-ext, C1-ext, R1-D2, V8-LV, CI.
   - **Carried failures:** V6/V6-sens, V10-A at dd976c5a, V7-R, the dd976c5a V5/V11/V7-scan-el, lv5b, P6 and the MF6 rule. The new-model V10-A is report-only (data seen).
   - **Claims and semantics:** no halo or measured-dose claim; pass:false and conformant:false are explained as non-criteria; the six-item operational merge condition is sound.
   - **Add (minor):**
     - The M2 table is a new pinned elastic table. It must carry the Amendment 16 qualification flags (`no_negative_density_in_domain`, `o16_mt5_sigma_not_c12_copy`, `o16_mt5_spectra_independent`, `shape_normalised`). P6-D and P7 are re-evaluated on it.
     - The pin change enters `src/` and the hashed set before the code freeze. Step 6 already follows the freeze; say why.
     - The claim "p+A shape agrees with Geant4 within the V11 tolerances" is conditional on the lv5d V11 pass, and it is a model-to-model statement, not a data statement.
5. **Suites and seeds (f).**
   - The arithmetic is verified:
     - V7-R2: 20531004 + 16000 + k gives 20547004–20547051, with the reference 20547052; the 7200 fallback gives 20547004–20547068.
     - V6-F at B = 20551004: r7 20558004–20558019, r8 20559004–20559011, r9 20560004–20560019.
     - Fixture: 20505000 + 16000 + k gives 20521000–20521003, inside the rehearsal range 20505000–20525999 and distinct from the retired 20521004.
   - **The misprint is confirmed.** Route rev 6, lines 566 and 633, print "20558004 … 20558015"; with k = 0..15 the range is 20558004–20558019.
   - **Important: the cross-base reservations are incomplete.** With bases 10 000 apart, lv r_index r ≡ hr r_index r − 10 for every r (route l.565). The draft reserves hr 0, 6 and 17–19, but lv 14 (V11), 15 (P7), 18 (V10-A cross-check) and 19 (X-ENDF, if seeded) coincide with hr r_index 4, 5, 8 and 9.
     - **Fix:** state the hr5d r_index set explicitly (hr5c used 1 and 3), and reserve every hr r_index whose lv counterpart draws seeds, or state that hr5d uses only r_index 1–3. Extend the seed-family test accordingly.
6. **Order and budget (g).**
   - One heavy job at a time, each step planned ≤ 1200 s, and V6-F at the evidence SHA: OK.
   - The ≈ 44.1 h total excludes the V6-F CUDA run, whose rate after the host slowdown has not been measured. Measure it with a timing-only step at the rehearsal seeds before the sizes commit, and add it to the budget.
   - **Important (V6 lesson):** V6-F is binding with no fallback, and it runs first. Pre-declare what happens if V6-F fails: whether lv5d/hr5d still run (archived science, ≈ 37 h) or are deferred, since a merge is then impossible.
   - Its coverage-minimum clause (≥ 40 points after the beam screen, else "not established") is the most likely failure mode. State it in the risk paragraph.
7. **Decision 0042 content: complete.** It covers the M2 lineage (SAID shared with the Hale p-p lineage; densities from electron scattering), the domain caveat below 100 MeV, the §5 contradiction verbatim, the rule, and the H3/H4 open findings with their reporting duty.
   - **Minor correction:** "σ_pp unchanged so the p-p channel stays independent of TOPAS" is true for the **level**, which is Hale NI at 150 MeV. Above 150 MeV the energy dependence S(E) is the Geant4 BGG p-p formula (decision 0041; Amendment 15 (c)2), shared with TOPAS. State that split instead of a blanket "independent".

**Verdict: REVISE.** No blocking issue: nothing is re-labelled, waived or relaxed, the plan of record is untouched, and the seeds are arithmetically correct. Fix four items before commit:
- **2(a):** pre-register the CONTINUE rule as necessary conditions on the exact V11 F band, the V11 gain and the V5 ΔIDD plateau integral;
- **2(b):** disclose the go/no-go nature of the diagnostic and the non-blindness of the lv5d V11;
- **5:** complete the cross-base reservations;
- **6:** add the V6-F CUDA timing and the V6-F failure branch.

Minor wording items: 1, 4 and 7.

---

# Round A20-2 review (Amendment 20 draft revision 2)

1. **Prior items: applied correctly.**
   - 2a: the bands are recomputed with the frozen `tost` floored-df t95, and I re-derived 63.1–136.9 % and 73.6–126.4 %. The gain bound 0.00360 / 0.00365 and the V5 ΔIDD bound 4.60 % / 6.76 % follow the same arithmetic. All are point-estimate conditions, stated as necessary and not sufficient.
   - 2b: the go/no-go and not-blind disclosure is written into four places.
   - 5: hr5d uses r_index 1–3, and hr 0, 4, 5, 6, 8, 9 and 17–19 are reserved. This is complete for the lv rows that draw seeds; hr 7 corresponds to lv 17, which draws none. The route rev 7 correction is included.
   - 6: the V6-F CUDA timing seed is 20514000 (= 20505000 + 9000, inside the rehearsal range). V6-F runs first, and the failure branch defers lv5d/hr5d. Budget: 44.2 / 54.3 / 7.0 / 0 h.
   - Minors applied.
   - **One gap remains (major, the V6 lesson again): V5 peak/plateau is not screened.**
     - At dd976c5a V5 failed peak/plateau at −3.66 % (CI [−4.39, −2.94]) and −4.89 % ([−5.93, −3.85]). With the archived half-widths the allowed point estimates are ±1.28 % and ±0.96 %.
     - V5's elastic part moves with the V11 F change, so the implied closure bands are about 61–127 % at 150 MeV and **≈ 80–120 % at 200 MeV**. The 200 MeV band is narrower than the V11 band (73.6–126.4 %), so a run with c ∈ [0.74, 0.80) at 200 MeV passes CONTINUE but is predicted to fail V5.
     - **Fix:** add condition (vi): predicted V5 peak/plateau = archived residual + (F_M2/F_base − 1), within ±(0.02 − t·SE_V5) at both energies, with the archived SE and df from step 21.
2. **The p-p level re-anchoring (b′)4: BLOCKING as written.**
   - **(a) The quantities do not match.**
     - SAID SGT is the **nuclear** total, Coulomb removed. The transported p-p density above the cut is Hale **NI**: nuclear plus Coulomb-nuclear interference, with the pure Coulomb part left to the MCS.
     - "Partial > total ⇒ unphysical" compares NI(θ > 16.26°) with σ_N,tot. Positive interference can make that inequality hold without any defect.
     - At 150 MeV the interference is small, so a real excess of about +6 % over SAID × f_cut is plausible there (Hale's fit ends at 150 MeV). But the argument has to be made with like quantities.
   - **(b) f_cut is ill-defined as worded.** "The fraction of the Hale-shape elastic cross section with θ_CM > 16.26°" needs the Hale density below the cut. There the NI term contains the interference term ∝ 1/θ², which is negative below about 12 MeV and not integrable to a meaningful total. That is the reason for the cut in the first place.
     - Only the nuclear-only term of the LTP=1 reconstruction gives a well-defined fraction. But nuclear-only level × NI shape is internally inconsistent.
   - **(c) Undisclosed effect: the rule raises σ_pp at ≤ 100 MeV.** From the acquired SAID SM16 pp SGT (`tmp/claude-1000/said2/1_SGT_5_5.txt`), with SAID × 0.96 (nuclear pp is nearly isotropic) against the table:

     | E | SAID × 0.96 | table | change |
     |---|---|---|---|
     | 20 MeV | 156.5 mb | 140.1 mb | **+12 %** |
     | 50 MeV | 56.6 mb | 49.8 mb | **+14 %** |
     | 100 MeV | 30.4 mb | 27.3 mb | **+11 %** |
     | 150 MeV | 24.7 mb | 26.3 mb | −6 % |
     | 200 MeV | 23.3 mb | 25.65 mb | −9 % |

     - This is because the Hale NI includes the destructive low-energy interference and SAID's nuclear total does not.
     - The draft's size table shows only 150/200 MeV. Along a 150 or 200 MeV path, most p-p events at depth have E < 150 MeV, so the net effect probably **increases** p-p removal and gain, the opposite of the stated aim. It would also move the D-E4 p-p ratios at 30–100 MeV (0.98–1.03) to about 1.1–1.15. The (c)5 prognosis ignores this.
   - **(d) Lineage and disclosure.**
     - SAID is independent of Geant4/TOPAS **code**. It shares **experimental-data** lineage with the BGG fit (world pp data, PDG) and with Hale.
     - The search for the defect was prompted by D-E2/D-E4/D-E6, written after those outcomes were seen. The draft calls (b′)4 "prospective" and "not a Geant4-informed tuning". The anchor value is indeed independent, but the timing is post-outcome and the discovery was Geant4-informed. Say both.
   - **Fix (one of):**
     - **(A) Recommended.** Anchor on a SAID **NI-equivalent**: acquire SAID DSG with Coulomb, and set σ_pp(E) = 2π∫_{|μ|≤0.96} (dσ_SAID,full − dσ_Mott) dμ (event-counting convention). This is the like-for-like quantity of the Hale NI at every energy. Report the comparison with Hale at ≤ 150 MeV as a check, and re-derive the size table at all energies.
     - **(B)** Restrict (b′)4 to E > 150 MeV, replacing the BGG-scaled extrapolation, with a continuity rule at 150 MeV.
     - **(C)** Drop (b′)4 from Amendment 20 and keep H3 open.

     In every option the diagnostic keeps M2 alone as the decomposition, and the full-energy size table must be stated.
3. **β_NN (b′)2.**
   - **The physics is correct.** Horiuchi Eq. (10) follows from σ_el = σ_tot² (1 + α²)/(16πβ) with σ_el = σ_tot below the pion threshold. Applying it per channel is equally valid for pp (identical particles; σ_el = σ_tot below threshold) and pn. Summing the pp and pn phase-shift functions with T_p and T_n is the standard optical-limit form (Abu-Ibrahim 2008), and with ρ_n = ρ_p it adds no freedom. The arithmetic checks: pp 0.186 fm², pn 0.196 fm².
   - **The classification is wrong (major).** This is **not** a trigger-1 situation: the plan-literal construction (Horiuchi Table I, isospin-averaged β, constant below 30 MeV) *is* constructible from the acquired preprint. So "construction impossible" did not occur, and the refinement is a **declared prospective revision of the §4.1 input rule**. It also replaces isospin-averaged σ_NN and α_NN with per-channel values.
     - It is allowed under the operator's response, it has no free parameter, and it comes before any M2 shape is computed, but it must be labelled as such.
     - **Fix:** relabel it, and compute the plan-literal variant (Table-I β, averaged NN) as a report-only decomposition in the diagnostic.
4. **Densities (b′)3.** By plan §4.2 the unavailability of the pinned source (de Vries 1987) is itself a trigger-1 situation for the pinned input. Allowing a substitute is the same kind of declared prospective revision as the β refinement.
   - The draft's asymmetry (an "exception" for β, "no exception" for densities while still substituting) is inconsistent. Treat both alike: these are declared revisions of the §4.1 inputs, and trigger 1, leading to M3, is reserved for "no acceptable input under the revised rule".
   - The post-preliminary 0.03 fm tolerance is disclosed with an independent rationale, which is acceptable. Freeze it in the committed amendment before any acquired values are compared.
5. **Seeds and budget: OK.**
   - 20991100–09, 20514000, 20521000–03, V6-F 20558004–19 / 20559004–11 / 20560004–19 and V7-R2 20547004–52 (fallback 20547068) are all asserted by the extended test.
   - The reservations are complete, and the budget branches are coherent.
6. **Remaining paths to a predicted failure, or to a reading as evidence.**
   - The V5 peak/plateau gap of item 1.
   - (b′)4 as written could cause an avoidable STOP: at 200 MeV the gain condition gets harder once σ_pp rises below 150 MeV.
   - Nothing reads as acceptance evidence; the necessary-not-sufficient wording and the disclosures are adequate.
   - **Minor:**
     - (b)5(ii) "H3 stays open until the V5/V11 re-qualification is evaluated" contradicts (b)5(iv) "none is closed by a re-qualification result". Keep (iv).
     - State that P6 item (3) remains an LA150 reconstruction check, unaffected by the transported level.

**Verdict: REVISE.**
- **Blocking:** (b′)4 as specified mismatches nuclear-only and NI quantities, leaves f_cut ill-defined, and raises σ_pp by 11–14 % at ≤ 100 MeV without disclosure. Fix it by (A), (B) or (C).
- **Major:**
  - add the V5 peak/plateau screen;
  - relabel the β and density departures as declared prospective input revisions, not trigger-1 resolutions, and add the plan-literal β variant to the report-only decomposition;
  - give the post-outcome and data-lineage disclosure for (b′)4.
- **Minor:** the (b)5 contradiction and P6 item (3).

---

# Round A20-3 review (Amendment 20 draft rev 3; decision-0042 draft rev 2)

**Disclosure, which must be carried into the amendment:** to verify item 2, this reviewer computed σ_anchor from the acquired SAID DSG scans (`tmp/claude-1000/said3/d_<T>_<seg>.txt`) and the SGT tables (`said2/1_SGT_5_5.txt`) on 2026-10-18, before Amendment 20 was committed. The values are reported below. They were seen before commit.

1. **The ramp, the units/frame check and the ln-interpolation rule: sound in form, and not readable as evidence.**
   - The [0.80, 1.05] bound is a coarse unit/frame guard. With the correct integral convention (item 2) the ratios are 0.86 / 0.88 / 0.90 / 0.92 / 0.93 at 100 / 125 / 150 / 175 / 200 MeV, so it passes.
   - ln(DSG − Mott) interpolation is fine on |μ| ≤ 0.96: the difference is positive there (DSG 3.7 vs Mott 0.31 mb/sr at 20°, 150 MeV).
   - The 125→150 MeV ramp is ad hoc but has no free parameter. It is moot, given item 2(c).
2. **Mott term, domain and units.**
   - **(a) BLOCKING: factor-2 convention error.**
     - The table's partial σ_pp is the event-counting **half-sphere** integral: 2π∫₀^{0.96} dσ/dΩ dμ (Amendment 14 addendum (4); table JSON "half sphere, 2 pi, once per event"; it reproduces 49.8 / 27.3 / 26.3 mb).
     - The draft's σ_anchor = 2π∫_{|μ|≤0.96}, over both hemispheres "with the same factor 2π", is **twice** that.
     - Computed: full-sphere 55.0 / 46.6 / 45.1 mb at 100 / 150 / 200 MeV, against half-sphere 27.26 / 23.22 / 22.47 mb. The full-sphere ratios to SGT are 1.74–1.86, so the [0.80, 1.05] check would block (b′)4 and wrongly record it as a "units/frame" failure.
     - **Fix:** σ_anchor = 2π∫₀^{0.96}(DSG − σ_Mott) dμ, equivalently π∫_{|μ|≤0.96}. This is a transcription correction, not a data-informed change.
   - **(b) Mott formula sign (major, transcription).**
     - ENDF-102 eq. 6.9: σ_c = 2η²/(k²(1−μ²))·[(1+μ²)/(1−μ²) + ((−1)^{2s}/(2s+1))·cos(η ln((1+μ)/(1−μ)))]. In θ-form this is (η/2k)²[1/s⁴ + 1/c⁴ **−** cos(η ln tan²(θ/2))/(s²c²)], the standard Mott exchange term for spin-½ particles.
     - The draft's "− 2·(−1/2)·cos" gives **+cos**, which is wrong. The numerical effect on σ_anchor is +0.016…0.058 mb (≤ 0.2 %), but the frozen text must carry the correct form. Make the builder reuse the P6-tested `coul` function.
     - Units: mb/sr in the CM frame is consistent; the ratio check then confirms it.
   - **(c) Consequence under the rules as written: the continuity check fails at 125 MeV.**
     - With the corrected convention, σ_anchor = 27.26 / 24.50 / 23.22 / 22.67 / 22.47 mb at 100 / 125 / 150 / 175 / 200 MeV.
     - The table σ_pp at the same energies is 27.30 / 22.94 / 26.26 / 25.80 / 25.65 mb.
     - The ratios: −0.1 % at 100 MeV (an excellent like-for-like check of the convention); **+6.8 % at 125 MeV, which exceeds the pre-declared 3 % bound**; −11.6 % at 150 MeV; −12.4 % at 200 MeV.
     - Per the pre-declared rule, the re-anchoring is **not applied**, and **V-b (M2 alone)** becomes the model under test.
     - **Do not re-tune the bound, the ramp or the check energies after these numbers.** Keep the rules as drafted with only the transcription fixes (a) and (b), and record the outcome. Alternatively, drop (b′)4 from Amendment 20 altogether.
     - Any redesigned p-p component would now be informed by SAID-vs-Hale data. That is construction data, not validation data, but it needs its own disclosed, reviewed amendment.
   - **(d) New finding to record (relevant to H3).** The transported Hale-NI σ_pp is **non-monotonic** at 100–150 MeV: 27.5 mb (99 MeV) → minimum 22.3 mb (130 MeV) → 26.3 mb (150 MeV).
     - The S(E) extension above 150 MeV inherits the upturn at the end of Hale's fit range. The SAID interference-inclusive equivalent is monotonic.
     - Table/SAID ≈ 1.13–1.14 at 150–200 MeV, which matches the D-E4 p-p rate excess (1.165 / 1.159) in sign and size.
     - H3 stays an open finding. This explains it; it does not close it.
3. **V5 screen bands: recomputed and confirmed.**
   - 150 MeV: archived CI [−4.386, −2.944] %, half-width 0.007211 (SE 0.004130, df 16, t 1.7459), allowed a = 0.012789. With PP_T/PP_i = 1.03804, the band is c ∈ [0.8580·0.02476/0.0348, 0.8580·0.05132/0.0348] = **61.1–126.5 %**.
   - 200 MeV: half-width 0.010386 (SE 0.005482, df 7, t 1.8946), allowed a = 0.009614. With PP_T/PP_i = 1.05143, the band is **80.7–120.2 %**.
   - Intersected with the V11 bands, the combined window is **63.1–126.5 % / 80.7–120.2 %**. Correct.
4. **Densities.**
   - The radius-constrained surrogates for N-14, Al-27 and P-31 are acceptable as a declared construction. The rms check is honestly declared void, and a diffuseness sensitivity is reported instead.
   - **No effect on the validation claims:** water is H and O, V10-A is carbon, and V6-F is water.
   - The Si-28 anomaly is recorded, and the passing 2pF set is used with provenance. Ca-40 3pF follows the plan's fallback order. The 1974/1987 single-family lineage is declared.
   - The policy note (published model parameters with provenance, not nuclear-data values) is a reasonable reading. Flag it explicitly as an interpretation of decision 0041's data rule for the operator or reviewer record.
   - The post-preliminary 0.03 fm tolerance is disclosed. OK.
5. **Decision 0042: complete, given item 2.** It covers lineage (SAID, Horiuchi, the density sources, IAEA, Barashenkov, Hale), input uncertainties, the §5 verdict verbatim, the open findings and what it does not claim. Fill in the shas at commit (6292c1b5…, 8b96f976…, 683d028d…). Two fixes:
   - The title and §2 assert "with SAID-anchored p-p level". Under the rules as written, the continuity check fails (item 2(c)), so make the title and §2 conditional: "p-p level re-anchoring per Amendment 20 (b′)4 if its checks pass; otherwise the level is unchanged".
   - Record finding 2(d) and the reviewer-computation disclosure.
6. **Remaining risks.**
   - No path remains for a predicted failure to proceed: the combined F window, gain, ΔIDD and V5 peak/plateau screens are all present and necessary.
   - Preservation is intact: no re-labelling, and the original failures are carried.
   - With (b′)4 inactive, the 200 MeV gain condition (ii) is the likely STOP point. The prognosis in (c)5 already says so.

**Verdict: REVISE.**
- **Blocking:** the factor-2 hemisphere convention of σ_anchor, 2(a).
- **Major:** the Mott sign, 2(b); disclose the reviewer computation and keep the frozen rules unchanged in response; record the consequence (continuity fails at 125 MeV, so V-b is the model under test) or drop (b′)4; record finding 2(d); make the decision-0042 title and §2 conditional.
- Everything else is confirmed.

---

# Round A20-4 review (Amendment 20 rev 4 against rev 3, and decision 0042 rev 3 against rev 2; checked by diff)

1. **Half-sphere anchor: applied.** σ_anchor = 2π∫₀^0.96(DSG − σ_Mott) dμ_CM is labelled a transcription correction and is applied consistently in the amendment ((b′)4 and the units text) and in decision §2.
2. **Mott sign: applied.** The exchange term is −cos(…) in both documents, and the builder reuses the P6-tested `coul` function.
3. **Continuity outcome: recorded without retuning.**
   - The table matches my computation exactly: σ_anchor 27.26 / 24.50 / 23.22 / 22.47 mb against the table's 27.30 / 22.94 / 26.26 / 25.65 mb; −0.1 % and +6.8 %.
   - The rule, the bound and the ramp are unchanged from rev 3. The reviewer's pre-commit sighting is disclosed, with who and when.
   - V-a is withdrawn, V-b (M2 alone) is the model under test, and V-c is kept. r(E) is report-only. The units/frame ratios 0.86–0.93 resolve the unit uncertainty.
4. **H3 lead: recorded consistently in both documents.**
   - The numbers check against the table npz: 22.31 mb at 130.0 MeV, 26.26 at 150 MeV, and table/SAID 1.131 / 1.142 at 150 / 200 MeV.
   - The separate investigation, "P6-D-type check of the Hale NI integration around 130–150 MeV", is named. "No re-qualification result closes H3" appears in both documents.
5. **Diff scope.** No other change between rev 3 and rev 4 beyond these items and their direct consequences:
   - (b)2 adopts M2 alone;
   - (b)5(ii) and 6;
   - the (b′)4 headings ("not applied", "had the change been applied");
   - V10 report-only with the p-p level unchanged;
   - the decision title, Supersedes, §2, §4 and Open findings;
   - all shas and ids filled in.

   The bands, seeds, gate list, order and budget are untouched. The decision and the amendment agree on adoption, the p-p lineage (unchanged; S(E) is BGG above 150 MeV) and the diagnostic variants.
6. **Three leftover inconsistencies in the amendment** (text only; minor, but the gate definition must be exact before freezing):
   - **(a)** l.68, (b)2: "…inputs per (b′)1–3, 5), The p-p level…". Replace with "…inputs per (b′)1–3, 5). The p-p level…".
   - **(b)** l.327, (c)6: "The adoption (M2 + (b′)4) is withdrawn." Replace with "The adoption of M2 is withdrawn."
   - **(c)** l.354, (e): "P6-D (evaluated on the new table, M2 + (b′)4)". Replace with "P6-D (evaluated on the new table, M2; p-p level unchanged)".

   l.238 ("…because the re-anchoring is applied after the reconstruction") sits under "Consequences (the rule as written; not applied)" and is acceptable as conditional text.

**Verdict: REVISE (minor; text only).** Apply the three one-line edits of item 6 verbatim. With those applied, no further review round is needed, and the package (Amendment 20 rev 4 and decision 0042 rev 3) is ready to be committed as written.

# Round A20-5 review (Amendment 20 rev 6 against the A20-4 rev 4 text; decision 0042 rev 5; evidence `m2-inputs/VERIFY.md`, `NOTES.txt`; decision 0041 at worktree HEAD)

1. **Preservation of the rev 4 commitments: OK, with every deviation recorded.**
   - The model under test changed (M2 → M3 → none), with change tables for rev 5 and rev 6. V-c is withdrawn. A new report-only seed block 20991200–09 is added to the seed-family test.
   - The diagnostic (c) is retained but not run. (d) lapses with its seeds unobserved. The (e) gate set is unchanged as text; its reduction is routed to the closure amendment.
   - The verdict "partially attributed", the recorded failures, the tolerances and the bases are untouched. The (b′)4 outcome and the H3 lead are kept.
2. **Trigger 3: applied as written; no reinterpretation.**
   - From VERIFY.md (independent implementation): O-16 σ_el^M2 = 299.4 / 223.7 / 165.4 mb against Barashenkov 352 / 159 / 102 mb, giving **0.851 / 1.407 / 1.622**. C-12: 229.1 / 167.0 / 122.1 against 253 / 125 / 84, giving 0.906 / 1.336 / 1.454.
   - O is outside [0.7, 1.4] at 2 of 3 energies, so the trigger fires; C is outside at 1 of 3.
   - The verified physics cause (free-NN optical limit; Glauber σ_tot +16–19 % against EXFOR 13569/13753; α_NN = 0 removes the excess) is recorded beside the outcome without overriding it. A medium correction is correctly refused as a post-outcome model choice.
   - **Important, precision and margin:** the decisive 150 MeV value is **1.407**, only 0.5 % above the bound. Print 1.407 (not "1.41") in the amendment and the decision, and state that the outcome hinges on this margin. The trigger stands as evaluated on M2 as defined, and no variant may be used to move it.
   - **Minor:** also record as a planning defect that the plan's §4.2 rationale ("optical-limit Glauber is expected to reproduce p+A σ_el … well within [the band]") was wrong at 150–200 MeV.
   - **Minor:** trigger 4 for C-12 at 150 MeV is "not evaluable (no minimum)". Under the frozen wording ("must lie within") that is not a pass. It is moot given trigger 3; record it as such.
3. **M3 impossibility: stated as an acquisition fact, with the defect attributed. OK, two wording fixes.**
   - (b‴)1 grounds "construction impossible" on the attempt table alone. (b‴)2 attributes the unverified energy coverage to the orchestrator as a planning defect and says it is not the ground.
   - **Fix (minor):** qualify "impossible" as "**not acquirable within the experiment's acquisition tool and the routes tried**". Name the remaining route explicitly, so the closure is not read as a physics judgment: an operator-assisted acquisition (the 256 MiB cap, or manual placement of MAT 625/825). Pursuing it would itself be an operator/infrastructure matter.
   - **Fix (minor):** decision 0042 §1 heading "fallback **ill-conceived**" asserts a physics judgment that is not established, since the coverage is *unverified*. Replace it with "fallback's energy coverage unverified (planning defect)".
4. **"No model change": follows from the plan, but not from the A20-4 text of (b)3.**
   - (b)3 as accepted in A20-4 says only "Triggers 1–4 lead to M3 as written". Rev 4 (c)6 defines "no model change" only for a STOP at the diagnostic. Neither covers M3 being construction-impossible.
   - The outcome does follow from the plan of record §4: M3 is the last fallback and M4 is rejected, so no model remains. Rev 5 (b″)4 pre-declared it before the M3 acquisition result was known.
   - **Fix (minor):** in (b‴)3 replace "This follows Amendment 20 (b)3 as revised" with "This follows plan-of-record §4 (M3 is the last fallback; M4 rejected) and the rule pre-declared in (b″)4 (revision 5, written before the M3 acquisition outcome); the A20-4 text of (b)3 did not address this case". Add this path to (c)6's "no model change" sentence.
5. **The closure amendment is not self-executing: OK.** It is "defined here; it requires operator approval before it is committed", with the reason given (disposition A reserves scope changes). Withdrawing the elastic depth-dose claim only narrows claims and is a consequence of the failed gates. **Important omissions in the proposal content,** which the operator needs for an informed choice:
   - **(a) Default-off has its own recorded V5 failure.** The non-elastic V3-005A configuration that would become the default failed V5 against TOPAS at lv5b 7aae5bb: peak/plateau **+12.5 % / +12.3 %**, ΔIDD −18.5 % / −14.7 % (acceptance table l.100). Elastic-on at dd976c5a: −3.7 % / −4.9 %. Both options carry a V5 failure, and elastic-on is the closer one. The "default off" recommendation must state this comparison.
   - **(b) V6-F configuration.** V6-F, the binding gate, transports the elastic channel through the D-1 table, so its reduced halo claim would apply to the elastic-on configuration. If elastic ships default-off, the default physics has no halo gate, and the V6-F claim refers to a non-default option with a known depth-dose deficiency. The proposal must state which configuration V6-F certifies and how the claim is worded.
6. **Decision 0042 rev 5: consistent with the amendment**, verified by diff. It covers the title, "Supersedes: nothing", M2/M3 records, the attempt table, the diagnostic not run, lineage, and "not acquired". Fix three things:
   - **(a)** The heading wording of item 3.
   - **(b)** "Dated sections: (None yet. The diagnostic outcome … and the re-qualification outcome are appended …)" is stale. Replace it with "None; no diagnostic or re-qualification in this state".
   - **(c)** Decision 0041's 2026-10-18 section records the disposition "adopt M2 …; one fresh V5 base at lv 20531004". That disposition is now superseded by the outcome. Add a dated one-line pointer in 0041 (2026-10-18, appended): "superseded by decision 0042: no model change (M2 failed trigger 3; M3 not acquirable); D-1 remains". Otherwise 0041 and 0042 disagree on the current model plan.

   The 0041 2026-10-17 verdict text is quoted verbatim and correctly in 0042.

**Verdict: REVISE (minor/important text fixes; no blocking issue).**
- Print 1.407 and its margin.
- Correct the (b‴)3 citation and the (c)6 path.
- Qualify "construction impossible" with the operator-assisted route, and remove "ill-conceived".
- Add to the closure proposal the default-off V5 comparison (+12.5 % / +12.3 % against −3.7 % / −4.9 %) and the V6-F configuration statement.
- Fix the stale dated-sections line in 0042 and add the supersession pointer in 0041.

Trigger application, preservation and non-self-execution are confirmed.

---

# Round A20-6 review (verification of the A20-5 findings; Amendment 20 rev 7 against rev 6, decision 0042 rev 6 against rev 5, and the proposed decision-0041 pointer; checked by diff)

1. **1.407 and the margin: applied.** The ratio is printed as 0.851 / 1.407 / 1.622 in the change table, in (b″)1 and in decision 0042 (Options and §1b). The 0.5 % margin and the "no variant may move it" statement are present. The band-rationale planning defect is recorded, and C-12 trigger 4 at 150 MeV is "not evaluable (no minimum), not a pass, moot".
2. **"No model change" citation: applied.**
   - (b‴)3 now cites plan §4 (M3 is the last fallback; M4 rejected) and (b″)4 (rev 5, written before the acquisition outcome), and says the A20-4 (b)3 text did not cover this case.
   - (c)6 names both paths.
3. **M3 wording: applied.** The amendment (b‴)1 and decision 0042 (title, §1 heading, §1 text) say "not acquirable within the experiment's tool and the routes tried", and name the operator-assisted route (raise the 256 MiB cap, or place MAT 625/825 manually). "Ill-conceived" is replaced by "energy coverage unverified (planning defect)".
4. **Closure proposal: applied, and the numbers check against the records.**
   - Default-off, lv5b 7aae5bb (acceptance archive table, l.100): peak/plateau +12.5 / +12.3 %, ΔIDD −18.5 / −14.7 %.
   - Elastic-on, lv5c dd976c5a step 21: −3.66 / −4.89 % → −3.7 / −4.9 %; ΔIDD +10.35 / +6.72 % → +10.3 / +6.7 %.
   - The V6-F statement (it certifies elastic-on only) and the claim wording are present.
   - The comparison remains valid at the current head: R1-D2 (B) reproduces the lv5b `r1-nuc` digest, so the elastic-off path is unchanged.
   - **Minor label fix:** the default-off row says "non-elastic V3-005A physics", but lv5b at 7aae5bb is slice-B physics: V3-005A plus secondary-proton non-elastic. Write "non-elastic physics (slice B, lv5b 7aae5bb)".
5. **Decision consistency: applied.**
   - 0042 "Dated sections" now reads "None. No diagnostic or re-qualification is run in this state."
   - The proposed **decision-0041 pointer** is accurate against the records: 1.407 / 1.622, band [0.7, 1.4], "not acquirable within … the routes tried", D-1 remains, and the fresh V5 base is not taken (Amendment 20 (b‴), (d)). Appending it as a dated paragraph at the end of 0041's 2026-10-18 section is consistent with the dated-entry convention.
   - **Minor grammar fix** in 0042 Consequences (l.226): "That amendment which requires operator approval because it removes V5 and V11 … implementation gates" should read "That amendment requires operator approval because it removes V5 and V11 from the selected gate set. The remaining gates are V6-F, V7-R2, V8-nuc and the implementation gates." (The final period is missing.)

The diff shows no other changes.

**Verdict: ACCEPT.** The five A20-5 findings are correctly applied. Apply the two one-line text fixes (the label in item 4 and the grammar in item 5) verbatim; no re-review is needed. With those, Amendment 20 rev 7, decision 0042 rev 6 and the decision-0041 dated pointer are ready to be committed as written.

# Round A21-1 review (Amendment 21 draft 5; evidence M2PRIME.md "Reading of ρ₁, ρ₂" and "R-C verification"; m2prime.py, m2prime2.py, cv_mc.py)

1. **Numbers against M2PRIME.md: all match, with one exception.**
   - **Exact matches:**
     - R-C σ_el/Barashenkov: C 0.595 / 0.976 / 1.129; O 0.595 / 1.062 / 1.286.
     - σ_tot deviations: C −20.3 / −1.9 / +0.5 %; O −19.3 / +0.1 / +3.9 %.
     - q_min/q_H: C 1.540/1.912 = 0.805 and 1.566/1.912 = 0.819; O 1.354/1.661 = 0.815.
     - P(θ_lab > θ_2mm) for O: 0.528 / 0.6035, printed 0.604.
     - R-B and R-A values and the Monte Carlo table: maximum |MC − form| = 4×10⁻⁴ (0.3004 against 0.3000 at 2E_F).
   - **σ_R/Barashenkov, re-derived:**
     - C: 241.9/222, 213.1/211, 200.1/210 = **1.090 / 1.010 / 0.953**.
     - O: σ_inel(p+O) = m_t − el = 645 − 352 = 293, 438 − 159 = 279, 378 − 102 = 276, giving **1.061 / 0.984 / 0.936**.
     - The m_t − el inversion is sound: it inverts decision 0041's σ_el = σ_tot(n+A) − σ_inel(p+A), with no clipping at 100–200 MeV.
   - **The exception (minor):** (c)7 says "C 43 and O 45 nodes" for the domain E ≥ 80 MeV. On the common pinned grid the count is **43 for every target** (45 is the count above 75.2 MeV = 2E_F(O)). Fix it to 43.
2. **The reading: settled correctly.** The quoted sentences support R-B as the paper's meaning of ρ₁ and ρ₂, and they support Eq. (14) as a bound-bound (nucleus–nucleus) fit. Because the paper itself points to Clementel–Villi for nucleon–nucleus scattering, R-C is a physics-grounded choice, not a band selection. I withdraw my R-A suggestion: it contradicts the text and has no free limit.
   - The Clementel–Villi factor is correct: f = 1 − 7x/5, and 1 − 7x/5 + (2x/5)(2 − 1/x)^{5/2} on E_F < E < 2E_F. It is continuous (0.3 at 2E_F) and zero at E_F. Identifying it with the flux-weighted ⟨(|k₀ − k₂|/k₀)Q⟩, verified by MC, is the right quantity for ρσ per unit length.
   - The conventions are standard: k_F = (3π²ρ/2)^{1/3}, E_F(0.17) = 38.39 MeV.
   - **Add a limitation (minor).** The derivation is nonrelativistic, and E is the bare kinetic energy. Using k₀²/2m from the relativistic momentum instead changes x by 7–9 % at 150–200 MeV: Δf ≈ +0.026, f 0.642 → 0.668 at 150 MeV and 0.731 → 0.757 at 200 MeV. Record it as a declared convention ambiguity, not adopted. Given the margins it does not affect trigger 3.
3. **Disclosure: complete and honest.** It gives the full sequence: R-B knife edge, then the reviewer's R-A objection, then the paper check, then R-C adopted with its outcomes in view. It discloses the non-blindness at trigger level, the timing of the domain choice, and the absence of a printed source.
   - **Add (important):** the expectation that the black disk kept below 80 MeV retains part of the known shape over-correction. D-E0 shows that p+O removal at E < 80 MeV is a non-negligible part of K. The diagnostic closure is therefore expected to fall short of a full-range M2′ (V-M2′ext is construction-level only). Stating this keeps a later STOP from being read as a failure of the physics alone.
4. **Keeping trigger 3 while the outcome is known: defensible now, but say what it is.**
   - My preliminary condition is met: the reading, the conventions, the domain and the numerical settings are fixed before the builder runs.
   - The R-C margins are robust: in band at 150 and 200 MeV by ≥ 8 % (O at 200 MeV is 1.286), and out at 100 MeV by 15 % for both targets. Implementation-level differences (≲ 1–2 %) cannot flip the outcome.
   - Consequence: the trigger no longer discriminates physics. Its remaining function is an **implementation check of the builder**.
   - **Fix (important):**
     - (i) State that trigger 3 carries no evidentiary weight for M2′.
     - (ii) Add a pre-declared **builder-versus-scratch reproduction record**: all trigger quantities, σ_R and q_min, with differences reported. A difference above 2 % goes to independent defect review.
     - (iii) Tighten the "no re-tuning" clause: a defect fix needs a reproducing test that does not depend on the trigger value, **and** it must not move the builder toward the scratch values except through the identified defect.

     With these, "no re-tuning except a reviewed defect fix" is tight enough.
   - The alternative, a report-only trigger with the diagnostic alone deciding, would also be acceptable. But keeping the frozen trigger is the more conservative and better-documented path.
   - **Physics caveat to record:** σ_el is 0.595 at 100 MeV, so the model over-blocks there, and M2′ is applied from 80 MeV upward. Only the normalised shape enters transport (the level is Barashenkov), but the 80–100 MeV shape comes from the regime where R-C is known to be worst. State it.
5. **Energy domain: acceptable** as a physics choice declared before the build. It does not affect the trigger nodes and keeps low-energy transport identical to D-1.
   - **Disclose precisely:** the pinned nodes adjacent to the boundary are 78.81 and 80.00 MeV, and the transport interpolates edges linearly in E between nodes. So 78.81–80.00 MeV is a blend region of black disk and M2′, not a point jump.
   - The node at exactly 80.00 MeV is included, and it is above 2E_F(C) = 79.5 MeV, so it uses the high branch. OK.
6. **Git projection: a problem with the committed scripts (important).**
   - `m2prime.py` and `m2prime2.py` (l.52) **hard-code the measured EXFOR means**: `TOT={'C-12':{100:492.5,150:341.5,200:293.5},'O-16':{100:645.0,150:443.0,200:375.0}}`. M2PRIME.md also prints them ("measured sigma_tot = EXFOR … means").
   - Committing either as-is violates (f): measured σ_tot excluded, EXFOR by entry number only.
   - **Fix:**
     - the scripts read the EXFOR values from the cache at run time (the draft-2 wording; restore it in (f));
     - the projected M2PRIME table drops the measured-value line;
     - it must not carry **both** absolute model σ_tot and the % deviations, because together they recover the measured values. Keep the deviations and drop absolute σ_tot, or the reverse.
   - Short quotations from arXiv:1004.2096 are fine under the arXiv licence; record the licence.
7. **Carry-over from earlier drafts.**
   - **Amendment 20 references (minor):** rev 7 (c), (e) and (f)8 refer to "the model under test (M3)", "the new M3 table" and "the M3 verification of (b″)4". Amendment 21 (e) must state that on CONTINUE every such reference reads as M2′: the M2′ table, with trigger 2–4 records in place of the M3 verification, plus P6-D, P7 and the Amendment 16 flags on it.
   - **Seeds:** add D4 seeds 20991300–09 to the seed-family test.
   - Nothing else conflicts. The bands, gates, bases, budget and order are unchanged, and STOP leads to the closure route requiring operator approval.

**Verdict: REVISE (minor/important, text and record hygiene; no blocking issue).** The required fixes:
- the node count (43);
- the relativistic-convention limitation;
- the expectation for below 80 MeV and the 80–100 MeV caveat;
- the trigger framing, the builder-versus-scratch reproduction record and the tightened re-tuning clause;
- the 78.81–80 MeV blend disclosure;
- the EXFOR values removed from the scripts and the projection;
- the Amendment 20 M3 → M2′ reference mapping and the D4 seeds in the test.

The physics definition (R-C, flux-weighted, conventions pinned) and the disclosure sequence are sound.

---

# Round A21-2 review (Amendment 21 draft 7 against draft 5; installed Amendment 20 rev 7 in the worktree acceptance plan; seed test)

1. **Node count: applied.** 43 for every target on the common grid.
2. **Trigger framing: applied.**
   - "Trigger 3 carries no evidentiary weight; implementation check".
   - The builder-vs-scratch reproduction record at six cells, with 2 % on σ_tot, σ_R and σ_el leading to independent defect review, and the outcome stays the builder's.
   - The re-tuning clause is tightened ("not toward the scratch values except through the identified defect") and now also names relativistic k₀²/2m. OK.
3. **Limitations (b)5–8: applied.**
   - The nonrelativistic Δf ≈ +0.026, with my numbers quoted correctly.
   - The black disk below 80 MeV, with the shortfall expectation and the warning against reading a STOP as a physics failure.
   - The 80–100 MeV over-blocking regime.
   - The 78.81–80.00 MeV blend.
4. **Git policy: mostly applied, one remaining leak path (important).**
   - The orchestrator's rule (no absolute σ of any variant; σ_tot pass/fail per cell; relative-difference reproduction record) and the script fix (cache reader, plus a test that no EXFOR value is a literal) close the direct leaks.
   - **Remaining leak.** (f) projects σ_R/Barashenkov and σ_el/Barashenkov ratios for **all** readings, including M2, variant (a). The Barashenkov arrays are committed in `src/ionmc/nuclear/bgg.py`, so σ_inel = m_t − el and σ_el are derivable. Therefore:
     - absolute σ_R^M2 and σ_el^M2 follow from the ratios;
     - σ_tot^M2 = σ_R + σ_el follows from those;
     - with the committed M2 deviations in Amendment 20 rev 7 (b″)1 (+18 / +19 % at 150 MeV, +16 / +18 % at 200 MeV), the measured n+C / n+O σ_tot then follows to about 0.5 %.
   - This is exactly the leak the orchestrator meant to prevent.
   - **Fix:** do not project σ_R/Barashenkov for variant (a) = M2. Its σ_el ratio alone is already committed and harmless. For any other variant that carries a printed % deviation anywhere in Git, apply the same rule. R-C, R-B, R-A and the other readings carry pass/fail only, so their ratio pairs are safe.
5. **Carry-over mapping (e): applied, but the scope is too broad, and there is a seed conflict (minor).** Against the installed text (acceptance plan):
   - **(a) Bare "M3" garbles two sentences.** The mapping lists bare "M3", which would also rewrite:
     - l.1115 "plan §4.2 forbids switching to M3 on such outcomes";
     - l.1136 "No M3 follows".

     Both must stay literal: M3 is the plan fallback, and M2′ has no fallback. Restrict the mapping to "the model under test (M3)", "V-M3", "the new M3 table" / "the M3 table", "the M3 verification of (b″)4" and "(b″)4 verification". State that every other "M3" stays literal.
   - **(b) D4 seeds.** Installed (c)1 (l.1069) and the seeds table (l.1247) give V-M3 the D4 seeds 20991100 + case. Under the mapping these would become M2′'s, which conflicts with A21's 20991300 + case.
     - Add: "for M2′, D4 uses 20991300–20991309. 20991100–20991109 stay reserved and unused (V-M3, never run)."
     - The seed test registers both blocks (l.57; set check at l.87), so only the text needs the statement.
   - l.1073 ("In D1–D7, 'M2' reads as the model under test (V-M3)") chains correctly to M2′ under the restricted mapping. OK.
6. **Seed test: OK.** 20991300–09 is registered in `tests/ionmc/test_v3_005d_seed_families.py` (l.57) and included in the used-set equality (l.87). I did not re-run the tests (read-only), so the "4 pass" result is taken from the report.

**Verdict: REVISE (two small text items; no blocking).**
- 4: drop σ_R/Barashenkov for M2 from the projection.
- 5: restrict the mapping to the listed phrases, keeping l.1115/1136 "M3" literal, and add the D4 seed statement.

Everything else from A21-1 is correctly applied. With these two edits applied verbatim, no further review round is needed.

# Round A21-3 review (Amendment 21 addendum (a) draft, file mtime 09:22:25; context: d4-execution-brief.md and de7 at 42db373)

The current file content is internally consistent. I found no sign of a partial edit; the status line matches de7 (triggers 2–4 pass, reproduction within 0.16 %).

1. **Prospectivity and strictness: correct in substance, but the strictness claim is overstated, and the (iv) outcome is predictable (important).**
   - True: no transport, D4-sampler or evaluator number exists, and de7 is a trigger-level record without D4/P values.
   - **The revised (iv) is predictable from committed numbers.** Amendment 21 (a)3/(a)6 already discloses the construction-level O-16 P(θ_lab > θ_2mm): M2′ 0.528 / 0.604 at 150 / 200 MeV, table 0.688 / 0.735, CHIPS 0.513 / 0.571. M2PRIME adds R-C 0.391 at 100 MeV, against table 0.622 and P_G4 ≈ 0.404 from de4 (1.541). So M2′/G4 ≈ 0.97 / 1.03 / 1.06, against table/G4 1.54 / 1.34 / 1.29. The revised (iv) at 100/150/200 MeV is therefore **expected to pass with large margins**, before any run.
   - **Strictness, per reading:**
     - (a) Literal reading (identical tables do not "move"): the original (iv) is a guaranteed STOP, and the revision is a **relaxation** from impossible to expected-pass.
     - (b) Sampling reading: neither rule dominates. The original passes with 2/3 high energies plus both low ones (≈ 25 %), which the revision forbids. The revision passes with 3/3 high energies and no low one, which the original forbids. The revision is stricter on the informative energies only.
   - So "stricter in proportion (100 % against 80 %)" is not true in every reading.
   - **Is the change justified?** Yes as a planning-defect repair: a clause evaluated where the model is identical to the reference by construction cannot test the model. It is not tuned to an unknown outcome. But its outcome is known in expectation, so the clause no longer screens anything.
   - **Fix:** replace the strictness sentence with: "Relative to the literal reading this removes a clause that fails by construction; relative to the sampling reading it is stricter on the model-informative energies. The outcome of the revised (iv) is predictable from the construction-level P values disclosed in Amendment 21 (a)3 and M2PRIME, so (iv) carries no screening weight, like trigger 3. The decisive, blind conditions are (i), (i-b), (ii) and (iii)."
   - Equally acceptable: declare (iv) report-only. The orchestrator chooses, but the disclosure is mandatory either way.
2. **The 30/50 MeV "identical by construction" listing weakens nothing else.**
   - (i), (i-b), (ii) and (iii) are evaluated at 150/200 MeV from transport, and (v) is the construction checks. None references the 30/50 MeV D4 cases.
   - A2 (water O-16 class, P at 2 mm) is the reading that matches the cited 1.29–1.67 (`de4-ratios` H2 exceedance_ratio_2mm 1.568 / 1.670 / 1.541 / 1.343 / 1.287). It is interpretive and not tuned.
   - A3 (exact frozen-rule values; the printed bands are roundings) is correct and makes the evaluator the authority.
   - A5 (the D1b multiplicativity and D3 additivity assumptions stated) is good practice.
3. **A4 is supported by the committed text, not a scope change.** Amendment 21 (e), installed at acceptance plan l.1485: "CONTINUE: Amendment 20 (d) (one fresh V5 base, lv 20531004), (e) …, (f) and (g) apply with M2′ as the model." Reading 20 (d)2's "only with M2" as M2′ follows directly from that clause.
4. **Execution rules: sound** (one run; re-execute once only after a crash with no published document, at the same seeds, both recorded; the verdict comes from the frozen evaluator). **Add (minor):** a run that publishes a document but is invalid (the run's `valid` flag false, or a non-zero failure counter) counts as executed. The screen verdict is then STOP ("not established"), never a re-run.
5. **Harness changes: complete against the brief.** The list covers:
   - the variant loader with table id and npz sha;
   - the R-B table;
   - reuse of the de3 base with a sha256 check against de6, the brief's preferred path;
   - the de5 seed-base and table-id parameters with the guard;
   - the de0 table-id argument;
   - the de8 evaluator reading only committed records, recomputing the bands with `t95`, with the rounding test.

   **Enforceability of the freeze (important):** "committed + Codex-reviewed + tested before the first run" is enforceable only if the run proves it. **Add:**
   - every D-E8 run document and the evaluator output record `git HEAD` and a clean-tree flag, plus the sha256 of `de3_probes.py`, `de5_sampling.py`, `de0_kernels.py` and `de8_diagnostic.py`;
   - the evaluator refuses to run on a dirty tree, or when HEAD does not contain the Codex-reviewed commit named in the record;
   - the operational definitions of the brief §1 (F and gain from the 20–60 mm plateau of 0.5 mm IDD; δ_E = Σ_{20–60 mm} ΔIDD × bin_g_cm2 from batch means; elastic-only runs; paired-jackknife SE reported, not used) are implemented in the evaluator and covered by its tests, so they are frozen with it.
6. **Disclosure header: good** (go/no-go, not a selection, stopping rule, lv5d not blind, construction outcomes and de7 known, executed once). **Add one sentence:** "Condition (iv) was revised after the de7 record (addendum (a)2) to the 100/150/200 MeV energies; its outcome there is predictable from the construction-level P values disclosed in Amendment 21 (a)3 and carries no screening weight."

**Verdict: REVISE (text and enforceability; no blocking issue).**
- Item 1: correct the strictness claim and add the (iv) predictability disclosure, or make (iv) report-only.
- Item 4: add the invalid-run rule.
- Item 5: make the freeze provable (HEAD, clean tree, script sha256s, evaluator refusal; operational definitions in the evaluator tests).
- Item 6: add the header sentence.

Items 2 and 3 are confirmed.

---

# Round A21-4 review (addendum (a) rev 2 against rev 1; checked by diff)

1. **Every A21-3 item is applied; the diff shows nothing else changed.** The changes:
   - (iv) is report-only, with the predictability statement (≈ 0.97 / 1.03 / 1.06 against 1.54 / 1.34 / 1.29) and "no screening weight".
   - The decisive blind conditions are (i), (i-b), (ii) and (iii); CONTINUE = (i), (i-b), (ii), (iii), (v).
   - "Stricter in proportion" is dropped. The literal-reading relaxation and the sampling-reading statement are given.
   - The A2 wording is adjusted ("reported").
   - A4 cites l.1485.
   - The operational definitions are now in the evaluator and its tests.
   - The freeze-provenance rule text and the invalid-run rule are added.
   - The header sentence is added.
2. **Refusal conditions: checkable offline as written, but three gaps.**
   - **Offline-checkable as written:**
     - dirty tree (`git status`);
     - HEAD containing the reviewed commit (`git merge-base --is-ancestor`);
     - the four script sha256 values against that commit (`git show <c>:<path>`);
     - run-document HEAD, clean-tree flag and sha256 against the manifest.
   - **Gaps (minor):**
     - **(a) Codex pass is not checkable offline.** Only the review *id* and the reviewed commit SHA are. State that the evaluator checks the id's presence and the commit SHA, and that the "passed" status is verified by the orchestrator at freeze time and recorded in the manifest and the commit message.
     - **(b) Add a table-identity refusal.** Refuse if any run document's elastic table id or npz sha256 differs from the de7-recorded M2′ table (`addd8179…`) for V-M2′, or from the de7 M2 and R-B records for the report-only variants. This guarantees that the screened table is the one whose triggers passed.
     - **(c) Add a seed refusal.** Refuse if any run seed differs from the registered set: harness 20990100 / 20990101; D4 20991300–20991309 (V-M2′) and 20991200–20991209 (V-M2). The evaluator output should also record the sha256 of every generated run document it consumed. The de3 base documents are already sha-checked against de6.

**Verdict: ACCEPT.** The three one-line additions of item 2 are to be inserted verbatim into the freeze-provenance rule text, with no further review round. With them, addendum (a) rev 2 is ready to be committed as written.

# Round D-E8 review (the executed go/no-go screen; worktree HEAD c509e351; record de8, manifest, attempt log, run documents and telemetry; read-only)

1. **Bands: they match the frozen rule and the plan.**
   - Recomputed independently:
     - a_F = 0.02 − 1.9432·0.002885 = 0.014394 and 0.02 − 1.8595·0.003838 = 0.012863. With F_i/F_T = 0.85678/0.89158 and 0.84407/0.88737, the closure windows are **0.6312–1.3688** and **0.7364–1.2636**.
     - a_PP = 0.01279 / 0.00961; a_gain = 0.00360 / 0.00365; a_I = 0.0460 / 0.0676.
   - These equal the record's band table and round to the printed plan bands. The t95 values are those of the frozen `t95` at the floored df (1.9432 / 1.8595 / 1.7459 / 1.8946 / 1.7291 / 1.7247 / 1.7207).
2. **The conditions are evaluated as frozen.**
   - Point estimates are compared against the exact bands (A3). CONTINUE = (i), (i-b), (ii), (iii), (v). (iv) is recorded and never part of the verdict. (v) comes from de7 by sha256 302d03a4….
   - Spot checks:
     - M2′ D1b at 150 / 200 MeV: (4.2406/4.4019)·(0.88009/0.8580) − 1 = −0.0118 and (3.8286/4.0255)·(0.87166/0.8455) − 1 = −0.0195. Matches.
     - D2 at 200 MeV = 0.0083 + Δg, so Δg_M2′ = −0.00368. Consistent.
   - The verdict is **STOP**: (i), (i-b) and (ii) fail at 200 MeV, and (iii) holds. The failing conditions are listed. The leak check passed.
3. **Provenance and once-only: satisfied.**
   - HEAD c509e351 with a clean tree, in the manifest and in every run document.
   - The four script sha256 values are in the manifest and the telemetry. The reviewed SHA equals HEAD.
   - Table ids and npz sha256 match de7: M2′ addd8179 / 9adbd8a0; M2 ad0d7757 / 3750727c; R-B 2a14cdcb / f9cc2438.
   - Seeds: CRN 20990100 / 20990101; D4 20991300–09 (M2′) and 20991200–09 (M2).
   - The append-only attempt log has 10 producers, each "started" then "published-valid" once, with no re-execution. Output mtimes show a single write of de8.
   - **The refused first evaluator call is compatible with the once-only rule.** Telemetry `20261018T082242Z-de8-evaluate` shows a refusal at argument validation ("--reviewed-head-review-id REVIEW-xxxxxxxx … is required"). That happens before any input is read or any quantity is computed, and nothing was written. The single evaluation is `…082335Z-de8-evaluate-2`.
   - **Minor (record hygiene):**
     - The refusal is visible only in telemetry. Reference it in the de8 README/decision text.
     - Record the full review id REVIEW-5949c5257eaf4c158cf7b8d9f0e12fe1 next to the 8-hex binding REVIEW-5949c525 the frozen evaluator accepts.
4. **Arm configuration: correct.**
   - All six arm documents: `elastic_only`, n = 1e5, 20 batches, CRN seeds by energy, valid, HEAD c509e351. Their sha256 match the manifest.
   - The reused base and EM-only documents match the committed `de6-report.md` input list (b42751ca / 5346b0bb / 04ba3904 / e2d093a1). The evaluator refuses a mismatch.
5. **D1 at 150 MeV, rounding gap: no issue.**
   - 0.6348 against the exact lower bound 0.63124 (printed 63.1 %) passes by 0.0035, judged on the exact value per A3.
   - Statistically the margin is far below the paired SE (0.023). The pass is fragile, but the rule is a point-estimate rule as frozen.
   - The verdict is robust: STOP follows from 200 MeV, where (i) fails by 0.13 (0.604 against 0.736, about 4 SE), (i-b) by 0.010 and (ii) by 0.001.
6. **What the record supports, and what it does not** (contrary evidence preserved).
   - **Supported:**
     - (a) The trigger-passing model M2′ does **not** close the 200 MeV V11 F gap within the screen (closure 0.60). The report-only M2 arm (free Glauber, failed trigger 3) **does** pass (i) and (i-b) at both energies (1.046 / 0.843). Within these candidates, the Pauli-blocked rate reduces the F closure, and the input-sanity trigger excluded the arm that closes F better. This is contrary evidence to record, but it is not a basis to revisit trigger 3 post hoc.
     - (b) **(ii) at 200 MeV fails for every p+A shape tested:**

       | Shape | d at 200 MeV | Δg |
       |---|---|---|
       | M2′ | 0.00462 | −0.00368 |
       | M2 | 0.00463 | −0.00367 |
       | R-B | 0.00463 | – |
       | D-E6 Geant4 CHIPS swap-pA | +0.00437 | −0.00393 |

       All exceed 0.00365, and the required Δg ≤ −0.00465 is not reached even by Geant4's own shape. So the 200 MeV gain condition is **not closable by the p+A shape alone** under the frozen D2 definition. This is model-independent with respect to the p+A shape.
     - (c) Two contributors are identifiable:
       - the p-p channel: the D-E6 post-hoc σ_pp arm gave a further Δg = −0.00287 at 200 MeV (H3, open);
       - the archived ionmc gain used by D2: g_base,arch = 0.02672 ± 0.00077 against the D-E3 base at independent seeds, 0.0249 ± 0.0010, a 1.4σ difference. With the D-E3 base the same Δg would give d ≈ 0.0028, a pass.
   - **Not supported:**
     - that (ii)@200 is a defect of the harness or evaluator: D2 is computed exactly as frozen, with the archived reference by design;
     - that M2′ specifically "fails the gain": every shape does;
     - that a lv5d V11 gain at fresh seeds would fail: the archived ionmc gain is one realisation;
     - that M2 would CONTINUE: it fails (ii) at 200 MeV too;
     - any statement about the physics validity of M2 vs M2′ beyond these screen quantities.

**Verdict: ACCEPT the D-E8 record** as a correctly executed, once-only, frozen-rule screen with verdict **STOP (established)**. The two record-hygiene notes of item 3 (telemetry reference for the refused call; the full review id) go into the README or decision text. Item 6's supported and unsupported statements should be carried into the decision text as written.

---

# Round D-E8-2 review (decision-0042 D-E8 section, README paragraph, Amendment 21 dated note, intervention-closure draft; both files mtime 10:28:07; checked against `de8/de8.json` and `de8-report.md`)

**(a) Decision section, README and dated note.** Most numbers trace to de8.json:
- D1 0.63485 / 0.60410 with SE 0.02335 / 0.02967; the exact lower bound 0.6312564 and margin 0.0036; (0.73640 − 0.60410)/0.02967 = 4.46;
- D1b −0.01184 / −0.01949; δ_E −0.164 / −0.070; D7 0.248 / 0.132 (M2 0.237 / 0.131);
- D5 ratios; D4 rows; gain_base 0.02491, so d ≈ 0.02123 − 0.01842 = 0.0028.

The disclosure header is verbatim. Provenance, attempt log, refused-call disclosure and full review id are complete. The contrary evidence does not re-open the verdict. Fix these items:
1. **D6 uses the as-frozen sign-error kernel (important).** The table reports K*_literal (M2′ 0.793 / 0.704; M2 0.639 / 0.588). Amendment 19 addendum (a)1 records K*_literal as a sign transcription error and K*_shifted as the meaningful form. Report K*_shifted from the record (M2′ **0.797 / 0.726**; M2 **0.623 / 0.606**), with K*_literal only as "as frozen". State that these are arm/table ratios.
2. **Mis-sourced P values (minor).** "The de7 O-16 P values were 0.528 / 0.604" is wrong: de7 carries no P values. They are the scratch values of Amendment 21 (a)3 and M2PRIME. Correct the source.
3. **Overstatement (important; my D-E8 item 6 rule).** "The 80 MeV domain of M2′ cost closure" is not supported. M2 differs from M2′ in **two** respects: no Pauli blocking, and full energy range. The record does not separate them: V-M2′ext is construction-level only, and D5 shows different in-domain shapes (tail ratio at 100 MeV 0.107 against 0.183). Replace with: "M2 differs from M2′ in Pauli blocking and in energy range; the record does not separate their contributions; a shortfall from the 80 MeV domain was anticipated (Amendment 21 (b)6)."
4. **Consequences not "exactly as committed" (minor).** Amendment 20 (b‴)4 says the proposal **must choose** between default-off and default-on with the limitation, and that the planner recommends default-off. Write "elastic channel shipped as a flagged option; default-off recommended (the closure proposal chooses)". The budget items match the committed (b‴)4 (6.15 + CUDA 0.16; 33.9; ≤ 3.28; ≈ 43.5 h).
5. **Derived values unlabelled (minor).** The V5 peak/plateau closure windows (0.6106–1.2653 / 0.8068–1.2015) are not record fields; they come from the frozen formula and the record's allowed |q|. They re-derive correctly. Label them "derived".

**(b) Intervention-closure draft** (current content is internally consistent; no sign of a conflicting edit):
- **Correct:** category experiment-boundary; the scientific state without a third-generation recommendation; options (A), (B) and (C) with (C) as default; the honest (B) prognosis (supported: the gain clause failed for every shape, including Geant4's).
- **Correct:** the scope change is identified (V5 and V11 removed from the Option II gate set), with both directives cited (disposition A of IR-20261013; Option II of IR-20261017).
- **Correct:** no autonomous scientific choice (model, calibration) is put to the operator. The default-off / on choice is legitimately a claims and scope matter. The item-6 statements, including the unsupported ones, are carried correctly.

Fix these items:
6. **Requirements status missing (important).** The operator's 2026-10-18 response says "This does not waive any final experiment requirement". (A) removes V5 and V11 and withdraws the elastic depth-dose claim. The request must state which requirement and ledger rows would then be recorded as **not met (not waived)**: S-PHYS-PROTON-NUCLEAR and V1-MUST-005 / 013 / 014, as Amendment 20 (b‴)4 item 5 names. It must also state that merging under (A) carries them as unmet. Otherwise the operator cannot see whether (A) conflicts with the earlier "no waiver" statement.
7. **V11 history (minor).** Removing V11 reverses part of the operator-selected 2B treatment, which added V11 as a gate. Say so explicitly.
8. **"All been tried" overstated (minor).** Write: "M2 failed its input check, M3 was not acquirable, and the declared second-generation M2′ failed the screen". M3 was not tried.
9. **Budget precision (minor).** Quote the committed values: V6-F ≈ 6.15 h CPU + CUDA (≈ 0.16 planned); remaining rows ≤ 3.28 h; total ≈ 43.5 h (Amendment 20 (b‴)4 item 6), not 6.2 / 3.3.

**Verdict: REVISE (text; no blocking issue).**
- (a): items 1 and 3 are important; 2, 4 and 5 are minor.
- (b): item 6 is important; 7–9 are minor.

With these applied verbatim, no further review round is needed for either document.

# Correction (2026-10-18, reviewer) to Rounds D-E6 and D-E8

Codex REVIEW-b4584380 is correct. Two statements in my earlier rounds go beyond the record:
- the D-E3/D-E6 review text: "the gain excess cannot be closed by any p+A shape";
- Round D-E8 item 6(b): "the 200 MeV gain condition is not closable by the p+A shape alone … model-independent with respect to the p+A shape".

The record supports only results for the **four p+A shapes actually tested**: M2′, M2, R-B, and the D-E6 Geant4 CHIPS swap. It contains no extremal bound and no search over admissible p+A angular distributions.

**The historical text stays as written, for provenance. Both statements are withdrawn as model-independent conclusions and replaced by:**

> "None of the four tested p+A shapes (M2′, M2, R-B, the D-E6 Geant4 CHIPS swap) closed the 200 MeV gain gap under the frozen D2 definition (d = 0.00462 / 0.00463 / 0.00463 / +0.00437 against 0.00365). The record does not establish that no p+A shape could close it. The p-p channel (H3) and the archived ionmc gain realisation used by D2 are candidate contributors; neither is shown to be necessary."

**What depended on the withdrawn claim:**
- **Verdicts and acceptances: nothing.**
  - The D-E8 record acceptance and its STOP verdict rest on the frozen rule applied to M2′ alone.
  - The D-E6 §5 outcome ("partially attributed") rests on clause 5 (the H3 and H4 exclusion columns), not on this claim.
  - H3's status (not excluded) rests on its own §2.2 signature (the rate ratio 1.165 / 1.159).
- **The contributor statements of D-E8 item 6(c) stand as written:** p-p, and the archived ionmc gain at 1.4σ. They were already phrased as identifiable, not necessary, contributors.
- **One wording I accepted must be tempered.** In Round D-E8-2 (b) I called the intervention's option-(B) prognosis "supported (the gain clause failed for every shape, including Geant4's)". The sentence "a shape change alone is not expected to pass it" should read: "the four tested shapes, including Geant4's, all failed the 200 MeV gain clause, so M3 passing it is not expected but not excluded".
- The D-E3-round recommendation that the swap-pp arm is essential is unaffected. It follows from the tested Geant4-shape residual and from H3's rate signature.

# Round D-E9-1 review (2026-10-19; investigation-de9-plan.md; read-only)

1. **No re-opening of D-E8, and no baseline substitution: satisfied in intent, with three wording and rule fixes.**
   - D-E8 is never re-run or re-evaluated, the verdicts stand, and the closure stays deferred. The "cannot justify" list (§3) is complete.
   - **(a) The A1 SE-calibration range is wrong (important).** The plan states [0.53, 1.47] as the "√(χ²₇/7), 90 % range". The correct values are:
     - **[0.556, 1.418]** for the 90 % range (95 % would be [0.491, 1.513]);
     - for the 150 MeV control, K = 4 (df 3), a different range: **[0.342, 1.614]** for 90 %.

     Fix both numbers now, before any run. Otherwise the range becomes a post-hoc degree of freedom.
   - **(b) The report-only d(μ̂) is a substituted baseline in all but name.** Keep it only as a sensitivity number with its SE, and **without any comparison flag against the 0.00365 band**. State: "d(μ̂) is never compared with the band, never a re-evaluation, and K is not extended."
   - The "|z| ≤ 2 consistent" labels are fixed and cannot produce a pass, because they are labels only. OK.
   - **(c) The Lead A sensitivity table** (D-E3 baseline 0.00281; pooled 0.00393) carries the same risk. Label each row "sensitivity; the frozen d = 0.00462 stands".
   - **The design finding is correct.** The (ii) band's t·SE (archived combined SE 0.000783) stands in for the future run's noise. The realisation error of g_arch inside D2's point estimate (SE 0.00077) is not in the rule.
2. **g_arch provenance: verified.** V11 compare step 16 gives the ionmc plateau gain at 200 MeV as 0.026722, SE 0.000772. The seed is 20481004 + 14000 + 1 = **20495005** (V11 r_index 14, k = 1). D2 reads it through `de0h5/h5.json`. All my sensitivity arithmetic reproduces: d ± 2SE gives 0.00307–0.00616; the pooled baseline is 0.026038 ± 0.000608, giving d = 0.00393; the prediction SE is 0.00086, a 1.1-SE excess.
3. **Decomposition: reproduced exactly from the generated D-E3 documents**, whose sha256 values are listed in the committed `de3/de3-summary.md`:

   | Quantity | 150 MeV (gain, SE) | 200 MeV (gain, SE) |
   |---|---|---|
   | all | 0.02207 (0.00079) | 0.02491 (0.00099) |
   | pp-only | 0.01368 (0.00057) | 0.01982 (0.00052) |
   | pA-only | 0.00809 (0.00035) | 0.00475 (0.00046) |

   - At 200 MeV p-p is 79.6 % of the elastic gain. Additivity residual: 0.00034 at 200 MeV and 0.00030 at 150 MeV.
   - The plan discloses it as planner-recomputed. OK.
   - Note: "pp-only" is the gain with σ_pA = 0, so it includes p-p primaries' path effects and secondaries. B2 is the right step to split these.
4. **Lead B.**
   - **(a) B0: well-posed after two clarifications (important).**
     - **Independence must be specified:** "an implementation that does not import `ionmc.nuclear.elastic_build`, `endf6` or the LAW=5 reconstruction functions; it reads the raw LA150 H-1 file". Otherwise it only re-tests interpolation.
     - **"Dip persists" must be defined before computing:** "an interior minimum of the half-sphere partial between 100 and 150 MeV, lower than both endpoint values by > 1 %".
     - The 1e-3 tolerance should be stated as relative, at nodes and at node midpoints.
     - The researcher's PP-LEVEL comparison (SAID/BGG) is a separate, report-only lineage check, and does not substitute for B0.
   - **(b) B2: possible src change (important).** A generation-resolved dose tally may need a change in `src/ionmc`, which is inside the rehearsal digest. That would invalidate the V7-R fixture (a host rebuild, about 0.7 h) and need a tally-neutrality (A16-type) proof. **Fix:** B2 uses only existing species/genealogy channels (V3-004 scoring). If a src change is unavoidable, state the digest consequence and the neutrality test, and do it through a dated note before the run.
   - **(c) B3 cannot become a tuned table, provided one more rule is added.** s ∈ {0.90, 0.95} is fixed, the monotone variant uses only the table's own data, nothing is adopted, and the variant tables are loader-refused. **Add:** "no inversion: the record never reports the s, or the level, that would close the residual; the slope is reported, the closing value is not computed". The monotone variant was designed after the dip was seen; it stays a sensitivity and is never a candidate model.
5. **The "physics reason" list is not fully outcome-free (important).**
   - Criterion 3 ("B2/B3 show the p-p level is quantitatively sufficient for the residual …") makes closing the TOPAS-defined residual a reason. That is outcome-based selection. **Remove it as a reason.** Sufficiency may be reported, never as a ground.
   - **Rephrase criterion 2** as: "B0 confirms the dip is in the evaluation, and an independent-lineage level (SAID-type, deterministic comparison) disagrees with the evaluation's level beyond the combined stated uncertainties". Monotonicity alone is a shape statement, not a level discrepancy.
   - Criterion 1 (builder or interpolation defect leads to the defect route) is correct.
   - Keep the closing sentence: any σ_pp change is a new model change with its own decision, pre-registration and re-qualification.
6. **Seeds, budget, records, disclosure.**
   - **Seeds:** the block 20992000–20992199 is disjoint from every used or reserved 2099xxxx seed (20990100/01, 20990200, 20991000–09, 20991100–09, 20991200–09, 20991300–09) and from the rehearsal and base families. It must be registered in the test before any run.
   - **Budget:** under 40 min, sandbox, one job. OK.
   - **Projection:** gains, ratios, z values and flags; B0 as ratios to the 125 MeV node with monotonicity flags; no absolute σ. OK.
   - **Disclosure:** complete for the numbers seen, including the planner-computed decomposition. **Add:** my A20-3 computation of σ_anchor (SAID) and the Hale-NI node values (27.30 / 22.94 / 26.26 / 25.65 mb), which bear directly on Lead B.
   - **Amendment rather than exploratory record:**
     - any src change for B2;
     - any B0 defect (the defect route, with an intended-change record);
     - any change of the screen design learned from Lead A (archived SE in D2's band);
     - any σ_pp change.

     Everything else can stay an exploratory record.

**Verdict: REVISE (no blocking issue). Fix before any run:**
- the ρ ranges ([0.556, 1.418] for K = 8; [0.342, 1.614] for K = 4);
- d(μ̂) and the Lead A sensitivity rows carry no band flag, and K is not extended;
- B0 independence and the dip definition;
- B2 restricted to existing channels, or a declared src and digest route;
- B3 no-inversion rule;
- physics-reason criterion 3 removed and criterion 2 rephrased;
- disclosure addition.

---

# Round D-E9-2 review (2026-10-19; investigation-de9-plan.md rev 2 against rev 1, checked by diff; PP-LEVEL.md)

**Every D-E9-1 item is applied:**
- the ρ ranges [0.556, 1.418] and [0.342, 1.614];
- no band comparison for any sensitivity row, each labelled "frozen d = 0.00462 stands", and K fixed;
- B0 independence (no import of elastic_build, endf6 or LAW=5; raw LA150 H-1; eqs 6.9 / 6.10 / 6.14 anew), "match" ≤ 1e-3 relative at nodes and midpoints, and the dip definition;
- B2 on existing channels only, otherwise dropped and recorded;
- the B3 no-inversion rule, with the monotone variant never a candidate;
- the physics-reason list without the outcome criterion and with the level discrepancy beyond combined uncertainties;
- the §3 amendment list;
- ratios in place of absolute values.

My table/σ_anchor ratios are reproduced correctly: 27.30/27.26 = 1.001, 22.94/24.50 = 0.936, 26.26/23.22 = 1.131, 25.65/22.47 = 1.142.

**(a) PP-LEVEL is reported as the researcher's: OK.** The plan presents it as "the researcher's verdict … undetermined whether data or reconstruction". **Add one sentence (minor):** "This verdict is not adopted as a plan conclusion; B0 decides between evaluation feature and builder defect." PP-LEVEL itself writes "DEFECT" in capitals, while the plan correctly leaves the evaluation-feature branch open.

**(b) B0 cannot authorise a table change: OK.**
- An evaluation feature leads to "S(E) carry-up recorded as a construction choice; any change needs its own amendment (§3)".
- A builder defect leads to the defect route and an amendment.
- Neither branch changes a table by itself, and §3/§3b forbid it.

**(c) Ratio-only presentation: OK for the Git rule as far as de9 adds anything (minor note).**
- table/BGG together with the committed BGG formula (P6 item (5); `bgg.py`) recovers absolute table levels. Those are already in committed Amendment 21 (b′)4 rev 4 text, so de9 adds no new leak.
- table/SGT together with those committed table values recovers SAID SGT points. If the projection rule treats SAID solution values as restricted, project table/SGT only at the minimum and the step (0.821 and 1.021, as the plan does) and not the full row. As drafted, the plan prints only those two. OK.

**(d) Like-for-like construction: the right comparison, and convention cannot produce the step.**
- SAID DSG minus the ENDF-102 Mott term, over the event-counting half sphere |μ_CM| ≤ 0.96, is the correct like-for-like quantity for Hale's NI = σ_e − σ_c partial: nuclear plus interference, pure Coulomb removed, same cut.
- A residual mismatch between SAID's Coulomb amplitude (relativistic phases, small electromagnetic corrections) and the ENDF Mott is of order 1 % of Mott at θ ≥ 16°. Mott is about 0.3 mb/sr against a DSG of about 3.7 mb/sr at 150 MeV, so the effect is ≲ 0.3 % of the integral. It cannot produce a 20 % step.
- **The step is also convention-free in the table's own energy dependence:** table/SGT goes from 0.859 at 140 MeV to 1.021 at 150 MeV, while SGT falls.
- **Two notes to record (minor):**
  - **(i) Two like-for-like computations differ by 1–2 %.** The researcher's table/anchor values (0.979 / 0.922 / 1.119 / 1.131 at 100 / 125 / 150 / 200 MeV) are systematically 1.2–2.2 % below mine (1.001 / 0.936 / 1.131 / 1.142). The likely cause is the Mott η and k definitions, or the angle interpolation (my simple trapezoid on the 1° grid against the ln-rule). Record both, cause unresolved. This does not affect the step and does not enter B0, which tests the table reconstruction, not SAID.
  - **(ii) The deviation pattern spans 125–150 MeV, not only the 150 MeV node.** table/SAID is about 0.92 at 125–130 MeV and about 0.95 at 140 MeV, then 1.12 at 150 MeV. So "entering at the 150 MeV node" is one description; the data also show a dip below SAID at 125–140 MeV. B0's dip test covers this; the plan's text should say so.

**(e) Length:** about 2140 words is acceptable for a plan record of this scope.

**Verdict: ACCEPT.** Insert the minor sentences of (a) and (d)(i)–(ii) verbatim before commit; no further review round is needed. The plan may be committed, the seed block registered, and B0 run first.

# Round D-E9-3 review (2026-10-19; D-E9 synthesis record at reviewed HEAD 7f688e23, REVIEW-622394ed; read-only)

1. **The reading rules are implemented as frozen.**
   - **A1.** K is 8 at 200 MeV and 4 at 150 MeV, as fixed. The calibration ranges are recorded per energy ([0.556, 1.418] and [0.342, 1.614]). ρ is 1.156 at 200 MeV and 0.513 at 150 MeV, both "per-run SE consistent".
   - **z values.** z_archived is 1.13 at 200 MeV and −0.16 at 150 MeV. z for the D-E3 baseline is −0.76 and −0.19. All are labelled "consistent with the realisation distribution", and no pass/fail is attached.
   - **Lead A table.** Every row carries `compared_with_band: false` and the label "sensitivity; frozen d = 0.00462 stands". The frozen row says "the frozen value; the STOP stands".
   - **d(μ̂) at 200 MeV.** The record gives 0.003627, with SE 0.000566 in the JSON. Neither the record nor the report states the D-E8 band value 0.00365, and neither compares d(μ̂) with it. This is compliant.
     - **Wording risk (minor).** The rendered report lists d(μ̂) without its SE, although the plan requires "with its SE" and the JSON has it. Add the SE column to the report table: 0.00057 for the A1 row, 0.00087 frozen, 0.00061 pooled, "–" for the fixed-offset rows. Then no reader can treat 0.00363 as distinguishable from the frozen 0.00462 or from any band. Decision and README text must not place 0.00363 next to 0.00365.
   - **B1, B2 and B3** are report-only. B3 reports deltas and slopes only, with no inversion and no relation to the residual. The monotone variant is labelled "never candidates".
   - **B0.** The pre-declared definitions are applied: match ≤ 1e-3 relative at every node and midpoint (max 3.6e-4), and dip > 1 % below both endpoints (minimum at 130 MeV, −18.3 % and −15.1 %). The pre-declared reading follows: "evaluation feature at the LA150 top nodes".
     - The pre-run look (the parser validated at the 100 MeV node only) is disclosed in the provenance.
     - The carry-up check confirms that the S(E)/S(150) propagation is exact.

2. **The numbers are internally consistent. I recomputed these from de9.json and the D-E3 documents:**
   - **B1 shares:**
     - 200 MeV: 0.01982 / (0.01982 + 0.00475) = 0.807; additivity gap (0.02491 − 0.02457)/0.02491 = 0.0137.
     - 150 MeV: 0.01368 / 0.02177 = 0.628; gap 0.0136.
     - These are the same gains I reproduced in round D-E9-1 from the hash-bound documents.
   - **B3 slopes at 200 MeV:**
     - 0.90–0.95: (−0.00081 + 0.00193)/ln(0.95/0.90) = 0.0207.
     - 0.95–1.00: 0.00081/ln(1/0.95) = 0.0158.
   - **A1:**
     - SE_μ = 0.0011577/√8 = 0.000409.
     - z_arch = (0.026722 − 0.025732)/√(0.00077² + 0.000409²) = 1.13.
     - z_DE3 uses the D-E3 baseline's own SE: (0.02491 − 0.02573)/√(0.00099² + 0.00041²) = −0.76. That is the sensible reading of "the same way".
     - d(μ̂) = 0.025732 − 0.00369 − 0.01842 = 0.00363.
   - **B2:** the shares sum to 1. The 200 MeV secondary-proton share of 0.651 is consistent with the p-p share of 0.807 in B1.

3. **Freeze and provenance: satisfied, with two minor gaps.**
   - **Code state.** reviewed_sha equals HEAD (7f688e23), the tree is clean, and review_id is REVIEW-622394ed. Script sha256 is recorded for de3_probes, de9_pp_level, de9_record and make_variant_tables.
   - **`documents_sha256` lists 42 documents:**

     | Group | Documents |
     |---|---|
     | Production runs | 32 (A1 24, B2 2, B3 6) |
     | Reused D-E3 documents (all, emonly, pp-only, pA-only at 150 and 200 MeV) | 8 |
     | D-E8 m2prime document at 200 MeV | 1 |
     | B0 | 1 |

     The reused hashes match the committed de6 and de3 lists.
   - **`legacy_fields_absent`** is plausible. The older D-E3 documents lack code, elastic and npz fields (and pp-only and pA-only also lack backend, precision, scoring and tallies). The D-E8 m2prime document lacks only `elastic`.
   - **The three refused synthesis attempts** before the single evaluation each wrote nothing, and each was fixed in the evaluator only (dated notes 5, 6 and the addendum). This is compatible with the once-only rule, because producers were never re-run.
   - **Gap (a): the manifest is missing from the results copy (minor).** `de9-manifest.json` was listed among the record files, but during this review it is present only in `validation/generated/elastic_attribution/de9/`, not in `validation/results/elastic_attribution/de9/`. Commit it with the record, or state where it lives.
   - **Gap (b): per-point series in de9-b0.json (minor).** The Git note says "no distributions or per-point arrays", but de9-b0.json carries per-point series as keyed dicts:
     - `carry_up_check/grid_150_200_over_150`, 51 values (S(E)/S(150) from the committed BGG formula);
     - `nodes`, 20 values, and `midpoints`, 19 values (reconstruction/table ratios ≈ 1).

     The content is harmless: formula ratios and agreement ratios, with no nuclear data and no absolute cross section. But it contradicts the stated projection. Trim to the summary, the named nodes and the levels relative to 125 MeV, or reword the note.
   - **Absolute values.** The leak check passes. No absolute cross section or nuclear-data value appears in the committed copies. `levels_over_125mev` are relative levels, as the plan specified.

4. **The physics-reason assessment is the right reading of the record.**
   - **Criterion 1 (builder or interpolation defect): excluded.** B0 reproduces the table within 3.6e-4 from the raw data, and the step is present in the evaluation.
   - **Criterion 2 (level discrepancy against an independent lineage beyond the combined uncertainties): not established.** The SAID outputs carry no uncertainties (Err = 0 in the scans), and the LA150 H-1 evaluation has no covariance in hand. A combined uncertainty cannot be formed.
   - **A decision section may legitimately state, as recorded facts and without going beyond "investigative lead":**
     - (i) The 130–150 MeV step of the p-p level is a feature of the LA150 H-1 evaluation as independently reconstructed, not a builder defect. The table's S(E) carry-up propagates the 150 MeV top-node level upward, which is a recorded construction choice; any change needs its own amendment.
     - (ii) At 200 MeV the p-p channel carries about 80 % of the ionmc elastic plateau gain (B1), and p-p secondary protons carry about 65 % (B2).
     - (iii) The report-only scaled-level and monotone-level variants show that the 200 MeV gain is sensitive to the p-p level above 150 MeV: slopes 0.016–0.021, monotone variant −0.0033. No variant is a candidate, and no inversion is made.
     - (iv) A1 finds the archived baseline consistent with the realisation distribution (z = 1.13), so it is not an outlier. The design finding stands: the frozen D2 prediction carries the baseline's own SE (prediction SE 0.00087 against a 0.00097 margin), which the screen's band did not include.
   - **It must not state:**
     - that the p-p level is wrong (criterion 2 is not established);
     - that a corrected level, or the A1 mean, would pass condition (ii) or V11;
     - anything about a fresh-seed lv5d result;
     - anything implying that the STOP or a failed row is affected.
   - "STOP stands; no amendment proposed" remains correct. Any σ_pp change or screen-design change requires its own reviewed amendment.

**Verdict: REVISE (minor; record hygiene only; no blocking issue).** The three items: the SE column for the Lead A rows in the report; `de9-manifest.json` in the committed copies; the per-point dicts in de9-b0.json trimmed, or the projection note reworded. With these applied verbatim, no further review round is needed, and the record may be committed and cited in decision 0042 within the bounds of item 4.

**Disposition (orchestrator, 2026-10-19, after the round).** (1) The Lead A table of the rendered report now carries the SE column (renderer only; the Markdown body was re-rendered from the unchanged `de9.json`, whose sha256 is unchanged; the D-E9 record was not re-evaluated). (2) `de9-manifest.json` is committed with the record (allowlisted in `de9_record.GIT_FILES`). (3) The per-point dicts of `de9-b0.json` (51 carry-up ratios of the committed BGG formula; 20 node and 19 midpoint reconstruction-over-table agreement ratios) are kept as produced: the B0 record is hash-bound in the D-E9 record (once-only), and the "no per-point arrays" sentence is the common projection note written into every record, so it cannot be reworded inside the record; the discrepancy is disclosed instead in the results README and in decision 0042 (formula and agreement ratios only; no nuclear data; no absolute cross section). No synthesis statistic, seed, document or reading rule changed.

*Addendum to the disposition (2026-10-19, Codex REVIEW-80ec0344 on cd028e3f, pass with two minor findings).* (i) The pooled-baseline row's `se` is the SE of the pooled baseline alone (the delta term is not propagated; its covariance with the D-E3 baseline is not formed): the rendered table now says so, as does decision 0042; the hash-bound `de9.json` keeps its field. (ii) The B0 per-point dicts: the Git projection `de9-b0.json` is the hash-bound once-only B0 document and is not trimmed; the Git policy of the results README carries a dated reconciliation (keyed agreement-ratio and formula-ratio series). No statistic, document, seed, setting or reading rule changes.
