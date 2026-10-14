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
