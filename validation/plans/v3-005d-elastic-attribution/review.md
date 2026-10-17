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
