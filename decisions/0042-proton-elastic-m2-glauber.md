# 0042 — No elastic model change in V3-005D: M2 (optical-limit Glauber) failed §4.2 trigger 3 on verified physics; M3 (ENDF neutron mirror) not acquirable within the tool and the routes tried; the black-disk table of decision 0041 remains the model; records preserved

- Status: **accepted** (2026-10-18; revision 6, independently reviewed in rounds A20-3 to A20-6 of `validation/plans/v3-005d-elastic-attribution/review.md`; committed together with Amendment 20 revision 7 of `validation/plans/v3-005-acceptance.md`).
- Date: 2026-10-18
- Task: V3-005B (scope name V3-005D), executed on the V3-005B branch; PR #82 unmerged.
- Affects: physical accuracy (p+A and p-p elastic channels), validation strategy, data provenance.
- Supersedes: **nothing.** Decision 0041, slice C, D-1 (black-disk p+A shape, radius rule; the p-p level) remains the transport model. The records made with its table (`elastic-proton-1fd24cff…`) stand.

## Problem

The hadronic elastic channel of the dd976c5a evidence SHA over-corrects the depth dose. Measured against TOPAS `g4h-elastic_HP`:
- V5 peak/plateau is −3.7 % / −4.9 % (criterion 2 %).
- V11 F is 0.857 / 0.844 against 0.892 / 0.887.

See decision 0041, "V3-005C outcome (2026-10-13)". The V3-005D attribution (plan of record `validation/plans/v3-005d-elastic-attribution.md`, sha256 fd9657df…; Amendment 19 and addenda (a)–(d)) was run to find the cause before any model change.

## Context

**§5 verdict, verbatim** from decision 0041, section 2026-10-17:

> **Partially attributed (frozen §5 rule, Amendment 19).** Clauses 1–4 hold at both energies: the base reproduces the archived V11 F; the closure residual of the frozen 2 × 2 swap is r = +0.0014 (SE_r 0.0027) at 150 MeV and +0.0051 (0.0037) at 200 MeV against the 0.010 threshold; the named component **H2 (p+A shape)** carries 100 % / 99 % of Δ(both) with a negligible interaction; its D-E5 signature holds (exceedance ratio ≥ 1.15 at 5 of 5 energies; tail-moment ratio 2.4–25). **Clause 5 fails:** **H3 (p-p)** does not meet its exclusion column (table/Geant4 p-p rate ratio 1.165 / 1.159 (ionmc table over Geant4, as `de4/de4-ratios-report.md` defines it) at 150 / 200 MeV, outside 1 ± 0.05; shape ratios ≈ 1.00); **H4 (recoil)** does not meet it at 150 MeV (recoil part 25 % of Δ_g against the 20 % bound; the recoil change occurs only in the shape-swap arms — recorded as a contrary reading, it does not change the verdict). H1, H5 and H6 are excluded. The post-hoc p-p σ arm (addendum (c)3, declared after D-E4) contributes 6 % / 13 % of Δ(all) (r_all −0.0007 / −0.0012); it is reported beside the verdict and changes no label. The original "sharp-edge tail" mechanism (H2a) is refuted for the peak deficit: the shape lever is the core width (D-E3: `pA-gauss` leaves F unchanged, `pA-s0.9` closes most of the gap; Geant4 CHIPS has the narrower core), while the tail affects the gain and the recoil. Contrary evidence preserved: r = +0.0051 at 200 MeV; the gain closure overshoots at 150 MeV and undershoots at 200 MeV; the peak-versus-gain ordering flips between the energies. V5, V11 and V10-A at dd976c5a stay failed. Exploratory, not evidence.

**The plan-of-record §4/§5 conflict, verbatim.**
- §5: "An 'established' attribution authorises a model-change decision (§4), never a passing record."
- §4, M2 row, status column: "adopted if §5 attributes H2".
- Decision 0041 (2026-10-17) recorded: "The two texts conflict for 'partially attributed'; they are NOT read as an authorisation. No model change is made on this record."

**Operator response to IR-20261017-223039-9DA564 (2026-10-18), verbatim:** "Scientific model selection, calibration design and fresh evaluation-base allocation remain autonomous responsibilities. You may revise internally authored plans prospectively, with independent review and preservation of original failures; my earlier instruction did not reserve those technical choices to me. For the integration boundary, I select Option II, including V6-F as a binding additional gate, together with the proposed V10-A/V11 and V7-R2 treatment. This does not waive any final experiment requirement or establish full halo validation. Perform the inexpensive candidate-model diagnostics before expensive qualification, retain the stability restrictions, and merge only after the revised applicable gates and independent review pass. Formally reconcile the intervention record with this response and resolve the outstanding review finding."

**Findings that bear on the choice** (records under `validation/results/elastic_attribution/`; exploratory):
- **The peak lever is the core width of the p+A shape.** D-E3: `pA-gauss` leaves F unchanged; `pA-s0.9` closes most of the gap. D-E6: the Geant4 CHIPS shape closes 96 % / 88 % of the observed gap at fixed σ.
- **The p-p level at E ≥ 150 MeV.** Only there is the table level set by the S(E) rule.
  - D-E2 (`de2/de2-report.md`): table/BGG 0.989 at 100 MeV, 1.101 at 150 and 200 MeV.
  - D-E4 (`de4/de4-ratios-report.md`): table/Geant4 p-p rate 0.977 / 0.997 / 1.032 / 1.165 / 1.159 at 30 / 50 / 100 / 150 / 200 MeV.
  - An earlier argument, that the table partial exceeds SAID's total at 150 MeV, is **withdrawn**: the Hale partial includes Coulomb–nuclear interference, while SAID SGT is nuclear-only (Amendment 20 (b′)4, revision 3).

## Options considered

1. **No model change** (keep D-1). V5 and V11 stay failed, and integration fails on V5 under the Option II gate set. Not chosen: the attribution names the shape, and the operator asked for candidate-model diagnostics.
2. **M2, the pre-committed optical-limit Glauber shape (plan §4.1, with declared prospective revisions).** Built and independently verified. It **failed §4.2 trigger 3** (O-16 σ_el^M2/σ_el^Barashenkov 1.407 / 1.622 at 150 / 200 MeV, outside [0.7, 1.4]), so it is **not adopted**; it is kept as a report-only sensitivity. The p-p re-anchoring of Amendment 20 (b′)4 was not applied by its pre-declared check.
3. **M2 with σ_pp scaled to the Geant4 rate** (the D-E6 post-hoc arm). Rejected: it would give the p-p level TOPAS/Geant4 lineage and amount to tuning to the reference.
4. **M3** (ENDF/B-VIII.0 neutron O-16 / C-12 MF4/MT2, isospin mirror). The plan's fallback, entered when M2 failed trigger 3. **Construction impossible:** the files (C-12 MAT 625, O-16 MAT 825) could not be acquired through the tool. **Planning defect:** the frozen row did not verify energy coverage; coverage to 135–200 MeV is unverified.
5. **M4** (Geant4 CHIPS shape). Rejected for gating: shared with TOPAS (plan §4). The D-E6 swap tables are diagnostics, never a model.
6. **Resolve H3 / H4 before any model change** (IR option 0B). Not chosen. H3 and H4 are carried as open findings (Amendment 20 (b)5). H3 now has a named lead and a separate investigation (Open findings).

## Decision

**Decision: no elastic model change.** Under the post-outcome rule of Amendment 20 (b), M2 was to be adopted. It was built, independently verified, and **failed plan §4.2 trigger 3** on verified physics (Amendment 20 (b″)). The fallback **M3 is construction-impossible**, and the frozen fallback row had not verified its energy coverage (Amendment 20 (b‴)1–2). The plan has no further fallback, so **the D-1 black-disk table of decision 0041 remains the transport model**, with the p-p level unchanged. **V5, V11 and V10-A at dd976c5a stay failed.** The go/no-go diagnostic is not run. The attempt is preserved as a record: the M2 construction record, the M3 definition, and the builder paths `--model m2` / `--model m3`; no production table is built. A second-generation candidate (M2′) is the subject of the separate Amendment 21 draft.

### 1. M3 for O-16 and C-12, attempted: not acquirable within the tool and the routes tried; fallback's energy coverage unverified (planning defect) (plan §4, M3 = (ii′); Amendment 20 (b″)4, (b‴)1–2)
- **Source.** ENDF/B-VIII.0 neutron sublibrary, n-O16 and n-C12, MF4/MT2 elastic angular distributions. They are used as the isospin-mirror shape for p+O and p+C at the same T_lab.
- **Conversion.**
  - Legendre, tabulated or mixed (LTT = 1 / 2 / 3), as given; lab frame converted to CM.
  - Evaluated on the pinned table energy grid with the ENDF interpolation law.
  - 256 equiprobable μ_CM bins per node.
  - Above the file's highest energy E_top, the E_top shape is held fixed in momentum transfer q.
- **Unchanged:** σ_el is Barashenkov; Coulomb stays in the MCS; the domains [E_min,shape, 250] MeV.
- **Other targets** (N, Al, Si, P, Ca) keep the D-1 black disk. The plan defines M3 for O and C only.
- **Verification before any transport:**
  - MF1 identity;
  - copy check: O-16 is not a C-12 copy above 20 MeV, relative difference > 1e-3 at > 50 % of the points at > 50 % of the shared energies;
  - normalisation ≤ 1e-6;
  - positivity at every node and midpoint;
  - monotone CDF;
  - E_top recorded.
- **Identification.** NNDC list https://www.nndc.bnl.gov/endf-b8.0/lists/neutrons.list: C-12 MAT 625 (LANL/ORNL, EVAL-AUG15, Hale, Young, Fu); O-16 MAT 825 (LANL, EVAL-DEC16, Hale, Paris, Young, Chadwick).
- **Acquisition: failed.** Verbatim attempt table:

  | Route | Result |
  |---|---|
  | NNDC neutron zip, 310 158 456 bytes (HEAD 200), ×2 | tool error, MAX_BYTES 256 MiB |
  | IAEA `n_0825_8-O-16.zip` | tool error; curl 403 challenge |
  | LANL T2 `…/ENDFB-VIII.0-neutron/O/16`, `/C/12` (guessed URLs) | tool error; curl SSL handshake failure |
  | NNDC Sigma | ≤ VII.1 only |
  | NNDC endf-releases | whole zips only |
  | OECD-NEA JANIS | session-based (HTTP 500 / 404) |
  | NNDC per-file guesses | 404 |

  "Construction impossible" rests on this acquisition failure alone, in the sense of not acquirable within the experiment's tool and the routes tried. The remaining route is operator-assisted acquisition (raising the 256 MiB cap, or manual placement of MAT 625 / 825). It is an operator or infrastructure matter and was not pursued.
- **Planning defect (author: the orchestrator).** The frozen row (ii′) did not verify energy coverage. Coverage of these evaluations to 135–200 MeV is **unverified**: open sources give only library-wide statements ("at least 20 MeV", IAEA-NDS-227; "0 < E < 20 MeV (higher for some)", Paris & Hale, CSEWG 2018).

### 1b. M2, optical-limit Glauber: built, failed §4.2 trigger 3; record only (no production table) (plan §4.1 with the Amendment 20 (b′) refinements)
- **Outcome.** σ_el^M2/σ_el^Barashenkov = C-12 0.906 / 1.336 / 1.454 and O-16 0.851 / 1.407 / 1.622 at 100 / 150 / 200 MeV (Si 0.78 / 1.34 / 1.71, Ca 0.79 / 1.26 / 1.58, rms-matched unfolding). **Trigger 3 fires for O.** The decisive 150 MeV value 1.407 is only 0.5 % above the bound. Plan §4.2's "well within" rationale for the band was wrong at 150–200 MeV (planning defect). An independent implementation reproduces the builder; the line-by-line review found no physics discrepancy.
- **Verified cause** (recorded; the outcome is unchanged):
  - free-NN optical-limit Glauber overpredicts the measured n+natC / n+natO σ_tot by +18 / +19 % at 150 MeV and +16 / +18 % at 200 MeV (EXFOR 13569, 13753);
  - Barashenkov matches those data within 3 %;
  - α_NN = 0 removes most of the excess;
  - the σ_R values agree with Abu-Ibrahim 2008, Table III, within 1 %.
- **No medium-corrected α_NN** is introduced.
- **Trigger-4 record.** HO Helm match R_H² = √350·c₁, s_pt² = s_ch² − r_p²/3. q_H = 1.912 fm⁻¹ (C-12) and 1.661 fm⁻¹ (O-16). O-16 q_min at 150 MeV is 0.83 q_H; C-12 at 150 MeV: "not evaluable (no minimum)", which is not a pass under the frozen wording and is moot given trigger 3; C-12 at 200 MeV is 0.83 q_H.
- **Fermi unfolding.** rms-matched (Shukla 2003, p. 7), validated against the exact unfolding at 100 MeV (σ_el −1.6 / −1.0 %, q_min +0.5 / 0.0 %).
- **Nodes below the SAID range** keep the black disk.
- **The construction, as specified:**
- **Amplitude.** F(q) = ik ∫ b db J₀(qb)[1 − exp(iχ(b))].
  - k = p_CM of the p+A system; q = 2k sin(θ_CM/2); dσ/dΩ_CM = |F|².
  - The **profile phases of pp and pn are summed separately** (Abu-Ibrahim et al., Phys. Rev. C 77 (2008) 034607, Eqs. (4)–(5)):

    χ(b) = ½ Σ_{c ∈ {pp, pn}} N_c σ_c(i + α_c) ∫ T_A(b′) f_c(b − b′) d²b′

    with the Gaussian NN profile f_c(b) = exp(−b²/(2β_c))/(2πβ_c). This is exact in the optical limit with ρ_n = (N/Z)ρ_p.
- **NN inputs.** σ_c and α_c come from SAID SM16 per channel at every grid node, interpolated linearly in T_lab. There is no hold-constant rule: the acquired 5–300 MeV range covers every node.
- **β_c** comes from the unitarity rule of Horiuchi et al. 2007, Eq. (10), β_c = (1 + α_c²) σ_c/(16π), applied to the SAID per-channel values. Check at 150 MeV: β_pp 0.186 fm² and β_pn 0.196 fm², against Abu-Ibrahim 2008 ≈ 0.18 / 0.196. This is a **declared prospective revision** of plan §4.1, not a §4.2 trigger (Amendment 20 (b′)2). **The plan-literal variant** (Horiuchi's isospin-averaged Table I, held constant below 30 MeV) is built and reported alongside as a sensitivity.
- **Densities** (declared prospective revision of plan §4.1 item 3; the plan-literal de Vries 1987 source is not acquirable because of a paywall). HO form for A ≤ 16 and 2pF for A ≥ 27 (fallback order 3pF → Fourier–Bessel → sum of Gaussians), with parameters from the open substitute sources of Amendment 20 (b′)3.
  - Each target's rms charge radius must agree with the IAEA charge-radius compilation to **≤ 0.03 fm**. The tolerance was chosen after a preliminary O-16 estimate was seen; disclosed in Amendment 20.
  - The parameters are **not** written here. They live in the builder's provenance JSON, keyed to the acquisition ids below (Amendment 20 (b′)3 policy note).
  - The 1974 and 1987 de Jager–de Vries compilations are treated as one compilation family.
  - Per target (Δ = substitute rms − IAEA rms):
    - C-12: HO, Shukla 2003, Δ −0.007 fm.
    - O-16: HO, Shukla 2003 + TGlauberMC, Δ +0.019 fm.
    - Si-28: 2pF, Mehndiratta & Shukla 2017, Δ +0.016 fm. Recorded anomaly: the de-Vries-1987-attributed 3pF and FB sets fail at −0.047 / −0.038 fm.
    - Ca-40: 3pF by the fallback order, TGlauberMC, Δ +0.004 fm.
    - N-14, Al-27, P-31: IAEA-radius-constrained surrogates (shape parameter from the neighbours, size fixed to the IAEA rms 2.5582 / 3.0610 / 3.1889 fm). The rms check is void by construction, so a diffuseness sensitivity is reported instead. These targets are not used in V5 / V11 (water) or V10-A (carbon).
- **Point density.** F_pt(q) = F_ch(q)/exp(−q² r_p²/6), with r_p = 0.8414 fm (CODATA 2018). The neutron charge radius and the Darwin–Foldy term are neglected (declared).
- **Level and normalisation.** The shape |F|² is normalised to unit integral over μ_CM and scaled to the **Barashenkov σ_el(p+A)**, which is unchanged (D-1 level; shared with TOPAS). Coulomb stays in the MCS model, with no Coulomb–nuclear interference. Domain [E_min,shape, 250] MeV per target, unchanged (Amendment 15 (b)).
- **Domain caveat.** The optical limit weakens below about 100 MeV: Pauli blocking and in-medium NN amplitudes are not modelled (plan §4.3).

### 2. p-p elastic level: unchanged (Amendment 20 (b′)4 not applied by its pre-declared check)
- **Pre-registered rule** (Amendment 20 (b′)4, revisions 3–4). At nodes E ≥ 150 MeV the level would have been

  σ_anchor(E) = 2π ∫₀^0.96 [DSG_SAID(E, μ) − σ_Mott(E, μ)] dμ_CM

  - half sphere, each event counted once (Amendment 14 addendum (4));
  - σ_Mott from the P6-tested `coul` function, with exchange term −cos(η ln tan²(θ/2))/(sin²(θ/2) cos²(θ/2));
  - a linear ramp of r = σ_anchor/σ_table from 125 MeV.

  It was **blocked if** |σ_anchor/σ_table − 1| > 3 % at 100 or 125 MeV.
- **Outcome.** Computed by the independent reviewer from the acquired scans, before commit:

  | T_lab | σ_anchor | table | continuity |
  |---|---|---|---|
  | 100 MeV | 27.26 mb | 27.30 mb | −0.1 % |
  | 125 MeV | 24.50 mb | 22.94 mb | **+6.8 %** |
  | 150 MeV | 23.22 mb | 26.26 mb | – |
  | 200 MeV | 22.47 mb | 25.65 mb | – |

  - The rule and the 3 % bound were fixed before this computation and are not changed. **The re-anchoring is not applied.**
  - The transported p-p level stays the table's: Hale NI ≤ 150 MeV, S(E) rule above, P6 item (3) unchanged.
  - The builder reports r(E) at every scan energy, **report-only**.
  - The units/frame ratios σ_anchor/SGT = 0.86–0.93 confirm the mb/sr CM assumption of the DSG scans.
- **Lineage of the p-p level:** unchanged from decision 0041 (level Hale NI ≤ 150 MeV; S(E) from the BGG formula above).

### 3. Table, qualification and re-evaluation
- 256 equiprobable μ_CM bins per node, the same table schema and loader, and a new table id. The builder records every acquisition id and the β, f_cut and density provenance.
- The Amendment 16 qualification flags are required: `no_negative_density_in_domain`, `o16_mt5_sigma_not_c12_copy`, `o16_mt5_spectra_independent`, `shape_normalised`.
- The plan §4.2 construction checks are evaluated once, before any transport: normalisation ≤ 1e-6; convergence ≤ 1e-3; σ_el^M2/σ_Barashenkov within [0.7, 1.4]; first minimum within [0.6, 1.0] × q_H.
- P6-D (CI) and P7 (lv5d r_index 15) are re-run on the new table. V10 (p-p level against PDG) is re-evaluated, report-only, and its stronger shared lineage (SAID/PDG) is stated.

### 4. Go/no-go diagnostic (Amendment 20 (c)): not run (no model under test); retained for Amendment 21
- **Pre-registered variants:**
  - **V-a:** withdrawn, because (b′)4 is not applied.
  - **V-M3: M3 for O-16 / C-12, p-p level unchanged — the model under test.** CONTINUE/STOP applies to it.
  - **V-M2:** M2, report-only sensitivity (D4 seeds 20991200 + case).
  - V-c: withdrawn.
- **CONTINUE iff all of these hold**, each necessary and none sufficient (Amendment 20 (c)2–3):
  - the V11 F closure lies within the exact frozen-rule band (63.1–136.9 % at 150 MeV, 73.6–126.4 % at 200 MeV), **and** within the V5 peak/plateau band (61.1–126.5 % at 150 MeV, 80.7–120.2 % at 200 MeV); effective window 63.1–126.5 % / 80.7–120.2 %;
  - the predicted V11 gain difference is ≤ 0.00360 / 0.00365;
  - the predicted V5 ΔIDD plateau integral is ≤ 4.60 % / 6.76 %;
  - the exceedance ratio moves toward 1 at ≥ 4 of 5 energies;
  - the §4.2 checks pass.
- In this state the screen has no model to test and is not run. Its bands and rule are retained for a later candidate (Amendment 21).
- **Not blind:** the screen is a go/no-go on this single model, not a selection among models. The lv5d V5 / V11 are not blind to its estimates; this is disclosed in every record.

## Lineage and provenance (all acquisitions role `construction`)

| Item | Acquisition id | Role in M2 |
|---|---|---|
| SAID start page (GWU DAC) | 6292c1b5e8ba4f40b213fcca6d6f1d0f | provenance of the service |
| SAID SM16, T_lab 10–300 MeV, step 5: pp SGT, np SGT, pp ALFA, np ALFA | efe6d363f28c423582df8d68d2efebd0, 6035179d45c046ad93445457ea48ee7d, 23ca28089b974a32b3ccd40d68665ecc, 761b99787a674a69a3478768b1cb9385 | first acquisition, kept on record |
| SAID SM16, T_lab 5–300 MeV, step 5 (same order) | 28e78745bfc74da98de99c76035f87ea, d4360764abe742c98e488b42c0c7b90a, c39d4c94600e455c98ed374dd0c17114, dab63152cbf64d8c97d0048c20987777 | σ_NN, α_NN; p-p level |
| SAID SM16, T_lab 5–30 MeV, step 1 (same order) | 932daa93b61b416a97455f92f44abbd7, 07803ab6aab2449ca940480325fa3655, b0206879d01941b2b0cb4d4d5a30fd34, 84a77e24e0444f47987b700899e684cb | low-energy nodes |
| Horiuchi, Suzuki, Abu-Ibrahim, Kohama, arXiv:nucl-th/0612029 v1 (Phys. Rev. C 75 (2007) 044607) | f56c546a11ae4dacafa35f6a614892e3 | Eq. (10) unitarity rule; Table I for the plan-literal sensitivity |
| SAID SM16 pp DSG, T_lab 100 (θ 1–60 / 61–120 / 121–179°) | e8973bddcd3e4c1d8142033cf99d883c / cb0dfcecb64e4aaaa5829a03ac349650 / 0ad01128a2404b37be349aba6a7ca269 | continuity check |
| SAID SM16 pp DSG, 125 MeV | 7639d5dc7c2f4baf9491cee446e1574f / a6ea1a216a334db3a94edd7a45dc52ac / 12797d52222349d997efffd0768eb96c | continuity check |
| SAID SM16 pp DSG, 150 MeV | 9b95bac9bc1247dab43156074583a5c1 / 3eadb73fcf07484eb392afd8d4ba19c2 / f406ab8af2aa4c4abd11238219dcb4ee | r(E), report-only |
| SAID SM16 pp DSG, 175 MeV | a3ece8cd3cfc42d6bb0aaa01193afe2d / 8909ef36e88f49578bfe5ca5f5c60968 / 14df69fa85c846a2999238d945ea7bca | r(E), report-only |
| SAID SM16 pp DSG, 200 MeV | c2c3d72963f24e518630bf33f870ae7b / e1b69d9c30a3490e9783a97c380c84c6 / 2ed37a4234d548099155ce9e9a4f8a48 | r(E), report-only |
| SAID SM16 pp DSG, 250 MeV | 0952913571f14abf84e492b4dbd724e2 / f7de116fbc3543d481200aa2ced8de21 / dbc89fb0adf248faa73c47790c82d83c | r(E), report-only |
| SAID SM16 pp DSG, 300 MeV | 537907c89e034dc7a47804fa9d4710e1 / e38519d84ca74a78ab317eac3afd74ba / 411c5c0ae6544649aec4263669473d87 | r(E), report-only |
| Shukla | dcb08bbcc4aa446e8808e1abfa5600e5 | density source (content to be stated with the researcher's final list) |
| TGlauberMC paper | 0d4d56435c234f5892e6a76bf5998393 | density-parameter reference |
| TGlauberMC 3.3.2b tarball (GPL; sha256 8b96f976b88c648e881fea400e0c775c3d157b311b72744df8a1359f4a62a110) | cdb9f72849a84e53a2b44e1e79f0b85f | **provenance only; no code reuse** |
| Mehndiratta & Shukla, arXiv:1702.07459 (sha256 d4187fe937228963ce19a5ec784fd98eabc1034182a7eefff639982b6a077e3e; arXiv licence) | a3d9b475e79948b8b4f2f5f2de7c9d82 | Si-28 2pF (de Jager 1974 lineage); Ca-40 3pF second source |
| Duda, Fourier–Bessel parameters | bf78ff8d0fc941f39cc79046a3dcab8a | density source (FB fallback form) |
| IAEA charge-radii CSV (sha256 683d028dae93235376ddf7a345f736561ee7c0dba070091912c034bdf740eee6) | 8dd779140e52459f81f506442382b4a8 | independent rms cross-check (≤ 0.03 fm) |
| de Vries–De Jager–De Vries 1987, DOI stub | 8e31f0ce1bf6466ebe2a128a10549571 | **explicitly NOT a source** (paywalled compilation, not acquired) |
| r_p = 0.8414 fm | CODATA 2018 (recorded constant) | proton-size unfolding |
| σ_el(p+A) | Geant4 11.4.2 Barashenkov data (decision 0041, unchanged) | p+A level (shared with TOPAS) |
| ENDF/B-VIII.0 neutron n-O16, n-C12 (MF4/MT2) | **not acquired** (construction impossible; Amendment 20 (b‴)1) | M3 (attempted) |
| EXFOR 13569 (Finlay 1993), 13753 (Abfalterer 2001), acquired (Amendment 14 (l)) | (existing acquisitions) | verification of the trigger-3 cause, report-only |
| Hale NI shape | ENDF/B-VIII.0 LA150 H-1 (decision 0041, unchanged) | p-p angular shape above the cut |

All ids and hashes are complete. The ENDF/B-VIII.0 neutron files were not acquired.

**Lineage statements:**
- **NN inputs.** SAID fits the world NN data, which includes the pre-1998 pp data of Hale's fit, so it is shared with the p-p shape lineage. It is not Geant4 or TOPAS lineage.
- **Densities.** From electron-scattering fits. Independent of the p+A elastic data and of TOPAS.
- **p+A level.** Shared with TOPAS: V5 and V11 test shape, p-p and implementation, not the p+A normalisation.
- **V10-A.** Data seen, so the new-model V10-A is report-only.

## Input uncertainties (declared)
- **SAID ALFA convention.** SAID does not state its Coulomb-removal and sign conventions in its own words. The sign is consistent with Ray 1979 and Abu-Ibrahim 2008.
- **Horiuchi erratum.** Phys. Rev. C 76 (2007) 039903 and the published Table I were not verifiable; arXiv v1 only.
- **SAID DSG units and angle interpolation.** The units were not stated in the output. **Resolved:** the half-sphere ratios σ_anchor/SGT = 0.86–0.93 lie inside [0.80, 1.05] and confirm mb/sr CM. The scans are used only for the report-only r(E). No Coulomb-only option exists on the SAID form, so the Mott term is analytic (P6-tested `coul`).
- **SAID header date.** The output header `SM16 0-3.0 GEV PP=50276/25348 NP=22832/12938 RAA (147)  5/03/07` carries a stale date field. It is unresolved and recorded verbatim.
- **Optical limit below about 100 MeV.** See the domain caveat in Decision §1.

## Open findings (Amendment 20 (b)5)
- **H3 (p-p).** The p-p level is unchanged. The p-p rate ratio against Geant4 is reported for the transported table.
  - **Lead, recorded 2026-10-18:** the table's Hale-NI partial σ_pp is non-monotonic, with a minimum of 22.3 mb at 130 MeV and 26.3 mb at 150 MeV.
  - The S(E) step above 150 MeV inherits the upturn, giving table/SAID ≈ 1.13 at 150–200 MeV. That matches the D-E4 p-p excess of 1.165 / 1.159.
  - So the lead is an anomaly in the table's own p-p level between 130 and 150 MeV (builder or Hale-NI evaluation).
  - It is to be examined by a separate pre-registered investigation, **"P6-D-type check of the Hale NI integration around 130–150 MeV"**, and reported in every re-qualification record.
  - **No re-qualification result closes H3.**
- **H4 (recoil).** The `elastic_recoil` deposit and its share of the plateau gain are reported per V5/V11 run.
- Neither finding is a gate.

## What this decision does not claim
- No validation of any new model; none is adopted. No statement about the p-p level beyond the recorded H3 lead. The diagnostic is not evidence.
- V5, V11 and V10-A at dd976c5a stay failed. V5 and V11 can pass only in the lv5d re-qualification at lv 20531004 under the frozen tolerances, and the new-model V10-A is report-only.
- No lateral-halo or measured-dose claim beyond the reduced V6-F scope (Amendment 18; operator 2026-10-18: "This does not … establish full halo validation").
- No clinical claim.

## Consequences
- A new elastic table, a builder extension (Glauber amplitude, SAID readers, the report-only r(E)) and tests are needed.
- The table pin enters `src/` and the hashed source set before the code freeze, so the rehearsal fixture is regenerated afterwards (Amendment 20 (f)).
- **No model change.** V5 and V11 cannot pass in this state, and the elastic depth-dose claim is withdrawn (route §6 S1/S2). Integration under Option II needs the closure amendment of Amendment 20 (b‴)4. The default-off configuration carries its own V5 failure (+12.5 % / +12.3 % peak/plateau at lv5b 7aae5bb, against −3.7 % / −4.9 % elastic-on), and the binding V6-F certifies the elastic-on configuration only. That amendment requires operator approval because it removes V5 and V11 from the selected gate set. The remaining gates are V6-F, V7-R2, V8-nuc and the implementation gates.

## Dated sections

None. No diagnostic or re-qualification is run in this state.
