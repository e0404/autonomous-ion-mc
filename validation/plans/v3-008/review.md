# V3-008 independent review log

# Round V3-008-1 review (2026-11-01; decision-0043 draft r2 and v3-008-acceptance draft r2; read-only; verbatim)

**Round V3-008-1 verdict: REVISE.** The drafts are scientifically sound and nothing in them is tuned to TOPAS or chosen from seen numbers. Five gaps should be fixed before freeze: a seed-sequence collision with the declared V3-005 continuation; an unquantified error term in δ_pred; P4 smoke runs that could expose V1/V2 data; full-physics TOPAS runs that exceed the 1500 s per-call limit; one lineage claim stated as fact before P2 confirms it. Nothing is blocking. I modified no file.

**1. Seeds (decision (i)): clear today, but not against the declared V3-005 continuation (important).**
- Checked against the registered families and clear: `tests/ionmc/test_v3_005d_seed_families.py`: CONSUMED_BASES 20261004–20491004; FRESH_BASES 20511004–20551004 with an envelope of base + 0…20999 (maximum 20572003); REHEARSAL 20505000–20525999; TOPAS 20270000 + 1000·row + k for rows 1–23 (maximum 20293003); 2098xxxx and 2099xxxx. develop c790df36, which carries no seed in 2056–2063xxxx. The fallback values: lv8 20581004–20587103 (base + 100·r ≤ 60 + k < 100); hr8 20591004–20597103; rehearsals 20600000–20609999; TOPAS 20619001–20619040; re-runs from 20621004 in steps of 10 000.
- The conflict. Committed V3-005 Amendment 20 (f)2 declares that "the sequence continues at 20561004". Under its own envelope rule, the next V3-005 bases 20561004 and 20571004 reserve up to 20582003 and 20592003. These overlap lv8 r_index 0–9 (20581004–20582003) and hr8 r_index 0–9.
- Required: one seed-base ledger shared by all tasks. Add a dated note on the V3-005 branch: "the V3-005 sequence continues at 20641004 (skipping 20561004–20632003, reserved for V3-008)". Make the V3-008 guard test assert disjointness from V3-005's declared next bases including their 21 000 envelope. Alternatively, shift V3-008 above 20592003.
- Minor. The TOPAS seeds 20619001–20619040 coincide numerically with V3-005's TOPAS formula at row 349. They are engine seeds and scientifically irrelevant, but the guard should assert disjointness from rows 1–23 and should reserve the formula range.

**2. Scientific soundness (question (a)): sound.**
- Stopping. The analytic Bethe default with z-explicit Barkas, Bloch and first-order Mott terms is sound. So is renaming the effective charge to the Barkas form (coefficient 125), with no numeric change. Omitting LS with a bound, checked by P3 against the 0.2 % range threshold, is appropriate. The optional ASTAR path for He is correctly classed as construction when selected. Rejecting the Geant4 construction as default, to keep V1 informative, is correct.
- MCS (d) and its selection rule. The physics (the collision-number change from z² in Ω and the Born screening factor) is real and theory-derived. The pre-declared order with a 3 % criterion on θ against a pre-defined Molière–Hanson reference, validated first on Gottschalk's proton tables to ≤ 1 %, is legitimate. Declare the near-tautology: (d) is built to approximate Molière–Hanson, so I11 is calibration, as the draft already says. Add: freeze c_Ω and F_scr as committed formulas before I11 is evaluated, with no adjustment after it.
- Straggling. Bohr with bare z and the existing sampler; charge-exchange and non-Bohr terms bounded. Fine.
- Cutoff and domains. The 2 MeV/u cutoff with residual ranges ≤ 0.075 mm, energy conserved, is fine. The domains are fine.
- Nothing is tuned to TOPAS or chosen from seen numbers. The He I2 bounds are derived from the SEEN decision-0038 maxima, and they are correctly labelled as a regression, not evidence.

**3. Acceptance plan (question (b)).**
- (a) Rows and tolerances. These are defensible from theory and source arguments, frozen before results, with unambiguous PASS/FAIL/STOP. They include V1 (R80 ≤ 0.35 % and ≤ 1.0 mm; width 10 %; peak 0.5 mm; peak/entrance 3 %), I3 (1.5 %/0.8 % rms), I9 (2 %) and I11 (3 %). STOP 3 correctly consumes V1/V2 after any model change and moves gating to V1-held at a fresh base. Add: V2 has no held-out analogue. If V2 fails and the model changes, the record must state that no independent lateral MC evidence remains, which is a gap.
- (b) The V2 δ_pred scheme is legitimate but incomplete (important). It is legitimate because it is fixed before any TOPAS output, derived from source formulas, and never recomputed (STOP 3). That closes the direct tuning path. It is incomplete because a Fermi–Eyges integration of the Urban θ0 is only an approximation of the Urban algorithm (step limitation, central-plus-tail sampling, lateral displacement). δ_pred carries a method error that the 5 % band ignores. Required wording (new pre-freeze step P2b): "P2b: the same δ_pred computation is applied to protons (z = 1) and compared, before any ion TOPAS run, with the archived V3-003 proton TOPAS lateral result (T15, already seen; a method check, not a model choice). The absolute discrepancy e_δ is recorded. Case B grades |ρ − 1| ≤ 5 % only if e_δ ≤ 2 %; if e_δ > 2 %, V2 falls to Case C (informative) for every ion whose |δ_pred| > 2 %." Also freeze in P2: the exact Urban formula, version and parameters used (θ0 expression, the y definition, Z-dependent c1/c2), committed with the source hash.
- (c) V3 is correctly declared. Sihver 1998 is the only candidate, P6 is a bounded access check, and failure is an EVIDENCE GAP leading to an intervention. He has no candidate and O none beyond Sihver; both are declared. The nuclear allowance is a tolerance addition, not a correction, which is correct. The I-value calibration lineage is declared.
- (d) The He V1 non-independence from I2 is adequately declared (shared ICRU 90 alpha array). Add one sentence: "For helium, V1 tests transport mechanics and straggling only; its range part is bounded by I2 and adds no independent stopping evidence."

**4. Lineage (question (c)): correct per the research facts, subject to P2.**
- Change "FACT" to "FACT per research note; confirmed by P2" for the GenericIon claims: LS ion model; ICRU 90 proton array × q_eff² + LS ΔL + Barkas above about 2 MeV/u; ICRU 73/90 ion tables only ≤ 2.5 MeV/u.
- P2 must confirm the use of `G4IonICRU73Data` inside `G4LindhardSorensenIonModel` (the energy limit at which ion tables are consulted) and record it in decision 0043 before freeze.
- If P2 finds the ion tables used above 2.5 MeV/u, I3's "largely independent of TOPAS" class becomes "shared" and V1's C/O lineage changes. Pre-declare that consequence.

**5. Binding rules (question (d)).**
- (a) Host-call limit (important). The full-physics TOPAS carbon throughput from V3-010 is about 38 ms/history. The V3 full-physics runs (opt4 + QGSP_BIC_HP + QMD, 1e5 primaries) would take about 3800 s per call, above the 1500 s limit. Required: "every REF case is split into parts of ≤ 1200 s planned (P4 throughput), each a separate host call with its own seed (seed table extended by part index)". The EM-only runs must be sized the same way once P4 measures throughput.
- (b) P4 smoke leakage (important). P4 smokes produce ion depth-dose and lateral output before freeze, which would unblind V1/V2. Required: "P4 smoke cases carry no dose, IDD or phase-space scorers (or their outputs are deleted unread by the runner); only the run log, EM-model printout and wall time are recorded".
- (c) Other rules are met: G4EMLOW goes to the cache only (cap 1 GiB) and is never committed, with only aggregate verdicts in Git; single process, one heavy job; a new worktree from v3/develop, with no v1/v2 reads; no clinical claim; the "materials.py theory rows only" domain statement.

**6. Pre-freeze checks P1–P8 (question (e)): add P2b, fix the order, add P9.**
- Order: P8 (repository facts, A16 baseline) and P7 (seeds, with the ledger of item 1) first. Then P2 and P2b (sources, Urban formula, δ_pred method error). Then P1, P3 and P5/P5b. Then P6. P4 (scorer-free smokes) last.
- Add P9: freeze record. The committed hashes of the P1, P2 and P3 scripts and their outputs go into decision 0043, so the freeze commit is provably after them, following the D-E8 freeze pattern.

**7. Implementation scope and STOP rules (question (f)).**
- 5–8 days is optimistic but achievable. It covers charge and mass plumbing through the python reference, the Warp kernels and CUDA, the (d) correction row, the loaders, the TOPAS bundles and the lv8/hr8 suites.
- Add a schedule checkpoint: "if implementation exceeds 10 working days, record the cause and re-plan; no gate changes".
- The STOP rules are adequate: R2 bit-identical; P1 selecting no candidate; V1/V2 failure handling with no tuning, δ_pred never recomputed, consumed bases preserved; evidence gaps leading to interventions.
- Add: "T12-ion and T1-ion failures are STOP (implementation defect) before any V-row is run".

**Exact wording changes required:**
1. Seeds (§C): add the shared ledger and the V3-005 continuation note (item 1). Make the guard assert disjointness from V3-005's declared next bases with their envelope, and from the TOPAS formula range.
2. §P: add P2b (item 3b) with the e_δ rule, freezing the Urban formula and parameters with the source hash. Add P9 (freeze record). Reorder as in item 6.
3. V2 row: append "Case B requires e_δ ≤ 2 % (P2b); otherwise Case C" and "V2 has no held-out analogue; a model change after V2 leaves lateral MC evidence as a declared gap".
4. P4: the scorer-free smoke wording (item 5b).
5. §B and §E: REF cases split into ≤ 1200 s parts with per-part seeds (item 5a).
6. Decision §6 and §B lineage: "FACT (research note), confirmed by P2" for the GenericIon claims, with the pre-declared consequence if ion tables are used above 2.5 MeV/u.
7. Decision §3: freeze c_Ω and F_scr before I11.
8. §D: T1-ion/T12-ion failure is STOP before V-rows; plus the 10-day schedule checkpoint.
9. V1 He: the transport-mechanics-only sentence (item 3d).

**Disposition (orchestrator, 2026-11-01).** All nine wording changes applied in draft r3 of decision 0043 and the acceptance plan (planner, string-checked). The seed ledger of finding 1 is committed on the V3-005B branch (06974d07: the V3-005 sequence continues at 20641004; 20572004–20632003 reserved for V3-008, the first value above the existing base envelopes, which the worker found overlapped the originally proposed 20561004 start). The drafts are committed as DRAFT r3; the pre-freeze checks P1–P9 follow; the freeze (revision 1) is a separate reviewed commit.

# Round V3-008-2 review (2026-11-01; decision-0043 r4 and v3-008-acceptance r4 against the pre-freeze records at 6e155b56; read-only; verbatim)

**Round V3-008-2 verdict: REVISE.** The r4 direction is sound: rest the MCS model on the qualified proton f_dM plus a theory-derived ion/proton relative mapping, and stop using an absolute water reference of about ±3 % uncertainty as a selector. But five things must change before freeze: the labelling of what was decided after seen numbers; two specification gaps in the exact-Ω mapping; the plan is left with no gating theory check of the charge dependence; the V2 premise depends on an unread P2 file; the order of the remaining pre-freeze checks. I modified no file.

**1. (a) Is the exact-Ω adoption a priori? Partly. The labelling must be exact (important).**
- Legitimately a priori: the model *structure*. The proton f_dM is qualified, theory supplies the z/A dependence through Ω, and the exact root replaces a linearisation the derivation itself introduced. r3 ranked (d) first, and the exact form is (d) without approximation.
- Not a priori, two things: (i) The exact-Ω variant was named after the first candidate table was seen. P1 OPEN 5 appears in the same record as the no-Fano table, which showed (d) failing marginally at x/R 0.1 through the linearisation. It was named before the P1b (Fano) tables, not "before the tables were seen". (ii) The change of selection rule (from "first candidate within 3 % of Molière–Hanson selects and gates" to "a priori form; reference informative") was made after the r3 rule reached STOP 2 on seen numbers. The reasons given are legitimate, but they arose together with the seen tables: Pb fails the frozen ≤ 1 % validation, water is unvalidated, the compound Fano rule is derived, not sourced.
- Required wording (decision §3, replacing the heading clause "declared a priori" and the bullet "Not selected by the seen tables"): "The model form is the exact version of r3's first-ranked candidate (d); the exact root was named in P1 OPEN 5 after the first (no-Fano) candidate table had shown the linearisation error, and before the P1b tables. **The selection rule was revised after the r3 rule reached STOP 2 on seen P1/P1b numbers**: the Molière–Hanson water reference failed its own pre-declared proton validation for Pb (recorded as FAILED, preserved), is not validated for water, and rests on a derived compound rule, so it is demoted from selector and gate to an informative report. The seen tables were not used to choose among forms; the predictability of the exact-Ω result against the reference from the seen decomposition is disclosed."
- Also record P1's ≤ 1.0 % proton validation as FAILED for Pb, not only as OPEN. Mark the rejection argument for (a) ("near-pass is an artefact …") as interpretation; the rejection itself rests on structure.

**2. Two specification gaps in the exact-Ω mapping (important).**
- (a) Ω or b. The Fano Δb cancels only at equal velocity. pv_eq lies at a different velocity than the ion's pv, so Δb(β_eq) ≠ Δb(β). State whether the root matches Ω (Fano excluded) or b = ln Ω + Δb (Fano included at the respective velocities). Recommended: match b. f_dM encodes the proton's full b. Stating Ω is acceptable only if it is declared and its consequence is bounded.
- (b) Dimensionality. Ω_p(pv, x/R) depends on depth fraction, and so on p1v1 as well as pv. The draft calls the mapping "a local function of (E, material, z, M)" and plans a 1-D row on the E grid. Either specify a 2-D dependence (pv, p1v1) and how it is tabulated or solved per step, keeping step independence; or state and bound the approximation that makes it 1-D (for example Ω evaluated at a fixed reference x/R). Commit both choices in the formula, script and test before P9.

**3. (b) Demoting I11 leaves no gating theory check of the charge dependence. Required minimum.**
- I11-map (CI, gating), already proposed and fine: the z = 1 path reproduces the proton model exactly (bit-identical or ≤ 1e-12); the root residual is ≤ 1e-10; the identity Ω_ion/Ω_p = A*/F_scr holds at equal β (existing test).
- I11-rel must be gating, not informative. The decision's own argument is that the relative ion/proton response is the quantity theory fixes reliably (Δb cancels at equal velocity, the compound-rule and absolute f_dM offsets drop out). If that check is only informative, the z-dependence has no theory gate at all. Required row text: "I11-rel (LV, gating): |(θ_ion,model/θ_p,model) / (θ_MH,ion/θ_MH,p) − 1| ≤ 2 % at every I11 point (equal velocity and x/R); failure is STOP (amendment), never a tolerance change. The predictability of this ratio from the seen P1/P1b decomposition is disclosed."
- I11-abs (informative): with the reference's uncertainty stated (Pb residual, water status, f_dM's proton offset).
- The ±5 % model uncertainty statement goes in the metadata and docs, as drafted.

**4. (c) V2 as informative: honest, but its premise depends on unread files (important).**
- The draft rests "V2 informative by default" on Urban MSC for ions (Urban's charge enters only the prefactor; core plus tail; step dependence). P2 itself records: the ion msc model as "Urban, INFERRED (strong)"; the `G4EmBuilder` path (`isWVI` default) as OPEN. If GenericIon or alpha actually get WentzelVI plus single Coulomb scattering, the structural argument changes.
- Required: close P2's follow-ups (`G4EmBuilder.cc`, `G4hMultipleScattering.cc`, `G4VMscModel` defaults) and archive the P4 run-time EM-table printout (msc model, step-limit type, RangeFactor, lateral displacement for GenericIon, alpha and proton) before freeze. Pre-declare both branches: "If Urban: V2 informative (Case C), the raw ratio and expected direction reported." / "If WentzelVI/single scattering: V2 is decided by an amendment before any ion TOPAS run."
- Add before declaring the lateral gap. A thin-slab angular-width check against a Geant4 single-scattering configuration (`G4CoulombScattering` / a `g4em-standard_SS`-type module, if OpenTOPAS 4.3 offers it for GenericIon) would be an independent MC test of the charge dependence, built on theory rather than a condensed-history model. Add a P-check: availability and cost, pre-registered before any run.
- Intervention wording for the gap (only after P6 and that check fail): "Unavailable input: independent evidence of the charge dependence of lateral spreading for He/C/O in water at 50–430 MeV/u, i.e. measured lateral dose widths or angular spreads with stated uncertainties (candidates Parodi 2013, Schwaab 2011, Tessonnier 2017: no open numeric tables found), or an accessible single-scattering MC reference. Request: provide access to such data, or accept that V3-008 integrates with ion lateral spreading supported only by theory-derived mapping and implementation checks (no independent evidence), recorded as a gap in the ledger." It asks for data access or gap acceptance only, not a model choice. "No held-out analogue" stays.

**5. (d) Lineage wording: correct per P2, with one consistency fix.** These are labelled to match the P2 record's FACT, INFERRED and OPEN tags: LS switch at 2.0 MeV scaled; the q_eff²(ΔL_LS + 2·Barkas·(z−1)/z) along-step term; ion tables at ≤ 2.0 MeV scaled and 2 < Z ≤ 80; "ICRU 90 proton array at equal velocity" per the research note, OPEN; the alpha array, OPEN. "No lineage trigger met" is correct, because 2.0 MeV scaled is below the 2.5 MeV/u trigger. Fix: the §6 lineage row for Molière–Hanson says "Urban MSC differs structurally … (P2 FACT)". The structure of θ0 is FACT, but that Urban is the ion model is INFERRED. Write "(θ0 form P2 FACT; Urban as the ion msc model P2 INFERRED, settled by P2 follow-up/P4)".

**6. (e) What remains before freeze, in this order:** 1. P2 follow-ups and P4 (read the remaining Geant4 files; run the scorer-free smokes with the EM-table and printParameters printouts archived; then the new single-scattering availability check). 2. Exact-Ω specification (Ω or b; 2-D or bounded 1-D; formula, script and the I11-map/I11-rel tests committed). 3. P3 (the LS bound). 4. P5/P5b (G4EMLOW span and erratum, libdEdx). 5. P1 completion ((i) the Pb attribution on ICRU 49 ranges; (ii) a water or compound θ_Hanson check if a source exists; P1 recorded as FAILED for Pb if (i) does not clear it). 6. P6 (measured C/O range, and the extended search for lateral data). 7. P9 (the freeze record, binding the hashes of every script, record and printout above). Nothing produces a qualification result.

**Exact wording changes required:** 1. Decision §3: the labelling text of item 1, and P1 recorded as FAILED for Pb. 2. Decision §3 and plan: specify Ω vs b and the 2-D or bounded-1-D form (item 2). 3. Plan: I11-rel gating at 2 %, STOP on failure, with the predictability disclosure; I11-abs informative with the reference uncertainty; I11-map unchanged. 4. Plan V2, and §P: close the P2 follow-ups and P4 printouts before freeze, with both branches pre-declared; add the single-scattering availability check; the intervention text of item 4. 5. Decision §6: the Urban FACT/INFERRED split. 6. Plan §P: the order of item 6.

**Disposition (orchestrator, 2026-11-01).** All six wording changes applied in draft r5 (planner): the selection statement and P1 FAILED-for-Pb record; the exact-Ω mapping matches b = ln Ω + Δb at each projectile's own velocities and is 2-D (pv, p1v1), solved per step (tabulation only with a ≤ 0.1 % test, new check P11); I11-map gating, I11-rel gating at 2 % with STOP, I11-abs informative with the reference uncertainty; V2 with both pre-declared branches, the P2 follow-ups and P4 printouts closing before freeze, the single-scattering availability check P10 (conditional row V5-SS), the lateral-gap intervention wording; the Urban FACT/INFERRED split; §P order 1–7. The Codex REVIEW-98315d4e code findings on the committed reference module are fixed separately.

# Round V3-008-3 review (2026-11-01; Codex REVIEW-b469dff8 on head b3eef90b, decision-0043 r5 and v3-008-acceptance r5; verdict changes_required; four important findings, verbatim; bracketed locations as relayed by the orchestrator)

1. "The canonical chronology remains false. This draft was first committed at 2632a46 and repeatedly revised while checks remained open, yet lines 24–25 [validation/plans/v3-008-acceptance.md] say it is committed as the task's first scientific commit only after every pre-freeze check closes. At this head P1/P3/P4/P9 remain open. Describe the draft commits accurately and reserve this statement for the eventual frozen revision."

2. "The prior lateral-evidence misclassification was not propagated through all canonical records. P11 [validation/plans/v3-008/p11-exact-omega-mapping.md:5] still lists V1-held as evidence for lateral charge dependence even though V1-held contains no lateral statistic, directly contradicting the corrected acceptance-plan statement. The plan's P6 row also still calls the lateral search OPEN and does not reference the new p6-lateral-data-search.md evidence-gap result. Reconcile all records before freeze."

3. "V5-SS can qualify while the actual Coulomb-model physics and lineage remain OPEN. Its run condition [acceptance plan ~line 150] checks only the alias, process list and cost, while the P10 record says screening charge dependence, electron treatment, form factor and even registration are unresolved. Nevertheless the 3% tolerance assumes reference-model differences are ≤1% and the row is classified as independent MC. Record the newly cited constructor source in P10, inspect the selected scattering model, justify the tolerance, and determine whether this is independent physical evidence or only shared-model implementation evidence before freeze."

4. "The exact-Ω construction extrapolates the empirical proton f_dM fit far outside its documented domain. Decision 0043 describes f_dM as based on proton results at E ≤ 300 MeV, but P11 [p11-exact-omega-mapping.md:66] reports artificial proton-equivalent states reaching 2.5 GeV and merely extends the Bethe table to 8 GeV. Extending the range table does not validate f_dM at those states; current checks are self-consistency or comparisons with the same Molière construction. Bound this extrapolation, replace it with an in-domain formulation, or add independent coverage across the claimed ion-energy endpoints."

**Disposition (planner, 2026-11-01; DRAFT r6, uncommitted at the time of writing).**
1. Accepted. The plan's revision history now gives the exact draft commit history: r3 at 2632a46, the task's first scientific commit but a draft; r5 at 4927958; r5 amendments at 3d2547a, 070b07f, 00af9a6 and b3eef90; r4 never committed separately. The false sentence is withdrawn and quoted. The "first scientific commit … after every pre-freeze check is closed" wording is replaced by a statement reserved for the frozen revision 1, a separate reviewed commit after §P closes. Decision 0043's status line carries the same history.
2. Accepted.
   - P11 status paragraph: corrected to "V5-SS (if its run condition holds) or measured data only", with a dated correction note quoting the superseded text.
   - Plan P6 row and §P order step 6: the lateral search is recorded as EVIDENCE GAP, citing `p6-lateral-data-search.md`; the r5 status is kept as history.
   - STOP 4 and decision 0043 §3 (lateral gap): the same.
   - `p6-measured-data-access.md`: dated addendum.
   - The remaining "V1-held" mentions in the records concern range only, or already state that V1-held has no lateral evidence.
3. Accepted, partly blocked (STOP on one sub-item).
   - P10 §7 records the constructor ([SS] 8aeaf91c…) and the Coulomb sources, read from cache: [CS] db354d80…, [eCS] fa5ddeb9…, [WOKVI] a64c94b1…/3a3f1db0….
   - New FACT, [CS] L120–148: `G4CoulombScattering` gives the proton `G4eCoulombScatteringModel` + `G4WentzelOKandVIxSection`, but alpha and GenericIon `G4IonCoulombScatteringModel`. The two sides of the V5-SS ratio use different models.
   - The proton path is documented: Born screening form, electron scattering, form factor, spin term.
   - The ion path (`G4IonCoulombScatteringModel`, `G4ionCoulombCrossSection`) and the `G4EmParameters.cc` defaults are not in the cache. They are listed for acquisition, and V5-SS registration is blocked until they are read.
   - The run condition gains (iii), the recorded model lineage of both sides, with pre-declared branches A/B.
   - The 3 % is re-justified for branch A (non-cancelling terms ≲ 1 %), is provisional until then, and is invalid under branch B.
   - The class becomes "held-out single-scattering MC: independent of the multiple-scattering approximation, not of the screened-Coulomb cross-section physics". It is not independent physical evidence of the charge dependence of screening.
   - Recommendation to the orchestrator: file the lateral-gap intervention whatever the V5-SS outcome.
4. Accepted: bounded and declared, with conditional independent coverage; model unchanged.
   - P11b (SEEN pre-freeze): f_dM is a function of (pv, p1v1). Matched b is in-domain only where p1v1_eq ≤ 527 MeV: 11/36 I11 end states, and 100 % of θ² extrapolated at every standard point and upper end (≤ 0.78 decades in lg pv).
   - No in-domain b-only formulation exists: b_ion reaches 17.06, against 14.95 for qualified protons.
   - A fail-closed clamp would remove most of the domain, which is a scope decision and not taken.
   - The extrapolation is bounded against Molière–Hanson proton theory: drift +1.2 to +1.5 pp from 300 MeV to 3 GeV. This is theory-anchored, not independent. It is added as a ≤ 1.5 % uncertainty component (±5 % kept, quadrature ≈ 4.9 %).
   - New report-only row I11-ext. I0 persists `fdm_extrapolated` and the p1v1_eq energy. I11-abs and I11-rel report the extrapolated fraction.
   - Independent coverage: V5-SS end cases He 250, C 430 and O 430 (cases 28–33, ≈ +1 h REF), subject to the finding-3 STOP.

**Disposition continuation (planner, 2026-11-01; finding 3, ion side; after the orchestrator acquired the ion-path sources; uncommitted at the time of writing).**
- **3 (continued). The STOP on the ion sub-item is lifted. Branch B is decided.**
  - **Sources read (P10 §7.4.1):** `G4IonCoulombScatteringModel.cc/.hh` (ef9b31d4…/a1a9a4f5…), `G4IonCoulombCrossSection.cc/.hh` (87989b93…/cd5d0b5c…), `G4EmParameters.cc` (aaa143c1…) and `G4eCoulombScatteringModel.hh` (68ce763e…).
  - **FACT, ion path:**
    - nuclear scattering only, no target electrons;
    - the universal screening radius a0 (z^0.23 + Z^0.23), which depends on the projectile charge;
    - the Born factor 1.13 + (10π Z z α)^0.04 · 3.76 α′², also for H;
    - no form factor, no spin term, exact CM kinematics.
  - **Defaults:** exponential form factor, ScreeningFactor 1, FactorForAngleLimit 1.
  - **Why branch B:** the ion path differs from the proton path in every item branch A required. A scratch transform (SEEN pre-freeze, aggregates) shows that the r6 ion/proton design would deviate by +3.9 to +5.1 % for an exact Molière–Hanson ionmc, so the proton is not a valid z = 1 reference.
  - **Redesign (V5-SS row, P10 §7.5.1):**
    - ion/ion ratios at a common velocity on the ion path: He/C and He/O at 150 and 250 MeV/u, C/O at 430 MeV/u;
    - 8 cases (seed cases reassigned before any run, §C);
    - an electron energy cut ≥ T_max.
  - **Tolerance:** 3 % per pair, kept and re-justified before any run.
    - Budget: bias ≤ 0.9 % (computed −0.17 to −0.69 %); ionmc against Molière–Hanson ≤ 2 % (I11-rel gate); statistics 0.7–1.0 %.
    - Combined: quadrature 2.4 %, conservative 3.0–3.1 %.
  - **Class:** shared-lineage, held-out single-scattering MC; not independent physical evidence of the charge dependence of screening.
  - **Consequence for finding 4:** the V5-SS end-case coverage of the f_dM extrapolation becomes differential only. The independent-coverage claim is withdrawn in decision 0043 §3 and in the plan rows P11/V5-SS (dated notes).
  - **Correction:** the r6 "electrons alone ≈ 6 %" is corrected to −3.3 to −4.1 %.
  - **Lateral-gap intervention:** recommendation with arguments for and against in P10 §7.6.2 (file at the freeze regardless). The orchestrator decides, including whether the verbatim clause "or an accessible single-scattering MC reference" needs an amendment.

# Round V3-008-4 review (2026-11-01; Codex REVIEW-e0b586fd; four findings, verbatim as relayed by the orchestrator; bracketed locations as relayed)

1. (important, `validation/scripts/ion_em/p3_driver.py:91`) "P3 does not evaluate the documented nominal-Z LS correction against ionmc's actual effective-charge corrections. `r_eff` uses q in LS, Bloch and Mott, while `r_bare` uses z in all three; neither computes LS(z) − Bloch(q) − Mott1(q), even though P2 records Geant4 indexing LS by nominal `Zin` while ionmc applies Bloch/Mott with q. At 2 MeV/u, combining the committed bare LS residuals with the current q/Bloch/Mott formulas gives approximately −0.0465 for C and −0.0947 for O, rather than the reported +0.00195 and +0.00232 stopping-number residuals. Recompute this mixed convention and its range effect, or scientifically justify excluding it, before claiming the ≤0.1% bound and closing P3."

**Disposition of finding 1 (orchestrator, 2026-11-01; computed by the P3 worker after the review, SEEN): ACCEPTED, and the reviewer's estimate is confirmed** (ΔL at 2 MeV/u: mixed −0.04649 for C, −0.09469 for O). The driver now evaluates the third variant r_mixed = ΔL_LS(z) − Bloch(q) − Mott1(q) (`p3-lindhard-sorensen.md` § Mixed charge convention; `p3-output.json` keys `mixed_charge_convention`, `verdict.*.mixed_convention`). Result: r_mixed/S = −1.22 % (C) / −2.56 % (O) at 2 MeV/u, |r/S| > 0.1 % below 4.69 (C) / 6.58 (O) MeV/u, standard points 0.046–0.057 % (C) / 0.080–0.095 % (O); CSDA-range shift from 2 MeV/u exceeds the pre-declared 0.2 % trigger for end energies ≤ 6.86 (C) / 11.9 (O) MeV/u (≈ 0.4 / 0.9 µm absolute), 0.046 / 0.080 % at 430 MeV/u (0.14 / 0.19 mm); He unaffected (≤ 0.008 %). Under the ionmc convention (every correction at q) the earlier result stands (≤ 0.097 % of S, ≤ 0.095 % of range). **Consequence under the pre-declared P3 rule (plan r3: "if it exceeds 0.2 % of the CSDA range at any domain energy, an LS option is added to `BetheOptions` and its default is decided in an amendment before freeze"): the trigger is MET for the mixed convention, so the LS option is added (implementation task, same branch) and its default is decided by a reviewed amendment after a bounded physics assessment of the charge convention for the higher-order terms of a partially stripped ion (bare z in close collisions vs effective q; Geant4 indexes LS by nominal Z with the q² prefactor; assessment in progress). P3 is re-opened as "trigger met (mixed convention) — option pending"; no tolerance is changed; I3's 6 % (2–10 MeV/u) and 1.5 % (≥ 10 MeV/u) tolerances are unaffected by a ≤ 2.6 % / ≤ 0.1 % term respectively.


2. "The canonical chronology is still contradictory. Line 366 [decisions/0043-ion-em-transport.md] correctly says the mapping, tests, and first evaluation were produced together and the required ordering was not kept, but this operative "Freezing" section [line ~386] still states that they are committed before evaluation. Replace this with the actual disclosed chronology so the frozen decision cannot make a false evidence-independence claim."

3. "The withdrawn independent-coverage claim remains in the canonical P11 record. This section [p11-exact-omega-mapping.md:138] still says V5-SS adds independent coverage of the f_dM extrapolation and leaves that coverage conditional on P10, whereas the completed branch-B analysis explicitly withdraws absolute independent coverage and the acceptance plan now requires an intervention for the gap. The same stale state remains in decision 0043's lineage/limitations summaries. Propagate the final branch-B conclusion through P11 and the canonical decision before freeze."

4. "The branch-B V5-SS design and 3% acceptance tolerance depend on `ss_budget.py`, `ss_check.py`, `ss_cost.py`, and `ss_signal.py`, but these are uncommitted session-scratch scripts; hashes and aggregate claims alone do not let an independent reviewer inspect or reproduce the calculation. This is consequential because it supplies the −0.17% to −0.69% bias, discrimination claims, and tolerance budget, whose stated conservative total already reaches 3.0–3.1% against a 3.0% gate. Commit a reproducible calculation with tests and bind it in P9 before freezing the row."

**Disposition (planner, 2026-11-01; at head ce4d3d6a, uncommitted at the time of writing).**
1. See the orchestrator's disposition of finding 1 above (accepted; trigger met for the mixed convention; LS option pending its amendment).
2. Accepted. The operative "Freezing" block of decision 0043 §3 is replaced by the actual chronology: the mapping, its b formulas, the evaluation choice, the I11-map/I11-rel tests and the first I11-rel/I11-abs evaluation were produced together at 070b07f0 (regenerated at 00af9a6c; P11b at b3eef90b); all P11 values are SEEN. The P9 hashes record what was evaluated and do not establish pre-registration. The freeze revision may claim pre-registration only for the rows whose values are unseen, listed there (I0, I1, I3, I3b (i), I3c, I4, I4b, I5, I6, I7, I8, I8b, I9, I10, I12–I16, C-res, T1-ion, T12-ion, T13-ion, V1, V1-held, V2, V4 and the V5-SS statistic); I11-rel/-abs/-ext/-map, I2, I3b (ii) and the P1/P1b tables are carried as disclosed SEEN. The r5 sentence is quoted there as superseded. Parallel dated notes: plan §P paragraph, P9 row, P11 row.
3. Accepted. P11b (c) item 5 and its "OPEN after P11b" line are marked withdrawn; new section P11c states the final branch-B conclusion (differential coverage only, ≤ 1.5 pp; shared lineage; the gap routed to the intervention filed at the freeze per the orchestrator decision in plan STOP 4 (a)–(c); declared limitation ≤ 1.5 %, I11-ext, `fdm_extrapolated`). The P11 status paragraph gets a dated note. Decision 0043: item 5 of the §3 P11b note, the r6 "End cases" bullet and the P10-outcome paragraph are marked superseded; §6 lineage rows (f_dM; Geant4 single scattering, ion path now FACT, branch B) and the limitations MCS bullet get dated notes. Grep over the records for "independent coverage", "end cases", "cases 28", "proton cases": the remaining hits are history lines or already state the withdrawal.
4. Accepted. The four scratch scripts (hashes re-checked, unchanged) are replaced as the record by `validation/scripts/ion_em/ss_budget.py` (sha256 b893c0420bb492368d89943a8bd8820752ae903375921f9005c7dd2bc8f00156) and `tests/ionmc/test_ion_em_ss_budget.py` (551c75838dbedbf2be0a7b643beb1cd4b5ce48d2d15e15f0f4883bf00be1cfcc): pure numpy functions with docstrings citing the equations and the Geant4 source lines, a `main()` printing the aggregates, no computation or I/O at import; 14 CI-fast tests (known limits of the transform, vanishing bias terms, SCRIPT-GENERATED golden aggregates, the budget and the 3 % rule). Every recorded aggregate is reproduced unchanged (bias −0.17 to −0.69 %, electron share −3.3 to −4.1 %, radius −2.3 to −3.0 %, He-pair signal −5.6 to −7.2 %, discrimination, cost 3991 s at 2 µs/event, quadrature 2.41 %, conservative 3.02–3.14 %). Bound in plan P9 and P10 §7.8 (files, sha256, command); the V5-SS budget sentence cites the module. The conservative total reaching the 3 % gate is unchanged and stays disclosed: the 3 % is justified by the quadrature rule used since r5, and the conservative variant reaches it only at the inconclusive limit σ = 1 %. The 1 mm SS smoke cost (≈ 23 ms per C 290 primary, as relayed) gives ≈ 1.0 h REF for the 8 cases (P10 §7.8, plan P10 row and §E).

# Round V3-008-5 review (2026-11-01; Codex REVIEW-6a14f4ec at head ea220a88; two findings, verbatim as relayed by the orchestrator)

1. (decisions/0043-ion-em-transport.md:163) "The P3 disposition is still contradictory. This decision says the triggered LS option is added, while `p3-lindhard-sorensen.md:93` says whether the trigger requires an option remains undecided, and the acceptance row labels the option pending. The complete diff contains no `BetheOptions` implementation or default-selection amendment. Make the canonical records consistently describe a pending prescription, or implement the option and complete the reviewed default amendment before freeze."

2. (validation/scripts/ion_em/ss_budget.py:677) "The newly codified V5-SS gate ignores statistical uncertainty whenever sigma is at most 1%: for example, its test explicitly passes a 2.99% deviation with 0.9% sigma. This contradicts the plan's scalar-comparison TOST rule (`|x-r| + 1.645 sigma <= tol`) and can pass a result whose confidence bound is far outside 3%. It is especially material because the committed conservative budget is already 3.02–3.14%. Apply the declared confidence-bound rule or explicitly revise and scientifically justify the acceptance method and tolerance before freezing V5-SS."

**Dispositions (planner, 2026-11-01; at head ea220a88, uncommitted at the time of writing; the LS implementation is by a separate worker in the same commit).**

1. Accepted. The second option was taken: the option is implemented and the default amendment completed, as the orchestrator decided.
   - Implementation: `BetheOptions.ls ∈ {"off", "eff", "bare_mixed"}`, z ≥ 2 only, in `src/ionmc/physics/lindhard_sorensen.py` and `src/ionmc/physics/stopping.py`, with tests `tests/ionmc/test_ion_em_ls_option.py`.
   - Default: `"off"`, by the dated decision 0043 §1 "Amendment (2026-11-01) — LS option". It states:
     - the reasons. The trigger came from Bloch(z) − Bloch(q) (−1.29 / −2.64 % of S at 2 MeV/u for C/O), not from the LS-beyond part (+0.05 / +0.06 %). No theory fixes the charge. Geant4 is internally inconsistent: the LS path uses nominal Z with the q² weight, while the Bethe–Bloch path evaluates Bloch and Mott at q_eff. The plausible band lies inside the declared low-energy model uncertainty;
     - the two alternatives and why they were not chosen: `"eff"` is defensible but changes the result by ≤ 0.1 %; `"bare_mixed"` is a hybrid that would remove a TOPAS difference by construction;
     - what could decide it later: a screened-ion PASS/DPASS or partial-wave calculation, or IAEA low-Z stopping data; never the I3 rows;
     - an updated uncertainty statement.
   - The physics record is `validation/plans/v3-008/p3b-ls-default.md`. It carries the sources, the decomposition table, the Geant4 line references, the "not verified" list, and the v11.4.2 `G4LindhardSorensenData.cc` hash 70ae7c6b…, identical to v11.2.0, which closes the P3 §7 item.
   - Consistency edits, all dated notes:
     - `p3-lindhard-sorensen.md`, at the "decision for the lead orchestrator" bullet and at §7;
     - the plan P3 row: "CLOSED (re-opened and closed 2026-11-01): option implemented, default off by amendment", with hashes bound at commit;
     - the plan rows I3b(ii) (three variants, report-only), I3c (`bare_mixed` allowed for the exploratory TOPAS-attribution diagnostic) and P9 (iii) (the new src module, stopping.py, the test and the record).

2. Accepted. The declared confidence-bound rule is applied, and the design is re-evaluated before any run.
   - **Rule.** `ss_budget.gate` now passes iff |dev| + 1.645 σ ≤ 3 %. σ > 1 % remains a separate INCONCLUSIVE (= FAIL) precondition. The test requires 2.99 % / 0.9 % to FAIL and −0.69 % / 0.30 % to PASS, plus boundary cases.
   - **Statistics.** The module now computes the statistics from the reference density itself: window moments to x⁴, x/y joint moments by the 2-D transform, and the delta method.
     - Per case at 2e4 primaries, x and y pooled: 0.34 % (window kurtosis 2.76–2.78, x/y correlation 0.04–0.05).
     - σ of the ratio statistic at 2e4 / 1e5: 0.53 %.
   - **Power (P(all five pass), Monte Carlo with correlated cases).**
     - At 2e4 / 1e5: 0.994 for a correct model (true deviation = the computed bias), but 0.254 if ionmc's ratio is off Molière–Hanson by the SEEN I11-rel spread.
   - **Design decision, pre-registered before any V5-SS run:**
     - TOPAS 6e4 primaries per case: 2.94 h REF, scaled from the measured smoke, ≤ 4 h, ≤ 2 parts per case;
     - ionmc ≥ 4e5 histories per case;
     - σ = 0.30 %. P(all pass) ≥ 0.9999 for a correct model, 0.918 at the SEEN spread, 0.006 at the I11-rel bound (a finding, not a false fail);
     - fallback 4e4 before any case is dropped.
   - **Tolerance and gated ratios.** 3 % and the five gated ratios are kept. Under TOST the tolerance covers only the systematic items, 2.90 % linear at the item bounds. The room left for ionmc's model term is ≈ 1.6–1.8 pp, above the SEEN 1.13 pp. 3.5 % was not adopted.
   - **Where recorded:**
     - P10 §7.5.2 (tables) and a §7.8 dated note;
     - plan V5-SS row and §E (dated notes), P9 (iv) (new module and test hashes bound at commit);
     - decision 0043 §3 (dated note).
   - The module's new golden values are SCRIPT-GENERATED regression values. The round V3-008-4 `budget()` and all earlier goldens are unchanged.

# Round FREEZE-1 (2026-11-01; planner; freeze revision prepared at f694b747; not a review verdict — the independent Codex and scientific reviews of the freeze commit follow)

**What revision 1 binds** (`validation/plans/v3-008/p9-freeze-record.md`):
- 53 repository files by sha256 at f694b747: the Molière–Hanson reference, the P11 mapping script, the P3 and V5-SS calculations (6); their tests and the seed guard (7); the src LS option, exact-Ω mapping, stopping and config, plus `pyproject.toml` (5); the pre-freeze records and research notes (14); the ten scorer-free smoke cases and the reference-case README (21).
- 34 external sources by sha256: 29 acquired (Geant4 v11.4.2 sources of §G 1/1b/1c, G4EMLOW 8.8 and four members, libdEdx and three members, Sihver 1998) and 5 scratch-referenced lineage files.
- 15 required ancestor commits (2632a46 … f694b74, with 6e155b5 added by the planner); the clean-tree and diff rule; the seed reservation; the smoke run ids (nine as recorded prefixes, full ids to be copied by the orchestrator).
- The pre-registration status per row, updated to f694b747: the I0 mapping sub-items and the I1 mapping-callable parity exercised by the 683480b tests are carried as SEEN, besides the I11 rows, I2, I3b (ii) and the P1/P1b tables.
- The statement that no transport or qualification row result exists at f694b747.
- Orchestrator decisions recorded as amendments: the single-scattering clause of the lateral-gap wording amended; the intervention filed at the freeze naming both gaps (charge dependence of lateral spreading; f_dM extrapolation); the f_dM extrapolation a declared limitation; He `alpha` naming with the BetheBloch lineage; the LS default off and the V5-SS branch B / TOST / 6e4 design (earlier rounds) bound. V2 = Case C recorded as the consequence of its pre-declared branch.

**Audit outcome** (plan §P "Freeze audit (revision 1)"):
- CLOSED: P8, P7, P2b, P3, P4 (with one record gap carried in FC-2), P10.
- CLOSED-with-gap: P11 (SEEN; f_dM extrapolation declared), P5/P5b (TO-VERIFY items declared under row I3's pre-declared rule; icru90 origin OPEN, lineage only).
- DECLARED: P6 (V3 measured-range EVIDENCE GAP; lateral-data EVIDENCE GAP).
- **OPEN, freeze conditions:** FC-1 = P1 completion (order step 5: Pb attribution on ICRU 49 ranges; water/compound θ_Hanson check) not executed; FC-2 = P2 follow-up residuals (`G4BetheBlochModel.cc/.hh` pin and ICRU 90 branch, X0/RadL, λ charge scaling, `G4VMscModel.hh` initialisers, icru90 silent-fallback check) and the G4_WATER material table (I-value, RadL) required by the P4 row but not recorded.
- Row-evaluability finding, reported and not changed: I3 pre-declares only the 75 eV case of the libdEdx water I-value; if ICRU73_NEW were not the revised water table, no attribution is pre-declared.
- No tolerance and no row definition changed. Freeze commit only after FC-1 and FC-2 are disposed (executed record, or dated orchestrator amendment) and `p9-freeze-record.md` §7 passes.

**FC-2 disposition (2026-11-01; planner, at 6bdf3897; record `validation/plans/v3-008/p2b-follow-ups.md`, sha256 `9c70eb1b9d49f0375af7976232817dbda925576eeb4b64e31aee61b33dfb5457`; not a review verdict):**
- Disposed by an executed record (option (a)): bounded source reads and a read of the archived smoke logs; no new run, no transport or scored output.
- `G4BetheBlochModel.cc/.hh` (5dc2a900…/cf05e98c…) and `G4NistMaterialBuilder.cc` (a1f8944d…) are pinned and added to `p9-freeze-record.md` §3. The record is added to §2 (counts 54 repository / 36 external / 90 total).
- ICRU 90 branch settled for G4_WATER + opt4: alpha uses the ICRU 90 alpha array (≤ 1 GeV, the whole He domain); C/O use the ICRU 90 proton array × q², plus the LS terms. The L323 `G4EmCorrections` path (Bloch/Mott at the effective charge, confirming the P3b reading) is reached only outside the ICRU 90 branch, e.g. in the SS constructor (`Use ICRU90 data 0`). The `p3b-ls-default.md` §5 item is resolved, and a possible z = 1 double count of ΔL_LS is a lineage item. The alpha switch at 7.9452 MeV equals 2 MeV × M_α/M_p.
- G4_WATER: I = 78 eV is FACT (compiled constant, which overrides the formula value). X0 = 36.083 g/cm² is derived from the source composition; ratio to ionmc's cross-check X0 = 1.0000000. P4 "material table archived" is corrected (no material table in any log). The optional `Ma/Verbosity` smoke is not required.
- msc for GenericIon/alpha: Minimal / 0.2 / DispFlag 1 / Skin 3 is FACT (printout + `G4EmParameters` defaults + opt4). The `G4VMscModel.hh` initialisers are moot, except `dtrl` (OPEN). The λ charge scaling is OPEN at the inline level. Both are declared lineage residuals (V2 Case C attribution only).
- icru90 z6/z8 read evidence is ABSENT (silent fallback confirmed in the source); declared limitation of the lineage record. No frozen row depends on it.
- Orchestrator decisions recorded with it:
  - (d) I3: the result stands with the attribution "reference-table version unresolved (P5b)" if ICRU73_NEW is the pre-revision table, and is never re-graded. This resolves the row-evaluability finding above.
  - (e) the V3 measured-range gap is bundled into the intervention filed at the freeze (one request, two unavailable inputs).
- No tolerance, row definition or seed changed. FC-1 is not addressed by this note.

**FC-1 disposition (orchestrator, 2026-11-01).** `p1c-completion.md` (e3cd3258…): Pb cleared by attribution to the range-table input; compound check passes; water not validated directly; no bound file changed; P1 CLOSED. With FC-2 (above) both freeze conditions are disposed; this commit is the freeze commit of revision 1 subject to its independent review.

# Round R2-TOPAS (2026-11-07; planner; plan revision 2 prepared at d0829601; not a review verdict — the independent Codex and scientific reviews of the revision 2 commit follow)

**Trigger.** The operator's answer to IR-20261102-010435-B59B2D (resolved 2026-11-07 20:37 UTC), recorded verbatim in
`reference_requests/operator-directives/2026-11-07-topas-validation-target.md` (sha256 ca079e9f…): agreement with the
specified TOPAS configuration (OpenTOPAS 4.3 / Geant4 11.4.2, seven modules) is the acceptance reference for depth-dose
distributions and lateral profiles; EM-only TOPAS may be intermediate; experimental validation deferred; the f_dM
extrapolation is a model-design decision assessed through the ion predictions; cases, observables, precision and
criteria chosen autonomously; history preserved.

**What revision 2 does** (design record `validation/plans/v3-008/r2-topas-target.md`; plan §B revision 2 block, §C,
§D STOP 4 resolution and STOP 5, §E, §F, Amendments 9–15; decision 0043 revision 2):
- Rows: TD-EM (R80 min(0.35 %, 1.0 mm), distal width 10 % where W ≥ 1 mm, peak 0.5 mm, peak/entrance 3 %, plateau
  groups 2 %), TD-FULL (R80 of the total IDD with A_frag 0.15 % He / 0.10 % C, O; peak 0.5 mm; the rest report-only),
  TL-EM and TL-FULL (core σ band −5 %/+40 %; f_dM trend D_ext ≤ 12 %; FULL on the generation-0 profiles), held-out
  TD-H/TL-H (He 100, C 150, O 200; once, after grading and any amendment), TF-DIAG report-only.
- Cases: He 50/150/250, C 80/290/430, O 100/350/430 MeV/u (the I11/P11b grid: low ends in the f_dM domain, the rest
  extrapolated), water; IDD and lateral in separate TOPAS runs; IDD bins 0.02–0.1 mm, projected lateral bins
  0.01–0.02 mm (deviations from the 0.5 mm sketch, justified in r2 §3).
- Statistics: four TOPAS parts per case and configuration (1e5 EM, 2e4 FULL), t-inflated part SE; ionmc 1e6 per case.
- Seeds: formula unchanged; cases 1–12 re-assigned; part blocks 0–3/4–7/8–11/12–15; 205 bundles (192 evidence parts,
  13 execution-only smokes at 20600411–20600423), each passing `validate_case`.
- Cost: REF ≈ 8.7 h planning (≈ 14.3 h pessimistic) with V5-SS kept (2.94 h); sizing by timing smokes before any
  scored run, fallback order pre-registered.
- Superseded (never run): V1, V1-held, V2, V4. V3 stays an EVIDENCE GAP under revision 1 (history). V5-SS, I3 and all
  A-rows unchanged.
- The lineage is stated: depth-dose rows are shared-lineage for stopping; the lateral reference (Urban for ions) is
  independent but structurally deficient, so the lateral rows are a coarse band and V5-SS stays the 3 %-level check of
  b(z).

**Bookkeeping (Codex REVIEW-ad85dff1 on d0829601):**
- (a) P9 counts: 53 + `p2b-follow-ups.md` + `p1c-completion.md` + `p1c_attribution.py` + `test_ion_em_p1c.py` +
  `p1c-aggregates.json` + 2 intervention files = **60** repository paths (the orchestrator's 2026-11-02 note omitted
  `p1c-completion.md` and read 59); external 44; total 104. Corrected by a dated note in `p9-freeze-record.md` §7.
- (b) The plan's P1 row cited the FC-1 record as e3cd3258… (its first committed content at 655d16a). Rebound by a dated
  note to the current hash a1a24cd34c664a5dbb74b732bebff9244a55b4c12bfa30f0bff9c0cb8ebed0a5; e3cd3258… (655d16a) and
  f0dc6ce9… (fa59c64) are superseded history. The FC-1 disposition paragraph of round FREEZE-1 above also cites
  e3cd3258…; it is history of that round and refers to the same superseded content.

**Not established (reported, not resolved):** whether TOPAS scoring bins limit steps (step-length smoke); the Urban
deficit and its energy trend (scratch estimates; the band and the D_ext budget rest on them); acceptance of the
`g4radioactivedecay` alias by OpenTOPAS 4.3 (timing smoke; rejection = STOP 5); the FULL cost; the cost of ionmc's fine
projected lateral tallies.

## Round R2-TOPAS-2 (2026-11-07; Codex REVIEW-7a7186a2 on b8c9836f)

**Finding 1 (p9-freeze-record.md:345), verbatim:** "The revision-2 freeze binding is internally inconsistent and cannot pass R2.6. The R2.1 whitelist omits three paths actually changed by this commit: `validation/scripts/ion_em/gen_r2_bundles.py`, `tests/ionmc/test_r2_bundles.py`, and `validation/reference_cases/README.md`. The README remains bound to its revision-1 hash `2c160fb5…`, but its current hash is `e16cc9d3…`. R2.5 also excludes the newly hash-bound generator, test, and intervention response, reporting 472 repository paths/516 total instead of 475/519 based on the record's own hash entries. Reconcile the whitelist, rebound hashes, counts, and R2.6 wording exactly."

**Disposition 1: accepted.** `p9-freeze-record.md` R2.1 whitelist now lists the generator, the test and the README (and the 410 bundle files, 26 of which changed in this round); the README is rebound in R2.2 to its current hash (`2c160fb5…` kept as superseded revision-1 history in §2.5); R2.5 states the counting rule (one count per distinct repository path with a current full sha256 entry; rebound paths count once; superseded hashes and the digest do not count) and the numbers: 60 + 415 = **475** repository paths, 44 external, **519** total; the generator, the test, `response.json` and the directive each have exactly one current hash entry (§2 dated note of 2026-11-07); R2.6 test 3 now checks the 475 paths, the digest and the whitelist exactly.

**Finding 2 (gen_r2_bundles.py:229), verbatim:** "The pre-implementation timing/step smokes claim that every scorer output is deleted unread, but their cases list full IDD/lateral files in `expected_outputs`. The controlled reference runner inventories and permanently archives the entire work directory and has no deletion mechanism, leaving these low-statistics reference observables readable/materializable before ionmc samples are locked. This defeats the stated evidence-blinding control. Use non-disclosing smoke outputs, enforce an auditable access quarantine/deletion mechanism, or lock ionmc predictions before these smokes run, and update the plan and bundle metadata accordingly."

**Disposition 2: accepted; first option (non-disclosing smoke outputs).** The 12 timing smokes and the step-length smoke carry no IDD or lateral scorer: one 1-bin `SmokeEdep` EnergyDeposit placeholder writing `smoke_unread.csv` (P4 convention), `case.json` scorers "none: timing/step smoke; the placeholder output is not an observable and is not read (P4 convention)", expected_outputs `["smoke_unread.csv"]`; particle, energy, physics, cuts, MaxStepNumber, seeds, primaries and the step smoke's verbosity unchanged. The 192 evidence bundles are byte-identical. The new bundle-set digest is `4a091ff6dafb250ac0c77431264d8f0c09594a8159d3a6f6d0d0b1e30c0fd373`; 26 file hashes changed (appendix A of `r2-topas-target.md`). The blinding control no longer relies on deletion: no reference observable exists before the ionmc samples are locked. `r2-topas-target.md` §6 step 2, §7, §4.3, §9 and the plan's order/blindness paragraph (dated note) are updated; the guard test now asserts that no smoke bundle has an IDD/lateral scorer and every evidence bundle has them. Consequence disclosed: the step-length smoke no longer carries the fine evidence bins, so it records primary step lengths in the unbinned geometry only and cannot decide whether those bins limit the step; that question stays OPEN (it was already listed as not established), and the bundles define the reference either way.
