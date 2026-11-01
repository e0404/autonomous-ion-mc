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
