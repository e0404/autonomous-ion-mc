# V3-008 acceptance criteria (revision 2 — TOPAS-agreement target, operator directive 2026-11-07; revision 1 frozen at the reviewed freeze commit of `p9-freeze-record.md`, kept as history)

> **Status: revision 2 (planner, 2026-11-07, at d0829601; round R2-TOPAS in `validation/plans/v3-008/review.md`). Revision 2 becomes the frozen revision by a reviewed commit; revision 1's freeze is superseded before any result.** Trigger: the operator's answer to IR-20261102-010435-B59B2D (resolved 2026-11-07 20:37 UTC), recorded verbatim in `reference_requests/operator-directives/2026-11-07-topas-validation-target.md` (sha256 ca079e9f0af990cdff9d8efde60d38ac55e2d0fc0fb9ef22355a81f4a2ab0b26): agreement with TOPAS (OpenTOPAS 4.3 / Geant4 11.4.2) is the acceptance reference for depth-dose distributions and lateral profiles; the final configuration is the seven-module list `g4em-standard_opt4`, `g4h-phy_QGSP_BIC_HP`, `g4decay`, `g4h-elastic_HP`, `g4stopping`, `g4ion-QMD`, `g4radioactivedecay` (the "specified TOPAS configuration", FULL); EM-only TOPAS may serve as intermediate validation; experimental validation is deferred; complementary checks are retained; the f_dM extrapolation is a model-design decision assessed through the ion predictions against the reference; existing results and historical verdicts are preserved. **What revision 2 changes** (design record `validation/plans/v3-008/r2-topas-target.md`; section Amendments, "Revision 2"): new gating rows TD-EM, TD-FULL, TL-EM, TL-FULL and the held-out TD-H/TL-H (§B, revision 2 block); V1, V1-held, V2 and V4 superseded (never run; no seed consumed); V3 stays an EVIDENCE GAP under revision 1 and is superseded as a target; STOP 4 resolved by the directive, new STOP 5 (§D); TOPAS case numbers re-assigned (§C); 205 TOPAS bundles added (bound at the revision-2 freeze, `p9-freeze-record.md` "Revision 2 additions"); REF planning total ≈ 8.7 h (§E). **Unchanged:** every A-row (I0–I16, I11-map/rel/abs/ext, T/R/C rows), I3/I3b/I3c, V5-SS (gating, complementary), the seed formulas, the statistics rules. At d0829601 no V3-008 ionmc ion-transport run, LV/HR qualification step or scored ion TOPAS run exists; revision 2 is pre-registered before any result. Also fixed in revision 2 (Codex REVIEW-ad85dff1 bookkeeping): the P1 row's FC-1 record hash (dated note in the P1 row) and the P9 repository-path count (`p9-freeze-record.md` §7, dated note).


> **Status: revision 1 (FROZEN at the reviewed freeze commit; prepared 2026-11-01 by the planner at f694b747, round FREEZE-1 in `validation/plans/v3-008/review.md`).** This revision freezes the V3-008 acceptance criteria together with decision 0043 ("frozen with the plan revision 1"). Freeze record: `validation/plans/v3-008/p9-freeze-record.md` (D-E8 pattern: bound sha256, required ancestor commits, clean-tree and diff rule, pre-registration status per row). **Freeze conditions.** The §P freeze audit (below the §P table) finds two items that this plan requires before freeze and that have no closing record at f694b747: **FC-1** (P1 completion, order step 5) and **FC-2** (P2 follow-up residuals and the P4 material printout). Revision 1 becomes frozen only through a reviewed commit that disposes both (by an executed, committed record or by a dated orchestrator amendment) and passes `p9-freeze-record.md` §7. Until then it is a freeze candidate, and no row may be run. No V3-008 transport or qualification result exists at f694b747. The DRAFT r6 status text below is kept as history.

> **Status: DRAFT r6 (2026-11-01, amended for Codex REVIEW-b469dff8; r5 was amended for REVIEW-12614ced), not yet frozen.** The acceptance criteria are frozen before any TRANSPORT or QUALIFICATION result; the pre-freeze construction/calibration computations P1/P1b (Molière–Hanson reference, candidate tables) and P11 (I11-rel/I11-abs values, evaluated at 070b07f0 and regenerated at 00af9a6c) exist, are disclosed and are not qualification evidence. This draft becomes revision 1 (frozen) by a reviewed commit after the pre-freeze checks of §P (steps 1–7, incl. P10 and P11) are executed; P4/P10 smoke runs carry no scorers. Independent review: rounds V3-008-1, V3-008-2 and V3-008-3 (Codex REVIEW-b469dff8) (`validation/plans/v3-008/review.md`). Seed bases per the shared ledger (V3-005B commit 06974d07): V3-008 reservation 20572004–20632003.


**Status.** DRAFT r6 (planner, 2026-11-01). Revision history:
- r2 was revised against `research-models.md` and `research-data-sources.md`.
- r3 applied the nine wording changes of review round V3-008-1.
- r4 folds in the committed pre-freeze records P1/P1b, P2, P7 and P8 (at 6e155b56). It changes:
  - §P status;
  - the I11 row (now informative) and the new rows I11-map and I11-rel;
  - V2 (informative by default);
  - the §B lineage;
  - STOP 2/3;
  - §G.
- r5 applies the six wording changes of review round V3-008-2:
  - I11-rel becomes gating; I11 becomes I11-abs (informative);
  - the exact-Ω specification (P11);
  - the V2 branches;
  - the single-scattering check (P10, row V5-SS);
  - the lateral-gap intervention wording;
  - the §P order (P2 follow-ups and P4 close before freeze).
- r6 (2026-11-01) applies the dispositions of Codex REVIEW-b469dff8 (round V3-008-3, `review.md`):
  - the commit history of this draft stated exactly (below); the "first scientific commit" statement reserved for the frozen revision;
  - the P6 lateral-data search recorded as EVIDENCE GAP in the P6 row, the §P order and STOP 4; V1-held/V3 carry no lateral evidence (P11 note corrected);
  - V5-SS: the single-scattering model lineage read (P10 §7): the proton and the ions use different Geant4 single-scattering models; registration BLOCKED until the ion-side model is read; run condition, branches A/B, evidence class and the status of the 3 % tolerance revised; end cases 28–33 added;
  - the f_dM domain extrapolation of the exact-Ω mapping quantified (P11b) and declared: new report-only row I11-ext; I0, I11-abs and I11-rel report the extrapolation.
- r6 continuation (2026-11-01, round V3-008-3 finding 3, after the ion-path sources were acquired; dated notes, r6 text kept):
  - the ion single-scattering path is read (P10 §7.4.1);
  - V5-SS is decided as **branch B**: ion/ion ratios at a common velocity on the ion path; the proton cases are dropped;
  - the 3 % is re-justified and kept (P10 §7.5.1); the f_dM-extrapolation coverage by V5-SS is reduced to differential;
  - lateral-gap recommendation in P10 §7.6.2.
- Round V3-008-4 (2026-11-01, Codex REVIEW-e0b586fd; dated notes, earlier text kept): the P11 chronology stated in decision 0043 §3 "Freezing" (pre-registration only for unseen rows; §P paragraph, P9 and P11 rows); the withdrawn independent coverage of the f_dM extrapolation propagated (§P, P11 row; P11c; decision 0043 §3, §6, limitations); the V5-SS budget committed as `validation/scripts/ion_em/ss_budget.py` with tests and bound in P9 (V5-SS row, P10 §7.8); the 1 mm SS smoke cost (P10 row, §E); finding 1: `p3-lindhard-sorensen.md` § Mixed charge convention.
- Round V3-008-5 (2026-11-01, Codex REVIEW-6a14f4ec; dated notes, earlier text kept): finding 1 — P3 closed with the LS option implemented and default "off" by the decision 0043 §1 amendment (`p3b-ls-default.md`; P3, I3b, I3c, P9 rows); finding 2 — the V5-SS gate is the plan's TOST rule (§ Statistics) and the statistics are raised to 6e4 TOPAS primaries / 4e5 ionmc histories per case before any run (V5-SS row, §E, P10 §7.5.2).
- **Revision 1 (FREEZE REVISION; planner, 2026-11-01, at f694b747; round FREEZE-1 in `review.md`; dated notes, earlier text kept):**
  - §P freeze audit: the status and closing record of every P-item; FC-1 and FC-2 named as freeze conditions; the handling of the TO-VERIFY items (P5b, I3) stated;
  - P9 freeze record `validation/plans/v3-008/p9-freeze-record.md`: 53 repository files and 34 external sources by sha256, 15 required ancestor commits, the clean-tree and diff rule, the pre-registration status updated to f694b747;
  - orchestrator decisions recorded as amendments (section Amendments): the single-scattering clause of the lateral-gap wording amended (STOP 4; decision 0043 §3); the intervention filed at the freeze naming both gaps (section "Intervention at the freeze (to be filed by the orchestrator)"); the f_dM extrapolation kept as a declared limitation; the He REF/V1 cases use `alpha` (= GenericIon(2,4)) with the BetheBloch lineage (P4). Already recorded in round V3-008-5 and bound here: the LS option with default off (decision 0043 §1 amendment); V5-SS branch B with the TOST rule and the 6e4-primary design;
  - consequence of a pre-declared branch recorded: V2 = Case C (report only; P2 F5 and P4);
  - no tolerance and no row definition changed.
- Nothing in §P has produced a transport or TOPAS result.
- **Revision 2 (planner, 2026-11-07, at d0829601; round R2-TOPAS in `review.md`; dated notes, earlier text kept):** the operator's directive of 2026-11-07 (IR-20261102-010435-B59B2D resolved 2026-11-07 20:37 UTC) makes agreement with the specified TOPAS configuration the acceptance reference; the TOPAS-agreement rows TD-EM, TD-FULL, TL-EM, TL-FULL, TD-H/TL-H and the report-only TF-DIAG are added (§B revision 2 block; design and justification in `validation/plans/v3-008/r2-topas-target.md`); V1, V1-held, V2, V4 superseded; V3 and the lateral gap recorded as revision-1 history and superseded as targets; STOP 4 resolved, STOP 5 added; the TOPAS case numbers re-assigned within the unchanged formula; §E re-estimated; the two REVIEW-ad85dff1 bookkeeping fixes. No A-row, I3 row or V5-SS definition, tolerance or seed changed. No result exists.
- **Commit history of this draft (exact; corrected in r6 for Codex REVIEW-b469dff8 finding 1).**
  - First committed with decision 0043 at 2632a46 (2026-11-01 12:44, DRAFT r3). That commit is the
    first scientific commit of task V3-008, but it is a draft, made while every pre-freeze check
    was open.
  - Revised by commits while pre-freeze checks remained open: 4927958 (DRAFT r5, 13:23), 3d2547a
    (13:29), 070b07f (14:06), 00af9a6 (16:52) and b3eef90 (17:10), all r5 amendments; r4 was never
    committed separately; r6 is this revision.
  - At b3eef90 P1, P3, P4 and P9 are open. The r3–r5 sentence "committed with decision 0043 as the
    first scientific commit of task V3-008, after independent review and after every pre-freeze
    check of section P is closed" did not describe these draft commits and is withdrawn.
  - **Reserved for the frozen revision:** revision 1 (frozen) is committed by a separate reviewed
    commit, after independent review and after every pre-freeze check of section P is closed; that
    commit carries the freeze record P9.
  - **Dated note 2026-11-01 (revision 1):** after b3eef90 the draft was further revised, with pre-freeze
    checks still open, at 41a1f8a (r6), a66867d, ce4d3d6, ea220a8 (round V3-008-4), daa1520 (round V3-008-5),
    683480b, 7fac30a and f694b74. Revision 1 is prepared at f694b74; the freeze commit is the reviewed commit
    named by `p9-freeze-record.md`, made once FC-1 and FC-2 are disposed.
- It is committed before any V3-008 TRANSPORT or QUALIFICATION result exists (no ionmc transport run, no LV/HR qualification step, no scored reference-engine run).
- Pre-freeze construction/calibration computations exist and are disclosed: P1/P1b (Molière–Hanson reference, candidate tables) and P11 (I11-rel/I11-abs values at 070b07f0, regenerated at 00af9a6c). They are SEEN and are not qualification evidence.
- The pre-freeze checks are construction work: theory computations, source reading and data
  inspection. They use no ionmc transport output and no TOPAS output.
- Later changes need a dated amendment with the scientific reason. Failures and consumed bases are
  preserved.

**Base and scope.**
- Base: `v3/develop` `c790df36`. The plan is independent of PR #82.
- Scope: 4He, 12C, 16O primaries, electromagnetic only (decision 0043).
- Domains: He 50–250, C 80–430, O 100–430 MeV/u.
- Evidence is claimed for water. The `materials.py` set gets theory rows only.

**Tiers.**

| Tier | Meaning |
|---|---|
| CI | GitHub, synthetic fixtures, warp CPU, ≤ 3 min added |
| LV | Local exact-SHA, single process, cached data. Data-dependent tests skip in CI and fail under `IONMC_REQUIRE_DATA=1` |
| HR | Host runner, CUDA |
| REF | v3 reference MCP, TOPAS |

**Statistics.** `statistics.NormalDist` plus Wilson–Hilferty, no scipy. Profile and scalar
comparisons reuse the frozen T12 machinery of `v3-003-acceptance.md` (amendments ⁹–²²): grouped
sparse bins, permutation/bootstrap-t calibration, and inconclusive counts as fail. TOST form: pass
iff |x − r| + 1.645 σ_comb ≤ tol.

**Notation.**
- R = CSDA range of the project table at the source energy.
- R80 and R20 are distal depths on 0.1 mm bins.
- **Standard points**: He 150, C 290, O 350 MeV/u. **Ends**: He 50/250, C 80/430, O 100/430 MeV/u.
- Numbers labelled "estimate" are planning values from hand arithmetic, not acceptance values.

## P. Pre-freeze checks (must be closed and recorded in decision 0043 before freezing)

P11: mapping, tests and the first I11-rel/I11-abs evaluation produced together at 070b07f0 (chronology in `p11-exact-omega-mapping.md`; SEEN, not blind). Genuinely held-out qualification evidence for the charge dependence of lateral spreading is V5-SS (if its run condition holds) or measured data only, never I11-rel. V1-held and V3 carry no lateral statistic and provide no lateral evidence. **r6:** measured data would be the only independent physical evidence, and P6 found none (EVIDENCE GAP, `p6-lateral-data-search.md`); V5-SS is held-out single-scattering MC evidence that shares the screened-Coulomb cross-section physics (evidence class in the V5-SS row; P10 §7). **Dated note 2026-11-01 (round V3-008-4, Codex REVIEW-e0b586fd findings 2–3):** "genuinely held-out qualification evidence … is V5-SS … or measured data" is qualified: under branch B, V5-SS is held-out but shared-lineage evidence, not independent physical evidence of the charge dependence of lateral spreading, and it covers the f_dM extrapolation only differentially; both gaps go into the intervention filed at the freeze (STOP 4, orchestrator decision; `p11-exact-omega-mapping.md` P11c). The P11 chronology is stated in decision 0043 §3 "Freezing" (corrected text): mapping, tests and the first evaluation were produced together, so the freeze claims pre-registration only for rows whose values are unseen. P5: G4EMLOW inspected (`p5-g4emlow-ion-stopping.md`): Geant4 files = lineage only; P5b (`p5b-libdedx-icru73.md`): libdEdx ICRU73_NEW identified as the I3 evaluation source (role exploratory until the row is frozen); whether it equals the 2009 Sigmund water revision and the I-value of its water table are TO-VERIFY.

Closed: P8, P7, P2 (source reading). P2b is very probably inapplicable. **Order before freeze (r5):**
1. P2 follow-ups, P4 and P10 (scorer-free smokes only).
2. P11: the exact-Ω specification, with the I11-map and I11-rel tests committed.
3. P3.
4. P5/P5b.
5. P1 completion.
6. P6 (measured C/O range; lateral-data search). r6: both recorded as EVIDENCE GAP (`p6-measured-data-access.md`, `p6-lateral-data-search.md`).
7. P9.

| # | Check | Outcome feeds |
|---|---|---|
| P8 | Repository facts. **CLOSED** (`validation/plans/v3-008/p8-repository-facts.md`): `A16_BASELINE` = f3a1dd62 (full SHA, fail-closed resolution); ASTAR object dd4a2e50… spans 1e-3 to 1e3 MeV total, i.e. ≤ 250 MeV/u; species ids 5 and 6 are free on this branch (PR #82 to be checked at integration); proton gates at `config.py:653–657` and `transport/tables.py:200–205` | R2, I2, I0 |
| P7 | Seed bases. **CLOSED**: ledger on the V3-005B branch (commit 06974d07, `V3_008_RESERVED = range(20572004, 20632004)`); guard test `tests/ionmc/test_v3_008_seed_families.py` committed | C |
| P2 | Geant4 v11.4.2 ion models from source. **Source reading CLOSED; follow-ups OPEN and required before freeze (order step 1)** (`p2-geant4-ion-models.md`; [URB] 22c5f52b…, [LS] e16ed6cd…, [OPT4] 04337300…). FACT: GenericIon = `G4ionIonisation` + `G4LindhardSorensenIonModel`; switch at 2.0 MeV scaled; above it the `G4BetheBlochModel` table + q_eff²(ΔL_LS + 2·Barkas·(z−1)/z), with no separate Bloch/Mott; ICRU 73/90 ion data only ≤ 2.0 MeV scaled and 2 < Z ≤ 80, never alpha; ion msc = `G4hMultipleScattering("ionmsc")`, model not set (Urban inferred); Urban θ0 = 13.6·\|charge\|·√y·invbetacp·(c1 + c2 ln y), y = t/X0, with no charge/β in the log (frozen with line references). The r3 lineage trigger ("ion tables above 2.5 MeV/u") is NOT met. **Follow-ups, required before freeze:** read `G4EmBuilder.cc` (40782f5c…) and `G4hMultipleScattering.cc` (a6643869…), both materialized in the v3-008 cache, and the `G4VMscModel` defaults (Urban for ions; WentzelVI for protons; alpha assignment); pin `G4BetheBlochModel.cc/.hh` (the ICRU 90 proton/alpha array claims stay "research note, partly confirmed" until then); read `G4IonICRU73Data.cc` (cf771a5a…) and `G4IonisParamMat.cc` (990c6f31…; Zeff of G4_WATER: θ0 differs by 2.5 % between 10/3 and 6.6) and X0 of G4_WATER; λ scaling and lateral displacement for ions (`G4VMscModel`, `G4EmTableUtil`) **Revision 1 (2026-11-01): follow-ups partly closed (`p2-geant4-ion-models.md` F0–F5 at 070b07f; P4 printouts); the residuals (`G4BetheBlochModel.cc/.hh` pin and its ICRU 90 branch, X0 of G4_WATER, λ charge scaling, `G4VMscModel.hh` initialisers) have no record → OPEN, freeze condition FC-2 (freeze audit below)** **Dated note 2026-11-01 (FC-2 disposition; `p2b-follow-ups.md`; Amendments, "Amendment (2026-11-01) — FC-2 disposition"):** `G4BetheBlochModel.cc/.hh` pinned (5dc2a900…/cf05e98c…); ICRU 90 branch settled for G4_WATER + opt4 (alpha ≤ 1 GeV: ICRU 90 alpha array; GenericIon: ICRU 90 proton array × q², then the LS terms; the L323 corrections are not reached; one INFERRED premise: the `G4ICRU90StoppingData` index/span); X0 of G4_WATER derived from the source composition (ratio to ionmc's cross-check X0 = 1.0000000); final msc parameters for GenericIon/alpha FACT (Minimal, 0.2, DispFlag 1, Skin 3); the `G4VMscModel.hh` initialisers moot except `dtrl` (OPEN, lineage); λ charge scaling OPEN at the inline level and declared a lineage residual of the TOPAS reference (V2 Case C attribution only). Follow-ups CLOSED for the freeze | §B lineage, V2 branch |
| P2b | **As defined in r3, very probably INAPPLICABLE**: opt4 protons use WentzelVI (INFERRED strong), so the archived T15 does not test an Urban method, and the r3 Fermi–Eyges δ_pred is not a well-defined Urban prediction (truncated core + tail; per-length width step dependent, factor 0.756 between 0.1 and 1 mm; TOPAS ion steps set by the bundle). Consequence (r5): V2 follows the pre-declared branches. **P2b′** (only if a later amendment reintroduces a graded Urban route): the Urban-emulator method is checked against a TOPAS configuration in which protons use Urban (source-verified constructor, e.g. opt3 if `ConstructCharged(..., false)`, TO VERIFY), same geometry and declared step model, run before any ion TOPAS run; e_δ ≤ 2 % is required for the extension | V2 |
| P4 | TOPAS settings, checked on smoke runs (execution only, not graded): `alpha` vs `GenericIon(2,4)` particle naming; `GenericIon(8,16)`; BeamEnergy = total kinetic energy; acceptance of a very large electron production cut so that deltas deposit locally (e.g. `d:Ph/Default/CutForElectron = 1 km`; exact parameter TO VERIFY); run-time printout confirming the LS ion model and the ICRU90 flag (`/process/em/printParameters` or verbose); throughput per 1e5 primaries P4 smoke cases carry no dose, IDD or phase-space scorers (or their outputs are deleted unread by the runner); only the run log, EM-model printout and wall time are recorded. The measured throughput sizes every REF case into parts of ≤ 1200 s planned (section E)  **r5: required before freeze (order step 1).** The scorer-free smokes archive the EM-table / `printParameters` printouts: msc model, step-limit type, RangeFactor and lateral displacement for GenericIon, alpha and proton; ionisation model names; the G4_WATER material table (RadL) **P4 progress 2026-11-01 (at 41a1f8ac; execution only):** proton (REF-188fa1d8…), alpha (REF-ed13a238…) and GenericIon(2,4) (REF-c764466c…) smokes ran without error (1e4 histories each, ≈ 1.3 s wall, EM-only opt4). Printout facts (Geant4 11.4.2 / OpenTOPAS 4.3): the EM parameters are printed at the first BeamOn by default (`Use ICRU90 data 1` with opt4 although `Ph/Default/ICRU90` was left at its default; tables to 10 GeV, so EMRangeMax is honoured); proton: msc WentzelVIUni (StepLim Minimal, Rfact 0.2, lateral displacement on) + CoulombScat (eCoulombScattering), hIoni Bragg < 2 MeV / BetheBloch; **alpha: msc UrbanMsc, ionIoni BraggIon < 7.9452 MeV / BetheBloch (NOT LindhardSorensen)**, StepFunction (0.1, 0.02 mm); **GenericIon: msc UrbanMsc, ionIoni LindhardSorensen 0–10 GeV**, StepFunction (0.1, 0.001 mm); nuclearStopping ICRU49 < 1 MeV for all three. Consequence: `alpha` and `GenericIon(2,4)` are different ionisation physics in TOPAS (same Urban msc); the He REF/V1 cases must fix the naming in the freeze revision (proposal: GenericIon(2,4) so He/C/O share the LS + Urban lineage; the alpha alternative is the BetheBloch lineage). `b:Ph/ListProcesses` raised no error and had no visible effect (the process table is printed by default). Not yet verified: the production-threshold table is not printed at default verbosity (1 km electron cut unconfirmed) and which process set each He case actually used (identical printouts) → verbose follow-up cases `he-water-150mevu-em-smoke-{alpha,genericion}-verbose` (5 histories, TrackingVerbosity 1, RunVerbosity 2, seeds 20600408/20600409); carbon, oxygen and the two SS smokes pending **P4 completed 2026-11-01 (all nine smokes at 41a1f8ac/a66867de; execution only):** carbon GenericIon(6,12) REF-5910656d… and oxygen GenericIon(8,16) REF-8031592a… ran without error or exception (1e4 histories, 2.1 / 2.2 s wall; same opt4 printout: GenericIon ionIoni LindhardSorensen, msc UrbanMsc). Verbose cases (5 histories, TrackingVerbosity 1, RunVerbosity 2; REF-c98e198f… alpha, REF-7537d7e4… GenericIon(2,4)): the couple table confirms the 1 km electron range cut is accepted and gives an e⁻ production threshold of 9.55 GeV in G4_WATER (no delta electrons; gamma 990 eV, e⁺ 56 keV, proton 5 keV at the 50 µm cut); **`GenericIon(2,4)` is mapped by Geant4/TOPAS to the predefined `alpha` particle** (the track printout reads `Particle = alpha` for both namings), so He always uses the alpha process set (BraggIon/BetheBloch + UrbanMsc), never LindhardSorensen — the naming question is moot and the He lineage is BetheBloch; the 600 MeV alpha primary takes ≈ 59 ionIoni-limited steps over its track, the first 16 mm long (StepFunction (0.1, 0.02 mm) with msc StepLim = Minimal), i.e. the Urban lateral displacement is applied over 16 mm steps at the entrance — recorded for the V2 branch and the P10 design. Throughput at 1e4 primaries, 300 mm water, EM-only opt4: p 1.2 s, He 1.3 s, C 2.1 s, O 2.2 s (≈ 0.1–0.2 ms per primary; no scorer). `b:Ph/ListProcesses` is accepted silently. P4 is CLOSED; the He REF/V1 cases use `alpha` (= GenericIon(2,4)) and the freeze revision records the BetheBloch lineage for He **Dated note 2026-11-01 (FC-2 disposition; `p2b-follow-ups.md`; Amendments, "Amendment (2026-11-01) — FC-2 disposition"):** **correction:** the smoke logs do not archive the G4_WATER material table (no density, I-value or RadL in any of the nine logs; the verbose cases print only the production-cut couple table). "Archived" is corrected to: source definition (`G4NistMaterialBuilder.cc` a1f8944d… L712–715: 1.0 g/cm³, H₂O by atom count, I = 78 eV as the compiled constant, which overrides the formula value, L230–235) + derived X0 (36.083 g/cm², Tsai form). Optional later smoke, not required: OpenTOPAS `i:Ma/Verbosity = 1` (reports each material definition; whether it prints RadL is unverified); no documented material-table printout parameter was found | V1, V2 branch, V3, host budget |
| P10 | **Single-scattering availability. CLOSED for coverage (source):** FACT, `G4EmStandardPhysicsSS.cc` v11.4.2 (cache 8aeaf91c…): L213–219 `particle = G4GenericIon::GenericIon(); G4ionIonisation* ionIoni = new G4ionIonisation(); … ionIoni->SetEmModel(new G4LindhardSorensenIonModel()); ph->RegisterProcess(ionIoni, particle); ph->RegisterProcess(new G4CoulombScattering(false), particle);` (no msc process for GenericIon); L222 `G4EmBuilder::ConstructChargedSS(hmsc);` → alpha `G4ionIonisation` + `G4CoulombScattering(false)` ([EMB] L158–160) and protons `G4hIonisation` + `G4CoulombScattering(false)` ([EMB] L212, L221–222); L107 `param->SetMscThetaLimit(0.0);`. GenericIon, alpha and protons all get pure single Coulomb scattering in the SS constructor. Design consequence (P10 record §3): the r5 design (10/50 mm, ≥ 1e5 primaries) is infeasible (≈ 10²–10³ h); V5-SS is redefined on the P10 design (1 mm slab, 2e4 primaries per case). **Remaining, as the V5-SS run condition (not a freeze blocker):** a scorer-free smoke confirming the OpenTOPAS 4.3 alias `g4em-standard_SS` and the process list (`CoulombScat`, no msc, for GenericIon, alpha, proton) and measuring the cost; lineage of the Coulomb model (screening z-dependence, electron treatment, form factor: `G4CoulombScattering.cc`, `G4eCoulombScatteringModel.cc`, `G4WentzelOKandVIxSection.cc`, OPEN) is declared in the V5-SS row. **r6 (Codex REVIEW-b469dff8; P10 §7): the Coulomb sources are read** ([CS] `G4CoulombScattering.cc` db354d80…, [eCS] `G4eCoulombScatteringModel.cc` fa5ddeb9…, [WOKVI] `G4WentzelOKandVIxSection.cc/.hh` a64c94b1…/3a3f1db0…). FACT [CS] L120–148: the SS constructor sets no model ([SS] L219; [EMB] L146–165, L221–229), and `G4CoulombScattering::InitialiseProcess` then selects `G4eCoulombScatteringModel` for the proton (mass < 1 GeV) but `G4IonCoulombScatteringModel` for alpha and GenericIon (mass > 1 GeV or type "nucleus"). **The ions and the z = 1 reference use different single-scattering models.** The ion-side model (`G4IonCoulombScatteringModel`, `G4ionCoulombCrossSection`) is not in the cache. Its screening, electron treatment and form factor decide whether the ratio statistic cancels model differences. **Reading them is required before freeze (STOP on this sub-item, P10 §7.4)**; the V5-SS branch (A or B) is then decided by the plan revision that freezes. **Dated note 2026-11-01 (continuation; P10 §7.4.1–7.7.1):** the ion path is read (`G4IonCoulombScatteringModel.cc/.hh` ef9b31d4…/a1a9a4f5…, `G4IonCoulombCrossSection.cc/.hh` 87989b93…/cd5d0b5c… (capital I; the "G4ionCoulombCrossSection" above is a misspelling), `G4EmParameters.cc` aaa143c1…, `G4eCoulombScatteringModel.hh` 68ce763e…). FACT: • nuclear scattering only (no target electrons); • the universal screening radius a0(z^0.23 + Z^0.23); • the Born factor 1.13 + corr·3.76 α′², with corr = (10π Z z α)^0.04, also for H; • no form factor, no spin term, exact CM kinematics. Defaults: NuclearFormfactorType exponential, ScreeningFactor 1, FactorForAngleLimit 1. The STOP on this sub-item is lifted. **Branch B** is decided (V5-SS row) **SS smokes 2026-11-01 (a66867de; REF-61e64b1b… carbon 1e3, REF-285075c5… proton 1e3; execution only):** the OpenTOPAS 4.3 alias `g4em-standard_SS` is accepted; printout: proton CoulombScat = `eCoulombScattering` (0 < θ < 180°, table 0–10 GeV), GenericIon and alpha CoulombScat = `IonCoulombScattering` (BuildTable = 0, on the fly), **no msc process for any of them**; ionisation hIoni Bragg/BetheBloch (p), ionIoni LindhardSorensen (GenericIon), with StepFunction (0.2, 0.1 mm) and 56-bin tables in the SS constructor. Cost: the 150 MeV proton full track (300 mm box) 0.61 s per primary (614 s / 1e3, no stuck tracks); **every one of the 1e3 carbon ions at 290 MeV/u was killed at TOPAS's `Ts/MaxStepNumber` = 1 000 000 after ≈ 1.25 s with ≈ 2530 of 3480 MeV remaining** (stderr: 1000 stuck tracks, 2.53e6 MeV lost), i.e. ≈ 1.25 µs per single-scattering step and > 1e6 steps per full carbon track. Consequence for the V5-SS design (1 mm slab): the step count per mm is far below the cap (order 1e4–1e5 per mm for C at the entrance; estimate), so the cap is not reached in a 1 mm slab, but every V5-SS case must set `Ts/MaxStepNumber` explicitly (≥ 1e7) and must fail closed on any killed track (stderr check, `Ts/KilledTrackMaxCount = 0` if the parameter exists — TO VERIFY); the 2e4-primary cost per case is bounded by ≈ 2e4 × (steps per mm) × 1.25 µs plus the 0.6 s-class proton-like overhead, consistent with the §E budget (re-check at the freeze with the measured step count of a 1 mm smoke, which is the next smoke to add). Run condition (i) alias and (ii) process list are satisfied; (iii) lineage recorded (branch B) **Dated note 2026-11-01 (round V3-008-4; 1 mm SS smoke, execution only; facts as relayed by the orchestrator, run not re-read by the planner):** case `carbon-water-290mevu-ss-slab1mm-smoke` (GenericIon(6,12), 3480 MeV, 100 histories, 1 mm G4_WATER, `g4em-standard_SS`, `Ts/MaxStepNumber` 1e8), REF-35df82fcd7723f6a9a0c-db361e98 at ce4d3d6a: no killed track; `IonCoulombScattering` for alpha and GenericIon; ≈ 23 ms wall per carbon primary. With 1.299e4 ion-path events per C 290 primary (`ss_budget.py`) this is ≈ 1.77 µs per event; the 8-case design costs ≈ 1.0 h REF (§E). The cost part of the run condition is met (≤ 4 h) | V5-SS, lateral gap |
| P11 | **Exact-Ω specification (order step 2):** commit the b-matching mapping (b = ln Ω + Δb_Fano at the respective velocities; 2-D dependence on (pv, p1v1); p1v1_eq/pv_eq = p1v1/pv), the evaluation choice (per-step root solve, or a (pv, p1v1) tabulation with a test bounding the θ interpolation error ≤ 0.1 %), and the I11-map and I11-rel tests, **before the mapping is evaluated against the Molière–Hanson reference**; hashes into P9. **Recorded** (`p11-exact-omega-mapping.md`): the mapping, the tests and the first I11-rel/I11-abs evaluation were produced together at 070b07f0, so the required order (mapping and tests committed before any evaluation) was NOT kept; the values are SEEN construction/calibration results, not qualification evidence, and a later rerun cannot make these points independent. I11-rel (equal-velocity proton, corrected after Codex REVIEW-42ab68dd) is evaluated pre-freeze: not blind; the freeze carries it as a disclosed gate. **r6, P11b (Codex REVIEW-b469dff8 finding 4; SEEN pre-freeze):** • f_dM is an explicit function of (pv, p1v1), not of b. • The equivalent proton state lies in the proton-qualified f_dM domain (p1v1_eq ≤ pv(300 MeV) = 527 MeV) only for He 50, C 80 and part of O 100 (11 of 36 I11 end states). • At every standard point and upper end, 100 % of θ² comes from extrapolated states: up to 0.37 decades in lg pv for He and 0.78 for C/O (≈ 2.5 GeV). • Decision (decision 0043 §3, note of 2026-11-01): the model is unchanged, because no in-domain formulation exists for this domain (b_ion up to 17.1 against 14.95 for qualified protons). • The extrapolation is declared and bounded against the Molière–Hanson proton theory (row I11-ext). It is carried in the ±5 % uncertainty, persisted in the metadata (I0) and covered by the V5-SS end cases, subject to the V5-SS STOP. **Dated note 2026-11-01 (continuation):** under V5-SS branch B the coverage is differential only (common-velocity ion/ion ratio, residual ≤ 1.5 pp, below the tolerance). The independent coverage of the extrapolation is withdrawn. It stays declared, theory-bounded (I11-ext) and carried in the ±5 % (P10 §7.5.1). **Dated note 2026-11-01 (round V3-008-4, Codex REVIEW-e0b586fd findings 2–3):** the requirement "… committed **before the mapping is evaluated against the Molière–Hanson reference**" at the head of this row is the r5 requirement; it was not met and is not claimed (decision 0043 §3 "Freezing", corrected text). The f_dM-extrapolation gap is routed to the lateral-gap intervention filed at the freeze (orchestrator decision (a), (c); P11c) | §3, I11-map, I11-rel, I11-ext |
| P3 | Committed numpy LS script: ΔL_LS − (Bloch + first-order Mott) for z = 6, 8, water, 2–430 MeV/u, tested against the Bloch limit at small z. Expected ≤ 0.1 % of S (estimate). If it exceeds 0.2 % of the CSDA range at any domain energy, an LS option is added (amendment). **CLOSED 2026-11-01** (`p3-lindhard-sorensen.md`, SEEN pre-freeze; threshold pre-declared in r3): max r/S 12C 0.060 %, 16O 0.097 %, 4He 0.010 % (protons: LS − Bloch ≤ 0.093 %); max ΔR/R (range from 2 MeV/u) 12C 0.060 %, 16O 0.095 %; the 0.2 % range trigger is not met and the ≤ 0.1 % expectation holds (16O with a small margin); no LS option. Sign positive everywhere (ionmc S slightly low, range slightly long, ≤ 0.19 mm for 16O at 430 MeV/u). Files `validation/scripts/ion_em/lindhard_sorensen.py` 3bfeb288…, `validation/scripts/ion_em/p3_driver.py` 50db856a…, `tests/ionmc/test_ion_em_lindhard_sorensen.py` d35ca4e3…, `validation/plans/v3-008/p3-output.json` 491355f2… (bound at P9). Exploratory, SEEN, not settled by source: the Geant4 LS table may enter `G4LindhardSorensenIonModel` with half the ionmc-convention magnitude (≈ 0.36 % of S for 16O at 430 MeV/u if so) — relevant to the I3c/V1 attribution only **Re-opened 2026-11-01 (Codex REVIEW-e0b586fd finding 1; SEEN):** the mixed charge convention r_mixed = ΔL_LS(z) − Bloch(q) − Mott1(q) (Geant4's indexing: LS at nominal Z, the ionmc corrections at q) gives r/S = −1.22 % (C) / −2.56 % (O) at 2 MeV/u, |r/S| > 0.1 % below 4.69 / 6.58 MeV/u, and a CSDA-range shift from 2 MeV/u above the 0.2 % trigger for end energies ≤ 6.86 (C) / 11.9 (O) MeV/u (≈ 0.4 / 0.9 µm; 0.046 / 0.080 % at 430 MeV/u); He ≤ 0.008 %. **The pre-declared trigger is MET for the mixed convention → the LS option is added to `BetheOptions` and its default is decided in an amendment before freeze** (physics assessment of the charge convention for the higher-order terms in progress; the ionmc all-at-q convention keeps the earlier ≤ 0.1 % result). Status: trigger met — option pending; no tolerance changed. **Dated note 2026-11-01 (round V3-008-5, Codex REVIEW-6a14f4ec finding 1). Status: CLOSED (re-opened and closed 2026-11-01): option implemented, default off by amendment.** `BetheOptions.ls ∈ {"off", "eff", "bare_mixed"}`, z ≥ 2 only (z = 1 unchanged), default `"off"` (all corrections at q), by decision 0043 §1 "Amendment (2026-11-01) — LS option"; physics assessment `validation/plans/v3-008/p3b-ls-default.md` (the trigger came from Bloch(z) − Bloch(q), −1.29 / −2.64 % of S at 2 MeV/u for C/O, the LS-beyond part being +0.05 / +0.06 %; the convention is not fixed by theory and is a declared model uncertainty at 2–10 MeV/u; no I3/V1 value used). Files: `src/ionmc/physics/lindhard_sorensen.py` sha256 800d1f0a13524af5bad3946032f894479f606ba9b1482bf3ea190eeb9b4d3a6a, `src/ionmc/physics/stopping.py` sha256 951f29acaae262255e1f0b6bc831d11118fe55025bf0a8cfb93c4cefad8db828, `tests/ionmc/test_ion_em_ls_option.py` sha256 c67966c8d2c3d9564032e33e666321dd11bd8d2386780afa0cb74850455ca84d, `validation/plans/v3-008/p3b-ls-default.md` sha256 01f75993b259a3b1c2d86a3ee88b345c048c1750f3963d005ed7f2dac7614a0b. P3 §7's "data files not shown identical" item is closed (v11.4.2 `G4LindhardSorensenData.cc` sha256 70ae7c6b… = the v11.2.0 file; GitHub-tag fidelity assumed). The sentence above "physics assessment … in progress" is superseded; no tolerance changed | decision §1, I3b |
| P5 | Data inspection of G4EMLOW 8.8 (`https://geant4-data.web.cern.ch/datasets/G4EMLOW.8.8.tar.gz`, 350 393 595 B, MD5 328330009df633f7e9b3a9f445745298), members `ion_stopping_data/icru73/z{6,8}_G4_WATER.dat` and `icru90/z6_G4_WATER.dat`. Checks: (i) energy span covers ≥ 2–430 MeV/u; (ii) units and per-nucleon axis; (iii) erratum provenance. The original ICRU 73 water used I = 67.2 eV and the erratum uses 78 eV; decide by the high-energy logarithmic slope against the Bethe L0 at two I values (aggregate verdict only, no values in Git) and by documentation; (iv) origin of `icru90/z6`. **P5b:** locate libdEdx (GPL) and, if it ships ICRU 73 (rev.) water files, use it as a cross-check or substitute. **P5 outcome (`p5-g4emlow-ion-stopping.md`):** archive acquired (sha256 b60cfd63…, exploratory); spans cover 2–430 MeV/u; the icru73 water files are the Sigmund-revised tables (README), so the erratum is resolved by documentation; origin of `icru90/z6` and `z8` OPEN (README silent; ICRU 90 has no oxygen table). The README states "DATA ARE *NOT* FOR COMMERCIAL USE AND MUST BE USED WITHIN GEANT4": the Geant4 files cannot be the independent I3 evaluation source and are a LINEAGE record only. libdEdx acquired (sha256 b941a0c3…, 962,747 B). **P5b (`p5b-libdedx-icru73.md`):** `data/raw/ICRU73_NEW.dat` (749768c1…), grid `icru73_newEng.dat` (8b5b809f…, 53 energies 0.025–1000 MeV/u), embedded `dedx_icru73new.h` (f6d088b2…; ions Z = 1–18, targets 276 liquid water and 104 air) identified as the I3 evaluation source; equality with the 2009 Sigmund revision and the I-value of its water table TO-VERIFY (libdEdx metadata lists I = 75.0 eV for target 276, a libdEdx material constant). **Oxygen/carbon lineage checks (P2 follow-up):** option 4 reads `ion_stopping_data/icru90/z8_G4_WATER.dat` (and `z6`) below the 2.0 MeV scaled switch and falls back silently when a file is absent; record existence, sha256 (z6 d9d4809d…, z8 863898af…), span/origin (0.00025–1000 MeV/u; origin OPEN) and a silent-fallback check (printout or log evidence that the file was read) for both **Dated note 2026-11-01 (FC-2 disposition; `p2b-follow-ups.md`; Amendments, "Amendment (2026-11-01) — FC-2 disposition"):** evidence ABSENT: the carbon, oxygen and carbon-SS logs contain no line showing the icru90 files were read; the loader is silent on a missing file (`RetrieveVector(ost, false)`, `G4IonICRU73Data.cc` L284, L404–416) and the LS fallback is Bragg × q² (`G4LindhardSorensenIonModel.cc` L281–291); the only read trace needs the loader's header-level `fVerbose > 0` (L417–419; whether an EM verbose command reaches it is OPEN). Declared limitation of the lineage record; no frozen row depends on it (I3 uses libdEdx; the files act on ≲ 0.03 mm of C/O range) | I3 |
| P1 | Molière–Hanson reference. **PARTLY CLOSED** (`p1-moliere-hanson-reference.md`, P1 and P1b). Done: the reference with the Fano correction is implemented; proton validation passes Be (0.10 %), Al (0.02 %) and Cu (0.21 %); **the ≤ 1.0 % proton validation FAILED for Pb** at x/R1 ≥ 0.9 (+1.10 %, +1.32 %; `xfail(strict=True)`, tolerance not loosened; recorded as FAILED and preserved). Findings carried into decision 0043 §3: water is not validated directly; the compound Fano rule is derived, not sourced; proton f_dM sits +2.5 to +3.5 % above the water reference; c_Ω ≈ 0.09–0.11 (r3's 0.13 was wrong); the r3 selection rule selects nothing (seen tables disclosed in §3). r4/r5 replace the selection rule by the exact-Ω mapping (the exact version of r3's first-ranked candidate; the selection rule was revised after the r3 rule reached STOP 2 on seen P1/P1b numbers — chronology and disclosure in decision 0043 §3) and demote the absolute reference comparison to the informative row I11-abs, with I11-map and I11-rel gating. **P1 completion (order step 5):** Pb attribution on ICRU 49 ranges; a water/compound θ_Hanson check if a source exists (Gottschalk polystyrene/Lexan/water values); if Pb is not cleared, P1 stays FAILED for Pb, as recorded. The reference-uncertainty statement goes into decision 0043. The P1 driver evaluates the exact I11 grid. Codex REVIEW-98315d4e: n_Ω is the full derivative of the proton Ω (or b) including the screening and Fano β-dependence; this matters for the superseded linearised (d) only. The Kanematsu ion form is moot **Revision 1 (2026-11-01): P1 completion has no record at f694b747 → OPEN, freeze condition FC-1 (freeze audit below; `p9-freeze-record.md` §0). The reference-uncertainty statement is in decision 0043 §3 ("order 3 %"). P1 stays FAILED for Pb** **FC-1 disposed 2026-11-01 (`p1c-completion.md` e3cd3258…):** P1 completion done — Pb CLEARED by attribution to the range-table input (ICRU 49 node path: ≤ 0.044 %; the committed xfail test unchanged), compound check within −0.57 to −0.04 % of Gottschalk's θ_Hanson for polystyrene/Lexan/Lucite/Nylon/Teflon, water not validated directly (b low by 0.040 ± 0.002 → ≤ 0.16 % in θ at I11-like B); I11-abs keeps its ≈ 3 % water reference uncertainty. P1 CLOSED **Dated note 2026-11-07 (revision 2; Codex REVIEW-ad85dff1 bookkeeping finding (b)):** the hash "e3cd3258…" above is the FC-1 record as first committed at 655d16a; it is superseded history. The record was revised in the freeze commits (fa59c64: f0dc6ce9…; d082960) and is bound at its current, final content: `validation/plans/v3-008/p1c-completion.md` sha256 a1a24cd34c664a5dbb74b732bebff9244a55b4c12bfa30f0bff9c0cb8ebed0a5 (= `p9-freeze-record.md` §2.4 dated note). The outcome stated above is unchanged | §3, I11, I11-rel, I11-map |
| P6 | Measured-range candidate: Sihver, Schardt, Kanai, Jpn. J. Med. Phys. 18 (1998) 1 (J-STAGE). Check open access and licence, whether Bragg-peak/R80 positions for C and O in water are tabulated with uncertainties, and whether beam-line water-equivalent offsets are given. Bounded investigation: two access approaches, ≤ 2 h. **r4: also search for measured ion lateral widths in water** (Parodi et al. 2013 J. Radiat. Res. (HIT); Tessonnier 2017; Schwaab 2011), numeric availability TO VERIFY; a negative result is recorded as the MCS evidence gap of STOP 4. **P6 outcome (`p6-measured-data-access.md`): the Sihver 1998 candidate FAILED** (measured Bragg curves normalised horizontally to the calculated peaks, figures only, oxygen beam is 18O) → V3 is an EVIDENCE GAP; the lateral-data search was OPEN at r5. **r6: the lateral-data search is recorded (`p6-lateral-data-search.md`, committed at b3eef90): EVIDENCE GAP.** No measured, independent, tabulated He/C/O lateral width or angular spread in water with stated uncertainties was found. The Tessonnier 2017 thesis has figures plus text-quoted, beamline-confounded ratios without uncertainties; Schwaab 2011 and Tessonnier 2017 PMB are paywalled; the "Phillips" ion data are of unknown provenance; Gottschalk's survey (arXiv:1204.4470 §5.3) concurs. The OPEN items left in that note (paywalled table content, the Bevalac Ne/Fe/U report) do not change the verdict | V3 (gap), MCS gap (the lateral-gap intervention is used only after P6 and V5-SS/P10 fail) |
| P9 | Freeze record (D-E8 pattern). It binds the committed sha256 listed in the P1 note (`validation/scripts/ion_em/__init__.py` 7962c4af…, `validation/scripts/ion_em/moliere_hanson.py` 5a99b3b8…, `tests/ionmc/test_ion_em_moliere_hanson.py` 789fd4b1…; the hashes are regenerated from the final committed files at every revision of this draft and bound verbatim at the freeze), the P11 mapping script and the I11-map/I11-rel tests, the P2 follow-up source hashes (40782f5c…, a6643869…, 990c6f31…, cf771a5a…), the P2 source hashes ([URB], [LS], [OPT4], [ICRU73]) and the P3 script and output. The freeze commit must have those commits as ancestors, the recorded hashes must match, and the tree must be clean. **Dated note 2026-11-01 (round V3-008-4, Codex REVIEW-e0b586fd findings 2 and 4):** (i) P9 also binds the V5-SS budget calculation `validation/scripts/ion_em/ss_budget.py` (b893c0420bb492368d89943a8bd8820752ae903375921f9005c7dd2bc8f00156) and `tests/ionmc/test_ion_em_ss_budget.py` (551c75838dbedbf2be0a7b643beb1cd4b5ce48d2d15e15f0f4883bf00be1cfcc) (P10 §7.8; command `uv run --extra dev python -I validation/scripts/ion_em/ss_budget.py`); the V5-SS row is frozen only with these hashes. (ii) The P9 hashes of the P11 mapping, script and tests record what was evaluated at 070b07f0/00af9a6c; they do not make I11-map/I11-rel/I11-abs/I11-ext pre-registered. The freeze revision claims pre-registration only for the rows whose values are unseen, listed in decision 0043 §3 "Freezing" (corrected text). **Dated note 2026-11-01 (round V3-008-5, Codex REVIEW-6a14f4ec findings 1 and 2):** (iii) P9 also binds the LS option: `src/ionmc/physics/lindhard_sorensen.py` sha256 800d1f0a13524af5bad3946032f894479f606ba9b1482bf3ea190eeb9b4d3a6a, `src/ionmc/physics/stopping.py` sha256 951f29acaae262255e1f0b6bc831d11118fe55025bf0a8cfb93c4cefad8db828, `tests/ionmc/test_ion_em_ls_option.py` sha256 c67966c8d2c3d9564032e33e666321dd11bd8d2386780afa0cb74850455ca84d and the record `validation/plans/v3-008/p3b-ls-default.md` sha256 01f75993b259a3b1c2d86a3ee88b345c048c1750f3963d005ed7f2dac7614a0b; (iv) the hashes of `ss_budget.py` and its test in (i) are superseded by the round V3-008-5 versions (TOST rule, window statistics, power): `validation/scripts/ion_em/ss_budget.py` sha256 13158795c32094609188b22ec64b2d612e3037a1d782b392801f5aa53af4e009, `tests/ionmc/test_ion_em_ss_budget.py` sha256 ca49735b5b00831aeaabd59d63faecfc23b98fe0d4247024ca963d7e9a77f097 (P10 §7.5.2) | freeze |

### Freeze audit (revision 1; planner, 2026-11-01, at f694b747)

Status words: **CLOSED**; **CLOSED-with-gap** (closed; a declared gap or limitation remains); **DECLARED** (an evidence gap or limitation carried into the ledger); **OPEN-FC** (required before freeze by this plan's text, no closing record at f694b747: freeze condition).

| # | Status at revision 1 | Closing record | Remaining or declared |
|---|---|---|---|
| P8 | CLOSED | `p8-repository-facts.md` (6e155b5) | species ids vs PR #82 are checked at integration (decision 0043 §4) |
| P7 | CLOSED | ledger commit 06974d07; guard test (sha256 in P9 §2.2); logged run in round FREEZE-1 | — |
| P2 | source reading CLOSED; follow-ups partly closed → **OPEN-FC (FC-2)**; **2026-11-01: FC-2 disposed** (`p2b-follow-ups.md`) | `p2-geant4-ion-models.md` F0–F5 (070b07f): GenericIon and alpha msc = Urban (FACT), no `G4CoulombScattering` for ions in opt4, Urban Zeff of water = 10/3, the loader reads `icru90/z{6,8}_G4_WATER.dat` for G4_WATER; P4 printouts: UrbanMsc for alpha and GenericIon, WentzelVIUni for protons | `G4BetheBlochModel.cc/.hh` not acquired or pinned through the reference tool and its ICRU 90 branch not settled (`p3b-ls-default.md` §2, §5): the "ICRU 90 proton/alpha array" claims stay "research note, partly confirmed", as the P2 row pre-declares; X0 (RadL) of G4_WATER; λ charge scaling; `G4VMscModel.hh` initialisers. By F5, X0, λ, lateral displacement and `dtrl` affect only the withdrawn Urban emulator, not V2 Case C. No gating row's definition or evaluation depends on FC-2; the lineage statements (§B) already assume the shared lineage, which is the conservative reading |
| P2b | CLOSED (inapplicable as defined; now FACT) | P4 proton printout: msc WentzelVIUni + CoulombScat | P2b′ only by a reviewed amendment |
| P4 | CLOSED, with one record gap carried in FC-2; **2026-11-01: record gap closed by correction** (`p2b-follow-ups.md` §2) | plan P4 row; nine smokes and the 1 mm SS smoke (case hashes and ids in P9 §2.5) | the r5 sentence requires archiving the G4_WATER material table (RadL); its content (I-value, RadL) is not recorded → FC-2. **He naming decided:** `alpha` (= GenericIon(2,4)), BetheBloch lineage (§B, revision 1 note) |
| P10 | CLOSED (coverage, lineage, branch B, cost) | `p10-single-scattering-availability.md` §1–7.8; SS smokes; 1 mm SS smoke REF-35df82fcd7723f6a9a0c-db361e98 | V5-SS run conditions, checked at run time and not freeze conditions (row text): production-cut table printout (electron energy cut ≥ T_max); `Ts/KilledTrackMaxCount` (TO VERIFY) or the stderr check; `Ts/MaxStepNumber` ≥ 1e7; measured cost of the 6e4 design ≤ 4 h |
| P11 | CLOSED-with-gap (SEEN, disclosed; the required order was not kept) | `p11-exact-omega-mapping.md` (070b07f, 00af9a6, b3eef90; P11c; P11d at 683480b) | f_dM extrapolation: DECLARED limitation (≤ 1.5 %, I11-ext, `fdm_extrapolated`), named in the intervention at the freeze |
| P3 | CLOSED (re-opened and closed) | `p3-lindhard-sorensen.md`, `p3-output.json`, `p3b-ls-default.md`; decision 0043 §1 amendment | convention band declared (decision 0043 §1 uncertainty statement) |
| P5 / P5b | CLOSED-with-gap | `p5-g4emlow-ion-stopping.md` (G4EMLOW = lineage only); `p5b-libdedx-icru73.md` (libdEdx ICRU73_NEW = I3 source) | DECLARED uncertainties, not freeze conditions (rule below): ICRU73_NEW = the 2009 water revision, TO-VERIFY; its water I-value, TO-VERIFY; block presence, by a loader test at implementation; origin of `icru90/z{6,8}` OPEN (lineage only), so the I3 clause "vs ICRU 90 for C if P5(iv) establishes its origin" is not triggered and I3 has no ICRU 90 carbon comparison; the silent-fallback run-time check → FC-2; **2026-11-01:** evidence absent, declared limitation of the lineage record (`p2b-follow-ups.md` §4) |
| P1 | PARTLY CLOSED → **OPEN-FC (FC-1)**; the P1 proton validation FAILED for Pb (preserved) | `p1-moliere-hanson-reference.md` (6e155b5, 4927958); reference-uncertainty statement in decision 0043 §3 ("order 3 %") | P1 completion (order step 5) not executed: Pb attribution on ICRU 49 ranges; the water/compound θ_Hanson check if a source exists. Required by the §P order and §F gate 1, hence a freeze condition. If it is disposed by amendment it is recorded as "not done" and Pb stays FAILED; only the I11-abs (report-only) wording depends on it |
| P6 | DECLARED (two EVIDENCE GAPs) | `p6-measured-data-access.md` (V3 gap); `p6-lateral-data-search.md` (lateral gap) | V3 row = EVIDENCE GAP; the lateral gap is filed at the freeze (section "Intervention at the freeze") |
| P9 | this revision | `validation/plans/v3-008/p9-freeze-record.md` | completed by the freeze commit (record §0 and §7) |

**Rule applied to TO-VERIFY items (revision 1).** A TO-VERIFY item blocks the freeze only if this plan names it as a pre-freeze check, or if a frozen row cannot be evaluated or graded without it. The two P5b items (revision identity; water I-value) do neither. The P5 row records P5b as done with these items open, and row I3 pre-declares their handling: the model is evaluated at 78 eV and never tuned, the bounded I-value contribution is stated, and a failure within it stays FAIL with the attribution, never re-graded. They are declared uncertainties of the I3 result. The same holds for I3's "reference uncertainty ≈ 1–2 % (TO VERIFY)": it is part of the tolerance defence; the tolerance is frozen as written, and the statement is reported with the result as unverified. **Disclosed limitation of row I3 (reported to the orchestrator, not resolved here):** the row pre-declares the handling of a 75 eV table only. If ICRU73_NEW were not the revised water table, I3 would still evaluate as written, but no attribution for that case is pre-declared. **Dated note 2026-11-01:** resolved by orchestrator decision (d) (row I3 note; Amendments): the result stands as evaluated with the attribution "reference-table version unresolved (P5b)" and is never re-graded.

## A. Verification rows

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| I0 | **Capability contract.** He/C/O sources run on python, warp-cpu and warp-cuda. Fail closed for: He-3 and any other ion; energy outside the domain or E0 + 6σ_E above it; `nuclear=True` with z ≥ 2; `nist-star` for C/O; ASTAR with a non-water material; LET medium ≠ water; unproducible species/generation tallies. Requested/effective configuration persists projectile, z, mass, `e_cut_ion_mev_per_u` and the MCS model id. **r6 (P11b):** it also persists the f_dM domain status of the source state: the kinetic energy of the proton-equivalent p1v1_eq at the source, and a flag `fdm_extrapolated` when that energy exceeds the proton-qualified 300 MeV. This is declared, not fail-closed, because the domain is set by decision 0043 §5. The proton b-slab tables must cover p1v1_eq (≥ 3 GeV); otherwise the run fails closed | every case raises `UnsupportedCombinationError` (or `valid=False` plus a named counter); fields present **Dated note 2026-11-01 (implementation, src/ionmc/physics/exact_omega.py e8427d23…, config.py 2ead8d09…):** `fdm_domain_status` gives the I0 fields: the source's proton-equivalent energy is DEFINED (orchestrator decision) as the maximum p1v1_eq kinetic energy over a depth scan at x/R ∈ {0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.97} over all geometry materials, nominal at E0 and a fail-closed value at E0 + 6σ_E (p1v1_eq is nearly depth-independent, the choice moves the value by a few percent); `fdm_extrapolated` = nominal value > 300 MeV (True for He 150/250, C 290/430, O 350/430 within 2 % of the P11b energies; False for protons); an ion run whose E0 + 6σ_E value exceeds the mapping proton-path edge (4000 MeV, see I11-map) or whose scanned mapping state fails raises UnsupportedCombinationError; `summary()['mcs_mapping']` and `['physics']['ls_requested'/'ls_effective']` persist the fields (effective = off for z = 1). `validate()` still rejects ion sources; the species-enabling commit calls `fdm_domain_status` | capability | CI |
| I1 | Shared functions with a charge argument (z ∈ {1, 2, 6, 8}), the velocity-scaled arguments and the exact-Ω (b-matching) mapping functions (root solve or bounded tabulation). U1 harness: twin vs Warp CPU f32/f64; the completeness test covers every new callable | decision 0001 classes as frozen in U1 ¹⁸ | backend parity | CI |
| I2 | **He table reproduction (SEEN, regression).** Bethe 4He at I = 75 eV vs ASTAR (75 eV), and Bethe at 78 eV vs ICRU 90 alpha, water (`compare_nist.py` extended). The I-value difference between ASTAR (75 eV) and the default model (78 eV) is stated in the result | S: max abs rel dev ≤ 0.30 % for e ≥ 10 MeV/u, ≤ 1.0 % for 2–10 MeV/u. CSDA range at 50/100/150/200/250 MeV/u: ≤ 0.20 % (ASTAR) and ≤ 0.25 % (integrated ICRU 90). Aggregates only in Git. Bounds = 1.5× the SEEN decision-0038 maxima, rounded up; derived from seen values, so I2 is a regression guard, not evidence | ion-specific tabulated (regression) | LV |
| I3 | **C/O table comparison.** Bethe (78 eV) vs ICRU 73 (rev.) water for C and O, and vs ICRU 90 for C if P5(iv) establishes its origin. Source: **libdEdx ICRU73_NEW, target 276 (liquid water), Z = 6 and 8, on its 53-point grid above 2 MeV/u** (P5b; role exploratory until this row is frozen). The Geant4 G4EMLOW files are lineage only — use restricted to Geant4 by the data set's README; not an evaluation source. **I-value:** the model is evaluated at I = 78 eV (not tuned). libdEdx lists I = 75.0 eV for target 276; whether the ICRU73_NEW water values were computed with 75 or 78 eV is TO-VERIFY (P5b). If 75 eV, the known, bounded contribution of the I-value difference is ≈ 0.4–0.6 % in S above 10 MeV/u and ≈ 0.4–0.5 % in CSDA range (estimate from L0 ∝ −ln I; model S lower and range longer at 78 eV); it is stated in the result, and a range-clause failure within that bound is recorded as FAIL with the attribution, never re-graded. **Dated note 2026-11-01 (orchestrator decision (d); recorded before any I3 evaluation):** if ICRU73_NEW is NOT the revised (2009) water table, i.e. it is the pre-revision table, the I3 result stands as evaluated (same tolerances, same 78 eV model) with the attribution "reference-table version unresolved (P5b)", and it is never re-graded | S, e ≥ 10 MeV/u: max abs rel dev ≤ 1.5 % and RMS ≤ 0.8 %. S, 2–10 MeV/u: ≤ 6 % (effective-charge model spread plus PASS low-energy behaviour). CSDA range at the ends and standard points ≤ 0.5 %. Defended from theory: omitted LS ≤ 0.1 % (P3); Barkas/shell parameterisation ≤ 0.5 % (alpha experience 0.2 % SEEN, scaled for z³ terms); reference uncertainty ≈ 1–2 % (TO VERIFY). **If P5 and P5b both fail:** row = EVIDENCE GAP, intervention per protocol; it cannot pass. O has ICRU 73 only | ion-specific tabulated (largely independent of TOPAS above the confirmed 2.0 MeV scaled-energy switch (≈ 1.98 MeV/u; P2)) | LV |
| I3b | Informative: (i) Bethe with each correction off in turn (Barkas, Bloch, Mott, effective charge) vs the I3 reference, range shift per term in mm; (ii) the P3 LS-minus-(Bloch + Mott₁) curve per energy. **Dated note 2026-11-01 (round V3-008-5):** (ii) also reports the three `BetheOptions.ls` variants `"off"` (default), `"eff"` and `"bare_mixed"` against the I3 reference (S per energy and CSDA range), report-only; they never select the default (decision 0043 §1 amendment) | report only | related model | LV |
| I3c | Exploratory attribution diagnostic: "TOPAS-like" C/O water stopping rebuilt from the registered ICRU 90 proton array × q_eff² + LS ΔL, vs ionmc Bethe; S and CSDA range | report only. Used only to attribute a V1 failure; never a pass criterion; never construction. **Dated note 2026-11-01 (round V3-008-5):** the LS ΔL of the rebuild may be taken from ionmc's `BetheOptions.ls = "bare_mixed"` (LS at the nominal Z with the q² weight, Geant4's indexing) for the TOPAS-attribution diagnostic; exploratory, never construction | exploratory | LV |
| I4 | **Velocity-scaling identity (theory).** With shell/Barkas/Bloch/Mott/effective charge off: S_ion(β) = z² S_p(β) and R_ion(β) = (M_ion/m_p)/z² · R_p(β) at every grid node 2–500 MeV/u, water and tissue, with the Tmax mass term evaluated exactly | rel ≤ 1e-6 (≤ 2e-3 without the Tmax correction, reported) | theory | CI |
| I4b | Range-straggling scaling, numpy quadrature σ_R² = ∫(dσ²/dx)/S³ dE with the implemented Bohr variance, corrections off (exact algebra) | σ_R/R(ion) / [σ_R/R(p) · (m_p/M)^½] − 1 within 1e-3 at the standard points | theory | CI |
| I5 | Range-table consistency for A > 1 (V3-003D rows D1–D3 re-pointed to ion tables) | round trip ≤ 1e-5; strictly monotone; D2 distance ≤ 0.5 µm | self-consistency | CI |
| I6 | Deterministic CSDA end depth, straggling and MCS off, at the standard points and ends (CI: He 100 MeV/u only). The residual range at E_cut is reported per ion against the decision's estimates (He ≈ 0.075, C ≈ 0.025–0.03, O ≈ 0.02–0.025 mm) | R(E0) − R(E_cut) − 1e-4 R ≤ z_end ≤ R(E0) + 1e-4 R | self-consistency | CI/LV |
| I7 | Energy balance with ion sources, all EM on | rel ≤ 1e-12 (f64, python), 1e-5 (f32); counters 0 | conservation | CI |
| I8 | Range straggling in transport: σ of the end depth (straggling on, MCS off), 1e5 histories warp-cpu per standard point | within 5 % of the I4b quadrature with corrections on | self-consistency | LV |
| I8b | Report-only: σ_R/R at the standard points vs σ_R/R(p at the same range)/√(M/m_p) and vs the literature statement σ_R/R(C) ≈ 0.3–0.4 % (Schardt et al. RMP 2010, TO VERIFY) | report; outside ±15 % is recorded as a finding (not a failure) | related model | LV |
| I9 | Lateral variance of deposits (MCS on, straggling off): 1 mm slabs, 0.1 mm lateral bins, Sheppard-corrected, z/R ∈ {0.25, 0.5, 0.75, 0.9}, standard points, 1e6 histories warp-cpu | within 2 % of Fermi–Eyges A2(z) with the implemented T_ion | self-consistency | LV |
| I10 | MCS integrator step independence (U5 analogue) for T_ion: s ∈ {0.01, 0.1, 0.5, 1, 2, 5} mm, x/R ∈ {0.05, 0.25, 0.5, 0.9}, standard points; negative control: per-step Highland | \|Δθ²/θ²\| ≤ 2e-3; the control varies ≥ 5 % | self-consistency | CI |
| I11-map | **Exact-Ω mapping implementation (deterministic):** at every grid node, water and the `materials.py` set, b_p(pv_eq, p1v1_eq) = b_ion(pv, p1v1) (b = ln Ω + Δb_Fano at the respective velocities; p1v1_eq/pv_eq = p1v1/pv) is satisfied on a (pv, p1v1) test set covering the domain, by the per-step root solve, or by the tabulation if chosen (then the θ interpolation error is ≤ 0.1 %); the python twin and the Warp row agree; z = 1, M = m_p gives pv_eq = pv exactly (the proton path is unchanged) | b residual ≤ 1e-10 absolute (float64 root); twin vs Warp per decision 0001; z = 1 identity exact **Dated note 2026-11-01 (implementation):** per-step root solve as declared (no tabulation): bracketed Illinois iteration on u = ln pv_eq with bisection fallback, cap 100 iterations (observed ≤ 19), convergence |g| ≤ 1e-13 (f64) / 5e-5 (f32) or step/bracket ≤ 1e-14 / 1e-6; status codes 0 converged, 1 invalid state, 2 no sign change (b_ion above the proton b range), 3 path too short, 4 not converged (non-zero → f_dM = 0, caller treats as invalid); mapping proton path = default-options Bethe proton table to 4000 MeV (1 MeV start, 200 points/decade, 722 nodes; 4 GeV because C/O at 430 MeV/u in Be reach 3.15/3.40 GeV); the 500 MeV proton transport tables and their hash unchanged; twin generated from the same source (`python_twin`), Warp CPU f64 vs twin 1e-12 rel / 1e-14 abs (bit-identical on the compared states), f32 within analytic budgets (≈ 1e-5 in pv_eq, 1e-6 in f_dM); z = 1 identity by a host-set `ion_is_proton` flag (not a value comparison); src b uses the transport-table layout (closed-form log-log range + inverse-range table) and differs from the reference script's b by ≈ 1.5e-6 (pv_eq/pv to ≈ 5 digits; test bound 1e-4). Tests tests/ionmc/test_ion_em_mapping_src.py e11b7b65… (I11-map residual He/C/O in water and the 12 materials.py materials; src vs reference; z = 1 identity; packing vs TransportTables; table edge; flag; LS metadata; I1 completeness and twin-vs-Warp f64/f32; a cuda-marked variant). Known limits: states below the 2 MeV/u cut in Pb give status 2 (b_ion above any proton b; outside the domain); f32 cancellation of t = R(E1) − R(E) for very thin slabs (parity set uses p1v1/pv > 1.05; the kernel keeps the birth step off this path); per-step cost at transport scale unmeasured | self-consistency | CI |
| I11-abs | **Informative: absolute MCS model vs Molière–Hanson ion theory** (P1 reference, committed script): θ_model(x)/θ_MH(x) − 1, x/R ∈ {0.1, 0.25, 0.5, 0.9}, standard points and ends, reported **together with the reference's own uncertainty** on the same grid (f_dM's proton offset vs the water reference; the Pb residual, P1 validation FAILED for Pb unless cleared by P1 completion; water not validated; compound Fano rule derived, not sourced). r6: also reported per point: the θ² fraction from proton-equivalent states above 300 MeV and the maximum extrapolation in lg pv (P11b) | **report only.** Not a pass criterion: the water reference carries about 3 % uncertainty (decision 0043 §3) and the model carries a ±5 % theory-derived uncertainty, stated in metadata and docs | source-model consistency (informative) | LV |
| I11-rel | **Evaluated pre-freeze; not blind; the freeze carries it as a disclosed gate** (mapping, tests and first evaluation were produced together at 070b07f0; SEEN construction/calibration values, not qualification evidence; chronology in the P11 note). I11-rel (LV, gating): \|(θ_ion,model/θ_p,model) / (θ_MH,ion/θ_MH,p) − 1\| ≤ 2 % at every I11 point (equal velocity, i.e. the proton at E_p = E_u·A·m_p/M, and equal x/R); failure is STOP (amendment), never a tolerance change. The predictability of this ratio from the seen P1/P1b decomposition is disclosed. Grid: the exact I11 grid of the P1 driver (standard points and ends × x/R ∈ {0.1, 0.25, 0.5, 0.9}). r6: reported per point with the extrapolated θ² fraction (P11b). For C/O the equal-velocity proton (E_p up to ≈ 430 MeV) is itself above 300 MeV, so I11-rel compares two extrapolated states; it does not test the extrapolation | ≤ 2 % at every point | theory consistency (relative; calibration lineage: same Molière theory as the mapping) | LV |
| I11-ext | **Report-only (r6, SEEN pre-freeze): extrapolation of f_dM beyond its proton-qualified domain (P11b).** • θ_p,model/θ_MH,p − 1 for protons in water at E_p ∈ {100, 200, 300, 430, 600, 1000, 2000, 3000} MeV × x/R ∈ {0.1, 0.25, 0.5, 0.9} (P1 reference; `exact_omega` functions). • Per I11 point: the θ² fraction from proton-equivalent states above 300 MeV and the extrapolation in lg pv. • SEEN at r6: +3.95/+4.01/+3.80/+2.53 % at 300 MeV and +5.44/+5.30/+5.01/+3.95 % at 3 GeV, so the drift beyond the qualified edge is ≤ 1.5 pp | report only. The drift beyond the edge enters the declared model uncertainty (decision 0043 §3). Not a gate: the values are seen and the check is theory-anchored, not independent | source-model consistency (informative; SEEN) | LV |
| I12 | Transport MCS step independence (T8 analogue): exit θ_rms of a 0.5·R slab and σ at 0.9 R, s_max ∈ {0.1, 0.5, 1, 5} mm, C 290 MeV/u, 1e6 histories warp-cpu | pairwise θ_rms ≤ 0.5 %, σ ≤ 1 %; θ_rms within 0.5 % of the I10 quadrature | falsification | LV |
| I13 | R80 vs CSDA (straggling on, MCS off), 1e5 histories, standard points and ends | \|R80/R − 1\| ≤ 0.2 % | related model | LV |
| I14 | IDD step independence (T9 analogue), C 290 MeV/u, s_max ∈ {0.1, 0.5, 1} mm, f_E ∈ {0.005, 0.02}, 1e6 histories each, all EM on | \|ΔR80\| ≤ 0.05 mm; IDD grouped χ² p > 0.001 and the max-statistic per amendment ¹⁷ | falsification | LV |
| I15 | Voxel-boundary scattering (reduced T14): transport voxels {1, 2} mm with half-voxel shifts, C 290 MeV/u, 1e6 histories | pairwise θ_rms ≤ 0.5 %, σ ≤ 1 % | falsification | LV |
| I16 | Primary-ion LET in a thin slab (deterministic, straggling/MCS off) vs S_w,ion(E_mid) from the ion's own water row (identity recorded) | rel ≤ 1e-3 | theory/self-consistency | CI |
| T1-ion | Trajectory parity python vs warp-cpu f64, all EM on: CI K = 16 He 100 MeV/u; LV K = 256 each for C 290 and O 350 | discrete columns identical for 100 % of steps; continuous within rtol = atol = 1e-10 | backend parity | CI/LV |
| T12-ion | Statistical parity: python (4·10³, 40 batches), warp-cpu f64 (1e6), warp-cuda f32 (1e6), C 290, all EM on. Observables: IDD and lateral profiles at z/R ∈ {0.25, 0.5, 0.9}; scalars R80, total deposit, σ_lat(0.5 R). Also the f32-vs-f64 probe, and He 150 / O 350 warp-cpu vs cuda | T12 rules as frozen (⁹–²²; inconclusive = fail) | backend parity | LV/HR |
| T13-ion | CUDA chunk invariance 2¹⁰ vs 2¹⁸, C 290 | bit-identical within each precision | self-consistency | HR |
| R1-ion | Repeatability, C 290 | bit-identical tallies, counters, grids, end states (python, warp-cpu) | self-consistency | CI/LV |
| R2 | **Proton regression:** A16 nuclear-off digests and the T1 proton trace fixture bit-identical to the develop baseline (P8); `A16_INTENDED_CHANGE` stays None | bit-identical; any difference = STOP | regression | CI/LV |
| C-res | Residual energy below E_cut deposited locally and tallied; `e_cut_ion ≥ 2 ×` table floor enforced. **Dated note 2026-11-09 (R2-STOP3):** for ion sources the residual is deposited along the table residual range (`continuous_residual_v1`, decision 0043 §4 amendment) and tallied in `cutoff` (the part beyond the world in `escaped`); `point_v1` keeps the one-point convention (proton default). The floor rule is unchanged | exact | capability | CI |

## B. Validation rows (independent references)

**TOPAS case definitions** (new bundles; template `topas/proton-water-150mev-lateral-emonly`):
- Physics: OpenTOPAS 4.3 / Geant4 11.4.2, `sv:Ph/Default/Modules = 1 "g4em-standard_opt4"`, G4_WATER (I = 78 eV).
- Geometry and beam: water box 120 × 120 × 400 mm; zero-size, zero-divergence monoenergetic pencil; BeamEnergy = total kinetic energy.
- Particles: `GenericIon(6,12)`, `GenericIon(8,16)`, alpha naming per P4. **Revision 1 (2026-11-01):** He cases use `alpha` (see the He lineage note below).
- Scorers:
  - IDD at 0.1 mm depth bins;
  - 3-D dose at 0.25 mm lateral × 1 mm depth over ±30 mm (binary);
  - EM-model printout.
- Variants:
  - **matched** (large electron production cut, deltas local, as in ionmc; P4) — gating;
  - **0.05 mm cut** (as the proton template) — informative.
- Statistics: 3 seeds × 1e5 primaries per evidence case, plus one smoke per ion (execution only).
- Parts: every REF case is split into parts of ≤ 1200 s planned (P4 throughput), each a separate host call with its own seed (seed table extended by part index, section C). The EM-only runs are sized the same way once P4 measures throughput. For scale: the V3-010 carbon full-physics throughput of ≈ 38 ms/history would exceed the 1500 s per-call limit at 1e5 primaries in one call.
- Directories: `validation/reference_cases/topas/{helium,carbon,oxygen}-water-<E>mevu-emonly[-matched|-cut005][-seedN|-smoke]`.

**Lineage, declared:**
- **C/O** stopping in TOPAS (P2 FACT): `G4LindhardSorensenIonModel`, switch at 2.0 MeV scaled (≈ 1.98 MeV/u); above it the `G4BetheBlochModel` table + q_eff²·(ΔL_LS + 2·Barkas·(z−1)/z), no separate Bloch/Mott; ICRU 73/90 ion data only ≤ 2.0 MeV scaled, 2 < Z ≤ 80. That the BetheBloch table is the ICRU 90 proton array at equal velocity is FACT (research note), partly confirmed (pin `G4BetheBlochModel.cc`). Shared with ionmc: Bethe theory, I = 78 eV, the Barkas F table.
- **He** stopping in TOPAS = the ICRU 90 alpha array: FACT (research note); P2 OPEN pending `G4EmBuilder.cc`/`G4BetheBlochModel.cc`. Planned as the same data as I2's reference, so the He range in V1 is NOT independent of I2. Alpha never reaches the ICRU 73 ion data (P2 FACT). **Revision 1 (2026-11-01; P4, orchestrator decision):** the He V1/V1-held cases use the particle name `alpha`. Geant4/TOPAS maps `GenericIon(2,4)` to `alpha` (P4 verbose smokes), so He always uses the alpha process set: ionIoni with BraggIon below 7.9452 MeV and `G4BetheBlochModel` above, UrbanMsc; never LindhardSorensen. The He lineage is therefore the BetheBloch lineage. Whether its table is the ICRU 90 alpha array stays "research note, partly confirmed" (FC-2), and the He range in V1 stays treated as NOT independent of I2 (unchanged; the conservative reading).
- Straggling: independent implementations. MSC: Urban (inferred) is structurally different from ionmc's model (no charge/β term in its log; truncated core + tail; step-length dependent); TOPAS lateral widths for C/O are expected 10–20 % narrower (scratch estimate, P2 §6.3, not frozen).
- The r3 pre-declared trigger (ion tables above 2.5 MeV/u; the confirmed switch is 2.0 MeV scaled, ≈ 1.98 MeV/u) is NOT met (P2): I3 keeps "largely independent of TOPAS above ≈ 2 MeV/u"; no amendment.

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| V1 | **[Superseded 2026-11-07 by revision 2 (TD-EM; tolerances (a)–(d) carried over unchanged); never run, no seed consumed.]** **TOPAS EM-only IDD (matched variant)** at the standard points (blind). ionmc: warp-cpu f64, 1e6 histories, all EM on | (a) \|R80_ionmc − R80_TOPAS\| ≤ 0.35 % of R80 and ≤ 1.0 mm (TOST; σ from seeds and batches). Defence: Bethe-78 vs integrated ICRU 90 proton range differs by up to 0.18 % (decision 0038, SEEN); effective-charge, LS and Barkas construction differences ≤ 0.1–0.2 % (estimate). (b) Distal width (R80 − R20) within 10 % (TOST); straggling theory spread is a few %, and Geant4 fluctuations are an independent implementation of Bohr-type theory. (c) Peak position within 0.5 mm. (d) Peak/entrance ratio within 3 % (TOST). For helium, V1 tests transport mechanics and straggling only; its range part is bounded by I2 and adds no independent stopping evidence. | independent MC (C/O: partly shared stopping lineage; He: shared stopping data) | REF + LV |
| V1-held | **[Superseded 2026-11-07 by revision 2 (held-out set TD-H/TL-H at He 100, C 150, O 200; same role); never run, no seed consumed.]** As V1 at **He 220, C 430, O 430, C 100 MeV/u**, run ONCE, only after V1 has been graded and after any model amendment | same criteria | independent MC (held-out) | REF + LV |
| V2 | **[Superseded 2026-11-07 by revision 2 (TL-EM, TL-FULL; the Case C finding (Urban for ions) is the basis of the TL band); never run, no seed consumed.]** **TOPAS EM-only lateral (matched variant):** Sheppard-corrected σ of the lateral dose profile in a ±10 mm window (window rule as T15 ¹⁰) at z/R80 ∈ {0.5, 0.9}, standard points. **Pre-declared branches, decided by the P2 follow-ups and the P4 printouts before freeze:** "If Urban: V2 informative (Case C), the raw ratio and expected direction reported." / "If WentzelVI/single scattering: V2 is decided by an amendment before any ion TOPAS run." Expected direction if Urban: TOPAS narrower (scratch magnitude 10–20 % for C/O, not frozen). The r3 Fermi–Eyges δ_pred and the r4 Urban-emulator extension are withdrawn; an extension can return only by a reviewed amendment before any ion TOPAS run | Case C: report only (raw ratio, \|z\|, direction). Any graded form only by the amendment of the second branch. V2 asserts at most consistency of the transported MCS with an independent MC, **not** an independent test of the z-dependence. V2 has no held-out analogue; a model change after V2 leaves lateral MC evidence as a declared gap. **Revision 1 (2026-11-01): branch decided — Case C.** P2 F5 (source: GenericIon and alpha get `G4hMultipleScattering` with no model set, hence Urban, and no `G4CoulombScattering` in opt4) and the P4 printouts (msc UrbanMsc for alpha and GenericIon) select the first pre-declared branch. V2 is informative: raw ratio, \|z\| and direction, report only. No printout contradicted it; no amendment | independent MC (informative if Urban) | REF + LV |
| V3 | **[Dated note 2026-11-07 (revision 2): stays recorded as an EVIDENCE GAP under revision 1 (history, not relabelled). Superseded as the validation target by the operator directive of 2026-11-07 (TD-EM/TD-FULL: agreement with the specified TOPAS configuration); experimental validation is deferred.]** **Measured range (C, O): EVIDENCE GAP.** The only open candidate, Sihver, Schardt, Kanai 1998 (HIMAC/GSI; J-STAGE scan), **FAILED** P6: the measured Bragg curves are normalised horizontally to the calculated peaks (no absolute measured R80/peak positions), the data are figures only, and the oxygen beam is 18O. No absolute measured R80/peak reference remains | EVIDENCE GAP → intervention path (inaccessible information: GSI Schardt 2007 precision Bragg curves, HIT Tessonnier 2017 / Krämer 2016 / Kurz 2012 numeric tables), bundled with the lateral gap if V5-SS is withdrawn. **Dated note 2026-11-01 (orchestrator decision (e)):** the measured-range gap is bundled into the intervention filed at the freeze, whatever V5-SS shows (one request, two unavailable inputs; section "Intervention at the freeze"). He: no candidate. The row cannot pass; it provides no range or lateral evidence | measured (gap) | — |
| V5-SS | **[Dated note 2026-11-07 (revision 2): kept unchanged as a gating complementary row (shared screened-Coulomb lineage); it is the only row resolving b(z) at the 3 % level, which the Urban-based TOPAS lateral rows cannot (r2-topas-target.md §5.3); cost 2.94 h within the REF budget.]** **Held-out lateral row (single scattering).** Run ONCE, after freeze, only if the run condition holds: GenericIon coverage by the SS constructor (P10: FACT, closed) and a scorer-free smoke confirming the OpenTOPAS 4.3 alias, the process list and a cost ≤ 4 h REF for this design. **r6, condition (iii):** the Coulomb-model lineage of both sides is recorded (P10 §7). The proton path is read (FACT). The ion path (`G4IonCoulombScatteringModel`, `G4ionCoulombCrossSection`) is read before freeze, and the branch is decided then. • **Branch A:** the ion model has the proton path's screening form and constant, its electron-scattering treatment (target electrons, angle limited by the energy-transfer cut) and its form factor. This design is kept. • **Branch B:** any of these differs. The proton is then not a valid z = 1 reference. The freezing revision either replaces the statistic by ion/ion ratios at a common velocity, all on the ion model (e.g. He, C, O at 150 or 250 MeV/u), with a tolerance justified then, or withdraws the row. • The smoke also archives the single-scattering model name per particle. Design (P10 §3): G4_WATER slab 1 mm, `g4em-standard_SS`, matched large electron cut (same for ion and proton), pencil beam, primaries scored by direction at the exit surface; cases He 150, C 290, O 350 MeV/u and protons at the matched velocities (T_p = (T/A)·m_p/u). **r6:** plus the end cases He 250, C 430, O 430 MeV/u and protons at their matched velocities (cases 28–33). They probe the most extrapolated f_dM states: p1v1_eq ≈ 0.8–2.5 GeV, against 0.46–1.9 GeV at the standard points (P11b). They are dropped first if the smoke-measured cost of the full design exceeds 4 h; 2e4 primaries per case in parts ≤ 1200 s. Statistic: ratio_ref = θ_SS,ion/θ_SS,p and ratio_model = θ_ionmc,ion/θ_ionmc,p at equal velocity and slab thickness, θ = projected rms in a ±2.5 θ window anchored on the ionmc prediction (fixed before the run), x and y pooled; ionmc: same slab, warp-cpu f64, ≥ 1e5 histories, same estimator; σ by bootstrap over primaries and parts. Seeds: TOPAS case numbers 22–27 (He 22, C 23, O 24, p@He-velocity 25, p@C-velocity 26, p@O-velocity 27; reassigned from the V3 full-physics cases, which never ran) and 28–33 (r6 end cases, §C), parts 0–3, seed = 20619000 + 100·part + case; ionmc seeds from the lv8 row block (r_index declared in `steps_v8.py`). Lineage: independent of Urban, of f_dM and of the Highland θ0; possibly shares the Molière-type screening (α′ term) with the mapping's Born factor (Coulomb-model lineage OPEN, P10), so V5-SS tests the mapping's realisation of b(z, M) and the angular integration rather than an independent screening theory; tests only the thin-slab state (b ≈ 9–11), declared limitation. **r6 lineage (P10 §7; FACT for the proton path):** • `G4WentzelOKandVIxSection` screens with screenZ = min(1.13 Z, 1.13 + 3.76 (zZα)²/β²) times a Thomas–Fermi-type radius, carrying a Geant4 factor 1 + exp(−0.001 Z²) and no Born term for hydrogen. This is the Born screening form of the mapping's χ_a. • It scatters off atomic electrons, with weight Z and the angle limited by min(cut, T_max). • It applies an exponential-type nuclear form factor (default type OPEN, `G4EmParameters.cc`) and a spin-½ factor. • The ion path is OPEN. • The thin-slab b is 10.0–11.9 (P11b). **Dated note 2026-11-01 (continuation of round V3-008-3 finding 3; P10 §7.4.1, §7.5.1). Condition (iii) is decided: BRANCH B.** • The ion path (`G4IonCoulombScatteringModel` + `G4IonCoulombCrossSection`, FACT) differs from the proton path in every item branch A required: screening radius (universal a0(z^0.23 + Z^0.23), no 1 + exp(−0.001 Z²)), Born factor (corr^0.04, also for H, no cap), electrons (none), form factor (none). • With the r6 design, an exact Molière–Hanson ionmc would deviate by +3.9 to +5.1 % (scratch, SEEN pre-freeze). So the proton is not a valid z = 1 reference. • **Design (supersedes the r6 cases and statistic above):** – cases He 150 (22), C 150 (25), O 150 (26), He 250 (28), C 250 (27), O 250 (31), C 430 (29), O 430 (30); 1 mm G4_WATER, `g4em-standard_SS`, 2e4 primaries per case; – the proton cases are dropped; cases 23, 24, 32 and 33 are retired unused (§C); – statistic ratio_ref = θ_SS,A/θ_SS,B and ratio_model = θ_ionmc,A/θ_ionmc,B at the same MeV/u (equal velocity; z/A = ½ for all three ions, so the ratio measures the b-dependence only); – gated pairs He/C 150, He/O 150, He/C 250, He/O 250 and C/O 430; – estimator, window, ionmc side and bootstrap as above. • **Electron production cut:** energy cut ≥ T_max of the fastest case (≈ 1.16 MeV at 430 MeV/u, kinematics), the same for all cases and checked in the smoke printout. No δ rays are produced, and the reference has no electron scattering at all. • **Run condition** additionally requires: the SS smoke archives `CoulombScat` with model `IonCoulombScattering` for GenericIon and alpha, no msc, the production-cut table, and the cost. Planning cost ≈ 1.1 h at 2 µs/event (≈ 3.3 h at 6 µs). Drop order if the full design exceeds 4 h: C/O 430 first, then the 250 MeV/u triplet. The minimum design is the 150 MeV/u triplet; if it exceeds 4 h, the row is withdrawn. • **Absolute θ_ionmc/θ_SS per ion** is report-only. The pre-registered expected reference-model factor is +3.1 to +4.4 % (P10 §7.5.1) | \|ratio_model/ratio_ref − 1\| ≤ 3 % at each ion, with the bootstrap σ of the ratio ≤ 1 % (else inconclusive = FAIL). Tolerance pre-registered now and justified: Molière–Hanson relative uncertainty (I11-rel bound 2 %; Born-screening, z-independent Fano D and nuclear size partly cancel in the ion/proton ratio) and model differences of the SS reference in the ratio (≲ 1 %) combined with the statistics of a ratio of two widths (2 × 0.5–0.7 % ≈ 0.7–1.0 %) give ≈ 2.4–2.6 %, rounded to 3 %. **r6: status of the 3 %.** It is provisional and cannot be frozen until the branch is decided. • **Branch A:** the cancellation argument is re-derived from the proton-path source. At equal velocity the screening constant, the electron share and the position of the form-factor cutoff relative to χ_a are the same for the ion and the proton, apart from the Born factor. The residual terms are: • the Born-factor shift of the form-factor logarithm (≲ 0.3 % in θ, estimate); • the electron part: Geant4 screens it with the Born-modified screenZ while ionmc's Fano term is z-independent (≈ 13 % electron share × ≈ 0.65 log shift for O ⇒ Δb ≈ 0.08 ⇒ ≲ 0.5 % in θ, estimate); • the proton's spin-½ term (< 0.1 % at these angles); • the hydrogen branch without a Born term (< 0.1 % in θ, estimate). The ≲ 1 % model term therefore holds, and 3 % stands. • **Branch B:** a mismatch in electron scattering alone would be ≈ 6 % in θ in water (the electron share of χ_c² is ≈ 13 %), larger than the tolerance. The absolute θ_ionmc/θ_SS is reported (informative). **Dated note 2026-11-01 (continuation): branch B. Tolerance set now from the source analysis, before any V5-SS run or SS smoke.** |ratio_model/ratio_ref − 1| ≤ 3 % at each of the five gated ion/ion pairs (all must pass), with bootstrap σ of the ratio ≤ 1 % (else inconclusive = FAIL). Budget (P10 §7.5.1; [inserted 2026-11-01, round V3-008-4 finding 4:] computed by the committed, tested module `validation/scripts/ion_em/ss_budget.py` (b893c0420bb492368d89943a8bd8820752ae903375921f9005c7dd2bc8f00156) and `tests/ionmc/test_ion_em_ss_budget.py` (551c75838dbedbf2be0a7b643beb1cd4b5ce48d2d15e15f0f4883bf00be1cfcc), P10 §7.8, which reproduces every figure below unchanged): • (1) reference-model bias, computed −0.17 to −0.69 % (radius −0.15 to −0.62, corr −0.06 to −0.19, missing electrons +0.07 to +0.59, window estimator −0.05 to −0.41), plus ≤ 0.2 pp for the approximations: ≤ 0.9 %; • (2) ionmc against Molière–Hanson in a common-velocity ion/ion ratio = (1 + r₁)/(1 + r₂), r = I11-rel per ion: ≤ 2 % by the I11-rel gate and the SEEN positive sign (SEEN spread ≤ 1.13 pp; thin-slab state assumed like the I11 grid, declared); • (3) statistics 0.7–1.0 %. Quadrature ≈ 2.4 %; conservative (bias linear) 3.0–3.1 %. **3 % kept.** Expected value of the statistic: −0.2 to −0.7 %. Discrimination: "b without z" fails every He pair (+6.4 to +8.4 %); "no Born term" fails He/O (−4.0 to −4.8 % with the bias) and sits at the tolerance for He/C. C/O 430 (signal ≈ −1 %) is a consistency check only. **If the run condition fails:** the row is withdrawn and becomes the declared EVIDENCE GAP with the lateral-gap intervention (decision 0043 §3), bundled with V3. **Dated note 2026-11-01 (round V3-008-5, Codex REVIEW-6a14f4ec finding 2; before any V5-SS run or scored smoke; P10 §7.5.2; numbers SCRIPT-GENERATED by `ss_budget.py` §(h), SEEN planning values).** • **Rule (supersedes "|ratio_model/ratio_ref − 1| ≤ 3 % … with bootstrap σ ≤ 1 %" above):** per gated ratio PASS iff \|ratio_model/ratio_ref − 1\| + 1.645 σ ≤ 3 % (the plan's scalar TOST form, § Statistics), σ = bootstrap SE of the ratio statistic over primaries and parts; σ > 1 % = INCONCLUSIVE (counts as FAIL); all five gated ratios must pass. The round V3-008-4 module rule, which passed e.g. 2.99 % at σ = 0.9 %, is withdrawn. • **Statistics:** relative SE of the windowed rms at 2e4 primaries 0.34 % per case (x and y pooled; window kurtosis 2.76–2.78, x/y correlation 0.04–0.05), ionmc at 1e5 0.14 %; σ of the ratio 0.53 % at 2e4 / 1e5. Power at 2e4 / 1e5 (true deviation = computed bias, five correlated ratios, Monte Carlo): P(all pass) 0.994 for a correct model, but 0.254 with ionmc's ratio off Molière–Hanson by the SEEN I11-rel spread (1.13 pp, adverse). • **Design (pre-registered now): TOPAS 6e4 primaries per case (≤ 2 parts of ≤ 1200 s; parts 0–1 of the seed formula), ionmc ≥ 4e5 histories per case.** σ of the ratio 0.30 %; P(all pass) ≥ 0.9999 (correct model), 1.000 (bias − 0.2 pp), 0.918 (SEEN model spread, adverse), 0.006 (model term at the I11-rel edge 2 pp, adverse: a finding, not a false fail). REF cost scaled from the 1 mm C 290 smoke: 10 596 s ≈ 2.94 h ≤ 4 h. If the measured cost of this design exceeds 4 h, N falls back to 4e4 (1.96 h; σ 0.36 %; P 1.000 / 0.751) before the drop order above. • **Tolerance 3 % and the five gated ratios kept.** Under TOST the tolerance covers the systematic items only (bias ≤ 0.9 % linear; room left for ionmc's model term ≈ 1.6–1.8 pp ≥ SEEN 1.13 pp); 3.5 % not adopted (not needed at this N; it would let "no Born" pass He/C 250). Module and test hashes: P9 (iv) | r5: independent MC (single scattering; held-out). **r6:** held-out single-scattering MC. • Independent of the multiple-scattering approximation: f_dM, the Highland θ0, the Molière angular integration and the b-matching. • NOT independent of the screened-Coulomb cross-section physics: the proton path has the same Born screening form; the ion path is OPEN. • It is evidence that the mapping realises the shared cross-section model, not independent physical evidence of the charge dependence of screening. **Dated note 2026-11-01 (branch B; P10 §7.5.1):** held-out single-scattering MC on the Geant4 ion path. • Independent of the MCS approximation. • Shares the screened Rutherford law and Molière's Born form (modified by corr). • Differs in the screening radius (universal, z-dependent) and the electrons (none); both are pre-computed as bias. • Evidence that the mapping realises the b(z)-dependence of the collision number in the thin slab. It is not independent physical evidence of the charge dependence of screening, and covers the f_dM extrapolation only differentially **Dated note 2026-11-01 (SS smokes, P10 row):** run conditions (i) alias and (ii) process list hold (`g4em-standard_SS`: `IonCoulombScattering` for GenericIon AND alpha, `eCoulombScattering` for the proton, no msc for any); since Geant4 maps `GenericIon(2,4)` to `alpha`, the He cases use the alpha process set, which in the SS constructor has the SAME `IonCoulombScattering` model as C/O (the one-ion-model premise of branch B holds for the scattering; only the ionisation model differs, BetheBloch vs LindhardSorensen, which does not deflect the primary with the 1 km electron cut). Every V5-SS case sets `Ts/MaxStepNumber` ≥ 1e7 and fails closed on any killed track (the 300 mm carbon SS smoke lost all tracks at the 1e6 default); the per-case cost is re-estimated from the 1 mm carbon SS smoke before the freeze (§E) | REF + LV |
| V4 | **[Superseded 2026-11-07 by revision 2: the FULL configuration carries the 0.05 mm cut; its delta-electron effect is covered by the generation-0 lateral observable and the report-only all-particle core (TL-FULL); never run.]** Informative: ionmc vs the TOPAS 0.05 mm-cut variant (delta transport), R80 and lateral σ | report only | independent MC | LV |


### Revision 2 (2026-11-07): TOPAS-agreement rows (operator directive; design and defence in `validation/plans/v3-008/r2-topas-target.md`)

Reference configurations (bundles under `validation/reference_cases/topas/*-r2-*`): **EM** (intermediate) =
`sv:Ph/Default/Modules = 1 "g4em-standard_opt4"`, `CutForAllParticles` 0.05 mm, `CutForElectron` 1 km;
**FULL** (final, the specified TOPAS configuration) = `sv:Ph/Default/Modules = 7 "g4em-standard_opt4"
"g4h-phy_QGSP_BIC_HP" "g4decay" "g4h-elastic_HP" "g4stopping" "g4ion-QMD" "g4radioactivedecay"`,
`CutForAllParticles` 0.05 mm. Both: G4_WATER, box 120 × 120 × L mm, pencil beam, He as `alpha`, `EMRangeMax` 10 GeV,
`Ts/MaxStepNumber` 1e7. Cases (water): main He 50/150/250, C 80/290/430, O 100/350/430 MeV/u (standard points and
ends; low ends in the f_dM domain, the rest extrapolated, P11b); held-out He 100, C 150, O 200. IDD runs and lateral
runs are separate TOPAS runs. **Lineage: shared-lineage comparisons for stopping** (TOPAS C/O = ICRU 90 proton array
× q² + LS, He = ICRU 90 alpha array; I = 78 eV and the Barkas F table shared with ionmc); MCS independent but
structurally deficient for ions in TOPAS (Urban: |charge| in the prefactor only, no z or β in the log, step
dependent; P2). All rows TOST (§ Statistics) with σ_comb from the four TOPAS parts (SD/√4 × t_{0.95,3}/1.645 = 1.430)
and the ionmc batches; inconclusive = FAIL. Observables, estimators and resolvability rules: r2-topas-target.md §4.

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| TD-EM | **TOPAS EM IDD, main set (9 cases).** TOPAS 4 × 2.5e4 primaries per case (IDD bins 0.02/0.1 mm); ionmc warp-cpu f64, 1e6 histories, 40 batches, same bin edges | per case, all of: (a) \|ΔR80\| ≤ max(min(0.35 % R80_T, 1.0 mm), 0.05 mm); (b) distal width R20 − R80 within 10 %, gated where W_T ≥ 1.0 mm (else report-only); (c) peak position within 0.5 mm; (d) peak-to-entrance ratio within 3 %; (e) plateau shape: laterally integrated energy per depth in groups (5/2/1 mm by range) from the entrance window to 0.8 R80_T, each group within 2 %. Defence: r2-topas-target.md §5.1 (stopping-construction difference ≤ 0.18 % C/O / ≤ 0.10 % He, SEEN decision 0038; ion terms ≤ 0.15 %; shared I-value, so neither the ±0.3 % I-value nor the 1–2 % table uncertainty enters) | shared-lineage MC agreement (stopping shared; straggling, end of range and transport independent) | REF + LV |
| TD-FULL | **TOPAS FULL IDD, main set (final).** TOPAS 4 × 5e3 per case; ionmc as TD-EM (EM-only) | per case: (a) \|ΔR80 (total IDD)\| ≤ TD-EM (a) tolerance + A_frag, A_frag = 0.15 % R80_T (He), 0.10 % R80_T (C, O) (= 0.25 σ_R, fragment-pedestal shift bound; SEEN proton support −0.11 σ_R); (c) peak position within 0.5 mm. **Report-only:** peak-to-entrance, plateau, distal width, tail beyond R80, absolute entrance dose, `IDDPrimary` R80, primary fluence, Δ_frag (no correction model is declared for nuclear attenuation) | shared-lineage MC agreement with a full-physics reference on EM-determined observables | REF + LV |
| TL-EM | **TOPAS EM lateral, main set.** Core σ (bin-integrated Gaussian ML fit within ±2σ, x/y projected, mean) in the 1 mm slab at z = f·R80_T,EM, f ∈ {0.25, 0.5, 0.75, 0.9}; TOPAS 4 × 2.5e4 (lateral bins 0.01/0.02 mm); point gated iff σ_T ≥ 3 h_x | (a) every gated point: r = σ_ionmc/σ_TOPAS with (r − 1) − 1.645σ ≥ −5 % and (r − 1) + 1.645σ ≤ +40 % (ionmc ±5 %; Urban deficit 10–30 %, P2 §6.3, ×1.04 low-β term); (b) **f_dM-extrapolation trend** per ion: \|D_ext\| + 1.645σ ≤ 12 %, D_ext = ⟨r⟩(upper end)/⟨r⟩(low end) − 1 over the gated depths in {0.75, 0.9} common to both (NOT EVALUABLE if none; recorded, not a pass). This is the criterion by which the f_dM extrapolation's adequacy is judged (with TL-FULL (b)) | independent MSC implementation, structurally deficient reference (coarse band; not a precision test of the charge dependence) | REF + LV |
| TL-FULL | **TOPAS FULL lateral, main set (final).** As TL-EM on the generation-0 (primary) projected profiles of the FULL run; TOPAS 4 × 5e3; same z_k | (a) and (b) as TL-EM. Report-only: all-particle core σ, halo fraction, ±20 mm second moments, σ_FULL/σ_EM | as TL-EM | REF + LV |
| TD-H / TL-H | **Held-out set** He 100, C 150, O 200 MeV/u (cases 10–12), all four configurations; run ONCE, only after the main TD/TL rows are graded and after any model amendment | TD-EM (a)–(e), TD-FULL (a), (c), TL-EM (a), TL-FULL (a) (no D_ext) | held-out shared-lineage MC agreement | REF + LV |
| TF-DIAG | Report-only TOPAS-internal diagnostics (no ionmc output): Δ_frag = R80_FULL − R80_EM, primary survival, FULL/EM plateau and peak ratios, σ_FULL,primary/σ_EM, halo fractions | report only; attribution rule in r2-topas-target.md §5.5 and STOP 5 | reference-internal | LV |

**Order and blindness (revision 2; revised 2026-11-07, round R2-TOPAS-3, Codex REVIEW-4c5b6803: blinding by ORDERING,
not by output suppression).** The controlled runner archives every TOPAS output permanently, smokes included, so
nothing is deleted or declared unread. Step 0: implementation and the A-rows (T1/T12 failures STOP first), then the
ionmc TD/TL samples of every main and held-out case are produced, their sealed records and sha256 committed in a
reviewed commit, the **sample seal commit** (its SHA and review id are recorded here by a dated seal note before any TOPAS
run). Step 1: timing and step-length smokes (representative, configuration-specific scorers identical to the
evidence configuration; outputs archived by the runner and read for sizing only after the seal; they disclose nothing
that can influence the sealed samples). Step 2: dated sizing note. Step 3: evidence runs, then TD-EM/TL-EM graded →
TD-FULL/TL-FULL graded → held-out set once (`r2-topas-target.md` §6 order). **No TOPAS run of the revision-2 set,
smoke or evidence, may start before the sample seal commit.** No ionmc parameter is fitted to TOPAS; no ionmc
prediction is recomputed after a TOPAS output of the same row has been seen. **Run validity** per TOPAS part: zero
killed/stuck tracks; the EM printout shows UrbanMsc for alpha/GenericIon and the bundled module list; the replicate
gate holds; a violating part is re-run at a re-run base declared by amendment and preserved. **Run validity against
the seal (R2-TOPAS-4, git ancestry and controlled records):** a revision-2 TOPAS run (smoke or evidence) is valid only if (i) the `code_sha` of its run record (the clean committed SHA the controlled runner executed) has the sample seal commit as an ancestor (`git merge-base --is-ancestor <seal> <code_sha>`), and (ii) the seal commit's independent review (review id in the dated seal note) finished before the run's controlled `reference_started` event (controlled-service records, never Git timestamps; `finished_at` of the run is not used). The seal commit is reviewed and pushed before any run (a run at an unpushed SHA is impossible: the runner snapshots the committed SHA). The evaluator `validation/scripts/ion_em/r2_evaluate.py` (implemented with the TD/TL rows) enforces exactly these two checks and fails closed when a record is missing; a run failing either check invalidates the TD/TL rows of its case, recorded as a FAIL of the blinding control and never re-graded.

**Ledger wording (revision 2):** achieved ion validation is described as **agreement with the specified TOPAS
configuration** (OpenTOPAS 4.3 / Geant4 11.4.2, the seven-module list; EM-only as intermediate), shared stopping
lineage stated; experimental validation is deferred; the lateral agreement is within the pre-registered coarse band;
the charge dependence of lateral spreading has shared-lineage single-scattering evidence (V5-SS) only.

**Dated seal note (2026-11-09, orchestrator; revision 2, plan section "Order" step 0 — the sample seal).** The
**sample seal commit** is `0c68b78ac191c850bf9101dbacf8eca12f4e9b1b` (branch
`task/v3-008-helium-carbon-and-oxygen-electromagnetic-transport-effective`, pushed 2026-11-09 ≈09:21 UTC, PR #84). It
carries the sealed record `validation/results/transport/v8r2-samples/sample-seal.json` (committed in its ancestor
`84b8b2cd`; schema `ionmc-v8r2-sample-seal-1`; 86 chunks of the twelve 1e6-history ionmc TD/TL samples produced by 76
signed host-runner runs at the frozen implementation SHA `1831f2149506b0bd403259f2d00981b835837e4d`; every step `pass`,
every diagnostic counter 0; one superseded and one failed run preserved) and the two review fixes that followed the record
(`f2141b03`: mapping-domain rule D1; `0c68b78a`: no mapping work with multiple scattering off; both leave full-coverage
results bit-identical, P9 freeze record amendment of 2026-11-09). Independent review of the seal commit: Codex
`REVIEW-9691e1b098e940d7bdfd1063c09325ee`, `head_sha` 0c68b78a…, status `passed`, verdict `pass`, `finished_at`
2026-11-09T09:15:57.707596+00:00 (preceded by REVIEW-16c63514 and REVIEW-f3361441, `changes_required`, both resolved).
Validation record `VAL-20261109-092143-547E49` (suite `v3-008-r2-sample-seal`, passed, `details.host_runs` = the 76 runs).
CI run 42592861874 at 0c68b78a: success (jobs docs, pre-commit and tests all success; completed 2026-11-09T09:30:30Z). Evaluator binding (`r2_evaluate.py`, plan record):

```json
{"seal_sha": "0c68b78ac191c850bf9101dbacf8eca12f4e9b1b", "seal_review_id": "REVIEW-9691e1b098e940d7bdfd1063c09325ee"}
```

No TOPAS run of the revision-2 set (smoke or evidence) exists at the time of this note; every such run must have a
`code_sha` descending from the seal commit and a controlled `reference_started` later than the review's `finished_at`
above (R2-TOPAS-4). The TD/TL statistics of `r2_evaluate.py` are implemented after this seal (plan section "Order").

## C. Seeds (fresh; registered before any run; shared seed-base ledger)

**Shared seed-base ledger.**
- One ledger of bases across tasks. The V3-005 sequence continues at **20641004**.
- **20572004–20632003 is reserved for V3-008.** This is the committed ledger on the V3-005 branch: commit
  06974d07, the dated note to Amendment 20 (f)2, and the guard constant
  `V3_008_RESERVED = range(20572004, 20632004)`. The block starts at the first value above the existing base
  envelopes; base 20551004 reserves up to 20572003. P7 confirms this note.

**V3-008 bases (the r2 fallback set, now declared):**
- lv8 qualification **20581004**; hr8 qualification **20591004**.
- Row seed `base + 100·r_index + k` (k < 100, r_index ≤ 60, fixed in `steps_v8.py`). Hence lv8 seeds lie in
  20581004–20587103 and hr8 seeds in 20591004–20597103.
- Rehearsals use only the family **2060xxxx** (20600000–20609999). They are never evidence.
- **TOPAS seeds:** `20619000 + 100·part + case`, with case 1–40 and part 0–19. Part 0 gives 20619001–20619040
  (the r2 table); the highest is 20620940.
  - Case numbers: V1/V2 matched 1–9; V1-held 10–21; V5-SS 22–27 (reassigned from the V3 full-physics cases, which never ran; V3 is an evidence gap); V5-SS end cases 28–33 (r6: He 250 = 28, C 430 = 29, O 430 = 30, p at the He 250 velocity = 31, p at the C 430 velocity = 32, p at the O 430 velocity = 33; previously unused, inside the declared formula range, no ledger change); cut variants 34–40.
  - **Dated note 2026-11-01 (V5-SS branch B; before any run):** case 25 = C 150, 26 = O 150, 27 = C 250, 31 = O 250 MeV/u (were protons at matched velocities). Cases 22 (He 150), 28 (He 250), 29 (C 430) and 30 (O 430) are unchanged. Cases 23 (C 290), 24 (O 350), 32 and 33 are retired unused. Nothing has run, so no seed is consumed. The formula and range are unchanged, so there is no ledger change.
  - Each part (≤ 1200 s planned) is a separate host call with its own seed.
  - **Dated note 2026-11-07 (revision 2; before any run; formula, ranges and guard test unchanged, so no ledger change).** Case numbers re-assigned (none of the old numbers ever ran; no seed consumed): **1–12** = the revision-2 TOPAS-agreement cases He 50 (1), He 150 (2), He 250 (3), C 80 (4), C 290 (5), C 430 (6), O 100 (7), O 350 (8), O 430 (9) and the held-out He 100 (10), C 150 (11), O 200 (12) (previously V1/V2 matched 1–9 and part of V1-held). Part blocks per configuration: IDD-EM parts 0–3, LAT-EM 4–7, IDD-FULL 8–11, LAT-FULL 12–15; parts 16–17 (IDD-FULL) and 18–19 (LAT-FULL) reserved for the pre-registered FULL sizing fallback (r2-topas-target.md §6). V5-SS keeps cases 22, 25–31 (parts 0–1). Retired unused and kept retired: 23, 24, 32, 33. Unassigned (use only by amendment): 13–21, 34–40. Revision-2 evidence seeds: 192 values, 20619001 (He 50 IDD-EM part 0) to 20620512 (O 200 LAT-FULL part 15), distinct from each other and from V5-SS. Execution-only smokes of revision 2 use the rehearsal family: 20600411–20600422 (timing) and 20600423 (step length); 20600401–20600410 are consumed by the P4/P10 smokes. **ionmc TD/TL samples:** lv8 rows r_index 30–41 = cases 1–12 (declared now; `steps_v8.py` adopts them), k = 0 for the single run per case, k = 1, 2, … for chunks of a split run; the held-out set after a model amendment uses the re-run base 20621004 with the same r_index.
- **Re-runs** take fresh bases from **20621004**, then 20631004, each declared in an amendment before use. The
  last base inside the reserved block is 20631004, whose `base + 100·r + k` seeds stay ≤ 20632003 only for
  r_index ≤ 9. A re-run needing more rows, or a third re-run base, requests a new block in the shared ledger.
  - **Dated note 2026-11-11 (amendment I9-REF, item 4k; before any re-run).** The third re-run base **20572004** (seeds
    20572004–20578103, full r_index ≤ 60; family `rerun3`) is declared for the lv8 re-run of I9, I12, I14 and I15
    (r_index 44–47). It is the first value of the V3-008 reservation, in the leading gap 20572004–20581003 that no
    family had used, disjoint from lv8/hr8/rehearsal/TOPAS/re-run families and from every V3-005 family (guard test
    extended); the sentence "a third re-run base requests a new block in the shared ledger" is superseded for this
    base because the reservation itself still holds it — no cross-task ledger change. 20631004 remains the declared
    base for re-runs with r_index ≤ 9.
  - **Dated note 2026-11-13 (results of the LV/HR block; after the runs).** Base **20572004** (family `rerun3`, lv8 I9/I12/I14/I15,
    r_index 44–47) and base **20591004** (hr8 T13-ion and T12-ion HR, r_index 54) are **consumed**: their full-scale
    results were observed (all PASS; dated results note of 2026-11-13). Base 20581004 (lv8 A-rows at f5f3a701, item 4d;
    I9 r 44 and I12 r 45 shard 0 at e602a583) was already consumed. Any further lv8/hr8 re-run takes a fresh declared base.
  - **Dated note 2026-11-13 (V5-SS ionmc samples; before any run).** Base **20631004** (family `rerun2`, declared for
    re-runs with r_index ≤ 9 and never drawn) is used for the FIRST run of the V5-SS ionmc thin-slab samples, row `v5ss`,
    r_index = case order index 0–7 (cases 22, 25, 26, 27, 28, 29, 30, 31), k = chunk index; seeds 20631004–20631704 step 100;
    it is consumed once the samples' statistics are sealed; a V5-SS re-run would need a new declared block (guard test
    extended). These r_index values index base 20631004 only (`steps_v8.SS_R_INDEX`), not the lv8 table `ROW_R_INDEX`.
- Never 2098xxxx or 2099xxxx. CI tests use fixed small seeds documented in the tests.

**Guard test** `tests/ionmc/test_v3_008_seed_families.py`. It asserts that every V3-008 seed is:
- inside 20572004–20632003 (`V3_008_RESERVED`), including every declared base (lv8 20581004, hr8 20591004, the rehearsal family 2060xxxx, re-run bases from 20621004);
- for TOPAS, generated by the part formula `20619000 + 100·part + case` (case 1–40, part 0–19), whose maximum (20620940) stays below 20621004;
- disjoint from the consumed v3/develop bases (20261004–20421004, families 2041xxxx/2043xxxx);
- disjoint from V3-005's fresh and consumed bases with their 21 000 envelope (`base + 0…20999`), including V3-005's
  declared next bases from 20641004 on (20641004, 20651004, …);
- disjoint from the V3-005 TOPAS formula range `20270000 + 1000·row + k` (row 1–23, k 1–3);
- outside 2098xxxx/2099xxxx.

The V3-005 constants are copied as literals with a cross-reference, because that branch is unmerged. A base whose
full-scale statistical results were observed is consumed, pass or fail.

## D. PASS / FAIL / STOP

- **PASS (task):** every gating row of A and B passes at one clean frozen SHA with the declared
  bases, or is an EVIDENCE GAP handled under STOP 4. Deferred multiprocess rows are DEFERRED
  (`conformant: false` while the single-process directive holds).
  **Dated note 2026-11-07 (revision 2):** the gating TOPAS-agreement rows are TD-EM, TD-FULL (a, c), TL-EM, TL-FULL and, once run, TD-H/TL-H; V5-SS stays gating. "EVIDENCE GAP handled under STOP 4" no longer applies to V3 and the lateral gap: the operator's directive resolved STOP 4 (below), and V3 stays recorded as a revision-1 EVIDENCE GAP (history), not as a gating row of revision 2.
- **FAIL:** recorded verbatim in the evidence table and in decision 0043 "Later validation
  outcome". It is never re-graded with a widened tolerance.
- **STOP** (halt; amendment and independent review before continuing):
  1. R2 differs.
  2. The exact-Ω (b-matching) mapping cannot be committed with passing I11-map before freeze, or I11-rel fails (> 2 % at any I11 point): STOP with amendment and re-review, never a tolerance change. (The r3 selection STOP is superseded: no selection remains.)
  3. V1 or V2 fails. Attribution with I3b/I3c and the P2 source reading is allowed. Any model change
     after seeing V1/V2 consumes them (preserved as failed). V1-held then becomes the sole gating
     independent-MC range evidence, at a fresh base with unchanged tolerances. V2 can fail only if a
     graded form is introduced by the amendment of its second branch. No reference prediction is
     recomputed after a TOPAS output has been seen, and no parameter is fitted to TOPAS.
     **Dated note 2026-11-07 (revision 2):** STOP 3 applies equally to TD-EM, TL-EM, TD-FULL and TL-FULL (attribution with I3b/I3c, TF-DIAG and the P2 source reading allowed). A model change after seeing any of them consumes all main TD/TL rows (preserved as graded); the held-out set TD-H/TL-H then becomes the sole gating TOPAS-agreement evidence, at the re-run ionmc base, with unchanged tolerances. A TL failure caused by the reference's structure (Urban outside the band) is recorded as FAIL with that attribution, never re-graded.
  4. Evidence gap (I3 data, V3 measured, ion MCS z-dependence): an intervention per protocol. **r6:** P6's lateral search has failed (EVIDENCE GAP). V5-SS shares the screened-Coulomb cross-section physics, so it is not independent physical evidence of the charge dependence of screening. The planner recommends filing the lateral-gap intervention before integration whatever the V5-SS outcome, with the V5-SS result and its class stated in the ledger. This is an orchestrator decision, to be recorded before freeze. r5 text: for the lateral gap it is used only after P6 and P10 fail, with the wording of decision 0043 §3 verbatim ("Unavailable input: independent evidence of the charge dependence of lateral spreading for He/C/O in water at 50–430 MeV/u, … recorded as a gap in the ledger."). Integration is allowed with the gap stated in the ledger; V2-ION stays "in progress". **Dated note 2026-11-01 (continuation):** under V5-SS branch B, the f_dM extrapolation also loses its absolute independent coverage. The arguments for and against filing regardless of V5-SS, and the question whether the verbatim clause "or an accessible single-scattering MC reference" needs an amendment, are in P10 §7.6.2. This is an orchestrator decision before freeze. **Orchestrator decision 2026-11-01 (before freeze):** (a) the lateral-gap intervention is filed at the freeze whatever V5-SS shows, naming as the unavailable input the independent evidence of the charge dependence of lateral spreading AND the absence of independent coverage of the f_dM extrapolation beyond the 300 MeV proton domain (declared and bounded only within the same Molière theory); V5-SS stays in the plan as shared-lineage single-scattering evidence (branch B) and its result and class go into the ledger; (b) the verbatim clause "or an accessible single-scattering MC reference" is amended in the reviewed freeze revision so that a reference sharing the screened-Coulomb cross-section lineage does not satisfy it; (c) the f_dM extrapolation is kept as a declared limitation (uncertainty component ≤ 1.5 %, I11-ext report row, `fdm_extrapolated` flag); no further coverage is claimed. **Revision 1 amendment (2026-11-01; orchestrator decision (b)).** The verbatim clause of the lateral-gap wording (decision 0043 §3), quoted: "or an accessible single-scattering MC reference". Amended wording: "or an accessible single-scattering MC reference whose screened-Coulomb cross-section physics is independent of the model used by the exact-Ω mapping; a reference that shares the screened Rutherford law and Molière's Born screening form (including the Geant4 single-scattering configuration of V5-SS, `G4IonCoulombScatteringModel`, branch B) does not satisfy this input". The intervention is filed at the freeze with this wording and both gaps (section "Intervention at the freeze (to be filed by the orchestrator)"). V5-SS stays a gating row; its result and class go into the ledger when it runs.
   **Dated resolution 2026-11-07 (revision 2; operator directive `reference_requests/operator-directives/2026-11-07-topas-validation-target.md`, IR-20261102-010435-B59B2D resolved 2026-11-07 20:37 UTC).** The operator's answer supersedes the gap-based STOP 4 for V3-008: measured data and an independent-physics single-scattering reference are no longer prerequisites; the acceptance reference for depth-dose distributions and lateral profiles is agreement with the specified TOPAS configuration (EM-only TOPAS as intermediate validation); experimental validation is deferred. Consequences: (i) V3 (measured C/O range) stays recorded as an EVIDENCE GAP under revision 1 (history; Sihver 1998 FAILED P6) and is superseded as a target by TD-EM/TD-FULL; (ii) the lateral gap (intervention item (1)) is recorded as history and superseded as a target by TL-EM/TL-FULL; the ledger states that independent physical evidence of the charge dependence of lateral spreading is part of the deferred experimental validation and that V5-SS is shared-lineage evidence; (iii) the f_dM extrapolation (item (2)) stays a declared model-design limitation; I11-ext stays report-only; its adequacy is judged by **TL-EM criterion (b) and TL-FULL criterion (b)** (the in-domain-to-extrapolated trend D_ext, \|D_ext\| + 1.645σ ≤ 12 %), with TL (a) at the extrapolated cases as the absolute check; (iv) nothing that has run is re-graded; the P1 Pb FAILED record (cleared by attribution), the P6 EVIDENCE GAPs and the SEEN rows stay as recorded.
  5. **(Revision 2, 2026-11-07.)** TD-FULL or TL-FULL fails while the corresponding EM row passes: STOP; attribution with TF-DIAG (Δ_frag against A_frag; σ_FULL,primary/σ_EM); the result stays FAIL, is never re-graded and goes to the orchestrator; integration then needs a reviewed amendment or an operator decision. A TOPAS run that rejects any of the seven modules is a STOP before any scored FULL run (the planner cannot change the operator's configuration).
- **(Revision 2, R2-TOPAS-3, 2026-11-07.)** A revision-2 TOPAS run (smoke or evidence) fails either check of the run-validity rule (R2-TOPAS-4: the sample seal commit is an ancestor of the run record's `code_sha`; the seal commit's review finished before the run's controlled `reference_started`; evaluator `r2_evaluate.py`, fail closed on a missing record): the TD/TL rows of its case are a FAIL of the blinding control (never re-graded); STOP and report to the orchestrator.
- **T1-ion and T12-ion failures are STOP (implementation defect) before any V-row is run.**
- Schedule checkpoint: if implementation exceeds 10 working days, record the cause and re-plan; no gate changes.
- Every failure, inconclusive result and consumed base is preserved.

## Intervention at the freeze (to be filed by the orchestrator)

**Dated note (2026-11-02, orchestrator) — FILED.** The intervention was filed through `request_human_intervention` on 2026-11-02 01:04 UTC as **IR-20261102-010435-B59B2D** (category inaccessible-information; notification sent, HTTP 200) at the freeze candidate 655d16a0; the record is committed verbatim on this branch at `reference_requests/interventions/IR-20261102-010435-B59B2D/` (README.md 86eba08d…, request.json 068cafd7…; commit cc2e4e26). Options offered: A measured lateral data, B independent-physics single-scattering reference, C accept the lateral and f_dM gaps, D measured C/O ranges with beam-line offsets, E accept the V3 EVIDENCE GAP. Integration is not blocked (STOP 4); the ledger wording of V2-ION and V3 follows the operator's answer, recorded as a dated resolution when it arrives.

**Dated resolution (2026-11-07, revision 2) — ANSWERED.** IR-20261102-010435-B59B2D was resolved by the operator on 2026-11-07 at 20:37 UTC. The answer is recorded verbatim in `reference_requests/operator-directives/2026-11-07-topas-validation-target.md` (sha256 ca079e9f…; the intervention's own resolution files are handled by the orchestrator). It selects none of the options A–E as offered: it changes the validation target to agreement with the specified TOPAS configuration and defers experimental validation. Processing: plan revision 2 (STOP 4 dated resolution; §B revision 2 block; `validation/plans/v3-008/r2-topas-target.md`) and decision 0043 revision 2. Ledger wording of V2-ION and V3: as in the §B "Ledger wording (revision 2)" paragraph. The text below is kept as the filed record.


Revision 1, 2026-11-01; orchestrator decisions (a)–(c) of STOP 4. The orchestrator files this request through the structured intervention route at the freeze commit, whatever V5-SS later shows. The filed record is copied verbatim onto the task branch.

**Category: inaccessible-information** (EXPERIMENT.md permitted intervention 1). The measured data exist: the Tessonnier et al. 2017 PMB and Schwaab et al. 2011 PMB lateral profiles in water, the HIT/GSI campaign data behind them, and the data "available from the authors" of Besuglow et al. 2022. They are paywalled or unpublished as numeric tables, so the agent system cannot access them (P6). "external-data" (permitted intervention 3) would mean a new measurement on operator-held equipment, which is not what is asked. The alternative of an independent-physics single-scattering reference is offered inside the same request.

**Wording to be filed** (the r5 verbatim wording of decision 0043 §3, amended in revision 1: the single-scattering clause per decision (b), and the second gap per decision (a)):

> Unavailable input: (1) independent evidence of the charge dependence of lateral spreading for He/C/O
> in water at 50–430 MeV/u, i.e. measured lateral dose widths or angular spreads with stated
> uncertainties (candidates Parodi 2013, Schwaab 2011, Tessonnier 2017: no open numeric tables
> found), or an accessible single-scattering MC reference whose screened-Coulomb cross-section
> physics is independent of the model used by the exact-Ω mapping; a reference that shares the
> screened Rutherford law and Molière's Born screening form (including the Geant4 single-scattering
> configuration of V5-SS, `G4IonCoulombScatteringModel`, branch B) does not satisfy this input;
> (2) independent coverage of the f_dM evaluation beyond its 300 MeV proton-qualified domain
> (proton-equivalent states up to about 2.5 GeV at the standard points and upper ends), which is
> declared and bounded only within the same Molière theory (row I11-ext). Request: provide access to
> such data or such a reference, or accept that V3-008 integrates with ion lateral spreading
> supported only by theory-derived mapping and implementation checks and by shared-lineage
> single-scattering MC evidence (V5-SS; no independent evidence), and with the f_dM extrapolation as
> a declared limitation (a ≤ 1.5 % component of the ±5 % model uncertainty), recorded as gaps in the
> ledger.

**Evidence state filed with it:**
- P6: measured C/O range = EVIDENCE GAP (Sihver, Schardt, Kanai 1998 failed: curves normalised horizontally to the calculated peaks, figures only, 18O). Measured lateral data = EVIDENCE GAP: no tabulated He/C/O lateral width or angular spread in water with stated uncertainties; the Tessonnier 2017 thesis has figures and text-quoted, beamline-confounded ratios without uncertainties; Schwaab 2011 and Tessonnier 2017 PMB are paywalled; the "Phillips" ion data are of unknown provenance; Gottschalk's survey (arXiv:1204.4470 §5.3) concurs.
- V5-SS (branch B; P10 §7.5.1–7.5.2): held-out single-scattering MC on the Geant4 ion path, independent of the multiple-scattering approximation (f_dM, Highland, Molière angular integration, b-matching). It shares the screened Rutherford law and Molière's Born form (modified by corr^0.04) and differs in the screening radius and the absent electrons (pre-computed bias). Class: shared-lineage evidence that the mapping realises the b(z)-dependence in the thin slab, not independent physical evidence of the charge dependence of lateral spreading. It covers the f_dM extrapolation only differentially (≤ 1.5 pp). It has not run at the freeze; its result and class are added to the ledger after the run.
- Theory and implementation rows: I11-rel (gating; SEEN, disclosed), I11-abs (informative), I11-ext (report-only; drift ≤ 1.5 pp from 300 MeV to 3 GeV against the Molière–Hanson proton theory; theory-anchored, not independent), I9, I10, I11-map, I12, I15.
- f_dM extrapolation: declared limitation, a ≤ 1.5 % component of the ±5 % model uncertainty, persisted as `fdm_extrapolated` and the source proton-equivalent energy (I0).
- Integration is allowed with the gaps stated in the ledger (STOP 4); V2-ION stays "in progress" for ion lateral evidence.

~~**Not part of this request:** the V3 measured-range gap (row V3 bundles it with the lateral gap only if V5-SS is withdrawn). Whether it is filed separately at the freeze is an orchestrator decision; none is recorded.~~

**Dated note 2026-11-01 (orchestrator decision (e)): the V3 measured-range gap is bundled into this request** (one request, two unavailable inputs: the lateral gap as worded above, and the measured range). Same category (inaccessible-information): the GSI/HIT measured Bragg curves exist but are not accessible as numeric data (P6). Additional wording to be filed after the lateral wording:

> Unavailable input (second): measured ranges of carbon and oxygen ions in water for row V3, i.e.
> Bragg-peak and/or R80 positions with stated uncertainties, together with the water-equivalent
> thickness of the beam-line elements upstream of the water (offsets) needed to compare them with a
> calculation, at energies within 80–430 MeV/u (C) and 100–430 MeV/u (O) (candidates: GSI precision
> Bragg curves (Schardt et al. 2007), HIT numeric tables (Tessonnier 2017, Krämer 2016, Kurz 2012)).
> The only open candidate, Sihver, Schardt, Kanai 1998, is unusable (P6: curves normalised
> horizontally to the calculated peaks, figures only, oxygen = 18O). Request: provide access to such
> data, or accept that V3-008 integrates with row V3 recorded as an EVIDENCE GAP in the ledger.

Depends on the answer, in addition: whether a later amendment adds a measured-range row (operator-provided data), pre-registered before use. Row V3 cannot pass at this revision either way.

**Decisions that depend on the answer:** the ledger wording for ion lateral spreading under V2-ION; whether a later amendment or task adds a measured-data row (operator-provided data) or an independent-physics reference row, each pre-registered before use. No frozen gate of this plan changes with the answer.

## E. Host budget (single process, one heavy job at a time; planning values)

Planning throughput (proton, 150 MeV, EM): python ≈ 60/s, warp-cpu ≈ 7·10³/s, CUDA ≈ 1.5·10⁵/s.
C at 290 MeV/u has a comparable step count, so these values are used until the rehearsal measures
ion throughput.

| Rows | Estimate |
|---|---|
| CI additions | ≤ 3 min |
| P1, P3, I2, I3, I3b, I3c, I11-abs, I11-rel (numpy) | < 10 min |
| I6, I8, I13 | ≈ 3 min |
| I9, I12, I14, I15 (≈ 17 × 1e6 warp-cpu) | ≈ 41 min, steps ≤ 1500 s |
| T1-ion LV | < 1 min |
| T12-ion LV | ≈ 9 min |
| T12/T13-ion HR | ≈ 3 min + compile |
| V1/V2/V1-held ionmc samples (≈ 9 × 1e6 warp-cpu) and V5-SS ionmc thin-slab samples (12 × 1e5, r6; 8 × 1e5 under branch B, dated note 2026-11-01; 8 × 4e5 under the round V3-008-5 TOST design, cost not measured) | ≈ 27 min (r6 estimate; not re-estimated for 8 × 4e5) |
| **TOPAS (REF)** | Every REF case is split into parts of ≤ 1200 s planned (P4 throughput), each a separate host call with its own seed (seed table extended by part index); the EM-only runs are sized the same way once P4 measures throughput (the V3-010 carbon full-physics throughput ≈ 38 ms/history would exceed the 1500 s per-call limit at 1e5 primaries). 3 smokes + 1 SS smoke (scorer-free); 9 matched evidence runs; 12 held-out; V5-SS 6 cases × 2e4 in parts (≈ 1 h, P10 estimate; withdrawn if the smoke shows > 4 h) plus the r6 end cases 28–33 (≈ 1 h more at the same per-mm cost; dropped first if the full design exceeds 4 h); the 6 V3 full-physics runs are dropped (V3 is an evidence gap); 7 cut variants. Ion throughput unknown (P4). Cap 16 h total, one run at a time; materialize, parse, delete per run |
| V5-SS under branch B (dated note 2026-11-01; P10 §7.5.1) | 8 ion cases × 2e4 instead of 6 + 6 (no proton cases): ≈ 1.1 h at 2 µs/event, ≈ 3.3 h at 6 µs; ion-path events per mm He 3.6e3–5.0e3, C 1.1e4–1.6e4, O 1.4e4–1.9e4; drop order C/O 430, then the 250 MeV/u triplet. **Dated note 2026-11-01 (round V3-008-4; measured):** the 1 mm C 290 SS smoke (REF-35df82fcd7723f6a9a0c-db361e98, as relayed) took ≈ 23 ms per primary = 1.299e4 events × ≈ 1.77 µs. Scaled per case by the events per primary (`ss_budget.cost_from_measured`): 126 s (He 250) to 674 s (O 150) per case at 2e4 primaries, one part ≤ 1200 s each; total ≈ 3.5e3 s ≈ 1.0 h for the 8 cases (upper estimate if part of the 23 ms is fixed overhead). Minutes per case, about one hour in total; ≤ 4 h holds and no case is dropped. **Dated note 2026-11-01 (round V3-008-5, V5-SS TOST design; before any run):** 6e4 primaries per case: 378 s (He 250) to 2022 s (O 150) per case, ≤ 2 parts of ≤ 1200 s, total 10 596 s ≈ 2.94 h ≤ 4 h (same scaling, upper estimate); fallback 4e4 (7064 s ≈ 1.96 h) before dropping cases. ionmc thin-slab samples 8 × 4e5 (cost not measured; expected minutes) |
| **Total** | LV ≈ 1.5 h; HR ≈ 0.25 h; REF ≤ 16 h |

All steps run with `--single-process --workers 1`, thread limits 1, logged per the 2026-10-12
operator rules, and each host call ≤ 1500 s.

**Dated note 2026-11-07 (revision 2; planning values, replaced by the timing smokes; r2-topas-target.md §8).** The V1/V1-held/V2/V4 TOPAS and ionmc items above are superseded. Revision-2 REF: IDD-EM ≈ 1.8 h, LAT-EM ≈ 0.33 h, IDD-FULL ≈ 2.9 h, LAT-FULL ≈ 0.5 h (12 cases each, 4 parts per case and configuration, every part ≤ 1200 s planned; largest planned part ≈ 470 s, C 430 IDD-FULL), 13 execution-only smokes ≈ 0.25 h, V5-SS 2.94 h (unchanged): **≈ 8.7 h** (≈ 14.3 h if FULL costs 25× EM instead of 10×), within the 16 h cap; 221 REF host calls. Sizing rule and fallback order before any scored run: r2-topas-target.md §6 step 3 (no case is dropped silently; exceeding the cap is an orchestrator decision). ionmc TD/TL samples: 12 cases × 1e6 histories warp-cpu f64 with three scoring grids, ≈ 0.5–1.5 h LV (fine-tally cost unmeasured), each host call ≤ 1500 s.

## F. Gates and integration criterion

1. Pre-freeze checks closed in the r5 order of section P (steps 1–7). Decision 0043 and this plan independently reviewed (scientific
   review + Codex) and committed before any transport or qualification result (the disclosed pre-freeze computations P1/P1b and P11 excepted; they are not qualification evidence).
2. Implementation, about 5–8 days:
   - charge/mass plumbing (tables, reference loop, kernels; no literal charge 1);
   - the exact-Ω (b-matching) mapping: per-step root solve on cumulative slab-integral tables, or a bounded (pv, p1v1) tabulation;
   - ion E_cut;
   - species append;
   - fail-closed rules;
   - ion water row for LET;
   - effective-charge docstring rename;
   - `compare_nist.py` extension;
   - C/O table loader (P5);
   - TOPAS bundles;
   - `steps_v8.py` and the run_suite `lv8`/`hr8` suites;
   - the seed guard.
   - **(revision 2, 2026-11-07)** the TD/TL analysis (estimators of r2-topas-target.md §4: R80/R20/peak/plateau groups, the core-width ML fit, D_ext, TOST with the t-inflated part SE) as a committed, tested script before any TD/TL grading; ionmc IDD and projected lateral scoring grids on the TOPAS bin edges; the revision-2 TOPAS bundles are committed with revision 2.
3. Tests: the CI rows, full `pytest tests/ionmc` via the logged single-process runner, pre-commit
   on tracked and new files, CI green at the exact head.
4. Codex review on the exact commit and the current base, re-reviewed until passed. Independent
   scientific review of the evidence summary (separate from Codex).
5. Exact-SHA local validation record (lv8, hr8, REF summaries) at the PR head; PR, CI, merge,
   retire.
6. **Integration criterion:**
   - all gating rows pass at the frozen SHA, or are declared evidence gaps with interventions filed;
   - R2 bit-identical;
   - ledger updated (V1-MUST-005, V2-ION: He/C/O EM implemented with evidence @SHA, gaps listed,
     nuclear in V3-009);
   - decision 0043 outcome appended from committed result files;
   - no NIST/ICRU/ENDF/Geant4-data values in Git;
   - no clinical claim.

## G. Minimal acquisition list (via `acquire_reference_data`; roles as stated)

| # | Item | URL | Size | Role |
|---|---|---|---|---|
| 1 | Geant4 v11.4.2 sources: `G4UrbanMscModel.cc`, `G4EmStandardPhysics_option4.cc`, `G4LindhardSorensenIonModel.cc`, `G4LindhardSorensenData.cc`, `G4BetheBlochModel.cc`, `G4IonICRU73Data.cc` | `https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/source/processes/electromagnetic/{standard,…}/src/…` and `…/physics_lists/constructors/electromagnetic/src/G4EmStandardPhysics_option4.cc` (exact msc path TO VERIFY) | ≈ 15–50 kB each | exploratory / lineage and δ_pred (P2) |
| 1b | Geant4 v11.4.2 follow-ups (P2): `G4hMultipleScattering.cc`, `G4EmBuilder.cc` (acquired, to read); `G4BetheBlochModel.cc/.hh`, `G4IonisParamMat.cc`, `G4VMscModel.cc/.hh`, `G4EmTableUtil.cc`, `G4EmParameters.cc`, `G4UrbanMscModel.hh`; `G4EmStandardPhysics_option3.cc` (P2b′ only) | raw.githubusercontent.com/Geant4/geant4/v11.4.2/source/… (paths in `p2-geant4-ion-models.md` §7) | ≈ 15–60 kB each | exploratory (lineage wording, V2 extension) |
| 1c | Geant4 v11.4.2 single-scattering sources (r6, acquired, exploratory / lineage): `G4CoulombScattering.cc` (db354d80…), `G4eCoulombScatteringModel.cc` (fa5ddeb9…), `G4WentzelOKandVIxSection.cc` (a64c94b1…) and `.hh` (3a3f1db0…); `G4EmStandardPhysicsSS.cc` (8aeaf91c…). **Still needed (P10 §7.4):** `G4IonCoulombScatteringModel.cc/.hh`, `G4ionCoulombCrossSection.cc/.hh`, `G4EmParameters.cc` (defaults of NuclearFormfactorType, ScreeningFactor, FactorForAngleLimit), `G4eCoulombScatteringModel.hh`. **Dated note 2026-11-01: acquired and read (P10 §7.4.1):** `G4IonCoulombScatteringModel.cc` ef9b31d4… (8 471 B), `.hh` a1a9a4f5… (5 956 B), `G4IonCoulombCrossSection.cc` 87989b93… (7 432 B), `.hh` cd5d0b5c… (5 536 B), `G4EmParameters.cc` aaa143c1… (43 246 B), `G4eCoulombScatteringModel.hh` 68ce763e… (6 844 B). Nothing is still needed for the V5-SS lineage. Outcome: branch B | `https://gitlab.cern.ch/geant4/geant4/-/raw/v11.4.2/source/processes/electromagnetic/standard/{src,include}/<file>`; `…/electromagnetic/utils/src/G4EmParameters.cc` | ≈ 5–40 kB each | exploratory (V5-SS lineage; Geant4 Software License) |
| 2 | G4EMLOW 8.8 (acquired, sha256 b60cfd63…) | `https://geant4-data.web.cern.ch/datasets/G4EMLOW.8.8.tar.gz` | 350 393 595 B | lineage only — use restricted to Geant4 by the data set's README; not an evaluation source |
| 3 | libdEdx source with data files (acquired, sha256 b941a0c3…, acquisition f062ca75…; GPL-3) | acquired | 962 747 B | exploratory until I3 is frozen; then evaluation (I3: ICRU73_NEW, target 276, Z = 6, 8; `p5b-libdedx-icru73.md`) |
| 4 | Sihver, Schardt, Kanai 1998, Jpn. J. Med. Phys. 18:1 (PDF, J-STAGE free access, sha256 68567feb…) | J-STAGE | 934 294 B | FAILED as a V3 reference (P6): curves normalised horizontally to the calculated peaks; figures only; 18O → V3 evidence gap |
| 5 | (optional) IAEA stopping DB tuples, C in H2O | via `https://www-nds.iaea.org/stopping/` | kB | exploratory: documents the measured-data limitation |

## Evidence archive

| Row(s) | Archive | Head SHA | Seed base | Result |
|---|---|---|---|---|
| TD-EM (main set, 9 cases) | `validation/results/transport/v8r2-topas/reports/TD-EM-case<1..9>.json`; signed run records `…/records/runs/<run_id>/` (index `…/records/index.json`, 36 IDD-EM parts) | runs at `ead583f1`; records and grading at `94bc948e` | TOPAS 20619000 + 100·part + case (parts 0–3); ionmc sample seal `0c68b78a` (REVIEW-9691e1b0) | **PASS 3** (He 150, C 290, O 350), **FAIL 6** (He 50, C 80, O 100: (d) PE; He 250, C 430, O 430: (c) z_peak; C 430, O 430 also (d)); see the dated results note of 2026-11-09 (STOP 3) |
| TL-EM (main set, 9 cases) | `validation/results/transport/v8r2-topas/reports/TL-EM-case<1..9>.json`; signed run records as above (36 LAT-EM parts) | runs at `ead583f1`; records and grading at `94bc948e` | TOPAS 20619000 + 100·part + case (parts 4–7); ionmc sample seal `0c68b78a` | **PASS 1** (C 80), **FAIL 8**: 25 of 36 points invalid by §4.2 fit non-convergence (2-cycle at a bin edge in single batches/parts), 1 point TOPAS σ undefined; every evaluable gated point inside the band (r = 1.17–1.31); D_ext NOT EVALUABLE (He, C, O) |
| TD-H = TD-EM + TD-FULL (held-out, cases 10–12) | `validation/results/transport/v8r2-topas/reports/TD-EM-case<10..12>.json`, `…/TD-FULL-case<10..12>.json`; signed run records `…/records/runs/<run_id>/` (index `…/records/index.json`, 24 IDD parts); host records `v8r2-samples/host-records/` | runs at `55e6a180`; records at `b18fa715`; grading at `1fc287cf` | TOPAS 20619000 + 100·part + case (parts 0–3, 8–11); ionmc held-out seal `ca21ec18` (REVIEW-557971cc), base 20621004 | **PASS 6/6** (TD-EM 3/3, TD-FULL 3/3); see the dated held-out results note of 2026-11-11 |
| TL-H = TL-EM + TL-FULL (held-out, cases 10–12) | `…/reports/TL-EM-case<10..12>.json`, `…/TL-FULL-case<10..12>.json`; signed run records as above (24 LAT parts) | runs at `55e6a180`; records at `b18fa715`; grading at `1fc287cf` | TOPAS as above (parts 4–7, 12–15); ionmc held-out seal `ca21ec18`, base 20621004 | **PASS 6/6** (TL-EM 3/3, TL-FULL 3/3); no point invalid under the §4.2 window-cycle rule; D_ext not defined |
| I9 (lv8, first run) | `validation/generated/RUN-20261111T132917Z-098c9e9a/transport/lv8-e602a583-28/28-i9-combine.txt` (host-run store and task worktree; signed host records RUN-20261111T124822Z-f978c3d3 … RUN-20261111T132114Z-876468a3 (shards 22–27), RUN-20261111T132917Z-098c9e9a (combine); not in Git) | `e602a583` | lv8 20581004 (r 44; consumed) | **FAIL** — variance/A2 = 1.733–1.740 at every point (He 150, C 290, O 350; f 0.25–0.9); attributed to the reference helper's placeholder 1/X_S = 0.0123 cm²/g (amendment I9-REF, item 4k); never re-graded; re-run as a fresh sample at base 20572004 |
| I12 (lv8, first run; abandoned) | shard 29 only (`RUN-20261111T133119Z-94ee45c3`, i12-theta-s01-s0, pass as a shard); no combine, no grade | `e602a583` | lv8 20581004 (r 45; shard 0 consumed) | ABANDONED after the I9 attribution (same reference helper); re-run at base 20572004 |
| I9 (lv8, re-run) | `validation/results/transport/lv8-712f484a-steps22-135-single-process.json` (step `28-i9-combine`); signed host records `validation/results/transport/v8-lvhr/host-records/` (shards 22–27 `RUN-20261111T140753Z-4ac10d13` … `RUN-20261111T144115Z-ce06b608`, combine `RUN-20261111T144917Z-4c7e600d`); step documents in the host-run store | `712f484a` | lv8 re-run base 20572004 (r 44; consumed) | **PASS** — variance/A2 0.9954–1.0012 at every point (He 150, C 290, O 350; f 0.25–0.9), tolerance 2 %; see the dated results note of 2026-11-13 |
| I12 (lv8) | same summary (step `65-i12-combine`); host records steps 29–64 and combine `RUN-20261111T214033Z-180cba48` | `712f484a` | 20572004 (r 45; consumed) | **PASS** — θ_rms deviation from quadrature ≤ 0.22 % (tol 0.5 %), pairwise θ spread 0.36 % (tol 0.5 %), σ(0.9 R) spread 0.44 % (tol 1 %) |
| I14 (lv8) | same summary (step `114-i14-combine`); host records steps 66–113 and combine `RUN-20261112T225438Z-a8fa755d` | `712f484a` | 20572004 (r 46; consumed) | **PASS** — max pairwise \|ΔR80\| 0.0037 mm (tol 0.05 mm); grouped IDD p 0.24–0.96, max \|z\| ≤ 3.36 < 4.52 for all five configurations vs s1-f002 |
| I15 (lv8) | same summary (step `135-i15-combine`); host records steps 115–134 and combine `RUN-20261113T015413Z-8c57c225` | `712f484a` | 20572004 (r 47; consumed) | **PASS** — pairwise θ_rms spread 0.29 % (tol 0.5 %), σ spread 0.97 % / 0.61 % at 0.5 R / 0.9 R (tol 1 %) |
| T13-ion (hr8) | `validation/results/transport/hr8-712f484a-single-process.json` (step `01-t13-ion-chunks-cuda-c290`); host record `RUN-20261113T015455Z-790a1181` | `712f484a` | hr8 20591004 (consumed) | **PASS** — chunks 2¹⁰ vs 2¹⁸ bit-identical (edep, tallies, counters) in float32 and float64 |
| T12-ion HR (hr8) | same summary (step `12-t12-ion-hr-combine`); host records steps 02–11 and combine `RUN-20261113T033810Z-176a1c8f` | `712f484a` | 20591004 (r 54; consumed) | **PASS** — all four pairs (C 290 cpu64/cuda32, cuda32/cuda64, He 150 cpu64/cuda32, O 350 cpu64/cuda32) pass the frozen T12 rules; array p ≥ 0.08, scalar \|z\| ≤ 2.05 (bound 3.5), total deposit identical to ≤ 8.7e-10 |

## Amendments

(none before revision 1)

**Revision 1 (2026-11-01; recorded before any V3-008 transport or qualification result; orchestrator decisions where stated):**
1. Lateral-gap intervention: the clause "or an accessible single-scattering MC reference" is amended (STOP 4; decision 0043 §3) so that a reference sharing the screened-Coulomb cross-section lineage does not satisfy it; the intervention is filed at the freeze whatever V5-SS shows and names both gaps (orchestrator decisions (a), (b); section "Intervention at the freeze").
2. The f_dM extrapolation is a declared limitation; no further coverage is claimed (orchestrator decision (c)).
3. He particle naming: `alpha` (= GenericIon(2,4)), BetheBloch lineage (P4; orchestrator decision; §B).
4. V2 = Case C by the pre-declared branch (P2 F5, P4); no amendment of V2 itself.
5. Pre-registration status updated to f694b747 (`p9-freeze-record.md` §6): the I0 sub-items and the I1 mapping-callable parity exercised by the 683480b tests are carried as SEEN.

**Amendment (2026-11-01) — FC-2 disposition** (planner record `validation/plans/v3-008/p2b-follow-ups.md`; recorded before any V3-008 transport or qualification result; no row definition, tolerance or seed changed):
6. FC-2 is disposed by the executed record `p2b-follow-ups.md` (option (a) of `p9-freeze-record.md` §0: bounded source reads and a read of the archived smoke logs, no new run). Pointer notes in rows P2, P4, P5 and in the freeze audit.
   - `G4BetheBlochModel.cc/.hh` pinned; the ICRU 90 branch settled for G4_WATER + option 4. The `p3b-ls-default.md` §5 item is resolved: the GenericIon base table is the ICRU 90 proton array × q², not the L323 corrections. A possible double count of the z = 1 part of ΔL_LS is a lineage item.
   - G4_WATER: I = 78 eV (source constant, FACT); X0 = 36.083 g/cm² derived from the source composition; ratio to ionmc's cross-check value 1.0000000.
   - P4 record corrected: the material table was never archived. "Archived" now means the source definition plus the derived X0; the optional `Ma/Verbosity` smoke is not required.
   - The final msc parameters for GenericIon/alpha are FACT. `dtrl` and the λ charge scaling are OPEN at the inline-header level and declared lineage residuals (V2 Case C attribution only).
   - icru90 z6/z8 read evidence absent: a declared limitation of the lineage record; no frozen row depends on it.
7. Orchestrator decision (d), I3 attribution: if ICRU73_NEW is the pre-revision water table, the I3 result stands as evaluated with the attribution "reference-table version unresolved (P5b)" and is never re-graded (row I3 note).
8. Orchestrator decision (e), V3 bundling: the V3 measured-range gap is bundled into the intervention filed at the freeze (one request, two unavailable inputs; section "Intervention at the freeze"; row V3 note).

Bound at revision 1 but recorded earlier as dated notes in the rows: the LS option with default "off" (decision 0043 §1 amendment, round V3-008-5); V5-SS branch B (r6 continuation), its TOST rule and the 6e4-primary design (round V3-008-5).

**Revision 2 (2026-11-07; operator directive of 2026-11-07, IR-20261102-010435-B59B2D resolved 2026-11-07 20:37 UTC; planner record `validation/plans/v3-008/r2-topas-target.md`; recorded before any V3-008 transport or qualification result; round R2-TOPAS):**
9. Validation target: agreement with the specified TOPAS configuration (OpenTOPAS 4.3 / Geant4 11.4.2, seven modules) for depth-dose distributions and lateral profiles; EM-only TOPAS as intermediate; experimental validation deferred (directive, verbatim file sha256 ca079e9f0af990cdff9d8efde60d38ac55e2d0fc0fb9ef22355a81f4a2ab0b26).
10. New rows TD-EM, TD-FULL, TL-EM, TL-FULL, TD-H/TL-H (gating as stated) and TF-DIAG (report-only); V1, V1-held, V2, V4 superseded (never run); V3 kept as a revision-1 EVIDENCE GAP and superseded as a target; V5-SS, I3 and every A-row unchanged.
11. STOP 4 resolved by the directive (dated resolution in §D); STOP 3 extended to the TD/TL rows; STOP 5 added (FULL fails while EM passes; module rejection).
12. f_dM extrapolation: unchanged model (model-design decision); I11-ext report-only; adequacy judged by TL-EM (b) and TL-FULL (b).
13. TOPAS case numbers re-assigned within the unchanged formula (§C); 205 bundles committed (192 evidence parts, 13 execution-only smokes); ionmc r_index 30–41 declared.
14. Bookkeeping (Codex REVIEW-ad85dff1): P1 row FC-1 record rebound to the current hash a1a24cd3… (e3cd3258… kept as superseded history); P9 §7 repository-path count corrected to 60 (dated note in `p9-freeze-record.md`).
15. Revision 2 is the frozen revision once committed by a reviewed commit; revision 1's freeze is superseded before any result (`p9-freeze-record.md`, "Revision 2 additions").

**Dated entry (2026-11-09; plan section "Order" step 2 — sizing note; planner).** The 13 step-1 smokes ran at 0330a82f, after the seal commit 0c68b78a and its review REVIEW-9691e1b0. All ran validly (exit 0, no killed or stuck track). UrbanMsc for alpha/GenericIon; the seven modules accepted, with QMDModel and radioactive decay shown in the printout.

Measured t_s and rule outcomes: `r2-topas-target.md` §8 "Dated sizing note (2026-11-09)".
- EM configurations and He IDD-FULL hold.
- Rule (i): C IDD-FULL N_part 3333 (13 332), O IDD-FULL 2529 (10 116), He LAT-FULL 4789 (19 156).
- Bundle regeneration (`NumberOfHistoriesInRun`) by reviewed commit before these runs. Seeds unchanged.
- **C and O LAT-FULL: STOP under rule (iii)** (7198 s / 9816 s per 5000-history part; even six parts reach < 1e4). Recorded and preserved.
- The step-length smoke shows that the scorer bins limit the primary step in IDD runs (97.5 %, 0.1 mm). The TL band is unchanged (§9 dated resolution).
- Projected REF without C/O LAT-FULL: ≈ 15.6 h (≈ 12.5 h with cost ∝ R), He LAT-FULL at 0.05 mm included; cap 16 h. Under R2-SIZE-1, ≈ 13.3 h plus LAT-FULL-A.

A **proposed** amendment, R2-SIZE-1 (LAT-FULL with `CutForElectron = 1 km`, the seven modules unchanged, three timing smokes 20600424–20600426 and a LAT step smoke 20600427, then re-sizing under the unchanged rule), is recorded in §8. It is **not in force** until independently reviewed and committed before any LAT-FULL run.

**Orchestrator decision (2026-11-09, same day, before any evidence run).** (1) R2-SIZE-1 option (A) is ADOPTED for the TL-FULL row of all three ions, helium included, so that the row and its energy-trend criterion (5.3 b) use one configuration (seven modules unchanged; `CutForElectron = 1 km` in the LAT-FULL bundles only; IDD-FULL keeps 0.05 mm); the ledger wording for the TL-FULL reference is "seven-module TOPAS configuration, delta electrons deposited locally". It enters into force only when the commit holding this entry and the regenerated bundles has passed its independent review; the C/O LAT-FULL STOP at 0.05 mm stays recorded as history. (2) The 16 h REF cap is a planning value of this plan (section E), not an operator directive; the decision whether to raise it is deferred until the three `ecut1km` timing smokes have measured t_s (threshold ≈ 85 ms at the upper ends); the V5-SS 4e4 fallback applies first if the cap is exceeded. (3) Rule (ii) (six parts) is not used; if it were, the σ_T inflation for K = 6 would have to be amended first. (4) Follow-up implementation: `gen_r2_bundles.py` per-(ion, configuration) history counts (C IDD-FULL 3333, O IDD-FULL 2529, He LAT-FULL 4789 at 0.05 mm are superseded for LAT-FULL by the post-smoke re-sizing), the 1 km cut for lat-full, smoke seeds 20600424–20600427, new bundle-set digest, `test_r2_bundles.py` updated; smoke evidence record `validation/results/transport/v8r2-topas/smokes.json`. **(2026-11-09, later the same day)** The bundles were regenerated accordingly (209 bundles, digest e51d3b3d…; `r2-topas-target.md` section 7 dated note and appendix A 2026-11-09 listing); R2-SIZE-1 enters into force when the commit holding this regeneration passes its independent review.

**Dated entry (2026-11-09, after the R2-SIZE-1 smokes; re-sizing note; planner).** The four smokes 20600424–20600427 ran at 39cf5478. The seal 0c68b78a is an ancestor. Exit 0, no killed or stuck track. Their run records are **unsigned schema-2 records**: the operator has not yet restarted the reference service after V3-INFRA-004. They are execution and sizing records only, not evidence. Details: `r2-topas-target.md` §8 "Dated re-sizing note (2026-11-09, after the R2-SIZE-1 smokes)".
- **LAT-FULL** (`CutForElectron` 1 km):
  - t_s He 45.2 / C 114.7 / O 114.3 ms; × 5000 = 226 / 574 / 572 s, all ≤ 1200 s.
  - The part plan holds: 5000 × 4 for every case. No bundle change.
- **Final parts (K = 4):** IDD-EM and LAT-EM 25 000; IDD-FULL He 5000 / C 3333 / O 2529; LAT-FULL 5000.
- **Budget.** Projected REF 16.25 h (13.12 h with cost ∝ R) against the 16 h cap.
  - The pre-registered V5-SS 4e4 fallback applies: −0.98 h, giving 15.27 h. No STOP.
  - Prospective option for the orchestrator: raise the cap (a planning value of §E, not an operator directive) to 17 h and keep V5-SS at 6e4. The planner recommends this. Either way, the choice is recorded before the first evidence run.
- **LAT step regime measured:** plateau median 0.63 mm, inside the 0.1–1 mm span. The TL band is unchanged (`r2-topas-target.md` §9).
- `Decay` is listed in the ListProcesses printout, so the executed module set, `g4decay` included, is verifiable.
- **Status of R2-SIZE-1 (a fact for the orchestrator to confirm).** R2-SIZE-1 is in force if and only if the commit chain carrying it, 3901421c … 39cf5478, with Codex REVIEW-841e82ae (pass), counts as its independent review.
  - The review id and its result are as relayed by the orchestrator. The planner has not seen that review record.
  - The orchestrator also states whether that review finished before the smokes started (≈ 11:39Z). The smokes are execution-only either way.

**Orchestrator decisions (2026-11-09, after the re-sizing note; before any evidence run).** (1) REF cap: raised from 16 h to **17 h** (a planning value of this plan, section E; not an operator directive); V5-SS keeps its 6e4 design and its 4e4 fallback is NOT invoked (projected REF total 16.25 h under the conservative scaling, 13.1 h with cost ∝ range; margin 0.75 h). Reason: V5-SS is the only row resolving b(z) at the 3 % level and the fallback would widen its statistical uncertainty by ≈ 1.22×; the overrun is 1.6 % of an upper projection. (2) R2-SIZE-1 is **in force**: the commit chain 3901421c → 39cf5478 (regenerated bundles with the 1 km electron cut for LAT-FULL and ListProcesses, frozen manifest 6d1ae3a0…) passed its independent review, Codex REVIEW-841e82aebece44679c0e9138b8aa361c (head 39cf5478, verdict pass, finished 2026-11-09T11:36:11Z); the four R2-SIZE-1 smokes started after that (first run 11:39:19Z) and are execution-only records in any case. (3) LAT-FULL stays at 5000 histories per part (rule holds: ≤ 574 s per part); the final part table of the re-sizing note is the evidence plan. (4) Evidence runs start only after the operator restart of the reference service (IR-20261109-112802-D6BAE3) and the publication of the reference-runner key (signed, timestamped run records; unsigned records are not evidence).

**Dated amendment (2026-11-09, Codex REVIEW-8429c50e; orchestrator) — run validity of the 17 smokes and the provisional status of the sizing decisions.** Plan rule R2-TOPAS-4 applies to every revision-2 TOPAS run, smoke or evidence, and (since V3-INFRA-004) requires the runner's signed, timestamped record. All 17 smokes (13 at 0330a82f, 4 at 39cf5478) were executed BEFORE the reference service was restarted with V3-INFRA-004, so their records are unsigned (schema 2) and their run validity under R2-TOPAS-4 is **NOT ESTABLISHED**. They are preserved as execution records (`validation/results/transport/v8r2-topas/smokes.json`, REF ids, wall times, log checks) and are never evidence. Consequently every statement derived from them is **PROVISIONAL**: the sizing note, the LAT-FULL re-sizing (N_part 5000 at the 1 km cut), the step-length findings (sections 8 and 9), the 17 h cap decision and the in-force status of R2-SIZE-1 (its bundle content is reviewed and frozen; its sizing confirmation is pending). Reconciliation procedure (pre-registered now): after the operator restarts the reference service and the `reference-runner` key is published (IR-20261109-112802-D6BAE3, key-publication task), the sizing-relevant smokes are RE-RUN once with signed records at the then-current committed SHA: the six EM timing smokes (20600411, 20600412, 20600415, 20600416, 20600419, 20600420), the three IDD-FULL timing smokes (20600413, 20600417, 20600421), the three R2-SIZE-1 LAT-FULL smokes (20600424–20600426) and the two step-length smokes (20600423, 20600427); the three 0.05 mm LAT-FULL timing smokes (20600414, 20600418, 20600422) are not re-run (that configuration was replaced by R2-SIZE-1 and their STOP record stands). The re-runs must pass `r2_io.run_checks` and the R2-TOPAS-4 chronology against the seal review; their `t_s` re-confirm the part plan under the unchanged section-6 rule: if every re-run part time stays ≤ 1200 s at N_part of the final part table, the provisional decisions become final by a dated confirmation note (same numbers or better); if any part exceeds 1200 s, the re-sizing note is redone from the signed values before any evidence run; if any re-run fails run validity, STOP (recorded). No evidence run starts before this confirmation. The step-length and ListProcesses diagnostics are likewise re-read from the signed re-runs (same seeds; identical results expected for the single-threaded engine, which is recorded as a reproducibility check, not a claim).

**Dated STOP record and amendment R2-SMOKE-LP (2026-11-09, later; orchestrator).** The pre-registered signed re-runs of the 14 smokes were made at fce29f49 after the publication of the `reference-runner` key (V3-INFRA-005, develop 1e93fe1e). All 14 records verify (signature, identity, request and frozen-manifest binding, output digest) and all 14 pass the R2-TOPAS-4 chronology (seal review finished 2026-11-09T09:15:57.707596Z, every `started_at` later). **Run validity (`r2_io.run_checks`): 4 VALID, 10 RUN-VALIDITY FAIL → STOP, recorded** (`validation/results/transport/v8r2-topas/smoke-reruns.json`, per-run records under `smoke-runs/<run_id>/`). Cause: the ten failing smokes (the original timing smokes 20600411–20600413, 20600415–20600417, 20600419–20600421 and the IDD-EM step smoke 20600423) use the thirteen original bundles that were deliberately kept byte-identical without `b:Ph/ListProcesses = "True"`, while the frozen-configuration rule now requires that line in every bundle and, for FULL bundles, the `Decay` process-list marker; stderr, histories, UrbanMsc and the hadronic markers of those runs are clean. The four VALID runs (the three R2-SIZE-1 LAT-FULL smokes 20600424–20600426 and the LAT-EM step smoke 20600427) stand. Resolution (prospective, effective after independent review of the commit carrying it): ten new execution-only smoke bundles, the ListProcesses variants of the ten failing ones (`...-lp-timing-smoke`, `carbon-water-290mevu-r2-idd-em-lp-steps-verbose-smoke`; seeds 20600428–20600437; identical to their originals except the process-list line and the seed; old bundles untouched; frozen manifest regenerated: 219 bundles), to be run once with signed records; the sizing confirmation then rests on the 10 lp runs plus the 4 VALID re-runs (all 14 must be VALID under `check-run`); otherwise STOP again. The measured re-run timings (all within 4 % of the unsigned runs; part times at the final part table: C IDD-FULL 1190 s, O IDD-FULL 1186 s, all others ≤ 575 s) are recorded but remain provisional until the lp runs confirm them.

**Dated correction (2026-11-09, Codex REVIEW-8286be8c; orchestrator).** The confirmation set of R2-SMOKE-LP cannot include the four VALID re-runs made at fce29f49: the frozen bundle manifest changed with the ten lp bundles (digest 6d1ae3a0… → 76d78566…) and the evaluator rejects any run whose manifest at `code_sha` differs from HEAD's (the rule is intentional and is not relaxed). Those four runs remain preserved records (VALID at their own manifest). The confirmation set is therefore the **14 runs made at the manifest-current head**: the ten lp smokes (20600428–20600437) and re-runs of the four R2-SIZE-1 / LAT-EM step smokes (20600424–20600427), all at commit 4cdc44ce or a descendant with the same manifest; all 14 must be VALID under `check-run` run at the committed head of this repository (the earlier `smoke-reruns.json` evaluator SHA f5512c4e was a throwaway clone used because records must be committed before `check-run` can read them — the regenerated summary records the real head). The `check-run` bundle-name grammar is extended to the lp names in the same commit. Any further manifest change before the evidence runs would again void the confirmation and require re-runs.

**Dated confirmation note (2026-11-09, late; orchestrator) — R2-SMOKE-LP confirmation set VALID; IDD-FULL re-sized; smoke binding per bundle.** The 14 runs of the confirmation set (ten lp smokes, seeds 20600428–20600437, nine at 4cdc44ce and one at 1141df44; re-runs of 20600424–20600427 at 1141df44) are **all VALID** under `check-run` at the committed head 111841e2 (`validation/results/transport/v8r2-topas/smoke-confirmation.json`, `confirmation_set_complete: true`; every record verifies, passes the R2-TOPAS-4 chronology against the seal review and all run-validity checks incl. the executed-module markers with `Decay` in the process list). Signed timings and part times at the final part table (this note re-sizes IDD-FULL, see below): he idd-em: t_s 7.0 ms, part 175 s at N = 25000 (REF-022775f8c6a8530c48aa-a50056d1); he lat-em: t_s 2.5 ms, part 62 s at N = 25000 (REF-2c6d31fa3f3ae83bce78-ea24d511); he idd-full: t_s 70.4 ms, part 352 s at N = 5000 (REF-1bee6c0ab49920cf6ff3-685816fc); carbon idd-em: t_s 6.4 ms, part 159 s at N = 25000 (REF-964a865042b11a34b6bc-688112f3); carbon lat-em: t_s 2.0 ms, part 49 s at N = 25000 (REF-706eae23ccef7407d213-b7cff919); carbon idd-full: t_s 364.8 ms, part 1167 s at N = 3200 (REF-494fafda858a0a1f2c2d-5e74e4ef); oxygen idd-em: t_s 5.1 ms, part 127 s at N = 25000 (REF-fda97285d2e8ab1adcb5-9699e0a1); oxygen lat-em: t_s 1.6 ms, part 40 s at N = 25000 (REF-ae7dbb3bd0e51a23a709-e348d838); oxygen idd-full: t_s 475.9 ms, part 1190 s at N = 2500 (REF-77d4d10efe5c035b544e-caf6cb7f); he lat-full: t_s 44.4 ms, part 222 s at N = 5000 (REF-100a453268dacee7ea3f-80664e3e); carbon lat-full: t_s 114.5 ms, part 573 s at N = 5000 (REF-c26d98d1bb1665012906-ef220531); oxygen lat-full: t_s 114.3 ms, part 572 s at N = 5000 (REF-c7ab0adc969cd4034512-f5df42e7). Two parts exceeded the 1200 s planning cap at the previous counts (C IDD-FULL 364.8 ms × 3333 = 1216 s; O IDD-FULL 475.9 ms × 2529 = 1204 s; run-to-run timing jitter ≈ 2 %), so section-6 rule (i) is applied: **C IDD-FULL 3200 per part (12 800 per case), O IDD-FULL 2500 per part (10 000 per case, the statistical minimum)**; all other counts unchanged (EM 25 000; He IDD-FULL 5000; LAT-FULL 5000 at the 1 km electron cut). The 32 IDD-FULL bundles are regenerated accordingly (frozen manifest regenerated; appendix A listing). The sizing, step-length (IDD 97.5 % bin-limited at 0.1 mm; LAT 83 % plane-limited, median 0.45 mm) and cap statements of the sizing and re-sizing notes are hereby **confirmed** with signed records; the provisional status is lifted. Amendment to the evaluator's manifest rule (prospective, reviewed with this commit): an EVIDENCE run is bound to the complete frozen manifest (manifest bytes identical at its `code_sha` and at HEAD, unchanged rule); a SMOKE run is bound to its own bundle's manifest entry (identical at `code_sha` and HEAD and equal to the signed request), so that re-sizing evidence bundles does not void execution-only smoke records; the whole-manifest digests are reported. The evidence runs may start once the commit carrying this note passes its independent review (R2-TOPAS-4 and the plan order then apply to every part).

### Dated results note and STOP 3 (2026-11-09, orchestrator; plan Order step 3: TD-EM and TL-EM graded)

**Evidence.** The 72 EM evidence parts (IDD-EM parts 0–3, LAT-EM parts 4–7, cases 1–9; 25 000 primaries each) ran
sequentially on 2026-11-09 18:10–19:47 UTC at `code_sha` ead583f1 (the clean committed head after REVIEW-5baeb516;
the seal commit 0c68b78a is an ancestor; REVIEW-9691e1b0 finished 09:15:57 UTC). Their signed schema-3 run records
(signer `reference-runner`, key_id 303cce4c129d387c) and artifact listings are committed under
`validation/results/transport/v8r2-topas/records/` (commit 94bc948e); every listed file was verified against its
sha256 and the signed `output_files_sha256` reproduces from the `work/` entries. `r2_evaluate.py check-run`:
72/72 VALID (binding, chronology, frozen manifest, module markers). The rows were then graded, once, with
`r2_evaluate.py grade` at evaluator HEAD 94bc948e (sealed ionmc samples: the 86 chunks of `sample-seal.json`);
the reports (statistics, record ids and sha256 only) are
`validation/results/transport/v8r2-topas/reports/<row>-case<N>.json`, bound by sha256 in the P9 freeze record.
**They are final: nothing in them is re-graded** (plan D, FAIL).

**TD-EM (gating; all five criteria per case).** 3 PASS, 6 FAIL.

| Case | Status | (a) R80 ionmc / TOPAS [mm] (stat ≤ tol) | (b) W ionmc / TOPAS [mm] (report-only) | (c) z_peak ionmc / TOPAS [mm] (stat ≤ 0.5) | (d) PE ionmc / TOPAS (stat ≤ tol) | (e) plateau groups (failing) |
|---|---|---|---|---|---|---|
| He 50 | **FAIL** | 22.49 / 22.55 (0.07 ≤ 0.08) PASS | 0.163 / 0.172 (-5.4 %) | 22.36 / 22.42 (0.06) PASS | 9.91 / 9.19 (0.74 ≤ 0.28) FAIL | 17 (0) PASS |
| He 150 | **PASS** | 159.22 / 159.57 (0.38 ≤ 0.56) PASS | 1.114 / 1.107 (+0.7 %) | 158.39 / 158.73 (0.36) PASS | 9.39 / 9.38 (0.03 ≤ 0.28) PASS | 25 (0) PASS |
| He 250 | **FAIL** | 382.54 / 383.36 (0.86 ≤ 1.00) PASS | 2.561 / 2.572 (-0.4 %) | 380.76 / 381.57 (0.91) FAIL | 9.13 / 9.00 (0.15 ≤ 0.27) PASS | 61 (0) PASS |
| C 80 | **FAIL** | 17.50 / 17.52 (0.02 ≤ 0.06) PASS | 0.072 / 0.087 (-17.4 %) | 17.43 / 17.43 (0.00) PASS | 11.47 / 10.35 (1.13 ≤ 0.31) FAIL | 13 (0) PASS |
| C 290 | **PASS** | 163.21 / 163.69 (0.49 ≤ 0.57) PASS | 0.621 / 0.644 (-3.6 %) | 162.68 / 163.16 (0.48) PASS | 11.11 / 10.97 (0.18 ≤ 0.33) PASS | 25 (0) PASS |
| C 430 | **FAIL** | 309.31 / 310.13 (0.84 ≤ 1.00) PASS | 1.130 / 1.189 (-5.0 %) | 308.50 / 309.24 (0.76) FAIL | 10.60 / 10.30 (0.34 ≤ 0.31) FAIL | 49 (0) PASS |
| O 100 | **FAIL** | 19.57 / 19.58 (0.00 ≤ 0.07) PASS | 0.066 / 0.079 (-16.8 %) | 19.50 / 19.48 (0.02) PASS | 11.10 / 9.96 (1.15 ≤ 0.30) FAIL | 15 (0) PASS |
| O 350 | **PASS** | 166.61 / 167.10 (0.49 ≤ 0.58) PASS | 0.538 / 0.562 (-4.3 %) | 166.14 / 166.62 (0.48) PASS | 11.37 / 11.20 (0.20 ≤ 0.34) PASS | 26 (0) PASS |
| O 430 | **FAIL** | 231.83 / 232.52 (0.71 ≤ 0.81) PASS | 0.729 / 0.768 (-5.1 %) | 231.23 / 231.89 (0.67) FAIL | 11.04 / 10.66 (0.41 ≤ 0.32) FAIL | 37 (0) PASS |

Findings: (a) R80 passes in all nine cases, but ionmc is systematically short: −0.29 % (He 50), −0.22 % (He 150),
−0.21 % (He 250), −0.09 % (C 80), −0.29 % (C 290), −0.27 % (C 430), −0.02 % (O 100), −0.29 % (O 350), −0.30 %
(O 430); at the standard points and upper ends this exceeds the §5.1 stopping-construction budget (≤ 0.18 % C/O, ≤ 0.10 % He) — a finding for the
I3b/I3c attribution. (c) z_peak fails at the three upper ends (ionmc 0.66–0.81 mm proximal at ranges 232–383 mm; the
absolute 0.5 mm tolerance is tighter there than the relative (a) tolerance, by design of revision 1 — not changed).
(d) PE fails at the three low-energy cases (ionmc/TOPAS − 1 = +7.8 % He 50, +10.8 % C 80, +11.4 % O 100) and at
C 430 (+2.8 %) and O 430 (+3.5 %); the report-only distal width W is 5–17 % narrower in ionmc at the same cases
(−5.4 %, −17.4 %, −16.8 %; −5.0 %, −5.1 %), i.e. less energy-loss straggling in ionmc than in TOPAS at low energy per
nucleon (and at the C/O ends), which raises the peak. (e) every plateau group passes in every case (local stopping
agrees within 2 % everywhere). No model is changed by this note.

**TL-EM (gating).** 1 PASS (C 80), 8 FAIL.

| Case | Status | f = 0.25 | f = 0.5 | f = 0.75 | f = 0.9 | (b) D_ext |
|---|---|---|---|---|---|---|
| He 50 | **FAIL** | σ_i/σ_T 1.381, report-only: REPORT-ONLY | σ_i/σ_T 1.259, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.222, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.206, gated: PASS | n/a |
| He 150 | **FAIL** | σ_i/σ_T 1.185, gated: PASS | σ_i/σ_T 1.176, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.178, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.181, gated: FAIL (invalid: a batch/part fit failed) | n/a |
| He 250 | **FAIL** | σ_i/σ_T 1.174, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.178, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.184, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.185, gated: FAIL (invalid: a batch/part fit failed) | FAIL (invalid: a point of D_ext is invalid) |
| C 80 | **PASS** | σ_i/σ_T 1.495, report-only: REPORT-ONLY (outside band) | σ_i/σ_T 1.369, report-only: REPORT-ONLY | σ_i/σ_T 1.301, gated: PASS | σ_i/σ_T 1.280, gated: PASS | n/a |
| C 290 | **FAIL** | σ_i/σ_T 1.275, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.262, gated: PASS | σ_i/σ_T 1.259, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.259, gated: FAIL (invalid: a batch/part fit failed) | n/a |
| C 430 | **FAIL** | σ_i/σ_T 1.265, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.270, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.274, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.277, gated: FAIL (invalid: a batch/part fit failed) | FAIL (invalid: a point of D_ext is invalid) |
| O 100 | **FAIL** | σ_i/σ_T n/a, gated: FAIL (invalid: fewer than 3 bins (or no signal) inside +-2 sigma) | σ_i/σ_T 1.354, report-only: REPORT-ONLY (invalid: a batch/part fit failed) | σ_i/σ_T 1.309, gated: PASS | σ_i/σ_T 1.286, gated: FAIL (invalid: a batch/part fit failed) | n/a |
| O 350 | **FAIL** | σ_i/σ_T 1.290, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.283, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.283, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.278, gated: FAIL (invalid: a batch/part fit failed) | n/a |
| O 430 | **FAIL** | σ_i/σ_T 1.282, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.284, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.287, gated: FAIL (invalid: a batch/part fit failed) | σ_i/σ_T 1.285, gated: FAIL (invalid: a batch/part fit failed) | FAIL (invalid: a point of D_ext is invalid) |

Findings and attribution (allowed under STOP 3; no re-grading): where a point is evaluable, r = σ_core,ionmc/
σ_core,TOPAS-EM = 1.17–1.31 at the gated points (every one inside the band [−5 %, +40 %], as expected for the Urban
deficit) and 1.35–1.50 at the shallow report-only points. 25 of the 36 points are invalid with "a batch/part fit
failed" and 1 with "fewer than 3 bins inside ±2σ" (TOPAS σ undefined at O 100, f = 0.25). A read-only diagnostic
(orchestrator scratch `diag_tl_fits.py`, `diag_trace.py`, reusing the committed estimator) shows the cause: the §4.2
fixed-point iteration of the window-truncated ML fit enters a **2-cycle when ±2σ crosses a bin edge** — e.g. He 50,
slab 11–12 mm, ionmc batch 20, x: the window alternates between 32 and 34 bins and σ between 0.082377 and
0.082592 mm (0.26 % apart, both ≈ the batch value), never meeting the 1e-4 relative criterion within 50 iterations.
It affects 2–8 % of single-batch/part fits (He 50: 3/80 and 2/80 at the f = 0.5/0.75 slabs; C 290: 6/80 and 5/80 at
f = 0.75/0.9 and one TOPAS part at f = 0.25); the pooled fits converge in 2–3 iterations. By the pre-registered §4.2
rule (non-convergence = invalid point, counted as FAIL for a gated point) **the verdicts stand as graded**; the
failure is attributed to the estimator's iteration, not to a lateral-width disagreement. Criterion (b) D_ext is NOT
EVALUABLE for He, C and O (a point of each pair invalid).

**STOP 3 (plan D).** TD-EM and TL-EM fail → halt; amendment and independent review before continuing. The evidence
loop was halted after the part in flight: the IDD-FULL parts already run (He 50 parts 8–11, He 150 parts 8–11,
He 250 parts 8–9; signed records in the orchestrator ledger, to be committed with the FULL set) are valid evidence
runs whose outputs have not been read; the remaining FULL and held-out parts wait for the amendment. V5-SS is
unaffected. The amendment (R2-STOP3) is written prospectively in a separate reviewed commit; it may change the
estimator's convergence rule only for rows not yet graded, must preserve these grades verbatim, and decides between
the two pre-registered branches of STOP 3 (no model change: the FAILs stand and go to the orchestrator/operator for
integration; model change after attribution: all main TD/TL rows are consumed and the held-out set becomes the sole
gating evidence at a re-run ionmc base, unchanged tolerances).

### Amendment R2-STOP3 (2026-11-09, orchestrator; prospective; in force once the commit carrying it and its implementation has passed independent review; nothing graded is re-graded)

**0. Scope.** Resolves STOP 3 of the dated results note above. Every grade of that note stands verbatim; no
tolerance is changed; no parameter is fitted to TOPAS; no TOPAS configuration, bundle or seed changes.

**1. Attribution (STOP 3 allows I3b/I3c-type table work and the P2 source reading; estimates are labelled).**
- **(a) TD-EM low-end (d) FAILs — dominant cause: ionmc's end-of-range point deposit.** Decision 0043 §4: at
  E ≤ E_cut = 2 MeV/u × A the whole residual energy is deposited at the current point
  (`src/ionmc/transport/kernels.py`, `deposit_point`, counter `t_cutoff`); the residual CSDA range at 2 MeV/u from
  ionmc's own table is 0.092 mm (He), 0.043 mm (C), 0.038 mm (O), against σ_R ≈ 0.05–0.13 mm at the low ends and a
  peak group of w_p = 5 h_z = 0.1 mm. TOPAS tracks ions continuously to rest (no cutoff deposit). The §5.1 (d)
  defence ("E_cut and the peak-group width are identical in both analyses", r2-topas-target.md §5.1) was therefore
  WRONG: the estimator is identical, the end-of-range physics is not; §5.1 (b) had identified the R20 effect only.
  Scratch forward model (orchestrator and planner, SEEN after grading, labelled estimate, not evidence: single-ion
  CSDA depth-dose from ionmc's committed stopping table, Gaussian range straggling from the first-order Bohr
  integral σ_R² = ∫(dΩ²/dx)/S³ dE, bin-averaged on the case grid, reduced with the committed `r2_statistics`
  depth estimators): with the point deposit it reproduces the sealed ionmc W within 1.3 % at all nine cases and the
  low-end PE within 0.5 % (so ionmc's straggling is the declared Bohr model, correctly implemented); replacing the
  point by a continuous residual deposit changes PE by −10.7 / −9.9 / −5.4 % (He 50 / C 80 / O 100), W by
  +7 / +11 / +14 % and R80 by +0.040 / +0.015 / +0.009 mm, against the observed ionmc − TOPAS differences PE
  +7.8 / +10.8 / +11.4 %, W −5.4 / −17.4 / −16.8 %, R80 −0.065 / −0.016 / −0.004 mm. **Residual (INFERRED):** TOPAS
  σ_R exceeds ionmc's Bohr value by ≈ 0–3 % (He), ≈ 4–9 % (C) and ≈ 4–5 % in W but up to ≈ 15 % in PE (O 100);
  ionmc uses the pure Bohr variance with the bare charge (code, `src/ionmc/physics/em.py`). **D-4 preliminary
  (2026-11-11; Geant4 v11.4.2 source reading by the physics researcher; the files were fetched outside the
  acquisition tool and must be re-acquired through `acquire_reference_data` before this reading is committed as a
  P2-type record):** opt4 gives alpha and GenericIon `G4IonFluctuations`, which hands over to
  `G4UniversalFluctuation` above ≈ 10·z·M/m_p MeV (≈ 20 / 60 / 79 MeV/u for He / C / O); with the 1 km electron cut
  the step's maximum transfer equals T_max, and the Urban ("Glandz") sampler of `G4UniversalFluctuation` is used
  unless the step's mean loss exceeds 10 T_max; its variance relative to Bohr rises with energy per nucleon (≈ 0.95
  at 80 MeV/u, ≈ 1.19 at 430 MeV/u, estimate) and depends on the step length; integrated at the 0.1 mm step regime
  of the IDD runs this gives σ_R(TOPAS)/σ_R(Bohr) ≈ 0.99–1.02 (He 150–250) and ≈ 1.03–1.06 (C/O 290–430) — the
  observed sign, ordering and size. The excess is therefore a step-dependent property of the reference's
  fluctuation sampler at high energy per nucleon, not an ion-charge effect and not physics missing from ionmc;
  at the held-out energies (He 100, C 150, O 200) it is ≈ 0–3 % (INFERRED). The charge-state (Yang) and
  relativistic terms of `G4IonFluctuations` act only below the handover and change σ_R by < 0.1 %. Upper-end (d) FAILs (C 430 +2.8 %, O 430 +3.5 %): E_cut contributes ≤ 0.3 %
  (estimate); the rest is straggling-type and partly unexplained.
- **(b) R80 systematically short by 0.2–0.3 % at the standard points and upper ends — on the reference side relative
  to its declared lineage.** ionmc's R80 equals its own table CSDA range within 0.03 % at R ≥ 150 mm (transport
  self-consistent). TOPAS-EM R80 exceeds ionmc's CSDA range by +0.19 % (He 150, He 250) and +0.26 … +0.29 %
  (C/O 290–430). Decision 0038 ("Later validation outcome") predicts the opposite sign for C/O (ionmc longer than
  integrated ICRU 90 by ≤ 0.175 %) and ≤ 0.07 % for He, so the "shared stopping lineage" statement of §B and
  r2-topas-target.md §2 holds only to ≈ 0.2–0.45 %; the §5.1 (a) budget is not violated by ionmc. Mechanism TO
  VERIFY (D-3: ICRU 90 arrays of the Geant4 source integrated to CSDA ranges at the nine energies; D-4: Geant4 ion
  stopping path). **D-4 preliminary (2026-11-11, same caveat as above):** in opt4 alpha uses the ICRU 90 alpha
  array directly (≤ 1 GeV, no along-step corrections); GenericIon uses the ICRU 90 proton array at equal velocity
  × effective charge² with the Lindhard–Sørensen correction and a Barkas term 2B(z−1)/z applied along the step
  (not in the range table); Geant4 has no ICRU 90 C/O arrays, so decision 0038's "integrated ICRU 90" C/O baseline
  must be re-examined (I3c rebuilds exactly this TOPAS-like construction and is the D-3 instrument). Candidate for
  the lengthening: steps losing ≤ 2 % of the energy use dE/dx at the pre-step energy (`G4VEnergyLossProcess`,
  limit 0.02 in `G4ionIonisation`), which lengthens the range by ≈ +0.03–0.06 % at 0.1 mm steps and up to
  ≈ +0.2 % at 1 mm steps (estimate); the LS + Barkas terms probably shorten it by ≈ 0.1–0.2 % (sign TO VERIFY);
  nuclear stopping, effective charge, interpolation and the R80-vs-CSDA definition are each ≤ 0.05 %. Geant4 ends
  ion tracks at 1 keV with a final-range step of 1 µm (ions) / 20 µm (alpha): no point deposit comparable to ionmc's
  2 MeV/u cutoff (ESTABLISHED). No ionmc stopping
  change: closing this gap would be fitting to the reference (plan D STOP 3). Low ends: ≈ 0.04 mm of the He 50
  deficit (0.065 mm) is the E_cut convention (estimate).
- **(c) z_peak FAILs at the three upper ends = (b) against the absolute 0.5 mm tolerance** (ESTABLISHED from the
  reports: z_peak shift / R80 shift = −0.81/−0.82 mm He 250, −0.74/−0.83 mm C 430, −0.66/−0.70 mm O 430). The
  absolute (c) tolerance is tighter than the relative (a) tolerance for R80 > 143 mm — a revision-1 design
  inconsistency, recorded, NOT changed (tolerances are never widened). The held-out ranges (≤ 78 mm) are not
  affected by it.
- **(d) TL-EM FAILs — estimator** (no self-consistent ±2σ window; results note). Every evaluable gated point lies in
  the band (r = 1.17–1.31); the lateral physics shows no disagreement at the band's resolution.
- Diagnostics to be committed as reviewed, report-only scripts before the held-out grading (none touches a sealed
  sample or any TOPAS output beyond the committed reports): **D-1** table CSDA range vs committed R80 per case;
  **D-2** the forward W/PE model (point vs continuous deposit, σ_R scan); **D-3** ICRU 90 array integration vs TOPAS
  R80; **D-4** P2 source reading (Geant4 11.4 ion fluctuation and ion stopping classes, via
  `acquire_reference_data`); **D-5** (optional, reference-internal) a scorer-free TOPAS smoke dumping the Geant4
  dE/dx and range tables for alpha and GenericIon(6,12)/(8,16) in G4_WATER, if TOPAS exposes the dump (TO VERIFY).
  D-3/D-4 inform the attribution only; they do not block the held-out set.

**D-4 record committed (2026-11-13, orchestrator; research specialist's reading; report-only, nothing re-graded).**
The ten Geant4 v11.4.2 source files of the preliminary D-4 reading (item 1, "D-4 preliminary") were re-acquired through
`acquire_reference_data` (role exploratory, Geant4 Software License; sha256 and acquisition ids in the record) and re-read:
`validation/plans/v3-008/d4-geant4-ion-fluctuation-stopping.md` (bound in the P9 record). Findings that change the
preliminary text: (i) the step regime of the TD rows is settled at 0.1 mm (plan step-length smoke, 97.5 % of IDD steps),
so the pre-step dE/dx candidate of item 1(b) contributes only ≈ +0.03–0.06 % (INFERRED) — the "up to ≈ +0.2 % at 1 mm
steps" clause of item 1(b) does not apply to these runs and is superseded by the record; (ii) alpha has no LS term
(G4BetheBlochModel branch), so the He R80 offset of +0.19–0.21 % is not an LS/Barkas effect; the leading remaining
candidate for the longer reference range at all three ions is the ICRU 90 base arrays of the Geant4 source against ionmc's
committed Bethe table (INFERRED; D-3 is the deciding instrument); (iii) the Bohr branch of the universal fluctuation
sampler has three conditions (with the 1 km cut still meanLoss ≥ 10·T_max), the handover charge and mass come from
`SetParticleAndCharge` before each step, and the preliminary line numbers of the LS along-step terms, the BetheBloch
alpha branch, the ICRU 90 vectors and the opt4 step functions are corrected; an archived TOPAS configuration printout
corroborates linLossLim 0.02, the step functions, LindhardSorensen for GenericIon and the 1 keV lowest energy. The
σ_R attribution of item 1(a) (Glandz-vs-Bohr variance rising with energy per nucleon; ≈ 0.98–1.02 at the held-out
energies) stands as INFERRED. Remaining TO VERIFY items and what resolves them are listed in the record's section 5.
No ionmc model, parameter or grade changes (plan D, STOP 3).

**2. Estimator (prospective; applies to rows not yet graded only).** r2-topas-target.md §4.2 receives the dated
"window-cycle" rule (exact text there): the 1e-4 relative convergence test is kept and applied first (every fit that
converged before returns a bit-identical σ and iteration count); if the window W_k equals a window used at an
earlier iteration j ≤ k − 2 (no self-consistent window exists), σ_core = arithmetic mean of the ML scales σ_j … σ_{k−1}
of the cycle's distinct windows, flagged with period and amplitude A = (max − min)/mean, valid iff A ≤ 2 % (cap fixed a priori: ≈ 2× the
single-part statistical SE of σ̂ at the smallest planned statistics and < 1/2 of the −5 % band edge); every other
non-termination, "fewer than 3 bins", non-contiguous window or bracket-edge ML stays invalid. Implementation
`validation/scripts/ion_em/r2_statistics.py` (`core_sigma`, `FIT_CYCLE_MAX_AMPLITUDE = 0.02`, `CoreFit.cycle_period`,
`CoreFit.cycle_amplitude`), per-point counts of cycled fits in the TL reports, tests
`tests/ionmc/test_r2_statistics_cycle.py`. It governs TD-H/TL-H (the only rows left to grade, item 3). The 18 EM
grades stay as graded. A **report-only re-evaluation** of TL-EM under the rule is permitted and will be committed at
`validation/results/transport/v8r2-topas/reports/report-only/TL-EM-case<N>-r2stop3.json` (`r2_evaluate.py
reeval-tl-em`, schema `ionmc-v3-008-r2-reeval-1`, no `status`/`verdict`, header "graded verdict stands"); it is an
attribution of the TL-EM FAILs, never a grade.

**3. Branch decision: B — model change after attribution, limited to the end-of-range deposition.**
- **Change (decision 0043 §4, dated amendment of the same date):** end-of-range model `continuous_residual_v1`
  (default for the ion sources He/C/O): when E ≤ E_cut the residual energy is deposited along the straight
  continuation of the current direction over ionmc's own table residual range R_tab(E), with the table's dE/ds
  (sub-stepped on the scoring grids), no straggling, no MCS and no secondaries on that ≤ 0.092 mm segment; the
  history ends as before (`code_cutoff`; `cutoff` counts the in-world residual, the part beyond the world boundary
  goes to the `escaped` tally; voxel-exact walk along the transport grid with up to 32 sub-steps per homogeneous
  leg, table inverse per sub-step — decision 0043 §4 dated precision of 2026-11-11). `point_v1` (the previous convention)
  stays selectable and is the proton default (`PhysicsOptions.end_of_range_model = "default"` resolves by charge;
  proton transport bit-identical). The requested and effective option are recorded in the run summary. No stopping
  change, no straggling change, no new parameter: nothing is fitted to TOPAS. Declared residual of the model
  (estimate): tail-profile dependence ≤ 2.3 % in PE and ≤ 0.008 mm in R80 at the low ends, ≤ 0.4 % / ≤ 0.004 mm
  at the held-out energies.
- **Consequences (plan D STOP 3, revision-2 note of 2026-11-07; r2-topas-target.md §5.4):** (i) all main TD/TL rows
  are **consumed**: TD-EM/TL-EM preserved as graded (3 + 1 PASS, 6 + 8 FAIL), TD-FULL/TL-FULL (main) never graded;
  (ii) **TD-H/TL-H** (cases 10–12: He 100, C 150, O 200 MeV/u; criteria TD-EM (a)–(e), TD-FULL (a), (c), TL-EM (a),
  TL-FULL (a); no D_ext) become the **sole gating TOPAS-agreement evidence**, with unchanged tolerances; (iii) the
  **new frozen SHA** is the reviewed commit carrying this amendment and the implementation (P9 freeze record
  rebound to it); the A-rows (CI, LV8/HR8 host rows, T1-ion, T12-ion) are re-run at that SHA before any held-out
  sample is produced (one-SHA rule, plan D PASS); (iv) the **held-out ionmc samples** are re-run at that SHA with
  the re-run base **20621004** (declared here before use; r_index 39–41 unchanged, `steps_v8.py --seed-base
  20621004`, seeds 20624904 + k, 20625004 + k, 20625104 + k; chunk plan as sealed: case 10: 4 × 2.5e5, case 11:
  4 × 2.5e5, case 12: 5 × 2e5 histories), their signed host records sealed in a new sample-seal record
  (`validation/results/transport/v8r2-samples/sample-seal-heldout-r2stop3.json`) by a reviewed commit, recorded in
  this plan by a dated held-out seal note (seal SHA + review id) **before any held-out TOPAS run**; the R2-TOPAS-4
  run-validity rule then binds the held-out TOPAS parts to that seal commit (ancestry; review finished before
  `reference_started`); the sealed main-set samples of `sample-seal.json` stay sealed (consumed, preserved); (v)
  **V5-SS** is not in the sample seal (its ionmc samples and TOPAS runs are not yet produced); it is produced at the
  new frozen SHA, unchanged design; (vi) **TOPAS:** the 48 held-out parts (cases 10–12, parts 0–15, seeds unchanged)
  are evidence; the 10 main IDD-FULL parts already run (He 50 parts 8–11, He 150 parts 8–11, He 250 parts 8–9) are
  preserved (records committed, outputs unread; TF-DIAG use only); the remaining 62 main FULL parts are **not
  required** and are deferred — orchestrator decision after the held-out grading whether to run them as TF-DIAG
  within the REF cap. **Loss:** TL (b) D_ext (the f_dM extrapolation trend, directive item 6) is never evaluated
  (TL-EM (b) NOT EVALUABLE and consumed); the adequacy of the f_dM extrapolation rests on TL-H (a) and I11-ext
  (declared limitation); the TOPAS-agreement claim narrows to He 100, C 150 and O 200 MeV/u plus the main EM rows
  that passed as graded.
- **Why B and not A.** The dominant cause of the low-end (d) FAILs is an identified ionmc model approximation whose
  effect the pre-registered defence misjudged; leaving it would knowingly spend the held-out set on the same defect.
  Branch A would end at the operator with gating FAILs that are partly an estimator defect. B is the pre-registered
  route for a model change after attribution. **Planning expectations (estimates, not predictions, recorded before
  the held-out runs):** TD-H (a)/(c) pass with margin (reference range excess ≈ 0.13–0.18 mm vs (a) tolerance
  0.19–0.27 mm; z_peak not binding at R ≤ 78 mm); TD-H (d) residual ≈ +1.2–2.1 % if TOPAS σ_R exceeds Bohr by 5 %
  — a FAIL there if the C/O excess is ≈ 10 % would be a legitimate finding against ionmc's straggling model, recorded
  as such.
- **The operator is informed through the ledger that** (i) the upper-end z_peak FAILs trace to a reference range
  0.15–0.45 % beyond its declared ICRU 90 lineage, which no non-fitted ionmc change can remove; (ii) the claim
  narrows to three interior energies without the f_dM trend test; (iii) the operator may override this branch choice
  (Branch A: integration decision on gating FAILs). No operator action is required to proceed.

**4. Budget (planning).** Host: 13 held-out chunk runs (≈ 3–4 h at the sealed per-chunk wall times) + the A-rows
(≈ 1–2 h); REF: 48 held-out parts (≈ 2–3 h at the measured t_s); EM evidence spent ≈ 1.5 h + 10 FULL parts ≈ 0.5 h;
all within the 17 h REF cap and the host budget of §E.

**4b. Gap recorded (2026-11-11, orchestrator; found while preparing the A-row re-run).** The LV- and HR-tier
verification rows of section A (I8, I8b, I9, I12, I13, I14, I15, T1-ion LV, T12-ion, T13-ion, R1-ion, R2 LV) were
never implemented as the `run_suite` suites `lv8`/`hr8` named in section F.2 and were never run; the CI-tier rows
ran in CI at every commit. Order step 0 ("implementation and the A-rows … then the ionmc TD/TL samples") was
therefore not followed for the LV/HR tier before the sample seal of 2026-11-09, an orchestrator omission. The
main-set samples are consumed under item 3 regardless. The suites are being implemented (`steps_v8.py`,
`run_suite.py`, `summarize.py`; lv8 base 20581004, hr8 base 20591004; r_index per step in `steps_v8.ROW_R_INDEX`,
42–55); implemented in this commit: lv8 `i6-lv-deterministic-end-depth` (r 51), `i8-range-straggling-{he150,c290,o350}`
(r 42; I8b reported in the same document), `i13-r80-vs-csda-{case}` (r 43, nine cases), `t1-ion-trace-parity-c290-o350`
(r 48, K = 256), `r1-ion-repeatability-c290` (r 50), `r2-a16-trace-fixture`; hr8 `t13-ion-chunks-cuda-c290` (r 55).
Implemented in the round-3 commit under item 4c as shard + combine steps (`steps_v8.py shard/combine`, lv5
pattern; partials sealed by content digest and bindings; fail closed on missing/altered/mis-bound partials): T12-ion
LV (r 49: python 2 × 2000 histories, warp-cpu f64 2 × 1e5), I9 (r 44: He 150, C 290, O 350, 2 × 1e5 each), I12
(r 45: θ_rms and σ at s_max 0.1/0.5/1/5 mm, shards of 2e4–1e5), I14 (r 46: six (s_max, f_E) configurations, shards
of 1.25e4–1e5), I15 (r 47: two voxel sizes × two shifts, shards of 5e4–1e5); hr8 T12-ion HR (r 54: cpu64 2e5 inside
hr8, cuda32 and cuda64 1e6, He 150 / O 350 cpu64 2e5 vs cuda32 1e6). Implementation definitions where the plan is open
(fixed now, before any of these rows runs): T12-ion uses the lv8 ion default step limit 1 mm on both backends (the
proton T12 used the 0.2 mm lateral bin); I9 gates |variance/A2 − 1| ≤ 2 % (the σ ratio is reported); I14 compares
each configuration with the reference (1 mm, 0.02) by the T9 IDD rules and checks |ΔR80| ≤ 0.05 mm over all pairs,
pooled shards rebatched to 40; I12 requires every θ_rms within 0.5 % of the quadrature at the exit plane; I15 gates
the pairwise checks only (θ/quadrature and σ/Fermi–Eyges reported); the f32-vs-f64 probe is cuda32 vs cuda64 at
1e6; CUDA samples use 100 batches, CPU/python 40. The ion Fermi–Eyges T_ion path is the I10 CI test's. Definitions fixed by the implementation
where the plan was open (documented in `steps_v8.py`): I8 end depth = depth at which E reaches E_cut (before the
residual deposit), compared with the I4b quadrature from E0 to E_cut (corrections on, Bohr variance with the bare
charge); I8b compares with a proton of equal CSDA range / √(M/m_p) and, for carbon, the 0.3–0.4 % literature band
(a deviation beyond ±15 % is a finding); I13 uses 0.1 mm IDD bins in a 1.1 R box and the `parity.r80_of` estimator;
every row runs warp-cpu float64 with the end-of-range default; R1-ion and T13-ion history counts are the
implementation's (1e4 per precision; 1e5 per chunk size). Order: the cheap A-rows and the STOP rows T1-ion and
T12-ion run at the frozen SHA BEFORE the held-out samples; the remaining LV rows (I9, I12, I14, I15) and hr8 run at
the same SHA before the PASS declaration (one-SHA rule), scheduled after the held-out grading because the host
executes one job at a time.

**4c. Host budget correction (2026-11-11, orchestrator decision; prospective).** Section E planned warp-cpu at
≈ 7·10³ histories/s. Measured: ≈ 118 histories/s for C 290 MeV/u, all EM on, 1 mm steps, one process (sandbox,
lv8 implementation), and 67–270 histories/s on the host for the sealed TD/TL chunks (He 50 268, He 150 149, He 250
84, C 80 227, C 290 112, C 430 67, O 100 206, O 350 98, O 430 79, He 100 206, C 150 196, O 200 172; from the sealed
`wall_s`), i.e. ≈ 60× slower than planned: the ≈ 17 × 1e6-history LV rows would cost ≈ 40 h of host time instead
of 41 min, and each 1e6 run needs ≥ 7 shards of ≤ 1500 s. Decision (no tolerance changes; the rows keep their
observables, geometries and criteria): the warp-cpu history counts of the LV falsification/self-consistency rows
are set to **2e5 per configuration** — I9 2e5 per standard point, I12 2e5 per s_max, I14 2e5 per (s_max, f_E),
I15 2e5 per voxel configuration, T12-ion warp-cpu f64 2e5 (python 4·10³ in 40 batches unchanged) — each as shards of
≤ 1500 s combined by the step; the CUDA rows (T12-ion cuda f32, the f32-vs-f64 probe, He 150 / O 350 cpu vs cuda,
T13-ion) keep 1e6 on the GPU where the planned throughput holds (to be measured by the first hr8 step; fallback 2e5
by the same rule). Statistical adequacy at 2e5 (planning): the relative SE of a projected rms over N histories is
≈ 1/√(2N) ≈ 0.16 %, so the pairwise 0.5 % θ_rms criterion of I12/I15 is resolved at ≈ 3 SE and the 1 % σ criterion
at ≈ 4 SE; the I9 2 % criterion at ≈ 10 SE; the I14 R80 SE ≈ σ_R/√N ≈ 1.5 µm against 0.05 mm, the grouped χ²
keeps its p > 0.001 rule at the reduced N; T12-ion's power is bounded by its python sample (4·10³) in any case.
Estimated host time: the sandbox rates at the smaller step limits and energy fractions are 3–6× lower than the
1 mm rate (C 290, all EM on, histories/s: 1 mm/0.02 → 111, 1 mm/0.005 → 35, 0.5 mm/0.02 → 78, 0.5 mm/0.005 → 34,
0.1 mm/0.02 → 23, 0.1 mm/0.005 → 18), so the sharded lv8 rows plan to ≈ 22 h of host time (T12-ion ≈ 1.1 h, I9
≈ 1.4 h, I12 ≈ 6.4 h, I14 ≈ 9.9 h, I15 ≈ 3.5 h; every shard ≤ 1000 s planned against the 1500 s timeout) plus ≈ 1 h
for the cheap rows and hr8's CPU parts ≈ 1.4 h; the held-out samples ≈ 4 h (13 chunks at the sealed rates). The §E
total is superseded accordingly; REF is unchanged. These rows gate the task PASS (one-SHA rule) but not the
held-out blinding, so they are scheduled after the held-out grading; the orchestrator informs the operator of the
≈ 22 h through the ledger and may ask for a decision if the host time is not available.

The report-only TL-EM re-evaluation of item 2 is committed at
`validation/results/transport/v8r2-topas/reports/report-only/TL-EM-case<N>-r2stop3.json` (evaluator HEAD e6f17bc2,
9 files): with the window-cycle rule every gated point of the main set lies inside the band (criterion (a) would
read PASS for 8 of 9 cases; O 100 keeps its undefined f = 0.25 point); 0–9 of 80 ionmc batch fits and 0–1 of 8 TOPAS
part fits per point used the cycle rule. The graded TL-EM verdicts stand.

**4d. Held-out samples produced (2026-11-11, orchestrator; before the held-out seal commit).** The A-rows of item 4b
(lv8 steps 01–21: I6-LV, I8, I13, T1-ion K = 256, R1-ion, R2, T12-ion LV) passed on the host at the frozen SHA
f5f3a701 (validation record VAL-20261111-032337-9F2C85; Codex REVIEW-ab510ca6; CI run 42943469701). The 13 held-out
chunks then ran on the host at f5f3a701 with `steps_v8.py td-tl-sample --seed-base 20621004` (r_index 39–41, seeds
20624905–20624908, 20625005–20625008, 20625105–20625109; case 10: 4 × 2.5e5, case 11: 4 × 2.5e5, case 12: 5 × 2e5
histories; every run succeeded with all counters zero; wall 4.38 h in total): `RUN-20261111T032108Z-0f2c091b` (case 10 k1), `RUN-20261111T034144Z-3d718129` (case 10 k2), `RUN-20261111T040210Z-bed1ad5f` (case 10 k3), `RUN-20261111T042235Z-eb5e4f67` (case 10 k4), `RUN-20261111T044303Z-10ad7312` (case 11 k1), `RUN-20261111T050431Z-9d8383b2` (case 11 k2), `RUN-20261111T052557Z-6ca8729d` (case 11 k3), `RUN-20261111T054739Z-333f15dd` (case 11 k4), `RUN-20261111T060908Z-fa584e14` (case 12 k1), `RUN-20261111T062857Z-6b98f3f0` (case 12 k2), `RUN-20261111T064831Z-daec24ca` (case 12 k3), `RUN-20261111T070805Z-71ec22bc` (case 12 k4), `RUN-20261111T072739Z-2ca81843` (case 12 k5). Their signed host records are
sealed in `validation/results/transport/v8r2-samples/sample-seal-heldout-r2stop3.json` (schema
`ionmc-v8r2-sample-seal-heldout-1`; `seed_base` 20621004, `npz_subdir` `r2-samples-heldout`, `implementation_sha`
f5f3a701; the main-set seal `sample-seal.json` stays as sealed). The evaluator reads held-out cases against this
record only (`r2_io.SealKind`; plan block keys `seal_sha_heldout` / `seal_review_id_heldout`; review record
`records/review-heldout.json`; fail closed on any mismatch). No held-out TOPAS output exists yet; the held-out seal
note follows in the commit after this seal commit's review, and the 48 held-out TOPAS parts start only after that.

**4e. Dated held-out seal note (2026-11-11, orchestrator; written after the independent review of the held-out seal
commit finished and BEFORE any held-out TOPAS run).** The held-out seal commit is `ca21ec18601c8558082f3f0d543c07d40204992b`
(seal record `validation/results/transport/v8r2-samples/sample-seal-heldout-r2stop3.json`, implementation f5f3a701);
its review REVIEW-557971cc64724e7684141671dfac057f passed (finished 2026-11-11T07:54:30.313805Z; record
`validation/results/transport/v8r2-topas/records/review-heldout.json`). Every held-out TOPAS run must start after
that time at a `code_sha` descending from the seal commit (R2-TOPAS-4 applied to the held-out seal; evaluator
`r2_io.SealKind heldout`).

```json
{"seal_sha_heldout": "ca21ec18601c8558082f3f0d543c07d40204992b", "seal_review_id_heldout": "REVIEW-557971cc64724e7684141671dfac057f"}
```

**4f. Held-out TOPAS evidence runs produced (2026-11-11, orchestrator; records commit, before any grade).** The 48
held-out parts (cases 10–12; per case IDD-EM parts 0–3 at 2.5e4, LAT-EM parts 4–7 at 2.5e4, IDD-FULL parts 8–11 at
the re-sized N_part 5000 / 3200 / 2500 for He / C / O, LAT-FULL parts 12–15 at 5000 histories; bundles frozen in
`frozen-bundles.json`) ran through the controlled reference runner in the pre-registered order (case 10, 11, 12; parts
0–15), all at `code_sha` 55e6a18041710ed1f97d2a6bb44b97ea21daaf12 (the seal-note commit, a descendant of the
held-out seal ca21ec18); the first run started 2026-11-11T08:02:33.8Z, after the held-out review finished
(07:54:30.3Z, item 4e), and the last finished 10:02:32.6Z; every run succeeded (exit 0, empty stderr, no missing
output); wall 0.24 h (He) + 0.55 h (C) + 0.88 h (O) = 1.67 h in total. Their signed records (`record.json`, signer
`reference-runner`, key 303cce4c129d387c) and verbatim artifact listings (`artifacts.json`) are committed under
`validation/results/transport/v8r2-topas/records/runs/<run_id>/` and indexed in `records/index.json` (130 entries =
82 main + 48 held-out). Run validity (`r2_evaluate.py check-run`, held-out chronology) and the grades of TD-H / TL-H
(rows TD-EM, TD-FULL, TL-EM, TL-FULL of cases 10–12, each once) follow in the next commit; nothing has been graded at
this commit.

**4g. Evaluator amendment — signed host-run records of the sample chunks (2026-11-11, orchestrator; Codex
REVIEW-2f0b4eec on the seal-note commit 55e6a180, finding 1 "important").** The review found that the grader bound the
ionmc samples only through the sealed npz hashes and the attestation *block* copied into the seal, without
re-verifying the host runner's signature over the complete run record. From this commit `r2_evaluate.py grade` and
`reeval-tl-em` (both seal kinds) require, for every sealed chunk of a graded case and before any number is read, the
complete signed host-run record committed at `validation/results/transport/v8r2-samples/host-records/<run_id>.json`
(read as committed at HEAD): the record must verify against the committed trust file with signer `host-runner`
(key 4d02f521f3fb9344); its attestation block must be the sealed `host_runner_attestation` (exact for the held-out
seal; the main seal carries a five-key subset — signer, key_id, signed_at_utc, payload_sha256, signature — which is
matched key by key, the signature itself binding the whole payload); its signed fields must state a clean successful
run of V3-008 at the seal's `implementation_sha` (`succeeded`, `exit_code` 0, not timed out, clean tree before and
after, no publication error) and agree with the seal (`run_id`, start/finish times, `stdout_sha256` =
`wrapper_stdout_sha256`, `output_manifest_sha256` = the chunk's and the seal's single `host_runs` entry); the signed
output manifest is reproduced from the published run directory (`validation/generated/<run_id>/`, sorted
`validation/generated/<path>\0<sha256>` lines) and must contain the chunk's npz with the sealed sha256. Any gap is an
ERROR (fail closed; `r2_io.verify_host_records`; tests in `tests/ionmc/test_r2_evaluate.py`). The complete records
live in the protected host-run store, which the agent cannot read; the orchestrator requested their verbatim copies
from the operator (IR-20261111-101840-EDCAE6, category inaccessible-information). Until the records of the 13 held-out
runs (and the 76 main runs, for the report-only re-evaluation) are committed, every held-out grade and every
re-evaluation fails closed with "no committed signed host-run record"; no grade is affected because none of cases
10–12 has been graded. Finding 2 (minor) is closed by an unconditional digest test of the D-3 aggregate
(`tests/ionmc/test_d3_range_attribution.py`, no ICRU cache needed).

**4h. Evaluator amendment, round 2 (2026-11-11, orchestrator; Codex REVIEW-e844de8a on b18fa715, finding 1 "important").**
The host-record binding of item 4g additionally binds the signed invocation (`argv`: the prescribed `steps_v8.py
td-tl-sample` call or the `python -c` chunk wrapper, with the chunk's case, k and `scale`, the kind's seed base and the
output directory; `r2_io.sample_invocation`, fail closed on any other shape) and the manifest-bound step document
(the step's published stdout, sha256 = the chunk's `step_stdout_sha256`, carrying exactly one step JSON document
whose case, k, seed, `seed_base`, `r_index`, histories, npz name/sha256, `git_sha`, `pass`, `valid`, all-zero
counters and effective n_histories/seed/cpu_workers agree with the seal; `r2_io.step_document`); negative tests for
wrong case, k, seed, scale, seed base, history count, counters, document and argv shapes. Finding 2 (minor): the
`validation/results/transport/README.md` paragraph now states that the host records are absent at this commit
(IR-20261111-101840-EDCAE6 open) and a cache-free inventory test (`tests/ionmc/test_v3_008_host_records.py`) verifies
every committed record against the seals and the trust file; its completeness test is a strict expected failure
until the records are committed. All 99 sealed chunks' published outputs satisfy the manifest and step-document
bindings at this commit (checked locally; the `argv` binding needs the records).

**4i. Host records committed; wrapper templates (2026-11-11, orchestrator; IR-20261111-101840-EDCAE6 resolved; Codex
REVIEW-84382c07 on eb29618d, finding "important").** The operator copied the complete signed host-run records of all
V3-008 host runs (106 files) into an agent-readable directory; the 89 records of the sealed sample runs (76 main
+ 13 held-out) each verify against the committed trust file (signer `host-runner`, key 4d02f521f3fb9344) and agree
with the seals (run id, `implementation_sha`, attestation block, start/finish times, stdout and output-manifest
digests; 0 problems) and are committed verbatim at `validation/results/transport/v8r2-samples/host-records/<run_id>.json`
(the 17 further records are A-row and other host runs, not part of the seals, and are not committed). The review found
that the round-2 invocation check (item 4h) matched the wrapper source by pattern and therefore tolerated additional
code; from this commit the signed `argv` of a sample run must be exactly `["python", "-c", <source>]` with the source
byte-identical to one of four reviewed wrapper templates (`validation/scripts/ion_em/r2_wrapper_templates.py`:
T1/T2 single-chunk main wrappers, T3 the multi-chunk main wrapper with the k list, T4 the held-out wrapper with
`--seed-base 20621004`; each template rendered with the run's sealed case, k value(s) and scale, so the k list of a
T3 run must equal the sealed chunks of that run), every other source fails closed. The 89 real records match the
templates 1 / 30 / 45 / 13 (T1–T4). The cache-free inventory test (`tests/ionmc/test_v3_008_host_records.py`) now
requires a record for every sealed run and checks each record's signature, seal binding and template match. With
the records committed the held-out grader is complete: TD-H / TL-H (rows TD-EM, TD-FULL, TL-EM, TL-FULL of cases
10–12) are graded once in the next step, at this evaluator, and reported in the following commit.

**4j. One-SHA rule applied to the remaining LV/HR rows (2026-11-11, orchestrator; prospective, before any of these
rows runs).** Item 4b scheduled the remaining lv8 rows (I9, I12, I14, I15; steps 22–135) and hr8 (T13-ion, T12-ion HR)
"at the same SHA" as the A-rows (f5f3a701). The controlled host runner executes the task worktree's current clean
HEAD only, and HEAD has moved past f5f3a701 by evaluator, test, record and plan commits only: `git diff --stat
f5f3a701 HEAD -- src validation/scripts/transport pyproject.toml uv.lock` is empty at the commit carrying this note,
i.e. the implementation and the suite scripts are byte-identical to the frozen implementation SHA. The remaining
rows therefore run at the HEAD of the commit carrying this note (recorded with every host run id in the host ledger
and in the dated results note that follows), with `--expected-sha` set to that HEAD; the sharded rows and their
combine steps run at that one HEAD (the partials bind the SHA), and no commit is made to the task branch while the
block is in flight. The tolerances, observables, history counts (item 4c) and seeds (lv8 base 20581004, hr8 base
20591004, r_index per `steps_v8.ROW_R_INDEX`) are unchanged. This satisfies the one-SHA rule in substance (one
implementation) and is recorded here because the letter of item 4b named f5f3a701.

**4k. I9 FAIL at e602a583 and its attribution; amendment I9-REF (2026-11-11, orchestrator; prospective for the
re-run; the FAIL is preserved).** The lv8 row I9 (lateral variance of deposits vs Fermi–Eyges A2 with the implemented
T_ion, |variance/A2 − 1| ≤ 2 %) ran at e602a583 (item 4j; shards 22–27, host runs RUN-20261111T124822Z-f978c3d3,
RUN-20261111T125420Z-301ed67a, RUN-20261111T130009Z-2aee46ff, RUN-20261111T130648Z-f00c3ca4,
RUN-20261111T131327Z-d3458298, RUN-20261111T132114Z-876468a3, 2e5 histories per standard point; combine 28 =
RUN-20261111T132917Z-098c9e9a, exit 1) and **FAILS**: variance/A2 = 1.733–1.740 (σ/√A2 = 1.317–1.319) for He 150,
C 290 and O 350 at every depth fraction 0.25–0.9 — a factor uniform in ion, energy and depth. **Attribution (read-only
code reading, no re-grading):** the row's reference helper `_IonPath` (`tests/ionmc/test_v3_008_rows_ci.py`, used by
`steps_v8.ion_path` for I9, I12 and the informative I15 ratios) evaluates the implemented ion scattering power with a
hard-coded inverse scattering length 1/X_S = 0.0123 cm²/g (a placeholder carried over from the ratio-only CI tests,
commit 53424cf), whereas the transport kernel uses the material value `inverse_scattering_length_cm2_per_g(WATER)` =
0.021331 cm²/g; the predicted ratio 0.021331/0.0123 = 1.7342 reproduces the twelve measured ratios to −0.08 … +0.35 %.
The transport is therefore consistent with the Fermi–Eyges prediction of its own T_ion to ≤ 0.35 % at 2e5 histories;
the FAIL is a defect of the verification row's reference constant, not of the implementation under test. (The TOPAS
TL ratios r = 1.18–1.31 are not a constant factor and keep their attribution of item 1(d): a model difference,
Gaussian core vs Geant4's single-scattering halo.) The proton references (`transport/mcs_checks.py`) use the material
value and are unaffected; I10 (CI) is ratio-only and passes with either constant; I14, T12-ion and T13-ion do not use
the helper. **Record:** I9 at e602a583 is FAIL in the evidence table (verbatim, never re-graded); the I12 shard 29
started at e602a583 (base 20581004) is preserved as an abandoned partial (no row graded from it); no further lv8 row is
run at e602a583. **Amendment I9-REF (reviewed with its implementation commit before any re-run):** (i) `_IonPath` and
the CI helpers take 1/X_S from the material table (water 0.021331 cm²/g) and a CI test pins the helper's constant to
the table value; the I9/I12/I15 definitions, tolerances, geometries, bins and history counts are unchanged; (ii) the
remaining lv8 rows I9, I12, I14, I15 run at the SHA of the amendment commit (one-SHA rule for the re-run block) with
the lv8 **re-run base 20572004** (declared here before use; section C amended: the second declared re-run base
20631004 supports r_index ≤ 9 only, while I9/I12/I14/I15 use r_index 44–47, and the leading gap 20572004–20581003 of
the V3-008 reservation was never assigned — the new family `rerun3` = 20572004 + 100·r + k, r ≤ 60, k < 100, lies
entirely inside the reservation and is disjoint from every other V3-008 and V3-005 family, so no cross-task ledger
block is requested; guard test `test_v3_008_seed_families.py` extended accordingly; `--seed-base 20572004`), i.e.
the re-run I9 is a fresh sample, not a re-grade of the consumed one; hr8
(T13-ion, T12-ion HR) runs at the same SHA with its untouched base 20591004; (iii) the A-rows of item 4d (steps 01–21
at f5f3a701, base 20581004) stand. The corrected lateral grid half-width (8·√A2 with the material constant) is the
only computational change of the rows.

Implementation of I9-REF (same commit): `IonPath` moved to the new helper module
`src/ionmc/transport/ion_mcs_checks.py` (verification helper only; `git diff f5f3a701 -- src/ionmc` shows this added
file and nothing else — no transport kernel, physics table or option changed), the material 1/X_S and density are its
defaults, the CI tests use the material value and pin it, `steps_v8.ion_path` imports the helper, and
`summarize.V8_RERUN_SEED_BASES` accepts the lv8 re-run base 20572004 as qualifying (plan decision of this item; the
consumed base 20581004 stays the base of the A-rows of item 4d).

**5. Order from here (pre-registered).** Implementation + estimator rule + this amendment in one reviewed commit (new
frozen SHA) → A-rows at that SHA → held-out samples (13 chunks) → held-out seal commit + review + dated seal note →
48 held-out TOPAS parts (IDD-EM, LAT-EM, IDD-FULL, LAT-FULL per case) → TD-H, TL-H graded once → V5-SS → report-only
TL-EM re-evaluation and D-1…D-4 (any time after this commit; they read only committed inputs).

### Dated results note — held-out set TD-H / TL-H (2026-11-11, orchestrator; amendment R2-STOP3 item 5: the held-out rows graded once)

**Evidence.** The 48 held-out evidence parts (cases 10 He 100, 11 C 150, 12 O 200; IDD-EM parts 0–3 and LAT-EM parts
4–7 at 2.5e4 primaries, IDD-FULL parts 8–11 at 5000 / 3200 / 2500, LAT-FULL parts 12–15 at 5000) ran sequentially on
2026-11-11 08:02–10:02 UTC at `code_sha` 55e6a180 (the held-out seal-note commit; the held-out seal ca21ec18 is an
ancestor; REVIEW-557971cc finished 07:54:30 UTC, item 4e). Their signed schema-3 run records (signer
`reference-runner`, key_id 303cce4c129d387c) and verbatim artifact listings are committed under
`validation/results/transport/v8r2-topas/records/` (commit b18fa715, item 4f); every listed file was verified
against its sha256 and the signed `output_files_sha256` reproduces from the `work/` entries. `r2_evaluate.py
check-run` at b18fa715: 48/48 VALID (binding, held-out chronology, frozen manifest, module markers). The sealed
ionmc samples are the 13 held-out chunks of `sample-seal-heldout-r2stop3.json` (implementation f5f3a701, re-run
base 20621004), each bound to its complete signed host-run record (items 4g–4i). The four rows were then graded,
once each, with `r2_evaluate.py grade` at evaluator HEAD 1fc287cf (2026-11-11 11:58–11:59 UTC; the grade
re-checks run validity, chronology, the replicate gate and the host-record binding before reading any number);
the reports (statistics, record ids and sha256 only) are `validation/results/transport/v8r2-topas/reports/<row>-case<N>.json`
for rows TD-EM, TD-FULL, TL-EM, TL-FULL and cases 10–12, bound by sha256 in the P9 freeze record. **They are final:
nothing in them is re-graded.**

**Verdict: TD-H PASS (TD-EM 3/3, TD-FULL 3/3) and TL-H PASS (TL-EM 3/3, TL-FULL 3/3)** — all twelve (case, row)
grades PASS; every gated criterion passes in every case.
**TD-EM (held-out):**

| Case | Status | (a) R80 ionmc / TOPAS [mm] (stat ≤ tol) | (b) W ionmc / TOPAS [mm] (report-only) | (c) z_peak ionmc / TOPAS [mm] (stat ≤ 0.5) | (d) PE ionmc / TOPAS (stat ≤ tol) | (e) plateau groups (failing) |
|---|---|---|---|---|---|---|
| He 100 | **PASS** | 77.95 / 78.12 (0.18 ≤ 0.27) PASS | 0.570 / 0.550 (+3.6 %) | 77.54 / 77.71 (0.19) PASS | 9.44 / 9.42 (0.04 ≤ 0.28) PASS | 30 (0) PASS |
| C 150 | **PASS** | 53.10 / 53.20 (0.11 ≤ 0.19) PASS | 0.215 / 0.216 (-0.4 %) | 52.86 / 53.01 (0.15) PASS | 10.95 / 11.06 (0.14 ≤ 0.33) PASS | 20 (0) PASS |
| O 200 | **PASS** | 65.47 / 65.57 (0.10 ≤ 0.23) PASS | 0.226 / 0.232 (-2.5 %) | 65.27 / 65.34 (0.07) PASS | 11.51 / 11.58 (0.09 ≤ 0.35) PASS | 25 (0) PASS |

**TD-FULL (held-out):**

| Case | Status | (a) R80 ionmc / TOPAS [mm] (stat ≤ tol) | (b) W ionmc / TOPAS [mm] (report-only) | (c) z_peak ionmc / TOPAS [mm] (stat ≤ 0.5) | (d) PE ionmc / TOPAS (stat ≤ tol) | (e) plateau groups (failing) |
|---|---|---|---|---|---|---|
| He 100 | **PASS** | 77.95 / 78.11 (0.21 ≤ 0.39) PASS | 0.570 / 0.580 (-1.7 %) | 77.54 / 77.66 (0.15) PASS | 9.44 / 6.65 (2.84 ≤ 0.20) REPORT-ONLY (report-only) | 30 (24) REPORT-ONLY (report-only) |
| C 150 | **PASS** | 53.10 / 53.21 (0.12 ≤ 0.24) PASS | 0.215 / 0.224 (-4.0 %) | 52.86 / 52.99 (0.13) PASS | 10.95 / 7.88 (3.29 ≤ 0.24) REPORT-ONLY (report-only) | 20 (15) REPORT-ONLY (report-only) |
| O 200 | **PASS** | 65.47 / 65.55 (0.09 ≤ 0.29) PASS | 0.226 / 0.260 (-13.0 %) | 65.27 / 65.32 (0.06) PASS | 11.51 / 8.18 (3.53 ≤ 0.25) REPORT-ONLY (report-only) | 25 (21) REPORT-ONLY (report-only) |

**TL-EM (held-out):**

| Case | Status | f = 0.25 | f = 0.5 | f = 0.75 | f = 0.9 | (b) D_ext |
|---|---|---|---|---|---|---|
| He 100 | **PASS** | σ_i/σ_T 1.208, gated: PASS | σ_i/σ_T 1.179, gated: PASS | σ_i/σ_T 1.181, gated: PASS; cycled fits 4 | σ_i/σ_T 1.178, gated: PASS; cycled fits 1 | n/a (held-out: no D_ext) |
| C 150 | **PASS** | σ_i/σ_T 1.311, report-only: REPORT-ONLY | σ_i/σ_T 1.267, gated: PASS | σ_i/σ_T 1.257, gated: PASS | σ_i/σ_T 1.253, gated: PASS; cycled fits 1 | n/a (held-out: no D_ext) |
| O 200 | **PASS** | σ_i/σ_T 1.320, report-only: REPORT-ONLY | σ_i/σ_T 1.280, gated: PASS | σ_i/σ_T 1.270, gated: PASS | σ_i/σ_T 1.266, gated: PASS; cycled fits 2 | n/a (held-out: no D_ext) |

**TL-FULL (held-out):**

| Case | Status | f = 0.25 | f = 0.5 | f = 0.75 | f = 0.9 | (b) D_ext |
|---|---|---|---|---|---|---|
| He 100 | **PASS** | σ_i/σ_T 1.210, gated: PASS | σ_i/σ_T 1.165, gated: PASS | σ_i/σ_T 1.140, gated: PASS; cycled fits 4 | σ_i/σ_T 1.127, gated: PASS; cycled fits 2 | n/a (held-out: no D_ext) |
| C 150 | **PASS** | σ_i/σ_T 1.318, report-only: REPORT-ONLY | σ_i/σ_T 1.285, gated: PASS | σ_i/σ_T 1.259, gated: PASS | σ_i/σ_T 1.252, gated: PASS; cycled fits 1 | n/a (held-out: no D_ext) |
| O 200 | **PASS** | σ_i/σ_T 1.308, report-only: REPORT-ONLY | σ_i/σ_T 1.279, gated: PASS | σ_i/σ_T 1.255, gated: PASS | σ_i/σ_T 1.258, gated: PASS; cycled fits 2 | n/a (held-out: no D_ext) |

**Findings (no re-grading; attribution only).** TD-EM: (a) R80 passes with ionmc again systematically short, by
−0.21 % (He 100), −0.20 % (C 150), −0.15 % (O 200) — the same sign and size as in the main set (−0.2 … −0.3 % at the
standard points), consistent with the D-3/D-4 attribution of item 1(b) (reference-side pre-step dE/dx and the
C/O stopping construction; no ionmc stopping change); (c) z_peak is 0.07–0.17 mm proximal, inside the absolute
0.5 mm tolerance at these ranges (53–78 mm); (d) PE agrees to +0.1 % (He 100), −1.0 % (C 150), −0.6 % (O 200),
i.e. the ionmc peak excess of the main low ends (+7.8 … +11.4 % at He 50, C 80, O 100 under the point deposit) is
absent at these energies under `continuous_residual_v1`, and the report-only distal width W is within +3.6 % /
−0.4 % / −2.5 % of TOPAS (main low ends: −5 … −17 %); (e) every plateau group passes. The held-out set contains no
case at the main set's low ends (50–100 MeV/u) or upper ends (430 MeV/u): the low-end PE/W behaviour and the
upper-end z_peak behaviour of the amended model are therefore not re-measured by this set, and the main-set
FAILs stand as graded (item 0). TD-FULL (gated (a), (c) only, §5.4): R80 and z_peak pass with the same −0.12 … −0.21 %
range offset; the report-only (d) and (e) show the expected signature of nuclear physics absent from ionmc — the
TOPAS-FULL peak is 28–30 % lower than the ionmc (EM-only) peak (ionmc/TOPAS − 1 = +42 %, +39 %, +41 %) and the plateau
rises along the depth in TOPAS (24/30, 15/20, 21/25 groups outside the 2 % band, all inconclusive-to-high), the
TF-DIAG attribution of §5.5 (V3-009 scope, no grade affected); W is −1.7 % / −4.0 % / −13.0 % (the O 200 FULL peak is
broadened by fragments). TL-EM and TL-FULL: at the ten gated points r = σ_core,ionmc / σ_core,TOPAS-EM is 1.178–1.280
(EM) and 1.127–1.285 (FULL), every point inside the band [−5 %, +40 %] of §5.3; the shallow f = 0.25 points of C 150
and O 200 are report-only by the gate-σ rule (r = 1.31–1.32); the dated §4.2 window-cycle rule (item 2) resolved
8 (EM) and 9 (FULL) batch/part fits that entered the 2-cycle at a bin edge — no point is invalid, and every
re-measured width lies in the same 1.13–1.32 band as the main-set re-evaluation (reports/report-only). D_ext is
not defined for the held-out set (§5.4).

**Consequence (plan D, revision-2 block; amendment R2-STOP3 item 3).** The held-out set was the sole gating
TOPAS-agreement evidence after the model change; it passes in full. The remaining gates of V3-008 before the
integration decision are unchanged: the lv8/hr8 rows not yet run (I9, I12, I14, I15 sharded and T12-ion HR, T13-ion;
items 4b–4c), V5-SS, the D-4 record, and the CI/review/validation gates of the reports commit. No clinical claim.

**Archive cross-check (2026-11-11, after the operator integrated the host-run record retrieval tools, develop 930b7a9).**
`list_host_run_records` lists 107 V3-008 host runs, 106 with a complete archived record (the first smoke run
RUN-20261107T224820Z-92f9f0be has none — a pre-existing gap, not a sealed run); the archive byte hashes of the
sampled sealed records (T1/T2/T3/T4 examples and one more held-out run) equal those of the operator's copies, and the
committed records remain bound by their verified signatures. No record was re-materialized or re-signed.

### Dated results note — lv8 re-run block (I9, I12, I14, I15) and hr8 (T13-ion, T12-ion HR) at 712f484a (2026-11-13, orchestrator; items 4j and 4k)

**Evidence.** The remaining lv8 rows ran as the pre-registered shard/combine steps 22–135 of `run_suite.py --suite lv8
--single-process --workers 1 --seed-base 20572004 --only <step>` (one host call per step; amendment I9-REF, item 4k:
re-run base 20572004, family `rerun3`, r_index 44–47) and the hr8 suite as its steps 01–12 (`--seed-base 20591004`,
r_index 54 and T13-ion), all at the clean HEAD `712f484a13dbee4a181d15addf127b1c58dc1cbc` (the amendment I9-REF commit;
`--expected-sha` set to it; no commit to the task branch while the block was in flight, item 4j). The block ran
2026-11-11 14:07 UTC – 2026-11-13 03:38 UTC (lv8 evidence wall 17.47 h, hr8 1.70 h; every evidence run `exit 0`,
`succeeded`, `committed_state_unchanged`, `dirty_before/after` false). Every sharded row's `PARTIAL <name> <sha256>`
stdout lines were taken from the signed host-runner records; each copied partial was re-hashed against that digest
before staging into the import archive whose manifest the combine step read (summaries `attestation.partials`, 110 lv8
+ 10 hr8, each with its `host_run_id`). **The complete signed host-runner record of every run of the block (136 files,
signer `host-runner`, key 4d02f521f3fb9344) is committed verbatim at
`validation/results/transport/v8-lvhr/host-records/<run_id>.json`** with a derived index (`v8-lvhr/index.json`); the
cache-free test `tests/ionmc/test_v3_008_lvhr_host_records.py` verifies every signature, requires the signed `argv` to
be byte-identical to the one wrapper template rendered with the run's suite, SHA, step, declared seed base and import
list, re-derives the index, and binds every partial and every step of the two committed combined summaries
(`validation/results/transport/lv8-712f484a-steps22-135-single-process.json`, `summarize.py --combine --subset` over
the 114 archives, `pass` true, `complete` false because steps 01–21 are the A-row block at f5f3a701 of item 4d;
`hr8-712f484a-single-process.json`, the complete hr8 suite, `pass` true) to the record that produced it. Both
summaries are `conformant: false` by code for the recorded reasons only (single-process diagnostic mode; combine steps
reading imported partials; the lv8 subset). The eight records of the preserved first I9 run and the abandoned I12
shard at e602a583 (item 4k) are committed in the same directory (roles `preserved_fail`, `preserved_fail_shard`,
`abandoned_partial`). Grades below are read from the combine step documents; **they are final and nothing is
re-graded.**

**Verdict: I9 PASS, I12 PASS, I14 PASS, I15 PASS (lv8); T13-ion PASS, T12-ion HR PASS (hr8).** With the A-rows of
item 4d (I6-LV, I8, I13, T1-ion, R1-ion, R2, T12-ion LV at f5f3a701) every LV/HR verification row of section A has
now passed, each at one clean SHA per block.

**I9 (lateral variance vs Fermi–Eyges A2 with the implemented T_ion; tolerance 2 %; 2e5 histories per standard point,
item 4c).** Shards 22–27: `RUN-20261111T140753Z-4ac10d13`, `RUN-20261111T141353Z-8b17b0cc`, `RUN-20261111T141948Z-b2cd0dee`,
`RUN-20261111T142631Z-a6420bf7`, `RUN-20261111T143323Z-f02c0468`, `RUN-20261111T144115Z-ce06b608`; combine 28
`RUN-20261111T144917Z-4c7e600d`. variance/A2 at z/R = 0.25 / 0.5 / 0.75 / 0.9: He 150 0.9993 / 0.9962 / 0.9963 / 0.9974;
C 290 1.0001 / 1.0012 / 1.0000 / 0.9997; O 350 0.9999 / 0.9954 / 0.9967 / 0.9995 — every point within 0.5 % of unity
(largest deviation −0.46 %, O 350 at 0.5 R), against the 2 % tolerance. **PASS.** The preserved FAIL at e602a583
(1.733–1.740, item 4k) is thereby confirmed as the reference-helper constant: the fresh sample with the material 1/X_S
reproduces A2 to ≤ 0.5 %, and no transport code changed between the two runs (`git diff e602a583 712f484a -- src/ionmc`
adds only the helper module `ion_mcs_checks.py`).

**I12 (transport MCS step independence, C 290, s_max ∈ {0.1, 0.5, 1, 5} mm; 8 × 2e5 histories).** Shards 29–64:
`RUN-20261111T144927Z-09cf2fc5` … `RUN-20261111T213722Z-420ae1d7` (36 runs, listed in `v8-lvhr/index.json`, steps 29–64);
combine 65 `RUN-20261111T214033Z-180cba48`. Exit θ_rms of the 0.5 R slab: 7.915 / 7.887 / 7.905 / 7.903 mrad at
s_max 0.1 / 0.5 / 1 / 5 mm against the I10 quadrature 7.898 mrad: largest deviation 0.22 % (tolerance 0.5 %), pairwise
spread 0.36 % (tolerance 0.5 %); σ at 0.9 R 0.8910 / 0.8949 / 0.8940 / 0.8932 mm, pairwise spread 0.44 % (tolerance
1 %). **PASS.** The step-57 host call `RUN-20261111T183304Z-bd523546` was interrupted by a WSL restart (no archived
record, no published output, nothing seen) and the step was re-launched as `RUN-20261111T205654Z-ecf93a1e` (same seed);
the interrupted call is preserved in the index as a run without record.

**I14 (IDD step independence, C 290, s_max ∈ {0.1, 0.5, 1} mm × f_E ∈ {0.005, 0.02}; 6 configurations × 2e5
histories, 48 shards).** Shards 66–113: `RUN-20261111T214043Z-4e67621f` … `RUN-20261112T224006Z-c165b084` (steps 66–113
in the index); combine 114 `RUN-20261112T225438Z-a8fa755d`. Reference s_max 1 mm / f_E 0.02; max pairwise |ΔR80| =
0.0037 mm (tolerance 0.05 mm; R80 162.953–162.957 mm over the six configurations); grouped IDD comparisons against the
reference (165 bins, studentized permutation, Bonferroni bound 4.52, no merged or unsupported bin): p = 0.761 / 0.960 /
0.946 / 0.799 / 0.240 and max |z| = 2.85 / 2.27 / 2.54 / 2.89 / 3.36 for s01-f0005 / s01-f002 / s05-f0005 / s05-f002 /
s1-f0005 (every p > 0.001, every max |z| below the bound, amendment ¹⁷). **PASS.** The step-70 host call
`RUN-20261111T222703Z-13954013` completed on the host (signed record: exit 0, `PARTIAL i14-s01-f0005-s4.json
030ec188…`) but, after an unexpected WSL restart, its published outputs were not visible from the sandbox; the step
was re-launched as `RUN-20261112T140322Z-964c767c`, which produced the byte-identical partial (same digest; same seed,
deterministic single-process transport) and is the run that was staged and imported. Both records are committed
(roles `evidence` and `superseded_duplicate`); the identity of the two partials is asserted by the inventory test.

**I15 (voxel-boundary scattering: transport voxels 1 and 2 mm, unshifted and half-shifted, C 290; θ 4 × 2 and σ
4 × (4, 4, 2, 2) shards × 2e5 = 1.6e6 histories).** Shards 115–134: `RUN-20261112T225448Z-24e2fcdc` …
`RUN-20261113T014227Z-3012a901` (steps 115–134 in the index); combine 135 `RUN-20261113T015413Z-8c57c225`. θ_rms
7.909 / 7.898 / 7.915 / 7.892 mrad for v1-sh0 / v1-sh1 / v2-sh0 / v2-sh1, pairwise spread 0.29 % (tolerance 0.5 %),
θ/quadrature 0.9992–1.0021; σ at 0.5 R / 0.9 R 0.3322 / 0.8967, 0.3299 / 0.8932, 0.3322 / 0.8971, 0.3291 / 0.8916 mm,
pairwise spread 0.97 % (0.5 R) and 0.61 % (0.9 R) (tolerance 1 %), σ/Fermi–Eyges 0.9988–1.0084 (informative). **PASS.**

**T13-ion (hr8 step 01, CUDA chunk invariance 2¹⁰ vs 2¹⁸, C 290, 1e5 histories; `RUN-20261113T015455Z-790a1181`).**
Energy deposition (1 602 000 voxels), tallies and every diagnostic counter bit-identical between the chunk sizes in
float32 and in float64; all counters 0. **PASS.** (An earlier launch of this step, `RUN-20261113T015423Z-d1d4986b`, was
rejected by `run_suite.py` before any step ran because the orchestrator passed `--only 1` instead of the zero-padded
`01`; exit 1, no output; the record is committed with role `rejected_launch`.)

**T12-ion HR (hr8 steps 02–11 shards, 12 combine `RUN-20261113T033810Z-176a1c8f`; 4.6e6 histories: c290-cpu64 2 × 2e5
(`RUN-20261113T020456Z-05ed7790`, `RUN-20261113T021924Z-1eaa5433`), c290-cuda32 1e6 (`RUN-20261113T023349Z-79c64988`),
c290-cuda64 1e6 (`RUN-20261113T023615Z-b875aca3`), he150-cpu64 2 × 2e5 (`RUN-20261113T023844Z-91c97e14`,
`RUN-20261113T024934Z-685d4820`), he150-cuda32 1e6 (`RUN-20261113T030022Z-482dcc1b`), o350-cpu64 2 × 2e5
(`RUN-20261113T030211Z-1ac4b788`, `RUN-20261113T031845Z-8a064208`), o350-cuda32 1e6 (`RUN-20261113T033518Z-2e8ccde1`)).**
All four frozen pairs pass every T12 rule (arrays IDD and lateral profiles at 0.25 / 0.5 / 0.9 R with the grouped χ²,
max-statistic and Bonferroni rules, no merged or unsupported bin; scalars R80 and σ_lat(0.5 R) by the z rule, bound
3.5; total deposit by the deterministic precision bound 1e-5): c290-cpu64 vs c290-cuda32 (studentized bootstrap-t,
40 / 100 batches): IDD p 0.373, max |z| 2.67 (bound 4.52); lat25 / lat50 / lat90 p 0.210 / 0.501 / 0.086; R80 z 0.36;
σ_lat z 0.08; total deposit relative difference 7.6e-10. c290-cuda32 vs c290-cuda64 (f32-vs-f64 probe, permutation):
IDD p 0.100, max |z| 3.02; lat p 0.759 / 0.289 / 0.719; R80 z −1.29; σ_lat z −2.05; deposit 7.6e-10. he150-cpu64 vs
he150-cuda32: IDD p 0.285, max |z| 3.30; lat p 0.504 / 0.080 / 0.459; R80 z −1.00; σ_lat z −1.04; deposit 4.2e-10.
o350-cpu64 vs o350-cuda32: IDD p 0.665, max |z| 2.25; lat p 0.102 / 0.342 / 0.732; R80 z −1.63; σ_lat z −1.29; deposit
8.7e-10. **PASS.** (The CUDA samples ran 100–155 s per 1e6 histories; the single-process warp-cpu shards 640–990 s per
2e5.)

**Findings (no re-grading).** (i) The I9 re-run closes the attribution of item 4k by measurement: the factor 1.734 was
the helper constant, the transport's lateral variance follows the Fermi–Eyges prediction of its own T_ion to ≤ 0.5 %.
(ii) The step-independence and voxel rows (I12, I14, I15) show no resolvable dependence of range, IDD shape, exit angle or
lateral width on the step controls or the voxel grid within their tolerances. (iii) Backend parity holds for He, C and
O between warp-cpu float64 and warp-cuda float32 and between CUDA precisions; T13-ion is bit-identical. (iv) The three
operational events of the block (two WSL restarts, one rejected launch) left no evidence gap: one interrupted call has
no record and was re-run, one completed call was superseded by a byte-identical re-run, one launch was rejected before
running; all are preserved in the committed index. **Consequence (plan D / F).** The LV/HR verification gate of
section A is closed. Remaining before the integration decision: V5-SS, the D-4 record, and the CI / review / validation
gates of this results commit. No clinical claim.

**Binding of the combine runs to their consumed partials (2026-11-13, orchestrator; Codex REVIEW-195f5f58 on a326cee3,
finding 1 "important").** The review found that the committed combined summaries were compared with the signed shard records
only through the unsigned summary itself, while the combine's signed stdout names only a mutable cache path. Closed without
re-running anything: the published outputs of the five combine runs, the T13-ion run and the preserved e602a583 I9 combine
(step file with the step document, `summary.json`, `manifest.txt`, `environment.txt`) are committed verbatim under
`validation/results/transport/v8-lvhr/published-outputs/`; the signed `output_manifest_sha256` of each record reproduces from
them, and the combine step document inside carries the sha256 of the partials manifest the step read and the content digest
it recomputed for every imported partial (steps_v8 combine fails closed on a digest mismatch) with the shard's host run id.
`tests/ionmc/test_v3_008_lvhr_host_records.py` now requires: signed manifest == reproduction; the signed step document's
partial set == the committed summary's set for that step; every digest == the `PARTIAL` line of the named shard record; the
T13 document's identity flags; the preserved I9 FAIL document's binding at e602a583. No grade, number or run changes.

### V5-SS execution note (2026-11-13, orchestrator; prospective; before any V5-SS ionmc sample or TOPAS run)

Executes row V5-SS (branch B) exactly as pre-registered (dated notes of 2026-11-01, rounds V3-008-4/-5; cost decision
2026-11-09: 6e4 TOPAS primaries per case, fallback not invoked). Nothing else of the pre-registered design changes: cases,
gated pairs, tolerance 3 %, TOST rule, σ > 1 % = INCONCLUSIVE, N (TOPAS 6e4 in parts 0–1, ionmc 4e5), estimator, lineage
wording. This note only fixes what the frozen text leaves open (D1–D6).

**Order of operations.** (i) One reviewed implementation commit (generator + bundles + SS manifest, sampler step, evaluator,
tests, this note, the section-C seed note); `git diff 712f484a -- src pyproject.toml uv.lock` empty (implementation of the
LV/HR block). (ii) Eight ionmc samples on the controlled host runner (one call per case, warp-cpu float64, single process,
4e5 histories, `--expected-sha` = HEAD of (i)). (iii) SS seal commit (reviewed, pushed): `ss_evaluate.py seal` →
`validation/results/transport/v8-ss/ss-sample-seal.json` + the signed host-run records, then the dated SS seal note with the
evaluator binding block. (iv) Scorer smoke (after (iii)). (v) 16 TOPAS evidence parts, one at a time, ≤ 1200 s planned each.
(vi) Records commit; `ss_evaluate.py check-run` on all 16. (vii) `ss_evaluate.py grade` once →
`validation/results/transport/v8-ss/ss-report.json`; dated results note; never re-graded.

**SS-SEAL-1 (mirror of R2-TOPAS-4).** A V5-SS TOPAS run (smoke or evidence) is valid only if (a) the `code_sha` of its run
record has the SS seal commit as an ancestor (`git merge-base --is-ancestor <ss_seal> <code_sha>`), and (b) the SS seal
commit's independent review (review id in the SS seal note) finished before the run's controlled `reference_started` event
(controlled-service records, never Git timestamps). `ss_evaluate.py` enforces both and fails closed when a record is missing;
a run failing either check invalidates the pairs that use its case (FAIL of the blinding control, never re-graded).

**Decisions.** D1 ionmc seeds ("ionmc seeds from the lv8 row block (r_index declared in `steps_v8.py`)"): base 20631004
(`RERUN_BASE_2`, declared, never drawn; lv8 20581004 is consumed), row `v5ss`, r_index = case order 0–7 (22, 25, 26, 27, 28,
29, 30, 31), k = chunk index → 20631004 + 100·r + k (`steps_v8.SS_SEED_BASE`, `SS_R_INDEX`); section-C note of 2026-11-13.
D2 window ("±2.5 θ window anchored on the ionmc prediction (fixed before the run)"): θ_anchor,c = untruncated pooled (x, y,
about zero) projected rms of the sealed ionmc sample; W_c = `ss_budget.WINDOW_K` · θ_anchor,c = 2.5 θ_anchor,c; it is ionmc's
own prediction, deterministic and sealed before any TOPAS output (P10 §7.5.1: a ±5 % anchor shift moves the pair biases ≤ 0.07
pp). D3 generator: new `validation/scripts/ion_em/gen_ss_bundles.py` with its own manifest
`validation/reference_cases/topas/v5ss-frozen-manifest.json`; the r2 generator and r2 manifest stay byte-identical (the r2
evaluator requires the r2 manifest at `code_sha` to equal HEAD's). D4 scorer smoke: the phase-space parameter names and the
ASCII layout are unverified, so one smoke `v5ss-helium-water-250mevu-scorer-smoke` (He 250, 200 primaries, rehearsal seed
**20600438**; **20600439** reserved for one re-smoke if a parameter name is wrong: regenerated bundles, new SS manifest,
reviewed commit) runs after the seal; its angles are never used. D5 bootstrap ("σ = bootstrap SE of the ratio statistic over
primaries and parts"): B = 2000; TOPAS resampled with replacement within each part, ionmc over histories, each (case, side)
independently; one dev replicate per b; σ = SD; analysis seed `ANALYSIS_SEED = 20261113` (numpy, not a transport seed).
D6 ionmc physics ("ionmc: same slab, warp-cpu f64"): the default all-EM configuration of the TD/TL samples
(`steps_v8.sample_config`: MCS on, straggling on, nuclear off, Bethe, `ls = "off"`), 40 batches, pencil beam at z = −0.001 mm,
`BoxPhantom((−60, −60, 0), (120, 120, 1.0), WATER)`, escape records.

**Statistic (both sides).** θ_w(W) = sqrt( Σ_i [θx_i² 1(|θx_i| ≤ W) + θy_i² 1(|θy_i| ≤ W)] / Σ_i [1(|θx_i| ≤ W) + 1(|θy_i| ≤
W)] ); TOPAS θ = atan2(u, w), w = +sqrt(1 − u² − v²) from the phase space; ionmc θ = atan2(u_x, u_z). Pair A/B:
dev = (θ_model,A/θ_model,B)/(θ_ref,A/θ_ref,B) − 1; verdict `ss_budget.gate(100·dev, 100·σ)`; row PASS iff all five PASS.

**Files.** Bundles `validation/reference_cases/topas/v5ss-<ion>-water-<E>mevu-case<NN>-part<P>/`; sampler `steps_v8.py
ss-sample` (npz `ionmc-ss-case<NN>-k<K>.npz`; wrapper `ss_wrapper_template.py`, argv `["python", "-c",
render_wrapper(case, sha, "v5ss-samples", "1")]`); evaluator `validation/scripts/ion_em/ss_evaluate.py`; records under
`validation/results/transport/v8-ss/` (`host-records/`, `ss-seal-review.json`, `topas-runs/<run_id>/`, `smoke-runs/`).

| case | ion | MeV/u | BeamEnergy | TOPAS seeds part 0 / 1 | ionmc seed (k = 0) |
|---|---|---|---|---|---|
| 22 | He (`alpha`) | 150 | 600 MeV | 20619022 / 20619122 | 20631004 |
| 25 | C (`GenericIon(6,12)`) | 150 | 1800 MeV | 20619025 / 20619125 | 20631104 |
| 26 | O (`GenericIon(8,16)`) | 150 | 2400 MeV | 20619026 / 20619126 | 20631204 |
| 27 | C | 250 | 3000 MeV | 20619027 / 20619127 | 20631304 |
| 28 | He (`alpha`) | 250 | 1000 MeV | 20619028 / 20619128 | 20631404 |
| 29 | C | 430 | 5160 MeV | 20619029 / 20619129 | 20631504 |
| 30 | O | 430 | 6880 MeV | 20619030 / 20619130 | 20631604 |
| 31 | O | 250 | 4000 MeV | 20619031 / 20619131 | 20631704 |

30 000 primaries per part. **Expected cost** (`ss_budget.cost_from_measured(0.023 s, 6e4)`, upper estimate): REF s per case
(per part) He 150 528 (264), C 150 1740 (870), O 150 2023 (1012), C 250 1454 (727), He 250 378 (189), C 430 1209 (604),
O 430 1514 (757), O 250 1752 (876); total 10 596 s = 2.94 h, every part ≤ 1200 s; ionmc ≈ 1–10 min per call (inferred).

**Dated SS seal note (2026-11-13, orchestrator; written after the independent review of the SS seal commit finished and
BEFORE any V5-SS TOPAS run, smoke or evidence).** The SS seal commit is `7fe6ab1d7b7eb18ab377156e0551057264fba40b` (seal record
`validation/results/transport/v8-ss/ss-sample-seal.json`, eight ionmc thin-slab samples at implementation SHA bef32873, seed base
20631004, windows W = 2.5 θ_anchor fixed; the eight signed host-run records committed; validation record VAL-20261113-061649-BF6904).
Its review REVIEW-5d9e8df442e94669aef5fe609794b608 passed (finished 2026-11-13T06:10:46.335228Z; record
`validation/results/transport/v8-ss/ss-seal-review.json`, copied verbatim from the review service). Binding rule SS-SEAL-1: every
V5-SS TOPAS run must start after that time at a `code_sha` descending from the seal commit; `ss_evaluate.py check-run` / `grade`
read the block below and fail closed otherwise. The scorer smoke (seed 20600438) and the 16 evidence parts follow.

```json
{"ss_seal_sha": "7fe6ab1d7b7eb18ab377156e0551057264fba40b", "ss_seal_review_id": "REVIEW-5d9e8df442e94669aef5fe609794b608"}
```

**Dated note (2026-11-13, orchestrator) — scorer smoke result and bundle regeneration (pre-declared route of the execution note).**
The scorer smoke `v5ss-helium-water-250mevu-scorer-smoke` (seed 20600438) ran at `code_sha` 677e7c27 after the SS seal review
(REF-740144c3826978e4ea78-ef394dd8, 2026-11-13 06:34:06–06:34:08 UTC, exit 0, 1.4 s): every PhaseSpace parameter was accepted;
the `.header` column list equals the evaluator's expected ten columns (X, Y, Z [cm], cosX, cosY, E [MeV], weight, PDG, the
negative-w flag, the first-scored flag); 200 of 200 primaries reached the surface and were scored once (PDG 1000020040, Z = 0.1 cm
in world coordinates, flags 0/1, cosines printed with six significant digits); `CoulombScat` for alpha and GenericIon with
`IonCoulombScattering` only and no `msc` block for any particle. Two log checks of `ss_evaluate.py` could not be satisfied by
this bundle: the log prints no production-cut couple table (so the e⁻ threshold ≥ 1.16 MeV cannot be read), and the global
"Registered Physics Processes" table lists an `msc` entry (for a particle other than the beam particle; no per-particle msc
block exists). Resolution, both execution/printout-level: (i) every V5-SS bundle gains `i:Ts/RunVerbosity = 2` (the setting of the
committed r2 verbose smoke, which prints the couple table with the e⁻ threshold); the 17 bundles and the SS manifest are
regenerated (seeds, physics, geometry, scorer unchanged), and the scorer smoke is re-run at the reserved seed **20600439**;
(ii) the evaluator's process-list rule becomes per beam particle (fail on an `msc: for <beam particle>` block, on a missing
`CoulombScat`, or on a model other than `IonCoulombScattering`; a global `msc` entry for another particle is recorded, not a
failure). The angles of the smokes are never used. No evidence run has started; the 16 evidence parts run at the regenerated
commit (a descendant of the seal commit, after its review).

**Dated note (2026-11-13, orchestrator) — re-smoke 20600439 and the couple-table rule.** The re-smoke
`v5ss-helium-water-250mevu-scorer-smoke` (seed 20600439) ran at `code_sha` ea10ce6a (REF-8358f3c4f2f32a8453f4-a7dbae17,
06:54:09–06:54:11 UTC, exit 0): identical scorer result (ten expected columns, 200/200 primaries scored once, alpha only,
Z = 0.1 cm, six-digit cosines), `CoulombScat`/`IonCoulombScattering` for alpha, no msc block for any particle, and the
production-cut couple table now printed: couple 1, `G4_WATER` (the slab), `Range cuts … e- 1 km`, `Energy thresholds … e-
9.54965 GeV` — above the required 1.16 MeV, so no δ electron can be produced in the slab; couple 0 is the vacuum world
(e⁻ threshold 990 eV, the floor; nothing interacts in vacuum). `ss_evaluate.py` had taken the minimum threshold over all
couples and therefore failed on the vacuum couple; the rule is corrected to read the couple of the slab material `G4_WATER`
(every other couple is recorded, not graded). Evaluator change only (reviewed commit); bundles, seeds and the seal are
unchanged; no evidence run has started.

**Dated note (2026-11-13, orchestrator; Codex REVIEW-8252be27 on 4505c89f, finding "important").** The comment rewording of
4505c89f also regenerated the scorer-smoke bundle after its re-smoke had run (REF-8358f3c4f2f32a8453f4-a7dbae17 at
ea10ce6a), so the evaluator's exact binding (the smoke's manifest entry at the run's `code_sha` must equal HEAD's) would
have rejected that smoke. Corrected without any new run: the generator now keeps the smoke bundle byte-identical to its
run-time version (input.txt sha256 eb9a4694…, legacy comment wording; a documented special case pinned by test), while
the 16 evidence bundles, which have not run, carry the clarified wording. The smoke of seed 20600439 therefore stays
the pre-registered run condition evidence; no third smoke seed is used.
