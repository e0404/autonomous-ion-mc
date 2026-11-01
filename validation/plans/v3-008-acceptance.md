# V3-008 acceptance criteria (DRAFT r5 — to be frozen before any result)

> **Status: DRAFT r5 (2026-11-01), not yet frozen.** Acceptance criteria are frozen BEFORE any result: this draft becomes revision 1 (frozen) by a reviewed commit after the pre-freeze checks of §P (steps 1–7, incl. the new P10 and P11) are executed; they produce no qualification result (P4 smoke runs carry no scorers). Independent review: rounds V3-008-1 and V3-008-2 (`validation/plans/v3-008/review.md`). Seed bases per the shared ledger (V3-005B commit 06974d07): V3-008 reservation 20572004–20632003.


**Status.** DRAFT r5 (planner, 2026-11-01). Revision history:
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
- Nothing in §P has produced a transport or TOPAS result.
- It is committed with decision 0043 as the first scientific commit of task V3-008, after
  independent review and after every pre-freeze check of section P is closed.
- It is committed BEFORE any V3-008 transport, LV, HR or reference-engine result exists.
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

Closed: P8, P7, P2 (source reading). P2b is very probably inapplicable. **Order before freeze (r5):**
1. P2 follow-ups, P4 and P10 (scorer-free smokes only).
2. P11: the exact-Ω specification, with the I11-map and I11-rel tests committed.
3. P3.
4. P5/P5b.
5. P1 completion.
6. P6 (measured C/O range; lateral-data search).
7. P9.

| # | Check | Outcome feeds |
|---|---|---|
| P8 | Repository facts. **CLOSED** (`validation/plans/v3-008/p8-repository-facts.md`): `A16_BASELINE` = f3a1dd62 (full SHA, fail-closed resolution); ASTAR object dd4a2e50… spans 1e-3 to 1e3 MeV total, i.e. ≤ 250 MeV/u; species ids 5 and 6 are free on this branch (PR #82 to be checked at integration); proton gates at `config.py:653–657` and `transport/tables.py:200–205` | R2, I2, I0 |
| P7 | Seed bases. **CLOSED**: ledger on the V3-005B branch (commit 06974d07, `V3_008_RESERVED = range(20572004, 20632004)`); guard test `tests/ionmc/test_v3_008_seed_families.py` committed | C |
| P2 | Geant4 v11.4.2 ion models from source. **Source reading CLOSED; follow-ups OPEN and required before freeze (order step 1)** (`p2-geant4-ion-models.md`; [URB] 22c5f52b…, [LS] e16ed6cd…, [OPT4] 04337300…). FACT: GenericIon = `G4ionIonisation` + `G4LindhardSorensenIonModel`; switch at 2.0 MeV scaled; above it the `G4BetheBlochModel` table + q_eff²(ΔL_LS + 2·Barkas·(z−1)/z), with no separate Bloch/Mott; ICRU 73/90 ion data only ≤ 2.0 MeV scaled and 2 < Z ≤ 80, never alpha; ion msc = `G4hMultipleScattering("ionmsc")`, model not set (Urban inferred); Urban θ0 = 13.6·\|charge\|·√y·invbetacp·(c1 + c2 ln y), y = t/X0, with no charge/β in the log (frozen with line references). The r3 lineage trigger ("ion tables above 2.5 MeV/u") is NOT met. **Follow-ups, required before freeze:** read `G4EmBuilder.cc` (40782f5c…) and `G4hMultipleScattering.cc` (a6643869…), both materialized in the v3-008 cache, and the `G4VMscModel` defaults (Urban for ions; WentzelVI for protons; alpha assignment); pin `G4BetheBlochModel.cc/.hh` (the ICRU 90 proton/alpha array claims stay "research note, partly confirmed" until then); read `G4IonICRU73Data.cc` (cf771a5a…) and `G4IonisParamMat.cc` (990c6f31…; Zeff of G4_WATER: θ0 differs by 2.5 % between 10/3 and 6.6) and X0 of G4_WATER; λ scaling and lateral displacement for ions (`G4VMscModel`, `G4EmTableUtil`) | §B lineage, V2 branch |
| P2b | **As defined in r3, very probably INAPPLICABLE**: opt4 protons use WentzelVI (INFERRED strong), so the archived T15 does not test an Urban method, and the r3 Fermi–Eyges δ_pred is not a well-defined Urban prediction (truncated core + tail; per-length width step dependent, factor 0.756 between 0.1 and 1 mm; TOPAS ion steps set by the bundle). Consequence (r5): V2 follows the pre-declared branches. **P2b′** (only if a later amendment reintroduces a graded Urban route): the Urban-emulator method is checked against a TOPAS configuration in which protons use Urban (source-verified constructor, e.g. opt3 if `ConstructCharged(..., false)`, TO VERIFY), same geometry and declared step model, run before any ion TOPAS run; e_δ ≤ 2 % is required for the extension | V2 |
| P4 | TOPAS settings, checked on smoke runs (execution only, not graded): `alpha` vs `GenericIon(2,4)` particle naming; `GenericIon(8,16)`; BeamEnergy = total kinetic energy; acceptance of a very large electron production cut so that deltas deposit locally (e.g. `d:Ph/Default/CutForElectron = 1 km`; exact parameter TO VERIFY); run-time printout confirming the LS ion model and the ICRU90 flag (`/process/em/printParameters` or verbose); throughput per 1e5 primaries P4 smoke cases carry no dose, IDD or phase-space scorers (or their outputs are deleted unread by the runner); only the run log, EM-model printout and wall time are recorded. The measured throughput sizes every REF case into parts of ≤ 1200 s planned (section E)  **r5: required before freeze (order step 1).** The scorer-free smokes archive the EM-table / `printParameters` printouts: msc model, step-limit type, RangeFactor and lateral displacement for GenericIon, alpha and proton; ionisation model names; the G4_WATER material table (RadL) | V1, V2 branch, V3, host budget |
| P10 | **Single-scattering availability (order step 1, with P4):** whether OpenTOPAS 4.3 / Geant4 11.4.2 offers a single-scattering configuration for GenericIon and alpha (`G4CoulombScattering`; a g4em-standard_SS-type module), its run cost for the V5-SS cases (from a scorer-free smoke), and, if available within 4 h REF, the pre-registration of V5-SS (observable, window, histories, seeds, final tolerance) **before any run** | V5-SS, lateral gap |
| P11 | **Exact-Ω specification (order step 2):** commit the b-matching mapping (b = ln Ω + Δb_Fano at the respective velocities; 2-D dependence on (pv, p1v1); p1v1_eq/pv_eq = p1v1/pv), the evaluation choice (per-step root solve, or a (pv, p1v1) tabulation with a test bounding the θ interpolation error ≤ 0.1 %), and the I11-map and I11-rel tests, **before the mapping is evaluated against the Molière–Hanson reference**; hashes into P9 | §3, I11-map, I11-rel |
| P3 | Committed numpy LS script: ΔL_LS − (Bloch + first-order Mott) for z = 6, 8, water, 2–430 MeV/u, tested against the Bloch limit at small z. Expected ≤ 0.1 % of S (estimate). If it exceeds 0.2 % of the CSDA range at any domain energy, an LS option is added (amendment). **OPEN** | decision §1, I3b |
| P5 | Data inspection of G4EMLOW 8.8 (`https://geant4-data.web.cern.ch/datasets/G4EMLOW.8.8.tar.gz`, 350 393 595 B, MD5 328330009df633f7e9b3a9f445745298), members `ion_stopping_data/icru73/z{6,8}_G4_WATER.dat` and `icru90/z6_G4_WATER.dat`. Checks: (i) energy span covers ≥ 2–430 MeV/u; (ii) units and per-nucleon axis; (iii) erratum provenance. The original ICRU 73 water used I = 67.2 eV and the erratum uses 78 eV; decide by the high-energy logarithmic slope against the Bethe L0 at two I values (aggregate verdict only, no values in Git) and by documentation; (iv) origin of `icru90/z6`. **P5b:** locate libdEdx (GPL) and, if it ships ICRU 73 (rev.) water files, use it as a cross-check or substitute | I3 |
| P1 | Molière–Hanson reference. **PARTLY CLOSED** (`p1-moliere-hanson-reference.md`, P1 and P1b). Done: the reference with the Fano correction is implemented; proton validation passes Be (0.10 %), Al (0.02 %) and Cu (0.21 %); **the ≤ 1.0 % proton validation FAILED for Pb** at x/R1 ≥ 0.9 (+1.10 %, +1.32 %; `xfail(strict=True)`, tolerance not loosened; recorded as FAILED and preserved). Findings carried into decision 0043 §3: water is not validated directly; the compound Fano rule is derived, not sourced; proton f_dM sits +2.5 to +3.5 % above the water reference; c_Ω ≈ 0.09–0.11 (r3's 0.13 was wrong); the r3 selection rule selects nothing (seen tables disclosed in §3). r4/r5 replace the selection rule by the exact-Ω mapping (the exact version of r3's first-ranked candidate; the selection rule was revised after the r3 rule reached STOP 2 on seen P1/P1b numbers — chronology and disclosure in decision 0043 §3) and demote the absolute reference comparison to the informative row I11-abs, with I11-map and I11-rel gating. **P1 completion (order step 5):** Pb attribution on ICRU 49 ranges; a water/compound θ_Hanson check if a source exists (Gottschalk polystyrene/Lexan/water values); if Pb is not cleared, P1 stays FAILED for Pb, as recorded. The reference-uncertainty statement goes into decision 0043. The P1 driver evaluates the exact I11 grid. Codex REVIEW-98315d4e: n_Ω is the full derivative of the proton Ω (or b) including the screening and Fano β-dependence; this matters for the superseded linearised (d) only. The Kanematsu ion form is moot | §3, I11, I11-rel, I11-map |
| P6 | Measured-range candidate: Sihver, Schardt, Kanai, Jpn. J. Med. Phys. 18 (1998) 1 (J-STAGE). Check open access and licence, whether Bragg-peak/R80 positions for C and O in water are tabulated with uncertainties, and whether beam-line water-equivalent offsets are given. Bounded investigation: two access approaches, ≤ 2 h. **r4: also search for measured ion lateral widths in water** (Parodi et al. 2013 J. Radiat. Res. (HIT); Tessonnier 2017; Schwaab 2011), numeric availability TO VERIFY; a negative result is recorded as the MCS evidence gap of STOP 4 | V3, MCS gap (the lateral-gap intervention is used only after P6 and P10 both fail) |
| P9 | Freeze record (D-E8 pattern). It binds the committed sha256 listed in the P1 note (`validation/scripts/ion_em/__init__.py` 7962c4af…, `validation/scripts/ion_em/moliere_hanson.py` 5a99b3b8…, `tests/ionmc/test_ion_em_moliere_hanson.py` 789fd4b1…; the hashes are regenerated from the final committed files at every revision of this draft and bound verbatim at the freeze), the P11 mapping script and the I11-map/I11-rel tests, the P2 follow-up source hashes (40782f5c…, a6643869…, 990c6f31…, cf771a5a…), the P2 source hashes ([URB], [LS], [OPT4], [ICRU73]) and the P3 script and output. The freeze commit must have those commits as ancestors, the recorded hashes must match, and the tree must be clean | freeze |

## A. Verification rows

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| I0 | **Capability contract.** He/C/O sources run on python, warp-cpu and warp-cuda. Fail closed for: He-3 and any other ion; energy outside the domain or E0 + 6σ_E above it; `nuclear=True` with z ≥ 2; `nist-star` for C/O; ASTAR with a non-water material; LET medium ≠ water; unproducible species/generation tallies. Requested/effective configuration persists projectile, z, mass, `e_cut_ion_mev_per_u` and the MCS model id | every case raises `UnsupportedCombinationError` (or `valid=False` plus a named counter); fields present | capability | CI |
| I1 | Shared functions with a charge argument (z ∈ {1, 2, 6, 8}), the velocity-scaled arguments and the exact-Ω (b-matching) mapping functions (root solve or bounded tabulation). U1 harness: twin vs Warp CPU f32/f64; the completeness test covers every new callable | decision 0001 classes as frozen in U1 ¹⁸ | backend parity | CI |
| I2 | **He table reproduction (SEEN, regression).** Bethe 4He at I = 75 eV vs ASTAR (75 eV), and Bethe at 78 eV vs ICRU 90 alpha, water (`compare_nist.py` extended). The I-value difference between ASTAR (75 eV) and the default model (78 eV) is stated in the result | S: max abs rel dev ≤ 0.30 % for e ≥ 10 MeV/u, ≤ 1.0 % for 2–10 MeV/u. CSDA range at 50/100/150/200/250 MeV/u: ≤ 0.20 % (ASTAR) and ≤ 0.25 % (integrated ICRU 90). Aggregates only in Git. Bounds = 1.5× the SEEN decision-0038 maxima, rounded up; derived from seen values, so I2 is a regression guard, not evidence | ion-specific tabulated (regression) | LV |
| I3 | **C/O table comparison.** Bethe (78 eV) vs ICRU 73 (rev.) water for C and O, and vs ICRU 90 for C if P5(iv) establishes its origin. Source per P5/P5b | S, e ≥ 10 MeV/u: max abs rel dev ≤ 1.5 % and RMS ≤ 0.8 %. S, 2–10 MeV/u: ≤ 6 % (effective-charge model spread plus PASS low-energy behaviour). CSDA range at the ends and standard points ≤ 0.5 %. Defended from theory: omitted LS ≤ 0.1 % (P3); Barkas/shell parameterisation ≤ 0.5 % (alpha experience 0.2 % SEEN, scaled for z³ terms); reference uncertainty ≈ 1–2 % (TO VERIFY). **If P5 and P5b both fail:** row = EVIDENCE GAP, intervention per protocol; it cannot pass. O has ICRU 73 only | ion-specific tabulated (largely independent of TOPAS above the confirmed 2.0 MeV scaled-energy switch (≈ 1.98 MeV/u; P2)) | LV |
| I3b | Informative: (i) Bethe with each correction off in turn (Barkas, Bloch, Mott, effective charge) vs the I3 reference, range shift per term in mm; (ii) the P3 LS-minus-(Bloch + Mott₁) curve per energy | report only | related model | LV |
| I3c | Exploratory attribution diagnostic: "TOPAS-like" C/O water stopping rebuilt from the registered ICRU 90 proton array × q_eff² + LS ΔL, vs ionmc Bethe; S and CSDA range | report only. Used only to attribute a V1 failure; never a pass criterion; never construction | exploratory | LV |
| I4 | **Velocity-scaling identity (theory).** With shell/Barkas/Bloch/Mott/effective charge off: S_ion(β) = z² S_p(β) and R_ion(β) = (M_ion/m_p)/z² · R_p(β) at every grid node 2–500 MeV/u, water and tissue, with the Tmax mass term evaluated exactly | rel ≤ 1e-6 (≤ 2e-3 without the Tmax correction, reported) | theory | CI |
| I4b | Range-straggling scaling, numpy quadrature σ_R² = ∫(dσ²/dx)/S³ dE with the implemented Bohr variance, corrections off (exact algebra) | σ_R/R(ion) / [σ_R/R(p) · (m_p/M)^½] − 1 within 1e-3 at the standard points | theory | CI |
| I5 | Range-table consistency for A > 1 (V3-003D rows D1–D3 re-pointed to ion tables) | round trip ≤ 1e-5; strictly monotone; D2 distance ≤ 0.5 µm | self-consistency | CI |
| I6 | Deterministic CSDA end depth, straggling and MCS off, at the standard points and ends (CI: He 100 MeV/u only). The residual range at E_cut is reported per ion against the decision's estimates (He ≈ 0.075, C ≈ 0.025–0.03, O ≈ 0.02–0.025 mm) | R(E0) − R(E_cut) − 1e-4 R ≤ z_end ≤ R(E0) + 1e-4 R | self-consistency | CI/LV |
| I7 | Energy balance with ion sources, all EM on | rel ≤ 1e-12 (f64, python), 1e-5 (f32); counters 0 | conservation | CI |
| I8 | Range straggling in transport: σ of the end depth (straggling on, MCS off), 1e5 histories warp-cpu per standard point | within 5 % of the I4b quadrature with corrections on | self-consistency | LV |
| I8b | Report-only: σ_R/R at the standard points vs σ_R/R(p at the same range)/√(M/m_p) and vs the literature statement σ_R/R(C) ≈ 0.3–0.4 % (Schardt et al. RMP 2010, TO VERIFY) | report; outside ±15 % is recorded as a finding (not a failure) | related model | LV |
| I9 | Lateral variance of deposits (MCS on, straggling off): 1 mm slabs, 0.1 mm lateral bins, Sheppard-corrected, z/R ∈ {0.25, 0.5, 0.75, 0.9}, standard points, 1e6 histories warp-cpu | within 2 % of Fermi–Eyges A2(z) with the implemented T_ion | self-consistency | LV |
| I10 | MCS integrator step independence (U5 analogue) for T_ion: s ∈ {0.01, 0.1, 0.5, 1, 2, 5} mm, x/R ∈ {0.05, 0.25, 0.5, 0.9}, standard points; negative control: per-step Highland | \|Δθ²/θ²\| ≤ 2e-3; the control varies ≥ 5 % | self-consistency | CI |
| I11-map | **Exact-Ω mapping implementation (deterministic):** at every grid node, water and the `materials.py` set, b_p(pv_eq, p1v1_eq) = b_ion(pv, p1v1) (b = ln Ω + Δb_Fano at the respective velocities; p1v1_eq/pv_eq = p1v1/pv) is satisfied on a (pv, p1v1) test set covering the domain, by the per-step root solve, or by the tabulation if chosen (then the θ interpolation error is ≤ 0.1 %); the python twin and the Warp row agree; z = 1, M = m_p gives pv_eq = pv exactly (the proton path is unchanged) | b residual ≤ 1e-10 absolute (float64 root); twin vs Warp per decision 0001; z = 1 identity exact | self-consistency | CI |
| I11-abs | **Informative: absolute MCS model vs Molière–Hanson ion theory** (P1 reference, committed script): θ_model(x)/θ_MH(x) − 1, x/R ∈ {0.1, 0.25, 0.5, 0.9}, standard points and ends, reported **together with the reference's own uncertainty** on the same grid (f_dM's proton offset vs the water reference; the Pb residual, P1 validation FAILED for Pb unless cleared by P1 completion; water not validated; compound Fano rule derived, not sourced) | **report only.** Not a pass criterion: the water reference carries about 3 % uncertainty (decision 0043 §3) and the model carries a ±5 % theory-derived uncertainty, stated in metadata and docs | source-model consistency (informative) | LV |
| I11-rel | I11-rel (LV, gating): \|(θ_ion,model/θ_p,model) / (θ_MH,ion/θ_MH,p) − 1\| ≤ 2 % at every I11 point (equal velocity and x/R); failure is STOP (amendment), never a tolerance change. The predictability of this ratio from the seen P1/P1b decomposition is disclosed. Grid: the exact I11 grid of the P1 driver (standard points and ends × x/R ∈ {0.1, 0.25, 0.5, 0.9}) | ≤ 2 % at every point | theory consistency (relative; calibration lineage: same Molière theory as the mapping) | LV |
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
| C-res | Residual energy below E_cut deposited locally and tallied; `e_cut_ion ≥ 2 ×` table floor enforced | exact | capability | CI |

## B. Validation rows (independent references)

**TOPAS case definitions** (new bundles; template `topas/proton-water-150mev-lateral-emonly`):
- Physics: OpenTOPAS 4.3 / Geant4 11.4.2, `sv:Ph/Default/Modules = 1 "g4em-standard_opt4"`, G4_WATER (I = 78 eV).
- Geometry and beam: water box 120 × 120 × 400 mm; zero-size, zero-divergence monoenergetic pencil; BeamEnergy = total kinetic energy.
- Particles: `GenericIon(6,12)`, `GenericIon(8,16)`, alpha naming per P4.
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
- **He** stopping in TOPAS = the ICRU 90 alpha array: FACT (research note); P2 OPEN pending `G4EmBuilder.cc`/`G4BetheBlochModel.cc`. Planned as the same data as I2's reference, so the He range in V1 is NOT independent of I2. Alpha never reaches the ICRU 73 ion data (P2 FACT).
- Straggling: independent implementations. MSC: Urban (inferred) is structurally different from ionmc's model (no charge/β term in its log; truncated core + tail; step-length dependent); TOPAS lateral widths for C/O are expected 10–20 % narrower (scratch estimate, P2 §6.3, not frozen).
- The r3 pre-declared trigger (ion tables above 2.5 MeV/u; the confirmed switch is 2.0 MeV scaled, ≈ 1.98 MeV/u) is NOT met (P2): I3 keeps "largely independent of TOPAS above ≈ 2 MeV/u"; no amendment.

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| V1 | **TOPAS EM-only IDD (matched variant)** at the standard points (blind). ionmc: warp-cpu f64, 1e6 histories, all EM on | (a) \|R80_ionmc − R80_TOPAS\| ≤ 0.35 % of R80 and ≤ 1.0 mm (TOST; σ from seeds and batches). Defence: Bethe-78 vs integrated ICRU 90 proton range differs by up to 0.18 % (decision 0038, SEEN); effective-charge, LS and Barkas construction differences ≤ 0.1–0.2 % (estimate). (b) Distal width (R80 − R20) within 10 % (TOST); straggling theory spread is a few %, and Geant4 fluctuations are an independent implementation of Bohr-type theory. (c) Peak position within 0.5 mm. (d) Peak/entrance ratio within 3 % (TOST). For helium, V1 tests transport mechanics and straggling only; its range part is bounded by I2 and adds no independent stopping evidence. | independent MC (C/O: partly shared stopping lineage; He: shared stopping data) | REF + LV |
| V1-held | As V1 at **He 220, C 430, O 430, C 100 MeV/u**, run ONCE, only after V1 has been graded and after any model amendment | same criteria | independent MC (held-out) | REF + LV |
| V2 | **TOPAS EM-only lateral (matched variant):** Sheppard-corrected σ of the lateral dose profile in a ±10 mm window (window rule as T15 ¹⁰) at z/R80 ∈ {0.5, 0.9}, standard points. **Pre-declared branches, decided by the P2 follow-ups and the P4 printouts before freeze:** "If Urban: V2 informative (Case C), the raw ratio and expected direction reported." / "If WentzelVI/single scattering: V2 is decided by an amendment before any ion TOPAS run." Expected direction if Urban: TOPAS narrower (scratch magnitude 10–20 % for C/O, not frozen). The r3 Fermi–Eyges δ_pred and the r4 Urban-emulator extension are withdrawn; an extension can return only by a reviewed amendment before any ion TOPAS run | Case C: report only (raw ratio, \|z\|, direction). Any graded form only by the amendment of the second branch. V2 asserts at most consistency of the transported MCS with an independent MC, **not** an independent test of the z-dependence. V2 has no held-out analogue; a model change after V2 leaves lateral MC evidence as a declared gap | independent MC (informative if Urban) | REF + LV |
| V3 | **Measured range (C, O)** from Sihver, Schardt, Kanai 1998 (HIMAC), only if P6 passes. ionmc EM-only R80/peak position vs measured, at the published energies inside the domain. Nuclear allowance per energy: \|ΔR80(TOPAS full physics − TOPAS EM-only)\|, from additional TOPAS full-physics runs (opt4 + QGSP_BIC_HP + QMD, same geometry, 2 seeds × 1e5) | \|ΔR80\| ≤ (stated measurement uncertainty incl. beam-line WET ⊕ 0.5 % R) + nuclear allowance, at every energy used. Lineage: I_water = 78 eV calibration lineage, declared (range-consistency evidence, not an I-value test). **If P6 fails: row = EVIDENCE GAP → intervention** (inaccessible information: GSI Schardt 2007 precision Bragg curves, HIT Tessonnier 2017 / Krämer 2016 / Kurz 2012 numeric tables). He: no candidate; declared gap | measured | LV |
| V5-SS | **Conditional on P10 (registered before any run):** thin-slab projected angular width of He 150, C 290, O 350 MeV/u behind 10 mm and 50 mm of water against a Geant4 single-scattering configuration (`G4CoulombScattering` / a g4em-standard_SS-type module for GenericIon, if OpenTOPAS 4.3 offers it), core width in a ±2.5 θ window, ≥ 1e5 primaries per case in parts ≤ 1200 s | proposed: \|θ_ionmc/θ_SS − 1\| ≤ 5 % (the declared model uncertainty) and \|z\| < 3; final observable, window, histories and seeds are fixed in P10 before any run. If P10 finds it unavailable or above 4 h REF: the row is withdrawn and the lateral-gap intervention applies | independent MC (single scattering; independent of Urban and of f_dM) | REF + LV |
| V4 | Informative: ionmc vs the TOPAS 0.05 mm-cut variant (delta transport), R80 and lateral σ | report only | independent MC | LV |

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
  - Case numbers: V1/V2 matched 1–9; V1-held 10–21; V3 full physics 22–33; cut variants 34–40.
  - Each part (≤ 1200 s planned) is a separate host call with its own seed.
- **Re-runs** take fresh bases from **20621004**, then 20631004, each declared in an amendment before use. The
  last base inside the reserved block is 20631004, whose `base + 100·r + k` seeds stay ≤ 20632003 only for
  r_index ≤ 9. A re-run needing more rows, or a third re-run base, requests a new block in the shared ledger.
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
  4. Evidence gap (I3 data, V3 measured, ion MCS z-dependence): an intervention per protocol. For the lateral gap it is used only after P6 and P10 fail, with the wording of decision 0043 §3 verbatim ("Unavailable input: independent evidence of the charge dependence of lateral spreading for He/C/O in water at 50–430 MeV/u, … recorded as a gap in the ledger."). Integration is allowed with the gap stated in the ledger; V2-ION stays "in progress".
- **T1-ion and T12-ion failures are STOP (implementation defect) before any V-row is run.**
- Schedule checkpoint: if implementation exceeds 10 working days, record the cause and re-plan; no gate changes.
- Every failure, inconclusive result and consumed base is preserved.

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
| V1/V2/V1-held/V3 ionmc samples (≈ 9 × 1e6 warp-cpu) | ≈ 22 min |
| **TOPAS (REF)** | Every REF case is split into parts of ≤ 1200 s planned (P4 throughput), each a separate host call with its own seed (seed table extended by part index); the EM-only runs are sized the same way once P4 measures throughput (the V3-010 carbon full-physics throughput ≈ 38 ms/history would exceed the 1500 s per-call limit at 1e5 primaries). 3 smokes; 9 matched evidence runs; 12 held-out; 6 full-physics (V3, only if P6 passes); 7 cut variants: ≈ 37 runs × 1e5. Ion throughput unknown (P4). Cap 16 h total, one run at a time; materialize, parse, delete per run |
| **Total** | LV ≈ 1.5 h; HR ≈ 0.25 h; REF ≤ 16 h |

All steps run with `--single-process --workers 1`, thread limits 1, logged per the 2026-10-12
operator rules, and each host call ≤ 1500 s.

## F. Gates and integration criterion

1. Pre-freeze checks closed in the r5 order of section P (steps 1–7). Decision 0043 and this plan independently reviewed (scientific
   review + Codex) and committed before any result.
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
| 2 | G4EMLOW 8.8 (members extracted by path after acquisition) | `https://geant4-data.web.cern.ch/datasets/G4EMLOW.8.8.tar.gz` | 350 393 595 B | evaluation (I3) **only after** the P5 span/erratum checks; otherwise exploratory |
| 3 | libdEdx source with data files | TO VERIFY (GitHub/SourceForge) | < 5 MB expected | evaluation alternative / cross-check (P5b) |
| 4 | Sihver, Schardt, Kanai 1998, Jpn. J. Med. Phys. 18:1 (PDF) | J-STAGE, TO VERIFY | ≈ 1–3 MB | evaluation, measured (V3), if P6 passes |
| 5 | (optional) IAEA stopping DB tuples, C in H2O | via `https://www-nds.iaea.org/stopping/` | kB | exploratory: documents the measured-data limitation |

## Evidence archive

| Row(s) | Archive | Head SHA | Seed base | Result |
|---|---|---|---|---|
| (filled after the runs) | | | | |

## Amendments

(none)
