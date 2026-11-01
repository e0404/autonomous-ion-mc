# V3-008 acceptance criteria (DRAFT r3 — to be frozen before any result)

> **Status: DRAFT r3 (2026-11-01), not yet frozen.** Acceptance criteria are frozen BEFORE any result: this draft becomes revision 1 (frozen) by a reviewed commit after the pre-freeze checks P1–P9 are executed (they produce no qualification result; P4 smoke runs carry no scorers). Until then no validation row runs. Independent review: round V3-008-1 (`validation/plans/v3-008/review.md`). Seed bases per the shared ledger (V3-005B commit 06974d07): V3-008 reservation 20572004–20632003.


**Status.** DRAFT r3 (planner, 2026-11-01). r2 was revised against `research-models.md` and
`research-data-sources.md`; r3 applies the nine wording changes of review round V3-008-1 (REVISE).
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

Execution order: P8, P7, P2, P2b, P1, P3, P5/P5b, P6, P4 (last: the only step that runs an engine), P9.

| # | Check | Outcome feeds |
|---|---|---|
| P8 | Repository facts: the A16 baseline constant and the CI digest mechanism on develop (expected `f3a1dd62`); ASTAR last row = 1000 MeV total (loader test); species-id merge risk with PR #82 | R2, I2, I0 |
| P7 | Seed bases (section C): register the V3-008 bases in the shared seed-base ledger and in the guard test before any run; confirm the dated note on the V3-005 branch (V3-005 continues at 20641004) | C |
| P2 | Acquire and read the Geant4 v11.4.2 sources (`G4UrbanMscModel.cc`, `G4EmStandardPhysics_option4.cc`, `G4LindhardSorensenIonModel.cc`, `G4BetheBlochModel.cc`, `G4IonICRU73Data.cc`, `G4LindhardSorensenData.cc`). Establish (i) the msc model assigned to GenericIon and alpha in opt4, and (ii) the θ0 formula and its z-dependence. From these, compute **δ_pred(ion, z/R80)**, the expected relative lateral-σ offset TOPAS vs ionmc at z/R80 ∈ {0.5, 0.9}, by numpy Fermi–Eyges integration of the Urban θ0 scattering power against the selected ionmc model. If the Urban formula cannot be determined from source, δ_pred is "undetermined" (see V2). **P2 freezes the exact Urban formula, version and parameters** (θ0 expression, y definition, Z-dependent c1/c2) together with the sha256 of each source file read. **P2 also confirms the energy limit at which `G4LindhardSorensenIonModel` consults `G4IonICRU73Data`** (research note: scaled energy ≤ 2.5 MeV, i.e. about 2.5 MeV/u). Pre-declared: if the ion tables are used above 2.5 MeV/u, the class of I3 becomes "shared with TOPAS" and the C/O lineage of V1 changes accordingly (amendment before freeze) | V2, I3, V1 |
| P2b | P2b: the same δ_pred computation is applied to protons (z = 1) and compared, before any ion TOPAS run, with the archived V3-003 proton TOPAS lateral result (T15, already seen; a method check, not a model choice). The absolute discrepancy e_δ is recorded. Case B grades |ρ − 1| ≤ 5 % only if e_δ ≤ 2 %; if e_δ > 2 %, V2 falls to Case C (informative) for every ion whose |δ_pred| > 2 %. | V2 |
| P1 | Committed numpy Molière–Hanson ion reference, as defined in decision 0043 §3. It is validated first against Gottschalk's proton θ_Hanson tables (Be/Al/Cu/Pb, cached): agreement ≤ 1.0 % at every tabulated point, because both are the same theory and 1 % covers integration/convention differences. Then the analytic derivation of c_Ω and F_scr for candidate (d) (formula plus derivation committed), and θ for candidates (d), (c), (b), (a) at the I11 grid. **The pre-declared selection (d) → (c) → (b) → (a), first within 3 %, picks the default.** If none passes: STOP and amend | I11, MCS default |
| P3 | Committed numpy LS script: ΔL_LS − (Bloch + first-order Mott) for z = 6, 8, water, 2–430 MeV/u, tested against the Bloch limit at small z. Expected ≤ 0.1 % of S (estimate). If it exceeds 0.2 % of the CSDA range at any domain energy, an LS option is added (amendment) | decision §1, I3b |
| P5 | Data inspection of G4EMLOW 8.8 (`https://geant4-data.web.cern.ch/datasets/G4EMLOW.8.8.tar.gz`, 350 393 595 B, MD5 328330009df633f7e9b3a9f445745298), members `ion_stopping_data/icru73/z{6,8}_G4_WATER.dat` and `icru90/z6_G4_WATER.dat`. Checks: (i) energy span covers ≥ 2–430 MeV/u; (ii) units and per-nucleon axis; (iii) erratum provenance. The original ICRU 73 water used I = 67.2 eV and the erratum uses 78 eV; decide by the high-energy logarithmic slope against the Bethe L0 at two I values (aggregate verdict only, no values in Git) and by documentation; (iv) origin of `icru90/z6`. **P5b:** locate libdEdx (GPL) and, if it ships ICRU 73 (rev.) water files, use it as a cross-check or substitute | I3 |
| P6 | Measured-range candidate: Sihver, Schardt, Kanai, Jpn. J. Med. Phys. 18 (1998) 1 (J-STAGE). Check open access and licence, whether Bragg-peak/R80 positions for C and O in water are tabulated with uncertainties, and whether beam-line water-equivalent offsets are given. Bounded investigation: two access approaches, ≤ 2 h | V3 |
| P4 | TOPAS settings, checked on smoke runs (execution only, not graded): `alpha` vs `GenericIon(2,4)` particle naming; `GenericIon(8,16)`; BeamEnergy = total kinetic energy; acceptance of a very large electron production cut so that deltas deposit locally (e.g. `d:Ph/Default/CutForElectron = 1 km`; exact parameter TO VERIFY); run-time printout confirming the LS ion model and the ICRU90 flag (`/process/em/printParameters` or verbose); throughput per 1e5 primaries P4 smoke cases carry no dose, IDD or phase-space scorers (or their outputs are deleted unread by the runner); only the run log, EM-model printout and wall time are recorded. The measured throughput sizes every REF case into parts of ≤ 1200 s planned (section E) | V1, V2, V3, host budget |
| P9 | Freeze record: the committed hashes of the P1, P2 (incl. P2b) and P3 scripts and their outputs are recorded in decision 0043, so the freeze commit is provably after them (D-E8 pattern: the freeze commit has those commits as ancestors, the recorded hashes match, the tree is clean) | freeze |

## A. Verification rows

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| I0 | **Capability contract.** He/C/O sources run on python, warp-cpu and warp-cuda. Fail closed for: He-3 and any other ion; energy outside the domain or E0 + 6σ_E above it; `nuclear=True` with z ≥ 2; `nist-star` for C/O; ASTAR with a non-water material; LET medium ≠ water; unproducible species/generation tallies. Requested/effective configuration persists projectile, z, mass, `e_cut_ion_mev_per_u` and the MCS model id | every case raises `UnsupportedCombinationError` (or `valid=False` plus a named counter); fields present | capability | CI |
| I1 | Shared functions with a charge argument (z ∈ {1, 2, 6, 8}), the velocity-scaled arguments and the (d) correction row. U1 harness: twin vs Warp CPU f32/f64; the completeness test covers every new callable | decision 0001 classes as frozen in U1 ¹⁸ | backend parity | CI |
| I2 | **He table reproduction (SEEN, regression).** Bethe 4He at I = 75 eV vs ASTAR (75 eV), and Bethe at 78 eV vs ICRU 90 alpha, water (`compare_nist.py` extended). The I-value difference between ASTAR (75 eV) and the default model (78 eV) is stated in the result | S: max abs rel dev ≤ 0.30 % for e ≥ 10 MeV/u, ≤ 1.0 % for 2–10 MeV/u. CSDA range at 50/100/150/200/250 MeV/u: ≤ 0.20 % (ASTAR) and ≤ 0.25 % (integrated ICRU 90). Aggregates only in Git. Bounds = 1.5× the SEEN decision-0038 maxima, rounded up; derived from seen values, so I2 is a regression guard, not evidence | ion-specific tabulated (regression) | LV |
| I3 | **C/O table comparison.** Bethe (78 eV) vs ICRU 73 (rev.) water for C and O, and vs ICRU 90 for C if P5(iv) establishes its origin. Source per P5/P5b | S, e ≥ 10 MeV/u: max abs rel dev ≤ 1.5 % and RMS ≤ 0.8 %. S, 2–10 MeV/u: ≤ 6 % (effective-charge model spread plus PASS low-energy behaviour). CSDA range at the ends and standard points ≤ 0.5 %. Defended from theory: omitted LS ≤ 0.1 % (P3); Barkas/shell parameterisation ≤ 0.5 % (alpha experience 0.2 % SEEN, scaled for z³ terms); reference uncertainty ≈ 1–2 % (TO VERIFY). **If P5 and P5b both fail:** row = EVIDENCE GAP, intervention per protocol; it cannot pass. O has ICRU 73 only | ion-specific tabulated (largely independent of TOPAS above 2.5 MeV/u) | LV |
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
| I11 | **MCS model vs Molière–Hanson ion theory** (P1 reference, committed script): θ_model(x) vs θ_MH(x), x/R ∈ {0.1, 0.25, 0.5, 0.9}, standard points and ends | max \|θ/θ_MH − 1\| ≤ 3 % (justification in decision §3: f_dM ±2 % ⊕ ion theory 1–2 %) | source-model reproduction (calibration role; NOT independent) | LV |
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
- FACT (research note), confirmed by P2: **C/O** stopping in TOPAS = ICRU 90 proton array at equal velocity × Geant4 q_eff² + LS ΔL + Geant4 Barkas above about 2 MeV/u. Shared with ionmc: Bethe theory, I = 78 eV, the Barkas F table.
- FACT (research note), confirmed by P2: **He** stopping in TOPAS = the ICRU 90 alpha array, the same as I2's reference, so the He range in V1 is NOT independent of I2.
- Straggling and MSC models are independent implementations.
- Pre-declared: if P2 finds the ion tables used above 2.5 MeV/u, I3 becomes "shared with TOPAS" and this C/O lineage statement is amended before freeze.

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| V1 | **TOPAS EM-only IDD (matched variant)** at the standard points (blind). ionmc: warp-cpu f64, 1e6 histories, all EM on | (a) \|R80_ionmc − R80_TOPAS\| ≤ 0.35 % of R80 and ≤ 1.0 mm (TOST; σ from seeds and batches). Defence: Bethe-78 vs integrated ICRU 90 proton range differs by up to 0.18 % (decision 0038, SEEN); effective-charge, LS and Barkas construction differences ≤ 0.1–0.2 % (estimate). (b) Distal width (R80 − R20) within 10 % (TOST); straggling theory spread is a few %, and Geant4 fluctuations are an independent implementation of Bohr-type theory. (c) Peak position within 0.5 mm. (d) Peak/entrance ratio within 3 % (TOST). For helium, V1 tests transport mechanics and straggling only; its range part is bounded by I2 and adds no independent stopping evidence. | independent MC (C/O: partly shared stopping lineage; He: shared stopping data) | REF + LV |
| V1-held | As V1 at **He 220, C 430, O 430, C 100 MeV/u**, run ONCE, only after V1 has been graded and after any model amendment | same criteria | independent MC (held-out) | REF + LV |
| V2 | **TOPAS EM-only lateral (matched variant):** Sheppard-corrected σ of the lateral dose profile in a ±10 mm window (window rule as T15 ¹⁰) at z/R80 ∈ {0.5, 0.9}, standard points. δ_pred from P2, committed before any TOPAS run. **Case A, \|δ_pred\| ≤ 2 % for an ion:** grade the raw ratio. **Case B, 2 % < \|δ_pred\| ≤ 10 %:** grade the corrected ratio ρ = σ_TOPAS / [σ_ionmc (1 + δ_pred)], report the raw ratio. δ_pred is derived from the Geant4 source formulas and the Molière reference only, never from TOPAS output. **Case C, δ_pred undetermined or \|δ_pred\| > 10 %:** V2 is informative only for that ion, and lateral spreading has no independent MC evidence (declared gap, intervention candidate) | Case A: \|σ_ionmc/σ_TOPAS − 1\| ≤ 5 % and \|z\| < 3 at both depths. Case B: \|ρ − 1\| ≤ 5 % and \|z\| < 3 at both depths. Either way V2 asserts consistency of the transported MCS with an independent MC after a source-derived model difference, **not** an independent test of the z-dependence. Case B requires e_δ ≤ 2 % (P2b); otherwise Case C. V2 has no held-out analogue; a model change after V2 leaves lateral MC evidence as a declared gap | independent MC | REF + LV |
| V3 | **Measured range (C, O)** from Sihver, Schardt, Kanai 1998 (HIMAC), only if P6 passes. ionmc EM-only R80/peak position vs measured, at the published energies inside the domain. Nuclear allowance per energy: \|ΔR80(TOPAS full physics − TOPAS EM-only)\|, from additional TOPAS full-physics runs (opt4 + QGSP_BIC_HP + QMD, same geometry, 2 seeds × 1e5) | \|ΔR80\| ≤ (stated measurement uncertainty incl. beam-line WET ⊕ 0.5 % R) + nuclear allowance, at every energy used. Lineage: I_water = 78 eV calibration lineage, declared (range-consistency evidence, not an I-value test). **If P6 fails: row = EVIDENCE GAP → intervention** (inaccessible information: GSI Schardt 2007 precision Bragg curves, HIT Tessonnier 2017 / Krämer 2016 / Kurz 2012 numeric tables). He: no candidate; declared gap | measured | LV |
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
  2. P1 selects no candidate.
  3. V1 or V2 fails. Attribution with I3b/I3c and the P2 source reading is allowed. Any model change
     after seeing V1/V2 consumes them (preserved as failed). V1-held then becomes the sole gating
     independent-MC range evidence, at a fresh base with unchanged tolerances. δ_pred is never
     recomputed after a TOPAS output has been seen, and no parameter is fitted to TOPAS.
  4. Evidence gap (I3 data, V3 measured, V2 case C): an intervention per protocol. Integration is
     allowed with the gap stated in the ledger; V2-ION stays "in progress".
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
| P1, P3, I2, I3, I3b, I3c, I11 (numpy) | < 10 min |
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

1. Pre-freeze checks P1–P9 closed, in the order of section P. Decision 0043 and this plan independently reviewed (scientific
   review + Codex) and committed before any result.
2. Implementation, about 5–8 days:
   - charge/mass plumbing (tables, reference loop, kernels; no literal charge 1);
   - the T_ion (d) correction row;
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
