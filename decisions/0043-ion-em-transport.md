# 0043 — Helium, carbon and oxygen electromagnetic transport (revision 1 — frozen with the plan revision 1)

> **Status: frozen with the plan revision 1** (`validation/plans/v3-008-acceptance.md`; prepared 2026-11-01 by the planner at f694b747, round FREEZE-1 in `validation/plans/v3-008/review.md`). The decision is frozen by the same reviewed freeze commit as the plan, named by `validation/plans/v3-008/p9-freeze-record.md`. That commit is made only after the freeze conditions FC-1 (P1 completion) and FC-2 (P2 follow-up residuals and the P4 material printout) of the plan's §P freeze audit are disposed; until then this revision is a freeze candidate. No V3-008 transport or qualification result exists at f694b747. The DRAFT r6 status text below is kept as history.

> **Status: DRAFT r6 (2026-11-01, amended for Codex REVIEW-b469dff8, round V3-008-3), not yet frozen.** Commit history (exact): r3 was committed with the plan at 2632a46 (the task's first scientific commit, a draft); r5 at 4927958, amended at 3d2547a, 070b07f, 00af9a6 and b3eef90 while pre-freeze checks were open; r4 was not committed separately. The frozen revision is a separate reviewed commit after §P closes. r5 status text follows: r3 was committed at 2632a465 and reviewed in round V3-008-1; r4/r5 fold in the pre-freeze records P8, P7, P2, P1/P1b (committed at 6e155b56) and the review rounds V3-008-1 and V3-008-2 (`validation/plans/v3-008/review.md`; the six V3-008-2 wording changes applied). The decision is frozen by a reviewed commit only after the pre-freeze checks of the acceptance plan (§P, steps 1–7) close. The decision and plan are frozen before any TRANSPORT or QUALIFICATION result; the pre-freeze construction/calibration computations P1/P1b (Molière–Hanson reference, candidate tables) and P11 (I11-rel/I11-abs values, evaluated at 070b07f0 and regenerated at 00af9a6c) exist, are disclosed and are not qualification evidence; the P2/P10 source readings and the P5/P5b/P6 data inspections also exist and are disclosed (amended for Codex REVIEW-12614ced).


- Status: DRAFT r6 (planner, 2026-11-01), not frozen. History:
  - r2 was revised against the two research notes;
  - r3 applied the wording changes of review round V3-008-1;
  - r4 folds in the committed pre-freeze records at 6e155b56: `validation/plans/v3-008/`
    `p1-moliere-hanson-reference.md`, `p2-geant4-ion-models.md` and `p8-repository-facts.md`.
  - Main r4 changes: §3 (MCS model and the status of I11/V2) and the Geant4 lineage wording.
  - r5 applies the six wording changes of review round V3-008-2: the §3 selection statement, the
    exact-Ω specification, gating I11-rel, the V2 branches, the lateral-gap intervention, and the
    §6 Molière–Hanson row.
  - Round V3-008-4 (Codex REVIEW-e0b586fd, 2026-11-01; dated notes, earlier text kept): the operative §3 "Freezing" text replaced by the actual chronology (finding 2); the withdrawn independent-coverage claim for the f_dM extrapolation propagated through §3, §6 and the limitations (finding 3); the V5-SS budget committed as `validation/scripts/ion_em/ss_budget.py` with tests (finding 4, P10 §7.8); finding 1 in `p3-lindhard-sorensen.md`.
  - Round V3-008-5 (Codex REVIEW-6a14f4ec, 2026-11-01; dated notes, earlier text kept): finding 1 — the LS option implemented with default "off" by the dated §1 amendment (record `validation/plans/v3-008/p3b-ls-default.md`), the uncertainty statement extended; finding 2 — V5-SS uses the plan's TOST rule, and the statistics are raised to 6e4 TOPAS primaries / 4e5 ionmc histories per case (§3 dated note).
  - **Revision 1 (freeze revision; planner, 2026-11-01, at f694b747; dated notes, earlier text kept):** status "frozen with the plan revision 1"; the single-scattering clause of the lateral-gap wording amended and the second gap (f_dM extrapolation) added, both by orchestrator decision (§3 "Lateral evidence gap", revision 1 amendment); the pre-registration list of §3 "Freezing" updated to f694b747; the P2 follow-up and P4 outcomes for the Geant4 ion models recorded (Context, dated note); the open P1 completion pointed to the freeze condition FC-1 (§3 P1/P1b findings). No model choice, tolerance or row definition changed.
  - r6 applies the dispositions of Codex REVIEW-b469dff8 (round V3-008-3): the draft commit history (above); the P6 lateral EVIDENCE GAP propagated; the V5-SS model lineage and class (P10 §7; registration blocked until the Geant4 ion single-scattering model is read); the f_dM domain extrapolation (§3, note of 2026-11-01; P11b). Not accepted. To be reviewed
  independently and committed together with `validation/plans/v3-008-acceptance.md` before any
  V3-008 transport or qualification result exists; the pre-freeze construction/calibration
  computations P1/P1b and P11 exist, are disclosed and are not qualification evidence. The pre-freeze checks of the plan (section P) must be closed and their
  outcomes written into this decision first.
- Task: V3-008 (main line, `validation/plans/v3-next/selection-2026-11-01.md`, rank 1).
- Base: `v3/develop` at `c790df36`. The unmerged V3-005B/C/D elastic work (PR #82) is neither an
  input nor a dependency.
- Affects: physical accuracy (ion range, straggling, lateral spread), capability contract,
  architecture (charge as a transport parameter), validation strategy, data provenance.
- Markers: **FACT** = checked in code or a primary source during planning/research;
  **estimate** = hand arithmetic on standard formulas, to be reproduced by a committed script;
  **TO VERIFY** = not checked; **OPEN (research needed)** = unresolved choice with a pre-declared
  rule.

## Problem

V2-ION requires helium and carbon transport with ion-specific stopping and range evidence and an
oxygen domain with ion-specific evidence for every claimed capability. The engine (decisions
0037–0041) transports z = 1 projectiles only. `config.validate` rejects any non-proton source
("ions arrive with V3-008"). `TransportTables.from_stopping_tables` raises for z ≠ 1. The shared EM
functions are called with charge 1 (`self.one` / `one_d` in `transport/reference.py`,
`transport/kernels.py`). This decision fixes, for 4He, 12C and 16O primaries without nuclear
interactions:
- stopping source, effective charge, uncertainty, grid;
- straggling;
- multiple Coulomb scattering (MCS);
- kinematics, cutoff and domains;
- capability boundary to V3-009;
- data roles and lineage.

## Context (FACT unless marked)

**Stopping model.** `physics/stopping.py` evaluates
`S/ρ = K q² (Z/A) β⁻² [L0 + q L1 + q² L2 + M]`:
- exact Tmax with the ion mass;
- Sternheimer δ;
- Bichsel shell correction, held at η = 0.13;
- Ashley–Ritchie–Brandt Barkas term, with the F table transcribed from Geant4 `G4EmCorrections`;
- Bloch term;
- first-order Mott term `½ π α β q` for z ≥ 2. This is the same term as Geant4's `π α β z` in
  Geant4's 2× bracket convention, so the two are consistent.

**Effective charge.** The code's effective charge `q = z(1 − exp(−125 β z^(−2/3)))` is the
**Barkas (1963) form with coefficient 125**. The code docstring calls it Pierce–Blann; Pierce &
Blann (1968) use 0.95·137.04 ≈ 130.2. This decision uses the correct name. The docstring is
corrected during implementation, with no numeric change. The two coefficients differ by 1.8 % in
q² for C at 2 MeV/u and by < 0.1 % above 10 MeV/u (arithmetic).

**Earlier stopping results (SEEN).** Decision 0038 compared the analytic alpha stopping above
10 MeV/u:
- vs ASTAR at I = 75 eV: max 0.203 %;
- vs ICRU 90 at I = 78 eV: max 0.200 %;
- alpha CSDA range: ≤ 0.10 %.

These numbers are SEEN, so helium table rows are regressions, not blind tests.

**Projectiles and species.**
- `physics/projectiles.py` defines ALPHA (CODATA), CARBON12 and OXYGEN16 (AME2020 minus Z m_e).
- `species.py` registers ids 0–4; 12C and 16O are not registered.
- `TransportTables` index by total kinetic energy `E = A e` (V3-005A), so the exact log-log range
  construction (V3-003D) carries over unchanged to A > 1.

**Scattering functions.** `bohr_variance` and `scattering_power_dm` already take a charge
argument. `f_dM` is a fit of T_Hanson/T_ICRU to Molière/Fano/Hanson **theory** for protons
(Be/Cu/Pb, E ≤ 300 MeV; arXiv:0908.1413). Gottschalk calls it "not explicitly material
dependent".

**Registered data.**
- `nist-astar-water-2005`: ASTAR water at I = 75 eV (ICRU 49), up to 1000 MeV total = 250 MeV/u
  (consistent with `icru90.py`; last row not re-read, TO VERIFY by a loader test).
- `geant4-icru90-stopping-11.4.2`: proton and alpha water arrays only.
- No C/O table is registered.

**What Geant4 11.4.2 `g4em-standard_opt4` (TOPAS) does for ions.** Source: P2 record, from v11.4.2
sources pinned by sha256. FACT unless marked.
- **GenericIon (C, O) ionisation** = `G4ionIonisation` + `G4LindhardSorensenIonModel` (FACT [OPT4]
  L316–324, default constructor). TOPAS's instantiation is confirmed at P4 from the EM printout.
- **The LS model switches at 2.0 MeV scaled energy** (`fElimit`, escaled = E·m_p/M), i.e.
  ≈ 1.985 MeV/u for 12C and ≈ 1.984 MeV/u for 16O. It does not switch at 2.5 MeV/u.
- **Above the switch:**
  - dE/dx is the `G4BetheBlochModel` table;
  - along each step it adds q_eff²·(ΔL_LS + 2·BarkasCorrection·(z−1)/z), with q_eff² from
    `G4EmCorrections::EffectiveChargeSquareRatio`;
  - there is **no separate Bloch or Mott term**: the LS ΔL stands in for the higher orders.
  - The claim that `G4BetheBlochModel` returns the **ICRU 90 proton array at equal velocity** for
    GenericIon in G4_WATER is only **partly confirmed**. It is FACT (research note, scratch read); it
    is not in the pinned P2 set. Pinning `G4BetheBlochModel.cc/.hh` is a follow-up.
- **ICRU 73/90 ion data** (`G4IonICRU73Data`) are consulted only for escaled ≤ 2.0 MeV and only for
  2 < Z ≤ 80. They are **never used for alpha.** The pre-declared r3 trigger, "ion tables used above
  2.5 MeV/u", is **not met**, so no lineage amendment is triggered.
- **Alpha** is configured in `G4EmBuilder::ConstructCharged`, not in [OPT4]. That alpha uses the ICRU
  90 alpha array up to 1 GeV total is FACT (research note) and OPEN in the pinned set; the follow-up
  reads `G4EmBuilder.cc` (acquired), `G4ionIonisation.cc` and `G4BetheBlochModel.cc`.
- **Ion and alpha multiple scattering** = `G4hMultipleScattering("ionmsc")` with no model set. Urban
  is INFERRED strong as the default; the follow-up reads `G4hMultipleScattering.cc` (acquired).
- **Protons:** opt4 protons use WentzelVI (INFERRED strong; same follow-up).
- **Dated note 2026-11-01 (revision 1; P2 follow-up F1–F5 at 070b07f, P4 smokes at 41a1f8a/a66867d/ce4d3d6):**
  - GenericIon ("ionmsc") and alpha ("msc") get `G4hMultipleScattering` with no model set, hence `G4UrbanMscModel`
    (FACT [HMSC] L76); no `G4CoulombScattering` for them in opt4. The P4 printouts show UrbanMsc for both and
    WentzelVIUni + CoulombScat for the proton. So V2 is Case C (plan V2 row, revision 1 note).
  - Alpha: ionIoni BraggIon below 7.9452 MeV, `G4BetheBlochModel` above (P4 printout). `GenericIon(2,4)` is mapped to
    `alpha` (P4 verbose smokes), so He never uses LindhardSorensen. The He REF/V1 cases use `alpha`, and the He
    lineage is the BetheBloch lineage (orchestrator decision; plan §B).
  - Still not settled: whether the `G4BetheBlochModel` table is the ICRU 90 proton/alpha array (the file was hashed
    from a scratch download in `p3b-ls-default.md` §2, not acquired, and its ICRU 90 branch was "not settled", §5).
    The claim stays "research note, partly confirmed"; it is part of the freeze condition FC-2. The lineage rows
    (§6) already assume the shared lineage.

**Reference engines.**
- TOPAS transports ions. The only existing ion TOPAS case is a full-physics QMD smoke run
  (`carbon-water-290mevu-smoke`), which is not EM-only and not usable as evidence.
- The EM-only proton case `proton-water-150mev-lateral-emonly` (modules = 1
  `"g4em-standard_opt4"`, cut 0.05 mm) is the template for new ion cases.
- FRED aborts on 12C fragmentation; an EM-only ion run is untested (OPEN, not pursued in V3-008).
- MCsquare is proton-only.

## 1. Stopping power

### Selected

- **Default for 4He, 12C, 16O: the existing analytic layer (`BetheStoppingSource`)**, with no
  formula change and no retuning; I_water = 78 eV (ICRU 90), configurable and recorded. Reasons:
  - it is z-explicit: the Barkas term scales as q, the Bloch and Mott terms with higher powers of
    the charge;
  - it shares one code path with the qualified proton engine;
  - it is independent of the tables used to evaluate it;
  - it satisfies the analytic-plus-tabulated architecture of decision 0038.
- **Energy grid:** 1–500 MeV/u, 200 points per decade. Transport tables are on `ln(A e)`.
- **Tabulated helium alternative:** `NistStarStoppingSource` with ASTAR, selectable. When selected
  it is construction data for that run (decision 0038 rule). Note that ASTAR is at I = 75 eV, while
  the default model uses 78 eV.
- **Lindhard–Sørensen (LS): omitted, with a bound.**
  - LS contains the Bloch term exactly. LS − (Bloch + first-order Mott) is the higher-order Mott
    series, O((αz)²).
  - Estimated bound: ≤ 0.05–0.1 % of S and ≤ 0.1 % of the CSDA range for z ≤ 8 at ≥ 10 MeV/u
    (estimate, `research-models.md` RQ-S1).
  - The first-order Mott term already present is ≈ 0.5 % (C, 290 MeV/u) to 0.7 % (O, 430 MeV/u)
    of S (estimate).
  - Finite nuclear size matters only at γ ≳ 10, so it is irrelevant here.
  - Pre-freeze check P3: a committed numpy script (complex lgamma for the Coulomb phases δ_k,
    tested against the Bloch limit) reproduces the bound. If the CSDA-range effect exceeds 0.2 % at
    any domain energy, an LS option is added to `BetheOptions` and its default is decided in an
    amendment before freeze. The I3 tolerance (1.5 % in S) does not depend on this term.
  - **P3 result (2026-11-01, SEEN pre-freeze; `validation/plans/v3-008/p3-lindhard-sorensen.md`):**
    the point-nucleus LS sum (Dirac–Coulomb phases, tail-corrected to 1e-11) minus the exact ionmc
    Bloch + Mott₁ forms gives max r/S = 0.060 % (12C), 0.097 % (16O), 0.010 % (4He) and a CSDA-range
    shift of at most 0.060 % (12C) and 0.095 % (16O) over 2–430 MeV/u, positive everywhere. The
    estimate held (16O with a small margin); the 0.2 % trigger is not met; no LS option is added. The
    second-order coefficient comes out at 2.0–2.1 (Morgan–Eby closed form ≈ 1.70; informative). An
    exploratory, source-unsettled observation that Geant4 may apply the LS table at half the
    ionmc-convention magnitude is recorded there for the I3c/V1 attribution only.
  - **P3 re-opened (2026-11-01, Codex REVIEW-e0b586fd finding 1; SEEN):** with the mixed charge
    convention (LS at the bare nuclear charge z, Bloch and Mott₁ at the effective charge q — the
    Geant4 indexing) the omitted term is −1.22 % (C) / −2.56 % (O) of S at 2 MeV/u and shifts the
    CSDA range from 2 MeV/u by more than 0.2 % for end energies up to ≈ 7 (C) / 12 (O) MeV/u
    (≈ 0.4 / 0.9 µm), by 0.046 / 0.080 % at 430 MeV/u; the variants coincide where q/z > 0.999
    (above 16 / 24 MeV/u). r_mixed ≈ Bloch(z) − Bloch(q): the convention decides whether the
    close-collision Bloch term sees the bare nucleus. The pre-declared trigger is met, so the LS
    option is added to `BetheOptions`; its default (off = all-at-q, ls-eff, ls-bare, or a
    documented mixture) is decided in a reviewed amendment after a bounded physics assessment.
  - **Amendment (2026-11-01) — LS option** (Codex REVIEW-6a14f4ec finding 1, round V3-008-5;
    record `validation/plans/v3-008/p3b-ls-default.md`; SEEN pre-freeze, no I3, V1 or other
    reference value seen or used). The option is implemented: `BetheOptions.ls ∈ {"off", "eff",
    "bare_mixed"}` (`src/ionmc/physics/lindhard_sorensen.py`, `src/ionmc/physics/stopping.py`;
    tests `tests/ionmc/test_ion_em_ls_option.py`). It acts for z ≥ 2 only (z = 1 is a no-op, so the
    qualified proton model is unchanged). `"eff"` replaces Bloch(q) + Mott₁(q) by ΔL_LS(q);
    `"bare_mixed"` replaces them by ΔL_LS(z), the Geant4/catima-lineage indexing. **Default: `"off"`**
    (the existing model, all corrections at q; within 0.1 % of S of `"eff"` everywhere).
    - Reasons. The trigger was met only through the charge argument of the Bloch term:
      Bloch(z) − Bloch(q) is −1.29 % (C) / −2.64 % (O) of S at 2 MeV/u, while the genuine
      LS-beyond-Bloch + Mott₁ part is +0.05 / +0.06 %. No theory fixes that charge for a partially
      stripped ion (PASS treats screening through a screened-Coulomb potential; empirical effective
      charges such as Barkas, ZBL and HBG were fitted with everything scaled by q², absorbing
      higher-order effects), and Geant4 is itself inconsistent: `G4LindhardSorensenIonModel` indexes
      LS by the nominal Z with the q² prefactor (L157–158, L296–300), while `G4EmCorrections` on the
      `G4BetheBlochModel` path (L323) evaluates Bloch and Mott at the effective charge (L158–160,
      L695, L718); catima uses the mixed form only above 10 MeV/u; no published rationale was found.
      The plausible band at 2 MeV/u (+0.05 to −3.2 % of S for C, +0.06 to −7.5 % for O, the lower end
      a z²-weighted bare close-collision term) lies inside the declared low-energy model uncertainty
      and acts on the last ≈ 1–1.5 mm of range with ≤ 1 µm absolute effect.
    - Alternatives not chosen. `"eff"` is equally defensible (more complete at high energy) but
      changes S and range by ≤ 0.1 % and would put a production LS generator on the default path
      before freeze for that change; it remains selectable. `"bare_mixed"` must not be the default:
      it is a hybrid (z in the correction, q² weight), not a derived result, and it would remove a
      difference from the TOPAS reference by construction, reducing V1 independence without a
      physics reason.
    - What could decide it later: a screened-ion calculation (PASS/DPASS, or a project-side
      partial-wave transport cross section for a screened C/O ion), or measured C/O stopping in
      low-Z targets at 2–10 MeV/u from the IAEA database. **Never** the libdEdx ICRU73_NEW I3 rows
      (that would turn the I3 evaluation into calibration).
    - Consequences. Standard points: LS on/off shortens the range by 0.046 % (C) / 0.080 % (O),
      13 % / 23 % of the V1 R80 tolerance of 0.35 % (the convention itself changes it by
      ≤ 0.0003 %); stated in the V1 attribution budget. I3 2–10 MeV/u: the convention moves S by up
      to 1.2 % (C) / 2.6 % (O) at 2 MeV/u (20 % / 43 % of the 6 % tolerance), a declared model
      uncertainty, not a used-up tolerance; ≥ 10 MeV/u ≤ 0.1 %. I3b(ii) reports the three variants
      (report-only); I3c may use `"bare_mixed"` for the TOPAS-attribution diagnostic (exploratory).
    - Geant4 data-file identity: `G4LindhardSorensenData.cc` at tag v11.4.2 has sha256 70ae7c6b…,
      identical to the v11.2.0 file used in P3 §7 (closes that item; GitHub-tag fidelity assumed).
- **Effective charge:** the Barkas form with coefficient 125, in the stopping power only.
  - Below 1 MeV/u no stopping power is evaluated, and transport deposits locally below
    E_cut = 2 MeV/u (section 4).
  - Low-velocity alternatives (Hubert–Bimbot–Gauvin 1990, Ziegler/Brandt–Kitagawa) would differ by
    a few % in q² at 2–5 MeV/u and by < 0.5 % above 10 MeV/u (TO VERIFY). That moves the end of
    range by ≤ 0.01–0.02 mm (estimate), below the 0.05 mm action threshold. No change; the bound is
    documented.

### Uncertainty statement (metadata and docs; not acceptance tolerances)

- I_water 78 ± 2 eV (ICRU 90) gives about ±0.3 % of range.
- Omitted LS term: ≤ 0.1 % of S (P3 computed: 16O 0.097 %, 12C 0.060 %; range ≤ 0.095 %).
- **Dated 2026-11-01 (LS-option amendment, §1):** the charge convention of the close-collision
  (Bloch/LS) correction for C/O is a declared model uncertainty at 2–10 MeV/u (and below
  16 / 24 MeV/u): +0.05 to −3.2 % (12C) and +0.06 to −7.5 % (16O) of S at 2 MeV/u, below 0.1 % of S
  above about 10 MeV/u, ≤ 1 µm in absolute range. Not resolvable from theory; the default `"off"`
  sits at the upper end of the band.
- Shell correction held at η = 0.13 and the Barkas parameterisation: alpha experience is ≤ 0.2 %
  above 10 MeV/u (SEEN). For C/O the provisional budget is 0.5 % in S above 10 MeV/u, to be replaced
  by the I3 result once it exists, and recorded as an outcome.
- Reference-table uncertainties: ICRU 90 (protons, alphas) and ICRU 73 (with the 2009 water
  erratum) quote uncertainties of about 1–2 % at therapeutic energies (TO VERIFY).

### Evaluation sources (never construction in default runs)

| Ion | Reference | Status | Lineage |
|---|---|---|---|
| 4He | NIST ASTAR water, I = 75 eV (compared with Bethe at 75 eV) | cached, pinned, SEEN | ICRU 49 / Bethe theory |
| 4He | ICRU 90 alpha water array (Geant4 11.4.2 source), compared at 78 eV | cached, pinned, SEEN | **the same array TOPAS uses for alpha** → helium TOPAS range (V1) is not independent of this row |
| 12C, 16O | ICRU 73 water tables (Sigmund-revised per the README): G4EMLOW 8.8 `ion_stopping_data/icru73/z6_G4_WATER.dat`, `z8_G4_WATER.dat` (P5: spans 0.025–1000 MeV/u, 53 points) | **lineage only — use restricted to Geant4 by the data set's README; not an evaluation source** | ICRU 73 (PASS/BEST theory). In TOPAS these files are consulted only at ≤ 2.0 MeV scaled energy (≈ 1.98 MeV/u) and only for 2 < Z ≤ 80 (P2 FACT) |
| 12C (16O?) | ICRU 90 carbon water table: G4EMLOW `icru90/z6_G4_WATER.dat` | **lineage only — use restricted to Geant4 by the data set's README; not an evaluation source**; origin OPEN (P5). **No ICRU 90 oxygen table exists** (ICRU 90 covers p, α, C) | ICRU 90 |
| 12C, 16O | **libdEdx ICRU73_NEW** (archive b941a0c3…, GPL-3; `data/raw/ICRU73_NEW.dat` 749768c1…, grid 53 energies 0.025–1000 MeV/u; target 276 liquid water) — the **I3 evaluation source** (P5b, `validation/plans/v3-008/p5b-libdedx-icru73.md`) | exploratory until I3 is frozen. TO-VERIFY: equality with the 2009 Sigmund water revision; the I-value of the water values (libdEdx metadata lists 75.0 eV for target 276; the model is compared at 78 eV and the difference is stated as a bounded contribution) | ICRU 73 (revised, as distributed by libdEdx); independent of TOPAS above ≈ 2 MeV/u |
| 12C, 16O | IAEA stopping database (v2026-09) | **Limitation (FACT, Montanari et al. arXiv:2402.03080):** water has measured data only for H, He, Li and C, and the carbon data lie near the stopping maximum (≲ 1 MeV/u). **No oxygen in water. No measured C/O water stopping at therapy energies exists publicly.** Not curable by intervention; declared | measured, outside the domain |

All reference values stay in the provenance cache. Git holds pins, hashes, aggregate ratios and
verdicts only (decision 0038 policy).

### Rejected alternatives

- **z²-scaled ASTAR/PSTAR for C/O.** It omits the Barkas term (∝ z³) and the Bloch/Mott terms
  (1–2 % at tens of MeV/u), and the scaled table would become construction data.
- **ICRU 73 (rev.) tables as the default C/O source.** Reasons, revised after the lineage research:
  - (i) The only candidate tabulated C/O reference would become construction data, leaving C/O
    stopping without a tabulated evaluation.
  - (ii) The files' span and erratum provenance are not yet established.
  - (iii) No ICRU 90 oxygen table exists.
  - (iv) Kept as an evaluation source, these tables are also the C/O check that is **most
    independent of TOPAS**, because Geant4 consults them only at ≤ 2.0 MeV scaled energy
    (≈ 1.98 MeV/u) and only for 2 < Z ≤ 80 (P2, FACT). Making them the default would turn that
    independent check into a shared-construction comparison.
  - The r1 reason ("TOPAS uses these data") was wrong above about 2 MeV/u and is withdrawn.
- **The Geant4 construction as default** (`G4BetheBlochModel` table, the ICRU 90 proton array at equal
  velocity per the research note, plus q_eff²·(ΔL_LS + 2·Barkas·(z−1)/z)). It would make ionmc's C/O stopping identical in construction to TOPAS's, so V1
  would test only transport mechanics. It is kept only as an **exploratory attribution diagnostic**
  ("TOPAS-like rebuild", row I3c), buildable from the already registered ICRU 90 proton array plus
  an LS ΔL computation. It has no evidence role.
- **Port of `G4LindhardSorensenIonModel`.** It is a fit table with shared lineage, and LS is below
  the action threshold (above).
- **SRIM and MSTAR.** Both are closed or executable-only software; their outputs are not versioned
  datasets, and MSTAR's high-energy part is ASTAR-anchored.

## 2. Energy-loss straggling

- **Selected:** the Bohr variance `σ² = (K/2)(Z/A) ρ x z² Tmax (1/β² − 1/2)` with the **bare charge
  z** and the exact ion Tmax, sampled by `bohr_gamma_v1` (decision 0039, amendment 3). No new
  sampler.
- **Validity:**
  - The ions are fully stripped (q/z > 0.99 above about 10–15 MeV/u).
  - The Vavilov κ = ξ/Tmax per 1 mm step is about 1 for C at 290 MeV/u (estimate), against about
    0.1 for 150 MeV protons. Ions therefore sit closer to the Gaussian regime: the Gamma shape is
    large and its skewness small, and the additive-shape step-independence argument carries over
    unchanged. Delta tails are irrelevant because deltas are deposited locally.
- **Charge-exchange straggling:** it acts only below a few MeV/u (TO VERIFY: Sigmund; ICRU 73),
  i.e. over ≲ 0.4 mm (C) or 0.3 mm (O) of range. It contributes ≪ 1 % of σ_R² (estimate) and is
  neglected. Using q² instead of z² would change only that segment and is rejected.
- **Non-Bohr effects** (shell/binding, correlation; Tschalär; Yang, O'Connor & Wang 1991): a few %
  in σ² (TO VERIFY), i.e. 1–2 % in σ_R. This is below the 10 % amendment trigger; documented.
- **Evidence intent:**
  - Exact scaling identity σ_R/R(ion) = (m_p/M)^½ σ_R/R(p) at equal velocity with corrections off
    (theory, I4b).
  - Transport vs quadrature (I8).
  - Report-only comparison with the literature statement σ_R/R(C) ≈ σ_R/R(p)/√12 ≈ 0.3–0.4 %
    (Schardt, Elsässer, Schulz-Ertner, RMP 82 (2010) 383; TO VERIFY), band ±15 % (I8b).
  - Distal width vs TOPAS (V1(b)).
  - Measured Bragg-peak widths mix beam energy spread and fragment tails, so they are informative
    only.

## 3. Multiple Coulomb scattering

### Theory used for the choice

These are the Molière–Bethe relations as written by Gottschalk (arXiv:1804.00022, 0908.1413) and
implemented in the committed P1 reference, `validation/scripts/ion_em/moliere_hanson.py`:
- χ_c² ∝ z² Σ w_i Z_i²/A_i · ∫dt/(pv)². This is Molière's Z² form, with the Fano correction to b for
  the atomic electrons.
- χ_a,i² = χ_0,i² (1.13 + 3.76 α_i′²), with χ_0 ∝ Z^(1/3)/p and α′ = z Z α/β.
- Ω = χ_c²/(1.167 χ_a²), with B − ln B = ln Ω + Δb_Fano.
- θ_0 = χ_c √(B − 1.2)/√2.

The ion/proton relation at equal velocity and equal x/R is exact, and is tested in P1:
- the ion slab is t_ion = (A*/z²) t_p;
- hence Ω_ion/Ω_p = A*/F_scr, where A* = M/m_p and F_scr is the Z_i²-weighted ratio of the Born
  screening factors;
- Δb_Fano depends on β only, so it cancels in this ratio.

### P1/P1b findings (recorded at 6e155b56; disclosed here before any model choice is frozen)

- **Reference validity.**
  - The Molière/Fano/Hanson reference reproduces Gottschalk's proton θ_Hanson for Be, Al and Cu to
    ≤ 0.21 %.
  - **For Pb, P1's pre-declared ≤ 1.0 % proton validation FAILED**: +1.10 % and +1.32 % at
    x/R1 ≥ 0.9. The failure is recorded as FAILED and preserved. Attribution on ICRU 49 ranges is
    still to do (P1 completion); clearing it later does not erase the failed record.
    **Revision 1 (2026-11-01):** at f694b747 P1 completion (plan §P order step 5) has no record. It is the freeze
    condition FC-1 of the plan's §P freeze audit (`p9-freeze-record.md` §0): executed before the freeze commit, or
    recorded "not done" by a dated orchestrator amendment. Pb stays FAILED unless cleared.
  - **Water is not validated directly**: Gottschalk's table has no compound block.
  - **The compound Fano rule** (electron term weighted by Z, nucleus by Z²; 1/(pv)² slab average of D)
    **is derived, not sourced.**
  - The proton f_dM sits **+2.5 to +3.5 %** (x/R 0.1–0.5) and +1.3 to +2.1 % (x/R 0.9) above this
    water reference. That is more than Gottschalk's stated ±2 % for f_dM. Either the derived compound
    rule or f_dM's water accuracy is off, and the available sources cannot decide which.
  - **Consequence:** the absolute water reference has its own uncertainty of order 3 %. It cannot
    arbitrate between models at a 3 % tolerance.
- **The r3 constant was wrong.** The r3 estimate c_Ω ≈ 0.13 was non-relativistic; P1 found
  n_Ω ≈ 0.90–1.12 and c_Ω ≈ 0.09–0.11. Codex REVIEW-98315d4e on 6e155b56 (being fixed in code) adds
  that **n_Ω must be the full derivative of the proton Ω (or b), including the β-dependence of the
  screening factor and of the Fano term**. This affects only the record of the superseded linearised
  (d); the exact root below needs no n_Ω.
- **The P1 driver evaluates the exact I11 grid**: ion-specific energies (standard points and ends) ×
  x/R ∈ {0.1, 0.25, 0.5, 0.9}.
- **Seen candidate predictions** (P1b, water, report only): max deviation from the Fano reference.

  | Candidate | Max deviation |
  |---|---|
  | (d), linearised c_Ω | about +6 % |
  | (c) slab proxy | +10 to +19 % |
  | (b) | +7 to +16 % |
  | (a), velocity scaling only | −3.4 % He (50 MeV/u), −2.9 % C, −2.5 % O |
  | (e) Gottschalk's ion f_I | −4 to −8 % |

- The r3 rule (first candidate within 3 %) selects nothing. That is STOP 2.

### Selected (r5): the exact-Ω mapping of the proton model

**Selection statement (review round V3-008-2, verbatim).** The model form is the exact version of
r3's first-ranked candidate (d); the exact root was named in P1 OPEN 5 after the first (no-Fano)
candidate table had shown the linearisation error, and before the P1b tables. The selection rule was
revised after the r3 rule reached STOP 2 on seen P1/P1b numbers: the Molière–Hanson water reference
failed its own pre-declared proton validation for Pb (recorded as FAILED, preserved), is not
validated for water, and rests on a derived compound rule, so it is demoted from selector and gate to
an informative report. The seen tables were not used to choose among forms; the predictability of the
exact-Ω result against the reference from the seen decomposition is disclosed.

**Model.**

```
T_ion(pv, p1v1; material) = f_dM(pv_eq, p1v1_eq) · (z E_s / pv)² / X_S
p1v1_eq / pv_eq = p1v1 / pv          (the memory term 1 − (pv/p1v1)² is unchanged)
pv_eq is the root of   b_p(pv_eq, p1v1_eq; material) = b_ion(pv, p1v1; z, M, material)
b = ln Ω + Δb_Fano   (Molière collision number of the slab traversed from p1v1 down to pv along the
                      CSDA path, with the Fano term evaluated at each projectile's own velocities)
```

**Matching quantity.**
- The root matches b, i.e. ln Ω plus the Fano correction, each evaluated at the respective
  projectile's velocities along its own slab.
- It does not match ln Ω alone: at the proton-equivalent state the velocities differ from the ion's,
  so Δb_Fano does not cancel.
- b_ion and b_p are computed with the same formulas as the committed P1 reference (Molière χ_c, χ_a
  with the Born factor, the Fano Δb, the 1/(pv)² slab average), from the project range tables.

**Two-dimensional dependence and step independence.**
- The mapping depends on (pv, p1v1): the ion's current pv and its initial p1v1.
- It is evaluated per step from the particle state. The current pv is the step's midpoint value, as
  in decision 0039. p1v1 is carried with the particle.
- So the per-step scattering power is a function of the particle state only. Step independence, the
  random hinge, the clamp and the birth-step average carry over.

**Evaluation.** Default: a per-step root solve.
- The 1-D monotone root in pv_eq uses tabulated b-slab integrals per material and projectile:
  cumulative ∫dt/(pv)² and ∫ ln χ_a² dt/(pv)² on the range grid, so the slab integral from p1v1 to pv
  is a difference of two table reads.
- **Allowed alternative:** a tabulation of pv_eq/pv on a (pv, p1v1) grid per material and ion, read
  by bilinear interpolation in (ln pv, ln p1v1). It is allowed only if a committed test bounds the
  interpolation error of the resulting θ (or T) by ≤ 0.1 % over the domain.
- Both the choice and its bound go into the formula, the script and the tests before P9.

**Physical content.**
- The Rutherford prefactor carries the ion's own z/(pv).
- The proton model f_dM carries the remaining dependence. Its absolute water calibration is
  evidenced for protons (decision 0039: U4/U4b source reproduction; T15 TOPAS lateral, status per
  the develop records, TO VERIFY).
- Theory supplies only the **relative** ion/proton change, through b. That relative quantity is the
  one Molière theory fixes most reliably: the compound-rule uncertainty and f_dM's absolute offset
  largely drop out of a ratio at equal velocity and x/R.

**Rejected.**
- (a): rejected **on structure**. It omits the z- and A-dependence of the collision number.
  - Interpretation (not the basis of the rejection): its seen near-pass is plausibly an artefact of
    a reference with about 3 % uncertainty combined with f_dM's opposite offset.
- (b), (c), (e): also rejected on structure. They lack the Born-screening term or are mis-scaled.
- Option (ii), an independent arbiter as selector: none is confirmed (P6, P10), and the TOPAS
  reference model is unsettled (P2 follow-ups / P4).

**Declared model uncertainty** (theory-derived, not validated): ±5 % on the ion θ at the I11 points,
recorded in result metadata and docs. Components, combined in quadrature (≈ 4–5 %):

| Component | Size |
|---|---|
| f_dM's water accuracy | ≤ 3.5 % |
| Born-screening form at α′ ≈ 0.5–0.7 | 1–2 % |
| z-independent Fano D (no Bloch-type electron terms) | ≤ 1 % |
| Nuclear size | ≤ 1–2 % |
| (r6) Extrapolation of f_dM beyond its proton-qualified pv (P11b: drift against Molière–Hanson ≤ 1.5 pp, protons 300 MeV → 3 GeV) | ≤ 1.5 % |

With the r6 component, the quadrature of the upper values is ≈ 4.9 %, so the ±5 % statement stands.

**Chronology (exact; Codex REVIEW-42ab68dd).** The mapping, the tests and the first I11-rel/I11-abs evaluation were produced together at 070b07f0, in one commit while the plan is unfrozen. Decision 0043 §3 and plan row P11 required the mapping and tests to be committed BEFORE the evaluation; that order was not kept. The values are SEEN construction/calibration results, not qualification evidence; a later rerun cannot make these points independent. I11-rel remains the gate as pre-registered (2 %, STOP on failure, tolerance never changed): evaluated pre-freeze; not blind; the freeze carries it as a disclosed gate. The first evaluation also took the proton at the ion's numeric MeV/u instead of at equal velocity; corrected (E_p = E_u A m_p / M) and regenerated, below.

**Dated note 2026-11-01 (r6; Codex REVIEW-b469dff8 finding 4): f_dM domain extrapolation (P11b).**
- f_dM is an explicit function of (pv, p1v1) and is qualified for protons ≤ 300 MeV, i.e. p1v1 ≤ 527 MeV.
- At fixed memory ratio, b_p is monotone in pv. Matched b is therefore in the qualified domain exactly when p1v1_eq ≤ 527 MeV. b-matching does not make the evaluation in-domain.
- SEEN pre-freeze, aggregates in `p11-exact-omega-mapping.md` §P11b:
  - Only He 50, C 80 and most of O 100 are in domain (11/36 I11 end states).
  - At every standard point and upper end, all of θ² comes from proton-equivalent states above 300 MeV: up to 0.8 GeV (He) and 2.5 GeV (C/O), i.e. up to 0.78 decades in lg pv.
- Extending the Bethe proton table to 8 GeV only made the root computable. It does not validate f_dM there.
- **Choice:** the model is kept unchanged. It is not modified after the seen numbers.
  - No in-domain formulation exists for this domain: a b-only form extrapolates in b (b_ion ≤ 17.06 against ≤ 14.95 for qualified protons).
  - A fail-closed clamp would remove most of the declared domain (§5). That is a scope decision outside this revision.
- **The extrapolation is:**
  1. declared;
  2. bounded against Molière–Hanson proton theory: θ_fdM/θ_MH drifts by +1.2 to +1.5 pp from 300 MeV to 3 GeV, smooth and logarithmic, as it does in-domain (row I11-ext, report-only, SEEN); this is theory-anchored, not independent;
  3. carried as a separate ≤ 1.5 % component of the declared ±5 %;
  4. persisted in result metadata (I0: the p1v1_eq energy of the source and `fdm_extrapolated`), with fail-closed proton b-tables to ≥ 3 GeV;
  5. covered, for the independent part, only through V5-SS. Its standard cases already map to extrapolated states (0.46–1.9 GeV); the r6 end cases He 250, C 430 and O 430 reach 0.8–2.5 GeV. This coverage is subject to the V5-SS STOP (P10 §7.4) and to its class (shared screened-Coulomb lineage, independent of f_dM).
     **[Superseded 2026-11-01 (round V3-008-4, Codex REVIEW-e0b586fd finding 3): item 5 is withdrawn.** Under V5-SS branch B the common-velocity ion/ion ratio covers the extrapolation only differentially (≤ 1.5 pp, below its 3 % tolerance), and V5-SS shares the screened-Coulomb lineage. There is no independent coverage. The gap is named in the lateral-gap intervention filed at the freeze (orchestrator decision 2026-11-01 (a), below), and the extrapolation stays a declared limitation (decision (c)).**]**
- I11-rel does not test the extrapolation: its C/O equal-velocity proton is itself above 300 MeV.

**Freezing (operative text corrected 2026-11-01, round V3-008-4, Codex REVIEW-e0b586fd finding 2; the superseded r5 bullet is quoted at the end of this block and is no longer operative).**
- **Actual chronology (disclosed; see "Chronology (exact)" above).** The exact-Ω mapping (`validation/scripts/ion_em/exact_omega.py`), its b formulas, its evaluation choice (the per-step root solve), the I11-map and I11-rel tests and the first I11-rel/I11-abs evaluation were produced together, in one commit at 070b07f0, while the plan was unfrozen. The equal-velocity correction was regenerated at 00af9a6c. P11b (I11-ext and the extrapolated θ² fractions) was computed at b3eef90b, after those values had been seen. The required order (mapping and tests committed before any evaluation) was not kept.
- **What is SEEN.** All P11 values (I11-rel, I11-abs, I11-ext, the P11b fractions) are SEEN construction/calibration results, not qualification evidence. The P9 hashes of the mapping, its script and its tests record what was evaluated; they do not establish pre-registration of these rows.
- **What the freeze revision may claim as pre-registered.** Only rows whose values have not been seen at the freeze commit (per the records at ce4d3d6a; no ionmc transport run, no scored reference-engine run, no LV/HR qualification step):
  - capability and implementation rows: I0, I1, I5, I7, I16, C-res, T1-ion, T12-ion, T13-ion (deterministic or statistical checks not yet run);
  - theory identities not yet evaluated on the implementation: I4, I4b;
  - stopping comparisons not yet made: I3 (libdEdx ICRU73_NEW against ionmc; P5/P5b inspected the data only), I3b (i), I3c;
  - transport rows: I6, I8, I8b, I9, I10, I12, I13, I14, I15;
  - reference rows: V1, V1-held, V2, V4, and the V5-SS statistic (its planning budget, expected bias and discrimination values are SEEN planning values, P10 §7.5.1 and §7.8; the V5-SS outcome is not).
- **Rows carried as disclosed SEEN, not pre-registered:** I11-rel (the gate, tolerance unchanged), I11-abs and I11-ext (report-only), I11-map (its deterministic identities ran at 070b07f0), I2 (regression on seen table values), I3b (ii) (the P3 curve), and the P1/P1b reference tables.
- Nothing is fitted to the reference, TOPAS or transport output. (Unchanged.)
- Superseded r5 text, quoted: "The exact-Ω mapping, its b formulas, its evaluation choice (root solve or a bounded tabulation) and the I11-map and I11-rel tests are committed **before the mapping is evaluated against the Molière–Hanson reference**. Their hashes go into P9."
- **Dated note 2026-11-01 (revision 1; lists updated from the records at ce4d3d6a to f694b747).** The src mapping
  implementation and its tests (683480b, rebound at f694b747) ran before the freeze. The following sub-items are
  therefore SEEN implementation checks, not pre-registered: in I0, the `fdm_domain_status` fields, the
  `fdm_extrapolated` flag at the standard and end cases, the fail-closed mapping-edge rule (4000 MeV) and the
  `ls_requested`/`ls_effective` persistence; in I1, the twin-vs-Warp CPU f64/f32 parity and the completeness check of
  the exact-Ω mapping callables; and the I11-map residuals of the src implementation. The remaining I0 and I1
  sub-items (species, the other fail-closed rules, the charge-argument functions) are unseen. The LS-option tests
  compare the three variants with each other, not with the I3 reference, so I3 and I3b (i) stay unseen. R2 is not yet
  evaluated against an ion-enabled build. The full per-row list is `p9-freeze-record.md` §6.

**Rows.** I11-map (CI, gating): implementation of the root. **I11-rel (LV, gating)**: the relative
check, failure STOP. **I11-abs** (informative): the absolute comparison, reported with the
reference's own uncertainty.

### Independent evidence for ion MCS and the TOPAS reference model (P2 findings; follow-ups close BEFORE freeze)

**Opt4 ion multiple scattering.**
- GenericIon gets `G4hMultipleScattering("ionmsc")` with no model set (FACT [OPT4] L161, L323).
- That this model is Urban is INFERRED, not yet confirmed.
- It is settled before freeze by:
  - reading `G4EmBuilder.cc` (40782f5c…) and `G4hMultipleScattering.cc` (a6643869…), both now
    materialized in the v3-008 cache;
  - the `G4VMscModel` defaults;
  - `G4IonisParamMat.cc` (990c6f31…) and `G4IonICRU73Data.cc` (cf771a5a…);
  - the P4 scorer-free smokes, with the EM-table / printParameters printouts archived: msc model,
    step-limit type, RangeFactor and lateral displacement for GenericIon, alpha and proton.

**Urban θ0** = 13.6 MeV·|charge|·√y·invbetacp·(c1 + c2 ln y), y = t/X0 (FACT [URB] L1023–1027,
L1115–1125).
- There is no charge, β or Born term in the logarithm.
- The θ0 core is truncated, and a power-law tail is added.
- The per-length width depends on the step length.
- TOPAS ion steps depend on the bundle design.
- The r3 Fermi–Eyges δ_pred is withdrawn, and P2b as defined is very probably inapplicable.
- **Scratch estimate** (not frozen): TOPAS C/O lateral widths are plausibly 10–20 % narrower than
  ionmc, if Urban.

**V2 branches (pre-declared, verbatim).**
- "If Urban: V2 informative (Case C), the raw ratio and expected direction reported."
- "If WentzelVI/single scattering: V2 is decided by an amendment before any ion TOPAS run."
- The r4 Urban-emulator extension is withdrawn from this plan. It can return only by a reviewed
  amendment before any ion TOPAS run.
- V2 has no held-out analogue: a model change after V2 leaves lateral MC evidence as a declared gap.

**Single-scattering reference (P10).** Before freeze, P10 checks whether a thin-slab angular-width
check against a Geant4 single-scattering configuration is available, and at what cost. Candidates:
`G4CoulombScattering`, or a g4em-standard_SS-type module if OpenTOPAS 4.3 offers it for GenericIon.
**P10 outcome (closed for coverage):** FACT, `G4EmStandardPhysicsSS.cc` v11.4.2 (cache 8aeaf91c…): L213–219 `particle = G4GenericIon::GenericIon(); G4ionIonisation* ionIoni = new G4ionIonisation(); … ionIoni->SetEmModel(new G4LindhardSorensenIonModel()); ph->RegisterProcess(ionIoni, particle); ph->RegisterProcess(new G4CoulombScattering(false), particle);` (no msc process for GenericIon); L222 `G4EmBuilder::ConstructChargedSS(hmsc);` → alpha `G4ionIonisation` + `G4CoulombScattering(false)` ([EMB] L158–160) and protons `G4hIonisation` + `G4CoulombScattering(false)` ([EMB] L212, L221–222); L107 `param->SetMscThetaLimit(0.0);`. The r5 thick-slab design was infeasible (P10 §3: ≈ 10²–10³ h). **V5-SS is the held-out lateral qualification row** on the P10 design: 1 mm water slab, 2e4 primaries per case, He 150, C 290, O 350 MeV/u and protons at matched velocity; statistic = the ion/proton ratio of the rms projected angular width compared with the model's ratio; pre-registered tolerance |ratio_model/ratio_ref − 1| ≤ 3 %; run ONCE after freeze, only if a scorer-free smoke confirms the OpenTOPAS 4.3 alias, the process list and a cost ≤ 4 h REF. V1-held and V3 provide no lateral evidence. **[Superseded 2026-11-01 (round V3-008-4): the cases and the ion/proton statistic of this paragraph are replaced by branch B (ion/ion ratios at a common velocity, 8 ion cases, no protons; the dated note below and the V5-SS row). V5-SS is shared-lineage single-scattering evidence, not independent physical evidence (orchestrator decision 2026-11-01 (a)–(b), below).]**

**r6 (2026-11-01; Codex REVIEW-b469dff8 finding 3; P10 §7): Coulomb-model lineage and evidence class.**
- **Model per particle.** With no model set by the constructor, `G4CoulombScattering` selects `G4eCoulombScatteringModel` + `G4WentzelOKandVIxSection` for the proton, and `G4IonCoulombScatteringModel` for alpha and GenericIon (FACT, [CS] L120–148).
- **Proton path** (FACT, P10 §7.3):
  - Wentzel screened Rutherford with Molière's Born screening, min(1.13 Z, 1.13 + 3.76 (zZα/β)²), times a Thomas–Fermi-type radius;
  - electron scattering up to min(cut, T_max);
  - an exponential-type nuclear form factor (default OPEN);
  - a spin-½ factor.
- **Ion path:** OPEN. The files are listed in P10 §7.4. **STOP on this sub-item:** the V5-SS registration, its branch (A: same physics on both sides, ratio design kept; B: proton invalid as the z = 1 reference, so ion/ion ratios at a common velocity or withdrawal) and the freezing of its 3 % tolerance wait for them. Under branch A the 3 % is re-justified from the source: the non-cancelling terms sum to ≲ 1 % (P10 §7.6).
- **Class:** held-out single-scattering MC.
  - Independent of the multiple-scattering approximation: f_dM, Highland, the Molière angular integration and the b-matching.
  - Not independent of the screened-Coulomb cross-section physics: the same Born screening form on the proton path.
  - It is evidence that the mapping realises the shared cross-section model, not independent physical evidence of the charge dependence of screening.
- **End cases** He 250, C 430 and O 430 (TOPAS cases 28–33) are added for the f_dM extrapolation coverage. **[Superseded 2026-11-01 (round V3-008-4, finding 3): under branch B He 250, C 430 and O 430 are cases 28–30 of the ion/ion design and give differential coverage only; cases 31–33 are re-assigned or retired (plan §C). No independent coverage of the extrapolation is claimed.]**

**Dated note 2026-11-01 (continuation of round V3-008-3 finding 3; P10 §7.4.1–7.7.1): the ion path is read; V5-SS branch B.**
- **Ion path (FACT).** `G4IonCoulombScatteringModel` + `G4IonCoulombCrossSection` (ef9b31d4…, 87989b93…; headers a1a9a4f5…, cd5d0b5c…):
  - nuclear scattering only; target electrons are not scattered;
  - the universal screening radius a0 (z^0.23 + Z^0.23), which depends on the projectile charge; no 1 + exp(−0.001 Z²);
  - the Born factor 1.13 + corr · 3.76 α′², with corr = (10π Z z α)^0.04, also for hydrogen and with no cap;
  - no nuclear form factor, no spin term;
  - exact CM kinematics, with β_rel = β_lab identically; the bare charge.
  - `G4EmParameters.cc` (aaa143c1…) defaults: NuclearFormfactorType exponential, ScreeningFactor 1, FactorForAngleLimit 1.
- **Branch B.** These differ from the proton path in screening, electrons and form factor, so the proton is not a valid
  z = 1 reference. With the r6 design, an exact Molière–Hanson ionmc would deviate by +3.9 to +5.1 % (planner
  scratch, SEEN pre-freeze, aggregates).
- **V5-SS is redesigned:**
  - ion/ion ratios at a common velocity, all on the ion path: He/C and He/O at 150 and 250 MeV/u, C/O at 430 MeV/u;
  - 8 cases × 2e4 primaries, 1 mm water;
  - an electron energy cut ≥ T_max, so no electron scattering of any kind;
  - the proton cases are dropped.
- **Tolerance** (set now from the source analysis, before any V5-SS run): 3 % per pair, kept, with σ ≤ 1 %.
  - Budget: reference-model bias ≤ 0.9 % (computed −0.17 to −0.69 %); ionmc against Molière–Hanson in the ion/ion
    ratio ≤ 2 % (I11-rel gate, SEEN spread ≤ 1.13 pp); statistics 0.7–1.0 %.
  - Combined: quadrature ≈ 2.4 %, conservative 3.0–3.1 %.
  - The He pairs exclude a b without z-dependence (+6.4 to +8.4 %), and the He/O pairs exclude a missing Born term.
  - **2026-11-01 (round V3-008-4, finding 4):** these figures are reproduced by the committed module `validation/scripts/ion_em/ss_budget.py` with `tests/ionmc/test_ion_em_ss_budget.py` (sha256 in P10 §7.8 and plan P9); every recorded aggregate was reproduced unchanged (quadrature 2.41 %, conservative 3.02–3.14 %).
- **Class:** held-out single-scattering MC on the ion path.
  - Independent of the MCS approximation.
  - Shares the screened Rutherford law and the Born form; differs in the radius and the electrons, both pre-computed
    as bias.
  - Not independent physical evidence of the charge dependence of screening.
- **Consequence for the f_dM-extrapolation note above (item 5).** The common-velocity ratio covers the extrapolation
  only differentially (≤ 1.5 pp, below the tolerance). Its independent coverage is withdrawn. The extrapolation stays
  declared, theory-bounded (I11-ext) and carried in the ±5 %; the absolute θ_ionmc/θ_SS (expected reference factor
  +3.1 to +4.4 %) is report-only.
- **Lateral gap.** The planner recommends filing the intervention at the freeze whatever the V5-SS outcome, with a
  ledger statement on the V5-SS class and on the f_dM extrapolation (P10 §7.6.2, arguments for and against).
  - The verbatim wording's clause "or an accessible single-scattering MC reference" is literally met by V5-SS, so
    whether it needs a reviewed amendment is an orchestrator decision before freeze.

**Dated note 2026-11-01 (round V3-008-5, Codex REVIEW-6a14f4ec finding 2; before any V5-SS run; SCRIPT-GENERATED numbers from `validation/scripts/ion_em/ss_budget.py` §(h), SEEN planning values).**
- **Rule.** The round V3-008-4 rule (PASS iff |dev| ≤ 3 % whenever σ ≤ 1 %) ignored σ and is withdrawn. V5-SS uses the
  plan's scalar TOST form per gated ratio: **PASS iff |ratio_model/ratio_ref − 1| + 1.645 σ ≤ 3 %**, σ the bootstrap
  standard error of the ratio statistic; σ > 1 % stays a separate precondition (INCONCLUSIVE, counted as FAIL).
- **Statistics (module).** Relative SE of the windowed rms from the window moments of the reference density: kurtosis in
  the ±2.5 θ_H window 2.76–2.78, x/y correlation 0.04–0.05; at 2e4 primaries 0.34 % per case (x and y pooled; 0.48 % with
  one projection); ionmc Gaussian at 1e5 histories 0.14 %. σ of the ratio statistic at the round V3-008-4 design
  (2e4 / 1e5): 0.53 %.
- **Power (true deviation per ratio = the computed bias −0.17 to −0.69 %, five correlated ratios, Monte Carlo 4e5
  draws).** At 2e4 / 1e5: P(all five pass) = 0.994 for a correct model, 0.982 with the bias shifted by the 0.2 pp
  approximation allowance, but 0.254 if ionmc's ion/ion ratio deviates from Molière–Hanson by the SEEN I11-rel spread
  (1.13 pp, adverse sign) — poor robustness.
- **Design decision (pre-registered now, before any V5-SS run): raise the statistics.** TOPAS 6e4 primaries per case
  (≤ 2 parts of ≤ 1200 s per case, parts 0–1 of the existing seed formula), ionmc ≥ 4e5 histories per case. σ of the
  ratio 0.30 %. P(all five pass) ≥ 0.9999 (correct model), 1.000 (+0.2 pp), 0.918 (SEEN model spread, adverse); a model
  term at the I11-rel edge (2 pp) fails (P = 0.006). REF cost scaled from the measured 1 mm C 290 smoke (23 ms per
  primary): 10 596 s ≈ 2.94 h ≤ 4 h (upper estimate). If the measured cost of this design exceeds 4 h, N falls back to
  4e4 (1.96 h; σ 0.36 %; 1.000 / 0.751) before the r6 drop order applies. The ionmc cost of 8 × 4e5 thin-slab histories is
  not measured (expected minutes).
- **Tolerance 3 % kept; five gated ratios kept.** Under TOST the statistics enter through 1.645 σ (0.49 % at the
  design), so the tolerance has to cover only the systematic items (bias ≤ 0.9 %, linear, plus ionmc's model term):
  the room left for ionmc's model term is ≈ 1.6 pp (worst-case bias) to 1.8 pp (computed bias), above the SEEN spread
  1.13 pp but below the I11-rel bound 2 pp. V5-SS is therefore deliberately tighter than the I11-rel band; a FAIL from a
  model term near 2 pp is a finding, not a false fail. 3.5 % was not adopted: it is not required by the budget at the
  pre-registered N and would let a model without the Born term pass He/C 250 (−2.66 % + 0.49 = 3.15 %).

**Lateral evidence gap.** It is used only after P6 (the measured lateral-data search) fails and V5-SS is withdrawn (its run condition fails). **r6:** the P6 lateral-data search is recorded as EVIDENCE GAP (`validation/plans/v3-008/p6-lateral-data-search.md`, committed at b3eef90). Because V5-SS shares the screened-Coulomb physics (above), no independent *physical* evidence of the charge dependence of screening can come from it. The planner recommends filing the intervention below before integration whatever the V5-SS outcome, with the V5-SS result and class stated in the ledger. This is an orchestrator decision before freeze; the intervention wording stays verbatim. **Orchestrator decision 2026-11-01 (before freeze):** (a) the lateral-gap intervention is filed at the freeze whatever V5-SS shows, naming as the unavailable input the independent evidence of the charge dependence of lateral spreading AND the absence of independent coverage of the f_dM extrapolation beyond the 300 MeV proton domain (declared and bounded only within the same Molière theory); V5-SS stays in the plan as shared-lineage single-scattering evidence (branch B) and its result and class go into the ledger; (b) the verbatim clause "or an accessible single-scattering MC reference" is amended in the reviewed freeze revision so that a reference sharing the screened-Coulomb cross-section lineage does not satisfy it; (c) the f_dM extrapolation is kept as a declared limitation (uncertainty component ≤ 1.5 %, I11-ext report row, `fdm_extrapolated` flag); no further coverage is claimed.
Intervention wording, verbatim:

> Unavailable input: independent evidence of the charge dependence of lateral spreading for He/C/O
> in water at 50–430 MeV/u, i.e. measured lateral dose widths or angular spreads with stated
> uncertainties (candidates Parodi 2013, Schwaab 2011, Tessonnier 2017: no open numeric tables
> found), or an accessible single-scattering MC reference. Request: provide access to such data, or
> accept that V3-008 integrates with ion lateral spreading supported only by theory-derived mapping
> and implementation checks (no independent evidence), recorded as a gap in the ledger.

**Revision 1 amendment (2026-11-01; orchestrator decisions (a) and (b) above; recorded before any V3-008 transport or
qualification result).** The wording above stays as the r5 record. Its clause, quoted verbatim:
"or an accessible single-scattering MC reference". Amended wording: "or an accessible single-scattering MC reference whose screened-Coulomb cross-section physics is independent of the model used by the exact-Ω mapping; a reference that shares the screened Rutherford law and Molière's Born screening form (including the Geant4 single-scattering configuration of V5-SS, `G4IonCoulombScatteringModel`, branch B) does not satisfy this input". The orchestrator files the intervention at
the freeze, whatever V5-SS shows, category inaccessible-information, with this wording (plan section "Intervention at
the freeze (to be filed by the orchestrator)"):

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

The sentence "It is used only after P6 … fails and V5-SS is withdrawn" at the head of this block is superseded by
decision (a). V5-SS stays a gating row (shared-lineage single-scattering evidence, branch B); its result and class go
into the ledger when it runs.

**Gating MCS rows:** I9, I10, I11-map, I11-rel, I12 and I15 (self-consistency, theory consistency
and falsification), and V5-SS (held-out, single scattering) if its run condition holds.

## 4. Kinematics, charge state, end of range, secondaries

- **Masses:** nuclear masses from `projectiles.py`. API energies are total kinetic energy in MeV
  (decision 0037); per-nucleon input only through an explicit helper, recorded as such.
- **Charge:** the bare nucleus z for MCS and straggling; the effective charge only in the stopping
  power.
- **End of range:** `E_cut = 2 MeV/u × A` (new field `e_cut_ion_mev_per_u = 2.0`, recorded), with
  local deposition below it; the table floor 1 MeV/u ≤ E_cut/2 rule is kept. Justification:
  - Residual CSDA ranges at 2 MeV/u are about 0.075 mm (He), 0.025–0.03 mm (C) and 0.02–0.025 mm
    (O) (estimates by velocity scaling).
  - As a fraction of the range at the lowest domain energy that is ≤ 0.34 % (He 50 MeV/u), 0.17 %
    (C 80) and 0.13 % (O 100).
  - Energy is conserved; only the last ≤ 0.075 mm of dose moves upstream by at most that distance.
    The effect on R80 is ≪ 0.01 mm against a straggling σ of 0.1–1 mm.
  - The estimates are reproduced from the tables in rows I6/C-res.
- **Deltas:** deposited locally (unrestricted S).
- **Nuclear interactions, fragments and secondary ions:** out of scope (V3-009). `nuclear=True` with
  z ≥ 2 raises `UnsupportedCombinationError`. The producible set is `{(ion, "primary")}`.
- **Species:** append `carbon12` (id 5) and `oxygen16` (id 6) after `alpha` (id 4), append-only.
  Merge risk with PR #82: TO VERIFY that it appends no species ids.
- **LET:** decision 0040 already defines S_w as the water stopping of the transported ion. The ion's
  own water row is built for primary-only fields; mixed-field LET is V3-009.

## 5. Supported domains (no clinical claim)

| Ion | Source energy (MeV/u) | Bounds |
|---|---|---|
| 4He | 50–250 | upper: end of the ASTAR/ICRU 90 alpha tables; lower: about 2 cm range |
| 12C | 80–430 | range about 2–30 cm water |
| 16O | 100–430 | the ledger's declared S-PHYS-OXYGEN entrance domain |

- Evidence is for water. Other `materials.py` materials go through the same analytic model with
  theory rows only, and are documented as such.
- Fail closed for: energies outside the domain (including E0 + 6σ_E); He-3 and other ions;
  `nist-star` for C/O; nuclear on.

**Dated note 2026-11-01 (§3 implementation).** The exact-Ω mapping is implemented in `src/ionmc/physics/exact_omega.py`
as the declared per-step root solve (bracketed Illinois on ln pv_eq, ≤ 19 iterations observed, float64 b residual
≤ 2e-12), with a dedicated proton path to 4000 MeV, fail-closed status codes, the `fdm_extrapolated` flag and the
source proton-equivalent energy defined as the maximum over a depth scan (plan I0/I11-map notes; P11d). The requested
and effective LS modes are recorded in the effective configuration outside the table hash (Codex REVIEW-68f985dd minor).
The per-step cost at transport scale is not yet measured.

## 6. Data roles and lineage (r4)

| Item | Role | Shared with Geant4/TOPAS? |
|---|---|---|
| Bethe model and corrections (ionmc) | construction | Bethe theory shared; the Barkas F table is transcribed from `G4EmCorrections` (shared code lineage; TOPAS adds its own Barkas term for ions) |
| I_water 78 eV (ICRU 90) | construction | shared (G4_WATER). ICRU 90's value was informed by measured Bragg-peak positions (GSI, Kumazaki 2007, Bichsel & Hiraoka 1992; TO VERIFY), so measured water ranges are calibration lineage |
| ASTAR (He) | evaluation | not used by TOPAS |
| ICRU 90 alpha array (He) | evaluation (I2: a regression guard on seen values, not evidence) | **identical to TOPAS alpha stopping**: FACT (research note); P2 OPEN pending `G4EmBuilder.cc`/`G4BetheBlochModel.cc`. For helium, V1 tests transport mechanics and straggling only |
| ICRU 73 (rev.) C/O: libdEdx ICRU73_NEW (I3 evaluation source); G4EMLOW icru73/icru90 files lineage only (Geant4-restricted use) | evaluation (libdEdx; exploratory until I3 frozen) / lineage (G4EMLOW) | the G4EMLOW files are consulted by TOPAS only at ≤ 2.0 MeV scaled (≈ 1.98 MeV/u), 2 < Z ≤ 80, never alpha (P2 FACT) → the libdEdx comparison is largely independent of TOPAS above |
| TOPAS C/O stopping: `G4BetheBlochModel` table (ICRU 90 proton array at equal velocity: research note, partly confirmed) + q_eff²(ΔL_LS + 2·Barkas·(z−1)/z), no separate Bloch/Mott, switch at 2.0 MeV scaled (P2 FACT) | the reference's construction; I3c attribution diagnostic | — (the r3 trigger "ion tables above 2.5 MeV/u" is not met; no amendment) |
| Molière–Hanson ion theory (numpy, P1) | theory for the relative ion/proton b mapping (construction); I11-rel gating (relative), I11-abs informative | Urban MSC differs structurally: no charge/β term in its log, truncated core + tail, step-length dependent (θ0 form P2 FACT; Urban as the ion msc model P2 INFERRED, settled by P2 follow-up/P4) |
| Gottschalk f_dM | construction (qualified for protons ≤ 300 MeV; evaluated at proton-equivalent states up to ≈ 2.5 GeV, extrapolation declared and bounded, §3 r6 note). **2026-11-01 (round V3-008-4, finding 3):** the extrapolation has no independent coverage. It is bounded only within the same Molière theory (I11-ext, ≤ 1.5 % component of the ±5 %), V5-SS covers it differentially only, and the gap is named in the lateral-gap intervention filed at the freeze (orchestrator decision 2026-11-01 (a), (c)) | no |
| Geant4 11.4.2 single scattering (`g4em-standard_SS`; V5-SS) | evaluation (held-out, conditional) | shares the screened-Coulomb cross-section physics (Born screening form; proton path FACT, ion path OPEN, P10 §7); independent of f_dM and of the multiple-scattering approximation. **2026-11-01 (round V3-008-4, finding 3):** ion path FACT (P10 §7.4.1), branch B: ion/ion ratios at a common velocity on `G4IonCoulombScatteringModel` (alpha and GenericIon); shares the screened Rutherford law and the Born form (with corr^0.04), differs in the universal z-dependent radius and the absent electrons (pre-computed bias, `validation/scripts/ion_em/ss_budget.py`, P10 §7.8). Shared-lineage evidence: it does not satisfy the lateral-gap input (orchestrator decision (b)) and covers the f_dM extrapolation only differentially |
| TOPAS EM-only ion runs | evaluation, independent MC (V1 gating; V2 informative by default) | stopping lineage partly shared (I-value, ICRU 90 proton/alpha arrays, Barkas); straggling independent; MSC model structurally different (Urban) |
| Sihver, Schardt, Kanai 1998 measured C/O depth-dose | **FAILED as a reference (P6):** curves normalised horizontally to the calculated peaks, figures only, oxygen = 18O | → V3 is an EVIDENCE GAP (intervention path, bundled with the lateral gap if V5-SS is withdrawn) |

## Expected tradeoffs and declared limitations

- **Helium** range evidence is tabulated (SEEN) plus TOPAS, which shares the ICRU 90 alpha array.
  **No measured He range evidence was found.**
- **Carbon and oxygen:** the tabulated check uses libdEdx ICRU73_NEW (P5b); the Geant4 G4EMLOW files are
  lineage only. Measured stopping at therapy energies does not exist publicly. The only measured-range
  candidate (Sihver 1998) FAILED P6, so measured range is an EVIDENCE GAP (intervention path).
- **MCS:** the ion model is the exact-Ω (b-matching) mapping of the proton model, with ±5 %
  theory-derived uncertainty. I11-rel (relative) gates. I11-abs is informative, because the water
  reference fails its Pb validation, is unvalidated for water and rests on a derived compound rule.
  V2 follows the pre-declared branches (Case C if Urban). Independent evidence for the z-dependence
  of lateral spreading rests on V5-SS (held-out, single scattering, run condition per P10) or measured data; otherwise it is a declared gap; the intervention wording is in
  §3. r6: measured data are an EVIDENCE GAP (P6). V5-SS shares the screened-Coulomb cross-section physics (P10 §7), and its registration waits for the ion-side Geant4 model. f_dM is evaluated beyond its proton-qualified domain at most ion energies: the extrapolation is declared and bounded, ≤ 1.5 % against Molière theory (§3 note of 2026-11-01).
  **Dated note 2026-11-01 (round V3-008-4, Codex REVIEW-e0b586fd finding 3; supersedes the V5-SS clauses of this bullet).** The ion-side model is read and V5-SS is branch B (ion/ion ratios at a common velocity, shared screened-Coulomb lineage). V5-SS is therefore not independent evidence for the z-dependence of lateral spreading, and it covers the f_dM extrapolation only differentially. Per the orchestrator decision of 2026-11-01 (§3), the lateral-gap intervention is filed at the freeze whatever V5-SS shows, naming both the missing independent evidence of the charge dependence of lateral spreading and the absence of independent coverage of the f_dM extrapolation above the 300 MeV proton domain. The extrapolation stays a declared limitation (≤ 1.5 % component, I11-ext, `fdm_extrapolated`); no further coverage is claimed.
- **Gaps** go to the intervention route. Nothing is relabelled, and no model is tuned to TOPAS.

## Validation strategy

Frozen in `validation/plans/v3-008-acceptance.md` (draft r6 alongside).

**Revision 1 (2026-11-01):** frozen in `validation/plans/v3-008-acceptance.md` revision 1, by the reviewed freeze
commit named in `validation/plans/v3-008/p9-freeze-record.md`, once FC-1 and FC-2 are disposed.
