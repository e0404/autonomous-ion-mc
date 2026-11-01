# 0043 — Helium, carbon and oxygen electromagnetic transport (DRAFT r3 for independent review)

> **Status: DRAFT r3 (2026-11-01), not yet frozen.** Written by the scientific planner, revised against the two research notes in `validation/plans/v3-008/research/`, independently reviewed in round V3-008-1 (REVISE; the nine required wording changes applied; `validation/plans/v3-008/review.md`). The decision is frozen by a reviewed commit only after the pre-freeze checks P1–P9 of the acceptance plan resolve the items marked OPEN / TO-VERIFY / "confirmed by P2". No result has been produced; nothing has run. Operator directive of 2026-11-01 and the reviewed selection `validation/plans/v3-next/selection-2026-11-01.md` (on the V3-005B branch) are the basis.


- Status: DRAFT r3 (planner, 2026-11-01). r2 was revised against `research-models.md` and
  `research-data-sources.md` in the same scratch directory; r3 applies the wording changes of review
  round V3-008-1. Not accepted. To be reviewed
  independently and committed together with `validation/plans/v3-008-acceptance.md` before any
  V3-008 result exists. The pre-freeze checks of the plan (section P) must be closed and their
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

**What Geant4 11.4.2 `g4em-standard_opt4` (TOPAS) does for ions** (FACT (research note), confirmed by P2):
- **GenericIon (C, O):**
  - The model is `G4LindhardSorensenIonModel`.
  - Above about 2 MeV/u (scaled energy E·m_p/M > 2 MeV), dE/dx is the ICRU 90 **proton** water
    array at equal velocity × Geant4 effective charge squared. This is `G4BetheBlochModel` with
    `UseICRU90Data`.
  - Along-step it adds the Lindhard–Sørensen ΔL (`G4LindhardSorensenData`) and a Barkas term
    `2·BarkasCorrection·(q−1)/q`.
  - ICRU 73/90 ion tables (`G4IonICRU73Data`, G4EMLOW `ion_stopping_data`) are read only for scaled
    energies ≤ 2.5 MeV, i.e. the last tens of µm of range. P2 must confirm the energy limit at which
    `G4LindhardSorensenIonModel` consults `G4IonICRU73Data`. Pre-declared: if the ion tables are used
    above 2.5 MeV/u, I3's class becomes "shared with TOPAS" and V1's C/O lineage changes (amendment
    before freeze).
- **Alpha:** the ICRU 90 alpha water array up to 1 GeV total.

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
- **Effective charge:** the Barkas form with coefficient 125, in the stopping power only.
  - Below 1 MeV/u no stopping power is evaluated, and transport deposits locally below
    E_cut = 2 MeV/u (section 4).
  - Low-velocity alternatives (Hubert–Bimbot–Gauvin 1990, Ziegler/Brandt–Kitagawa) would differ by
    a few % in q² at 2–5 MeV/u and by < 0.5 % above 10 MeV/u (TO VERIFY). That moves the end of
    range by ≤ 0.01–0.02 mm (estimate), below the 0.05 mm action threshold. No change; the bound is
    documented.

### Uncertainty statement (metadata and docs; not acceptance tolerances)

- I_water 78 ± 2 eV (ICRU 90) gives about ±0.3 % of range.
- Omitted LS term: ≤ 0.1 % of S (estimate; P3).
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
| 12C, 16O | ICRU 73 water tables, if they carry the 2009 erratum (I = 78 eV): G4EMLOW 8.8 `ion_stopping_data/icru73/z6_G4_WATER.dat`, `z8_G4_WATER.dat` | OPEN, pre-freeze checks P5: file span must cover ≥ 2–430 MeV/u; erratum vs original-ICRU-73 provenance (the original water tables used I = 67.2 eV; the difference is a few %); licence = Geant4 Software Licence | ICRU 73 (PASS/BEST theory). In TOPAS these files act only ≤ 2.5 MeV/u, so above that the row is **largely independent of TOPAS's stopping** (shared: Bethe theory, I-value) |
| 12C (16O?) | ICRU 90 carbon water table: G4EMLOW `icru90/z6_G4_WATER.dat` | OPEN P5: whether it holds ICRU 90 published carbon values or a Geant4 recomputation. **No ICRU 90 oxygen table exists** (ICRU 90 covers p, α, C), so `z8` is of unclear origin and not used | ICRU 90 |
| 12C, 16O | libdEdx data files (ICRU 73 with erratum, per memory) | OPEN P5b: location, licence (GPL), file provenance | cross-check of the G4EMLOW files |
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
    independent of TOPAS**, because Geant4 uses them only below about 2.5 MeV/u. Making them the
    default would turn that independent check into a shared-construction comparison.
  - The r1 reason ("TOPAS uses these data") was wrong above about 2 MeV/u and is withdrawn.
- **The Geant4 construction as default** (ICRU 90 proton array at equal velocity × Geant4 q_eff² +
  LS ΔL + Barkas). It would make ionmc's C/O stopping identical in construction to TOPAS's, so V1
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

All Molière–Bethe relations as cited by Gottschalk and Lynch & Dahl (constants TO VERIFY against
Lynch & Dahl, NIM B 58 (1991) 6):
- χ_c² ∝ z² Σ w_i Z_i(Z_i+1)/A_i · t/(pv)².
- χ_a² = χ_0² (1.13 + 3.76 α'²), with χ_0 ∝ Z^(1/3)/p and the Born parameter α' = z Z α/β.
- Collision number Ω = χ_c²/(1.167 χ_a²) ∝ z² t/[β² (1.13 + 3.76 α'²)].
- Hanson width θ_0 = χ_c √(B − 1.2)/√2, with B − ln B = ln Ω.

Consequence at equal velocity and equal x/R:
- The ion thickness is t_ion = (A/z²) t_p.
- Hence Ω_ion/Ω_p = A/F_scr. F_scr is the Z(Z+1)-weighted ratio of the Born screening factors of
  ion and proton over the target elements (oxygen dominates water: 72 of 76).
- Estimated F_scr in water: about 1.15 (He, 150 MeV/u), 1.9 (C, 290), 2.5 (O, 350).
- The ion Molière width therefore exceeds pure z/(pv) proton scaling by about **+5 % (He) and +8 %
  (C, O)** (estimate).

### Selected: candidate (d), velocity-scaled differential Molière with a Molière collision-number correction

```
T_ion = f_dM,ion · (z E_s / pv)² / X_S
f_dM,ion = f_dM(pv*, p1v1*) + c_Ω · ln( A / F_scr(β, z; material) )
pv* = pv · m_p / M_ion,   p1v1* = p1v1 · m_p / M_ion
```

- `pv` is the ion's own momentum × velocity in the Rutherford prefactor (exact χ_c² ∝ z²/(pv)²).
- The memory term `1 − (pv/p1v1)²` is invariant under the velocity scaling. So the random hinge,
  the clamp, the birth-step average and U5-type step independence carry over by construction: the
  correction is a local function of (β, material).
- `c_Ω` is the sensitivity of f_dM to ln Ω at fixed x/R, derived **analytically from the f_dM proton
  fit and Molière theory**: the lg pv coefficient 0.232 together with Ω ∝ (pv)^≈0.77 at fixed x/R,
  so c_Ω ≈ 0.13 (estimate).
- F_scr comes from the Born factors above.
- Both are committed as a formula with its derivation (pre-freeze check P1). **c_Ω and F_scr are
  frozen as committed formulas before I11 is evaluated; no adjustment after.** Their hashes enter the
  P9 freeze record. **Nothing is fitted to TOPAS or to any transport output.**
- Implementation: the per-material correction is precomputed on the transport energy grid (a table
  row in `TransportTables`, inside the hash) and read like S(E).
- **Pre-declared selection order** after the I11 reference exists: (d), (c), (b), (a). The first
  candidate with max |θ/θ_MH − 1| ≤ 3 % at the I11 grid points for all three ions is the default.
  If none passes, STOP and amend.
- Expected outcomes, declared now (estimates):

  | Candidate | Expected θ vs Molière–Hanson |
  |---|---|
  | (d) | ≤ 1–2 % |
  | (c) Highland/Lynch–Dahl with z² in the logarithm, differential form (Kanematsu 2008, NIM B 266, 5056; ion form TO VERIFY) | +1 to +3 % |
  | (b) naive T_dM with the ion pv | +2 to +4 % |
  | (a) velocity scaling only (the r1 choice) | −5 % (He), −7 to −8 % (C/O): expected to FAIL |

- **Why the tolerance is 3 % on θ** (projected rms):
  - f_dM reproduces θ_Hanson for protons within about ±2 % (Gottschalk; exact band TO VERIFY).
  - Ion-specific theory uncertainty is about 1–2 %, from the Born-screening form at α' ≈ 0.5–0.7,
    the Fano electron term and the nuclear form factor (core-width effect ≤ 1–2 %, estimate).
  - In quadrature that is ≈ 2.5–3 %; tightening below 3 % is not defensible.
- **Role:** the Molière–Hanson ion reference is **calibration/construction** for this choice. Row
  I11 is declared calibration: for candidate (d), whose correction is derived from the same Molière
  theory, it is close to a tautology. It checks the implementation and the derivation, not the physics,
  and it is not independent evidence.
- **Pre-declared Molière–Hanson reference definition** (P1):
  - χ_c with Z(Z+1) and the Fano electron term with its β-dependence (as Gottschalk 1993);
  - χ_a per element with 1.13 + 3.76 α'²;
  - ln χ_a averaged along the slab with weight 1/(pv)²;
  - point nucleus, with the nuclear-size effect reported separately.
  - The code is validated first against Gottschalk's proton θ_Hanson tables (Be/Al/Cu/Pb, already
    acquired).
- **Rejected:**
  - per-step Highland: step-size dependent;
  - full Molière angular sampling with single-scattering tails: a large kernel change, and the
    Gaussian-core width is the target; the tail deficit stays a documented limitation as for protons
    (T15 window rule ¹⁰);
  - (a) as default: it misses the collision-number change, which is physics, not a fit artefact.

### Independent evidence for MCS and the TOPAS reference model

- TOPAS opt4 multiple scattering for GenericIon/alpha is believed to be `G4UrbanMscModel`. Its
  `ComputeTheta0` is believed to be a Highland-like `13.6·|q|·√y/(βcp)·(c1 + c2 ln y)` **without
  z² in the logarithm** (TO VERIFY: pre-freeze check P2 reads the v11.4.2 source and the opt4 msc
  assignment).
- If confirmed, TOPAS ion widths are expected to be about 5–8 % (C/O) and about 5 % (He) narrower
  than a Molière-consistent ionmc (estimate). V2 could then fail on the reference model's account.
- The method is checked on protons first (P2b): the same δ_pred computation for z = 1 is compared
  with the archived V3-003 proton TOPAS lateral result (T15, already seen), and the discrepancy e_δ is
  recorded. Case B is available only if e_δ ≤ 2 %. P2 freezes the exact Urban formula, version and
  parameters with the source hashes.
- The plan therefore **pre-declares δ_pred**, the expected TOPAS/ionmc lateral-σ offset per ion and
  depth. It is computed before any TOPAS run, from the Urban formula as read in the source and from
  the I11 reference, by a theory-to-theory numpy Fermi–Eyges integration. V2 grades ionmc against
  TOPAS corrected by δ_pred and reports the raw ratio alongside.
- V2 then asserts that the transport implements the selected MCS model consistently with an
  independent MC code once the known, source-derived model difference is removed. **It does not
  independently test the z-dependence**, which rests on Molière theory (calibration role). V2 has no
  held-out analogue: a model change after V2 leaves lateral MC evidence as a declared gap.
- Measured lateral data for ions in water (Parodi et al. 2013 J. Radiat. Res.; Tessonnier 2017;
  Schwaab 2011) are the only independent arbiter. No public numeric set is confirmed, so this is a
  declared evidence gap and intervention candidate. Wong et al. 1990 (U on Cu) is not applicable.

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

## 6. Data roles and lineage (r3)

| Item | Role | Shared with Geant4/TOPAS? |
|---|---|---|
| Bethe model and corrections (ionmc) | construction | Bethe theory shared; the Barkas F table is transcribed from `G4EmCorrections` (shared code lineage; TOPAS adds its own Barkas term for ions) |
| I_water 78 eV (ICRU 90) | construction | shared (G4_WATER). ICRU 90's value was informed by measured Bragg-peak positions (GSI, Kumazaki 2007, Bichsel & Hiraoka 1992; TO VERIFY), so measured water ranges are calibration lineage |
| ASTAR (He) | evaluation | not used by TOPAS |
| ICRU 90 alpha array (He) | evaluation (I2: a regression guard on seen values, not evidence) | **identical to TOPAS alpha stopping** — FACT (research note), confirmed by P2; for helium, V1 tests transport mechanics and straggling only |
| ICRU 73 (rev.) / ICRU 90 C/O files (if P5 passes) | evaluation | used by TOPAS only ≤ 2.5 MeV/u → largely independent above |
| TOPAS C/O stopping (ICRU 90 proton × q_eff,G4² + LS ΔL + G4 Barkas) — FACT (research note), confirmed by P2 | the reference's construction; I3c attribution diagnostic | — (if P2 finds ion tables used above 2.5 MeV/u: I3 becomes shared, V1 C/O lineage amended) |
| Molière–Hanson ion theory (numpy) | calibration of the MCS choice (d) | Urban MSC differs (independent implementation, possibly different z-dependence; P2) |
| Gottschalk f_dM | construction | no |
| TOPAS EM-only ion runs | evaluation, independent MC | stopping lineage partly shared (I-value, ICRU 90 proton/alpha arrays, Barkas); straggling and MSC independent |
| Sihver, Schardt, Kanai 1998 measured C/O depth-dose (if acquirable) | evaluation, measured | I-value calibration lineage; pre-ICRU 90 |

## Expected tradeoffs and declared limitations

- **Helium** range evidence is tabulated (SEEN) plus TOPAS, which shares the ICRU 90 alpha array.
  **No measured He range evidence was found.**
- **Carbon and oxygen:** a tabulated check depends on P5. Measured stopping at therapy energies does
  not exist publicly. Measured range depends on one candidate (Sihver 1998); oxygen has no other.
- **MCS:** the z-dependence is theory-calibrated (d). The independent MC check corrects for a
  source-derived reference offset, and measured lateral data are an evidence gap.
- **Gaps** go to the intervention route. Nothing is relabelled, and no model is tuned to TOPAS.

## Validation strategy

Frozen in `validation/plans/v3-008-acceptance.md` (draft r3 alongside).
