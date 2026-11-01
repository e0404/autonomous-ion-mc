# 0043 — Helium, carbon and oxygen electromagnetic transport (DRAFT r5 for independent review)

> **Status: DRAFT r5 (2026-11-01), not yet frozen.** r3 was committed at 2632a465 and reviewed in round V3-008-1; r4/r5 fold in the pre-freeze records P8, P7, P2, P1/P1b (committed at 6e155b56) and the review rounds V3-008-1 and V3-008-2 (`validation/plans/v3-008/review.md`; the six V3-008-2 wording changes applied). The decision is frozen by a reviewed commit only after the pre-freeze checks of the acceptance plan (§P, steps 1–7) close. No transport or qualification result has been produced; the pre-freeze computations P1/P1b (reference module, candidate tables) and the P2 source reading exist and are disclosed.


- Status: DRAFT r5 (planner, 2026-11-01), not frozen. History:
  - r2 was revised against the two research notes;
  - r3 applied the wording changes of review round V3-008-1;
  - r4 folds in the committed pre-freeze records at 6e155b56: `validation/plans/v3-008/`
    `p1-moliere-hanson-reference.md`, `p2-geant4-ion-models.md` and `p8-repository-facts.md`.
  - Main r4 changes: §3 (MCS model and the status of I11/V2) and the Geant4 lineage wording.
  - r5 applies the six wording changes of review round V3-008-2: the §3 selection statement, the
    exact-Ω specification, gating I11-rel, the V2 branches, the lateral-gap intervention, and the
    §6 Molière–Hanson row. Not accepted. To be reviewed
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
| 12C, 16O | ICRU 73 water tables, if they carry the 2009 erratum (I = 78 eV): G4EMLOW 8.8 `ion_stopping_data/icru73/z6_G4_WATER.dat`, `z8_G4_WATER.dat` | OPEN, pre-freeze checks P5: file span must cover ≥ 2–430 MeV/u; erratum vs original-ICRU-73 provenance (the original water tables used I = 67.2 eV; the difference is a few %); licence = Geant4 Software Licence | ICRU 73 (PASS/BEST theory). In TOPAS these files are consulted only at ≤ 2.0 MeV scaled energy (≈ 1.98 MeV/u) and only for 2 < Z ≤ 80 (P2 FACT), so above that the row is **largely independent of TOPAS's stopping** (shared: Bethe theory, I-value) |
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

**Freezing.**
- The exact-Ω mapping, its b formulas, its evaluation choice (root solve or a bounded tabulation)
  and the I11-map and I11-rel tests are committed **before the mapping is evaluated against the
  Molière–Hanson reference**. Their hashes go into P9.
- Nothing is fitted to the reference, TOPAS or transport output.

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
If it is available within the REF budget, row V5-SS is pre-registered before any run.

**Lateral evidence gap.** It is used only after P6 (the measured lateral-data search) and P10 fail.
Intervention wording, verbatim:

> Unavailable input: independent evidence of the charge dependence of lateral spreading for He/C/O
> in water at 50–430 MeV/u, i.e. measured lateral dose widths or angular spreads with stated
> uncertainties (candidates Parodi 2013, Schwaab 2011, Tessonnier 2017: no open numeric tables
> found), or an accessible single-scattering MC reference. Request: provide access to such data, or
> accept that V3-008 integrates with ion lateral spreading supported only by theory-derived mapping
> and implementation checks (no independent evidence), recorded as a gap in the ledger.

**Gating MCS rows:** I9, I10, I11-map, I11-rel, I12 and I15 (self-consistency, theory consistency
and falsification), and V5-SS if registered.

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

## 6. Data roles and lineage (r4)

| Item | Role | Shared with Geant4/TOPAS? |
|---|---|---|
| Bethe model and corrections (ionmc) | construction | Bethe theory shared; the Barkas F table is transcribed from `G4EmCorrections` (shared code lineage; TOPAS adds its own Barkas term for ions) |
| I_water 78 eV (ICRU 90) | construction | shared (G4_WATER). ICRU 90's value was informed by measured Bragg-peak positions (GSI, Kumazaki 2007, Bichsel & Hiraoka 1992; TO VERIFY), so measured water ranges are calibration lineage |
| ASTAR (He) | evaluation | not used by TOPAS |
| ICRU 90 alpha array (He) | evaluation (I2: a regression guard on seen values, not evidence) | **identical to TOPAS alpha stopping**: FACT (research note); P2 OPEN pending `G4EmBuilder.cc`/`G4BetheBlochModel.cc`. For helium, V1 tests transport mechanics and straggling only |
| ICRU 73 (rev.) / ICRU 90 C/O files (if P5 passes) | evaluation | consulted by TOPAS only at ≤ 2.0 MeV scaled (≈ 1.98 MeV/u), 2 < Z ≤ 80, never alpha (P2 FACT) → largely independent above |
| TOPAS C/O stopping: `G4BetheBlochModel` table (ICRU 90 proton array at equal velocity: research note, partly confirmed) + q_eff²(ΔL_LS + 2·Barkas·(z−1)/z), no separate Bloch/Mott, switch at 2.0 MeV scaled (P2 FACT) | the reference's construction; I3c attribution diagnostic | — (the r3 trigger "ion tables above 2.5 MeV/u" is not met; no amendment) |
| Molière–Hanson ion theory (numpy, P1) | theory for the relative ion/proton b mapping (construction); I11-rel gating (relative), I11-abs informative | Urban MSC differs structurally: no charge/β term in its log, truncated core + tail, step-length dependent (θ0 form P2 FACT; Urban as the ion msc model P2 INFERRED, settled by P2 follow-up/P4) |
| Gottschalk f_dM | construction | no |
| TOPAS EM-only ion runs | evaluation, independent MC (V1 gating; V2 informative by default) | stopping lineage partly shared (I-value, ICRU 90 proton/alpha arrays, Barkas); straggling independent; MSC model structurally different (Urban) |
| Sihver, Schardt, Kanai 1998 measured C/O depth-dose (if acquirable) | evaluation, measured | I-value calibration lineage; pre-ICRU 90 |

## Expected tradeoffs and declared limitations

- **Helium** range evidence is tabulated (SEEN) plus TOPAS, which shares the ICRU 90 alpha array.
  **No measured He range evidence was found.**
- **Carbon and oxygen:** a tabulated check depends on P5. Measured stopping at therapy energies does
  not exist publicly. Measured range depends on one candidate (Sihver 1998); oxygen has no other.
- **MCS:** the ion model is the exact-Ω (b-matching) mapping of the proton model, with ±5 %
  theory-derived uncertainty. I11-rel (relative) gates. I11-abs is informative, because the water
  reference fails its Pb validation, is unvalidated for water and rests on a derived compound rule.
  V2 follows the pre-declared branches (Case C if Urban). Independent evidence for the z-dependence
  of lateral spreading is a declared gap, unless P6 or P10 provide it; the intervention wording is in
  §3.
- **Gaps** go to the intervention route. Nothing is relabelled, and no model is tuned to TOPAS.

## Validation strategy

Frozen in `validation/plans/v3-008-acceptance.md` (draft r4 alongside).
