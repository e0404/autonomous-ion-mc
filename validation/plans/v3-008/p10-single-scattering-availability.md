# P10: availability and cost of a Geant4 single-scattering thin-slab reference (pre-freeze check, V3-008)

Status: pre-freeze check P10 of `validation/plans/v3-008-acceptance.md` (DRAFT r5). Written
2026-11-01 by the research specialist at worktree head `3d2547ad`. Documentation and source
reading plus hand arithmetic only. No TOPAS run, no smoke and no transport output. Nothing here is
an acceptance value. The plan's §P requires a scorer-free smoke for the cost figure; it has **not**
been run, so P10 is **not closed**.

Markers:
- **FACT**: read in a pinned Geant4 file (sha256 in `p2-geant4-ion-models.md` §0 and §F0), quoted
  with line numbers.
- **LIT**: public documentation or a secondary source, cited with URL; not version-pinned to the
  engine.
- **INFERRED**: arithmetic, or recalled code that is not pinned.
- **OPEN**: the settling file or run is named.

## 1. Is a single-scattering EM configuration available in OpenTOPAS 4.3?

- **LIT:** the OpenTOPAS module table ("The table lists the aliases supplied directly by
  OpenTOPAS") lists `g4em-standard_SS | G4EmStandardPhysicsSS`, next to `g4em-standard_WVI`,
  `g4em-standard_GS` and `g4em-standard_opt0…opt4`
  (<https://opentopas.readthedocs.io/en/main/parameters/physics/modular.html>, read 2026-11-01).
  - The page documents v4.3.0 features: case-insensitive lookup, and Geant4 registry names
    accepted. It does not state its version.
  - So the alias's presence in the installed 4.3 build is LIT, not FACT. It is confirmed by a
    smoke run (section 4) that the engine does not reject.
- **FACT, Geant4 11.4.2 [EMB] `G4EmBuilder.cc`:** a single-scattering builder exists and gives
  alpha single Coulomb scattering with no msc process. [EMB] L146–165 (`ConstructIonEmPhysicsSS`):
  ```
  part = G4Alpha::Alpha();
  ph->RegisterProcess(new G4ionIonisation(), part);
  ph->RegisterProcess(new G4CoulombScattering(false), part);
  ```
  - d, t and He3 are the same (L150–164).
  - It is reached from `ConstructChargedSS` (L301–343, L343).
  - **GenericIon is not configured in `G4EmBuilder`'s SS path** (L146–165, L301–358). It must be
    configured in `G4EmStandardPhysicsSS::ConstructProcess` itself.
- **LIT (version-mismatched):** a Geant4 10.4 source listing of `G4EmStandardPhysicsSS.cc` (found
  by web search, not fetched) registers `G4ionIonisation` + `G4CoulombScattering` for alpha, He3
  and GenericIon (<https://nuclear.korea.ac.kr/~ejungwoo/geant4/G4EmStandardPhysicsSS_8cc_source.html>).
  A Lund Geant4 EM lecture describes the SS list as single elastic scattering for all particles,
  for validation and verification
  (<https://indico.lucas.lu.se/event/932/contributions/2744/attachments/718/1372/EMPhysics1.pdf>).
- **OPEN (11.4.2):**
  - whether `G4EmStandardPhysicsSS` calls `ConstructChargedSS`;
  - what it registers for GenericIon (process and scattering model);
  - any EM parameters it sets.
  - **Acquire:**
    - `source/physics_lists/constructors/electromagnetic/src/G4EmStandardPhysicsSS.cc` (+ `include/G4EmStandardPhysicsSS.hh`);
    - for the model physics, `source/processes/electromagnetic/standard/src/G4CoulombScattering.cc` (default model, `isCombined = false` meaning, table flag for nuclei), `.../standard/src/G4eCoulombScatteringModel.cc` and `.../standard/src/G4WentzelOKandVIxSection.cc` (screening, and whether the projectile charge enters the Molière α′ term; nuclear form factor; treatment of atomic electrons versus the ionisation cut);
    - `.../standard/src/G4IonCoulombScatteringModel.cc`, if the SS constructor uses it for GenericIon.
- INFERRED (recalled, not pinned):
  - `G4CoulombScattering(false)` means pure single scattering (no msc combination, no lower angular
    limit);
  - the default model is `G4eCoulombScatteringModel` (Wentzel screened Rutherford with a nuclear
    form factor);
  - for particle type "nucleus" no cross-section table is built (on-the-fly cross sections, a cost
    driver).
- Alternative if SS proves too expensive: `g4em-standard_WVI` (`G4EmStandardPhysicsWVI`; LIT
  alias). That is WentzelVI msc plus single scattering above a θ limit. It is cheaper but **not**
  independent of a condensed msc model; its ion assignment is OPEN
  (`.../constructors/electromagnetic/src/G4EmStandardPhysicsWVI.cc`). It is not recommended as a
  V5-SS substitute without an amendment.

**Availability verdict:**
- alpha: available in Geant4 11.4.2 (FACT, conditional on the SS constructor calling
  `ConstructChargedSS`);
- GenericIon: LIT (10.4) / OPEN (11.4.2 file above);
- the OpenTOPAS alias: LIT.
- All three are settled together by acquiring the SS constructor file and by one scorer-free smoke
  whose process list shows `CoulombScat` (and no `msc`/`ionmsc`) for GenericIon and alpha.

## 2. Independence and lineage (to be declared in V5-SS)

- Independent of Urban, of the Highland-type θ0 and of Gottschalk's f_dM: the angular distribution
  is built from individual elastic events.
- **Not fully independent of the screening theory.** If `G4WentzelOKandVIxSection` uses a
  Molière-type screening with the projectile charge in α′ (INFERRED, OPEN, file above), it shares
  the screening ingredient (1.13 + 3.76 α′²) with the exact-Ω mapping's Born factor (decision 0043
  §3). V5-SS then tests the angular integration and the mapping's realisation of b(z, M), not an
  independent determination of the z-dependence of screening. The lineage must be stated in the
  row.
- The treatment of atomic electrons matters at the few-% level. Scratch: the electron share of
  χ_c² in water ≈ 13 %, i.e. ≈ 6 % of θ. Geant4 splits electron scattering between the Coulomb
  model (below the ionisation cut) and delta production (above it; whether the primary is deflected
  there is OPEN, `G4BetheBlochModel.cc` `SampleSecondaries`). The electron production cut is
  therefore part of the pre-registered configuration. Recommended: the matched large cut, so that
  all electron scattering stays inside the Coulomb model, with the same cut for the ion and the
  proton case.
- **Recommended observable form:** the relative ratio
  (θ_ion,SS/θ_p,SS) / (θ_ion,ionmc/θ_p,ionmc) at equal velocity and slab thickness, plus the
  absolute ratio θ_ionmc/θ_SS as a second statistic. Rationale: the SS model's own absolute offsets
  for water (screening constant, electron term) largely cancel in the ion/proton ratio, as in
  I11-rel. The proton SS case is cheap (section 3). Which statistic gates is the orchestrator's
  decision at pre-registration.

## 3. Cost estimate (INFERRED; to be replaced by the smoke measurement)

**Events per history.**
- Single scattering samples about Ω elastic events (Molière collision number) per slab.
- Scratch computation (Bethe/Gottschalk constants, constant energy, water; script in the session
  scratchpad, not committed):

| Case | Ω per 1 mm | per 10 mm | per 50 mm |
|---|---|---|---|
| p 150 MeV | 2.8e3 | 2.8e4 | 1.4e5 |
| He 150 MeV/u | 1.0e4 | 1.0e5 | 5.1e5 |
| C 290 MeV/u | 3.4e4 | 3.4e5 | 1.7e6 |
| O 350 MeV/u | 4.3e4 | 4.3e5 | 2.1e6 |

- Ω_ion/Ω_p at equal velocity and thickness ≈ z²/F_scr (≈ 19 for C), so ion SS is 4–15× costlier
  per mm than proton SS.

**Cost per step.**
- The project's measured TOPAS EM-only proton throughput is ≈ 6 min per 1e5 histories, i.e.
  ≈ 3.6 ms/history: 150 MeV, full range, two 3-D scorers
  (`validation/results/reference/README.md`, REF-de1e8f41…, REF-104a8234…, REF-ad4a98a0…).
- Assuming ~1–2·10³ steps per such history (proton steps plus delta electrons plus voxel
  boundaries), that is **≈ 1–3 µs per Geant4 step** (INFERRED, factor-3 uncertain). On-the-fly
  nucleus cross sections may raise it.
- The V3-010 carbon full-physics figure of ≈ 38 ms/history (plan §E) is not comparable: a full
  range, hadronic.

**Estimated wall time at 2 µs/step** (no scoring overhead beyond one exit-plane scorer):

| Slab | He 150 | C 290 | O 350 | p (×3 matched velocities) |
|---|---|---|---|---|
| 1 mm, 1e5 primaries | ≈ 2.0e3 s | ≈ 6.8e3 s | ≈ 8.5e3 s | ≈ 0.2e3 s each |
| 1 mm, 2e4 primaries | ≈ 400 s | ≈ 1.4e3 s | ≈ 1.7e3 s | < 100 s each |
| 10 mm, 1e5 primaries | ≈ 2.0e4 s | ≈ 6.8e4 s | ≈ 8.5e4 s | – |
| 50 mm, 1e5 primaries | ≈ 1e5 s | ≈ 3.4e5 s | ≈ 4.3e5 s | – |

**Consequences (INFERRED).**
- The r5 draft V5-SS (10 mm and 50 mm, ≥ 1e5 primaries, three ions) costs ≈ 10²–10³ h. It is
  **far above the 4 h REF budget** and cannot be registered as drafted.
- A **1 mm slab with ≈ 2e4 primaries per case** costs ≈ 1 h of REF for He/C/O plus three proton
  cases, in parts of ≤ 1200 s (C and O: 2 parts each). It fits in 4 h even if the per-step cost is
  3× higher, with fewer primaries if needed.
- Statistics: the window-trimmed rms width from N primaries has a relative standard error of order
  1/√(2N), inflated by trimming: ≈ 0.5–0.7 % at N = 2e4. That is adequate for a 5 % tolerance with
  |z| < 3.
- A thin slab keeps the energy loss negligible: ΔE/E ≈ 0.1 % for C 290 MeV/u at 1 mm (INFERRED).
  ionmc's per-step model is then compared at one well-defined (pv, p1v1 ≈ pv) state with b ≈ 9–11.
  The deeper-slab b-range (b ≈ 13–14.5) is out of reach within budget. That is a declared
  limitation, not a reason to relax the tolerance.

## 4. Scorer-free smoke (the measurement P10 still needs; execution only)

- **Bundle:** OpenTOPAS 4.3, `sv:Ph/Default/Modules = 1 "g4em-standard_SS"`, the matched electron
  cut, `d:Ph/Default/EMRangeMax` set explicitly ≥ 1.1 GeV (see P2 §F6), G4_WATER slab 1 mm thick,
  wide lateral extent, vacuum world, pencil beam, `GenericIon(6,12)` 3480 MeV.
- **No scorers at all** (no dose, no phase space, no fluence). Only:
  - the run log;
  - `/process/list` and `/process/em/printParameters` (or EM verbose 1);
  - wall time.
- Two history counts (e.g. 200 and 1000) separate initialisation from per-history cost.
- One run each for alpha (`alpha` naming per P4) and the proton is optional at the same size.
- **Archived:** the process list per particle (confirming `CoulombScat` for GenericIon and alpha and
  no msc), the model names, the maximum kinetic energy, and the wall times. No angular information
  exists in the output, so nothing is seen before pre-registration.
- **Seeds:** smokes use no evidence seed. If the runner requires one, use a rehearsal-family seed
  (2060xxxx) declared as such. TOPAS evidence seeds come from the case table only after
  registration.

## 5. Pre-registration content (to be fixed before any scored run; proposals only)

1. **Cases:**
   - He 150, C 290, O 350 MeV/u;
   - protons at the matched velocities (≈ 150, 290 and 350 MeV per nucleon-equivalent kinetic
     energy, i.e. T_p = (T/A)·m_p/u);
   - 1 mm G4_WATER slab. A 2 mm or 5 mm second thickness only if the smoke shows ≤ 1 µs/step.
2. **Scorer:** a phase space of primaries (generation filter, the projectile's Z) at the slab exit
   surface, direction cosines only. The projected angles θx = atan(u/w) and θy are pooled.
3. **Statistic:**
   - the projected rms in the ±2.5 θ window, iterated as in the plan (window anchored on the ionmc
     prediction, fixed before the run);
   - a bootstrap σ over primaries and parts.
4. **ionmc side:** the same slab and energy, warp-cpu f64, ≥ 1e5 histories, the same estimator.
5. **Tolerance (proposed):** |θ_ionmc/θ_SS − 1| ≤ 5 % (the declared ±5 % model uncertainty) and
   |z| < 3. For the relative ratio of section 2: the same 5 %, or a tighter value if pre-registered
   with a reason.
6. **Seeds and parts:** histories per part ≤ 1200 s at the smoke-measured rate. TOPAS seeds are
   formula-generated from case numbers not yet assigned in the r5 table (cases 1–40 are taken).
   Taking new case numbers requires a seed-ledger note.
7. **Lineage statement:** per section 2.
8. **Withdrawal rule:** if the smoke shows that the alias is missing, that GenericIon gets no
   single scattering, or a cost > 4 h for the minimum design (He/C/O at 1 mm, 1e4 primaries each),
   V5-SS is withdrawn. The lateral-gap intervention then applies after P6, as the plan states.

## 6. Status summary

| Item | Status |
|---|---|
| OpenTOPAS alias `g4em-standard_SS` → `G4EmStandardPhysicsSS` | LIT (docs, version not stated); smoke confirms |
| Alpha → `G4CoulombScattering(false)`, no msc, in Geant4 11.4.2 SS builder | FACT [EMB] L158–160 (conditional on the SS constructor calling `ConstructChargedSS`: OPEN) |
| GenericIon single scattering in 11.4.2 | LIT (10.4 listing) / OPEN: acquire `G4EmStandardPhysicsSS.cc` |
| Coulomb model physics (screening z-dependence, electrons, form factor) | OPEN: `G4CoulombScattering.cc`, `G4eCoulombScatteringModel.cc`, `G4WentzelOKandVIxSection.cc` |
| Cost | INFERRED: r5 design (10/50 mm, 1e5) infeasible; 1 mm × 2e4 ≈ 1 h REF; smoke measurement OPEN |
| V5-SS registration | not yet possible; requires the smoke and the SS constructor reading, then a reviewed plan revision |

## 7. Addendum 2026-11-01 (r6, Codex REVIEW-b469dff8 finding 3): constructor and Coulomb-model lineage

Planner, at worktree head `b3eef90b`. Source reading only (cat/grep/sed of the pinned files); no run. Sections 1–6 above
are the record as written at `3d2547ad` and are not changed. The "OPEN (11.4.2)" items of section 1 and the first two
lines of the section 6 table are superseded as follows.

### 7.1 Pinned sources (Geant4 v11.4.2; Geant4 Software License; role exploratory / lineage)

| Tag | File | sha256 | Bytes |
|---|---|---|---|
| [SS] | `physics_lists/constructors/electromagnetic/src/G4EmStandardPhysicsSS.cc` | `8aeaf91ce73c520e09e482d64a2e6a3b6c090124a9c14451d07b0b6eed9d9767` | 8 018 |
| [CS] | `processes/electromagnetic/standard/src/G4CoulombScattering.cc` | `db354d80c0a6efd278abf7793144bac1e5f3939a8d842ccd123014ea9ca0b754` | 7 415 |
| [eCS] | `processes/electromagnetic/standard/src/G4eCoulombScatteringModel.cc` | `fa5ddeb923b221c6c086593d0a006fd1f21502b5d460692117ba7e52e2f8c9af` | 10 716 |
| [WOKVI] | `processes/electromagnetic/standard/src/G4WentzelOKandVIxSection.cc` | `a64c94b1adf953beba1289575134514ca1261587e06b15b3286461418339a24c` | 14 478 |
| [WOKVI-h] | `processes/electromagnetic/standard/include/G4WentzelOKandVIxSection.hh` | `3a3f1db069d55ae06bd1c2e84c84c15a64dac2362371ed4a7dc783362184090d` | 7 445 |

[EMB] is `G4EmBuilder.cc` (40782f5c…, P2 record). The files are cached as `.ionmc-cache/reference/<sha256>`.

### 7.2 Constructor (FACT)

- [SS] L104–110: `param->SetDefaults()`, `SetMscThetaLimit(0.0)`, `SetUseMottCorrection(true)` (for e±), Auger and PIXE on.
- [SS] L213–219: GenericIon gets `G4ionIonisation` with `G4LindhardSorensenIonModel`, and `new G4CoulombScattering(false)`.
  No msc process is registered and no scattering model is set.
- [SS] L222 → [EMB] `ConstructChargedSS` (L301–343):
  - protons: `G4hIonisation` + `G4CoulombScattering(false)`, via `ConstructLightHadronsSS`, [EMB] L205–229;
  - alpha: `G4ionIonisation` + `G4CoulombScattering(false)`, via `ConstructIonEmPhysicsSS`, [EMB] L146–165.
  - In neither case is a model set.
- This settles the first two items of section 6: GenericIon and alpha single scattering are FACT. The OpenTOPAS alias stays LIT until the smoke.

### 7.3 Which single-scattering model each particle gets (FACT, [CS])

- [CS] L58–81: `G4CoulombScattering(false)`, process name `CoulombScat`, `isCombined = false`.
- [CS] L103–118: q2Max = 0.5 (FactorForAngleLimit · ħc/fm)². With `isCombined = false` the cross section is of type `fEmDecreasing`, splined. `MscThetaLimit = 0` gives no lower angular limit ([CS] L160–176: MinPrimaryEnergy = 0).
- **[CS] L120–148: the default model depends on the particle.**
  - `if (mass > CLHEP::GeV || p->GetParticleType() == "nucleus")`: no lambda table is built, `ion = true`, and the model is `G4IonCoulombScatteringModel`.
  - Otherwise the model is `G4eCoulombScatteringModel(isCombined)`.
  - Consequence: the **proton** (938 MeV, type "baryon") gets `G4eCoulombScatteringModel`; **alpha** (3.73 GeV, type "nucleus") and **GenericIon** (type "nucleus") get `G4IonCoulombScatteringModel`.
  - **The ions and the z = 1 reference of V5-SS use different single-scattering models.** The ratio design of the plan assumed one model on both sides.

**Proton path: `G4eCoulombScatteringModel` + `G4WentzelOKandVIxSection` (FACT).**

| Aspect | Source | Content |
|---|---|---|
| Cross section | [eCS] ComputeCrossSectionPerAtom; [WOKVI] SetupTarget | Wentzel (screened Rutherford). kinFactor = 2π (m_e r_e)² Z q² /(β² p²): the projectile charge q² and the target Z. The nuclear part is weighted by Z (×Z in the transport form), the electron part by 1 |
| Screening | [WOKVI] InitialiseA, SetupTarget | Thomas–Fermi radius a0 = m_e/0.88534 (× α), ScreenRSquare[Z] = ½ ScreeningFactor α² a0² (1 + exp(−0.001 Z²)) Z^(2/3). For Z > 1 and mass > 1 MeV: screenZ = min(1.13 Z, 1.13 + 3.76 Z² α² q²/β²) · ScreenRSquare[Z]/p². This is Molière's Born form with α′ = zZα/β, the form of the mapping's χ_a. For Z = 1 there is no Born term (screenZ = ScreenRSquare[1]/p²). The factor 1 + exp(−0.001 Z²) (≈ 1.94 for O) is a Geant4-specific constant |
| Electrons | [WOKVI] ComputeMaxElectronScattering; [eCS] SampleSecondaries | Scattering off atomic electrons is included, with the same screenZ, up to 1 − cos θ = min(cut, T_max) m_e/p². Above the cut, deflection by delta production belongs to the ionisation model (unread). With the matched large cut, all electron scattering stays in the Coulomb model, up to the kinematic limit. A sampled event is either electron or nuclear by the cross-section ratio |
| Form factor | [WOKVI] Initialise, SampleSingleScattering | `fNucFormfactor = G4EmParameters::NuclearFormfactorType()`: exponential 1/(1 + A_ff z1)², Gaussian or flat by rejection; FormFactor[Z] = 6.937e-6 MeV⁻² · A^0.54 (Butkevich et al., NIM A 488 (2002) 282). The default type is OPEN (`G4EmParameters.cc`, not pinned) |
| Spin / relativistic | [WOKVI] SampleSingleScattering | For spin ½ (the proton): rejection (1 − z1 β²/2 + (π/2) α Z β √(z1 β²/2)(2 − z1))·fm/(1 + z1 factD). At these angles (z1 ≈ θ²/4 ≲ 1e-4) the term is < 0.1 %. Spin-0 projectiles have factB = 0 |
| Recoil | [eCS] SampleSecondaries | First-order recoil kinetics; the recoil is deposited locally below the cut. Angles are in the lab via the momentum transfer |
| Limits | [CS] L150–154 | emin/emax from EmParameters and the model's own limits; no projectile-charge limit is visible in the read files |

### 7.4 Ion path: OPEN, STOP on this sub-item (files needed)

These were not in the cache and were not guessed. They decide whether branch A or B of the V5-SS row applies: the
screening radius (target-only Thomas–Fermi, or a projectile-dependent radius such as Lindhard or ZBL); whether the Born
factor uses z; whether target electrons are scattered at all; the form factor; the charge state.

- `source/processes/electromagnetic/standard/src/G4IonCoulombScatteringModel.cc`
- `source/processes/electromagnetic/standard/include/G4IonCoulombScatteringModel.hh`
- `source/processes/electromagnetic/standard/src/G4ionCoulombCrossSection.cc`
- `source/processes/electromagnetic/standard/include/G4ionCoulombCrossSection.hh`
- `source/processes/electromagnetic/utils/src/G4EmParameters.cc`: the defaults of `NuclearFormfactorType`, `ScreeningFactor` and `FactorForAngleLimit`, for both paths.
- `source/processes/electromagnetic/standard/include/G4eCoulombScatteringModel.hh`: defaults of `fixedCut`, `recoilThreshold`, `cosThetaMax` (proton path; low priority).

URL pattern: `https://gitlab.cern.ch/geant4/geant4/-/raw/v11.4.2/<path above>`.

### 7.5 What "independent" means for V5-SS (classification, r6)

- **Independent of** the multiple-scattering approximation used by ionmc:
  - Gottschalk's f_dM fit;
  - the Highland/Rossi prefactor and X_S;
  - the Molière angular expansion and Hanson width;
  - the b-matching root;
  - the per-step Gaussian sampling.
  The angular distribution is built from individual elastic events.
- **Not independent of** the screened-Rutherford cross-section physics. The proton path uses Molière's Born screening form (1.13 + 3.76 α′²) with a Thomas–Fermi-type radius, which is exactly the ingredient through which the exact-Ω mapping gets its charge dependence (beyond z²). The ion path is OPEN.
- **Classification:** held-out single-scattering MC evidence for how the mapping realises the shared screened-Coulomb model, including the angular integration and the thin-slab b. It is **not** independent physical evidence of the charge dependence of screening. That would need measurement, and P6 found none (EVIDENCE GAP).

### 7.6 The 3 % tolerance (r6)

At equal velocity, for the ion/proton ratio:
- **Cancels**, if both sides use the same model (branch A): the Thomas–Fermi constant and the Geant4 factor 1 + exp(−0.001 Z²); the electron share; the position of the nuclear-size cutoff and of the electron kinematic limit relative to χ_a. Both scale as 1/p, so in units of the screening angle they are velocity-only quantities.
- **Does not cancel:**
  - the Born-factor (z-dependent) shift of the form-factor and electron-limit logarithms (≲ 0.3 % and ≲ 0.5 % in θ, estimates);
  - the proton's spin-½ term (< 0.1 %);
  - the hydrogen branch without a Born term (< 0.1 %).
  - Sum ≲ 1 %: the pre-registered 3 % stands under branch A.
- **Branch B:** with different models on the two sides, nothing cancels by construction. For example, if the ion model omits electron scattering, θ_ion,SS falls by ≈ 6 % in water (the electron share of χ_c² is ≈ 13 %), and a 3 % test of the ratio would fail for a reason unrelated to ionmc. Under branch B, the proton cannot serve as the z = 1 reference. The freezing revision either uses ion/ion ratios at a common velocity, all on the ion model, or withdraws the row.

### 7.7 Status (r6)

| Item | Status |
|---|---|
| GenericIon / alpha / proton single scattering in 11.4.2 SS | FACT ([SS], [EMB]) |
| Model per particle | FACT ([CS] L120–148): proton `G4eCoulombScatteringModel`; alpha and GenericIon `G4IonCoulombScatteringModel` |
| Proton-path physics | FACT (table above); form-factor default type OPEN (`G4EmParameters.cc`) |
| Ion-path physics | **OPEN: STOP on this sub-item** (files in 7.4) |
| V5-SS registration | blocked until 7.4 is read and the branch is decided in the freezing plan revision |
