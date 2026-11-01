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

#### 7.4.1 Addendum 2026-11-01 (continuation of round V3-008-3 finding 3): the ion path read; the STOP on this sub-item is lifted

Planner, at worktree head `41a1f8ac` (DRAFT r6). Source reading only (cat/grep/sed of the cached files; role
exploratory / lineage; Geant4 Software License). Section 7.4 above is the record of the STOP and is not changed.
**File-name correction:** the class and file are `G4IonCoulombCrossSection` (capital I), not
`G4ionCoulombCrossSection` as written in 7.4 and in the plan.

| Tag | File (Geant4 v11.4.2) | sha256 | Bytes |
|---|---|---|---|
| [ICS] | `processes/electromagnetic/standard/src/G4IonCoulombScatteringModel.cc` | `ef9b31d40e99d72545dd06d329e13d7453bc20dea0a7681ce6ecc929e28d9296` | 8 471 |
| [ICS-h] | `processes/electromagnetic/standard/include/G4IonCoulombScatteringModel.hh` | `a1a9a4f5eb1dd0f323e5eac08873f435e87b739ff705c5b441caf4bb2262e4c0` | 5 956 |
| [IXS] | `processes/electromagnetic/standard/src/G4IonCoulombCrossSection.cc` | `87989b93748a8026c7b7031098166fa24c2f6d9c00979671bd0d9649e321e2f4` | 7 432 |
| [IXS-h] | `processes/electromagnetic/standard/include/G4IonCoulombCrossSection.hh` | `cd5d0b5cbb4dbc2963ededc5519d110ac1dd9042d38a9a8a0742051c46928b84` | 5 536 |
| [EMP] | `processes/electromagnetic/utils/src/G4EmParameters.cc` | `aaa143c1cb34b544f0d67fa6a7242d2b3fbd86ebfeb1bf4bc5c9fcca2973909e` | 43 246 |
| [eCS-h] | `processes/electromagnetic/standard/include/G4eCoulombScatteringModel.hh` | `68ce763e099bc31a8dfc67609ad39dcb7c070719a3a446f3971ed044255f109e` | 6 844 |

**Ion path, `G4IonCoulombScatteringModel` + `G4IonCoulombCrossSection` (FACT unless marked).**
- Origin: Consolandi 2010, "from G4eCoulombScatteringModel & G4CoulombScatteringModel"; reference Boschini et al.,
  arXiv:1011.4822 ([ICS] L33–48, [ICS-h] L43–52). The paper is LIT, not read.
- **Nuclear scattering only. Target electrons are not scattered.** [ICS] L125–148: the cross section per atom is
  `ioncross->NuclearCrossSection()`; the cut and emax arguments are unnamed (L129). [IXS] has no electron term
  (L186–219). The production cut enters only as the threshold for creating a recoil nucleus ([ICS] L220–236; below it
  the recoil energy is deposited locally as non-ionising).
- **Rutherford factor:** [IXS] L80–81, L190: 2π (m_e r_e)² Z² z²/(β² p_CM²), weight Z² per atom. z² is the squared
  PDG charge of the ion definition, i.e. the bare charge, no effective charge ([IXS-h] L146–156).
- **Kinematics:** [IXS] L107–135, relativistic reduced mass (Martynenko–Faustov): p_CM = p·M_t/E_cm,
  invbeta2 = 1 + μ²/p_CM². Because μ/p_CM = m/p, β_rel equals the lab β identically. The angle is sampled in the CM and
  boosted to the lab with exact two-body kinematics ([ICS] L184–211); the recoil energy is removed from the primary
  (L215–248). Target mass: m_p for Z = 1, the NIST atomic mass for Z > 1 in the cross section (L137–143), the isotope's
  nuclear mass in the sampling (L166–171). In momentum transfer this is the lab screened Rutherford law. At the V5-SS
  angles the CM treatment changes nothing: the lab deflection limit on hydrogen (≈ m_p/M ≥ 0.06 rad) is far outside
  the window (≤ 3 mrad).
- **Screening radius:** [IXS] L68, L167–182: x = a0 (z^0.23 + Z^0.23) with a0 = m_e c²/0.88534, ScreenRSquare = α² x².
  This is the universal (Ziegler–Biersack–Littmark-type) screening length 0.88534 a_B/(z^0.23 + Z^0.23). It
  **depends on the projectile charge**. The Thomas–Fermi Z^(1/3) is used only for the `G4Proton` definition (L177–179).
  There is no ScreeningFactor and no 1 + exp(−0.001 Z²) factor.
- **Born factor:** [IXS] L153–155: screenZ = ½ ScreenRSquare/p_CM² · (1.13 + corr · 3.76 Z² z² α²/β²), with
  corr = (10π Z z α)^0.04 (0.97–1.11 for He, C, O on H and O). There is no min(1.13 Z, …) cap, and the Born term also
  applies to hydrogen. The `heavycorr` argument is unused (L139–140).
- **Angular limits:** [IXS] L70–72, L132–133: cosθ_min = 1 and cosθ_max = −1, fixed; [ICS] L77 cosThetaMin = 1.0 has
  no setter ([ICS-h] L71–136). The PolarAngleLimit set by [CS] L152 and FactorForAngleLimit (q2Max) are not read. The
  law is integrated from 0 to θ_CM = π.
- **No nuclear form factor** (NuclearFormfactorType is not read) and **no spin term** (spin is stored, [IXS-h]
  L151–152, but unused in [IXS]).
- **Tables and limits:** [CS] L124–127: no lambda table for ions; the cross section is recomputed on the fly (cost:
  smoke). The energy limits come from [CS] L150–154 and [EMP] L146–147 (0.1 keV to 100 TeV); the G4VEmModel defaults
  are not read. There is no charge limit.

**Defaults ([EMP], `G4EmParameters::Initialise`, L118 ff.):** L157 FactorForAngleLimit = 1.0; L158 MscThetaLimit = π
(the SS constructor sets 0, [SS] L107); L166 ScreeningFactor = 1.0; L177 NuclearFormfactorType = fExponentialNF.
This closes the OPEN form-factor default of the proton path (7.3): it is exponential. [eCS-h] L147–151:
cosThetaMin 1, cosThetaMax −1, recoilThreshold 0, fixedCut −1, so the couple's electron production cut is used
([eCS] L186, L228); the L82 default `combined = true` is overridden by [CS] (`isCombined = false`).

**Term-by-term comparison** (mapping = `validation/scripts/ion_em/moliere_hanson.py`, decision 0043 §3):

| Term | Proton path ([eCS] + [WOKVI]) | Ion path ([ICS] + [IXS]) | ionmc mapping (Molière–Fano) |
|---|---|---|---|
| Screening radius | Thomas–Fermi a0 Z^(1/3), a0 = m_e/0.88534; × (1 + exp(−0.001 Z²)) for Z ≥ 2 (≈ 1.94 for O); ScreeningFactor 1 | universal a0 (z^0.23 + Z^0.23): projectile-dependent; no extra factor | Thomas–Fermi α m_e Z^(1/3)/0.885 |
| Born factor | min(1.13 Z, 1.13 + 3.76 α′²) for Z ≥ 2; none for H | 1.13 + corr · 3.76 α′², corr = (10π Z z α)^0.04; also for H; no cap | 1.13 + 3.76 α′² for every element |
| Electrons | weight Z per atom, same screenZ, up to 1 − cos θ = min(cut, T_max) m_e/p² | **none** | Z² in χ_c plus Fano Δb (z-independent at fixed β) |
| Form factor | exponential (default, FACT) | none | none (nuclear size is in the ±5 % statement) |
| Spin | spin-½ factor | none | none |
| Kinematics | lab, first-order recoil | CM, exact two-body | lab |

**Result:** the ion path differs from the proton path in every item that branch A required to be equal: screening
radius and constant, electron treatment and form factor. It also differs from the mapping in the radius (z^0.23
dependence), in corr and in the electrons.

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

#### 7.5.1 Addendum 2026-11-01 (continuation): branch B decided; quantification, design, tolerance and class

**Decision: branch B** of run condition (iii). The rule pre-registered in r6 ("Branch B: any of these differs")
applies on the source facts of 7.4.1 alone: the proton cannot serve as the z = 1 reference. The numbers below show
it is also quantitatively necessary: with branch A's design, an ionmc that equalled Molière–Hanson exactly would
miss the 3 % tolerance.

**Method** (SEEN pre-freeze; planner scratch, not committed; aggregates only; no ionmc and no TOPAS output used):
- The exact small-angle transform of the screened-Rutherford single-scattering law:
  Ω(η) = Σ_i 2 χ_c,i² ∫ θ (1 − J0(ηθ))/(θ² + a_i²)² dθ, in closed form through K1, with cut-off electron terms by
  quadrature.
- The projected distribution by Fourier inversion, and its rms in ±2.5 θ_H. θ_H is the Molière–Hanson width with
  Molière Z² plus the Fano Δb, as in `moliere_hanson.py`; it stands in for ionmc's Gaussian.
- 1 mm water at constant energy.
- Scripts in the session scratchpad, `v3-008/p10b/`: `ss_budget.py` (e50e3450…), `ss_check.py` (ff8be412…),
  `ss_cost.py` (2e378188…), `ss_signal.py` (1dc8b5eb…).
- Checks:
  - 1 − x K1(x): the series and the integral agree at the switch, and the small-x asymptote holds to 1e-7 relative;
  - η-grid refinement changes the pair biases by < 0.001 pp;
  - a ±5 % shift of the window anchor changes them by ≤ 0.07 pp;
  - P(|x| < 50 θ_H) = 0.99997.

**(a) Absolute windowed width, ion path relative to the mapping's Molière–Fano physics** (ions at the r6 cases):

| Case | b | no electrons | corr^0.04 | universal radius | ion path total | Molière_w / Gaussian_w |
|---|---|---|---|---|---|---|
| He 150 | 10.31 | −3.7 % | −0.03 % | −2.3 % | −6.0 % | +3.1 % |
| C 290 | 11.62 | −3.6 % | −0.14 % | −2.8 % | −6.5 % | +2.8 % |
| O 350 | 11.86 | −3.6 % | −0.19 % | −3.0 % | −6.7 % | +2.8 % |
| He 250 / C 430 / O 430 | 10.05–11.84 | −3.7 to −4.0 % | −0.02 to −0.18 % | −2.4 to −3.0 % | −6.4 to −6.8 % | +2.8 to +3.2 % |

Proton path at the matched velocities, relative to Molière–Fano: −2.7 to −2.8 %. Without its electron term it would be
−6.5 to −7.5 %; the electron limit 2 m_e/M ≈ 1.09 mrad is comparable to θ_H(p) = 0.72–1.93 mrad.
**Branch-A design** (ion/proton at equal velocity): the expected ratio_model/ratio_ref − 1 for an exact Molière–Hanson
ionmc is +3.9 % (He 150), +4.8 % (C 290), +5.1 % (O 350), +4.3 % (He 250), +5.0 % (C 430) and +5.1 % (O 430).
That is above 3 % before any ionmc error. **Correction to 7.6:** the r6 estimate "electron omission alone ≈ 6 %"
overstated the electron term. Removing the Fano electrons gives −3.3 to −4.1 %; the total of −6.0 to −6.8 % arises
only together with the universal radius.

**(b) Branch-B design (ion/ion at a common velocity, all on the ion path).**
- At equal velocity z/A = ½ for 4He, 12C and 16O, so χ_c is common and the ratio measures only the b-dependence:
  z² in the collision number, and the Born term.
- The universal radius and corr enter only through the difference of their logarithms. The missing electrons and the
  window estimator enter only through the different B of the two ions.

| Pair (MeV/u) | Signal ratio_ref − 1 | Expected bias ratio_model/ratio_ref − 1 | radius / corr / electrons / estimator | Alternative "b without z" | Alternative "no Born term" |
|---|---|---|---|---|---|
| He/C 150 | −5.6 % | −0.49 % | −0.45 / −0.14 / +0.42 / −0.32 | +6.4 % | −2.8 % |
| He/O 150 | −6.1 % | −0.69 % | −0.62 / −0.19 / +0.47 / −0.35 | +7.3 % | −4.1 % |
| He/C 250 | −6.5 % | −0.38 % | −0.41 / −0.13 / +0.52 / −0.36 | +7.3 % | −2.3 % |
| He/O 250 | −7.2 % | −0.57 % | −0.57 / −0.18 / +0.59 / −0.41 | +8.4 % | −3.4 % |
| C/O 430 | −1.0 % | −0.17 % | −0.15 / −0.06 / +0.09 / −0.06 | +1.2 % | −1.0 % |

(C/O at 290 and 350: bias −0.19 and −0.18 %, signal ≈ −1 %; not in the design.) The alternative columns are
Molière–Hanson ratio changes for a model without the z-dependence of b, or without the Born term.

**Design (branch B).**
- **Cases.** TOPAS `g4em-standard_SS`, 1 mm G4_WATER, pencil beam, 2e4 primaries per case in parts ≤ 1200 s:
  - He 150 (case 22), C 150 (25), O 150 (26);
  - He 250 (28), C 250 (27), O 250 (31);
  - C 430 (29), O 430 (30).
  - Cases 23 (C 290), 24 (O 350), 32 and 33 are retired unused; nothing has run. The seed formula is unchanged.
  - The proton cases are dropped.
- **Gated ratios (5):** He/C and He/O at 150 and at 250 MeV/u, and C/O at 430 MeV/u.
- **Estimator (unchanged):** projected rms in ±2.5 θ, the window anchored on the ionmc prediction per particle, x and
  y pooled. ionmc: the same slab, warp-cpu f64, ≥ 1e5 histories. Bootstrap σ over primaries and parts.
- **Electron production cut:** an energy cut ≥ T_max of the fastest case (430 MeV/u: 2 m_e c² β²γ²/(1 + 2γ m_e/M + …)
  ≈ 1.16 MeV, kinematics), the same for every case and checked in the smoke's production-cut printout. Then no δ ray
  is produced, and the reference contains no electron scattering of any kind. The Coulomb part has none (FACT); the
  ionisation part has no secondaries. Whether G4BetheBlochModel deflects the primary at δ production stays unread, and
  this choice makes the question irrelevant.
- **Smoke archive:** for GenericIon and alpha, the process `CoulombScat` with model `IonCoulombScattering`
  ([ICS-h] L75–76) and no msc. Proton smokes are no longer needed for V5-SS.
- **Cost (INFERRED):**
  - ion-path elastic events per primary per mm: He 3.6e3–5.0e3, C 1.1e4–1.6e4, O 1.4e4–1.9e4. That is 0.43–0.59 of
    the Molière Ω of §3, because the universal radius screens more strongly;
  - 8 cases × 2e4 ≈ 4.0e3 s at 2 µs/event (≈ 1.1 h), ≈ 3.3 h at 6 µs/event;
  - the on-the-fly cross section ([CS] L124–127) may raise the per-event cost, which the smoke measures.
- **Drop order** if the smoke-measured full design exceeds 4 h:
  1. C/O 430;
  2. the 250 MeV/u triplet.
  - The minimum design is the 150 MeV/u triplet: He/C and He/O at 150. If it exceeds 4 h, the row is withdrawn
    (§5 item 8, with the minimum redefined).
- **Absolute θ_ionmc/θ_SS per ion: report-only.** The expected reference-model factor θ_Gauss,w(θ_H)/θ_SS,w is
  +3.1 to +4.4 % over the eight cases (scratch, above). It is recorded now so that the report is not interpreted
  post hoc. It comes before ionmc's own offset from Molière–Hanson (I11-abs: +1.85 to +5.32 % on the I11 grid).

**Tolerance (set now from the source analysis, before any V5-SS run or smoke).**

| Component | Size |
|---|---|
| (1) Reference-model bias. Computed: −0.17 to −0.69 %. Added: ≤ 0.2 pp for the estimate's approximations (window anchor, small angle, one effective Fano species, constant energy) | ≤ 0.9 % |
| (2) ionmc against Molière–Hanson in a common-velocity ion/ion ratio. This equals (1 + r₁)/(1 + r₂), with r the I11-rel deviation of each ion against the same equal-velocity proton. The I11-rel gate (≤ 2 %) and the SEEN sign (all r > 0; spread 0.69–1.82 %) bound it. Assumption, declared: the thin-slab state (not on the I11 grid, which has x/R ≥ 0.1) behaves like the grid | ≤ 2 % (SEEN spread ≤ 1.13 pp) |
| (3) Statistics: σ of the ratio of two windowed widths at 2e4 each (bootstrap); inconclusive, counted as FAIL, if > 1 % | 0.7–1.0 % |

- Quadrature, the rule used since r5: √(0.9² + 2² + 1²) ≈ 2.4 %.
- Conservative variant, with the bias added linearly to the quadrature of (2) and (3): 3.0–3.1 %.
- **|ratio_model/ratio_ref − 1| ≤ 3 % at each of the five gated ratios (all must pass), with the bootstrap σ of the
  ratio ≤ 1 %. The 3 % is kept.** It is justified by the quadrature budget; the conservative variant reaches it only
  at the inconclusive limit σ = 1 %.
- The expected bias is negative (−0.2 to −0.7 %) and is recorded as the expected value.
- **Discrimination at 3 %:**
  - a model without the z-dependence of b fails every He pair (6.4–8.4 % plus the bias);
  - a model without the Born term fails the He/O pairs (−4.0 to −4.8 % including the bias) and sits at the tolerance
    for He/C (−2.7 to −3.3 %);
  - C/O 430 discriminates neither (≈ 1 %). It is a consistency check of the z 6 → 8 scaling at the most extrapolated
    f_dM states, and is sensitive only to errors > 3 %.
- Risk: five correlated ratios (He is shared) raise the chance of a false FAIL. The planning estimate is ≈ 5 % per
  ratio if the realised model term is at its SEEN spread.

**Class (branch B).** Held-out single-scattering MC on a Geant4 model different from the proton path.
- **Independent of** the multiple-scattering approximation: f_dM, Highland, the Molière angular integration, the
  b-matching and the per-step Gaussian.
- **Shares with the mapping** the screened Rutherford law (Z² z²/(β² p²)) and Molière's Born form 1.13 + 3.76 α′²,
  modified by corr.
- **Differs from it** in the screening radius (universal, z-dependent) and in the electrons (none). These are
  pre-computed above as bias.
- **What it is evidence for:** the mapping realises the b(z)-dependence of the Molière collision number (z² at fixed
  velocity; the Born term at the He/O level) in the thin-slab state.
- **What it is not:** independent physical evidence of the charge dependence of screening. The reference's own
  screening z-dependence (the universal radius) is another model, not data.
- **Consequence for P11b** (decision 0043 §3, note on the f_dM extrapolation, item 5):
  - Under branch B, V5-SS covers the extrapolation only differentially. The common-velocity ratio cancels the shared
    equal-velocity-proton factor. What remains is the difference in f_dM/Molière–Hanson drift between the two ions'
    equivalent states, ≤ 1.5 pp over 300 MeV to 3 GeV (I11-ext), below the tolerance.
  - The "independent coverage by the end cases" is therefore withdrawn. The extrapolation stays declared,
    theory-bounded (I11-ext) and carried in the ±5 %. The absolute V5-SS values are report-only context.

#### 7.5.2 Addendum 2026-11-01 (round V3-008-5, Codex REVIEW-6a14f4ec finding 2): TOST rule, statistics and power; design decision before any run

**Finding.** The round V3-008-4 module rule ("PASS iff |dev| ≤ 3 % whenever σ ≤ 1 %") ignored σ: it passed a 2.99 %
deviation at σ = 0.9 % (upper confidence bound 4.47 %). That contradicts the plan's scalar-comparison rule (§ Statistics:
pass iff |x − r| + 1.645 σ_comb ≤ tol). The tolerance text of 7.5.1 ("≤ 3 % … with the bootstrap σ of the ratio ≤ 1 %")
is superseded for the rule; the 7.5.1 numbers are unchanged.

**Rule (V5-SS, per gated ratio).** PASS iff |ratio_model/ratio_ref − 1| + 1.645 σ ≤ 3 %; σ = bootstrap SE of the ratio
statistic over primaries and parts; σ > 1 % stays a separate precondition (INCONCLUSIVE, counted as FAIL), kept because it
labels a statistically weak run as such instead of as a model failure. All five gated ratios must pass. `ss_budget.gate`
implements it; the test now requires the 2.99 % / 0.9 % example to FAIL and the expected −0.69 % at σ 0.30 % to PASS.

**Statistics (module §(h); SCRIPT-GENERATED, SEEN planning values).** The windowed moments P, M2, M4 of the reference
density (Geant4 ion path) in ±2.5 θ_H follow from the same transform (x⁴ kernel); the x/y joint moments of one primary from
the 2-D transform exp(−Ω(|k|)); the delta method gives the relative SE of the windowed rms (module docstring). Checks: the
moments of a Gaussian characteristic function reproduce the truncated-Gaussian closed form (1e-9) and the joint moments
factorise (1e-6); the untruncated limits 1/√(2N) (one projection) and 1/√(4N) (pooled) are reproduced.

| Quantity (8 design cases) | Value |
|---|---|
| P(window), reference | 0.972–0.976 |
| Kurtosis in the window (Gaussian truncated at 2.5 σ: 2.62) | 2.76–2.78 |
| Correlation of the x and y window terms of one primary | 0.04–0.05 |
| Relative SE of the rms at 2e4 primaries: one projection / x and y pooled | 0.48 % / 0.34 % |
| ionmc Gaussian, relative SE at 1e5 histories, pooled | 0.14 % |
| σ of the ratio statistic: 2e4 / 1e5 (round V3-008-4 design), pooled / one projection | 0.53 % / 0.73 % |
| σ: 6e4 / 4e5 (pre-registered below), pooled / one projection | 0.30 % / 0.41 % |

The planner's estimate "√2/√(2N) ≈ 0.5 % at 2e4" is reproduced (0.53 %, incl. the ionmc side). The budget item (3)
"0.7–1.0 %" of 7.5.1 corresponds to one projection per primary; the design pools x and y.

**Power: P(all five ratios pass) under the TOST rule at tol 3 %** (true deviation per ratio = the computed bias of
7.5.1, −0.49 / −0.69 / −0.38 / −0.57 / −0.17 %, plus the stated shift; Monte Carlo over the eight reference and eight ionmc
cases, which correlates ratios sharing a case; 4e5 draws, seed 20261101; σ taken as known):

| Scenario (shift of every ratio) | 2e4 / 1e5 | 2e4 / 4e5 | 4e4 / 4e5 | **6e4 / 4e5** |
|---|---|---|---|---|
| S0 correct model (0) | 0.994 | 0.998 | 1.000 | **1.000** |
| S1 −0.2 pp (approximation allowance of the bias) | 0.982 | 0.991 | 1.000 | **1.000** |
| S2 −1.33 pp (S1 + SEEN I11-rel spread 1.13 pp, adverse) | 0.254 | 0.319 | 0.751 | **0.918** |
| S3 −2.2 pp (S1 + I11-rel bound 2 pp, adverse) | 0.0006 | 0.0007 | 0.003 | **0.006** |
| REF cost (scaled from the 23 ms/primary C 290 smoke) | 0.98 h | 0.98 h | 1.96 h | **2.94 h** |

(One projection instead of pooling, S0: 0.83 at 2e4 / 1e5, 0.9999 at 6e4 / 4e5.) False-fail probability for a correct
model: 0.6 % at 2e4 / 1e5, < 0.01 % at 6e4 / 4e5. Per-ratio values are printed by the module.

**Decision (pre-registered 2026-11-01, before any V5-SS run or scored smoke): raise N.**
- TOPAS 6e4 primaries per case: 378 s (He 250) to 2022 s (O 150) per case, ≤ 2 parts of ≤ 1200 s (parts 0–1 of the seed
  formula), total 10 596 s ≈ 2.94 h ≤ 4 h (upper estimate; same scaling as 7.8). ionmc ≥ 4e5 histories per case (cost not
  measured; a 1 mm slab, expected minutes).
- If the measured cost of this design exceeds 4 h, N falls back to 4e4 (1.96 h) before the case drop order of 7.5.1.
- **Tolerance 3 % kept; five gated ratios kept.** Under TOST the statistics enter through 1.645 σ (0.49 % at the design),
  so the tolerance covers only the systematic items: `ss_budget.budget_systematic` gives bias + model = 2.90 % linear
  (2.19 % quadrature) at the item bounds. The room left for ionmc's own model term is 3 − 0.9 − 0.49 ≈ 1.6 pp at the bias
  bound (≈ 1.8 pp at the computed bias), above the SEEN I11-rel spread (1.13 pp) but below the I11-rel bound (2 pp).
  V5-SS is therefore deliberately tighter than the I11-rel band: a FAIL caused by a model term near 2 pp (S3) is a
  finding, not a false fail.
- Not chosen: 3.5 % (a revision the budget does not need at this N; it would let a model without the Born term pass
  He/C 250, −2.66 % + 0.49 = 3.15 %); fewer gated ratios (C/O 430 has the smallest bias and does not limit the power).

### 7.6 The 3 % tolerance (r6)

At equal velocity, for the ion/proton ratio:
- **Cancels**, if both sides use the same model (branch A): the Thomas–Fermi constant and the Geant4 factor 1 + exp(−0.001 Z²); the electron share; the position of the nuclear-size cutoff and of the electron kinematic limit relative to χ_a. Both scale as 1/p, so in units of the screening angle they are velocity-only quantities.
- **Does not cancel:**
  - the Born-factor (z-dependent) shift of the form-factor and electron-limit logarithms (≲ 0.3 % and ≲ 0.5 % in θ, estimates);
  - the proton's spin-½ term (< 0.1 %);
  - the hydrogen branch without a Born term (< 0.1 %).
  - Sum ≲ 1 %: the pre-registered 3 % stands under branch A.
- **Branch B:** with different models on the two sides, nothing cancels by construction. For example, if the ion model omits electron scattering, θ_ion,SS falls by ≈ 6 % in water (the electron share of χ_c² is ≈ 13 %), and a 3 % test of the ratio would fail for a reason unrelated to ionmc. Under branch B, the proton cannot serve as the z = 1 reference. The freezing revision either uses ion/ion ratios at a common velocity, all on the ion model, or withdraws the row.

#### 7.6.1 Addendum 2026-11-01 (continuation): the tolerance under branch B

- The branch-A argument of 7.6 does not apply: branch B is decided (7.5.1).
- The 3 % is re-justified for the ion/ion design from the source analysis (7.5.1, budget table) and is kept. It is
  set now, before any V5-SS run or SS smoke, and recorded as such.
- The electron figure of 7.6 is corrected in 7.5.1: −3.3 to −4.1 %, not ≈ 6 %.

#### 7.6.2 Addendum 2026-11-01 (continuation): whether to file the lateral-gap intervention regardless of V5-SS (recommendation only; nothing filed)

**For filing at the freeze, whatever V5-SS shows:**
- Under branch B, V5-SS is shared-lineage evidence. It uses the same screened Rutherford law and the same Born form
  as the mapping. Its own screening z-dependence (the universal radius) is a different model, not data.
- P6 found no measured He/C/O lateral data (EVIDENCE GAP). So the charge dependence of lateral spreading has no
  independent physical support, whether V5-SS passes or fails.
- Branch B also removes V5-SS's absolute coverage of the f_dM extrapolation (7.5.1). That was part of the
  finding-4 disposition, so a second, related gap now exists regardless of the outcome.
- The protocol requires an authentic intervention for missing essential external evidence after a bounded
  investigation, and P6 has completed that investigation.
- Filing early avoids a serial block late in integration, and the operator's answer (data access or gap acceptance)
  does not depend on the V5-SS result.

**Against:**
- The r5 sequence uses the intervention only after P6 and V5-SS fail. Filing earlier changes a pre-registered order.
  That is harmless before freeze, but it needs a dated amendment.
- The verbatim wording offers "or an accessible single-scattering MC reference" as an acceptable input. V5-SS
  literally is one, so the request may read as already partly met unless the ledger states its shared-lineage class.
- If V5-SS passes, an early request may look like over-claiming the gap.
- It costs operator attention.

**Recommendation:**
- File at the freeze, regardless of V5-SS, with the verbatim wording.
- Add a ledger statement that:
  - V5-SS is shared-lineage, implementation-level evidence on the ion path;
  - the f_dM extrapolation above 300 MeV proton-equivalent has no independent coverage.
- Because the verbatim text cannot be edited here, the orchestrator also decides whether its clause "or an accessible
  single-scattering MC reference" needs a reviewed amendment, e.g. "an accessible single-scattering MC reference
  independent of the screened-Coulomb model used by the mapping".

### 7.7 Status (r6)

| Item | Status |
|---|---|
| GenericIon / alpha / proton single scattering in 11.4.2 SS | FACT ([SS], [EMB]) |
| Model per particle | FACT ([CS] L120–148): proton `G4eCoulombScatteringModel`; alpha and GenericIon `G4IonCoulombScatteringModel` |
| Proton-path physics | FACT (table above); form-factor default type OPEN (`G4EmParameters.cc`) |
| Ion-path physics | **OPEN: STOP on this sub-item** (files in 7.4) |
| V5-SS registration | blocked until 7.4 is read and the branch is decided in the freezing plan revision |

#### 7.7.1 Status addendum 2026-11-01 (continuation)

| Item | Status |
|---|---|
| Ion-path physics | FACT (7.4.1): nuclear only, no electrons, universal z-dependent screening radius, Born factor with corr^0.04, no form factor, no spin term, exact CM kinematics |
| Proton-path form-factor default | FACT ([EMP] L177): exponential. ScreeningFactor 1.0 (L166), FactorForAngleLimit 1.0 (L157) |
| Branch of run condition (iii) | **B** (7.5.1): the proton is not a valid z = 1 reference; ion/ion ratios at a common velocity on the ion path |
| Tolerance | 3 % kept, re-justified for branch B (quadrature 2.4 %, conservative 3.0–3.1 %); set before any run |
| f_dM extrapolation coverage by V5-SS | differential only; the absolute coverage is withdrawn (7.5.1) |
| V5-SS registration | unblocked on lineage. It still needs the SS smoke (alias, process and model list, production-cut printout, cost) and the reviewed freezing revision |
| Lateral-gap intervention | recommendation in 7.6.2; an orchestrator decision before freeze |

### 7.8 Addendum 2026-11-01 (round V3-008-4, Codex REVIEW-e0b586fd finding 4): reproducible calculation

Planner, at head `ce4d3d6a` (files uncommitted at the time of writing; bound at the commit that adds them). The
session-scratch scripts of 7.5.1 (`ss_budget.py` e50e3450…, `ss_check.py` ff8be412…, `ss_cost.py` 2e378188…,
`ss_signal.py` 1dc8b5eb…; hashes re-checked, unchanged) are replaced as the record by one committed module with tests.
Same constants, same numerics (η grid 200 001 nodes to 40/θ_H, window 2.5 θ_H, J0 and K1 quadratures). Changes: pure
functions with docstrings citing the equations and the Geant4 source lines ([CS] L120–148; [ICS] L125–148; [IXS] L68,
L80–81, L153–155, L167–182, L186–219; [WOKVI] functions as in 7.3); no computation or I/O at import (the G(S) table is
built lazily); J0 is evaluated in row chunks (identical arithmetic per row); a cost scaling from a measured time per
primary is added. Every number is SEEN pre-freeze, a planning value, not qualification evidence.

| File | sha256 | Bytes |
|---|---|---|
| `validation/scripts/ion_em/ss_budget.py` | `b893c0420bb492368d89943a8bd8820752ae903375921f9005c7dd2bc8f00156` | 32 261 |
| `tests/ionmc/test_ion_em_ss_budget.py` | `551c75838dbedbf2be0a7b643beb1cd4b5ce48d2d15e15f0f4883bf00be1cfcc` | 8 445 |

**Command** (numpy only, ≈ 10 s, ≈ 0.1 GB): `uv run --extra dev python -I validation/scripts/ion_em/ss_budget.py`
prints (a) the absolute windowed widths, (b) the superseded branch-A table, (c) the branch-B pair signals, biases and
their decomposition, (d) the discrimination signals, (e) the convergence and window-anchor checks, (f) the cost table,
(g) the budget. Tests: `pytest tests/ionmc/test_ion_em_ss_budget.py` (14 tests, ≈ 5 s).

**Reproduction (module output against the values recorded in 7.5.1 from the scratch scripts): every recorded aggregate
is reproduced unchanged at its recorded precision.**

| Quantity | Recorded (7.5.1) | Module |
|---|---|---|
| Gated pair bias | −0.17 to −0.69 % | −0.17 to −0.69 % (−0.1706 to −0.6936) |
| Decomposition per pair (radius / corr / electrons / estimator) | table of 7.5.1 | identical to two decimals |
| Electron share (no electrons, 13 cases) | −3.3 to −4.1 % | −3.26 to −4.14 % (the 8 design cases: −3.26 to −4.04 %) |
| Universal radius | −2.3 to −3.0 % | −2.32 to −3.00 % |
| He-pair signal ratio_ref − 1 | −5.6 to −7.2 % | −5.58 to −7.22 %; C/O 430 −1.02 % |
| "b without z" / "no Born" (He pairs) | +6.4 to +8.4 % / −2.3 to −4.1 % | +6.43 to +8.39 % / −2.28 to −4.11 % |
| Branch A, exact Molière–Hanson ionmc | +3.9 to +5.1 % | +3.90 to +5.14 % |
| Expected reference factor θ_Gauss,w/θ_SS,w − 1 (8 cases) | +3.1 to +4.4 % | +3.12 to +4.35 % |
| Grid refinement / window anchor ±5 % | < 0.001 pp / ≤ 0.07 pp | 0.000 pp / ≤ 0.066 pp |
| P(\|x\| < 50 θ_H) | 0.99997 | 0.99997 |
| Events per primary per mm; ratio to Ω₀ | He 3.6e3–5.0e3, C 1.1e4–1.6e4, O 1.4e4–1.9e4; 0.43–0.59 | the same |
| Cost, 8 × 2e4 at 2 µs/event | ≈ 4.0e3 s (1.1 h); 3.3 h at 6 µs | 3991 s (1.11 h); 3.33 h |
| Budget: quadrature / conservative | ≈ 2.4 % / 3.0–3.1 % | 2.41 % / 3.02–3.14 % |

**Tests (CI-fast).** 1 − x K1(x) at the series/integral switch against 1 − 2 K1(2) and its small-x asymptote; J0
(zero and J0(1)); the closed-form Ω against direct quadrature of its defining integral (η a = 0.3, 1, 4); Ω's small-η
limit equals Molière's (χ_c² η²/4)(−ln(a² η²/4) + 1 − 2 C_E) and its large-η limit the collision number χ_c²/a²; the
windowed transform of a pure Gaussian characteristic function reproduces the closed-form truncated Gaussian (rms and
probability, 1e-6); the bias terms vanish exactly for equal models and for identical ions, and the decomposition
telescopes to the total; the Born-only z-dependence of the Thomas–Fermi variant; equal-velocity kinematics; golden
aggregates for the five gated ratios (signal, bias and the four terms; SCRIPT-GENERATED by the module, tolerance
0.002 pp), the discrimination signals and the cost table; and a regression that the budget gives quadrature 2.4 % and
conservative 3.0–3.1 % and that the rule (PASS iff |dev| ≤ 3 % with σ ≤ 1 %, else FAIL or INCONCLUSIVE = FAIL)
evaluates as recorded. The golden values are a regression guard on SEEN planning values, not evidence.

**Cost from the 1 mm smoke (dated note 2026-11-01).** Facts as relayed by the orchestrator (the planner did not re-read
the run): `carbon-water-290mevu-ss-slab1mm-smoke` (GenericIon(6,12), 3480 MeV, 100 histories, 1 mm G4_WATER,
`g4em-standard_SS`, `Ts/MaxStepNumber` 1e8), REF-35df82fcd7723f6a9a0c-db361e98 at `ce4d3d6a`: no killed track; `IonCoulombScattering` for
alpha and GenericIon; wall time ≈ 23 ms per carbon primary. Computation: the module gives 1.299e4 ion-path events per
C 290 primary in the slab, so 23 ms corresponds to ≈ 1.77 µs per event (within the 2 µs planning value). Scaling by the
events per primary of each design case: 126 s (He 250) to 674 s (O 150) per case at 2e4 primaries, i.e. one part
≤ 1200 s per case, and **≈ 3.5e3 s ≈ 1.0 h for the 8 cases** (an upper estimate if part of the 23 ms is fixed
per-run overhead). Under the 4 h condition with a factor ≈ 4 margin; the drop order is not triggered. Minutes per case,
about an hour in total, not several hours.

**Dated note 2026-11-01 (round V3-008-5, Codex REVIEW-6a14f4ec finding 2).** The module and test of the table above are
superseded by the round V3-008-5 versions (sha256 bound at commit, plan P9 (iv)). Changes: `gate` implements the TOST rule
(7.5.2); new functions `budget_systematic`, `windowed_moments_char`, `windowed_joint_char`, `rms_rel_se`,
`gauss_window_moments`, `case_statistics`, `pair_sigma`, `tost_pass_probability`, `ref_cost_parts`; constants `Z_TOST`,
`N_REF_PRIMARIES` = 6e4, `N_MODEL_HISTORIES` = 4e5 (the round V3-008-4 values kept as `*_R6`); `main` prints §(h). Every
earlier function, constant and golden value is unchanged (the round V3-008-4 `budget` is kept for the record). Tests: 20
(the rule test replaced; window kernels, Gaussian moments and factorisation, SE limits, window statistics and TOST power/
cost goldens added; SCRIPT-GENERATED regression values, not evidence).
