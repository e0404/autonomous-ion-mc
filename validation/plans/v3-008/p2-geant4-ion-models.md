# P2: Geant4 11.4.2 ion models read from source (pre-freeze check, V3-008)

Status: source reading for pre-freeze check P2 of `validation/plans/v3-008-acceptance.md` (DRAFT r3).
Worktree head `2632a465`. Written 2026-11-01 by the research specialist. This note contains no
transport output, no TOPAS output and no tabulated reference data values. The Urban `ComputeTheta0`
coefficients below are source-code constants. The plan requires P2 to freeze them, so they are
quoted with line references. The tabulated arrays (`celectron`, `cpositron`, `sig0`, `hecorr`) are
cited by line only. A script that needs them reads them from the pinned cached file.

Markers:
- **FACT (quoted)**: read in a pinned file below, with the line number and verbatim code.
- **INFERRED**: derived from FACT lines by algebra, or recalled from Geant4 code that is *not* in
  the pinned set. Each INFERRED item names the file that would make it FACT.
- **OPEN**: not settled. The settling file or run-time evidence is named.

## 0. Provenance (sources read)

All files are Geant4 tag v11.4.2, Geant4 Software Licence (attribution required; quoting short
excerpts with attribution is permitted; nothing here is redistributed as data). Read with plain
tools only, never compiled or executed.

| Short name | File (upstream path under `source/`) | sha256 | Bytes | Lines | Acquisition id | Materialized |
|---|---|---|---|---|---|---|
| [URB] | `processes/electromagnetic/standard/src/G4UrbanMscModel.cc` | `22c5f52b3a1cfc8dd99dfd07929bed702cb2a35a0d3a7c727450bfbfc9b69558` | 41 382 | 1155 | `81d32d16…` | yes, `.ionmc-cache/reference/<sha256>` |
| [LS] | `processes/electromagnetic/standard/src/G4LindhardSorensenIonModel.cc` | `e16ed6cd1d132183684ef47948b6f8fe36526a48663a6dc1878d5bb8330eb12c` | 16 490 | 447 | `0a4b7456…` | yes |
| [ICRU73] | `processes/electromagnetic/standard/src/G4IonICRU73Data.cc` | `cf771a5a44ea936c56593b7dd2d67bb4157345eb316815f60880986404fe08d2` | 13 701 | – | `33881db243e64160be539305e8782dd7` | **no**: blocked until the worktree is clean, so not read here |
| [OPT4] | `physics_lists/constructors/electromagnetic/src/G4EmStandardPhysics_option4.cc` | `04337300cc119a4fd893d515059cc75cef9f042849b183848643a4c4bf104c89` | 11 966 | 334 | `4100c163…` | yes |

- The sha256 and byte counts of the three materialized files were re-computed here and match.
- The coordinator assigned acquisition ids to files in the order listed above. Only the [ICRU73] id
  was given in full.

## 1. Option 4 assignments for GenericIon and alpha

### Ionisation

- FACT [OPT4] L316–324. GenericIon gets `G4ionIonisation` with the LS model unless an external dE/dx
  is requested:
  ```
  particle = G4GenericIon::GenericIon();
  G4ionIonisation* ionIoni = new G4ionIonisation();
  if (fUseExternalDEDX) {
    ionIoni->SetEmModel(new G4IonParametrisedLossModel());
  } else {
    ionIoni->SetEmModel(new G4LindhardSorensenIonModel());
  }
  ph->RegisterProcess(hmsc, particle);
  ph->RegisterProcess(ionIoni, particle);
  ```
  - `fUseExternalDEDX` is set only when the constructor receives a non-empty name. L137:
    `if (!nam.empty()) { fUseExternalDEDX = true; }`.
  - INFERRED: the default argument is an empty name (needs `G4EmStandardPhysics_option4.hh`).
  - OPEN for TOPAS: how TOPAS instantiates `g4em-standard_opt4`. Settled by the P4 run-time EM
    printout: model name "LindhardSorensen" versus "ParamICRU73".
- The LS model sets no explicit energy limits in [OPT4]. INFERRED: it covers the whole ion energy
  range, and its internal switch is `fElimit` (section 3).
- FACT [OPT4] L131: `param->SetUseICRU90Data(true);`
- FACT [OPT4] L122–124: step functions (ionisation step limitation, which affects the step pattern):
  ```
  param->SetStepFunctionMuHad(0.1, 50*CLHEP::um);
  param->SetStepFunctionLightIons(0.1, 20*CLHEP::um);
  param->SetStepFunctionIons(0.1, 1*CLHEP::um);
  ```
  INFERRED: "LightIons" applies to alpha (and d, t, He3) and "Ions" to GenericIon (needs
  `processes/electromagnetic/utils/src/G4EmParameters.cc` and `.../utils/src/G4VEnergyLossProcess.cc`).
- **Alpha is not configured in [OPT4].** It is handled inside `G4EmBuilder::ConstructCharged(hmsc, pnuc)`
  (FACT L328).
  - INFERRED: alpha gets `new G4hMultipleScattering("ionmsc")` and a `G4ionIonisation` with default
    models (Bragg-type below ≈ 2 MeV·M/m_p, `G4BetheBlochModel` above).
  - OPEN: needs `physics_lists/constructors/electromagnetic/src/G4EmBuilder.cc` and
    `processes/electromagnetic/standard/src/G4ionIonisation.cc`.

### Multiple scattering

- FACT [OPT4] L161: `G4hMultipleScattering* hmsc = new G4hMultipleScattering("ionmsc");`
  - No `SetEmModel` is called on `hmsc` anywhere in [OPT4].
  - The same object is registered for GenericIon (L323) and passed to `ConstructCharged` (L328).
- INFERRED: the default model of `G4hMultipleScattering` is `G4UrbanMscModel` (in
  `InitialiseProcess`: `if(nullptr == EmModel(0)) { SetEmModel(new G4UrbanMscModel()); }`). The
  project's own `docs/research/em-physics.md` L79 also says "Urban MSC for GenericIon" [NV].
  - OPEN: needs `processes/electromagnetic/standard/src/G4hMultipleScattering.cc`. Also settled by
    the P4 printout (`msc: for GenericIon … UrbanMsc`).
- Msc parameters set in [OPT4] (FACT L125–129):
  ```
  param->SetUseMottCorrection(true); // use Mott-correction for e-/e+ msc gs
  param->SetMscStepLimitType(fUseSafetyPlus); // for e-/e+ msc gs
  param->SetMscSkin(3);              // error-free stepping for e-/e+ msc gs
  param->SetMscRangeFactor(0.08);    // error-free stepping for e-/e+ msc gs
  param->SetMuHadLateralDisplacement(true);
  ```
  - FACT: [OPT4] does not call `SetMscMuHadStepLimitType`, `SetMscMuHadRangeFactor`,
    `SetLateralDisplacementAlg96` or `SetMscLambdaLimit`, so their library defaults apply.
- **Step-limit type, RangeFactor and lateral displacement for ions and alpha: OPEN.** Two code paths
  compete; both are recalled, not read:
  - (a) `G4EmTableUtil::PrepareMscProcess`, for particle type "nucleus" (alpha and GenericIon): `fMinimal`,
    `latDisplacement = false`, `facrange = 0.2`.
  - (b) `G4VMscModel::InitialiseParameters(p)`, called at [URB] L139. For non-e± particles it takes
    `MscMuHadStepLimitType()` (default believed `fMinimal`), `MscMuHadRangeFactor()` (default
    believed 0.2) and `MuHadLateralDisplacement()`, which [OPT4] sets to `true`.
  - Which one wins decides the lateral displacement flag for ions.
  - Files: `processes/electromagnetic/utils/src/G4EmTableUtil.cc`, `.../utils/src/G4VMscModel.cc`,
    `.../utils/src/G4VMultipleScattering.cc`, `.../utils/src/G4EmParameters.cc`.
  - Run-time evidence: the Geant4 EM table printout per msc process ("RangeFactor= …, stepLimType: …,
    latDisp: …") for GenericIon, alpha and proton. **Recommend that P4 archives these lines.**
- Skin: `SetMscSkin(3)` is used by [URB] only in the `fUseDistanceToBoundary` branch (FACT L494–589,
  `skindepth = skin*stepmin` L515). INFERRED irrelevant for ions if they run `fMinimal`.
- FACT [URB] L661–672: if the ion step-limit type is `fMinimal` (the `else` branch), msc limits the
  step only after a geometry boundary, to `facrange*max(range, lambda0)`:
  ```
  if (stepStatus == fGeomBoundary)
    {
      tlimit = (currentRange > lambda0)
        ? facrange*currentRange : facrange*lambda0;
  ```
  - At track start `tlimit = … = geombig` (FACT L424).
  - INFERRED: λ0 ≫ R for ions, so msc never limits ion steps in practice.
  - **Steps are then set by** the ionisation step function, delta production (cut-dependent),
    geometry and scoring boundaries, and any user max step. Section 6 shows why this matters.
- FACT [OPT4] L331: `G4EmModelActivator mact(GetPhysicsName());`
  - INFERRED: it applies only per-region model overrides requested through G4EmParameters (none in
    the template bundle).
- **Protons (needed for P2b): not in the pinned set.**
  - INFERRED (strong): `G4EmBuilder::ConstructCharged(hmsc, pnuc)` with its default `isWVI = true`
    gives protons `G4WentzelVIModel` + `G4CoulombScattering`, not Urban.
  - The project's `docs/research/validation-evidence.md` L53 states "opt4: WentzelVI for protons [K]",
    and `docs/research/em-physics.md` L79 says the same.
  - OPEN until `G4EmBuilder.cc` is read or the archived T15 TOPAS log is inspected.
  - Consequence for P2b: section 6.4.

## 2. G4UrbanMscModel: θ0 and the angular distribution

### 2.1 θ0 (central width), FACT [URB] L982–1028

```
G4double invbetacp = (kinEnergy+mass)/(kinEnergy*(kinEnergy+2.*mass));
if(currentKinEnergy != kinEnergy) {
  invbetacp = std::sqrt(invbetacp*(currentKinEnergy+mass)/
                        (currentKinEnergy*(currentKinEnergy+2.*mass)));
}
G4double y = trueStepLength/currentRadLength;
...                                   // positron-only branch L995–1021
static const G4double c_highland = 13.6*CLHEP::MeV;
G4double theta0 = c_highland*std::abs(charge)*std::sqrt(y)*invbetacp;
// correction factor from e- scattering data
theta0 *= (msc[idx]->coeffth1+msc[idx]->coeffth2*G4Log(y));
```

- **y = t/X0** is the true path length over the material radiation length. FACT L993, with
  `currentRadLength = couple->GetMaterial()->GetRadlen()` at L850.
- The log is the **natural log of y**. There is no lg(pv), no β and **no charge in the logarithm**
  (FACT L1027).
- `invbetacp = (E+M)/(E(E+2M)) = 1/(pβc)` (exact algebra, INFERRED).
  - With a post-step energy E2 ≠ E1 it is the geometric mean of 1/(pβc) at E1 (pre) and E2 (post)
    (FACT L988–992).
  - E2 is set in `SampleScattering` (FACT L794–800):
    - E2 = E(R − t) from the range table if t > R·dtrl;
    - E2 = E1 − t·dE/dx(E1) if t > 0.01 R;
    - otherwise E2 = E1.
  - `dtrl` is a G4VMscModel member; INFERRED default 0.05 (needs `G4VMscModel.cc`/`.hh`).
- **Charge dependence of θ0: linear, `std::abs(charge)` in the prefactor only** (FACT L1024).
  - INFERRED: `charge` is the PDG charge of the track's particle definition, i.e. the bare z,
    not an effective charge. It is set by `SetParticle(track->GetDynamicParticle()->GetDefinition())`
    in `StartTracking` (FACT L420).
  - OPEN: the inline `SetParticle` body (`mass = p->GetPDGMass(); charge = p->GetPDGCharge()/eplus;
    chargeSquare = charge*charge`) is in `processes/electromagnetic/standard/include/G4UrbanMscModel.hh`.
- **Mass dependence:** only through `invbetacp` (exact kinematics) and the electron-equivalent energy
  of the cross section (section 2.3).
  - FACT: no branch in [URB] tests for ions; a case-insensitive search for "ion/GenericIon/nucleus/
    alpha" finds only comments.
  - The only mass branches are `mass < masslimite` (0.6 MeV, e± only; L83, L473, L599, L630),
    `mass > electron_mass_c2` (L296) and `currentKinEnergy < mass` (L712).
- Z-dependent coefficients, FACT L1115–1125 (per material-cuts couple):
  ```
  G4double Zeff = aCouple->GetMaterial()->GetIonisation()->GetZeffective();
  ...
  G4double lnZ = G4Log(Zeff);
  // correction in theta0 formula
  G4double w = G4Exp(lnZ/6.);
  G4double facz = 0.990395+w*(-0.168386+w*0.093286);
  msc[j]->coeffth1 = facz*(1. - 8.7780e-2/Zeff);
  msc[j]->coeffth2 = facz*(4.0780e-2 + 1.7315e-4*Zeff);
  ```
- Small-step rule, FACT L872–881:
  - θ0 = √(t/tsmall)·θ0(tsmall) if t ≤ tsmall = min(tlimitmin, lambdalimit).
  - With `fMinimal`, `tlimitmin` stays 10·tlimitminfix = 0.1 nm (FACT L89, L96, L427); it is
    recomputed only in the other three step-limit branches (L516, L608, L638).
  - INFERRED: irrelevant for ion steps.

### 2.2 Tail and sampling, FACT [URB] L834–978 and L1127–1132

- Tail coefficients:
  ```
  G4double Z13 = w*w;
  msc[j]->coeffc1 = 2.3785    - Z13*(4.1981e-1 - Z13*6.3100e-2);
  msc[j]->coeffc2 = 4.7526e-1 + Z13*(1.7694    - Z13*3.3885e-1);
  msc[j]->coeffc3 = 2.3683e-1 - Z13*(1.8111    - Z13*3.2774e-1);
  msc[j]->coeffc4 = 1.7888e-2 + Z13*(1.9659e-2 - Z13*2.6664e-3);
  ```
- Tail parameter (L900–909):
  - `u = exp(ln(tau)/6)` with τ = t/λ (λ averaged over the step, L841–846);
  - `xx = ln(lambdaeff/X0)`;
  - `xsi = coeffc1 + u*(coeffc2 + coeffc3*u) + coeffc4*xx`, floored at 1.9.
- The sampled cosθ is a mixture (L921–975):
  - a **core** `cth = 1 + ln(ea + r·(1−ea))·x` with `x ≈ θ0²` (L893–897) and `ea = exp(−xsi)`. This is a
    Gaussian in θ with projected σ = θ0, **truncated at θ² ≤ 2·xsi·θ0²**;
  - a **power-law tail** with exponent `c = xsi`, weight `1 − prob`;
  - with probability `1 − qprob` (only when qprob < 1) an **isotropic** direction.
  - `qprob = xmeanth/(prob·xmean1 + (1−prob)·xmean2)` with `xmeanth = exp(−τ)`.
  - So only ⟨cosθ⟩ is tied to λ (the transport mean free path); the core width is θ0.
- Simple scattering fallbacks (outside the ion regime, INFERRED): `SimpleScattering()` if
  relloss > 0.5, θ0 > π/6 or xmean1 ≤ 0.999·xmeanth (L867–869, L891, L934).

### 2.3 Charge in the transport cross section (λ), FACT [URB] L172–414

- Electron-equivalent energy (L296–303). INFERRED by exact algebra: this is the electron with the
  same pβ as the projectile.
- Screened-Rutherford branch, used only if that energy ≤ 10 MeV (L316–320):
  `sigma *= chargeSquare*atomicNumber*atomicNumber/(beta2*bg2);`
- High-energy branch (> 10 MeV, which every V3-008 ion and proton case reaches): L399–409
  **overwrite** `sigma` from the `sig0`/`hecorr` tables, **without `chargeSquare`**:
  ```
  c1 = bg2lim*sig0[iZ]*(1.+hecorr[iZ]*(beta2-beta2lim))/bg2;
  ...
  sigma = c1+ratZ*(c2-c1) ;
  ```
  The low-energy factor `sigma *= (1.+0.30/(1.+std::sqrt(1000.*eKineticEnergy)));` follows (L411).
- **OPEN:** whether the z² of an ion is applied to λ outside [URB]. The candidates are the
  `G4VMscModel::GetTransportMeanFreePath` table scaling (`pFactor`, base-particle tables for
  GenericIon) or `G4VMultipleScattering` for `isIon`.
  - Files: `processes/electromagnetic/utils/include/G4VMscModel.hh`, `.../utils/src/G4VMscModel.cc`,
    `.../utils/src/G4VMultipleScattering.cc`, `.../utils/src/G4EmTableUtil.cc`.
  - λ does not enter θ0. It enters xsi (through u and xx), qprob and the true-to-geometric path
    conversion.
  - Scratch emulation (section 6.3): the two λ hypotheses (unscaled vs ×z²) change the core rms per
    1 mm step by about 1 % for C 290 MeV/u.

## 3. G4LindhardSorensenIonModel: stopping construction

- Energy switch, FACT [LS] L82: `fElimit = 2.0*CLHEP::MeV;`, applied to the scaled energy
  `escaled = kinEnergy*pRatio` with `pRatio = CLHEP::proton_mass_c2/mass` (L160, L193, L239, L265).
  - INFERRED (arithmetic with nuclear masses): the switch is at E/A ≈ 2.0·(M/A)/m_p ≈ 1.985 MeV/u
    for 12C and ≈ 1.984 MeV/u for 16O.
- **Above fElimit (escaled > 2 MeV)**, FACT L240–242 (tables) and L292–300 (along step):
  - The tabulated dE/dx is `fBBModel->ComputeDEDXPerVolume(mat, p, kinEnergy, cutEnergy)`, with
    `fBBModel = new G4BetheBlochModel()` (L81).
  - Along step:
    ```
    const G4double q2 = corr->EffectiveChargeSquareRatio(p, mat, e);
    ...
    G4double deltaL0 = 2.0*corr->BarkasCorrection(p, mat, e)*(charge-1.)/charge;
    G4double deltaL = lsdata->GetDeltaL(Zin, gam);
    res = eloss +
      CLHEP::twopi_mc2_rcl2*q2*eDensity*(deltaL+deltaL0)*length/beta2;
    ```
  - `e` is the mid-step energy `max(E − eloss/2, E/2)` (L263).
  - `charge` and `Zin` are the bare PDG charge (L157–158).
  - The effective charge enters through `G4EmCorrections::EffectiveChargeSquareRatio` (L135, L267),
    which is also pushed into the Bragg and BB models (L136–137).
  - INFERRED: `G4EmCorrections` delegates to `G4ionEffectiveCharge` (Ziegler-type; bare charge above
    a Z-dependent energy). OPEN: `processes/electromagnetic/utils/src/G4EmCorrections.cc` and
    `G4ionEffectiveCharge.cc`.
  - **Bloch and Mott:** no separate Bloch or Mott term is added in [LS] (FACT: neither is called).
    The LS ΔL (`G4LindhardSorensenData`, tabulated in (Z, γ)) is the higher-order-in-z correction
    that stands in for them. There is no "high-order corrections" switch in [LS]: `CorrectionsAlongStep`
    is overridden, so `G4BetheBlochModel`'s own ion corrections are not called through this path.
- **Whether `G4BetheBlochModel` returns the ICRU 90 proton array at equal velocity × q² for
  GenericIon in G4_WATER: not visible in the pinned set. OPEN.**
  - The research note `research/research-data-sources.md` §0 reports it as FACT from a scratch read
    of `G4BetheBlochModel.cc`.
  - To pin it, acquire `processes/electromagnetic/standard/src/G4BetheBlochModel.cc` and
    `include/G4BetheBlochModel.hh` (for `fAlphaTlimit` and any upper limit of the ICRU 90 branch),
    plus the ICRU 90 data holder `G4ICRU90StoppingData.cc` (path as in the registered
    `geant4-icru90-stopping-11.4.2` record).
- **At or below fElimit**, FACT L240–241 and L271–291:
  - The tabulated dE/dx is `G4BraggModel` (the proton Bragg parameterisation, L80) × q².
  - Along step, for `Z > 2 && Z <= 80` only, it is replaced by the ion data:
    ```
    if (escaled <= fElimit) {
      // data from ICRU73 or ICRU90
      if(Z > 2 && Z <= 80) {
        res = fIonData->GetDEDX(mat, Z, escaled, G4Log(escaled));
    ```
    - The restricted-loss correction is added when cut < tmax (L282–286), then the result is
      multiplied by the step length.
    - If no data exist (`res == 0`) it falls back to `res = eloss*q2/chargeSquareRatio` (L289–290).
  - **Alpha (Z = 2) never reaches `fIonData`** (FACT L273), and in opt4 alpha does not use this
    model at all (section 1).
- **The exact energy limit below which `G4IonICRU73Data` is consulted is escaled ≤ 2.0 MeV
  (`fElimit`, FACT L82, L271), for 2 < Z ≤ 80.**
  - It is not 2.5 MeV. The value 2.5 MeV is, per the research note, the loader's resampling upper
    bound (`fEmax`), unverified here (section 4).
  - The call into the loader is `fIonData->GetDEDX(mat, Z, escaled, G4Log(escaled))` (L274), after
    `fIonData->Initialise()` (L124). `fIonData = new G4IonICRU73Data()` is created once (L119).
- Delta production above the cut uses BB/Bragg cross sections switched at the same fElimit
  (L193–196). There is also a projectile form factor (L161–170, L363–383). INFERRED: this matters
  only for the cut variant (V4).

## 4. G4IonICRU73Data loader: OPEN, pending materialization of `cf771a5a…`

Everything in this section is **OPEN, pending materialization of `cf771a5a44ea…08d2`** (acquisition
`33881db243e64160be539305e8782dd7`). The research note (`research/research-data-sources.md` §0)
reports, from a scratch read of the same sha256 (to be re-verified by a short follow-up):
- the files are resampled only over 0.025–2.5 MeV scaled (`fEmin`, `fEmax`);
- `$G4LEDATA/ion_stopping_data/icru90/z<Z>_<material>.dat` for Z ≤ 18 in G4_AIR/G4_GRAPHITE/G4_WATER,
  "density coefficient 0.998 for water";
- otherwise `icru73/z<Z>_<material>.dat`, and element files `z<Z>_<Ztarget>.dat`.

To verify in the follow-up:
1. The directory choice (icru90 vs icru73) for Z = 6 and Z = 8 in G4_WATER.
2. `fEmin`/`fEmax` and the number of resampling points.
3. The G4_WATER special case and the meaning of the 0.998 factor (a density ratio or a stopping
   scale).
4. Any Z_eff² or q² scaling inside `GetDEDX`. (FACT [LS] L274 passes `escaled` and `Z` and adds no
   q² to `res`, which suggests the files are per-ion stopping at the scaled energy.)
5. Any erratum or revision comment.

Whatever the loader does, by section 3 it acts only for escaled ≤ 2.0 MeV. Its span cannot widen the
energy range in which TOPAS uses these files. The P5 span question (whether the *files* cover
2–430 MeV/u) is a G4EMLOW data question, separate from the loader's resampling window.

## 5. Verdicts on the decision's claims, and the lineage consequence

| Claim (decision 0043 Context; plan §B lineage) | Verdict | Basis |
|---|---|---|
| GenericIon (C, O) stopping model is `G4LindhardSorensenIonModel` | **CONFIRMED** (for the default constructor; TOPAS instantiation confirmed at P4) | FACT [OPT4] L316–322, L137 |
| "ICRU 90 proton array × q_eff² + LS ΔL + Barkas above ~2 MeV/u" | **PARTLY CONFIRMED** | CONFIRMED from [LS]: the switch at escaled = 2.0 MeV (≈ 1.985 MeV/u C, 1.984 O); `G4BetheBlochModel` above it; effective charge from `G4EmCorrections::EffectiveChargeSquareRatio`; LS ΔL from `G4LindhardSorensenData`; Barkas term `2·BarkasCorrection·(z−1)/z`, both multiplied by q_eff² (not z²). NOT confirmed from the pinned set: that the BB dE/dx for GenericIon in G4_WATER is the **ICRU 90 proton array** at equal velocity (OPEN, `G4BetheBlochModel.cc`; FACT only in the research note's scratch read). Wording correction: there is no separate Bloch or Mott term; the LS ΔL replaces them |
| "ICRU 73/90 ion tables only ≤ 2.5 MeV/u" | **CONFIRMED in substance, with a tighter limit**: ion tables are consulted only for escaled ≤ **2.0 MeV** (≈ 1.98 MeV/u for C/O) and only for 2 < Z ≤ 80, i.e. never for alpha | FACT [LS] L82, L265, L271–274 |
| Alpha = ICRU 90 alpha array | **OPEN** (not in the pinned set) | needs `G4EmBuilder.cc`, `G4ionIonisation.cc`, `G4BetheBlochModel.cc`, `G4BraggIonModel.cc` |
| Opt4 msc for GenericIon/alpha is `G4UrbanMscModel` | **INFERRED** (strong), OPEN until `G4hMultipleScattering.cc` / `G4EmBuilder.cc` or the P4 printout | [OPT4] L161, L323, L328 (no model set) |
| Urban θ0 is Highland-like `13.6·|q|·√y/(βcp)·(c1 + c2 ln y)` without z² in the log | **CONFIRMED** | FACT [URB] L1023–1027; y = t/X0; c1, c2 = coeffth1/2(Zeff of the material) |

**Lineage consequence pre-declared in decision 0043 (I3 class, V1 C/O lineage).**
- The trigger was "ion tables used above 2.5 MeV/u". It is **not met**: the tables are used only at
  ≤ 2.0 MeV scaled.
- Therefore **no amendment is triggered**. I3 keeps the class "ion-specific tabulated (largely
  independent of TOPAS above ~2 MeV/u)", and the V1 C/O lineage statement stands.
- Two wording corrections are recommended for the freeze text:
  - "≤ 2.5 MeV/u" becomes "≤ 2.0 MeV scaled energy (≈ 1.98 MeV/u)";
  - "+ Geant4 Barkas" becomes "+ Geant4 Barkas term 2·B·(z−1)/z, both ΔL terms weighted by q_eff²".
- The "ICRU 90 proton array" part of the lineage stays "FACT (research note)" until
  `G4BetheBlochModel.cc` is pinned.

## 6. P2b preparation: frozen Urban θ0, and what δ_pred needs

### 6.1 θ0 as a python function (frozen from [URB] L982–1028, L1115–1125)

```
def urban_theta0(z, M, E1, t, E2=None, X0=X0_WATER, Zeff=ZEFF_WATER):
    # z: bare charge (|PDG charge|), M: projectile mass [MeV], E1: pre-step kinetic energy [MeV],
    # t: true step length [same unit as X0], E2: kinetic energy used by SampleScattering (L794-800)
    inv = (E1 + M) / (E1 * (E1 + 2*M))                    # 1/(p beta c) [1/MeV]
    if E2 is not None and E2 != E1:
        inv = sqrt(inv * (E2 + M) / (E2 * (E2 + 2*M)))     # geometric mean (code evaluates E2 first; symmetric)
    y = t / X0
    w = exp(log(Zeff) / 6)
    facz = 0.990395 + w*(-0.168386 + w*0.093286)
    c1 = facz * (1 - 8.7780e-2/Zeff)
    c2 = facz * (4.0780e-2 + 1.7315e-4*Zeff)
    return 13.6 * abs(z) * sqrt(y) * inv * (c1 + c2*log(y))   # [rad]
```

Constants and their status:
- 13.6 MeV, the `facz`/c1/c2 constants, the √y·ln y form and the geometric-mean rule: **FACT**.
- **X0_WATER (G4_WATER radiation length): INFERRED 36.083 cm.** This is the Geant4 material printout
  value as recalled.
  - OPEN: confirm from the P4 run log (material table printout) or `source/materials/src/G4Material.cc`.
  - Sensitivity: θ0 ∝ X0^(−1/2)·(c1 + c2 ln(t/X0)).
- **ZEFF_WATER = `G4IonisParamMat::GetZeffective()` for G4_WATER: OPEN.**
  - INFERRED, recalled: the atom-number-density-weighted mean Z, (2·1 + 8)/3 = 10/3.
  - Settled by `source/materials/src/G4IonisParamMat.cc`.
  - Sensitivity: Zeff = 10/3 gives c1 = 0.8996113, c2 = 0.0382116. Zeff = 6.6 (electron-weighted)
    gives c1 = 0.9223267, c2 = 0.0391877, i.e. **θ0 is +2.5 % higher**. This must be closed before
    δ_pred is frozen.
- Particle values: z = bare charge (INFERRED, see 2.1). M = Geant4 ion mass (`G4NucleiProperties`;
  the nuclear masses from `projectiles.py` agree to ≪ 1e-4, which is negligible).

Test vectors from my scratch implementation. They are for checking an implementation, not
acceptance values (Zeff = 10/3, X0 = 36.083 cm, E2 = E1, nuclear masses p 938.272, α 3727.379,
12C 11174.862, 16O 14895.079 MeV):

| Case | t | θ0 (rad) |
|---|---|---|
| p 150 MeV | 1 mm | 1.729130e-3 (Zeff 6.6: 1.772622e-3) |
| 4He 150 MeV/u | 1 mm | 8.649455e-4 |
| 12C 290 MeV/u | 1 mm | 4.724671e-4 |
| 12C 290 MeV/u | 0.1 mm | 1.299207e-4 |
| 16O 350 MeV/u | 1 mm | 3.995826e-4 |
| 12C 80 MeV/u | 0.5 mm | 1.067592e-3 |
| 12C 290 MeV/u, E2 = 0.99·E1 | 1 mm | 4.746029e-4 |

### 6.2 Why "Fermi–Eyges integration of the Urban θ0 scattering power" is not sufficient (INFERRED, method finding)

1. **θ0 is not the projected rms of an Urban step.**
   - The core is a Gaussian truncated at θ ≤ √(2·xsi)·θ0 (xsi ≈ 2.4–2.5 for these cases), mixed
     with a power-law tail (weight 1 − prob ≈ 13–15 %).
   - In scratch sampling of [URB] L834–978, the projected rms inside ±5θ0 is ≈ 1.13–1.16·θ0 per
     step. The full rms is 1.26–1.5·θ0, with heavy-tail noise.
   - The lateral σ graded by V2 is a ±10 mm-window second moment of the multi-step sum. It is
     therefore neither θ0-based nor full-rms-based.
2. **Strong step-length dependence.**
   - θ0²/t ∝ (c1 + c2 ln(t/X0))². Between t = 0.1 mm and 1 mm this per-length "scattering power"
     changes by a factor 0.756.
   - Urban does not restore additivity for heavy particles: qprob ≥ 1 in practice (τ ~ 1e-9 to
     1e-6), so ⟨cosθ⟩ is not pinned to exp(−t/λ).
   - Ion steps in TOPAS are not set by msc (section 1, `fMinimal` branch). They are set by the
     ionisation step function (0.1·R + 1 µm for GenericIon, i.e. cm-scale), by delta production
     (cut-dependent: the matched and 0.05 mm variants differ), by the IDD/3-D scoring bins (0.1 mm and
     1 mm depth, if TOPAS realises them as step-limiting boundaries; OPEN, TOPAS behaviour) and by
     any `MaxStepSize`.
   - **δ_pred is therefore a property of the TOPAS bundle design, not of Urban alone.**
3. **Lateral displacement on or off (OPEN, section 1)** matters only for long steps. With steps
   ≲ 0.5 mm the discretisation bias is ≲ 1 % in σ.
4. λ for ions (section 2.3, OPEN) shifts xsi, hence the in-window rms, by ≈ 1 %.

**Recommendation for the P2 method (orchestrator decision; plan amendment before freeze):**
- Compute δ_pred with a committed numpy **emulator of the Urban per-step sampler** ([URB]
  L834–978 + L982–1028, coefficients from L1115–1132, `sig0`/`hecorr` read from the pinned file),
  run over the ion's slowing-down with a **declared TOPAS step-length model**.
- Score it with the same ±10 mm Sheppard-corrected σ estimator as V2.
- Do not use a Gaussian Fermi–Eyges integral of θ0.
- Freeze the step-length model from the V2 bundle geometry. Best: make V2's TOPAS step length
  explicit (a dedicated V2 run with only the 3-D scorer and a fixed `MaxStepSize` in water, TOPAS
  parameter name TO VERIFY).
- At P4, record the primary-ion step-length statistics from a few histories with tracking verbose
  (execution only, no dose).

### 6.3 Scratch magnitude estimate (estimate only; not a δ_pred value; not frozen)

Scratch emulation: constant energy, Zeff = 10/3, λ ×z², no lateral displacement, 8 cm of water.
It gives the Urban lateral σ (3-MAD window) relative to a Highland-slab Fermi–Eyges value:

| Case | 0.1 mm steps | 0.5 mm steps | 1 mm steps |
|---|---|---|---|
| 12C 290 MeV/u | 0.835 | 0.876 | 0.890 |
| p 150 MeV | 0.857 | 0.885 | 0.895 |

- The proton figure is consistent with the literature statement in `docs/research/validation-evidence.md`
  L12 (Makarova, Gottschalk & Sauerwein: Geant4 Urban about 8 % low for low-Z targets).
- Decision 0043 expects the Molière ion width to exceed z/(pv) scaling by +5 % (He) and +8 % (C/O).
  Combining the two, **|δ_pred| for C/O plausibly lies at 10–20 %** (TOPAS narrower). This is well
  outside the decision's "5–8 %" expectation and likely in **Case C (|δ_pred| > 10 %)**.
  He is plausibly near the 10 % boundary.
- The proper committed computation (6.2) must replace this estimate. The figure is given only
  because it bears on whether Case B is reachable at all.

### 6.4 P2b (proton e_δ) applicability: OPEN, likely invalid as defined

- If opt4 protons use WentzelVI + single scattering (INFERRED strong; section 1), the archived T15
  proton TOPAS width was **not** produced by Urban.
- An Urban-based proton δ_pred compared with T15 would then not check the δ_pred method. A small
  e_δ would be coincidental, and a large one uninformative.
- Options for the orchestrator, after `G4EmBuilder.cc` (or the T15 run log's msc printout) settles
  the proton model:
  - **(A)** Declare P2b inapplicable. e_δ is undetermined, so under the plan's rule Case B is
    unavailable and V2 is Case C for every ion with |δ_pred| > 2 %. This is honest and costs
    nothing; lateral MC evidence becomes a declared gap.
  - **(B)** Replace the T15 comparison with a proton TOPAS run that uses Urban for protons. Candidates:
    a Geant4 EM constructor where protons get Urban (for example opt3, if
    `G4EmStandardPhysics_option3.cc` shows `ConstructCharged(..., false)`; TO VERIFY) or another
    TOPAS-level setting. Use the same geometry and step design as the ion V2 bundle, run before any
    ion TOPAS run. Cost: ≈ 3 × 6 min of REF time and one amendment. It tests the emulator and
    step-model method on the reference engine's actual Urban code path.
- Recommendation: settle the proton model first. If WentzelVI is confirmed, prefer (B) if a
  source-verified Urban-for-protons configuration exists that leaves water, X0 and Zeff unchanged;
  otherwise (A).

## 7. Companion files to acquire (exact upstream paths, tag v11.4.2)

| Point settled | Path |
|---|---|
| Loader internals (section 4) | `source/processes/electromagnetic/standard/src/G4IonICRU73Data.cc` (acquired, not materialized) |
| Default msc model of `G4hMultipleScattering` | `source/processes/electromagnetic/standard/src/G4hMultipleScattering.cc` |
| Alpha, proton and d/t assignments; `isWVI` | `source/physics_lists/constructors/electromagnetic/src/G4EmBuilder.cc` (+ `include/G4EmBuilder.hh`) |
| Ion step limit, RangeFactor, lateral displacement | `source/processes/electromagnetic/utils/src/G4EmTableUtil.cc`, `.../utils/src/G4VMscModel.cc`, `.../utils/include/G4VMscModel.hh`, `.../utils/src/G4VMultipleScattering.cc`, `.../utils/src/G4EmParameters.cc` |
| `SetParticle`, inline helpers of Urban | `source/processes/electromagnetic/standard/include/G4UrbanMscModel.hh` |
| Zeff of G4_WATER used by Urban | `source/materials/src/G4IonisParamMat.cc` |
| X0 of G4_WATER | `source/materials/src/G4Material.cc` (or the P4 material printout) |
| ICRU 90 proton/alpha branch for GenericIon/alpha | `source/processes/electromagnetic/standard/src/G4BetheBlochModel.cc` + `include/G4BetheBlochModel.hh`; alpha low-energy: `.../standard/src/G4BraggIonModel.cc`; `.../standard/src/G4ionIonisation.cc` |
| Effective charge | `source/processes/electromagnetic/utils/src/G4EmCorrections.cc`, `.../utils/src/G4ionEffectiveCharge.cc` (directory TO VERIFY) |
| Low-energy dE/dx table below fElimit | `source/processes/electromagnetic/standard/src/G4BraggModel.cc` (informative) |
| LS ΔL table | `source/processes/electromagnetic/standard/src/G4LindhardSorensenData.cc` (already in the plan list; needed for I3c only) |
| Proton-Urban alternative for P2b (B) | `source/physics_lists/constructors/electromagnetic/src/G4EmStandardPhysics_option3.cc` |

Priority for closing P2/P2b:
1. `G4EmBuilder.cc`, `G4hMultipleScattering.cc`, `G4IonisParamMat.cc`;
2. `G4EmTableUtil.cc`, `G4VMscModel.cc/.hh`, `G4EmParameters.cc`, `G4UrbanMscModel.hh`;
3. `G4BetheBlochModel.cc/.hh`, `G4IonICRU73Data.cc` (lineage wording only).

Run-time evidence at P4 (execution only) can replace several of these:
- the msc block per particle (model name, RangeFactor, stepLimType, latDisp);
- the ionisation model names;
- the G4_WATER material table (RadL).

## P2 follow-up (2026-11-01)

Written 2026-11-01 by the research specialist at worktree head `3d2547ad`. Source reading only: no
transport output, no TOPAS output, no tabulated reference values. Markers as in section 0 (FACT
(quoted), INFERRED, OPEN). Each OPEN item names the file or run-time line that settles it.

### F0. Provenance (sha256 re-computed here; all match)

| Short name | File (upstream path under `source/`, tag v11.4.2) | sha256 | Bytes | Lines | Acquisition id |
|---|---|---|---|---|---|
| [EMB] | `physics_lists/constructors/electromagnetic/src/G4EmBuilder.cc` | `40782f5c5bc0e66f9b4b635a08d8b0c236b8576ad5b885afca04b679bb18d3d1` | 15 843 | 466 | `3e362c8b…` |
| [HMSC] | `processes/electromagnetic/standard/src/G4hMultipleScattering.cc` | `a664386984f22898ca8bcaf0f38326479e1e463de88446d41fd465dac83f50ea` | 4 288 | 106 | `78bd2972…` |
| [VMSC] | `processes/electromagnetic/utils/src/G4VMscModel.cc` | `87b9f203542f2412b3bf327b2cc41eb0c5295a92fabf703009b66f7d37b97a69` | 8 224 | 221 | `514d42e3…` |
| [IPM] | `materials/src/G4IonisParamMat.cc` | `990c6f315a3b6951253eb251024ba439f8edf4e151994c33991510645807849d` | 21 112 | 622 | `4199285b…` |
| [ICRU73] | `processes/electromagnetic/standard/src/G4IonICRU73Data.cc` | `cf771a5a44ea936c56593b7dd2d67bb4157345eb316815f60880986404fe08d2` | 13 701 | 445 | `33881db2…` |

The upstream paths are the coordinator's acquisition paths; the files do not state their own
paths. Geant4 Software Licence; short excerpts are quoted with attribution; nothing is
redistributed as data.

### F1. The msc model for GenericIon and alpha in option 4, and their parameters

**Default model of `G4hMultipleScattering`: Urban. FACT [HMSC] L54–58, L73–79:**
```
G4hMultipleScattering::G4hMultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName)
{
  isInitialized = false;
}
...
void G4hMultipleScattering::InitialiseProcess(const G4ParticleDefinition*)
{
  if(isInitialized) { return; }
  if(nullptr == EmModel(0)) { SetEmModel( new G4UrbanMscModel() ); }
  AddEmModel(1, EmModel(0));
  isInitialized = true;
}
```
- FACT: the constructor and `InitialiseProcess` set **no** step-limit type, RangeFactor, skin or
  lateral-displacement value. These come from the model (F2) or from base-class code (OPEN, F2).
- FACT [HMSC] L85–87: the process prints its parameters at initialisation:
  `"      RangeFactor= " … ", stepLimType: " … ", latDisp: "`. **P4 should archive these lines** for
  "ionmsc" (GenericIon) and for the alpha msc process.

**GenericIon: Urban. FACT.**
- [OPT4] L161 creates `new G4hMultipleScattering("ionmsc")`. [OPT4] never calls `SetEmModel` on it.
  [OPT4] L323 registers it for GenericIon. By [HMSC] L76 the model is `G4UrbanMscModel`.
- [OPT4] registers no `G4CoulombScattering` for GenericIon (L316–325: only `hmsc`, `ionIoni`,
  `pnuc`).
- The same "ionmsc" object goes to `ConstructCharged` (L328). There it is registered for deuteron
  and triton ([EMB] L123–129) and, if `isHEP`, for heavy hadrons and hyper-nuclei ([EMB] L287–297).

**Alpha: Urban, through its own msc process (not "ionmsc"). FACT [EMB] L131–136:**
```
  part = G4Alpha::Alpha();
  ph->RegisterProcess(new G4hMultipleScattering(), part);
  ph->RegisterProcess(new G4ionIonisation(), part);
  if( nucStopping != nullptr ) {
    ph->RegisterProcess(nucStopping, part);
  }
```
- No `SetEmModel` call, so the model is Urban ([HMSC] L76). No `G4CoulombScattering` for alpha.
- He3 is identical (L138–143).
- INFERRED: the default process name is "msc" (`G4hMultipleScattering.hh`). P4 should therefore
  search the printout for the alpha block under "msc", not "ionmsc".
- Alpha ionisation is `G4ionIonisation` with **no model set** ([EMB] L133). The default models
  (Bragg-ion type below the switch, `G4BetheBlochModel` above) are INFERRED. They are settled by
  `processes/electromagnetic/standard/src/G4ionIonisation.cc` and the P4 printout. The claim
  "alpha = ICRU 90 alpha array" stays OPEN until `G4BetheBlochModel.cc` and
  `G4BraggIonModel.cc` are pinned.

**Path from option 4 into `G4EmBuilder` (FACT).**
- [OPT4] L328: `G4EmBuilder::ConstructCharged(hmsc, pnuc);` is called with two arguments.
- [EMB] L232–234: `ConstructCharged(G4hMultipleScattering* hmsc, G4NuclearStopping* nucStopping, G4bool isWVI)`.
- [EMB] L284: `ConstructIonEmPhysics(hmsc, nucStopping);`
- The default value of `isWVI` is in `include/G4EmBuilder.hh`: OPEN (INFERRED `true`).
- **`isWVI` does not affect alpha or GenericIon** (FACT): it is used only for muons (L243–244),
  pi/K (L272–275) and p/pbar (L278).
- Protons get `G4WentzelVIModel` + `G4CoulombScattering` if `isWVI` (FACT [EMB] L190–195, L203);
  otherwise Urban. The proton model therefore stays INFERRED (WentzelVI) until
  `G4EmBuilder.hh` or the P4 proton printout settles it. This does not change P2b's status
  ("very probably inapplicable").

**Nuclear stopping (FACT).**
- [OPT4] L134 sets `SetMaxNIELEnergy(1*CLHEP::MeV)`, and L164–169 then create `G4NuclearStopping`.
- It is registered for GenericIon (L325) and alpha ([EMB] L134–136).
- INFERRED: it acts only below 1 MeV kinetic energy. ionmc deposits locally below 2 MeV/u, so this
  is irrelevant to V1/V2.

**Step-limit type, RangeFactor and lateral displacement for GenericIon and alpha: partly FACT,
final value OPEN.**
- FACT [VMSC] L96–113, called from [URB] L139:
  ```
  void G4VMscModel::InitialiseParameters(const G4ParticleDefinition* part)
  {
    if(IsLocked()) { return; }
    G4EmParameters* param = G4EmParameters::Instance();
    if(std::abs(part->GetPDGEncoding()) == 11) {
      steppingAlgorithm = param->MscStepLimitType();
      facrange = param->MscRangeFactor();
      latDisplasment = param->LateralDisplacement();
    } else {
      steppingAlgorithm = param->MscMuHadStepLimitType();
      facrange = param->MscMuHadRangeFactor();
      latDisplasment = param->MuHadLateralDisplacement();
    }
    skin = param->MscSkin();
    facgeom = param->MscGeomFactor();
    facsafety = param->MscSafetyFactor();
    lambdalimit = param->MscLambdaLimit();
  }
  ```
- Consequences for the model path (FACT):
  - [OPT4] L126 (`fUseSafetyPlus`) and L128 (`SetMscRangeFactor(0.08)`) apply to e± only.
  - Ions and alpha take the MuHad values.
  - [OPT4] L129 `SetMuHadLateralDisplacement(true)` gives `latDisplasment = true` for them on
    this path.
  - `skin = 3` ([OPT4] L127) is copied for every particle. [URB] uses it only in the
    `fUseDistanceToBoundary` branch.
- The default values of `MscMuHadStepLimitType()` and `MscMuHadRangeFactor()` are OPEN
  (`processes/electromagnetic/utils/src/G4EmParameters.cc`, `SetDefaults`). INFERRED: `fMinimal`
  and 0.2.
- OPEN, competing path: a process-level override for particle type "nucleus" (INFERRED: `fMinimal`,
  `latDisp = false`, factor 0.2), in `processes/electromagnetic/utils/src/G4EmTableUtil.cc`
  (`PrepareMscProcess`) or `.../utils/src/G4VMultipleScattering.cc`. Whether it reaches the model
  after `InitialiseParameters` is not visible in the pinned set.
- FACT [URB] L485–491: even with `latDisplasment = true`, Urban switches displacement off for a step
  when `currentRange*doverrb < presafety`, i.e. far from boundaries.
- **Run-time settlement (P4):**
  - the [HMSC] L85–87 process line;
  - the per-model `G4VMscModel::DumpParameters` line, FACT [VMSC] L124–126:
    `"StepLim=" … " Rfact=" … " Gfact=" … " Sfact=" … " DispFlag:" … " Skin=" … " Llim=" … " mm"`.
- None of this changes the model identity, and none of it affects the V2 branch.

### F2. `G4VMscModel` defaults and the charge-square question ("λ gets z²")

- FACT [VMSC] L59–68, constructor: `lambdalimit(1.*CLHEP::mm)`, `geomMin(1.e-6*CLHEP::mm)`,
  `geomMax(1.e50*CLHEP::mm)`, `steppingAlgorithm(fUseSafety)`,
  `dedx = 2.0*CLHEP::MeV*CLHEP::cm2/CLHEP::g`. Unless the model is locked, `InitialiseParameters`
  (F1) overwrites the step limit, facrange, latDisplasment, skin, facgeom, facsafety and
  lambdalimit.
- OPEN: the in-class initialisers of `facrange`, `latDisplasment`, `skin` and `dtrl` (P2 §2.1
  INFERRED `dtrl` = 0.05). Settled by `processes/electromagnetic/utils/include/G4VMscModel.hh`.
- **Charge square in [VMSC] (FACT, whole file read).**
  - q² appears only in the fallbacks used when no ionisation process is attached:
    - `x = dedx*q*q` (L146–147, L161–162);
    - `localrange = kinEnergy/(dedx*q*q*…)` (L177–178, L195–196);
    - `e -= (localrange - range)*dedx*q*q*…` (L214–215).
  - **No charge scaling of the transport cross section λ is applied in `G4VMscModel.cc`.**
- FACT [URB] L399–409 (P2 §2.3): the Urban high-energy branch has no `chargeSquare`.
  - INFERRED arithmetic: the electron with the same pβ is above 10 MeV for He/C/O at all energies
    ≥ 2 MeV/u (C at 2 MeV/u: pβc ≈ 48 MeV). So this branch is the one used.
- **Still OPEN:** whether a z² enters outside these files, through:
  - the inline `GetTransportMeanFreePath` ([URB] L448, L726, L842 call it), in `G4VMscModel.hh`;
  - the λ table built for GenericIon and reused for C/O, in `G4EmTableUtil.cc` and
    `G4VMultipleScattering.cc`.
- Impact:
  - λ does not enter θ0;
  - it shifts the in-window rms by ≈ 1 % (P2 §6.3 scratch);
  - it is irrelevant to the V2 branch.

### F3. `G4IonICRU73Data`: directories, naming, span, the water special case, erratum

**Data directory (FACT L158–161):** `G4EmParameters::GetDirLEDATA() + "/ion_stopping_data/"`.
INFERRED: this is `$G4LEDATA`.

**Directory choice for G4_WATER under option 4 (FACT).**
- [OPT4] L131 sets `SetUseICRU90Data(true)`.
- [ICRU73] L88: `namesICRU90[3] = {"G4_AIR","G4_GRAPHITE","G4_WATER"}`.
- [ICRU73] L195–198: a match calls `ReadMaterialData(mat, densityCoef[j], true)`.
- [ICRU73] L261–264 and L283:
  ```
  if(useICRU90 && Z <= 18) {
    ost << "90";
  } else {
    ost << "73";
  ...
  ost << "/z" << Z1 << "_" << name << ".dat";
  ```
- Hence for C and O in G4_WATER, Geant4 11.4.2 + option 4 reads
  **`ion_stopping_data/icru90/z6_G4_WATER.dat` and `icru90/z8_G4_WATER.dat`**. It never reads
  `icru73/z{6,8}_G4_WATER.dat`. For Z ≤ 18, Z1 = Z and there is no Z² rescaling (L259–262).
- **If `icru90/z8_G4_WATER.dat` is absent from G4EMLOW 8.8:**
  - `RetrieveVector(ost, false)` returns nullptr without a warning (L404–416);
  - `GetDEDX` returns 0 (L147);
  - LS falls back to Bragg × q² (P2 §3, [LS] L289–290).
  - **P5 should add `icru90/z8_G4_WATER.dat` to the member list** and record whether it exists.

**The ICRU 73 named-material list is inoperative (FACT by C++ semantics).**
- [ICRU73] L56–87: the 31 initialisers of `namesICRU73[31]` contain no commas (checked: 0 commas in
  L57–87). Adjacent string literals concatenate, so `namesICRU73[0]` is one concatenated string and
  `namesICRU73[1..30]` are empty.
- The name match at L208–219 therefore never succeeds.
- Consequence: without the ICRU 90 flag, G4_WATER would be built from element files
  `icru73/z{Z}_{1,8}.dat` by mass-fraction additivity (L222–238, L310–336, L385–389), not from
  `icru73/z{Z}_G4_WATER.dat`.
- For V3-008 this means the `icru73/z{6,8}_G4_WATER.dat` files, if used as the I3 reference, are
  not read by TOPAS at any energy.
- INFERRED: this is an upstream code defect. It does not affect option 4 + G4_WATER.

**Energy span (FACT).**
- L100–102: `fEmin = 0.025*CLHEP::MeV; fEmax = 2.5*CLHEP::MeV;`, with
  `fNbins = fNbinsPerDecade*G4lrint(std::log10(fEmax/fEmin))` (2 decades). The value of
  `fNbinsPerDecade` is OPEN (`include/G4IonICRU73Data.hh`).
- Each file is resampled onto this log grid (L432–437).
- Below fEmin: `(*v)[0]*std::sqrt(e/fEmin)` (L148–149).
- It is queried only at escaled ≤ 2.0 MeV (P2 §3, [LS] L271–274). The loader window therefore says
  nothing about the files' own span; that is the P5 (i) data question.

**Units and the water scale factor (FACT L89, L198, L286–288).**
```
const G4double densityCoef[3] = {0.996, 1.025, 0.998};
...
const G4double fact = coeff *
  mat->GetDensity() * CLHEP::MeV * 1000 * CLHEP::cm2 / CLHEP::g;
v->ScaleVector(CLHEP::MeV, fact);
```
- G4_WATER file values are multiplied by **0.998** × density × 1000 MeV cm²/g.
- INFERRED: the file values are in MeV cm²/mg on a MeV (proton-scaled) energy axis; P5 (ii)
  confirms this.
- The meaning of 0.998 is not documented in the file (OPEN). It is not a plain density ratio
  (graphite takes 1.025).
- It acts only at ≤ 2.0 MeV scaled, i.e. on ≲ 0.03 mm of residual C/O range: negligible for V1.

**Erratum (FACT negative).**
- The file has no comment on an erratum, a revision or an I-value (grep for "errat", "revis",
  "67.2", "78": none).
- The only provenance comment is "data from ICRU73 or ICRU90" ([LS] L272).
- The erratum question stays with P5 (iii) (data and documentation).

**Effective-charge scaling (FACT).**
- For ICRU 90 files with Z ≤ 18 there is none.
- For ICRU 73 projectiles missing from `zdat`, the nearest tabulated Z1 is rescaled by
  (Z/Z1)² (L265–278, L297–303). C and O are in `zdat` (L91).

### F4. `G4IonisParamMat`: Zeff for Urban, X0, and the I value of G4_WATER

**Zeff used by Urban (`GetZeffective()`, [URB] L1115): FACT [IPM] L412–448:**
```
      const G4double weight = theAtomicNumDensityVector[iel];
      norm += weight;
      z += element->GetZ() * weight;
...
    z *= norm;
...
  fZeff = z;
```
- This is the **atom-number-density-weighted mean Z**.
- For H2O (atom ratio 2:1, any H2O definition): **Zeff = (2·1 + 1·8)/3 = 10/3 = 3.3333**.
- Hence c1 = 0.8996113 and c2 = 0.0382116 in P2 §6.1 (the 10/3 column; the 6.6 alternative is
  rejected).
- INFERRED: `GetZeffective()` is the inline accessor returning `fZeff` (`materials/include/G4IonisParamMat.hh`).
- FACT L382–385: a different, mass-fraction-weighted Zeff (≈ 7.22 for water) is used only in
  `ComputeFluctModel` (energy-loss fluctuations), not by Urban msc.

**X0 (radiation length) is not in this file. FACT.**
- The `fX0density` symbol in [IPM] (L99, L104, L372, L473) is the Sternheimer density-effect
  parameter, not the radiation length. The two must not be confused.
- `currentRadLength = couple->GetMaterial()->GetRadlen()` ([URB] L850) comes from
  `materials/src/G4Material.cc` (`ComputeRadiationLength`) with
  `materials/src/G4Element.cc` (Tsai factors). **OPEN** until these are pinned or the P4 material
  printout ("RadL") is archived.
- INFERRED cross-check: Tsai's formula with Geant4's light-element constants gives
  X0(H2O) ≈ 36.08 g/cm², i.e. 36.08 cm at the G4_WATER density 1.0. This is consistent with P2's
  36.083 cm.

**I value (new finding; for the lineage wording).**
- FACT [IPM] L556–619: a chemical formula `"H_2O"` maps to 75.0 eV. The name is at L581, index 38;
  the value at L606, index 38. This overrides the density-effect data.
- INFERRED: G4_WATER carries the formula "H_2O" and is then set to 78 eV by
  `G4NistMaterialBuilder` via `SetMeanExcitationEnergy` (FACT L456–474 corrects the density
  effect when this happens).
- The plan's "G4_WATER (I = 78 eV)" stays INFERRED until the P4 material printout (or
  `materials/src/G4NistMaterialBuilder.cc`) shows 78 eV.

### F5. V2 branch

**Applies: "If Urban: V2 informative (Case C), the raw ratio and expected direction reported." —
FACT at the Geant4-source level.**
- In Geant4 11.4.2 `G4EmStandardPhysics_option4`, GenericIon and alpha each get a
  `G4hMultipleScattering` with no model set, hence `G4UrbanMscModel`.
- Neither gets a `G4CoulombScattering` process (F1: [OPT4] L161, L316–325; [EMB] L131–136;
  [HMSC] L76).
- The second branch (WentzelVI/single scattering) is excluded for ions by source.
- Residual confirmation is at the TOPAS level only. The P4 scorer-free printout must show the Urban
  msc model for GenericIon ("ionmsc") and alpha, and no CoulombScat process for them, which
  confirms that OpenTOPAS's `g4em-standard_opt4` alias instantiates the unmodified constructor and
  that no region model activator applies.
- A contradicting P4 printout is a STOP for this branch assignment (amendment before any ion TOPAS
  run). It is not a re-interpretation.
- The remaining OPEN items (step-limit type, RangeFactor, latDisp, λ z², `dtrl`, X0) affect only
  an Urban emulator, which r5 withdrew. They do not affect Case C.

### F6. Note for P4 (outside P2 scope; flagged, not investigated further)

- The OpenTOPAS documentation shows `d:Ph/Default/EMRangeMax = 600. MeV # maximum for EM tables`
  and states that Geant4 ≥ 11.3 requires at least 600 MeV.
- If OpenTOPAS applies such a default, the alpha tables (in alpha kinetic energy; INFERRED from
  `G4ionIonisation`) would end below He 150–250 MeV/u = 600–1000 MeV total.
- P4 should archive the `printParameters` maximum kinetic energy. The He bundles should set
  `d:Ph/Default/EMRangeMax` ≥ 1.1 GeV explicitly, unless the printout shows a larger value.
- GenericIon tables (scaled energy ≤ 433 MeV for O 430 MeV/u) are inside 600 MeV.
