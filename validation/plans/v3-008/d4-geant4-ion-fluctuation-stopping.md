# D-4: Geant4 11.4.2 ion fluctuation and ion stopping paths read from source (V3-008, amendment R2-STOP3 item 1)

## 0. Status

Report-only diagnostic for attribution (plan `validation/plans/v3-008-acceptance.md`, amendment R2-STOP3 item 1,
diagnostic D-4). Dated 2026-11-13, research specialist, worktree head `ed5f1d4f`. **No ionmc change follows from this
record and nothing is re-graded.** Every grade of the main TD/TL rows and of the held-out rows TD-H/TL-H stands as
recorded. Plan D, STOP 3: no ionmc parameter is fitted to the reference.

This record replaces the preliminary D-4 reading of 2026-11-11, which was made from files fetched outside the
acquisition tool. The source files below were re-acquired through `acquire_reference_data` and re-read.

No transport output and no scored TOPAS output was read. The only run-time item cited is the EM configuration printout
(process/model table) of an archived TOPAS stdout (section 2.6). No Geant4 tabulated data value and no stopping-power
or range value is transcribed. Source constants (thresholds, model parameters) are quoted with file:line where they
are code. Untrusted files were read with `cat`/`sed`/`grep` only, never compiled, executed or imported.

Labels:
- **ESTABLISHED**: read in a pinned file below, at file:line.
- **INFERRED**: derived by algebra or planning arithmetic, or from physics literature not in the pinned set.
- **TO VERIFY**: open; the settling evidence is named.

## 1. Pinned sources

All files are Geant4 tag **v11.4.2** (GitHub), Geant4 Software License, role exploratory. They are materialized
verbatim at `.ionmc-cache/reference/<sha256>` in the V3-008 worktree. sha256, byte count and line count were
recomputed on 2026-11-13 and all match. **Every line citation below refers to exactly these files.** Paths are
upstream paths under `source/`.

| Tag | Path | sha256 | Bytes | Lines | Acquisition id |
|---|---|---|---|---|---|
| [IF] | `processes/electromagnetic/standard/src/G4IonFluctuations.cc` | `a888247514685c7d49b2ca5d5a0b9fc12dc995a4b4895e39152307f5c1cfc65e` | 15 346 | 439 | `eaf2c86693e845b3ad2774b93b0d4d52` |
| [IFh] | `processes/electromagnetic/standard/include/G4IonFluctuations.hh` | `962a89f116074caf779092591e7dd1042ee6ed15eb806c2555859f480b8bbff0` | 4 063 | 122 | `14d1646b50834ebda79ce113a41942cb` |
| [UF] | `processes/electromagnetic/standard/src/G4UniversalFluctuation.cc` | `c4b1f80c16e3668996801f18b812aa2ec508dd10517128c88e8ee90be787b151` | 9 346 | 274 | `4f5598c8de2b48e6b7bc9c7a6c14e44c` |
| [UFh] | `processes/electromagnetic/standard/include/G4UniversalFluctuation.hh` | `1ff4cbe1b76ae9b8ede40ed9e9c28ba0805122d6890755fe92ee07029570094c` | 5 692 | 166 | `23279fcc45034ac8a2291d62365a35ea` |
| [II] | `processes/electromagnetic/standard/src/G4ionIonisation.cc` | `4b41d8e10b70ba9808d6fd6753b4685fc50a1222d73123d617d2c99f661e2f13` | 6 779 | 179 | `3cc2b75128094351920fb39fbe683a8b` |
| [ESU] | `processes/electromagnetic/standard/src/G4EmStandUtil.cc` | `c9bad251d39e1422c20789d8c0e10d92ca41f6183ac11d41c23cfdd821838d52` | 2 764 | 65 | `8d6fc8bdba2a47bdb25f3a87aad386d6` |
| [VELP] | `processes/electromagnetic/utils/src/G4VEnergyLossProcess.cc` | `2f53a58c5a0c0feb4cd3a5fef3c7464fe8220939f00cba6e6768ba2fa8f116e5` | 52 184 | 1511 | `56e7034106f54131b871943a5392efed` |
| [IEC] | `processes/electromagnetic/utils/src/G4ionEffectiveCharge.cc` | `891854a97b97f038f9a0fa88b9f2f4ed139ea080e8edf9bd349cfd93031ae4c5` | 7 225 | 188 | `1954c41a50c64064bba3e97d0c20b108` |
| [I90] | `materials/src/G4ICRU90StoppingData.cc` | `afd1479497a605d9675c1fe2a3870ddb3f305110fb2d76e13bb600fa3ab281cd` | 9 639 | 175 | `a45512277d5547d69244f16daad75401` |
| [LSD] | `processes/electromagnetic/standard/src/G4LindhardSorensenData.cc` | `70ae7c6ba3ef990b7005240bd2b60e4040a3fb60bb1c9e4dbf829c109180eda9` | 10 838 | 212 | `11afcec583d046bc8939ec79e746150c` |

These earlier registered files are also present in the cache (hash re-verified) and are cited by line:

| Tag | Path | sha256 | Bytes | Lines | Registered in |
|---|---|---|---|---|---|
| [OPT4] | `physics_lists/constructors/electromagnetic/src/G4EmStandardPhysics_option4.cc` | `04337300cc119a4fd893d515059cc75cef9f042849b183848643a4c4bf104c89` | 11 966 | 334 | `p2-geant4-ion-models.md` §0 |
| [EMB] | `physics_lists/constructors/electromagnetic/src/G4EmBuilder.cc` | `40782f5c5bc0e66f9b4b635a08d8b0c236b8576ad5b885afca04b679bb18d3d1` | 15 843 | 466 | `p2-geant4-ion-models.md` F0 |
| [BB] | `processes/electromagnetic/standard/src/G4BetheBlochModel.cc` | `5dc2a9002bff25f239a6402d3976a1f285f071abf4b41585b76f275e4539e01a` | 18 688 | 492 | `p2b-follow-ups.md` §0 |
| [LS] | `processes/electromagnetic/standard/src/G4LindhardSorensenIonModel.cc` | `e16ed6cd1d132183684ef47948b6f8fe36526a48663a6dc1878d5bb8330eb12c` | 16 490 | 447 | `p2-geant4-ion-models.md` §0 |
| [IPM] | `materials/src/G4IonisParamMat.cc` | `990c6f315a3b6951253eb251024ba439f8edf4e151994c33991510645807849d` | 21 112 | 622 | `p2-geant4-ion-models.md` F0 |
| [EMP] | `processes/electromagnetic/utils/src/G4EmParameters.cc` | `aaa143c1cb34b544f0d67fa6a7242d2b3fbd86ebfeb1bf4bc5c9fcca2973909e` | 43 246 | 1627 | `p10-single-scattering-availability.md` §7 |

## 2. Fluctuation path for alpha and GenericIon under opt4

### 2.1 Which fluctuation model (ESTABLISHED)
- [II] L136–138 install `G4EmStandUtil::ModelOfFluctuations(true)` when no fluctuation model is set.
- [ESU] L52–56: `fDummyFluctuation` gives the dummy model; otherwise `isIon` gives `new G4IonFluctuations()`.
- So [OPT4] L133 `SetFluctuationType(fUrbanFluctuation)` does not reach alpha or GenericIon. That setting acts only
  in the `isIon == false` branches ([ESU] L57–60).
- Alpha and GenericIon both use `G4ionIonisation`: alpha through [EMB] (`new G4ionIonisation()` at L133, L140, L159,
  L163; which of these lines is the alpha registration was not re-read, TO VERIFY), GenericIon at [OPT4] L316–322
  (P2 §1).
- The same fluctuation object serves both energy-loss models of the process ([II] L154, L165).
- Fluctuations are on: `lossFluctuation = true` ([EMP] L120).

### 2.2 Handover G4IonFluctuations → G4UniversalFluctuation (ESTABLISHED code; INFERRED arithmetic)
- [IF] L76 `parameter(10.0*CLHEP::MeV/CLHEP::proton_mass_c2)`.
- [IF] L117–119: `if(dp->GetKineticEnergy() > parameter*charge*particleMass) return uniFluct->SampleFluctuations(...)`.
  Here `charge` is the bare PDG charge ([IF] L94, L432) and `particleMass` is the PDG mass of the track's
  definition ([IF] L93, L431).
- So the threshold is T > 10 MeV · z · M/m_p, i.e. T/A > 10 · z · (M/A)/m_p MeV.
  - INFERRED (nuclear masses): ≈ 19.9 MeV/u (He), 59.6 MeV/u (C), 79.4 MeV/u (O).
- For ions, [VELP] L608–616 call `fluc->SetParticleAndCharge(track.GetParticleDefinition(), q2)` before each step,
  with q2 the model's effective-charge-square ratio. [IF] L426–436 store the definition, mass and bare charge and
  forward q2 to Universal (L436, [UF] L259–272).
- Below the handover, [IF] uses its own dispersion:
  - Bohr form, [IF] L182–183: `(tmax/beta2 - 0.5*tcut)*twopi_mc2_rcl2*length*n_e*effChargeSquare`.
  - Multiplied by `fac_cut = 1 + (fac-1)·2 m_e β²/(tmax(1-β²))` (L200–204). INFERRED: with tcut = tmax for a
    heavy ion, fac_cut ≈ fac.
  - `fac = s1·q_eff²/z² + s2` (L401).
  - s1 is the Geissel et al. NIM B 195 (2002) factor (L406–422). It is replaced by the Yang tabulation only for
    `beta2 < 3*theBohrBeta2*Z` (L226), where theBohrBeta2 = 50 keV/m_p (L77) and Z = Z_eff of the material (L189).
    INFERRED: that is a few hundred keV/u, far below every energy here.
  - s2 is the Q. Yang et al. NIM B61 (1991) charge-exchange term. For ions `factor = charge·(charge/Z)^(1/3)`
    (L372), and the "ions in solids" parameter set is used for a non-gas medium (L383–386).
  - The old heavy-ion correction is commented out (L192–197).
  - Big-fraction widening applies when meanLoss > 0.2·E (L131–138; `minFraction = 0.2`, [IFh] L111).
  - Sampling below the handover: Gaussian truncated to [0, 2·meanLoss] if meanLoss/σ ≥ 2 (L146–151), Gamma if
    > 0.1 (L154–157), uniform otherwise (L160–161).
- INFERRED magnitudes (from the preliminary planning estimate, not recomputed): s2 ≈ 1 % for C near 50 MeV/u and
  ≈ 0.1 % for He near 15 MeV/u; s1 − 1 < 1e-3. Below the handover the variance is therefore Bohr × (1 + ≲ 1 %).
  This changes σ_R by < 0.1 %.

### 2.3 Electron cut of 1 km: tcut = tmax (ESTABLISHED code; cut value from the bundles)
- [VELP] L868–872: `tmax = currentModel->MaxSecondaryKinEnergy(dynParticle); tcut = std::min(cut, tmax);`, then
  `SampleFluctuations(couple, dp, tcut, tmax, length, eloss)`.
- The TOPAS bundles set CutForElectron = 1 km. The bundle request text records an e- production threshold of
  9.55 GeV in G4_WATER, from the P4 verbose smoke. That is orders of magnitude above tmax (≤ ≈ 1.2 MeV at 430 MeV/u).
  Hence **tcut = tmax**.
- In [UF] the Bohr branch then reduces to the Bohr variance with q_eff² (next subsection). In [IF] the `− 0.5 tcut`
  term likewise gives (tmax/β² − tmax/2).

### 2.4 Universal: Bohr branch versus Urban "Glandz" branch (ESTABLISHED)
- Gaussian/Bohr branch only if all three conditions hold ([UF] L126–127):
  - `particleMass > electron_mass_c2`;
  - `meanLoss >= minNumberInteractionsBohr*tcut`, with `minNumberInteractionsBohr = 10` ([UFh] L112);
  - `tmax <= 2.*tcut`.
  With tcut = tmax this is **meanLoss ≥ 10 · tmax**.
- In the Bohr branch: σ² = `(tmax/beta2 - 0.5*tcut)*twopi_mc2_rcl2*length*chargeSquare*n_e` (L129–130), with
  chargeSquare = q2 from SetParticleAndCharge (L271). Gaussian truncated to [0, 2·meanLoss] if sn ≥ 2 (L134–140),
  Gamma otherwise (L143–147). INFERRED: this equals ionmc's Bohr variance (decision 0043, bare charge) wherever
  q_eff = z.
- Otherwise the Glandz sampler runs ([UF] L152–167, `SampleGlandz` L172–239). Inputs:
  - e0 = `GetEnergy0fluct()` = 10 eV ([IPM] L406), e1 = I (L158, L179);
  - rate = 0.56, fw = 4, a0 = 42, nmaxCont = 8 ([UFh] L114–117).
  - **Width correction for small cuts:** `scaling = min(1 + 0.5 keV/tcut, 1.5)`; the mean loss is divided by it
    and the sample multiplied by it (L162–163, L167). The variance therefore gains a factor `scaling`. INFERRED:
    ≤ 1.003 at ≥ 80 MeV/u with tcut = tmax.
  - Excitation, single level: a1 = meanLoss·(1 − rate)/e1 (L182), rescaled by fw (L183–190). Sampled as Poisson or
    Gauss ([UFh] L128–141).
  - Ionisation: a3 = rate·meanLoss·(tcut − e0)/(e0·tcut·ln(tcut/e0)) (L193–194). Energies follow 1/E² on
    [α·e0, tcut], Poisson in number (L221–233). The part below α·e0 is approximated by a Gaussian for a3 > 8
    (L212–218).
- **Variance of the Glandz branch (INFERRED algebra on these lines):**
  σ² ≈ scaling · [rate · meanLoss · (tmax − e0)/ln(tmax/e0) + (1 − rate) · meanLoss · fw · I].
  - The Gaussian sub-part below α·e0 uses the fixed-number variance n·Var(E), not the compound-Poisson n·⟨E²⟩
    (L217). INFERRED: the deficit is below 1e-3 of the total, because the variance is carried by the few
    collisions near tmax.
  - Mean preserved: rate·meanLoss + (1 − rate)·meanLoss.
  - The distribution is skewed toward high losses: of order 10 Poisson collisions per 0.1 mm step, with 1/E²
    energies up to tmax.

### 2.5 Glandz versus Bohr per unit length (INFERRED planning arithmetic, `python3 -I`, 2026-11-13)
- Method: plain Bethe stopping (I = 78 eV, Z/A = 0.5551, no shell or density correction), nuclear masses,
  tcut = tmax.
- Ratio: R_var = [0.56·S·tmax/ln(tmax/e0) + 0.44·4·I·S] / Bohr. It depends on energy per nucleon only, not on z,
  since q² cancels between S and the Bohr term.
- Minimum step for the Bohr branch: 10·tmax/S. It scales as 1/z².

| E/A [MeV/u] | 80 | 100 | 150 | 200 | 250 | 290 | 430 |
|---|---|---|---|---|---|---|---|
| R_var (Glandz/Bohr) | 0.952 | 0.969 | 1.009 | 1.047 | 1.082 | 1.108 | 1.190 |
| min step for Bohr branch, He [mm] | 0.54 | 0.80 | 1.65 | 2.73 | 4.01 | 5.16 | 9.97 |
| same, C [mm] | 0.06 | 0.09 | 0.18 | 0.30 | 0.45 | 0.57 | 1.11 |
| same, O [mm] | 0.03 | 0.05 | 0.10 | 0.17 | 0.25 | 0.32 | 0.62 |

### 2.6 Run-time corroboration (configuration printout only; ESTABLISHED for the archived run)
- Source: the Geant4 EM table printed in the stdout of an archived V3-008 TOPAS run (OpenTOPAS 4.3, opt4 EM-only).
  Run REF-b73b308392747b407aba-884b2a64, stdout sha256 `04b4bd6580cc27e88b732b095be4a876d476df750464b7f83cc5263b7d9a303f`,
  lines 94–103 and 284–306. Only the configuration block was read.
- The same `linLossLim= 0.02` lines occur in the execution-only timing smokes of revision 2 (grep only).
- What the printout shows:
  - `Use ICRU90 data 1`; "Lowest muon/hadron kinetic energy 1 keV"; "Type of fluctuation model for leptons and
    hadrons Urban".
  - ionIoni for GenericIon: `StepFunction=(0.1, 0.001 mm) … fluct: 1, linLossLim= 0.02`, model `LindhardSorensen`
    over Emin 0 eV – Emax 10 GeV.
  - ionIoni for alpha: `StepFunction=(0.1, 0.02 mm) … linLossLim= 0.02`, models `BraggIon` below and `BetheBloch`
    above 7.9452 MeV.
  - dE/dx and range tables from 100 eV to 10 GeV in 160 bins.
- This closes P2's INFERRED item "LightIons applies to alpha" (the step function of [OPT4] L123–124).
- The printout does not name the fluctuation model. The source path above names it.

## 3. Stopping path (alpha and GenericIon above ≈ 2 MeV/u)

### 3.1 Alpha (ESTABLISHED)
- [II] L141–142: `if (pdg == 1000020040) SetEmModel(new G4BraggIonModel())` below eth.
- [II] L129: eth = 2 MeV · M/m_p. The constructor default at L83 is overwritten. `G4BetheBlochModel` is used above
  eth ([II] L157–165).
- [BB] L273: `if(nullptr != fICRU90 && kineticEnergy < fProtonTlimit)`, with fProtonTlimit = 10 GeV (L79).
- [BB] L284–285: for alpha with `kineticEnergy <= fAlphaTlimit` (1 GeV, L78) the dE/dx is
  `GetElectronicDEDXforAlpha(...)`, the **ICRU 90 alpha array** (P2b §1: branch L273–302 and the direct return).
- [BB] L343: `CorrectionsAlongStep`: `if(isAlpha) { return; }`. No along-step corrections for alpha.
- [I90] L82: ICRU 90 materials G4_AIR, G4_WATER, G4_GRAPHITE.
- [I90] L164–169: vectors built as `G4PhysicsFreeVector(n, e[0], e[n-1], true)`. That the last argument is the
  spline flag is INFERRED (needs `G4PhysicsFreeVector.hh`).
- He 250 MeV/u is ≈ 1.0 GeV, at the 1 GeV alpha limit. INFERRED: relevant only for the first millimetres of that
  one case.

### 3.2 GenericIon: C and O (ESTABLISHED, except where marked)
- **Table.** `G4LindhardSorensenIonModel` over the whole range (runtime printout §2.6; [OPT4] L316–322). Above
  escaled = 2 MeV ([LS] L82) the table is `fBBModel->ComputeDEDXPerVolume` (P2 §3, [LS] L240–242).
  - For GenericIon this is the **ICRU 90 proton array at the proton-scaled energy × chargeSquareRatio**
    ([BB] L273 and L290–293 per P2b §1; the ion kinetic energy stays < 10 GeV in every case).
  - C/O tracks reuse the GenericIon tables ([II] L117–126: base particle GenericIon) with
    massRatio = M_base/M ([VELP] L256–261, L544–548).
- **Effective charge** ([IEC] L66, L108, L113): bare charge when T·m_p/M > z · 20 MeV, i.e. above ≈ 40 (He),
  119 (C) and 159 (O) MeV/u (INFERRED arithmetic). Below that, a Ziegler-type q_eff applies.
  - [VELP] L608–612 update q² every step for `isIon` particles.
- **Along-step LS terms, not in the table.** [VELP] L853–861 call `CorrectionsAlongStep` for `isIon` after the
  table-based loss. [LS] L296–300 then compute:
  - `deltaL0 = 2·BarkasCorrection(p, mat, e)·(charge−1)/charge`;
  - `deltaL = lsdata->GetDeltaL(Zin, gam)`;
  - `res = eloss + twopi_mc2_rcl2·q2·n_e·(deltaL + deltaL0)·length/beta2`.
  - `e` is the mid-step energy max(E − eloss/2, E/2) ([LS] L263).
  - [LSD] L71 `GetDeltaL(Z, gamma)` interpolates a tabulated (Z, γ) set (InitialiseData L115 ff.; values not
    transcribed).
- **Linear-loss limit and pre-step dE/dx.**
  - [II] L80 `SetLinearLossLimit(0.02)`, which sets `actLinLossLimit` ([VELP] L1427–1433), so the default 0.01
    ([EMP] L154; [VELP] L238) is not applied.
  - [VELP] L813–814: `eloss = length * GetDEDXForScaledEnergy(preStepScaledEnergy, …)`, i.e. dE/dx at the
    **pre-step** energy.
  - Only if eloss > 0.02 · E_pre ([VELP] L824–828) is the loss recomputed exactly from the range table:
    E_pre − E(R − length).
  - The LS terms use the mid-step energy (above), so the pre-step bias concerns the table part only.

### 3.3 End of track (ESTABLISHED)
- Step limit [VELP] L567–572: x = R·dRoverRange + finR·(1 − dRoverRange)·(2 − finR/R) if R > finR, else x = R.
  - Step functions: (0.1, 1 µm) for ions, (0.1, 20 µm) for alpha ([OPT4] L123–124; runtime §2.6).
  - `Use cut as a final range` = 0 (printout).
- [VELP] L797–807: `if (length >= fRange || preStepKinEnergy <= lowestKinEnergy)` the whole pre-step energy is
  deposited along that step and the energy is set to 0.
- [VELP] L263–265: lowestKinEnergy = LowestMuHadEnergy for M ≥ 1 MeV. This is 1 keV ([EMP] L151; printout).
  [OPT4] L118 changes only the electron value.
- [VELP] L908–911: a post-step energy ≤ 1 keV is absorbed into the step.
- INFERRED: the residual is deposited continuously within the final ≤ 1 µm (C/O) or ≤ 20 µm (alpha), with no
  point deposit. ionmc's `continuous_residual_v1` (decision 0043 §4 amendment) deposits along its own table residual
  range below E_cut = 2 MeV/u. This is the same continuous-to-rest convention to within the table difference.
- **Nuclear stopping** is registered for alpha and GenericIon, `ICRU49NucStopping` up to 1 MeV (printout §2.6;
  [OPT4] L168).

## 4. Implications for the two attribution items

### 4.1 Item 1(a): σ_R, W and PE excess of the reference
- ESTABLISHED (§2.2–2.4): above 20/60/79 MeV/u (He/C/O) the reference samples Urban's Glandz model whenever the
  step's mean loss is below 10 · tmax, and Bohr only above it.
- ESTABLISHED by an earlier record: the step regime of the IDD runs. Plan line 533 records the step-length smoke
  result: the scorer bins limit the primary step in the IDD runs (97.5 %, 0.1 mm).
- INFERRED (§2.5): at 0.1 mm steps the Glandz branch applies:
  - He: above ≈ 50–80 MeV/u;
  - C: above ≈ 100–110 MeV/u;
  - O: above ≈ 150 MeV/u.
  There the variance per unit length relative to Bohr is 0.95 (80 MeV/u) to 1.19 (430 MeV/u).
- **Sign and ordering (INFERRED):**
  - The reference's σ_R exceeds Bohr when the slowing-down spends its straggling-weighted path at E/A ≳ 150 MeV/u.
    The weight is 1/S³, so it is concentrated near the end of range, where the ratio is ≤ 1 or the Bohr branch
    applies.
  - The excess grows with initial energy per nucleon, not with z, and shrinks with longer steps (more steps fall
    in the Bohr branch).
  - The shape is skewed (§2.4), which widens a distal-edge-fitted width beyond its variance.
- **Magnitudes (INFERRED planning estimates; preliminary scratch emulation 2026-11-11, not committed):**
  - Method: integral σ_R² = ∫(dΩ²/dx)/S³ dE, plain Bethe S, per-energy branch choice at 0.1 mm steps, and Bohr × 1.01
    below the handover.
  - Result: σ_R(G4)/σ_R(Bohr) ≈ 0.99 (He 150), 1.02 (He 250), 1.03 (C/O 290), 1.04 (C/O 350), 1.06 (C/O 430).
  - Against the item 1(a) residuals (He ≈ 0–3 %, C/O ≈ 4–9 %) this has the right sign and ordering. It accounts for
    most of the size at C/O 290–430. The remaining gap (to the ≈ 4–9 % and the O 100 PE ≈ 15 %) is unexplained
    by this mechanism and is TO VERIFY. The candidates are the skewness and the ICRU 90/LS-based S in place of
    plain Bethe.
- **Held-out energies (INFERRED):** He 100, C 150, O 200 lie at R_var ≈ 0.97 / 1.01 / 1.05 at entry, and below
  ≈ 1 along most of the weighted path. Expected σ_R(G4)/σ_R(Bohr) ≈ 0.98–1.02.
  - The held-out results note records W(ionmc)/W(TOPAS) − 1 = +3.6 % / −0.4 % / −2.5 % and PE within +0.1 … −1.0 %
    (TD-EM). This is consistent with a near-Bohr reference at these energies; it is no confirmation of the
    mechanism.
- **ionmc side:** no change. ionmc keeps the declared Bohr model with the bare charge (decision 0043). The reference
  excess is a step-dependent property of the reference's sampler, not physics absent from ionmc. Adopting it would
  be fitting to the reference.

### 4.2 Item 1(b): TOPAS R80 longer by +0.19 … +0.29 % (main set) and +0.15 … +0.21 % (held-out)

Sign convention: "+" means the mechanism makes the TOPAS range longer.

| # | Mechanism | Evidence | Sign | Size (INFERRED unless stated) |
|---|---|---|---|---|
| 1 | Pre-step dE/dx for steps losing ≤ 2 % ([VELP] L813–828, [II] L80) | code ESTABLISHED; 0.1 mm IDD step regime ESTABLISHED (plan L533) | + | ≈ +0.03 … +0.06 % at 0.1 mm (scratch emulation 2026-11-11); +0.14 … +0.23 % would need ≈ 1 mm steps, which the IDD runs do not have |
| 2 | Base stopping array: ICRU 90 alpha (He), ICRU 90 proton × q² at equal velocity (C/O), versus ionmc's own Bethe-with-corrections table | code ESTABLISHED (§3.1–3.2); difference size TO VERIFY (D-3) | TO VERIFY | must carry ≈ +0.13 … +0.16 % for He (no LS terms for alpha, [BB] L343), against decision 0038's ≤ 0.07 % |
| 3 | LS ΔL + Barkas 2B(z−1)/z, along the step, C/O only ([LS] L296–300) | code ESTABLISHED; sign TO VERIFY | INFERRED − | Physics (Lindhard–Sørensen 1996): at β ≈ 0.65–0.73 the Mott-type part (∝ zαβ) dominates the negative Bloch part (∝ (zα/β)²), so ΔL > 0; the Barkas term is positive for z > 0. Both raise S, i.e. the C/O range is shorter by ≈ 0.1–0.2 %. The ICRU 90 proton base already contains the z = 1 corrections, so ΔL is the z-dependent excess. **Not ESTABLISHED:** the tabulated [LSD] values were not read for sign |
| 4 | q_eff < z below 40/119/159 MeV/u ([IEC] L113) | switch ESTABLISHED | + | < 0.01 % |
| 5 | Nuclear stopping ≤ 1 MeV | registration ESTABLISHED (printout) | − | < 0.01 % |
| 6 | Table interpolation: 160 bins, 100 eV–10 GeV; ICRU 90 vectors with spline flag ([I90] L169) | INFERRED | ? | ~1e-4 |
| 7 | Mean-loss bias of the fluctuation samplers: truncations symmetric about the mean ([UF] L140, [UFh] L157; [IF] L151) | code ESTABLISHED | ≈ 0 | second order in σ, ≲ 1e-4 |
| 8 | Detour (projected < path), R80 versus CSDA | INFERRED | − / ≈ 0 | ≤ 0.05 %, cannot lengthen |

- **Net (INFERRED):** at the established 0.1 mm step regime, mechanism 1 explains only ≈ 0.03–0.06 % of the
  +0.15 … +0.29 %. The He offset (+0.19 % main, +0.21 % held-out) has no LS term and must sit mostly in mechanism
  2 (alpha array versus ionmc's table) or in R80 extraction.
- For C/O, mechanism 3 (if negative) makes mechanism 2 for the proton array still larger. This needs ≈ +0.3–0.4 %
  for C/O.
- This revises the preliminary note, which named mechanism 1 as the leading candidate on the assumption of 0.5–1 mm
  steps.
- The ICRU 90 base arrays are therefore the leading candidate. D-3 (integration of the [I90] arrays to CSDA ranges at
  the case energies, with and without the LS terms) is the discriminating calculation.

## 5. What remains TO VERIFY, and what resolves it

1. **Sign and size of LS ΔL at Z = 6, 8 and γ ≈ 1.15–1.46.**
   - Read `G4LindhardSorensenData::InitialiseData` ([LSD] L115 ff.) for the (Z, γ) grid and the sign. Report only
     the sign and the ratio to L; do not transcribe values.
   - Read `G4EmCorrections::BarkasCorrection` (not in the pinned set; `G4EmCorrections.cc`, scratch hash
     `a357d076…` in `p9-freeze-record.md` §3) for the Barkas sign and normalisation.
2. **Size of mechanism 2.** D-3: integrate the ICRU 90 alpha and proton (× q², at equal velocity) arrays of [I90]
   for G4_WATER to CSDA ranges at the case energies, and compare with ionmc's table ranges. Also run D-5, the
   optional scorer-free TOPAS dump of the Geant4 dE/dx and range tables for alpha and GenericIon(6,12)/(8,16) in
   G4_WATER, if TOPAS exposes it. The dump would directly show the table range, including the base-particle
   scaling but excluding the along-step LS terms.
3. **Mechanism 1 with measured steps.** Re-run the pre-step emulation with the step-length smoke's measured step
   distribution and ICRU 90 S in place of plain Bethe.
4. **σ_R magnitude.** Recompute σ_R(G4)/σ_R(Bohr) with S from the ICRU 90 arrays (D-3 instrument), the measured
   step distribution and the Glandz skewness (sampled, not variance-only). Compare the predicted W/PE shift with
   the committed D-2 forward model. This decides whether the C/O 4–9 % and the O 100 PE ≈ 15 % residuals are fully
   attributed.
5. **isIon for alpha.** `G4EmTableUtil::CheckIon` (not acquired) decides whether alpha takes the `isIon` paths
   ([VELP] L544–549, L608–616, L854, L912). If it does not, the fluctuation model keeps alpha's PDG charge from
   `InitialiseMe` ([IF] L90–97). This leaves the He threshold of §2.2 unchanged, as long as `InitialiseMe` runs
   before the first step; where it is called was not read (`G4VEmModel`/`G4EmModelManager`).
6. **Alpha registration line in [EMB].** Identify which of L133/L140/L159/L163 is alpha. The runtime printout
   already establishes the alpha model set.
7. **`GetDEDXForScaledEnergy` factor.** That it multiplies by fFactor = q²·(density factor) is INFERRED;
   `G4VEnergyLossProcess.hh` was not acquired.

**None of this changes an ionmc model, parameter or grade.** Plan D, STOP 3: closing the reference's range offset
or adopting its fluctuation sampler would be fitting to the reference. Both remain attributions recorded as
reference-side properties.

## 6. Differences from the preliminary note of 2026-11-11

- **Step regime.** The preliminary left the TOPAS step length open and treated 0.5–1 mm steps as plausible. The
  plan's step-length smoke result (plan L533: 0.1 mm, 97.5 %) settles the IDD regime as 0.1 mm. Consequence:
  mechanism 1 drops from "leading candidate" to ≈ +0.03–0.06 %, and the base ICRU 90 arrays (D-3) become the
  leading candidate for the range offset. The plan text at item 1(b) ("up to ≈ +0.2 % at 1 mm steps") is
  superseded in emphasis, not in fact.
- **Glandz width correction.** The preliminary variance formula omitted the small-cut width correction ([UF]
  L162–163, L167, factor ≤ 1.003 here) and the fixed-number Gaussian sub-part (L212–218, < 1e-3). Both are
  negligible; the formula now names them.
- **Bohr branch condition.** It has three conditions ([UF] L126–127), not one: particleMass > m_e,
  meanLoss ≥ 10·tcut and tmax ≤ 2·tcut. Bohr uses Gamma sampling for sn < 2 (L143–147).
- **Handover bookkeeping.** `charge` and `particleMass` at [IF] L117 come from `SetParticleAndCharge` (L426–436),
  called pre-step by [VELP] L608–616 for ions; the preliminary cited only L117–118. For alpha this depends on
  `isIon` (TO VERIFY item 5).
- **[II] eth.** The constructor sets 2 MeV (L83); L129 overwrites it with 2 MeV·M/m_p. The runtime printout
  confirms 7.9452 MeV for alpha.
- **Line corrections:**
  - [LS] along-step lines are L296–300 (preliminary and P2: L292–300 block; formula unchanged).
  - [BB] alpha branch L284–285 (preliminary L273–296 range).
  - [I90] vector construction L164–169 (preliminary L156–173).
  - [IF] Yang ion branch L372 with the solid set L383–386.
  - [OPT4] SetStepFunction lines L121–124 (ions L124, light ions L123).
  - All other preliminary line numbers ([IF] L76, L117, L182, L200, L226, L401, L412; [UF] L127, L129, L182,
    L194; [UFh] L112, L114–117; [II] L80, L129, L137, L141–142, L158; [ESU] L53–56; [VELP] L263, L571, L797,
    L813, L824, L854, L869, L909–911; [IEC] L66, L113; [IPM] L406; [BB] L78, L79, L343; [LS] L82, L263;
    [I90] L82; [EMP] L151, L154) were re-verified in the pinned files and stand.
- **New runtime corroboration.** The archived EM printout (§2.6) establishes the step functions per particle,
  linLossLim 0.02, the model sets and the table binning. The preliminary had these as source-only or INFERRED.
- **New held-out estimates.** R_var at 100 MeV/u (0.969) and 200 MeV/u (1.047), and the held-out comparison of
  §4.1.
- **Unchanged and still not established:** the LS ΔL + Barkas sign, still TO VERIFY. The preliminary called it
  "most likely raises the stopping"; this record keeps it INFERRED from theory, not from [LSD].
