# P2 follow-up, part 2: disposition of freeze condition FC-2 (V3-008, record of 2026-11-01)

Status: pre-freeze source reading and a read of the archived smoke logs (no new run), written 2026-11-01 by the
planner at task head `6bdf3897` (freeze CANDIDATE revision 1). It disposes freeze condition FC-2 of
`validation/plans/v3-008/p9-freeze-record.md` §0. The file name is the one prescribed by the orchestrator; it is
**not** the P2b check of the plan (the Urban-emulator applicability check, CLOSED as inapplicable). No transport
output, no scored TOPAS output and no tabulated stopping, range or scattering reference value appears here. Geant4
source constants are quoted with the line number (Geant4 Software Licence; short excerpts with attribution; nothing
is redistributed as data). Untrusted sources were read with `cat`/`grep`/`sed` only, never compiled or executed.

Markers as in `p2-geant4-ion-models.md` §0: **FACT** (read in a pinned file or archived log, with line), **INFERRED**
(algebra, or code outside the pinned set; the settling file is named), **OPEN**.

## 0. Sources

| Tag | File (Geant4 v11.4.2, upstream path under `source/`) | sha256 | Lines | Status |
|---|---|---|---|---|
| [BB] | `processes/electromagnetic/standard/src/G4BetheBlochModel.cc` | `5dc2a9002bff25f239a6402d3976a1f285f071abf4b41585b76f275e4539e01a` | 492 | acquired, exploratory, 2026-11-01 (same hash as the `p3b-ls-default.md` §2 scratch download) |
| [BB-h] | `processes/electromagnetic/standard/include/G4BetheBlochModel.hh` | `cf05e98c6b56aa07e227800ccff169b74e2c724d56ce5753fcab886be8fac58e` | 169 | acquired, exploratory, 2026-11-01 |
| [NMB] | `materials/src/G4NistMaterialBuilder.cc` | `a1f8944dd63f4c9915c2205ee2ca5ec0ffade15ac0ca8afe38d53597c2671abc` | 2002 | acquired by the orchestrator 2026-11-01 (lineage); identified by its header (L30 `File name: G4NistMaterialBuilder`) |
| [LS] [URB] [ICRU73] [OPT4] [EMB] [HMSC] [VMSC] [IPM] [EMP] | as in `p9-freeze-record.md` §3 | e16ed6cd… 22c5f52b… cf771a5a… 04337300… 40782f5c… a6643869… 87b9f203… 990c6f31… aaa143c1… | – | pinned earlier |
| [EMC], [IEC] | `G4EmCorrections.cc`, `G4ionEffectiveCharge.cc` (session scratch, `p3b-ls-default.md` §2) | `a357d076…`, `891854a9…` (full values in `p9-freeze-record.md` §3) | – | scratch, lineage |

The three new hashes were recomputed here with `sha256sum` and match the cache file names. Archived smoke logs read
(stdout only; execution-only smokes, no scorer output exists): carbon `REF-5910656de47adfc817ef-e043426a`, oxygen
`REF-8031592a201b7ddf5a48-7a2c2ae5`, carbon SS `REF-61e64b1b80f37f99f9e5-02b74b6c`, alpha verbose
`REF-c98e198face1e783fde0-4f13dcfb`. These are the full run ids as materialized in `.ionmc-cache/reference-runs/`.

## 1. G4BetheBlochModel: the ICRU 90 branch, the corrections path, the effective charge, the low-energy limit

**ICRU 90 branch (FACT [BB] L104–107, L273–302).**
- L104–107: `if(G4EmParameters::Instance()->UseICRU90Data() && (pname == "proton" || pname == "GenericIon" || pname == "alpha")) { fICRU90 = nist->GetICRU90StoppingData(); }`.
- [OPT4] L131 `SetUseICRU90Data(true)`; the opt4 smokes print `Use ICRU90 data 1` (carbon/oxygen stdout L89).
- L273: `if(nullptr != fICRU90 && kineticEnergy < fProtonTlimit)`, with `fProtonTlimit(10*CLHEP::GeV)` (L79) and
  `fAlphaTlimit(1*CLHEP::GeV)` (L78). L278: `iICRU90 = fICRU90->GetIndex(baseMaterial)`.
- If the material is indexed (L280):
  - alpha (L283–289): `kineticEnergy <= fAlphaTlimit` → `GetElectronicDEDXforAlpha(iICRU90, kineticEnergy)`;
    above 1 GeV → `GetElectronicDEDXforProton(iICRU90, e)*chargeSquareRatio` at `e = T·m_p/M`;
  - GenericIon and protons (L290–293): `GetElectronicDEDXforProton(iICRU90, e)*chargeSquareRatio` at the
    proton-scaled energy;
  - L294 `dedx *= material->GetDensity()`; L295–298 the restricted-loss subtraction only if `cutEnergy < tmax`;
  - L300 `if(dedx > 0.0) { return dedx; }`: the general Bethe formula (L304–320) and the high-order corrections
    (L323) are then **not** evaluated.
- Which materials `GetIndex` accepts, and the data span of the proton and alpha arrays, are in
  `G4ICRU90StoppingData.cc` (not pinned). INFERRED: G4_WATER is indexed (the same three ICRU 90 materials,
  `{"G4_AIR","G4_GRAPHITE","G4_WATER"}`, are listed in [ICRU73] L88), and `fProtonTlimit`/`fAlphaTlimit` mark the
  array ends. If water were not indexed, or an array returned 0, the Bethe path below would apply.

**Consequences for the V3-008 cases in G4_WATER under option 4 (FACT from the lines above, conditional on the
INFERRED indexing).**
- **alpha** (He cases, 50–250 MeV/u = 200–1000 MeV ≤ `fAlphaTlimit`): above the BraggIon/BetheBloch switch the
  stopping is the ICRU 90 **alpha** array, multiplied by density only. `chargeSquareRatio` enters only the
  restricted subtraction, which is inactive with the 1 km electron cut (e⁻ production threshold 9.55 GeV in the
  verbose smokes, P4). `CorrectionsAlongStep` returns immediately for alpha (L342–343). The decision 0043 §6 claim
  "ICRU 90 alpha array (He) identical to TOPAS alpha stopping" is **CONFIRMED at the source level** (from "research
  note"), with the one INFERRED premise above.
- **GenericIon** (C, O): the LS model owns a `G4BetheBlochModel` ([LS] L81) and initialises it with the GenericIon
  definition ([LS] L102), so L104–107 set `fICRU90`. Above 2.0 MeV scaled ([LS] L240–242) the base table is the
  ICRU 90 **proton** array at equal velocity × q² (L291–292). The LS `CorrectionsAlongStep` then adds
  q²·(ΔL_LS(Z) + 2·Barkas·(z−1)/z) ([LS] L294–300). The decision's "ICRU 90 proton array at equal velocity: research
  note, partly confirmed" is **CONFIRMED at the source level** (same premise).
- **Resolution of `p3b-ls-default.md` §5 ("not settled whether the G4BetheBlochModel base table for GenericIon
  already contains the charge-1 Bloch/Mott terms (L323) or the ICRU 90 proton array")**: for G4_WATER under option 4
  it is the ICRU 90 proton array (L291–292, returning at L300); L323 is not reached. INFERRED: an evaluated proton
  stopping table implicitly contains the z = 1 Barkas and Bloch contributions. The LS path adds the Barkas top-up
  (z−1)/z, which is consistent with that, and the full ΔL_LS(Z), whose z = 1-equivalent part may be counted twice.
  The size of this possible double count is unquantified; it is a lineage item of the reference (I3c attribution,
  V1 attribution budget), not a freeze condition.
- **Single-scattering constructor:** the SS smoke prints `Use ICRU90 data 0` (carbon SS stdout L89). There
  `fICRU90` stays null. The GenericIon base table above 2 MeV scaled is then the general Bethe formula with the
  L323 high-order corrections at q, and the [ICRU73] loader builds G4_WATER from element files (P2 F3: the
  named-material list is inoperative). This affects only V5-SS lineage (an angular/lateral ratio statistic), not a
  stopping row.

**Corrections path (FACT; confirms the P3b reading).** [BB] L323 `dedx += corr->HighOrderCorrections(p,material,kineticEnergy,cutEnergy);`
→ [EMC] `HighOrderCorrections`: `sum = 2.0*(Barkas + Bloch) + Mott`, scaled by `q2` ([EMC] L172–198). `q2` comes from
`SetupKinematics` ([EMC] L158–160: `if(charge > 1.5) { charge = effCharge.EffectiveCharge(p,mat,kinEnergy); } q2 = charge*charge;`).
This path is reached only when the ICRU 90 branch does not return (see above).

**Effective charge for alpha (FACT).** [BB] L129–136: `chargeSquareRatio = corr->EffectiveChargeSquareRatio(p, mat, kinEnergy)`,
which delegates to `G4ionEffectiveCharge` ([EMC] L826, L846). [IEC] L104–114: the ZBL 1985 effective charge, with the
bare charge returned when `reducedEnergy > effCharge*energyHighLimit` (`energyHighLimit = 20.0*CLHEP::MeV`, L66),
i.e. above 40 MeV proton-scaled for helium. [IEC] L119–180 is the helium branch (`Zi == 2`, ZBL polynomial,
L122–123 coefficients). For alpha in water under option 4, the effective charge multiplies the stopping only
above 1 GeV (outside the He domain), so it is inert for the He V1 cases (INFERRED from the alpha branch above).

**Low-energy limit (FACT printout; INFERRED code).** [BB] L84 `SetLowEnergyLimit(2.0*CLHEP::MeV)`. The alpha
printout reads `BraggIon : Emin= 0 eV Emax=7.9452 MeV` and `BetheBloch : Emin=7.9452 MeV Emax= 10 GeV` (carbon stdout
L300–301, identical in every opt4 smoke). 2 MeV × M_α/M_p = 2 × 3727.379/938.272 = 7.9452 MeV (scratch arithmetic). The
model switch is therefore 2 MeV at proton-equivalent velocity (≈ 2 MeV/u). INFERRED: the scaling is applied by
`G4ionIonisation` (not pinned). Below the switch, `G4BraggIonModel` (not pinned) applies over the last ≈ 2 MeV/u of
the He track. GenericIon prints a single `LindhardSorensen` model over 0–10 GeV (carbon stdout L284); its internal
Bragg/BetheBloch switch at 2.0 MeV scaled is P2 FACT ([LS] L82, L240).

**Disposition.** The pin is done ([BB], [BB-h] acquired). The ICRU 90 branch is settled for G4_WATER + option 4,
with one INFERRED premise (the `G4ICRU90StoppingData` material index and span). `p3b` §5 item 5 is resolved as above;
a possible double count of the z = 1 part of ΔL_LS is recorded as a lineage item. No row definition or tolerance
changes. He: V1 tests transport mechanics and straggling only (as already stated in decision 0043 §6).

## 2. G4_WATER material: compiled definition, radiation length, and the record gap of P4

**Definition (FACT [NMB]).**
- L712: `AddMaterial("G4_WATER", 1.0, 0, 78., 2, kStateLiquid);` gives density 1.0 g/cm³, a mean excitation energy of
  78 eV (source constant), two components and the liquid state.
- L713–714: `AddElementByAtomCount("H", 2); AddElementByAtomCount("O", 1);`.
- L715: `chFormulas[nMaterials-1] = "H_2O";`.
- `BuildMaterial` (L221–236) first sets I from the chemical formula (`FindMeanExcitationEnergy`, i.e. the [IPM]
  "H_2O" → 75 eV entry of P2 F4). L230–233 (`// If exists, NIST DB data always overwrites other data`;
  `if (ionPotentials[i] > 0.0) { exc1 = ionPotentials[i]; }`) then overwrite it, and L234–235 call
  `SetMeanExcitationEnergy(exc1)`.
- **Hence G4_WATER carries I = 78 eV.** This settles P2 F4 ("stays INFERRED until … G4NistMaterialBuilder.cc shows
  78 eV"): FACT at the source level. It matches ionmc's construction value (decision 0043 §6, "I_water 78 eV",
  shared lineage).

**Radiation length (derived; X0 is not a compiled constant).** Geant4 computes X0 at material construction from the
composition: `G4Material::ComputeRadiationLength` sums, over elements, the atom number density times the per-atom
Tsai factor of `G4Element` (`ComputeLradTsaiFactor`, `ComputeCoulombFactor`). Those two files are not pinned, so the
Geant4 code form is INFERRED. The physics form is Tsai's (PDG RPP):

  1/X0 = (N_A/M_H2O)·[2·f(1) + 1·f(8)], f(Z) = 4 α r_e² Z [Z (L_rad(Z) − f_c(Z)) + L′_rad(Z)],

- f_c(Z) = a²[1/(1+a²) + 0.20206 − 0.0369 a² + 0.0083 a⁴ − 0.002 a⁶], with a = αZ;
- L_rad, L′_rad = 5.31, 6.144 for Z = 1 (Tsai's light-element values) and ln(184.15 Z^−1/3), ln(1194 Z^−2/3) for Z = 8;
- composition H₂O by atom count ([NMB] L713–714); density 1.0 g/cm³ ([NMB] L712).

CALC (scratch script `v3-008/fc2/x0_water.py`, sha256 in §6; numpy; standard atomic weights H 1.00794, O 15.9994,
which are not read from Geant4's element builder): **X0(G4_WATER) = 36.083 g/cm² = 36.083 cm** at density 1.0. It
agrees with the INFERRED value of P2 §6.1 / F4.

**Comparison with ionmc.** ionmc's materials module (`src/ionmc/materials.py`) holds no radiation length. The only X0
in `src` is `radiation_length_g_cm2(WATER)` in `src/ionmc/transport/mcs_checks.py` L153–173: the same Tsai/PDG form,
used by the generalised-Highland cross-check; the transport MCS (exact-Ω mapping of f_dM) does not use X0. **Ratio
ionmc / Geant4-form X0 = 1.0000000** (identical form and inputs). Sensitivity to the element weights: +0.001 g/mol on
H changes X0 by 1.1e-4 relative, so a difference in Geant4's element-builder masses would be negligible.

Where X0 matters: only the Urban θ0 of the TOPAS reference (y = t/X0, P2 §2.1), i.e. the V2 Case C lineage. The
withdrawn Urban emulator would also have needed it.

**Correction of the P4 record (dated note, 2026-11-01).** The P4 row (r5 sentence) says that the scorer-free smokes
archive "the G4_WATER material table (RadL)". That is not the case.
- The archived stdout files contain the EM parameter block, the per-process model tables and, in the verbose smokes
  only, the production-cut couple table (`Material : G4_WATER` with range cuts and energy thresholds; alpha verbose
  stdout L535–554).
- They contain **no material table**: no density, I-value or radiation length. A grep for "RadL", "Rad. len",
  "ExcitationEnergy" and "I mean" in all nine materialized smoke logs returns 0 hits.
- "Archived" is therefore corrected to: **source definition ([NMB] L712–715, L221–236) + derived X0 (this section)**.
- Optional later smoke (not required, not a freeze condition): OpenTOPAS documents `i:Ma/Verbosity = 0 # Set to 1 to
  report each time a material is defined` (OpenTOPAS defaults page,
  https://opentopas.readthedocs.io/en/main/parameters/defaults.html, read 2026-11-01). Whether that report prints the
  Geant4 material table (density, I, RadL) is not established. No documented OpenTOPAS parameter that prints the
  Geant4 material table was found. A Geant4 UI command passed through TOPAS (e.g. `/material/g4/printMaterial`) is not
  documented on the pages read and is untested.

**Disposition.** I-value: FACT (78 eV, source). X0: derived from the source composition, with ratio 1 to ionmc's
cross-check value. The P4 record gap is closed by correcting the record, not by a new run.

## 3. λ charge scaling for ions and the `G4VMscModel` initialisers

**What the P2 row and F1/F2 ask.**
- "λ scaling and lateral displacement for ions (`G4VMscModel`, `G4EmTableUtil`)": whether a z² enters the Urban
  transport mean free path λ for C/O (GenericIon tables), and the lateral-displacement setting.
- F2: the in-class initialisers of `facrange`, `latDisplasment`, `skin` and `dtrl` in `G4VMscModel.hh`, and the final
  step-limit type and RangeFactor (the competing "nucleus" override path).

**Final msc parameters for GenericIon and alpha (FACT, run time and source).**
- Printout (carbon stdout L274–277 and L290–293; identical in the oxygen and He smokes):
  `UrbanMsc : Emin= 0 eV Emax= 10 GeV` / `StepLim=Minimal Rfact=0.2 Gfact=2.5 Sfact=0.6 DispFlag:1 Skin=3 Llim=1 mm`.
  This is the per-model `DumpParameters` line ([VMSC] L124–126), printed after initialisation, so it gives the values
  in force. It settles F1's OPEN "competing path": no override to other values reaches the model.
- The source of each value:
  - `steppingAlgorithm` = `MscMuHadStepLimitType()` = `fMinimal` (default [EMP] L176; the printout's "Type of msc
    step limit algorithm for muons/hadrons 0");
  - `facrange` = `MscMuHadRangeFactor()` = 0.2 ([EMP] L161);
  - `facgeom` = 2.5 ([EMP] L162);
  - `facsafety` = 0.6 ([EMP] L164);
  - `lambdalimit` = 1 mm ([EMP] L165; [VMSC] L61);
  - `skin` = 3 ([EMP] L163 default 1.0, set to 3 by [OPT4] L127 for all particles);
  - `latDisplasment` = `MuHadLateralDisplacement()` = true ([EMP] L127 default **false**, set **true** by [OPT4] L129;
    printout "Msc lateral displacement for muons and hadrons 1");
  - "Urban msc model lateral displacement alg96 1" ([EMP] L126).
  All of these are copied by [VMSC] `InitialiseParameters` (L96–113, P2 F1), which overwrites any in-class
  initialiser. So the `G4VMscModel.hh` initialisers of `facrange`, `latDisplasment`, `skin`, `facgeom` and
  `steppingAlgorithm` are **moot** for GenericIon and alpha.
- **`dtrl` is not overwritten and not printed.** Its initialiser is in `G4VMscModel.hh`, which is not in the cache:
  OPEN (INFERRED 0.05, P2 §2.1). It sets the step fraction below which [URB] uses the linear energy-loss
  approximation in the true-to-geometric conversion ([URB] L708, L795). It affects only the Urban path-length
  conversion of the TOPAS reference (V2 Case C lineage).
- Displacement switch-off far from boundaries: [URB] L485–491 (P2 F1). The P4 step-length note (16 mm entrance steps
  for alpha) stands.

**λ charge scaling (from the cached sources: partly FACT, the decisive inline OPEN).**
- [URB] uses the ion's own definition: `SetParticle(track->GetDynamicParticle()->GetDefinition())` (L420), then
  `GetTransportMeanFreePath(particle, currentKinEnergy, …)` (L448, also L726, L842).
- In `ComputeCrossSectionPerAtom` the charge enters only the screened-Rutherford branch (L320
  `sigma *= chargeSquare*atomicNumber*atomicNumber/(beta2*bg2)`). The high-energy branch used by every V3-008 case
  (electron-equivalent energy > 10 MeV) overwrites `sigma` without `chargeSquare` (L399–409; P2 §2.3).
- [VMSC] applies no charge factor to λ (P2 F2).
- Whether the λ table built for GenericIon (charge 1) is rescaled for a C/O ion, by charge or by a scaled energy,
  is decided by the inline `G4VMscModel::GetTransportMeanFreePath` (`G4VMscModel.hh`), the inline
  `G4UrbanMscModel::SetParticle` (`G4UrbanMscModel.hh`) and the table builder (`G4VMultipleScattering.cc` /
  `G4EmTableUtil.cc`). **None of these is cached → OPEN.**
- Impact (P2 §6.3, F2, F5): λ does not enter θ0. It shifts the Urban in-window rms by about 1 % (scratch estimate)
  and affects only an Urban emulator, which r5 withdrew, and the attribution text of V2 (Case C, informative). No
  frozen row is graded on it.

**Disposition.** Final msc parameters: FACT (printout + source). In-class initialisers: moot, except `dtrl` (OPEN,
lineage only). λ charge scaling: OPEN at the inline level and **declared a lineage residual** of the TOPAS reference
(V2 Case C attribution). Settling files: `G4VMscModel.hh`, `G4UrbanMscModel.hh`, `G4VMultipleScattering.cc`.

## 4. icru90 z6/z8 silent-fallback check (P5 row)

**What happens when the file is absent (FACT [ICRU73], [LS]).**
- [ICRU73] L284 calls `RetrieveVector(ost, false)`.
- With `warn = false`, a file that cannot be opened produces neither a `G4Exception` nor a printout (L404–416). The
  vector stays null and `GetDEDX` returns 0 (L147).
- [LS] L281–291: if `res > 0.0` the data are used; otherwise `res = eloss*q2/chargeSquareRatio` ("simplified
  correction"), i.e. the Bragg-model loss rescaled to the current q².
- The only positive trace of a successful read is `"File " … " is opened by G4IonICRU73Data"`, printed only if
  `fVerbose > 0` ([ICRU73] L417–419). The other loader printouts need `fVerbose > 1` (L171–243).
- `fVerbose` is declared and initialised in `G4IonICRU73Data.hh` (not cached). The `.cc` never sets it from
  `G4EmParameters`, so whether `/process/em/verbose` reaches it is OPEN.
- The LS model's own printouts at these points are commented out ([LS] L275–279, L301–306, L309–313). The model
  has no verbose path.

**Archived logs (FACT).** grep for `icru90`, `ICRU90`, `ion_stopping`, `G4LEDATA`, `LindhardSorensen` in the carbon,
oxygen and carbon-SS stdout files:
- carbon/oxygen: `Use ICRU90 data 1` (L89); `Directory in G4LEDATA for fluorescence data files fluor` (L127);
  `LindhardSorensen : Emin= 0 eV Emax= 10 GeV` (L284); `ICRU73QO` lines for other particles (L317, L381, L445, L509).
- carbon SS: `Use ICRU90 data 0` (L89); `LindhardSorensen` (L258).
- No line names a file under `ion_stopping_data` or records that it was opened. stderr is empty for the carbon and
  oxygen smokes.
- The data directory is `TOPAS_G4_DATA_DIR = /opt/reference/g4data` (stdout L6). Its content is provisioned with the
  engine and was not inspected by the planner. The files exist in the G4EMLOW 8.8 archive (P5: z6 `d9d4809d…`,
  z8 `863898af…`).

**Conclusion: evidence ABSENT.** The archived logs cannot show whether `icru90/z{6,8}_G4_WATER.dat` were read, and
the code is silent on failure.
- What would produce the evidence:
  - (i) a scorer-free smoke in which `G4IonICRU73Data::fVerbose > 0`. That requires a setter or a link to the EM
    verbose level, which is OPEN (`G4IonICRU73Data.hh`); `/process/em/verbose 2` alone is not shown to reach it;
  - (ii) an engine-side listing with sha256 of `$G4LEDATA/ion_stopping_data/icru90/z{6,8}_G4_WATER.dat` in the
    provisioned data directory. This shows presence, not reading.
- The LS model has no own verbose.
- **Recorded as a declared limitation of the lineage record, not a blocker.**
  - No frozen row depends on it. I3 uses libdEdx ICRU73_NEW, not the Geant4 files.
  - The files act only at ≤ 2.0 MeV scaled, i.e. on the last ≲ 0.03 mm of C/O range (P2 F3), far below the V1 R80
    tolerance (0.35 %, ≈ 0.8–1.1 mm).
  - I3c (TOPAS attribution diagnostic) is exploratory.
  - The V5-SS ratio is a lateral statistic.
- **Hence no row is affected.**

## 5. Summary of the FC-2 disposition

| FC-2 sub-item | Outcome |
|---|---|
| Pin `G4BetheBlochModel.cc/.hh`; settle its ICRU 90 branch | DONE: pinned; branch settled for G4_WATER + opt4 (alpha → ICRU 90 alpha array ≤ 1 GeV; GenericIon → ICRU 90 proton array × q² + LS terms); one INFERRED premise (`G4ICRU90StoppingData` index/span); p3b §5 item resolved; possible z = 1 double count = lineage item |
| X0 (RadL) of G4_WATER | DONE: derived from the source composition, 36.083 g/cm²; ratio to ionmc's cross-check X0 = 1.0000000 |
| G4_WATER I-value | DONE: 78 eV, FACT ([NMB] L712, L230–235) |
| P4 "material table (RadL) archived" | CORRECTED: no material table is in the logs; archived = source definition + derived X0; optional `Ma/Verbosity` smoke (content unverified) |
| λ charge scaling; lateral displacement for ions | lateral displacement FACT (on; printout + [OPT4] L129); λ z²: OPEN at the inline level → declared lineage residual (V2 Case C attribution only) |
| `G4VMscModel.hh` initialisers | moot (overwritten; final values printed), except `dtrl` OPEN (lineage only) |
| icru90 z6/z8 silent-fallback check | evidence ABSENT; silent fallback confirmed in source; declared limitation of the lineage record; no row depends on it |

No row definition, tolerance or seed changes. No transport, qualification or scored reference result was produced
or read.

## 6. Scratch artefacts (session scratch, not committed)

- `v3-008/fc2/x0_water.py`, sha256 `f6005e0bae99b230856900947c7a30c73196919458e0a8baf698679514eb18cd`. Output: X0 Geant4-form 36.0830 g/cm², ionmc 36.0830 g/cm², ratio
  1.00000000, weight sensitivity 1.11e-04.
