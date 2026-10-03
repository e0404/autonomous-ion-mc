# Validation evidence for p / He / C / O transport: sources, reference engines, evaluation matrix

Status: the research was cut short at the coordinator's request. Each claim is labelled:
- **[V]**: verified in this session (I fetched the data, the API or the paper text).
- **[K]**: from domain knowledge, not re-verified here. Check it before freezing any decision.
- **[NV]**: not verified. It could be wrong.

## 0. Key findings

1. **The Geant4 validation database geant-val.cern.ch has a public JSON API with *measured* carbon fragmentation numbers [V].** It holds Haettner et al. (12C at 400 MeV/u in water, GSI). The data cover fragment yields N/N0 for Z=1–5 at water depths of 59, 159, 258, 279, 288, 312 and 347 mm, with errors. There are also angular distributions and energy spectra at 0–8°. It also holds EXFOR-derived inelastic cross sections for p+C, p+O and 12C+C, the Barashenkov compilation, and Toshito 2007 (12C charge-changing cross sections in water). This is the strongest public numeric evidence for carbon attenuation, build-up and the fragment tail. Only 400 MeV/u is present; the 200 MeV/u Haettner data are not in the API.
2. **Gottschalk et al. 2015 (177 MeV proton halo) has about 300 *absolute* dose points in a table [V].** The values are in MeV/g/p, at depths z = 1.5–~22 cm and radial offsets r = 0–10 cm (arXiv:1409.1938, Table 1, extractable text). It is the best public proton measurement for both depth dose and the nuclear halo.
3. **Proton MCS has measured θ0 data [V].** Gottschalk 1993 (Go93) measured 158.6 MeV protons through 14 materials in 115 material/thickness combinations. Makarova, Gottschalk & Sauerwein (arXiv:1610.01279) tabulate a subset. They show that Molière/Fano/Hanson theory matches the measurements to <1% on average. They also show that Geant4 Urban MSC is about 8% low for low-Z targets and WentzelVI about 4% low. **For MCS, TOPAS is therefore weaker evidence than Hanson theory plus Go93.**
4. **I found no public numeric measured Bragg curves or lateral profiles for He, C or O in water.** I searched for HIT/Tessonnier 2017, the GSI Steidl/Schardt precision Bragg curves and Zenodo/Figshare. HIT's Tessonnier 2017 is not open access and has no data statement [V]. For He and O the measured numeric evidence is therefore limited to nuclear cross sections (Horst 2019 for He, Zeitlin 2011 for O) plus ion-specific stopping tables. A reference request should be prepared (section 5).
5. **The reference engines share stopping-power lineage.** TOPAS uses G4 ICRU90/ICRU73 data [K], and MCsquare and FRED are proton-tuned and PSTAR-like [K/NV]. Agreement among the engines on range is therefore *not* independent of the I-value choice. Range claims need a separately documented I-value and measured anchors.

## 1. Measured and tabulated data with usable numbers

| Topic | Source | Where the numbers are | Category | Access / licence |
|---|---|---|---|---|
| p depth dose plus halo, absolute | Gottschalk, Cascio, Daartz, Wagner, PMB 60 (2015) 5627; [arXiv:1409.1938](https://arxiv.org/abs/1409.1938) | Table 1: log10(dose/(MeV/g/p)) on a z × r grid, 177 MeV, r = 0..10 cm. Appendix D gives the beam parameters (σx, emittance) needed to model the source. | measured | arXiv non-exclusive licence. Numbers are facts and can be used with citation. The BGware.ZIP links in the paper return 502/403 [V]. |
| p MCS θ0 | Go93: Gottschalk et al., NIM B 74 (1993) 467, [doi:10.1016/0168-583X(93)95944-Z](https://doi.org/10.1016/0168-583X(93)95944-Z) | Full table in Go93 (paywalled). Subset in [arXiv:1610.01279](https://arxiv.org/abs/1610.01279), Table 1: Be, polystyrene, C, Lexan, …, Pb, U, with θ0 exptl / Hanson / Urban / Wentzel. | measured plus independent theory | arXiv subset is usable now. The full Go93 table needs the paywalled paper (intervention category 1). |
| p scattering-power theory | Gottschalk, Med Phys 37 (2010) 352; [arXiv:0908.1413](https://arxiv.org/abs/0908.1413); "Techniques of proton radiotherapy: transport theory", [arXiv:1204.4470](https://arxiv.org/abs/1204.4470) | Formulas plus Tables 4/5 (θ for 8 models) | independent theory | arXiv |
| p/He stopping and CSDA range | NIST PSTAR/ASTAR (ICRU 49, I_water = 75 eV), https://physics.nist.gov/PhysRefData/Star/Text/ | HTML form (POST) output tables | ion-specific tabulated | US Government work, public. Cache the POST body and bytes. |
| Water I-value update | ICRU 90 (2016): I_water = 78 ± 2 eV [V via search]. New proposal of 79.4 eV: de Vera et al., [arXiv:2608.15368](https://arxiv.org/abs/2608.15368) [V abstract]. | ICRU 90 report is paywalled. The G4 G4ICRU90StoppingData class contains the p/α tables for water, air and graphite [K]. | tabulated / theory | ICRU is paywalled. Geant4 licence allows reuse but lineage is shared with TOPAS. |
| C/O stopping | ICRU 73 (+2009 erratum, I = 78 eV), MSTAR (Paul) [K] | report / program | ion-specific tabulated | paywalled / program licence [NV]. The libdEdx package bundles PSTAR/ASTAR/MSTAR/ICRU73 [K, licence NV]. |
| p depth dose, independent MC | Berger, NISTIR 5226 (1993), PTRAN, "Penetration of proton beams through water I", https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir5226.pdf | Tables of depth dose for 50–250 MeV in water [K: tables exist; not opened] | independent MC (older nuclear data) | US Government, public |
| p nonelastic cross sections | EXFOR (IAEA NDS, https://www-nds.iaea.org/exfor/); geant-val record inspire_id=-7 (p+C, p+O, p+Al, p+Ca, 12C+C, 12C+Al) [V]; Barashenkov compilation (inspire_id=-2) [V]; ICRU 63 [K]; Carlson ADNDT 63 (1996) 93 [K] | JSON API / EXFOR retrieval | measured | EXFOR is free. geant-val states no licence [V]. Cite the original experiments. |
| 12C thick-target fragmentation in water | Haettner et al., PMB 58 (2013) 8265, [doi:10.1088/0031-9155/58/23/8265](https://doi.org/10.1088/0031-9155/58/23/8265) | **geant-val** `GET https://geant-val.cern.ch/api/getExpPlotsByInspireId?inspire_id=-5`. 384 records: 5 "fragment yield" curves (Z = 1..5 vs depth, N/N0 ± err), 7 yield-by-Z sets, 34 angular, 152 energy spectra; 400 MeV/u only [V]. Example: Z=2 at 279 mm gives 0.570 ± 0.011. | measured | No licence on geant-val. Treat it as citing published experimental facts and keep it in the provenance cache. |
| 12C charge-changing cross sections in water | Toshito et al., PRC 75 054606 (2007) (geant-val inspire_id 747302) [V] | API | measured | as above |
| 12C thin-target double-differential cross sections, 95 and 50 MeV/u, on H, C, O, Al, Ti | Dudouet et al., PRC 88 024606 (2013); Divay et al., PRC 95 044602 (2017); https://hadrontherapy-data.in2p3.fr/ [V: site lists E600, Zero Degree and 50 MeV/u experiments, "free access"] | Website (file format not verified) | measured | free access; licence not stated |
| 12C neutron yields, 62 MeV/u | INFN-LNS (geant-val -12) [V] | API | measured | as above |
| 16O (and 14N, 20Ne, 24Mg) fragmentation, 290–1000 MeV/u on C, CH2, Al, … | Zeitlin et al., PRC 83 034909 (2011), [arXiv:1102.2848](https://arxiv.org/abs/1102.2848) [V] | Tables in the paper. H-target values come from CH2 − C subtraction. | measured | arXiv |
| 12C fragmentation, 290/400 MeV/u | Zeitlin et al., PRC 76 014911 (2007) [K] | Tables | measured | APS paywall; arXiv version NV |
| 4He charge- and mass-changing cross sections on H, C, O, Si, 70–220 MeV/u | Horst et al., PRC 99 014603 (2019), [doi:10.1103/PhysRevC.99.014603](https://doi.org/10.1103/PhysRevC.99.014603) [V: exists] | Tables in the paper [NV] | measured | APS; possibly open access [NV] |
| Total reaction cross-section database (space/therapy) | Luoni et al. 2021, [arXiv:2105.11981](https://arxiv.org/abs/2105.11981) [V: exists] | Compiled σ_R database (p, He, C, O on many targets) [NV content] | measured compilation | arXiv / data licence NV |
| He/C/O Bragg curves, peak widths, FWHM at HIT | Tessonnier et al., PMB 62 (2017) 3958, [doi:10.1088/1361-6560/aa6516](https://doi.org/10.1088/1361-6560/aa6516) | Not open access at LMU, no data statement [V]. Figures only. | measured (inaccessible numbers) | intervention candidate |
| p, 3He, 7Li, 12C, 16O precision Bragg curves, 100–400 MeV/u (GSI) | Steidl, Schardt, Weber, Krämer, GSI report / DPG abstracts [V: existence] | Numbers not retrieved. repository.gsi.de is behind a bot challenge [V]. | measured | intervention candidate |
| p 67.5 MeV depth dose (UCSF) | Faddegon et al., Med Phys 42 (2015) 4199 | Supplementary availability NV | measured | NV |
| LET measured | Microdosimetry (TEPC: Kase 2006/2013, Martino 2010; silicon microdosimeters; FNTD: Granville 2016) [K] | Mostly figures. These give y-spectra, not LETd, and converting to LETd needs a model. | related model | Not recommended as an acceptance anchor |

**Missing measured numeric data:**
- No public numeric measured Bragg curves for He, C or O in water.
- No public numeric lateral pencil-beam profiles for He, C or O.
- No public Haettner 200 MeV/u data.
- No directly usable measured LETd for any species.

## 2. Reference engines: roles and native configuration

### TOPAS 4.3 / Geant4 11.4.2 (p, He, C, O; the only engine covering all four species with fragmentation)
- **Physics [K]:** use `Geant4_Modular` with `g4em-standard_opt4`, `g4h-phy_QGSP_BIC_HP` (or `_BIC`), `g4decay`, `g4ion-binarycascade`, `g4h-elastic_HP` and `g4stopping`. This is the TOPAS default. Run variants with the ion module replaced by **`g4ion-QMD`** and **`g4ion-inclxx`**; the spread between them is a model-uncertainty band.
  - Bolst et al. (NIM A 869, 2017) and G4-Med (Arce et al., Med Phys 2021 and 2025) found that BIC, QMD and INCL++ differ by tens of percent in fragment yields and angular distributions against Haettner [K; G4-Med confirmed to use Haettner for the BP and fragment tests [V]].
  - Geant4 ≥11.1 improved p and C Bragg peaks through ICRU90 low-energy proton stopping and the Lindhard–Sørensen ion model [V, G4-Med 2025].
  - Use a range cut of ≤0.05–0.1 mm, or 1 mm when dose-only scoring at ≥1 mm voxels. Enable `/process/em/printParameters` and archive the log to confirm the MSC model (opt4: WentzelVI for protons [K]) and ICRU90 use.
- **Scorers:**
  - Dose and energy: `DoseToWater`, `DoseToMedium` and `EnergyDeposit` on TsBox/TsCylinder grids.
  - Fluence: `Fluence` and `EnergyFluence`, plus `SurfaceTrackCount` on a surface with `OnlyIncludeParticlesGoingIn`.
  - Species and charge filters: `OnlyIncludeParticlesNamed`, `OnlyIncludeParticlesOfAtomicNumber`, `OnlyIncludeParticlesOfCharge` and `OnlyIncludeIfIncidentParticlesNamed`. Use them for Z-resolved yields comparable to Haettner (count Z-filtered particles crossing planes at the Haettner depths, within an angular cone) [K].
  - Energy binning: `EBins`, `EBinMin` and `EBinMax` give species-resolved spectra [K].
  - LET: `ProtonLET` (TsScoreProtonLET) gives dose-weighted LETd by default, or track-weighted. It covers protons only (primary plus secondary) and uses unrestricted electronic stopping power from the pre-step energy, following Granville & Sawakuchi 2015, with options such as `MaxScoredLET` [K]. **For He/C/O, do not rely on a built-in mixed-field LET [NV whether 4.3 has one].** Instead, score per-species fluence-energy spectra per voxel and compute LETd = Σ Φ_i(E)·S_i(E)² / Σ Φ_i(E)·S_i(E) offline with a declared stopping table. This keeps the estimator definition controlled and identical across codes.
- **Uncertainty:** `Report` options are Sum, Mean, Histories, Count_In_Bin, Second_Moment, Variance and Standard_Deviation [K]. Whether Standard_Deviation is per-history or of the mean must be checked in the 4.3 docs [NV]. **Recommendation:** run ≥10 independent seeds or batches and compute the batch standard error of the mean yourself. This avoids definitional ambiguity and matches the project's batch method.
- **Output:** CSV with `#` header lines (scorer, bins, units) followed by rows `ix, iy, iz, value...` [K]. Parse the header. Do not assume the column order.
- **Cost [estimate, measure in calibration]:** about 1–3 ms per 200 MeV proton history per thread with opt4 and 0.1 mm cuts, i.e. roughly 2–5 min per 1e5 per thread. 290 MeV/u 12C with fragmentation costs about 10–50× more per history, i.e. roughly 0.5–3 h per 1e5 per thread. Run multithreaded if the service allows (`Ts/NumberOfThreads`).
- **Role:**
  - Independent MC for p (depth dose, halo, LETd) and for He/C/O (attenuation, fragment yields, mixed LET, tail).
  - Not authoritative for MCS θ0 below about 4% (Makarova 2017).
  - Its ion nuclear models have documented 10–30% discrepancies, so use the model-variant band, not one list.

### MCsquare (e0404 binary): protons only
- **Physics:** Souris, Lee & Sterpin, Med Phys 43 (2016) 2584 [K]. Class II condensed history and Rossi-type MCS. ICRU-63-type nuclear cross sections. Transports secondary p, d and α; heavier fragments are deposited locally; neutrons and γ escape [K].
- **Config keys:** from the smoke test, `Simulate_Nuclear_Interactions` and `Simulate_Secondary_*`. Additional keys to verify against the bundled `config` template or the binary's help [NV]:
  - `LET_MHD_Output True` and `LET_Calculation_Method StopPow|DepositedEnergy`;
  - `Compute_stat_uncertainty True`, `Stat_uncertainty <target>`;
  - `Energy_ASCII_Output`, `Dose_Sparse_Output`, `Beamlet_Mode`, `Export_batch_dose` (or a similar batch key), `Max_Num_Primaries`.
- **Monoenergetic point beam [K/NV]:** write a custom BDL file with one or two nominal energies:
  - MeanEnergy = NominalEnergy and EnergySpread = 0;
  - one Gaussian (Weight1 = 1) with SpotSize ≈ 1e-3 mm, Divergence 0 and Correlation 0.
  - First test whether exactly 0 is accepted (σ = 0 may trigger divide-by-zero in the emittance parametrisation), so use 1e-3 mm / 1e-6 rad.
  - Set the nozzle-to-isocentre distances so that the source plane lies at the phantom surface.
  - Verify the result by scoring a thin entrance slab: lateral σ ≈ 0 and peak width consistent with pure straggling.
- **Role:** second independent proton MC (different MCS, nuclear and code lineage from Geant4) for depth dose, R80, lateral profiles in water/heterogeneities and proton LETd. It provides no evidence for ions.

### FRED 3.76 (CPU, -nogpu)
- Proton physics: Schiavi et al., PMB 62 (2017) 7482 [K]. The Fred group has published carbon-ion extensions [K]; **whether the 3.76 binary transports He/C/O and with which nuclear model was NOT verified.**
- Action: run `fred -h` and inspect the manual shipped with the binary. Smoke-test a pencil beam with `particle = C12` (or the documented equivalent) and **inspect the log for actual ion transport and nuclear settings** before assigning any ion role.
- Output: `out/score/*.mhd` [V from smoke]. LET scoring in Fred exists for protons [K]; exact keys NV.
- **Provisional role:** a third proton MC. Use it for ions only after verification.

## 3. Evaluation matrix

How to read the tolerances: they are proposals to freeze before qualification. They are taken from published agreement levels; where I could not verify a number it is marked [K]. Tolerances apply to the difference after the reference's own uncertainty is accounted for (section 4).

| Capability | Observable(s), unit, normalisation | Best evidence (category) | Proposed tolerance and justification | If sources disagree |
|---|---|---|---|---|
| **p stopping / range** | R80, R90 of the IDD of a monoenergetic pencil beam in water (mm), 70–230 MeV; d(R)/dE | PSTAR CSDA with I = 75 eV and ICRU 90 with I = 78 eV (tabulated); TOPAS and MCsquare (independent MC); Gottschalk 2015 Table 1 at 177 MeV (measured) | **R80 vs CSDA with the *same* I-value: ≤0.3% or 0.3 mm**, whichever is larger. R80 ≈ R_CSDA for a pristine peak to ~0.1–0.2% [K: Bortfeld 1997]. Between-code R80: ≤0.5 mm. The 75→78 eV choice alone shifts range by ~0.5–0.7% [estimate]. | Declare the I-value. Treat 75 vs 78 eV as a documented systematic band, not as an error. Measurement decides between them. |
| **p straggling** | Distal 80–20% fall-off (mm); peak-to-plateau ratio | Bohr/Tschalär theory (independent theory); TOPAS and MCsquare (MC) | Fall-off within 5% or 0.1 mm of the MC mean. Peak/plateau within 2%. | Theory first for an energy-loss-only case (turn nuclear off in all codes). |
| **p MCS** | θ0 (mrad) after slabs; σ_x(z) in water (mm) | Go93 θ0 plus Hanson (measured plus theory); Preston–Koehler / Fermi–Eyges for σ_x(z) (theory) | θ0 within 2% of Hanson for thin-to-thick slabs. Hanson matches Go93 to <1% on average [V]. σ_x(z) within 3% or 0.2 mm. | Prefer Hanson/Go93 over TOPAS: G4 WentzelVI is ~4% low and Urban ~8% low for low-Z targets [V]. |
| **p nuclear attenuation / halo** | Primary fluence vs depth (N/N0); absolute dose in MeV/g/p on Gottschalk's z × r grid | EXFOR/ICRU63 σ_nonel for p+O and p+H (measured); Gottschalk 2015 Table 1 (measured); TOPAS BIC (MC; Hall et al. 2016 report agreement over 5 orders of magnitude [V abstract]) | Core (r = 0): 2%. Halo (r ≥ 1 cm): 15% per point. The fit residuals in Gottschalk 2015 are 9% (MI) and 15% (MD) rms [V]. Attenuation slope: within 3% of σ-based prediction. | Check the beam model first (σx and energy from Appendix D), then the nuclear cross section. |
| **p LETd** | LETd (keV/µm) along the axis; definition: dose-weighted, unrestricted electronic stopping, primaries plus secondary protons, with the energy-cut handling stated | Analytic: LETd ≈ S(E(z)) in the plateau for primaries (theory); TOPAS ProtonLET and MCsquare LET (MC) | Plateau: 3%. Bragg peak and distal: 10% or 1 keV/µm, with the estimator definition matched. Published estimator choices (Cortés-Giraldo & Carabe PMB 2015; Granville & Sawakuchi PMB 2015) change distal LETd strongly [K]. | Recompute all codes with the *same* offline estimator from fluence spectra before declaring a disagreement. |
| **He range** | R80 for 4He in water (mm) | ASTAR (ICRU 49, α) with an I-value correction (tabulated); TOPAS (MC); GSI precision Bragg curves (not retrieved) | 0.3% or 0.3 mm vs table with the same I | Measured HIT/GSI numbers needed (intervention). |
| **He nuclear** | Charge-changing σ (mb); He fluence attenuation vs depth; secondary Z=1 build-up | Horst 2019 (measured); TOPAS BIC/QMD/INCL (MC) | σ within combined uncertainty or 5%. Attenuation slope within 5%. | Use the measured σ. Report the model band. |
| **C range / Bragg peak** | R80 (mm) and peak/plateau, 12C 100–430 MeV/u | ICRU 73 / MSTAR (tabulated); TOPAS (MC); measured BP positions at GSI/HIT (not public) | Range: 0.3% or 0.3 mm vs table with the same I. Peak/plateau: 5% vs TOPAS (sensitive to the nuclear model). | Intervention for measured curves. |
| **C attenuation / fragments / tail** | Primary N/N0 vs depth; Z = 1..5 yields N/N0 at 59–347 mm; angular distributions; tail-to-peak dose ratio | **Haettner 400 MeV/u via geant-val** (measured); Toshito 2007 and Zeitlin 2007 charge-changing σ (measured); Dudouet 95 MeV/u thin target (measured); TOPAS model band (MC) | Primary attenuation: 3%. Z = 1, 2 yields: 15%. Z = 3–5 yields: 25%, or within combined experimental error. G4-Med/Bolst-type comparisons show tens-of-percent model spreads [K]. Tail dose: within the TOPAS BIC/QMD/INCL envelope plus 10%. | Measured data win. If the implementation sits inside the measurement error but outside the TOPAS band, report TOPAS model inadequacy and do not tune to TOPAS. |
| **C mixed-field LET** | LETd_total and species-resolved LETd (keV/µm) | Offline LETd from species fluence spectra in TOPAS (MC); ICRU 73 stopping (tabulated) | Plateau: 5%. Peak: 10%. Tail: 20%, with the same estimator. | Same estimator first. Species decomposition next. |
| **O range / attenuation** | R80; σ_cc (mb); N/N0 | ICRU 73 / MSTAR (tabulated); **Zeitlin 2011** (16O on C and H via CH2−C, ≥290 MeV/u) (measured); Luoni 2021 σ_R database; TOPAS (MC) | Range: 0.3%/0.3 mm. σ: within combined uncertainty or 5%. Fragment claims only with a documented domain. | Restrict claims to the energies and targets that have measurements. |

## 4. Statistical comparison method

- **Primary test: pointwise standardized residuals.** z_i = (x_i − r_i) / sqrt(σ_x,i² + σ_r,i² + σ_sys,i²). Use batch standard errors with ≥10 batches; for an MC reference, the reference's own batch standard error.
  - Acceptance needs both:
    - an **equivalence (bias) test**: |x − r| ≤ tolerance + 2σ_comb in a frozen ROI. This is TOST-style and avoids "passing" because of noise.
    - a **consistency test**: the χ²/ndf of z over the ROI is consistent with 1 (p > 0.01), and coverage of the 1σ intervals is about 68%.
  - Because the dose fall-off is steep, evaluate high-gradient regions in *position* (distance-to-agreement) rather than dose.
- **Range metrics.**
  - Fit or spline-interpolate the distal edge; avoid raw bin values. Use a cubic or monotone spline, or fit the Bortfeld curve.
  - Report R90, R80, R50 and the 80–20 fall-off.
  - Bootstrap over batches to get σ(R80). This is typically 0.02–0.05 mm for 1e6 protons at 0.1–0.2 mm bins [estimate].
  - Depth bins must be ≤ 1/5 of the fall-off width: ≤0.2 mm for protons, ≤0.1 mm for carbon near the peak. Report how the bin centre is defined.
- **Gamma index.** Useful as a summary for 3D/2D clinical-like cases: 1%/1 mm for MC-vs-MC water, 2%/2 mm for MC-vs-measurement, with a 10% low-dose threshold. It must not be the sole acceptance test, because noise in the evaluated distribution inflates the pass rate and it hides systematic range shifts. Always pair it with R80 and z-tests. Compute gamma on a denoised or high-statistics reference, and state the local or global normalisation.
- **Normalisation choice is part of the claim.**
  - Absolute per-primary (MeV/g/p, Gy·cm²) comparisons are the most discriminating. Use them for MC-vs-MC and for Gottschalk 2015.
  - Normalising to the peak hides nuclear-attenuation errors; normalising to the plateau hides peak-height errors. State which one is used and test the other as a secondary.
  - Haettner yields are normalised to primaries incident (N/N0) and cut on acceptance angle. The simulation must apply the same angular acceptance and detector distance.
- **Common pitfalls.**
  - Comparing IDD versus central-axis dose, or finite-radius chamber IDD (Bragg peak chamber radius truncates the halo by a few percent, per Gottschalk 2015).
  - Different I-values.
  - Comparing a "dose to water" scorer with "dose to medium".
  - Different secondary-particle cuts for LET.
  - Statistical correlation across voxels from the same histories, so use batch-level, not voxel-independent, resampling for integral metrics.
  - Treating reference MC noise as zero.
  - Re-tuning to the evaluation set. Freeze the calibration and evaluation roles before running.

## 5. Gaps and proposed intervention or reference requests

After bounded public search, the following measured **numbers** were not obtained. Each needs a structured intervention (category 1: inaccessible information) once two access approaches have been documented as failed.
1. **HIT Bragg curves, peak positions, FWHM and lateral profiles for 1H, 4He, 12C, 16O** (Tessonnier et al. 2017, doi 10.1088/1361-6560/aa6516, and the companion MC-validation paper PMB 62 (2017) 6579). Request the operator to supply the paper tables, supplementary files, or digitised curves with a stated digitisation uncertainty.
2. **GSI precision Bragg-curve measurements (Steidl/Schardt)** for p, 3He, 7Li, 12C and 16O at 100–400 MeV/u. The repository is behind a bot challenge [V]. Request the report PDF.
3. **Full Go93 θ0 table** (paywalled NIM B), to supplement the arXiv subset.
4. **Haettner 200 MeV/u data and the attenuation table** from PMB 58 8265 (only 400 MeV/u is in geant-val).
5. **Horst 2019 tables** if the APS version is not accessible.
6. **If He/O ion transport in FRED is unavailable**, the only independent ion MC is Geant4/TOPAS. Consider a plug-and-play **FLUKA** reference package (FLUKA needs licence registration, i.e. operator action). It would give an ion nuclear-model lineage independent of Geant4, for 12C at 290 MeV/u and 16O at 300 MeV/u depth dose and Z-resolved yields.

Until these arrive, He and O claims should be scoped to:
- range against ion-specific tables, plus the TOPAS model band;
- nuclear cross sections against Horst (He) and Zeitlin (O).

Release records must state explicitly that measured depth-dose evidence for He and O is pending or unavailable.

## 6. Immediate actions for the lead
- Acquire through the provenance cache: geant-val `api/getExpPlotsByInspireId?inspire_id=-5` (Haettner), `=-7` (EXFOR σ), `=747302` (Toshito) and `=-2` (Barashenkov); arXiv 1409.1938, 1610.01279 and 1102.2848 (for table extraction); PSTAR/ASTAR POST outputs; NISTIR 5226.
- Smoke-test the reference engines:
  - FRED with a C12 primary;
  - MCsquare with a custom zero-spread, tiny-spot BDL;
  - TOPAS with species-filtered fluence spectra and the BIC/QMD/INCL variants.
- Measure per-history cost to fix the run budgets.
- Freeze the matrix above, with tolerances and normalisations, in a decision record before running any qualification comparisons.
