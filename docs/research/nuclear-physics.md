# Nuclear interaction physics for IonMC v3: research report

Status: delivered early because the coordinator asked for it. Each item has one of these tags:
- **[V]** I checked it against a primary source or the reference source code during this session.
- **[M]** It comes from memory or the secondary literature and I did not check it in this session. Verify it before you freeze it in a decision record.
- Physics labels: **EST** = established physics, **APPROX** = approximation, **ENG** = engineering choice.

Downloaded source files, for provenance only. They are not project inputs: `.../scratchpad/research/{tp3621,tp209726,nucfrg2_tp3533,nucfrg3}.pdf/.txt` and the Geant4 `G4Tripathi*`, `G4IonsSihver/Kox/Shen` `.cc` files.

---

## 0. Executive recommendations

| Topic | Recommendation | Confidence |
|---|---|---|
| sigma_R, ion-nucleus (He, C, O on H, C, N, O, P, Ca) | Tripathi universal model (TP-3621), plus the light-system branch (TP-1999-209726) when either partner has A<=4. Compare with Sihver-1993 and Kox. Calibrate a target-specific scale factor only against measured data that is held out from evaluation. | High that it is adequate (~5-10%) |
| sigma_nonelastic, proton | A smooth fit or table to measured p+12C and p+16O data (Carlson 1996, EXFOR). Use Tripathi-light as the independent cross-check. ENDF/B-VIII.0 proton sublibrary (LA150) as the tabulated source <=150 MeV. | Medium-high |
| Proton secondaries | Event-level model from ENDF/B-VIII.0 MF6 energy-angle yields (p, d, t, 3He, alpha, recoil; Kalbach-Mann angular form). Transport p, d and alpha; deposit recoils locally; neutrons and gammas escape and are tallied. | Medium |
| He/C/O fragmentation | A **table-driven semi-empirical fragmentation kernel**. Per (projectile, target element, energy) it stores partial production cross sections sigma(Z_F, A_F), calibrated to measured charge-changing partial cross sections. Kinematics from Goldhaber/Morrissey. Target fragments are deposited locally and neutrons escape. Option: generate the tables with abrasion-ablation (NUCFRG2-style, public) and correct them to data. | Medium |
| Validation | Primary (Z=P) attenuation vs depth, per-Z fragment build-up vs depth, distal tail dose, fragment angular distributions and lateral halo. Use measured data first (Haettner, Schall, Toshito, Zeitlin, Tessonnier, Rovituso, Horst) and Geant4-QMD/BIC/INCLXX as a second line with a spread band. | High |
| Correction to the brief | At **400 MeV/u only ~25-35% of 12C survive to the Bragg peak**, not 50-70%. At 200 MeV/u ~65-70% survive. | Medium-high (my calculation agrees with [M] literature) |

---

## 1. Total reaction (inelastic) cross sections

### 1.1 Tripathi universal parameterization [V]
Sources (both NASA reports are US-government works, so the formulas are freely usable):
- R.K. Tripathi, F.A. Cucinotta, J.W. Wilson, *Universal Parameterization of Absorption Cross Sections*, NASA TP-3621 (Jan 1997). https://ntrs.nasa.gov/citations/19970011098. Journal version: NIM B 117 (1996) 347, doi:10.1016/0168-583X(96)00331-X.
- *Universal Parameterization of Absorption Cross Sections: Light Systems*, NASA/TP-1999-209726. https://ntrs.nasa.gov/citations/20000004258. Journal version: NIM B 155 (1999) 349, doi:10.1016/S0168-583X(99)00479-6.
- Neutron extension: NIM B 129 (1997) 11 and NASA TP-3656.

**EST/APPROX (empirical, physically motivated):**

sigma_R = pi r0^2 (A_P^{1/3} + A_T^{1/3} + delta_E)^2 (1 - R_c B/E_cm) X_m

- r0 = 1.1 fm. E_cm is the c.m. kinetic energy in MeV; Geant4 uses invariant-mass kinematics.
- B = 1.44 Z_P Z_T / R [MeV, R in fm].
- R = r_P + r_T + 1.2 (A_P^{1/3} + A_T^{1/3}) / E_cm^{1/3}.
- r_i = 1.29 r_rms,i. Use measured rms radii (de Vries 1987 ADNDT 36:495, charge radii with the nucleon form factor removed).
- delta_E = 1.85 S + 0.16 S / E_cm^{1/3} - C_E + 0.91 (A_T - 2Z_T) Z_P / (A_T A_P). For light systems with A_P > A_T, swap P and T in the isospin term.
- S = A_P^{1/3} A_T^{1/3} / (A_P^{1/3} + A_T^{1/3}).
- C_E = D [1 - exp(-E/T1)] - 0.292 exp(-E/792) cos(0.229 E^{0.453}). Here E is the lab energy in MeV/u.

The D and T1 parameters:
- **Heavy systems (TP-3621):** D = 1.75 (rho_AP + rho_AT)/(rho_C12 + rho_C12), with hard-sphere densities rho = A / (4/3 pi r_i^3). T1 = 40. Note that Geant4's `G4TripathiCrossSection` hard-codes D = 1.75 and uses r_rms = 0.6 x 1.36 fm x A^{1/3}. Both are simplifications of the paper.
- **Proton in TP-3621:** D = 2.05.
- **Lithium:** D/3.
- **Light systems (TP-1999):**
  - n+X: T1 = 18, D = 1.85 + 0.16/(1+exp((500-E)/200)).
  - p+X: T1 = 23, same D as n+X.
  - d+X: T1 = 23, D = 1.65 + 0.1/(1+exp((500-E)/200)).
  - 3He+X: T1 = 40, D = 1.55.
  - 4He+X: D = 2.77 - 8.0e-3 A_T + 1.8e-5 A_T^2 - 0.8/(1+exp((250-E)/G)).
  - (T1, G) for 4He: (40, 75) in general. Special targets: alpha+alpha (40, 300), Be (25, 300), N (40, 500), Al (25, 300), Fe (40, 300).
- **Coulomb multiplier R_c (Table 2, TP-1999):** p+d 13.5, p+3He 21, p+4He 27, p+Li 2.2, d+d 13.5, d+4He 13.5, d+C 6.0, 4He+Ta/Au 0.6. Otherwise R_c = 1. TP-3621 instead raised the barrier for p+4He by x27 and for p+12C by x3.5.
- **X_m (optical-model multiplier, used for neutrons and light systems):** X_m = 1 - X1 exp(-E/(X1 S_L)). X1 = 2.83 - 3.1e-2 A_T + 1.7e-4 A_T^2 (n+4He: X1 = 5.2). S_L = 1.2 + 1.6 [1 - exp(-E/15)]. For charged projectiles above ~30 MeV/u, X_m is about 1.
- **Reference implementation:** Geant4 `G4TripathiLightCrossSection.cc` (tag v10.7.4; removed from current master; Geant4 licence). Its radii come from `G4WilsonRadius` and include an extra sqrt(5/3) factor. This only affects the Coulomb term, which is negligible above ~20 MeV/u.

**My implementation check (TP-3621/TP-1999 formulas, de Vries radii; per element, mb):**

| Proj | E [MeV/u] | sigma(+H) | sigma(+C) | sigma(+O) | lambda_water [cm] |
|---|---|---|---|---|---|
| p | 100 / 200 | — | 270 / 242 | 325 / 295 | 92 / 102 |
| 4He | 100 / 200 | 96 / 93 | 384 / 353 | 466 / 432 | 46 / 48 |
| 12C | 100 / 200 / 400 | 270 / 242 / 239 | 955 / 867 / 870 | 1086 / 994 / 999 | 18.4 / 20.2 / 20.3 |
| 16O | 200 / 400 | 295 / 292 | 994 / 999 | 1134 / 1140 | 17.4 / 17.3 |

How these compare with expectations:
- C+C at ~870 mb agrees with measurements (Kox 1987; Jaros 1978 ~ 850-880 mb) [M].
- p+C at 240-270 mb looks ~5-15% high compared with my recollection of ICRU 63/Carlson (~220-230 mb) [M]. **Check this against data before adopting Tripathi for protons.**
- The 4He+C value (~350-380 mb) needs checking against Horst 2017/2019. The 4He D-sigmoid makes sigma rise above 250 MeV/u, which is outside the clinical He range (<=~250 MeV/u).

### 1.2 Alternatives [V for formulas from Geant4 source]
- **Sihver et al. 1993**, PRC 47:1225, doi:10.1103/PhysRevC.47.1225. sigma = pi r0^2 [A_P^{1/3} + A_T^{1/3} - b0 (A_P^{-1/3} + A_T^{-1/3})]^2 with r0 = 1.36 fm and b0 = 1.581 - 0.876 (A_P^{-1/3} + A_T^{-1/3}). [M] The proton form uses b0 = 2.247 - 0.915 (1 + A_T^{-1/3}). It is energy-independent and valid above ~100 MeV/u. The same paper gives **partial cross-section** formulas, which is useful for section 3.
- **Kox et al. 1987**, PRC 35:1678. sigma = pi R_int^2 (1 - B_c/E_cm). R_int = r0(A_P^{1/3} + A_T^{1/3}) + r0[a A_P^{1/3} A_T^{1/3}/(A_P^{1/3} + A_T^{1/3}) - c(E)] + 5(A_T - 2Z_T) Z_P/(A_P A_T) [fm]. r0 = 1.1, a = 1.85, B_c = Z_P Z_T / (1.3 (A_P^{1/3} + A_T^{1/3})). Geant4's c(E) fit is c = 2 - 10/(log10 E)^5 for log10 E > 1.5.
- **Shen et al. 1989**, NPA 491:130. A Kox-like form with an added E_cm^{-1/3} surface term. This was the Geant4 ion default before Glauber-Gribov.
- **Bradt-Peters 1950**, Phys Rev 77:54. sigma = pi r0^2 (A_P^{1/3} + A_T^{1/3} - b)^2 with r0 ~1.2-1.4 fm and b ~ 0.8-1. It has no energy dependence; use it only as a sanity check.
- **Glauber-Gribov** (`G4ComponentGGNuclNuclXsc`, Grichine). This is the current Geant4 ion default [M]. It makes TOPAS a non-identical cross-section lineage to Tripathi, which is good for independence.
- **Luoni et al. 2021**, New J Phys 23:101201. A total nuclear reaction cross-section database for space and therapy, with data and a model comparison [M]. Good for evidence assembly; check the data files and licence.

### 1.3 Proton-nucleus data access
- **ICRU Report 63 (2000).** Copyrighted and sold by ICRU/SAGE; **not free** [M]. Its proton data come from the LA150 evaluation (Chadwick et al., Nucl Sci Eng 131 (1999) 293). LA150 is distributed in the free **ENDF/B-VIII.0 proton sublibrary** at https://www.nndc.bnl.gov/endf-b8.0/ [M]. It contains p+12C, 14N, 16O, Ca and others up to 150 MeV, including MF3 nonelastic and MF6 double-differential yields. Above 150 MeV: TENDL-2023 (p up to 200 MeV; licence statement to check) or JENDL-5 (p up to 200 MeV for light nuclei) [M].
- **Carlson 1996**, ADNDT 63:93, doi:10.1006/adnd.1996.0010. A compilation of proton reaction and total cross sections up to 1 GeV. It is a paywalled journal, but the factual numbers can be used with citation.
- **EXFOR** (IAEA NDS, free with citation). Web entry points: https://www-nds.iaea.org/exfor/ and the new interface https://nds.iaea.org/exfor/ (x4 search by reaction, e.g. `6-C-12(P,NON),,SIG`, `8-O-16(P,NON),,SIG`, `6-C-12(6-C-12,X),,SIG`). [M] The IAEA also distributes the full EXFOR library as a downloadable bulk and SQLite "X4Pro" package, and on GitHub (IAEA-NDS). **I did not verify exact URL patterns.** Recommended approach: download the bulk X4Pro/EXFOR master once into the provenance cache, hash it, and query it offline. This avoids fragile servlet URLs.
- **Recommended choice:** for protons, fit to measured data (Carlson/EXFOR/ENDF MF3, expected about +/-3-5%). Keep Tripathi-light as an independent check (~5-10%). For He/C/O, use Tripathi with He from the light branch (~5-10% vs data at 100-400 MeV/u) and cross-check with Sihver and Kox. Required evidence: a residual plot against EXFOR points for p+C, p+O, He+C, He+O, C+H, C+C, C+O, O+C, O+O with calibration and evaluation splits recorded.

---

## 2. Proton nuclear secondaries

- **What fast codes do [M, verify against the papers]:**
  - *MCsquare* (Souris, Lee, Sterpin, Med Phys 43 (2016) 2447, doi:10.1118/1.4945036). Nonelastic cross sections are based on ICRU 63. It transports secondary p, d and alpha with sampled energy and angle. Recoils are deposited locally, and neutrons and gammas are discarded (energy escapes).
  - *VMCpro* (Fippel & Soukup, Med Phys 31 (2004) 2263, doi:10.1118/1.1769631). Elastic p-p and p-16O, plus nonelastic p-16O fitted to ICRU 63. Secondary protons are generated with a simplified energy-angle model, and a fixed fraction of energy goes to alphas and recoils (local) and to neutral particles (escape). **I did not retrieve the explicit Fippel-Soukup formulas.** Do not implement it from memory.
  - *gPMC* (Jia et al. 2012 PMB 57:7783; Giantsoudi 2015) follows VMCpro.
  - *FRED* (Schiavi et al. 2017 PMB 62:7482) uses a parameterised nuclear model fitted to FLUKA.
  - Wan Chan Tseung et al. 2015 (Med Phys 42:2967) uses Geant4-derived nonelastic event tables.
- **Energy budget [M]:** Seltzer, NISTIR 5221 (1993), free: https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir5221.pdf. For ~150 MeV p+16O, about 60-65% of the energy lost in nonelastic events goes to charged secondaries (mostly protons, with surprisingly large alpha and recoil fractions). About 35-40% goes to neutrons and gammas, which essentially escape a patient-sized phantom, plus a binding Q of ~ a few %. In water, ~1%/cm of primaries undergo nonelastic events (my calculation: lambda ~ 92-102 cm, so ~20% loss over 26 cm at 200 MeV) [EST].
- **Recommended first implementation [ENG]:**
  1. sigma_ne(E) per element from §1.3.
  2. At each event, sample secondary multiplicities and (E, mu) from preprocessed **ENDF/B-VIII.0 MF6 tables** for O, C, N (Ca and P by A-scaling) at 20-150 MeV. Extend to 250 MeV with TENDL/JENDL or a scaled continuation, labelled APPROX. The MF6 continuum uses the Kalbach-Mann form f(mu; E, E') = a/(2 sinh a) [cosh(a mu) + r sinh(a mu)] [EST, ENDF-6 manual LAW=1 LANG=2]. Preprocess offline into compact CDF tables so the runtime stays numpy-only.
  3. Transport secondary p (and d, t via their stopping power scaled by velocity). Deposit alpha, 3He and A>4 recoils locally, which is APPROX (ranges < ~mm for most). Tally neutron and gamma energy as escaped.
  4. Enforce per-event energy and charge bookkeeping, sampling correlated multiplicities from mean yields. This is APPROX, because uncorrelated sampling preserves only means.
- **Discriminating validation:**
  - The nuclear halo in integral depth dose and lateral profiles. Gottschalk et al. 2014, "Nuclear halo of a 177 MeV proton beam in water", arXiv:1409.1938 (open; contains measured and parameterised data).
  - The fraction of integral depth dose (IDD) missing from narrow detectors.
  - Reference MC: TOPAS QGSP_BIC_HP vs MCsquare.

---

## 3. Projectile fragmentation for He/C/O

### 3.1 Options
| Option | Description | Pros | Cons | Evidence class |
|---|---|---|---|---|
| (a) Semi-empirical partial cross sections | Sihver 1993 (PRC 47:1225); Silberberg-Tsao (ApJS 1973 onward; covers proton targets, scaled); Webber 1990 (PRC 41:566, formula for H targets) | Cheap and transparent. Directly calibratable to measured sigma_cc(Z_F). | Isotope resolution is crude. These systematics target cosmic-ray energies (>~200 MeV/u). | independent theory/empirical |
| (b) Abrasion-ablation | Bowman-Swiatecki-Tsang geometric abrasion (LBL-2908, 1973) plus evaporation. **NUCFRG2**: Wilson et al., NASA TP-3533 (1995), public at https://ntrs.nasa.gov/citations/19960003438 [V downloaded]. **NUCFRG3**: Adamczyk et al. NIM A 678 (2012) 21, and NASA report https://ntrs.nasa.gov/citations/20200005597 [V]. NUCFRG source is not publicly released [M]. **However, Geant4 contains a public NUCFRG2-derived implementation (`G4WilsonAbrasionModel` + `G4WilsonAblationModel`, Geant4 licence) [V file paths].** | Physics-based (Z, A) yields with energy dependence. Can be generated offline into tables. | Poor for very light projectiles (4He) and light fragments (H, He) [M]. | independent theory (but shares lineage with Tripathi/NASA) |
| (c) Tabulated yields from Geant4/FLUKA | Run the reference engine once and table the fragment yields | Realistic correlations | **The same engine cannot then serve as independent validation.** PROTOCOL says "a data table from a reference engine is not a run of that engine", so classify as related model. Licence: FLUKA output use terms. | related model |

**Recommendation [ENG]:** use (a) and (b) to build per-(projectile, target element) tables on an energy grid of 50-450 MeV/u in ~50 MeV/u steps.
- Normalise the sum of fragment-production probabilities to the Tripathi sigma_R.
- Constrain the charge-changing partials to measured data, with a recorded calibration set (e.g. Zeitlin 2007 and Webber 1990 for 12C; Zeitlin 2011 for 16O; Horst 2017/2019 for 4He).
- Evaluate on held-out thick-target data (Haettner 2013, Schall 1996, Tessonnier 2017).
- Use Geant4 QMD/INCLXX/BIC only as second-line independent MC with an inter-model spread.

### 3.2 Kinematics (projectile frame -> lab boost)
- **Fragment velocity [EST, APPROX]:** the fragment keeps approximately the projectile velocity. Mean parallel momentum downshift (Morrissey 1989, PRC 39:460) [M]: <Delta p_par> ~ -8 MeV/c x Delta A, in the projectile frame, with Delta A = A_P - A_F. Result: T_F/A_F ~ T_P/A_P x (1 - small). This is <~1% for Delta A <= 3 at 400 MeV/u.
- **Parallel width (Goldhaber 1974, PLB 53:306) [EST, APPROX]:** sigma_par^2 = sigma0^2 A_F (A_P - A_F)/(A_P - 1), with sigma0 ~ 86-90 MeV/c (Greiner et al. 1975 PRL 35:152) [M]. Equivalent form: sigma0^2 = p_F^2/5 with p_F ~ 200-220 MeV/c. Morrissey gives sigma ~ 87 MeV/c sqrt(Delta A) [M] (the same in the small-Delta A limit).
- **Transverse width [M]:** sigma_perp^2 = sigma_par^2 + sigma_D^2 A_F (A_F - 1)/[A_P (A_P - 1)], with orbital deflection sigma_D ~ 200 MeV/c (Van Bibber et al. 1979 PRL 43:840).
- **Angles:** sample p_perp,x and p_perp,y as Gaussians and set theta ~ p_perp/p_lab.
  - Worked example: 400 MeV/u C has p/A = 951 MeV/c. For B-11, sigma_theta ~ 90/10460 ~ 9 mrad. For protons from C, sigma_theta ~ 90-100 mrad (~5-6 deg). For alphas, ~ 25-30 mrad.
  - This ordering (H broad, He intermediate, Z>=3 narrow) is what Haettner and Gunzert-Marx measured [M].
  - Light fragments (Z=1) also have a non-Gaussian high-p_perp tail from participant/coalescence processes [M]. Add a second, wider component (e.g. 10-20% weight, ~2-3x width) **calibrated on Dudouet 2013 or Haettner angular data**. Mark it APPROX.
- **Target fragments:** low-energy (a few MeV/u) recoils of C, N and O targets deposit locally (ENG/APPROX, range < ~0.1 mm). Target-like protons and alphas with E > ~10 MeV may be transported (optional).
- **Neutrons:** escape. Tally their energy as escaped (V2-NUM requirement).
- **Energy bookkeeping:** use Q-value per channel from mass tables (AME2020, public). Deficit = binding plus neutrons/gammas, recorded explicitly.

---

## 4. Helium and oxygen specifics

**4He channels [EST]:**
- 4He + X -> 3He + n, 3H + p, d + d / d + p + n, p + p + n + n, plus target fragments.
- Mass-changing without charge change is about 3He production, so **sigma(3He-like) ~ sigma_MCC - sigma_CCC**. This gives a discriminating observable from Horst data.
- Measured data [M, numbers in tables in the papers; verify]:
  - **Horst et al. 2017**, PRC 96:024624, doi:10.1103/PhysRevC.96.024624. 4He+12C charge- and mass-changing cross sections, 80-220 MeV/u.
  - **Horst et al. 2019**, PRC 99:014603, doi:10.1103/PhysRevC.99.014603. 4He on H, C, O, Si targets, 70-220 MeV/u.
  - **Rovituso et al. 2017**, PMB 62:1310, doi:10.1088/1361-6560/aa5302. 120 and 200 MeV/u 4He in water and PMMA at HIT: H and He fragment yields, angular and energy spectra.
  - **Aricò et al. 2019** [M; exact citation not verified]: 4He fragment characterisation at HIT with Timepix.
  - **Tessonnier et al. 2017**, PMB 62:3958 (measured 1H/4He/12C/16O depth dose and lateral profiles at HIT) and PMB 62:6579 (MC verification) [M].
  - Norbury et al. 2020, Front Phys 8:565954, reviews He cross-section data gaps [M].

**16O [M]:**
- Zeitlin et al. 2011, PRC 83:034909: 14N, 16O, 20Ne, 24Mg at 290-1000 MeV/u, charge-changing totals and partials on C, CH2, etc.
- Webber 1990, PRC 41:520/533/547/566: 16O on H, He, C at ~ 600 MeV/u.
- Schall et al. 1996, NIM B 117:221 (5 <= Z <= 10, thick water absorbers, attenuation).
- Sihver et al. 1998, Jpn J Med Phys 18:1 (HIMAC depth dose for C, O, Ne).
- Tessonnier 2017 (HIT O depth dose).
- FOOT collaboration 16O+C data at 200-400 MeV/u (2022+, check).
- Böhlen 2010 (PMB 55:5833) is a **carbon** benchmark (Haettner data), not oxygen.

**Proposed oxygen supported domain [ENG]:**
- 16O in water and soft tissue (H, C, N, O, P, Ca), 100-430 MeV/u.
- Claimed capabilities: primary attenuation; charge-changing partial yields (Z=1-7) and their depth build-up; IDD including the distal tail; dose-averaged LET of the mixed field with stated model uncertainty.
- Not claimed unless validated: isotope-resolved fragment spectra, double-differential fragment spectra, and claims below 100 MeV/u entrance energy (slow fragments are still transported).
- Unsupported combinations must fail closed.

---

## 5. Carbon validation data [M unless noted; verify access]

| Dataset | Content | Numeric availability |
|---|---|---|
| Haettner, Iwase, Krämer, Kraft, Schardt 2013 PMB 58:8265, doi:10.1088/0031-9155/58/23/8265 | 200 and 400 MeV/u 12C in water: Z-resolved build-up/attenuation vs depth, angular distributions, energy spectra | Figures. Some tables in Haettner's 2006 KTH MSc thesis and in Haettner et al. 2006 Radiat Prot Dosim 122:485. Widely digitised (Böhlen 2010, G4-Med). Digitising is acceptable with a recorded digitisation uncertainty. |
| Schall et al. 1996 NIM B 117:221 | Charge-changing attenuation in water for Z=5-10 | Table in paper |
| Toshito et al. 2007 PRC 75:054606 | 12C on water and polycarbonate, 200-400 MeV/u, total and partial sigma_cc | Tables |
| Zeitlin et al. 2007 PRC 76:014911 | 12C at 290 and 400 MeV/u on elemental targets, total and partial sigma_cc | Tables |
| Webber et al. 1990 PRC 41:520-566 | Charge- and mass-changing and isotopic cross sections on H, He, C | Tables |
| Golovkov et al. 1997 (Adv. Hadrontherapy, Excerpta Medica) | 12C fragmentation in water | Paper |
| Matsufuji et al. 2003 PMB 48:1605; 2005 PMB 50:3393 | Fragment yields in water (HIMAC), spatial distribution | Figures and tables |
| Gunzert-Marx et al. 2008 NJP 10:075003 (open access) | 200 MeV/u C in water: secondary fragments and dose contributions | Open-access figures and tables |
| Dudouet et al. 2013 PRC 88:024606; Divay et al. 2017 PRC 95:044602 | 95 and 50 MeV/u 12C thin targets (H, C, O, Al, Ti): double-differential cross sections | **Public numeric database at http://hadrontherapy-data.in2p3.fr** [M, verify] |
| Kurz/Tessonnier/HIT and GSI measured IDDs | Bragg curves with tails | Figures. HIT base data are not public [M]. |

Licensing: numeric experimental results are facts. Use them with citation and record the provenance (source URL, retrieval time, hash, digitisation method). Do not redistribute publisher PDFs.

---

## 6. Validation design (discriminating observables)

| Observable | What it falsifies | Expected magnitude | Suggested tolerance for a semi-empirical model [ENG, freeze before qualification] |
|---|---|---|---|
| Z=P survival vs depth (charge-changing; Haettner counts all Z=6, including 11C and 10C) | sigma_R/sigma_cc scale and energy dependence | 12C: ~65-70% at the BP for 200 MeV/u (8.6 cm). **~25-35% for 400 MeV/u (27.5 cm)**. My Tripathi calculation: 65% and 26% (mass-changing; charge-changing is slightly higher). 16O 250 MeV/u: ~50%. 4He 200 MeV/u: ~58%. Proton 200 MeV: ~78%. | 3% absolute in surviving fraction |
| Per-Z fragment build-up vs depth (Z=1..P-1) | Partial cross sections and secondary attenuation | H and He dominate by number (He about 0.2-0.4 per primary near the BP at 400 MeV/u [M]) | 20-30% for Z=1,2; 30-50% for Z=3-5 (inter-code spread levels) |
| Distal tail dose 1-10 cm beyond the BP | Fragment yields x ranges | 12C pristine peak at 270-400 MeV/u: tail ~5-10% of peak (~15-25% of entrance dose), decreasing gradually [M]. SOBP tail ~10-15% of plateau. He: ~1% of peak. O: larger than C. | 15-20% relative or 1% of peak absolute |
| Fragment angular distributions (Z=1,2) at depth | Goldhaber widths and wide component | H: FWHM of several degrees; He: ~2-3 deg [M] | 20-30% on width |
| Lateral halo and low-dose envelope | Fragment angles plus proton secondaries | Few % of the central dose at 1-3 cm off-axis | Shape check against TOPAS spread |
| Energy spectra of fragments at fixed depth | Velocity shift and widths | Peak near the primary velocity | Position within ~5% |
| Thin-target sigma_cc partials vs E | Calibration and held-out separation | — | 10-15% |

Published inter-code and model-data disagreement [M]:
- Böhlen 2010 (PMB 55:5833): FLUKA and Geant4 vs Haettner give yield deviations of ~10-30%, larger for angular tails.
- Dudouet 2014 (PRC 89:064615): Geant4 models vs 95 MeV/u double-differential data differ by up to factor ~2 at large angles, with integrated yields within ~20-50%.
- Bolst et al. 2017 (NIM A 869:68) and Arce et al. 2021 (Med Phys 48:19, G4-Med): QMD is best overall for Haettner build-up within ~10-20%.

So a semi-empirical model cannot credibly claim better than ~10-20% on fragment yields. It must achieve ~2-5% on primary attenuation and on dose in the plateau and peak.

---

## 7. TOPAS / Geant4 11.4 as reference [M]
- **TOPAS default modular list:** `g4em-standard_opt4`, `g4h-phy_QGSP_BIC_HP`, `g4decay`, `g4ion-binarycascade`, `g4h-elastic_HP`, `g4stopping`.
- **For C and O, run at least two ion models:** `g4ion-qmd` (G4IonQMDPhysics) and BIC. Optionally add `g4ion-inclxx` (INCL++, light-ion projectiles A<=18). Use the spread as reference-model uncertainty.
- **Model performance:** literature (Bolst 2017, Arce 2021, Böhlen 2010, Dudouet 2014) finds QMD best for fragment yields and angular distributions at 200-400 MeV/u. BIC tends to misestimate light-fragment yields and angular widths. All Geant4 models show tens-% discrepancies for large-angle H/He.
- **Cross sections:** the Geant4 ion inelastic cross section defaults to Glauber-Gribov, so it is independent of Tripathi.
- **Logs and cuts:** inspect the physics list actually printed in logs, the production cuts (e.g. 0.05-0.1 mm) and the step limits. Geant4 11.x changed some de-excitation and ion defaults, so record the version and data sets (G4ENSDFSTATE, G4PARTICLEXS, G4TENDL).
- **Other engines:** MCsquare is proton-only. FRED's ion capability needs checking. So TOPAS is the main independent ion MC; measured data must carry the release claims.

---

## 8. Open items for the orchestrator
1. Retrieve and verify Fippel-Soukup 2004, Souris 2016 and Seltzer NISTIR 5221 before implementing proton secondaries. Prefer the ENDF MF6 route.
2. Run a Tripathi residual analysis against EXFOR. Watch the p+C overestimate (5-15%?).
3. Verify access to hadrontherapy-data.in2p3.fr and the EXFOR bulk download. Plan to digitise Haettner 2013, with digitisation uncertainty.
4. Decide the calibration/evaluation split: thin-target partials for calibration, thick-target Haettner/Schall/Tessonnier/Rovituso for evaluation.
