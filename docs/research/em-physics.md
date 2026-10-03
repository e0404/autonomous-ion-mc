# EM condensed-history physics for therapeutic ions: decision report

Scope: p, He, C, O ions from 1 to about 500 MeV/u in water and tissue. Runtime dependencies are numpy and warp only.
Labels: **[EST]** established physics; **[APPROX]** approximation; **[ENG]** engineering choice; **[NV]** not verified in this session (it comes from memory or secondary sources, so check it before you rely on it).
Notation: T = kinetic energy; T/A or MeV/u = energy per nucleon; z = projectile charge; S = mass electronic stopping power in MeV cm2/g.

*[Lead redaction note, 2026-10-03: exact PSTAR/ASTAR values that the agent quoted as examples were removed from this note; NIST SRD 124 data are copyrighted, all rights reserved, and the project commits no NIST values (decision 0038). Earlier revisions of this file in the repository history contain the quoted values; see the intervention record IR-20261003-192301-6CF85A.]*

## 0. Recommendations at a glance
| Topic | Recommendation | Label |
|---|---|---|
| Stopping (p, He) | Build our own Bethe tables, with shell, Barkas, Bloch and Mott corrections and I_water = 78 eV, inside the package. Validate them against NIST PSTAR/ASTAR (I = 75 eV, auto-downloaded) and against the ICRU 90 water arrays embedded in Geant4 source. | ENG |
| Stopping (C, O) | Use the same Bethe core with the actual z. Use Barkas (z^3) and Bloch (z^4) terms, not plain z^2 scaling. Use effective charge only below about 10 MeV/u. Validate against ICRU 73 (+ errata) tables shipped in G4EMLOW, TOPAS runs and ATIMA/catima. | ENG/EST |
| Energy loss along a step | Use the range-table inversion E1 = R^-1(R(E0) - rho*s). | ENG (exact in CSDA) |
| Straggling | Gaussian with Bohr variance (relativistic form). Switch to a Gamma distribution when mean/sigma < 2, as G4IonFluctuations does. Yang correction optional. | APPROX |
| MCS | Gaussian core from a step-size-independent scattering power (differential Highland or Molière-fitted). Apply it with a random hinge. An optional single-scattering tail is a fidelity switch. | APPROX/ENG |
| Steps | Limit to min(voxel boundary, energy-loss fraction 1-2%, range fraction, s_max). Below an end-of-range cutoff of about 1-2 MeV/u, deposit the remaining energy locally. | ENG |
| delta electrons | Deposit locally; there is no electron transport. Use unrestricted S. | APPROX |
| LET | Compute LET_t = sum(l*S)/sum(l) and LET_d = sum(l*S^2)/sum(l*S), with S = unrestricted electronic stopping in water at the mean step energy, summed over all ions including fragments. Also report the epsilon-weighted variant. | ENG |

## 1. Electronic stopping power

### 1a. Free data sources
| Source | Coverage | Access (tested 2026-10-03 unless marked) | Terms |
|---|---|---|---|
| **NIST PSTAR** (SRD 124, doi:10.18434/T4NC7P) | Protons, 74 materials, 1 keV-10 GeV, 133 default energies. Equivalent to ICRU 49, which uses I_water = 75 eV. | **Verified working:** `POST https://physics.nist.gov/cgi-bin/Star/apdata.pl` (urlencoded or multipart) with fields `prog=PSTAR`, `matno=276` (liquid water; 104 = dry air, 119 = compact bone ICRU, 223 = PMMA, etc. from the `<select name="matno">` list on PSTAR.html), `ShowDefault=on`, `NumofEnergies=0`, `character=space`, `electronic=on`, `nuclear=on`, `total=on`, `csda=on`, `project=on`, `detour=on`. The response is plain text: 8 header lines, then columns T [MeV], S_el, S_nuc, S_tot [MeV cm2/g], CSDA and projected range [g/cm2], and the detour factor. GET is rejected ("Request Invalid"). The HTML form endpoint is `/cgi-bin/Star/ap_table.pl`, which requires `GraphType=None` and accepts custom energies in the `Energies` textarea. Python stdlib `urllib` is enough to download. | NIST SRD: **copyright claimed** under the Standard Reference Data Act (15 USC 290e). Attribution and the notice "Copyright protection on this compilation of data has been secured by the Secretary of the U.S. Department of Commerce on behalf of the United States" are required (https://www.nist.gov/open/license). **Auto-download to a cache is fine; do not commit the tables to Git.** |
| **NIST ASTAR** | 4He, same materials, 1 keV-1 GeV **total** kinetic energy (that is, up to 250 MeV/u), 122 energies | Same POST with `prog=ASTAR`. Verified (example values redacted: NIST SRD 124 data are use-only and are not reproduced in this repository). | Same as PSTAR |
| **ICRU 90 water/air/graphite (p, alpha)** | Proton arrays 1 keV-10 GeV, alpha arrays 1 keV-1 GeV, I_water = 78 eV | Hard-coded in Geant4 `source/materials/src/G4ICRU90StoppingData.cc` (arrays `T0_proton`/`e1_proton`, where index 1 is G4_WATER). For example, at 200 MeV the ICRU 90 value is S = 4.470 MeV cm2/g (the corresponding PSTAR value is redacted: NIST SRD 124 data are use-only). | Geant4 Software License (permissive, attribution required). The ICRU 90 report itself is copyrighted, so the provenance of the numbers needs a recorded legal judgement [NV]. |
| **Geant4 G4EMLOW 8.8** `ion_stopping_data/` | ICRU 73 (+ errata) tables: projectile Z = 3-80 for 31 named NIST materials (including G4_WATER, G4_MUSCLE_STRIATED_ICRU, G4_BONE_COMPACT_ICRU, G4_ADIPOSE_TISSUE_ICRP, G4_LUCITE, G4_AIR) and element targets. There is also an `icru90/` variant for Z <= 18 in G4_WATER, G4_AIR and G4_GRAPHITE. | URL https://cern.ch/geant4-data/datasets/G4EMLOW.8.8.tar.gz (Geant4 11.4). Paths, verified from the reader source `G4IonICRU73Data.cc`/`G4IonStoppingData.cc`: `G4EMLOW8.8/ion_stopping_data/icru73/z{Zion}_{G4_MATNAME}.dat`, `.../icru90/z{Zion}_G4_WATER.dat`, and element files `z{Zion}_{Ztarget}.dat`. The format is G4PhysicsFreeVector ASCII (edge min/max/n, then energy-value pairs). Energy is per nucleon in MeV, and dE/dx is in MeV cm2/mg (the reader multiplies by 1000 x density). **The tarball listing was not completed in this session** (the download of about 300 MB was still running), so confirm the file names and the energy span of each file [NV]. Note that current Geant4 uses these tables only between 0.025 and 2.5 MeV/u (G4IonICRU73Data fEmin/fEmax) and uses Bethe + Lindhard-Sorensen above that. | Geant4 Software License for the dataset [NV: confirm a LICENSE file in the tarball]. This is a data table from a reference engine, **not** an engine run (PROTOCOL.md). |
| ICRU 73 report / ICRU 90 report | C and O in water, 0.025-1000 MeV/u | Purchase only | Copyrighted; do not redistribute. |
| SRIM | All ions | Closed-source Windows executable; cannot be auto-downloaded or run under numpy+warp | Output is usable for comparison; the redistribution terms are unclear [NV]. **Not recommended.** |
| libdEdx (PSTAR/ASTAR/MSTAR/ICRU73 tables), ATIMA/catima/pycatima (GSI, Lindhard-Sorensen + effective charge) | All ions | Public repositories | Copyleft licences likely (GPL/AGPL) [NV]. Use them only as offline independent-theory references and never vendor them. |

Lineage caution [EST]: PSTAR, ASTAR, ICRU 73 and ICRU 90 all rest on Bethe theory above about 1 MeV/u. Agreement with them is *ion-specific tabulated* evidence, not *measured* evidence. ICRU 73 water was originally computed with I = 67.2 eV and later corrected by errata (Sigmund, Schinner and Paul, J. ICRU 5(1), 2009) [NV]. Record which version the G4EMLOW file reflects.

### 1b. Analytical model (use it to generate the tables in the package) [EST, with APPROX corrections]
S_el/rho = K z^2 (Z/A) beta^-2 [ L0 + z L1 + z^2 L2 ] (MeV cm2/g)

- K = 4 pi N_A r_e^2 m_e c^2 = 0.307075 MeV cm2/mol. Water: Z/A = 0.555087 mol/g. Mixtures follow Bragg additivity of Z/A; I comes from ln I = sum w_i (Z/A)_i ln I_i / (Z/A).
- L0 = 1/2 ln(2 m_e c^2 beta^2 gamma^2 T_max / I^2) - beta^2 - delta/2 - C/Z
- T_max = 2 m_e c^2 beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2), with m_e c^2 = 0.51099895 MeV.
- Density effect delta (Sternheimer): water uses x0 = 0.2400, x1 = 2.8004, C-bar = 3.5017, a = 0.09116, m = 3.4773 [NV: PDG/Sternheimer 1984 table]. For T/A up to 500 MeV/u, x = log10(beta*gamma) is at most 0.09 < x0, so **delta = 0 in water across the whole therapeutic domain**. Keep the term for generality.
- Shell correction C/Z: it matters below about 20 MeV/u (several % at about 2 MeV). Options are the Bichsel parameterisation (ICRU 49) or Geant4's K/L-shell model (`G4EmCorrections::ShellCorrection`). It affects only the last few mm of range.
- Barkas term z L1: use Ashley-Ritchie-Brandt / Lindhard as tabulated in ICRU 49, or the Geant4 `BarkasCorrection` implementation [APPROX].
- Bloch term z^2 L2 = -y^2 sum_{n>=1} 1/[n(n^2+y^2)], with y = z alpha/beta and alpha = 1/137.036 [EST].
- Mott term (Geant4 form) = 0.5 * pi alpha beta z [APPROX]. For C/O above about 100 MeV/u, the Lindhard-Sorensen correction (Phys. Rev. A 53, 2443, 1996) is the complete replacement for Bloch + Mott. Expect differences below 1% [NV].
- Effective charge below about 10 MeV/u: Pierce-Blann/Barkas z_eff = z [1 - exp(-125 beta z^-2/3)]. For carbon at 10 MeV/u, z_eff/z = 0.996; at 1 MeV/u it is about 0.83. Ziegler/Geant4 (`G4ionEffectiveCharge.cc`) has separate He and heavy-ion fits [APPROX]. This affects only the last ~0.4 mm of carbon range (the CSDA range of 10 MeV/u carbon in water is about 0.4 mm), so **its effect on range is below dose-voxel resolution**.
- Below about 1 MeV/u, Bethe breaks down (L can become negative). Splice to the PSTAR/ASTAR low-energy shape, or simply deposit locally at the cutoff [ENG]. A 1 MeV proton has a range of only 24 um in water.
- Ion stopping by "z^2-scaled proton table" is an approximation. It misses the (z-1) L1 and (z^2-1) L2 differences, which are about 1-2% at tens of MeV/u for C/O [NV magnitude]. Use the explicit formula with the real z.

### 1c. I-value of water [EST]
ICRU 90 (2016) gives I_water = 78 +/- 2 eV; ICRU 49 and PSTAR use 75 eV. This session integrated the electronic stopping from 1 MeV upward (log-log interpolation) using the ICRU 90 arrays from Geant4 against PSTAR:

*(Per-energy ICRU 90 − PSTAR range differences computed by the agent were redacted from this note on 2026-10-03; the committed aggregate comparison in `validation/results/stopping/` and decision 0038 carry the project's own evaluation. The agent's finding was that the 75 → 78 eV change shifts the 200 MeV proton range by about +1 mm, of order 0.4 %.)*

A ±2 eV uncertainty in I corresponds to about ±0.3% of range. Recommendation: default to 78 eV for water and make I a configurable, provenance-recorded parameter. Validate twice: 75 eV against PSTAR/ASTAR (this should reproduce them to about 0.1-0.2%), and 78 eV against the ICRU 90 arrays. This is the discriminating test of whether the implementation is correct or only calibrated.

## 2. Energy-loss straggling
- Vavilov parameter kappa = xi/T_max, with xi = (K/2)(Z/A) rho x z^2/beta^2 [MeV]. kappa > ~10 gives Gaussian (Bohr) straggling; kappa < 0.01 gives Landau [EST]. Example: a 200 MeV proton crossing 1 mm of water has xi = 0.027 MeV, T_max = 0.48 MeV and kappa = 0.055, so a single step is Vavilov-like. Carbon has a z^2 = 36 times larger kappa and is near-Gaussian.
- Bohr variance with the relativistic factor and no delta production: sigma^2 = 2 pi r_e^2 m_e c^2 n_el z^2 x (T_max/beta^2)(1 - beta^2/2) ≈ K m_e c^2 rho (Z/A) z^2 gamma^2 (1 - beta^2/2) x. Here K m_e c^2 = 0.15692 MeV^2 cm2/mol, so water gives ≈ 0.0871 z^2 gamma^2 (1 - beta^2/2) MeV^2/cm [EST].
- Why a Gaussian per step is acceptable for dose [APPROX]: delta rays are deposited locally, so the energy loss in each step is deposited in the same voxel. Range straggling is a sum over many steps, so by the central limit theorem only the mean and variance per step need to be right. The Landau tail of a single step does not reach the dose distribution.
- Practical model: copy the logic of G4IonFluctuations (verified from source) [ENG]. If mean/sigma >= 2, sample a Gaussian truncated to (0, 2*mean). If 0.1 < mean/sigma < 2, sample a Gamma distribution with n = (mean/sigma)^2. Otherwise sample uniformly on (0, 2*mean). There is also an optional Yang et al. (NIM B 61, 149, 1991) charge-state and low-velocity factor. Geant4 uses UniversalFluctuation (Urban) above T > 10 MeV x z x (M/m_p), but Urban is not needed when deltas are local.
- Validation: range straggling sigma_R ≈ 0.012 R^0.935 cm for protons (Bortfeld 1997), which is about 1.0-1.1% of R. For ions, sigma_R/R scales roughly as 1/sqrt(A) at equal range [APPROX]. Test that the dose is independent of step size when the step is refined.

## 3. Multiple Coulomb scattering
| Model | Formula / notes | Use |
|---|---|---|
| Highland/PDG (Lynch & Dahl, NIM B 58, 6 (1991), doi:10.1016/0168-583X(91)95671-Y) | theta0 = (13.6 MeV/(beta c p)) z sqrt(x/X0) [1 + 0.038 ln(x z^2/(X0 beta^2))]; 11% accuracy for 1e-3 < x/X0 < 100. Water X0 = 36.08 g/cm2. Gottschalk et al. 1993 use 14.1 MeV and (1 + (1/9) log10(x/X0)) (NIM B 74, 467, doi:10.1016/0168-583X(93)95944-Z), checked against measured 158.6 MeV proton data. | Validation only; see the pitfall below |
| Differential Highland / differential Molière (Gottschalk, Med Phys 37, 352 (2010), arXiv:0908.1413) | Scattering power T = f(pv, p1v1)(E_s/pv)^2/X_S, with E_s = 15.0 MeV and X_S = "scattering length" (water about 46.9 g/cm2 [NV]). f_dH = 0.5244 + 0.1975 lg(1 - (pv/p1v1)^2) + 0.2320 lg(pv) - 0.0098 lg(pv) lg(1 - (pv/p1v1)^2), with pv in MeV [NV: coefficients from memory; check against arXiv:0908.1413]. It is local and step-additive. | **Recommended core** |
| Molière (Bethe form) | chi_c^2 = 0.157 Z(Z+1) z^2 x/(A (pv)^2) (MeV^2, x in g/cm2). B is solved from B - ln B = ln(chi_c^2/(1.167 chi_a^2)), and theta_1/e = chi_c sqrt(B - 1.2). Includes the single-scattering tail. | Reference/validation; a "tail" fidelity option |
| Goudsmit-Saunderson | Exact angular distribution for a given step length, built from screened Rutherford/Mott cross sections. | Overkill for ions; Geant4 uses it for electrons only |

- **Pitfall [EST]:** the log term in Highland is not additive. Applying Highland per step makes the total angle depend on step size, so smaller steps give too little scattering. Use a step-independent scattering power, or compute increments Delta(theta^2) = theta0^2(s0 + s) - theta0^2(s0) along the integrated path, the generalised Highland form with ∫ dx/(pv)^2 [EST/ENG]. Make this an explicit V2-NUM falsification test: compare 0.1, 1 and 5 mm steps.
- Scaling: theta0 is proportional to z/(pv). At equal velocity, pv is proportional to A, so theta0 is proportional to z/A: about 1/2 of the proton value for He, C and O at the same MeV/u. At equal water range, the lateral spread at end of range is about 2% of R for protons, about 1% for He and about 0.6% for C [APPROX; NV magnitudes].
- Lateral displacement within a step [EST]: the Fermi-Eyges second moments for one step are y = s theta0 (z1/sqrt(12) + z2/2) and theta = z2 theta0 (PDG). A **random hinge** (move u*s, deflect, move (1-u)*s, with u uniform) reproduces <y^2> = s^2 theta0^2/3 and <y theta> = s theta0^2/2 exactly, so no explicit correlation sampling is needed. Hinge positions must respect voxel boundaries: truncate at the boundary, then re-hinge.
- Single-scattering tail: at the 1e-3-1e-2 level it adds to the low-dose envelope and field-size factors. For protons the halo is dominated by nuclear secondaries (Gottschalk et al., PMB 60, 5627, 2015 [NV doi]). Offer a Gaussian core as the default and a "Molière/Rutherford tail" fidelity flag (Poisson number of single scatters above theta_cut, with screened-Rutherford sampling).
- What other codes use: TOPAS/Geant4 opt4 uses WentzelVI + single Coulomb scattering for protons and Urban MSC for GenericIon [NV]. MCsquare (Souris et al., Med Phys 43, 1700, 2016) uses a Gaussian based on the Rossi formula with a correction plus random hinge [NV]. FRED (Schiavi et al., PMB 62, 7482, 2017) uses a Gaussian core plus single-scattering tail [NV]. gPMC (Jia et al., PMB 57, 7783, 2012) [NV].

## 4. Condensed-history step control [ENG unless noted]
- Energy loss: tabulate CSDA R_m(E) (g/cm2) and its inverse per material on a log grid with about 100-200 points per decade (in float64 when building the tables). The mean energy after a step is E1 = R_m^-1(R_m(E0) - rho*s). In a homogeneous medium this is exact in CSDA and independent of step size [EST], which removes the S(E)*s bias at the Bragg peak. Deposit (E0 - E1) plus the fluctuation, and evaluate S for LET at the mean step energy.
- Limits: s = min(d_boundary, s_eloss, s_range, s_max), where
  - s_eloss caps the energy-loss fraction at 1% (accurate) to 5% (fast) of T. It matters mainly for MCS and straggling, not for the mean loss.
  - s_range follows the Geant4 step function. For R > rho_f, s = alpha R + rho_f (1 - alpha)(2 - rho_f/R), with defaults alpha = 0.2, rho_f = 0.1 mm for protons and alpha = 0.1, rho_f = 0.02-0.05 mm for ions [NV defaults].
  - s_max is about 1 voxel.
  Expose a "fidelity" preset that fixes these values and record the effective values.
- Voxels: use DDA/Siddon voxel walking with a step truncated at each boundary, and recompute material quantities on crossing. Alternatively, a water-equivalent path-length step with stopping-power-ratio scaling (MCsquare-like) is faster but gives an approximate MCS-energy coupling. Score on separate scoring grids by splitting track segments at score-grid boundaries; otherwise LET and dose depend on alignment, which V2-NUM requires you to test.
- Class I vs class II: with local delta deposition, a pure continuous loss with unrestricted S plus straggling is enough. Class II (explicit deltas above a cut) only adds value if electrons are ever transported. Keep the interface open for it.
- Cutoff: when T/A < E_cut (about 1-2 MeV/u, residual range ≤ 0.07 mm for protons), deposit the rest of the energy locally at the current point (or spread it over the residual range). Expose E_cut as a parameter.

## 5. delta electrons: local deposition [APPROX]
- Kinematics: T_max = 0.48 MeV for 200 MeV protons, 0.62 MeV for 250 MeV protons and 1.16 MeV for 430 MeV/u carbon. The corresponding CSDA electron ranges in water are about 1.7 mm, 2.3 mm and 5 mm (ESTAR [NV values]).
- Energy fraction carried above a threshold T_c: f ≈ [ln(T_max/T_c) - beta^2(1 - T_c/T_max)] / [ln(2 m_e c^2 beta^2 gamma^2 T_max/I^2) - 2 beta^2]. For 200 MeV protons, f ≈ 8% above 100 keV (electron range about 0.14 mm) and ≈ 2% above 300 keV (about 0.8 mm). These electrons are forward-peaked (cos theta ≈ sqrt(T(T_max + 2mc^2)/(T_max(T + 2mc^2)))).
- Consequences: inside a uniform medium, charged-particle equilibrium cancels the transport error. Residual errors appear in the first few mm (entrance build-up of order 1% or less) and at large density steps (air cavities, lung and bone interfaces), with errors expected at the about 1% level within about 1-2 mm of interfaces [APPROX; estimate, not measured]. Make this a documented domain limitation and check it with a TOPAS run that compares production cuts (for example, 1 mm vs 1 um) at an air/water interface.
- LET implication: with local deposition, use **unrestricted** electronic stopping (LET_inf) for the LET definition, consistent with the deposited energy. Restricted LET_Delta would need a delta cut, which this code does not have. State this in the metadata.

## 6. LET estimators [ENG, grounded in the literature]
Literature: Cortés-Giraldo & Carabe, PMB 60, 2645 (2015), doi:10.1088/0031-9155/60/7/2645; Granville & Sawakuchi, PMB 60, N283 (2015), doi:10.1088/0031-9155/60/14/N283; Guan et al., Med Phys 42, 6234 (2015) [NV]. These studies show that estimators using epsilon_i/l_i are step-size dependent and biased by fluctuations and the treatment of secondary electrons. Estimators that use tabulated S at the mean step energy are stable.
- Track-averaged: LET_t(v) = sum_i l_i S_i / sum_i l_i
- Dose-averaged: LET_d(v) = sum_i l_i S_i^2 / sum_i l_i S_i. This is the expectation of the epsilon-weighted form when epsilon_i ≈ S_i l_i, and it has lower variance. Also offer the variant LET_d^eps = sum epsilon_i S_i / sum epsilon_i for cross-checking.
- S_i: unrestricted electronic stopping power of the transported ion at the mean step energy. Default to **water** ("LET_d in water"), with a material option.
- Mixed fields: sum over primaries **and** all transported secondary ions (p, d, t, 3He, alpha, Li-O fragments). Provide species-resolved partial sums (numerator and denominator per species) so users can recombine them. Nuclear-recoil or locally deposited fragment energy that is not transported needs an explicit rule: exclude it from LET and record that it was excluded [ENG]. For carbon, the fragments dominate LET_t beyond the peak, and LET_d there is set by heavy fragments.
- Uncertainty: LET is a ratio estimator. Use batch means on the numerator and denominator (delta method), and mark it undefined where the denominator is zero.

## 7. Validation targets
| Quantity | Reference | Evidence class |
|---|---|---|
| CSDA range p | PSTAR at 100/150/200/250 MeV (values redacted: NIST SRD 124 data are use-only; the table is downloaded by the user) | ion-specific tabulated |
| CSDA range He | ASTAR at 400/600/800 MeV (100/150/200 MeV/u) (values redacted, as above) | ion-specific tabulated |
| C, O range and S | ICRU 73 (errata) via G4EMLOW; ATIMA/catima; TOPAS runs | tabulated / independent MC |
| I-value sensitivity | 75 vs 78 eV gives +0.4% range (section 1c) | self-consistency + tabulated |
| Bragg-Kleeman | R = alpha E^p, alpha = 0.0022 cm MeV^-p, p = 1.77 (water, protons; Bortfeld, Med Phys 24, 2024 (1997), doi:10.1118/1.598116) | related model (about 1-2%) |
| Bragg curve shape | Bortfeld 1997 analytic curve (includes nuclear terms) | related model only |
| Range straggling | sigma_R ≈ 0.012 R^0.935 cm (Bortfeld) | related model |
| MCS | Gottschalk 1993 measured theta0 tables (158.6 MeV p, many materials and thicknesses: digital numbers in the paper); Highland/Molière theory | measured / theory |
| Measured ion depth dose and lateral profiles | Tessonnier et al., PMB 62, 6579 (2017) (HIT p/He/C/O: tabulated range, peak width, fall-off and lateral Gaussian parameters; full curves are not verified as public) [NV doi]; Schwaab et al., PMB 56, 7493 (2011) (lateral tails of p/C, HIT) [NV]; Haettner et al., PMB 58, 8265 (2013) (12C 200/400 MeV/u fragment yields and angular distributions in water, tabulated) [NV doi] | measured |

No clean public digital measured He/C/O Bragg-curve dataset was confirmed in this session. If a search of Zenodo, geant-val.cern.ch and PTCOG/HIT supplements, within the PROTOCOL budget, fails to find one, file a reference request and run TOPAS/FRED/MCsquare through the reference MCP.

## Uncertainties and discriminating evidence
1. **ICRU 73 version and span in G4EMLOW:** read `ion_stopping_data/icru73/z6_G4_WATER.dat` and compare it with Bethe(78 eV) + Barkas + Bloch at 10/100/400 MeV/u. If the file reflects the pre-errata data (I = 67.2 eV), expect a few % at low energy.
2. **Barkas/Bloch implementation for C/O:** compare the computed carbon CSDA range at 400 MeV/u (about 27.5 cm in water [NV]) with ICRU 73, catima and TOPAS. A z^2-scaling bug shows up as a difference of 1 mm or more.
3. **MCS step dependence:** compute lateral sigma at the Bragg peak with 0.1/1/5 mm steps. Per-step Highland fails; the differential form must be flat. Then compare against TOPAS opt4 and Gottschalk's measured theta0.
4. **delta local deposition:** run TOPAS with high vs low production cuts in water/air/bone slabs and quantify the interface error domain.
5. **LET estimator choice:** refine the step and grid, and compare the l*S^2 and epsilon*S estimators. The difference should go to zero within statistics in the validated domain.
6. **Legal questions:** the attribution for ICRU 90 numbers redistributed via Geant4 source, the G4EMLOW license file, and the libdEdx/catima licences are all [NV]. Record them as provenance decisions before ingesting any of these.

Sources: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html ; https://physics.nist.gov/cgi-bin/Star/apdata.pl ; https://www.nist.gov/open/license ; https://catalog.data.gov/dataset/nist-stopping-power-range-tables-for-electrons-protons-and-helium-ions-srd-124 ; https://geant4.web.cern.ch/download/11.4.0.html ; Geant4 source (github.com/Geant4/geant4, master, read 2026-10-03): G4IonICRU73Data.cc, G4IonStoppingData.cc, G4ICRU90StoppingData.cc, G4IonFluctuations.cc, G4EmCorrections.cc, G4ionEffectiveCharge.cc ; https://arxiv.org/abs/0908.1413 ; https://heibib.ub.uni-heidelberg.de/search/Record/1571674160
