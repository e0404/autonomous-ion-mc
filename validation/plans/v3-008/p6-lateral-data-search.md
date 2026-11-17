# P6 part 2 — measured ion lateral-spreading data search (V3-008), research note 2026-11-01 (not evidence)

Labels: FACT = fetched and verified this session; LIT = from secondary/search snippets or memory; OPEN = not checked.
Scope: independent evidence of the CHARGE DEPENDENCE of lateral spreading for He/C/O in water, 50–430 MeV/u.
Nothing written to any worktree. Scratch downloads (outside the repository) under
.../scratchpad/v3-008/p6dl/{parodi,kanematsu,gottschalk,tessonnier,minibeam,dpg}/ (sha256 below).
Stopped at the coordinator's turn limit: items marked OPEN were not searched.

## 0. Key framing finding (FACT)
Gottschalk, "Techniques of Proton Radiotherapy: Transport Theory", arXiv:1204.4470 (2.44 MB, sha256 3ddf35b7…),
§5.3 "Experimental Tests" of heavy-ion beam spreading: "A number of transverse scans of heavy ion beams vs. depth are
available. Unfortunately they tend to be motivated entirely by clinical applications (beam commissioning), thus
dominated by ion-optical effects or beam contamination [Schwaab 2011] and not too useful as tests of the basic theory.
After moderate effort we have found only the following three": Preston & Koehler (protons only), "Phillips" (p, He, C;
provenance doubtful, see 2), Wong et al. 1990 (238U in Cu; angular, not water). Conclusion (p. ~43): "There are few
experimental tests of beam spreading in matter." This independent expert survey supports the evidence-gap verdict.

## 1. Candidates

| # | Candidate | What / ions / energies / medium | Uncertainties | Numbers tabulated? | URL / licence / size | Lineage | Status |
|---|---|---|---|---|---|---|---|
| 1 | Parodi, Mairani, Sommerer 2013, J. Radiat. Res. 54 Suppl 1 i91 (doi 10.1093/jrr/rrt051, PMC3700515) | FLUKA-based double-Gauss TPS parametrisation for p and C at HIT (20 energies each), water target 12 cm upstream of isocentre, incl. confidential nozzle geometry; validated against Schwaab 2011 measurements (Fig. 1, p 157.53 MeV, C 299.94 MeV/u, ~1.5 cm and ~16.5 cm). No He, no O | none stated for the parametrisation | NO: figures only (Figs 1–6), no tables | Europe PMC full-text XML https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3700515/fullTextXML (39 905 B, sha256 e9cb45c7…); CC BY 3.0 | MODEL-DERIVED (FLUKA), includes beamline → not independent | FACT |
| 2 | Schwaab, Brons, Fieres, Parodi 2011, PMB 56:7813 (doi 10.1088/0031-9155/56/24/009) | Measured lateral profiles (line scans, pinpoint IC block) of p and C at HIT, several energies/depths, water and air; low-dose envelope focus | OPEN | Not OA (Europe PMC: OA N, no suppl.). Tables: OPEN (not accessible); Parodi 2013 and arXiv:2407.09851 use them as plotted data | IOP paywall | MEASURED, independent of Geant4, but beamline-dominated (Gottschalk) | FACT (access), OPEN (content) |
| 3 | Tessonnier, Mairani, Brons, Haberer, Debus, Parodi 2017, PMB 62:3958 (doi 10.1088/1361-6560/aa6516) | Measured lateral profiles in water at HIT for p, He, C, O (same campaign as thesis, row 4) | OPEN in paper | Not OA (Europe PMC: OA N, hasSuppl N; LMU epub 53991: "Volltext nicht verfügbar") | IOP paywall | MEASURED, independent | FACT (access) |
| 4 | Tessonnier, T., LMU dissertation 2017, "Treatment of low-grade meningiomas with protons and helium ions", ch. 3 + App. A | MEASURED lateral profiles in water (PTW pinpoint 0.03 cm³ IC block, 24 chambers, line scans of 101 spots, experimental room, entrance window at isocentre), from ~16 mm to BP+30 mm; p 80.90/157.43/221.06 MeV, He 82.33/158.08/220.51, C 150.42/299.94/430.10, O 176.88/356.20/430.32 MeV/u (matched ranges except O highest). Single-Gauss σ0 (core >20 %) and double-Gauss σ1/σ2/w; σcorr = sqrt(σ² − σ(0 cm)²), σ(0 cm) by linear extrapolation | No uncertainty stated for σcorr or the ratios (CCE etc. elsewhere); O foci out of tolerance at the time | NO table of σ(z): Fig. 3.6 (σ0corr, σ1corr, σ2corr, w vs depth/BP) and Fig. A.3 plotted only. Numeric TEXT statements (p. 44, App. A p. 157): at medium energy, width relative to C: p 3.2×, He 1.7×, O 1/1.1×; scaling p 0.9·sqrt(Mc/Mp), He sqrt(Mc/MHe), O 1.05·sqrt(Mc/MO) (factor F ≈ sqrt(M) ± 10 %); σ0corr at BP larger than C by [2.3, 3.4, 5.1] mm (p) and [0.8, 1.3, 1.8] mm (He) for low/medium/high energy; O−C < 0.3 mm except highest energy | https://edoc.ub.uni-muenchen.de/20602/7/Tessonnier_Thomas.pdf (6 861 659 B, sha256 70df3593…); LMU edoc, free access, licence not checked (no CC seen; OPEN) | MEASURED, independent of Geant4/TOPAS; BUT σcorr contains nozzle/air-induced angular divergence (only the entrance width is subtracted), i.e. ion-optics- and beamline-confounded; Gottschalk's caveat applies | FACT |
| 5 | Besuglow et al. 2022, Front. Phys. 9:797354 (doi 10.3389/fphy.2021.797354) | He at HIT, air-gap widths (measured + FLUKA, Tables 1–2), water profiles at 107.93/158.08/220.51 MeV/u near maximum dose (OCTAVIUS 2D array, Fig. 6) | ±0.8 mm MC-vs-measurement agreement | Water: NOT tabulated (figures only); data "available from the authors" | CC BY | measured (water figure), He only | FACT (via WebFetch summary) |
| 6 | Krämer et al. 2016, Med. Phys. 43:1995 (doi 10.1118/1.4944593) | He 200 MeV/u at HIT; 3-D dose in extended targets (SOBP), fragments; "slight underestimation of the lateral spread" by TRiP98 | — | No pencil-beam σ(z); Europe PMC: CC BY, subscription link | Wiley | measured, extended fields: not an MCS test | FACT (abstract) |
| 7 | "Phillips (LBL)" data, via Kraft 1990, Hollmark et al. 2004 PMB 49:3247, Kanematsu 2008 NIM B 266:5056 (arXiv:0806.0106) and 2009 PMB (arXiv:0810.1390) | RMS displacement vs depth in water: p, He (~215 MeV/u), C (~420 MeV/u), R0 ≈ 29.4–29.7 cm, plus σy0(R0) | none | Plots only (no tables) in all citing works; Gottschalk: original LBL report not found, Phillips himself suspects data were "collated … or actually just obtained from the treatment planning code data files"; Hollmark's beam-size correction implies an unrealistic 1.2 mm beam | arXiv PDFs (sha256 60ee41ad…, 515daf1e…) | PROVENANCE UNKNOWN, possibly model-derived → not usable as independent evidence | FACT |
| 8 | Wong et al. 1990, Med. Phys. 17:163 (doi 10.1118/1.596538) | Projected angles of 650 MeV/u 238U behind 0.19 X0 Cu (Bevalac, Si position-sensitive detectors) | not explicit | OPEN (paywall) | Wiley | measured; out of scope (U in Cu, not He/C/O in water) | LIT (abstract/Gottschalk) |
| 9 | LBL/Bevalac "Multiple Coulomb scattering of high-energy heavy charged particle beams used in biology and medicine" (INIS ntp0a-gn489; Rad. Res.?) | Lateral displacement + angular distributions, Ne, Fe, U on Al, Cu, Pb, WATER | OPEN | OPEN (INIS 403/402; OSTI unreachable via proxy) | — | measured, but Ne/Fe/U only (no He/C/O); would test charge scaling at higher Z | LIT |
| 10 | Matsufuji et al. 2005 PMB 50:3393 (HIMAC, 290 MeV/u? C in thick water) | Fragment angular distributions (ΔE–E telescope); primary deflection reportedly consistent with Molière | OPEN | Not OA (Europe PMC) | IOP | measured; fragment-focused | LIT |
| 11 | GSI Schardt angular distributions of C on water (cited in DPG 2009 ST6 abstract) | used to calibrate a fragment (Goldhaber/Morrissey) model | — | not published as data (abstract mention only) | https://www.dpg-verhandlungen.de/2009/muenchen/st6.pdf (75 786 B) | fragments, not primary MCS | FACT (abstract text) |
| 12 | Mirandola et al. 2015 Med. Phys. 42:5287 (CNAO) | Spot sizes in AIR (EBT3), p 0.7–2.2 cm, C 0.4–1.1 cm FWHM | — | not water σ(z) | paywall | commissioning | LIT |
| 13 | MedAustron (Resch 2019 protons OA; Jia et al. 2025 carbon GATE-RTion vs measured transverse profiles) | protons in water (MCS model validation); carbon: lateral profiles in air + field-size factors | OPEN | OPEN | PMC (bot challenge) | commissioning | LIT |
| 14 | HIMAC/QST (Inaniwa triple Gaussian, IBIC 2014 wepf23, NIRS 2021 multi-ion He/C/O/Ne) | beam-model parameters from measurement+MC | — | OPEN | — | model-fitted | LIT |
| 15 | geant-val | Earlier repo note: Haettner entries = fragment yields only; no ion lateral-profile test found by search | — | — | — | — | LIT/OPEN |
| 16 | Kempe & Brahme 2010; Hollmark 2004; Kanematsu 2008/2009 | theory/semi-empirical; their ion "data" = Phillips (row 7) | — | — | — | theory | FACT (Kanematsu), LIT (Kempe) |
| 17 | ICRU 93 / PTCOG | no MCS tables known | — | — | — | — | LIT/OPEN |
| 18 | arXiv:2407.09851 (minibeams p/He/C/O, Geant4) | Geant4 validated vs Schwaab 2011 C 350.84 MeV/u, p 182.66 MeV in water (plateau 35.5 mm); states "experimental data on the lateral profiles … are lacking"; initial width fitted at entrance | — | data read from Schwaab plots | arXiv (2.49 MB, sha256 6590326c…) | Geant4 — not independent | FACT |
| 19 | Haberer et al. 1993 (NIM A 330:296) | scanning system; no MCS data | — | — | — | — | LIT |

## 2. Proton alternative (matched velocity / relative check)
- Gottschalk et al. 1993 NIM B 74:467 (158.6 MeV p, thin targets incl. water; already pinned) and Gottschalk 2015 halo
  (arXiv:1409.1938, 177 MeV, tabulated): measured, independent, tabulated → test the proton (z = 1) baseline at the
  velocity of He/C/O ≈ 158–177 MeV/u. FACT (known to project) / LIT.
- Preston & Koehler (Harvard, unpublished, 127–158 MeV p in water/Al; universal curve): available only as the facsimile
  figure in Gottschalk arXiv:1204.4470 Fig. 23 — plotted only. FACT.
- Limitation: proton data constrain the z = 1 baseline (and thereby f(z=1)) but cannot by themselves test the charge
  dependence; a ratio test needs ion data at matched velocity/range, which is the missing item.
- Within Tessonnier (row 4), p and ions were measured in the same setup: the p/C, He/C, O/C width ratios are
  measured relative numbers (text-quoted, no uncertainty, beamline-confounded) — the only same-facility charge-dependence
  numbers found.

## 3. Verdict
**EVIDENCE GAP (for a graded/gating row).** No measured, independent, TABULATED dataset of He/C/O lateral width vs
depth in water (or thin-water angular widths) with stated uncertainties was found. Measured He/C/O-in-water data exist
(Tessonnier 2017 PMB, paywalled; Tessonnier 2017 LMU thesis, free) but only as figures plus a few text-quoted ratios,
without uncertainties and confounded by beamline/ion-optics divergence (Gottschalk's assessment). The "Phillips" ion
data are of unknown provenance (possibly TPS files). Thin-target heavy-ion angular data (Wong 1990; Bevalac Ne/Fe/U
in water) are outside He/C/O or inaccessible.
→ Per decision 0043 §3, the lateral-gap intervention wording applies once P10 has also failed.

Optional, non-gating (orchestrator's decision; must be pre-registered as INFORMATIVE before any ionmc lateral run):
an informative relative check against the thesis' text-quoted measured ratios at medium energy (σ_rel to C:
p 3.2, He 1.7, O 1/1.1; BP differences He−C 1.3 mm, p−C 3.4 mm, O−C < 0.3 mm), with ionmc run as an ideal pencil
beam in water, declared as beamline-confounded with an assumed ≥ 10 % ratio uncertainty (from the thesis'
"F ≈ sqrt(M) ± 10 %"). Acquisition: https://edoc.ub.uni-muenchen.de/20602/7/Tessonnier_Thomas.pdf, 6.86 MB,
sha256 70df35938f57acbf6d88922593e8cb9612a0e0f83701499f181244fc9850428a, cite only (licence OPEN, no CC seen),
numbers transcribed with page references (p. 31–33 methods, p. 44, p. 157). It cannot discriminate the ~5 % model
uncertainty and must not be counted as held-out qualification evidence.

## 4. OPEN (not done before stop)
Schwaab 2011 and Tessonnier 2017 PMB table content (paywalled); Bevalac Ne/Fe/U-in-water report (OSTI/INIS
unreachable); HIMAC multi-ion papers; MedAustron Jia 2025 full text; geant-val direct API query; thesis licence;
Heidelberg/GSI theses (Schwaab Diplom, Weber 1996 GSI report).
