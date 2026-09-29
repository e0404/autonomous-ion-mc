# Nuclear interaction models and measured data (research note, 2026-09-30)

Source: physics-researcher subagent; formulas in section 1 were checked against
Geant4 v10.4.3 source and the open MCsquare source; sections 3–5 and all
dataset locations are **unverified** recollections pending first-hand lookup.

## 1. Total reaction cross sections (verified against Geant4 source)

Tripathi universal (NASA TP-3621, 1997; NIM B 117 (1996) 347):
σ_R = π r0² (A_P^⅓ + A_T^⅓ + δ_E)² (1 − R_c B/E_cm) X_m, r0 = 1.1 fm;
δ_E = 1.85 S + 0.16 S/E_cm^⅓ − C_E + 0.91 (A_T − 2Z_T) Z_P/(A_T A_P);
S = A_P^⅓ A_T^⅓/(A_P^⅓ + A_T^⅓);
C_E = D (1 − e^(−E/T1)) − 0.292 e^(−E/792) cos(0.229 E^0.453), E in MeV/u lab,
E_cm the c.m. kinetic energy in MeV; B = 1.44 Z_P Z_T/R,
R = r_P + r_T + 1.2 (A_P^⅓ + A_T^⅓)/E_cm^⅓, r_i = 1.29 r_rms,i (Geant4 uses
r_rms ≈ 0.6·1.36·A^⅓ fm). Heavy systems: T1 = 40, X_m = 1, R_c = 1; D = 1.75
(2.05 for protons; alpha: 2.77 − 8e-3 A_T + 1.8e-5 A_T² − 0.8/(1 + e^((250−E)/75))).

Tripathi light systems (NASA TP-1999-209726; NIM B 155 (1999) 349) for p, n, d,
³He, ⁴He as projectile or target: X_m = 1 − X1 e^(−E/(X1 S_L)),
X1 = 2.83 − 3.1e-2 A + 1.7e-4 A² (A = heavier partner),
S_L = 1.2 + 1.6 (1 − e^(−E/15)); p/H partner: T1 = 23,
D = 1.85 + 0.16/(1 + e^((500−E)/200)); ⁴He: T1 = 40, G = 75 (Z_T = 7: T1 = 40,
G = 500), D = 2.77 − 8e-3 A + 1.8e-5 A² − 0.8/(1 + e^((250−E)/G));
R_c = 13.5 (p+d), 21 (p+³He), 27 (p+⁴He), 2.2 (p+Li), 6.0 (d+¹²C), else 1.
Geant4 reference: `G4TripathiLightCrossSection.cc`, `G4TripathiCrossSection.cc`
(https://github.com/Geant4/geant4/tree/v10.4.3/source/processes/hadronic/cross_sections/src).

Sihver (PRC 47 (1993) 1225): σ = π r0² [A_P^⅓ + A_T^⅓ − b0 (A_P^−⅓ + A_T^−⅓)]²,
r0 = 1.36 fm, b0 = 1.581 − 0.876 (A_P^−⅓ + A_T^−⅓); energy independent
(Geant4 applies only above 100 MeV/u). Kox (PRC 35 (1987) 1678) and Shen
(NPA 491 (1989) 130) forms also recorded in the subagent transcript.

Plausibility anchors (unverified): ¹²C+¹²C σ_R ≈ 850–900 mb at 200–900 MeV/u;
¹²C+p ≈ 220–240 mb. Beam attenuation experiments count the Z = 6 primary and
therefore measure the charge-changing cross section σ_cc < σ_R.

Protons: MCsquare (Apache-2.0, https://gitlab.com/openmcsquare/MCsquare)
samples ICRU 63 tables (`Materials/<element>/ICRU_Nuclear_{elastic,inelastic}.dat`,
7–249 MeV, p/d/α multiplicities, 13-angle spectra, local recoil fraction) in
`src/compute_nuclear_interaction.c`; pp elastic macroscopic cross section
Σ_pp/ρ = 0.315 E^(−1.126) + 3.78e-6 E cm²/g (E > 10 MeV), isotropic in c.m.
ICRU 63 tables are ICRU-copyright evaluated data; open alternatives (unverified):
ENDF/B proton sublibrary, JENDL/HE. Fippel & Soukup 2004 (Med Phys 31 2263) is
paywalled; its formulas were not recovered from open sources.

## 3. Fragmentation models (unverified)

EPAX is unsuitable for ¹²C/¹⁶O. Abrasion–ablation (Bowman–Swiatecki–Tsang;
NUCFRG2) gives σ(Z,A) at ≈20–30% accuracy. Goldhaber (PLB 53 (1974) 306)
momentum widths σ∥² = σ0² A_F (A_P − A_F)/(A_P − 1), σ0 ≈ 90 MeV/c empirical;
adequate for Z ≥ 3 forward fragments, underpredicts wide-angle light fragments.
Recommended GPU hybrid: tabulated/parameterised partial cross sections fitted to
measured charge-changing data; Goldhaber kinematics for Z ≥ 3; energy–angle
parameterisation (moving source) or tabulation for Z ≤ 2; per-event energy balance.

## 4. Measured evidence candidates (unverified availability)

Haettner 2006/2013 (¹²C 200/400 MeV/u in water: fragment yields vs depth);
Gunzert-Marx 2008 NJP (200 MeV/u spectra behind 12.78 cm water, likely open);
Toshito 2007 PRC (charge-changing σ, 200–400 MeV/u ¹²C on water);
Golovchenko 2002 PRC; Dudouet 2013 PRC (95 MeV/u double-differential);
Divay 2017 PRC (50 MeV/u); FOOT collaboration (C/O 200–400 MeV/u);
Rovituso 2017 PMB and Horst 2019 PRC (⁴He); Bragg curves: Sihver, Schardt &
Kanai 1998; Tessonnier 2017 PMB (H/He/C/O at HIT); Kurz 2012 PMB (¹⁶O);
Matsufuji 2003/2005 PMB (fragment fluence).

## 5. Proton evidence (unverified)

Bortfeld 1997 analytic Bragg curve (related theory, shared stopping lineage);
Gottschalk 1993 NIM B 74 467 measured MCS in many materials (tables);
no open numerically tabulated clinical proton depth-dose dataset confirmed.

## Recommendation summary

Tripathi 1997 + 1999 for all ion/fragment–target pairs; tabulated proton
non-elastic data; transport all charged fragments; declare He and O domains
separately; never use the same measured set for calibration and evaluation;
note that FLUKA/Geant4 were tuned on most of these measured sets.
