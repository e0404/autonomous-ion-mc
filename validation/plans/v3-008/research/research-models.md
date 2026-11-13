> **Research note (physics-researcher, 2026-11-01; not evidence).** Statements are labelled FACT / LIT / TO-VERIFY / OPEN by the author; estimates are hand estimates, nothing was run. Items marked TO-VERIFY or OPEN are resolved by the pre-freeze checks of `validation/plans/v3-008-acceptance.md` before the plan is frozen.

# V3-008 research: ion EM models (RQ-M1, RQ-S1, RQ-S2, RQ-ST1)

Status: interim report, written on coordinator instruction before the investigation was finished.
Nothing in the repository was modified. No simulations were run. All numbers below marked "estimate"
come from hand arithmetic on standard formulas, not from code. They must be reproduced by the I11/I3b
scripts before freeze.

Markers: FACT = checked in a primary source or in code during this pass. LIT = literature statement
whose source exists and was located, but whose detail was only partly read. TO-VERIFY = from memory or
from a snippet only. OPEN = not resolved.

Sources actually read in this pass:
- the repository code `src/ionmc/physics/{stopping,em,scattering,projectiles}.py` (worktree v3-005b);
- decisions 0038/0039 (grep);
- Gottschalk, arXiv:0908.1413 (local LaTeX copy in scratchpad `v3-005c-downloads/arxiv-0908.1413/x/ScatPowerV2.tex`);
- the PDG 2024 review "Passage of particles through matter", text extracted to `scratchpad/v3-008/dl-pdg/pdg.txt`;
- Kats, Halenka, Voskresenskaya, arXiv:2008.04819 (text extracted to `dl-kats/kats.txt`, read pp. 1–5);
- the MARS15 paper arXiv:physics/0505043 (LaTeX, `dl-mars/x/`);
- Geant4 v11.2.0 `G4EmCorrections.cc`, `G4LindhardSorensenIonModel.cc`, `G4LindhardSorensenData.cc` (`dl-g4/`).
  `G4UrbanMscModel.cc` was downloaded but not read.

## 0. Facts about the existing code that the answers rely on

- FACT: the proton MCS model is Gottschalk's differential Molière power,
  `T_dM = f_dM(pv, p1v1) (z E_s/pv)^2 / X_S` with `E_s = 15 MeV`. Here
  `f_dM = 0.5244 + 0.1975 lg(1-(pv/p1v1)^2) + 0.2320 lg pv - 0.0098 lg pv lg(1-(pv/p1v1)^2)`.
  - `X_S` is Gottschalk's scattering length (Bragg-additive).
  - The angle is applied as a Gaussian with a random hinge and a birth-step average (decision 0039,
    `em.py`).
  - f_dM is a fit of T_Hanson/T_IC to Molière/Fano/Hanson **theory** for protons (Be, Cu, Pb;
    E ≤ 300 MeV). Gottschalk (sec. TdM) states it is "not explicitly material dependent".
  - `z` already enters the prefactor only.
- FACT: the Mott term in `stopping.py` is `0.5 π α β q`. Ahlen's first-order Mott term, written in the
  4πr_e² convention of `S = K q² (Z/A) β⁻² [L…]`, is the same quantity.
  - Geant4 `G4EmCorrections::MottCorrection` returns `π α β z` and adds it inside a 2× bracket
    (`sum = 2(Barkas + Bloch) + Mott`, multiplied by `2π r_e² mc² n_e q²/β²`, lines 190/242/715).
    That is the same term in the other convention, so the two implementations are consistent.
- FACT (labelling): the "Pierce–Blann" charge in the code is `z(1 − exp(−125 β z^(−2/3)))`. This is
  the Barkas (1963) form. Pierce & Blann (1968) use `0.95 v/(v0 z^(2/3))`, i.e. a coefficient of
  0.95·137.04 ≈ 130.2 instead of 125.
  - Gottschalk's footnote calls the 125 form "Barkas formula" and cites Kraft and Barkas.
  - Recommend the decision names it "Barkas/Pierce–Blann-type (coefficient 125)".
  - The two coefficients differ by about 1.8 % in q² for C at 2 MeV/u (estimate) and by < 0.1 % above
    10 MeV/u.

## RQ-M1 — MCS for He/C/O at 50–430 MeV/u in water (decisive)

### Physics: how the Molière quantities scale with z, A and β

All of the following is FACT as standard Molière/Bethe 1953 theory; the explicit forms are taken from
Lynch & Dahl as cited by Gottschalk and PDG (TO-VERIFY against Lynch & Dahl 1991, NIM B 58, 6).

- **Characteristic angle:** `χ_c² = 0.1569 z² Σ w_i Z_i(Z_i+1)/A_i · t / (pv)²`. It depends on the ion
  only through z²/(pv)².
- **Screening angle:** `χ_a² = χ_0² (1.13 + 3.76 α'²)`, with:
  - `χ_0 = ħ/(p a_TF) ∝ Z^(1/3)/p`. This depends on **p**, not pv.
  - `α' = z Z α/β`, the Born parameter.
  - Lynch & Dahl's equivalent form is `2.007e-5 Z^(2/3)(1 + 3.34 α'²)/p²`.
- **Collision number:** `Ω = χ_c²/(1.167 χ_a²) ∝ z² t / [β² (1.13 + 3.76 α'²)]`. B follows from
  `B − ln B = ln Ω`. The Hanson width is `θ_0 = χ_c √(B − 1.2)/√2` (Gottschalk eq. thetaHanson,
  FACT).
- **What this means for the ion width.** Ω is mass-independent, and for an ion at the same velocity and
  the same x/R the thickness is `t_ion = (A/z²) t_p`. Hence `Ω_ion/Ω_p = A / F_scr`. Here F_scr is the
  ratio of the Born screening factors (1.13 + 3.76 α'²) of ion and proton, weighted by Z(Z+1).
  Water's Z(Z+1) sum is dominated by O: 72 of 76.
- **Born-parameter screening factor for an O target (estimate):**
  - C at 290 MeV/u: α' = 0.54 → F_scr ≈ 1.9.
  - O at 350 MeV/u: α' = 0.68 → F_scr ≈ 2.5.
  - He at 150 MeV/u: α' = 0.23 → F_scr ≈ 1.15.
  - Proton: α' ≈ 0.09 → factor 1.03.
- **Resulting shift in ln Ω (estimate):** about +1.8 for C, +1.9 for O and +1.25 for He relative to the
  velocity-matched proton. This gives ΔB ≈ +1.4 to +2.0. With proton B ≈ 12–14 for therapeutic water
  slabs (TO-VERIFY), the ion width is larger than z/(pv)-scaled proton theory by about +8 % (C, O) and
  +5–6 % (He).
- **Projectile screening.** The ions are fully stripped in the transported domain: q/z > 0.99 above
  about 10–15 MeV/u (FACT from the formula). Projectile electrons therefore play no role. Screening is
  by target electrons only, which is already in χ_a.
- **Nuclear size / form factor.** For C at 290 MeV/u after about 10 cm, the core momentum transfer is
  p·χ_c ≈ 20 MeV/c (estimate). The form-factor scale ħc/R_ch(O) is about 70 MeV/c. The form factor
  therefore cuts the single-scattering tail beyond roughly 3χ_c. The effect on the Gaussian (1/e,
  Hanson) core width is expected at the ≤ 1–2 % level (estimate, OPEN). It grows with z and becomes
  important only for very heavy ions. Wong et al. 1990 (Med. Phys. 17, 163, DOI 10.1118/1.596538)
  measured MCS of 650 MeV/u U on Cu (0.19 X0) against Fermi, Molière and modified Fermi theory: FACT
  that the paper exists and what it measured; its conclusions are OPEN (full text not read). It is not
  water and not light ions, so it is weak evidence for V3-008.
- **Atomic-electron (Fano) term.** In momentum-transfer terms the maximum transfer to an electron,
  2 m_e c βγ, is mass-independent. Relative to the ion's larger p·χ_c, the electron contribution
  (13 % of χ_c² in water) becomes less "logarithmically efficient". The resulting change is a −1 to
  −2 % effect on θ² relative to treating electrons like nuclei (estimate, TO-VERIFY). This belongs in
  the θ_MH reference computation (use the Fano form with β-dependence, as Gottschalk 1993 does for
  protons).

### Candidates and expected deviation from Molière–Hanson (estimates; to be confirmed by I11)

Two sensitivities underlie the table:
- f_dM in protons responds by about 0.13 per unit ln Ω at fixed x/R. This comes from the lg pv
  coefficient 0.232 and Ω ∝ (pv)^≈0.77 at fixed x/R (estimate).
- Molière gives d ln(B−1.2)/d ln Ω ≈ 1.08/(B−1.2) ≈ 0.09.

| Candidate | What it does to the log term | Expected θ vs θ_MH |
|---|---|---|
| (a) velocity-scaled T_dM (planner's first choice) | treats the ion as the velocity-matched proton: Ω shift 0 | **−5 % (He), −7 to −8 % (C, O)**. Expected to FAIL the 3 % rule |
| (b) naive T_dM with the ion pv | `lg pv_ion = lg pv* + lg A`: Ω shift ≈ ln A, misses F_scr | about +2 to +4 % (over-predicts by roughly ½ ln F_scr / (B−1.2)) |
| (c) Highland/Lynch–Dahl with z² in the log (Kanematsu-style differential form) | ln(x z²/X0β²) at fixed x/R shifts by ln A, misses F_scr | similar to (b): +1 to +3 %. PDG 2024 (FACT, `pdg.txt` l. 686–690): the formula "takes into account the p and z dependence quite well at small Z, but for large Z and small x the β-dependence is not well represented" |
| (d) T_dM with a Molière Ω correction: velocity-scaled arguments plus an additive term `Δf = c · ln(A/F_scr(β, z))` (or, equivalently, pv_arg chosen so that f_dM sees the Molière Ω) | matches Molière Ω by construction | ≤ 1–2 %, limited by f_dM's own fit error. The constant c (≈ 0.13) and the F_scr weighting must be derived from Molière theory, never from TOPAS |
| (iv) Urban/Geant4 | `G4UrbanMscModel::ComputeTheta0` is believed to use a Highland-like `13.6·|charge|·√y/(βcp)·(c1 + c2 ln y)` with y = step/X0 and **no z² in the log** (TO-VERIFY; file is in `dl-g4/` but was not read) | if confirmed, Urban (TOPAS) ion widths sit on the (a) side, i.e. several % below Molière |

Lineage notes:
- Kanematsu 2008 (NIM B 266, 5056) "treats heavy ions as well as protons" with T_dH (FACT: Gottschalk
  l. 70). Its exact ion form is TO-VERIFY.
- Kanematsu 2009 (PMB 54 N67, arXiv:0810.1390) gives a water-referenced T_LD ∝ 1/R_W for protons and
  heavy ions (LIT). It is semi-empirical (it calibrates E_s down); do not use it as the reference.

### Recommendation

1. Keep the planner's structure: a numpy Molière–Hanson ion reference (I11), validated first against
   Gottschalk's proton θ_Hanson tables. Pre-declare **the reference definition**:
   - χ_c with Z(Z+1) and the Fano electron term;
   - χ_a with the Born term `1.13 + 3.76 α'²` per element;
   - log-averaged χ_a along the slab, weighted by 1/(pv)², as in Gottschalk 1993;
   - point nucleus, with the nuclear-size effect reported separately as an uncertainty.
2. Reorder the pre-declared selection to **(d), (c), (b), (a)**, or keep the rule "first candidate
   within the bound" with (a) listed last. My estimate is that (a) misses by 5–8 % for C/O. If the
   order is kept, the expected outcome is that (a) fails and the rule falls through: acceptable, but
   predictable now.
   - (d) is the physically correct generalisation of the existing model. It keeps step independence,
     the random hinge, the clamp, the birth-step average and U5 by construction, because it only
     changes the argument of f_dM.
   - Its correction constant comes from theory: c from the f_dM proton fit sensitivity, F_scr from
     Molière. The correction should be committed as a formula with derivation, not fitted.
3. **Tolerance for I11: 3 % on θ (projected rms) is defensible from theory, not TOPAS.**
   - f_dM reproduces θ_Hanson for protons within about ±2 % (Gottschalk Fig. ThP; TO-VERIFY the exact
     band).
   - The extra theory uncertainty for ions is about 1–2 %, from the Born-screening form at α' ≈ 0.5–0.7,
     the Fano/electron term and the nuclear form factor (estimates).
   - Adding these in quadrature gives about 2.5–3 %.
   - Do not tighten below 3 %. A failure of every candidate at 3 % should STOP, as the plan says.
4. **V2 (TOPAS lateral, 5 %) needs attention before freeze (OPEN, important).**
   - If Urban's ion θ0 indeed lacks the z²/Born-dependent log term, TOPAS will be about 5–8 % narrower
     in angle than a Molière-consistent ionmc for C/O. In σ_x at 0.5–0.9 R that is roughly the same
     fraction.
   - V2 could then fail because the reference model is wrong, not ionmc.
   - Before freeze, read `G4UrbanMscModel::ComputeTheta0` and the opt4 ion msc assignment for
     GenericIon/alpha (Urban vs WentzelVI). Pre-declare the expected TOPAS/Molière offset per ion, and
     state in the plan that V2 is evidence of shared-Highland lineage, not an independent test of the
     z-dependence.
   - Measured water lateral data for C (e.g. Parodi et al. 2013 J. Radiat. Res. double-Gaussian
     parametrisation at HIT, located; Schwaab et al. 2011 not located; Tessonnier et al. 2017) are the
     only independent arbiter. Their availability as numbers is OPEN.
   - Kempe & Brahme 2010 was not checked (OPEN).

## RQ-S1 — Lindhard–Sørensen vs Bloch + first-order Mott for C and O

Formula and literature:
- FACT: the point-nucleus LS correction (Kats et al. eq. 10; MARS15 eq. LS) is
  `ΔL_LS = Σ_k [ k/η² (k−1)/(2k−1) sin²(δ_k − δ_{k−1}) + k/η² (k+1)/(2k+1) sin²(δ_{−k} − δ_{−k−1})
  + k/(4k²−1) · 1/(γ²k² + η²) − 1/k ] + β²/2`, with η = αz/β and relativistic Coulomb phases δ_k.
- I checked by hand that in the non-relativistic limit the bracket sums to
  `−η² Σ 1/(k(k² + η²))`, i.e. the Bloch term exactly. LS therefore contains Bloch, and
  LS − Bloch = the exact (all-order) Mott correction (Kats: LS and the exact Mott–Bloch correction
  "coincide up to the seventh significant digit", FACT from the abstract).
- FACT (Kats p. 2, citing Eby & Morgan): exact Mott calculations "demonstrated the significance of
  Mott's correction to the Bethe−Bloch formula for incident nuclei with Z ≥ 20".

Size estimates:
- The first-order Mott term is 0.5 π α z β / L. For C at 290 MeV/u that is 0.0446/8.7 ≈ **0.5 % of S**;
  for O at 430 MeV/u about **0.7 %** (estimate). This is already in the code.
- The residual LS − (Bloch + first-order Mott) is the higher-order Mott series, O((αz)²) with
  coefficients of order 1. (αz)² is 0.0019 (C) and 0.0034 (O); divided by L ≈ 7–9 this gives roughly
  **≤ 0.05–0.1 % of S** across 80–430 MeV/u (estimate, TO-VERIFY numerically).
- The integrated CSDA-range effect is ≤ 0.1 %, below the planner's 0.2 % trigger.
- Finite nuclear size enters only at γ ≳ 10 (MARS15, FACT), so it is irrelevant here.

Context:
- ICRU 73 high-energy tables were computed with the program BEST, and the report states that the LS
  program and the tables rest on first principles (LIT, core.ac.uk snippet; TO-VERIFY whether BEST
  applies LS).
- The ICRU 90 carbon method is OPEN.
- Geant4's `G4LindhardSorensenIonModel` adds a tabulated `ΔL(Z, γ)·Z^(2/3)` above its low-energy
  limit. It is a fit table; I did not decode its normalisation, so do not use it as a reference
  (shared lineage anyway).

Recommendation:
- Omit LS and document the bound "|LS − (Bloch + Mott₁)| ≤ 0.1 % of S, ≤ 0.1 % of R for z ≤ 8 at
  ≥ 10 MeV/u (estimate)". Keep the I3 tolerance at 1.5 % in S; it is not consumed by this term.
- Before freeze, confirm with a small numpy script (complex lgamma for the δ_k; about 1e4 terms plus a
  tail; test against the Bloch limit at small z). Report it as an I3b informative row "LS minus
  Bloch+Mott₁ per energy". If the script shows > 0.2 % in range, the planner's option rule applies.

## RQ-S2 — effective charge below 10 MeV/u and the 2 MeV/u cutoff

Effective charge (estimates from the formula):
- Above 10 MeV/u: q/z ≥ 0.99 for C and ≥ 0.98 for O.
- At 2 MeV/u: q/z ≈ 0.92 (C, q² −16 %), 0.87 (O, q² −24 %) and 0.994 (He).
- Alternative forms (Hubert, Bimbot & Gauvin 1990, ADNDT 46, 1; Ziegler/Brandt–Kitagawa;
  TO-VERIFY, not read) differ from the Barkas-125 form by a few % in q² at 2–5 MeV/u for C/O and by
  < 0.5 % above 10 MeV/u (TO-VERIFY). The 125 vs 130.2 coefficient alone gives 1.8 % in q² at 2 MeV/u
  for C (FACT, arithmetic).

Residual range at the 2 MeV/u cutoff (velocity scaling R_ion = (A/z²) R_p(β), with PSTAR
R_p(2 MeV) ≈ 0.0075 g/cm² TO-VERIFY, corrected upward by about 1/q̄²):
- He ≈ 0.075 mm; C ≈ 0.025–0.03 mm; O ≈ 0.02–0.025 mm (estimates; the decision's numbers are
  consistent).
- As a fraction of the CSDA range at the lowest domain energy (R_p(50/80/100 MeV) ≈ 22/52/77 mm,
  memory):
  - He 50 MeV/u: R ≈ 22 mm → 0.34 %;
  - C 80 MeV/u: R ≈ 17 mm → 0.17 %;
  - O 100 MeV/u: R ≈ 19 mm → 0.13 %.

What the cutoff error is:
- The cutoff is not a range error. Energy is conserved and deposited locally. The only effect is
  that the last ≤ 0.03–0.075 mm of dose is moved upstream by at most that distance.
- On R80 the effect is ≪ 0.01 mm, because the distal edge is set by straggling σ of order 0.1–1 mm.
  I6/C-res already cover it.
- Effective-charge model spread from 20 to 2 MeV/u: R_C(20 MeV/u) ≈ 1.4 mm, of which the 2–10 MeV/u
  part is about 0.4 mm. A 2–4 % q² spread there moves the end of range by ≤ 0.01–0.02 mm (estimate).
  This is below the decision's 0.05 mm threshold, so no change is needed.
- Recommendation: document it; make the I3 2–10 MeV/u band (6 %) informative or keep it as is (it is
  defensible: model spread plus PASS/ICRU 73 low-energy behaviour).

## RQ-ST1 — energy-loss straggling for ions

Bohr variance and charge state:
- Bohr variance with bare z and the relativistic factor Tmax(1 − β²/2) (FACT, `em.py`) is the right
  leading model in the domain: ions are fully stripped above about 10 MeV/u.
- Charge-exchange straggling arises only where capture and loss are active, i.e. below a few MeV/u for
  C/O (TO-VERIFY: Sigmund, "Particle Penetration and Radiation Effects"; ICRU 73 straggling chapter).
- The range interval below 10 MeV/u is ≈ 0.4 mm for C and ≈ 0.3 mm for O. Its contribution to σ_R²
  is ≪ 1 % of the total variance (estimate), so it is negligible for σ_R.
- Using q² instead of z² would lower the variance only in that last segment: a negligible change,
  consistent with the planner's rejection.

Other deviations from Bohr:
- Shell/binding corrections (Chu; Lindhard–Scharff) and bunching or correlation effects for light ions
  in water at ≥ 50 MeV/u are expected at the few-% level in σ² (TO-VERIFY: Yang, O'Connor & Wang 1991
  NIM B 61, 149; Tschalär). That is about 1–2 % in σ_R, below the planner's 10 % amendment trigger.

Sampler validity:
- Vavilov κ = ξ/Tmax per 1 mm step (estimate):
  - C 290 MeV/u: ξ ≈ 0.73 MeV, Tmax ≈ 0.74 MeV → κ ≈ 1;
  - proton 150 MeV: κ ≈ 0.1.
- κ ∝ z² at equal velocity, so ions sit closer to the Gaussian regime than protons. mean/σ is large, so
  the Gamma shape k = (mean/σ)² is large and its skewness 2/√k is small.
- The `bohr_gamma_v1` sampler (exact Bohr mean and variance, additive shape) remains valid. Its
  step-independence argument carries over unchanged. Delta tails are irrelevant because deltas are
  deposited locally.

What the plan rows should assert:
- I4b: the scaling identity σ_R/R(ion) = (m_p/M)^½ σ_R/R(p) at equal velocity, corrections off.
  This is exact algebra (FACT), so 1e-3 is fine.
- I8: within 5 % of the quadrature (self-consistency), fine.
- Add an informative row comparing σ_R/R at the standard points with the literature statement
  "σ_R/R(C) ≈ σ_R/R(p)/√12 ≈ 0.3–0.4 %" (Schardt, Elsässer & Schulz-Ertner, RMP 82 (2010) 383,
  TO-VERIFY). Use a ±15 % band, report only.
- V1(b) distal width within 10 % vs TOPAS is defensible:
  - Geant4 ion fluctuation models (G4IonFluctuations/Universal) are independent of the ionmc
    implementation but share Bohr theory.
  - The theory spread in σ_R is a few %.
  - Beam energy spread is zero in both codes for this case.
- Measured carbon Bragg-peak widths in water mix beam energy spread and fragment tails into the
  width; use them only as informative evidence.

## OPEN list (not verified in this pass)

1. Numerical I11 Molière–Hanson ion computation. All RQ-M1 percentages are hand estimates (±3 %).
2. Content of `G4UrbanMscModel::ComputeTheta0` and the opt4 msc model for GenericIon/alpha. This is
   decisive for interpreting V2.
3. Wong et al. 1990 conclusions; Kempe & Brahme 2010; Kanematsu 2008 ion form of T_dH; Lynch & Dahl
   constants.
4. Measured lateral σ of C/He/O in water, usable as numbers (Parodi 2013 JRR, Schwaab 2011,
   Tessonnier 2017).
5. Numerical LS − (Bloch + Mott₁) for z = 6, 8 (script). Whether ICRU 73 BEST and ICRU 90 C tables
   include LS.
6. Hubert et al. 1990 effective-charge form and its difference from Barkas-125 in water.
7. Literature size of non-Bohr straggling (Tschalär, Yang 1991) and charge-exchange straggling for C/O
   in water.
8. Proton B range of 12–14 for therapeutic water slabs (sets the size of the RQ-M1 correction).
