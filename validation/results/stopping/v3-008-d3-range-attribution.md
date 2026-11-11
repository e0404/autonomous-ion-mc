# V3-008 D-3: range attribution of the TOPAS-EM R80 excess (report-only, exploratory, 2026-11-11)

Script `validation/scripts/ion_em/d3_range_attribution.py`, aggregates `v3-008-d3-range-attribution.json` (input and script
sha256 inside). Amendment R2-STOP3 item 1(b). Never a pass criterion; nothing fitted; no reference values stored. All ranges from
E_cut = 2 MeV/u, ionmc exact log-log quadrature, 400 points/decade. The C/O "TOPAS-like" stopping is an APPROXIMATION of the
Geant4 construction (S = S_ICRU90,p(equal velocity) q_eff^2 (1 + (dL_LS + dL_Barkas)/L_p); LS at the bare Z minus the proton's
own Bloch + Mott; Barkas = (q - 1) L1 of ionmc's form; exact Geant4 forms not verified, D-4/D-5 pending). All columns in % of
R_ionmc (Bethe-78, all options, as transported); residual = TOPAS excess - construction (He: ICRU 90 alpha direct; C/O: LS+Barkas
primary) - pre-step (2 % rule). Held-out energies have no TOPAS-EM R80 (no residual).

| ion | E (MeV/u) | R_ionmc (mm) | ICRU90 alpha direct (He) | TOPAS-like plain | TOPAS-like half-LS+Barkas | TOPAS-like LS+Barkas (primary) | pre-step (2 % rule) | pre-step (pure Euler) | TOPAS R80 / ionmc R80 | residual (primary) |
|---|---|---|---|---|---|---|---|---|---|---|
| he | 50 | 22.449 | +0.055 | - | - | - | +0.208 | +0.528 | +0.290 | +0.028 |
| he | 100 | 77.890 | +0.051 | - | - | - | +0.094 | +0.185 | - | - |
| he | 150 | 159.170 | +0.076 | - | - | - | +0.055 | +0.100 | +0.220 | +0.088 |
| he | 250 | 382.527 | +0.077 | - | - | - | +0.027 | +0.047 | +0.215 | +0.110 |
| carbon | 80 | 17.485 | - | -0.078 | -0.095 | -0.026 | +0.232 | +0.709 | +0.091 | -0.115 |
| carbon | 150 | 53.067 | - | +0.076 | -0.023 | -0.081 | +0.120 | +0.275 | - | - |
| carbon | 290 | 163.193 | - | +0.132 | -0.028 | -0.170 | +0.052 | +0.101 | +0.294 | +0.413 |
| carbon | 430 | 309.302 | - | +0.236 | +0.046 | -0.131 | +0.031 | +0.057 | +0.267 | +0.368 |
| oxygen | 100 | 19.553 | - | -0.120 | -0.134 | -0.055 | +0.215 | +0.643 | +0.020 | -0.140 |
| oxygen | 200 | 65.448 | - | +0.155 | +0.000 | -0.113 | +0.102 | +0.233 | - | - |
| oxygen | 350 | 166.605 | - | +0.252 | +0.019 | -0.192 | +0.049 | +0.099 | +0.291 | +0.433 |
| oxygen | 430 | 231.820 | - | +0.330 | +0.073 | -0.165 | +0.038 | +0.074 | +0.300 | +0.427 |

Other C/O variants and residuals are in the JSON (`barkas_only`, `halfls_barkas` = dL_LS x 0.5, the unverified half-magnitude
reading of the Geant4 LS table).

## Conclusions (labelled estimates, not evidence)

* **He (alpha uses the ICRU 90 alpha array directly, ESTABLISHED by D-4 preliminary):** the construction difference is +0.05 to
  +0.08 % (ionmc range shorter than integrated ICRU 90; decision 0038's <= 0.07 % confirmed). Explained: +0.076 / +0.077 of the
  +0.22 / +0.215 % excess (He 150 / 250), plus the pre-step estimate +0.055 / +0.027 %. **Unexplained (estimate): +0.09 / +0.11 %**
  (about 0.1 % = 0.15 / 0.4 mm), i.e. about 40-50 % of the He excess.
* **C/O 290-430 (+0.27 ... +0.30 %):** the explained part depends entirely on the unverified ion terms. ICRU 90 proton x q^2
  alone lengthens the TOPAS-like range by +0.13 (C 290) ... +0.33 % (O 430) relative to ionmc, the observed sign and size
  (residual after pre-step +0.0 ... +0.1 %). With the LS + Barkas terms as constructed here the sign flips to -0.13 ... -0.19 %
  (shorter) and the unexplained part grows to +0.37 ... +0.43 %; with half the LS magnitude it is +0.19 ... +0.27 %. The
  LS/Barkas terms enter at +/-0.2 ... 0.5 % (as large as the effect), so the sign of the net construction difference is NOT
  determined by this work; the decision needs the Geant4 dE/dx table (D-5) or the exact G4LindhardSorensenIonModel/Barkas forms (D-4).
  Robust: the ICRU 90 proton base + q^2 scaling alone points the same way as the observation; the pre-step lengthening adds only
  +0.03 ... +0.05 % at >= 290 MeV/u.
* **Pre-step lengthening (0.1 mm steps, 2 % rule):** +0.03 ... +0.06 % for R >= 150 mm (the amendment's estimate holds); the
  unconditional Euler bound is +0.05 ... +0.10 %. At the low ends (+0.2 %) it is dominated by the last millimetres and is not
  meaningful for R80.
* **Low ends (He 50, C 80, O 100):** residuals +0.03 / -0.12 / -0.14 % are inside the R80-vs-CSDA and E_cut-convention
  ambiguity (the sub-2 MeV/u residual range is 0.04-0.09 mm = 0.2-0.5 % of these ranges and is not in the integral).
* Net: nothing here supports changing ionmc's stopping; at most ~0.1 % (He) to ~0.3-0.4 % (C/O, only if the Geant4 ion terms are
  as large as the primary variant) remains unattributed.
