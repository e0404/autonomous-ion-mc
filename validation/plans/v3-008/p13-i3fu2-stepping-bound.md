# P13: stepping-bias bound of the I3-FU-1 step-2 thin-slab extrapolation (research record, 2026-11-16)

Status: research record written by the physics researcher at the orchestrator's request after Codex REVIEW-c4f3fdd1
(finding 2) on amendment I3-FU-1 revision 3 asked for "a quantitative bound derived and independently verified from the
actual Geant4 algorithm". It is NOT evidence and changes no grade; it supports the declared 0.2 % residual-stepping term
of the step-2 method systematic (plan rule 2a) and the refined along-step emulation of `i3fu2_evaluate.py` (revision 4).
Labels FACT (Geant4 11.4.2 source line read) / ESTIMATE (model arithmetic with the committed ionmc model; script
`validation/plans/v3-008/p13-scripts/p13.py`, model-only, reads no reference data) / OPEN as stated by the author. The
Geant4 files fetched for this record are identified by the sha256 of the fetched bytes in §1; the acquisition tool returned
an opaque error for Geant4 source URLs on 2026-11-16 (three attempts), so these files are NOT in the provenance cache —
their line citations are bound to the recorded hashes only (disclosed limitation; a later successful acquisition can
re-bind them). No Geant4 data values are reproduced. Not a clinical claim.


Scope: amendment I3-FU-1 step 2 "Quantity" bullet / rule 2a (plan `validation/plans/v3-008-acceptance.md`, heading
"### Amendment I3-FU-1"). Question: residual stepping bias of r0 = (f2 r(f1) - f1 r(f2))/(f2 - f1), f1 = 3 %, f2 = 4.8 %,
for a slab traversed in one ionIoni step (range-table branch). Researcher record, 2026-11-16; no worktree edits, no
simulation; model arithmetic only. Labels: FACT (source line read), ESTIMATE (model arithmetic), OPEN.
Status: complete except the items in section 7 (OPEN), which the coordinator asked to close out at the turn limit.

## 1. Sources (Geant4 v11.4.2)
- Cached (`.ionmc-cache/reference/<sha256>`, hashes as D-4 section 1): [VELP] G4VEnergyLossProcess.cc `2f53a58c...`;
  [LS] G4LindhardSorensenIonModel.cc `e16ed6cd...`; [BB] G4BetheBlochModel.cc `5dc2a900...`; [IEC]
  G4ionEffectiveCharge.cc `891854a9...`; [OPT4] `04337300...`; [EMP] G4EmParameters.cc `aaa143c1...`.
- Fetched for P13 (raw.githubusercontent.com/Geant4/geant4/v11.4.2/source/processes/electromagnetic/... unless noted;
  full sha256 of the fetched bytes, computed by the orchestrator on the researcher's copies; NOT in the provenance
  cache — see the status paragraph):
  | Tag | Path (under source/) | sha256 |
  |---|---|---|
  | [VELPh] | processes/electromagnetic/utils/include/G4VEnergyLossProcess.hh | 23e9dd285fbae71c1a11e682542ca790d745a68580a732cd4cb7962adf9769e4 |
  | [LTB] | processes/electromagnetic/utils/src/G4LossTableBuilder.cc | 71c24723f3eaddd14642f2a38e714151998b04afa1c183e98a54d9a8fd4b9aed |
  | [LTBh] | processes/electromagnetic/utils/include/G4LossTableBuilder.hh | 5dd6fcd762c2af045c9a8b7ffc22e82ddeb36ce33b5771301b1d26b06c8b7549 |
  | [LTM] | processes/electromagnetic/utils/src/G4LossTableManager.cc | 80211e548c33d656e44e041362bdb63a1e4a4f4e444f1fad6adc9cdb563aabf0 |
  | [PVi] | global/management/include/G4PhysicsVector.icc | 7cb2f12424df6598ff6ba320666dbd29510d1a727fbdac6a49d3fbac0047ae8f |
  | [PVc] | global/management/src/G4PhysicsVector.cc | 6d1437c9d23ed63ecbabc96ed09d8fd68789e5f4bb0575542dbd9756d66de943 |
  | [PVh] | global/management/include/G4PhysicsVector.hh | 275ad147aca9f75a65a44160ca4884b0f866c30e391eadac71cd9a8ebd37adc9 |
  | [BR] | processes/electromagnetic/standard/src/G4BraggModel.cc | ceb5158a646dc1d3ae3a0be345eae8f363bfcf964615b8c7c06e9623f05be858 |
  | [BRI] | processes/electromagnetic/standard/src/G4BraggIonModel.cc | 77fd6989ae2ea79fe124781a8c6d918e57763d7a4665a8fe8fe073d19c31f1dd |
  | [ETU] | processes/electromagnetic/utils/src/G4EmTableUtil.cc | 41a2f7ba0a8739e36028dbe036cb15b063bb318aa07e8037a3b8307e3e30eadc |
  | [EMC] | processes/electromagnetic/utils/src/G4EmCorrections.cc | a357d07617816f4e3ef1c7558890cf4966ac55bf85ba4ad2392713b10b52d3c7 |
  | [LSh] | processes/electromagnetic/standard/include/G4LindhardSorensenIonModel.hh | f888df1a2e5e93556fa80e0d0e92d758c11310a13eec35c1760a824d4fec4fe6 |
  | [VEM] | processes/electromagnetic/utils/src/G4VEmModel.cc | 623d8fbb92c6b3b5d7619baa603a6e7061dd9e80147f0502036cc117c58ffb51 |
- Model arithmetic: `validation/plans/v3-008/p13-scripts/p13.py` (ionmc at worktree HEAD 2aaf1b67; model-only).

## 2. The algorithm as executed (FACT unless marked)
1. Pre-step quantities, once per step in PostStepGPIL: preStepScaledEnergy = E_pre * massRatio [VELP] L596, massRatio =
   M_GenericIon/M_ion set per track [VELP] L546; q2 = ChargeSquareRatio(track) at the **pre-step kinetic energy**
   [VELP] L609 (-> [LS] L135 `corr->EffectiveChargeSquareRatio(p, mat, kinEnergy)`); fFactor = q2, reduceFactor =
   1/(fFactor*massRatio) [VELP] L612.
2. Range: fRange = reduceFactor * R_base(E_s) from the spline log-vector [VELPh] L646 (cached per energy/couple, L643).
3. AlongStepDoIt: full stop only if length >= fRange [VELP] L797; short-step estimate eloss = length * fFactor *
   S_base(E_s) [VELP] L813 with [VELPh] L602; **branch**: if eloss > E_pre * linLossLimit [VELP] L824 then
   x = (fRange - length)/reduceFactor [VELP] L826 and eloss = E_pre - E_base^-1(x)/massRatio [VELP] L827 via the inverse
   range vector [VELPh] L681-692 (ScaledKinEnergyForLoss). linLossLimit = 0.02 (D-4 section 3.2, [II] L80; runtime
   printout `linLossLim= 0.02`, D-4 section 2.6).
   => the table part is the exact solution (up to table error) of dE/ds = -q2(E_pre) * S_base(E * massRatio):
   **q2 frozen at the pre-step energy, base stopping integrated along the step.**
4. Ions only: CorrectionsAlongStep(mat, p, **preStepKinEnergy**, cut, length, eloss) [VELP] L856, after the table loss.
   [LS] L259 return if eloss >= E_pre; L263 e = max(E_pre - eloss/2, E_pre/2) (mid-step by the table loss); L267 q2 at e;
   above escaled 2 MeV: L296 deltaL0 = 2*Barkas(e)*(z-1)/z, L297 deltaL = LS(Z, gamma(e)), L299-300
   res = eloss + 2pi mc2 re2 * q2(e) * n_e * (deltaL + deltaL0) * length / beta2(e); L308 guard (res > E_pre or
   2 res < eloss -> unchanged). So: table part with q2(E_pre); additive higher-order term by a one-point midpoint rule
   with q2 and beta at the mid-step energy. No q2 re-scaling of the table part for GenericIon (contrast [BB] L346-L354:
   BB rescales by q2(mid)/q2(pre) only if eloss >= 5 % E_pre; BB returns for alpha, [BB] L343).
5. Table content: LS table for GenericIon = [BR]/[BB] ICRU90 proton array at E_s (x q2 = 1 for charge 1, [IEC] L95-96
   `effCharge <= 1.5`) on both sides of 2 MeV ([LS] L240-242; [BR] L644-647 ICRU90 for water below; [BB] L273-293 above):
   continuous at the 2 MeV switch for G4_WATER (cut = 1 km >= tmax, so no cut term).
6. Tables: 20 bins/decade [OPT4] L119 (default 7, [EMP] L168), 100 eV-10 GeV = 160 bins (printout, D-4 2.6); dE/dx log
   vector with spline (`spline = true` [VELPh] L516, passed at [VELP] L430-432); range by midpoint quadrature with 100
   sub-intervals per bin over the spline dE/dx [LTB] L188-250; inverse range = free vector (x = range, y = energy) with
   spline [LTB] L263-285 (`splineFlag = true` [LTBh] L112, no SetSplineFlag in [LTM], built at [LTM] L878-879);
   spline = not-a-knot cubic in **linear** x and y ([PVc] L228-290, L312-372; evaluation [PVi] L125-147).
7. q_eff in Geant4 ([IEC]): bare above E_s > z * 20 MeV ([IEC] L66, L113), i.e. above ~119 (C) / ~159 (O) MeV/u;
   below, a Brandt-Kitagawa/Ziegler form for Z > 2 and a Ziegler He form for Z = 2.

## 3. Exact expression and expansion
Notation per point: E = E_pre (total kinetic), S(E) = q2(E) S_b(E) (table part), C(E) = higher-order additive term,
F = S + C, ell = path. Exact: dE/ds = -F. Algorithm: Delta_alg = Delta_tab + C(E - Delta_tab/2) * ell with
Delta_tab solving dE/ds = -q2(E) S_b, q2 frozen at E. Expansion (E'' = F F', E''' = -(F''F^2 + F'^2 F)):
Delta = F ell - F F' ell^2/2 + (F'' F^2 + F'^2 F) ell^3/6 + ..., same for Delta_tab with G = q2(E_pre) S_b.
With f = ell F/E (loss fraction; the generator uses ell S_M(E0)/E0 = 3 %), lambda = dln q2/dlnE, sigma = dlnS/dlnE:
- g(f) = Delta_alg/Delta_exact - 1 = a f + b f^2 + O(f^3), a = lambda/2 * (S/F) (the C-term midpoint is exact to O(f));
  b collects lambda*sigma, lambda^2, dlambda/dlnE and the midpoint-rule curvature of C (order (C/F) f^2/24, ~1e-5 rel).
- Richardson residual: R = (f2 g(f1) - f1 g(f2))/(f2 - f1) = -b f1 f2 + O(f^3); f1 f2 = 1.44e-3.
- a and b obtained numerically from g at f = 0.5 % and 1 %; g(f1), g(f2), R by RK4 (400 steps) on ionmc's default model
  (Bethe, I = 78 eV, Barkas-125 q, shell+Barkas+Bloch+Mott). Two emulations: (A) "frozen_all" = the committed evaluator's
  emulation (whole S scaled by q2(E_pre)/q2(E)); (B) "g4split" = section 2 item 4 literally: leading Bethe part with
  frozen q2, higher-order part (Barkas+Bloch+Mott of ionmc as proxy for deltaL + deltaL0) at the table mid-step energy.

## 4. Results per point (ESTIMATE; ionmc model; % values)
| ion E0 [MeV/u] | dlnq2/dlnE | dlnS/dlnE | a (B) | b (B) | g(f1) B | g(f2) B | R (A) | R (B) |
|---|---|---|---|---|---|---|---|---|
| C 2.4 | +0.192 | -0.536 | +0.116 | +0.152 | +0.360 | +0.590 | (A not printed; ~ -0.018) | -0.0240 |
| C 3   | +0.153 | -0.587 | +0.093 | +0.128 | +0.288 | +0.473 | -0.0158 | -0.0202 |
| C 5   | +0.080 | -0.691 | +0.050 | +0.076 | +0.155 | +0.256 | -0.0094 | -0.0121 |
| C 10  | +0.022 | -0.767 | +0.016 | +0.029 | +0.049 | +0.082 | -0.0032 | -0.0046 |
| C 50  | +0.000 | -0.789 | -0.000 | +0.001 | -0.000 | -0.001 | -0.0000 | -0.0001 |
| O 2.4 | +0.266 | -0.454 | +0.168 | +0.209 | +0.520 | +0.851 | -0.0241 | **-0.0328** |
| O 3   | +0.222 | -0.507 | +0.140 | +0.183 | +0.434 | +0.712 | -0.0211 | -0.0287 |
| O 5   | +0.132 | -0.627 | +0.085 | +0.122 | +0.264 | +0.434 | -0.0142 | -0.0192 |
| O 10  | +0.048 | -0.732 | +0.033 | +0.056 | +0.104 | +0.171 | -0.0062 | -0.0089 |
| O 50  | +0.001 | -0.786 | +0.001 | +0.002 | +0.002 | +0.003 | -0.0001 | -0.0004 |
| He 10 | +0.000 | -0.802 | -0.000 | +0.000 | -0.000 | -0.000 | -0.0000 | -0.0000 |
| He 50 | +0.000 | -0.792 | -0.000 | -0.000 | -0.001 | -0.002 | -0.0000 | +0.0001 |
(7, 15, 20, 30 MeV/u computed as well; all |R| <= 0.014 %.) Checks: a ~ lambda/2 (e.g. O 2.4: lambda/2 = 0.133 vs a(A)
0.134); R ~ -b f1 f2 (O 2.4: -0.0302 % vs -0.0328 %, the rest is O(f^3)). Max |R| = 0.033 % (O 2.4, emulation B);
emulation A (the committed evaluator) underestimates R and Delta_pred = g(f2) - g(f1) by ~25-35 % at 2.4-5 MeV/u
(O 2.4: Delta_pred 0.262 % (A) vs 0.332 % (B)). Matches SR-3's -0.019/-0.024 % for emulation A at 2.4 MeV/u.

## 5. Neglected contributions, bounded separately
1. **Effective-charge model (ESTIMATE).** Geant4 uses [IEC] (Ziegler/Brandt-Kitagawa), ionmc Barkas-125; R scales
   ~ with lambda*(sigma, lambda). Allowance: factor 2 on the q-driven part -> |R_q| <= 0.066 % (O 2.4). OPEN: evaluate
   lambda of the [IEC] formula for water (needs G4_WATER Fermi energy / Zeff from G4IonisParamMat; not transcribed here).
2. **Mid-step correction evaluation (FACT mechanism, ESTIMATE size).** Included in emulation B; the one-point midpoint
   error is O((C/F) f^2 dln^2C/dlnE^2 / 24) < 1e-5 relative; its R contribution < 0.001 %. The [LS] L263 mid-step uses
   the table loss (not the corrected loss): shift of the evaluation energy by (C/F) f/2 ~ 1e-3 f -> < 0.001 % in R.
3. **Table interpolation (FACT structure, ESTIMATE size) - the largest non-physics term.** Emulated in p13.py: G4 node grid
   (100 eV * 10^(k/20)), G4 not-a-knot spline ported from [PVc] L312-372 / [PVi] L125-147, range by the [LTB] L230-249
   quadrature, inverse range spline in linear (R, E); smooth proxy S_b = ionmc proton Bethe (>= 1 MeV, power-law below).
   Single-step relative errors eps(f1) = -0.006...-0.014 %, eps(f2) = -0.002...-0.008 %; the extrapolation amplifies
   them (errors of fixed absolute size scale like 1/f; Richardson weights f2/f1/(f2-f1) etc.): |R_int| <= 0.028 % at the
   nine C/O energies (scaled energy E0 * 1.0073), <= 0.028 % over a 120-point scan 2.3-60 MeV (95th pct 0.026 %).
   Mechanism: inverse-range spline of E ~ R^0.57 with node ratio ~1.22 in R; bound h^4 scaling, 5/384 h^4|E''''| ~ 3e-5 E.
   The real ICRU90 proton array (itself spline-interpolated data, [I90] L164-169) may carry more high-order structure
   than the proxy: allowance factor 3 -> |R_int| <= 0.085 %. 2 MeV switch: continuous for water (section 2 item 5), no kink.
4. **msc path-length excess (ESTIMATE).** Energy loss uses the true step length (step.GetStepLength(), [VELP] L~773).
   Urban msc (P2 F1) converts z -> t ~ -lambda1 ln(1 - z/lambda1): t/z - 1 ~ z/(2 lambda1) ~ <theta^2>_exit/4; Highland
   estimate for C 2.4 MeV/u (slab ~2 um, pbeta ~ 58 MeV, theta0 ~ 2-3 mrad) gives ~5e-6; it is linear in f, so the
   extrapolation removes it; residual < 1e-5 %. Lateral displacement does not change the energy loss of a one-step
   traversal. OPEN (not read): [URB] ComputeTrueStepLength / tlimitminfix2 lines and the ion step-limit type
   (P2 F1 INFERRED fMinimal); not needed for the bound because the step-count smokes enforce one step.
5. **Energy convention (FACT, negligible).** E_s = E * M_GenericIon/M_ion ([VELP] L546, L596) = same-velocity scaling;
   E0 (MeV/u) vs E_s differ by 0.73-0.76 %; g depends on log-slopes only, so the shift changes R by < 1 % of itself.
   Velocity-mapping differences between ionmc and Geant4 enter rho (the stopping ratio), not the stepping bias.
6. **Branch condition (FACT).** The test is length * q2(E_pre) * S_base(E_s) > 0.02 E_pre ([VELP] L813, L824), i.e. the
   pre-step **table** loss without the LS term. The design sets ell S_M(E0) = 3 % E0; Geant4's table S differs from S_M by
   a few % (I3 residuals), so the estimate is 2.8-3.2 % >= 2 % at f1 (margin ~1/3); f2 a fortiori. The short-step
   branch (forward Euler on the pre-step dE/dx, q² per step) biases by about σ·δ/2 per step with σ = dlnS/dlnE, i.e. −0.04 to −0.06 % for the revision-5 fine-step runs (δ ≈ 0.15 %; SR-3 fourth pass), not −a f; it is not taken by the one-step slabs. Step-count smokes verify the single step; OPEN: they do not log the
   branch; no direct runtime evidence beyond the arithmetic.
7. **Energy dependence of the stopping ratio itself (ESTIMATE).** r(f) also contains rho(f) = <S_G4>/<S_M> - 1 over the
   slab; its quadratic part ~ d^2 rho/dlnE^2 f^2 / 6 with |d^2 rho/dlnE^2| <~ 0.02 -> < 0.001 %; not stepping.
8. **Not stepping, out of scope:** mean preservation of the fluctuation models (covered by the fluctuation-off smoke /
   control (iii)); the 1 keV/fRange end-of-track rules ([VELP] L797, L908-911) are never reached in a 5 % slab.

## 6. Bound on the residual after extrapolation (per point, all C/O/He points)
| contribution | best estimate | allowance (bound) | basis |
|---|---|---|---|
| q2 frozen + LS midpoint (emulation B) | <= 0.033 % | 0.066 % (x2 for q_eff model) | FACT mechanism ([VELP] L609-612, L824-827, [LS] L263-300); ESTIMATE size |
| table interpolation | <= 0.028 % | 0.085 % (x3 for real data) | FACT structure ([LTB], [PVc], [PVi], [OPT4] L119); ESTIMATE size |
| midpoint curvature, msc, convention, rho curvature | < 0.003 % | 0.005 % | ESTIMATE |
| **linear sum** | **<= 0.06 %** | **<= 0.16 %** | |
Signs: at 2.4 MeV/u both leading terms are negative (R_B -0.033 %, R_int -0.028 %): best-estimate sum -0.06 %.
**Conclusion (ESTIMATE):** the residual stepping bias after extrapolation is <= 0.1 % at the best-estimate level at every
point (max 0.06 %, O 2.4 MeV/u) and <= 0.2 % with the stated allowances (0.16 %). The <= 0.1 % premise holds only with
modest margin at 2.4-3 MeV/u if the allowances are taken at face value; the run-time check on R_model alone (emulation A)
covers only ~0.024 % of the ~0.16 % allowance.
Source-verified (FACT): the algorithm in section 2 (q2 at the pre-step energy, range-table inversion with frozen q2,
branch condition and its quantity, LS term at mid-step with q2(mid), alpha without along-step correction, spline and
grid settings, continuity of the base table at 2 MeV in water). Model estimates (ESTIMATE): every number in sections
4-6 (ionmc q_eff and slopes as proxy for Geant4's; proxy proton stopping for the interpolation emulation; Highland msc).

## 7. OPEN items and the discriminating evidence
- O1: Geant4 [IEC] dlnq2/dlnE for C/O in water (replaces the x2 allowance); needs G4IonisParamMat Fermi energy/Zeff. **ESTIMATE (SR-3 fifth pass, 2026-11-16): the Ziegler-type Geant4 charge ([IEC], v_F ≈ 0.97) against ionmc's Barkas-125 form gives dln q²/dlnE 0.190 vs 0.192 (C) and 0.258 vs 0.266 (O) at 2.4 MeV/u, |Δλ| ≤ 0.01; the revision-6 gate e₁ ≈ 0.0155·Δλ fails at |Δλ| ≈ 0.1.**
- O2: interpolation with the actual ICRU90 proton array (replaces the x3 allowance); needs the data values, which this
  record does not transcribe; alternatively a one-step vs many-step Geant4 comparison (MaxStep 0.1 um in the slab) at
  2.4 and 3 MeV/u would measure (frozen-q2 + interpolation) jointly at f1, f2 (statistics: SE of Delta E <~ 0.02 %).
- O3: [URB] ComputeTrueStepLength lines and the GenericIon msc step-limit type (bound already negligible).
- O4: emulation A (committed evaluator) vs B: B implements [LS] L263-300 literally; consider adopting B for Delta_pred
  and R_model, or widen the stated budget; decision for the orchestrator.
- O5: register the fetched sources (hashes above) before any of their line citations enter Git. **Status 2026-11-16 (orchestrator): NOT met — the citations are in Git bound to the full sha256 listed in §1 without a provenance-cache entry (the acquisition tool errored); rebinding pending.**
