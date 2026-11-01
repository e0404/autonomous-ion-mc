# P3b (V3-008): default of the LS option — charge convention of the higher-order stopping terms, record of 2026-11-01

Status: pre-freeze physics assessment, SEEN. It follows the P3 re-opening (`p3-lindhard-sorensen.md`, section "Mixed charge convention", Codex REVIEW-e0b586fd finding 1) and was written on 2026-11-01 before the acceptance-plan freeze. It is the record behind the decision 0043 §1 "Amendment (2026-11-01) — LS option" and the closing of Codex REVIEW-6a14f4ec finding 1 (round V3-008-5, `review.md`). It is not qualification evidence. No I3, V1 or other reference value was used. The text is the planner's bounded assessment (session scratch `v3-008/p3b/ls-default-assessment.md`), copied with its sources, numbers and open items. Changes against the scratch note: the option field name follows the implementation (`BetheOptions.ls`); the decomposition table gains the 3, 16 (C) / 3, 24 (O) MeV/u rows and the z²-weighted column from the same script output; the Geant4 line references quote the source lines; the outcome is stated in section 6.

Inputs:
- the P3 record (`validation/plans/v3-008/p3-lindhard-sorensen.md`, mixed-convention section) and `p3-output.json`;
- decision 0043 §1, the P2 record and RQ-S1;
- `src/ionmc/physics/stopping.py` (at ea220a88, before the option was added);
- `validation/scripts/ion_em/lindhard_sorensen.py`.

Markers:
- **SAYS**: what the source states (read in this pass).
- **CODE**: read in source code (untrusted files, read with grep/sed only, never compiled or executed).
- **INF**: inference.
- **CALC**: scratch arithmetic. Scripts: `v3-008/p3b/scripts/conv.py` (sha256 1467711fe4e2997370debad59cf3c2512c35ccb0dc290bd6239d9af5f9b964c1) and `scripts/rng.py` (d4386b3233def2984d68cd23e38b10e9cd008e37bcb09458c8934beea80245b4). Output: `conv-output.txt` (a19f62e826717898d33425656b14c7b31ade8ab8c585978d1cb6699619130ba7). Model: ionmc, I = 78 eV. Aggregates only. The scripts are session scratch and are not committed. They use the committed `lindhard_sorensen.py` functions and the ionmc Bloch/Mott forms.

## 0. What the trigger actually measured (CALC, decomposition)

r_mixed = [LS(z) − Bloch(z) − Mott1(z)] + [Bloch(z) − Bloch(q)] + [Mott1(z) − Mott1(q)]. All values in % of S:

| Ion, E (MeV/u) | q/z ionmc (Barkas-125) | q/z ZBL (Geant4 formula, vF = 1 assumed) | LS−B−M at z | Bloch(z) − Bloch(q) | Mott(z) − Mott(q) | r_mixed | r_eff | Bloch(q) itself | z²-weighted LS(z) (§1) |
|---|---|---|---|---|---|---|---|---|---|
| C 2 | 0.916 | 0.912 | +0.05 | −1.29 | +0.01 | −1.22 | +0.05 | −9.1 | −3.18 |
| C 3 | 0.952 | 0.947 | +0.06 | −0.52 | +0.01 | −0.46 | +0.05 | −6.3 | −1.15 |
| C 5 | 0.980 | 0.976 | +0.06 | −0.14 | +0.00 | −0.07 | +0.06 | −3.8 | −0.23 |
| C 10 | 0.996 | 0.994 | +0.06 | −0.01 | +0.00 | +0.05 | +0.06 | −1.8 | +0.03 |
| C 16 | 0.999 | 0.998 | +0.06 | −0.00 | +0.00 | +0.06 | +0.06 | −1.1 | +0.06 |
| O 2 | 0.871 | 0.868 | +0.06 | −2.64 | +0.02 | −2.56 | +0.06 | −13.1 | −7.52 |
| O 3 | 0.918 | 0.914 | +0.08 | −1.25 | +0.02 | −1.16 | +0.07 | −9.5 | −3.12 |
| O 5 | 0.960 | 0.955 | +0.09 | −0.41 | +0.01 | −0.31 | +0.08 | −6.0 | −0.83 |
| O 10 | 0.989 | 0.987 | +0.10 | −0.06 | +0.00 | +0.04 | +0.09 | −3.1 | −0.02 |
| O 24 | 0.9991 | 0.9985 | +0.10 | −0.00 | +0.00 | +0.09 | +0.10 | −1.2 | +0.09 |

Conclusions (INF from CALC):
- The part that is genuinely "LS beyond Bloch + Mott1" is +0.05 to +0.10 % of S at every energy and for every charge choice.
- The trigger fired only through Bloch(z) − Bloch(q). The open question is the **charge argument of the Bloch term**, which is a large correction (−9 % (C) / −13 % (O) of S at 2 MeV/u). It is not LS physics.
- Geant4's ZBL effective charge is within about 0.5 % of ionmc's Barkas-125 q at 2–10 MeV/u (vF = 1 assumed; the Fermi velocity of water in Geant4 was not read). So the spread between the mixed and the all-at-q conventions is the same size in Geant4 as in ionmc.
- A fully "bare close-collision" treatment would weight the close-collision correction with z², not q² (§1). In ionmc normalisation that term is (z/q)² LS(z) − Bloch(q) − Mott1(q):
  - −3.2 % (C) / −7.5 % (O) of S at 2 MeV/u;
  - −0.2 % / −0.8 % at 5 MeV/u;
  - below 0.1 % in magnitude above 10 MeV/u.
- The plausible convention band at 2 MeV/u is therefore about +0.05 to −3.2 % (C) and +0.06 to −7.5 % (O) of S. All variants coincide above 16 (C) / 24 (O) MeV/u.
- Region concerned:
  - ionmc CSDA range from 2 to 16 MeV/u (C) is 0.97 mm, and from 2 to 24 MeV/u (O) 1.53 mm;
  - to 10 MeV/u it is 0.41 / 0.32 mm;
  - the range effect of the mixed convention is about 0.4 / 0.9 µm absolute (P3).

## 1. Physics of the charge convention

### What the sources say (SAYS)

- Sigmund & Schinner, "The Bloch correction, key to heavy-ion stopping", J. Appl. Phys. 128, 100903 (2020), https://doi.org/10.1063/5.0015478 (read via the AIP page):
  - the Bloch correction "picks up contributions from close collisions". It is "a property of Coulomb scattering" and is written with the atomic number Z1 (bare-ion curves);
  - partially stripped ions are treated (in PASS) "by replacing the Coulomb interaction with the target electrons by a pertinent screened-Coulomb potential". Frozen-charge results show that "the dependence on the ion charge q1 ... is much weaker" than on Z1;
  - against the effective charge: "despite its name, the effective charge is not a charge". Replacing Z1 by q1 "would predict a vanishing stopping cross section for a neutral ion". The effective-charge factor "primarily expresses the transition from the Bethe to the Bohr stopping regime", i.e. it has traditionally absorbed the Bloch/Bohr (higher-order) physics near the Bragg peak;
  - the article gives no prescription "Bloch at q" or "Bloch at z" for use inside a Bethe model with an effective charge.
- Sigmund & Schinner, Eur. Phys. J. D 12, 425 (2000), https://doi.org/10.1007/s100530070004 (abstract only):
  - the binary theory uses an effective potential that "holds equally well for dressed as for stripped ions" and avoids a formal division into close and distant collisions;
  - PASS is the code behind ICRU 73 (secondary statements in search results; the ICRU 73 report itself was not read; icru.org/?p=88 gives only the scope, 0.025–1000 MeV/u and Z1 = 3–18).
- Lindhard & Sørensen, Phys. Rev. A 53, 2443 (1996): a theory for a point (or finite) **bare** nucleus. It was not read directly (as in P3). The formula used rests on Kats et al., arXiv:2008.04819, and MARS15, arXiv:physics/0505043.
- ZBL and Hubert–Bimbot–Gauvin are empirical effective-charge schemes, S_ion = (γ Z1)² S_ref at equal velocity, fitted to stopping data. The Z1³/Z1⁴ effects of the ion are absorbed into the fitted charge fraction. This is SAYS for ZBL, from the Geant4 source comments citing ZBL 1985 and Ziegler–Manoyan 1988. HBG 1990 (ADNDT 46, 1) was not read in this pass; INF from its known form.

### Inference (INF)

- In a screened-ion picture, collisions with an impact parameter inside the projectile's electron cloud see nearly the bare charge z, and distant ones see the net charge (the Brandt–Kitagawa idea).
  - The Bloch, LS and Mott corrections live at impact parameters of order max(ħ/mv, z e²/mv²). For C at 2 MeV/u that is about 0.04–0.06 Å, inside the projectile K-shell radius (about 0.09 Å).
  - On this argument the close-collision correction should use a charge near z, with weight z².
  - That is neither the ionmc convention (q, q²) nor the Geant4/catima convention (z in the correction, q² weight).
- But ionmc's q (Barkas 1963, coefficient 125) is an **empirical stopping effective charge**, not a mean ionic charge.
  - Like ZBL's, it was fitted in a framework where everything is scaled by q².
  - Putting a bare-charge Bloch term (or a z²-weighted one) on top of an empirical q risks counting the low-velocity higher-order physics twice.
  - "All at q" is the convention in which the effective-charge parameterisations were built. It is internally consistent, but it is not a theory of the screened ion either.
- A theory-based answer needs the screened-ion problem solved directly: a PASS-type calculation (screened-Coulomb potential plus charge-state averaging) or a partial-wave transport cross section for a screened potential. Neither the q nor the z convention is a derived result.
- **This cannot be decided from theory alone at the 1–3 % level.** It sits inside the 2–5 MeV/u model uncertainty already declared: the shell correction held at η = 0.13, a large Barkas term, and effective-charge alternatives that differ by a few % in q² (decision 0043 §1).

## 2. What Geant4 does (CODE; v11.4.2 from the GitHub tag; untrusted, read with grep/sed)

Research downloads to the session scratchpad (`v3-008/p3b/dl-*`, not committed, not acquired through the reference-data tool; role lineage, exploratory). sha256 of each file as downloaded:

| File (tag v11.4.2) | sha256 | Note |
|---|---|---|
| `G4LindhardSorensenIonModel.cc` | e16ed6cd1d132183684ef47948b6f8fe36526a48663a6dc1878d5bb8330eb12c | identical to the pinned P2 cache file [LS] |
| `G4EmCorrections.cc` | a357d07617816f4e3ef1c7558890cf4966ac55bf85ba4ad2392713b10b52d3c7 | |
| `G4BetheBlochModel.cc` | 5dc2a9002bff25f239a6402d3976a1f285f071abf4b41585b76f275e4539e01a | |
| `G4ionEffectiveCharge.cc` | 891854a97b97f038f9a0fa88b9f2f4ed139ea080e8edf9bd349cfd93031ae4c5 | |
| `G4LindhardSorensenData.cc` | 70ae7c6ba3ef990b7005240bd2b60e4040a3fb60bb1c9e4dbf829c109180eda9 | identical to the v11.2.0 SEEN download of P3 §7 |

- `G4LindhardSorensenIonModel.cc`:
  - L157–158: `charge = particle->GetPDGCharge()*inveplus;` and `Zin = G4lrint(std::abs(charge));`, i.e. the nominal Z;
  - L135 and L267: `q2 = corr->EffectiveChargeSquareRatio(...)`;
  - L296: `deltaL0 = 2.0*corr->BarkasCorrection(p, mat, e)*(charge-1.)/charge;`, with the nominal charge in the ratio;
  - L297: `deltaL = lsdata->GetDeltaL(Zin, gam);`, the LS table at the nominal Z;
  - L299–300: `res = eloss + CLHEP::twopi_mc2_rcl2*q2*eDensity*(deltaL+deltaL0)*length/beta2;`, so both terms carry the q2 weight;
  - applied for escaled > fElimit = 2 MeV (L82, L271);
  - no comment in the file explains the choice of Zin.
- `G4EmCorrections.cc` is the other Geant4 ion path (`G4BetheBlochModel` L323 `dedx += corr->HighOrderCorrections(p,material,kineticEnergy,cutEnergy);`; `ComputeIonCorrections` / `IonHighOrderCorrections` L218–283):
  - L158–160: `charge = p->GetPDGCharge()*inveplus; if(charge > 1.5) { charge = effCharge.EffectiveCharge(p,mat,kinEnergy); } q2 = charge*charge;`;
  - L695: Bloch with `G4double y2 = q2/ba2;`;
  - L718: Mott `return CLHEP::pi*CLHEP::fine_structure_const*beta*charge;`;
  - in this path **Bloch and Mott are evaluated at the effective charge**, which is the ionmc convention;
  - L250–283: a matching term (`rest/e` at the ICRU73 transition energy) removes the junction offset.
- `G4ionEffectiveCharge.cc`:
  - L104–180: the ZBL fractional charge with a Brandt–Kitagawa-type `chargeCorrection`;
  - L66 (`energyHighLimit = 20.0*CLHEP::MeV`) and L113: the bare charge is returned above reduced energy Z × 20 MeV.
- `G4LindhardSorensenData.cc` at tag v11.4.2 has sha256 70ae7c6b… (full hash in the table above). That is **identical to the v11.2.0 SEEN download** used in P3 §7. It closes the P3 §7 item "the v11.2.0 vs v11.4.2 data files are not shown to be identical", on the assumption that the GitHub tag equals the release. The LS-model file from the same source matched the pinned cache hash.
- Summary (INF): Geant4 is not self-consistent across its own ion models.
  - The LS path is "mixed": LS at Z, a q² weight, applied from about 2 MeV/u.
  - The Bethe–Bloch ion path is "all at q".
  - The likely reason is engineering: one (Z, γ) table indexed by an integer Z, with the charge-state effect taken by the existing q² prefactor.
  - No published rationale was found.
- ATIMA lineage (catima, AGPL, read only; `calculations.cpp` sha256 ff08e35799c06cc7b5d0a6c95a73171482c0aac3be6120200ae059f35a3ecbe5, `catima.cpp` 2ef012e5d38a1f7da1fc5e4c04effda1eaa16d0bafdeedc45b33254efec7bcf7):
  - `bethek_dedx_e` multiplies `(f2*barkas + LS - delta/2)` by `zp_eff^2`, with `precalculated_lindhard(p)` at the bare p.Z: also mixed;
  - catima.cpp L33–43: this is used only above 10 MeV/u. 10–30 MeV/u is blended with the SRIM-type `sezi_dedx_e`, and pure LS-Bethe applies above 30 MeV/u;
  - so in that lineage the mixed convention is not applied where q/z deviates strongly from 1.

## 3. Recommendation (as assessed)

**Add the option, default `off`.** `off` is the current model: Bloch and Mott1 at q, no LS. It equals `eff` within 0.1 % of S everywhere. Values:
- `off` (default);
- `eff`: ΔL_LS(q) replaces Bloch(q) + Mott1(q);
- `bare_mixed`: ΔL_LS(z) replaces Bloch(q) + Mott1(q). This is the Geant4/catima-lineage convention, for sensitivity and attribution only.

For z = 1 the option is a no-op, so the qualified proton model is unchanged.

Reasons:
1. The trigger fired on the charge argument of the Bloch term (§0), not on omitted LS physics. Within ionmc's own convention P3 shows r ≤ 0.1 % of S and ≤ 0.095 % of R.
2. None of the conventions is derived theory for a partially stripped ion (§1).
   - `bare_mixed` is a hybrid: z in the correction with a q² weight.
   - A consistent bare-close-collision treatment would be about 2.6 times (C) / 2.9 times (O) larger at 2 MeV/u, with the same sign.
   - The all-at-q convention is the one in which empirical effective charges (Barkas, ZBL, HBG) were constructed, and the one Geant4 itself uses in its Bethe–Bloch ion path.
3. The convention spread lies inside the already declared low-energy model uncertainty: ≤ 1.2 % (C) / 2.6 % (O) of S at 2 MeV/u, < 0.1 % above 4.7 / 6.6 MeV/u. It acts on the last 1–1.5 mm of range, with an absolute effect ≤ 1 µm.
4. Adopting Geant4's convention as the default would remove a difference from the TOPAS reference by construction. That reduces V1 independence without a physics reason.
5. `eff` as default is defensible.
   - It is more complete at high energy and shortens the range by 0.05–0.08 % at the standard points.
   - It costs a production LS table generator in `src`, with long-double partial-wave sums and tests, before freeze, for a change ≤ 0.1 %.
   - It is a second choice, not the recommendation.

Uncertainty that theory cannot resolve: the correct charge weighting of the close-collision (Bloch/LS) correction for C/O at 2–16/24 MeV/u. The band is about +0.05 to −3 % (C) / −7.5 % (O) of S at 2 MeV/u, and it falls below 0.1 % by about 10 MeV/u. Discriminating evidence:
- (a) a screened-ion calculation. Either the PASS/DPASS program (free from SDU; availability not checked), or a project-side partial-wave transport cross section for a screened C/O ion (Hartree–Fock or Thomas–Fermi 1s cloud, charge-state weighted; non-relativistic suffices below 25 MeV/u);
- (b) measured C/O stopping in low-Z targets at 2–10 MeV/u from the IAEA database. No water data exist, and shell and Barkas effects confound the comparison.

**Do not** use the libdEdx ICRU73_NEW I3 rows to choose: that would turn the I3 evaluation into calibration.

## 4. Consequences for I3, V1 and range at the standard points

- Standard points (50–430 MeV/u): q/z > 0.999, so the convention choice changes the range by ≤ 0.0003 % (P3: 0.0456 vs 0.0458 % for C at 430). Only "LS or not" matters there.
  - `eff` or `bare_mixed` vs `off` shortens the range by 0.046 % (C, 0.14 mm of 309 mm) and 0.080 % (O, 0.19 mm of 232 mm).
  - Against the V1 R80 tolerance of 0.35 % (≈ 1.08 mm C, 0.81 mm O), that is 13 % (C) / 23 % (O) of the tolerance, and below 2 % of it for He (0.005 %).
  - This comes from LS on/off, not from the charge convention. It is not decisive. It belongs in the V1 attribution budget, together with the unsettled possibility (P3 §7) that Geant4 applies half the LS term.
- I3 (libdEdx ICRU73_NEW):
  - 2–10 MeV/u, 6 % tolerance: the default choice moves ionmc S by up to 1.2 % (C) / 2.6 % (O) at 2 MeV/u, i.e. up to 20 % / 43 % of the tolerance at the lowest points. That falls to 0.07 % / 0.31 % at 5 MeV/u. With `off`, this spread is a declared model uncertainty, not a used-up tolerance.
  - ≥ 10 MeV/u, 1.5 %: ≤ 0.1 % either way.
  - The proton I3/V1 rows are unaffected because the option acts only for z ≥ 2. For protons q = z and LS − Bloch ≤ 0.093 % of S.
- Transport observables at clinical energies are unaffected by the convention.
  - The difference sits in the last ≤ 1.5 mm of range, the part below the energies where the variants coincide (16 MeV/u for C, 24 MeV/u for O), and amounts to about 1 µm in absolute range.
  - Its effect on the Bragg-peak shape at ≥ 1 mm scoring resolution is negligible (INF).

## 5. Not verified

- Lindhard & Sørensen 1996, ICRU 73 (its text, PASS settings and any LS use at high energy), the 2009 ICRU 73 erratum, HBG 1990 and ZBL 1985 were not read. Statements about them here come from secondary sources or from the known form of the formulas.
- The PASS treatment of the inverse-Bloch term for screened ions (which charge it uses) was not confirmed from a primary text.
- No published comparison quantifying the q-vs-z convention for C/O in water at 2–10 MeV/u was found. Geant4 validation papers (Lechner et al. 2010; Ivanchenko et al.) were not located or read in this pass.
- The water Fermi velocity in Geant4 (`GetFermiEnergy`) was not read. The ZBL q/z column assumes vF = 1 v0.
- It was not settled whether the G4BetheBlochModel base table for GenericIon already contains the charge-1 Bloch/Mott terms (L323) or the ICRU 90 proton array (earlier branch, ~L290–300). A possible small double count in Geant4 is unquantified.
- GitHub-tag fidelity to the Geant4 11.4.2 release tarball is assumed. One matching hash supports it.
- DPASS availability and licence were not checked.

## 6. Outcome (2026-11-01; orchestrator decision)

- The option is implemented in the same commit that adds this record: `BetheOptions.ls ∈ {"off", "eff", "bare_mixed"}`, acting for z ≥ 2 only, default `"off"`.
  - Implementation: `src/ionmc/physics/lindhard_sorensen.py` and `src/ionmc/physics/stopping.py`.
  - Tests: `tests/ionmc/test_ion_em_ls_option.py`.
  - The sha256 values are bound at commit in the acceptance plan rows P3 and P9.
- The default is set by the decision 0043 §1 amendment "Amendment (2026-11-01) — LS option", on the reasons in §3.
- I3b(ii) reports the three variants (report-only). I3c may use `bare_mixed` for the TOPAS-attribution diagnostic (exploratory).
- No I3, V1 or other reference value was seen or used for this choice.
