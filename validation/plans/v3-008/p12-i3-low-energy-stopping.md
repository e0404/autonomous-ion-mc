# P12: research note on the row-I3 deficit (C/O stopping below 50 MeV/u vs ICRU 73 revised), 2026-11-15

Status: research record written by the physics researcher at the orchestrator's request after row I3 FAILED (plan results
note "rows I2 and I3", 2026-11-15). It is NOT evidence and changes no grade; it is the input to a future prospective
amendment (plan D). Labels FACT / ESTIMATE / OPEN as stated by the author; estimates use only the committed ionmc model.
Per-energy deviations against the reference table are given as band values only (stopping-results policy); the model-only
scripts of the estimates are committed under `validation/plans/v3-008/p12-scripts/` (`decomp.py`, `transport.py`; they read no
reference table and need only the committed ionmc model); the band diagnostic against the table itself is NOT committed
(a hypothesis until the follow-up amendment commits a tested reproducer). Not a clinical claim.

# V3-008 row I3: explaining the 2–50 MeV/u deficit (research note, 2026-11-15)

Labels: FACT (code or cited document), ESTIMATE (my calculation with the committed model; scripts `decomp.py` and `transport.py` are in `validation/plans/v3-008/p12-scripts/`, single process), OPEN. No reference table values are reproduced here.

## 1. What ICRU73_NEW is
- FACT [3]: Geant4 `G4WaterStopping.cc` says: "29.04.2009 A.Ivantchenko added revised data for G4WATER, provided by Prof. P. Sigmund, Univ. Southern Denmark, in the framework of the ESA Technology Research Programme (ESA contract 21435/08/NL/AT)". The table covers Z = 3–18 and 26 on a 0.025–1000 MeV/u grid. That is the same span as libdEdx `icru73_newEng.dat` (P5b found the grid identical to the G4 icru73 files).
- FACT [2,6,9]: the original ICRU 73 water table assumed I = 67.2 eV. Paul [6] pointed out the error using carbon-ion range data. The Errata and Addenda (P. Sigmund, J. ICRU 5(1)) recommend I = 78 eV and give revised water tables. An independent carbon stopping measurement agreed best with the Errata data at 78 eV [9].
- ESTIMATE: ICRU73_NEW is the 2009 Sigmund/Schinner revision, computed with I = 78 eV. Three things support this: (i) the provenance and grid above; (ii) ionmc at 78 eV agrees within 0.1 % above 50 MeV/u, whereas 75 eV would give a +0.4–0.6 % offset; (iii) the libdEdx "75 eV" is a material constant, not the table's I-value. Byte-level equality with the G4WaterStopping arrays is OPEN. It could be checked by a hash comparison that prints no values.
- FACT, general [14]: ICRU 73 values for Z ≥ 3 come from PASS (Sigmund–Schinner binary theory). OPEN: whether the 2009 water revision used PASS over the whole 2–20 MeV/u band. The primary text was not read.
- Was the revision mainly a low-energy change? No. Its documented driver is the I-value, so it acts mainly at high energy. Burigo & Greilich [9] report 1–2 % plateau differences when the unrevised table is used.

## 2. Is a 5–13 % deficit at 2–10 MeV/u expected?
- FACT [8]: Paul & Schinner 2006 compared stopping programs with experiment for heavy ions over 0.025–1000 MeV/u. MSTAR and SRIM have standard deviations of about 6–8 %. The original ICRU 73 had mean deviations and standard deviations of 8–23 %, partly from the water error. Data for C/O in liquid water at 2–20 MeV/u are sparse. OPEN: there is no published band-specific uncertainty for ICRU 73-revised water. For an expected uncertainty I use about 3–5 % at 2–10 MeV/u and about 1–2 % at 10–20 MeV/u (judgement).
- FACT [4]: Geant4 does not use raw Bethe plus effective charge just above 2 MeV/u. `G4EmModelManager::FillDEDXVector` multiplies the upper model by `(1.0 + del/e)` with `del = (dedx1/dedx2 - 1.0)*elow`. This forces continuity with the ICRU 73-based Bragg-type model at the 2 MeV-scaled switch and decays as 1/E. The fact that the main reference code needs this patch is indirect evidence that Bethe + effective charge falls short at a few MeV/u. Lechner et al. [13] added full ICRU 73 tables to Geant4 for this energy range.
- ESTIMATE: the observed deficit falls roughly as 1/E: C band maxima 10.1 % (2–10 MeV/u) and 2.5 % (≥ 10 MeV/u, attained in the 10–20 band), band means 1.8 % at 10–20 and 0.3 % at 20–50. It is negative, with the model low. Its size grows faster than z: the O/C ratio is 1.4 at 2–5, 2.2 at 10–20 and 2.7 at 20–50 MeV/u. A deficit of this size and sign at 2–10 MeV/u is not surprising for this model class. The 10–50 MeV/u part (1.8–3.9 % in the band mean, 5 % maximum for O in the 10–20 MeV/u band) is larger than I would expect from the charge state alone (see §3).

## 3. Decomposition of the ionmc model (ESTIMATE, committed code, water 78 eV)
Each correction term is given as a percentage of the full bracket L. "q-loss" is 1 − q²/z².

| E MeV/u | C q/z | C q-loss | C Barkas | C Bloch | C shell | O q/z | O q-loss | O Barkas | O Bloch | O shell |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | .916 | 16.1 | +4.2 | −9.1 | −1.3 | .871 | 24.2 | +5.5 | −13.1 | −1.3 |
| 3 | .952 | 9.4 | +2.7 | −6.3 | −1.2 | .918 | 15.7 | +3.5 | −9.5 | −1.2 |
| 5 | .980 | 4.0 | +1.4 | −3.8 | −1.0 | .960 | 7.8 | +1.8 | −6.1 | −1.0 |
| 10 | .996 | 0.8 | +0.5 | −1.8 | −0.8 | .989 | 2.1 | +0.7 | −3.1 | −0.8 |
| 20 | .9996 | 0.09 | +0.2 | −0.9 | −0.4 | .998 | 0.34 | +0.3 | −1.5 | −0.4 |
| 50 | 1.000 | 0.00 | +0.06 | −0.3 | −0.1 | .9999 | 0.01 | +0.08 | −0.6 | −0.1 |

The Mott term is +0.1–0.4 %. `ls="bare_mixed"` changes S by −1.2 % (C) and −2.6 % (O) at 2 MeV/u, which is the wrong direction, and by less than 0.1 % above 7 MeV/u. This is consistent with the observation that the LS options do not help.

What would have to change to close the gap (ESTIMATE):
- **Charge only.** S would need to rise by about 9–14 % at 2–5 MeV/u. That needs q higher by about 4.6 % (C) and 6.5 % (O), and by about 0.9 % / 2 % at 10–20 MeV/u. For C above about 10 MeV/u this would push q/z above 1. A charge-only fix is therefore ruled out above about 10 MeV/u.
- **Above about 10 MeV/u** the gap has to come from z-dependent higher-order terms or their interplay with the shell term. These terms are of the same order: Bloch −0.9 to −3 %, Barkas +0.2 to 0.7 %. Examples of changes that would close it: a Barkas term 3–5 times larger, or a Bloch term about 40–60 % smaller in magnitude, at 10–50 MeV/u. Both look physically implausible as single-term errors. PASS treats Barkas and Bloch together and non-perturbatively, including projectile screening. OPEN which term is the cause.
- **Discriminating test (no fit).** ICRU73_NEW also holds Z = 3–10 in water. At fixed velocity, compute table/model − 1 against z. Scaling with (1 − q²/z²) points to the charge state; scaling with z points to Barkas; scaling with z² points to Bloch. Running the same test on the existing proton and He rows (I1/I2 against ICRU 90) isolates the shell term, which is the same for every z.

Alternative charge parameterisations (ESTIMATE, q/z):

| E | C Barkas-125 | C Ziegler/G4 (+ZM) | C S-G 2001 mean | O Barkas-125 | O Ziegler/G4 | O S-G |
|---|---|---|---|---|---|---|
| 2 | .916 | .915 | .946 | .871 | .870 | .902 |
| 5 | .980 | .977 | .997 | .960 | .956 | .984 |
| 10 | .996 | .995 | 1.002* | .989 | .987 | 1.000 |

- Ziegler/G4: formula verbatim from `G4ionEffectiveCharge` [5], with the Ziegler–Manoyan screening correction and v_F = 0.94–1. It is 0.2–0.5 % below Barkas-125, so switching to it would lower S further by up to about 1 %.
- S-G 2001 [11]: coefficients recalled from secondary sources and NOT checked against the primary. Values above 1 (*) show it is being used outside its intended scope. It is a mean equilibrium charge, not a stopping effective charge. It would raise q² by 7 % (C) / 7 % (O) at 2 MeV/u and 3.5 % / 5 % at 5 MeV/u, about half of the deficit, and nothing above 10 MeV/u.
- Hubert–Bimbot–Gauvin 1990 [10]: OPEN. The formula was not retrieved because the PDF was too large. It is an empirical stopping-charge fit to about 600 data points for 16O–238U projectiles at 3–90 MeV/u, so carbon falls outside its projectile set.

## 4. Transport relevance (ESTIMATE; range from 2 MeV/u; Bohr straggling forward model)
- C 290: R = 163.2 mm. R(<20 MeV/u) = 1.45 mm (0.9 %); R(<50 MeV/u) = 7.5 mm (4.6 %).
- O 350: R = 166.6 mm. R(<20) = 1.11 mm (0.7 %); R(<50) = 5.7 mm (3.4 %).
- Range shift from a stopping increase below 20 MeV/u:
  - +5 %: −0.07 mm (C), −0.05 mm (O).
  - +10 %: −0.13 mm (C), −0.10 mm (O).
  - Using the observed band-mean profile: −0.06 mm (C 290, −0.035 %), −0.09 mm (O 350, −0.055 %), −0.06 mm (C 80, −0.33 %), −0.09 mm (O 100, −0.47 %).
  - This is consistent in sign and size with the reported maxima of 0.26 / 0.39 % over 80–430 MeV/u, which are probably reached at the low-energy end.
- Bragg peak. The energy deposited below 20 MeV/u is fixed; only its path is compressed, by about 0.1 mm. In my crude model (0.1 mm groups, phase noise about ±2 %):
  - C 290 / O 350 (σ_R about 0.4–0.5 mm): peak-to-entrance (PE) +1 to +3 %, peak shift ≤ one group.
  - C 80 / O 100 (σ_R about 0.06 mm): PE +3 to +10 %.
- Relation to TD-EM:
  - The I3 deficit cannot explain the low-end (d) FAILs. ionmc's PE is already about 11 % higher than TOPAS (C 80 11.47 vs 10.35; O 100 11.10 vs 9.96), and raising ionmc's stopping would increase PE further.
  - If TOPAS/Geant4 carries ICRU 73-like stopping at 2–20 MeV/u (through the `(1+del/e)` blend), the point-deposit excess is larger than observed, not smaller.
  - The R80 pattern points the same way. ionmc−TOPAS is −0.29 % for the upper beams but −0.09 % (C 80) and −0.02 % (O 100). The differences, +0.20 % / +0.27 %, match the band-profile prediction of +0.29 % / +0.41 %.
  - This undercuts the TD-EM defence that "stopping is shared" at low energy.
  - Discriminating evidence: Geant4 dE/dx for C/O in G4_WATER at 2–20 MeV/u, printed by the P4 smoke. Also confirm that `FillDEDXVector` behaves the same way in 11.4.2.

## 5. Follow-up options (for the orchestrator)
**(a) Different effective charge.**
- Ziegler/G4 goes the wrong way. S-G recovers at most about half of the 2–10 MeV/u gap and nothing at ≥ 10 MeV/u. Hubert is OPEN.
- It cannot meet the frozen ≥ 10 MeV/u gate.
- Evidence consumed: it is a model change after the results were seen. I3 stays FAIL as graded. Under STOP 3, the main TD/TL rows are consumed and TD-H/TL-H become the gating rows. It also changes He (q acts for z ≥ 2), so I2, the P3 LS bound and the range rows must be re-run.
- Fresh evidence needed: a held-out stopping reference that has not been seen. One candidate is ICRU 90 carbon in water. ICRU 90 does include carbon ions [9,12], contrary to the brief, but has no O and may share lineage with ICRU 73. Others are IAEA-database experimental points, plus TD-H/TL-H.

**(b) Tabulated low-energy stopping.** Use a table at 2 to about 20–50 MeV/u and blend into Bethe, Geant4-style or with a fixed match energy.
- Sources:
  - libdEdx ICRU73_NEW: GPL-3. As a runtime model input, the generated tables would carry GPL obligations. That needs a decision 0038-style check.
  - The same Sigmund data appear to be embedded in Geant4 source (`G4WaterStopping.cc`) under the Geant4 licence, which is permissive with attribution. That may be a cleaner licensing route, OPEN until a licence review.
  - G4EMLOW: restricted to Geant4 use.
  - SRIM: not redistributable.
  - MSTAR: licence OPEN.
- Problem: once ICRU73_NEW becomes model input, I3 is circular. Fresh independent evidence would be needed (ICRU 90 C, experimental data, TD-H/TL-H). The rows consumed are the same as in (a).
- This option is the most likely to reproduce TOPAS's low end.

**(c) Declared limitation.**
- Every beam passes through 2–20 MeV/u, so excluding an energy domain is not meaningful. What can be declared is an uncertainty: S up to −10 % (C) / −13.5 % (O) in the 2–5 MeV/u band, falling about as 1/E. Its consequences are R ≤ 0.13 mm (≤ 0.5 % at R ≈ 20 mm) and Bragg-peak PE of a few %.
- No model change and no rows consumed. I3 stays FAIL.
- It needs a plan-D amendment to the capability statement and a check against REQUIREMENTS for ion-specific stopping evidence.

Recommendation (non-binding): first run the cost-free z-scan attribution (§3) and the TOPAS dE/dx print. Then choose between (b) with a licence review and (c). (a) alone cannot meet I3.

## Sources
1. Lühr A. et al., Int. J. Radiat. Biol. 88 (2012) 209 (libdEdx). https://libdedx.github.io
2. ICRU, Errata and Addenda: ICRU Report 73, J. ICRU 5(1). https://www.icru.org/?p=97
3. Geant4 G4WaterStopping.cc. https://geant4.kek.jp/lxr/source/processes/electromagnetic/standard/src/G4WaterStopping.cc
4. Geant4 v11.2.0 G4EmModelManager.cc. https://raw.githubusercontent.com/Geant4/geant4/v11.2.0/source/processes/electromagnetic/utils/src/G4EmModelManager.cc
5. Geant4 v11.2.0 G4ionEffectiveCharge.cc (same path, utils/src).
6. Paul H., NIM B 255 (2007) 435. https://research.jku.at/en/publications/the-mean-ionization-potential-of-water-and-its-connection-to-the-/
7. IAEA stopping database, statistical analyses. https://www-nds.iaea.org/stopping/statisticalAnalyses (not retrievable, HTTP 402)
8. Paul H., Schinner A., NIM B 247 (2006) 166. https://ui.adsabs.harvard.edu/abs/2006NIMPB.247..166P (abstract via search summary)
9. Burigo L., Greilich S., arXiv:1812.07877 (ICRU 90 carbon; Errata at 78 eV). https://arxiv.org/abs/1812.07877
10. Hubert F., Bimbot R., Gauvin H., ADNDT 46 (1990) 1. doi:10.1016/0092-640X(90)90001-Z
11. Schiwietz G., Grande P.L., NIM B 175–177 (2001) 125. doi:10.1016/S0168-583X(00)00583-8
12. ICRU Report 90 (2016). https://www.icru.org/?p=107
13. Lechner A., Ivanchenko V., Knobloch J., NIM B 268 (2010) 2343. doi:10.1016/j.nimb.2010.04.008
14. Sigmund P., Schinner A., Eur. Phys. J. D 12 (2000) 425 (binary theory; recalled, not re-read).
