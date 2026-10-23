# PP-LEVEL — the pinned table's p-p level against SAID and BGG (physics researcher, 2026-10-19; report-only; ratios and verdict wording only; obtained after the D-E8 STOP; no table of record changes)

Projected record committed with validation/plans/v3-005d-investigation-de9.md (D-E9 Lead B). The helper scripts were not kept; the ratios are recomputable from the pinned table in the cache, the SAID SM16 SGT/DSG acquisitions (pins in src/ionmc/data/elastic_m2_inputs.json) and ionmc.nuclear.bgg / elastic_glauber.

# PP-LEVEL - p-p level of the pinned elastic table (research agent, 2026-10-19; report-only, scratch; ratios only)
Table: pin 1fd24cff (H-1 row = Hale LA150 LAW=5 LTP=1 NI = sigma_e - sigma_c over |mu_CM| <= 0.96, half sphere, ratio
interpolation between LAW=5 nodes, <= 150 MeV; above 150 MeV sigma_NI(150) * S(E)/S(150), S = Geant4/BGG p-p formula).
References: SAID SM16 pp SGT (m2-official/said); SAID SM16 pp DSG scans (m2-official/dsg) minus the ENDF-102 Mott density
over the same |mu_CM| <= 0.96 cut ("anchor", like-for-like with the table); BGG p-p = ionmc.nuclear.bgg.pp_bgg_mb.
Scripts: $TMPDIR pp.py / pp2.py (not kept; recomputable from the cache and the scratch SAID files).

## (1) Ratios
 E MeV | table/SAID SGT | table/BGG | table/SAID(DSG-Mott, cut) | SAID(DSG-Mott,cut)/SGT
  30   | 0.854 | 0.903 |   -                  | -
  50   | 0.844 | 0.947 |   -                  | -
 100   | 0.863 | 0.989 | 0.979                | 0.881
 125   | 0.828 | 0.920 | 0.922                | 0.898
 130   | 0.821 | 0.905 | ~0.91 (interpolated) | -
 140   | 0.859 | 0.936 | ~0.95 (interpolated) | -
 150   | 1.021 | 1.101 | 1.119                | 0.913
 175   | 1.043 | 1.101 | 1.128                | 0.925
 200   | 1.057 | 1.101 | 1.131                | 0.934
 250   | 1.068 | 1.101 | 1.129                | 0.946
 Non-monotonic feature: table/SGT falls to its minimum 0.821 at 130 MeV, then rises 0.824 (130.7) -> 0.859 (140) -> 1.021 (150);
 d sigma_table/dE changes sign from negative to positive at ~130 MeV and back to negative above 150 MeV (S(E) scaling);
 SAID SGT and SAID(DSG-Mott, cut) are monotonically falling, their ratio smooth (0.881 -> 0.946, 100 -> 250 MeV).
 Above 150 MeV table/BGG = 1.1008 at every node by construction (the 150 MeV node level is propagated by S(E)/S(150)).
## (2) Geant4 reference of the D-E4 p-p ratio
 D-E4 applies the same 16.26 deg CM cut to Geant4; Geant4 accepted fraction f = kept/(on-curve) = 0.961/0.962/0.964/0.955/0.959.
 Prediction "table / (BGG * f)": 0.940 / 0.985 / 1.026 / 1.153 / 1.148 at 30/50/100/150/200 MeV
 vs D-E4 0.977(6) / 0.997(9) / 1.032(12) / 1.165(15) / 1.159(15): agreement within 1 sigma at 100-200 MeV, 1.3 sigma at 50 MeV;
 30 MeV off by ~6 sigma (unexplained here; slab energy loss / acceptance - not investigated).
 => the Geant4 p-p rate behaves as BGG sigma_pp (G4BGGNucleonElasticXS, Z=1: 1.0115 x NS low-energy formula; inferred from the
 agreement, Geant4 source not re-read). The D-E4 1.165/1.159 is the table's +10 % step over BGG that enters at 150 MeV.
 BGG itself is 0.93-0.96 of SAID SGT at 150-250 MeV (full sphere, nuclear only).
## (3) Origin
 "Hale-NI partial": LA150/ENDF-B-VIII.0 H-1 MF6 MT2 LAW=5 LTP=1 (Hale R-matrix/phase-shift evaluation, nuclear amplitude with
 Coulomb-nuclear interference), NI = sigma_e - sigma_c (Mott subtracted, interference kept), integrated over |mu_CM| <= 0.96 on the
 half sphere. A free-pp elastic cross section and its cut nuclear+interference integral fall smoothly with T in this range
 (SAID: SGT and DSG-Mott both monotone); a 20 % rise between 130 and 150 MeV has no physical counterpart (no threshold below
 pi production at 280 MeV).
## (4) Verdict
 The p-p level excess at >= 150 MeV is a DEFECT of the pinned construction, not a Geant4 convention difference: against the
 like-for-like SAID reference (same cut, same Mott subtraction) the table goes from -8 % (125 MeV) to +12 % (150 MeV) and stays
 +13 % to 250 MeV through the S(E) scaling of the 150 MeV node. Undetermined: whether the anomaly is in the LA150 H-1 LAW=5 data at
 its top nodes or in the reconstruction (discriminator: independent LAW=5 NI reconstruction at the 130/140/150 MeV nodes).
 The interference/cut convention cannot explain it (it is present at 100-125 MeV, where the table is -2 to -8 % vs SAID).
