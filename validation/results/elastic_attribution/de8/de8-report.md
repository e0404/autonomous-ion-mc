Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E8: go/no-go screen of the model under test M2' (Amendment 20 (c), Amendment 21 (e))

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de8_diagnostic --reviewed-head-review-id REVIEW-5949c525 --reviewed-sha c509e351a4c402cf99f2993e2416d3560347d5ff`

Input files (sha256):

- `validation/results/transport/lv5c-dd976c5-single-process.json`  c35a62d45b8ea98b67031e072aa424a9355b5db3a94b04ccd182047e8e1c708b
- `validation/results/elastic_attribution/de7/de7.json`  302d03a4c77bbdaebb5db17b20fa632f4e306630305beb17745391f05b1a4d73
- `validation/results/elastic_attribution/de6/de6-report.md`  44496c43b6f245efdfc33a96a2c672ca168c5c9587cf215f6a18b6351fe99c8b
- `validation/results/elastic_attribution/de4/de4-ratios.json`  696c6dcd7fb780f2c020dbc74d9fc31994ac0f95e302413e9fb87ea0dfa0d47c
- `validation/results/elastic_attribution/de5/de5.json`  f170c1d382bbbe8a9226759b558638ac12428312a9b40fbd3a461caaab1434d2

<!-- body below is rendered from the Git projection of the JSON record -->

Go/no-go screen, Amendment 20 (c) as carried by Amendment 21 (e) and its addendum (a); report-only, never evidence. It uses the V11 reference (archived TOPAS F and gain), the V5 reference (TOPAS ΔIDD and peak/plateau) and the Geant4 thin-target ratios to decide whether the single pre-committed model M2′ is transported at all. It is not a selection among models: M2′ has no free parameter, plan §4.2 forbids switching to M3 on such outcomes, and STOP leads to "no model change". It is a stopping rule on gating quantities. The lv5d V11 and V5 runs are therefore not blind to these estimates; they replicate them on independent seeds with the frozen comparators. The M2′ construction outcomes were known before registration (Amendment 21 (a)3) and the de7 record before this addendum; the screen's transport quantities were not. Condition (iv) was revised after de7 to 100 / 150 / 200 MeV and made report-only: its outcome there is predictable from the disclosed construction-level P values. The decisive blind conditions are (i), (i-b), (ii) and (iii). Executed once; recorded whichever way it falls.

## Verdict: **STOP**

Failing conditions: (i) at 200 MeV, (i-b) at 200 MeV, (ii) at 200 MeV

condition (iv) is recorded, never part of the verdict; the 30 / 50 MeV cases are identical to the pinned table by construction for M2' and R-B.

## Bands (recomputed from the frozen rule)

| E | F closure window | allowed |q| (i-b) | allowed |d_g| (ii) | allowed |p| (iii) |
|---|---|---|---|---|
| 150 | 0.631 to 1.369 | 0.01279 | 0.00360 | 0.0460 |
| 200 | 0.736 to 1.264 | 0.00961 | 0.00365 | 0.0676 |

## Arm m2  (report-only)

Table id `ad0d77578f29cde5dd9c4747d847cfcc586ec8bfb0643544aaa8addcfeb82034`; documents valid: True.

| E | D1 closure (SE) | D1b q (SE) | D2 d_g (SE) | D3 p (SE) | (i) | (i-b) | (ii) | (iii) |
|---|---|---|---|---|---|---|---|---|
| 150 | 1.046 (0.033) | 0.0042 (0.0013) | -0.00076 (0.00036) | 0.0113 (0.0051) | pass | pass | pass | pass |
| 200 | 0.843 (0.034) | -0.0079 (0.0017) | 0.00463 (0.00040) | 0.0345 (0.0036) | pass | pass | fail | pass |

## Arm m2prime  (model under test)

Table id `addd817961ca24ebde185698ece4b38440a16cfdcf3b21dfad9416a089519b47`; documents valid: True.

| E | D1 closure (SE) | D1b q (SE) | D2 d_g (SE) | D3 p (SE) | (i) | (i-b) | (ii) | (iii) |
|---|---|---|---|---|---|---|---|---|
| 150 | 0.635 (0.023) | -0.0118 (0.0009) | -0.00088 (0.00036) | 0.0097 (0.0050) | pass | pass | pass | pass |
| 200 | 0.604 (0.030) | -0.0195 (0.0015) | 0.00462 (0.00039) | 0.0344 (0.0035) | fail | fail | fail | pass |

## Arm rb  (report-only)

Table id `2a14cdcb93db416f30db0b0fb3bf1607af414590fadd92ea0fe3c06f802ec221`; documents valid: True.

| E | D1 closure (SE) | D1b q (SE) | D2 d_g (SE) | D3 p (SE) | (i) | (i-b) | (ii) | (iii) |
|---|---|---|---|---|---|---|---|---|
| 150 | 0.611 (0.022) | -0.0128 (0.0009) | -0.00086 (0.00036) | 0.0100 (0.0050) | fail | pass | pass | pass |
| 200 | 0.586 (0.031) | -0.0204 (0.0015) | 0.00463 (0.00039) | 0.0345 (0.0035) | fail | fail | fail | pass |

## D4 (report-only for (iv)): P(theta_lab > theta_2mm), arm / Geant4

Table `ad0d77578f29cde5dd9c4747d847cfcc586ec8bfb0643544aaa8addcfeb82034`: moves toward 1 at 5 of 5; in-domain 3 of 3: True.

| E | arm/table | table/G4 | arm/G4 | toward 1 | identical by construction |
|---|---|---|---|---|---|
| 100 | 0.6595 | 1.5405 | 1.0160 | True | False |
| 150 | 0.7841 | 1.3427 | 1.0528 | True | False |
| 200 | 0.8325 | 1.2867 | 1.0711 | True | False |
| 30 | 0.3382 | 1.5676 | 0.5302 | True | False |
| 50 | 0.4462 | 1.6703 | 0.7453 | True | False |

Table `addd817961ca24ebde185698ece4b38440a16cfdcf3b21dfad9416a089519b47`: moves toward 1 at 3 of 5; in-domain 3 of 3: True.

| E | arm/table | table/G4 | arm/G4 | toward 1 | identical by construction |
|---|---|---|---|---|---|
| 100 | 0.6270 | 1.5405 | 0.9659 | True | False |
| 150 | 0.7665 | 1.3427 | 1.0291 | True | False |
| 200 | 0.8222 | 1.2867 | 1.0579 | True | False |
| 30 | 1.0002 | 1.5676 | 1.5679 | False | True |
| 50 | 1.0028 | 1.6703 | 1.6751 | False | True |

## (v) construction record: holds (de7 sha256 `302d03a4c77bbdaebb5db17b20fa632f4e306630305beb17745391f05b1a4d73`)

Provenance: HEAD `c509e351a4c402cf99f2993e2416d3560347d5ff`, clean tree True.

Leak check (no absolute cross section, no measured value): passed True.
