Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E9 record: B0 / B1 / B2 / B3 / A1 and the Lead A sensitivity table

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de9_record`

Input files (sha256):

- `validation/results/transport/lv5c-dd976c5-single-process.json`  c35a62d45b8ea98b67031e072aa424a9355b5db3a94b04ccd182047e8e1c708b

<!-- body below is rendered from the Git projection of the JSON record -->

D-E9: exploratory, report-only; leads identified after the D-E8 STOP; not a re-evaluation. The D-E8 STOP stands; no band, seed or table of record changes.

## B0: independent p-p level reconstruction: match True, dip persists True; reading: evaluation feature at the LA150 top nodes (the independent reconstruction reproduces the step).

## B1: decomposition of the elastic plateau gain

| E | gain all | pp-only | pA-only | pp share of sum | additivity gap / all |
|---|---|---|---|---|---|
| 150 | 0.02207 | 0.01368 | 0.00809 | 0.628 | 0.013 |
| 200 | 0.02491 | 0.01982 | 0.00475 | 0.807 | 0.013 |

## B3: report-only variants (gain change against all, CRN; never candidates)

| variant | E | delta gain | SE |
|---|---|---|---|
| pp-level-0.90-ge150 | 150 | -0.00002 | 0.00001 |
| pp-level-0.90-ge150 | 200 | -0.00193 | 0.00019 |
| pp-level-0.95-ge150 | 150 | -0.00001 | 0.00001 |
| pp-level-0.95-ge150 | 200 | -0.00081 | 0.00009 |
| pp-level-monotone | 150 | -0.00028 | 0.00007 |
| pp-level-monotone | 200 | -0.00327 | 0.00027 |

Slopes dg / d ln s of the scaled family: 150 MeV 0.00023 (0.90-0.95), 0.00012 (0.95-1.00); 200 MeV 0.02074 (0.90-0.95), 0.01583 (0.95-1.00).

## A1: realisations of the pinned-table gain

| E | K | mean | SD | SE mean | rho | SE calibration | z archived | z D-E3 |
|---|---|---|---|---|---|---|---|---|
| 150 | 4 | 0.02222 | 0.00041 | 0.00021 | 0.513 | per-run SE consistent | -0.16 (consistent with the realisation distribution) | -0.19 (consistent with the realisation distribution) |
| 200 | 8 | 0.02573 | 0.00116 | 0.00041 | 1.156 | per-run SE consistent | 1.13 (consistent with the realisation distribution) | -0.76 (consistent with the realisation distribution) |

## Lead A sensitivity table at 200 MeV (report-only; never compared with the band)

| reading | d | SE | label |
|---|---|---|---|
| A1 realisation mean (K fixed) | 0.00363 | 0.00057 | sensitivity; frozen d = 0.00462 stands |
| archived baseline - 2 SE | 0.00307 | – | sensitivity; frozen d = 0.00462 stands |
| archived baseline + 2 SE | 0.00616 | – | sensitivity; frozen d = 0.00462 stands |
| D-E3 baseline | 0.00281 | – | sensitivity; frozen d = 0.00462 stands |
| the prediction as frozen | 0.00462 | 0.00087 | the frozen value; the STOP stands |
| inverse-variance pooled baseline | 0.00393 | 0.00061 | sensitivity; frozen d = 0.00462 stands |

## B2: share of the elastic plateau gain

- 150 MeV: secondary-proton share 0.407, recoil share 0.094, remainder 0.499
- 200 MeV: secondary-proton share 0.651, recoil share 0.054, remainder 0.295

Leak check (no absolute cross section): passed True.
