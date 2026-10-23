Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E6 swap test: summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de6_swap --de6-summary`

Input files (sha256):

- `validation/generated/elastic_attribution/de3/de3-150-all.json`  b42751ca39f13efe566bfe29ff85ca2f0940910ee98e7619f16b9f57bae11478
- `validation/generated/elastic_attribution/de3/de3-150-emonly.json`  5346b0bb447a33e38d7462f1c9c52974f0af87c077f2bf7264990fd940cd536a
- `validation/generated/elastic_attribution/de3/de3-150-swap-all.json`  dade944f58cae97039f23022df11d9fa6ef2ce644554a101211a694418b4dab3
- `validation/generated/elastic_attribution/de3/de3-150-swap-both.json`  a88f003e12b3a77c572cd747f528baff0004f05199419757d6fbb0910d86e93b
- `validation/generated/elastic_attribution/de3/de3-150-swap-pA.json`  4393b7534b54eab45d9bd5efc79c3f64f3bc172a7a1394760cc7a30780678685
- `validation/generated/elastic_attribution/de3/de3-150-swap-pp-sigma.json`  5ad8444a1e2c46e9ab8e406331bc829dc544370370a7c62e1c6cde12d2e28e4f
- `validation/generated/elastic_attribution/de3/de3-150-swap-pp.json`  b93d3dd796a8e05ea2d85d3fda5b76379e397fdadf775dd919a6005bf98160a3
- `validation/generated/elastic_attribution/de3/de3-200-all.json`  04ba39046251948c0edce8e3554dfa76fbcd61542d56a3985cafd4b5ed4a183f
- `validation/generated/elastic_attribution/de3/de3-200-emonly.json`  e2d093a1a724c303a47c5e7e55a7afca9375087b2bfd3df1c909fce996a47a91
- `validation/generated/elastic_attribution/de3/de3-200-swap-all.json`  27e513f5faf22995068cf877e22c6df859d47a3f22c574003d8ace21ce6ca8bc
- `validation/generated/elastic_attribution/de3/de3-200-swap-both.json`  112e017db2a8cc57c0a53e687015f885e0869f897cf29efac1619e958576edf4
- `validation/generated/elastic_attribution/de3/de3-200-swap-pA.json`  3a9ed03773b5c940ce01cc4a158bb4eafc41e388eea6301497d56d5cf1bb1a63
- `validation/generated/elastic_attribution/de3/de3-200-swap-pp-sigma.json`  f953bef0d4b795efe67d42e8de402b554ef29b45e79f3af55523dc6ecac28530
- `validation/generated/elastic_attribution/de3/de3-200-swap-pp.json`  06768ea7de370e016f3ef6d577dc7769b0e30b79fcbe1a4d2b67ea49c046ace7
- `validation/results/elastic_attribution/de0h5/h5.json`  40421563c3c55383d413de82af8775b36cb0a28f569a62435703762730ea76a7
- `validation/generated/elastic_attribution/de6/cdf-error.json`  f2af13a2891ebbb6065c91bab6daa843d9f5c68e5771e06897f0d574c3b54698
- `validation/generated/elastic_attribution/de6/de6-lineage.json`  e54dd5775a88dd911ca0a02eb9a921973e461d6deaafbe7a5cefb5672baf748d
- `validation/generated/elastic_attribution/de4/de4-ratios.json`  7807cb09ffb8dfba73a561979bd366f9e92c58a83363b254eaf71049b30b2365

<!-- body below is rendered from the Git projection of the JSON record -->

Clauses 1-3 of plan section 5 hold at both energies (rule of record, 1 - F): **True**; (-ln F sensitivity): True. Clauses 4 and 5 are not decided here.

Archived V11 quantities: committed H5 record `validation/results/elastic_attribution/de0h5/h5.json` (sha256 40421563c3c55383d413de82af8775b36cb0a28f569a62435703762730ea76a7), source RUN-20261012T142957Z-a8e41216 / lv5c-dd976c5-16/16-v11-compare.txt at dd976c5a (sha256 1f81c774d5209761fce23f757e6f6c2d4f731c533be9b31d0e017a253175f6e2); values digest f9bd66b4c3f572df4c3a4ae75e63dbeef3372485b9421c93d62874a254b0d1f7; archive cross-checked (values to 1e-12, source sha256 equal to the pin).

## 150 MeV

Clause 1: base F 0.8580 +- 0.0014 against archived 0.8568: difference +0.0013 (2 SE = 0.0042), table identity bit for bit: True, holds: True.

### 1 - F (rule of record; Delta = F_X - F_base, positive closes)

| quantity | value | SE |
|---|---|---|
| Delta(pA) | +0.0335 | 0.0009 |
| Delta(pp) | -0.0000 | 0.0001 |
| Delta(both) | +0.0334 | 0.0009 |
| interaction I | -0.0000 | 0.0000 |
| observed (archived V11) | +0.0348 | 0.0026 |
| closure residual r | +0.0014 | SE_r 0.0027 (swap 0.0009, CDF 0.0002) |
| threshold max(0.010, 2 SE_r) | 0.0100 | clause 2 holds: True |
| main effects pA / pp (fraction of Delta(both)) | 1.00 / -0.00 | named pA, |I|/|Delta(both)| 0.00, clause 3 holds: True |

Extended set (addendum (c)3; post-hoc arms declared after D-E4):

| quantity | value | SE |
|---|---|---|
| Delta(pA) | +0.0335 | 0.0009 |
| Delta(pp) | -0.0000 | 0.0001 |
| Delta(ppS) | +0.0022 | 0.0006 |
| Delta(all) | +0.0355 | 0.0011 |
| I_all = Delta(all) - sum of main effects | -0.0001 | 0.0005 |
| closure residual r_all | -0.0007 | SE_r 0.0028 (swap 0.0011, CDF 0.0002) |
| threshold max(0.010, 2 SE_r) | 0.0100 | clause 2 (extended) holds: True |
| main effects pA shape / pp shape / pp sigma (fraction of Delta(all)) | 0.94 / -0.00 / 0.06 | |I_all|/|Delta(all)| 0.00 |

### -ln F (sensitivity)

| quantity | value | SE |
|---|---|---|
| Delta(pA) | +0.0382 | 0.0010 |
| Delta(pp) | -0.0000 | 0.0001 |
| Delta(both) | +0.0382 | 0.0010 |
| interaction I | -0.0000 | 0.0000 |
| observed (archived V11) | +0.0398 | 0.0029 |
| closure residual r | +0.0016 | SE_r 0.0031 (swap 0.0010, CDF 0.0003) |
| threshold max(0.010, 2 SE_r) | 0.0100 | clause 2 holds: True |
| main effects pA / pp (fraction of Delta(both)) | 1.00 / -0.00 | named pA, |I|/|Delta(both)| 0.00, clause 3 holds: True |

Extended set (addendum (c)3; post-hoc arms declared after D-E4):

| quantity | value | SE |
|---|---|---|
| Delta(pA) | +0.0382 | 0.0010 |
| Delta(pp) | -0.0000 | 0.0001 |
| Delta(ppS) | +0.0025 | 0.0007 |
| Delta(all) | +0.0405 | 0.0012 |
| I_all = Delta(all) - sum of main effects | -0.0002 | 0.0005 |
| closure residual r_all | -0.0007 | SE_r 0.0032 (swap 0.0012, CDF 0.0003) |
| threshold max(0.010, 2 SE_r) | 0.0100 | clause 2 (extended) holds: True |
| main effects pA shape / pp shape / pp sigma (fraction of Delta(all)) | 0.94 / -0.00 / 0.06 | |I_all|/|Delta(all)| 0.01 |

Clause 3 of the frozen 2x2: **holds under both metrics**. Outcome label: **closure holds on the frozen 2x2 (clause 2)**.

### Gain and peak-versus-gain ordering (descriptive)

| X | Delta_g | SE | relative peak change | relative gain change | gain - peak (SE) | ordering |
|---|---|---|---|---|---|---|
| pA | +0.00694 | 0.00039 | +0.236 | +0.314 | +0.079 (0.013) | tail-like |
| pp | -0.00020 | 0.00002 | -0.000 | -0.009 | -0.009 (0.001) | core-like |
| both | +0.00674 | 0.00039 | +0.235 | +0.306 | +0.070 (0.013) | tail-like |
| ppS | +0.00121 | 0.00023 | +0.015 | +0.055 | +0.040 (0.007) | tail-like |
| all | +0.00813 | 0.00046 | +0.250 | +0.369 | +0.119 (0.012) | tail-like |

Observed gain reduction (ionmc - TOPAS, archived): +0.00585.

### H4: recoil channel (share of the observed gain excess; the section 2.2 exclusion bound is 20 %)

| arm | recoil MeV/primary | recoil plateau contribution | reduction vs base | share of gain excess | >= 20 % |
|---|---|---|---|---|---|
| all | 0.2678 | 0.00200 | +0.00000 | +0.00 | False |
| both | 0.1091 | 0.00054 | +0.00146 | +0.25 | True |
| pA | 0.1091 | 0.00054 | +0.00146 | +0.25 | True |
| pp | 0.2677 | 0.00200 | -0.00000 | -0.00 | False |
| ppS | 0.2686 | 0.00204 | -0.00003 | -0.01 | False |
| swap_all | 0.1090 | 0.00054 | +0.00146 | +0.25 | True |

## 200 MeV

Clause 1: base F 0.8455 +- 0.0027 against archived 0.8441: difference +0.0015 (2 SE = 0.0075), table identity bit for bit: True, holds: True.

### 1 - F (rule of record; Delta = F_X - F_base, positive closes)

| quantity | value | SE |
|---|---|---|
| Delta(pA) | +0.0380 | 0.0013 |
| Delta(pp) | +0.0002 | 0.0001 |
| Delta(both) | +0.0382 | 0.0013 |
| interaction I | +0.0001 | 0.0001 |
| observed (archived V11) | +0.0433 | 0.0034 |
| closure residual r | +0.0051 | SE_r 0.0037 (swap 0.0013, CDF 0.0002) |
| threshold max(0.010, 2 SE_r) | 0.0100 | clause 2 holds: True |
| main effects pA / pp (fraction of Delta(both)) | 0.99 / 0.00 | named pA, |I|/|Delta(both)| 0.00, clause 3 holds: True |

Extended set (addendum (c)3; post-hoc arms declared after D-E4):

| quantity | value | SE |
|---|---|---|
| Delta(pA) | +0.0380 | 0.0013 |
| Delta(pp) | +0.0002 | 0.0001 |
| Delta(ppS) | +0.0058 | 0.0008 |
| Delta(all) | +0.0445 | 0.0016 |
| I_all = Delta(all) - sum of main effects | +0.0005 | 0.0009 |
| closure residual r_all | -0.0012 | SE_r 0.0038 (swap 0.0016, CDF 0.0002) |
| threshold max(0.010, 2 SE_r) | 0.0100 | clause 2 (extended) holds: True |
| main effects pA shape / pp shape / pp sigma (fraction of Delta(all)) | 0.85 / 0.00 / 0.13 | |I_all|/|Delta(all)| 0.01 |

### -ln F (sensitivity)

| quantity | value | SE |
|---|---|---|
| Delta(pA) | +0.0440 | 0.0015 |
| Delta(pp) | +0.0002 | 0.0002 |
| Delta(both) | +0.0442 | 0.0015 |
| interaction I | +0.0001 | 0.0001 |
| observed (archived V11) | +0.0500 | 0.0038 |
| closure residual r | +0.0058 | SE_r 0.0041 (swap 0.0015, CDF 0.0003) |
| threshold max(0.010, 2 SE_r) | 0.0100 | clause 2 holds: True |
| main effects pA / pp (fraction of Delta(both)) | 0.99 / 0.00 | named pA, |I|/|Delta(both)| 0.00, clause 3 holds: True |

Extended set (addendum (c)3; post-hoc arms declared after D-E4):

| quantity | value | SE |
|---|---|---|
| Delta(pA) | +0.0440 | 0.0015 |
| Delta(pp) | +0.0002 | 0.0002 |
| Delta(ppS) | +0.0069 | 0.0010 |
| Delta(all) | +0.0513 | 0.0018 |
| I_all = Delta(all) - sum of main effects | +0.0003 | 0.0010 |
| closure residual r_all | -0.0013 | SE_r 0.0043 (swap 0.0018, CDF 0.0003) |
| threshold max(0.010, 2 SE_r) | 0.0100 | clause 2 (extended) holds: True |
| main effects pA shape / pp shape / pp sigma (fraction of Delta(all)) | 0.86 / 0.00 / 0.13 | |I_all|/|Delta(all)| 0.01 |

Clause 3 of the frozen 2x2: **holds under both metrics**. Outcome label: **closure holds on the frozen 2x2 (clause 2)**.

### Gain and peak-versus-gain ordering (descriptive)

| X | Delta_g | SE | relative peak change | relative gain change | gain - peak (SE) | ordering |
|---|---|---|---|---|---|---|
| pA | +0.00393 | 0.00042 | +0.246 | +0.158 | -0.088 (0.018) | core-like |
| pp | +0.00023 | 0.00002 | +0.001 | +0.009 | +0.008 (0.001) | tail-like |
| both | +0.00417 | 0.00042 | +0.247 | +0.167 | -0.080 (0.018) | core-like |
| ppS | +0.00287 | 0.00027 | +0.038 | +0.115 | +0.077 (0.009) | tail-like |
| all | +0.00694 | 0.00052 | +0.288 | +0.279 | -0.009 (0.019) | indeterminate |

Observed gain reduction (ionmc - TOPAS, archived): +0.00830.

### H4: recoil channel (share of the observed gain excess; the section 2.2 exclusion bound is 20 %)

| arm | recoil MeV/primary | recoil plateau contribution | reduction vs base | share of gain excess | >= 20 % |
|---|---|---|---|---|---|
| all | 0.3445 | 0.00155 | +0.00000 | +0.00 | False |
| both | 0.1274 | 0.00039 | +0.00116 | +0.14 | False |
| pA | 0.1274 | 0.00039 | +0.00116 | +0.14 | False |
| pp | 0.3444 | 0.00155 | +0.00000 | +0.00 | False |
| ppS | 0.3448 | 0.00155 | +0.00000 | +0.00 | False |
| swap_all | 0.1272 | 0.00039 | +0.00116 | +0.14 | False |

## Section 2.2 signature inputs (from the D-E4 ratios record; inputs, not decisions)

| quantity (30 / 50 / 100 / 150 / 200 MeV) | values |
|---|---|
| H1 rate_ratio_O | 1.005 / 1.006 / 1.006 / 1.013 / 1.013 |
| H1 rate_ratio_C | 1.025 / 1.014 / 1.031 / 1.019 / 1.008 |
| H2 exceedance_ratio_2mm | 1.568 / 1.670 / 1.541 / 1.343 / 1.287 |
| H2 exceedance_ratio_3mm | 1.341 / 2.277 / 1.904 / 1.557 / 1.453 |
| H2 theta_med_ratio | 1.069 / 1.196 / 1.378 / 1.338 / 1.348 |
| H2 T_ratio | 2.407 / 7.207 / 14.196 / 16.409 / 24.965 |
| H3 pp_rate_ratio | 0.977 / 0.997 / 1.032 / 1.165 / 1.159 |
| H3 pp_T_slow_ratio | 1.047 / 1.021 / 1.013 / 0.991 / 1.022 |
| H3 pp_exceedance_ratio_2mm | 1.068 / 1.032 / 1.008 / 0.996 / 1.000 |

H1: O rate ratio outside 1 +- 0.10 at 0 energies; all O and C within 1 +- 0.05: True. H2: exceedance ratio >= 1.15 at 5 energies, T ratio >= 1.5 at 5, all within 1 +- 0.10: False. H3: any p-p rate or T_slow ratio outside 1 +- 0.10: True; all within 1 +- 0.05: False.

## Finite-sample CDF error on the predicted F (addendum (c)4): bootstrap sd, tail-edge systematic, leave-one-energy-out (150 MeV value, proxy at 200), quadrature

| E | (a) bootstrap sd | (b) tail edge | (c) leave-one-out | total |
|---|---|---|---|---|
| 150 | 0.00019 | +0.00009 | +0.00011 | 0.00024 |
| 200 | 0.00020 | +0.00003 | +0.00011 | 0.00023 |
