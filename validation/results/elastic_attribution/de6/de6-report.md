Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E6 swap test: summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de6_swap --de6-summary`

Input files (sha256):

- `validation/generated/elastic_attribution/de3/de3-150-all.json`  26fa08a6f50af52a72fff8ac13428bb07434e77db3399f379cd010a35f8cd2bf
- `validation/generated/elastic_attribution/de3/de3-150-emonly.json`  ccbc5da52e57416246124405d1b8206c89de2d1667c0e1de900f53afa9ad4d56
- `validation/generated/elastic_attribution/de3/de3-150-swap-all.json`  994929f33b5cd2b2803427a66772d06edbec47dccbcdb8d7b79e89c730767d1a
- `validation/generated/elastic_attribution/de3/de3-150-swap-both.json`  1d422f15814a3333b2894ab01e76fcea0412aff546a109f19ce0f828ed69db01
- `validation/generated/elastic_attribution/de3/de3-150-swap-pA.json`  60b7c2496d040cf1ff8ab17458c7bc7572a5ad75f7454d53697aaa1892c5798c
- `validation/generated/elastic_attribution/de3/de3-150-swap-pp-sigma.json`  ec93d6e4dea2eb375908c74df158b1d75f39f9c022d1a439291a800c311ef188
- `validation/generated/elastic_attribution/de3/de3-150-swap-pp.json`  046dab13892e0edcabd4e6116b83d1b42bafaed15dbbddff0e3112ba64f4eb67
- `validation/generated/elastic_attribution/de3/de3-200-all.json`  81956f90ccdb89380bc26240be370d91afa11df9a889aa658e734e71fbdd1538
- `validation/generated/elastic_attribution/de3/de3-200-emonly.json`  09361d6649f64994ea88b98d219c211a75b6ad63711a62c96faaf83f60424ad3
- `validation/generated/elastic_attribution/de3/de3-200-swap-all.json`  9c3465491f19db9c2caba250c91ac4ed63ea25101d7368da79c870533aa558f9
- `validation/generated/elastic_attribution/de3/de3-200-swap-both.json`  f3f5eac1ce65e3be09c9f31bb0996328a51b40a320e9548dd162a8bb720da568
- `validation/generated/elastic_attribution/de3/de3-200-swap-pA.json`  ca3dc9c7c1459de37a0af26ac3f17355076fb9ed9057535a68581f516cd43c59
- `validation/generated/elastic_attribution/de3/de3-200-swap-pp-sigma.json`  bd616145e1d3d2a10fc101144e94d98d546c4535714468389dce6740bc4fd3ce
- `validation/generated/elastic_attribution/de3/de3-200-swap-pp.json`  cc1b3a375c669f1dea23bea4afe78e115a4bd849ab377f80a534fd8814b53cbc
- `validation/generated/elastic_attribution/de6/cdf-error.json`  588898b2cc961195eae7b2e838c21c804d91a8481fcee5314846400086596a86
- `validation/generated/elastic_attribution/de4/de4-ratios.json`  7807cb09ffb8dfba73a561979bd366f9e92c58a83363b254eaf71049b30b2365

<!-- body below is rendered from the Git projection of the JSON record -->

Clauses 1-3 of plan section 5 hold at both energies (rule of record, 1 - F): **True**; (-ln F sensitivity): True. Clauses 4 and 5 are not decided here.

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
| 150 | 0.00018 | +0.00009 | +0.00011 | 0.00023 |
| 200 | 0.00018 | +0.00003 | +0.00011 | 0.00021 |
