Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E0: deterministic kernels of the elastic channel

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de0_kernels`

Input files (sha256):

- `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.npz`  1028bb67f57bc96edafe6d3dbcce532013edbbd31b9b80861239b09d9ad4202a
- `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.json`  ab5e2a357ab53060a84182607d93fe4373e50ce60622450af74728fd4e495e77
- `validation/generated/RUN-20261012T142957Z-a8e41216/transport/lv5c-dd976c5-16/16-v11-compare.txt`  1f81c774d5209761fce23f757e6f6c2d4f731c533be9b31d0e017a253175f6e2
- `validation/generated/RUN-20261012T142343Z-081123c3/transport/lv5c-dd976c5-03/03-v10-a.txt`  05499b3a30a8695544f153dea9cb7367c744c2ea43e3e4e225a2417eb7397cf9
- `.ionmc-cache/transport-samples/dd976c5/12/samples/v11-150-emel.json`  3f3aca4c5c32c26175372bf1a7fcbfc4d0b3eb095e830e664c4928da8632fd6a
- `.ionmc-cache/transport-samples/dd976c5/13/samples/v11-150-emonly.json`  cb46ec341d0334740ea0e585f5ea9c0f084d60a883c4097c80a5b070837c0b6f
- `.ionmc-cache/transport-samples/dd976c5/14/samples/v11-200-emel.json`  e9002536a2580e3f12f6c66bb5644cf1de050718f1682f7be1695dab890f89b3
- `.ionmc-cache/transport-samples/dd976c5/15/samples/v11-200-emonly.json`  1bb5e78e387cdd8528895d6c22b196fe1ba3f792db34d72a4bbcbbd81cc7866b
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emelastic-seed1/case.json`  e7d6d7174b6d8ba04170db366dd8e3ea0ecb5228439861d61427f357a002f005
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emelastic-seed1/input.txt`  1195d133a0743260840d59fd0b65627e3a0a88333b18b344383e016229cc847d
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emelastic-seed2/case.json`  af11728c5a575cd115e3a015a30f2d6df753d9837a8b64aa758d135e993c0932
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emelastic-seed2/input.txt`  c8b63cc6966f0725b4827aadc67fece1125e232a128689eeba3773a0c81797fa
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emelastic-seed3/case.json`  63eae4445d66daef9d03e3f44fb96c171e8dfbcda4d98f97b2b35359c8c0b370
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emelastic-seed3/input.txt`  ea0993fccaae8ba88584a17a18a7489c2f348a7b01b8f0f2969c02d939e240ce
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emonly-seed1/case.json`  4710de3ebf9d395e47c920bc76837818e75311f1850a6526b0dd77a8f067ccaa
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emonly-seed1/input.txt`  f82a6328e70f7049650201234c228265bafdc3d458f454ec68f13ce482a924de
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emonly-seed2/case.json`  54d0eb6511e8b4d302b22ef299e9eac001a5abcd68b7730dfca35e1bf1358c4e
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emonly-seed2/input.txt`  24bc2bec7722ffcc54c3c6a42c6ec7616ce6a0ae4474d9ae628e94e68bd28ed5
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emonly-seed3/case.json`  32ac048798eb1a0f3a0bba1748901b922ea2aa4ba97c370953809bb47abcae92
- `validation/reference_cases/topas/proton-water-150mev-idd-r20-emonly-seed3/input.txt`  cacf94431e94aaf1be0ce0961d5fbcd4b31e2bd86610cb04a4719b78eaf3a030
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emelastic-seed1/case.json`  7e8eebade4e18a0e50edcbc9cfa2e114fb80110e1caed45281273f8caf4672ea
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emelastic-seed1/input.txt`  841265268397cb83af9025f84ad1f95f9206c45eff8fb92f0b151f023c1d2c62
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emelastic-seed2/case.json`  44ebacad2b0c0673fc83574becf30cd044e80205ea0aefbf968311dc7fdd79a6
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emelastic-seed2/input.txt`  dff6141ed90c8f92dc556bf7f080c3be35c705186c54a556b534e9ad406cf168
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emelastic-seed3/case.json`  aa761d568a41e37716917b91bc920a0fc969fa0dbd39818fbc8bb3e0634f6c58
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emelastic-seed3/input.txt`  258714219c98105f033289dfdfcfef4fd5aff907729a989a87be734200c7140a
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emonly-seed1/case.json`  8fe4d5dd6bd11337e6d8b0fcf3a86974835d21d6b4503c2a6d26b8636acd5cdb
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emonly-seed1/input.txt`  bf07f7152859dcfaff455577818290ee910393d5c5ff21c9eb09a6efcf80d6e5
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emonly-seed2/case.json`  cfbc24016f6d7d4edba6da36cd03a08a83c6edad161efc6a4e3afe03b6f94bed
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emonly-seed2/input.txt`  d8830a0cdb85ef17f8846b7bbc66a283e9ee141fd7ca6f912963c4757881284d
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emonly-seed3/case.json`  a4285a57b427cc72fc2c4f4f60c6ca9e45492645e2fa9956338a5aa3fa0c087a
- `validation/reference_cases/topas/proton-water-200mev-idd-r20-emonly-seed3/input.txt`  b25768a8ea5b15db3e5fef3badb945f96a61e73b31fde11890d665fe1e691b22

Definitions: module docstring of `de0_kernels.py` (plan section 2.1). Water, Bethe CSDA range (the source of the V5/V11 runs), straight-line path, exact two-body kinematics, no straggling, 2048 equiprobable quantile points per event, 300 depth nodes. The plan's G wording is under-specified; the operational form is in the docstring.


## Dated note (2026-10-14): sign of the frozen K* formula

The plan (section 2.1, frozen by Amendment 19 (a)) writes K* = int n sigma <1 - IDD_EM(z_pk - Delta)/IDD_EM(z_pk)> dz. This formula has a sign error: a proton whose range is shortened by Delta has the depth-dose curve IDD_EM(z + Delta), so its contribution at z_pk is IDD_EM(z_pk + Delta) (distal fall-off); IDD_EM(z_pk - Delta) is the proximal rise, i.e. a range *extended* by Delta. The plan's own prose ('moves that proton's contribution ... to z_pk - Delta') describes the shifted form. The as-frozen form is computed and labelled 'as frozen (sign error)'; **no conclusion is drawn from it**. The shifted form is the physically meaningful kernel (F predicted 0.853 / 0.832 against 0.857 / 0.844 observed in V11; peak removal about 10 % high). In D-E3 both forms are to be reported and validated against the shifted one. The frozen plan is not edited. The kernel is adequate for ranking and sensitivity to about 10 %, not for a 0.010 closure.

G caveat: G (operational form, see the docstring) omits the higher stopping power of the proton that continues after a p-p event (E_fast = E - T_slow); that term is positive. The column `G_total incl. post-p-p stopping` adds the post-event stopping power of every continuing proton at z' (all events, exact residual range); its size is in `post_event_stopping_term` of kernels.json. G is still the plan's incomplete descriptive quantity and enters only the descriptive peak-versus-gain ordering.

## Kernels in water (p+O and p-p; per primary)

| E | shape | P_el(O) | P_el(pp) | K_1mm | K_2mm | K_3mm | K_5mm | K* as frozen (sign error; no conclusion) | K* shifted (meaningful) | K*_gauss (shifted form) | G_total | G_path | G_recoil | G_pp_slow(up.) | G_total incl. post-p-p stopping |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 50 | table | 0.07256 | 0.01204 | 0.04641 | 0.02629 | 0.01701 | 0.01011 | - | - | - | - | - | - | - | - |
| 50 | disk | 0.07256 | 0.01204 | 0.04641 | 0.02629 | 0.017 | 0.0101 | - | - | - | - | - | - | - | - |
| 100 | table | 0.1689 | 0.02494 | 0.1345 | 0.09715 | 0.07415 | 0.04916 | - | - | - | 0.007725 | 0.0009738 | 0.002343 | 0.004408 | 0.02046 |
| 100 | disk | 0.1689 | 0.02494 | 0.1345 | 0.09714 | 0.07413 | 0.04914 | - | - | - | 0.007753 | 0.001003 | 0.002341 | 0.004408 | 0.02052 |
| 150 | table | 0.2271 | 0.03618 | 0.1987 | 0.1553 | 0.1255 | 0.08949 | 0.08366 | 0.1356 | 0.1511 | 0.01336 | 0.003081 | 0.002015 | 0.008261 | 0.01936 |
| 150 | disk | 0.2271 | 0.03618 | 0.1987 | 0.1553 | 0.1255 | 0.08947 | 0.08365 | 0.1356 | 0.1511 | 0.01339 | 0.00311 | 0.002014 | 0.008261 | 0.0194 |
| 150 | Rx1.2 | 0.2271 | 0.03618 | 0.1794 | 0.132 | 0.1035 | 0.07275 | 0.07311 | 0.1148 | 0.1304 | 0.01207 | 0.002174 | 0.001633 | 0.008261 | 0.01705 |
| 150 | s0.5 | 0.2271 | 0.03618 | 0.1904 | 0.1431 | 0.1114 | 0.074 | 0.07417 | 0.1225 | 0.1391 | 0.01103 | 0.001909 | 0.0008575 | 0.008261 | 0.01555 |
| 150 | s0.9 | 0.2271 | 0.03618 | 0.178 | 0.127 | 0.09483 | 0.05936 | 0.06552 | 0.107 | 0.1242 | 0.01022 | 0.001403 | 0.0005577 | 0.008261 | 0.01451 |
| 150 | gauss | 0.2271 | 0.03618 | 0.2032 | 0.1613 | 0.1308 | 0.0911 | 0.08277 | 0.1401 | 0.1558 | 0.01107 | 0.001876 | 0.0009341 | 0.008261 | 0.01546 |
| 200 | table | 0.2611 | 0.0519 | 0.2485 | 0.2028 | 0.1701 | 0.1279 | 0.08962 | 0.1506 | 0.1608 | 0.02051 | 0.002064 | 0.001522 | 0.01692 | 0.02404 |
| 200 | disk | 0.2611 | 0.0519 | 0.2485 | 0.2028 | 0.1701 | 0.1279 | 0.0896 | 0.1506 | 0.1608 | 0.02047 | 0.002023 | 0.001524 | 0.01692 | 0.02397 |
| 200 | Rx1.2 | 0.2611 | 0.0519 | 0.2279 | 0.1762 | 0.1435 | 0.1061 | 0.07866 | 0.1286 | 0.1382 | 0.02018 | 0.001998 | 0.001254 | 0.01692 | 0.02357 |
| 200 | s0.5 | 0.2611 | 0.0519 | 0.2394 | 0.1889 | 0.1536 | 0.1092 | 0.07783 | 0.1337 | 0.1446 | 0.019 | 0.001536 | 0.000539 | 0.01692 | 0.02213 |
| 200 | s0.9 | 0.2611 | 0.0519 | 0.226 | 0.1707 | 0.1342 | 0.09117 | 0.06838 | 0.1164 | 0.1272 | 0.01863 | 0.001356 | 0.0003501 | 0.01692 | 0.02173 |
| 200 | gauss | 0.2611 | 0.0519 | 0.2527 | 0.2088 | 0.1759 | 0.1305 | 0.08699 | 0.1536 | 0.1643 | 0.01903 | 0.001528 | 0.0005803 | 0.01692 | 0.02216 |

Per target (table shape): K_w of p+O alone and p-p alone, and the conditional removal probability K_w/P_el.

| E | target | K_1mm | K_2mm | K_3mm | K_5mm | K_2mm/P_el | K_3mm/P_el |
|---|---|---|---|---|---|---|---|
| 50 | O-16 | 0.03556 | 0.01772 | 0.01019 | 0.005649 | 0.2442 | 0.1404 |
| 50 | H-1 | 0.01085 | 0.008566 | 0.006823 | 0.004463 | 0.7115 | 0.5668 |
| 50 | C-12 | - | - | - | - | 0.3387 | 0.2 |
| 100 | O-16 | 0.1108 | 0.07579 | 0.0547 | 0.0326 | 0.4487 | 0.3238 |
| 100 | H-1 | 0.02373 | 0.02136 | 0.01945 | 0.01656 | 0.8564 | 0.7798 |
| 100 | C-12 | - | - | - | - | 0.529 | 0.3961 |
| 150 | O-16 | 0.1631 | 0.1211 | 0.09302 | 0.05992 | 0.533 | 0.4096 |
| 150 | H-1 | 0.03562 | 0.03426 | 0.0325 | 0.02957 | 0.947 | 0.8982 |
| 150 | C-12 | - | - | - | - | 0.6129 | 0.493 |
| 200 | O-16 | 0.1968 | 0.1519 | 0.1202 | 0.08064 | 0.5817 | 0.4604 |
| 200 | H-1 | 0.05165 | 0.05089 | 0.04987 | 0.04725 | 0.9807 | 0.9609 |
| 200 | C-12 | - | - | - | - | 0.666 | 0.5542 |

## Pointwise shape moments at the entrance energy (table; per interaction)

| E | target | <1-mu_CM> | theta_rms,CM | theta_med,lab | theta_2mm | P(th>th_2mm) | theta_3mm | P(th>th_3mm) | T (3 th_med) | <T_other> MeV | P(D>2mm) | P(D>3mm) | first min (bin density) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 50 | O-16 | 0.1592 | 35.6 | 20.10 | 24.38 | 0.3628 | 29.98 | 0.2246 | 0.08043 | 0.9056 | 0.4331 | 0.2852 | 56.6 |
| 50 | C-12 | 0.1783 | 36.8 | 22.47 | 24.38 | 0.4424 | 29.98 | 0.2944 | 0.07443 | 1.298 | 0.5312 | 0.3823 | 62.6 |
| 50 | H-1 | 1 | 96.6 | 30.76 | 24.38 | 0.7012 | 29.98 | 0.5264 | 0 | 13.32 | 0.9756 | 0.9385 | - |
| 100 | O-16 | 0.1176 | 30.0 | 15.76 | 13.03 | 0.6221 | 15.98 | 0.4907 | 0.06675 | 1.365 | 0.6802 | 0.5596 | 41.9 |
| 100 | C-12 | 0.1568 | 35.7 | 18.55 | 13.03 | 0.71 | 15.98 | 0.5972 | 0.08463 | 2.324 | 0.7695 | 0.6738 | 48.4 |
| 100 | H-1 | 1 | 97.1 | 29.91 | 13.03 | 0.9365 | 15.98 | 0.8818 | 0 | 26.01 | 1 | 1 | - |
| 150 | O-16 | 0.09283 | 26.8 | 12.44 | 9.11 | 0.6885 | 11.16 | 0.5718 | 0.05934 | 1.647 | 0.7393 | 0.6353 | 33.4 |
| 150 | C-12 | 0.1145 | 30.0 | 14.62 | 9.11 | 0.7637 | 11.16 | 0.667 | 0.06638 | 2.591 | 0.8145 | 0.7344 | 40.0 |
| 150 | H-1 | 1 | 97.3 | 29.57 | 9.11 | 0.9844 | 11.16 | 0.9531 | 0 | 38.71 | 1 | 1 | - |
| 200 | O-16 | 0.08216 | 25.5 | 10.68 | 7.10 | 0.7349 | 8.70 | 0.6299 | 0.05621 | 1.98 | 0.7798 | 0.688 | 28.5 |
| 200 | C-12 | 0.09761 | 27.6 | 12.53 | 7.10 | 0.7998 | 8.70 | 0.7153 | 0.06033 | 2.996 | 0.8442 | 0.7754 | 34.1 |
| 200 | H-1 | 1 | 97.3 | 29.27 | 7.10 | 1 | 8.70 | 0.9883 | 0 | 51.62 | 1 | 1 | - |

## Pointwise: p+O shape families

| E | shape | <1-mu_CM> | theta_rms,CM | theta_med,lab | P(th>th_2mm) | P(th>th_3mm) | T | P(D>2mm) | P(D>3mm) |
|---|---|---|---|---|---|---|---|---|---|
| 50 | table | 0.1592 | 35.6 | 20.10 | 0.3628 | 0.2246 | 0.08043 | 0.4331 | 0.2852 |
| 50 | disk | 0.1594 | 35.7 | 20.10 | 0.3628 | 0.2246 | 0.08057 | 0.4331 | 0.2852 |
| 50 | Rx1.2 | 0.131 | 31.8 | 16.89 | 0.2476 | 0.1475 | 0.07401 | 0.311 | 0.186 |
| 50 | s0.5 | 0.1041 | 27.4 | 18.62 | 0.3008 | 0.1621 | 0.03097 | 0.374 | 0.2217 |
| 50 | s0.9 | 0.06965 | 21.7 | 16.82 | 0.2246 | 0.09912 | 0.005282 | 0.2964 | 0.1514 |
| 50 | gauss | 0.1096 | 27.4 | 21.12 | 0.3989 | 0.2529 | 0.002558 | 0.4688 | 0.3188 |
| 100 | table | 0.1176 | 30.0 | 15.76 | 0.6221 | 0.4907 | 0.06675 | 0.6802 | 0.5596 |
| 100 | disk | 0.1176 | 30.0 | 15.76 | 0.6221 | 0.4907 | 0.06672 | 0.6802 | 0.5596 |
| 100 | Rx1.2 | 0.1028 | 28.4 | 13.25 | 0.5112 | 0.3691 | 0.06574 | 0.5791 | 0.4414 |
| 100 | s0.5 | 0.06021 | 20.3 | 14.25 | 0.561 | 0.416 | 0.01596 | 0.6265 | 0.4917 |
| 100 | s0.9 | 0.03933 | 16.2 | 12.70 | 0.4814 | 0.3281 | 0.001798 | 0.5537 | 0.4067 |
| 100 | gauss | 0.06662 | 21.2 | 16.36 | 0.6431 | 0.5161 | 0.00119 | 0.6987 | 0.5835 |
| 150 | table | 0.09283 | 26.8 | 12.44 | 0.6885 | 0.5718 | 0.05934 | 0.7393 | 0.6353 |
| 150 | disk | 0.09217 | 26.6 | 12.44 | 0.6885 | 0.5718 | 0.0586 | 0.7393 | 0.6353 |
| 150 | Rx1.2 | 0.0776 | 24.5 | 10.42 | 0.5879 | 0.4531 | 0.05293 | 0.6504 | 0.5249 |
| 150 | s0.5 | 0.03804 | 16.0 | 11.18 | 0.6338 | 0.501 | 0.01035 | 0.6919 | 0.5728 |
| 150 | s0.9 | 0.02472 | 12.8 | 10.00 | 0.5645 | 0.4194 | 0.001308 | 0.6304 | 0.4966 |
| 150 | gauss | 0.041 | 16.5 | 12.78 | 0.7021 | 0.5889 | 0.0006158 | 0.751 | 0.6509 |
| 200 | table | 0.08216 | 25.5 | 10.68 | 0.7349 | 0.6299 | 0.05621 | 0.7798 | 0.688 |
| 200 | disk | 0.08178 | 25.3 | 10.68 | 0.7349 | 0.6299 | 0.05585 | 0.7798 | 0.688 |
| 200 | Rx1.2 | 0.06722 | 22.9 | 8.93 | 0.644 | 0.5181 | 0.04805 | 0.7007 | 0.5869 |
| 200 | s0.5 | 0.02793 | 13.7 | 9.54 | 0.6846 | 0.5635 | 0.007518 | 0.7368 | 0.6304 |
| 200 | s0.9 | 0.01814 | 11.0 | 8.54 | 0.6216 | 0.4868 | 0.0009539 | 0.6821 | 0.5605 |
| 200 | gauss | 0.03004 | 14.1 | 10.89 | 0.7446 | 0.6426 | 0.0004466 | 0.7881 | 0.6992 |

## Reconciliation of the planning expectation (Amendment 14 (f): 0.985 / 0.986) with the observed ratio

Observed V5 ratio ionmc/TOPAS peak/plateau: 0.9634 / 0.9511 (step 22, -3.66 % / -4.89 %); V11 F ratio below. Kernel prediction F = (1 - K*_shifted)/(1 + G), O + p-p, table shape. The as-frozen K* is shown for the record only (sign error).

| quantity | 150 MeV | 200 MeV |
|---|---|---|
| F ionmc (V11, archived) | 0.8568 | 0.8441 |
| F TOPAS (V11 rows 10/11, archived) | 0.8916 | 0.8874 |
| F ionmc / F TOPAS | 0.961 | 0.9512 |
| peak removal ionmc = 1 - F(1+gain) | 0.1243 | 0.1334 |
| peak removal TOPAS = 1 - F(1+gain) | 0.09392 | 0.09629 |
| planning expectation ratio | 0.985 | 0.986 |
| planning-implied F ionmc (= ratio x F TOPAS) | 0.8782 | 0.8749 |
| gap: planning-implied F - observed F ionmc | 0.02143 | 0.03088 |
| K* as frozen (sign error), table (no conclusion) | 0.08366 | 0.08962 |
| K* shifted, table (O + pp) | 0.1356 | 0.1506 |
| G total, table / incl. post-p-p stopping | 0.0134 / 0.0194 | 0.0205 / 0.0240 |
| gain ionmc observed | 0.0221 | 0.02672 |
| F predicted, table (shifted K*) | 0.853 | 0.8323 |
| F predicted / F observed ionmc | 0.9956 | 0.9861 |
| peak removal predicted / observed | 1.091 | 1.129 |
| R_O table (fm) / planning R_O (fm) | 2.681 / 2.6 | 2.687 / 2.64 |
| RADIUS: F predicted at planning R_O (disk) | 0.8496 | 0.8302 |
| RADIUS: dF vs table | -0.003457 | -0.00215 |
| SHAPE: s (fm) giving the planning theta_rms,CM (17 / 14 deg) | 0.4379 | 0.4749 |
| SHAPE: F predicted at that s | 0.8654 | 0.8489 |
| SHAPE: dF vs table | 0.01239 | 0.01661 |
| radius share of the gap (planning F - observed F) | -0.1613 | -0.06962 |
| shape share of the gap (planning F - observed F) | 0.5781 | 0.5377 |
| shape share of the gap to the kernel's own table prediction | 0.4918 | 0.3895 |

## Reviewer's crude removal estimate (150 / 200 MeV, p+O only): cited against recomputed

Cited (review.md B.1, 150 MeV): moment, removal 2 mm, 3 mm. Recomputed with the same definitions (`crude_reviewer_estimate`: Bragg-Kleeman, straight line, step threshold) and with the plan definitions (project Bethe range, exact kinematics, table-interpolated radius, K_w of the O channel).

| E | shape | cited moment / 2 mm / 3 mm | crude recomputed | plan K_w O-only 2 mm / 3 mm | relative change vs sharp: crude 2mm / 3mm / moment | plan 2mm / 3mm |
|---|---|---|---|---|---|---|
| 150 | sharp | 0.0305/0.119/0.091 | 0.0305 / 0.1186 / 0.0909 | 0.1210 / 0.0930 | +0.0% / +0.0% / +0.0% | +0.0% / +0.0% |
| 150 | s0.5 | 0.0182/0.107/0.077 | 0.0182 / 0.1065 / 0.0769 | 0.1088 / 0.0789 | -10.2% / -15.4% / -40.2% | -10.1% / -15.2% |
| 150 | s0.9 | 0.0123/0.091/0.061 | 0.0123 / 0.0906 / 0.0606 | 0.0927 / 0.0623 | -23.6% / -33.3% / -59.7% | -23.4% / -33.0% |
| 150 | Rx1.2 | 0.025/0.096/0.069 | 0.0250 / 0.0956 / 0.0693 | 0.0978 / 0.0710 | -19.4% / -23.8% / -17.9% | -19.2% / -23.6% |
| 150 | gauss | 0.0196/0.125/0.096 | 0.0196 / 0.1246 / 0.0962 | 0.1270 / 0.0983 | +5.0% / +5.8% / -35.6% | +4.9% / +5.7% |
| 200 | sharp | - | 0.0319 / 0.1494 / 0.1180 | 0.1519 / 0.1202 | +0.0% / +0.0% / +0.0% | +0.0% / +0.0% |
| 200 | s0.5 | - | 0.0178 / 0.1356 / 0.1017 | 0.1380 / 0.1037 | -9.2% / -13.8% / -44.2% | -9.1% / -13.7% |
| 200 | s0.9 | - | 0.0118 / 0.1176 / 0.0826 | 0.1198 / 0.0844 | -21.3% / -30.1% / -63.0% | -21.1% / -29.8% |
| 200 | Rx1.2 | - | 0.0261 / 0.1230 / 0.0918 | 0.1253 / 0.0937 | -17.6% / -22.2% / -18.1% | -17.5% / -22.1% |
| 200 | gauss | - | 0.0190 / 0.1554 / 0.1238 | 0.1579 / 0.1260 | +4.0% / +4.9% / -40.5% | +4.0% / +4.8% |

## H6: TOPAS rows 10/11 (emelastic) against the frozen EM-only cases (emonly)

`input.txt` differs only in the header comment, the seed and the module list for all 6 pairs: **True**.

Differing lines of the seed-1 pairs (emelastic `+` against emonly `-`):

- proton-water-150mev-idd-r20-emelastic-seed1 vs proton-water-150mev-idd-r20-emonly-seed1:
  - `-# Native OpenTOPAS 4.3 input. V3-005B row V5 absolute integral depth dose (r = 20 cm), EM-ONLY physics; no validation claim.`
  - `+# Native OpenTOPAS 4.3 input. V3-005C plan Amendment 14 (f) row 10 (emelastic, 150 MeV) for acceptance row V11: EM + hadronic elastic only (no inelas`
  - `-i:Ts/Seed = 20273001`
  - `+i:Ts/Seed = 20280001`
  - `-sv:Ph/Default/Modules = 1 "g4em-standard_opt4"`
  - `+sv:Ph/Default/Modules = 2 "g4em-standard_opt4" "g4h-elastic_HP"`
- proton-water-200mev-idd-r20-emelastic-seed1 vs proton-water-200mev-idd-r20-emonly-seed1:
  - `-# Native OpenTOPAS 4.3 input. V3-005B row V5 absolute integral depth dose (r = 20 cm), EM-ONLY physics; no validation claim.`
  - `+# Native OpenTOPAS 4.3 input. V3-005C plan Amendment 14 (f) row 11 (emelastic, 200 MeV) for acceptance row V11: EM + hadronic elastic only (no inelas`
  - `-i:Ts/Seed = 20274001`
  - `+i:Ts/Seed = 20281001`
  - `-sv:Ph/Default/Modules = 1 "g4em-standard_opt4"`
  - `+sv:Ph/Default/Modules = 2 "g4em-standard_opt4" "g4h-elastic_HP"`

`case.json` differs in the descriptive text, row/mode/modules and the seed (not scored by the engine). Hashed inputs are listed above.
