Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E3 probes: summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de3_probes --summary`

Input files (sha256):

- `validation/generated/elastic_attribution/de3/de3-150-all.json`  81c883b94d67022108b2ff540cf53bf5ecdb8e8b8d937193946988ee273544ac
- `validation/generated/elastic_attribution/de3/de3-150-emonly.json`  63c7af6022849a8981c54c3267ce6c11ed0db66525856c133181941f0a607762
- `validation/generated/elastic_attribution/de3/de3-150-pA-Rx1.2.json`  4652a77fc9001f6cc91820a6995a748dbc7ca1410ea7848251f370ee2a098f5e
- `validation/generated/elastic_attribution/de3/de3-150-pA-gauss.json`  23ddcde95d89cc722299dd7f2794bddbb1484a7ffef1d72230d3b30880240adb
- `validation/generated/elastic_attribution/de3/de3-150-pA-only.json`  37d0e9becf78430a534fddfecd8ce3c036e35006486ce32880e201d4f4d6012e
- `validation/generated/elastic_attribution/de3/de3-150-pA-s0.5.json`  9b8a168ea558dc974d24ba8fe4f845d2d313e01ff299306db555cd00b96e765c
- `validation/generated/elastic_attribution/de3/de3-150-pA-s0.9.json`  d12783d51341c2449fb3a4ca6c11fdd50f4ce4ac2b5a50a010a12ac29dbcef80
- `validation/generated/elastic_attribution/de3/de3-150-pA-sigma0.8.json`  e5380f9705ad1b4bafad85d56b5e11f1d59eb3eb6f064ae6f184266f3f9020b4
- `validation/generated/elastic_attribution/de3/de3-150-pinned.json`  9d4aef601df5b5217c0d2b7872b8888c98d4ae5e54ecd6dd48401e402d176edb
- `validation/generated/elastic_attribution/de3/de3-150-pp-only.json`  a321da0ee2372a7fed115fdd7569db0c9f477c9db010c667660e3fc3658db27c
- `validation/generated/elastic_attribution/de3/de3-200-all.json`  a93432513cc74a1954fce65ccaab7ad6a8f613dd1212659b09ac506ed1b71579
- `validation/generated/elastic_attribution/de3/de3-200-emonly.json`  a85cf343e6b3557d2a26e31ea47a46f1a0d28e648526100aa1c663e7fec724e9
- `validation/generated/elastic_attribution/de3/de3-200-pA-Rx1.2.json`  93a75e752dabc90a97a4f14fb623d3d21cfbd2501d460af809d527b5b6c6805e
- `validation/generated/elastic_attribution/de3/de3-200-pA-gauss.json`  09e3e40fdc3ca654a93a843df4d22d8e8ee4cb8f3d70797efba4767007a1f692
- `validation/generated/elastic_attribution/de3/de3-200-pA-only.json`  af60be8492efa17d12afaa49cd71edfae919c5addc32df711f641284f80d1196
- `validation/generated/elastic_attribution/de3/de3-200-pA-s0.5.json`  df418c707122746aebe0613cd1bccb3c04853b8b5b0001749fc8e6b5369a198c
- `validation/generated/elastic_attribution/de3/de3-200-pA-s0.9.json`  e9e82a7c28ad19d8a501ed9e89e44045d217a6bfc3bccddef7aed0b4503947b6
- `validation/generated/elastic_attribution/de3/de3-200-pA-sigma0.8.json`  0ef609b4ef8d63e618ae4a9752a4040638bd6d3d84658d149175a661b180234a
- `validation/generated/elastic_attribution/de3/de3-200-pinned.json`  bee346f88e3b64cf05379110cb1f963ca8654693e15edf60c335a4861be23df8
- `validation/generated/elastic_attribution/de3/de3-200-pp-only.json`  f807e208d195793a7f49ca48bdfcf151f667c0afbb1c7346e51b835589ca4524

<!-- body below is rendered from the Git projection of the JSON record -->

| E | variant | F | gain | peak removal | K* shifted (pred) | pred/meas removal | F pred (meas. gain) | recoil MeV | R80 mm | R80 EM-only | lateral rms mm (EM+el / EM) | valid |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 150 | all | 0.8580 +- 0.0014 | 0.0221 +- 0.0008 | 0.1230 | 0.1356 | 1.102 | 0.8458 | 0.2678 | 158.323 | 158.456 | 5.515 / 1.540 | True |
| 150 | pA-Rx1.2 | 0.8757 +- 0.0016 | 0.0202 +- 0.0007 | 0.1067 | 0.1148 | 1.076 | 0.8677 | 0.2197 | 158.321 | 158.456 | 5.193 / 1.540 | True |
| 150 | pA-gauss | 0.8592 +- 0.0016 | 0.0163 +- 0.0006 | 0.1268 | 0.1401 | 1.105 | 0.8461 | 0.1667 | 158.324 | 158.456 | 4.841 / 1.540 | True |
| 150 | pA-only | 0.8986 +- 0.0014 | 0.0081 +- 0.0003 | 0.0942 | 0.1028 | 1.092 | 0.8900 | 0.2680 | 158.325 | 158.456 | 4.467 / 1.540 | True |
| 150 | pA-s0.5 | 0.8744 +- 0.0015 | 0.0162 +- 0.0006 | 0.1115 | 0.1225 | 1.098 | 0.8636 | 0.1499 | 158.314 | 158.456 | 4.690 / 1.540 | True |
| 150 | pA-s0.9 | 0.8872 +- 0.0015 | 0.0154 +- 0.0006 | 0.0992 | 0.1070 | 1.079 | 0.8794 | 0.1040 | 158.311 | 158.456 | 4.299 / 1.540 | True |
| 150 | pA-sigma0.8 | 0.8763 +- 0.0012 | 0.0203 +- 0.0007 | 0.1058 | 0.1150 | 1.087 | 0.8673 | 0.2123 | 158.344 | 158.456 | 5.138 / 1.540 | True |
| 150 | pinned | 0.8580 +- 0.0014 | 0.0221 +- 0.0008 | 0.1230 | 0.1356 | 1.102 | 0.8458 | 0.2678 | 158.323 | 158.456 | 5.515 / 1.540 | True |
| 150 | pp-only | 0.9547 +- 0.0010 | 0.0137 +- 0.0006 | 0.0323 | 0.0327 | 1.014 | 0.9542 | 0.0000 | 158.449 | 158.456 | 3.645 / 1.540 | True |
| 200 | all | 0.8455 +- 0.0027 | 0.0249 +- 0.0010 | 0.1334 | 0.1506 | 1.129 | 0.8287 | 0.3445 | 260.435 | 260.750 | 7.301 / 2.468 | True |
| 200 | pA-Rx1.2 | 0.8651 +- 0.0030 | 0.0239 +- 0.0010 | 0.1143 | 0.1286 | 1.126 | 0.8510 | 0.2831 | 260.444 | 260.750 | 7.119 / 2.468 | True |
| 200 | pA-gauss | 0.8457 +- 0.0028 | 0.0217 +- 0.0009 | 0.1359 | 0.1536 | 1.130 | 0.8285 | 0.2006 | 260.429 | 260.750 | 7.269 / 2.468 | True |
| 200 | pA-only | 0.9044 +- 0.0027 | 0.0048 +- 0.0005 | 0.0913 | 0.1025 | 1.122 | 0.8932 | 0.3452 | 260.463 | 260.750 | 4.890 / 2.468 | True |
| 200 | pA-s0.5 | 0.8638 +- 0.0029 | 0.0215 +- 0.0010 | 0.1176 | 0.1337 | 1.137 | 0.8480 | 0.1806 | 260.414 | 260.750 | 7.002 / 2.468 | True |
| 200 | pA-s0.9 | 0.8785 +- 0.0029 | 0.0212 +- 0.0009 | 0.1029 | 0.1164 | 1.132 | 0.8652 | 0.1262 | 260.419 | 260.750 | 6.786 / 2.468 | True |
| 200 | pA-sigma0.8 | 0.8614 +- 0.0031 | 0.0240 +- 0.0009 | 0.1179 | 0.1301 | 1.103 | 0.8495 | 0.2742 | 260.499 | 260.750 | 7.148 / 2.468 | True |
| 200 | pinned | 0.8455 +- 0.0027 | 0.0249 +- 0.0010 | 0.1334 | 0.1506 | 1.129 | 0.8287 | 0.3445 | 260.435 | 260.750 | 7.301 / 2.468 | True |
| 200 | pp-only | 0.9349 +- 0.0010 | 0.0198 +- 0.0005 | 0.0466 | 0.0481 | 1.033 | 0.9334 | 0.0000 | 260.728 | 260.750 | 6.190 / 2.468 | True |

## Derived (D-E3 review): composition and gap closure

| E | F_all | F_pA x F_pp | F_all - F_pA F_pp (paired jackknife) | removal additivity pA + pp - all (paired jackknife) | F_TOPAS (reconstructed) | gap closed: sigma0.8 / Rx1.2 / s0.5 / s0.9 / gauss |
|---|---|---|---|---|---|---|
| 150 | 0.8580 | 0.8578 | +0.00024 +- 0.00144 | +0.00343 +- 0.00139 | 0.8929 | 0.52 / 0.51 / 0.47 / 0.84 / 0.03 |
| 200 | 0.8455 | 0.8455 | +0.00005 +- 0.00146 | +0.00451 +- 0.00150 | 0.8889 | 0.36 / 0.45 / 0.42 / 0.76 / 0.00 |

## Findings recorded with the independent D-E3 review (exploratory; nothing enters the section 5 rule)

- **Composition:** 150 MeV: F_all - F_pA F_pp = +0.00024 +- 0.00144 (consistent with multiplicative within the paired uncertainty); removal additivity +0.00343 +- 0.00139 (NOT consistent with additive); 200 MeV: F_all - F_pA F_pp = +0.00005 +- 0.00146 (consistent with multiplicative within the paired uncertainty); removal additivity +0.00451 +- 0.00150 (NOT consistent with additive). Where it is consistent, D-E6 effects are to be read on -ln F; on F itself the product term would consume part of the 30 % interaction budget of section 5.
- **Gap closure toward TOPAS F** (F_TOPAS reconstructed from the archived V11 Delta_pk): pA-s0.9 closes most of the gap, pA-s0.5 and pA-Rx1.2 about half; the shape (core width) is the lever on F.
- **Contrary to a pure tail story:** pA-gauss removes the large-angle tail yet leaves F unchanged while lowering the gain; the peak deficit is driven by the core width, the plateau gain by the tail.
- **Contrary evidence on p-p (H3 stays live):** the pp-only plateau gain alone is of the order of TOPAS's total elastic gain (84 % / 108 % per the review, from the archived V11 gain); no p+A shape variant closes the gain excess, so the swap-pp arm of D-E6 is essential.
- **K* model:** the shifted K* overpredicts the measured peak removal by 8-14 % for every p+A variant and by 1-3 % for pp-only, the bias anticipated in the D-E0 dated note (adequate for ranking, not for a 0.010 closure).

## Provenance

- recorded git HEAD c7c8b662c268b2bcf8237a317a34369cd5905c52; recorded status of validation/scripts/elastic_attribution: dirty: M validation/scripts/elastic_attribution/de3_probes.py;  M validation/scripts/elastic_attribution/de4_parse.py;  M validation/scripts/elastic_attribution/de5_sampling.py
- script sha256 recorded at summary time: common.py 29171edbda46e951, de0_kernels.py 38e5ec2dbe285af1, de3_probes.py 257ae085db885c5e, make_variant_tables.py beafad7d90f4a4d6
- script sha256 recorded for the run (staged copies, if any): none recorded
- table files sha256: elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.json ab5e2a357ab53060, elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.npz 1028bb67f57bc96e
- backend warp-cpu float64, single process True; seeds {'150': 20990100, '200': 20990101}; runs (seed, n, batches, table id, wall s) in `runs` of de3-summary.json
- the run files were produced with the scripts whose hashes are recorded in the logged run metadata (log directory names not recorded in this record)
