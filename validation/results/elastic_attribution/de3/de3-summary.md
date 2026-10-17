Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E3 probes: summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de3_probes --summary`

Input files (sha256):

- `validation/generated/elastic_attribution/de3/de3-150-all.json`  a3443dcaef1007198425460274cae4266ed61a0b35602306208368cfad274d53
- `validation/generated/elastic_attribution/de3/de3-150-emonly.json`  284b081cdb1517c2647b33169ed5e1ae08457d4563de135b77aa7ec550832876
- `validation/generated/elastic_attribution/de3/de3-150-pA-Rx1.2.json`  1beb630a733e4fd17f12564e7c6b66dc7caa33c9c231394ae3fb5c88b0528260
- `validation/generated/elastic_attribution/de3/de3-150-pA-gauss.json`  6525184c6c3c687a7c5455375d8b33ef5f1e01aabb2104ec6b9d7589f4df2288
- `validation/generated/elastic_attribution/de3/de3-150-pA-only.json`  e3a859b3e59aef608a3343d67e9b6bba9269efa7f7c84648eff218616406431a
- `validation/generated/elastic_attribution/de3/de3-150-pA-s0.5.json`  be577ccc52c2e7efd9aea0f648c6187b5dba160825843591d15fa1a4320f41d2
- `validation/generated/elastic_attribution/de3/de3-150-pA-s0.9.json`  0997bfeb5fbedad92dfc41134532f9215933b06025c9e92d9108b700ed3efdb1
- `validation/generated/elastic_attribution/de3/de3-150-pA-sigma0.8.json`  c184b276ebee230573dd1cfb29694f562066b8daca3b95ebe524042bcb25df5a
- `validation/generated/elastic_attribution/de3/de3-150-pinned.json`  3abc5c3978faa69232af72b0359d8b8378a6b4032738a7d975c19df32233bda1
- `validation/generated/elastic_attribution/de3/de3-150-pp-only.json`  8ad70e1a67cee8d1d9361a4ccc7078c531d92850336e8bad349272ecbc9100f6
- `validation/generated/elastic_attribution/de3/de3-200-all.json`  1358a0b1addeb161e0ed79a859e8cd10ac17584c034bd38ab59096408be0aae1
- `validation/generated/elastic_attribution/de3/de3-200-emonly.json`  29fbcdceb47d0a839bb6c0753f8b9ed6c8292c428b6f899723b580ed04f9b3cc
- `validation/generated/elastic_attribution/de3/de3-200-pA-Rx1.2.json`  a7801a7fed917121dd3b79c321f836a26140dd0564bc4a68e29317378634296d
- `validation/generated/elastic_attribution/de3/de3-200-pA-gauss.json`  4f45f0faffc54f44a40067855064bf3f30913f3aa498e1daa02bb9d83ccbb1c5
- `validation/generated/elastic_attribution/de3/de3-200-pA-only.json`  3bc5fbb8d1d4dad11c98d6c2673b7013e72e1c4d5a626dd3d9d04ae4b072943b
- `validation/generated/elastic_attribution/de3/de3-200-pA-s0.5.json`  1fd760865a7d22c47cd1e4919e245a975bd6d97c68cb77a7f352caad926b40e1
- `validation/generated/elastic_attribution/de3/de3-200-pA-s0.9.json`  21792b3c460174b335e2de1f8b2f139b61e2bb014c3bf6c8262073c0605a5e3d
- `validation/generated/elastic_attribution/de3/de3-200-pA-sigma0.8.json`  9d3c5aff1eceed5806d8fb007730ed1a7b4f2d36ba8dee7fd01e52962f521d57
- `validation/generated/elastic_attribution/de3/de3-200-pinned.json`  23b897ff060fab214d388f64a62e8c2bd0c74a4e07295ea6af242ed2441bb1db
- `validation/generated/elastic_attribution/de3/de3-200-pp-only.json`  fb4f6fe3bdf088d780efbb4e05c9f2699468c5bced34ea4fa4569be795f7e50e

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

- recorded git HEAD d2bce93c1f521b89ef7fd1176a9fe7a5ac927a9e; recorded status of validation/scripts/elastic_attribution: dirty: M validation/scripts/elastic_attribution/de3_probes.py
- script sha256 recorded at summary time: common.py 29171edbda46e951, de0_kernels.py 38e5ec2dbe285af1, de3_probes.py c4e01676fc3e2840, make_variant_tables.py beafad7d90f4a4d6
- script sha256 recorded for the run (staged copies, if any): none recorded
- table files sha256: elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.json ab5e2a357ab53060, elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.npz 1028bb67f57bc96e
- backend warp-cpu float64, single process True; seeds {'150': 20990100, '200': 20990101}; runs (seed, n, batches, table id, wall s) in `runs` of de3-summary.json
- the run files were produced with the scripts whose hashes are recorded in the logged run metadata: 20261015T113620Z-de3-emonly-150, 20261015T113649Z-de3-pinned-150, 20261015T113747Z-de3-all-150, 20261015T113844Z-de3-pp-only-150, 20261015T113941Z-de3-pA-only-150, 20261015T114038Z-de3-pA-sigma0.8-150, 20261015T114134Z-de3-pA-Rx1.2-150, 20261015T114231Z-de3-pA-s0.5-150, 20261015T114328Z-de3-pA-s0.9-150, 20261015T114424Z-de3-pA-gauss-150, 20261015T114521Z-de3-emonly-200, 20261015T114600Z-de3-pinned-200, 20261015T114716Z-de3-all-200, 20261015T114832Z-de3-pp-only-200, 20261015T114947Z-de3-pA-only-200, 20261015T115101Z-de3-pA-sigma0.8-200, 20261015T115216Z-de3-pA-Rx1.2-200, 20261015T115331Z-de3-pA-s0.5-200, 20261015T115446Z-de3-pA-s0.9-200, 20261015T115601Z-de3-pA-gauss-200
