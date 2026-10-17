Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis: combined per-case summaries

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --combine`

Input files (sha256):

- `validation/generated/elastic_attribution/de4/de4-case0.json`  d9aedb8c59a7fcdfbdb8dd312c3a24a0e4b3bc6025ac46daf562782797a8e03c
- `validation/generated/elastic_attribution/de4/de4-case1.json`  7fdaa24c0b1fdc220eb14bcc22ae5429193bf0c32c4061cd6cb8b1067daa92fb
- `validation/generated/elastic_attribution/de4/de4-case2.json`  b630857ee99c9c97284835370e0f1419540c4eab82e341f89574ed19320cb807
- `validation/generated/elastic_attribution/de4/de4-case3.json`  fbf53800b4a3a519a71d303b4d434bd351b31588ab5cff3742ddfd76af4c7d24
- `validation/generated/elastic_attribution/de4/de4-case4.json`  413119bae1c2729326f77f447197236a0ca4fa60508b1b1e23157c7f82db3c3c
- `validation/generated/elastic_attribution/de4/de4-case5.json`  b7fec574b07f4d1440fa6243266eb1c015730c2e30730f9a54ee7e7261908d49
- `validation/generated/elastic_attribution/de4/de4-case6.json`  af2a6149bd636b8e674a99a6bab70f33f6036cfa922a24b317fa9f9cb58eea79
- `validation/generated/elastic_attribution/de4/de4-case7.json`  d5114e13e47ea5f57cbd296a49ef26aeaa46581c4973fc95aa657f6124fe2a0b
- `validation/generated/elastic_attribution/de4/de4-case8.json`  3c029c0ec3dddd1013fe552c6e6ea6a6b8b89d1c1a453dd6a19eb01c7cd75c51
- `validation/generated/elastic_attribution/de4/de4-case9.json`  fefd70f6be67d43af36f1a67d5e54700706371b851ad28a768ed1fd88f15bcff

<!-- body below is rendered from the Git projection of the JSON record -->

## Case 0

Case 0 (row 14: water, 30 MeV, 0.2 cm, 20000000 primaries, 20025531 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 19833617 primaries removed (fraction 0.99168); 20000000 primary records, 25531 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 25325); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 24103 | 0.001205 +- 7.8e-06 | 30.42 | 0.2027 | 0.0002 | 0.0000 | 7.749 |
| O-16 | 140322 | 0.007016 +- 1.9e-05 | 22.52 | 0.1108 | 0.0722 | 0.0345 | 0.486 |

p-p events: 25325 with a secondary proton; 234 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 988 removed by the 16.26 deg CM cut; 24103 kept (class H-1).

Classification: 166383 deflected primaries -> 141063 p+A, 25320 p-p; off both kinematic curves: 731 p+A-pool tracks excluded (fraction 5.182e-03; 1071 including primaries of p-p events); raw conflict tracks 527 (3.167e-03; crossover 1.86 deg), of which 25 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 2.39e+02, fraction < 1e-3 0.9790.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 1

Case 1 (row 15: water, 50 MeV, 0.2 cm, 20000000 primaries, 20014038 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 19877083 primaries removed (fraction 0.99385); 20000000 primary records, 14038 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 13968); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 13351 | 0.0006676 +- 5.8e-06 | 30.44 | 0.6790 | 0.5140 | 0.0000 | 13.047 |
| O-16 | 108606 | 0.00543 +- 1.6e-05 | 16.81 | 0.2171 | 0.0988 | 0.0111 | 0.425 |

p-p events: 13968 with a secondary proton; 84 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 533 removed by the 16.26 deg CM cut; 13351 kept (class H-1).

Classification: 122917 deflected primaries -> 108937 p+A, 13980 p-p; off both kinematic curves: 337 p+A-pool tracks excluded (fraction 3.093e-03; 483 including primaries of p-p events); raw conflict tracks 865 (7.037e-03; crossover 1.85 deg), of which 9 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 2.33e+02, fraction < 1e-3 0.9848.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 2

Case 2 (row 16: water, 100 MeV, 1 cm, 4000000 primaries, 4007563 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3945499 primaries removed (fraction 0.98637); 4000000 primary records, 7563 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 7473); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 7078 | 0.00177 +- 2.1e-05 | 29.61 | 0.9294 | 0.8751 | 0.0000 | 25.727 |
| O-16 | 46790 | 0.0117 +- 5.4e-05 | 11.44 | 0.4044 | 0.2575 | 0.0047 | 0.418 |

p-p events: 7473 with a secondary proton; 134 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 261 removed by the 16.26 deg CM cut; 7078 kept (class H-1).

Classification: 54501 deflected primaries -> 47021 p+A, 7480 p-p; off both kinematic curves: 233 p+A-pool tracks excluded (fraction 4.954e-03; 451 including primaries of p-p events); raw conflict tracks 784 (1.439e-02; crossover 1.83 deg), of which 3 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 8.86e+02, fraction < 1e-3 0.9553.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 3

Case 3 (row 17: water, 150 MeV, 1 cm, 4000000 primaries, 4006429 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3972614 primaries removed (fraction 0.99315); 4000000 primary records, 6429 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 6370); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 6028 | 0.001507 +- 1.9e-05 | 29.96 | 0.9877 | 0.9580 | 0.0000 | 39.083 |
| O-16 | 20974 | 0.005243 +- 3.6e-05 | 9.30 | 0.5132 | 0.3676 | 0.0036 | 0.434 |

p-p events: 6370 with a secondary proton; 55 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 287 removed by the 16.26 deg CM cut; 6028 kept (class H-1).

Classification: 27386 deflected primaries -> 21015 p+A, 6371 p-p; off both kinematic curves: 39 p+A-pool tracks excluded (fraction 1.856e-03; 179 including primaries of p-p events); raw conflict tracks 570 (2.081e-02; crossover 1.81 deg), of which 6 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 8.27e+02, fraction < 1e-3 0.9631.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 4

Case 4 (row 18: water, 200 MeV, 1 cm, 4000000 primaries, 4006265 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3980303 primaries removed (fraction 0.99508); 4000000 primary records, 6265 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 6214); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 5919 | 0.00148 +- 1.9e-05 | 28.70 | 1.0000 | 0.9856 | 0.0000 | 50.508 |
| O-16 | 13464 | 0.003366 +- 2.9e-05 | 7.93 | 0.5711 | 0.4340 | 0.0022 | 0.430 |

p-p events: 6214 with a secondary proton; 41 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 254 removed by the 16.26 deg CM cut; 5919 kept (class H-1).

Classification: 19697 deflected primaries -> 13484 p+A, 6213 p-p; off both kinematic curves: 18 p+A-pool tracks excluded (fraction 1.335e-03; 142 including primaries of p-p events); raw conflict tracks 505 (2.564e-02; crossover 1.78 deg), of which 10 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 3.33e+02, fraction < 1e-3 0.9733.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 5

Case 5 (row 19: graphite, 30 MeV, 0.2 cm, 20000000 primaries, 20000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 19654435 primaries removed (fraction 0.98272); 20000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 340490 | 0.01702 +- 2.9e-05 | 20.77 | 0.0776 | 0.0471 | 0.0310 | 0.546 |

p-p events: 0 with a secondary proton; 0 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 345565 deflected primaries -> 345291 p+A, 274 p-p; off both kinematic curves: 5051 p+A-pool tracks excluded (fraction 1.462e-02; 5051 including primaries of p-p events); raw conflict tracks 1622 (4.694e-03; crossover 1.88 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 6

Case 6 (row 20: graphite, 50 MeV, 0.2 cm, 20000000 primaries, 20000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 19752729 primaries removed (fraction 0.98764); 20000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 245121 | 0.01226 +- 2.5e-05 | 16.09 | 0.1903 | 0.0806 | 0.0097 | 0.514 |

p-p events: 0 with a secondary proton; 0 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 247271 deflected primaries -> 247154 p+A, 117 p-p; off both kinematic curves: 2131 p+A-pool tracks excluded (fraction 8.618e-03; 2131 including primaries of p-p events); raw conflict tracks 2126 (8.598e-03; crossover 1.87 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 7

Case 7 (row 21: graphite, 100 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3900305 primaries removed (fraction 0.97508); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 98446 | 0.02461 +- 7.8e-05 | 11.55 | 0.4123 | 0.2627 | 0.0045 | 0.559 |

p-p events: 0 with a secondary proton; 0 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 99695 deflected primaries -> 99586 p+A, 109 p-p; off both kinematic curves: 1231 p+A-pool tracks excluded (fraction 1.235e-02; 1231 including primaries of p-p events); raw conflict tracks 1640 (1.645e-02; crossover 1.85 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 8

Case 8 (row 22: graphite, 150 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3950527 primaries removed (fraction 0.98763); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 49206 | 0.0123 +- 5.5e-05 | 9.49 | 0.5262 | 0.3833 | 0.0031 | 0.589 |

p-p events: 0 with a secondary proton; 0 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 49473 deflected primaries -> 49448 p+A, 25 p-p; off both kinematic curves: 263 p+A-pool tracks excluded (fraction 5.316e-03; 263 including primaries of p-p events); raw conflict tracks 1238 (2.502e-02; crossover 1.83 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 9

Case 9 (row 23: graphite, 200 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3966450 primaries removed (fraction 0.99161); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 33435 | 0.008359 +- 4.6e-05 | 8.23 | 0.5992 | 0.4625 | 0.0023 | 0.608 |

p-p events: 0 with a secondary proton; 0 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 33550 deflected primaries -> 33536 p+A, 14 p-p; off both kinematic curves: 112 p+A-pool tracks excluded (fraction 3.338e-03; 112 including primaries of p-p events); raw conflict tracks 1071 (3.192e-02; crossover 1.80 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
