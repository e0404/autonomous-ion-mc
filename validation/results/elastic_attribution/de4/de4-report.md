Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis: combined per-case summaries

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --combine`

Input files (sha256):

- `validation/generated/elastic_attribution/de4/de4-case0.json`  13fa9594a6f719c49835f5100c70b633bc2d4b5ea09cdf279fcb85462c3da254
- `validation/generated/elastic_attribution/de4/de4-case1.json`  3a2f98403e38bc87ae0f69598aee2ce223a5e4fad5cae5cc4da2d0c4ccab41bc
- `validation/generated/elastic_attribution/de4/de4-case2.json`  19d23663b416e790b2a052ebeea6e41c4cf82f8229e7d199ded2f5cac5c237a7
- `validation/generated/elastic_attribution/de4/de4-case3.json`  73eb2c4de870a6e71f1518f43e06e889a91716ca7081ff7b311ffb27de74f543
- `validation/generated/elastic_attribution/de4/de4-case4.json`  ed46a0d437c8e236c95dac99a87f6ff06394beda0242fc5ee0ecc8c70465a742
- `validation/generated/elastic_attribution/de4/de4-case5.json`  1d88619b05066999f54b31819136d0ca558e170b2917046201a0720075088a37
- `validation/generated/elastic_attribution/de4/de4-case6.json`  4a1e759037af65453736095ecfe652e061d1d4467042020037675613453fa926
- `validation/generated/elastic_attribution/de4/de4-case7.json`  95d0d7ec5237c322005d018d0cd78fe1a26855270303db4f143a051e1f02eef7
- `validation/generated/elastic_attribution/de4/de4-case8.json`  3e11200e2bec655758bc8a7852cf540827550013db87e84f8e1e40085431af30
- `validation/generated/elastic_attribution/de4/de4-case9.json`  7fc4cec8d2d92b3b9a3d01a38ec7a7cc66dcfa4b9335da1103c4fcf1f70e4c04

<!-- body below is rendered from the Git projection of the JSON record -->

## Case 0

Case 0 (row 14: water, 30 MeV, 0.2 cm, 20000000 primaries, 20025531 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 19833617 primaries removed (fraction 0.99168); 20000000 primary records, 25531 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 25325); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 23782 | 0.001189 +- 7.7e-06 | 30.44 | 0.2031 | 0.0000 | 0.0000 | 7.789 |
| O-16 | 140322 | 0.007016 +- 1.9e-05 | 22.52 | 0.1108 | 0.0722 | 0.0345 | 0.486 |

p-p events: 25325 with a secondary proton; 556 excluded (off the p+H curve or not exactly two proton records); 987 removed by the 16.26 deg CM cut; 23782 kept (class H-1).

Classification: 166383 deflected primaries -> 141063 p+A, 25320 p-p; off both kinematic curves: 731 p+A-pool tracks excluded (fraction 5.182e-03; 1071 including primaries of p-p events); raw conflict tracks 527 (3.167e-03; crossover 1.86 deg), of which 25 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 2.39e+02, fraction < 1e-3 0.9790.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 1

Case 1 (row 15: water, 50 MeV, 0.2 cm, 20000000 primaries, 20014038 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 19877083 primaries removed (fraction 0.99385); 20000000 primary records, 14038 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 13968); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 13179 | 0.000659 +- 5.7e-06 | 30.51 | 0.6817 | 0.5161 | 0.0000 | 13.111 |
| O-16 | 108606 | 0.00543 +- 1.6e-05 | 16.81 | 0.2171 | 0.0988 | 0.0111 | 0.425 |

p-p events: 13968 with a secondary proton; 256 excluded (off the p+H curve or not exactly two proton records); 533 removed by the 16.26 deg CM cut; 13179 kept (class H-1).

Classification: 122917 deflected primaries -> 108937 p+A, 13980 p-p; off both kinematic curves: 337 p+A-pool tracks excluded (fraction 3.093e-03; 483 including primaries of p-p events); raw conflict tracks 865 (7.037e-03; crossover 1.85 deg), of which 9 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 2.33e+02, fraction < 1e-3 0.9848.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 2

Case 2 (row 16: water, 100 MeV, 1 cm, 4000000 primaries, 4007563 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3945499 primaries removed (fraction 0.98637); 4000000 primary records, 7563 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 7473); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 6805 | 0.001701 +- 2.1e-05 | 29.75 | 0.9311 | 0.8773 | 0.0000 | 25.969 |
| O-16 | 46790 | 0.0117 +- 5.4e-05 | 11.44 | 0.4044 | 0.2575 | 0.0047 | 0.418 |

p-p events: 7473 with a secondary proton; 415 excluded (off the p+H curve or not exactly two proton records); 253 removed by the 16.26 deg CM cut; 6805 kept (class H-1).

Classification: 54501 deflected primaries -> 47021 p+A, 7480 p-p; off both kinematic curves: 233 p+A-pool tracks excluded (fraction 4.954e-03; 451 including primaries of p-p events); raw conflict tracks 784 (1.439e-02; crossover 1.83 deg), of which 3 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 8.86e+02, fraction < 1e-3 0.9553.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 3

Case 3 (row 17: water, 150 MeV, 1 cm, 4000000 primaries, 4006429 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3972614 primaries removed (fraction 0.99315); 4000000 primary records, 6429 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 6370); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 5798 | 0.00145 +- 1.9e-05 | 30.17 | 0.9884 | 0.9600 | 0.0000 | 39.514 |
| O-16 | 20974 | 0.005243 +- 3.6e-05 | 9.30 | 0.5132 | 0.3676 | 0.0036 | 0.434 |

p-p events: 6370 with a secondary proton; 294 excluded (off the p+H curve or not exactly two proton records); 278 removed by the 16.26 deg CM cut; 5798 kept (class H-1).

Classification: 27386 deflected primaries -> 21015 p+A, 6371 p-p; off both kinematic curves: 39 p+A-pool tracks excluded (fraction 1.856e-03; 179 including primaries of p-p events); raw conflict tracks 570 (2.081e-02; crossover 1.81 deg), of which 6 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 8.27e+02, fraction < 1e-3 0.9631.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 4

Case 4 (row 18: water, 200 MeV, 1 cm, 4000000 primaries, 4006265 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3980303 primaries removed (fraction 0.99508); 4000000 primary records, 6265 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 6214); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 5749 | 0.001437 +- 1.9e-05 | 28.84 | 1.0000 | 0.9863 | 0.0000 | 50.969 |
| O-16 | 13464 | 0.003366 +- 2.9e-05 | 7.93 | 0.5711 | 0.4340 | 0.0022 | 0.430 |

p-p events: 6214 with a secondary proton; 222 excluded (off the p+H curve or not exactly two proton records); 243 removed by the 16.26 deg CM cut; 5749 kept (class H-1).

Classification: 19697 deflected primaries -> 13484 p+A, 6213 p-p; off both kinematic curves: 18 p+A-pool tracks excluded (fraction 1.335e-03; 142 including primaries of p-p events); raw conflict tracks 505 (2.564e-02; crossover 1.78 deg), of which 10 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 3.33e+02, fraction < 1e-3 0.9733.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 5

Case 5 (row 19: graphite, 30 MeV, 0.2 cm, 20000000 primaries, 20000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 19654435 primaries removed (fraction 0.98272); 20000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 340490 | 0.01702 +- 2.9e-05 | 20.77 | 0.0776 | 0.0471 | 0.0310 | 0.546 |

p-p events: 0 with a secondary proton; 0 excluded (off the p+H curve or not exactly two proton records); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 345565 deflected primaries -> 345291 p+A, 274 p-p; off both kinematic curves: 5051 p+A-pool tracks excluded (fraction 1.462e-02; 5051 including primaries of p-p events); raw conflict tracks 1622 (4.694e-03; crossover 1.88 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 6

Case 6 (row 20: graphite, 50 MeV, 0.2 cm, 20000000 primaries, 20000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 19752729 primaries removed (fraction 0.98764); 20000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 245121 | 0.01226 +- 2.5e-05 | 16.09 | 0.1903 | 0.0806 | 0.0097 | 0.514 |

p-p events: 0 with a secondary proton; 0 excluded (off the p+H curve or not exactly two proton records); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 247271 deflected primaries -> 247154 p+A, 117 p-p; off both kinematic curves: 2131 p+A-pool tracks excluded (fraction 8.618e-03; 2131 including primaries of p-p events); raw conflict tracks 2126 (8.598e-03; crossover 1.87 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 7

Case 7 (row 21: graphite, 100 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3900305 primaries removed (fraction 0.97508); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 98446 | 0.02461 +- 7.8e-05 | 11.55 | 0.4123 | 0.2627 | 0.0045 | 0.559 |

p-p events: 0 with a secondary proton; 0 excluded (off the p+H curve or not exactly two proton records); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 99695 deflected primaries -> 99586 p+A, 109 p-p; off both kinematic curves: 1231 p+A-pool tracks excluded (fraction 1.235e-02; 1231 including primaries of p-p events); raw conflict tracks 1640 (1.645e-02; crossover 1.85 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 8

Case 8 (row 22: graphite, 150 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3950527 primaries removed (fraction 0.98763); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 49206 | 0.0123 +- 5.5e-05 | 9.49 | 0.5262 | 0.3833 | 0.0031 | 0.589 |

p-p events: 0 with a secondary proton; 0 excluded (off the p+H curve or not exactly two proton records); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 49473 deflected primaries -> 49448 p+A, 25 p-p; off both kinematic curves: 263 p+A-pool tracks excluded (fraction 5.316e-03; 263 including primaries of p-p events); raw conflict tracks 1238 (2.502e-02; crossover 1.83 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.


## Case 9

Case 9 (row 23: graphite, 200 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3966450 primaries removed (fraction 0.99161); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 33435 | 0.008359 +- 4.6e-05 | 8.23 | 0.5992 | 0.4625 | 0.0023 | 0.608 |

p-p events: 0 with a secondary proton; 0 excluded (off the p+H curve or not exactly two proton records); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 33550 deflected primaries -> 33536 p+A, 14 p-p; off both kinematic curves: 112 p+A-pool tracks excluded (fraction 3.338e-03; 112 including primaries of p-p events); raw conflict tracks 1071 (3.192e-02; crossover 1.80 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
