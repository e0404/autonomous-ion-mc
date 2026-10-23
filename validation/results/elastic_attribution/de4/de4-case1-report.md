Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 1: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-2f91c432adc6f10f1d00-3bf76819/work .ionmc-cache/reference-runs/REF-fff91ef2d473e092ab93-86e6ad20/work .ionmc-cache/reference-runs/REF-cda480cdcddf8e1acc50-8e677ef4/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-2f91c432adc6f10f1d00-3bf76819/work/ps.header`  c46ddd5ca7e71410c9b85e8745ea7f349fec4b7a5ed11921018410740cd0f33d
- `.ionmc-cache/reference-runs/REF-2f91c432adc6f10f1d00-3bf76819/work/ps.phsp`  108dc7ea71f137baca01e169fefbb0e20c65124f1b99593ec6a35c4433ef8dba
- `.ionmc-cache/reference-runs/REF-2f91c432adc6f10f1d00-3bf76819/work/input.txt`  6cc9291b96b3feea70face37a8008376d8d2b5a3a9ebcd04e97206834fb5193f
- `.ionmc-cache/reference-runs/REF-fff91ef2d473e092ab93-86e6ad20/work/ps.header`  30ea5322d7381e24fe3d1c2a6c0528247530371f9892b15c9dad6eef76d513e4
- `.ionmc-cache/reference-runs/REF-fff91ef2d473e092ab93-86e6ad20/work/ps.phsp`  edc0552dfec405e82751adc5bb7cc94b71d8503136659796a6b295c3df787bec
- `.ionmc-cache/reference-runs/REF-fff91ef2d473e092ab93-86e6ad20/work/input.txt`  b2ea0b4b4a9dd56480a28359c591fa3f2fed17f5b9957fc1ba6bf39952206828
- `.ionmc-cache/reference-runs/REF-cda480cdcddf8e1acc50-8e677ef4/work/ps.header`  045acc317a2177bb4d96a43f8cd1fd0ed271ffb36bd4b9a573e257587c8e6208
- `.ionmc-cache/reference-runs/REF-cda480cdcddf8e1acc50-8e677ef4/work/ps.phsp`  ad7c7c8858b3c40d065ff54dbfcb710558eeb253a4e2c600316d23fc93866aa5
- `.ionmc-cache/reference-runs/REF-cda480cdcddf8e1acc50-8e677ef4/work/input.txt`  b4a2006c43b2db9cc3485f8dbe376806e6e480793a9bf1936dfbc6c101974b45

<!-- body below is rendered from the Git projection of the JSON record -->

Case 1 (row 15: water, 50 MeV, 0.2 cm, 20000000 primaries, 20014038 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 19877083 primaries removed (fraction 0.99385); 20000000 primary records, 14038 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 13968); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 13351 | 0.0006676 +- 5.8e-06 | 30.44 | 0.6790 | 0.5140 | 0.0000 | 13.047 |
| O-16 | 108606 | 0.00543 +- 1.6e-05 | 16.81 | 0.2171 | 0.0988 | 0.0111 | 0.425 |

p-p events: 13968 with a secondary proton; 84 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 533 removed by the 16.26 deg CM cut; 13351 kept (class H-1).

Classification: 122917 deflected primaries -> 108937 p+A, 13980 p-p; off both kinematic curves: 337 p+A-pool tracks excluded (fraction 3.093e-03; 483 including primaries of p-p events); raw conflict tracks 865 (7.037e-03; crossover 1.85 deg), of which 9 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 2.33e+02, fraction < 1e-3 0.9848.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
