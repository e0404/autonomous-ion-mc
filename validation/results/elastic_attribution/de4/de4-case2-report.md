Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 2: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-bf3f98356e973aede9ac-1a91e885/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-bf3f98356e973aede9ac-1a91e885/work/ps.header`  a59d90a44f67ee582b8c732a3b004f4c4fe75be45a9ecfc3a5103acb33f4feb2
- `.ionmc-cache/reference-runs/REF-bf3f98356e973aede9ac-1a91e885/work/ps.phsp`  7de07171974c2b188b4e849d4d6cdedd4cab73bbb8333d1a5837da99759dd446
- `.ionmc-cache/reference-runs/REF-bf3f98356e973aede9ac-1a91e885/work/input.txt`  21da81599245e6a6ea624e393045af00756c46e6ff33afa800c7c2d6e99dae0a

<!-- body below is rendered from the Git projection of the JSON record -->

Case 2 (row 16: water, 100 MeV, 1 cm, 4000000 primaries, 4007563 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3945499 primaries removed (fraction 0.98637); 4000000 primary records, 7563 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 7473); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 6805 | 0.001701 +- 2.1e-05 | 29.75 | 0.9311 | 0.8773 | 0.0000 | 25.969 |
| O-16 | 46790 | 0.0117 +- 5.4e-05 | 11.44 | 0.4044 | 0.2575 | 0.0047 | 0.418 |

p-p events: 7473 with a secondary proton; 415 excluded (off the p+H curve or not exactly two proton records); 253 removed by the 16.26 deg CM cut; 6805 kept (class H-1).

Classification: 54501 deflected primaries -> 47021 p+A, 7480 p-p; off both kinematic curves: 233 p+A-pool tracks excluded (fraction 4.954e-03; 451 including primaries of p-p events); raw conflict tracks 784 (1.439e-02; crossover 1.83 deg), of which 3 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 8.86e+02, fraction < 1e-3 0.9553.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
