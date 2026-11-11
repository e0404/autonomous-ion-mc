Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 8: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-92208aebedb8ce3aa598-b665cc82/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-92208aebedb8ce3aa598-b665cc82/work/ps.header`  3618858e040439b5fb40aff158a209b5e2aba0c674f37a9e828737952359bb00
- `.ionmc-cache/reference-runs/REF-92208aebedb8ce3aa598-b665cc82/work/ps.phsp`  085ff8e4ec1bda425fe1db47142ce5af3df0a3a4bc367e0a2e865b471619b021
- `.ionmc-cache/reference-runs/REF-92208aebedb8ce3aa598-b665cc82/work/input.txt`  b6c07c31049061140d0ce134b6af4068106d98ae9576b432bae8a333357d5732

<!-- body below is rendered from the Git projection of the JSON record -->

Case 8 (row 22: graphite, 150 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3950527 primaries removed (fraction 0.98763); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 49206 | 0.0123 +- 5.5e-05 | 9.49 | 0.5262 | 0.3833 | 0.0031 | 0.589 |

p-p events: 0 with a secondary proton; 0 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 49473 deflected primaries -> 49448 p+A, 25 p-p; off both kinematic curves: 263 p+A-pool tracks excluded (fraction 5.316e-03; 263 including primaries of p-p events); raw conflict tracks 1238 (2.502e-02; crossover 1.83 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
