Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 9: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-937006e321a3e4822984-d1a0df23/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-937006e321a3e4822984-d1a0df23/work/ps.header`  4169c9ace57971eba1e874ae08bb8f4eb590133b3d36912532cf1b8521481c0c
- `.ionmc-cache/reference-runs/REF-937006e321a3e4822984-d1a0df23/work/ps.phsp`  e333d68e3e801250d79428ee65067e4b91f57c423c1ee19189784faca8dc021f
- `.ionmc-cache/reference-runs/REF-937006e321a3e4822984-d1a0df23/work/input.txt`  0ce90ec75465a77ebd7ebfcdda101929446147463051c4c85c0273f3e68a119e

<!-- body below is rendered from the Git projection of the JSON record -->

Case 9 (row 23: graphite, 200 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3966450 primaries removed (fraction 0.99161); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 33435 | 0.008359 +- 4.6e-05 | 8.23 | 0.5992 | 0.4625 | 0.0023 | 0.608 |

p-p events: 0 with a secondary proton; 0 excluded (off the p+H curve or not exactly two proton records); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 33550 deflected primaries -> 33536 p+A, 14 p-p; off both kinematic curves: 112 p+A-pool tracks excluded (fraction 3.338e-03; 112 including primaries of p-p events); raw conflict tracks 1071 (3.192e-02; crossover 1.80 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
