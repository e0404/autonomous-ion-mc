Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 4: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-8c7b3198f2262cf40909-cd63ec2a/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-8c7b3198f2262cf40909-cd63ec2a/work/ps.header`  3297d4f1818bea9a47b4593662d50e8b086eaa6645c9275c5ed2199c7d607f82
- `.ionmc-cache/reference-runs/REF-8c7b3198f2262cf40909-cd63ec2a/work/ps.phsp`  46e53fa578a98bfbb5795bba94183aeb7a90a6bf18195656af5990d26019d4ca
- `.ionmc-cache/reference-runs/REF-8c7b3198f2262cf40909-cd63ec2a/work/input.txt`  1b9b5c5995e5086575c1b7bba7baed0c3fb4df27ad675d9bf27dce0ffd1bfeac

<!-- body below is rendered from the Git projection of the JSON record -->

Case 4 (row 18: water, 200 MeV, 1 cm, 4000000 primaries, 4006265 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3980303 primaries removed (fraction 0.99508); 4000000 primary records, 6265 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 6214); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 5919 | 0.00148 +- 1.9e-05 | 28.70 | 1.0000 | 0.9856 | 0.0000 | 50.508 |
| O-16 | 13464 | 0.003366 +- 2.9e-05 | 7.93 | 0.5711 | 0.4340 | 0.0022 | 0.430 |

p-p events: 6214 with a secondary proton; 41 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 254 removed by the 16.26 deg CM cut; 5919 kept (class H-1).

Classification: 19697 deflected primaries -> 13484 p+A, 6213 p-p; off both kinematic curves: 18 p+A-pool tracks excluded (fraction 1.335e-03; 142 including primaries of p-p events); raw conflict tracks 505 (2.564e-02; crossover 1.78 deg), of which 10 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 3.33e+02, fraction < 1e-3 0.9733.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
