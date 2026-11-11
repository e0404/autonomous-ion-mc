Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 7: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-87e1645504b0e38176a9-94553de7/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-87e1645504b0e38176a9-94553de7/work/ps.header`  8bebb312a63027f9ed5bcee251da56c789a91ae2ab7855fc08f9c8d51725e405
- `.ionmc-cache/reference-runs/REF-87e1645504b0e38176a9-94553de7/work/ps.phsp`  61e8fd2b89e69d03486d6125166ba822e3d236a4aa3baaa7787cfe6a2d1e8057
- `.ionmc-cache/reference-runs/REF-87e1645504b0e38176a9-94553de7/work/input.txt`  04e8ecc3f6ca8bfe923db85995b13b9725862b7175bf3e284ae2b13ab0a7cdc5

<!-- body below is rendered from the Git projection of the JSON record -->

Case 7 (row 21: graphite, 100 MeV, 1 cm, 4000000 primaries, 4000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 3900305 primaries removed (fraction 0.97508); 4000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 98446 | 0.02461 +- 7.8e-05 | 11.55 | 0.4123 | 0.2627 | 0.0045 | 0.559 |

p-p events: 0 with a secondary proton; 0 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 99695 deflected primaries -> 99586 p+A, 109 p-p; off both kinematic curves: 1231 p+A-pool tracks excluded (fraction 1.235e-02; 1231 including primaries of p-p events); raw conflict tracks 1640 (1.645e-02; crossover 1.85 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
