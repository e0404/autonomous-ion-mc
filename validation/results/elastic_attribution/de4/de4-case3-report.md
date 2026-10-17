Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 3: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-bc42d09720c8d824c5cf-e71c94ad/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-bc42d09720c8d824c5cf-e71c94ad/work/ps.header`  0ffbfac845f553d4d4c46ace79920c6b56f8dca284115f363149ca910bebc720
- `.ionmc-cache/reference-runs/REF-bc42d09720c8d824c5cf-e71c94ad/work/ps.phsp`  76fdae7011c81a39e4e7270191325cfabe51d3d8cc985ff575526985cca854ee
- `.ionmc-cache/reference-runs/REF-bc42d09720c8d824c5cf-e71c94ad/work/input.txt`  8aa87f5db1a77ae962b2a99f53417b9c175d684b6611138bfbbf807b0d93fcea

<!-- body below is rendered from the Git projection of the JSON record -->

Case 3 (row 17: water, 150 MeV, 1 cm, 4000000 primaries, 4006429 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 3972614 primaries removed (fraction 0.99315); 4000000 primary records, 6429 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 6370); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 6028 | 0.001507 +- 1.9e-05 | 29.96 | 0.9877 | 0.9580 | 0.0000 | 39.083 |
| O-16 | 20974 | 0.005243 +- 3.6e-05 | 9.30 | 0.5132 | 0.3676 | 0.0036 | 0.434 |

p-p events: 6370 with a secondary proton; 55 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 287 removed by the 16.26 deg CM cut; 6028 kept (class H-1).

Classification: 27386 deflected primaries -> 21015 p+A, 6371 p-p; off both kinematic curves: 39 p+A-pool tracks excluded (fraction 1.856e-03; 179 including primaries of p-p events); raw conflict tracks 570 (2.081e-02; crossover 1.81 deg), of which 6 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 8.27e+02, fraction < 1e-3 0.9631.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
