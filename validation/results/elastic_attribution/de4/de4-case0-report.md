Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 0: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-579ee8d7c1bd2869cf86-c8aa5c87/work .ionmc-cache/reference-runs/REF-14940501a696f1e68846-bd77bc8d/work .ionmc-cache/reference-runs/REF-80139f68b6ba7fc8c05c-752c737f/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-579ee8d7c1bd2869cf86-c8aa5c87/work/ps.header`  5aa2e0aaedf625c967eb1fb4fd9466596aff859f2186277c12de781a8dee35a4
- `.ionmc-cache/reference-runs/REF-579ee8d7c1bd2869cf86-c8aa5c87/work/ps.phsp`  f454032b4661b409a0a2b9589c99b8b54a9f01238d28ca1e2f4e3d6b0e297e5f
- `.ionmc-cache/reference-runs/REF-579ee8d7c1bd2869cf86-c8aa5c87/work/input.txt`  fa9bbef7bdf1feb3d712f5e4889527e5196d10b298c7fc5f80cadaf2bce81482
- `.ionmc-cache/reference-runs/REF-14940501a696f1e68846-bd77bc8d/work/ps.header`  94ccb040a25fbfb77d169abe024392a784dbf7614ca44dfc62fa190bc75d4efa
- `.ionmc-cache/reference-runs/REF-14940501a696f1e68846-bd77bc8d/work/ps.phsp`  f467281bdf44d19f6e86d9348ef525da1629dcac515f570fb2234b43265ac6a6
- `.ionmc-cache/reference-runs/REF-14940501a696f1e68846-bd77bc8d/work/input.txt`  4bbfc3e3b130679fcb37eff9bdd82991702689607ae5e2160b163b1f26ff0ea3
- `.ionmc-cache/reference-runs/REF-80139f68b6ba7fc8c05c-752c737f/work/ps.header`  f930f453597d3968ecd2cbff6af981003b755459a2f40787cf1c1c0a026131ec
- `.ionmc-cache/reference-runs/REF-80139f68b6ba7fc8c05c-752c737f/work/ps.phsp`  c77392b91a179366ec24328036fb3cb08997005dd52fc07db66187a538ba05d4
- `.ionmc-cache/reference-runs/REF-80139f68b6ba7fc8c05c-752c737f/work/input.txt`  ac0028336fb23f0d4294aaa5c041d84e10d21c3b4b4b3627c655bd7582a37e94

<!-- body below is rendered from the Git projection of the JSON record -->

Case 0 (row 14: water, 30 MeV, 0.2 cm, 20000000 primaries, 20025531 proton records, target O-16).

KE filter (E > E0 - 1e-4 MeV): 19833617 primaries removed (fraction 0.99168); 20000000 primary records, 25531 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 25325); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| H-1 | 23782 | 0.001189 +- 7.7e-06 | 30.44 | 0.2031 | 0.0000 | 0.0000 | 7.789 |
| O-16 | 140322 | 0.007016 +- 1.9e-05 | 22.52 | 0.1108 | 0.0722 | 0.0345 | 0.486 |

p-p events: 25325 with a secondary proton; 556 excluded (off the p+H curve or not exactly two proton records); 987 removed by the 16.26 deg CM cut; 23782 kept (class H-1).

Classification: 166383 deflected primaries -> 141063 p+A, 25320 p-p; off both kinematic curves: 731 p+A-pool tracks excluded (fraction 5.182e-03; 1071 including primaries of p-p events); raw conflict tracks 527 (3.167e-03; crossover 1.86 deg), of which 25 p-p by the secondary rule. Secondary protons on the p+H curve: max residual 2.39e+02, fraction < 1e-3 0.9790.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
