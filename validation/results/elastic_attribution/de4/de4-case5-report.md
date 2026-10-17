Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 5: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-ed7bb5ebe8d6e5c36ea0-c34e534a/work .ionmc-cache/reference-runs/REF-180ebf0cf5858cb478c1-183bb0a5/work .ionmc-cache/reference-runs/REF-8b170b44753aa87093b9-341267cf/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-ed7bb5ebe8d6e5c36ea0-c34e534a/work/ps.header`  aad790ff65c74432b7f868e672a80bc8fdfea3e22870f55b8f31a87f9ec2a3e7
- `.ionmc-cache/reference-runs/REF-ed7bb5ebe8d6e5c36ea0-c34e534a/work/ps.phsp`  aa4787144142483a75069a8a8a9649caf73c87a5b6b3e0e7a1fbb30603b95e1f
- `.ionmc-cache/reference-runs/REF-ed7bb5ebe8d6e5c36ea0-c34e534a/work/input.txt`  b6b30dd0e917d6a0b3318ca6cc3386f5cf4e5504cdd8a50160a23bb37ba501bd
- `.ionmc-cache/reference-runs/REF-180ebf0cf5858cb478c1-183bb0a5/work/ps.header`  595d5c83d5e4a673891f72427515eccfde603a136f1e4323e4872c703f377329
- `.ionmc-cache/reference-runs/REF-180ebf0cf5858cb478c1-183bb0a5/work/ps.phsp`  6b8543712f18795f8f8826952772921f77cdf4521c845238d03be49db194478c
- `.ionmc-cache/reference-runs/REF-180ebf0cf5858cb478c1-183bb0a5/work/input.txt`  1f1328d2fd056e22da89e34f0038b3fb26b5f4edfe6816af507cd8e99205322f
- `.ionmc-cache/reference-runs/REF-8b170b44753aa87093b9-341267cf/work/ps.header`  dcdf53492dce22fec144a5a9202a5c020371279667d4ebf74d31de528cc5a39e
- `.ionmc-cache/reference-runs/REF-8b170b44753aa87093b9-341267cf/work/ps.phsp`  95470b664ef2680ef686e76e5e0dbaf6d9500d2a09c3c1813da14b359ccef781
- `.ionmc-cache/reference-runs/REF-8b170b44753aa87093b9-341267cf/work/input.txt`  391cc2b73aa32c96207106d268fa21cd174177e97effc8dc49a8e609ea3653f0

<!-- body below is rendered from the Git projection of the JSON record -->

Case 5 (row 19: graphite, 30 MeV, 0.2 cm, 20000000 primaries, 20000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 19654435 primaries removed (fraction 0.98272); 20000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event is accepted when its fast proton lies on the p+H curve. Selection rule: fast-proton-on-curve-v2.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 340490 | 0.01702 +- 2.9e-05 | 20.77 | 0.0776 | 0.0471 | 0.0310 | 0.546 |

p-p events: 0 with a secondary proton; 0 excluded (fast proton off the p+H curve: p-p followed by a second interaction); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 345565 deflected primaries -> 345291 p+A, 274 p-p; off both kinematic curves: 5051 p+A-pool tracks excluded (fraction 1.462e-02; 5051 including primaries of p-p events); raw conflict tracks 1622 (4.694e-03; crossover 1.88 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
