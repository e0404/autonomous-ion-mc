Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 analysis, case 6: TOPAS phase space -> per-class summary

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --parts .ionmc-cache/reference-runs/REF-cbd81eb04a582c623dbf-78402c95/work .ionmc-cache/reference-runs/REF-57494cfafb768a80ce07-0466d071/work .ionmc-cache/reference-runs/REF-11dd54c6c3c5470b58cb-835b0b5c/work`

Input files (sha256):

- `.ionmc-cache/reference-runs/REF-cbd81eb04a582c623dbf-78402c95/work/ps.header`  32bb4eaed753bea24be221db7dc75efef86f2a1972fbc73d926f0311fc688126
- `.ionmc-cache/reference-runs/REF-cbd81eb04a582c623dbf-78402c95/work/ps.phsp`  334d1eb5ecda7532314f78c758e54afc5c5afa16f662859083bb814ee19f8583
- `.ionmc-cache/reference-runs/REF-cbd81eb04a582c623dbf-78402c95/work/input.txt`  093e075dc315f0d501d05dca44c73dfbb65307e97d400c6814c5a1e3285d6bd4
- `.ionmc-cache/reference-runs/REF-57494cfafb768a80ce07-0466d071/work/ps.header`  94a4ded33286f36c6a5f4d294e2993fb186b7b1a7e38d452e0f6886b0881691b
- `.ionmc-cache/reference-runs/REF-57494cfafb768a80ce07-0466d071/work/ps.phsp`  8701a82076b6c2ec3acc3e70e2d9d3ce0f036d1aac1c7dda0a5fb2796d992a93
- `.ionmc-cache/reference-runs/REF-57494cfafb768a80ce07-0466d071/work/input.txt`  9f7d6b3e9bb09423b2364bd7f97c974d92f3f9fbf8a3acf7be62153b0c451521
- `.ionmc-cache/reference-runs/REF-11dd54c6c3c5470b58cb-835b0b5c/work/ps.header`  1008cd450d729d81f10f5ee92de5e270c4af059a7c14da28704ef70ba558fa7f
- `.ionmc-cache/reference-runs/REF-11dd54c6c3c5470b58cb-835b0b5c/work/ps.phsp`  19cd67f442bc1f0a69bba21d1094bd724da08b2a7cb3c2ed80063f14f87b885c
- `.ionmc-cache/reference-runs/REF-11dd54c6c3c5470b58cb-835b0b5c/work/input.txt`  d2aa3604bd57f93ee33e402854989d6963e2bb79237be4640faaa1d991548c64

<!-- body below is rendered from the Git projection of the JSON record -->

Case 6 (row 20: graphite, 50 MeV, 0.2 cm, 20000000 primaries, 20000000 proton records, target C-12).

KE filter (E > E0 - 1e-4 MeV): 19752729 primaries removed (fraction 0.98764); 20000000 primary records, 0 secondary proton records (parent != 0). Secondary definition: event has a proton record with parent ID != 0 (cross-check: events with >1 proton record = 0); a p-p event needs exactly two proton records.

| class | events | rate/primary +- SE | theta_lab median | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV |
|---|---|---|---|---|---|---|---|
| C-12 | 245121 | 0.01226 +- 2.5e-05 | 16.09 | 0.1903 | 0.0806 | 0.0097 | 0.514 |

p-p events: 0 with a secondary proton; 0 excluded (off the p+H curve or not exactly two proton records); 0 removed by the 16.26 deg CM cut; 0 kept (class H-1).

Classification: 247271 deflected primaries -> 247154 p+A, 117 p-p; off both kinematic curves: 2131 p+A-pool tracks excluded (fraction 8.618e-03; 2131 including primaries of p-p events); raw conflict tracks 2126 (8.598e-03; crossover 1.87 deg), of which 0 p-p by the secondary rule. Secondary protons on the p+H curve: max residual nan, fraction < 1e-3 nan.

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted.
