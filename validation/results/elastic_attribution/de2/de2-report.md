Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E2: p-p audit

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de2_pp`

Input files (sha256):

- `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.npz`  1028bb67f57bc96edafe6d3dbcce532013edbbd31b9b80861239b09d9ad4202a
- `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.json`  ab5e2a357ab53060a84182607d93fe4373e50ce60622450af74728fd4e495e77
- `.ionmc-cache/ionmc-data/objects/41268d804bd1fc3f165c88a36f055df5f65cfad3c7acdccdf06c3dd442b73dfd`  41268d804bd1fc3f165c88a36f055df5f65cfad3c7acdccdf06c3dd442b73dfd

Table sigma = Hale NI above 16.26 deg CM (above 150 MeV anchor x BGG ratio); BGG formula = 1.0115 (23 + 50 sqrt(ln(0.73/p)^7)) mb = total elastic, a different quantity. sigma_ET = sigma <T_slow>. Angles: the continuing (faster) proton. No download was made; PDG = the cached rpp2022 file.


| E | sigma table | builder formula | BGG | table/BGG | <T_slow> MeV | <T_slow>/E | sigma_ET (mb MeV) | th_fast lab p10/p50/p90 | th_fast CM p50 | P(th>th_2mm) | P(th>th_3mm) | PDG within 15 % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 140.12 | 140.12 | 157.36 | 0.890 | 5.44 | 0.272 | 762 | 17.0 / 31.4 / 42.3 | 63.1 | 0.000 | 0.000 | none |
| 50 | 49.77 | 49.77 | 52.53 | 0.947 | 13.32 | 0.266 | 663 | 16.3 / 30.8 / 41.9 | 62.2 | 0.701 | 0.526 | none |
| 100 | 27.30 | 27.30 | 27.61 | 0.989 | 26.01 | 0.260 | 710 | 15.1 / 29.9 / 41.5 | 61.1 | 0.937 | 0.882 | none |
| 150 | 26.26 | 26.26 | 23.86 | 1.101 | 38.71 | 0.258 | 1017 | 14.2 / 29.6 / 41.1 | 61.0 | 0.984 | 0.953 | none |
| 200 | 25.65 | 25.65 | 23.30 | 1.101 | 51.62 | 0.258 | 1324 | 14.0 / 29.3 / 40.8 | 61.0 | 1.000 | 0.988 | none |

## p-p events per primary along the water CSDA path

| E0 | P(pp) table | P(pp) table, E>50 MeV | P(pp) BGG formula (E >= 12.5) | transferred energy (MeV/primary) |
|---|---|---|---|---|
| 100 | 0.0249 | 0.0129 | 0.0263 | 0.349 |
| 150 | 0.0362 | 0.0259 | 0.0379 | 0.767 |
| 200 | 0.0519 | 0.0436 | 0.0517 | 1.558 |

PDG file: 154 points; lowest: 325 MeV (0.846 GeV/c) 23.2 mb; 418 MeV (0.979 GeV/c) 24.0 mb; 433 MeV (1.000 GeV/c) 25.8 mb; 440 MeV (1.010 GeV/c) 24.0 mb; 460 MeV (1.037 GeV/c) 22.0 mb
