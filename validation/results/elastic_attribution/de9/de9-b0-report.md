Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E9 B0: independent reconstruction of the pinned p-p level, 100-150 MeV

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de9_pp_level`

Input files (sha256):

- `.ionmc-cache/ionmc-data/objects/27bcafb89cf0444c53c6b9f3dd17618c62e2e6b3694f70d31c4f502399f103f7`  27bcafb89cf0444c53c6b9f3dd17618c62e2e6b3694f70d31c4f502399f103f7

<!-- body below is rendered from the Git projection of the JSON record -->

D-E9 B0: exploratory; leads identified after the D-E8 STOP; not a re-evaluation. Independent reconstruction of the pinned p-p level; ratios only.

Match (relative difference <= 0.001 at every node and midpoint of 100-150 MeV): **True** (max |rec/table - 1| = 3.61e-04).

Dip (interior minimum more than 1% below both endpoints): **True** (minimum at 130 MeV; -0.1830 against 100 MeV, -0.1505 against 150 MeV).

Reading (pre-declared): **evaluation feature at the LA150 top nodes (the independent reconstruction reproduces the step)**.

| node (MeV) | reconstruction / table | within 1e-3 |
|---|---|---|
| 100 | 1.000000 | True |
| 130 | 1.000000 | True |
| 140 | 1.000000 | True |
| 150 | 1.000000 | True |

Quadrature convergence (GL 1024 / Simpson 4001 / Simpson 20001 against GL 2048 at 100 / 125 / 150 MeV): max |relative difference| 2.08e-12; accepted (< 1e-06): **True**.

Carry-up check above 150 MeV (rule: sigma(E) = sigma_rec(150 MeV) x S(E)/S(150 MeV), S = the table's own bgg.pp_bgg_mb): all nodes within 0.001: **True** (max |carried/table - 1| = 8.32e-13).

| node (MeV) | carried / table | within |
|---|---|---|
| 200 | 1.000000 | True |
| 250 | 1.000000 | True |

Level relative to the 125 MeV node:

| node (MeV) | reconstruction | table |
|---|---|---|
| 100 | 1.1905 | - |
| 125 | 1.0000 | - |
| 130 | 0.9727 | - |
| 140 | 0.9872 | - |
| 150 | 1.1450 | - |
