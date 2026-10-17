Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E6 lineage manifest of the swap tables

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de6_swap --build`

Input files (sha256):

- `validation/generated/elastic_attribution/de4/cdf/de4-case0-cdf.npz`  a393aa37d1c96439e8a42c99ab7575941f1e8269c5bff9f3044b525004c7d23f
- `validation/generated/elastic_attribution/de4/cdf/de4-case1-cdf.npz`  fb087211abd63b7c261a233ecd2f2b256e3c3dd070c73796ee39c152f705681f
- `validation/generated/elastic_attribution/de4/cdf/de4-case2-cdf.npz`  f1e0387efd51625f9c336f1bab6b0784e13d7a04a5d2a7a780ea81ecfc2118a0
- `validation/generated/elastic_attribution/de4/cdf/de4-case3-cdf.npz`  0f1a0c3d76e4e2e2605b1c8315d42348497d4aca59419360a59541d1523dd4d7
- `validation/generated/elastic_attribution/de4/cdf/de4-case4-cdf.npz`  a2ac4acba11d8cef9c27cb5a49115cc9d47945754824f3d5acc0c31d731b45b2

<!-- body below is rendered from the Git projection of the JSON record -->

table ids and npz hashes recomputed from the table files; CDF hashes from the D-E4 CDF files; bound by de6_swap --de6-summary

| variant | table id | npz sha256 | classes | CDF-derived |
|---|---|---|---|---|
| all | 08a10dd42ec6c913e26673e3770e778e81330ac73f2756a4ea403ac629c2cb50 | 84c371e94faffd56cb7b494befdd258527ba4e6fdb2faac4d24a3ad82ac221fc | - | False |
| swap-all | 69d3581511aa78d01bec5f09fe1538668e7d2369a11aec4253410e3ce4d03d2c | 177b6839af4397d0e657a148aee8507a5a8e1a9fb1cf9ac0e7c70c37e0bdf644 | O-16, H-1 | True |
| swap-both | 1082dfc68cf0f73ba523bae50c3cf8692c24e3dde726abb480ca15f98c2e97d6 | 09e37ca2c934444ab816ddd189935ca9591e8a8d9ccdca9d437e991839fc34ad | O-16, H-1 | True |
| swap-pA | 89090af08a53541c0fdb4294052237de6ca1739dafad2f41ba98ab9227d47cdc | 2cbe0feb6a533fd82920a844813a507efbb6e70c1dc223b8bab8703ba28ec0b5 | O-16 | True |
| swap-pp | 0366f9d116c0cc335bfe3e3157c72b4be04f59baf8eee86ac464e90eec09335f | 5f3981ec6576019e9abaade11a9b3859332f55811e7f51b6ab69125a533e9b41 | H-1 | True |
| swap-pp-sigma | cbb4495f9c88d24c200580966a535bc6b331da7ed806bdfe468f92697cf3be64 | 9d31435af82cf44e583d759a6d205e15665be62e1f36a79f9384445ccdc2add3 | - | True |

| case | row | E (MeV) | primaries | CDF npz sha256 | events O-16 / H-1 (stored) |
|---|---|---|---|---|---|
| 0 | 14 | 30 | 20000000 | a393aa37d1c96439e8a42c99ab7575941f1e8269c5bff9f3044b525004c7d23f | 140322 / 48206 |
| 1 | 15 | 50 | 20000000 | fb087211abd63b7c261a233ecd2f2b256e3c3dd070c73796ee39c152f705681f | 108606 / 26702 |
| 2 | 16 | 100 | 4000000 | f1e0387efd51625f9c336f1bab6b0784e13d7a04a5d2a7a780ea81ecfc2118a0 | 46790 / 14156 |
| 3 | 17 | 150 | 4000000 | 0f1a0c3d76e4e2e2605b1c8315d42348497d4aca59419360a59541d1523dd4d7 | 20974 / 12056 |
| 4 | 18 | 200 | 4000000 | a2ac4acba11d8cef9c27cb5a49115cc9d47945754824f3d5acc0c31d731b45b2 | 13464 / 11838 |
