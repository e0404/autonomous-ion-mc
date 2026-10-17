Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E3 base check at 200 MeV: pinned table against the base copy

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de3_probes --check-base --reuse-runs --energy 200`

Input files (sha256):

- `validation/generated/elastic_attribution/de3/de3-200-pinned.json`  81cdc0c19ae5a52ab522046a00d7f865c872086887aedf5588a2a06f82f2c2ab
- `validation/generated/elastic_attribution/de3/de3-200-all.json`  81956f90ccdb89380bc26240be370d91afa11df9a889aa658e734e71fbdd1538

<!-- body below is rendered from the Git projection of the JSON record -->

Energy 200 MeV, n = 100000, 20 batches, seed 20990101; pinned table 1fd24cff, base copy 08a10dd4.

**Bit for bit: True** (depth-dose batches True, recoil channel True, counters True).

| array | sha256 pinned | sha256 base copy |
|---|---|---|
| idd_batches | 1647c1bc7d9df59cd7509cf397139c89e1afff3b833bd38968adeabe5aaf6705 | 1647c1bc7d9df59cd7509cf397139c89e1afff3b833bd38968adeabe5aaf6705 |
| recoil_idd_batches | 32224351a7a571e7078245ac908e9e169cc2ffe1c13328aee18526ad9cd36e2d | 32224351a7a571e7078245ac908e9e169cc2ffe1c13328aee18526ad9cd36e2d |
