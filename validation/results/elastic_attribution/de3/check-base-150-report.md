Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E3 base check at 150 MeV: pinned table against the base copy

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de3_probes --check-base --reuse-runs --energy 150`

Input files (sha256):

- `validation/generated/elastic_attribution/de3/de3-150-pinned.json`  ba2c506e40f78f7fe24ab221c7dffbb5511b5e57f67aec5839316aee5dd7f6bb
- `validation/generated/elastic_attribution/de3/de3-150-all.json`  b42751ca39f13efe566bfe29ff85ca2f0940910ee98e7619f16b9f57bae11478

<!-- body below is rendered from the Git projection of the JSON record -->

Energy 150 MeV, n = 100000, 20 batches, seed 20990100; pinned table 1fd24cff, base copy 08a10dd4.

**Bit for bit: True** (depth-dose batches True, recoil channel True, counters True).

| array | sha256 pinned | sha256 base copy |
|---|---|---|
| idd_batches | 8df1969be0973e6f4b154a362845d232b46a0e5a99da00bf671e303d5b9e57da | 8df1969be0973e6f4b154a362845d232b46a0e5a99da00bf671e303d5b9e57da |
| recoil_idd_batches | f7585fb9c0da2fa3833efe8193bf81d02020da1effa1e980a6acac6b5ad50122 | f7585fb9c0da2fa3833efe8193bf81d02020da1effa1e980a6acac6b5ad50122 |
