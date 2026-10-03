# RNG measurement scripts (decision 0037)

Scripts used on 2026-10-03 (Warp 1.17.0, numpy 2.5.3, CPU device) to measure
the properties of Warp's built-in generator and to check the Philox4x32-10
construction. Run them with a Python environment containing `warp-lang` and
`numpy` (for example `uv run python validation/scripts/rng/<script>`).

| Script | Purpose | Result obtained on 2026-10-03 |
|---|---|---|
| `rng_overlap.py N L` | Marks every 32-bit state visited by `wp.rand_init(seed, i)` followed by `L` draws for `N` histories in a bitmap and counts draws that revisit a state already used by another history | N=1e5, L=1000: 1.1 % of draws; N=1e6, L=1000: 10.7 %; N=1e6, L=2000: 20.0 % |
| `rng_seed_dupes.py` | Counts histories whose first draws are bit-identical between seeds `s` and `s+1` | 260 of 1e6 histories identical for seed pairs (1,2), (42,43), (1000,1001) |
| `philox_lib.py` | Philox4x32-10 as a precision-generic `@wp.func` factory (`make_philox(real)`) | — |
| `philox.py` | Checks the kernel implementation against the Random123 known-answer vectors and against a pure-Python integer implementation on 20 000 random counter/key blocks | all known-answer vectors match; 20 000/20 000 blocks bit-identical |

These scripts are measurement provenance for the decision, not part of the
`ionmc` package; they are archived unchanged. The package implementation is
`ionmc.rng.philox` with its tests in `tests/ionmc/test_rng.py` (task V3-003A).

Note (decision 0039): `philox_lib.py` maps a word to a uniform with the 24-bit form
`((w >> 8) + 0.5) * 2**-24` in both precisions. In float32 this form rounds to exactly 1.0
for `w >> 8 = 2**24 - 1`; the package uses `((w >> 9) + 0.5) * 2**-23` for float32 and keeps
the 24-bit form for float64 (regression test `test_u6c_u01_strictly_inside_unit_interval_in_both_precisions`).
