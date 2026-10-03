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
`ionmc` package. The package implementation of Philox and its tests are
delivered by the transport task (V3-003).
