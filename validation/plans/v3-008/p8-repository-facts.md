# V3-008 P8: repository facts (read-only record)

Worktree `/home/wahln/aiprojects/ion-mc-worktrees/v3-008`, HEAD `2632a46560af7f607e590896c7727406c9aeb61c`,
merge-base with v3/develop `c790df369c03d1d775c843d4f1e13ab9337137e5`. No design content; facts only.

## (a) A16 baseline constant

```
$ grep -n "A16_BASELINE\s*=" validation/scripts/transport/steps_v4.py
112:A16_BASELINE = "f3a1dd62ea2f57a3f4c07999935f317b3c044871"
$ git rev-parse f3a1dd62
f3a1dd62ea2f57a3f4c07999935f317b3c044871
$ git log --oneline -1 f3a1dd62
f3a1dd6 V3-003D: range-table consistency - exact log-log quadrature CSDA range and re-qualification (#80)
```

Result: the constant equals `f3a1dd62ea2f...` (the V3-003D merge commit), a full 40-character SHA. The test
`tests/ionmc/test_run_suite.py:974-979` (`test_a16_has_no_recorded_exception_and_the_baseline_is_the_v3_003d_merge`)
asserts the 40-hex form and `startswith("f3a1dd62")`. `step_a16` (steps_v4.py:1191-1278) resolves
`A16_BASELINE^{commit}` in git (line 1204) and fails closed if the commit is absent. `A16_INTENDED_CHANGE` in
`run_suite.py` is `None` (docstring near lines 436-458).

## (b) ASTAR file in the task cache and registry

Registry entry, `src/ionmc/data/registry.py` lines 82-95: id `nist-astar-water-2005`, version
`ASTAR 2.0.1 (2005), liquid water (matno 276)`, POST to `https://physics.nist.gov/cgi-bin/Star/apdata.pl`
(body `prog=ASTAR&matno=276&ShowDefault=on&NumofEnergies=0&...`), pinned
sha256 `dd4a2e504e19beca1f893ca1e2a7d75629bcabc04e002d8e62ef176f3b69b46f`, 8568 bytes, parser `nist_star`,
role `evaluation`. The registry entry declares no energy range (the `Dataset` fields carry none); the range is
that of the file.

Cache (`/home/wahln/.cache/ionmc-experiment/experiment-v3/ionmc-data`): manifest
`manifests/nist-astar-water-2005.json` (retrieved_at 2026-10-03T18:29:05+00:00) and object
`objects/dd4a2e50...b69b46f`.

```
$ sha256sum objects/dd4a2e504e19beca1f893ca1e2a7d75629bcabc04e002d8e62ef176f3b69b46f
dd4a2e504e19beca1f893ca1e2a7d75629bcabc04e002d8e62ef176f3b69b46f   (matches the pinned value)
$ tail -1 <that object>
1.000E+03 1.557E+01 4.848E-03 1.558E+01 3.812E+01 3.811E+01 0.9997
$ head (first data row)
1.000E-03 9.891E+01 2.282E+02 3.271E+02 3.273E-06 1.706E-06 0.5212
```

Result: ASTAR (alpha in liquid water, matno 276) spans 1.000E-03 to 1.000E+03 MeV total kinetic energy; last row
1000 MeV, total stopping power 15.58 MeV cm2/g, CSDA range 38.12 g/cm2. Columns: kinetic energy, electronic,
nuclear, total stopping power, CSDA range, projected range, detour factor. For comparison the cached PSTAR
object (`f5fc6482...`) ends at 1.000E+04 MeV. The loader is `src/ionmc/data/nist_star.py`
(`parse_star_text` line 57, `load_star_table` line 112); no loader test of the ASTAR last row was run for this record.

## (c) Species table and ids 5 and 6

`src/ionmc/species.py` (not `physics/`): `MAX_ION_SPECIES = 64`, `PSEUDO_BASE = 64`, `N_SPECIES_SLOTS = 65`
(lines 28-30). `_ION_ENTRIES` (lines 50-56) holds, append-only:

| id | projectile |
|---|---|
| 0 | PROTON |
| 1 | DEUTERON |
| 2 | TRITON |
| 3 | HELIUM3 |
| 4 | ALPHA |

Pseudo-species: id 64 `nuclear_local`. Ids 5 and 6 are free on this branch (ion block 0-63; nothing at 5-63).
The unmerged V3-005B branch (PR #82) is NOT inspected from this worktree: whether it claims ids 5 or 6 is
"to be checked at integration".

## (d) Implementation entry points (proton-only gates)

```
$ grep -n "projectile != PROTON" src/ionmc/config.py
653:    if src.projectile != PROTON:
$ grep -n "if projectile.z != 1" src/ionmc/transport/tables.py
200:        if projectile.z != 1:
```

- `src/ionmc/config.py:653-657`, inside `validate()`: `if src.projectile != PROTON: raise fail(f"projectile
  {src.projectile!r} is not supported: the transport engine implements exactly the canonical proton {PROTON!r}
  (ions arrive with V3-008)")`. This rejects every non-proton source, including deuteron and alpha.
- `src/ionmc/transport/tables.py:200-201`, in `TransportTables.from_stopping_tables`: `if projectile.z != 1:
  raise ValueError("transport tables are implemented for z = 1 projectiles only")`. A second guard at lines
  204-205: `if projectile == PROTON and tables[0].projectile.z != 1: raise ValueError("... protons only (z = 1 ions)")`.
