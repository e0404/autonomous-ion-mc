# V3-006 acceptance plan — self-describing persisted results, reports, configuration files, `ionmc run`, requested/effective capability contract

> **Status: freeze candidate, revised 2026-11-23 after Codex REVIEW-bd983b92 of the S0 commit 00b5e3de
> (changes_required; findings 11–15 of the pre-freeze review log, resolved in place). The revised S0 commit is
> re-reviewed; the freeze is the first S0 commit whose Codex review passes (planner, 2026-11-23, base `v3/develop`
> 95937553 = the V3-013 squash merge, PR #92).** Frozen together with decision 0045
> (`decisions/0045-persisted-results-and-capability-contract.md`) by that commit, after an independent
> architecture/scientific review, before any V3-006 code exists. Later changes are dated amendments; earlier text is
> kept. Not a clinical claim.

Markers: **FACT** (read at 95937553; file:line in decision 0045, Context), **ESTIMATE**, **TO VERIFY**. Row classes:
*regression*, *contract* (fail closed), *round-trip*, *determinism*, *self-consistency*, *workflow*. Tiers: CI
(GitHub, small, warp-cpu, single process), LV (logged local exact-SHA validation), HR (controlled host runner, CUDA).
**No row is physics evidence**; the evidence category of every row is self-consistency (PROTOCOL taxonomy).
Single-process directive: `cpu_workers = 1` in every LV/HR run; multiprocess checks are `multiprocess`-marked,
deferred and recorded. Orchestrator decisions O1–O10 (2026-11-23) are adopted as stated in decision 0045.

## A. Requirement mapping (ledger rows advanced; exact evidence each needs)

V3-006 may move a row to "implemented" only (code merged, evidence pending); never "evidenced" or "satisfied" without
conformant qualification.

| Ledger row | Suite | What V3-006 delivers | Evidence required at the exact integration head | Remaining after V3-006 |
|---|---|---|---|---|
| V2-OUTPUT | S-PROV-OUTPUT | library bundle writer/reader; every V2-OUTPUT clause mapped to a JSON pointer or array | OUT-0..OUT-10, PROV-1..PROV-8, REP-1..REP-4 PASS (CI + LV; OUT-2 and PROV-7 on HR) | influence-matrix metadata (V3-007); conformant qualification |
| V2-CAP | S-CAP-CONTRACT | declared == enforced sweep; reconciliation; honoured check; requested AND effective persisted | CAP-1a, CAP-1c, CAP-2..CAP-7 PASS (CI); CAP-1b PASS (LV, CUDA cases HR) | V3-005B, V3-007, V3-009 extend the matrix |
| V2-USER | S-WORKFLOW-USER | TOML config, `ionmc run/check/capabilities/verify/report`, documentation | CLI-1..CLI-6, CFG-1..CFG-7, DOC-1 PASS (advanced, not closed) | clean-install qualification incl. data acquisition, multi-beamlet, influence output |
| V1-MUST-019 Separation of transport from data source | S-CAP-CONTRACT | selection by registered id in configs; identity persisted (no new mechanism) | CFG-4, PROV-3 | qualification |
| V1-MUST-020 No large physics datasets in Git | S-PROV-DATA | no table values in bundles or committed fixtures | OUT-8, REG-4 | qualification |
| V1-MUST-025 Provenance | S-PROV-OUTPUT | `provenance.code/data/models` incl. the decision-0044 block and nuclear table id | PROV-1..PROV-4 | qualification |
| V1-MUST-026 Offline reuse | S-PROV-DATA, S-WORKFLOW-USER | CLI runs from the cache with no network | CLI-4 | qualification |
| V1-MUST-027 Explicit random control | S-NUM-UNCERTAINTY | seed in config, persisted, reproduction from the bundle | OUT-3 | coverage rows (V3-011) |
| V1-MUST-030 Statistical uncertainty estimation | S-NUM-UNCERTAINTY | uncertainty definition, validity masks, undefined = NaN persisted | OUT-7 | coverage (V3-011) |
| V1-MUST-037 Explicit physical units | S-PROV-OUTPUT | unit table and per-array units | OUT-5 | qualification |
| V1-MUST-038 Coordinate-system documentation | S-PROV-OUTPUT | grid convention persisted; edges reproducible; docs page | OUT-6, DOC-1 | qualification |

Not addressed and not claimed: V1-MUST-009/010 (V3-007), V2-ION nuclear/fragments (V3-009), warp nuclear (V3-005B),
V2-RELEASE, any physics row, any evidence record.

## B. Design decisions (frozen by decision 0045; summary)

1. **Result schema** `ionmc.result/1`: directory bundle `result.json` (strict JSON, authoritative metadata, file
   hashes) + `arrays.npz` (deterministic zip writer, per-array `content_sha256`) + `report.md` + `plots/*.svg`
   (derived, byte-reproducible). `schema_version` "1.0"; the reader accepts only known versions and rejects unknown
   fields at every level. `kind = "transport"` (`"influence"` reserved).
2. **Provenance fields**: `code` (runtime Git tree oid of the package; commit/dirty from a wheel build stamp or the
   read-only `git` fallback with `GIT_*` stripped, else `unknown`; O2; build stamp per decision 0045 §3),
   `environment` (versions and devices; no hostname/user/home), `data` (dataset ids, archive/member sha256, per-run
   role), `models` (stopping id, version, `spec_sha256`, `evidence_at_run` verbatim as at run time; other model names;
   nuclear table id), `run` (backend, device, precision, seed, RNG identity, histories, batches, workers, chunking,
   non-deterministic timings, `output_path`, `cache_dir`, `cache_dir_source`), `statistics`, `accounting`, `units`.
   Absolute paths only at the frozen path set P of decision 0045 §2 (O7/O8).
3. **Requested/effective**: `requested` = canonical JSON of the config; `effective` = `EffectiveConfig.summary()`
   byte-for-byte (O4) incl. the decision-0044 `stopping` block; `effective_extra` holds new information (materials,
   grid convention, profile axis, resolutions).
4. **Capability contract**: declared (`capabilities()`, derived from enforcing objects) = enforced (sweep); the frozen
   field-by-field reconciliation table of decision 0045 §4 (equal / derived / unused / requested-only); result-level
   honoured check; no negotiation or fallback.
5. **Configuration**: TOML (`tomllib`), strict, `ionmc.config/1`; `physics.stopping_model` required (absence →
   `ConfigError`); files by path + mandatory sha256, relative paths resolved against the config directory and recorded
   absolute (O8).
6. **CLI**: `ionmc run | check | capabilities | verify | report`; the frozen exception → exit-code table of decision
   0045 §7 (0 success, 1 dataset/cache, 2 usage/config incl. unknown model id, 3 invalid result, 4 unsupported
   combination, 5 backend unavailable, 6 result format, 70 defect; subclasses before bases).
7. **Batches** persisted by default; `persist_batches = false` recorded (O3). **Profile axis** = dominant axis of the
   source direction, recorded (O9). **Golden fixture**: one synthetic KB-scale bundle (O10).


8. **Input manifest**: immutable `InputManifest` from the loader, passed to `save`; API runs record `source = "api"`
   with `unavailable` file pointers; `save` refuses changed files. **Round trip**: canonical semantic comparison over
   the reconstructible subset; the rest is recorded `not-reconstructible` (decision 0045 §6).
9. **Untrusted bundles**: `load` validates ZIP directory, catalogue and NPY headers before allocation and streams
   decompression under configurable byte budgets with safe defaults (decision 0045 §1).

## C. Seeds and determinism rules (declared before any run)

- CI tests use fixed small seeds written in the tests (not evidence).
- LV/HR workflow runs use the new family **`out6`: seed = 20700004 + 100·r_index + k**, k ≤ 99, r_index ≤ 49.
  Extent 20700004–20705003. r_index 0 = CLI water box proton (python, float64); 1 = same (warp-cpu, float32);
  2 = same (warp-cuda, float32, HR); 3 = heterogeneous voxel slab proton with LET and a synthetic lookup tally
  (warp-cpu); 4 = carbon box (warp-cpu); 5 = nuclear proton (python; derived nuclear table in the cache);
  6 = reproduction-from-bundle runs (OUT-3); 7 = `icru90-water-v1` alpha box (warp-cpu; cache required);
  all others reserved. Consumed once drawn.
- **Collision check (2026-11-23, every tracked text file):** no decimal seed literal lies in 20660000–20760000
  (hexadecimal digests contain such digit runs and are not seeds). The nearest declared families are V3-008 rerun2
  (≤ 20632003), the V3-005 continuation base 20641004 and the next base 20651004 with the 21 000 envelope
  (≤ 20672003), V3-013 `topas13`/`ion13` (20610001–20618103), 2098xxxx/2099xxxx and 21000000; `out6` is disjoint
  from all of them. The V3-005 fresh sequence (plan v3-005 "Seeds", step 10 000, envelope `base + 0…20999`) is
  open-ended: the windows of its bases 20681004, 20691004 and 20701004 reach `out6`; those of 20671004 (≤ 20692003)
  and 20711004 do not. The three bases are excluded by the dated note of 2026-11-23 in plan v3-005 (shared seed-base
  ledger); V3-005B restates the exclusion in its seed guard when it declares its bases.
- Guard: `tests/ionmc/test_v3_006_seed_families.py` (pure python) asserts the extent, disjointness from every family
  of `test_v3_008_seed_families.py` (V3-008 and `V3_005_FAMILIES`) and of `test_v3_013_seed_families.py`, the
  exclusion of the three V3-005 bases above, and that this section's literals match its constants.
- Determinism: (a) `save` of the same `Result` twice gives byte-identical files; (b) two runs with identical config
  and seed on the same backend/precision/host give identical `content_sha256` for every array, byte-identical
  `report.md` and SVGs (reports print neither timings nor path values), identical `files` hashes, and a `result.json`
  equal apart from exactly `/run/timings` and the pointers of the path set P (decision 0045 §2); (c) reports and plots
  are a pure function of the bundle. No statistical acceptance row exists in V3-006.

## D. Ordered steps (per step: commit → Codex review → exact-SHA validation record → push → CI)

Operator diagnostic rules apply to every step (binding, 2026-10-12): one resource-intensive job globally (no sandbox
pytest while a host job runs; no concurrent heavy sub-agents); `IONMC_SINGLE_PROCESS=1` and thread limits
(OMP/OPENBLAS/MKL/NUMEXPR = 1); every pytest run is a logged run `python -X faulthandler -u -m pytest -v
--capture=tee-sys --tb=long` with `PYTHONFAULTHANDLER=1 PYTHONUNBUFFERED=1`, complete unfiltered output in
`telemetry/test-logs/<UTC>-<label>/{command.txt, meta.json, log.txt}` outside the repository; after a crash the rules
for re-runs and narrowing apply. The worktree stays clean while a Codex review runs.

| Step | Content | Rows that must pass at the step's commit |
|---|---|---|
| S0 | decision 0045, this plan, ledger notes, V3-013 integration notes; no code | review only |
| S1 | `ionmc.configio`: `to_dict`/`from_dict` of every config dataclass, TOML loader, `InputManifest`, `ConfigError` | CFG-1..CFG-7, REG-1..REG-3 |
| S2 | `ionmc.results`: deterministic npz writer, `save`/`load`, schema 1.0, allowlist, invariants, resource budgets and adversarial fixtures; seed guard test | OUT-0, OUT-1, OUT-2 (CI part), OUT-4..OUT-10, PROV-8, seed guard |
| S3 | `ionmc.provenance`: tree oid, hatchling build-stamp hook (wheel only; `.gitignore` guard), read-only git fallback, environment, path set P, data/model aggregation; `environment.py` docstring amended | PROV-1..PROV-6 |
| S4 | contract: derived `capabilities()` entries, sweep, targeted rule pairs, `reconcile`, honoured check | CAP-1a, CAP-1c, CAP-2..CAP-7 |
| S5 | `ionmc.svg`, `ionmc.report`, `regenerate_report`, synthetic golden fixture | REP-1..REP-4 |
| S6 | CLI subcommands and exit codes | CLI-1..CLI-3, CLI-6 (CI part) |
| S7 | docs (results format, configuration reference, CLI, workflow), ledger rows | DOC-1, docs build |
| S8 | LV/HR runs (`out6`), one at a time: OUT-2/OUT-3 on three backends, CLI-4 offline, CLI-5, CLI-6, PROV-7 and REG-1 on the host, CAP-1b (LV, CUDA cases HR) | all rows |

### Discriminating rows (defined before implementation; all gating)

**OUT — persisted results**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| OUT-0 | Clause map: a table in the test module maps every V2-OUTPUT clause (incl. `Result.valid`, every `GridResult.valid`, backend, device, RNG identity) to a JSON pointer or array name; each exists and is non-empty or carries an explicit status | every clause resolved | CI |
| OUT-1 | In-process round trip `load(save(r))`: (a) proton box python; (b) voxel slab with LET, spectrum and synthetic lookup tallies, warp-cpu float64; (c) invalid result with `allow_invalid_result`; (d) a TOML run with its `InputManifest` (paths and hashes at `/input_manifest/*`); (e) the same config built through the API without a manifest (`source = "api"`, every file pointer `unavailable`) | arrays bitwise equal (dtype, shape, NaN positions); JSON values equal | CI |
| OUT-2 | **Clean-process round trip**: bundle written by one process, read by `python -I` (package path injected explicitly) which recomputes every digest and prints a canonical dump | dump equals the writer's | CI; LV/HR on all three backends |
| OUT-3 | Reproduction from the bundle: `from_dict(requested)` + same seed re-run (python, warp-cpu float64) | every `content_sha256` and `effective` equal | CI (small), LV (`out6` r_index 6) |
| OUT-4 | **Unknown-field fail-closed** and other negative controls, each separate: unknown key at top level and nested; missing required key; `schema_version` "1.1" and "2.0"; flipped byte in `arrays.npz`; altered member with recomputed file hash; NaN literal; array shape ≠ grid shape; extra npz member; missing file hash entry | each raises `ResultFormatError` naming the defect | CI |
| OUT-5 | Units: every catalogue entry has `units` in the vocabulary; every JSON quantity carries its unit in its key or a `units` sibling | exact | CI |
| OUT-6 | Coordinates: edges from origin/spacing/shape equal `ScoringGrid.edges_mm` bitwise; convention string present; `material_index`/density equal the input | exact | CI |
| OUT-7 | Undefined values: `defined_mask` False ⇒ NaN, never 0; JSON undefined scalars carry `status`; a writer patched to write 0 there is rejected by `load` | exact | CI |
| OUT-8 | **No table values**: npz names ⊆ allowlist; no member equals bitwise any stopping/range/LET-water/nuclear table array of the run; injecting `tables.range_mm` is refused by `save` | exact | CI |
| OUT-9 | Determinism: `save` twice → identical bytes; a different zlib level changes the file hash but not the content digests (documents which identity is authoritative) | exact | CI |
| OUT-10 | **Untrusted-bundle budgets** (synthetic KB-scale adversarial fixtures): NPY header longer than the bound; header disagreeing with the catalogue; duplicate ZIP member; more members than catalogue entries; declared shape whose byte size exceeds the budget (e.g. 2**40 elements); member whose declared uncompressed size exceeds the budget; compressed member that inflates beyond its declared size; oversized `result.json`; plus a positive control at the default budgets | each raises `ResultFormatError` before allocating the declared size (peak RSS increase bounded by 64 MiB, measured with `resource.getrusage`); the positive control loads | CI |

**PROV — provenance**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| PROV-1 | `package_tree_oid` of the checkout equals `git rev-parse HEAD:src/ionmc` on a clean tree (git run by the test) | equal | CI, LV |
| PROV-2 | A temporary package copy with one changed byte gives another oid and `git_dirty = true`; with no stamp and no git: `git_commit = null`, `status = "unknown"`; with `GIT_DIR`/`GIT_WORK_TREE` set to a foreign repository the fallback still reports the package's own repository | exact | CI |
| PROV-3 | `provenance.models.stopping` equals `effective.stopping` (`model_id_effective`, `version`, `spec_sha256`, `evidence_at_run`) for `bethe-ionmc-v1` and, cache present, `nist-pstar-astar-water-2005` and `icru90-water-v1` (two records: `consistency-only`, `topas-agreement-shared-lineage:V-C1:PASS`); a bundle written before a hypothetical registry append keeps its records on `load` | exact | CI (default), LV (tabulated) |
| PROV-4 | `provenance.data` = union of `data_used` pins, the nuclear table id (nuclear runs) and lookup sha256; a removed pin is rejected by `load` | exact | CI, LV |
| PROV-5 | Build stamp: `uv build --wheel` from a clean checkout in the LV environment (not the host runner) gives a wheel whose `ionmc/_build_info.json` carries the commit and a tree oid equal to the runtime oid of the installed wheel; the checkout stays clean (`git status --porcelain` empty, no stamp under `src/ionmc`); a build without a reachable repository stamps `commit = null`, `status = "unknown"`; `[build-system] requires` and `uv.lock` are unchanged by S3 | exact | LV (TO VERIFY: offline build environment from the uv cache) |
| PROV-6 | Privacy: absolute paths only at the pointers of the path set P (decision 0045 §2); outside P no string contains `Path.home()`, and no path component (split at `/` and `\`) and no whole non-path string equals `getpass.getuser()` or `os.uname().nodename` (no substring matching; a user named `root` does not hit `mcs_mapping.evaluation = "root"`, which is not a path); `cache_dir_source` recorded; reports print no P value | exact | CI |
| PROV-7 | Host runner: an `out6` r_index 2 bundle written under the controlled host runner (read-only snapshot at `/workspace`, `GIT_WORK_TREE` set, `HOME=/tmp/home`) at a clean committed SHA has `code.source = "git"`, `git_dirty = false`, `git_commit` = the run's code SHA and a tree oid equal to `git rev-parse <sha>:src/ionmc`; any other outcome is a FAIL handled by prospective amendment | exact | HR |
| PROV-8 | Input manifest: a TOML run records `source = "toml"`, `base_dir`, the config file and every referenced file with path and sha256; an API run records `source = "api"` and `unavailable` file pointers; editing a referenced file between load and `save` makes `save` raise `ConfigError` (CLI exit 2, no bundle); a manifest `content_sha256` that disagrees with the config arrays is refused | exact | CI |

**CAP — capability contract**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| CAP-1a | Cache-independent sweep: species × energy-domain class {below, inside, upper edge, above} × backend × precision × nuclear × tallies × trace, for `bethe-ionmc-v1` and its LS variants against `validate`; tabulated models only through their declared domains (no cache) | declared == enforced over these dimensions: supported ⇔ `validate` passes, unsupported ⇔ raises before transport; zero third outcomes | CI (CUDA cases as `BackendUnavailableError` when absent) |
| CAP-1b | Cache-backed sweep of `nist-pstar-astar-water-2005` and `icru90-water-v1` (species × energy-domain class × backend × precision × nuclear × {water, PMMA}) against `validate` | as CAP-1a; a missing cache object is a FAIL, never a skip | LV; HR (CUDA cases) |
| CAP-1c | Targeted rule pairs: one positive and one negative case for every other `validate` rule (geometry kind and voxel constraints, scoring-grid count and shape, source sigmas, cutoffs, end-of-range and MCS modes, `mcs_mapping_eval`, worker and chunk limits, memory budget, diagnostics incl. trace); a checklist of the rules (`validate` source order) is part of the test module and reviewed | each positive passes, each negative raises before transport | CI |
| CAP-2 | Declaration is derived: monkeypatching an enforcing constant changes `capabilities()` identically | exact | CI |
| CAP-3 | Reconciliation: (a) every row of the frozen table of decision 0045 §4 item 3 is exercised and recorded with its class and rule; (b) `reconcile` over the CFG-1 generated config set yields zero differences outside the table for every config that passes `validate`, and every effective pointer of the table resolves; (c) a synthetic unlisted difference (patched effective field) raises `UnsupportedCombinationError`; (d) the API behaviour of decision 0044 is preserved: neither or both of `stopping`/`stopping_model` still fail at `PhysicsOptions` construction | exact | CI |
| CAP-4 | Discriminating dispatch: for tallies, nuclear, lookups, multiple scattering, straggling and stopping model, patching the consumer to drop the feature makes `Simulation(...)` or `save` raise | each raises | CI |
| CAP-5 | Honoured check: a `Result` lacking one requested tally, or carrying nuclear tallies on a non-nuclear request, is refused by `save` | exact | CI |
| CAP-6 | Requested and effective persisted separately (`end_of_range_model = "default"`, LS `eff` on a proton, `max_steps` unset) | exact | CI |
| CAP-7 | **Unsupported combination refuses with a stable exit code**: the V3-005A C1 list and the V3-013 DOM/ROLE rows through `ionmc check` and `ionmc run` give the code of the frozen table of decision 0045 §7 (domain/role/combination 4; DOM-6 missing or mismatched cache object and parse failure 1; DOM-7 unknown id 2) with no bundle; warp + nuclear exits 4 until V3-005B | exact | CI |

**CFG — configuration files**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| CFG-1 | Round trip over the reconstructible subset (decision 0045 §6): for generated configs covering every field of every config dataclass except user `stopping` objects and in-memory arrays, `to_dict(from_dict(to_dict(c))) == to_dict(c)` with arrays compared by `content_sha256`, and `validate(...).summary()` equal | exact | CI |
| CFG-2 | Documented TOML examples give the same `validate(...).summary()` as the equivalent API constructions | equal | CI |
| CFG-3 | Unknown key, misspelt table, wrong type, missing `nuclear`, missing `physics.stopping_model` (error names that key), missing sha256, sha256 mismatch: `ConfigError` with the key path; `..` paths resolve and are recorded absolute with sha256 | exact | CI |
| CFG-4 | Stopping by `stopping_model` id only; a `stopping` table is rejected; unknown id raises | exact | CI |
| CFG-5 | Built-in material names resolve to the identical `Material`; inline definitions checked as by the API | exact | CI |
| CFG-6 | No environment variable overrides a physics/run field (`IONMC_*` set in the test) | exact | CI |
| CFG-7 | Not-reconstructible records: a config with a user `stopping` object and one with an API-created lookup table serialise to `{"status": "not-reconstructible", "reason": ...}` at the affected pointer; `from_dict` of either raises `ConfigError` naming it; no default is substituted | exact | CI |

**REP — reports and plots**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| REP-1 | **Byte-identical regeneration**: `regenerate_report(load(bundle))` in a clean process equals the files written at run time | bytes equal | CI, LV |
| REP-2 | Golden synthetic fixture (KB-scale, synthetic arrays, no table values, no stopping data) and its committed report/SVGs | bytes equal | CI |
| REP-3 | Report content: tree oid, commit or "unknown", validity banner, requested/effective table, every `evidence_at_run` record verbatim, data identities, statistics definition, energy balance, counters, missing items, "Not a clinical claim"; the regex `\bvalidated\b` never matches the generated text (report minus the verbatim evidence records), and a registry guard asserts it matches no evidence record text in `STOPPING_EVIDENCE` | exact | CI |
| REP-4 | SVGs parse with `xml.etree`; no NaN coordinates; undefined voxels are gaps; profile axis = dominant source-direction axis, recorded | exact | CI |

**CLI and workflow**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| CLI-1 | `ionmc run example.toml` → bundle; `ionmc verify` exit 0; `ionmc report --check` exit 0 | exact | CI |
| CLI-2 | Every line of the exception → exit-code table of decision 0045 §7 is provoked separately (incl. `ConfigError` vs its base `UnsupportedCombinationError`, the `__cause__` chain to a cache error, `TransportLimitError` 3 without a bundle and `allow_invalid_result` 3 with a `valid = false` bundle); no partial bundle on any failure | exact | CI |
| CLI-3 | `ionmc capabilities` JSON equals `capabilities()` | exact | CI |
| CLI-4 | Offline: `out6` runs with network denied on the host runner (default model) and from the populated cache (tabulated, nuclear) | complete; bundles verify | HR, LV |
| CLI-5 | `out6` r_index 0–2: bundles verify and differ only in backend/device/precision fields and array contents | exact | LV, HR |
| CLI-6 | r_index 4 (carbon), 5 (nuclear proton) and 7 (`icru90-water-v1`) carry their stopping/nuclear identities | exact | LV |
| DOC-1 | Results-format reference checked against the schema constants; configuration reference lists every key; docs build passes | exact | CI |

**REG — no regression and repository hygiene**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| REG-1 | A16 digests and the T1 trace fixture (asserted identity) | bit-identical | CI, LV |
| REG-2 | `EffectiveConfig.summary()` of a fixed config set equals the values at 95937553 | byte-identical JSON | CI |
| REG-3 | `Result` fields and `capabilities()` keys present at 95937553 still present (additions allowed) | exact | CI |
| REG-4 | **No table values in Git**: no file added by V3-006 contains a physics table value (fixtures synthetic) | exact | CI |
| REG-5 | **Single process only**: every V3-006 test and `ionmc run` of the LV/HR set runs with `cpu_workers = 1` under `IONMC_SINGLE_PROCESS=1`; no test spawns worker processes | exact | CI, LV |

## E. Independent review items

1. Codex review of every step's exact commit (lifecycle; a Codex approval is not scientific validation).
2. Independent architecture review (reviewer agent, read-only, not the author) of S0 and of the final head: (a) OUT-0
   completeness against the V2-OUTPUT text; (b) reader strictness: attempt a bundle that loads but misrepresents
   (swapped grids, renamed arrays); (c) CAP-1a/1b/1c completeness against `capabilities()` and the `validate` rules;
   (d) licence: no node-derived ICRU 90/NIST/ENDF values in bundles or fixtures; (e) privacy of provenance fields
   (O7); (f) no report text implies validation or clinical use, and evidence records stay verbatim; (g) the
   `environment.py` amendment and `GIT_*` handling (O2).
3. Contrary evidence is preserved: a row failing after the freeze is amended prospectively, never deleted.

## F. Integration gates

- **STOP**: any REG difference; any change of a stopping spec, `spec_sha256`, table hash or evidence record.
- **Integration** requires every row of section D PASS at the exact head (CI rows in GitHub CI and in a logged local
  run; LV/HR rows with exact-SHA records), S7 documentation, ledger rows set to "implemented" citing the rows (no
  "evidenced"/"satisfied"), decision 0045 "Later validation outcome" filled, and the standard lifecycle gates: passing
  independent Codex review of the exact head and current integration base, exact-SHA local validation record, push,
  PR, exact-head CI green, controlled merge, worktree retired with its branch preserved.
- Deferred, recorded, non-gating under the single-process directive: bundles from `cpu_workers > 1`.

## G. Acquisition list

None. Tabulated-model and nuclear LV rows use datasets already in the provenance cache; no download.

## Pre-freeze review log (2026-11-23; findings recorded, resolutions applied in place before the freeze)

Independent pre-freeze review (reviewer agent, read-only) of decision 0045 and this plan; findings 1–10:

1. **Blocking.** The requested → effective allowlist missed resolutions that `validate` makes (unit direction,
   `mcs_mapping_eval`, ion cutoff per u × A, unused nuclear/straggling/MCS/chunk options, added deuteron species and
   water row, LS-variant id on z = 1) → frozen field-by-field table with the classes equal / derived / unused /
   requested-only (decision 0045 §4 item 3; Context); row CAP-3 requires zero unlisted differences over the CFG-1 set.
2. Numbering of the V3-005 seed amendment → handled by the orchestrator (dated note of 2026-11-23 in plan v3-005);
   the §C refers to it and obliges V3-005B to restate the exclusion in its seed guard.
3. **Important.** O7 vs O8: absolute input and config paths contradicted the privacy rule; the default cache is not a
   user choice; a substring scan for the user name false-positives → frozen path set P of JSON pointers, scan at path
   components only, `cache_dir_source` recorded (decision 0045 §2); row PROV-6.
4. **Important.** Exit codes ambiguous (`ConfigError` subclasses `UnsupportedCombinationError`; `select()` raises
   unknown id and missing cache alike; invalid-result code unstated) → loader pre-checks and a frozen exception →
   exit-code table with subclass precedence (decision 0045 §7); rows CLI-2, CAP-7.
5. **Important.** Build stamp under-specified (sdist path without `.git`, host without hatchling, stamp not ignored)
   → `uv build --wheel` from a clean checkout in the local environment, stamp mapped into the wheel only, `.gitignore`
   guard, pre-declared `commit = null` outcome, no new build dependency asserted (decision 0045 §3); row PROV-5.
6. **Important.** Determinism (b) excluded too little → reports never print paths or timings; the excluded pointers
   are `run.timings` and the path set P (decision 0045 §2, §5); §C.
7. **Minor.** Context stale (derived species sections, `cli.py:99–100`, `Result` fields) → Context corrected;
   `GridResult.valid` added to OUT-0 clause map.
8. **Minor.** The `validated` regex collided with verbatim evidence records → regex scoped to generated text, plus a
   registry guard (decision 0045 §5); row REP-3.
9. **Minor.** Host-runner provenance unasserted → row PROV-7 (HR: `code.source = "git"`, clean tree).
10. **Minor.** "No integer literal" was literally false (hex digests) → "no decimal seed literal" with the verified
    window facts (§C).

Codex REVIEW-bd983b92 (REVIEW-bd983b92e86c4f14bdadde38fd9a2619, 2026-11-23) of the S0 commit 00b5e3de:
changes_required, five important findings, resolved in place:

11. **Important.** The reconciliation row "neither `stopping` nor `stopping_model` → `bethe-ionmc-v1`" contradicted
    decision 0044 and `select()` (exactly one is required; absence fails at `PhysicsOptions` construction) → row
    deleted; TOML configs require `physics.stopping_model` (absence → `ConfigError` naming the key); the API
    behaviour of decision 0044 is unchanged and preserved by CAP-3 (decision 0045 §4, §6; rows CFG-3, CAP-3).
12. **Important.** An exact config round trip over every field is impossible (user `stopping` objects, in-memory
    lookup and geometry arrays, `eq=False` classes) → canonical semantic comparison (`to_dict` equality, arrays by
    `content_sha256`) over the declared reconstructible subset; the rest round-trips to an explicit
    `not-reconstructible` record (decision 0045 §6; rows CFG-1, CFG-7).
13. **Important.** `save` could not produce the input-file pointers (`Result` and `LookupTable` keep no paths) →
    immutable `InputManifest` produced by the loader and passed to `save`; API runs record `source = "api"` with
    explicit `unavailable` file pointers; `save` refuses a manifest whose hashes no longer match (decision 0045 §2,
    §6; rows OUT-1 (d)–(e), PROV-8).
14. **Important.** CAP-1 over-claimed (dimensions omitted; no CI cache for tabulated models) → CAP-1a (CI,
    cache-independent matrix), CAP-1b (LV/HR, cache-backed tabulated models) and CAP-1c (one positive and one
    negative case for every other `validate` rule); claim restated as "declared == enforced over the enumerated
    dimensions"; the global "zero undeclared combinations" claim removed (decision 0045 §4 item 2).
15. **Important.** `load`/`ionmc verify` lacked resource-exhaustion defences for untrusted bundles → pre-allocation
    validation of the ZIP directory, catalogue and NPY headers, streamed decompression under configurable byte
    budgets with safe defaults (decision 0045 §1); row OUT-10 with synthetic KB-scale adversarial fixtures in S2.

## Status

- 2026-11-23 (planner, S0): freeze candidate; no code. Base 95937553 (V3-013 integrated). Freeze note:
  `<placeholder: S0 commit SHA, Codex review id, architecture review>`. Not a clinical claim.
- 2026-11-23 (planner): revised after Codex REVIEW-bd983b92 of 00b5e3de (findings 11–15, resolved in place); the
  revised S0 commit is re-reviewed; the freeze is the first S0 commit whose Codex review passes. Not a clinical claim.
