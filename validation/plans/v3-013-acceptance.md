# V3-013 acceptance plan — explicit selection of stopping-power models and tabulated sources

> **Status: DRAFT freeze candidate (S0, planner, 2026-11-17, base `v3/develop` a9563412; revised 2026-11-17 after
> the pre-freeze reviews, section H).** Frozen together with
> decision 0044 (`decisions/0044-stopping-source-selection.md`) by one reviewed commit (independent Codex review and an
> independent scientific reviewer) before any V3-013 code, transport run, TOPAS run or comparison output exists. Later
> changes are dated amendments; earlier text is kept. No outcome changes the default or any registered spec within this
> task. Not a clinical claim.

Markers: **FACT** (read in the repository at a9563412), **ESTIMATE**, **TO VERIFY**. Row classes: *regression*,
*self-consistency*, *capability* (fail closed), *consistency* (input-table reproduction; never evidence),
*shared-lineage MC agreement*, *report-only*. Every repository fact quoted here is recorded with file and line in
decision 0044, Context (Q4–Q8).

## A. Rows and pre-declared tolerances

### A.1 Regression and identity of the default (gating; any difference = STOP)

| Row | What | Tolerance | Class | Where |
|---|---|---|---|---|
| R2 | A16 nuclear-off digests (`a16_digest.py`, all specs used by plan v3-008 R2) and the T1 proton trace fixture (`test_transport_diagnostic.py`, asserted table identity) against `v3/develop` a9563412; `A16_INTENDED_CHANGE` stays `None` | bit-identical | regression | CI + LV |
| R1-ion | Repeatability C 290 MeV/u (tallies, counters, grids, end states) on python and warp-cpu, default model selected through `stopping_model="bethe-ionmc-v1"` and through `stopping=BetheStoppingSource()` | bit-identical between the two selections and between repeats | regression | CI + LV |
| ID-DEF | `TransportTables.identity` and `sha256` for proton, alpha, carbon12, oxygen16 in water and in `materials.PMMA`, default model: identity dicts equal to the values computed at a9563412 (portable); `sha256` equal on the same machine (libm-dependent, decision 0044 Q4) | identical | regression | CI (identity), LV (sha256, same host) |
| ID-NIST | As ID-DEF for `nist-pstar-astar-water-2005` (proton, alpha; water) against the `"nist-star"` source at a9563412, cached tables required (skip-with-record in CI if the cache is absent; run in LV) | identical | regression | LV |
| ID-GUARD | `DEFAULT_STOPPING_MODEL_ID == "bethe-ionmc-v1"`, its `spec_sha256` pinned; `STOPPING_DEFAULT_INTENDED_CHANGE is None`; a synthetic change of the default without a record naming a decision fails (mirrors the A16 intended-change tests of `test_run_suite.py`) | exact | regression | CI |
| ID-DIAG | The sources of `parity.py`, `ion_mcs_checks.py`, `nuclear/build.py` D6 and `a16_digest.py` equal `resolve("bethe-ionmc-v1")` (dataclass equality) | exact | regression | CI |
| ID-USER | `stopping=BetheStoppingSource()` resolves to `bethe-ionmc-v1`; `BetheStoppingSource(e_min_per_u=1.5)`, `BetheStoppingSource(e_max_per_u=20)` and a source with `barkas=False` resolve to `model_id_effective = "unregistered"` with their complete options in the identity block, never to a registered id (decision 0044 §2) | exact | regression | CI |

### A.2 Consumer consistency (gating; synthetic test-only model, no reference data)

A test-only model `test-bethe-x1.05` (registered in the test module only; tabulated kind; S = 1.05 × the default's S on
the default grid, same start-range rule scaled consistently) is selected by `stopping_model` in a small water geometry
(python and warp-cpu backends, float64, fixed seed from §C r_index 10). Every length parameter of the test-model run
(voxel spacing, box length, source start position, `z_clip`, `max_step_mm`, the residual-range fraction in mm, the
cutoff residual length) is scaled
by 1/1.05 so that both runs take the same step sequence in range units (scientific review I1).

| Row | What | Tolerance | Class |
|---|---|---|---|
| CC-R | CSDA range at 10 energies (proton, alpha) = default range / 1.05 | ≤ 1e-12 relative | self-consistency |
| CC-T | Mean depth of the cutoff points (proton 100 MeV, alpha 100 MeV/u; MCS and straggling off; 2 000 histories; lengths scaled as above) = default mean depth / 1.05, same seed, on python and on warp-cpu | ≤ 1e-9 relative (deterministic straight paths with identical step sequences in range units; a larger residual is a code defect to be attributed, not a tolerance to relax) | self-consistency |
| CC-LET | LET water row (decision 0040; `_water_table`) of the run = 1.05 × the default's row at every node | ≤ 1e-12 relative | self-consistency |
| CC-RES | Unit level: the `continuous_residual_v1` residual-range segment (`_residual_range` / `_deposit_residual` of `transport/reference.py`) at fixed energies (E_cut and four energies below it; alpha and proton) = default / 1.05 | ≤ 1e-12 relative | self-consistency |
| CC-ID | Identity block names `test-bethe-x1.05`; every table identity of the run (materials, water row) names the same model; `TransportTables.sha256` differs from the default's | exact | self-consistency |
| CC-NEG | Monkeypatching one consumer at a time (material tables, water row, deuteron tables, exact-Ω ion rows) to the default source is detected by the single-model identity assertion of `validate` | each patch raises | capability |
| CC-HASH | New-model metadata carry `model_id`/`version` (hash differs); pre-existing models' metadata byte-identical; the (metadata, projectile) → effective-model-id map is injective over the registry (z = 1 LS variants resolve to the default) | exact | regression |

### A.3 Domains and roles (gating; every case raises `UnsupportedCombinationError` naming model, species, material, energy and domain)

| Row | Case |
|---|---|
| DOM-1 | `icru90-water-v1` with carbon12, oxygen16 or deuteron |
| DOM-2 | `icru90-water-v1` or `nist-pstar-astar-water-2005` with any non-water geometry material (structural rule, not the name) |
| DOM-3 | `icru90-water-v1` alpha with `E0 + 6σ_E` above 250 MeV/u; proton above 500 MeV |
| DOM-4 | `icru90-water-v1` with requested water I ≠ 78 eV |
| DOM-5 | `nuclear=True` with `icru90-water-v1` (no deuteron table) and with `nist-pstar-astar-water-2005` (existing rule) |
| DOM-6 | Missing cache object, archive or member hash mismatch, parse failure |
| DOM-7 | Unknown `model_id`; both or neither of `stopping`/`stopping_model` given; a `stopping` object without `spec`/`identity()` |
| ROLE-1 | A spec whose `DataPin.use_role` is `evaluation-only` (synthetic pin shaped like libdEdx ICRU73_NEW) or `lineage-only` (shaped like a G4EMLOW file) cannot be resolved |
| ROLE-2 | `ExactOmegaMapping.build`/`fdm_domain_status` with a missing ion table raise (no Bethe default in run paths) |
| ROLE-3 | z = 1 with `bethe-ionmc-v1+ls-eff`: `model_id_requested = "bethe-ionmc-v1+ls-eff"`, `model_id_effective = "bethe-ionmc-v1"`, tables bit-identical to the default; for z ≥ 2 requested = effective |

### A.4 Persisted identity (gating)

| Row | What | Tolerance |
|---|---|---|
| PER-1 | No library writer exists (decision 0044 Q5): `EffectiveConfig.summary()` is written with the JSON writer of `steps_v8.py` (or its V3-013 successor), read back in a **clean process** (`python -I`, fresh interpreter); the reloaded `stopping` block resolves (`resolve(model_id_effective)`) to a source whose `spec_sha256`, data pins, `content_sha256` and table identities equal the block | exact |
| PER-2 | Ion runs record `mcs_proton_path_stopping = "bethe-ionmc-v1"`; runs with `icru90-water-v1` record `circular_rows` (the ICRU 90 part of I2 and the decision-0038 ICRU 90 comparisons) | exact |
| PER-3 | Per-run data pins (Codex REVIEW-88ef9d57 finding 3): a proton run of `nist-pstar-astar-water-2005` has `data_used` == [the PSTAR pin, `role_in_this_run = "construction"`, table `("proton", water)`] and `data_declared` == [PSTAR, ASTAR]; an alpha run has `data_used` == [the ASTAR pin, construction, `("alpha", water)`]; a run of `icru90-water-v1` with a proton lists the single ICRU 90 pin with table `("proton", water)`; negative control: a block whose `data_used` lists a pin of a projectile not transported is rejected by the block validator used by PER-1 | exact (S3 for the NIST cases, S4 for ICRU 90) |

### A.5 `icru90-water-v1` consistency checks (gating for implementation correctness; **not evidence**)

| Row | What | Tolerance |
|---|---|---|
| CC-NODE | The transport table's S at every fine-grid node (200/decade) equals the log-log interpolant of the parsed ICRU 90 array evaluated there (proton, alpha); the chord error of the fine table at the ICRU nodes themselves (off the fine grid) is reported as an aggregate only (max, RMS) | ≤ 1e-12 relative at the fine nodes; ICRU-node error report-only |
| CC-INT | Interpolation error between nodes: log-log vs a monotone cubic in log-log space (both from the array only), max and RMS per decade | report-only (aggregates in Git) |
| CC-I2 | On the common fine grid (`compare_nist.py` uses `log_grid(1, min(500, last node))`, the same nodes), the node arrays R(E_i) − R(E_0) (cumulative exact log-log increments, E_0 = 1 MeV/u) of `icru90-water-v1` and of the I2 "integrated ICRU 90" construction (`_icru_range_table`) are compared node by node (alpha and proton); the start-range offset (own-array integral vs the I2 NIST start) is reported; off-node `range_at` values are not compared (the two constructions differ by their start offset inside the ln R interpolation); at the first node both arrays are exactly 0 (Codex REVIEW-8e40febe minor 2) | exact equality at the first node (both 0.0), ≤ 1e-12 relative at every later node; offset report-only |
| CC-SPAN | Loader test: counts and coverage of both grids (57 / 49 nodes; window inside the span); no value printed | exact |

### A.6 Configuration-specific validation of `icru90-water-v1` (graded once, after the seal of §D)

**V-C1 — TOPAS agreement, helium (gating for the model's evidence record; shared-lineage transport agreement,
never a validation of the stopping physics).**
- Reference: TOPAS 4.3 / Geant4 11.4.2, `g4em-standard_opt4` only, the TD-EM configuration of plan v3-008 (water,
  1 km electron cut, beam and scoring as TD-EM, IDD bins per the sizing table below), 4 parts × 2.5e4 primaries per case.
- Cases (fresh, held-out; decision 0044 Q6): **He 75, 125 and 200 MeV/u** (cases 1, 2, 3 of §C). None was run in any
  V3-008 row; all lie inside the He domain, away from the domain ends whose behaviour was not re-measured after
  R2-STOP3, and away from He 220 MeV/u (declared for the superseded, never-run V1-held row).
- Sizing by the §3 rules of r2-topas-target.md (scientific review I3): R_plan = CSDA range of the default model in
  water (`BetheStoppingSource()`, rounded to 0.1 mm, computed before any run); L = 1.3 R_plan rounded up to 50 mm;
  h_z = 0.02 / 0.05 / 0.1 mm for R < 30 / 30–120 / ≥ 120 mm; entrance window [1, 3] mm (R ≥ 30 mm); plateau group
  g = 1 / 2 / 5 mm for R < 30 / 30–120 / ≥ 120 mm; IDD scorer only (no lateral scorers).

| case | He MeV/u | total MeV | R_plan (mm) | box L (mm) | h_z (mm) | entrance window (mm) | plateau group g (mm) |
|---|---|---|---|---|---|---|---|
| 1 | 75 | 300 | 46.7 | 100 | 0.05 | [1, 3] | 2 |
| 2 | 125 | 500 | 115.7 | 200 | 0.05 | [1, 3] | 2 |
| 3 | 200 | 800 | 261.9 | 350 | 0.1 | [1, 3] | 5 |

- ionmc: `stopping_model = "icru90-water-v1"`, warp-cpu float64, 1e6 histories, 40 batches, default end-of-range model
  for ions (`continuous_residual_v1`), all other options default.
- Criteria and estimators: TD-EM (a)–(e) with the tolerances, the estimators (r2-topas-target.md §4.1), the TOST rule
  (§5: PASS iff |x − r| + 1.645 σ_comb ≤ tol; inconclusive = FAIL) and the run-validity rules (§6), as in force at
  a9563412 including every R2-STOP3 amendment; evaluated with the V3-008 evaluator functions where their interfaces
  allow, otherwise with a sealed V3-013 evaluator that imports them.
- Lineage statement (binding wording of any result): TOPAS He stopping above 7.9452 MeV is the same ICRU 90 alpha
  array (decision 0043 §6, FC-2); V-C1 is evidence for transport, straggling and end-of-range consistency of this
  configuration with the reference, **not** for stopping physics. Residual construction differences: table resampling
  and interpolation (ionmc log-log on 200/decade; Geant4 its own dE/dx table binning), start range below 1 MeV/u.
- Outcome per case: PASS / FAIL. The evidence registry of decision 0044 (separate from the immutable spec) receives
  one record: `topas-agreement-shared-lineage:V-C1:PASS` only if all three cases PASS, otherwise `failed:V-C1`;
  never `validated` (Codex REVIEW-0f80fad8 finding 2; scientific review I5). The record is appended once and never
  re-graded; the spec, its hash and the table hashes are unchanged either way.

**V-C1-A — the default at the same energies (report-only).** The paired default-model samples of V-C2 graded with the
same criteria. They never change a V3-008 verdict, the default or any spec. Its results are committed in full and
recorded in the ledger whatever the outcome.

**V-C2 — selection reaches transport (gating for integration, §F).**
- Prediction, computed and committed in the seal before any run (model-only quantities, no reference value):
  ΔR_pred(E0) = R_CSDA,icru90(E0) − R_CSDA,default(E0) in mm of water for the three V-C1 energies (Gaussian-averaged
  over σ_E if the TD-EM source has σ_E > 0).
- Measurement: ΔR80 = R80(icru90) − R80(default) from paired ionmc samples with **identical seeds** (common random
  numbers, declared in §C), same estimator as TD-EM (a); σ by delete-one-batch jackknife over the 40 paired batches.
- PASS per energy iff |ΔR80 − ΔR_pred| + 1.645 σ ≤ max(0.02 mm, 0.10 |ΔR_pred|). Justification (ESTIMATE): the
  model-dependence of the R80 − CSDA offset through the Bohr-variance weighting is ≲ 0.005 mm; detour and binning
  terms cancel in the pair. The tolerance does not depend on any TOPAS output.
- Informativeness precondition (scientific review I2), evaluated and committed in the seal before any run: at each
  energy |ΔR_pred| ≥ 2 · (0.02 mm + 1.645 σ_exp), with σ_exp = √2 · SE_1e6 and SE_1e6 = 0.01 mm · √(1e5/1e6), the
  R80 standard error of the V3-008 TD-EM He budget (r2-topas-target.md §6) scaled to 1e6 histories, with no credit
  for the common random numbers (ESTIMATE; the precondition is therefore reproducible from these literals). Where
  the precondition holds, V-C2 is gating at that energy. Where it fails, V-C2 at that energy is non-informative, is
  still run and reported, and the gate for that energy is V-C2-D alone.

**V-C2-D — selection reaches transport, deterministic (gating for integration at every energy).**
- For the three energies, MCS and straggling off, python and warp-cpu float64, 2 000 histories per model, ionmc
  seeds of §C r_index 15. Per model m the telescoped end point Z_m = mean over histories of
  (z_end + R_m(E_end)), where z_end is the cutoff depth and E_end ≤ E_cut the energy at which the history ended
  (the last step's overshoot below E_cut and the residual below E_cut are thereby removed; second-pass finding N1).
  The prediction ΔR_pred,D(E0) = R_CSDA,icru90(E0) − R_CSDA,default(E0) at E0 exactly (no σ_E) is computed and
  committed in the seal.
- PASS per energy and backend iff (i) for each model |Z_m − R_m(E0)| ≤ 1e-4 · R_m(E0) (round trip of the forward
  and inverse tables along the step sequence; residual ≲ 1e-5 relative, scientific review I1; ESTIMATE) and
  (ii) |(Z_icru90 − Z_default) − ΔR_pred,D| ≤ 1e-4 · R_CSDA,default(E0).
- Informativeness: |ΔR_pred,D| ≥ 5e-4 · R_CSDA,default(E0) at every energy is checked in the seal; if it fails at an
  energy, the seal records it and the declared substitute gate there is the synthetic-model row CC-T (both
  backends) together with CC-ID and PER-1 evaluated for `icru90-water-v1` (the identity block and every table
  identity of the run name `icru90-water-v1`), so that the selected model's own path to transport is still tested.

**V-C3 — table spread (report-only).** `icru90-water-v1` vs ASTAR (alpha) and PSTAR (proton) stopping and CSDA range
over the common window, aggregates only (max, RMS; decision 0038 already reported this spread, SEEN).

## B. What is NOT claimed; measured-data gaps

- No claim that any model is TOPAS-equivalent; V-C1 is a shared-lineage agreement claim for one configuration.
- No physical validation of ICRU 90 stopping by V-C1, CC-NODE or CC-I2; no claim for protons with `icru90-water-v1`
  beyond consistency (no proton row is pre-registered; a later amendment may add one with fresh seeds).
- No C/O alternative; the default's C/O low-energy limitation and the I3 FAIL stand.
- Measured He/C/O ranges, measured lateral data and measured C/O stopping at therapy energies: evidence gaps
  (experimental validation deferred by the operator directive of 2026-11-07); not relabelled.
- No evidence transfers between models. Not a clinical claim.

## C. Seeds (declared before any run)

- Rule source: plan v3-008 §C and `tests/ionmc/test_v3_008_seed_families.py` (reservation
  `range(20572004, 20632004)`). Used inside it: rerun3 20572004–20578103, lv8 20581004–20587103, hr8 20591004–20597103,
  rehearsal8 20600000–20609999, topas8 20619001–20620940 (incl. the consumed shadow range 20619013–20620521), rerun1
  20621004–20627103, i3fu2 20627213–20628721, rerun2 20631004–20632003. Free gap used here: **20610000–20618999**.
- **TOPAS family `topas13`: seed = 20610000 + 100·part + case**, case 1 = He 75, 2 = He 125, 3 = He 200 MeV/u; parts
  0–3 = V-C1 evidence; parts 4–7 = frozen repeats (part + 4), drawn only on a run-validity repeat trigger of
  r2-topas-target §6; part 19 = execution smokes (2 000 primaries; never evidence); every other (part, case) with
  part ≤ 19, case ≤ 40 reserved. Extent 20610001–20611940.
- **ionmc family `ion13`: seed = 20612004 + 100·r_index + k**, k ≤ 99 (sample chunk/part index as in the V3-008 r2
  sample layout), r_index ≤ 60: r_index 0–2 = He 75/125/200 V-C1/V-C2 samples, used for **both** `icru90-water-v1`
  and the default (common random numbers, V-C2); r_index 10 = A.2 test seed; r_index 11–14 = LV R1-ion/consumer runs; r_index 15 = V-C2-D;
  others reserved. Extent 20612004–20618103.
- Consumed by V3-008 and never re-drawn: every family listed above. Consumed by V3-013 once drawn: `topas13` parts used,
  `ion13` r_index used.
- Ownership: the gap 20610000–20618999 is transferred from the V3-008 reservation to V3-013 by a dated note in plan
  v3-008 §C (same commit; scientific review M3); the V3-008 guard is unchanged.
- Guard: a new pure-python test (`tests/ionmc/test_v3_013_seed_families.py`) asserts both families lie inside the
  reservation, are mutually disjoint and disjoint from every V3-008 and V3-005 family of the existing guard, and that
  this section's literals match its constants.

## D. Order of execution and chronology

1. **Freeze (S0):** decision 0044, this plan, ledger rows, decision 0043 note — one commit; independent Codex review
   and scientific review; no code before it.
2. **S1** registry with the default (A.1 rows); **S2** single selection point, consumers, domains, identity
   (A.2–A.4); **S3** `nist-pstar-astar-water-2005` (ID-NIST); **S4** `icru90-water-v1` (A.5). Each: commit → Codex
   review on the exact commit until passed → exact-SHA local validation record → push → CI at the exact head.
3. **Seal:** one reviewed commit containing the V-C evaluator (`validation/scripts/stopping/v13_evaluate.py`, with
   tests on synthetic inputs), the three TOPAS bundles (parts 0–7, 19) with their hashes, the ionmc sample step
   definitions, the committed V-C2/V-C2-D prediction file with both informativeness preconditions evaluated, and
   the pins `spec_sha256`, `content_sha256` and table identities of both models plus the tree hash of `src/ionmc`
   at the seal. **Sealed closure (Codex REVIEW-c661e903 finding 2 and REVIEW-88ef9d57 findings 1–2; the V3-008
   I3-FU-1 mechanism):** the seal commit also carries a frozen manifest
   (`validation/reference_cases/topas/v13-frozen-manifest.json`) listing every sealed path: the evaluator, every
   analysis module it imports (its transitive closure inside the repository), `sealed_launch.py`, `pyproject.toml`,
   `uv.lock`, `infrastructure/keys/trusted_attestation_keys.json`, the environment-digest script
   `validation/scripts/stopping/env_digest.py` and the sealed environment record (below), every TOPAS bundle file
   (`input.txt`, `case.json`) and the ionmc sample step definitions, plus the `src/ionmc` tree hash. The manifest
   does not list itself (a file cannot contain its own digest; Codex REVIEW-0e4576e7 finding 1): it is bound by its
   blob at the signed seal SHA, and its sha256 is additionally written into the pin note below. **Binding to the
   reviewed seal SHA:** the seal commit's signed Codex review record (passed, `head_sha` = seal) is committed
   verbatim at `validation/results/stopping/v3-013/seal-review.json` by the next commit, and the plan carries a
   one-line dated pin note `seal review <id> at <sha>; manifest sha256 <digest>` in the form of the V3-008 I3-FU-1
   pins (that commit is itself reviewed and signed). **Sealed environment (finding 2):** `env_digest.py` (stdlib only)
   computes one digest per component — the interpreter executable, the stdlib tree (`sysconfig` `stdlib` and
   `platstdlib`, extension modules included), the complete installed distributions of `numpy` and `warp` (every
   file listed in each distribution's `RECORD` via `importlib.metadata`, plus every sibling native-library
   directory such as `numpy.libs` and every `*.libs` directory next to them; Codex REVIEW-fd8c722d) — as the
   sha256 over the sorted list of (relative path, file sha256), plus the `uv.lock` sha256; the seal file
   `validation/results/stopping/v3-013/seal-environment.json` holds these digests as recorded by a controlled
   host-run record (`run_host_validation` of `env_digest.py` on the execution environment, made before the seal
   commit and cited by id). **Warp generated code (REVIEW-fd8c722d):** every ionmc evidence run on a Warp backend
   starts with a fresh, empty kernel cache (the sample runner sets `WARP_CACHE_PATH` to a new empty directory
   inside the run's output directory before importing `warp`; a non-empty cache at start is a run-validity
   failure), and at its end records the sha256 of every generated shared object and source file in that cache
   (`warp_modules`) in its run record; the evaluator requires that every evidence run of the same backend lists
   identical `warp_modules` hashes and records them in the result, and marks a run invalid otherwise. The evaluator (a) verifies the seal review record
   against the sealed trusted-keys file (blob at the seal SHA) and the pinned seal SHA, (b) for every sealed path
   compares the on-disk bytes AND the blob at the pinned seal SHA (`git cat-file -p <seal>:<path>`, stdlib
   subprocess, read-only) with the manifest hash, and the on-disk manifest with its blob at the seal SHA and with
   the pin-note digest — the evaluator's own file and `sealed_launch.py` included — so that evaluator and
   manifest cannot be changed together after outcomes while keeping the old seal review, (c) verifies the
   `sealed_launch.py` launch record and recomputes the environment digests with the sealed `env_digest.py`,
   stopping on any difference from `seal-environment.json`; every ionmc evidence run computes the same digests at
   its start (the sample runner imports `env_digest.py`) and writes them into its run record, and the evaluator
   marks the run invalid on any difference, (d) for every evidence run (ionmc and TOPAS) requires that the blobs
   of every sealed path and the `src/ionmc` tree at the run's `code_sha` are identical to the seal's (so sample
   definitions, scoring, extraction, the validation scripts and the trust map cannot be changed at the execution
   SHA and restored later; equality of blobs, not ancestry), and (e) for every TOPAS run compares the executed
   `input.txt`/`case.json` hashes recorded in the host-run record with the sealed bundle hashes. Any mismatch marks
   the run invalid (TOPAS, ionmc) or stops the evaluation (seal files, environment, trust map).
   A change to any sealed file after the seal is a new seal (dated amendment, reviewed before any further run;
   seen samples consumed).
4. **ionmc samples first** (scientific review B1; the step-0 rule of r2-topas-target.md §6): the V-C2-D runs and the
   V-C1/V-C2 samples of both models (`icru90-water-v1` and the default, common random numbers) are produced at the
   seal commit; their records (sample digests, identity blocks) are committed and reviewed **before any TOPAS run,
   smokes included**. The evaluator fails closed unless every ionmc evidence run's identity block equals the seal
   pins exactly (equality, not ancestry), its `code_sha` has the seal as an ancestor, and the blobs of every sealed
   path and the `src/ionmc` tree at its `code_sha` are identical to the seal's (step 3 (d); second-pass finding N3
   and Codex REVIEW-88ef9d57 finding 2: neither transport code nor sample definitions, scoring or extraction can
   differ at the execution SHA).
5. **TOPAS smokes** (execution only, part 19), then **TOPAS evidence runs**, one heavy job at a time. A TOPAS run is
   valid only if the seal review and the ionmc-sample commit's review both finished before the run's controlled
   `reference_started` event (the V3-008 R2-TOPAS-4 rule; controlled records, not Git timestamps), and only if
   its executed input hashes equal the sealed bundle hashes (step 3). The evaluator enforces these checks and
   fails closed.
6. **Evaluate once** at a committed SHA; commit the reports and run records; dated results note in this plan; ledger
   update. A defect found after the seal is handled by a prospective dated amendment, reviewed before any re-run;
   seen samples are consumed (the I3-FU-1 revision-17 precedent).

## E. Compute and effort (ESTIMATES)

- S1–S4 implementation: ≈ 3–4 working days including review rounds; tests only (A rows: < 1 h LV on the host).
- TOPAS V-C1: 3 cases × 4 parts × 2.5e4 He primaries, EM-only: ≈ 0.5–1.5 h (scaled from the V3-008 timing smokes,
  TO VERIFY by the part-19 smokes).
- ionmc: 6 samples (3 energies × 2 models) × 1e6 histories, warp-cpu f64: ≈ 1.5–3 h.
- Evaluation: < 0.5 h. Total compute ≈ 2.5–5 h, single process.

## F. Gates and integration criterion

- **STOP**: any R2, R1-ion, ID-DEF or ID-GUARD difference; any A.2–A.4 failure that is not a test defect.
- **Integration of V3-013** requires: all A.1–A.4 rows PASS, A.5 gating rows PASS, V-C2-D PASS for all three energies on both backends and V-C2 PASS at every energy where its
  precondition holds,
  V-C1 graded (PASS or FAIL, never INCONCLUSIVE-by-omission: a run-validity failure after the frozen repeats is a FAIL),
  the documentation step S6 (default statement of decision 0044 §6, registry table, THIRD_PARTY_NOTICES update),
  the standard lifecycle gates (Codex review, exact-SHA validation record, CI, controlled merge).
- **No outcome changes the default or any spec within this task.** A FAIL is recorded and sets the model's evidence
  status; a revised alternative is a new `model_id`/`version` with fresh, held-out evidence pre-registered by amendment.
- A V-C2-D or V-C2 FAIL blocks integration (selection does not reach transport as predicted) and triggers attribution by
  prospective amendment.

## G. Acquisition list

- S1–S4 and V-C1/V-C2/V-C3: **none new** (ICRU 90 file afd14794…, PSTAR f5fc6482…, ASTAR dd4a2e50… are cached and
  pinned).
- Only if a later task opens the low-energy C/O composite (S5b): `G4WaterStopping.cc` at Geant4 v11.4.2 via
  `acquire_reference_data` (role exploratory until a decision), followed by a written licence and provenance review
  (Geant4 licence; data provenance stated in the file) before any construction role. Optional, for lineage only:
  `G4ICRU90StoppingData.hh` and `G4BraggModel.cc` v11.4.2 (Q7 TO VERIFY items).
- libdEdx ICRU73_NEW as construction: only after an operator authorization (GPL-3); not requested by this plan.

## H. Pre-freeze review log (2026-11-17)

Codex REVIEW-fd8c722d432349a5b2f2c2a9b4c90837 at cbfe16bb: changes_required, one important finding, applied:
numpy sibling native libraries and Warp's generated kernel cache were outside the digest → complete installed
distributions (RECORD plus `*.libs`), fresh Warp kernel cache per evidence run with every generated module
hashed into the run record and required identical across runs of a backend (§D step 3).
Codex REVIEW-0e4576e72e644515a7fffa05e9b21efd at 61ab457e: changes_required, two important findings, applied:
(1) the manifest could not list its own digest → the manifest is bound by its blob at the signed seal SHA and by
its sha256 in the reviewed pin note; (2) runtime and trust map not bound → sealed environment record
(`env_digest.py`: interpreter, stdlib, numpy and warp trees, `uv.lock`) recorded by a controlled host run before
the seal, recomputed by every ionmc evidence run and by the evaluator; the trusted-keys file is a sealed path
(§D step 3).
Codex REVIEW-88ef9d57a5324a7286c5439282864c1c at 338db616: changes_required, three important findings, applied:
(1) the manifest and closure were not bound to the reviewed seal SHA → committed seal review record, pin note,
blob-level verification at the seal SHA of every sealed path including the manifest, the evaluator and
`sealed_launch.py`, launch record and environment recorded (§D step 3); (2) ionmc runs were bound only by the
`src/ionmc` tree hash → blobs of every sealed path at the execution SHA must equal the seal's (§D steps 3–4);
(3) `data_declared`/`data_used` had no discriminating row → PER-3 with per-projectile exact assertions and a
negative control (A.4).
Codex REVIEW-c661e903895c44eb9b7d41fb7d95f404 at 625ddab0: changes_required, two important findings, applied
(decision 0044 review log): per-run data pins `data_declared`/`data_used` (0044 §4) and the sealed closure
with manifest, self-verification and executed-input hash checks (§D step 3, step 5).
Codex REVIEW-8e40febed0a04536ab61b1c354190a4d at 92cf8929: **pass** (signed), two minor findings: the V-C2-D end
point (same defect as N1 below; resolved by the telescoped end point) and the undefined relative error at the
first CC-I2 node (resolved: exact equality there). Second pass (verification of 92cf8929 by the same reviewer):
all first-pass findings and both first-pass Codex findings
RESOLVED; new findings N1 (important: ΔR_pred,D compared with a cutoff depth that is R(E0) − R(E_end) plus the
last-step overshoot, so the V-C2-D row could fail without a defect → telescoped end point Z_m = mean(z_end +
R_m(E_end)) with a per-model round-trip gate), N2 (r_index 11 assigned twice → V-C2-D uses r_index 15), N3 (pinned
`src/ionmc` tree hash never checked → evaluator checks it), N4 (σ_exp scaling ambiguous → literals fixed, no
common-random-number credit), N5 (the CC-T substitute did not exercise `icru90-water-v1` → CC-ID/PER-1 for
`icru90-water-v1` added to the substitute). Applied in place before the freeze.

Decision 0044, "Pre-freeze review log", records the Codex review REVIEW-0f80fad83827437582b3def7bf34b9b7 at
14099895 (changes_required: evidence status inside the immutable spec; `validated:V-C1` overstatement) and the
independent scientific review (B1 run order; I1 tolerances of CC-T/CC-RES/CC-I2/CC-NODE; I2 V-C2 informativeness;
I3 bin sizes; I4 unregistered sources; I5 status label; I6 decision-0038 amendment; M1–M6) with their resolutions.
In this plan the resolutions are: row ID-USER (A.1); the length-scaled, two-backend CC-T, the unit-level CC-RES and
the (metadata, projectile) injectivity (A.2); the fine-node CC-NODE and the node-array CC-I2 (A.5); the per-case
sizing table, the evidence-record wording, V-C1-A full commitment, the V-C2 precondition and the deterministic
V-C2-D (A.6); the seed-gap transfer (§C); ionmc samples before any TOPAS run with exact seal pins (§D); the gates
(§F). All applied before the freeze; no run exists.

## Preserved obligations (verbatim from the V3-008 integration note, `validation/plans/v3-008-acceptance.md`)

"**Unresolved qualification obligations preserved:** row I3 (C/O stopping vs the revised ICRU 73 table) FAIL;
measured-data rows are evidence gaps (experimental validation deferred, directive 2026-11-07); carbon/oxygen nuclear
interactions, fragmentation, secondary ions and distal fragment dose are V3-009; the low ends and 430 MeV/u ends of the
domains were not re-measured after the R2-STOP3 model change; lineage identity of ICRU73_NEW with the 2009 revision
OPEN. Not a clinical claim."

## Status

- 2026-11-17 (planner, S0): draft freeze candidate; no code, run or comparison output exists for V3-013. Base
  a9563412 (V3-008 integrated, PR #84). Open before freeze: independent Codex review and scientific review of decision
  0044 and this plan. Decision points reported to the orchestrator: Q3 (MCS proton path fixed and recorded), the use of
  the V3-008 reservation gap 20610000–20618999 for the V3-013 seed families, and V-C1-A as report-only.
- 2026-11-17 (orchestrator, S0 revision): Codex REVIEW-0f80fad8 (changes_required) and the independent scientific
  review applied in place (section H); still a freeze candidate pending the re-review.
- 2026-11-17 (orchestrator, S2a amendment A-1; prospective, before any V-C run). (a) DOM-7 is implemented as frozen with
  one explicit hook: a `stopping` object whose type `spec_of` does not support is refused unless it declares itself by
  a method `stopping_spec()` returning a user spec (`decision = "user"`, recorded as `model_id_effective =
  "unregistered"`, never a registered id); the duck-typed test-only sources of `test_config_validation.py` and
  `test_scoring_transport.py` declare it. (b) The `E0 + 6σ_E` part of the window check of decision 0044 §3 (row DOM-3)
  fails closed at `validate` for registered models; for unregistered user sources the pre-existing V3-003 runtime rule
  is retained unchanged (`source_energy_out_of_range` counter, result invalid; `test_transport_reference.py::
  test_source_energy_out_of_range_is_counted`): such a run cannot extrapolate silently because its result is invalid,
  and the table-range, LET-water-row and deuteron checks still apply to it. Reason: a frozen row contradicted a
  tested pre-existing behaviour; recorded here rather than silently narrowed.
- 2026-11-17 (orchestrator): **frozen at 7c139ca1** (Codex REVIEW-65665aa0ec074fa285c5e65da13646be, pass, signed; scientific review two passes applied); S1 starts from this commit.
- 2026-11-17 (orchestrator, S5a): seal material part 1 committed (bundles manifest digest aa2f24e436f5cfb22b6bc0920367bfe7d2677bb00ea05cd9a05d407cd1466b83, predictions file sha256 f30ac1f7265ace91cbf341d337a98e5103da68b9adda558d032d79042f733663); no run yet.
- 2026-11-18 (orchestrator, S5b): seal material part 2 committed (evaluator sha256 77dba18a3d902260aa557a399587c6d36ee07286434b6745f0401eebb38917a3, env_digest sha256 c49f3106ea994e25df01b6cee215e2e0436fbd4d3a82d9f6c83e43984df2bc10); seal-environment.json pending the controlled host run.
- 2026-11-18 (orchestrator, seal candidate): `seal-environment.json` from the controlled host run RUN-20261118T120550Z-168e03a6
  (`env_digest.py --write` through `sealed_launch.py` at 4130d569); the host numpy digest (269e58d9…) differs from the
  sandbox's (f98a87d5…) while interpreter, stdlib, warp and uv.lock agree, so the seal pins are host-produced:
  `vc-predictions.json` replaced by the output of RUN-20261118T120833Z-627be664 (`steps_v13.py predict --cache-dir`
  at 6abda9eb; two earlier attempts RUN-20261118T120658Z-7ce8d7c7 and RUN-20261118T120706Z-13e646d1 failed at
  import, before any computation) — every value, table identity and table sha256 equals the sandbox file; only the
  pinned script sha256 changed (the `--cache-dir` option). Frozen manifest `v13-frozen-manifest.json` written
  (80 sealed paths, `src/ionmc` tree 68db1d3b…). This commit is the seal candidate; the seal is the signed Codex
  review of its exact SHA (pending the codex-worker service, IR-20261117-151541-ED9CA9). No evidence run exists.
- 2026-11-18 (orchestrator, gate records of the implementation steps; rows A.1–A.5 run by the lightweight test suite in
  CI and in the logged local runs): S1 3a23f6f3 → 8a3791bc (VAL-20261117-120611-2F373F, VAL-20261117-151354-0F8D19; CI
  run 44112144023 green at 8a3791bc after the D-3 digest regeneration; Codex REVIEW-2a71314b pass at 3a23f6f3), S2a
  dd8929dd (VAL-20261117-154631-2B79E5), S2b afe42b30 (VAL-20261117-155710-EFCF01; CI run 44120214581 green), S3
  983eb1dd (VAL-20261117-161952-D1A9AE; CI run 44124558686 green), S4 f524fe90 (VAL-20261117-163333-730832; CI run
  44127163573 green), S5a 5248730e (VAL-20261117-164538-0703D8; CI run 44129446092 green), S5b 4130d569
  (VAL-20261118-120527-2228BE), seal candidate f1dfa55d (VAL-20261118-121113-9DDF7E). Codex reviews of 8a3791bc and
  later are pending the codex-worker service (IR-20261117-151541-ED9CA9); the same-host sha256 part of row ID-DEF
  and the full R1-ion repeatability are local-validation items still to be run on the host runner before integration.
- 2026-11-18 (orchestrator, re-seal): the S6 documentation commit 368cad59 edited `src/ionmc/THIRD_PARTY_NOTICES.md`,
  which lies inside the pinned `src/ionmc` tree (68db1d3b… → 5be8b2b8…), so the f1dfa55d seal candidate's pins went
  stale before any review (found by the host test run RUN-20261118T121617Z-a9e14c06 and by `gen_v13_manifest.py
  --check`). `vc-predictions.json` re-derived on the host at 29eacc1b (RUN-20261118T122425Z-6065f566; every value,
  identity and table sha256 unchanged; only the tree pin changed) and the frozen manifest regenerated. This commit is
  the seal candidate; the seal is the signed Codex review of its exact SHA (pending IR-20261117-151541-ED9CA9). Rule
  adopted: nothing under `src/ionmc` (notices included) changes after the seal material is generated; a change means a
  new seal. No evidence run exists.
- 2026-11-18 (orchestrator, seal candidate 3): the host runner exports `GIT_WORK_TREE=/workspace` (diagnostic
  RUN-20261118T123758Z-3ee984a3), which redirected the `git -C <repo>` calls of the manifest generator and the
  evaluator when run on another repository (the synthetic-repo tests failed on the host, RUN-20261118T122616Z-a4ccfae3
  and RUN-20261118T123658Z-f24e84de; the real workspace checks were correct: tree 5be8b2b8…, clean). Both helpers now
  drop inherited `GIT_*` variables; PER-1's clean-process reader injects the package path explicitly (`python -I`
  ignores PYTHONPATH and the host venv has no site-packages install). Two sealed files changed → manifest regenerated.
  This commit is the seal candidate; the seal is the signed Codex review of its exact SHA (pending the service).
- 2026-11-18 (orchestrator, seal candidate 4 after the integration-base merge): the operator fast-forwarded the branch to
  the prepared merge commit 795de2a1 (parents 2d9e90e1 and `v3/develop` 8ca9d4fc = PR #93; IR-20261118-134659-22DDF8,
  resolved), as the review service's ancestry check requires. The sealed infrastructure closure changed with PR #93
  (`host_inputs.py`, `storage_guard.py` new; `host_snapshot.py`, `host_runner.py` changed), so the frozen manifest is
  regenerated; `src/ionmc` tree 5be8b2b8… and `vc-predictions.json` unchanged. **Row R1-ion (LV) PASS** at 2d9e90e1 on
  both backends (same `src/ionmc` tree as this commit): python 100 histories RUN-20261118T130803Z-6d0338d9 and
  warp-cpu 10 000 histories RUN-20261118T131049Z-5051c1b1 — A1/A2 (`stopping_model`) and B1/B2 (`stopping` object)
  bit-identical in every tally, end state, counter, energy balance, table block and identity block (`selected_by`
  differs as declared); results committed under `validation/results/stopping/v3-013/lv13-r1ion/`, the six signed
  host-run records of the seal and LV runs under `host-records/` (all verified against the trust file). Host test run
  RUN-20261118T124118Z-1e5bc05b at 2d9e90e1: the eleven V3-013 test files, 165 passed (host-marked ID-DEF sha256
  included). This commit is the seal candidate; the seal is the signed Codex review of its exact SHA.
- 2026-11-18 (orchestrator, seal candidate 5): Codex REVIEW-388775b2 at 9c12855d (changes_required) — registered
  tabulated identities were not bound to the live array content and the spec freeze was shallow — fixed at 4426833e
  (content digests in the specs, read-only arrays, `spec_of` recomputes from the live arrays; recursive freeze; docs
  path). The digest-pinned aggregate files (D-3, `cc-int.json`, `v-c3.json`) were regenerated with unchanged numerics.
  Because `src/ionmc` and the `icru90-water-v1` spec hash changed, `vc-predictions.json` was re-derived on the host at
  4426833e (RUN-20261118T162043Z-84864809): every range, prediction and gating entry unchanged; the pins (spec_sha256
  of `icru90-water-v1`, `src/ionmc` tree) changed. Frozen manifest regenerated. This commit is the seal candidate.
- 2026-11-18 (orchestrator, seal candidate 6): Codex REVIEW-c8c734b6 at 68f5b178 (changes_required, one important
  finding: a replaced start-range callable or a subclass overriding `table()` kept the registered identity) fixed at
  9aba6381 (registered recognition only for the exact composite and inner classes; the ICRU start-range values bound
  into the spec as `options["start_range_sha256"]`; NIST's start range is inside the digested `csda_range`).
  `cc-int.json`/`v-c3.json` regenerated (one input-hash line each, numerics unchanged). `vc-predictions.json`
  re-derived on the host at 9aba6381 (RUN-20261118T163350Z-73072441, record committed and verified): every value and
  gating entry unchanged; pins changed (`icru90-water-v1` spec_sha256 7b6cd304… → 4bc54369…, `src/ionmc` tree →
  5b129987…). Frozen manifest regenerated (82 sealed paths). This commit is the seal candidate.
- 2026-11-18 (orchestrator, seal candidate 7): Codex REVIEW-58783c86 at 8d1d150e (changes_required: an overridable trust
  predicate and a conditional start-range callable could keep the registered identity) fixed at eaaf6cea (trust decided
  inside `spec_of` by exact types of the composite and inner objects; the ICRU start rule is a dedicated class installed
  only by `build_icru90_source` and bound by object identity). Spec hashes unchanged; `cc-int.json`/`v-c3.json`
  regenerated (one input-hash line, numerics unchanged). `vc-predictions.json` re-derived on the host at eaaf6cea
  (RUN-20261118T164500Z-24846584, record committed and verified): only the `src/ionmc` tree pin changed (→ 7b753045…).
  Frozen manifest regenerated (82 sealed paths). This commit is the seal candidate.
- 2026-11-18 (orchestrator, seal candidate 8): Codex REVIEW-77a533e9 at 19b0323a (changes_required: per-instance method
  shadowing on the non-slotted frozen dataclasses) fixed at b17b5be1 (every trusted source class slotted; `spec_of`
  rejects instances with a `__dict__` and classes whose behaviour-defining attributes differ from the import-time
  snapshot). Spec hashes unchanged; the digest-pinned aggregates (D-3, `cc-int.json`, `v-c3.json`) regenerated with
  unchanged numerics. `vc-predictions.json` re-derived on the host at b17b5be1 (RUN-20261118T165917Z-bbc33f71, record
  committed and verified): only the `src/ionmc` tree pin changed (→ b65e6861…). Frozen manifest regenerated (82 sealed
  paths). This commit is the seal candidate.
- 2026-11-18 (orchestrator, seal candidate 9): Codex REVIEW-b01b9fa4 at 4fd38abb (changes_required: mutable start-rule
  slot state and a time-of-check/time-of-use gap between selection and table construction) fixed at 26da5997 (frozen,
  slotted start rule with read-only exact arrays; table digests computed from the arrays consumed inside `table()`;
  `validate()` recomputes the identity after every table is built and fails closed on any change). Spec hashes
  unchanged; D-3, `cc-int.json`, `v-c3.json` regenerated with unchanged numerics. `vc-predictions.json` re-derived on
  the host at 26da5997 (RUN-20261118T171417Z-226fd2bd, record committed and verified): only the `src/ionmc` tree pin
  changed (→ a6378d86…). Frozen manifest regenerated (82 sealed paths). This commit is the seal candidate.
- 2026-11-18 (orchestrator, seal candidate 10): Codex REVIEW-bcd2acbc at 3926190f (changes_required: the NIST table
  provider was not required to be the exact `StarTable` and NIST tables escaped the post-build content comparison)
  fixed at e8e01808 (exact nested providers for NIST and ICRU; every NIST table carries the digest of the arrays
  consumed inside `table()` outside its metadata; `validate()` compares every table of a registered tabulated model
  with the live and registered digests and fails closed). Spec hashes unchanged; D-3 regenerated (numerics
  unchanged). `vc-predictions.json` re-derived on the host at e8e01808 (RUN-20261118T172835Z-ca2c1648, record
  committed and verified): only the `src/ionmc` tree pin changed (→ 6ac381de…). Frozen manifest regenerated (82
  sealed paths). This commit is the seal candidate.
- 2026-11-18 (orchestrator, seal candidate 11): Codex REVIEW-dd6310ee at 80c24847 (changes_required: lossy normalisation
  of stateful construction fields and no consumed-output binding for analytic tables) fixed at a3c30226 (exact
  primitive types for every trusted field; `validate()` rebuilds every table of a registered selection from a fresh
  `resolve()` and requires bit-for-bit equality — the general closure of this finding class, about 6 ms per
  validation). Spec hashes and the digest-pinned aggregates unchanged. `vc-predictions.json` re-derived on the host
  at a3c30226 (RUN-20261118T175042Z-852732ab, record committed and verified): only the `src/ionmc` tree pin changed
  (→ 5cbd23a2…). Frozen manifest regenerated (82 sealed paths). This commit is the seal candidate.
- 2026-11-18 (orchestrator, seal candidate 12): Codex REVIEW-6483e470 at 8df08b18 (changes_required: the rebuild check
  ignored the selected cache directory and compared cache-local provenance) fixed at 5e0f29be (the cache directory is
  carried on the resolved source, the selection, the identity block and every `data_used` entry; the rebuild reads the
  same cache; `retrieved_at` is excluded from the metadata comparison). Because that directory is run provenance, the
  sealed identity pins carry it as the placeholder `<run-cache-dir>` and the evaluator normalises a run's block the
  same way while requiring a real path (01f40fa8). `cc-int.json`/`v-c3.json` regenerated (numerics unchanged).
  `vc-predictions.json` re-derived on the host at 01f40fa8 (RUN-20261118T180956Z-e7f93e4c; the intermediate host run
  RUN-20261118T180655Z-f7fed0c8 at 5e0f29be produced the un-normalised file that 01f40fa8 carried transitionally;
  both records committed and verified): values, gating and spec hashes unchanged; `src/ionmc` tree pin → 740e60bb….
  Frozen manifest regenerated (82 sealed paths). This commit is the seal candidate.
- 2026-11-18 (orchestrator, seal candidate 13): Codex REVIEW-f5d6f67e at aed1f072 (changes_required). Finding 2
  (the selected cache directory must be an absolute resolved path) fixed at d7b149b6 (`cc-int.json`/`v-c3.json`
  regenerated, numerics unchanged). Finding 1 (the normaliser turned the analytic default's null `cache_dir` into
  the placeholder and the evaluator then required a directory for every block, which would have invalidated every
  default-model evidence run) was fixed at 9a99846a by the orchestrator under the operator's explicit authorization
  (IR-20261118-182151-FD0238, option A) after the implementation worker's edit of the two sealed scripts had been
  denied by the permission classifier: both normalisers replace only string directories by `<run-cache-dir>` and
  keep null as null; the evaluator requires a real directory (non-empty, not the placeholder) exactly where the pin
  carries the placeholder and null where the pin is null, at the top level and per `data_used` entry; regression
  tests for every case. `vc-predictions.json` re-derived on the host at 9a99846a (RUN-20261118T193652Z-e14c09b0,
  record committed and verified): values, gating and spec hashes unchanged; `bethe-ionmc-v1` pins now carry
  `cache_dir: null`; script sha256 → 9bf72756…; `src/ionmc` tree pin → 54d226da…. Frozen manifest regenerated
  (82 sealed paths, sha256 3bc6abe9…). This commit is the seal candidate; no evidence run before its Codex review.
- 2026-11-18 (orchestrator, SEAL — §D step 3 pin note): seal review REVIEW-0e6ca486200b47e2a958376470ad50f5 at ad73d94eff7ce25dbb1db04b00f1136603627e03; manifest sha256 3bc6abe916eabfa2ba7ced9adfbf8926eeb32e0b910214aef754f2f67a1c3b39
  (signed codex-worker record, verdict pass, no findings, finished 2026-11-18T19:44:23Z; committed verbatim as
  `validation/results/stopping/v3-013/seal-review.json` and verified against the committed trust file). The sealed
  closure is the 82 paths of the frozen manifest plus the `src/ionmc` tree 54d226da… at the seal SHA; every evidence
  run must carry identical blobs and that tree. Next, in order (§D steps 4–6): this pin commit's own Codex review;
  the 18 ionmc samples (one controlled host run at a time, each started after this review finished); the samples
  commit and its review (second pin note `samples review <id> at <sha>`); the three TOPAS smokes; the evidence parts
  0–3 of the three cases; the sealed evaluator. No evidence run has been started. Not a clinical claim.
- 2026-11-18 (orchestrator, timing smokes after the seal; contrary evidence to §E): three `vc-sample` smokes at the
  non-evidence seed k = 99 (seeds 20612103 / 20612303), after the pin review REVIEW-7db57441 finished (19:59:31 UTC):
  case 1 `icru90-water-v1` 10 000 histories 52.0 s (RUN-20261118T200117Z-173deba7); case 3 `bethe-ionmc-v1` 10 000
  histories 93.2 s (RUN-20261118T200243Z-7502b1a2) and 50 000 histories 437 s (RUN-20261118T200438Z-5cb9e01e);
  all exit 0, `pass` true, every counter zero, identity blocks equal to the sealed pins, environment digests equal
  to `seal-environment.json`, identical `warp_modules`. Rate ≈ 4.5 ms/history (75 MeV/u) and ≈ 8.6 ms/history
  (200 MeV/u) on warp-cpu float64 with two threads, so one frozen 1e6-history sample needs ≈ 75 / ≈ 110 / ≈ 143 min
  for cases 1/2/3 — the §E estimate (1.5–3 h for all six) was ≈ 5× too low, and cases 2–3 exceed the controlled
  host runner's cap of 7200 s per run. No evidence seed (k = 0) has been used. Execution is paused pending the
  operator's decision IR-20261118-201459-7CFC95 (raise the runner cap, or a prospective amendment to chunked
  samples with a new seal). Not a clinical claim.
