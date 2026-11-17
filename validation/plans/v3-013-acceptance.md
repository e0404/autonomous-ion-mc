# V3-013 acceptance plan — explicit selection of stopping-power models and tabulated sources

> **Status: DRAFT freeze candidate (S0, planner, 2026-11-17, base `v3/develop` a9563412).** Frozen together with
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

### A.2 Consumer consistency (gating; synthetic test-only model, no reference data)

A test-only model `test-bethe-x1.05` (registered in the test module only; tabulated kind; S = 1.05 × the default's S on
the default grid, same start-range rule scaled consistently) is selected by `stopping_model` in a small water geometry
(python backend, float64, fixed seed from §C r_index 10).

| Row | What | Tolerance | Class |
|---|---|---|---|
| CC-R | CSDA range at 10 energies (proton, alpha) = default range / 1.05 | ≤ 1e-12 relative | self-consistency |
| CC-T | Mean depth of the cutoff points (proton 100 MeV, alpha 100 MeV/u; MCS and straggling off; 2 000 histories) = default mean depth / 1.05, same seed | ≤ 1e-9 relative (deterministic straight paths) | self-consistency |
| CC-LET | LET water row (decision 0040; `_water_table`) of the run = 1.05 × the default's row at every node | ≤ 1e-12 relative | self-consistency |
| CC-RES | `continuous_residual_v1` segment length for an alpha history = default / 1.05 | ≤ 1e-12 relative | self-consistency |
| CC-ID | Identity block names `test-bethe-x1.05`; every table identity of the run (materials, water row) names the same model; `TransportTables.sha256` differs from the default's | exact | self-consistency |
| CC-NEG | Monkeypatching one consumer at a time (material tables, water row, deuteron tables, exact-Ω ion rows) to the default source is detected by the single-model identity assertion of `validate` | each patch raises | capability |
| CC-HASH | New-model metadata carry `model_id`/`version` (hash differs); pre-existing models' metadata byte-identical; the metadata → model-id map is injective over the registry | exact | regression |

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

### A.5 `icru90-water-v1` consistency checks (gating for implementation correctness; **not evidence**)

| Row | What | Tolerance |
|---|---|---|
| CC-NODE | Resampled S reproduces the parsed ICRU 90 array at every node inside the window (proton, alpha) | ≤ 1e-12 relative |
| CC-INT | Interpolation error between nodes: log-log vs a monotone cubic in log-log space (both from the array only), max and RMS per decade | report-only (aggregates in Git) |
| CC-I2 | CSDA range above 1 MeV/u, R(E) − R(1 MeV/u), equals the I2 "integrated ICRU 90" construction of `compare_nist.py` (`_icru_range_table`) minus its value at 1 MeV/u (same interpolation and quadrature) at 50/75/100/125/150/200/250 MeV/u (alpha) and 100/150/200/250 MeV (proton); the start-range offset (own-array integral vs the I2 NIST start) is reported | ≤ 1e-12 relative; offset report-only |
| CC-SPAN | Loader test: counts and coverage of both grids (57 / 49 nodes; window inside the span); no value printed | exact |

### A.6 Configuration-specific validation of `icru90-water-v1` (graded once, after the seal of §D)

**V-C1 — TOPAS agreement, helium (gating for the model's evidence status `validated:V-C1`).**
- Reference: TOPAS 4.3 / Geant4 11.4.2, `g4em-standard_opt4` only, the TD-EM configuration of plan v3-008 (water,
  1 km electron cut, beam and scoring as TD-EM, IDD bins 0.02/0.1 mm), 4 parts × 2.5e4 primaries per case.
- Cases (fresh, held-out; decision 0044 Q6): **He 75, 125 and 200 MeV/u** (cases 1, 2, 3 of §C). None was run in any
  V3-008 row; all lie inside the He domain, away from the domain ends whose behaviour was not re-measured after
  R2-STOP3, and away from He 220 MeV/u (declared for the superseded, never-run V1-held row).
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
- Outcome per case: PASS / FAIL. The model's status is `validated:V-C1` only if all three cases PASS; otherwise
  `failed:V-C1` (recorded, never re-graded).

**V-C1-A — the default at the same energies (report-only).** The paired default-model samples of V-C2 graded with the
same criteria. They never change a V3-008 verdict, the default or any spec.

**V-C2 — selection reaches transport (gating for integration, §F).**
- Prediction, computed and committed in the seal before any run (model-only quantities, no reference value):
  ΔR_pred(E0) = R_CSDA,icru90(E0) − R_CSDA,default(E0) in mm of water for the three V-C1 energies (Gaussian-averaged
  over σ_E if the TD-EM source has σ_E > 0).
- Measurement: ΔR80 = R80(icru90) − R80(default) from paired ionmc samples with **identical seeds** (common random
  numbers, declared in §C), same estimator as TD-EM (a); σ by delete-one-batch jackknife over the 40 paired batches.
- PASS per energy iff |ΔR80 − ΔR_pred| + 1.645 σ ≤ max(0.02 mm, 0.10 |ΔR_pred|). Justification (ESTIMATE): the
  model-dependence of the R80 − CSDA offset through the Bohr-variance weighting is ≲ 0.005 mm; detour and binning
  terms cancel in the pair. The tolerance does not depend on any TOPAS output.

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
  and the default (common random numbers, V-C2); r_index 10 = A.2 test seed; r_index 11–14 = LV R1-ion/consumer runs;
  others reserved. Extent 20612004–20618103.
- Consumed by V3-008 and never re-drawn: every family listed above. Consumed by V3-013 once drawn: `topas13` parts used,
  `ion13` r_index used.
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
   definitions and the committed V-C2 prediction file. A run (smoke or evidence) is valid only if its `code_sha` has
   the seal commit as an ancestor and the seal review finished before the run's controlled `reference_started` /
   host-run start event (the V3-008 R2-TOPAS-4 rule; controlled records, not Git timestamps). The evaluator enforces
   both checks and fails closed.
4. **Smokes** (execution only), then **evidence runs** (TOPAS, then ionmc), one heavy job at a time.
5. **Evaluate once** at a committed SHA; commit the reports and run records; dated results note in this plan; ledger
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
- **Integration of V3-013** requires: all A.1–A.4 rows PASS, A.5 gating rows PASS, V-C2 PASS for all three energies,
  V-C1 graded (PASS or FAIL, never INCONCLUSIVE-by-omission: a run-validity failure after the frozen repeats is a FAIL),
  the documentation step S6 (default statement of decision 0044 §6, registry table, THIRD_PARTY_NOTICES update),
  the standard lifecycle gates (Codex review, exact-SHA validation record, CI, controlled merge).
- **No outcome changes the default or any spec within this task.** A FAIL is recorded and sets the model's evidence
  status; a revised alternative is a new `model_id`/`version` with fresh, held-out evidence pre-registered by amendment.
- A V-C2 FAIL blocks integration (selection does not reach transport as predicted) and triggers attribution by
  prospective amendment.

## G. Acquisition list

- S1–S4 and V-C1/V-C2/V-C3: **none new** (ICRU 90 file afd14794…, PSTAR f5fc6482…, ASTAR dd4a2e50… are cached and
  pinned).
- Only if a later task opens the low-energy C/O composite (S5b): `G4WaterStopping.cc` at Geant4 v11.4.2 via
  `acquire_reference_data` (role exploratory until a decision), followed by a written licence and provenance review
  (Geant4 licence; data provenance stated in the file) before any construction role. Optional, for lineage only:
  `G4ICRU90StoppingData.hh` and `G4BraggModel.cc` v11.4.2 (Q7 TO VERIFY items).
- libdEdx ICRU73_NEW as construction: only after an operator authorization (GPL-3); not requested by this plan.

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
