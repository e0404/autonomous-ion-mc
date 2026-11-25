# V3-007 acceptance plan — beamlet-resolved transport, sparse dose influence matrices, per-beamlet uncertainty, persisted influence bundles

> **Status: FROZEN at the S0 commit (2026-11-25) after independent review (summary in the Status section) — to be
> bound to its exact SHA by the S0 commit's Codex review.** Base `v3/develop` 3f3736a6 = the V3-006 squash merge,
> PR #95. Frozen together with decision 0046 **before any V3-007 code exists**. Later changes are dated amendments;
> earlier text is kept; failures are recorded and never re-graded. Not a clinical claim.

Markers: **FACT** (read at 3f3736a6), **ESTIMATE**, **TO VERIFY**, **PROPOSAL** (frozen by decision 0046). Row
classes: *regression*, *contract* (fail closed), *identity* (exact algebraic identity of the estimator),
*round-trip*, *determinism*, *statistical*, *measurement*, *workflow*. Tiers: **CI** (GitHub, small, warp-cpu, single
process, < 30 min for the whole tests job), **LV** (logged local exact-SHA validation), **HR** (controlled host runner,
the only CUDA path, RTX A6000). Evidence category of every row: **self-consistency** or **backend-parity** (PROTOCOL
taxonomy). **No row is physics evidence**: beamlet resolution re-partitions the same transported histories; the
physics of each history is unchanged and keeps the evidence status it has (proton EM: decisions 0038/0039; proton
nuclear: python reference only, decision 0041 with V3-005A diagnostic-mode grade; ions: decision 0043/0044 grades).
Single-process directive: `cpu_workers = 1` in every row; `IONMC_SINGLE_PROCESS=1`, thread limits 1.

## A. Requirement mapping (ledger rows advanced; exact evidence each needs)

V3-007 may move a row to "implemented" (code merged, evidence pending) only; never "evidenced"/"satisfied" without
conformant qualification (V3-012/V3-014+).

| Ledger row | Suite | What V3-007 delivers | Evidence required at the exact integration head | Remaining after V3-007 |
|---|---|---|---|---|
| V1-MUST-009 Beamlet-resolved scoring | S-CAP-PLANNING | spot-list source `BeamletSet`; per-history beamlet id derived from the global history index; per-beamlet energy/dose accumulators and per-beamlet energy balance on python / warp-cpu / warp-cuda | BEAM-1..BEAM-6, INF-1..INF-7 PASS (CI + LV); PAR-1..PAR-3 PASS (LV; CUDA rows HR) | qualification; beamlet-resolved LET/lookup channels (declared unsupported, §B.7) |
| V1-MUST-010 Dose influence matrices | S-CAP-PLANNING | sparse CSC influence matrix (Gy per primary of the beamlet), explicit threshold with per-beamlet discarded accounting; `matvec`, `rmatvec`, `to_dense`, `plan_dose`, `plan_uncertainty` helpers (numpy only) | SPA-1..SPA-5, INF-2, UNC-1..UNC-4, PER-1..PER-10 PASS | qualification; ion nuclear influence (V3-009); warp nuclear influence (V3-005B) |
| V1-MUST-001 Treatment-planning orientation | S-CAP-PLANNING | batched beamlet execution, device context reused across beamlets and across runs of one `Simulation` | INF-3, REUSE-1, PERF-1 (measurement) | qualification; frozen targets (V3-012) |
| V2-OUTPUT (clause "Influence matrices must retain geometry, beamlet mapping, units, history and normalization metadata") | S-PROV-OUTPUT | `ionmc.result` schema 1.1, `kind = "influence"`: beamlet table, history offsets, grid, units, normalisation, threshold, discarded accounting, per-beamlet uncertainty and balance; 1.0 bundles stay readable | PER-1..PER-10, OUT-0-INF clause map PASS (CI; PER-2 LV/HR on three backends) | conformant qualification; then V2-OUTPUT has no open clause of its own beyond qualification |
| V2-USER (clause "multi-beamlet and influence output") | S-WORKFLOW-USER | TOML `[source] type = "beamlets"` with inline `[[source.beamlets]]` or a spot file `{file, sha256}`, `[influence]`; `ionmc run/check/verify/report` for influence bundles; docs | CLI-1..CLI-5, DOC-1 PASS | clean-install qualification with data acquisition; heavy-ion workflow (V3-009/V3-010) |
| V2-CAP | S-CAP-CONTRACT | declared == enforced extended by the dimensions source kind {pencil, beamlets} × influence {off, on}; reconciliation rows for every new field; honoured check for influence | CAP-1..CAP-6 PASS (CI; CAP-1-LV/HR for cache-backed models and CUDA) | V3-005B, V3-009 extend the matrix |
| V1-MUST-030 Statistical uncertainty estimation | S-NUM-UNCERTAINTY | per-entry batch standard error; independence of beamlets by construction (disjoint counters) | UNC-1, UNC-4 | V3-011 coverage suite |
| V1-MUST-027 / V1-MUST-028 Random control, parallel-safe streams | S-NUM-UNCERTAINTY | no counter reuse across beamlets (decision 0037 encoding, `c0` never reset) | BEAM-4, INF-6 negative control, INF-7 | V3-011 |
| V2-NUM (clauses "seed/history/batch/backend changes", "history increases must improve statistics") | S-NUM-FALSIFICATION | per-beamlet history scaling, chunk/group invariance, backend parity of influence columns | UNC-3, INF-3, PAR-1..PAR-3 | V3-011 owns the adversarial suite |
| V1-MUST-034/035/036, V2-PERF | S-PERF-WORKLOADS | **measurement only** (PERF-1, PERF-2), labelled "EM-only, not the performance.json influence workload" | measurement recorded with hardware/config/method | frozen targets and full-physics workload (needs V3-005B; V3-012) |
| SHOULD: sparse representations, batched beamlet execution, planning-oriented reuse, planning-aware uncertainty, batch-based estimation | S-CAP-PLANNING / S-NUM-UNCERTAINTY | as above | SPA-*, INF-3, REUSE-1, UNC-* | **FACT: the ledger has no SHOULD rows**; S7 adds a "SHOULD coverage" note section (no new status vocabulary) |

Not addressed and not claimed: beamlet-resolved LET / lookup / species channels (fail closed), warp nuclear (V3-005B),
ion nuclear and fragments (V3-009), multi-process beamlet runs (deferred), any physics row, any evidence record,
performance targets.

## B. Design decisions (frozen by decision 0046; summary)

1. **Beamlet model.** `Beamlet(position_mm, direction, kinetic_energy_mev, energy_sigma_mev=0, lateral_sigma_mm=0,
   n_histories, weight=None, label=None)`; `BeamletSet(projectile, beamlets: tuple[Beamlet, ...])` is a second
   source type accepted by `SimulationConfig.source` (union with `PencilBeamSource`). One projectile per set. `weight`
   (planning weight, primaries) is metadata for `plan_dose`; it **never** changes sampling. Limits (fail closed):
   1 ≤ J ≤ 4 096 beamlets (ESTIMATE-justified: ≤ 400 B per beamlet inline in `requested` incl. a 64-character label
   ≈ 1.6 MiB at the limit, inside the 0045 `max_json_bytes` and token budgets; the frozen workload has 100),
   `label` optional, ≤ 64 characters, printable (`str.isprintable()`, no control characters) and unique within the
   set (A-1 item 3), every `n_histories` ≥ B and a multiple of `n_batches` B, Σ n_j ≤ 2³² (decision 0037),
   `RunOptions.n_histories` = Σ_j n_j (`validate` refuses a mismatch with `UnsupportedCombinationError` naming both
   values; TOML: absent → derived, equal → accepted, different → `ConfigError`; A-1 item 1), every beamlet
   inside the projectile's energy domain with the same
   ±6σ rules as a pencil source, mixed projectiles rejected. The existing `validate` constraints stay (FACT at
   3f3736a6): 1 ≤ N ≤ 2³², B ≥ 2, N % B == 0, accumulators `B·ΣV·8·workers ≤ memory_budget_bytes` (default 2³¹) and
   the per-batch quanta capacity check (`per_batch·e_cap/QUANTUM < MAX_QUANTA`).
   **Envelope rules** (a `BeamletSet` is consumed by code written for one source through its envelope): for the
   quanta capacity check `e_cap = max_j (E_j + 6σ_j)`; `max_steps` is derived from the envelope E_hi;
   `fdm_domain_status` and the σ ≤ fraction·E0 rule are evaluated **per beamlet**, not on the envelope.
2. **History mapping.** Beamlet j owns global history indices `[s_j, s_j + n_j)`, `s_j = Σ_{i<j} n_i` in spot-list
   order (decision 0037: "beamlet membership is a range of history indices", `c0` never reset). Batch = `h mod B`;
   since every s_j and n_j is a multiple of B, each beamlet contributes exactly n_j/B histories to every batch.
   Source sampling (purpose 1) reads beamlet j's parameters; transport streams are untouched.
3. **Execution (all backends, one kernel).** Per-history beamlet id `j = searchsorted(s, h, "right") − 1` on a small
   offsets array; per-beamlet source parameters in arrays (float64 energies/positions, `real` direction as today).
   Accumulators `int64[G·B, V_inf]` for the influence grid, row `(j − j0)·B + batch`; G = beamlet group size,
   chosen deterministically as the largest G with `G·B·V_inf·8 + existing accumulators ≤ memory_budget_bytes`
   (or user-set `influence.group_size`), recorded. Launch chunks may span beamlets inside a group. Integer
   accumulation ⇒ per-beamlet arrays are **bit-identical for any G and any chunk size** within a precision. After
   each group: one device→host copy, per-beamlet reduction and sparsification on the host, accumulators zeroed. A
   pencil source is the J = 1 case of the same code path (bit identity guarded by REG-1). The beamlet id is
   range-asserted (0 ≤ j < J) in every backend; it can never wrap.
4. **Influence representation.** Canonical CSC by beamlet: `indptr int64[J+1]`, `voxel int64[nnz]` (flat C-order
   index into the influence grid, strictly increasing within a column), `energy_mev`, `energy_std_mev`, `dose_gy`,
   `dose_std_gy` float64[nnz], `n_nonzero_batches int32[nnz]`, batch values `float64[B, nnz]` (persisted by default,
   O3; opt-out recorded). Units: MeV and Gy **per primary of beamlet j**. Helpers in numpy only.
5. **Thresholding.** Default `threshold_rel = 0.0`: lossless (every voxel with a nonzero integer accumulator is
   stored). A user threshold τ ∈ (0, 1) drops entries with `dose_gy < τ · max_v dose_gy[:, j]` per beamlet; it is a
   requested field, reconciled `equal`, persisted, and the per-beamlet discarded energy (MeV/primary), discarded
   fraction of in-grid energy, discarded entry count and maximum discarded dose are persisted arrays and reported.
   No statistical (relative-uncertainty) threshold (biased).
6. **Per-beamlet uncertainty.** Batch standard error per entry (decision 0040 `reduce_batches`, hpb_j = n_j/B);
   `defined`/`sparse` masks as for grids. Beamlets use disjoint counter ranges and are statistically independent,
   so for plan weights x: `Var(D_v) = Σ_j x_j² σ_vj²` (`plan_uncertainty`); no cross-beamlet covariance is persisted
   or needed. `plan_uncertainty` excludes entries dropped by a threshold; the persisted per-beamlet
   `influence__discarded_var_gy2` (Σ σ² over dropped entries) quantifies the omission (stated in the docs).
   Per-beamlet energy balance (initial, deposited, escaped, truncated, cutoff, unaccounted, in-grid,
   outside, quantization) is reduced exactly from the per-history tally rows.
7. **Capability contract.** Declared and enforced together: beamlets supported for every species/backend/precision
   the pencil source supports (EM), and with nuclear on python only (warp + nuclear already exits 4 until V3-005B).
   **Fail closed:** influence + any `TallyRequest` or lookup (beamlet-resolved channels not implemented); beamlets
   with `cpu_workers > 1` (single-process directive; untested path); mixed projectiles; n_j not a multiple of B;
   influence grid name unknown; τ outside [0, 1); influence requested with a pencil source; `persist_batches` without
   influence. Spot list without `[influence]` = composite run (allowed, no tallies in V3-007).
8. **Persistence.** `ionmc.result` **schema 1.1** (MINOR: optional blocks), `kind = "influence"` (reserved by 0045):
   new blocks `beamlets` (table, offsets, hashes of the spot list), `influence` (grid name, representation, units,
   normalisation "per primary of beamlet j", threshold, group size, discarded summary); arrays under the safe-name
   prefix `influence__*`. `EffectiveConfig.summary()` is byte-identical for every pencil configuration (REG-2); for a
   `BeamletSet` its `source` block is the beamlet-set summary and an `influence` key is added. For a `BeamletSet`,
   `unit_direction` is null-with-status
   (`not-applicable: beamlets`), a per-beamlet direction array `beamlet__direction` is persisted, and the report's
   profile axis is taken from beamlet 0 and labelled so (every reader of the source, F1 list in BEAM-2, is covered).
   The reader accepts 1.0 and 1.1; tier-1 binds every new invariant (PER-4). The composite grids of the same run are
   persisted as
   ordinary grids. **V3-006 row OUT-4's unknown-version control "1.1" moves to "1.2" by a dated note in plan v3-006**
   (intent preserved, not a re-grade).
9. **CLI.** `ionmc run` with `[source] type = "beamlets"`, inline `[[source.beamlets]]` or `spots = {file, sha256}`
   (CSV with a fixed header, read by `numpy.loadtxt` on bytes read once, recorded in the `InputManifest` as role
   `spots`; bounds before parsing: `lstat` size ≤ 1 MiB, line count ≤ 4 097, sha256 over exactly the bytes read;
   violation → `ConfigError`, exit 2; A-1 item 3; `[run] n_histories` is optional and derived = Σ n_j, A-1 item 1),
   `[influence]
   grid, threshold_rel, group_size, persist_batches`. `ionmc report`: beamlet table (id, label, E, position,
   n_histories, nnz, discarded fraction, median relative uncertainty of entries ≥ 50 % of the column maximum) and SVGs
   (composite depth dose; per-beamlet nnz and discarded fraction vs beamlet index). Exit codes unchanged (0045 §7).
10. **Reuse.** Inside a run: tables, materials, geometry, grids, nuclear setup and the compiled kernel are built
    once; per beamlet only source parameters and accumulator rows change. Across runs: a `Simulation` caches its
    device context; `Simulation.run(beamlets=subset)` re-runs a subset with the **original** offsets (identical
    streams, so columns reproduce bit-identically; results of runs with the same seed are never merged — `merge` is not
    offered in V3-007). **Subset semantics (A-1 item 4):** the result's `requested` is the full configuration plus the
    requested-only field `run.beamlet_subset` (sorted unique ids; carried into `effective.run.beamlet_subset`);
    `run.n_histories` = Σ_{j ∈ subset} n_j while `beamlets.offsets` stay the original full-set offsets; composite
    grids, balances, accounting and statistics are those of the transported subset only (normalised per transported
    primary); the matrix keeps J columns, untransported ones empty (nnz = 0, `beamlet__transported = False`,
    per-beamlet balances NaN with mask); `save` persists `influence.subset` with `kind = "influence"`; a subset
    bundle and a full bundle are never merged.

### Amendments

**Amendment A-1 (2026-11-25, S0, before any code; Codex REVIEW-1652a7c1 on 9805dfdb, changes_required; four
important findings, all adopted).** No code exists and the plan is bound only by the S0 commit whose review passes,
so the frozen text above is edited in place; this amendment records what changed and why.

1. **History-count invariant (finding 1).** For a `BeamletSet`, `RunOptions.n_histories` MUST equal Σ_j n_j.
   `validate` refuses a mismatch (`UnsupportedCombinationError`, message naming both values). In TOML
   (`ionmc.config/1`) `[run] n_histories` is optional for `type = "beamlets"`: absent → derived = Σ n_j
   (reconciliation row `derived: sum-of-beamlets`); present and different → `ConfigError` (exit 2); present and equal
   → `equal`. The persisted `run.n_histories`, capacity checks, composite normalisation and the transport loop
   therefore all see N = Σ n_j. Rows: BEAM-1, CAP-3, PER-1.
2. **Complete reader inventory (finding 2).** BEAM-2's static list and the decision's Context gain the readers
   `config.py:348`, `471–498` (the ion-mapping path; its failure branch at 498 reads the singular
   `source.kinetic_energy_mev` and MUST be converted to per-beamlet values), `526`, `564`, `588` and `1016`. Rule:
   every reader, including exceptional/failure paths, is in the list; the test's grep check (`config.source`,
   `.source.`, `source.kinetic_energy_mev`, `source.unit_direction`, `source.projectile`,
   `source.energy_sigma_mev`, `source.lateral_sigma_mm`, `source.position_mm`) must find no reader outside the
   list; a converted failure branch gets a test that provokes it with a BeamletSet (BEAM-1 twin). Rows: BEAM-1,
   BEAM-2.
3. **Resource bounds before allocation (finding 3).** `label` is optional, ≤ 64 characters, printable and unique
   within the set; the per-beamlet JSON size is therefore bounded (ESTIMATE ≤ 400 B incl. a 64-character label →
   ≤ 1.6 MiB at J = 4 096, inside the 0045 budgets, re-stated). Spot file: before parsing, `ionmc` checks `lstat`
   size ≤ 1 MiB and, after reading the bytes once, line count ≤ 4 097 (header + J ≤ 4 096) and the sha256 of exactly
   those bytes against the declared value, then parses the fixed header with `numpy.loadtxt` on the in-memory bytes;
   any violation → `ConfigError` (exit 2) naming the bound. Rows: BEAM-1, CLI-2, PER-10.
4. **Subset-run semantics (finding 4).** As stated in §B.10: `run.beamlet_subset` is a requested-only field
   (reconciliation row `requested-only`); `run.n_histories` = Σ_{j ∈ subset} n_j; offsets stay the original;
   composite results are those of the transported subset; the matrix keeps J columns with untransported columns
   empty; `save` persists `influence.subset`, `kind = "influence"`; tier 1 binds Σ_{transported} n_j =
   `run.n_histories` and "untransported columns are empty"; subset and full bundles are never merged (no merge API).
   Rows: INF-4, PER-4, CAP-3.

Not a clinical claim.

## C. Seeds and determinism rules (declared before any run)

- CI tests use fixed small seeds written in the tests (not evidence).
- LV/HR runs use the new family **`inf7`: seed = 20800000 + 100·r + k**, k ≤ 99, r ≤ 49 → extent 20800000–20804999.
  r assignment: 0 = INF/SPA reference set, python float64; 1 = same, warp-cpu float64; 2 = same, warp-cpu float32;
  3 = same, warp-cuda float32 (HR); 4 = UNC-2 replicas (k = 0, 1), warp-cpu float32; 5 = UNC-2 replicas, warp-cuda
  (HR); 6 = UNC-3 4× histories; 7 = nuclear python beamlets (CAP/BEAM-6); 8 = carbon beamlets warp-cpu; 9 = PERF
  runs (HR; k = repeat·10 + mode); 10 = CLI workflow; 11 = reproduction from bundle; others reserved. Consumed once
  drawn.
- **Collision check (executed in S0, worktree at 3f3736a6):** `git grep -nE
  '\b207[5-9][0-9]{4}\b|\b208[0-5][0-9]{4}\b|\b2086000[0-9]\b' -- ':!*.npz'`
  over every tracked file. Hits: `uv.lock:250` (a digit run inside a package URL hash, not a seed) and
  `validation/plans/v3-006-acceptance.md:90` (the V3-006 collision-check sentence quoting its own window
  20660000–20760000); the V3-005 fresh-sequence bases 20781004, 20791004, 20801004 appear in no tracked file (they are
  arithmetic members of an open-ended sequence, not literals). Conclusion: no decimal seed literal in
  20750000–20860000, so `inf7` (20800000–20804999) has no literal collision; the arithmetic exclusion is the dated
  note below. **FACT from plan v3-006 §C:** the V3-005 fresh sequence (base step 10 000, envelope base +
  0…20999) is open-ended; the windows of bases **20781004, 20791004 and 20801004** reach `inf7`. They are excluded by
  a dated note in plan v3-005 (shared seed-base ledger, 2026-11-25), as V3-006 did for its window; V3-005B restates the
  exclusion.
- Guard `tests/ionmc/test_v3_007_seed_families.py` (pure python): extent, disjointness from the V3-006/V3-008/V3-013
  families and `V3_005_FAMILIES`, the three exclusions, literals match constants.
- Determinism (exact, within one backend/precision/host): (a) per-beamlet arrays independent of G and
  `chunk_histories` (INF-3); (b) a J = 1 set equals the pencil run (INF-1); (c) re-running a subset reproduces its
  columns (INF-4); (d) bundles: as V3-006 §C (content digests equal; `result.json` equal apart from `/run/timings` and
  path set P). **Backend parity is statistical** (requirements: no bitwise cross-backend claim), except python vs
  warp-cpu float64 trajectory parity (decision 0037 trajectory class).

## D. Ordered steps (per step: commit → Codex review → exact-SHA validation record → push → CI)

Operator diagnostic rules apply to every step: one heavy job globally; logged pytest runs (`python -X faulthandler
-u -m pytest -v --capture=tee-sys --tb=long`, logs in `telemetry/test-logs/<UTC>-<label>/`); worktree clean while a
Codex review runs; delete generated LV/HR bundles after their digests and summaries are committed (data-retention
rule; REG-4: no bundle in Git).

| Step | Content | Rows that must pass at the step's commit |
|---|---|---|
| S0 | decision 0046, this plan, seed note in plan v3-005, OUT-4 note in plan v3-006, ledger "in progress" notes; no code | review only |
| S1 | `Beamlet`, `BeamletSet`, `InfluenceOptions`, `validate` rules, source-consumer audit (every `config.source` reader uses the beamlet envelope or per-beamlet values), summary/reconcile rows | BEAM-1..BEAM-3, CAP-1 (CI part), CAP-2..CAP-4, REG-2, REG-3 |
| S2 | reference backend: per-history beamlet id, per-beamlet accumulators, per-beamlet balance | BEAM-4..BEAM-6, INF-1, INF-2, INF-6, INF-7 (python) |
| S3 | Warp driver: device context, offsets/parameter arrays, grouped accumulators, host sparsification | INF-1..INF-5, INF-7 (warp-cpu), REG-1 |
| S4 | `ionmc.influence`: CSC container, threshold accounting, helpers, per-beamlet statistics | SPA-1..SPA-5, UNC-1, UNC-4 |
| S5 | schema 1.1: writer/reader/tier-1 bindings/tier-2 re-derivation, golden 1.1 fixture, v3-006 OUT-4 control moved | PER-1, PER-3..PER-10, OUT-0-INF, REG-5 |
| S6 | CLI/TOML/spot file, report rows and SVGs, `Simulation` reuse and subset runs | CLI-1..CLI-3, REUSE-1 (CI part), INF-4 |
| S7 | docs (results format 1.1, configuration, CLI, workflow, architecture/transport beamlet section), ledger notes incl. SHOULD coverage section | DOC-1, docs build `--strict` |
| S8 | LV/HR evidence driver `validation/scripts/influence/inf7_run.py` (pattern of `out6_run.py`), one heavy job at a time | PER-2, PAR-1..PAR-3, UNC-2, UNC-3, CLI-4, CLI-5, CAP-1-LV/HR, REG-1-HR, PERF-1, PERF-2 |

### Discriminating rows (defined before implementation; all gating unless marked *measurement*)

Reference configuration **R-INF** (used where a row says "R-INF"): water box 60×60×200 mm, influence grid 2 mm
(30×30×100), protons, 4 beamlets at (−10, 0), (+10, 0), (0, −10), (0, +10) mm lateral, energies 70 / 90 / 110 /
130 MeV, σ_E = 0.5 MeV, σ_lat = 3 mm, direction +z, B = 10, n_j = 2 000 (CI: 200, B = 4). Default physics,
`bethe-ionmc-v1`, nuclear off unless stated.

**BEAM — beamlet model and mapping**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| BEAM-1 | Config validation: each fail-closed rule of §B.1 and §B.7 provoked separately | `UnsupportedCombinationError` (or `ConfigError` for TOML) before transport, message names the field; each positive twin passes; includes (A-1) `n_histories` ≠ Σ n_j refused (message names both values), absent TOML value derived, equal value accepted; label of 65 characters / with a control character / duplicate refused; the converted ion-mapping failure branch (`config.py` 498) provoked with a BeamletSet and naming the beamlet | CI |
| BEAM-2 | Source-consumer audit: a test enumerates every reader of `config.source` (static list checked against `grep` of `src/ionmc`, reviewed; FACT at 3f3736a6, complete after A-1 item 2: `config.py` (348; 471–498 ion-mapping path, its failure branch at 498 converted to per-beamlet values; 526; 564; 588; 1016; max_steps E_hi 755–756; channel compile 807; nuclear bounds 921–943; the decision-0044 stopping-model window check 1035–1064; ion domain 1127–1153; sigma rules and ion mapping 1237–1294 incl. `fdm_domain_status(proj, mats, e0, σ)` at 1287; stopping identity species 594–611; summary 624; `unit_direction=src.unit_direction()` 1351), `reference.py` (projectile 193–260, sampling 840–860), `warp_driver.py` (ctl source fields 121–135, `ctl.dir` at 135), `kernels.py` (source sampling), `results.py` (producible 403; profile axis from `unit_direction` 725 and 1626; binding 1668–1675; stopping-window check on E0 + 6σ 1806; tally E_hi bound 2093), `_result_shapes.py` 335–345, `contract.py` 78–87, 308–320, 701 (`SOURCE_TYPES`), 939, `simulation.py` 736 (`capabilities()["sources"]`), `report.py` 335, `configio.py` 590 and 1162) and asserts for a beamlet set the envelope rules: `max_steps`, channel `E_hi`, nuclear energy bounds, ion mapping domain and energy-sigma rules use min/max over beamlets | every consumer covered; a beamlet with E above the domain is refused although the mean energy is inside; every reader, including exceptional/failure paths, is in the list: a grep check (`config.source`, `.source.`, `source.kinetic_energy_mev`, `source.unit_direction`, `source.projectile`, `source.energy_sigma_mev`, `source.lateral_sigma_mm`, `source.position_mm`) finds no reader outside it | CI |
| BEAM-3 | Offsets: `s_j`, batch membership and per-batch count n_j/B for unequal n_j | exact | CI |
| BEAM-4 | Stream disjointness: the (c0, c3) pairs used by source sampling over all beamlets of R-INF are unique and equal `range(Σ n_j)` (instrumented reference run) | exact | CI |
| BEAM-5 | Source sampling: the sampled energy/position of the first history of every beamlet equals a direct evaluation of the shared source functions with that beamlet's parameters and counter | bitwise (python), float32 class (warp) | CI |
| BEAM-6 | Nuclear python beamlets (r = 7, 2 beamlets, 150/200 MeV, nuclear on): per-beamlet balance closes (`source − deposited − escaped − truncated − residual`, nuclear signed terms, rel 1e-9); secondaries' deposits land in their primary's beamlet (sum identity INF-2 holds with nuclear) | exact identities | CI (small), LV |

**INF — identities of the beamlet estimator**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| INF-1 | J = 1 `BeamletSet` vs the equal `PencilBeamSource`, same seed/backend/precision | every grid array, balance and counter **bitwise** equal; influence column (τ = 0) equals the composite grid's nonzero entries bitwise | CI (python, warp-cpu), LV (all), HR (cuda) |
| INF-2 | **Composite identity**: R-INF with influence vs the same spot list without influence | (a) Σ_j integer accumulators `[B, V]` of the columns equal the composite run's integer accumulators **bitwise**; (b) float64: composite `energy_mev` vs Σ_j (n_j/N)·column energy (τ = 0) within relative 1e-12 per voxel (only float64 summation-order rounding; J ≤ 100 terms); (c) composite grid of the influence run bitwise equal to the composite-only run | CI, LV, HR |
| INF-3 | Group and chunk invariance: G ∈ {1, 2, J}; CI: J = 4, n_j = 600 (B = 4), `chunk_histories` ∈ {2¹⁰, 2¹¹}, so chunk boundaries fall inside beamlets; LV/HR keep R-INF with `chunk_histories` ∈ {2¹⁰, 2¹⁴} | per-beamlet integer accumulators and balances bitwise equal across all variants | CI (warp-cpu), LV, HR (cuda) |
| INF-4 | Column independence: `run(beamlets=[2])` and `run(beamlets=[0, 3])` reproduce those columns of the full run | bitwise; the reconciliation table of the subset result has exactly one `requested-only` row (`beamlet_subset`), its accounting closes on the transported histories (`run.n_histories` = Σ_{j ∈ subset} n_j), untransported columns are empty (nnz = 0, `beamlet__transported = False`), `save`/`load` round-trips it and `verify` exits 0 (A-1 item 4) | CI |
| INF-5 | Trajectory parity python vs warp-cpu float64 for the first 3 histories of beamlets 0, 1 and J−1 (across beamlet boundaries) | decision-0037 trajectory class (existing T1 tolerances) | CI |
| INF-6 | **Negative controls** (patched in the test): (i) `c0` reset per beamlet → BEAM-4 and INF-4 fail; (ii) beamlet id off by one at a boundary, or an accumulator row error with correct source parameters → INF-2(a) still passes; INF-7 fails (BEAM-5 and INF-1 catch only a wrong source assignment; shows INF-2 alone is not sufficient); (iii) accumulator row stride B±1 → INF-2(a) fails | each detected by the named row | CI |
| INF-7 | **Per-beamlet closure**: Σ_v column_j quanta·Q + outside_j + quantization_j = deposited_j, reduced from the tally rows of h ∈ [s_j, s_j + n_j), for every j; negative controls "accumulator row off by one" and "columns permuted" must FAIL this row | exact (rel 1e-12) | CI, LV, HR |

**SPA — sparse representation and threshold**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| SPA-1 | τ = 0 round trip dense ↔ CSC (`to_dense`, `from_dense`) | bitwise; nnz = count of nonzero integer accumulators summed over batches | CI |
| SPA-2 | Canonical form: indptr monotone, voxel strictly increasing per column, dtypes, no explicit zeros at τ = 0 | exact | CI |
| SPA-3 | Threshold accounting, τ ∈ {1e-4, 1e-2, 0.1}: every kept entry ≥ τ·max_j, every dropped < τ·max_j; dense − kept = discarded per beamlet (energy, rel 1e-12); reported fraction = discarded/in-grid energy; counts exact; τ = 0 reports exactly 0 | exact | CI, LV |
| SPA-4 | Helpers: `matvec(w)` vs `to_dense() @ w`, `rmatvec(y)` vs `to_dense().T @ y`, `plan_dose(weights)` | relative 1e-12 | CI |
| SPA-5 | Threshold is never silent: τ appears in requested/effective/report; changing τ changes `config_sha256`; a bundle whose discarded arrays are removed fails `load` | exact | CI |

**UNC — per-beamlet uncertainty**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| UNC-1 | Entry std = batch standard error recomputed from the persisted `[B, nnz]` batch values; masks per decision 0040 rules | relative 1e-9 (reader tolerance of 0045 A-1) | CI |
| UNC-2 | *statistical*, **coverage between independent replicas** (r = 4, k = 0/1 compared; HR r = 5): R-INF variant with J = 16 beamlets (4 energies × 4 positions), n_j = 10 000, B = 20. Entry set S_j is fixed by an independent pilot replica (r = 4, k = 2; HR r = 5, k = 2): defined, n_nonzero = B and relative uncertainty ≤ 0.05 in the pilot; the compared replicas are k = 0 and k = 1. z = Δ/√(s₁² + s₂²); per-beamlet coverage c_j(1) = fraction |z| ≤ 1, c_j(2) = fraction |z| ≤ 2. Reference p(1), p(2) = P(\|T_ν\| ≤ 1, 2) with ν = 2(B − 1) = 38 (computed by numpy quadrature in the test) | for both levels L ∈ {1, 2}: \|mean_j c_j(L) − p(L)\| ≤ max(3·sd_j(c_j(L))/√J, 0.03) (beamlets independent; 0.03 = declared margin for non-normal low-count batches, fixed now). **Power control:** recomputing with s scaled by 0.8 must FAIL the same criterion; if it passes, the row is INCONCLUSIVE (gating fail) | LV, HR |
| UNC-3 | *statistical*, history scaling (r = 6): R-INF J = 4, n_j = 4 000 (k = 0) vs 16 000 (k = 1) (B = 20); entry set S fixed by an independent pilot replica (r = 6, k = 2; same criteria as UNC-2) | median over S of s(16 000)/s(4 000) ∈ [0.45, 0.55] (expected 0.5) | LV |
| UNC-4 | `plan_uncertainty(x)` = √(Σ_j x_j² σ_vj²) exactly (rel 1e-12); docstring and results-format page state the independence argument (disjoint counters, decision 0037) | exact | CI |

**PER — persisted influence bundles (schema 1.1)**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| OUT-0-INF | Clause map: geometry, beamlet mapping (offsets, labels, parameters), units, histories per beamlet, normalisation, threshold, discarded accounting, uncertainty definition → JSON pointer/array | every clause resolved | CI |
| PER-1 | In-process round trip `load(save(r))` for R-INF (python, warp-cpu f64), τ = 0 and τ = 1e-2, `persist_batches` on/off (off recorded `not-persisted`); includes a bundle whose `requested.run` carries the derived N = Σ n_j (A-1 item 1) | arrays bitwise, JSON equal | CI |
| PER-2 | Clean-process round trip (`python -I`) | canonical dump equal | CI; LV/HR on three backends |
| PER-3 | Version handling: the frozen V3-006 1.0 golden fixture loads unchanged; a 1.1 golden synthetic fixture (KB-scale) loads and regenerates its report byte-identically; "1.2"/"2.0" refused; `kind = "influence"` with schema 1.0 refused | exact | CI |
| PER-4 | Tier-1 negative controls, one each: indptr non-monotone or ≠ nnz; voxel index ≥ V or unsorted or duplicate; nnz arrays of unequal length; negative std; batch values not reproducing mean/std; discarded accounting inconsistent; beamlet table ≠ `requested`; offsets ≠ cumulative n_j; columns swapped with balances unchanged and digests recomputed (fails through the INF-7 identity); `Σ_j (n_j/N)·column + discarded` ≠ composite in-grid energy (rel 1e-9); influence block missing for `kind = "influence"`; subset bundle with a non-empty untransported column; Σ_{transported} n_j ≠ `run.n_histories` | each `ResultFormatError` naming the defect | CI |
| PER-5 | Budgets: a catalogue declaring nnz = 2⁴⁰ or a 3-D batch array beyond `max_array_bytes` refused before allocation (V3-006 OUT-10 method, RSS bound 64 MiB) | exact | CI |
| PER-6 | Reproduction from the bundle (r = 11): `from_dict(requested)` + seed → every content digest equal | exact | CI (small), LV |
| PER-7 | Tier 2 `verify_effective` re-derives `effective` for a beamlet bundle | exact | CI (default model), LV |
| PER-8 | No table values (allowlist extended only by `influence__*`, `beamlet__*`) | exact | CI |
| PER-9 | Honoured check: influence requested ⇒ influence arrays present; absent request ⇒ none | exact | CI |
| PER-10 | Bundle at the J limit (J = 4 096, n_j = B = 2, warp-cpu) saves and loads under the default budgets; `requested` + `beamlets` stay within `max_json_bytes` and the token budget; the size is computed with 64-character labels at J = 4 096 and asserted within `max_json_bytes` (A-1 item 3) | exact | CI |

**CAP — capability contract**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| CAP-1 | V3-006 sweep (CAP-1a/1b dimensions) × source kind {pencil, beamlets} × influence {off, on} | declared == enforced; zero third outcomes | CI; LV (tabulated); HR (CUDA) |
| CAP-2 | Declaration derived (monkeypatched constant changes the `beamlets` entries of `capabilities()["sources"]`; no new top-level key) | exact | CI |
| CAP-3 | Reconciliation table extended (every new field equal/derived/unused/requested-only; `group_size` None → derived `memory-budget`; `run.n_histories` absent in TOML → derived `sum-of-beamlets`; `run.beamlet_subset` → `requested-only`, carried into `effective.run.beamlet_subset`); zero unlisted differences | exact | CI |
| CAP-4 | Dispatch: patching the driver to ignore beamlet ids or the threshold makes `Simulation`/`save` raise | exact | CI |
| CAP-5 | Exit codes through `ionmc check/run` for every §B.7 refusal (4), TOML errors (2) | exact | CI |
| CAP-6 | warp backends + beamlets + nuclear → exit 4 until V3-005B, message names V3-005B's absence; python + nuclear + beamlets accepted | exact | CI |

**PAR — backend parity of influence columns**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| PAR-1 | python f64 (r = 0, n_j = 2 000) vs warp-cpu f32 (r = 2): per-beamlet in-grid energy, z = Δ/√(σ₁²+σ₂²) | \|z\| ≤ 4 for every beamlet (J = 4) | LV |
| PAR-2 | warp-cpu f32 vs warp-cuda f32 (r = 2/3): per-beamlet in-grid energy \|z\| ≤ 4 and entry-wise coverage as UNC-2 | both | HR |
| PAR-3 | warp-cpu f32 vs warp-cuda f32, same seed: per-beamlet integer accumulators | report-only (not gating; float32 device math differs); recorded | HR |

**CLI, reuse, performance, regression, docs**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| CLI-1 | `ionmc run beamlets.toml` (inline) and with a spot file → bundle; `verify` 0; `report --check` 0 | exact | CI |
| CLI-2 | Spot file: sha256 mismatch → 2; malformed header/row → 2 naming the line; recorded in `InputManifest` role `spots`; a file above 1 MiB (sparse file, nothing read) and one above 4 097 lines refused with exit 2 before parsing, naming the bound; the sha256 is computed over the bytes read once (A-1 item 3) | exact | CI |
| CLI-3 | Report content: beamlet table, threshold line, discarded summary, "Not a clinical claim.", no `\bvalidated\b` in generated text | exact | CI |
| CLI-4 | Offline host run (r = 10) on warp-cuda; bundle verifies | exact | HR |
| CLI-5 | r = 0..3 bundles verify and differ only in backend/device/precision fields and array contents | exact | LV, HR |
| REUSE-1 | A second `run()` on the same `Simulation` performs no table construction and no geometry/table upload (instrumented counters), results bitwise equal to a fresh object's | exact | CI; HR measurement of setup time |
| PERF-1 | *measurement* (HR, r = 9, 3 repeats): EM-only proton, influence-workload geometry (120×120×300 mm, 1 mm transport, 2 mm influence grid), 100 beamlets × 10 000 histories, multiple energies/positions; modes grouped (G auto) vs sequential (G = 1, chunk boundary at every beamlet); metrics of performance.json incl. assembly and save time and bundle size; batch-count scaling B ∈ {2, 10, 40} recorded (performance.json) | recorded with hardware/config/method; **no target, no pass/fail**; labelled "not the performance.json influence workload (nuclear unavailable on Warp until V3-005B)" | HR |
| PERF-2 | *measurement*: scaling beamlets 1/10/100 at 10 000 histories each, B = 20 | recorded | HR |
| REG-1 | A16 digests and T1 trace fixture | bitwise | CI, LV, HR |
| REG-2 | `EffectiveConfig.summary()` of the V3-006 fixed config set | byte-identical (CI tolerance rule of plan v3-006 note 2026-11-23) | CI |
| REG-3 | `Result` fields and `capabilities()` keys of 3f3736a6 still present | exact | CI |
| REG-4 | No table values and no LV/HR bundle in Git | exact | CI |
| REG-5 | V3-006 rows stay green (the OUT-4 control moved to "1.2" by dated note) | exact | CI |
| DOC-1 | Results-format 1.1 page checked against schema constants; configuration reference lists every new key; docs build `--strict` | exact | CI |

CI budget: every CI row uses ≤ 4 beamlets × ≤ 200 histories on python/warp-cpu (ESTIMATE: < 4 min added); heavy rows
are LV/HR only.

## E. Independent review items

1. Codex review of every step's exact commit (lifecycle; approval is not scientific validation).
2. Independent reviewer (read-only, not the author) of S0 and of the final head: (a) the source-consumer audit is
   complete (any reader of the single source not converted silently uses beamlet 0); (b) INF-2 cannot pass with a
   wrong beamlet assignment unless INF-1/BEAM-5 catch it (INF-6 shows this); (c) UNC-2 power and the 0.03 margin;
   (d) threshold accounting can never be silent; (e) schema 1.1 strictness: a swap that is not coherent with the
   per-beamlet balances must fail tier 1 (PER-4 row: columns
   swapped, balances unchanged, digests recomputed → `ResultFormatError` via the INF-7 identity); a coherent forgery is
   detectable only by PER-6 (reproduction); (f) no LET/lookup silently accepted with influence; (g) privacy
   rules of 0045 (spot files appear by basename and sha256 only); (h) the plan makes no physics or clinical claim; (i)
   final-head items the S0 reviewer could not verify: (c2) the CUDA
   float32 lateral basis stays bitwise equal for J = 1 with per-beamlet arrays (the INF-1 HR row decides); (f2) tier 1
   recomputes `config_sha256` for API runs (results.py:2291 seen); (h2) the CI time ESTIMATE (< 4 min) and the cost of
   `validate` for LS ions (time locally with `--durations` before push), and disjointness from the V3-008/V3-013 seed
   families (the guard test covers it).
3. Contrary evidence preserved: a failing row is amended prospectively, never deleted or re-graded.

## F. Integration gates

- **STOP**: any REG difference; any change of a stopping spec, table hash, evidence record or of a V3-006 frozen
  outcome other than the dated OUT-4 note; any gating row FAIL without a prospective amendment reviewed
  independently.
- **Integration** requires every gating row of §D PASS at the exact head (CI rows in GitHub CI and in a logged local
  run; LV rows with exact-SHA records; HR rows from one controlled host run at the clean head, signed record
  verified), PERF-1/PERF-2 recorded, S7 documentation, ledger rows V1-MUST-009/010 set to "implemented" and the
  V2-OUTPUT/V2-USER/V2-CAP/V1-MUST-001/030 notes citing rows and exact SHAs (no "evidenced"/"satisfied"), decision
  0046 "Later validation outcome" filled, and the lifecycle gates: passing independent Codex review of the exact head
  and current base, exact-SHA local validation record, push, PR, exact-head CI green, controlled merge, worktree
  retired with branch preserved. "Not a clinical claim." on every evidence document.
- With V3-007 integrated the M1 exit criterion is **not** declared reached by V3-007 alone: proton nuclear runs on
  the python backend only (V3-005B deferred), which the checkpoint note must state.
- Deferred, recorded, non-gating: multi-process beamlet runs; full-physics influence workload on Warp (V3-005B);
  beamlet-resolved LET/lookup channels.

## G. Acquisition list

None. Cache-backed CAP-1 rows use datasets already cached; no download.

## Status

- 2026-11-25 (planner): draft; no code; open questions O-1..O-12 raised. Not a clinical claim.
- 2026-11-25 (orchestrator decisions): all planner recommendations adopted, see decision 0046.
- 2026-11-25 (independent pre-freeze review, Opus reviewer; verdict: freeze with listed changes, every finding
  adopted and applied in this text and in decision 0046): FACT corrections F1 (the reader list of `config.source` in
  BEAM-2 and the decision's Context was incomplete; the decision-0044 stopping-model window check at
  `config.py` 1035–1064 had been mislabelled "nuclear bounds"; `report.py` 335, not 332) and F2 (`validate`
  constraints now complete, §B.1). Findings: (1) INF-6(ii) did not hold: new row INF-7 and the beamlet-id range
  assertion (§B.3); (2) INF-3 CI parameters changed so chunk boundaries fall inside beamlets; (3) E.2(e) replaced by
  the tier-1/PER-6 split and the PER-4 "columns swapped" control; (4) the beamlet limit lowered from 65 536 to
  4 096 (0045 budgets) and row PER-10 added; (5) UNC-2/UNC-3 entry sets fixed by an independent pilot replica;
  (6) direction/summary semantics for beamlets (§B.8); (7) envelope rules (§B.1); (8) PERF-1 batch-count scaling and
  the `capabilities()["sources"]` wording of CAP-2; (9) `influence__discarded_var_gy2`. Reviewer's unverified items
  carried into §E.2(i). Not a clinical claim.
- 2026-11-25 (S0 fix 1, Codex REVIEW-1652a7c1 on 9805dfdb, changes_required; four important findings adopted, see
  Amendment A-1 after §B): (1) `n_histories` = Σ n_j invariant; (2) complete `config.source` reader inventory incl.
  the ion-mapping failure branch; (3) label and spot-file bounds before allocation; (4) subset-run semantics. Rows
  BEAM-1, BEAM-2, CAP-3, CLI-2, INF-4, PER-1, PER-4, PER-10 updated. Not a clinical claim.
