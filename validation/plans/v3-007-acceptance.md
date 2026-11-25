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
| V1-MUST-001 Treatment-planning orientation | S-CAP-PLANNING | batched beamlet execution, device context reused across beamlets and across runs of one `Simulation` | INF-3, REUSE-1, MEM-1, PERF-1 (measurement) | qualification; frozen targets (V3-012) |
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
   `label` optional, matching `^[A-Za-z0-9_.-]{1,64}$` and unique within the set (A-1 item 3, grammar A-5 item 3),
   every `n_histories` ≥ B and a multiple of `n_batches` B, Σ n_j ≤ 2³² (decision 0037),
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
   with 1 ≤ G ≤ J (an explicit `influence.group_size` outside [1, J] is refused). `memory_budget_bytes` bounds the
   accumulator payload of a run by the exact formula below (A-3 item 2, A-4 item 2, narrowed by A-5 item 1).
   Automatic G = the largest value in [1, J] satisfying the bound; if G = 1 does not fit,
   `validate` refuses before any allocation (`UnsupportedCombinationError` naming the budget and the requirement at
   G = 1); an explicit G that exceeds the bound is refused the same way; the effective G is
   recorded. Launch chunks may span beamlets inside a group. Integer
   accumulation ⇒ per-beamlet arrays are **bit-identical for any G and any chunk size** within a precision. After
   each group: one device→host copy, per-beamlet reduction and sparsification on the host, accumulators zeroed. A
   pencil source is the J = 1 case of the same code path (bit identity guarded by REG-1). The beamlet id is
   range-asserted (0 ≤ j < J) in every backend; it can never wrap.

   **Allocation rule (A-4 item 2; FACT: array names read in `transport/warp_driver.py` and `transport/reference.py`
   at 3f3736a6).** `estimate(backend, B, ΣV, chunk, V_inf, G) ≤ memory_budget_bytes`, where `estimate` is
   `ionmc.config.allocation_estimate` (implemented in S1). **Scope (A-5 item 1):** it bounds the ACCUMULATOR PAYLOAD
   of a run, exactly: the grid accumulators `B·ΣV·8`, the per-chunk tally rows and counter rows of ONE chunk (chunk =
   min(`chunk_histories`, length of the CURRENT group's history range), as `run_warp_range` sizes them), the influence
   accumulators `G·B·V_inf·8` on the device and one simultaneous host copy of them during transfer (python path: no
   second copy). It does NOT bound, and the documentation (S7, configuration reference) says so: diagnostic arrays
   (traces, end positions, escape records; bounded by their own limits such as `trace_histories`), the retained
   CSC/batch output of earlier groups and the assembled `Result` (bounded at `save` by the 0045 `max_array_bytes` /
   `max_total_bytes` budgets, which refuse an oversized bundle, and measured by PERF-1/PERF-2 as bundle size),
   Python object overhead, and the Warp kernel/module cache. Automatic G is a payload rule, not a guarantee about
   process RSS. Symbols: B = `n_batches`; S = ΣV, the total voxels of all scoring grids
   (`total_vox`); n_g = number of scoring grids; V_inf = voxels of the influence grid; G = group size; N = Σ_j n_j;
   chunk = min(`RunOptions.chunk_histories`, max(N, 1)) (`chunk` in `run_warp_range`). Row widths are the named
   constants and functions of the code, never re-typed numbers (the values at 3f3736a6 are in parentheses).
   Scope: the rule applies to `BeamletSet` runs (with or without influence; V_inf = 0 and G = 1 without);
   pencil configurations keep the existing `B·ΣV·8·workers ≤ memory_budget_bytes` check unchanged (REG-2, REG-3):

   - tally-row width W_t, warp: `N_FIXED_TALLIES` (= `len(TALLY_NAMES)`, 6) + 2·n_g + `channel_columns(eff)`
     (`n_cols` of `run_warp_range`); python: `N_FIXED_TALLIES` + 2·n_g + `channel_columns(eff)` +
     `len(NUCLEAR_TALLY_NAMES)` (6) when nuclear is on (`tally_rows` of `reference.py`, 580–584);
   - counter-row width W_c, warp: `len(COUNTER_NAMES)` (9) + `len(ION_COUNTER_NAMES)` (2) for an ion
     (`counter_names`, `counter_rows` of `run_warp_range`); python: `len(self.counter_names)` of `reference.py`
     (237–243): `COUNTER_NAMES` + `ION_COUNTER_NAMES` for an ion, `COUNTER_NAMES` + `NUCLEAR_COUNTER_NAMES` (3) when
     nuclear is on;
   - F = bytes of every other array the backend holds during the run (geometry `material_index` int32[S_geo] and
     float64 densities, the uploaded tables, grid descriptors and offset arrays, the per-beamlet parameter and offset
     arrays, the non-diagnostic dummy `end_state`/`end_code`/`trace_*` arrays); the function sums the `nbytes` of
     the same arrays, so F follows the code instead of a typed constant.

   Warp paths (cpu and cuda), peak = persistent device arrays + the largest of three phases:
   `P = 8·B·S` (`edep` int64[B, S]) `+ 8·chunk·W_t` (`tally_rows` float64) `+ 4·chunk·W_c` (`counter_rows` int32)
   `+ F + 8·G·B·V_inf` (device influence accumulators int64[G·B, V_inf]); chunk phase
   `X_chunk = 8·chunk·W_c` (`crows.astype(int64)`) plus, on cuda only (`numpy()` is a view on cpu),
   `8·chunk·W_t + 4·chunk·W_c` (host copies of the rows); group phase `X_group = 8·G·B·V_inf` (ONE simultaneous host
   dense copy of one group during the device→host transfer); end phase `X_end = 2·8·B·S` (`edep_all` and the
   `edep_grids` copies). `estimate_warp = P + max(X_chunk, X_group, X_end)`.

   Python reference path (whole range N in one call, existing behaviour, now counted): `estimate_py = 8·B·S`
   (`self.edep`) `+ 8·N·W_t` (`tally_rows`) `+ 4·N·W_c` (`counter_rows`) `+ 89·N` (`end_pos`, `end_dir`,
   `end_ctrl_sum` float64[N, 3], `end_energy`, `end_ctrl` float64[N], `end_code` int8[N]: 24 + 24 + 24 + 8 + 8 + 1
   bytes per history, allocated unconditionally) `+ F + 8·G·B·V_inf` (host accumulators, no second copy)
   `+ 8·N·W_c` (`counter_rows.astype(int64)` in the reduction). Python lists and float objects of the exact
   summation (`exact_components`) are not arrays and are not counted; their size is bounded by one column of one
   chunk or range at a time (S1 states the bound in the docstring and MEM-1 keeps the 1.5 factor honest). For
   configurations WITHOUT a `BeamletSet` the `validate` rule stays as at 3f3736a6 (`B·ΣV·8·workers`), so REG-1..REG-3
   are untouched; the formula governs beamlet runs with influence.
   Group selection: G = the largest value in [1, J] with `estimate ≤ memory_budget_bytes`; G = 1 over budget →
   refused before allocation, naming the budget and `estimate` at G = 1; an explicit G outside [1, J] or over the
   budget → refused. Diagnostics are refused with influence (§B.10), so the diagnostic arrays are not in the formula.
4. **Influence representation.** Canonical CSC by beamlet: `indptr int64[J+1]`, `voxel int64[nnz]` (flat C-order
   index into the influence grid, strictly increasing within a column), `energy_mev`, `energy_std_mev`, `dose_gy`,
   `dose_std_gy` float64[nnz], `n_nonzero_batches int32[nnz]`, batch values `float64[B, nnz]` (persisted by default,
   O3; opt-out recorded). Units: MeV and Gy **per primary of beamlet j**. Helpers in numpy only.
5. **Thresholding.** Default `threshold_rel = 0.0`: lossless (every voxel with a nonzero integer accumulator is
   stored). A user threshold τ ∈ (0, 1) drops entries with `dose_gy < τ · max_v dose_gy[:, j]` per beamlet; it is a
   requested field, reconciled `equal`, persisted, and the per-beamlet discarded energy (MeV/primary), discarded
   fraction of in-grid energy, discarded entry count and maximum discarded dose are persisted arrays and reported.
       No statistical (relative-uncertainty) threshold (biased). **Zero-deposition columns (A-2 item 4):** a transported
    beamlet with no deposit in the influence grid (outside or non-intersecting source) is legal: nnz = 0, no
    threshold applies (column maximum 0), per-beamlet `discarded_energy_mev`,
    `discarded_count`, `discarded_max_dose_gy`, `discarded_var_gy2` are 0 and `discarded_fraction` is UNDEFINED,
    persisted as NaN with the mask `influence__discarded_fraction_defined` (decision 0040), at τ = 0 AND at every τ
    > 0 (A-3 item 4: at τ = 0 the fraction is exactly 0 for every column with positive in-grid energy); the report
    distinguishes only "defined" from "undefined (no deposition in the influence grid)" and prints
    that string and "undefined (no entries)" (median relative uncertainty); the
    per-beamlet energy balance still closes (outside the grid / escaped). **Composite identity under a
    threshold (A-4 item 3):** `Σ_j (n_j/N)·(kept_j + discarded_j)` = composite in-grid energy (rel 1e-9), with
    `kept_j` = Σ_v column_j energy and `discarded_j` the per-beamlet discarded energy; per-beamlet quantities are
    per primary of beamlet j and carry the weight n_j/N only inside this composite sum.
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
   The reader accepts 1.0 and 1.1; tier-1 binds every new invariant (PER-4), including the thresholded composite
   identity `Σ_j (n_j/N)·(kept_j + discarded_j)` = composite in-grid energy (A-4 item 3).
   The composite grids of the same run are
   persisted as
   ordinary grids. **V3-006 row OUT-4's unknown-version control "1.1" moves to "1.2" by a dated note in plan v3-006**
   (intent preserved, not a re-grade).
9. **CLI.** `ionmc run` with `[source] type = "beamlets"`, inline `[[source.beamlets]]` or `spots = {file, sha256}`
   (CSV with a fixed header, read by `numpy.loadtxt` on bytes read once, recorded in the `InputManifest` as role
   `spots`; read bound to one descriptor (A-2 item 3): opened with `O_NOFOLLOW | O_RDONLY | O_CLOEXEC` (a symbolic
    link is refused, `ConfigError` "symbolic link"), `fstat` on the descriptor must report a regular file of size
    ≤ 1 MiB,
    exactly `st_size` bytes are read from THAT descriptor (the V3-006 `_read_exact` pattern: a shrink fails, nothing
    beyond `st_size` is read), a SECOND `fstat` on the same descriptor must then report the same size and inode (any
    change → `ConfigError` "spot file changed while being read"; A-3 item 3), the sha256 is computed over those
    bytes, line count ≤ 4 097, then `numpy.loadtxt`
    parses the in-memory bytes; violation → `ConfigError`, exit 2; A-1 item 3; `[run] n_histories` is optional and
    derived = Σ n_j, A-1 item 1),
   `[influence]
   grid, threshold_rel, group_size, persist_batches`. `ionmc report`: beamlet table (id, label, E, position,
   n_histories, nnz, discarded fraction, median relative uncertainty of entries ≥ 50 % of the column maximum) and SVGs
   (composite depth dose; per-beamlet nnz and discarded fraction vs beamlet index). Exit codes unchanged (0045 §7).
10. **Reuse.** Inside a run: tables, materials, geometry, grids, nuclear setup and the compiled kernel are built
    once; per beamlet only source parameters and accumulator rows change. Across runs: a `Simulation` caches its
    device context. **Diagnostic column API (A-3 item 1; supersedes the A-1 item 4 / A-2 item 1 subset runs):**
    `Simulation.influence_columns(beamlets: Sequence[int]) -> ColumnSet` re-transports exactly the histories of the
    named beamlets with the ORIGINAL offsets (bit-identical streams) and returns a `ColumnSet`: the CSC columns of
    those beamlets with their batch values, std, nnz, discarded accounting and per-beamlet balances, plus the beamlet
    ids and the configured offsets. A `ColumnSet` is NOT a `Result`: no composite grids, no statistics block, no
    traces (a configuration with `requested_trace_count > 0` is refused), and `save` refuses it (`TypeError`); it
    serves column-independence checks (INF-4) and the V3-011 hooks only.
    **Validity (A-5 item 2):** a `ColumnSet` carries the transport counters of its run (the same counter block a
    `Result` has) and is only ever returned when every invalidating counter is zero; otherwise the API raises
    `TransportLimitError` (or the existing invalidating error) unconditionally. `allow_invalid_result` is one of the
    refused fields: a configuration with it set is refused before transport (allowlist rule below).
    Invalid subsets (empty, duplicate ids, ids
    outside [0, J), non-integers) are refused by the API (`UnsupportedCombinationError`). Every bundle is a full run:
    `RunOptions.n_histories` == Σ_{all j} n_j is the only history count, and results of runs with the same seed are
    never merged (`merge` is not offered in V3-007).
    **`influence_columns` is fully fail-closed (A-4 item 1).** The API requires an effective configuration with a
    `BeamletSet` AND `InfluenceOptions`. It is refused (`UnsupportedCombinationError`) for a pencil source, for a
    configuration without `InfluenceOptions`, and for ANY requested diagnostic: `requested_trace_count > 0` (the
    config field `diagnostics.trace_histories`), `track_end_positions`, `escape_records`, and every other
    diagnostic or tally request that a `Result` would carry. The test enumerates the diagnostic and tally fields of
    `RunOptions`/`SimulationConfig` (the `DiagnosticsOptions` dataclass, `tallies`, lookups) from the dataclass
    fields at the time of S6, so a later-added diagnostic cannot be silently discarded: any field outside the
    explicit allowlist of "fields `influence_columns` honours" is refused.

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
4. **[superseded by A-3] Subset-run semantics (finding 4).** [superseded by A-3] As stated in §B.10:
   `run.beamlet_subset` is a requested-only field
   (reconciliation row `requested-only`); `run.n_histories` = Σ_{j ∈ subset} n_j; offsets stay the original;
   composite results are those of the transported subset; the matrix keeps J columns with untransported columns
   empty; `save` persists `influence.subset`, `kind = "influence"`; tier 1 binds Σ_{transported} n_j =
   `run.n_histories` and "untransported columns are empty"; subset and full bundles are never merged (no merge API).
   [superseded by A-3: subset runs are the diagnostic `influence_columns` API and produce no `Result`.]
      Rows: INF-4, PER-4, CAP-3.

**Amendment A-2 (2026-11-25, S0, before any code; Codex REVIEW-15a83983 on 0901763, changes_required; four
important findings, all adopted).** The frozen text is edited in place as for A-1; where it conflicts with A-1 item 4
(`run.n_histories` = Σ_{transported} n_j) or the A-1 item 3 spot-file rule, A-2 supersedes.

1. **[superseded by A-3] Configured vs. transported history counts (finding 1).** [superseded by A-3]
   `RunOptions.n_histories` == Σ_{all j} n_j always, also
   in the `requested` configuration of a subset result (so `from_dict(requested)` validates and
   `run(beamlets=requested.run.beamlet_subset)` reproduces it in a clean process). The transported count is the
   separate persisted field `run.histories_transported` (= Σ_{j ∈ subset} n_j; equal to `n_histories` without a
   subset); composite grids, balances, accounting and statistics are normalised per transported primary and tier 1
   binds `histories_transported` = Σ_{transported} n_j. `config_sha256` hashes `requested` and therefore differs
   between the full run and a subset run. Invalid subsets (empty, duplicate ids, ids outside [0, J), non-integers)
   are refused before transport (`UnsupportedCombinationError`). Rows: §B.10, INF-4, PER-4, PER-6.
   [superseded by A-3: only the configured count `RunOptions.n_histories` == Σ_{all j} n_j remains.]
2. **AST-based source-consumer audit (finding 2).** BEAM-2's grep rule is replaced by an exhaustive `ast` audit with
   alias following and a static allowlist (`envelope` / `per-beamlet` / `refused-in-S1` / `not-a-source`) that is
   the inventory. Additions: `config.py` 844, 1163, 1214–1229, `warp_driver.py` 187; `configio.py` 1162 moves to
   `not-a-source`. The allowlist is re-derived (and the decision's Context updated) whenever S2/S3 add readers.
   Row: BEAM-2.
3. **Spot-file read bound to a descriptor (finding 3).** §B.9 and CLI-2: `O_NOFOLLOW | O_RDONLY | O_CLOEXEC`,
   `fstat` regular file ≤ 1 MiB, exactly `st_size` bytes read from that descriptor, sha256 over those bytes, line
   count ≤ 4 097, `numpy.loadtxt` on the in-memory bytes; four adversarial controls (symlink, sparse oversized,
   rename between open and read, growth after `fstat`; the growth control is replaced by A-3 item 3).
4. **Zero-deposition transported columns (finding 4).** Legal: nnz = 0, `beamlet__transported = True`
   [superseded by A-3: no such field; every column is transported], discarded
   counts/energy/max dose/variance 0, `discarded_fraction` undefined (NaN + mask
   `influence__discarded_fraction_defined`), report strings "undefined (no deposition in the influence grid)" and
   "undefined (no entries)"; the per-beamlet energy balance still closes; UNC-2/UNC-3 entry sets exclude such
   beamlets. New variant R-INF-Z. Rows: §B.5, SPA-3, PER-1, PER-4, CLI-3, UNC-2/UNC-3.

**Amendment A-3 (2026-11-25, S0, before any code; Codex REVIEW-1d19386a on dd5ad774, changes_required; six important
findings, all adopted).** The frozen text is edited in place as for A-1/A-2; superseded sentences of A-1/A-2 stay and
are marked "[superseded by A-3]". Design decision behind items 1-3: subset runs no longer produce a `Result` or a
bundle; they become a diagnostic column API.

1. **Subset runs are a diagnostic column API, not a result (findings 1, 2, 3).** §B.10: `Simulation.influence_columns(
   beamlets: Sequence[int]) -> ColumnSet` as above; `run.histories_transported`, `run.beamlet_subset`,
   `influence.subset`, `beamlet__transported` and the "untransported columns" notion are removed everywhere. The
   configured invariant `RunOptions.n_histories` == Σ_{all j} n_j is the only history count; every bundle is a full
   run; the PER-4 identity `Σ_j (n_j/N)·(kept_j + discarded_j) = composite in-grid energy` [corrected by A-4 item 3: the
   discarded term is weighted by n_j/N] holds with N = Σ_j n_j. Rows:
   INF-4, PER-4, PER-6 (full-run variant only; r = 11, k = 1 released, marked unused), CAP-3 (no `requested-only` row
   for subsets), PER-1, CLI-3, REUSE-1 unchanged.
2. **Fail-closed influence grouping (finding 4).** §B.3: 1 ≤ G ≤ J; the bound counts every simultaneous accumulator [superseded by A-5 item 1: accumulator-payload budget];
   automatic G = largest feasible in [1, J]; G = 1 not fitting, or an explicit G over the bound, is refused by
   `validate` before allocation. Row: BEAM-1.
3. **Spot-file size change (finding 5).** §B.9: a second `fstat` on the same descriptor after the read must report
   the same size and inode, else `ConfigError` "spot file changed while being read". Row: CLI-2 (growth control).
4. **τ = 0 and zero-deposition columns (finding 6).** §B.5, SPA-3: at τ = 0 the discarded fraction is exactly 0 for
   every column with positive in-grid energy and UNDEFINED (NaN + mask false) for a column with no in-grid
   deposition, at every τ; R-INF-Z is exercised at τ = 0 AND τ = 0.1 in SPA-3, PER-1, PER-4 and CLI-3.

**Amendment A-4 (2026-11-25, S0, before any code; Codex REVIEW-a52098eb on af4aebc6, changes_required; four
important and two minor findings, all adopted).** The frozen text is edited in place as for A-1..A-3.

1. **`influence_columns` is fully fail-closed (finding 1).** §B.10: effective configuration with a `BeamletSet` AND
   `InfluenceOptions` required; refused for a pencil source, for no `InfluenceOptions` and for ANY requested
   diagnostic (`requested_trace_count > 0`, `track_end_positions`, `escape_records`, every other diagnostic/tally
   field a `Result` would carry); the test enumerates the dataclass fields at S6 and refuses every field outside the
   explicit "honoured" allowlist. Row INF-4 gains the capability controls (pencil source, no `InfluenceOptions`, each
   diagnostic field); the allowlist test fails when a config field is neither honoured nor refused.
2. **Exact backend-aware allocation formula and an allocation audit (finding 2).** §B.3 gives the per-backend formulas
   (`allocation_estimate`, S1) with the symbol definitions and the row widths read from the driver constants
   (`N_FIXED_TALLIES`, `channel_columns`, `COUNTER_NAMES`, `ION_COUNTER_NAMES`, `NUCLEAR_TALLY_NAMES`,
   `NUCLEAR_COUNTER_NAMES`); warp: device accumulators `G·B·V_inf·8` plus ONE host dense copy of one group, per-chunk
   tally AND counter rows; python reference path: tally and counter rows for the whole range N plus host accumulators
   `G·B·V_inf·8`, no second copy. New row MEM-1 (instrumented allocation audit; estimate within a factor 1.5 of the
   measured peak, a looser estimate is also a fail). BEAM-1 keeps the arithmetic group-selection cases. Steps: S1
   (function), S2/S3 (audit), S8 (LV warp-cpu, HR warp-cuda). §A evidence column of V1-MUST-001 gains MEM-1.
3. **Thresholded composite identity (finding 3).** PER-4 (and the tier-1 binding of §B.5/§B.8): `Σ_j (n_j/N)·
   (kept_j + discarded_j)` = composite in-grid energy (rel 1e-9); the reader control uses τ = 1e-2 AND unequal n_j
   (R-INF with n_j ∈ {2 000, 4 000, 2 000, 4 000}) so an implementation that does not weight `discarded_j` by n_j/N
   fails. SPA-3's "dense − kept = discarded per beamlet" is unchanged (per-beamlet quantities need no weight).
4. **INF-4 at both thresholds, every field (finding 4).** At τ = 0 and τ = 1e-2 every `ColumnSet` field of the named
   beamlets is bitwise equal to the full run's (INF-4 row).
5. **Minor (finding 5).** Decision 0046 §8: the `Simulation.run(beamlets=subset)` sentence is marked "[superseded by
   A-3: `Simulation.influence_columns`]"; the ledger says "INF-4 (diagnostic column API)" instead of "subset runs".
6. **Minor (finding 6).** SPA-3: "τ = 0 reports exactly 0 for every column with positive in-grid energy; a transported
   zero-deposition column reports undefined (NaN + mask false) at every τ" — one sentence; no unconditional "exactly
   0" remains in the row.


**Amendment A-5 (2026-11-25, S0, before any code; Codex REVIEW-54588747 on a783ac8a, changes_required; three
important findings and one minor, all adopted).** The frozen text is edited in place as for A-1..A-4.

1. **The memory rule is an accumulator-payload budget (finding 1).** §B.3: `memory_budget_bytes` bounds the
   accumulator payload only (grid accumulators, tally and counter rows of one chunk, influence accumulators plus one
   host copy); the plan says explicitly what it does NOT bound (diagnostic arrays, retained CSC/batch output and the
   assembled `Result`, Python object overhead, the Warp kernel/module cache); automatic G is a payload rule, not an
   RSS guarantee; S7 documents this. MEM-1 measures exactly this payload class, with the estimate computed for the
   current group's range.
2. **`influence_columns` fails unconditionally on an invalid transport (finding 2).** §B.10: a `ColumnSet` carries
   the transport counters and is returned only when every invalidating counter is zero, else `TransportLimitError`
   (or the existing invalidating error); `allow_invalid_result` is a refused field. INF-4 gains the discriminating
   control (`max_steps` patched to truncate beamlet 1: `influence_columns([1])` raises; `run()` with
   `allow_invalid_result` yields an invalid `Result`).
3. **Label grammar (finding 3).** §B.1: `^[A-Za-z0-9_.-]{1,64}$`, unique; the same grammar for the spot-file label
   column; refused by `validate` (`UnsupportedCombinationError`) and the spot parser (`ConfigError`). Rows BEAM-1,
   CLI-2, CLI-3, PER-4 gain the controls; the report test asserts no `<`, `[`, `](` or `|` from labels.
4. **Thresholded controls must discard something (minor finding 4).** INF-4 (τ = 1e-2) and the PER-4 unequal-n_j
   control assert `discarded_count_j ≥ 1` and `discarded_energy_j > 0` for EVERY beamlet j (R-INF, n_j ∈ {2 000,
   4 000, 2 000, 4 000}, σ = 3 mm; TO VERIFY at S4: if a column has nothing below τ the test fails and the
   configuration is amended prospectively, never the assertion).

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
  `chunk_histories` (INF-3); (b) a J = 1 set equals the pencil run (INF-1); (c) `influence_columns` reproduces its
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
| S1 | `Beamlet`, `BeamletSet`, `InfluenceOptions`, `validate` rules, `ionmc.config.allocation_estimate` (§B.3, A-4 item 2), source-consumer audit (every `config.source` reader uses the beamlet envelope or per-beamlet values), summary/reconcile rows | BEAM-1..BEAM-3, CAP-1 (CI part), CAP-2..CAP-4, REG-2, REG-3 |
| S2 | reference backend: per-history beamlet id, per-beamlet accumulators, per-beamlet balance | BEAM-4..BEAM-6, INF-1, INF-2, INF-6, INF-7, MEM-1 (python) |
| S3 | Warp driver: device context, offsets/parameter arrays, grouped accumulators, host sparsification | INF-1..INF-5, INF-7, MEM-1 (warp-cpu), REG-1 |
| S4 | `ionmc.influence`: CSC container, threshold accounting, helpers, per-beamlet statistics | SPA-1..SPA-5, UNC-1, UNC-4 |
| S5 | schema 1.1: writer/reader/tier-1 bindings/tier-2 re-derivation, golden 1.1 fixture, v3-006 OUT-4 control moved | PER-1, PER-3..PER-10, OUT-0-INF, REG-5 |
| S6 | CLI/TOML/spot file, report rows and SVGs, `Simulation` reuse and the `influence_columns` diagnostic API | CLI-1..CLI-3, REUSE-1 (CI part), INF-4 |
| S7 | docs (results format 1.1, configuration, CLI, workflow, architecture/transport beamlet section), ledger notes incl. SHOULD coverage section | DOC-1, docs build `--strict` |
| S8 | LV/HR evidence driver `validation/scripts/influence/inf7_run.py` (pattern of `out6_run.py`), one heavy job at a time | PER-2, PAR-1..PAR-3, UNC-2, UNC-3, CLI-4, CLI-5, CAP-1-LV/HR, MEM-1 (LV warp-cpu, HR warp-cuda), REG-1-HR, PERF-1, PERF-2 |

### Discriminating rows (defined before implementation; all gating unless marked *measurement*)

Reference configuration **R-INF** (used where a row says "R-INF"): water box 60×60×200 mm, influence grid 2 mm
(30×30×100), protons, 4 beamlets at (−10, 0), (+10, 0), (0, −10), (0, +10) mm lateral, energies 70 / 90 / 110 /
130 MeV, σ_E = 0.5 MeV, σ_lat = 3 mm, direction +z, B = 10, n_j = 2 000 (CI: 200, B = 4). Default physics,
`bethe-ionmc-v1`, nuclear off unless stated. Fifth variant **R-INF-Z** (A-2 item 4): R-INF plus one beamlet aimed
outside the influence grid (zero-deposition column), used by SPA-3, PER-1, PER-4, CLI-3.

**BEAM — beamlet model and mapping**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| BEAM-1 | Config validation: each fail-closed rule of §B.1 and §B.7 provoked separately | `UnsupportedCombinationError` (or `ConfigError` for TOML) before transport, message names the field; each positive twin passes; includes (A-1) `n_histories` ≠ Σ n_j refused (message names both values), absent TOML value derived, equal value accepted; label of 65 characters / duplicate refused / (A-5 item 3) a label with whitespace, `<`, `[`, `(`, `|`, a control character or non-ASCII refused (grammar `^[A-Za-z0-9_.-]{1,64}$`); the converted ion-mapping failure branch (`config.py` 498) provoked with a BeamletSet and naming the beamlet; (A-3 item 2) explicit `group_size` 0 / J + 1 / over-budget refused; automatic G = 1 over budget refused; a budget that admits exactly G = 2 yields effective group_size 2 (`derived: memory-budget`), the budget computed from `allocation_estimate` (A-4 item 2) | CI |
| BEAM-2 | Source-consumer audit: a test enumerates every reader of `config.source` (static list checked against `grep` of `src/ionmc`, reviewed; FACT at 3f3736a6, complete after A-1 item 2: `config.py` (348; 471–498 ion-mapping path, its failure branch at 498 converted to per-beamlet values; 526; 564; 588; 1016; max_steps E_hi 755–756; channel compile 807; nuclear bounds 921–943; the decision-0044 stopping-model window check 1035–1064; ion domain 1127–1153; sigma rules and ion mapping 1237–1294 incl. `fdm_domain_status(proj, mats, e0, σ)` at 1287; stopping identity species 594–611; summary 624; `unit_direction=src.unit_direction()` 1351), `reference.py` (projectile 193–260, sampling 840–860), `warp_driver.py` (ctl source fields 121–135, `ctl.dir` at 135), `kernels.py` (source sampling), `results.py` (producible 403; profile axis from `unit_direction` 725 and 1626; binding 1668–1675; stopping-window check on E0 + 6σ 1806; tally E_hi bound 2093), `_result_shapes.py` 335–345, `contract.py` 78–87, 308–320, 701 (`SOURCE_TYPES`), 939, `simulation.py` 736 (`capabilities()["sources"]`), `report.py` 335, `configio.py` 590; A-2 item 2 additions `config.py` 844, 1163, 1214–1229, `warp_driver.py` 187; `configio.py` 1162 is a lookup-file source, not `SimulationConfig.source`, listed as `not-a-source`) and asserts for a beamlet set the envelope rules: `max_steps`, channel `E_hi`, nuclear energy bounds, ion mapping domain and energy-sigma rules use min/max over beamlets | every consumer covered; a beamlet with E above the domain is refused although the mean energy is inside; exhaustive AST audit (A-2 item 2, replaces the grep rule): a test parses every module of `src/ionmc` with `ast`, follows attribute accesses rooted at `.source` of a config / effective-config value AND every alias binding in the same scope (`src = cfg.source`, `source = eff.source`; tuple/keyword passing is recorded as a reader at the call site) and asserts every (module, line, attribute) is in the static allowlist annotated `envelope` / `per-beamlet` / `refused-in-S1` / `not-a-source` (the allowlist is the inventory); any access outside it fails; the allowlist is re-derived (and the decision's Context updated) whenever S2/S3 add readers; exceptional/failure paths included | CI |
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
| INF-4 | Column independence (A-3 item 1, A-4 items 1 and 4): `influence_columns([2])` and `influence_columns([0, 3])`, each at τ = 0 AND τ = 1e-2 | at both thresholds, for the beamlets named, every `ColumnSet` field — column voxel indices, energy/dose values, std, batch values, nnz, discarded energy / fraction (incl. the undefined mask) / count / max dose / variance, and the per-beamlet balance — is bitwise equal to the full run's; invalid subsets (empty, duplicate, out of range, non-integer) refused (`UnsupportedCombinationError`); `save(ColumnSet)` refused (`TypeError`); **capability controls (A-4 item 1):** pencil source, no `InfluenceOptions` and each diagnostic field (traces, end positions, escape records, every other diagnostic/tally field enumerated from the dataclass) → refused (`UnsupportedCombinationError`); the allowlist test fails when a config field is neither honoured nor refused; (A-5 item 2) with `max_steps` patched to force step truncation in beamlet 1, `influence_columns([1])` raises and returns nothing, while the same configuration through `run()` with `allow_invalid_result` yields an invalid `Result` (existing behaviour), and `allow_invalid_result` set in the configuration is refused by `influence_columns`; (A-5 item 4) at τ = 1e-2, `discarded_count_j ≥ 1` and `discarded_energy_j > 0` for EVERY beamlet j of the configuration (R-INF, n_j ∈ {2 000, 4 000, 2 000, 4 000}, σ = 3 mm; TO VERIFY at S4, never weaken the assertion) | CI |
| INF-5 | Trajectory parity python vs warp-cpu float64 for the first 3 histories of beamlets 0, 1 and J−1 (across beamlet boundaries) | decision-0037 trajectory class (existing T1 tolerances) | CI |
| INF-6 | **Negative controls** (patched in the test): (i) `c0` reset per beamlet → BEAM-4 and INF-4 fail; (ii) beamlet id off by one at a boundary, or an accumulator row error with correct source parameters → INF-2(a) still passes; INF-7 fails (BEAM-5 and INF-1 catch only a wrong source assignment; shows INF-2 alone is not sufficient); (iii) accumulator row stride B±1 → INF-2(a) fails | each detected by the named row | CI |
| INF-7 | **Per-beamlet closure**: Σ_v column_j quanta·Q + outside_j + quantization_j = deposited_j, reduced from the tally rows of h ∈ [s_j, s_j + n_j), for every j; negative controls "accumulator row off by one" and "columns permuted" must FAIL this row | exact (rel 1e-12) | CI, LV, HR |

**SPA — sparse representation and threshold**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| SPA-1 | τ = 0 round trip dense ↔ CSC (`to_dense`, `from_dense`) | bitwise; nnz = count of nonzero integer accumulators summed over batches | CI |
| SPA-2 | Canonical form: indptr monotone, voxel strictly increasing per column, dtypes, no explicit zeros at τ = 0 | exact | CI |
| SPA-3 | Threshold accounting, τ ∈ {1e-4, 1e-2, 0.1}: every kept entry ≥ τ·max_j, every dropped < τ·max_j; dense − kept = discarded per beamlet (energy, rel 1e-12); reported fraction = discarded/in-grid energy; counts exact; τ = 0 reports exactly 0 for every column with positive in-grid energy; a transported zero-deposition column reports undefined (NaN + mask false) at every τ, with counts 0, which R-INF-Z exercises at τ = 0 AND τ = 0.1 (A-2 item 4, A-3 item 4, A-4 item 6); per-beamlet quantities are per primary of beamlet j and need no n_j/N weight | exact | CI, LV |
| SPA-4 | Helpers: `matvec(w)` vs `to_dense() @ w`, `rmatvec(y)` vs `to_dense().T @ y`, `plan_dose(weights)` | relative 1e-12 | CI |
| SPA-5 | Threshold is never silent: τ appears in requested/effective/report; changing τ changes `config_sha256`; a bundle whose discarded arrays are removed fails `load` | exact | CI |

**UNC — per-beamlet uncertainty**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| UNC-1 | Entry std = batch standard error recomputed from the persisted `[B, nnz]` batch values; masks per decision 0040 rules | relative 1e-9 (reader tolerance of 0045 A-1) | CI |
| UNC-2 | *statistical*, **coverage between independent replicas** (r = 4, k = 0/1 compared; HR r = 5): R-INF variant with J = 16 beamlets (4 energies × 4 positions), n_j = 10 000, B = 20. Entry set S_j is fixed by an independent pilot replica (r = 4, k = 2; HR r = 5, k = 2; beamlets without entries are simply excluded, A-2 item 4): defined, n_nonzero = B and relative uncertainty ≤ 0.05 in the pilot; the compared replicas are k = 0 and k = 1. z = Δ/√(s₁² + s₂²); per-beamlet coverage c_j(1) = fraction |z| ≤ 1, c_j(2) = fraction |z| ≤ 2. Reference p(1), p(2) = P(\|T_ν\| ≤ 1, 2) with ν = 2(B − 1) = 38 (computed by numpy quadrature in the test) | for both levels L ∈ {1, 2}: \|mean_j c_j(L) − p(L)\| ≤ max(3·sd_j(c_j(L))/√J, 0.03) (beamlets independent; 0.03 = declared margin for non-normal low-count batches, fixed now). **Power control:** recomputing with s scaled by 0.8 must FAIL the same criterion; if it passes, the row is INCONCLUSIVE (gating fail) | LV, HR |
| UNC-3 | *statistical*, history scaling (r = 6): R-INF J = 4, n_j = 4 000 (k = 0) vs 16 000 (k = 1) (B = 20); entry set S fixed by an independent pilot replica (r = 6, k = 2; same criteria as UNC-2) | median over S of s(16 000)/s(4 000) ∈ [0.45, 0.55] (expected 0.5) | LV |
| UNC-4 | `plan_uncertainty(x)` = √(Σ_j x_j² σ_vj²) exactly (rel 1e-12); docstring and results-format page state the independence argument (disjoint counters, decision 0037) | exact | CI |

**PER — persisted influence bundles (schema 1.1)**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| OUT-0-INF | Clause map: geometry, beamlet mapping (offsets, labels, parameters), units, histories per beamlet, normalisation, threshold, discarded accounting, uncertainty definition → JSON pointer/array | every clause resolved | CI |
| PER-1 | In-process round trip `load(save(r))` for R-INF (python, warp-cpu f64), τ = 0 and τ = 1e-2, `persist_batches` on/off (off recorded `not-persisted`); includes a bundle whose `requested.run` carries the derived N = Σ n_j (A-1 item 1) and an R-INF-Z bundle with a zero-deposition column at τ = 0 and τ = 0.1 (A-2 item 4, A-3 item 4) | arrays bitwise, JSON equal | CI |
| PER-2 | Clean-process round trip (`python -I`) | canonical dump equal | CI; LV/HR on three backends |
| PER-3 | Version handling: the frozen V3-006 1.0 golden fixture loads unchanged; a 1.1 golden synthetic fixture (KB-scale) loads and regenerates its report byte-identically; "1.2"/"2.0" refused; `kind = "influence"` with schema 1.0 refused | exact | CI |
| PER-4 | Tier-1 negative controls, one each: indptr non-monotone or ≠ nnz; voxel index ≥ V or unsorted or duplicate; nnz arrays of unequal length; negative std; batch values not reproducing mean/std; discarded accounting inconsistent; beamlet table ≠ `requested`; offsets ≠ cumulative n_j; columns swapped with balances unchanged and digests recomputed (fails through the INF-7 identity); `Σ_j (n_j/N)·(kept_j + discarded_j)` ≠ composite in-grid energy (rel 1e-9; N = Σ_j n_j; `kept_j` = Σ_v column_j energy, `discarded_j` the per-beamlet discarded energy; reader control at τ = 1e-2 with unequal n_j, R-INF with n_j ∈ {2 000, 4 000, 2 000, 4 000}, so an implementation that does not weight `discarded_j` by n_j/N fails; A-4 item 3); influence block missing for `kind = "influence"`; a label violating the §B.1 grammar (e.g. `<img src=x>` edited into `result.json`; A-5 item 3); for the τ = 1e-2 unequal-n_j control, `discarded_count_j ≥ 1` and `discarded_energy_j > 0` are asserted for EVERY beamlet j (A-5 item 4); zero column with nonzero discarded count or a defined fraction (A-2 item 4; R-INF-Z at τ = 0 and τ = 0.1, A-3 item 4) | each `ResultFormatError` naming the defect | CI |
| PER-5 | Budgets: a catalogue declaring nnz = 2⁴⁰ or a 3-D batch array beyond `max_array_bytes` refused before allocation (V3-006 OUT-10 method, RSS bound 64 MiB) | exact | CI |
| PER-6 | Reproduction from the bundle (r = 11): `from_dict(requested)` + seed → every content digest equal (full run only; k = 1 released, marked unused, A-3 item 1) | exact | CI (small), LV |
| PER-7 | Tier 2 `verify_effective` re-derives `effective` for a beamlet bundle | exact | CI (default model), LV |
| PER-8 | No table values (allowlist extended only by `influence__*`, `beamlet__*`) | exact | CI |
| PER-9 | Honoured check: influence requested ⇒ influence arrays present; absent request ⇒ none | exact | CI |
| PER-10 | Bundle at the J limit (J = 4 096, n_j = B = 2, warp-cpu) saves and loads under the default budgets; `requested` + `beamlets` stay within `max_json_bytes` and the token budget; the size is computed with 64-character labels at J = 4 096 and asserted within `max_json_bytes` (A-1 item 3) | exact | CI |

**CAP — capability contract**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| CAP-1 | V3-006 sweep (CAP-1a/1b dimensions) × source kind {pencil, beamlets} × influence {off, on} | declared == enforced; zero third outcomes | CI; LV (tabulated); HR (CUDA) |
| CAP-2 | Declaration derived (monkeypatched constant changes the `beamlets` entries of `capabilities()["sources"]`; no new top-level key) | exact | CI |
| CAP-3 | Reconciliation table extended (every new field equal/derived/unused/requested-only; `group_size` None → derived `memory-budget`; `run.n_histories` absent in TOML → derived `sum-of-beamlets`; no `requested-only` row exists for subsets, A-3 item 1); zero unlisted differences | exact | CI |
| CAP-4 | Dispatch: patching the driver to ignore beamlet ids or the threshold makes `Simulation`/`save` raise | exact | CI |
| CAP-5 | Exit codes through `ionmc check/run` for every §B.7 refusal (4), TOML errors (2) | exact | CI |
| CAP-6 | warp backends + beamlets + nuclear → exit 4 until V3-005B, message names V3-005B's absence; python + nuclear + beamlets accepted | exact | CI |

**PAR — backend parity of influence columns**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| PAR-1 | python f64 (r = 0, n_j = 2 000) vs warp-cpu f32 (r = 2): per-beamlet in-grid energy, z = Δ/√(σ₁²+σ₂²) | \|z\| ≤ 4 for every beamlet (J = 4) | LV |
| PAR-2 | warp-cpu f32 vs warp-cuda f32 (r = 2/3): per-beamlet in-grid energy \|z\| ≤ 4 and entry-wise coverage as UNC-2 | both | HR |
| PAR-3 | warp-cpu f32 vs warp-cuda f32, same seed: per-beamlet integer accumulators | report-only (not gating; float32 device math differs); recorded | HR |

**MEM — allocation audit (A-4 item 2)**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| MEM-1 | Instrumented allocation audit: the accumulator-payload class of §B.3 (A-5 item 1: grid accumulators, tally and counter rows of one chunk, influence accumulators and one host copy; NOT diagnostics, retained output, Python overhead or the kernel cache) is measured as the peak over chunks of the sum of the byte sizes of its arrays, with the estimate computed for the current group's range; every such array (numpy `ndarray.nbytes` of arrays created by the driver and `wp.array` sizes, recorded through a test hook on the allocation helpers) is instrumented and its peak sum is at most `allocation_estimate(...)` (§B.3); for the python path and warp-cpu the peak is measured per chunk; the estimate must not exceed the measured peak by more than a factor 1.5 (an estimate that is too loose is also a fail, recorded) | measured peak ≤ estimate ≤ 1.5 × measured peak, for each of the python, warp-cpu f64 and f32 variants, with and without ions/nuclear where the backend supports them, and for G ∈ {1, 2, J} | CI (python, warp-cpu); LV (warp-cpu); HR (warp-cuda) |

**CLI, reuse, performance, regression, docs**

| Row | What | Pass iff | Tier |
|---|---|---|---|
| CLI-1 | `ionmc run beamlets.toml` (inline) and with a spot file → bundle; `verify` 0; `report --check` 0 | exact | CI |
| CLI-2 | Spot file: sha256 mismatch → 2; malformed header/row → 2 naming the line; recorded in `InputManifest` role `spots`; a spot-file label violating the §B.1 grammar (whitespace, `<`, `[`, `(`, `|`, non-ASCII) → `ConfigError`, exit 2 (A-5 item 3); a file above 1 MiB (sparse file, nothing read) and one above 4 097 lines refused with exit 2 before parsing, naming the bound; the sha256 is computed over the bytes read once (A-1 item 3); adversarial controls (A-2 item 3): symlink to a valid file → refused; oversized sparse regular file → refused by `fstat` before any read; file replaced at its name between open and read (the test patches the open helper to rename after opening) → bytes and sha256 are those of the opened descriptor; file grown after the first `fstat` → refused by the second `fstat` (the prefix is never accepted; A-3 item 3) | exact | CI |
| CLI-3 | Report content: beamlet table, threshold line, discarded summary, the strings "undefined (no deposition in the influence grid)" and "undefined (no entries)" for an R-INF-Z zero column at τ = 0 and τ = 0.1 (A-2 item 4, A-3 item 4), "Not a clinical claim.", no `\bvalidated\b` in generated text; the generated Markdown contains no `<`, `[`, `](` or `|` originating from labels, and a report regenerated from a bundle whose `result.json` label was edited to `<img src=x>` fails `load` (tier 1 re-validates the grammar), so the report is never produced (A-5 item 3) | exact | CI |
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
- 2026-11-25 (S0 fix 2, Codex REVIEW-15a83983 on 0901763, changes_required; four important findings adopted, see
  Amendment A-2 after A-1): (1) configured vs. transported history counts, subset validation; (2) AST-based
  source-consumer audit; (3) descriptor-bound spot-file read; (4) zero-deposition columns and R-INF-Z. Rows BEAM-2,
  CLI-2, CLI-3, INF-4, PER-1, PER-4, PER-6, SPA-3, UNC-2, UNC-3 updated. Not a clinical claim.
- 2026-11-25 (S0 fix 3, Codex REVIEW-1d19386a on dd5ad774, changes_required; six important findings adopted, see
  Amendment A-3 after A-2): (1) subset runs become the diagnostic `influence_columns` API (no `Result`, no bundle);
  (2) fail-closed influence grouping; (3) second `fstat` on the spot file; (4) τ = 0 and zero-deposition columns.
  Rows BEAM-1, INF-4, PER-1, PER-4, PER-6, CAP-3, SPA-3, CLI-2, CLI-3 updated. Not a clinical claim.
- 2026-11-25 (S0 fix 4, Codex REVIEW-a52098eb on af4aebc6, changes_required; four important and two minor findings
  adopted, see Amendment A-4 after A-3): (1) fully fail-closed `influence_columns`; (2) exact backend-aware allocation
  formula and the new row MEM-1; (3) thresholded composite identity with the n_j/N weight and the unequal-n_j reader
  control; (4) INF-4 at both thresholds, every field; (5) superseded-sentence marking in decision 0046 and the ledger
  wording; (6) SPA-3 τ = 0 sentence. Rows BEAM-1, INF-4, MEM-1, PER-4, SPA-3 updated. Not a clinical claim.
- 2026-11-25 (S0 fix 5, Codex REVIEW-54588747 on a783ac8a, changes_required; three important findings and one minor
  adopted, see Amendment A-5 after A-4): (1) the memory rule is an accumulator-payload budget with an explicit "does
  not bound" list; (2) `influence_columns` raises unconditionally on an invalid transport; (3) label grammar; (4)
  thresholded controls assert nonzero discards. Rows BEAM-1, CLI-2, CLI-3, INF-4, MEM-1, PER-4 updated. Not a
  clinical claim.
