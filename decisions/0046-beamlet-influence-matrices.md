# 0046 — Beamlet-resolved transport, sparse dose influence matrices and persisted influence bundles

- Status: Accepted (frozen with plan `validation/plans/v3-007-acceptance.md` at the S0 commit, 2026-11-25). Later
  changes are dated amendments. Not a clinical claim.
- Amended 2026-11-25 (S0 fix 6, Codex REVIEW-0d4bb5f6 on 8db925e7): (1) one consistent accumulator-payload
  contract (§2); (2) `Beamlet.weight` contract (§1).
- Amended 2026-11-25 (S0 fix 5, Codex REVIEW-54588747 on a783ac8a): (1) the memory rule is an accumulator-payload
  budget (not a process-RSS guarantee); (2) `influence_columns` raises unconditionally on an invalid transport;
  (3) label grammar `^[A-Za-z0-9_.-]{1,64}$`. See plan Amendment A-5.
- Amended 2026-11-25 (S0 fix 4, Codex REVIEW-a52098eb on af4aebc6): (1) `influence_columns` fully fail-closed; (2) exact
  backend-aware allocation formula and audit row MEM-1; (3) thresholded composite identity weighted by n_j/N;
  (4) INF-4 at both thresholds; minor: superseded-sentence marking, SPA-3 τ = 0 sentence. See plan Amendment A-4.
- Amended 2026-11-25 (S0 fix 3, Codex REVIEW-1d19386a on dd5ad774): (1) subset runs become the diagnostic
  `Simulation.influence_columns` API (no `Result`, no bundle); (2) fail-closed influence grouping; (3) second `fstat`
  on the spot file; (4) τ = 0 and zero-deposition columns. See plan Amendment A-3.
- Amended 2026-11-25 (S0 fix 1, Codex REVIEW-1652a7c1 on 9805dfdb): (1) `n_histories` = Σ n_j for a `BeamletSet`;
  (2) reader inventory completed; (3) label and spot-file bounds; (4) subset-run semantics. See plan Amendment A-1.
- Amended 2026-11-25 (S0 fix 2, Codex REVIEW-15a83983 on 0901763): (1) configured `n_histories` = Σ_all j n_j versus
  the persisted `run.histories_transported` [superseded by A-3];
  (3) descriptor-bound spot-file read; (4) zero-deposition columns. See
  plan Amendment A-2.

- Limits: 1 ≤ J ≤ 4 096 beamlets (ESTIMATE-justified: ≈ 300 B per beamlet inline in `requested` ≈ 1.2 MiB, within the
  0045 `max_json_bytes` and token budgets; the frozen workload has 100). Envelope rules: quanta capacity
  `e_cap = max_j (E_j + 6σ_j)`; `max_steps` from the envelope E_hi; `fdm_domain_status` and the σ ≤ fraction·E0 rule per
  beamlet.
- Task: V3-007. Base: `v3/develop` 3f3736a6 (V3-006, decision 0045).
- Affects: architecture (source model, execution scheduling, memory layout), public API (`BeamletSet`,
  `InfluenceOptions`, `ionmc.influence`), persisted format (`ionmc.result` 1.1), capability contract, performance.
  **No physics, table, RNG encoding, stopping-model or evidence change.**
- Builds on: 0037 (shared physics, Philox encoding: "beamlet membership is a range of history indices", `c0` never
  reset; int64 fixed-point accumulators; numpy + warp only), 0040 (batch reduction, masks, undefined = NaN), 0041
  (nuclear python only), 0044 (no silent switch), 0045 (bundle, versioning MAJOR.MINOR, `kind = "influence"` reserved,
  reconciliation table, tiers 1/2, exit codes, O3 batch persistence, path set P).

## Problem

V1-MUST-009/010 require scoring attributable to the originating beamlet and efficient construction of dose influence
matrices; SHOULDs ask for sparse representations, batched beamlet execution without re-initialising immutable data,
planning-oriented reuse and planning-aware uncertainty; V2-OUTPUT requires influence matrices to retain geometry,
beamlet mapping, units, history and normalisation metadata; V2-CAP forbids influence flags becoming silent no-ops.

## Context (FACT at 3f3736a6)

- `SimulationConfig.source` is one `PencilBeamSource` (`config.py:306–319`, `sources.py`). Complete list of its
  readers at 3f3736a6: `config.py` (max_steps E_hi 755–756; channel compile 807; nuclear bounds 921–943; the
  decision-0044 stopping-model window check 1035–1064; ion domain 1127–1153; sigma rules and ion mapping 1237–1294
  incl. `fdm_domain_status(proj, mats, e0, σ)` at 1287; stopping identity species 594–611; summary 624;
  `unit_direction=src.unit_direction()` 1351), `reference.py` (projectile 193–260, sampling 840–860),
  `warp_driver.py` (ctl source fields 121–135, `ctl.dir` at 135), `kernels.py` (source sampling), `results.py`
  (producible 403; profile axis from `unit_direction` 725 and 1626; binding 1668–1675; stopping-window check on
  E0 + 6σ 1806; tally E_hi bound 2093), `_result_shapes.py` 335–345, `contract.py` 78–87, 308–320, 701
  (`SOURCE_TYPES`), 939, `simulation.py` 736 (`capabilities()["sources"]`), `report.py` 335, `configio.py` 590 and
  1162. Amended list (REVIEW-1652a7c1): `config.py` also 348, 471–498 (ion-mapping path; its failure branch at 498
  reads the singular `source.kinetic_energy_mev` and is converted to per-beamlet values), 526, 564, 588 and 1016.
  Amended again (REVIEW-15a83983): `config.py` also 844, 1163, 1214–1229 and `warp_driver.py` 187; `configio.py` 1162
  is a lookup-file source, not `SimulationConfig.source` (kept in the audit allowlist as `not-a-source`); the
  inventory is the allowlist of the AST-based audit of plan row BEAM-2 and is re-derived whenever S2/S3 add readers.
  Every reader, including exceptional paths, is in the list; the plan's grep check finds none outside it.
- `run_warp_range(eff, h0, h1, device)` rebuilds kernel handle, uploads tables, geometry and grids, and allocates
  `edep int64[B, ΣV]` on every call; source parameters are scalars in the `ctl` struct; launches chunks of
  `chunk_histories`.
- `validate` requires 1 ≤ N ≤ 2³², B ≥ 2, N % B == 0, accumulators `B·ΣV·8·workers ≤ memory_budget_bytes`
  (default 2³¹) and the per-batch quanta capacity check (`per_batch·e_cap/QUANTUM < MAX_QUANTA`); for beamlets
  `e_cap = max_j (E_j + 6σ_j)`.
- Batch membership `h mod B`; per-history tally rows are reduced exactly on the host.
- Warp + nuclear raises (V3-005B unmerged, deferred by the operator).
- `experiment/v3/performance.json` workload `influence`: 100 beamlets × 10 000 primaries, multiple energies/
  positions, proton-3d planning volume, "same complete proton physics" (incl. nuclear and transported secondaries),
  including source generation, transport, beamlet uncertainty, assembly and persisted output.
- ESTIMATE: influence grid of that workload at 2 mm = 60×60×150 = 540 000 voxels; dense per-beamlet batch
  accumulators for all 100 beamlets at B = 20 would need 100·20·540 000·8 B ≈ 8.6 GB.

## Candidate approaches and selection

### 1. What a beamlet is

| Option | Assessment |
|---|---|
| (a) **`BeamletSet(projectile, beamlets)`; `Beamlet` = position, direction, E, σ_E, σ_lat, n_histories, optional weight and label** (selected) | reuses the validated pencil-beam sampling per beamlet; arbitrary directions (no axis privileged); one projectile keeps one table set and one kernel variant |
| (b) per-beamlet full `SimulationConfig`s | duplicates physics/geometry per beamlet; invites inconsistent physics between columns |
| (c) histories sampled proportionally to plan weights | uneven column uncertainty; weights change during optimisation, so the matrix would have to be recomputed — rejected; `weight` is metadata only |
| (d) mixed projectiles in one set | two table sets/kernels per run; no planning need in M1 — rejected (fail closed) |

Normalisation: columns in MeV and Gy **per primary of beamlet j**; plan dose `D = A·x` with x in primaries.
`Beamlet.weight` (A-6 item 2): `None` or a finite float ≥ 0 (primaries in the plan), persisted in the beamlet
table and `requested`; `plan_dose`/`plan_uncertainty` use the stored weights when no vector is given and refuse
when any is `None`; an explicit vector (length J, finite, ≥ 0) overrides; weights never change sampling, columns or
statistics (plan row WGT-1).

Amendment (REVIEW-1652a7c1): for a `BeamletSet`, `RunOptions.n_histories` MUST equal Σ_j n_j; `validate` refuses a
mismatch, and TOML `[run] n_histories` is optional (absent → derived Σ n_j, different → `ConfigError`); this is
the CONFIGURED count and the only history count (A-3: every bundle is a full run; the A-2 transported-count field
`run.histories_transported` is
[superseded by A-3]). `label` is optional, matches `^[A-Za-z0-9_.-]{1,64}$` (A-5; no whitespace or
Markdown/HTML metacharacters, also for the spot-file label column) and is unique within the set. A spot file is
opened with
`O_NOFOLLOW | O_RDONLY | O_CLOEXEC` (a symbolic link is refused), `fstat` on the descriptor must report a regular
file ≤ 1 MiB, exactly `st_size` bytes are read from that descriptor, a second `fstat` on the same descriptor must
report the same size and inode (else `ConfigError` "spot file changed while being read", A-3), sha256 over those
bytes, ≤ 4 097 lines (violation → `ConfigError`).

### 2. Execution and batching

| Option | Assessment |
|---|---|
| (a) sequential beamlets, one `run_warp_range` per beamlet | re-uploads tables/geometry per beamlet (violates the batched-execution SHOULD); a 10 000-history launch fills only a small fraction of the A6000's resident threads (ESTIMATE ≈ 8 %), so throughput collapses |
| (b) **one device context per run; per-history beamlet id from the global history index; source parameters in per-beamlet arrays; accumulators `int64[G·B, V]` for a group of G beamlets; launch chunks may span beamlets** (selected) | full occupancy; memory bounded by G; per-beamlet results bit-identical for any G and chunk size because integer accumulation is associative and each (beamlet, batch, voxel) cell receives the same integer increments; the pencil source is J = 1 of the same path (bit identity guarded by A16/T1) |
| (c) all beamlets at once with dense `[J·B, V]` | 8.6 GB for the workload; does not scale — it is (b) with G = J when it fits |
| (d) on-device sparse/hash accumulation | non-deterministic insertion order, atomics on hash slots, high complexity; reconsider in V3-014 only with measurements |
| (e) per-history (history-by-history) variance instead of batches | needs per-history scratch per voxel or a second moment accumulator with float rounding; departs from decision 0040 — rejected |

The kernel change is confined to source sampling (parameters read from arrays indexed by the beamlet id) and the
accumulator row (`(j − j0)·B + batch`). The beamlet id is `searchsorted(offsets, h, right) − 1` on a device array
(J ≤ 4 096: ≤ 12 comparisons); the id is range-asserted (0 ≤ j < J) in every backend and can never wrap. G = largest
group in [1, J] fitting `memory_budget_bytes` as an ACCUMULATOR-PAYLOAD budget (A-3; exact since A-4 item 2; narrowed
by A-5 item 1): the rule is `allocation_estimate(backend, B, ΣV, chunk, V_inf, G) ≤ memory_budget_bytes`,
implemented in S1 and bounding exactly the accumulator-payload class (A-6 item 1): (i) grid accumulators
`8·B·S`, (ii) influence accumulators `8·G·B·V_inf`, (iii) the tally rows (width `N_FIXED_TALLIES + 2·n_g +
channel_columns`, + `len(NUCLEAR_TALLY_NAMES)` on python with nuclear) and counter rows (width `len(COUNTER_NAMES)` +
`len(ION_COUNTER_NAMES)` for ions, + `len(NUCLEAR_COUNTER_NAMES)` on python with nuclear) a backend holds at one time
(one chunk on warp, the whole range N on python), (iv) on warp-cuda one host copy of the current group's influence
accumulators during transfer; chunk = min(`chunk_histories`, longest group range that G implies). The formulas are
in plan §B.3. Not in the class and not bounded: geometry/tables/descriptors/parameter arrays, end-state arrays and
end-phase `edep` copies, reduction temporaries, diagnostic arrays, the retained CSC/batch output of earlier groups
and the assembled `Result` (bounded at `save` by the 0045 budgets), Python object overhead, the Warp kernel/module
cache; automatic G is a payload rule, not an RSS guarantee. MEM-1 instruments exactly the class (i)–(iv) and
compares its peak with the estimate of the same configuration (measured ≤ estimate ≤ 1.5 × measured).
 the audit row MEM-1 measures the allocations and fails
an estimate that is below the measured peak or more than a factor 1.5 above it. Pencil configurations keep the
`validate` rule of 3f3736a6; G = 1 not fitting, an explicit G outside [1, J] or over the
bound is refused by `validate` before allocation (`UnsupportedCombinationError`); recorded as effective
(`derived: memory-budget`).

### 3. Influence matrix representation

| Option | Assessment |
|---|---|
| (a) dense `[V, J]` | 432 MB float64 per quantity for the workload, mostly zeros |
| (b) COO unsorted | not canonical; content digest depends on order |
| (c) **canonical CSC by beamlet (indptr, voxel, values), voxel strictly increasing per column** (selected) | assembly is column-wise by construction; canonical ⇒ deterministic digests; `A·x` and `Aᵀ·y` in numpy (`np.add.at`/`bincount`); users with scipy wrap the arrays with no copy (documented recipe; scipy is not a dependency) |
| (d) CSR by voxel | requires a global transpose at assembly; no advantage for column-wise production |

### 4. Thresholding

Selected: default `threshold_rel = 0` (lossless: every voxel with a nonzero integer sum). Optional τ ∈ (0, 1) per
beamlet relative to the column maximum dose; requested and effective, in `config_sha256`, reported; per-beamlet
discarded energy, fraction, entry count and maximum discarded dose persisted and bound by the reader (dense identity
`Σ_j (n_j/N)(kept_j + discarded_j) = composite in-grid energy` (rel 1e-9; `kept_j` = Σ_v column_j energy, discarded
energy weighted by n_j/N like the kept energy, A-4 item 3. Reader control: τ = 1e-2 with unequal n_j (R-INF,
n_j ∈ {2 000, 4 000, 2 000, 4 000}), so an implementation that does not weight `discarded_j` by n_j/N fails).
Rejected: absolute thresholds (unit- and
normalisation-dependent), thresholds on relative uncertainty (selects on noise ⇒ bias), silent default thresholds.

### 5. Per-beamlet uncertainty

Batch standard error per entry (decision 0040, hpb_j = n_j/B; every s_j and n_j a multiple of B so each batch has
exactly n_j/B histories of beamlet j). Beamlets use disjoint counter ranges (0037) and are statistically
independent; hence `Var(Σ_j x_j d_vj) = Σ_j x_j² σ_vj²` and no cross-beamlet covariance needs persisting.
Zero-deposition columns (A-2, A-3): a beamlet without deposit in the influence grid has nnz = 0, no threshold
applies, its discarded counts, energy, maximum dose and variance are 0 and its discarded fraction is undefined (NaN
with an explicit mask, decision 0040) at every τ including 0, while at τ = 0 the fraction is exactly 0 for every
column with positive in-grid energy; such beamlets have no uncertainty entries and the report distinguishes only
"defined" from "undefined (no deposition in the influence grid)".
`plan_uncertainty` excludes entries dropped by a threshold; the persisted
per-beamlet `influence__discarded_var_gy2` (Σ σ² over dropped entries) quantifies the omission. Batch values
per entry `[B, nnz]` are persisted by default (O3, tier-1 recomputation), opt-out recorded. Per-beamlet energy
balances are reduced exactly from per-history tally rows. Coverage is tested by independent replicas (plan UNC-2);
the broader coverage suite stays with V3-011.

### 6. Persistence

| Option | Assessment |
|---|---|
| (a) **`ionmc.result` schema 1.1, `kind = "influence"`** (selected) | 0045 reserved the kind and defined MINOR = optional blocks; reuses provenance, requested/effective, tiers, budgets, path safety, reports; one reader |
| (b) separate `ionmc.influence/1` | duplicates the whole provenance/contract machinery and a second strict reader; two formats to keep consistent |
| (c) influence arrays inside a 1.0 bundle | 1.0 readers reject unknown members (fail closed by design) — not possible without a version bump |

The composite grids of the run are persisted as ordinary grids in the same bundle. For a `BeamletSet`, `unit_direction`
is null-with-status (`not-applicable: beamlets`), a per-beamlet direction array `beamlet__direction` is persisted and
the report's profile axis is taken from beamlet 0 and labelled so. A species axis on the influence accumulators
(V3-009 fragment dose per beamlet) needs its own memory rule and a MINOR schema version. `EffectiveConfig.summary()`
stays
byte-identical for pencil configurations (O4, REG-2); for a `BeamletSet` its `source` block is the beamlet-set
summary (projectile, J, envelope, sha256 of the canonical beamlet table) and an `influence` key is added. The V3-006
negative control "schema 1.1 refused" moves to "1.2" by a dated note (intent unchanged). No bundle enters Git (REG-4).

### 7. Capability contract (declared and enforced together)

Supported: beamlets for every species/backend/precision the pencil source supports with EM physics; nuclear with
beamlets on python. Fail closed: influence + tallies or lookups (beamlet-resolved channels not implemented: a
composite LET silently labelled per-beamlet would be a no-op of the influence flag); warp + nuclear (unchanged, until
V3-005B); `cpu_workers > 1` with beamlets (single-process directive; untested path must not be declared supported);
mixed projectiles; n_j not a multiple of B; unknown influence grid; τ ∉ [0, 1); influence options with a pencil
source. Every refusal is `UnsupportedCombinationError` (exit 4) or `ConfigError` (exit 2).

### 8. Reuse

Within a run everything except source parameters and accumulator rows is built once. Across runs a `Simulation` keeps
its device context (tables, geometry, grids, kernel, nuclear setup); `Simulation.run(beamlets=subset)` re-runs a
subset with the original offsets (bit-identical columns) [superseded by A-3: `Simulation.influence_columns`]. No merge
of runs sharing a seed is offered.

Amendment (REVIEW-1d19386a, A-3; supersedes the REVIEW-1652a7c1 / REVIEW-15a83983 subset-run semantics): subset
runs produce no `Result` and no bundle. `Simulation.influence_columns(beamlets)` re-transports exactly the named
beamlets with the original offsets and returns a `ColumnSet` (CSC columns with batch values, std, nnz, discarded
accounting, per-beamlet balances, beamlet ids and configured offsets). A `ColumnSet` has no composite grids, no
statistics block and no traces; `save` refuses it (`TypeError`); it serves column-independence checks (INF-4) and the
V3-011 hooks. Fail closed (A-4 item 1, REVIEW-a52098eb): the API requires a `BeamletSet` AND `InfluenceOptions` and is
refused for a pencil source, for no `InfluenceOptions` and for ANY requested diagnostic (`requested_trace_count > 0`,
`track_end_positions`, `escape_records`, every other diagnostic/tally field a `Result` would carry); the test
enumerates the config dataclass fields at S6 and refuses every field outside the explicit "honoured" allowlist.
Validity (A-5, REVIEW-54588747): a `ColumnSet` carries the transport counters and is returned only when every
invalidating counter is zero, else `TransportLimitError` unconditionally; `allow_invalid_result` is a refused field.
Invalid subsets (empty, duplicate, out of range,
non-integer) are refused (`UnsupportedCombinationError`). Every bundle is a full run with
`RunOptions.n_histories` = Σ_{all j} n_j; the fields `run.histories_transported`, `run.beamlet_subset`,
`influence.subset` and `beamlet__transported` do not exist.

## Rejected alternatives (summary)

scipy.sparse / HDF5 (dependency outside numpy + warp, 0037); a separate influence schema; dense matrices; on-device
hashing; default thresholds; weight-proportional sampling; per-beamlet re-initialisation; resetting `c0` per beamlet
(stream reuse across beamlets, contradicts 0037); declaring LET influence or multi-process beamlets without tests.

## Uncertainty

No physical uncertainty is introduced. Software risks: kernel change in source sampling (guarded by REG-1 A16/T1
bitwise identity and INF-1); host-side sparsification cost and device→host transfer per group (measured, PERF-1);
bundle size with batch values (ESTIMATE for the workload: nnz ≈ 5·10⁶, batch values ≈ 0.8 GB uncompressed at B = 20;
recorded, opt-out explicit).

## Validation strategy

Frozen in `validation/plans/v3-007-acceptance.md` (rows BEAM, INF, SPA, UNC, PER, CAP, PAR, CLI, REUSE, PERF, REG,
DOC); exact identities defined before implementation (J = 1 equals pencil; Σ columns equals composite at the integer
level; per-beamlet closure INF-7; group/chunk invariance; column independence), negative controls showing the
identities discriminate, replica
coverage with a power control, three-backend parity, measurements without targets. No row is physics evidence.

## Expected tradeoffs

- Memory per group G·B·V·8 bytes; small G lowers occupancy only through more groups (each group re-launches).
- The searchsorted per history costs a few comparisons; measured in PERF-1.
- The full-physics influence workload of performance.json cannot run on Warp until V3-005B; V3-007 measures EM-only.

## Preserved obligations

Decisions 0037/0040/0041/0044/0045 unchanged except the additive schema 1.1; all V3-005/V3-008/V3-013/V3-006
verdicts as graded; single-process directive and diagnostic rules; V3-005B and V3-009 not claimed.

## Later validation outcome

(to be filled after S8)
