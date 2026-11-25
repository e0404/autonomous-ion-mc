# 0046 — Beamlet-resolved transport, sparse dose influence matrices and persisted influence bundles

- Status: Accepted (frozen with plan `validation/plans/v3-007-acceptance.md` at the S0 commit, 2026-11-25). Later
  changes are dated amendments. Not a clinical claim.
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
  1162.
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
group fitting `memory_budget_bytes` with the existing accumulators,
or user-set; recorded as effective (`derived: memory-budget`).

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
`Σ_j (n_j/N)(kept_j + discarded_j) = composite in-grid energy`). Rejected: absolute thresholds (unit- and
normalisation-dependent), thresholds on relative uncertainty (selects on noise ⇒ bias), silent default thresholds.

### 5. Per-beamlet uncertainty

Batch standard error per entry (decision 0040, hpb_j = n_j/B; every s_j and n_j a multiple of B so each batch has
exactly n_j/B histories of beamlet j). Beamlets use disjoint counter ranges (0037) and are statistically
independent; hence `Var(Σ_j x_j d_vj) = Σ_j x_j² σ_vj²` and no cross-beamlet covariance needs persisting.
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
subset with the original offsets (bit-identical columns). No merge of runs sharing a seed is offered.

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
