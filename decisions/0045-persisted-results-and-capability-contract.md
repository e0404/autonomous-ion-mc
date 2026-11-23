# 0045 — Persisted result bundle, configuration files, `ionmc run` and the requested/effective capability contract

- Status: freeze candidate, fourth revision 2026-11-23: revised after Codex REVIEW-bd983b92 of the S0 commit 00b5e3de
  (findings 11–15), Codex REVIEW-0f7aad46 of 00989a8b (findings 16–19), Codex REVIEW-c8397a61 of 034e3fb8 (findings
  20–22) and Codex REVIEW-34a0a3e8 of 84662725 (findings 23–24); each changes_required, resolved in place. The revised
  S0 commit is re-reviewed; the freeze is the first S0 commit whose Codex review passes (together with
  `validation/plans/v3-006-acceptance.md`, after an independent architecture/scientific review), before any V3-006
  code exists (planner). Later changes are dated amendments; earlier text is kept.
- Task: V3-006 (program plan task table: self-describing persisted results, human-readable reports and SVG plots,
  configuration files, `ionmc run` CLI, fail-closed capability contract with requested/effective configuration).
- Base: `v3/develop` 95937553 (the V3-013 squash merge, PR #92). V3-006 is implemented on its task branch from this
  base (orchestrator decision O1); no V3-013 seal binds this branch any more.
- Affects: architecture (public I/O API, CLI), reproducibility, provenance, capability contract, licensing of
  distributed outputs. **No physics, table, RNG, stopping-model or transport change.**
- Builds on: decision 0037 (numpy + warp only; NPZ + JSON sidecars; Markdown reports; SVG by an internal writer; unit
  and coordinate tables; Philox stream encoding), decision 0039 (fail-closed `validate`,
  `UnsupportedCombinationError`, `BackendUnavailableError`, `TransportLimitError`), decision 0040 (channels,
  undefined = NaN + mask, ratio statistics), decision 0041 (nuclear table id), decision 0044 §4–§5 (stopping identity
  block, `evidence_at_run`, `data_declared`/`data_used`, no silent switch, no environment override).
- Amends: the docstring rule of `src/ionmc/environment.py` "git is never executed" (§3, orchestrator decision O2):
  one read-only, timeout-bounded `git` query with every inherited `GIT_*` variable removed becomes the documented
  fallback of the code identity. Nothing else in decisions 0037–0044 is amended.
- Markers: **FACT** (read at 95937553, file:line), **ESTIMATE**, **TO VERIFY**, **PROPOSAL** (a design choice frozen
  by this decision). Not a clinical claim.

## Problem

V2-OUTPUT requires that, outside the producing process, a user recovers the code version/SHA and dirty state,
requested and effective physics, species/source, geometry/material definitions, transport and scoring grids with
coordinates, units, external data identifiers/hashes, model/reference versions, actual histories, batches, seeds/RNG
identity, normalisation, uncertainty definition/validity, escape/truncation, backend and relevant environment;
missing/undefined quantities are explicit; machine-readable data and human-readable summaries/profiles/plots/reports
exist; reopening is tested in a separate clean process. V2-CAP requires every public path to honour requested physics
or raise before returning, discoverable and documented support, persisted requested AND effective configuration and
discriminating dispatch tests. V2-USER requires a documented configuration/CLI workflow whose results are inspectable
without test code. None of this exists as a library product today.

## Context

**FACT (read at 95937553).**
- No library result writer exists. `Result` (`src/ionmc/simulation.py:306–353`) holds `valid`, `requested_config`,
  `effective_config`, `backend`, `device`, `precision`, `seed`, `rng`, `n_histories`, `n_batches`, `grids`
  (`GridResult`, `simulation.py:273–303`: mean/std energy and dose, `relative_uncertainty` with NaN where undefined,
  `defined_mask`, `mass_g`, `n_nonzero_batches`, `sparse_mask`, `low_batch_count`, `batch_energy_mev`, `valid`,
  `quantities`), `energy_balance`, `counters`, `timings`, `environment`, `diagnostics`, `transport_report`,
  `channel_raw`, `lookups`.
- `Simulation.__init__` calls `validate` before any transport (`simulation.py:710–714`); `run` raises
  `TransportLimitError` carrying the invalid result unless `allow_invalid_result` (`simulation.py:716–747`).
- `EffectiveConfig.summary()` (`src/ionmc/config.py:607–698`) carries backend, precision, seed, histories, batches,
  workers, chunking, source, physics (with `end_of_range_model`/`_effective`, `ls_requested`/`ls_effective`,
  `e_cut_effective_mev`), `mcs_mapping`, geometry (origin/spacing/shape, material **names only**, `z_exit_mm`),
  scoring grids, table identities and sha256, the stopping identity block `summary()["stopping"]`
  (`stopping_identity()`, `config.py:590–605`; decision 0044 §4), RNG and tallies; nuclear runs add `nuclear` and
  `e_cut_deuteron_mev` (`config.py:694–697`). It is bound by frozen V3-008 step documents
  (`validation/scripts/ion_em/r2_io.py`) and by the V3-013 identity records, so it is a compatibility surface.
- `validate` resolves more than the draft allowlist named: a non-unit `direction` becomes `unit_direction`
  (`src/ionmc/sources.py:17, 43–47`); `mcs_mapping_eval` becomes `"root"` and unused for an ion with multiple
  scattering off and is ignored for protons (`config.py:160–164, 448–453, 1262–1265`); the cutoff is `e_cut_mev`
  for a proton and `e_cut_ion_mev_per_u × A` for an ion (`config.py:523–530`); `e_cut_deuteron_mev` and
  `nuclear_table_id` are unused with nuclear off; `straggling_model`/`mcs_model` are unused when their switch is off;
  `chunk_histories` is unused on the python backend; the stopping identity adds the deuteron species for nuclear
  runs and a water material row when tallies are requested (`config.py:594–598`); an LS-variant id on z = 1 has the
  default as effective id (`src/ionmc/physics/stopping_selection.py:105`, `effective_model_id`).
- `select()` (`stopping_selection.py:112–144`) runs at `PhysicsOptions` construction and wraps every resolution error
  (unknown id, missing or mismatched cache object, parse failure) into `UnsupportedCombinationError` with the cause
  kept (`raise ... from exc`).
- V3-013 is integrated: the evidence registry (`src/ionmc/physics/stopping_evidence.py`) is append-only; the model
  `icru90-water-v1` now has two records, `consistency-only` and `topas-agreement-shared-lineage:V-C1:PASS` (graded at
  9f1fa8ec); the default `bethe-ionmc-v1` has `as-graded:0043`, the LS variants `unvalidated`, the NIST model
  `consistency-only`. No status is named `validated`. `evidence_at_run` therefore depends on the code version and
  must be persisted as it was at run time, never recomputed at load time.
- `capabilities()` (`simulation.py:618–630`, `_capabilities_base` 665–704) derives the tally section, the stopping
  model section (`_stopping_model_capabilities`, 632–663) and the species sections (`species_details`,
  `species_roles`, `unsupported_species`, with the ion energy domains from `ION_DOMAIN_MEV_PER_U`; 467–485 and
  following) from the enforcing objects; the proton energy domain, backends, geometry, sources and trace entries are
  free text. No test proves declared == enforced over the whole matrix.
- `describe_environment()` (`src/ionmc/environment.py:48–65`) reads `_build_info.json`, which nothing generates, so
  every run records `git_sha = None`; its docstring states "git is never executed". `.gitignore` does not list
  `_build_info.json`.
- `resolve_cache_dir` (`src/ionmc/data/cache.py:37–48`): explicit argument, else `IONMC_CACHE_DIR`, else
  `~/.cache/ionmc` (a default, not a user choice).
- CLI (`src/ionmc/cli.py`): `version`, `info`, `notices`, `data {list,fetch,verify,path,import,build}`; argparse
  errors and `data` without an action return 2 (`cli.py:99–100, 124–131`), `data` errors return 1
  (`cli.py:118–121`); no `run`. Entry point `ionmc = "ionmc.cli:main"` (`pyproject.toml:27–28`); build backend
  hatchling 1.32.4; `requires-python >= 3.12`.
- Validation scripts persist ad hoc (`np.savez_compressed` + JSON step document with `summary()` and the npz sha256).
- `VoxelGeometry.z_exit_mm` (`src/ionmc/geometry.py:28–52`) is an optional world exit plane; a source has an
  arbitrary direction and no axis is privileged (decision 0037, "Coordinates and grids").
- Host runner (`infrastructure/host_runner/host_runner.py:38, 179–195`): read-only bind of the task worktree at
  `/workspace`, `GIT_WORK_TREE=/workspace`, `HOME=/tmp/home`; Python 3.12.0 with numpy, warp and pytest only (no
  hatchling), no network.

**ESTIMATE / TO VERIFY.** `np.savez*` stamps zip members with the current local time, so npz file bytes are not
reproducible (ESTIMATE from CPython `zipfile`; TO VERIFY by row OUT-9). A hatchling build hook can write the build
stamp without a new build dependency (TO VERIFY by row PROV-5).

## Selected (PROPOSAL)

### 1. Result bundle `ionmc.result/1` (module `src/ionmc/results.py`)

A directory `<name>/`:

| File | Content | Authoritative? |
|---|---|---|
| `result.json` | strict JSON (`allow_nan=False`, sorted keys, UTF-8, two-space indent, one final `\n`): `schema = "ionmc.result"`, `schema_version = "1.0"`, `kind = "transport"` (`"influence"` reserved for V3-007), `valid`, the blocks of §2, the `arrays` catalogue and `files` (sha256 and size of every other file) | yes |
| `arrays.npz` | every array, written by a deterministic writer (`zipfile`, fixed `date_time = (1980, 1, 1, 0, 0, 0)`, fixed member order, `ZIP_DEFLATED` level 6, `np.lib.format.write_array`; read with `allow_pickle=False`) | yes |
| `report.md` | human-readable report (§5); a pure function of `result.json` + `arrays.npz` | derived |
| `plots/*.svg` | profiles (§5); a pure function of the same | derived |

- Catalogue entry per array: `name`, `dtype`, `shape`, `axes`, `grid` (scoring-grid name or `null`), `units`
  (decision 0037 vocabulary plus `"1"`, `"bool"`, `"count"`), `definition`, `content_sha256` = sha256 over
  `dtype.str | shape | C-order bytes` (compression-independent). Integrity = file sha256 AND every content digest.
- Persisted per grid: `energy_mev`, `energy_std_mev`, `dose_gy`, `dose_std_gy`, `relative_uncertainty`,
  `defined_mask`, `mass_g`, `n_nonzero_batches`, `sparse_mask`, `batch_energy_mev`, every `QuantityResult` array;
  geometry `material_index` and, if set, the density array; requested diagnostics.
- **Batch arrays are persisted by default** (orchestrator decision O3); `persist_batches = false` is an explicit,
  recorded choice (`{"status": "not-persisted", "reason": "user option persist_batches=false"}`), never a silent
  omission.
- **Never persisted:** stopping, range, LET-water and nuclear table arrays and any array derived node by node from
  licensed data (NIST, ICRU 90/Geant4, ENDF, AME). Only identities and digests (§2). Enforced by an array-name
  allowlist in the writer (decisions 0038/0044 licence roles, V1-MUST-020). User lookup tables are referenced by
  name, metadata and sha256; their values are not embedded.
- Undefined values: NaN in arrays with their masks; in JSON `{"value": null, "status": "undefined", "reason": ...}`;
  never 0, never omitted.
- Reader `load(path) -> LoadedResult` verifies schema, version, every hash, the allowlist and the cross-field
  invariants and **fails closed** (`ResultFormatError`) on an unknown field at any level, a missing required field, an
  unknown version, a hash mismatch, a NaN in a JSON number position, an array whose shape disagrees with its grid, or
  an npz member outside the catalogue. It never reconstructs transport internals and never recomputes
  `evidence_at_run`.
- **ZIP preflight on raw bytes (before `zipfile.ZipFile` is constructed, which builds one `ZipInfo` per member).**
  `load` checks the file size of `arrays.npz` against `max_total_bytes`, locates the end-of-central-directory record
  within the last 65 557 bytes, reads the entry count and central-directory size and offset from it (and from the
  ZIP64 end record if present), and refuses the archive unless entry count ≤ catalogue size, central-directory size ≤
  `max_central_directory_bytes` (default 1 MiB), offsets lie inside the file and ZIP64 values respect the same
  budgets; only then is the archive opened.
- **Untrusted bundles (resource budgets).** Before any array is read, `load` validates the ZIP central directory and
  the catalogue: member count ≤ catalogue size, no duplicate or unexpected member names, every member's declared
  uncompressed size and the catalogue's `shape` × `dtype.itemsize` within a byte budget derived from the catalogue and
  bounded by `max_array_bytes` and `max_total_bytes`, shapes with at most 4 dimensions and no negative or overflowing
  extent; each NPY header is parsed with `np.lib.format.read_magic` and `read_array_header_1_0`/`_2_0` (header length
  bounded) and must equal the catalogue entry before data are read; decompression is streamed with a byte budget and
  stops at the declared size; `result.json` is read with a size limit. Budgets are keyword arguments of `load` (and
  options of `ionmc verify`/`report`) with safe defaults (`max_total_bytes` = 8 GiB, `max_array_bytes` = 4 GiB,
  `max_json_bytes` = 16 MiB); before `json.loads`, a single pass over the raw bytes of `result.json` (no allocation
  proportional to content) counts every JSON token (objects, arrays, strings, numbers, `true`/`false`/`null`) and
  every object member against one total-token budget derived from the schema (default ≤ 2 000 000) and tracks the
  nesting depth (≤ 64); any violation raises `ResultFormatError` before allocation.
- **Safe names and bundle-relative paths.** Every derived file and npz member name comes from a fixed safe-name
  mapping of the user's grid, tally and lookup names: characters `[a-z0-9_-]` after lower-casing, length ≤ 64,
  collisions resolved by an index suffix; the mapping (name → file/member) is recorded in the catalogue. Every `files`
  entry and npz member name is bundle-relative and either single-level or `plots/<safe>.svg`: never absolute, no `..`
  or `.` segments, no other separator, no NUL. The writer refuses anything else and never follows or creates symlinks;
  the reader checks each entry with `os.lstat` (regular file, not a symlink) and `os.path.realpath` containment in the
  bundle directory before opening it, and refuses (`ResultFormatError`) otherwise.
- Writer `save(result, path, *, manifest=None, persist_batches=True, overwrite=False)` writes into a temporary sibling
  directory and renames (no partial bundle), refuses an existing path unless `overwrite`, writes an invalid result
  only when the run allowed it (`valid = false` and a banner in the report).
- Versioning: `"MAJOR.MINOR"`; MINOR adds optional blocks, MAJOR breaks; a reader accepts exactly the versions it
  knows; a frozen synthetic 1.0 fixture stays readable by later readers (orchestrator decision O10).

### 2. Blocks of `result.json`

- `requested`: canonical JSON of the `SimulationConfig` (§6), including `stopping_model` or, for a user `stopping`
  object (API only), `{"stopping_object": <identity>, "reconstructible": false}`.
- `input_manifest`: the `InputManifest` of §6 (config file and every referenced file by resolved absolute path and
  sha256, `base_dir`, `source`; orchestrator decision O8); for API runs `source = "api"` and every file pointer
  `{"status": "unavailable", "reason": "constructed in process"}`.
- `effective`: `EffectiveConfig.summary()` **byte-for-byte unchanged** (orchestrator decision O4), including the
  decision-0044 `stopping` block as produced at run time; plus `effective_extra` outside it: material definitions
  (composition, density, I requested/effective, Sternheimer parameters) of every geometry material, array references
  of `material_index`/density, the grid convention string (decision 0037: origin = corner of voxel (0,0,0), C order
  `[ix, iy, iz]`, mm), the profile axis of §5 and the resolution list of §4.
- `provenance`: `code` (§3); `environment` (python, platform, numpy, warp versions, device names, CUDA driver if any;
  no hostname, user name or home path); `data` (every dataset id, archive and member sha256 and per-run role,
  collected from `stopping.data_used`, the nuclear setup and the table identities; lookup sha256); `models` (stopping
  `model_id_effective`, `version`, `spec_sha256`, `evidence_at_run` verbatim; straggling, MCS, end-of-range and
  delta-electron model names; nuclear table id).
- `run`: backend, device, precision, production flag, seed, RNG identity (Philox4x32-10, counter/key encoding version,
  u01 mapping), `n_histories`, `n_batches`, histories per batch, `cpu_workers`, `chunk_histories`; `timings` labelled
  `"deterministic": false`; `output_path`; `cache_dir` and `cache_dir_source` (`"explicit" | "IONMC_CACHE_DIR" |
  "default"`, the latter being `~/.cache/ionmc`).
- **Path pointers (set P, frozen; orchestrator decisions O7/O8).** The only JSON pointers that may hold an absolute
  path, and therefore a home directory or user name: `/run/output_path`, `/run/cache_dir`, `/input_manifest/base_dir`,
  `/input_manifest/config_file/path`, `/input_manifest/files/<i>/path`, and every `cache_dir` key inside
  `/effective/stopping` and `/effective/nuclear` (the decision-0044 and nuclear blocks as `summary()` produces them).
  They are recorded as user-supplied or user-environment paths. Outside P no string contains the home path, and no
  path component (split at `/` and `\`) and no whole non-path string equals the user name or the hostname. Reports and
  plots never print a P value (input files appear by basename and sha256, caches by dataset id).
- `statistics`: normalisation ("per primary history"), uncertainty definition (batch standard error, B batches,
  history h in batch h mod B; ratio quantities by the delta method of decision 0040), `low_batch_count`, the
  definitions of `defined_mask`/`sparse_mask`.
- `accounting`: every `energy_balance` field (escaped, truncated, residual), every counter, `valid` and the reasons.
- `units`: the decision-0037 unit table as data.

### 3. Code identity (orchestrator decision O2)

- `code.package_tree_oid`: computed at run time from the imported `ionmc` package directory as a **Git tree
  object id** (SHA-1 tree/blob encoding, regular-file modes, excluding `__pycache__`, `*.pyc` and the build stamp);
  equal to `git rev-parse <commit>:src/ionmc` exactly when the code that ran equals that commit's package tree. Needs
  no git binary; works for wheel, checkout and host snapshot. Authoritative identity and dirty test.
- `code.git_commit`, `code.git_dirty`, `code.source` (`"stamp" | "git" | "unavailable"`): (a) a build stamp
  `ionmc/_build_info.json` present in an installed wheel: `git_dirty = runtime oid != stamped oid`; (b) else, if the
  package lies in a Git work tree, one read-only query (`git rev-parse HEAD` and `HEAD:src/ionmc`) with every `GIT_*`
  variable removed, `cwd` = the package directory and a timeout; `git_dirty = HEAD tree oid != runtime oid`; (c) else
  `git_commit = null`, `git_dirty = null`, `status = "unknown"` and the reason.
- **Build stamp.** Written by a hatchling custom build hook (built into hatchling; no new build dependency:
  `[build-system] requires` and `uv.lock` stay unchanged, `pyproject.toml` gains only the hook table) when a wheel is
  built with `uv build --wheel` from a clean checkout in the local validation environment (never on the host runner,
  which has no hatchling and no network). The hook computes the tree oid and one read-only `git rev-parse HEAD`
  (`GIT_*` stripped), writes the stamp into a temporary directory and maps it into the wheel with `force_include`;
  it never writes into `src/ionmc` (S3 also adds `src/ionmc/_build_info.json` to `.gitignore` as a guard). A build
  without a reachable repository (for example a wheel built from an sdist) stamps `commit = null`,
  `status = "unknown"`; this outcome is pre-declared, never an error.
- `ionmc_version` from `_version.py`. Apart from this query the library runs no subprocess and never uses the network.

### 4. Capability contract: declared, enforced, requested, effective

1. **Declared** = `capabilities()`; every entry naming a supported or rejected combination is derived from the
   constant or predicate that `validate` enforces (already so for tallies, stopping models and the ion species
   domains), extended to the proton energy domain, backends × precision × (nuclear, tallies, trace, workers), geometry
   and source types. A free-text entry either becomes derived or is labelled `"informative"`. Enabling a combination
   later (V3-005B warp nuclear, V3-009 ion nuclear, V3-007 beamlets) changes declaration and enforcement together.
2. **Enforced = declared over the enumerated dimensions**: a sweep enumerates the finite declared matrix (species ×
   energy-domain class × backend × precision × nuclear × tallies × trace × stopping model) plus boundary energies and
   asserts declared-supported ⇒ `validate` succeeds; declared-unsupported ⇒ `UnsupportedCombinationError` or
   `BackendUnavailableError` before transport; no third outcome. Tabulated models are swept against `validate` only
   where their cache is present (local and host tiers); without cache only their declared domains are checked. Every
   remaining `validate` rule (geometry kind and voxel constraints, scoring-grid count and shape, source sigmas,
   cutoffs, end-of-range and MCS modes, worker and chunk limits, diagnostics) has a targeted positive/negative pair.
   No global claim about combinations outside the enumerated dimensions is made.
3. **Requested → effective reconciliation** (`reconcile(requested, effective) -> tuple[Resolution, ...]`, run inside
   `Simulation.__init__` after `validate`, persisted in `effective_extra.resolutions`). Every requested field has one
   row of the frozen table below with one of four classes: **equal** (effective value must equal the requested one),
   **derived** (a documented deterministic rule, recorded with the rule name), **unused** (not applicable for this
   configuration; recorded with the reason, never silently dropped) and **requested-only** (an execution option with
   no effective counterpart; persisted in `requested`). Any difference outside the table raises
   `UnsupportedCombinationError("requested X, effective Y")` before transport. A requested feature without its
   effective implementation object (tallies → channel plan, nuclear → nuclear setup, lookups → lookup channels,
   multiple scattering → MCS model, straggling → straggling model) raises.

   | Requested field | Effective pointer | Class and rule |
   |---|---|---|
   | `source.projectile`, `position_mm`, `kinetic_energy_mev`, `energy_sigma_mev`, `lateral_sigma_mm` | `/effective/source/*` | equal |
   | `source.direction` | `/effective/unit_direction` | derived `normalise` (d / ‖d‖; `sources.py:43–47`); equal if already unit |
   | `geometry` (`BoxPhantom`) | `/effective/geometry/*` | derived `box-to-voxel` (one-voxel `VoxelGeometry`) |
   | `geometry` (`VoxelGeometry`) origin, spacing, shape, materials, `z_exit_mm`; `material_index`, density | `/effective/geometry/*`; `/effective_extra/geometry/*` (content sha256) | equal |
   | `scoring[i]` | `/effective/scoring/<i>` | equal |
   | `physics.stopping` / `stopping_model` | `/effective/stopping/model_id_effective` | derived: LS-variant id on z = 1 → `bethe-ionmc-v1` (`effective_model_id`); `stopping` object (API only) → `spec_of` id or `unregistered`; otherwise equal. Exactly one of the two is required (decision 0044 §2, `select()`); absence fails at `PhysicsOptions` construction and is not a resolution |
   | stopping identity species and materials | `/effective/stopping` (species, tables) | derived: deuteron added for nuclear runs; water row added when tallies are requested and no geometry material is water (`config.py:594–598`) |
   | stopping option `ls` | `/effective/physics/ls_effective` | derived `ls-z1-off` (z = 1 → `off`); equal for z ≥ 2 |
   | `physics.nuclear`, `straggling`, `multiple_scattering`, `delta_electrons` | `/effective/physics/*` | equal |
   | `physics.straggling_model` | `/effective/physics/straggling_model` | equal if `straggling`; unused otherwise |
   | `physics.mcs_model` | `/effective/physics/mcs_model` | equal if `multiple_scattering`; unused otherwise |
   | `physics.mcs_mapping_eval` | `/effective/mcs_mapping/evaluation` | ion with MCS on: equal; ion with MCS off: derived `mcs-off-root` (`"root"`, `unused = true`); proton: unused |
   | `physics.e_cut_mev` | `/effective/physics/e_cut_effective_mev` | proton: equal; ion: unused (must stay at its default 2.0) |
   | `physics.e_cut_ion_mev_per_u` | `/effective/physics/e_cut_effective_mev` | ion: derived `per-u-times-A`; proton: unused |
   | `physics.e_cut_deuteron_mev`, `nuclear_table_id` | `/effective/physics/e_cut_deuteron_mev`, `/effective/nuclear` | nuclear on: equal; nuclear off: unused |
   | `physics.end_of_range_model` | `/effective/physics/end_of_range_model_effective` | `"default"` → derived `point_v1` (z = 1) / `continuous_residual_v1` (z ≥ 2); otherwise equal |
   | `physics.max_step_mm`, `max_energy_loss_fraction`, `range_alpha`, `range_rho_f_mm`, `short_step_fraction`, `truncated_hinge_diagnostic` | `/effective/physics/*` | equal |
   | `run.backend`, `precision`, `seed`, `n_histories`, `n_batches`, `cpu_workers` | `/effective/*` | equal; `production` derived (float32 Warp only) |
   | `run.max_steps` | `/effective/max_steps`, `/effective/max_steps_origin` | `None` → derived `computed`; otherwise equal (`"user"`) |
   | `run.chunk_histories` | `/effective/chunk_histories` | Warp backends: equal; python: unused |
   | `run.allow_invalid_result`, `worker_timeout_s`, `memory_budget_bytes` | — | requested-only |
   | `diagnostics.*`, `tallies[i]`, `lookups[i]` | `/effective/tallies`, diagnostics arrays, lookup sha256 | equal (each request present in the compiled plan) |

4. **Result-level honoured check**: `save` and `load` verify that every requested feature left its trace (tally arrays
   for every request, nuclear tallies iff nuclear, trace arrays iff requested); a missing trace is an error, never a
   thinner bundle.
5. No negotiation, downgrade or fallback of backend, precision or model (decisions 0039, 0044 §5).

### 5. Reports and plots

- `report.md`: identity (tree oid, commit or "unknown", dirty, version, schema), validity banner, requested vs
  effective table with resolutions, source/geometry/materials (input files by basename and sha256), physics and
  stopping model with every `evidence_at_run` record reproduced verbatim (status, row, plan, graded-at SHA), data
  identities and roles, run/RNG, statistics definition, per-grid summary, energy balance and counters, missing-items
  list, and the line "Not a clinical claim." The report never prints timings or a value of the path set P (§2). The
  generated text (the report with the verbatim evidence records removed) never matches the regex `\bvalidated\b`; no
  evidence status is paraphrased; a `topas-agreement-shared-lineage:*` record is followed by the note "shared
  stopping-data lineage; transport agreement only".
- SVG (internal writer `src/ionmc/svg.py`): per grid the depth dose profile ± 1σ band, the depth-integrated energy
  profile and one lateral profile at the depth of maximum; per LET tally the LET_d/LET_t depth profile over defined
  voxels only (gaps where undefined). Fixed canvas, `format(x, ".6g")`, no timestamps, random ids or embedded fonts.
- Profile axis (orchestrator decision O9): the grid axis of the largest |component| of the source unit direction
  (ties → lowest axis index), recorded in `effective_extra.profile_axis`; no axis is privileged by the engine.
- `regenerate_report(bundle)` reproduces `report.md` and every SVG byte-identically; `save` calls the same function.

### 6. Configuration files `ionmc.config/1` (orchestrator decision O5)

- TOML input read with stdlib `tomllib`; canonical JSON (the `to_dict` form) is what the bundle persists as
  `requested`. Strict: unknown table or key, wrong type or missing required key → `ConfigError` (a subclass of
  `UnsupportedCombinationError`) naming the key path (e.g. `physics.nuclear`).
- Tables: `[source]`, `[geometry]` (`type = "box"` or `"voxel"` with `material_index = {file, sha256}` and optional
  `density = {file, sha256}`), `[[materials]]` (built-in name or inline definition checked as by the API),
  `[[scoring]]`, `[physics]` (every `PhysicsOptions` field except the `stopping` object; `stopping_model` is
  **required**: absence → `ConfigError` naming `physics.stopping_model`; the API keeps decision 0044's rule
  unchanged), `[run]`, `[diagnostics]`, `[[tallies]]`, `[[lookups]]` (`file`, `sha256`), `[output]` (`path`,
  `persist_batches`). `nuclear` has no default, as in the API.
- Relative paths resolve against the config file's directory (`..` allowed); the resolved absolute path and the
  sha256 are recorded; the sha256 is mandatory for every referenced file (O8). No environment variable overrides a
  physics or run field (decision 0044 §5).
- **`InputManifest`** (immutable; produced by the loader): `source = "toml" | "api"`, `base_dir`, `config_file` (path,
  sha256) and `files` (one entry per referenced file: role `material_index | density | lookup:<name>`, resolved
  absolute path, file sha256, `content_sha256` of the array and `item_sha256` = sha256 of the canonical JSON of every
  serialised field of the item: for a `LookupTable` all fields incl. `citation`, `license`, `source`, `resampling` and
  the array content digest; for a geometry array its dtype, shape, content digest and the geometry fields that
  describe it). It also carries `config_sha256`, the sha256 of the canonical JSON of `to_dict(config, manifest=m)`
  computed by the loader with its own manifest `m`, which binds it to the complete configuration. It is a side
  carrier: `Simulation` and `SimulationConfig` are unchanged; the CLI passes it to `save(..., manifest=...)`.
  `InputManifest.verify(config)` and `save` recompute `config_sha256` as the sha256 of
  `to_dict(result.requested_config, manifest=m)` with the same manifest (deterministic by construction) and refuse any
  mismatch (`ConfigError`, exit 2, no bundle): a manifest of another configuration or a configuration mutated after
  loading is never attributed to these files. `save` re-hashes every listed file and checks every `content_sha256`
  against the arrays of the requested config; a mismatch raises `ConfigError` (exit 2, no bundle). Without a manifest
  `save` records `source = "api"` with `unavailable` file pointers.
- **Round-trip semantics.** Equality of configurations is canonical and semantic: `to_dict` equality, arrays compared
  by `content_sha256` (`VoxelGeometry` and `LookupTable` use identity equality, `eq=False`). Reconstruction API:
  `from_dict(d, *, manifest: InputManifest | None = None)`. Serialisation takes the manifest: `to_dict(config, *,
  manifest: InputManifest | None = None)`. With a manifest, every array (voxel `material_index`, density, lookup
  table) whose `item_sha256` (computed from the in-memory object as above) equals that of a manifest file entry is
  emitted as the reference `{"ref": "<role>", "sha256": ..., "content_sha256": ...}` (matching by item digest, so
  numerically equal data with other provenance or resampling metadata never match; an in-process item identical in
  every serialised field matches the same reference by construction); without a manifest, or for an array without a
  matching entry, it is emitted as `{"status": "not-reconstructible", "content_sha256": ...}`. `from_dict` reloads it
  through the manifest entry of that role (path from the manifest, file re-hashed, content digest checked); without a
  manifest, or with a manifest lacking the role or with a different digest, it raises `ConfigError` naming the
  pointer. The exact round trip `to_dict(from_dict(to_dict(c, manifest=m), manifest=m), manifest=m) == to_dict(c,
  manifest=m)` is guaranteed for the reconstructible subset: every field except user `stopping` objects and arrays
  that exist only in memory (API-created lookup tables, geometry arrays without a file); API configs without
  file-backed arrays need no manifest. The excluded fields serialise to `{"status": "not-reconstructible", "reason":
  ...}` and `from_dict` refuses them (`ConfigError`), never substituting a default.


### 7. CLI

`ionmc run CONFIG.toml [--out DIR] [--overwrite]` (load → validate → reconcile → run → save → report), `ionmc check
CONFIG.toml` (load, validate and reconcile; print the effective summary; no transport), `ionmc capabilities
[--nuclear]` (JSON), `ionmc verify DIR` (load with every check), `ionmc report DIR [--check]` (regenerate; `--check`
compares bytes). Existing subcommands and their codes unchanged. Before constructing `PhysicsOptions` the loader
checks the `stopping_model` id against `STOPPING_MODELS` (`ConfigError`) and verifies every cache object the selected
model and the nuclear table need for this run (`cache.verify`; failure → `DatasetError`, a new class that is not an
`UnsupportedCombinationError`). **Exception → exit code (frozen; first matching line wins, so subclasses precede their
bases):**

| Order | Condition | Exit | Bundle |
|---|---|---|---|
| 1 | argparse usage error | 2 | none |
| 2 | `ConfigError` (TOML syntax, unknown table/key, wrong type, missing key, unknown `stopping_model` id, sha256 mismatch of a config-referenced input file) | 2 | none |
| 3 | `DatasetError`; or an `UnsupportedCombinationError` whose `__cause__` chain contains `cache.IntegrityError`, `OfflineError` or `FileNotFoundError` (a missing or mismatched cache object or a dataset parse failure reached through `select()`) | 1 | none |
| 4 | `UnsupportedCombinationError` (domain, role, combination, reconciliation) | 4 | none |
| 5 | `BackendUnavailableError` | 5 | none |
| 6 | `TransportLimitError` (invalid result, `allow_invalid_result = false`) | 3 | none |
| 7 | run completed with `valid = false` (`allow_invalid_result = true`) | 3 | written, `valid = false`, banner |
| 8 | `ResultFormatError` (`save`, `verify`, `report --check`) | 6 | none (temporary directory removed) |
| 9 | any other exception (a defect) | 70 | none; traceback printed |
| — | success | 0 | written |

No partial bundle on any failure. `run` executes in one process (`cpu_workers = 1` unless the user sets it; the
single-process directive governs every validation run).

## Rejected alternatives

| Topic | Alternative | Reason rejected |
|---|---|---|
| Container | HDF5 / Zarr / netCDF | runtime dependency outside numpy + warp (decision 0037; host venv) |
| Container | single `.npz` with embedded JSON | not inspectable without code (V2-USER) |
| Container | pickle of `Result` | executes code on load; not self-describing; Python-version-coupled |
| Container | pure JSON with arrays | size of batch arrays; lossless only via repr |
| Arrays | `np.savez_compressed` | non-reproducible zip timestamps; file hash not a content identity |
| Schema | lenient reader ignoring unknown fields | silent loss of meaning (V2-OUTPUT, V2-CAP) |
| Schema | extend `EffectiveConfig.summary()` | bound by frozen V3-008 step documents and V3-013 identity records (O4) |
| Code identity | build stamp only | absent for checkout and host-snapshot runs, which produce the evidence runs |
| Code identity | runtime `git` only | absent for wheels; `GIT_*` pitfalls; dirtiness only via status |
| Code identity | caller-supplied `--code-sha` | unverifiable |
| Contract | negotiation / automatic downgrade | silent switch (V2-CAP, decision 0044 §5) |
| Contract | documentation-only capability table | drifts from enforcement; untestable |
| Config | YAML | new dependency |
| Config | JSON as hand-written input | no comments; JSON stays the persisted canonical form only |
| Config | Python scripts as configs | executes code; not declarative |
| Report | HTML / matplotlib | dependency; decision 0037 fixes Markdown + internal SVG |
| Bundle | persist stopping/LET-water tables | licence roles (0038/0044), V1-MUST-020; identities suffice |
| Batches | omit batch arrays by default | lossy; re-analysis and parity need them (O3) |
| Paths | strip every absolute path | breaks the decision-0044 identity (`cache_dir`) and O8; O7 restricts absolute paths to the frozen path set P (§2) |

## Uncertainty

No physical uncertainty is introduced or changed. Software uncertainties: the reproducibility of deflate output across
zlib versions (the content digest, not the npz file hash, is the cross-host identity; row OUT-9 documents it); Warp
version strings and device names on hosts without CUDA; the hatchling hook behaviour (PROV-5, TO VERIFY).

## Data roles

Bundles carry dataset identities, hashes and per-run roles exactly as the decision-0044 block and the nuclear setup
state them; never table values. The roles of decisions 0038, 0041 and 0044 are unchanged. User lookup data are
referenced by hash. Evidence statuses are reproduced verbatim from `evidence_at_run`; no V3-006 row produces physics
evidence and no evidence record is appended by V3-006.

## Validation strategy

Frozen in `validation/plans/v3-006-acceptance.md` (rows OUT, PROV, CAP, CFG, REP, CLI, DOC, REG). Every row is a
software, capability, round-trip, determinism or self-consistency row; **no row is physics evidence**. Discriminating
tests are defined before implementation: clean-process round trip, unknown-field fail-closed, unsupported combination
refused with a stable exit code, byte-identical report regeneration, no table values in Git, single process only.

## Expected tradeoffs

- Batch arrays can be large (ESTIMATE: 20 batches × 1e6 voxels × 8 B ≈ 160 MB per grid before compression); the
  opt-out is recorded.
- The tree-oid computation reads every package file at run start (ESTIMATE < 50 ms for ≈ 80 files).
- A strict reader makes every extension (V3-007 influence, V3-009 fragments, V3-012 resources) a MINOR bump with
  reader support in the same task.
- `requested` is not reconstructible for runs that pass a user `stopping` object; recorded.
- One documented exception to "no environment details": the frozen path set P (§2; O7/O8).

## Preserved obligations

- Decision 0044 in full: no silent switch, no environment override, the default and every registered spec unchanged,
  evidence statuses never named `validated`, no evidence transfer between models.
- V3-013 results as graded (V-C1 PASS, V-C2 PASS, V-C2-D PASS, V-C1-A case 3 FAIL recorded); V3-008 verdicts incl.
  row I3 FAIL and I3-FU-1 INCONCLUSIVE; the V3-008 integration note's unresolved qualification obligations.
- Single-process directive and the operator diagnostic rules (one heavy job, logged pytest runs); multiprocess
  checks stay deferred and recorded.
- V1-MUST-009/010 (beamlets, influence matrices), V3-005B (warp nuclear) and V3-009 (ion nuclear) are not claimed.

## Pre-freeze review log (2026-11-23; findings recorded, resolutions applied in place before the freeze)

Independent pre-freeze review (reviewer agent, read-only) of this decision and the plan; findings 1–10:

1. **Blocking.** The requested → effective allowlist missed resolutions that `validate` makes (unit direction,
   `mcs_mapping_eval`, ion cutoff per u × A, unused nuclear/straggling/MCS/chunk options, added deuteron species and
   water row, LS-variant id on z = 1) → frozen field-by-field table with the classes equal / derived / unused /
   requested-only (§4 item 3; Context); plan row CAP-3 requires zero unlisted differences over the CFG-1 set.
2. Numbering of the V3-005 seed amendment → handled by the orchestrator (dated note of 2026-11-23 in plan v3-005);
   the plan §C refers to it and obliges V3-005B to restate the exclusion in its seed guard.
3. **Important.** O7 vs O8: absolute input and config paths contradicted the privacy rule; the default cache is not a
   user choice; a substring scan for the user name false-positives → frozen path set P of JSON pointers, scan at path
   components only, `cache_dir_source` recorded (§2); plan row PROV-6.
4. **Important.** Exit codes ambiguous (`ConfigError` subclasses `UnsupportedCombinationError`; `select()` raises
   unknown id and missing cache alike; invalid-result code unstated) → loader pre-checks and a frozen exception →
   exit-code table with subclass precedence (§7); plan rows CLI-2, CAP-7.
5. **Important.** Build stamp under-specified (sdist path without `.git`, host without hatchling, stamp not ignored)
   → `uv build --wheel` from a clean checkout in the local environment, stamp mapped into the wheel only, `.gitignore`
   guard, pre-declared `commit = null` outcome, no new build dependency asserted (§3); plan row PROV-5.
6. **Important.** Determinism (b) excluded too little → reports never print paths or timings; the excluded pointers
   are `run.timings` and the path set P (§2, §5); plan §C.
7. **Minor.** Context stale (derived species sections, `cli.py:99–100`, `Result` fields) → Context corrected;
   `GridResult.valid` added to the plan's OUT-0 clause map.
8. **Minor.** The `validated` regex collided with verbatim evidence records → regex scoped to generated text, plus a
   registry guard (§5); plan row REP-3.
9. **Minor.** Host-runner provenance unasserted → plan row PROV-7 (HR: `code.source = "git"`, clean tree).
10. **Minor.** "No integer literal" was literally false (hex digests) → "no decimal seed literal" with the verified
    window facts (plan §C).

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
    dimensions"; the global "zero undeclared combinations" claim removed (§4 item 2).
15. **Important.** `load`/`ionmc verify` lacked resource-exhaustion defences for untrusted bundles → pre-allocation
    validation of the ZIP directory, catalogue and NPY headers, streamed decompression under configurable byte
    budgets with safe defaults (decision 0045 §1); row OUT-10 with synthetic KB-scale adversarial fixtures in S2.

Codex REVIEW-0f7aad46 (REVIEW-0f7aad4697ae4f6b9679d724571922d2, 2026-11-23) of the revised S0 commit 00989a8b:
changes_required, four important findings, resolved in place:

16. **Important.** File-backed arrays could not be rebuilt by `from_dict(to_dict(c))` (the manifest is a side
    carrier; `VoxelGeometry` keeps no path, `LookupTable` only its file hash) → explicit reconstruction API
    `from_dict(d, *, manifest=None)`: file-backed arrays serialise as references (`ref`, `sha256`, `content_sha256`)
    and are reloaded through the manifest's file entries with both digests checked; a missing manifest or role raises
    `ConfigError` naming the pointer; in-memory arrays stay not-reconstructible (decision 0045 §6; rows CFG-1, CFG-7).
17. **Important.** Nothing bound the manifest to the complete configuration → `config_sha256` (canonical JSON of
    `to_dict(config)`) in the manifest, recomputed by `InputManifest.verify(config)` and `save`; a mismatch raises
    `ConfigError` (exit 2, no bundle) (decision 0045 §6; row PROV-8 with a swapped-manifest and a post-load-mutation
    negative test).
18. **Important.** CAP-1b omitted tallies and trace for the tabulated models → CAP-1b varies tallies (none / LET /
    full set with a synthetic lookup) and trace per model, backend and precision, plus cross-model cases of the
    water-row path (`config.py:594–598`); the enumerated dimensions are stated (row CAP-1b).
19. **Important.** `zipfile.ZipFile` builds one `ZipInfo` per member on open, before any member-count check → a
    raw-byte preflight of the end-of-central-directory record bounds file size, entry count and central-directory
    size (ZIP64 included) before the archive is opened (decision 0045 §1; row OUT-10 with a programmatically built
    oversized archive).

Codex REVIEW-c8397a61 (REVIEW-c8397a617d184374a1aa585b772da729, 2026-11-23) of the revised S0 commit 034e3fb8:
changes_required, three important findings, resolved in place (CAP-1b and the ZIP preflight accepted):

20. **Important.** `to_dict(config)` could not tell file-backed arrays from identical in-process arrays, so it could
    neither emit the references nor reproduce `config_sha256` → the manifest is passed to serialisation:
    `to_dict(config, *, manifest=None)` emits a reference for every array whose `content_sha256` matches a manifest
    entry and `{"status": "not-reconstructible", "content_sha256": ...}` otherwise; `config_sha256` is computed and
    recomputed with the same manifest (decision 0045 §6; rows CFG-1, CFG-7, PROV-8).
21. **Important.** A `result.json` within the size limit could still allocate far beyond its size in `json.loads`
    → default `max_json_bytes` lowered to 16 MiB and a single-pass structural pre-scan of the raw bytes (container,
    string and depth budgets) before decoding (decision 0045 §1; row OUT-10).
22. **Important.** Names from grids, tallies and lookups could become traversal paths, and the reader trusted the
    `files` catalogue → fixed safe-name mapping for every derived file and npz member, bundle-relative single-level
    or `plots/<safe>.svg` entries only, no symlinks, `realpath` containment and `lstat` checks in the reader, refusal
    in the writer (decision 0045 §1; new row OUT-11).

Codex REVIEW-34a0a3e8 (REVIEW-34a0a3e812f549e5b50aea78e47ac571, 2026-11-23) of the revised S0 commit 84662725:
changes_required, two important findings, resolved in place (path safety and the rest accepted):

23. **Important.** Matching a lookup to a manifest entry by `content_sha256` ignored its provenance fields
    (`citation`, `license`, `source`, `resampling`; `src/ionmc/lookup.py:133–145`) → every manifest file entry carries
    `item_sha256` over every serialised field of the item, and a reference is emitted only when the in-memory item's
    `item_sha256` equals the entry's; otherwise the item is not-reconstructible (decision 0045 §6; rows CFG-7,
    PROV-8).
24. **Important.** The JSON pre-scan did not count numbers, literals or object members → every JSON token and every
    object member counts against one schema-derived total-token budget (default 2 000 000) besides the depth budget
    (decision 0045 §1; row OUT-10 with a scalar-heavy fixture).

## Later validation outcome

(empty until graded)

## Amendment A-1 (2026-11-23; after Codex REVIEW-982a8c25 of 9e1eda65): the reader's binding contract, enumerated

The clause "cross-field invariants" of §1 is replaced by the following enumerated contract (a prospective amendment:
no row outcome has been observed; the S2a rows are implementation checks). `load` never opens a file outside the bundle
and never runs transport; everything below is computed from `result.json`, `arrays.npz` and the library code.

1. **Re-derivation of `effective`.** When `requested` contains no lookup-table reference, `load` rebuilds the
   `SimulationConfig` from `requested` (file-backed geometry arrays are taken from the persisted `material_index` /
   `density` members, whose content digests the manifest entries already bind) and calls `validate(config)`;
   `effective` must equal `EffectiveConfig.summary()` of that configuration byte-for-byte (canonical JSON). The
   re-derivation needs the data cache for tabulated stopping models (`IONMC_CACHE_DIR` or the bundle's recorded
   `cache_dir_source` rule); when the needed cache object is absent, `load` fails closed (`ResultFormatError`
   "effective not re-derivable: dataset … not cached") unless called with `rederive=False`, which records
   `effective_binding = "typed-only"` in the `LoadedResult` and is never used by `ionmc verify`.
2. **Lookup configurations.** When `requested` references lookup tables (values never persisted), `effective`
   is bound field by field: every field that does not depend on lookup values is compared with the re-derived
   configuration (built with the lookup tables replaced by placeholder tables of the persisted provenance, shape
   and units), and the lookup-dependent channel-plan fields (`k`, `bound_per_history`, `bounds`) are bound through
   the lookup `content_sha256`/`item_sha256` of the manifest only; this limit is recorded in `LoadedResult`.
3. **Stopping identity.** The persisted `effective.stopping` block is checked with
   `stopping_selection.validate_identity_block` against the registry at load time (spec and table identities; the
   `evidence_at_run` list is persisted verbatim and never recomputed, as §1 says).
4. **Accounting identities (relative tolerance 1e-9 on sums; clarified 2026-11-23 to match the implementation).** Per grid: `in_grid + quantization + outside
   = deposits`; `in_grid_mev[i] = sum(energy_mev[i]) × n_histories`; closure: `source − deposited − escaped −
   truncated − residual = 0` with the documented signed nuclear terms; nuclear results persist the exact key set
   of the engine's nuclear energy balance; counters consistent with `valid`.
5. **Array invariants (tolerances stated in the reader).** Standard errors finite and ≥ 0; `dose_gy = energy_mev ×
   MeV→Gy / mass_g` where `mass_g > 0` (NaN where it is 0); `dose_std_gy` the same transform of `energy_std_mev`;
   `relative_uncertainty = energy_std / energy` where defined (NaN otherwise); `defined_mask ⇔ n_nonzero_batches >
   0`; `sparse_mask ⇔ n_nonzero_batches < low-batch threshold`; when batch arrays are persisted, their mean and
   batch standard error reproduce `energy_mev` and `energy_std_mev` (relative tolerance 1e-9); quantity companion
   arrays follow the decision-0040 ratio rules (NaN wherever their mask is False).
6. Everything else in `effective` and `provenance` is bound by items 1–3; nothing in the bundle is "typed-only"
   except what item 2 records for lookup configurations and (clarified 2026-11-23) a `requested` configuration
   with a not-reconstructible user `stopping` object (§6; API only), which cannot be rebuilt and loads with
   `effective_binding = "typed-only"` and a recorded limit; `ionmc verify` reports that limit. Rows OUT-4/OUT-7 gain the corresponding negative controls
   (every identity and invariant above broken one at a time → `ResultFormatError`).
