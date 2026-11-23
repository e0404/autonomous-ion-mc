# Results format reference

Schema `ionmc.result`, version `1.0`, kind `transport`. A reader accepts exactly the versions it
knows; a MINOR change adds optional blocks and a MAJOR change breaks the format. The constants
below are those of `ionmc.results`; the documentation test `tests/ionmc/test_v3_006_docs.py`
checks this page against them. See the [overview](index.md) for the concepts.

## `result.json`

Strict JSON: `NaN` and infinity are never written, keys are sorted, UTF-8, two-space indent and one
final newline. Top-level keys:

| Key | Content |
|---|---|
| `schema` | `ionmc.result` |
| `schema_version` | `1.0` |
| `kind` | `transport` (`influence` is reserved) |
| `valid` | whether the run produced a valid result |
| `requested` | canonical JSON of the configuration as given |
| `input_manifest` | the config file and every referenced file by resolved path and sha256, `base_dir`, `source` (`toml` or `api`) |
| `effective` | `EffectiveConfig.summary()`, byte for byte |
| `effective_extra` | `materials`, `geometry_arrays`, `grid_convention`, `profile_axis`, `resolutions`, `resolutions_status`, `lookup_plan` |
| `provenance` | `code` (package tree object id, commit, dirty flag), `environment`, `data` (dataset ids and hashes), `models` (including `evidence_at_run`, verbatim) |
| `run` | backend, device, precision, production flag, seed, RNG identity, `n_histories`, `n_batches`, `histories_per_batch`, `cpu_workers`, `chunk_histories`, `timings` (labelled non-deterministic), `output_path`, `cache_dir`, `cache_dir_source` |
| `statistics` | normalisation (per primary history), uncertainty definition, mask definitions, the status of batch arrays and quantity batch arrays, one record per grid |
| `accounting` | `energy_balance`, `counters`, `channel_counters`, `valid`, `reasons` |
| `units` | the units table below, as data |
| `arrays` | the array catalogue |
| `files` | `{name: {"sha256", "size"}}` of every other file |

`effective` is never changed by the persistence layer; `effective_extra` holds what the bundle
adds around it. `cache_dir_source` is one of `explicit`, `IONMC_CACHE_DIR`, `default` or `not-used`.

## Array catalogue

Each entry of `arrays` has the fields `name`, `dtype`, `shape`, `axes`, `grid` (the scoring grid name or
`null`), `units`, `definition`, `content_sha256` and `origin`. `content_sha256` is the sha256 of
`dtype.str | shape | C-order bytes`, independent of compression. `origin` has `kind`
(`grid`, `quantity`, `geometry` or `diagnostic`) and the keys that locate the array (`array`,
`grid`, `name`, `quantity`, `index`). Array order is C order `[ix, iy, iz]`; the grid origin is the
corner of voxel (0,0,0); lengths are in mm. An `.npz` member is named like the array plus `.npy`.

Names are checked against an allowlist; an array outside it is refused by the writer and by the
reader. In the patterns, `<g>` is the safe name of the scoring grid and `<q>` the safe name of the
tally (see [safe names](#files-and-safe-names)).

### Grid arrays: `grid__<g>__<name>`

| Name | Units | Definition |
|---|---|---|
| `energy_mev` | MeV | mean deposited energy per primary history |
| `energy_std_mev` | MeV | standard error of the mean of `energy_mev` (batch statistics) |
| `dose_gy` | Gy | dose per primary history (NaN where the voxel mass is 0) |
| `dose_std_gy` | Gy | standard error of the mean of `dose_gy` |
| `relative_uncertainty` | 1 | std / mean (NaN where the mean is 0) |
| `defined_mask` | bool | positive mean energy and positive voxel mass |
| `mass_g` | g | exact voxel mass |
| `n_nonzero_batches` | count | batches with a nonzero sum |
| `sparse_mask` | bool | defined voxels with fewer than B/2 nonzero batches |
| `batch_energy_mev` | MeV | per-batch energy estimates E_b / (N/B); axes `batch, x, y, z`; only with persisted batches |

### Tally arrays: `grid__<g>__q__<q>__<name>`

For every tally the quantity result is stored as `mean`, `std`, `defined_mask`, `n_nonzero` and
`rounding_bound`. With persisted batches the per-batch numerator and denominator sums are stored as
`batch_num` and `batch_den`, so that means and standard errors can be recomputed. The batch arrays of
a spectrum tally (`fluence_spectrum`) have five dimensions
(batch, x, y, z, energy bin); every other array has at most four (`MAX_NDIM`, `MAX_NDIM_SPECTRUM_BATCH`
below). Ratio quantities (`let_t`, `let_d`, `let_d_eps`, `lookup_dose_avg`) are NaN wherever
their mask is false.

### Geometry arrays: `geometry__<name>`

| Name | Units | Definition |
|---|---|---|
| `material_index` | 1 | index of the voxel material in the geometry material list |
| `density_g_cm3` | g/cm3 | per-voxel density override (only if set) |

### Diagnostic arrays: `diag__<name>`

Present only when requested with `[diagnostics]`.

| Name | Units | Definition |
|---|---|---|
| `end_position_mm` | mm | track end position |
| `end_code` | 1 | track end reason code |
| `end_direction` | 1 | track end unit direction |
| `control_residual` | 1 | control residual per track |
| `control_displacement` | 1 | control displacement per track |
| `end_energy_mev` | MeV | track end kinetic energy |
| `escape_history` | count | history index of an escaping particle |
| `escape_position_mm` | mm | escape position |
| `escape_direction` | 1 | escape unit direction |
| `escape_energy_mev` | MeV | escape kinetic energy |
| `trace_end_history` | count | traced history index |
| `trace_end_code` | 1 | end code of a traced history |
| `trace_end_energy_mev` | MeV | end energy of a traced history |

Any other name starting with `trace_` (the per-step trace arrays) is accepted as a diagnostic.

## Units

The `units` block maps each quantity class to its unit; the unit vocabulary of the arrays is that
set plus the entries `MeV/mm` and `channel-unit`.

| Class | Unit |
|---|---|
| energy | `MeV` |
| dose | `Gy` |
| length | `mm` |
| mass | `g` |
| density | `g/cm3` |
| let | `keV/um` |
| fluence | `mm^-2` |
| counts | `count` |
| dimensionless | `1` |
| mask | `bool` |

The full unit vocabulary is `MeV`, `Gy`, `g`, `mm`, `mm^-2`, `keV/um`, `MeV/mm`, `g/cm3`, `1`,
`bool`, `count` and `channel-unit`.

## Reader budgets

`results.load` and the `verify` and `report` commands refuse a bundle that exceeds a budget before any
allocation proportional to it. The defaults are the following module constants; each can be lowered or
raised with a keyword of `load` or the option of the same name of the commands.

| Constant | Value | Applies to | Option |
|---|---|---|---|
| `MAX_TOTAL_BYTES` | 8589934592 | `arrays.npz` and all arrays together (8 GiB) | `--max-total-bytes` |
| `MAX_ARRAY_BYTES` | 4294967296 | one array (4 GiB) | `--max-array-bytes` |
| `MAX_JSON_BYTES` | 16777216 | `result.json` (16 MiB) | `--max-json-bytes` |
| `MAX_CENTRAL_DIRECTORY_BYTES` | 1048576 | the zip central directory (1 MiB) | `--max-central-directory-bytes` |
| `MAX_JSON_TOKENS` | 2000000 | JSON tokens of `result.json` | `--max-json-tokens` |
| `MAX_JSON_DEPTH` | 64 | JSON nesting depth | `--max-json-depth` |
| `MAX_DERIVED_BYTES` | 67108864 | `report.md` and each `plots/*.svg`, each and together (64 MiB) | `--max-derived-bytes` |

Fixed limits (not options): `MAX_NDIM` = 4 dimensions per array, `MAX_NDIM_SPECTRUM_BATCH` = 5 for the
batch arrays of a spectrum tally, `MAX_NPY_HEADER_BYTES` = 1024 for an `.npy` header and
`MAX_INPUT_FILE_BYTES` = 268435456 for one input file read by `verify --rederive` (bounded by its recorded
size). Before the archive is opened the reader locates the end-of-central-directory record within the
last 65557 bytes and checks entry count, directory size and offsets against these budgets; each
`.npy` header must equal its catalogue entry before data are read, and decompression is streamed with
a byte budget.

## Files and safe names

`files` lists `report.md`, every `plots/<name>.svg` and `arrays.npz`, each with `sha256` and `size`.
The fixed file names are `arrays.npz` and `report.md`; plots are `plots/<safe>.svg`. Every name
is bundle-relative, single-level or `plots/<safe>.svg`: never absolute, no `..`, no other separator.
The reader checks each entry with `lstat` (a regular file, not a symbolic link) and refuses
otherwise; the writer never follows or creates symbolic links.

Safe names come from the user's grid, tally and lookup names: lower case, characters `[a-z0-9_-]` only
(other runs become `-`), at most 20 characters for a grid or tally component and at most 64 for a file,
with an index suffix on collisions. The mapping name to member is recorded in the catalogue.

## What tier 1 binds

`load` binds, from `result.json` and `arrays.npz` alone: the run and the run-level effective fields to
`requested`; the array catalogue and the tally channel plan to the requested tallies, grids and lookup
records (the per-channel `f_max` of a lookup tally only through the lookup digests); material definitions
and `mass_g` from the persisted density; the stopping identity block against the registry; every
effective field that does not depend on tables, lookups or a user stopping object, through the
reconciliation table; the accounting identities (relative tolerance 1e-9); and the array invariants
(dose from energy and mass, masks from `n_nonzero_batches`, means and standard errors from the batch
arrays with relative tolerance 1e-9).

What it does not bind is returned as `LoadedResult.binding_limits` (and printed by `ionmc verify`):
the table-derived parts of `effective` (`tables`, the multiple-scattering mapping, the scattering
length, the step bound, the scoring pieces), the per-projectile digests of the ICRU 90 model and of a
user stopping object, the `evidence_at_run` records (format-checked and persisted verbatim, never
recomputed), the table-dependent bounds of lookup tallies, the nuclear table summary except its id, and,
when `persist_batches = false`, the recomputation of means and standard errors. Those parts are protected by
the file hashes only; `ionmc verify --rederive` rebuilds them under explicit trust. Lookup tables are
excluded from the tier-2 comparison when only placeholder tables are available:
`effective.tallies.channels[kind == 'FE'].k`, `effective.tallies.channels[kind == 'FE'].bound_per_history`,
`effective.tallies.bounds.f_max`, `effective.tallies.lookups[*].content_sha256` and
`effective_extra.lookup_plan[*].f_max`.

## Writing

`save(result, path, manifest=None, persist_batches=True, overwrite=False)` writes into a temporary
sibling directory and publishes it with one atomic rename; the destination never holds a partial
bundle. An existing destination is refused unless `overwrite=True`. `save` never deletes anything:
with `overwrite=True` the new bundle is exchanged atomically with the old one, and the previous bundle
is kept as the sibling `<path>.previous-<token>` (`ionmc run --overwrite` reports it on stderr;
removing it is the user's decision). If a write fails after its temporary directory exists, that
directory is kept as `<path>.failed-<token>` and the exception carries a note naming it, which the
CLI prints. A result with `valid = false` is written only when the run allowed it. Without the
required rename primitives (Linux `renameat2`) the write fails closed.

Not a clinical claim.
