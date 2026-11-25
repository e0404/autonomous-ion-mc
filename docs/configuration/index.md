# Configuration reference

A run is described by a TOML file of schema `ionmc.config`, version `1.0`, read with the standard
library `tomllib`. `ionmc check CONFIG.toml` loads and validates it without transport; `ionmc run
CONFIG.toml` runs it. The reading is strict: an unknown table or key, a wrong type or a missing
required key is a `ConfigError` (exit code 2) that names the key path, for example `physics.nuclear`.
The Python API takes the same fields; the TOML loader only adds the file layer (see the
[workflow](../workflow.md) and the [CLI](../cli.md)).

Documented examples (small, synthetic, used by the tests) are in `tests/ionmc/data/v3_006/`:
`proton_box.toml` (python backend, box phantom, `[output]`), `proton_box_warp_material.toml` (warp-cpu
float32, inline material, LET tally), `voxel_slab_lookup.toml` (voxel geometry from `.npy` files and a
lookup table), `carbon_box.toml` (a carbon-12 beam), `icru90_box.toml` (the ICRU 90 stopping model; needs
its cached dataset) and `nuclear_box.toml` (nuclear interactions; needs a built nuclear table, its
`@TABLE_ID@` placeholder must be replaced).

## Top level

| Key | Type | Required | Meaning |
|---|---|---|---|
| `schema` | string | yes | must be `ionmc.config` |
| `schema_version` | string | yes | must be `1.0` |
| `source` | table | yes | the beam, see `[source]` |
| `geometry` | table | yes | the phantom, see `[geometry]` |
| `scoring` | array of tables | yes | at least one scoring grid |
| `physics` | table | yes | physics options |
| `run` | table | yes | execution options |
| `materials` | array of tables | no | named or inline materials |
| `diagnostics` | table | no | optional diagnostic records |
| `tallies` | array of tables | no | extra scored quantities |
| `lookups` | array of tables | no | user lookup tables |
| `influence` | table | no | dose influence matrix request of a `[source] type = "beamlets"` run (decision 0046; keys `grid`, `threshold_rel`, `group_size`, `persist_batches`; documented in full by V3-007 step S7) |
| `output` | table | no | where and what to persist |

No environment variable overrides a physics or run field. Floats must be finite; a boolean is never accepted where a number is expected.

## [source]

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `projectile` | string or table | yes | | a registered species name (`proton`, `deuteron`, `triton`, `helium3`, `alpha`, `carbon12` or `oxygen16`; which of them can be transported is reported by `ionmc capabilities`) or an inline definition (see below) |
| `position_mm` | 3 floats | yes | | start position of the pencil beam, mm |
| `direction` | 3 floats | yes | | beam direction; normalised (recorded as a derived value) |
| `kinetic_energy_mev` | float | yes | | mean kinetic energy, MeV |
| `energy_sigma_mev` | float | no | 0.0 | Gaussian energy spread, MeV |
| `lateral_sigma_mm` | float | no | 0.0 | Gaussian lateral spread, mm |

## Inline projectile definition

An inline table in place of the registered name; it is checked against the registry.

| Key | Type | Required | Meaning |
|---|---|---|---|
| `name` | string | yes | species name |
| `symbol` | string | yes | chemical symbol |
| `z` | integer | yes | charge number |
| `a` | integer | yes | mass number |
| `mass_mev` | float | yes | mass energy, MeV |

## [geometry]

`type` selects the form. Both forms need `type`.

### Geometry type box

| Key | Type | Required | Meaning |
|---|---|---|---|
| `type` | string | yes | `box` |
| `lower_mm` | 3 floats | yes | lower corner, mm |
| `size_mm` | 3 floats | yes | edge lengths, mm |
| `material` | string | yes | a built-in material or one defined in `[[materials]]` |

A box is recorded as a derived one-voxel voxel geometry (`box-to-voxel`).

### Geometry type voxel

| Key | Type | Required | Meaning |
|---|---|---|---|
| `type` | string | yes | `voxel` |
| `origin_mm` | 3 floats | yes | corner of voxel (0,0,0), mm |
| `spacing_mm` | 3 floats | yes | voxel size, mm |
| `shape` | 3 integers | yes | `[nx, ny, nz]`, C order `[ix, iy, iz]` |
| `materials` | array of strings | yes | the material list that `material_index` refers to |
| `z_exit_mm` | float or absent | no | optional exit plane |
| `material_index` | `{file, sha256}` | yes | `.npy` file of integer indices |
| `density` | `{file, sha256}` | no | `.npy` file of per-voxel densities, g/cm3 (stored as `density_g_cm3` in the API) |

## [[materials]]

A table `name = "water"` alone selects a built-in material. Otherwise the table is an inline definition,
checked exactly as the API checks a `Material`.

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `name` | string | yes | | material name |
| `density_g_cm3` | float | yes | | density, g/cm3 |
| `mass_fractions` | table | yes | | element symbol to mass fraction, for example `{ H = 0.11, O = 0.89 }` |
| `I_eV` | float or absent | no | none (the material model computes it) | mean excitation energy, eV |
| `sternheimer` | table | no | | density-effect parameters, see below |
| `source` | string | no | empty | free-text origin of the definition |

## Sternheimer parameters

Inside `[materials.sternheimer]` (all required together): `x0`, `x1`, `cbar`, `a`, `m`.

## [[scoring]]

At least one grid; the number of grids is limited. Grid arrays are written in C order `[ix, iy, iz]`.

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `origin_mm` | 3 floats | yes | | corner of voxel (0,0,0), mm |
| `spacing_mm` | 3 floats | yes | | voxel size, mm |
| `shape` | 3 integers | yes | | voxels per axis |
| `name` | string | no | dose | grid name (appears in the array names and the report) |

## [physics]

`nuclear` has no default and `stopping_model` is required: a missing `stopping_model` is a
`ConfigError` naming `physics.stopping_model`. Possible model ids: `bethe-ionmc-v1`,
`bethe-ionmc-v1+ls-eff`, `bethe-ionmc-v1+ls-bare_mixed`, `nist-pstar-astar-water-2005` and
`icru90-water-v1` (the tabulated ones need their cached datasets, see the [workflow](../workflow.md)).
A user stopping object exists only in the API.

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `nuclear` | bool | yes | | nuclear interactions (proton source and python backend only) |
| `stopping_model` | string | yes | | registered stopping model id |
| `straggling` | bool | no | true | energy-loss straggling |
| `straggling_model` | string | no | bohr_gamma_v1 | straggling model (`bohr_gamma_v1` or `bohr_gauss_clamped_gamma_v1`) |
| `multiple_scattering` | bool | no | true | multiple Coulomb scattering |
| `mcs_model` | string | no | differential_moliere | scattering model |
| `delta_electrons` | string | no | local | delta-electron model (`local`) |
| `e_cut_mev` | float | no | 2.0 | proton cutoff energy, MeV (must stay 2.0 for ions) |
| `max_step_mm` | float | no | 1.0 | step length bound, mm |
| `max_energy_loss_fraction` | float | no | 0.02 | bound on mean energy loss per step, fraction of the kinetic energy |
| `range_alpha` | float | no | 0.2 | range step function parameter |
| `range_rho_f_mm` | float | no | 0.1 | range step function parameter, mm |
| `short_step_fraction` | float | no | 0.01 | short steps use the midpoint rule below this fraction of the residual range |
| `truncated_hinge_diagnostic` | bool | no | false | a diagnostic, not a physics model |
| `e_cut_deuteron_mev` | float | no | 4.0 | deuteron cutoff for nuclear runs, MeV |
| `nuclear_table_id` | string | no | none | nuclear table id (needed with `nuclear = true`) |
| `e_cut_ion_mev_per_u` | float | no | 2.0 | cutoff of an ion source, MeV per nucleon |
| `mcs_mapping_eval` | string | no | table | evaluation of the ion scattering mapping (`table` or `root`) |
| `end_of_range_model` | string | no | default | `default`, `point_v1` or `continuous_residual_v1` (ions only) |

## [run]

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `backend` | string | yes | | `python`, `warp-cpu` or `warp-cuda` |
| `precision` | string | yes | | `float32` or `float64` |
| `seed` | integer | yes | | 64-bit unsigned random seed |
| `n_histories` | integer | yes | | primary histories |
| `n_batches` | integer | no | 20 | interleaved batches (history h is in batch h mod B) |
| `cpu_workers` | integer | no | 1 | worker processes (python and warp-cpu only) |
| `max_steps` | integer | no | computed | per-history step bound; the computed bound is recorded |
| `allow_invalid_result` | bool | no | false | write a result that hit a transport limit, marked `valid = false` |
| `worker_timeout_s` | float | no | none | worker timeout, s |
| `memory_budget_bytes` | integer | no | 2147483648 | bound of the energy accumulators |
| `chunk_histories` | integer | no | 262144 | histories per Warp launch (power of two, at least 1024) |

## [diagnostics]

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `track_end_positions` | bool | no | false | persist the end position and reason of every track |
| `escape_records` | bool | no | false | persist records of escaping particles |
| `trace_histories` | integer | no | 0 | per-step trace of the first N histories |

## [[tallies]]

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `name` | string | yes | | tally name |
| `grid` | string | yes | | name of a `[[scoring]]` grid |
| `quantity` | string | yes | | one of `edep`, `dose`, `fluence`, `let_t`, `let_d`, `let_d_eps`, `lookup_sum`, `lookup_dose_avg`, `fluence_spectrum` |
| `species` | array of strings | no | all | restrict to species |
| `generation` | string | no | all | `all`, `primary` or `secondary` |
| `lookup` | string | no | none | name of a `[[lookups]]` table (for `lookup_sum`, `lookup_dose_avg`) |
| `energy_edges_mev_per_u` | array of floats | no | none | bin edges of a `fluence_spectrum` |
| `let_medium` | string | no | water | only `water` is supported |
| `dose_reference` | string | no | medium | only `medium` is supported |

## [[lookups]]

User lookup tables are referenced by file, never embedded in the bundle.

| Key | Type | Required | Meaning |
|---|---|---|---|
| `file` | string | yes | JSON table file; a relative path resolves against the directory of the config file |
| `sha256` | string | yes | sha256 of the file, mandatory |

## [output]

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `path` | string | no | none | bundle directory; `ionmc run --out` overrides it; a run needs one of the two |
| `persist_batches` | bool | no | true | write the batch arrays; `false` is recorded explicitly, never silent |

## Input files and the manifest

Every referenced file (`material_index`, `density`, each lookup) needs its `sha256`. A relative path
resolves against the directory of the config file (`..` is allowed). The loader records the resolved
absolute path and the hashes in the `InputManifest` (`source = "toml"`, `base_dir`, the config file and
one entry per file with the role `material_index`, `density` or `lookup:<name>`); it also records
`config_sha256`, the digest of the canonical configuration, which binds the manifest to the complete
configuration. A mismatching hash is a `ConfigError`. The manifest is persisted as `input_manifest` in
the bundle; the report names input files by basename and sha256 only.

## Requested and effective

`requested` is the canonical JSON of the configuration; `effective` is the engine's resolved
summary. The reconciliation table gives every requested field one row and one class: **equal**
(effective equals requested), **derived** (a deterministic, named rule), **unused** (not applicable
here, recorded with the reason) or **requested-only** (an execution option without an effective
counterpart). Any difference outside the table raises `UnsupportedCombinationError` before transport.
The table is `ionmc.contract.RECONCILIATION_TABLE`; each run records the outcome of every row in
`effective_extra.resolutions`.

| Row id | Requested fields | Class and rule |
|---|---|---|
| `source` | `source.projectile`, `source.position_mm`, `source.kinetic_energy_mev`, `source.energy_sigma_mev`, `source.lateral_sigma_mm` | equal `equal` |
| `direction` | `source.direction` | derived `normalise`; equal `equal` |
| `geometry-box` | `geometry` | derived `box-to-voxel` |
| `geometry-voxel` | `geometry` | equal `equal` |
| `scoring` | `scoring` | equal `equal` |
| `stopping-model` | `physics.stopping_model`, `physics.stopping` | equal `equal`; derived `ls-variant-z1-default`; derived `user-object-id` |
| `stopping-identity` | `physics.nuclear`, `tallies` | derived `deuteron-and-water-row` |
| `ls` | `physics.stopping.options.ls` | derived `ls-z1-off`; equal `equal` |
| `physics-flags` | `physics.nuclear`, `physics.straggling`, `physics.multiple_scattering`, `physics.delta_electrons` | equal `equal` |
| `straggling-model` | `physics.straggling_model` | equal `equal`; unused `straggling-off` |
| `mcs-model` | `physics.mcs_model` | equal `equal`; unused `mcs-off` |
| `mcs-mapping` | `physics.mcs_mapping_eval` | equal `equal`; derived `mcs-off-root`; unused `proton-has-no-mapping` |
| `e-cut` | `physics.e_cut_mev` | equal `equal`; unused `ion-uses-per-u` |
| `e-cut-ion` | `physics.e_cut_ion_mev_per_u` | derived `per-u-times-A`; unused `proton-uses-total` |
| `nuclear-options` | `physics.e_cut_deuteron_mev`, `physics.nuclear_table_id` | equal `equal`; unused `nuclear-off` |
| `end-of-range` | `physics.end_of_range_model` | equal `equal`; derived `default-by-z` |
| `physics-numerics` | `physics.max_step_mm`, `physics.max_energy_loss_fraction`, `physics.range_alpha`, `physics.range_rho_f_mm`, `physics.short_step_fraction`, `physics.truncated_hinge_diagnostic` | equal `equal` |
| `run` | `run.backend`, `run.precision`, `run.seed`, `run.n_histories`, `run.n_batches`, `run.cpu_workers` | equal `equal`; derived `float32-warp` |
| `max-steps` | `run.max_steps` | equal `user`; derived `computed` |
| `chunk` | `run.chunk_histories` | equal `equal`; unused `python-backend` |
| `run-only` | `run.allow_invalid_result`, `run.worker_timeout_s`, `run.memory_budget_bytes` | requested-only `execution-option` |
| `diagnostics` | `diagnostics` | equal `arrays-present` |
| `tallies` | `tallies` | equal `in-channel-plan` |
| `lookups` | `lookups` | equal `in-channel-plan` |

A feature that was requested but has no effective implementation object (tallies without a channel plan,
nuclear without a nuclear setup, lookups without lookup channels, multiple scattering without a scattering
model, straggling without a straggling model) is refused. `ionmc capabilities` prints the declared
capability matrix; what a configuration outside it does is described in the [workflow](../workflow.md).

Not a clinical claim.
