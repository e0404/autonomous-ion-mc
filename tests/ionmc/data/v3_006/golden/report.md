# ionmc result report

## Identity

- schema: ionmc.result 1.0 (kind: transport)
- code version: 0.1.0.dev0
- package tree object id: 0123456789abcdef0123456789abcdef01234567
- commit: unknown
- dirty: unknown
- code identity status: unknown (synthetic golden fixture: not produced by a code version)

## Validity

**Result status: valid** (no transport-limit counter is nonzero).

## Requested versus effective configuration

Resolution status: computed: contract.reconcile(requested, effective). Classes: equal, derived, unused, requested-only (decision 0045 section 4).

| row | class | rule | requested | effective | reason |
|---|---|---|---|---|---|
| source | equal | equal | source.projectile = proton; source.position_mm = [0, -12, 2]; source.kinetic_energy_mev = 20; source.energy_sigma_mev = 0; source.lateral_sigma_mm = 0 | /source/projectile = proton; /source/z = 1; /source/a = 1; /source/mass_mev = 938.272; /source/position_mm = [0, -12, 2]; /source/kinetic_energy_mev = 20; /source/energy_sigma_mev = 0; /source/lateral_sigma_mm = 0 | - |
| direction | derived | normalise | source.direction = [0.2, 0.95, 0.1] | /unit_direction = [0.204926, 0.973399, 0.102463] | - |
| geometry-box | derived | box-to-voxel | geometry = {lower_mm: [-2, -8, 0], material: water, size_mm: [4, 14, 4], type:... | /geometry/origin_mm = [-2, -8, 0]; /geometry/spacing_mm = [4, 14, 4]; /geometry/shape = [1, 1, 1]; /geometry/materials = [water]; /geometry/z_exit_mm = undefined | - |
| scoring | equal | equal | scoring = [{name: dose, origin_mm: [-2, -8, 0], shape: [2, 8, 2], spacing_mm:... | /scoring = [{name: dose, origin_mm: [-2, -8, 0], shape: [2, 8, 2], spacing_mm:... | - |
| stopping-model | equal | equal | physics.stopping_model = bethe-ionmc-v1; physics.stopping = not given | /stopping/model_id_requested = bethe-ionmc-v1; /stopping/model_id_effective = bethe-ionmc-v1 | - |
| stopping-identity | derived | deuteron-and-water-row | physics.nuclear = false; tallies = [{dose_reference: medium, energy_edges_mev_per_u: undefined, genera... | /stopping/species_transported = [proton]; /tables/materials = [{I_eV_effective: 78, I_eV_requested: 78, content_sha256: undefined...; /tables/water_row = {I_eV_effective: 78, content_sha256: undefined, dataset_id: undefin... | - |
| ls | derived | ls-z1-off | physics.stopping.options.ls = not given | /physics/ls_requested = off; /physics/ls_effective = off | - |
| physics-flags | equal | equal | physics.nuclear = false; physics.straggling = true; physics.multiple_scattering = true; physics.delta_electrons = local | /physics/nuclear = false; /physics/straggling = true; /physics/multiple_scattering = true; /physics/delta_electrons = local | - |
| straggling-model | equal | equal | physics.straggling_model = bohr_gamma_v1 | /physics/straggling_model = bohr_gamma_v1 | - |
| mcs-model | equal | equal | physics.mcs_model = differential_moliere | /physics/mcs_model = differential_moliere | - |
| mcs-mapping | unused | proton-has-no-mapping | physics.mcs_mapping_eval = table | - | - |
| e-cut | equal | equal | physics.e_cut_mev = 2 | /physics/e_cut_effective_mev = 2 | - |
| e-cut-ion | unused | proton-uses-total | physics.e_cut_ion_mev_per_u = 2 | - | - |
| nuclear-options | unused | nuclear-off | physics.e_cut_deuteron_mev = 4; physics.nuclear_table_id = undefined | - | - |
| end-of-range | derived | default-by-z | physics.end_of_range_model = default | /physics/end_of_range_model_effective = point_v1 | - |
| physics-numerics | equal | equal | physics.max_step_mm = 1; physics.max_energy_loss_fraction = 0.02; physics.range_alpha = 0.2; physics.range_rho_f_mm = 0.1; physics.short_step_fraction = 0.01; physics.truncated_hinge_diagnostic = false | /physics/max_step_mm = 1; /physics/max_energy_loss_fraction = 0.02; /physics/range_alpha = 0.2; /physics/range_rho_f_mm = 0.1; /physics/short_step_fraction = 0.01; /physics/truncated_hinge_diagnostic = false | - |
| run | equal | equal | run.backend = python; run.precision = float64; run.seed = 11; run.n_histories = 8; run.n_batches = 4; run.cpu_workers = 1 | /backend = python; /precision = float64; /seed = 11; /n_histories = 8; /n_batches = 4; /cpu_workers = 1; /production = false | - |
| max-steps | derived | computed | run.max_steps = undefined | /max_steps = 952; /max_steps_origin = computed | - |
| chunk | unused | python-backend | run.chunk_histories = 262144 | /chunk_histories = 262144 | - |
| run-only | requested-only | execution-option | run.allow_invalid_result = false; run.worker_timeout_s = undefined; run.memory_budget_bytes = 2147483648 | - | - |
| diagnostics | equal | arrays-present | diagnostics = {escape_records: false, trace_histories: 0, track_end_positions: fa... | - | - |
| tallies | equal | in-channel-plan | tallies = [{dose_reference: medium, energy_edges_mev_per_u: undefined, genera... | /tallies/channels = [{bound_per_history: 10952.9, class_mask: 1, generation: [0, 214748...; /tallies/quantities = [{definition: sum(l S^2)/sum(l S), S in water, denominator: 1, kind... | - |
| lookups | equal | in-channel-plan | lookups = [] | - | - |

## Source, geometry and materials

- source: a = 1, energy_sigma_mev = 0, kinetic_energy_mev = 20, lateral_sigma_mm = 0, mass_mev = 938.272, position_mm = [0, -12, 2], projectile = proton, z = 1
- source unit direction: [0.204926, 0.973399, 0.102463]
- geometry: origin [-2, -8, 0] mm, spacing [4, 14, 4] mm, shape [1, 1, 1], materials [water], exit plane z = undefined
- grid convention: origin = corner of voxel (0,0,0); C order [ix, iy, iz]; lengths in mm (decision 0037)
- profile axis: 1 (y), the grid axis of the largest absolute component of the source unit direction
- scoring grid dose: origin [-2, -8, 0] mm, spacing [2, 2, 2] mm, shape [2, 8, 2]

Input files (basename and sha256 only):

- manifest source: api; configuration sha256: 00a6403c912182ce1b486032cc3ec0df38fcd7b6864c63adcee776b7241990d0
- no referenced input files

Materials:

| name | density g/cm3 | I eV requested / effective | mass fractions | source |
|---|---|---|---|---|
| water | 1.0 | 78 / 78 | {H: 0.111898, O: 0.888102} | G4NistMaterialBuilder v11.4.2: G4_WATER; I and Sternheimer parameters as cited in module docstring |

## Physics and stopping model

- physics options: delta_electrons = local, e_cut_effective_mev = 2, e_cut_ion_mev_per_u = 2, e_cut_mev = 2, end_of_range_model = default, end_of_range_model_effective = point_v1, ls_effective = off, ls_requested = off, max_energy_loss_fraction = 0.02, max_rejection_attempts = 64, max_step_mm = 1, mcs_model = differential_moliere, multiple_scattering = true, nuclear = false, range_alpha = 0.2, range_rho_f_mm = 0.1, scattering_E_s_mev = 15, short_step_fraction = 0.01, straggling = true, straggling_model = bohr_gamma_v1, truncated_hinge_diagnostic = false
- stopping model: requested bethe-ionmc-v1, effective bethe-ionmc-v1, version 1, spec sha256 3aa2d5b987c77afb228f4c015e579be7c9935fbbf1cfdb6accda44ccd9b9cc26
- stopping kind: analytic; interpolation: analytic on a log grid; range exact-loglog-quadrature-v1; range construction: exact-loglog-quadrature-v1
- straggling: bohr_gamma_v1; multiple scattering: differential_moliere; end of range: point_v1; delta electrons: local; nuclear table: undefined

Stopping-evidence records at run time (verbatim, as persisted; never recomputed):

```evidence
{"graded_at_sha": null, "plan": "validation/plans/v3-008-acceptance.md", "row": null, "spec_sha256": "3aa2d5b987c77afb228f4c015e579be7c9935fbbf1cfdb6accda44ccd9b9cc26", "status": "as-graded:0043"}
```

## Data identities and roles

- stopping data used: none (analytic model)
- stopping data content sha256: undefined
- table sha256: 20d74f7cab0b7411b227dbbf303f294feedf73f769ea09dc245187c331c5dd19; table stopping source: bethe
- table material water (proton): source bethe, dataset undefined, content sha256 undefined

## Run and random numbers

- backend: python, device: cpu, precision: float64, production: false
- seed: 11; histories: 8 in 4 batches (2 per batch); cpu workers: 1; chunk histories: 262144
- random numbers: counter = (history, genealogy id, block, purpose), generator = philox4x32-10, key = seed low/high 32 bits, u01 = ((w >> 8) + 0.5) * 2**-24
- timings are recorded in result.json only (not deterministic) and are not printed here

## Statistics definition

- normalisation: per primary history
- uncertainty: batch standard error of the mean over B batches (history h in batch h mod B); ratio quantities by the delta method (decision 0040)
- defined_mask: positive mean energy and positive voxel mass
- sparse_mask: defined voxels with fewer than B/2 nonzero batches
- batch arrays: {status: persisted}
- quantity batch arrays: {reason: the result carries no channel accumulators, status: not-persisted}

## Per-grid summary

### Grid dose

- shape: [2, 8, 2]; valid: true; low batch count: true
- defined voxels: 27; sparse voxels: 0
- total energy per history: 29.0659 MeV
- maximum dose per history: 4.70238e-08 Gy
- maximum relative uncertainty (defined voxels): 0.0876863 1
- voxel of maximum dose (indices): [0, 5, 0]
- tally edep_excluded_from_let (edep): defined voxels 32 of 32
- tally let (let_d): defined voxels 28 of 32
- tally scoring_pieces (count): defined voxels 32 of 32
- tally scoring_pieces_local (count): defined voxels 32 of 32

## Energy balance and counters

- closure_residual_mev: 0 MeV
- cutoff_mev: 0.5 MeV
- escaped_mev: 0.25 MeV
- in_grid_mev: [232.528] MeV
- initial_mev: 233.778 MeV
- outside_mev: [1] MeV
- quantization_mev: [] MeV
- step_deposit_mev: 233.028 MeV
- truncated_mev: 0 MeV
- unaccounted_mev: 0 MeV
- counter accumulator_overflow: 0
- counter energy_inversion: 0
- counter genealogy_overflow: 0
- counter queue_overflow: 0
- counter scoring_pieces_overflow: 0
- counter source_energy_out_of_range: 0
- counter stall: 0
- counter step_truncation: 0
- counter straggling_rejection: 0
- channel counter lookup_out_of_domain: 0
- channel counter path_bound_exceeded: 0

## Missing items

- quantity batch arrays: {reason: the result carries no channel accumulators, status: not-persisted}
- input files and configuration file: constructed in process (unavailable)
- stopping, range, LET-water and nuclear table values: identities only; values are never persisted

Not a clinical claim.
