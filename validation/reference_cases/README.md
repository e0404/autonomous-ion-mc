# Reference-engine case bundles

Native inputs for the registered independent Monte Carlo engines (OpenTOPAS 4.3 / Geant4 11.4.2,
MCsquare e0404, FRED 3.76 CPU). Each directory holds `case.json` (schema checked by
`infrastructure/experiment_v3/reference.py::validate_case`) and the native input files. The
bundles must be committed before `run_reference_calculation` is called on them. No bundle is a
physics validation by itself; evidence requires the frozen criteria of the requirement ledger.

## Layout

```
topas/proton-water-150mev[-smoke]/       150 MeV p, QGSP_BIC_HP + opt4, dose 1 mm column + 2 mm 3D + LETd
topas/carbon-water-290mevu-smoke/        12C 3480 MeV total, QMD, Edep + primary-carbon Fluence (1/mm2) vs depth
topas/proton-water-150mev-fine[-seed2|-seed3]/   0.5 mm depth bins, 1e5 histories (evidence-grade candidates)
topas/proton-water-150mev-lateral[-seed2|-seed3|-smoke]/  T15 lateral case: IDD 0.5 mm + binary 3-D dose 240x240x300 (0.5 x 0.5 x 1 mm)
fred/proton-water-150mev-fine[-seed2|-seed3]/    1 x 1 x 0.5 mm voxels, 1e5 primaries (evidence-grade candidates)
mcsquare/proton-water-150mev-fine[-seed2|-seed3]/ 0.5 mm depth bins, 1e5 primaries (evidence-grade candidates)
mcsquare/proton-water-150mev[-smoke]/    150 MeV p, hand-written BDL (sigma 1.0 mm), 2 mm water CT, dose + LET
fred/proton-water-150mev[-smoke]/        150 MeV p, 1 mm water phantom, dose + LETd
fred/carbon-water-290mevu-smoke/         12C 290 MeV/u, 100 primaries: tests ion + nuclear support of the build
```

Common setup: water 120 x 120 x 300 mm, entrance face at depth 0, monoenergetic pencil on the
axis, seed 20261003, one thread. TOPAS and FRED use a zero-size, zero-divergence source. The
committed MCsquare BDL uses sigma 1.0 mm and divergence 1e-6 rad (correlation 0, energy spread 0 %)
because smaller or zero-width configurations placed all primaries outside the geometry.

## Roles

* Full cases (20000 histories, no suffix): currently classified as exploratory (see "Evidence
  status and fine cases" below); they are candidates for independent MC evidence towards suites
  S-PHYS-PROTON-EM, S-PHYS-PROTON-NUCLEAR and S-PHYS-LET once evidence-grade configurations
  (multiple seeds, >= 1e5 primaries per batch, <= 0.5 mm bins) are run.
* `-smoke` cases (100 to 200 histories): execution, timing and output-format checks only.
  The carbon smoke cases additionally test whether TOPAS (QMD) and FRED transport 12C and the
  energy convention; they provide no ion evidence.

## Engine roles established by the smoke runs

* TOPAS (Geant4): protons and ions (12C runs with QMD; primary-carbon Fluence scoring).
* MCsquare: protons only.
* FRED 3.76 CPU: protons only. 12C is recognised but aborts with "fragmentation of C12 ...
  not implemented"; `fred/carbon-water-290mevu-smoke` is kept as a documented negative result.
  `lTracking_nuc = t` enables FRED's nuclear elastic and inelastic modules for protons only.

## Materializing outputs

After a run, use the controlled tools, not shell copies:

1. `list_reference_artifacts` for the run id: exact relative paths, sizes and sha256.
2. `materialize_reference_artifacts` copies all or selected files to
   `.ionmc-cache/reference-runs/<run-id>/`.
3. Analyse locally with `ionmc.reference` (`read_topas_csv`, `read_metaimage`), which returns
   numpy arrays plus metadata. MetaImage arrays are `(nz, ny, nx)`; TOPAS arrays are `(nx, ny, nz)`.
   Raw logs stay archived; do not paste file contents into conversations.

## Items to verify on first (smoke) run

See `case.json` text and the task report: TOPAS ion BeamEnergy convention and filter parameter
names; MCsquare beam axis/isocentre placement, output file names for
LET and dose normalisation; FRED particle name for 12C, mhd spacing units, LETd output name.

## MCsquare beam model and output conventions

Sources: MCsquare master, gitlab.com/openmcsquare/MCsquare (`src/compute_beam_model.c`,
`src/compute_scoring.c`, `src/compute_simulation.c`, `BDL/BDL_default_DN.txt`).

* SpotSize is a standard deviation (sigma, mm) and Divergence a sigma (rad); the code builds the
  covariance matrix [[s^2, c s d], [c s d, d^2]] and samples with the square roots of its eigenvalues.
  Our BDL uses SpotSize 1.0 mm, Divergence 1e-6 rad, Correlation 0, EnergySpread 0 %: the beam is
  not point-like (sigma 1 mm), the smallest verified working size.
* Dose.mhd is dose per delivered proton (sum of deposited energy / simulated primaries / voxel
  volume, normalisation 1.0), not scaled by meterset or ProtonsMU. The source converts to Gy only
  for DVH with 1.602176e-19 * 1000 * N_delivered, so the file unit is inferred to be eV/g per
  proton: confirm the magnitude before use.
* Gantry 0: the beam enters at the y_max face and travels in -y, so depth = (y_max - y); the
  isocentre is the CT centre (60, 150, 60) mm.

### Beam-model diagnostics (150 MeV, 320 primaries; "outside" = all primaries generated outside)

| Variant | Run | Result |
| --- | --- | --- |
| default BDL verbatim | REF-73894625 | OK |
| spread 0 | REF-29377f19 | OK |
| 3 energies only | REF-88b26422 | OK |
| ProtonsMU in scientific notation | REF-fd4476b5 | OK |
| SpotSize 0.001 alone | REF-c1359df6 | OK |
| Divergence 1e-6 alone | REF-6447ed2d | OK |
| Correlation 0 alone | REF-75d2ccdf | FAIL (all outside) |
| point spot (0.001, 1e-6, 0) | REF-281e8414 | FAIL |
| 0.5 mm, 1e-3 rad, corr 0 | REF-bb2e900a | FAIL |
| 1.0 mm, 1e-6 rad, corr 0 | REF-62136da2 | OK (used in the cases) |

The diagnostic case directories were removed from the tree after these runs.

## Evidence runs (commit c6bc047, 20000 histories)

| Engine | Run id |
| --- | --- |
| TOPAS | REF-247e0913d2a55767d65f-432cc9ae |
| FRED | REF-ef9079a711bbbc7346a4-625cee73 |
| MCsquare | REF-e62d4ae032ee7d48700d-608a98c2 |

Reproduce the depth-dose comparison (peak depth, R90, R80, distal 80-20 width, pairwise
differences, output file hashes; no raw curves) after materializing the runs:

```
uv run python validation/scripts/reference/compare_depth_dose.py \
  --runs .ionmc-cache/reference-runs/<run-id>... --output <out>.json
```

### Evidence status and fine cases

The three 20000-history runs above were produced at commit c6bc047, where the `rationale` in the
case definitions read only "independent Monte Carlo evidence for suite
S-PHYS-PROTON-EM/-NUCLEAR/-LET" (no exploratory qualifier). The orchestrator reclassified them as
exploratory afterwards, as a conservative reclassification (single seed, 1-2 mm bins, no
uncertainty), and added the qualifier "exploratory reference; evidence-grade configuration
pending" to the current case definitions. The `inputs/case.json` hashes in the run manifests
refer to the c6bc047 files, not the current ones.
Evidence-grade runs use several committed seed variants per case (`-fine`, `-fine-seed2`,
`-fine-seed3`, and the `-lateral...` variants; seeds 20261003 / 20261004 / 20261005), 1e5 primaries
per run and depth bins of at most 0.5 mm: `topas/proton-water-150mev-fine` (ZBins 600),
`mcsquare/proton-water-150mev-fine` (CT 60 x 600 x 60 at 2 x 0.5 x 2 mm; the anisotropic spacing
was accepted by MCsquare e0404 in the archived runs) and `fred/proton-water-150mev-fine`
(120 x 120 x 600 voxels at 1 x 1 x 0.5 mm). Their batch analysis is implemented in
`validation/scripts/reference/compare_batches.py` (replicate gate on the full executed
configuration except the seed, per-engine mean / SD / SE / 95 % t-interval) and the archived results
are listed in `validation/results/reference/README.md`. `compare_depth_dose.py` (single runs) always
labels its results exploratory and computes no uncertainty; neither script accepts detached
metadata to change a status.

Library code: `ionmc.reference.runs` (`load_run`, `depth_dose`, per-engine extraction, fail closed)
and `ionmc.reference.metrics`. Metrics use bin-centre depths, argmax for the peak (resolution is
the depth bin width of the run: 1 mm TOPAS/FRED and 2 mm MCsquare for the older exploratory runs,
0.5 mm for all fine batch runs) and linear interpolation on the distal side.

## Evidence-grade batch protocol (V3-010B)

`run_reference_calculation` accepts no seed or history override (arguments: `task_id`,
`case_path`, `engine`, `gpu`, `timeout_seconds`), so the seed is varied by committed variant
directories. Each directory differs from its siblings only in the native seed line and
`case.json` `seeds`: base `-fine` / `-lateral` = 20261003, `-seed2` = 20261004, `-seed3` =
20261005 (TOPAS `i:Ts/Seed`, MCsquare `RNG_Seed`, FRED `-rseed` in `case.json` `arguments`). Every
run is reproducible from its committed bundle; the manifested `inputs/case.json` records the seed.

* Depth-dose batches (1e5 primaries, depth bins 0.5 mm): TOPAS `-fine`, MCsquare `-fine`
  (anisotropic 2 x 0.5 x 2 mm CT: if e0404 rejects it, a 1 mm isotropic fall-back is not
  evidence-grade and must be reported), FRED `-fine` (new; 120 x 120 x 600 voxels of
  1 x 1 x 0.5 mm; the 1 mm case was too coarse), each with `-seed2`, `-seed3`.
* Lateral case (T15): `topas/proton-water-150mev-lateral`. The 3-D DoseToMedium scorer covers the
  whole 120 x 120 mm field (240 x 240 bins of 0.5 mm) and 300 depth bins of 1 mm, output
  `Binary` (`dose3d.bin` + `dose3d.binheader`, Sum only, 17.28e6 bins, 138 MB if float64). The full
  field is scored, not +-30 mm: the scorer shares the `Water` component with the 0.5 mm IDD (a
  child scoring box would remove its volume from the IDD), and the halo then lies inside the grid. The frozen T15 text names 0.2 mm lateral bins;
  240 x 240 bins of 0.2 mm over the field would be 6.7x larger (about 0.9 GB per run), so 0.5 mm is used
  (Sheppard-corrected, tested at sigma = 2-6 mm; accepted by the orchestrator).
  A second identical 3-D scorer `Dose3DPrimary` (`dose3d_primary.bin`) adds the TOPAS filter
  `OnlyIncludeParticlesOfGeneration = "Primary"` (generation-0 particles, no nuclear-secondary
  dose); it is INFORMATIVE only and not the EM-only reference (the EM-only all-particle run is, see
  below), because it also excludes delta-electron dose; the full-physics all-particle scorer serves
  the comparison after nuclear physics exists. The parameter `OnlyIncludeParticlesOfGeneration =
  "Primary"` was accepted by TOPAS 4.3 in the smoke run (REF-95ef7515c835ec8c348f-e1a37fc0) and is
  used by all archived lateral runs (the binheader records the filter). Outputs
  and `compare_batches.py` report sigma_lat for both (`lateral`, `lateral_primary`). Dose3D binary layout and precision are not documented in the
  repository: `ionmc.reference.parsers.parse_topas_binary` infers float32/float64 from the file
  size and `runs.dose_3d` requires the lateral mean of the 3-D dose to reproduce the IDD (1e-4 of the
  peak bin), which fails closed on a wrong memory order. Run the `-lateral-smoke` case first.
* Batch statistics: `validation/scripts/reference/compare_batches.py` (>= 2 runs per engine,
  distinct seeds read from the manifested `inputs/case.json` and cross-checked against the
  manifested native input, equal histories >= 1e5 and bins <= 0.5 mm; otherwise it exits
  non-zero). Output status `batched` means only that the protocol holds. Results are evidence only
  after the frozen criteria of the acceptance plan have been evaluated against them; the physics list
  (QGSP_BIC_HP + opt4) and I-values (TOPAS G4_WATER 78 eV, FRED 75 eV, MCsquare not printed) label
  every result.
* Lateral estimator (`ionmc.reference.metrics.lateral_variance_1d`): per 1 mm slab at z = f R80
  (f = 0.5, 0.9; R80 from the IDD of the same run), x and y projections, second central moment
  of the dose profile minus Sheppard h^2/12, mean of the two projections, over the full field and
  a +-20 mm window. NO pedestal subtraction (Monte Carlo dose has no additive background and
  the tails are physical; negative dose is rejected). The full-field value includes the tails by
  definition; the +-20 mm value truncates a Gaussian core by a relative sigma^2 bias of
  -1.5e-5 (sigma 4 mm), -1.1e-3 (5 mm), -1.0e-2 (6 mm) (formulas in the module notes). The standard
  error is the spread across seeds.
* Observables and cases for T15. `topas/proton-water-150mev-lateral-emonly[-seed2|-seed3|-smoke]`
  is identical to the lateral cases except for the physics list
  (`sv:Ph/Default/Modules = 1 "g4em-standard_opt4"`: no hadronic, ion or decay modules); its
  ALL-particle scorer (`dose3d`) is the T15 reference for EM-only backends. The full-physics
  lateral cases serve the comparison after nuclear physics exists. The generation-filtered
  `dose3d_primary` scorer (in every lateral case) is INFORMATIVE only: it excludes the dose of EM
  secondaries (delta electrons), whereas ionmc deposits all electronic energy loss locally.
  Even the EM-only reference differs from ionmc by Geant4's delta-electron transport above the
  production cut (CutForAllParticles 0.05 mm, electron cut energy roughly 50 keV in water): Tmax at
  150 MeV is 0.35 MeV (CSDA range about 1 mm, [unverified] ESTAR value quoted from memory), the
  energy fraction in such deltas is a few percent with mean displacement well below 0.5 mm,
  bounding the extra variance by about 0.01 mm^2 against sigma^2 >= 2.9 mm^2, i.e. <= 0.4 % on
  sigma^2 and <= 0.2 % on sigma (range argument, not computed).
* The comments inside the lateral TOPAS `input.txt` files were corrected after the archived runs
  (no 'background zone'; the primary-generation scorer is informative only). Comments do not
  change the executed physics; the run manifests hash the input copies as executed (source commits
  0d0bb49 / a17f47a), which carried the old comment text.
* Replicate gate (`compare_batches.py`): within a group all runs must have identical engine
  identity (request.json engine block, runner/OS/sandbox hashes, clean commit), identical case.json
  (minus `seeds` and `rationale`, the single `-rseed <int>` pair of FRED removed; zero, several or malformed `-rseed` options are rejected), identical native input after deleting only the
  seed line(s) and identical hashes of all other manifested inputs (CT, BDL, Plan), and declared
  histories equal to the native input, FRED `-nprim` and the engine's own run summary. TOPAS
  EM-only runs form their own group `topas-emonly`.

Smoke finding (REF-95ef7515c835ec8c348f-e1a37fc0, 200 histories, 2.6 s): the filter
`OnlyIncludeParticlesOfGeneration = "Primary"` is accepted (recorded in the binheader as
"Filtered by"); binary files are float64, x fastest ("F" order), dose = energy / voxel mass, so
the lateral mean of the 3-D dose equals the IDD dose of the same slab (agreement 2e-8).
