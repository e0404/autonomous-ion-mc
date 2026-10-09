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
topas/proton-water-{150,200}mev-idd-r20[-emonly]-seed{1,2,3}/  V3-005B row V5: IDD over a 400 x 400 mm water box, 0.5 mm bins (full / EM-only physics)
mcsquare/proton-water-{150,200}mev-idd-r20-{on,off}-seed{1,2,3}/  V3-005B row V5: same grid, nuclear interactions on / off
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

## V5 reference cases (V3-005B, acceptance row V5; Amendment 7 (b))

Absolute integral depth dose of ionmc at r = 20 cm against TOPAS (gating) and MCsquare (report-only),
150 and 200 MeV, nuclear on and off. 24 committed bundles (12 TOPAS, 12 MCsquare), 1e5 primaries,
0.5 mm depth bins, three seeds per configuration; each `case.json` carries a `v5` block (`energy_mev`,
`mode`, `nuclear`, `depth_bins`) used by `validation/scripts/reference/compare_idd_v5.py` to group runs.
Seeds are `20270000 + 1000 * row + k` (k = 1, 2, 3; distinct per configuration so that on/off and full/EM-only
runs are statistically independent); row = TOPAS full 150 / 200 = 1 / 2, TOPAS EM-only = 3 / 4,
MCsquare on = 5 / 6, off = 7 / 8 (e.g. TOPAS full 150 MeV seeds 20271001-20271003; MCsquare off 200 MeV
20278001-20278003). Directories: `topas/proton-water-{E}mev-idd-r20-seed{k}` (full: opt4 + QGSP_BIC_HP + HP
elastic + stopping + binary-cascade ions + decay, as the V3-010B cases), `topas/proton-water-{E}mev-idd-r20-emonly-seed{k}`
(`g4em-standard_opt4` only), `mcsquare/proton-water-{E}mev-idd-r20-{on,off}-seed{k}`.

| Item | TOPAS | MCsquare |
| --- | --- | --- |
| Geometry | G4_WATER 400 x 400 x 320 mm (150 MeV) / 330 mm (200 MeV), ZBins 640 / 660 | CT 5 x 640 (660) x 5 voxels of 80 x 0.5 x 80 mm (CT.raw 64-66 kB zeros) |
| Lateral extent | one lateral bin over the whole 400 x 400 mm box (contains r = 20 cm; corners r > 200 mm are summed, negligible for a pencil beam) | same (5 x 5 coarse voxels summed) |
| Source | zero-size pencil, BeamEnergySpread 0 (as V3-010B) | BDL sigma 1.0 mm, divergence 1e-6 rad (as V3-010B; 190/200/210 MeV rows for 200 MeV) |
| Nuclear on/off | full vs EM-only module list | `Simulate_Nuclear_Interactions` and the three `Simulate_Secondary_*` flags True / False |
| Native unit | Gy per run (1e5 histories) | file value, inferred eV/g per primary |
| IDD conversion | MeV = Gy x 0.08 kg / 1.602176634e-13 J/MeV per bin; IDD = MeV / (1e5 x 1 g/cm3 x 0.05 cm) | IDD = sum_xz (value x 1e-6 MeV/g) x (8 cm x 8 cm) (density cancels) |

Output unit of the comparison: MeV/(g/cm^2)/primary. The script cuts every curve to the depth range
of the ionmc grid (1.1 R(CSDA)), refuses an engine run whose in-grid total is outside 0.85-1.02 of the
beam energy (unit or geometry error), requires >= 3 runs with distinct seeds per (engine, energy, nuclear) group, and
verifies each ionmc partial against its `content_sha256`. The comparison is:

```
uv run python validation/scripts/reference/compare_idd_v5.py --ionmc-dir <partials dir> \
  --topas-runs <12 TOPAS run dirs> --mcsquare-runs <12 MCsquare run dirs> --output <out>.json
```

Unverified until the engines run: (a) MCsquare e0404 accepting 80 mm lateral voxels and the
`Simulate_Nuclear_Interactions False` switch (fallback: 20 mm voxels, 21 x N x 21, 1 MB each; if nuclear cannot be switched
off, commit only the "on" cases); (b) the eV/g unit of `Dose.mhd` (the in-grid total check
catches a wrong factor); (c) the lateral-voxel invariance of the MCsquare depth-dose.

### Exploratory physics-attribution cases (V3-005B C18b; not evidence)

Four single-seed TOPAS cases at 150 MeV (1e5 primaries, geometry, source, 0.5 mm IDD scorer, cuts and
`NumberOfThreads` identical to `proton-water-150mev-idd-r20-seed1`) probe the V5 peak discrepancy seen in the
C18a dry run: the nuclear-on/off IDD ratio agrees between TOPAS and ionmc up to about 6 mm before the Bragg
peak, but TOPAS peak(full)/peak(EM-only) = 0.769 against ionmc 0.872 (TOPAS full-physics peak 0.5 mm proximal
of the EM-only peak, R80 0.19 mm shorter). ionmc has no hadronic elastic scattering (p-p and p-nucleus; decision
0041 limitation 7, V3-005C). They differ from the full case only in `Ph/Default/Modules` (and the seed). Each
`case.json` has `"role": "exploratory"` and an `exploratory` block instead of a `v5` block, so
`compare_idd_v5.py` never groups them with the evidence runs. They carry no validation claim and no
batch statistics (one seed each).

| Case (`topas/proton-water-150mev-idd-r20-...`) | Seed | Modules | Purpose |
| --- | --- | --- | --- |
| `x1-no-hadron-elastic` | 20279001 | opt4, QGSP_BIC_HP, stopping, ion-binarycascade, decay | full minus `g4h-elastic_HP`: is the peak deficit hadronic elastic? |
| `x2-emonly-plus-elastic` | 20279002 | opt4, `g4h-elastic_HP` | EM-only plus elastic without inelastic |
| `x3-no-ion-physics` | 20279003 | opt4, QGSP_BIC_HP, elastic_HP, decay | full minus `g4ion-binarycascade` and `g4stopping` (control) |
| `x4-no-hp-neutrons` | 20279004 | opt4, `g4h-phy_QGSP_BIC`, elastic_HP, stopping, ion-binarycascade, decay | control; the full and EM-only inputs already share every EM setting (same `g4em-standard_opt4`, `CutForAllParticles = 0.05 mm`, no step limits), so X4 swaps `QGSP_BIC_HP` for `QGSP_BIC` (no HP neutrons) |

**Outcomes (exploratory, not evidence).** Derived locally from `work/idd_dose.csv` of each materialized run with the
`compare_idd_v5.py` loader and `curve_metrics` (0.5 mm IDD in MeV cm⁻¹ per primary; plateau = mean over 20–60 mm;
peak depth is the bin centre; total = integral over the 320 mm grid, MeV per primary; R80 distal, interpolated).
Baselines are the committed frozen runs `proton-water-150mev-idd-r20[-emonly]-seed{1,2,3}` (REF-cec79def…, REF-e4f252d8…,
REF-f306bdc1… full; REF-7072862d…, REF-278f82f0…, REF-b9387263… EM-only), seed-averaged curves.

| Case (run id) | Peak | Peak depth (mm) | Plateau | Peak/plateau | R80 (mm) | Total (MeV) | Peak / EM-only | Plateau / EM-only | Total / EM-only |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| EM-only (3 seeds) | 37.371 | 157.25 | 6.0891 | 6.137 | 158.818 | 149.9996 | 1 | 1 | 1 |
| full (3 seeds) | 28.732 | 156.75 | 6.5271 | 4.402 | 158.636 | 144.905 | 0.7688 | 1.0719 | 0.9660 |
| X1 REF-e1fcaab4cb3a1d312991-36cabeda | 31.899 | 157.25 | 6.4270 | 4.963 | 158.776 | 144.390 | 0.8536 | 1.0555 | 0.9626 |
| X2 REF-1a17d04f66bf3f4861be-e248591a | 33.975 | 156.75 | 6.1922 | 5.487 | 158.608 | 150.002 | 0.9091 | 1.0169 | 1.0000 |
| X3 REF-fb22060d63a6943329fc-ded17160 | 28.699 | 156.75 | 6.5302 | 4.395 | 158.640 | 144.906 | 0.7679 | 1.0724 | 0.9660 |
| X4 REF-52f701f35ff7d797669e-b6824a24 | 28.778 | 156.75 | 6.5244 | 4.411 | 158.625 | 144.938 | 0.7701 | 1.0715 | 0.9663 |

Per-seed spread of the baselines (sample standard deviation, n = 3): full peak 0.106 (per-seed peak/EM-only 0.7655,
0.7706, 0.7703), plateau 0.010, peak/plateau 0.017, R80 0.003 mm, total 0.028 MeV; EM-only peak 0.061, plateau 0.0013,
peak/plateau 0.011, R80 0.015 mm, total 0.0007 MeV. Peak depth is 156.75 mm (full) and 157.25 mm (EM-only) in every
seed, one bin apart. Reading: removing hadronic elastic (X1) restores 0.085 of the 0.231 peak deficit of the full
case, and elastic alone without inelastic (X2) lowers the EM-only peak by 0.091; the product X1×X2 = 0.776 is close to
the full value 0.769 (multiplicative reading only, X1 and X2 are single seeds, so the seed spread of about 0.003 in
the ratio applies). X3 and X4 are indistinguishable from full within that spread. Elastic scattering is thus a large
but not the only contributor to the peak deficit; the attribution is to be tested by V3-005C.
