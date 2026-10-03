# Reference-engine case bundles

Native inputs for the registered independent Monte Carlo engines (OpenTOPAS 4.3 / Geant4 11.4.2,
MCsquare e0404, FRED 3.76 CPU). Each directory holds `case.json` (schema checked by
`infrastructure/experiment_v3/reference.py::validate_case`) and the native input files. The
bundles must be committed before `run_reference_calculation` is called on them. No bundle is a
physics validation by itself; evidence requires the frozen criteria of the requirement ledger.

## Layout

```
topas/proton-water-150mev[-smoke]/       150 MeV p, QGSP_BIC_HP + opt4, dose 1 mm column + 2 mm 3D + LETd
topas/carbon-water-290mevu-smoke/        12C 3480 MeV total, QMD, Edep + surviving-12C surface counts
mcsquare/proton-water-150mev[-smoke]/    150 MeV p, hand-written zero-width BDL, 2 mm water CT, dose + LET
fred/proton-water-150mev[-smoke]/        150 MeV p, 1 mm water phantom, dose + LETd
fred/carbon-water-290mevu-smoke/         12C 290 MeV/u, 100 primaries: tests ion + nuclear support of the build
```

Common setup: water 120 x 120 x 300 mm, entrance face at depth 0, monoenergetic zero-size
zero-divergence pencil on the axis, seed 20261003, one thread.

## Roles

* Evidence cases (20000 histories, no suffix): independent MC evidence for suites
  S-PHYS-PROTON-EM, S-PHYS-PROTON-NUCLEAR and S-PHYS-PROTON-LET.
* `-smoke` cases (100 to 200 histories): execution, timing and output-format checks only.
  The carbon smoke cases additionally test whether TOPAS (QMD) and FRED transport 12C and the
  energy convention; they provide no ion evidence.

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
names; MCsquare beam axis/isocentre placement, zero-width BDL acceptance, output file names for
LET and dose normalisation; FRED particle name for 12C, mhd spacing units, LETd output name.
