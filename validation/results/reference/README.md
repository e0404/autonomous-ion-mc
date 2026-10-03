# Reference-engine comparison results

Aggregate metrics only (no dose curves), produced by
`validation/scripts/reference/compare_depth_dose.py` from materialized runs. Every consumed file
(and `inputs/case.json`) is verified against `transfer-manifest.json` (size and sha256) before
analysis, and the file hashes recorded in the JSON are the digests of the bytes actually read.
Raw engine outputs remain archived by the controlled reference service.

## Evidence status

- `2026-10-03-proton-150mev-water-depth-dose-exploratory.json`: 150 MeV proton pencil beam in a
  120x120x300 mm water box, 20000 primaries per engine (TOPAS REF-247e0913d2a55767d65f-432cc9ae,
  FRED REF-ef9079a711bbbc7346a4-625cee73, MCsquare REF-e62d4ae032ee7d48700d-608a98c2; source
  commit c6bc047). Peak depth, R80, R90 and distal 80-20 width of the laterally integrated depth
  dose and their pairwise differences.
  The runs were produced at c6bc047 under the original case rationale (no exploratory qualifier)
  and were reclassified as exploratory afterwards (conservative reclassification); their manifested
  `inputs/case.json` are the c6bc047 files.
  **Exploratory diagnostic only**: single seed, no statistical uncertainty, depth bins of
  1-2 mm against a 2.4-2.7 mm distal falloff. It is not acceptance evidence and does not
  validate ionmc. The engines also use different water I-values (Geant4 11 G4_WATER 78 eV; FRED
  prints 75 eV; MCsquare's value is not printed), an expected contribution to the sub-2 mm range
  spread.

Evidence-grade runs use several seeded runs per engine (committed `-fine`, `-fine-seed2`,
`-fine-seed3`, `-lateral...` variants, at least 1e5 primaries and depth bins of at most 0.5 mm).
`validation/scripts/reference/compare_batches.py` groups runs by engine, reads each seed from the
manifested `inputs/case.json`, requires at least 2 runs per engine with distinct seeds (otherwise
exit non-zero, no exploratory fallback) and writes per-engine mean, sample SD, standard error and
95 % t-interval of R80, peak depth, R90, the distal 80-20 width and (TOPAS lateral case)
the Sheppard-corrected lateral sigma at z/R80 = 0.5 and 0.9. Its status `batched` is not a grade:
results are evidence only after the frozen criteria of the acceptance plan (e.g. T15) have been
evaluated against them. `compare_depth_dose.py` still labels its single-run results exploratory.
Both scripts record git HEAD of the analysis code and a dirty flag; `--code-sha` must equal HEAD.
