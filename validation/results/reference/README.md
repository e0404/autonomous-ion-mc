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

## Archived batch results (2026-10-04, analysis code 0d0bb49, clean)

Generated with `compare_batches.py` from the materialized runs listed below (1e5 primaries each,
seeds 20261003 / 20261004 / 20261005, source commit 0d0bb49 for every case bundle); every consumed
file was verified against its run manifest. Aggregates only; **not graded** against any frozen
criterion — these files are the reference side of later evaluations (V3-003B T15, V3-005 E1 and the
S-PHYS-PROTON-EM/-NUCLEAR suites), which must be run and graded against the plan, not read off here.

- `2026-10-04-proton-150mev-water-depth-dose-batches.json` — 0.5 mm depth bins; TOPAS fine
  (REF-b759e74cd1fcb17cad84-1499c6b5, REF-bd3d07f341ae59e10a94-dec77a96,
  REF-ed4b3dcaf979eeededeb-ae27667a; about 3 min each), MCsquare fine (REF-fe5685071ef3aedb34e5-ebdb0776,
  REF-a6c5b24b445dc41682ff-641dc0d1, REF-c9f99750f5a4b8ac2b74-40ffbcfd; 3.6 s each, anisotropic
  2 x 0.5 x 2 mm CT accepted by e0404), FRED fine (REF-51defb32216edf600549-089532cc,
  REF-5aa54bead7c692ff8b50-7ef21059, REF-1214fc77a5472063f551-b6eeaf5b; 10 s each). Per-engine
  R80 (mean +- SE over the three seeds): TOPAS 158.63 +- 0.01 mm, MCsquare 157.96 +- 0.00 mm,
  FRED 157.46 +- 0.01 mm; distal 80-20 widths 2.32 / 2.26 / 2.26 mm. The inter-engine range spread
  (1.2 mm) is far larger than the seed-to-seed uncertainty and is systematic; the different water
  I-values (Geant4 78 eV, FRED 75 eV, MCsquare unprinted) are one known contribution. No claim is
  made here about which engine is closer to the truth.
- `2026-10-04-proton-150mev-water-lateral-batches.json` — TOPAS lateral case
  (REF-c005bde2015f55dfa3ac-d542c739, REF-83114c1459cf799954c6-1c5d7e4d,
  REF-2e65bdfe7e29fe0d0df3-9b24464b; 6.5-12 min each): 3-D dose on 0.5 mm lateral x 1 mm depth bins,
  all-particle and primary-generation scorers (filter `OnlyIncludeParticlesOfGeneration = "Primary"`,
  accepted by TOPAS; primary dose is about 88 % of the total). Sheppard-corrected lateral sigma at
  z/R80 = 0.5 and 0.9, both for the +-20 mm window (core width) and the full +-60 mm field (second
  moment including the halo): e.g. primary-only, +-20 mm: 1.70 +- 0.02 mm and 3.59 +- 0.02 mm; the
  full-field values are 2-3 mm larger because 0.4-2 % of the dose lies outside +-20 mm. The
  primary-only scorer still contains Geant4's hadron-elastic deflections of primaries, which an
  EM-only ionmc run lacks (plan clarification 6).
