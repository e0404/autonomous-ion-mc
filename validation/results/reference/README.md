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
the Sheppard-corrected lateral sigma (second moment of the dose profile, no pedestal
subtraction; full field and +-20 mm window) at z/R80 = 0.5 and 0.9. Replicates must share the
complete executed configuration except the seed (engine identity, case, normalized native input,
auxiliary inputs, histories; see the case README). Its status `batched` is not a grade:
results are evidence only after the frozen criteria of the acceptance plan (e.g. T15) have been
evaluated against them. `compare_depth_dose.py` still labels its single-run results exploratory.
Both scripts record git HEAD of the analysis code and a dirty flag; `--code-sha` must equal HEAD.

## Archived batch results (2026-10-04, analysis code a17f47a, clean)

Generated with `compare_batches.py` from the materialized runs listed below (1e5 primaries each,
seeds 20261003 / 20261004 / 20261005; every run passed the full-configuration replicate gate and
the manifest hash verification). Aggregates only; **not graded** against any frozen criterion —
these files are the reference side of later evaluations (V3-003B T15, V3-005 E1 and the
S-PHYS-PROTON-EM/-NUCLEAR suites), which must be run and graded against the plan, not read off here.
An earlier version of the lateral file (analysis code 0d0bb49) used a pedestal-subtracted estimator
and was superseded; its values are not comparable and are not kept.

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
- `2026-10-04-proton-150mev-water-lateral-batches.json` — two TOPAS groups with 3-D dose on
  0.5 mm lateral x 1 mm depth bins over the full +-60 mm field, each with an all-particle and a
  primary-generation scorer (the latter informative only, see plan clarification 8):
  `topas` = full physics (REF-c005bde2015f55dfa3ac-d542c739, REF-83114c1459cf799954c6-1c5d7e4d,
  REF-2e65bdfe7e29fe0d0df3-9b24464b; 6.5-12 min each) and `topas-emonly` = electromagnetic modules
  only, `g4em-standard_opt4`, no hadronic/elastic/decay processes in the TOPAS process listing
  (REF-de1e8f4142b8b17411fc-747e54b9, REF-104a823447078de5446a-a9ecef08,
  REF-ad4a98a0c3601c281445-f0a4ffb1; about 6 min each). Lateral sigma is the Sheppard-corrected
  second moment of the dose profile over the stated window, no pedestal subtraction. Values
  (mm, mean +- SE over seeds) at z/R80 = 0.5 and 0.9: full physics, all particles, +-20 mm window
  2.93 +- 0.01 and 4.42 +- 0.01 (full field 5.11 / 6.11); EM-only, all particles, +-20 mm
  1.154 +- 0.002 and 3.019 +- 0.008 (full field 1.18 / 3.12); EM-only primary-only differs from
  EM-only all-particle by about 0.0008 mm at both depths (paired over seeds, roughly six standard
  errors of the paired difference): the delta-electron contribution is detectable but negligible
  against the 3 % T15 tolerance. The large difference between full physics and EM-only is the nuclear halo (secondary
  particles) plus hadron-elastic deflection of primaries. R80: EM-only 158.81 +- 0.00 mm versus
  full physics 158.61 +- 0.01 mm. The EM-only all-particle values are the T15 reference for the
  EM-only ionmc backends; the acceptance window is +-20 mm (plan clarification 10), and the
  full-field second moment is informative only (ionmc's Gaussian MCS lacks the single-scattering
  tail, so its full-field sigma is expected to fall below the reference).
