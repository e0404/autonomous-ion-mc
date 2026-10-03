# Reference-engine comparison results

Aggregate metrics only (no dose curves). Produced by
`validation/scripts/reference/compare_depth_dose.py` at clean task SHA `67d59f8ba74a53493b0cb453d25c917ea0bdcf3b`
from the materialized runs listed inside each file (run ids, engine, source
SHA of the committed case bundle, histories, output-file SHA-256). The raw
engine outputs remain archived by the controlled reference service.

- `2026-10-03-proton-150mev-water-depth-dose.json`: 150 MeV proton pencil
  beam in a 120x120x300 mm water box, 20000 primaries per engine (TOPAS
  REF-247e0913d2a55767d65f-432cc9ae, FRED REF-ef9079a711bbbc7346a4-625cee73,
  MCsquare REF-e62d4ae032ee7d48700d-608a98c2): peak depth, R80, R90 and
  distal 80-20 width of the laterally integrated depth dose and their pairwise
  differences. These are independent-Monte-Carlo reference values for later
  evaluation of ionmc; they are not validation of ionmc by themselves. Note
  the engines use different water I-values (Geant4 11 G4_WATER 78 eV; FRED
  prints 75 eV in its log; MCsquare's value is not printed), which is one
  expected source of the sub-2 mm range spread.
