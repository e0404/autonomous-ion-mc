# Stopping-power comparison results

Output of `validation/scripts/stopping/compare_nist.py`, run with the cached NIST
PSTAR/ASTAR water tables and the ICRU 90 water arrays (Geant4 source). The JSON contains
aggregates only: per comparison the number of points, the energy range and the maximum and
RMS relative deviation for E >= 2 and E >= 10 MeV/u; for CSDA ranges the maximum absolute
relative deviation (Bethe 75 eV vs NIST, Bethe 78 eV vs integrated ICRU 90) and the
model-only Bethe 78 - 75 eV range shift per energy. It contains no NIST or ICRU 90 table
values, per-energy deviations or ratios (field `redistribution_note`), because NIST SRD 124
data must not be reproducible from Git. Dataset identities and hashes are recorded under
`provenance`. The result file is regenerated at a clean SHA by the orchestrator.

Reproduce: `ionmc data fetch` the three datasets, then
`uv run python validation/scripts/stopping/compare_nist.py --output-dir OUT --cache-dir CACHE --code-sha $(git rev-parse HEAD)`.

Archived run: produced at clean task SHA `b634e2841a9379339b3b60cd100294904ea01829`.

**Dated note (2026-11-15, orchestrator; scientific review SR-1 finding 8).** `v3-008-d3-range-attribution.json` stored a per-case ratio of the directly integrated ICRU 90 alpha range to the ionmc range next to the ionmc absolute ranges, from which the ICRU 90 ranges could be reconstructed. The field is removed (`redaction_note` in the file); the D-3 attribution (construction difference ≤ 0.08 %, the step-regime terms) is unchanged and stays in the markdown record. `d3_range_attribution.py` no longer writes the field.
