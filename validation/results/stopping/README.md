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
