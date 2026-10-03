# Stopping-power comparison results

Output of `validation/scripts/stopping/compare_nist.py` executed at clean task
SHA `7a2d817085a1094bbb8cc7574391dd8d77c8fb49` (V3-002) with the cached NIST PSTAR/ASTAR water tables and the
ICRU 90 water arrays; dataset identities and hashes are recorded inside the
JSON under `provenance`. The NIST tables themselves are not redistributed.
Decision 0038 copies its validation-outcome tables from this file.

Reproduce: `ionmc data fetch` the three datasets, then
`uv run python validation/scripts/stopping/compare_nist.py --output-dir OUT --cache-dir CACHE --code-sha $(git rev-parse HEAD)`.
