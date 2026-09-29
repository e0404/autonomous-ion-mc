# Architecture

The intended architecture is recorded in decision
[0038](../generated/decisions/0038-v2-architecture-and-backend-strategy.md):
shared physics functions executed by a Python reference path and by Warp
CPU/CUDA generic kernels, tables that separate transport from data sources,
and planning-oriented batched execution.

## Implemented so far

| Module | Status |
|---|---|
| `ionmc.provenance` | Exact code identity (live Git or build-time stamp). |
| `ionmc.cli` | `ionmc version`, `ionmc data fetch/list/cache-dir`, `ionmc stopping`. |
| `ionmc.species`, `ionmc.materials` | Species (nuclear masses, MeV/u conventions) and materials (composition, density, I values). |
| `ionmc.physics.stopping`, `ionmc.physics.tables` | Corrected Bethe layer and blended stopping/range tables with provenance (decision 0040). |
| `ionmc.data.icru90` | Shipped ICRU 90 water tables (p, He, C). |
| `ionmc.data` | Content-addressed download cache with SHA-256 verification, provenance records and offline reuse. |
| `ionmc.data.nist_star` | NIST PSTAR/ASTAR (SRD 124, ICRU 49) acquisition and parsing. |

Transport, geometry, sources, scoring, data and results modules are added by
subsequent tasks; this page lists only implemented behavior.
