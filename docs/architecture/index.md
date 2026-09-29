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
| `ionmc.geometry`, `ionmc.sources`, `ionmc.scoring` | Voxel/box/slab geometries, pencil beams with spot/energy/angular spread, scoring grids with exact overlap masses and batch tallies. |
| `ionmc.transport.step_physics`, `ionmc.transport.shared` | Single-source step physics bound to Python math (reference) or Warp builtins (kernels). |
| `ionmc.transport.reference` | Float64 reference backend (per-particle loop). |
| `ionmc.config`, `ionmc.simulation`, `ionmc.results`, `ionmc.runconfig` | Capability contract, `run()`, self-describing result persistence, JSON run configuration. |
| `ionmc.data` | Content-addressed download cache with SHA-256 verification, provenance records and offline reuse. |
| `ionmc.data.nist_star` | NIST PSTAR/ASTAR (SRD 124, ICRU 49) acquisition and parsing. |

Warp CPU/CUDA backends, nuclear interactions, LET and influence matrices are added by
subsequent tasks; this page lists only implemented behavior.
