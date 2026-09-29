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
| `ionmc.cli` | `ionmc version`. |

Transport, geometry, sources, scoring, data and results modules are added by
subsequent tasks; this page lists only implemented behavior.
