# V3 program plan and task sequence

This page is the lead orchestrator's living plan for the v3 autonomous run. It
records intended milestones and task IDs so that progress, checkpoints and
scope can be traced in Git. It describes *planned* work; implemented behaviour
is documented in the architecture, physics and validation sections only once
it exists.

Governing inputs: `EXPERIMENT.md`, `REQUIREMENTS.md`,
`experiment/v3/{PROTOCOL,REQUIREMENTS,AGENTS}.md`, `experiment/v3/performance.json`
and `experiment/v3/requirements-index.json`. The requirement ledger that maps
every MUST requirement to planned evidence is `validation/requirement-ledger.md`.

## Operating rules applied to every task

1. One isolated task worktree per `V3-nnn` task, branch preserved after merge.
2. Tests and local validation run before commit; commit through the controlled
   commit tool; mandatory independent Codex review of the exact commit; fix,
   recommit, re-review until passed.
3. Exact-SHA local validation record, push, pull request, lightweight CI,
   squash merge into `v3/develop`, retire worktree.
4. Scientific evidence is recorded separately from review approval. Reference
   engine execution success is never physics validation.
5. Acceptance observables and tolerances are frozen (committed) before the
   corresponding qualification results are produced.

## Milestones

| Milestone | Content | Exit criterion |
|---|---|---|
| M0 Scaffold | Installable `ionmc` package, tooling, CI, plan, ledger, baseline architecture decision | V3-001 merged |
| M1 Proton core | External data layer, materials, geometry, scoring grids, proton EM transport on all three backends, uncertainty, LET, nuclear interactions, persisted outputs, CLI, influence matrices | Proton workflow runs end-to-end on python / warp-cpu / warp-cuda with analytic and tabulated validation |
| M2 Heavy ions | He/C/O stopping and scattering, reaction cross-sections, projectile fragmentation with transported fragments, species-resolved dose and mixed-field LET, oxygen supported domain | Ion-specific range/attenuation/fragment evidence collected |
| M3 Independent evidence | Native TOPAS/MCsquare/FRED case bundles, measured datasets with provenance, comparison tooling, numerical falsification suite | Evidence records for every physics suite; contradictions preserved |
| M4 Calibration and freeze | Performance harness, reference workload measurements, reviewed release plan with frozen targets and acceptance criteria | `validation/release-plan.json` frozen at a committed SHA |
| M5 Optimisation and qualification | Performance work within frozen targets, qualification runner, full exact-SHA qualification, promotion to `v3/main`, tag and GitHub release | `release_ready` true and promotion complete |

## Task sequence (planned IDs; may be refined, never renumbered after use)

| Task | Scope | Depends on |
|---|---|---|
| V3-001 | Package scaffold, uv tooling, CI update, program plan, requirement ledger, decision 0037 (architecture and conventions baseline) | – |
| V3-002 | External data layer: dataset registry, NIST PSTAR/ASTAR acquisition with checksum and provenance, cache location control, offline reuse, material definitions, physics table construction (analytic Bethe layer and tabulated layer behind one interface) | V3-001 |
| V3-003 | Geometry (homogeneous and voxel grids), independent scoring grid, sources with arbitrary direction, proton electromagnetic condensed-history transport (CSDA, straggling, multiple scattering) as shared Warp functions; reference, Warp CPU and Warp CUDA backends; batch uncertainty; dose and energy deposition scoring; backend parity tests | V3-002 |
| V3-004 | Extensible scoring: track- and dose-averaged LET, fluence, species-resolved tallies, external biological lookup tallies with provenance | V3-003 |
| V3-005 | Proton nuclear interactions: attenuation, secondary protons/deuterons/alphas, local deposition and escape accounting, conservation tests | V3-003 |
| V3-006 | Self-describing persisted results, human-readable reports and SVG plots, configuration files, `ionmc run` CLI, fail-closed capability contract with requested/effective configuration | V3-004, V3-005 |
| V3-007 | Beamlet-resolved transport and sparse dose influence matrices with per-beamlet uncertainty and persisted output | V3-006 |
| V3-008 | Helium, carbon and oxygen electromagnetic transport: effective-charge stopping, ion straggling and scattering, ion-specific tabulated range evidence | V3-003 |
| V3-009 | Nuclear reaction cross-sections and projectile fragmentation for He/C/O with transported secondary ions, distal fragment dose and mixed-field LET; the oxygen supported domain is declared explicitly and every claimed oxygen capability (range, attenuation, fragment build-up, distal dose, mixed-field LET) receives its own ion-specific evidence in suite `S-PHYS-OXYGEN`; unsupported oxygen combinations fail closed | V3-005, V3-008 |
| V3-010 | Reference evidence: native TOPAS/MCsquare/FRED case bundles and runs, measured data acquisition, parsers, comparison metrics (range metrics, gamma, z-tests), evidence records | V3-006, V3-009 |
| V3-011 | Adversarial numerical falsification suite and uncertainty coverage tests | V3-007, V3-009 |
| V3-012 | Performance measurement harness, calibration measurements of the three workload families, frozen release plan and targets (reviewed before merge) | V3-010, V3-011 |
| V3-013+ | Optimisation within frozen targets; qualification runner producing the exact-SHA report; full qualification; promotion | V3-012 |

Research supporting these tasks is delegated to scientific agents; their reports
are preserved under `docs/research/` once used for a decision.

## Checkpoints

Checkpoints are appended here as milestones complete (date, task, develop SHA,
summary).

- 2026-10-03, M0, V3-001 merged as `c4e63dd4` (PR #72): package scaffold,
  tooling, CI, program plan, requirement ledger, decision 0037, research
  reports under `docs/research/`, reproducible RNG and Warp architecture
  measurement archive. Ten independent review rounds were needed; the
  recurring theme was that every quoted number must be reproducible from
  committed scripts and archived raw outputs at a clean SHA.
