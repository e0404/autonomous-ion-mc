# Requirement ledger (v3)

Every MUST requirement listed in `experiment/v3/requirements-index.json` is mapped
to the planned qualification suite(s) and the evidence categories the lead
intends to use. Statuses are updated as tasks merge; nothing here is evidence.
The frozen `validation/release-plan.json` supersedes the suite mapping once it
exists. Evidence categories follow `experiment/v3/PROTOCOL.md`: measured,
ion-specific-tabulated, independent-monte-carlo, independent-theory,
related-model, backend-parity, self-consistency.

## Planned suites

| Suite | Category and content |
|---|---|
| `S-PHYS-PROTON-EM` | physics: proton stopping, range, straggling, lateral scattering |
| `S-PHYS-PROTON-NUCLEAR` | physics: proton nuclear attenuation and secondaries |
| `S-PHYS-LET` | physics: track- and dose-averaged LET including mixed fields |
| `S-PHYS-HELIUM` | physics: helium range and fragmentation |
| `S-PHYS-CARBON` | physics: carbon range, attenuation, fragment build-up, distal dose |
| `S-PHYS-OXYGEN` | physics: oxygen range and attenuation in its supported domain |
| `S-NUM-FALSIFICATION` | numerical: step/grid/precision/orientation/seed probes |
| `S-NUM-PARITY` | numerical: python vs warp-cpu vs warp-cuda statistical consistency |
| `S-NUM-UNCERTAINTY` | numerical: batch uncertainty validity and coverage |
| `S-CAP-CONTRACT` | capability: fail-closed dispatch, discoverable support, requested/effective config |
| `S-CAP-PLANNING` | capability: beamlets, influence matrices, scoring grids, biological lookups |
| `S-PROV-DATA` | provenance: dataset acquisition, versioning, checksums, caching, offline reuse |
| `S-PROV-OUTPUT` | provenance: self-describing persisted results reopened in a clean process |
| `S-PERF-WORKLOADS` | performance: frozen workload targets and resource metrics |
| `S-WORKFLOW-USER` | workflow: clean-install documented end-to-end user workflow |
| `S-WORKFLOW-QUALITY` | workflow: automated tests, CI, pre-commit, docs build |

## Requirements

| ID | Title | Planned suite(s) | Intended evidence | Status |
|---|---|---|---|---|
| V1-MUST-001 | Ion-therapy treatment planning orientation | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-002 | Python-first implementation | S-WORKFLOW-QUALITY | self-consistency | planned |
| V1-MUST-003 | Reference Python execution | S-NUM-PARITY | backend-parity | planned |
| V1-MUST-004 | NVIDIA Warp acceleration | S-NUM-PARITY | backend-parity | planned |
| V1-MUST-005 | Ion transport | S-PHYS-PROTON-EM, S-PHYS-HELIUM, S-PHYS-CARBON, S-PHYS-OXYGEN | ion-specific-tabulated, independent-monte-carlo, measured | planned |
| V1-MUST-006 | Homogeneous reference geometries | S-PHYS-PROTON-EM | ion-specific-tabulated | planned |
| V1-MUST-007 | Voxelized geometries | S-CAP-PLANNING, S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-008 | Material assignment | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-009 | Beamlet-resolved scoring | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-010 | Dose influence matrices | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-011 | Changeable scoring grids | S-CAP-PLANNING, S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-012 | Extensible scoring | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-013 | Absorbed dose | S-PHYS-PROTON-EM | independent-monte-carlo | planned |
| V1-MUST-014 | Energy deposition | S-PHYS-PROTON-EM | independent-monte-carlo | planned |
| V1-MUST-015 | LET scoring | S-PHYS-LET | independent-monte-carlo | planned |
| V1-MUST-016 | External biological lookup data | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-017 | Analytical physics layer | S-PHYS-PROTON-EM | independent-theory | planned |
| V1-MUST-018 | External/tabulated data layer | S-PROV-DATA | ion-specific-tabulated | planned |
| V1-MUST-019 | Separation of transport from data source | S-CAP-CONTRACT | self-consistency | planned |
| V1-MUST-020 | No large physics datasets in Git | S-PROV-DATA | self-consistency | planned |
| V1-MUST-021 | Automatic acquisition | S-PROV-DATA | self-consistency | planned |
| V1-MUST-022 | Local caching | S-PROV-DATA | self-consistency | planned |
| V1-MUST-023 | Dataset versioning | S-PROV-DATA | self-consistency | planned |
| V1-MUST-024 | Integrity verification | S-PROV-DATA | self-consistency | planned |
| V1-MUST-025 | Provenance | S-PROV-OUTPUT | self-consistency | planned |
| V1-MUST-026 | Offline reuse | S-PROV-DATA, S-WORKFLOW-USER | self-consistency | planned |
| V1-MUST-027 | Explicit random control | S-NUM-UNCERTAINTY | self-consistency | planned |
| V1-MUST-028 | Parallel-safe random streams | S-NUM-UNCERTAINTY, S-NUM-PARITY | self-consistency | planned |
| V1-MUST-029 | Cross-backend statistical consistency | S-NUM-PARITY | backend-parity | planned |
| V1-MUST-030 | Statistical uncertainty estimation | S-NUM-UNCERTAINTY | self-consistency | planned |
| V1-MUST-031 | Validated physical behavior | all S-PHYS-* | measured, ion-specific-tabulated, independent-monte-carlo | planned |
| V1-MUST-032 | Conservation checks | S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-033 | Numerical convergence | S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-034 | High-throughput execution | S-PERF-WORKLOADS | self-consistency | planned |
| V1-MUST-035 | GPU-oriented architecture | S-PERF-WORKLOADS | self-consistency | planned |
| V1-MUST-036 | Reproducible benchmarking | S-PERF-WORKLOADS | self-consistency | planned |
| V1-MUST-037 | Explicit physical units | S-PROV-OUTPUT | self-consistency | planned |
| V1-MUST-038 | Coordinate-system documentation | S-PROV-OUTPUT | self-consistency | planned |
| V1-MUST-039 | Independent validation | all S-PHYS-* | measured, ion-specific-tabulated, independent-monte-carlo | planned |
| V1-MUST-040 | Backend consistency | S-NUM-PARITY | backend-parity | planned |
| V1-MUST-041 | Automated testing | S-WORKFLOW-QUALITY | self-consistency | planned |
| V1-MUST-042 | CI integration | S-WORKFLOW-QUALITY | self-consistency | planned |
| V1-MUST-043 | Pre-commit framework | S-WORKFLOW-QUALITY | self-consistency | planned |
| V1-SCOPE | Complete MUST scope-control priorities | all suites | all categories | planned |
| V2-ION | independently qualified therapeutic ions | S-PHYS-HELIUM, S-PHYS-CARBON, S-PHYS-OXYGEN, S-PHYS-LET | ion-specific-tabulated, independent-monte-carlo, measured | planned |
| V2-CAP | fail-closed capability contract | S-CAP-CONTRACT | self-consistency | planned |
| V2-NUM | adversarial numerical falsification | S-NUM-FALSIFICATION, S-NUM-UNCERTAINTY | self-consistency | planned |
| V2-EVID | physical evidence lineage and sufficiency | all S-PHYS-* | measured, ion-specific-tabulated, independent-monte-carlo | planned |
| V2-BIO | external biological lookup functionality | S-CAP-PLANNING | self-consistency | planned |
| V2-OUTPUT | self-describing persisted scientific results | S-PROV-OUTPUT | self-consistency | planned |
| V2-USER | complete clean-install scientific workflow | S-WORKFLOW-USER | self-consistency | planned |
| V2-PERF | frozen application-level performance/resource gates | S-PERF-WORKLOADS | self-consistency | planned |
| V2-RELEASE | complete qualification and promotion | all suites | all categories | planned |

Status vocabulary: planned, implemented (code merged, evidence pending),
evidenced (exact-SHA evidence recorded), satisfied (full qualification pass).
