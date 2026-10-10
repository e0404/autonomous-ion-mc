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
| `S-PHYS-PROTON-NUCLEAR` | physics: proton nuclear attenuation and secondaries (decision 0041; frozen rows `validation/plans/v3-005-acceptance.md`: slice A P1-P5, N1, V1, V1b informative, V2, V3, V4, V4b, V9, R1, X1, C1, D6, E1; slice B V2b, V5-V8; amended 2026-10-07: Amendments 1-3 of the plan, D6 both tiers failed and local alpha deposition kept under a declared ceiling; slice B evaluated 2026-10-09 at 7aae5bb (V3-005B): V8-LV, R1-nuc, V2b, V7 scan/shift/f32 pass, V5 and V7 replicate coverage FAIL, V8 python:cpu64 nuc_local inconclusive; V6, D2, E1-B moved to V3-005C; V3-005C C2: P6 FAIL for the H-1 evaluation below 12.53 MeV (frozen item (4), 542 negative nodes/midpoints, recorded; P6-D revised-domain row pre-registered, Amendment 15); elastic table 25361f81… (C2c; C2b build f932d6ae…, identical npz) built, P6 items (1)-(3),(5)-(7) pass; P6-D and E-shape listed in the Amendment 15 row-table addendum) safeguard (d)1 MF6 rule: fail as written (zap3007 14/30, zap5012 10/30); superseded by Amendment 16 spectra evidence: pass (KS D 0.042-0.230 at 9 of 9, all > 0.02; table 1fd24cff…); |
| `S-PHYS-LET` | physics: track- and dose-averaged LET including mixed fields |
| `S-PHYS-HELIUM` | physics: helium range, primary attenuation, fragment (Z=1,2) yields and distal dose |
| `S-PHYS-CARBON` | physics: carbon range, attenuation, fragment build-up, distal dose |
| `S-PHYS-OXYGEN` | physics: oxygen range, primary attenuation, charge-resolved fragment build-up, distal fragment dose and mixed-field LET within the documented supported domain (16O in water/tissue, 100-430 MeV/u entrance energy); any oxygen capability outside this domain must fail closed and is not claimed |
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
| V1-MUST-005 | Ion transport | S-PHYS-PROTON-EM, S-PHYS-PROTON-NUCLEAR, S-PHYS-HELIUM, S-PHYS-CARBON, S-PHYS-OXYGEN | ion-specific-tabulated, independent-monte-carlo, measured | in progress (proton EM implemented V3-003A/B; proton non-elastic nuclear V3-005A python reference backend with diagnostic-mode evidence @b84fdf3 (VAL-20261008-064717-C8A3EE; every lv5 row passes; non-conformant by code: single-process directive, deferred partition row, imported partials), Warp CPU/CUDA nuclear kernels and secondary-proton interactions V3-005B @7aae5bb (`validation/results/transport/lv5b-7aae5bb-single-process.json`, `validation/results/transport/hr5-7aae5bb-single-process.json`; non-conformant by code: single-process directive, imported partials, deferred partition row; validation records `failed` at the evidence commit following 7aae5bb1): V8-LV trajectory parity pass, R1-nuc pass, V2b pass; V5 vs TOPAS FAIL on peak/plateau (+12.5 % / +12.3 %) and the nuclear on-off plateau integral (−18.5 % / −14.7 %), pass on plateau, R80 and total deposit, hypothesis: the absent hadronic elastic channel (to be tested by V11/V5 in V3-005C; the exploratory X cases show elastic is not the only contributor); evidence grade per `validation/plans/v3-005-acceptance.md`; V3-005C C2: P6 FAIL for the H-1 evaluation below 12.53 MeV (frozen item (4), 542 negative nodes/midpoints, recorded; P6-D revised-domain row pre-registered, Amendment 15); elastic table 25361f81… (C2c; C2b build f932d6ae…, identical npz) built, P6 items (1)-(3),(5)-(7) pass (data layer only, no elastic transport result); He/C/O planned) safeguard (d)1 MF6 rule: fail as written (zap3007 14/30, zap5012 10/30); superseded by Amendment 16 spectra evidence: pass (KS D 0.042-0.230 at 9 of 9, all > 0.02; table 1fd24cff…); |
| V1-MUST-006 | Homogeneous reference geometries | S-PHYS-PROTON-EM | ion-specific-tabulated | planned |
| V1-MUST-007 | Voxelized geometries | S-CAP-PLANNING, S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-008 | Material assignment | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-009 | Beamlet-resolved scoring | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-010 | Dose influence matrices | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-011 | Changeable scoring grids | S-CAP-PLANNING, S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-012 | Extensible scoring | S-CAP-PLANNING | self-consistency | implemented (V3-004); diagnostic-mode evidence @fccb687 (base 20381004) and @637ef82 after the V3-003D range-table change (lv4/hr4 archives at base 20401004; non-conformant: single-process directive 2026-10-07, a15-workers deferred); conformant exact-SHA qualification pending |
| V1-MUST-013 | Absorbed dose | S-PHYS-PROTON-EM, S-PHYS-PROTON-NUCLEAR | independent-monte-carlo, measured | implemented (V3-004); diagnostic-mode self-consistency evidence @fccb687 and @637ef82 (non-conformant, see V1-MUST-012); nuclear contribution implemented on the python backend (V3-005A: nuclear_local and secondary-p/d dose, rows V3/X1/E1 pass, diagnostic-mode evidence @b84fdf3 (VAL-20261008-064717-C8A3EE; every lv5 row passes; non-conformant by code: single-process directive, deferred partition row, imported partials); gating independent MC V5 evaluated in V3-005B @7aae5bb (`validation/results/transport/lv5b-7aae5bb-single-process.json`): FAIL on peak/plateau and nuclear on-off plateau integral, pass on total deposit, plateau and R80 (hypothesis: the absent hadronic elastic channel (to be tested by V11/V5 in V3-005C; the exploratory X cases show elastic is not the only contributor)); measured V6 moved to V3-005C); conformant qualification and independent Monte Carlo (T15) pending |
| V1-MUST-014 | Energy deposition | S-PHYS-PROTON-EM, S-PHYS-PROTON-NUCLEAR | independent-monte-carlo, self-consistency | implemented (V3-004); diagnostic-mode self-consistency evidence @fccb687 and @637ef82 (non-conformant, see V1-MUST-012); nuclear energy bookkeeping implemented (V3-005A: nuclear_local / alpha_local / escaped neutron and gamma / binding / signed imbalance tallies, balance row V3 at 1e-12 with independent AME recomputation of the binding term, diagnostic-mode evidence @b84fdf3 (VAL-20261008-064717-C8A3EE; every lv5 row passes; non-conformant by code: single-process directive, deferred partition row, imported partials); warp backends V3-005B @7aae5bb: V8-LV trajectory parity incl. the nuclear ledger, V5 total deposit within 1 % of TOPAS (−0.30 % / −0.54 %), `validation/results/transport/lv5b-7aae5bb-single-process.json`); conformant qualification and independent Monte Carlo (T15) pending |
| V1-MUST-015 | LET scoring | S-PHYS-LET | independent-monte-carlo | implemented (V3-004); diagnostic-mode theory/self-consistency evidence @fccb687 (A1-A9, A11, A14-A16) and @637ef82 (non-conformant, see V1-MUST-012); conformant qualification and independent Monte Carlo (TOPAS ProtonLET, V3-010) pending |
| V1-MUST-016 | External biological lookup data | S-CAP-PLANNING | self-consistency | implemented (V3-004); diagnostic-mode evidence @fccb687 (A10, A12; synthetic fixture only, no clinical tables; non-conformant, see V1-MUST-012); conformant qualification pending |
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
| V1-MUST-029 | Cross-backend statistical consistency | S-NUM-PARITY | backend-parity | planned (nuclear-on V3-005B @7aae5bb, `validation/results/transport/hr5-7aae5bb-single-process.json`, diagnostic mode, non-conformant: V8 statistical parity passes for every CPU/CUDA pair; python:cpu64 idd and sec_p pass, nuc_local inconclusive by the T12 sparse-profile rule at 2.4e4 python histories) |
| V1-MUST-030 | Statistical uncertainty estimation | S-NUM-UNCERTAINTY | self-consistency | planned (nuclear secondary estimators V3-005B @7aae5bb, `validation/results/transport/lv5b-7aae5bb-single-process.json`: V7 relative-SE slopes pass; V7 replicate coverage FAILS: escaped neutral energy fails the paired gate (lower bound 0.6445 < 0.6464), sec_p and nuclear_local pass both gates; re-test pre-registered for V3-005C) |
| V1-MUST-031 | Validated physical behavior | all S-PHYS-* | measured, ion-specific-tabulated, independent-monte-carlo | planned |
| V1-MUST-032 | Conservation checks | S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-033 | Numerical convergence | S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-034 | High-throughput execution | S-PERF-WORKLOADS | self-consistency | planned |
| V1-MUST-035 | GPU-oriented architecture | S-PERF-WORKLOADS | self-consistency | planned |
| V1-MUST-036 | Reproducible benchmarking | S-PERF-WORKLOADS | self-consistency | planned |
| V1-MUST-037 | Explicit physical units | S-PROV-OUTPUT | self-consistency | planned |
| V1-MUST-038 | Coordinate-system documentation | S-PROV-OUTPUT | self-consistency | planned |
| V1-MUST-039 | Independent validation | all S-PHYS-* | measured, ion-specific-tabulated, independent-monte-carlo | planned |
| V1-MUST-040 | Backend consistency | S-NUM-PARITY | backend-parity | planned (nuclear-on V3-005B @7aae5bb: V8-LV python vs warp-cpu f64 trajectory parity, identical events and genealogy, continuous difference 0.0 within 1e-10, `validation/results/transport/lv5b-7aae5bb-single-process.json`; V7 f32 vs f64 on CUDA pass, `validation/results/transport/hr5-7aae5bb-single-process.json`) |
| V1-MUST-041 | Automated testing | S-WORKFLOW-QUALITY | self-consistency | planned |
| V1-MUST-042 | CI integration | S-WORKFLOW-QUALITY | self-consistency | planned |
| V1-MUST-043 | Pre-commit framework | S-WORKFLOW-QUALITY | self-consistency | planned |
| V1-SCOPE | Complete MUST scope-control priorities | all suites | all categories | planned |
| V2-ION | independently qualified therapeutic ions | S-PHYS-HELIUM, S-PHYS-CARBON, S-PHYS-OXYGEN, S-PHYS-LET | ion-specific-tabulated, independent-monte-carlo, measured | planned |
| V2-CAP | fail-closed capability contract | S-CAP-CONTRACT, S-PHYS-PROTON-NUCLEAR | self-consistency | implemented for nuclear runs (V3-005A row C1, diagnostic-mode evidence @b84fdf3 (VAL-20261008-064717-C8A3EE; every lv5 row passes; non-conformant by code: single-process directive, deferred partition row, imported partials): unsupported element, E > 250 MeV, non-proton source, nist-star, warp backend with nuclear, missing/stale/mis-pinned table, unproducible species, every nuclear overflow raises or invalidates; V3-005B: warp backends support nuclear runs, the warp case of C1 is replaced by the forced-counter cases on warp-cpu (Amendment 7(g)), secondaries above the 250 MeV table domain are transported without nuclear sampling (Amendment 10(a))) |
| V2-NUM | adversarial numerical falsification | S-NUM-FALSIFICATION, S-NUM-UNCERTAINTY | self-consistency | implemented (V3-003B T9/T14 negative controls, V3-004 A7/A9, V3-003D D1-D4 with trapezoid negative controls); diagnostic-mode evidence @637ef82 (non-conformant, see V1-MUST-012); V3-005A thinning step-independence probes P4 and V2-probe; V3-005B @7aae5bb (`validation/results/transport/lv5b-7aae5bb-single-process.json`): V2b (s_max 0.1 vs 1.0 mm) pass, V7 scan and grid shift/refinement pass, V7 replicate coverage FAIL (escaped neutral energy, paired gate) — uncertainty-coverage evidence of row V7 not established; conformant qualification pending |
| V2-EVID | physical evidence lineage and sufficiency | all S-PHYS-* | measured, ion-specific-tabulated, independent-monte-carlo | in progress (decision 0041 data-role table: LA150 and AME2020 construction; EXFOR from 1997 on evaluation, informative V1b with the pre-declared LA150 p+C deficiency; pre-1997 EXFOR and geant-val report-only; MCsquare shared sigma lineage; TOPAS independent MC; V3-005C C2b: PDG rpp2022 p-p compilation role corrected from construction to evaluation (report-only V10, use id 770debc97dbd4d00b4cb2e6e44a6a8b0), S(E) above 150 MeV constructed from the Geant4 BGG p-p formula only) |
| V2-BIO | external biological lookup functionality | S-CAP-PLANNING | self-consistency | implemented (V3-004); diagnostic-mode evidence @fccb687 (lookup tallies with provenance; synthetic fixture only; non-conformant, see V1-MUST-012); conformant qualification pending |
| V2-OUTPUT | self-describing persisted scientific results | S-PROV-OUTPUT | self-consistency | planned |
| V2-USER | complete clean-install scientific workflow | S-WORKFLOW-USER | self-consistency | planned |
| V2-PERF | frozen application-level performance/resource gates | S-PERF-WORKLOADS | self-consistency | planned |
| V2-RELEASE | complete qualification and promotion | all suites | all categories | planned |

Status vocabulary: planned, in progress (task open, nothing merged as evidence), implemented (code merged, evidence pending),
evidenced (exact-SHA evidence recorded), satisfied (full qualification pass).
