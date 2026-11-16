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
| `S-PHYS-PROTON-NUCLEAR` | physics: proton nuclear attenuation and secondaries (decision 0041; frozen rows `validation/plans/v3-005-acceptance.md`: slice A P1-P5, N1, V1, V1b informative, V2, V3, V4, V4b, V9, R1, X1, C1, D6, E1; slice B V2b, V5-V8; amended 2026-10-07: Amendments 1-3 of the plan, D6 both tiers failed and local alpha deposition kept under a declared ceiling) |
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
| V1-MUST-005 | Ion transport | S-PHYS-PROTON-EM, S-PHYS-PROTON-NUCLEAR, S-PHYS-HELIUM, S-PHYS-CARBON, S-PHYS-OXYGEN | ion-specific-tabulated, independent-monte-carlo, measured | in progress (proton EM implemented V3-003A/B; proton non-elastic nuclear V3-005A python reference backend with diagnostic-mode evidence @b84fdf3 (VAL-20261008-064717-C8A3EE; every lv5 row passes; non-conformant by code: single-process directive, deferred partition row, imported partials), Warp kernels V3-005B; evidence grade per `validation/plans/v3-005-acceptance.md`; **He/C/O electromagnetic transport implemented in V3-008** (decision 0043; PR #84 open, NOT integrated: row I3 FAIL, plan F.6 unmet). Evidence (runs at f5f3a701: section-A rows I6-LV, I8, I13, T1-ion LV, R1-ion, R2, T12-ion LV, PASS; 712f484a: I9, I12, I14, I15, T13-ion, T12-ion HR, PASS; TOPAS main set run at ead583f1 and graded at 94bc948e; held-out set run at 55e6a180, records b18fa715, graded at 1fc287cf; bef32873 / 769f2dd6: V5-SS PASS 5/5 on block 2 graded at 833ccdcd, shared-lineage single scattering; grades bound in the plan evidence archive). **Rows I2 and I3 run once at 7d0e796b (2026-11-15): I2 PASS (He vs ASTAR / ICRU 90, regression guard); I3 FAIL (C/O vs libdEdx ICRU73_NEW (revised ICRU 73 lineage, identity with the 2009 revision OPEN): stopping gates fail: max |dev| 10.1 / 13.5 % (C / O) at 2–10 MeV/u against 6 %, 2.51 / 5.05 % at ≥ 10 MeV/u against 1.5 %; the CSDA-range gate passes (0.26 / 0.39 %); the attribution of the deficit to energies below ≈ 50 MeV/u is a hypothesis WITHDRAWN by the pre-declared rule of amendment I3-FU-1 step 1 (committed band diagnostic at 5999f789: nine of ten quoted band quantities match at 0.1 pp, the oxygen ≥ 50 MeV/u rounded mean does not; the committed diagnostic records the band-wise deficit; z-scan: [10,20) MeV/u z²-like, other bands not attributed to a single term); step 2 (TOPAS thin-slab measurement) graded once at c223f5e5: INCONCLUSIVE on two instrument defects (stderr banner rule, single-bin EDep parser; nothing measured), revision 16 re-run pending; preserved, never re-graded) — a gating FAIL: plan F.6 is not met and V3-008 is NOT integrated at this head; follow-up per plan D (research record, prospective amendment with fresh evidence) (`validation/results/stopping/v3-008-i2/`, `v3-008-i3/`).** TOPAS-agreement main set (2026-11-09): TD-EM 3 PASS / 6 FAIL, TL-EM 1 PASS / 8 FAIL, consumed and preserved under amendment R2-STOP3; held-out TD-H/TL-H PASS 12/12 at He 100 / C 150 / O 200 MeV/u — the TOPAS-agreement claim is limited to these three interior energies, the f_dM energy-trend test D_ext was not evaluable (R2-STOP3 item 3; the operator may override the Branch-B route). Gaps (decision 0043 §3): V3 measured C/O range EVIDENCE GAP; no independent evidence of the charge dependence of lateral spreading and no independent coverage of the f_dM extrapolation (intervention IR-20261102-010435-B59B2D, resolved 2026-11-07 by the operator directive: experimental validation deferred); ion nuclear interactions not implemented → V3-009) |
| V1-MUST-006 | Homogeneous reference geometries | S-PHYS-PROTON-EM | ion-specific-tabulated | planned |
| V1-MUST-007 | Voxelized geometries | S-CAP-PLANNING, S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-008 | Material assignment | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-009 | Beamlet-resolved scoring | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-010 | Dose influence matrices | S-CAP-PLANNING | self-consistency | planned |
| V1-MUST-011 | Changeable scoring grids | S-CAP-PLANNING, S-NUM-FALSIFICATION | self-consistency | planned |
| V1-MUST-012 | Extensible scoring | S-CAP-PLANNING | self-consistency | implemented (V3-004); diagnostic-mode evidence @fccb687 (base 20381004) and @637ef82 after the V3-003D range-table change (lv4/hr4 archives at base 20401004; non-conformant: single-process directive 2026-10-07, a15-workers deferred); conformant exact-SHA qualification pending |
| V1-MUST-013 | Absorbed dose | S-PHYS-PROTON-EM, S-PHYS-PROTON-NUCLEAR | independent-monte-carlo, measured | implemented (V3-004); diagnostic-mode self-consistency evidence @fccb687 and @637ef82 (non-conformant, see V1-MUST-012); nuclear contribution implemented on the python backend (V3-005A: nuclear_local and secondary-p/d dose, rows V3/X1/E1 pass, diagnostic-mode evidence @b84fdf3 (VAL-20261008-064717-C8A3EE; every lv5 row passes; non-conformant by code: single-process directive, deferred partition row, imported partials); gating independent MC V5 and measured V6 in V3-005B); conformant qualification and independent Monte Carlo (T15) pending |
| V1-MUST-014 | Energy deposition | S-PHYS-PROTON-EM, S-PHYS-PROTON-NUCLEAR | independent-monte-carlo, self-consistency | implemented (V3-004); diagnostic-mode self-consistency evidence @fccb687 and @637ef82 (non-conformant, see V1-MUST-012); nuclear energy bookkeeping implemented (V3-005A: nuclear_local / alpha_local / escaped neutron and gamma / binding / signed imbalance tallies, balance row V3 at 1e-12 with independent AME recomputation of the binding term, diagnostic-mode evidence @b84fdf3 (VAL-20261008-064717-C8A3EE; every lv5 row passes; non-conformant by code: single-process directive, deferred partition row, imported partials)); conformant qualification and independent Monte Carlo (T15) pending |
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
| V2-ION | independently qualified therapeutic ions | S-PHYS-HELIUM, S-PHYS-CARBON, S-PHYS-OXYGEN, S-PHYS-LET | ion-specific-tabulated, independent-monte-carlo, measured | in progress (V3-008, PR #84 open, not integrated — row I3 FAIL, plan F.6 unmet: He/C/O electromagnetic transport; ion-specific stopping evidence I13 (R80 vs CSDA, PASS) committed, **I2 (He table reproduction) PASS and I3 (C/O vs libdEdx ICRU73_NEW, revised ICRU 73 lineage per P5b, identity with the 2009 revision OPEN) FAIL at 7d0e796b — committed result: the stopping gates fail with max |dev| 10.1 / 13.5 % at 2–10 MeV/u and 2.5 / 5.1 % at ≥ 10 MeV/u, the range gate passes; the finer attribution (deficit concentrated below ≈ 50 MeV/u) is a hypothesis WITHDRAWN by the pre-declared rule of amendment I3-FU-1 step 1 (committed band diagnostic at 5999f789: nine of ten quoted band quantities match at 0.1 pp, the oxygen ≥ 50 MeV/u rounded mean does not; the committed diagnostic records the band-wise deficit; z-scan: [10,20) MeV/u z²-like, other bands not attributed to a single term); step 2 (TOPAS thin-slab measurement) graded once at c223f5e5: INCONCLUSIVE on two instrument defects (stderr banner rule, single-bin EDep parser; nothing measured), revision 16 re-run pending; every beam traverses that region, so it is a model limitation, not an excludable domain; gating FAIL, V3-008 not integrated at this head (plan F.6), follow-up per plan D**; lateral scattering (exact-Ω mapping; V5-SS single-scattering PASS 5/5, shared lineage) and TOPAS-agreement evidence restricted to the held-out energies He 100 / C 150 / O 200 MeV/u (EM-only intermediate; main set consumed with FAILs under R2-STOP3); documented supported domains per species (decision 0043 §5; the low ends and the 430 MeV/u ends were not re-measured after the R2-STOP3 model change); measured-data rows are EVIDENCE GAPs (experimental validation deferred, operator directive 2026-11-07); carbon/oxygen nuclear interactions, fragmentation, secondary ions and distal fragment dose NOT implemented → V3-009; the release cannot claim ion nuclear capability) |
| V2-CAP | fail-closed capability contract | S-CAP-CONTRACT, S-PHYS-PROTON-NUCLEAR | self-consistency | implemented for nuclear runs (V3-005A row C1, diagnostic-mode evidence @b84fdf3 (VAL-20261008-064717-C8A3EE; every lv5 row passes; non-conformant by code: single-process directive, deferred partition row, imported partials): unsupported element, E > 250 MeV, non-proton source, nist-star, warp backend with nuclear, missing/stale/mis-pinned table, unproducible species, every nuclear overflow raises or invalidates) |
| V2-NUM | adversarial numerical falsification | S-NUM-FALSIFICATION, S-NUM-UNCERTAINTY | self-consistency | implemented (V3-003B T9/T14 negative controls, V3-004 A7/A9, V3-003D D1-D4 with trapezoid negative controls); diagnostic-mode evidence @637ef82 (non-conformant, see V1-MUST-012); V3-005A thinning step-independence probes P4 and V2-probe (s_max 0.1 mm transport probe V2b in V3-005B) in progress; conformant qualification pending |
| V2-EVID | physical evidence lineage and sufficiency | all S-PHYS-* | measured, ion-specific-tabulated, independent-monte-carlo | in progress (decision 0041 data-role table: LA150 and AME2020 construction; EXFOR from 1997 on evaluation, informative V1b with the pre-declared LA150 p+C deficiency; pre-1997 EXFOR and geant-val report-only; MCsquare shared sigma lineage; TOPAS independent MC) |
| V2-BIO | external biological lookup functionality | S-CAP-PLANNING | self-consistency | implemented (V3-004); diagnostic-mode evidence @fccb687 (lookup tallies with provenance; synthetic fixture only; non-conformant, see V1-MUST-012); conformant qualification pending |
| V2-OUTPUT | self-describing persisted scientific results | S-PROV-OUTPUT | self-consistency | planned |
| V2-USER | complete clean-install scientific workflow | S-WORKFLOW-USER | self-consistency | planned |
| V2-PERF | frozen application-level performance/resource gates | S-PERF-WORKLOADS | self-consistency | planned |
| V2-RELEASE | complete qualification and promotion | all suites | all categories | planned |

Status vocabulary: planned, in progress (task open, nothing merged as evidence), implemented (code merged, evidence pending),
evidenced (exact-SHA evidence recorded), satisfied (full qualification pass).
