# V3-003D acceptance criteria (frozen before any result)

Status: frozen with the first commit of task V3-003D (decision 0039, design entry of 2026-10-07),
before any V3-003D test, local-validation (LV) or host-runner (HR) result exists. Later changes
require a new dated entry in the amendments section with the scientific reason.

Scope: the CSDA range table of the proton transport engine is replaced by the exact integral of the
log-log interpolated stopping power, with the closed-form range between the nodes
(`range_construction = exact-loglog-quadrature-v1`; decision 0039, V3-003D design entry). The
range-carrying step is rejected. The V3-003 and V3-004 qualifications are repeated at the changed
code under unchanged criteria (amendments 24 of `v3-003-acceptance.md` and 6 of
`v3-004-acceptance.md`).

Tiers: CI = GitHub (small, Warp CPU, single process); LV = local exact-SHA validation on the CPU;
HR = controlled host runner with CUDA. Evidence classes follow `experiment/v3/PROTOCOL.md`.

## General rules

- Every D row is graded on the analytic Bethe water source unless stated otherwise.
- Negative controls use a test-only monkeypatch of
  `ionmc.physics.stopping.exact_loglog_range_increments` to the trapezoid formula (the tables are
  data, so it applies to both backends).
- `E_cut` = 2 MeV.
- The numbers marked "predicted" are informative expectations from a pure-Python mock of the
  deterministic step loop on the real Bethe water tables of 1e1f73b2 (`f_short` = 0.01,
  `max_energy_loss_fraction` = 0.02, range step function 0.2 / 0.1 mm, planes every `s_max`). They
  are not pass criteria. The mock reproduced the recorded A4b failures of the trapezoid table
  (5.0e-4, 1.01e-3, 1.29e-3 against 5.3e-4, 1.04e-3, 1.32e-3). If an observation differs from a
  prediction it is reported as such and the pass criterion alone decides.

## Acceptance rows

| # | Check | Pass iff | Class | Tier |
|---|---|---|---|---|
| D1 | Table against an independent quadrature: water and Al (Bethe) and NIST water (when cached; else skipped in CI and run in LV) | (a) Each nodal increment `r_mass[i+1] - r_mass[i]` equals the 32-point Gauss-Legendre quadrature in `u = ln E` of `a E / S_interp(E)` (S from `stopping_mass`, an independent code path) within 1e-12 relative. (b) `range_g_cm2` at 10^4 log-uniform random energies equals `R_0` plus the Gauss-Legendre quadrature from `E_min` within 1e-12 relative. (c) `StoppingTable.csda_range_g_cm2` equals the `TransportTables` nodes within 1e-12 (same grid). (d) Control: the trapezoid increments deviate by at least 1e-5 relative in at least 90 % of the intervals | theory | CI |
| D1b | Round trip and monotonicity (U2 strengthened; the strengthening is informative, U2 at 1e-5 stays graded) | `\|Rinv(R(E))/E - 1\| <= 1e-6` and `\|R(Rinv(r))/r - 1\| <= 2e-6` over the table; `R` and `Rinv` strictly monotone | self-consistency | CI |
| D2 | Branch consistency, engine, deterministic 150 MeV (straggling and MCS off), slab geometry with planes every `s_max` as in A4b, trace of one history; python and warp-cpu float64; `s_max` in {1, 0.5, 0.25, 0.1} mm | For every step boundary `(z_b, E_b)` with `E_b > E_cut`: `\|R(E_b) - (R(E0) - rho z_b / 10)\|` <= 0.5 um (residual-range distance; relative energy agreement of 1e-5 at every depth is not attainable for any remedy, because `dE/E = dr/(p r)` with `p` about 1.7). Control (trapezoid): the maximum is >= 1.0 um at each `s_max`. Predicted: <= 0.23 um (new), 1.6 / 3.2 / 4.1 / 4.7 um (control) at 1 / 0.5 / 0.25 / 0.1 mm | falsification (numerical) | CI (python at 0.1 mm if <= 60 s, else LV) |
| D2b | Step-size independence of the deterministic end depth (the T2 quantity, projected to `E_cut` by the table) | The spread over the four `s_max` is <= 0.5 um. Control (trapezoid): spread >= 2 um. Predicted: 0.1 um (HEAD 3.1 um) | falsification | CI |
| D3 | V3-004 A4b with the table CSDA energy `Rinv(R(E0) - rho z)` as reference (formerly non-gating, footnote 1 of the V3-004 plan) | Worst relative error <= 1e-4 at `s_max` = 1, 0.5 and 0.25 mm. The transported-energy A4b and A4 stay <= 1e-4 (both references now gating). Control: the trapezoid table gives > 1e-4. Predicted: 7.1e-5 / 5.5e-5 / 3.3e-5 (new), 5e-4 / 1.0e-3 / 1.3e-3 (control). If the 1 mm margin (1.4x) fails, the bound is not widened: the observation is recorded as contrary evidence and a follow-up (short-branch midpoint energy from the table) is opened | falsification | CI |
| D4 | Range shift against the previous head (deterministic table numbers) | `R_new(150 MeV)` = 158.6248 +- 0.0002 mm (shift -5.1 um); `\|dR(100 MeV)\|` and `\|dR(200 MeV)\|` <= 10 um; the deterministic transported end-depth shift at 150 MeV is within +-0.5 um of -3.41 / -1.74 / -0.88 / -0.35 um for `s_max` 1 / 0.5 / 0.25 / 0.1 mm (the old values come from running the 1e1f73b2 formula in the test with the trapezoid monkeypatch). Statistical R80 shifts (T5, T9) are reported, informative, `\|dR80\|` <= 10 um expected | regression (intended change) | CI / LV (informative part) |
| D5 | Re-qualification of V3-003 under the frozen criteria (T1-LV 256 x 150 MeV, T2, T4, T5 to T10, T12, T13 chunks, T14, U rows, R1) at the frozen head | The `lv` and `hr` suites pass at the new base 20391004 with every V3-003 criterion and amendments 1 to 23 unchanged; T1 CI passes with the regenerated fixture (amendment 24) | per row | CI / LV / HR |
| D6 | Re-qualification of V3-004 (A1 to A16; A16 per amendment 6; A4b with both references gating) | The `lv4` and `hr4` suites pass at the new base 20401004 with every criterion unchanged | per row | CI / LV / HR |
| D7 | Execution mode | All LV/HR runs execute in single-process diagnostic mode (`run_suite.py --workers 1`, `SINGLE_PROCESS_ENV`). `t13-workers` and `a15-workers` are archived as DEFERRED, not evidence; archives are non-conformant; the ledger entry is "implemented; diagnostic evidence; conformant qualification pending" | process | - |

Predicted shifts of the physics (informative expectations, not pass criteria):

| quantity | expected change |
|---|---|
| Table `R(150 MeV)` | 158.62994 to 158.62484 mm (-5.10 um); the exact-node log-log value is 158.62460 mm |
| `R(100 MeV)`, `R(200 MeV)` | about -2.7 um, about -8.2 um (nodal; the closed form adds at most +0.5 um) |
| Deterministic transported end depth, 150 MeV | -3.41 / -1.74 / -0.88 / -0.35 um for `s_max` 1 / 0.5 / 0.25 / 0.1 mm (100 MeV at 1 mm: -2.68 um; 200 MeV: -3.41 um) |
| End-depth spread over `s_max` | 3.1 um to at most 0.1 um |
| Transported `E(z)` near the end (z = 158.02 mm, 1 mm steps) | 6.6212 to about 6.6000 MeV (-0.3 %) |
| Round trip `eta' = R(Rinv(r))/r - 1` | at most 7.2e-7 (r < 0.01 g/cm2), at most 1e-7 above 1 g/cm2 |
| `compare_nist.py` range deviations (V3-002, informative re-run) | change by at most 3.6e-5 relative |
| Source hashes (`content_sha256`, `source_sha256`) | unchanged (they hash the inputs); `TransportTables.sha256` changes |

## Seeds

- **Previous bases (consumed).** The V3-003 bases 20261004 to 20341004 and the V3-004 bases 20351004
  (rehearsal), 20361004, 20371004 and 20381004 (consumed by observation at the previous head; their
  archives belong to the previous head and are not evidence for this task).
- **V3-003 suites (`lv`, `hr`): base 20391004.** Qualification base of this task.
- **V3-004 suites (`lv4`, `hr4`): base 20401004.** Distinct from the V3-003 base so that the samples
  of the two plans are never correlated.
- **Rehearsals** of V3-003D use the 2041xxxx family only (for example 20411004). They are never
  qualification evidence.
- **Consumption rule (one list for both plans).** A base whose full-scale statistical results have
  been observed (including by an investigation, a rehearsal at that base or a failed qualification)
  is consumed and never reused. The next unused value of 20391004, 20401004, 20411004, ... (skipping
  any rehearsal value used) is frozen by a dated amendment before any result at the changed code is
  observed. CI tests use fixed small seeds, documented in the tests, and are not qualification
  evidence.
- `validation/scripts/transport/a16_digest.py` uses the fixed digest seed 20351004; it is not a
  statistical base and is unchanged.
- Tolerances and history counts do not change when a base is re-frozen.

Runner changes (implementation phase C2):

- `run_suite.py`: `QUALIFICATION_SEED_BASE = 20391004`, `V4_QUALIFICATION_SEED_BASE = 20401004`; the
  consumed tuples gain 20341004 (lv/hr, observed at the old head) and 20381004 (lv4/hr4); docstring
  and help text are updated.
- `steps.py`: `DEFAULT_SEED_BASE = 20391004`; `steps_v4.py`: its default as well.
- The hashed `SOURCE_FILES` of every suite gain `validation/plans/v3-003d-acceptance.md`.
- `steps_v4.py a16` gets `--mode {regression,intended-change}`; the suite default is `regression`.
  The V3-003D exception is the explicit record `A16_INTENDED_CHANGE` in `run_suite.py` (task
  V3-003D, baseline a524f209, identity field `range_construction`, baseline value absent, new value
  `exact-loglog-quadrature-v1`); while it is set, `lv4` passes it to `a16 --mode intended-change`,
  which verifies at run time that the baseline tree's table identity is the baseline value and the
  tree under test's the new value (fail closed otherwise; after the merge the baseline already
  carries the new value, so the next task's first commit deletes the record and advances
  `A16_BASELINE`, amendment 6). The verified identities are written to the step output.
- The `lv` pytest step `pytest-warp-cpu-t2-t4-t11-c1-t13` gains `tests/ionmc/test_range_quadrature.py`
  with `IONMC_REQUIRE_NIST=1` and `IONMC_CACHE_DIR=<workspace>/.ionmc-cache`: the NIST-water case
  of D1 cannot skip in LV (a missing cache fails; only the missing-cache condition is tolerated
  elsewhere, every other loader error fails). D3 also asserts the A4 table-CSDA error (3.4e-5).

## Single-process diagnostic mode

Operator directive 2026-10-07 (amendment 23 of the V3-003 plan, amendment 4 context of the V3-004
plan): until the host instability is understood, every LV/HR run of this task executes in single-process diagnostic
mode with `--workers 1` and unchanged histories, seeds and criteria. The multiprocessing-specific
checks (`t13-workers` of the V3-003 suites, `a15-workers` of the V3-004 suites) are archived as
deferred and are not evidence. All V3-003D qualification archives are therefore non-conformant
diagnostic evidence; the ledger records "implemented; diagnostic evidence; conformant qualification
pending".

## Evidence classes and scope

- D1, D4: theory and regression of an intended change (deterministic bounds, never z-scores).
- D1b: self-consistency. D2, D2b, D3: numerical falsification, with a trapezoid negative control.
- D5, D6: re-qualification of the existing plans, each row under its own class.
- D7: process. No row is independent-Monte-Carlo evidence. Nothing here is a clinical claim.

## Evidence archives (which archive satisfies which row)

To be filled from committed result files (qualification bases 20391004 for `lv`/`hr`, 20401004 for
`lv4`/`hr4`).

| Rows | Archive | Head SHA | Run id | Status |
|---|---|---|---|---|
| D1, D1b, D2, D2b, D3, D4 (CI part), T1 CI, U rows | CI (GitHub) and the pytest steps of `lv` (incl. D1 NIST water, required) and `lv4` | pending | pending | pending |
| D4 informative R80 shifts, D5 (LV rows), T1-LV | `lv` | pending | pending | pending |
| D5 (HR rows: T12, CUDA) | `hr` | pending | pending | pending |
| D6 (LV rows), D3 A4b | `lv4` | pending | pending | pending |
| D6 (HR rows) | `hr4` | pending | pending | pending |
| D7 | all archives, ledger entry | pending | pending | pending |

## Amendments

(None yet.)

## Documentation sweep list for phase C2

Places that state the trapezoid construction, or hard-code range values that a shift of about 5 um
could break, found by `grep -rn` for "trapezoid", "158.6", "range_g_cm2" and "a524f209" on
1e1f73b2 (the phase C1 commit changes none of them).

States the trapezoid construction (to be rewritten to the exact quadrature):

- `src/ionmc/physics/stopping.py`: `build_table` docstring and body (line 346 ff., `cumsum(0.5 (f1 + f0) dln)` at line 359), `StoppingTable` docstring (line 252, "interpolated in ln E, ln S and ln R"; add the construction and the between-node error), `BetheStoppingSource` and `NistStarStoppingSource` docstrings.
- `src/ionmc/transport/tables.py`: module docstring lines 5 to 10 ("trapezoid rule in ln E"), `from_stopping_tables` (lines 180 to 210, the re-integration from `t.csda_range_g_cm2[0]`).
- `docs/physics/stopping-power.md` line 46 ff. ("The CSDA range is the trapezoid integral").
- `docs/physics/em-transport.md` lines 12 to 36 (`dR = E d ln E / S`, the `ln R` table, the round-trip statement of line 37).
- `docs/architecture/transport.md` lines 166 to 178 (the finding paragraph and the inverse-range table density) and the A16 row at line 645 ("digests against `git archive a524f209`").
- `validation/scripts/stopping/compare_nist.py` line 121 (comment "the same trapezoid integral as for the analytic tables").
- `tests/ionmc/test_scoring_transport.py` lines 207 and 321 (docstring and comment about the trapezoid bias of the table CSDA comparison; the table comparison becomes gating), 852 (a test-local range integral "trapezoid in ln E"); line 291 uses `np.trapezoid` for an unrelated energy average (keep).
- `tests/ionmc/test_transport_diagnostic.py` lines 65 to 86 (T1 fixture provenance and the 8e-11 midpoint allowance; amendment 24).
- `decisions/0040-let-and-extensible-scoring.md` follow-up note (link to the 0039 design entry; done in C1).

Uses of range values to check against a 5 um shift (assertion tolerances, not code changes unless a
check fails):

- `tests/ionmc/test_tables_scattering.py` lines 76 to 99: `range_g_cm2` against `source_table.range_at` with rel 1e-4 (5 um of 158 mm is 3e-5: passes, check the 1e-13 round-trip at line 99 against the new inverse), the finite-difference check at line 89.
- `tests/ionmc/test_transport_reference.py` lines 43, 44, 71, 212, 213 and `tests/ionmc/test_transport_warp.py` lines 122, 123, 151: depths computed from the table itself (self-consistent; check expected-energy tolerances).
- `tests/ionmc/test_scoring_transport.py` lines 210, 228, 260, 1053: table-derived reference energies and path bounds (A4, A4b, the capacity-bound test).
- `tests/ionmc/test_stopping.py` lines 85 to 88, 187, 257 to 285: `StoppingTable` validation fixtures built from `csda_range_g_cm2` (the metadata and hash of tables change).
- `tests/ionmc/test_transport_scoring.py` line 34: `depth = 1.3 * 158.6` (geometric margin, harmless).
- `validation/scripts/transport/a16_digest.py` line 61: `depth = 1.1 * 158.6` (geometric, harmless); line 5 and `steps_v4.py` line 111 (`A16_BASELINE = "a524f209"`: kept; amendment 6 of the V3-004 plan; the first task after the V3-003D merge sets it to the merge SHA).
- `validation/scripts/transport/steps.py` lines 75 and 287: `TABLES.range_g_cm2` for the R80/T2 geometry and the cut-off range (table-derived, self-consistent).
- `src/ionmc/transport/parity.py` line 955 and `src/ionmc/config.py` line 380: table-derived range for geometry and the path bound (the path bound `B_L` of decision 0040 changes by about 3e-5 relative; check that the capacity proof margin is unaffected).
- `src/ionmc/transport/mcs_checks.py` lines 79, 243 to 309, `validation/scripts/transport/mcs_checks.py` line 97 and `docs/physics/em-transport.md` line 150, `docs/research/*.md`: 158.6 MeV is the Gottschalk beam energy, not a range; the U4 check (`theta_dM` by quadrature along the CSDA path) uses `range_g_cm2` through `r1_g`; its numbers shift by about 3e-5 relative (check the stated deviations of `test_tables_scattering.py` line 206 ff.).
- `docs/architecture/transport.md` line 645 and `validation/scripts/transport/steps_v4.py` / `a16_digest.py`: every use of `a524f209` (re-grep in C2 and justify each remaining hit).
- Archived summaries under `validation/results/transport/` that record the old table sha or R values belong to the old head and are left unchanged.
