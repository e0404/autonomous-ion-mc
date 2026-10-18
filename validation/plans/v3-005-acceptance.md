# V3-005 acceptance criteria (frozen before any result)

Status: frozen with the second commit of task V3-005A (the first, 9448d5e, was the A16 housekeeping: record deleted, baseline advanced), before any V3-005 test result, local-validation (LV) or host-runner (HR) result exists. This commit includes decision 0041 (`decisions/0041-proton-nuclear-interactions.md`) and no nuclear code or data result.

- Base: v3/develop at `f3a1dd62` (after V3-003A/B, V3-004, V3-003D, V3-010B).
- Rows of both slices (A and B) are frozen now. Rows marked **slice B** are evaluated by task V3-005B.
- Later changes require a new dated entry in the amendments section, stating the scientific reason.

**Scope.** Proton non-elastic nuclear interactions in water and soft tissue for 1–250 MeV protons (decision 0041), on the electromagnetic engine of decision 0039 and the scoring of decision 0040.
- Slice A covers the data pipeline, the shared nuclear functions, and the nuclear branch of the python reference backend.
- Slice B covers the Warp kernels, the generation queue, p-p elastic scattering (D2), and kernel parity.

**Decision references.** 0037 (purposes, amended), 0038 (data roles), 0039 (nuclear row, amended), 0040 (species, hook, capacity bounds), 0041 (this physics).

**Tiers.**

| Tier | Meaning |
|---|---|
| CI | GitHub CI: synthetic fixtures, no nuclear data. Data-dependent tests skip in CI but fail under `IONMC_REQUIRE_DATA=1`. |
| LV | Local exact-SHA validation on the CPU, with hash-verified tables staged in `.ionmc-cache/`. |
| HR | Controlled host runner (slice B only). |

**Evidence classes** follow `experiment/v3/PROTOCOL.md`. Statistics use `statistics.NormalDist` plus Wilson–Hilferty, without scipy.

## Notation

- **Σ(E)**: the runtime macroscopic non-elastic cross section of a material [cm²/g], composed from the per-element σ of the built table (ENDF/B-VIII.0 MF3/MT5 ≤ 150 MeV, Tripathi-shape extension above).
- **Σ_ENDF(E)**: the same composition evaluated with the native ENDF TAB1 interpolation laws, independently in numpy.
- **TOST form**: pass iff |x − r| + 1.645 σ_comb ≤ tol. σ_comb comes from the batch standard errors (20 batches unless stated).
- **S(z)**: primary survival. It is the primary fluence (V3-004 `fluence`, `generation="primary"`) at depth z divided by its value in the first voxel.
- **ρ_f**: the float rounding term, 1e-12 for python/float64.
- **Counters 0**: every transport, channel and nuclear counter is zero (`valid=True`).

## Criteria

| # | Check | Pass iff | Class | Tier | Slice |
|---|---|---|---|---|---|
| P1 | ENDF-6, EXFOR and AME2020 parsers on hand-written synthetic fixtures: TEXT/CONT/HEAD/LIST/TAB1/TAB2; INT 1–5; NR > 1; floats without `E`; MF6 LAW=1 LANG=1 and LANG=2 with NA=1; LCT=3; EXFOR units MB/B and MEV vs MEV/A; AME `#` estimated values | Exact equality with the hand-derived values (integers exact, floats bitwise) | self-consistency | CI | A |
| P2 | Samplers: Kalbach μ for a ∈ {0.01, 1, 5} × r ∈ {0, 0.5, 1}; Poisson for λ ∈ {0.1, 1, 4}; **A1 (2026-10-07):** floor+Bernoulli (`multiplicity_round`) for λ ∈ {0.1, 1.0, 2.4}, exact two-point frequencies; inverse-CDF E′ on a synthetic 64-bin table; 10⁶ draws each from a fixed CI seed | χ² goodness of fit (≥ 50 equiprobable bins, ≥ 5 expected per bin) with p > 0.001 for every case; the Kalbach first moment within 4 σ | self-consistency | CI | A |
| P3 | Per-event ledger closure (**A1, 2026-10-07**; was: per-event conservation), 10⁵ events from synthetic tables, float64 | ~~max \|ΔE\| ≤ 1e-9 MeV and max \|Δp\| ≤ 1e-9 MeV/c (four-vector of p + target vs products + residual + E\*); ΔZ = ΔA = 0 for every event; no event with m_r < M_r accepted~~ **A1 (2026-10-07):** \|T₁ − Σ T_lab − T_r − nuclear_binding − Δ\| ≤ 1e-9 MeV for every event (Δ = the event's `nuclear_imbalance`); ΔZ = ΔA = 0 for every event; no accepted event with a residual outside AME2020. Per-event energy and momentum conservation are no longer claimed. | conservation | CI | A |
| P4 | Thinning harness (shared funcs via twins, through the reference step glue): constant Σ and Σ(E) ∝ E^−0.5 in a synthetic medium; s_max ∈ {0.1, 1} mm; 2·10⁵ particles each | \|S_MC(x) − exp(−∫Σ dx)\| ≤ 0.002 abs at 10 depths, and z = \|S_MC − S\|/σ ≤ 4 | self-consistency | CI | A |
| P5 | Every `make_nuclear` function: python twin vs a Warp-CPU float64 test kernel on recorded inputs (the U1 completeness harness covers every function) | Decision 0001 classes: discrete outputs equal; continuous outputs within 1e-12 relative | backend parity | CI | A |
| N1 | Runtime Σ_mass (~~uniform ln E grid, ≥ 50 points/decade~~ **A3 (2026-10-07):** union grid of all ENDF MF3/MT5 and MF6/MT5 energies with a uniform ln E grid ≥ 50 points/decade, lin-lin in E) vs Σ_ENDF at every ENDF node and every node midpoint, 1–150 MeV; water, ICRU tissue (`materials.py` set) | rel ≤ 1e-3 at every point | self-consistency | LV | A |
| V1 | **Source-model reproduction:** the composed per-element σ (after surrogate scaling and composition) vs ENDF MF3/MT5 per element and for the compositions, 1–150 MeV; plus the extension's continuity at 150 MeV | rel ≤ 1e-3 at every ENDF node; extension: σ(150⁺)/σ(150⁻) − 1 within 1e-6 | self-consistency (source model) | LV | A |
| V1b | **Informative:** σ_nonel of p+C (and p+O where present) vs EXFOR raw entries published 1997 or later. The set is D0356 (Auce et al. 2005) plus the entries in the committed manifest. Weighted mean ratio per window [20,40), [40,70), [70,110), [110,160), [160,250] MeV, with the PDG scale factor. | **Not gating.** Ratios and 1 σ archived. Pre-declared known deficiency: LA150 p+C ≈ 227 mb at 100 MeV vs 245–275 mb measured (up to about −18 %). A window without independent data is recorded as an evidence gap. | measured | LV | A |
| V2 | Primary survival S(z) at 100, 150 and 200 MeV in water, MCS and straggling off, default steps (s_max 1 mm, f_E 0.02); N = 2.4·10⁵ per energy (sharded, see Seeds); 1 mm depth bins | TOST, tol 0.003 abs, vs numpy exp(−∫ρΣ dz) along the deterministic CSDA path (same tables), at every 10 mm depth until the CSDA end minus 5 mm | self-consistency | LV | A |
| V2-probe | 150 MeV, same seed (common random numbers): s_max 0.5 mm vs 1.0 mm, and f_E 0.005 vs 0.02; N = 10⁵ each | \|S_a(z) − S_b(z)\| ≤ 0.002 at every 10 mm depth. The row is valid only if σ_Δ(z) ≤ 7e-4 at the deepest depth; otherwise it is inconclusive, which counts as **fail** and is re-run with more shards and a new base | falsification | LV | A |
| V2b | As V2-probe for s_max 0.1 mm vs 1.0 mm on warp-cpu float64, N = 10⁶ | as V2-probe | falsification | LV | **slice B** |
| V3 | Energy balance with nuclear tallies: `initial = step_deposit + cutoff + nuclear_local + escaped + nuclear_escaped_neutron + nuclear_escaped_gamma + nuclear_binding + nuclear_imbalance + truncated + unaccounted` (**A1, 2026-10-07:** `nuclear_imbalance` added; it is signed and never deposited); the shipped per-grid identity of `EnergyBalance` (V3-004; fields `in_grid_mev`, `quantization_mev`, `outside_mev`), `in_grid + quantization + outside = step_deposit + cutoff`, extended with the nuclear-local destination: `in_grid + quantization + outside = step_deposit + cutoff + nuclear_local` per grid; `nuclear_binding` = Σ(Σm_out + M_r − m_p − M_t) recomputed from AME masses per event | relative residual ≤ 1e-12 (python, float64); counters 0. CI: 200 histories at 150 MeV. LV: 2·10⁴ at 150 MeV and 10⁴ at 250 MeV in water / tissue / bone slabs | conservation | CI, LV | A |
| V4 | Event sampler vs the ENDF table: C-12 and O-16 at 100 and 150 MeV, 10⁵ events each (offline, no transport) | Mean yields y_s (s ∈ {n, p, d, α, γ}) within 1 % + 3 σ of ENDF. Per-event Σ y⟨E′_CM⟩ within 5 % of ENDF (C-12 at 150 MeV: 106 MeV). Distortion table archived (~~accepted attempts, rejection fraction~~ **A1 (2026-10-07):** exact P_accept, λ/y, Δ̄_lab, sd(Δ_lab), Δ̄_CM, mean \|Σp_CM\|). | self-consistency | LV | A |
| V4b | **A1 (2026-10-07), regression guard set after the C7 experiments, not validation.** Same events as V4 | Per-species ⟨E′_CM⟩ ratio sampled/ENDF within 1 % for n, p, d and within 6 % for α; γ report-only. \|Δ̄_lab\| ≤ 0.12 E_avail. Exact P_accept ≥ 0.99 at C-12/O-16 100/150 MeV. Per-node ledger closure as P3 on every event. **A4 (2026-10-07):** the ⟨E′_CM⟩ limits carry the statistical term of V4: within 1 % + 3·sem (n, p, d) and 6 % + 3·sem (α), sem from the 10⁵ events. | self-consistency | LV | A |
| V5 | Absolute IDD (r = 20 cm) at 150 and 200 MeV, nuclear on and off, vs TOPAS (opt4, QGSP_BIC_HP, I = 78 eV) and MCsquare (V3-010B batches) | Plateau within 2 %; R80 ≤ 0.5 mm; peak/plateau within 2 %; plateau-integrated ΔIDD (on − off) within 10 % of TOPAS; total deposit within 1 % | independent MC (MCsquare shares σ lineage) | LV | **slice B** |
| V6 | Gottschalk 2015 `tbl:dmlg` grid (177 MeV, Appendix D beam), MeV/g/p | r = 0: 2 %; r ≥ 1 cm: 15 % per point, ≥ 90 % of points pass, median ≤ 10 %, no r-trend (sign test p > 0.01) | measured | LV | **slice B** |
| V7 | Secondary estimators (secondary-p dose, `nuclear_local`, escaped neutral energy) at N = 10⁴…10⁶; grid shift and refinement; f32/f64; batches | Relative-SE slope −0.5 ± 0.05; 1 σ coverage 68 ± 3 % | falsification | LV | **slice B** |
| V8 | Trajectory parity python vs warp-cpu f64 with nuclear on; statistical parity including CUDA | Identical events, species and genealogy ids; continuous columns 1e-10; T12-style statistical tests | backend parity | CI (K = 8) / LV / HR | **slice B** |
| V9 | Deuteron transport tables vs the velocity-scaling identity S_d(E) = S_p(E m_p/m_d), R_d(E) = (m_d/m_p) R_p(E m_p/m_d), E ≥ 2 MeV/u, Bethe source, water and tissue | rel ≤ 1e-3 for S and R at every grid node | theory | CI | A |
| R1 | `nuclear=False` regression = row A16 of `v3-004-acceptance.md`, re-pointed to baseline `f3a1dd62` (amendment 6 there; `A16_BASELINE` already advanced by the housekeeping commit 9448d5e). The lv4 suite runs the gated regression mode (`--mode regression`). | A16 digest (tallies, counters, grids, channel accumulators, residuals, trace) bit-identical to f3a1dd62 at digest seed 20351004, python and warp-cpu f64. **The V3-005A intended-change set for `nuclear=False` is empty: `A16_INTENDED_CHANGE` stays `None` and no `A16_INTENDED_CHANGE_005A` exists.** T1 passes unchanged. The digest (`a16_digest.py`) hashes the grid batch energies, every field of `asdict(energy_balance)`, `counters.as_dict()`, `diagnostics` and `valid`, and the comparison fails on any differing field set, so for `nuclear=False` the nuclear tallies and the counters `majorant_violation`, `nuclear_rejection_limit`, `nuclear_conservation` must not appear in those structures (they exist only in conditional blocks with `nuclear=True`; `TALLY_NAMES`, `COUNTER_NAMES` and `TRACE_COLUMNS` are unchanged). The capability report (`ionmc.simulation.capabilities()`) is **not** part of the A16 digest (neither `a16_digest.py` nor `steps_v4.py` reads it); its changes (the `physics.nuclear` flag and a `nuclear` section) lie outside the digest and are tested separately (row C1 and the capability tests). Other changes outside the digest: the new `ionmc data` import command, registry entries and EXFOR manifest, the nuclear data/table loaders, the deuteron `TransportTables`, new modules and tests. | regression | CI, LV | A |
| X1 | Paired on/off runs, same seed, 150 MeV, 2 × 2·10⁴ histories, all EM physics on | The peak/plateau ratio (plateau = mean IDD over 20–60 mm) decreases with nuclear on, and the paired z of the difference (batch-paired) > 3 | related model | LV | A |
| C1 | Fail-closed configuration and runtime. Every case below either raises `UnsupportedCombinationError` before transport, or ends with `valid=False` and the named counter > 0: unsupported element (no evaluated or surrogate entry); E0 + 6σ_E > 250 MeV; non-proton source; `nist-star` stopping with nuclear; any warp backend with nuclear; missing table; stale table (re-hash mismatch); source pin mismatch; a tally for an unproducible species or generation; and forced `queue_overflow`, `genealogy_overflow`, `nuclear_rejection_limit`, `majorant_violation`, `path_bound_exceeded` | each case raises or is `valid=False` with that counter | capability | CI | A |
| D6 | α local-deposition gate (decision 0041). G = share of the MF6/MT5 α energy (water-weighted, lab frame) carried by α with ASTAR water CSDA range > 0.1 mm, at 150 and 250 MeV. D = (α energy carried by α with range > 0.1 mm per event) × P_event(E) / E, with P_event the non-elastic probability over the full CSDA path in water at that beam energy (computed from the built Σ along the deterministic CSDA path). The builder computes G and D and writes both numbers to the table JSON; the row reads them and the archive carries both. | **Tier 1 (as ratified):** G(150) ≤ 0.05 and G(250) ≤ 0.05 → α deposited locally, no further condition. **Tier 2 (pre-declared fallback, applies only if tier 1 fails):** local deposition is still accepted iff D ≤ 1e-3 AND the 99.9th-percentile lab α energy has an ASTAR water CSDA range ≤ 2 mm at both energies. Passing via tier 2 is recorded in decision 0041 as a documented approximation that forbids sub-millimetre statements about α dose until V3-008 (z = 2 transport), and the ledger status says so. **Neither tier:** the task stops and escalates to the orchestrator, which decides between blocking α-dependent claims and re-scoping (no α transport in slice A). Both tiers were fixed before any G or D value was computed. | theory / model-domain | LV | A |
| E1 | **Exploratory:** python IDD (150 MeV, 10⁵ histories) vs V3-010B TOPAS/MCsquare via `compare_depth_dose.py` | Report only | independent MC | LV | A |

## Seeds

The consumed bases (shared list; `validation/results/transport/README.md`, `v3-003d-acceptance.md` Seeds) are 20261004 … 20401004 and the rehearsal family 2041xxxx. V3-005 declares:

- **lv5 (slice A) qualification base: 20421004.** Row r and shard k use the seed `20421004 + 1000·r_index + k`. `r_index` is fixed in `steps_v5.py` (V2-100: 1, V2-150: 2, V2-200: 3, V2-probe-s05: 4, V2-probe-fE: 5, V3-LV: 6, X1: 7, E1: 8, V4: 9). Shards of one row are independent samples and are pooled as batches.
- **Rehearsals** of V3-005A use the 2043xxxx family only. They are never qualification evidence.
- **Re-runs.** A failed or inconclusive qualification base is consumed. The next fresh base is the next unused value of 20441004, 20451004, … The new base is declared in an amendment *before* it is observed.
- **Slice B** (lv5b/hr5) declares its bases in an amendment at the start of V3-005B, from the same fresh sequence.
- **CI tests** use fixed small seeds documented in the tests; they are not qualification evidence. The A16/R1 digest seed is 20351004 (fixed, not a sample).

## Execution under the single-process directive (temporary)

- Every LV step runs with `cpu_workers=1` as one compute job, with `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1` and `IONMC_SINGLE_PROCESS=1`.
- Steps are run as `run_suite.py --only <step> --step-timeout 3300`.
- A row that would exceed 3300 s is **sharded** (more steps of fewer histories, distinct shard seeds). Its history count and criterion never decrease.
- The worker records the measured nuclear-on throughput before the first lv5 run. The planning value at commit f3a1dd6 is 77.9 hist/s at 150 MeV and 68.8 hist/s at 250 MeV (EM only, no tallies, python, one process).

## Deferred multiprocess checks

The temporary single-process directive is an execution mode, not a property of the engine: the python pool path is implemented for nuclear runs (the per-history stack is history-local and the partitioning is unchanged). The lv5 suite runs `cpu_workers = 1`. These checks are archived as **DEFERRED**. They are not passed and are not part of a conformant qualification:

- 1-vs-N-worker partition invariance of nuclear runs (python pool);
- any `cpu_workers > 1` LV row.

A summary carrying them is `pass` only for the executed rows and is `conformant: false` until the operator lifts the directive and the deferred rows run.

## Evidence archive

| Row(s) | Archive (`validation/results/transport/...`) | Head SHA | Seed base | Result |
|---|---|---|---|---|
| P1–P5, V3-CI, V9, R1, C1 | CI run; R1 also in `lv5-b84fdf3-single-process.json` | b84fdf3 | fixed CI seeds | pass |
| N1, V1, V1b, D6 | `lv5-b84fdf3-single-process.json` | b84fdf3 | n/a (deterministic) | pass (D6: tiers failed, ceiling holds); non-conformant by code |
| V2 (100/150/200), V2-probe | `lv5-b84fdf3-single-process.json` | b84fdf3 | 20421004 | pass (imported partials → non-conformant by code, Amendment 5) |
| V3-LV, V4, V4b, X1, E1 | `lv5-b84fdf3-single-process.json` | b84fdf3 | 20421004 | pass (E1 exploratory); non-conformant by code (diagnostic mode) |
| V2b, V5–V8 | slice B | — | declared in B | — |
| V8-LV, R1-nuc, V5, V2b, V7 (scan, shift, replicates) — consumed base | `lv5b-b40d812-single-process.json` | b40d812 | 20441004 (consumed, Amendment 11) | not evidence: V7 replicate step failed under the pooled rule (232/360 = 0.644); V5 fail (peak/plateau +12.6 % / +11.9 %, ΔIDD −20.8 % / −13.4 % vs TOPAS; total deposit, plateau, R80 pass; κ not used); other executed steps pass |
| V8-LV, R1-nuc | `lv5b-7aae5bb-single-process.json` (steps 02, 03) | 7aae5bb | 20471004 | pass (trajectory parity: maximum continuous difference 0.0, tolerance 1e-10; A16 bit-identical, empty intended-change set); non-conformant by code (diagnostic mode) |
| V5 | `lv5b-7aae5bb-single-process.json` (steps 04–08) | 7aae5bb | 20471004 | **fail** vs TOPAS: peak/plateau +12.5 % / +12.3 %, ΔIDD −18.5 % / −14.7 %; plateau, R80, total deposit (−0.30 % / −0.54 %) pass; κ not used; attributed to the absent hadronic elastic channel (V3-005C); imported partials → non-conformant by code |
| V2b | `lv5b-7aae5bb-single-process.json` (steps 09–11) | 7aae5bb | 20471004 | pass (max abs diff 8.7e-6 vs 2e-3, σ_Δ 1.4e-5 vs 7e-4, 1e6 per variant); imported partials → non-conformant by code |
| V7 N-scan, grid shift/refinement, f32/f64 (warp-cpu) | `lv5b-7aae5bb-single-process.json` (steps 12, 13) | 7aae5bb | 20471004 | pass (slopes −0.495 / −0.518 / −0.480; shift 5.4e-4, refinement 1.4e-13) |
| V7 replicate coverage (Amendment 13) | `lv5b-7aae5bb-single-process.json` (steps 14–24) | 7aae5bb | 20471004 | **fail**: `escaped_neutral` paired gate (Clopper–Pearson lower bound 0.6445 < 0.6464; single gate passes); `sec_p` and `nuclear_local` pass both gates; calibration step 23 pass; no further base (Amendment 11); imported partials → non-conformant by code |
| 1-vs-N worker partition | `lv5b-7aae5bb-single-process.json` (step 25) | 7aae5bb | — | deferred (single-process directive) |
| V8 statistical parity | `hr5-7aae5bb-single-process.json` (steps 01–06) | 7aae5bb | 20451004 | all CPU/CUDA pairs pass; python:cpu64 `idd`, `sec_p` pass; python:cpu64 `nuc_local` **inconclusive** (T12 sparse-profile rule: 60/161 bins supported, hull mass 0.259 < 0.5); imported partials → non-conformant by code |
| V7 f32 vs f64 (CUDA) | `hr5-7aae5bb-single-process.json` (step 07) | 7aae5bb | 20451004 | pass (\|z\| ≤ 2.55 of 3.0) |
| V6, D2, E1-B | — | — | — | re-assigned to V3-005C (Amendment 7(a)) |

**Correction note (2026-10-10, Codex REVIEW-7d749e6a).** In the V5 row of the table above, "attributed to the absent hadronic elastic channel (V3-005C)" is to be read as "hypothesis: the absent hadronic elastic channel (to be tested by V11/V5 in V3-005C; the exploratory X cases show elastic is not the only contributor)". The row text is left unchanged.

**Execution of Amendments 11–13 (recorded 2026-10-09; the amendment texts above are unchanged).**
- Amendment 11: the re-run at base 20471004 was executed at head 7aae5bb (archive `lv5b-7aae5bb-single-process.json`) with the replicate design that Amendment 13 later fixed; row V7 failed there, so no further base is taken and the failure is handed to V3-005C as declared.
- Amendment 12: superseded by Amendment 13 before any V7 observation at 20471004; its 300-replicate rule was never executed.
- Amendment 13: executed at 7aae5bb (steps 14–24, 7200 replicates, calibration step passed); `sec_p` and `nuclear_local` passed both gates, `escaped_neutral` failed the paired gate (m 0.6586, lower bound 0.6445 < 0.6464) and passed the single-interval gate, so row V7 fails.

## Amendments

**Amendment 1 (2026-10-07) — event model (rows P2, P3, V3, V4; new V4b).** Reason: the C6/C7 build (table `99b51c2c…`, 100 points/decade, 10⁵ events/node/iteration, 1534 s) showed that the frozen event model cannot reproduce the ENDF mean yields. Contrary evidence, verbatim from the C7 record: "ENDF Σ y⟨E'⟩ at 150 MeV is 103–107 MeV against ≈138–140 MeV CM kinetic energy minus ≈15–20 MeV separation energies. With λ = ENDF yield about half of all attempts are rejected and post-acceptance yields are 30–50 % low (O-16 ratios 0.57–0.97, C-12 0.53–0.99). Raising λ pushes the fraction of events exhausting 64 attempts above 1 % (converged C-12 150 MeV solution exhausts 26 %; exhaustion sets nuclear_rejection_limit and invalidates a run). Below 60–80 MeV yields are unreachable; at 20–30 MeV λ saturates at the cap 20. … Converged nodes: C-12, N-14 none; O-16 ≥ 125 MeV; Al-27, Si-28, P-31, Ca-40 ≥ 80 MeV. Independent-sample residuals (n,p,d,α,γ): C-12 100 MeV [−0.20,−0.19,−0.21,−0.21,−0.03]; C-12 150 MeV [−0.12,−0.11,−0.14,−0.04,−0.005] (accepted fraction 0.991); O-16 100 MeV [−0.027,−0.034,−0.07,+0.043,−0.008]; O-16 150 MeV [0.010,0.000,−0.010,−0.005,−0.007]." Planner experiments (λ = y, 10⁵ events per case): independent sampling of the inclusive spectra overshoots the available energy in 44–64 % of events (CM Δ sd 89–145 MeV); every per-event-conserving repair loses 34–64 % of Σ y⟨E'⟩ (one-factor KE rescale 0.50–0.66, rejection 0.36–0.48, sequential truncation equivalent); momentum by difference gives a residual recoil of 15–38 MeV against 3.1–3.6 MeV evaluated; the C7 design deposited E* ≈ 0.5 E_avail (46–70 MeV) locally against 9–12 MeV implied by the evaluation. Per-event energy and momentum conservation is therefore withdrawn (decision 0041 §3, amendment 2026-10-07). New model: floor+Bernoulli multiplicities n_s = ⌊λ_s⌋ + [u < λ_s − ⌊λ_s⌋] with λ solved at build time by exact enumeration; acceptance only on the existence of the residual nuclide in AME2020; the residual receives the ENDF mean heavy-recoil energy (local); E* is not deposited (de-excitation is already in the ENDF γ); the signed per-event difference is tallied as `nuclear_imbalance` and never deposited. Rejected alternatives (orchestrator): a momentum-balancing residual (its lab kinetic energy would be negative for the typical net CM momentum, i.e. the configuration is unphysical); a per-node product-energy scale κ making Δ̄_lab = 0 (would break the V4 Σ y⟨E'⟩ 5 % criterion by ≈ −8 %; it stays the pre-registered fallback if V5 (slice B) fails on total deposit); 256 inverse-CDF bins (×4 table memory). Row changes: P2 adds the floor+Bernoulli case; P3 becomes ledger closure; V3 adds `nuclear_imbalance`; V4 tolerances unchanged, distortion table redefined; new row V4b whose tolerances were set after the planner experiments and are regression guards. Expected from those experiments: Σ y⟨E'⟩ +0.4…+0.7 %; ⟨E'⟩ α +3…+5 %, γ +14…+24 % (64-bin uniform-in-bin sampling); Δ̄_lab/E_avail −0.070…−0.098 (products carry more energy than the event has; in water at 150 MeV ≈ 1.5 % of the beam energy). The C7 table and its numbers stay archived as superseded evidence.

**Amendment 2 (2026-10-07) — row D6 outcome.** Both tiers failed on the C7 table (150 MeV: G 0.5718, D 3.78e-3, 99.9th-percentile α 60.3 MeV, range 0.258 g/cm²; 250 MeV: G 0.8038, D 1.10e-2, 117 MeV, 0.852 g/cm²; P_event 0.160 / 0.321). The planner's estimate with the amended event model is higher, because the rejection sampler suppressed the hard α tail (150: G 0.672, D 5.55e-3, 99.9th-percentile 128 MeV → 0.995 g/cm²; 250: G 0.850, D 1.42e-2, 213 MeV → 2.51 g/cm²). Orchestrator decision (pre-declared stop rule): α stay locally deposited as a documented approximation; D6 is recomputed with the amended sampler and archived; the decision stands only if D ≤ 2e-2 and the 99.9th-percentile α range ≤ 3 g/cm² at both energies (ceiling set with the estimate above known), otherwise the task stops again. Every nuclear run reports the tally `nuclear_alpha_local`. No statement about α dose, α LET, or dose structure attributable to nuclear secondaries on a scale below 1 cm (150 MeV) / 2.5 cm (250 MeV) is made until V3-008 re-runs D6 with transported α (paired local-vs-transported IDD and lateral profiles at 150/250 MeV). Straight-line CSDA α placement was rejected as a re-scope that duplicates V3-008.

**Amendment 3 (2026-10-07) — rows N1, V1 grid.** The uniform-ln E table (100 points/decade) failed N1/V1: above 10 MeV max rel error 1.35e-3 (water), 0.8–1.3e-3 tissues; at threshold kinks up to 0.15 (muscle 4.8 MeV), water 0.019 at 7 MeV; V1 per target up to 9.3e-2 (Si-28 at 3 MeV); extension continuity 1.5e-10 (pass). MF3/MT5 of every target is one lin-lin region (INT=2; nodes C-12 33, N-14 143, O-16 32, Al-27 180, Si-28 170, P-31 191, Ca-40 194), so the runtime grid becomes the union of all ENDF MF3/MT5 and MF6/MT5 energies in [1,150] MeV with a uniform ln E grid of ≥ 50 points/decade over [1,250] MeV (646 nodes at 50/decade), interpolated lin-lin in E with a fixed-step binary search (decision 0041 §1 amended). Tolerances unchanged. Relaxing the tolerance was rejected: it would hide genuine table errors.

**Amendment 4 (2026-10-07) — C7b build findings (rows V4b, N1/V1 domain, capacity bound).** The first C7b build (table `12b1ac85…`, 604 grid nodes, 12.8 s) met N1 (2.8e-15), V1 (2.2e-16), V4 (yield ratios 0.988–1.003 within 1 % + 3σ; Σ y⟨E′⟩ + recoil ratios 1.002–1.009) and the D6 ceiling (150 MeV: G 0.668, D 5.48e-3, 99.9th-percentile α range 0.988 g/cm² count-weighted / 1.109 energy-weighted; 250 MeV: G 0.849, D 1.420e-2, 2.529 / 2.827 g/cm²; tiers 1 and 2 false, ceiling true). Three findings needed orchestrator decisions, recorded before the rebuild: (a) the residual-existence rule rejected the full break-up (Z_r, A_r) = (0, 0) (p + C-12 → p + 3α), so the exact P_accept of C-12 fell to 0.256 at 17 MeV and the build failed its own P_accept ≥ 0.5 rule; the empty residual with M_r = 0 is physically valid and is always allowed (decision 0041 §3). (b) Al-27 and P-31 at 1 MeV and Ca-40 at 3.44 MeV have σ > 0 where the evaluation tabulates zero yields (the first MF6 energy equals the threshold with yield 0; the first positive total yield is at 2.0 MeV for Al-27, 1.307 MeV for P-31 and 4.0 MeV for Ca-40); between the MT5 threshold and that first positive-yield energy every product's yield is held at its value there (distribution rows were already constant below their first energy; recorded per target as `yield_extended_below_mev` in the table JSON), so that N1/V1 are evaluated at every point and no σ is zeroed. The affected σ is 13 mb (Al-27) and 2 mb (P-31) at 1 MeV, where a proton stops within tens of µm. (c) V4b as written had no statistical term although V4 has one; the first build gave an n ⟨E′⟩ ratio of 1.0122 for O-16 at 100 MeV against the bare 1 % limit (the other 19 species/case ratios for n, p, d were within 1 %; α 1.031–1.045; γ 1.14–1.25 report-only); the row now reads 1 % + 3·sem (n, p, d) and 6 % + 3·sem (α); the raw 1.0122 is preserved here. (d) The per-event energy sum bound T_sum,max = 1030 MeV exceeded the 500 MeV stopping table, so the path bound is re-derived as B_L = 1.25·(R(E_hi) + Σ_s N_s,max·R_s(T_lab,max,s)) with per-species counts and single-particle energies (no range extrapolation); the per-particle energy bound is max_s T_lab,max,s ≤ 500 MeV (decision 0041 §5 amended). **Rebuild outcome (flagless, table `2b8d94cf4a82ef8a5119ddb8c84d4a2c11bd3511eb6672bd210562730cc0cf1a`, 604 nodes, 12.8 s, P_accept ≥ 0.9864, 0 non-converged λ nodes):** N1 2.8e-15, V1 2.2e-16, continuity 1.48e-10; V4 pass (yield ratios 0.988–1.003, Σ y⟨E′⟩ + recoil 1.002–1.009); V4b pass under the amended rule (n/p/d ratios 0.9923–1.0122 with sem 0.0025–0.0063; α 1.031–1.045; γ 1.14–1.25 report-only; Δ̄_lab/E_avail −0.072…−0.101; P(Δ<0) 0.51; mean |Σp_CM| 388–515 MeV/c; ledger closure ≤ 4.0e-12 MeV); D6 as above, both tiers false, ceiling true; `transport_energy_bound_mev` 262.08 (p n_max 3 at 262.08 MeV, d n_max 1 at 246.23 MeV), B_L(250 MeV, water) = 2489 mm with the proton range standing in for the deuteron term. **Rebuilt 2026-10-08 after Codex review as `dfee19d303c7fdd28d2080ebbb37a01b27aa6569cca846ad317a2d91a7a437a3` (builder version 3; identical physics, same N1/V1/V4/D6 numbers; new per-history energy bound 3154.4 MeV for the accumulators, per-particle bound 262.08 MeV unchanged).** **Rebuilt again 2026-10-08 after the second Codex review as `3bcf146e38dd2b5581bd1ff245c127a7d059e6789421a5d3c6760974f2784504` (builder version 4, schema 3: the id now covers the npz digest; qualification flags and bounds live in the npz; identical physics and numbers).** **Rebuilt 2026-10-08 after the third Codex review as `00e8031d5f67704f256e904ebab89b0c2a75845bcf1887c29e1518bffaafe5fb` (builder version 5: the id is the sha256 of the canonical sidecar, which carries the npz digest; identical physics and numbers).**

**Amendment 5 (2026-10-08) — conformance of sharded rows under the single-process directive (Codex REVIEW-3a1a1c8b finding 1).** The V2 and V2-probe shards are run as separate host-runner invocations (3600 s limit) and combined from imported partials. The code cannot verify the protected host-runner records, so a combine that used any imported partial is `conformant: false` by code with the reason `imported_partials_unverified_by_code`, whatever the manifest says; its `pass` verdict on the row criterion is still evaluated and archived. The orchestrator's validation record lists the host-runner run ids whose protected stdout carries the `PARTIAL <name> <sha256>` lines matching the manifest; that record, not the code, is the attestation. Sharded rows become conformant only when all shards and the combine run in one invocation, which requires the operator to lift the single-process directive (the lv5 summary is already non-conformant for the deferred partition row). History counts and criteria are unchanged.

**Amendment 6 (2026-10-08) — slice-B seeds.** lv5b qualification base **20441004**, hr5 base **20451004**. Row r and shard k use `base + 1000·r_index + k`, with r_index fixed in `steps_v5b.py` (throughput 1, V8-LV 2, R1-nuc 3, V5-150-on 4, V5-150-off 5, V5-200-on 6, V5-200-off 7, V2b 8, V7-scan 9, V7-shift 10, V7-replicates 11 (replicate j uses shard k = j), V6 12 and E1-B 13 reserved for V3-005C; hr5: V8-stat 1, V7-f32 2, V6-CUDA 3 reserved for V3-005C). Rehearsals of V3-005B use the **2046xxxx** family only and are never evidence. A failed or inconclusive base is consumed; re-runs take the next unused value of 20471004, 20481004, …, declared in an amendment before observation. The κ fallback re-run (if triggered) takes the next unused value. CI seeds are fixed in the tests.

**Amendment 7 (2026-10-08) — slice-B scope and decision rules, fixed before any slice-B result.** (a) *Split:* V3-005B evaluates V2b, V5, V7, V8 and R1/C1 on the warp backends; V6, the p-p elastic model (D2) with its new rows and the E1-B exploratory re-run are evaluated by **V3-005C** (dated amendment there; new rows P6, P7, V10, V11 are proposed, not yet frozen). (b) *V5 gating engine:* TOPAS (opt4, QGSP_BIC_HP, I = 78 eV) gates; MCsquare (shared σ lineage; R80 already 0.67 mm from TOPAS in V3-010B) is report-only. New reference runs: TOPAS full and EM-only, MCsquare nuclear on and off, at 150 and 200 MeV, r = 20 cm, 1e5 primaries, 0.5 mm bins, 3 seeds each. (c) *κ fallback rule:* with D = (dep_ionmc − dep_TOPAS)/dep_TOPAS per energy (TOST, 1 %): if V5 total deposit passes at both energies, κ is not used; if it fails with D > 0 and |D| ≤ 1.5 % at either energy, the per-node product-energy scale κ (Δ̄_lab = 0; expected ≈ −0.5 % deposit) is built, V5 is re-run on a fresh base and κ is adopted iff V5 then passes, with V4's Σ y⟨E′⟩ criterion amended to a documented failure (≈ −8 %); if it fails with D < 0 or |D| > 1.5 %, κ cannot be the remedy and the task stops. (d) *Secondary-proton non-elastic interactions* are included (python and warp together, after V8 passes on the frozen slice-A model) with an A16-style intended-change record for nuclear-on outputs (expected 0.0037 / 0.0096 events per history at 150 / 200 MeV); deuteron interactions stay excluded (0.0009 / 0.0024 per history); p + C/O nuclear elastic scattering stays deferred, so no lateral-halo claim is made in V3-005B. (e) *V7 coverage protocol:* ≥ 30 replicates × ≥ 10 bins (≥ 300 intervals) on warp-cpu f64; f32 vs f64 on CUDA in hr5. (f) *R1 extension:* until the secondary-interaction change lands, the python nuclear-on outputs must be bit-identical to the lv5 archive at b84fdf38. (g) *C1:* the case "any warp backend with nuclear raises" is replaced by the forced-counter cases on warp-cpu (same counters as python). (h) *CUDA:* the nuclear kernel caps the chunk at 2^14 histories (4 KiB stack scratch per thread). (i) *D9* evaluates Tripathi-light and Geant4 BGG/Barashenkov σ against the post-1997 EXFOR set offline with the pre-registered window rule; BGG is report-only (lineage with TOPAS). Tolerances of the frozen rows are unchanged.

**Amendment 8 (2026-10-08) — V7 and V8 execution details declared before any slice-B result.** V7 N-scan uses 100 batches per N (slope standard error ≈ 0.02 against the frozen ±0.05); the grid-shift case passes if the estimator totals differ by ≤ 1e-3 relative and the refinement case by ≤ 1e-6 relative (new tolerances for the two layout cases the frozen row names without a number); replicates ≥ 30 × ≥ 10 bins (Amendment 7). Because the warp kernels assign batches as `h mod n_batches`, V7 and the hr5 V8 samples run contiguous history blocks as partials of a two-batch configuration and sum the rows (verified exact against a full run at rehearsal scale). hr5 V8 statistical parity: python 2.4e4 histories in two shards, warp samples 1e6 each; pairs python:cpu64, cpu64:cuda32, cpu64:cuda64, cuda32:cuda64 with the T12 comparison. V2b runs 1e6 per variant in two shards with common random numbers. lv5b step numbering follows `steps_v5b.py` (V2b = steps 08–10; the V3-005C rows keep r_index 12, 13 and hr5 3). Tolerances of the frozen rows are otherwise unchanged.

**Amendment 9 (2026-10-08) — V5 execution details declared before any V5 result.** Reference cases: TOPAS full (six V3-010B modules, I = 78 eV, zero-size pencil beam) and EM-only (`g4em-standard_opt4` only), MCsquare nuclear on and off (`Simulate_Nuclear_Interactions` and the three `Simulate_Secondary_*` flags), 150 and 200 MeV, 1e5 primaries, three independent seeds per configuration (`20270000 + 1000·row + k`), scoring the IDD in 0.5 mm depth bins over a 400 × 400 mm square (contains the r = 20 cm circle; the corners carry negligible dose) in a 320 / 330 mm deep water box; the ionmc side is the lv5b V5 step (warp-cpu f64, 20 batches, same transverse extent). `compare_idd_v5.py`: ionmc metrics from the pooled batch curve with a delete-one jackknife SE (df 19), engine metrics as the mean over seeds (df n−1), TOST on the 90 % interval with Welch df; TOPAS EM-only serves as the "off" reference for ΔIDD; every curve is cut to the ionmc depth range (1.1 R) so that "total deposit" is the same quantity; a run whose in-grid total lies outside 0.85–1.02 of the beam energy is refused as a unit error. κ rule clarification: if one energy fails with D > 0 and |D| ≤ 1.5 % while the other fails for another reason (D < 0 or |D| > 1.5 %), the task stops (κ is built only when every failing energy qualifies). MCsquare remains report-only (Amendment 7).

**Amendment 10 (2026-10-08) — Codex REVIEW-7516c131 at 88986f17, before any slice-B qualification result.** (a) *Table domain for secondaries:* a secondary proton whose kinetic energy exceeds the nuclear table's 250 MeV endpoint is transported electromagnetically without nuclear sampling until its energy is inside the table domain (no clamping to the endpoint); both backends identically; counted in the nuclear diagnostics. (b) *V5 lineage:* `compare_idd_v5.py` derives the energy and the full/EM-only or nuclear on/off label of every engine run from the fingerprint-verified normalized input (`compare_batches.run_fingerprint`: native history counts, output summaries, request/manifest agreement, clean source SHA, runtime identity, input equivalence, physics-module list), never from `case.json`; replicates must share the verified configuration; the comparator exits non-zero on a failed gating verdict. (c) *V5 enforced step:* lv5b gains the step `v5-compare` which runs the comparator on the four ionmc partials and the reference runs present in the snapshot cache and whose verdict is row V5's verdict (the four `v5-ionmc` steps remain producers). (d) *D9 scope:* the Geant4 BGG/Barashenkov evaluation of Amendment 7(i) is deferred to V3-005C: its parameterisation is not in the repository or the cache; the acquisition to make (Geant4 11.4.x `G4ComponentBarNucleonNucleusXsc.cc`, Geant4 Software License, role exploratory/report-only because of its lineage with TOPAS) is recorded; D9 in V3-005B covers LA150 versus Tripathi-light only. Tolerances unchanged.

**Amendment 11 (2026-10-08, revised the same day after Codex REVIEW-db4474cd, before the re-run) — V7 replicate coverage: base 20441004 consumed and the re-run declared.** At base 20441004 (host run RUN-20261008T132316Z-11c02312, code b40d8121) the V7 replicate-coverage step failed: 232 of 360 one-sigma intervals covered the 1e6-history reference, coverage 0.644 against the frozen 0.68 ± 0.03 (every other executed lv5b step passed; the archive is preserved as `validation/results/transport/lv5b-b40d812-single-process.json` with `pass = false`). *Uncertainty, cluster-aware.* The 360 indicators are not independent: the 12 depth bins of a replicate are scored from the same histories and every replicate is compared with the same finite reference, so the replicate is the independent unit. The step did not record per-replicate counts at the time; a sandbox reproduction of the consumed run through the step's command line (same code path, base 20441004, 30 replicates, reference shard 30; bitwise identical to the host run: 232/360 with the same per-bin counts) gives per-replicate covered bins [9, 10, 10, 5, 6, 8, 10, 5, 7, 5, 10, 6, 10, 8, 10, 5, 9, 10, 9, 5, 7, 8, 6, 6, 7, 8, 7, 9, 7, 10], cluster standard error 0.028 and 90 % interval [0.596, 0.692] (t₂₉), which contains both the frozen target 0.68 and the nominal coverage of one-sigma intervals built from the sample standard error of 20 batches, P(|t₁₉| ≤ 1) = 0.670; the reference contributes a variance ratio of 0.009 (its own standard error is negligible). A second sandbox diagnostic, run by mistake at the V3-003 base 20391004 (the module default when the step is called without its command line; neither a slice-B base nor the declared re-run base, never evidence) gave 224/360 = 0.622 with cluster standard error 0.029 and per-replicate counts [9, 7, 8, 8, 5, 7, 6, 6, 11, 5, 9, 6, 6, 9, 8, 6, 7, 10, 7, 6, 9, 9, 6, 8, 12, 7, 11, 6, 5, 5]; it is disclosed because it was observed, and it is not used. Both observations lie below the nominal 0.670 (by 0.9 and 1.6 cluster standard errors; pooled 0.633 ± 0.020, 1.8 standard errors below), so the failure is recorded as *unexplained*, not as fluctuation: a small under-coverage of the batch standard error for the secondary-proton dose (skewed batch distributions at 500 histories per batch, the nearly empty distal bin at the end of range, which is bin 12 with reference 2e-4 and 17/30 and 12/30 hits) is a hypothesis, not a finding. *Decision (frozen before the re-run).* Per Amendment 6 the base 20441004 is consumed; the whole lv5b suite is re-run at the next base, **20471004**, at the code SHA that carries this amendment (`V5B_QUALIFICATION_SEED_BASE`; 20441004 joins the consumed list the code refuses). For the re-run the V7 replicate count is 90 (Amendment 8 fixes "≥ 30"; 90 × 1e4 histories, the same 1e6 reference, 12 bins, 20 batches), which lowers the expected cluster standard error to about 0.016 (0.028 × √(30/90)); the step now records per-replicate counts, the cluster standard error, the 90 % interval, the nominal coverage and the reference variance ratio. The gate is unchanged: pooled coverage within 0.68 ± 0.03, at least 300 intervals and 10 used bins. Calibration of that gate, stated so that the outcome can be read: with a cluster standard error of 0.016 the gate rejects a correct estimator (true coverage 0.670) on the low side with probability about 0.11, and passes an estimator with true coverage 0.640 with probability about 0.27. Row V7 is judged by the re-run alone; the consumed run, the diagnostics and the pooled coverage are reported next to it. If the re-run fails, row V7 fails, the task records the scientific failure (V2-NUM uncertainty-coverage evidence not established for the secondary-proton estimators) without a further base, and V3-005C receives a pre-registered re-test with a replicate-level criterion. V5 is re-measured with the new seeds; the consumed run's V5 verdict (fail on peak-over-plateau and ΔIDD, pass on total deposit, plateau and R80, κ not used) is preserved as a second independent observation of the attributed failure. hr5 keeps base 20451004.

**Amendment 12 (2026-10-08, after Codex REVIEW-d3f216ca, before any V7 observation at base 20471004) — V7 replicate coverage judged by a replicate-level equivalence test.** The frozen row V7 states "1 σ coverage 68 ± 3 %" for the batch standard error of the secondary estimators; executed as a point-estimate gate it passes an estimator with true coverage 0.640 with probability about 0.27 (Amendment 11), so a pass would not establish the claim. The exact expectation of the implemented interval rule (|mean − reference| ≤ sample standard error of 20 batches) for a correct estimator is P(|t₁₉| ≤ 1) = 0.670, not the Gaussian 0.683; the shared 1e6-history reference (variance ratio 0.009, Amendment 11) lowers it by at most 0.003. The row is therefore judged, for the re-run at base 20471004 and the code SHA carrying this amendment, by the following pre-registered rule, which is stricter in inference than the frozen gate and centred on the exact nominal: 300 replicates of 1e4 histories (reference shard 300; seeds inside the row's block), 12 bins, 20 batches; f_j = covered bins of replicate j / bins used; m = mean(f_j), s = sd(f_j)/√300, t = t₀.₉₅,₂₉₉; **pass iff the 90 % interval [m − t·s, m + t·s] lies within [0.640, 0.700]** (two one-sided tests at α = 0.05 against the margins 0.670 ± 0.03), with at least 300 intervals and 10 used bins as before. Calibration (from the measured replicate spread sd(f_j) = 0.155, expected s = 0.009): an estimator with true coverage 0.640 or 0.700 is accepted with probability at most 0.05; an estimator at the nominal 0.670 passes with probability about 0.91; one at the consumed run's 0.644 fails with probability above 0.95. The former gate (pooled coverage within 0.68 ± 0.03) is computed and reported as the legacy point gate; it no longer decides. A pass establishes the V2-NUM uncertainty-coverage evidence for the secondary-proton estimators on warp-cpu float64 at 150 MeV; a fail is row V7's failure (no further base, Amendment 11). This is a criterion change recorded before observation (the scientific-event service being unavailable, this paragraph and the decision 0041 addendum are the record); the consumed run at 20441004 remains a failure under both rules.

**Amendment 13 (2026-10-08, after Codex REVIEW-be30e621, REVIEW-873ab9cd, REVIEW-edb9970b, REVIEW-f1b32614, REVIEW-852837a3, REVIEW-22711e6f, REVIEW-97a3377e and the calibration study, before any V7 observation at base 20471004) — V7 replicate coverage for all three frozen estimators, judged by a reference-free paired design with exact finite-sample bounds.** Supersedes the execution details of Amendment 12 (the margin ±0.03 and the total false-acceptance level ≤ 0.05 stay; the centre is the exact nominal coverage of the implemented rule, as Amendment 12 defined it). *Estimators:* the secondary-proton dose profile (`sec_p`, 12 depth bins), the local nuclear deposit profile (`nuclear_local`, 12 bins, `nuc_local_dose` channel) and the escaped neutral energy (`escaped_neutral`, scalar), as the frozen row names them. *Replicates:* eight shard steps `v7-rep-s0..s7`, each one warp-cpu float64 run of 9e6 histories at 150 MeV in the coarse 12-bin geometry evaluated in 18 000 contiguous blocks of 500 histories (block estimates as in the N-scan; every history has its own stream, so blocks are independent samples; the effective configuration used to run the blocks takes the smallest batch count that satisfies the fixed-point accumulator guard, 8 for a shard); replicate j of a shard is the contiguous group of 20 blocks, giving 900 replicates of 1e4 histories per shard and R = 7200; per replicate and bin the mean and the sample standard error over its 20 blocks are stored in a hash-bound, attested partial (shard seeds k = 0..7 of r_index 11). *Paired rule, no reference:* the replicates form P = 3600 fixed disjoint pairs, replicate j of shard s with replicate j of shard s + 4 (s = 0..3); for pair p and used bin b, hit_pb = |mean_jb − mean_kb| ≤ √(sem_jb² + sem_kb²), **except that a pair-bin in which either replicate has a zero standard error (identical blocks, in practice all zero) is a miss**: an undefined standard error cannot certify coverage (V2-NUM); such pair-bins are counted and reported. f_p = hits_p / bins used. The pairs are independent and identically distributed by construction (disjoint pairs of independent replicates; nothing is shared and no reference enters), so no reference error, box, bootstrap or shift has to be bounded. *Bins used* (profiles), fixed by the independent 1e6-history `v7-rep-ref` run (seed k = 8, 20 blocks) before any replicate is looked at: ref_mean_b > 0, ref_spread_b > 0, and the expected per-replicate relative standard error r_b = (SEM_ref,b / ref_mean_b)·√(n_ref / n_rep) < 0.5 (rel-SEM_ref < 0.05 at n_ref = 1e6, n_rep = 1e4): a replicate's batch standard error is a usable uncertainty only where its expected relative error is well below 100 %; in sparser bins most 1e4-history replicates record nothing and their standard error is undefined. The cut 0.5 lies in the wide gap seen in rehearsals between the resolvable bins (r_b ≤ 0.15) and the distal bin (r_b ≈ 1.0), so no bin's inclusion is decided by reference noise; the mask, r_b and the excluded bins with their reason are recorded; at least 10 used bins are required (fail closed). The scalar has one bin. *Nominal and region:* under Gaussian batch means with correct standard errors the paired statistic (X_j − X_k)/√(s_j² + s_k²) is Student t with 38 degrees of freedom, so the nominal one-sigma paired coverage is c₀ = P(|t₃₈| ≤ 1) = 0.6764 and the equivalence region is **[0.6464, 0.7064]** (c₀ ± 0.03). *Link to single-interval coverage:* with a standard-error scale κ, single-interval coverage is P(|t₁₉| ≤ κ) and paired coverage P(|t₃₈| ≤ κ), both strictly increasing; the paired region corresponds to κ ∈ [0.939, 1.065], i.e. to a single-interval coverage band [0.6405, 0.6997], so the paired gate tests the same calibration of the batch standard error as the frozen row's "68 ± 3 %"; the single-interval coverage against the 1e6 reference (m_ref1e6) is reported as a diagnostic. *Single-interval gate (added after Codex REVIEW-25b77cf0 and corrected after REVIEW-0af7808c, before any observation):* in addition to the paired gate, each estimator is judged on the coverage of its own reported mean ± standard-error interval against a held-out reference: H = shards 6–7 (1800 replicates, 1.8e7 histories) gives μ_H,b and SEM_H,b = sd(H means)/√1800; for the 5400 evaluation replicates of shards 0–5 and a candidate true mean μ_H,b + δ·SEM_H,b, hit_jb(δ) = |mean_jb − μ_H,b − δ·SEM_H,b| ≤ sem_jb (the reported half-width only; a zero sem_jb or SEM_H,b is a miss); g_j(δ) = hits_j(δ) / bins used (the same reference-fixed mask); given H the g_j(δ) at any fixed δ are exactly i.i.d. Nominal (one-sample, Gaussian batch means): c₁ = P(|t₁₉| ≤ 1) = 0.6701, region [0.6401, 0.7001]. Reference error: the box over the used bins with half-width Z_box·SEM_H,b, Z_box = max(Z_t, Z_boot) with Z_t the Student-t quantile at 1 − 0.01/(2B) with 1799 degrees of freedom and Z_boot the empirical 99 % quantile of the simultaneous standardized error over 5000 bootstrap resamples of the actual 7200 replicate means (seeded from the run seed base, recorded); the exact per-bin extrema of the coverage over the box by the interval-count sweep give m_lo and m_hi (and K_min, K_max for the scalar's hit count). Bounds, valid uniformly over the box by monotonicity: if the box contains the true shift δ*, then K_min ≤ K(δ*) and m_lo ≤ m(δ*) ≤ m_hi, so a one-sided bound applied to K_min / m_lo (lower) and K_max / m_hi (upper) is at least as conservative as the same bound at the unknown δ*, where the hits are i.i.d. given H: scalar → Clopper–Pearson at α_tost = 0.04 per side on K_min and K_max; profiles → the fixed distribution-free Hoeffding radius ε_H = √(ln(1/α_tost)/(2·5400)) = 0.0173 (no variance term), lower = m_lo − ε_H, upper = m_hi + ε_H. **The single-interval gate passes iff lower ≥ 0.6401 and upper ≤ 0.7001**, with at least 300 replicates and, for profiles, at least 10 used bins. What is exact: the sweep, the monotonicity argument and the bounds conditional on H. What rests on an assumption: the 0.99 level of the box needs Gaussian marginals of the H mean (a mean of 36 000 blocks; t-exact under Gaussianity and checked on the actual replicate means by the bootstrap). The paired gate carries the finite-sample, assumption-free guarantee. **An estimator passes iff both gates pass; row V7 passes iff all three estimators pass.** Calibration of the two gates (same Monte Carlo, each gate asserted separately at its own margins with κ calibrated on its own true coverage): Gaussian profiles — paired, single and row false acceptance 0 at every margin (upper bounds 0.005), row power 1.000 / 0.958; Gaussian scalar — single-gate false acceptance 0.009 / 0.020 (upper bounds ≤ 0.025), row ≤ 0.004, paired gate at its exact value 0.0387 (the Monte Carlo a consistency check), row power 0.947, joint ≈ 0.90; Gamma shape 2 and zero-inflated shape 2 (two bins excluded): false acceptance ≤ 0.007 (upper bounds ≤ 0.021), but the single gate's power at the nominal is only 0.14 for these skewed surrogates — the cost of the distribution-free radius plus the box on a half-region of 0.03 — stated here as a known limitation: a correctly calibrated but strongly skewed tally can fail the single-interval gate, and such a failure is reported with the single-interval coverage itself; Gamma shape 0.5 independent (single-interval coverage 0.638 below the band) is rejected by the single gate as intended; the correlated shape-0.5 surrogates have no resolvable bins and fail closed. *Finite-sample bounds, α_tost = 0.04 per side, exact for i.i.d. bounded variables with no other assumption:* scalar (f_p ∈ {0, 1}): one-sided Clopper–Pearson binomial bounds on the hit count of 3600 pairs; profiles (f_p ∈ [0, 1]): the empirical Bernstein bound ε = √(2V ln(2/α_tost)/P) + 7 ln(2/α_tost)/(3(P − 1)) (Maurer & Pontil 2009) with V the sample variance of the f_p; lower/upper = CP bounds or m ∓ ε. *Rule (step `v7-rep-combine`, consumes the nine partials, fail closed on any missing or invalid one or unequal shard replicate counts):* **the paired gate passes iff lower ≥ 0.6464 and upper ≤ 0.7064**, with at least 300 pairs and, for profiles, at least 10 used bins; an estimator passes iff the paired gate and the single-interval gate (below) both pass; **row V7 passes iff all three estimators pass.** A true paired coverage outside the region is accepted with probability at most 0.04 (for the scalar exactly 0.0386 / 0.0387 at the two margins, computed from the Clopper–Pearson acceptance set 2378 ≤ hits ≤ 2494), without any distributional assumption. The former pooled point gate 0.68 ± 0.03 is reported only. *Calibration as exact-SHA local evidence (gate-specific, current implementation):* step `pytest-v7-rep-calibration` runs (`tests/ionmc/test_v7_coverage.py`, marker `calibration`, excluded from the lightweight CI) the exact scalar false-acceptance computation of the paired gate (0.0386 / 0.0387 at its two margins, from the Clopper–Pearson acceptance set 2378 ≤ hits ≤ 2494 of 3600; this exact value is the scalar's paired-gate calibration and the Monte Carlo below is a consistency check at four standard errors) and the Monte Carlo of the two-gate procedure exactly as implemented (3600 simulated pairs; 5400 evaluation and 1800 held-out replicates; a simulated 1e6 reference for the bin mask; the degenerate-miss rule; 500 bootstrap resamples inside the Monte Carlo), each gate asserted on its own at its own margins with κ calibrated on its own true coverage (paired coverage at c₀ ± 0.03; single-interval coverage with sem_j alone at c₁ ± 0.03), false acceptance as the 95 % Clopper–Pearson upper bound ≤ 0.05. Results: Gaussian 12-bin profiles (between-bin correlation 0 and 0.8) — paired, single and row false acceptance 0 at every margin (upper bounds 0.005); power at the nominal: paired gate 1.000 / 0.983, single gate 1.000 / 0.972, row 1.000 / 0.958. Gaussian scalar — single-gate false acceptance 0.009 / 0.020 (upper bounds 0.012 / 0.024), row false acceptance ≤ 0.004 at all four margins; power: paired gate 0.958, single gate 0.988, row 0.947; joint row pass over the three estimators ≈ 0.90. Gamma-distributed blocks (sparse-tally surrogates): shape 2 independent — false acceptance 0 at all margins, paired-gate power 1.0, single-gate and row power 0.14; shape 2 zero-inflated (two bins of 90 % zeros excluded by the criterion, 10 bins used) — false acceptance ≤ 0.007 (upper bounds ≤ 0.021), paired-gate power 1.0, single-gate and row power 0.14; shape 0.5 independent — single-interval coverage 0.638 (below its band) while the paired coverage 0.660 is inside: the single gate rejects, row power 0 as intended; shape 2 correlated — paired 0.644 and single 0.604 at κ = 1, both gates reject; the correlated shape-0.5 surrogates have no resolvable bins and fail closed. The low single-gate power on strongly skewed surrogates is the stated cost of the distribution-free radius plus the reference box on a half-region of 0.03: a correctly calibrated but strongly skewed tally can fail the single-interval gate; such a failure is reported together with the single-interval coverage itself and the per-bin skewness, and is a V7 failure under this amendment. The lv5b suite has 25 steps (14–21 shards, 22 reference for the bin set and the cross-check, 23 calibration, 24 combine, 25 the deferred partition). The seeds k = 0..8 are unused; the superseded chains at 1790dc99, 15a4c640, 2554be89, df641ded, 273180af, 2f0d2f47, e31e15d7, 2c53b4e7, 353d14f5 and c5807962 executed steps 01–13 at most (one 0.8 s guard failure of step 14 at 353d14f5 before any transport; no V7 result) and are cross-checks, not evidence. A failed re-run is row V7's failure without a further base (Amendment 11).

**Amendment 14 (2026-10-10, revised the same day after Codex REVIEW-7d749e6a, before any V3-005C result; written and executed on the V3-005B task branch) — V3-005C scope: hadronic elastic scattering (D2), the V5 re-run with a pre-declared attribution ladder, first execution of V6, a pre-registered V7 re-test of `escaped_neutral`, and a conclusive python `nuc_local` parity sample.** Nothing in this amendment changes a frozen row, a frozen tolerance or the text of Amendments 1–13; it adds rows, execution details and seeds. Decision 0041, section "Slice C", records the physics choices and the rejected alternatives.

*Revisions (2026-10-10, after Codex REVIEW-7d749e6ae7d643258ced384524f8e20f at 20c234f; orchestrator decisions; nothing had been observed):*
1. V6: the z ≤ 0.95·R0 exclusion at r = 0 is withdrawn. Every r = 0 point gates at the frozen 2 %, and the distal-edge sensitivity becomes an informative annotation reported with the verdict. E₀ is set once from R0 through the range table, before any dose comparison ((g)).
2. V6 lineage: V6 is classified as evaluation with shared beam lineage (calibration-type inputs). The r ≥ 1 cm halo gates only under a V6-sens condition fixed now ((g), row V6-sens).
3. V10 becomes report-only, because S(E) shares lineage with PDG/BGG. The evidence gap for p-p above 150 MeV is recorded ((c)3, row V10).
4. The V10-A data become related-model evidence (shared/related lineage). V10-A stays gating as a model-consistency test, and the plan states that no fully independent p+O elastic angular data exist ((c)1, (c)4, row V10-A, (e)1).
5. V7-R: passing the calibration step is part of the pass condition ((h)).
6. The "attributed" wording for the V5 failure is corrected to "hypothesis" in the requirement ledger and the results README, and by dated correction notes in the archive table above and in the decision 0041 slice-B outcome.

*(a) Operator conditions (binding).* The operator resolved IR-20261009-032741-C590CF on 2026-10-10T13:48:52Z. The response is recorded verbatim in `reference_requests/interventions/IR-20261009-032741-C590CF/response.json`: "Regarding IR-20261009-032741-C590CF: Keep PR #82 unmerged under interpretation B. Continue V3-005C on top of the V3-005B task branch. Address the dose discrepancies, investigate the escaped-neutral uncertainty calibration, and obtain sufficient statistics for the inconclusive Python/Warp comparison. Preserve all existing failed and inconclusive records. Pre-register follow-up validation before generating new results; do not relax acceptance thresholds or repeat seeds merely to obtain a pass. Missing elastic scattering is a hypothesis to test, not an assumed complete explanation. Merge only after the applicable scientific and implementation gates pass. Do not create an implementation-only passing record to bypass the failed scientific gate." The consequences are as follows.
1. PR #82 stays unmerged. All V3-005C work is committed on `task/v3-005b-…`, and every task tool (commit, review, local validation, host runner, reference calculation, data acquisition) runs under task id **V3-005B**. "V3-005C" names the scope only.
2. The failed and inconclusive records stay as archived, and no later result re-labels them. These are: `lv5b-b40d812` (consumed); lv5b-7aae5bb V5 fail and V7 fail (the `escaped_neutral` paired gate); hr5-7aae5bb python:cpu64 `nuc_local` inconclusive. A pass of a new row below is reported next to them, never in place of them.
3. Every new tolerance below is fixed now, from the model-error and data-uncertainty estimates stated with it, before any model output is compared with the data it is judged against. The provisional-tolerance mechanism of the planning draft ("freeze V10-A in a later Amendment 14b, inflate by the data normalisation") is rejected, because it would choose tolerances after the data were in hand.
4. No base is re-used. A failed or inconclusive base is consumed (Amendment 6). Every fresh base is declared before it is observed ((j)).
5. The hadronic-elastic attribution of the V5 failure is tested by V11, the V5 re-run and the ladder in (f). It is not assumed.
6. The task is merged only if the scientific gates (V5, V6 with the status fixed by V6-sens, V7-R including its calibration step, V8-nuc, V10-A) and the implementation gates (P6, P7, P5-ext, C1-ext, R1-D2, V8-LV, CI) pass. A summary that passes only the implementation rows is not a passing record of this task, and none is created.

*(b) Physics under test (decision 0041 "Slice C"; orchestrator decision D-1).*
- **p+A (C, N, O, Ca and the decision-0041 surrogates), 1–250 MeV:** σ_el(E) = BGG/Barashenkov max(σ_tot(n+A) − σ_inel(p+A), 0).
  - The Barashenkov arrays are hand-transcribed from Geant4 11.4.2 `G4BarashenkovData.hh`. The rules are linear interpolation in E and A^{2/3} interpolation for untabulated Z, with the BGG low-energy rule below 14 MeV. The formulas are hand-re-implemented from `G4BGGNucleonElasticXS` / `G4NucleonNuclearCrossSection::GetElasticCrossSection`, and the Geant4 Software License notice is recorded.
  - Shared lineage: TOPAS uses the same σ_el.
  - Angular shape: the nuclear-only black-disk form dσ/dΩ_CM ∝ |2J₁(qR)/(qR)|², q = 2p_CM sin(θ_CM/2). R comes from π(R + ƛ)² = σ_nonel(E) of the transport's own table (LA150 MT5 ≤ 150 MeV, Tripathi-shape extension above), so no elastic data enter the shape.
  - The form is sampled over all angles from an inverse-CDF in μ_CM per energy node. Coulomb scattering stays entirely in the MCS model, and no Coulomb–nuclear interference is modelled. The recoil nucleus is deposited locally (`nuclear_local`).
- **p-p (H-1):**
  - ≤ 150 MeV: the LA150 Hale LAW=5 LTP=1 (LIDP=1) reconstruction (ENDF-102 eqs 6.9/6.10/6.14, b/sr CM, factor 2π on integration). The tables interpolate the ratio σ_e/σ_c and never P_NI. The transported density is NI = σ_e − σ_c over θ_CM ≥ 16.26° (|μ_CM| ≤ 0.96, the NJOY `umin` convention).
  - Above 150 MeV: σ(E) = σ_NI(150)·S(E)/S(150), with S = 1.0115·(23 + 50·√(ln(0.73/p)^7)) mb (BGG/PDG systematics, p_lab in GeV/c). The shape is held fixed in μ_CM.
  - Kinematics are exact relativistic two-body. The faster proton continues as the primary, and the slower one is transported as a secondary proton.
  - Hale's π·b₀ is 14–16 % above the PDG σ_el(pp) at 100–150 MeV. This excess is recorded, not corrected.
- **Sampling:** Σ_tot = Σ_nonel + Σ_el in the existing thinning (one candidate per step, fail-closed `majorant_violation` unchanged). The channel is chosen by Σ_el/Σ_tot. The EM random-number streams are identical with nuclear on and off. The particle bound of 32 and the CUDA chunk ≤ 2^14 are unchanged, and elastic adds at most one secondary per event.

*(c) Lineage and data roles, per channel.* Decision 0041 sets the rule: measured data published before 1997 share lineage with LA150 and are exploratory or report-only; data from 1997 on are evaluation data. The rule is applied per channel as follows.
1. **Elastic angular shape (p+A).** No elastic data enter its construction directly. The form is parameter-free, and R comes from the transport's own σ_nonel. The pre-1997 p+16O and p+12C angular distributions (list in (e)) are nevertheless **related-model evidence (shared/related lineage)**, for two reasons:
   - LA150's σ_nonel (optical-model/GNASH evaluation) may have been adjusted to the same elastic measurements. Through R it sets the angular scale qR, which is the whole shape of a parameter-free form, not one minor scalar.
   - The σ_el level that multiplies the shape comes from Barashenkov's compilation (1989/1993), which may contain the same measurements.

   V10-A therefore stays **gating as a model-consistency test** under this label. **No fully independent p+O elastic angular data exist**: every p+16O angular dataset found in EXFOR for 75–250 MeV predates 1997.
2. **Elastic level σ_el(p+A).** This level is gated by post-1997 data only, through the proxy σ_tot(n+A) − σ_nonel(p+A). For C the proxy is Abfalterer 2001 (EXFOR 13753-007, n+C σ_tot, 5.3–559 MeV) minus Auce 2005 (D0356-002, p+C σ_nonel, 81–≈180 MeV) and Ingemarsson 1999 (O0579-003, p+C σ_R at 65.5 MeV, 295.5 ± 7.7 mb).
   - For O no post-1997 n+O total cross section exists in EXFOR (13753 has no O subentry). The O level is therefore **report-only**: Ingemarsson 1999 (O0579-004, p+16O σ_R at 65.5 MeV, 365 ± 15 mb) minus Finlay 1993 (13569-005, n+O σ_tot, pre-1997).
   - This is an **evidence gap**, stated in decision 0041: O carries about 89 % of the nuclear-elastic mass in water.
3. **p-p.**
   - Hale's R-matrix fit (FEB98, 0–150 MeV) lists no data, so every p-p measurement at ≤ 150 MeV published before February 1998 is treated as in the fit.
   - V10 compares with the PDG compilation **above** 150 MeV, outside Hale's range. It is **report-only**: the extrapolation rule S(E) is itself PDG/BGG-derived (shared lineage), and the 160 MeV point was seen before registration. **Evidence gap:** no p-p elastic data independent of the Hale fit and of the PDG compilation were identified above 150 MeV.
   - The p-p dσ/dΩ entries at 140–260 MeV are report-only: Mahjour-Shafiei 2004 O1226-002 (190 MeV, post-1997); Avan 1984 O0145-007; Chamberlain 1954 C2637-002 and 1951 C2606-003; Cassels 1956 O2057-002.
   - The p-p channel is independent of TOPAS, which uses BGG pp.
4. **Consequence for the reference engines.** σ_el(p+A) is the same Barashenkov quantity in TOPAS. On the p+A elastic normalisation, V5 and V11 are therefore **not independent**: they test transport, angular shape and implementation. No fully independent evidence for the elastic channel exists in this scope. What remains is:
   - V10-A, gating as a model-consistency test with related lineage; its C level is checked against post-1997 data.
   - V10 (p-p), report-only.
   - The V6 halo at r ≥ 1 cm, which has shared beam lineage and gates only under the V6-sens condition of (g).

*(d) Rows frozen now.* These rows enter the plan-row accounting test of `tests/ionmc/test_run_suite.py` ("a new plan row needs a step or a declaration") through the explicit allowlist `PLAN_ROWS_PENDING_STEPS` until their lv5c/hr5c steps exist ((m)).

| # | Check | Pass iff | Class / tier | Justification |
|---|---|---|---|---|
| P6 | Data layer, deterministic: (1) MF6 LAW=5 parser for LTP=1 (LIDP 0/1; NW = 4NL+3 / 3NL+3) and LTP=12 (NW = 2NL) on hand-written fixtures; (2) eqs 6.9/6.10/6.14 in b/sr CM, identical-particle (spin ½) interference coefficient −1/2, phase sign exp(+iη ln …) as NJOY `coul`, integrated cross sections × 2π; (3) LA150 reconstruction values: O-16 100 MeV NI(θ_CM ≥ 5°) 123.8 mb and σ_e(θ_CM ≥ 5°) 185.6 mb; H-1 NI(θ_CM ≥ 16.26°) 49.8 / 27.3 / 26.3 mb at 50 / 100 / 150 MeV; (4) ratio interpolation, and **failure on any negative transported density at any node or node midpoint**; (5) BGG: pp 27.61 mb (100 MeV) and 23.86 mb (150 MeV) from the formula with p_lab 0.4446 / 0.5513 GeV/c; p+O 352 / 159 / 102 / 92 mb at 100 / 150 / 200 / 250 MeV; p+C 253 / 125, N 282 / 132, Ca 912 / 400 mb at 100 / 150 MeV; (6) shape normalisation ∫dσ/dΩ dΩ = σ_el and the inversion residual of R(σ_nonel); (7) two-body kinematics closure per event | Integers exact. Parser floats bitwise. Formula re-evaluations 1e-12 relative. LA150 reconstruction values within 1e-3 relative. Research-quoted BGG values within half a unit of the last quoted digit. Kinematics \|ΔE\| ≤ 1e-9 MeV and \|Δp\| ≤ 1e-9 MeV/c per event. Zero negative densities. | self-consistency, CI (data cases skip in CI and fail under `IONMC_REQUIRE_DATA=1`) | The tolerances are set by float precision and by the precision to which the reference numbers were recorded. The round-1 values come from an independent numpy reconstruction (scratch evidence, not repository code). The two quoted pp values at 150 MeV (23.9 in round 1, 23.8 in round 2) are one number, 23.857 mb, recomputed here from the formula. |
| P7 | Sampler vs the built tables, offline: H-1 at 20 / 50 / 100 / 150 / 200 MeV; C-12 and O-16 at 20 / 50 / 100 / 150 / 200 / 250 MeV; 1e5 events per case; μ_CM histogram in ≥ 50 equiprobable bins of the table CDF; σ_el table vs direct builder evaluation at every node and node midpoint | χ² p > 0.001 in every case; σ relative error ≤ 1e-3 at every point | self-consistency, LV, r_index 15 | As P2/N1. With 17 cases, a correct sampler fails at least one case with probability 1.7 %. |
| P5-ext | Every new `@wp.func` twin vs Warp-CPU float64 on recorded inputs (U1 harness). Event-level parity of complete elastic events on recorded inputs: channel choice, μ_CM, outgoing direction and energy, recoil energy, and the recoil proton's state for p-p | Discrete outputs equal; continuous outputs within 1e-12 relative | backend parity, CI | As P5. |
| V10 | p-p level above Hale's range: the transported σ_pp(E) = σ_NI,θ≥16.26°(150 MeV)·S(E)/S(150) at **every** PDG rpp2022 `pp_elastic` point with 0.556 < p_lab ≤ 0.729 GeV/c (T = 152.4–249.9 MeV) | **Report only.** The per-point ratios and the flag \|model/data − 1\| ≤ 0.08 are recorded; there is no verdict | measured, report-only (shared lineage: S(E) is PDG/BGG-derived), LV, deterministic | PDG point errors are 3–5 % plus about 3 % normalisation. On the model side, S(E) changes σ by ≤ 4 % from 150 to 250 MeV, BGG is 0–5 % below the PDG points at 94–160 MeV, and the excluded forward cone (< 16.26° CM) is part of the model error. **Disclosure:** the 160 MeV point (p 0.57 GeV/c, 26.1 ± 1.0 mb) was seen in research (model 26.3 mb at 150 MeV), so the row was not blind at registration. No p-p elastic data independent of the Hale fit and of the PDG compilation were identified above 150 MeV (evidence gap). Points at 94–145 MeV are report-only (possibly in Hale's fit). **V10-pp-b** (report-only): dσ/dΩ_CM(90°) at 50 / 100 / 150 MeV vs SAID SP07 (shared data lineage with Hale); dropped and recorded if SAID is unreachable. **V10-pp-c** (report-only): p-p dσ/dΩ shape at 140–260 MeV vs the EXFOR entries of (c)3. |
| V10-A | p+A nuclear elastic vs measured data (entries and execution in (e)): (a) dσ/dΩ_CM at θ_CM = 10°, 15°, 20° for every gating angular dataset; (b) the partial integral from 10° to the data's first minimum; (c) the σ_el level vs σ_tot(n+C) − σ_nonel(p+C) at the post-1997 energies | (a) every evaluated point has \|model/data − 1\| ≤ **0.30**; (b) every evaluated dataset has \|model/data − 1\| ≤ **0.20**; (c) every gating point has \|model/proxy − 1\| ≤ **0.12**; plus the minimum coverage of (e). Otherwise the row fails, or it is "not executed" when the coverage is missing | related-model evidence (shared/related lineage), gating as a model-consistency test, LV, deterministic on the built table; r_index 18 seeds a sampled-histogram cross-check (report-only) | (a) A diffraction form vs data before the first minimum is typically within 20–30 % (round-2 model-error estimate for this option), and the data normalisation is 3–5 %. The tolerance sits at the upper end of that range, so with about 20–30 points a correct model may still fail one; this is accepted and stated. (b) Integration averages the per-angle scatter: ⟨θ²⟩ within about 15 %, σ_el ±10 %, normalisation 3–5 %. (c) The Barashenkov table quality is ±10 %. The proxy is a difference that is 35–55 % of σ_tot at 65–180 MeV, so σ_tot 1–2 % and σ_R 3–5 % propagate to 5–8 %. The quadrature sum is 11–13 %, so a single-point fail is weak evidence of a > 12 % model error; this is stated, not adjusted. |
| X-ENDF | Model nuclear dσ/dΩ_CM at 20–40° vs LA150 NI (ratio interpolation) for C-12, N-14, Ca-40 at 50 / 100 / 150 MeV | Report only | related evaluation, LV, r_index 19 (sampled cross-check) | Not like-for-like, because CNI is present in ENDF and absent in the model. O-16 is excluded: its MT2 is a copy of C-12 (decision 0041 data finding). |
| V11 | Elastic-only effect, ionmc vs TOPAS, at 150 and 200 MeV. ionmc warp-cpu float64: EM+elastic (internal switch `elastic_only`, non-elastic off) and EM-only, with common random numbers, 1e5 histories each, 20 batches, the V5 geometry (r = 20 cm). TOPAS: EM+elastic (rows 10/11 of (f), 3 seeds) and the frozen V5 EM-only cases. F = (peak/plateau)_EM+el / (peak/plateau)_EM | \|F_ionmc/F_TOPAS − 1\| within 0.02, and the plateau gain (plateau_EM+el/plateau_EM − 1) within 0.005 absolute. Both use the V5 TOST on the 90 % interval with Welch df (Amendment 9) | independent MC (shape and implementation only; σ_el shared), LV, r_index 14 | 2 % equals V5's peak/plateau tolerance, so an elastic-model difference that V11 accepts cannot by itself fail V5. The expected plateau gain is +1.2 / +1.4 % (estimate in (f)), so 0.5 % absolute is about 40 % of the effect. The TOPAS X2 single seed gave F = 0.894 at 150 MeV. |
| C1-ext | Fail-closed cases: elastic table missing, stale (re-hash mismatch) or unqualified; an element without elastic data; E above the 250 MeV domain; a forced Σ_tot majorant violation; a negative built density (the builder refuses); missing provenance, licence or lineage metadata; p-p kinematics closure failure; an LA150 O-16 MT2 offered as a construction input | Raises `UnsupportedCombinationError` (or a build error) before transport, or `valid=False` with the named counter > 0 | capability, CI | As C1. |
| R1-D2 | Intended-change record for the elastic change, A16-style. (A) `nuclear=False`: the A16 digest is bit-identical to f3a1dd62 and the intended-change set is empty (row R1 unchanged). (B) `nuclear=True` with Σ_el ≡ 0 (internal switch): bit-identical to the V3-005B nuclear-on outputs (python primary-only digest `3edfce9c…` = b84fdf38, and the C13 secondary-proton record), provided the channel choice takes the free candidate-block slot (purpose 2). If C4 shows that the slot is not free, (B) is replaced, in that commit and before any lv5c result, by an intended-change record of the RNG keying plus V8-LV parity. (C) The record for nuclear-on outputs, with directions declared now: peak/plateau down, plateau up, lateral spread up, \|Δ total deposit\| ≤ 0.1 % of E₀ | (A, B) bit-identical; (C) record present, reviewed and matching every declared direction | regression, CI/LV, r_index 3 | (B) shows that elastic touches nothing else. The directions in (C) follow from the estimate in (f): Δpeak −11 / −15 %, Δplateau +1.2 / +1.4 %. Elastic conserves energy except the locally deposited recoil. |
| E1-B | Exploratory python IDD at 150 MeV (1e5 histories) vs TOPAS and MCsquare | Report only | –, LV, r_index 13 | As E1 (slice A). |
| D9-BGG | The acquired Barashenkov σ_inel(p+A) and σ_el(p+A) arrays vs the post-1997 EXFOR set, windows of V1b / Amendment 7(i) | Report only | related model (TOPAS lineage), LV, deterministic | Quantifies the shared lineage of a construction input of (b). |
| V6-sens | V6 beam sensitivity with common random numbers, N_V6 histories per variant: σ_x, σ_θ and the x–θ correlation each varied by ±20 % (six variants), plus R0 ±0.5 mm through E₀ (report-only) | **Status condition for V6, fixed now.** The r ≥ 1 cm halo part of V6 gates only if, for every gating r ≥ 1 cm point and each of the six beam variants, \|D_var/D_nom − 1\| + 2σ_Δ < 0.15 (the per-point tolerance). Otherwise V6 as a whole becomes report-only, and this is recorded. The R0 variants are reported only | –, LV, r_index 20 | The beam parameters p1 and p5–p7 come from the paper's 25-parameter fit to the same `tbl:dmlg` data, and A0/A1 are fit artefacts of a 2-D model. The halo counts as a physics test only if it is insensitive to these calibration-type inputs at the scale of its tolerance. |
| V7-R | Re-test of the V7 coverage with the bootstrap-t interval for `escaped_neutral` (h) | (h) | falsification, LV, r_index 16 (diagnostic 17) | (h) |
| V2b-el | Frozen row V2b re-executed with the elastic channel on (warp-cpu float64, s_max 0.1 vs 1.0 mm, 1e6 per variant, common random numbers) | As V2b (frozen, unchanged) | falsification, LV, r_index 8 | Elastic adds a thinned channel to the step machinery that V2b tests; the pre-elastic lv5b pass stays archived. |
| V7-scan-el | Frozen row V7, N-scan and f32/f64 on warp-cpu, re-executed with the elastic channel on (Amendment 8 design: 100 batches per N) | As V7 N-scan (relative-SE slope −0.5 ± 0.05; f32/f64 as lv5b step 12) | falsification, LV, r_index 9 | Elastic changes every secondary tally; the pre-elastic lv5b pass stays archived. |
| V7-shift-el | Frozen row V7, grid shift and refinement, re-executed with the elastic channel on | As Amendment 8 (≤ 1e-3 and ≤ 1e-6 relative) | falsification, LV, r_index 10 | As V7-scan-el. |
| V8-nuc | python:cpu64 `nuc_local` T12 parity with the (i) sample, as part of the full hr5c V8 statistical parity | T12 conclusive and passing under its unchanged rule. An inconclusive verdict is not a pass, and no larger sample is taken afterwards | parity, HR, hr5c r_index 1 | (i) |

The frozen rows V5, V6, V8 and R1 are evaluated again with their tolerances unchanged. V2b and the V7 N-scan, grid shift and refinement are re-executed with elastic on as V2b-el, V7-scan-el and V7-shift-el, with the frozen tolerances; the pre-elastic lv5b results at 7aae5bb stay archived.

*(e) V10-A execution.*
1. **Gating angular datasets (role "related-model evidence (shared/related lineage)", (c)1; no fully independent p+O elastic angular data exist).**
   - p+16O: Kelly 1989 C0063-002 (135 MeV); Kelly 1990 C0141-002 (179.9 MeV); Seifert 1993 C0148-002 (201.4 MeV); Glover 1985 C0550-002 (200 MeV).
   - p+12C: Strauch & Titus 1956 O0274-009 (96 MeV); Rolland 1966 O0553-008 / O0553-002 (75 / 150 MeV); Meyer 1983 C0055-002 (every energy of its EN column in 120–200 MeV); Meyer 1981 C0057-002 (200 MeV); Meyer 1988 C1420-002 (250 MeV).
   - EXFOR has no p+16O entry for Rolland 1966 (O0553 holds C, Ca, Bi and Ce) and no entry for Gerstein 1957, so neither is used.
2. **Gating level points ((c)2).**
   - Each Auce 2005 D0356-002 energy, and Ingemarsson 1999 O0579-003 at 65.5 MeV. Abfalterer 2001 13753-007 σ_tot(n+C) is interpolated linearly in E between its two neighbouring points, and its errors are propagated.
   - O: Ingemarsson O0579-004 minus Finlay 13569-005, report-only.
3. **Frames.**
   - Laboratory-frame data are converted to CM with relativistic two-body elastic kinematics (AME2020 target masses) and the Jacobian.
   - Ratio-to-Rutherford data are converted with the Rutherford formula of item 4.
   - Units follow the P1 EXFOR parser.
4. **Model value.** Nuclear Bessel dσ/dΩ_CM (σ_el × normalised shape) plus the unscreened relativistic point-Rutherford dσ/dΩ_CM, added incoherently.
   - This is the transport's content: Coulomb in the MCS, nuclear in the elastic channel.
   - The Rutherford share is reported per point. The estimate for p+O at 135 MeV is about 10 % at 10°, ≤ 4 % at 15° and ≤ 3 % at 20°.
   - CNI is not modelled. Its effect is part of the model error that the 30 % covers.
5. **Measured value at an angle.** Log-linear interpolation of dσ/dΩ in θ_CM between the two data points that bracket the angle, both within 3° of it. Otherwise the point is not evaluated, and this is recorded.
   - θ_min,1 is the first measured point above 10° that is lower than both of its neighbours.
   - Points with θ ≥ 0.8·θ_min,1 are not evaluated in (a), because the model error near and beyond a minimum is a factor ≥ 2. They are reported. If a dataset has no measured minimum, the model's first minimum sets this cut, and (b) is not evaluated for that dataset.
6. **Partial integral (b).** The trapezoid in θ of 2π sin θ dσ/dΩ over the measured points in [10°, θ_min,1], with interpolated endpoints. The model is integrated on the same nodes with the same rule, so the quadrature error cancels.
7. **Minimum coverage.** ≥ 12 evaluated (a) points; at least one p+16O dataset at 130–140 MeV and one at 175–205 MeV; at least one p+12C dataset; ≥ 3 (b) datasets; ≥ 4 (c) points.
   - If the coverage is not met, V10-A is reported as **not executed (evidence missing)**. D2 is then not claimed as validated for p+A, and the summary stays `conformant: false` with a deferred row. This is never a pass.
8. **Uncertainties.** The data's statistical errors are reported and are not added to the frozen tolerances. The r_index-18 cross-check (1e6 sampled events per dataset energy, histogrammed at the data angles) is report-only.
9. **Fallback (ii′) trigger (orchestrator D-3).** The ENDF/B-VIII.0 **neutron** O-16 / C-12 MF4/MT2 elastic shape replaces the Bessel shape only if all of the following hold:
   - V10-A (a) or (b) fails for p+O;
   - V10-A (c) passes;
   - P7 passes.

   Never because of V5 or V11 alone. In that case:
   - The (ii) failure is a recorded scientific failure.
   - A new decision is taken.
   - The neutron file is acquired, and its MF1 is checked so that n-O16 above 20 MeV is not itself a C-12 copy.
   - All V10-A results for (ii′) are **exploratory**, because the data have been seen. Gating for the shape then rests on V6 and V11.
   - The re-run takes the next declared base ((j)).

*(f) V5 re-run, expectation, attribution ladder and reference requests.*
- **Re-run.** V5 is re-run at lv5c r_index 4–7 with the frozen tolerances, frozen reference runs and frozen comparator. TOPAS gates and MCsquare is report-only (Amendments 7(b), 9, 10).
- **Pre-registered expectation** (planning estimates, not criteria).
  - At 150 MeV: ionmc/TOPAS peak/plateau 0.985 (range 0.95–1.03); ΔIDD about −2.5 % (−8 … +3 %).
  - At 200 MeV: 0.986 at δ = 5 mm (range 0.90–0.99).
  - The basis is a straight-line CSDA estimate with Barashenkov σ_el and P(p-O) 0.21 / 0.24 per primary, valid to about ±4 % on the factor. The TOPAS X-case route gives 0.998 / 0.996.
  - Plateau, R80 and total deposit stay passing, and |Δ total deposit| ≤ 0.1 %.
- **Attribution ladder.** It applies only if V5 fails. An attribution is recorded as an attribution, and V5 stays failed.
  - Rung 0 (implementation): P6, P7, P5-ext, R1-D2 and V8-LV must pass. If one fails, the defect is fixed, the fix is reviewed, and one fresh base follows.
  - Rung 1 (elastic-model difference): V11 at the failing energy. If V11 fails in the direction of the V5 residual, the failure is attributed to the p+A angular shape (radius rule of D-1). It is a recorded scientific failure, and (ii′) is entered only through its trigger in (e).
  - Rung 2 (V11 passes, ionmc above TOPAS by 2–5 % on peak/plateau): the inclusive-yield energy excess (`nuclear_imbalance`, −0.7 % of E₀) and α local deposition (0.6 % of E₀).
    - Test: the non-elastic factor F_ne = (peak/plateau)_nuclear on / (peak/plateau)_elastic only, ionmc (V5 "on" vs V11) against TOPAS (full vs rows 10/11).
    - Attributed iff the 90 % interval of F_ne,ionmc/F_ne,TOPAS − 1 lies beyond 0.02 on the side of the residual. The depth-resolved difference in the last 10 mm before R80 is reported.
  - Rung 3 (nuclear-recoil local deposition): considered only if rungs 1–2 exclude.
    - Attributed iff the recorded recoil energy deposited in [peak − 10 mm, R80] is ≥ 50 % of the residual's energy equivalent.
    - Otherwise the failure is recorded as **unexplained**.
  - **Stop rule:** if peak/plateau fails **low** (ionmc below TOPAS by more than 2 %) or by more than 5 % in either direction, the task stops. This is a recorded scientific failure, escalated to the orchestrator, with no tolerance change.
  - Exactly one fresh base per reviewed fix. Each such fix is recorded as an intended change. No further base is taken without a new reviewed physics reason.
- **Reference-calculation requests.** These run inside the task through `run_reference_calculation`, with the same engine access as Amendment 9 (TOPAS, Geant4 11.4.2), committed case directories, the Amendment 9 geometry and scoring (400 × 400 mm, 0.5 mm bins, 1e5 primaries) and the fingerprint verification of Amendment 10(b).
  - Seeds are 20270000 + 1000·row + k, k = 1..3. Row 9 is consumed by the exploratory X1–X4 cases (20279001–20279004), so the new rows start at 10:
    - row 10: EM+elastic at 150 MeV (`g4em-standard_opt4` + `g4h-elastic_HP`, as X2);
    - row 11: EM+elastic at 200 MeV;
    - row 12: full minus hadronic elastic at 200 MeV (as X1);
    - row 13: full minus hadronic elastic at 150 MeV (three-seed replacement of the single-seed X1).
  - Rows 10/11 serve V11 and rung 2. Rows 12/13 give the full/X1 elastic factor, report-only.

*(g) V6 execution (frozen row V6 unchanged).*
- **Data.** Gottschalk 2015, arXiv:1409.1938 source (`halo.tex`), acquired through `acquire_reference_data` and pinned by sha256 (scratch copy `e5a270e7…23d6`; a mismatch is recorded).
  - Role: **evaluation (measured)** for `tbl:dmlg`: 46 depths z = 1.533–25.000 cm × 10 radii (0–8 and 10 cm), log₁₀(MeV/g/p). 'M' marks a missing point, which is excluded and counted.
  - Role: **construction for the beam only** for `tbl:parms` (p1, p2, p5–p7), never for the dose comparison.
  - The arXiv licence field is recorded at acquisition. A difference from the PMB 60 (2015) 5627 table, if found, is recorded, and the arXiv table stays the pinned one.
- **Lineage and classification.** V6 is **evaluation with shared beam lineage (calibration-type inputs)**: the beam parameters p1 (R0) and p5–p7 come from the paper's 25-parameter MD fit to the same `tbl:dmlg` data.
  - The r = 0 depth dose and the beam-size-sensitive core are shared-lineage evidence.
  - The r ≥ 1 cm halo is the physics test. It gates only under the V6-sens condition (row V6-sens): if V6-sens fails that condition, V6 as a whole becomes report-only, and this is recorded.
  - The source contains no independent entrance-beam measurement (σ_x ≈ 0.6 cm "as determined by fits", l.147), and we did not obtain one.
- **Beam at the water entrance.** The MD fit transports the moments in water only (l.728–732). Its A-moments are defined at the tank entrance and absorb the 0.008 cm Kapton window, the 150 cm He bag, the 30 cm air gap and the monitor. These are therefore **not** simulated: simulating them would double count.
  - **E₀ procedure (pre-declared; computed once, before any dose comparison; ionmc never fits the dose data).** R0 = p1 = 21.321 cm is Gottschalk's d80, which for a pristine beam equals the mean projected range.
    - E₀ solves R_CSDA(E₀) = R0 in the project's water range table (the transport's stopping power, I = 78 eV, the V3-003D log-log quadrature), by bisection to 1e-6 relative, with no detour-factor correction.
    - It is computed deterministically in C8 and frozen in Amendment 14a with its value. The expected value is 177 ± 0.5 MeV; a value outside 176–178 MeV is recorded as a range-model finding and still used.
    - The energy spread below uses dR/dE from the same table.
    - ionmc's own d80 at that E₀ is reported next to R0 and is never used to adjust E₀.
  - Gaussian energy spread: the beam contributes 0.8 % of R0 in range (σ_R,beam = 0.1706 cm, converted with dR/dE at E₀), on top of the model's own straggling. Gottschalk's total is σ_sem/R0 = p2 = 1.372 %, with Janni straggling 1.11 %.
  - Per transverse plane, a bivariate Gaussian with σ_x = p5 = 5.948 mm and σ_θ = p6 = 31.55 mrad. The fit parameter p7 = −0.1462 gives f = −(1 + p7) = −0.8538 (halo.tex l.776–780) and A1 = −√(1 − f²)·√(A0A2) = −9.770e-3 cm·rad.
  - The **x–θ correlation coefficient is −0.521** (converging; not −0.1462, which is p7). The emittance is π·160.2 mm·mrad.
- **Scoring.** Axisymmetric dose on a 0.5 mm (r, z) grid, in MeV/g/p, with absolute normalisation per primary (no fitted scale; the monitor uncertainty of ±3 % is reported, not added).
  - At each data (r, z) the dose is averaged over the Exradin T1 internal volume, modelled as a 0.4 cm diameter sphere (0.053 cm³ against a sphere's 0.034 cm³; the averaging scale, not the volume, is what matters).
- **Criteria execution.**
  - r = 0: **every** r = 0 point gates at the frozen 2 %, including the distal edge. A failure there is a recorded failure, not an excluded region.
  - *Informative annotation, reported next to the verdict:* σ_sem = 1.372 % of R0 = 0.293 cm. For a Gaussian distal edge (φ/σ relative to the local dose), the relative dose gradient is about 1.2 % per 0.1 mm at d80 and 2.8 % per 0.1 mm at d50. So a range difference of a few tenths of a millimetre (for example CSDA vs projected range) uses most or all of the 2 % tolerance there. At z ≤ 0.95·R0, a 0.1 mm shift changes the dose by ≲ 0.2 %. Each r = 0 residual is reported with its range-shift equivalent.
  - r ≥ 1 cm: 15 % per point, ≥ 90 % of points pass, median |residual| ≤ 10 %.
  - "No r-trend": for every depth with points in both groups, the sign of median residual(4–10 cm) − median residual(1–3 cm). Two-sided binomial sign test over the depths, p > 0.01.
- **History count rule.** N_V6 is the smallest whole number of 1e6-history shards that gives a rehearsal-estimated relative SE ≤ 3 % at every r ≥ 1 cm point whose measured dose is ≥ 1e-4 of the r = 0 maximum, and ≤ 0.5 % at every gating r = 0 point. It is estimated from a 2050xxxx rehearsal at the evidence code SHA, fixed in Amendment 14a and never decreased. The expectation is ≥ 1e7.
- **Backends.** V6 is evaluated on warp-cpu float64 (lv5c r_index 12, sharded so that each step takes ≤ 1500 s) and on warp-cuda float32 (hr5c r_index 3) at N_V6. Both must pass, and a disagreement between them is a V6 failure.

*(h) V7-R: re-test of the `escaped_neutral` uncertainty calibration (orchestrator D-4).*
- **Design.** As Amendment 13: 150 MeV, warp-cpu float64, 12-bin geometry, blocks of 500 histories, replicates of 20 blocks, 7200 replicates, an independent 1e6-history reference.
  - Change, justified now: the lv5b 9e6-history shards took 1308–1315 s without elastic, which leaves 13 % headroom against the 1500 s step limit, and elastic adds candidate events.
  - The shards therefore become **16 × 4.5e6 histories** (450 replicates each). Pairs are replicate j of shard s with replicate j of shard s + 8 (s = 0..7, 3600 pairs). The held-out set H is shards 12–15 (1800 replicates), and the evaluation set is shards 0–11 (5400). The reference is k = 16.
  - Every count of Amendment 13 is unchanged. Shard seeds are k = 0..15 of r_index 16.
- **Recording.** Each shard additionally writes an npz sidecar `escaped_neutral_block_sums[replicate, block]` (float64, 450 × 20). Its sha256, `batch_histories` and `blocks_per_replicate` go into the hash-bound JSON partial, and the per-replicate mean/sem are kept. Only `escaped_neutral` is recorded at block level, because only it failed (1.15 MB in total).
- **Estimators and gates (pass rule, stated before any observation).**
  - **`escaped_neutral`:** the bootstrap-t (percentile-t) interval, judged by the two-gate structure of Amendment 13. Per replicate there are B = 1999 resamples of its 20 block sums, t*_b = (m*_b − m̂)/se*_b, and the interval is [m̂ − t*_(1−α/2)·se, m̂ − t*_(α/2)·se] at α = 0.3173 (nominal 0.6827).
    - Bootstrap streams are seeded by `numpy.random.SeedSequence([20481004, 16, shard, replicate])`, and `[…, 1000 + s, j]` for pairs.
    - Paired gate: for pair (j, k), a two-sample bootstrap-t interval for μ_j − μ_k, with t* = ((m*_j − m*_k) − (m̂_j − m̂_k))/√(se*_j² + se*_k²) and B = 1999. A hit is a pair whose interval contains 0. A degenerate se (zero) is a miss. The bound is Clopper–Pearson at α_tost 0.04 per side.
    - Single-interval gate: the Amendment 13 held-out reference H with its t/bootstrap box and the exact sweep. A hit is a replicate whose bootstrap-t interval contains μ_H + δ·SEM_H. The bound is Clopper–Pearson on K_min / K_max at α_tost 0.04.
    - Both gates use the region **[0.6527, 0.7127]** (0.6827 ± 0.03, the margin of the frozen row). ≥ 300 pairs / replicates.
    - The 20-block Student-t interval of `escaped_neutral` is evaluated under the Amendment 13 rule and is **report-only**. It cannot turn the Amendment 13 failure into a pass.
  - **`sec_p` and `nuclear_local`:** the Amendment 13 two-gate rule, unchanged, at the new base. They are re-tested because elastic changes these tallies, and row V7-R needs all three estimators.
  - **V7-R passes iff the calibration step `pytest-v7r-calibration` passes at the evidence SHA (requirements below) and all three estimators pass both of their gates.** The gates alone never make V7-R pass. If the calibration step fails, V7-R fails with "uncertainty-coverage evidence not established".
- **Calibration** (exact SHA, step `pytest-v7r-calibration`, part of the V7-R pass condition). The requirements hold for each estimator and each surrogate:
  - false acceptance of the estimator's two gates jointly ≤ **0.02** at each margin (95 % Clopper–Pearson upper bound of the calibration Monte Carlo);
  - each gate alone ≤ 0.04 (its α_tost construction);
  - power ≥ **0.9** at the nominal.
  1. Surrogate from the archived 7200 lv5b replicates (7aae5bb partials, (mean, sem) pairs): the empirical studentised quantiles are 2.5 % / 97.5 % = −2.228 / +2.005 (t₁₉: ±2.093; two-sided rejection 5.18 %). The 15.87 / 84.13 % quantiles are computed deterministically from the same archive in C8. A parametric skewed block distribution (Gamma or zero-inflated, as in Amendment 13) is fitted so that its simulated 20-block t distribution reproduces those four quantiles.
  2. Surrogate from the new block sums: a rehearsal at 2050xxxx with ≥ 2 shards of 4.5e6 at the evidence code SHA. Replicates are drawn as 20 blocks from the pooled rehearsal blocks, with the pooled mean as the truth.
  3. Margins are reached by a scale κ on the interval half-widths, as in Amendment 13.
  4. If the rehearsal calibration misses a requirement, the replicate count or B may be increased in Amendment 14a before the evidence run, never after and never decreased. If the calibration step at the evidence SHA misses a requirement, **V7-R fails** ("uncertainty-coverage evidence not established"), whatever the gates give.
- **Diagnostic (exploratory, r_index 17).**
  - b ∈ {5, 10, 20} blocks per replicate, by exact aggregation of the recorded block sums.
  - b ∈ {40, 80} only if C6 shows that 250- and 125-history sub-block sums can be emitted **without changing the RNG keying** (otherwise they are dropped and the drop is recorded). r_index 17 is used only if such a separate sub-block run is needed.
  - The per-replicate excess kurtosis is reported. There is no pass or fail.
- **Outcome.**
  - Pass: V2-NUM coverage evidence for `escaped_neutral` with the bootstrap-t estimator. The reported uncertainty of that scalar switches to it, recorded as an intended change.
  - Fail: a recorded scientific failure without a further base. The 20-block SEM of `escaped_neutral` is then documented as under-covering. The Amendment 13 failure stands in either case.
- **Literature.** The Hall 1992 and Efron–Tibshirani 1993 statements are paraphrases until their page references are verified. Only DiCiccio & Efron 1996 has been verified verbatim.

*(i) V8-nuc sample (`nuc_local`).*
- **T12 rule unchanged** (`ionmc.transport.parity`). A bin is supported iff ≥ max(2, ⌈n_batches/2⌉) batches are non-empty in each sample. Support therefore depends on the histories per batch h, not on N.
- **Python sample.** **4 shards × 20 batches × 1500 histories = 1.2e5** (hr5c r_index 1, shards k = 0..3). The cpu64 side is unchanged (100 × 1e4).
  - The occupancy model reproduces the archived verdict (0.259 exact rule vs 0.253 model at h = 600). It gives a singleton mass fraction of 0.98–0.99 at h = 1500 and P(conclusive) ≈ 1.00, against 0.00 at h = 600. Elastic adds recoil deposits and can only raise the occupancy.
  - Planning time: 3e4 histories per shard at 60 hist/s (hr5 python shards: 1.2e4 in about 200 s) is about 500 s, or about 650 s with a 1.3× elastic overhead, below 1500 s.
  - If the rehearsal shows > 1500 s, the sample becomes 8 shards × 10 batches × 1500 (same h, same 80 batches, same N).
- **Records.** The hr5 inconclusive record is preserved, and the new verdict is reported next to it.

*(j) Seeds and r_index (orchestrator D-5).*
- **Bases.** The seed of row r and shard k is base + 1000·r_index + k (Amendment 6). lv5c qualification base **20481004**: the next unused value of the Amendment 6 sequence; 20441004 and 20471004 are consumed. hr5c base **20491004**, declared here.
- **Rehearsal family.** **2050xxxx**, restricted to **20505000–20509999** so that it stays disjoint from lv5c r_index 19 and 20 (20500004 + k and 20501004 + k). Rehearsals are never evidence. The planning numpy seed 20500001 of `estimate_elastic_bgg.py` was not a transport seed.
- **Fresh bases after a reviewed fix** ((e)9, (f)), declared now: lv **20511004** and hr **20521004** (20491004 is hr5c and 2050xxxx is the rehearsal family), then 20531004, 20541004, … Each further base is declared in an amendment before it is observed.
- **lv5c r_index:**

  | r_index | Row |
  |---|---|
  | 1 | throughput |
  | 2 | V8-LV |
  | 3 | R1-D2 |
  | 4–7 | V5 (150 on / 150 off / 200 on / 200 off) |
  | 8 | V2b-el |
  | 9 | V7-scan-el |
  | 10 | V7-shift-el |
  | 11 | unused |
  | 12 | V6 |
  | 13 | E1-B |
  | 14 | V11 |
  | 15 | P7 |
  | 16 | V7-R |
  | 17 | V7-R diagnostic |
  | 18 | V10-A cross-check |
  | 19 | X-ENDF |
  | 20 | V6-sens |

  V10, D9-BGG and the V10-A table comparison are deterministic.
- **hr5c r_index:** 1 V8-stat (including V8-nuc), 2 V7-f32, 3 V6-CUDA. The draft's optional r_index 4 (V7-R on CUDA) is dropped, because its seeds would equal lv5c V11's.
- **Seed overlaps, stated because they exist.**
  - lv5c r_index 1 (throughput, 20482004 + k) equals the consumed lv5b V7-replicate seeds. The throughput step yields no row verdict.
  - lv5c r_index 12/13 equal hr5c r_index 2/3 (V6-LV with V7-f32, E1-B with V6-CUDA).
  - lv5c r_index 11 is left unused because it equals hr5c r_index 1 (V8-stat).
  - lv5c r_index 10 (V7-shift-el, 20491004 + k) shares only the base value of hr5c, which no hr5c row uses (r_index ≥ 1).
- **Rule against seed overlap.** No statistic combines two samples whose seed sets intersect. Every step records its seed list, and every combine/compare step checks that the sets are disjoint and that none is a consumed evidence seed (fail closed).
- **TOPAS seeds** as in (f).

*(k) Execution constraints (unchanged in substance).*
- The single-process directive applies: `cpu_workers=1`, OMP/OPENBLAS/MKL/NUMEXPR = 1, `IONMC_SINGLE_PROCESS=1`. The summary is `conformant: false` with the deferred partition row until the directive is lifted.
- Every host step is sized for ≤ 1500 s of wall time, because the MCP client aborts at 1800 s idle (RUN-20261008T201947Z-6a3036f4). Rows are sharded, and history counts and criteria never decrease.
- Partials come only from the protected host-runner stdout (attested manifest, content sha256; Amendment 5). Steps run only through the CLI (`run_suite.py --suite lv5c|hr5c --only <step>`).
- The per-process nuclear device cache (C37) also holds the elastic tables, so there is no per-block rebuild. A fault handler is set in every step.
- CUDA chunks are ≤ 2^14 histories.

*(l) Acquisitions (`acquire_reference_data`; sha256 at acquisition; `record_scientific_event` for each failure; roles via `assign_reference_data_role`).*
- **Tool role vs plan labels.** The data tool accepts only the roles {construction, calibration, evaluation, exploratory}. The acquisitions were made with these tool roles:
  - **evaluation:** the pre-1997 angular entries (C0063, C0141, C0148, C0550, O0274, O0553, C0055, C0057, C1420), 13753, O0579 and the arXiv:1409.1938 e-print;
  - **exploratory:** 13569 and the p-p entries (O1226, O0145, C2637, C2606, O2057);
  - **construction:** the six Geant4 files.

  The finer labels of this amendment ("related-model evidence (shared/related lineage)", "report-only", "evaluation with shared beam lineage (calibration-type inputs)") live in the acquisition rationale, the registry entries and the ledger text, not in the tool's role field. Where the two differ, the label of this amendment decides how a row is judged. Acquisition ids and hashes are registered in the repository by step C1b.
- **EXFOR (CC BY 4.0).**
  - Entries: C0063, C0141, C0148, C0550, O0274, O0553, C0055, C0057, C1420 (angular: related-model evidence; tool role evaluation); 13753, O0579 (evaluation); 13569 (report-only); O1226, O0145, C2637, C2606, O2057 (p-p, report-only). D0356 is already pinned.
  - Source: the IAEA raw master `https://www-nds.iaea.org/nrdc/exfor-master/entry/<x>/<entry>.txt` (the path pattern used for C1862 in `src/ionmc/data/exfor_manifest.json`). The recorded fallback is the GitHub mirror `https://raw.githubusercontent.com/IAEA-NDS/exfor_master/main/exforall/<first three characters>/<entry>.x4`. The source actually used is recorded per entry.
  - The mirror sha256 values from planning (for example C0063 `bb27336c…`, 13753 `a254d71c…`) are cross-checks only. If the two sources give different parsed DATA blocks, this is recorded as a scientific event. No intervention is planned.
- **Geant4 11.4.2 sources** (Geant4 Software License, notice recorded verbatim): `G4BarashenkovData.hh`, `G4BGGNucleonElasticXS.cc`, `G4NucleonNuclearCrossSection.cc`, `G4ComponentBarNucleonNucleusXsc.cc`, `G4HadronNucleonXsc.cc`, `G4HadronElasticPhysics.cc`.
  - Role: construction (σ_el(p+A); S(E) for pp above 150 MeV) and the D9-BGG lineage record. Numbers are hand-transcribed and no code is copied.
- **PDG rpp2022 `pp_elastic.dat`** (evaluation, V10).
- **arXiv:1409.1938 source** ((g)).
- **ENDF/B-VIII.0 neutron O-16/C-12**: only if (ii′) triggers.
- **SAID SP07**: optional, report-only.

*(m) Work breakdown (each commit goes through task-commit, `start_codex_review` / `inspect_codex_review`, fixes and re-review until passed; a Codex approval is not scientific validation).*
- **C0:** this amendment and the decision 0041 slice-C section.
- **C1:** acquisitions (l) with registry roles and licences.
- **C2:** the data layer — LAW=5 parser and Hale reconstruction, BGG σ_el and S(E), radius rule and Bessel shape, the builder with provenance. It also includes two checks:
  - that the LA150 O-16 MT5 used for R is not itself a C-12 copy (a copy is a recorded finding before any table is used);
  - that the MCS scattering power contains no nuclear term that elastic would double count.

  Gate: P6.
- **C3:** the planning estimate is re-issued with the built tables (deterministic) and recorded in decision 0041 before transport.
- **C4:** the python reference — Σ_tot thinning, channel branch, elastic events, recoil push, `nuclear_local`, balance. Gates: P7 and R1-D2 (B).
- **C5:** Warp twins, kernels and device cache. Gates: P5-ext and V8-LV.
- **C6:** C1-ext, the R1-D2 record, the V7-R sidecar and the sub-block keying check, the V8-nuc configuration, `steps_v5c.py` (lv5c/hr5c suites, seed guards of (j)), and the plan-row accounting of (d). Gate: CI.
  - Until the lv5c/hr5c steps exist, the Amendment 14 rows are carried in the explicit allowlist `PLAN_ROWS_PENDING_STEPS` of `tests/ionmc/test_run_suite.py`. The plan-row accounting check skips exactly these ids, and a companion test requires every listed id to be a plan row. Each id is removed when its step lands, and **the set must be empty by C9**, before the V3-005C evidence commit.
- **C7:** TOPAS rows 10–13 through `run_reference_calculation`, with fingerprint checks.
- **C8:** E₀ computed deterministically from the range table ((g)); a 2050xxxx rehearsal (N_V6, V7-R calibration on both surrogates, shard timing), then **Amendment 14a** (sizes, E₀ and quantiles only; nothing decreases).
- **C9:** lv5c at a reviewed exact head, base 20481004. The deterministic rows run first (V10, V10-A, X-ENDF, D9-BGG). Then come R1-D2, V8-LV, P7, V11, V5, V2b-el, V7-scan-el, V7-shift-el, V6, V6-sens, V7-R with its diagnostic and calibration, and E1-B.
  - If the (ii′) trigger fires after a seeded step has run, the base is consumed.
- **C10:** hr5c at the same head, base 20491004 (V8-stat with V8-nuc, V7-f32, V6-CUDA). Then archives, README tables, the decision outcome, the ledger, and an exact-SHA `record_local_validation`.
  - PR #82 is merged only after the gates of (a)6 pass.

**Amendment 14 addendum (2026-10-10, C2 data findings, before any transport result).** (1) *p-p lower energy bound.* The LAW=5 reconstruction of the LA150 H-1 evaluation gives a negative nuclear-plus-interference density at the 16.26° CM cut below about 20 MeV (−5.6 mb/sr at 10 MeV; positive at 20 MeV and above). The elastic table builder therefore records E_min,pp, the lowest grid node from which every node and midpoint has a non-negative transported density at and above the cut, sets σ_pp = 0 below E_min,pp, and writes the value and the first-negative diagnostics into the table JSON; the P6 negativity rule applies at and above E_min,pp. The value is reported in the C2 record. (2) *BGG below 14 MeV.* The Geant4 rule below 14 MeV (the 14 MeV value times the Coulomb factor of `G4NuclearRadii` with the `r0[Z]` table and `aeff[Z]` of `G4IsotopeList.hh`) is reproduced as Geant4 implements it; the two files are acquired as construction data (sha256 53cae34b…, d343b827…). The decision-0041 planning value of about 1060 mb for p+O at 10 MeV was a constant-below-14 MeV estimate; the table JSON records σ_el at 10 MeV per target for comparison. (3) *PDG p-p compilation* (rpp2022-pp_elastic.dat, sha256 41268d80…) is acquired as construction data for the above-150 MeV scaling ratio S(E); it is the same compilation as the report-only V10 comparison (shared lineage, (c)3). (4) *Conventions recorded by the builder:* target set and indices aligned with the nuclear table (seven LA150 targets plus H-1, the same element mapping and (A/A_ref)^{2/3} scaling); the p-p integral counts each event once (half sphere 0 ≤ μ_CM ≤ cos 16.26°, factor 2π), which reproduces the research values (49.4/27.4/26.9 mb above 10° CM at 50/100/150 MeV); union grid of the H-1 ENDF and Barashenkov nodes with adaptive midpoint refinement at 1–14 MeV; 256 equiprobable inverse-CDF bins per node (the nuclear table uses 64) to resolve the diffraction minima. (5) *Correction of a transcribed test value (2026-10-10, before any result).* Row P6 (5) quotes p+Ca σ_el = 912 mb at 100 MeV; `G4BarashenkovData.hh` at 0.1 GeV gives n-total 1340 mb and p-inelastic 573 mb, i.e. σ_el = 767 mb (912 mb is the 90 MeV node, 1500 − 588: a hand-reading slip in the research note). The P6 test value is 767 mb; the row text is left unchanged and read with this correction. The 150 MeV value (400 mb) stands. The builder also records σ_el at 10 MeV per target (C 899.5, N 965.6, O 1021.1, Ca 1281.5 mb by the Geant4 rule) next to the decision-0041 planning value of about 1060 mb for p+O. (6) *Shape lower bound (2026-10-10, C2 builder finding, before any transport result).* The black-disk radius R(E) = √(σ_nonel(E)/π) has no real value where the transport's σ_nonel is zero near the non-elastic threshold, so the builder records per target E_min,shape, the lowest grid node from which every node has a valid radius, sets σ_el = 0 below it and writes the values into the table JSON (first build: C 6.0, N 3.25, O 6.58, Al 2.59, Si 3.97, P 2.48, Ca 5.74 MeV; protons below about 7 MeV have a residual range in water under 0.5 mm). E_min,pp of (1) was 12.53 MeV in the first build (negative density at 12.5 and 12 MeV recorded). The radius from the LA150 σ_nonel at 150 MeV is R_O = 2.68 fm against the planning value 2.60 fm (recorded, not corrected). R uses the transport's own σ_nonel including its extension above 150 MeV (decision 0041), as (b) states.

**Amendment 15 (2026-10-10, after Codex review of the C2 builder (findings F1-F4), orchestrator decisions D1-D4, before any transport result) — P6 failure record, pre-registered elastic model domain, construction-source correction and the two frozen C2 checks.** Nothing in this amendment changes a frozen tolerance, the text of a frozen row, the text of Amendments 1-14 or of the Amendment 14 addendum. Addendum items (1), (3) and (6) are to be read with sections (a), (c) and (b) below. No ionmc transport with elastic scattering has been run; the TOPAS runs for rows 10-13 (step C7) that are in progress are reference data, are not ionmc results and depend on none of these decisions. Decision 0041 "Slice C" receives a dated paragraph that points here.

*(a) Row P6 item (4): recorded failure for H-1, and a separately justified p-p domain (F1, D1).*
1. **Record.** Row P6 item (4) was frozen as "failure on any negative transported density at any node or node midpoint". The LA150 H-1 LAW=5 LTP=1 reconstruction (ENDF-102 eqs 6.9/6.10/6.14) gives a negative nuclear-plus-interference density at the 16.26° CM cut (|μ_CM| = 0.96) below 12.53 MeV: −0.022 mb/sr at 12.5 MeV, −0.75 mb/sr at 12 MeV, −5.6 mb/sr at 10 MeV. **P6 item (4) is therefore FAILED for the H-1 evaluation below E_min,pp = 12.53 MeV.** The failure stays in the record. Addendum (1), which applied the negativity rule only at and above E_min,pp, is a domain change made after the observation; it is not a pass of the frozen item and is superseded in that respect by this section. Ledger wording (`docs/requirements/…` ledger and results README): "P6: fail for H-1 below E_min,pp (frozen item (4); negative LA150 Hale NI density at the 16.26° cut, −0.022 / −0.75 / −5.6 mb/sr at 12.5 / 12 / 10 MeV); pass above (revised domain, Amendment 15 (a))". Archive table: a row "P6 item (4), H-1 < E_min,pp | table JSON `elastic-proton-<id>.json`, field `negative_density_below_e_min_pp` | C2b head | n/a (deterministic) | **fail**, recorded; revised-domain row P6-D below".
2. **Revised model domain (pre-registered, separately justified).** The p-p elastic channel is defined on [E_min,pp, 250] MeV with σ_pp = 0 below. E_min,pp is the lowest grid node from which every node and node midpoint up to 150 MeV has a non-negative density at and above the cut (first build 12.53 MeV). The justification does not rest on the failure:
   - (i) *Separability.* The Coulomb part of p-p scattering is transported by the MCS model. Below E_min,pp the LA150 Hale LTP=1 representation gives no non-negative nuclear-plus-interference density at the 16.26° cut, and the model of Amendment 14 (b) has no Coulomb-nuclear interference term. The nuclear p-p part below E_min,pp therefore cannot be separated from the Coulomb part without an additional, unvalidated model, and it is **omitted**. No claim is made that it is small: the domain revision is justified by this separability argument and by the quantified impact of (ii)-(iii) only. The reconstruction's σ_pp(15 MeV) and σ_pp(20 MeV) above the cut are reported by C2b from the table JSON (`sigma_barn[0]` at the nearest nodes and direct evaluation).
   - (ii) *Residual range.* A proton at 12.53 MeV has a CSDA range in water of about 0.19 g/cm² (ICRU 90/PSTAR level; Bragg-Kleeman 0.193 g/cm²). C3 records the value from the project's own range table (Bethe source, water) in the re-issued planning estimate. Every omitted p-p event displaces energy by at most this residual range, about 1.9 mm (not sub-millimetre).
   - (iii) *Omitted events per 150 MeV history (the bound of record),* N_pp = ∫ n_H σ_pp(E) / S(E) dE over E < E_min,pp. Here n_H = 6.69e22 cm⁻³ (two H atoms per water molecule). σ_pp below E_min,pp cannot be taken from the reconstruction, because its density is negative there; σ_pp is therefore bounded by S-wave unitarity, dσ/dΩ_CM ≤ 1/k². Over the transported half sphere μ_CM ≤ 0.96 this gives σ_pp ≤ 397 mb at 12.53 MeV, growing as 1/E below. With the residual range of (ii) this gives **N_pp ≤ 1.2e-2 omitted p-p events per 150 MeV history**, carrying **≤ 0.15 MeV** in total (≤ 0.1 % of the beam energy). The displacement of that energy is bounded by the 1.9 mm residual range (0.193 g/cm²). This bound stands as the record. C3 re-evaluates it with the project's range table and records those values; it does not replace the bound by a reconstruction value.
3. **Fail-closed bound.** A build with E_min,pp > 15 MeV fails (`BuildError`), because (ii)/(iii) were evaluated at 12.53 MeV and a residual range of 0.27 g/cm² at 15 MeV is the largest one this justification covers. A build with E_min,pp > 150 MeV already fails (existing `BuildError`).
4. **New row P6-D (revised domain; the P6 row text is unchanged).** Pass iff: (1) every node and node midpoint in [E_min,pp, 150] MeV, plus the fixed 150 MeV shape used above 150 MeV, has a density ≥ 0 over |μ_CM| ≤ 0.96. This is evaluated in a pass separate from the one that determines E_min,pp, so the flag is not true by construction. (2) σ_pp = 0, `valid` = 0 at every node below E_min,pp, and the loader asserts both. (3) E_min,pp within the bound of 3. (4) Synthetic fixtures: a density made negative at one in-domain midpoint fails the build, and a density negative only below the domain is recorded as `negative_density_below_e_min_pp = true` without failing. Tolerances: exact (sign tests), as in P6. P6 items (1)-(3) and (5)-(7) are unaffected and keep their tolerances.
5. **Merge gate.** P6 (frozen) is recorded as FAILED for H-1 below 12.53 MeV and stays failed in the ledger. P6-D, the pre-registered revised-domain row, replaces P6 in the merge-gate list of Amendment 14 (a)6. **This is a model-domain revision after a recorded data-quality failure, not a tolerance change; the recorded P6 failure is preserved and disclosed in the PR and the merge-gate report.** No operator intervention is raised for it; the operator reviews the PR at merge time.

*(b) p+A elastic domain declared, not silent (F2, D2).*
1. **Domain.** The p+A elastic channel is defined per target on [E_min,shape(target), 250] MeV (first build: C 6.0, N 3.25, O 6.58, Al 2.59, Si 3.97, P 2.48, Ca 5.74 MeV). Below it σ_el = 0. This is a **declared model limitation**: BGG σ_el there is large (O 1021 mb at 10 MeV, addendum (5)) and is not transported. Addendum (6) is read with this section.
2. **Interfaces.** The table JSON carries `elastic_domain = {target: [e_min_mev, 250.0]}` for all eight targets (H-1 from (a)2) and, per target, `sigma_bgg_below_domain_max_mb` (max BGG σ_el below the domain). `ElasticTable.elastic_domain -> dict[str, tuple[float, float]]` reads it, and the loader rejects a table where σ ≠ 0 or `valid` ≠ 0 below the domain (stale/unqualified, fail closed). Transport (python and Warp, both counted identically) adds the diagnostic counter `elastic_below_domain`. It is incremented once per proton, at the first step whose post-step energy is below max over the material's targets with BGG σ_el > 0 of E_min,shape. A second counter, `pp_below_domain`, is incremented once per proton at the first step whose post-step energy is below E_min,pp in a material that contains hydrogen. Both counters are diagnostics: they never set `valid = False`, are reported with every nuclear run, and exist only with elastic on, so R1/A16 for `nuclear=False` is untouched. They are part of the P5-ext/V8-LV parity of the counters.
3. **Quantified impact (before transport).** N_el,O = n_O σ_el R(E_min,shape,O) with n_O = 3.34e22 cm⁻³, σ_el ≤ 1.1 b and R(6.58 MeV) ≈ 0.062 g/cm² gives 2.3e-3 events per history in water. Each event concerns a proton with ≤ 0.6 mm residual range, so the dose displacement per event is < 0.6 mm. H contributes no p+A elastic. C3 recomputes ∫ n_t σ_BGG/S dE below each target's domain for water and every `materials.py` tissue with the project range table and records it in decision 0041.
4. **Rejected alternative.** A low-energy angular shape below E_min,shape (isotropic nuclear, optical model) cannot be validated with the data in scope ((c)1 of Amendment 14: no independent p+O angular data) and would add an unvalidated model. Rejected.
5. **Fail-closed rule.** An element with E_min,shape > 10 MeV fails the build (`BuildError`). Justification: the impact bound of 3 scales with R(E_min,shape), and R(10 MeV) ≈ 0.13 g/cm² is the largest residual range accepted.
6. **E-shape (exploratory, report-only, never evidence).** One python run at 150 MeV in water with 2e4 histories per variant: σ_el below E_min,shape off (the model) versus forced isotropic in μ_CM with BGG σ_el (a test-only builder option whose table id differs and which the loader refuses outside this step). Both variants use the same seeds from the rehearsal family: base **20509500**, shard k → 20509500 + k (inside 20505000-20509999, outside the C8 range, which C8 is told to avoid). Reported: ΔIDD(z) in the last 10 mm with batch SE, R80 and R90 shift, and the `elastic_below_domain` count. With about 45 expected affected events in 2e4 histories the comparison is an upper-bound check. Its result changes no row and no tolerance.

*(c) PDG p-p compilation: role and construction source (F3, D3).*
1. **Finding.** The builder hash-verifies `pdg-rpp2022-pp-elastic` and includes it in the table identity (`SOURCE_IDS`), but it never parses the file. S(E) above 150 MeV is computed entirely from the BGG formula `pp_bgg_mb` = 1.0115·(23 + 50·√(ln(0.73/p)^7)) mb, hand-re-implemented from Geant4 11.4.2 `G4HadronNucleonXsc` / `G4BGGNucleonElasticXS`.
2. **Correction.** The construction source of S(E) is the Geant4 BGG p-p formula (registry ids `geant4-g4hadronnucleonxsc-cc-11.4.2` and `geant4-g4bggnucleonelasticxs-cc-11.4.2`, already construction). The PDG rpp2022 compilation is **evaluation (report-only V10, shared lineage with the BGG formula, which Geant4 fitted to PDG-type data)**, as Amendment 14 (l) originally listed it. Addendum (3) is read with this correction.
3. **Implementation (C2b).** `pdg-rpp2022-pp-elastic` is removed from `SOURCE_IDS`, so the table id changes. `assign_reference_data_role` appends role `evaluation` to acquisition b742ae43…, and the registry entry gets `role="evaluation"`, with description and lineage stating that it was registered as construction in error and why (no value is used in construction). Tests: the id is absent from `SOURCE_IDS` and from `info["sources"]`; no module under `ionmc.nuclear.elastic_*` references the dataset id; P6 (5) pp values (27.61 / 23.86 mb) are unchanged.

*(d) The two frozen C2 checks, implemented before the table is used (F4, D4).*
1. **O-16 MF3/MT5 and MF6/MT5 are not a C-12 copy.** The existing `o16_copy_finding` tests MT2 and stays as the documented MT2 finding (C1-ext). It is not the frozen check. New `o16_mt5_copy_check(c12, o16)`:
   - Compare MF3/MT5 σ and the MF6/MT5 yields of n, p, d, α and γ at every shared incident energy in [max threshold, 150] MeV. If fewer than 10 energies are shared, also compare on the union of both node sets by lin-lin interpolation (INT=2, Amendment 3).
   - **Not a copy** iff the relative difference is > 1e-3 at > 50 % of the compared nodes, for σ and for each yield with non-zero mean. An exact copy gives 0, and the MT2 copy was ≤ 1e-4.
   - **In addition**, the coefficient of variation of the ratio σ_O/σ_C (MF3/MT5) across the shared nodes must be > 1e-3, so that a scaled copy (for example (16/12)^{2/3}·σ_C) is detected. Both conditions must hold.
   - A copy raises `BuildError` and is recorded with `record_scientific_event`. Such a copy would also affect the V3-005A non-elastic table, which has the same MT5 lineage.
   - The result goes to the table JSON as `o16_mt5_finding` (n compared, fraction above 1e-3, median and max relative difference, ratio CV).
   - Tests: a synthetic exact copy and a synthetic scaled copy both raise; the real data pass (data case: skip in CI, fail under `IONMC_REQUIRE_DATA=1`).
2. **The MCS scattering power has no hadronic term.**
   - The MCS model is the differential Molière scattering power (Gottschalk, Med. Phys. 37 (2010) 352): T_dM = f_dM(pv, p₁v₁)·(z E_s/pv)²/X_S, with E_s = 15 MeV (`physics/em.py` `scattering_power_dm`) and 1/X_S = α N_A r_e² Σ w_i (Z_i²/A_i){2 ln(33219 (A_i Z_i)^{−1/3}) − 1} (`physics/scattering.py`).
   - Its only inputs are the electromagnetic constants α, r_e, N_A, E_s and the material's Z, A and weight fractions. X_S is a screened-Coulomb (Molière) quantity, and the Gaussian-core power has no nuclear-elastic tail, so the new channel does not double count.
   - Recorded in `docs/architecture/transport.md`, section "Nuclear interactions" (elastic paragraph).
   - Tests:
     - the transitive in-package import closure of `ionmc.physics.scattering` and `ionmc.physics.em` contains no `ionmc.nuclear.*` and no `ionmc.data.*` module;
     - the parameter list of `scattering_power_dm` is exactly (pv_mev, p1v1_mev, z, inv_rho_xs_cm2_g, rho_g_cm3);
     - the existing Gottschalk value test stays.
   - C2b confirms against the paper that f_dM is calibrated to Molière (Coulomb-only) theory, not to measured MCS data. If it is calibrated to measured data, this is recorded as a finding before C4.
3. **Qualification flags renamed (they say what they check).** `QUALIFICATION_FIELDS = ("no_negative_density_in_domain", "o16_mt5_not_c12_copy", "shape_normalised")`. `no_negative_density_in_domain` is P6-D (1)-(2), computed independently.
   - The table JSON adds the diagnostics `negative_density_below_e_min_pp` (bool) with `first_negative_nodes`, `o16_mt2_copy_recorded` (the former `o16_copy_asserted`, now a finding, not a qualification) and `o16_mt5_finding`.
   - `SCHEMA` and `BUILDER_VERSION` are bumped, so the loader rejects every earlier table as stale.

*(e) Work step C2b (fix round, before C3).* Builder, loader, registry and documentation changes of (a)-(d) (the transport counters `elastic_below_domain` and `pp_below_domain` land with C4/C5 and their parity tests); the tests of P6-D, (c)3 and (d)1-2; one rebuild, whose table id, E_min,pp, E_min,shape, `elastic_domain`, σ_pp at 15/20 MeV and `o16_mt5_finding` are reported in the C2b record. The commit goes through task-commit under task id V3-005B, then `start_codex_review` / `inspect_codex_review`, fixes and re-review until passed (a Codex approval is not scientific validation). Gate: P6 items (1)-(3), (5)-(7), P6-D, (d)1-2 and CI. C3 then adds the deterministic impact numbers of (a)2(ii)-(iii) and (b)3. No tolerance changes; no seed is consumed (E-shape uses the rehearsal family).

**Amendment 15 addendum (2026-10-10, row table; after Codex REVIEW-54fe3b7b findings F1 and F3, before any transport result).** Amendment 15 (a)4 and (b)6 introduced the rows P6-D and E-shape in prose only, so the plan-row accounting of `tests/ionmc/test_run_suite.py` (which reads Markdown table rows) did not see them. They are listed here in the format of the criteria table; the prose of Amendment 15 stays the binding detail, and nothing in this addendum changes a frozen tolerance, a frozen row or the text of Amendments 1-15. Both ids enter the allowlist `PLAN_ROWS_PENDING_STEPS` until their steps exist.

| # | Check | Pass iff | Class | Tier | Slice |
|---|---|---|---|---|---|
| P6-D | Revised p-p elastic domain [E_min,pp, 250] MeV (Amendment 15 (a)2-(a)4; the frozen row P6 is unchanged and stays FAILED for H-1 below E_min,pp) | A sign-test pass separate from the one that determines E_min,pp, on the finer μ grid of the builder, finds a transported nuclear-plus-interference density ≥ 0 over \|μ_CM\| ≤ 0.96 at every node and node midpoint at and above E_min,pp (up to 150 MeV, plus the fixed 150 MeV shape used above); σ_pp = 0 and `valid` = 0 at every node below E_min,pp, asserted by the loader; the builder fails (`BuildError`) if E_min,pp > 15 MeV; the synthetic fixtures of (a)4 (4) behave as stated. Exact (sign tests). | self-consistency | CI | C |
| E-shape | Exploratory low-energy shape sensitivity (Amendment 15 (b)6): python, 150 MeV, water, 2e4 histories per variant, σ_el below E_min,shape off (the model) vs forced isotropic in μ_CM with BGG σ_el (test-only builder option); seeds of the rehearsal family, base 20509500, shard k → 20509500 + k | Report-only, never evidence: ΔIDD(z) in the last 10 mm with batch SE, R80 and R90 shift and the `elastic_below_domain` count are reported; the result changes no row and no tolerance | exploratory | report-only | C |

*Correction to Amendment 15 (a)2(iii) and (a)3 (2026-10-10, Codex REVIEW-54fe3b7b F1).* The estimate in (a)2(iii) used an S-wave unitarity "bound", dσ/dΩ_CM ≤ 1/k². That inequality bounds a single partial wave, not the full nuclear-plus-interference amplitude, so it is not a bound on the omitted p-p events and is **withdrawn** (with it the figures ≤ 1.2e-2 events and ≤ 0.15 MeV per history, and the table JSON field `pp_omitted_events_bound_per_history_150mev` of the C2b build f932d6ae…, which carries 1.095e-2). It is replaced by a data-based measure that the builder records: the absolute nuclear-plus-interference correction M(E) = 2π ∫ |ρ_NI(E, μ_CM)| dμ_CM over the transported half sphere 0 ≤ μ_CM ≤ 0.96, evaluated from the LA150 Hale reconstruction at the grid nodes in [1 MeV, E_min,pp], and integrated with n_H = 6.69e22 cm⁻³ along the residual path of a proton slowing down in water from E_min,pp to 1 MeV, with the path taken from the project's own water range table. The residual range below 1 MeV is reported separately and is not extrapolated. Because |ρ_NI| counts negative and positive density alike, M measures the size of the part of the evaluation that is omitted; it is not a cross section that the model would transport. Values (C2c build): **4.775e-3 events per 150 MeV history; 9.23e-3 MeV.** The 15 MeV build limit of (a)3 is a pre-registered ceiling on the revised domain; it is not derived from this estimate or from the withdrawn bound. The separability argument of (a)2(i) and the residual-range statement of (a)2(ii) are unchanged.

*Amendment 15 addendum (2026-10-10, (d)1 specification revision, before any transport result).* The MF6/MT5 part of the O-16 not-a-C-12-copy check of (d)1 was evaluated as written with every compared incident energy (30 shared MF3/MT5 nodes, 7 to 150 MeV) in the denominator and a node with both yields zero counting as not differing. As written it FAILS on the LA150 data: zap 3007 differs from C-12 at 14/30 = 0.467 of the compared energies (16 both-zero nodes) and zap 5012 at 10/30 = 0.333 (20 both-zero nodes); this result is preserved in the table JSON as `mf6_rule_as_frozen` with `passes: false`. The rule is revised as a specification defect, since a node at which both materials have zero yield for a product cannot discriminate a copy from an independent evaluation: for each product the denominator is the number of INFORMATIVE compared energies (at least one material nonzero; relative difference 1.0 where exactly one yield is zero), a product gates only if it has at least 10 informative energies and is otherwise report-only, and every gating product must differ by more than 1e-3 at more than 50 % of its informative energies; the sigma half and all tolerances (1e-3, 50 %, CV 1e-3) are unchanged. The revised evaluation is recorded as `mf6_rule_revised`; `o16_mt5_not_c12_copy` is set from it (sigma test AND every gating product). Under the revised rule every gating product differs at fraction 1.0 (informative counts 10 to 30; zap 5012 has exactly 10 and zap 3007 has 14). Relabelling of the omitted p-p correction (Codex REVIEW-f2767757 F2): the values 4.775e-3 and 9.23e-3 MeV above are total-variation diagnostics of the signed nuclear-plus-interference correction (JSON fields `pp_omitted_ni_correction_total_variation_per_history_150mev` and `pp_omitted_ni_correction_weighted_total_variation_mev`), event-equivalent magnitudes, not expected numbers of physical events and not energy transferred; the signed correction can cancel and the physical sigma_pp below E_min,pp is not available from the evaluation. C2d table of record id 045c31a2b5681b2d8b718c14645d562ecb787c5c81f0871b24bf092d4004d95c (builder `ionmc-elastic-proton-builder-4`).

**Amendment 16 (2026-10-10, after Codex REVIEW-67f87f7f on the C2d head 6d8f125; before any transport result) — the frozen MF6 copy rule stays failed; independently pre-registered follow-up evidence for the O-16 MT5 independence question.**
1. *Record.* The Amendment 15 (d)1 MF6 per-product rule, evaluated as frozen with all compared incident energies in the denominator, FAILS on the LA150 data for two sparse products (zap3007 14/30, zap5012 10/30; 16 and 20 nodes with zero yield in both evaluations). The "informative-node" revision of the Amendment 15 addendum of the same day was written after that outcome was observed; it therefore does **not** qualify the table. It stays in the table JSON as a report-only analysis (`mf6_rule_revised`) and the as-frozen result stays recorded (`mf6_rule_as_frozen`, `passes: false`). The qualification flag `o16_mt5_not_c12_copy` is replaced by two flags: `o16_mt5_sigma_not_c12_copy` (the σ half of (d)1 exactly as frozen: relative difference > 1e-3 at > 50 % of the compared MF3/MT5 energies and ratio CV > 1e-3; observed fraction 1.0, CV 0.058) and `o16_mt5_spectra_independent` (item 2). The MF6 per-product rule as frozen is recorded as a failed safeguard specification (a both-zero node cannot discriminate a copy from an independent evaluation) and is not used for qualification.
2. *Follow-up evidence, pre-registered before computation.* The question the safeguard asks — whether the O-16 MF6/MT5 product data are an independent evaluation rather than a C-12 copy — is answered on a facet not yet examined: the product **energy spectra**. For each product in {n, p, α} and each incident energy in {50, 100, 150} MeV (nine combinations), the normalised angle-integrated secondary-energy distribution of O-16 and of C-12 is reconstructed from MF6/MT5 (LAW=7 or the law the file uses; the same reader as the nuclear table builder) and compared by the Kolmogorov–Smirnov distance D between the two cumulative distributions on the union E′ grid. An exact copy gives D = 0 at every combination; independent evaluations of different nuclei differ in Q-values, level densities and Coulomb barriers, so D well above 0.02 is expected. **Pass iff D > 0.02 at at least eight of the nine combinations and D > 0 at all nine.** The nine D values are recorded in the table JSON (`o16_mt5_spectra_finding`) and the flag `o16_mt5_spectra_independent` is set from this rule only. A fixture with copied spectra must fail; a fixture with spectra differing only by a 1e-6 scale must fail (D = 0 after normalisation); the real data are evaluated once, after this amendment is committed.
3. *Qualification.* The elastic table qualifies iff `no_negative_density_in_domain`, `o16_mt5_sigma_not_c12_copy`, `o16_mt5_spectra_independent` and `shape_normalised` are all true; the failed MF6 per-product rule is carried in the record and in the ledger ("safeguard (d)1 MF6 rule: fail as written; superseded by Amendment 16 evidence"). No tolerance of any scientific row changes. If item 2 fails, the elastic table does not qualify and the task stops for a reviewed physics reason.
4. *Rows.* `X-MT5-SPEC` (self-consistency, CI, slice C) is the row of item 2; it enters the plan-row accounting through the pending-steps allowlist until its CI test lands (it is deterministic and lands in C2e).

| # | Check | Pass iff | Class | Tier | Slice |
|---|---|---|---|---|---|
| X-MT5-SPEC | O-16 vs C-12 MF6/MT5 product energy spectra (n, p, α at 50/100/150 MeV) | KS distance D > 0.02 at ≥ 8 of 9 combinations and D > 0 at all 9 (Amendment 16) | self-consistency | CI | C |

**Amendment 17 (2026-10-10, before any lv5c/hr5c step and before the C8 rehearsal) — definitional clarifications of Amendments 13–16 found while specifying the step producers (C5/C6/C8 briefs). Nothing below changes a tolerance, a history count, a seed or a gate direction; every item fixes a reading of existing text so that the steps are fully determined before data are seen. Where two readings exist, the chosen one is stated and the verdict under the other reading is also reported (contrary evidence is preserved). Amendment 14a remains reserved for the post-rehearsal sizes (N_V6, E₀, the studentised quantiles, E1-B and V8-nuc shard counts); nothing there decreases.**

*(a) V7-R (Amendment 14 (h)).*
1. "Block sums" are the true sums of the escaped neutral energy (neutron + gamma components) over the 500 histories of a block, in MeV; the sidecar states `units: "MeV per 500-history block"`. The JSON replicate `mean` equals `mean(sums)/500` (asserted to 1e-14 relative). The bootstrap-t interval and both gates are scale-invariant, so this is a units declaration.
2. Order statistics of the B = 1999 resampled t*: sorted ascending, 1-based; t*_(α/2) = the ⌊(B + 1)·α/2⌋-th = 317th value and t*_(1 − α/2) = the ⌈(B + 1)(1 − α/2)⌉-th = 1683rd value (the outer, conservative pair).
3. A resample with se*_b = 0: t*_b = +∞ if m*_b > m̂, −∞ if m*_b < m̂, 0 if equal; the count of such resamples is recorded per replicate. (The replicate-level rule "a degenerate se is a miss" is unchanged.)
4. The literal `SeedSequence([20481004, 16, shard, replicate])` and `[20481004, 16, 1000 + s, j]` use the run's seed base, which is 20481004 at evidence; a rehearsal base therefore gives different streams. The generator is numpy `PCG64` via `default_rng`; the numpy version is recorded.
5. Calibration requirements per estimator and surrogate: the two surrogates of (h) (the archive-fitted skewed block distribution and the rehearsal block sums) are `escaped_neutral` surrogates and apply to its bootstrap-t gates. For `sec_p` and `nuclear_local` the rule is the unchanged Amendment 13 rule, and the requirement "power ≥ 0.9 at the nominal" is evaluated on the Amendment 13 Gaussian 12-bin profile surrogates (the sampling distribution of these profile means, bootstrap-checked in Amendment 13). The Amendment 13 skewed-Gamma surrogate's single-gate power of 0.14 stays a recorded limitation of the distribution-free radius. Contrary reading preserved: under the literal reading "each estimator and each surrogate", the row cannot pass for `sec_p`/`nuclear_local` irrespective of the data; the calibration step reports the verdict under both readings, and the row verdict uses the reading stated here.
6. The calibration step may be sharded into `pytest-v7r-calibration-s{k}` (a fixed pytest parameter split, listed in the suite) if one step would exceed 1500 s; the row's calibration condition then reads "every shard passed at the evidence SHA".
7. "At the evidence code SHA" for the rehearsal means: the sha256 of the hashed source set (`source_file_list("lv5c")`) is identical to that of the evidence commit; the rehearsal document records that digest. A commit that changes only plans, decisions, fixtures, README or tests outside the hashed set keeps the rehearsal valid.
8. Sub-block sums (b ∈ {40, 80}): decided in Amendment 14a from the rehearsal timing; not emitted inside the r_index-16 shards; if a separate r_index-17 run does not fit the time budget, the diagnostic is recorded as dropped for b ∈ {40, 80} and reported for b ∈ {5, 10, 20} only.
9. Defect found and fixed before any lv5c step: `steps_v5b.batch_estimates` addressed the nuclear tally block as the last six columns, which is wrong when the three elastic columns follow it; the fixed code addresses columns by name, and a bit-identity test shows the frozen lv5b/hr5 (elastic off) path unchanged.

*(b) lv5c/hr5c suites (Amendment 14 (j), (m)).*
1. R1-D2 (C) measurement rule: paired common-random-number runs at r_index 3 (seed 20484004 + k), warp-cpu float64, 150 and 200 MeV, elastic on versus elastic off, on the V5 IDD geometry (peak/plateau, plateau, total deposit) and on a Cartesian 1 mm lateral grid of half-width 50 mm; "lateral spread" is the dose-weighted rms lateral radius in the depth slab [0.45, 0.55]·R80 of the elastic-off run. A declared direction "matches" iff the paired batch-mean difference has the declared sign with |z| ≥ 3 (batch standard error, 20 batches); |Δ total deposit| ≤ 0.1 % of E₀ is a point criterion with its 90 % interval reported. The comparison with the archived lv5b V5 partials at 7aae5bb1 (independent seeds) is report-only.
2. hr5c r_index 1 shard indices: python k = 0..3 (fixed in Amendment 14 (i)); cpu64 k = 4, cuda32 k = 5, cuda64 k = 6.
3. E-shape (Amendment 15 (b)6) is an lv5c step `e-shape-exploratory` with the absolute seeds 20509500 + k (a documented exception to `base + 1000·r + k`), tagged `exploratory`; its verdict is "executed" and it never enters a gate.
4. Stop rule (Amendment 14 (f)): when `v5-attribution` detects the trigger it writes its document, exits non-zero and the remaining steps of the C9 order still run (they are independent evidence and are archived); no fresh base is taken for V5, no tolerance changes, and the orchestrator raises an intervention. "The task stops" refers to fix-and-rerun, not to the archiving of the remaining rows.
5. Rung 3 of the attribution ladder needs the elastic p + A recoil deposit resolved in depth: the V5 "on" producers score it through a dedicated pseudo-species channel of the recoil deposit (tally-neutral; A16 and the R1-D2 (B) identity unchanged); it is added before any lv5c step (C6).
6. E1-B: 1e5 python histories are sharded as `e1-b-s{k}` with the shard count fixed in Amendment 14a from the throughput step (≥ 2), the total never below 1e5.
7. P6 leaves `PLAN_ROWS_PENDING_STEPS` by mapping to its CI test; the ledger keeps its recorded failure; P6-D gates (Amendment 15 (a)5).
8. Frozen-suite pin (C3b): `steps_v5.nuc_config` carries `elastic=False`; lv5c/hr5c build their configurations through `steps_v5c.nuc_config_c` (elastic on, the pinned table id) and `nuc_config_c_elastic_only` (row V11).

*(c) V6 / V6-sens (Amendment 14 (g)).*
1. Scoring observable: at every data point (r_i, z_j) the dose is the overlap-volume-weighted mean, over the 0.4 cm diameter sphere centred at the point, of the dose of a 0.5 mm × 0.5 mm axisymmetric (r, z) annular tally (dose = energy deposit / (ρ · annulus volume), MeV/g per primary); the overlap weights are computed by ≥ 8³ sub-samples per cell and recorded. The tally implementation (annular tally or equivalent Cartesian strips) is an implementation choice proved tally-neutral (A16 unchanged); the observable is fixed here.
2. r = 0 rule: point criterion |D_ionmc/D_meas − 1| ≤ 0.02 at every non-missing r = 0 depth, with the 90 % interval and the range-shift equivalent reported per point.
3. Sign-test ties (zero difference of the two medians at a depth): excluded from the test, counted and reported.
4. V6-sens nominal: a nominal run at r_index 20 (k shards) with the same seeds as the six beam variants and the two R0 variants (common random numbers: the variants scale the same normal deviates); σ_Δ is the paired batch standard error of the ratio. The r_index-12 V6 sample is the gating sample and is not reused for V6-sens.
5. "x–θ correlation varied by ±20 %" varies the correlation coefficient ρ_xθ = −0.521 (to −0.417 and −0.625) at fixed σ_x and σ_θ.
6. Gating halo set: every non-missing r ≥ 1 cm point (the frozen row), for the per-point, 90 %, median and V6-sens criteria; the ≥ 1e-4-of-maximum subset (the N_V6 sizing set) is reported alongside, not gated separately.
7. If the V6-sens condition fails, V6 as a whole is report-only (literal Amendment 14 (g)); the r = 0 verdict is still computed and reported under its 2 % rule, and the merge gate (a)6 treats V6 as "status report-only, recorded".
8. Backend disagreement: V6 fails if the lv5c (warp-cpu f64) and hr5c (warp-cuda f32) verdicts differ on any gate component (r = 0, halo, overall); the per-point paired z between the backends is reported, with any |z| > 4 listed as a finding.
9. Beam source capability: `PencilBeamSource` gains an angular spread and an x–θ correlation (per plane, bivariate Gaussian) with defaults that draw no random number, so that every existing configuration stays bit-identical (A16, R1-D2 (A)/(B)); the direction convention is (tan θx, tan θy, 1) normalised. If the source RNG layout cannot accommodate the two extra deviates without shifting an existing stream, the row is blocked and recorded before any run.
10. Role disclosure: the `tbl:parms` beam parameters are construction inputs of the V6 beam (Amendment 14 (g)); the construction role is appended to the arXiv:1409.1938 acquisition's role record and to the registry description; the dose table stays evaluation and is read only by the comparator (isolation test).

*(d) Row–step mapping.* The pending-row allowlist is reduced in the order: CI rows (P6, P6-D, X-MT5-SPEC, P5-ext, C1-ext) → lv5c/hr5c registrations (R1-D2, P7, V2b-el, V7-scan-el, V7-shift-el, E1-B, E-shape, V8-nuc) → V7-R, V11, V6-sens/V6 → V10, V10-A, X-ENDF, D9-BGG → empty before C9.

**Amendment 14a (2026-10-11, after the C8 rehearsal at base 20505000 and before any lv5c/hr5c evidence step) — post-rehearsal sizes, E₀ and the studentised quantiles. Nothing here decreases a count, changes a tolerance, a seed base or a gate direction; every value is fixed from the rehearsal records named here and is then frozen. The rehearsal ran at the reviewed head c847cc75 (Codex REVIEW-e1853ef7 pass); the hashed-source digest of the evidence commit is the one recorded by the V7-R rehearsal fixture (Amendment 17 (a)7).**

*(a) V6 / V6-sens (Amendment 14 (g), Amendment 17 (c)).*
1. E₀ = **177.84075737 MeV** (bisection to 1e-6 relative on R_CSDA(E₀) = 21.321 cm in the project's water range table; the exact-inverse cross-check differs by 2.14e-7 relative), σ_E = **0.8233 MeV** from dR/dE at E₀. Finding recorded, E₀ still used: E₀ lies inside the 176–178 MeV band but outside the 177 ± 0.5 MeV expectation; the range-inversion residual 1.4503e-6 of the transport table is a report-only finding (decision 0041, revised 2026-10-11). ionmc's own d80 at E₀ is reported next to R0 by `v6-compare` and never used to adjust E₀.
2. Rehearsal shard `v6-s0` (RUN-20261011T030353Z-0b25db97, 1e6 histories, warp-cpu float64, 184 s, valid): the halo clause of the history-count rule (153 points with measured dose ≥ 1e-4 of the r = 0 maximum) gives **88** shards (worst relative SE 28.1 % at z = 22.0 cm, r = 3 cm). The r = 0 clause (34 gating points) is satisfied with ≤ 18 shards up to z = 21.8 cm but needs 116, 651 and 8910 shards at z = 22.0, 22.2 and 22.4 cm and cannot be satisfied at z = 22.6 cm (model dose zero in 1e6 histories). The rule therefore defines no finite N_V6.
3. Rehearsal shards `v6-sens-nominal-s0` (RUN-20261011T030829Z-537a9857) and `v6-sens-st-hi-s0` (RUN-20261011T031141Z-88a536e0), 185 s each, common random numbers: of the 261 gating halo points of Amendment 17 (c)6, 32 had zero model dose in both 1e6-history runs (z ≥ 20 cm; measured 6e-6 to 4e-5 of the r = 0 maximum), so the V6-sens ratio was undefined there in the rehearsal; among the 229 points with a finite paired SE the worst needs 345 shards per variant by SE scaling (SE ∝ N^−1/2); on the 153-point sizing subset the worst paired SE 0.42 gives k = 32 shards per variant.
4. Consequence (revised 2026-10-11 after Codex REVIEW-7b4d96c9 and REVIEW-ea1441ef, before any V6 evidence shard): the fail-closed condition already defined applies, justified by a resource bound. For the r = 0 clause, SE scaling from the rehearsal (SE ∝ N^−1/2) needs 8910 shards (8.9e9 histories) at z = 22.4 cm, and at z = 22.6 cm no depositing history was observed in 1e6, so the depositing rate is below 3.0e-6 per history (95 % upper bound) and the ≥ 4e4 depositing histories that a 0.5 % relative SE requires need ≥ 1.3e10 histories; at the measured 5437 hist/s these are ≥ 19 and ≥ 28 days of host time for one sample. The rule therefore defines no feasible N_V6, and **V6 fails closed at the sizing stage**: `V6_SHARDS` stays `None`, the steps `v6-unsized` (lv5c and hr5c) record this reason with the rehearsal run ids and the bound, and rows V6 and V6-sens are recorded as **failed — not sized under the pre-registered rule** (no shard is run, no reading is amended, no tolerance changes). For V6-sens the same scaling gives ≥ 345 shards per variant on the finite points, and at the 32 points undefined in the rehearsal the rate bound above means that even 180 depositing histories per point (a paired SE of 0.075 at best) need ≥ 6e7 histories per variant, with the SE scaling far above that; the sensitivity condition could therefore not be established within any feasible budget either. The intervention IR-20261011-031924-C205C8, raised earlier the same day before this resolution, is retained as historical provenance only; the scientific course is fixed here autonomously and is not conditional on an operator preference.
5. Contrary-evidence note: while selecting the points of the sizing rule the orchestrator displayed measured and model mean doses side by side for the rehearsal shard; no residual, interval or verdict was computed and nothing was changed as a consequence; the evidence bases 20481004 and 20491004 are unobserved.

*(b) V7-R (Amendment 14 (h), Amendment 17 (a)).*
1. Archive studentised quantiles (7aae5bb lv5b partials, 7200 replicates; fixture `tests/ionmc/fixtures/v7r/lv5b-7aae5bb-escaped-neutral-replicates.json`, sha256 14306bbc…): 2.5 % / 15.87 % / 84.13 % / 97.5 % = **−2.2282 / −1.1016 / +1.0360 / +2.0047**; the skewed block surrogate is the zero-inflated Gamma fitted to these four values (rehearsal fit k = 5.0703, π = 0.21667, max |Δquantile| 0.0377, seed 20505000).
2. Rehearsal block-sum surrogate: shards `v7r-s0` and `v7r-s1` at base 20505000 run at the evidence hashed-source digest, committed as `tests/ionmc/fixtures/v7r/rehearsal/` with `PROVENANCE.json` (host run ids, partial digests, digest, seed base, code SHA) in the commit that follows the last change of a hashed source file; a CI test loads the committed fixture in required mode against the current digest whenever the directory exists, so a later hashed-source change without a re-run is caught before the suite. The shard `v7r-s0` of RUN-20261011T013718Z-1fc90c28 (992b9ec, 758 s for 4.5e6 histories) is timing evidence only (worktree edited during the run).
3. Calibration sizing (Amendment 14 (h)4). Rehearsal of the archive-fitted surrogate at the default 200 replicates per margin (RUN-20261011T032408Z-08dd0427 at c847cc75, 1018 s): joint false acceptance 3/200 (single_high) and 2/200 (paired_low) → CP95 upper 0.038 / 0.031 > 0.02; power 187/200 = 0.935 → CP95 lower 0.899 < 0.9; paired gate 9/200 and 6/200, single gate 6/200 at their own margins. The rehearsal misses the requirements at n = 200, and the replicate count is raised to **3200 per margin** (`run_suite.V7R_CAL_REPS`, exported as `IONMC_V7R_CAL_REPS` and verified in every calibration step header; B = 1999 and the 5000 bootstrap checks unchanged; never decreased). At n = 3200 the joint requirement is met iff ≤ 50 of 3200 are accepted at a margin (observed rate 0.0156, CP95 upper 0.0197 ≤ 0.02; 51 gives 0.0201) and the power requirement iff ≥ 2909 of 3200 at the nominal (CP95 lower ≥ 0.9); the bounds table of `v7r.py` carries the n = 3200 row.
4. Gate-alone requirement (revised 2026-10-11 after Codex REVIEW-7b4d96c9; the first version of this item, which made the gate-alone clause a construction statement, is withdrawn as a post-rehearsal weakening). The criteria stay exactly as implemented and reviewed before the rehearsal (C5): joint — CP95 upper ≤ 0.02 at every non-nominal margin; **each gate alone — CP95 upper ≤ 0.04 of its Monte Carlo false-acceptance rate at its own margin** (at n = 3200: ≤ 109 of 3200 accepted); power — CP95 lower ≥ 0.9 at the nominal. The observation that a gate built at α_tost = 0.04 has a false-acceptance rate near 0.04 by construction, so that this clause may be unprovable by Monte Carlo at any n (rehearsal: paired 9/200 and 6/200, single 6/200), is recorded as a limitation and reported next to the verdict; it does not change the gate. If the clause misses at the evidence SHA, V7-R fails ("uncertainty-coverage evidence not established") and the finding is preserved.
5. Execution. One Monte Carlo replicate costs ≈ 4.86 s with all five margins (nominal, paired_low, paired_high, single_low, single_high) evaluated from one shared bootstrap, plus ≈ 42 s per shard for the fit and κ (the rehearsal and a sandbox smoke shard). The calibration is therefore sharded into `pytest-v7r-calibration-s{k}`: shard = (surrogate, part), 15 parts per surrogate (214/213 replicates, ≈ 1080 s each; s0–s14 archive-fitted, s15–s29 rehearsal block sums; data and bootstrap streams seeded by `SeedSequence([20507000, surrogate_idx, part])` and its bootstrap analogue), every shard evaluating all five margins; the Amendment 13 Gaussian/skewed cases of `sec_p`/`nuclear_local` and the readings report form the last shard s30 (308 s measured, RUN-20261011T034222Z-16d64100). `V7R_CALIBRATION_SHARDS` = 31. Each shard writes a sealed count record; the consumer `v7r-calibration-combine` requires every record present and consistent (same fit and κ within a surrogate, parts summing to 3200 per margin, replicate count 3200), sums the counts per (surrogate, margin) and applies the requirements of items 3–4 (the gate-alone CP95 upper bound included); "the calibration step passes at the evidence SHA" means every shard passed and the combine verdict is pass, and `v7r-combine` requires that status partial.
6. The calibration is deterministic (fixed fit/boot/data/κ seeds 20505000–20508000 and the part seeds above, independent of the suite seed base), so the C9 run at the evidence SHA is the calibration of record; a miss there fails V7-R ("uncertainty-coverage evidence not established", Amendment 14 (h)4) and no further count increase follows.
7. Sub-block sums b ∈ {40, 80} (Amendment 17 (a)8): the keying check passed (`test_sub_block_run_range_250_plus_250_equals_the_500_history_block_bit_for_bit`), so emission is feasible without an RNG change, but a separate r_index-17 run of the 16 shards (≈ 16 × 810 s) does not fit the C9 host budget; b ∈ {40, 80} are **dropped and recorded**; the diagnostic reports b ∈ {5, 10, 20} only.
8. Rehearsal seed reading (recorded, not a change): the rehearsal shards use base 20505000 with the shard seed base + 16000 + k; the evidence seeds are 20481004 + 16000 + k as frozen.

*(c) Shard counts (Amendment 17 (b)6, Amendment 14 (i), (k)).*
1. `E1B_SHARDS` = **3** (python, 150 MeV; 3 × 33340 = 100020 histories, never below 1e5; 33340 because a run's history count must be a multiple of its 20 batches). From the throughput step at c847cc75 (RUN-20261011T034756Z-a78b4deb): python 59.4 hist/s at 150 MeV → ≈ 560 s per shard against the 1500/1.25 = 1200 s effective limit; E1-B was the only undeclared plan (≥ 2 needed).
2. V8-nuc: the python sample stays **4 shards × 20 batches × 1500 histories** (hr5c r_index 1, k = 0..3); the rehearsal shard `v8-stat-python-s0` (RUN-20261011T031454Z-2c505136) took 545 s ≤ 1500 s, so the 8 × 10 fallback of Amendment 14 (i) is not taken.
3. V7-R shards: 16 × 4.5e6 confirmed (813 s planned per shard at the measured block rate 5537 hist/s; 758 s measured at 992b9ec).
4. E-shape (exploratory): 2 × 2e4 python histories in one step (≈ 674 s planned), both variants of Amendment 15 (b)6 at the same absolute seed 20509500 + k.
5. The throughput step at the evidence head must report every declared shard plan as fitting the 1200 s effective limit (`shard_table_ok`); a plan that does not fit stops C9 before the first evidence step, and no size decreases afterwards.

*(d) Records carried into the results README (not plan changes).* The frozen lv5b `r1-nuc` digest no longer reproduces at HEAD because of the added `elastic_mev` field (schema only; every other field identical field by field); lv5c R1-D2 (B) reproduces 3edfce9c… with the field dropped. `CONSUMED_R_INDICES` is range(0, 12) (lv5b r_index 12/13 were never run). R1-D2 (C) seeds: IDD runs k = 0, 1; lateral runs k = 2, 3 of r_index 3. Rung-3 residual energy equivalent (Amendment 14 (f)): Σ over [TOPAS peak − 10 mm, TOPAS R80] of (IDD_on − IDD_TOPAS,full) × bin mass per primary; attributed iff the elastic-recoil deposit (species channel `elastic_recoil`) in that window is ≥ 50 % of it.

**Amendment 18 (2026-10-11, after the operator's response to IR-20261011-031924-C205C8 and before any V6-F observation) — follow-up validation design V6-F for the Gottschalk halo, the merge-gate treatment of "failed — not sized", and the record of the intervention.**

*(a) Operator response (verbatim, recorded through `resolve_human_intervention` on 2026-10-11):* "Considering IR-20261011-031924-C205C8, I acknowledge the subsequent autonomous resolution in Amendmentment 14a. Preserve V6/V6-sens as failed at the sizing stage and do not spend the proposed 20–41 hours on a run that cannot satisfy the frozen requirements. Continue independent work, keeping the applicable merge gates enforced. Develop and independently review a computationally feasible follow-up validation design, explicitly documenting any change in scope or criteria and preserving the original failure. Reconcile this intervention's open status with the subsequent decision, recording this acknowledgment accurately." The intervention is resolved; the autonomous resolution of Amendment 14a (a)4 stands unchanged.

*(b) Row V6-F (pre-registered follow-up; the full text is `validation/plans/v3-005c-v6f-followup.md`, revision 2, with the independent review record appended).* V6-F is a NEW row. It never replaces V6 or V6-sens in the ledger, the results tables or merge gate (a)6; the original verdict "failed — not sized under the pre-registered rule" (Amendment 14a (a)4) stays recorded verbatim. Summary of the frozen design (the file governs in any difference):
1. Observable, beam, E₀, scoring and the sphere averaging of Amendment 14 (g) / 17 (c), unchanged. Every retained tolerance is the frozen one: r = 0 point criterion 2 %; halo 15 % per point, ≥ 90 % of points, median |residual| ≤ 10 %; r-trend sign test p > 0.01.
2. Point set fixed from the frozen inputs, the measured table and the paper's aura-term parameters p23–p25 only (never from model output; the aura parameters are a declared extension of the `tbl:parms` construction role, recorded before execution): depth cut z_F = 21.6 cm (the deepest data depth ≤ R0 + σ_sem = 21.614 cm); r = 0: the 29 non-missing depths z ≤ z_F gate, the five at 21.8–22.6 cm are report-only (D1); halo: non-missing r ≥ 1 cm points with z ≤ z_F and measured dose ≥ 1e-4 of the r = 0 maximum (145 points), minus the points where the paper's own aura term (`tbl:parms` p23–p25) exceeds 5 % of the measured dose (41 points removed; 104 gating points at r = 1–5 cm; D7, a new scope difference: the transport deposits no neutral energy, so those points test the missing aura physics, not the charged halo); the aura share is reported for every point and the 145-point verdict with its pre-declared low bias is a contrary reading.
3. Beam-sensitivity screen (V6-F's own per-point screen, D3; row V6-sens itself stays unchanged and failed): a halo point stays in the gating set iff |D_var − D_nom| / D_meas + 2σ′ / D_meas < 0.15 for every one of the six beam variants (σ′ the paired common-random-number SE of the difference); the model-normalised screen is a contrary reading. Coverage minimum (D8): ≥ 40 gating points, ≥ 3 radii and ≥ 12 sign-test depths after excluding ties, otherwise the halo verdict is "not established (beam-sensitive)", never a pass. Expected after the screen: ≈ 58 points, 16 depths (planning estimate from the beam parameters and the measured table).
4. Sizes: N_F = 16 shards × 1e6 histories per backend (rule value 12 with a 1.25 factor on the rehearsal SE; expected worst relative SE 1.7 % halo, 0.34 % r = 0); sensitivity runs k = 8 shards × 7 runs (nominal + six beam variants, common random numbers; the R0 ± 0.5 mm variants are not run); 72 CPU shards × ≈ 184 s ≈ 3.7 h, every step far below the binding 1200 s effective limit of Amendment 14a (c)5; the CPU/CUDA backend-agreement rule of Amendment 17 (c)8 is part of the pass condition.
5. Execution: ONCE, in a later task that makes the lateral-halo claim, gating there, on that task's fresh base B declared before any observation, with seeds B + 7000 + k (gating, k = 0..15), B + 8000 + k (sensitivity, k = 0..7, shared by the seven runs) and B + 9000 + k (CUDA, k = 0..15); r_index ≤ 9 avoids the cross-base collision of bases spaced 10000 apart; B + 7000 must exceed the largest rehearsal seed actually used (20525999) or a seed-by-seed disjointness proof is recorded; no seed of the original V6 (r_index 12/20 of 20481004, 3 of 20491004, never run) is reused. V6-F is NOT executed in V3-005C (its blindness is weaker than a row frozen before the Amendment 14a (a)5 display, and no V6-F code enters the lv5c/hr5c hashed source set).
6. Contrary-evidence notes carried from the design process: the designer saw the sensitivity SE under both normalisations at one gating point, which bounds the rehearsal model-to-measured dose ratio at (21.6 cm, 3 cm) by ≲ 0.37 (0.037/0.0988), at that point only; no rule was changed afterwards; revision-1 per-point tables are quarantined and never committed or archived next to the revision-2 table; recorded in the design file.
7. Independent review (record appended to the design file): first review of revision 1 — "approve with listed changes", three blocking findings B1–B3 resolved in revision 2 as items 2, 3 and 5 above; re-check of revision 2 — "no blocking issue; approve with listed changes" (R2-I1 quarantine done; R2-I2 is item (c) below; R2-m1–m5 applied in the design file); check of item (c) — approved with the wording changes applied here.

*(c) Merge gate (a)6 and "failed — not sized".* Amendment 14 (a)6 reads: "The task is merged only if the scientific gates (V5, V6 with the status fixed by V6-sens, V7-R including its calibration step, V8-nuc, V10-A) and the implementation gates (P6, P7, P5-ext, C1-ext, R1-D2, V8-LV, CI) pass" (P6 replaced by P6-D per Amendment 15 (a)5). Amendment 17 (c)7 covers only the report-only status. A row recorded as "failed — not sized under the pre-registered rule" is a failed scientific gate: it is not a pass and not report-only, and V6-F cannot stand in for it. V3-005C therefore cannot be merged on V6 under (a)6. V10-A is likewise a recorded failure: step `v10-a` failed at 37e87ecf (RUN-20261011T180918Z-e5cdfbe1, partial `v10-a-status` de6e2599…; the (e)9 fallback trigger was evaluated by `v10-a-trigger`, RUN-20261011T181515Z-4e0f6041, and did not fire because clause (c) fails), the step is deterministic and is re-recorded at the evidence head of C9. Both failures are preserved, the implementation gates are evaluated and recorded regardless, and no passing record is created. Default without an operator decision: NOT merged, PR #82 stays open. After C9/C10 the orchestrator raises an ambiguous-boundary intervention that proposes no relaxation, re-labelling or replacement of any gate and asks only how the recorded failures are carried (the V3-005B interpretation-B precedent); any answer is recorded verbatim and changes no row verdict.

*(d) Rehearsal-validity digest (Amendment 17 (a)7, Amendment 14a (b)2) — implementation finding.* The implementation of the "evidence code SHA" digest hashed the full suite source identity, which includes `validation/plans/*.md` and `decisions/*.md`, so every documentation commit after the C8 rehearsal made the committed V7-R rehearsal fixture stale and forced two 4.5e6-history host shards (rebuilt at 6df0e554, 1a86a08c and df596559 before this was found). Amendment 17 (a)7 states that a commit changing only plans, decisions, fixtures, README or tests outside the hashed set keeps the rehearsal valid; the rehearsal block sums depend on code and data only. The rehearsal digest is therefore computed over the hashed source set minus `validation/plans/`, `decisions/`, `docs/` and README files (C9b); the suite source identity recorded in `environment.txt` and checked by `summarize.attest()` is unchanged and still covers every file. Nothing else changes; the fixture is regenerated once at the head carrying this change.

*(e) Addendum to (d) (2026-10-13, C12).* The CI job of the pull request has no data cache; thirty-six slice-C tests that open the nuclear or elastic table lacked the data-availability skip of the slice-B tests and failed in CI at the PR head 4ce33832 (run 38355773758) while passing locally. The skip guards are added (tests only; in required mode the tests still run). Because `tests/` is not part of the shard computation either, it is excluded from the rehearsal-validity digest together with the documentation set of (d); the suite source identity and `summarize.attest()` are unchanged. The fixture is regenerated once more at the head carrying this change.

**Amendment 19 (2026-10-14, before any V3-005D D-E step is run; after the operator's disposition A of IR-20261013-164616-DEFAF4) — V3-005D elastic attribution (exploratory, report-only).** Rows D-E0 to D-E6 are exploratory and never evidence. They change no frozen row, tolerance, base or verdict; V5, V11 and V10-A at dd976c5a stay failed.

(a) **Plan of record:** `validation/plans/v3-005d-elastic-attribution.md`, sha256 `fd9657df047d2e4ef00ad4c0a55c4d390e395abbd342ab200ff80b2964fc7854`, committed in the same commit as this amendment, before D-E0. Supporting already-seen scripts are in `validation/plans/v3-005d-elastic-attribution/`. **Kernels and signatures:** as §2.1–§2.2 of that file. Definitions K_w (w = 1, 2, 3, 5 mm), K*, G, θ_w (w = 2, 3 mm), T (3·θ_med). Thresholds as in the "identified" / "excluded" columns.

(b) **D-E3 probes:** warp-cpu float64, single process, 1e5 histories, 20 batches, V5 geometry r = 20 cm, 150 and 200 MeV. Variants:
- `all`;
- `pp-only` (σ_el = 0 for every target ≠ H-1);
- `pA-only` (σ_pp = 0);
- `pA-sigma0.8` (σ_el × 0.8 for every target ≠ H-1);
- `pA-Rx1.2` (black disk at 1.2·R(E));
- `pA-s0.5` and `pA-s0.9` (sharp form × exp(−q²s²), s in fm);
- `pA-gauss` (exp(−q²R²/4));
- EM-only.

256 equiprobable μ_CM bins per node as in the pinned table; normalisation exact; σ unchanged except where stated.

(c) **D-E4:** rows 14–23 as in §3.2 of the plan of record. Water rows 14–18 and graphite rows 19–23 at 30 / 50 / 100 / 150 / 200 MeV. Primaries per case: **2e7 at 30 and 50 MeV with a 2 mm target, and 4e6 at 100, 150 and 200 MeV with a 1 cm target.** Module `g4h-elastic_HP` only (declared fallback in §3.2); KE filter E₀ − 1e-4 MeV; kinematic-curve classification; seed 20270000 + 1000·row + 1 (k = 2 only if a case is split in two equal halves).

(d) **D-E5:** 1e6 table-sampled events per case.

(e) **D-E6:** swaps `swap-pA` (O CDF from water rows 14–18; graphite rows 19–23 used for D-E5 signatures only), `swap-pp` (H-1 from rows 14–18 above 16.26° CM), `swap-both`; q-scaled interpolation; ionmc shape below 30 MeV; same sizes and seeds as (b).

(f) **Seeds:** investigation family 2099xxxx, starting at **20990100**. D-E3 and D-E6 use one seed per energy: **20990100** (150 MeV) and **20990101** (200 MeV), shared by every variant and by EM-only (CRN). D-E5 uses **20991000 + case index** (0–9, ordered as rows 14–23). The family is disjoint from:
- every consumed base;
- the declared fresh bases 20511004 / 20521004 / 20531004 / 20541004 / 20551004 and their r_index blocks (envelope base + 0 … 20999, maximum 20572003); the re-qualification bases of record are lv 20531004 / hr 20541004 (route document §4.4; 20511004 / 20521004 are retired unused, declared before any observation);
- the rehearsal seeds used (20505000–20525999);
- E-shape 20509500 + k;
- the TOPAS rows.

Disjointness is asserted by a test under `tests/` only (outside the rehearsal digest). The `run_suite.py` seed-guard extension (hashed, l.1303–1313) is deferred to the re-qualification commit, which changes hashed code anyway. The D-E runs are executed by `validation/scripts/elastic_attribution/` and never through `run_suite.py` or an evidence step. Disclosure: numpy seed 20990001 was used once in a planner scratch calculation (analytic V7-R binomial study, no transport, `scripts/candidates2.py`); it lies below the family start and is not used.

(g) **The §5 rule** (rules 1–5, thresholds 0.010, 2·SE, 70 %, 30 %, the signature thresholds, the conditional `swap-sigma` / `swap-all`) and **the §4 pre-commitment** (M2 with the §4.1 inputs, conventions and rules; M3 only on the §4.2 triggers) are frozen with this amendment. No §4 or §5 parameter changes after D-E0 is run.

(h) **Already seen and declared:** the V10-A per-point ratios and first minima, the built-table moments, the shape-moment and removal-kernel estimates of rev 1 / the review (crude; recomputed in D-E0), the reviewer's exceedance and median ratios per shape, and P(pp).

(i) **Timing:** D-E0 to D-E3 and D-E5 to D-E6 are sandbox-only, run one at a time and not concurrently with any host job. D-E4 consists of 10 reference-engine calls, each ≤ 1500 s, one at a time, preceded by a 1e5-primary smoke case.

A later model change goes through a new decision and a pre-registered re-qualification at the declared fresh bases of the integration route (lv 20531004, hr 20541004; V6-F, if the operator selects it, at 20551004) with **unchanged** tolerances.

Supporting plans (planning documents, not evidence): `validation/plans/v3-005d-v7r-calibration-repair.md` (V7-R2 design; its pre-registration is a separate later amendment and its host execution waits for the operator's budget decision) and the integration route (committed with this amendment). Seed disjointness is asserted by `tests/ionmc/test_v3_005d_seed_families.py`.

*(a) Addendum to Amendment 19 (2026-10-14, after D-E0–D-E2 were run and independently reviewed; before D-E3).* Recorded, not changed: (1) the K* formula in §2.1 of the plan of record writes IDD_EM(z_pk − Δ); its prose ("an event with deficit Δ moves that proton's contribution … to z_pk − Δ" read against the physical statement that a range deficit moves the peak contribution to a shallower position of the *shifted* curve) is a sign transcription error. The frozen literal form stays computed and reported "as frozen"; no conclusion is drawn from it. The form IDD_EM(z_pk + Δ) is reported beside it as the physically meaningful kernel (D-E0 record: F predicted 0.853 / 0.832 against the archived 0.857 / 0.844). K* validates D-E3 only and never enters the §5 rule, so no identification threshold is affected. (2) The reviewer's round-2/3 statement that the archived V10-A first minima at 17.17° / 14.99° show "a core about 2× too wide" is refuted by D-E1: those minima are duplicate-angle artefacts (0.74σ / 0.48σ against their neighbours; the data fall monotonically there). The refuted prior is recorded with this date; the §4.2 M3 first-minimum trigger is evaluated against the plan's frozen minimum criterion only. V10-A's verdict at dd976c5a is unchanged (its 13.74° cut was set by that spurious minimum, which is recorded as a note on the row, not a change). (3) Any shape parameter fitted to V11 or to Geant4 in D-E0 (the "TOPAS-equivalent s / R×" diagnostics) is data-seen, quarantined in `validation/results/elastic_attribution/de0/quarantine-data-seen.json`, and is never cited or used as a model parameter. (4) D-E0–D-E2 are exploratory records; they change no row, tolerance, base or verdict.

*(b) Addendum to Amendment 19 (2026-10-14, after D-E3, the H5 check and the D-E5 table side were run and independently reviewed; before any D-E4 reference run and before D-E6).* Recorded, not changed, unless stated:
1. **D-E3 / H5 records** (`validation/results/elastic_attribution/de3/`, `de0h5/`; logged runs in the experiment telemetry): the base copy of the pinned table reproduces it bit for bit (array digests recorded); `all` reproduces the archived V11 F within 0.6 σ at the same 1e5-history size; H5 (MCS overlap) is **excluded** (EM-only peak/plateau ionmc/TOPAS within ±2 % at both energies). Findings preserved with the review: F_all is consistent with F_pA × F_pp within the paired-jackknife uncertainty of the contrast (value ± SE in the D-E3 summary; Codex REVIEW-a7c99b45 required the uncertainty), i.e. the composition is multiplicative to the precision of the probes; the shape variant `pA-s0.9` closes 84 % / 76 % of the archived V11 gap, `pA-s0.5` and `pA-Rx1.2` about half, `pA-gauss` none on F (contrary to a pure tail story: the peak deficit follows the core width, the plateau gain the tail); the pp-only plateau gain alone is of the order of TOPAS's total elastic gain, so **H3 stays live** and the `swap-pp` arm of D-E6 is essential; the shifted K* overpredicts peak removal by 8–14 % (p+A) and 1–3 % (p-p). The committed records come from a second execution of all 19 runs with the committed tooling (portable scratch cache); every run reproduced the first execution (untracked scripts, preserved in the experiment telemetry) bit for bit (depth-dose and recoil batches, counters, lateral rms batches). D-E3 never enters §5.
2. **D-E4 preconditions** (read-only verification of the engine, 2026-10-14): (i) OpenTOPAS has no filter on the *current* kinetic energy at the scoring surface (only initial-KE filters), so every exiting proton is recorded; the plan's expectation of "≤ 5e4 records" per case is wrong, the record count is ≈ the number of primaries, and the E > E₀ − 1e-4 MeV cut is applied in the analysis instead. (ii) The reference runner limits every output file, stdout included, to 512 MiB; stdout carries one line per history, so `i:Ts/ShowHistoryCountAtInterval` is set, and at 2e7 primaries a binary phase space (≈ 1 GB estimated) exceeds the limit. **Split rule, declared now and applied only after the smoke case measures bytes per record:** a case whose estimated phase-space size exceeds 400 MiB is split into k equal parts, k = ⌈size / 400 MiB⌉, with seeds 20270000 + 1000·row + j, j = 1 … k (this supersedes the "k = 2 only" wording of (c); the j = 1 seed is the one declared there); the parts are concatenated before analysis and treated as one case. (iii) Whether the module list `g4h-elastic_HP` alone is accepted is undecidable without a run; the §3.2 fallback of the plan of record applies if it is rejected. (iv) The smoke case `validation/reference_cases/topas/v3-005d-de4-smoke-water-150mev-1cm` (1e5 primaries, seed **20299001**, inside the 2027xxxx block but outside the row seeds 20284001–20293001 and unused elsewhere) is exploratory and report-only: its outputs are used for format, size and timing only, never for a signature. The 10 cases are frozen after it, with the counts of (c) unchanged.
3. **D-E5 classification rule** (review of the table side): the two kinematic curves lie within 1e-3 of each other up to θ_lab ≈ 1.8° (the plan assumed 0.5°), so 0.4–1.9 % of events are "conflict" tracks; excluding them would bias the Geant4/ionmc ratios shape-dependently. Rule, applied identically to the table side and to the Geant4 phase space: a conflict track is assigned to p+A unless its event contains a secondary proton (then p-p); the raw conflict fraction stays reported as a diagnostic. The first-minimum, Helm and q-slope descriptors of D-E5 are **recorded, not interpreted** (the frozen 1° estimator picks noise dips in a 256-bin piecewise-uniform table; the 3σ variant misses the J₁ zero).
4. **D-E6 interaction sensitivity (declared before D-E6 is run):** the §5 rule of record is unchanged (Δ on 1 − F, thresholds 0.010 / 2·SE / 70 % / 30 %). Because D-E3 showed F to be multiplicative, the same main effects and interaction are also reported on −ln F as a pre-declared sensitivity. If clause 3 of §5 (70 % / 30 %) holds under one metric and fails under the other, clause 3 is recorded as **"inconclusive (metric-dependent)"**, never as satisfied. No threshold is relaxed.
5. **Git policy reconciliation (Codex REVIEW-a7c99b45):** the committed exploratory records follow the explicit allowlist stated in `validation/results/elastic_attribution/README.md` (dimensionless ratios, per-primary probabilities and rates, kinematic statistics, kernel integrals, counts, flags, provenance; never absolute cross-section values, tabulated or sampled distributions, PDG/EXFOR points). The D-E3 scratch cache is portable (`validation/generated/elastic_attribution/de3-cache/`, git-ignored) and the 18 D-E3 runs are repeated with the committed tooling and compared bit for bit against the preserved first-run files.
6. Everything in this addendum is exploratory; no row, tolerance, base or verdict changes.

*(c) Addendum to Amendment 19 (2026-10-15, written AFTER the D-E4 reference runs and the D-E4/D-E5 comparison were seen and independently reviewed; before any D-E6 run).* This addendum changes the scope of D-E6 and says so; the frozen §5 rule of record and its thresholds are unchanged and are evaluated first.
1. **D-E4 record** (`validation/results/elastic_attribution/de4/`): 18 OpenTOPAS 4.3 / Geant4 11.4.2 runs (rows 14–23; 4e6 primaries at 100–200 MeV, 2e7 as three seed parts at 30/50 MeV per (b)2; module `g4h-elastic_HP` alone was accepted — proton model hElasticCHIPS with BarashenkovGlauberGribov cross sections; G4_C taken at the Geant4 NIST density 2.0 g/cm³, not printed in the log), parsed one case at a time with every part bound to the frozen case (histories, seeds, module list, geometry), phase spaces deleted after parsing (disk-retention rule), empirical μ_CM CDFs extracted for D-E6. Findings preserved with the review: **H1 (σ level) excluded** — p+A rates table/Geant4 1.005–1.031 at every energy, as expected from the shared Barashenkov σ_el, so D-E4 tests shape only; **H2 signature holds** — the table's p+A angular distribution is wider than CHIPS at every energy (θ_lab median ratio 1.07–1.61, exceedance ratios 1.3–3.4, tail moment 2.4–26×), nothing contradicts H2a; **H3 identified on rate** — the p-p rate ratio is ≈ 1.00 at 30–50 MeV, 1.07 at 100 MeV and ≈ 1.15–1.21 at 150/200 MeV (shape ratios ≈ 1.00), consistent with the D-E2 finding (table σ_pp = 1.10 × BGG at 150/200 MeV) and with the share of CHIPS p-p below the 16.26° cut; **H4 borderline** — the recoil-energy ratio ≈ 4–5 could carry ≈ 24 % of the gain excess at 150 MeV against the 20 % exclusion bound; evaluated in D-E6 through the recoil channel.
2. **p-p selection rule corrected (both sides):** the review found that "exactly two on-curve proton records" excluded 2.2–4.6 % of Geant4 p-p events that are p-p followed by a second (p+O) interaction in the target; the table side has no such events. Rule now: a p-p event is accepted when its fast proton lies on the p+H curve (the recoil record is used for the secondary flag only); the excluded fraction stays reported. The D-E4 record is regenerated under this rule before D-E6 (the 150 MeV p-p rate ratio is expected to drop to ≈ 1.16).
3. **D-E6 arms extended (post-hoc, declared here):** because H3 is identified on rate and `swap-pp` exchanges only the p-p shape, the frozen 2 × 2 (base, swap-pA, swap-pp, swap-both) cannot measure the p-p rate effect. Two arms are added: `swap-pp-sigma` (ionmc shapes; σ_pp scaled per energy by the Geant4/table p-p rate ratio after the corrected rule, interpolated linearly in E between the D-E4 energies, factor 1 below 30 MeV) and `swap-all` (both Geant4 shapes and the scaled σ_pp). Seeds, sizes and CRN as in (b)/(e). Evaluation order: (i) the frozen §5 rule on the frozen 2 × 2 is evaluated and recorded first, exactly as written; (ii) the extended set is reported beside it: Δ(all), the three main effects (pA shape, pp shape, pp σ), the interaction I_all = Δ(all) − ΣΔ(main), and the closure r_all = Δ_pk,obs − Δ(all) with SE_r per clause 2, on 1 − F and on −ln F ((b)4). If clause 2 fails on the frozen 2 × 2 but holds on the extended set, the outcome is recorded as **"closure reached only with the post-hoc p-p σ arm (declared after D-E4)"**, never as "established" under the frozen rule; thresholds (0.010, 2·SE, 70 %, 30 %) are unchanged.
4. **CDF edges and finite-sample error:** the end edges of the empirical equiprobable binning are the sample minimum and maximum of μ_CM (not −1/+1; the review found the −1 edge inflates the swapped ⟨1 − μ⟩ by ≈ 11.5 % at 150 MeV). The finite-sample CDF error of §3.2 is estimated as the quadrature sum of (a) a 200-resample bootstrap of the empirical CDF per energy propagated through the D-E0 kernels K* and G, (b) the tail-edge systematic (sample-extreme edge versus the next quantile), and (c) a leave-one-energy-out check of the q-scaled interpolation at 150 MeV (predicted from 100 and 200 MeV; applied to 200 MeV as a flagged proxy); it enters SE_r of clause 2. The bootstrap uses the per-energy seeds **20990200 + E**, i.e. 20990350 (150 MeV) and 20990400 (200 MeV) (base 20990200; investigation family 2099xxxx, unused elsewhere; the earlier wording "seed 20990200" named the base only).
5. Everything in this addendum is exploratory; no row, tolerance, base or verdict changes.

*(d) Addendum to Amendment 19 (2026-10-17, after D-E6 was run and independently reviewed).* **§5 outcome: "partially attributed"** — clauses 1–4 hold at both energies (frozen 2 × 2: r = +0.0014 / +0.0051, SE_r 0.0027 / 0.0037, threshold 0.010; H2 (p+A shape) carries 100 % / 99 % of Δ(both); D-E5 signature holds); clause 5 fails (H3 not excluded: table/Geant4 p-p rate ratio 1.165 / 1.159; H4 not excluded at 150 MeV: recoil part 25 % vs the 20 % bound, arising only in the shape-swap arms — contrary reading recorded). The −ln F sensitivity ((b)4) agrees with the rule of record on every clause. The post-hoc arms ((c)3) change no label (closure already holds on the frozen 2 × 2). Full wording in decision 0041 (dated section 2026-10-17); the independent review rounds D-E3–D-E6 are in `validation/plans/v3-005d-elastic-attribution/review.md`. **No model change is authorised by this record** (§5 "established" clause; the §4 M2 wording conflicts and is not read as an authorisation); the question goes to the operator with the batched V3-005D intervention. V5, V11, V10-A at dd976c5a stay failed; nothing in Amendment 19 changes a row, tolerance, base or verdict.

**Amendment 20 (2026-10-18; revision 7, independently reviewed in rounds A20 to A20-6 of `validation/plans/v3-005d-elastic-attribution/review.md`; committed together with decision 0042) — operator response to IR-20261017-223039-9DA564: model rule under partial attribution, M2 inputs and construction, pre-registered go/no-go diagnostic, Option II gate set, seeds and budget; outcome "no model change" (M2 failed §4.2 trigger 3; M3 not acquirable within the tool and the routes tried).** The revision tables below are the review history; earlier revisions are kept in the session scratchpad and are not of record.

**Revision 2 (2026-10-18)** answers review round A20 in `validation/plans/v3-005d-elastic-attribution/review.md` (REVISE, nothing blocking). Revision 1 is kept as `amendment-20-draft.rev1.bak`.

| Review item | Severity | Change in this revision |
|---|---|---|
| 2a | major | CONTINUE becomes three **necessary** point-estimate conditions computed with the exact frozen Welch-t rule (V11 F band 63.1–136.9 % / 73.6–126.4 %; V11 gain; V5 ΔIDD plateau integral). The arithmetic is shown. An honest prognosis is added: a Geant4-like shape fails the 200 MeV gain condition (item (c)) |
| 2b | important | The diagnostic is declared a go/no-go on the single pre-committed model. The lv5d V11 is not blind to it; the disclosure goes into the amendment and the diagnostic record header |
| 5 | important | hr5d uses r_index 1–3 only; every other hr r_index is reserved (lv r ↔ hr r − 10). The route's V6-F range misprint 20558015 → 20558019 becomes a required one-line correction of the route (revision 7) in the same commit |
| 6 | important | A measured-first V6-F CUDA timing step is added. V6-F runs first in the evidence run, and the failure branch is pre-declared: lv5d/hr5d are deferred, not run |
| minor | – | Rule (b) is worded as a post-outcome rule. The M2 table carries the Amendment 16 qualification flags, with P6-D and P7 re-evaluated on it. The σ_pp lineage is split: the level is independent of TOPAS, S(E) above 150 MeV is shared BGG |
| new | – | Section (b′) "M2 inputs": SAID SM16 record with the header verbatim and the convention uncertainty; **β_NN declared refinement** (the Horiuchi unitarity rule applied to per-channel SAID at every node, pp and pn phases summed; supersedes the plan wording); density substitution rule (≤ 0.03 fm rms cross-check, chosen after a preliminary estimate, disclosed; if no substitute passes, plan §4.2 trigger 1 → M3 as the plan says) |
| new (d) *(superseded by revision 3, item 1)* | – | **(b′)4: the p-p level is re-anchored to SAID SM16 σ_tot × f_cut.** The defect: the table's partial 26.3 mb exceeds SAID's total 25.72 mb at 150 MeV. Size ≈ −6 % / −9 % at 150 / 200 MeV. BGG is no longer used for p-p. The diagnostic tests M2 + (b′)4, with M2 alone as a report-only decomposition. V10 is re-evaluated; H3 stays open |
| new (e) | – | **(b′)5: SAID SM16 re-acquired for 5–300 MeV** (8 acquisition ids; consistent with the 10–300 tables at all 59 common energies), linear interpolation in T_lab, no hold-constant rule |

**Revision 7 (2026-10-18)** applies the text fixes of review round A20-5. Revision 6 is kept as `amendment-20-draft.rev6.bak`.

| Item | Change in this revision |
|---|---|
| 1 | The O-16 value is printed as 1.407 at 150 MeV (0.851 / 1.407 / 1.622), with its 0.5 % margin above the bound. The plan's band rationale is recorded as wrong at 150–200 MeV. C-12 trigger 4 at 150 MeV is "not evaluable (no minimum)", moot |
| 2 | "No model change" now cites plan §4 (M3 last fallback, M4 rejected) plus (b″)4, not (b)3. The same path is added to (c)6 |
| 3 | M3 is "not acquirable within the tool and the routes tried"; the operator-assisted route (raising the 256 MiB cap, or manual placement of MAT 625 / 825) is named |
| 4 | The closure proposal now carries the default-off V5 comparison (+12.5 / +12.3 % against −3.7 / −4.9 %) and states that V6-F certifies the elastic-on configuration only |

**Revision 6 (2026-10-18)** records the outcome of the M2/M3 chapter. Revision 5 is kept as `amendment-20-draft.rev5.bak`.

| Item | Change in this revision |
|---|---|
| M3 | Construction impossible: the ENDF/B-VIII.0 neutron files (C-12 MAT 625, O-16 MAT 825) could not be acquired. The 310 158 456-byte zip exceeds the 256 MiB tool cap, and every per-material route failed (verbatim attempt table in (b‴)1) |
| Planning defect | Recorded, author the orchestrator: the frozen row (ii′) did not verify the energy coverage of the neutron files; coverage to 135–200 MeV is unverified (open sources state only library-wide "at least 20 MeV") |
| Outcome | **"No model change"**: V5 / V11 / V10-A stay failed; the D-1 table remains the model; the M2 / M3 records and the `--model m2/m3` code are preserved; no production table is built; the diagnostic (c) is not run; (d) lapses |
| Integration | Gates still executable without a model change: V6-F, V7-R2, V8-nuc and the implementation gates. V5 / V11 cannot pass, and the elastic depth-dose claim is withdrawn. Merge needs a closure amendment, defined in (b‴)4, which **requires operator approval** because it reduces the operator-selected gate set. Budget ≈ 43.5 h |
| Successor | Amendment 21 (separate draft): candidate M2′ |

**Revision 5 (2026-10-18)** follows from building and independently verifying the pre-committed M2 (scratch; no transport). Revision 4 is kept as `amendment-20-draft.rev4.bak`.

| Item | Change in this revision |
|---|---|
| Trigger 3 | New (b″). σ_el^M2/σ_el^Barashenkov for O-16 is 0.851 / **1.407 / 1.622** at 100 / 150 / 200 MeV, so **trigger 3 fires**; C-12 is outside the band at 1 of 3. The trigger is applied as written. Verified cause: a free-NN optical-limit limitation (σ_tot +16–19 % against EXFOR 13569 / 13753; α_NN = 0 removes it), not wrong inputs. Disclosed: seen before commit |
| Model under test | **M3** (ENDF/B-VIII.0 neutron O-16 / C-12 MF4/MT2, isospin mirror) for O and C; other targets keep the black disk. Defined: conversion to 256 μ_CM bins, q-scaled hold above E_top, MF1 identity, copy check, normalisation, positivity, CDF checks. If no open host serves the files, the outcome is "no model change" |
| M2 | Kept as the report-only sensitivity V-M2 (D4 seeds 20991200 + case). No medium-corrected α_NN. Construction record kept: trigger-4 Helm definition (q_H C 1.912, O 1.661 fm⁻¹; O q_min 0.83 q_H; C not evaluable at 150 MeV); rms-matched Fermi unfolding; black-disk shape below the SAID range |
| Diagnostic | V-M3 is the model under test, V-M2 is report-only, V-c is withdrawn. Seeds, bands, gates and budget are unchanged |

**Revision 4 (2026-10-18)** answers review round A20-3 (REVISE). Revision 3 is kept as `amendment-20-draft.rev3.bak`.

| Review item | Severity | Change in this revision |
|---|---|---|
| 1 | blocking | σ_anchor uses the half-sphere integral 2π∫₀^0.96 (Amendment 14 addendum (4)); transcription correction |
| 2 | major | Mott exchange term: −cos(…); the builder reuses the P6-tested `coul` function |
| 3 | major | **Outcome recorded:** σ_anchor 27.26 / 24.50 / 23.22 / 22.47 mb against the table's 27.30 / 22.94 / 26.26 / 25.65 mb at 100 / 125 / 150 / 200 MeV. Continuity is +6.8 % at 125 MeV, above 3 %, so (b′)4 is **not applied**. The values were seen by the reviewer before commit, after the rule was fixed in revision 3; the rule is not changed. V-a is withdrawn; **V-b (M2 alone) is the model under test**; V-c stays. r(E) is reported, report-only. Units/frame ratios 0.86–0.93 resolve the DSG unit uncertainty |
| 4 | major | New H3 lead recorded: the table's Hale-NI partial is non-monotonic (22.3 mb at 130 MeV, 26.3 mb at 150 MeV); the S(E) step inherits it (table/SAID ≈ 1.13 at 150–200 MeV, matching D-E4). Separate pre-registered investigation named: "P6-D-type check of the Hale NI integration around 130–150 MeV" |

**Revision 3 (2026-10-18)** answers review round A20-2 (REVISE). Revision 2 is kept as `amendment-20-draft.rev2.bak`.

| Review item | Severity | Change in this revision |
|---|---|---|
| 1 | blocking | (b′)4 rewritten. The "partial > total" argument is withdrawn: the Hale partial includes interference, SAID SGT is nuclear-only. Re-anchoring is restricted to **E ≥ 150 MeV**, with a linear ramp of the factor from 125 MeV. The anchor is SAID pp DSG minus the analytic Mott term (formula stated), integrated over \|μ_CM\| ≤ 0.96 like the table. Pre-declared continuity check: more than 3 % at 100 or 125 MeV means no re-anchoring, and M2 alone becomes the model under test. Full disclosure: SAID shares data lineage with BGG and Hale, and the component was written after D-E2/D-E4/D-E6. Size not assumed. Pre-registered diagnostic variants V-a / V-b |
| 2 | major | V5 peak/plateau screen added from the frozen rule (step 21): closure **61.1–126.5 %** at 150 MeV and **80.7–120.2 %** at 200 MeV. The effective F window becomes 63.1–126.5 % / 80.7–120.2 % |
| 3 | major | β_NN and densities relabelled as **declared prospective revisions**, not trigger 1. The plan-literal β variant (Horiuchi Table I, constant below 30 MeV) is built and reported as a sensitivity (V-c). Densities: the plan-literal source is unavailable (paywall); a failed substitute means trigger 1 → M3 |
| 4 | minor | (b)5(ii)/(iv) reconciled: no re-qualification result closes H3/H4; that needs a new pre-registered test. P6 item (3) stays a reconstruction check |
| densities | – | (b′)3 rewritten with the per-target decisions. The 1974/1987 compilations are treated as one family. C-12 / O-16 HO; Si-28 2pF (Mehndiratta & Shukla), with the de-Vries-1987 Si-28 sets failing the check as a recorded anomaly; Ca-40 3pF by the fallback order. N-14 / Al-27 / P-31 use IAEA-radius-constrained surrogates (rms check void by construction; diffuseness sensitivity; not in V5/V11/V10-A). The density-parameter policy note is declared for review |

---

> **Amendment 20 (2026-10-18; written after the operator's response to IR-20261017-223039-9DA564 and after the §5 outcome of Amendment 19 (d); before any new transport, reference or calibration result; the M2 table was later built in scratch and failed §4.2 trigger 3, revision 5) — model decision for the recorded outcome "partially attributed", M2 inputs, the pre-registered M2 diagnostic (go/no-go), one fresh V5 base, the revised merge-gate set (Option II), the lv5d/hr5d suites and seeds, and the order of work.**
>
> **What this amendment does not change:**
> - no frozen row, frozen tolerance, recorded verdict or consumed base;
> - the text of Amendments 1–19 or their addenda;
> - the plan of record `validation/plans/v3-005d-elastic-attribution.md` (sha256 fd9657df…, manifest `validation/plans/v3-005d-sha256.txt`).
>
> **[new scope]** marks an added or replaced rule. **[unchanged]** marks a restatement.

## (a) Operator response (verbatim) and its reading

The response to IR-20261017-223039-9DA564 is recorded in `reference_requests/interventions/IR-20261017-223039-9DA564/response.json`:

"Scientific model selection, calibration design and fresh evaluation-base allocation remain autonomous responsibilities. You may revise internally authored plans prospectively, with independent review and preservation of original failures; my earlier instruction did not reserve those technical choices to me. For the integration boundary, I select Option II, including V6-F as a binding additional gate, together with the proposed V10-A/V11 and V7-R2 treatment. This does not waive any final experiment requirement or establish full halo validation. Perform the inexpensive candidate-model diagnostics before expensive qualification, retain the stability restrictions, and merge only after the revised applicable gates and independent review pass. Formally reconcile the intervention record with this response and resolve the outstanding review finding."

Reading:
1. **Autonomous.** Each of these is taken autonomously, through this independently reviewed amendment:
   - the model decision (b);
   - the diagnostic rule (c);
   - V7-R2 "as written" (`validation/plans/v3-005d-v7r-calibration-repair.md`, revision 4);
   - the bases (f);
   - the revision of Amendment 17 (b)4 for one V5 base (d).
2. **Operator-selected.** Option II of `validation/plans/v3-005d-integration-route.md` (§2.2), with 2B (V10-A report-only, V11 added) and 3A (V7-R2 in V7-R's slot).
3. **Not waived or established.** No final experiment requirement is waived. No full halo validation is established: V6-F is reduced scope (`validation/plans/v3-005c-v6f-followup.md` §10).

## (b) Model decision rule for the recorded outcome [new scope; post-outcome rule]

1. **Status of this rule.** It was written **after** the §5 outcome was known (Amendment 19 (d)) and is adopted under the operator's response above.
   - It authorises the M2 adoption **for this outcome**. It does not re-interpret §5 of the plan of record.
   - §5 authorises a model change only for "established". §4 says "M2 adopted if §5 attributes H2". This conflict, recorded in decision 0041 (2026-10-17), is settled for this case only.
   - The recorded outcome stays **"partially attributed"** and is never re-labelled.
2. **Rule.** Clauses 1–4 of §5 hold at both energies and the named component is **H2 (p+A shape)**, so the orchestrator adopts the pre-committed **M2** (optical-limit Glauber, plan §4.1, inputs per (b′)1–3, 5). The p-p level is unchanged: the re-anchoring of (b′)4 is not applied, by its pre-declared check. The hypotheses that clause 5 did not exclude become **open findings**.
   - *Revision 5:* M2 was built and failed §4.2 trigger 3 for O-16. By the plan's fallback rule, **M3 is the model under test** ((b″)). M2 is kept as a report-only sensitivity.
   - Facts (records under `validation/results/elastic_attribution/`): r = +0.0014 / +0.0051 against 0.010; H2 carries 100 % / 99 % of Δ(both) (`de6/de6-report.md`).
   - Open findings: H3 at both energies; H4 at 150 MeV.
3. **M3 triggers of plan §4.2** [unchanged]: construction impossible; numerical verification fails; σ_el level outside [0.7, 1.4]; first minimum outside [0.6, 1.0] × q_H. They are evaluated once on the M2 construction, before any transport with M2.
   - The β_NN rule (b′)2 and the density substitution (b′)3 are **declared prospective revisions** of plan §4.1, not trigger-1 cases. The plan-literal β construction could be built, and it is reported as a sensitivity.
   - If the density substitution of (b′)3 fails its cross-check, no qualified density source exists, because the plan-literal source is unavailable. That is trigger 1 (construction impossible), and the fallback M3 follows with its own pre-committed definition, inputs and §4.2 conditions.
   - Triggers 1–4 lead to M3 as written. M3 then passes the same pre-registered screen (c) before any host step.
4. **The new table** (M3 for O-16 / C-12 since revision 5; (b″)) is a new pinned elastic table.
   - It must carry the Amendment 16 qualification flags: `no_negative_density_in_domain`, `o16_mt5_sigma_not_c12_copy`, `o16_mt5_spectra_independent`, `shape_normalised`.
   - **P6-D and P7 are re-evaluated on it** (P6-D in CI; P7 in lv5d, r_index 15).
   - The pin change enters `src/` and the hashed source set **before the code freeze**. The rehearsal fixture of (f)5 therefore has to be regenerated after it.
5. **What "open finding" obliges** [new scope]:
   - **(i) Reporting.** Every lv5d V5 and V11 record, the V6-F record, the decision 0042 outcome and the results README carry an "Open findings H3 / H4" block next to the verdicts. Its contents:
     - the p-p rate ratio, ionmc table / Geant4 = 1.165 / 1.159 at 150 / 200 MeV (`de4/de4-ratios-report.md`, "ratio = ionmc table / Geant4");
     - the post-hoc σ_pp arm's share, 6 % / 13 % of Δ(all) (`de6/de6-report.md`);
     - the recoil share of the plateau-gain change.
   - **(ii) H3.** σ_pp is **not** corrected. The re-anchoring of (b′)4 is not applied (revision 4), so the p-p level is the table's own: Hale NI ≤ 150 MeV, S(E) rule above. The p-p rate ratio against Geant4 is reported. **New lead (revision 4), recorded under the open finding H3:**
     - The table's Hale-NI partial σ_pp is **non-monotonic**: a minimum of **22.3 mb at 130 MeV** and **26.3 mb at 150 MeV**.
     - The S(E) step above 150 MeV inherits that upturn, giving table/SAID ≈ 1.13 at 150–200 MeV (σ_anchor values above). That matches the D-E4 p-p excess of 1.165 / 1.159.
     - So the lead for H3 is an **anomaly in the table's own p-p level between 130 and 150 MeV** (builder or Hale-NI evaluation), not a missing SAID re-anchoring.
     - It stays an open finding, for a separate pre-registered investigation: **"P6-D-type check of the Hale NI integration around 130–150 MeV"**. It is reported in every re-qualification record.
     - Any σ_pp change needs its own decision and amendment.
   - **(iii) H4.** The `elastic_recoil` deposit (MeV per primary) and its share of the plateau gain are reported per V5/V11 run at both energies.
   - **(iv)** No open finding is a gate. **No re-qualification result closes an open finding:** V5 / V11 / V6-F records report H3 and H4, they do not resolve them. Closing H3 or H4 needs a new pre-registered test, for example a D-E4/D-E5-type thin-target comparison of the transported table.
6. **Decision 0042** records:
   - the M2 choice and lineage (SAID NN shared with the Hale p-p lineage; densities per (b′)3), the β_NN refinement of (b′)2, the SAID convention uncertainty of (b′)1, the p-p re-anchoring rule of (b′)4 and its outcome (not applied by its pre-declared check), the H3 lead (non-monotonic Hale-NI level at 130–150 MeV) and the SAID range of (b′)5;
   - the domain caveat below about 100 MeV (plan §4.3);
   - the §5 contradiction verbatim from decision 0041 (2026-10-17);
   - this rule;
   - the go/no-go disclosure of (c)4.

## (b′) M2 inputs (acquisitions of record; declared refinements) [new scope]

**(b′)1. σ_NN and α_NN: GWU SAID solution SM16** (Workman et al., Phys. Rev. C 94 (2016) 065203; pp to 3 GeV, np to 1.3 GeV).
- **Acquisitions:** ids **efe6d363f28c423582df8d68d2efebd0**, **6035179d45c046ad93445457ea48ee7d**, **23ca28089b974a32b3ccd40d68665ecc**, **761b99787a674a69a3478768b1cb9385** (pp and np; SGT and ALFA). Role `construction`, retrieved 2026-10-18. The SP07 fallback of plan §4.1 was not needed.
- **Quantities:**
  - SGT is the nuclear total cross section (mb).
  - ALFA is Re/Im of the spin-averaged nuclear forward amplitude, Coulomb removed.
- **Coverage.** All 59 rows, T_lab 10–300 MeV in 5 MeV steps, are present and physical:
  - pp SGT 349.6 → 24.17 mb, minimum 23.98 mb at 255–280 MeV;
  - np SGT 945.5 → 35.64 mb;
  - pp ALFA maximum 1.846 at 95 MeV, 1.623 at 150 MeV;
  - np ALFA changes sign between 10 and 15 MeV.
- **Output header, recorded verbatim:** `SM16 0-3.0 GEV PP=50276/25348 NP=22832/12938 RAA (147)  5/03/07`. Its date field "5/03/07" is stale relative to SM16 (2016) and is **unresolved**. It is recorded as is, together with the retrieval date.
- **Input uncertainty, declared.** SAID does not state its exact Coulomb-removal and sign conventions in its own words. The sign is consistent with Ray, Phys. Rev. C 20 (1979) 1857 and with Abu-Ibrahim et al., Phys. Rev. C 77 (2008) 034607. This is recorded in decision 0042 as an input uncertainty of M2.
- **Below 10 MeV.** M2 runs from E_min,shape (C 6.0, O 6.58 MeV; Amendment 15 (b)1). The SAID tables were re-acquired for 5–300 MeV; see (b′)5.

**(b′)2. β_NN: declared prospective revision of plan §4.1** (orchestrator decision, under review with this amendment; revision 3: not a trigger-1 case).
- **Finding.** The acquired Horiuchi et al. 2007 preprint (arXiv:nucl-th/0612029, v1, id **f56c546a11ae4dacafa35f6a614892e3**) does not supply the input that plan §4.1 assumed:
  - its Table I is isospin-averaged ("here we use the mean values"; no pp/pn columns);
  - it covers E = 30–1000 MeV;
  - for E < 300 MeV its β is **not an independent datum**: it is fixed by the paper's unitarity rule, Eq. (10), β = (1 + α²)·σ_tot/(16π), chosen so that σ_el = σ_tot.
  - The erratum Phys. Rev. C 76 (2007) 039903 and the published Table I could not be verified (arXiv v1 only).
- **Rule (supersedes plan-of-record §4.1 item 1 and item 2's β wording).** That wording was "Horiuchi 2007 β_NN … held constant below the lowest tabulated energy".
  - M2 applies **the same unitarity rule** to the acquired SAID per-channel σ_NN and α_NN, **at every grid node**, separately for pp and pn.
  - It **sums the separate pp and pn profile phases**: Abu-Ibrahim et al., Phys. Rev. C 77 (2008) 034607, Eqs. (4)–(5). This is exact in the optical limit with ρ_n ∝ ρ_p, the declared approximation of plan §4.1 item 4.
  - No hold-constant rule is needed for β below 30 MeV.
- **Check at 150 MeV**, from the SAID values above. With σ in fm² and 16π = 50.27:

  | Channel | Arithmetic | β_NN | Abu-Ibrahim 2008 (interpolated) |
  |---|---|---|---|
  | pp | (1 + 1.623²) × 2.572 / 50.27 | **0.186 fm²** | ≈ 0.18 fm² |
  | pn | (1 + 0.9435²) × 5.214 / 50.27 | **0.196 fm²** | ≈ 0.196 fm² |

- **Cross-check, report only.** Horiuchi's averaged Table I is reported as a cross-check at its tabulated energies (for example β = 0.195 fm² at 150 MeV, mean). Differences are recorded, not used.
- **Status (revision 3).** This is a **declared prospective revision**, not a §4.2 trigger-1 case: the plan-literal construction (Horiuchi's isospin-averaged Table I, held constant below its lowest energy of 30 MeV) **can be built**. Reason for the revision: the source itself fixes β for E < 300 MeV by the unitarity rule, so applying that rule to per-channel SAID inputs is the source's construction without isospin averaging, and without the constant extension below 30 MeV. **The plan-literal variant is built and reported alongside as a sensitivity** (diagnostic (c), report-only). It never enters CONTINUE/STOP.

**(b′)3. Densities: declared prospective revision of plan §4.1 item 3 (substitution rule).** The plan-literal source, de Vries–De Jager–De Vries 1987 (At. Data Nucl. Data Tables 36, 495), is not acquirable (paywall).
- **Lineage of the substitution (declared).** The de Jager–de Vries–de Vries 1974 compilation (At. Data Nucl. Data Tables 14, 479) and the 1987 compilation are **one electron-scattering compilation family**. Open, citable sources that reproduce parameters from either count as substitutes. Reason: the 1987 compilation updates and extends the 1974 one with the same functional forms. Lineage: electron scattering only, independent of the p+A elastic data and of TOPAS.
- **Not a source.** The id 8e31f0ce1bf6466ebe2a128a10549571 is a DOI stub and is never used as one.
- **Validation (tabulated targets).** The substitute's rms charge radius must agree with the independently acquired IAEA charge-radius compilation (CSV, id 8dd779140e52459f81f506442382b4a8) to **≤ 0.03 fm**.
  - Disclosure: this tolerance was chosen after a preliminary, unverified estimate was seen (O-16 Δ ≈ 0.019 fm).
  - Independent rationale: electron-scattering compilations scatter at the 0.01–0.03 fm level.
  - The tolerance is not changed after comparison.
- **Per-target decisions** (Δ = substitute rms − IAEA rms; values from the researcher's per-target results):

  | Target | Form | Source (second source) | Δ [fm] | Notes |
  |---|---|---|---|---|
  | C-12 | HO | Shukla 2003, Table I (Fourier–Bessel confirmation, Duda) | −0.007 | – |
  | O-16 | HO | Shukla 2003 (TGlauberMC) | +0.019 | – |
  | Si-28 | 2pF, c = 3.106 fm, d = 0.542 fm | Mehndiratta & Shukla, arXiv:1702.07459, Table 1, de Jager 1974 lineage (Shukla 2003 row, Δ +0.028) | +0.016 | **Recorded provenance anomaly:** the de-Vries-1987-attributed Si-28 3pF set (TGlauberMC) and the Fourier–Bessel set (Duda) **fail** the cross-check at −0.047 and −0.038 fm. They are not used |
  | Ca-40 | **3pF** | TGlauberMC (de Vries 1987 lineage; source line 2090) (Mehndiratta 3pF, Δ +0.003) | +0.004 | The plan's preferred 2pF is not tabulated, so the plan's fallback order (3pF next) applies |

- **Targets with no open tabulated entry (N-14, Al-27, P-31): declared IAEA-radius-constrained surrogate.**
  - **Construction:**
    - N-14: the preferred HO form, with α = the mean of the C-12 and O-16 HO α.
    - Al-27 and P-31: the preferred 2pF form, with the Si-28 2pF diffuseness d.
    - The size parameter (a, or c) is fixed so that the charge rms radius equals the IAEA value exactly: N-14 2.5582 fm, Al-27 3.0610 fm, P-31 3.1889 fm.
  - **The rms cross-check is void by construction** for these targets.
  - **Declared check instead: diffuseness sensitivity.** The shape parameter is varied between the two neighbours' values (N-14: the C-12 and O-16 α; Al-27 and P-31: the Si-28 d and that of the next tabulated neighbour). The M2 shape change at 150 and 200 MeV (⟨1 − μ_CM⟩ and the first-minimum position) is reported.
  - **Why the construction is needed.** Naive nearest-nucleus scaling, radius × (A/A_ref)^{1/3}, fails the 0.03 fm rule for N-14 (+0.035 / +0.042 from C-12 / O-16) and Al-27 (+0.055). It passes for P-31 (+0.013); the constrained construction is used there too, for uniformity.
  - **Scope.** These targets do not enter V5 or V11 (water: H, O) or V10-A (carbon). They enter only the non-water `materials.py` tissues and bone.
  - This supersedes, for densities, the surrogate radius scaling of decision 0041.
- **Policy note, declared for review.** The density parameters (a, α, c, z, w) are **published model parameters of open sources, not nuclear-data values**. They are stored in a provenance JSON of the builder, together with the acquisition ids: Shukla dcb08bbcc4aa446e8808e1abfa5600e5; Mehndiratta & Shukla a3d9b475e79948b8b4f2f5f2de7c9d82 (sha256 d4187fe9…7e3e; arXiv licence); TGlauberMC 0d4d56435c234f5892e6a76bf5998393 / cdb9f72849a84e53a2b44e1e79f0b85f; Duda bf78ff8d0fc941f39cc79046a3dcab8a; IAEA 8dd779140e52459f81f506442382b4a8. Decision 0041's data rule ("absolute values are not committed") is read as not applying to such published model parameters. This reading is declared for review.
- **Proton-size unfolding** as planned: F_pt(q) = F_ch(q)/exp(−q² r_p²/6), r_p = 0.8414 fm (CODATA 2018).
- **If a tabulated target has no substitute passing the cross-check,** M2 is "construction impossible" (plan §4.2 trigger 1), and the fallback M3 follows exactly as the plan of record says, passing the same screen (c). The radius-constrained surrogates cannot fail the rms check, so the sensitivity report is their check.

**(b′)4. p-p elastic level at E ≥ 150 MeV re-anchored to a SAID interference-inclusive equivalent (second model component; orchestrator decision; declared prospective revision, under review with this amendment)** [new scope; revision 3; **not applied** by its pre-declared check, revision 4]
- **Withdrawn (revision 3).** The revision-2 argument "the table's partial (26.3 mb at 150 MeV) exceeds SAID's total (25.72 mb), which is unphysical" is **withdrawn**. The comparison is not like-for-like:
  - the table's transported NI density includes the Coulomb–nuclear interference (Amendment 14 (b): NI = σ_e − σ_c);
  - SAID SGT is the nuclear-only total.
- **Reason for the component, stated within its domain.** Only at **E ≥ 150 MeV** is the table level set by the builder's S(E) rule ("Hale NI at 150 MeV × BGG ratio", Amendment 14 (b), Amendment 15 (c)1). There the records show:
  - the level at 1.10 × the BGG formula (D-E2, `de2/de2-report.md`);
  - a p-p rate excess over Geant4 of 1.165 / 1.159 at 150 / 200 MeV (D-E4, `de4/de4-ratios-report.md`).

  At ≤ 100 MeV, table / Geant4 is 0.977 / 0.997 / 1.032 at 30 / 50 / 100 MeV, and **the table level is kept**.
- **Anchor quantity.** At each grid node E ≥ 150 MeV:

  σ_anchor(E) = 2π ∫_{0}^{0.96} [ DSG_SAID(E, μ) − σ_Mott(E, μ) ] dμ_CM

  This is the half-sphere integral, each event counted once (Amendment 14 addendum (4)). *(Revision 4, transcription correction: revision 3 wrote a both-hemisphere integral, which doubles it.)*

  - **DSG_SAID** is the SAID SM16 pp differential cross section as SAID outputs it: nuclear + Coulomb + interference, CM frame, mb/sr.
  - **σ_Mott** is the pure Coulomb cross section of two identical spin-½ particles. It is the Coulomb term of ENDF-102 eq. 6.9, the same function the LAW=5 reconstruction uses (NJOY `coul` convention, P6 item (2)):

    σ_Mott(μ) = (η/(2k))² [ 1/sin⁴(θ/2) + 1/cos⁴(θ/2) − cos(η ln tan²(θ/2)) / (sin²(θ/2) cos²(θ/2)) ]

    *(Revision 4: the exchange term of spin-½ identical particles in ENDF-102 eq. 6.9 is −cos(…). Revision 3 printed +cos. The **builder reuses the P6-tested `coul` function** of the existing LAW=5 code; it does not re-implement the formula.)*

    with η = α/β, where β is the relative velocity of the two protons (relativistic, as the reconstruction defines it), k = p_CM/ħ, θ = θ_CM, and the identical-particle spin-½ interference coefficient −1/2.
  - The SAID form has no Coulomb-only or nuclear-only option, so the Mott term is computed analytically.
  - The integration functional, the half sphere 0 ≤ μ_CM ≤ 0.96 with factor 2π and each identical-particle event counted once, is the one the builder applies to the Hale NI density (the P6 item (3) convention; Amendment 14 addendum (4)). σ_anchor and the table partial are therefore like-for-like interference-inclusive quantities.
- **Inputs: SAID SM16 pp DSG angle scans** at 100, 125, 150, 175, 200, 250 and 300 MeV, role `construction`, acquired.
  - URL template: `https://gwdac.phys.gwu.edu/cgi-bin/go3nn?sl=SM16&rt=1&ot=DSG&iv=A&il=<lo>&ii=1&iu=<hi>&fv=T&fn=<T>&jpeg=n&u=190&l=210`.
  - The CGI caps 60 points per request, so each energy is three 1°-step segments (θ_CM 1–60, 61–120, 121–179), 179 rows in all.
  - Header on all: `SM16 0-3.0 GEV PP=50276/25348 NP=22832/12938 RAA (147)  5/03/07`.
  - The output is the full DSG, including Coulomb and interference.

  | T_lab [MeV] | θ 1–60° | θ 61–120° | θ 121–179° |
  |---|---|---|---|
  | 100 | e8973bddcd3e4c1d8142033cf99d883c | cb0dfcecb64e4aaaa5829a03ac349650 | 0ad01128a2404b37be349aba6a7ca269 |
  | 125 | 7639d5dc7c2f4baf9491cee446e1574f | a6ea1a216a334db3a94edd7a45dc52ac | 12797d52222349d997efffd0768eb96c |
  | 150 | 9b95bac9bc1247dab43156074583a5c1 | 3eadb73fcf07484eb392afd8d4ba19c2 | f406ab8af2aa4c4abd11238219dcb4ee |
  | 175 | a3ece8cd3cfc42d6bb0aaa01193afe2d | 8909ef36e88f49578bfe5ca5f5c60968 | 14df69fa85c846a2999238d945ea7bca |
  | 200 | c2c3d72963f24e518630bf33f870ae7b | e1b69d9c30a3490e9783a97c380c84c6 | 2ed37a4234d548099155ce9e9a4f8a48 |
  | 250 | 0952913571f14abf84e492b4dbd724e2 | f7de116fbc3543d481200aa2ced8de21 | dbc89fb0adf248faa73c47790c82d83c |
  | 300 | 537907c89e034dc7a47804fa9d4710e1 | e38519d84ca74a78ab317eac3afd74ba | 411c5c0ae6544649aec4263669473d87 |

  - **Input uncertainty: units.** The output does not state its units. They are assumed to be mb/sr in the CM frame. Verification rule, applied by the builder before σ_anchor is used, at every scan energy: the ratio σ_anchor(E)/SGT_SAID(E) must lie in **[0.80, 1.05]**. σ_anchor is the half-sphere interference-inclusive partial, so it is expected below the nuclear total. A unit error (μb, fm², b) would move the ratio by a factor of 10 or more, and a frame error (lab instead of CM) by a factor of order 2–4. The bound is a planner bound for the unit and frame check, not an accuracy claim, and it is fixed before any scan value is examined. A failure at any scan energy blocks (b′)4 entirely, and the diagnostic runs M2 alone (V-b).
  - **Input uncertainty: interpolation in angle.** The 1° sampling must be interpolated onto the builder's μ grid. Rule: linear interpolation in θ_CM of ln(DSG − σ_Mott) where the difference is positive. If the difference is not positive, the scan cannot be used on that interval (σ_anchor is not formed at that energy), and this is recorded.
- **Rule.**
  - At nodes ≥ 150 MeV: σ_pp,table(θ_CM > 16.26°; E) = σ_anchor(E), interpolated linearly in T_lab between the scan energies. The Hale angular shape (the fixed 150 MeV shape, Amendment 14 (b)) is unchanged; only the level changes.
  - Between 125 and 150 MeV the multiplicative factor r(E) = σ_anchor/σ_table is interpolated linearly from 1 at 125 MeV to r(150), so that no step appears. Below 125 MeV, r = 1.
- **Continuity check, pre-declared; the scans at 100 and 125 MeV are report-only for the level.** σ_anchor is also computed at 100 and 125 MeV. **If |σ_anchor/σ_table − 1| > 3 % at 100 or at 125 MeV, the re-anchoring is not applied at all.** The model under test of (c) is then M2 alone, and this is recorded. The 3 % bound is fixed now, before any DSG scan is seen.
- **Outcome of the check (revision 4; recorded, the rule is not changed).**
  - **Who computed it, and when.** The independent reviewer computed σ_anchor (corrected half-sphere convention, corrected Mott sign) from the acquired SAID scans while verifying revision 3, i.e. **before this amendment is committed.** These values were therefore seen before commit.
  - **Rule fixed first.** The rule, the 3 % bound and the ramp were fixed in revision 3, before that computation. They are **not** retuned. The outcome is recorded.

  | T_lab | σ_anchor | table | continuity |
  |---|---|---|---|
  | 100 MeV | 27.26 mb | 27.30 mb | −0.1 % |
  | 125 MeV | 24.50 mb | 22.94 mb | **+6.8 % > 3 %** |
  | 150 MeV | 23.22 mb | 26.26 mb | – |
  | 200 MeV | 22.47 mb | 25.65 mb | – |

  - **By the pre-declared rule, (b′)4 is NOT applied.** The p-p level of the transported table stays as it is: Hale NI ≤ 150 MeV, S(E) rule above.
  - The builder still computes r(E) = σ_anchor/σ_table at every scan energy and reports it, **report-only**.
  - **The DSG unit uncertainty is resolved.** The half-sphere units/frame ratios σ_anchor/SGT_SAID lie at 0.86–0.93, inside [0.80, 1.05], which confirms the mb/sr CM assumption.
- **Expected sign and size.** **Not assumed.** They are stated in the record from the acquired numbers once available. The revision-2 estimates (≈ −6 % / −9 %) are withdrawn.
- **Disclosure.**
  - SAID shares **experimental-data lineage** with the BGG formula (fitted to PDG-type NN data) and with Hale's fit (pre-1998 pp data). It does not share lineage with Geant4 models.
  - This component was written **after** D-E2, D-E4 and D-E6 were seen.
  - It is not a Geant4-scaled level: the D-E6 post-hoc arm stays rejected as a model.
- **Lineage, had the change been applied** (it is not; revision 4):
  - The p-p level at E ≥ 150 MeV comes from SAID; the BGG p-p formula is no longer used there.
  - Below 125 MeV the level is Hale NI as before.
  - The p+A σ_el level stays Barashenkov / BGG, shared with TOPAS.
  - V10 is re-evaluated, report-only, with its stronger shared lineage stated.
- **Consequences** (the rule as written; not applied, per revision 4).
  - Part of decision 0042 and of the diagnostic (c).
  - P6-D, P7 and the Amendment 16 flags are re-run on the new table. **P6 item (3) stays a reconstruction check** of the LA150 values, because the re-anchoring is applied after the reconstruction.
  - H3 per (b)5(ii)/(iv).

**(b′)5. SAID energy range: re-acquired 5–300 MeV** [new scope]
- **Re-acquisition succeeded.** Role `construction`. The header on all tables is `SM16 0-3.0 GEV PP=50276/25348 NP=22832/12938 RAA (147)  5/03/07`. The tables are byte-reproducible.

  | Range and step | pp SGT | np SGT | pp ALFA | np ALFA |
  |---|---|---|---|---|
  | T_lab 5–300 MeV, step 5 | **28e78745bfc74da98de99c76035f87ea** | **d4360764abe742c98e488b42c0c7b90a** | **c39d4c94600e455c98ed374dd0c17114** | **dab63152cbf64d8c97d0048c20987777** |
  | T_lab 5–30 MeV, step 1 | **932daa93b61b416a97455f92f44abbd7** | **07803ab6aab2449ca940480325fa3655** | **b0206879d01941b2b0cb4d4d5a30fd34** | **84a77e24e0444f47987b700899e684cb** |

- **Consistency.** The 5–300 tables equal the earlier 10–300 tables ((b′)1) at all 59 common energies. The step-1 tables equal the step-5 tables at 5, 10, …, 30 MeV.
- **First rows at 5 MeV:** pp SGT 690.6 mb; np SGT 1635 mb; pp ALFA 0.7241; np ALFA −0.2913.
- **Use.** M2 and (b′)4 use the 5–300 step-5 tables, with the step-1 tables for 5–30 MeV. Interpolation is linear in T_lab, with **no hold-constant rule**: every M2 node (≥ E_min,shape = 6.0 MeV for C) lies inside the acquired range. The earlier 10–300 ids of (b′)1 stay on record.
- The ALFA sign and Coulomb convention remains the input uncertainty recorded in (b′)1.
- Below E_min,pp, σ_pp = 0 as before. (b′)4 uses SAID only at and above 12.53 MeV.

## (b″) M2 construction outcome: §4.2 trigger 3 fires; the fallback M3 is the model under test [new scope; revision 5]

**(b″)1. Facts.** The M2 builder (`src/ionmc/nuclear/elastic_glauber.py`, scratch build) was reviewed line by line: no physics discrepancy was found. An independent implementation (numeric 3-D form factor, J₀ by its integral representation, Fourier-space folding) reproduces it.
- **σ_el^M2 / σ_el^Barashenkov.** These are the committed P6 item (5) BGG values for C and O.

  | Target | 100 MeV | 150 MeV | 200 MeV | Unfolding |
  |---|---|---|---|---|
  | C-12 | 0.906 | 1.336 | **1.454** | HO, exact |
  | O-16 | 0.851 | **1.407** | **1.622** | HO, exact |
  | Si-28 | 0.78 | 1.34 | 1.71 | rms-matched |
  | Ca-40 | 0.79 | 1.26 | 1.58 | rms-matched |

- **Trigger 3** (band [0.7, 1.4] at ≥ 2 of 100 / 150 / 200 MeV, for O or for C) **fires for O** (1.407 at 150 MeV and 1.622 at 200 MeV). C is outside the band at 1 of 3 energies (1.454 at 200 MeV), so it alone would not trigger.
  - **Margin, stated:** the decisive 150 MeV value **1.407 is only 0.5 % above the bound 1.4**. The outcome hinges on this margin. The trigger stands as evaluated on M2 as defined, and no variant may be used to move it.
  - **Planning defect, recorded:** plan §4.2's rationale for the band (optical-limit Glauber "is expected to reproduce p+A σ_el at this energy to well within it") was wrong at 150–200 MeV.
- **Cause, verified.** The trigger text says a deviation of this size "indicates wrong inputs". That is **not** the cause here: the reaction cross sections agree with Abu-Ibrahim et al. 2008, Table III, within 1 %. The cause is a real limitation of free-NN optical-limit Glauber at 150–200 MeV:
  - Glauber σ_tot overpredicts the measured n+natC and n+natO σ_tot by +18 / +19 % at 150 MeV and +16 / +18 % at 200 MeV (Finlay 1993, EXFOR 13569; Abfalterer 2001, EXFOR 13753; both acquired, Amendment 14 (l)).
  - Barashenkov matches those data within 3 %.
  - With α_NN = 0 the excess largely disappears. It comes from the real part of the free NN phase; medium effects are not modelled.
- **Disclosure.** The trigger outcome was seen before this amendment is committed, from the scratch build and its verification. No transport with M2 or M3 has been run.

**(b″)2. Decision: the frozen trigger is applied as written.** M2 fails §4.2 trigger 3, so the plan's fallback **M3** is the model under test.
- The verified cause (a model limitation, not wrong inputs) is recorded next to the outcome. It does not change the outcome.
- A medium-corrected α_NN is **not** introduced: it would be a new model choice after the outcome was seen.
- **M2 is kept as a report-only sensitivity** (V-M2), never a gate, with the recorded limitation. Its σ_el level is Barashenkov, as for every table; the shape only is M2's.

**(b″)3. M2 construction record, kept for the report-only sensitivity.**
- **Trigger 4 for HO densities, definition recorded.**
  - The Helm parameters come from a q⁴-cumulant match: R_H² = √350·c₁ (invariant under Gaussian unfolding) and s_pt² = s_ch² − r_p²/3. That gives q_H = 1.912 fm⁻¹ for C-12 and 1.661 fm⁻¹ for O-16.
  - O-16 at 150 MeV: q_min = 1.380 fm⁻¹ = 0.83 q_H. Passes.
  - C-12 has **no true minimum at 150 MeV**, because the real part fills it. Trigger 4 there is recorded as **"not evaluable (no minimum)"**: under the frozen wording ("must lie within") that is not a pass. It is **moot**, because trigger 3 already decided. At 200 MeV, q_min = 1.596 fm⁻¹ = 0.83 q_H.
- **Proton-size unfolding for Fermi forms.** Where the plan-literal unfolding F_ch/exp(−q² r_p²/6) diverges (β < r_p²/3), the declared rule is the rms-matched point Fermi density (Shukla 2003, p. 7 prescription).
  - Validated at 100 MeV against the exact unfolding: σ_el −1.6 % / −1.0 % and q_min +0.5 % / 0.0 % for Si / Ca.
- **Nodes below the lowest SAID energy** keep the pinned black-disk shape. The counts are recorded.
- **β of the plan-literal variant.** The independent review notes that `HoriuchiBeta` (builder lines 266–269) applies the isospin-mean β to both channels with the per-channel σ and α. This concerned only the plan-literal sensitivity V-c, which is withdrawn (see (c)1).

**(b″)4. M3 construction** (plan of record §4, M3 = (ii′); Amendment 14 (e)9) [the plan's definition, made executable].
- **Targets:** O-16 and C-12 only (plan §4 table; (e)9: "replaces the Bessel shape" for p+O and p+C).
  - **Source:** the ENDF/B-VIII.0 **neutron** sublibrary, MAT for n-O16 and n-C12, MF4/MT2 elastic angular distributions, used as the isospin-mirror shape for p+O and p+C at the same T_lab.
  - **Other targets** (N-14, Al-27, Si-28, P-31, Ca-40): the plan of record defines M3 for O and C only, so they **keep the pinned black-disk shape** (decision 0041, D-1). N, Al, Si, P and Ca do not enter V5/V11 (water) or V10-A (carbon).
- **Conversion.**
  - MF4 may be given as Legendre coefficients (LTT = 1), tabulated (LTT = 2) or mixed (LTT = 3), and is read as given. Frame: CM; a lab-frame representation (LCT = 1) is converted to CM.
  - The shape is evaluated at the pinned table energy grid with the ENDF interpolation law. It is converted to **256 equiprobable μ_CM bins per node**: the same table schema and loader, and a new table id.
  - The domains [E_min,shape, 250] MeV are unchanged (Amendment 15 (b)).
  - **Above the highest MF4 energy E_top of a file**, the shape at E_top is held **fixed in momentum transfer q** (q-scaled to each node's p_CM), as D-E6 interpolated the Geant4 CDFs (Amendment 19 (c)). Declared now; E_top is recorded per file.
  - **σ_el stays Barashenkov**: only the shape is M3's. Coulomb stays in the MCS.
- **M3 verification** (§4.2-type, before any transport; each failure is a recorded construction failure):
  1. **MF1 identity.** MAT, NLIB = ENDF/B, version VIII.0, neutron sublibrary, as recorded at acquisition.
  2. **Copy check** ((e)9: "n-O16 above 20 MeV is not itself a C-12 copy"), the rule of Amendment 15 (d)1 transposed:
     - Compare the MF4/MT2 angular distributions of O-16 and C-12 at every shared incident energy above 20 MeV, evaluated on 181 equally spaced μ_CM points.
     - **Not a copy** iff the relative difference exceeds 1e-3 at > 50 % of the points at > 50 % of the shared energies.
  3. **Normalisation** ≤ 1e-6 per node, after conversion.
  4. **Positivity.** A non-negative angular distribution at every node and node midpoint (P6-D style). A negative value is a BuildError.
  5. **Monotone CDF and finite first moment** at every node.
  6. **Energy coverage.** E_top ≥ 150 MeV, or the q-scaling range is recorded as a limitation.
- **Lineage.**
  - ENDF/B-VIII.0 neutron evaluations of O-16 / C-12 are related to the proton channel by isospin symmetry (an isospin-mirror approximation, of the order of 5–10 %, plan §4).
  - They are independent of TOPAS's proton elastic model (CHIPS). Geant4's neutron high-precision data (below 20 MeV) are not used by TOPAS for protons.
  - **M3's V10-A results are exploratory**, because the data have been seen ((e)9).
- **Acquisition.** In progress. Ids are filled in when they exist.
  - **Pre-declared:** if no open host serves the ENDF/B-VIII.0 neutron O-16 and C-12 files, M3 is "construction impossible". The plan of record has no further fallback, so **the program closes as "no model change"**: V5 and V11 are not re-run, and the outcome is recorded.
  - If M3 is built and passes its verification, it passes the screen (c) as the model under test.

## (b‴) Outcome of the M2/M3 chapter: **no model change** [new scope; revision 6]

**(b‴)1. M3 is construction-impossible (recorded; grounded on the acquisition failure alone).** The ENDF/B-VIII.0 neutron evaluations were identified on the NNDC list (https://www.nndc.bnl.gov/endf-b8.0/lists/neutrons.list):
- **C-12: MAT 625** (LANL/ORNL, EVAL-AUG15, Hale, Young, Fu);
- **O-16: MAT 825** (LANL, EVAL-DEC16, Hale, Paris, Young, Chadwick).

They could not be acquired through the experiment's acquisition tool. Attempts, verbatim:

| Route | Result |
|---|---|
| NNDC neutron sublibrary zip, 310 158 456 bytes (HEAD 200), tried twice | tool error: exceeds MAX_BYTES 256 MiB (`infrastructure/experiment_v3/data.py`) |
| IAEA `n_0825_8-O-16.zip` | tool error; curl: 403 (bot challenge) |
| LANL T2 `…/ENDFB-VIII.0-neutron/O/16` and `…/C/12` (guessed URLs) | tool error; curl: SSL handshake failure |
| NNDC Sigma | offers libraries ≤ ENDF/B-VII.1 only |
| NNDC endf-releases | whole-library zips only |
| OECD-NEA JANIS | session-based (HTTP 500 / 404) |
| NNDC per-file URL guesses | 404 |

By plan §4.2 trigger 1 and (b″)4, M3 is **"construction impossible"**, in the precise sense of **not acquirable within the experiment's acquisition tool and the routes tried.** One route remains: **operator-assisted acquisition**, either raising the tool's 256 MiB cap or placing the MAT 625 / 825 files manually. Pursuing it would be an operator or infrastructure matter. It is named here so that the closure is not read as a physics judgment on M3.

**(b‴)2. Planning defect of the frozen plan, recorded (author: the orchestrator).** The angular shape matters at 135–200 MeV, the V10-A data range and the V5/V11 energies. **Whether the ENDF/B-VIII.0 neutron C-12 / O-16 evaluations extend to 135–200 MeV is unverified.** Open sources state only library-wide coverage:
- "at least 20 MeV" for most evaluations (IAEA-NDS-227, POINT2018 documentation);
- "0 < E < 20 MeV (higher for some)" (Paris & Hale, CSEWG 2018 slides).

**The frozen fallback row (ii′) did not verify energy coverage** (Amendment 14 (e)9; plan §4). Had the files covered only a lower range, (b″)4's q-scaled hold above E_top would have transported a low-energy shape at 150–200 MeV. This omission is recorded as an orchestrator planning defect. It changes no verdict, and it is not the ground for "construction impossible"; that ground is the acquisition failure.

**(b‴)3. Outcome: "no model change".** This follows **plan-of-record §4** (M3 is the last fallback; M4 is rejected) and the rule pre-declared in **(b″)4** (revision 5, written before the M3 acquisition outcome was known). The A20-4 text of (b)3 did not address this case.
- **V5, V11 and V10-A at dd976c5a stay failed.** The D-1 black-disk table of decision 0041 (`elastic-proton-1fd24cff…`) **remains the transport model.** The p-p level is unchanged as well ((b′)4 not applied).
- **Records preserved:**
  - the M2 construction record ((b″)1–3: trigger 3 fired on verified physics; the trigger-4 record; the rms-matched unfolding);
  - the M3 definition and its failure ((b″)4, (b‴)1–2);
  - the builder code paths `--model m2` and `--model m3`, kept as records of the attempt. **No production table is built** from either.
- **The diagnostic of (c) is not run**, because there is no model under test. (d) lapses unused: no V5 base is taken, and lv 20531004 r_index 4–7 and 14 stay unobserved.
- A second-generation candidate, if pursued, needs its own prospective amendment (draft: Amendment 21, M2′).

**(b‴)4. Integration in the "no model change" state.**
- **(a)6 gates that remain executable without a model change** (on the existing D-1 table, at a new evidence SHA):
  - **V6-F**, binding, executed once at base 20551004. Its observable includes the elastic channel as transported by the D-1 table.
  - **V7-R2** in V7-R's slot.
  - **V8-nuc**, in the hr5d V8 statistical parity.
  - **The implementation gates:** P6-D, P7, P5-ext, C1-ext, R1-D2, V8-LV, CI, on the pinned D-1 table.
- **(a)6 gates that cannot pass:** **V5 and V11.** No model change means no fresh V5 base ((d) lapses), and V11 at the D-1 table is the recorded dd976c5a failure.
- **Consequence.** Under the gate set of (e) as selected by the operator (Option II), **the task cannot merge in this state.** **The elastic depth-dose claim is withdrawn** (route §6 S1/S2): the elastic channel has no validated depth-dose claim.
- **Closure amendment, defined here; it requires operator approval before it is committed.** It would remove V5 and V11 from the operator-selected gate set. That is a change of required validation scope, which disposition A of IR-20261013-164616-DEFAF4 reserves to the operator. Its content:
  1. **V5 and V11** carried as recorded failures (dd976c5a; lv5b V5), with the attribution record (Amendment 19 (d)) and this chapter's outcome next to them.
  2. **The elastic channel** ships as an explicitly flagged physics option with a declared limitation: "elastic depth-dose effect over-corrected against TOPAS by 3.7–4.9 % in peak/plateau; no validated depth-dose claim". The proposal must choose between **default off** and default on with the limitation reported in every result's capability block (V2-CAP). The planner recommends default off: the non-elastic V3-005A configuration then remains the default physics.
     - **Both options carry a recorded V5 failure.** The comparison the operator needs:

       | Configuration | Archive | V5 peak/plateau vs TOPAS (150 / 200 MeV) | ΔIDD (150 / 200 MeV) |
       |---|---|---|---|
       | Default off (non-elastic physics (slice B, lv5b 7aae5bb)) | lv5b 7aae5bb (plan archive table) | **+12.5 % / +12.3 %** | −18.5 % / −14.7 % |
       | Elastic on (D-1 table) | lv5c dd976c5a | **−3.7 % / −4.9 %** | +10.3 % / +6.7 % |

       Elastic-on is the closer configuration. The default-off recommendation rests on not shipping an unvalidated model by default, not on a better depth dose, and it must be presented with this comparison.
     - **Configuration certified by V6-F.** The binding V6-F transports the elastic channel through the D-1 table (elastic **on**). Its reduced-scope halo claim therefore certifies **the elastic-on configuration only**. If elastic ships default-off, the default physics has no halo gate, and the V6-F claim refers to a non-default option with a known depth-dose deficiency.
     - Claim wording: "V6-F (reduced scope) passed for the elastic-on configuration (D-1 table); not established for the default (elastic-off) configuration"
  3. **The gates that remain** (list above) must pass at one evidence SHA, under the operational merge condition of (e).
  4. **Claims:**
     - no elastic depth-dose claim;
     - no lateral-halo claim beyond V6-F's reduced scope;
     - no measured-dose claim from V6;
     - V10-A report-only.
  5. **Ledger rows** (S-PHYS-PROTON-NUCLEAR, V1-MUST-005/013) state the withdrawn claim.
  6. **Budget:**

     | Item | Host hours |
     |---|---|
     | V6-F | ≈ 6.15 CPU + CUDA (≈ 0.16 planned) |
     | V7-R2 block (includes the fixture) | ≈ 33.9 |
     | lv5d / hr5d rows without V5 / V11 | ≤ 3.28 (route §4.2 remainder, an upper bound) |
     | **Total** | **≈ 43.5** |

     The V5 / V11 steps are not run.

## (c) Pre-registered diagnostic of the model under test: a go/no-go screen [new scope; **not run** in the "no model change" state of (b‴); retained for any later candidate, see Amendment 21]

It runs once, after the table of the model under test (M3, revision 5) is built and has passed its verification ((b″)4), and before any host step. It is report-only and never evidence. It decides only whether the model under test is used at all: **CONTINUE or STOP.**

1. **Quantities.** Pre-registered variants:
   - **(V-a) withdrawn** (revision 4).
   - **(V-M3) M3 for O-16 and C-12** (pinned black disk for the other targets), p-p level unchanged: **the model under test** (revision 5). CONTINUE/STOP applies to it. D4 seeds 20991100 + case.
   - **(V-M2) M2, p-p level unchanged: report-only sensitivity** (revision 5; formerly V-b). Never a gate. The trigger-3 limitation is recorded. D4 seeds **20991200 + case** (new, 2099xxxx family, asserted unused by the seed-family test).
   - **(V-c) withdrawn** (revision 5): it was a sensitivity of M2, which is now itself report-only.

   In D1–D7, "M2" reads as "the model under test" (V-M3); the same quantities are reported for V-M2. All runs use the D-E3 harness: warp-cpu f64, V5 geometry r = 20 cm, 1e5 histories, 20 batches, the table through the exploratory-table route.
   - Common random numbers with the archived D-E3 base, seeds **20990100 / 20990101** (the declared D-E3 CRN seeds, reused deliberately for CRN).
   - Paired-jackknife SE throughout.

   | # | Quantity | Definition |
   |---|---|---|
   | D1 | F closure | c_E = (F_M2 − F_base) / Δ_obs,E, with Δ_obs = +0.0348 / +0.0433 (`de6/de6-report.md`; archived V11 values in `de0h5/h5.json`, bound to step 16) |
   | D1b | predicted V5 peak/plateau relative difference | q_E = PP_ionmc,arch × (F_M2/F_base) / PP_TOPAS,arch − 1, with PP_ionmc,arch = 4.2406 / 3.8286 and PP_TOPAS,arch = 4.4019 / 4.0255 (step 21, lv5c archive) and F_base = 0.8580 / 0.8455 (`de6/de6-report.md`). Assumes the V5 "on" peak/plateau scales with the elastic factor, the multiplicative composition shown in D-E3 (Amendment 19 (b)1) |
   | D2 | predicted V11 gain difference | d_g,E = (g_base,arch + (g_M2 − g_base)) − g_TOPAS,arch. Archived ionmc gain 0.02210 / 0.02672 and TOPAS gain 0.01625 / 0.01842 (`de0h5/h5.json`) |
   | D3 | predicted V5 ΔIDD plateau-integral relative difference | p_E = (I_ionmc,arch + δ_E) / I_TOPAS,arch − 1. I_ionmc,arch = 1.9332 / 2.2848 and I_TOPAS,arch = 1.7520 / 2.1409 are the step-21 values in `validation/results/transport/lv5c-dd976c5-single-process.json`. δ_E is the change of the plateau-integrated IDD of the EM+elastic harness run, M2 minus base, at the CRN seeds, in the V5 ΔIDD units (same geometry and plateau window as V5) |
   | D4 | exceedance ratio M2 / Geant4 at 2 mm, 30–200 MeV | (M2 / table) from the D-E5 sampler (1e6 events per case; seeds **20991100 + case**, new, in the 2099xxxx family) × (table / Geant4) from the committed `de4/de4-ratios.json` |
   | D5–D7 | report-only | tail-moment and θ_lab median ratios; D-E0 kernels K* and G; the recoil share (H4) |

2. **CONTINUE iff all of the following hold** (each is necessary; together they are **not sufficient** for any gate):
   - **(i)** D1 lies inside the exact V11 F band at both energies.
   - **(i-b)** |q_E| ≤ 0.02 − t·SE_PP,E at both energies, the exact V5 peak/plateau band; equivalently, D1 inside 61.1–126.5 % at 150 MeV and 80.7–120.2 % at 200 MeV. Combined with (i), the effective F-closure window is **63.1–126.5 %** at 150 MeV and **80.7–120.2 %** at 200 MeV.
   - **(ii)** |d_g,E| ≤ 0.005 − t·SE_g,E at both energies.
   - **(iii)** |p_E| ≤ 0.10 − t·SE_I,E at both energies.
   - **(iv)** D4 moves toward 1, i.e. |ratio_M2/G4 − 1| < |ratio_table/G4 − 1|, at ≥ 4 of 5 energies. Current table / Geant4: 1.29–1.67 at 5 of 5 (review round D-E6, item 2).
   - **(v)** the §4.2 construction checks passed.

   Failing any one means **STOP**.
3. **Bands, from the frozen rule.**
   - The frozen rule (`validation/scripts/reference/compare_idd_v5.py`, `tost`, Amendment 9) passes iff the 90 % interval d ± t95·SE lies inside (−tol, tol). t95 is the one-sided 95 % Student-t quantile at the **floored** Welch–Satterthwaite df (`t95()`, table `_T95`).
   - Each band below is the set of point estimates whose interval would lie inside the tolerance **if the new run had the archived SE**. Archived SE and df are from step 16 (V11) and step 21 (V5) of the lv5c archive.

   | Clause | E | archived SE | df (floored) | t95 | allowed \|estimate\| | result |
   |---|---|---|---|---|---|---|
   | V11 F (tol 0.02, relative) | 150 | 0.002885 | 6.48 (6) | 1.9432 | 0.02 − 0.005607 = 0.01439 | closure band **63.1–136.9 %** |
   | V11 F | 200 | 0.003838 | 8.92 (8) | 1.8595 | 0.02 − 0.007136 = 0.01286 | closure band **73.6–126.4 %** |
   | V5 peak/plateau (tol 0.02, relative) | 150 | 0.004130 | 16.07 (16) | 1.7459 | 0.02 − 0.007211 = 0.01279 | closure band **61.1–126.5 %** (archived −3.66 %) |
   | V5 peak/plateau | 200 | 0.005482 | 7.41 (7) | 1.8946 | 0.02 − 0.010386 = 0.00961 | closure band **80.7–120.2 %** (archived −4.89 %) |
   | V11 gain (tol 0.005, absolute) | 150 | 0.000809 | 19.83 (19) | 1.7291 | 0.005 − 0.001399 = **0.00360** | – |
   | V11 gain | 200 | 0.000783 | 20.08 (20) | 1.7247 | 0.005 − 0.001351 = **0.00365** | – |
   | V5 ΔIDD plateau integral (tol 0.10, relative) | 150 | 0.03093 | 16.89 (16) | 1.7459 | 0.10 − 0.0540 = **0.0460** | archived +0.1035 (fail) |
   | V5 ΔIDD | 200 | 0.01881 | 21.24 (21) | 1.7207 | 0.10 − 0.0324 = **0.0676** | archived +0.0672 (pass, marginal) |

   For V11 F the closure band is c ∈ [(F_T(1 − a) − F_i)/(F_T − F_i), (F_T(1 + a) − F_i)/(F_T − F_i)], where a is the allowed |estimate| and F_i, F_T are the archived ionmc and TOPAS F (0.85678 / 0.89158 at 150 MeV, 0.84407 / 0.88737 at 200 MeV). For V5 peak/plateau, c ∈ [((1 − a)·PP_T/PP_i − 1)·F_base/Δ_obs, ((1 + a)·PP_T/PP_i − 1)·F_base/Δ_obs], with the D1b values.

   A run outside any band is already predicted to fail its clause, so it does not proceed. A run inside every band is not predicted to pass: the new SE, the M2 table's own V5 terms and the fresh seeds decide at lv5d.
4. **Disclosure: go/no-go, not blind** [new scope].
   - The screen uses the V11 reference (the archived TOPAS F and gain), the V5 reference (TOPAS ΔIDD) and the Geant4 thin-target ratios to decide whether the single pre-committed model is transported at all.
   - It is **not a selection among models**. M2 has no free parameter, plan §4.2 forbids switching to M3 on such outcomes, and STOP leads to "no model change".
   - It is, however, a stopping rule on gating quantities. **The lv5d V11 and V5 runs are therefore not blind to the diagnostic's estimates.** They replicate D1–D3 on independent seeds with the frozen comparators.
   - This disclosure is written into:
     - the diagnostic record header;
     - decision 0042;
     - the lv5d V5 / V11 records;
     - the results README.
5. **Honest prognosis from the recorded arms** (D-E6, `de6/de6-report.md`; report-only, not a rule). Observed gain excess: 0.00585 / 0.00830.

   | Shape arm | Δ_g at 150 MeV | residual at 150 MeV | Δ_g at 200 MeV | residual at 200 MeV |
   |---|---|---|---|---|
   | Geant4 CHIPS p+A shape (swap-pA) | +0.00694 | −0.00109 | +0.00393 | **+0.00437** |
   | both shapes (swap-both) | +0.00674 | −0.00089 | +0.00417 | **+0.00413** |
   | all, with the post-hoc σ_pp arm | +0.00813 | −0.00228 | +0.00694 | +0.00136 |

   - **At 200 MeV a shape as narrow as Geant4's would fail condition (ii)** (0.0041–0.0044 > 0.00365). With the archived TOPAS reference, the V11 gain clause at 200 MeV would also fail. Closing it in the recorded arms required the σ_pp change as well, i.e. H3, which (b)5(ii) does not authorise.
   - M2 passes (ii) at 200 MeV only if its shape lowers the 200 MeV gain by at least 0.00830 − 0.00365 = 0.00465, about 0.0005 more than the CHIPS shape did.
   - **STOP is therefore a probable outcome.** It costs minutes, and it saves the host time of (g). If STOP is reached on (ii) at 200 MeV, the record says so. A σ_pp change would need its own reviewed decision and amendment, with the lineage consequence of (b)5(ii) (the level would move to Geant4/TOPAS lineage). It is not decided here.
   - For (iii) at 150 MeV, a rough planner estimate (assumptions: a 20–60 mm plateau window and the archived plateau level 6.59) puts a CHIPS-like gain reduction at δ ≈ −0.17. The required shift is ≤ −0.1007. This is not a rule.
6. **STOP-and-record.**
   - The record (table id of the model under test, D1–D7, the failing condition) goes into decision 0042 and the results README.
   - The adoption of the model under test is withdrawn. No M3 follows, because §4.2 does not allow it on screen outcomes. "No model change" closure applies **on STOP at the diagnostic** for the model under test, **and** when no model remains under plan-of-record §4 (M3 the last fallback, M4 rejected), as pre-declared in (b″)4 and recorded in (b‴)3.
   - The outcome is **"no model change"**: V5 and V11 are not re-run; (d) lapses unused; lv 20531004 r_index 4–7 and 14 stay unused.
   - Under the gate set of (e), integration fails on V5/V11. That is a recorded outcome.
   - A different model, or a σ_pp change, needs a new reviewed physics reason and a new amendment.
7. **Cost.** Sandbox, single process, with no host job and no other heavy job running: two harness runs at 15–43 s each (`wall_s` in `de3/de3-summary.json`), plus 10 sampler cases. Minutes.

## (d) One fresh V5 base — revision of Amendment 17 (b)4 [new scope; **lapses unused** in the state of (b‴)]

1. **Background.** Amendment 17 (b)4 applied after the stop rule fired at dd976c5a. The intervention it required was raised (IR-20261017-223039-9DA564) and answered (a).
2. **Exactly one** fresh V5 base: **lv 20531004**, V5 at lv5d r_index 4–7.
   - Conditions: only after CONTINUE, and only with M2.
   - **Unchanged:** tolerances, reference runs, comparator, attribution ladder and stop rule (Amendment 14 (f); Amendments 9, 10).
3. **Preserved:** the V5 failures at b40d812, 7aae5bb and dd976c5a, and the dd976c5a stop-rule record.
4. **No second base** without a new amendment. If the stop rule fires again, 17 (b)4 applies as written and the remaining steps are archived.

## (e) Merge gate (a)6 for V3-005B/C/D under Option II [new scope; in the state of (b‴), V5 and V11 cannot pass, so see (b‴)4 for the closure amendment, which requires operator approval]

**Scientific gates, at one evidence SHA:**

| Gate | Status and conditions |
|---|---|
| V5 | lv 20531004, M2 |
| **V11** | **added**; frozen tolerances; lv5d r_index 14 |
| **V7-R2** | **in V7-R's slot**; repair plan revision 4, measured joint-power rule. V7-R stays recorded as failed |
| V8-nuc | hr5d |
| **V6-F** | **binding**, executed once at base 20551004; v6f-followup revision 2 unchanged, including the coverage minimum |

**Implementation gates:** P6-D (on the new M3 table; p-p level unchanged), P7 (on the new M3 table, lv5d r_index 15), the M3 verification of (b″)4, P5-ext, C1-ext, R1-D2, V8-LV, CI; the M3 table's Amendment 16 qualification flags.

**Executed, not gates:** V2b-el, V7-scan-el, V7-shift-el, E1-B, V10 (report-only; p-p level unchanged), X-ENDF, D9-BGG, E-shape.

**Carried as recorded failures** (shown next to every pass):
- V6 / V6-sens: failed — not sized, re-recorded by `v6-unsized`;
- V10-A at dd976c5a; the new-model V10-A is report-only (data seen);
- V7-R;
- the dd976c5a V5 / V11 / V7-scan-el; lv5b V5 and V7; P6 (H-1 below 12.53 MeV); the MF6 rule.

**Claims.**
- No lateral-halo claim beyond V6-F's reduced scope. No measured-dose claim from V6.
- The p+A shape statement "agrees with Geant4 `g4h-elastic_HP` within the V11 tolerances" is made **only if** the lv5d V11 passes. It is a model-to-model statement, not a data statement.
- σ_el is Barashenkov, shared with TOPAS; p+C at 100 MeV is +17 % above the post-1997 proxy (recorded deficiency).

**Archive semantics.** lv5d `pass: false` (V6 rows); `conformant: false` on both (single-process directive; imported partials; deferred partition step). Neither flag is the merge criterion.

**Merge condition: "revised applicable gates pass"**, operationally, at the PR head:
1. every gate above passes at the evidence SHA;
2. `record_local_validation` passes at that SHA, and a truthful head-check record passes at the PR head;
3. the independent review of the closure and Codex pass;
4. CI is green;
5. the recorded failures are listed;
6. no implementation-only passing record exists.

## (f) Suites, bases and seeds [new scope]

1. **Suites.**
   - `lv5d` / `hr5d` are new. The lv5c / hr5c definitions are kept.
   - **lv5d** is the lv5c row set with: V7-R → V7-R2 (`v7r-kappa`, 78 calibration parts, 48 transport shards, reference, combine, diagnostic); the V6-F CPU rows (gating and sensitivity) added; V6 / V6-sens re-recorded by `v6-unsized`.
   - **hr5d**: V8-stat (python × 4, cpu64, cuda32, cuda64, compare), V7 f32/f64 CUDA, `v6-unsized` / `v6-cuda-compare`, and the V6-F CUDA replicate and timing step.
2. **Bases.**

   | Base | Use |
   |---|---|
   | lv **20531004** | lv5d rows |
   | hr **20541004** | hr5d rows |
   | **20551004** | V6-F rows (a row-specific base, like the absolute E-shape seeds of Amendment 17 (b)3) |

   20511004 / 20521004 are **retired unused**: declared before observation, no seed drawn, no repetition. The sequence continues at 20561004 and skips 2098xxxx / 2099xxxx.
3. **lv5d r_index** (as lv5c):

   | r_index | Row | Seeds drawn |
   |---|---|---|
   | 1 | throughput | yes |
   | 2 | V8-LV | yes |
   | 3 | R1-D2 | yes |
   | 4–7 | V5 | yes |
   | 8 | V2b-el | yes |
   | 9 | V7-scan-el | yes |
   | 10 | V7-shift-el | yes |
   | 11 | unused | – |
   | 12 | V6 | no seed |
   | 13 | E1-B | yes |
   | 14 | V11 | yes |
   | 15 | P7 | yes |
   | 16 | V7-R2 transport | yes |
   | 17 | V7-R2 diagnostic | no seed |
   | 18 | V10-A cross-check | yes |
   | 19 | X-ENDF sampled cross-check | yes |
   | 20 | V6-sens | no seed |

4. **hr5d uses r_index 1, 2 and 3 only:** 1 V8-stat; 2 V7-f32; 3 V6-CUDA, unsized, so no seed is drawn. Bases 10 000 apart make lv r ≡ hr r − 10 at equal k. Therefore **every other hr r_index at 20541004 is reserved and unused**, in particular:

   | hr r_index | Coincides with | Treatment |
   |---|---|---|
   | 0 | lv 10 | reserved |
   | 4 | lv 14 (V11) | reserved |
   | 5 | lv 15 (P7) | reserved |
   | 6 | lv 16 (V7-R2 transport, 20547004–20547052; 7200 fallback to 20547068) | reserved |
   | 8 | lv 18 | reserved |
   | 9 | lv 19 | reserved |
   | 17–19 | V6-F r_index 7–9 at 20551004 | reserved |
   | 1–3 | lv 11–13 | the existing pattern of Amendment 14 (j): lv 11 is unused; lv 12 and hr 3 draw no seed; lv 13 (E1-B) coincides with hr 3, which draws none |

5. **Seeds of record:**

   | Item | Seeds |
   |---|---|
   | V7-R2 transport | 20547004–20547051, reference 20547052 (7200 fallback: k = 0..63, reference k = 64 → 20547068); repair plan revision 4 §5 |
   | V7-R2 κ / calibration roots | 20980000 / 20982000 / 20984000 |
   | V6-F gating CPU (r_index 7, k = 0..15) | **20558004–20558019** |
   | V6-F sensitivity (r_index 8, k = 0..7) | **20559004–20559011** |
   | V6-F CUDA (r_index 9, k = 0..15) | **20560004–20560019** |
   | Diagnostic D4 | 20991100–20991109 (V-M3); 20991200–20991209 (V-M2, report-only) |

   **Required one-line correction in the same commit:** `validation/plans/v3-005d-integration-route.md` §4.4 (table row V6-F) and Appendix A print "20558004 … 20558015". This becomes "20558004–20558019" in route **revision 7**.
6. **Rehearsal, timing and size re-check (before the evidence SHA):**
   - **V7-R fixture** regenerated at the evidence digest as 4 × 2.25e6 histories at rehearsal seeds **20521000–20521003** (20505000 + 16000 + k; digest rule of Amendments 17 (a)7, 18 (d)–(e)). Never evidence.
   - **V6-F CUDA timing step**, timing only, measured first: one 1e6-history CUDA float32 shard of the V6-F beam at rehearsal seed **20514000** (rehearsal base 20505000, r_index 9, k = 0; unused by the recorded rehearsals, route §4.4 scan). It fixes the CUDA shard sizes (≤ 1200 s planned) and the budget line.
   - **Joint-power re-check** (repair plan revision 4, §3.8) on the rebuilt pool, sandbox: 5400/0.005 iff ≥ 378/400 for both surrogates; else 7200/0.005 iff ≥ 378/400; else operator.
   - Then the sizes addendum (14a-style; nothing decreases) is committed.
7. **Seed-family test** (`tests/ionmc/test_v3_005d_seed_families.py`), extended to assert:
   - all seeds of item 5, the timing seed 20514000 and the fixture seeds;
   - disjointness from `consumed_seed_set()`, the recorded rehearsal seeds and 2098xxxx / 2099xxxx;
   - every reserved hr r_index of item 4;
   - checks (a)–(e) of v6f-followup §7.
8. **Implementation items:**
   - lv5d / hr5d qualification-base constants, `summarize.seed_blockers`, `DEFAULT_SEED_BASES`, the V6-F row base;
   - `run_suite.full_step_names`;
   - `rehearsal_preflight` for lv5d;
   - the single-process workers = 1 force list;
   - row mapping and pending rows;
   - a V7-R2 reference index distinct from every shard index (the constant `V7R_REF_INDEX = 16` collides with shard 16);
   - `declared_plans()` entries, or a timing rule, for the rep-based κ / calibration parts, with one short calibration part measured on the slow host;
   - calibration `PARTIAL` lines;
   - the M3 table pin and its Amendment 16 flags and (b″)4 verification, with P6-D in CI and P7 in lv5d;
   - V6-F code in the hashed source set;
   - all of the above before the code freeze.

## (g) Order of work, failure branch and budget (one heavy job at a time; stability restrictions retained) [new scope]

| # | Step | Where | Cost (source) |
|---|---|---|---|
| 1 | Acquisitions, **all complete**: SAID SM16 (10–300 and 5–300 MeV tables, pp DSG scans 100–300 MeV); Horiuchi preprint; density sources (Shukla, Mehndiratta & Shukla, TGlauberMC paper and tarball, Duda) and the IAEA radii (ids in (b′)) | MCP | done |
| 2 | Decision 0042 | repo | – |
| 3 | M2 builder, §4.2 checks, Amendment 16 flags, tests, Codex | sandbox | – |
| 4 | M2 diagnostic (c); **STOP/CONTINUE record** committed, with the disclosure header | sandbox | minutes |
| 5 | On CONTINUE: implementation (f)8, Codex, CI, sandbox suite in two logged halves; then code freeze | sandbox | – |
| 6 | Fixture rebuild (host), fixture commit | host | ≈ 0.7 h (4 × 2.25e6 at the slow rate; route §4.2) |
| 7 | V6-F CUDA timing step | host | one shard; planning ≈ 37 s (hr5c step 06: 1e6 CUDA f32 histories in 22 s at 150 MeV × 1.672); measured, then recorded |
| 8 | Joint-power re-check | sandbox | ≈ 2.2 h; ≈ 3 h more if 7200 |
| 9 | Sizes addendum, independent review and Codex at the **evidence SHA** | – | – |
| 10 | **V6-F first:** 72 CPU shards, then 16 CUDA shards; then `v6f-compare` | host | CPU ≈ 6.15 h (72 × 184 s × 1.672; v6f-followup §8.3). CUDA from step 7 (planning ≈ 16 × 37 s ≈ 0.16 h) |
| 11 | Only if V6-F passes: lv5d (throughput first), then hr5d | host | ≈ 37.2 h (≈ 47.3 h at 7200; route §4.2) |
| 12 | Closure | – | – |

**Failure branch: V6-F fails, or its halo verdict is "not established (beam-sensitive)"** [new scope]. The coverage minimum (≥ 40 gating points, ≥ 3 radii, ≥ 12 sign-test depths after the screen; planning expectation ≈ 58 points and 16 depths, Amendment 18 (b)3) is the most likely failure mode.
- **Merge is impossible under Option II.** V6-F is executed once; there is no fall-back (route §6 S7).
- **lv5d / hr5d are deferred, not run.**
  - Reason: no lv5d / hr5d outcome can make the merge possible.
  - The operator directs "inexpensive diagnostics before expensive qualification", and the V6 lesson applies: do not spend host time that cannot change a gate.
  - Their seeds stay unobserved.
  - They may be executed later at the same evidence SHA without a new amendment, as archived science, if the orchestrator records a reason. They never count toward a merge on this route.
- **Closure wording:** "V3-005 integration blocked: binding gate V6-F failed (or: not established) at <SHA>. V5 / V11 / V7-R2 not executed (deferred; seeds unobserved). Recorded failures preserved. PR #82 stays open. No passing record."
- Any change of the gate set afterwards is a new, post-hoc operator decision.

**Budget:**

| Branch | Host hours | Sandbox |
|---|---|---|
| Full, 5400/0.005 | ≈ 0.7 + 0.01 (timing) + 6.15 + CUDA (≈ 0.16 planned) + 37.2 ≈ **44.2 h** | minutes (diagnostic) + ≈ 2.2 h (joint power) |
| Full, 7200/0.005 | ≈ **54.3 h** | – |
| V6-F fails | ≈ **7.0 h** (0.7 + 6.15 + CUDA) | – |
| STOP at the diagnostic | **0 h** | – |

D-E4 is done and is not re-budgeted.

## (h) Unchanged, restated

- The frozen tolerances of V5, V11, V6-F, V8, R1-D2 and P7.
- The V7-R2 clauses (joint ≤ 0.02, gate alone ≤ 0.04, power ≥ 0.9, R = 3200, B = 1999).
- The single-process directive and the operator diagnostic rules.
- The crash rule (re-execute once at the same seeds a step that published no document; never re-run a step with a published verdict).
- The Amendment 19 records stay exploratory.

**Amendment 21 (2026-10-18; draft 8, independently reviewed in rounds A21-1 and A21-2 of `validation/plans/v3-005d-elastic-attribution/review.md`; registered before any builder table, transport or diagnostic of the candidate exists; record files in `validation/plans/v3-005d-m2prime/`) — second-generation elastic candidate M2′: optical-limit Glauber with the Clementel–Villi Pauli-blocking factor (reading R-C), applied at table nodes with E ≥ 80 MeV; trigger 3 of Amendment 20 unchanged as the sole input-sanity gate.** The change tables below are the review history. Earlier drafts are kept in the session scratchpad and are not part of the record.

**Changes in draft 8 (review round A21-2, applied verbatim).** (f): σ_R/Barashenkov is withheld for M2 and for every variant with a printed % deviation in Git; the residual-inference disclosure is added. (e): the M3 → M2′ mapping is restricted to the listed phrases, with the D4 seed statement.

**Changes in draft 7.** Git policy (f): no absolute σ_tot, σ_R or σ_el of any variant is projected. This follows the elastic-attribution projection rule; note σ_R + σ_el = σ_tot. The builder-vs-scratch record is projected as relative differences only.

**Changes in draft 6 (review round A21-1).**
- (1) The node count is 43 for both targets.
- (2) Trigger 3 is framed as an implementation check. A builder-vs-scratch reproduction record is added, and the re-tuning clause is tightened.
- (3) New limitations: the nonrelativistic convention; the black disk kept below 80 MeV; the 80–100 MeV shape; the 78.81–80.00 MeV blend region.
- (4) Git policy: no measured EXFOR values anywhere. σ_tot is compared only as pass/fail per cell; the scripts read the cache at run time.
- (5) An explicit Amendment 20 M3 → M2′ reference mapping on CONTINUE.

**Changes in draft 5.**
- The Clementel–Villi factor is pinned in full, with both branches and the flux-weighted definition. It is verified by an independent Monte Carlo derivation (`validation/plans/v3-005d-m2prime/M2PRIME-projected.md`, section "R-C verification", script `cv_mc.py`), because no open printed source was obtained.
- The conventions are pinned.
- **Energy domain:** M2′ applies at table nodes with E ≥ 80 MeV only. Below that the black disk is kept. The extension down to 5 MeV is a report-only sensitivity.
- The disclosure records when the domain was chosen.

**Changes in draft 4.**
- **The model changes from Bertulani & De Conti Eq. (14) to the Clementel–Villi single-Fermi-sphere Pauli factor** (reading R-C). The reasons are the paper's own text and the physics, recorded in `validation/plans/v3-005d-m2prime/M2PRIME-projected.md`, section "Reading of ρ₁, ρ₂", with the script `m2prime2.py`.
- Eq. (14) under reading R-B becomes a report-only sensitivity.
- The disclosure quotes the R-C scratch numbers and gives, in order, the sequence of readings considered.
- Trigger 3 of Amendment 20 is unchanged: the sole gate, both targets, ≥ 2 of 3 energies.
- Everything else is as in draft 3.

---

**What this amendment does not change:**
- Amendments 1–20 or any recorded verdict;
- the plan of record `validation/plans/v3-005d-elastic-attribution.md` (sha256 fd9657df…);
- the "no model change" outcome of Amendment 20 (b‴), which stands unless M2′ passes every step below;
- V5, V11 and V10-A at dd976c5a, which stay failed.

**[new scope]** throughout.

## (a) Disclosure (mandatory)

1. **Timing.** M2′ was defined:
   - after M2 failed plan §4.2 trigger 3 (O-16 1.407 at 150 MeV, 0.5 % above the bound; Amendment 20 (b″));
   - after D-E2–D-E6 were seen (Amendment 19 (a)–(d));
   - after the cause of the M2 failure was verified: free-NN optical-limit Glauber overpredicts the measured n+C / n+O σ_tot by +16–19 % at 150–200 MeV (EXFOR 13569, 13753).
2. **The sequence of readings considered, in order.** The form, the reading and the outcomes below were **all seen before registration**.
   1. **R-B:** Bertulani & De Conti, Phys. Rev. C 81 (2010) 064603, Eq. (14), with ρ₁ and ρ₂ the local proton and neutron densities and ρ̃ = ρ_tot/ρ₀. This was draft 2's model. Scratch σ_el/Barashenkov: C 0.717 / 1.099 / 1.235; O 0.699 / 1.182 / 1.396. O-16 sat on a knife edge (margins −0.001 and +0.004).
   2. **R-A objection** (independent reviewer): read ρ₁ = 0 for the free projectile. **Rejected.** It contradicts the paper's text, which defines ρ₁ and ρ₂ as the proton and neutron densities. It is also unphysical: it has no free limit, and it removes about 60 % of the free-NN Glauber σ_tot (M2PRIME, R-A against M2).
   3. **Paper check** (acquired text, arXiv:1004.2096):
      - Eq. (14) is a fit to the Pauli-blocking geometry for **two bound nucleons**: two Fermi spheres, with the bound-bound geometry of the paper's Eq. (7).
      - The paper applies it with an "effective local density … where the two densities (one from the projectile and the other from the target) cross each other". That is a nucleus–nucleus application.
      - It does **not** apply Eq. (14) to a free nucleon on nuclear matter.
      - For nucleon–nucleus scattering it names Clementel & Villi as the standard: "an analytical expression for the scattering of a nucleon on a nucleon Fermi gas … widely used to describe Pauli-blocking in nucleon-nucleus scattering".
   4. **R-C: the Clementel–Villi single-Fermi-sphere factor, adopted as M2′.** The knife edge of R-B is recorded as the **reason the reading was re-examined**. It was not a selection on the band: the switch is grounded in the paper's text (step 3), and R-C was adopted with its outcomes in view, as disclosed here.
3. **R-C outcomes known before registration** (scratch implementation `validation/plans/v3-005d-m2prime/scripts/m2prime2.py`):

   | Quantity | C-12 (100 / 150 / 200 MeV) | O-16 (100 / 150 / 200 MeV) | Reading |
   |---|---|---|---|
   | σ_el^M2′/σ_el^Barashenkov | **0.595** / 0.976 / 1.129 | **0.595** / 1.062 / 1.286 | trigger 3: 100 MeV fails for both, 150 and 200 MeV pass with robust margins. Each target is inside the band at 2 of 3, so trigger 3 would pass |
   | σ_tot against measured n+natC / n+natO, ±10 % band (report-only; pass/fail per cell, (f)) | fail / pass / pass | fail / pass / pass | ±10 % at ≥ 2 of 3: passes; at all 3: fails |
   | σ_R/σ_inel^Barashenkov (report-only) | 1.090 / 1.010 / 0.953 | 1.061 / 0.984 / 0.936 | – |
   | Trigger 4: q_min/q_H | 0.805 at 150 MeV; 0.819 at 200 MeV | 0.815 at 150 MeV | passes ([0.6, 1.0]) |

   **This amendment is therefore NOT a blind pre-registration at the trigger level.**
4. **Expectation:** trigger 3 passes for both targets. Trigger 2 (numerics) is not yet evaluated on the builder.
5. **What remains blind:**
   - the transport go/no-go diagnostic of Amendment 20 (c): D1–D3 transport quantities at the CRN seeds, and D4 at the new seeds 20991300 + case;
   - every gate.
6. **The diagnostic is not naive either.**
   - The medium correction barely changes the angular shape. O-16 P(θ_lab > θ_2mm) is 0.528 / 0.604 at 150 / 200 MeV under R-C, against M2's 0.540 / 0.612, the black disk's 0.688 / 0.735, and Geant4 CHIPS's 0.513 / 0.571 (CHIPS for information only).
   - The diagnostic therefore effectively tests "Glauber shape + Pauli-blocked rate" against the black disk.
   - D-E3's `pA-s0.9` probe (84 % / 76 % gap closure) is prior evidence in the same direction.
6b. **The energy domain and the form of the factor.**
   - The energy domain of (c)7 (M2′ only at nodes with E ≥ 80 MeV) was chosen **after** the branch analysis of the Clementel–Villi factor was seen (`M2PRIME-projected.md`, "R-C verification") and **before** any M2′ build or trigger evaluation.
   - It is a physics choice: only the high-energy branch is exercised, inside the eikonal validity domain. **It is not tuned to the band:** the trigger nodes at 100, 150 and 200 MeV lie above 80 MeV and are unaffected by it.
   - **No open printed source** for the full Clementel–Villi factor, including its low-energy branch, was obtained. Hussein–Rego–Bertulani 1991 and Kikuchi–Kawai 1968 were not acquirable; Kelly (nucl-th/9607001) and Chen–Sammarruca–Bertulani (arXiv:1304.6096) were checked and do not contain it.
   - The form of (c)2 is therefore pinned from the researcher's **independent Monte Carlo derivation record**, not from a printed equation, and this amendment says so. The derivation used 4×10⁶ samples per point and agrees with the form within 4×10⁻⁴ on both branches; both branches give 0.3000 at E = 2E_F, and the factor is zero at E = E_F.
7. **Alternatives that remain open, named.**
   - **(i)** M3, the plan-of-record fallback: not excluded on physics grounds. It was "not acquirable within the experiment's tool and the routes tried" (Amendment 20 (b‴)1). An **operator-assisted acquisition** (raising the 256 MiB tool cap, or manual placement of ENDF/B-VIII.0 MAT 625 / 825) would make it constructible. Its coverage to 135–200 MeV is unverified.
   - **(ii)** The closure route without a model change (Amendment 20 (b‴)4).
   - M2′ is the autonomous candidate. The M3 route is an operator or infrastructure matter.
8. **Selection on physics grounds.** The standard Pauli-blocking correction for nucleon–nucleus scattering, as the cited paper names it. **Never** selected on V5, V11, D-E5 or D-E6 outcomes.
   - Rejected on physics grounds (M2PRIME):
     - (b2), ρ̃ = 2ρ/ρ₀;
     - R-A;
     - α scaling (no published prescription);
     - χ scaling (a fit without physics basis);
     - α_NN = 0;
     - the Horiuchi 2007 Table I set: it is the isospin-averaged free amplitude (at 150 MeV σ 3.845 against 3.893 fm², α 1.245 against 1.168), with M2's own β rule.
   - **The E + |U| refinement is explicitly not adopted:** it would be tuning after the outcome is seen.

## (b) Known limitations, recorded now

1. **100 MeV over-correction.** With the bare lab energy E, the factor ignores the potential-energy gain of the projectile inside the nucleus. R-C over-corrects at 100 MeV: σ_el/Barashenkov is 0.595 for both targets, and σ_tot falls outside the ±10 % band below the measured value (pass/fail record, (f)). The E + |U| refinement is not adopted ((a)8).
2. **Single Fermi sphere.** Clementel–Villi treats a free nucleon on a Fermi gas, with isotropic NN scattering and E_F from the local total density (N = Z assumed in E_F; the isospin asymmetry is not resolved).
3. **Low-energy behaviour and domain.**
   - At the interior peak point densities (C 0.179, O 0.165 fm⁻³), 2E_F is 79.5 / 75.2 MeV and E_F is 39.7 / 37.6 MeV.
   - Below about 75–80 MeV the low branch is active in part of the nucleus. Below about 38–40 MeV the interior is fully blocked (f = 0), and the phase comes from the surface only.
   - On the pinned grid (1–250 MeV, 400 nodes), C has 122 and O 120 nodes in [5 MeV, 2E_F), and C 95 and O 91 nodes in [5, E_F].
   - These regimes are **outside the model's domain**: see (c)7.
4. **Optical limit.** The optical-limit approximation weakens below about 100 MeV (plan §4.3). The real part of the free NN amplitude (α) is kept. The factor acts through σ in χ = ½σ(i + α)…, so it scales the real and the imaginary phase together.
5. **Nonrelativistic derivation (declared convention ambiguity, not adopted).** The factor is derived nonrelativistically, with E the bare kinetic energy. Using k₀²/2m from the relativistic momentum instead changes x by 7–9 % at 150–200 MeV, i.e. Δf ≈ +0.026 (f 0.642 → 0.668 at 150 MeV, 0.731 → 0.757 at 200 MeV; review round A21-1). It does not affect trigger 3 at the known margins.
6. **The black disk below 80 MeV keeps part of the over-correction.** D-E0 shows that p+O removal below 80 MeV is a non-negligible part of the kernel K. With the black disk kept there ((c)7), **the diagnostic closure is expected to fall short of what a full-range M2′ would give** (V-M2′ext is construction-level only). A later STOP is therefore not to be read as a failure of the in-medium physics alone.
7. **The 80–100 MeV shape comes from the regime where R-C is worst.** At 100 MeV σ_el/Barashenkov is 0.595: the model over-blocks there. Only the normalised shape enters transport (the level is Barashenkov), but the shape at 80–100 MeV comes from that regime.
8. **Blend region at the domain boundary.** The pinned nodes next to 80 MeV are 78.81 and 80.00 MeV. The transport interpolates the bin edges linearly in E between nodes, so **78.81–80.00 MeV is a blend of the black-disk and M2′ shapes**, not a point jump. The 80.00 MeV node is included and lies above 2E_F(C) = 79.5 MeV, so it uses the high branch.

## (c) M2′ construction (pinned)

1. **Amplitude and phase.** As M2 (Amendment 20 (b′)): pp and pn phases summed (Abu-Ibrahim et al. 2008, Eqs. (4)–(5)), with σ_c evaluated locally:

   χ(b) = ½ Σ_{c ∈ {pp, pn}} (i + α_c) ∫ d²b′ f_c(b − b′) ∫ dz σ_c^free(E) · P_CV(E; ρ_tot(r)) · ρ_c(r),  with r = √(b′² + z²)

2. **Clementel–Villi factor**, the flux-weighted definition. Lineage: Clementel & Villi 1955; Hussein, Rego & Bertulani, Phys. Rep. 201 (1991) 279. **Pinned from the independent derivation record** (`M2PRIME-projected.md`, "R-C verification"; `cv_mc.py`); no printed source was obtained ((a)6b).

   With x = E_F/E:
   - f = 1 − 7x/5 for E ≥ 2E_F;
   - f = 1 − 7x/5 + (2x/5)(2 − E/E_F)^{5/2} for E_F < E < 2E_F;
   - f = 0 for E ≤ E_F.

   **Definition.** f is the flux-weighted, per-path-length suppression ⟨(|k₀ − k₂|/k₀)·Q⟩, the quantity that enters the Glauber attenuation ρσ per unit length. It is **not** the plain blocked fraction ⟨Q⟩, which differs: 0.2711 against 0.3000 at E = 2E_F.

   **Conventions:**
   - E_F = (ħc)² k_F²/(2 m_p c²), k_F = (3π² ρ_tot/2)^{1/3}, with ρ_tot = ρ_p + ρ_n the local point density in fm⁻³;
   - ħc = 197.3269804 MeV fm and m_p c² = 938.27208816 MeV, giving E_F(0.17 fm⁻³) = 38.392 MeV;
   - **E is the bare proton lab kinetic energy**;
   - the same factor applies to pp and pn.
3. **Inputs.**
   - σ_c^free and α_c: SAID SM16 (Amendment 20 (b′)1, (b′)5), free α.
   - β_c: the free unitarity rule (Horiuchi Eq. (10)) with the free σ_c and α_c.
   - Densities: Amendment 20 (b′)3 with the unfolding rules of (b″)3.
4. **Level and targets.** The shape |F|² is normalised and scaled to the **Barashenkov σ_el(p+A)**, unchanged. Coulomb stays in the MCS, and the p-p level is unchanged. All seven targets are built; only O and C enter V5, V11 and V10-A. Nodes outside the domain of (c)7 keep the pinned black-disk shape.
5. **Report-only sensitivity variant V-RB:** Eq. (14) of Bertulani & De Conti under reading R-B (ρ₁, ρ₂ the local proton and neutron densities; ρ̃ = ρ_tot/ρ₀; ρ₀ = 0.17 fm⁻³).
   - It is built and reported at construction level (trigger quantities, σ_tot, σ_R, q_min, P) and in the D-E3 harness at the CRN seeds.
   - It is never a gate, and D4 is not run for it.
6. **Sources and acquisition.**
   - arXiv:1004.2096 (the paper check) is acquired through `acquire_reference_data`, role `construction`.
   - The Clementel–Villi factor is pinned from the committed derivation record ((c)2; (f)), **not** from a printed source, and decision 0042 states this. No construction-impossible condition attaches to the printed source.
   - If an open printed source is obtained later, it is compared with the pinned form and the comparison is recorded. A disagreement goes to independent review before any use.
7. **Energy domain** (orchestrator decision; physics-based; declared before any M2′ build or trigger evaluation, (a)6b).
   - **M2′ replaces the black-disk shape and rate only at table nodes with E ≥ 80 MeV.** That is about 2E_F at the peak density of C and O, so only the high branch is exercised, inside the eikonal validity domain. On the common pinned grid this is **43 nodes for every target**.
   - **Nodes below 80 MeV keep the pinned black-disk table for all targets.** The low-energy transport is therefore unchanged from decision 0041's D-1 table.
   - **Disclosed:** the shape changes at the domain boundary, with the 78.81–80.00 MeV node interval acting as a blend region ((b)8). The σ_el level is Barashenkov on both sides.
   - **Report-only sensitivity V-M2′ext:** M2′ extended down to 5 MeV, built and reported at construction level only. It is never the model.
   - **Other targets** (Al, Si, P, Ca): their 2E_F thresholds at peak density are expected near 77 MeV but are **unverified** (recorded). The common 80 MeV boundary applies to all targets.

## (d) Triggers

- **Trigger 3, unchanged from Amendment 20 / plan §4.2, the sole input-sanity gate.** σ_el^M2′/σ_el^Barashenkov within **[0.7, 1.4]** at ≥ 2 of 100 / 150 / 200 MeV, **for C-12 and for O-16**.
  - **Framing (round A21-1).** With the outcome known ((a)3) and robust margins (in band at 150 and 200 MeV by ≥ 8 %, out at 100 MeV by 15 % for both targets), **trigger 3 carries no evidentiary weight for M2′.** Its remaining function is an **implementation check of the builder**.
  - It is evaluated on the **builder's official construction only**. The scratch numbers of (a)3 are never used for it.
  - **No re-tuning,** in particular no switch to E + |U|, no relativistic k₀²/2m and no change of reading after the builder's numbers are seen. The only exception is a defect fix that meets two conditions: it comes from an independent review with a reproducing test that does not depend on the trigger value, **and** it does not move the builder toward the scratch values except through the identified defect.
  - **Builder-vs-scratch reproduction record** (pre-declared). At the six trigger cells (C and O at 100 / 150 / 200 MeV), σ_tot, σ_R, σ_el, σ_el/Barashenkov and q_min of the builder are compared with the scratch values of `m2prime2.py`, and the differences are reported. **A difference above 2 % in σ_tot, σ_R or σ_el goes to independent defect review** before the trigger outcome is recorded. The trigger outcome is the builder's, whatever the review finds, unless a defect satisfying the clause above is fixed.
  - **The outcome is accepted whichever side it lands on.**
- **Triggers 2 and 4, unchanged:**
  - trigger 2: normalisation ≤ 1e-6; b-, z- and q-integration convergence ≤ 1e-3;
  - trigger 4: q_min within [0.6, 1.0] × q_H (q_H from Amendment 20 (b″)3). C-12 at 150 MeV is "not evaluable" if it has no true minimum, and 200 MeV is used instead.
- **No new σ_tot trigger** is introduced.
- **Report-only sanity records** (never gates):
  - σ_tot^M2′ against the EXFOR 13569 / 13753 average within a ±3 % energy window, recorded as **pass/fail per cell** against the ±10 % band, with both readings (≥ 2 of 3 and all 3). No % deviation is printed ((f));
  - σ_R^M2′ against the Barashenkov σ_inel(p+A).
- **If any trigger fails, M2′ is withdrawn.** There is no fallback. "No model change" stands, and the closure route of Amendment 20 (b‴)4 follows (operator approval required).

## (e) Diagnostic, gates, seeds, budget, order (as in draft 3)

- **Diagnostic.** Amendment 20 (c), with every band unchanged (V11 F, V5 peak/plateau, V11 gain, V5 ΔIDD, D4 direction; the go/no-go and non-blind disclosures). Variants:
  - **V-M2′ (R-C), the model under test:** D-E3 CRN seeds 20990100 / 20990101; **D4 seeds 20991300 + case**.
  - **V-M2, report-only:** D4 seeds 20991200 + case.
  - **V-RB, report-only:** construction level and harness only.
  - **V-M2′ext, report-only:** the 5 MeV extension of (c)7, construction level only.
- **CONTINUE:** Amendment 20 (d) (one fresh V5 base, lv 20531004), (e) (the Option II gate set, with V5 and V11), (f) and (g) apply with M2′ as the model. Decision 0042 receives a dated section adopting M2′.
- **Reference mapping on CONTINUE (explicit; restricted, round A21-2).** In Amendment 20 revision 7 ((c), (e), (f)8), **only** the following phrases read as **M2′**:
  - "the model under test (M3)";
  - "V-M3";
  - "the new M3 table" / "the M3 table";
  - "the M3 verification of (b″)4" and "(b″)4 verification".

  Their meaning then becomes: the M2′ table; the trigger 2–4 records plus the builder-vs-scratch reproduction record of (d), in place of the M3 verification; and P6-D, P7 and the Amendment 16 flags evaluated on the M2′ table.
  - **Every other "M3" stays literal.** In particular "plan §4.2 forbids switching to M3 on such outcomes" and "No M3 follows" keep M3 as the plan's fallback; M2′ has no fallback.
  - **D4 seeds.** For M2′, D4 uses **20991300–20991309**. The block 20991100–20991109, given to V-M3 in the installed plan's (c) and (f) seed lines, **stays reserved and unused** (V-M3 was never run).
- **Seed-family test:** D4 seeds 20991300–20991309 are registered.
- **STOP:** "no model change" stands, and the closure route of Amendment 20 (b‴)4 follows.
- **Seeds and bases:** unchanged (lv 20531004, hr 20541004, V6-F 20551004; V7-R2 20547004–20547052).
- **Budget:** minutes of sandbox time (acquisition, build, triggers, diagnostic); host time on CONTINUE only, ≈ 44.2 h plus V6-F CUDA (Amendment 20 (g)).
- **Order:**
  1. this amendment, reviewed independently;
  2. acquisitions;
  3. builder `--model m2prime`, with the factor and the conventions pinned from the derivation record ((c)2) and the 80 MeV domain ((c)7);
  4. tests and Codex review;
  5. triggers;
  6. diagnostic;
  7. the STOP/CONTINUE record.

  No host step happens before CONTINUE.

## (f) Git policy for the M2PRIME record (revised in draft 6)

The researcher's assessment (all readings; the sections "Reading of ρ₁, ρ₂" and "R-C verification") is committed as a projected table in a plan-record directory (path decided at commit, for example `validation/plans/v3-005d-m2prime/`).
- **Projected (revised in drafts 7–8, following the projection rule already in force for the elastic-attribution records: the `common.py` allowlists):**
  - σ_el/Barashenkov and σ_R/Barashenkov ratios, **with these exceptions (round A21-2):**
    - **σ_R/Barashenkov is not projected for M2 (variant (a))**, nor for any variant that has a printed % deviation from measurement anywhere in Git. The Barashenkov arrays are committed in `src/ionmc/nuclear/bgg.py`, so σ_R and σ_el ratios together with Amendment 20's M2 percentages would recover the measured σ_tot to about 0.5 %.
    - For the constructions fitted to the measured σ_tot ((c) α scaling, (d) χ scaling, (e) Pauli + α scaling), both ratios are withheld.
    - R-C (M2′), R-B, R-A and (b2) carry pass/fail only and keep both ratios;
  - the σ_tot comparison as **pass/fail per cell** against the ±10 % band;
  - q_min and q_min/q_H;
  - P(θ_lab > θ_2mm);
  - the Pauli factors;
  - the Monte Carlo derivation table of the Clementel–Villi factor;
  - the defining sentences quoted from arXiv:1004.2096, with the arXiv licence recorded.
- **Not projected:**
  - **no absolute σ_tot, σ_R or σ_el for ANY variant** (M2, M2′, V-RB, V-M2′ext, or any other reading). Absolute values together with a % deviation, or σ_R + σ_el = σ_tot, would reveal the measured values;
  - measured σ_tot values (M2PRIME's "measured sigma_tot … means" line is dropped);
  - EXFOR points;
  - any % deviation from measurement.

  EXFOR is referenced by entry number only (13569, 13753).
- **Residual inference, disclosed and accepted.** For the non-fitted rows, the projected ratios multiplied by the public Barashenkov values (Geant4 `G4BarashenkovData`, committed in `bgg.py`) rebuild the model cross sections. Together with the pass/fail flags this does not give the measured value. It is accepted as a disclosed residual: the rule forbids redistributing database values, not approximate inference from public tables.
- **Builder-vs-scratch reproduction record ((d)).** The 2 % rule is evaluated on the absolute values at run time. The record is **projected as relative differences only**.
- Amendment 20 revision 7 and decision 0042 keep their reviewed M2 percentages. They are derived quantities, and nothing absolute accompanies them in Git.
- **Scripts.** `m2prime.py`, `m2prime2.py` and `cv_mc.py` are committed outside the rehearsal digest. **They must read the measured EXFOR σ_tot from the cache at run time.** The hard-coded `TOT = {...}` table of `m2prime.py` / `m2prime2.py` (line 52) is replaced by a cache reader before commit, and a test asserts that no EXFOR value is a literal in the committed scripts.
- **Record files** (Git-safe projection; scripts outside the rehearsal digest):
  - `validation/plans/v3-005d-m2prime/M2PRIME-projected.md`
  - `validation/plans/v3-005d-m2prime/scripts/exfor_tot.py`
  - `validation/plans/v3-005d-m2prime/scripts/m2prime.py`
  - `validation/plans/v3-005d-m2prime/scripts/m2prime2.py`
  - `validation/plans/v3-005d-m2prime/scripts/cv_mc.py`
  - `validation/plans/v3-005d-m2prime/scripts/make_projection.py`
  - `validation/plans/v3-005d-m2prime/scripts/olg.py`
