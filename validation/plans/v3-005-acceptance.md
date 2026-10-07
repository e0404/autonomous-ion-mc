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
| P2 | Samplers: Kalbach μ for a ∈ {0.01, 1, 5} × r ∈ {0, 0.5, 1}; Poisson for λ ∈ {0.1, 1, 4}; inverse-CDF E′ on a synthetic 64-bin table; 10⁶ draws each from a fixed CI seed | χ² goodness of fit (≥ 50 equiprobable bins, ≥ 5 expected per bin) with p > 0.001 for every case; the Kalbach first moment within 4 σ | self-consistency | CI | A |
| P3 | Per-event conservation, 10⁵ events from synthetic tables, float64 | max \|ΔE\| ≤ 1e-9 MeV and max \|Δp\| ≤ 1e-9 MeV/c (four-vector of p + target vs products + residual + E\*); ΔZ = ΔA = 0 for every event; no event with m_r < M_r accepted | conservation | CI | A |
| P4 | Thinning harness (shared funcs via twins, through the reference step glue): constant Σ and Σ(E) ∝ E^−0.5 in a synthetic medium; s_max ∈ {0.1, 1} mm; 2·10⁵ particles each | \|S_MC(x) − exp(−∫Σ dx)\| ≤ 0.002 abs at 10 depths, and z = \|S_MC − S\|/σ ≤ 4 | self-consistency | CI | A |
| P5 | Every `make_nuclear` function: python twin vs a Warp-CPU float64 test kernel on recorded inputs (the U1 completeness harness covers every function) | Decision 0001 classes: discrete outputs equal; continuous outputs within 1e-12 relative | backend parity | CI | A |
| N1 | Runtime Σ_mass (uniform ln E grid, ≥ 50 points/decade) vs Σ_ENDF at every ENDF node and every node midpoint, 1–150 MeV; water, ICRU tissue (`materials.py` set) | rel ≤ 1e-3 at every point | self-consistency | LV | A |
| V1 | **Source-model reproduction:** the composed per-element σ (after surrogate scaling and composition) vs ENDF MF3/MT5 per element and for the compositions, 1–150 MeV; plus the extension's continuity at 150 MeV | rel ≤ 1e-3 at every ENDF node; extension: σ(150⁺)/σ(150⁻) − 1 within 1e-6 | self-consistency (source model) | LV | A |
| V1b | **Informative:** σ_nonel of p+C (and p+O where present) vs EXFOR raw entries published 1997 or later. The set is D0356 (Auce et al. 2005) plus the entries in the committed manifest. Weighted mean ratio per window [20,40), [40,70), [70,110), [110,160), [160,250] MeV, with the PDG scale factor. | **Not gating.** Ratios and 1 σ archived. Pre-declared known deficiency: LA150 p+C ≈ 227 mb at 100 MeV vs 245–275 mb measured (up to about −18 %). A window without independent data is recorded as an evidence gap. | measured | LV | A |
| V2 | Primary survival S(z) at 100, 150 and 200 MeV in water, MCS and straggling off, default steps (s_max 1 mm, f_E 0.02); N = 2.4·10⁵ per energy (sharded, see Seeds); 1 mm depth bins | TOST, tol 0.003 abs, vs numpy exp(−∫ρΣ dz) along the deterministic CSDA path (same tables), at every 10 mm depth until the CSDA end minus 5 mm | self-consistency | LV | A |
| V2-probe | 150 MeV, same seed (common random numbers): s_max 0.5 mm vs 1.0 mm, and f_E 0.005 vs 0.02; N = 10⁵ each | \|S_a(z) − S_b(z)\| ≤ 0.002 at every 10 mm depth. The row is valid only if σ_Δ(z) ≤ 7e-4 at the deepest depth; otherwise it is inconclusive, which counts as **fail** and is re-run with more shards and a new base | falsification | LV | A |
| V2b | As V2-probe for s_max 0.1 mm vs 1.0 mm on warp-cpu float64, N = 10⁶ | as V2-probe | falsification | LV | **slice B** |
| V3 | Energy balance with nuclear tallies: `initial = step_deposit + cutoff + nuclear_local + escaped + nuclear_escaped_neutron + nuclear_escaped_gamma + nuclear_binding + truncated + unaccounted`; the shipped per-grid identity of `EnergyBalance` (V3-004; fields `in_grid_mev`, `quantization_mev`, `outside_mev`), `in_grid + quantization + outside = step_deposit + cutoff`, extended with the nuclear-local destination: `in_grid + quantization + outside = step_deposit + cutoff + nuclear_local` per grid; `nuclear_binding` = Σ(Σm_out + M_r − m_p − M_t) recomputed from AME masses per event | relative residual ≤ 1e-12 (python, float64); counters 0. CI: 200 histories at 150 MeV. LV: 2·10⁴ at 150 MeV and 10⁴ at 250 MeV in water / tissue / bone slabs | conservation | CI, LV | A |
| V4 | Event sampler vs the ENDF table: C-12 and O-16 at 100 and 150 MeV, 10⁵ events each (offline, no transport) | Mean yields y_s (s ∈ {n, p, d, α, γ}) within 1 % + 3 σ of ENDF. Per-event Σ y⟨E′_CM⟩ within 5 % of ENDF (C-12 at 150 MeV: 106 MeV). Distortion table archived (accepted attempts, rejection fraction). | self-consistency | LV | A |
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
| P1–P5, V3-CI, V9, R1, C1 | CI run | — | fixed CI seeds | — |
| N1, V1, V1b, D6 | `lv5-<sha>-...json` | — | n/a (deterministic) | — |
| V2 (100/150/200), V2-probe | `lv5-<sha>-...json` | — | 20421004 | — |
| V3-LV, V4, X1, E1 | `lv5-<sha>-...json` | — | 20421004 | — |
| V2b, V5–V8 | slice B | — | declared in B | — |

## Amendments

(none)
