# 0041 — Proton non-elastic nuclear interactions: data, sampling, secondaries and bookkeeping

- Status: accepted (plan 2026-10-03, slice-A brief 2026-10-04, plan delta 2026-10-07)
- Date: 2026-10-07
- Task: V3-005A (data pipeline, shared functions, python reference backend); V3-005B (Warp kernels, generation queue, p-p elastic scattering, kernel parity)
- Affects: physical accuracy, numerical accuracy, reproducibility, architecture, validation strategy, data provenance

## Problem

The engine of decision 0039 transports protons with electromagnetic physics only, and `nuclear=True` raises. Clinical proton dose needs the non-elastic nuclear channel:

- the attenuation of primaries, about 1 % per cm of water;
- the secondary protons and deuterons that it creates;
- the local deposition of heavy fragments;
- the energy that leaves as neutrons and γ.

The requirements are V1-MUST-005, V1-MUST-013 and V1-MUST-014, with fail-closed capability rules. This decision fixes the data source, the sampling algorithm, the secondary model, the bookkeeping and the validation, before any nuclear result exists.

## Context

- **Primary evidence.** The ENDF/B-VIII.0 proton sublibrary was parsed during planning: sha256 `27bcafb8…03f7`, 14 112 244 B. It contains the LA150 evaluation (Chadwick et al., Nucl. Sci. Eng. 131 (1999) 293; library: Brown et al., Nucl. Data Sheets 148 (2018) 1).
  - C-12, N-14, O-16, P-31 and Ca-40 carry only MF1/451, MF3/{2,5} and MF6/{2,5}, with EMAX = 150 MeV. Al-27 is also evaluated.
  - Non-elastic is MF3/**MT5**.
  - MF6 light products (n, p, d, α) use LAW=1, LANG=2 (Kalbach–Mann) with NA=1. Residuals and γ use LANG=1. LCT=3.
  - There are **no t and no ³He products**.
  - H-1 has elastic data only.
- **MT5 values (mb):**

  | | 20 MeV | 100 MeV | 150 MeV |
  |---|---|---|---|
  | C-12 | 441 | 227 | 222 |
  | O-16 | 535 | 297 | 295 |

- **Measured data.** Auce et al. 2005 (EXFOR D0356) gives p+C 279/275/237 mb at 81/100/119 MeV. geant-val (EXFOR-derived) gives 245 ± 7 mb (99.1 MeV) and 275 ± 21 mb (100 MeV).
  - LA150 was fitted to data published before 1997, so it shares lineage with them.
  - The EXFOR Data Explorer API is behind a bot challenge; raw EXFOR master entries are used.
- **Masses.** AME2020 (Wang et al., Chin. Phys. C 45 (2021) 030003), sha256 `e8599c6d…3307`.
- **Engine facts.**
  - Decision 0039: the hinge step, the ramp-apportioned deposit, Gamma straggling and float64 bookkeeping.
  - Decision 0040: the species registry (pseudo-species `nuclear_local`, id 64), the hook `score_piece(…, species, gen, cls)`, and the runtime-checked capacity bounds.
  - Decision 0037: Philox counters `(history, genealogy, block, purpose)`.

## Options considered

| Topic | Options | Selected | Reason |
|---|---|---|---|
| σ_nonel source | ENDF/B-VIII.0 MT5 (LA150); fit to EXFOR; Tripathi-light; Geant4 BGG | **LA150 MT5 ≤ 150 MeV** | Consistent with the MF6 yields used for the secondaries (one evaluation for σ and products). Fitting to EXFOR would decouple σ from the yields and consume the independent evaluation set. |
| 150–250 MeV | TENDL-2023 / JENDL-5 (stop at 200 MeV, licences unverified); constant σ; Tripathi shape | **σ_ENDF(150)·σ_TL(E)/σ_TL(150)**, spectra stretched by E_avail(E)/E_avail(150) | Smooth and anchored. Measured p+C is nearly flat over 158–231 MeV. The fraction of events above 150 MeV is recorded. |
| Unevaluated elements | raise; A-scaling from a neighbour | Na, Mg → Al-27; S, Cl → P-31; K, Ar → Ca-40, scaled by (A/A_ref)^{2/3} and recorded; anything else raises | Tissue mass fractions ≤ ~1 %. Fail closed otherwise. |
| Hydrogen | none; p-p elastic | Non-elastic 0 below the pion threshold. **p-p elastic (D2) in V3-005B.** | p-p elastic is strong-interaction scattering that T_dM does not cover. A claims no lateral halo. |
| Interaction sampling | Σ(E_mid)ρs per step; distance sampling with the inverse range; **majorant thinning** | **Majorant thinning with remaining mean free paths**, candidate at the step end | Exact for any step size given Σ ≤ Σ̂. It composes with the hinge and the step limits, and needs no position sampling. |
| Products | parameterised (Fippel–Soukup style); Geant4-derived event tables; **ENDF MF6 event model** | **ENDF MF6**: Poisson multiplicities (λ adjusted offline), inverse-CDF E′_CM, exact Kalbach μ, residual by four-momentum difference | Uses only the construction data. Gives conservation per event by construction. Independent of the reference engines. |
| α | transport (z = 2 machinery, V3-008); local | **Local deposition**, subject to the two-tier gate D6 (below) | Mean α energy 4.5–6 MeV; range tens of µm. |
| Secondaries' nuclear interactions | include; omit | **Omitted in A** (declared, below) | No deuteron data. Secondary-proton interactions are moved to V3-005B with an estimate. |

## Decision

### 1. Cross sections

- **Per-element σ(E).** For C-12 (natural C; C-13 at 1.1 % neglected), N-14, O-16, Al-27, P-31 and Ca-40, σ(E) comes from MF3/MT5. Surrogates are as above.
- **Hydrogen** contributes 0.
- **Runtime rows.** The builder writes per-material rows on a uniform ln E grid with ≥ 50 points per decade (*amended 2026-10-07, C7 outcome:* on the union grid of all ENDF MF3/MT5 and MF6/MT5 energies in [1,150] MeV with a uniform ln E grid ≥ 50 points/decade over [1,250] MeV, interpolated lin-lin in E, looked up by the fixed-step binary search `grid_locate`; the builder fails if any MF3/MT5 or MF6 yield TAB1 has an interpolation law other than INT=2 in [1,150] MeV; where σ > 0 while the evaluation tabulates zero yields (Al-27, P-31 at 1 MeV; Ca-40 at 3.44 MeV) every product's yield is held at its value at the first energy with positive total yield (2.0 / 1.307 / 4.0 MeV) down to the MT5 threshold, recorded per target in the JSON (`yield_extended_below_mev`); see acceptance-plan Amendments 3 and 4):
  - Σ_mass(E) = N_A Σ_el w_el σ_el/A_el [cm²/g];
  - cumulative target fractions;
  - the majorant Σ̂(E) = 1.02 × max Σ over [E(1 − 2f_E − 0.01), E], or over [0, E] for an end-of-range step.
- **Linear cross section.** The linear Σ is ρ_voxel Σ_mass, scaled with density like the stopping tables.

### 2. Sampling (thinning)

- **Birth.** A primary samples n_λ = −ln u at birth.
- **Step limit.** `select_step_nuclear` adds the limit d_nuc = n_λ/(ρΣ̂(E₀)) as reason 4; geometry still wins ties.
- **After the unchanged EM step** (hinge, straggling, ramp deposit):
  - n_λ is reduced by ρΣ̂(E₀)·s_act.
  - If the step was nuclear-limited, leg 2 was not truncated, and E₁ > E_cut, a candidate occurs at the post-step point.
  - The candidate is accepted with probability Σ(E₁)/Σ̂(E₀). Otherwise it is fictitious and n_λ is resampled.
- **Fail closed.** Σ(E₁) > Σ̂(E₀) increments `majorant_violation` and invalidates the result.
- **Window factor.** *(corrected 2026-10-08, Codex REVIEW-c5c9149b finding 5)* The factor 1.02 widens the ratified window [E(1 − 2f_E − 0.01), E] (derived for the clamped Gaussian) so that the bulk of the Gamma straggling distribution lies inside it; no finite window covers the unbounded Gamma tail. Exactness of the thinning therefore rests on the fail-closed rule below, which is evaluated after every step of the primary, not only at candidates: any step whose post-step Σ exceeds the majorant used for its optical-depth decrement invalidates the result.
- **Termination.** The parent proton is terminated at an accepted event; the yields include the leading proton.

### 3. Event model

- **Target and rows.** Select the target from the cumulative fractions. Interpolate its table rows at E₁ (above 150 MeV: the stretched 150 MeV rows).
- **Multiplicities.** Independent Poisson(λ_s) for s ∈ {n, p, d, α, γ}, with n ≤ 16. λ_s is fixed-point adjusted at build time, so that the post-acceptance mean yields match ENDF within 1 %.
- **Kinematics.**
  - Draw E′_CM by inverse CDF (64 equiprobable bins, linear within a bin).
  - Draw the Kalbach μ with a from the Kalbach (1988) systematics (NA=1) and r from MF6; φ is uniform.
  - LCT=3 is taken as the p + target CM.
- **Residual.**
  - (Z_r, A_r) = (Z_t + 1 − Σz, A_t + 1 − Σa).
  - Its four-momentum is P_tot − Σp_i.
  - Accept iff the nuclide has an AME2020 mass and m_r ≥ M_r. Then E* = m_r − M_r.
- **Attempts.** At most 64 attempts. Exhaustion increments `nuclear_rejection_limit`, tallies the energy as `unaccounted` and invalidates the result.
- **Lab frame.** Boost everything with the p + target CM velocity. Energy, momentum, charge and baryon number are then conserved exactly per event, to floating-point rounding.
- **Disposition** (the α row and the residual row are read with the amendment below):

  | Product | Treatment |
  |---|---|
  | p and d | Transported as generation-1 secondaries from a per-history LIFO stack (at most 32 particles per history). Each carries position, direction, T, species, genealogy id, generation, the parent's geometry voxel indices (never recomputed by floor) and p₁v₁ = pv(T). |
  | α, residual kinetic energy, E* | Deposited at the event point through the scoring hook with `length = 0`, species `nuclear_local` (id 64), generation 1 and class local. They therefore appear in dose and in `edep_excluded_from_let`, never in LET, fluence or lookups. |
  | n and γ | Tallied as `nuclear_escaped_neutron` / `nuclear_escaped_gamma`. |
  | Mass difference | Σm_out + M_r − m_p − M_t (from masses, signed) is tallied as `nuclear_binding`. |

- **Deuterons.**
  - They use their own `TransportTables` from `stopping.table(material, DEUTERON)` (Bethe with q = z = 1), on the total-kinetic-energy axis, and a deuteron water row for LET.
  - E_cut,d = 4 MeV.
  - The EM twins take the deuteron mass and charge.
  - Deuterons have no nuclear interactions.
  - `nist-star` with `nuclear=True` raises (no deuteron table).

*Amendment 2026-10-07 (C7 outcome).* The event model above failed V4 before qualification (acceptance-plan Amendment 1 quotes the numbers verbatim; the C7 table `99b51c2c…` is archived as superseded evidence). It is replaced by: multiplicities n_s = ⌊λ_s⌋ + [u_s < λ_s − ⌊λ_s⌋] (shared function `multiplicity_round`, cap 16 kept, same RNG slots) with λ solved at build time by exact enumeration of the ≤ 32 outcomes so that the post-acceptance mean yields equal ENDF (tolerance 1e-6, nodes not reaching 1e-3 recorded `converged=false`); acceptance iff the residual (Z_r, A_r) has an AME2020 mass or is the empty residual (0, 0) with M_r = 0 (full break-up into light products, e.g. p + C-12 → p + 3α; added 2026-10-07 after the first C7b build, acceptance-plan Amendment 4) (≤ 64 attempts, exhaustion rule unchanged; the builder fails if the exact P_accept < 0.5 at any node or midpoint, so the runtime exhaustion probability is ≤ 2⁻⁶⁴); products sampled independently as before; the residual receives the ENDF mean heavy-recoil energy Σ_r y_r⟨E_r⟩ (formally CM under LCT=3, deposited without a boost; lin-lin in E, stretched above 150 MeV) and is deposited locally with the α; E* is not deposited (residual de-excitation is represented by the ENDF γ yields); the per-event difference Δ = T₁ + m_p + M_t − Σ E_lab − M_r − T_r is tallied as `nuclear_imbalance` (signed, never deposited), so that T₁ = Σ T_lab + T_r + nuclear_binding + Δ exactly. Energy and momentum are conserved on average to the extent the evaluation is, not per event; single histories may deposit more than their initial energy. Reason: with inclusive LA150 yields and spectra any per-event-conserving scheme removes 34–64 % of the product energy and the by-difference residual recoil is 5–10× the evaluated one (planner experiments 2026-10-07). Rejected: a momentum-balancing residual (negative lab kinetic energy for the typical net CM momentum); a per-node product-energy scale κ (breaks the V4 5 % criterion; pre-registered fallback if V5 fails); 256 inverse-CDF bins (×4 memory). The "Products" row of the options table is read with this amendment.

### 4. Random numbers (amends decision 0037)

- The purpose table becomes: 0 transport (EM), 1 source sampling, **2 nuclear** (`PURPOSE_NUCLEAR`; formerly reserved). The name `PURPOSE_RESERVED` is kept as an alias of the same value 2, so existing imports keep working.
- Each primary has a nuclear block counter `nuc_c2` on purpose 2, with this fixed layout:
  - birth `(u_nλ, ·, ·, ·)`;
  - candidate `(u_accept, u_nλ, u_target, ·)`;
  - event attempt: 2 multiplicity blocks, then per particle `(u_bin, u_frac, u_branch, u_μ)` and `(u_φ, ·, ·, ·)`.
- Secondaries use purpose 0 with their own genealogy id (decision 0037 encoding: child b of generation g gets parent + b·32^g).
- EM streams are therefore identical with nuclear on and off, which gives low-variance paired differences.

### 5. Bookkeeping and fail-closed rules

- **Balance.** `initial = step_deposit + cutoff + nuclear_local + escaped + nuclear_escaped_neutron + nuclear_escaped_gamma + nuclear_binding + truncated + unaccounted` (*amended 2026-10-07:* `+ nuclear_imbalance`, signed, never deposited; the conditional block also carries `nuclear_alpha_local`, the α part of `nuclear_local`).
- **Grid identity.** `in_grid + quantization + outside = step_deposit + cutoff + nuclear_local`.
- **Layout.**
  - The nuclear tallies and the counters `majorant_violation`, `nuclear_rejection_limit` and `nuclear_conservation` form blocks that exist only with `nuclear=True`.
  - `TALLY_NAMES`, `COUNTER_NAMES` and the trace columns are unchanged.
  - The existing `genealogy_overflow` and `queue_overflow` guard the stack.
- **Capacity bounds** (re-derivation required by decision 0040 section 4):
  - The per-history path bound becomes B_L = 1.25 · mixed_path_bound(E_hi + Q⁺). Q⁺ is the largest positive Q of the reachable channels.
  - The argument: superadditivity of the convex CSDA range, plus R_d ≤ R_p at equal total energy.
  - The scored path sums over all particles of a history.
  - The E bound is E_hi + Q⁺. The piece-count bound is multiplied by the particle bound 32.
  - *Amended 2026-10-07 (no per-event conservation; second form after the first C7b build, Amendment 4):* per species s ∈ {p, d} the builder records N_s,max = max over grid nodes of ⌈λ_s⌉ (≤ 16) and T_lab,max,s (top bin edge boosted with μ = 1) in the table JSON (`transport_path_bound_terms`), and `transport_energy_bound_mev` = max_s T_lab,max,s (per-particle bound, asserted ≤ the stopping-table maximum). B_L = 1.25 · (mixed_path_bound(E_hi) + Σ_s N_s,max · mixed_path_bound_s(T_lab,max,s)). Two energy bounds (clarified 2026-10-08 after Codex review REVIEW-d8700e0e, finding 1): the per-particle bound `transport_energy_bound_mev` = max over transported species of T_lab,max,s governs the stopping-table coverage; the per-history energy bound that sizes the fixed-point accumulators (E/ES/FE) is `history_energy_bound_mev` = max(E_hi, Σ over ALL product species s ∈ {n, p, d, α, γ} of N_s,max · T_lab,max,s + T_r,max), because per-event conservation is withdrawn and a history may deposit more than its initial energy; the builder writes both. The earlier single-term form R(T_sum,max) with T_sum,max = 1030 MeV would have required extrapolating the stopping table beyond 500 MeV and is withdrawn.
  - `path_bound_exceeded` still backstops these at runtime.
- **Rejected before transport** (`validate()`):
  - an unsupported element;
  - E₀ + 6σ_E > 250 MeV; *(added 2026-10-08, Codex REVIEW-c5c9149b finding 1)* at runtime every sampled source energy above 250 MeV increments `source_energy_out_of_range` and invalidates the result, because the Gaussian source is unbounded;
  - a non-proton source;
  - `nist-star`;
  - any warp backend (until V3-005B);
  - a missing, stale or mis-pinned table (*2026-10-08, Codex REVIEW-c5c9149b finding 2:* the qualification flags, convergence flags and capacity bounds are stored inside the npz; *2026-10-08, REVIEW-7d0fb725 finding 1:* the table id is the sha256 of the canonical JSON of the whole sidecar without the id field, which contains `npz_sha256`, so the pinned id authenticates every sidecar field (element rows, target Z/A, gate numbers) and the npz bytes; the loader recomputes the npz digest and the id before using anything);
  - *(added 2026-10-08, Codex finding 3)* an unqualified table: the loader refuses a table whose JSON records `gate_d6.ceiling_pass` false or any non-converged λ node (`NuclearTableUnqualifiedError`); the builder still writes such tables as evidence but they cannot be transported with.
- **Producible species.** With `nuclear=True` the producible set is {(proton, primary), (proton, secondary), (deuteron, secondary), (nuclear_local, secondary)}.

### 6. Data pipeline and roles

- `ionmc data build nuclear-proton` reads the hash-pinned sources from the cache. It writes `<cache>/derived/nuclear-proton-<id>.npz` plus JSON. *(Amended 2026-10-08, Codex REVIEW-7d0fb725 / REVIEW-3a1a1c8b:)* id = sha256 of the canonical JSON (sorted keys, compact separators) of the whole JSON sidecar without the `table_id` field; the sidecar carries `npz_sha256`, the source hashes, the builder version and the canonical options, so the id authenticates the npz bytes and every sidecar field. The original formula sha256(source hashes, builder version, canonical options) is superseded.
- The JSON records surrogates, the extension method, the λ residuals, the gate values G and D, and units.
- Loading uses `np.load(allow_pickle=False)`, re-hashes the npz, recomputes the id from the sidecar, checks the source pin, derives the qualification flags and capacity bounds from the npz arrays (refusing an unqualified table), cross-checks the sidecar, and freezes the arrays.
- Git holds no raw or derived nuclear data. It holds only the registry entries, the EXFOR manifest (IDs, roles, hashes; no values), synthetic fixtures and aggregate statistics.

| Data | Role | Licence / citation |
|---|---|---|
| ENDF/B-VIII.0 protons (LA150) MF3/MF6, sha256 `27bcafb89cf0444c53c6b9f3dd17618c62e2e6b3694f70d31c4f502399f103f7` | construction | No licence text (free NNDC download); cite Brown 2018, Chadwick 1999 |
| AME2020 `mass_1.mas20.txt`, sha256 `e8599c6d7f724fac91934e59f1b9de8fb8f63e820f4b39456b790665ed2a3307` | construction | Citation only; Wang et al. 2021 |
| Tripathi light-system parameterisation (NASA TP-1999-209726) | construction (shape above 150 MeV); report-only theory ≤ 150 MeV | US Government work |
| NIST ASTAR water | construction of the D6 gate only | NIST SRD, public |
| EXFOR D0356 (Auce et al., PRC 71 (2005) 064606; p+C, p+Ca), sha256 `2ef17fb10aaaeb5096dc92f9c3175f51f2252163eb4a83aa15669027371747bf`, and further raw entries published 1997 or later in the manifest | evaluation, informative (V1b) | CC BY 4.0; cite the experiment and NDS 120 (2014) 272 |
| EXFOR C1862 (Slaus et al., PRC 12 (1975) 1093; p+C, p+Be, p+O), sha256 `157aca7e1b8f6cb6814a5fa99fa749823721af6829cbf84712a56fce159ae4ef` | exploratory, report-only (published before 1997: shared lineage with LA150) | CC BY 4.0; cite the experiment and NDS 120 (2014) 272 |
| geant-val EXFOR-derived inelastic data (−7), sha256 `fa7ac90fd259f728e5948c71b7a3636bec7b7d9756abeb60837c7083dc2a4c49` | exploratory, report-only (shared lineage with LA150) | no licence stated |
| MCsquare (ICRU 63 = LA150 lineage) | shared lineage for σ, not independent | — |
| TOPAS QGSP_BIC_HP | independent MC (slice B gating, slice A exploratory) | — |

### 7. α local-deposition gate D6

Definitions. G is the share of the MF6/MT5 α energy (water-weighted, lab frame) carried by α with ASTAR water CSDA range > 0.1 mm, at 150 and 250 MeV. D = (α energy carried by α with range > 0.1 mm per event) × P_event(E) / E, where P_event is the non-elastic probability over the full CSDA path in water at that beam energy, computed from the built Σ along the deterministic CSDA path. The builder computes G and D and writes both to the table JSON; both are archived.

Both tiers below are fixed now, before any G or D value exists.

- **Tier 1 (as ratified):** G(150) ≤ 0.05 and G(250) ≤ 0.05. Then α are deposited locally with no further condition.
- **Tier 2 (fallback, only if tier 1 fails):** local deposition is still accepted iff D ≤ 1e-3 AND the 99.9th-percentile lab α energy has an ASTAR water CSDA range ≤ 2 mm at both energies. Passing via tier 2 is recorded in the Outcome as a documented approximation. It forbids sub-millimetre statements about α dose until V3-008 (z = 2 transport), and the ledger status says so.
- **Neither tier:** the task stops. The orchestrator decides between blocking α-dependent claims and re-scoping. There is no α transport in slice A.

*Outcome and amendment 2026-10-07.* Neither tier passed on the C7 table (150 MeV: G 0.5718, D 3.78e-3, 99.9th-percentile α range 0.258 g/cm²; 250 MeV: G 0.8038, D 1.10e-2, 0.852 g/cm²; the planner's estimate with the amended sampler: G 0.672 / 0.850, D 5.55e-3 / 1.42e-2, 0.995 / 2.51 g/cm²). Orchestrator decision: local deposition is kept as a documented approximation, valid only if the recomputed D ≤ 2e-2 and the 99.9th-percentile α range ≤ 3 g/cm² at both energies (ceiling set with the estimate known); every nuclear run reports `nuclear_alpha_local`; no statement about α dose, α LET or nuclear-secondary dose structure below 1 cm (150 MeV) / 2.5 cm (250 MeV) until V3-008 re-runs D6 with transported α. Straight-line CSDA α placement was rejected as a re-scope duplicating V3-008. See acceptance-plan Amendment 2.

### 8. Amendment to decision 0039

The row "Nuclear flag: `nuclear=True` raises until V3-005" becomes:

> `nuclear=True` is supported on the python backend from V3-005A (this decision). Warp backends raise until V3-005B.

The nuclear step limit is reason 4 of `select_step_nuclear`; `select_step` (reasons 0–3) is unchanged.

### 9. Multiprocessing

The python pool path is implemented for `nuclear=True` like any other run: the per-history stack is history-local and the partitioning is unchanged. The current operator directive (single-process execution) is an execution mode of the validation, not a property of the engine, and is not baked into `validate()`. The lv5 suite runs `cpu_workers = 1`; the 1-vs-N partition-invariance check of nuclear runs is listed as **deferred** in the acceptance plan until the directive is lifted.

## Approximations and known limitations

1. **LA150 p+C deficiency** (pre-declared before any comparison). σ_LA150(p+C, 100 MeV) = 227 mb, against 245–275 mb measured (Auce 2005; geant-val), i.e. up to about −18 %.
   - At about −10 % in σ_C, the effect on water attenuation is small, because O dominates (water σ is about 89 % O by mass).
   - The p+O evaluation has almost no independent data between 65 and 250 MeV (evidence gap).
   - *Planned remedy:* V3-005B evaluates Tripathi-light and the Geant4 BIC/BGG σ_inel as an alternative data role against the post-1997 EXFOR set. Pre-registered fallback D9: refit on pre-1997 entries only, evaluate on the unchanged post-1997 set.
2. **No t and ³He** (absent in LA150; physically about 0.1–0.2 per event [M]). TOPAS comparisons will show the gap.
3. **α deposited locally.** *(amended 2026-10-07)* Both D6 tiers failed (section 7); local deposition is kept by orchestrator decision as a documented approximation under the D ≤ 2e-2 / 3 g/cm² ceiling, with the claim restrictions of section 7 until V3-008. About 0.6 % (150 MeV) / 1.4 % (250 MeV) of the beam energy in water is α energy deposited at the interaction point instead of along α tracks of up to ≈ 1 / 2.5 cm (99.9th percentile).
4. *(amended 2026-10-07)* **No per-event energy or momentum conservation.** Planner estimate of the mean `nuclear_imbalance` per event at 100/150 MeV: C-12 −6.5/−10.5 MeV, O-16 −8.2/−13.7 MeV (products carry more energy than available; 7–10 % of E_avail; sd 50–78 MeV); in water at 150 MeV about 1.5 % of the beam energy appears in non-elastic products (≈ 65 % of it in transported p/d). Every run reports the sum. E* is not deposited; the residual carries only the ENDF mean recoil.
5. **Uncorrelated floor+Bernoulli multiplicities** preserve mean yields only (minimum-variance choice), not fluctuations or correlations. Microdosimetry is out of scope.
6. **No nuclear interactions of secondaries.** Estimate [M]: a secondary proton (mean about 40 MeV) undergoes a non-elastic event with P ≈ 1.5 %. That is ≈ 0.3 % of 150 MeV histories, carrying ≲ 0.1 % of the energy. V3-005B quantifies it.
7. **No p-p elastic and no nuclear elastic scattering** on C/O. No lateral-halo claim is made until V3-005B (V6).
8. **LCT=3** is interpreted as the p + target CM. *(amended 2026-10-07)* E′_CM is uniform within each of 64 equiprobable bins: ⟨E′⟩ is biased by +3…+5 % for α and +14…+24 % for γ at 100–150 MeV (n, p, d within 0.2 %); the γ bias affects only the escaped-γ tally. Residual-existence rejection changes yields by < 1 % after the exact λ correction (V4 records it).
9. **Kalbach separation energies** use the Kalbach (1988) systematics formula, not AME2020 masses (decision at C4, 2026-10-06); the a parameter is insensitive to this at the per-mille level.

## Validation strategy (frozen before results exist)

`validation/plans/v3-005-acceptance.md` is committed with this decision.

- **Slice A rows:** P1–P5, N1, V1, V1b (informative), V2, V2-probe, V3, V4, V9, R1 (= A16 re-pointed to f3a1dd62, empty intended-change set), X1, C1, D6 (two tiers), E1 (exploratory).
- **Slice B rows:** V2b and V5–V8.
- **Gating physics evidence:** the depth-dose and attenuation comparison with the reference engines (V5, the V3-010B batches) and the measured comparisons of V5/V6 and V3-011. V1b is informative only.

## Outcome

- 2026-10-07, C7 table `99b51c2c9001bb6f41a87e2092ca70172c14ee2d207f5e4774ff7460409dcdd1` (superseded): D6 neither tier (section 7); N1/V1 failed on the uniform grid (Amendment 3); the λ fixed point did not converge for C-12/N-14 at any node (Amendment 1). All three led to dated amendments before any transport result existed.
- 2026-10-07, C7b table `2b8d94cf4a82ef8a5119ddb8c84d4a2c11bd3511eb6672bd210562730cc0cf1a` (builder/schema v2, 604 nodes, 12.8 s): D6 recomputed — 150 MeV G 0.6680, D 5.481e-3, 99.9th-percentile α 127.1 MeV / 0.988 g/cm² (energy-weighted 1.109); 250 MeV G 0.8492, D 1.420e-2, 213.9 MeV / 2.529 g/cm² (energy-weighted 2.827); tiers 1 and 2 not passed; the section-7 ceiling (D ≤ 2e-2, range ≤ 3 g/cm²) holds, so local α deposition stands as the declared approximation. N1 2.8e-15, V1 2.2e-16, continuity 1.5e-10. V4 pass; V4b pass under Amendment 4 (raw n ratio 1.0122 for O-16 at 100 MeV, sem 0.0036). Mean `nuclear_imbalance` per event −6.6/−10.8 MeV (C-12) and −8.2/−14.2 MeV (O-16) at 100/150 MeV. Capacity terms: per-particle bound 262.08 MeV, B_L(250 MeV, water) 2489 mm.
- 2026-10-08: after Codex review REVIEW-d8700e0e (six findings fixed: per-history energy bound, majorant check on every step, loader refuses unqualified tables, independent AME recomputation of `nuclear_binding`, hash-bound shard partials, consistent capability report) the table was rebuilt with builder version 3 as `dfee19d303c7fdd28d2080ebbb37a01b27aa6569cca846ad317a2d91a7a437a3` (604 nodes, 11.6 s; same physics and gate numbers as `2b8d94cf…`, which is superseded); `history_energy_bound_mev` 3154.4 MeV (2 n × 253.3 + 3 p × 262.1 + 1 d × 246.2 + 3 α × 257.4 + 3 γ × 278.9 + recoil 6.3), per-particle bound 262.08 MeV; loader qualification: D6 ceiling true, 0 non-converged nodes.
- 2026-10-08: after Codex review REVIEW-c5c9149b (five findings fixed: runtime rejection of source energies above 250 MeV, table id covering the npz digest with qualification flags and bounds stored in the npz, orchestrator-attested shard manifests, consistent capability report, corrected window-factor text) the table was rebuilt with builder version 4 / schema 3 as `3bcf146e38dd2b5581bd1ff245c127a7d059e6789421a5d3c6760974f2784504` (604 nodes, 11.7 s; identical physics and gate numbers; `dfee19d3…` and `2b8d94cf…` superseded and now refused as stale).
- 2026-10-08: after Codex review REVIEW-7d0fb725 (four findings fixed: whole sidecar bound to the table id, attested partials manifest with host-runner run ids and an attestation block in the combined summary, data-layer docs refreshed, capability wording) the table was rebuilt with builder version 5 / schema 3 as `00e8031d5f67704f256e904ebab89b0c2a75845bcf1887c29e1518bffaafe5fb` (identical physics and gate numbers; all earlier ids superseded and refused as stale).
- 2026-10-08, lv5 at head `b84fdf389e8bae3b3820b73497086d071c3325ad` (single-process diagnostic mode, seed base 20421004, VAL-20261008-064717-C8A3EE, archive `validation/results/transport/lv5-b84fdf3-single-process.json`): V2 pooled 2.4e5 histories per energy — max |ΔS| 5.4e-4 (100 MeV), 8.3e-4 (150), 8.7e-4 (200), pass; V2-probe conclusive (max |Δ| 5.1e-5 / 2.4e-5); V3-LV balance residual ≤ 6.2e-16 with the independent AME recomputation of `nuclear_binding` within 6.2e-9 MeV on six runs; V4/V4b pass at 1e5 events per case; X1 pass (all fail-closed counters 0); E1 exploratory 150 MeV: peak 156.75 mm, R80 158.46 mm, 80–20 fall-off 2.27 mm, mean `nuclear_imbalance` −1.05 MeV per history (0.7 % of the beam energy appears in products, inside the planner's estimate), α deposited locally 0.88 MeV per history; R1 field-by-field identity with f3a1dd62. Every executed row passes; the summary is non-conformant by code (single-process directive, deferred partition row, imported partials). Contrary evidence carried: D6 tiers 1 and 2 failed (ceiling holds); the 64-bin ⟨E′⟩ biases; the mean energy excess of the inclusive sampling; the LA150 p+C deficiency (V1b 0.894 at 70–110 MeV).
- Open for slice B (V3-005B): warp kernels and parity, TOPAS comparison V5/V6 (deciding test for the inclusive-sampling energy excess; pre-registered κ fallback), secondary-proton nuclear interactions, multiprocess rows once the directive is lifted.

## Slice B (V3-005B) — ratified 2026-10-08, before any slice-B result

Acceptance-plan Amendments 6 (seeds) and 7 (scope and decision rules) are the binding record. In brief: Warp CPU/CUDA nuclear kernels keep the per-history LIFO stack of 32 in a per-thread scratch array (CUDA chunk ≤ 2^14) in a separate kernel variant so that the nuclear=False path stays textually unchanged; the event sampler becomes a `@wp.func` built from the existing twins with event-level parity on recorded inputs (P5 extended); V8 trajectory parity python vs warp-cpu f64 and statistical parity including CUDA; secondary-proton non-elastic interactions added after V8 on the frozen slice-A model (intended-change record); V5 gated by TOPAS with the frozen κ fallback rule; V2b, V7, D9 (report-only BGG). p-p elastic (D2), V6 (Gottschalk) and the exploratory E1-B are re-assigned to V3-005C. Planner evidence (python, water, 150 / 200 MeV, single process): deposited energy 96.19 / 94.37 % of E₀; escaped n 1.42 / 2.35 %, escaped γ 0.25 / 0.24 %, binding 2.84 / 3.68 %, imbalance −0.72 / −0.65 %; P_event 0.160 / 0.239; secondary protons 0.231 / 0.390 per history (mean 26 / 34 MeV); expected secondary-proton non-elastic events 0.0037 / 0.0096 per history (0.12 / 0.32 % of E₀); secondary step work 12 / 19 % of primary step work (CUDA SIMT efficiency ≈ 0.70 with the in-kernel stack; a two-wave generation queue would be only 4–17 % cheaper).

*Addendum 2026-10-08 (slice B, Codex REVIEW-7516c131):* secondaries above the 250 MeV table endpoint are transported without nuclear sampling until inside the domain (no extrapolation); the V5 comparator takes every engine configuration from the fingerprint-verified normalized input; the BGG part of D9 is deferred to V3-005C (plan Amendment 10).

*Addendum 2026-10-08 (slice B, V7 replicates):* the lv5b qualification base 20441004 is consumed by a failed V7 replicate-coverage step (0.644 against 0.68 ± 0.03; cluster-aware 90 % interval [0.596, 0.692] from a bitwise sandbox reproduction; a second diagnostic at a non-slice base gave 0.622; both below the nominal 0.670, recorded as unexplained); the suite is re-run at base 20471004 with 90 replicates and the step records replicate-level statistics; tolerances unchanged (plan Amendment 11, revised after Codex REVIEW-db4474cd before the re-run). The consumed run's V5 verdict (fail, attributed to the absent hadronic elastic channel) is preserved.

*Addendum 2026-10-08 (slice B, V7 rule):* before any V7 observation at base 20471004, row V7 is judged by a replicate-level equivalence test (300 replicates; 90 % t interval of the mean replicate coverage within [0.640, 0.700], centred on the exact nominal 0.670 of the implemented 20-batch interval rule; false acceptance ≤ 0.05 at the margins, power ≈ 0.91 at the nominal); the pooled point gate 0.68 ± 0.03 is reported only (plan Amendment 12, Codex REVIEW-d3f216ca).

*Addendum 2026-10-08 (slice B, V7 rule, final):* row V7 covers all three frozen estimators (sec_p and nuclear_local as 12-bin profiles, escaped_neutral as a scalar) with 7200 replicates from eight sharded 9e6-history runs and two gates per estimator: (1) a reference-free paired gate — 3600 fixed disjoint pairs, hit iff |Δmean| ≤ √(sem_j² + sem_k²) with a zero standard error a miss, nominal c₀ = P(|t₃₈| ≤ 1) = 0.6764, region c₀ ± 0.03, exact Clopper–Pearson/empirical-Bernstein bounds at α 0.04, finite-sample and assumption-free; (2) a single-interval gate — each of the 5400 evaluation replicates' reported mean ± standard error against the held-out mean of shards 6–7 swept over a t/bootstrap confidence box on that reference (exact per-bin worst case, monotonicity argument), nominal c₁ = P(|t₁₉| ≤ 1) = 0.6701, region c₁ ± 0.03, Clopper–Pearson (scalar) / fixed Hoeffding radius 0.0173 (profiles) at α 0.04 — the box level assumes Gaussian marginals of a mean of 36 000 blocks, bootstrap-checked; profile bins fixed from the independent 1e6 run by a resolvability criterion (expected per-replicate relative standard error < 0.5). Calibrated by an exact-SHA step, per gate (Amendment 13): paired gate — false acceptance exactly 0.0386 / 0.0387 at the two margins for the scalar (≤ 0.04 by construction) and 0 (upper bounds 0.005) for the Gaussian profiles; single-interval gate — ≤ 0.020 (scalar) and ≤ 0.001 (profiles); row — ≤ 0.004 at every margin; Gamma-block and zero-inflated surrogates ≤ 0.007; power at the nominal for Gaussian-like tallies 0.96–1.00 (paired gate), 0.97–1.00 (single gate), 0.95–1.00 (row); for strongly skewed Gamma surrogates the single gate's power is only 0.14 (a stated limitation of the distribution-free radius); row passes iff all three estimators pass both gates (plan Amendment 13, Codex REVIEW-be30e621/873ab9cd/edb9970b/f1b32614/852837a3/22711e6f/97a3377e/25b77cf0/0af7808c/219d9e08).

### Slice B outcome

- 2026-10-09, lv5b at head `7aae5bb126924c27f0a61622384b23623c2877e3` (Codex REVIEW-a9f890415d8b4fc0a9e9774b622a96b7 passed; single-process diagnostic mode, seed base 20471004, archive `validation/results/transport/lv5b-7aae5bb-single-process.json`, run ids in `validation/results/transport/README.md`; validation records keyed to the evidence commit that follows, status `failed`, listed by `inspect_local_validation V3-005B`). **Two rows fail their pre-registered criteria (V5, V7); every other executed row passes.** The summary is non-conformant by code (single-process directive; imported partials in steps 08, 11, 24; deferred partition step 25).
  - **V8-LV** (step 02): python vs warp-cpu f64 with nuclear on, 256 histories at 150 MeV and 96 nuclear-dense histories at 100 MeV (σ × 40; 185 events, 279 secondaries in total): identical events, species and genealogy, maximum continuous difference 0.0 (tolerance 1e-10), counters 0 — pass.
  - **R1-nuc** (step 03): A16 `nuclear=False` regression against f3a1dd62 bit-identical with an empty intended-change set, T1 pass; python nuclear-on primary-only digest equal to the b84fdf38 digest; the secondary-proton intended-change record (C13) matches — pass.
  - **V5** (step 08, TOPAS gating): **fail.** peak/plateau +12.5 % (150 MeV) / +12.3 % (200 MeV) against 2 %; ΔIDD plateau integral (on − off) −18.5 % / −14.7 % against 10 %; pass on plateau (−0.97 % / −1.56 %), R80 (−0.18 / +0.01 mm) and total deposit (−0.30 % / −0.54 %). κ is not used (total deposit passes at both energies, Amendment 7(c)). MCsquare (report-only) also fails peak/plateau and ΔIDD and, unlike TOPAS, plateau (−1.8 % / −1.7 %) and R80 (+0.48 / +0.64 mm); total deposit agrees. *Attribution (pre-declared, Amendments 7(d) and 11; limitation 7):* the missing hadronic elastic channel (p-p and p-nucleus), present in the TOPAS physics list; the exploratory TOPAS cases X1–X4 (`validation/reference_cases/README.md`, not evidence) were built to test it. Outcomes at 150 MeV, peak relative to the seed-averaged TOPAS EM-only run (full 0.769, per-seed 0.766–0.771): X1 (full minus hadronic elastic) 0.854, X2 (EM-only plus elastic) 0.909, product X1×X2 0.776, controls X3 0.768 and X4 0.770; elastic scattering is a large but not the only contributor (table and run ids in `validation/reference_cases/README.md`, single seed each, not evidence). The consumed run at base 20441004 (b40d8121, independent samples) gave the same verdict pattern (peak/plateau +12.6 % / +11.9 %, ΔIDD −20.8 % / −13.4 %, total deposit −0.29 % / −0.58 %). The attribution is a hypothesis for V3-005C to test, not a pass.
  - **V2b** (steps 09–11): s_max 0.1 vs 1.0 mm on warp-cpu f64, 1e6 per variant: max |Δ| 8.7e-6 against 2e-3, σ_Δ at the deepest depth 1.4e-5 (σ_max 7e-4), conclusive — pass.
  - **V7 N-scan and f32/f64 on warp-cpu** (step 12): relative-SE slopes −0.495 (escaped neutral), −0.518 (`nuclear_local`), −0.480 (secondary-p dose) within −0.5 ± 0.05; f32 vs f64 |z| ≤ 0.68 — pass. **V7 grid shift / refinement** (step 13): 5.4e-4 (≤ 1e-3) and 1.4e-13 (≤ 1e-6) — pass.
  - **V7 replicate coverage** (steps 14–24, Amendment 13): **row fails.** `sec_p` passes both gates (paired 0.6723, bounds [0.6626, 0.6821]; single 0.671, [0.6499, 0.6909]; 11 bins, bin 11 excluded with r_b 1.14); `nuclear_local` passes both (paired 0.6756, [0.6665, 0.6848]; single 0.6664, [0.6435, 0.6868]; 11 bins, bin 11 excluded with r_b 1.19); `escaped_neutral` **fails the paired gate** (m 0.6586, Clopper–Pearson lower bound 0.6445 < 0.6464) while its single-interval gate passes (0.6557, [0.6420, 0.6676]). No degenerate pair-bins; the legacy point gate passes for all three; the calibration step 23 passed. Per Amendment 11 no further base is taken: the V2-NUM uncertainty-coverage evidence of row V7 is not established. *Hypothesis (not a finding):* heavy-tailed 500-history block sums of the escaped neutral energy make the 20-block standard error under-cover by about 1.7 % for this scalar (the skewness of its 1e4-history replicate means is only 0.05).
  - **Deferred:** the multiprocess partition row (step 25) under the single-process directive; the Geant4 BGG/Barashenkov part of D9 (Amendment 10(d); the LA150 vs Tripathi-light part ran in C17: Tripathi-light not preferred, LA150 p+C −10.7 % at 70–110 MeV, no model adopted); V6 and the p-p elastic model D2 (V3-005C, Amendment 7(a)).
- 2026-10-09, hr5 at the same head (seed base 20451004, CUDA, archive `validation/results/transport/hr5-7aae5bb-single-process.json`):
  - **V8 statistical parity** (T12): every CPU/CUDA pair (cpu64:cuda32, cpu64:cuda64, cuda32:cuda64) passes for `idd`, `nuc_local`, `sec_p` and the four scalars (|z| ≤ 2.26 against 3.5). python:cpu64 `idd` (p 0.54) and `sec_p` (p 0.51) pass; `nuc_local` is **inconclusive** by the frozen T12 sparse-profile rule at the pre-registered 2.4e4 python histories (60 of 161 bins supported, 0.259 of the hull's mean deposit < 0.5; the supported bins give p = 0.75, the `nuclear_local` total z = −1.06). Recorded as inconclusive, not as a pass; the binding python ↔ warp evidence is the V8-LV trajectory parity above.
  - **V7 f32 vs f64 on CUDA:** pass (z −1.66, −2.16, −2.55 against 3.0; 1e6 histories in 10.7 s float32 and 15.0 s float64). All three float32 means are lower than float64 by 1.0–1.2 %, recorded as an observation.
- **Execution incidents (chain ba8751a0, superseded):** two native SIGSEGV of V7 shards (RUN-20261008T205127Z-ccc40171 after 691 s, RUN-20261008T213342Z-dbb8a88c after 318 s; the re-run RUN-20261008T210343Z-dccacb06 with the first crash's seed passed), traced in a sandbox reproduction with the fault handler to the per-block rebuild of the packed nuclear device (`build_event_model`); fixed by the per-process device cache (C37, 41e4bc36): bitwise-identical estimators, about 29 % less wall time per shard. The MCP client's 1800 s idle limit aborted the 1838 s `v7-rep-s0` call (RUN-20261008T201947Z-6a3036f4) while the host run completed; its digest came from the published step document, for that superseded chain only. All partials of the 7aae5bb1 chain came from the protected stdout.
- **Inherited by V3-005C:** D2 (p-p elastic) and p-nucleus elastic scattering, then V5 and V6 re-run on a fresh base; a pre-registered V7 re-test with a replicate-level criterion that tests the escaped-neutral hypothesis above; a larger python sample for the python:cpu64 `nuc_local` profile parity; the BGG part of D9; the multiprocess rows once the operator lifts the single-process directive.

- *Correction note 2026-10-10 (Codex REVIEW-7d749e6a):* the words "attributed to the absent hadronic elastic channel" (slice-B addendum on V7 replicates) and "Attribution" in the V5 bullet above are to be read as "hypothesis: the absent hadronic elastic channel (to be tested by V11/V5 in V3-005C; the exploratory X cases show elastic is not the only contributor)". The original text is left unchanged.

## Slice C (V3-005C, executed on the V3-005B branch) — ratified 2026-10-10 (revised the same day after Codex REVIEW-7d749e6a), before any slice-C result

Acceptance-plan Amendment 14 is the binding record for rows, tolerances, seeds and execution. This section records the physics decisions (orchestrator D-1 to D-6) and the alternatives that were rejected.

The operator resolved IR-20261009-032741-C590CF on 2026-10-10. Under that resolution:
- PR #82 stays unmerged, and V3-005C continues on this branch under task id V3-005B.
- All failed and inconclusive records are preserved.
- Follow-up validation is pre-registered before new results are generated, and no threshold is relaxed or seed repeated to obtain a pass.
- In the operator's words, "missing elastic scattering is a hypothesis to test, not an assumed complete explanation".
- The task merges only after the scientific and the implementation gates pass.

### D-1: model for hadronic elastic scattering (D2)

1. **p + nucleus (C, N, O, Ca and the surrogates of section 1), 1–250 MeV.**
   - The integrated cross section is σ_el(E) = max(σ_tot(n+A) − σ_inel(p+A), 0) from the Barashenkov tables, as used by Geant4 `G4BGGNucleonElasticXS` / `G4NucleonNuclearCrossSection::GetElasticCrossSection`. It uses linear interpolation in E, A^{2/3} interpolation for untabulated Z, and the BGG rule below 14 MeV.
   - The numbers are hand-transcribed from Geant4 11.4.2 `G4BarashenkovData.hh` and the formulas are re-implemented by hand. No code is copied. The Geant4 Software License notice and per-file sha256 are recorded.
   - **Lineage:** TOPAS (QGSP_BIC_HP with `g4h-elastic_HP`) uses the same σ_el, so V5 and V11 cannot test this normalisation.
   - The angular distribution is the nuclear-only black-disk form |2J₁(qR)/(qR)|², with q = 2p_CM sin(θ_CM/2). The radius comes from π(R + ƛ)² = σ_nonel(E) of the transport's own non-elastic table (LA150 MT5, Tripathi-shape extension above 150 MeV). The form is parameter-free, and no elastic data enter the shape.
   - It is tabulated as σ_el(E) plus an inverse-CDF in μ_CM per energy node and sampled over all angles.
   - Planning values: R_O = 2.58 / 2.60 / 2.64 fm and θ_rms(p+O) = 21 / 17 / 14° at 100 / 150 / 200 MeV. σ_el(p+O) = 1060 / 817 / 352 / 159 / 102 / 92 mb at 10 / 50 / 100 / 150 / 200 / 250 MeV.
2. **Coulomb scattering** stays entirely in the multiple-scattering model, and no Coulomb–nuclear interference is modelled for p+A. This is the same split as in the gating engine (hadronic elastic plus EM Coulomb).
3. **Recoil nucleus.** It is deposited locally via `nuclear_local` (species 64, class local, excluded from LET). ⟨T_r⟩ is about 0.55–0.8 MeV at 100–150 MeV for O, and its range in tissue is ≪ 10 µm.
4. **p-p (H-1).**
   - ≤ 150 MeV: the LA150 Hale R-matrix evaluation (LAW=5, LTP=1, LIDP=1), reconstructed with ENDF-102 eqs 6.9/6.10/6.14. Units are b/sr CM, so integrated cross sections carry 2π, and the identical-particle interference coefficient is −1/2.
   - The tables interpolate the ratio σ_e/σ_c in μ and between energy nodes, never P_NI. The transported density is NI = σ_e − σ_c over θ_CM ≥ 16.26° (|μ_CM| ≤ 0.96, the NJOY `umin` convention), where NI ≥ 0 at every checked node.
   - Above 150 MeV: σ(E) = σ_NI(150)·S(E)/S(150), with S the BGG/PDG systematics 1.0115·(23 + 50·√(ln(0.73/p)^7)) mb, and the shape held fixed in μ_CM.
   - Kinematics are exact relativistic two-body. The faster proton continues as the primary, and the slower one (⟨T⟩ ≈ T/4, 37.9 MeV at 150 MeV) is transported as a secondary proton.
   - The p-p channel is independent of TOPAS (BGG pp).
   - Hale's π·b₀ is 14–16 % above the PDG σ_el(pp) at 100–150 MeV, while the transported NI above 16° agrees with PDG to 1–7 %. The excess is recorded, not corrected, because it is a property of the evaluation.
5. **Sampling.** Σ_tot = Σ_nonel + Σ_el in the existing thinning, with the channel chosen by Σ_el/Σ_tot. `majorant_violation` stays fail-closed, the EM streams are identical with nuclear on and off, and the bounds of section 5 (32 particles, CUDA chunk ≤ 2^14) are unchanged. Elastic adds at most one secondary per event.

### Rejected alternatives

| Option | Reason |
|---|---|
| (i) LA150 LTP=12 nuclear-plus-interference for p+A | **Data finding:** above about 20 MeV the O-16 MT2 data (MF3 and the full MF6 P_NI tables, 24–150 MeV) are a numerical copy of C-12. MF3 agrees to ≤ 1e-5 relative. The backward floor σ_e/σ_c = 0.4152 at 100 / 120 / 150 MeV equals 1 − σ_c(Z=6)/σ_c(Z=8) exactly, i.e. an NI clipped against the carbon Coulomb amplitude and re-used with the oxygen one. The MF1 texts give no hint of copying. ENDF/B-VIII.1 changed only ⁴He in the proton sublibrary. O carries about 89 % of the nuclear-elastic mass in water. Beyond the copy: lin-lin P_NI goes negative between 5° and 10° (C-12, Ca-40), and CNI above 5° would double-count with the multiple-scattering model. Kept only as the report-only cross-check X-ENDF (C-12, N-14, Ca-40). |
| (ii′) shape from the ENDF/B-VIII.0 **neutron** O-16 / C-12 MF4/MT2 | The isospin-mirror error is 5–10 %. It needs a new acquisition and an MF1 copy check. It is the **pre-registered fallback**, entered only if V10-A (a) or (b) fails for p+O while V10-A (c) and P7 pass (Amendment 14 (e)9). Its V10-A results would be exploratory. |
| (iii) global optical model plus a partial-wave solver | The verification burden is high. Koning–Delaroche is not valid for A < 24. The Weppner 2009 light-nucleus fit contains Rolland 1966, which would remove the independence of V10-A. |
| TENDL-2023 p-O016 | Not reachable from the sandbox (proxy 502), so its content is unverified. |
| Gaussian form exp(−q²R²/4) | It is the small-q limit of the adopted form and has no diffraction minimum. Kept as a report-only sensitivity. |
| R from the Barashenkov σ_inel, or the CHIPS t-slope | Either would make the angular shape share TOPAS lineage too. |
| p-p cut at 5° or 10° CM | NI is negative below about 14° CM at ≤ 100 MeV (down to −11.2 mb/sr at 50 MeV), so the density would not be positive-definite. |
| Local deposition of the p-p recoil proton | Its ranges reach centimetres, and it carries about 1–1.5 % of E₀. |
| Correcting Hale's π·b₀ excess | It is a property of the evaluation, and the transported NI already agrees with PDG. |
| CNI for p+A | Not in the gating engine; Coulomb is owned by the multiple-scattering model; the ENDF CNI is unusable for O. |
| Recording block sums of all V7 bins (28.8 MB) | Only `escaped_neutral` failed. `sec_p` and `nuclear_local` passed both Amendment 13 gates. |

### D-2 to D-6 (summary; details in Amendment 14)

- **D-2:** V10-A is **gating** for any p+A elastic validation claim. Its tolerances are frozen now: ±30 % per angle at 10° / 15° / 20° CM, 20 % on the partial integral, 12 % on the C level. EXFOR is acquired from the IAEA raw master with the GitHub mirror as the recorded fallback, and no intervention is needed.
  - **Data role (revised after Codex REVIEW-7d749e6a).** The pre-1997 p+16O / p+12C angular data are **related-model evidence (shared/related lineage)**: Kelly 1989 / 1990, Seifert 1993, Glover 1985, Meyer 1981 / 1983 / 1988, Rolland 1966, Strauch & Titus 1956. Two lineage paths connect them to the model:
    - LA150's σ_nonel, which sets R and hence the whole angular scale qR, may have been adjusted to the same elastic measurements;
    - Barashenkov's compilation, which sets the level, may contain them.
  - V10-A stays gating as a **model-consistency test**. **No fully independent p+O elastic angular data exist.**
  - **Level.** The σ_el level is gated only by post-1997 data for C (Abfalterer 2001 n+C σ_tot minus Auce 2005 / Ingemarsson 1999 p+C σ_R). For O, only Ingemarsson 1999 (65.5 MeV) minus the pre-1997 Finlay 1993 n+O σ_tot exists, so it is report-only. **The p+O elastic level is therefore an evidence gap.**
- **D-3:** the (ii′) trigger above is adopted, and its results are exploratory.
- **D-4:** the V7 re-test records the 20 block sums per replicate for `escaped_neutral` only and judges it with the bootstrap-t interval under the two-gate structure of Amendment 13. `sec_p` and `nuclear_local` keep their Amendment 13 rule at the new base. The Amendment 13 failure stands whatever the re-test gives.
- **D-5:** seeds are lv5c 20481004 and hr5c 20491004; rehearsals use 20505000–20509999; later fresh bases are 20511004 (lv) and 20521004 (hr).
- **D-6:** all work happens on this branch.

### Additional limitations (slice C)

1. The p+A elastic **level** is shared with TOPAS. For O it is not independently tested (D-2).
2. The black-disk radius from σ_nonel is about 2.6 fm for O, about 20 % below typical strong-absorption radii. The angular distribution may therefore be too wide; V10-A (a) and V11 discriminate this. At 10–50 MeV, p+O elastic (σ_el ≈ 0.8–1.1 b, θ_rms 28–90°) is outside the diffraction regime (compound elastic). Its dose effect is small because the residual ranges are short, but it raises the event count to P(p-O) ≈ 0.21–0.24 per primary at 150–200 MeV.
3. No CNI. Coulomb remains in the multiple-scattering model only.
4. V10 is report-only. S(E) is PDG/BGG-derived, the PDG point at 160 MeV was seen during research, and no p-p elastic data independent of the Hale fit and of the PDG compilation were identified above 150 MeV (evidence gap).
5. V6 has shared beam lineage: p1 and p5–p7 come from the paper's fit to the same `tbl:dmlg` data. No independent entrance-beam measurement exists in the source, and none was obtained. Every r = 0 point gates at 2 %, and the r ≥ 1 cm halo gates only under the V6-sens condition (Amendment 14 (g)).

### Slice C work breakdown

C0 (Amendment 14 and this section) to C10, as listed in Amendment 14 (m). Every commit goes through the Codex review loop until it passes. Validation is recorded at the exact SHA. A Codex approval is not scientific validation.

### Revised 2026-10-10 (C2b/C2c)

These bullets record the slice-C revisions after Codex REVIEW-c856e576 and REVIEW-54fe3b7b (plan Amendment 15 and its 2026-10-10 row-table addendum). The text above is left unchanged; where it differs from these bullets, these bullets apply. No ionmc transport with elastic scattering had been run when they were written.

- **Model domains (supersedes "1–250 MeV" in D-1 1 and the p-p range implied in D-1 4).** The p+A elastic channel is defined per target on [e_min_shape(target), 250] MeV, with σ_el = 0 below: C 6.0, N 3.252, O 6.585, Al 2.588, Si 3.971, P 2.475, Ca 5.736 MeV. This is a declared model limitation (BGG σ_el below these energies is large and is not transported; Amendment 15 (b)). The p-p channel is defined on [12.532, 250] MeV (E_min,pp), with σ_pp = 0 below. The table JSON carries `elastic_domain` for all eight targets; the loader rejects a table with a non-zero cross section below a target's domain.
- **P6 failure (supersedes "NI ≥ 0 at every checked node" in D-1 4).** The frozen P6 item (4) negativity check FAILED for the LA150 H-1 evaluation below 12.532 MeV: 542 nodes and node midpoints have a negative transported nuclear-plus-interference density at the 16.26° CM cut, the last one −0.0219 mb/sr at 12.5 MeV (table JSON `p6.first_negative_nodes`, `negative_density_below_e_min_pp = true`). The failure is preserved and is not converted into a pass. P6-D is the pre-registered revised-domain row (Amendment 15 (a)4 and the row-table addendum); P6 items (1)-(3) and (5)-(7) are unaffected.
- **Omitted p-p correction below E_min,pp (total-variation diagnostic).** The S-wave unitarity "bound" quoted in Amendment 15 (a)2(iii) is not a bound on the full amplitude and is withdrawn. The record is the total variation of the signed nuclear-plus-interference correction to Rutherford scattering: M(E) = 2π ∫ |ρ_NI| dμ_CM over the half sphere above the cut, from the Hale reconstruction at the grid nodes in [1 MeV, E_min,pp], integrated with n_H = 6.69e22 cm⁻³ over the residual path in water (project range table) down to 1 MeV: 4.775e-3 per 150 MeV history (`pp_omitted_ni_correction_total_variation_per_history_150mev`) and 9.23e-3 MeV (`pp_omitted_ni_correction_weighted_total_variation_mev`, recoil-energy weighted). These are event-equivalent magnitudes, not an expected number of physical events and not energy transferred; the signed correction can cancel and the physical σ_pp below E_min,pp is not available from the evaluation (JSON `semantics`). The residual range below 1 MeV is reported separately, not extrapolated. The 15 MeV build ceiling on E_min,pp is pre-registered, not derived from this measure.
- **S(E) above 150 MeV (supersedes "BGG/PDG systematics" in D-1 4 and "S(E) is PDG/BGG-derived" in limitation 4).** The scaling S(E) is constructed from the Geant4 p-p BGG formula only (`G4HadronNucleonXsc` / `G4BGGNucleonElasticXS`, 11.4.2). The PDG rpp2022 p-p compilation is evaluation data (report-only V10, shared lineage with the BGG fit); it was registered as construction in error, never parsed by the builder, and is removed from the table identity. Role corrected through `assign_reference_data_role`, use id 770debc97dbd4d00b4cb2e6e44a6a8b0.
- **O-16 findings.** MF3/MT5 sigma is not a C-12 copy: the relative difference exceeds 1e-3 at a fraction 1.0 of the 30 shared MF3/MT5 nodes in 7-150 MeV (required > 0.5; median 0.22, maximum 0.28), CV of σ_O/σ_C 0.058 (required > 1e-3). MF6/MT5 yields, rule as frozen (every compared energy in the denominator, both-zero nodes not differing; JSON `mf6_rule_as_frozen`): FAILS, `passes: false`, because two products fall below 50 % owing to both-zero nodes, zap 3007 at 14/30 = 0.467 (16 both-zero nodes) and zap 5012 at 10/30 = 0.333 (20 both-zero nodes); preserved as a recorded finding. Revised rule (Amendment 15 addendum, (d)1 specification revision before any transport result; JSON `mf6_rule_revised`): the denominator is the informative energies (at least one material nonzero), a product gates with at least 10 informative energies, else report-only, tolerances unchanged; every product gates and differs at fraction 1.0 (informative counts: n 19, d 22, zap 3007 14, zap 5009 24, zap 5012 10, zap 7012 17, zap 7015 23, zap 9016 19, all other products 30; medians n 0.117, p 0.082, d 0.268, α 0.688, γ 0.307; the proton minimum is 1.2e-3). `o16_mt5_not_c12_copy` is set from the revised rule. The MT2 copy finding stands as a recorded finding (`o16_mt2_copy_recorded`; MF3 maximum relative difference 9.09e-5 over 34 nodes from 24 MeV, tolerance 1e-4); MT2 is never a construction input.
- **MCS has no hadronic term (Amendment 15 (d)2).** `test_mcs_scattering_power_has_no_hadronic_term` passes: the transitive in-package import closure of `ionmc.physics.em` and `ionmc.physics.scattering` contains no `ionmc.nuclear.*` module, the argument list of `scattering_power_dm` is the electromagnetic one, and the value matches the EM formula. Caveat: the closure contains the EM-only data modules `ionmc.data.cache`, `ionmc.data.nist_star` and `ionmc.data.registry` (through the unit constant `K_MEV_CM2_MOL` of `ionmc.physics.stopping`); `docs/architecture/transport.md` records this. Provenance of f_dM: decision 0039 and `docs/physics/em-transport.md` cite Gottschalk 2010 without stating its calibration target; that it is calibrated to Molière (Coulomb-only) theory is unverified here (source text not in the local cache) and stays an open finding to resolve before C4.
- **Diagnostic counters reserved.** `elastic_below_domain` and `pp_below_domain` (Amendment 15 (b)2) are reserved for the transport steps C4/C5; they never set `valid = False` and are part of the P5-ext/V8-LV counter parity.
- **Table of record at the C2b build:** id f932d6ae21dd24135a8f98f8dcb7c8357dbf292e3c7f6ce183e64fb194a8f584 (npz sha256 68c765f5…, 400 nodes, schema `ionmc-elastic-proton-table-2`). C2c build id 25361f8125c9ae03ed4e2b1808ab5ba90a85b6738a57c5cc7c83ca0173a3cfcc (builder -3). **C2d build (table of record):** id 045c31a2b5681b2d8b718c14645d562ecb787c5c81f0871b24bf092d4004d95c, builder `ionmc-elastic-proton-builder-4`, identical npz arrays (sha256 68c765f5…); only the JSON diagnostics changed (total-variation relabelling, MF6 as-frozen and revised records).
