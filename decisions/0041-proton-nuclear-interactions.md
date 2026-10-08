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
