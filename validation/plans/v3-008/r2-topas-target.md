# R2: TOPAS-agreement validation target (V3-008 acceptance plan revision 2; planner record, 2026-11-07)

Prepared 2026-11-07 by the planner at task head `d0829601` (revision 1 frozen at the reviewed freeze commit). This
record is the design behind plan revision 2 (`validation/plans/v3-008-acceptance.md`, section "Revision 2") and
decision 0043 revision 2. It is written **before any result**: at `d0829601` the V3-008 code has produced no ionmc
ion-transport run, no LV/HR qualification step and no scored ion reference-engine (TOPAS) run. Nothing in this record
was run to write it. The planning ranges below come from the project Bethe table (ionmc's own construction model,
`BetheStoppingSource()` default options, water, rounded to 0.1 mm); they size boxes and bins and are not reference
values. No NIST/ICRU/ENDF value is quoted.

## 0. Trigger: the operator's answer to IR-20261102-010435-B59B2D

- Directive, recorded verbatim by the orchestrator:
  `reference_requests/operator-directives/2026-11-07-topas-validation-target.md`
  (sha256 `ca079e9f0af990cdff9d8efde60d38ac55e2d0fc0fb9ef22355a81f4a2ab0b26`). The intervention was resolved on
  **2026-11-07 at 20:37 UTC** (resolution record under `reference_requests/interventions/IR-20261102-010435-B59B2D/`,
  handled by the orchestrator; not bound by this planner record).
- Content, as relevant to the plan:
  1. the acceptance reference for depth-dose distributions and lateral profiles is **agreement with TOPAS**;
  2. the final reference configuration is OpenTOPAS 4.3 / Geant4 11.4.2 with
     `sv:Ph/Default/Modules = 7 "g4em-standard_opt4" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4h-elastic_HP" "g4stopping"
     "g4ion-QMD" "g4radioactivedecay"` (here: the **FULL** configuration, "the specified TOPAS configuration");
  3. simpler TOPAS physics, including EM-only, may be used for intermediate validation (here: the **EM**
     configuration);
  4. cases, materials, energies, observables, statistical precision and acceptance criteria are chosen and justified
     autonomously;
  5. complementary checks (libdEdx I3 and similar) are retained;
  6. the proton-equivalent (f_dM) extrapolation is a model-design decision, assessed through the ion predictions
     against the reference;
  7. experimental validation is deferred. Achieved validation is described as agreement with the specified TOPAS
     configuration. Existing results and historical verdicts are preserved.

## 1. What revision 2 changes and keeps

| Item | Revision 2 |
|---|---|
| V1, V1-held, V2, V4 (TOPAS EM-only rows of revision 1) | **Superseded** by the TD/TL rows below. Never run; no seed consumed; their text stays in the plan as history. V1's tolerances (a)–(d) are carried into TD-EM unchanged |
| V3 (measured C/O range) | Stays recorded as an **EVIDENCE GAP under revision 1** (P6: Sihver 1998 FAILED). Superseded as a validation target by the TD rows (operator directive); the gap is not relabelled as closed. Experimental validation is deferred |
| Lateral evidence gap (P6 lateral search; intervention item (1)) | Recorded as history. Superseded as a validation target by the TL rows (agreement with the specified TOPAS configuration). Independent physical evidence of the charge dependence of lateral spreading remains absent; that is the deferred experimental validation, stated in the ledger wording |
| f_dM extrapolation (intervention item (2); I11-ext) | Model-design decision (unchanged model). I11-ext stays report-only. Adequacy is judged by **TL-EM (b) and TL-FULL (b)** (section 5.3) |
| V5-SS | **Kept unchanged** (gating, complementary, shared screened-Coulomb lineage; 2.94 h REF fits the budget, section 8). It is the only row that resolves the charge dependence of the scattering (b(z)) at the 3 % level; the TOPAS opt4 lateral rows cannot (section 5.3) |
| I0–I16, I11-map/rel/abs/ext, T1/T12/T13/R1/R2/C-res, I3/I3b/I3c | Unchanged |
| Seeds | Formula `20619000 + 100·part + case` and the guard test unchanged; case numbers re-assigned (section 7) |
| REF budget cap | 16 h unchanged; planning total ≈ 8.7 h (section 8) |

## 2. Reference configurations and lineage

| Config | Physics | Cuts | Role |
|---|---|---|---|
| **EM** (intermediate) | `sv:Ph/Default/Modules = 1 "g4em-standard_opt4"` | `CutForAllParticles` 0.05 mm, `CutForElectron` 1 km (P4 verbose smoke: e⁻ threshold 9.55 GeV in G4_WATER, no delta electrons), `EMRangeMax` 10 GeV | like-for-like test of ionmc's EM-only transport (ionmc has no ion nuclear model: I0 fails closed for `nuclear=True` with z ≥ 2) |
| **FULL** (final; the operator's configuration) | the seven modules above, verbatim | `CutForAllParticles` 0.05 mm (default of the proton full-physics cases), `EMRangeMax` 10 GeV | final acceptance reference, on observables that an EM-only code can be expected to reproduce (section 5); everything nuclear-dominated is report-only |

Both: G4_WATER (I = 78 eV), water box 120 × 120 × L mm, zero-size zero-divergence monoenergetic pencil on the axis
from z = −0.001 mm, `BeamEnergy` = total kinetic energy (P4 convention), He as `alpha` (P4: `GenericIon(2,4)` maps
to `alpha`), C/O as `GenericIon(6,12)` / `GenericIon(8,16)`, one thread, `Ts/MaxStepNumber` 1e7, vacuum world
400 × 400 × 1200 mm.

**Lineage (declared; these are SHARED-LINEAGE comparisons for stopping):**
- C/O stopping in TOPAS (P2, FC-2 FACT): ICRU 90 proton array at equal velocity × q_eff² + q_eff²(ΔL_LS +
  2·Barkas·(z−1)/z) above 2.0 MeV scaled. He: the ICRU 90 alpha array (≤ 1 GeV), BraggIon below 7.9452 MeV. Shared
  with ionmc: Bethe theory, I = 78 eV, the Barkas F table (transcribed from `G4EmCorrections`). For He the TOPAS
  stopping is the same data as I2's reference. R80 agreement is therefore agreement between two constructions of the
  same Bethe-type stopping at the same I-value, plus independent straggling, end-of-range and MCS-detour
  implementations. It is not independent stopping evidence; I3 (libdEdx ICRU73_NEW) stays the largely independent
  C/O stopping check above 2 MeV/u.
- MCS in TOPAS (P2 F1–F5, P4 FACT): Urban (`G4UrbanMscModel`) for alpha and GenericIon in opt4, also in FULL (opt4 is
  its EM module); WentzelVI + single Coulomb for protons. Urban is structurally different from ionmc's exact-Ω
  mapping of f_dM (independent implementation), but its θ0 carries |charge| only in the prefactor, with no charge or
  β in the logarithm, a truncated core plus tail, and a step-length dependent per-length width (P2 §2.1, §2.2, §6.2).
  Scratch estimate, not frozen (P2 §6.3): Urban/Highland-FE 0.835–0.89 per step 0.1–1 mm for C 290 MeV/u; with the
  Molière ion excess of +5 % (He) to +8 % (C/O) over z/(pv) scaling, TOPAS ion lateral widths are expected
  ≈ 10–30 % narrower than ionmc's. The TL bands are built on this (section 5.3).
- FULL adds fragmentation (QMD), hadron elastic, decay and radioactive decay. ionmc has none of these (V3-009).

## 3. Cases (water only)

Main set = the plan's standard points and ends, i.e. the exact I11/P11b grid. Reasons: (i) it spans the declared
domain (decision 0043 §5) at both ends; (ii) the f_dM domain status of every case is known pre-freeze from P11b
(SEEN theory, not TOPAS): the low ends are in the proton-qualified f_dM domain (He 50, C 80; O 100 mostly), the
standard points and upper ends are 100 % extrapolated (up to 0.8 GeV He, 2.5 GeV C/O), which is what the f_dM
assessment (TL (b)) needs; (iii) the shortest ranges give the tightest absolute R80 tolerance and test the
end-of-range handling (E_cut 2 MeV/u, P3 convention band); (iv) the upper ends maximise the fragmentation effect in
FULL and the extrapolation in TL. Held-out set (run once, section 6): interior energies with no I11 state.

| # | Ion | MeV/u | Particle | Total MeV | R_plan (mm) | Box L (mm) | IDD bin h_z (mm) | Lateral bin h_x (mm) | Set |
|---|---|---|---|---|---|---|---|---|---|
| 1 | He | 50 | alpha | 200 | 22.5 | 50 | 0.02 | 0.01 | main (low end, f_dM in domain) |
| 2 | He | 150 | alpha | 600 | 159.3 | 250 | 0.1 | 0.02 | main (standard point) |
| 3 | He | 250 | alpha | 1000 | 382.6 | 500 | 0.1 | 0.02 | main (upper end) |
| 4 | C | 80 | GenericIon(6,12) | 960 | 17.5 | 50 | 0.02 | 0.01 | main (low end, in domain) |
| 5 | C | 290 | GenericIon(6,12) | 3480 | 163.2 | 250 | 0.1 | 0.02 | main (standard point) |
| 6 | C | 430 | GenericIon(6,12) | 5160 | 309.3 | 450 | 0.1 | 0.02 | main (upper end) |
| 7 | O | 100 | GenericIon(8,16) | 1600 | 19.6 | 50 | 0.02 | 0.01 | main (low end, mostly in domain) |
| 8 | O | 350 | GenericIon(8,16) | 5600 | 166.6 | 250 | 0.1 | 0.02 | main (standard point) |
| 9 | O | 430 | GenericIon(8,16) | 6880 | 231.9 | 350 | 0.1 | 0.02 | main (upper end) |
| 10 | He | 100 | alpha | 400 | 78.0 | 150 | 0.05 | 0.02 | held-out |
| 11 | C | 150 | GenericIon(6,12) | 1800 | 53.1 | 100 | 0.05 | 0.02 | held-out |
| 12 | O | 200 | GenericIon(8,16) | 3200 | 65.5 | 100 | 0.05 | 0.02 | held-out |

Rules: L = 1.3 R_plan rounded up to 50 mm (distal edge and part of the FULL fragment tail inside the box);
h_z = 0.02 mm (R < 30 mm), 0.05 mm (30–120 mm), 0.1 mm (≥ 120 mm), so that σ_R/h_z ≥ 3 with the planning straggling
σ_R/R ≈ 0.60 % (He), 0.35 % (C), 0.30 % (O) (proton σ_R/R ≈ 1.1–1.2 % scaled by (m_p/M)^½; estimate, I4b form);
h_x = 0.01 mm (R < 30 mm) else 0.02 mm, so that the expected TOPAS core σ at z/R = 0.25 (≈ 0.03 mm for C 80 up to
≈ 0.5 mm for He 250; scaling estimate from the equal-velocity proton, σ_ion ≈ σ_p·(z/A)·1.05–1.08·(R_ion/R_p)^{3/2})
spans ≥ 3 bins.

**Deviation from the orchestrator's sketch, with reasons.** (i) IDD bins 0.02–0.1 mm, not 0.5 mm: the distal 80–20
width of C/O is ≈ 0.1 mm (low ends) to ≈ 1.5 mm (upper ends), so 0.5 mm bins put the distal edge inside one or two
bins, and linear interpolation of bin-averaged dose then compresses or garbles a true R80 shift (for an edge much
sharper than the bin, the apparent shift is between 0.2 and 0.8 of the true shift depending on the grid phase).
(ii) Projected 1-D lateral profiles at 0.01–0.02 mm, not a 0.5 mm 3-D grid: the expected ion core widths are
0.03–3.5 mm, most of them below 0.5 mm. (iii) Separate IDD and lateral runs (section 4.3).

## 4. Observables and estimators (identical algorithm for ionmc and TOPAS, on identical bin edges)

### 4.1 Depth dose
- E(z): laterally integrated energy deposit per primary per unit depth (MeV/mm) on the h_z grid over the full
  120 × 120 mm field. TOPAS: `IDD` dose × slab mass (120 × 120 mm × h_z at 1.0 g/cm³) / histories. ionmc: IDD
  scoring grid (1 × 1 × L/h_z) with the same edges.
- Peak: E rebinned to w_p = 5 h_z; D_peak and z_peak from a parabola through the maximal group and its two
  neighbours.
- R80 (R20): the first depth distal of z_peak at which E falls to 80 % (20 %) of D_peak, by linear interpolation on
  the h_z grid.
- Distal width W = R20 − R80.
- Entrance: D_ent = mean of E over [1, 3] mm (R ≥ 30 mm) or [0.5, 1.5] mm (R < 30 mm); PE = D_peak/D_ent.
- Plateau groups: E averaged over groups of g = 5 mm (R ≥ 120 mm), 2 mm (30–120 mm), 1 mm (< 30 mm) from the
  entrance-window start to 0.8·R80_TOPAS (last partial group dropped).

### 4.2 Lateral profiles
- Depths z_k = f_k·R80_TOPAS-EM of the same case, f_k ∈ {0.25, 0.5, 0.75, 0.9}; the 1 mm depth slab containing z_k.
  The same absolute slabs are used for ionmc and for both TOPAS configurations.
- Projected profiles: x profile (summed over the full y extent) and y profile, h_x bins, from the `LatX`/`LatY`
  scorers (EM) and `LatXPrimary`/`LatYPrimary` (FULL, generation-0 filter).
- **Core width σ_core**: binned maximum-likelihood fit of a centred Gaussian, integrated over each bin, to the
  profile truncated to |x| ≤ 2σ (window-truncated normal likelihood), iterated from the Sheppard-corrected second
  moment in ±20 mm until σ changes by < 1e-4 relative (≤ 50 iterations; non-convergence = invalid point, counted as
  FAIL for a gated point). σ_core = mean of the x and y fits. Robust to a broad pedestal and to the halo; no rms over
  the full profile.
- Report-only: the Sheppard-corrected second moment in ±20 mm; the halo fraction (dose outside ±3σ_core over the
  full field); in FULL also the all-particle core σ.
- Resolvability rule (fixed now, reference-only): a (case, depth) point is gated iff σ_core,TOPAS-EM ≥ 3 h_x; else
  report-only.

### 4.3 Why IDD and lateral runs are separate TOPAS runs
Whether TOPAS scoring bins act as step-limiting boundaries is OPEN (P2 §6.2). Urban's per-length width depends on the
step length (factor 0.756 between 0.1 and 1 mm per step, P2 §6.2). Fine depth bins (0.02–0.1 mm) in the lateral run
could therefore narrow the TOPAS lateral width by design. The lateral runs carry only the 1 mm-depth projected
scorers; the IDD runs carry only depth scorers (laterally integrated, insensitive to the lateral width). The
step-length smoke (section 6) records which regime applies; either way the bundle defines the reference.

## 5. Rows, criteria and tolerances (all TOST: PASS iff |x − r| + 1.645 σ_comb ≤ tol; inconclusive = FAIL)

σ_comb = √(σ_T² + σ_i²). σ_T: sample SD over the K = 4 TOPAS parts / √4, multiplied by t_{0.95,3}/1.645 = 1.430
(small-sample inflation). σ_i: ionmc batch SE (40 batches). Bands: section 5.3.

### 5.1 TD-EM (gating; intermediate; per case, every criterion must pass)

| | Observable | Tolerance | Defence |
|---|---|---|---|
| (a) | R80 | max(min(0.35 % R80_T, 1.0 mm), 0.05 mm) | V1's value kept. Budget (relative to range): stopping-construction difference, C/O ≤ 0.18 % (Bethe-78 vs integrated ICRU 90 proton range, decision 0038, SEEN) + ion terms ≤ 0.15 % (effective charge fully stripped above ≈ 10 MeV/u; LS vs Bloch + Mott ≤ 0.1 %, P3; Barkas construction and the possible z = 1 double count of ΔL_LS, estimate; the P3 convention band ≤ 1 µm); He ≤ 0.10 % (ionmc He vs the ICRU 90 alpha array CSDA range, decision 0038, SEEN); straggling shape relative to the mean range ≤ 0.02 %; MCS detour (Urban narrower) ≤ 0.01 %; same-grid interpolation ≤ 0.01 %. Linear ≤ 0.36 % (C/O) / 0.14 % (He); quadrature ≈ 0.24 %. The I-value ±0.3 % and the 1–2 % table uncertainties (I3, TO VERIFY) do not enter: both codes use 78 eV and Bethe-type tables; they are uncertainties against nature, not between the codes. The 1.0 mm cap binds only at He 250 and C 430 (≈ 0.26 % / 0.32 % of R80), inside the budget. The 0.05 mm floor is decision 0043's action threshold; it binds nowhere (smallest tolerance ≈ 0.061 mm, C 80) |
| (b) | W = R20 − R80 | 10 % of W_T | V1 (b). Gated only where W_T ≥ 1.0 mm (reference-only rule): below it, ionmc's local deposition below E_cut (residual range ≤ 0.075 mm, decision 0043 §4) can move R20 by up to half that, ≥ 3.75 % of W, and the grid resolves W only to a few bins. Else report-only |
| (c) | z_peak | 0.5 mm | V1 (c) |
| (d) | PE ratio | 3 % | V1 (d): straggling-theory spread is a few % on the peak height; E_cut and the peak-group width are identical in both analyses |
| (e) | plateau groups, E_i/E_T per group | 2 % per group, every group | New (IDD shape). Local stopping difference ≤ 0.35 % (0.2 % table + 0.15 % ion terms); the R80 tolerance projects to ≤ 0.76 % at 0.8 R (d ln E/dz ≈ 0.435/(R − z) for S ∝ E^−0.77); both deposit deltas locally (1 km cut); EM-only fluence = 1 per primary in both. Linear ≈ 1.1–1.3 %; statistics negligible (≲ 0.05 %) |

Report-only: the ratio curve E_i/E_T over [0, R80_T + 5σ_R]; R90; the C/O 430 equal-velocity R80 ratio.

### 5.2 TD-FULL (final; gating on (a) and (c) only)

EM-only ionmc against the FULL all-particle IDD.

| | Observable | Tolerance | Status |
|---|---|---|---|
| (a) | R80 of the **total** (all-particle) IDD | TD-EM (a) tolerance + A_frag, A_frag = 0.15 % R80_T (He), 0.10 % R80_T (C, O) | **gating** |
| (c) | z_peak (total IDD) | 0.5 mm | **gating** |
| — | PE, plateau ratio, W, the tail dose beyond R80_T + 3σ_R, the absolute entrance dose ratio | — | report-only (fluence attenuation by nuclear reactions changes them by up to tens of %; no correction model is declared) |
| — | R80 of `IDDPrimary`, primary fluence vs depth, Δ_frag = R80_T,FULL − R80_T,EM (TOPAS-internal) | — | report-only (attribution) |

**Why the total IDD, and A_frag.** The operator's target is the depth-dose distribution, i.e. the total dose; it is
also the clinically used curve. Smooth fluence attenuation exp(−z/λ) shifts a Gaussian-convolved edge by ≈ −σ_R²/λ
(≪ 0.01 mm). A fragment dose pedestal b near the peak, relative to the primary peak P, moves the 80 % point to the
primary's (80 − 20 b/P) % level if b is the same on both sides of the edge: +0.07 σ_R (b/P = 0.1) to +0.20 σ_R
(b/P = 0.3), deeper; if the pedestal ends before the edge the shift is shallower by a similar amount. A_frag =
0.25 σ_R covers |b/P| up to ≈ 0.35 in either direction, i.e. 0.15 % R80 (He) and 0.09 → 0.10 % (C), 0.075 → 0.10 %
(O). Disclosed SEEN support (V3-003 proton archive, `validation/results/reference/README.md`): the 150 MeV proton
R80 moves by −0.2 mm (≈ −0.11 σ_R) from EM-only to full physics in TOPAS. The `IDDPrimary` R80 is report-only because
the generation filter also removes the delta-electron dose (0.05 mm cut), so it is not like-for-like either.

### 5.3 TL-EM and TL-FULL (lateral)

TL-EM uses the EM all-particle profiles (no deltas with the 1 km cut, so all-particle = primary-ion profile).
TL-FULL uses the FULL **generation-0 (primary) profiles**: with the 0.05 mm cut, delta electrons (T_max up to
≈ 1.2 MeV at 430 MeV/u) and fragments broaden the all-particle profile around a 0.03–3.5 mm ion core, whereas the
generation-0 profile is the surviving primaries' track distribution weighted by restricted stopping, i.e. the
quantity ionmc transports. The all-particle core σ and the halo are report-only in FULL.

| | Statistic | Tolerance | Status |
|---|---|---|---|
| (a) | r = σ_core,ionmc/σ_core,TOPAS per gated (case, depth) | band: (r − 1) − 1.645σ ≥ −5 % **and** (r − 1) + 1.645σ ≤ +40 % | **gating**, every gated point, TL-EM and TL-FULL |
| (b) | **f_dM-extrapolation trend** per ion: D_ext = ⟨r⟩_upper end / ⟨r⟩_low end − 1, ⟨r⟩ = mean of r over the gated depths among f ∈ {0.75, 0.9} common to both cases | \|D_ext\| + 1.645σ ≤ 12 % | **gating**, TL-EM and TL-FULL; NOT EVALUABLE (recorded, not a pass) for an ion whose low-end case has no gated point at 0.75/0.9 |
| — | D_std = ⟨r⟩_standard/⟨r⟩_low − 1; r at the non-gated points; ±20 mm second moments; halo fraction; FULL all-particle core; σ_FULL,primary/σ_EM (TOPAS-internal); C/O 430 equal-velocity σ ratio ionmc vs TOPAS | — | report-only |

**Band (a), defended.** Lower edge −5 % = ionmc's declared ±5 % model uncertainty (decision 0043 §3), with no
lineage allowance on that side: every identified structural difference of Urban for ions makes TOPAS narrower.
Upper edge: 1.08/0.835 (largest Molière ion excess over z/(pv) scaling, smallest Urban/Highland-FE factor of P2
§6.3) × 1.05 (ionmc) × 1.04 (Urban lacks the β term, ≈ 0.038·Δln(1/β²) relative to its ≈ 0.9 core factor, larger at
low β) ≈ 1.41 → +40 %. What the band catches: MCS off or not applied to ions; a wrong power of z (z² instead of z in
θ: × 2, 6, 8); momentum per nucleon instead of total (× A); any ionmc width more than 5 % below TOPAS. What it does
not resolve: the z- and β-dependence of the Molière logarithm (6–8 % between He and C/O at equal velocity, V5-SS
discrimination) and the ≤ 1.5 % f_dM extrapolation component. **The achieved lateral claim is therefore "agreement
with the specified TOPAS configuration within the pre-registered band", not a precision test of the charge
dependence.** V5-SS remains the 3 %-level check of b(z) (shared screened-Coulomb lineage).

**Criterion (b): where the f_dM extrapolation is judged (answer to directive item 6).** The low-end cases map to
in-domain proton states (P11b), the upper ends to fully extrapolated ones (0.8 GeV He; 2.5 GeV C/O). Within one
ion, D_ext compares ionmc/TOPAS at extrapolated states with ionmc/TOPAS at in-domain states, so f_dM's absolute
water calibration and the energy-independent part of the Urban deficit cancel. Budget: f_dM extrapolation ≤ 1.5 %
(declared, I11-ext); Born-screening model trend with α′ = zZα/β 1–2 %; Urban's missing β term changes between the
ends by 0.038·ln(β²_up/β²_low) ≈ 5.1 % (He 50 → 250), 4.8 % (C 80 → 430), 4.0 % (O 100 → 430) in θ0 relative to its
≈ 0.9 core factor; Urban's per-step log with TOPAS steps scaling with the range ≈ 2 % (estimate). Linear ≈ 10.5 %;
with statistics (σ(D_ext) ≈ 1 %) → 12 %. Expected sign negative (≈ −4 to −7 %: the Urban deficit shrinks with
energy). A PASS shows that the extrapolation produces no energy-dependent error larger than ≈ 5–10 % relative to the
in-domain cases, at the resolution this reference allows; it cannot confirm the declared ≤ 1.5 % component, which
stays declared, theory-bounded (I11-ext) and carried in the ±5 %. A FAIL is a finding against the extrapolation (or
the declared Urban-trend estimate), recorded and attributed, never re-graded.

**Rejected: the ion/proton ratio at equal velocity against TOPAS.** TOPAS opt4 gives protons WentzelVI + single
Coulomb and ions Urban, so the TOPAS ratio σ_ion/σ_p measures the model switch, not the charge dependence; it
carries the full Urban deficit and cancels only f_dM's absolute calibration, which the band already covers. The
equal-velocity protons of the upper ends (≈ 434 MeV) would also need ≈ 1 m boxes. Recommendation (not adopted, needs
a reviewed amendment before any ion TOPAS lateral output is read): an additional intermediate configuration in which
ions use a Wentzel-type msc (e.g. `G4EmStandardPhysicsWVI`, alias and ion assignment TO VERIFY by source read and a
scorer-free smoke) would give an ion lateral reference with the z-dependent Born screening and could tighten (a) and
(b).

### 5.4 Held-out set (TD-H, TL-H): cases 10–12

Same criteria: TD-EM (a)–(e), TD-FULL (a), (c), TL-EM (a), TL-FULL (a) (no D_ext). Run ONCE, only after the main TD/TL
rows are graded and after any model amendment. It plays V1-held's role: after a model change the main rows are
consumed and the held-out set becomes the sole gating TOPAS-agreement evidence, at fresh ionmc seeds, unchanged
tolerances.

### 5.5 Report-only TOPAS-internal diagnostics (TF-DIAG)

Δ_frag per case; primary-ion survival (fluence) at R80; the FULL/EM ratio of the plateau and of PE; σ_FULL,primary /
σ_EM per depth; halo fractions. They involve no ionmc output and serve attribution only. **Attribution rule
(pre-registered):** if TD-FULL (a) fails while TD-EM (a) passes and |Δ_frag| > A_frag, the failure is attributed to
nuclear physics absent from ionmc (V3-009 scope); it stays FAIL, is never re-graded, and goes to the orchestrator
(plan STOP 5).

## 6. Statistics, order, blindness and run validity

**N per case** (TOPAS: K = 4 parts, separate host calls, each with its own seed; ionmc: one run per case).

| Config | TOPAS per part × K | Total | Driving criterion |
|---|---|---|---|
| IDD-EM | 2.5e4 × 4 | 1e5 | (d) PE 3 %: per-primary relative SD of the peak-group deposit ≈ 1, so σ_T(PE) ≈ 1.43/√1e5 ≈ 0.45 %; 1.645 σ_comb ≈ 0.75 % leaves ≈ 2.2 % for systematics. R80: SE ≈ 1.4 σ_R/√N (edge-quantile argument) ≤ 0.01 mm everywhere, far below 0.05 mm (R80 SE ≤ 0.05 mm needs only N ≥ (1.4 σ_R/0.05 mm)² ≈ 4.2e3 at He 250) |
| LAT-EM | 2.5e4 × 4 | 1e5 | core σ SE ≈ 0.25 % (x, y pooled); cheap |
| IDD-FULL | 5e3 × 4 | 2e4 | R80/z_peak only gated; survivors to the peak ≳ 50 %; R80 SE ≤ 0.02 mm |
| LAT-FULL | 5e3 × 4 | 2e4 | generation-0 core σ SE ≈ 0.6–0.8 % (survivors) |

ionmc: warp-cpu f64, 1e6 histories per case, 40 batches, all EM on, default options (LS off, E_cut 2 MeV/u), one
run with three scoring grids (IDD 1 × 1 × L/h_z; LatX and LatY projected, h_x, 1 mm depth, lateral extent at least
±20 mm), same bin edges as TOPAS. Seeds section 7.

**Order (pre-registered; revision 2 of this block, R2-TOPAS-3, 2026-11-07: blinding by ORDERING, not by output
suppression).** The controlled reference runner inventories and permanently archives the whole work directory and
has no deletion mechanism, so every TOPAS output, smoke outputs included, is readable once the run has happened. The
blinding control is therefore the order of events, with the ionmc sample seal as a hard precondition of ANY TOPAS
run of the revision-2 set:
0. **Sample seal (hard precondition of every TOPAS run of this set, smokes included).** Revision 2 frozen by a
   reviewed commit (this record, the bundles and the P9 additions bound); implementation, CI/LV/HR rows, T1-ion and
   T12-ion (their failures stop everything before any V/T row); then the ionmc TD/TL samples of EVERY main and
   held-out case produced at the frozen implementation SHA, their sealed records and sha256 committed in a reviewed
   commit, the **sample seal commit**. Its SHA and its review id are recorded in the plan (`v3-008-acceptance.md`) by a
   dated seal note before any TOPAS run; the seal commit is reviewed and pushed before any run. No TOPAS run of the revision-2 set (smoke or evidence) may start before that
   commit exists.
1. **Timing and step-length smokes** (rehearsal seeds 20600411–20600423; representative, configuration-specific
   scorers: each smoke bundle carries the SAME scorer block as its evidence configuration, so the scorer boundaries
   and bins that can limit TOPAS steps are exercised and the measured cost is the real cost): twelve timing smokes
   (He 250, C 430, O 430 × IDD-EM, LAT-EM, IDD-FULL, LAT-FULL; 1e3 primaries EM, 200 FULL) and one step-length smoke
   (C 290 IDD-EM configuration, 3 histories, `TrackingVerbosity` 1). Recorded: wall time, module acceptance (in
   particular `g4radioactivedecay` and `g4ion-QMD` in the seven-module list), process list (UrbanMsc for
   alpha/GenericIon), killed tracks, whether scorer bins limit the primary step. The runner archives the smoke
   outputs; they may be read for sizing only after the seal; they disclose nothing that can influence the sealed
   samples (sealed before the first run).
2. **Sizing note** (dated, after the seal, before any evidence run): measured s/primary t_s. The part plan holds if
   t_s × N_part ≤ 1200 s at the upper end of each ion (the most expensive case of each ion and config). Otherwise,
   in this order: (i) reduce N_part, keeping the totals ≥ IDD-EM 6e4, LAT-EM 2e4, IDD-FULL 1e4, LAT-FULL 1e4 (the
   statistical minima of the table above); (ii) IDD-FULL may add parts 16–17 and LAT-FULL parts 18–19 of the seed
   formula (K = 6 for that config; new bundles committed by the dated note); (iii) otherwise STOP (amendment). If the
   projected REF total exceeds the 16 h cap after (i), V5-SS takes its own pre-registered 4e4 fallback, then STOP
   (orchestrator decision on the budget). No case is dropped silently.
3. **Evidence runs** of the TOPAS parts, then TD-EM, TL-EM graded; then TD-FULL, TL-FULL; then V5-SS (independent of
   order); then the held-out set (once). The ionmc samples are never recomputed after any TOPAS output of the same
   row has been seen.

**Run validity (each TOPAS part; a violation invalidates the part, which is re-run at a fresh re-run base declared
by amendment, the invalid part preserved):** zero killed or stuck tracks in stderr; `NumberOfHistoriesInRun`
completed; the EM printout shows UrbanMsc for alpha/GenericIon and the module list as bundled; the replicate gate of
`compare_batches.py` (identical configuration except the seed) holds within each group.

**Run validity against the seal (R2-TOPAS-3; rule revised in R2-TOPAS-4 to git ancestry and controlled records).** A revision-2 TOPAS run (smoke or evidence) is valid only if BOTH of the following hold. (i) **Ancestry:** the `code_sha` of its run record (the clean committed SHA that the controlled reference runner executed; every run record carries it) has the sample seal commit as an ancestor, i.e. `git merge-base --is-ancestor <seal> <code_sha>` succeeds (the seal commit itself counts). (ii) **Review before start:** the seal commit has passed its independent review (its review id is recorded in the plan's dated seal note), and that review's `finished_at` precedes the run's controlled `reference_started` event (the runner's own start record). Both comparisons use records of the controlled services, never Git timestamps (a Git timestamp does not prove when a commit existed), and the run's `finished_at` is not used (a run begun before sealing but completed afterwards must not pass). The sample seal commit must be reviewed and pushed before any run; the runner executes a clean committed SHA; pushing the seal commit before any run is required by the gate sequence (push → CI → run), not enforced by the runner. The evaluator, `validation/scripts/ion_em/r2_evaluate.py` (to be implemented with the TD/TL rows), enforces exactly these two checks and fails closed when either record is missing (seal SHA, seal review id or its `finished_at`, run `code_sha`, `reference_started`). A run failing either check invalidates the TD/TL rows of its case: FAIL of the blinding control, never re-graded.

## 7. Seeds and bundles

TOPAS seed = `20619000 + 100·part + case` (formula, ranges and guard test unchanged). Case numbers (re-assigned;
none of the old numbers was ever run, so no seed was consumed): **1–12** = the cases of section 3 (previously V1/V2
matched 1–9 and V1-held 10–12). Part blocks per configuration: **IDD-EM parts 0–3, LAT-EM 4–7, IDD-FULL 8–11,
LAT-FULL 12–15**; parts 16–19 reserved for the FULL fallback of section 6. V5-SS keeps cases 22, 25–31 (parts 0–1).
Retired unused and kept retired: 23, 24, 32, 33. Unassigned (any use needs an amendment): 13–21, 34–40.
Examples: C 290 IDD-EM part 0 = 20619005; O 430 LAT-FULL part 15 = 20620509; highest revision-2 seed 20620512.
The 192 evidence seeds are distinct from each other and from the V5-SS seeds.

Bundle names (`validation/reference_cases/topas/`): `<ion>-water-<E>mevu-r2-<config>-part<p>` with ion ∈ {he,
carbon, oxygen} and config ∈ {idd-em, lat-em, idd-full, lat-full}: 12 cases × 4 configs × 4 parts = **192**.
Smokes: `<ion>-water-<E>mevu-r2-<config>-timing-smoke` (12; seeds 20600411–20600422 in the order He 250, C 430, O 430
× idd-em, lat-em, idd-full, lat-full) and `carbon-water-290mevu-r2-idd-em-steps-verbose-smoke` (20600423). Total
**205** bundles; every `case.json` passes `infrastructure/experiment_v3/reference.py::validate_case`, the seed in
`case.json` equals the native seed line, and `histories` equals `NumberOfHistoriesInRun` (checked 2026-11-07).
Hashes: appendix A.

**Dated note (2026-11-09, R2-SIZE-1 regeneration; orchestrator).** After the sizing note the bundle set is **209** bundles / 418 files: the 48 LAT-FULL evidence bundles carry `CutForElectron = 1 km` (all ions), the 32 carbon/oxygen IDD-FULL bundles carry the rule-(i) history counts (C 3333, O 2529 per part), and four execution-only smokes were added (`he-water-250mevu-r2-lat-full-ecut1km-timing-smoke` 20600424, `carbon-water-430mevu-r2-lat-full-ecut1km-timing-smoke` 20600425, `oxygen-water-430mevu-r2-lat-full-ecut1km-timing-smoke` 20600426, `carbon-water-290mevu-r2-lat-em-steps-verbose-smoke` 20600427); the three 0.05 mm LAT-FULL timing smokes already run stay byte-identical as the record of those runs. The 96 EM bundles and the 16 He IDD-FULL bundles are byte-identical to the frozen set. Bundle-set digest `e51d3b3d1bf00e51daadf889f6541bf8bd6ba8aea6453eee132110e07ff0c4d6` (appendix A, 2026-11-09 listing). The LAT-FULL `NumberOfHistoriesInRun` (5000) is re-sized by a further dated note after the ecut1km smokes. (2026-11-09, ListProcesses; Codex REVIEW-6b4731a7) The 192 evidence bundles and the four smokes 20600424–20600427 now request the per-particle process list (`b:Ph/ListProcesses = "True"`, `input.txt` only), so that the executed module set, `g4decay` included, is verifiable from the log; the thirteen smokes already run are unchanged; bundle-set digest `6d1ae3a06211831cfb7e8924c65d6d8fe31fae7781443ff84c077204045ea790` (appendix A, ListProcesses listing).

ionmc seeds: lv8 row block `20581004 + 100·r_index + k`, r_index **30–41** = cases 1–12 (declared now; `steps_v8.py`
must adopt them), k = 0 for the single run, k = 1, 2, … for chunks if a run is split to stay ≤ 1500 s per host call.
The held-out set after a model amendment uses the re-run base 20621004 with the same r_index.

## 8. Cost (planning values; replaced by the timing smokes)

Per-primary cost, assuming scorer bins limit steps (upper estimate): EM IDD ≈ 0.15 ms + (R/h_z) × 3 µs; EM lateral
≈ (R/1 mm + 2 σ_end/h_x) × 4 µs; FULL = EM × 10 (C, O) or × 5 (He) (proton full/EM lateral ≈ 1.1–2 in the V3-003
archive; carbon full ≈ 38 ms/history in V3-010, configuration not recorded here; fragments cross every bin).

| Item | Estimate | Largest part |
|---|---|---|
| IDD-EM, main + held-out (1e5 each) | 5250 + 1220 s ≈ 1.8 h | He 250: 2.5e4 × 11.6 ms ≈ 290 s |
| LAT-EM (1e5 each) | 1050 + 135 s ≈ 0.33 h | He 250 ≈ 73 s |
| IDD-FULL (2e4 each) | 8500 + 1960 s ≈ 2.9 h | C 430: 5e3 × 94 ms ≈ 470 s (1175 s at × 25) |
| LAT-FULL (2e4 each) | 1640 + 200 s ≈ 0.5 h | C 430 ≈ 95 s |
| 13 smokes | ≈ 0.25 h | — |
| V5-SS (unchanged) | 2.94 h | ≤ 1200 s per part |
| **REF total** | **≈ 8.7 h** (≈ 14.3 h if FULL costs × 25 EM) | cap 16 h |
| ionmc TD/TL samples (12 × 1e6, warp-cpu) | ≈ 0.5–1.5 h LV (fine tallies unmeasured) | ≤ 1500 s per call |

Host calls: 192 evidence parts + 13 smokes + 16 V5-SS parts = 221 REF calls.

**Dated note (2026-11-08, sizing/mapping amendment).** The ionmc timing smoke RUN-20261107T224820Z-92f9f0be (1e4 histories per case, root-solve mapping) measured per-case ms/history: case 1 7.66, 2 12.26, 3 20.43, 4 8.01, 5 13.39, 6 20.69, 7 8.40, 8 14.44, 9 17.87, 10 9.30, 11 9.41, 12 9.96; projected 42.2 h for the twelve 1e6-history samples with the root solve (the planning value above, 0.5–1.5 h, is superseded). The exact-Omega root solve was ≈ 90 % of that cost (39.5 µs per step). Amendment: mapping evaluation = tabulation (decision 0043 section 3, "Evaluation", the allowed alternative), chosen before any revision-2 ionmc sample or TOPAS evidence run; reason: measured cost. Table: `pv_eq/pv` per (material, ion) on a uniform 128 x 128 grid in (ln pv, ln(p1v1/pv)) (a shear of the decision's (ln pv, ln p1v1) plane), nodes computed by the unchanged root solve `pv_equivalent` at validation (only the nodes of cells meeting the band of the run's initial states pv(E0 - 6 sigma_E) ... pv(E0 + 6 sigma_E); every other node invalid, a state there fails closed as `mcs_mapping_failure`), bilinear interpolation (`exact_omega.pv_equivalent_table`, shared by the python twin and the Warp kernel); option `PhysicsOptions.mcs_mapping_eval` = "table" (default for z >= 2) or "root" (the per-step root solve, retained); the summary records the evaluation, table sizes, axes, domain and the sha256 of the table rows. Bound test `tests/ionmc/test_ion_mcs_table.py::test_table_variance_error_bound` (He, C, O in water and skeletal muscle; >= 2000 random states plus edge states per ion and material, upper and lower domain end, thin-slab edge included; bound |var_table/var_root - 1| <= 1e-3): measured maximum 2.05e-5 (C in skeletal muscle; He 4.4e-6 / 5.5e-6, C 1.83e-5 / 2.05e-5, O 1.84e-5 / 2.04e-5 in water / muscle); transport-level counterpart `test_v3_008_rows_ci.py::test_i11_map_table_against_root_at_transport_states` max 6.8e-6. Measured MCS cost per step (warp-cpu float64, single thread, MCS-on minus MCS-off, 2000 histories, CI configuration): He 100 MeV/u 0.16 us (root 39.85 us), C 150 MeV/u 0.17 us (root 39.57 us). The model (decision 0043 section 3) is unchanged; the root solve remains available as an option. The revision-2 sample cost with the table is not yet measured at the case configurations; it is re-measured by a timing smoke before the samples. (2026-11-08, Codex REVIEW-d3a17ebf) The bound is in addition certified deterministically at build for every produced table (`exact_omega.certify_pv_eq_table`: the root solve at 5 interior points per solved cell, the centre and +-1/4 of the cell in both axes; per-step variance of the first and later steps; a table above TABLE_VARIANCE_BOUND = 1e-3 fails closed; the maximum, cells and points are recorded in the summary), measured maxima at the upper domain end (128 x 128, 628 cells, 3140 points per table) He / C / O: water 4.11e-6 / 1.97e-5 / 1.98e-5, skeletal muscle 4.02e-6 / 2.20e-5 / 2.20e-5, over all eleven built-in materials and a user-defined material at most 4.44e-6 / 2.87e-5 / 3.21e-5 (`test_ion_mcs_table.py::test_certification_every_builtin_and_user_material`). (2026-11-08, Codex REVIEW-9529fd35) The five-sample certification is replaced by a domain-wide analytic bound per solved cell (`exact_omega.certify_pv_eq_table`): the tensor-product bilinear error bound hx^2/8 sup|q_xx| + hy^2/8 sup|q_yy| of q = pv_eq/pv, second derivatives from second differences of the root solve on a 2x-refined grid (3 x 3 stencil per cell, maximum over the cell and its solved neighbours, which covers cells where a second derivative changes sign), times safety factor 2 (q is smooth in ln pv and ln r, P11; the factor covers the variation of the second derivative within a cell; empirically, over the 70 tables below the sampled cross-check passes everywhere with factor 1 and fails in every table with 0.75, so factor 2 leaves a per-cell margin of at least 2), converted exactly into the relative per-step variance error (only pv_eq differs between the paths and f_dM is linear in log10 pv_eq at fixed memory term l1: var_table/var_root - 1 = (c2 + c3 l1) dl2 / f, maximized at the corners of each cell's (l1, l2) rectangle) for later steps and for the first step over every physical end point pv_end in (0, min(pv_mid, (1 - 1e-4) p1v1)] (the first-step ratio is endpoint-dependent, c3 != 0); 5 x 5 interior samples per cell must stay below each cell's analytic bound, else the build fails closed, as it does when the bound exceeds 1e-3; scope: states of a particle that has lost less than 1e-4 of its initial pv (TABLE_MIN_LOSS) are outside the relative contract, because the memory term drives the f_dM fit to its zero there (e.g. He at 50 MeV/u in aluminium, copper, lead, compact bone and the user material, where the fit is clamped to 0 for both paths); measured analytic bound / sampled maximum, upper and lower domain end: water He 1.10e-5 / 4.1e-6 and 1.36e-5 / 3.2e-6, C 4.38e-5 / 1.97e-5 and 1.47e-5 / 7.1e-6, O 4.37e-5 / 1.98e-5 and 1.94e-5 / 9.3e-6; skeletal muscle He 1.09e-5 / 4.0e-6 and 1.32e-5 / 3.1e-6, C 4.90e-5 / 2.20e-5 and 1.46e-5 / 7.0e-6, O 4.91e-5 / 2.20e-5 and 1.94e-5 / 9.3e-6; largest analytic bound over the eleven built-in materials and a user material (70 tables) He 1.36e-5 (water, lower), C 6.37e-5, O 7.20e-5 (beryllium, upper), at most 7.2 % of the 1e-3 threshold; oxygen at 100 MeV/u in copper and lead cannot be tabulated (a node has no root; fails closed before certification, root option remains); build plus certification 0.5-1.6 s per (material, ion) table. (2026-11-08, Codex REVIEW-ee53867c) The table is demoted to an accelerator and the two preceding certification sentences are superseded: mapping evaluation = exact root solve on every step in both modes (the same Illinois iteration, tolerance |residual| <= 1e-13 or step/bracket <= 1e-14 in ln pv_eq, float64, and status codes as the root path); with mcs_mapping_eval = "table" the bilinear pv_eq/pv value only supplies the starting bracket (one end at the table value, the other max(10 x the table's build-time interpolation estimate, 1e-5) away on the side of the sign change, widened once x10, else the full proton-path bracket, counted in the diagnostic counter mcs_table_bracket_miss, which never invalidates a result); no interpolation bound is claimed (the build-time estimate is a performance diagnostic) and the TABLE_MIN_LOSS scope is withdrawn; the decision-0043 test of the table evaluation is the per-step root guarantee (`test_ion_mcs_table.py::test_accelerated_root_equals_root`: He, C, O at both domain ends in the eleven built-in materials and a user material, >= 2000 random plus edge states each, later and first steps with both end-point extremes: equal status, pv_eq within 1e-8 of the root path (observed <= 1.43e-9), variance within 1e-8 times the f_dM sensitivity (observed <= 5.3e-10), 0 bracket misses in 427 704 evaluations; `test_v3_008_rows_ci.py::test_i11_map_table_against_root_at_transport_states`: residual <= 1e-10, observed <= 9.4e-14); oxygen at 100 MeV/u in copper and lead still cannot be tabulated (build fails closed, root option); measured MCS cost (warp-cpu float64, single thread, MCS-on minus MCS-off, 2000 histories, CI configuration): He 100 MeV/u 21.00 us/step (root 40.04), C 150 MeV/u 21.65 us/step (root 39.63), 4.58 / 4.84 ms per history against 8.66 / 8.79; Illinois iterations per step on random table-domain states: accelerated mean 3.8-5.7 (max 26), root 7.8-8.7 (max 33); build 0.11-0.16 s (water) and 0.21-0.31 s (skeletal muscle) per (material, ion) table.

**Dated sizing note (2026-11-09, plan section 6 "Order" step 2; planner, after the sample seal and the step-1 smokes, before any evidence run).** The 13 smokes ran on 2026-11-09 at `code_sha` 0330a82f (descendant of the seal commit 0c68b78a; every run started after REVIEW-9691e1b0 finished at 09:15:57Z): exit 0, stderr 0 bytes, no "killed"/"stuck" line, "Total number of histories" = `NumberOfHistoriesInRun` in all 13. Only stdout/stderr were materialized; no scorer output was read. Ledger of the runs (orchestrator scratch, `smokes/ledger.jsonl`, `printout-checks.md`); t_s = TOPAS "Execution: Real" / histories, single thread. This time includes the physics-table build at BeamOn (a fixed per-run cost; "Initialization" is ≤ 0.21 s), so t_s is slightly high, i.e. conservative.

| Ion (upper end) | IDD-EM | LAT-EM | IDD-FULL | LAT-FULL |
|---|---|---|---|---|
| He 250 (case 3) | 7.16 ms (REF-2e6aa910, 1000) | 2.46 ms (REF-41e5cd18, 1000) | 71.72 ms (REF-2062230d, 200) | 250.54 ms (REF-893f1fba, 200) |
| C 430 (case 6) | 6.52 ms (REF-ffb564df, 1000) | 2.00 ms (REF-2f873a5c, 1000) | 360.03 ms (REF-cdce6317, 200) | 1439.68 ms (REF-c8dac078, 200) |
| O 430 (case 9) | 5.08 ms (REF-83f7c733, 1000) | 1.62 ms (REF-709c1045, 1000) | 474.32 ms (REF-793fe683, 200) | 1963.19 ms (REF-171fec5d, 200) |

The step-length smoke (C 290 IDD-EM, 3 histories, `TrackingVerbosity` 1) is REF-ab109b03. Full run ids are in the ledger.

**Rule outcomes (t_s × N_part against 1200 s at the upper end).**
- IDD-EM: 179 / 163 / 127 s (He / C / O). The plan holds: 25 000 × 4.
- LAT-EM: 61 / 50 / 41 s. The plan holds: 25 000 × 4.
- IDD-FULL:
  - He: 359 s, holds at 5000 × 4.
  - C: 1800 s, so rule (i) applies: **N_part = 3333** (= floor(1200 s / t_s)), total 13 332 ≥ 1e4.
  - O: 2372 s, so rule (i) applies: **N_part = 2529**, total 10 116 ≥ 1e4.
- LAT-FULL:
  - He: 1253 s, so rule (i) applies: **N_part = 4789**, total 19 156 ≥ 1e4.
  - C: 7198 s. N_part ≤ 833 gives 3332 with four parts and 4998 with six (rule (ii), parts 18–19), both < 1e4. **STOP under rule (iii).**
  - O: 9816 s. N_part ≤ 611 gives 2444 / 3666, both < 1e4. **STOP under rule (iii).**

**The LAT-FULL STOP for carbon and oxygen is recorded as a STOP of the plan as frozen. It is preserved and is not re-graded. No case is dropped.** Floor(1200/t_s) leaves ≥ 25 % to the 1500 s per-call limit, against an estimated 3–7 % statistical uncertainty of t_s from 200 histories.

**N_part for the lower-energy cases.** Every case of an ion uses the upper-end N_part of its configuration. No scaling argument is used to raise N_part. The upper end is the most expensive case under each of the three scalings below (all ratios ≤ 1), so every part stays ≤ 1200 s.
- (a) IDD-EM: cost ∝ R/h_z (number of depth bins). This is supported by the step smoke (97.5 % of primary steps are bin-limited) and by the measured 1.9–2.2 µs per bin at the three upper ends.
- (b) LAT-EM: cost ∝ R × (0.02 mm/h_x).
- (c) FULL with the 0.05 mm cut: cost ∝ N_δ/h. N_δ = ∫ z²/β² (1/T_cut − 1/T_max) dx is the delta-electron count over the track (T_cut ≈ 58 keV assumed; not printed in these logs), and h is the scorer bin.
  - Ratio to the upper end, He 50/150/100: 0.24/0.49/0.50; C 80/290/150: 0.43/0.61/0.50; O 100/350/200: 0.67/0.78/0.76.
  - Projected largest lower-energy part: O 350 IDD-FULL ≈ 930 s.

**Bundle consequence.** The committed bundles carry `NumberOfHistoriesInRun` = 5000 for IDD-FULL and LAT-FULL. Before any of these evidence runs, the values must be changed by a **reviewed bundle regeneration**: `gen_r2_bundles.py`, the new bundle-set digest and appendix A, and the P9 rebinding by the orchestrator. Parts keep their numbers and seeds (`20619000 + 100·part + case`). The EM configurations are unchanged.

**FULL cost is delta-electron transport, not nuclear (measured).** The excess cost per mm of range, (t_FULL − t_EM)/R, scales across He, C and O as z²/β², the delta-production rate:
- excess / (z²/β²): IDD 0.0160 / 0.0169 / 0.0168 µs·mm⁻¹ per unit; LAT 0.061 / 0.069 / 0.070.
- A geometric nuclear cross-section scaling, (A^{1/3} + 16^{1/3})², would give He:C:O ≈ 1 : 1.37 : 1.51. Measured: 1 : 6.8 : 12.0 (IDD) and 1 : 7.2 : 13.0 (LAT).
- A two-term least-squares fit (z²/β² plus nuclear) returns a nuclear term of ≈ 0 (slightly negative). The nuclear share is therefore not resolved, and is small at these t_s.
- Mechanism: in the FULL logs, e⁻ msc is GoudsmitSaunderson with StepLim = SafetyPlus, Skin = 3. Every delta electron above the cut is tracked through the scorer planes: 0.1 mm depth planes (IDD); 0.02 mm x and y planes plus 1 mm depth planes (LAT), which are two parallel worlds, `Water_6000x1xL` and `Water_1x6000xL`.

**Module and process checks (all 13 logs).**
- UrbanMsc for alpha and GenericIon in EM and FULL (GenericIon: StepLim = Minimal, Rfact 0.2, Emax 10 GeV). Step functions printed: light ions (0.1, 0.02 mm), general ions (0.1, 0.001 mm).
- All six FULL logs show the hInelastic QGSP_BIC_HP thresholds block and the "Radioactive Decay Physics Parameters" block. This resolves the §9 `g4radioactivedecay` risk: accepted.
- **QMD acceptance has positive evidence.** The HADRONIC PROCESSES SUMMARY printed at the default verbosity lists `QMDModel: 100 MeV/n ---> 6 GeV/n` for ionInelastic (GenericIon), alphaInelastic, He3Inelastic and the deuteron, next to Binary Light Ion Cascade (≤ 110 MeV/n) and FTFP (≥ 5.99 GeV/n). It appears 5 times in each of the six FULL logs. This corrects the orchestrator's printout check ("QMD acceptance not printed").
- Report-only fact: GenericIon lists no hadElastic process (only ionInelastic and RadioactiveDecay); alpha has hadElastic (hElasticLHEP).

**Step-length result (resolves the §4.3/§9 risk for IDD runs).** REF-ab109b03, C 290 MeV/u IDD-EM, 3 histories:
- 5046 primary steps (1682 per history; R/h_z = 1632).
- 4919 are limited by the parallel-world scorer volume `Water_1x1x2500` (97.5 %), 118 by ionIoni, 3 by CoupledTransportation (the entrance), none by msc.
- Median step 0.1000 mm, mean 0.0976 mm.

TOPAS scoring bins therefore limit the primary step: the IDD runs transport the ion with ≈ h_z steps. This affects only laterally integrated quantities (MCS detour ≤ 0.01 % of R80, §5.1 (a) budget), which is the reason for the separate runs (§4.3).

LAT runs were **not** measured. Their 0.02 mm planes are parallel to the beam and are crossed only through lateral motion, every ≈ h_x/|θ| of path (estimate: ≈ 2–4 mm at the upper ends at mid-depth, shorter at the low ends with h_x = 0.01 mm). The 1 mm depth planes cap the step. The expected LAT primary step is therefore ≈ 0.3–1 mm, not 0.02 mm. That is inside the 0.1–1 mm span on which the TL band was built (smallest Urban/Highland-FE factor 0.835 at 0.1 mm, P2 §6.3). **The TL band is unchanged.** Measuring the LAT step regime is part of the proposed smoke set below (seed 20600427).

**Projected REF total under rule (i)** (scalings (a)–(c); R from §3; the R-only scaling in parentheses):

| Item | Hours |
|---|---|
| IDD-EM | 1.19 (0.95) |
| LAT-EM | 0.32 (0.31) |
| IDD-FULL (He 5000, C 3333, O 2529 × 4 × 4 cases) | 8.55 (5.79) |
| LAT-FULL He at 0.05 mm, N_part 4789 (superseded if the proposal below is adopted) | 2.32 (2.24) |
| V5-SS | 2.94 |
| 13 smokes run | 0.26 |
| **Subtotal without C/O LAT-FULL** | **≈ 15.6 h (12.5 h)**, cap 16 h |

C/O LAT-FULL at 0.05 mm with ≥ 1e4 histories would add ≥ 2 × 4 × 1e4 × (1.44 or 1.96 s) × Σ scale ≈ 25–30 h. That is beyond any reading of the cap.

**PROPOSED amendment R2-SIZE-1 (dated 2026-11-09; NOT IN FORCE until independently reviewed and committed before any LAT-FULL run).**

Options for resolving the C/O LAT-FULL STOP without dropping a case:
- **(A) LAT-FULL with `CutForElectron = 1 km`, the cut already used in the EM configuration. Recommended.** The seven modules stay verbatim; `CutForAllParticles` 0.05 mm stays for all other particles; IDD-FULL keeps 0.05 mm, so TD-FULL is unchanged.
  - Gating observable. The TL-FULL observable is the generation-0 core σ (§5.3). Delta electrons are not generation 0, and the primary's msc and step functions do not depend on the electron cut: `Use cut as a final range` = 0 in the logs.
  - Restricted vs total loss. Only the primary's deposit changes, from restricted to total loss. At one depth that is a nearly uniform lateral weight (all primaries have nearly the same β), so σ_core is unchanged to ≪ 1 % (estimate).
  - Like-for-like gain. The deposit becomes like-for-like with ionmc's local deposition, and the TOPAS-internal ratio σ_FULL,primary/σ_EM isolates the nuclear and hadronic effects (identical EM cuts).
  - Lost (report-only items only). The delta-electron contribution (T_max ≈ 1.2 MeV at 430 MeV/u, ranges up to ≈ 5 mm) to the FULL all-particle core σ, the halo fraction and the ±20 mm second moments. The delta production of secondary fragments is also lost. These become "seven modules, delta electrons deposited locally" quantities and must be labelled so in the ledger.
  - Cost. The z²/β² scaling above predicts that most of the LAT-FULL excess disappears; the nuclear part is not resolved by the data. Planning bracket for t_s at the upper ends: 3–50 ms. **To be measured.**
- (B) Coarser lateral bins or a narrower extent for LAT-FULL. The ionmc samples hold `latx`/`laty` on 120/h_x bins from −60 mm with their edges, so exact unions (k·h_x) and sub-ranges are re-binnable; nothing else is.
  - Cost gain. A narrower extent hardly lowers the cost, because the crossings happen near the core. Coarsening by k = 5 (0.1 mm) cuts the plane density by 5, but the GS SafetyPlus cost per electron near planes is not proportional. A further 3–5× reduction (C ≈ 300–500 ms) is uncertain and may still fail (iii) for O.
  - Science cost. The TL-FULL bins and resolvability rule would differ from TL-EM: f = 0.25 points at the C/O upper ends may become report-only. The primary step regime would also change between TL-EM and TL-FULL, confounding σ_FULL/σ_EM.
  - Not recommended.
- (B′) A depth-windowed fine scorer (fine x/y bins only in 1 mm-edged windows around f·R80). It is exact on ionmc's 1 mm slabs and gives ≈ 5–10× fewer crossings. But it needs either a data-dependent bundle (z_k depends on R80_TOPAS-EM) or windows wide enough for the planning-R uncertainty, plus eight extra parallel worlds. Not recommended.
- (C) Keep the 0.05 mm cut and fine bins and accept the budget: ≈ 25–30 h more (projection above) at ≤ 1200 s per part. That means ≈ 70–100 more host calls and a REF total ≈ 40 h. This is a budget decision for the orchestrator; it is not justified while the gating observable does not need delta-electron transport.
- (D) Option from the plan's own structure: drop the all-particle `LatX`/`LatY` from LAT-FULL. This gains nothing: the cost is plane crossings, not scoring. (A) is the plan-internal answer, because §5.3 already chose the generation-0 profile and declared the all-particle FULL profile report-only.

**Recommended procedure (A); bounded.**
1. Reviewed bundle regeneration (follow-up implementation task; list below):
   - all 48 LAT-FULL evidence bundles (12 cases × parts 12–15, same names and seeds) gain `d:Ph/Default/CutForElectron = 1 km`, with comment and `case.json` "cuts"/"physics" text saying so. Applied to every case, He included, so that TL-FULL and D_ext use one configuration per row;
   - IDD-FULL `NumberOfHistoriesInRun` is set to C 3333 and O 2529 (He 5000) for all four cases of each ion;
   - new execution-only smokes in the rehearsal family:
     - `he-water-250mevu-r2-lat-full-ecut1km-timing-smoke` (**20600424**);
     - `carbon-water-430mevu-r2-lat-full-ecut1km-timing-smoke` (**20600425**);
     - `oxygen-water-430mevu-r2-lat-full-ecut1km-timing-smoke` (**20600426**), each 200 histories with the LAT-FULL scorer block;
     - `carbon-water-290mevu-r2-lat-em-steps-verbose-smoke` (**20600427**, 3 histories, `TrackingVerbosity` 1, LAT-EM scorer block; measures the LAT primary step regime, execution only).
   - Total 209 bundles, 418 files, new digest. 20600424–20600427 appear nowhere in the repository except as the end of the current `SMOKE_RANGE`.
2. Run the four smokes, one at a time. Each must pass the run-validity and seal checks (R2-TOPAS-4). Read stdout only.
3. Dated re-sizing note: apply the unchanged §6 step-2 rule to LAT-FULL-A (N_part 5000, parts 12–15; then (i) with total ≥ 1e4; then (ii) parts 18–19; then (iii) STOP), with the scalings above for the projection. If the projected total exceeds 16 h: V5-SS takes its 4e4 fallback (−0.98 h), then STOP for the orchestrator's budget decision.
4. Frozen and unchanged: the modules, the TD-FULL configuration and cut, all tolerances and bands, the seeds formula, the part blocks, K = 4, the sealed ionmc samples, and the evaluator's estimators. The 0.05 mm LAT-FULL smokes, their C/O STOP and the He rule-(i) outcome stay recorded as history.

**Budget with (A).** Projected total ≈ 13.3 h + LAT-FULL-A. Break-even: with the conservative scaling, the cap holds for t_s(A) ≤ ≈ 85 ms at the upper ends (2e4 per case). The bracket upper bound (50 ms) gives ≈ 1.6 h, so ≈ 14.9 h ≤ 16 h; at ≤ 10 ms, ≈ 13.6 h. The V5-SS fallback is projected not to be needed.

The 16 h cap is a planning value of this plan (acceptance plan §E). It is not an operator directive: the operator's directives fix single-process execution, one heavy job, logging and the seven modules, and set no REF wall-clock cap. Raising it is the orchestrator's prospective decision with review.

**Decisions for the orchestrator.**
- (1) Adopt R2-SIZE-1 (A) after independent review. The cut is the plan's own choice: the directive fixes only the module list. The ledger wording becomes "lateral agreement with the specified seven-module TOPAS configuration with delta electrons deposited locally (CutForElectron 1 km)".
- (2) Whether He LAT-FULL also moves to (A). Recommended: yes, one configuration per row.
- (3) Whether the cap is raised, should the smoke exceed ≈ 85 ms.
- (4) Whether K = 6 under rule (ii) needs a stated σ_T inflation factor. §5 fixes t_{0.95,3}/1.645 for K = 4; with K = 6 it would be t_{0.95,5}/1.645 = 1.226, which the plan does not state. It is not needed under the current outcomes.

**Follow-up implementation task (not done here).**
- `validation/scripts/ion_em/gen_r2_bundles.py`:
  - CONFIGS histories per (ion, config): IDD-FULL C 3333, O 2529; LAT-FULL per the re-sizing note;
  - the `CutForElectron = 1 km` line and its texts for lat-full;
  - SMOKE_SEEDS plus the three ecut1km smokes and a second verbose smoke (C 290 lat-em, 20600427);
  - the "part p of 4" rationale stays;
  - BUNDLE_SET_DIGEST.
- `tests/ionmc/test_r2_bundles.py`:
  - `SMOKE_RANGE = range(20600411, 20600428)`, count 209 and 418 files, RECORDED_DIGEST;
  - a new test: lat-full bundles carry `CutForElectron = 1 km` and idd-full ones do not;
  - a new test: `NumberOfHistoriesInRun` per (ion, config) equals the sizing table and the totals reach the minima;
  - the smoke-scorer test must accept the `-ecut1km-timing-smoke` stems.
- Appendix A: a new dated listing; the old one is kept.
- `p9-freeze-record.md` rebinding: orchestrator.

**Dated re-sizing note (2026-11-09, after the R2-SIZE-1 smokes; plan section 6 "Order" step 2 applied to LAT-FULL under R2-SIZE-1; planner, before any evidence run).** The four smokes 20600424–20600427 ran on 2026-11-09 at `code_sha` 39cf5478 (the seal commit 0c68b78a is an ancestor; ledger `finished_at` 11:39–11:41Z, each run ≤ 24 s, i.e. after REVIEW-9691e1b0 finished at 09:15:57Z): exit 0, stderr 0 bytes, no "killed"/"stuck" line, "Total number of histories" = `NumberOfHistoriesInRun` in all four. **Their run records are unsigned schema-2 records:** the operator has not yet restarted the controlled reference service after V3-INFRA-004 (signed, timestamped run records, ed05934). They are execution and sizing records only, not evidence, and no evidence run may use such a record. Only stdout/stderr were materialized; no scorer output was read. t_s = TOPAS "Execution: Real" / histories, single thread, physics-table build included (conservative).

| LAT-FULL, `CutForElectron = 1 km` (upper end) | t_s | was at 0.05 mm | t_s × 5000 |
|---|---|---|---|
| He 250 (case 3), REF-fe3eb38a06c93dfcd4da-4dbe8e46, 200 histories, 9.035 s | 45.2 ms | 250.5 ms | 226 s |
| C 430 (case 6), REF-b07c8cb0f27e0da26d60-9903bf4c, 200 histories, 22.941 s | 114.7 ms | 1439.7 ms | 574 s |
| O 430 (case 9), REF-9433e888c8a566eecc8b-77d7385e, 200 histories, 22.860 s | 114.3 ms | 1963.2 ms | 572 s |

**Rule outcome (unchanged section 6 rule, N_part 5000, parts 12–15).** He 226 s, C 574 s, O 572 s, all ≤ 1200 s. **The part plan holds for LAT-FULL at 5000 × 4 for every ion; rule (i) and (ii) are not used.** The He rule-(i) value 4789 (0.05 mm) is superseded. **The committed LAT-FULL `NumberOfHistoriesInRun` stays 5000: no bundle change.** The C/O LAT-FULL STOP at 0.05 mm stays recorded as history.

**Forecast check (preserved).** The cost fell by 5.5× (He), 12.6× (C) and 17.2× (O). For C and O it is above the planning bracket of the proposal (3–50 ms) and above its break-even (≈ 85 ms). The proposal's statement "the V5-SS fallback is projected not to be needed" is therefore not borne out under the conservative scaling (projection below).

**What the remaining cost is (report-only).** Excess over LAT-EM per mm of range, (t_FULL-A − t_LAT-EM)/R: He 0.112, C 0.364, O 0.486 ms/mm, ratio 1 : 3.3 : 4.3. This lies between the geometric nuclear cross-section ((A^{1/3} + 16^{1/3})², 1 : 1.37 : 1.51) and cross-section × projectile mass number (1 : 4.1 : 6.0). It is consistent with nuclear secondaries (fragments, nucleons, gammas) now dominating the cost; no attribution beyond that is claimed.

**N_part for the lower-energy cases.** Every case uses the upper-end N_part (5000). Scaling (b) of the earlier note (cost ∝ R × 0.02 mm/h_x; for the delta-free LAT-FULL-A this is the conservative choice: the nuclear-interaction count per primary is ∝ R and the 0.01 mm bins of the low ends double the plane density) gives ratios to the upper end He 0.12/0.42/0.20, C 0.11/0.53/0.17, O 0.17/0.72/0.28 (low end/standard/held-out). Largest lower-energy part: O 350 ≈ 411 s.

**Final part table of the evidence runs (K = 4, seeds `20619000 + 100·part + case`, every case of an ion uses the same N_part).**

| Config (parts) | He (cases 1, 2, 3, 10) | C (4, 5, 6, 11) | O (7, 8, 9, 12) | Total per case |
|---|---|---|---|---|
| IDD-EM (0–3) | 25 000 | 25 000 | 25 000 | 1e5 |
| LAT-EM (4–7) | 25 000 | 25 000 | 25 000 | 1e5 |
| IDD-FULL (8–11), 0.05 mm | 5 000 | 3 333 | 2 529 | 2e4 / 13 332 / 10 116 |
| LAT-FULL (12–15), `CutForElectron` 1 km | 5 000 | 5 000 | 5 000 | 2e4 |

Parts 16–19 stay unused. These are exactly the committed bundles (frozen manifest `validation/results/transport/v8r2-topas/frozen-bundles.json`, digest 6d1ae3a0…).

**LAT step-length result (REF-8f7842c1487b5df484aa-31f28612, C 290 MeV/u LAT-EM, 3 histories, `TrackingVerbosity` 1; parallel worlds `Water_6000x1x250` and `Water_1x6000x250`).** The limiting plane of each primary (Track ID 1, Parent ID 0) step was classified by the replica change of the scorer touchables printed with the step.
- **999 primary steps, 333 per history** (296–362). The IDD run had 1682 per history.
- Limits:
  - parallel-world planes 82.9 % (828):
    - 1 mm depth planes 48.9 % (489). They are reported as `Water_1x6000x250` or `Water_6000x1x250`, which tie at the same z.
    - x planes 22.0 % (220), mostly reported as `CoupledTransportation`.
    - y planes 11.9 % (119), reported as `Water_1x6000x250`.
  - ionIoni 16.5 % (165), all at z ≥ 155 mm (kinetic energy ≤ 630 MeV, ≈ 53 MeV/u; the last ≈ 8 mm);
  - msc 0;
  - the entrance and end steps 0.6 %.
  - The pencil starts on the x = 0 and y = 0 bin boundaries. As a result, the first 1 mm step of each history changes the x, y and z bins at once, and is followed by a ≈ 2e-9 mm y-plane step (report-only; the ionmc samples use the same edges).
- Step length, all primary steps: median 0.446 mm, mean 0.491 mm, max 1.000 mm. 25 % of the steps are full 1 mm steps.
- Plateau (z ≤ 150 mm), 245 steps per history:
  - median 0.634 mm, mean 0.614 mm, quartiles 0.31 / 0.63 / 1.00 mm;
  - path-weighted geometric-mean step 0.74 mm;
  - 0.6 % of the path is in steps < 0.1 mm.
- At the TL depths (± 5 mm around 41, 82, 122, 147 mm, i.e. f ≈ 0.25–0.9 of R ≈ 163 mm), median step 0.93, 0.78, 0.52, 0.35 mm. The steps shorten with depth because the lateral-plane crossings increase as the angular spread grows.
- Mean path between plane crossings per lateral axis: 2.2 mm (x) and 4.1 mm (y); 1.4 mm for both axes together. This is in line with the earlier per-axis estimate of ≈ 2–4 mm.

Comparison with the IDD run: there, 97.5 % of steps are limited by the 0.1 mm depth planes, median 0.1000 mm. **The LAT runs transport the primary with steps ≈ 0.3–1 mm, ≈ 5× longer than the IDD runs.** This confirms the earlier expectation (≈ 0.3–1 mm). Consequences:
- The LAT steps lie inside the 0.1–1 mm span on which the TL band (a) was built (Urban/Highland-FE 0.835–0.89, P2 §6.3; the band uses the smallest factor). **The TL band is unchanged.**
- LAT-FULL uses the same two parallel worlds (C 430 ecut1km log: `Water_6000x1x450`, `Water_1x6000x450`, all four lateral scorers). TL-EM and TL-FULL therefore share the step regime, and the TOPAS-internal ratio σ_FULL,primary/σ_EM is not confounded by it.
- The step regime differs between IDD and LAT runs by design (§4.3). It enters the IDD runs only through laterally integrated quantities.

**Process list (ListProcesses).** The C 430 ecut1km log (REF-b07c8cb0…) prints the "Registered Physics Processes" block with `Decay` (log line 73) and `RadioactiveDecay`, besides ionInelastic, alphaInelastic and hadElastic. It also shows UrbanMsc for GenericIon (StepLim = Minimal), five QMDModel lines and the Radioactive Decay Physics Parameters block. The executed module set, `g4decay` included, is now verifiable from the printout of every run that requests the list.

**Projected REF total (same scaling assumptions as the 2026-11-09 sizing note).**
- t_s from the smokes at the upper end of each ion, scaled to the other cases:
  - IDD-EM ∝ R/h_z;
  - LAT-EM and LAT-FULL-A ∝ R × 0.02 mm/h_x;
  - IDD-FULL (0.05 mm cut) ∝ N_δ/h_z with the ratios of that note.
- In parentheses: cost ∝ R only.
- V5-SS as designed (10 596 s).
- Smokes: the wall time actually used.

| Item | Hours |
|---|---|
| IDD-EM (25 000 × 4 × 12) | 1.19 (0.95) |
| LAT-EM (25 000 × 4 × 12) | 0.32 (0.31) |
| IDD-FULL (He 5000, C 3333, O 2529 × 4 × 4 cases) | 8.55 (5.79) |
| LAT-FULL-A (5000 × 4 × 12) | 2.97 (2.87) |
| V5-SS (6e4 per case) | 2.94 |
| 17 smokes run | 0.28 |
| **REF total** | **16.25 h (13.12 h)**, cap 16 h |

**Pre-registered order applied.**
- The conservative projection exceeds the 16 h cap by 0.25 h (1.6 %).
- Per section 6 step 2, V5-SS takes its pre-registered 4e4 fallback (7064 s ≈ 1.96 h): **saving 0.98 h → projected REF total 15.27 h (12.15 h), 0.73 h under the cap.**
- **No STOP and no budget decision remains under the plan as frozen.** No case is dropped.
- Under the R-only scaling the cap holds without the fallback (13.12 h).

**Framing for the orchestrator (a prospective option, not required by the plan).**
- The fallback is triggered by a 0.25 h overrun of a projection. Its two scalings span 13.1–16.3 h. Its largest term, IDD-FULL, rests on the N_δ/h scaling (8.55 h against 5.79 h).
- The 16 h cap is a planning value of this plan (acceptance plan §E), not an operator directive.
- Raising it to **17 h** before any evidence run would keep V5-SS at its 6e4 design (16.25 h projected, 0.75 h margin), at the price of ≈ 1 h more REF time.
- Keeping 16 h applies the fallback (V5-SS at 4e4, its own pre-registered minimum before any case drop).
- Planner recommendation: **raise the cap to 17 h**. V5-SS is the only row that resolves the charge dependence b(z) at the 3 % level (§1, §5.3). The fallback widens its statistical uncertainty by √1.5 ≈ 1.22, which is not worth trading for a 1.6 % overrun of an upper projection. If no cap decision is recorded before the first V5-SS run, the fallback applies as written.
- Either choice is recorded by a dated entry before the first evidence run.
- Host calls: 192 evidence parts, 17 smokes already run, and the V5-SS parts (16 at the 6e4 design).

## 9. Not established, risks

- Whether TOPAS scoring bins limit steps (OPEN; step-length smoke). It moves the Urban width by up to ≈ 6 % (P2
  §6.3) inside the TL band, and the cost.
  **Dated resolution (2026-11-09, sizing note §8):** YES for IDD runs. REF-ab109b03: 97.5 % of primary steps are limited by the 0.1 mm depth-scorer planes (median step 0.1000 mm). This affects only laterally integrated quantities. LAT runs are not measured. Their primary step is expected at ≈ 0.3–1 mm (1 mm depth planes; 0.02 mm lateral planes crossed only through lateral motion), inside the 0.1–1 mm span of the band. TL band unchanged. Proposed LAT step smoke 20600427 (R2-SIZE-1). Cost consequence: the FULL cost is dominated by delta electrons crossing the scorer planes, which leads to the C/O LAT-FULL STOP (§8).
  **Dated resolution (2026-11-09, after the R2-SIZE-1 smokes; re-sizing note §8):** LAT runs measured (REF-8f7842c1487b5df484aa-31f28612, C 290 MeV/u LAT-EM, 3 histories):
  - 333 primary steps per history (IDD run: 1682).
  - Step limits: scorer planes 82.9 % (1 mm depth 48.9 %, x 22.0 %, y 11.9 %); ionIoni 16.5 %, only in the last ≈ 8 mm; msc 0.
  - Plateau step: median 0.63 mm, mean 0.61 mm, max 1.000 mm; 0.6 % of the path is in steps < 0.1 mm. Median at the four TL depths: 0.93–0.35 mm.
  - This is inside the 0.1–1 mm span of the band: **TL band unchanged.**
  - TL-EM and TL-FULL use the same two parallel worlds and therefore the same step regime.
  - The risk is closed for both run types.
- The Urban deficit for ions (10–30 %) and its energy trend are scratch estimates (P2 §6.3; the Urban emulator was
  withdrawn in r5). If the actual deficit exceeds the band (e.g. at the low ends with very short steps), TL (a)
  fails on the reference's structure; that is recorded as FAIL with attribution, never re-graded.
- `g4radioactivedecay` acceptance by OpenTOPAS 4.3 is unverified (timing smoke; rejection = STOP, the operator's
  configuration cannot be changed by the planner).
  **Dated resolution (2026-11-09):** accepted. All six FULL smokes print the Radioactive Decay Physics Parameters block, and the hadronic summary lists QMDModel (100 MeV/n – 6 GeV/n) for GenericIon and alpha inelastic. Exit 0 (§8 sizing note).
- FULL cost (fragments, neutrons crossing bins) is uncertain by a factor of several; the sizing rule handles it.
  **Dated note (2026-11-09):** measured. The cost is dominated by delta electrons, not fragments (§8). The rule gives (i) for C/O IDD-FULL and He LAT-FULL, and a STOP (iii) for C/O LAT-FULL. Proposed resolution R2-SIZE-1 is not in force. **(2026-11-09, later)** R2-SIZE-1 is in force (acceptance plan, orchestrator decisions of 2026-11-09 after the re-sizing note; Codex REVIEW-841e82ae); the measured 1 km-cut costs and the final part table are in the section 8 re-sizing note.
- ionmc fine lateral tallies (0.01–0.02 mm projected grids) are within the scoring design (grids independent of the
  transport grid, up to 4 grids) but their per-step cost is unmeasured.
- σ_R/R planning values and the FULL A_frag derivation are estimates; A_frag is frozen as a number now.
- `validation/reference_cases/README.md` (bound in P9 revision 1) is not updated by this record; a layout entry for
  the r2 bundles can be added by the orchestrator in a later commit.

## Appendix A. Bundle hashes (sha256, `sha256sum`, at preparation)

The listing below is generated by `sha256sum` over every `case.json` and `input.txt` of the 205 r2 bundles, sorted
by path. Bundle-set digest = sha256 of this listing block's text (lines joined with "\n", final newline included).

Bundle-set digest: `a8a1c27d0177d38e763344248753f0cb9792bbf2277158622b41bcbf9d8f423a` (410 files).

```
0d05020abee11f2167ef97b0cbcf77c088c9342e3899a76b43e113b8b3cb9054  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part0/case.json
a739636613115825667cf652568d4fa549f57779e96fd91eb5333a2313c28f6d  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part0/input.txt
e88ac3408f05dc31e2c4a3d35b2d59c5de77017e0f0ce50e7921c4c26f386639  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part1/case.json
1e2cc994f70282560c4a5371a979aaf1a83a584c5c8878210e7a9959d5513007  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part1/input.txt
acc80727a6c32e4542725ce742e4eee7170f9bd4202df5109f30785b8c2afcc0  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part2/case.json
d7dd335a82bcf952173a2fdd491174b71168e93197357036e69af0306ec54930  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part2/input.txt
71a168f5a2b9362f9e83882db15f325ed11abd4ddc731444b535816773073769  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part3/case.json
be3f4c371f65e99f07901b2fefffe9a50833e0672b6829dacb432b22f9e16c15  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part3/input.txt
385bd83d6976aef52914c0fb93601e01af1bb3ff97a1a2afa6a28653ae555b0c  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part10/case.json
b3fff7e00264469e5d6536dc40e3ad0ec48a90c4908f78e50d786271aa3a66ac  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part10/input.txt
792c7680b4c45e0e9042f009a5735a999b660411d751c3354d659035454f9518  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part11/case.json
410f5a4ff2f9be2f6350566547219fc60e46ff2a8b1230ca1429941f5bd3882e  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part11/input.txt
b9f6153d312601cbe2d02674eb28cbee2ba610aef75a41ed5b4ce70c19e4b60e  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part8/case.json
3e5c843cf0ce95f3ac2064b8de8e38d1b1c4c1a77e2fff40b646abdd5bae5e37  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part8/input.txt
56b7309ed8e72c3ae28397d71c4d7d5132992168f6680259d65bd82963a94845  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part9/case.json
2e73d8d68e79e70dd3c060909c30731b66a27e3909d8889fb5c1deda390b5882  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part9/input.txt
c7da3edb0d3f786124007891b6ed731100e8fde449f29c7f37d84b9c50b6c317  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part4/case.json
bc0df5ddeeae5f69cf771b4792cc7d9cb50343482a49be06a5550cd509a0074e  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part4/input.txt
c42d407bccea453fcdeb3bb09100a3944280f59674ed52d997fe00ecdf1d1e6e  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part5/case.json
c062cc352472cd4b369c3c6addcb09abc28c33726d0abc9396aeed1bee5eecde  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part5/input.txt
021bfab9f143b21afe4b00615f887b17987aa2b4428e5c5f74794d6bb476cb43  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part6/case.json
b1c8d594a1a4491cf7d7b10337d2b29b7046f3127822f025a8a618567244c498  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part6/input.txt
d92698482c6d8ad40c753d04f884dc000da71604e8e8dfc0cbf2c58a2446be98  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part7/case.json
8742b333bfb5bcdc8038edcee488f4b7d8a68a32480fa053f5b195ff339621a7  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part7/input.txt
65265dd880a8420cd1b8be927b03ff86f035c225c29493d53fe6f48b190cef76  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part12/case.json
12637dd774241a42f5013a291b78f87c4bbb535f0fca55d93939482322862e19  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part12/input.txt
f07166666c28b3d7e631ab4fef55bd0448b8d44020d9da7badfda85a156faedd  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part13/case.json
cfe0345ebf44cd5c11a64cde23a34fe35e25f7732eb05b28e256a5fde9834d59  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part13/input.txt
c1dffd7662e186900e039515933136b35a9e7ceaa59e8a161b237cc905aacece  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part14/case.json
6824cef85804663bd120e4db19f840f90550491ce801d1c302c7c81de465d92d  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part14/input.txt
eda50eb0de16b6282337ebc5d6980465ccb3fcdf70bd5649d34390115527bef6  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part15/case.json
9b6b4fee49afbb5b16625cdd8e48a1b23a83df61691abd444ef3adea08f3a738  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part15/input.txt
d5b03125879e78930de5c9a25ddb6bb01eb408412e1fa2180642a9316ada44e7  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part0/case.json
41ab2a9038fa6dae2954b35c470095dcdcfe168bc66db211c27ab9b7ae8bc6bd  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part0/input.txt
fae636cbfb79689d33fa612e7cb46473646a5dd5475f7c591252601514528f58  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part1/case.json
5f23c64c8111c012908b8274eff923d2374b42e7fc3a60f855421c4cca2d26e1  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part1/input.txt
7b3a868ff939f8687c73ea930b003435637a33563437bbfc6db19ed693f9a4e8  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part2/case.json
c76ae4a6b82167ba92b80599f2e64df264a4ccbe296662a792b2c40a2c32527b  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part2/input.txt
64f1945b14ec9000446675870263a8ca2e782bec4a0e24f1186b9c4ed922fbbf  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part3/case.json
678eb414cb8394995a4e638d9f97c44ac75c7a894b929d8b1b0a2785852eb128  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part3/input.txt
a35c696049bf3139520693c381791a74689752551d9216ff4133e7e4aaf63f11  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-steps-verbose-smoke/case.json
a5f178bfc114719b14fd9032715896fbdfce0215f507372c5bfedea087aab827  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-steps-verbose-smoke/input.txt
a56caefcece0d2978149b85f2a446aa35cb93ef27c0c0edcce930d732b2973df  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part10/case.json
b30bf9a79778d615d19909dd5542e0bb7a982cb4dacb8803fe81fb0d7ac2487e  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part10/input.txt
9735c08b68c351212aecbfb928448acd580a9882202a4069a4d9955efc7eca23  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part11/case.json
5290146582dcda1fbfa3be3d52b1695ad7e79ce2fa658fe68fb21f63d47663f1  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part11/input.txt
f4c71ddd081ca5c0fa277ef1a3139ed1688d20f130c9603dae036bb87852d19d  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part8/case.json
3c2f031329bb4ec27bd5274a55ac803a0b0ec84cc8843fa2895317461361f54b  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part8/input.txt
cccd91a5398831b4f80d02d79bed6802da428cbd5a7c406736e2f06629fc78e4  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part9/case.json
24ff771c354ca397ed4da86753f6efd1cdb860ea5cd03b5c1cff2b33b746a4ca  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part9/input.txt
5d46532272b0fc9d0d92765066a7bc09180b4573abc4ebc63bf483401f1942ab  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part4/case.json
47a389e66b4948114bade7fd250ecea7f55b396353b4a3f3969c6fca0e7e3b4b  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part4/input.txt
543753d65dc76196f31643f42333ddc31b33a84172f93100c319aa7c88f86ca2  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part5/case.json
d31c877d60de66b47ebd69309360706e404e1729ec9412a80dd3fbb3850e46d1  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part5/input.txt
d62d36eedd486aaddd4701a4343db59c1168743f18088f3d9042124fd109b4aa  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part6/case.json
e0c90ea4a99116e6059fb5ca2ec7dd3ead937fd6a5fef052b9ff8b8fc07a9124  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part6/input.txt
5d3d82e1e3a94e7f27229121bd4a74a7ff496fad41df93558a665cfe180a10c6  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part7/case.json
41ded92c17eb2b825ec6097bdb4e6a636abb7ab261513a316385c600e181ed44  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part7/input.txt
723026c8003caf91a5c59de07e59a4cb1fad611511ca75a47c3c080215bda8d1  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part12/case.json
cae4e0053181be33425975c12dec78d6c361320578a16703b8ecdd372ff923c1  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part12/input.txt
cae10efa493a9a5fcc219dab92cf87405e8adaa7f59a5e0fc3826a7b0e0c97a6  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part13/case.json
17b0ed9cf035138542364db3a51ee72f59d8b65a3f17adba304c3bbb6675aa28  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part13/input.txt
91f2f21dd257a7ce5e6144e4af31a7777a8ecbf094750f66e333faacb405f677  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part14/case.json
c0e6c453443b4d48ee1755910fb7fd9be4fe0846435ca7a6ac5d31fb9628845d  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part14/input.txt
e95890d03f160588cd29bec6f7aa315d54a817ebb4bc1b4c779988990fb1787a  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part15/case.json
4c82f0f51ca026e2f3ab3daef3e02644eb930410982b0944bbc445baa2e45fbe  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part15/input.txt
a0ff91b6a1911fa35744f9b0d15850cf87a82617be71de75a9bb216950f18f63  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part0/case.json
e8e4bbd231a76ed89c79bdfab7a79d6bff295603f7e781729b782cdfe54854d7  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part0/input.txt
02477f2e202070672b638ab738a1fea8c38b70bdc127fdaab43f62f417ab2197  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part1/case.json
d15bfd749939cb5fb9b7b0683779f1245f3aa6bdfff9ccb8f42a9c2b65a9004e  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part1/input.txt
5915e4a195f550f6fd80fad02354b33f134948c0a89769a7e4fd2fdf2eae8b00  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part2/case.json
09dbd66c89f97e49731a82ee789826ffb784212a5d4d5f08504a48d04f54d6ee  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part2/input.txt
f5e750785c6f3af7c023348235441a02d079dc7653ded2379f27f0d3f8e61e1b  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part3/case.json
41851c95bcf228deaeb364c8a7aaa4aca2b5d507f691d452919905251aecf58e  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part3/input.txt
bc1b51aaaedcea7b5c9faef52dbfbbb2eeb8d145aada80a6953ebabc1a1980e1  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-timing-smoke/case.json
e764f08fd7aa63b82b206194fb08f657da0ea0038112fdf9694d5971728dd140  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-timing-smoke/input.txt
382c3c91647efdc1314c68a1285cb9948fc833d91cd0542db6e22be7b48ec39b  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part10/case.json
8ed41c2e9e1d469fd75671de40622a98c37a003d758506d376ab6ef48ff796b7  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part10/input.txt
fe0014d599403d6f14358d5e5ff68e715c2ee4ecc9a8b8b90f0bd661f6f7ef8e  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part11/case.json
3d346adeee9624d394f5182545afa233b7bd0d0c5405d3ebe51a100f41463292  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part11/input.txt
bf16150b9470495d1a640cac3c5d1b3812fd2130dbe48465f5d1c545593c4e7c  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part8/case.json
8de503e2d54df10a1eae4b8b4e3cc5eb46ab48e143ce592a3780b42e822e9d03  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part8/input.txt
069c258f44848b84485bc606ce811c0d8853321a5320850d1d3cfd69df288e0c  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part9/case.json
68ef690e282c4b9e5ee8a4acaa82bde1195e0431ed1ea3bab0b40208ab974adc  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part9/input.txt
abe48fa2c7ed2bb16b94304d145249442383b7f15e3460e28f6cf0c1f43343bd  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-timing-smoke/case.json
62291644f88c80cb77ed3a50b3282276994053dbe47b20d7287b7edc09d4ef78  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-timing-smoke/input.txt
a7975cb46a40e8a5eb5b36a352db39aceedf2b0e7d7145708cc565fa452a7d9b  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part4/case.json
a32e18e87d30f184c065d0f93b4622deb2e38e4d7bf26748637d9b8d8c84b178  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part4/input.txt
ae855029b11701a7585f3282b2106141f360573563672d6ec744d41d049d0687  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part5/case.json
1acaeb718f82700479620d4c083e8bbac28c520bd2f9593dda85f50ed171c652  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part5/input.txt
4dc458b7e061532809730f0321c6430c253b330bf6c19d3da6de88fafdf66d4d  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part6/case.json
25ac6fab3bad4f91cb014e514648488913e5e03c5c2216a0f7d77647a6d372ef  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part6/input.txt
882defbf9104e46740e4326faf7d04e62989d09e2e0e35ee275dae2368786464  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part7/case.json
b73c788c01307a4ba3de511039ebe4c89b332c614d217915a758d4f3b0688caa  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part7/input.txt
34bd5eff92da781bc36b5336ee36665d2013e0d0f4621bc61d3921294e7df790  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-timing-smoke/case.json
374cff2aea2a3537ea1390f2632b998b2fe2a01fbe1275bd345cdf1f04b38f0e  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-timing-smoke/input.txt
3c9f46063fb4da211212854ae7cdccd7de02a064d137fe6cc0ca9a4b7cc21f9f  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part12/case.json
2df4287de9921460651af5b571711a482a5c21b4dd2338c4e6de5bd180757cd5  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part12/input.txt
f5ba25d1b555c04225945f9f759dacf51d3e84120b327c41856c848a1f1c6964  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part13/case.json
bf082d2aefe5f7a7342b2ed932a06c025575933bbacf4bcfe50aa48b4501a1c5  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part13/input.txt
ae0677d3cfdc6a2326d7e28f1da3919520fefe6bd1bf232f7fe5e2fec32fb74d  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part14/case.json
edf0ef74cf84741376417691f06e9f8ba8598150a2294b103f8cc9538c5876c6  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part14/input.txt
ca6be186ab24f2ce2650b51663509b190e0019e195ef3aa2382a028f8cccabd2  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part15/case.json
d54a3a6ba350191e4fd0651ae0de2bd5402b3591b4c030a5c13cbdb2017455fb  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part15/input.txt
9f5fffea2438ac0ac65a5a1c25cdfd2bf53b852fd9eff4bcabf2fe8773a82f6f  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-timing-smoke/case.json
bd7058eec0b55604119c1de3ea75a29d0cf541b4712a6e00dfd7b1a1665cf9b5  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-timing-smoke/input.txt
8896b767f450df43302d1582c079786be4577be705d7193eb0cb35e1369a44a7  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part0/case.json
9e5593efa371613ad14482748506b7b3f4ce71bdedb40bd5b40030a923c4630d  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part0/input.txt
d5ff09091bd7e48d982696dfac73f650f1e00673c9ad7a786a05530b62016ecf  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part1/case.json
4d36de807841fc6ddc741088fecc0202481eac8fbedf8e75c97b318946361a6f  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part1/input.txt
324dab6c6e9374d791c313a5b2d374f4c50d9114b6fada53411d8c0b3e669687  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part2/case.json
ff70d7b84442f825d69d4bb16941d0bb17f755234c675512d0c7aa70ed5e25ec  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part2/input.txt
0d4f2924cb29cf6901784fe4a745ae010c7fa86a61af594bce0d12f5ec7238f2  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part3/case.json
3bd234105b63338e1a36c028a9683c7b8abd08d985d5bab664c84987198b9763  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part3/input.txt
427e37f613bf6d31ec7e4c56c15404b744b3b114e67d1eeb62b9faa5fb40218b  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part10/case.json
df5f292b17996188ef0a67e36ac0cae8123ea66fb95f5d58785e70d3d3377dd1  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part10/input.txt
8600df7cd0e29ff864e9180a8d0d30abe6050466ff90d8676a6116cc4c205105  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part11/case.json
4948a141f094d47c451a9454d18b007a687c53f3580dc5c08964dae5cf139be1  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part11/input.txt
439a28c283ef4a4664028e439c33dac19fc8d35f3142649c8c226769ceab686f  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part8/case.json
076cbe352ec413a7d73af07e35d9025d91085bc0533a18d4489a4930ec369260  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part8/input.txt
2f96b810b712afa4dfa739929ce9ad94a3a8e0f54dc76d71e1259f90f125e3d7  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part9/case.json
729ff41655b7ef9933d5d8b676381f46311a225759de37c30450b5a7816b2fd6  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part9/input.txt
0eceef04c7b75b1d16d029a197580c0bb05799ee2fb353b58cdc68cee21fbccf  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part4/case.json
1ba50fd5c3d82da4d0fa5bf16d8e15bc7fb398a5c96d9f07902af9deb594b63c  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part4/input.txt
10d54312445b76fecb09e6d441bb2abe85d1d4fd4e53e56bffd9d753cdc05b39  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part5/case.json
e72c3be25b2645be2d0b648a8f2046f3a6c5ac33f07835ccf1746e792a599dc5  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part5/input.txt
44603556587b7de581d254d40852e9a071d6020ccdfb6bb5cdae7167e4c91765  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part6/case.json
358b7ca92806b83eea3ce294f5a696e725c01a952b4117f2eaa2492ffff0e8ea  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part6/input.txt
db2420538316cfd33bc14dfb6f31d7e07b4b7eee04ba97533fd0b25a443ceaa1  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part7/case.json
3c0a68125d67f85be2ea0d36fdad1564741f9aeb4a83aff0188f3f8ccce64663  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part7/input.txt
6c1aba73096cb62f37a0ce7f4af2fcd8c43024c1c5a2cc88ef137447d5dd7968  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part12/case.json
cc2431bd6ef5a7ba288b6172061c1e8e5538783d6aa487695f93a0d23a6fe6fa  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part12/input.txt
20ae45d4b17838c125cfa8994341baaaecb15a21150b6aa193838ab7098a92e7  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part13/case.json
321c3ae8d514e48c0cd4aca94c661b9b5503942e91d7e316fbfd8b69041daf41  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part13/input.txt
32b4c1b4e0655a580232df270033cc673e880c3b08becfebf50007985ad58001  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part14/case.json
9ccec95dbda406d01b36b376214c0f66ce1a80cf56df4b160e8dbec99c5ca5eb  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part14/input.txt
e9f621c272b336a6cac7058fbaa8d81f78fd24c3820a0f225db9c866fd0f6b78  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part15/case.json
c3d1b3736ac50012de60d29e2f8756ab41cdefcae2318544fe2278a5e2b2fbce  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part15/input.txt
e6447aef899b7d35ec548fcc4ce8245efe6ff184a4b499da43fb793a93368a6c  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part0/case.json
6defdbe3534608b02c88ba83262c8eeeda9db5225d7320098ac9b07d72bae3d9  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part0/input.txt
bf6c8095d9c65d0ba3c99417f1b28615082054aa85870aebeb033ac461aa1980  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part1/case.json
bc772dcb64f05e76dac60e6e12e4a577da3677b40f28df513ce07972a4ced254  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part1/input.txt
469b0eb051b9288ef9086745b4b51ffed7dc157c63de0fb39ef9f6a42a890d28  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part2/case.json
f0caf317650f2ef126a74b680c7bf044df0e10b4c7c8b5366fc4a9c0e39cf889  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part2/input.txt
c6fd29ca3bada5b9b58eab5c401f77f2410f66d7910ac0f8dd25ac8f54b98266  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part3/case.json
7a2932b3cb7016857803f1ffe4e4cee97e1ca6bb4d227eaefddb291953f79714  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part3/input.txt
ba5670ff0611233edd15747ab1b08455f068e31bcf4e48622b97b411416228ed  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part10/case.json
b6b331f4bfff5d1584d485f8f5fdccc884d7be3252673a1b318d523228ec8d3b  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part10/input.txt
0a81392df95171929d3b2c13664594ba9af499052f450999581f3609e717eb32  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part11/case.json
e7e31ed8793c2692399d392560200eec39a5e34526cbaf6803d0442c583db5a9  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part11/input.txt
272fa219e1e773878fee7085cd87acad5a053a3b999651dec4fc379e7a34ae5f  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part8/case.json
87577171f3f74b8b62e424379dc3a816fc48b32ec0b3a5bd9bc0c7a1f9ae5add  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part8/input.txt
5e5968bf3cc4b3267502fd37fa914c0eb100786da8d296b0b8dd2639e4e08cf7  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part9/case.json
89a11ef7817f062004c31880b0187b4cac68b65903a0faef6405350e60fe581e  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part9/input.txt
66cf863b9ea6727b6ec744850be81972170b89df350afde397515730d68446cd  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part4/case.json
5fd7328c07fbe70abcce9f65935d76bd40d4d65316f5f0c50ecc1c13b8373a42  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part4/input.txt
87b8e2a181a577d60f5205b56131af9e229c9378b6c23ed5bad11905268b47f1  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part5/case.json
b65df798485fdc7bc2ef183808baf18e4efa11bcff333f4bd0b4e30fa1a9d022  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part5/input.txt
0c3ad6dab6384fc3d2c8aabeb3057fe5600521adf3fac0b8608394b175d7aaab  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part6/case.json
29bc3bbcc54d0acd2bae41fb335020e7ab50163aa8f144fece64df0ab893df5f  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part6/input.txt
9994dea2d25b729edf3d4a7c3d754d471cc7fe089e7c3137493fc53461b9355e  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part7/case.json
5f517f0a19b79f9a6aaca65012e2f8d24aa917a3b05be4113201a779808fb2b2  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part7/input.txt
46d6e1253ff7a34cd3b8067793fdaac8405066e7e1e2cfcef2920ce502006af7  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part12/case.json
2b8df4e05b5964f77477ba9ebcf08f3c58e878d72ad4b7bf0dd1a6652c412233  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part12/input.txt
73012481f754d9e683cab88980154cfca02612f0f8a228337131b7948ea66d99  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part13/case.json
9ef5b1d937d1c15e83b8b94fb47e4cc2ca22a9597b7c8fd4e8277a091e4ca420  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part13/input.txt
5cf4e2f8c4d5fa32f9fa7342ec541230fc7c42be65a564dd4a56748e4c62583b  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part14/case.json
b9e8b5222c4ba1ae5c925d7e8c9f360e9923233e0db3f6db181105b670d2ccb8  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part14/input.txt
f81e7c81386d58246b9dfd8c085f62a3955cb022f7c605346f725fde42780c10  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part15/case.json
99fee2eeac111d48849831c69b17bc09ed713a98d060d4258ad0d6c6743ec80c  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part15/input.txt
1e0ced17112f261f4c1a91d69cc01df9cd82bd77a27dc2a5da7c8fbd385c799e  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part0/case.json
cea1e2eb832f4ded7579e28ceee70235e35cd33bc689d139b0eadced1b921c22  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part0/input.txt
85d57721675b2fd893f3f3db23ab5509dae036ccd19983643700be000f333bc2  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part1/case.json
e2201f92cbe75f96424864140ffba06d6df4bd6207d1d0d0f92e2daa56b7d422  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part1/input.txt
33653c422076dd5a978ceba9e329f8b24c2e5c06db0d80dc527ba831dcee61d5  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part2/case.json
42f78199174978a773f40fb15ab3517afb2de572c1dafef38e1a5c155ddc498f  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part2/input.txt
f8c9e5d9947923229cc2b10f054221af6f3304f843a4e168f92c7e4c584de4ab  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part3/case.json
c709d059800e442f6316e8d113d8d098360c273617affc29139357a87f6272f1  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part3/input.txt
23209d8eb7264fe388b1a735965fdb5e66bdd517c2ab3b5c2b497d310d732aa6  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part10/case.json
79bbe5761fa017b852b58e7a6dc370f5bdf345cc74e562e0e74285622c29aec3  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part10/input.txt
375984531b4c98e8858c06717af4194443e460590d51d53ad001248eb9f7b236  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part11/case.json
70aefa0feaa2cb27e752cce2ebfe6cfebd600b343ee2d8e8055365a516ce8b86  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part11/input.txt
cb91130d3533244534798416f5f60d82937dd109d6e93ec9c77f771c9b637a99  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part8/case.json
21918eeca4a155f2874e532e46acae0e4fad04d87d2ce61975bbadfa408c76bb  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part8/input.txt
1e96ba2fb61863eaef9799dc4de7995c1b6b06679303370a25a33a2b75967fb2  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part9/case.json
faa64686350a2727ca18e3c92bc8dceb9844d37c40c89dfce3641525b91143be  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part9/input.txt
f45d7bb098d77897465eebf300d03a8a5f1cf6b81af5d746d3327332fcb25e18  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part4/case.json
56357b81bd1f98dca22530e6c9eb5b10fe444ee40361f8e206fd0dcbfb3e69ef  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part4/input.txt
6d4a798e09457fa30b523bd97497a55442c05f97a0fa749c34ef851ca31d8d4b  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part5/case.json
f58a07c73a4bb1cb0869bb5cefbbf98461b90599d96f27b516f4f6e24728fc29  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part5/input.txt
2807e73e12eb891f818f913c44c629abd9db509278b36f72b797bb06c5e20089  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part6/case.json
55214145802c4c93f4a73e4df43b8f2adf3919ef8505b327773df9bf25c3c45f  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part6/input.txt
0f02c3f02d0c99d7190885486cd36b1636dbdb6e92da39867bba65b1a5d98207  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part7/case.json
c8b2fb4d20c7a1250e2a7ac08052a1787863f10392c7483761e37a79814f817b  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part7/input.txt
e0cd5fe59842d0960f2e23669b6cfaace795abae0923d9f85342459d7b67e6eb  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part12/case.json
c1b8ff0d8c8f22ce3fafe19e345cdb39f890ba49d3d5d6ce60d41361bf458f2e  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part12/input.txt
fc427b3d5fd5b6dc3da3f187eaf86accf81a81823a31103e7d77af9fb469227e  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part13/case.json
3ae4f8a15012be447992661abe22811f491ae960f49a4cbebab871662e79bb62  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part13/input.txt
1c1ad505231dfe3fa5538f8e9fea7bc3b1bd11bc9063c8e294e604df0c030ba1  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part14/case.json
8f3934a5058142e578111bd2ec17c336900ffdb84825bb16f793b7599400a63f  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part14/input.txt
7a733f215b57fa097dfef5bd357d4bf5c222333af233a320ac0f95de0e8f4a31  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part15/case.json
0ffdd4c9e01f3114dce3da3cb1b75a9d003563a73f105bd5f9e2c2b4cb2bb11f  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part15/input.txt
92adb9f38bd96668c16504fbf7dc344ccaa19ddb05f81a98edbca9016d0d7213  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part0/case.json
4fec2628b39cc38490170f02b178c721904958a692d7398ea14d8355ce20a0f5  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part0/input.txt
2d5411857198f7e56759a9dbdd937bc6fae9d4f738838e3dca2ec76d90089454  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part1/case.json
86f5b6fdfb42381fb52bb9d6fb77449d5a31f2f27a35119b83127301ad36d734  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part1/input.txt
8295591a58a2b97157c6a7f7e57f93c60c6df5fde7177786a1a31f3d80b3c6ef  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part2/case.json
56924c7ac145bab6ab5eedcf7069ba424f63121e839e6a8d2297ce086651d517  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part2/input.txt
13778fdbeaa1d5e6561036b0264fb8f7e9dfd555529b1f6a380497aa10d29df6  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part3/case.json
8ea7551e3e9671358228aceef9f6a689ad7f3891fbd9c61f6d7bec67f65d5812  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part3/input.txt
ec898d5de552402ebd5b72ba8f3f35b8f940dd8789c6e80c8d017df1cf8764d7  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-timing-smoke/case.json
a3135dd305cc0bbb45ee3ea9798c3ad3b25dc47bfeb8894251341e7ce4ddec9c  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-timing-smoke/input.txt
79dddeeba7d850d126605408c29b4377f04ffde2a4ef07ce240177bd91554ff8  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part10/case.json
c2fe5a7f39d99789984ad814394ae1e22b975cd284f4ccdc99fe5f3cada8382b  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part10/input.txt
77729bda3bbafc47426e5f6438dd787d12990d95227ff9cc6d0c082eeb200324  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part11/case.json
e8902a4f38c21babf64fb3facfacdef7c39ef820e6c68d3809ee3d96f5cd4f6a  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part11/input.txt
9a22747713c613bad7789c5d9a79410893fa8ff59d09494f1869aa2b55aaac34  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part8/case.json
bb07276391aa51fedcdc6ce5203d4359cd5bd091fce938cc4937b46d02c934dc  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part8/input.txt
bc0ce4406c3feca53725218fbb16955adcec0900ff665bf9a7e5928b83edc4af  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part9/case.json
c0de13d48b54de58baa278b7c1fb6aefe8e6ef3c9670185545b7abcea7dec201  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part9/input.txt
9e0a3c4be4a24868ba33f0dc75b9ca68a613e33b10ea2a70039cdfe1f7b3f084  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-timing-smoke/case.json
e54e544be6562257ea6acaf0cefb61bde0a2ba9970a2eee5289b3fc2428a67cb  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-timing-smoke/input.txt
3fe6bd332baa2c9565d6983830eef7f36510840e5cf1d622ee4dc333fe140e03  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part4/case.json
48e0e76b9c51d72b6a946b1d1f90b2f499c896e1067108a00fc5c5b40d4ce83c  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part4/input.txt
88e4c8fc89e21739e266aa272477cdd04eaefdc7de90e4846da30591a867b81c  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part5/case.json
3d0cec3801dee8e2449695c13e2644c2ffa07ea24d4b910227b642dae232b47a  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part5/input.txt
454fd2d8480ee20c0c2ba54ad631035373008d0fec8f127af3a46ff16de8712e  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part6/case.json
15cf79945b42ce82d2cf0dcaf801d7bff7703bddd6ddc30c23f57baee9aa5199  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part6/input.txt
0a13f5c3bfd3583f6dcc35ffe19fff1a5e24ad50e30354bd771ab04c8f1fc467  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part7/case.json
d2e25dd1caf7873ffefd0b2518c7230585bc10077e9a7c80ac3f63b45f1a51fd  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part7/input.txt
87b01fa63601b547ad839337b0823b4c419826244030847d2fd743cbab6f99f2  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-timing-smoke/case.json
627556ce583cbf3655b786f317bda3326b58388d49f085d653adb807a4de22d2  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-timing-smoke/input.txt
39c4b2134a1cc0c40d09bbb2175a9f95a85844bf07b318ef096d24d65ee74533  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part12/case.json
0224f5e5577501d1b1378912ae83e605da9d7da3930b1894520c839c6eaae13c  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part12/input.txt
6905e1399c87b3a9ed8e4a0640b2e51bb7b0b20b12c4207d26f6b62517b75f65  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part13/case.json
bc9993746ee6a780422f4570afe189b3f09e40c0bd7ca8705d9cbffaf3344230  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part13/input.txt
544bb3037ebf6c9e094b3d2eb2dfa7ba0cab0544feee9759e997d503f6fe643e  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part14/case.json
9b7a6ce6450e98b79fdd948468a6e114c468d53456c2dff2f3144f3ac3f6de4c  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part14/input.txt
9b881e9f4072e3bd66620ff9374fc74da09e5b5e2ef0a20eb4c8d343d6b1b291  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part15/case.json
af94ad4e0d9ce2e90d629d452c936c164b25324f4b11ea4bdc628669aed3c7ed  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part15/input.txt
127f4db956c58bd6f078e1105d59e9cfae1afa1500a7e859107c074a37b08a62  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-timing-smoke/case.json
09bd3a0143ba1ffd76f764f113f2aaa58e322511701e92967240c472c29bbc41  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-timing-smoke/input.txt
083a1b9f64f5220823d6698377dc45636718f81568fdc0d62d8b4e7d480cdbe8  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part0/case.json
e23af6b48664aa8c5ff36d8a5bf181746c6e1598c5e6b6e5d3eb0a51e38b79b6  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part0/input.txt
7b997ff97d90f6f6b0e8a2c2cd07c7db7f91580381ed0e9bd978c74e64cf854d  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part1/case.json
040e7f1d7ba312a116ddbad7af2dbf686d1cdb12d61899ad0339a15d031d6b89  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part1/input.txt
5df63253f8ce43b52660cf43ae27513402c814746842db868930275fe6b5f8e0  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part2/case.json
68afa8fc75b5e89fa0770594df56f93204c419d16904d8c2faabf7844010efc4  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part2/input.txt
e65672e9004586c2e816027e911ae19d8686a7b27b9c5789e8a86783491c937b  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part3/case.json
6372c3fad8be6f2818f560ff9d750839aaa18c3c70fb2658884a998eeba85a90  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part3/input.txt
b16d8ff417c4ce79af00ac7e29d30f336c1c77505b923f619a06a01fdb59d337  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part10/case.json
7baaf4aafd10c0f59ff8d099646f3aeca545af2df5bbf6c465f4265c3b1e6e3e  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part10/input.txt
20e46f14ac94c9b80b9a973bc30cf090467cebba53b54fad7068f45c5c553d94  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part11/case.json
221a2d1de324dcc4c521984f03b42910e6f32e1e863b3dbc2a55b3fdfe568699  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part11/input.txt
f6b526ca7ac989b480b7876b54adf99c420af1ab683dce646b4aa0ca4cd47ae9  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part8/case.json
34336fc74f1066ba4c0dbb8d10878babda5ff6c4aba6623941bbb6e800dd5a6e  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part8/input.txt
d2990f59bd43b5b0310f8eaf3561af41047dd153a16e54679f98a49ebca0584e  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part9/case.json
05128714c972be33ef845f53dfa17a1eff94962b4d7aab1bb8ea4950cdc2acbe  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part9/input.txt
61d45b0ce5d6a4739ce74d3ddc7fcb5832143d6b562fff497c5db25d3378bfe0  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part4/case.json
11d0053a4f0edf6de9b94f52fe265e96941a07af7c9a8c7ff79f429191fa8828  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part4/input.txt
586728ee86844bfa3de2e63ae2621a02aa5e63acf41362c4eada69839ff9fd5e  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part5/case.json
47b4a65752ba972ed951ffeeb44c340c6c767a8a8a82a34c7a679d5c3631d863  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part5/input.txt
4e139de318b54f8635ddb8b703183080fd7aa1f1de1238e35e6326dd54fd77f4  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part6/case.json
f25d066fa8755af38828568d1be9c8c5aa6cc47bf5faf6d5b5ee8b3c1be294e4  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part6/input.txt
0cbbf9f87322165e07e1f98510c16e81eabfba730ac8d26e4ae1d0eb22f2c43f  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part7/case.json
64c9c41477cca6f4aa98e26131356cc830a16c759716a5466ea7c51d9422f9f4  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part7/input.txt
62db4e7415ab2649888aaae7a5fa097b56258d53932a233b1cb84fd0a494da1a  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part12/case.json
99adba9b7ea26d8819ccd7def7cdf7646aa518133e6757a8d027d51092a8da32  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part12/input.txt
e0d4498d11f02db8c2713990766d727702db5f68f5f550072617ba2ee71098a8  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part13/case.json
9aec4f529f7ab658e555dc11d48173f065ab10e47400581aade735ddada07178  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part13/input.txt
561d79247cbec3b3bad0df8ebba248db396f6a4aac04a4298a1f6d4eb2709c92  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part14/case.json
7ae6d449b5fbc3adc2e05291a3ed1b4f539ae5a54aab1d46c6f3ec956e2dc4c1  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part14/input.txt
8bc414acf2899e1c1fc1e44b55edd5f3223766b87e24633ad749a5c16cc05a8c  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part15/case.json
0a5bda33f22ab0bcc99d4e73a3da75768b235fe7538d9ea4e6b92aa6f7c506d2  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part15/input.txt
1c3031316cfb3ba292591115d4d18366cc86cd5279eb11f301baa49a485ee0a6  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part0/case.json
6965f8c54e7ea579fe13161226272aaf7ade2186f2a3a6773b09b5c35ed0cd9e  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part0/input.txt
ca40ddec322ac51ee333d6ff74e74035dac0a1e68e7b35b6aee0aa3eb874a4e5  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part1/case.json
e0a2b94b95f78bc9b3bf756d82c91c8228f2b8261bdbc8ab1ec601f9f1b711a2  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part1/input.txt
f240d6b20b307a9b9b5b4acfa026a2d80b0b56a96f32bf52cd4b41e92175bb55  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part2/case.json
b71d226becfc7b44aa06a846d000711630e7f76fe16d0ae265fdfdc979904589  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part2/input.txt
bf5f6d989bc489424187b7370a0ffe0025fcdf8be825cfc68fb12873f8a31380  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part3/case.json
2464bcfa7536bcfdd66d4d8c25f5b39d78a9a8b6842dd57ba671c0bc1798c1f2  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part3/input.txt
4da71f38a7310a71aab43905f81ae8e821e6fdc629c81ad3d242094d0886413b  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part10/case.json
e1af23e91ed161a6f9226954fe4bb19b9dc8c1472e42afc067a5d671f6197a57  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part10/input.txt
91741227b97639c2d05cbfded4d2357f56a9e7fb819b108628ac1ab8905998aa  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part11/case.json
a31df7ca0d88f541de1a41e7bae49b4035b449fec427ddffbbf148190b4e663a  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part11/input.txt
13ed25930844658ec8286dd219fdd56c65af48bc301593fcb067a97e4f7dcbb7  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part8/case.json
e4a60b95cb76a701f3883073f1002796a65a9ab1bb3cfde8e2a4be10b65615a4  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part8/input.txt
4f1f6ebb45614fb77b7fe25a8b9f1d4564f403af346f15d63c657219cd567083  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part9/case.json
1f9977c1dcebf666330076ccb6dc0b456b10032fc3a81b254d51f0619ec6b717  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part9/input.txt
a9085b2ae68636a91dc10bd8463c2be4b19c74765be35e8ea6f3db5451480207  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part4/case.json
0ac301998308e436c9b1d22d051c5c1732a4c2adf812e9aaf71bf8b15aaf578f  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part4/input.txt
42914eb9c800dc4ac54b02b5f052cb38a5f7b8afd1e1748677497dc4c7435713  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part5/case.json
dead92255425f4fcfbb8d89009b7ae04c208469fd0d7cb256aa55e492728da0d  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part5/input.txt
1b2406a4cb7080a318af64c74f614e350b48efdc1a0810b25af02160e439e86f  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part6/case.json
1025d20712d71d4a7911f4a352c44048856f6b71673fca20c217dbc6a6e7daad  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part6/input.txt
734ae6f02236c96aae4385ab593e7c7ac17731e79bcb93250f272e6c312d89f3  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part7/case.json
e57f47c980ba5528e4fc4524c7c69f55d768e9a655f1a9e43d946c8c33c572d1  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part7/input.txt
960e69d4c4a329b89f35c620af8da0f5595546d3e081067946a38231c57c7e69  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part12/case.json
b2bf2ab0c1b189561564a3d16fee80e05e34a27d2bbacd9660cb5e1d236a0481  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part12/input.txt
8b1709cd95f0d14ddfdb93eb4a1e09b821bf6c6283ef59122af3b0e0b17e80c5  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part13/case.json
ce2c7b28fdcdd824c99c02f3e164794c1fcaabd63eb029e3119f5a722ec61c0b  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part13/input.txt
24ce96a3f4991dc49f7c05bad271ba3df132c7af21f560709eced3c17563dc08  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part14/case.json
fd6183a25d396ed9bbe646520fd0a2e690c6f400a7b43da31fee6abad6ee7674  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part14/input.txt
9b4f9b7cfdfdfb522fa7502e48a55899d0a79353f52d77a84088d504c4e7ae77  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part15/case.json
52bbf738d7713222ae3715148277c8312dce9a4926ed7e298b41f69c11022f3b  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part15/input.txt
343218486a7c086fe01031f29411aeab40fd6e074c9785d731798479d117d7e4  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part0/case.json
29a73f36400c34042e8dc209e68e13e735dab49999398722f5db58be2419f163  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part0/input.txt
474f51ae87da6bcddcdfe16ad3deb4fefce4bb7a99c15d526d57e9625688f476  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part1/case.json
ce1d39a0583a44ea8e66ec096b75e6c3a1e3fc9d9b9ad783a31bcaec08ae9dcc  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part1/input.txt
bf662a6d5e7ba07403251f1a325c8f11e62cd35f0bdcbe94c659aac7adf4b2d3  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part2/case.json
beb73019728522b21bf1d3d105a0bb509fbe9a28ba3c8206af8c3905b61b12c3  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part2/input.txt
6d285d28d1828ae3b24cbc5430f39b6100733147e93fe26c3cd9e01c34f8cdca  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part3/case.json
76dc2d3ba916a762a44ae14cea41b9f32ed5be7673298a4d020009f545101fcc  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part3/input.txt
f5688d86418b42b707ced92fd14d4b0b03fa2dfbcd2069bee4b0f060d2697b0c  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part10/case.json
b893f5d4cd4e1ba84c945565f171ba63d1b3cd06ef882d5bd7f49199c77be009  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part10/input.txt
24db7f2fdbff071741aa689cc1bcbe56b32e9c04ca8de6cd201e6a32db70c617  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part11/case.json
6fe3c234d7b3bb7dd8690e23e4bce430d4fbc3f6311ed7e476bcecbf2946e3b2  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part11/input.txt
3bbcbacebb2fba1198646380dd784405c409325f307a695d727ee051097d5196  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part8/case.json
951c21807eedf1a5071f17f2b82ada599b4d66da772fdf81b77aea4f0b817da9  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part8/input.txt
bba05697f8657131f1f5c74d24b49ec4f58fa11118002181e1f3706e46f11f26  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part9/case.json
66d93fd802c9bb0a9f3e0e7059e4402f8de4d447a20e157b2681bb34e8b2e8e6  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part9/input.txt
274295a5ba8e2c6014e91194def95851c7aaebc5042649cb65b4d530d9dfb8cc  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part4/case.json
c61ee72fd26c7235d22e6cf93a019f83adfe8a96b19d06cbe54287a09466508b  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part4/input.txt
d8609358f7c4b8bba2ba38ce8423dc8cefec4344842de93d9e154e9cfe73eabb  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part5/case.json
768223c09509c890f0256842e2a4fedf192038d0d0f8f37a4b5cb5e71bbc82c0  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part5/input.txt
474fb7dc51ffc4d6a103bf4052f958e1e1a6fe53028566284655db48792810a9  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part6/case.json
2bf4301045719145d804180f733eb619325e20ffb86dcb8c6f61956c8151f74e  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part6/input.txt
01e4fa14b97d1176e9ee8f441ca1259f8a399a1b3e41725969144761bda16305  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part7/case.json
cb165b53184f0f77b20c4745c65508406f184fea16ca1c1c4af7a183a40566b0  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part7/input.txt
c053e4f2d0381ab866e0720b348dd732c71edb0714e1eb1466a586a88f240b48  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part12/case.json
774f336aa4c5b6d4a8126b8b9673a40aac4dc2883f12b56e9f49f59de7196bb5  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part12/input.txt
e509bd15a06dabce41994a64d7efc8a66d02969ccd6e72669a086c40b3a81666  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part13/case.json
0388682e6f73a9f01945213a410f6aaeb5ed53883394ac8a0496492a0f8940be  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part13/input.txt
5515879912517e52f8accfbb8037936ee482a195ee9801be88c0e5ff1c484916  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part14/case.json
43d08d2524149433af5a5048d0b21411c074c287adba8327bba7bce2981c7243  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part14/input.txt
aa74ed8a3cf01615272046f2d5f1e52c997de44f35645f15820a300bff2be33b  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part15/case.json
a57e79fadab5cfb3576bb1faa724439c15f8050c949c7a62f981cd7989ff0ae8  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part15/input.txt
6768d4e0590b3e0d14be84261e423a3878e3dc0f740c79dd5a2170e0bb182445  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part0/case.json
c785e0702e8700fd2e156ca1f2a40b9d4b48bc2428d0a09246e1ae2dbedb662a  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part0/input.txt
aad74c0ad6c21fab40a7faf53e052c8163dae7a490346927353629531aea70aa  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part1/case.json
d10fc9b24fc47dce8914490c9e5848dfd0146fc67c983247c1b3cd4db8d57841  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part1/input.txt
4088dacfd4bc728e4979267cf785555b8fe244edff04f24eef551ae5242df5eb  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part2/case.json
d74c95a0895eaba02c1c12e4b87e1be21832199a178d2486d0a2e4580c69c898  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part2/input.txt
b11821d31776a41df4147490f30f69b6325f19195c8053026e5b0719751dc9ab  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part3/case.json
a89f0e3ef0181eae9e59d3bd32c9a7adc183d680b2ce690f79f61af5054c65f3  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part3/input.txt
063249503278226ff8d9bf9e2f4400c31ebe2828f583a3a107c33d2a2fbdefa3  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part10/case.json
feadc045c123b263d4d8748516a6c0738d45d6c34f09b800fd2eae95023b96a1  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part10/input.txt
197829acf4984edcdcf4526896980c709aa73de3c5649d2f36c27de9660ddee9  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part11/case.json
c654a691c912443071923e50d8ba08378a38f75fc03c566626096bcb020cc3a2  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part11/input.txt
c8cd4a29768cd66e766f23ce44f0997a9adf0a83950233da347aedc4a5deffa7  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part8/case.json
b7c257e474483f113485857ac0bd02e6240c4f92250ce858139b9a25c6e7dd18  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part8/input.txt
94e06587c4628811ee0202816be30d9c9799374d8992ad7c4284a5b49950d9d9  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part9/case.json
b08aa10f4712570c4ddb3dd49ba2060d6b39287d3ba8165b6d1e3e7d87de6afc  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part9/input.txt
b6d7b4aa3070a34aa7b2b388711238e3ce3622fbaeff72ed9351a4159bc986b4  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part4/case.json
cc3ccef356d76eb76e3c962f96b9fb5c3116d00107e94442351e5c47f8d6d723  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part4/input.txt
681966a05a222bd2c81c9447d8f8dd40f73744b8660f052ead134f01402dd4ad  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part5/case.json
debeac1fe8d08484c6a28672b683c1f0f502a68a4a507c5e90e61d967319ec44  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part5/input.txt
3cdd490e5529beabeddc3df4a19c426b8ef8b88bd9257c94899dcd2d8145cd12  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part6/case.json
456eeaa023c2754ac78035c06c99ea2feae274e48c460a55eb7ff757823636f3  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part6/input.txt
42672b34c49f000b2d4c18f1ef8d8bccff017918d622af0ab46a028e751ab5b9  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part7/case.json
31b969cab317d04b0120cc104ee3d729b3ba3acd72b00745fd0b196654e34ab8  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part7/input.txt
46633778435a5d502c1daca6936ff6d24e610ea3ac32c51474d0995c328bc819  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part12/case.json
58c6be816a8487da597a6e39fb8774a31b90f4baf32a06378eb9f96dd979d852  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part12/input.txt
5a746275cc40cb78c8cecdfaa6b881b7f870afad17dee19f45ee5e5f8a232920  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part13/case.json
2ac2d94d7b0928e59fa4d923b555e7a819f9c98c72d55406c4d7d36ac278943a  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part13/input.txt
3bbfa8b6f3f77294c21853d99f62494591cde9e1da304a139e978d8fa56b00d8  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part14/case.json
fd426f60353812261ce40ae914d1e2051b6c733cc08ddd41d765f416538b20ca  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part14/input.txt
96bab17658bcb5e71b660f17b17ad2c2619e82bfd1f98f7a8e6efa5e0e9de46b  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part15/case.json
74f6aea9db6a419a1ab59007d5d9f2c0c24af17abfa181085791a9b61144ebb4  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part15/input.txt
a1d4655a0aac8ec99a71e83f68b1415856532be166cdd8d91f20da0d774cf6bd  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part0/case.json
79af9b20a7c911eef0a48e6394cad6a9ced5d7a5a9352d26b6f377e5164371eb  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part0/input.txt
6f676ee12a45640278f2b32d0aa4a9ebcfda3de58a0eef4b723e1789b979cc6a  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part1/case.json
0c8c3a64c3f9deb96c6e404d9cb3bbd2601394c06ba6ba58f936869ae9c88cc9  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part1/input.txt
081840fbaa91b1e39cf21909d502fe02ade812172725b989a5e1f4926ef4d8b7  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part2/case.json
1493b90e2acd54d0108d4cd0afc5d9c6ca5379dac03274f143e69233cf12d0a8  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part2/input.txt
3d8a8a17492ef44ee5b7621867978512e451c09a164432d7c93f5c99cbfb3983  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part3/case.json
f617f28e14cd625196b7c4d27d68adc14e4731ccd45e6562db746c011bd3891e  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part3/input.txt
1dc45b8271737544010e69af341646bb22b8df58aabc8bfb1293442c19b72860  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-timing-smoke/case.json
99e874b301d296d390e828581213f500b8be61d08d8f30102876aefe4e40fd2b  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-timing-smoke/input.txt
50c4c7b2a05fc5c32e0122d4a32cce7528abdc3c70a361c99e9e6c1ab8fddd9e  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part10/case.json
52e7bf1a5ba389f7a830828d995b0722859ecb2ccd614f1aab258d3e976de789  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part10/input.txt
c5831bb9708c659d1447468c664378ce3f4f7bfa3fe99283604dd8e6f6b584e1  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part11/case.json
1233d80725147c49a5381090d1bbb01197a9b48494b8432dfbc289db29e024ac  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part11/input.txt
5006a3af39b6a8e7dc4d0cc55bcbc4bf56f3b3839cdf01561b5dff35bd2cdc79  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part8/case.json
e6433b95072be6428a9f7e5d0c7460c93adff2a6f43e534d1a145402b602fb10  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part8/input.txt
66bb11cd4ebe2df13f1ad44a32b4ff24d3719bb3de49c496b6aaba69c74148e5  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part9/case.json
0c3fb51e8dc280920bc4af39713a3705f5d9b013fbecb5b5dc802ceed43d8007  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part9/input.txt
498747f20ceb861496849755436f6635fc876f7b08f87d29a508937a7c003988  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-timing-smoke/case.json
05a521e18813326b6d671eab0d12e35b9d6f70a599962920ad6a0a26659c8523  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-timing-smoke/input.txt
b6723eece668406022f671b4bcf1eb77646ccf47bfd5cb667edaa659dcf38ebd  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part4/case.json
b3b29c17517892940e928bccca0a5638b3e86b07bd77e69317a298faf96b4bc8  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part4/input.txt
0466de97357aba8b671745dbb541a71a301d59bf7970145c4c0ad066ae8cadb6  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part5/case.json
bf43bed492bd6be6ab0bf519c503503579e47dbafbc695ce4759a04db9da8ee9  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part5/input.txt
b220748eaa26a7a96d99602417ef79510e984d9d9f472addbed162a8ca1a3f88  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part6/case.json
b63656311e6fb7b98bddaebf8300583b84c4dffbd878cdd83154882809279394  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part6/input.txt
0da7b1f42cbe5120a643c975929ccb661ce91e033ce09024b013da234e1c1573  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part7/case.json
92381bbf13b2f4f904ea93245d6002e1d17fccc43ffb84d9ff80e7a18b155479  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part7/input.txt
90fd8557849085efa4db3ebad8e40c88d034fca31096ad6105572267649038ea  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-timing-smoke/case.json
1be32efd866a1c8f886a4389295c190ab13667bc4310e37815485e33c48d6006  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-timing-smoke/input.txt
941e3f46bd855c23dbe2a6002929c913f557a3ad86ae08422a8f309aa9297258  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part12/case.json
aeda25a264763054b34228efe950ecd664b2aa54dd8ad50496fc41b80697c00c  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part12/input.txt
ca657d73e9c73afcf7402c9073179b6c840851d5f491fb67e518479170ddd3ec  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part13/case.json
049facb62b200d6ceaa3641a75371a0d3de51b5e4c52b30d2b084a6bdb990f16  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part13/input.txt
77b4ae2d9b8e509ceae0850d1c831b30eaf3ce951b3ca19ee1dffe5ba2b621f2  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part14/case.json
cb6f5ef65f3c454240722ac5e25150c7ae945145df9fb1f72bbf2d4f93738eb9  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part14/input.txt
2d68e025e788bf6e683eb220ea3dbf68c91706575a6fb2390a88b5fa5b8441e1  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part15/case.json
e1e7902363c3c5735fec1ebe6adfe75a7e982a158c7791a3fb5e8f5c96769b39  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part15/input.txt
2a9a10221447d70c4a0f2dd0298421c1bc1a181300cd99005319ca1d7fb8e232  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-timing-smoke/case.json
1793e74862e37bb6d89b8cd2733cb61e5d4fd8a8621a5e62e40c4e6a7283b610  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-timing-smoke/input.txt
```

**Dated note (2026-11-08, ionmc sample chunk plan; orchestrator, before any sealed sample).** Timing smoke of the
table-accelerated exact root (RUN-20261108T043826Z-968614a1 at c5838411, 1e4 histories per case, counters all 0,
`mcs_table_bracket_miss` 0 in every case): per-case ms/history and the resulting plan for the twelve 1e6-history
samples, split into chunks so that every host call stays below 1500 s (plan section 7: chunks are independent runs
with `k = 1 .. n`, seed `20581004 + 100 r_index + k`, histories `1e6 / n`, combined per batch by the evaluator):

| case | ms/history | hours at 1e6 | chunks n | --scale | histories per chunk | s per chunk |
|---|---|---|---|---|---|---|
| 1 | 4.52 | 1.26 | 4 | 1/4 | 250000 | 1130 |
| 2 | 7.03 | 1.95 | 5 | 1/5 | 200000 | 1406 |
| 3 | 12.35 | 3.43 | 10 | 1/10 | 100000 | 1235 |
| 4 | 4.54 | 1.26 | 4 | 1/4 | 250000 | 1135 |
| 5 | 9.15 | 2.54 | 8 | 1/8 | 125000 | 1144 |
| 6 | 15.18 | 4.22 | 20 | 1/20 | 50000 | 759 |
| 7 | 4.98 | 1.38 | 4 | 1/4 | 250000 | 1245 |
| 8 | 10.46 | 2.91 | 8 | 1/8 | 125000 | 1308 |
| 9 | 12.83 | 3.56 | 10 | 1/10 | 100000 | 1283 |
| 10 | 5.04 | 1.40 | 4 | 1/4 | 250000 | 1260 |
| 11 | 5.26 | 1.46 | 4 | 1/4 | 250000 | 1315 |
| 12 | 5.96 | 1.66 | 5 | 1/5 | 200000 | 1192 |

86 host calls, 27.0 h of single-process compute (root solve: 42.2 h). Every chunk's `npz_sha256`, its JSON step document
and the signed host-runner record (run id, output manifest) are the sealed sample artefacts; the npz arrays themselves
stay in the published run directories and the host-run store (not in Git). The chunk runs start at this commit's SHA or
a descendant; the sample seal commit (plan section "Order") follows the last chunk.

**Dated note (2026-11-09, ionmc sample chunks complete; orchestrator, before any TOPAS run).** All 86 chunks of the
2026-11-08 chunk plan ran at the frozen implementation SHA 1831f2149506b0bd403259f2d00981b835837e4d (cases 1–12,
1e6 histories each, 26.1 h of single-process compute, every step `pass`, every diagnostic counter 0 including
`mcs_table_bracket_miss`). The sealed record is `validation/results/transport/v8r2-samples/sample-seal.json`
(schema `ionmc-v8r2-sample-seal-1`; 76 host runs, case 6 two chunks per call; structure checked by
`tests/ionmc/test_v3_008_sample_seal.py`). Two non-sample runs are preserved in the record: the superseded first
chunk at 20021d2a (uncompressed npz) and the first attempt of case 12 k 5 (RUN-20261109T074122Z-9838df64), whose
step passed but whose signed record reads `succeeded: false` because the orchestrator edited the task worktree while
the run was in flight (`dirty_after`); the clean-tree re-run (RUN-20261109T080127Z-1011e276) reproduced the identical
npz sha256. The commit holding the record is the **sample seal commit**; its SHA and review id are recorded in
`v3-008-acceptance.md` by the dated seal note once its review has passed, before any TOPAS run.

**Dated note (2026-11-09, seal recorded; orchestrator).** Sample seal commit `0c68b78ac191c850bf9101dbacf8eca12f4e9b1b` (the record commit 84b8b2cd plus two review fixes: D1 mapping-domain rule, MCS-off gating; full-coverage results bit-identical), Codex REVIEW-9691e1b098e940d7bdfd1063c09325ee passed (finished 2026-11-09T09:15:57Z), VAL-20261109-092143-547E49, CI run 42592861874 success; dated seal note in `v3-008-acceptance.md` (Revision 2 section). No TOPAS run of this set has been started.

## Appendix A (2026-11-09, after R2-SIZE-1 regeneration). Bundle hashes of the changed and new files

Regenerated by `validation/scripts/ion_em/gen_r2_bundles.py` under amendment R2-SIZE-1 (orchestrator decision
2026-11-09, `v3-008-acceptance.md` "Amendments"): 209 bundles, 418 files. The listing above (205 bundles, digest
`a8a1c27d…`) is kept as the record at preparation; it is superseded for the files listed below. Every file not listed
here is byte-identical to its entry above: all IDD-EM and LAT-EM bundles, the He IDD-FULL bundles and the 13 smokes of
the preparation set (the three 0.05 mm LAT-FULL timing smokes included). Changed: the 48 LAT-FULL evidence bundles
(`d:Ph/Default/CutForElectron = 1 km` with its comment, `case.json` "physics", "cuts" and "rationale": "seven-module
TOPAS configuration, delta electrons deposited locally"; `NumberOfHistoriesInRun` 5000 unchanged pending the LAT-FULL
re-sizing) and the 32 carbon and oxygen IDD-FULL evidence bundles (`NumberOfHistoriesInRun` 3333 and 2529; cut
unchanged at 0.05 mm). New: `he-water-250mevu-r2-lat-full-ecut1km-timing-smoke` (20600424),
`carbon-water-430mevu-r2-lat-full-ecut1km-timing-smoke` (20600425), `oxygen-water-430mevu-r2-lat-full-ecut1km-timing-smoke`
(20600426), each 200 histories, and `carbon-water-290mevu-r2-lat-em-steps-verbose-smoke` (20600427, 3 histories,
`TrackingVerbosity` 1). Seeds, names and part numbers of all existing bundles are unchanged.

Bundle-set digest of the full 418-file set (same definition as above; `gen_r2_bundles.py --check`):
`e51d3b3d1bf00e51daadf889f6541bf8bd6ba8aea6453eee132110e07ff0c4d6` (418 files). Changed and new files (168, sorted by
path):

```
69b372ca81649397bacaa2eea1d454590c2f54db5a78fd5ae0bec02723b384bd  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part10/case.json
30134e926b95eabfd62929602294052860acc5209af01949046a16209d1c57c7  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part10/input.txt
87d61e149162cd8d82b96dd8069f72c9d04540ca9f25100e6ee1304b5db8bc75  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part11/case.json
303a09dea3adab91234893938286e31f0e27b8884a9ccb62a58f3a292c55e887  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part11/input.txt
14a2493d681ab431ae3abbfda5e8b3022c2d3e68c28d420efe440105a64324b7  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part8/case.json
93f22c5e7080e86bf54731cdd1a01ebed20414935038f8ff5572486aabc0b27f  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part8/input.txt
7dbe8f0c03393ee02feca35ad22a1e34a75a0b93bb94caabd7ec79ee982aeceb  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part9/case.json
a6d5be7388abf75e99362dd2d69fb624de3e9f277ac0f2a67de309db7f87e466  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part9/input.txt
35f95a048545ee22685ea9dc0dc0e73669ca6745cfd5cb3f0688d05bbce93714  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part12/case.json
552b6ba828f7ca93251b4217089a456cee0c738758b9413829437cbe418298dd  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part12/input.txt
d367f8c3ffa4d04c937fa4f77bd265259b1cdaad1b4ebe76970cb0810c736bbe  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part13/case.json
e8cf3bb2391eb2b586a0de11e8ea3236da08f7a252ce71ea80172b6745bb8388  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part13/input.txt
b94a65733e495ae40f1636d1ec08857d90d988b64a97b0eed84f5318910743f5  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part14/case.json
7464bc7af8a947c909202f8b0675915cacfb7b0b702c4afdf62848373b6d3b70  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part14/input.txt
dc2d28158fd1234a8d3657bf783fbb8b6ee16f92bd450d3e7e3bbc25453bbe4b  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part15/case.json
aa70add26790d6a920f84c02113c297ab187dabf169483a8dd3f59682a39efa7  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part15/input.txt
9fe4bc4f8116d08879879e31d1206b06d76f5afbf2c9f10d817f935f1df31ff0  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part10/case.json
73371d0288f2686fa9030c74c46d07b36a1b34f82822fd83c4946684e0aee2d1  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part10/input.txt
35621ca3f3572d5c8ff7915771643db86828f968122ff4da0196f45f92b19b3c  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part11/case.json
ddc21343942e8c37149fb9b54adada1d6521c48fdcd01209287eccbc71fa1d75  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part11/input.txt
4a036dd7417c78254fa77d46a40111e94061058e795b3ceea66d7554c37e746e  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part8/case.json
eef8c91e7be1583094ea57f0cc0fc176222b4967ad765daafc1e694b6e1f6939  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part8/input.txt
6ef9b550217264e8a91f19dfc12f415a180c56627927b78b05fb15cda583f046  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part9/case.json
a19b2ad89f0cc2ac4e2bd7cda968cdd8ed32cb2c65d7f7bf00aeff202cfcc177  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part9/input.txt
812c1335cd298ae102af825b80d781409b83d0dc612fd48c74b83637ec74aead  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-steps-verbose-smoke/case.json
4e9e5179cf0b3404a7174ca44ee9efd97460f0c927537a49ad89150efbf7feec  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-steps-verbose-smoke/input.txt
b60a3214703d64f780abaa1ce03e461d1b0cc4e819dbcfc08ab02155de8291c8  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part12/case.json
c77fccff3bc7d06ad39da616f91231c14320e06bea2316ca07696928dc40b097  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part12/input.txt
92200c9e0403c457a8be92f77da3131adbde598d61dbc116e2e6a139aadd0af4  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part13/case.json
dcf1b60544f0f40e3c7e58b056f1358c6c68f3c198d9b30c77fdf53b630cb858  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part13/input.txt
e9194f420474e12a9ebf7916b291586b18be86e87905e3bd4469f36d25f93d86  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part14/case.json
07b5772c76f87e6eb2fa9406957d15f64b2f991f5c737946b08a2d7b38eeaeff  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part14/input.txt
b4f94ebcb21f33ec322e80432745456a0d000cf07deca25e4ccaed683e625842  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part15/case.json
7be4fa8fece28d2cb273c3ddb0e0ef082e25112f75138c61f61ba5d1bd9d5ee8  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part15/input.txt
e2aee57bb7ca86dfda93ad4fb18eb6e9fe613a1ec41ac8cac19fd3bad7d7c330  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part10/case.json
0aacc7bbfae3a7ee6dbd9e8977aa0c62a5a537fa3ddf5c42c72fcbee8ff48626  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part10/input.txt
ba45181639e5b9e62ecca9ce8ba81a0d66a9ac51a703f63bc360e519243d1982  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part11/case.json
7c123b0036dcc55d5519d52a6ade9c80049bb8b2b491ff8ba61599176e628f35  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part11/input.txt
ac1c50c7375740955cb20cf19e49504d8f2b3870557950086f54cf48f7fc795a  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part8/case.json
346d9eaf9ca0679831bcb036e4b0d4675f040de25f11a43e7d350dc1615fdfc3  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part8/input.txt
a97aa605deefb3caed825d1ce55e0434d72bf5a6c37b52cc5f507aa0e80a6e21  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part9/case.json
f6decee8115cd47fb8bcc0242710ab4b78fd64bbce10f7e69f81bf99aac56d82  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part9/input.txt
94f255cd117727a4370bec056c4cef18944d670abcdca31ba8adc6794e0cdc16  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-ecut1km-timing-smoke/case.json
92c029746c421a6afc02bc159108818ae4a3b3679e59c043d05d2984ca393517  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-ecut1km-timing-smoke/input.txt
05ca3279c373a805e39625311618d943a92ce27b36610865b7b0610721bdaf61  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part12/case.json
f889ef2f3b0241a32b2fca2805aaf1c43318047f27ee60b38df348c32c5344b4  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part12/input.txt
c4cc85745daa08f58294413aa0fa27baf37b65e55be595c1323171f017abdd58  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part13/case.json
b71d73656884dd69334fbdd8931448c1527b34ef28ba27ba0e5d86dad7736531  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part13/input.txt
26674a35b36b3e1977bb9429f02f4fa5de1e3f32a459d15512104e2d51f0c52c  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part14/case.json
37ca7d9a62fb4694fa69fab1a8f2f93391aae9d6e4b930762414d6363549e475  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part14/input.txt
052828170907ac32de2ced3cc2cabcc0ba0bbb3cbf4edf15138198502f129673  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part15/case.json
32d432837bf4fc891ca0ba5cdc2f2c4602e4355b268e3528bdc04cabf1ef64a2  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part15/input.txt
102ac5714974fdf54727d3a9c83d8998d135230f9463e1e02ba0b8fe942e3a58  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part10/case.json
8845ab0afd37efcb721239d9cf4aa413980a96553f3526c75e49c0b7f9a94792  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part10/input.txt
1e38110abd91adcd0c0ea326c53b087944e2be0b186999f65dcf25bc49a62236  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part11/case.json
eb553d4c3abecc35ee42173a2f179f01cef91f900d30e1be24a9c4aa60d338b3  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part11/input.txt
e29506453ede430f65a3ae18d73d216acb04bc3520abdeb903fe44765295dc01  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part8/case.json
a09b92eac0adfa6a314566ccf7a4c0ee378e4224c81e6c01c82dfdc7f30c1d58  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part8/input.txt
852123dfd6f76e1611dae0284d7d8309d7f870ecdd2609a23cbb6050344a8e6c  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part9/case.json
f29a4bab41349e2a45bbd7b6e6d2b5af6d25703c2b32dfb3e439a42ad2c4c463  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part9/input.txt
c2138f993be883ad0bc27939df2c85d3bdcf0d8e0cd3226c2739a6307086b08b  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part12/case.json
6b447d344c9c4b53738d02064d2055c51bfa2096f1e6da54c87c0124f80723ff  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part12/input.txt
4dc9a0287eceedfce14dfc5478f62214803f9861d282293e684c4d4b358fea31  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part13/case.json
08e2bfb25521600f36a24ddf0b9594be823dcc46e956b216f1993c0e059d0802  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part13/input.txt
fc23363f19b735be9cb2ce2fd81c65e20feea261a302e637ef7c404a792d6a26  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part14/case.json
87132a396171964b96d2717d0dac47217a6c2400e6933412921102030acfff5f  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part14/input.txt
31cd3c4767c1202a4fce73c8d80826254404d198c24e95951a0c88ea495c7279  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part15/case.json
b467031099dc0b5f7cd88d71bfc22a2bfa30858c4d7db8939deb0c2701896705  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part15/input.txt
f8bdafbf349167631b47a64a83f1fcb8fd7656d175d6a0ca29a0c66b3be2b76f  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part12/case.json
b79a4c263bd87041adf8f1c58bdbf22711ef610887d6eca7384f1e57d825e2ea  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part12/input.txt
6e87f71d10b54f73919652c24971fef46bb26f03dd4d17a16ced53b9fb6f523a  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part13/case.json
fed37c409b1c2063ed77d696e3d121dafc912ed45ad73e55101330f19575e099  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part13/input.txt
4ea50c5f4748a795985f19dcf2c07d8a659b1250d2e176294314807f3910cddf  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part14/case.json
25dc070991f68d6bdcbafd93c57a05b3807137bc0fccea584bb29300c33482ca  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part14/input.txt
64ec2dbf18f1531aab91a93d48bc65d345975fdb6758a1dedda8edc1fc83f6d5  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part15/case.json
98d6d2bddaf612f86ef19a70dec7c672f20db5cc1a79ad96b332319803ca9836  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part15/input.txt
31448a00fb9a5698fcef3a8a77eb756fcc9edf1e3c39205c11745bbca4ba7369  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part12/case.json
3810571d861e88cef4a206b88b81e240a9b7741749706749e151db0b842e4cb1  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part12/input.txt
6cbafab0a8900352cb0dadeb097d90caf731f358ee9fd7dcf63419ad86a27c7d  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part13/case.json
1081d21f4dd964eaeffe72c2c638c7d415723164a71043be57b4f5dba970c843  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part13/input.txt
9d0b8c4bd5a6f698795febfc7729882e8a4b61e3a1ae3a3e4729521a19d33b66  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part14/case.json
cee8e63728bd5a288ee0cb69f6f837513044303cda20cbe499a7b3416f03e618  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part14/input.txt
f4b38340d1cfe7dad5c0436baeee6808aadedbe4373dbf1b70247bcbdc93e936  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part15/case.json
a662cf1ad08eaff127ccdc30af56ad78b7243b2f04cb2a5ff53d8cdf0e9a19ab  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part15/input.txt
6567ac23198f5b65a0405a6209200409daec0504c4907d40c203b0408369ddc6  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-ecut1km-timing-smoke/case.json
00428331f3a074b949db04fb054ddafd6c758e41aa4654b263a6148c39676483  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-ecut1km-timing-smoke/input.txt
09f0deac621215099c2a1edc55b848abae067672d454d1139d887d07bf8068f6  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part12/case.json
227800e52a0397dfdbd786dac793a4e03f9f9e65e1af07d791de34991f992d23  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part12/input.txt
0750965220cd177b57e5e7488b989314241c58995f56c17aa2b7524d5b19ae0b  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part13/case.json
b835f9fc706299515865394350dc0dd2b8be160e39c330fd7b0c4494d4b93f3e  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part13/input.txt
e9c06ebb688ce74def1786fd8c726c51fe4f7941adedbc98a95e39db7c8bd2b5  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part14/case.json
4e6f217518d4937518491ee8600fd065b09755a0549d4a4a6d829d8d3e193102  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part14/input.txt
12ff3cfaa2314475beb5173ed24a0d706adaa06c9a39f3758f3123f9636b2b60  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part15/case.json
49f864dcc508aebd0890d55cad0c340aef2fd0fc508b400b7ef73dbe09201be0  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part15/input.txt
afd38361fa9af822fbb0497a21690d663062d560d66fef6d66c5ea0cefec7f9f  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part12/case.json
34167692fe9152dbd467a1ba46b0818ba2ff5bdccd3980a4a88883f19f98c647  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part12/input.txt
987a05f6fe0fb5b5ca5abcdcb0f0b2c2075bc7093befe669873f9e960bc52ca8  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part13/case.json
91b80814445fffd2176b95797411b5747d0dae6039ff0ed818710b93ead5f96b  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part13/input.txt
cfde3edf3abc235ed7f440857540681f722fe4d3714c06cb101d279a9f67f8d2  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part14/case.json
237052b3f92f6e7704f415cf567e9d6f4514108dadac82d7e56d935aa8591fac  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part14/input.txt
b3cac6dcdf1c0698837ef180cfa4e0ef404710cfc41ef53fad4ed9a668106673  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part15/case.json
37d7c40dc7f79805db4a2e7694f05d82dfab6b92674be6351c282df613843672  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part15/input.txt
c41c01239800422e9e76404852cd6e92be9e0eb7b06a2ae03731c35ad4e0df70  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part10/case.json
881ea299529fb58490a272bb657094da0519de51eaf8775f8adcc22caa3df5cd  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part10/input.txt
16a70de56bbf97cd337caafe0e6e884e31773b2a6d7ddd6490215edc666443b8  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part11/case.json
9435b6923c3c29a133344c9c598eeeb2871e3aa7c2d4df54f6ae5a7a2f279dc1  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part11/input.txt
215a069a02b96e237b303d930e45249f2bfc1ca45f68376474cdaeadc1a44acf  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part8/case.json
d191e9bf51f688460ea16d13c350f53ef0462314a9ee607ec01c06d8b91bcec4  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part8/input.txt
947a6030dd5de2a78312d036c3bcbfacb366354868f3edbfc6ad1aca5d122f55  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part9/case.json
5ae1f96c494433f4e64b7a626cd32b729970fbd9b6dcc1108227a67891105742  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part9/input.txt
284891f546fa72db685ec222fa8ee01966373bc93385b4bdd5f8f1577a5e70de  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part12/case.json
22e37a70decff663deb8800a9302790aece287dfb30278dcfdba46c875c0d67e  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part12/input.txt
5c3095620d507ce9c3e8b2678cc3f169c44d79354aacfe0e18870fcecf4ed018  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part13/case.json
334f3bee308b6618a122b38c9e4dbeb9761df6ad3fade2f2eb0bf9b82c8fa9cf  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part13/input.txt
c81c40f935018617753fee351a93407088c2b6158c346934cd8faff8b3607cc8  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part14/case.json
3e3c35b10dd523095b1aac274b5a04bdd74402875902b0c8518c9740452512ef  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part14/input.txt
da1cc86f3d1f0dda447644c3abddd6b317bbaeb89278c533ec92e5f12051aba7  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part15/case.json
d9b353852a26716ea2162bf35b1c6addd611bb86219b20ac62b9cc2451b2c12d  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part15/input.txt
b7eea448db342d6e77182e767fb5a22d9ee3e02611286c52083042f23283a1ec  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part10/case.json
840251556fd4ae0ddd91068781d3f773ddf3b6290c2e43d42c195a6b8f478365  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part10/input.txt
341e125ac2c13b9a349235fe8df092b0f3e1c441c1a4abd3cdfcda5c9c2e0188  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part11/case.json
35efd84a8e167ec33ac0d7f1eac1dafbce14ae26b44ff046cb894facd1e1903c  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part11/input.txt
426cef2b188379dbc8af7ff26b3fd5920db78bd579ea3701007feb049de8b4ae  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part8/case.json
eb6c44b71f17e6331b8f6a46358646e04365f6e5278b10a02419d998652e008c  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part8/input.txt
65997bb78dce9a628fddae3a329d4b2c9f9985ebce8cca76ca8591bfcab5d103  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part9/case.json
aff8ba2f801efdf4d2d8b9435beff42a6b1363c0eeceab3b45d3dd8bc93342d5  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part9/input.txt
6fc11c19824a31fb679c929428c030041d5d8cde1e389f311110108d1930a5b8  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part12/case.json
3b34bae82d7d7ee984600f8987a9c4988129859a8039adc1517a1537f3ce7ea0  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part12/input.txt
3cf7446965a19ea8ca4d1cfc7b6668fdcbdda1c94a57629da3ffba3e8d15f25d  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part13/case.json
4dd87d66fe41b797fe52a34f57d27e159e2a0abda9be7a0988f88a86c210fe70  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part13/input.txt
a66a004044fe594fb1568b596355bdcc8a9382832f2ef68a72b12b2105dc6699  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part14/case.json
06b190a441d4579eb29d20f98ebba3bcbe6a61b89db43caca5634aedf133f45b  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part14/input.txt
0c05b5971dfecda91b7109a96ad7b2a713822a1b8942c009258040e9cc556432  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part15/case.json
b2582f85d211695ab69bc0d9b83f9e54f5d8865f76a17aac1ef6e07a5b4d8da3  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part15/input.txt
e0683f847b78cd347de80f86ce76e8d15d2269c98af53277d07092a6c9eb65de  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part10/case.json
19b2d615a01adf884c77de668afc92aeb4462d18de34c7fb6fb753219900fe10  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part10/input.txt
5d9226cd775a719d6ca357796e6b7d21a7c3d13e3334e545ba36258743b03642  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part11/case.json
127a8abd0bd82a576930ab7b8bdabf5a8548645d715e2ae6b8c2d995fb068afb  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part11/input.txt
946b1f23c846a7b185b5ed647f880aed22d2a47e069409066babb4bed232743a  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part8/case.json
ad9506cc46a9daa01737977b466237bfb10b863447f24dd3d39cab2ccbd53a7f  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part8/input.txt
aaab61bae9a44bb4f44c9013518fc64d5546a3a59593c656eece1f13dc52ebd2  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part9/case.json
677cd2aa02ea4b10695fef1134e60baa9096c0cf71d0c9282055e8e65f58fc1a  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part9/input.txt
f1fc83b8cc0475e96a016b5887205635b156ca138f798f84a24d8039ba5e3b1f  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part12/case.json
a3b9b76e2f8ee561fe3c06a4dffdfb77dd821fcbd6b70b402b10a3761dc92108  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part12/input.txt
290a60af8843f23070182a46eb25ddb1c785bfd8719ac9dcedb48e76a21654c9  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part13/case.json
e6b01b6cb9f83964301cdef8a216652e8ef3e5bd8f3f3dee26d31356606f0cf4  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part13/input.txt
b25edd224d28162fd9db13e8269290887ba89589ac3c61c73e9fc22710f0b512  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part14/case.json
83bc1bdb1db8b3409753f2c4f11225dc3e97080b8d109b59d2c1e32c7848f4ed  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part14/input.txt
6f985691a2948673971608e6f04f60b0cba8659ab91f11ab86cca6824f811b95  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part15/case.json
bbc55a1d4c338324acae2b785ecccb003b315d3dfe25d282c56aa53d60905257  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part15/input.txt
48eea0640cf36461e4ad1fd4efda8cc0591bd94519b803b5b92f7dedd43afd58  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part10/case.json
b0222794cd5972a39133eac30a67a827f48fbdc60014eccdf30903dfe078740a  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part10/input.txt
09da0bca46e11de3d702c2ef85caadcf13e029997a2f5f0f2032d2f2b1aa1530  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part11/case.json
8cede8061fa99d9f21655eb05fc8b60031d441b151cec7bebef8b0cc64d874d5  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part11/input.txt
e2dd451a41002321fc2c6f16877870f629e5ba6cf4743d87e5d061c8f7a4a588  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part8/case.json
05c86ab8e9a7fc778693d760886aeffc7d0d32aaf680688bbe962ee4fd1c969d  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part8/input.txt
effc044b69c2e0ef1855e2b5e9596c08514a95ee874fea71b1c3e51b11fd3f2d  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part9/case.json
662e00bd79b0a6dfc5f3061204564373e95fa82424136480ac10cbc430a9463d  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part9/input.txt
f81f382d3ee31d505d08ea30b63506a3055df72e231a5be88c7d2bcf7928ddac  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-ecut1km-timing-smoke/case.json
95ca1a53958a836ac8ff39b96dbc78645aec5d3c2f7e982a9f80f9a68e273781  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-ecut1km-timing-smoke/input.txt
49bc7f362ad9ab052c2d2c33e41aa5f2aee0df7ee45ac1fc4124c6a527818f7f  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part12/case.json
0310fd56082c3387bc0010a72441be28ba14e2b46724f1db76824889926ebf12  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part12/input.txt
45e044e8e7e7f27d6a20a004339c1a29f10e073ade95fd2e1da478648b42f315  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part13/case.json
9db2a3bfbf3fbc6188a145272925a32a2ff74dd19c97bab537cf139612467fe4  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part13/input.txt
244df15aadd4e35fdd52245448dae4c7224dac2449eafbca2eee25a89edc9bf8  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part14/case.json
f116a59c9b9837c1ce77565161ec45ee2c1b78c9d685ba1238fa35d29ea8b311  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part14/input.txt
4a0f0df594fab8c33e010b1f4c2d90f9ac9ff6d2e5f0fbb5b750798d3ccc25e5  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part15/case.json
78dd1aef946e7cba6dd34575d17dd8cca2d0f7aad69f5826aa30230a271e4684  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part15/input.txt
```

## Appendix A (2026-11-09, ListProcesses). Bundle hashes of the changed files

Regenerated by `validation/scripts/ion_em/gen_r2_bundles.py` after Codex REVIEW-6b4731a7: every evidence bundle
(192) and the four R2-SIZE-1 smokes (20600424–20600427) gain, in `input.txt` only, the line
`b:Ph/ListProcesses = "True"` with a one-line comment (the per-particle process list makes the executed module set,
`g4decay` included, verifiable from the log). No `case.json` changed. The thirteen smokes already run
(20600411–20600423) are byte-identical. Every file not listed here keeps its hash from the latest listing above
(the R2-SIZE-1 listing, else the preparation listing). 209 bundles, 418 files.

Bundle-set digest of the full 418-file set (same definition as above; `gen_r2_bundles.py --check`):
`6d1ae3a06211831cfb7e8924c65d6d8fe31fae7781443ff84c077204045ea790` (418 files). Changed files (196 `input.txt`, sorted
by path):

```
cf44bda6a4f06d72a62cc740f9b9c0ed9961ad0660d64164c631d3d89f5b3528  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part0/input.txt
cb260c7230f3e02f3040aaa368b86a334856699e0c37218b5e5f33e559f845eb  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part1/input.txt
f2a6df410c7fcdf3106931a02dd9f95da9a7afd5561ffa4f1c39cf8420e19514  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part2/input.txt
1987df7558f8665603d1e95689bd0f9b4bdd28a5e8eebcd57379fafd01fa1a3e  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-em-part3/input.txt
9b3585687cc7551c78f16f867a89507efda95d846af85a4919661958270dd567  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part10/input.txt
753fb9b15d343548a8d59f0c14a948312347cf8f25b77226f3b7bd74e34bf713  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part11/input.txt
a2cc9eec25614e34a0c3fb2db3843da7f0ec1ac63ffae47e327b37231ec17f46  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part8/input.txt
0e1cdad501248a1840197e77f817878b15778f39a0c7a6ed1e21ab698d62284e  validation/reference_cases/topas/carbon-water-150mevu-r2-idd-full-part9/input.txt
fa496402cf11e6b5d7f204e82cc43e81dbd41ddf0903178ac233737abeb79952  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part4/input.txt
d69b330e259d1c5ba753b68b8ad4ab3ec96e35275f77617b93df87f4f0369301  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part5/input.txt
b3aba111786e1462c2da766595e341f949216f13a560b17f311818af0dda07b4  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part6/input.txt
8cad04a593b0562c7bf892fea88cfa6152403b29ee82dc4e8516ca81c9a788d2  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-em-part7/input.txt
541547dae768ab2bdb28f7a6a6c37f92aeb2b9a6e559941f1d98dce9440c71bd  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part12/input.txt
4b105281c1cd48791572d9456d86818dcb037d48084582be2f93acdbe29674ad  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part13/input.txt
198661c784944ed0999051dc3dcf62024a7c19bbf97d4649025c522873cfae8c  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part14/input.txt
d2168e8bc57b784ed2f7d9fa9083c0ed4b55c2af9babea9d6b7ac4b710821d82  validation/reference_cases/topas/carbon-water-150mevu-r2-lat-full-part15/input.txt
bcfc235e65a1944fadfe8d4273674453b48a627a3547b03ff517e1d2336681f4  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part0/input.txt
8881c2713590cd3294cdbaa3c5b55f251fe1fde95c0ea04e0ec6558d2dd46c96  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part1/input.txt
121e61da216bb4c5668d9b96c7eff6330760296df59f84a70740cb08987aa0a8  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part2/input.txt
798535a6f4c7bf0df84e1d59c2051dca8419153fdb957279a9bfc1dc06825d6b  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-part3/input.txt
72d683a7fd10d078faf08329c6f1d84600237211c62d15733c16f070b203feaf  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part10/input.txt
389688b03fb882216adf59aa1c26604d5c9d2820716351b9df11a33826eaf7b7  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part11/input.txt
6eeff6b47575ed2e9be54635cd8a369967320fb7c6e43a70d39b754c7c5aad94  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part8/input.txt
0378dfcfce5bb502c5f0ff43d97e80668ff6ca1d8d317caed29fbe62842e63f4  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-full-part9/input.txt
e1f0729eaa3a7a257b12ea0748fb4c54271183d0a4533cd5547747eac0c25035  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part4/input.txt
373bd257845f1835cc27b0eac6cdec8a5d7f8a8b9c72aa68b41b3b4c9e8eeaa2  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part5/input.txt
0565c8684fb0e8f9166d5eda4574430d70677f5f7b66e85bc4afd3b219431f61  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part6/input.txt
4aa5678f08104a3be4b721979314eb122d7ed46557a6609e227168396245ce46  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-part7/input.txt
37f21fee0a6e64b3a90e9c67e0865a72bab4f6764234cce41db210872fa76ddb  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-em-steps-verbose-smoke/input.txt
fcbe2efe1926b1c33bf04a8e4f8b2eee8dae7eceb362ed8f6123f617b7ac2450  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part12/input.txt
2455aaf2b1477e148d463d2d98ccf6ce81a6f6f5eccc38c125de5fdf8862f35d  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part13/input.txt
896e5604063ddcd78a9a59ec9b62813cec3bd54d9dc29b9bbe633713cd203c74  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part14/input.txt
53cf1aff3074dd53627e4cbcb31941b5a0fbac3aaf141abb443f2568ac8ab7c0  validation/reference_cases/topas/carbon-water-290mevu-r2-lat-full-part15/input.txt
91c8115f2d1728cc414e81b3e2e6889e828ad023dd83cff46cfa13fd6c72d195  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part0/input.txt
9fff9235f52ffcbb591c5d06c85fcaf25c67c591aea1f54d87a938258b6ae9d2  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part1/input.txt
da2c7e6033d1746209764904b81b9b2d79c407154340c2768a8d01a8a2d878f1  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part2/input.txt
9e446a7b9129a6e45965dbee8e834e233885ebbb59a3fd1a1e1495d68e0640dc  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-part3/input.txt
4b0f039293051565b2410d110ef6f4fa43ece4120a975549779927f3500dfda9  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part10/input.txt
577ae4032545b1e6c0d4398d85c78c34778c7464e749d2d5d655f6c101a1cfd5  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part11/input.txt
d4ee9dc17d92e703da8e1b92cdadd9b3171b595d0d13629904ca094fad708a1c  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part8/input.txt
4fc2c7fced547797981a0dcf2548f2241b42ae58e2363069df799365c8f3c137  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part9/input.txt
b54a7004237aae51790e7da4913418f19bcf8660152f794443294b120832b690  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part4/input.txt
b0fcd4c0cb9a735136da6185e30fa8c2e3e4a93fac549f0d4efd431460e766ff  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part5/input.txt
baf28aba40f3a9810fe783149796a9e8d644c3bddb8e6b44ce790ac7b4afd597  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part6/input.txt
2282e1070ef546a90a377f02814430a363eb79b2b080d6779e73a355ea104b11  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part7/input.txt
f66cd544dcca03d43a90bbb1984a8e3e7267d611e613888e853caff50b2c2ddd  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-ecut1km-timing-smoke/input.txt
80f48e1ce57f1ddbf115b7e8199376bf642d696e828fd6d4f354e9659de8a2aa  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part12/input.txt
12c1b849771207b3fe2e9554c3096d05a7fa2f0327499af52b8a7755d19dccbb  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part13/input.txt
c3ea15f972c9239ca813bcc4a8bab94e3ad8b843f682df538bb1dd43ab86a148  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part14/input.txt
8b4199a641b4ac9f055cf1b60ae4a869ee886bae27b4179ccc0db7f20c8b085a  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part15/input.txt
060205d1a56bc1b3ad095c59267dfbd6c44c44fca6009f2900af6c9e27783958  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part0/input.txt
f0d3d41f6d82b70ba26caf4570ce8fe115fd0a5570f84d10011dd5ca8239ed92  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part1/input.txt
72bffa39310a802b9e71fbf014939d95b5626d3e878eacf7ea12a047a8061b8a  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part2/input.txt
cdec42dd7165f2814499d6aea08e99867ba6c539f61c6daae64fca55e6fbb32f  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-em-part3/input.txt
57b1dd03b456fb4541ec77d02f927c5c1146c62d999fec62ccb20a6f1c7214a9  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part10/input.txt
23f7dcdcd9d922218cfec9e7912475e81f8334261980b7f2d03a79974f693c16  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part11/input.txt
e5cbaa389692e8ef1c0e68bdaed0962465fe34ae139c888fac02636bfb33304d  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part8/input.txt
d71085a117e26d27fb7bc2d1ea79df4c289246d6d2d92a11fa10197ee86bdb44  validation/reference_cases/topas/carbon-water-80mevu-r2-idd-full-part9/input.txt
6bfe3a7e126531373582fa81a0e4985600618443712ed7da80ae4ea4fabfbbc7  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part4/input.txt
27b2b95403f421c6b0fa649e7b5106377fc7e8d253c1528d75a3deb759fd12cd  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part5/input.txt
3b500f94caec8770ea3682f7bf26dee659fa2878f0e60422687536bc3ba160b1  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part6/input.txt
ad322cb9719da9ca221b3a7731a6c326e9467f93953f54f57b73d9d07d081bf6  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-em-part7/input.txt
5c82e4e4edcbbc988736f5277a42cff7aa42683f39dba529e092467375ca031f  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part12/input.txt
05109fc654d8e053e8890febeae3f98b55593e94fdc60bc818adaa0cb7a3274e  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part13/input.txt
b31ae3da28b29ad373e2cee9bd98b3369e7728080285d134770a4ad394b3fb52  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part14/input.txt
85fb8cafe5d248a4018b43fa868572b8e194dec9ea980d007569199ae77d66fb  validation/reference_cases/topas/carbon-water-80mevu-r2-lat-full-part15/input.txt
8334385b487b83eab7f69f46b75f64aa72711dc78b878b540b787ede7f33bb94  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part0/input.txt
d076b6747b2ba88f13ca39cd9e57c2cf1ab94aa60a9810aa1e7019a53421be68  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part1/input.txt
0d4c6befa5a99f6f7f68d26f1d57de62c64b0e765e782dc8ae26c3bc52f0a970  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part2/input.txt
89cf59d7862d045989571ca3a1a1729c1c09b5b8c98c1cd449a485487341d925  validation/reference_cases/topas/he-water-100mevu-r2-idd-em-part3/input.txt
da2bb9cb14fcd5347e1c2c79288ad2978f5d1bd219551af6a24283d13e694db6  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part10/input.txt
b070ac00af19d8522670102d1d0d7fc3e378cc62b0a86edcb65e0f0e2632a122  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part11/input.txt
9a83749133826aee08d14f315fd65bca3af456f8ba7d686bc94ac23e549c583e  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part8/input.txt
111ff3fcc594dc244eb271510cca8a313d348b0f86c3204872ce671c4ff905f0  validation/reference_cases/topas/he-water-100mevu-r2-idd-full-part9/input.txt
91c781dcf87b8dfc7def9ed3e91441e73845d51a7e8b60fb3a9dc2482698f6e1  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part4/input.txt
13276a32e985b0bb29668c3de05866265c925a063183a3e4567bd5266925e4b8  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part5/input.txt
3f249288075938f5f0d331d26a7590fb94a6062d7e79214f0258be816eb7c5c2  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part6/input.txt
7c37c5a5da95501544eb8311a77805bacbf52eb448cb2fbe2c04b75d75cb5b49  validation/reference_cases/topas/he-water-100mevu-r2-lat-em-part7/input.txt
17d1bef3c713c51580a1fa8462816e76d93b1cea81755c68293268d00f204a1a  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part12/input.txt
762889ab41f798b093d6c993a68151b3f1760da476fb907147803382ea39012a  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part13/input.txt
16b0d66bd4803b91372a16873d16ca4497576a2397422f22110e07e3cc4493e9  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part14/input.txt
88bd5e81ff177a7bdae0d83401e3e490e6e82d0d120ce5e91a029738d7c879ab  validation/reference_cases/topas/he-water-100mevu-r2-lat-full-part15/input.txt
11c49a5e733cdab0157abd6ca6bfdb2b95511458dc06d7bf23386a839eebc079  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part0/input.txt
814bbdd83bc5c6aa8d8ab98bcf88cfa71fd1ed7037ee0d60ec478b4efc9d8ca1  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part1/input.txt
8eb0c0e8c2c386f6c578a2a99766e67d082d77de123bbb8a4cc05f0a8f1c875c  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part2/input.txt
1be015e4cbeb1dc0ce0c8b1f17af948b4faef3c51fb9cf8297126d1896d93703  validation/reference_cases/topas/he-water-150mevu-r2-idd-em-part3/input.txt
7a86b789e223091bbfe3a499659bb972add935ce2973e8b4ef063ae7b0005f8a  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part10/input.txt
4d13f90673a7368f6e4037975cebbc9beaadefd5a2a45649058450cedc1a757c  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part11/input.txt
da6b025c7fb4a3008ba55d8900be0bd80c528fb4311bbf464d9ea2a007913768  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part8/input.txt
0f825308e094c20c543c37e2f4ea1b4826eab9b2c8f2bdbe1ca98edcce52a5e3  validation/reference_cases/topas/he-water-150mevu-r2-idd-full-part9/input.txt
5c81b06b84b6de936d785afc19898452bf610ee55b7a300114012de360dfdc09  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part4/input.txt
7e9f94ee3ae61a595a19bd73aeca8dbfce6ebb25027b4f8cfa83736df43c38b7  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part5/input.txt
c58c35d77097b0ed1e0aab498c52e8d2f26c524be82a164072654bc05e322841  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part6/input.txt
0e4589ba7e48c49da5b27e003bab18c5cabdb51c342e4b55391e4d99a6c54cee  validation/reference_cases/topas/he-water-150mevu-r2-lat-em-part7/input.txt
4d06c05e6dfea6fb183c5b0736109f4729ce7f6279584d334e952202d0c9549b  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part12/input.txt
450e2163bf2a7170344195fa681019010428156fcb34ee8789a856e59e1fc811  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part13/input.txt
c5ab6f1377e96509fb31b4b2be04808c0338c367f4729cfdb5bc1bbbba309901  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part14/input.txt
edb95621fe2f9e041f1b3b3aec5a596a52cb6f3fba2af22b050e6a6a6304c799  validation/reference_cases/topas/he-water-150mevu-r2-lat-full-part15/input.txt
252d99f4c9011c87cb2a34cdb3f3dd0f17d829c0a7b635b4ec6d3af2d5a085f0  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part0/input.txt
4a0d282deb4d611e9115fac727d84a543959fbcd852c398fe98880b3ea314c95  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part1/input.txt
79e758aa985af052fe077f7a91837a54b99f4fe7a3ade750aa75c93ddeee94d2  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part2/input.txt
88e9642e1ea740d01b6c633074f526cc8bb5a59bdebba2506e5ba00d2e494f6c  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-part3/input.txt
83cf5e5fa8f1b222ac38bc2626f06a67b1bb9b8850d763711d03ec626eaef33b  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part10/input.txt
9284cf6c7be05e74439ef26b1e02b5e8e79311ff1b97aa212e2b40172fc7de07  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part11/input.txt
69944ff96c33ebcea4381522284788f9e6243eb611befeec10425e050864c92f  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part8/input.txt
0016ce9aea5b52dfe82e8aa40fc52002b8a6d8457d2f8c836e3e7ab2c0c04ee6  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part9/input.txt
95f073539d7115fa89b5be457a07eb1febcca5411a7d34c93f7948b785c10937  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part4/input.txt
67ad45437c5ee7728efb8add6ac5193a3069c860019c1c391a93c024d39a7af4  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part5/input.txt
63399f9df331e41769a3e5a1a229fb76b70cdc99a1a0808c4592c69ec142668a  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part6/input.txt
3f06effbbea19dedfebe62c2023fa15ce7b1a1d169e1a92bb17f428792f6a3f0  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part7/input.txt
88ebffbe8c7bbef064fc324adc4a71347f411582ce6111e55c836da578370bb5  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-ecut1km-timing-smoke/input.txt
29d121b47b8044b9baae489ffca123b5f3b442dd287718d51fe557d09f547767  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part12/input.txt
643c4662743da8dc61e387fa502c39c44737431f46256ba8508de843b240331a  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part13/input.txt
3c8e04f7da93a5bcddcabb5b493cf1392c70ec53b8d19f1980921c63c1a96834  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part14/input.txt
259cef9309724796c2294e030ac4a49018e2ac3ac18eccb02be3db84ba4491bf  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part15/input.txt
93b0df2c7baec99bc583bd1b9918167e03e84c17d988f4b7462dddba4cfbee00  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part0/input.txt
3bbf046af0b3754e402942a9ca79af23d1e56a1c30a7fbbc22392cef09def090  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part1/input.txt
93b4ec5b4d37ccc4e3a834a2868f38793de14dea0dc51e6440e6823d0570dd08  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part2/input.txt
bd6502caa971dcf050a8baf314ddeab7e07977234b3bd9b5126f42b64e074dc6  validation/reference_cases/topas/he-water-50mevu-r2-idd-em-part3/input.txt
8ef7c84223f146a23b7d7bb34bb57b8fa8660e31c768396f4b2afe7c70cfbd9c  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part10/input.txt
0363c7105826f36b89a9b04602acd2a6d2c18590639ab64af84b3f1d2fc4da2d  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part11/input.txt
aa0c404c0390bcfb1ef5da43e69669b08d9e2bde2f197e2eb54354c5e1c6200f  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part8/input.txt
6789f74ca78201c9c5d2c60cb6c56f1eaec49c7840592ac8a7ec4fdee9322760  validation/reference_cases/topas/he-water-50mevu-r2-idd-full-part9/input.txt
5ced92db423af1699d05f9b076b6dcec2d324c459ae772a1ca3d99218ebe7bb2  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part4/input.txt
ffd3bad235ca868ac25a94a861445b93f3340971e61766ba2c926d24f3d840b4  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part5/input.txt
7280f3660c12838f24aa79be3bfa3db7aa342361e6fe5a5686dd0a8e4f35dded  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part6/input.txt
10d7b9a8f372f2fa4748b7d9af7f721d0d3162f5f821f890e81d6d7d6e8ef221  validation/reference_cases/topas/he-water-50mevu-r2-lat-em-part7/input.txt
df6c9b4abd83cf01694d96ff429fad9f8da47d9b5059885c420bd49167cf6b53  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part12/input.txt
ef088ddba5a74c5ddca352c57ef604ad4864b4a63a946bd3cf1d77da24c2af71  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part13/input.txt
e83f5ef441de1555c9dc9092746772616b9aabe608e75893c27a7d95a84b0fa1  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part14/input.txt
26d9906ed9459fdd08015115e7803fad455e2ae970e58fccb9e69a5595da4207  validation/reference_cases/topas/he-water-50mevu-r2-lat-full-part15/input.txt
0f369d36dc7c17e2253ff6032b505e15744f4165e75da79c1dde03f34f094dd5  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part0/input.txt
fea806f16d2949a6bb2a5deae85563a250d3777caa2ea7ec014351f8e39ed275  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part1/input.txt
29302a7aa18588782cdd5049238fd6d4746a3537fa9bb3f43976b53264524e2d  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part2/input.txt
4c7ef7d5ef3ef3391a153837ae283a736625809adb98bcfbc707dcd55b7fac75  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-em-part3/input.txt
ed7471475e1344a4e4b0ea6014d428e59f0b82b5f198e6e6b29b122a71e507a3  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part10/input.txt
0cd59521d2ec231a25ea0744522727a13383e4e3acf8acb57bca5e03ad7c5d48  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part11/input.txt
c93c24915f68d568f3b4263664b54bc3c185b04bcbea9d06727df5d86873fe79  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part8/input.txt
6d8d10f7ceda306a2170e04fae8b0e71e33043c835dcda308c6f1670ea9cdb0a  validation/reference_cases/topas/oxygen-water-100mevu-r2-idd-full-part9/input.txt
4279ba95a47af99f06e8ce5d6211feb01199b8469368ab5c47853e4191e8ec00  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part4/input.txt
f450cb7d78ceb788b5392bb86396b7b402986dc015cdccf1354e2ebcbba5bcfd  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part5/input.txt
95162255797bcd85266d2f5334f56c6e7b97f9b154d9da911c50ec57e4a38e0a  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part6/input.txt
f0732f5181171c6687f16fbea86c55e1633d00e7b8710b85df72f78b3705dd0d  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-em-part7/input.txt
71d7b9c694c18bb76b7f79f5880789a347c48bbb5274126f44ce4a15f4e065e6  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part12/input.txt
e0f1ee6fd49b37bbd41085d34f4075f30ba1ec8c18f5000a201a527cc52fa581  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part13/input.txt
34de47677fbcd820b2ec364a7bd9f70fb0c7d1b6abe90672750cfb390aa067ce  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part14/input.txt
cce57bc927c6d47d9247d7eba60959062762ca33d91b83b8e5b0c1473ea83d18  validation/reference_cases/topas/oxygen-water-100mevu-r2-lat-full-part15/input.txt
0c90f7480c420e18359deb770cc7ff6620a3b5854be40e930661bcaf1bbafc94  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part0/input.txt
3ba88761ac631d5883ebed3a89acd009e57b260053df350e2d1f8265c7ff0940  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part1/input.txt
4160c96b0630fa76fcc9ded08e24e8b4086bf40c528d36545b60c19730ffa58c  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part2/input.txt
d5d901b47b817403407b9aef6c78408da962b82a2fdf81e8b27adbd60d3981b1  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-em-part3/input.txt
faa1b146bfecdd18e405db87f1c508fa5269f21d0698fb7d7ed3ba7b4b4f3b37  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part10/input.txt
563064739f6ae2c0993b4764cbb10f6e908fa8f35d4c10c189f1d7216f1cd6af  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part11/input.txt
58276c8c0b0cb0b04fe1712f57b5f81bcd6d93e51058d92850aaebad3512fabf  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part8/input.txt
37283f34eb3e190ac0bae47150b23faf577066bdebb5949a6dbdd92f8f8e571f  validation/reference_cases/topas/oxygen-water-200mevu-r2-idd-full-part9/input.txt
4f1c63dccd75bd981c8e35b8982c273d86de95bfed74e11c6fc1938210121856  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part4/input.txt
0d652ed857313d0d7d3f0763cea1cee3d3513b612eb77e097d31997ae0fc6f47  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part5/input.txt
b7615f2cbea5e2c00c3e1725456e2da3cb631abdda8f3e50714851408beb40ab  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part6/input.txt
75eaad89c46a9c3deb84f5d5635d6b66e42418b8df5f1e1c15ca23bea7733df5  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-em-part7/input.txt
cc3c11b0190eef3e076f609ea3926dd7c69f02df7a64c09a23c1179ed72aa306  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part12/input.txt
022555a253ad4a2c551446ddf29ac5ea7526ea87521be160e38f7c83e3e0fa18  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part13/input.txt
298c8844983f3f9862b5c7eeff519005a5e7a2fe5a74f7487b9c3fa29b4e7e9f  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part14/input.txt
fa0c83cb9b63a8fa3c74ed339a26b1651e58ca723eb155950623829cac466a6c  validation/reference_cases/topas/oxygen-water-200mevu-r2-lat-full-part15/input.txt
0e984c85d4b403ec571bc33f8126e05b93d039609f8ef86d3a0abdb66a4d53d8  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part0/input.txt
dba46b90ad5acbc3fb5020026c5dd92f93f794c9cca8056852c03aa678f232e6  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part1/input.txt
9cbc96e42fd9894037d0ddbdd7890f9df0734d192a0bbcadb0052babe337c5e9  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part2/input.txt
59fb7337d8847f7e0cdad091105099495e69a27debde54838f422d6bbc210df7  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-em-part3/input.txt
ea1aad9897e5adba2c18bb85e068aab9d99e25219180bc7679e45994f8f5c0b0  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part10/input.txt
d8b00049592623da1d85b3c4469e1bac6493661503c4f0df7889ce1218426a53  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part11/input.txt
bd349cf772050455b707688caf51cf39b88941b9a21536beb8967741f5633e93  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part8/input.txt
aa22c8f9cdbec38395c319b968471b1316aff4e9a7d2d89c4d8c2c306b6b266c  validation/reference_cases/topas/oxygen-water-350mevu-r2-idd-full-part9/input.txt
be6b90ae5a8f63272fbde95da42f4d85551aa95b7f3dc0c9112ac2d3ea33bc65  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part4/input.txt
8dc7921d99234a033c2b6965771e61f1406887f57afe0471932d9f6f7cdd9505  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part5/input.txt
1efa9afad66d7b301bbd36a3505545f49822841d89e44d72dfab9dea8394dd80  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part6/input.txt
43ea338f29a8501d5a920dc0e8c26aff94a209356803a94cfae203c1cb9adb87  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-em-part7/input.txt
21fc11a1d864b74e98969f8515b0acadc4b8adeb184b2f8197c590f8b88abaef  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part12/input.txt
056357acecc0f1ae30e0d13e46111a3a3665e55919848dc6acca424e90d105f6  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part13/input.txt
414f2ed514687fc3198fa6a97d9ce9ec754a8e5c6328dd3c944462f34ce19a20  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part14/input.txt
aabc8fb1b1beca52f9dcddc50c3a21a358816b9bf33e191ae11b93e9c87368c9  validation/reference_cases/topas/oxygen-water-350mevu-r2-lat-full-part15/input.txt
81b8ddbea43835da50c7d9a5c317f434c55f7396b635b91b2c0391734b57089d  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part0/input.txt
e80147ab991bf3b335086c88ee5caf40435a6ad77cffce292da45d17d0b17b0e  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part1/input.txt
b540fa5ca3c1f83ec45c0714896be45e76dc7bb2cb0b5d17218c2489471f255a  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part2/input.txt
d1bba3f7f04fba46509d82ccff29a34ac7b156372bc65e67e75128e87ee9eb88  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-part3/input.txt
74d9f6056005ce1060260e2e4e2344fb78e017deb0bed13094275fe21c94a5fb  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part10/input.txt
943ad022cad47c4e4c99625467f22eab3af23a34b11eccc629c77704e50c4829  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part11/input.txt
df52e8373829274b50a47ce20f0183bcd77571e56909e408f9766bf7f6523852  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part8/input.txt
9e88d3b39815a1ea4680863f453afeef7624f58bf6eba26c4b047e59dc2088b4  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part9/input.txt
a17a0909cdea1015291e4c8dff705df1a1f2d77a609db638b653152f1b00a5a6  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part4/input.txt
797d1c1b51cb327521918a970797ae65d4294cfc0b507916bfccb2e88d6beeb7  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part5/input.txt
cdef636b3c1364f69755cd58bd4b1cae82a2c80cab815f85f91c0a83e473e2fb  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part6/input.txt
4f02baec9368b5867ade352f769ac5299e5e99cb76d634c646e4738ec89ad69f  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part7/input.txt
21bbac3a4881ad9516eafd8cff5748dce5dbfb7f1ff93104e42d06e8fd8e6175  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-ecut1km-timing-smoke/input.txt
c5adf48eb5166ddbfc30dc7c9a2b755c0c105150ec26b66d5077ba1cd9cf5be6  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part12/input.txt
9ef37716845f846a740c342babf49c0add319b4058b97446918f47d8614a57b1  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part13/input.txt
dc5e57b51e6ec490f68e2487cdb9f28d4f5c92f68164e179bfa9e49ee7284307  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part14/input.txt
1de7a06bafe88677c4d41894d8e4fd74379ab0dc8e5bdf1662e3cae090d34e9e  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part15/input.txt
```
