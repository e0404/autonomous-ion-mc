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

**Order (pre-registered).**
1. Revision 2 frozen by a reviewed commit (this record, the bundles and the P9 additions bound).
2. Execution-only smokes (rehearsal seeds 20600411–20600423; every scorer output deleted unread, P4 pattern):
   twelve timing smokes (He 250, C 430, O 430 × IDD-EM, LAT-EM, IDD-FULL, LAT-FULL; 1e3 primaries EM, 200 FULL) and
   one step-length smoke (C 290 IDD-EM configuration, 3 histories, `TrackingVerbosity` 1). Recorded: wall time,
   module acceptance (in particular `g4radioactivedecay` and `g4ion-QMD` in the seven-module list), process list
   (UrbanMsc for alpha/GenericIon), killed tracks, and whether scorer bins limit the primary step.
3. Sizing note (dated, before any scored run): measured s/primary t_s. The part plan holds if t_s × N_part ≤ 1200 s
   at the upper end of each ion (the most expensive case of each ion and config). Otherwise, in this order:
   (i) reduce N_part, keeping the totals ≥ IDD-EM 6e4, LAT-EM 2e4, IDD-FULL 1e4, LAT-FULL 1e4 (the statistical
   minima of the table above); (ii) IDD-FULL may add parts 16–17 and LAT-FULL parts 18–19 of the seed formula
   (K = 6 for that config; new bundles committed by the dated note); (iii) otherwise STOP (amendment). If the
   projected REF total exceeds the 16 h cap after (i), V5-SS takes its own pre-registered 4e4 fallback, then STOP
   (orchestrator decision on the budget). No case is dropped silently.
4. Implementation, CI/LV/HR rows, T1-ion and T12-ion (their failures stop everything before any V/T row).
5. ionmc samples of the TD/TL rows at the frozen implementation SHA; their result files and hashes are recorded in the
   evidence archive **before any scored TOPAS output of these rows is materialized or read**. TOPAS runs may execute
   earlier, but their outputs stay archived unread until then.
6. TD-EM, TL-EM graded; then TD-FULL, TL-FULL; then V5-SS (independent of order); then the held-out set (once).

**Run validity (each TOPAS part; a violation invalidates the part, which is re-run at a fresh re-run base declared
by amendment, the invalid part preserved):** zero killed or stuck tracks in stderr; `NumberOfHistoriesInRun`
completed; the EM printout shows UrbanMsc for alpha/GenericIon and the module list as bundled; the replicate gate of
`compare_batches.py` (identical configuration except the seed) holds within each group.

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

## 9. Not established, risks

- Whether TOPAS scoring bins limit steps (OPEN; step-length smoke). It moves the Urban width by up to ≈ 6 % (P2
  §6.3) inside the TL band, and the cost.
- The Urban deficit for ions (10–30 %) and its energy trend are scratch estimates (P2 §6.3; the Urban emulator was
  withdrawn in r5). If the actual deficit exceeds the band (e.g. at the low ends with very short steps), TL (a)
  fails on the reference's structure; that is recorded as FAIL with attribution, never re-graded.
- `g4radioactivedecay` acceptance by OpenTOPAS 4.3 is unverified (timing smoke; rejection = STOP, the operator's
  configuration cannot be changed by the planner).
- FULL cost (fragments, neutrons crossing bins) is uncertain by a factor of several; the sizing rule handles it.
- ionmc fine lateral tallies (0.01–0.02 mm projected grids) are within the scoring design (grids independent of the
  transport grid, up to 4 grids) but their per-step cost is unmeasured.
- σ_R/R planning values and the FULL A_frag derivation are estimates; A_frag is frozen as a number now.
- `validation/reference_cases/README.md` (bound in P9 revision 1) is not updated by this record; a layout entry for
  the r2 bundles can be added by the orchestrator in a later commit.

## Appendix A. Bundle hashes (sha256, `sha256sum`, at preparation)

The listing below is generated by `sha256sum` over every `case.json` and `input.txt` of the 205 r2 bundles, sorted
by path. Bundle-set digest = sha256 of this listing block's text (lines joined with "\n", final newline included).

Bundle-set digest: `9bd6f8023b4aa4b0bebbb21c5e237ea5a71c136dbc5243dbef5ef772c54a212f` (410 files).

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
905c8d3a717b8469ddaf9e446a66805accfcf1c805da05176d82fad25b5bc00c  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-steps-verbose-smoke/case.json
562bdb351570489176010ae44b684222918913e00220bc12dded719cb556ad15  validation/reference_cases/topas/carbon-water-290mevu-r2-idd-em-steps-verbose-smoke/input.txt
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
1bacfca42afc4cc334041b9d1dbaf8ea9dc338cefb9ff12a95480d1942713020  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-timing-smoke/case.json
a4099b16e5c81a1ffaf629744c07f9300b6acecbd43c0bfe30b3471852059823  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-em-timing-smoke/input.txt
382c3c91647efdc1314c68a1285cb9948fc833d91cd0542db6e22be7b48ec39b  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part10/case.json
8ed41c2e9e1d469fd75671de40622a98c37a003d758506d376ab6ef48ff796b7  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part10/input.txt
fe0014d599403d6f14358d5e5ff68e715c2ee4ecc9a8b8b90f0bd661f6f7ef8e  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part11/case.json
3d346adeee9624d394f5182545afa233b7bd0d0c5405d3ebe51a100f41463292  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part11/input.txt
bf16150b9470495d1a640cac3c5d1b3812fd2130dbe48465f5d1c545593c4e7c  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part8/case.json
8de503e2d54df10a1eae4b8b4e3cc5eb46ab48e143ce592a3780b42e822e9d03  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part8/input.txt
069c258f44848b84485bc606ce811c0d8853321a5320850d1d3cfd69df288e0c  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part9/case.json
68ef690e282c4b9e5ee8a4acaa82bde1195e0431ed1ea3bab0b40208ab974adc  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-part9/input.txt
bc323427602bd313f4467145ce070a2e4f4689d45fc7e236bcf9a2633cd4f490  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-timing-smoke/case.json
4fe475021f99f4bbc66a709d8e32daf1410190b7481a7cd8b27323d8f4736043  validation/reference_cases/topas/carbon-water-430mevu-r2-idd-full-timing-smoke/input.txt
a7975cb46a40e8a5eb5b36a352db39aceedf2b0e7d7145708cc565fa452a7d9b  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part4/case.json
a32e18e87d30f184c065d0f93b4622deb2e38e4d7bf26748637d9b8d8c84b178  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part4/input.txt
ae855029b11701a7585f3282b2106141f360573563672d6ec744d41d049d0687  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part5/case.json
1acaeb718f82700479620d4c083e8bbac28c520bd2f9593dda85f50ed171c652  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part5/input.txt
4dc458b7e061532809730f0321c6430c253b330bf6c19d3da6de88fafdf66d4d  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part6/case.json
25ac6fab3bad4f91cb014e514648488913e5e03c5c2216a0f7d77647a6d372ef  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part6/input.txt
882defbf9104e46740e4326faf7d04e62989d09e2e0e35ee275dae2368786464  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part7/case.json
b73c788c01307a4ba3de511039ebe4c89b332c614d217915a758d4f3b0688caa  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-part7/input.txt
460dd249ad154e6aff58d67963ef5bbc5f6256c84c1491e3b813e6d6815ad4e0  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-timing-smoke/case.json
0879efff78ae63e13389e72feb5921b504ae4476ef4b37ef07653d198d01754a  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-em-timing-smoke/input.txt
3c9f46063fb4da211212854ae7cdccd7de02a064d137fe6cc0ca9a4b7cc21f9f  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part12/case.json
2df4287de9921460651af5b571711a482a5c21b4dd2338c4e6de5bd180757cd5  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part12/input.txt
f5ba25d1b555c04225945f9f759dacf51d3e84120b327c41856c848a1f1c6964  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part13/case.json
bf082d2aefe5f7a7342b2ed932a06c025575933bbacf4bcfe50aa48b4501a1c5  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part13/input.txt
ae0677d3cfdc6a2326d7e28f1da3919520fefe6bd1bf232f7fe5e2fec32fb74d  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part14/case.json
edf0ef74cf84741376417691f06e9f8ba8598150a2294b103f8cc9538c5876c6  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part14/input.txt
ca6be186ab24f2ce2650b51663509b190e0019e195ef3aa2382a028f8cccabd2  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part15/case.json
d54a3a6ba350191e4fd0651ae0de2bd5402b3591b4c030a5c13cbdb2017455fb  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-part15/input.txt
ae64307d48ee0391850f14130fac3dbca80d15b1a78dcb64eb1b865943091b68  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-timing-smoke/case.json
7adc35980c353584c5bce11e06c8c2f9c5c43a7be074473a0f4ff1a3da9fd4a2  validation/reference_cases/topas/carbon-water-430mevu-r2-lat-full-timing-smoke/input.txt
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
2b13ddf6982c251c687b6758d563edea5ae33fac68c6f73dd9d4d7d1eb2d181e  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-timing-smoke/case.json
38a85507b8a5ca3a5bc917bf68e593bf466e5d59240fa1679788f26c7bcaa26b  validation/reference_cases/topas/he-water-250mevu-r2-idd-em-timing-smoke/input.txt
79dddeeba7d850d126605408c29b4377f04ffde2a4ef07ce240177bd91554ff8  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part10/case.json
c2fe5a7f39d99789984ad814394ae1e22b975cd284f4ccdc99fe5f3cada8382b  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part10/input.txt
77729bda3bbafc47426e5f6438dd787d12990d95227ff9cc6d0c082eeb200324  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part11/case.json
e8902a4f38c21babf64fb3facfacdef7c39ef820e6c68d3809ee3d96f5cd4f6a  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part11/input.txt
9a22747713c613bad7789c5d9a79410893fa8ff59d09494f1869aa2b55aaac34  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part8/case.json
bb07276391aa51fedcdc6ce5203d4359cd5bd091fce938cc4937b46d02c934dc  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part8/input.txt
bc0ce4406c3feca53725218fbb16955adcec0900ff665bf9a7e5928b83edc4af  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part9/case.json
c0de13d48b54de58baa278b7c1fb6aefe8e6ef3c9670185545b7abcea7dec201  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-part9/input.txt
1839a15d54ad2f1ba884be45119a60490f7a3715aff213e7e017fc2644cfdac5  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-timing-smoke/case.json
7e0140484fe773fde1766319b315c118d3c33558a455525388ee58400cb7f5b6  validation/reference_cases/topas/he-water-250mevu-r2-idd-full-timing-smoke/input.txt
3fe6bd332baa2c9565d6983830eef7f36510840e5cf1d622ee4dc333fe140e03  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part4/case.json
48e0e76b9c51d72b6a946b1d1f90b2f499c896e1067108a00fc5c5b40d4ce83c  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part4/input.txt
88e4c8fc89e21739e266aa272477cdd04eaefdc7de90e4846da30591a867b81c  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part5/case.json
3d0cec3801dee8e2449695c13e2644c2ffa07ea24d4b910227b642dae232b47a  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part5/input.txt
454fd2d8480ee20c0c2ba54ad631035373008d0fec8f127af3a46ff16de8712e  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part6/case.json
15cf79945b42ce82d2cf0dcaf801d7bff7703bddd6ddc30c23f57baee9aa5199  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part6/input.txt
0a13f5c3bfd3583f6dcc35ffe19fff1a5e24ad50e30354bd771ab04c8f1fc467  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part7/case.json
d2e25dd1caf7873ffefd0b2518c7230585bc10077e9a7c80ac3f63b45f1a51fd  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-part7/input.txt
2c8d25150fd9d7019df2aeb99adc02c80f2dbd9647fc80b784f8eb6bb23f1fce  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-timing-smoke/case.json
489777cf893d47d3c6cde16ded7e802099e84109957d6f781b32037382d96500  validation/reference_cases/topas/he-water-250mevu-r2-lat-em-timing-smoke/input.txt
39c4b2134a1cc0c40d09bbb2175a9f95a85844bf07b318ef096d24d65ee74533  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part12/case.json
0224f5e5577501d1b1378912ae83e605da9d7da3930b1894520c839c6eaae13c  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part12/input.txt
6905e1399c87b3a9ed8e4a0640b2e51bb7b0b20b12c4207d26f6b62517b75f65  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part13/case.json
bc9993746ee6a780422f4570afe189b3f09e40c0bd7ca8705d9cbffaf3344230  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part13/input.txt
544bb3037ebf6c9e094b3d2eb2dfa7ba0cab0544feee9759e997d503f6fe643e  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part14/case.json
9b7a6ce6450e98b79fdd948468a6e114c468d53456c2dff2f3144f3ac3f6de4c  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part14/input.txt
9b881e9f4072e3bd66620ff9374fc74da09e5b5e2ef0a20eb4c8d343d6b1b291  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part15/case.json
af94ad4e0d9ce2e90d629d452c936c164b25324f4b11ea4bdc628669aed3c7ed  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-part15/input.txt
9749d2a050c800a409b4363a4f2f7238591ecf6ad631b78a297b76f6b0697e5e  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-timing-smoke/case.json
15516290aca2d55452dd0fc7067598ebc9460079a5dbfdf7b7d2878566e3d112  validation/reference_cases/topas/he-water-250mevu-r2-lat-full-timing-smoke/input.txt
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
0c6074f48f23dac17e03a6e9e5430ad6327b5ba488cfd741262b2a711b9449c3  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-timing-smoke/case.json
f45c3917f36572262d0f4dfda2ac468562e1880c9732c9840cfa83cfda940151  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-em-timing-smoke/input.txt
50c4c7b2a05fc5c32e0122d4a32cce7528abdc3c70a361c99e9e6c1ab8fddd9e  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part10/case.json
52e7bf1a5ba389f7a830828d995b0722859ecb2ccd614f1aab258d3e976de789  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part10/input.txt
c5831bb9708c659d1447468c664378ce3f4f7bfa3fe99283604dd8e6f6b584e1  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part11/case.json
1233d80725147c49a5381090d1bbb01197a9b48494b8432dfbc289db29e024ac  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part11/input.txt
5006a3af39b6a8e7dc4d0cc55bcbc4bf56f3b3839cdf01561b5dff35bd2cdc79  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part8/case.json
e6433b95072be6428a9f7e5d0c7460c93adff2a6f43e534d1a145402b602fb10  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part8/input.txt
66bb11cd4ebe2df13f1ad44a32b4ff24d3719bb3de49c496b6aaba69c74148e5  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part9/case.json
0c3fb51e8dc280920bc4af39713a3705f5d9b013fbecb5b5dc802ceed43d8007  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-part9/input.txt
3c9ef78a621686cf2e27dcc0b2c462d11566d6a77dde2d7cbbd97dafa5df1e18  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-timing-smoke/case.json
95a4277bfce6162c027360f38f54284c0105a75a39dc2b68985edce833a50e5e  validation/reference_cases/topas/oxygen-water-430mevu-r2-idd-full-timing-smoke/input.txt
b6723eece668406022f671b4bcf1eb77646ccf47bfd5cb667edaa659dcf38ebd  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part4/case.json
b3b29c17517892940e928bccca0a5638b3e86b07bd77e69317a298faf96b4bc8  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part4/input.txt
0466de97357aba8b671745dbb541a71a301d59bf7970145c4c0ad066ae8cadb6  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part5/case.json
bf43bed492bd6be6ab0bf519c503503579e47dbafbc695ce4759a04db9da8ee9  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part5/input.txt
b220748eaa26a7a96d99602417ef79510e984d9d9f472addbed162a8ca1a3f88  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part6/case.json
b63656311e6fb7b98bddaebf8300583b84c4dffbd878cdd83154882809279394  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part6/input.txt
0da7b1f42cbe5120a643c975929ccb661ce91e033ce09024b013da234e1c1573  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part7/case.json
92381bbf13b2f4f904ea93245d6002e1d17fccc43ffb84d9ff80e7a18b155479  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-part7/input.txt
9dc44d37cf64acf87251e6e4b7c156fdd216bb2f89435fb09490f8dc383ab396  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-timing-smoke/case.json
21e4387166ceb61bc120a1de230fc04a8e0c06b684105f4420c8a992b4cb31b3  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-em-timing-smoke/input.txt
941e3f46bd855c23dbe2a6002929c913f557a3ad86ae08422a8f309aa9297258  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part12/case.json
aeda25a264763054b34228efe950ecd664b2aa54dd8ad50496fc41b80697c00c  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part12/input.txt
ca657d73e9c73afcf7402c9073179b6c840851d5f491fb67e518479170ddd3ec  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part13/case.json
049facb62b200d6ceaa3641a75371a0d3de51b5e4c52b30d2b084a6bdb990f16  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part13/input.txt
77b4ae2d9b8e509ceae0850d1c831b30eaf3ce951b3ca19ee1dffe5ba2b621f2  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part14/case.json
cb6f5ef65f3c454240722ac5e25150c7ae945145df9fb1f72bbf2d4f93738eb9  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part14/input.txt
2d68e025e788bf6e683eb220ea3dbf68c91706575a6fb2390a88b5fa5b8441e1  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part15/case.json
e1e7902363c3c5735fec1ebe6adfe75a7e982a158c7791a3fb5e8f5c96769b39  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-part15/input.txt
b874a716cba12e8b42b1ddc5f8287036a4a41ccf03c21358a68828abbe8ca75a  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-timing-smoke/case.json
909b10422554ebdd2061e4cb8b7c80d061ea99aa7b8dfa68bdfd8232645af795  validation/reference_cases/topas/oxygen-water-430mevu-r2-lat-full-timing-smoke/input.txt
```
