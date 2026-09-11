# 0027 — Helium-4 transport via effective-charge scaling of the proton table

- Status: accepted
- Date: 2026-09-11
- Task: DEV-022
- Affects: multi-ion architecture, stopping-power data, validation strategy

## Problem

`REQUIREMENTS.md` requires **ion transport** (MUST) and a **multi-ion architecture**
(SHOULD: helium, carbon, oxygen); Stage 5 opens multi-ion transport with milestone
**V5** (helium/carbon Bragg curves versus published data). This decision adds the
first non-proton ion — **helium-4** (alpha, z = 2, A = 4) — reusing the existing
proton CSDA transport, which is already species-agnostic (the `Particle` carries
`charge`, `rest_energy_mev`, `mass_number`, and the engine threads them into Bohr
straggling ∝ z² and Highland MCS ∝ z).

## Decision

1. **Equal-velocity z² stopping-power scaling.** At the same velocity β (same
   energy per nucleon), a bare ion's **mass** stopping power is
   `S_ion(E_ion) = z² · S_p(E_p)` with `E_p = E_ion · (m_p / m_ion)` — the
   first-Born (Bethe) z² scaling; only the projectile charge differs. It acts on
   mass stopping power (MeV·cm²/g), so the **same water table applies directly**.
   For helium-4, `E_p = E_He · MP_OVER_MHE` with
   `MP_OVER_MHE = m_p / m_He = 938.27208816 / 3727.3794066 = 0.2517243`
   (this exact mass ratio, not `E/4`, keeps the z² scaling exact to 4 significant
   figures across the therapeutic range; `E/4` leaves a systematic ~0.5 % error).

2. **Bare effective charge `z = 2`.** Barkas electron-pickup lowers `z_eff` below
   2 only in the sub-2 MeV/u tail (the last ~0.03 mm of a 15.8 cm range), shifting
   the Bragg peak by < 0.02 % of range — far below the ~1 mm range-straggling
   width. The first implementation uses the constant `z_eff = z = 2` (factor 4),
   documented as valid for `E/A ≳ 2 MeV/u` with the residual deposited locally at
   the table floor; a Barkas `z_eff(β)` is a later, purely table-build refinement.

3. **Build-time table transform (kernels unchanged).**
   `scale_ion_stopping_table(proton_table, particle)` builds the ion `StoppingTable`
   once: `E_ion,i = E_p,i · (m_ion/m_p)`, `S_ion,i = z² · S_p,i`, with the PCHIP
   slopes and cumulative CSDA range recomputed on the ion grid (`prepare_stopping_
   table`). Transport then runs the *identical* `linear_stopping_power` /
   `midpoint_energy_loss` / `energy_loss_step_length` with the ion table and
   `particle = ALPHA`; no hot-path change. Straggling (via `charge` → variance ∝ z²)
   and Highland MCS (via `charge` and `m_He` in `pv` → θ₀ ∝ z/A at equal velocity)
   are already generic. The proton table floor of 0.5 MeV maps to `E_He ≈ 1.99 MeV`
   (0.5 MeV/u), so the constant `z = 2` scaling is actually applied down to
   0.5 MeV/u — into the Barkas-pickup band — but that spans only the last
   ~0.03–0.07 mm of a 15.8 cm range, so the residual error is negligible; below the
   floor the residual is deposited locally exactly as the proton path handles its
   own floor.

4. **The range identity.** Integrating `R = ∫dE/S` under the scaling gives
   `R_ion(E_ion) = (m_ion/(z²·m_p)) · R_p(E_ion·m_p/m_ion)` exactly. For helium
   `m_He/(z²·m_p) = 0.99315`, so a helium-4 ion has essentially the **same range in
   g/cm² as a proton of the same energy per nucleon**: a 600 MeV (150 MeV/u) helium
   beam ranges to **≈ 15.86 g/cm² ≈ 15.8 cm** in water — the depth of a 150 MeV
   proton. This is the primary Bragg-curve validation target.

## Validation (`validation/v5_helium.py`, pre-registered per decision 0001)

- **stopping_vs_bethe** — the scaled-PSTAR helium table agrees with the *independent*
  analytic Bethe model (`AnalyticStoppingPower(WATER, ALPHA)`, which reproduces PSTAR to
  < 0.1 % for protons) at several energies spanning `E/A = 10–250 MeV/u` to ≲ 1.5 %
  (the two share ICRU-49 physics; the z² scaling breaks only in the low-energy
  pickup tail).
- **range_identity** — the built table's CSDA range satisfies
  `R_He(E) = (m_He/(z²·m_p))·R_p(E·m_p/m_He)` to the interpolation tolerance, and
  `R_He(600 MeV) ≈ 15.86 g/cm²`.
- **helium_bragg_range** — a 600 MeV helium beam transported (deterministic) stops
  at ≈ 158 mm (the depth-dose R80 within a bin of the table range), the Bragg-curve
  target.
- **cross-backend** — reference vs Warp CPU (and CUDA) helium depth dose agree
  (deterministic total to the float32 budget), confirming the species-agnostic
  kernels transport helium consistently.

## Consequences

- Helium-4 is transported end-to-end (Bragg curve in water) by reusing the proton
  table, opening the multi-ion architecture with no hot-path change.
- **Deferred:** nuclear fragmentation and the fragment dose tail beyond the Bragg
  peak, ion-specific nonelastic cross-sections, carbon/oxygen (higher z, stronger
  fragmentation), species-resolved scoring, and the low-energy shell/Barkas/Bloch
  corrections and `z_eff(β)` beyond the constant `z = 2`.
