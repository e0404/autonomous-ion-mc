# 0028 — Carbon-12 and oxygen-16 primary transport and the z² scaling accuracy

- Status: accepted
- Date: 2026-09-11
- Task: DEV-023
- Affects: multi-ion architecture, validation strategy

## Problem

Decision `0027` added helium-4 via equal-velocity z² scaling of the proton
stopping table (`scale_ion_stopping_table`, which is species-generic). The
**multi-ion architecture** (SHOULD) names carbon and oxygen, and milestone **V5**
requires carbon Bragg curves. This decision generalises primary-ion transport to
**carbon-12** (z=6, A=12) and **oxygen-16** (z=8, A=16), and — because the pure z²
scaling omits the higher-order z³ (Barkas) and z⁴ (Bloch) terms that grow with the
projectile charge — **quantifies the scaling's accuracy** against the independent
analytic Bethe model (which includes those corrections).

## Decision

1. **No new transport code.** `scale_ion_stopping_table(proton_table, particle)`
   (decision `0027`) already builds the ion table for any `Particle` — `E_ion =
   E_p·(m_ion/m_p)`, `S_ion = z²·S_p` — and the kernels are species-agnostic
   (`ALPHA`, `CARBON_12`, `OXYGEN_16` are defined in `ionmc.particles`). Carbon and
   oxygen are transported exactly as helium: build the scaled table, run with the
   corresponding particle. This closes the *Multi-ion architecture* breadth: all
   three clinical ions (He, C, O) are transportable.

2. **z² scaling accuracy is Z-dependent and bounded.** Compared with the
   independent Bethe model (`AnalyticStoppingPower(WATER, ion)`), the pure z²
   scaling agrees to (worst case, at the 10 MeV/u end-of-range end where the
   higher-order z terms are largest): **helium ~0.3 %, carbon ~1.1 %, oxygen
   ~2.0 %**; across the therapeutic plateau (50–400 MeV/u) it is < 0.5 % for all
   three. The error grows with z (the omitted Barkas ∝ z³ and Bloch ∝ z⁴ terms),
   so the pre-registered Bethe cross-check tolerances are ion-specific: helium
   ≤ 1.5 %, carbon ≤ 2 %, oxygen ≤ 3 %. This documents the model's fidelity limit;
   a per-ion table (ASTAR/ICRU 73) or an explicit Barkas/Bloch term would tighten
   it and is a later refinement.

3. **Fragmentation is a critical, explicit caveat for carbon/oxygen.** Unlike
   helium (fragmentation ~2–3 % of dose), carbon and oxygen fragment strongly; the
   **fragment dose tail beyond the Bragg peak is a defining clinical feature** and
   is **not** modelled here. The primary-ion Bragg curve produced by this task is
   therefore the primary-particle CSDA curve *without* the fragment tail, valid for
   the range/peak position but not the distal tail. Fragmentation is the next
   Stage-5 task and the remaining V5 gate.

## Validation (`validation/v5_carbon.py`, pre-registered per decision 0001)

- **stopping_vs_bethe** — the scaled carbon (and oxygen) table matches the
  independent Bethe model within the ion-specific tolerance (carbon ≤ 2 %, oxygen
  ≤ 3 %) over E/A = 10–400 MeV/u, and the Z-dependent worst-case error is recorded.
- **range_identity** — the carbon/oxygen CSDA range obeys
  `R_ion(E) = (m_ion/(z²·m_p))·R_p(E·m_p/m_ion)` to round-off (exact by
  construction).
- **carbon_bragg_range** — a 290 MeV/u carbon beam (3480 MeV) transported
  (deterministic) stops at ≈ 163 mm (R_C(290 MeV/u) ≈ 16.3 g/cm²), the clinical
  carbon range, conserving energy.
- **cross-backend** — reference vs Warp CPU (and CUDA) carbon depth dose agree to
  the float32 budget.

## Consequences

- Carbon-12 and oxygen-16 are transportable end-to-end (primary Bragg curve), and
  the z² scaling's Z-dependent accuracy is quantified against physics.
- **Deferred:** nuclear fragmentation and the fragment dose tail (the next Stage-5
  task, essential for a clinically faithful carbon Bragg curve), ion-specific
  nonelastic cross-sections, species-resolved scoring, per-ion ASTAR/ICRU tables,
  and explicit Barkas/Bloch corrections.
