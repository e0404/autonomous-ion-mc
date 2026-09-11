# 0029 — Bounded carbon-12 nuclear fragmentation and the distal dose tail

- Status: accepted
- Date: 2026-09-11
- Task: DEV-024
- Affects: multi-ion physics, validation strategy (V5 fragment-tail gate)

## Problem

Carbon and oxygen ions **fragment strongly**: ~50 % of a 290 MeV/u carbon beam
undergoes a nonelastic nuclear reaction before the Bragg peak (Haettner 2013), and
the lighter charged fragments (H, He, Li–B) — which at the same velocity have a
**longer range** than the primary carbon — deposit a characteristic **dose tail
beyond the Bragg peak**. This distal tail is the defining clinical feature that
distinguishes an ion depth dose from a proton's, and reproducing it is the last
**V5** milestone gate. DEV-023 transports carbon *primaries* only (no tail); this
decision adds a bounded fragmentation model that produces the tail with the correct
approximate magnitude and reach.

## Decision

1. **Carbon nonelastic attenuation — constant cross-section.** The ¹²C total
   reaction cross-section in water is ~energy-independent over the therapeutic
   range: `σ_R ≈ 1.4 barn per water molecule` (Zeitlin et al.; the ~1.2 b
   charge-changing value plus isotope-changing channels). With the water molecular
   density `n_mol = ρ·N_A/M_water = 3.343e22/cm³` this gives a macroscopic reaction
   rate `Σ = n_mol·σ_R ≈ 0.047 /cm` (mean free path ~21 cm), so the primary carbon
   **survival** is `S(z) = exp(−Σ·z)`. This reproduces the published fragmenting
   fractions (`1−S(R)`): ~53 % at 290 MeV/u (16.3 cm), ~72 % at 400 MeV/u (27 cm,
   the measured ~70 %).

2. **Three representative charged fragment species.** Per reaction the model emits
   **proton** (reuse `PROTON`; lumps H = p,d,t), **alpha** (reuse `ALPHA`; lumps
   He), and **boron-11** (new `BORON_11`, z=5, A=11; lumps the Z=3–5 heavies), with
   default multiplicities `N_H=2.0, N_He=0.7, N_B=0.35` per reaction (grounded in
   the measured partial charge-changing cross-sections; H/He set the far tail,
   boron the near-peak shoulder). These are the calibration handle for the
   tail-to-peak magnitude.

3. **Same-velocity forward emission.** In the peripheral (abrasion-ablation)
   picture a projectile fragment is a spectator continuing at the **beam velocity**,
   so at a reaction vertex with residual carbon energy `E_C` a fragment of mass `A_f`
   takes `E_f = A_f·(E_C/12)` (same energy per nucleon) emitted **forward (+z)**.
   The Goldhaber momentum spread gives only a few-percent range spread (below the
   straggling width), and emission angles are small, so forward-only is the right
   default for the 1-D depth-dose tail (a Gaussian angular halo is a 3-D later
   refinement).

4. **Deterministic orchestration reusing the multi-ion transport.**
   `carbon_fragmentation_depth_dose` (a standalone orchestration, no transport-kernel
   change): transport the carbon primary (nuclear off) for `D₀(z)`; apply the
   analytic survival to get the attenuated primary `D_p(z) = D₀(z)·S(z)`; the
   reactions in each depth bin (`N·[S(edge_lo)−S(edge_hi)]`) emit, per species, a
   forward fragment at `E_f = A_f·E_C(z)/12` (with `E_C(z)` from inverting the
   carbon CSDA range), transported by that species' z²-scaled table (decisions
   `0027`/`0028`) into the depth-dose grid. The total is
   `D = D_p + Σ_species D_fragment`. Because the primary and every fragment species
   transport through the *existing* validated CSDA drivers (reference and Warp
   CPU/CUDA), and the fragmentation itself is deterministic numpy, the model runs
   on both backends and inherits their parity. Energy is booked as
   `escaped = E₀ − D.sum() ≥ 0` (the transported fragments carry
   `Σ N_s·A_s/12 ≈ 0.72` of each reaction's energy; the remainder — target
   fragments, neutrons, binding, transverse momentum, and fragments leaving the
   grid — escapes).

## Validation (`validation/v5_fragmentation.py`, pre-registered per decision 0001)

- **fragment_tail_present** — for a 290 MeV/u carbon beam the **integrated
  distal-dose fraction** (the fraction of the total deposited dose landing more than
  a fixed margin distal to the peak) is **10–25 %** at +10 mm (canonical ~16 %),
  decreasing with margin, and is ~0 for a primary-only (no-fragmentation) run. This
  integrated metric is the gate because it is **resolution-robust** (invariant to the
  depth-bin width to < 0.01 across 0.5–4 mm bins); the single-bin *tail-to-peak*
  point ratio (~15 % only near 2–4 mm bins) is a diagnostic, since the model omits
  the beam energy spread that sets real peak sharpness.
- **tail_reach** — fragment dose stays `> 1 %` of the peak out to `≥ 1.5×` the
  carbon range, with non-zero dose extending toward `~2.5–3×` (H/He fragments).
- **primary_survival** — the surviving primary fraction at the peak matches
  `S(R) ≈ 0.50` at 290 MeV/u (and ~0.28 fragmented-away consistent with the
  cross-section).
- **energy_conservation** — reconciled independently rather than tautologically:
  the injected fragment KE reconstructed from the reaction weights, multiplicities
  and residual carbon energy matches the transported value (rel < 1e-6) and is
  deposited in full (rel < 1e-4, the grid contains the fragments); only then does
  `D.sum() + escaped = E₀` close, with `escaped ≥ 0`.
- **cross-backend** — reference vs Warp CPU (and CUDA) total depth dose (primary +
  fragments) agree to the float32 budget.

## Consequences

- Carbon depth dose now shows the characteristic fragment tail with the right
  approximate magnitude and reach, closing the V5 fragment-tail gate.
- **Deferred:** energy/angular-resolved fragment spectra, secondary fragmentation,
  neutron/gamma production, the full isotopic cocktail (d,t,³He,individual Li/Be/B),
  target fragmentation, species-resolved LET of the tail, the lateral fragment halo,
  and a discriminating TOPAS/Geant4 fragment-resolved reference (a `reference_
  requests/` item to lock the multiplicity/normalization).
