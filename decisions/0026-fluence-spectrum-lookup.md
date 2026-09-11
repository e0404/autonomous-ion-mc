# 0026 — Energy-resolved fluence-spectrum scoring and lookup-table accumulation

- Status: accepted
- Date: 2026-09-11
- Task: DEV-021
- Affects: scoring, treatment-planning/biological post-processing, validation

## Problem

The final **V4** milestone gate is "lookup-table accumulation reproduces offline
post-processing on scored spectra", and `REQUIREMENTS.md` lists **fluence scoring**
and **energy-resolved scoring** (SHOULDs). This needs (a) a scored proton
fluence-vs-energy spectrum and (b) an on-the-fly accumulation through a per-energy
lookup table that provably equals the same quantity computed by post-processing the
scored spectrum. It reuses the per-step quantities the dose (`0021`) and LET
(`0023`) scorers already form.

## Decision

1. **Track-length fluence spectrum (`FluenceSpectrum`).** The standard Monte-Carlo
   track-length estimator (Kellerer/Chilton; ICRU 85; Geant4 `G4PSCellFlux`,
   TOPAS `Fluence`): in a scoring region of volume `V`, `Φ = Σ_i w_i·ℓ_i / V`. For
   the energy-resolved spectrum, each step's weighted track length `w_i·ℓ_i`
   (`ℓ_i = s`, the step length in mm) is histogrammed into the bin of its
   **step-mean energy** `E_mid = E − ½·de` — the *same* energy the LET Method C
   scorer uses (`0023`). The scored array is the **raw** un-normalised histogram
   `counts[k] = Σ_{steps in k} w_i·ℓ_i` (mm); differential fluence
   `Φ(E)[cm⁻²·MeV⁻¹] = counts·100/(V_mm³·ΔE)` (100 = mm→cm on `ℓ` and mm³→cm³ on
   `V`), optionally `/N_histories` for per-primary fluence. Bins are linear,
   `n_bins = 160` over `0–160 MeV` (ΔE = 1 MeV > the per-step Δe ≈ 0.5 MeV at the
   `max_fraction=0.02`, `max_step_mm=1` defaults, so a step never straddles a bin
   and the midpoint binning is second-order and convergent as ΔE→0).

2. **Exact-identity lookup gate.** A per-energy lookup `w_tab[k]` (precomputed once
   at bin centres) drives an on-the-fly accumulator
   `A_gate = Σ_i (w_i·ℓ_i)·w_tab[k_i]`, `k_i = bin(E_mid_i)`. Because both the
   accumulator and the offline post-processing `A_offline = Σ_k counts[k]·w_tab[k]`
   read the **same** `w_tab` and the **same** bin rule, they are identical term for
   term — a **round-off-level identity** (float64 accumulators, `~4e-15` reference,
   within the float32 budget cross-backend), not a statistical one. This is the V4
   "reproduces offline post-processing" gate. Accumulators are float64 even in the
   float32 kernel (matching `0021`).

3. **Physical lookup `w = S_lin` reproduces deposited energy.** With
   `w_tab[k] = linear_stopping_power(center[k], ρ_ref)`, `A ≈ Σ_i ℓ_i·S_lin(E_mid_i)
   ≈ Σ_i de_i` = the total electronic energy deposited by the stepped transport, so
   the fluence spectrum folded with a stopping-power lookup recovers the total dose
   — a conservation cross-check (TOPAS's `EnergyFluence` pattern). Checked in
   **deterministic mode** (straggling/scattering off, so `de` is not stochastic);
   the terminal energy-cut residual is deposited as a point with no track step, so
   it is accounted separately (`total = A + Σ terminal residual`). The relative
   discretisation error is `~1e-3` (midpoint + bin quantisation) and shrinks with
   finer bins/steps — a convergence check.

4. **Backend-parity-friendly.** `counts`, `A_gate`, and the (optional) exact-E_mid
   `A_exact` are additive, strictly-positive float64 accumulations (no
   cancellation), so they inherit the dose scorer's cross-backend methodology
   (`0022`): the deterministic per-bin/scalar comparison is the tight gate. The bin
   rule is `k = floor((E_mid − E_lo)/ΔE)`, half-open bins, identical on both paths;
   out-of-range `E_mid` is dropped identically on both.

## Validation (`validation/v4_fluence.py`, pre-registered per decision 0001)

- **lookup_reproduces_offline** — `A_gate` equals the offline `Σ_k counts·w_tab`
  to round-off (`≤ 1e-12` reference), the V4 gate.
- **stopping_lookup_reproduces_energy** — with `w = S_lin`, `A_gate + terminal`
  equals the total deposited energy in deterministic mode to the binning tolerance
  (`~1e-3`), shrinking as `n_bins` grows (binning convergence).
- **fluence_sanity** — `counts ≥ 0`; the spectrum's support lies within
  `[E_cut, E_beam]`; the total track length `Σ counts` equals the summed step
  lengths.
- **cross-backend** — the deterministic per-bin `counts` and `A_gate` agree
  reference-vs-Warp-CPU and CUDA-vs-CPU to the float32 budget.

## Consequences

- A proton fluence-vs-energy spectrum is scored, and a lookup-table accumulator
  reproduces its offline post-processing exactly and recovers deposited energy —
  closing the last V4 milestone gate and satisfying fluence/energy-resolved
  scoring.
- **Deferred:** per-voxel/spatially-resolved spectra, angular/direction-resolved
  fluence, secondary-/heavier-species spectra, restricted (delta-cutoff) lookups,
  exact `1/S` intra-bin sub-splitting, log-spaced bins, and multiple simultaneous
  scoring regions.
