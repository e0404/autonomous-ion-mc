# 0023 — Dose-averaged LET (LET_d) scoring

- Status: accepted
- Date: 2026-09-11
- Task: DEV-018
- Affects: scoring, treatment-planning capability, validation strategy

## Problem

Ion-therapy planning and biological modelling need **LET** (linear energy
transfer), and `REQUIREMENTS.md` makes **dose-averaged LET** a MUST (*LET
scoring*: "at minimum, dose-averaged LET should be supported"). It is also the
next named validation item of milestone **V4**: "LET estimators versus published
proton LET-in-water data computed with the *same* estimator convention … and
versus analytical limits." This builds on the lab-frame 3-D dose scorer
(decision `0021`): LET_d is scored on a grid co-registered with `DoseGrid3D`.

## Decision

1. **Estimator — unrestricted dose-averaged electronic LET, "Method C".** Follow
   Cortés-Giraldo & Carabe (2015, PMB 60:2645) Method C. For each transport step
   `i` with start energy `E_i`, mean (unstraggled) midpoint energy loss
   `Δe_i`, and deposited energy `ε_i`:

   ```
              Σ_i  ε_i · L_i
     LET_d =  ───────────────        (den > 0)
                Σ_i  ε_i
   ```

   - `L_i = linear_stopping_power(E_mid_i, ρ_i)` — the tabulated **unrestricted
     electronic** linear stopping power (`ionmc.physics.transport.
     linear_stopping_power` → `tabulated_mass_stopping_power`) at the step-mean
     energy `E_mid_i = E_i − ½·Δe_i` (the *same* midpoint energy the RK2 step
     `midpoint_energy_loss` already forms), in **MeV/mm**. The evaluation energy
     is clamped to the table floor `table_e[0]` so the terminal deposit near the
     stopping point never extrapolates below the table.
   - `ε_i = w·de` — **bit-identical** to the energy already deposited into the
     dose grid at the step's lab midpoint (decision `0021`); weighting by the
     actual (possibly straggled) `de` keeps `den = Σ ε_i` exactly equal to the
     scored dose energy, while `L_i` stays on the physical `S(E)` curve at the
     unstraggled `E_mid`.

   **Why S(E), not `ε/Δl`.** The "deposited-energy-over-length" estimator
   (Method A, `L_i = ε_i/Δl_i`) is contaminated by spurious high-LET components
   wherever a step is truncated — and the grid path **clips every step to the
   nearest voxel face** (`_face_dist`, decision `0020`) and applies Gaussian
   straggling to `de`. Cortés-Giraldo & Carabe show Method A is not
   voxel-size-convergent (≈1.8× too high at the plateau at 0.2 mm vs 2.0 mm),
   whereas Method C (tabulated `S(E_mid)`) is stable against simulation
   parameters. Method C is immune to face-clipping and straggling by
   construction.

   **Unrestricted is the consistent choice.** The engine transports no delta
   rays (CSDA + local deposition), so there is no escaping energy to exclude;
   unrestricted `S_el` is both the community standard for "LET" comparisons and
   the only choice consistent with the transport model.

2. **Units — keV/µm, LET to medium.** Report LET_d in **keV/µm**. Because the
   engine's linear stopping power is already MeV/mm and **1 MeV/mm ≡ 1 keV/µm**
   exactly, the accumulated MeV/mm value *is* the LET in keV/µm — no conversion
   factor. `L_i` uses the local material stopping power (LET **to medium**); in
   water / homogeneous validation this equals LET to water. A water-normalised
   "LET to water" variant for heterogeneous media is a future option.

3. **Two co-registered accumulators, ratio on read.** LET scoring carries a
   second grid `let_num` (`Σ ε_i·L_i`, MeV/mm) parallel to the dose energy grid
   `dose` (`Σ ε_i`, the denominator). Both are deposited at the *same* lab
   midpoint as the dose, once per step, plus the terminal energy-cut deposit
   (`w·e` at `L` evaluated at the clamped residual energy — a small distal
   contributor that carries high LET, so including it preserves the distal-LET
   rise). `LET_d = let_num / dose` where `dose > 0`, with a **low-dose reporting
   mask** (`DoseGrid3D.let_d_kev_um(let_num, dose, min_dose_frac)` returns 0
   where the voxel dose is below a small fraction of the peak, since the distal
   falloff — highest LET_d — is exactly where dose→0 and the ratio is noisiest).

4. **Backend parity.** `let_num` and `dose` are two additive, strictly positive
   accumulations (no cancellation), so they inherit the dose scorer's
   cross-backend methodology (decision `0001`) and face-flip decorrelation
   caveat (decision `0020`): deterministic (scattering-off) per-voxel parity is
   the tight discriminating gate; under scattering the per-voxel ratio
   decorrelates as dose does and only aggregate/deterministic checks are tight.
   Because LET_d is a *ratio* of two correlated grids, correlated face-flip
   perturbations partially cancel, so its parity is at least as good as dose.

## Validation (`validation/v4_let.py`, pre-registered per decision 0001)

- **thin_voxel_analytic** — a monoenergetic proton crossing a thin voxel gives
  LET_d → `linear_stopping_power(E, ρ)` at that energy; check against the
  committed NIST PSTAR water values (keV/µm) at several energies (1, 10, 100,
  150, 250 MeV) to the interpolation tolerance.
- **distal_let_peak** — for a therapeutic pencil beam the LET_d rises with depth
  and its maximum lies **distal to** the Bragg dose maximum; entrance LET_d for
  ~150 MeV ≈ 0.4–0.5 keV/µm, distal region a few keV/µm.
- **voxel_size_convergence** — at a plateau depth (slowly varying energy) the
  per-voxel LET_d is stable across dose-voxel z-sizes (0.5/1/2 mm agree to a few
  percent), demonstrating the voxel-size convergence that motivates Method C over
  the face-clip-sensitive `ε/Δl` estimator.
- **dose_weight_consistency** — the LET_d denominator grid equals the scored
  dose energy to round-off (same ε_i).
- **cross-backend** — deterministic reference-vs-Warp-CPU LET_d per voxel to the
  float32 budget, and CPU-vs-CUDA deterministic per voxel; stochastic per-voxel
  ratio recorded as a diagnostic (decision `0022`'s corrected cross-backend
  story).

## Consequences

- A per-voxel dose-averaged LET grid is scored alongside dose on the reference
  and Warp CPU/CUDA paths, validated against published proton LET-in-water.
- **Deferred:** track-averaged / fluence-averaged LET_t (a third `Σ ℓ_i·L_i`
  accumulator; SHOULD, not MUST), restricted LET with a delta-ray cutoff
  (needs delta-ray transport), species-resolved LET (needs secondary/heavy-ion
  transport), energy-resolved LET spectra and dose-mean lineal energy, and a
  water-normalised LET-to-water variant for heterogeneous media.
