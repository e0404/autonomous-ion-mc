# Condensed-history methods for fast ion-therapy Monte Carlo (research note, 2026-09-30)

Source: physics-researcher subagent. Verified: gPMC (Jia 2012, PMC4474737),
goCMC (Qin 2017, PMC5730973), Holmes review (arXiv:2404.14543), Gottschalk
scattering power (arXiv:0908.1413), Kanematsu (arXiv:0806.0106, 0810.1390),
FRED carbon (De Simoni 2022, Front Oncol). Other items **unverified**.

## Stopping and steps
- gPMC: Bethe-Bloch restricted stopping power in water (Te,min 100 keV), other
  materials via tabulated ratio to water; step = min(voxel, hard event, 2 mm),
  ΔE/E ≤ 25%; mean loss via range-table inversion; cutoff 0.5 MeV local.
- goCMC: same scheme for carbon; **tuned I_water = 69.2 eV to match Geant4**
  (a calibration to another code, not physics).
- Review: min(d_vox, d_hard, d_max) universal; deltas deposited locally.
- Recommendation: CSDA range-table inversion for mean loss; tabulated
  unrestricted stopping powers (deltas not transported ⇒ unrestricted);
  defaults voxel boundary, d_max = 1 mm, ΔE/E ≤ 0.10, configurable; cutoff
  0.5 MeV/u with residual energy deposited locally and counted.

## Straggling
- gPMC/goCMC: Gaussian Bohr. FRED: Gaussian/Vavilov/Landau by κ (unverified).
- κ ≈ 0.06 for a 200 MeV proton in 1 mm water: per-step losses are not Gaussian,
  but variances add and the depth-dose is governed by the total variance.
- Recommendation: σ² = 2π r_e² m_e c² n_e z_eff² T_up/β² (1 − β²/2) s;
  Gaussian when κ ≥ 10 else moment-matched Gamma; clamp ΔE ≤ E and count clamps.

## Multiple Coulomb scattering
- Highland applied per step with quadrature addition is step-size dependent
  (Urban Geant4 test: −1.1% → +10.9% as steps shrink 1000×); goCMC re-fitted
  E_s = 19.8 + 0.0023 E MeV to compensate.
- Remedies: corrected Rossi (no log term), generalized Highland on the
  accumulated path integral, Gottschalk differential Molière
  T_dM = f_dM (E_s/pv)²/X_S, E_s = 15.0 MeV,
  f_dM = 0.5244 + 0.1975 lg(1 − (pv/p₁v₁)²) + 0.2320 lg(pv) − 0.0098 lg(pv) lg(1 − (pv/p₁v₁)²).
- Ions: θ ∝ z/(pv); Kanematsu agrees within ~2% for p/He/C in water; FRED
  carbon uses an empirical Highland scale factor 1.29–1.43.
- Water tests barely discriminate MCS models; use high-Z/mixed slabs
  (Gottschalk 1993) for discrimination.

## Effective charge / corrections (unverified magnitudes)
C and O fully stripped above a few MeV/u; charge pickup affects only the last
micrometres; I = 75 → 78 eV shifts proton range ≈1%; shell ≈1% below 10 MeV/u;
Barkas/Bloch 0.1–1% for carbon; density effect ≲0.3% to 430 MeV/u.

## LET estimators
Cortés-Giraldo & Carabe PMB 60 (2015) 2645; Granville & Sawakuchi PMB 60 (2015)
N283. Recommended: LET_d = Σ S(Ē)² ℓ / Σ S(Ē) ℓ, LET_t = Σ S(Ē) ℓ / Σ ℓ with
unrestricted electronic stopping power in water, transported charged hadrons
only, species-resolved numerators/denominators, NaN where undefined.

## GPU patterns
gPMC: particle-per-thread batches, secondary stack, XORWOW, float atomics,
batch uncertainty. goCMC: float32 with per-batch flush to reduce accumulation
error. Recommendation: generation-based kernels with compacted secondary
queues; float32 per-batch scoring with float64 accumulation; ≥10 batches;
counter-based RNG keyed by (seed, beamlet, history); dense per-beamlet scratch
grids compressed to sparse Dij (100 beamlets × 60×60×150 float32 ≈ 216 MB).

## Falsification probes suggested
Step refinement (ΔE/E 0.25→0.01, d_max 2→0.05 mm); high-Z slab MCS; LET grid
shift/step invariance; float32 vs float64 at high history counts and far-off
sources; batch count and seeds; I-value sensitivity; clamp/truncation counters.
