# 0006 — Analytical electronic stopping-power model, I-value policy and V0 acceptance criteria

- Status: accepted
- Date: 2026-09-10
- Task: DEV-002
- Affects: physical accuracy, validation strategy, scientific interpretation, reproducibility

## Problem

`REQUIREMENTS.md` requires an analytical physics layer alongside tabulated
data. The first quantity is the electronic stopping power of protons (and
later ions) in water and tissue-like materials. The model form, its
corrections, the mean excitation energy policy, the validity range and the
acceptance criteria against reference data must be fixed **before** the
comparison is run (review finding B2 on decision 0002).

## Context and evidence

- Reference data: NIST PSTAR (ICRU 49, I(water) = 75 eV) is *theory*
  (Bethe with shell, Barkas, Bloch and density corrections) above 0.5 MeV,
  fitted to experiment below; quoted uncertainty 1–2 % (elements) and
  1–4 % (compounds) at high energy (https://physics.nist.gov/PhysRefData/Star/Text/programs.html).
  NIST SRD tables are copyrighted; MCsquare redistributes PSTAR-derived
  per-material tables under Apache-2.0 (https://github.com/e0404/MCsquare,
  commit `211eefe6eaf2b8572d196d17f546f35ffb0ae0cf`). Two independent
  fetches of `Materials/Water/PSTAR_Stop_Pow.dat` agreed on all 20 native
  grid values now kept in `src/ionmc/reference_data/pstar_water.py`.
- Nuclear stopping of protons in water is < 10⁻⁴ of the electronic one
  above 1 MeV, so the table's totals equal electronic values at four
  significant figures.
- Density effect: Sternheimer–Berger–Seltzer (1984) parameters for liquid
  water as tabulated in Geant4 `G4DensityEffectData.cc` (`G4_WATER`:
  Cbar 3.5017, x0 0.2400, x1 2.8004, a 0.09116, m 3.4773, δ0 = 0; the often
  misquoted 0.097 is the fit error). For protons below ~900 MeV, δ = 0 in
  water.
- Shell correction: Barkas–Berger empirical formula
  `C(I, η) = A(η) I² + B(η) I³` with the coefficients of Geant4
  `G4IonisParamElm.cc` (0.422377, 0.0304043, −0.00038106 × 10⁻⁶;
  3.858019, −0.1667989, 0.00157955 × 10⁻⁹; I in eV, η = βγ), applied per
  element with ICRU 37 elemental I and combined as `Σ n_k C_k / N_e`
  (Geant4 `ShellCorrectionSTD`, `G4IonisParamMat.cc`); valid for
  η ≥ 0.13, frozen below the velocity of an 8 MeV proton and tapered
  logarithmically to zero at the velocity of a 2 MeV proton (Geant4).
- Barkas correction: Ashley–Ritchie–Brandt with the 47-point `F(W)` table,
  `b` per element and the 1.29 normalisation of Geant4
  `G4EmCorrections::BarkasCorrection` (cites Ashley & Ritchie, Phys. Rev. B
  5 (1972) 2393, and ICRU 49); `X = β²/(α² Z)`, `W = b/√X`,
  `L1 = 1.29 F(W)/(√(Z X) X)`, combined over elements with electron weights
  (Geant4 uses atom weights; the difference is < 0.01 % of S for water).
- Bloch correction: exact series `L2 = −y² Σ_j 1/(j(j² + y²))`, `y = zα/β`
  (Geant4 `BlochCorrection`), 16 terms plus integral tail.
- Mean excitation energy of water: 75.0 eV (ICRU 37/49, PSTAR) versus
  78 ± 2 eV (ICRU 90); PDG 2024 lists 79.7 eV. Δ ln I = 0.039 changes S by
  −0.4 to −1.1 % over 1–400 MeV.

## Candidate approaches

1. **Bethe + all four corrections with per-element shell/Barkas additivity
   (selected).**
2. Bethe + density effect only (as RayStation MC reportedly does):
   simpler; 0.6 % low at 10 MeV, 3 % at 2 MeV against PSTAR — inadequate
   for a stated 1 % criterion at 10 MeV.
3. Geant4-style split (Bragg parameterisation below 2 MeV, Bethe above):
   the low-energy parameterisation is itself a tabulated fit; deferred to the
   tabulated layer (DEV-003), which is the architecturally clean home for it.

## Selected approach

Approach 1, implemented as shared source in `src/ionmc/physics/stopping.py`
with the Python API `ionmc.stopping_power.AnalyticStoppingPower`. Each
correction is individually switchable (`Corrections`), and every result
carries provenance (I value and source, density-effect set, enabled
corrections, execution path).

**I-value policy.** `I` is a per-material, provenance-tagged parameter.
Liquid water defaults to 75.0 eV (`WATER = WATER_ICRU49`) so that the
reference comparison is I-consistent; `WATER_ICRU90` (78 eV) is provided and
its offset is reported separately. Density-effect parameter sets are bound
to the I value they were fitted with; a material with a different I drops
them (negligible below 900 MeV/u in condensed media).

**Validity claim.** Proton kinetic energies 10–1000 MeV (and the same
velocity range for ions). Below 10 MeV the analytic layer is degraded and
documented as such; the tabulated layer (DEV-003) covers it.

**CSDA range.** Integral of 1/S over ln E by composite Simpson (200
intervals, quadrature error < 10⁻⁶) from a 1 MeV floor; the neglected
residual range below 1 MeV is ≈ 0.0025 g/cm² for protons in water.

## Acceptance criteria (fixed before the comparison)

| comparison | criterion |
|---|---|
| analytic vs PSTAR table, E ≥ 10 MeV (17 points) | every point within **1.0 %** |
| analytic vs PSTAR table, 5 / 2 / 1 MeV | regression guards 2 % / 5 % / 10 % (not an accuracy claim) |
| CSDA range vs recalled PSTAR values 7.718 / 15.77 / 25.96 / 37.94 g/cm² (100–250 MeV) | 0.5 % guard; *informational* until the full table is acquired (DEV-003) |
| R(200) − R(100) vs recalled 18.242 g/cm² | 0.2 % |
| numpy binding vs Python binding | bitwise equal |
| float64 Python vs float32 Warp (CPU, CUDA) | S: rtol 1e-5; range: rtol 2e-5 (decision 0005) |
| Warp CPU vs CUDA (decision 0001 mixed criterion) | S kernel, class *transcendental*: rtol 4e-6, atol 1e-6; range kernel, class *iterative accumulation*: rtol 1e-5, atol 1e-6 |

## Rationale

The selected model is the same physics that generated the reference table
above its theory threshold, so agreement at the sub-percent level is the
correct expectation and any larger deviation is diagnostic of an
implementation error rather than a model limitation. Per-element shell and
Barkas additivity is the physically consistent form and keeps the shared
functions scalar (two material constants for the shell term, small arrays
for Barkas). Switchable corrections serve the configurable-fidelity
requirement and make each term's magnitude testable.

## Expected tradeoffs

- The shell-correction formula is empirical and freezes below 8 MeV; the
  analytic layer is not usable alone for the last ~1 mm of a proton track.
- Barkas needs a 47-point table and a per-element loop inside the shared
  function; cheap for table generation, irrelevant for transport, which will
  use precomputed tables.
- Water's I = 75 eV default follows the reference data rather than the
  latest recommendation; users must opt in to ICRU 90.

## Validation strategy

`tests/ionmc/test_stopping.py` (reference paths, every criterion above
except Warp), `tests/ionmc/test_warp_paths.py` (Warp paths, skipped without
warp) and `validation/v0_stopping_power.py` executed on the host runner for
the exact task SHA; results recorded below and in the DEV-002 validation
record.

## Validation outcome

*Reference paths (sandbox, Python 3.12, numpy 1.26; also Python 3.12.3 /
numpy 2.5.3 on the host):* maximum deviation from the PSTAR table for
E ≥ 10 MeV is **0.075 %** (all 17 points low by 0.01–0.07 %); 5 MeV +0.90 %,
2 MeV +3.46 %, 1 MeV +5.20 % (documented degradation). CSDA ranges
7.7169 / 15.7760 / 25.9623 / 37.9443 g/cm² versus recalled 7.718 / 15.77 /
25.96 / 37.94 (−0.015 %, +0.038 %, +0.009 %, +0.011 %); R(200) − R(100) =
18.2455 g/cm² versus 18.242 (+0.02 %). numpy and Python bindings bitwise
equal. Correction magnitudes at 10 MeV: shell −0.78 %, Barkas +0.30 %,
Bloch −0.05 %, density 0; ICRU 90 offset −1.14 % (1 MeV) to −0.43 %
(400 MeV).

*Warp paths (host runner):* see the DEV-002 validation record; to be
appended here.
