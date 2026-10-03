# 0038 — Electronic stopping power model, I-value and data roles

- Status: accepted
- Date: 2026-10-03
- Task: V3-002
- Affects: physical accuracy (range), reproducibility, validation strategy, data provenance

## Problem

Ion range in water is the single most consequential quantity of a therapy
Monte Carlo. The project must choose how electronic stopping powers are
obtained for protons, helium, carbon and oxygen in water and tissue-like
materials, which mean excitation energy of water is used, and which external
datasets play which evidential role, so that validation does not merely
reproduce its own inputs.

## Context

- `REQUIREMENTS.md` requires both an analytical physics layer and an
  external/tabulated data layer behind one interface, with versioned,
  checksummed, cached and provenance-tracked external data.
- Research report `docs/research/em-physics.md` (2026-10-03) established:
  - NIST PSTAR/ASTAR (SRD 124) can be downloaded as text by HTTP POST and
    carry an I_water of 75 eV (ICRU 49 basis); NIST SRD terms require
    attribution and the tables must not be redistributed in Git.
  - ICRU Report 90 recommends I_water = 78 ± 2 eV; its proton and alpha
    water stopping arrays are embedded in the Geant4 source file
    `G4ICRU90StoppingData.cc` (Geant4 Software License).
  - The 75 → 78 eV change shifts the 200 MeV proton range by about +1.1 mm
    (+0.4 %); ±2 eV corresponds to about ±0.3 % of range.
  - PSTAR, ASTAR, ICRU 73/90 and Bethe theory share one theoretical lineage
    above ~1 MeV/u, so agreement with them is ion-specific tabulated evidence,
    not measured evidence.
- The controlled host runner has no network; everything needed at validation
  time must be cached beforehand.

## Candidate approaches

1. **Tabulated layer only** (interpolate PSTAR/ASTAR; scale by z² for C/O).
   Rejected: no analytical layer, z²-scaling omits the Barkas and Bloch
   terms (1–2 % for C/O at tens of MeV/u), and the tables would be both
   construction and evaluation data.
2. **Analytical layer only** (Bethe-Bloch with corrections). Rejected: the
   requirements demand a tabulated data layer as well, and a single lineage
   gives no cross-check.
3. **Analytical default evaluated against independent tables, tabulated
   source selectable (selected).**

## Selected approach

### Model (analytical layer)

Mass electronic stopping power
`S/ρ = K z² (Z/A) β⁻² [L₀ + z L₁ + z² L₂]` with

- `K = 0.307075 MeV cm² mol⁻¹`, Bragg additivity for Z/A and ln I of
  compounds;
- `L₀ = ½ ln(2 mₑc² β²γ² T_max / I²) − β² − δ/2 − C/Z`, exact `T_max`
  including the mass ratio terms; Sternheimer density effect with
  per-material parameters (zero in water below 500 MeV/u, kept for
  generality); shell correction `C/Z` from the Bichsel/ICRU 49 parameterisation
  in η = βγ, evaluated per target element with the elemental I-values and
  combined as an electron-weighted sum (so overriding a compound's I does
  not alter its shell correction); below the validity limit η = 0.13 the
  correction is held at its boundary value, which affects only the last
  ≈ 0.1 mm of range;
- `L₁` Barkas term in the Ashley–Ritchie–Brandt form, implemented as the
  tabulated function and shell parameters used by Geant4
  (`G4EmCorrections::BarkasCorrection`, read from the v11.4.2 source),
  summed over the target elements;
- `L₂` Bloch term `−y² Σ 1/[n(n²+y²)]`, `y = zα/β`; Mott term for z ≥ 2 as
  in Geant4 (`½ π α β z`), both explicitly z-dependent so that He/C/O are not
  z²-scaled protons;
- effective charge for z ≥ 2 (Pierce–Blann form; its departure from z is
  confined to low energies);
- below `E_low = 1 MeV/u` no stopping power is evaluated: transport deposits
  the residual energy locally (documented domain boundary).

Default `I_water = 78 eV` (ICRU 90), configurable per material and recorded
in every result. Tables are built in float64 on a logarithmic energy grid
(per nucleon), stored per material: `S_el(E)`, CSDA range `R(E)` by
integration of `1/S`, and the inverse `E(R)`.

### Tabulated layer

Datasets are declared in a registry with id, version, URL (and POST body),
pinned SHA-256, licence and citation. `ionmc data fetch` downloads to a
user-controlled cache (`IONMC_CACHE_DIR`), verifies the hash, and stores a
manifest; offline reuse re-verifies the hash. Parsers turn NIST STAR text
into the same `StoppingTable` structure the analytical layer produces, so
transport is indifferent to the source (`stopping_source = "bethe" |
"nist-star"`).

### Legal assessment of the data sources

- NIST PSTAR/ASTAR: Standard Reference Data; use with attribution is
  permitted. The full tables are cached locally and not redistributed in
  Git; the test suite contains short excerpts (about a dozen rows) with the
  NIST copyright notice and citation, which the lead judges to be within
  the attribution-based use permitted by NIST's licence statement
  (https://www.nist.gov/open/license). The notices are collected in
  `THIRD_PARTY_NOTICES.md`.
- ICRU 90 water arrays: the numbers are obtained from a Geant4 source file
  distributed under the Geant4 Software License, which permits use and
  modification with attribution. The ICRU report itself is copyrighted; the
  project uses the values only as a cached evaluation reference, cites ICRU
  Report 90 and Geant4, and does not redistribute the ICRU report. The
  Barkas-correction table implemented in the analytic layer is transcribed
  from Geant4's `G4EmCorrections` and the test suite contains short
  excerpts of the ICRU 90 arrays; both are covered by the Geant4 Software
  License, whose text and required acknowledgment are included in
  `THIRD_PARTY_NOTICES.md`. This is the lead's assessment; it is recorded
  here because the research report flagged it as unverified.

### Evidence roles

| Dataset | Role | Claim it supports |
|---|---|---|
| NIST PSTAR water (SRD 124) | evaluation | Bethe implementation at I = 75 eV reproduces ICRU 49 stopping and CSDA range |
| NIST ASTAR water | evaluation | same for alpha particles |
| ICRU 90 water arrays (Geant4 source) | evaluation | default I = 78 eV tables |
| ICRU 73/90 ion tables (if acquired) or TOPAS runs | evaluation | carbon/oxygen stopping and range |

When `stopping_source = "nist-star"` is selected, PSTAR/ASTAR become
construction data for that run; the result metadata records this and such
runs cannot be evaluated against the same table.

## Expected tradeoffs

- The ICRU 90 I-value makes proton ranges ~0.4 % longer than PSTAR-based
  codes; this is deliberate and documented. Comparisons with MCsquare, FRED
  and TOPAS must state each engine's I-value.
- The analytical model is expected to agree with PSTAR/ASTAR within ≈ 1 %
  above 10 MeV/u and to degrade below; residual-range consequences are
  sub-voxel.
- No measured stopping data are used; measured evidence enters through
  depth-dose and range comparisons in later suites.

## Validation strategy

- Unit tests against transcribed anchor values (a handful of PSTAR/ASTAR
  points) at I = 75 eV.
- Validation script comparing full PSTAR/ASTAR tables (cached) with the
  analytical layer at 75 eV, and ICRU 90 arrays at 78 eV, over 2–500 MeV/u,
  reporting relative deviations and CSDA range differences.
- Later: carbon/oxygen range against TOPAS and tabulated ion data (V3-008).

## Later validation outcome

Comparison script `validation/scripts/stopping/compare_nist.py` executed at
the clean task SHA `d3945d1f54e26c24b82e8737109f2a5c2e2034b2` against the cached NIST PSTAR/ASTAR water tables
and the ICRU 90 water arrays. Every number below is copied from the committed
aggregate result file `validation/results/stopping/2026-10-03-compare-nist.json`
(dataset hashes inside; the file contains aggregates only, no reference
values or per-energy ratios).

| Comparison of electronic stopping power | Points (E ≥ 10 MeV/u) | Max abs. rel. deviation | RMS rel. deviation | Max abs. rel. deviation, E ≥ 2 MeV/u |
|---|---|---|---|---|
| proton, Bethe (I = 75 eV) vs NIST | 34 | 0.079 % | 0.028 % | 1.835 % |
| proton, Bethe (I = 78 eV) vs ICRU 90 | 14 | 0.077 % | 0.032 % | 1.871 % |
| proton, ICRU 90 vs NIST (reference tables themselves) | 25 | 0.701 % | 0.430 % | 1.009 % |
| alpha, Bethe (I = 75 eV) vs NIST | 35 | 0.203 % | 0.114 % | 0.711 % |
| alpha, Bethe (I = 78 eV) vs ICRU 90 | 13 | 0.200 % | 0.114 % | 0.712 % |
| alpha, ICRU 90 vs NIST (reference tables themselves) | 13 | 0.716 % | 0.558 % | 0.984 % |

CSDA range (water):

| Projectile | Energies (MeV/u) | Max abs. rel. deviation, Bethe 75 eV vs NIST | Max abs. rel. deviation, Bethe 78 eV vs integrated ICRU 90 | Bethe 78 − 75 eV shift (mm) per energy |
|---|---|---|---|---|
| proton | 100, 150, 200, 250 | 0.073 % | 0.175 % | +0.42, +0.81, +1.30, +1.85 |
| alpha | 100, 150, 200 | 0.101 % | 0.067 % | +0.42, +0.82, +1.30 |

Interpretation:

- Above 10 MeV/u the analytic layer reproduces the ICRU 49 (PSTAR/ASTAR) and
  ICRU 90 electronic stopping powers within the tabulated maxima (below
  0.1 % for protons and 0.21 % for alpha particles at the respective
  I-values) and CSDA ranges within 0.2 %. Between 2 and 10 MeV/u the
  deviation grows to about 2 % (protons) because the shell correction is
  held at its validity boundary; the residual-range consequence is below
  0.1 mm, as anticipated in the selected approach.
- The reference tables themselves differ: ICRU 90 stopping powers lie up to
  about 1 % below PSTAR/ASTAR over the therapeutic range (table row "ICRU 90
  vs NIST"), more than the pure I-value change the analytic model produces,
  so the analytic model at 78 eV reproduces ICRU 90 stopping to 0.1 % while
  its ranges exceed the integrated ICRU 90 ranges by up to 0.18 % at the
  highest proton energy. This is preserved as a known model-table
  difference, not tuned away.
- The ICRU 90 alpha arrays in the Geant4 source are indexed by total alpha
  kinetic energy; this reading makes the ICRU 90/ASTAR ratio consistent with
  the proton ratio and with the expected I-value effect.
- These are ion-specific tabulated comparisons sharing Bethe-theory lineage
  with the model; they establish implementation correctness, not measured
  physical accuracy. Measured and independent Monte Carlo range evidence is
  collected in later suites.
