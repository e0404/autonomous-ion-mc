# 0040 — Electronic stopping power: analytic layer, tabulated data and lineage

- Status: accepted
- Date: 2026-09-30
- Task: V2-003
- Affects: physical accuracy, validation strategy, dependency footprint, reproducibility

## Problem

Ion ranges in tissue are set by electronic stopping power at the 0.1–1 %
level. The protocol requires an analytic physics layer, an external tabulated
layer, transport independent of the source, automatic acquisition with
provenance and offline reuse, and ion-specific evidence per species. The
research record (`validation/research/2026-09-30-stopping-power-data-sources.md`)
found no licence-clean, machine-readable heavy-ion table except the public
NIST addendum implementing ICRU Report 90 (protons, alphas, carbon ions in
water, air, graphite), while Geant4's G4EMLOW ion tables carry a
"use within Geant4" restriction and were rejected.

## Evidence

- NIST PSTAR/ASTAR (ICRU 49, I(water) = 75 eV) are reachable through a form
  POST and parse into 7-column text tables; acquisitions
  `de5f4f25…` (PSTAR water) and `fed7209e…` (ASTAR water), role *evaluation*.
- The NIST ICRU 90 addendum PDF (sha256 `e8340e3d…`, acquisition
  `7e9d36b1…`) contains tables A.9/A.12/A.15 for water at I = 78 eV; the
  values were transcribed by text extraction and shipped as
  `ionmc/data/icru90_water.json` (8 kB) with the source hash.
- The corrected Bethe formula reproduces PSTAR/ASTAR at I = 75 eV within
  0.1 % (p, 8–500 MeV) / 0.25 % (He, 8–250 MeV/u) and ICRU 90 at I = 78 eV
  within 0.12 % (p) / 0.25 % (He) / 0.6 % (C, ≥ 16 MeV/u) / 0.3 % (C,
  ≥ 30 MeV/u); it is 1.5 % low for carbon at 10 MeV/u and 9 % low at
  2 MeV/u. Effective-charge scaling of the proton table reproduces the ICRU 90
  helium and carbon tables within ±2 % above 1 MeV/u.
- I = 75 → 78 eV changes proton stopping by −0.45 % at 100 MeV and range by
  +0.5 %, as expected.

## Decision

1. Mean excitation energies follow ICRU 90 (water 78 eV, air 85.7, graphite
   81) and ICRU 37/49 elsewhere; no I value is tuned to any Monte Carlo code.
2. Tables blend a tabulated low-energy source with the analytic layer above
   species-dependent windows (8–16 MeV/u for z ≤ 2, 16–32 MeV/u for z ≥ 3):
   ICRU 90 for water (p, He, C direct; other ions by effective-charge scaling
   of the proton table), PSTAR/ASTAR for materials with verified NIST codes,
   and the ICRU 90 water shape scaled by the analytic ratio for tissues
   without tables (ICRP muscle, soft tissue, lung; calcium is absent from
   PSTAR so bone uses its own NIST tissue table). The fallback is recorded in
   the provenance and affects only the last millimetre of range.
3. The analytic layer stays the source above the window even where an ICRU 90
   table exists, so that the ICRU 90 values above the window remain
   *held-out* evidence for the formalism (tested in
   `tests/ionmc/test_stopping_physics.py`).
4. Evidence lineage for PHY-STOPPING-RANGE: for water, the same-I reference
   is ICRU 90 (ion-specific tabulated, construction below the window,
   evaluation above it); PSTAR/ASTAR (ICRU 49) are a second, independent
   evaluation whose known −0.45 % stopping / +0.5 % range offset from the
   I value is documented rather than absorbed. For other materials
   PSTAR/ASTAR share the ICRU 37/49 I values and are the held-out reference.
   Measured ranges and TOPAS remain the downstream evidence for transport.

## Tradeoffs

- The Barkas and Mott terms are approximate; the residual −0.3…−0.45 %
  for carbon at 100–200 MeV/u is accepted within the frozen 3 % table
  tolerance and could be improved with a Lindhard–Sørensen implementation.
- Shipping the ICRU 90 water table in the package trades strict "no data in
  Git" for offline availability of the single most important dataset; it is
  8 kB of key data with provenance, analogous to physical constants.
- Materials without tables rely on a water-shape approximation below
  16 MeV/u.

## Validation

`tests/ionmc/test_stopping_physics.py` (analytic vs ICRU 90 above the window,
continuity, range integration, inverse tables, fallbacks) and
`tests/ionmc/test_nist_star.py::test_real_nist_fetch` (local, network). Ranges
of transported particles are validated in later tasks against these tables
(PHY-STOPPING-RANGE) and independent references.
