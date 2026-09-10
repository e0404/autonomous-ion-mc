# 0004 — Internal units and coordinate conventions

- Status: accepted
- Date: 2026-09-10
- Task: DEV-002
- Affects: numerical accuracy (unit-conversion errors), maintainability, API clarity, reproducibility

## Problem

`REQUIREMENTS.md` requires explicit physical units at every API boundary, an
explicit coordinate convention, and a consistent internal unit system. The
first physics code (stopping power and ranges) needs these fixed before any
transport code builds on them.

## Context

- Treatment-planning inputs (CT grids, dose grids, spot positions) are
  conventionally in millimetres; ICRU stopping-power tables are in
  MeV cm²/g and ranges in g/cm²; mass densities in g/cm³; mean excitation
  energies in eV.
- Warp kernels run in float32 (decision 0005); a unit system that keeps
  typical magnitudes near unity (energies in MeV, lengths in mm, densities in
  g/cm³) avoids needless dynamic-range loss.
- Ions are expressed either as total kinetic energy or per nucleon
  (MeV/u); mixing them is a classic error source.

## Candidate approaches

1. Pure SI (J, m, kg/m³). Rejected: every input and reference table needs
   conversion, and float32 magnitudes become awkward (energies ~1e-11 J).
2. Pure CGS with MeV (cm, g/cm³). Consistent with the tables, but conflicts
   with the mm convention of every planning system and CT.
3. **MeV, mm, g/cm³ with explicit helpers for the CGS table units
   (selected).**

## Selected approach

Internal units (documented in `src/ionmc/units.py`):

| quantity | unit |
|---|---|
| kinetic energy | MeV (total, not per nucleon; helpers convert) |
| rest energy | MeV |
| length, position | mm |
| mass density | g/cm³ |
| mass stopping power | MeV cm²/g (as in ICRU 49 / PSTAR) |
| linear stopping power | MeV/mm |
| mass (CSDA) range | g/cm² |
| linear range | mm |
| mean excitation energy | eV in material definitions; converted to MeV inside physics functions |
| angle | rad |
| charge | elementary charge |

Variables and parameters not in the internal unit carry a suffix
(`_cm`, `_ev`, `_g_per_cm2`). Physics functions receive `I` in MeV so that all
energies inside a formula share one unit.

Coordinates: right-handed Cartesian in mm; geometry objects define their
origin; beam and patient frames are related by explicit rigid transforms
introduced with the geometry stage. No axis is privileged by the physics
layer.

Physical constants: CODATA 2018 values in `src/ionmc/constants.py`. The
Bethe prefactor `K` is computed from the constants (0.3070749… MeV cm²/mol)
rather than hard-coded to the PDG rounding 0.307075; the difference is
3 × 10⁻⁷ relative and invisible at the 4-significant-figure level of the
reference tables.

## Rationale

The chosen units minimise conversions at the two boundaries that matter most
(planning geometry in mm; stopping-power and range tables in CGS-based ICRU
units) and keep float32 magnitudes benign.

## Expected tradeoffs

- Linear stopping power in MeV/mm and mass stopping power in MeV cm²/g
  differ by the density and a factor 10; the helper functions in
  `ionmc.units` are the only sanctioned conversion path.
- Total kinetic energy (not per nucleon) as the primary variable means ion
  APIs must convert explicitly; helpers are provided.

## Validation strategy

Unit tests of the conversion helpers and of the constant `K` against the PDG
value; the stopping-power validation against PSTAR (decision 0006)
exercises the MeV/eV/cm conventions end to end.

## Later validation outcome

Recorded with task DEV-002's validation record.
