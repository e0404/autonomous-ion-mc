# P5b: libdEdx ICRU 73 (revised) tables as the I3 evaluation source (pre-freeze check, V3-008)

Status: pre-freeze record for check P5b of `validation/plans/v3-008-acceptance.md`, written 2026-11-01.
- The facts below were verified by the orchestrator from the materialized libdEdx archive. The planner
  recomputed the three member hashes from the archive.
- No stopping value is reproduced here, and none is committed (decision 0038 policy).
- **The values were NOT compared with the Geant4 G4EMLOW files.** The G4EMLOW README restricts their use
  to Geant4, so they are a lineage record only (P5).
- No comparison with the ionmc model has been made.

## Source

| Item | Value |
|---|---|
| Archive | libdEdx source with data files, sha256 `b941a0c30793daddc4d119a2cc664bfbdd430554c3a4788be04e97ff6fcc0863`, 962 747 B, acquisition `f062ca75e28842e3b1c5ca0f1f0e2a72` |
| Licence | GPL-3. The data are used as a cached evaluation reference outside Git, as decision 0038 does for NIST data |
| Citation | Lühr, Toftegaard, Kantemiris, Hansen, Bassler, Int. J. Radiat. Biol. 88 (2012) 209 (libdEdx; mentions a "new and old version" of ICRU 73) |

## Members (sha256 recomputed from the archive)

| Member | sha256 | Bytes | Content (FACT, orchestrator) |
|---|---|---|---|
| `data/raw/ICRU73_NEW.dat` | `749768c1c2de475c928cc2c50a354903577139c11debd4b5948b035ca3a3d87e` | 21 260 | 36 blocks with header `#target:ion:53` |
| `data/raw/icru73_newEng.dat` | `8b5b809f5c0cc4f6018f9faa89e9f73c7b9b81a3c588f98bef556712d4c77dd1` | 198 | 53 energies, 0.025–1000 MeV/u. The grid is identical to that of the Geant4 icru73 water files (grid identity only; no values compared) |
| `src/data/embedded/dedx_icru73new.h` | `f6d088b233e6584363ee09db6b6d828dde0a3c2fe62850e2a35285e83ff78459` | 35 877 | ions 1–18 (= Z); targets 276 and 104; energies in MeV/nucl; stopping in MeV cm²/g |

## Metadata (`dedx_metadata.h`, FACT, orchestrator)

- Target 276: density 1.000 g/cm³, I = 75.0 eV (liquid water, NIST material id 276).
- Target 104: density 1.20479e-3 g/cm³, I = 85.7 eV (dry air).
- These are libdEdx's material constants. Whether the ICRU73_NEW water values were *computed* with
  I = 75 eV or with the 78 eV of the 2009 revision is **TO-VERIFY**.

## Identification

- libdEdx README: `DEDX_ICRU73` = "ICRU Report 73 (2005)"; `DEDX_ICRU73_OLD` = "older parametrization".
- Whether ICRU73_NEW equals the 2009 Sigmund water revision ("Errata and Addenda for ICRU Report 73",
  J. ICRU 5(1), 2009) is **TO-VERIFY**. To check: the libdEdx documentation and commit history, and
  Lühr et al. 2012.

## Consequence for row I3 (as written in the plan)

- **What is compared:** the ionmc analytic model at I = 78 eV (not tuned) against libdEdx ICRU73_NEW,
  target 276, for Z = 6 and 8, on the 53-point grid above 2 MeV/u.
- **Tolerances:** unchanged (S: 1.5 % max / 0.8 % RMS above 10 MeV/u, 6 % for 2–10 MeV/u; CSDA range
  0.5 %).
- **I-value difference:** if the table was computed with 75 eV, the known, bounded contribution is
  about 0.4–0.6 % in S above 10 MeV/u and about 0.4–0.5 % in CSDA range (estimate from L0 ∝ −ln I).
  It is stated in the result. A range-clause failure within that bound is recorded as FAIL with the
  attribution and is never re-graded.
- **Role:** exploratory until I3 is frozen, evaluation afterwards. Only aggregates go into Git.
- **Independence:** TOPAS consults ICRU 73/90 ion data only at ≤ 2.0 MeV scaled energy (P2), so the
  comparison is largely independent of TOPAS above about 2 MeV/u. The shared lineage is Bethe theory and
  the I-value.

## OPEN

1. Equality of ICRU73_NEW water with the 2009 Sigmund revision.
2. The I-value used for the ICRU73_NEW water values (75 vs 78 eV).
3. Whether target 276 holds liquid water for both Z = 6 and Z = 8 across the full 53-point grid
   (block presence, checked by a loader test without printing values).
