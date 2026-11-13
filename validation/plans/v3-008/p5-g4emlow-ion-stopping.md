# V3-008 pre-freeze check P5: G4EMLOW 8.8 ion stopping data (record, 2026-11-01)

Nothing is frozen by this note. No table values are committed; only hashes, spans and point counts.
Facts below are orchestrator-verified from the acquired archive.

## Acquisition

| Item | Value |
|---|---|
| Archive | `G4EMLOW.8.8.tar.gz` |
| sha256 | `b60cfd63176f5d16107e2a25b35b235155032d1735d749670ca50fede12624cf` |
| Size | 350,393,595 bytes |
| Acquisition id | `3ee0c1f90a224ebbbc5b058433911f92` |
| Role | exploratory |
| Location | task cache `.ionmc-cache/reference/b60cfd63…` |

Four members were extracted read-only into the scratchpad. Nothing is committed from them.

## Members (format: header "Emin Emax N", then N rows of energy per nucleon [MeV] and mass stopping power [MeV cm2/mg], per the README)

| Member (`G4EMLOW8.8/ion_stopping_data/...`) | sha256 | Points | Span [MeV/u] | Points >= 2 MeV/u | Points in [50, 500] |
|---|---|---|---|---|---|
| `icru73/z6_G4_WATER.dat` | `444be013f716f5701051a5d4acf1a5ee7ebbdff7b6e2ebc78ca07563f9820159` | 53 | 0.025-1000 | 32 | 12 |
| `icru73/z8_G4_WATER.dat` | `6d8255d47aca417bd10d19684350203c0276d3df828886fbdccaea8cdb19bd67` | 53 | 0.025-1000 | 32 | 12 |
| `icru90/z6_G4_WATER.dat` | `d9d4809dcbacdb74e5e6c069709ef848c9d715b097f4ff4db24ef9969e727d75` | 243 | 0.00025-1000 | 170 | 70 |
| `icru90/z8_G4_WATER.dat` | `863898af3649947df532c6eeb244b2658503c0c7fa3882f81a504278793717f6` | 234 | 0.00025-1000 | 161 | 70 |

`ion_stopping_data/README` sha256 begins `a4251df8f2d4af17…`.

## Outcome of the checks

- (i) Span: all four files cover 2-430 MeV/u (to 1000 MeV/u).
- (ii) Units and axis: energy per nucleon [MeV], mass stopping power [MeV cm2/mg] (README).
- (iii) Erratum: the README states that the icru73 water files are the "revised tables provided by
  P. Sigmund/Univ. of Southern Denmark (revised tables are: G4_WATER (all considered ions) ...)".
  The ICRU 73 water erratum is therefore resolved in the icru73 files, by documentation. The
  high-energy slope test against Bethe at two I values was not run.
- (iv) Origin of `icru90/z6` and `icru90/z8`: **OPEN.** The README does not describe the icru90
  directory. ICRU 90 itself published no oxygen table (research note), so `icru90/z8` is a
  Geant4-side construction of unknown origin.

**Terms.** The README states: "DATA ARE *NOT* FOR COMMERCIAL USE AND MUST BE USED WITHIN GEANT4".

**Consequence.** Under the stated terms the G4EMLOW tables may not be used outside Geant4. They
cannot be the project's independent evaluation source for row I3. They remain a LINEAGE record of
what TOPAS reads (below the 2.0 MeV scaled switch; see P2). The acquisition role stays exploratory.

## P5b (OPEN, not done here)

libdEdx (GPL-3; github APTG/libdedx master tarball): sha256
`b941a0c30793daddc4d119a2cc664bfbdd430554c3a4788be04e97ff6fcc0863`, 962,747 bytes, acquisition id
`f062ca75e28842e3b1c5ca0f1f0e2a72`. It embeds ICRU 73 (revised) heavy-ion tables. Its C-12 and O-16
liquid-water tables must be located and identified. Until that is done the I3 evidence source is
pending.

## Lineage checks to run at P5/P2 closure (acceptance plan row P5)

Because option 4 reads `icru90/z8_G4_WATER.dat` and `icru90/z6_G4_WATER.dat` below the 2.0 MeV scaled
switch and falls back silently when a file is absent: record existence, sha256, span/origin and a
silent-fallback check for both files (above).
