# P6 — measured-data access (V3-008), research note 2026-11-01 (not evidence)

Labels: FACT = verified this session by fetching; LIT = from memory/secondary source; OPEN = not checked.
Stopped early on coordinator instruction (turn limit). Nothing written to any worktree.

## 1. Sihver, Schardt, Kanai, Jpn. J. Med. Phys. 18(1) (1998) 1-21

Access (FACT)
- DOI 10.11323/jjmp1992.18.1_1 -> https://www.jstage.jst.go.jp/article/jjmp1992/18/1/18_1/_article
  Note: J-STAGE journal code is `jjmp1992`, NOT `jjmp` (the `jjmp` series starts at vol. 20; that is why the
  earlier lookup failed). Found via CiNii Research (crid 1390001204672683904).
- Article page states "JOURNAL FREE ACCESS"; received 24 Jun 1997, accepted 1 Jun 1998, on J-STAGE 24 Sep 2012.
  No Creative Commons licence on the page (grep: 0 hits) -> free to read, copyright with the publisher
  (Japan Society of Medical Physics); redistribution of the PDF is NOT licensed; extracting a few numbers
  with citation for validation is the normal scholarly use (licence: free access, all rights reserved).
- PDF: https://www.jstage.jst.go.jp/article/jjmp1992/18/1/18_1/_pdf/-char/en
  Content-Disposition filename 18_1.pdf; 934 294 bytes; sha256
  68567feb58bf3c2a669193e87d1f0ce41f86c21dcd36959f5ffd7c390bfcda55; PDF 1.4, 21 pages.
  Scratch copy only: $TMPDIR/sihver/pdf/18_1.pdf (TMPDIR = /home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000).
- The PDF is a SCAN: 21 CCITT G4 page images, no text layer (FACT). Numbers must be read from images
  (transcription by eye, or OCR); any extraction tool chain must state this. Pages rendered with
  `mutool draw -r 110` to $TMPDIR/sihver/png/p01..p21.png.

Content (FACT for pp. 1-7, 10-11 inspected; rest OPEN)
- Measurements at GSI SIS, water target: 12C 195, 270, 330 MeV/u; **18O** 200, 300, 400 MeV/u; 20Ne 670 MeV/u.
  RIKEN: 12C and 20Ne at 135 MeV/u (depth-dose and dose-average LET). (abstract, p.1)
- **The oxygen beam is 18O, not 16O.** ionmc's O species is 16O. Using 18O needs either an 18O mass/species in
  ionmc (range scales ~A/Z^2 at equal MeV/u, i.e. 18O range ~ 9/8 of 16O plus small corrections) or a
  documented scaling — an approximation that must be declared. C (12C) is directly usable.
- Paper is mainly a 1-D model paper (Bethe/Salamon-Heinrich code, fragmentation model).
- Table 1 (p.6) = CALCULATED ranges (g/cm^2) for 270 MeV/u 12C, 400 MeV/u 18O, 670 MeV/u 20Ne from six codes
  (ATIMA, Steward, Salamon/Heinrich, TRIM92, Bichsel RST, Sihver-Kanai) with I = 75 / 67.5 / 80 eV.
  NOT measured data. Fig. 1 = calculated range vs I. Fig. 2/3 (pp.7,10) = charge-changing cross sections
  and attenuation (Schardt et al.), not ranges.
- Whether measured Bragg-peak positions/R80 with uncertainties are TABULATED (pp. 12-21): **OPEN** — pages
  rendered but not inspected before the stop. Indicator against tabulation (LIT/secondary, FACT that the
  sentence exists): Pshenichnov et al. arXiv:0704.2519 acknowledge "Dr. Dieter Schardt ... for providing us
  the data tables for the depth-dose distributions of carbon ions", suggesting the paper gives figures only.
  Next step (≤ 15 min): view $TMPDIR/sihver/png/p12..p21.png for tables of measured peak positions,
  uncertainty statements, and beam-line WET/offset statements (ripple filters, windows, ionisation chambers).

Lineage (FACT/LIT)
- Measured data at GSI/RIKEN, independent of Geant4/TOPAS and of ICRU 90 as data. LIT/FACT (secondary):
  Pshenichnov 2008 (arXiv:0704.2519) set I_water = 77 eV because it "better describes the set of available
  data on depth-dose distributions measured with carbon-ion beams (Sihver 1998; Haettner 2006)" -> these data
  are part of the I-value calibration lineage (consistent with plan V3 row wording). Whether ICRU 90 cites
  Sihver 1998 directly: OPEN.

## 2. Measured ion lateral widths in water (Parodi 2013, Schwaab 2011, Tessonnier 2017, others)

- NOT investigated in this session (stopped). Status remains as in research-data-sources.md:
  Tessonnier 2017 PMB 62:3958 paywalled, figures only (repo note, earlier FACT); GSI report behind a bot
  challenge; Haettner geant-val = fragment yields only (no lateral dose). Parodi 2013 J. Radiat. Res. (likely
  OA, Oxford/JRR CC BY-NC, LIT), Schwaab 2011, Kempe & Brahme, Haberer, HIMAC/CNAO/MedAustron
  tabulated sigma/FWHM: all OPEN.
- (3) p+water lateral data at matched velocity: OPEN (Gottschalk 2015 halo arXiv:1409.1938 is already
  known to the project — proton, 177 MeV, tabulated dose vs radius; usable only for a relative check).

## 3. Verdicts (provisional)

- V3 (measured C/O range): Sihver 1998 is ACCESSIBLE (free, 934 kB scanned PDF, URL above) — access part of P6
  PASSES. Usability NOT yet established: (a) tabulation of measured peak/R80 positions with uncertainties OPEN
  (pp. 12-21); (b) O is 18O (species mismatch; C usable as is, O only via an added 18O species or a declared
  scaling); (c) scanned PDF -> manual transcription with provenance. If pp. 12-21 show only figures: V3 =
  EVIDENCE GAP (figure digitisation would be a weaker, declared alternative; decide in the plan), He: no
  candidate (gap).
- Lateral (V5-SS / declared gap): no data found in this session -> EVIDENCE GAP not yet established by search
  (search not performed); per decision 0043 §3 the intervention applies only after P6 and P10 fail, so the
  lateral part of P6 must be completed (bounded search) before invoking it.

## Orchestrator addendum (2026-11-01): pages 12–15 of Sihver, Schardt, Kanai 1998 read

FACT (read from the rendered pages p12–p15 of the J-STAGE PDF, sha256 68567feb58bf3c2a669193e87d1f0ce41f86c21dcd36959f5ffd7c390bfcda55; scanned, no text layer): section III describes the GSI depth-dose measurements (195, 270, 330 MeV/u 12C; 200, 300, 400 MeV/u 18O; 670 MeV/u 20Ne; water column 2–30 cm, entrance wall reproducible to better than 50 µm) and the RIKEN measurements (135 MeV/u 12C and 20Ne). Section IV states: "In Fig. 5-12, the measured curves were always normalized in the horizontal direction so that the measured Bragg peaks coincide with the calculated ones." Figures 5–8 show the measured and calculated Bragg curves as plots only; no table of measured peak positions or ranges with uncertainties appears on pages 12–15, and the horizontal normalisation removes the absolute range information from the figures. The oxygen beam is 18O, not 16O.

**V3 verdict: EVIDENCE GAP for the measured-range row.** This candidate cannot provide absolute measured ranges in water for C-12 (or O-16); it remains usable only as a shape (peak-width, fragment-tail) illustration, informative at most. Row V3 therefore goes to the intervention path of the plan unless another public measured-range source is found (none in hand: Tessonnier 2017 paywalled, GSI reports behind a bot challenge, Haettner = fragment yields). The lateral-data search (P6 part 2) is still OPEN; the intervention for the lateral gap is used only after P6 and P10 close.

## Orchestrator/planner addendum (2026-11-01, r6, Codex REVIEW-b469dff8 finding 2): lateral search outcome

The text above calls the lateral-data search OPEN; that was its status when written. The lateral part of P6 is recorded in `p6-lateral-data-search.md` (committed at b3eef90): **EVIDENCE GAP**. No measured, independent, tabulated He/C/O lateral widths or angular spreads in water with stated uncertainties were found. This record is not otherwise changed.
