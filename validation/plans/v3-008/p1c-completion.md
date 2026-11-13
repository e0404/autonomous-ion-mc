# V3-008 P1 completion (order step 5, freeze condition FC-1): record draft, 2026-11-01

Scope: (1) attribution of the P1b Pb residual on ICRU 49 ranges; (2) a compound/water θ_Hanson
check. The committed reference `validation/scripts/ion_em/moliere_hanson.py` (sha256 `5a99b3b8…`)
is used unmodified. Variants (Fano −u_in, Δb off, compound rules, node mapping) are runtime
monkeypatches in scratch scripts only. All numbers are ratios or differences. No NIST table value
is reproduced. The 1.0 % tolerance is not changed.

## Inputs and provenance

- Gottschalk, arXiv:0908.1413 (ScatPowerV2.tex), from the cached object
  `.ionmc-cache/ionmc-data/objects/67fb1478…` (sha256 checked by `gottschalk_tex`). The tables used:
  - `tbl:theta0Single`: x/R1, ρx, E and θ_Hanson for Be/Al/Cu/Pb at 158.6 MeV, with ρR1 in each header.
  - `tbl:B`: the B = 4.5 targets, which include water and air.
  - `tbl:modStack` rows 1 (Pb only) and 12 (Lexan only), at 230 MeV.

  The text states that ICRU 49 (MIXED.RET) is used throughout, with a cubic spline in ln R(ln T) (l. 411, 417).
- Makarova, Gottschalk, Sauerwein, arXiv:1610.01279, e-print sha256 `f8aedc66…0804175`, `main.tex`:
  - `tab:results`: the Go93 measured θ0 and the LOOKUP θ0(Hanson), Molière/Fano/Hanson with MIXED = ICRU 49 ranges (Janni 82 for Nylon).
  - `tab:materials`: the Go93 compositions.
- NIST PSTAR (ICRU 49) CSDA ranges, used as an independent cross-check only. Downloaded 2026-11-01
  to the scratchpad and not cached. The materials are Pb 082, Be 004, Al 013, Cu 029 and polystyrene 226. Only ratios are reported.

## 1. Pb attribution: CLEARED (input difference in the range table, not the reference physics)

**Range table used by the P1/P1b evaluation.**
- The test's `_path` uses the analytic ionmc `BetheStoppingSource()` defaults: Pb I = 823 eV, shell, Barkas, Bloch and Mott corrections, and no Sternheimer data for Pb (δ = 0).
- Each point is evaluated at x = (x/R1)·R1(ours).
- Gottschalk used ICRU 49.

**Range comparison.** All values are ours/reference − 1.

| | Pb | Cu | Al | Be |
|---|---|---|---|---|
| R1(ours)/ρR1(Gottschalk) | **+0.846 %** | −0.179 % | −0.003 % | +0.038 % |
| R1(ours)/R1(PSTAR) | +0.822 % | −0.241 % | −0.026 % | −0.008 % |
| R1(PSTAR)/ρR1(Gottschalk) | +0.024 % | +0.062 % | +0.023 % | +0.046 % |

**Residual range at Gottschalk's E, Pb.** Ours/Gottschalk − 1 is +0.85 % at 158 MeV, +0.93 % (149 MeV),
+1.46 % (105 MeV), **+3.12 % (41 MeV)** and **+2.56 % (20 MeV)**. PSTAR/Gottschalk is within ±0.05 % at every node.
Gottschalk's range-energy columns are therefore ICRU 49 to ≤ 0.06 %. ionmc's Bethe Pb table is the outlier.

**Sensitivity of θ to depth near the end of range.**
- d ln θ / d ln x on a fixed path is 0.79, 2.27 and 5.66 at x/R1 = 0.5, 0.9 and 0.97 (Pb).
- Be/Al/Cu are almost the same.
- A residual-range difference of a few per cent therefore becomes a θ difference of order 1 % at x/R1 ≥ 0.9.

**(a) Range table.** θ_MH/θ_Hanson − 1 for Pb, in %. "ICRU 49 path" means the committed reference
evaluated at x = (x/R1)·ρR1(G) with E(x) from:
- K: the ionmc table × a smooth κ(E) through Gottschalk's (E, ρR1 − ρx) nodes;
- S: a natural cubic spline of ln R(ln T) through Gottschalk's nodes only (his method).

Two node sets were used for each.

| x/R1 | 0.001 | 0.01 | 0.1 | 0.2 | 0.5 | 0.9 | 0.97 | max abs |
|---|---|---|---|---|---|---|---|---|
| P1b baseline (ionmc Bethe ranges) | +0.51 | +0.48 | +0.49 | +0.51 | +0.61 | **+1.10** | **+1.32** | 1.32 FAIL |
| Gottschalk's absolute ρx on our path | −0.08 | +0.07 | +0.00 | −0.01 | −0.05 | −0.76 | −2.97 | (mixed input, shows the sensitivity) |
| ICRU 49 path, K (all nodes) | +0.009 | +0.001 | +0.001 | +0.000 | −0.011 | +0.004 | +0.029 | **0.029** |
| ICRU 49 path, S (all nodes) | +0.009 | −0.002 | +0.017 | +0.004 | +0.001 | −0.044 | −0.013 | **0.044** |
| PSTAR path (spline ln R(ln T)) | +0.024 | +0.015 | +0.013 | +0.013 | +0.002 | −0.022 | −0.020 | **0.024** |

The maximum over the four ICRU 49 path variants is 0.044 %.

Controls on the same ICRU 49/PSTAR paths, max |dev|:
- Be: 0.080 % / 0.105 % (the 0.001 point);
- Al: 0.029 % / 0.030 %;
- Cu: 0.030 % / 0.043 %.

Be, Al and Cu stay passing; Cu's baseline −0.21 % becomes ≤ 0.04 %.

Further Pb points, all ICRU 49:
- Makarova LOOKUP 31.566 g/cm², x/R1 0.875, on the PSTAR path: **+0.02 %**.
  On the ionmc Bethe path at the same absolute thickness it is −0.49 %.
- Go93-thin 0.029 and 0.907 g/cm²: +0.01 % and −0.00 %.
- 0908.1413 `tbl:modStack` row 1 (230 MeV, 6.429 g/cm²): +0.01 %.
  Our exit energy is −0.11 MeV relative to the table.

**Attribution.** All of the Pb residual (the flat +0.5 % and the growth to +1.32 % at 0.97) comes from
the range-energy input:
- The flat part comes from x = (x/R1)·R1(ours) with R1(ours) 0.85 % too long. θ ∝ x^≈0.55 at small x.
- The growth comes from the end-of-range residual-range excess of +2.6 % to +3.1 %, amplified by d ln θ/d ln x ≈ 2.3 to 5.7.

**(b) Fano form for high Z.**
- For Pb, Δb = 0.048 to 0.061 along the slabs, against b ≈ 5.7 to 12.3.
- Switching Δb off changes the deviation by −0.54 % (thin) to −0.19 % (0.97).
- −u_in from 3.6 to 6.0 moves it by −0.15 % to +0.11 %.
- Bethe Z(Z+1) adds +0.18 % to +0.47 %.
- None of these produces the end-of-range growth (all are flat or shrinking toward 0.97).
- With the ICRU 49 path the committed Fano/Z² form matches to ≤ 0.044 %, so no Fano change is indicated.

**(c) Compound rule.** Not applicable: Pb is one element, and Δb reduces to D/Z (tested).

**(d) Quadrature and interpolation.** Varying the following changes the deviations by at most 8e-7 (Pb) and 1.2e-6 (all):
- Gauss–Legendre n from 50 to 3200;
- uniform GL with 4000 nodes instead of the clustered map;
- the Bethe table at 1000 instead of 200 points/decade.

**Verdict.** Pb is CLEARED by attribution. With Gottschalk's own ICRU 49 range-energy relation, the
committed reference reproduces his Pb θ_Hanson to ≤ 0.044 % at every point (≤ 0.024 % on PSTAR).
The P1b failure (+1.10 %, +1.32 %) is a documented input difference: the ionmc analytic Bethe Pb
range is +0.82 % above ICRU 49 at 158.6 MeV and up to +3.1 % in residual range at 20 to 41 MeV. It
is not caused by the Molière/Fano/Hanson physics, the Fano form or the numerics.

The 1.0 % tolerance is unchanged. The committed test, as it stands (ionmc Bethe ranges), still
exceeds 1.0 % at x/R1 ≥ 0.9, so its `xfail(strict=True)` remains a true statement. How the record
and test are re-bound (for example a P1 test on the ICRU 49 node path) is the orchestrator's decision.

**Side finding, outside P1.** ionmc's analytic Bethe Pb stopping gives ranges 0.8 to 3 % longer than
ICRU 49 at 20 to 160 MeV. For Cu, R1 is −0.24 %. This is relevant to any high-Z use of the Bethe tables.

## 2. Water/compound θ_Hanson check

**Search trail.**
- Go93 (NIM B 74 (1993) 467) is paywalled and not read.
- Its 14 materials are Be, polystyrene, C, Lexan, Nylon, Lucite, Teflon, Al, Cu, Zn, brass, Sn, Pb and U (arXiv:1610.01279 l. 61, `tab:materials`). **There is no water target.**
- arXiv:1610.01279 `tab:results` gives, for selected Go93 points, the measured θ0 and the LOOKUP θ0(Hanson) (Molière/Fano/Hanson, ICRU 49).
- The Go93 1σ errors are "not shown" (l. 170), so the measurement uncertainty is not accessible.
- The arXiv:0908.1413 measured polystyrene/Pb data exist only as figures (ExptFigPoly/Lead.eps), not tables.

**Compound check against Gottschalk's θ_Hanson** (theory vs theory, the P1 analogue; bound ≤ 1.0 %).
The table gives ours/LOOKUP − 1 with our Bethe path (ICRU 49 I-values).

| Material (g/cm²) | x/R1 | ours/LOOKUP − 1 | Bethe-Z(Z+1) variant |
|---|---|---|---|
| Polystyrene 0.347 | 0.02 | −0.16 % | +1.49 % |
| Polystyrene 15.751 | 0.89 | −0.20 % (PSTAR path −0.16 %) | +4.00 % |
| Lexan, C16H14O3, 0.094 / 1.455 | 0.005 / 0.08 | −0.15 % / −0.10 % | +0.39 / +2.36 % |
| Lexan, Go93 listed composition (H 0.074), same points | | −0.57 % / −0.54 % | |
| Lucite 0.366 / 1.449 | 0.02 / 0.08 | −0.18 % / −0.15 % | +1.48 / +2.34 % |
| Nylon (Go93 composition, Janni ranges) 0.093 / 3.010 | 0.005 / 0.18 | −0.25 % / −0.16 % | +0.34 / +2.79 % |
| Teflon 0.055 (constant E, loss < 0.3 MeV) | 0.003 | −0.05 % | +0.20 % |
| Carbon 0.316 / 1.616 (element) | | +0.01 % / +0.02 % | |
| Lexan C16H14O3, 230 MeV, 30.082 g/cm² (`tbl:modStack` row 12) | 0.87 | −0.04 % (exit-energy-pinned −0.06 %) | |

- Every compound point is within 1.0 %. The maximum is 0.25 % with stoichiometric compositions and 0.57 % with the Go93-listed Lexan composition.
- The C16H14O3 Lexan agreement at 230 MeV (−0.04 %, exit energy within 0.04 MeV) indicates that LOOKUP uses stoichiometric polycarbonate.
- The comparison discriminates between compound rules. Against the committed rule (Δb = Σ w Z/A D / Σ w Z²/A):
  - an electron-weighted Δb average would shift θ by +4.1 % (polystyrene) and +4.6 % (Lucite);
  - a mass-weighted average would shift it by +1.8 % and +2.2 %.

  The derived compound Fano rule is therefore consistent with Gottschalk's implementation to about 0.2 %. This closes OPEN 7 for H/C/N/O/F compounds.
- Informative only (no uncertainty available): ours/measured − 1 tracks LOOKUP/measured − 1 within 0.6 %. Examples:
  - polystyrene 15.751: −0.34 % vs −0.14 %;
  - Lucite: −0.57 / +0.26 % vs −0.39 / +0.41 %.

**Water.** No measured water θ0 exists in Go93, and no tabulated water θ_Hanson at a usable precision was found.
The only water θ_Hanson is 0908.1413 `tbl:B` (158.6 MeV, 2.35 µm, B = 4.5, θ = 0.051 mrad, two digits):
- ours/G − 1 = −1.3 %, against a rounding half-interval of ±1.0 %. This is not usable against a 1.0 % bound.
- As a b check, our B at that thickness is 4.448. Δb = −0.040 ± 0.002 (thickness rounding).
  The same check gives +0.004 for Be, −0.006 for Cu, +0.005 for Pb, and +0.022 ± 0.021 for air.
- A b shift of −0.040 changes water θ by −0.13 % to −0.16 % at 150 MeV, x/R 0.1 to 0.9 (B ≈ 14.5 to 17).

Outcome: **water is not validated directly (no accessible measured or precise θ_Hanson source).** It is
supported indirectly by the H-containing compounds (≤ 0.25 %) and the `tbl:B` b check (≈ −0.15 % at I11-like B).

## Not verified

- The Go93 measurement uncertainties.
- The Go93 paper itself.
- The exact compositions inside LOOKUP; the Lexan stoichiometry is inferred from the agreement.
- The meaning of the `tbl:modStack` E_2 column for row 12 (taken as the exit energy).
- The PSTAR downloads are not in the hash-pinned repository cache.

**Note on the test cache path.** The test's default `GOTTSCHALK_CACHE` is `.ionmc-cache/reference`, but the object is
in `.ionmc-cache/ionmc-data/objects`. The Pb/P1 tests therefore skip unless `IONMC_CACHE_DIR=.ionmc-cache/ionmc-data`.

## Committed provenance (2026-11-02)

Codex REVIEW-100e7fb6 finding 2: the scratch calculation is now committed and reproducible.
The Pb clearing remains an ATTRIBUTION (input difference in the range table), not a pass of the
committed Pb test, which keeps its strict xfail on ionmc Bethe ranges.

- Module `validation/scripts/ion_em/p1c_attribution.py` (sha256 `5ec7482a0604b78129e38ca59bdd0b29acccdbd1c423aa398ba33093661b500a`). Run
  `python validation/scripts/ion_em/p1c_attribution.py [--out FILE]`: it reproduces the aggregates
  of this record and writes `validation/plans/v3-008/p1c-aggregates.json` (ratios only, no PSTAR
  values). The committed `moliere_hanson.py` is used read-only (variants are explicit arguments of
  `slab_theta`, equal to `MH.moliere_slab` at the defaults; tested). The JSON equals the scratch
  `p1c-aggregates.json` of the researcher exactly (701 numbers, identical).
- Tests `tests/ionmc/test_ion_em_p1c.py` (sha256 `d95597c79d74e17e1cf5ef20e2436f719716a359527239435afc9f8f6b07f13e`): parsers on the cached inputs, Pb baseline
  (+1.10 / +1.32 %), Pb on the Gottschalk-node spline path (<= 0.05 %), on PSTAR (<= 0.03 %), the
  Be/Al/Cu controls (<= 0.11 %), R1 ratios, Fano/quadrature controls (<= 1.2e-6), the compound
  LOOKUP ratios (+-0.03 pp), the excluded compound Fano rules (+4.1 % / +1.8 %) and the water b check.
- Hash-pinned inputs (cache; looked up by sha256 in `$IONMC_CACHE_DIR[/objects]`,
  `.ionmc-cache/reference`, `.ionmc-cache/ionmc-data/objects`; hash verified at read):
  - `67fb1478e51534064f0e5363b1ee4160da2a37a600f51897b56dd74ea22ba5ed`: arXiv:0908.1413 e-print
    (ScatPowerV2.tex; `tbl:theta0Single`, `tbl:B`, `tbl:modStack`).
  - `934fd31644eb97734a0f4c9a756a40d49012d3c126ef2b4391e088be7a092aa4`: arXiv:1610.01279 PDF
    (Table 1; text via `mutool draw -F txt` or `pdftotext`; the researcher used the LaTeX e-print,
    the PDF has the same table and the parsed values agree).
  - NIST PSTAR POST responses (public domain; plain text, material name printed in the response;
    CSDA range column parsed at run time; same hashes as the researcher's scratch files), POST to
    `https://physics.nist.gov/cgi-bin/Star/apdata.pl` in the form of the registry `_STAR_BODY`
    (`prog=PSTAR&matno=NNN&ShowDefault=on&NumofEnergies=0&character=space&electronic=on&nuclear=on&total=on&csda=on&project=on&detour=on`;
    body inferred from the registry template, the acquisition record is authoritative):
    - Be (004) `0fcbb24d9b7f07af7898ac3954802d08010d22339e195ee7cd18fe2638380757`
    - Al (013) `3fe6459ec4abde6052bcb63f427ea856bbc35dd8506014191145b4eb9253db5e`
    - Cu (029) `986e8577db37120c405f6686394c1250c3e7f1d6af7e88547b964de2e99e0063`
    - Pb (082) `11dbcd1b014d926ef3e5b8bdd2bf9e539c153b5a4c8d51a62d81c34f5b0a05a7`
    - polycarbonate (219) `8974efe7242af76f1df4fb696a7b957c696f6ca7de711e01c3c68a2f40c04040` (parser test only)
    - PMMA (223) `1e13828ffc3dee1e5748d2d922e5ecf2c583380903bb787f84fbe2886dc549e9` (parser test only)
    - polystyrene (226) `dc65ded99e0cb9b455482404642054c628f0639fd75ed8ec5ebc0e5e68f981c1`
- LOGDIRs (`~/.local/share/ionmc-experiment/experiment-v3/telemetry/test-logs/`):
  `20261102T011850Z-p1c-final` (module run writing the aggregates JSON, exit 0),
  `20261102T012449Z-p1c-tests2` (p1c, moliere_hanson, ls_option, seed_families: 68 passed, 2 xfailed),
  `20261102T012723Z-p1c-stopping` (`test_stopping.py -k "not real_"`: 49 passed).
- The test cache-path default of `test_ion_em_moliere_hanson.py` now also searches
  `.ionmc-cache/ionmc-data/objects` (required-mode behaviour unchanged), so the P1 proton tests run
  without `IONMC_CACHE_DIR`. This resolves the "PSTAR downloads are not in the hash-pinned cache"
  item of "Not verified" above.

**Dated note (2026-11-02, Codex REVIEW-2c2b8078 minor):** the aggregates writer now emits standard JSON (`allow_nan=False`; NaN/inf → `null`, documented by the `_null_means` key; `bethe_b_minus_b_G` is undefined for that configuration); the artifact was regenerated by the committed module (logged run 20261102T013457Z-p1c-regen) and the hashes above are updated accordingly.
