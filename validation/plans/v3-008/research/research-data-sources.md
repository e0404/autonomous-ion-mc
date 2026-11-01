> **Research note (physics-researcher, 2026-11-01; not evidence).** Statements are labelled FACT / LIT / TO-VERIFY / OPEN by the author; estimates are hand estimates, nothing was run. Items marked TO-VERIFY or OPEN are resolved by the pre-freeze checks of `validation/plans/v3-008-acceptance.md` before the plan is frozen.

# V3-008 data-source research (RQ-D1, RQ-D2, RQ-D3, RQ-V1) — research notes, 2026-11-01

Scope: bounded investigation, stopped early on the coordinator's instruction. Markers: **FACT** = checked
this session against a primary source (Geant4 v11.4.2 source on GitHub, a live index page, or a repository
file); **LIT** = literature statement from memory or a secondary source, not re-read; **TO-VERIFY** =
plausible but not checked; **OPEN** = not answered.
Nothing was downloaded into the repository. Scratch fetches (Geant4 .cc files, the IAEA arXiv PDF) are in
`$TMPDIR` only.

## 0. Headline finding: TOPAS ion stopping lineage differs from what the draft says (FACT)

The decision 0043 draft (sections 1 and 6, and the V1 row) says TOPAS uses the ICRU 73 C/O tables. In
Geant4 11.4.2 `g4em-standard_opt4` that is wrong above about 2 MeV/u:

- `G4EmStandardPhysics_option4.cc` (v11.4.2) calls `param->SetUseICRU90Data(true)`. GenericIon ionisation
  uses `G4LindhardSorensenIonModel` unless `fUseExternalDEDX` is set; in that case it uses
  `G4IonParametrisedLossModel`. The default is LS. (FACT)
- `G4LindhardSorensenIonModel.cc`: `fElimit = 2.0 MeV` on the proton-scaled energy (`escaled = E·m_p/M`,
  so about 2 MeV/u).
  - Below fElimit it uses `G4IonICRU73Data::GetDEDX`, i.e. ICRU 90 / ICRU 73 ion tables.
  - Above fElimit, dE/dx comes from `G4BetheBlochModel`. Along-step it adds the Lindhard–Sørensen ΔL
    (`G4LindhardSorensenData`) and a Barkas term `2·BarkasCorrection·(q−1)/q`. (FACT)
- `G4BetheBlochModel.cc` (v11.4.2), when `UseICRU90Data` is set and the particle is GenericIon in an ICRU 90
  material (G4_WATER, G4_AIR, G4_GRAPHITE): `dedx = S_p,ICRU90(E·m_p/M) × chargeSquareRatio`. That is the
  ICRU 90 **proton** array at equal velocity, scaled by the Geant4 effective charge squared. For alpha it
  uses the ICRU 90 **alpha** array up to 1 GeV total (`fAlphaTlimit`). (FACT)
- `G4IonICRU73Data` resamples its files only over 0.025–2.5 MeV (scaled) (`fEmin`, `fEmax`). (FACT)
  - Directory: `$G4LEDATA/ion_stopping_data/icru90/z<Z>_<G4 material>.dat` for Z ≤ 18 in
    G4_AIR/G4_GRAPHITE/G4_WATER (density coefficient 0.998 for water).
  - Otherwise `icru73/z<Z>_<material>.dat`; element files `z<Z>_<Ztarget>.dat`.
- **Consequence:** the TOPAS EM-only C/O range lineage in water is ICRU 90 proton stopping + Geant4
  effective charge + Geant4 LS ΔL + Geant4 Barkas.
  - ionmc's Barkas F table is transcribed from G4EmCorrections, so that part is shared code lineage.
  - I = 78 eV is shared.
  - The ICRU 73/90 C/O tables enter TOPAS only below about 2 MeV/u (the last tens of µm of range).
  - So an ICRU 73 (rev.) or ICRU 90 C/O table comparison (row I3) is *more* independent of TOPAS than the
    draft assumes. The draft's rejection reason "copy the stopping data TOPAS uses" for a tabulated C/O
    default does not hold above 2 MeV/u.
- **Attribution aid for STOP 3:** ionmc can rebuild "TOPAS-like" C/O water stopping from data it already has
  (registered `geant4-icru90-stopping-11.4.2` proton array) plus LS ΔL. No new data needed.
- Scratch copies with SHA-256: `G4IonICRU73Data.cc` = cf771a5a…08d2 (13 701 B), fetched from
  `https://raw.githubusercontent.com/Geant4/geant4/v11.4.2/source/processes/electromagnetic/standard/src/`.
- Sibling files, same tag:
  - `.../standard/src/G4LindhardSorensenIonModel.cc`
  - `.../standard/src/G4BetheBlochModel.cc`
  - `source/physics_lists/constructors/electromagnetic/src/G4EmStandardPhysics_option4.cc`
- The repository's V1 rationale (TOPAS EM-only) should declare this lineage. The V1 case file should enable
  `/process/em/printParameters` or verbose output to confirm at run time that the LS model and ICRU90 flag
  are active (TO-VERIFY at run time).

## RQ-D1 Public C/O stopping tables in liquid water, 1–500 MeV/u

### (a) ICRU 73 + 2009 water erratum; (b) ICRU 90 carbon

- **Geant4 G4EMLOW** (data set for Geant4 11.4.2 = **G4EMLOW 8.8**; FACT from
  `cmake/Modules/G4DatasetDefinitions.cmake` at v11.4.2, MD5 `328330009df633f7e9b3a9f445745298`).
  - URL: `https://cern.ch/geant4-data/datasets/G4EMLOW.8.8.tar.gz`, which redirects to
    `https://geant4-data.web.cern.ch/datasets/G4EMLOW.8.8.tar.gz`. HEAD gives Content-Length
    **350 393 595 B (≈ 334 MiB)**. (FACT)
  - It contains `ion_stopping_data/icru90/` (Z = 3–18 in water, air, graphite) and `ion_stopping_data/icru73/`
    (31 ICRU 73 materials incl. G4_WATER, G4_WATER_VAPOR, plus element files). Existence is inferred from
    the loader code. (FACT for paths in code; file presence and contents TO-VERIFY)
  - Format: `G4PhysicsFreeVector::Retrieve(ascii)` text. Values appear to be in MeV cm²/mg (the loader
    multiplies by density·MeV·1000·cm²/g); energy is per nucleon (TO-VERIFY).
  - **Energy span of the files: OPEN.** The loader uses only ≤ 2.5 MeV/u. Older Geant4 docs say the
    ICRU 73 data go "up to 1 GeV/nucleon" (G4IonParametrisedLossModel, since 9.3; LIT, secondary search
    result). If the files hold the full ICRU 73 span (0.025–1000 MeV/u), they are a usable I3 reference.
    If truncated, they are not.
  - Whether `icru73/z6_G4_WATER.dat` carries the **2009 erratum** values (I = 78 eV) or the original ICRU 73
    water (which used I = 67.2 eV): TO-VERIFY by comparison with literature values. The size of the
    difference (a few % at therapy energies) makes this easy to tell.
  - The `icru90/z6_G4_WATER.dat` file: it is OPEN whether it holds ICRU 90 published carbon values or a
    Geant4-team recomputation for Z = 3–18 (ICRU 90 tabulates only p, α, C).
  - Licence: Geant4 Software Licence (redistribution allowed with notice). Fine for the cache; values still
    stay out of Git by decision 0041.
  - Cost: 334 MiB download, within the 1 GiB cap (V3-INFRA-001), for perhaps a few kB of needed files. The
    provenance tool would pin the tarball hash and extract the members by path. No smaller official
    per-file mirror was found (search for `z6_G4_WATER.dat` returned nothing).
  - Alternative cheap probe: some third-party Geant4 data mirrors exist (KEK LXR shows source only).
    Not checked.
- **ICRU 73 erratum document**: Sigmund, Schinner, Paul, "Errata and Addenda for ICRU Report 73", J. ICRU
  5(1) (2009). LIT. Public PDF availability TO-VERIFY. H. Paul's former site (exphys.jku.at/stopping) and
  the IAEA "legacy" pages may host it. Not checked.
- **ICRU 90 carbon tables**: carbon ions in water, air and graphite are tabulated in ICRU 90 (LIT; an
  arXiv paper, 1812.07877, states ICRU 90 covers p, α and C ions for water, graphite and air). The report is
  paywalled (OUP J ICRU). It is not in `G4ICRU90StoppingData.cc` (FACT, decision 0038 and the parser). The
  only public machine-readable candidate is G4EMLOW `icru90/z6_G4_WATER.dat`, with provenance OPEN as
  above.
  - For O there is no ICRU 90 table (LIT). G4EMLOW `icru90/z8_*` exists per the loader's Z ≤ 18 rule, but
    its origin is OPEN.
- **libdEdx** (Lühr, Toftegaard, Kantemiris, Hansen, Bassler, Int. J. Radiat. Biol. 88 (2012) 209; LIT).
  It bundles PSTAR/ASTAR/MSTAR/ICRU 49/ICRU 73 (with the water erratum, per memory) and Bethe. Licence GPL
  (LIT). Repository URL and the current host (SourceForge / GitHub "APTG/libdedx") TO-VERIFY; not checked
  in this pass. If the data files are in the source tree, this is a far smaller acquisition than
  G4EMLOW (likely < 5 MB).
- **NIST**: PSTAR/ASTAR cover only protons and helium. NIST has no ion (Z > 2) stopping program. LIT,
  consistent with the registry; not re-fetched.

### (c) IAEA stopping database / MSTAR

- IAEA "Electronic Stopping Power of Matter for Ions", `https://www-nds.iaea.org/stopping/` (FACT, live).
  - Current version **2026-09** (25 Sep 2026); bulk file `StoppingPower_v2026-09.zip`. A HEAD on the guessed
    download path returned 403 (wrong path or bot block); exact URL and size TO-VERIFY.
  - Has an API page (`/stopping/api`, content JS-rendered, not parsed), per-tuple CSV/TXT downloads
    (`/stopping/download/tuples`, `/stopping/download/refs`), and a `/stopping-legacy` branch.
  - Terms of use: `https://www.iaea.org/about/terms-of-use` (FACT link; content not read).
- Content for water (Montanari et al., arXiv:2402.03080, CC BY 4.0; FACT from the paper text):
  - "Only four ions in water feature energy loss experimental data: H, He, Li, and C."
  - **No O in water.**
  - He/Li data are about 40 years old. "Carbon in water is the system less examined."
  - The recent liquid-water C data (Baek et al., Phys. Rev. E 2020) are around the stopping maximum.
  - **Consequence:** measured C stopping in water exists only far below the 2 MeV/u transport floor
    (energies TO-VERIFY, expected ≲ 1 MeV/u). It is useless for I3 / range evidence.
- MSTAR (H. Paul & A. Schinner): fits to measured data, Z = 3–18. Its span is 0.001–1000 MeV/u per memory
  (LIT, TO-VERIFY). It is distributed as a program (Windows/DOS executable; IAEA "Computer Programs" page),
  not as a table. Its high-energy part is anchored on ASTAR scaling (LIT), so it is not independent of
  ICRU 49 theory. MSTAR is software → same rule as SRIM.
- (d) SRIM/TRIM: free, closed-source Windows software (J.F. Ziegler). Its outputs are a computation, not a
  versioned dataset. The project rule (decision 0043 draft) already rejects it. Keep rejected. LIT.

## RQ-D2 ASTAR coverage (FACT from repository; NIST page not re-fetched)

- The registered `nist-astar-water-2005` POST (`apdata.pl`, matno 276, `NumofEnergies=0` = default grid)
  returns electronic, nuclear and total stopping, CSDA range, projected range and detour factor.
- The `icru90.py` docstring states the alpha grid "ends at 1000 MeV, the same span as NIST ASTAR (total
  kinetic energy …, up to 250 MeV/u)".
- ASTAR default grid is 1 keV–1000 MeV total kinetic energy (LIT; consistent with the docstring).
- ASTAR water uses I = 75 eV (ICRU 49).
- **The draft's He upper bound of 250 MeV/u is therefore consistent.** I did not re-read the cached file's
  last row (to avoid unnecessary exposure). TO-VERIFY by a header-free line count in the loader test.

## RQ-D3 Independent (non-ICRU-lineage) measured stopping or range of C/O in water

- Measured *stopping* at therapy energies: none public in the IAEA database (see above). FACT for the
  database scope as described by Montanari et al.; O absent.
- Measured *ranges* (Bragg-peak positions) of C/O in water: GSI (Schardt et al., "Precision Bragg-curve
  measurements for light-ion beams in water", GSI Sci. Rep. 2007; LIT).
  - These and the proton/ion range data of Kumazaki et al. 2007 and Bichsel & Hiraoka 1992 (LIT) informed
    the ICRU 90 choice I_water = 78 eV.
  - **Any such range data is calibration lineage for I = 78 eV** (shared with Geant4/TOPAS G4_WATER and
    with ionmc's default). It is measured evidence of range *consistency*, not an independent I-value test.
    The draft already states this correctly.
- Water-equivalent plastics, CT-calibration ranges (Jäkel et al., Med. Phys. 28 (2001) 701; LIT):
  relative ranges (WEPL); not usable for absolute water range.
- **No public, numeric, non-ICRU-lineage C/O stopping or range dataset was found in this pass.** OPEN
  (bounded search incomplete; see "not checked").

## RQ-V1 Measured ion Bragg curves / ranges in water

- The repository research note `docs/research/validation-evidence.md` (FACT, read this session) already
  records:
  - "No public numeric measured Bragg curves or lateral profiles for He, C or O in water".
  - Tessonnier et al. PMB 62 (2017) 3958 is not open access, has no data statement, and its numbers are in
    figures only.
  - The GSI Steidl/Schardt report repository is behind a bot challenge.
  - Haettner 200 MeV/u is not public.
  - Intervention candidates are listed there.
- **Haettner geant-val** (`https://geant-val.cern.ch/api/getExpPlotsByInspireId?inspire_id=-5`): per that
  note it holds 384 records, all at 400 MeV/u.
  - The records are fragment yields N/N0 vs depth (Z = 1–5), yield-by-Z sets, angular distributions and
    energy spectra.
  - **No depth-dose curve is listed.** (FACT as recorded in the repo note; the file is NOT in the project
    cache, checked: neither `ionmc-data` nor the V3-005A/C backups contain it. Only `inspire_id=-7` is
    cached.)
  - Primary-carbon attenuation is a V3-009 observable, not V3-008.
- HIMAC: Sihver, Schardt, Kanai, Jpn. J. Med. Phys. 18 (1998) 1, "Depth-dose distributions of high-energy
  carbon, oxygen and neon beams in water" (LIT). J-STAGE is likely open access. Bragg-peak positions are
  possibly tabulated. **TO-VERIFY: this is the most promising public candidate for C and O measured ranges
  in water.** It is ICRU/I-value lineage-neutral as data, but pre-ICRU 90 and part of the I-value
  literature.
- HIT: Krämer et al. 2016 (He; Med. Phys.), Kurz et al. 2012 PMB 57 5017 (O) (LIT). Open-access status
  and numeric tables TO-VERIFY. They are probably figures only.
- ICRU 93 (2016) / PTCOG: no public numeric Bragg-curve tables known (LIT).
- Geant4 G4-Med (Arce et al. 2021/2025) uses Haettner Bragg data for its carbon Bragg-peak test (LIT, repo
  note). Whether geant-val exposes the Bragg curve as a separate inspireId is OPEN.
- **Independence:** every measured water range dataset is calibration lineage for I = 78 eV (see RQ-D3).
  None is independent of ICRU 90 for the I-value. They are independent for straggling, distal shape and
  the z-dependent stopping terms.

## Existing project TOPAS / reference cases (FACT, repository)

- `validation/reference_cases/topas/carbon-water-290mevu-smoke`: **full physics**: `g4em-standard_opt4`,
  `g4h-phy_QGSP_BIC_HP`, `g4ion-QMD`, `g4h-elastic_HP`, `g4stopping`, `g4decay`.
  - Settings: cut 0.05 mm, 200 histories, seed 20261003, `GenericIon(6,12)` with BeamEnergy 3480 MeV total.
  - Execution smoke only, not EM-only, not usable for V1.
- The only EM-only TOPAS cases are protons (`proton-water-150mev-lateral-emonly*`, module list = 1
  `"g4em-standard_opt4"`, cut 0.05 mm). That module string is the template for the ion EM-only bundles.
- No He or O TOPAS case exists. `validation/results/reference/` holds proton results only. FRED has a
  `carbon-water-290mevu-smoke` case (the draft says FRED aborts on 12C fragmentation; EM-only ion FRED
  untested, OPEN).
- Lineage of a TOPAS EM-only ion run: see section 0. The alpha run uses the ICRU 90 alpha array up to
  1 GeV total, i.e. the same data the project registered (shared; helium V1 range is not independent of
  the I2 reference).

## Recommended minimal acquisition list

| # | Item | URL | Size | Role | Note |
|---|---|---|---|---|---|
| 1 | Geant4 11.4.2 ion stopping model sources (G4LindhardSorensenIonModel.cc, G4BetheBlochModel.cc, G4IonICRU73Data.cc, G4EmStandardPhysics_option4.cc, G4LindhardSorensenData.cc) | raw.githubusercontent.com/Geant4/geant4/v11.4.2/source/... (paths in section 0) | ≈ 15–50 kB each | exploratory / attribution (lineage declaration; V1 failure attribution) | Geant4 licence; no table values in Git |
| 2 | G4EMLOW 8.8 tarball → extract `ion_stopping_data/icru73/z6_G4_WATER.dat`, `z8_G4_WATER.dat`, `icru90/z6_G4_WATER.dat`, `z8_G4_WATER.dat` (+ `z2`? not present: Z ≥ 3) | https://geant4-data.web.cern.ch/datasets/G4EMLOW.8.8.tar.gz (MD5 328330009df633f7e9b3a9f445745298) | 350 393 595 B | evaluation (I3) **only if** the files span ≥ 2–500 MeV/u and the water-erratum provenance is established | first verify the span; if truncated at 2.5 MeV/u, I3 has no table |
| 3 | libdEdx source (if it ships ICRU 73 rev. water tables as files) | TO-VERIFY (GitHub/SourceForge) | expected < 5 MB | evaluation (I3) alternative / cross-check of #2 | GPL; data provenance must trace to ICRU 73 errata |
| 4 | Sihver, Schardt, Kanai 1998 Jpn J Med Phys 18:1 (J-STAGE PDF) | TO-VERIFY | ≈ 1–3 MB | evaluation (V3 measured range, C and O) | pre-ICRU 90; I-value calibration lineage |
| 5 | (optional) IAEA stopping DB tuple CSV for C in H2O | via `https://www-nds.iaea.org/stopping/` tuples/API | kB | exploratory (documents that no ≥ 2 MeV/u data exist) | IAEA terms of use |

## Rows without independent evidence (intervention candidates)

- **V3 (measured range) for He, C, O**: no public numeric dataset confirmed. Sihver 1998 (C, O) is the only
  open candidate, unverified.
  - If it fails: intervention for the GSI Schardt 2007 report and HIT Tessonnier 2017 numeric tables
    (inaccessible information, EXPERIMENT.md case 1).
  - He has no measured range candidate at all.
- **I3 (C/O table)**: depends on the G4EMLOW file span and erratum provenance (OPEN). **O has no ICRU 90
  table**; any O reference is ICRU 73 (rev.) or a Geant4 file of unclear origin.
- **Measured C/O stopping ≥ 2 MeV/u in water**: none exists publicly (IAEA scope). This is not curable by
  intervention short of new experiments; declare it as a limitation.
- **He range evidence**: I2 (ASTAR/ICRU 90) and TOPAS He (ICRU 90 alpha array) share data lineage.
  Independent He evidence would need measured data (none found).

## Not checked in this pass (stopped on instruction)

- Actual contents and span of the G4EMLOW `ion_stopping_data` files.
- libdEdx location and licence.
- Public PDF of the ICRU 73 errata.
- J-STAGE availability of Sihver 1998.
- Open access for Krämer 2016 and Kurz 2012.
- IAEA bulk zip URL and size.
- Whether geant-val has a separate carbon Bragg-curve record.
- MSTAR energy span.
- The ASTAR file's last energy (no re-read).
