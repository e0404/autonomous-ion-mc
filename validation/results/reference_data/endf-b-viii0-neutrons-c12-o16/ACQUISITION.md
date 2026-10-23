# ENDF/B-VIII.0 neutron sublibrary: acquisition, identity and coverage record (C-12 MAT 625, O-16 MAT 825)

Exploratory acquisition for the V3-005D M3 candidate. **No construction, calibration or evaluation use.** This page is rendered from `acquisition.json` by `validation/scripts/reference_data/endf_verify.py`; no absolute cross-section value, Legendre coefficient or per-point distribution is recorded (decision 0041).

## Acquisition

- URL `https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_neutrons.zip`; 310,158,456 bytes; sha256 `a1cced11710eb36690c5ad7b64d414465c5849f5428904049c864e9beca28894`; retrieved 2026-10-23T18:37:16Z; acquisition id `7490952096b342c98f533fab81d1b1aa`; role `exploratory`; acquisition code at `c790df36`; 561 zip entries.
- Authorisation: response of 2026-10-19 to IR-20261018-091819-903058 (bounded, acquisition-first Option B)
- Authorisation: IR-20261019-000144-2C50C8 resolved 2026-10-23T18:38:00Z (V3-INFRA-001 integrated via PR #83 at c790df369c03d1d775c843d4f1e13ab9337137e5)
- Citation: D. A. Brown et al., Nucl. Data Sheets 148 (2018) 1-142, doi:10.1016/j.nds.2018.02.001; license basis: public NNDC distribution of the US evaluated nuclear data library.

## Provenance of this record

- Verifier `validation/scripts/reference_data/endf_verify.py` sha256 `09795c35d0b126371c30cb2704da63a52a4db4f94ab01909bf7bfec448ca70b3`; git HEAD `76c17e4a1fdd2116e6ff6e4c1837aa1c58cb770c`; working tree clean: False.
- Command: `python -m validation.scripts.reference_data.endf_verify --zip .ionmc-cache/reference/a1cced11710eb36690c5ad7b64d414465c5849f5428904049c864e9beca28894 --member ENDF-B-VIII.0_neutrons/n-006_C_012.endf --member ENDF-B-VIII.0_neutrons/n-008_O_016.endf --changelog ENDF-B-VIII.0_neutrons/CHANGELOG.txt --acquisition-meta validation/results/reference_data/endf-b-viii0-neutrons-c12-o16/acquisition-meta.json --out validation/results/reference_data/endf-b-viii0-neutrons-c12-o16`
- Source archive sha256 `a1cced11710eb36690c5ad7b64d414465c5849f5428904049c864e9beca28894`.

## Identity of the members

| member | MAT | ZA | AWR | bytes | sha256 | evaluation (MF1/MT451) |
|---|---|---|---|---|---|---|
| `n-006_C_012.endf` | 625 | 6012 | 11.8936 | 3,805,168 | `daa7901e0067f9ac6d9b2392514ba00e71e37d00c8a8ecbbde0620237815f981` | 6-C - 12 LANL,ORNL  EVAL-AUG15 G.M. Hale, P.G. Young, C.Y. Fu / NDS 148, 1 (2018)    DIST-FEB18 |
| `n-008_O_016.endf` | 825 | 8016 | 15.8575 | 41,873,568 | `aabef9c710c132b4e4d8028e58eff2e542d61d34b93823b5f1199649bff7f781` | 8-O - 16 LANL       EVAL-DEC16 Hale,Paris,Young,Chadwick / CH96AB,CA98,Ch99,Yo01 DIST-FEB18 |

| member | NLIB | NVER | NFOR | NSUB | NMOD | LREL | EMAX (MeV) | NWD | NXC |
|---|---|---|---|---|---|---|---|---|---|
| `n-006_C_012.endf` | 0 | 8 | 6 | 10 | 0 | 0 | 20 | 79 | 57 |
| `n-008_O_016.endf` | 0 | 8 | 6 | 10 | 3 | 0 | 30 | 963 | 192 |

## Coverage (MF3/MT2 and MF4/MT2)

| member | MF3/MT2 range (MeV), points | LTT / LI / LCT | Legendre part | tabulated part |
|---|---|---|---|---|
| `n-006_C_012.endf` | 1e-05 eV to 150 MeV, 1077 | 3 / 0 / 2 | 1e-05 eV to 20 MeV, NE 702, NL <= 8 | 20 MeV to 150 MeV, NE 17, NP <= 91, INT [4] |
| `n-008_O_016.endf` | 1e-05 eV to 150 MeV, 2323 | 3 / 0 / 2 | 1e-05 eV to 30 MeV, NE 1622, NL <= 10 | 30 MeV to 150 MeV, NE 14, NP <= 73, INT [4] |

- `n-006_C_012.endf` Legendre-part incident energies >= 20 MeV (MeV): 20
- `n-006_C_012.endf` tabulated-part incident energies (MeV): 20, 23, 27, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150
- `n-008_O_016.endf` Legendre-part incident energies >= 20 MeV (MeV): 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30
- `n-008_O_016.endf` tabulated-part incident energies (MeV): 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150

## Per-energy statistics of the tabulated f(mu) (CM; aggregates only)

### `n-006_C_012.endf`

| E (MeV) | NP | INT | norm | <1 - mu> | P(theta > 10 deg) | P(theta > 20 deg) | local minima forward of 130 deg (deg) |
|---|---|---|---|---|---|---|---|
| 20 | 91 | [4] | 1.000000 | 0.21915 | 0.8955 | 0.6497 | 64 |
| 23 | 91 | [4] | 1.000000 | 0.18066 | 0.8818 | 0.6110 | 64, 126 |
| 27 | 91 | [4] | 1.000000 | 0.14890 | 0.8648 | 0.5662 | 118 |
| 30 | 91 | [4] | 1.000000 | 0.13105 | 0.8519 | 0.5340 | 114 |
| 35 | 91 | [4] | 1.000000 | 0.10885 | 0.8302 | 0.4833 | 108 |
| 40 | 91 | [4] | 1.000000 | 0.09272 | 0.8080 | 0.4360 | 104 |
| 50 | 91 | [4] | 1.000000 | 0.07152 | 0.7643 | 0.3547 | none |
| 60 | 91 | [4] | 1.000000 | 0.05870 | 0.7227 | 0.2917 | 128 |
| 70 | 91 | [4] | 1.000000 | 0.05056 | 0.6852 | 0.2462 | 84, 122 |
| 80 | 91 | [4] | 1.000000 | 0.04547 | 0.6544 | 0.2168 | 46, 78, 116 |
| 90 | 91 | [4] | 1.000000 | 0.04271 | 0.6335 | 0.2021 | 42, 76, 110 |
| 100 | 91 | [4] | 1.000000 | 0.04206 | 0.6272 | 0.2015 | 40, 72, 106 |
| 110 | 73 | [4] | 1.000000 | 0.02952 | 0.5773 | 0.1166 | none |
| 120 | 73 | [4] | 1.000000 | 0.02680 | 0.5469 | 0.0954 | none |
| 130 | 73 | [4] | 1.000000 | 0.02450 | 0.5172 | 0.0779 | none |
| 140 | 73 | [4] | 1.000000 | 0.02251 | 0.4879 | 0.0636 | none |
| 150 | 73 | [4] | 1.000000 | 0.02075 | 0.4590 | 0.0518 | none |

### `n-008_O_016.endf`

| E (MeV) | NP | INT | norm | <1 - mu> | P(theta > 10 deg) | P(theta > 20 deg) | local minima forward of 130 deg (deg) |
|---|---|---|---|---|---|---|---|
| 30 | 73 | [4] | 1.000000 | 0.11993 | 0.8291 | 0.4820 | 52.5 |
| 35 | 73 | [4] | 1.000000 | 0.10030 | 0.8097 | 0.4388 | none |
| 40 | 73 | [4] | 1.000000 | 0.08676 | 0.7915 | 0.4009 | none |
| 50 | 73 | [4] | 1.000000 | 0.06762 | 0.7497 | 0.3282 | none |
| 60 | 73 | [4] | 1.000000 | 0.05279 | 0.7081 | 0.2609 | none |
| 70 | 73 | [4] | 1.000000 | 0.04341 | 0.6685 | 0.2072 | none |
| 80 | 73 | [4] | 1.000000 | 0.03694 | 0.6308 | 0.1644 | none |
| 90 | 73 | [4] | 1.000000 | 0.03220 | 0.5946 | 0.1305 | none |
| 100 | 73 | [4] | 1.000000 | 0.02855 | 0.5600 | 0.1036 | none |
| 110 | 73 | [4] | 1.000000 | 0.02564 | 0.5266 | 0.0823 | none |
| 120 | 73 | [4] | 1.000000 | 0.02325 | 0.4944 | 0.0656 | none |
| 130 | 73 | [4] | 1.000000 | 0.02123 | 0.4632 | 0.0525 | none |
| 140 | 73 | [4] | 1.000000 | 0.01949 | 0.4328 | 0.0422 | 122.5 |
| 150 | 73 | [4] | 1.000000 | 0.01796 | 0.4032 | 0.0341 | 115, 122.5 |

## sigma_el ratio (MF3/MT2 over the committed BGG/Barashenkov elastic formula)

| member | Z | 100MeV | 150MeV |
|---|---|---|---|
| `n-006_C_012.endf` | 6 | 1.0356 | 0.9040 |
| `n-008_O_016.endf` | 8 | 0.9659 | 0.9308 |

## Covariance and other sections

- `n-006_C_012.endf`: MF [1, 2, 3, 4, 5, 6, 12, 14, 33]; covariance MF [33]; MF33 MT [1, 2, 4, 5, 28, 51, 102, 103, 104, 107]; MF34 MT []; MF6/MT2 present: False; no MF34/MT2.
- `n-008_O_016.endf`: MF [1, 2, 3, 4, 6, 12, 13, 14, 33, 34]; covariance MF [33, 34]; MF33 MT [1, 2, 4, 5, 16, 22, 23, 28, 32, 41, 44, 45, 102, 103, 104, 105, 107, 108, 112, 600, 601, 602, 603, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 749, 800, 801, 802, 803]; MF34 MT [2]; MF6/MT2 present: False; MF34/MT2 NMT1=1; MT1=2 NL=1 NL1=1: L=1/L1=1 LB [5] E 1e-05 eV-7.005 MeV.

## MF1 comment lines mentioning 'retain' (verbatim)

- `n-006_C_012.endf` comment line 33: ` to 150 MeV (Chadwick and Young) was not retained in this file.`
- `n-008_O_016.endf` comment line 126: `it used to be before. To retain continuity, I introduced a region`

## CHANGELOG quotes

`ENDF-B-VIII.0_neutrons/CHANGELOG.txt`: 99,474 bytes, sha256 `45933f47ceb3c02c402df25382515beebfdf65c2f0a6510a6bd37927ba1ae47f`, 1969 lines. Verbatim lines mentioning the members:

- line 14 (n-008_O_016.endf): `    * n-008_O_016.endf: Make 1/v portion of capture cross section log-log`
- line 32 (n-006_C_012.endf): `    * n-006_C_013.endf, n-006_C_012.endf: Ensure total cross section is`
- line 143 (n-008_O_016.endf): `    * n-008_O_016.endf: Add covariance for many channels.  We suspect format problems`
- line 202 (n-006_C_012.endf): `    * n-006_C_012.endf: Revised from Hale & Paris: natC evaluation restored from 20-150 MeV.`
- line 262 (n-008_O_016.endf): `    * n-008_O_016.endf: Covariances added`
- line 323 (n-006_C_012.endf): `    * n-006_C_012.endf: Added covariances`
- line 434 (n-008_O_016.endf): `    * n-008_O_016.endf:`
- line 886 (n-008_O_016.endf): `    * n-008_O_016.endf: new evaluation of capture cross section based on the R-matrix (up to`
- line 1272 (n-006_C_012.endf): `    * n-006_C_012.endf: 2 fixes:`
- line 1352 (n-008_O_016.endf): `* Bug fixes in several files (n-008_O_016.endf, n-006_C_013.endf, n-092_U_235.endf,`
- line 1368 (n-008_O_016.endf): `    * n-008_O_016.endf: CIELO 16O file from G. Hale, et al. (LANL)`

## Findings

- Both files: MF3/MT2 (elastic cross section) from 1e-5 eV to 150 MeV; MF4/MT2 with LTT = 3 (Legendre below, tabulated above a switch energy), LCT = 2 (centre-of-mass). Tabulated f(mu) uses INT 4 (log-y, lin-mu) throughout; the statistics below honour it, with the exceedance probabilities interpolated exactly to the threshold angle.
- C-12: Legendre part to 20 MeV (NL <= 8; its only incident energy >= 20 MeV is 20 MeV); tabulated f(mu) at 17 incident energies from 20 to 150 MeV (NP <= 91). O-16: Legendre part to 30 MeV (NL <= 10; the incident energies 20-30 MeV listed under it are Legendre entries); tabulated at 14 incident energies from 30 to 150 MeV (NP <= 73).
- The project's energies of interest are 100, 150 and 200 MeV and the elastic table grid runs to 250 MeV: 100 and 150 MeV are inside the evaluated domain (150 MeV is the LAST node), 200 MeV and the grid above 150 MeV are OUTSIDE it. No extrapolation is authorised without scientific justification (operator, 2026-10-19).
- O1 RESOLVED BY RECORD: CHANGELOG line 202 reads 'n-006_C_012.endf: Revised from Hale & Paris: natC evaluation restored from 20-150 MeV.' The C-12 20-150 MeV content is therefore the natural-carbon evaluation restored in the VIII.0 cycle; the MF1 comment 'was not retained in this file' is stale text from the August 2015 revision (the tape label is dated 2015-08-29).
- Header anomaly: MF1/MT451 EMAX is 2.0e7 eV (C-12) and 3.0e7 eV (O-16) while the data extend to 1.5e8 eV. Recorded as is; not corrected.
- No MF6/MT2 in either file. No MF4/MT2 shape covariance above 20 MeV in C-12 (MF33 only; no MF34 section). O-16: MF34/MT2 parsed (one block, MT1 = 2, NL = 1, NL1 = 1, i.e. Legendre order L = 1 only) and its energy grid ends at 7.005 MeV, below the 30 MeV switch to the tabulated part, so no MF4/MT2 shape covariance is carried above 7.005 MeV in O-16 (stricter than the 30 MeV previously assumed).
- These are NEUTRON evaluations. Their use as a mirror for p+A elastic angular distributions (M3) is a scientific hypothesis to be assessed (isospin mirror validity, Coulomb-nuclear interference, isovector terms) before any plan; this record establishes identity, provenance and coverage only.

- Evaluation comment, `n-006_C_012.endf`: MF1/MT451: 'ENDF/B-VII MOD 1 Evaluation, August 2015, G. M. Hale and P. G. Young (LANL)'; R-matrix analysis up to 6.45 MeV; 'At neutron energies above 6.45 MeV, the evaluation was joined to the existing ENDF/B-VII.0 evaluation for natural carbon, which is essentially version VI.1 (Fu, Axton, and Perey)'; 'All cross sections were kept the same as the VII.0 file in the energy range 6.47 <= E <= 20 MeV'; 'The version VI.3 extension to 150 MeV (Chadwick and Young) was not retained in this file.' This last sentence is stale text from the August 2015 revision (see finding O1 and the CHANGELOG quote).
- Evaluation comment, `n-008_O_016.endf`: MF1/MT451: 'EVAL-DEC16 Hale,Paris,Young,Chadwick'; R-matrix below 6.25 MeV, GNASH/optical-model evaluation 6.25-30 MeV; the previous MOD 2 (LA150-type, Chadwick et al.) evaluation retained from 30 to 150 MeV ('the previous evaluation, which extends to 150 MeV'); above 30 MeV MF6/MT5 carries the non-elastic emission and 'MF=4 MT=2 Elastic Angular Distributions' are listed among the non-zero data above 30 MeV; covariances 'continuing to 150 MeV'

**Status:** identity, provenance and energy coverage VERIFIED by the committed verifier; scientific suitability recorded in decision 0042
