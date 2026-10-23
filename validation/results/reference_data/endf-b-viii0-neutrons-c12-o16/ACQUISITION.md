# ENDF/B-VIII.0 neutron sublibrary — acquisition, identity and coverage record (C-12 MAT 625, O-16 MAT 825)

Exploratory acquisition for the V3-005D M3 candidate (operator-authorised, acquisition-first). **No construction, calibration or evaluation use.** The JSON next to this file is the record; this page summarises it. No cross-section value, angular-distribution point or Legendre coefficient is recorded (decision 0041).

## Acquisition

- URL `https://www.nndc.bnl.gov/endf-b8.0/zips/ENDF-B-VIII.0_neutrons.zip`; 310,158,456 bytes; sha256 `a1cced11710eb36690c5ad7b64d414465c5849f5428904049c864e9beca28894`; retrieved 2026-10-23T18:37:16Z; acquisition id `7490952096b342c98f533fab81d1b1aa`; role `exploratory`; acquisition code at `c790df36` (V3-INFRA-001, PR #83, cap 1 GiB); 561 zip entries.
- Authorisation: operator response of 2026-10-19 to IR-20261018-091819-903058; IR-20261019-000144-2C50C8 resolved 2026-10-23.

## Identity of the two members

| member | MAT | ZA | AWR | bytes | sha256 | evaluation (MF1/MT451) |
|---|---|---|---|---|---|---|
| `n-006_C_012.endf` | 625 | 6012 | 11.89365 | 3,805,168 | `daa7901e0067f9ac6d9b2392514ba00e71e37d00c8a8ecbbde0620237815f981` | 6-C - 12 LANL,ORNL  EVAL-AUG15 G.M. Hale, P.G. Young, C.Y. Fu / NDS 148, 1 (2018)    DIST-FEB18 |
| `n-008_O_016.endf` | 825 | 8016 | 15.85751 | 41,873,568 | `aabef9c710c132b4e4d8028e58eff2e542d61d34b93823b5f1199649bff7f781` | 8-O - 16 LANL       EVAL-DEC16 Hale,Paris,Young,Chadwick / CH96AB,CA98,Ch99,Yo01 DIST-FEB18 |

NLIB 0 (ENDF/B), NVER 8, NFOR 6, NSUB 10 (incident neutron) in both.

## Coverage (MF3/MT2 and MF4/MT2)

| member | MF3/MT2 range | MF4/MT2 LTT / LCT | Legendre part | tabulated part | tabulated incident energies (MeV) |
|---|---|---|---|---|---|
| `n-006_C_012.endf` | 1e-05–150 MeV (1077 points) | 3 / 2 | to 20 MeV, NE 702, NL ≤ 8 | 20–150 MeV, NE 17, NP ≤ 91 | 20, 23, 27, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150 |
| `n-008_O_016.endf` | 1e-05–150 MeV (2323 points) | 3 / 2 | to 30 MeV, NE 1622, NL ≤ 10 | 30–150 MeV, NE 14, NP ≤ 73 | 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150 |

## Findings

- Both files: MF3/MT2 (elastic cross section) from 1e-5 eV to 150 MeV; MF4/MT2 with LTT = 3 (Legendre below, tabulated above a switch energy), LCT = 2 (centre-of-mass).
- C-12: Legendre to 20 MeV (NL <= 8); tabulated f(mu) at 17 incident energies from 20 to 150 MeV (NP <= 91). O-16: Legendre to 30 MeV (NL <= 10); tabulated at 14 incident energies from 30 to 150 MeV (NP <= 73).
- The project's energies of interest are 100, 150 and 200 MeV and the elastic table grid runs to 250 MeV: 100 and 150 MeV are inside the evaluated domain (150 MeV is the LAST node), 200 MeV and the grid above 150 MeV are OUTSIDE it. No extrapolation is authorised without scientific justification (operator, 2026-10-19).
- Header anomaly: MF1/MT451 EMAX is 2.0e7 eV (C-12) and 3.0e7 eV (O-16) while the data extend to 1.5e8 eV. Recorded as is; not corrected.
- No MF6/MT2 in either file. Covariance files present: C-12 MF33; O-16 MF33 and MF34 (the suitability assessment states whether MF4/MT2 above 20/30 MeV carries any covariance).
- These are NEUTRON evaluations. Their use as a mirror for p+A elastic angular distributions (M3) is a scientific hypothesis to be assessed (isospin mirror validity, Coulomb-nuclear interference, isovector terms) before any plan; this record establishes identity, provenance and coverage only.

## Evaluation comment blocks (quoted for provenance)

- `n-006_C_012.endf`: MF1/MT451: 'ENDF/B-VII MOD 1 Evaluation, August 2015, G. M. Hale and P. G. Young (LANL)'; R-matrix analysis up to 6.45 MeV; 'At neutron energies above 6.45 MeV, the evaluation was joined to the existing ENDF/B-VII.0 evaluation for natural carbon, which is essentially version VI.1 (Fu, Axton, and Perey)'; 'All cross sections were kept the same as the VII.0 file in the energy range 6.47 <= E <= 20 MeV'; 'The version VI.3 extension to 150 MeV (Chadwick and Young) was not retained in this file.' (sic: the file nevertheless carries MF3/MT1-2 and tabulated MF4/MT2 data to 150 MeV; provenance of the 20-150 MeV part is NOT established by the comment block and is an open item of the suitability assessment)
- `n-008_O_016.endf`: MF1/MT451: 'EVAL-DEC16 Hale,Paris,Young,Chadwick'; R-matrix below 6.25 MeV, GNASH/optical-model evaluation 6.25-30 MeV; the previous MOD 2 (LA150-type, Chadwick et al.) evaluation retained from 30 to 150 MeV ('the previous evaluation, which extends to 150 MeV'); above 30 MeV MF6/MT5 carries the non-elastic emission and 'MF=4 MT=2 Elastic Angular Distributions' are listed among the non-zero data above 30 MeV; covariances 'continuing to 150 MeV'

**Status:** identity and provenance VERIFIED; energy coverage VERIFIED (<= 150 MeV); scientific suitability NOT YET ASSESSED; no construction use.
