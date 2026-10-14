Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E1: pinned p+A shape and shape families against the V10-A data

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de1_shapes`

Input files (sha256):

- `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.npz`  1028bb67f57bc96edafe6d3dbcce532013edbbd31b9b80861239b09d9ad4202a
- `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.json`  ab5e2a357ab53060a84182607d93fe4373e50ce60622450af74728fd4e495e77
- `validation/generated/RUN-20261012T142343Z-081123c3/transport/lv5c-dd976c5-03/03-v10-a.txt`  05499b3a30a8695544f153dea9cb7367c744c2ea43e3e4e225a2417eb7397cf9

Model = pinned table (analytic black disk of the table radius, plus unscreened Rutherford, incoherent; the V10-A model value). Families: Amendment 19 (b) at their a-priori parameters (no fit). Partial integrals: trapezoid over the data nodes from 10 deg to the window end (10-30, 10-40 deg, and the last data angle), nuclear + Rutherford. The R fit is **data-seen**, a diagnostic, never a model parameter.


## Minima: measured against the table

| dataset | target | E | theta_min,1 (V10-A rule) | significant? | table first min: bin-density / 1-deg bins / analytic disk (total) / J1 zero | Ru share at data min | R_table | R implied by V10-A min (3.8317/q) | significant minima (deg, merged duplicates, >=10 % and >=3 sigma below both neighbours) |
|---|---|---|---|---|---|---|---|---|---|
| C0063002@135 | O-16 | 135 | 17.17 | False | 35.2 / 33.5 / 34.1 / 34.0 | 0.045 | 2.66 | 5.21 | none |
| C0141002@179.9 | O-16 | 179.9 | 14.99 | False | 29.2 / 26.5 / 28.9 / 28.8 | 0.059 | 2.68 | 5.12 | none |
| C0148002@201.4 | O-16 | 201.4 | - | - | 28.3 / 26.5 / 27.2 / 27.1 | - | 2.69 | - | none |
| C0550002@200 | O-16 | 200 | 79.15 | True | 28.5 / 27.5 / 27.3 / 27.2 | 0.213 | 2.69 | 0.99 | 79.2 |
| O0274009@96 | C-12 | 96 | 10.86 | False | 52.0 / 44.5 / 50.8 / 50.7 | 0.173 | 2.22 | 10.02 | none |
| O0553008@75 | C-12 | 75 | 11.89 | False | 54.3 / 51.5 / 54.4 / 54.3 | 0.153 | 2.36 | 10.40 | none |
| O0553002@150 | C-12 | 150 | - | - | 40.0 / 36.5 / 38.8 / 38.7 | - | 2.27 | - | none |
| C0055002@122 | C-12 | 122 | 138.80 | True | 44.0 / 42.5 / 44.4 / 44.3 | 0.012 | 2.22 | 0.89 | 138.8 |
| C0055002@160 | C-12 | 160 | 105.40 | False | 38.2 / 34.5 / 37.5 / 37.4 | 0.016 | 2.27 | 0.91 | 138.9 |
| C0057002@200 | C-12 | 200 | 93.77 | True | 34.1 / 31.5 / 33.0 / 32.9 | 0.066 | 2.28 | 0.88 | 93.8 |
| C1420002@250 | C-12 | 250 | 31.56 | False | 29.2 / 27.5 / 29.0 / 28.9 | 0.127 | 2.29 | 2.11 | 56.5, 58.0, 60.0, 85.8 |

**Frozen criterion** (plan: dip in >= 2 adjacent distinct-angle points, deeper than 2 sigma_stat, model Rutherford share < 20 %), applied to every dataset; the `significant` column above is my stricter (>= 10 % and >= 3 sigma) cross-check.

| dataset | V10-A rule theta_min,1 | drop vs lower neighbour (sigma, duplicates not merged) | V10-A cut 0.8 theta_min | frozen criterion met at that minimum | frozen-criterion minima anywhere (deg; depth sigma; Ru share) |
|---|---|---|---|---|---|
| C0063002@135 | 17.17 | 0.74 | 13.74 | False | none |
| C0141002@179.9 | 14.99 | 0.48 | 11.99 | False | none |
| C0148002@201.4 | - | - | - | - | none |
| C0550002@200 | 79.15 | 11.27 | 63.32 | False | 79.2 (10.7; 0.21; not met) |
| O0274009@96 | 10.86 | - | 8.69 | False | 10.9 (-; 0.17; not met) |
| O0553008@75 | 11.89 | - | 9.51 | False | none |
| O0553002@150 | - | - | - | - | none |
| C0055002@122 | 138.80 | 3.46 | 111.04 | True | 138.8 (10.4; 0.01; met) |
| C0055002@160 | 105.40 | 0.21 | 84.32 | True | 105.4 (2.2; 0.02; met); 138.9 (5.3; 0.01; met) |
| C0057002@200 | 93.77 | 0.61 | 75.02 | True | 93.8 (6.1; 0.07; met); 129.6 (3.4; 0.01; met) |
| C1420002@250 | 31.56 | 0.34 | 25.25 | False | 31.6 (0.8; 0.13; not met); 34.8 (0.4; 0.03; not met); 35.9 (1.5; 0.02; not met); 60.0 (18.9; 0.02; met) |

**Finding.** The V10-A minima of C0063 @ 135 MeV (17.17 deg) and C0141 @ 179.9 MeV (14.99 deg) are duplicate-angle artefacts (a 0.74 sigma / 0.48 sigma wiggle between two runs at the same angle), not diffraction minima: the data fall monotonically over 10-38 deg. The recorded V10-A evaluation used them to set the (a) cut 0.8 theta_min,1 = 13.74 deg (C0063; 11.99 deg for C0141), which excluded the 15 and 20 deg (a) points of C0063 and C0141 (`03-v10-a.txt`: "theta >= 0.8 theta_min,1 (measured: 13.74 deg)"). This is a defect of the V10-A minimum rule (duplicate angles not merged), recorded here as a finding; the V10-A verdict is unchanged and nothing is re-evaluated. The statement 'first minima about 2x too wide / R about 4.9 fm' rested on these minima and is withdrawn.


## Partial sigma_tr, window 10-30 deg (model / data; families at a-priori parameters)

| dataset | range (deg) | data sigma_tr (mb) | model | disk | Rx1.2 | s0.5 | s0.9 | gauss | model sigma/data sigma | data partial sigma / table sigma_el | Ru-only / data |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C0063002@135 | 10-30 | 2.64 | 2.103 | 2.103 | 1.419 | 2.125 | 1.857 | 2.685 | 1.667 | 0.38 | 0.120 |
| C0141002@179.9 | 10-30 | 1.54 | 1.477 | 1.477 | 1.109 | 1.489 | 1.286 | 2.159 | 1.375 | 0.40 | 0.119 |
| C0148002@201.4 | 10-29 | 0.995 | 1.749 | 1.749 | 1.305 | 1.753 | 1.491 | 2.646 | 1.597 | 0.31 | 0.144 |
| C0550002@200 | 10-28 | 1.14 | 1.526 | 1.526 | 1.080 | 1.533 | 1.307 | 2.275 | 1.401 | 0.36 | 0.124 |
| O0274009@96 | 10-28 | 7 | 1.271 | 1.271 | 1.259 | 1.368 | 1.341 | 1.253 | 1.097 | 0.54 | 0.047 |
| O0553008@75 | 10-28 | 60.2 | 0.218 | 0.218 | 0.222 | 0.234 | 0.238 | 0.211 | 0.187 | 3.10 | 0.009 |
| O0553002@150 | 10-27 | 1.59 | 2.357 | 2.357 | 1.873 | 2.419 | 2.119 | 2.571 | 1.847 | 0.34 | 0.087 |
| C0055002@122 | 10-29 | 3.15 | 1.939 | 1.939 | 1.684 | 2.024 | 1.843 | 2.015 | 1.498 | 0.43 | 0.068 |
| C0055002@160 | 10-27 | 1.54 | 2.161 | 2.161 | 1.659 | 2.207 | 1.903 | 2.415 | 1.689 | 0.37 | 0.080 |
| C0057002@200 | 10-29 | 1.06 | 2.023 | 2.023 | 1.348 | 2.009 | 1.636 | 2.568 | 1.601 | 0.37 | 0.080 |
| C1420002@250 | 10-30 | 1.28 | 1.086 | 1.086 | 0.769 | 1.072 | 0.850 | 1.626 | 0.938 | 0.56 | 0.045 |

## Partial sigma_tr, window 10-40 deg (model / data; families at a-priori parameters)

| dataset | range (deg) | data sigma_tr (mb) | model | disk | Rx1.2 | s0.5 | s0.9 | gauss | model sigma/data sigma | data partial sigma / table sigma_el | Ru-only / data |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C0063002@135 | 10-39 | 2.92 | 2.013 | 2.013 | 1.710 | 2.006 | 1.726 | 2.982 | 1.655 | 0.39 | 0.133 |
| C0141002@179.9 | 10-39 | 1.7 | 1.802 | 1.802 | 1.635 | 1.609 | 1.242 | 2.197 | 1.433 | 0.41 | 0.135 |
| C0148002@201.4 | 10-40 | 1.11 | 2.468 | 2.468 | 2.009 | 2.041 | 1.453 | 2.623 | 1.728 | 0.32 | 0.168 |
| C0550002@200 | 10-39 | 1.27 | 2.020 | 2.020 | 1.750 | 1.736 | 1.276 | 2.320 | 1.493 | 0.37 | 0.146 |
| O0274009@96 | 10-39 | 8.6 | 1.814 | 1.814 | 1.380 | 1.816 | 1.557 | 2.024 | 1.297 | 0.59 | 0.050 |
| O0553008@75 | 10-37 | 75.6 | 0.308 | 0.308 | 0.253 | 0.313 | 0.285 | 0.324 | 0.222 | 3.37 | 0.009 |
| O0553002@150 | 10-38 | 1.82 | 2.566 | 2.566 | 1.827 | 2.527 | 2.072 | 3.561 | 1.936 | 0.35 | 0.101 |
| C0055002@122 | 10-39 | 3.54 | 2.417 | 2.417 | 1.666 | 2.377 | 1.958 | 3.049 | 1.654 | 0.45 | 0.078 |
| C0055002@160 | 10-37 | 1.73 | 2.290 | 2.290 | 1.640 | 2.261 | 1.848 | 3.210 | 1.748 | 0.38 | 0.092 |
| C0057002@200 | 10-39 | 1.13 | 2.087 | 2.087 | 1.916 | 1.996 | 1.573 | 3.067 | 1.616 | 0.38 | 0.097 |
| C1420002@250 | 10-39 | 1.36 | 1.358 | 1.358 | 1.246 | 1.154 | 0.822 | 1.695 | 0.986 | 0.56 | 0.053 |

## Partial sigma_tr, window 10 deg to the last data angle (model / data; families at a-priori parameters)

| dataset | range (deg) | data sigma_tr (mb) | model | disk | Rx1.2 | s0.5 | s0.9 | gauss | model sigma/data sigma | data partial sigma / table sigma_el | Ru-only / data |
|---|---|---|---|---|---|---|---|---|---|---|---|
| C0063002@135 | 10-92 | 3.23 | 4.767 | 4.767 | 3.579 | 2.646 | 1.718 | 2.922 | 1.895 | 0.39 | 0.190 |
| C0141002@179.9 | 10-90 | 1.85 | 3.698 | 3.698 | 3.074 | 1.944 | 1.242 | 2.113 | 1.605 | 0.41 | 0.193 |
| C0148002@201.4 | 10-59 | 1.16 | 3.214 | 3.214 | 2.890 | 2.228 | 1.448 | 2.559 | 1.817 | 0.32 | 0.204 |
| C0550002@200 | 10-126 | 1.34 | 5.368 | 5.368 | 4.384 | 2.113 | 1.332 | 2.310 | 1.706 | 0.37 | 0.241 |
| O0274009@96 | 10-95 | 9.91 | 3.273 | 3.273 | 2.410 | 2.173 | 1.466 | 2.768 | 1.447 | 0.60 | 0.070 |
| O0553008@75 | 10-44 | 79.7 | 0.346 | 0.346 | 0.250 | 0.342 | 0.296 | 0.396 | 0.233 | 3.41 | 0.009 |
| O0553002@150 | 10-74 | 2 | 4.461 | 4.461 | 3.458 | 2.988 | 1.980 | 3.773 | 2.137 | 0.36 | 0.134 |
| C0055002@122 | 10-156 | 3.9 | 6.054 | 6.054 | 4.961 | 2.850 | 1.879 | 3.557 | 1.880 | 0.45 | 0.127 |
| C0055002@160 | 10-156 | 1.9 | 6.640 | 6.640 | 5.323 | 2.774 | 1.803 | 3.427 | 2.028 | 0.39 | 0.157 |
| C0057002@200 | 10-157 | 1.21 | 6.774 | 6.774 | 5.644 | 2.454 | 1.571 | 3.043 | 1.903 | 0.38 | 0.161 |
| C1420002@250 | 10-102 | 1.42 | 3.380 | 3.380 | 2.640 | 1.336 | 0.826 | 1.671 | 1.127 | 0.57 | 0.082 |

The `full` window is dominated by the model's high-angle diffraction lobes (model 100-600x above the data at 60-90 deg for the O sets), so the 10-30 / 10-40 deg windows are the like-for-like comparison. Datasets with `data partial sigma / table sigma_el` above 1 (the data integral alone exceeds the table sigma_el) are not credible against the table: **O0553008 @ 75 MeV (p+C, Rolland 1966)** has 5800 mb/sr at 10 deg against a model 694 mb/sr; its EXFOR subentry gives MB/SR with a 5 % digital error and no other factor, so the cause is not an ingestion unit error that can be identified here (it is a source/transcription or frame issue of that subentry); it is reported but excluded from summaries and nothing was changed.


## Model / data - 1 at the V10-A angles (nuclear + Rutherford)

| dataset | 10 deg | 15 deg | 20 deg | 30 deg |
|---|---|---|---|---|
| C0063002@135 | +0.12 @10.8 | +0.45 @15.1 | +1.48 @19.4 | +1.11 @30.0 |
| C0141002@179.9 | +0.14 @10.4 | +0.37 @15.0 | +1.23 @20.7 | -0.59 @29.8 |
| C0148002@201.4 | +0.41 @9.8 | +0.73 @14.7 | +1.93 @19.6 | +1.66 @30.3 |
| C0550002@200 | +0.07 @10.8 | +0.50 @15.1 | +1.03 @19.3 | +0.92 @30.0 |
| O0274009@96 | +0.05 @9.8 | -0.17 @14.8 | +0.19 @19.4 | +2.72 @30.5 |
| O0553008@75 | -0.88 @10.2 | -0.84 @15.2 | -0.82 @19.1 | - |
| O0553002@150 | +0.19 @10.8 | +0.81 @16.4 | +1.43 @19.2 | +4.42 @30.1 |
| C0055002@122 | -0.03 @10.1 | +0.19 @15.5 | +0.81 @19.9 | +6.39 @30.7 |
| C0055002@160 | +0.08 @11.0 | +0.32 @14.3 | +1.94 @20.8 | +2.99 @30.6 |
| C0057002@200 | +0.09 @8.9 | +0.77 @15.5 | +1.76 @19.9 | +0.23 @30.8 |
| C1420002@250 | -0.34 @9.9 | +0.02 @15.4 | +0.80 @20.0 | -0.58 @29.8 |

## One-parameter R diagnostic (DATA-SEEN; never a model parameter)

Black-disk shape with sigma_el and the Rutherford term fixed; least squares in log(dsigma) on the stated CM window (R scanned 1-8 fm in 0.01 fm steps). The windows 10-30 and 10-40 deg are fixed a priori; `to_V10A_first_minimum` ends at the V10-A rule minimum and is only as meaningful as that minimum.

| dataset | window | R_fit (fm) | R_table (fm) | R_fit/R_table | rms log resid at R_fit | at R_table |
|---|---|---|---|---|---|---|
| C0063002@135 | 10-30deg | 3.27 | 2.66 | 1.23 | 0.475 | 1.025 |
| C0063002@135 | 10-40deg | 2.94 | 2.66 | 1.11 | 0.940 | 1.035 |
| C0063002@135 | to_V10A_first_minimum | 3.87 | 2.66 | 1.46 | 0.210 | 0.354 |
| C0141002@179.9 | 10-30deg | 2.79 | 2.68 | 1.04 | 0.265 | 0.450 |
| C0141002@179.9 | 10-40deg | 2.41 | 2.68 | 0.90 | 0.738 | 0.882 |
| C0141002@179.9 | to_V10A_first_minimum | 3.48 | 2.68 | 1.30 | 0.130 | 0.219 |
| C0148002@201.4 | 10-30deg | 2.97 | 2.69 | 1.11 | 0.577 | 0.737 |
| C0148002@201.4 | 10-40deg | 7.79 | 2.69 | 2.90 | 0.911 | 1.534 |
| C0148002@201.4 | to_V10A_first_minimum | 2.97 | 2.69 | 1.11 | 0.577 | 0.737 |
| C0550002@200 | 10-30deg | 2.91 | 2.69 | 1.08 | 0.488 | 0.603 |
| C0550002@200 | 10-40deg | 8.00 | 2.69 | 2.98 | 1.050 | 1.333 |
| C0550002@200 | to_V10A_first_minimum | 7.60 | 2.69 | 2.83 | 3.165 | 3.636 |
| O0274009@96 | 10-30deg | 3.08 | 2.22 | 1.39 | 0.161 | 0.405 |
| O0274009@96 | 10-40deg | 2.71 | 2.22 | 1.22 | 0.354 | 0.758 |
| O0274009@96 | to_V10A_first_minimum | - | - | - | - | - |
| O0553008@75 | 10-30deg | 3.09 | 2.36 | 1.31 | 1.564 | 1.710 |
| O0553008@75 | 10-40deg | 2.75 | 2.36 | 1.16 | 1.500 | 1.565 |
| O0553008@75 | to_V10A_first_minimum | - | - | - | - | - |
| O0553002@150 | 10-30deg | 2.92 | 2.27 | 1.29 | 0.490 | 1.106 |
| O0553002@150 | 10-40deg | 2.56 | 2.27 | 1.13 | 0.838 | 1.155 |
| O0553002@150 | to_V10A_first_minimum | 2.92 | 2.27 | 1.29 | 0.490 | 1.106 |
| C0055002@122 | 10-30deg | 3.06 | 2.22 | 1.38 | 0.323 | 0.843 |
| C0055002@122 | 10-40deg | 2.85 | 2.22 | 1.28 | 0.627 | 1.148 |
| C0055002@122 | to_V10A_first_minimum | 7.53 | 2.22 | 3.39 | 3.676 | 4.639 |
| C0055002@160 | 10-30deg | 2.81 | 2.27 | 1.24 | 0.458 | 1.037 |
| C0055002@160 | 10-40deg | 2.59 | 2.27 | 1.14 | 0.795 | 1.059 |
| C0055002@160 | to_V10A_first_minimum | 7.76 | 2.27 | 3.42 | 3.702 | 4.832 |
| C0057002@200 | 10-30deg | 2.71 | 2.28 | 1.19 | 0.585 | 0.956 |
| C0057002@200 | 10-40deg | 2.20 | 2.28 | 0.97 | 1.106 | 1.149 |
| C0057002@200 | to_V10A_first_minimum | 7.43 | 2.28 | 3.26 | 3.259 | 4.065 |
| C1420002@250 | 10-30deg | 2.38 | 2.29 | 1.04 | 0.519 | 0.564 |
| C1420002@250 | 10-40deg | 2.06 | 2.29 | 0.90 | 0.982 | 1.213 |
| C1420002@250 | to_V10A_first_minimum | 2.31 | 2.29 | 1.01 | 0.554 | 0.557 |
