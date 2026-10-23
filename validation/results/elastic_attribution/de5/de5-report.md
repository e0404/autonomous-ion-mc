Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E5 (table side): the pinned elastic table sampled in the D-E4 frame

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de5_sampling`

Input files (sha256):

- `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.npz`  1028bb67f57bc96edafe6d3dbcce532013edbbd31b9b80861239b09d9ad4202a
- `.ionmc-cache/ionmc-data/derived/elastic-proton-1fd24cff23d3f7e72aabcda2af8bc70ad47b81755f94a2ca55890cb6ae52c591.json`  ab5e2a357ab53060a84182607d93fe4373e50ce60622450af74728fd4e495e77

1e6 events per case from the P7 sampler (`ElasticTable.sample_mu`), seeds SeedSequence(20991000 + case), thin target at the beam energy, exact two-body kinematics, p-p cut at 16.26 deg CM, D-E4 kinetic-energy filter. Signature definitions: module docstring (plan 2.2). Thresholds theta_w use the Bethe range of water.
<!-- body below is rendered from the Git projection of the JSON record -->

## Per case and class

| case (row) | material | E | class | rate/primary | share | KE cut removed | theta_med,lab | P(th>th_2mm) | P(th>th_3mm) | T | <T_other> MeV | first min (deg; 3-sigma) | q-slope | Helm R_eff / s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 (14) | water | 30 | H-1 | 0.001177 | 0.143 | 0.00e+00 | 31.22 | 0.2165 | 0.0000 | 0.0000 | 8.111 | -; - | - | - |
| 0 (14) | water | 30 | O-16 | 0.007053 | 0.857 | 1.77e-04 | 24.07 | 0.1736 | 0.0969 | 0.0829 | 0.659 | 3.5; - | -2.51 | 3.70 / 0.00 |
| 1 (15) | water | 50 | H-1 | 0.0006654 | 0.109 | 0.00e+00 | 30.77 | 0.7010 | 0.5268 | 0.0000 | 13.321 | -; - | - | - |
| 1 (15) | water | 50 | O-16 | 0.005462 | 0.891 | 1.67e-04 | 20.11 | 0.3626 | 0.2250 | 0.0802 | 0.905 | 49.5; - | 3.65 | 3.35 / 0.00 |
| 2 (16) | water | 100 | H-1 | 0.001825 | 0.134 | 0.00e+00 | 29.98 | 0.9365 | 0.8832 | 0.0000 | 26.073 | -; - | - | - |
| 2 (16) | water | 100 | O-16 | 0.01177 | 0.866 | 1.41e-04 | 15.76 | 0.6229 | 0.4902 | 0.0667 | 1.365 | 40.5; - | -3.03 | 2.95 / 0.00 |
| 3 (17) | water | 150 | H-1 | 0.001756 | 0.249 | 0.00e+00 | 29.56 | 0.9842 | 0.9528 | 0.0000 | 38.712 | -; - | - | - |
| 3 (17) | water | 150 | O-16 | 0.005315 | 0.751 | 1.52e-04 | 12.45 | 0.6891 | 0.5722 | 0.0597 | 1.655 | 31.5; - | -2.31 | 2.95 / 0.00 |
| 4 (18) | water | 200 | H-1 | 0.001715 | 0.335 | 0.00e+00 | 29.26 | 1.0000 | 0.9878 | 0.0000 | 51.611 | -; - | - | - |
| 4 (18) | water | 200 | O-16 | 0.00341 | 0.665 | 1.56e-04 | 10.69 | 0.7348 | 0.6306 | 0.0561 | 1.978 | 29.5; - | -3.26 | 3.00 / 0.00 |
| 5 (19) | graphite | 30 | C-12 | 0.01745 | 1.000 | 1.13e-04 | 27.93 | 0.2617 | 0.1416 | 0.0973 | 1.093 | 3.5; - | -2.40 | 2.20 / 1.30 |
| 6 (20) | graphite | 50 | C-12 | 0.01243 | 1.000 | 1.17e-04 | 22.46 | 0.4420 | 0.2941 | 0.0742 | 1.296 | 5.5; - | -2.73 | 3.05 / 0.00 |
| 7 (21) | graphite | 100 | C-12 | 0.02537 | 1.000 | 7.70e-05 | 18.57 | 0.7100 | 0.5982 | 0.0843 | 2.322 | 44.5; - | 1.92 | 2.55 / 0.00 |
| 8 (22) | graphite | 150 | C-12 | 0.01253 | 1.000 | 8.20e-05 | 14.61 | 0.7631 | 0.6664 | 0.0668 | 2.598 | 39.5; - | -2.90 | 2.55 / 0.00 |
| 9 (23) | graphite | 200 | C-12 | 0.008423 | 1.000 | 7.80e-05 | 12.54 | 0.7998 | 0.7152 | 0.0603 | 2.995 | 32.5; - | -1.84 | 2.45 / 0.00 |

## Thresholds, p-p cut and the classification rule

Conflict tracks (both kinematic residuals < 1e-3; the p+H and p+A curves agree to 1e-3 below the crossover angle listed) are not excluded: p-p iff the event has a secondary proton, else p+A.

| case | theta_2mm | theta_3mm | class | theta_CM min (p-p) | crossover (deg) | raw conflict | misclassified | off curve | rate/primary after rule |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 39.11 | 48.40 | H-1 | 16.26 | 1.86 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.001177 |
| 0 | 39.11 | 48.40 | O-16 | - | 1.86 | 3.87e-03 | 0.00e+00 | 0.00e+00 | 0.007052 |
| 1 | 24.38 | 29.98 | H-1 | 16.26 | 1.85 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.0006654 |
| 1 | 24.38 | 29.98 | O-16 | - | 1.85 | 5.68e-03 | 0.00e+00 | 0.00e+00 | 0.005461 |
| 2 | 13.03 | 15.98 | H-1 | 16.27 | 1.83 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.001825 |
| 2 | 13.03 | 15.98 | O-16 | - | 1.83 | 9.10e-03 | 0.00e+00 | 0.00e+00 | 0.01176 |
| 3 | 9.11 | 11.16 | H-1 | 16.26 | 1.81 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.001756 |
| 3 | 9.11 | 11.16 | O-16 | - | 1.81 | 1.42e-02 | 0.00e+00 | 0.00e+00 | 0.005314 |
| 4 | 7.10 | 8.70 | H-1 | 16.26 | 1.78 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.001715 |
| 4 | 7.10 | 8.70 | O-16 | - | 1.78 | 1.90e-02 | 0.00e+00 | 0.00e+00 | 0.003409 |
| 5 | 39.11 | 48.40 | C-12 | - | 1.88 | 2.96e-03 | 0.00e+00 | 0.00e+00 | 0.01745 |
| 6 | 24.38 | 29.98 | C-12 | - | 1.87 | 4.75e-03 | 0.00e+00 | 0.00e+00 | 0.01243 |
| 7 | 13.03 | 15.98 | C-12 | - | 1.85 | 6.80e-03 | 0.00e+00 | 0.00e+00 | 0.02537 |
| 8 | 9.11 | 11.16 | C-12 | - | 1.83 | 1.07e-02 | 0.00e+00 | 0.00e+00 | 0.01253 |
| 9 | 7.10 | 8.70 | C-12 | - | 1.80 | 1.40e-02 | 0.00e+00 | 0.00e+00 | 0.008423 |

The first-minimum, q-slope and Helm descriptors are recorded, not interpreted: the frozen 1-degree estimator picks noise dips in a 256-bin piecewise-uniform table, and the 3-sigma variant misses the J1 zero.

Absolute cross sections and angle / energy histograms are not committed (decision 0041); they are in the generated record. Geant4 ratios (ionmc / Geant4) are produced later by `signature_ratios` from the D-E4 summary.
