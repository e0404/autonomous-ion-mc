Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E4 / D-E5: table / Geant4 signature ratios

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de4_parse --ratios`

Input files (sha256):

- `validation/generated/elastic_attribution/de4/de4.json`  168f232539352f7e75d4db48ae111418fef284b5ae34ae4bff29e77b9eaa3d05
- `validation/generated/elastic_attribution/de5/de5.json`  969c30d8dfa9110f98d5cbe773e88606f6f9bc816c44c8b678b3c166db662594

<!-- body below is rendered from the Git projection of the JSON record -->

ratio = ionmc table / Geant4; rates per primary after the frozen KE filter, off-curve exclusion and classification

## Case 0

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| H-1 | 24103 | 0.977 +- 0.006 | 1.068 +- 0.014 | 0.000 +- 0.000 | 1.026 | - | 1.047 | - | - | - |
| O-16 | 140322 | 1.005 +- 0.003 | 1.568 +- 0.012 | 1.341 +- 0.013 | 1.069 | 2.407 | 1.357 | 0.636 | - | 0.957 |

## Case 1

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| H-1 | 13351 | 0.997 +- 0.009 | 1.032 +- 0.006 | 1.025 +- 0.009 | 1.011 | - | 1.021 | - | - | - |
| O-16 | 108606 | 1.006 +- 0.003 | 1.670 +- 0.010 | 2.277 +- 0.021 | 1.196 | 7.207 | 2.130 | 14.143 | - | -1.103 |

## Case 2

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| H-1 | 7078 | 1.032 +- 0.012 | 1.008 +- 0.003 | 1.009 +- 0.005 | 1.013 | - | 1.013 | - | - | - |
| O-16 | 46790 | 1.006 +- 0.005 | 1.541 +- 0.009 | 1.904 +- 0.015 | 1.378 | 14.196 | 3.266 | 1.373 | - | 0.620 |

## Case 3

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| H-1 | 6028 | 1.165 +- 0.015 | 0.996 +- 0.001 | 0.995 +- 0.003 | 0.987 | - | 0.991 | - | - | - |
| O-16 | 20974 | 1.013 +- 0.007 | 1.343 +- 0.009 | 1.557 +- 0.014 | 1.338 | 16.409 | 3.810 | 1.105 | - | 0.424 |

## Case 4

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| H-1 | 5919 | 1.159 +- 0.015 | 1.000 +- 0.000 | 1.002 +- 0.002 | 1.020 | - | 1.022 | - | - | - |
| O-16 | 13464 | 1.013 +- 0.009 | 1.287 +- 0.010 | 1.453 +- 0.014 | 1.348 | 24.965 | 4.602 | 1.157 | - | - |

## Case 5

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| C-12 | 340490 | 1.025 +- 0.002 | 3.373 +- 0.020 | 3.009 +- 0.023 | 1.345 | 3.136 | 2.002 | 0.636 | - | 0.796 |

## Case 6

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| C-12 | 245121 | 1.014 +- 0.002 | 2.322 +- 0.010 | 3.648 +- 0.025 | 1.396 | 7.668 | 2.519 | 0.647 | - | 0.626 |

## Case 7

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| C-12 | 98446 | 1.031 +- 0.003 | 1.722 +- 0.007 | 2.277 +- 0.012 | 1.608 | 18.718 | 4.151 | 1.219 | - | -0.434 |

## Case 8

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| C-12 | 49206 | 1.019 +- 0.005 | 1.450 +- 0.006 | 1.739 +- 0.010 | 1.540 | 21.840 | 4.409 | 1.339 | - | 0.484 |

## Case 9

| class | Geant4 events | rate_per_primary_after_rule | P_theta_lab_gt_theta_2mm | P_theta_lab_gt_theta_3mm | theta_lab_median_deg | T_tail_moment | mean_t_other_mev | first_minimum_deg | first_minimum_3sigma_deg | qslope |
|---|---|---|---|---|---|---|---|---|---|---|
| C-12 | 33435 | 1.008 +- 0.006 | 1.335 +- 0.006 | 1.546 +- 0.009 | 1.523 | 26.314 | 4.926 | 1.140 | - | 0.241 |
