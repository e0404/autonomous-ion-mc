Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E7: official M2 / M2' builder records, Git-safe projection

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de7_m2prime_record --m2 DIR --m2prime DIR --scratch-text FILE`

Input files (sha256):

- `/home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/-home-wahln-aiprojects-ion-mc/3057d970-cf12-4711-8c0a-ae0cbe21cb5d/scratchpad/m2-official/m2/elastic-proton-ad0d77578f29cde5dd9c4747d847cfcc586ec8bfb0643544aaa8addcfeb82034.json`  10afbe559a65a85249b231997bf81e3e1ff31272565eed8854d9513d12df1b3e
- `/home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/-home-wahln-aiprojects-ion-mc/3057d970-cf12-4711-8c0a-ae0cbe21cb5d/scratchpad/m2-official/m2prime/elastic-proton-addd817961ca24ebde185698ece4b38440a16cfdcf3b21dfad9416a089519b47.json`  18434dace7271bf9ba0e1c12842784107e173fe9905bb36ee0286611b7d50e72
- `/home/wahln/.local/share/ionmc-experiment/experiment-v3/tmp/claude-1000/-home-wahln-aiprojects-ion-mc/3057d970-cf12-4711-8c0a-ae0cbe21cb5d/scratchpad/m2-official/rb/elastic-proton-2a14cdcb93db416f30db0b0fb3bf1607af414590fadd92ea0fe3c06f802ec221.json`  65ddb9f647c0af85d8673542219b9e989786f571076545401312d33e35aec6d7

Relative differences against the researcher's scratch numbers only; no absolute cross section, no measured value, no % deviation. sigma_R / Barashenkov is withheld for M2 (Amendment 21 (f)).
<!-- body below is rendered from the Git projection of the JSON record -->

Builder records of the official M2 and M2' constructions (Amendments 20 (b''), 21). Ratios, pass/fail flags, counts and hashes only.

## M2

Table id `ad0d77578f29cde5dd9c4747d847cfcc586ec8bfb0643544aaa8addcfeb82034`; model `m2-optical-glauber`, variant `V-b`, Pauli variant `None`, Pauli E_min - MeV, beta `unitarity`, pp level `none`, unfold `rms-matched`, inputs unpinned: False.

sigma_R / Barashenkov: withheld for M2 (Amendment 21 (f)): the ratios with the Amendment 20 M2 percentages would recover the measured sigma_tot

| target | E (MeV) | sigma_el / Barashenkov | sigma_tot in +-10 % band | q_min (1/fm) | q_min / q_H | Pauli factor (peak) |
|---|---|---|---|---|---|---|
| Al-27 | 100 | 0.761 | - | 1.068 | - | - |
| Al-27 | 150 | 1.265 | - | 1.193 | - | - |
| Al-27 | 200 | 1.669 | - | 1.225 | - | - |
| C-12 | 100 | 0.905 | pass | - | - | - |
| C-12 | 150 | 1.336 | fail | - | - | - |
| C-12 | 200 | 1.454 | fail | 1.596 | 0.835 | - |
| Ca-40 | 100 | 0.788 | - | 0.885 | - | - |
| Ca-40 | 150 | 1.264 | - | 0.968 | - | - |
| Ca-40 | 200 | 1.582 | - | 1.004 | - | - |
| N-14 | 100 | 0.930 | - | 1.383 | - | - |
| N-14 | 150 | 1.479 | - | 1.487 | - | - |
| N-14 | 200 | 1.442 | - | 1.499 | - | - |
| O-16 | 100 | 0.851 | pass | 1.285 | 0.774 | - |
| O-16 | 150 | 1.407 | fail | 1.379 | 0.830 | - |
| O-16 | 200 | 1.622 | fail | 1.397 | 0.841 | - |
| P-31 | 100 | 0.772 | - | 1.000 | - | - |
| P-31 | 150 | 1.315 | - | 1.106 | - | - |
| P-31 | 200 | 1.677 | - | 1.143 | - | - |
| Si-28 | 100 | 0.782 | - | 1.044 | - | - |
| Si-28 | 150 | 1.345 | - | 1.157 | - | - |
| Si-28 | 200 | 1.710 | - | 1.186 | - | - |

sigma_tot against the measured values, pass/fail per cell (band +-10 %):

- Al-27: at least 2 of 3: -; all 3: -
- C-12: at least 2 of 3: fail; all 3: fail
- Ca-40: at least 2 of 3: -; all 3: -
- N-14: at least 2 of 3: -; all 3: -
- O-16: at least 2 of 3: fail; all 3: fail
- P-31: at least 2 of 3: -; all 3: -
- Si-28: at least 2 of 3: -; all 3: -

Trigger 2: max 0.000244 against 0.001, passed: True, nodes evaluated 68. Component maxima: b_step 0.000037, q_cutoff 0.000000, q_step 0.000244, radial_step 0.000000, s_grid 0.000002.

Trigger 3 (>= 2 of 3 in [0.7, 1.4], C-12 and O-16): **fails**.

| target | E (MeV) | ratio | in band | margin to nearest bound | excess over upper bound |
|---|---|---|---|---|---|
| C-12 | 100 | 0.9054 | pass | 0.2054 | -0.3533 |
| C-12 | 150 | 1.3365 | pass | 0.0635 | -0.0454 |
| C-12 | 200 | 1.4543 | fail | -0.0543 | 0.0388 |
| C-12 | verdict | 2 of 3 in band | pass | | |
| O-16 | 100 | 0.8508 | pass | 0.1508 | -0.3923 |
| O-16 | 150 | 1.4072 | fail | -0.0072 | 0.0052 |
| O-16 | 200 | 1.6218 | fail | -0.2218 | 0.1584 |
| O-16 | verdict | 1 of 3 in band | fail | | |

Trigger 4 (q_min / q_H in [0.6, 1.0]): **pass**.

- C-12: q_H 1.912 1/fm; 150 MeV q_min not evaluable (no true minimum), q_min/q_H -, in band -; 200 MeV q_min 1.596, q_min/q_H 0.835, in band pass; verdict cell 200_mev: pass
- O-16: q_H 1.661 1/fm; 150 MeV q_min 1.379, q_min/q_H 0.830, in band pass; 200 MeV q_min 1.397, q_min/q_H 0.841, in band pass; verdict cell 150_mev: pass

Pins:
- inputs JSON sha256 `d4ac274f7b8aaea6df915a06e0735e1b73ea8075c072de22f6badd0235b52747`; Horiuchi JSON sha256 `None`; inputs unpinned: False
- SAID NP_ALFA: `47404925af3cbcf8d98613ef06dc72ca254ef667fc75510bc5d4c427acf72685`, `db95e4dcdf76e12e4fdb480981ef4c353c232d84a23dc2de47b4453f17d13136`, `14399f8deb631857b389f3399b9d7df08b789ebf6fedc7014c6151dc09d960c4`
- SAID NP_SGT: `eec6f944a1a9b85800fbeabd07b716f697a2db0edb5abf5852b7221614652290`, `d3f35fe2b873dc39f8c3c6f9de1efdb916aaec356221e23d0c5debe65a1f120d`, `d4d955c71d4a4a3a198d968e884ff62534b9f4381345ef514be32aa669cd410a`
- SAID PP_ALFA: `b2ca5813878a1c7a9c7e1ce808e8971dbbafaa0e7375240b2d58f6c60e53e252`, `f7eeda20f311db64f38c6e6b09d7b97385bc2ae009cdea0ccceecb9a3d83c7fd`, `4bae553e8a06e82f02b8032ca648e50631b5cee6f8437c316061e577701d8c8a`
- SAID PP_SGT: `508d47be0adcc78aa5546d50d60cec8410d6b85a6c94e8860368ca173a8b1b37`, `b72360a2c3813bd523a57e529ba964b1c30a806c02963b6a941efbd7e7a6bcb8`, `d93b3f4d0c7ae46b941c9ae82292238384624af416d1acccbf4e2146f478954c`
- DSG files: none

## M2prime

Table id `addd817961ca24ebde185698ece4b38440a16cfdcf3b21dfad9416a089519b47`; model `m2prime-pauli-blocked`, variant `M2prime`, Pauli variant `clementel-villi`, Pauli E_min 80.0 MeV, beta `unitarity`, pp level `none`, unfold `rms-matched`, inputs unpinned: False.

| target | E (MeV) | sigma_el / Barashenkov | sigma_R / Barashenkov | sigma_tot in +-10 % band | q_min (1/fm) | q_min / q_H | Pauli factor (peak) |
|---|---|---|---|---|---|---|---|
| Al-27 | 100 | 0.593 | 1.110 | - | 1.113 | - | 0.364 |
| Al-27 | 150 | 1.017 | 1.002 | - | 1.198 | - | 0.576 |
| Al-27 | 200 | 1.370 | 0.945 | - | 1.221 | - | 0.682 |
| C-12 | 100 | 0.595 | 1.091 | fail | 1.443 | 0.755 | 0.503 |
| C-12 | 150 | 0.976 | 1.011 | pass | 1.538 | 0.805 | 0.669 |
| C-12 | 200 | 1.129 | 0.954 | pass | 1.567 | 0.819 | 0.752 |
| Ca-40 | 100 | 0.660 | 1.076 | - | 0.932 | - | 0.421 |
| Ca-40 | 150 | 1.061 | 1.028 | - | 0.986 | - | 0.614 |
| Ca-40 | 200 | 1.337 | 0.986 | - | 1.009 | - | 0.711 |
| N-14 | 100 | 0.628 | 1.006 | - | 1.367 | - | 0.546 |
| N-14 | 150 | 1.095 | 0.954 | - | 1.453 | - | 0.697 |
| N-14 | 200 | 1.128 | 0.923 | - | 1.480 | - | 0.773 |
| O-16 | 100 | 0.596 | 1.062 | fail | 1.282 | 0.772 | 0.599 |
| O-16 | 150 | 1.062 | 0.986 | pass | 1.356 | 0.816 | 0.732 |
| O-16 | 200 | 1.286 | 0.938 | pass | 1.382 | 0.832 | 0.799 |
| P-31 | 100 | 0.619 | 1.080 | - | 1.049 | - | 0.388 |
| P-31 | 150 | 1.075 | 0.991 | - | 1.120 | - | 0.592 |
| P-31 | 200 | 1.391 | 0.946 | - | 1.145 | - | 0.694 |
| Si-28 | 100 | 0.612 | 1.093 | - | 1.087 | - | 0.399 |
| Si-28 | 150 | 1.083 | 0.989 | - | 1.161 | - | 0.599 |
| Si-28 | 200 | 1.407 | 0.942 | - | 1.183 | - | 0.699 |

sigma_tot against the measured values, pass/fail per cell (band +-10 %):

- Al-27: at least 2 of 3: -; all 3: -
- C-12: at least 2 of 3: pass; all 3: fail
- Ca-40: at least 2 of 3: -; all 3: -
- N-14: at least 2 of 3: -; all 3: -
- O-16: at least 2 of 3: pass; all 3: fail
- P-31: at least 2 of 3: -; all 3: -
- Si-28: at least 2 of 3: -; all 3: -

Trigger 2: max 0.000038 against 0.001, passed: True, nodes evaluated 33. Component maxima: b_step 0.000038, q_cutoff 0.000000, q_step 0.000033, radial_step 0.000000, s_grid 0.000002.

Trigger 3 (>= 2 of 3 in [0.7, 1.4], C-12 and O-16): **passes**.

| target | E (MeV) | ratio | in band | margin to nearest bound | excess over upper bound |
|---|---|---|---|---|---|
| C-12 | 100 | 0.5955 | fail | -0.1045 | -0.5747 |
| C-12 | 150 | 0.9761 | pass | 0.2761 | -0.3028 |
| C-12 | 200 | 1.1288 | pass | 0.2712 | -0.1937 |
| C-12 | verdict | 2 of 3 in band | pass | | |
| O-16 | 100 | 0.5955 | fail | -0.1045 | -0.5746 |
| O-16 | 150 | 1.0622 | pass | 0.3378 | -0.2413 |
| O-16 | 200 | 1.2860 | pass | 0.1140 | -0.0814 |
| O-16 | verdict | 2 of 3 in band | pass | | |

Trigger 4 (q_min / q_H in [0.6, 1.0]): **pass**.

- C-12: q_H 1.912 1/fm; 150 MeV q_min 1.538, q_min/q_H 0.805, in band pass; 200 MeV q_min 1.567, q_min/q_H 0.819, in band pass; verdict cell 150_mev: pass
- O-16: q_H 1.661 1/fm; 150 MeV q_min 1.356, q_min/q_H 0.816, in band pass; 200 MeV q_min 1.382, q_min/q_H 0.832, in band pass; verdict cell 150_mev: pass

Pauli: variant `clementel-villi`, scale 1.0, E_min 80.0 MeV, rho_0 0.17 1/fm^3. Nodes with the factor applied per target: {'Al-27': 43, 'C-12': 43, 'Ca-40': 43, 'N-14': 43, 'O-16': 43, 'P-31': 43, 'Si-28': 43}; black-disk nodes below E_min: {'Al-27': 122, 'C-12': 116, 'Ca-40': 117, 'N-14': 122, 'O-16': 112, 'P-31': 122, 'Si-28': 122}; nodes below the SAID range (black disk kept): {'Al-27': 61, 'C-12': 0, 'Ca-40': 0, 'N-14': 24, 'O-16': 0, 'P-31': 67, 'Si-28': 12}.

Pins:
- inputs JSON sha256 `d4ac274f7b8aaea6df915a06e0735e1b73ea8075c072de22f6badd0235b52747`; Horiuchi JSON sha256 `None`; inputs unpinned: False
- SAID NP_ALFA: `47404925af3cbcf8d98613ef06dc72ca254ef667fc75510bc5d4c427acf72685`, `db95e4dcdf76e12e4fdb480981ef4c353c232d84a23dc2de47b4453f17d13136`, `14399f8deb631857b389f3399b9d7df08b789ebf6fedc7014c6151dc09d960c4`
- SAID NP_SGT: `eec6f944a1a9b85800fbeabd07b716f697a2db0edb5abf5852b7221614652290`, `d3f35fe2b873dc39f8c3c6f9de1efdb916aaec356221e23d0c5debe65a1f120d`, `d4d955c71d4a4a3a198d968e884ff62534b9f4381345ef514be32aa669cd410a`
- SAID PP_ALFA: `b2ca5813878a1c7a9c7e1ce808e8971dbbafaa0e7375240b2d58f6c60e53e252`, `f7eeda20f311db64f38c6e6b09d7b97385bc2ae009cdea0ccceecb9a3d83c7fd`, `4bae553e8a06e82f02b8032ca648e50631b5cee6f8437c316061e577701d8c8a`
- SAID PP_SGT: `508d47be0adcc78aa5546d50d60cec8410d6b85a6c94e8860368ca173a8b1b37`, `b72360a2c3813bd523a57e529ba964b1c30a806c02963b6a941efbd7e7a6bcb8`, `d93b3f4d0c7ae46b941c9ae82292238384624af416d1acccbf4e2146f478954c`
- DSG files: none

## RB

Table id `2a14cdcb93db416f30db0b0fb3bf1607af414590fadd92ea0fe3c06f802ec221`; model `m2prime-pauli-blocked`, variant `M2prime`, Pauli variant `bdc-pn`, Pauli E_min 80.0 MeV, beta `unitarity`, pp level `none`, unfold `rms-matched`, inputs unpinned: False.

| target | E (MeV) | sigma_el / Barashenkov | sigma_R / Barashenkov | sigma_tot in +-10 % band | q_min (1/fm) | q_min / q_H | Pauli factor (peak) |
|---|---|---|---|---|---|---|---|
| Al-27 | 100 | 0.676 | 1.164 | - | 1.136 | - | 0.562 |
| Al-27 | 150 | 1.113 | 1.033 | - | 1.212 | - | 0.708 |
| Al-27 | 200 | 1.475 | 0.966 | - | 1.227 | - | 0.781 |
| C-12 | 100 | 0.717 | 1.162 | fail | 1.505 | 0.787 | 0.658 |
| C-12 | 150 | 1.099 | 1.053 | pass | 1.566 | 0.819 | 0.772 |
| C-12 | 200 | 1.235 | 0.983 | pass | 1.580 | 0.826 | 0.829 |
| Ca-40 | 100 | 0.727 | 1.119 | - | 0.928 | - | 0.601 |
| Ca-40 | 150 | 1.142 | 1.056 | - | 0.986 | - | 0.734 |
| Ca-40 | 200 | 1.423 | 1.006 | - | 1.011 | - | 0.801 |
| N-14 | 100 | 0.750 | 1.069 | - | 1.404 | - | 0.687 |
| N-14 | 150 | 1.228 | 0.991 | - | 1.472 | - | 0.792 |
| N-14 | 200 | 1.231 | 0.950 | - | 1.490 | - | 0.844 |
| O-16 | 100 | 0.699 | 1.123 | fail | 1.307 | 0.787 | 0.723 |
| O-16 | 150 | 1.182 | 1.022 | pass | 1.370 | 0.825 | 0.816 |
| O-16 | 200 | 1.396 | 0.963 | pass | 1.389 | 0.836 | 0.862 |
| P-31 | 100 | 0.697 | 1.129 | - | 1.057 | - | 0.579 |
| P-31 | 150 | 1.170 | 1.021 | - | 1.126 | - | 0.719 |
| P-31 | 200 | 1.491 | 0.966 | - | 1.148 | - | 0.789 |
| Si-28 | 100 | 0.695 | 1.144 | - | 1.101 | - | 0.586 |
| Si-28 | 150 | 1.183 | 1.019 | - | 1.171 | - | 0.724 |
| Si-28 | 200 | 1.512 | 0.963 | - | 1.189 | - | 0.793 |

sigma_tot against the measured values, pass/fail per cell (band +-10 %):

- Al-27: at least 2 of 3: -; all 3: -
- C-12: at least 2 of 3: pass; all 3: fail
- Ca-40: at least 2 of 3: -; all 3: -
- N-14: at least 2 of 3: -; all 3: -
- O-16: at least 2 of 3: pass; all 3: fail
- P-31: at least 2 of 3: -; all 3: -
- Si-28: at least 2 of 3: -; all 3: -

Trigger 2: max 0.000038 against 0.001, passed: True, nodes evaluated 33. Component maxima: b_step 0.000038, q_cutoff 0.000000, q_step 0.000036, radial_step 0.000000, s_grid 0.000002.

Trigger 3 (>= 2 of 3 in [0.7, 1.4], C-12 and O-16): **passes**.

| target | E (MeV) | ratio | in band | margin to nearest bound | excess over upper bound |
|---|---|---|---|---|---|
| C-12 | 100 | 0.7166 | pass | 0.0166 | -0.4881 |
| C-12 | 150 | 1.0992 | pass | 0.3008 | -0.2149 |
| C-12 | 200 | 1.2347 | pass | 0.1653 | -0.1181 |
| C-12 | verdict | 3 of 3 in band | pass | | |
| O-16 | 100 | 0.6987 | fail | -0.0013 | -0.5010 |
| O-16 | 150 | 1.1817 | pass | 0.2183 | -0.1559 |
| O-16 | 200 | 1.3960 | pass | 0.0040 | -0.0028 |
| O-16 | verdict | 2 of 3 in band | pass | | |

Trigger 4 (q_min / q_H in [0.6, 1.0]): **pass**.

- C-12: q_H 1.912 1/fm; 150 MeV q_min 1.566, q_min/q_H 0.819, in band pass; 200 MeV q_min 1.580, q_min/q_H 0.826, in band pass; verdict cell 150_mev: pass
- O-16: q_H 1.661 1/fm; 150 MeV q_min 1.370, q_min/q_H 0.825, in band pass; 200 MeV q_min 1.389, q_min/q_H 0.836, in band pass; verdict cell 150_mev: pass

Pauli: variant `bdc-pn`, scale 1.0, E_min 80.0 MeV, rho_0 0.17 1/fm^3. Nodes with the factor applied per target: {'Al-27': 43, 'C-12': 43, 'Ca-40': 43, 'N-14': 43, 'O-16': 43, 'P-31': 43, 'Si-28': 43}; black-disk nodes below E_min: {'Al-27': 122, 'C-12': 116, 'Ca-40': 117, 'N-14': 122, 'O-16': 112, 'P-31': 122, 'Si-28': 122}; nodes below the SAID range (black disk kept): {'Al-27': 61, 'C-12': 0, 'Ca-40': 0, 'N-14': 24, 'O-16': 0, 'P-31': 67, 'Si-28': 12}.

Pins:
- inputs JSON sha256 `d4ac274f7b8aaea6df915a06e0735e1b73ea8075c072de22f6badd0235b52747`; Horiuchi JSON sha256 `None`; inputs unpinned: False
- SAID NP_ALFA: `47404925af3cbcf8d98613ef06dc72ca254ef667fc75510bc5d4c427acf72685`, `db95e4dcdf76e12e4fdb480981ef4c353c232d84a23dc2de47b4453f17d13136`, `14399f8deb631857b389f3399b9d7df08b789ebf6fedc7014c6151dc09d960c4`
- SAID NP_SGT: `eec6f944a1a9b85800fbeabd07b716f697a2db0edb5abf5852b7221614652290`, `d3f35fe2b873dc39f8c3c6f9de1efdb916aaec356221e23d0c5debe65a1f120d`, `d4d955c71d4a4a3a198d968e884ff62534b9f4381345ef514be32aa669cd410a`
- SAID PP_ALFA: `b2ca5813878a1c7a9c7e1ce808e8971dbbafaa0e7375240b2d58f6c60e53e252`, `f7eeda20f311db64f38c6e6b09d7b97385bc2ae009cdea0ccceecb9a3d83c7fd`, `4bae553e8a06e82f02b8032ca648e50631b5cee6f8437c316061e577701d8c8a`
- SAID PP_SGT: `508d47be0adcc78aa5546d50d60cec8410d6b85a6c94e8860368ca173a8b1b37`, `b72360a2c3813bd523a57e529ba964b1c30a806c02963b6a941efbd7e7a6bcb8`, `d93b3f4d0c7ae46b941c9ae82292238384624af416d1acccbf4e2146f478954c`
- DSG files: none

## Builder-vs-scratch reproduction (M2', six trigger cells; relative differences, builder / scratch - 1)

| target | E (MeV) | sigma_tot | sigma_R | sigma_el | sigma_el / Barashenkov | q_min | max (tot, R, el) | within 2 % |
|---|---|---|---|---|---|---|---|---|
| C-12 | 100 | +0.0009 | +0.0013 | +0.0003 | +0.0008 | +0.0010 | 0.0013 | pass |
| C-12 | 150 | +0.0009 | +0.0013 | +0.0001 | +0.0001 | -0.0012 | 0.0013 | pass |
| C-12 | 200 | +0.0011 | +0.0016 | +0.0002 | -0.0002 | +0.0004 | 0.0016 | pass |
| O-16 | 100 | +0.0010 | +0.0015 | +0.0001 | +0.0009 | +0.0007 | 0.0015 | pass |
| O-16 | 150 | +0.0011 | +0.0014 | -0.0000 | +0.0002 | +0.0016 | 0.0014 | pass |
| O-16 | 200 | +0.0011 | +0.0014 | -0.0002 | -0.0000 | +0.0006 | 0.0014 | pass |

2 % rule: **all six cells within 2 % in sigma_tot, sigma_R, sigma_el**.

Leak check (no absolute cross section in the projection): passed True, 20886 absolute values compared.
