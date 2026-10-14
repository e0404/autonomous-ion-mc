Exploratory, report-only (Amendment 19); not evidence; changes no row.

# D-E0 / H5: EM-only peak/plateau, ionmc against TOPAS (archived V11 EM-only inputs)

Command (from the worktree root, `.venv` python, single thread): `python -m validation.scripts.elastic_attribution.de0_h5`

Input files (sha256):

- `validation/generated/RUN-20261012T142957Z-a8e41216/transport/lv5c-dd976c5-16/16-v11-compare.txt`  1f81c774d5209761fce23f757e6f6c2d4f731c533be9b31d0e017a253175f6e2
- `.ionmc-cache/transport-samples/dd976c5/13/samples/v11-150-emonly.json`  cb46ec341d0334740ea0e585f5ea9c0f084d60a883c4097c80a5b070837c0b6f
- `.ionmc-cache/transport-samples/dd976c5/15/samples/v11-200-emonly.json`  1bb5e78e387cdd8528895d6c22b196fe1ba3f792db34d72a4bbcbbd81cc7866b
- `.ionmc-cache/reference-runs/REF-278f82f0abda74861681-8a56e570/work/idd_dose.csv`  b6e7a374f62dae58dfe6c6b75e0841e6376704d0eaf2b6cb54f07c5931b9178e
- `.ionmc-cache/reference-runs/REF-7072862d365da2cf275b-475cdb44/work/idd_dose.csv`  76acde29effcb10f067868a393257c1d89daa3448d83e7686205181ca40301cd
- `.ionmc-cache/reference-runs/REF-b9387263d9a5d9db0593-7c20f8b2/work/idd_dose.csv`  81a9d8e18ae7384728da33403366f5231227b14704195f8813c532dcdcce1a54
- `.ionmc-cache/reference-runs/REF-095c91aab5a7a6d92ebe-7e80de3c/work/idd_dose.csv`  cd23c00ed5eee4c39882830da83e08dcefd4701d220349ec323c79d7c70ed961
- `.ionmc-cache/reference-runs/REF-c10bdf1824e76827a539-b1ded06c/work/idd_dose.csv`  37de1c34504970153951921b931383fcf488e7b9eaea855ba74b12c08b13fedb
- `.ionmc-cache/reference-runs/REF-f3ed032bc11595bea382-3d206d0e/work/idd_dose.csv`  973cf76d0ec5477c90e37585501c0088b35ac7048ea6455340bf0b989115c0d5

Plateau window 20-60 mm and all metric definitions are those of the V11 compare step (module docstring). Peak/plateau and plateau: ionmc batch-jackknife, TOPAS mean of the 3 frozen seeds. 90 % intervals by the compare step's Welch df. Only ratios and metrics are recorded.
<!-- body below is rendered from the Git projection of the JSON record -->

| E | peak/plateau ionmc / TOPAS (90 % CI) | plateau ratio | R80 ionmc - TOPAS (mm) | within 2 % (point / CI) | direction | verdict |
|---|---|---|---|---|---|---|
| 150 | 1.0086 [1.0048, 1.0124] | 1.0030 [1.0027, 1.0033] | -0.334 [-0.361, -0.308] | True / True | amplifies Delta_pk (ionmc EM-only peak/plateau above TOPAS) | H5 excluded (EM-only within +-2 %) |
| 200 | 1.0111 [1.0041, 1.0181] | 0.9995 [0.9992, 0.9997] | -0.340 [-0.390, -0.290] | True / True | amplifies Delta_pk (ionmc EM-only peak/plateau above TOPAS) | H5 excluded (EM-only within +-2 %) |

| E | EM-only contribution to F_ionmc/F_TOPAS - 1 if uncorrected (= 1/ratio - 1) | V11 F_ionmc/F_TOPAS - 1 (Delta_pk, archived) | ionmc peak/plateau (EM-only) | TOPAS peak/plateau (EM-only) |
|---|---|---|---|---|
| 150 | -0.0085 | -0.0390 | 6.1902 +- 0.0117 | 6.1374 +- 0.0064 |
| 200 | -0.0110 | -0.0488 | 6.4611 +- 0.0120 | 6.3900 +- 0.0171 |

Rule (plan H5): excluded if the EM-only peak/plateau agree within +-2%; H5 counts only if the difference exceeds 2 % and amplifies Delta_pk. F is a ratio, so a common EM difference cancels to first order.
