# V3-006 step S8 evidence: the `out6` workflow runs (plan `validation/plans/v3-006-acceptance.md`)

Driver: `validation/scripts/results/out6_run.py` (cases in `out6_cases.py`, templates in
`validation/scripts/results/out6/`). Every file here is a small JSON record of identities, digests, exit codes
and verdicts; no bundle, array or table value is committed (REG-4). The summaries and verdict files written by the
driver carry only paths relative to their run directory and no value of the path set P (decision 0045 §2). The one
exception is the signed host-runner record `hr/host-records/RUN-*.json`: it is the runner's own document, names
the absolute run, snapshot and output paths of the host (`/home/...`, `/workspace/...`) and is committed unaltered
because its signature covers those bytes.

## `lv/` — local-validation tier (this sandbox, `IONMC_LV_TIER=1`, one process, thread limits 1)

Produced at the clean commit 154fc2da (every summary's `code` block: `source = git`, `git_dirty = false`,
`git_commit = 154fc2da…`, `package_tree_oid = 4970c9a4…` = `git rev-parse 154fc2da:src/ionmc`) by the logged run
`20261125T100614Z-out6-lv-evidence-154fc2da` (`telemetry/test-logs`). Earlier sets of the same rows at 48b782bf,
cebd9109 and 154fc2da (all PASS) were superseded by the reader and driver corrections of the final review rounds.

| file | row(s) | verdict |
|---|---|---|
| `r0k0.summary.json` … `r7k0.summary.json`, `r6k1.summary.json` | OUT-2 (LV part: `ionmc run`, clean-process `ionmc verify` and `ionmc report --check` in `python -I` children, exit codes 0), CLI-4 (LV part: r5 nuclear and r7 `icru90-water-v1` from the populated cache, no network in this sandbox), CLI-6 inputs | RUN PASS for r0, r1, r3, r4, r5, r6 (k = 0 and k = 1), r7 |
| `compare.json` | CLI-5 (LV part, r0 and r1; r2 and the three-way comparison are HR) | CLI-5-LV PASS |
| `identities-carbon.json`, `identities-nuclear.json`, `identities-icru90.json` | CLI-6 (r4 carbon, r5 nuclear table id and npz digest, r7 `icru90-water-v1` dataset identity) | PASS, PASS, PASS |
| `reproduce-r6k0.json`, `reproduce-r6k1.json` | OUT-3 (`configio.from_dict(requested)` + the same seed on python/float64 and warp-cpu/float64: every array content digest and the effective block equal) | PASS, PASS |

Seeds: `out6` family, seed = 20700004 + 100·r_index + k (k = 0; r6 k = 1 for the warp-cpu original of OUT-3).

Preserved failure: the first LV attempt at 1f371352 (logged run `20261125T084343Z-out6-lv-evidence`) passed every
RUN, CLI-5-LV and CLI-6 row but reported **OUT-3 FAIL** for r6k0 and r6k1 (child exit code 1): the driver started
the reproduction child with its working directory set to the output directory and passed the original bundle as a
repository-relative path, which did not resolve there. The driver was corrected at 48b782bf (the child receives the
resolved path; regression test with a relative path) and the whole set was regenerated at that commit. No tolerance
and no row definition changed.

## `hr/` — host-runner tier (controlled GPU host, read-only snapshot, no network)

Filled by the controlled host run RUN-20261125T100705Z-1cae10a5 at the clean commit 154fc2da (signed record in
`hr/host-records/`, verified with `validation/scripts/transport/attest_verify.py`): RUN-r2 (warp-cuda, `cuda:0`,
float32), RUN-r0 and RUN-r1 on the host, CLI-5 (three-way r0/r1/r2) PASS, PROV-7 PASS, OUT-2-HR PASS, CLI-4-HR PASS,
REG-1-HR PASS and CAP-1b-HR PASS (90 contract tests with CUDA present, 0 skipped). See `hr/README.md`.

Not a clinical claim.
