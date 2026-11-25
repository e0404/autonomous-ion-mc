# V3-006 step S8 evidence: the `out6` workflow runs (plan `validation/plans/v3-006-acceptance.md`)

Driver: `validation/scripts/results/out6_run.py` (cases in `out6_cases.py`, templates in
`validation/scripts/results/out6/`). Every file here is a small JSON record of identities, digests, exit codes
and verdicts; no bundle, array or table value is committed (REG-4). Paths inside the records are relative to the
run directory; no value of the path set P appears (decision 0045 §2).

## `lv/` — local-validation tier (this sandbox, `IONMC_LV_TIER=1`, one process, thread limits 1)

Produced at the clean commit 48b782bf (every summary's `code` block: `source = git`, `git_dirty = false`,
`git_commit = 48b782bf…`, `package_tree_oid = 2be42e2a…` = `git rev-parse 48b782bf:src/ionmc`) by the logged run
`20261125T084654Z-out6-lv-evidence-48b782bf` (`telemetry/test-logs`).

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

Filled by the controlled host run of `out6_run.py host` (r2 on warp-cuda/float32 with OUT-2-HR, the three-way
CLI-5 comparison r0/r1/r2 on the host, PROV-7 of the r2 bundle, CAP-1a/CAP-1b with CUDA present, the trace baseline
and the A16 digest for REG-1); see `hr/README.md` once present. Rows not covered by `lv/`: OUT-2 HR, CLI-4 HR
(network denied by the host sandbox), CLI-5 (three-way), PROV-7, REG-1 host, CAP-1b HR.

Not a clinical claim.
