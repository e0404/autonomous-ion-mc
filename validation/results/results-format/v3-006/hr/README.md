# `hr/` — host-runner tier of the V3-006 `out6` evidence

Controlled host run **RUN-20261125T100705Z-1cae10a5** (`host-records/RUN-20261125T100705Z-1cae10a5.json`, the
signed record of the host runner; signer `host-runner`, key id 4d02f521f3fb9344; verify with
`validation/scripts/transport/attest_verify.py`) at the clean commit 154fc2da (`dirty_before = false`,
`dirty_after = false`, `committed_state_unchanged = true`, exit code 0, 33 s). The snapshot is mounted read-only at
`/workspace`; the run has no network; `HOME=/tmp/home`; every child single-process with thread limits 1;
`IONMC_REQUIRE_DATA=1` and `IONMC_LV_TIER=1` (a missing cache object or a skipped test is a failure, never a skip).
Driver entry: `validation/scripts/results/out6_run.py host --out validation/generated/out6-hr --cache-dir
/workspace/.ionmc-cache/ionmc-data --sha 154fc2da… --tree-oid 4970c9a4…` (the full argv is in the record).

| file | row(s) | verdict |
|---|---|---|
| `r2k0.summary.json` | OUT-2 HR (r2 on warp-cuda, device `cuda:0`, float32: `ionmc run`, clean-process `ionmc verify` and `ionmc report --check` exit 0), CLI-4 HR (network denied by the host sandbox) | RUN-r2 PASS |
| `r0k0.summary.json`, `r1k0.summary.json` | the r0 (python/float64) and r1 (warp-cpu/float32) bundles of the same box on the host for the three-way comparison | RUN-r0 PASS, RUN-r1 PASS |
| `compare.json` | CLI-5 (r0/r1/r2 bind to their verified bundles; they differ only in the backend/device/precision/production pointers and their effects: RNG `u01` mapping, run and chunk reconciliation rows, config digest, array contents and derived-file hashes, accounting and statistics values, timings, output path; the formula seeds) | CLI-5 PASS |
| `prov7.json` | PROV-7 (the r2 bundle written under the host runner: `code.source = git`, `git_dirty = false`, `git_commit` = 154fc2da, `package_tree_oid` = `git rev-parse 154fc2da:src/ionmc` = 4970c9a4…) | PASS |
| `host-summary.json` | OUT-2-HR (r2 exit codes AND `test_out2_clean_process_round_trip`: 1 passed, 0 skipped), CLI-4-HR, REG-1-HR (`test_default_trace_is_unchanged_from_the_stored_baseline`: 1 passed; the four A16 digest tests: 4 passed; 0 skipped), CAP-1b-HR (`tests/ionmc/test_v3_006_contract.py` with CUDA present: 90 passed, 0 skipped, so the CAP-1a CUDA cells and the cache-backed CAP-1b matrix ran on the host) | PASS, PASS, PASS, PASS |

The pytest summary lines of the host children (from the run's published logs):
`out2 → 1 passed`, `contract → 90 passed`, `trace → 1 passed`, `a16 → 4 passed`; no skips, failures or errors.
Every summary records identities, digests, exit codes and verdicts only; the bundles stay in the runner-published
output directory (ignored by Git) and are not committed (REG-4).

Not a clinical claim.

The same host run also re-derived the V3-013 predictions file (`validation/scripts/stopping/steps_v13.py predict`,
`vc-predictions.json` sha256 61df725f…): it differs from the previously committed file only in the `src/ionmc` tree
pin (c5b33495… → 4970c9a4…, the tree of the V3-006 head 154fc2da); every value, identity and table hash is
unchanged (see the dated note of 2026-11-25 in `validation/plans/v3-013-acceptance.md`).
