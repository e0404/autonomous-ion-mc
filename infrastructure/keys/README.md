# Trusted attestation keys

`trusted_attestation_keys.json` lists the Ed25519 public keys of the protected
MCP services (`validation-manager`, `host-runner`). Trust is on first use: the
lists are empty on merge. After the servers are restarted, the operator (or the
orchestrator in a reviewed commit) adds each published public key
(`python infrastructure/experiment_v3/attest.py --print-public-key <signer>`) as
`{"key_id": ..., "public_key": ..., "added": "YYYY-MM-DD", "note": ...}`.
A record whose key is not listed is "unverified", never "verified".

## Enforcement

The local-validation gate (`infrastructure/validation/local_validation.py`)
verifies every record against this file and the signed task/SHA identity. Only
`verified` records count as passing; unverified and invalid ones are rejected
with a reason.

**Bootstrap is restricted to the migration task.** The constant
`BOOTSTRAP_TASKS = {"V3-INFRA-002"}` is the only place where unverified records
(no trusted `validation-manager` key published yet) are still counted, reported
as "unverified (no trusted key published yet)". For every other task id,
unverified records are rejected whatever the trust map state, with
"attestation required: no trusted validation-manager key published (operator
runbook) or record not verified". Consequently, **after the V3-INFRA-002 change
merges, every other task's merge gate is blocked until the operator publishes
the keys**; the runbook (SETUP.md) publishes both keys in one sitting and then
re-records validation. Invalid (tampered, wrong-signer, wrong-identity) records
are rejected in all modes.

Validation records that name a host run (`details.host_run_id` or
`details.host_run_ids`) additionally require the signed host-runner record
`<host-runs root>/<task>/<RUN-id>/record.json` (via
`infrastructure/validation/host_run_records.py`): present, verified by a trusted
`host-runner` key (bootstrap-tolerated only for the migration task while the
list is empty), bound to the validation record's task id and SHA and to its own
run id, with `exit_code` 0 and `timed_out` false. A missing, unverified or
invalid run record rejects the validation record. Transport consumers use
`require_verified` from `validation/scripts/transport/attest_verify.py`, which
is strict by default and takes expected task/SHA/run id/manifest values.
