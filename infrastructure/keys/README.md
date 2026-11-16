# Trusted attestation keys

`trusted_attestation_keys.json` lists the Ed25519 public keys of the protected
MCP services (`validation-manager`, `host-runner`, and optionally
`reference-runner`, V3-INFRA-004, which signs `REF-...` reference-run
records). Trust is on first use: the lists are empty on merge. The gate
requires the first two signers; `reference-runner` may be absent or `[]` until
its key-publication commit (records it signed are then `unverified`), and once
listed it obeys the same rules (exactly one key in a key-publication branch,
evidenced by a verified successful reference run referenced from a passing
validation record in `details.reference_runs`; see `experiment/v3/SETUP.md`,
"Signed reference-run records"). A signer that the integration copy lists
(any signer, optional ones included) may never be removed or emptied by a
branch; it can only be rotated. `codex-worker` (V3-INFRA-006, signs
`REVIEW-...` Codex review records) is optional in the same way; see "Publishing
the codex-worker key" below. After the servers are restarted, the operator (or the
orchestrator in a reviewed commit) adds each published public key
(`python infrastructure/experiment_v3/attest.py --print-public-key <signer>`) as
`{"key_id": ..., "public_key": ..., "added": "YYYY-MM-DD", "note": ...}`.
A record whose key is not listed is "unverified", never "verified".

**Publishing the first keys (key-publication mode).** The merge gate reads
this file from the integration checkout, so the first key-publication commit
cannot be verified by the integration copy. The gate therefore evaluates a task
branch whose only change against the integration base is this file, listing
exactly one well-formed (non-small-order) key for both signers, against the
branch's own copy (`trust_source: "key-publication"`), with full enforcement,
no bootstrap tolerance, and signed evidence required for every published key
(a verified passing validation record that references a verified successful
host run of that task and SHA). The same rule governs every later change of
this file (rotation = replace, evidenced by the new keys): a branch that changes
it together with anything else, or lists more than one key per signer, fails
the gate (`trust_problem`) even when the integration copy already trusts a key.
The CI test `test_committed_trust_file_is_empty_or_strictly_published` accepts
both committed states (empty, or exactly one valid key per signer). See
`experiment/v3/SETUP.md` (Signed attestations) for the exact rules.

**Publishing the codex-worker key (V3-INFRA-006).** The codex-worker MCP
service creates `<share root>/attestation/codex-worker.key` on first use, when
the first review finishes after the operator restarted the service with the
V3-INFRA-006 change, and signs every finished review as
`<state>/reviews/<REVIEW-id>/record.json`. Until a `"codex-worker"` entry is
published here (absent or `[]`), those records are "unverified (key not
trusted)"; the local-validation gate does not require the signer and stays
green. Publication procedure: (a) after the restart, let one review finish (a
discovery review; its record is unverified and not evidence) and read the
`attestation.key_id` and `attestation.public_key` of its signed record
(`inspect_codex_review` returns the record with its attestation block); check
key_id = sha256(public key)[:16] and cross-check with `python
infrastructure/experiment_v3/attest.py --print-public-key codex-worker` run by
the operator outside the agent sandbox. (b) Make a pure key-publication branch
whose only change adds exactly that entry under `"codex-worker"` (`added`, and a
`note` naming the discovery REVIEW record). (c) Record validation at that
commit's SHA as for every key publication, and obtain a Codex review of the
same SHA: the gate then requires, next to the evidence of the other published
keys, a verified, passed review record of this task and SHA signed by the
published key (`key_publication_evidence` gains a `codex-worker` entry). (d)
Push, CI, merge. Reviews finished before the restart have no `record.json` and
are returned with `attestation: null`. Verify a record with `python -m
infrastructure.validation.review_records <record.json>` (exit 0 only when
verified).

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

Validation records that name host runs do so with
`details.host_runs = [{"run_id": "RUN-...", "output_manifest_sha256": "<64 hex
or null>"}]`; `record_local_validation` fills a missing digest from the signed
run record itself (a missing run record marks the validation record
`host_run_manifest_lookup: "missing"`, which the gate rejects). Each entry
additionally requires the signed host-runner record
`<host-runs root>/<task>/<RUN-id>/record.json` (via
`infrastructure/validation/host_run_records.py`): present, verified by a trusted
`host-runner` key (bootstrap-tolerated only for the migration task while the
list is empty), bound to the validation record's task id and SHA and to its own
run id, carrying a signed `output_manifest_sha256` equal to the claimed one
(null must match null; a record without the field is invalid), and succeeded:
`succeeded` true, `exit_code` 0, `timed_out` false,
`committed_state_unchanged` true and `dirty_after` false (when present), no
`output_publication_error`. The legacy `details.host_run_id` / `host_run_ids`
forms carry no manifest claim and are rejected outside the bootstrap task. A
missing, unverified or invalid run record rejects the validation record. Transport consumers use
`require_verified` from `validation/scripts/transport/attest_verify.py`, which
is strict by default and takes expected task/SHA/run id/manifest values (manifest default
`OMIT` = no expectation, explicit `None` = signed null); its `bootstrap=True` is accepted only for
an `expected_task_id` in the mirrored `BOOTSTRAP_TASKS` while the signer has no trusted key.
