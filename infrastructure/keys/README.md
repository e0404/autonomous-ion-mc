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
verifies every record against this file and the signed task/SHA identity.
Once at least one `validation-manager` key is listed, only `verified` records
count; unverified and invalid ones are rejected with a reason. While the list
is empty (bootstrap, before the operator publishes the key) records are
reported as "unverified (no trusted key published yet)" and tolerated; invalid
(tampered or wrong-signer) records are rejected in both modes. Publishing the
key therefore ends the bootstrap exception: re-record validation afterwards.
Transport consumers call `require_verified` from
`validation/scripts/transport/attest_verify.py` for `host-runner` records.
