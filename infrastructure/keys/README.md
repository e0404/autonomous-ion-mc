# Trusted attestation keys

`trusted_attestation_keys.json` lists the Ed25519 public keys of the protected
MCP services (`validation-manager`, `host-runner`). Trust is on first use: the
lists are empty on merge. After the servers are restarted, the operator (or the
orchestrator in a reviewed commit) adds each published public key
(`python infrastructure/experiment_v3/attest.py --print-public-key <signer>`) as
`{"key_id": ..., "public_key": ..., "added": "YYYY-MM-DD", "note": ...}`.
A record whose key is not listed is "unverified", never "verified".
