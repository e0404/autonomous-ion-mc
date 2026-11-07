"""Stdlib-only attestation verifier for repository-side consumers.

Verify half of infrastructure/experiment_v3/attest.py, duplicated so that
validation scripts (run in the host sandbox without infrastructure imports) can
check "attestation" blocks on VAL-/RUN- records. Keep in sync with attest.py
(tests assert both agree). Trusted keys: infrastructure/keys/trusted_attestation_keys.json.

Status values: "verified" (valid signature by a listed key of the expected
signer), "unverified" (no attestation, or signing key not listed), "invalid"
(attestation present but malformed, tampered, or wrongly signed).
"""

import hashlib
import json
from pathlib import Path

SCHEMA = "ionmc-attestation-1"
ALG = "ed25519"
TRUSTED_KEYS_PATH = (
    Path(__file__).resolve().parents[3] / "infrastructure" / "keys" / "trusted_attestation_keys.json"
)

_P = 2**255 - 19
_L = 2**252 + 27742317777372353535851937790883648493
_D = -121665 * pow(121666, _P - 2, _P) % _P
_I = pow(2, (_P - 1) // 4, _P)


def _sha512(data):
    return hashlib.sha512(data).digest()


def _inv(x):
    return pow(x, _P - 2, _P)


def _recover_x(y, sign):
    if y >= _P:
        return None
    x2 = (y * y - 1) * _inv(_D * y * y + 1) % _P
    if x2 == 0:
        return None if sign else 0
    x = pow(x2, (_P + 3) // 8, _P)
    if (x * x - x2) % _P != 0:
        x = x * _I % _P
    if (x * x - x2) % _P != 0:
        return None
    if (x & 1) != sign:
        x = _P - x
    return x


_GY = 4 * _inv(5) % _P
_GX = _recover_x(_GY, 0)
_G = (_GX, _GY, 1, _GX * _GY % _P)  # extended coordinates (X, Y, Z, T)
_IDENTITY = (0, 1, 1, 0)


def _add(p, q):
    a = (p[1] - p[0]) * (q[1] - q[0]) % _P
    b = (p[1] + p[0]) * (q[1] + q[0]) % _P
    c = 2 * p[3] * q[3] * _D % _P
    d = 2 * p[2] * q[2] % _P
    e, f, g, h = b - a, d - c, d + c, b + a
    return (e * f % _P, g * h % _P, f * g % _P, e * h % _P)


def _mul(s, p):
    q = _IDENTITY
    while s > 0:
        if s & 1:
            q = _add(q, p)
        p = _add(p, p)
        s >>= 1
    return q


def _equal(p, q):
    return (p[0] * q[2] - q[0] * p[2]) % _P == 0 and (p[1] * q[2] - q[1] * p[2]) % _P == 0


def _compress(p):
    zi = _inv(p[2])
    x, y = p[0] * zi % _P, p[1] * zi % _P
    return int.to_bytes(y | ((x & 1) << 255), 32, "little")


def _decompress(s):
    if len(s) != 32:
        return None
    y = int.from_bytes(s, "little")
    sign = y >> 255
    y &= (1 << 255) - 1
    x = _recover_x(y, sign)
    if x is None:
        return None
    return (x, y, 1, x * y % _P)


def verify(message, signature, public_key):
    """Return True iff signature is a valid Ed25519 signature (never raises)."""
    try:
        if len(public_key) != 32 or len(signature) != 64:
            return False
        point_a = _decompress(public_key)
        point_r = _decompress(signature[:32])
        if point_a is None or point_r is None:
            return False
        s = int.from_bytes(signature[32:], "little")
        if s >= _L:
            return False
        k = int.from_bytes(_sha512(signature[:32] + public_key + message), "little") % _L
        return _equal(_mul(s, _G), _add(point_r, _mul(k, point_a)))
    except Exception:
        return False


def _signed_message(record, signer, signed_at_utc, key_id):
    body = {
        "payload": {k: v for k, v in record.items() if k != "attestation"},
        "signer": signer,
        "signed_at_utc": signed_at_utc,
        "key_id": key_id,
    }
    return json.dumps(
        body, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def _canonical_payload(record):
    payload = {k: v for k, v in record.items() if k != "attestation"}
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def load_trusted(path=None):
    """Return {signer: {key_id: public_key_hex}} from the trust file."""
    data = json.loads(Path(path or TRUSTED_KEYS_PATH).read_text(encoding="utf-8"))
    return {
        signer: {e["key_id"]: e["public_key"] for e in entries}
        for signer, entries in data.items()
    }


def verify_record(record, expected_signer, trusted=None):
    """Return {"status": verified|unverified|invalid, "reason": str, "signer", "key_id"}."""
    block = record.get("attestation")
    if not isinstance(block, dict):
        return {"status": "unverified", "reason": "no attestation"}
    try:
        if block.get("schema") != SCHEMA or block.get("alg") != ALG:
            raise ValueError("unsupported schema or algorithm")
        public_key = bytes.fromhex(block["public_key"])
        signature = bytes.fromhex(block["signature"])
        key_id, signer = block["key_id"], block["signer"]
        signed_at = block["signed_at_utc"]
        if signer != expected_signer:
            raise ValueError(f"signer {signer!r} != expected {expected_signer!r}")
        if key_id != hashlib.sha256(public_key).hexdigest()[:16]:
            raise ValueError("key_id does not match public key")
        payload = _canonical_payload(record)
        if hashlib.sha256(payload).hexdigest() != block["payload_sha256"]:
            raise ValueError("payload hash mismatch")
        if not verify(_signed_message(record, signer, signed_at, key_id), signature, public_key):
            raise ValueError("invalid signature")
    except (KeyError, ValueError, TypeError) as exc:
        return {"status": "invalid", "reason": str(exc)}
    keys = (load_trusted() if trusted is None else trusted).get(expected_signer, {})
    if keys.get(key_id) != block["public_key"]:
        return {"status": "unverified", "reason": "signing key is not trusted",
                "signer": signer, "key_id": key_id}
    return {"status": "verified", "reason": "ok", "signer": signer, "key_id": key_id}


class AttestationRequired(Exception):
    """A record was rejected by require_verified."""


def require_verified(
    record,
    expected_signer,
    trust=None,
    *,
    expected_task_id=None,
    expected_git_sha=None,
    expected_run_id=None,
    expected_output_manifest_sha256=None,
    bootstrap=False,
):
    """Raise AttestationRequired unless ``record`` may be consumed as evidence.

    Strict by default: only a ``verified`` record passes; ``invalid`` and
    ``unverified`` raise. ``bootstrap=True`` (explicit; for the migration task
    only) additionally tolerates ``unverified`` while the trust map lists no
    key for ``expected_signer``; invalid always raises.

    Identity binding (replay protection): each ``expected_*`` that is not None
    must equal the signed record's field, else AttestationRequired. The git
    SHA is read from ``git_sha`` (VAL- records) or ``sha`` (RUN- records).
    Returns the verify_record result dict.

    Production consumer: infrastructure/validation/local_validation.py (via
    infrastructure/validation/host_run_records.py) for RUN- records.
    """
    trusted = load_trusted() if trust is None else trust
    result = verify_record(record, expected_signer, trusted)
    if result["status"] == "invalid":
        raise AttestationRequired(f"invalid attestation: {result['reason']}")
    git_sha = record.get("git_sha", record.get("sha"))
    expected = (
        ("task_id", expected_task_id, record.get("task_id")),
        ("git_sha", expected_git_sha, git_sha),
        ("run_id", expected_run_id, record.get("run_id")),
        (
            "output_manifest_sha256",
            expected_output_manifest_sha256,
            record.get("output_manifest_sha256"),
        ),
    )
    for name, want, got in expected:
        if want is not None and got != want:
            raise AttestationRequired(
                f"signed identity mismatch: {name} is {got!r}, expected {want!r}"
            )
    if result["status"] != "verified":
        if bootstrap and not trusted.get(expected_signer):
            return result
        raise AttestationRequired(
            f"{result['status']} record rejected (attestation required for "
            f"{expected_signer!r}): {result['reason']}"
        )
    return result
