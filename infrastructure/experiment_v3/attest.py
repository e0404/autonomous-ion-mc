"""Ed25519 attestations for records produced by the protected MCP services.

Pure standard library (hashlib, os, json): the MCP servers run in an
environment with no third-party crypto package. Ed25519 follows the RFC 8032
section 5.1 reference algorithm. It is NOT constant time (Python big integers);
acceptable here because the signing key never leaves a single-user host and no
remote party can time the signer. Do not reuse this module for multi-tenant use.

A signed record carries an additive "attestation" block. The signed message is
the canonical JSON of {"payload": record without that block, "signer",
"signed_at_utc", "key_id"}, so signer identity and signing time are covered.
payload_sha256 is the hash of the canonical payload alone. Trust is decided by
looking up the public key in a per-signer trusted-key map
({signer: {key_id: public_key_hex}}, as loaded by load_trusted).

The verify half is duplicated in validation/scripts/transport/attest_verify.py
for consumers that cannot import infrastructure; keep them in sync.
"""

import argparse
import hashlib
import json
import os
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

SCHEMA = "ionmc-attestation-1"
ALG = "ed25519"
SIGNERS = ("validation-manager", "host-runner")

# ---------------------------------------------------------------- Ed25519 --

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
    return (p[0] * q[2] - q[0] * p[2]) % _P == 0 and (
        p[1] * q[2] - q[1] * p[2]
    ) % _P == 0


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


def _expand(seed):
    if len(seed) != 32:
        raise ValueError("Ed25519 private key seed must be 32 bytes")
    h = _sha512(seed)
    a = int.from_bytes(h[:32], "little")
    a &= (1 << 254) - 8
    a |= 1 << 254
    return a, h[32:]


def public_key_from_seed(seed):
    a, _ = _expand(seed)
    return _compress(_mul(a, _G))


def sign(message, private_key_seed):
    """Return the 64-byte Ed25519 signature of message."""
    a, prefix = _expand(private_key_seed)
    pub = _compress(_mul(a, _G))
    r = int.from_bytes(_sha512(prefix + message), "little") % _L
    big_r = _compress(_mul(r, _G))
    k = int.from_bytes(_sha512(big_r + pub + message), "little") % _L
    s = (r + k * a) % _L
    return big_r + int.to_bytes(s, 32, "little")


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
        k = (
            int.from_bytes(_sha512(signature[:32] + public_key + message), "little")
            % _L
        )
        return _equal(_mul(s, _G), _add(point_r, _mul(k, point_a)))
    except Exception:
        return False


# ------------------------------------------------------------- keys/records --


class AttestationError(Exception):
    pass


@dataclass(frozen=True)
class KeyPair:
    seed: bytes
    public_key: bytes

    @property
    def key_id(self):
        return key_id_of(self.public_key)


@dataclass(frozen=True)
class VerifiedAttestation:
    signer: str
    key_id: str
    public_key: str
    signed_at_utc: str
    payload_sha256: str


def key_id_of(public_key):
    return hashlib.sha256(public_key).hexdigest()[:16]


def keypair_from_seed(seed):
    return KeyPair(seed=seed, public_key=public_key_from_seed(seed))


def generate_keypair():
    return keypair_from_seed(os.urandom(32))


def canonical_payload(record):
    """Canonical JSON bytes of the record without its attestation field."""
    payload = {k: v for k, v in record.items() if k != "attestation"}
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def signed_message(record, signer, signed_at_utc, key_id):
    """Bytes actually signed: canonical JSON binding payload, signer, time, key."""
    body = {
        "payload": {k: v for k, v in record.items() if k != "attestation"},
        "signer": signer,
        "signed_at_utc": signed_at_utc,
        "key_id": key_id,
    }
    return json.dumps(
        body, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def load_trusted(path=None):
    """Return {signer: {key_id: public_key_hex}} from the trust file."""
    default = (
        Path(__file__).resolve().parents[1] / "keys" / "trusted_attestation_keys.json"
    )
    data = json.loads(Path(path or default).read_text(encoding="utf-8"))
    return {
        s: {e["key_id"]: e["public_key"] for e in entries}
        for s, entries in data.items()
    }


def sign_record(record, key, signer):
    if signer not in SIGNERS:
        raise AttestationError(f"unknown signer: {signer!r}")
    payload = canonical_payload(record)
    signed = {k: v for k, v in record.items() if k != "attestation"}
    signed_at = datetime.now(UTC).isoformat()
    message = signed_message(record, signer, signed_at, key.key_id)
    signed["attestation"] = {
        "schema": SCHEMA,
        "alg": ALG,
        "signer": signer,
        "key_id": key.key_id,
        "public_key": key.public_key.hex(),
        "signed_at_utc": signed_at,
        "payload_sha256": hashlib.sha256(payload).hexdigest(),
        "signature": sign(message, key.seed).hex(),
    }
    return signed


def verify_record(record, trusted, expected_signer=None):
    """Verify ``record``; ``trusted`` is ``{signer: {key_id: public_key_hex}}``.

    Raises ``AttestationError`` on any failure.
    """
    block = record.get("attestation")
    if not isinstance(block, dict):
        raise AttestationError("record has no attestation")
    if block.get("schema") != SCHEMA or block.get("alg") != ALG:
        raise AttestationError("unsupported attestation schema or algorithm")
    try:
        public_key = bytes.fromhex(block["public_key"])
        signature = bytes.fromhex(block["signature"])
        key_id, signer = block["key_id"], block["signer"]
        claimed_hash = block["payload_sha256"]
        signed_at = block["signed_at_utc"]
    except (KeyError, ValueError, TypeError) as exc:
        raise AttestationError(f"malformed attestation: {exc}") from exc
    if expected_signer is not None and signer != expected_signer:
        raise AttestationError(f"signer {signer!r} != expected {expected_signer!r}")
    if key_id != key_id_of(public_key):
        raise AttestationError("key_id does not match public key")
    payload = canonical_payload(record)
    if hashlib.sha256(payload).hexdigest() != claimed_hash:
        raise AttestationError("payload hash mismatch")
    if not verify(
        signed_message(record, signer, signed_at, key_id), signature, public_key
    ):
        raise AttestationError("invalid signature")
    if trusted.get(signer, {}).get(key_id) != block["public_key"]:
        raise AttestationError("signing key is not trusted")
    return VerifiedAttestation(
        signer, key_id, block["public_key"], str(signed_at), claimed_hash
    )


def key_dir():
    """Protected per-experiment attestation directory (monkeypatched in tests)."""
    from infrastructure.experiment_v3.paths import SHARE_ROOT

    return SHARE_ROOT / "attestation"


def _write_public(directory, signer, key):
    pub = directory / f"{signer}.pub"
    if not pub.exists():
        pub.write_text(key.public_key.hex() + "\n", encoding="utf-8")
        os.chmod(pub, 0o644)


def load_or_create_key(signer, directory=None):
    """Load (or create with mode 0600) the signer's seed; ensure the .pub file."""
    if signer not in SIGNERS:
        raise AttestationError(f"unknown signer: {signer!r}")
    directory = Path(directory) if directory is not None else key_dir()
    directory.mkdir(parents=True, exist_ok=True)
    os.chmod(directory, 0o700)
    path = directory / f"{signer}.key"
    if path.exists():
        if path.stat().st_mode & 0o077:
            raise AttestationError(f"{path} must not be accessible by group/others")
        key = keypair_from_seed(bytes.fromhex(path.read_text(encoding="utf-8").strip()))
    else:
        key = generate_keypair()
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(key.seed.hex() + "\n")
    _write_public(directory, signer, key)
    return key


def sign_with_service_key(record, signer):
    return sign_record(record, load_or_create_key(signer), signer)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--print-public-key", choices=SIGNERS, required=True)
    args = parser.parse_args(argv)
    key = load_or_create_key(args.print_public_key)
    print(
        json.dumps(
            {
                "signer": args.print_public_key,
                "key_id": key.key_id,
                "public_key": key.public_key.hex(),
            }
        )
    )
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    raise SystemExit(main())
