"""Public reference-data acquisition with content-addressed provenance and use roles."""

from __future__ import annotations

import ipaddress
import socket
import urllib.parse
import urllib.request
import uuid
from pathlib import Path

from infrastructure.experiment_v3.common import STATE, event, file_hash, now, write_json

ROLES = {"construction", "calibration", "evaluation", "exploratory"}
MAX_BYTES = 256 * 1024**2


def public_url(url):
    parsed = urllib.parse.urlsplit(url)
    if (
        parsed.scheme != "https"
        or not parsed.hostname
        or parsed.username
        or parsed.password
        or parsed.fragment
    ):
        raise ValueError(
            "Only public HTTPS URLs without credentials/fragments are supported"
        )
    if parsed.port not in (None, 443):
        raise ValueError("Only HTTPS port 443 is supported")
    if parsed.hostname.lower() == "localhost":
        raise ValueError("Local network destinations are forbidden")
    try:
        addresses = socket.getaddrinfo(parsed.hostname, 443, type=socket.SOCK_STREAM)
    except OSError as exc:
        raise ValueError("Public hostname cannot be resolved") from exc
    if not addresses or any(
        not ipaddress.ip_address(item[4][0]).is_global for item in addresses
    ):
        raise ValueError("Private/local network destinations are forbidden")
    return url


class PublicRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        public_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def assign_role(content_hash, role, reason, *, sha, state=None):
    if role not in ROLES or not reason.strip():
        raise ValueError("A valid use role and scientific rationale are required")
    if len(content_hash) != 64 or any(
        x not in "0123456789abcdef" for x in content_hash
    ):
        raise ValueError("Invalid content SHA256")
    state = Path(state or STATE)
    obj = state / "data" / "objects" / content_hash
    if not obj.is_file() or file_hash(obj) != content_hash:
        raise ValueError("Missing or corrupted cached object")
    record = {
        "content_sha256": content_hash,
        "role": role,
        "reason": reason,
        "code_sha": sha,
        "at": now(),
    }
    use_id = uuid.uuid4().hex
    write_json(state / "data" / "uses" / f"{use_id}.json", record)
    event("data_role_assigned", sha=sha, details=record, state=state)
    return record | {"use_id": use_id}


def acquire(
    url,
    license_basis,
    citation,
    role,
    rationale,
    *,
    sha,
    expected_sha256=None,
    state=None,
    post_body=None,
):
    if (
        not license_basis.strip()
        or not citation.strip()
        or not rationale.strip()
        or role not in ROLES
    ):
        raise ValueError(
            "License/access basis, citation, use role and rationale required"
        )
    public_url(url)
    state = Path(state or STATE)
    temp = state / "data" / "incoming" / uuid.uuid4().hex
    temp.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        data=post_body.encode() if post_body is not None else None,
        headers={
            "User-Agent": "IonMC-reference-research/2",
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    # Use configured public egress proxies only when they contain no credentials.
    # These are trusted service environment settings, never agent request fields.
    proxies = {}
    for scheme, proxy in urllib.request.getproxies().items():
        parsed = urllib.parse.urlsplit(proxy)
        if (
            scheme in ("http", "https")
            and parsed.scheme in ("http", "https")
            and parsed.hostname
            and not parsed.username
            and not parsed.password
        ):
            proxies[scheme] = proxy
    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler(proxies), PublicRedirect()
    )
    try:
        with opener.open(request, timeout=30) as response, temp.open("wb") as out:
            final_url = public_url(response.url)
            total = 0
            while block := response.read(65536):
                total += len(block)
                if total > MAX_BYTES:
                    raise ValueError("Download exceeds 256 MiB limit")
                out.write(block)
        content_hash = file_hash(temp)
        if expected_sha256 is not None and content_hash != expected_sha256:
            raise ValueError("Upstream checksum mismatch")
        dest = state / "data" / "objects" / content_hash
        dest.parent.mkdir(parents=True, exist_ok=True)
        temp.replace(dest)
    finally:
        temp.unlink(missing_ok=True)
    record = {
        "schema_version": 2,
        "url": url,
        "final_url": final_url,
        "method": "POST" if post_body is not None else "GET",
        "request_body": post_body,
        "sha256": content_hash,
        "bytes": total,
        "retrieved_at": now(),
        "license_basis": license_basis,
        "citation": citation,
        "upstream_sha256": expected_sha256,
        "code_sha": sha,
    }
    acquisition_id = uuid.uuid4().hex
    write_json(state / "data" / "acquisitions" / f"{acquisition_id}.json", record)
    assign_role(content_hash, role, rationale, sha=sha, state=state)
    event(
        "data_acquired",
        sha=sha,
        details={"acquisition_id": acquisition_id, "sha256": content_hash},
        state=state,
    )
    return record | {"acquisition_id": acquisition_id}
