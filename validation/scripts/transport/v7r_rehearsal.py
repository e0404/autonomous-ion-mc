"""Loader and verifier of the committed V7-R rehearsal fixture (``tests/ionmc/fixtures/v7r/rehearsal/``).

Shared by the CI tests (``tests/ionmc/test_v7r_coverage.py``) and the builder
``v7r_rehearsal_fixture.py``, so that a fixture is produced and verified by the same code. The
requirements enforced are listed in the module docstring of ``test_v7r_coverage.py`` (Amendment 14a
(b)2, Amendment 17 (a)7). ``steps`` is the ``steps_v5c`` module and ``run_suite`` the ``run_suite``
module (passed in, so that this file imports neither).
"""

# ruff: noqa: E501

from __future__ import annotations

import json
from pathlib import Path
from types import ModuleType

import numpy as np

EXPECTED_SHARD = {
    "n": 4_500_000,
    "replicates": 450,
    "batch_histories": 500,
    "blocks_per_replicate": 20,
}
SEED_BASE_RANGE = (20505000, 20509499)


class RehearsalError(Exception):
    """The rehearsal fixture violates at least one requirement; ``str`` lists every violation."""


def load_rehearsal(
    directory: Path, steps: ModuleType, run_suite: ModuleType, digest: str | None = None
) -> tuple[np.ndarray, dict]:  # type: ignore[type-arg]
    """Verified pooled block sums ``[replicates, 20]`` of the rehearsal in ``directory`` and a
    record of what was verified, or ``RehearsalError`` listing every violation (module docstring).
    Seed reading (Amendment 14 (j)): the rehearsal BASE lies in 20505000-20509499 and shard seeds
    are ``base + 16000 + k`` (20521000+ for base 20505000). These coincide with the declared future
    hr5c base 20521004 only as a number, never as a used seed (hr5c seeds are
    ``base + 1000 r + k``, r >= 1).
    ``digest`` overrides the current hashed-source digest (tests only)."""
    bad: list[str] = []
    v7r = steps.v7r
    prov_path = directory / "PROVENANCE.json"
    prov: dict = {}  # type: ignore[type-arg]
    if not prov_path.is_file():
        bad.append(f"PROVENANCE.json missing in {directory}")
    else:
        try:
            prov = json.loads(prov_path.read_text())
        except ValueError as exc:
            bad.append(f"PROVENANCE.json unreadable: {exc}")
    if prov:
        for key in (
            "host_run_ids",
            "partial_digests",
            "hashed_source_digest",
            "seed_base",
            "code_sha",
        ):
            if key not in prov:
                bad.append(f"PROVENANCE.json lacks {key}")
        ids = prov.get("host_run_ids")
        if not (isinstance(ids, list) and ids and all(isinstance(i, str) and i for i in ids)):
            bad.append("PROVENANCE.json host_run_ids must be a non-empty list of run ids")
        if not (
            isinstance(prov.get("seed_base"), int)
            and SEED_BASE_RANGE[0] <= prov["seed_base"] <= SEED_BASE_RANGE[1]
        ):
            bad.append(
                f"PROVENANCE.json seed_base {prov.get('seed_base')!r} outside {SEED_BASE_RANGE}"
            )
        if "head_sha" in prov and prov["head_sha"] != prov.get("code_sha"):
            bad.append(
                f"PROVENANCE.json head_sha {prov['head_sha']!r} differs from code_sha {prov.get('code_sha')!r}"
            )
        have = digest if digest is not None else run_suite.hashed_source_digest("lv5c")
        if prov.get("hashed_source_digest") != have:
            bad.append(
                "hashed_source_digest differs from the current hashed source set of lv5c "
                "(the rehearsal is not at the evidence code, Amendment 17 (a)7)"
            )
    parts = sorted(directory.glob("v7r-s*.json"))
    if len(parts) < 2:
        bad.append(f"{len(parts)} rehearsal partial(s) in {directory}, need >= 2 shards")
    recorded = prov.get("partial_digests") if isinstance(prov.get("partial_digests"), dict) else {}
    if prov and {p.name for p in parts} != set(recorded):
        bad.append(f"partial_digests names {sorted(recorded)} differ from the partials present")
    blocks, shards, seeds = [], [], []
    for path in parts:
        nm = path.name
        try:
            doc = json.loads(path.read_text())
        except ValueError as exc:
            bad.append(f"{nm}: unreadable: {exc}")
            continue
        dig = steps.v5.content_digest(doc)
        if doc.get("content_sha256") != dig:
            bad.append(f"{nm}: content_sha256 does not reproduce")
        if nm in recorded and recorded[nm] != dig:
            bad.append(f"{nm}: differs from the digest recorded in PROVENANCE.json")
        if doc.get("valid") is not True:
            bad.append(f"{nm}: valid is not true")
        if doc.get("reduced") is not False:
            bad.append(f"{nm}: reduced is not false")
        for key, want in EXPECTED_SHARD.items():
            if doc.get(key) != want:
                bad.append(f"{nm}: {key} is {doc.get(key)!r}, expected {want}")
        shards.append(doc.get("shard"))
        seeds.append(doc.get("seed"))
        sb = doc.get("seed_base")
        if not (isinstance(sb, int) and SEED_BASE_RANGE[0] <= sb <= SEED_BASE_RANGE[1]):
            bad.append(f"{nm}: seed_base {sb!r} outside {SEED_BASE_RANGE}")
        if prov and sb != prov.get("seed_base"):
            bad.append(f"{nm}: seed_base differs from PROVENANCE.json")
        if (
            isinstance(sb, int)
            and isinstance(doc.get("shard"), int)
            and (doc.get("seed") != sb + 1000 * v7r.V7R_R_INDEX + doc["shard"])
        ):
            bad.append(f"{nm}: seed is not seed_base + 1000 * 16 + shard")
        if prov and doc.get("git_sha") != prov.get("code_sha"):
            bad.append(f"{nm}: git_sha differs from code_sha of PROVENANCE.json")
        try:
            blocks.append(steps.load_v7r_sidecar(doc, path))
        except SystemExit as exc:
            bad.append(f"{nm}: sidecar rejected: {exc}")
    if len(set(shards)) != len(shards):
        bad.append(f"shard indices not unique: {shards}")
    if len(set(seeds)) != len(seeds):
        bad.append(f"seeds not unique: {seeds}")
    if bad:
        raise RehearsalError("; ".join(bad))
    info = {"shards": len(parts), **{k: prov[k] for k in
            ("seed_base", "code_sha", "host_run_ids", "hashed_source_digest")}}  # fmt: skip
    return np.concatenate(blocks), info
