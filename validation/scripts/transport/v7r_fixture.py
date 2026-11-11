"""Derives the committed calibration fixture of V7-R from the archived lv5b V7 replicate partials
(V3-005C C5; plan Amendment 14 (h): "Surrogate from the archived 7200 lv5b replicates (7aae5bb
partials, (mean, sem) pairs)").

    python validation/scripts/transport/v7r_fixture.py [--archive validation/generated] \
        [--out tests/ionmc/fixtures/v7r/lv5b-7aae5bb-escaped-neutral-replicates.json]

Reads ``v7-rep-s0.json`` .. ``v7-rep-s7.json`` (shards, 900 replicates each, seeds 20482004-20482011)
and ``v7-rep-ref.json`` (seed 20482012) from the host archives ``<archive>/RUN-*/transport/
lv5b-7aae5bb-*/samples/``. It refuses a document whose ``git_sha`` does not start with ``7aae5bb1``,
whose ``content_sha256`` does not reproduce, that is reduced or invalid, or whose seed/shard/count is
not the expected one, and writes a deterministic JSON file (sorted keys, no timestamps): the 7200
``escaped_neutral`` replicate means and SEMs in shard order (per primary, 20 blocks of 500 histories)
with the provenance of every source (host run directory, file sha256, ``content_sha256``, seed).
Pure numpy-free standard library; no ``ionmc`` import.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
GIT_SHA_PREFIX = "7aae5bb1"
SCHEMA = "ionmc-v7r-lv5b-escaped-neutral-replicates-1"
DEFAULT_OUT = REPO / "tests" / "ionmc" / "fixtures" / "v7r" / "lv5b-7aae5bb-escaped-neutral-replicates.json"
SHARDS, REPS, FIRST_SEED = 8, 900, 20482004


def content_digest(doc: dict[str, Any]) -> str:
    """As ``steps_v5.content_digest`` (canonical JSON of the document without ``content_sha256``)."""
    norm = json.loads(json.dumps(doc))
    norm.pop("content_sha256", None)
    return hashlib.sha256(json.dumps(norm, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _find(archive: Path, name: str) -> Path:
    hits = sorted(glob.glob(str(archive / "RUN-*" / "transport" / "lv5b-7aae5bb-*" / "samples" / name)))
    if len(hits) != 1:
        raise SystemExit(f"{name}: found {len(hits)} copies under {archive} (need exactly 1)")
    return Path(hits[0])


def derive(archive: Path) -> dict[str, Any]:
    means: list[float] = []
    sems: list[float] = []
    sources = []
    for k in range(SHARDS + 1):
        name = f"v7-rep-s{k}.json" if k < SHARDS else "v7-rep-ref.json"
        path = _find(archive, name)
        raw = path.read_bytes()
        doc = json.loads(raw)
        if not str(doc.get("git_sha", "")).startswith(GIT_SHA_PREFIX):
            raise SystemExit(f"{name}: git_sha {doc.get('git_sha')!r} is not {GIT_SHA_PREFIX}...")
        if doc.get("content_sha256") != content_digest(doc):
            raise SystemExit(f"{name}: content_sha256 does not reproduce")
        if doc.get("reduced") or not doc.get("valid"):
            raise SystemExit(f"{name}: reduced or invalid")
        if doc.get("seed") != FIRST_SEED + k or (k < SHARDS and doc.get("shard") != k):
            raise SystemExit(f"{name}: unexpected seed/shard")
        if k < SHARDS:
            est = doc["estimators"]["escaped_neutral"]
            m = [row[0] for row in est["mean"]]
            s = [row[0] for row in est["sem"]]
            if len(m) != REPS or len(s) != REPS:
                raise SystemExit(f"{name}: {len(m)} replicates, expected {REPS}")
            means += m
            sems += s
        sources.append({
            "name": name, "host_run_dir": path.parents[3].name, "file_sha256": hashlib.sha256(raw).hexdigest(),
            "content_sha256": doc["content_sha256"], "seed": doc["seed"], "git_sha": doc["git_sha"],
            "n": doc["n"], "seed_base": doc["seed_base"],
        })  # fmt: skip
    return {
        "schema": SCHEMA, "estimator": "escaped_neutral", "git_sha_prefix": GIT_SHA_PREFIX,
        "replicates": SHARDS * REPS, "blocks_per_replicate": 20, "batch_histories": 500,
        "units": "MeV per primary (mean and SEM of 20 blocks of 500 histories)",
        "order": "shards 0..7 in order, 900 replicates each", "sources": sources,
        "mean": means, "sem": sems,
    }  # fmt: skip


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--archive", default=str(REPO / "validation" / "generated"))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    a = ap.parse_args(argv)
    doc = derive(Path(a.archive))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, sort_keys=True, indent=0) + "\n")
    print(f"{out} sha256 {hashlib.sha256(out.read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
