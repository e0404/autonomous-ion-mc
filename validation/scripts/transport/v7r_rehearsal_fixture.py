"""Builds the committed V7-R rehearsal fixture ``tests/ionmc/fixtures/v7r/rehearsal/`` from the host archive
(V3-005C C8d; plan Amendment 14a (b)2 and Amendment 17 (a)7).

    uv run --extra dev python validation/scripts/transport/v7r_rehearsal_fixture.py \
        --archive validation/generated --run-ids RUN-A RUN-B ... --code-sha <40 hex> \
        [--seed-base 20505000] [--out tests/ionmc/fixtures/v7r/rehearsal] [--force]

The rehearsal shards are produced by the controlled generation path of ``run_suite.py``: a run of the suite
``lv5c`` whose ``--only`` selection consists solely of ``v7r-s{k}`` steps and whose ``--seed-base`` lies in
the rehearsal family 20505000-20509499 is exempt from the fixture preflight. This script scans the host run
directories named by ``--run-ids`` (``<archive>/RUN-*/transport/*/samples/`` and ``<archive>/RUN-*/transport/*/``)
for ``v7r-s{k}.json`` and the npz sidecar bound in each partial, copies them byte-identically, requires every
partial's ``git_sha`` to equal ``--code-sha``, writes ``PROVENANCE.json`` (sorted keys: ``host_run_ids``,
``partial_digests``, ``hashed_source_digest`` of the CURRENT tree, ``seed_base``, ``code_sha``, ``head_sha``, ``clean_tree``, numpy and python
versions, ``created_from``) and then runs ``v7r_rehearsal.load_rehearsal``, the verification of the CI test,
exiting non-zero on any violation (nothing is left in ``--out`` in that case). The builder refuses unless the
repository HEAD is exactly ``--code-sha`` and ``git status --porcelain`` (untracked included) is empty, so the
recorded digest is that of the commit's hashed sources. The fixture can only be produced at the code whose
hashed-source digest it records. An existing ``--out`` is refused unless ``--force``.
"""

from __future__ import annotations

import argparse
import json
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_OUT = REPO / "tests" / "ionmc" / "fixtures" / "v7r" / "rehearsal"
DEFAULT_SEED_BASE = 20505000


def _fail(msg: str) -> None:
    raise SystemExit(f"v7r_rehearsal_fixture: {msg}")


def _git(*args: str) -> str:
    """Output of ``git <args>`` in ``REPO`` (module-level so tests can point it at a temporary repository)."""
    res = subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        _fail(f"git {' '.join(args)} failed: {res.stderr.strip()}")
    return res.stdout


def require_clean_head(code_sha: str) -> str:
    """Refuses unless HEAD is exactly ``code_sha`` and the tree is clean (untracked files included), so the
    recorded hashed-source digest is the digest of that commit's hashed sources. Returns HEAD."""
    head = _git("rev-parse", "HEAD").strip()
    if head != code_sha:
        _fail(f"repository HEAD {head} differs from --code-sha {code_sha}")
    dirty = _git("status", "--porcelain", "--untracked-files=all").strip()
    if dirty:
        _fail(f"working tree is not clean (the digest would not be that of {code_sha}):\n{dirty}")
    return head


def scan_partials(archive: Path, run_ids: list[str]) -> list[Path]:
    """The ``v7r-s{k}.json`` partials of the named host runs (samples directory or suite directory)."""
    hits: list[Path] = []
    for rid in run_ids:
        run = archive / rid
        if not run.is_dir():
            _fail(f"host run directory {run} not found")
        found = sorted({*run.glob("transport/*/samples/v7r-s*.json"), *run.glob("transport/*/v7r-s*.json")})
        found = [p for p in found if p.name.removeprefix("v7r-s").removesuffix(".json").isdigit()]
        if not found:
            _fail(f"no v7r-s{{k}}.json partial under {run}")
        hits += found
    return hits


def build(archive: Path, run_ids: list[str], code_sha: str, seed_base: int, out: Path, force: bool) -> dict:  # type: ignore[type-arg]
    import numpy as np

    if out.exists() and not force:
        _fail(f"{out} exists (use --force to replace it)")
    if len(code_sha) != 40 or any(c not in "0123456789abcdef" for c in code_sha):
        _fail("--code-sha must be the full 40-hex-digit lowercase SHA")
    head = require_clean_head(code_sha)
    hits = scan_partials(archive, run_ids)
    names = [p.name for p in hits]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        _fail(f"shard partial(s) {dup} found more than once under {run_ids}")
    sys.path.insert(0, str(HERE))
    import run_suite
    import steps_v5c
    import v7r_rehearsal

    # staged outside the repository: a directory inside tests/ionmc would enter the hashed-source digest
    stage = Path(tempfile.mkdtemp(prefix="rehearsal-stage-"))
    try:
        digests: dict[str, str] = {}
        created_from: list[str] = []
        for path in hits:
            doc = json.loads(path.read_bytes())
            if doc.get("git_sha") != code_sha:
                _fail(f"{path}: git_sha {doc.get('git_sha')!r} differs from --code-sha {code_sha}")
            shutil.copyfile(path, stage / path.name)
            created_from.append(str(path.relative_to(archive)))
            digests[path.name] = doc.get("content_sha256", "")
            side = (doc.get("sidecar") or {}).get("file")
            if not isinstance(side, str) or side != Path(side).name:
                _fail(f"{path}: no plain sidecar file name in the partial")
            src = next((c for c in (path.parent / side, path.parent.parent / side) if c.is_file()), None)
            if src is None:
                _fail(f"{path}: sidecar {side} not found next to the partial or in its parent directory")
            shutil.copyfile(src, stage / side)
            created_from.append(str(src.relative_to(archive)))
        prov = {
            "host_run_ids": sorted(run_ids), "partial_digests": digests,
            "hashed_source_digest": run_suite.hashed_source_digest("lv5c"), "seed_base": seed_base,
            "code_sha": code_sha, "head_sha": head, "clean_tree": True, "numpy_version": np.__version__,
            "python_version": platform.python_version(), "created_from": sorted(created_from),
        }  # fmt: skip
        (stage / "PROVENANCE.json").write_text(json.dumps(prov, sort_keys=True, indent=1) + "\n")
        try:
            _, info = v7r_rehearsal.load_rehearsal(stage, steps_v5c, run_suite)
        except v7r_rehearsal.RehearsalError as exc:
            _fail(f"verification failed: {exc}")
        if force and out.exists():
            shutil.rmtree(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(stage, out)
        return info
    finally:
        shutil.rmtree(stage, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("--archive", required=True, help="host archive root, e.g. validation/generated")
    ap.add_argument("--run-ids", nargs="+", required=True, help="host run directories (RUN-...) to scan")
    ap.add_argument("--code-sha", required=True, help="full git sha every partial must carry")
    ap.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--force", action="store_true", help="replace an existing fixture directory")
    a = ap.parse_args(argv)
    info = build(Path(a.archive), a.run_ids, a.code_sha, a.seed_base, Path(a.out), a.force)
    print(f"{a.out}: verified, {info['shards']} shards, hashed_source_digest {info['hashed_source_digest']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
