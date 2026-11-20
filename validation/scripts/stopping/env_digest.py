# ruff: noqa: E501
"""Environment digests of the V3-013 seal (stdlib only; plan ``v3-013-acceptance.md`` section D step 3).

Usage (the accepted launch is through ``sealed_launch.py``; direct launch also works for a local
dry run, the record then carries ``"launch": null``)::

    .venv/bin/python -I -S -B -X pycache_prefix="$(mktemp -d)" \
        validation/scripts/stopping/sealed_launch.py validation/scripts/stopping/env_digest.py \
        --write validation/results/stopping/v3-013/seal-environment.json
    ... env_digest.py --check validation/results/stopping/v3-013/seal-environment.json

One digest per component, each the sha256 over the sorted list of (relative path, file sha256)
(:func:`env_seal.tree_digest`), except the interpreter and ``uv.lock`` (the file sha256):

* ``interpreter``: the resolved real executable (not the venv symlink);
* ``stdlib``: every file of the ``sysconfig`` ``stdlib`` directory and, when different, of
  ``platstdlib`` (extension modules included; ``site-packages``, ``dist-packages``, ``test`` and
  ``__pycache__`` excluded, exactly the walker of ``env_seal``);
* ``numpy`` and ``warp``: the COMPLETE installed distributions (``numpy``, ``warp-lang``): every file
  listed in the distribution's ``RECORD`` (``importlib.metadata``) plus every file of every
  ``*.libs`` directory next to the package in site-packages (vendored native libraries such as
  ``numpy.libs``). The distributions are located by metadata, not imported;
* ``uv_lock``: the sha256 of the repository's ``uv.lock``.

The record holds the digests, the versions, ``sys.version``, ``platform.platform()`` and the number
of files per component; it holds no absolute path. ``environment_digests()`` returns that record
(used by ``steps_v13.py`` and the sealed evaluator). :func:`warp_kernel_cache_attestation` hashes
every generated source and shared object of a Warp kernel cache directory.
"""

from __future__ import annotations

import argparse
import builtins
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SCHEMA = "ionmc-v13-environment-digests-1"
COMPONENTS = ("interpreter", "stdlib", "numpy", "warp", "uv_lock")
#: keys of the record that must be equal for two environments to be the same
COMPARE_KEYS = ("schema", "digests", "versions", "file_counts", "sys_version", "platform")
DISTRIBUTIONS = {"numpy": "numpy", "warp": "warp-lang"}
WARP_CACHE_SCHEMA = "ionmc-v13-warp-kernel-cache-1"
LAUNCH_RECORD = "__ionmc_sealed_launch__"


def _env_seal() -> Any:
    """``env_seal.py`` next to this file (stdlib only at import), shared with the I3-FU-1 seal."""
    mod = sys.modules.get("env_seal")
    if mod is not None and Path(getattr(mod, "__file__", "")).resolve() == HERE / "env_seal.py":
        return mod
    spec = importlib.util.spec_from_file_location("env_seal", HERE / "env_seal.py")
    if spec is None or spec.loader is None:
        raise ImportError(HERE / "env_seal.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["env_seal"] = mod
    spec.loader.exec_module(mod)
    return mod


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def component_entries(repo: Path | None = None) -> dict[str, dict[str, str]]:
    """{component: {relative path: sha256}} of every component (the interpreter and ``uv.lock``
    hold one entry each)."""
    es = _env_seal()
    repo = REPO if repo is None else Path(repo)
    stdlib = {f"stdlib/{k}": v for k, v in es.directory_files(es.stdlib_dir()).items()}
    plat = es.platstdlib_dir()
    if plat != es.stdlib_dir():
        stdlib.update({f"platstdlib/{k}": v for k, v in es.directory_files(plat).items()})
    lock = repo / "uv.lock"
    if lock.is_symlink() or not lock.is_file():
        raise es.EnvironmentSealError(f"{lock} is not a regular file")
    return {
        "interpreter": {"python": _sha256_file(Path(sys.executable).resolve())},
        "stdlib": stdlib,
        "numpy": es.distribution_files(DISTRIBUTIONS["numpy"]),
        "warp": es.distribution_files(DISTRIBUTIONS["warp"]),
        "uv_lock": {"uv.lock": _sha256_file(lock)},
    }


def environment_digests(repo: Path | None = None) -> dict[str, Any]:
    """The environment record (see the module docstring): ``digests`` per component,
    ``file_counts``, ``versions``, ``sys_version`` and ``platform``."""
    es = _env_seal()
    entries = component_entries(repo)
    digests = {}
    for name in COMPONENTS:
        e = entries[name]
        digests[name] = (
            next(iter(e.values())) if name in ("interpreter", "uv_lock") else (es.tree_digest(e))
        )
    return {
        "schema": SCHEMA,
        "digests": digests,
        "file_counts": {name: len(entries[name]) for name in COMPONENTS},
        "versions": {
            "python": platform.python_version(),
            **{n: importlib.metadata.version(d) for n, d in DISTRIBUTIONS.items()},
        },
        "sys_version": sys.version,
        "platform": platform.platform(),
    }


def differences(a: Mapping[str, Any], b: Mapping[str, Any]) -> list[str]:
    """Names of the components / keys in which two environment records differ (empty: equal);
    a missing or malformed record differs in every key."""
    out: list[str] = []
    for key in COMPARE_KEYS:
        va, vb = a.get(key), b.get(key)
        if key in ("digests", "file_counts", "versions"):
            if not isinstance(va, Mapping) or not isinstance(vb, Mapping):
                out.append(key)
                continue
            out += [f"{key}.{k}" for k in sorted(set(va) | set(vb)) if va.get(k) != vb.get(k)]
        elif va != vb:
            out.append(key)
    return out


def warp_kernel_cache_attestation(cache_dir: Path | str) -> dict[str, Any]:
    """sha256 of every file (generated source and shared object) under the Warp kernel cache
    ``cache_dir`` (relative path -> sha256), their count and the tree digest."""
    root = Path(cache_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"Warp kernel cache directory {root} does not exist")
    files = {
        p.relative_to(root).as_posix(): _sha256_file(p)
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }
    return {
        "schema": WARP_CACHE_SCHEMA,
        "count": len(files),
        "digest": _env_seal().tree_digest(files),
        "files": files,
    }


def _launch_summary() -> dict[str, Any] | None:
    rec = getattr(builtins, LAUNCH_RECORD, None)
    if not isinstance(rec, dict):
        return None
    keys = ("origin", "launcher_sha256", "target", "target_sha256", "head_sha", "target_committed")
    return {k: rec.get(k) for k in keys}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Environment digests of the V3-013 seal")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", metavar="PATH", help="compute the record and write it")
    mode.add_argument("--check", metavar="SEAL", help="recompute and compare with a written record")
    a = ap.parse_args(argv)
    if a.write:
        rec = {**environment_digests(), "launch": _launch_summary()}
        out = Path(a.write)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(rec, indent=1, sort_keys=True) + "\n")
        print(f"wrote {out}")
        for name in COMPONENTS:
            print(f"  {name}: {rec['digests'][name]} ({rec['file_counts'][name]} files)")
        return 0
    try:
        seal = json.loads(Path(a.check).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        print(f"cannot read {a.check}: {exc}", file=sys.stderr)
        return 2
    if not isinstance(seal, dict):
        print(f"{a.check} is not a JSON object", file=sys.stderr)
        return 2
    diff = differences(environment_digests(), seal)
    if diff:
        print(f"environment differs from {a.check} in: {', '.join(diff)}", file=sys.stderr)
        return 1
    print("environment matches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
