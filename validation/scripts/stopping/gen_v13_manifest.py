# ruff: noqa: E501
"""Frozen manifest of the V3-013 seal (stdlib only; plan ``v3-013-acceptance.md`` section D step 3).

Usage::

    python -I -B validation/scripts/stopping/gen_v13_manifest.py --write [--out PATH]
    python -I -B validation/scripts/stopping/gen_v13_manifest.py --check [--manifest PATH]

Writes ``validation/reference_cases/topas/v13-frozen-manifest.json``: the sha256 of every sealed
path and ``src_ionmc_tree`` = ``git rev-parse HEAD:src/ionmc`` (``src/ionmc`` itself is covered by
the tree hash). The manifest does not list itself (a file cannot contain its own digest): it is
bound by its blob at the signed seal SHA and by its sha256 in the plan's pin note. Sealed paths:

* the sealed evaluator ``v13_evaluate.py`` and its transitive in-repository import closure
  (the AST of every file is walked for ``import`` / ``from ... import`` under ``validation/`` and
  ``infrastructure/``, plus string constants naming a repository script (``*.py`` or a name with an underscore), because the evaluator
  loads some scripts by path; ``ionmc`` is covered by the tree hash);
* ``sealed_launch.py``, ``env_digest.py`` (and ``env_seal.py`` through the closure),
  ``pyproject.toml``, ``uv.lock``, ``infrastructure/keys/trusted_attestation_keys.json``;
* ``validation/results/stopping/v3-013/seal-environment.json`` (when present: the controlled
  host-run record of ``env_digest.py``), ``vc-predictions.json``, ``steps_v13.py``;
* every file of the V3-013 TOPAS bundles (``*-v13-*`` directories) and their S5a manifest.

``--check`` recomputes the manifest and fails (exit 1) on any difference from the file.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SCHEMA = "ionmc-v13-frozen-manifest-1"
MANIFEST_REL = "validation/reference_cases/topas/v13-frozen-manifest.json"
EVALUATOR_REL = "validation/scripts/stopping/v13_evaluate.py"
LAUNCHER_REL = "validation/scripts/stopping/sealed_launch.py"
ENV_DIGEST_REL = "validation/scripts/stopping/env_digest.py"
STEPS_REL = "validation/scripts/stopping/steps_v13.py"
TRUST_REL = "infrastructure/keys/trusted_attestation_keys.json"
RESULTS_DIR = "validation/results/stopping/v3-013"
ENV_SEAL_REL = f"{RESULTS_DIR}/seal-environment.json"
PREDICTIONS_REL = f"{RESULTS_DIR}/vc-predictions.json"
BUNDLES_REL = "validation/reference_cases/topas"
BUNDLES_MANIFEST_REL = f"{BUNDLES_REL}/v13-bundles-manifest.json"
FIXED = (
    EVALUATOR_REL, LAUNCHER_REL, ENV_DIGEST_REL, STEPS_REL, TRUST_REL, PREDICTIONS_REL,
    "pyproject.toml", "uv.lock",
)  # fmt: skip
OPTIONAL = (ENV_SEAL_REL,)
SCOPES = ("validation", "infrastructure")
_STEM = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(repo: Path, *args: str) -> str:
    proc = subprocess.run(  # noqa: S603 - fixed argv, no shell
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=False
    )
    if proc.returncode != 0:
        raise SystemExit(f"git {' '.join(args)} failed: {proc.stderr.strip()}")
    return proc.stdout.strip()


def _py_files(repo: Path) -> list[Path]:
    out = []
    for scope in SCOPES:
        out += [
            p
            for p in sorted((repo / scope).rglob("*.py"))
            if "__pycache__" not in p.parts and not p.is_symlink()
        ]
    return out


def _dotted_candidates(repo: Path, dotted: str) -> list[Path]:
    """Files of the in-repository module ``dotted`` (``infrastructure.a.b``) and the
    ``__init__.py`` of its parent packages; ``ionmc`` is covered by the tree hash."""
    parts = dotted.split(".")
    if parts[0] != "infrastructure":
        return []
    found = []
    for n in range(1, len(parts) + 1):
        base = repo.joinpath(*parts[:n])
        if (base / "__init__.py").is_file():
            found.append(base / "__init__.py")
    leaf = repo.joinpath(*parts).with_suffix(".py")
    if leaf.is_file():
        found.append(leaf)
    return found


def _references(repo: Path, path: Path, stems: dict[str, list[Path]]) -> set[Path]:
    """In-repository files named by the imports and the string constants of ``path``."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except SyntaxError as exc:
        raise SystemExit(f"cannot parse {path}: {exc}") from exc
    found: set[Path] = set()

    def add_module(name: str) -> None:
        found.update(_dotted_candidates(repo, name))
        if "." not in name:
            found.update(stems.get(name, []))

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                add_module(a.name)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            add_module(node.module)
            for a in node.names:
                if node.module.split(".")[0] == "infrastructure":
                    add_module(f"{node.module}.{a.name}")
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            s = node.value
            if s.startswith("infrastructure.") and _STEM.fullmatch(s.replace(".", "_")):
                add_module(s)
                continue
            stem = s[:-3] if s.endswith(".py") else s
            if _STEM.fullmatch(stem) and (s.endswith(".py") or "_" in stem):
                found.update(stems.get(stem, []))
    found.discard(path)
    return found


def import_closure(repo: Path, roots: list[str]) -> list[str]:
    """Repository-relative paths of ``roots`` and everything they (transitively) import or name."""
    stems: dict[str, list[Path]] = {}
    for p in _py_files(repo):
        stems.setdefault(p.stem, []).append(p)
    seen: set[Path] = set()
    todo = [repo / r for r in roots]
    while todo:
        p = todo.pop()
        if p in seen:
            continue
        if not p.is_file():
            raise SystemExit(f"closure root {p} does not exist")
        seen.add(p)
        todo += sorted(_references(repo, p, stems) - seen)
    return sorted(p.relative_to(repo).as_posix() for p in seen)


def sealed_paths(repo: Path) -> list[str]:
    """Every path of the manifest (sorted, unique)."""
    out = set(import_closure(repo, [EVALUATOR_REL, ENV_DIGEST_REL, STEPS_REL]))
    out.update(FIXED)
    out.update(rel for rel in OPTIONAL if (repo / rel).is_file())
    for d in sorted((repo / BUNDLES_REL).glob("*-v13-*")):
        if d.is_dir():
            out.update(p.relative_to(repo).as_posix() for p in sorted(d.iterdir()) if p.is_file())
    if (repo / BUNDLES_MANIFEST_REL).is_file():
        out.add(BUNDLES_MANIFEST_REL)
    missing = sorted(r for r in out if not (repo / r).is_file() or (repo / r).is_symlink())
    if missing:
        raise SystemExit(f"sealed paths that are not regular files: {missing}")
    return sorted(out)


def build(repo: Path, require_committed: bool = False) -> dict:
    """The manifest document of ``repo`` (hashes of the on-disk bytes, tree hash of HEAD)."""
    paths = sealed_paths(repo)
    dirty = _git(repo, "status", "--porcelain", "--", "src/ionmc")
    if dirty:
        raise SystemExit(
            f"src/ionmc has uncommitted changes; its tree hash would not match:\n{dirty}"
        )
    if require_committed:
        pending = _git(repo, "status", "--porcelain", "--", *paths)
        if pending:
            raise SystemExit(f"sealed paths with uncommitted changes:\n{pending}")
    return {
        "schema": SCHEMA,
        "plan": "validation/plans/v3-013-acceptance.md section D step 3",
        "note": "the manifest does not list itself; it is bound by its blob at the signed seal SHA "
        "and by its sha256 in the plan pin note",
        "src_ionmc_tree": _git(repo, "rev-parse", "HEAD:src/ionmc"),
        "file_count": len(paths),
        "files": {r: sha256_file(repo / r) for r in paths},
    }


def text(doc: dict) -> str:
    return json.dumps(doc, indent=1, sort_keys=True) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Frozen manifest of the V3-013 seal")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--out", default=None, help="--write: output file (default: the manifest path)")
    ap.add_argument("--manifest", default=None, help="--check: manifest file (default: the path)")
    ap.add_argument("--repo", default=str(REPO), help="repository root (tests)")
    ap.add_argument(
        "--require-committed", action="store_true", help="refuse uncommitted sealed paths"
    )
    a = ap.parse_args(argv)
    repo = Path(a.repo).resolve()
    doc = build(repo, a.require_committed)
    new = text(doc)
    if a.write:
        out = Path(a.out) if a.out else repo / MANIFEST_REL
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(new)
        print(
            f"wrote {out}: {doc['file_count']} sealed paths, src/ionmc tree {doc['src_ionmc_tree']}"
        )
        print(f"manifest sha256 {hashlib.sha256(new.encode()).hexdigest()}")
        return 0
    path = Path(a.manifest) if a.manifest else repo / MANIFEST_REL
    try:
        old = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        print(f"cannot read {path}: {exc}", file=sys.stderr)
        return 2
    bad = []
    if old.get("src_ionmc_tree") != doc["src_ionmc_tree"]:
        bad.append("src_ionmc_tree")
    of, nf = old.get("files", {}), doc["files"]
    bad += [f"{k} ({'missing' if k not in nf else 'unlisted' if k not in of else 'hash'})"
            for k in sorted(set(of) | set(nf)) if of.get(k) != nf.get(k)]  # fmt: skip
    if path.read_text(encoding="utf-8") != new:
        bad.append("file text (formatting or other fields)") if not bad else None
    if bad:
        print(f"{path} differs from the recomputed manifest in: {bad}", file=sys.stderr)
        return 1
    print(f"manifest matches: {doc['file_count']} sealed paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
